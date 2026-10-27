#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
indicators.py —— 路线无关光化学指示因子引擎 (五路线屋脊实验核心)
  输入: <pref>_conc.csv + <pref>_irr.csv + 机制 .fac (FORMULA 段 + 计量)
  输出: 末 24h 窗口均值指示因子 dict
    第一类 (浓度比): FNR, H2O2_HNO3, O3_NOy, HCHO_NOy, O3_NOz, NOy_ppb
    第二类 (过程量, 由 IRR 逐反应组装): Q, L_NOx, L_NOx_Q, P_H2O2, P_HNO3,
           PH2O2_PHNO3, OPE, P_O3, P_NOz
  NOy 定义: 从 .fac FORMULA 段自动取含 N 物种 (M3 单一事实源的直接收益),
           N2O5 计双 N; NOz = NOy − NO − NO2 (−N2O5 快平衡项不计产率)
  用法:
    from indicators import route_indicators
    d = route_indicators(pref, fac)          # 窗口均值
    ts = route_timeseries(pref, fac)         # 逐时序列 (进阶诊断)
"""
import csv
import re
import numpy as np

PPB = 2.46e10


def _fix_fortran_exp(tok):
    if "E" not in tok and "e" not in tok:
        m = re.search(r"(\d)-(\d{3})$", tok)
        if m:
            return tok[:m.start(1) + 1] + "E" + tok[m.start(1) + 1:]
    return tok


def _load_csv(path):
    rows = list(csv.reader(open(path)))
    hdr = [h.strip().split()[0] for h in rows[0]]   # 剥行内残留空白
    data = np.array([[float(_fix_fortran_exp(x)) for x in r] for r in rows[1:] if r])
    return hdr, data


def _parse_fac(fac_path):
    """物种序 / FORMULA(含 N 物种与 N 数) / 反应表 (速率表达式, lhs, rhs)"""
    species, noy_w, rxns = [], {}, []
    in_vars = in_form = False
    for raw in open(fac_path):
        s = raw.strip()
        if not s or s.startswith("*"):
            continue
        if s == "VARIABLE":
            in_vars = True; continue
        if s.startswith("FORMULA"):
            in_form = True; continue
        body = s.rstrip(";").strip()
        if in_vars:
            species += body.split()
            if s.endswith(";"):
                in_vars = False
            continue
        if in_form:
            if body:
                tok = body.split()
                if tok[0] in species:
                    for e in tok[1:]:
                        mm = re.match(r"^([A-Z][a-z]?)(\d+)$", e)
                        if mm and mm.group(1) == "N":
                            noy_w[tok[0]] = int(mm.group(2))
            if s.endswith(";"):
                in_form = False
            continue
        if s.startswith("%"):
            b, eq = body[1:].split(":", 1)     # body 已剥 ';' — 尾分号毁所有名字匹配
            lhs, rhs = eq.split("=", 1)
            rxns.append((b.strip(), lhs.strip(), rhs.strip()))
    return species, noy_w, rxns


def _rate_tokens(expr):
    """速率表达式 → (命名速率引用列表); 供 L_NOx/Q 语义匹配"""
    return re.findall(r"\b[A-Za-z_]\w*\b", expr)


class Mechanism:
    """IRR → 过程量的映射 (与 fac_batch_omp 预算表同一语义, Python 版)"""

    def __init__(self, fac_path):
        self.species, self.noy_w, rxns = _parse_fac(fac_path)
        self.sid = {s: i for i, s in enumerate(self.species)}
        self.nrx = len(rxns)
        self.prod = [[] for _ in rxns]     # (isp, coef)
        self.reac = [[] for _ in rxns]
        for k, (_, lhs, rhs) in enumerate(rxns):
            for side, store in ((lhs, self.reac), (rhs, self.prod)):
                for tok in side.split("+"):
                    tok = tok.strip()
                    if not tok:
                        continue
                    mm = re.match(r"^(\d+\.?\d*)\s+(\w+)$", tok)
                    c, nm = (float(mm.group(1)), mm.group(2)) if mm else (1.0, tok)
                    if nm in self.sid:
                        store[k].append((self.sid[nm], c))
        # 预算: Q (自由基源) / L_NOx / P(H2O2) / P(HNO3) 的反应集合
        # Q 匹配 = lhs/rhs 名多重集 (FACSIMILE 反应物在计量侧; 与 Fortran 预算表同语义)
        self.q_rx = []      # (k, coef)
        self.lnox_rx = []
        self.pnoz_rx = []   # 产物含 N 的 NOz 物种反应
        q_patterns = [
            ((["O1D", "H2O"], ["OH", "OH"]), 2.0),
            ((["HONO"], ["NO", "OH"]), 1.0),
            ((["H2O2"], ["OH", "OH"]), 2.0),
            ((["HCHO"], ["HO2", "HO2", "CO"]), 2.0),
            ((["CH3OOH"], ["CH3O", "OH"]), 1.0),
            ((["ISOPOOH"], ["ISOPBO", "OH"]), 1.0),
        ]
        for k, (_, lhs, rhs) in enumerate(rxns):
            lset = sorted(t.strip() for t in lhs.split("+") if t.strip())
            rset = sorted(t.strip() for t in rhs.split("+") if t.strip())
            for (lp, rp), c in q_patterns:
                if lset == sorted(lp) and rset == sorted(rp):
                    self.q_rx.append((k, c))
                    break
            if self._is_lnox(lhs, rhs):
                self.lnox_rx.append(k)
        noz_spec = [s for s, w in self.noy_w.items() if s not in ("NO", "NO2", "N2O5")]
        for k in range(self.nrx):
            for isp, c in self.prod[k]:
                if self.species[isp] in noz_spec:
                    self.pnoz_rx.append((k, c))
                    break

    @staticmethod
    def _is_lnox(lhs, rhs):
        """NOx→储库: NO2+OH→HNO3 / NO+OH→HONO (名多重集, 双反应物)"""
        l = [t.strip() for t in lhs.split("+") if t.strip()]
        r = [t.strip() for t in rhs.split("+") if t.strip()]
        pairs = [(["NO2", "OH"], ["HNO3"]), (["NO", "OH"], ["HONO"])]
        for lp, rp in pairs:
            if sorted(l) == sorted(lp) and sorted(r) == sorted(rp):
                return True
        return False

    def budgets(self, peff_row):
        """一行 IRR (nr,) → Q, L_NOx, P_H2O2, P_HNO3, P_O3, P_NOz"""
        q = sum(c * peff_row[k] for k, c in self.q_rx)
        ln = sum(peff_row[k] for k in self.lnox_rx)
        ph = sum(cf * peff_row[k] for k in range(self.nrx) for i, cf in self.prod[k]
                 if self.species[i] == "H2O2")
        phn = sum(cf * peff_row[k] for k in range(self.nrx) for i, cf in self.prod[k]
                  if self.species[i] == "HNO3")
        po3 = sum(cf * peff_row[k] for k in range(self.nrx) for i, cf in self.prod[k]
                  if self.species[i] == "O3")
        pnz = sum(c * peff_row[k] for k, c in self.pnoz_rx)
        return q, ln, ph, phn, po3, pnz


def route_timeseries(pref, fac, mech=None):
    """逐时指示因子 (末 24h); 返回 dict of arrays"""
    m = mech or Mechanism(fac)
    chdr, C = _load_csv(f"{pref}_conc.csv")
    ihrd, I = _load_csv(f"{pref}_irr.csv")   # I[:,0] = h — budgets 前必须剥掉
    sid = m.sid
    idx = {s: chdr.index(s) for s in m.species if s in chdr}
    T = C.shape[0]
    out = {k: [] for k in ["FNR", "H2O2_HNO3", "O3_NOy", "HCHO_NOy", "O3_NOz",
                           "NOy_ppb", "Q", "L_NOx", "L_NOx_Q", "P_H2O2_HNO3",
                           "OPE", "O3_ppb"]}
    noy_arr = np.zeros(T)
    for t in range(T):
        noy = sum(w * max(C[t, idx[s]], 0.0) for s, w in m.noy_w.items() if s in idx)
        noy_arr[t] = noy
        noz = noy - max(C[t, idx.get("NO", 0)], 0.0) - max(C[t, idx.get("NO2", 0)], 0.0)
        q, ln, ph, phn, po3, pnz = m.budgets(I[t, 1:])   # off-by-one 修复: 剥 h 列
        r = lambda a, b: a / b if b > 1e-30 else 0.0
        out["FNR"].append(r(C[t, idx["HCHO"]], C[t, idx["NO2"]]))
        out["H2O2_HNO3"].append(r(C[t, idx["H2O2"]], C[t, idx["HNO3"]]))
        out["O3_NOy"].append(r(C[t, idx["O3"]], noy))
        out["HCHO_NOy"].append(r(C[t, idx["HCHO"]], noy))
        out["O3_NOz"].append(r(C[t, idx["O3"]], noz))
        out["NOy_ppb"].append(noy / PPB)
        out["Q"].append(q); out["L_NOx"].append(ln)
        out["L_NOx_Q"].append(r(ln, q))
        out["P_H2O2_HNO3"].append(r(ph, phn))
        out["OPE"].append(r(po3, pnz))
        out["O3_ppb"].append(C[t, idx["O3"]] / PPB)
    return {k: np.array(v) for k, v in out.items()}


def route_indicators(pref, fac, mech=None, pm_only=False):
    """末 24h 窗口均值 (pm_only=True: 12–18h, Sillman 口径)

    pm 窗口的小时取自 conc CSV 的 h 列 (h % 24; dtout=3600 时 h 即小时),
    不假设窗口起点 —— --days 3 (窗口自 h48 起) 与默认 5 天 (h96 起) 均正确。
    """
    ts = route_timeseries(pref, fac, mech)
    d = {}
    sel = None
    if pm_only:
        _, C = _load_csv(f"{pref}_conc.csv")     # C[:,0] = h 列 (输出步序号)
        hrs = C[:, 0] % 24.0
        sel = np.array([12.0 <= h <= 18.0 for h in hrs])
    for k, v in ts.items():
        if pm_only:
            v = v[sel]
        d[k] = float(np.mean(v)) if v.size else float("nan")
    return d
