#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
indicators.py —— 路线无关光化学指示因子引擎 (五路线屋脊实验核心)
  输入: <pref>_conc.csv + <pref>_irr.csv + 机制 .fac (FORMULA 段 + 计量)
  输出: post-spinup 全部输出行均值指示因子 dict
    官方指示因子 5 族 —— 判据、图件、run_B 诊断只用这一组:
      浓度比  FNR = HCHO/NO2 ; H2O2_HNO3
      过程量  L_NOx_Q ; P_H2O2_HNO3 ; OPE = P(O3)/P(NOz)   (由 IRR 逐反应组装)
    状态量    NOy_ppb —— 25 物种含 N 加和；不是比值，用于跨情景可比与观测对照
    备用扩展  O3_NOy ; HCHO_NOy ; O3_NOz —— 仍计算、仍在 batch 24 列内，
              但不进判据/图件（与官方 5 族信息重叠，历史上属 Sillman 全家福）
    绝对量    Q ; L_NOx ; O3_ppb
  NOy 定义: 从 .fac FORMULA 段自动取含 N 物种 (M3 单一事实源的直接收益),
           N2O5 计双 N; NOz = NOy − NO − NO2 (−N2O5 快平衡项不计产率)
  用法:
    from indicators import route_indicators
    d = route_indicators(pref, fac)          # 窗口均值
    ts = route_timeseries(pref, fac)         # 逐时序列 (进阶诊断)

    O3 收支 (论文 Eq4-6) 用 Mechanism.p_o3_gd(peff_row)：
      G=HO2+NO+ΣRO2+NO ; D=O1D+H2O+OH+NO2+OH+O3+HO2+O3+Σalkene+O3 ; P=G-D
"""
import csv
import re
import numpy as np

PPB = 2.46e10

# V2.5+ 冻结机制的 RO2 物种（Eq(5) 的 Σ RO2_i+NO）
# 机制改版新增 RO2 时必须同步更新；p_o3_gd() 启动时校验物种存在
PEROXY_SPECIES = frozenset({
    "HO2", "CH3O2", "CH3CO3", "ISOPBO2", "ISOPNO3", "MVKO2", "MACRO2",
    "OLERO2", "ALKRO2", "ARO1O2", "ARO2O2", "KETRO2", "TERPRO2",
    "CSLRO2", "LALKRO2", "MTBERO2", "ESTERO2", "SESQRO2", "ETHO2",
    "ALD2O2", "IOLERO2",
})
ALKENE_SPECIES = frozenset({"OLE", "ISOP", "TERP", "SESQ", "ETH", "IOLE"})


def _side_terms(side):
    """反应侧 → {物种: 系数}; 如 '0.90 NO2 + 0.10 ISOPNIT'"""
    d = {}
    for tok in side.split("+"):
        tok = tok.strip()
        if not tok:
            continue
        mm = re.match(r"^(\d+\.?\d*)\s+(\w+)$", tok)
        c, nm = (float(mm.group(1)), mm.group(2)) if mm else (1.0, tok)
        d[nm] = d.get(nm, 0.0) + c
    return d


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
    """物种序 / FORMULA(N、C 数) / 反应表 (速率表达式, lhs, rhs)"""
    species, noy_w, carbon_w, rxns = [], {}, {}, []
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
                        if not mm:
                            continue
                        if mm.group(1) == "N":
                            noy_w[tok[0]] = int(mm.group(2))
                        elif mm.group(1) == "C":
                            carbon_w[tok[0]] = int(mm.group(2))
            if s.endswith(";"):
                in_form = False
            continue
        if s.startswith("%"):
            b, eq = body[1:].split(":", 1)     # body 已剥 ';' — 尾分号毁所有名字匹配
            lhs, rhs = eq.split("=", 1)
            rxns.append((b.strip(), lhs.strip(), rhs.strip()))
    return species, noy_w, carbon_w, rxns


def _rate_tokens(expr):
    """速率表达式 → (命名速率引用列表); 供 L_NOx/Q 语义匹配"""
    return re.findall(r"\b[A-Za-z_]\w*\b", expr)


class Mechanism:
    """IRR → 过程量的映射 (与 fac_batch_omp 预算表同一语义, Python 版)"""

    def __init__(self, fac_path):
        self.species, self.noy_w, self.carbon_w, rxns = _parse_fac(fac_path)
        self.rxns = rxns
        self.sid = {s: i for i, s in enumerate(self.species)}
        self.nrx = len(rxns)
        self.prod = [[] for _ in rxns]     # (isp, coef)
        self.reac = [[] for _ in rxns]
        for k, (_, lhs, rhs) in enumerate(rxns):
            for side, store in ((lhs, self.reac), (rhs, self.prod)):
                for nm, c in _side_terms(side).items():
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
        # NOz 单一定义: NOy 减去本元组（与 fac_batch_omp 的 NOz 口径一致）
        self.nox_res = ("NO", "NO2", "N2O5")
        noz_spec = [s for s in self.noy_w if s not in self.nox_res]
        for k in range(self.nrx):
            for isp, c in self.prod[k]:
                if self.species[isp] in noz_spec:
                    self.pnoz_rx.append((k, c))
                    break

        # O3 收支 (论文 Eq4-6): G = HO2+NO + Σ RO2+NO; D = 论文五类汇
        missing = PEROXY_SPECIES - set(self.species)
        if missing:
            raise ValueError(f"PEROXY_SPECIES 不在机制物种表: {sorted(missing)}")
        self.g_ho2_rx, self.g_ro2_rx = [], []
        self.d_o1d_h2o_rx, self.d_ohno2_rx = [], []
        self.d_oho3_rx, self.d_ho2o3_rx, self.d_alkene_o3_rx = [], [], []
        self.r_o3no_rx, self.r_o3no2_rx, self.r_o3hv_rx, self.o3_loss_rx = [], [], [], []
        for k, (rate, lhs, rhs) in enumerate(rxns):
            ld, rd = _side_terms(lhs), _side_terms(rhs)
            if ld.get("NO", 0) > 0 and rd.get("NO2", 0) > 0:
                # 一条 RO2+NO 最多氧化 1 个 NO；额外 NO2 视为硝酸根释放
                coef = min(rd["NO2"], 1.0)
                if "HO2" in ld:
                    self.g_ho2_rx.append((k, coef))
                elif any(sp in ld for sp in PEROXY_SPECIES):
                    self.g_ro2_rx.append((k, coef))
            if set(ld) == {"O1D", "H2O"} and set(rd) == {"OH"}:
                self.d_o1d_h2o_rx.append(k)
            if set(ld) == {"OH", "NO2"} and set(rd) == {"HNO3"}:
                self.d_ohno2_rx.append(k)
            if set(ld) == {"OH", "O3"} and set(rd) == {"HO2", "O2"}:
                self.d_oho3_rx.append(k)
            if set(ld) == {"HO2", "O3"} and set(rd) == {"OH", "O2"}:
                self.d_ho2o3_rx.append(k)
            if "O3" in ld and any(sp in ld for sp in ALKENE_SPECIES):
                self.d_alkene_o3_rx.append(k)
            if set(ld) == {"O3", "NO"} and set(rd) == {"NO2"}:
                self.r_o3no_rx.append(k)
            if set(ld) == {"O3", "NO2"} and set(rd) == {"NO3"}:
                self.r_o3no2_rx.append(k)
            if "O3" in ld:
                self.o3_loss_rx.append(k)
                if rate in ("J<1>", "J<2>"):
                    self.r_o3hv_rx.append(k)

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

    def p_o3_gd(self, peff_row):
        """论文 Eq(4)-(6): G、D_paper、P_paper; O3+NO/NO2/hv 为独立诊断

        IRR 行是窗口积分，故返回值也是“窗口累计量”；调用方负责 /dtout/PPB。
        """
        g_ho2 = sum(c * peff_row[k] for k, c in self.g_ho2_rx)
        g_ro2 = sum(c * peff_row[k] for k, c in self.g_ro2_rx)
        d_o1d = sum(peff_row[k] for k in self.d_o1d_h2o_rx)
        d_ohno2 = sum(peff_row[k] for k in self.d_ohno2_rx)
        d_oho3 = sum(peff_row[k] for k in self.d_oho3_rx)
        d_ho2o3 = sum(peff_row[k] for k in self.d_ho2o3_rx)
        d_alkene = sum(peff_row[k] for k in self.d_alkene_o3_rx)
        d_paper = d_o1d + d_ohno2 + d_oho3 + d_ho2o3 + d_alkene
        g_total = g_ho2 + g_ro2
        r_o3no = sum(peff_row[k] for k in self.r_o3no_rx)
        r_o3no2 = sum(peff_row[k] for k in self.r_o3no2_rx)
        r_o3hv = sum(peff_row[k] for k in self.r_o3hv_rx)
        r_total = sum(peff_row[k] for k in self.o3_loss_rx)
        r_other = r_total - r_o3no - r_o3no2 - r_o3hv - d_oho3 - d_ho2o3 - d_alkene
        return {
            "G_HO2NO": g_ho2, "G_RO2NO": g_ro2, "G_total": g_total,
            "D_O1D_H2O": d_o1d, "D_OH_NO2": d_ohno2, "D_OH_O3": d_oho3,
            "D_HO2_O3": d_ho2o3, "D_alkene_O3": d_alkene, "D_paper": d_paper,
            "P_paper": g_total - d_paper,
            "R_O3_NO": r_o3no, "R_O3_NO2": r_o3no2, "R_O3_hv": r_o3hv,
            "R_O3_other": r_other, "R_O3_total": r_total,
        }


def route_timeseries(pref, fac, mech=None):
    """post-spinup 全部输出行的逐时指示因子; 返回 dict of arrays"""
    m = mech or Mechanism(fac)
    chdr, C = _load_csv(f"{pref}_conc.csv")
    ihrd, I = _load_csv(f"{pref}_irr.csv")   # I[:,0] = h — budgets 前必须剥掉
    sid = m.sid
    idx = {s: chdr.index(s) for s in m.species if s in chdr}
    T = C.shape[0]
    # 顺序即 run_E CSV 列序: 官方 5 族 → 状态量 NOy → 备用 3 比值 → 绝对量
    out = {k: [] for k in ["FNR", "H2O2_HNO3", "L_NOx_Q", "P_H2O2_HNO3", "OPE",
                           "NOy_ppb",
                           "O3_NOy", "HCHO_NOy", "O3_NOz",
                           "Q", "L_NOx", "O3_ppb"]}
    noy_arr = np.zeros(T)
    for t in range(T):
        noy = sum(w * max(C[t, idx[s]], 0.0) for s, w in m.noy_w.items() if s in idx)
        noy_arr[t] = noy
        noz = noy - sum(max(C[t, idx[s]], 0.0) for s in m.nox_res if s in idx)
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


def route_indicators(pref, fac, mech=None, pm_only=False, dtout=3600.0):
    """post-spinup 全部输出行的窗口均值 (pm_only=True: 12–18h, Sillman 口径)

    conc CSV 的 h 列是**输出步序号**（fac_v3 写 step），当地小时 =
    h × dtout/3600 % 24，故 dtout 由调用方从 runtime 的 MODEL 行传入
    （缺省 3600 = 逐时输出，此时 h 与小时数值相同）。
    不假设窗口起点 —— 窗口自 h48 起还是 h96 起都正确。
    """
    ts = route_timeseries(pref, fac, mech)
    d = {}
    sel = None
    if pm_only:
        _, C = _load_csv(f"{pref}_conc.csv")     # C[:,0] = h 列 (输出步序号)
        hrs = (C[:, 0] * dtout / 3600.0) % 24.0
        sel = np.array([12.0 <= h <= 18.0 for h in hrs])
    for k, v in ts.items():
        if pm_only:
            v = v[sel]
        d[k] = float(np.mean(v)) if v.size else float("nan")
    return d
