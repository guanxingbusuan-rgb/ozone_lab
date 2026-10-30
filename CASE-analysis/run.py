#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CASE-analysis/run.py —— 个例案例分析（单一入口：模拟 + 全部 CSV 存档）

三类过程分析合一（六账账本 / 反应级归因 / 前体物 RIR）的自完备程序：
一个「个例」= 本目录 config_override/ 的一份完整配置快照。

自包含: 不 import 任何编排层代码；共享仅"模型本体"四件套
        （facsimile_impl/fac_v3、inputproc/cfgbuild.py、机制 .fac、
         analysis_lib/csvio 与 indicators）。
输入: config_override/ 完整配置集（本目录即唯一输入源；直接改这里即换个例）
输出（本目录 output/）:
  in_case/ + runtime_case.dat         配置快照与编译产物（可追溯的第一事实源）
  case_{conc,ipr,irr}.csv             基线 fac_v3 原始输出（h 列 = 输出步序号）
  ipr_ledger.csv                      六账长表（h × 物种 × P/L/E/D/B/X）
  irr_top.csv / hourly_rxn_top.csv    反应级净贡献（默认 O3）
  rir_mean.csv / rir_diurnal.csv      前体物 RIR（论文 Eq4–8，P(O3) 加权）
  rir_species_{mean,diurnal}.csv      逐物种 RIR（仅 --per-species）
  p_o3no_baseline.csv                 基线 O3 收支（Eq4–6）
  runs.csv                            基线+RIR 扰动运行记录
  run.log / run_meta.txt              运行日志与元数据

用法:
  python CASE-analysis/run.py                       # 基线 + 六账 + IRR + RIR(-10%)
  python CASE-analysis/run.py --no-rir              # 只做基线/六账/IRR（快）
  python CASE-analysis/run.py --pct 20 --groups NOx,CO --window 8-18
  python CASE-analysis/run.py --species HCHO --per-species
  python CASE-analysis/run.py --no-assim        # 关闭观测同化（自由运行基线）
"""
import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "output")
SRC = os.path.join(HERE, "config_override")
FACV3 = os.path.join(ROOT, "facsimile_impl", "fac_v3")
CFGBUILD = os.path.join(ROOT, "inputproc", "cfgbuild.py")
MECH_PATH = os.path.join(ROOT, "mechanism", "mechanism_3_full_chain.fac")
GEN = os.path.join(ROOT, "facsimile_impl", "gen_v3")

sys.path.insert(0, os.path.join(ROOT, "analysis_lib"))
from csvio import PPB, read_table, write_csv          # noqa: E402
from indicators import Mechanism                      # noqa: E402

LOG_F = None


def log(msg=""):
    print(msg)
    if LOG_F:
        LOG_F.write(msg + "\n")
        LOG_F.flush()


# ── 前体物分组（与 D3 原定义一致）──
VOC_CLASSES = ["ALK", "OLE", "ARO", "ALD", "KET", "OVOC", "ISOP", "TERP", "DMS"]
GROUPS = {
    "NOx": ["NO", "NO2"], "CO": ["CO"],
    "ALK": ["ALK", "LALK"], "OLE": ["OLE", "ETH", "IOLE"],
    "ARO": ["ARO1", "ARO2", "CSL"], "ALD": ["HCHO", "CH3CHO", "ALD2"],
    "KET": ["KET"], "OVOC": ["MEOH", "ETOH", "MTBE", "ESTER"],
    "ISOP": ["ISOP"], "TERP": ["TERP", "SESQ"], "DMS": ["DMS"],
}
VOC_SPECIES = sorted({s for g in VOC_CLASSES for s in GROUPS[g]})
GROUPS["TVOC"] = VOC_SPECIES
PER_SPECIES = sorted(set(VOC_SPECIES) | {"NO", "NO2", "CO"})
BUDGET_KEYS = ["G_HO2NO", "G_RO2NO", "G_total",
               "D_O1D_H2O", "D_OH_NO2", "D_OH_O3", "D_HO2_O3", "D_alkene_O3",
               "D_paper", "P_paper", "R_O3_NO", "R_O3_NO2", "R_O3_hv",
               "R_O3_other", "R_O3_total"]


# ── 观测同化约束（数据同化 bc 通道）──────────────────────────
# 输入约定: config_override/obs/ 目录下
#   obs_<物种名>.txt      物种直测序列（头部 # unit: ppb；两列 t[s] 值）
#   obsclass_<类>.txt     类别总量序列（ppbC），按权重拆到类内物种
#   tvoc_obs.txt          TVOC 总量序列（ppbC），按基线排放谱拆到全部 VOC
# 类别 → 机制物种映射（权重 = in_emission 基线排放谱归一，随快照存档可追溯）:
#   烷烃 ALKANES → ALK 0.60 / LALK 0.40
#   烯烃 ALKENES → OLE 0.70 / ETH 0.15 / IOLE 0.15
#   芳烃 AROMATICS → ARO1 0.51 / ARO2 0.40 / CSL 0.09
#   OVOC → KET 0.38 / HCHO 0.15 / ALD2 0.15 / MEOH 0.12 / ETOH 0.08 /
#          MTBE 0.05 / ESTER 0.04 / CH3CHO 0.03
#   炔烃 / 卤代烃: 机制 V2.5+ 无对应物种 → 不可约束（观测仅作验证对比）
# 优先级: 物种直测 > 类别拆分 > TVOC 拆分（高一级命中后，低一级跳过该物种，防双计）。
OBS_DIR = os.path.join(SRC, "obs")
CLASS_WEIGHTS = {
    "ALKANES":   {"ALK": 0.60, "LALK": 0.40},
    "ALKENES":   {"OLE": 0.70, "ETH": 0.15, "IOLE": 0.15},
    "AROMATICS": {"ARO1": 0.51, "ARO2": 0.40, "CSL": 0.09},
    "OVOC":      {"KET": 0.38, "HCHO": 0.15, "ALD2": 0.15, "MEOH": 0.12,
                  "ETOH": 0.08, "MTBE": 0.05, "ESTER": 0.04, "CH3CHO": 0.03},
}
TVOC_WEIGHTS = {"ALK": 0.267, "OLE": 0.253, "LALK": 0.178, "ARO1": 0.173,
                "ARO2": 0.136, "ETH": 0.054, "IOLE": 0.054, "ISOP": 0.035,
                "CSL": 0.030, "KET": 0.025, "HCHO": 0.011, "MEOH": 0.009,
                "ALD2": 0.005, "ETOH": 0.005, "MTBE": 0.002, "DMS": 0.001,
                "ESTER": 0.001, "TERP": 0.001, "SESQ": 0.001, "CH3CHO": 0.001}


def _read_obs_series(path):
    """观测序列: 两列 t[s] 值（值单位 ppb 或 ppbC）；跳过 # 注释"""
    t, v = [], []
    for ln in open(path, encoding="utf-8"):
        s = ln.split("#")[0].split()
        if len(s) < 2:
            continue
        t.append(float(s[0])); v.append(float(s[1]))
    return t, v


def _write_bc_series(cfg_dir, sp, t, v, src_note):
    dst = os.path.join(cfg_dir, "constraints", "obs_" + sp + ".txt")
    with open(dst, "w", encoding="utf-8") as f:
        f.write("# 观测约束序列 (来源: " + src_note + ")\n# unit: ppb\n# t[s]  " + sp + "\n")
        for tt, vv in zip(t, v):
            f.write(str(int(round(tt))).rjust(9) + "  " + format(vv, "10.4f") + "\n")
    return dst


def prepare_observations(cfg_dir):
    """扫描 config_override/obs/ → 生成逐物种观测序列 + 改写 bc.config（同化装配）。

    返回 (constrained 物种列表, 装配报告行)。无 obs/ 目录时返回空（自由运行）。"""
    if not os.path.isdir(OBS_DIR):
        return [], []
    direct, classes, tvoc = {}, {}, None
    for fn in sorted(os.listdir(OBS_DIR)):
        if fn.startswith("obsclass_") and fn.endswith(".txt"):
            classes[fn[len("obsclass_"):-4].upper()] = os.path.join(OBS_DIR, fn)
        elif fn.startswith("obs_") and fn.endswith(".txt"):
            direct[fn[len("obs_"):-4].upper()] = os.path.join(OBS_DIR, fn)
        elif fn == "tvoc_obs.txt":
            tvoc = os.path.join(OBS_DIR, fn)
    mech_species = set(Mechanism(MECH_PATH).species)

    per_species = {}
    report = []
    for name, path in sorted(direct.items()):
        if name not in mech_species:
            report.append("  [obs跳过] " + name + ": 不在机制 VARIABLE 段")
            continue
        t, v = _read_obs_series(path)
        per_species[name] = (t, v, "直测 obs_" + name + ".txt")
    for cname, path in sorted(classes.items()):
        w = CLASS_WEIGHTS.get(cname)
        if not w:
            report.append("  [obs跳过] obsclass_" + cname + ": 未注册的类别")
            continue
        t, v = _read_obs_series(path)
        for sp, wt in w.items():
            if sp in per_species:      # 物种直测优先，防双计
                continue
            per_species[sp] = (t, [x * wt for x in v], "类别 " + cname + " x" + format(wt, ".2f"))
        report.append("  [obs类别] " + cname + ": 拆分至 " + str(len(w)) + " 物种（ppbC 权重和=1）")
    if tvoc:
        t, v = _read_obs_series(tvoc)
        n_used = 0
        for sp, wt in TVOC_WEIGHTS.items():
            if sp in per_species or sp not in mech_species:
                continue
            per_species[sp] = (t, [x * wt for x in v], "TVOC x" + format(wt, ".3f"))
            n_used += 1
        report.append("  [obsTVOC] 拆分至 " + str(n_used) + " 个未覆盖 VOC 物种")

    for sp, (t, v, note) in sorted(per_species.items()):
        _write_bc_series(cfg_dir, sp, t, v, note)
    bc = os.path.join(cfg_dir, "bc.config")
    kept = [ln for ln in open(bc, encoding="utf-8").read().splitlines()
            if ln.strip().startswith("TAU_MIX") or ln.strip().startswith("#")
            or not ln.strip()]
    if not any("TAU_MIX" in ln for ln in kept):
        kept.insert(0, "TAU_MIX = 86400")
    for sp in sorted(per_species):
        kept.append(sp + " = file:constraints/obs_" + sp + ".txt")
    open(bc, "w", encoding="utf-8").write("\n".join(kept) + "\n")
    report.append("  [obs装配] " + str(len(per_species)) + " 个物种约束 → bc.config（τ=TAU_MIX 全局）")
    return sorted(per_species), report


# ── 自包含驱动：快照/编译/运行/窗口 ──
def prepare_case(cfg_dir):
    if not os.path.isdir(SRC):
        raise SystemExit(f"✗ 缺输入目录 {SRC}（个例 = 该目录的完整配置集）")
    if os.path.exists(cfg_dir):
        shutil.rmtree(cfg_dir)
    shutil.copytree(SRC, cfg_dir, ignore=shutil.ignore_patterns("__pycache__"))
    return cfg_dir


def compile_runtime(cfg_dir, rt_out):
    r = subprocess.run([sys.executable, CFGBUILD, cfg_dir, rt_out, "--mech", MECH_PATH],
                       capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(rt_out):
        log(r.stdout[-1200:] + r.stderr[-600:])
        raise SystemExit("✗ cfgbuild 编译失败")
    return rt_out


def runtime_window(rt):
    for ln in open(rt):
        t = ln.split()
        if len(t) >= 2 and t[0] == "MODEL":
            d, dt = float(t[1]), float(t[2])
            sp = float(t[5]) if len(t) >= 6 else 0.0
            if not (d > sp):
                raise SystemExit(f"✗ 模拟天数 {d:g} ≤ SPINUP {sp:g}（无 post-spinup 输出）")
            return d, dt, sp
    raise SystemExit(f"✗ runtime 缺 MODEL 行: {rt}")


def run_fac(scen, pref, rt):
    t0 = time.time()
    r = subprocess.run([FACV3, scen, pref, rt], cwd=os.path.dirname(FACV3),
                       capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(f"{pref}_conc.csv"):
        log(r.stdout[-800:] + r.stderr[-600:])
        raise SystemExit(f"✗ fac_v3 失败: {scen} → {pref}")
    m = re.search(r"nst=\s*(\d+)", r.stdout)
    return {"nst": int(m.group(1)) if m else -1,
            "elapsed_s": time.time() - t0,
            "line": next((x.strip() for x in r.stdout.splitlines() if "FAC-V3 ok" in x), "")}


def scale_emissions(cfg_dir, species, scale):
    """把 emissions.config 中指定物种的基线列乘 scale（保留单位行与注释）"""
    p = os.path.join(cfg_dir, "emissions.config")
    out = []
    for ln in open(p, encoding="utf-8"):
        body, sep, com = ln.partition("#")
        t = body.split()
        if t and t[0] in species and len(t) >= 2:
            body = body.replace(t[1], f"{float(t[1]) * scale:.6e}", 1)
        out.append(body + (sep + com if sep else ""))
    open(p, "w", encoding="utf-8").write("".join(out))


def emitted_species():
    out = set()
    for ln in open(os.path.join(SRC, "emissions.config")):
        s = ln.split("#")[0].split()
        if s and s[0].lower() != "unit":
            out.add(s[0])
    return out


# ── 分析核 ──
# fac_v3 的 CSV 首列 h 是**输出步序号**（oz_fac_v3.f90 写 step），不是小时：
# 第 h 步的行代表后向窗 [h·dtout/3600 − dtout/3600, h·dtout/3600]（小时）。
# dtout=3600 时两者数值相同，所以只有改输出步长才会暴露混用。统一走下面三个函数。
def win_hours(hstep, dtout):
    """→ (窗口起点小时, 窗口终点小时)"""
    end = hstep * dtout / 3600.0
    return end - dtout / 3600.0, end


def win_hour(hstep, dtout):
    """窗口起点所在的当地小时 0-23（逐时输出→24 个唯一标签；子小时输出→同段多窗共一标签）"""
    return int(win_hours(hstep, dtout)[0] // 1) % 24


def win_day(hstep, dtout):
    """窗口起点落在模拟的第几天（1-based）"""
    return int(win_hours(hstep, dtout)[0] // 24.0) + 1


def export_ipr(pref, outdir, dtout):
    hp, DP = read_table(f"{pref}_ipr.csv")
    species, j = [], 1
    while j < len(hp) and hp[j].split("_")[0] != "fnox":
        species.append(hp[j].split("_")[0])
        j += 6
    rows = []
    for r in DP:
        hh = int(r[0])
        hl = win_hour(hh, dtout)              # 本窗起点的当地小时
        for si, sp in enumerate(species):
            b = 1 + 6 * si
            rows.append([hh, hl, sp] + [f"{r[b + q]:.6e}" for q in range(6)])
    p = write_csv(os.path.join(outdir, "ipr_ledger.csv"),
                  ["h", "hlocal", "species", "P", "L", "E", "D", "B", "X"], rows,
                  note="个例分析·六账 | P机制产 L机制耗 E排放 D干湿沉 B同化+稀释 X BLH膨胀 | "
                       "单位 molec/cm3，每行=本窗积分（勿再 trapz）| "
                       f"h=输出步序号（非小时）, dtout={dtout:g}s; "
                       "hlocal=本窗起点当地小时 0-23")
    return p, species, len(DP)


def export_irr(pref, outdir, species, mech, dtout):
    jSP = mech.species.index(species) + 1
    reac = [ln.split() for ln in open(os.path.join(GEN, "oz.reac"))]
    prod = [ln.split() for ln in open(os.path.join(GEN, "oz.prod"))]
    eqs = [l.strip().lstrip("%").strip() for l in open(MECH_PATH) if l.startswith("%")]
    coef = np.zeros(len(reac))
    for k in range(len(reac)):
        pr, rr = prod[k], reac[k]
        c = sum(float(pr[2 + 2 * t]) for t in range(int(pr[0]))
                if int(pr[1 + 2 * t]) == jSP)
        c -= sum(1.0 for t in range(int(rr[0])) if int(rr[1 + t]) == jSP)
        coef[k] = c
    hi, DI = read_table(f"{pref}_irr.csv")
    I = np.array([r[1:1 + len(reac)] for r in DI], dtype=float)
    W = I * coef
    s = W.sum(axis=0)
    rows = [[k + 1, f"{s[k]:.5e}", eqs[k][:70] if k < len(eqs) else "?"]
            for k in np.argsort(-np.abs(s))[:12]]
    p1 = write_csv(os.path.join(outdir, "irr_top.csv"),
                   ["r_index", f"net_{species}_contribution", "reaction"], rows,
                   note=f"个例分析·IRR | {species} 净贡献 top-12 | 正=净生成 负=净消耗 | 窗口积分求和")
    rh = []
    for i, r in enumerate(DI):
        hlocal = win_hour(int(r[0]), dtout)
        for rank, k in enumerate(np.argsort(-np.abs(W[i]))[:4], 1):
            rh.append([int(r[0]), hlocal, rank, k + 1, f"{W[i][k]:.4e}",
                       eqs[k][:64] if k < len(eqs) else "?"])
    p2 = write_csv(os.path.join(outdir, "hourly_rxn_top.csv"),
                   ["h", "hlocal", "rank", "r_index", f"net_{species}", "reaction"], rh,
                   note=f"个例分析·IRR | 每窗 {species} 净贡献 top-4 | h=输出步序号, "
                        f"hlocal=本窗起点当地小时, dtout={dtout:g}s")
    return p1, p2


def _parse_window(s):
    if s.lower() == "all":
        return 0.0, 23.0
    lo, hi = (float(x) for x in s.split("-", 1))
    if not (0.0 <= lo <= hi <= 23.0):
        raise SystemExit("--window 需为 all 或 0-23 内的 lo-hi")
    return lo, hi


def _s_value(row, chdr, species, m, ppbc):
    tot, found = 0.0, False
    for sp in species:
        if sp not in chdr:
            continue
        c = row[chdr.index(sp)] / PPB
        if ppbc:
            c *= m.carbon_w.get(sp, 1.0)
        tot += c
        found = True
    return tot if found else np.nan


def _p_ppbh(m, irr_row, dtout):
    return m.p_o3_gd(irr_row)["P_paper"] * 3600.0 / (dtout * PPB)


def _calc_rows(group, species, chdr, Cb, Ib, Cp, Ip, m, dtout, lo, hi, tmin_h=0.0):
    ppbc = group not in ("NOx", "CO", "NO", "NO2")
    rows = []
    for t in range(len(Ib)):
        hstep = Ib[t, 0]
        if win_hours(hstep, dtout)[0] < tmin_h - 1e-9:
            continue            # 跨 spinup 边界的首窗不参与诊断（与 plot.py 的剔除同规则）
        hour = float(win_hour(hstep, dtout))      # 本窗起点当地小时
        if not (lo - 1e-9 <= hour <= hi + 1e-9):
            continue
        pb, pp = _p_ppbh(m, Ib[t, 1:], dtout), _p_ppbh(m, Ip[t, 1:], dtout)
        sb = _s_value(Cb[t], chdr, species, m, ppbc)
        sp = _s_value(Cp[t], chdr, species, m, ppbc)
        ds = sb - sp
        ok = (np.isfinite(pb) and np.isfinite(pp) and np.isfinite(sb)
              and np.isfinite(sp) and pb > 0 and sb > 0 and ds / sb > 1e-12)
        rows.append({"group": group, "day": win_day(hstep, dtout),
                     "hour": hour, "P_base": pb, "RIR":
                     ((pb - pp) / pb) / (ds / sb) if ok else np.nan,
                     "valid": bool(ok), "basis": "ppbC" if ppbc else "ppb"})
    return rows


def _weighted(rows):
    v = [r for r in rows if r["valid"]]
    w = sum(r["P_base"] for r in v)
    if not v or w <= 0:
        return np.nan, w, 0
    return sum(r["RIR"] * r["P_base"] for r in v) / w, w, len(v)


def rir_tables(rows, order):
    """RIR 行池 → (逐时表, 全窗表)。s_basis 取该池计算时真正用的口径，不重新推断。"""
    d, w = [], []
    for g in order:
        grp = [r for r in rows if r["group"] == g]
        for hh in range(24):
            sel = [r for r in grp if int(round(r["hour"])) % 24 == hh]
            rir, wsp, n = _weighted(sel)
            if n:
                d.append([g, hh, f"{rir:.6e}", f"{wsp:.6e}",
                          len({r["day"] for r in sel if r["valid"]}), n,
                          sel[0]["basis"]])
        rir, wsp, n = _weighted(grp)
        if n:
            w.append([g, f"{rir:.6e}", f"{wsp:.6e}", n, grp[0]["basis"]])
    return d, w


def export_rir(base_pref, outdir, rt, args, runs):
    days, dtout, spinup = runtime_window(rt)
    m = Mechanism(MECH_PATH)
    emitted = emitted_species()
    chdr, Cb = read_table(f"{base_pref}_conc.csv")
    _, Ib = read_table(f"{base_pref}_irr.csv")
    Cb, Ib = np.array(Cb, float), np.array(Ib, float)
    scale = 1.0 - args.pct / 100.0
    names = list(GROUPS) if args.groups.lower() == "all" else \
        [x.strip() for x in args.groups.split(",") if x.strip()]
    bad = [x for x in names if x not in GROUPS]
    if bad:
        raise SystemExit(f"✗ 未知组 {bad}; 可用 {list(GROUPS)}")
    lo, hi = _parse_window(args.window)

    def one(group, sp_list):
        sp_list = [s for s in sp_list if s in emitted]
        if not sp_list:
            log(f"  [skip] {group}: emissions.config 无该组物种")
            return []
        cfg = os.path.join(OUT, "rir_runs", f"in_{group}")
        prepare_case(cfg)
        prepare_observations(cfg)          # 同化约束与基线一致（唯一差异=前体物减排）
        scale_emissions(cfg, sp_list, scale)
        rt_g = compile_runtime(cfg, os.path.join(OUT, "rir_runs", f"rt_{group}.dat"))
        pref = os.path.join(OUT, "rir_runs", group)
        r = run_fac(args.scen, pref, rt_g)
        runs.append([group, len(sp_list), f"{scale:.4f}", r["nst"], f"{r['elapsed_s']:.2f}"])
        _, Cp = read_table(f"{pref}_conc.csv")
        _, Ip = read_table(f"{pref}_irr.csv")
        Cp, Ip = np.array(Cp, float), np.array(Ip, float)
        if len(Cp) != len(Cb) or not np.allclose(Cp[:, 0], Cb[:, 0]):
            raise SystemExit(f"✗ {group} 与基线 h 列不一致")
        log(f"  [{group}] {len(sp_list)} 物种 × {scale:.3f} → nst={r['nst']}")
        return _calc_rows(group, sp_list, chdr, Cb, Ib, Cp, Ip, m, dtout, lo, hi,
                          tmin_h=spinup * 24.0)

    grp_rows, sp_rows = [], []
    for g in names:
        grp_rows += one(g, GROUPS[g])
    if args.per_species:
        for sp in PER_SPECIES:
            sp_rows += one(sp, [sp])
    if not grp_rows and not sp_rows:
        raise SystemExit("✗ RIR 无可计算样本")

    ndays_post = len({win_day(h, dtout) for h in Ib[:, 0]
                      if win_hours(h, dtout)[0] >= spinup * 24.0 - 1e-9})
    note = (f"个例分析·RIR | 平均RIR Eq(8)（P(O3) 加权）| -{args.pct:g}% 排放扰动 | "
            f"window={args.window} | P=G-D(Eq4-6) | S=模拟浓度 | dtout={dtout:g}s | "
            f"post-spinup {ndays_post} 天（n_days=每日加权样本数）")
    if ndays_post < 2:
        log(f"  ⚠ post-spinup 仅 {ndays_post} 天 → Eq(8) 的逐时加权退化：每个 hour 只有 "
            f"1 个样本，RIR_avg 就是该小时的 Eq(7) 单点值。要多日加权，需 "
            f"post-spinup ≥ 2 天（延长 model.parameters 天数或减小 SPINUP），"
            f"且气象/观测序列必须覆盖到新天数")

    paths = []
    if grp_rows:
        dd, ww = rir_tables(grp_rows, names)
        paths.append(write_csv(
            os.path.join(outdir, "rir_diurnal.csv"),
            ["group", "hour", "RIR_avg", "P_weight_sum", "n_days",
             "n_samples", "s_basis"], dd, note=note))
        paths.append(write_csv(
            os.path.join(outdir, "rir_mean.csv"),
            ["group", "RIR_mean", "P_weight_sum", "n_samples", "s_basis"],
            ww, note=note))
    if sp_rows:
        dd, ww = rir_tables(sp_rows, PER_SPECIES)
        note_sp = note + " | 逐物种：一次只扰一个排放物种（与组扰动不同池，勿混比）"
        paths.append(write_csv(
            os.path.join(outdir, "rir_species_diurnal.csv"),
            ["species", "hour", "RIR_avg", "P_weight_sum", "n_days",
             "n_samples", "s_basis"], dd, note=note_sp))
        paths.append(write_csv(
            os.path.join(outdir, "rir_species_mean.csv"),
            ["species", "RIR_mean", "P_weight_sum", "n_samples", "s_basis"],
            ww, note=note_sp))

    # 基线 O3 收支（同样剔除跨 spinup 边界的首窗）
    f = 3600.0 / (dtout * PPB)
    brow = []
    for t in range(len(Ib)):
        hstep = Ib[t, 0]
        if win_hours(hstep, dtout)[0] < spinup * 24.0 - 1e-9:
            continue
        d = m.p_o3_gd(Ib[t, 1:])
        brow.append([hstep, win_day(hstep, dtout), float(win_hour(hstep, dtout)),
                     *[d[k] * f for k in BUDGET_KEYS], "ppb/h"])
    paths.append(write_csv(
        os.path.join(outdir, "p_o3no_baseline.csv"),
        ["h", "day", "hour", *BUDGET_KEYS, "unit"], brow,
        note=f"个例分析·O3收支 Eq(4)-(6) | unit=ppb/h | h=输出步序号, "
             f"hour=本窗起点当地小时 0-23（dtout={dtout:g}s）| "
             f"post-spinup 完整窗（跨边界首窗已剔除）| "
             f"全窗输出（RIR 的 window 仅作用于 RIR 加权）"))
    return paths


def main(argv):
    ap = argparse.ArgumentParser(description="个例案例分析（六账 + IRR + RIR）")
    ap.add_argument("--scen", default="SCEN_BAU", help="情景名（默认 SCEN_BAU）")
    ap.add_argument("--species", default="O3", help="IRR 归因目标物种（默认 O3）")
    ap.add_argument("--pct", type=float, default=10.0, help="RIR 前体物减排百分比")
    ap.add_argument("--groups", default="all", help="all 或逗号分隔组名")
    ap.add_argument("--window", default="all", help="RIR 当地小时窗（all 或 lo-hi；all=0-23 全逐时，夜间 P≈0 自动剔除）")
    ap.add_argument("--per-species", action="store_true", help="附加逐物种 RIR")
    ap.add_argument("--no-rir", action="store_true", help="跳过 RIR（快速只做六账/IRR）")
    ap.add_argument("--no-assim", action="store_true",
                    help="关闭观测同化（忽略 config_override/obs/，自由运行基线）")
    a = ap.parse_args(argv[1:])

    global LOG_F
    os.makedirs(OUT, exist_ok=True)
    LOG_F = open(os.path.join(OUT, "run.log"), "w", encoding="utf-8")
    t0 = time.time()
    log(f"═══ 个例案例分析（情景 {a.scen}；六账 + IRR({a.species})"
        + ("；RIR 跳过" if a.no_rir else f" + RIR(-{a.pct:g}%)") + "）═══")

    cfg = prepare_case(os.path.join(OUT, "in_case"))
    constrained, obs_report = ([], []) if a.no_assim else prepare_observations(cfg)
    for ln in obs_report:
        log(ln)
    if constrained:
        log("  [同化] 观测约束物种 (" + str(len(constrained)) + "): " + ", ".join(constrained))
    rt = compile_runtime(cfg, os.path.join(OUT, "runtime_case.dat"))
    days, dtout, spinup = runtime_window(rt)
    log(f"  案例: 天数 {days:g} / spinup {spinup:g} / dtout {dtout:g}s")

    pref = os.path.join(OUT, "case")
    rb = run_fac(a.scen, pref, rt)
    log(f"  基线: {rb['line']}")

    m = Mechanism(MECH_PATH)
    p_ipr, species, nh = export_ipr(pref, OUT, dtout)
    p1, p2 = export_irr(pref, OUT, a.species, m, dtout)
    log(f"  六账: {len(species)} 物种 × {nh} 步; IRR top 已导出")
    paths = [p_ipr, p1, p2]
    runs = [["base", "-", "1.0", rb["nst"], f"{rb['elapsed_s']:.2f}"]]
    if not a.no_rir:
        paths += export_rir(pref, OUT, rt, a, runs)
    p_runs = write_csv(os.path.join(OUT, "runs.csv"),
                       ["case", "n_species", "scale", "nst", "elapsed_s"], runs,
                       note="个例分析 | 基线与 RIR 扰动运行记录")
    paths.append(p_runs)

    md5 = hashlib.md5(open(MECH_PATH, "rb").read()).hexdigest()[:12]
    rir_desc = ("跳过" if a.no_rir else
                f"-{a.pct:g}% window={a.window} groups={a.groups}"
                + (" +per-species" if a.per_species else ""))
    with open(os.path.join(OUT, "run_meta.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"个例案例分析 运行元数据\n"
                 f" 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                 f" 情景: {a.scen} | 天数 {days:g} / spinup {spinup:g} / dtout {dtout:g}s\n"
                 f" 机制 md5={md5}\n"
                 f" 输入: {os.path.relpath(SRC, ROOT)}（完整配置集，本目录即唯一输入源）\n"
                 f" 同化约束: " + (str(len(constrained)) + " 物种: " + ", ".join(constrained)
                 if constrained else "关闭（--no-assim 或无 obs/）") + "\n"
                 f" IRR 物种: {a.species} | RIR: {rir_desc}\n"
                 f" 耗时: {time.time() - t0:.1f} s\n")
    log(f"\n产物: {os.path.relpath(OUT, ROOT)}/（{len(paths)} 个 CSV + run.log/meta）")
    LOG_F.close()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
