#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CASE-analysis/run.py —— 个例案例分析（单一入口：模拟 + 全部 CSV 存档）

把原 experiments/run 下的 D1（六账账本）/ D2（反应级归因）/ D3（前体物 RIR）
合并为一个自完备程序：一个「个例」= 本目录 config_override/ 的一份完整配置快照。

自包含: 不 import experiments/run；共享仅"模型本体"四件套
        （facsimile_impl/fac_v3、inputproc/cfgbuild.py、机制 .fac、
         analysis_lib/csvio 与 indicators）。
输入: config_override/ 完整配置集（缺省已由根 configs/ 种子化；直接改这里即换个例）
输出（本目录 output/）:
  in_case/ + runtime_case.dat         配置快照与编译产物（可追溯的第一事实源）
  case_{conc,ipr,irr}.csv             基线 fac_v3 原始输出
  ipr_ledger.csv                      六账长表（h × 物种 × P/L/E/D/B/X）
  irr_top.csv / hourly_rxn_top.csv    反应级净贡献（默认 O3）
  rir_mean.csv / rir_diurnal.csv      前体物 RIR（论文 Eq4–8，P(O3) 加权）
  p_o3no_baseline.csv                 基线 O3 收支（Eq4–6）
  runs.csv                            基线+RIR 扰动运行记录
  run.log / run_meta.txt              运行日志与元数据

用法:
  python CASE-analysis/run.py                       # 基线 + 六账 + IRR + RIR(-10%)
  python CASE-analysis/run.py --no-rir              # 只做基线/六账/IRR（快）
  python CASE-analysis/run.py --pct 20 --groups NOx,CO --window 8-18
  python CASE-analysis/run.py --species HCHO --per-species
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


# ── 分析核（原 D1/D2/D3 逻辑，去 experiments/run 依赖）──
def export_ipr(pref, outdir):
    hp, DP = read_table(f"{pref}_ipr.csv")
    species, j = [], 1
    while j < len(hp) and hp[j].split("_")[0] != "fnox":
        species.append(hp[j].split("_")[0])
        j += 6
    rows = []
    for r in DP:
        hh = int(r[0])
        for si, sp in enumerate(species):
            b = 1 + 6 * si
            rows.append([hh, sp] + [f"{r[b + q]:.6e}" for q in range(6)])
    p = write_csv(os.path.join(outdir, "ipr_ledger.csv"),
                  ["h", "species", "P", "L", "E", "D", "B", "X"], rows,
                  note="个例分析·六账 | P机制产 L机制耗 E排放 D干湿沉 B同化+稀释 X BLH膨胀 | "
                       "单位 molec/cm3，每行=本小时窗积分（勿再 trapz）")
    return p, species, len(DP)


def export_irr(pref, outdir, species, mech):
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
        for rank, k in enumerate(np.argsort(-np.abs(W[i]))[:4], 1):
            rh.append([int(r[0]), rank, k + 1, f"{W[i][k]:.4e}",
                       eqs[k][:64] if k < len(eqs) else "?"])
    p2 = write_csv(os.path.join(outdir, "hourly_rxn_top.csv"),
                   ["h", "rank", "r_index", f"net_{species}", "reaction"], rh,
                   note=f"个例分析·IRR | 每小时 {species} 净贡献 top-4")
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


def _calc_rows(group, species, chdr, Cb, Ib, Cp, Ip, m, dtout, lo, hi):
    ppbc = group not in ("NOx", "CO", "NO", "NO2")
    rows = []
    for t in range(len(Ib)):
        hour = (Ib[t, 0] * dtout / 3600.0) % 24.0
        if not (lo - 1e-9 <= hour <= hi + 1e-9):
            continue
        pb, pp = _p_ppbh(m, Ib[t, 1:], dtout), _p_ppbh(m, Ip[t, 1:], dtout)
        sb = _s_value(Cb[t], chdr, species, m, ppbc)
        sp = _s_value(Cp[t], chdr, species, m, ppbc)
        ds = sb - sp
        ok = (np.isfinite(pb) and np.isfinite(pp) and np.isfinite(sb)
              and np.isfinite(sp) and pb > 0 and sb > 0 and ds / sb > 1e-12)
        rows.append({"group": group, "day": int(Ib[t, 0] * dtout / 86400.0) + 1,
                     "hour": float(hour), "P_base": pb, "RIR":
                     ((pb - pp) / pb) / (ds / sb) if ok else np.nan,
                     "valid": bool(ok), "basis": "ppbC" if ppbc else "ppb"})
    return rows


def _weighted(rows):
    v = [r for r in rows if r["valid"]]
    w = sum(r["P_base"] for r in v)
    if not v or w <= 0:
        return np.nan, w, 0
    return sum(r["RIR"] * r["P_base"] for r in v) / w, w, len(v)


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
        return _calc_rows(group, sp_list, chdr, Cb, Ib, Cp, Ip, m, dtout, lo, hi)

    all_rows = []
    for g in names:
        all_rows += one(g, GROUPS[g])
    if args.per_species:
        for sp in PER_SPECIES:
            all_rows += one(sp, [sp])
    if not all_rows:
        raise SystemExit("✗ RIR 无可计算样本")

    note = (f"个例分析·RIR | 平均RIR Eq(8)（P(O3) 加权）| -{args.pct:g}% 排放扰动 | "
            f"window={args.window} | P=G-D(Eq4-6) | S=模拟浓度 | dtout={dtout:g}s")
    drows, wrows = [], []
    for g in names:
        for hh in range(24):
            sel = [r for r in all_rows if r["group"] == g and int(round(r["hour"])) % 24 == hh]
            rir, w, n = _weighted(sel)
            if n:
                drows.append([g, hh, f"{rir:.6e}", f"{w:.6e}",
                              len({r["day"] for r in sel if r["valid"]}), n,
                              "ppbC" if sel[0]["basis"] == "ppbC" else "ppb"])
        rir, w, n = _weighted([r for r in all_rows if r["group"] == g])
        if n:
            wrows.append([g, f"{rir:.6e}", f"{w:.6e}", n,
                          "ppbC" if GROUPS and g not in ("NOx", "CO") else "ppb"])
    p_d = write_csv(os.path.join(outdir, "rir_diurnal.csv"),
                    ["group", "hour", "RIR_avg", "P_weight_sum", "n_days",
                     "n_samples", "s_basis"], drows, note=note)
    p_m = write_csv(os.path.join(outdir, "rir_mean.csv"),
                    ["group", "RIR_mean", "P_weight_sum", "n_samples", "s_basis"],
                    wrows, note=note)
    # 基线 O3 收支
    f = 3600.0 / (dtout * PPB)
    brow = []
    for t in range(len(Ib)):
        hour = (Ib[t, 0] * dtout / 3600.0) % 24.0
        if lo - 1e-9 <= hour <= hi + 1e-9:
            d = m.p_o3_gd(Ib[t, 1:])
            brow.append([Ib[t, 0], int(Ib[t, 0] * dtout / 86400.0) + 1, hour,
                         *[d[k] * f for k in BUDGET_KEYS], "ppb/h"])
    p_b = write_csv(os.path.join(outdir, "p_o3no_baseline.csv"),
                    ["h", "day", "hour", *BUDGET_KEYS, "unit"], brow,
                    note=f"个例分析·O3收支 Eq(4)-(6) | unit=ppb/h | window={args.window}")
    return p_d, p_m, p_b


def main(argv):
    ap = argparse.ArgumentParser(description="个例案例分析（六账 + IRR + RIR）")
    ap.add_argument("--scen", default="SCEN_BAU", help="情景名（默认 SCEN_BAU）")
    ap.add_argument("--species", default="O3", help="IRR 归因目标物种（默认 O3）")
    ap.add_argument("--pct", type=float, default=10.0, help="RIR 前体物减排百分比")
    ap.add_argument("--groups", default="all", help="all 或逗号分隔组名")
    ap.add_argument("--window", default="8-18", help="RIR 当地小时窗（all 或 lo-hi）")
    ap.add_argument("--per-species", action="store_true", help="附加逐物种 RIR")
    ap.add_argument("--no-rir", action="store_true", help="跳过 RIR（快速只做六账/IRR）")
    a = ap.parse_args(argv[1:])

    global LOG_F
    os.makedirs(OUT, exist_ok=True)
    LOG_F = open(os.path.join(OUT, "run.log"), "w", encoding="utf-8")
    t0 = time.time()
    log(f"═══ 个例案例分析（情景 {a.scen}；六账 + IRR({a.species})"
        + ("；RIR 跳过" if a.no_rir else f" + RIR(-{a.pct:g}%)") + "）═══")

    cfg = prepare_case(os.path.join(OUT, "in_case"))
    rt = compile_runtime(cfg, os.path.join(OUT, "runtime_case.dat"))
    days, dtout, spinup = runtime_window(rt)
    log(f"  案例: 天数 {days:g} / spinup {spinup:g} / dtout {dtout:g}s")

    pref = os.path.join(OUT, "case")
    rb = run_fac(a.scen, pref, rt)
    log(f"  基线: {rb['line']}")

    m = Mechanism(MECH_PATH)
    p_ipr, species, nh = export_ipr(pref, OUT)
    p1, p2 = export_irr(pref, OUT, a.species, m)
    log(f"  六账: {len(species)} 物种 × {nh} 步; IRR top 已导出")
    paths = [p_ipr, p1, p2]
    runs = [["base", "-", "1.0", rb["nst"], f"{rb['elapsed_s']:.2f}"]]
    if not a.no_rir:
        paths += list(export_rir(pref, OUT, rt, a, runs))
    p_runs = write_csv(os.path.join(OUT, "runs.csv"),
                       ["case", "n_species", "scale", "nst", "elapsed_s"], runs,
                       note="个例分析 | 基线与 RIR 扰动运行记录")
    paths.append(p_runs)

    md5 = hashlib.md5(open(MECH_PATH, "rb").read()).hexdigest()[:12]
    with open(os.path.join(OUT, "run_meta.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"个例案例分析 运行元数据\n"
                 f" 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                 f" 情景: {a.scen} | 天数 {days:g} / spinup {spinup:g} / dtout {dtout:g}s\n"
                 f" 机制 md5={md5}\n"
                 f" 输入: {os.path.relpath(SRC, ROOT)}（完整配置集，不回退根 configs）\n"
                 f" IRR 物种: {a.species} | RIR: "
                 f"{'跳过' if a.no_rir else f'-{a.pct:g}% window={a.window} groups={a.groups}'}\n"
                 f" 耗时: {time.time() - t0:.1f} s\n")
    log(f"\n产物: {os.path.relpath(OUT, ROOT)}/（{len(paths)} 个 CSV + run.log/meta）")
    LOG_F.close()
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
