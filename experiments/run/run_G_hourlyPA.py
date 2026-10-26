#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_G_hourlyPA.py —— 实验 G：逐小时过程分析（CMAq-PA 式，无绘图）

对一个 fac_v3 前缀，逐小时给出：
  ① 六账占比（P/L/E/D/B/X 归一化 %）——回答"这一小时谁主导"
  ② O3 的反应级 top 通道（每小时取净贡献 top-4）
  ③ IPR↔IRR 双路径一致性（两组独立累积互证）

输出:
  output/run_G/hourly_ledger.csv     h × species × 六账（原始积分）
  output/run_G/hourly_share.csv      h × species × 六账占比（%）
  output/run_G/hourly_rxn_top.csv    h × top4 反应（净 O3 贡献）
  output/run_G/consistency.csv       IPR 与 IRR 对 O3 的逐时对比 + 最大相对差

用法:
  python experiments/run/run_G_hourlyPA.py                    # 现跑 BAU
  python experiments/run/run_G_hourlyPA.py output/run_A/BAU   # 复用
"""
import os
import sys
import numpy as np
import csv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (ROOT, require_prereq, out_dir, run_fac_v3, read_table, write_csv,
                     manifest, summary_table, log, col, _fnum)

FOCUS = ["O3", "HNO3", "H2O2", "CO", "PAN", "NO2"]
PROCS = list("PLEDBX")


def main(argv):
    log("═══ 实验 G：逐小时过程分析 ═══")
    require_prereq()
    d = out_dir("run_G")

    if len(argv) > 1 and os.path.exists(f"{argv[1]}_conc.csv"):
        pref = argv[1]; log(f"  复用: {pref}")
    else:
        pref = os.path.join(d, "BAU")
        run_fac_v3("SCEN_BAU", pref, verbose=True)

    # ── 读三份 CSV ──
    hp, DP = read_table(f"{pref}_ipr.csv")           # 首列 h
    hi, DI = read_table(f"{pref}_irr.csv")
    hc, DC = read_table(f"{pref}_conc.csv")

    hours = [int(r[0]) for r in DP]
    n = len(hours)

    # 物种顺序（从 IPR 表头每 6 列一个）
    sp_names = []
    j = 1
    while j < len(hp) and hp[j].split("_")[0] not in ("fnox",):
        sp_names.append(hp[j].split("_")[0]); j += 6
    # 机制物种序（用于 conc 列对应）
    sys.path.insert(0, os.path.join(ROOT, "diagnostics"))
    from indicators import Mechanism
    mech = Mechanism(os.path.join(ROOT, "mechanism", "mechanism_3_full_chain.fac"))
    names = mech.species

    # ── ① 逐小时六账（原始 + 占比）──
    rows_raw, rows_share = [], []
    for i, hh in enumerate(hours):
        for sp in FOCUS:
            if sp not in sp_names:
                continue
            cols = [hp.index(f"{sp}_{p}") for p in PROCS]
            v = [DP[i][c] for c in cols]
            tot = sum(abs(x) for x in v) or 1.0
            rows_raw.append([hh, sp] + [f"{x:.6e}" for x in v])
            rows_share.append([hh, sp] + [f"{100.0 * x / tot:.2f}" for x in v])
    p1 = write_csv(os.path.join(d, "hourly_ledger.csv"),
                   ["h", "species"] + PROCS, rows_raw,
                   note="实验G | 逐小时六账积分 [molec/cm3] | P机制产 L机制耗 E排放 D干湿沉 B同化+稀释 X BLH膨胀")
    p2 = write_csv(os.path.join(d, "hourly_share.csv"),
                   ["h", "species"] + [f"{p}_pct" for p in PROCS], rows_share,
                   note="实验G | 六账占比 (%) | 每行按 |六账| 之和归一")

    # ── ② 逐小时反应级 top-4（对 O3 的净贡献）──
    gen = os.path.join(ROOT, "facsimile_impl", "gen_v3")
    jO3 = names.index("O3") + 1
    reac = [ln.split() for ln in open(f"{gen}/oz.reac")]
    prod = [ln.split() for ln in open(f"{gen}/oz.prod")]
    mech_lines = [l.strip().lstrip("%").strip()
                  for l in open(os.path.join(ROOT, "mechanism", "mechanism_3_full_chain.fac"))
                  if l.startswith("%")]
    coef = np.zeros(len(reac))
    for k in range(len(reac)):
        c = 0.0
        pr = prod[k]; npr = int(pr[0])
        for t in range(npr):
            if int(pr[1 + 2 * t]) == jO3:
                c += float(pr[2 + 2 * t])
        rr = reac[k]; nr = int(rr[0])
        for t in range(nr):
            if int(rr[1 + t]) == jO3:
                c -= 1.0
        coef[k] = c

    rows_top = []
    for i, hh in enumerate(hours):
        contrib = coef * np.array(DI[i][1:1 + len(reac)])
        order = np.argsort(-np.abs(contrib))[:4]
        for rank, k in enumerate(order, 1):
            rows_top.append([hh, rank, k + 1, f"{contrib[k]:.4e}",
                             mech_lines[k][:64] if k < len(mech_lines) else "?"])
    p3 = write_csv(os.path.join(d, "hourly_rxn_top.csv"),
                   ["h", "rank", "r_index", "net_O3", "reaction"], rows_top,
                   note="实验G | 每小时 O3 净贡献 top-4 反应 | 正=生成 负=消耗")

    # ── ③ IPR↔IRR 双路径一致性（对 O3，分产/耗两路核对）──
    #
    # ★ 关键: IPR_L 已是**负值**（机制耗），而 IRR 侧算出的"耗"若是正数，
    #   直接相加会得到 2 倍关系（曾误判为 2.027 倍不一致）。
    #   正确做法是**分别**比对产（P）与耗（L）两条路径。
    #
    # 路径 A(IPR): O3_P / O3_L 两列
    # 路径 B(IRR): 由 gen_v3 计量表 Σ coef_k·I_k，正向(产)/负向(耗) 分别累加
    coefP = np.maximum(coef, 0.0)          # 只取产 O3 的反应
    coefL = np.minimum(coef, 0.0)          # 只取耗 O3 的反应
    rows_c, worst = [], 0.0
    for i, hh in enumerate(hours):
        I_row = np.array(DI[i][1:1 + len(reac)])
        ipr_P = DP[i][hp.index("O3_P")]
        ipr_L = DP[i][hp.index("O3_L")]
        irr_P = float(np.dot(coefP, I_row))
        irr_L = float(np.dot(coefL, I_row))     # 负值
        relP = abs(ipr_P - irr_P) / max(abs(ipr_P), abs(irr_P), 1e-30)
        relL = abs(ipr_L - irr_L) / max(abs(ipr_L), abs(irr_L), 1e-30)
        worst = max(worst, relP, relL)
        rows_c.append([hh, f"{ipr_P:.6e}", f"{irr_P:.6e}", f"{relP:.3e}",
                       f"{ipr_L:.6e}", f"{irr_L:.6e}", f"{relL:.3e}"])
    p4 = write_csv(os.path.join(d, "consistency.csv"),
                   ["h", "IPR_P", "IRR_P", "rel_P", "IPR_L", "IRR_L", "rel_L"], rows_c,
                   note=("实验G | O3 产/耗双路核对: IPR 六账 vs IRR 计量重算 | "
                         "两独立累积路径应逐位一致(rel<1e-6)"))

    log(f"\n  IPR↔IRR 最大相对差 = {worst * 100:.6f}%  "
        f"{'✓ 两累积路径一致' if worst < 1e-5 else '⚠ 检查对位'}")

    # 终端摘要：O3 逐小时占比（每 3 小时取样）
    sample = [i for i, hh in enumerate(hours) if hh % 3 == 0]
    summary_table("实验 G · O₃ 六账占比 (%) — 每 3 小时",
                  ["h", "P", "L", "E", "D", "B", "X"],
                  [[hours[i]] + [next(float(r[2 + k]) for r in rows_share
                                      if int(r[0]) == hours[i] and r[1] == "O3")
                                 for k in range(6)] for i in sample])

    mp = manifest(d, "G_hourlyPA", {"prefix": pref, "hours": [hours[0], hours[-1]]},
                  [p1, p2, p3, p4, f"{pref}_conc.csv", f"{pref}_irr.csv", f"{pref}_ipr.csv"])
    log(f"\n产物: {os.path.relpath(d, ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
