#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_G_consist.py —— 实验 G：IPR↔IRR 双路径一致性互证（物种参数化，无绘图）

**这是唯一同时读 _ipr.csv 与 _irr.csv 的脚本**（互证的定义使然）：
  路径 A(IPR): <pref>_ipr.csv 的 {SP}_P / {SP}_L 列（fac_v3 求解器内累积）
  路径 B(IRR): gen_v3 计量表 Σ coef_k·I_k（由 _irr.csv 独立重算）
两路径来自独立累积，应逐位一致（rel < 1e-6）。

★ 关键: IPR_L 已是**负值**（机制耗），而 IRR 侧算出的"耗"也是负值，
  直接相加会得到 2 倍关系（曾误判为 2.027 倍不一致）。
  正确做法是**分别**比对产（P）与耗（L）两条路径。

输出:
  output/run_G_consist/consistency.csv   IPR 与 IRR 对目标物种的逐时对比 + 最大相对差

用法:
  python experiments/run/run_G_consist.py                     # O3，现跑 BAU
  python experiments/run/run_G_consist.py HCHO                # 指定物种
  python experiments/run/run_G_consist.py O3 output/run_A/BAU # 指定物种 + 复用前缀
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _runner import (ROOT, require_prereq, out_dir, run_fac_v3, read_table,
                     write_csv, manifest, log)

sys.path.insert(0, os.path.join(ROOT, "analysis_lib"))
from indicators import Mechanism


def main(argv):
    log("═══ 实验 G：IPR↔IRR 双路径互证 ═══")
    require_prereq()

    # ── 参数解析: [物种] [前缀] 顺序任意 ──
    mech = Mechanism(os.path.join(ROOT, "mechanism", "mechanism_3_full_chain.fac"))
    species, pref_arg = None, None
    for a in argv[1:]:
        if os.path.exists(f"{a}_conc.csv"):
            pref_arg = a
        elif a.upper() in mech.species:
            species = a.upper()
        else:
            log(f"✗ 无法识别参数 '{a}'（既非已有输出前缀，也非机制物种）")
            return 1
    species = species or "O3"
    jSP = mech.species.index(species) + 1                # gen 计量表 1-based 物种号

    if pref_arg:
        pref = pref_arg
        log(f"  复用已有输出: {pref}")
        src_tag = "reused"
    else:
        d = out_dir("run_G_consist")
        pref = os.path.join(d, "BAU")
        r = run_fac_v3("SCEN_BAU", pref, verbose=True)
        src_tag = "fresh"
    d = out_dir("run_G_consist")
    log(f"  目标物种: {species} (gen #{jSP})")

    # ── 计量系数（与 run_D_irr 同法）: 产 +coef / 耗 −1 ──
    gen = os.path.join(ROOT, "facsimile_impl", "gen_v3")
    reac = [ln.split() for ln in open(os.path.join(gen, "oz.reac"))]
    prod = [ln.split() for ln in open(os.path.join(gen, "oz.prod"))]
    coef = np.zeros(len(reac))
    for k in range(len(reac)):
        c = 0.0
        pr = prod[k]
        for t in range(int(pr[0])):
            if int(pr[1 + 2 * t]) == jSP:
                c += float(pr[2 + 2 * t])
        rr = reac[k]
        for t in range(int(rr[0])):
            if int(rr[1 + t]) == jSP:
                c -= 1.0
        coef[k] = c

    # ── 读两份 CSV（互证的两条独立路径）──
    hp, DP = read_table(f"{pref}_ipr.csv")
    hi, DI = read_table(f"{pref}_irr.csv")
    if f"{species}_P" not in hp:
        log(f"✗ IPR 表头无 {species}_P 列（检查机制是否含该物种）")
        return 1

    # 分产/耗两路: 只取产 SP 的反应 / 只取耗 SP 的反应
    coefP = np.maximum(coef, 0.0)          # 只取产 SP 的反应
    coefL = np.minimum(coef, 0.0)          # 只取耗 SP 的反应
    rows_c, worst = [], 0.0
    for i in range(min(len(DP), len(DI))):
        I_row = np.array(DI[i][1:1 + len(reac)])
        ipr_P = DP[i][hp.index(f"{species}_P")]
        ipr_L = DP[i][hp.index(f"{species}_L")]
        irr_P = float(np.dot(coefP, I_row))
        irr_L = float(np.dot(coefL, I_row))     # 负值
        relP = abs(ipr_P - irr_P) / max(abs(ipr_P), abs(irr_P), 1e-30)
        relL = abs(ipr_L - irr_L) / max(abs(ipr_L), abs(irr_L), 1e-30)
        worst = max(worst, relP, relL)
        rows_c.append([int(DP[i][0]), f"{ipr_P:.6e}", f"{irr_P:.6e}", f"{relP:.3e}",
                       f"{ipr_L:.6e}", f"{irr_L:.6e}", f"{relL:.3e}"])
    p1 = write_csv(os.path.join(d, "consistency.csv"),
                   ["h", "IPR_P", "IRR_P", "rel_P", "IPR_L", "IRR_L", "rel_L"], rows_c,
                   note=(f"实验G | {species} 产/耗双路核对: IPR 六账 vs IRR 计量重算 | "
                         "两独立累积路径应逐位一致(rel<1e-6)"))

    log(f"\n  IPR↔IRR 最大相对差 = {worst * 100:.6f}%  "
        f"{'✓ 两累积路径一致' if worst < 1e-5 else '⚠ 检查对位'}")
    mp = manifest(d, "G_consist",
                  {"source": src_tag, "prefix": pref, "species": species},
                  [p1, f"{pref}_conc.csv", f"{pref}_irr.csv", f"{pref}_ipr.csv"])
    log(f"\n产物: {os.path.relpath(d, ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
