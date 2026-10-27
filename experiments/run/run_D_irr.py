#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_D_irr.py —— 实验 D（IRR 侧）：反应级归因（物种参数化，无绘图）

只读 <pref>_irr.csv + 计量表 gen_v3/oz.reac|oz.prod + 机制 .fac（物种序/反应式）
—— **不依赖** _ipr.csv / _conc.csv（过程账本见 run_D_ipr.py；互证见 run_G_consist.py）

对目标物种 SP，净贡献 = Σ_k coef_k(SP)·I_k，其中 coef 由计量表计算
（产物取 +系数，反应物取 −1），I_k 为该小时窗 IRR 积分。**窗口积分直接求和**
（每行已是本窗口积分，trapz 会二次积分）。

输出:
  output/run_D_irr/irr_top.csv          全程累计净贡献 top-12 反应
  output/run_D_irr/hourly_rxn_top.csv   逐小时净贡献 top-4 反应

用法:
  python experiments/run/run_D_irr.py                     # O3，现跑 BAU
  python experiments/run/run_D_irr.py HCHO                # 指定物种
  python experiments/run/run_D_irr.py O3 output/run_A/BAU # 指定物种 + 复用前缀
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _runner import (ROOT, require_prereq, out_dir, run_fac_v3, read_table,
                     write_csv, summary_table, log)

sys.path.insert(0, os.path.join(ROOT, "analysis_lib"))
from indicators import Mechanism

TOP_CUM, TOP_HR = 12, 4


def main(argv):
    log("═══ 实验 D·IRR：反应级归因 ═══")
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
        tmpdir = None
    else:
        # 自跑 BAU：原始模拟输出落系统临时目录（读完即删），
        # 保证 output/run_D_irr/ 只留分析 CSV 产物
        import tempfile
        tmpdir = tempfile.mkdtemp(prefix="oz_D_irr_")
        pref = os.path.join(tmpdir, "BAU")
        run_fac_v3("SCEN_BAU", pref, verbose=True)
    d = out_dir("run_D_irr")
    log(f"  目标物种: {species} (gen #{jSP})")

    # ── 计量系数: 每反应对 SP 的净计量（产 +coef / 耗 −1）──
    gen = os.path.join(ROOT, "facsimile_impl", "gen_v3")
    reac = [ln.split() for ln in open(os.path.join(gen, "oz.reac"))]
    prod = [ln.split() for ln in open(os.path.join(gen, "oz.prod"))]
    mech_lines = [l.strip().lstrip("%").strip()
                  for l in open(os.path.join(ROOT, "mechanism",
                                             "mechanism_3_full_chain.fac"))
                  if l.startswith("%")]
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

    # ── IRR 逐窗口读取（首列 h，其后为反应列）──
    hi, DI = read_table(f"{pref}_irr.csv")
    I = np.array([r[1:1 + len(reac)] for r in DI], dtype=float)
    contrib_w = I * coef                                  # (T, nr) 逐窗口净贡献

    # ── ① 全程累计 top-N（窗口积分直接求和，不做 trapz）──
    contrib_sum = contrib_w.sum(axis=0)
    order = np.argsort(-np.abs(contrib_sum))[:TOP_CUM]
    rows = []
    for k in order:
        eq = mech_lines[k] if k < len(mech_lines) else "?"
        rows.append([k + 1, f"{contrib_sum[k]:.5e}", eq[:70]])
    summary_table(f"实验 D·IRR · {species} 反应级 top-{TOP_CUM}"
                  f"（净贡献 = 计量×Σ窗口IRR）",
                  ["r#", f"net_{species}", "reaction"], rows)
    p1 = write_csv(os.path.join(d, "irr_top.csv"),
                   ["r_index", f"net_{species}_contribution", "reaction"], rows,
                   note=(f"实验D·IRR | {species} 净贡献 top-{TOP_CUM} | "
                         "正值=净生成, 负值=净消耗 | 窗口积分求和（勿 trapz）"))

    # ── ② 逐小时 top-N ──
    hours = [int(r[0]) for r in DI]
    rows_h = []
    for i, hh in enumerate(hours):
        ordh = np.argsort(-np.abs(contrib_w[i]))[:TOP_HR]
        for rank, k in enumerate(ordh, 1):
            rows_h.append([hh, rank, k + 1, f"{contrib_w[i][k]:.4e}",
                           mech_lines[k][:64] if k < len(mech_lines) else "?"])
    p2 = write_csv(os.path.join(d, "hourly_rxn_top.csv"),
                   ["h", "rank", "r_index", f"net_{species}", "reaction"], rows_h,
                   note=f"实验D·IRR | 每小时 {species} 净贡献 top-{TOP_HR} 反应 | 正=生成 负=消耗")

    mp_log = [os.path.relpath(p1, ROOT), os.path.relpath(p2, ROOT)]
    log(f"\n产物: {mp_log[0]}\n      {mp_log[1]}")
    if tmpdir:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
