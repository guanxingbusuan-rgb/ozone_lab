#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_D_ipr.py —— 实验 D（IPR 侧）：六账过程账本，逐小时 × 全物种（无绘图）

只读 <pref>_ipr.csv —— **单一输入依赖**（不碰 _conc.csv / _irr.csv / 计量表 /
机制文件）。反应级归因见 run_D_irr.py；IPR↔IRR 互证见 run_G_consist.py；
逐窗口闭合校验方法见旧版 experiments/ipr_closure_check.py。

输出:
  output/run_D_ipr/ipr_ledger.csv   长表: h × 物种 × 六账（P L E D B X，小时窗积分）
                                    覆盖机制全部物种；占比 % 可由 |六账| 归一化重算

用法:
  python experiments/run/run_D_ipr.py                     # 现跑 BAU
  python experiments/run/run_D_ipr.py output/run_A/BAU    # 复用已有输出前缀
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _runner import (ROOT, require_prereq, out_dir, run_fac_v3, read_table,
                     write_csv, summary_table, log)

PROCS = list("PLEDBX")


def main(argv):
    log("═══ 实验 D·IPR：六账过程账本（逐小时 × 全物种）═══")
    require_prereq()
    d = out_dir("run_D_ipr")

    if len(argv) > 1 and os.path.exists(f"{argv[1]}_ipr.csv"):
        pref = argv[1]
        log(f"  复用已有输出: {pref}")
        tmpdir = None
    else:
        # 自跑 BAU：原始模拟输出落系统临时目录（读完即删），
        # 保证 output/run_D_ipr/ 只留 ipr_ledger.csv 一个产物
        import tempfile
        tmpdir = tempfile.mkdtemp(prefix="oz_D_ipr_")
        pref = os.path.join(tmpdir, "BAU")
        run_fac_v3("SCEN_BAU", pref, verbose=True)

    # ── 读 _ipr.csv: 表头每 6 列一个物种，每行 = 本小时窗积分 ──
    hp, DP = read_table(f"{pref}_ipr.csv")
    species, j = [], 1
    while j < len(hp) and hp[j].split("_")[0] not in ("fnox",):
        species.append(hp[j].split("_")[0])
        j += 6

    # ── 长表: h, species, P, L, E, D, B, X ──
    rows = []
    for i, r in enumerate(DP):
        hh = int(r[0])
        for si, sp in enumerate(species):
            base = 1 + 6 * si
            rows.append([hh, sp] + [f"{r[base + p]:.6e}" for p in range(6)])

    p1 = write_csv(os.path.join(d, "ipr_ledger.csv"),
                   ["h", "species", "P", "L", "E", "D", "B", "X"], rows,
                   note=("实验D·IPR | 六账过程账本（长表，逐小时 × 全物种）| "
                         "P机制产 L机制耗 E排放 D干湿沉 B同化+稀释 X BLH膨胀 | "
                         "单位 molec/cm3，每行=本小时窗积分（勿再 trapz）| "
                         "占比% = 该行账目/Σ|六账|，全程总量 = 按物种对列求和"))

    # 终端摘要：O3 逐小时占比（每 3 小时取样）
    if "O3" in species:
        si = species.index("O3")
        sample = []
        for r in DP:
            if int(r[0]) % 3 == 0:
                v = [r[1 + 6 * si + p] for p in range(6)]
                tot = sum(abs(x) for x in v) or 1.0
                sample.append([int(r[0])] + [f"{100.0 * x / tot:.2f}" for x in v])
        summary_table("实验 D·IPR · O₃ 六账占比 (%) — 每 3 小时",
                      ["h", "P", "L", "E", "D", "B", "X"], sample)

    log(f"\n产物: {os.path.relpath(p1, ROOT)}  "
        f"({len(species)} 物种 × {len(DP)} 小时 = {len(rows)} 行)")
    if tmpdir:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
