#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_E_indicators.py —— 实验 E：指示因子工具箱（12+ 指标，任意路线通用，无绘图）

对指定前缀（fac_v3 / batch / mcm 输出）批量计算指示因子，落 CSV。
**不做反演**（反演见 run_H）。

输出:
  output/run_E/<tag>_indicators.csv    窗口均值指标（12 个）
  output/run_E/<tag>_timeseries.csv    逐时指标序列（24 行 × 12 列）
  output/run_E/compare.csv             多前缀横向对比（若给了多个）

用法:
  python experiments/run/run_E_indicators.py output/run_A/BAU
  python experiments/run/run_E_indicators.py output/run_A/BAU output/run_A/LNNOX
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (ROOT, require_prereq, out_dir, write_csv, manifest,
                     summary_table, log)

sys.path.insert(0, os.path.join(ROOT, "diagnostics"))
from indicators import Mechanism, route_indicators, route_timeseries

MECH = os.path.join(ROOT, "mechanism", "mechanism_3_full_chain.fac")


def analyze(pref, tag, d, mech):
    log(f"  [{tag}] {pref}")
    win = route_indicators(pref, MECH, mech=mech)            # 末 24h 窗口均值
    ts = route_timeseries(pref, MECH, mech=mech)             # 逐时序列
    keys = list(win.keys())
    p1 = write_csv(os.path.join(d, f"{tag}_indicators.csv"),
                   ["indicator", "value"], [[k, f"{win[k]:.6g}"] for k in keys],
                   note=f"实验E | {tag} | 末24h窗口均值 | 12 指标")
    n = len(ts.get("FNR", []))
    rows = []
    for i in range(n):
        rows.append([96 + i] + [f"{ts[k][i]:.6g}" for k in keys])
    p2 = write_csv(os.path.join(d, f"{tag}_timeseries.csv"), ["h"] + keys, rows,
                   note=f"实验E | {tag} | 逐时指标（末24h）")
    return win, keys, [p1, p2]


def main(argv):
    prefs = argv[1:] or [os.path.join("output", "run_D", "BAU")]
    log(f"═══ 实验 E：指示因子（{len(prefs)} 个前缀）═══")
    require_prereq()
    d = out_dir("run_E")
    mech = Mechanism(MECH)

    results, all_products = {}, []
    for p in prefs:
        base = os.path.abspath(p)
        if not os.path.exists(f"{base}_conc.csv"):
            log(f"  ✗ 跳过（无 {base}_conc.csv）: {p}")
            continue
        tag = os.path.basename(base.rstrip("/"))
        win, keys, prods = analyze(base, tag, d, mech)
        results[tag] = win
        all_products += prods

    if not results:
        log("✗ 无可分析前缀")
        return 1

    # 横向对比
    keys = list(next(iter(results.values())).keys())
    rows = [[k] + [f"{results[t][k]:.6g}" for t in results] for k in keys]
    p = write_csv(os.path.join(d, "compare.csv"), ["indicator"] + list(results), rows,
                  note="实验E | 多前缀指标对比 | 末24h均值")
    all_products.append(p)

    summary_table("实验 E · 指示因子对比",
                  ["indicator"] + list(results),
                  [[k] + [round(results[t][k], 3) for t in results] for k in keys])

    log("\n  判读: FNR<1→VOC-limited; >2→NOx-limited; L_NOx/Q>0.5→NOx 饱和;"
        " OPE 低→NOx 饱和")
    mp = manifest(d, "E_indicators", {"prefixes": prefs}, all_products)
    log(f"\n产物: {os.path.relpath(d, ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
