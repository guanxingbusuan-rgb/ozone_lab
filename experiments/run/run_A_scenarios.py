#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_A_scenarios.py —— 实验 A：单箱基线与六情景政策沙盘（只跑 + 落 CSV，无绘图）

输入假设: runtime_v3.dat 已编译完成（本脚本不调 cfgbuild）
输出:
  output/run_A/<SCEN>_conc.csv     逐时 74 物种（末 24 h）
  output/run_A/<SCEN>_irr.csv      131 反应小时窗积分
  output/run_A/<SCEN>_ipr.csv      六账小时窗积分
  output/run_A/summary.csv         情景汇总（O3max/末24h均值/关键物种/nst/耗时）

用法:
  python experiments/run/run_A_scenarios.py                # 六情景全跑
  python experiments/run/run_A_scenarios.py BAU LNNOX      # 指定情景
  OZONE_OUT=/data/out python .../run_A_scenarios.py        # 换输出根目录
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (ROOT, require_prereq, out_dir, run_fac_v3, last24_mean,
                     write_csv, manifest, summary_table, log, SCEN_ALL)


def main(argv):
    scen = [s.upper().replace("SCEN_", "") for s in argv[1:]] or SCEN_ALL
    log(f"═══ 实验 A：单箱六情景（{len(scen)} 个）═══")
    require_prereq()
    d = out_dir("run_A")

    rows, products = [], []
    for s in scen:
        r = run_fac_v3(f"SCEN_{s}", os.path.join(d, s))
        pref = r["prefix"]
        o3_24 = last24_mean(f"{pref}_conc.csv", "O3")
        hcho24 = last24_mean(f"{pref}_conc.csv", "HCHO")
        no2_24 = last24_mean(f"{pref}_conc.csv", "NO2")
        no24 = last24_mean(f"{pref}_conc.csv", "NO")
        rows.append([s, round(r["o3max"], 2) if r["o3max"] == r["o3max"] else "nan",
                     round(o3_24, 2), round(hcho24, 3), round(no2_24, 3),
                     round(no24, 2), r["nst"], round(r["elapsed_s"], 2)])
        products += [f"{pref}_conc.csv", f"{pref}_irr.csv", f"{pref}_ipr.csv"]
        log(f"    {s:6s} O3max={r['o3max']:7.2f} ppb  末24h O3={o3_24:6.2f}  "
            f"FNR≈{hcho24 / max(no2_24, 1e-9):6.2f}  ({r['elapsed_s']:.2f}s)")

    summary_table("实验 A · 六情景汇总（ppb）",
                  ["scenario", "O3max", "O3_24h", "HCHO_24h", "NO2_24h", "NO_24h", "nst", "t_s"],
                  rows)
    p = write_csv(os.path.join(d, "summary.csv"),
                  ["scenario", "O3max_ppb", "O3_last24h_ppb", "HCHO_last24h_ppb",
                   "NO2_last24h_ppb", "NO_last24h_ppb", "nst", "elapsed_s"], rows,
                  note="实验A 单箱六情景汇总 | O3max 来自 fac_v3 stdout，末24h均值来自 conc CSV")
    products.append(p)
    mp = manifest(d, "A_scenarios", {"scenarios": scen}, products)
    log(f"\n产物: {os.path.relpath(d, ROOT)}/  ({len(products)} 个 CSV + manifest.json)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
