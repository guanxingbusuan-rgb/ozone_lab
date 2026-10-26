#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_C_ridge.py —— 实验 C：屋脊的跨路线/跨机制对比（无绘图）

对同一 fnox 序列，分别用 74sp 工程机制（fac_v3 / batch）与 29sp MCM 子集（mcm_run）
跑，比较各自 O3max 屋脊位置与指示因子穿越点。

输入假设: runtime_v3.dat 已编译；mcm_run 已构建
输出:
  output/run_C/ridge_batch.csv    batch 路线（从 run_B 或现跑）
  output/run_C/ridge_mcm.csv      MCM 路线（逐 fnox）
  output/run_C/ridge_compare.csv  两路线屋脊位置/一致性对比
  output/run_C/mcm_<fnox>_{conc,irr,ipr}.csv   MCM 逐点输出

用法:
  python experiments/run/run_C_ridge.py            # 13 点 logspace
  python experiments/run/run_C_ridge.py 9          # 9 点
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (ROOT, require_prereq, out_dir, run_batch, run_mcm, read_table,
                     write_csv, manifest, summary_table, log, col, last24_mean)

FN_DEFAULT = np.logspace(-1.0, 1.0, 13)


def main(argv):
    npt = int(argv[1]) if len(argv) > 1 else 13
    FN = np.logspace(-1.0, 1.0, npt)
    log(f"═══ 实验 C：屋脊跨路线对比（{npt} 点, fnox {FN[0]:.3g}–{FN[-1]:.3g}）═══")
    require_prereq(need_batch=True, need_mcm=True)
    d = out_dir("run_C")

    # ── 路线 1：batch（9×9 网格，取 fvoc≈1 一列即可；或小网格专扫）──
    log("\n[1/2] batch 74sp …")
    # 用 ng=9 一遍拿到 fnox 轴（几何网格 0.3–7.69）
    bp = os.path.join(d, "batch9.csv")
    run_batch(9, bp, verbose=True)
    hdr, data = read_table(bp)
    fvoc = col(hdr, data, "fvoc")
    j1 = float(fvoc[np.argmin(np.abs(fvoc - 1.0))])
    sel = np.where(np.abs(fvoc - j1) < 1e-3)[0]
    sel = sel[np.argsort(col(hdr, data, "fnox")[sel])]
    b_fn = col(hdr, data, "fnox")[sel]
    b_o3 = col(hdr, data, "O3max")[sel]
    b_fnr = col(hdr, data, "FNR_mean")[sel]
    write_csv(os.path.join(d, "ridge_batch.csv"),
              ["fnox", "O3max", "FNR_mean"],
              [[f"{b_fn[k]:.5g}", f"{b_o3[k]:.6g}", f"{b_fnr[k]:.6g}"] for k in range(len(b_fn))],
              note=f"实验C | batch 74sp | fvoc≈{j1:.4g} 切面")

    # ── 路线 2：mcm_run（逐 fnox，逐点跑）──
    log("\n[2/2] mcm 29sp …")
    rows_m, prods = [], []
    for fn in FN:
        pref = os.path.join(d, f"mcm_{fn:.3f}")
        r = run_mcm(5, pref, fnox=fn, fvoc=1.0)
        o3 = last24_mean(f"{pref}_conc.csv", "O3")
        hcho = last24_mean(f"{pref}_conc.csv", "HCHO")
        no2 = last24_mean(f"{pref}_conc.csv", "NO2")
        fnr = hcho / max(no2, 1e-9)
        rows_m.append([f"{fn:.5g}", f"{o3:.6g}", f"{hcho:.6g}", f"{no2:.6g}", f"{fnr:.6g}"])
        prods += [f"{pref}_conc.csv", f"{pref}_irr.csv", f"{pref}_ipr.csv"]
        log(f"    fnox={fn:8.4f}  O3={o3:7.2f} ppb  FNR={fnr:6.2f}")
    write_csv(os.path.join(d, "ridge_mcm.csv"),
              ["fnox", "O3_last24h", "HCHO_last24h", "NO2_last24h", "FNR"], rows_m,
              note="实验C | mcm 29sp | fvoc=1.0 | 末24h均值")
    m_fn = np.array([float(r[0]) for r in rows_m])
    m_o3 = np.array([float(r[1]) for r in rows_m])
    m_fnr = np.array([float(r[4]) for r in rows_m])

    # ── 对比 ──
    ib = int(np.nanargmax(b_o3)); im = int(np.nanargmax(m_o3))
    b_band = [b_fn[k] for k in range(len(b_fnr)) if 1 <= b_fnr[k] <= 2]
    m_band = [m_fn[k] for k in range(len(m_fnr)) if 1 <= m_fnr[k] <= 2]
    rows_c = [
        ["batch 74sp", f"{b_fn[ib]:.4g}", f"{b_o3[ib]:.2f}",
         f"{b_band[0]:.3g}–{b_band[-1]:.3g}" if b_band else "—",
         f"{'重合' if b_band and b_band[0] <= b_fn[ib] <= b_band[-1] else '不重合'}"],
        ["mcm 29sp", f"{m_fn[im]:.4g}", f"{m_o3[im]:.2f}",
         f"{m_band[0]:.3g}–{m_band[-1]:.3g}" if m_band else "—",
         f"{'重合' if m_band and m_band[0] <= m_fn[im] <= m_band[-1] else '不重合'}"],
    ]
    summary_table("实验 C · 跨路线屋脊对比",
                  ["route", "ridge_fnox", "O3max", "FNR_1-2_band", "verdict"], rows_c)
    log("\n  注: 两路线排放谱不同 → 屋脊位置本就不同；有意义的是各自『屋脊 vs 自身 FNR 穿越』。")
    p = write_csv(os.path.join(d, "ridge_compare.csv"),
                  ["route", "ridge_fnox", "O3max_ppb", "FNR_1-2_band_fnox", "verdict"], rows_c,
                  note="实验C | 屋脊位置 vs 指示因子穿越一致性")
    mp = manifest(d, "C_ridge", {"npoints": npt, "fn_range": [FN[0], FN[-1]]},
                  [os.path.join(d, "ridge_batch.csv"), os.path.join(d, "ridge_mcm.csv"), p] + prods)
    log(f"\n产物: {os.path.relpath(d, ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
