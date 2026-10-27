#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_B_surface.py —— 实验 B：O3–NOx–VOC 响应面（4D 参数扫描）+ 指示因子（无绘图）

输入假设: runtime_v3.dat 已编译
输出:
  output/run_B/batch.csv              原始 24 列（fnox,fvoc,jfac,dt + 19 响应量 + nst）
  output/run_B/indicator_matrix_<var>.csv   每个响应量转成 fnox×fvoc 矩阵（便于直接画/读）
  output/run_B/slices.csv             fnox 切面（fvoc≈1）逐点序列
  （终端另打印体制诊断：屋脊 / FNR 1-2 带 / L_NOx-Q 穿 0.5 / OPE 谷）

用法:
  python experiments/run/run_B_surface.py                     # 9×9 二维
  python experiments/run/run_B_surface.py 13                  # 13×13
  python experiments/run/run_B_surface.py 9 0.6,0.8,1.0,1.2   # 加光解维
  python experiments/run/run_B_surface.py 9 1.0 0,5,10        # 加温度维
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _runner import (ROOT, require_prereq, out_dir, run_batch, read_table,
                     write_csv, manifest, summary_table, log, col)

INDICATORS = ["O3max", "O3mean", "FNR_mean", "H2O2HNO3", "O3NOy", "HCHONOy",
              "O3NOz", "NOy_ppb", "LnOxQ", "PH2O2PHNO3", "OPE",
              "FNR_pm", "H2O2HNO3_pm", "O3NOy_pm", "HCHONOy_pm", "O3NOz_pm", "NOy_pm"]


def main(argv):
    ng = int(argv[1]) if len(argv) > 1 else 9
    jf = [float(x) for x in argv[2].split(",")] if len(argv) > 2 else None
    dts = [float(x) for x in argv[3].split(",")] if len(argv) > 3 else None

    log(f"═══ 实验 B：响应面 {ng}×{ng}"
        f"{' × jfac'+str(jf) if jf else ''}{' × dT'+str(dts) if dts else ''} ═══")
    require_prereq(need_batch=True)
    d = out_dir("run_B")
    pref = os.path.join(d, "batch.csv")

    r = run_batch(ng, pref, jfac_list=jf, dt_list=dts, verbose=True)
    products = [pref, f"{pref}_ipr.csv", f"{pref}_irr.csv"]

    hdr, data = read_table(pref)
    ncase = len(data)
    log(f"  案例数 {ncase} · 列数 {len(hdr)} · 耗时 {r['elapsed_s']:.1f}s")

    # ── 每个响应量 → fnox×fvoc 矩阵（仅当只有 fnox/fvoc 两轴时）──
    have = [v for v in INDICATORS if v in hdr]
    if not jf and not dts:
        fn = sorted({round(x, 4) for x in col(hdr, data, "fnox")})
        fv = sorted({round(x, 4) for x in col(hdr, data, "fvoc")})
        idx = {(round(a, 4), round(b, 4)): k
               for k, (a, b) in enumerate(zip(col(hdr, data, "fnox"),
                                              col(hdr, data, "fvoc")))}
        for var in have:
            vals = col(hdr, data, var)
            rows = []
            for a in fn:
                rows.append([a] + [f"{vals[idx[(a, b)]]:.6g}" for b in fv])
            p = write_csv(os.path.join(d, f"indicator_matrix_{var}.csv"),
                          ["fnox\\fvoc"] + [f"{b:.4g}" for b in fv], rows,
                          note=f"实验B | 响应量 {var} | 行=fnox 列=fvoc | 空值=NaN")
            products.append(p)
        log(f"  矩阵文件 {len(have)} 个 (每个响应量一个)")

    # ── fnox 切面（fvoc 最接近 1 的那一列）──
    fvoc_all = col(hdr, data, "fvoc")
    j1 = float(fvoc_all[np.argmin(np.abs(fvoc_all - 1.0))])
    sel = np.where(np.abs(fvoc_all - j1) < 1e-3)[0]
    sel = sel[np.argsort(col(hdr, data, "fnox")[sel])]
    rows = [[f"{col(hdr, data, 'fnox')[k]:.4g}"] +
            [f"{col(hdr, data, v)[k]:.6g}" if v in hdr else "nan" for v in have]
            for k in sel]
    p = write_csv(os.path.join(d, "slices.csv"), ["fnox"] + have, rows,
                  note=f"实验B | fvoc≈{j1:.4g} 切面 | 用于看屋脊与指标穿越")
    products.append(p)

    # ── 终端摘要：屋脊与 FNR 穿越 ──
    o3 = col(hdr, data, "O3max")[sel]
    fnox = col(hdr, data, "fnox")[sel]
    try:
        fnr = col(hdr, data, "FNR_mean")[sel]
    except KeyError:
        fnr = np.full_like(o3, np.nan)
    iridge = int(np.nanargmax(o3))
    summary_table(f"实验 B · fvoc≈{j1:.3f} 切面（屋脊定位）",
                  ["fnox", "O3max", "FNR"],
                  [[round(float(fnox[k]), 3), round(float(o3[k]), 2),
                    round(float(fnr[k]), 3)] for k in range(len(fnox))])
    log(f"\n  ★ 屋脊: fnox={fnox[iridge]:.3f} 处 O3max={o3[iridge]:.2f} ppb")
    in_band = [k for k in range(len(fnr)) if 1.0 <= fnr[k] <= 2.0]
    if in_band:
        log(f"  ★ FNR 落 1–2 带: fnox ∈ [{fnox[in_band[0]]:.3f}, {fnox[in_band[-1]]:.3f}]"
            f"  → 与屋脊 {'重合 ✓' if in_band[0] <= iridge <= in_band[-1] else '不重合'}")

    # 体制诊断补充（原 scan_suite.py 的两项独有输出，已并入本实验）
    try:
        lq = col(hdr, data, "LnOxQ")[sel]
        i_lq = next((k for k in range(len(lq) - 1) if lq[k] < 0.5 <= lq[k + 1]), None)
        if i_lq is not None:
            log(f"  ★ L_NOx/Q 穿 0.5 (Kleinman): fnox {fnox[i_lq]:.3f} → {fnox[i_lq + 1]:.3f}")
    except KeyError:
        pass
    try:
        ope = col(hdr, data, "OPE")[sel]
        i_ope = int(np.nanargmin(ope))
        log(f"  ★ OPE 谷: fnox={fnox[i_ope]:.3f} ({ope[i_ope]:.2f})  [谷=NOx 饱和标志]")
    except KeyError:
        pass

    mp = manifest(d, "B_surface",
                  {"ng": ng, "jfac": jf, "dt": dts, "cases": ncase},
                  products + [mp] if False else products)
    log(f"\n产物: {os.path.relpath(d, ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
