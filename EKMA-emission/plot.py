#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
EKMA-emission/plot.py —— S1 响应面可视化（读 output/ekma_emission_grid.csv）

图面:
  横轴 = NMVOC 排放扰动 [%]（0 = 无扰动，200 = +200%）
  纵轴 = NOx  排放扰动 [%]（0 = 无扰动，200 = +200%）
  色标 = O3max [ppb]（末 24 步窗口最大，不是 8h 日最大）；黑实线 = O3max 等值线；白虚线 = 屋脊（每列 VOC 的 O3 峰 NOx）

产物（output/）:
  ekma_emission_S1.png         主图（响应面 + 屋脊）
  ekma_emission_slices.png     NOx 切片线图

屋脊轨迹 ekma_emission_ridge.csv 与派生表 ekma_emission_derived.csv 由 run.py 写出；
本脚本只读 output/ekma_emission_grid.csv，不产出 CSV。

用法:
  python plot.py                （在 EKMA-emission/ 目录内）
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "output")

sys.path.insert(0, os.path.join(ROOT, "analysis_lib"))
from csvio import read_table, write_csv, fcol  # noqa: E402


def main():
    src = os.path.join(OUT, "ekma_emission_grid.csv")
    if not os.path.exists(src):
        print(f"✗ 缺 {os.path.relpath(src, ROOT)} —— 先跑 python EKMA-emission/run.py",
              file=sys.stderr)
        return 1
    hdr, D = read_table(src)
    x = fcol(hdr, D, "pert_voc_pct")      # NMVOC
    y = fcol(hdr, D, "pert_nox_pct")      # NOx
    z = fcol(hdr, D, "O3max_ppb")

    xs = np.array(sorted(set(np.round(x, 4))))
    ys = np.array(sorted(set(np.round(y, 4))))
    if len(xs) * len(ys) != len(x):
        print(f"✗ 不是完整正交网格（{len(xs)}×{len(ys)}={len(xs)*len(ys)} vs {len(x)} 行）",
              file=sys.stderr)
        return 1
    # 键序必须与 x/y 一致: x=pert_voc, y=pert_nox → idx[(voc, nox)]
    idx = {(round(a, 4), round(b, 4)): k for k, (a, b) in enumerate(zip(x, y))}
    Z = np.array([[z[idx[(a, b)]] for b in ys] for a in xs])   # Z[voc_idx, nox_idx]

    # ── 屋脊：每个 VOC 列取 O3max 峰对应的 NOx 扰动 ──
    ridge_nox, ridge_o3 = [], []
    for j, xv in enumerate(xs):
        col = Z[j, :]
        if np.isfinite(col).any():
            k = int(np.nanargmax(col))
            ridge_nox.append(ys[k]); ridge_o3.append(col[k])
    # ── 主图 ──
    fig, ax = plt.subplots(figsize=(9.2, 7.6))
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    lv = 18
    cs = ax.contourf(X, Y, Z, levels=lv, cmap="turbo")
    cb = plt.colorbar(cs, ax=ax, pad=0.02)
    cb.set_label("O3max [ppb]  (max over the last 24 output steps)", fontsize=10)

    zlo, zhi = float(np.nanmin(Z)), float(np.nanmax(Z))
    clv = [v for v in MaxNLocator(7).tick_values(zlo, zhi) if zlo < v < zhi]
    if clv:
        c2 = ax.contour(X, Y, Z, levels=clv, colors="k", linewidths=0.85)
        ax.clabel(c2, fmt="%g", fontsize=8, inline=True)

    ax.plot(xs, ridge_nox, "w--", lw=2.4, label="O3 ridge")
    ax.plot(xs, ridge_nox, "k-", lw=1.2, label=None)
    ax.scatter([0], [0], marker="+", s=180, color="k", zorder=6,
               label="baseline (0%, 0%)")

    # 高低 NOx 两侧的体制提示（定性箭头文字，不依赖判据切点）
    imax = int(np.nanargmax(Z))
    jj, kk = np.unravel_index(imax, Z.shape)
    ax.annotate(f"max O3 = {Z[jj, kk]:.1f} ppb\n(NMVOC +{xs[jj]:.0f}%, NOx +{ys[kk]:.0f}%)",
                xy=(xs[jj], ys[kk]), xytext=(0.98, 0.97), textcoords="axes fraction",
                fontsize=9, color="w", ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="k", alpha=0.35, ec="none"))

    ax.set_xlabel(f"NMVOC emission perturbation [%]  (0 = baseline, +{xs[-1]:g}% = x{1 + xs[-1] / 100:.2f})",
                  fontsize=11)
    ax.set_ylabel(f"NOx emission perturbation [%]  (0 = baseline, +{ys[-1]:g}% = x{1 + ys[-1] / 100:.2f})",
                  fontsize=11)
    ridge_note = ""
    if ridge_nox and max(ridge_nox) >= ys[-1] - 1e-6:
        n_edge = sum(1 for v in ridge_nox if v >= ys[-1] - 1e-6)
        ridge_note = (f"\n⚠ ridge at the +{ys[-1]:g}% NOx edge for {n_edge}/{len(ridge_nox)} "
                      f"VOC columns — true peak may lie beyond the sampled range")
    ax.set_title("S1 · EKMA emission response surface — O3max [ppb]\n"
                 f"{len(xs)}×{len(ys)} = {len(x)} orthogonal cases "
                 f"(0–{xs[-1]:g}% in {xs[1]-xs[0]:g}% steps)\n"
                 "fixed: mechanism V2.5+ · met/IC from config_override/ (effective input: "
                 "output/runtime_emission.dat, see run_meta.txt) · "
                 "response window = last 24 output steps (last day)"
                 + ridge_note,
                 fontsize=10)
    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.15, color="w")

    # ── 附加图: O3max–NOx 切片（每个选定 NMVOC 水平一条线，屋脊=各线峰） ──
    fig2, ax2 = plt.subplots(figsize=(9.2, 6.2))
    picks = [v for v in (0, 25, 50, 100, 150, 200) if v in xs]
    cmap = plt.get_cmap("viridis")
    for q, v in enumerate(picks):
        j = int(np.where(xs == v)[0][0])
        col = np.array([Z[j, :]], dtype=float)[0]
        ax2.plot(ys, col, "-o", ms=2.6, lw=1.6,
                 color=cmap(q / max(len(picks) - 1, 1)),
                 label=f"NMVOC +{v:g}%")
        k = int(np.nanargmax(col))
        ax2.plot(ys[k], col[k], "*", ms=11, color=cmap(q / max(len(picks) - 1, 1)),
                 mec="k", mew=0.6)
    ax2.set_xlabel("NOx emission perturbation [%]  (0 = baseline)", fontsize=11)
    ax2.set_ylabel("O3max [ppb]", fontsize=11)
    ax2.set_title("S1 slices · O3max vs NOx at selected NMVOC levels\n"
                  "star = local peak (the ridge point for that NMVOC column; "
                  "off-right star = peak beyond sampled max)", fontsize=11)
    ax2.grid(alpha=0.3); ax2.legend(fontsize=9)
    png2 = os.path.join(OUT, "ekma_emission_slices.png")
    fig2.tight_layout(); fig2.savefig(png2, dpi=170, bbox_inches="tight"); plt.close(fig2)

    png = os.path.join(OUT, "ekma_emission_S1.png")
    fig.tight_layout()
    fig.savefig(png, dpi=170, bbox_inches="tight")
    plt.close(fig)

    print(f"  O3max {zlo:.2f} – {zhi:.2f} ppb")
    print(f"  产物: {os.path.relpath(png, ROOT)}")
    print(f"        {os.path.relpath(png2, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
