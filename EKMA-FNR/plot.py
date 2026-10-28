#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
EKMA-FNR/plot.py —— S2 FNR 响应面可视化（读 output/ekma_fnr_grid.csv）

图面（IC 扰动口径，轴 = 控制变量 → 铺满）:
  横轴 = 初始 NO₂ 扰动 [%]（0 = config_override 的 NO₂ 基线）—— 分母
         ic1 是 NOx 组：NO 与 NO₂ 同因子缩放，NO₂:NO 比例恒定（图上标出实际比例）
  纵轴 = 初始 HCHO 扰动 [%]（0 = config_override 的 HCHO 基线）—— 分子
  基线 ppb 数值从 output/fnr_ic_base.csv 读取，换基线后标签自动跟随
  范围/步长由 run.py 的 --pmax/--step 决定（缺省 0–1000% / 100%）
  色标 = O3max [ppb]；黑实线 = O3max 等值线；白虚线 = FNR 等值线（1/2/4，Sillman）
  白点虚线 = 屋脊（每列 HCHO 的 O₃ 峰对应 NO₂ 扰动）

产物（output/）:
  ekma_fnr_S2.png              主图
  ekma_fnr_slices.png          O₃max–NO₂ 切片线图

屋脊轨迹 ekma_fnr_ridge.csv 与派生表 ekma_fnr_derived.csv 由 run.py 写出；
本脚本只读 output/ekma_fnr_grid.csv，不产出 CSV。

用法:
  python plot.py                （在 EKMA-FNR/ 目录内）
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
from csvio import read_table, read_rows, write_csv, fcol  # noqa: E402


def main():
    src = os.path.join(OUT, "ekma_fnr_grid.csv")
    if not os.path.exists(src):
        print(f"✗ 缺 {os.path.relpath(src, ROOT)} —— 先跑 python EKMA-FNR/run.py",
              file=sys.stderr)
        return 1
    hdr, D = read_table(src)
    x = fcol(hdr, D, "pert_no2_pct")       # NO₂  ← ic1
    y = fcol(hdr, D, "pert_hcho_pct")      # HCHO ← ic2
    z = fcol(hdr, D, "O3max_ppb")
    fr = fcol(hdr, D, "FNR_mean")

    xs = np.array(sorted(set(np.round(x, 4))))   # NO₂ 轴
    ys = np.array(sorted(set(np.round(y, 4))))   # HCHO 轴
    if len(xs) * len(ys) != len(x):
        print(f"✗ 不是完整正交网格（{len(xs)}×{len(ys)} vs {len(x)} 行）", file=sys.stderr)
        return 1
    idx = {(round(a, 4), round(b, 4)): k for k, (a, b) in enumerate(zip(x, y))}
    Z = np.array([[z[idx[(a, b)]] for b in ys] for a in xs])   # Z[no2_idx, hcho_idx]
    R = np.array([[fr[idx[(a, b)]] for b in ys] for a in xs])

    # ── 基线 ppb 从本轮 IC 存档读取（换 config_override 基线后轴标签自动跟随） ──
    base = {"NO2": float("nan"), "HCHO": float("nan"), "NO": float("nan")}
    ic_path = os.path.join(OUT, "fnr_ic_base.csv")
    if os.path.exists(ic_path):
        ih, irows = read_rows(ic_path)
        for row in irows:
            sp = row[ih.index("species")]
            if sp in base:
                base[sp] = float(row[ih.index("value_cm3")]) / 2.46e10
    b_no2, b_hcho, b_no = base["NO2"], base["HCHO"], base["NO"]
    def ppb(v):
        return f"{v:g} ppb" if v == v else "baseline n/a"
    ratio = f"{b_no2 / b_no:g}" if (b_no == b_no and b_no and b_no2 == b_no2) else "n/a"

    # ── 屋脊: 每个 HCHO 列取 O₃max 峰对应的 NO₂ 扰动 ──
    ridge_no2, ridge_o3 = [], []
    for j, yv in enumerate(ys):
        col = Z[:, j]
        if np.isfinite(col).any():
            k = int(np.nanargmax(col))
            ridge_no2.append(xs[k]); ridge_o3.append(col[k])
    # ── 主图（x = 初始 NO₂ 扰动 ← 分母, y = 初始 HCHO 扰动 ← 分子）──
    fig, ax = plt.subplots(figsize=(9.2, 7.6))
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    cs = ax.contourf(X, Y, Z, levels=18, cmap="turbo")
    cb = plt.colorbar(cs, ax=ax, pad=0.02)
    cb.set_label("O3max [ppb]", fontsize=11)

    zlo, zhi = float(np.nanmin(Z)), float(np.nanmax(Z))
    clv = [v for v in MaxNLocator(7).tick_values(zlo, zhi) if zlo < v < zhi]
    if clv:
        c2 = ax.contour(X, Y, Z, levels=clv, colors="k", linewidths=0.85)
        ax.clabel(c2, fmt="%g", fontsize=8, inline=True)

    # FNR 等值线（轴=控制变量 → 比值不再是射线，改画等值线）
    rin = [g for g in (1.0, 2.0, 4.0) if np.nanmin(R) < g < np.nanmax(R)]
    if rin:
        cr = ax.contour(X, Y, R, levels=rin, colors="w", linestyles="--",
                        linewidths=1.8)
        ax.clabel(cr, fmt="FNR=%g", fontsize=8, inline=True)

    ax.plot(ridge_no2, ys, "w.", ms=9, mec="k", mew=0.7, label="O3 ridge")
    ax.scatter([0], [0], marker="+", s=180, color="k", zorder=6,
               label="baseline (0%, 0%)")

    imax = int(np.nanargmax(Z))
    jj, kk = np.unravel_index(imax, Z.shape)
    ax.annotate(f"max O3 = {Z[jj, kk]:.1f} ppb\n(HCHO +{ys[kk]:.0f}%, NO₂ +{xs[jj]:.0f}%)",
                xy=(xs[jj], ys[kk]), xytext=(0.98, 0.97), textcoords="axes fraction",
                fontsize=9, color="w", ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="k", alpha=0.35, ec="none"))

    ax.set_xlabel(f"initial NO₂ perturbation [%]  (baseline {ppb(b_no2)}, 0 = no perturbation)",
                  fontsize=11)
    ax.set_ylabel(f"initial HCHO perturbation [%]  (baseline {ppb(b_hcho)}, 0 = no perturbation)",
                  fontsize=11)
    ax.set_title("S2 · EKMA-FNR response surface — O3max [ppb]\n"
                 f"{len(xs)}×{len(ys)} = {len(x)} orthogonal IC cases "
                 f"(0–{xs[-1]:g}% in {xs[1]-xs[0]:g}% steps)\n"
                 "fixed: mechanism V2.5+ · met/IC from config_override/ (effective input: "
                 "output/runtime_fnr.dat, see run_meta.txt) · "
                 f"emissions OFF · NO scales with the NOx group "
                 f"(NO2:NO = {ratio}:1 fixed) · O₃ not perturbed",
                 fontsize=11)
    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.15, color="w")

    png = os.path.join(OUT, "ekma_fnr_S2.png")
    fig.tight_layout()
    fig.savefig(png, dpi=170, bbox_inches="tight")
    plt.close(fig)

    # ── 切片图: O₃max vs NO₂（选定 HCHO 水平）──
    fig2, ax2 = plt.subplots(figsize=(9.2, 6.2))
    picks = [v for v in (0, 100, 200, 300, 400) if v in ys]
    cmap = plt.get_cmap("viridis")
    for q, v in enumerate(picks):
        j = int(np.where(ys == v)[0][0])
        colv = Z[:, j]
        ax2.plot(xs, colv, "-o", ms=2.6, lw=1.6,
                 color=cmap(q / max(len(picks) - 1, 1)), label=f"HCHO +{v:g}%")
        k = int(np.nanargmax(colv))
        ax2.plot(xs[k], colv[k], "*", ms=11,
                 color=cmap(q / max(len(picks) - 1, 1)), mec="k", mew=0.6)
    ax2.set_xlabel(f"initial NO₂ perturbation [%]  (baseline {ppb(b_no2)})", fontsize=11)
    ax2.set_ylabel("O3max [ppb]", fontsize=11)
    ax2.set_title("S2 slices · O3max vs initial NO₂ at selected HCHO levels\n"
                  "star = local peak (ridge point)", fontsize=11)
    ax2.grid(alpha=0.3); ax2.legend(fontsize=9)
    png2 = os.path.join(OUT, "ekma_fnr_slices.png")
    fig2.tight_layout(); fig2.savefig(png2, dpi=170, bbox_inches="tight")
    plt.close(fig2)

    print(f"  O3max {zlo:.2f} – {zhi:.2f} ppb")
    print(f"  产物: {os.path.relpath(png, ROOT)}")
    print(f"        {os.path.relpath(png2, ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
