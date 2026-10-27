#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fig_response_v3.py —— O3 排放空间响应面 3 面板图 (屋脊 / FNR / 屋脊切片)
数据: output/batch_v3_indicators.csv (或旧 7 列 batch_v3.csv — 取前几列)
用法: python fig_response_v3.py [csv] [out.png]
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

csv_path = sys.argv[1] if len(sys.argv) > 1 else "../output/batch_v3_indicators.csv"
out_png = sys.argv[2] if len(sys.argv) > 2 else "response_surface_v3.png"

lines = [l for l in open(csv_path) if l.strip() and not l.startswith("fnox")]
recs = [[float(x) for x in l.split(",")] for l in lines]
arr = np.array(recs)
fx = sorted(set(arr[:, 0])); fv = sorted(set(arr[:, 1]))
ni, nj = len(fx), len(fv)

def grid(col):
    Z = np.full((ni, nj), np.nan)
    for r in arr:
        Z[fx.index(r[0]), fv.index(r[1])] = r[col]
    return Z

O3, F = grid(4), grid(8)   # 24 列: O3max=4, FNR=8
X, Y = np.meshgrid(fv, fx)
ridge_fv, ridge_fx = [], []
for j in range(nj):
    col = O3[:, j]
    if np.isfinite(col).any():
        ridge_fx.append(fx[int(np.nanargmax(col))]); ridge_fv.append(fv[j])

fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.2))
im0 = axes[0].contourf(X, Y, O3, levels=14, cmap="turbo")
axes[0].contour(X, Y, O3, levels=[27, 28], colors="white", linewidths=1.5, linestyles="--")
plt.colorbar(im0, ax=axes[0], label="O3max [ppb]")
axes[0].set_title("(a) O3 response (white=ridge)")
im1 = axes[1].contourf(X, Y, np.log10(np.maximum(F, 1e-6)), levels=14, cmap="RdYlBu_r")
cs1 = axes[1].contour(X, Y, F, levels=[1, 2, 4], colors="k", linewidths=1.2)
axes[1].clabel(cs1, fmt="FNR=%g")
plt.colorbar(im1, ax=axes[1], label="log10 FNR")
axes[1].set_title("(b) FNR = HCHO/NO2 indicator")
j1 = min(range(nj), key=lambda j: abs(fv[j] - 1.0))
axes[2].plot(fx, O3[:, j1], "o-", color="#C0392B", lw=2)
axes[2].axvline(1.5, ls="--", color="gray")
axes[2].text(1.55, np.nanmax(O3[:, j1])*0.92, "ridge (NOx-saturated)", fontsize=8)
axes[2].set_title("(c) O3-NOx ridge: NOx-limited to saturated")
for ax in axes:
    ax.set_xlabel("f_VOC"); ax.set_ylabel("f_NOx")
    ax.set_xscale("log"); ax.set_yscale("log")
axes[2].set_ylabel("O3max [ppb]"); axes[2].set_yscale("linear")
plt.tight_layout()
plt.savefig(out_png, dpi=170, bbox_inches="tight")
print(f"{out_png} 完成")
