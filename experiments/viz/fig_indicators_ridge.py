#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fig_indicators_ridge.py —— 指示因子全家族 × O3 屋脊 9 面板图
数据: output/batch_v3_indicators.csv (fac_batch_omp 81 案例, 16 列)
用法: python fig_indicators_ridge.py [csv] [out.png]
面板: O3max / FNR / H2O2-HNO3 / O3-NOy / HCHO-NOy / NOy / L_NOx-Q / P(H2O2)-P(HNO3) / OPE
每面板: 黑色等值线 = Sillman 1995 / Kleinman 2005 阈值; 白虚线 = O3 屋脊 (每 fvoc 列 O3max 峰)
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

csv_path = sys.argv[1] if len(sys.argv) > 1 else "../output/batch_v3_indicators.csv"
out_png = sys.argv[2] if len(sys.argv) > 2 else "indicators_vs_ridge.png"

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

# 24 列格式: 0fnox 1fvoc 2jfac 3dt 4O3max 5O3mean 6HCHOmax 7NO2mean
#            8FNR 9H2O2HNO3 10O3NOy 11HCHONOy 12O3NOz 13NOy 14LnOxQ 15PH2O2PHNO3 16OPE 17-22pm 23nst
O3 = grid(4)
X, Y = np.meshgrid(fv, fx)
ridge_fv, ridge_fx = [], []
for j in range(nj):
    col = O3[:, j]
    if np.isfinite(col).any():
        ridge_fx.append(fx[int(np.nanargmax(col))]); ridge_fv.append(fv[j])

panels = [
    ("(a) O3max [ppb]",        O3,       None,            "turbo",    False),
    ("(b) FNR = HCHO/NO2",     grid(8),  [1.0, 2.0, 4.0], "RdYlBu_r", False),
    ("(c) H2O2/HNO3",          grid(9),  [0.4, 0.8],      "RdYlBu_r", False),
    ("(d) O3/NOy",             grid(10),  [10, 15],        "RdYlBu_r", False),
    ("(e) HCHO/NOy",           grid(11),  [0.28, 0.32],    "RdYlBu_r", False),
    ("(f) NOy [ppb]",          grid(13), None,            "viridis",  False),
    ("(g) L_NOx/Q (Kleinman)", grid(14), [0.5],           "coolwarm", False),
    ("(h) P(H2O2)/P(HNO3)",    grid(15), [0.15, 0.5],     "RdYlBu_r", True),
    ("(i) OPE = P(O3)/P(NOz)", grid(16), None,            "viridis",  True),
]
fig, axes = plt.subplots(3, 3, figsize=(16.5, 13.5))
for ax, (title, Z, levels, cmap, logc) in zip(axes.flat, panels):
    ZZ = np.log10(np.maximum(Z, 1e-9)) if logc else Z
    im = ax.contourf(X, Y, ZZ, levels=14, cmap=cmap)
    plt.colorbar(im, ax=ax, label="log10" if logc else "")
    if levels:
        cs = ax.contour(X, Y, Z, levels=levels, colors="k", linewidths=1.3)
        ax.clabel(cs, fmt="%g", fontsize=7)
    ax.plot(ridge_fv, ridge_fx, "w--", lw=2.4, label="O3 ridge")
    ax.plot(ridge_fv, ridge_fx, "k-", lw=1.1)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("f_VOC"); ax.set_ylabel("f_NOx")
    ax.set_xscale("log"); ax.set_yscale("log")
axes.flat[0].legend(loc="upper left", fontsize=8, framealpha=0.9)
fig.suptitle("Ozone Ridge vs Indicator Family — Emission-Space Response Surface\n"
             "mechanism_3 (53 sp, 96 rxn) · black contours = Sillman-1995 / Kleinman-2005 "
             "thresholds · OPE = P(O3)/P(NOz) with HNO3 single-counted", fontsize=12)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig(out_png, dpi=170, bbox_inches="tight")
print(f"{out_png} 完成")
