#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""fig_ridge5.py —— 五路线指示因子屋脊实验图 (2×2, 归一化)"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
while not os.path.isdir(os.path.join(HERE, 'analysis_lib')):
    _up = os.path.dirname(HERE)
    if _up == HERE:
        raise SystemExit('cannot locate project root (analysis_lib/ not found)')
    HERE = _up
sys.path.insert(0, f"{HERE}/analysis_lib")
from indicators import route_timeseries  # noqa

# R13: Fortran 三位指数(3.07932-174)兜底 — 共享实现
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) if _os.path.basename(_os.path.dirname(_os.path.abspath(__file__))) == "experiments" else _os.path.dirname(_os.path.abspath(__file__)), "analysis_lib"))
from csvfloat import fix_fortran_exp as _fixexp, fnum as _fnum

FN = np.logspace(-1.5, 1.5, 15)
PANELS = [
    ("mcm_run 29sp MCM emission-space", "mc_{:.2f}", f"{HERE}/mcm_city/input/mechanism_mcm_ch4.fac"),
]
KEYS = ["FNR", "H2O2_HNO3", "O3_NOy", "O3_NOz", "OPE"]
COLORS = dict(zip(KEYS, ["#E74C3C", "#8E44AD", "#2980B9", "#27AE60", "#F39C12"]))

fig, ax = plt.subplots(figsize=(9, 5.2))
for (title, pat, fac), axp in zip(PANELS, [ax]):
    data = {}
    for f in FN:
        pref = f"{HERE}/output/ridge5/" + pat.format(f)
        try:
            data[f] = route_timeseries(pref, fac)
        except Exception:
            pass
    fns = sorted(data)
    for k in KEYS:
        v = np.array([data[f][k].mean() for f in fns])
        rng = v.max() - v.min()
        vn = (v - v.min()) / rng if rng > 0 else v * 0 + 0.5
        axp.plot(fns, vn, "o-", color=COLORS[k], ms=3, label=k, lw=1.4)
    o3 = np.array([data[f]["O3_ppb"].mean() for f in fns])
    axp.plot(fns, (o3 - o3.min()) / (o3.max() - o3.min()), "k-", lw=2.6, alpha=0.65, label="O3")
    ir = int(np.argmax(o3))
    axp.axvline(fns[ir], color="k", ls="--", lw=1)
    axp.text(fns[ir] * 1.1, 0.04, f"ridge\n{fns[ir]:.2f}", fontsize=7)
    axp.set_xscale("log"); axp.set_ylim(-0.05, 1.1)
    axp.set_title(title, fontsize=10)
    axp.set_xlabel("f_NOx"); axp.set_ylabel("normalized")
    axp.legend(fontsize=6, ncol=2)
fig.suptitle("Ridge5 experiment: indicator transitions vs O3 ridge per route "
             "(dashed = O3 ridge; 74sp batch panel omitted - see indicators_vs_ridge.png)", fontsize=11)
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig(f"{HERE}/figures/ridge5_routes.png", dpi=160, bbox_inches="tight")
print("figures/ridge5_routes.png 完成")
