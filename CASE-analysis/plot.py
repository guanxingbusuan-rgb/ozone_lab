#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CASE-analysis/plot.py —— 个例分析图（只读 output/ 的 CSV，不出数据）

三张图（缺对应 CSV 时跳过）:
  output/case_ipr.png   目标物种六账随时间（P/L/E/D/B/X + 净值）
  output/case_irr.png   IRR 反应级净贡献 top-12（横条）
  output/case_rir.png   前体物 RIR（组均值横条 + 逐时曲线）

用法:
  python CASE-analysis/plot.py [--species O3]
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "output")
sys.path.insert(0, os.path.join(ROOT, "analysis_lib"))
from csvio import read_rows, fcol                   # noqa: E402

PROCS = ["P", "L", "E", "D", "B", "X"]


def plot_ipr(sp="O3"):
    p = os.path.join(OUT, "ipr_ledger.csv")
    if not os.path.exists(p):
        return None
    h, R = read_rows(p)
    sps = [r[h.index("species")] for r in R]
    if sp not in sps:
        print(f"  [skip] ipr_ledger 无 {sp}"); return None
    rows = [r for r in R if r[h.index("species")] == sp]
    hr = np.array([float(r[h.index("h")]) for r in rows])
    vals = {k: np.array([float(r[h.index(k)]) for r in rows]) for k in PROCS}
    net = sum(vals.values())
    fig, ax = plt.subplots(1, 2, figsize=(12.6, 4.6))
    for k, mk in zip(PROCS, ["-", "-", "--", "-", "-", "-"]):
        ax[0].plot(hr, vals[k], mk, lw=1.4, label=k)
    ax[0].set_yscale("symlog", linthresh=1e6)
    ax[0].set_xlabel("output step h"); ax[0].set_ylabel("window integral [molec/cm3]")
    ax[0].set_title(f"({sp}) six accounts — P mech prod, L mech loss, E emis,"
                    " D dep, B nudging+dilution, X BLH")
    ax[0].legend(fontsize=8, ncol=2); ax[0].grid(alpha=.25)
    ax[1].plot(hr, net, "k-", lw=1.6)
    ax[1].axhline(0, color="grey", lw=.8)
    ax[1].set_xlabel("output step h"); ax[1].set_ylabel("net [molec/cm3]")
    ax[1].set_title(f"({sp}) net = P+L+E+D+B+X")
    ax[1].grid(alpha=.25)
    out = os.path.join(OUT, "case_ipr.png")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    return out


def plot_irr():
    p = os.path.join(OUT, "irr_top.csv")
    if not os.path.exists(p):
        return None
    h, R = read_rows(p)
    net = fcol(h, R, h[[i for i, n in enumerate(h) if n.startswith("net_")][0]])
    lab = [f"r{r[h.index('r_index')]}  {r[h.index('reaction')][:46]}" for r in R]
    order = np.argsort(np.abs(net))
    fig, ax = plt.subplots(figsize=(10.4, 5.2))
    ax.barh(np.arange(len(net)), net[order],
            color=["#c0392b" if v < 0 else "#2471a3" for v in net[order]])
    ax.set_yticks(np.arange(len(net))); ax.set_yticklabels([lab[i] for i in order], fontsize=7.5)
    ax.axvline(0, color="k", lw=.8)
    ax.set_xlabel("net contribution (stoich × Σ window IRR)")
    ax.set_title("Reaction-level attribution — top-12 by |net| (blue = production, red = loss)")
    ax.grid(axis="x", alpha=.25)
    out = os.path.join(OUT, "case_irr.png")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    return out


def plot_rir():
    pm = os.path.join(OUT, "rir_mean.csv")
    if not os.path.exists(pm):
        return None
    h, R = read_rows(pm)
    g = [r[h.index("group")] for r in R]
    v = fcol(h, R, "RIR_mean")
    fig, ax = plt.subplots(1, 2, figsize=(13.2, 4.8))
    ax[0].barh(np.arange(len(g)), v,
               color=["#1f77b4" if x >= 0 else "#d62728" for x in v])
    ax[0].set_yticks(np.arange(len(g))); ax[0].set_yticklabels(g, fontsize=8)
    ax[0].axvline(0, color="k", lw=.8)
    ax[0].set_xlabel("RIR (P(O3)-weighted avg, Eq.8)")
    ax[0].set_title("Precursor RIR — group mean")
    ax[0].grid(axis="x", alpha=.25)
    pd_ = os.path.join(OUT, "rir_diurnal.csv")
    if os.path.exists(pd_):
        h2, R2 = read_rows(pd_)
        for name in g[:5]:
            rows = [r for r in R2 if r[h2.index("group")] == name]
            if rows:
                ax[1].plot([float(r[h2.index("hour")]) for r in rows],
                           [float(r[h2.index("RIR_avg")]) for r in rows],
                           "-o", ms=3, lw=1.3, label=name)
        ax[1].set_xlabel("hour of day"); ax[1].set_ylabel("RIR_avg")
        ax[1].set_title("RIR diurnal — top-5 groups")
        ax[1].legend(fontsize=8); ax[1].grid(alpha=.25)
    out = os.path.join(OUT, "case_rir.png")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="O3")
    a = ap.parse_args()
    made = [p for p in (plot_ipr(a.species), plot_irr(), plot_rir()) if p]
    if not made:
        print("✗ output/ 无可绘 CSV —— 先跑 python CASE-analysis/run.py", file=sys.stderr)
        return 1
    for p in made:
        print(os.path.relpath(p, ROOT), "完成")
    return 0


if __name__ == "__main__":
    sys.exit(main())
