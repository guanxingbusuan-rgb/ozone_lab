#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CASE-analysis/plot.py —— 个例分析图（只读 output/ 的 CSV，不出数据）

五张图（缺对应 CSV 时跳过）:
  output/case_ipr.png   六过程三板：逐时堆叠柱 + 六线 + 净值（后向窗 0–23）
  output/case_irr.png   IRR 反应级净贡献 top-12（横条）
  output/case_rir.png   前体物 RIR（组均值横条 + 逐时曲线）
  output/case_o3_budget.png  O₃ 收支堆叠柱（Eq4–6，图 4 同款）
  output/case_rir_fig5.png   图 5 同款 RIR 分组柱图（逐时 0–23，夜间自动剔除）

用法:
  python CASE-analysis/plot.py [--species O3]
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "output")
sys.path.insert(0, os.path.join(ROOT, "analysis_lib"))
from csvio import read_rows, fcol, PPB               # noqa: E402

PROCS = ["P", "L", "E", "D", "B", "X"]
IPR_NAME = {"P": "mech prod", "L": "mech loss", "E": "emission",
            "D": "deposition", "B": "nudging+dilution", "X": "BLH growth-dil"}
IPR_COLR = {"P": "#1f77b4", "L": "#d62728", "E": "#2ca02c",
            "D": "#9467bd", "B": "#ff7f0e", "X": "#7f7f7f"}


def plot_ipr(sp="O3"):
    p = os.path.join(OUT, "ipr_ledger.csv")
    if not os.path.exists(p):
        return None
    h, R = read_rows(p)
    sps = [r[h.index("species")] for r in R]
    if sp not in sps:
        print(f"  [skip] ipr_ledger 无 {sp}"); return None
    rows = [r for r in R if r[h.index("species")] == sp]
    # hlocal 由 run.py 按 dtout 折算（= 本窗起点当地小时）；h 是**输出步序号**。
    # 首个输出步(h=最小)的窗口跨 spinup 边界，剔除 —— 余下恰为一整日循环。
    hh = [int(r[h.index("h")]) for r in rows]
    keep = [i for i, x in enumerate(hh) if x > min(hh)]
    hr = np.array([float(rows[i][h.index("hlocal")]) for i in keep])
    if len(np.unique(hr)) != len(hr):
        print(f"  [注意] {sp} 的 hlocal 有重复（dtout<3600 时一日有多窗）——"
              f"本图按逐时窗设计，子小时输出请自行聚合")
    o = np.argsort(hr); hr = hr[o]
    vals = {k: np.array([float(rows[i][h.index(k)]) for i in keep])[o] / PPB
            for k in PROCS}
    net = sum(vals.values())

    fig, ax = plt.subplots(1, 3, figsize=(19.2, 5.2))
    # A: 六过程逐时堆叠柱（正过程向上 / 负过程向下）+ 净值折线（图 4 同构）
    w = 0.8
    bup = np.zeros(len(hr)); bdn = np.zeros(len(hr))
    for k in PROCS:
        vp, vn = np.maximum(vals[k], 0), np.minimum(vals[k], 0)
        ax[0].bar(hr, vp, w, bottom=bup, color=IPR_COLR[k], edgecolor="k", lw=.4)
        bup = bup + vp
        ax[0].bar(hr, vn, w, bottom=bdn, color=IPR_COLR[k], edgecolor="k", lw=.4)
        bdn = bdn + vn
    ax[0].plot(hr, net, "k-o", ms=3.5, lw=1.6, zorder=5)
    ax[0].axhline(0, color="k", lw=.8)
    ax[0].set_xticks(range(24))
    ax[0].set_xticklabels([f"{i:02d}" for i in range(24)], fontsize=8)
    ax[0].set_xlabel("hour of window start")
    ax[0].set_ylabel("process integral [ppb/h]  (1-h backward window)")
    ax[0].set_title(f"({sp}) six-process hourly budget — up: P/E/B+, down: L/D/X−")
    ax[0].legend(handles=[Patch(facecolor=IPR_COLR[k], edgecolor="k",
                                label=f"{k} {IPR_NAME[k]}") for k in PROCS]
                          + [Line2D([0], [0], color="k", marker="o", label="net")],
                 fontsize=7, ncol=2, loc="upper left")
    ax[0].grid(axis="y", alpha=.2)
    # B: 六条过程线（symlog）
    for k in PROCS:
        ax[1].plot(hr, vals[k], "-o", ms=2.5, lw=1.3, color=IPR_COLR[k], label=k)
    ax[1].set_yscale("symlog", linthresh=1.0)
    ax[1].set_xlabel("hour of window start"); ax[1].set_ylabel("[ppb/h]")
    ax[1].set_title(f"({sp}) six accounts by process")
    ax[1].legend(fontsize=7, ncol=2); ax[1].grid(alpha=.25)
    # C: 净值 = 六账合计（≈ ΔC(h−1h→h)，后向窗闭合 ≤2%，见 README §6.4）
    ax[2].plot(hr, net, "k-o", ms=3, lw=1.6)
    ax[2].axhline(0, color="grey", lw=.8)
    ax[2].set_xlabel("hour of window start"); ax[2].set_ylabel("net [ppb/h]")
    ax[2].set_title(f"({sp}) net = P+L+E+D+B+X ≈ ΔC (closure ≤2%)")
    ax[2].grid(alpha=.25)
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


# ── 图 4 风格：O₃ 收支堆叠柱状图（Eq4–6 + P(O₃) 折线，0–23 逐时）──
def plot_o3_budget():
    p = os.path.join(OUT, "p_o3no_baseline.csv")
    if not os.path.exists(p):
        return None
    h, R = read_rows(p)
    hr = np.array([float(r[h.index("hour")]) for r in R])
    G1 = fcol(h, R, "G_RO2NO")           # RO₂+NO（白）
    G2 = fcol(h, R, "G_HO2NO")           # HO₂+NO（灰）
    # D(O₃) 五项（Eq6），全部为正数（绘制时画在负方向）
    D1 = fcol(h, R, "D_HO2_O3")          # HO₂+O₃（深灰）
    D2 = fcol(h, R, "D_OH_O3")           # OH+O₃（中灰）
    D3 = fcol(h, R, "D_OH_NO2")          # OH+NO₂（浅灰）
    D4 = fcol(h, R, "D_O1D_H2O")         # O¹D+H₂O（淡灰）
    D5 = fcol(h, R, "D_alkene_O3")       # 烯烃+O₃（淡）
    P  = fcol(h, R, "P_paper")           # P(O₃) 折线

    fig, ax = plt.subplots(figsize=(13.2, 6.2))
    w = 0.7
    # G 堆叠（正向）：RO₂+NO 在最外，HO₂+NO 在内
    ax.bar(hr, G1, w, color="white", edgecolor="k", lw=0.6, label="RO$_2$+NO→NO$_2$+RO$_2$")
    ax.bar(hr, G2, w, bottom=G1, color="#b0b0b0", edgecolor="k", lw=0.6,
           label="HO$_2$+NO→NO$_2$+OH")
    # D 堆叠（负向）：按 Eq6 顺序
    ax.bar(hr, -D1, w, color="#404040", edgecolor="k", lw=0.6, label="HO$_2$+O$_3$→OH+O$_2$")
    ax.bar(hr, -D2, w, bottom=-D1, color="#808080", edgecolor="k", lw=0.6,
           label="OH+O$_3$→HO$_2$+O$_2$")
    ax.bar(hr, -D3, w, bottom=-(D1+D2), color="#c0c0c0", edgecolor="k", lw=0.6,
           label="NO$_2$+OH→HNO$_3$")
    ax.bar(hr, -D4, w, bottom=-(D1+D2+D3), color="#e0e0e0", edgecolor="k", lw=0.6,
           label="O$^1$D+H$_2$O→2OH")
    ax.bar(hr, -D5, w, bottom=-(D1+D2+D3+D4), color="#f0f0f0", edgecolor="k", lw=0.6,
           label="alkenes+O$_3$")
    # P(O₃) 折线
    ax.plot(hr, P, "k-o", ms=4, lw=1.8, label="P(O$_3$) = G−D")
    ax.axhline(0, color="k", lw=.8)
    ax.set_xlabel("hour of day", fontsize=11)
    ax.set_ylabel("D(O$_3$) and G(O$_3$) [ppb/h]", fontsize=11)
    ax.set_title("O$_3$ budget — Eq(4)–(6) stacked bars (G up / D down) + P(O$_3$) line\n"
                 "0–23 h hourly, 4-day case (post-spinup day)", fontsize=11.5)
    ax.legend(fontsize=7, ncol=2, loc="upper left", framealpha=0.9)
    ax.grid(axis="y", alpha=.2)
    ax.set_xticks(range(24))
    ax.set_xticklabels([f"{v:02d}" for v in range(24)], fontsize=8)
    out = os.path.join(OUT, "case_o3_budget.png")
    fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    return out

# ── 图 5 风格：前体物 RIR 分组柱图（2h 时窗 × 类别，0-23h 完整逐时）──
def plot_rir_fig5():
    pd_ = os.path.join(OUT, "rir_diurnal.csv")
    if not os.path.exists(pd_):
        return None
    h, R = read_rows(pd_)
    groups = sorted({r[h.index("group")] for r in R})
    hours = sorted({int(float(r[h.index("hour")])) for r in R})
    dat = {}
    for r in R:
        key = (r[h.index("group")], int(float(r[h.index("hour")])))
        dat[key] = float(r[h.index("RIR_avg")])
    # 类别颜色映射（图 5 风格）
    colors = {"NOx": "#1a1a1a", "TVOC": "#2c3e50", "CO": "#7f8c8d",
              "ALK": "#e74c3c", "OLE": "#e67e22", "ARO": "#f1c40f",
              "ALD": "#2ecc71", "KET": "#1abc9c", "OVOC": "#3498db",
              "ISOP": "#9b59b6", "TERP": "#34495e", "DMS": "#95a5a6"}
    fig, ax = plt.subplots(figsize=(14.8, 6.4))
    w = 0.8 / len(groups)
    for gi, g in enumerate(groups):
        xs = [hh + gi * w for hh in hours]
        vs = [dat.get((g, hh), 0) for hh in hours]
        ax.bar(xs, vs, w, label=g, color=colors.get(g, "#bdc3c7"),
               edgecolor="grey", lw=0.4)
    ax.axhline(0, color="k", lw=.8)
    ax.set_xlabel("hour of day", fontsize=11)
    ax.set_ylabel("RIR [%]", fontsize=11)
    ax.set_title("Figure-5 style: precursor RIR by group, hourly 0–23\n"
                 "(night hours auto-excluded: P(O$_3$) ≈ 0 → RIR undefined)",
                 fontsize=11.5)
    ax.legend(fontsize=7.5, ncol=4, loc="upper left")
    ax.grid(axis="y", alpha=.2)
    ax.set_xticks(hours)
    out = os.path.join(OUT, "case_rir_fig5.png")
    fig.tight_layout(); fig.savefig(out, dpi=160, bbox_inches="tight"); plt.close(fig)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="O3")
    a = ap.parse_args()
    made = [p for p in (plot_ipr(a.species), plot_irr(), plot_rir(),
                        plot_o3_budget(), plot_rir_fig5()) if p]
    if not made:
        print("✗ output/ 无可绘 CSV —— 先跑 python CASE-analysis/run.py", file=sys.stderr)
        return 1
    for p in made:
        print(os.path.relpath(p, ROOT), "完成")
    return 0


if __name__ == "__main__":
    sys.exit(main())
