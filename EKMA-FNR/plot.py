#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
EKMA-FNR/plot.py —— S2 FNR 响应面可视化（读 output/ekma_fnr_grid.csv）

轴序与 S1（EKMA-emission）一致：VOC 在横轴、NOx 在纵轴。
图面（IC 扰动口径，轴 = 控制变量 → 铺满）:
  横轴 = 初始 HCHO 扰动 [%]（0 = config_override 的 HCHO 基线）—— 分子
  纵轴 = 初始 NO₂ 扰动 [%]（0 = config_override 的 NO₂ 基线）—— 分母
         ic1 是 NOx 组：NO 与 NO₂ 同因子缩放，NO₂:NO 比例恒定（图上标出实际比例）
  基线 ppb 数值从 output/fnr_ic_base.csv 读取，换基线后标签自动跟随
  范围/步长由 run.py 的 --pmax/--step 决定（缺省 0–1000% / 100%）
  色标 = O3max [ppb]；黑实线 = O3max 等值线；白虚线 = FNR 等值线（1/2/4，Sillman）
  白点虚线 = 屋脊（每列 HCHO 的 O₃ 峰对应 NO₂ 扰动）

产物（output/）:
  ekma_fnr_S2.png              主图
  ekma_fnr_slices.png          O₃max–NO₂ 切片线图

屋脊轨迹 ekma_fnr_ridge.csv、派生表 ekma_fnr_derived.csv 等**全部 CSV 由 run.py 写出**；
本脚本只读 output/ 的 CSV、只写 PNG。站点定位模式（--obs-single / --obs-more）
的候选带与逐站明细**只打印到终端**，不落 CSV。

用法:
  python plot.py                                              （在 EKMA-FNR/ 目录内）
  python plot.py --obs-single "HCHO=8.0,NO2=4.5,O3=30"        # 单站点定位（ppb）
  python plot.py --obs-more stations.csv                      # 多站点定位（CSV: station,HCHO,NO2,O3）
      观测口径: O₃=日最大, HCHO/NO₂=日均（与模拟窗口均值同口径）
      输出: 各站最佳匹配坐标 (pert_hcho, pert_no2) + 体制判定 + 候选带（终端表）
      图面标记三档: 1 站 = 红星 + 8 近邻候选**包络框+叉点** + 坐标/体制标注;
                   2–5 站 = 红星 + 候选叉点 + 编号: 站名 (体制); >5 站 = 红星 + 只标编号
                   （三档的完整候选明细始终打印在终端；与 EKMA-emission 同规则）
      ⚠ 标记直接画进 ekma_fnr_S2.png（覆盖同名图）；要无标记版就重跑不带 --obs-* 的本脚本
"""
import os
import sys

import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "output")

sys.path.insert(0, os.path.join(ROOT, "analysis_lib"))
from csvio import read_table, read_rows, fcol  # noqa: E402


def _parse_obs(s):
    """'HCHO=8.0,NO2=4.5,O3=30' → {'HCHO':8.0,'NO2':4.5,'O3':30.0}（ppb）"""
    out = {}
    for tok in s.split(","):
        k, _, v = tok.partition("=")
        out[k.strip().upper()] = float(v)
    need = {"HCHO", "NO2", "O3"}
    miss = need - set(out)
    if miss:
        raise SystemExit(f"✗ --obs-single 缺 {sorted(miss)}（格式: HCHO=ppb,NO2=ppb,O3=ppb）")
    extra = set(out) - need
    if extra:
        print(f"  ⚠ --obs-single 忽略未支持的键: {sorted(extra)}（当前匹配空间仅 HCHO/NO₂/O₃）")
    return {k: out[k] for k in need}


def _regime(fnr):
    if fnr != fnr:
        return "n/a"
    return "VOC-limited" if fnr < 1.0 else ("transition" if fnr <= 2.0 else "NOx-limited")


def locate_station(obs, grid_csv, state_csv):
    """状态匹配定位: (HCHO,NO2,O₃)_obs → 最近模拟情景 → 站点坐标与体制。

    匹配空间 = 每案例的窗口状态量 (NO2mean, HCHOmean, O3max)——三者都是可观测口径
    （NO₂/HCHO 日均、O₃ 日最大）。返回 dict（最佳匹配/候选带/越界检查）。"""
    hg, G = read_table(grid_csv)
    hs, S = read_table(state_csv)
    if len(G) != len(S):
        raise SystemExit("✗ grid 与 state 行数不一致（行序对齐失败）")
    no2m = fcol(hs, S, "NO2mean_ppb"); hchom = fcol(hs, S, "HCHOmean_ppb")
    o3m = fcol(hs, S, "O3max_ppb")
    pn = (fcol(hs, S, "fnox") - 1.0) * 100.0          # NO₂ 轴
    ph = (fcol(hs, S, "fvoc") - 1.0) * 100.0          # HCHO 轴
    finit = fcol(hg, G, "FNR_init"); fmean = fcol(hg, G, "FNR_mean")
    o3g = fcol(hg, G, "O3max_ppb")

    q = np.array([obs["NO2"], obs["HCHO"], obs["O3"]])
    cols = np.c_[no2m, hchom, o3m]
    in_range = [bool(cols[:, k].min() <= q[k] <= cols[:, k].max()) for k in range(3)]
    M, SD = cols.mean(0), cols.std(0)
    P = (cols - M) / SD
    qq = (q - M) / SD
    d = np.linalg.norm(P - qq, axis=1)
    k = min(8, len(d))
    idx = np.argsort(d)[:k]
    best = int(idx[0])

    # 屋脊侧判读: 最佳匹配 HCHO 列的屋脊 NO₂
    rh, RD = read_rows(os.path.join(OUT, "ekma_fnr_ridge.csv"))
    rx = fcol(rh, RD, "pert_hcho_pct"); rn = fcol(rh, RD, "ridge_pert_no2_pct")
    ridge_at = float(np.interp(ph[best], rx, rn)) if len(rx) > 1 else float("nan")

    cands = sorted([(pn[i], ph[i], finit[i], fmean[i], o3g[i], d[i]) for i in idx])
    return {"obs": obs, "dist": float(d[best]), "best": best,
            "no2_pct": float(pn[best]), "hcho_pct": float(ph[best]),
            "fnr_init": float(finit[best]), "fnr_mean": float(fmean[best]),
            "o3_sim": float(o3g[best]), "ridge_no2": ridge_at,
            "in_range": in_range, "cands": cands, "k": k}


def report_station(loc):
    o = loc["obs"]
    print("  ── 站点定位（状态匹配: HCHO/NO₂/O₃ 三观测 → 最近模拟情景）──")
    print(f"  观测: HCHO={o['HCHO']:g} ppb(日均)  NO₂={o['NO2']:g} ppb(日均)  "
          f"O₃={o['O3']:g} ppb(日最大)")
    labels = ["NO₂", "HCHO", "O₃"]
    for name, ok in zip(labels, loc["in_range"]):
        print(f"    范围检查 {name}: {'✓ 在模拟状态包络内' if ok else '⚠ 超出包络 → 外推，结果不可靠'}")
    print(f"  匹配距离: {loc['dist']:.3f}σ（越小 = 该实验越能解释你的观测日）")
    b = loc
    print(f"  最佳匹配情景: 初始 HCHO {b['hcho_pct']:+.0f}% / 初始 NOx {b['no2_pct']:+.0f}%"
          f"（相对基线）→ 模拟 O₃max={b['o3_sim']:.1f} ppb")
    print(f"  体制判定（三重校验）:")
    print(f"    ① FNR_init = {b['fnr_init']:.2f} → {_regime(b['fnr_init'])}"
          f"（初始组成比，对应图上 FNR 等值线口径的初始版）")
    print(f"    ② FNR_mean = {b['fnr_mean']:.2f} → {_regime(b['fnr_mean'])}（窗口均浓度比）")
    ridge = b["ridge_no2"]
    if ridge == ridge:
        side = "屋脊左侧 → NOx-limited（加 NO₂ 增 O₃）" if b["no2_pct"] < ridge else                "屋脊右侧 → 滴定区（加 NO₂ 降 O₃）"
        print(f"    ③ 屋脊侧   = 该 HCHO 列屋脊在 NO₂ {ridge:+.0f}%，站点 {b['no2_pct']:+.0f}% → {side}")
    regs = sorted({_regime(c[2]) for c in loc["cands"]})
    if len(regs) > 1:
        print(f"  ⚠ 候选带横跨多体制 {regs} → 定位不确定度大，请报告区间而非单点")
    else:
        print(f"  ⚠ 候选带体制一致（{regs[0]}）→ 定位稳健")


def _annot_off(ax, x, y, dx=7.0, dy=7.0):
    """标注偏移避让：右半区向左标（防压右下角图例），底部抬高（防压 x 轴）"""
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    fx = (x - x0) / max(x1 - x0, 1e-9)
    fy = (y - y0) / max(y1 - y0, 1e-9)
    ox = -dx if fx > 0.55 else dx
    oy = dy * 3.0 if fy < 0.22 else dy
    return ox, oy, ("right" if fx > 0.55 else "left")


def draw_station(ax, loc, many, q, solo):
    """站点标记三档（与 EKMA-emission 同规则）：
    1 站 = 红星 + 8 近邻候选**包络框** + 全标注; 2–5 站 = 红星 + 候选**叉点** + 紧凑标注;
    >5 站 = 红星 + 只标编号（明细在终端候选表）。"""
    b = loc
    nm = "" if b["station"] == "single" else b["station"]
    cx = [c[1] for c in b["cands"]]          # HCHO% 轴
    cy = [c[0] for c in b["cands"]]          # NO₂% 轴
    ax.scatter([b["hcho_pct"]], [b["no2_pct"]], marker="*", s=320 if solo else 260,
               color="red", edgecolors="k", linewidths=0.8, zorder=8,
               label="station (matched)" if q == 1 else None)
    if many:
        ax.annotate(str(q), (b["hcho_pct"], b["no2_pct"]),
                    textcoords="offset points", xytext=(6, 6),
                    fontsize=9, color="k", zorder=9)
        return
    ax.scatter(cx, cy, marker="x", s=48, color="red", linewidths=1.1, zorder=7,
               label="8-NN candidates" if q == 1 else None)
    if not solo:
        ox, oy, cha = _annot_off(ax, b["hcho_pct"], b["no2_pct"])
        ax.annotate(f"{q}: {nm} ({_regime(b['fnr_init'])})",
                    (b["hcho_pct"], b["no2_pct"]),
                    textcoords="offset points", xytext=(ox, oy),
                    fontsize=8.5, color="k", zorder=9, ha=cha,
                    bbox=dict(boxstyle="round,pad=0.25", fc="w", alpha=0.85, ec="red"))
        return
    if len(cx) > 1 and (max(cx) > min(cx) or max(cy) > min(cy)):
        ax.add_patch(plt.Rectangle((min(cx), min(cy)), max(cx) - min(cx),
                                   max(cy) - min(cy), fill=False,
                                   edgecolor="red", lw=1.0, ls=":", zorder=6))
    ax.annotate(f"station: HCHO {b['hcho_pct']:+.0f}% / NO₂ {b['no2_pct']:+.0f}%\n"
                f"regime: {_regime(b['fnr_init'])} (FNR_init={b['fnr_init']:.2f})",
                xy=(b["hcho_pct"], b["no2_pct"]),
                xytext=(0.02, 0.02), textcoords="axes fraction",
                fontsize=9, color="k", ha="left", va="bottom", zorder=9,
                bbox=dict(boxstyle="round,pad=0.3", fc="w", alpha=0.85, ec="red"))


def main():
    ap = argparse.ArgumentParser(description="S2 FNR 响应面可视化（含可选站点定位）")
    ap.add_argument("--obs-single", default=None, metavar='"HCHO=a,NO2=b,O3=c"',
                    help='单站点定位: 三观测 ppb（O₃=日最大, HCHO/NO₂=日均）')
    ap.add_argument("--obs-more", default=None, metavar="stations.csv",
                    help="多站点定位: CSV 列 station,HCHO,NO2,O3（一行一站, ppb）")
    args = ap.parse_args()
    if args.obs_single and args.obs_more:
        print("✗ --obs-single 与 --obs-more 二选一（一次只跑一种站点输入）", file=sys.stderr)
        return 2
    src = os.path.join(OUT, "ekma_fnr_grid.csv")
    if not os.path.exists(src):
        print(f"✗ 缺 {os.path.relpath(src, ROOT)} —— 先跑 python EKMA-FNR/run.py",
              file=sys.stderr)
        return 1
    hdr, D = read_table(src)
    stations = []                                    # [(name, obs_dict), ...]
    if args.obs_more:
        if not os.path.exists(args.obs_more):
            print(f"✗ 站点文件不存在: {args.obs_more}", file=sys.stderr)
            return 2
        import csv as _csv
        _csv.field_size_limit(16 * 1024 * 1024)
        rdr = _csv.DictReader(open(args.obs_more, encoding="utf-8-sig"))
        need = {"station", "hcho", "no2", "o3"}
        if not need <= {k.strip().lower() for k in (rdr.fieldnames or [])}:
            print(f"✗ {args.obs_more} 表头须含 station,HCHO,NO2,O3", file=sys.stderr)
            return 2
        for row in rdr:
            r = {k.strip().lower(): v.strip() for k, v in row.items()}
            if not r.get("station"):
                continue
            try:
                obs = {"HCHO": float(r["hcho"]), "NO2": float(r["no2"]),
                       "O3": float(r["o3"])}
            except (ValueError, KeyError):
                print(f"  ⚠ 跳过站点 {r.get('station')!r}: HCHO/NO2/O3 有非数值", file=sys.stderr)
                continue
            stations.append((r["station"], obs))
        if not stations:
            print("✗ 站点文件无有效行", file=sys.stderr)
            return 2
    if args.obs_single:
        stations.append(("single", _parse_obs(args.obs_single)))
    locs = [dict(station=nm, **locate_station(obs, src,
              os.path.join(OUT, "ekma_fnr_grid_raw.csv_state.csv")))
            for nm, obs in stations]
    bad = [nm for nm, _ in stations if not nm.isascii()]
    if bad:
        print(f"  ⚠ 站名含非 ASCII {bad} → 图上会显示缺字方块（图件不配 CJK 字体）；"
              f"建议站名用 ASCII，中文明细仍在终端表", file=sys.stderr)
    x = fcol(hdr, D, "pert_hcho_pct")     # HCHO ← ic2（分子）
    y = fcol(hdr, D, "pert_no2_pct")      # NO₂  ← ic1（分母，NOx 组）
    z = fcol(hdr, D, "O3max_ppb")
    fr = fcol(hdr, D, "FNR_mean")

    xs = np.array(sorted(set(np.round(x, 4))))   # HCHO 轴
    ys = np.array(sorted(set(np.round(y, 4))))   # NO₂ 轴
    if len(xs) * len(ys) != len(x):
        print(f"✗ 不是完整正交网格（{len(xs)}×{len(ys)} vs {len(x)} 行）", file=sys.stderr)
        return 1
    # 键序必须与 x/y 一致: x=pert_hcho, y=pert_no2 → idx[(hcho, no2)]
    idx = {(round(a, 4), round(b, 4)): k for k, (a, b) in enumerate(zip(x, y))}
    Z = np.array([[z[idx[(a, b)]] for b in ys] for a in xs])   # Z[hcho_idx, no2_idx]
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
    for j, xv in enumerate(xs):
        col = Z[j, :]
        if np.isfinite(col).any():
            k = int(np.nanargmax(col))
            ridge_no2.append(ys[k]); ridge_o3.append(col[k])
    # ── 主图（x = 初始 HCHO 扰动 ← 分子, y = 初始 NO₂ 扰动 ← 分母）──
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

    ax.plot(xs, ridge_no2, "w--", lw=2.4, label="O3 ridge")
    ax.plot(xs, ridge_no2, "k-", lw=1.2, label=None)
    ax.scatter([0], [0], marker="+", s=180, color="k", zorder=6,
               label="baseline (0%, 0%)")

    imax = int(np.nanargmax(Z))
    jj, kk = np.unravel_index(imax, Z.shape)
    ax.annotate(f"max O3 = {Z[jj, kk]:.1f} ppb\n(HCHO +{xs[jj]:.0f}%, NO₂ +{ys[kk]:.0f}%)",
                xy=(xs[jj], ys[kk]), xytext=(0.98, 0.97), textcoords="axes fraction",
                fontsize=9, color="w", ha="right", va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="k", alpha=0.35, ec="none"))

    ax.set_xlabel(f"initial HCHO perturbation [%]  (baseline {ppb(b_hcho)}, 0 = no perturbation)",
                  fontsize=11)
    ax.set_ylabel(f"initial NO₂ perturbation [%]  (baseline {ppb(b_no2)}, 0 = no perturbation)",
                  fontsize=11)
    ridge_note = ""
    if ridge_no2 and max(ridge_no2) >= ys[-1] - 1e-6:
        n_edge = sum(1 for v in ridge_no2 if v >= ys[-1] - 1e-6)
        ridge_note = (f"\n! ridge at the +{ys[-1]:g}% NO2 edge for {n_edge}/{len(ridge_no2)} "
                      f"HCHO columns - true peak may lie beyond the sampled range")
    ax.set_title("S2 · EKMA-FNR response surface — O3max [ppb]\n"
                 f"{len(xs)}×{len(ys)} = {len(x)} orthogonal IC cases "
                 f"(0–{xs[-1]:g}% in {xs[1]-xs[0]:g}% steps)\n"
                 "fixed: mechanism V2.5+ · met/IC from config_override/ (effective input: "
                 "output/runtime_fnr.dat, see run_meta.txt) · "
                 f"emissions OFF · NO scales with the NOx group "
                 f"(NO2:NO = {ratio}:1 fixed) · O₃ not perturbed"
                 + ridge_note,
                 fontsize=11)
    if locs:
        many, solo = len(locs) > 5, len(locs) == 1
        for q, L in enumerate(locs, 1):
            draw_station(ax, L, many, q, solo)
    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax.grid(alpha=0.15, color="w")

    # ── 切片图: O₃max vs NO₂（选定 HCHO 水平）──
    fig2, ax2 = plt.subplots(figsize=(9.2, 6.2))
    picks = [v for v in (0, 100, 200, 300, 400) if v in xs]
    cmap = plt.get_cmap("viridis")
    for q, v in enumerate(picks):
        j = int(np.where(xs == v)[0][0])
        colv = Z[j, :]
        ax2.plot(ys, colv, "-o", ms=2.6, lw=1.6,
                 color=cmap(q / max(len(picks) - 1, 1)), label=f"HCHO +{v:g}%")
        k = int(np.nanargmax(colv))
        ax2.plot(ys[k], colv[k], "*", ms=11,
                 color=cmap(q / max(len(picks) - 1, 1)), mec="k", mew=0.6)
    ax2.set_xlabel(f"initial NO₂ perturbation [%]  (baseline {ppb(b_no2)})", fontsize=11)
    ax2.set_ylabel("O3max [ppb]", fontsize=11)
    ax2.set_title("S2 slices · O3max vs initial NO₂ at selected HCHO levels\n"
                  "star = local peak (ridge point)", fontsize=11)
    ax2.grid(alpha=0.3); ax2.legend(fontsize=9)
    png2 = os.path.join(OUT, "ekma_fnr_slices.png")
    fig2.tight_layout(); fig2.savefig(png2, dpi=170, bbox_inches="tight")
    plt.close(fig2)

    png = os.path.join(OUT, "ekma_fnr_S2.png")
    fig.tight_layout()
    fig.savefig(png, dpi=170, bbox_inches="tight")
    plt.close(fig)

    print(f"  O3max {zlo:.2f} – {zhi:.2f} ppb")
    print(f"  产物: {os.path.relpath(png, ROOT)}")
    print(f"        {os.path.relpath(png2, ROOT)}")
    if locs:
        # ── 站点状态匹配候选：全部走终端表（plot.py 只出 PNG，不写 CSV）──
        note_st = "; ".join(f"{nm}(HCHO={o['HCHO']:g},NO2={o['NO2']:g},O3={o['O3']:g})"
                            for nm, o in stations)
        print(f"  ── 匹配空间 = NO2mean/HCHOmean/O3max 标准化欧氏距离；"
              f"regime 按 FNR_init: <1 VOC-limited / 1–2 transition / >2 NOx-limited ──")
        print(f"  观测状态: {note_st}")
        for L in locs:
            head = L["station"] if len(locs) > 1 else "obs"
            print(f"\n  [{head}] 最近匹配距离 {L['dist']:.2f} σ")
            print(f"    {'rank':>4s} {'距离σ':>7s} {'NO₂%':>7s} {'HCHO%':>7s} "
                  f"{'FNR_init':>9s} {'FNR_mean':>9s} {'O3max':>8s}  regime")
            for i, c in enumerate(L["cands"], 1):
                print(f"    {i:>4d} {c[5]:7.4f} {c[0]:7.1f} {c[1]:7.1f} "
                      f"{c[2]:9.4f} {c[3]:9.4f} {c[4]:8.3f}  {_regime(c[2])}")
        # ── 终端报告 ──
        if len(locs) == 1:
            report_station(locs[0])
        else:
            print("  ── 多站点定位汇总（状态匹配）──")
            print(f"  {'站点':<12s} {'距离σ':>7s} {'HCHO%':>7s} {'NO₂%':>7s} "
                  f"{'FNR_init':>8s}  体制")
            for L in locs:
                print(f"  {L['station']:<12.12s} {L['dist']:7.2f} "
                      f"{L['hcho_pct']:7.0f} {L['no2_pct']:7.0f} "
                      f"{L['fnr_init']:8.2f}  {_regime(L['fnr_init'])}"
                      + ("  ⚠ 候选跨多体制" if len({_regime(c[2]) for c in L['cands']}) > 1 else ""))
            regs = {_regime(L['fnr_init']) for L in locs}
            print(f"  体制分布: " + ", ".join(f"{r}×{sum(1 for L in locs if _regime(L['fnr_init'])==r)}"
                                              for r in sorted(regs)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
