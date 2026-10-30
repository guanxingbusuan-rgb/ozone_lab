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

屋脊轨迹 ekma_emission_ridge.csv、派生表 ekma_emission_derived.csv 等**全部 CSV 由 run.py 写出**；
本脚本只读 output/ 的 CSV、只写 PNG。站点反查模式（--obs-single / --obs-more）的
候选带与逐站明细**只打印到终端**，不落 CSV。

用法:
  python plot.py                                                    （在 EKMA-emission/ 目录内）
  python plot.py --obs-single "O3=142,NO2=6.5"                      # 单站点反查（ppb）
  python plot.py --obs-single "O3=142,NO2=6.5,HCHO=8.0" --with-hcho # 三观测反查（VOC 轴更准）
  python plot.py --obs-more stations.csv --with-hcho                # 多站点（CSV: station,O3,NO2[,HCHO]）
      观测口径: O₃=日最大, NO₂/HCHO=窗口均值（与模拟统计窗同口径, 末 24 输出步）
      输出: 反查坐标 (pert_voc_pct, pert_nox_pct) + 屋脊侧判读 + 每站 8 近邻候选（终端表）
      图面标记三档: 1 站 = 红星 + 8 近邻候选**包络框+叉点** + 坐标/屋脊侧/match σ 标注;
                   2–5 站 = 红星 + 候选叉点 + 编号: 站名 (VOC% / NOx%); >5 站 = 红星 + 只标编号
                   （三档的完整候选明细始终打印在终端；与 EKMA-FNR 同规则）
      ⚠ 标记直接画进 ekma_emission_S1.png（覆盖同名图）；要无标记版就重跑不带 --obs-* 的本脚本
      图面标注一律纯 ASCII（DejaVu 无 CJK）；中文明细只走终端
"""
import argparse
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
from csvio import read_table, fcol  # noqa: E402


def _parse_obs(s, with_hcho):
    """'O3=a,NO2=b[,HCHO=c]' → dict（ppb）；键与 --with-hcho 决定的匹配空间一致"""
    out = {}
    for tok in s.split(","):
        k, _, v = tok.partition("=")
        out[k.strip().upper()] = float(v)
    need = {"O3", "NO2"} | ({"HCHO"} if with_hcho else set())
    miss = need - set(out)
    if miss:
        raise SystemExit(f"✗ --obs-single 缺 {sorted(miss)}（当前匹配空间: O3,NO2"
                         + (",HCHO" if with_hcho else "") + "）")
    extra = set(out) - need
    if extra:
        print(f"  ⚠ 忽略未支持的键: {sorted(extra)}")
    return out


def locate_station(obs, grid_csv, state_csv, with_hcho):
    """状态匹配定位: (O₃,NO₂[,HCHO])_obs → 最近模拟情景 → (VOC%, NOx%) 坐标。

    匹配空间 = 每案例的 (NO2mean_ppb, O3max_ppb[, HCHOmean_ppb])——三者皆可观测
    （NO₂ 站点常规、O₃ 站点常规、HCHO 站点 DNPH/卫星柱经 sat2surf）。
    与 FNR §11.7 同构；差异: 本反查的答案落在**排放政策空间**（与 §9.3 清单法同轴可比）。"""
    hg, G = read_table(grid_csv)
    hs, S = read_table(state_csv)
    if len(G) != len(S):
        raise SystemExit("✗ grid 与 state 行数不一致（行序对齐失败）")
    no2m = fcol(hs, S, "NO2mean_ppb"); o3m = fcol(hs, S, "O3max_ppb")
    pv = fcol(hg, G, "pert_voc_pct"); pn = fcol(hg, G, "pert_nox_pct")
    o3g = fcol(hg, G, "O3max_ppb")
    if with_hcho:
        hchom = fcol(hs, S, "HCHOmean_ppb")
        cols = np.c_[no2m, hchom, o3m]
        q_raw = np.array([obs["NO2"], obs["HCHO"], obs["O3"]])
    else:
        cols = np.c_[no2m, o3m]
        q_raw = np.array([obs["NO2"], obs["O3"]])
    in_range = [bool(cols[:, k].min() <= q_raw[k] <= cols[:, k].max())
                for k in range(cols.shape[1])]
    M, SD = cols.mean(0), cols.std(0)
    P = (cols - M) / SD
    qq = (q_raw - M) / SD
    d = np.linalg.norm(P - qq, axis=1)
    kk = min(8, len(d))
    idx = np.argsort(d)[:kk]
    best = int(idx[0])
    rh, RD = read_table(os.path.join(OUT, "ekma_emission_ridge.csv"))
    rx = fcol(rh, RD, "pert_voc_pct"); rn = fcol(rh, RD, "ridge_pert_nox_pct")
    ridge_at = float(np.interp(pv[best], rx, rn)) if len(rx) > 1 else float("nan")
    cands = [(float(d[i]), float(pv[i]), float(pn[i])) for i in idx]   # (距离, VOC%, NOx%)
    return {"obs": obs, "dist": float(d[best]), "best": best,
            "voc_pct": float(pv[best]), "nox_pct": float(pn[best]),
            "o3_sim": float(o3g[best]), "ridge_nox": ridge_at,
            "in_range": in_range, "cands": cands, "with_hcho": with_hcho}


def _regime_side(voc, nox, ridge_nox):
    """屋脊侧判读（左=NOx-limited, 右=滴定区）；ridge_nox 为 nan 时返回未知"""
    if ridge_nox != ridge_nox:
        return "n/a"
    return "NOx-limited（屋脊左侧，加 NOx 增 O₃）" if voc < ridge_nox else            "NOx-saturated（屋脊右侧，加 NOx 降 O₃）"


def report_station(loc):
    o = loc["obs"]; b = loc
    print("  ── 站点定位（状态匹配 → 排放政策空间坐标）──")
    print(f"  观测: O₃={o['O3']:g} ppb(日最大)  NO₂={o['NO2']:g} ppb(均值)"
          + (f"  HCHO={o['HCHO']:g} ppb(均值)" if loc["with_hcho"] else ""))
    dims = ["NO₂", "O₃"] + (["HCHO"] if loc["with_hcho"] else [])
    for name, ok in zip(dims, loc["in_range"]):
        print(f"    范围检查 {name}: " + ("✓ 在模拟包络内" if ok else "⚠ 超出包络 → 外推，不可靠"))
    print(f"  匹配距离: {b['dist']:.3f}σ（越小 = 该实验越能解释你的观测状态）")
    print(f"  反查坐标: NMVOC {b['voc_pct']:+.0f}% / NOx {b['nox_pct']:+.0f}%"
          f"（模拟 O₃max={b['o3_sim']:.1f} ppb）")
    side = _regime_side(b["voc_pct"], b["nox_pct"], b["ridge_nox"])
    print(f"  屋脊侧   : {side}（该 VOC 列屋脊在 NOx {b['ridge_nox']:+.0f}%）")
    if b["with_hcho"]:
        print("  （三观测反查: HCHO 维参与匹配 → VOC 轴歧义已压缩）")
    else:
        print("  ⚠ 仅两观测: VOC 轴歧义 ±40–50%（EKMA 固有）。加 HCHO 观测（--with-hcho）可压缩。")


def _side_short(voc, nox, ridge_nox):
    """屋脊侧的**图用短标签**（纯 ASCII —— 图件不放中文，DejaVu 无 CJK 会成豆腐块）"""
    if ridge_nox != ridge_nox:
        return "ridge n/a"
    return ("ridge-left: NOx-limited" if voc < ridge_nox
            else "ridge-right: NOx-saturated")


def _annot_off(ax, x, y, dx=9.0, dy=9.0):
    """标注偏移避让：右半区向左标（防压右下角图例），底部抬高（防压 x 轴与色标）"""
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    fx = (x - x0) / max(x1 - x0, 1e-9)
    fy = (y - y0) / max(y1 - y0, 1e-9)
    ox = -dx if fx > 0.55 else dx
    oy = dy * 3.0 if fy < 0.22 else dy
    return ox, oy, ("right" if fx > 0.55 else "left")


def draw_station(ax, loc, many, q, solo):
    """站点标记：红星 = 反查落点。分三档（与 EKMA-FNR 同规则）：
    >5 站 → 只标编号；2–5 站 → 编号 + 紧凑标注；单站 → 再加 8 近邻候选带（虚线框）+ 全标注。"""
    b = loc
    nm = "" if b["station"] == "single" else b["station"]
    ax.scatter([b["voc_pct"]], [b["nox_pct"]], marker="*", s=280, color="red",
               edgecolors="k", linewidths=0.8, zorder=8,
               label="station (matched)" if q == 1 else None)
    if many:
        ax.annotate(str(q), (b["voc_pct"], b["nox_pct"]),
                    textcoords="offset points", xytext=(6, 6),
                    fontsize=9, color="k", zorder=9)
        return
    cx = [c[1] for c in b["cands"]]          # VOC% 轴
    cy = [c[2] for c in b["cands"]]          # NOx% 轴
    ax.scatter(cx, cy, marker="x", s=48, color="red", linewidths=1.1, zorder=7,
               label="8-NN candidates" if q == 1 else None)
    if not solo:
        ox, oy, cha = _annot_off(ax, b["voc_pct"], b["nox_pct"])
        ax.annotate(f"{q}: {nm} ({b['voc_pct']:+.0f}% / {b['nox_pct']:+.0f}%)",
                    (b["voc_pct"], b["nox_pct"]),
                    textcoords="offset points", xytext=(ox, oy),
                    fontsize=8.5, color="k", zorder=9, ha=cha,
                    bbox=dict(boxstyle="round,pad=0.25", fc="w", alpha=0.85, ec="red"))
        return
    if len(cx) > 1 and (max(cx) > min(cx) or max(cy) > min(cy)):
        ax.add_patch(plt.Rectangle((min(cx), min(cy)), max(cx) - min(cx),
                                   max(cy) - min(cy), fill=False,
                                   edgecolor="red", lw=1.0, ls=":", zorder=6))
    ax.annotate(f"station: NMVOC {b['voc_pct']:+.0f}% / NOx {b['nox_pct']:+.0f}%\n"
                f"{_side_short(b['voc_pct'], b['nox_pct'], b['ridge_nox'])}"
                f" · match {b['dist']:.2f}σ"
                + ("" if b["with_hcho"] else "\n2-obs: VOC axis ambiguous ±40–50%"),
                xy=(b["voc_pct"], b["nox_pct"]),
                textcoords="offset points", xytext=(10, -34),
                fontsize=8.5, color="k", ha="left", zorder=9,
                bbox=dict(boxstyle="round,pad=0.3", fc="w", alpha=0.85, ec="red"))


def main():
    ap = argparse.ArgumentParser(description="S1 响应面可视化（含可选站点反查定位）")
    ap.add_argument("--obs-single", default=None, metavar='"O3=a,NO2=b[,HCHO=c]"',
                    help="单站点反查: 观测 ppb（O₃=日最大, NO₂=窗口均值; HCHO 可选第三维）")
    ap.add_argument("--obs-more", default=None, metavar="stations.csv",
                    help="多站点反查: CSV 列 station,O3,NO2[,HCHO]（一行一站, ppb）")
    ap.add_argument("--with-hcho", action="store_true",
                    help="把 HCHO 纳入匹配空间（第三维, 破 VOC 轴歧义; 需各站提供 HCHO 观测）")
    args = ap.parse_args()
    if args.obs_single and args.obs_more:
        print("✗ --obs-single 与 --obs-more 二选一（一次只跑一种站点输入）", file=sys.stderr)
        return 2
    src = os.path.join(OUT, "ekma_emission_grid.csv")
    if not os.path.exists(src):
        print(f"✗ 缺 {os.path.relpath(src, ROOT)} —— 先跑 python EKMA-emission/run.py",
              file=sys.stderr)
        return 1
    hdr, D = read_table(src)
    locs = []
    if args.obs_single or args.obs_more:
        state = os.path.join(OUT, "ekma_emission_grid_raw.csv_state.csv")
        if not os.path.exists(state):
            print("✗ 缺状态侧车 _state.csv（先跑 run.py）", file=sys.stderr)
            return 1
        stations = []
        if args.obs_more:
            import csv as _csv
            _csv.field_size_limit(16 * 1024 * 1024)
            rdr = _csv.DictReader(open(args.obs_more, encoding="utf-8-sig"))
            need = {k.strip().lower() for k in (rdr.fieldnames or [])}
            base_need = {"station", "o3", "no2"} | ({"hcho"} if args.with_hcho else set())
            if not base_need <= need:
                print(f"✗ {args.obs_more} 表头缺 {sorted(base_need - need)}", file=sys.stderr)
                return 2
            for row in rdr:
                r = {k.strip().lower(): v.strip() for k, v in row.items()}
                if not r.get("station"):
                    continue
                try:
                    o = {"O3": float(r["o3"]), "NO2": float(r["no2"])}
                    if args.with_hcho:
                        o["HCHO"] = float(r["hcho"])
                except (ValueError, KeyError):
                    print(f"  ⚠ 跳过站点 {r.get('station')!r}: 非数值", file=sys.stderr)
                    continue
                stations.append((r["station"], o))
            if not stations:
                print("✗ 站点文件无有效行", file=sys.stderr)
                return 2
        if args.obs_single:
            stations.append(("single", _parse_obs(args.obs_single, args.with_hcho)))
        locs = [dict(station=nm, **locate_station(o, src, state, args.with_hcho))
                for nm, o in stations]
        bad = [nm for nm, _ in stations if not nm.isascii()]
        if bad:
            print(f"  ⚠ 站名含非 ASCII {bad} → 图上会显示缺字方块（图件不配 CJK 字体）；"
                  f"建议站名用 ASCII，中文明细仍在终端表", file=sys.stderr)
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
    if locs:
        many, solo = len(locs) > 5, len(locs) == 1
        for q, L in enumerate(locs, 1):
            draw_station(ax, L, many, q, solo)
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
    if locs:
        wh = "--with-hcho" if args.with_hcho else ""
        space = ("匹配空间=NO2mean/O3max" + ("+HCHOmean" if args.with_hcho else "") +
                 " 标准化欧氏距离 " + wh + "; 坐标=排放政策空间 VOC%/NOx%")
        obs_txt = "; ".join(f"{nm}(O3={o['O3']:g},NO2={o['NO2']:g})"
                            for nm, o in stations)
        print(f"  ── {space} ──")
        print(f"  观测状态: " + (f"{len(locs)} 站; " if len(locs) > 1 else "obs: ") + obs_txt)
        # ── 终端报告 ──
        if len(locs) == 1:
            report_station(locs[0])
        else:
            print("  ── 多站点反查汇总 ──")
            print(f"  {'站点':<12s} {'距离σ':>7s} {'VOC%':>7s} {'NOx%':>7s}  屋脊侧")
            for L in locs:
                side = _regime_side(L["voc_pct"], L["nox_pct"], L["ridge_nox"])
                print(f"  {L['station']:<12.12s} {L['dist']:7.2f} {L['voc_pct']:7.0f} "
                      f"{L['nox_pct']:7.0f}  {side[:12]}")
        # ── 每站 8 近邻候选（plot.py 只写 PNG，明细不落 CSV）──
        for L in locs:
            head = L["station"] if len(locs) > 1 else "obs"
            print(f"\n  [{head}] 最近匹配 {L['dist']:.2f} σ · 8 近邻候选")
            print(f"    {'rank':>4s} {'距离σ':>7s} {'VOC%':>7s} {'NOx%':>7s}")
            for i, (dd, v, n) in enumerate(L["cands"], 1):
                print(f"    {i:>4d} {dd:7.4f} {v:7.1f} {n:7.1f}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
