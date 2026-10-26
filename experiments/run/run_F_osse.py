#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_F_osse.py —— 实验 F：观测系统模拟实验 OSSE（同化技能评估，无绘图）

双胞胎法：真值运行（排放×1.4）→ 采样加噪成"观测" → 用不同 τ 同化 → 评估 RMSE 技能。

输入假设: runtime_v3.dat 已编译（本脚本会为每个 τ 生成变体 runtime 并调用 cfgbuild）
输出:
  output/run_F/obs_samples.csv      合成观测（时间 × 物种，ppb，含噪声设定）
  output/run_F/skill_rmse.csv       每条通路 × 每个 τ 的 RMSE
  output/run_F/skill_delta.csv      RMSE 相对无约束基线的改善
  output/run_F/<tau>h_{conc,irr,ipr}.csv   各 τ 的同化运行输出

用法:
  python experiments/run/run_F_osse.py
  python experiments/run/run_F_osse.py 0.3       # 改观测噪声 (相对标准差)
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (ROOT, require_prereq, out_dir, run_fac_v3, clone_input,
                     read_table, col, write_csv, manifest, summary_table, log, PY)

OBS_SPECIES = ["O3", "HCHO", "ISOP"]
# ★ 关键: 序列的 t 必须是"自模拟起点秒"。conc CSV 只含末 24h (h=96..120),
#   故绝对小时 = 96..119。写成相对 0..23 会被 ser_interp 当作"全在序列之前",
#   外延钳制到末值 → 同化曲线被拉平（曾踩过）。
OBS_TIMES_H = [96 + i for i in range(24)]          # 绝对小时 96..119（逐时 24 点）
TAUS_H = [2, 6, 24, 72]
TRUTH_SCALE = {"NO": 1.4, "NO2": 1.4, "ISOP": 1.4}  # 真值：NOx/ISOP 排放 ×1.4


def rmse(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    return float(np.sqrt(np.mean((a[m] - b[m]) ** 2))) if m.any() else np.nan


def main(argv):
    noise = float(argv[1]) if len(argv) > 1 else 0.10
    log(f"═══ 实验 F：OSSE 双胞胎同化（噪声 {noise:.0%}）═══")
    require_prereq()
    d = out_dir("run_F")

    # ── ① 真值运行（排放 ×1.4）──
    log("\n[1/4] 真值运行 (NO/NO2/ISOP 排放 ×1.4) …")
    rt_truth = clone_input(os.path.join(d, "in_truth"), emissions_scale=TRUTH_SCALE)
    r = run_fac_v3("SCEN_BAU", os.path.join(d, "truth"), rt=rt_truth, verbose=True)

    hdr, data = read_table(f"{os.path.join(d, 'truth')}_conc.csv")
    h = col(hdr, data, "h")
    truth = {sp: col(hdr, data, sp) / 2.46e10 for sp in OBS_SPECIES}

    # ── ② 采样加噪 → 合成观测 ──
    log(f"\n[2/4] 采样 {len(OBS_TIMES_H)} 个时刻 (噪声 σ={noise:.0%}, seed=42) …")
    rng = np.random.default_rng(42)
    obs = {}
    rows_obs = []
    for sp in OBS_SPECIES:
        vals = []
        for th in OBS_TIMES_H:
            k = int(np.argmin(np.abs(h - th)))
            v = float(truth[sp][k])
            vals.append(v * (1.0 + rng.normal(0, noise)))
        obs[sp] = (OBS_TIMES_H, vals)
        rows_obs.append([sp] + [f"{v:.4f}" for v in vals])
        log(f"    {sp:5s}: 真值 {np.mean([truth[sp][int(np.argmin(np.abs(h - t)))] for t in OBS_TIMES_H]):.3f} ppb"
            f" → 观测均值 {np.mean(vals):.3f} ppb")
    write_csv(os.path.join(d, "obs_samples.csv"),
              ["species"] + [f"t={t}h" for t in OBS_TIMES_H], rows_obs,
              note=f"实验F | 合成观测 (ppb) | 真值排放 ×1.4 | 相对噪声 σ={noise}")

    # ── ③ 无约束基线 ──
    log("\n[3/4] 无约束基线 …")
    rt_free = clone_input(os.path.join(d, "in_free"))
    run_fac_v3("SCEN_BAU", os.path.join(d, "free"), rt=rt_free, verbose=True)
    hdr_f, data_f = read_table(f"{os.path.join(d, 'free')}_conc.csv")
    hf = col(hdr_f, data_f, "h")
    free = {}
    for sp in OBS_SPECIES:
        v = col(hdr_f, data_f, sp) / 2.46e10
        free[sp] = [float(v[int(np.argmin(np.abs(hf - t)))]) for t in OBS_TIMES_H]

    base_rmse = {sp: rmse(free[sp], obs[sp][1]) for sp in OBS_SPECIES}

    # ── ④ τ 扫描同化 ──
    log("\n[4/4] τ 扫描同化 …")
    rows_skill, rows_delta = [], []
    products = []
    for th in TAUS_H:
        log(f"    τ={th:2d}h …")
        rt_c = clone_input(os.path.join(d, f"in_{th}h"), tau=th * 3600.0, obs=obs)
        pref = os.path.join(d, f"{th}h")
        run_fac_v3("SCEN_BAU", pref, rt=rt_c)
        products += [f"{pref}_conc.csv", f"{pref}_irr.csv", f"{pref}_ipr.csv"]

        hdr_c, data_c = read_table(f"{pref}_conc.csv")
        hc = col(hdr_c, data_c, "h")
        row = [th]
        drow = [th]
        for sp in OBS_SPECIES:
            v = col(hdr_c, data_c, sp) / 2.46e10
            sim = [float(v[int(np.argmin(np.abs(hc - t)))]) for t in OBS_TIMES_H]
            e = rmse(sim, obs[sp][1])
            row.append(f"{e:.4f}")
            drow.append(f"{(1 - e / base_rmse[sp]) * 100:.1f}%" if base_rmse[sp] > 0 else "nan")
        rows_skill.append(row)
        rows_delta.append(drow)

    # 无约束行也放进表（τ=∞）
    rows_skill.insert(0, ["free"] + [f"{base_rmse[sp]:.4f}" for sp in OBS_SPECIES])
    rows_delta.insert(0, ["free"] + ["0.0%" for _ in OBS_SPECIES])

    header = ["tau_h"] + OBS_SPECIES
    write_csv(os.path.join(d, "skill_rmse.csv"), header, rows_skill,
              note="实验F | 同化技能 RMSE [ppb] | free=无约束基线")
    write_csv(os.path.join(d, "skill_delta.csv"), header, rows_delta,
              note="实验F | RMSE 相对 free 的改善 (%) | 正值=改善, 负值=劣化")
    summary_table("实验 F · 同化技能 RMSE（ppb）", header, rows_skill)
    summary_table("实验 F · 相对无约束改善（%）", header, rows_delta)

    best = {}
    for sp in OBS_SPECIES:
        cand = [(TAUS_H[i], float(rows_skill[i + 1][1 + OBS_SPECIES.index(sp)]))
                for i in range(len(TAUS_H))]
        best[sp] = min(cand, key=lambda x: x[1])
    log("\n  最优 τ: " + " · ".join(f"{sp}={best[sp][0]}h (RMSE {best[sp][1]:.2f})"
                                   for sp in OBS_SPECIES))
    log("  注: NO2 未纳入站点约束 —— 其受光稳态支配，强 nudging 会破坏 PSS（见 obs_osse.py 通路③）")

    mp = manifest(d, "F_osse", {"noise": noise, "taus_h": TAUS_H,
                                "obs_species": OBS_SPECIES,
                                "truth_scale": TRUTH_SCALE},
                  [os.path.join(d, "obs_samples.csv"),
                   os.path.join(d, "skill_rmse.csv"),
                   os.path.join(d, "skill_delta.csv")] + products)
    log(f"\n产物: {os.path.relpath(d, ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
