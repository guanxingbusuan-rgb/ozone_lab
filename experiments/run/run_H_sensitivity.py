#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_H_sensitivity.py —— 实验 H：OAT 敏感性（**不含伴随反演**，无绘图）

对 9 个排放物种各 ±δ 扰动（中心差分），构造雅可比 J（观测向量 × 参数），
并给出每个参数的**归一化灵敏度** S = (ΔO3max/O3max)/(ΔE/E)。

⚠ 本脚本**不做反演**（SVD 反演另属伴随/反演类，未包含）。
   若要反演，见 experiments/sens_adj.py ②段。

输出:
  output/run_H/jacobian.csv         J: 观测维 × 参数（中心差分）
  output/run_H/sensitivity.csv      参数 × {基线值, 扰动值, ΔO3max, 归一化灵敏度}
  output/run_H/oat_runs.csv         每次扰动的运行记录（rc/nst/耗时）

用法:
  python experiments/run/run_H_sensitivity.py          # 默认 ±10%
  python experiments/run/run_H_sensitivity.py 5        # ±5%
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (ROOT, require_prereq, out_dir, run_fac_v3, clone_input,
                     read_table, col, write_csv, manifest, summary_table, log)

# OAT 扰动的排放物种（覆盖主要源类）
# V2.4: 原 9 个 + 新增 10 个排放类（长链烷烃/乙烯/内烯烃/高级醛/甲酚/倍半萜/醇/醚/酯）
# 分组含义见 input_v3/emissions.config 与 templates/SPECIES_INPUT_MATRIX.md
EMIT = [
    # 无机 + 一次醛
    "NO", "NO2", "CO", "HCHO",
    # 生物源
    "ISOP", "TERP", "SESQ",
    # V2.2 集总类
    "ALK", "OLE", "ARO1", "ARO2", "KET",
    # V2.4 新增类
    "LALK", "ETH", "IOLE", "ALD2", "CSL", "MEOH", "ETOH", "MTBE", "ESTER",
]
OBS_SPECIES = ["O3", "HCHO", "NO2", "CO"]
OBS_TIMES_H = list(range(96, 120))


def obs_vec(pref):
    """观测向量: 4 物种 × 24 时刻（ppb），按 CSV 行对齐"""
    hdr, data = read_table(f"{pref}_conc.csv")
    h = col(hdr, data, "h")
    y = []
    for sp in OBS_SPECIES:
        v = col(hdr, data, sp) / 2.46e10
        for t in OBS_TIMES_H:
            y.append(float(v[int(np.argmin(np.abs(h - t)))]))
    return np.array(y)


def main(argv):
    pct = float(argv[1]) if len(argv) > 1 else 10.0
    f_up, f_dn = 1.0 + pct / 100.0, 1.0 - pct / 100.0
    log(f"═══ 实验 H：OAT 敏感性（±{pct:g}%，{len(EMIT)} 个排放物种）═══")
    log("  ⚠ 本实验只做 OAT 梯度，不含伴随/反演")
    require_prereq()
    d = out_dir("run_H")

    # 基线
    log("\n[基线] …")
    rt0 = clone_input(os.path.join(d, "in_base"))
    b = run_fac_v3("SCEN_BAU", os.path.join(d, "base"), rt=rt0, verbose=True)
    y0 = obs_vec(os.path.join(d, "base"))
    o3_base = b["o3max"]

    # OAT 中心差分
    log(f"\n[OAT ± {pct:g}%] …")
    J = np.zeros((len(y0), len(EMIT)))
    rows, runs = [], []
    for j, sp in enumerate(EMIT):
        rp = run_fac_v3("SCEN_BAU", os.path.join(d, f"p_{sp}"),
                        rt=clone_input(os.path.join(d, f"in_p_{sp}"),
                                       emissions_scale={sp: f_up}))
        rm = run_fac_v3("SCEN_BAU", os.path.join(d, f"m_{sp}"),
                        rt=clone_input(os.path.join(d, f"in_m_{sp}"),
                                       emissions_scale={sp: f_dn}))
        yp, ym = obs_vec(os.path.join(d, f"p_{sp}")), obs_vec(os.path.join(d, f"m_{sp}"))
        J[:, j] = (yp - ym) / (f_up - f_dn)        # ∂y/∂(E/E0)
        # 归一化灵敏度（以 O3max 计）
        d_o3 = rp["o3max"] - rm["o3max"]
        S = (d_o3 / o3_base) / ((f_up - f_dn) / 2.0) / 2.0 if o3_base else np.nan
        rows.append([sp, f"{f_dn:.3f}", f"{f_up:.3f}",
                     f"{rm['o3max']:.3f}", f"{rp['o3max']:.3f}",
                     f"{d_o3:+.3f}", f"{S:+.4f}", f"{np.linalg.norm(J[:, j]):.4f}"])
        runs.append([sp, rp["nst"], f"{rp['elapsed_s']:.2f}",
                     rm["nst"], f"{rm['elapsed_s']:.2f}"])
        log(f"    {sp:5s}: O3max {rm['o3max']:7.2f} → {rp['o3max']:7.2f} "
            f"(Δ={d_o3:+6.2f}, S={S:+.3f}, ||J||={np.linalg.norm(J[:, j]):.2f})")

    p1 = write_csv(os.path.join(d, "jacobian.csv"),
                   [f"obs{i}" for i in range(J.shape[0])],
                   [[f"{x:.6e}" for x in J[:, j]] for j in range(len(EMIT))],
                   note=("实验H | ∂y/∂(E/E0) 中心差分 | 行=参数(" + ",".join(EMIT) +
                         ") 列=观测(4物种×24时刻, ppb)"))
    p2 = write_csv(os.path.join(d, "sensitivity.csv"),
                   ["species", f"E_down(×{f_dn:g})", f"E_up(×{f_up:g})",
                    "O3max_down", "O3max_up", "dO3max", "S_normalized", "norm_J_col"],
                   rows, note=f"实验H | OAT ±{pct:g}% | S=(ΔO3max/O3max)/(ΔE/E)")
    p3 = write_csv(os.path.join(d, "oat_runs.csv"),
                   ["species", "nst_up", "t_up_s", "nst_down", "t_down_s"], runs,
                   note=f"实验H | 每次扰动的运行记录（共 {2 * len(EMIT)} 次 fac_v3）")
    summary_table(f"实验 H · OAT 灵敏度（±{pct:g}%）",
                  ["species", f"O3max(-{pct:g}%)", f"O3max(+{pct:g}%)", "ΔO3max", "S", "||J||"],
                  [[r[0], r[3], r[4], r[5], r[6], r[7]] for r in rows])
    log(f"\n  共 {2 * len(EMIT)} 次运行；总耗时 "
        f"{sum(float(r[2]) + float(r[4]) for r in runs):.1f}s")
    mp = manifest(d, "H_sensitivity",
                  {"perturb_pct": pct, "species": EMIT,
                   "obs_species": OBS_SPECIES, "runs": 2 * len(EMIT)},
                  [p1, p2, p3])
    log(f"\n产物: {os.path.relpath(d, ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
