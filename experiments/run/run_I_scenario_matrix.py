#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_I_scenario_matrix.py —— 实验 I：情景矩阵（多情景 × 多路线，方向一致性，无绘图）

对 6 个政策情景，在 fac_v3（74sp）与 mcm_run（29sp）上分别跑，给出：
  · 各路线 × 情景的 O3max 矩阵
  · 相对 BAU 的 Δ% 矩阵
  · 方向一致性（是否与 fac_v3 同号）

输出:
  output/run_I/matrix_o3max.csv     route × scenario 的 O3max (ppb)
  output/run_I/matrix_delta.csv     route × scenario 的 Δ% (vs 该路线 BAU)
  output/run_I/direction.csv        方向一致性判定表
  output/run_I/mcm/<SCEN>_*.csv     MCM 各情景输出

用法:
  python experiments/run/run_I_scenario_matrix.py
  python experiments/run/run_I_scenario_matrix.py --no-mcm     # 只跑 fac_v3
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _runner import (ROOT, require_prereq, out_dir, run_fac_v3, run_mcm,
                     write_csv, manifest, summary_table, log, last24_mean, SCEN_ALL)

# mcm_run 只支持排放缩放（fnox/fvoc），无 JFAC/DT —— 仅映射可比情景
MCM_MAP = {
    "BAU": (1.0, 1.0),
    "LNNOX": (0.3, 1.0),
    "HVOC": (1.0, 3.0),
    "ISOPX": (0.5, 2.0),
}


def main(argv):
    no_mcm = "--no-mcm" in argv
    log(f"═══ 实验 I：情景矩阵（{len(SCEN_ALL)} 情景 × "
        f"{'1' if no_mcm else '2'} 路线）═══")
    require_prereq(need_mcm=not no_mcm)
    d = out_dir("run_I")

    routes = {"fac_v3 74sp": {}, "mcm 29sp": {}}
    prods = []

    # ── fac_v3：全 6 情景 ──
    log("\n[1] fac_v3 74sp …")
    for s in SCEN_ALL:
        r = run_fac_v3(f"SCEN_{s}", os.path.join(d, f"v3_{s}"), verbose=True)
        routes["fac_v3 74sp"][s] = r["o3max"]
        prods += [os.path.join(d, f"v3_{s}_conc.csv")]

    # ── mcm：仅 4 个可比情景 ──
    if not no_mcm:
        log("\n[2] mcm 29sp（仅 fnox/fvoc 可比情景）…")
        for s, (fn, fv) in MCM_MAP.items():
            pref = os.path.join(d, "mcm", s)
            run_mcm(5, pref, fnox=fn, fvoc=fv)
            o3 = last24_mean(f"{pref}_conc.csv", "O3")
            # mcm conc 的单位是 molec/cm3 → last24_mean 已除 PPB
            routes["mcm 29sp"][s] = o3
            prods += [f"{pref}_conc.csv"]
            log(f"    {s:6s} fnox={fn:.1f} fvoc={fv:.1f} → O3(末24h)={o3:6.2f} ppb")

    # ── 矩阵 ──
    scen_order = [s for s in SCEN_ALL]
    hdr = ["route"] + scen_order
    rows_max, rows_delta, rows_dir = [], [], []
    for rt, vals in routes.items():
        if not vals:
            continue
        base = vals.get("BAU", np.nan)
        rows_max.append([rt] + [f"{vals.get(s, np.nan):.3f}"
                                if s in vals else "—" for s in scen_order])
        rows_delta.append([rt] + [f"{(vals[s] / base - 1) * 100:+.1f}%"
                                  if s in vals and base == base and base else "—"
                                  for s in scen_order])

    ref = routes["fac_v3 74sp"]
    for rt, vals in routes.items():
        if rt == "fac_v3 74sp":
            continue
        row = [rt]
        for s in scen_order:
            if s not in vals or s == "BAU":
                row.append("—")
                continue
            dref = ref.get(s, np.nan) - ref.get("BAU", np.nan)
            dcur = vals[s] - vals.get("BAU", np.nan)
            same = (np.sign(dref) == np.sign(dcur))
            row.append(f"{'✓' if same else '✗'} ({dcur:+.2f} vs {dref:+.2f})")
        rows_dir.append(row)

    p1 = write_csv(os.path.join(d, "matrix_o3max.csv"), hdr, rows_max,
                   note="实验I | O3max [ppb] | fac_v3 为完整机制；mcm 仅支持排放缩放情景")
    p2 = write_csv(os.path.join(d, "matrix_delta.csv"), hdr, rows_delta,
                   note="实验I | 相对本路线 BAU 的 Δ%")
    p3 = write_csv(os.path.join(d, "direction.csv"), ["route"] + scen_order, rows_dir,
                   note="实验I | 方向一致性 (相对 fac_v3): ✓同向 ✗反向；括号内为 Δppb")
    summary_table("实验 I · O3max 矩阵（ppb）", hdr, rows_max)
    summary_table("实验 I · Δ% vs BAU", hdr, rows_delta)
    if rows_dir:
        summary_table("实验 I · 方向一致性（vs fac_v3）", ["route"] + scen_order, rows_dir)

    log("\n  注: mcm_run 只接受 fnox/fvoc → HOT/CLOUD 无法在其上表达（标 —）")
    log("      mcm ISOPX 方向与 fac_v3 相反是真实的跨机制分歧（非 bug）")
    mp = manifest(d, "I_scenario_matrix",
                  {"scenarios": scen_order, "routes": list(routes),
                   "mcm_skip": no_mcm},
                  [p1, p2, p3] + prods)
    log(f"\n产物: {os.path.relpath(d, ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
