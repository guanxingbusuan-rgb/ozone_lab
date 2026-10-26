#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_all.py —— 全部实验运行总驱动（类似 CMAQ 的 run.cctm.csh 角色）

只跑模拟 + 落 CSV（**不含任何绘图、不含伴随反演**）。
对应关系:
  run.cctm.csh   只管"跑 CMAQ 化学输运并写输出"，绘图是另外的 ncl/python 脚本
  run_all.py     只管"跑 ozone_lab 各项模拟并写 output/run_X/*.csv"，绘图见 figures/

用法:
  python experiments/run/run_all.py                # 全跑（A–I，除反演）
  python experiments/run/run_all.py A B D          # 只跑指定实验
  python experiments/run/run_all.py --list         # 列出实验清单
  OZONE_OUT=/data/out python .../run_all.py        # 换输出根

输出总览:
  output/run_A/  六情景单箱（summary.csv + 各情景三件 CSV）
  output/run_B/  响应面（batch.csv + 19 个 indicator_matrix_*.csv + slices.csv）
  output/run_C/  跨路线屋脊（ridge_batch/mcm/compare.csv）
  output/run_D/  IPR/IRR 过程诊断（ipr_ledger/closure/irr_top.csv）
  output/run_E/  指示因子（12 指标 + 对比）
  output/run_F/  OSSE（obs_samples/skill_rmse/skill_delta.csv）
  output/run_G/  逐小时过程（hourly_ledger/share/rxn_top/consistency.csv）
  output/run_H/  OAT 敏感性（jacobian/sensitivity.csv）
  output/run_I/  情景矩阵（matrix_o3max/delta/direction.csv）
  output/run_M/  气象序列驱动（summary/attribution/o3_timeseries.csv）
  output/run_all_manifest.json  全局清单（时间/耗时/产物）
"""
import os
import sys
import json
import time
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from _common import ROOT, OUTROOT, PY, log, ensure_dir

EXPERIMENTS = {
    "A": ("run_A_scenarios.py", "六情景单箱基线/政策沙盘"),
    "B": ("run_B_surface.py", "O3–NOx–VOC 响应面（4D 参数扫描）"),
    "C": ("run_C_ridge.py", "屋脊跨路线/跨机制对比"),
    "D": ("run_D_process.py", "IPR 六账 + IRR 反应级过程诊断"),
    "E": ("run_E_indicators.py", "指示因子工具箱（12 指标）"),
    "F": ("run_F_osse.py", "观测系统模拟实验 OSSE"),
    "G": ("run_G_hourlyPA.py", "逐小时过程分析 + IPR↔IRR 互证"),
    "H": ("run_H_sensitivity.py", "OAT 敏感性（±10% 中心差分）"),
    "I": ("run_I_scenario_matrix.py", "情景 × 路线矩阵 + 方向一致性"),
    "M": ("run_M_met_series.py", "气象观测序列驱动（站点逐时）+ 单变量归因"),
}

ORDER = ["A", "B", "D", "G", "E", "F", "H", "M", "I", "C"]   # C 依赖 B，放最后


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("-")]
    if "--list" in argv:
        log("实验清单:")
        for k in sorted(EXPERIMENTS):
            log(f"  {k}  {EXPERIMENTS[k][1]:36s}  → experiments/run/{EXPERIMENTS[k][0]}")
        log("\n注: 未包含『伴随/反演』（按需求排除；另见 experiments/sens_adj.py）")
        return 0

    todo = [a.upper() for a in args] or ORDER
    bad = [t for t in todo if t not in EXPERIMENTS]
    if bad:
        log(f"✗ 未知实验: {bad}；用 --list 查看")
        return 1

    log("═" * 62)
    log("ozone_lab 实验总驱动（只跑模拟 + CSV，无绘图）")
    log(f"输出根目录: {OUTROOT}")
    log(f"计划执行: {todo}")
    log("═" * 62)

    results = {}
    t_all = time.time()
    for tag in todo:
        script, desc = EXPERIMENTS[tag]
        log(f"\n▶▶▶ 实验 {tag}: {desc}")
        t0 = time.time()
        r = subprocess.run([PY, os.path.join(HERE, script)],
                           cwd=ROOT, capture_output=True, text=True)
        dt = time.time() - t0
        ok = r.returncode == 0
        results[tag] = {"ok": ok, "elapsed_s": round(dt, 2),
                        "script": script, "desc": desc}
        # 透传（缩进）脚本输出
        for ln in (r.stdout or "").splitlines():
            log("  " + ln)
        if not ok:
            for ln in (r.stderr or "").splitlines()[-8:]:
                log("  [stderr] " + ln)
            log(f"✗ 实验 {tag} 失败（{dt:.1f}s）")
        else:
            log(f"✓ 实验 {tag} 完成（{dt:.1f}s）")

    total = time.time() - t_all
    log("\n" + "═" * 62)
    log("执行汇总")
    log("═" * 62)
    for tag in todo:
        r = results[tag]
        log(f"  {tag}  {'✓' if r['ok'] else '✗'}  {r['elapsed_s']:7.1f}s  {r['desc']}")
    log(f"\n总计 {total:.1f}s · 输出在 {OUTROOT}/run_*/")

    p = os.path.join(OUTROOT, "run_all_manifest.json")
    with open(p, "w") as f:
        json.dump({"timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                   "total_s": round(total, 2),
                   "outroot": OUTROOT,
                   "experiments": results,
                   "note": "只跑模拟+CSV；绘图见 figures/；未含伴随反演"}, f,
                  ensure_ascii=False, indent=2)
    log(f"\n清单: {os.path.relpath(p, ROOT)}")
    return 0 if all(r["ok"] for r in results.values()) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
