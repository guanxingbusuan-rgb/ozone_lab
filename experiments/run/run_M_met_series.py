#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_M_met_series.py —— 实验 M：气象观测序列接入（站点逐时气象驱动）

回答的问题:
  "ozone_lab 能不能像监测站点约束那样，把气象变量也按时间序列喂进去？"
  能 —— 6 个气象量 (TEMP/PRES/RH/BLH/RAIN/JFAC) 都可以是 file: 序列,
  与 BC 浓度序列走同一套 SERIES 机制 (cfgbuild 编译进 runtime.dat)。

本实验用"同一平均强迫、只改时间结构"的对照设计, 把气象时序的影响分离出来:
  M0_const   六变量全取各自序列的 120 h 均值 (常数)      ← 对照基准
  M1_all     六变量全部走逐时序列                        ← 真实驱动
  M2_*..M7_* 每次只让一个变量走序列, 其余保持均值          ← 归因

输出 (全部 CSV, 无绘图):
  output/run_M/<run>_conc.csv       末 24 h 逐时 74 物种 (fac_v3 原始输出)
  output/run_M/summary.csv          每 run 的 O3 峰值/均值/日振幅/FNR/耗时
  output/run_M/o3_timeseries.csv    h + 各 run 的 O3 逐时值 (直接可画)
  output/run_M/attribution.csv      各变量"时序化"相对常数的 O3 变化量

用法:
  python experiments/run/run_M_met_series.py
  python experiments/run/run_M_met_series.py --const-only   # 只跑基准
  OZONE_OUT=/data/out python experiments/run/run_M_met_series.py
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (ROOT, INPUT, PY, CFGBUILD, MECH, require_prereq, out_dir,
                     run_fac_v3, clone_input, read_table, col, write_csv,
                     manifest, summary_table, log, run)

# 变量 → (参数名, 序列文件名, 原生单位, 该量的物理身份)
VARS = [
    ("TEMP", "constraints/met_TEMP_demo.txt", "K",    "气温 → Arrhenius + 数密度 M"),
    ("PRES", "constraints/met_PRES_demo.txt", "hPa",  "气压 → 数密度 M (三体反应)"),
    ("RH",   "constraints/met_RH_demo.txt",   "%",    "湿度 → H2O → O(1D)+H2O 供 OH"),
    ("BLH",  "constraints/met_BLH_demo.txt",  "m",    "边界层 → 体积/稀释/沉降/膨胀"),
    ("RAIN", "constraints/met_RAIN_demo.txt", "mm/h", "降水 → 湿沉降 k_wet"),
    ("JFAC", "constraints/met_JFAC_demo.txt", "1",    "光解衰减 → 全部 J 的乘子"),
]


def read_series_raw(path):
    """读两列 't v'（跳过注释），返回 (t, v) 原生单位 numpy 数组"""
    t, v = [], []
    for ln in open(path):
        s = ln.split("#")[0].strip()
        if not s:
            continue
        p = s.split()
        t.append(float(p[0]))
        v.append(float(p[1]))
    return np.array(t), np.array(v)


def write_meteo(dst, series_vars, means):
    """写 meteo.config: series_vars 中的走 file:, 其余写常数均值"""
    lines = ["# 实验 M 自动生成 —— 气象驱动配置",
             "# 序列型变量走 file:, 其余固定为对应序列的 120 h 均值"]
    for k, rel, unit, _ in VARS:
        if k in series_vars:
            lines.append(f"{k:<7s}= file:{rel}")
        else:
            lines.append(f"{k:<7s}= {means[k]:.6g}")
    lines.append("RADIAT = 0          # 当前无消费者 (cfgbuild 会告警), 留位给在线光解")
    open(os.path.join(dst, "meteo.config"), "w").write("\n".join(lines) + "\n")


def analyze(pref):
    """从 conc CSV（末 24 h 逐时）取 O3 峰值/均值/日振幅 + FNR 代理"""
    hdr, data = read_table(f"{pref}_conc.csv")
    if not data:
        return {}
    o3 = col(hdr, data, "O3") / 2.46e10
    hcho = col(hdr, data, "HCHO") / 2.46e10
    no2 = col(hdr, data, "NO2") / 2.46e10
    fnr = hcho / np.maximum(no2, 1e-12)
    return {"o3": o3, "h": col(hdr, data, "h"),
            "o3_max": float(np.nanmax(o3)), "o3_mean": float(np.nanmean(o3)),
            "o3_amp": float(np.nanmax(o3) - np.nanmin(o3)),
            "fnr": float(np.nanmean(fnr))}


def main(argv):
    const_only = "--const-only" in argv
    log("═══ 实验 M：气象观测序列驱动（站点逐时）═══")
    require_prereq()
    d = out_dir("run_M")

    # ── ① 读 6 条序列的均值 (作为"常数基准"的等价强迫) ──
    means, raw = {}, {}
    for k, rel, unit, desc in VARS:
        p = os.path.join(INPUT, rel)
        if not os.path.exists(p):
            log(f"✗ 缺少样例序列 {rel} —— 请先确认 input_v3/constraints/ 完整")
            return 1
        t, v = read_series_raw(p)
        raw[k] = (t, v)
        means[k] = float(np.mean(v))
        log(f"  {k:<5s} {len(t):3d} 点  [{v.min():.4g}, {v.max():.4g}] {unit:<5s} "
            f"均值 {means[k]:.4g}   {desc}")
    if len({len(v[0]) for v in raw.values()}) != 1:
        log("⚠ 六条序列的点数不一致 —— 插值仍可用, 但对照不如等长严格")

    # ── ② 构造 run 清单 ──
    runs = [("M0_const", set())]
    if not const_only:
        runs.append(("M1_all", {k for k, *_ in VARS}))
        for k, *_ in VARS:
            runs.append((f"M{2 + [v[0] for v in VARS].index(k)}_{k}_only", {k}))

    # ── ③ 逐个编译 + 运行 ──
    rows, products, o3cols = [], [], {}
    h_axis = None
    for name, series_vars in runs:
        dst = os.path.join(d, f"in_{name}")
        rt = clone_input(dst)              # 复制 input_v3 并用原 meteo 编译一次
        write_meteo(dst, series_vars, means)   # 改成本 run 的气象驱动
        run([PY, CFGBUILD, dst, rt, "--mech", MECH])   # 内容变了 → 必须重编

        pref = os.path.join(d, name)
        r = run_fac_v3("SCEN_BAU", pref, rt=rt)
        a = analyze(pref)
        if not a:
            log(f"  ✗ {name} 无输出")
            return 1
        o3cols[name] = a["o3"]
        if h_axis is None:
            h_axis = a["h"]
        rows.append([name, "+".join(sorted(series_vars)) or "(常数)",
                     round(a["o3_max"], 2), round(a["o3_mean"], 2),
                     round(a["o3_amp"], 2), round(a["fnr"], 3),
                     r["nst"], round(r["elapsed_s"], 2)])
        products += [f"{pref}_conc.csv", f"{pref}_irr.csv", f"{pref}_ipr.csv"]
        log(f"    {name:<18s} O3max={a['o3_max']:7.2f}  均值={a['o3_mean']:6.2f}  "
            f"日振幅={a['o3_amp']:5.2f}  FNR={a['fnr']:5.2f}  ({r['elapsed_s']:.2f}s)")

    # ── ④ summary.csv ──
    p = write_csv(os.path.join(d, "summary.csv"),
                  ["run", "met_series_vars", "O3max_ppb", "O3_mean_last24h_ppb",
                   "O3_diurnal_amp_ppb", "FNR_mean", "nst", "elapsed_s"], rows,
                  note="实验M 气象序列驱动 | 对照设计: 同一 120h 均值, 只改时间结构\n"
                       "O3max 来自 fac_v3 stdout; 均值/振幅/FNR 来自末 24h conc CSV")
    products.append(p)

    # ── ⑤ o3_timeseries.csv (逐时可画) ──
    if o3cols:
        names = [n for n, _ in runs]
        tsrows = [[int(h_axis[i])] + [round(float(o3cols[n][i]), 3) for n in names]
                  for i in range(len(h_axis))]
        p = write_csv(os.path.join(d, "o3_timeseries.csv"),
                      ["h"] + names, tsrows,
                      note="实验M 末 24h O3 逐时值 (ppb) | 首列 h 为模拟小时")
        products.append(p)

    # ── ⑥ attribution.csv (各变量时序化的净效应) ──
    if not const_only:
        base = [r for r in rows if r[0] == "M0_const"][0]
        allr = [r for r in rows if r[0] == "M1_all"][0]
        attr = []
        for r in rows:
            if r[0] in ("M0_const", "M1_all"):
                continue
            attr.append([r[1].replace("_only", ""),
                         round(r[2] - base[2], 2), round(r[3] - base[3], 2),
                         round(r[4] - base[4], 2),
                         round((r[2] - base[2]) / max(abs(allr[2] - base[2]), 1e-9) * 100, 1)])
        attr.sort(key=lambda x: -abs(x[1]))
        p = write_csv(os.path.join(d, "attribution.csv"),
                      ["variable", "dO3max_ppb", "dO3mean_ppb", "dO3amp_ppb",
                       "share_of_full_effect_%"], attr,
                      note="实验M 归因 | 单变量时序化 − 常数基准; 末列 = 该变量占"
                           "『六变量全时序化』总效应的百分比 (可正可负, 求和不等于100%)")
        products.append(p)
        summary_table("实验 M · 单变量时序化归因（相对常数基准）",
                      ["variable", "dO3max", "dO3mean", "dO3amp", "share%"], attr)

    # ── ⑦ 收尾 ──
    summary_table("实验 M · 全部 run",
                  ["run", "series", "O3max", "O3mean", "O3amp", "FNR", "nst", "t_s"], rows)
    mp = manifest(d, "M_met_series",
                  {"const_baseline_means": {k: round(v, 6) for k, v in means.items()},
                   "runs": [n for n, _ in runs],
                   "design": "同一120h均值、只改时间结构的对照"},
                  products)
    log(f"\n产物: {os.path.relpath(d, ROOT)}/  ({len(products)} 个 CSV + manifest.json)")
    log("提示: o3_timeseries.csv 可直接喂 matplotlib; 站点格式见 input_v3/constraints/met_*_demo.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
