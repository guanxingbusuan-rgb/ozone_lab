#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
scen5.py —— 五路线光化学情景矩阵实验
  情景集 (政策语义):
    BAU    基线
    LNNOX  NOx ×0.3            (减排情景)
    HVOC   VOC ×3.0            (的增长情景)
    ISOPX  NOx ×0.5 + VOC ×2.0 (联合情景)
    HOT    ΔT = +5K            (温度情景; 仅 fac_v3 支持)
    CLOUD  JFAC ×0.6           (云情景; 仅 fac_v3 支持)
    IC×2   初值 ×2 (config 级; 仅 fac_v3 — 情景表无 IC 因子, 缺口文档化)
    BG+50% 背景 O₃ ×1.5 (config 级; 同上)
  路线 × 支持空间:
    fac_v3   74sp | 情景表全 6 因子 + config 级扩展
    mcm_run  29sp | fnox/fvoc CLI
    batch    74sp | 网格最近点 (LNNOX/HVOC/ISOPX 可读)
  输出: 路线×情景 O₃max 响应矩阵 (%Δ) + 方向一致性
用法: python3 experiments/scen5.py
"""
import os
import shutil
import subprocess
import sys
import numpy as np

# R13: Fortran 三位指数(3.07932-174)兜底 — 共享实现
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) if _os.path.basename(_os.path.dirname(_os.path.abspath(__file__))) == "experiments" else _os.path.dirname(_os.path.abspath(__file__)), "diagnostics"))
from csvfloat import fix_fortran_exp as _fixexp, fnum as _fnum

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = sys.executable
PPB = 2.46e10
WORK = f"{HERE}/output/scen5"


def sh(cmd, cwd=None):
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    if r.returncode != 0:
        print(r.stdout[-300:], r.stderr[-300:])
        raise SystemExit(f"失败: {cmd}")


def o3max(pref):
    import csv

    rows = list(csv.reader(open(f"{pref}_conc.csv")))
    hdr = [h.strip() for h in rows[0]]
    i = hdr.index("O3")
    return max(_fnum(r[i]) for r in rows[1:] if r) / PPB


def make_cfg(tag, emis_scale=None, ic_scale=None, bg_scale=None):
    src, dst = f"{HERE}/input_v3", f"{WORK}/in_{tag}"
    os.makedirs(f"{dst}/constraints", exist_ok=True)
    for f in os.listdir(src):
        p = os.path.join(src, f)
        if os.path.isfile(p):
            shutil.copy(p, dst)
    for f in os.listdir(f"{src}/constraints"):
        shutil.copy(f"{src}/constraints/{f}", f"{dst}/constraints")
    if emis_scale:
        lines = []
        for ln in open(f"{src}/emissions.config"):
            s = ln.split("#")[0].split()
            if s and s[0] in emis_scale:
                s[1] = f"{float(s[1]) * emis_scale[s[0]]:.5g}"
                lines.append("  ".join(s) + "\n")
            else:
                lines.append(ln)
        open(f"{dst}/emissions.config", "w").writelines(lines)
    if ic_scale:
        lines = []
        for ln in open(f"{src}/initialConcentrations.config"):
            s = ln.split("#")[0].split()
            if s and s[0] in ic_scale:
                s[1] = f"{float(s[1]) * ic_scale[s[0]]:.5g}"
                lines.append("  ".join(s) + "\n")
            else:
                lines.append(ln)
        open(f"{dst}/initialConcentrations.config", "w").writelines(lines)
    if bg_scale:
        lines = []
        for ln in open(f"{src}/bc.config"):
            s = ln.split("#")[0].split("=")
            if len(s) == 2 and s[0].strip() == "O3":
                lines.append(f"O3 = {float(s[1]) * bg_scale:.4g}\n")
            else:
                lines.append(ln)
        open(f"{dst}/bc.config", "w").writelines(lines)
    rt = f"{WORK}/rt_{tag}.dat"
    sh([PY, f"{HERE}/src/cfgbuild.py", dst, rt,
        "--mech", f"{HERE}/mechanism/mechanism_3_full_chain.fac"])
    return rt


def main():
    os.makedirs(WORK, exist_ok=True)
    M = {}   # M[route][scen] = O3max ppb 或 None (不支持)

    # ---- fac_v3: 情景表 6 情景 + config 级 2 扩展 ----
    print("运行 fac_v3 (74sp, 情景表 + config 扩展)...")
    M["fac_v3 74sp"] = {}
    for scen in ["SCEN_BAU", "SCEN_LNNOX", "SCEN_HVOC", "SCEN_ISOPX",
                 "SCEN_HOT", "SCEN_CLOUD"]:
        sh([f"{HERE}/facsimile_impl/fac_v3", scen, f"{WORK}/v3_{scen}",
            f"{HERE}/runtime_v3.dat"], cwd=f"{HERE}/facsimile_impl")
        M["fac_v3 74sp"][scen.replace("SCEN_", "")] = o3max(f"{WORK}/v3_{scen}")
    rt = make_cfg("v3ic", ic_scale={s: 2.0 for s in
                  ["NO", "NO2", "O3", "CO", "HCHO", "CH4", "ISOP", "SO2"]})
    M["fac_v3 74sp"]["IC×2"] = o3max(run_facv3_cfg(rt, f"{WORK}/v3_ic"))
    rt = make_cfg("v3bg", bg_scale=1.5)
    M["fac_v3 74sp"]["BG+50%"] = o3max(run_facv3_cfg(rt, f"{WORK}/v3_bg"))

    # ---- fac_run/kpp_run 已退役 (15sp 层) ----

    # ---- mcm_run (29sp, fnox/fvoc CLI) ----
    print("运行 mcm_run (29sp)...")
    M["mcm_run 29sp"] = {}
    for tag, fn, fv in [("BAU", 1.0, 1.0), ("LNNOX", 0.3, 1.0),
                        ("HVOC", 1.0, 3.0), ("ISOPX", 0.5, 2.0)]:
        sh([f"{HERE}/mcm_city/mcm_run", "5", f"{WORK}/mc_{tag}", f"{fn}", f"{fv}"],
           cwd=f"{HERE}/mcm_city")
        for suf in ("conc", "irr", "ipr"):
            src = f"{HERE}/mcm_city/mcm_city_{suf}.csv"
            if os.path.exists(src):
                os.replace(src, f"{WORK}/mc_{tag}_{suf}.csv")
        M["mcm_run 29sp"][tag] = o3max(f"{WORK}/mc_{tag}")
        M["mcm_run 29sp"]["HOT"] = None
        M["mcm_run 29sp"]["CLOUD"] = None
        M["mcm_run 29sp"]["IC×2"] = None
        M["mcm_run 29sp"]["BG+50%"] = None

    # ---- batch (74sp 网格最近点) ----
    print("读取 batch 网格最近点...")
    lines = [l for l in open(f"{HERE}/output/batch_v3_indicators.csv")
             if l.strip() and not l.startswith("fnox")]
    grid = [[_fnum(x) for x in l.split(",")] for l in lines]
    M["batch 74sp"] = {"BAU": near(grid, 1.0, 1.0), "LNNOX": near(grid, 0.3, 1.0),
                       "HVOC": near(grid, 1.0, 3.0), "ISOPX": near(grid, 0.5, 2.0),
                       "HOT": None, "CLOUD": None, "IC×2": None, "BG+50%": None}

    # ---- 响应矩阵 ----
    scens = ["LNNOX", "HVOC", "ISOPX", "HOT", "CLOUD", "IC×2", "BG+50%"]
    print("\n===== 五路线情景矩阵: O₃max [ppb] =====")
    print(f"{'路线':18s} " + " ".join(f"{s:>9s}" for s in ["BAU"] + scens))
    for route, row in M.items():
        cells = []
        for s in ["BAU"] + scens:
            v = row.get(s)
            cells.append(f"{v:9.2f}" if v else "      —  ")
        print(f"{route:18s} " + " ".join(cells))

    print("\n===== 响应方向一致性 (%Δ vs 各路线 BAU; ✓=与 fac_v3 同向) =====")
    ref = M["fac_v3 74sp"]
    hdr = f"{'路线':18s} " + " ".join(f"{s:>10s}" for s in scens)
    print(hdr)
    for route, row in M.items():
        base = row.get("BAU")
        cells = []
        for s in scens:
            v = row.get(s)
            if v is None or not base:
                cells.append(f"{'—':>10s}")
                continue
            pct = (v / base - 1) * 100
            refpct = (ref[s] / ref["BAU"] - 1) * 100 if ref.get(s) else 0
            agree = "✓" if (pct * refpct > 0) or abs(refpct) < 1 else "✗"
            cells.append(f"{pct:+9.1f}{agree}")
        print(f"{route:18s} " + " ".join(cells))
    print(f"\n数据目录: {WORK}/")


def run_facv3_cfg(rt, pref):
    sh([f"{HERE}/facsimile_impl/fac_v3", "SCEN_BAU", pref, rt],
       cwd=f"{HERE}/facsimile_impl")
    return pref


def near(grid, fn, fv):
    row = min(grid, key=lambda r: (r[0] - fn) ** 2 + (r[1] - fv) ** 2)
    return row[4]


if __name__ == "__main__":
    main()
