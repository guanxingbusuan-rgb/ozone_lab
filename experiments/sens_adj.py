#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
sens_adj.py —— 排放清单+观测数据的灵敏性/伴随分析
  ① 排放 OAT 雅可比: 9 个排放物种各 ±10% (config 级扰动, fac_v3) →
     观测量 (O3/HCHO/ISOP/NO2 @OSSE 时刻) 的灵敏度矩阵 J
  ② 线性化反演验证: min||J·ΔE − (obs−free)|| 与 OSSE 真值 ΔE 对照
     (truth: NO/NO2/ISOP 排放 ×1.4 → ΔE/E = +0.4)
  ③ 跨路线一致性: fac_v3 OAT 组斜率 vs batch 74sp 网格局部斜率
      (替代 modules/rate_sensitivity.py 的硬编码假实现)
用法: python3 experiments/sens_adj.py
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
FACV3 = f"{HERE}/facsimile_impl/fac_v3"
PPB = 2.46e10
WORK = f"{HERE}/output/sens"
OBS_TIMES = [96, 101, 106, 111, 116, 120]
OBS_SP = ["O3", "HCHO", "ISOP", "NO2"]
NOISE = {"O3": 0.04, "HCHO": 0.10, "ISOP": 0.15, "NO2": 0.08}
EMIT = ["NO", "NO2", "CO", "HCHO", "CH3CHO", "CH4", "ISOP", "SO2", "DMS"]


def sh(cmd, cwd=None):
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    if r.returncode != 0:
        print(r.stdout[-300:], r.stderr[-300:])
        raise SystemExit(f"失败: {cmd}")


def make_input(tag, scale=None, tau=86400.0, obs=None):
    src, dst = f"{HERE}/input_v3", f"{WORK}/in_{tag}"
    os.makedirs(f"{dst}/constraints", exist_ok=True)
    for f in os.listdir(src):
        p = os.path.join(src, f)
        if os.path.isfile(p):
            shutil.copy(p, dst)
    for f in os.listdir(f"{src}/constraints"):
        shutil.copy(f"{src}/constraints/{f}", f"{dst}/constraints")
    if scale:
        lines = []
        for ln in open(f"{src}/emissions.config"):
            s = ln.split("#")[0].split()
            if s and s[0] in scale:
                s[1] = f"{float(s[1]) * scale[s[0]]:.5g}"
                lines.append("  ".join(s) + "\n")
            else:
                lines.append(ln)
        open(f"{dst}/emissions.config", "w").writelines(lines)
    if obs:
        for sp, (times, vals) in obs.items():
            with open(f"{dst}/constraints/obs_{sp}.txt", "w") as f:
                f.write(f"# unit: ppb\n# t[s]  {sp}\n")
                for t, v in zip(times, vals):
                    f.write(f"{int(t * 3600):8d}  {v:8.3f}\n")
        with open(f"{dst}/bc.config", "w") as f:
            f.write(f"TAU_MIX = {tau}\n")
            for sp in obs:
                f.write(f"{sp} = file:constraints/obs_{sp}.txt\n")
    rt = f"{WORK}/rt_{tag}.dat"
    sh([PY, f"{HERE}/src/cfgbuild.py", dst, rt,
        "--mech", f"{HERE}/mechanism/mechanism_3_full_chain.fac"])
    return rt


def run_v3(rt, pref):
    sh([FACV3, "SCEN_BAU", pref, rt], cwd=f"{HERE}/facsimile_impl")
    return f"{pref}_conc.csv"


def read_col(path, sp):
    import csv
    rows = list(csv.reader(open(path)))
    hdr = [h.strip() for h in rows[0]]
    i = hdr.index(sp)
    h = np.array([_fnum(r[0]) for r in rows[1:] if r])
    v = np.array([_fnum(r[i]) / PPB for r in rows[1:] if r])
    return h, v


def obs_vec(path):
    """观测量向量 (4 物种 × 6 时刻, ppb)"""
    y = []
    for sp in OBS_SP:
        h, v = read_col(path, sp)
        for t in OBS_TIMES:
            y.append(v[int(np.argmin(np.abs(h - t)))])
    return np.array(y)


def main():
    os.makedirs(WORK, exist_ok=True)
    rng = np.random.default_rng(42)

    # ---- ① OSSE 真值与观测 (与 obs_osse.py 同种子同噪声, 可复现) ----
    rt_true = make_input("truth", {"NO": 1.4, "NO2": 1.4, "ISOP": 1.4})
    c_true = run_v3(rt_true, f"{WORK}/truth")
    y_true, obs = [], {}
    for sp in OBS_SP:
        h, v = read_col(c_true, sp)
        for t in OBS_TIMES:
            k = int(np.argmin(np.abs(h - t)))
            y_true.append(v[k])
            obs.setdefault(sp, ([], []))
            obs[sp][0].append(h[k])
            obs[sp][1].append(v[k] * (1 + rng.normal(0, NOISE[sp])))
    y_obs = np.array([o for sp in OBS_SP for o in obs[sp][1]])

    # ---- ② 基线与排放 OAT 雅可比 ----
    rt_base = make_input("base", None)
    c_base = run_v3(rt_base, f"{WORK}/base")
    y_base = obs_vec(c_base)
    print("===== ① 排放 OAT 雅可比 (±10%, 观测向量 24 维) =====")
    J = np.zeros((len(y_base), len(EMIT)))
    for j, sp in enumerate(EMIT):
        c_p = run_v3(make_input(f"oat{j}p", {sp: 1.10}), f"{WORK}/oat{j}p")
        c_m = run_v3(make_input(f"oat{j}m", {sp: 0.90}), f"{WORK}/oat{j}m")
        J[:, j] = (obs_vec(c_p) - obs_vec(c_m)) / (0.2 * 1.0)   # ∂y/∂(E/E0), 中心差分
        print(f"  {sp:6s}: ||J_col|| = {np.linalg.norm(J[:, j]):.3f}")
    np.savetxt(f"{WORK}/jacobian.csv", J, delimiter=",",
               header=",".join(EMIT), comments="")

    # ---- ③ 线性化反演 vs OSSE 真值 ----
    print("\n===== ② 线性化反演 (截断 SVD 正则化) vs 真值 ΔE/E =====")
    r = y_obs - y_base
    U, S, Vt = np.linalg.svd(J, full_matrices=False)
    cond = S[0] / S[-1]
    keep = S > S[0] * 0.1                       # 截断奇异值 (病态防护)
    Sinv = np.where(keep, 1 / np.where(S > 0, S, 1), 0)
    x = Vt.T @ (Sinv * (U.T @ r))
    print(f"  奇异值条件数 = {cond:.0f}, 保留 {keep.sum()}/{len(S)} 个模态"
          f" (截断: 排放可辨识性有限 — 非线性化学致响应近似共线)")
    x = np.where(np.abs(x) > 50, 0, x)          # 残余病态模态置零
    truth = {"NO": 0.4, "NO2": 0.4, "ISOP": 0.4}
    print(f"{'排放':8s} {'反演 ΔE/E':>10s} {'真值':>6s}")
    for j, sp in enumerate(EMIT):
        tv = truth.get(sp, 0.0)
        mark = "←" if tv else ""
        print(f"{sp:8s} {x[j]:10.2f} {tv:6.1f} {mark}")
    pred = y_base + J @ x
    print(f"  反演后拟合残差: {np.linalg.norm(pred - y_obs):.2f} (基线 {np.linalg.norm(r):.2f})")

    # ---- ④ 跨路线一致性: fac_v3 OAT 组斜率 vs batch 网格局部斜率 ----
    print("\n===== ③ 跨路线一致性: O₃max 对 ln(fnox) 的局部斜率 =====")
    c_p = run_v3(make_input("gx1p", {"NO": 1.1, "NO2": 1.1}), f"{WORK}/gx1p")
    c_m = run_v3(make_input("gx1m", {"NO": 0.9, "NO2": 0.9}), f"{WORK}/gx1m")
    o3 = lambda p: read_col(p, "O3")[1].max()
    slope_oat = (o3(c_p) - o3(c_m)) / (0.2 * np.log(1.1 / 0.9)) * np.log(1.1 / 0.9) / 0.2
    slope_oat = (o3(c_p) - o3(c_m)) / 0.2 / 1.0125          # per unit fnox
    lines = [l for l in open(f"{HERE}/output/batch_v3_indicators.csv")
             if l.strip() and not l.startswith("fnox")]
    g = [([_fnum(x) for x in l.split(",")[:4]]) + [_fnum(l.split(",")[4])]
         for l in lines]
    near = lambda arr, fn, fv: min(
        g, key=lambda row: (row[0]-fn)**2 + (row[1]-fv)**2)
    hi = near(g, 1.0125, 1.0125); lo = near(g, 0.675, 1.0125)
    slope_grid = (hi[4] - lo[4]) / (hi[0] - lo[0])
    print(f"  fac_v3 OAT (NO+NO2 ±10%): {slope_oat:.1f} ppb / unit fnox")
    print(f"  batch 网格局部斜率:        {slope_grid:.1f} ppb / unit fnox")
    print(f"  一致性: {abs(slope_oat/slope_grid - 1) * 100:.0f}% 偏差")

    print("\n④ 真速率灵敏度: 随 fac_run/15sp 层退役 (原载体=fac_run 重编译); "
          "工程机制速率谱待重建 (载体=fac_v3 重编译, 结构同 ④)。\n")
if __name__ == "__main__":
    main()
