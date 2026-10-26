#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
obs_osse.py —— 观测约束 OSSE (Observing System Simulation Experiment)
  双胞胎设计: 真值运行 (NOx/ISOP 排放 ×1.4) 采样加噪 → 观测文件 → 基础配置约束
  运行 → 技能指标。三条同化通路 (物理语义各异):
    ① 地面 O₃ 站点 → 张弛 nudging (τ 扫描 2/6/24/72h)
    ② VOCs 组份站点 (HCHO/ISOP) → 同①
    ③ 卫星 NO₂ 柱 (转化后) → 【不做直接 nudge】— NO₂ 受光稳态支配, 强直接
       约束破坏 PSS (本轮实验实测 h→1e-79 崩溃); 正确通路 = 自上而下排放反演
       (E_NOx ×= mean(obs)/mean(free), 经典排放反演一步法)
  分析: RMSE 改善 / τ 强度曲线 / IPR B 账同化倾向。
  路线能力: fac_v3 / fac_batch_omp (chem_phys_mod bc_s) 支持观测约束;
  fac_run/mcm_run 仅常数背景 (机制验证定位); kpp 纯机制无约束 (定位使然)。
用法: python3 experiments/obs_osse.py
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
WORK = f"{HERE}/output/osse"
OBS_TIMES = [96, 101, 106, 111, 116, 120]       # 技能评估窗内 (小时)
NOISE = {"O3": 0.04, "HCHO": 0.10, "ISOP": 0.15, "NO2": 0.08}
STATION = ["O3", "HCHO", "ISOP"]


def sh(cmd, cwd=None):
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
    if r.returncode != 0:
        print(r.stdout[-400:], r.stderr[-400:])
        raise SystemExit(f"命令失败: {cmd}")


def make_input(tag, emissions_scale, tau, obs=None):
    """生成输入目录: base 配置 + 排放缩放 (truth/反演) + 约束文件 (constrained)"""
    src, dst = f"{HERE}/input_v3", f"{WORK}/in_{tag}"
    os.makedirs(f"{dst}/constraints", exist_ok=True)
    for f in os.listdir(src):
        p = os.path.join(src, f)
        if os.path.isfile(p):
            shutil.copy(p, dst)
    for f in os.listdir(f"{src}/constraints"):
        shutil.copy(f"{src}/constraints/{f}", f"{dst}/constraints")
    if emissions_scale:
        lines = []
        for ln in open(f"{src}/emissions.config"):
            s = ln.split("#")[0].split()
            if s and s[0] in emissions_scale:
                s[1] = f"{float(s[1]) * emissions_scale[s[0]]:.4g}"
                lines.append("  ".join(s) + "\n")
            else:
                lines.append(ln)
        open(f"{dst}/emissions.config", "w").writelines(lines)
    if obs:
        for sp, (times, vals) in obs.items():
            with open(f"{dst}/constraints/obs_{sp}.txt", "w") as f:
                f.write(f"# unit: ppb\n# station: OSSE-synth  species: {sp}\n")
                f.write("# t[s]  value\n")
                for t, v in zip(times, vals):
                    f.write(f"{int(t * 3600):8d}  {v:8.3f}\n")
        with open(f"{dst}/bc.config", "w") as f:
            f.write(f"# OSSE 站点约束 (τ={tau}s)\nTAU_MIX = {tau}\n")
            for sp in obs:
                f.write(f"{sp} = file:constraints/obs_{sp}.txt\n")
    rt = f"{WORK}/rt_{tag}.dat"
    sh([PY, f"{HERE}/src/cfgbuild.py", dst, rt,
        "--mech", f"{HERE}/mechanism/mechanism_3_full_chain.fac"])
    return rt


def run_facv3(rt, pref):
    sh([FACV3, "SCEN_BAU", pref, rt], cwd=f"{HERE}/facsimile_impl")
    return f"{pref}_conc.csv"


def read_conc(path, sp):
    import csv
    rows = list(csv.reader(open(path)))
    hdr = [h.strip() for h in rows[0]]
    i = hdr.index(sp)
    h = [_fnum(r[0]) for r in rows[1:] if r]
    v = [_fnum(r[i]) / PPB for r in rows[1:] if r]
    return np.array(h), np.array(v)


def main():
    os.makedirs(WORK, exist_ok=True)
    rng = np.random.default_rng(42)

    # ---- ① 真值运行 (NOx/ISOP 排放 ×1.4) ----
    rt_true = make_input("truth", {"NO": 1.4, "NO2": 1.4, "ISOP": 1.4}, 86400)
    c_true = run_facv3(rt_true, f"{WORK}/truth")

    # ---- ② 真值采样 + 噪声 → 观测 ----
    obs = {}
    for sp, form in [("O3", "地面 O₃ 站点"), ("HCHO", "VOCs 组份站点"),
                     ("ISOP", "VOCs 组份站点"), ("NO2", "卫星转化 (TROPOMI 13:30)")]:
        h, v = read_conc(c_true, sp)
        times, vals = [], []
        for t in OBS_TIMES:
            k = int(np.argmin(np.abs(h - t)))
            noise = 1.0 + rng.normal(0, NOISE[sp])
            times.append(h[k])
            vals.append(max(v[k] * noise, 1e-4))
        obs[sp] = (np.array(times), np.array(vals))
        print(f"  观测 {sp} ({form}): " +
              " ".join(f"h{int(t)}={v:.2f}" for t, v in zip(times, vals)))

    def skill(path, sp, times, vals):
        h, v = read_conc(path, sp)
        e = [v[int(np.argmin(np.abs(h - t)))] - o for t, o in zip(times, vals)]
        return float(np.sqrt(np.mean(np.square(e))))

    # ---- ③ 站点物种: 自由 vs τ 扫描 nudging ----
    rt_free = make_input("free", None, 86400)
    c_free = run_facv3(rt_free, f"{WORK}/free")
    print("\n===== 通路①② 技能: RMSE [ppb] (站点物种 nudge) =====")
    print(f"{'τ':>6s} " + " ".join(f"{sp:>8s}" for sp in STATION))
    free_rmse = {sp: skill(c_free, sp, *obs[sp]) for sp in STATION}
    print(f"{'free':>6s} " + " ".join(f"{free_rmse[sp]:8.2f}" for sp in STATION))
    rmse = {}
    for tau_h in (2, 6, 24, 72):
        rt = make_input(f"c{tau_h}", None, tau_h * 3600.0,
                        obs={sp: obs[sp] for sp in STATION})
        try:
            c = run_facv3(rt, f"{WORK}/c{tau_h}")
        except SystemExit:
            print(f"{tau_h:5d}h  运行失败 (强 nudge 破坏光稳态)")
            continue
        rmse[tau_h] = {sp: skill(c, sp, *obs[sp]) for sp in STATION}
        print(f"{tau_h:5d}h " + " ".join(f"{rmse[tau_h][sp]:8.2f}" for sp in STATION))

    print("\n----- 过冲指标 (O₃max, 强 nudge 在锚点间振荡的隐藏代价) -----")
    h_o3, _ = read_conc(c_free, "O3")
    print(f"  free: O₃max = {read_conc(c_free, 'O3')[1].max():.1f} ppb")
    for tau_h in (2, 6, 24, 72):
        cpath = f"{WORK}/c{tau_h}_conc.csv"
        if os.path.exists(cpath):
            print(f"  τ={tau_h:3d}h: O₃max = {read_conc(cpath, 'O3')[1].max():.1f} ppb")

    print("\n----- RMSE 改善 (自由 → 最优 τ) -----")
    for sp in STATION:
        cands = {t: r[sp] for t, r in rmse.items() if sp in r}
        if not cands:
            continue
        bt = min(cands, key=cands.get)
        print(f"  {sp:5s}: {free_rmse[sp]:.2f} → {cands[bt]:.2f} ppb "
              f"(τ={bt}h, 改善 {(1 - cands[bt]/free_rmse[sp])*100:+.0f}%)")

    # ---- ④ 卫星 NO2: 自上而下排放反演 (不直接 nudge) ----
    print("\n===== 通路③ 卫星 NO2: 自上而下排放反演 =====")
    t_n, v_n = obs["NO2"]
    h_f, v_f = read_conc(c_free, "NO2")
    mo = float(np.mean(v_n))
    mf = float(np.mean([v_f[int(np.argmin(np.abs(h_f - t)))] for t in t_n]))
    scale = mo / mf
    print(f"  mean(obs)={mo:.2f} / mean(free)={mf:.2f} → E_NOx × {scale:.2f}")
    rt_inv = make_input("inv", {"NO": scale, "NO2": scale}, 86400)
    c_inv = run_facv3(rt_inv, f"{WORK}/inv")
    inv_rmse = skill(c_inv, "NO2", t_n, v_n)
    free_n = skill(c_free, "NO2", t_n, v_n)
    print(f"  NO2 RMSE: free {free_n:.2f} → 反演 {inv_rmse:.2f} ppb "
          f"(改善 {(1 - inv_rmse/free_n)*100:+.0f}%)")
    # 反演对 O3 的连带效应 (NOx 反演改变化学)
    o3_free = skill(c_free, "O3", *obs["O3"])
    o3_inv = skill(c_inv, "O3", *obs["O3"])
    print(f"  O3 RMSE 连带: {o3_free:.2f} → {o3_inv:.2f} ppb")

    # ---- ⑤ IPR B 账同化倾向 (τ=6h) ----
    print("\n===== 同化倾向 (IPR B 账, τ=6h 约束运行, 末 24h 窗) =====")
    ipr_path = f"{WORK}/c6_ipr.csv"
    if os.path.exists(ipr_path):
        import csv
        import re
        rws = list(csv.reader(open(ipr_path)))
        hdr = [h.strip() for h in rws[0]]

        def fexp(tok):
            if "E" not in tok and "e" not in tok:
                m = re.search(r"(\d)-(\d{3})$", tok)
                if m:
                    return tok[:m.start(1) + 1] + "E" + tok[m.start(1) + 1:]
            return tok
        D = np.array([[float(fexp(x)) for x in r] for r in rws[1:] if r])
        for sp in STATION:
            if f"{sp}_B" not in hdr:
                continue
            ib = hdr.index(f"{sp}_B")
            m = D[1:, ib].mean()          # 排除首窗 (spin-up)
            rate = m / 3600.0                  # 窗口积分 → 平均速率
            tot = np.abs(D[1:, [hdr.index(f"{sp}_{p}") for p in "PLEDBX"]]).mean(axis=0).sum() / 3600
            print(f"  {sp:5s} B 账均速率={rate:+.3e} molec/cm³/s "
                  f"({'拉向观测' if abs(rate) > 1e3 else '≈0 (已贴合)'}), "
                  f"|B| 占过程总账 {abs(rate)/max(tot,1e-30)*100:.0f}%")
    else:
        print("  (τ=6h 运行不存在 — 跳过)")
    print(f"\n数据目录: {WORK}/")


if __name__ == "__main__":
    main()
