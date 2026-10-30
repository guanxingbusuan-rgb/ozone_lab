#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
EKMA-emission/run.py —— S1 排放正交扰动响应面（EKMA-emission）

实验定义（用户口径）:
  · 固定: 机制 .fac（V2.5+ 冻结） + 本目录 config_override/ 中的气象/初值
  · 变量: emissions.config 的 NOx 排放组 × NMVOC 排放组，正交扰动
  · 扰动范围: 默认 0%（无扰动 = baseline）→ 625%（+625%），**25% 步长**
      → 因子 1.00, 1.25, …, 7.25 共 26 点/轴
  · 网格: 默认 26 × 26 = 676 个独立 0-D 箱案例（fac_batch_omp 并行）
  · 时长: **模拟 4 天 / SPINUP 3 天**（响应统计 = 末 24 步 = 第 4 天；
    batch 的 SPINUP 仅作 ≥24 步闸门，聚合窗固定为最后一天）

对应物理量:
  x 轴 = E_VOC / E_VOC,BAU   （NMVOC 排放组扰动 %）
  y 轴 = E_NOx / E_NOx,BAU   （NOx 排放组扰动 %）
  色标 = O3max [ppb]（末 24 步窗口最大值）

输出（本目录 output/，全部由本文件写出）:
  runtime_emission.dat            本轮有效输入记录（cfgbuild 产物）
  ekma_emission_grid_raw.csv      fac_batch_omp 原始 24 列（默认 26×26 = 676 行）
  ekma_emission_grid_raw.csv_{ipr,irr,state}.csv   执行器侧车（**平均倾向** [molec/cm3/s]，
                                 与 fac_v3 的窗口积分口径差 dtout 倍）
  ekma_emission_grid.csv          规整表（pert_nox_pct, pert_voc_pct, 指标列, VNr_group/VNr_nmv）
  ekma_emission_ridge.csv         每 VOC 列的屋脊点
  ekma_emission_derived.csv       PR_bau / PR_vs_baseNOx / VNr_ridge_* / ridge_at_edge
  run_meta.txt                    运行元数据（网格/时长/耗时/失败/降级/氮守恒/侧车单位）

用法:
  python EKMA-emission/run.py                # 默认 0–625% / 25%（26×26 = 676 案例）
  python EKMA-emission/run.py --pmax 200 --step 5
  python EKMA-emission/run.py --from path/to/ekma_emission_grid_raw.csv   # 只重导出
"""
import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "output")
BATCH = os.path.join(ROOT, "facsimile_impl", "fac_batch_omp")
MECH = os.path.join(ROOT, "mechanism", "mechanism_3_full_chain.fac")
CFGBUILD = os.path.join(ROOT, "inputproc", "cfgbuild.py")
OVERRIDE_DIR = os.path.join(HERE, "config_override")   # 本实验唯一输入目录（含 constraints/）
REQUIRED_INPUT = ("model.parameters",)                 # 本入口直接读取的文件（缺失即报错）
DAYS = 4.0        # 模拟天数（响应统计 = 末 24 步 = 第 4 天）
SPINUP = 3.0      # SPINUP 丢弃天数（须 < DAYS；batch 聚合窗仍为末 24 步）

sys.path.insert(0, os.path.join(ROOT, "analysis_lib"))
from csvio import read_table, write_csv, fcol  # noqa: E402

# 情景组口径（与 modules/runtime_dat_mod.f90::build_groups 单源一致）
NOX_G = {"NO", "NO2"}
NONVOC = {"SO2", "DMS", "MSA"}          # 在 E_VOC 组内但不是 NMVOC（硫）
VOC_G = {"HCHO", "CH3CHO", "ISOP", "DMS", "MVK", "MACR", "MSA", "SO2",
         "ALK", "OLE", "ARO1", "ARO2", "KET", "TERP", "CSL", "LALK",
         "MEOH", "ETOH", "MTBE", "ESTER", "SESQ", "ETH", "IOLE", "ALD2"}


def model_window(rt_path):
    """读 runtime 的 MODEL 行 → (days, spinup)；缺失返回 (nan, nan)"""
    if not rt_path or not os.path.exists(rt_path):
        return float("nan"), float("nan")
    for ln in open(rt_path):
        t = ln.split()
        if len(t) >= 2 and t[0] == "MODEL":
            d = float(t[1]); sp = float(t[5]) if len(t) >= 6 else 0.0
            return d, sp
    return float("nan"), float("nan")


def emis_totals(rt_path):
    """解析 runtime 的 SPECIES+EMIS 段 → (N0, V0_group, V0_nmv) 24h 均值
    排放总量 [molec/cm3/s]。BOX 日型 24h 均值 = E0*(0.1+0.225*amp)；FLAT = E0。
    V0_group = build_groups 的 E_VOC 组（含 SO2/DMS/MSA）；
    V0_nmv   = 剔除 SO2/DMS/MSA 后的非硫 VOC（更接近"NMVOC"）。"""
    lines = open(rt_path).read().splitlines()
    names, emis, sec, k = None, {}, "", 0
    for ln in lines:
        if ln.startswith("#"):
            continue
        t = ln.split()
        if not t:
            continue
        if t[0] == "SPECIES":
            names = t[1:]
        elif t[0] == "EMIS":
            sec, k = "EMIS", 0
        elif t[0] in ("DRY", "WET", "DIL", "SCEN", "SERIES", "PHOT", "MET", "BC",
                      "IC", "MODEL", "RUN"):
            sec = ""
        elif sec == "EMIS" and t[0] == "E" and names is not None:
            e0, shp, amp = float(t[1]), t[2], float(t[3])
            emis[names[k]] = e0 * ((0.1 + 0.225 * amp) if shp == "BOX" else 1.0)
            k += 1
    n0 = sum(v for s_, v in emis.items() if s_ in NOX_G)
    v0 = sum(v for s_, v in emis.items() if s_ in VOC_G)
    v_nmv = sum(v for s_, v in emis.items() if s_ in VOC_G - NONVOC)
    return n0, v0, v_nmv


def prepare_configs(cfg_dir):
    """把 config_override/ 装配到 cfg_dir（临时目录），再固定 4 天模拟 / 3 天 spinup。

    cfg_dir 只服务于本次 cfgbuild：有效输入的最终记录是 output/runtime_emission.dat，
    不在 output/ 留配置快照目录（config_override/ 本身就是用户编辑入口与输入存档）。"""
    shutil.copytree(OVERRIDE_DIR, cfg_dir,
                    ignore=shutil.ignore_patterns("__pycache__"))
    print(f"  [配置] 输入目录: {OVERRIDE_DIR}")
    mp = os.path.join(cfg_dir, "model.parameters")
    lines = open(mp).read().splitlines()
    lines[0] = f"{DAYS:g}        # 模拟天数（EKMA-emission: 4 天）"
    if len(lines) >= 5:
        lines[4] = f"{SPINUP:g}        # SPINUP 丢弃天数 [d]（响应统计 = 末 24 步 = 第 4 天）"
    else:
        lines.append(f"{SPINUP:g}        # SPINUP 丢弃天数 [d]")
    open(mp, "w").write("\n".join(lines) + "\n")
    return cfg_dir


def compile_runtime(cfg_dir, rt_out):
    """cfgbuild: cfg_dir → rt_out（Fortran 顺序表控运行时）"""
    r = subprocess.run([sys.executable, CFGBUILD, cfg_dir, rt_out, "--mech", MECH],
                       capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(rt_out):
        print(r.stdout[-1500:], r.stderr[-500:], file=sys.stderr)
        raise SystemExit("✗ cfgbuild 编译失败")
    return rt_out


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmax", type=float, default=625.0,
                    help="最大扰动百分比（缺省 625）")
    ap.add_argument("--step", type=float, default=25.0,
                    help="扰动步长百分比（缺省 25）")
    ap.add_argument("--from", dest="from_raw", default=None,
                    help="复用指定 raw CSV（绝对/相对仓库根；不装配配置、"
                         "不编译、不跑模拟，仅重导出 grid/derived）")
    a = ap.parse_args(argv[1:])
    os.makedirs(OUT, exist_ok=True)
    dt = None
    nfail_txt = "沿用原 raw（--from 复用模式）"
    guard_txt = "沿用原 raw（未跑模拟，无从判断）"
    rt_i = None

    if a.from_raw:
        # ── 复用模式: 不重建配置/不编译/不跑模拟；基线总量取自产出该 raw 的 runtime ──
        raw = a.from_raw if os.path.isabs(a.from_raw) else os.path.join(ROOT, a.from_raw)
        if not os.path.exists(raw):
            print(f"✗ --from 文件不存在: {raw}", file=sys.stderr)
            return 1
        cands = [os.path.join(os.path.dirname(raw), "runtime_emission.dat"),
                 os.path.join(OUT, "runtime_emission.dat")]
        rt_i = next((q for q in cands if os.path.exists(q)), None)
        if rt_i is None:
            print("  ⚠ 找不到 runtime_emission.dat → VNr 列留空（无法确定基线总量）",
                  file=sys.stderr)
        else:
            print(f"  [复用] 基线总量取自 {os.path.relpath(rt_i, ROOT) if rt_i.startswith(ROOT) else rt_i}")
        print(f"═══ EKMA-emission [复用] "
              f"{os.path.relpath(raw, ROOT) if raw.startswith(ROOT) else raw}"
              f" → grid/derived（不装配配置、不编译、不跑模拟）═══")
    else:
        if a.pmax <= 0 or a.step <= 0 or a.step > a.pmax:
            print("✗ --pmax/--step 非法", file=sys.stderr)
            return 2
        if not os.path.isdir(OVERRIDE_DIR):
            print(f"缺少完整输入目录: {OVERRIDE_DIR}（本程序唯一输入源，无全局配置层）",
                  file=sys.stderr)
            return 1
        miss = [f for f in REQUIRED_INPUT if not os.path.isfile(os.path.join(OVERRIDE_DIR, f))]
        if miss:
            print(f"✗ {OVERRIDE_DIR} 缺必要输入文件: {' '.join(miss)}", file=sys.stderr)
            return 1
        for q in (BATCH, CFGBUILD, MECH):
            if not os.path.exists(q):
                print(f"缺 {q}\n  → 先 make all；配置由本入口自动编译", file=sys.stderr)
                return 1
        n = int(round(a.pmax / a.step)) + 1          # 含 0% 与 pmax%
        fmin, fmax = 1.0, 1.0 + a.pmax / 100.0       # 0% → 1.00
        # ── ① 临时装配配置 + 口径改写（4 天 / 3 天 spinup）+ 编译 ──
        cfgroot = tempfile.mkdtemp(prefix="ekma_emission_cfg_")
        cfgi = os.path.join(cfgroot, "cfg")       # copytree 要求目标目录尚不存在
        try:
            prepare_configs(cfgi)
            rt_i = compile_runtime(cfgi, os.path.join(OUT, "runtime_emission.dat"))
        finally:
            shutil.rmtree(cfgroot, ignore_errors=True)
        print("  [配置] 装配完成 → output/runtime_emission.dat（有效输入记录，不在 output 留快照）")

        # ── ② 批量模拟 ──
        raw = os.path.join(OUT, "ekma_emission_grid_raw.csv")
        print(f"═══ EKMA-emission: {n}×{n} = {n*n} 案例  "
              f"(扰动 0–{a.pmax:g}% / 步长 {a.step:g}% ⇒ 因子 {fmin:.2f}–{fmax:.2f}, "
              f"线性网格, {DAYS:g} 天 / spinup {SPINUP:g} 天) ═══")
        t0 = time.time()
        # 参数: ng outf rtfile jfac dt fmin fmax EMIS g2 LIN
        cmd = [BATCH, str(n), raw, rt_i, "1.0", "0.0",
               f"{fmin:.6g}", f"{fmax:.6g}", "EMIS", "1", "LIN"]
        r = subprocess.run(cmd, cwd=os.path.join(ROOT, "facsimile_impl"),
                           capture_output=True, text=True)
        dt = time.time() - t0
        tail = [l for l in r.stdout.splitlines()
                if "FAC-BATCH" in l or "FATAL" in l or "警告" in l or "氮守恒" in l]
        print("  " + (" | ".join(tail) if tail else "batch done"))
        _mg = [l.strip() for l in r.stdout.splitlines() if "[氮守恒 FAIL]" in l]
        guard_txt = (f"违例 {len(_mg)} 例: " + " ; ".join(_mg[:3])) if _mg else "0 违例"
        _mf = re.search(r"失败 (\d+)", r.stdout)
        _md = re.search(r"降级 (\d+)", r.stdout)
        nfail_txt = (f"失败 {_mf.group(1) if _mf else '?'} / "
                     f"降级 {_md.group(1) if _md else '?'}")
        if r.returncode != 0:
            print(r.stdout[-800:], r.stderr[-800:], file=sys.stderr)
            print(f"  ✗ 批量模拟非零退出（{nfail_txt}；{guard_txt}）"
                  f"—— 网格不完整，拒绝出表", file=sys.stderr)
            return 1
        if _md and int(_md.group(1)) > 0:
            print(f"  ⚠ 有 {_md.group(1)} 个降级案例（窗口不完整）—— 已记入 run_meta")

    # ── 规整：因子 → 扰动百分比（网格范围/步长由 raw 数据推断，两种模式同源）──
    import csv as _csv
    _csv.field_size_limit(16 * 1024 * 1024)
    hdr, D = read_table(raw)
    fx, fv = fcol(hdr, D, "fnox"), fcol(hdr, D, "fvoc")
    keep = ["O3max", "O3mean", "HCHOmax", "NO2mean", "FNR_mean",
            "H2O2HNO3", "O3NOy", "HCHONOy", "O3NOz", "NOy_ppb", "LnOxQ",
            "PH2O2PHNO3", "OPE", "nst"]
    keep = [k for k in keep if k in hdr]
    cols = [fcol(hdr, D, k) for k in keep]
    pn = (fx - 1.0) * 100.0          # NOx 扰动 %
    pv = (fv - 1.0) * 100.0          # VOC 扰动 %
    ux = sorted({round(v, 4) for v in pn})
    uy = sorted({round(v, 4) for v in pv})
    sx = (ux[1] - ux[0]) if len(ux) > 1 else 0.0
    sy = (uy[1] - uy[0]) if len(uy) > 1 else 0.0
    if rt_i is not None:
        N0, V0g, V0n = emis_totals(rt_i)   # 基线总量取自实际产出该 raw 的 runtime
        vnr_g = [(fv[i] * V0g) / (fx[i] * N0) for i in range(len(fx))]
        vnr_n = [(fv[i] * V0n) / (fx[i] * N0) for i in range(len(fx))]
    else:
        N0 = V0g = V0n = float("nan")
        vnr_g = vnr_n = [None] * len(fx)
    rows = []
    for i in range(len(fx)):
        rows.append([f"{pn[i]:.1f}", f"{pv[i]:.1f}"] +
                    [f"{c[i]:.6g}" for c in cols] +
                    [(f"{vnr_g[i]:.4g}" if vnr_g[i] is not None else ""),
                     (f"{vnr_n[i]:.4g}" if vnr_n[i] is not None else "")])
    grid_csv = os.path.join(OUT, "ekma_emission_grid.csv")
    out_names = [("O3max", "O3max_ppb"), ("O3mean", "O3mean_ppb")]
    hdr_out = (["pert_nox_pct", "pert_voc_pct"] +
               [dict(out_names).get(k, k) for k in keep] +
               ["VNr_group", "VNr_nmv"])
    write_csv(grid_csv, hdr_out, rows,
              note="\n".join([
                  f"EKMA-emission | S1 排放正交扰动响应面（{len(ux)}×{len(uy)}={len(rows)} 案例）",
                  f"NOx 扰动 0–{ux[-1]:g}% / 步长 {sx:g}% ⇒ 因子 1.00–{1.0 + ux[-1] / 100.0:.2f}",
                  f"VOC 扰动 0–{uy[-1]:g}% / 步长 {sy:g}% ⇒ 因子 1.00–{1.0 + uy[-1] / 100.0:.2f}",
                  "pert_nox_pct: NOx 排放组 (E_NOx = NO,NO2) 相对基线的扰动 %  ← y 轴",
                  "pert_voc_pct: E_VOC 排放组相对基线的扰动 %  ← x 轴",
                  "VNr_group = (1+VOC%)/(1+NOx%) × V0_group/N0 —— **摩尔比**（集总类无单一摩尔质量，不可当质量比）",
                  "VNr_nmv   = 同上但 V0 剔除 SO2/DMS/MSA 的硫（更接近 NMVOC）",
                  "固定条件: mechanism V2.5+ / 本目录 config_override 的气象与初值",
                  f"时长: 模拟 {DAYS:g} 天 / SPINUP {SPINUP:g} 天",
                  "O3max_ppb = **末 24 个输出步窗口内最大值**（= 第 4 天，不是 8h 日最大）；"
                  "其余列为同窗口均值/比值"]))

    # ── 派生表: PR（峰值率）与 VNr（排放 VOC/NOx 比）——每 NMVOC 列一行 ──
    o3 = fcol(hdr, D, "O3max")
    o3_at = {(round((fx[i] - 1.0) * 100.0, 1), round((fv[i] - 1.0) * 100.0, 1)): o3[i]
             for i in range(len(fx))}
    o3bau = o3_at.get((0.0, 0.0), float("nan"))
    by_voc = {}
    for i in range(len(fx)):
        by_voc.setdefault(round((fv[i] - 1.0) * 100.0, 1), []).append(
            (round((fx[i] - 1.0) * 100.0, 1), o3[i]))
    drows = []
    for voc in sorted(by_voc):
        lst = sorted(by_voc[voc])
        k = max(range(len(lst)), key=lambda t: lst[t][1])
        nox_r, o3pk = lst[k]
        o3_base_nox = dict(lst).get(0.0, float("nan"))   # 0% 扰动 = 基线排放, 非零排放
        vg = ((1.0 + voc / 100.0) * V0g / ((1.0 + nox_r / 100.0) * N0)
              if rt_i is not None else float("nan"))
        vn = ((1.0 + voc / 100.0) * V0n / ((1.0 + nox_r / 100.0) * N0)
              if rt_i is not None else float("nan"))
        edge = 1 if nox_r >= ux[-1] - 1e-9 else 0     # 屋脊顶到采样上界 = 峰值未覆盖
        drows.append([f"{voc:.1f}", f"{o3pk:.4f}", f"{o3pk / o3bau:.4f}",
                      f"{o3pk / o3_base_nox:.4f}", f"{nox_r:.1f}",
                      f"{vg:.4f}" if vg == vg else "",
                      f"{vn:.4f}" if vn == vn else "",
                      edge, f"{o3_base_nox:.4f}"])
    write_csv(os.path.join(OUT, "ekma_emission_derived.csv"),
              ["pert_voc_pct", "O3peak_ppb", "PR_bau", "PR_vs_baseNOx",
               "ridge_pert_nox_pct", "VNr_ridge_group", "VNr_ridge_nmv",
               "ridge_at_edge", "O3_at_baseNOx_ppb"], drows,
              note="\n".join([
                  "EKMA-emission | 派生量: 峰值率 PR 与屋脊点组成比 VNr（每 VOC 列）",
                  "O3peak_ppb = 该 VOC 扰动列沿 NOx 轴的 O3max 峰值（屋脊点）",
                  "PR_bau  = 峰值率① = O3peak / O3max(基线 0%,0%)",
                  "PR_vs_baseNOx = 峰值率② = O3peak / O3max(同列 NOx=0% 端; 0% = 基线)",
                  "⚠ 0% 扰动 = 基线排放（不是零排放）→ 该比值是「从基线到峰值」的放大，列名为历史命名；PR_bau 才是「相对 (0%,0%) 基线」的峰值比",
                  "VNr_ridge_group = 屋脊点 (1+VOC%)/(1+NOx%) × V0_group/N0（**摩尔比**，含硫）",
                  "VNr_ridge_nmv   = 同上但剔除 SO2/DMS/MSA",
                  "ridge_at_edge=1 → 该列峰值落在 NOx 采样上界，真实峰可能在界外（外推无效）",
                  (f"基线总量: NOx N0={N0:.4g}, V0_group={V0g:.4g}, V0_nmv={V0n:.4g} molec/cm3/s"
                   if rt_i is not None else "基线总量: 未知（复用模式未找到 runtime）"),
                  "经典 EKMA: 屋脊处组成比近似常数 = 临界排放比；本节按归档数据复核其恒定性"]))

    # ── 屋脊表（run.py 是唯一 CSV 导出者；plot.py 只画图）──
    rid = []
    for uv in uy:
        col = [(a, c) for a, b, c in zip(pn, pv, o3) if round(b, 1) == uv]
        if col:
            noxv, o3v = max(col, key=lambda t: t[1])
            rid.append([f"{uv:.1f}", f"{noxv:.1f}", f"{o3v:.4f}"])
    p_ridge = write_csv(os.path.join(OUT, "ekma_emission_ridge.csv"),
                        ["pert_voc_pct", "ridge_pert_nox_pct", "ridge_O3max_ppb"],
                        rid,
                        note="EKMA-emission | 屋脊：每 VOC 扰动列的 O3max 峰对应 NOx 扰动")

    mech_md5 = (hashlib.md5(open(MECH, "rb").read()).hexdigest()[:12]
                if os.path.exists(MECH) else "n/a")
    _d, _sp = model_window(rt_i) if rt_i else (float("nan"), float("nan"))
    with open(os.path.join(OUT, "run_meta.txt"), "w") as f:
        f.write(f"EKMA-emission 运行元数据\n 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                f" 模式: {'复用 raw（--from，未跑模拟）' if a.from_raw else '完整运行'}\n"
                + (f" 来源 raw: {raw}\n" if a.from_raw else "")
                + f" 网格: {len(ux)}×{len(uy)} = {len(rows)} 案例 "
                  f"(NOx 0–{ux[-1]:g}%/{sx:g}%, VOC 0–{uy[-1]:g}%/{sy:g}%)\n"
                f" 时长: {_d:g} 天 / SPINUP {_sp:g} 天（响应统计 = 末 24 步 = 最后一天）\n"
                f" 耗时: {'—（复用模式未计时）' if dt is None else f'{dt:.1f} s'}\n"
                f" 机制: {os.path.relpath(MECH, ROOT)} md5={mech_md5}\n"
                f" 基线总量 runtime: "
                f"{os.path.relpath(rt_i, ROOT) if rt_i and rt_i.startswith(ROOT) else rt_i}\n"
                f" 输入目录: {OVERRIDE_DIR}\n"
                f" 有效输入记录: output/runtime_emission.dat"
                f"（IC/EMIS/MODEL/MET 解析值；不在 output 留配置快照目录）\n"
                f" 失败/降级: {nfail_txt}\n"
                f" 氮守恒不变量: {guard_txt}\n"
                f" 侧车单位: {os.path.basename(raw)}_ipr.csv/_irr.csv = 末 24 步**平均倾向** "
                f"[molec/cm3/s]；fac_v3 的三件套是窗口**积分** [molec/cm3]，差 dtout 倍\n")

    o3 = fcol(hdr, D, "O3max") if "O3max" in hdr else []
    print(f"  产物: {os.path.relpath(grid_csv, ROOT)}, {os.path.relpath(p_ridge, ROOT)}")
    if len(o3):
        print(f"  O3max 范围 {min(o3):.2f} – {max(o3):.2f} ppb")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
