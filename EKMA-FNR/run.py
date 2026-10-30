#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
EKMA-FNR/run.py —— S2 FNR 响应面（初始 HCHO × 初始 NOx 正交扰动，EKMA-FNR 口径）

实验定义（用户确认口径）:
  · 固定: 机制 .fac（V2.5+ 冻结）+ 气象（config_override 的逐时序列，仅光解/物理环境生效）
  · 排放: **关闭**（emissions.config 清空 → 纯初值化学）
  · 变量: **初始 NOx**（ic1 = NOX 组：**NO 与 NO₂ 同因子缩放 → 保持 NO₂:NO 固定比例**）
          × **初始 HCHO 浓度**（只缩 HCHO）
  · 基线: 全部取自 config_override/initialConcentrations.config
          （当前试验基线 NO 1 / NO₂ 3 / HCHO 6 ppb ⇒ NO₂:NO = 3:1；O₃ 不缩放）
  · 扰动范围: 默认 0%（无扰动）→ 1000%（+1000%），**100% 步长**
      → 因子 1.00, 2.00, …, 11.00 共 11 点/轴（--pmax/--step 可改）
  · 网格: 默认 11 × 11 = 121 个独立 0-D 箱案例，每案例 1 天纯初值化学

与 S1 (EKMA-emission) 的区别:
  S1 扰动的是其 config_override/emissions.config（持续排放）；本实验扰动的是
  本目录 config_override/initialConcentrations.config（初始浓度），且排放关闭——
  这是 Dodge 1977 经典 EKMA 的"初始浓度空间"口径。

自包含性: 本目录不 import 任何编排层代码；配置读取/编译/批量调用
全部在本文件内实现。输入 = 本目录 config_override/（仓库无全局配置层）。
共享依赖仅: 机制 .fac、inputproc/cfgbuild.py（编译器）、
facsimile_impl/fac_batch_omp（执行器）、analysis_lib/csvio（CSV 读写）。

输出（本目录 output/，全部由本文件写出）:
  runtime_fnr.dat                本轮有效输入记录（cfgbuild 产物）
  fnr_ic_base.csv                基线初值存档（含每物种的缩放角色）
  ekma_fnr_grid_raw.csv          fac_batch_omp 原始 24 列（n² 行）
  ekma_fnr_grid_raw.csv_{ipr,irr,state}.csv   执行器侧车（**平均倾向** [molec/cm3/s]，
                                 与 fac_v3 的窗口积分口径差 dtout 倍）
  ekma_fnr_grid.csv              规整表（pert_hcho_pct, pert_no2_pct, 指标列, FNR_init）
  ekma_fnr_ridge.csv             每 HCHO 列的屋脊点（NO₂ 扰动 + O3max）
  ekma_fnr_derived.csv           PR_bau / PR_vs_baseNO2 / FNR_init_ridge 等派生量
  run_meta.txt                   运行元数据（网格/耗时/失败/降级/氮守恒/侧车单位）

用法:
  python EKMA-FNR/run.py                # 默认 0–1000% / 100%（11×11 = 121 案例）
  python EKMA-FNR/run.py --pmax 100 --step 20
  python EKMA-FNR/run.py --from path/to/ekma_fnr_grid_raw.csv   # 只重导出 grid
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
CFGBUILD = os.path.join(ROOT, "inputproc", "cfgbuild.py")
MECH = os.path.join(ROOT, "mechanism", "mechanism_3_full_chain.fac")
sys.path.insert(0, os.path.join(ROOT, "analysis_lib"))
from csvio import read_table, read_rows, write_csv, fcol  # noqa: E402

PY = sys.executable

DAYS = 1.0            # 纯初值化学: 1 天（排放关闭）
EMPTY = "# EKMA-FNR: 排放关闭（纯初值化学）\n"
OVERRIDE_DIR = os.path.join(HERE, "config_override")
# ic1 轴掩膜：NOX = NO + NO₂ 同因子缩放（保持 IC 基线里的 NO₂:NO 比例不变）
IC1_SPEC = "NOX"
IC2_SPEC = "HCHO"
# 本程序直接读取的输入文件（缺失即报错：无根配置回退、无内置常量兜底）
REQUIRED_INPUT = ("initialConcentrations.config", "model.parameters")


def prepare_configs(cfg_dir):
    """把 config_override/ 装配到 cfg_dir（临时目录），再固定排放关闭及 1 天／0 spinup。

    cfg_dir 只服务于本次 cfgbuild：有效输入的最终记录是 output/runtime_fnr.dat，
    不在 output/ 留配置快照目录（config_override/ 本身就是用户编辑入口与输入存档）。"""
    shutil.copytree(OVERRIDE_DIR, cfg_dir,
                    ignore=shutil.ignore_patterns("__pycache__"))
    print(f"  [配置] 输入目录: {OVERRIDE_DIR}")
    ic_source = "config_override/initialConcentrations.config"
    with open(os.path.join(cfg_dir, "initialConcentrations.config")) as f:
        ic_lines_used = f.read().splitlines()
    print(f"  [IC 基线] {ic_source}")
    open(os.path.join(cfg_dir, "emissions.config"), "w").write(EMPTY)
    mp = os.path.join(cfg_dir, "model.parameters")
    lines = open(mp).read().splitlines()
    lines[0] = f"{DAYS:g}        # 模拟天数（EKMA-FNR: 1 天纯初值化学）"
    if len(lines) >= 5:
        lines[4] = "0.0        # SPINUP 丢弃天数 [d]（1 天模拟须 0）"
    else:
        lines.append("0.0        # SPINUP 丢弃天数 [d]")
    open(mp, "w").write("\n".join(lines) + "\n")
    return cfg_dir, ic_source, ic_lines_used


def ic_baseline_ppb(lines, csv_path):
    """→ (HCHO_ppb, NO2_ppb, 来源)。优先用本轮 config_override 装配值，否则读归档 fnr_ic_base.csv。"""
    def to_ppb(v, unit):
        return v / 2.46e10 if unit in ("cm3", "") else (v * 1e3 if unit == "ppm" else v)

    if lines:
        d = {}
        for ln in lines:
            t = ln.split("#")[0].strip()
            if not t:
                continue
            q = t.split()
            if len(q) < 2:
                continue
            unit = q[2].lower() if len(q) > 2 else "cm3"
            d[q[0]] = to_ppb(float(q[1]), unit)
        if "HCHO" in d and "NO2" in d:
            return d["HCHO"], d["NO2"], "本轮 config_override 装配"
    if os.path.exists(csv_path):
        h, rows = read_rows(csv_path)
        d = {r[h.index("species")]: float(r[h.index("value_cm3")]) / 2.46e10
             for r in rows}
        if "HCHO" in d and "NO2" in d:
            return d["HCHO"], d["NO2"], "fnr_ic_base.csv（归档）"
    print("  ⚠ 找不到 HCHO/NO₂ 基线 → FNR_init 记 nan（检查 config_override/ 或归档）",
          file=sys.stderr)
    return float("nan"), float("nan"), "缺失"


def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmax", type=float, default=1000.0,
                    help="最大扰动百分比（缺省 1000；250/500 实测窗内无完整屋脊或降级更多，见 README §8.6）")
    ap.add_argument("--step", type=float, default=100.0,
                    help="扰动步长百分比（缺省 100）")
    ap.add_argument("--from", dest="from_raw", default=None,
                    help="复用指定 raw CSV（绝对路径或相对仓库根；不重建配置、不编译、"
                         "不跑模拟，仅重导出 grid）")
    a = ap.parse_args(argv[1:])
    from_raw = a.from_raw

    os.makedirs(OUT, exist_ok=True)
    if from_raw:
        raw = from_raw if os.path.isabs(from_raw) else os.path.join(ROOT, from_raw)
        if not os.path.exists(raw):
            print(f"✗ --from 文件不存在: {raw}", file=sys.stderr)
            return 1
    else:
        if a.pmax <= 0 or a.step <= 0 or a.step > a.pmax:
            print("✗ --pmax/--step 非法", file=sys.stderr)
            return 2
        if not os.path.isdir(OVERRIDE_DIR):
            print(f"缺少完整输入目录: {OVERRIDE_DIR}（本程序无全局配置回退）", file=sys.stderr)
            return 1
        miss = [f for f in REQUIRED_INPUT if not os.path.isfile(os.path.join(OVERRIDE_DIR, f))]
        if miss:
            print(f"✗ {OVERRIDE_DIR} 缺必要输入文件: {' '.join(miss)}", file=sys.stderr)
            return 1
        for p in (BATCH, CFGBUILD, MECH):
            if not os.path.exists(p):
                print(f"✗ 缺 {p}", file=sys.stderr)
                return 1
        n = int(round(a.pmax / a.step)) + 1
        fmin, fmax = 1.0, 1.0 + a.pmax / 100.0   # 0% → 1.00 ; 250% → 3.50

    # ── ① 临时装配配置 + 编译（有效输入 = output/runtime_fnr.dat，不在 output 留快照目录）──
    #    --from 复用模式跳过本步：不装配、不编译、不跑模拟
    rt = os.path.join(OUT, "runtime_fnr.dat")
    mech_md5 = (hashlib.md5(open(MECH, "rb").read()).hexdigest()[:12]
                if os.path.exists(MECH) else "n/a")
    if not from_raw:
        cfgroot = tempfile.mkdtemp(prefix="ekma_fnr_cfg_")
        cfgi = os.path.join(cfgroot, "cfg")       # copytree 要求目标目录尚不存在
        cfgi, ic_source, ic_lines_used = prepare_configs(cfgi)
        r = subprocess.run([PY, CFGBUILD, cfgi, rt, "--mech", MECH],
                           capture_output=True, text=True)
        shutil.rmtree(cfgroot, ignore_errors=True)
        if r.returncode != 0:
            print(r.stdout[-800:], r.stderr[-800:], file=sys.stderr)
            print("✗ cfgbuild 失败（装配目录已清理；有效输入请查 output/runtime_fnr.dat）",
                  file=sys.stderr)
            return 1
        _ic_rows = []
        for ln in ic_lines_used:
            ln0 = ln.split("#")[0].strip()
            if not ln0:
                continue
            parts = ln0.split()
            sp, val = parts[0], float(parts[1])
            unit = parts[2].lower() if len(parts) > 2 else "cm3"
            fac = {"cm3": 1.0, "ppb": 2.46e10, "ppm": 2.46e13}.get(unit)
            if fac is None:
                print(f"  ✗ IC 基线未知单位 '{unit}' ({sp})", file=sys.stderr)
                return 2
            role = ("perturbed by ic1 (NOX: NO+NO2 same factor)" if sp in ("NO", "NO2")
                    else "perturbed by ic2 (HCHO)" if sp == IC2_SPEC else "fixed baseline")
            _ic_rows.append([sp, f"{val * fac:.6g}", unit, role])
        write_csv(os.path.join(OUT, "fnr_ic_base.csv"),
                  ["species", "value_cm3", "unit_as_written", "scaling"],
                  _ic_rows,
                  note="EKMA-FNR | 基线初值（ic1=NOX 组：NO 与 NO₂ 同因子缩放，"
                       "保持 NO₂:NO 比例; ic2=HCHO; O₃ 不缩放）\n"
                       f"网格: 0–{a.pmax:g}% / {a.step:g}%\n"
                       f"基线来源: {ic_source}\n"
                       "value_cm3 = 按 unit_as_written 换算后的 molec/cm3（ppb ×2.46e10）")
        ic_note = "本目录 IC 基线存档已随本次运行重建"
    else:
        # 复用: IC 基线存档随 raw 一起取（源目录有则复制，无则保留本地并明示）
        src_ic = os.path.join(os.path.dirname(raw), "fnr_ic_base.csv")
        dst_ic = os.path.join(OUT, "fnr_ic_base.csv")
        if os.path.abspath(src_ic) == os.path.abspath(dst_ic):
            ic_note = "IC 基线存档即本地 output 内（--from 指向本目录 raw，无需复制）"
        elif os.path.exists(src_ic):
            shutil.copy(src_ic, dst_ic)
            ic_note = "IC 基线存档随 raw 复制"
        else:
            ic_note = "源目录无 fnr_ic_base.csv → 本地存档沿用（未更新，基线可能不符）"

    # ── ② 批量模拟: n×n，IC 模式，物种对 NOX:HCHO（ic1 同因子缩 NO+NO₂；--from 时整段跳过）──
    dt = None
    r = None
    nfail_txt = "沿用原 raw（--from 复用模式）"
    guard_txt = "沿用原 raw（未跑模拟，无从判断）"
    if not from_raw:
        raw = os.path.join(OUT, "ekma_fnr_grid_raw.csv")
        print(f"═══ EKMA-FNR: {n}×{n} = {n*n} 案例  (NOx(=NO+NO₂)/HCHO 初值扰动 "
              f"0–{a.pmax:g}% / {a.step:g}% ⇒ 因子 {fmin:.2f}–{fmax:.2f}, "
              f"线性网格, 排放关闭, {DAYS:g} 天) ═══")
        t0 = time.time()
        cmd = [BATCH, str(n), raw, rt, "1.0", "0.0",
               f"{fmin:.6g}", f"{fmax:.6g}", "IC", f"{IC1_SPEC}:{IC2_SPEC}", "LIN"]
        r = subprocess.run(cmd, cwd=os.path.join(ROOT, "facsimile_impl"),
                           capture_output=True, text=True)
        dt = time.time() - t0
        # 执行器的关键行一律捕获：违例原始行 [氮守恒 FAIL] 不含 FAC-BATCH/FATAL 关键词，
        # 不显式列进过滤器就会被丢掉（计数本身已进收尾行并触发 stop 1）
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
    else:
        print(f"═══ EKMA-FNR [复用] "
              f"{os.path.relpath(raw, ROOT) if raw.startswith(ROOT) else raw}"
              f" → 仅重导出 grid/派生表（不重建配置、不编译、不跑模拟）═══")

    # ── ③ 规整: 因子 → 扰动百分比（网格范围/步长由 raw 数据推断，两种模式同源）──
    hdr, D = read_table(raw)
    fx, fv = fcol(hdr, D, "fnox"), fcol(hdr, D, "fvoc")
    keep = ["O3max", "O3mean", "HCHOmax", "NO2mean", "FNR_mean", "H2O2HNO3",
            "O3NOy", "HCHONOy", "O3NOz", "NOy_ppb", "LnOxQ", "PH2O2PHNO3",
            "OPE", "FNR_pm", "H2O2HNO3_pm", "O3NOy_pm", "HCHONOy_pm",
            "O3NOz_pm", "NOy_pm", "nst"]
    keep = [k for k in keep if k in hdr]
    cols = [fcol(hdr, D, k) for k in keep]
    H0, N0i, base_src = ic_baseline_ppb(
        None if from_raw else ic_lines_used, os.path.join(OUT, "fnr_ic_base.csv"))
    print(f"  [基线] HCHO={H0:g} ppb, NO₂={N0i:g} ppb（{base_src}）")
    pn = (fx - 1.0) * 100.0          # NO₂ 扰动 %
    ph = (fv - 1.0) * 100.0          # HCHO 扰动 %
    ux = sorted({round(v, 4) for v in pn})
    uy = sorted({round(v, 4) for v in ph})
    sx = (ux[1] - ux[0]) if len(ux) > 1 else 0.0
    sy = (uy[1] - uy[0]) if len(uy) > 1 else 0.0
    rows = [[f"{ph[i]:.1f}", f"{pn[i]:.1f}"] +
            [f"{c[i]:.6g}" for c in cols] +
            [f"{(1.0 + (fv[i] - 1.0)) * H0 / ((1.0 + (fx[i] - 1.0)) * N0i):.4f}"]
            for i in range(len(fx))]
    grid_csv = os.path.join(OUT, "ekma_fnr_grid.csv")
    names = {"O3max": "O3max_ppb", "O3mean": "O3mean_ppb"}
    write_csv(grid_csv,
              ["pert_hcho_pct", "pert_no2_pct"] + [names.get(k, k) for k in keep] + ["FNR_init"],
              rows,
              note="\n".join([
                  f"EKMA-FNR | S2 初始浓度正交扰动响应面（{len(ux)}×{len(uy)}={len(rows)} 案例）",
                  f"NO₂ 扰动 0–{ux[-1]:g}% / 步长 {sx:g}% ⇒ 因子 1.00–{1.0 + ux[-1] / 100.0:.2f}"
                  f"（ic1 = NOX 组：NO 与 NO₂ 同因子，NO₂:NO 比例恒定）",
                  f"HCHO 扰动 0–{uy[-1]:g}% / 步长 {sy:g}% ⇒ 因子 1.00–{1.0 + uy[-1] / 100.0:.2f}",
                  f"pert_hcho_pct: 初始 HCHO 相对基线 ({H0:g} ppb) 的扰动 %  ← 横轴（VOC 轴，FNR 分子）",
                  f"pert_no2_pct:  初始 NO₂ 相对基线 ({N0i:g} ppb) 的扰动 %   ← 纵轴（NOx 轴，FNR 分母）"
                  f"（= NOx 组因子；NO 同步缩放至同比例）",
                  f"基线来源: {base_src}（FNR_init = (1+HCHO%)/(1+NO₂%) × H0/N0）",
                  "O₃ 不缩放（因变量）；其余物种保持基线",
                  "排放关闭；1 天纯初值化学；气象 = config_override 的逐时序列（仅第 1 天生效）",
                  "FNR_mean = 窗口内 HCHO/NO₂ 均值比（Sillman 判据）"]))

    # ── ④ 屋脊与派生表（run.py 是唯一 CSV 导出者；plot.py 只画图）──
    Z = {(round(a, 4), round(b, 4)): c
         for a, b, c in zip(pn, ph, fcol(hdr, D, "O3max"))}
    o3_00 = Z.get((0.0, 0.0), float("nan"))
    ridge_no2, drows = [], []
    for yv in uy:
        col = [(x, Z[(x, yv)]) for x in ux if (x, yv) in Z]
        if not col:
            continue
        no2_r, o3pk = max(col, key=lambda t: t[1])
        o3_base = dict(col).get(0.0, float("nan"))
        ridge_no2.append([f"{yv:.1f}", f"{no2_r:.1f}", f"{o3pk:.4f}"])
        drows.append([f"{yv:.1f}", f"{o3pk:.4f}", f"{o3pk / o3_00:.4f}",
                      f"{o3pk / o3_base:.4f}", f"{no2_r:.1f}",
                      f"{(1.0 + yv / 100.0) * H0 / ((1.0 + no2_r / 100.0) * N0i):.4f}",
                      f"{o3_base:.4f}"])
    p_ridge = write_csv(os.path.join(OUT, "ekma_fnr_ridge.csv"),
                        ["pert_hcho_pct", "ridge_pert_no2_pct", "ridge_O3max_ppb"],
                        ridge_no2,
                        note="EKMA-FNR | 屋脊：每 HCHO 扰动列的 O3max 峰对应 NO₂ 扰动")
    p_der = write_csv(os.path.join(OUT, "ekma_fnr_derived.csv"),
                      ["pert_hcho_pct", "O3peak_ppb", "PR_bau", "PR_vs_baseNO2",
                       "ridge_pert_no2_pct", "FNR_init_ridge", "O3_at_baseNO2_ppb"],
                      drows,
                      note="\n".join([
                          "EKMA-FNR | 峰值率 PR 与初始组成比（每 HCHO 列一行）",
                          "PR_bau = O3peak / O3max(基线 0%,0%)；"
                          "PR_vs_baseNO2 = O3peak / O3max(同列 NO₂=0% 端; 0% = 基线初值)",
                          f"FNR_init_ridge = (1+HCHO%)/(1+NO₂%) × H0/N0，"
                          f"H0={H0:g}/N0={N0i:g} ppb（{base_src}）",
                          "O3max = 末 24 个输出步窗口最大（1 天模拟即全天）"]))

    with open(os.path.join(OUT, "run_meta.txt"), "w") as f:
        f.write(f"EKMA-FNR 运行元数据\n 时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
                f" 模式: {'复用 raw（--from，未跑模拟）' if from_raw else '完整运行'}\n"
                + (f" 来源 raw: {raw}\n" if from_raw else "")
                + f" 网格: {len(ux)}×{len(uy)} = {len(rows)} 案例 "
                  f"(NO₂ 0–{ux[-1]:g}%/{sx:g}%, HCHO 0–{uy[-1]:g}%/{sy:g}%)\n"
                f" 耗时: {'—（复用模式未计时）' if dt is None else f'{dt:.1f} s'}\n"
                f" 机制 md5={mech_md5}\n"
                f" IC 缩放: ic1=NOX（NO+NO₂ 同因子，保持 NO₂:NO 比例）, ic2=HCHO\n"
                f" IC 基线存档: {ic_note}\n"
                f" FNR 基线: HCHO={H0:g} ppb, NO₂={N0i:g} ppb（{base_src}）\n"
                f" 输入目录: {OVERRIDE_DIR}\n"
                f" 有效输入记录: {rt}（IC/EMIS/MODEL/MET 解析值；不在 output 留配置快照目录）\n"
                f" 排放: 关闭 | 时长: {DAYS:g} 天 | spinup: 0\n"
                f" 失败案例: {nfail_txt}\n"
                f" 氮守恒不变量: {guard_txt}\n"
                f" 侧车单位: {os.path.basename(raw)}_ipr.csv/_irr.csv = 末 24 步**平均倾向** "
                f"[molec/cm3/s]；fac_v3 的三件套是窗口**积分** [molec/cm3]，差 dtout 倍\n")

    o3 = fcol(hdr, D, "O3max")
    print(f"  网格: {len(ux)}×{len(uy)} = {len(rows)} 案例 "
          f"(NO₂ 0–{ux[-1]:g}%/{sx:g}%, HCHO 0–{uy[-1]:g}%/{sy:g}%)")
    print(f"  产物: {os.path.relpath(grid_csv, ROOT)}, "
          f"{os.path.relpath(p_ridge, ROOT)}, {os.path.relpath(p_der, ROOT)}")
    if len(o3):
        print(f"  O3max 范围 {min(o3):.2f} – {max(o3):.2f} ppb")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
