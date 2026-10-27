#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
_runner.py —— run/ 入口脚本的**仿真编排库**（是库不是入口；只跑模拟，不含绘图与分析）

分层（2026-10-27 拆分）:
  · 本模块 = 编排层：路径常量 / 前置检查 / 三个驱动封装 / 输入变体 clone_input /
    manifest 溯源 / 终端简表。
  · CSV 读写与 Fortran 数值兜底 = analysis_lib/csvio.py（**唯一实现**）。
    本模块把 read_table/col/last24_mean/write_csv/PPB 原样转发，脚本 import 无需改动。

设计约定:
  · 输入假设 = "runtime_v3.dat 已由 cfgbuild 编译完成"（本模块只做前置检查，不编译）
  · 输出统一落到 output/run_<X>/ 下，全是 CSV/TXT，可直接喂给 matplotlib/pandas/Excel
  · 绝不删除任何已有数据（幂等覆盖本实验自己的文件，不碰 output/ 其它内容）
  · 每个脚本结束时写 manifest.json（运行时间/参数/产物清单）便于溯源

用法（任何 run_X.py 顶部，按需取用）:
    from _runner import (ROOT, OUTROOT, require_prereq, run_fac_v3,
                         read_table, col, last24_mean, write_csv, manifest)
"""
import os
import sys
import json
import time
import subprocess
import numpy as np

# ── 路径（相对本文件定位，可被环境变量覆盖）──
RUN_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(RUN_DIR))          # ozone_lab/
OUTROOT = os.environ.get("OZONE_OUT", os.path.join(ROOT, "output"))

PY = sys.executable
CONDA_ENV = os.environ.get("CONDA_ENV", os.path.join(os.path.expanduser("~"),
                                                    "miniconda3", "envs", "venv"))
if not os.path.exists(PY) and os.path.exists(os.path.join(CONDA_ENV, "bin", "python")):
    PY = os.path.join(CONDA_ENV, "bin", "python")

FACV3 = os.path.join(ROOT, "facsimile_impl", "fac_v3")
BATCH = os.path.join(ROOT, "facsimile_impl", "fac_batch_omp")
MCM = os.path.join(ROOT, "mcm_city", "mcm_run")
CFGBUILD = os.path.join(ROOT, "inputproc", "cfgbuild.py")
MECH = os.path.join(ROOT, "mechanism", "mechanism_3_full_chain.fac")
RT = os.path.join(ROOT, "runtime_v3.dat")
INPUT = os.path.join(ROOT, "configs")

PPB = 2.46e10

sys.path.insert(0, os.path.join(ROOT, "analysis_lib"))
# CSV 唯一实现（analysis_lib/csvio.py），此处转发以保持 run_* 的 import 面不变
from csvio import PPB, col, last24_mean, read_table, write_csv  # noqa: E402

SCEN_ALL = ["BAU", "LNNOX", "HVOC", "HOT", "CLOUD", "ISOPX"]


# ══════════════════════ 基础设施 ══════════════════════
def log(msg=""):
    print(msg, flush=True)


def ensure_dir(d):
    os.makedirs(d, exist_ok=True)
    return d


def out_dir(name):
    """output/run_<name>/ 并确保存在"""
    return ensure_dir(os.path.join(OUTROOT, name))


def require_prereq(need_batch=False, need_mcm=False, need_rt=True):
    """前置检查：只检查，不自动编译（输入 dat 由用户/上一步准备）"""
    problems = []
    if need_rt and not os.path.exists(RT):
        problems.append(f"缺少 {RT}\n      → 先运行: make runtime_v3.dat")
    if not os.path.exists(FACV3):
        problems.append(f"缺少 {FACV3}\n      → 先运行: make all")
    if need_batch and not os.path.exists(BATCH):
        problems.append(f"缺少 {BATCH}\n      → 先运行: make all")
    if need_mcm and not os.path.exists(MCM):
        problems.append(f"缺少 {MCM}\n      → 先构建 mcm_city/ (见 environment/ENVIRONMENT.md)")
    if problems:
        log("✗ 前置条件不满足:")
        for p in problems:
            log("    - " + p)
        raise SystemExit(1)
    if need_rt:
        sz = os.path.getsize(RT)
        log(f"  [前置OK] runtime_v3.dat ({sz} B) · {os.path.basename(FACV3)}"
            + (f" · {os.path.basename(BATCH)}" if need_batch else "")
            + (f" · {os.path.basename(MCM)}" if need_mcm else ""))


def run(cmd, cwd=None, timeout=None, check=True):
    """执行命令，返回 (rc, stdout, stderr, elapsed)"""
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd, timeout=timeout)
    dt = time.time() - t0
    if check and r.returncode != 0:
        log(f"  ✗ 命令失败 (rc={r.returncode}): {' '.join(map(str, cmd))}")
        log("    stdout: " + (r.stdout or "")[-500:])
        log("    stderr: " + (r.stderr or "")[-500:])
        raise SystemExit(1)
    return r.returncode, r.stdout, r.stderr, dt


# ══════════════════════ 求解器封装 ══════════════════════
def _parse_facv3(out):
    """从 fac_v3 stdout 解析 O3max / nst"""
    o3max, nst = np.nan, -1
    for ln in out.splitlines():
        if "O3max=" in ln:
            try:
                o3max = float(ln.split("O3max=")[1].split()[0])
            except (IndexError, ValueError):
                pass
            for tok in ln.replace(",", " ").split():
                if tok.startswith("nst="):
                    try:
                        nst = int(tok[4:])
                    except ValueError:
                        pass
    return o3max, nst


def run_fac_v3(scenario, pref, rt=RT, verbose=False):
    """跑一个 fac_v3 情景。pref 为输出前缀（绝对路径或相对 facsimile_impl/）"""
    ensure_dir(os.path.dirname(os.path.abspath(pref)) or ".")
    rc, out, err, dt = run([FACV3, scenario, pref, rt], cwd=os.path.join(ROOT, "facsimile_impl"))
    o3max, nst = _parse_facv3(out)
    if verbose:
        log(f"    {scenario}: O3max={o3max:.2f} ppb, nst={nst}, {dt:.2f}s")
    return {"scenario": scenario, "prefix": pref, "o3max": o3max, "nst": nst,
            "elapsed_s": dt, "rc": rc}


def run_batch(ng, pref, rt=RT, jfac_list=None, dt_list=None, verbose=False):
    """跑 fac_batch_omp 参数扫描。返回 stdout 供解析"""
    cmd = [BATCH, str(ng), pref, rt]
    if jfac_list:
        cmd.append(",".join(str(x) for x in jfac_list))
    if dt_list:
        if not jfac_list:
            cmd.append("1.0")
        cmd.append(",".join(str(x) for x in dt_list))
    rc, out, err, dt = run(cmd, cwd=os.path.join(ROOT, "facsimile_impl"))
    if verbose:
        line = [l for l in out.splitlines() if "FAC-BATCH" in l]
        log("    " + (line[0].strip() if line else "batch done"))
    return {"prefix": pref, "ng": ng, "jfac": jfac_list, "dt": dt_list,
            "elapsed_s": dt, "rc": rc, "stdout": out}


def run_mcm(days, pref, fnox=1.0, fvoc=1.0):
    """跑 mcm_run（MCM 对照）。pref 可含路径（已支持 mkdir -p）"""
    ensure_dir(os.path.dirname(os.path.abspath(pref)) or ".")
    rc, out, err, dt = run([MCM, str(days), pref, f"{fnox:.3f}", f"{fvoc:.3f}"],
                           cwd=os.path.join(ROOT, "mcm_city"))
    return {"prefix": pref, "elapsed_s": dt, "rc": rc, "stdout": out}


# ══════════════════════ 输入变体 ══════════════════════
def clone_input(dst, emissions_scale=None, bc_lines=None, tau=None, obs=None,
                extra_constraints=None, days=None, dtout=None):
    """复制 configs → dst，并按需改写排放/bc/约束，然后 cfgbuild 编译。

    参数:
      emissions_scale : {species: factor}
      bc_lines        : 覆盖整个 bc.config 的行列表（如 ["TAU_MIX = 21600", "O3 = file:..."]）
      tau             : 仅设 TAU_MIX（保留原 bc 物种行）
      obs             : {species: (times_h, vals_ppb)} → 生成序列文件并登记为约束
      extra_constraints : {dst_name: source_path} 复制到 constraints/
      days            : 模拟天数；走 cfgbuild --days 覆盖，**不改 model.parameters 文件**
      dtout           : 输出步长 [s]；一般不需改（见 GUIDE §3⑫ 坑①）
    返回: runtime_dat_path
    """
    import shutil
    ensure_dir(dst)
    ensure_dir(os.path.join(dst, "constraints"))
    for f in os.listdir(INPUT):
        p = os.path.join(INPUT, f)
        if os.path.isfile(p):
            shutil.copy(p, dst)
    for f in os.listdir(os.path.join(INPUT, "constraints")):
        shutil.copy(os.path.join(INPUT, "constraints", f),
                    os.path.join(dst, "constraints", f))

    if extra_constraints:
        for name, src in extra_constraints.items():
            shutil.copy(src, os.path.join(dst, "constraints", name))

    if emissions_scale:
        lines = []
        for ln in open(os.path.join(INPUT, "emissions.config")):
            s = ln.split("#")[0].split()
            if s and s[0] in emissions_scale:
                s[1] = f"{float(s[1]) * emissions_scale[s[0]]:.6g}"
                lines.append("  ".join(s) + "\n")
            else:
                lines.append(ln)
        open(os.path.join(dst, "emissions.config"), "w").writelines(lines)

    if bc_lines is not None:
        open(os.path.join(dst, "bc.config"), "w").write("\n".join(bc_lines) + "\n")
    elif tau is not None:
        txt = open(os.path.join(INPUT, "bc.config")).read()
        out = []
        hit = False
        for ln in txt.splitlines():
            if ln.strip().startswith("TAU_MIX"):
                out.append(f"TAU_MIX = {tau}")
                hit = True
            else:
                out.append(ln)
        if not hit:
            out.insert(0, f"TAU_MIX = {tau}")
        open(os.path.join(dst, "bc.config"), "w").write("\n".join(out) + "\n")

    if obs:
        for sp, (times_h, vals) in obs.items():
            fn = f"obs_{sp}.txt"
            with open(os.path.join(dst, "constraints", fn), "w") as f:
                f.write(f"# unit: ppb\n# station: exp-synth  species: {sp}\n# t[s]  value\n")
                for t, v in zip(times_h, vals):
                    f.write(f"{int(round(t * 3600)):10d}  {v:9.4f}\n")
        # 在现有 bc.config 上登记（★ 保留背景边界常数行，只把被约束物种改成 file:）
        src_bc = os.path.join(dst, "bc.config")
        lines = open(src_bc).read().splitlines() if os.path.exists(src_bc) else []
        kept = [ln for ln in lines
                if not any(ln.strip().startswith(sp + " ") or ln.strip().startswith(sp + "=")
                           or ln.strip().startswith(sp + "\t") for sp in obs)]
        if not any(ln.strip().startswith("TAU_MIX") for ln in kept):
            kept.insert(0, "TAU_MIX = 86400")
        for sp in obs:
            kept.append(f"{sp} = file:constraints/obs_{sp}.txt")
        open(src_bc, "w").write("\n".join(kept) + "\n")

    rt = os.path.join(os.path.dirname(os.path.abspath(dst)), os.path.basename(dst) + "_rt.dat")
    cmd = [PY, CFGBUILD, dst, rt, "--mech", MECH]
    if days is not None:
        cmd += ["--days", f"{days:g}"]
    if dtout is not None:
        cmd += ["--dtout", f"{dtout:g}"]
    run(cmd)
    return rt


# ══════════════════════ 溯源 ══════════════════════
def manifest(exdir, name, params, products, extra=None):
    """写 manifest.json（溯源用：实验名/时间戳/根路径/机制/参数/产物清单）"""
    m = {
        "experiment": name,
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "root": ROOT,
        "runtime_dat": RT,
        "mechanism": MECH,
        "params": params,
        "products": [os.path.relpath(p, ROOT) if os.path.isabs(p) else p
                     for p in products if p],
    }
    if extra:
        m.update(extra)
    p = os.path.join(exdir, "manifest.json")
    with open(p, "w") as f:
        json.dump(m, f, ensure_ascii=False, indent=2)
    return p


def summary_table(title, header, rows, width=None):
    """终端打印简表（不落盘；落盘请用 write_csv）"""
    log("")
    log(f"═══ {title} ═══")
    if not rows:
        log("  (无数据)")
        return
    cells = [[str(h) for h in header]] + [[_fmt(x) for x in r] for r in rows]
    ncol = len(header)
    w = [max(len(c[i]) for c in cells) for i in range(ncol)]
    if width:
        w = [min(a, width) for a in w]
    for j, c in enumerate(cells):
        log("  " + "  ".join(str(c[i])[:w[i]].rjust(w[i]) for i in range(ncol)))
        if j == 0:
            log("  " + "  ".join("-" * w[i] for i in range(ncol)))


def _fmt(x):
    if isinstance(x, float):
        if np.isnan(x):
            return "nan"
        if abs(x) >= 1e4 or (x != 0 and abs(x) < 1e-3):
            return f"{x:.3e}"
        return f"{x:.3f}"
    return str(x)
