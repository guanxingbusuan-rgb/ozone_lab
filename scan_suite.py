#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
scan_suite.py —— ozone_lab 排放/气象响应面套件 (P1: batch 引擎)
  引擎: fac_batch_omp (111sp mechanism_3 V2.4, 线程级并行, 16→24 列指示因子+IPR)
  ① 运行 9×9 (f_NOx×f_VOC) 响应面 → output/batch_v3_indicators.csv (+_ipr.csv)
  ② 体制诊断表: 屋脊位置 / FNR 过渡带 / L_NOx-Q 穿越 / OPE 谷
  ③ 图: make figures (fig_indicators_ridge / fig_response_v3)
旧 15sp fac_run 逐案例引擎已退役 (响应面能力由 batch 全面覆盖;
fac_run/kpp_run 已退役, 响应面/冒烟由 batch+fac_v3 全面接管)。
"""
import os
import subprocess
import sys
import numpy as np

# R13: Fortran 三位指数(3.07932-174)兜底 — 共享实现
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) if _os.path.basename(_os.path.dirname(_os.path.abspath(__file__))) == "experiments" else _os.path.dirname(_os.path.abspath(__file__)), "diagnostics"))
from csvfloat import fix_fortran_exp as _fixexp, fnum as _fnum

HERE = os.path.dirname(os.path.abspath(__file__))
BATCH = f"{HERE}/facsimile_impl/fac_batch_omp"
OUT = f"{HERE}/output/batch_v3_indicators.csv"
RT = f"{HERE}/runtime_v3.dat"
PPB = 2.46e10


def main():
    ng = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    os.makedirs(f"{HERE}/output", exist_ok=True)
    r = subprocess.run([BATCH, str(ng), OUT, RT], cwd=f"{HERE}/facsimile_impl")
    if r.returncode != 0:
        sys.exit("batch 失败 (含部分格点求解失败 — 见 stdout 的 [FATAL])")
    lines = [l for l in open(OUT) if l.strip() and not l.startswith("fnox")]
    recs = [[_fnum(x) for x in l.split(",")] for l in lines]
    arr = np.array(recs)
    fnox_u = sorted(set(arr[:, 0]))
    fvs = sorted(set(arr[:, 1]))
    fv_target = min(fvs, key=lambda v: abs(v - 1.0))
    col = sorted([r for r in arr if abs(r[1] - fv_target) < 1e-6], key=lambda r: r[0])
    o3 = [r[4] for r in col]
    iridge = int(np.argmax(o3))
    fnr = [r[8] for r in col]
    lq = [r[14] for r in col]
    ope = [r[16] for r in col]
    print("\n===== 体制诊断 (f_VOC=1 切面, 111sp mechanism_3 V2.4) =====")
    print(f"  屋脊: f_NOx = {col[iridge][0]:.2f} (O3max {o3[iridge]:.1f} ppb)")
    it = next((i for i in range(len(fnr) - 1) if fnr[i] > 2 >= fnr[i + 1]), None)
    if it is not None:
        print(f"  FNR 穿 2.0: f_NOx {col[it][0]:.2f} → {col[it + 1][0]:.2f} (Sillman 过渡带)")
    il = next((i for i in range(len(lq) - 1) if lq[i] < 0.5 <= lq[i + 1]), None)
    if il is not None:
        print(f"  L_NOx/Q 穿 0.5 (Kleinman): f_NOx {col[il][0]:.2f} → {col[il + 1][0]:.2f}")
    iop = int(np.argmin(ope))
    print(f"  OPE 谷: f_NOx = {col[iop][0]:.2f} ({ope[iop]:.1f})")
    print(f"  数据: {OUT} (+_ipr.csv) — 图: make figures")


if __name__ == "__main__":
    main()
