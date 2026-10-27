#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
diagnostics.py —— CSV 契约读入 + IPR 六账积分 (analysis 引擎的数据底座)
消费现役三路线的统一输出: <prefix>_conc.csv + <prefix>_irr.csv + <prefix>_ipr.csv
反应列序 r1..rN = .fac 反应行序 (与 gen_v3 计量表/机制行序一致)
"""
import numpy as np
import csv

# 单源化 (2026-10-27): Fortran 三位指数兜底与 PPB 常数统一在 csvio.py
try:                          # 包式导入
    from .csvio import PPB, fix_fortran_exp as _fix_fortran_exp
except ImportError:           # 扁平导入 (analysis_lib 在 sys.path 时)
    from csvio import PPB, fix_fortran_exp as _fix_fortran_exp


def load(prefix):
    conc = list(csv.reader(open(f"{prefix}_conc.csv")))[1:]
    irr = list(csv.reader(open(f"{prefix}_irr.csv")))[1:]
    h = np.array([float(r[0]) for r in conc])
    C = np.array([[float(_fix_fortran_exp(x)) for x in r[1:]] for r in conc])   # (T, ns)
    irr_h = np.array([float(r[0]) for r in irr])
    I = np.array([[float(_fix_fortran_exp(x)) for x in r[1:]] for r in irr])    # (T, nr)
    return h, C, irr_h, I

# ── IPR: 过程贡献分解 (路线无关, 六本账) ──
IPR_PROC = ["P", "L", "E", "D", "B", "X"]   # 机制产/耗/排放/干湿沉/BC+稀释/BLH膨胀

def ipr_from_fac(prefix):
    """读 <prefix>_ipr.csv → 逐物种六账总积分 [molec/cm3]

    契约: fac_v3/batch 的 _ipr.csv 每行**已是该小时窗口的精确积分**
          (源码: ipr_int - ipr_w = 本小时窗积分)。因此总积分 = **逐行求和**,
          绝不能再做 np.trapz(那样是二次积分 —— 曾有的 bug)。

    返回 (species, {sp: {proc: value}})，value = Σ_窗口 该过程积分

    注: 用总积分做"闭合校验"会被大项对消放大误差(P/L 各~1e13 抵消到~1e10,
        1.5% 的逐窗偏差可致总和符号翻转)。闭合校验请用**逐窗口比值中位**,
        见 experiments/ipr_closure_check.py。
    """
    rows = list(csv.reader(open(f"{prefix}_ipr.csv")))
    hdr = [h.strip().split("_")[0] for h in rows[0]]
    data = np.array([[float(_fix_fortran_exp(x)) for x in r] for r in rows[1:] if r])
    species, out = [], {}
    j = 1
    while j < len(hdr) and hdr[j] != "fnox":
        sp = hdr[j]
        species.append(sp)
        out[sp] = {}
        for pi, proc in enumerate(IPR_PROC):
            out[sp][proc] = float(np.sum(data[:, j + pi]))
        j += 6
    return species, out

