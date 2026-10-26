#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
diagnostics.py —— ozone_lab 诊断库: 指示因子 + IRR/IPR + family + 化学动力账本
消费三实现的统一输出: <prefix>_conc.csv + <prefix>_irr.csv
反应序号 = ozone_kpp.eqn 的 K1..K25 (与 .fac 行序一致)
"""
import numpy as np
import csv

PPB = 2.46e10
# 反应索引 (1-based K1..K25)

def _fix_fortran_exp(tok):
    """gfortran 3 位指数省略 E 的遗留格式 (8.27280-108) → 标准科学计数"""
    import re as _re
    if "E" not in tok and "e" not in tok:
        m = _re.search(r"(\d)-(\d{3})$", tok)
        if m:
            return tok[:m.start(1) + 1] + "E" + tok[m.start(1) + 1:]
    return tok


def load(prefix):
    conc = list(csv.reader(open(f"{prefix}_conc.csv")))[1:]
    irr = list(csv.reader(open(f"{prefix}_irr.csv")))[1:]
    h = np.array([float(r[0]) for r in conc])
    C = np.array([[float(_fix_fortran_exp(x)) for x in r[1:]] for r in conc])   # (T, ns)
    irr_h = np.array([float(r[0]) for r in irr])
    I = np.array([[float(_fix_fortran_exp(x)) for x in r[1:]] for r in irr])    # (T, nr)
    return h, C, irr_h, I

# ============ ① 光化学指示因子 (Sillman 1995 框架的教学版) ============

# ============ ② IRR: 积分反应速率 (CMAQ IRR 教学版) ============


# ============ ③ IPR: 过程贡献分解 (路线无关, 六本账梯形积分) ============
IPR_PROC = ["P", "L", "E", "D", "B", "X"]   # 机制产/耗/排放/干湿沉/BC+稀释/BLH膨胀

def ipr_from_fac(prefix):
    """读 <prefix>_ipr.csv → 逐物种六账总积分 [molec/cm3]

    契约: fac_v3/batch 的 _ipr.csv 每行**已是该小时窗口的精确积分**
          (源码: ipr_int - ipr_w = 本小时窗积分)。因此总积分 = **逐行求和**,
          绝不能再做 np.trapz(那样是二次积分 —— 曾有的 bug)。

    返回 (species, {sp: {proc: value}})，value = Σ_窗口 该过程积分

    注: 用总积分做"闭合校验"会被大项对消放大误差(P/L 各~1e13 抵消到~1e10,
        1.5% 的逐窗偏差可致总和符号翻转)。闭合校验请用**逐窗口比值中位**,
        见 experiments/ipr_closure_check.py 与 run_D_process.py。
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

# ============ ④ family 追踪 (NOx/HOx/ROx/Cl 教学版) ============

# ============ ⑤ 化学动力账本 (Day 15 黑炭式: 每物种来源/去向) ============
RXN_FOR = {
    # 每物种的 IRR 通道 (来源/去向), 由机制表静态生成
    "O3": (["HO2_NO", "CH3O2_NO", "O_O2", "J_O3P"],   # 来源: NO 氧化循环节点
          ["J_O3D", "O3_NO", "O3_NO2"]),              # 去向: 光解+滴定
    "HCHO": (["CH3O_O2", "MOOH_OH", "J_MOOH"], ["J_HCHO", "HCHO_OH"]),
    "HNO3": (["NO2_OH"], ["J_HNO3"]),
    "H2O2": (["HO2_HO2"], ["J_H2O2", "H2O2_OH"]),
}

