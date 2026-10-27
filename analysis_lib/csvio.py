#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
csvio.py —— CSV 契约读写的**唯一实现**（analysis_lib 层；纯 I/O，不启动仿真、不写实验目录）

职责边界:
  · 本模块只做"CSV 文本 ↔ 数值"的两件事: 读（含 gfortran 三位指数省略 E 的兜底解析）
    与写（可带 # 注释头）。
  · 不 import 任何仿真/编排模块（experiments/run/_runner.py 等）——依赖方向恒为单向。

与旧实现的关系（2026-10-27 单源化）:
  · 原 experiments/run/_lib.py 的 _fnum/read_table/col/last24_mean/write_csv → 本模块
  · 原 analysis_lib/csvfloat.py 的 fix_fortran_exp/fnum/load_csv_rows → 本模块
    （csvfloat.py 降为兼容转发层）
  · 原 analysis_lib/diagnostics.py 的本地 _fix_fortran_exp → 本模块
    （diagnostics.py 改为 import，函数名保留为别名）

用法（analysis_lib 已在 sys.path 时）:
    from csvio import read_table, col, last24_mean, write_csv, PPB
"""
import csv
import os
import re

import numpy as np

# 单位换算常数（1013.25 hPa / 298.15 K；全仓唯一定义点）
PPB = 2.46e10

_PAT_NEG = re.compile(r"^(-?[0-9.]+)-([0-9]{3})$")
_PAT_POS = re.compile(r"^(-?[0-9.]+)\+([0-9]{3})$")


def fix_fortran_exp(tok):
    """'3.07932-174' → '3.07932E-174'; '-9.75647-278' → '-9.75647E-278'; 其余原样"""
    if not isinstance(tok, str):
        return tok
    s = tok.strip()
    if "E" in s or "e" in s:
        return s
    m = _PAT_NEG.match(s)
    if m:
        return f"{m.group(1)}E-{m.group(2)}"
    m = _PAT_POS.match(s)
    if m:
        return f"{m.group(1)}E+{m.group(2)}"
    return s


def fnum(tok):
    """修复后转 float；空串/不可解析 → nan（不抛异常，便于整表解析）"""
    s = str(tok).strip()
    if not s:
        return np.nan
    try:
        return float(s)
    except ValueError:
        pass
    try:
        return float(fix_fortran_exp(s))
    except ValueError:
        return np.nan


def load_csv_rows(path):
    """读 CSV, 返回 (header, rows_of_float) —— 全部数值列已修复格式"""
    with open(path) as f:
        rd = list(csv.reader(f))
    if not rd:
        return [], []
    hdr = [h.strip() for h in rd[0]]
    out = []
    for r in rd[1:]:
        if not r or not any(x.strip() for x in r):
            continue
        out.append([fnum(x) for x in r])
    return hdr, out


def read_table(path):
    """读 CSV → (header[list[str]], data[list[list[float]]])，首列常为 h 或坐标"""
    rows = [r for r in csv.reader(open(path)) if r and any(x.strip() for x in r)]
    if not rows:
        return [], []
    hdr = [h.strip() for h in rows[0]]
    data = [[fnum(x) for x in r] for r in rows[1:]]
    return hdr, data


def col(hdr, data, name):
    """按列名取 numpy 数组"""
    if name not in hdr:
        raise KeyError(f"列 '{name}' 不存在; 可用列: {hdr[:8]}...")
    i = hdr.index(name)
    return np.array([r[i] for r in data], dtype=float)


def last24_mean(path, species, ppb=True):
    """conc CSV 末 24 行(小时) 的物种均值；sp 缺列返回 nan"""
    hdr, data = read_table(path)
    if species not in hdr or not data:
        return np.nan
    v = col(hdr, data, species)[-24:]
    return float(np.nanmean(v) / (PPB if ppb else 1.0))


def write_csv(path, header, rows, note=None):
    """写 CSV（可带 # 注释头）；自动建父目录"""
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w", newline="") as f:
        if note:
            for ln in note.splitlines():
                f.write(f"# {ln}\n")
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(r)
    return path
