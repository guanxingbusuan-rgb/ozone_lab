#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
csvfloat.py —— Fortran CSV 数值解析共享实现 (R13)

背景: gfortran 的 ES 编辑符在 |指数| >= 100 时因字段宽度不足**静默省略 E**
      （如 3.07932-174），而 ES15.5E3 可强制保留。本模块统一兜底旧/新格式。

用法:
    from csvfloat import fix_fortran_exp, load_csv_rows
"""
import re

_PAT = re.compile(r"^(-?[0-9.]+)-([0-9]{3})$")


def fix_fortran_exp(tok):
    """'3.07932-174' → '3.07932E-174'; '-9.75647-278' → '-9.75647E-278'; 其余原样"""
    if not isinstance(tok, str):
        return tok
    s = tok.strip()
    if "E" in s or "e" in s:
        return s
    m = _PAT.match(s)
    if m:
        return f"{m.group(1)}E-{m.group(2)}"
    # 正指数无 E 的形态 (如 1.00000+200)
    m = re.match(r"^(-?[0-9.]+)\+([0-9]{3})$", s)
    if m:
        return f"{m.group(1)}E+{m.group(2)}"
    return s


def fnum(tok):
    """修复后转 float"""
    return float(fix_fortran_exp(tok))


def load_csv_rows(path):
    """读 CSV, 返回 (header, rows_of_float) —— 全部数值列已修复格式"""
    import csv
    with open(path) as f:
        rd = list(csv.reader(f))
    hdr = [h.strip() for h in rd[0]]
    out = []
    for r in rd[1:]:
        if not r or not any(x.strip() for x in r):
            continue
        out.append([fnum(x) for x in r])
    return hdr, out
