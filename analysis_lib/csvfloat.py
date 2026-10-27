#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
csvfloat.py —— 兼容转发层（**实现已单源化到 csvio.py**，2026-10-27）

保留原因: 历史代码仍以 `from csvfloat import fix_fortran_exp, fnum` 引用
（tests/functional_test.py、experiments/viz/fig_ridge5.py）。
新代码请直接 `from csvio import ...`。

背景（原始说明）: gfortran 的 ES 编辑符在 |指数| >= 100 时因字段宽度不足
**静默省略 E**（如 3.07932-174），ES15.5E3 可强制保留。csvio.fix_fortran_exp 统一兜底。
"""
try:                          # 包式导入 (from analysis_lib.csvfloat import ...)
    from .csvio import PPB, fix_fortran_exp, fnum, load_csv_rows
except ImportError:           # 扁平导入 (analysis_lib 在 sys.path 时)
    from csvio import PPB, fix_fortran_exp, fnum, load_csv_rows

__all__ = ["PPB", "fix_fortran_exp", "fnum", "load_csv_rows"]
