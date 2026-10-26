#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ipr_closure_check.py —— fac_v3 IPR 分解完备性验证
  契约: <pref>_ipr.csv 行 h=k 的六账 = 窗口 (k-1, k] 的精确积分 (ONE_STEP 梯形)
  检验: Σ六账[i] ≈ C[i] - C[i-1] (同一窗口) → 比值中位 ≈ 1
用法: python ipr_closure_check.py <prefix>
"""
import sys, csv
import numpy as np
HERE = __file__.rsplit("/experiments/", 1)[0]
sys.path.insert(0, f"{HERE}/diagnostics")
import diagnostics as dg

pref = sys.argv[1] if len(sys.argv) > 1 else "/tmp/v3t2"
rows = list(csv.reader(open(f"{pref}_ipr.csv")))
hdr = [h.strip() for h in rows[0]]
D = np.array([[float(dg._fix_fortran_exp(x)) for x in r] for r in rows[1:] if r])
h, C, _, _ = dg.load(pref)
chdr = [c.strip().split()[0] for c in list(csv.reader(open(f"{pref}_conc.csv")))[0]]
n = min(D.shape[0], C.shape[0])
print("分解完备性 (窗口对位: D[i]=窗(i-1,i] ↔ C[i]-C[i-1]):")
for sp in ["O3", "NO2", "HNO3", "H2O2", "CO", "PAN", "ISOP", "CH3O2", "OH"]:
    if sp not in chdr:
        continue
    ic = chdr.index(sp) - 1                       # dg.load C 剥 h 列
    cols = [hdr.index(f"{sp}_{p}") for p in ["P", "L", "E", "D", "B", "X"]]
    six = D[:, cols].sum(axis=1)
    npair = min(len(six), C.shape[0]) - 1             # ipr 行 i ↔ C[i]-C[i-1]
    six_a = six[1:1+npair]
    dC = C[1:1+npair, ic] - C[:npair, ic]
    big = np.abs(six_a) > 1e6
    if big.sum() < 3:
        print(f"  {sp:6s} 信量不足 (快物种/近零) — 跳过")
        continue
    ratio = np.median(dC[big] / six_a[big])
    flag = "✓" if abs(ratio - 1) < 0.05 else ("~" if abs(ratio - 1) < 0.15 else "✗")
    print(f"  {sp:6s} 中位比值 = {ratio:.4f} {flag}  ({big.sum()} 窗)")
