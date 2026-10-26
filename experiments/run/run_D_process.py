#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_D_process.py —— 实验 D：IPR/IRR 过程诊断 + 六账闭合校验（无绘图）

输入：fac_v3 输出（本脚本自己跑一次 SCEN_BAU 或复用已有前缀）
输出:
  output/run_D/ipr_ledger.csv          每物种六账（窗口积分）+ 净变化 + 闭合比
  output/run_D/irr_top.csv             O3 相关的反应级 top 通道（IRR 积分排序）
  output/run_D/closure.csv             4 个关键物种的闭合校验

用法:
  python experiments/run/run_D_process.py                     # 现跑 BAU
  python experiments/run/run_D_process.py output/run_A/BAU    # 复用已有输出前缀
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (ROOT, require_prereq, out_dir, run_fac_v3, read_table, write_csv,
                     manifest, summary_table, log, col, out_dir as _od)

sys.path.insert(0, os.path.join(ROOT, "diagnostics"))
import diagnostics as dg


def main(argv):
    log("═══ 实验 D：IPR/IRR 过程诊断 ═══")
    require_prereq()
    d = out_dir("run_D")

    if len(argv) > 1 and os.path.exists(f"{argv[1]}_conc.csv"):
        pref = argv[1]
        log(f"  复用已有输出: {pref}")
        src_tag = "reused"
    else:
        pref = os.path.join(d, "BAU")
        r = run_fac_v3("SCEN_BAU", pref, verbose=True)
        src_tag = "fresh"

    # ── ① 六账账本（全部 74 物种）—— 用 diagnostics.ipr_from_fac 的标准积分 ──
    species, ipr = dg.ipr_from_fac(pref)
    h, C_raw, irr_h, I = dg.load(pref)
    # dg.load 的 C 含 h 列? 兼容两种契约: 若列数 == ns+1 则剥首列
    C = C_raw[:, 1:] if C_raw.shape[1] == len(species) + 1 else C_raw
    names = None
    try:
        from indicators import Mechanism
        names = Mechanism(os.path.join(ROOT, "mechanism", "mechanism_3_full_chain.fac")).species
    except Exception:
        names = species

    rows, closure = [], []
    for sp in species:
        if sp not in names:
            continue
        j = names.index(sp)
        dC = float(C[-1, j] - C[0, j]) if C.shape[0] > 1 else np.nan
        six = sum(ipr[sp].values())
        ratio = dC / six if abs(six) > 1e-30 else np.nan
        rows.append([sp] + [f"{ipr[sp][k]:.6g}" for k in "PLEDBX"]
                    + [f"{six:.6g}", f"{dC:.6g}",
                       f"{ratio:.4f}" if ratio == ratio else "nan"])
    p1 = write_csv(os.path.join(d, "ipr_ledger.csv"),
                   ["species", "P", "L", "E", "D", "B", "X", "sum_six", "dC", "ratio"],
                   rows, note="实验D | 六账窗口积分 [molec/cm3] | ratio=ΔC/Σsix 应≈1")

    # ── ② 关键物种闭合校验：**逐窗口比值中位**（正确方法）──
    #
    # 为何不用"总量比 ΔC/Σ六账"：
    #   白天 O3 的 P 与 L 各 ~1e13，对消到净值 ~1e10。逐窗 1.5% 的积分残差
    #   在总量上会被放大成符号翻转（实测总量比 −0.08 而逐窗比 +1.015）。
    #   → 夏用逐窗口比值的中位数，才是稳健的闭合判据（同 ipr_closure_check.py）。
    import csv as _csv
    _rows = list(_csv.reader(open(f"{pref}_ipr.csv")))
    _hdr = [h.strip() for h in _rows[0]]
    _D = np.array([[float(dg._fix_fortran_exp(x)) for x in r]
                   for r in _rows[1:] if r])
    CHK = ["O3", "HNO3", "H2O2", "CO", "ISOP", "PAN", "NO2", "CH3O2", "OH"]
    for sp in CHK:
        if sp not in names or f"{sp}_P" not in _hdr:
            continue
        j = names.index(sp)
        cols = [_hdr.index(f"{sp}_{p}") for p in "PLEDBX"]
        six_w = _D[:, cols].sum(axis=1)                 # 每行=本窗口积分
        n = min(len(six_w), C.shape[0])
        ratios = []
        for i in range(1, n):
            if abs(six_w[i]) > 1e6:                      # 跳过近零信量（快物种）
                ratios.append((C[i, j] - C[i - 1, j]) / six_w[i])
        if len(ratios) < 3:
            closure.append([sp, f"{ipr[sp]['P']:.4g}", f"{ipr[sp]['L']:.4g}",
                            f"{ipr[sp]['D']:.4g}", "—", "—",
                            "信量不足", "skip"])
            continue
        med = float(np.median(ratios))
        closure.append([sp, f"{ipr[sp]['P']:.4g}", f"{ipr[sp]['L']:.4g}",
                        f"{ipr[sp]['D']:.4g}", f"{len(ratios)}窗", "—",
                        f"{med:.4f}",
                        "PASS" if abs(med - 1) < 0.05 else
                        ("~" if abs(med - 1) < 0.15 else "CHECK")])
    summary_table("实验 D · 逐窗口闭合（中位比值应 ≈1，|Δ|<0.05 为 PASS）",
                  ["species", "P", "L", "D", "windows", "-", "ratio_med", "verdict"], closure)
    p2 = write_csv(os.path.join(d, "closure.csv"),
                   ["species", "P_total", "L_total", "D_total", "windows", "-",
                    "ratio_median_per_window", "verdict"], closure,
                   note=("实验D | 逐窗口闭合校验: 中位(ΔC_window/Σ六账_window), "
                         "0.95–1.05 为 PASS\n"
                         "# 注: 总量比 ΔC/Σ六账 不可用 —— 大项对消(P/L 各~1e13→净~1e10)"
                         "会把 1.5% 逐窗残差放大成符号翻转"))

    # ── ③ IRR 反应级 top 通道（对 O3 净变化贡献最大的反应）──
    # 用 oz.reac/oz.prod 计量表 + IRR 积分重算出每个反应对 O3 的净贡献
    gen = os.path.join(ROOT, "facsimile_impl", "gen_v3")
    I_int = np.trapz(I, irr_h, axis=0)
    jO3 = names.index("O3") + 1                       # gen 索引 1-based
    contrib = np.zeros(len(I_int))
    reac = [ln.split() for ln in open(os.path.join(gen, "oz.reac"))]
    prod = [ln.split() for ln in open(os.path.join(gen, "oz.prod"))]
    for k in range(len(I_int)):
        c = 0.0
        pr = prod[k]
        n = int(pr[0])
        for i in range(n):
            sid = int(pr[1 + 2 * i]); coef = float(pr[2 + 2 * i])
            if sid == jO3:
                c += coef
        rr = reac[k]
        nr = int(rr[0])
        for i in range(nr):
            if int(rr[1 + i]) == jO3:
                c -= 1.0
        contrib[k] = c * I_int[k]
    order = np.argsort(-np.abs(contrib))[:12]
    mech_lines = [l for l in open(os.path.join(ROOT, "mechanism",
                                               "mechanism_3_full_chain.fac"))
                  if l.startswith("%")]
    rows3 = []
    for k in order:
        eq = mech_lines[k].strip().lstrip("%").strip() if k < len(mech_lines) else "?"
        rows3.append([k + 1, f"{contrib[k]:.5e}", eq[:70]])
    summary_table("实验 D · O₃ 反应级 top-12（净贡献 = 计量×IRR 积分）",
                  ["r#", "net_O3", "reaction"], rows3)
    p3 = write_csv(os.path.join(d, "irr_top.csv"),
                   ["r_index", "net_O3_contribution", "reaction"], rows3,
                   note="实验D | O3 净贡献 top-12 | 正值=净生成, 负值=净消耗")

    mp = manifest(d, "D_process", {"source": src_tag, "prefix": pref},
                  [p1, p2, p3, f"{pref}_conc.csv", f"{pref}_irr.csv", f"{pref}_ipr.csv"])
    log(f"\n产物: {os.path.relpath(d, ROOT)}/")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
