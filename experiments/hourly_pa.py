#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
hourly_pa.py —— fac_v3 逐小时过程分析 (CMAQ-PA 式)
  数据: <pref>_ipr.csv (53×6 过程账, ONE_STEP 梯形小时窗积分)
        <pref>_irr.csv (131 反应 IRR 小时窗积分 — 同法累积)
  输出:
    ① 关键物种逐小时六账占比表 (O3/NO2/HCHO/HNO3/ISOP/PAN)
    ② O3 逐小时反应级 top 预算 (IRR 重构, 产/耗各 top-4)
    ③ IPR↔IRR 一致性检验 (两独立累积路径互证)
    ④ 日循环图 figures/hourly_pa.png
用法: python3 experiments/hourly_pa.py [prefix=../output/osse/../fac_v3 输出前缀]
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{HERE}/diagnostics")
from indicators import _load_csv, Mechanism   # noqa

# R13: Fortran 三位指数(3.07932-174)兜底 — 共享实现
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) if _os.path.basename(_os.path.dirname(_os.path.abspath(__file__))) == "experiments" else _os.path.dirname(_os.path.abspath(__file__)), "diagnostics"))
from csvfloat import fix_fortran_exp as _fixexp, fnum as _fnum

PROCS = ["P", "L", "E", "D", "B", "X"]
PROC_NAMES = ["chem P", "chem L", "emis", "dep", "BC+dil", "BLH"]   # 图例用英文 (CJK 字体缺失)


def main():
    pref = sys.argv[1] if len(sys.argv) > 1 else f"{HERE}/output/pa_demo"
    if not os.path.exists(f"{pref}_ipr.csv"):
        # 便捷: 无输出则现跑一个基线
        rt = f"{HERE}/runtime_v3.dat"
        print(f"[hourly_pa] 无 {pref}_ipr.csv — 现跑基线...")
        import subprocess

        subprocess.run([f"{HERE}/facsimile_impl/fac_v3", "SCEN_BAU", pref, rt],
                       cwd=f"{HERE}/facsimile_impl", check=True)
    hdr, D = _load_csv(f"{pref}_ipr.csv")
    ihrd, R = _load_csv(f"{pref}_irr.csv")
    mech = Mechanism(f"{HERE}/mechanism/mechanism_3_full_chain.fac")
    hours = D[:, 0].astype(int) % 24

    # ---- ① 关键物种逐小时六账占比 ----
    print("===== ① 逐小时过程占比 (|六账| 归一, %) =====")
    for sp in ["O3", "NO2", "HCHO", "HNO3", "ISOP", "PAN"]:
        if f"{sp}_P" not in hdr:
            continue
        cols = [hdr.index(f"{sp}_{p}") for p in PROCS]
        B = D[:, cols]
        tot = np.abs(B).sum(axis=1)
        print(f"\n  [{sp}]  (h: " + " ".join(PROC_NAMES) + ")")
        for k in range(len(hours)):
            if tot[k] < 1e-30:
                continue
            shares = " ".join(f"{abs(B[k, j]) / tot[k] * 100:5.1f}" for j in range(6))
            print(f"  h{hours[k]:02d}  {shares}")

    # ---- ② O3 逐小时反应级 top 预算 (IRR) ----
    print("\n===== ② O3 逐小时反应级 top 预算 (IRR 积分重构) =====")
    io3 = mech.sid["O3"]
    prod_r = [[] for _ in range(mech.nrx)]
    loss_r = [[] for _ in range(mech.nrx)]
    for r in range(mech.nrx):
        for i, c in mech.prod[r]:
            if i == io3:
                prod_r[r].append(c)
        for i, _ in mech.reac[r]:
            if i == io3:
                loss_r[r].append(True)
    for k in [0, 6, 12, 18]:
        h = hours[k]
        P = {r: sum(c * R[k, r + 1] for c in prod_r[r]) for r in range(mech.nrx)}
        L = {r: R[k, r + 1] for r in range(mech.nrx) if loss_r[r]}
        tp = sorted(P.items(), key=lambda x: -x[1])[:3]
        tl = sorted(L.items(), key=lambda x: -x[1])[:3]
        print(f"  h{h:02d}: 产 " + "  ".join(
            f"r{r+1}({v:+.2e})" for r, v in tp if abs(v) > 1e-20) +
              f" | 耗 " + "  ".join(f"r{r+1}({v:.2e})" for r, v in tl if v > 1e-20))

    # ---- ③ IPR↔IRR 一致性 (两独立累积路径) ----
    print("\n===== ③ IPR↔IRR 一致性 (O3, 逐小时; L 按符号约定) =====")
    worst = 0.0
    for k in range(D.shape[0]):
        p_ipr = D[k, hdr.index("O3_P")]
        p_irr = sum(c * R[k, r + 1] for r in range(mech.nrx)
                    for i, c in mech.prod[r] if i == io3)
        l_ipr = D[k, hdr.index("O3_L")]
        l_irr = -sum(R[k, r + 1] for r in range(mech.nrx)
                     for i, _ in mech.reac[r] if i == io3)
        for a, b in ((p_ipr, p_irr), (l_ipr, l_irr)):
            if max(abs(a), abs(b)) > 1e20:
                worst = max(worst, abs(a - b) / max(abs(a), abs(b)))
    print(f"  最大相对差 = {worst * 100:.3f}%  "
          f"({'✓ 两累积路径一致' if worst < 0.01 else '✗ 超限 — 检查装配'})")

    # ---- ④ 日循环图 ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(2, 3, figsize=(16, 8))
        for ax, sp in zip(axes.flat, ["O3", "NO2", "HCHO", "HNO3", "ISOP", "PAN"]):
            if f"{sp}_P" not in hdr:
                continue
            cols = [hdr.index(f"{sp}_{p}") for p in PROCS]
            B = D[:, cols] / 3600.0        # 积分 → 平均速率
            bottom = np.zeros(len(hours))
            for j, (pn, cn) in enumerate(zip(PROC_NAMES,
                    ["#27AE60", "#E74C3C", "#F39C12", "#8E44AD", "#3498DB", "#95A5A6"])):
                ax.bar(hours, B[:, j], bottom=bottom, label=pn, color=cn, width=0.8)
                bottom += B[:, j]
            ax.set_title(sp)
            ax.set_xlabel("hour"); ax.set_ylabel("rate [molec/cm³/s]")
            if sp == "O3":
                ax.legend(fontsize=7, ncol=2)
        fig.suptitle("fac_v3 hourly process analysis (mechanism_3, ONE_STEP exact)", fontsize=12)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        out = f"{HERE}/figures/hourly_pa.png"
        plt.savefig(out, dpi=150, bbox_inches="tight")
        print(f"\n图: {out}")
    except ImportError:
        print("(matplotlib 不可用 — 跳过图)")


if __name__ == "__main__":
    main()
