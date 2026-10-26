#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ipr_irr5.py —— 五路线 IPR/IRR 过程诊断分析
  ① 收支闭合验证: 各路线 IPR 六账 vs 数值 dC/dt (长寿命物种; 快物种频闪已文档化)
  ② 求解器无关性: fac_v3 (fcmix-DENSE 单箱) vs fac_batch_omp (C-API OpenMP)
     同情景 131 反应 IRR 逐项对照 —— 旗舰双胞胎的交叉验收
  ③ 过程分解 (CMAQ-PA 式): O3/HCHO/HNO3 的六账贡献占比 per 路线
用法: python experiments/ipr_irr5.py
"""
import os
import sys
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{HERE}/diagnostics")
import diagnostics as dg
from indicators import _fix_fortran_exp   # noqa
import csv

PPB = 2.46e10
LONG_LIVED = ["O3", "HNO3", "H2O2", "CO"]


def load_ipr(prefix):
    rows = list(csv.reader(open(f"{prefix}_ipr.csv")))
    hdr = [h.strip() for h in rows[0]]
    data = np.array([[float(dg._fix_fortran_exp(x)) for x in r] for r in rows[1:] if r])
    species, acc = [], {}
    j = 1
    while j + 6 <= len(hdr) and hdr[j] != "fnox":
        sp = hdr[j].rsplit("_", 1)[0]          # "NO_P" → "NO"
        cols = [hdr.index(f"{sp}_{p}") for p in ["P", "L", "E", "D", "B", "X"]]
        acc[sp] = data[:, cols]
        species.append(sp)
        j += 6
    return species, acc, data[:, 0]


def closure(prefix, label, sp_list=None):
    """IPR 六账窗口均值 vs (C_end − C_start)/Δt — 长寿命物种比值应≈1"""
    _, acc, _ = load_ipr(prefix)
    h, C, _, _ = dg.load(prefix)
    chdr = [c.strip().split()[0] for c in
            list(csv.reader(open(f"{prefix}_conc.csv")))[0]]
    # 契约注意: dg.load 的 C 已剥 h 列, chdr 含 h → index 需 −1
    print(f"  [{label}]")
    for sp in (sp_list or LONG_LIVED):
        if sp not in acc or sp not in chdr:
            print(f"    {sp:6s} — 无列 (跳过)")
            continue
        ic = chdr.index(sp) - 1
        dC = (C[-1, ic] - C[0, ic]) / (h[-1] - h[0]) / 3600.0   # /s
        six = acc[sp].mean(axis=0).sum()
        ratio = dC / six if abs(six) > 1e-30 else float("nan")
        print(f"    {sp:6s} dC/dt={dC:+.3e}  Σ六账={six:+.3e}  比值={ratio:.3f}")


def irr_cross(pref_a, pref_b, nr, label):
    """131 反应 IRR 逐项对照 (同情景双求解器)"""
    _, _, ha, Ia = dg.load(pref_a)
    _, _, hb, Ib = dg.load(pref_b)
    ra, rb = Ia[-12:].mean(axis=0), Ib[-12:].mean(axis=0)   # 末 12h 均值
    scale = np.maximum(np.abs(ra), np.abs(rb))
    ok = scale > 1e-3          # 忽略 ~0 速率 (夜间光解等)
    rel = np.where(ok, np.abs(ra - rb) / scale, 0.0)
    big = np.where(ok & (np.maximum(np.abs(ra), np.abs(rb)) > 1e4))[0]
    worst = big[np.argsort(-rel[big])][:5] if len(big) else []
    med = np.median(rel[ok]) * 100
    print(f"  [{label}] 有效反应 {ok.sum()}/{nr}, 相对差中位 {med:.2f}%")
    for k in worst:
        print(f"    r{k+1:3d}: fac_v3={ra[k]:+.3e} batch={rb[k]:+.3e} rel={rel[k]*100:.1f}%")


def process_share(prefix, label, sp_list):
    """六账贡献占比 (窗口均值, |归一|) — CMAQ-PA 式过程分解"""
    species, acc, _ = load_ipr(prefix)
    procs = ["机制P", "机制L", "排放", "干湿沉", "BC+稀释", "BLH膨胀"]
    print(f"  [{label}]")
    for sp in sp_list:
        if sp not in acc:
            continue
        m = acc[sp].mean(axis=0)
        tot = np.abs(m).sum()
        if tot < 1e-30:
            continue
        shares = "  ".join(f"{p}:{v / tot * 100:5.1f}%" for p, v in zip(procs, m))
        print(f"    {sp:6s} {shares}")




def closure_snap(prefix, sp, label):
    """逐快照闭合: six(t) vs 前向差分 dC/dt(t) — 按昼夜定位错账"""
    species, acc, t_h = load_ipr(prefix)
    h, C, _, _ = dg.load(prefix)
    chdr = [c.strip().split()[0] for c in
            list(csv.reader(open(f"{prefix}_conc.csv")))[0]]
    if sp not in acc or sp not in chdr:
        print(f"  [{label}] {sp}: 无列")
        return
    ic = chdr.index(sp) - 1   # dg.load C 剥 h 列 → chdr 索引 −1
    six = acc[sp].sum(axis=1)
    dC = (C[1:, ic] - C[:-1, ic]) / 3600.0     # 前向差分 @t_k
    print(f"  [{label}] {sp} 逐快照 (hour, dC/dt, Σ六账, 差):")
    for k in range(len(dC)):
        hh = h[k] % 24
        tag = "昼" if 6 <= hh <= 18 else "夜"
        print(f"    h={int(h[k]):3d}({tag})  dC/dt={dC[k]:+.3e}  six={six[k]:+.3e}  diff={dC[k]-six[k]:+.3e}")

def main():
    print("===== ① IPR 收支闭合 (长寿命物种, 比值≈1) =====")
    closure(f"{HERE}/output/ridge5/mc_1.00", "mcm_run 29sp MCM")
    closure("/tmp/v3ipr", "fac_v3 74sp runtime 通路")
    # batch 用 ipr csv (窗口均值直接给, 与快照均值口径一致)
    print("  [fac_batch_omp 74sp] 窗口均值六账量级抽查:")
    rows = list(csv.reader(open("/tmp/irr5.csv_ipr.csv")))
    hdr = [h.strip() for h in rows[0]]
    for r in rows[1:]:
        f = [x.strip() for x in r]
        if abs(float(f[0]) - 1.0125) < 0.01 and abs(float(f[1]) - 1.0125) < 0.01:
            for sp in ("O3", "HNO3"):
                cols = [hdr.index(f"{sp}_{p}") for p in ["P", "L", "E", "D", "B", "X"]]
                six = sum(float(f[c]) for c in cols)
                print(f"    {sp:6s} Σ六账均值={six:+.3e}/s")
            break

    print("\n===== ①b 逐快照闭合 (错账定位) =====")
    closure_snap(f"{HERE}/output/ridge5/mc_1.00", "O3", "mcm O3 (178× 疑点)", )

    print("\n===== ② 求解器无关性: fac_v3 vs batch 131 反应 IRR =====")
    print("  [跨求解器对照已随 15sp 层退役; batch(74sp C-API) 自身 IRR↔IPR 一致性见 hourly_pa]")

    print("\n===== ③ 过程分解 (CMAQ-PA 式, 窗口均值占比) =====")
    process_share("/tmp/v3ipr", "fac_v3 74sp", ["O3", "NO2", "HCHO", "HNO3", "ISOP", "PAN"])
    process_share(f"{HERE}/output/ridge5/mc_1.00", "mcm 29sp", ["O3", "NO2", "HCHO", "HNO3"])


if __name__ == "__main__":
    main()
