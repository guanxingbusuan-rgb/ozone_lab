#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ridge5.py —— 五路线光化学指示因子屋脊实验
  路线 × 情景空间:
    [已退役] fac_run  15sp 教学 (2026-10 退役, 代码级守卫 run_fac_run)
    [已退役] kpp_run  15sp EKMA (2026-10 退役, kpp_impl/ 整目录删除)
    mcm_run      29sp MCM  | 排放空间 | CVODE-DENSE
    fac_batch_omp 74sp 工程 V2.2 | 排放空间 | CVODE-DENSE (复用 9×9 网格 fvoc≈1 列)
    fac_v3       = batch 同机制 (单箱; 不重复实验)
  指标: 逐路线 O3 屋脊 fnox* 与各指示因子过渡 fnox (区间中位穿越), 一致性 = |Δfnox|/fnox*
  网格: fnox ∈ logspace(-1, 1, 13), fvoc = 1
用法: python experiments/ridge5.py [--fast]
"""
import os
import subprocess
import sys
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{HERE}/diagnostics")
from indicators import route_indicators   # noqa: E402

# R13: Fortran 三位指数(3.07932-174)兜底 — 共享实现
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))) if _os.path.basename(_os.path.dirname(_os.path.abspath(__file__))) == "experiments" else _os.path.dirname(_os.path.abspath(__file__)), "diagnostics"))
from csvfloat import fix_fortran_exp as _fixexp, fnum as _fnum

PPB = 2.46e10
FN = np.logspace(-1.5, 1.5, 15)      # 0.032 .. 31.6 (mcm 屋脊在 0.1 之下, kpp 在 10 之上)
OUTD = f"{HERE}/output/ridge5"


def run_fac_run(fn):
    raise SystemExit("fac_run 已退役")


def run_mcm(fn):
    pref = f"{OUTD}/mc_{fn:.2f}"
    subprocess.run([f"{HERE}/mcm_city/mcm_run", "5", pref, f"{fn:.3f}", "1.0"],
                   capture_output=True, cwd=f"{HERE}/mcm_city")
    # v2.3.1 起 mcm_run 直接写 out_pref (旧版写死 mcm_city_* 后改名 — 已修)
    return route_indicators(pref, f"{HERE}/mcm_city/input/mechanism_mcm_ch4.fac")


def load_batch():
    """batch 9×9 的 fvoc≈1 列 → {(fnox): indicators}"""
    import csv

    rows = list(csv.reader(open(f"{HERE}/output/batch_v3_indicators.csv")))
    hdr = rows[0]
    ix = {h: i for i, h in enumerate(hdr)}
    out = {}
    for r in rows[1:]:
        f = [_fnum(x) for x in r]
        if abs(f[1] - 1.0125) < 0.01 and abs(f[2] - 1.0) < 1e-6 and abs(f[3]) < 1e-6:
            out[f[0]] = {
                "FNR": f[ix["FNR_mean"]], "H2O2_HNO3": f[ix["H2O2HNO3"]],
                "O3_NOy": f[ix["O3NOy"]], "HCHO_NOy": f[ix["HCHONOy"]],
                "O3_NOz": f[ix["O3NOz"]], "NOy_ppb": f[ix["NOy_ppb"]],
                "L_NOx_Q": f[ix["LnOxQ"]], "OPE": f[ix["OPE"]],
                "O3_ppb": f[ix["O3max"]],
            }
    return out


def ridge_fnox(o3_by_fn):
    """O3 屋脊 = argmax(O3)"""
    fns = sorted(o3_by_fn)
    o3 = [o3_by_fn[f]["O3_ppb"] for f in fns]
    return fns[int(np.argmax(o3))], fns, o3


def indicator_transition(o3_by_fn, key):
    """指示因子过渡 = 穿越 (min+max)/2 中位的 fnox (线性插值)"""
    fns = sorted(o3_by_fn)
    v = np.array([o3_by_fn[f][key] for f in fns])
    lo, hi = v.min(), v.max()
    mid = (lo + hi) / 2
    for i in range(len(v) - 1):
        a, b = v[i], v[i + 1]
        if (a - mid) * (b - mid) <= 0 and a != b:
            frac = (mid - a) / (b - a)
            return fns[i] * (1 - frac) + fns[i + 1] * frac
    return None


def main():
    os.makedirs(OUTD, exist_ok=True)

    routes = {}
    print("运行 mcm_run (29sp MCM 排放空间)...")
    routes["mcm_run 29sp 排放"] = {f: run_mcm(f) for f in FN}
    print("加载 fac_batch_omp (74sp 排放空间, 9×9 fvoc≈1)...")
    routes["batch 74sp 排放"] = load_batch()

    print("\n===== 五路线屋脊实验: O3 屋脊 vs 指示因子过渡 =====")
    hdr_keys = ["FNR", "H2O2_HNO3", "O3_NOy", "O3_NOz", "L_NOx_Q", "OPE"]
    fmin, fmax = FN.min(), FN.max()
    print(f"{'路线':22s} {'O3屋脊fnox':>10s} " + " ".join(f"{k:>11s}" for k in hdr_keys))
    results = {}
    for name, data in routes.items():
        rf, fns, o3 = ridge_fnox(data)
        rf_b = " *" if rf <= fmin * 1.01 or rf >= fmax * 0.99 else ""
        trans = {}
        for k in hdr_keys:
            if k == "L_NOx_Q":
                v = [data[f][k] for f in fns]
                tt = None
                for i in range(len(v) - 1):
                    if (v[i] - 0.5) * (v[i + 1] - 0.5) <= 0 and v[i] != v[i + 1]:
                        frac = (0.5 - v[i]) / (v[i + 1] - v[i])
                        tt = fns[i] * (1 - frac) + fns[i + 1] * frac
                        break
            elif k == "OPE":
                fns_ = sorted(data)
                ope = [data[f][k] for f in fns_]
                tt = fns_[int(np.argmin(ope))]
            else:
                tt = indicator_transition(data, k)
            if tt is not None and (tt <= fmin * 1.01 or tt >= fmax * 0.99):
                tt = None       # 网格边界值 = 伪迹, 不参与一致性
            trans[k] = tt
        results[name] = (rf, trans)
        row = f"{name:22s} {rf:10.2f}{rf_b} "
        row += " ".join(f"{(trans[k] if trans[k] else float('nan')):11.2f}" for k in hdr_keys)
        print(row)

    print("\n===== 一致性 (|过渡fnox − O3屋脊fnox| / 屋脊fnox) =====")
    print(f"{'路线':22s} " + " ".join(f"{k:>11s}" for k in hdr_keys))
    for name, (rf, trans) in results.items():
        row = f"{name:22s} "
        for k in hdr_keys:
            tt = trans[k]
            row += f"{(abs(tt - rf) / rf if tt else float('nan')):11.2f} "
        print(row)

        print("    mcm 29sp 为真实 MCM 速率; batch 74sp 为工程机制 V2.2 (含 isoprene/卤素).")
    print(f"    数据目录: {OUTD}/")


if __name__ == "__main__":
    main()
