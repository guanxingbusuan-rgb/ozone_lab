#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
sat2surf.py —— 卫星柱浓度 → 地表等效 ppb 约束序列 (0-D 箱模型的离线转换器)

背景: ozone_lab 是 0-D 箱体, 无气柱概念 —— 引擎不直读卫星柱浓度。
      卫星数据的正确用法 = 离线转换为表面等效体积混合比, 再走 bc.config 张弛约束。

转换公式 (直观法为主准则, 混合比法仅作参考核对):
  ② 直观法 ★主准则★:
       表面等效浓度 [molec/cm3] = 柱浓度 [molec/cm2] / BLH[cm]
       表面等效混合比 [ppb]     = 上式 / PPB
     理由: 只依赖 PPB (1 ppb 数密度, 由 P/kT 直接得出) + BLH 实测,
           物理上就是"柱浓度摊到混合层厚度 = 层平均浓度"的定义式, 假设最少。
  ① 混合比法 (参考):
       混合比 = 柱/空气柱 × 1e9 / (BLH/等效柱高)
     注意: 2.15e25 与 8.74 km 实为同一"标准大气柱"假设的两种表述,
           比直观法多一层假设, 故仅作交叉参考。
  两法分歧时: **输出警告, 不中断**; 一律以直观法结果为准写出文件。
  ⚠ 历史 bug: 旧模板写 柱/(空气柱×混合层分数) 且标 [ppb] —— 那是无量纲
    摩尔混合比, 漏 ×1e9, 差 9 个量级。本工具以量级自检 + 参考核对防复发。

输入格式 (注释头 + 三列数据; --selftest 为三条已知值的换算核对):
  # satellite: TROPOMI/S5P  product: OFFL  L2  v: 2.4.0
  # species: NO2_tropo
  # qa_filter: qa_value > 0.75
  # t[YYYY-MM-DDTHH:MM]  column[molec/cm2]  qa
  2024-06-01T13:30   4.2e15   0.92
  ...

输出格式 (直接可登记到 bc.config):
  # unit: ppb
  # source: sat2surf.py ...
  # t[s]  NO2
        0    1.71
    86400   15.45

用法:
  python inputproc/sat2surf.py <raw.raw> <out.txt> [选项]
选项:
  --blh M            混合层高度常数 [m] (缺省 1000)
  --blh-file PATH    BLH 序列 (两列 t[s] BLH[m]), 在过境时刻线性插值 (优先于 --blh)
  --qa-min 0.75      QA 阈值, 低于则剔除
  --species NAME     覆盖头部 species (无白名单; 凡机制 VARIABLE 物种皆可,
                     常用: NO2_tropo|HCHO_tropo|SO2_tropo|CO|CH4|O3)
  --t0 YYYY-MM-DD    模拟起点日期 (缺省 = 首条观测当日 00:00)
  --max-points 128   Fortran 端上限; 超限自动等间隔抽稀 (会打印警告)
  --selftest         运行内置校验后退出

退出码: 0=成功 (含警告情形) 1=输入错误 (无法产出结果)
"""
import sys, re, argparse
from datetime import datetime, timedelta

# R13: PPB 常量单源 (与 cfgbuild/unit_factor 一致)
# 稳健导入: 直接用脚本路径时 sys.path[0] 已是 inputproc/; 但被 importlib/其它 cwd 加载时
# 需要显式把本文件所在目录加入 sys.path, 否则 ModuleNotFoundError。
import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)
from cfgbuild import PPB

AIR_COL = 2.15e25      # 标准大气空气柱 [molec/cm2]
# 空气柱对应厚度: 由 AIR_COL/PPB 解析派生 (勿手写 8.74e5 —— 四舍五入会与直观式
# 产生 1.9e-5 相对差, 破环两法互证). 值 = 2.15e25/(PPB×1e9) = 8.73984e5 cm
AIR_SCALE = AIR_COL / (PPB * 1e9)   # [cm] ≈ 8.7398 km


def parse_dt(s):
    for fmt in ("%Y-%m-%dT%H:%M", "%Y-%m-%d %H:%M", "%Y%m%d%H", "%Y-%m-%d"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            continue
    raise ValueError(f"无法解析时间 '{s}' (支持 YYYY-MM-DDTHH:MM / YYYYMMDDHH / YYYY-MM-DD)")


def interp_blh(t_sec, ts, vs):
    """BLH 序列线性插值, 端点钳制"""
    if not ts:
        return None
    if t_sec <= ts[0]:
        return vs[0]
    if t_sec >= ts[-1]:
        return vs[-1]
    for i in range(len(ts) - 1):
        if t_sec <= ts[i + 1]:
            f = (t_sec - ts[i]) / max(ts[i + 1] - ts[i], 1e-9)
            return vs[i] + f * (vs[i + 1] - vs[i])
    return vs[-1]


def convert(col, blh_m, air_col=AIR_COL, ppb=PPB, scale_cm=None):
    """返回 (ppb_直观法_主准则, ppb_混合比法_参考)

    直观法 (主准则): col / BLH / ppb
        只依赖 PPB + 实测 BLH —— "柱摊到混合层厚度"的定义式, 假设最少。

    混合比法 (仅参考): col/air_col × 1e9 / (BLH/scale_cm)
        scale_cm 缺省由 air_col 派生 (⇒ 与直观法代数等价);
        若显式传入不同 scale_cm (如假设柱混合到 2 km), 两法将产生真实分歧 ——
        此时工具仅告警, 输出一律采用直观法。
    """
    blh_cm = blh_m * 100.0
    surf = col / blh_cm                       # molec/cm3
    ppb_b = surf / ppb                        # ★ 主准则: 直观法
    scale = scale_cm if scale_cm else air_col / (ppb * 1e9)
    frac = blh_cm / scale                     # 混合层分数
    ppb_a = col / air_col * 1e9 / frac         # 参考: 混合比法
    return ppb_b, ppb_a


def selftest():
    """内置校验: 三条已知值"""
    ok = True
    cases = [
        (3.8e16, 1000.0, 15.45),   # TROPOMI 城市 NO2 柱 + 1 km BLH
        (4.2e15, 1000.0, 1.71),
        (1.0e16, 500.0, 8.13),
    ]
    print("== sat2surf 自检 ==")
    for col, blh, expect in cases:
        b, a = convert(col, blh)
        rel = abs(b - a) / max(abs(b), 1e-30)
        # 主准则判定: 直观法与期望值一致 (两法相对差仅提示, 不作为失败条件)
        good = abs(b - expect) / expect < 0.01
        ok &= good
        flag = "✓" if good else "✗"
        note = "" if rel < 1e-9 else f"  [警告] 两法分歧 {rel:.1e}, 以直观法为准"
        print(f"  柱={col:.2e} BLH={blh:g}m → 直观 {b:.2f} ppb (期望 {expect}) {flag}"
              f" | 参考: 混合比 {a:.2f} ppb (差 {rel:.1e}){note}")
    # 反证: 漏 ×1e9 的旧式
    wrong = 3.8e16 / (AIR_COL * (1e5 / AIR_SCALE))
    print(f"  旧式(漏 1e9) = {wrong:.3e} → 无量纲混合比, 非 ppb; ×1e9 = {wrong*1e9:.2f} ppb")
    print("结果:", "全部通过" if ok else "存在失败")
    return 0 if ok else 2


def main():
    ap = argparse.ArgumentParser(add_help=True)
    ap.add_argument("raw", nargs="?")
    ap.add_argument("out", nargs="?")
    ap.add_argument("--blh", type=float, default=1000.0)
    ap.add_argument("--blh-file")
    ap.add_argument("--qa-min", type=float, default=0.75)
    ap.add_argument("--species")
    ap.add_argument("--t0")
    ap.add_argument("--max-points", type=int, default=128)
    ap.add_argument("--air-col", type=float, default=AIR_COL,
                    help="参考法用的空气柱 [molec/cm2], 缺省 2.15e25 (不影响主准则)")
    ap.add_argument("--scale-height-km", type=float, default=None,
                    help="参考法独立假设的等效柱高 [km]; 缺省由 air-col 派生(与主准则等价)。"
                         "显式给出不同值可产生真实分歧, 工具仅告警并以直观法输出")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if not a.raw or not a.out:
        ap.print_help()
        return 1

    # ---- 读原始文件 ----
    species, rows, qa_thr = None, [], a.qa_min
    with open(a.raw, encoding="utf-8") as f:
        for ln in f:
            s = ln.strip()
            if not s:
                continue
            if s.startswith("#"):
                m = re.search(r"species:\s*([A-Za-z0-9_]+)", s)
                if m:
                    species = m.group(1)
                m = re.search(r"qa_value\s*>\s*([0-9.]+)", s)
                if m:
                    qa_thr = float(m.group(1))
                continue
            p = s.replace(",", " ").split()
            if len(p) < 2:
                continue
            try:
                dt = parse_dt(p[0])
                col = float(p[1])
                qa = float(p[2]) if len(p) > 2 else None
            except ValueError as e:
                print(f"[跳过] 无法解析行: {s} ({e})")
                continue
            rows.append((dt, col, qa))

    if not rows:
        print(f"[错误] {a.raw}: 无有效数据行", file=sys.stderr)
        return 1
    species = (a.species or species or "NO2").replace("_tropo", "")
    rows.sort(key=lambda r: r[0])

    # ---- QA 过滤 ----
    kept = [r for r in rows if (r[2] is None or r[2] >= qa_thr)]
    n_drop = len(rows) - len(kept)
    if n_drop:
        print(f"[QA] 阈值 {qa_thr}: 剔除 {n_drop} 条低质量过境 (保留 {len(kept)})")
    if not kept:
        print(f"[错误] QA 过滤后无数据", file=sys.stderr)
        return 1

    # ---- BLH 序列 ----
    blh_ts, blh_vs = [], []
    if a.blh_file:
        for ln in open(a.blh_file):
            s = ln.split("#")[0].strip()
            if not s:
                continue
            p = s.split()
            blh_ts.append(float(p[0]))
            blh_vs.append(float(p[1]))
        if any(v <= 0 for v in blh_vs):
            print(f"[错误] {a.blh_file}: BLH 含非正值", file=sys.stderr)
            return 1
        print(f"[BLH] 序列 {a.blh_file}: {len(blh_ts)} 点, 范围 {min(blh_vs):g}–{max(blh_vs):g} m")

    # ---- 时间轴 ----
    t0 = parse_dt(a.t0) if a.t0 else kept[0][0].replace(hour=0, minute=0, second=0)
    times = [(r[0] - t0).total_seconds() for r in kept]
    if any(times[i] <= times[i - 1] for i in range(1, len(times))):
        print("[错误] 时间非严格递增 (检查重复过境)", file=sys.stderr)
        return 1

    # ---- 转换: 直观法为主准则; 混合比法仅参考核对 (分歧只警告, 不中断) ----
    out_rows, max_rel = [], 0.0
    for (dt, col, qa), t in zip(kept, times):
        blh = interp_blh(t, blh_ts, blh_vs) if blh_ts else a.blh
        b, aa = convert(col, blh, air_col=a.air_col,
                        scale_cm=(a.scale_height_km * 1e5) if a.scale_height_km else None)
        rel = abs(b - aa) / max(abs(aa), 1e-30)
        max_rel = max(max_rel, rel)
        out_rows.append((t, b, blh, col))     # 一律采用直观法 (主准则)
    if max_rel > 1e-9:
        print(f"[警告] 两法分歧: 最大相对差 {max_rel:.2e} > 1e-9 "
              f"— 以直观法 (柱/BLH/PPB, 假设最少) 为准继续; "
              f"混合比法仅供参考, 不中断", file=sys.stderr)

    # ---- 量级自检 ----
    vals = [r[1] for r in out_rows]
    lo, hi = min(vals), max(vals)
    if hi > 500 or lo < 0.001:
        print(f"[警告] 输出范围 {lo:.3g}–{hi:.3g} ppb 异常 (典型: NO2 5–20, HCHO 1–10, SO2 0.1–5); "
              f"检查 BLH 假设与柱单位", file=sys.stderr)

    # ---- 点数限制 (Fortran 端 ≤128) ----
    stride = 1
    if len(out_rows) > a.max_points:
        stride = (len(out_rows) + a.max_points - 1) // a.max_points
        out_rows = out_rows[::stride]
        print(f"[警告] {len(kept)} 点 > 上限 {a.max_points}, 等间隔抽稀 ×{stride} → {len(out_rows)} 点")

    # ---- 写出 ----
    lines = [
        "# unit: ppb",
        f"# source: sat2surf.py {a.raw} | species={species} | "
        f"QA>={qa_thr} | BLH={'file:' + a.blh_file if a.blh_file else f'{a.blh:g}m'}",
        f"# method: 柱/BLH/PPB (直观法=主准则); 混合比法参考"
        + (f" (air_col={a.air_col:.3g}, scale={a.scale_height_km:g}km)"
           if a.scale_height_km else "")
        + f" max_rel={max_rel:.1e}"
        + ("" if max_rel <= 1e-9 else "  [警告: 两法分歧, 已按直观法输出]"),
        f"# t[s]  {species}",
    ]
    for t, ppb, blh, col in out_rows:
        lines.append(f"{int(t):8d} {ppb:9.3f}")
    open(a.out, "w").write("\n".join(lines) + "\n")

    print(f"\n== 转换结果 → {a.out} ==")
    print(f"  {'过境时刻':20s} {'柱[molec/cm2]':>14s} {'BLH[m]':>8s} {'表面[ppb]':>10s}")
    for (dt, col, qa), t, (tt, ppb, blh, cc) in zip(kept, times, out_rows):
        print(f"  {dt:%Y-%m-%dT%H:%M}   {cc:14.3e} {blh:8.0f} {ppb:10.2f}")
    print(f"\n  范围 {lo:.2f}–{hi:.2f} ppb | 点数 {len(out_rows)} | 两法最大相对差 {max_rel:.1e}")
    print(f"  下一步: bc.config 加行  {species} = file:{a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
