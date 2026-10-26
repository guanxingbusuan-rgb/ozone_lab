#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
cfgbuild.py —— ozone_lab 配置编译器
  input/ (12 类人类可读配置) → runtime.dat (Fortran 顺序表控平面文件)
  编译期完成: 物种存在性校验 / 单位换算(ppb→molec) / 序列单调性校验 /
  情景表达式静态解析 / 光解 l/m/n 表映射 / 排放日型展开
用法: python cfgbuild.py [input_dir] [out_runtime.dat]
"""
import os, re, sys
from pathlib import Path

PPB = 2.46e10

def unit_factor(u):
    """浓度单位 → molec/cm³ 换算因子 (cm3|ppb|ppm); 单一事实源"""
    f = {"cm3": 1.0, "ppb": PPB, "ppm": 1e3 * PPB}.get(str(u).lower())
    if f is None:
        raise ValueError(f"未知单位 '{u}' (cm3|ppb|ppm)")
    return f

def strip_comment(line):
    """去除 # 注释, 返回净内容"""
    return line.split("#")[0].strip()

def read_kv(path):
    """读 'KEY = value' 型配置"""
    out = {}
    for ln in open(path):
        s = strip_comment(ln)
        if not s or "=" not in s: continue
        k, v = [x.strip() for x in s.split("=", 1)]
        out[k] = v
    return out

def read_series(path):
    """读 't v' 两列时间序列, 校验单调
    单位指令: 文件头 '# unit: cm3|ppb|ppm' (缺省 cm3=molec/cm3) — 模板约定"""
    t, v = [], []
    unit = "cm3"
    for ln in open(path):
        s = strip_comment(ln)
        if not s:
            c = ln.strip()
            if c.startswith("#") and "unit:" in c:
                u = c.split("unit:", 1)[1].split()[0].lower()   # 首 token (容忍行内注释)
                if u not in ("cm3", "ppb", "ppm"):
                    raise ValueError(f"{path}: 未知单位 '{u}' (cm3|ppb|ppm)")
                unit = u
            continue
        p = s.split()
        t.append(float(p[0])); v.append(float(p[1]))
    f = unit_factor(unit)
    v = [x * f for x in v]
    for i in range(1, len(t)):
        if t[i] <= t[i-1]:
            raise ValueError(f"{path}: 时间非单调 @ t={t[i]}")
    if len(t) > 128:
        raise ValueError(f"{path}: 序列 {len(t)} 点 > Fortran 端上限 128 "
                         f"(runtime_dat_mod NSER_MAX=8, 每序列 128 点)")
    return t, v

def read_ic(path):
    out = {}
    for ln in open(path):
        s = strip_comment(ln)
        if not s: continue
        p = s.split()
        val = float(p[1])
        unit = p[2].lower() if len(p) > 2 else "cm3"
        out[p[0]] = val * unit_factor(unit)
    return out

def parse_fac_species(fac_path):
    """从 .fac 的 VARIABLE 段读物种列表 (支持多行, 以 ; 结束)"""
    species = []
    in_vars = False
    for ln in open(fac_path):
        s = ln.strip()
        if s == "VARIABLE":
            in_vars = True; continue
        if in_vars:
            ended = s.endswith(";")
            species += s.rstrip(";").split()
            if ended:
                break
    return species

def parse_photolysis_table(path):
    """读 MCM photolysis-rates 表: j l m n name tau"""
    entries = {}
    for ln in open(path):
        s = strip_comment(ln)
        if not s or not s[0].isdigit(): continue
        p = s.split()
        j = int(p[0])
        l_ = float(p[1].replace('D','E'))
        m_ = float(p[2].replace('D','E'))
        n_ = float(p[3].replace('D','E'))
        entries[j] = (l_, m_, n_)
    return entries

def compile_cfg(input_dir="input_v3", out_path="runtime_v3.dat", mech_override=None):
    errs = []
    inp = Path(input_dir)

    # ============ ① 机制文件 ============
    # 优先找 mechanism_*.fac (用户可放多个, 按字母序取第一个)
    facs = [Path(mech_override)] if mech_override else sorted(inp.glob("mechanism_*.fac"))
    facs = [f for f in facs if f.exists()]
    if not facs:
        errs.append("input/ 下无 mechanism_*.fac")
        fac = None; species = []
    else:
        fac = facs[0]
        species = parse_fac_species(fac)
    sid = {s: i for i, s in enumerate(species)}
    ns = len(species)

    # ============ ② 光解参数表 ============
    phot_files = sorted(inp.glob("photolysis-rates_*"))
    phot_tab = {}
    if phot_files:
        phot_tab = parse_photolysis_table(phot_files[0])

    # ============ ③ 光解配置 ============
    ph_cfg = read_kv(inp / "photolysis.config")
    jmode = ph_cfg.get("MODE", "MCM")
    jlat = float(ph_cfg.get("LAT", "45.0"))
    jdoy = float(ph_cfg.get("DOY", "172.0"))
    # JMAP 显式语法: "4=4,13=D4.0e-5" (facJ=slot 或 facJ=D<正午直接值>);
    # 缺省 = 恒等映射; 兼容旧 MAPPING 位置语法
    jmap_raw = ph_cfg.get("JMAP", "")
    if jmap_raw:
        mapping = []           # 延迟: 长度由机制 J 索引最大值定
        jmap_pairs = {}
        for tok in jmap_raw.split(","):
            k, v = tok.split("=")
            jmap_pairs[int(k)] = v
    else:
        jmap_pairs = {}
        mapping = [int(x) for x in ph_cfg.get("MAPPING", "").split(",")] if \
            ph_cfg.get("MAPPING") else []

    # ============ ④ IC ============
    ic = read_ic(inp / "initialConcentrations.config")
    y0 = [0.0]*ns
    for k, v in ic.items():
        if k not in sid:
            errs.append(f"IC 未知物种 '{k}' (不在机制 VARIABLE 段)")
        else:
            y0[sid[k]] = v

    # ============ ⑤ 气象 ============
    METS = ["TEMP","PRES","RH","RADIAT","BLH","RAIN","JFAC"]
    met_defaults = {"TEMP":298.15, "PRES":1013.25, "RH":50.0, "RADIAT":0.0,
                    "BLH":1000.0, "RAIN":0.0, "JFAC":1.0}
    mets = read_kv(inp / "meteo.config")
    met_c, met_s = {}, {}
    series = {}
    sidx = {}  # 序列名→下标

    def register_series(name, t, v):
        if name in sidx: return sidx[name]
        if len(sidx) >= 8:
            raise ValueError(
                f"序列总数 {len(sidx)+1} > Fortran 端上限 8 (runtime_dat_mod NSER_MAX=8; "
                f"气象+BC 共享同一序列池) — 请合并或降采样: '{name}' 注册被拒绝")
        sidx[name] = len(sidx) + 1
        series[sidx[name]] = (name, t, v)

    for k in METS:
        val = mets.get(k)
        if val is None:
            met_c[k] = met_defaults[k]; met_s[k] = 0
        elif val.startswith("file:"):
            fn = inp / val[5:]
            t, v = read_series(fn)
            if k == "BLH" and min(v) <= 0:
                raise ValueError(f"{fn}: BLH 序列含非正值 (min={min(v):g}) — "
                                 f"物理账②④ 中 vd/BLH、dB/BLH 将除零")
            register_series(k, t, v)
            met_s[k] = sidx[k]; met_c[k] = 0.0
        else:
            met_c[k] = float(val); met_s[k] = 0

    # ============ ⑥ BC ============
    bcc = read_kv(inp / "bc.config")
    tau_mix = float(bcc.pop("TAU_MIX", "86400"))
    bc_c = [0.0]*ns; bc_s = [0]*ns
    for k, v in bcc.items():
        if k not in sid:
            errs.append(f"BC 未知物种 '{k}'")
            continue
        if v.startswith("file:"):
            nm = "bc_" + k
            t, v2 = read_series(inp / v[5:])
            register_series(nm, t, v2)
            bc_s[sid[k]] = sidx[nm]
        else:
            parts = v.split()
            try:
                bc_c[sid[k]] = float(parts[0]) * unit_factor(parts[1] if len(parts) > 1 else "cm3")
            except ValueError as e:
                errs.append(f"BC '{k}': {e} (写法: 值 [单位], 如 '40 ppb')")

    # ============ ⑦ 排放 ============
    emis = [(0.0, "FLAT", 0.0)] * ns
    emit_unit = "molec/cm3/s"
    for ln in open(inp / "emissions.config"):
        c = ln.split("#", 1)[1].strip() if "#" in ln else ""
        if c.lower().startswith("unit:"):
            eu = c.split(":", 1)[1].strip().lower()
            if eu in ("ppb/h", "ppbh"):
                emit_unit = "ppb/h"
            elif eu in ("molec/cm3/s", "molec/cm³/s", "cm3/s"):
                emit_unit = "molec/cm3/s"
            else:
                errs.append(f"排放未知单位指令 '{eu}' (molec/cm3/s | ppb/h)")
            continue
    for ln in open(inp / "emissions.config"):
        s2 = strip_comment(ln)
        if not s2: continue
        p = s2.split()
        if p[0].upper() == "UNIT": continue
        if p[0] not in sid:
            errs.append(f"排放未知物种 '{p[0]}'")
            continue
        e0 = float(p[1]) * (PPB / 3600.0 if emit_unit == "ppb/h" else 1.0)
        amp = float(p[3]) if len(p) > 3 else 0.0
        shp = p[2] if len(p) > 2 else "FLAT"
        if e0 < 0:
            errs.append(f"排放 '{p[0]}': 基线速率 {e0:g} < 0 (负排放=非物理汇)")
        if shp not in ("BOX", "FLAT"):
            errs.append(f"排放 '{p[0]}': 日型 '{shp}' 非法 (只支持 BOX|FLAT, 大小写敏感; "
                        f"否则运行时会被静默按 FLAT 处理)")
        if shp == "BOX" and not (0.0 <= amp <= 1.0):
            errs.append(f"排放 '{p[0]}': BOX 振幅 {amp:g} ∉ [0,1] "
                        f"(amp<0 → 白天负排放; amp>1 → 峰值超基线无依据)")
        emis[sid[p[0]]] = (e0, shp, amp)

    # ============ ⑧ 干沉降 ============
    dry = [(0.0, 0)] * ns
    for ln in open(inp / "drydep.config"):
        s = strip_comment(ln)
        if not s: continue
        p = s.split()
        if p[0] not in sid:
            errs.append(f"干沉降未知物种 '{p[0]}'")
            continue
        dry[sid[p[0]]] = (float(p[1]), 1 if len(p)>2 and p[2]=="on" else 0)

    # ============ ⑨ 湿沉降 ============
    wet = [(0.0, 0.0)] * ns
    for ln in open(inp / "wetdep.config"):
        s = strip_comment(ln)
        if not s: continue
        p = s.split()
        if p[0] not in sid:
            errs.append(f"湿沉降未知物种 '{p[0]}'")
            continue
        wet[sid[p[0]]] = (float(p[1]), float(p[2]))

    # ============ ⑩ 稀释 ============
    dil = read_kv(inp / "dilution.config").get("DIL", "0.0")
    dil_v = float(dil) if not dil.startswith("file:") else 0.0

    # ============ ⑪ 情景 ============
    scen = []
    for ln in open(inp / "scenarios.config"):
        s = strip_comment(ln)
        if not s: continue
        name_part, _, mod_part = s.partition("=")
        head = name_part.strip()
        if not mod_part.strip(): continue
        d = dict(E_NOx=1.0, E_VOC=1.0, E_CO=1.0, E_CH4=1.0, JFAC=1.0, DT=0.0)
        for m in mod_part.split():
            mm = re.match(r"(E_NOx|E_VOC|E_CO|E_CH4|JFAC)\*([\d.eE+-]+)|(DT)([+-][\d.]+)", m)
            if not mm:
                errs.append(f"情景 '{head}': 无法解析 '{m}'")
                continue
            if mm.group(3):
                d["DT"] += float(mm.group(4))
            else:
                d[mm.group(1)] *= float(mm.group(2))
        scen.append((head, d))

    # ============ ⑫ 模型参数 ============
    mp = []
    for ln in open(inp / "model.parameters"):
        s = strip_comment(ln)
        if s:
            mp.append(float(s.split()[0]))
        if len(mp) == 4: break

    # ============ R16: 死配置告警 (RADIAT 当前无消费者) ============
    if abs(met_c.get("RADIAT", 0.0)) > 0.0 and met_s.get("RADIAT", 0) == 0:
        print("[cfgbuild] ⚠ 注意: meteo.config 的 RADIAT ≠ 0, 但当前无任何模块消费该量 "
              "(保留给 B2 在线光解集成位) — 不影响本次模拟结果", file=sys.stderr)

    # ============ 校验汇总 ============
    if errs:
        print("[cfgbuild] 配置错误:")
        for e in errs:
            print(f"   ✗ {e}")
        sys.exit(1)

    # ============ 光解参数化 l/m/n 映射 ============
    # 机制 J 索引最大值 (从 .fac 扫描)
    fac_txt = fac.read_text() if fac else ""
    jmax = max([int(i) for i in re.findall(r"J<(\d+)>", fac_txt)] or [1])
    if not mapping:
        mapping = list(range(1, jmax + 1))
    nj = max(jmax, len(mapping))
    jlmn = []
    for pos, slot in enumerate(mapping, 1):
        v = jmap_pairs.get(pos, str(slot))
        if v.startswith("D"):
            jlmn.append((float(v[1:]), 0.79, 0.33))   # 直接值: 有机物典型 m/n
        else:
            jlmn.append(phot_tab.get(int(v), (0.0, 1.0, 0.0)))

    # ============ 光解完备性校验 (防静默零光解) ============
    # 机制 J<slot> 使用的每个槽位必须有非零 l (表参数或 JMAP 直接值)
    _fac_txt = fac.read_text() if fac else ""
    fac_slots = sorted(set(int(i) for i in re.findall(r"J<(\d+)>", _fac_txt)))
    missing_slots = [s for s in fac_slots if s > len(jlmn) or jlmn[s - 1][0] == 0.0]
    if missing_slots:
        errs.append(f"机制使用光解槽位无参数 (J={missing_slots}) — "
                    f"补 photolysis-rates 表或 JMAP D 直接值")
    if phot_files and not phot_tab:
        errs.append(f"{phot_files[0].name} 存在但未解析出任何槽位")

    if errs:
        print("[cfgbuild] 配置错误:")
        for e in errs:
            print(f"   ✗ {e}")
        sys.exit(1)

    # ============ 编译输出 ============
    L = []
    w = L.append
    w(f"# runtime.dat 编译自 {input_dir}/ — Fortran 顺序表控读, 勿手改")
    w(f"# 机制: {fac.name if fac else 'NONE'} ({ns} 物种)")
    w(f"RUN {ns} {nj} {len(series)} {len(scen)}")
    w("SPECIES " + " ".join(species))
    w(f"MODEL {mp[0]:.4g} {mp[1]:.4g} {mp[2]:.4g} {mp[3]:.4g}")
    w("IC")
    w("  " + " ".join(f"{y:.6e}" for y in y0))
    w("MET")
    for k in METS:
        if met_s[k] == 0:
            w(f"  C {met_c[k]:.6g}")
        else:
            w(f"  S {met_s[k]}")
    w("BC")
    w(f"  {tau_mix:.6g} " + " ".join(f"{x:.6e}" for x in bc_c) + " " +
      " ".join(str(x) for x in bc_s))
    w("PHOT")
    w(f"  {jmode} {jlat:.4g} {jdoy:.6g}")
    for j, (l_, m_, n_) in enumerate(jlmn, 1):
        w(f"  J {j} {l_:.6e} {m_:.6g} {n_:.6g}")
    w("EMIS")
    for (e0, shp, amp) in emis:
        w(f"  E {e0:.6e} {shp} {amp:.4g}")
    w("DRY")
    for (v, f_) in dry:
        w(f"  D {v:.6g} {f_}")
    w("WET")
    for (a, b) in wet:
        w(f"  W {a:.6e} {b:.6g}")
    w("DIL")
    w(f"  {dil_v:.6g}")
    w("SCEN")
    for name, d in scen:
        w(f"  {name} {d['E_NOx']:.4g} {d['E_VOC']:.4g} {d['E_CO']:.4g} "
          f"{d['E_CH4']:.4g} {d['JFAC']:.4g} {d['DT']:+.4g}")
    w("SERIES")
    for idx in sorted(series):
        nm, t, v = series[idx]
        w(f"  {nm} {len(t)}")
        w("   " + " ".join(f"{x:.6g}" for x in t))
        w("   " + " ".join(f"{x:.6g}" for x in v))
    Path(out_path).write_text("\n".join(L) + "\n")
    print(f"[cfgbuild] {input_dir}/ → {out_path}")
    print(f"  机制: {fac.name} ({ns} 物种, {len(phot_tab)} J 参数)")
    print(f"  序列: {list(sidx.keys())}")
    print(f"  情景: {[s[0] for s in scen]}")
    return 0

if __name__ == "__main__":
    args = sys.argv[1:]
    mech_override = None
    if "--mech" in args:
        i = args.index("--mech")
        mech_override = args[i + 1]
        args = args[:i] + args[i + 2:]
    inp = args[0] if args else "input_v3"
    out = args[1] if len(args) > 1 else "runtime_v3.dat"
    compile_cfg(inp, out, mech_override=mech_override)
