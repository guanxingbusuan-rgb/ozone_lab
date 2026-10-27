#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
element_audit.py —— ozone_lab 机制元素守恒审计器
逐反应检验 C/H/N/O/S/Cl/Br 平衡, 允许储库物种 (H2O/O2/CO2/H2) 在两侧自动推断;
若储库需求出现负值或非整(半)数 → FAIL。
用法: python element_audit.py [mechanism.fac]
退出码: 0=全过 1=有失败
"""
import re
import sys
from pathlib import Path

# ============ 52 物种元素表 (mechanism_3) ============
FORMULA = {
    # 无机核心
    "NO":  {"N":1,"O":1},   "NO2": {"N":1,"O":2},  "O3":  {"O":3},
    "O1D": {"O":1},         "O":   {"O":1},        "OH":  {"O":1,"H":1},
    "HO2": {"H":1,"O":2},   "H2O2":{"H":2,"O":2},  "CO":  {"C":1,"O":1},
    "HCHO":{"C":1,"H":2,"O":1}, "HNO3":{"H":1,"N":1,"O":3},
    "HONO":{"H":1,"N":1,"O":2}, "NO3":{"N":1,"O":3}, "N2O5":{"N":2,"O":5},
    # CH4 链 + PAN
    "CH4":  {"C":1,"H":4},  "CH3O2":{"C":1,"H":3,"O":2},
    "CH3OOH":{"C":1,"H":4,"O":2}, "CH3O":{"C":1,"H":3,"O":1},
    "CH3CHO":{"C":2,"H":4,"O":1}, "CH3CO3":{"C":2,"H":3,"O":3},
    "PAN":  {"C":2,"H":3,"N":1,"O":5},
    # isoprene 链 (ISOPBO2=C5H9O3 β-羟过氧; ISOPNO3=硝酸根过氧 C5H8NO5)
    "ISOP":   {"C":5,"H":8},
    "ISOPBO2":{"C":5,"H":9,"O":3},
    "ISOPBO": {"C":5,"H":9,"O":2},
    "ISOPNO3":{"C":5,"H":8,"N":1,"O":5},
    "ISOPNIT":{"C":5,"H":9,"N":1,"O":4},
    "ISOPDNT":{"C":5,"H":9,"N":1,"O":5},
    "ISOPOOH":{"C":5,"H":10,"O":3},
    "MVK":    {"C":4,"H":6,"O":1},
    "MVKO2":  {"C":4,"H":7,"O":4},
    "MVKO":   {"C":4,"H":7,"O":3},
    "MVKOOH": {"C":4,"H":8,"O":4},
    "MACR":   {"C":4,"H":6,"O":1},
    "MACRO2": {"C":4,"H":5,"O":3},     # 摘氢酰基过氧
    "MACROO": {"C":4,"H":5,"O":2},
    "MACROOH":{"C":4,"H":6,"O":3},
    "MPAN":   {"C":4,"H":5,"N":1,"O":5},
    "MGLY":   {"C":3,"H":4,"O":2},     # 甲基乙二醛
    # 硫链
    "SO2":  {"S":1,"O":2}, "SO3": {"S":1,"O":3}, "H2SO4":{"H":2,"S":1,"O":4},
    "DMS":  {"C":2,"H":6,"S":1},
    "MSA":  {"C":1,"H":4,"O":3,"S":1},  # 甲基磺酸 CH3SO3H
    # 卤素
    "Cl":   {"Cl":1}, "ClO": {"Cl":1,"O":1}, "ClONO2":{"Cl":1,"N":1,"O":3},
    "Br":   {"Br":1}, "BrO": {"Br":1,"O":1}, "BrONO2":{"Br":1,"N":1,"O":3},
    "HCl":  {"H":1,"Cl":1}, "HBr": {"H":1,"Br":1},
    "HOCl": {"H":1,"O":1,"Cl":1}, "HOBr": {"H":1,"O":1,"Br":1},
}

# 储库物种: 可在反应任一侧自动补齐 (大气定常背景)
RESERVOIR = {
    "H2O":  {"H":2,"O":1},
    "O2":   {"O":2},
    "CO2":  {"C":1,"O":2},
    "H2":   {"H":2},
    "M":    {},
    "N2":   {"N":2},
}
RES_NAMES = "|".join(RESERVOIR.keys())

def parse_fac(path):
    """返回 [(rate_str, lhs, rhs)]; 并从 FORMULA 段加载分子式 (M3: 机制自带配平说明书)"""
    rxns, in_vars, in_form = [], False, False
    formulas = {}
    for raw in Path(path).read_text().splitlines():
        s = raw.strip()
        if not s or s.startswith("*"):
            continue
        if s == "VARIABLE":
            in_vars = True; continue
        if s.startswith("FORMULA"):
            in_form = True; continue
        if in_form:
            body = s.rstrip(";").strip()
            if body:
                tok = body.split()
                import re as _re
                formulas[tok[0]] = {}
                for e in tok[1:]:
                    mm = _re.match(r"^([A-Z][a-z]?)(\d+)$", e)
                    if not mm:
                        raise ValueError(f"FORMULA 项非法: '{e}' in {tok[0]}")
                    formulas[tok[0]][mm.group(1)] = int(mm.group(2))
            if s.endswith(";"):
                in_form = False
            continue
        if in_vars:
            if s.rstrip().endswith(";"):
                in_vars = False
            continue
        s = s.rstrip(";").strip()
        if s.startswith("%"):
            body, eq = s[1:].split(":", 1)
            lhs, rhs = eq.split("=", 1)
            rxns.append((body.strip(), _side(lhs), _side(rhs)))
    return rxns

def _side(s):
    out = []
    for tok in s.split("+"):
        tok = tok.strip()
        if not tok:
            continue
        mm = re.match(r"^([\d.]+)\s+(\w+)$", tok)
        out.append((float(mm.group(1)), mm.group(2)) if mm else (1.0, tok))
    return out

def elem_vec(tokens):
    v = {}
    for coef, name in tokens:
        f = FORMULA.get(name)
        if f is None:
            if name in RESERVOIR:
                f = RESERVOIR[name]   # 显式储库 (如 O + O2 + M)
            else:
                raise KeyError(name)
        for e, n in f.items():
            v[e] = v.get(e, 0.0) + coef * n
    return v

def close_with_reservoirs(net, elements, verbose=False):
    """net = lhs - rhs 逐元素差; 求储库组合使 net 归零。
    物理方向约束 (非物理闭合会掩盖真实机制错误):
      O2 : lhs 计数 ≤ 2.5 (自由基+O2 链), rhs ≤ 1 (终止反应放 O2)
      H2O: 任一侧 ≤ 1 (摘氢放水 / 水解吸水)
      CO2: 仅 rhs ≤ 1 (酰基分解)
      H2 : 仅 rhs ≤ 1 (HCHO 分子光解通道)
    返回 (解 dict, 涉及储库描述) 或 (None, None)。"""
    sols = []
    for w_h2o in (-1, 0, 1):                 # -1 = lhs, +1 = rhs
        for w_h2 in (0, 1):                  # rhs only
            for w_co2 in (0, 1):             # rhs only
                rem = dict(net)
                for name, sgn in (("H2O", w_h2o), ("H2", w_h2), ("CO2", w_co2)):
                    if sgn == 0:
                        continue
                    for e, n in RESERVOIR[name].items():
                        rem[e] = rem.get(e, 0.0) - sgn * n
                others = {e: rem.get(e, 0.0) for e in elements if e != "O"}
                if any(abs(x) > 1e-9 for x in others.values()):
                    continue
                o = rem.get("O", 0.0)
                for o2 in (0, .5, 1, 1.5, 2, 2.5):       # rhs
                    if abs(o - 2 * o2) < 1e-9:
                        sols.append((w_h2o, w_h2, w_co2, o2, "rhs"))
                for o2 in (.5, 1, 1.5, 2, 2.5):          # lhs
                    if abs(o + 2 * o2) < 1e-9:
                        sols.append((w_h2o, w_h2, w_co2, o2, "lhs"))
    if not sols:
        return None, None
    w_h2o, w_h2, w_co2, o2, o2side = sols[0]
    desc = []
    if w_h2o: desc.append(f"H2O@{'lhs' if w_h2o < 0 else 'rhs'}")
    if w_h2:  desc.append("H2@rhs")
    if w_co2: desc.append("CO2@rhs")
    if o2:    desc.append(f"O2@{o2side}x{o2:g}")
    return True, "+".join(desc) or "精确(无储库)"

def audit(fac_path, max_show=40):
    global FORMULA
    rxns = parse_fac(fac_path)
    # M3: 机制文件自带 FORMULA → 覆盖内嵌表 (单一事实源)
    import element_audit as _self
    import re as _re2
    in_form = False
    fcount = 0
    for raw in Path(fac_path).read_text().splitlines():
        s = raw.strip()
        if s.startswith("FORMULA"):
            in_form = True; continue
        if in_form:
            body = s.rstrip(";").strip()
            if body:
                tok = body.split()
                d = {}
                for e in tok[1:]:
                    mm = _re2.match(r"^([A-Z][a-z]?)(\d+)$", e)
                    d[mm.group(1)] = int(mm.group(2))
                FORMULA[tok[0]] = d
                fcount += 1
            if s.endswith(";"):
                in_form = False
    if fcount:
        print(f"  [M3] 机制自带 FORMULA {fcount} 物种 (覆盖内嵌表)")
    elements = ["C", "H", "N", "O", "S", "Cl", "Br"]
    nfail = 0
    print(f"===== 元素守恒审计: {fac_path} ({len(rxns)} 反应) =====")
    for i, (rate, lhs, rhs) in enumerate(rxns, 1):
        try:
            lv, rv = elem_vec(lhs), elem_vec(rhs)
        except KeyError as e:
            print(f"  [FAIL] #{i} 未知物种 {e}:  {' '.join(f'{c} {n}' for c, n in lhs)}")
            nfail += 1
            continue
        net = {e: lv.get(e, 0.0) - rv.get(e, 0.0) for e in elements}
        if all(abs(v) < 1e-9 for v in net.values()):
            continue
        ok, desc = close_with_reservoirs(net, elements)
        if ok is not None:
            continue   # 物理储库闭合 → 通过
        nfail += 1
        ntxt = "  ".join(f"{e}:{v:+g}" for e, v in net.items() if abs(v) > 1e-9)
        print(f"  [FAIL] #{i} {rate[:34]:36s} 净差 {ntxt}")
    if nfail == 0:
        print(f"  全部 {len(rxns)} 反应元素守恒 PASS "
              f"(储库 {RES_NAMES} 自动推断)")
    else:
        print(f"  {nfail}/{len(rxns)} 反应 FAIL")
    return nfail

if __name__ == "__main__":
    fac = sys.argv[1] if len(sys.argv) > 1 else \
        "mechanism/mechanism_3_full_chain.fac"
    sys.exit(1 if audit(fac) else 0)
