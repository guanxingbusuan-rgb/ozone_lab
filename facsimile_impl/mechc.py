#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mechc.py —— Day 35: FACSIMILE 子集 → Fortran 代码生成 (AtChem2/Day28 思想独立精简版)
生成: gen/oz_mech.f90 (q/p 速率代码) + gen/oz.reac / gen/oz.prod (计量数据)
支持: 命名速率互相引用、J<idx> 光解索引、TEMP/M/O2 环境、@幂、三反应物行
"""
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

@dataclass
class Mech:
    species: list = field(default_factory=list)
    rates: list = field(default_factory=list)      # (name, expr)
    reactions: list = field(default_factory=list)  # (rate_expr, lhs[(c,name)], rhs[(c,name)])
    ro2_pool: list = field(default_factory=list)   # RO2 池成员 (单成员时表达式可替换)

def parse_fac(path):
    m = Mech()
    lines = Path(path).read_text().splitlines()
    in_vars = False
    in_form = False
    for raw in lines:
        s = raw.strip()
        if not s or s.startswith("*"):
            continue
        if s == "VARIABLE":
            in_vars = True; continue
        # R15: 段标记可能行内连写 ("... CH3O2 ;FORMULA"), 先归一化再判定
        _seg = s
        if ";" in s:
            _tail = s.split(";", 1)[1].strip()
            if _tail.startswith("FORMULA"):
                _seg = "FORMULA"          # 行内紧跟 FORMULA 段头
        if _seg == "FORMULA":
            # 若本行还含 VARIABLE 物种 (如 "...  CH3O2 ;FORMULA"), 先收物种
            if in_vars:
                _heads = s.split(";", 1)[0].split()
                if _heads:
                    m.species += _heads
                in_vars = False
            in_form = True; continue
        if in_form:
            s2 = s.lstrip(";").strip()
            if s2 == "" or s2.startswith("FORMULA"):
                continue                      # R15: 段间重复声明, 忽略
            if s.endswith(";"):
                in_form = False
            continue
        if in_vars:
            clean = s.strip()
            # R15: 分号可在行内任意位置 (MCM 风格 "... CH3O2 ;FORMULA" —— 段间无换行)
            if ";" in clean:
                head = clean.split(";", 1)[0]
                new_sp = head.split()
                in_vars = False
            else:
                new_sp = clean.split()
            if new_sp:
                dup = sorted(set(new_sp) & set(m.species)) \
                      + sorted({x for x in new_sp if new_sp.count(x) > 1})
                if dup:
                    raise SystemExit(f"[mechc] VARIABLE 重复物种: {sorted(set(dup))} "
                                     "(ns 虚增 + 索引错位 — 修复机制文件)")
                m.species += new_sp
            continue
        s = s.rstrip(";").strip()
        if s.startswith("%"):
            body, eq = s[1:].split(":", 1)
            lhs, rhs = eq.split("=", 1)
            m.reactions.append((body.strip(), _side(lhs), _side(rhs)))
        else:
            name, expr = s.split("=", 1)
            name = name.strip()
            if name == "RO2":
                m.ro2_pool = expr.split()
                continue
            if any(nm == name for nm, _ in m.rates):
                raise SystemExit(f"[mechc] 重复速率定义: {name} "
                                 "(后者静默覆盖前者 — 旧 Br+O3 同类 bug)")
            m.rates.append((name, expr.strip()))
    return m

def _side(s):
    out = []
    for tok in s.split("+"):
        tok = tok.strip()
        if not tok:
            continue
        mm = re.match(r"^(\d+\.?\d*)\s+(\w+)$", tok)
        out.append((float(mm.group(1)), mm.group(2)) if mm else (1.0, tok))
    return out

RESERVED = {"TEMP", "M", "O2", "EXP", "LOG10", "H2O", "N2"}

def to_f(expr, m, jn):
    e = expr.replace("@", "**")
    e = re.sub(r"\*\*-(\d+\.?\d*)", r"**(-\1)", e)
    e = re.sub(r"J<(\d+)>", r"JARR(\1)", e)      # FACSIMILE 光解索引 → JARR 元素
    def repl(mm):
        t = mm.group(0)
        if t == "RO2":
            if len(m.ro2_pool) != 1:
                raise ValueError(f"RO2 池多成员 ({m.ro2_pool})")
            return f"Y({m.species.index(m.ro2_pool[0]) + 1})"   # 池浓度 (替换值不再重扫)
        if t == "N2":
            return "N2C"                       # 背景 N2 常数 (gen 内 parameter)
        if t in m.species:
            return f"Y({m.species.index(t) + 1})"   # 速率表达式中的物种浓度 (RO2 池展开)
        if t in RESERVED:
            return t
        if t in [n for n, _ in m.rates]:
            return f"q({[n for n,_ in m.rates].index(t)+1})"
        return t
    return re.sub(r"[A-Za-z_][A-Za-z0-9_]*", repl, e)

def convert(fac, gen_dir):
    m = parse_fac(fac)
    ns, nr, nq = len(m.species), len(m.reactions), len(m.rates)
    jmax = [int(i) for i in re.findall(r"J<(\d+)>",
            Path(fac).read_text())] or [1]
    jn = [f"J{i}" for i in range(1, max(jmax) + 1)]
    sid = {s: i + 1 for i, s in enumerate(m.species)}
    # 三反应物上限由机制行决定
    maxr = max(len(l) for _, l, _ in m.reactions)
    maxp = max(len(r) for _, _, r in m.reactions)
    L = []
    L.append(f"! ===== 自动生成: mechc.py < {fac} =====")
    L.append("module oz_mech")
    L.append("  implicit none")
    L.append(f"  integer, parameter :: ns={ns}, nr={nr}, nq={nq}, nj={len(jn)}, mre={maxr}, mpe={maxp}")
    L.append("  real(8), parameter :: N2C = 2.46d19   ! 背景 N2 (FACSIMILE 恒定物种)")
    L.append("contains")
    L.append("  subroutine update_rates(TEMP, M, O2, H2O, JARR, Y, q, p)")
    L.append(f"    real(8), intent(in)  :: TEMP, M, O2, H2O, JARR(nj), Y(ns)")
    L.append("    real(8), intent(out) :: q(nq), p(nr)")
    for i, (name, expr) in enumerate(m.rates):
        L.append(f"    q({i+1}) = {to_f(expr, m, jn)}   ! {name}")
    for k, (expr, lhs, rhs) in enumerate(m.reactions):
        cm = " : " + "+".join(n for _, n in lhs) + " = " + "+".join(n for _, n in rhs)
        L.append(f"    p({k+1}) = {to_f(expr, m, jn)}   !{cm}")
    L.append("  end subroutine update_rates")
    L.append("end module oz_mech")
    Path(gen_dir, "oz_mech.f90").write_text("\n".join(L) + "\n")
    # 计量数据只记演化物种: RESERVED(O2/M/H2O)浓度恒定, 由速率表达式承载
    with open(Path(gen_dir, "oz.reac"), "w") as f:
        for _, lhs, _ in m.reactions:
            idx = [sid[n] for _, n in lhs if n in sid]
            f.write(f"{len(idx)}  " + "  ".join(map(str, idx)) + "\n")
    with open(Path(gen_dir, "oz.prod"), "w") as f:
        for _, _lhs, rhs in m.reactions:   # 元组序=(rate,lhs,rhs): 第二位是反应物!
            toks = []
            for c, n in rhs:
                if n in sid:
                    toks += [str(sid[n]), f"{c:g}"]
            f.write(f"{len(toks)//2}  " + "  ".join(toks) + "\n")
    print(f"[mechc] {ns} 物种, {nr} 反应, {nq} 命名速率, nj={jn[-1][1:]} 光解 → {gen_dir}/")

if __name__ == "__main__":
    fac = sys.argv[1] if len(sys.argv) > 1 else "mechanism/ozone.fac"
    gen = sys.argv[2] if len(sys.argv) > 2 else "gen"
    Path(gen).mkdir(exist_ok=True)
    convert(fac, gen)
