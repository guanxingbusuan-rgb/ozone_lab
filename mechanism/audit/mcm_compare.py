#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mcm_compare.py —— 与 MCM v3.3.1 官方导出逐反应比对速率系数

用途:
  本机制（mechanism_3_full_chain.fac）的速率常数大量来自 MCM/JPL，但手工
  录入难免偏差。本工具用官方导出（mechanism/reference/mcm_v3.3.1_export.fac）
  做数值校验：对每条能按"反应物集合"唯一匹配上的反应，在标准条件
  (T=298.15 K, M=2.46e19 molec/cm3) 下计算两边速率常数并报告相对偏差。

用法:
  python mechanism/audit/mcm_compare.py                      # 全部报告，按偏差排序
  python mechanism/audit/mcm_compare.py --threshold 0.05     # 只显示 >5% 偏差
  python mechanism/audit/mcm_compare.py --quiet              # 仅汇总（供 CI）

退出码: 0 = 最大偏差 < 阈值; 1 = 存在超阈值项（便于调用方按码判断）

注意:
  · 只比"反应物集合完全一致"的反应；集总类（OLE/ALK/ARO1/TERP…）在 MCM
    里没有对应物种，天然不参与。
  · 光解反应（J<n>）不比数值（MCM 的 J 槽位要靠光解表，不在这里）。
  · 三体/Troe 反应两边都用 M 或 N2 表达式，按同条件求值即可。
"""
import os
import re
import sys
import argparse
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OURS = os.path.join(ROOT, "mechanism", "mechanism_3_full_chain.fac")
MCM = os.path.join(ROOT, "mechanism", "reference", "mcm_v3.3.1_export.fac")

# 标准比对条件
T_REF = 298.15
M_REF = 2.46e19
O2_REF = 0.2095 * M_REF
N2_REF = 0.7805 * M_REF
H2O_REF = 0.0          # 干燥，避免水汽增强因子掩盖主参数差异

# 集总物种 → MCM 代表分子（单一事实源：mcm_compare 关口② 与 jslot_audit 共用）
LUMP_MAP = {
    "ALK": "C3H8", "OLE": "C3H6", "ARO1": "TOLUENE", "ARO2": "MXYL",
    "KET": "CH3COCH3", "TERP": "APINENE", "SESQ": "BCARY",
    "ISOP": "C5H8", "CH4": "CH4", "LALK": "NC6H14",
    "ETH": "C2H4", "IOLE": "CBUT2ENE", "ALD2": "C2H5CHO",
    "MEOH": "CH3OH", "ETOH": "C2H5OH", "ESTER": "ETHACET",
    "MGLY": "MGLYOX", "GLY": "GLYOX",
}


# ══════════════════ .fac 解析 ══════════════════
def parse_fac(path):
    """返回 (species, reactions, rates) 三元组

    species  : 物种名列表（VARIABLE 段）
    reactions: [(rate_token, frozenset(lhs_names), tuple(rhs_names))]
    rates    : {name: expr} 命名速率
    """
    species, reactions, rates = [], [], {}
    seg = ""
    for raw in open(path):
        s = raw.strip()
        if not s or s.startswith("*"):
            continue
        if s == "VARIABLE":
            seg = "VAR"; continue
        if s.startswith("FORMULA"):
            seg = "FORM"; continue
        if seg == "VAR":
            body = s.rstrip(";").strip()
            species += body.split()
            if s.endswith(";"):
                seg = ""
            continue
        if seg == "FORM":
            if s.endswith(";"):
                seg = ""
            continue
        s2 = s.rstrip(";").strip()
        if not s2:
            continue
        if s2.startswith("%"):
            if ":" not in s2 or "=" not in s2:
                continue
            body, eq = s2[1:].split(":", 1)
            lhs_s, rhs_s = eq.split("=", 1)

            def side(x):
                out = []
                for tok in x.split("+"):
                    tok = tok.strip()
                    if not tok:
                        continue
                    p = tok.split()
                    out.append(p[1] if len(p) == 2 else p[0])
                return out

            reactions.append((body.strip(), frozenset(side(lhs_s)), tuple(side(rhs_s))))
        elif "=" in s2:
            name, expr = s2.split("=", 1)
            name = name.strip()
            # 去行内 ! 注释, 再去分号 (顺序重要: 注释前的 ; 才是语句终结符)
            expr = expr.split("!", 1)[0].strip().rstrip(";").strip()
            if name == "RO2":
                continue
            rates[name] = expr
    return species, reactions, rates


# ══════════════════ 表达式求值 ══════════════════
# 只匹配"不在数字内部"的标识符: 排除 6.0E-34 里被误认的 E
_TOK = re.compile(r"(?<![\d.])[A-Za-z_][A-Za-z0-9_]*")


def to_py(expr):
    """MCM/FACSIMILE 表达式 → Python 可 eval 形式

    关键: D 指数只在"数字之后紧跟 D 再跟数字"时才换成 E
    —— 不能全局替换，否则 KO1DH2O 会被改成 KO1EH2O（真实踩过的坑）。
    """
    e = expr.replace("@", "**")
    e = re.sub(r"(?<=\d)[Dd]([+-]?\d+)", r"e\1", e)     # 3.45D-12 → 3.45e-12
    return e


def evaluate(expr, table, depth=0, cache=None):
    """递归求值：先把命名速率展开，再交给 eval

    table: {name: expr} 命名速率表（本机制或 MCM 各自的）
    """
    if cache is None:
        cache = {}
    if depth > 20:
        raise ValueError("命名速率展开过深（循环引用？）")

    def repl(m):
        name = m.group(0)
        if name in ("TEMP", "M", "O2", "N2", "H2O", "EXP", "LOG10", "LOG"):
            return name
        if name in cache:
            return repr(cache[name])
        if name in table:
            val = evaluate(table[name], table, depth + 1, cache)
            cache[name] = val
            return repr(val)
        raise KeyError(name)

    py = _TOK.sub(repl, to_py(expr))
    env = {
        "TEMP": T_REF, "M": M_REF, "O2": O2_REF, "N2": N2_REF, "H2O": H2O_REF,
        "EXP": __import__("math").exp,
        "LOG10": __import__("math").log10,
        "LOG": __import__("math").log,
    }
    return eval(py, {"__builtins__": {}}, env)


# ══════════════════ 主流程 ══════════════════
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=0.10,
                    help="报告阈值（相对偏差），默认 0.10")
    ap.add_argument("--quiet", action="store_true", help="仅汇总")
    args = ap.parse_args(argv)

    for p in (OURS, MCM):
        if not os.path.exists(p):
            print(f"✗ 缺文件: {p}")
            return 2

    our_sp, our_rx, our_rates = parse_fac(OURS)
    mcm_sp, mcm_rx, mcm_rates = parse_fac(MCM)

    # MCM 按反应物集合索引
    midx = defaultdict(list)
    for rate, lhs, rhs in mcm_rx:
        midx[lhs].append((rate, rhs))

    # 第二道关口用模块级 LUMP_MAP（集总物种 → MCM 代表分子）：
    # 纯无机比对只覆盖 57 条，集总物种的速率曾被漏查，
    # 实测抓出 9 处错误（SESQ+O3 偏慢 100x、MGLY+OH 漏 EXP 项偏慢 7x 等）。
    PHOT = re.compile(r"^J<\d+>$")

    rows, struct, skipped = [], [], 0
    # 本机制可能把同一反应物集合拆成多条（如 HO2+HO2 的双通道 A/B），
    # 与 MCM 的"分通道之和"比较才公平 → 先按反应物集合聚合本机制速率。
    our_sum = defaultdict(float)
    our_expr = defaultdict(list)
    for rate, lhs, rhs in our_rx:
        if PHOT.match(rate):
            continue
        try:
            our_sum[lhs] += evaluate(rate, our_rates)
            our_expr[lhs].append(rate)
        except Exception:
            skipped += 1
    _grouped = set()
    for rate, lhs, rhs in our_rx:
        if PHOT.match(rate) or lhs in _grouped:
            continue
        _grouped.add(lhs)
        v_our = our_sum[lhs]              # 同集合多条已求和
        rate = "+".join(our_expr[lhs])    # 展示用
        hits = midx.get(lhs)
        if not hits:
            continue
        # 分两类比较（关键，否则全是假阳性）:
        #  (A) 产物集合完全一致的 MCM 分支 → 严格数值比对（可判定录入错误）
        #  (B) 只有反应物集合相同、产物不同 → 属"集总化/分支取舍"，
        #      单独归类，不判 FAIL（本机制用单产物代表 MCM 的多产物分支）
        exact = [(r, h) for r, h in ((None, h) for h in hits)
                 if frozenset(h[1]) == frozenset(rhs)]
        v_exact, used_exact = 0.0, []
        for _r, (mrate, mrhs) in exact:
            try:
                v = evaluate(mrate, mcm_rates)
            except Exception:
                continue
            if v > 0:
                v_exact += v
                used_exact.append(mrate)
        # 全部命中的速率和（用于 (B) 类说明）
        v_all, used_all = 0.0, []
        for mrate, mrhs in hits:
            try:
                v = evaluate(mrate, mcm_rates)
            except Exception:
                continue
            if v > 0:
                v_all += v
                used_all.append(mrate)

        if v_exact > 0:
            rows.append({
                "our_rate": rate, "lhs": " + ".join(sorted(lhs)),
                "v_our": v_our, "v_mcm": v_exact,
                "mcm_rate": (" | ".join(dict.fromkeys(used_exact)))[:60],
                "n_branch": len(used_exact), "rel": abs(v_our - v_exact) / v_exact,
                "kind": "exact",
            })
        elif v_all > 0:
            struct.append({
                "lhs": " + ".join(sorted(lhs)), "v_our": v_our,
                "v_mcm": v_all, "our_rhs": " + ".join(rhs),
                "mcm_rhs": " | ".join(" + ".join(h[1]) for h in hits)[:70],
                "n_branch": len(hits),
            })

    rows.sort(key=lambda r: -r["rel"])
    over = [r for r in rows if r["rel"] > args.threshold]

    if not args.quiet:
        print("=" * 78)
        print(f"MCM v3.3.1 速率系数比对  (T={T_REF} K, M={M_REF:.2e})")
        print(f"本机制 {len(our_rx)} 反应 / MCM {len(mcm_rx)} 反应")
        print("=" * 78)
        print(f"\n可比对 {len(rows)} 条（跳过量纲不符 {skipped} 条）")
        if rows:
            print(f"最大偏差 {rows[0]['rel'] * 100:.1f}%  "
                  f"中位偏差 {sorted(r['rel'] for r in rows)[len(rows) // 2] * 100:.2f}%")
        print(f"\n超过 {args.threshold:.0%} 阈值的 {len(over)} 条:\n")
        print(f"{'反应物':30s} {'本机制':>12s} {'MCM合计':>12s} {'偏差':>8s} 分支  表达式")
        print("-" * 78)
        for r in over:
            print(f"{r['lhs'][:30]:30s} {r['v_our']:12.3e} {r['v_mcm']:12.3e} "
                  f"{r['rel'] * 100:7.1f}% {r.get('n_branch', 1):4d}  {r['our_rate'][:20]}")
            print(f"{'':32s}{'':13s}{'':13s}{'':9s}{'':5s} MCM: {r['mcm_rate']}")

    # ══ 第二道关口: 集总物种 → MCM 代表分子（跨机制可映射比对）══
    lump_rows = []
    for rate, lhs, rhs in our_rx:
        if PHOT.match(rate) or not (lhs & set(LUMP_MAP)) or not (1 <= len(lhs) <= 2):
            continue
        mapped = frozenset(LUMP_MAP.get(x, x) for x in lhs)
        if not mapped <= set(mcm_sp):
            continue
        hits2 = midx.get(mapped)
        if not hits2:
            continue
        try:
            v_our2 = evaluate(rate, our_rates)
        except Exception:
            continue
        v_sum2, used2 = 0.0, []
        for mr, mrhs in hits2:
            try:
                v2 = evaluate(mr, mcm_rates)
            except Exception:
                continue
            if v2 > 0:
                v_sum2 += v2
                used2.append(mr)
        if v_sum2 <= 0:
            continue
        lump_rows.append({
            "lhs": " + ".join(sorted(lhs)),
            "mapped": " + ".join(sorted(mapped)),
            "v_our": v_our2, "v_mcm": v_sum2,
            "rel": abs(v_our2 - v_sum2) / v_sum2,
            "our_rate": rate, "n_branch": len(used2),
        })
    lump_rows.sort(key=lambda r: -r["rel"])
    lump_over = [r for r in lump_rows if r["rel"] > args.threshold]

    if not args.quiet and lump_rows:
        print(f"\n{'=' * 78}")
        print(f"集总物种跨机制比对 {len(lump_rows)} 条（映射到 MCM 代表分子）")
        print("=" * 78)
        print(f"{'本机制':26s} {'MCM 代表':22s} {'本值':>11s} {'MCM':>11s} {'偏差':>8s}")
        print("-" * 78)
        for r in lump_over:
            print(f"{r['lhs'][:26]:26s} {r['mapped'][:22]:22s} "
                  f"{r['v_our']:11.3e} {r['v_mcm']:11.3e} {r['rel'] * 100:7.1f}%")
        if not lump_over:
            print(f"  (全部在 {args.threshold:.0%} 内，最大 "
                  f"{lump_rows[0]['rel'] * 100:.1f}%)")

    if not args.quiet and struct:
        print(f"\n{'=' * 78}")
        print(f"结构性差异 {len(struct)} 条（集总化/分支取舍 —— 非录入错误）")
        print("=" * 78)
        print("本机制用【单产物】代表 MCM 的【多分支】，速率按产物集合不完全对应，")
        print("故不做数值判定；列出以便核对产物方向是否符合设计意图。\n")
        for r in struct:
            print(f"  {' + '.join(r['lhs'].split(' + '))[:28]:28s} "
                  f"本={r['v_our']:.3e}  MCM合计={r['v_mcm']:.3e} ({r['n_branch']} 支)")
            print(f"    本: → {r['our_rhs'][:60]}")
            print(f"    MCM: → {r['mcm_rhs']}")

    # ══ 第三道关口: 文档一致性（mechanism/README 反应表 vs .fac 实际值）══
    doc_bad = []
    _rd_path = os.path.join(ROOT, "mechanism", "README.md")
    if os.path.exists(_rd_path):
        _rd = open(_rd_path).read()
        _defs = {k: v.strip() for k, v in
                 re.findall(r"^([A-Z][A-Za-z0-9_]*)\s*=\s*([^;]+);",
                            open(OURS).read(), re.M)}
        _rows = re.findall(r"`([A-Z][A-Za-z0-9_]*)` = ([^|`\n]+)", _rd)
        _norm = lambda x: x.replace(" ", "").replace("D", "e").replace("d", "e").upper()
        for _n, _v in _rows:
            _v = _v.strip().rstrip("| ").strip()
            if _n in _defs:
                if (_norm(_v) not in _norm(_defs[_n])
                        and _norm(_defs[_n]) not in _norm(_v)):
                    doc_bad.append((_n, _v, _defs[_n]))
    if not args.quiet and doc_bad:
        print(f"\n{'=' * 78}")
        print(f"文档漂移 {len(doc_bad)} 处（mechanism/README 与 .fac 不一致）")
        print("=" * 78)
        for _n, _v, _a in doc_bad:
            print(f"  {_n:12s} README={_v[:26]:26s} 实际={_a[:30]}")

    max_rel = max([r["rel"] for r in rows + lump_rows] or [0])
    ok = max_rel <= args.threshold and not doc_bad
    print(f"\n[{'PASS' if ok else 'FAIL'}] 速率系数一致性 "
          f"({len(rows)} 条无机 + {len(lump_rows)} 条集总映射, "
          f"最大偏差 {max_rel * 100:.1f}% vs 阈值 {args.threshold:.0%})"
          if (rows or lump_rows) else "[SKIP] 无可比对反应")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
