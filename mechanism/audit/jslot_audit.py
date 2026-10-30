#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
jslot_audit.py —— 光解槽位语义审计（机制质检第三维度）

为什么单独立一个维度:
  mcm_compare.py 比的是热反应速率数值；光解 J<n> 的 l/m/n 完全由**槽号**决定，
  而槽号是 MCM 官方命名空间里的"外键"。借错槽 = 借了别的分子的 l/m/n，
  速率数值比对对此毫无感知（J 值根本不参与求值）。V2.5 曾因此迁移 7 处
  （见 mechanism/V25_ASSESSMENT.md §2.2）。本脚本把那次人工逐槽核对固化为断言。

四条不变量:
  ① 语义包含 —— 机制占用的每个 MCM 编号槽，其反应物（剔储库 + 集总映射后）
     必须 ⊆ MCM 官方该槽的反应物分子集合。
  ② 自定义槽封闭 —— MCM 无对应分子的槽必须显式落在白名单 CUSTOM 里，
     且 MCM 该槽确实为空（否则就是借用了别人占着的槽），
     且 photolysis.config 的 JMAP 给了 D<正午值> 直达参数。
  ③ JMAP 卫生 —— `photolysis.config` 的 JMAP 里每个槽都必须真的被机制引用
     （无死参数），且**只允许**给白名单里的自定义槽（不得悄悄覆盖 MCM 官方 l/m/n）。
  ④ 文档一致 —— mechanism/README.md 每条 J<n> 行的槽号/反应物/产物/参数来源
     必须与 .fac + photolysis.config 实际一致（§5 反应表与 §6 光解方案表都查）。

判定口径（防假阳/假阴的三个设计）:
  · 只取 MCM 中速率 token **恰为** `J<n>`（无产率乘子）的反应作为该槽的官方
    分子归属。MCM 大量用 `J<4>*0.14*0.4` 这类写法把 J4 的**日变化形状**借给
    别的分子；若把这些也算进槽归属，不变量①会退化成永真。
  · 反应物侧的 O2/H2O/CO2/H2/N2/M 是本机制的记账储库（MCM 不写在反应物侧），
    比对前剔除；储库口径与 element_audit.py 共用同一张表。
  · J34/J41 这类"通用槽"在 MCM 里本就归属上百个分子，①在那里天然宽松；
    报告会打印该槽的 MCM 分子数，让读者知道这条断言的实际强度。

用法:
  python mechanism/audit/jslot_audit.py            # 详细报告
  python mechanism/audit/jslot_audit.py --quiet    # 仅结论
退出码: 0 = 全部通过; 1 = 有 FAIL; 2 = 缺输入文件
"""
import argparse
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from element_audit import RESERVOIR                       # noqa: E402
from mcm_compare import LUMP_MAP, MCM, OURS, parse_fac    # noqa: E402

ROOT = os.path.dirname(os.path.dirname(HERE))
README = os.path.join(ROOT, "mechanism", "README.md")
PHOTCFG = next((q for q in (
    os.path.join(ROOT, "configs", "photolysis.config"),
    os.path.join(ROOT, "CASE-analysis", "config_override", "photolysis.config"),
    os.path.join(ROOT, "EKMA-emission", "config_override", "photolysis.config"),
) if os.path.exists(q)), os.path.join(ROOT, "configs", "photolysis.config"))

# MCM 官方导出里没有对应分子 → photolysis.config JMAP 走文献正午值。
# 白名单封闭: 新增自定义槽必须显式写在这里，否则不变量② 判 FAIL。
CUSTOM = {
    25: "BrO", 26: "ClONO2", 27: "BrONO2", 28: "HOCl",
    29: "HOBr", 30: "N2O5", 45: "ClO",
}

# 光解专属代理（MCM 写显式异构体，本机制写集总名）。
# ISOPOOH 是异戊二烯-OOH 的集总；MCM 六个异构体（ISOPA/B/C/D/34/N-OOH）
# 全部走 J41，取 ISOPBOOH 作代表（与 ISOP+OH → ISOPBO2 的主路代表一致）。
PHOT_PROXY = {"ISOPOOH": "ISOPBOOH"}
PROXY = dict(LUMP_MAP, **PHOT_PROXY)

JSLOT = re.compile(r"^J<(\d+)>")
ROW = re.compile(r"^\|(.+)\|$")


def proxy(names):
    """机制物种名 → MCM 代表分子集合（剔除记账储库）"""
    return frozenset(PROXY.get(n, n) for n in names if n not in RESERVOIR)


def mcm_slot_molecules(mcm_rx):
    """→ (官方归属, 借形状) 两张 {slot: set(反应物分子)}"""
    exact, borrowed = {}, {}
    for tok, lhs, _rhs in mcm_rx:
        m = JSLOT.match(tok)
        if not m:
            continue
        s = int(m.group(1))
        tgt = exact if tok.strip() == f"J<{s}>" else borrowed
        tgt.setdefault(s, set()).update(lhs)
    return exact, borrowed


def parse_jmap(path):
    m = re.search(r"^JMAP\s*=\s*([^#\n]+)", open(path, encoding="utf-8").read(), re.M)
    out = {}
    for tok in (m.group(1).split(",") if m else []):
        tok = tok.strip()
        if tok:
            k, v = tok.split("=", 1)
            out[int(k)] = v.strip()
    return out


def parse_rxn(text):
    def side(x):
        out = []
        for tok in x.split("+"):
            tok = tok.strip()
            if tok:
                p = tok.split()
                out.append(p[1] if len(p) == 2 else p[0])
        return out
    lhs, rhs = text.split("=", 1)
    return side(lhs), side(rhs)


def rxn_of(cells):
    """取该行反引号里的反应式（允许反引号后跟注释文字）"""
    for c in reversed(cells):
        m = re.search(r"`([^`]*=[^`]*)`", c)
        if m and "J<" not in m.group(1):
            return m.group(1)
    return None


def scan_readme(path, known):
    """扫出 README 里所有光解表行 → [(行号, 槽号, 反应式, 参数来源列 or None)]

    两种表都收: §5 反应表把 `J<n>` 写在速率列; §6 光解方案表把裸槽号写在首列、
    末列是参数来源。只保留反应式里每个物种都是机制物种或记账储库的行
    （滤掉无关表格）。
    """
    allowed = set(known) | set(RESERVOIR)
    rows = []
    for ln, raw in enumerate(open(path, encoding="utf-8").read().splitlines(), 1):
        m = ROW.match(raw.strip())
        if not m:
            continue
        cells = [c.strip() for c in m.group(1).split("|")]
        rx = rxn_of(cells)
        if rx is None:
            continue
        jm = re.search(r"`J<(\d+)>`", " ".join(cells))
        if jm:
            slot = int(jm.group(1))
            src = None
        elif cells[0].isdigit() and len(cells) >= 6:
            slot, src = int(cells[0]), cells[-1]
        else:
            continue
        lhs, rhs = parse_rxn(rx)
        if not lhs or not all(n in allowed for n in lhs + rhs):
            continue
        rows.append((ln, slot, rx, src))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet", action="store_true", help="仅结论")
    args = ap.parse_args(argv)

    for p in (OURS, MCM, README, PHOTCFG):
        if not os.path.exists(p):
            print(f"✗ 缺文件: {p}")
            return 2

    our_sp, our_rx, _ = parse_fac(OURS)
    _msp, mcm_rx, _ = parse_fac(MCM)
    exact, borrowed = mcm_slot_molecules(mcm_rx)
    jmap = parse_jmap(PHOTCFG)
    known = set(our_sp)

    phot = [(int(JSLOT.match(t).group(1)), t, lhs, rhs)
            for t, lhs, rhs in our_rx if JSLOT.match(t)]
    by_slot = {}
    for s, _t, lhs, rhs in phot:
        by_slot.setdefault(s, []).append((lhs, rhs))

    fails, notes = [], []

    def fail(msg):
        fails.append(msg)

    if not args.quiet:
        print("=" * 78)
        print("光解槽位语义审计  (机制 J<n> vs MCM v3.3.1 官方槽归属)")
        print(f"机制 {len(phot)} 条光解 / {len(by_slot)} 个不同槽位 / "
              f"最高 J{max(by_slot)};  MCM 官方归属槽 {len(exact)} 个")
        print("=" * 78)
        print(f"{'槽':>4} {'机制反应物':22s} {'→ MCM 代表':22s} 判定")
        print("-" * 78)

    # ══ ①② 语义包含 / 自定义槽封闭 ══
    for s in sorted(by_slot):
        for lhs, _rhs in by_slot[s]:
            mine = proxy(lhs)
            if s in CUSTOM:
                if exact.get(s):
                    fail(f"J{s} 声明为自定义槽，但 MCM 官方该槽已归属 "
                         f"{sorted(exact[s])[:6]} → 借用了被占槽位")
                if CUSTOM[s] not in mine:
                    fail(f"J{s} 白名单声明分子 {CUSTOM[s]}，机制反应物为 {sorted(mine)}")
                if not jmap.get(s, "").startswith("D"):
                    fail(f"J{s} 是自定义槽但 JMAP 缺 D<正午值> 直达参数 "
                         f"(实际 {jmap.get(s, '无')!r})")
                verdict = f"自定义槽（MCM 无基准，JMAP {jmap.get(s, '?')}）"
            else:
                if not exact.get(s):
                    fail(f"J{s} 不在自定义白名单，且 MCM 官方该槽无归属反应 "
                         f"→ 应加入 CUSTOM 或改用正确槽号")
                    verdict = "无 MCM 归属（未声明）"
                elif not mine <= exact[s]:
                    bad = sorted(mine - exact[s])
                    hint = "; ".join(
                        f"{b} 在 MCM 走 J{sorted(k for k, v in exact.items() if b in v)}"
                        if any(b in v for v in exact.values()) else f"{b} 在 MCM 无光解"
                        for b in bad)
                    fail(f"J{s} 语义错位: 机制反应物 {sorted(mine)} ⊄ MCM 该槽 "
                         f"{len(exact[s])} 个分子 —— {hint}")
                    verdict = f"✗ 错位（MCM 该槽 {len(exact[s])} 分子）"
                else:
                    verdict = f"✓ MCM 槽（该槽 {len(exact[s])} 分子）"
                    if borrowed.get(s):
                        verdict += f"，另有 {len(borrowed[s])} 个分子借该槽形状"
            if not args.quiet:
                print(f"J{s:<3} {' + '.join(sorted(lhs)):22s} "
                      f"{' + '.join(sorted(mine)):22s} {verdict}")

    for s, mol in sorted(CUSTOM.items()):
        if s not in by_slot:
            notes.append(f"自定义槽白名单 J{s}={mol} 未被机制引用（陈旧条目）")

    # ══ ③ JMAP 无死参数 / 不覆盖 MCM 官方槽 ══
    dead = sorted(set(jmap) - set(by_slot))
    if dead:
        fail(f"photolysis.config JMAP 有未被机制引用的死槽位 J{dead} —— "
             f"删掉或补对应光解反应")
    override = sorted(set(jmap) - set(CUSTOM))
    if override:
        fail(f"JMAP 覆盖了 MCM 官方槽 J{override} —— 官方 l/m/n 被手写值悄悄替换"
             f"（V2.5 前 J13 就是这样偏快 4.4×）；JMAP 只允许给 CUSTOM 白名单里的槽")

    # ══ ④ README 光解行与 .fac / JMAP 一致 ══
    rows = scan_readme(README, known)
    doc_pairs = {(ds, proxy(parse_rxn(dtxt)[0])) for _ln, ds, dtxt, _src in rows}
    for s, _t, lhs, _rhs in sorted(phot):
        if (s, proxy(lhs)) not in doc_pairs:
            fail(f".fac 的 J{s} `{' + '.join(lhs)}` 在 README 光解表里没有对应行"
                 f"（§5 反应表或 §6 光解方案表被删/被改槽号）")
    for ln, s, text, src in rows:
        lhs, rhs = parse_rxn(text)
        if s not in by_slot:
            where = sorted(k for k, v in by_slot.items()
                           if any(proxy(l) == proxy(lhs) for l, _ in v))
            hint = f"—— .fac 里该反应物走 J{where}" if where else "—— .fac 无此光解"
            fail(f"README:{ln} 声明 J{s} `{text}`，但 .fac 未使用该槽位 {hint}")
            continue
        cands = by_slot[s]
        hit = next((r for l, r in cands if proxy(l) == proxy(lhs)), None)
        if hit is None:
            actual = " | ".join(" + ".join(sorted(proxy(l))) for l, _ in cands)
            fail(f"README:{ln} J{s} 反应物 {' + '.join(lhs)} 与 .fac 不符（实际 {actual}）")
            continue
        if frozenset(hit) != frozenset(rhs):
            fail(f"README:{ln} J{s} 产物 {' + '.join(rhs)} 与 .fac 不符"
                 f"（实际 {' + '.join(sorted(set(hit)))}）")
        # 参数来源列: §6 惯例以 "MCM 表…" / "JMAP 直达…" 开头（括注里可能提到
        # 旧 JMAP 值，故只看开头，不做子串匹配）
        if src and src.startswith("JMAP") != (s in jmap):
            fail(f"README:{ln} J{s} 参数来源写作 "
                 f"{'JMAP 直达' if src.startswith('JMAP') else 'MCM 表'}，"
                 f"实际{'在' if s in jmap else '不在'} JMAP")

    if not args.quiet:
        print("-" * 78)
        print(f"README 光解表行扫描: {len(rows)} 行（§5 反应表 + §6 光解方案表）")
    # 失败原因与陈旧条目在 --quiet 下也要打印：只留结论行不够定位问题
    for n in notes:
        print(f"  ⚠ {n}")
    if fails:
        print(f"\nFAIL {len(fails)} 处:")
        for f in fails:
            print(f"  ✗ {f}")

    ok = not fails
    n_mcm = len([s for s in by_slot if s not in CUSTOM])
    n_cus = len([s for s in by_slot if s in CUSTOM])
    print(f"\n[{'PASS' if ok else 'FAIL'}] 光解槽位语义 "
          f"({len(phot)} 条光解 / {len(by_slot)} 槽; "
          f"MCM 官方槽 {n_mcm} 个逐槽核对, 自定义槽 {n_cus} 个白名单封闭, "
          f"README {len(rows)} 行"
          + ("一致)" if ok else f", {len(fails)} 处不符)"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
