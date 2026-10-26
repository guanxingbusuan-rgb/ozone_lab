#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mechc_c.py —— #5 优化: 机制 → C 代码生成 (绕过 Fortran/f2py 链)
  mechc.py 的 C 版: .fac → C 速率函数 → gcc -shared → Python ctypes 直调
  消除 f2py 三步构建开销 + Fortran-C ABI 转换层
用法: python mechc_c.py <fac> <out.c> [--compile]
"""
import re, sys, subprocess, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "facsimile_impl"))
from mechc import parse_fac, to_f, RESERVED

def to_c(expr, m, jn=()):
    """FACSIMILE → C: @→pow, D→E, EXP→exp, 命名速率→q[i], J<n>→J[n-1]"""
    expr = expr.split("!", 1)[0].strip()   # 剥离行内 ! 注释
    rate_names = [n for n, _ in m.rates]
    def repl(mm):
        tok = mm.group(0)
        if tok in rate_names:
            return f"q[{rate_names.index(tok)}]"
        if tok in ("TEMP", "M", "O2", "H2O"):
            return tok
        if tok == "EXP":
            return "exp"
        if tok == "LOG10":
            return "log10"
        return tok
    # 先处理 J<n>
    e = re.sub(r"J<([0-9]+)>", lambda g: "J[" + str(int(g.group(1))-1) + "]", expr)
    # FACSIMILE 的 10@x 是"10 的 x 次方", 必须先于 @→** 统一转成 pow(10,x)
    # (否则会变成 10**(x), C 的 ** 不是运算符 → 编译失败; 主机制 Troe 全用此写法)
    e = re.sub(r"10@\(([^)]*)\)", r"pow(10,(\1))", e)
    # @ → ** 再转 pow
    e = e.replace("@", "**")
    # D 记法 → E (含无符号指数: 1.90D17 → 1.90E17)
    e = re.sub(r"(\d)[Dd]([+-]?\d)", r"\1E\2", e)
    # 标识符替换
    e = re.sub(r"[A-Za-z_][A-Za-z0-9_]*", repl, e)
    return _pow_to_c(e)


def _pow_to_c(e):
    r"""把 Python 风格 a**b 转成 C 的 pow(a,b)。

    必须做"平衡括号反向扫描"而不是正则 —— 因为左操作数可能是
    `(...)` 形式（Troe 的 `(log10(q[3])/(1+(...)**2))`），
    正则 `\w+\*\*` 只能吃掉标识符，漏掉 `)**2` 导致生成非法 C。
    """
    while "**" in e:
        i = e.index("**")
        # 反向扫描左操作数
        j = i
        depth = 0
        while j > 0:
            ch = e[j - 1]
            if ch == ")":
                depth += 1
            elif ch == "(":
                if depth == 0:
                    break
                depth -= 1
            elif depth == 0 and not (ch.isalnum() or ch in "._[]"):
                break
            j -= 1
        left = e[j:i]
        # 正向扫描右操作数
        k = i + 2
        depth = 0
        while k < len(e):
            ch = e[k]
            if ch == "(":
                depth += 1
            elif ch == ")":
                if depth == 0:
                    break
                depth -= 1
            elif depth == 0 and not (ch.isalnum() or ch in "._[]+-"):
                break
            k += 1
        right = e[i + 2:k]
        e = e[:j] + f"pow({left},{right})" + e[k:]
    return e

def gen_c(fac_path, out_c):
    m = parse_fac(fac_path)
    ns, nr, nq = len(m.species), len(m.reactions), len(m.rates)
    rate_names = [n for n, _ in m.rates]

    L = []
    w = L.append
    w('#include <math.h>')
    w(f'// ===== 自动生成: mechc_c.py < {os.path.basename(fac_path)} =====')
    w(f'// {ns} species, {nr} reactions, {nq} named rates')
    w(f'#define NS {ns}')
    w(f'#define NR {nr}')
    w(f'#define NQ {nq}')
    w('')
    w('// 速率更新: 环境参数 + 浓度 → q(命名速率) + p(反应速率)')
    w('void update_rates_c(double TEMP, double M, double O2, double H2O,')
    w(f'                    const double *J, const double *Y, double *q, double *p) {{')
    w('    (void)O2; (void)H2O; (void)Y;  /* 避免未使用告警 */')
    # 命名速率
    for i, (name, expr) in enumerate(zip(rate_names, [e for _, e in m.rates])):
        w(f'    q[{i}] = {to_c(expr, m)};')
    # 反应速率
    for k, (expr, lhs, rhs) in enumerate(m.reactions):
        ce = to_c(expr, m)
        # 引用 q 的地方替换为 q[i-1]
        for i, name in enumerate(rate_names):
            ce = re.sub(rf'\bq\({i+1}\)', f'q[{i}]', ce)
        w(f'    p[{k}] = {ce};')
    w('}')
    w('')
    # resid: 通用稀疏累积 (与 Fortran 同构)
    w('// resid: 排放+沉降+背景松弛+化学 通用累积')
    w(f'void resid_c(int nr, const int *nreact, const int (*reac_id)[2],')
    w(f'             const int *nprod, const int (*prod_id)[3], const double (*prod_coef)[3],')
    w(f'             const double *p_eff, const double *E, const double *kdep,')
    w(f'             const double *bcv, double tau_mix, double *ydot) {{')
    w(f'    for (int i = 0; i < {ns}; i++) ydot[i] = E[i];')
    w(f'    // 化学累积')
    w(f'    for (int k = 0; k < NR; k++) {{')
    w(f'        for (int i = 0; i < nreact[k]; i++)')
    w(f'            ydot[reac_id[k][i] - 1] -= p_eff[k];')
    w(f'        for (int i = 0; i < nprod[k]; i++)')
    w(f'            ydot[prod_id[k][i] - 1] += prod_coef[k][i] * p_eff[k];')
    w(f'    }}')
    w(f'}}')

        # 后处理: ** → pow()
    text = "\n".join(L)
    # (expr)**num → pow(expr, num)
    text = re.sub(r"\(([^)]+)\)\*\*([+-]?[\d.]+)", r"pow(\1,\2)", text)
    # var**num → pow(var, num)
    text = re.sub(r"(\w+)\*\*(\d+)", r"pow(\1,\2)", text)
    text = re.sub(r"(\w+)\*\*([+-]?[\d.]+)", r"pow(\1,\2)", text)
    open(out_c, "w").write(text + "\n")
    print(f"[mechc_c] {ns} species, {nr} reactions → {out_c}")
    return 0

if __name__ == "__main__":
    fac = sys.argv[1] if len(sys.argv) > 1 else "mechanism/mechanism_3_full_chain.fac"
    out = sys.argv[2] if len(sys.argv) > 2 else "output/mech_rates.c"
    os.makedirs(os.path.dirname(out), exist_ok=True)
    gen_c(fac, out)
    # 自动编译
    cc = os.path.expanduser("~/gcc-13.2.0/bin/gcc")
    dylib = out.replace(".c", ".dylib")
    r = subprocess.run([cc, "-O3", "-fPIC", "-shared", out, "-o", dylib, "-lm"],
                       capture_output=True, text=True)
    if r.returncode == 0:
        print(f"[mechc_c] 编译 → {dylib}")
    else:
        print(f"[mechc_c] 编译失败: {r.stderr[:200]}")
