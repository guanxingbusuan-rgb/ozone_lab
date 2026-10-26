#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
functional_test.py —— ozone_lab 工程级收官功能测试 ("""
import os
import sys
import subprocess
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)
sys.path.insert(0, "diagnostics")
import shutil
import tempfile
import diagnostics as dg

PASS = []
def check(name, ok, detail=""):
    PASS.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")

PY = sys.executable          # 用当前解释器（脱敏 + 更可移植）
PPB = 2.46e10
FS = np.array([0.3, 0.45, 0.67, 1.0, 1.5, 2.0, 2.5, 2.8, 3.0])

# ================= 编译期负向用例 =================
# 用完整的 input_v3 做底（保证能走到物种校验那一步），只注入一个假物种。
# 旧写法只复制一个 .fac 到空目录 → cfgbuild 先因缺 photolysis.config 崩溃，
# returncode 非 0 就"PASS"，其实从未校验过物种 —— 假绿。现改为断言错误信息本身。
bad_dir = tempfile.mkdtemp()
for f in os.listdir("input_v3"):
    p = os.path.join("input_v3", f)
    if os.path.isfile(p):
        shutil.copy(p, bad_dir)
shutil.copytree("input_v3/constraints", os.path.join(bad_dir, "constraints"))
open(os.path.join(bad_dir, "initialConcentrations.config"), "a").write("XXX_UNKNOWN 1.0e9\n")
r = subprocess.run([PY, "src/cfgbuild.py", bad_dir, "/tmp/bad_rt.dat", "--mech",
                    "mechanism/mechanism_3_full_chain.fac"],
                   capture_output=True, text=True)
_out = (r.stdout or "") + (r.stderr or "")
check("未知物种被拒绝 (编译期)", r.returncode != 0 and "XXX_UNKNOWN" in _out,
      "错误信息含该物种名" if "XXX_UNKNOWN" in _out else "未命中物种校验 (可能假绿)")
shutil.rmtree(bad_dir, ignore_errors=True)

# 反向守卫: 气象序列若误写浓度单位指令, 必须报错而非静默 ×2.46e10
bad_dir = tempfile.mkdtemp()
for f in os.listdir("input_v3"):
    p = os.path.join("input_v3", f)
    if os.path.isfile(p):
        shutil.copy(p, bad_dir)
shutil.copytree("input_v3/constraints", os.path.join(bad_dir, "constraints"))
with open(os.path.join(bad_dir, "meteo.config"), "w") as _f:
    _f.write("TEMP = file:constraints/met_TEMP_demo.txt\n"
             "BLH  = file:constraints/met_BLH_demo.txt\n")
_metp = os.path.join(bad_dir, "constraints/met_TEMP_demo.txt")
open(_metp, "w").write("# unit: ppb\n" + open(_metp).read())
r = subprocess.run([PY, "src/cfgbuild.py", bad_dir, "/tmp/bad_rt2.dat", "--mech",
                    "mechanism/mechanism_3_full_chain.fac"],
                   capture_output=True, text=True)
_out = (r.stdout or "") + (r.stderr or "")
check("气象序列拒绝浓度单位指令 (防静默 ×2.46e10)", r.returncode != 0 and "气象序列" in _out,
      "cfgbuild 拒绝并说明原因")
shutil.rmtree(bad_dir, ignore_errors=True)


# ---------- FAC-BATCH-OMP 指示因子批量自检 (6×6) ----------
import subprocess, csv, io
HERE2 = os.path.dirname(os.path.abspath(__file__))
# FT0 假绿防线: 生成码物种数必须与 .fac VARIABLE 同步 (陈旧二进制 → 此处硬失败)
import re as _re
_mech = open(f"{HERE2}/../mechanism/mechanism_3_full_chain.fac").read()
_n_sp = len(set(_re.search(r"VARIABLE(.*?)FORMULA", _mech, _re.S).group(1)
                .replace(";", " ").split()))
_gen_src = open(f"{HERE2}/../facsimile_impl/gen_v3/oz_mech.f90").read()
_ns = int(_re.search(r"::\s*ns=(\d+)", _gen_src).group(1))
_nr_gen = int(_re.search(r"\bnr=(\d+)", _gen_src).group(1))  # ns/nr 同行, nr 前无 ::
_n_rx = len([l for l in _mech.split("\n") if l.strip().startswith("%")])
print(f"== FT0 生成码一致性: gen ns={_ns}/nr={_nr_gen} vs 机制物种 {_n_sp}/反应 {_n_rx} ==")
assert _ns == _n_sp, f"gen_v3 陈旧! ns={_ns} != 机制物种 {_n_sp} — 先 make v3"
assert _nr_gen == _n_rx, f"gen_v3 陈旧! nr={_nr_gen} != 机制反应 {_n_rx} — 先 make v3"
print("   [PASS] 生成码与机制同步 (物种+反应双守卫)")

# R10-P0 情景缩放完备性: HVOC 必须放大全部 E_VOC 组物种 (含 ALK/OLE/ARO/KET/TERP/SO2)
import subprocess as _sp, csv as _csv
def _eacc(pref, sp):
    rows = list(_csv.reader(open(f"{pref}_ipr.csv"))); hd = rows[0]; j = hd.index(f"{sp}_E")
    return sum(float(r[j]) for r in rows[1:] if r[j].strip())
_sp.run([f"{HERE2}/../facsimile_impl/fac_v3", "SCEN_BAU", "/tmp/ft_p0_bau",
         f"{HERE2}/../runtime_v3.dat"], capture_output=True, cwd=f"{HERE2}/../facsimile_impl")
_sp.run([f"{HERE2}/../facsimile_impl/fac_v3", "SCEN_HVOC", "/tmp/ft_p0_hvoc",
         f"{HERE2}/../runtime_v3.dat"], capture_output=True, cwd=f"{HERE2}/../facsimile_impl")
_bad = []
for _sp_name in ("OLE", "ALK", "ARO1", "ARO2", "KET", "TERP", "SO2", "ISOP", "HCHO"):
    _a = _eacc("/tmp/ft_p0_bau", _sp_name); _b = _eacc("/tmp/ft_p0_hvoc", _sp_name)
    if not (_a > 0 and 2.7 < _b/_a < 3.3):
        _bad.append(f"{_sp_name}:{_b/_a:.2f}")
check("R10-P0 HVOC 缩放全部 E_VOC 组物种 (×3)", not _bad,
      f"未缩放/异常: {_bad}" if _bad else "9 物种全部 ≈×3.03")

# R11-P0 等价路径防线: cfg JFAC=0.8 + CLI 1.0  应等于  cfg 1.0 + CLI 0.8 (容差 1%)
# (旧 bug: eval_J_jf 与 eval_J 二次相乘 → 差 8.07%)
import tempfile as _tf, shutil as _sh
_eq = _tf.mkdtemp(prefix="r11eq_")
for _f in ("emissions.config", "initialConcentrations.config", "meteo.config",
           "photolysis.config", "photolysis-rates_MCM", "bc.config", "drydep.config",
           "wetdep.config", "dilution.config", "scenarios.config", "model.parameters"):
    _sh.copy(f"{HERE2}/../input_v3/{_f}", f"{_eq}/{_f}")
_sh.copytree(f"{HERE2}/../input_v3/constraints", f"{_eq}/constraints")
_m = open(f"{_eq}/meteo.config").read().replace("JFAC   = 1.0", "JFAC   = 0.8")
open(f"{_eq}/meteo.config", "w").write(_m)
subprocess.run([PY, f"{HERE2}/../src/cfgbuild.py", _eq, f"{_eq}/rt.dat",
                "--mech", f"{HERE2}/../mechanism/mechanism_3_full_chain.fac"],
               capture_output=True)
def _o3max_csv(p):
    rows = list(csv.reader(open(p)))
    for r in rows[1:]:
        if abs(float(r[0]) - 1.0125) < 0.01 and abs(float(r[1]) - 1.0125) < 0.01:
            return float(r[4])
    return None
for _tag, _rt, _jf in (("A", f"{HERE2}/../runtime_v3.dat", "0.8"),
                       ("B", f"{_eq}/rt.dat", "1.0")):
    subprocess.run([f"{HERE2}/../facsimile_impl/fac_batch_omp", "4",
                    f"{_eq}/{_tag}.csv", _rt, _jf, "0"],
                   capture_output=True, cwd=f"{HERE2}/../facsimile_impl")
_a, _b = _o3max_csv(f"{_eq}/A.csv"), _o3max_csv(f"{_eq}/B.csv")
if _a and _b:
    _d = abs(_a - _b)/_a*100
    check("R11-P0 jfac 等价路径 (cfg0.8×CLI1.0 == cfg1.0×CLI0.8)", _d < 1.0,
          f"A={_a:.2f} B={_b:.2f} 差 {_d:.3f}% (旧 bug 为 8.07%)")
else:
    check("R11-P0 jfac 等价路径", False, "无法取得 O3max")
_sh.rmtree(_eq, ignore_errors=True)

# R12 数据链防线: Makefile 的 batch 输出名 必须 == figures 依赖名
_mk = open(f"{HERE2}/../Makefile").read()
import re as _re2
_batch_out = _re2.search(r"fac_batch_omp 9 \.\./(output/\S+?\.csv)", _mk)
_fig_dep = _re2.search(r"^figures: (\S+)", _mk, _re2.M)
if _batch_out and _fig_dep:
    _bo, _fd = _batch_out.group(1), _fig_dep.group(1)
    check("R12 batch 输出名 == figures 依赖名 (数据链不断)", _bo == _fd,
          f"batch→{_bo}  figures←{_fd}" + (" (数据链断裂!)" if _bo != _fd else ""))
else:
    check("R12 batch/figures 数据链", False, "无法解析 Makefile")

# R13 数值格式防线: 极痕量物种必须带 E (ES15.5E3 强制 3 位指数)
import re as _re3, glob as _glob
_bad_fmt = 0; _worst = None
for _csvf in _glob.glob("/tmp/ft_p0_bau_conc.csv"):
    _txt = open(_csvf).read()
    for _m in _re3.finditer(r"(-?[0-9]\.[0-9]+)([-+])([0-9]{3})\b", _txt):
        _bad_fmt += 1
        _worst = _m.group(0)
check("R13 输出格式无三位无E指数 (ES15.5E3)", _bad_fmt == 0,
      f"发现 {_bad_fmt} 处" + (f", 如 {_worst}" if _worst else ""))

# R13 解析健壮性: csvfloat 必须能吃四种形态
import sys as _sys3
_sys3.path.insert(0, f"{HERE2}/../diagnostics")
from csvfloat import fix_fortran_exp as _fx3
_cases = [("3.07932-174", 3.07932e-174), ("-9.75647-278", -9.75647e-278),
          ("5.00000E-99", 5.0e-99), ("1.00000+200", 1.0e200)]
_ok3 = all(abs(float(_fx3(a)) - b) <= abs(b)*1e-15 for a, b in _cases)
check("R13 csvfloat 四种指数形态解析", _ok3, f"{len(_cases)} 组")

# R16 配置项消费防线: MODE=DIRECT 与 MCM 必须产生不同结果 (防死开关)
import tempfile as _tf16, shutil as _sh16
_eq16 = _tf16.mkdtemp(prefix="r16mode_")
for _f in ("emissions.config","initialConcentrations.config","meteo.config",
           "photolysis.config","photolysis-rates_MCM","bc.config","drydep.config",
           "wetdep.config","dilution.config","scenarios.config","model.parameters"):
    _sh16.copy(f"{HERE2}/../input_v3/{_f}", f"{_eq16}/{_f}")
_sh16.copytree(f"{HERE2}/../input_v3/constraints", f"{_eq16}/constraints")
_pc = open(f"{_eq16}/photolysis.config").read().replace("MODE = MCM", "MODE = DIRECT")
open(f"{_eq16}/photolysis.config", "w").write(_pc)
subprocess.run(["python3", f"{HERE2}/../src/cfgbuild.py", _eq16, f"{_eq16}/rt.dat",
                "--mech", f"{HERE2}/../mechanism/mechanism_3_full_chain.fac"],
               capture_output=True)
def _o3_of(pref):
    rr = list(csv.reader(open(f"{pref}_conc.csv"))); hh = [x.strip() for x in rr[0]]
    jj = hh.index("O3")
    def _ff(v):
        v = v.strip()
        if "E" in v or "e" in v: return float(v)
        m = _re3.match(r"^(-?[0-9.]+)-([0-9]{3})$", v)
        return float(v[:-4]+"E-"+v[-3:]) if m else float(v)
    return max(_ff(r[jj]) for r in rr[1:] if r[jj].strip())/2.46e10
for _tag, _rt in (("MCM", f"{HERE2}/../runtime_v3.dat"), ("DIRECT", f"{_eq16}/rt.dat")):
    subprocess.run([f"{HERE2}/../facsimile_impl/fac_v3", "SCEN_BAU",
                    f"{_eq16}/{_tag}", _rt], capture_output=True,
                   cwd=f"{HERE2}/../facsimile_impl")
_m, _d = _o3_of(f"{_eq16}/MCM"), _o3_of(f"{_eq16}/DIRECT")
check("R16 MODE=DIRECT 与 MCM 产生不同结果 (非死开关)", abs(_m-_d)/_m > 0.005,
      f"MCM={_m:.2f} DIRECT={_d:.2f} 差 {abs(_m-_d)/_m*100:.2f}%")
_sh16.rmtree(_eq16, ignore_errors=True)

r = subprocess.run([f"{HERE2}/../facsimile_impl/fac_batch_omp", "9",
                    "/tmp/ft_batch_ind.csv", f"{HERE2}/../runtime_v3.dat"],
                   capture_output=True, text=True, cwd=f"{HERE2}/../facsimile_impl")
rows = list(csv.reader(open("/tmp/ft_batch_ind.csv")))
# R10 尾逗号 bug 防线: ipr/irr 侧文件的表头列数必须等于数据列数
for _side in ("_ipr", "_irr"):
    _p = f"/tmp/ft_batch_ind.csv{_side}.csv"
    try:
        _l = open(_p).read().split("\n")
        _nh = len(_l[0].split(",")); _nd = len(_l[1].split(","))
        check(f"R10 batch{_side} 表头/数据列数一致", _nh == _nd,
              f"{_nh} vs {_nd}" + (" (尾逗号空列 bug)" if _nh != _nd else ""))
    except FileNotFoundError:
        pass
hdr, data = rows[0], [r for r in rows[1:] if r]
ok = len(hdr) == 24 and len(data) == 81
check("FAC-BATCH 24 列 × 81 案例", ok, f"列={len(hdr)} 行={len(data)}")
d = {(round(float(r[0]),2), round(float(r[1]),2)): r for r in data}
fnox_l = sorted({float(r[0]) for r in data})
j1 = min((k[1] for k in d), key=lambda v: abs(v-1.0))
col = sorted([d[(round(fn,2), j1)] for fn in fnox_l], key=lambda r: float(r[0]))
o3 = [float(r[4]) for r in col]; fnr = [float(r[8]) for r in col]
lq = [float(r[14]) for r in col]
ir = int(np.argmax(o3))
drop = 1 - min(o3[ir:]) / max(o3[ir:])          # 峰后最大回撤
check("屋脊存在 (先升后峰后回撤≥20%)",
     max(o3) > o3[0]*1.5 and ir < len(o3)-1 and drop >= 0.20,
     f"O3 {o3[0]:.1f}→{max(o3):.1f}@{ir}→{o3[-1]:.1f} 回撤{drop*100:.0f}%")
check("FNR 随 NOx 单调下降", fnr[-1] < fnr[0]*0.3, f"{fnr[0]:.2f} → {fnr[-1]:.2f}")
# ---------- 回归哨兵: 基准格点金值 (抓数值漂移, 如 OPE 双计类 bug) ----------
SENT = [  # (fnox, 列, 期望, 容差) — 列为 0 基 CSV 索引 (24 列格式); V2.2 哨兵
    ("2.28", 4, 132.10, 0.10),   # O3max 屋脊峰 (燃料补齐后上移 4.4×)
    ("3.42", 4, 93.68, 0.10),    # O3max 饱和区
    ("1.01", 8, 3.539, 0.15),    # FNR 转折带
    ("2.28", 16, 6.487, 0.15),   # OPE 饱和区谷 (双计 bug 的教训位)
]
for fn, c, exp, tol in SENT:
    row = next(r for r in data if abs(float(r[0]) - float(fn)) < 0.01 and abs(float(r[1]) - 1.0125) < 0.01)
    v = float(row[c])
    check(f"哨兵 fnox={fn} 列{c}={exp}", abs(v - exp) <= tol*exp, f"{v:.3f} ∈ {exp}±{tol:.0%}")
# ---------- IPR 分解完备性 (物理层补偿检验: 六账积分 vs ΔC) ----------
import subprocess as _sp
_r = _sp.run([f"{HERE2}/../facsimile_impl/fac_v3", "SCEN_BAU", "/tmp/ft_ipr",
              f"{HERE2}/../runtime_v3.dat"], capture_output=True,
             cwd=f"{HERE2}/../facsimile_impl")
_rows = list(csv.reader(open("/tmp/ft_ipr_ipr.csv")))
_hdr2 = [h.strip() for h in _rows[0]]
_D = np.array([[float(dg._fix_fortran_exp(x)) for x in r] for r in _rows[1:] if r])
_crows = list(csv.reader(open("/tmp/ft_ipr_conc.csv")))
_ch = [c.strip() for c in _crows[0]]
_C = np.array([[float(dg._fix_fortran_exp(x)) for x in r] for r in _crows[1:] if r])
_ratios = []
for _sp in ("O3", "HNO3", "H2O2", "CO"):
    _ic = _ch.index(_sp)      # _C 含 h 列 (未剥), 无需 -1
    _cols = [_hdr2.index(f"{_sp}_{p}") for p in "PLEDBX"]
    _six = _D[:, _cols].sum(axis=1)
    _npair = min(len(_six), _C.shape[0]) - 1
    _dC = np.full(len(_six), np.nan)
    _dC[1:1+_npair] = _C[1:1+_npair, _ic] - _C[:_npair, _ic]
    _big = np.abs(_six) > 1e7
    for _k in range(len(_six)):
        if _big[_k] and not np.isnan(_dC[_k]) and abs(_six[_k]) > 1e-30:
            _ratios.append(_dC[_k] / _six[_k])   # 逐窗比值池 (跳过 nan/零信量)
_med = float(np.median(_ratios)) if _ratios else float("nan")
check("IPR 分解完备性 (六账积分 vs ΔC, 比值中位 ∈ [0.9,1.1])",
      0.9 <= _med <= 1.1, f"中位比值 {_med:.3f} ({len(_ratios)} 样本=4 物种×时窗)")

check("L_NOx/Q ∈ [0,1.5] 全格点 (NOx 饱和区合法 >1)", all(0 <= x <= 1.5 for x in lq),
     f"{min(lq):.3f}–{max(lq):.3f}")

# ================= FT-SAT: sat2surf 卫星转换器 (量纲 bug 防线) =================
print("== FT-SAT sat2surf 卫星柱→表面 ppb 转换 ==")
import importlib.util as _ilu
_spec = _ilu.spec_from_file_location("sat2surf", os.path.join(HERE, "src", "sat2surf.py"))
_s2s = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_s2s)

# ① 内置自检 (3 组已知值 + 两法互证 + 旧式反证)
_r = subprocess.run([sys.executable, "src/sat2surf.py", "--selftest"],
                    capture_output=True, text=True)
check("sat2surf 内置自检 (3 组已知值/两法互证/旧式反证)", _r.returncode == 0,
      _r.stdout.strip().splitlines()[-1] if _r.stdout else "")

# ② 量纲防线: 修正式与直观式必须一致到机器精度; 旧式(漏 1e9)必须差 9 量级
_b, _a = _s2s.convert(3.8e16, 1000.0)
_rel = abs(_b - _a) / _a
_old = 3.8e16 / (_s2s.AIR_COL * (1e5 / _s2s.AIR_SCALE))
_decades = abs(np.log10(_old * 1e9 / _old))     # 应为 9
check("sat2surf 默认两法等价 (相对差 <1e-12) 且旧式确实差 9 量级",
      _rel < 1e-12 and abs(_decades - 9) < 1e-6,
      f"两法差 {_rel:.1e}, 旧式漏 1e9 = {_decades:.0f} 个量级, 3.8e16→{_b:.2f} ppb")


_tmp = tempfile.mkdtemp()
# ③ 端到端: 造原始柱 → 转换 → bc 登记 → cfgbuild → SERIES 值 = ppb×2.46e10
try:
    _raw = os.path.join(_tmp, "demo.raw")
    open(_raw, "w").write(
        "# species: NO2_tropo\n"
        "# t  column  qa\n"
        "2024-06-01T13:30  4.2e15  0.92\n"
        "2024-06-02T13:30  3.8e16  0.88\n"
        "2024-06-03T13:30  5.1e16  0.40\n")     # 第 3 条 QA<0.75 应被剔除
    _surf = os.path.join(_tmp, "surf.txt")
    _r = subprocess.run([sys.executable, "src/sat2surf.py", _raw, _surf],
                        capture_output=True, text=True)
    _lines = [l for l in open(_surf) if l.strip() and not l.startswith("#")]
    _ok_tool = _r.returncode == 0 and len(_lines) == 2      # QA 剔除后 2 点
    _vals = [float(l.split()[1]) for l in _lines]
    _ok_val = abs(_vals[0] - 1.71) < 0.02 and abs(_vals[1] - 15.45) < 0.05

    # 接入 cfgbuild
    _inp = os.path.join(_tmp, "input")
    shutil.copytree("input_v3", _inp)
    shutil.copy(_surf, os.path.join(_inp, "constraints", "sat_NO2.txt"))
    open(os.path.join(_inp, "bc.config"), "w").write(
        "TAU_MIX = 21600\nO3 = 6.2e11\nCO = 1.5e12\nCH4 = 4.5e13\n"
        "NO2 = file:constraints/sat_NO2.txt\n")
    _rt = os.path.join(_tmp, "rt.dat")
    _r2 = subprocess.run([sys.executable, "src/cfgbuild.py", _inp, _rt,
                          "--mech", "mechanism/mechanism_3_full_chain.fac"],
                         capture_output=True, text=True)
    _txt = open(_rt).read() if os.path.exists(_rt) else ""
    _ok_cfg = _r2.returncode == 0 and "bc_NO2" in _txt
    # SERIES 值应 = ppb × 2.46e10
    _ser = []
    _inf = False
    for l in _txt.splitlines():
        if l.strip().startswith("bc_NO2"):
            _inf = True; continue
        if _inf and l.strip() and l.strip()[0].isdigit():
            _ser.append(l.split())
            if len(_ser) >= 2: break
    _fit = False
    if len(_ser) >= 2 and len(_ser[1]) >= 2:
        _v1 = float(_ser[1][0]); _fit = abs(_v1 / (1.71 * 2.46e10) - 1) < 0.01
    check("sat2surf 端到端 (QA剔除/数值正确/cfgbuild SERIES 单位换算)",
          _ok_tool and _ok_val and _ok_cfg and _fit,
          f"2 点 {_vals}, SERIES 首值 {_ser[1][0] if len(_ser)>1 else 'N/A'}")
    # ④ HCHO 产品同样支持 (无白名单): 转换 + 量级合理
    _raw_h = os.path.join(_tmp, "hcho.raw")
    open(_raw_h, "w").write(
        "# species: HCHO_tropo\n"
        "# t  column  qa\n"
        "2024-06-01T13:30  8.5e15  0.88\n"
        "2024-06-02T13:30  1.8e16  0.86\n")
    _surf_h = os.path.join(_tmp, "hcho_surf.txt")
    _r3 = subprocess.run([sys.executable, "src/sat2surf.py", _raw_h, _surf_h],
                         capture_output=True, text=True)
    _lh = [l for l in open(_surf_h) if l.strip() and not l.startswith("#")]
    _vh = [float(l.split()[1]) for l in _lh]
    # HCHO 典型 1-10 ppb; 且 species 名应正确识别为 HCHO
    _ok_h = (_r3.returncode == 0 and len(_vh) == 2
             and all(1.0 <= v <= 10.0 for v in _vh)
             and "HCHO" in open(_surf_h).read())
    check("sat2surf 支持 HCHO 产品 (无物种白名单 + 量级 1-10 ppb)",
          _ok_h, f"HCHO 转换 {['%.2f' % v for v in _vh]} ppb")

    # ⑤ 分歧策略: 显式给不同 scale-height → 仅警告不中断, 且输出仍用直观法
    _demo_raw = os.path.join(_tmp, "demo.raw")
    with open(_demo_raw, "w") as _f:
        _f.write("# t column[molec/cm2] qa\n")
        _f.write("2026-07-01T12:00 3.8e16 0.95\n")
        _f.write("2026-07-01T13:00 7.6e16 0.95\n")
    _dv = os.path.join(_tmp, "diverge.txt")
    _r4 = subprocess.run([sys.executable, "src/sat2surf.py", _demo_raw,
                          _dv, "--scale-height-km", "2.0"],
                         capture_output=True, text=True)
    _dv_txt = open(_dv).read() if os.path.exists(_dv) else ""
    _dv_vals = [float(l.split()[1]) for l in _dv_txt.splitlines()
                if l.strip() and not l.startswith("#")]
    _ok_dv = (_r4.returncode == 0                       # 不中断
              and "警告" in _r4.stderr                   # 有警告
              and "直观法" in _r4.stderr                 # 指明主准则
              and len(_dv_vals) == 2                     # 仍产出
              and abs(_dv_vals[0] - 15.45) < 0.1         # 用直观法值(参考法会高 2.19 倍)
              and "已按直观法输出" in _dv_txt)           # 文件留痕
    check("sat2surf 两法分歧时仅警告且按直观法继续 (exit 0 + 产物)",
          _ok_dv,
          f"exit={_r4.returncode}, 值={['%.2f' % v for v in _dv_vals]}, "
          f"stderr={'有警告' if '警告' in _r4.stderr else '无'}")

    # ⑥ 选择性输入: 极简配置 (IC 3 个/排放 2 个/沉降 2 个/BC 1 个) 必须能跑
    _mn = os.path.join(_tmp, "minimal")
    shutil.copytree("input_v3", _mn)
    open(os.path.join(_mn, "initialConcentrations.config"), "w").write(
        "O3 3.0e11 cm3\nNO 5.0e8 cm3\nCO 2.0e12 cm3\n")
    open(os.path.join(_mn, "emissions.config"), "w").write(
        "NO 2.0e7 BOX 0.9\nISOP 2.0e6 BOX 0.9\n")
    open(os.path.join(_mn, "drydep.config"), "w").write(
        "O3 0.8 on\nHNO3 2.0 on\n")
    open(os.path.join(_mn, "wetdep.config"), "w").write("")
    open(os.path.join(_mn, "bc.config"), "w").write("TAU_MIX = 86400\nO3 = 6.2e11\n")
    open(os.path.join(_mn, "meteo.config"), "w").write(
        "TEMP = 298.15\nPRES = 1013.25\nRH = 60\nBLH = 1000.0\n"
        "RAIN = 0.0\nJFAC = 1.0\n")
    _mrt = os.path.join(_tmp, "minimal_rt.dat")
    _r5 = subprocess.run([sys.executable, "src/cfgbuild.py", _mn, _mrt,
                          "--mech", "mechanism/mechanism_3_full_chain.fac"],
                         capture_output=True, text=True)
    _r6 = subprocess.run([os.path.join(HERE, "facsimile_impl", "fac_v3"),
                          "SCEN_BAU", os.path.join(_tmp, "minimal_run"), _mrt],
                         capture_output=True, text=True,
                         cwd=os.path.join(HERE, "facsimile_impl"))
    check("选择性输入: 极简配置 (IC3/排放2/沉降2/BC1) 可跑通",
          _r5.returncode == 0 and _r6.returncode == 0 and "FAC-V3 ok" in _r6.stdout,
          f"cfgbuild rc={_r5.returncode}, fac_v3 rc={_r6.returncode}")

    # ⑦ 选择性输入: 全空配置 (全部清空) 也必须合法
    _em = os.path.join(_tmp, "empty")
    shutil.copytree(_mn, _em)
    for _f in ("initialConcentrations.config", "emissions.config",
               "drydep.config", "wetdep.config"):
        open(os.path.join(_em, _f), "w").write("")
    open(os.path.join(_em, "bc.config"), "w").write("TAU_MIX = 86400\n")
    _ert = os.path.join(_tmp, "empty_rt.dat")
    _r7 = subprocess.run([sys.executable, "src/cfgbuild.py", _em, _ert,
                          "--mech", "mechanism/mechanism_3_full_chain.fac"],
                         capture_output=True, text=True)
    _r8 = subprocess.run([os.path.join(HERE, "facsimile_impl", "fac_v3"),
                          "SCEN_BAU", os.path.join(_tmp, "empty_run"), _ert],
                         capture_output=True, text=True,
                         cwd=os.path.join(HERE, "facsimile_impl"))
    check("选择性输入: 全空配置 (IC/排放/沉降全清空) 仍合法可跑",
          _r7.returncode == 0 and _r8.returncode == 0 and "FAC-V3 ok" in _r8.stdout,
          f"cfgbuild rc={_r7.returncode}, fac_v3 rc={_r8.returncode}")

finally:
    shutil.rmtree(_tmp, ignore_errors=True)

# ================= FT-CFG: runtime.dat 与 input_v3 一致性 (mtime 陷阱防线) =================
print("== FT-CFG runtime_v3.dat 与配置一致性 ==")
# 陷阱: 同一秒内改配置又改回 → make 误判 up-to-date → runtime.dat 保留中间值
# 防线: 用 cfgbuild 重新生成一份到临时文件, 与仓库内 runtime_v3.dat 逐字节比对
_tmp_rt = os.path.join(tempfile.mkdtemp(), "rt_check.dat")
try:
    _rc = subprocess.run([sys.executable, "src/cfgbuild.py", "input_v3", _tmp_rt,
                          "--mech", "mechanism/mechanism_3_full_chain.fac"],
                         capture_output=True, text=True)
    _fresh = open(_tmp_rt).read() if os.path.exists(_tmp_rt) else ""
    _repo = open("runtime_v3.dat").read() if os.path.exists("runtime_v3.dat") else ""
    check("runtime_v3.dat 与 input_v3 当前内容一致 (非陈旧)",
          _rc.returncode == 0 and _fresh == _repo,
          "逐字节一致" if _fresh == _repo else "⚠ 仓库内 runtime_v3.dat 陈旧, 请 make runtime_v3.dat")
finally:
    shutil.rmtree(os.path.dirname(_tmp_rt), ignore_errors=True)

n_pass = sum(1 for x in PASS if x)
print(f"\n===== functional_test: {n_pass}/{len(PASS)} PASS =====")
sys.exit(0 if n_pass == len(PASS) else 1)
