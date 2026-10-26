# ozone_lab 实验操作指南（EXPERIMENTS_GUIDE）

> **定位**：本文是"动手做实验"的**完整教学级操作手册**——每个实验都给出
> **科学问题 → 前置条件 → 逐步命令 → 预期输出（含实测样例）→ 产物逐列解读 →
> 判读要点 → 变异玩法 → 常见坑**八段式。零基础按 §1→§4 顺序读，之后按需跳章节。
>
> **四份文档分工**：
> [README.md](README.md)（项目是什么）· [OPERATIONS_MANUAL.md](OPERATIONS_MANUAL.md)
> （架构与速查）· **本文（实验怎么做，手把手）** ·
> [figures/README.md](figures/README.md)（机制图集与七机制对比）。
>
> **适用版本**：2026-10-26 · mechanism_3 **V2.3（74sp/131rx/41 光解）** ·
> 现役三路线（fac_v3 / fac_batch_omp / mcm_run）·
> **本文所有命令与数字均为当日实测**（可复现）。

---

## 目录

- [0. 实验运行脚本（experiments/run/）](#0-实验运行脚本experimentsrun--先读这一节)
  - [0.3 十个运行脚本](#03-十个运行脚本ai--m-不含伴随反演)
  - [0.3b 气象序列输入（实验 M）](#03b-气象序列输入实验-m-的完整打法)
- [1. 环境与首次构建](#1-环境与首次构建)
- [2. 化学基础（读实验前必读）](#2-化学基础读实验前必读)
- [3. 数据如何进程序（12 类配置详解）](#3-数据如何进程序12-类配置详解)
- [4. 路线选择与编译边界](#4-路线选择与编译边界)
- [5. 实验 A：单箱基线与政策情景](#5-实验-a单箱基线与政策情景)
- [6. 实验 B：响应面与指示因子屋脊](#6-实验-b响应面与指示因子屋脊)
  - [6.9 响应面能力全景（能跑多少种类 EKMA 曲面）](#69-响应面能力全景能跑多少种类-ekma-曲面)
  - [6.9.6 经典响应面清单：能实现 4 款，缺 3 款](#696-经典响应面清单能实现-4-款缺-3-款)
  - [6.9.7 FNR 响应面专项（最核心的一张图）](#697-fnr-响应面专项最核心的一张图)
- [6.10 四大分析能力状态核实](#610-四大分析能力状态核实观测约束--敏感性--伴随--响应面)
- [6.11 情景分析与过程分析](#611-情景分析与过程分析两项都是核心能力)
- [7. 实验 C：屋脊的跨路线对比](#7-实验-c屋脊的跨路线对比)
- [8. 实验 D：IPR/IRR 过程诊断](#8-实验-diprirr-过程诊断)
- [9. 实验 E：指示因子工具箱](#9-实验-e指示因子工具箱)
- [10. 实验 F：观测约束 OSSE](#10-实验-f观测约束-osse)
- [11. 实验 G：逐小时过程分析](#11-实验-g逐小时过程分析)
- [12. 实验 H：灵敏度/伴随分析](#12-实验-h灵敏度伴随分析)
- [13. 实验 I：情景矩阵](#13-实验-i情景矩阵)
- [14. 六类可选输入数据文件接入操作（templates/ 逐一详解）](#14-六类可选输入数据文件接入操作templates-逐一详解)
- [15. 输出契约（含真实样例行）](#15-输出契约含真实样例行)
- [16. 指标解读手册（推导 + 实测 + 冲突案例）](#16-指标解读手册推导--实测--冲突案例)
- [17. 质检四关卡](#17-质检四关卡机制代码变更后必过)
- [18. 故障排查速查表](#18-故障排查速查表)
- [19. 性能与耗时预算](#19-性能与耗时预算)
- [20. Python API 速查](#20-python-api-速查)
- [附录：最小示例数据集](#附录最小示例数据集)

---

## 0. 实验运行脚本（experiments/run/ —— 先读这一节）

> **新增**：`experiments/run/` 是"**只跑模拟 + 落 CSV**"的运行脚本集，
> 类比 CMAQ 的 `run.cctm.csh`（只跑输运写输出，绘图另用 NCL）。
> 它把 9 个实验封装成一条命令，输出统一落到 `output/run_*/`。

### 0.1 output/ 与 figures/ 的准确关系（常见误解）

| 目录 | 是什么 | 会不会被删 |
|---|---|---|
| **`output/`** | **实验产出的 CSV 数据**（每个实验一个 `run_X/` 子目录） | 只有 `make clean` 会删；`run_*.py` 只写自己的子目录 |
| **`figures/gen_all.py` / `gen_cartoon.py`** | **机制示意图**（教学图，纯手绘，**不读任何模拟结果**） | — |
| **`figures/fig_response_v3.py` / `fig_indicators_ridge.py`** | **工程图**（**读 `output/` 的 CSV** 才画） | — |

**所以你的理解是对的**：`figures/` 里的**工程图确实是通过输出文件可视化的**；
但 `figures/` 里还有一类**机制示意图**，它们与模拟结果无关。

### 0.2 输出会不会被自动删除？——不会（除 make clean）

```bash
make clean           # ⚠ 会 rm -rf output（唯一危险命令）
make clean-build     # ✓ 新增：只删二进制+生成码，保留所有实验结果
make clean-output    # ✓ 新增：只删 output/（需输入 y 确认）
```

做实验期间请用 `clean-build`。也可把结果放仓库外：
```bash
OZONE_OUT=$HOME/ozone_results python experiments/run/run_all.py
```

### 0.3 十个运行脚本（A–I + M，**不含伴随反演**）

| 脚本 | 实验 | 产物 | 耗时 |
|---|---|---|---|
| `run_A_scenarios.py` | 六情景政策沙盘 | `run_A/<SCEN>_{conc,irr,ipr}.csv` + `summary.csv` | ~3 s |
| `run_B_surface.py` | 响应面（4D） | `run_B/batch.csv` + **17 个 `indicator_matrix_*.csv`** + `slices.csv` | ~15 s |
| `run_C_ridge.py` | 屋脊跨路线对比 | `run_C/ridge_{batch,mcm,compare}.csv` | ~16 s |
| `run_D_process.py` | IPR 六账 + IRR 反应级 | `run_D/ipr_ledger.csv`（74×6）`closure.csv` `irr_top.csv` | ~1 s |
| `run_E_indicators.py` | 12 指标工具箱 | `run_E/<tag>_indicators.csv` `_timeseries.csv` `compare.csv` | <1 s |
| `run_F_osse.py` | OSSE 同化技能 | `run_F/obs_samples.csv` `skill_rmse.csv` `skill_delta.csv` | ~4 s |
| `run_G_hourlyPA.py` | 逐小时过程 + IPR↔IRR | `run_G/hourly_{ledger,share,rxn_top}.csv` `consistency.csv` | ~1 s |
| `run_H_sensitivity.py` | OAT 敏感性（±10%） | `run_H/jacobian.csv` `sensitivity.csv` `oat_runs.csv` | ~11 s |
| `run_I_scenario_matrix.py` | 情景 × 路线矩阵 | `run_I/matrix_{o3max,delta}.csv` `direction.csv` | ~4 s |
| **`run_M_met_series.py`** | **气象观测序列驱动 + 单变量归因** | `run_M/summary.csv` `attribution.csv` `o3_timeseries.csv` | ~6 s |

**总驱动**：`run_all.py`（顺序 A→B→D→G→E→F→H→M→I→C，全跑 **78 s**）

```bash
make all && make runtime_v3.dat                  # 前置
python experiments/run/run_all.py                # 全跑十个
python experiments/run/run_all.py A D            # 只跑指定
python experiments/run/run_all.py --list         # 列出清单
```

### 0.3b 气象序列输入（实验 M 的完整打法）

**能。** 6 个气象量都可以写 `file:` 走时间序列，机制与 BC 浓度约束完全同源
（`cfgbuild` 编译进 `runtime.dat` 的 `SERIES` 段）。

```bash
# 1) 已有 6 条 121 点逐时样例
ls input_v3/constraints/met_*_demo.txt
#    met_TEMP_demo.txt  met_PRES_demo.txt  met_RH_demo.txt
#    met_BLH_demo.txt   met_RAIN_demo.txt  met_JFAC_demo.txt

# 2) 切换驱动方式：改 input_v3/meteo.config（文件末尾已备好注释块，取消注释即可）
#    TEMP = file:constraints/met_TEMP_demo.txt
#    ...

# 3) 一键跑完整对照实验（常数基准 / 六变量全序列 / 六个单变量）
python experiments/run/run_M_met_series.py
```

**实验 M 的对照设计**（关键：把"平均强迫"与"时间结构"分离）：

| run | 气象驱动 | 作用 |
|---|---|---|
| `M0_const` | 六变量全取各自 120 h **均值** | 对照基准（同一平均强迫，无时间结构） |
| `M1_all` | 六变量**全部**逐时序列 | 真实驱动 |
| `M2..M7_*_only` | 每次只让**一个**变量走序列 | 归因 |

**实测结果**（`output/run_M/attribution.csv`，O3 峰值相对常数基准）：

| 变量 | ΔO3max (ppb) | 说明 |
|---|---|---|
| BLH | **−45.4** | 主导项：BLH 序列均值 918 m，但夜间塌到 154 m 使积累与稀释的时序严重不对称 |
| TEMP | **+16.8** | 次之：日均温相同时，日间高温时段恰好与光解峰值重合 → 放大 O3 |
| JFAC | +2.6 | 云盖日削弱光解，但晴天日的增强部分抵消 |
| RH | −1.2 | 改变 H2O → OH 初生，量级小 |
| RAIN | +0.2 | 5 天仅 6 h 降水，湿沉降窗口太窄 |
| PRES | 0.00 | 5 天 ≤1 hPa，经 M 的效应可忽略（实务中取常数即可） |

⚠ **单变量效应不可加和**（−45.4 + 16.8 + 2.6 − 1.2 + 0.2 = −27.0，而 `M1_all` 实测 −33.1）。
非线性来自 O3 对 BLH/TEMP 的耦合（温度影响速率常数又影响数密度，BLH 影响稀释又影响
沉降）。所以 `attribution.csv` 的末列"占总效应百分比"会出现 **−137% 和 +50.6%**——
这不是 bug，是"单个因子的贡献超过总量"的典型标志，写论文时应当写成
"效应不可加，需用响应面/全因子设计"而不是硬凑成 100%。

**三条硬规则**（`cfgbuild` 编译期强制）：

1. **气象序列里绝不能写 `# unit: ppb` 这类浓度单位指令行**。那条指令的作用是"乘以
   浓度换算因子"，写进气象文件会把 298.15 K 静默变成 7.33e12（×2.46e10）。
   现在 `cfgbuild` 遇到会直接报错（`read_series(allow_unit=False)`），
   并有专门的回归用例守卫。
2. 序列**不循环**，超出首末时刻钳制到端点值 → 想覆盖 5 天就写满 120 h。
3. 每序列 ≤128 点，且与 BC 序列**共享上限 8 的池子**。启用 6 个气象序列后，
   BC 侧最多再用 2 个 `file:` 约束。

### 0.4 输入假设（与你的要求一致）

> 所有 `run_*.py` **假设 `runtime_v3.dat` 已由 cfgbuild 编译完成**——
> 只做**前置检查**（存在性/大小），不自动编译。需要改数据时先手工：
> `make runtime_v3.dat`（或 `python src/cfgbuild.py input_v3 ...`）。

例外：`run_F_osse.py`（需为每个 τ 生成变体）、`run_H_sensitivity.py`（需为每个 ±扰动生成
变体）会在**自己输出目录内**建 `in_*/` 并调 cfgbuild —— 但它们不触碰仓库里的 `input_v3/`。

### 0.5 详细说明

完整用法、输出契约、mcm 前置构建、与旧脚本的对应关系：
**[experiments/run/README.md](experiments/run/README.md)**

---

## 1. 环境与首次构建

### 1.1 依赖清单

> 完整版（版本号、路径、依赖分层、从零复现、自检脚本）见 [ENVIRONMENT.md](ENVIRONMENT.md)；
> 一键体检：`sh check_env.sh`

| 依赖 | 位置 | 版本 | 用途 | 缺失后果 |
|---|---|---|---|---|
| gfortran / gcc | `~/gcc-13.2.0/bin` | 13.2.0 | 编译求解核 | `make` 失败 |
| Python | miniconda venv | 3.12 | cfgbuild/诊断/实验脚本 | 脚本无法运行 |
| numpy + matplotlib | 同 venv | — | 数值/绘图 | 诊断与图失败 |
| SUNDIALS CVODE | `~/soft/atmenv/cvode` | 2.7.0 | ODE 求解（fcmix 与 C API） | 链接失败 |
| KPP | `~/soft/src/KPP-3.2.0` | 3.2.0 | **仅历史**（KPP 路线已退役） | 无影响 |

`Makefile` 已固化上述路径（`PY`/`FC`/`CVODELIB` 变量），一般无需改。

### 1.2 首次构建

```bash
cd ozone_lab
make all     # 构建 fac_v3 + fac_batch_omp 两个工程可执行体
make test    # 全链自检（见 1.3）
```

构建链内部做了什么（理解它有助排错）：

```
mechanism/mechanism_3_full_chain.fac   (74sp/131rx 机制唯一事实源)
      │ mechc.py  ← 只在 .fac 变更时触发
      ▼
facsimile_impl/gen_v3/oz_mech.f90      (生成码: ns/nr/nq 常量 + update_rates + 计量表)
      + gen_v3/oz.reac  (每反应的反应物索引)
      + gen_v3/oz.prod  (每反应的产物索引+系数)
      │ gfortran
      ▼
fac_v3 / fac_batch_omp                 (两个二进制)
      + modules/runtime_dat_mod.f90    (runtime.dat 解析 + 气象/光解求值)
      + modules/chem_phys_mod.f90      (物理四账 + ATOL 调优, 双胞胎共享)
```

### 1.3 首次检查清单（4 勾选项，全过才继续）

| # | 检查 | 期望输出 |
|---|---|---|
| 1 | `make smoke` | `[PASS] FAC-V3 冒烟 (conc 26 行, O3末=…)` + `1/1 PASS` |
| 2 | `python tests/functional_test.py` | `functional_test: 27/27 PASS`（含 4 金值哨兵 + 7 条 R10-R13 防线 + sat2surf 5 项 + 2 条负向守卫） |
| 3 | `python audit/element_audit.py mechanism/mechanism_3_full_chain.fac` | `全部 131 反应元素守恒 PASS` |
| 4 | `ls facsimile_impl/fac_v3 facsimile_impl/fac_batch_omp` | 两个二进制存在 |

❌ 任何一项失败 → 先看 §18 故障排查，不要继续做实验（垃圾进垃圾出）。

### 1.4 目录地图（每个目录的职责）

| 目录/文件 | 职责 | 你会不会改它 |
|---|---|---|
| `mechanism/` | 机制唯一事实源（.fac）；`teaching/` 存 15sp 教学存档档 | 扩机制时改 |
| `input_v3/` | **12 类生效配置** + `constraints/` 序列 | **最常改这里** |
| `templates/` | 5 类数据落位模板（人读格式合同） | 只读（照它落位） |
| `runtime_v3.dat` | cfgbuild 产物（Fortran 唯一输入） | 不手改（自动生成） |
| `facsimile_impl/` | 工程主线两个可执行体 + 生成码 | 改物理逻辑时 |
| `mcm_city/` | MCM 对照子树的独立机制与二进制 | 机制扩链后跑一次 |
| `modules/` | 共享 Fortran 模块（runtime_dat_mod / chem_phys_mod） | 改物理/解析时 |
| `src/` | Python 编译前端（cfgbuild / mechc / mechc_c / fac2kpp） | 改输入格式时 |
| `diagnostics/` | 诊断库（diagnostics.py + indicators.py） | 改指标定义时 |
| `audit/` | 元素审计与 C 证人 | 加机制时跑 |
| `tests/` | 冒烟 + 功能测试 | 加哨兵时 |
| `experiments/run/` | **10 个实跑脚本**（A–I + M，只落 CSV 不绘图） | 加实验时 |
| `experiments/*.py` | 旧版实验脚本（含终端报告/绘图/反演） | 一般不用 |
| `output/` | 所有实验产物（已 gitignore） | 清理用 `make clean-build` |
| `audit_reports/` | 历次自审过程报告（已 gitignore，结论已进文档） | 只读 |
| `figures/` | 机制图集（SCI 版 + 漫画版）与工程图 | 改图时 |

### 1.5 一次最小完整工作流（跑通它你就入门了）

```bash
cd ozone_lab
make all && make test                       # ① 构建 + 自检
cd facsimile_impl
./fac_v3 SCEN_BAU ../output/my_first ../runtime_v3.dat   # ② 跑一个情景（0.2–0.9 s）
cd ..
python3 - <<'PY'                             # ③ 看指标
import sys; sys.path.insert(0, "diagnostics")
from indicators import route_indicators
d = route_indicators("output/my_first", "mechanism/mechanism_3_full_chain.fac")
for k in ("O3_ppb", "FNR", "L_NOx_Q", "OPE"): print(f"{k:10s} {d[k]:.3f}")
PY
# 期望: O3_ppb≈58.4  FNR≈3.43  L_NOx_Q≈2.01  OPE≈5.32
```

---

## 2. 化学基础（读实验前必读）

不懂机制就看不懂输出。本节用最短篇幅建立**看 CSV 与指标所需的全部化学知识**。

### 2.1 三大循环与两条 ★ 通道（一切的基础）

**第一循环·NOx 光稳态（Leith，时间尺度 ~100 s）**

```
NO₂ + hν (J₁) → NO + O(³P)        光解，O₃ 生成入口
O(³P) + O₂ + M → O₃               三体结合，τ≈μs
O₃ + NO → NO₂ + O₂                滴定
```
净效果为零——只锁定比值 **[O₃] = J₁[NO₂]/(k[NO])**。所以"只有 NOx 和阳光"永远造不出净 O₃。

**第二循环·HOx（O₃ 净生成通道 1）**

```
O₃ + hν (J₂) → O(¹D) + O₂          ~5–10% 走此通道
O(¹D) + H₂O → 2 OH                 OH 唯一强源（另 ~50% 被 N₂/O₂ 猝灭）
CO + OH (+O₂) → CO₂ + HO₂          载体
HO₂ + NO → NO₂ + OH  ★             把 NO 翻成 NO₂ 但不消耗 O₃
```

**第三循环·ROx（O₃ 净生成通道 2）**

```
VOC + OH → RO₂                     燃料进入
RO₂ + NO → RO + NO₂  ★             第二净生成通道
RO + O₂ → HCHO + HO₂               碎渣变自由基
HCHO + hν (J₄) → 2 HO₂ + CO        一张 HCHO 光解 = 两个 HO₂
```

**为什么 ★ 能净产生 O₃**：`HO₂+NO` 把 NO 氧化为 NO₂ 且**没有动用 O₃**；随后 NO₂ 光解
产生 O → O+O₂ → O₃，于是净多了一个 O₃。走滴定边（O₃+NO）再走光解则是净零。
**P(O₃) = k[HO₂][NO] + k[RO₂][NO]** —— 这就是 IRR 里的净化学生成。

### 2.2 mechanism_3 V2.3 十段链条（对应图 H 的十个进料管）

| 链 | 名称 | 关键内容 | 涉及物种（举例） |
|---|---|---|---|
| ① | 无机核心 + HONO 晨爆 | J₁–J₇ 光解 + 18 热反应；Leith 光稳态；OH/HO₂ 循环；HONO 夜间积累→日出爆发 | NO NO₂ O₃ O¹D O OH HO₂ H₂O₂ CO HNO₃ HONO |
| ② | NO₃/N₂O₅ 夜化学 | O₃+NO₂→NO₃；NO₃+NO₂⇌N₂O₅（MCM KMT03/04 **Troe 逐字拷贝**）；水解→HNO₃；**白天 J₅/J₆ 光解闭环** | NO₃ N₂O₅ |
| ③ | CH₄ + PAN 储库 | CH₄→CH₃O₂→CH₃O→HCHO；CH₃CHO→CH₃CO₃；CH₃CO₃+NO₂⇌PAN（Ea≈112 kJ/mol） | CH₄ CH₃O₂ CH₃OOH CH₃O CH₃CHO CH₃CO₃ PAN |
| ④ | ISOP 二代链 | ISOP→ISOPBO₂→三分支（NO★/HO₂/NO₃）→MVK/MACR→MPAN/MGLY→CH₃CO₃ 回灌 | ISOP ISOPBO₂ ISOPNO₃ ISOPNIT ISOPDNT ISOPOOH MVK MACR MPAN MGLY |
| ⑤ | 硫链 | DMS+OH(日)/NO₃(夜)→SO₂→SO₃→H₂SO₄；+MSA 夜通道 | SO₂ SO₃ H₂SO₄ DMS MSA |
| ⑥ | 卤素闭环 | X+O₃→XO；XO+NO→X+NO₂；XO+NO₂→XONO₂ / XO+HO₂→HOX → hν 再生 | Cl ClO ClONO₂ Br BrO BrONO₂ HCl HBr HOCl HOBr |
| ⑦ | V2.2 人为集总 | ALK/OLE/ARO1/ARO2/KET 五类；**OLE+O₃ 臭氧化**；芳烃→GLY+MGLY | ALK OLE ARO1 ARO2 KET GLY |
| ⑧ | V2.2 支撑化学 | **PNA 温度计**（HO₂+NO₂⇌PNA）；N₂O₅ 白天光解；OH+HNO₃ 再循环 | PNA |
| ⑨ | V2.3 烷基硝酸盐 | RONO₂/OLENIT/ISOPNIT 的 OH 氧化通道（NOx 再循环） | RONO₂ OLENIT |
| ⑩ | V2.3 单萜迷你链 | TERP + OH/O₃/NO₃ 三入口 → TERPRO₂/TERPNIT/TERPOOH → SOA 前体 | TERP TERPRO₂ TERPNIT TERPOOH TERPCHO |

### 2.3 一天中的化学时钟（对照 CSV 的 h 列读）

| 时段 | 主导过程 | CSV 里你会看到 |
|---|---|---|
| 00–06 | 夜间：NO 滴定 O₃、NO₃/N₂O₅ 储库、HONO 地表积累 | O₃ 下降、NO 高、N₂O₅ 升、HONO 升 |
| 06–09 | 日出：HONO 光解晨爆 → OH 骤升；BLH 抬升稀释 | OH 尖峰、HONO 骤降、O₃ 被稀释掉一点 |
| 09–15 | 白天主战场：自由基循环全速，O₃ 攀升 | O₃ 单调升、HO₂/RO₂ 高、P(O₃) 峰 |
| 15–18 | 峰值后：PAN 累积、H₂O₂/HNO₃ 累积 | O₃ 见顶、PAN 峰 |
| 18–24 | 傍晚：交通 NO 滴定，光解停止 | O₃ 快速降、NO 升、NO₃ 起 |

### 2.4 物种分类（决定你该关注哪一列）

| 类别 | 成员 | 特征 |
|---|---|---|
| **自由基**（稳态，寿命 μs–s） | OH HO₂ CH₃O₂ RO₂ 家族 O O¹D | 浓度极小（1e4–1e8），**看比值不看绝对值** |
| **无机核心** | NO NO₂ O₃ HNO₃ HONO NO₃ N₂O₅ PNA | 主指标来源 |
| **有机燃料与产物** | CH₄ CO HCHO CH₃CHO ISOP ALK OLE ARO1 ARO2 KET TERP MVK MACR | 排放 + 降解 |
| **有机储库** | PAN MPAN RONO₂ OLENIT ISOPNIT ISOPDNT TERPNIT | 长寿命、可远距离输送 |
| **氢过氧化物** | H₂O₂ CH₃OOH ISOPOOH MVKOOH MACROOH ALKOOH OLEOOH TERPOOH | 自由基汇（干沉降 0.3–0.5 cm/s） |
| **硫/卤素** | SO₂ SO₃ H₂SO₄ DMS MSA Cl ClO ClONO₂ Br BrO BrONO₂ HCl HBr HOCl HOBr | 海洋场景核心 |

---

## 3. 数据如何进程序（12 类配置详解）

### 3.1 五段管线（程序只认 runtime_v3.dat）

```
templates/*.template          ① 格式合同（人读；程序零直读, grep 可证）
     │ 人工落位（唯一的"手工"环节）
     ▼
input_v3/constraints/*.txt    ② 两列 "t[s] 值" + "# unit: cm3|ppb|ppm"
     │ 在 config 中登记（bc.config / meteo.config 的 file: 行）
     ▼
input_v3/*.config (12 类)     ③ 生效配置（活数据）
     │ src/cfgbuild.py 编译：物种存在性/单位/单调/点数≤128/序列≤8/槽位完备
     ▼
runtime_v3.dat                ④ Fortran 唯一输入（零文本解析）
     │ modules/runtime_dat_mod.f90 → chem_phys_mod 物理四账
     ▼
fac_v3 / fac_batch_omp        ⑤ 输出 CSV → diagnostics/indicators
```

#### 源码级验证：谁在真正读文件（逐 open 点，带行号）

**`templates/` 零代码引用**——可自行验证：

```bash
grep -rn "templates" --include="*.py" --include="*.f90" --include="*.sh" \
     --include="Makefile" src/ modules/ facsimile_impl/ mcm_city/ tests/ audit/ experiments/
# 输出为空 → 模板从未被任何程序打开
grep -rn "\.template" src/ modules/ facsimile_impl/ tests/ experiments/
# 输出为空 → 也没有 *.template 通配读取
```

**`src/cfgbuild.py` 的全部文件读取点**（11 个显式 + 2 个 glob = 12 类）：

| 行号 | 读取内容 | 用途 |
|---|---|---|
| :109 | `glob("mechanism_*.fac")` | 机制（物种表/速率/反应） |
| :121 | `glob("photolysis-rates_*")` | 光解 l/m/n 参数表 |
| :127 | `read_kv("photolysis.config")` | MODE/LAT/DOY/JMAP |
| :158 | `read_kv("meteo.config")` | 气象七量（可含 `file:`） |
| :178 | `read_series(...)` | **meteo 的 `file:` 序列本体** |
| :188 | `read_kv("bc.config")` | TAU_MIX + 常数/`file:` 登记 |
| :197 | `read_series(...)` | **BC 的 `file:` 序列本体** |
| :210, :221 | `emissions.config` | 排放（两遍：先扫 `unit:` 指令，再读数） |
| :244 | `drydep.config` | 干沉降 vd |
| :255 | `wetdep.config` | 湿沉降 a/b |
| :265 | `read_kv("dilution.config")` | DIL |
| :270 | `scenarios.config` | 情景因子表达式 |
| :290 | `model.parameters` | 天数/步长/rtol/atol |
| :379 | `write(runtime.dat)` | **唯一输出** |

**Fortran 端只有一个输入文件**：

| 位置 | 语句 | 说明 |
|---|---|---|
| `modules/runtime_dat_mod.f90:71` | `open(newunit=u, file=path, status="old", action="read")` | 只打开 `runtime_v3.dat` |
| `facsimile_impl/oz_fac_v3.f90:107` | `call read_runtime(trim(rtfile))` | rtfile 默认 `runtime_v3.dat`（CLI 第 3 参可覆盖） |

**结论（三点，均可 grep 复现）**：
1. **模板是格式合同**：`templates/*.template` 不被任何程序读取（零引用、零通配）
2. **真正的数据入口只有两条**：`input_v3/*.config`（12 类，cfgbuild 读）与
   `input_v3/constraints/*.txt`（被 config 的 `file:` 指向后由 cfgbuild 读）
3. **运行时零文本解析**：Fortran 只读 cfgbuild 产出的 `runtime_v3.dat`

### 3.2 十二类配置逐个详解（含当前实测值）

#### ① `mechanism_*.fac` —— 机制（唯一事实源）

```
VARIABLE  <物种名列表，空格分隔，分号结尾>
<速率名> = <表达式>             例: KHO2NO = 8.6D-12*EXP(270/TEMP)
RO2 = <过氧自由基列表>          例: RO2 = CH3O2 CH3CO3 ISOPBO2
% <速率名|J<n>|数值表达式> : <反应物> = <产物> ;
FORMULA  <物种 元素计数>         例: HCHO C1 H2 O1   (审计与 NOy 判定的事实源)
```
- `J<n>` = 光解槽位（1–41），槽位必须有非零 l 参数（cfgbuild 校验）
- 元素审计读 FORMULA 段；indicators 的 NOy 判定也从这里自动提取含 N 物种

#### ② `initialConcentrations.config` —— 初始浓度

```
# 初值: 物种 值 单位(cm3|ppb|ppm)
NO      5.0e8    cm3          # = 20 ppb
NO2     2.0e8    cm3
O3      3.0e11   cm3
CO      2.0e12   cm3
CH4     4.5e13   cm3
HCHO    2.0e9    cm3
PAN     5.0e9    cm3
ISOP    2.5e9    cm3
```
- 未列出的物种 = 0；单位可写 `ppb`（自动 ×2.46e10）
- IC 只影响起始，5 天模拟下前 1–2 天是"建立期"，末 24 h 才代表准稳态

#### ③ `meteo.config` —— 气象七量（常数或序列）

```
TEMP   = 298.15                        # K
PRES   = 1013.25                       # hPa
RH     = 60                            # %
RADIAT = 0                             # W/m²（保留，当前不参与 J）
BLH    = file:constraints/met_BLH.txt  # m（序列, 必须 >0）
RAIN   = file:constraints/met_RAIN.txt # mm/h（序列）
JFAC   = 1.0                           # 光解缩放（云订正 0–1）
```
- 7 个量都支持 `file:` 序列
- BLH 用于：干沉降 `vd/BLH`、BLH 膨胀稀释 `−y·dB/BLH`（仅增长段）
- **BLH 序列最小值必须 > 0**（cfgbuild 强制）
- **六变量全序列的现成样例**：`constraints/met_{TEMP,PRES,RH,BLH,RAIN,JFAC}_demo.txt`
  （各 121 点 / 0–120 h 逐时）。切换方式与配套实验 M 见 **§0.3b**。

#### ④ `constraints/*.txt` —— 序列数据本体

```
# t[s]  BLH[m]        ← 注释行（第一行可有列名）
      0    300        ← 自模拟起点秒；严格递增
  21600    300
  28800   1200
  43200   1500
```
- 两列空白分隔；行内 `#` 后为注释
- **单位指令行 `# unit: ppb` 只对浓度序列（BC 约束）有意义**。
  气象序列**不允许**出现该行——它的作用是"乘以浓度换算因子"，
  误写会把 298.15 K 静默变成 7.33e12。`cfgbuild` 现已编译期拒绝
  （`read_series(allow_unit=False)`），并有回归用例守卫。
- 上限：**每序列 ≤128 点、全池 ≤8 条**（气象 + BC 共享）
- 外延：超出时间轴取端点值（**不是外推**）；缺测：直接省行，自动线性插值

**现有序列文件清单**（`input_v3/constraints/`）：

| 文件 | 点数 | 用途 |
|---|---|---|
| `met_BLH.txt` / `met_RAIN.txt` | 7 | 默认启用的轻量序列（5 天快照） |
| `met_{TEMP,PRES,RH,BLH,RAIN,JFAC}_demo.txt` | 各 121 | 站点逐时样例（实验 M 用） |
| `bc_O3_demo.txt` | 24 | 地面 O₃ 站点逐时 → BC 张弛约束 |
| `bc_HCHO_demo.txt` | 24 | HCHO 站点/卫星 → BC 张弛约束 |
| `bc_NO2_sat_demo.txt` | 24 | S5P TROPOMI NO₂ 柱 → BC 张弛约束 |


#### ⑤ `bc.config` —— 边界/观测约束（张弛同化）

```
TAU_MIX = 86400                              # 张弛时间常数 [s]
O3  = 6.2e11                                 # 常数直接写
CO  = 1.5e12
CH4 = 4.5e13
# O3 = file:constraints/o3_station_MYSITE.txt   ← 序列约束（推荐）
```
物理含义：`dy/dt += (观测插值 − y)/TAU_MIX`——箱体被"拉向"观测。
- TAU_MIX 越小 = 约束越强（86400 s≈1 天=弱；7200 s=2 h=强，见实验 F）
- 序列优先、常数兜底（`chem_phys_mod` 物理账③）；每物种独立选择

#### ⑥ `photolysis.config` + `photolysis-rates_MCM` —— 光解

```
MODE = MCM            # MCM 参数化（SZA 驱动）或 DIRECT（查表）
LAT  = 40.0           # 纬度（决定 SZA 日变化）
DOY  = 172.0          # 年积日（172≈夏至）
JMAP = 13=D4.0e-5,14=D1.0e-5,...   # 槽位→l 参数（速率表未定义时的直给）
```
- J 计算式：`J = JFAC·l·cosᶻᵐ·exp(−n/cosχ)`；`cosχ≤1e-3 → J=0`（夜间）
- 槽位 1–41，机制用到 23 个；MAPPING 缺省恒等

#### ⑦ `emissions.config` —— 排放（四列格式）

```
# 物种   基线速率        日型  振幅
NO       2.157e+07     BOX   0.9
NO2      1.078e+06     BOX   0.9
CO       1.126e+07     BOX   0.9
HCHO     5.815e+04     BOX   0.9
CH4      2000          FLAT  0.0
ISOP     2.037e+06     BOX   0.9
SO2      4e+05         BOX   0.9
DMS      1e+04         FLAT  0.0
ALK      2.622e+07     BOX   0.9
OLE      2.134e+07     BOX   0.9
ARO1     1.199e+07     BOX   0.9
ARO2     7.992e+06     BOX   0.9
KET      1.5e+06       BOX   0.9
TERP     1.0e5         BOX   0.75
```
- 单位：**molec·cm⁻³·s⁻¹**（箱体单位体积排放率）
- `BOX`：06–18 时 `E0×(0.1+0.9·amp·sin²(π(h−6)/12))`，**夜间固定 0.1×E0**
- `FLAT`：全天 `E0×1.0`；振幅 ∈ [0,1]（强制）；日型大小写敏感
- 未列出 = 零排放

#### ⑧ `dilution.config` —— 稀释
```
DIL = 0.0     # [1/s]，整体一阶稀释
```

#### ⑨ `drydep.config` —— 干沉降
```
O3      0.8    on      # vd [cm/s] on/off；k = vd/BLH
HNO3    2.0    on
HCHO    0.3    on
H2SO4   2.0    on
OLEOOH  0.5    on      # 氢过氧化物同族 0.3–0.5
```
- **新物种若无化学汇，必须在此配沉降**（OLEOOH 曾因此无界累积）

#### ⑩ `wetdep.config` —— 湿沉降
```
HNO3    5.0e-5  0.8    # Λ = a·RAIN^b
```

#### ⑪ `scenarios.config` —— 政策情景
```
SCEN_BAU   = E_NOx*1 E_VOC*1 E_CO*1 E_CH4*1 JFAC*1
SCEN_LNNOX = E_NOx*0.3 E_VOC*1 E_CO*1 E_CH4*1 JFAC*1
SCEN_HVOC  = E_NOx*1 E_VOC*3 E_CO*1 E_CH4*1 JFAC*1
SCEN_HOT   = E_NOx*1 E_VOC*1 E_CO*1 E_CH4*1 JFAC*1 DT+5
SCEN_CLOUD = E_NOx*1 E_VOC*1 E_CO*1 E_CH4*1 JFAC*0.6
SCEN_ISOPX = E_NOx*0.5 E_VOC*2 E_CO*1 E_CH4*1 JFAC*1
```
- 只支持乘法因子与 `DT+n`（加法，温度偏移）
- 组映射（`runtime_dat_mod.build_groups`）：`E_NOx`→NO/NO₂；
  `E_VOC`→HCHO/CH₃CHO/ISOP/DMS/MVK/MACR/MSA/SO₂/ALK/OLE/ARO1/ARO2/KET；
  `E_CO`→CO；`E_CH4`→CH₄

#### ⑫ `model.parameters` —— 模拟控制
```
5.0        # 模拟天数
3600.0     # 输出步长 [s]（=1 h）
1.0e-6     # RTOL
1.0e3      # ATOL 基线（主量按 chem_phys_mod 调优表覆盖）
```

### 3.3 排放单位换算完整推导（三种来源）

单位统一为 **molec·cm⁻³·s⁻¹**。

**(a) 从 ppb/h 换算**（站点/箱模型常用）

1 ppb 在 298 K、1 atm 下 = 2.46×10¹⁰ molec/cm³。"每小时涨 1 ppb"即每秒涨 1/3600 ppb：

```
E [molec/cm³/s] = C_ppb_per_hour × 2.46e10 / 3600 = C × 6.833e6
```
**例**：NO 20 ppb/h（早高峰城市典型）→ `20 × 6.833e6 = 1.37e8`。

**(b) 从质量清单换算**（t/yr 落到箱体）

```
E = Q[t/yr] × 1e6 / M[g/mol] × 6.022e23 / (3.154e7 [s/yr] × A[cm²] × H[cm])
```
**例**：NO₂ 1000 t/yr，A=100 km²=1e12 cm²，H=1 km=1e5 cm：
```
E = 1000×1e6/46 × 6.022e23 / (3.154e7 × 1e12 × 1e5) = 4.15e9
```

**(c) 从面通量换算**（mol/m²/s）

```
E = F[mol/m²/s] × 6.022e23 / (1e4 [cm²/m²] × H[cm])
```
**例**：1e-7 mol/m²/s、H=1 km → `1e-7×6.022e23/(1e4×1e5) = 6.02e7`。

**(d) 合理性自检（稳态法）**

稳态浓度 ≈ E × τ。NO 的 τ 由 O₃ 滴定决定：`τ = 1/(k[O₃]) ≈ 1/(1.8e-14×7e11) ≈ 79 s`。
E=2.157e7 → SS ≈ 1.7e9 ≈ **70 ppb**，与城市量级吻合 ✓。
若 SS 偏离 1–3 个量级，回查单位（最常见：忘记 /3600 或漏 H）。

### 3.4 五类模板落位（templates/ 是合同，不是输入）

| 模板 | 你提供什么 | 落位到 | 登记到 |
|---|---|---|---|
| `surface_o3_station` | 站点 O₃ 序列（ppb） | `constraints/o3_station_<站名>.txt` | `bc.config` 一行 |
| `surface_aqi_station` | AQI 宽表（μg/m³） | 先离线拆分换算成各物种序列 | 同上（PM、NO 不接） |
| `surface_voc_station` | VOC 物种序列 | `constraints/<voc>_station_*.txt` | `bc.config` |
| `emissions_inventory` | 排放清单（换算成 molec/cm³/s） | **直接编辑** `emissions.config` | — |
| `satellite_column` | 卫星柱浓度（特殊语义） | **需离线转表面等效 ppb** | `bc.config` |

> ⚠ 卫星模板是唯一"引擎不直读"的例外：0-D 箱无气柱概念，必须先离线反演成
> 表面等效 ppb（OSSE 实验 F 演示了这一决策）。

### 3.5 数据硬约束（违反即 cfgbuild 拒绝，错误信息会指路）

| 约束 | 值 | 违反表现 |
|---|---|---|
| 每序列点数 | ≤ 128 | `序列 N 点 > Fortran 端上限 128` |
| 全池序列数 | ≤ 8 | `序列总数 9 > Fortran 端上限 8` |
| 时间轴 | 严格递增、正值 | `时间非单调 @ t=…` |
| 单位指令 | `# unit: cm3\|ppb\|ppm` | `未知单位 'x'` |
| BLH 序列 | 全 > 0 | `BLH 序列含非正值` |
| 排放基线 | ≥ 0 | `基线速率 < 0` |
| 排放日型 | 严格 `BOX`/`FLAT` | `日型 'box' 非法` |
| BOX 振幅 | [0, 1] | `振幅 ∉ [0,1]` |
| 物种名 | ∈ 机制 VARIABLE | `排放未知物种 / IC 未知物种` |

### 3.6 编译边界（改什么才需要编译）

| 你改了什么 | 命令 | 实测耗时 |
|---|---|---|
| 12 类 config / constraints / 排放 / 情景 / IC | `make runtime_v3.dat` | **0.37 s** |
| 只换情景名跑 | 直接 `./fac_v3 SCEN_X …` | **0.19–0.87 s** |
| `.fac`（加删物种/反应/改速率） | `make v3`（mechc 重生成 + gfortran） | **~8 s** |
| Fortran 逻辑（物理/组映射/诊断） | `make v3` | ~8 s |

**数据 ≠ 代码**：边界画在 `.fac` 上。make 按时间戳判断（配置改了代码没变 →
`Nothing to be done`）。

---

## 4. 路线选择与编译边界

### 4.1 三条现役路线的定位

| | ① fac_v3 | ② fac_batch_omp | ③ mcm_run |
|---|---|---|---|
| **一句话** | 显微镜：单点问"为什么" | 望远镜：曲面问"在哪转折" | 体检器：无机块有没有被带偏 |
| 机制 | V2.3 74sp/131rx（gen_v3） | 同左（共享生成码） | MCM 29sp/71rx（独立子树） |
| 内核 | fcmix CVODE-BDF+DENSE | **C-ABI 直连 CVODE + OpenMP** | fcmix CVODE-BDF+DENSE |
| 输入 | runtime.dat 全通路 + 观测同化 | 同左 + 网格参数 | CLI 四参（days/pref/fnox/fvoc） |
| 一次运行 | 1 情景 | 81+ 格点 | 1 案例 |
| 输出 | 逐时 74 物种 + 131 IRR + 六账 | 24 列指示因子（末 24h 均值） | 29 物种 + 71 IRR + 六账 |
| 生产角色 | ✓ 主力 | ✓ 批量 | ✗ 对照 |

### 4.2 决策树

```
问"为什么/过程/政策响应" → ① fac_v3
问"整个曲面的形状/转折点" → ② fac_batch_omp
问"机制改完无机块还正常吗" → ③ mcm_run（扩链后跑一次）
问"∂O₃/∂kᵢ 逐反应灵敏度" → ① + experiments/rate_sensitivity.py
```

### 4.3 显微镜 vs 望远镜（同一个内核，两个窗口）

- **内核四样全同**：机制（gen_v3）· 物理（chem_phys_mod 四账）·
  数值（CVODE-BDF+DENSE）· 输入基线（runtime_v3.dat）
- **差别只在接口层**：情景空间 / 并行方式 / 输出深度
- **互证容差 ±10~15%**：同格点不是逐位相同（batch 有 jfac/ΔT 覆写路径与
  IPR 均值窗差异）；"逐位一致"仅指 batch 自身并行/串行对照
- **标准接力**：广度找现象（②）→ 深度查原因（①）→ 广度验稳健（②）

---

## 5. 实验 A：单箱基线与政策情景

### 5.1 科学问题

一团代表性城市空气，在给定排放/气象下：**O₃ 涨到多少？政策减谁有效？**

### 5.2 前置条件

- `make test` 全绿
- `input_v3/` 的 12 类配置已就绪
- 想改政策：先看 `scenarios.config` 现有 6 个情景

### 5.3 逐步命令

```bash
# 步骤 1: 检查排放与情景配置
grep -v "^#" input_v3/emissions.config      # 15 个物种的排放
grep "^SCEN" input_v3/scenarios.config      # 6 个情景定义

# 步骤 2: 编译数据（改过任何 config 后必须；没改也幂等）
make runtime_v3.dat
# 期望: 情景: ['SCEN_BAU', 'SCEN_LNNOX', ...] （打印情景列表）

# 步骤 3: 跑六个情景（在 facsimile_impl/ 下）
cd facsimile_impl
for S in BAU LNNOX HVOC HOT CLOUD ISOPX; do
  echo "=== SCEN_$S ==="
  ./fac_v3 SCEN_$S ../output/expA_$S ../runtime_v3.dat
done
```

### 5.4 预期输出（实测样例）

```
FAC-V3 情景: SCEN_BAU (fNOx= 1.00, fVOC= 1.00)
FAC-V3 ok (  0.70 s, nst=5329, O3max= 101.61 ppb)
```
- 第一行：匹配到的情景名与两个因子
- 第二行：CPU 时间、总积分步数（nst）、**全程 O₃ 峰值**（ppb）
- `nst≈5300` 是健康值（若 >2 万说明数值困难，查 §18）

### 5.5 产物逐列解读

**(a) `output/expA_BAU_conc.csv`** —— 逐时浓度

```
h,NO,NO2,O3,O1D,O,OH,HO2,H2O2,CO,HCHO,HNO3,...
 96,  5.83904E+07,  4.92624E+10,  9.75556E+11, ...
 97,  6.89242E+07,  5.08632E+10,  8.48995E+11, ...
```
- 第 1 列 `h` = 小时序号（**从 96 起 = 第 4 天 0 时**，只写末 24 h）
- 单位 molec/cm³；换 ppb 除 2.46e10
- **可能看到极小负值**（如 `O = -3.1E-43`）——BDF 的浮点噪声，10⁻⁴³ 量级等于零，正常
- 读法：`O3` 列除以 2.46e10 得 ppb；`OH` 列看日变化形状（不看绝对值）

**(b) `output/expA_BAU_irr.csv`** —— 逐反应速率
```
h,r1,r2,...,r131
```
`r<k>` 对应 `.fac` 的第 k 个反应行（见 `gen_v3/oz.reac`/`oz.prod` 索引）。

**(c) `output/expA_BAU_ipr.csv`** —— 六本账
```
h,NO_P,NO_L,NO_E,NO_D,NO_B,NO_X,NO2_P,...
```
每物种六账：`P`机制产 / `L`机制耗 / `E`排放 / `D`干湿沉 / `B`BC+稀释 / `X`BLH膨胀。
**这是 fac_v3 独有的深度诊断**——回答"这个物种的变化由哪个过程主导"。

### 5.6 V2.3 六情景实测与判读

| 情景 | 因子 | O₃max (ppb) | 末24h均值 (ppb) | 判读 |
|---|---|---|---|---|
| BAU | 全 1 | 101.61 | 58.43 | 基线 |
| LNNOX | NOx×0.3 | 53.9 | 32.8 | **减 NOx 大幅降 O₃ → 处于屋脊右侧（NOx 饱和）** |
| HVOC | VOC×3 | 101.4 | 58.2 | 加 VOC 几乎无变化 → VOC 不是瓶颈 |
| HOT | T+5K | 112.6 | 63.5 | +5 K 涨 11 ppb：PAN 分解释放 NOx + 速率温度依赖 |
| CLOUD | JFAC×0.6 | 84.3 | 48.2 | 光解减 40% → O₃ 降 17% |
| ISOPX | NOx×0.5+VOC×2 | 71.3 | 42.2 | 组合，LNNOX 主导 |

**三个判读结论**：
1. **本配置处 NOx 饱和区**（LNNOX 降 47%、HVOC 无感）→ 政策应减 NOx；
2. **温度敏感**（+5 K → +11 ppb）→ 气候变暖会加剧本区域 O₃；
3. **光解是第一驱动力**（JFAC×0.6 → −17%）。

**反例检查**：若你的数据跑出 LNNOX 反升 → 你在屋脊左侧（NOx-limited），
体制已变，回到实验 B 重新定位。

### 5.7 变异玩法

| 想研究 | 怎么做 |
|---|---|
| 排放的绝对值影响 | 改 `emissions.config` 基线（注意 §3.3 换算） |
| 新政策组合 | `scenarios.config` 加 `SCEN_MINE = E_NOx*0.7 E_CO*0.5 JFAC*1` → 重编译 → 跑 |
| 初始条件敏感性 | 改 `initialConcentrations.config`（前 2 天建立期会吸收大部分差异） |
| 想看全程而非末 24h | 把 `oz_fac_v3.f90` 的 `step >= nsteps-24` 改 `step >= 0`，`make v3` |
| 站点观测同化 | 见 §14 配方 1 |
| 延长模拟 | `model.parameters` 第 1 行改天数 |

### 5.8 常见坑

| 坑 | 症状 | 解 |
|---|---|---|
| 忘记重编译数据 | 结果与上次完全一致 | `make runtime_v3.dat` |
| 在错误目录跑 | `Cannot open file 'gen_v3/oz.reac'` | `cd facsimile_impl` |
| 期望看到第 1 天 | CSV 只有 24 行 | 设计如此（只写末 24h） |
| 把负的浮点尾数当错误 | `O=-3E-43` | 10⁻⁴³ = 零，正常 |
| 用 CSV 首末 24h 比"累积" | 两者是同一窗口 | 用 IPR 分账或改输出条件 |

---

## 6. 实验 B：响应面与指示因子屋脊

### 6.1 科学问题

把 "NOx 排放 × VOC 排放" 铺成曲面：**O₃ 屋脊（转折点）在哪？
哪条指示因子能提前指示体制？**

### 6.2 前置条件

- 实验 A 已跑通（同一批配置）
- `make v3`（需 fac_batch_omp 二进制）

### 6.3 逐步命令

```bash
# 步骤 1: 默认 9×9 网格（几何网格 0.3×1.5^(i-1)，0.30–7.69）
make batch
# 内部: v3 → runtime_v3.dat → fac_batch_omp 9 ../output/batch_v3.csv ../runtime_v3.dat

# 步骤 2: 加密网格（如 13×13）
cd facsimile_impl
./fac_batch_omp 13 ../output/batch_13.csv ../runtime_v3.dat

# 步骤 3: 出图
cd .. && make figures
# → figures/response_surface_v3.png + figures/indicators_vs_ridge.png
```

### 6.4 预期输出

运行 81 案例，4 线程 ~47 s，零失败。产物 `output/batch_v3.csv`（曲面原始值）
与 `output/batch_v3_indicators.csv`（24 列）。

### 6.5 产物逐列解读（24 列）

| 列 | 含义 | 用途 |
|---|---|---|
| 1–4 | `fnox, fvoc, jfac, dt` | 网格坐标 |
| 5–8 | `O3max, O3mean, HCHOmax, NO2mean` | 地面量 |
| 9 | `FNR_mean` | 体制判据主指标 |
| 10–13 | `H2O2HNO3, O3NOy, HCHONOy, O3NOz` | 浓度比族 |
| 14 | `NOy_ppb` | NOx 存量（17 物种加权） |
| 15–17 | `LnOxQ, PH2O2PHNO3, OPE` | 过程量族 |
| 18–23 | `*_pm` 系列（12–18h 正午窗变体） | Sillman 原始口径 |
| 24 | `nst` | 步数（数值健康诊断） |

### 6.6 V2.3 实测与判读

- **屋脊 fnox = 2.28**（O₃max 125.5 ppb），右侧回撤 90% → 强烈 NOx 饱和
- **FNR 过渡带（1–2）与屋脊位置一致**（一致性 0.67）→ 指示因子方向性验证通过
- 画图读法：取 `fvoc≈1` 的行，画 O₃-fnox 切面 → 峰即屋脊；叠 FNR 曲线看穿越点

**判读要点**：屋脊右侧（fnox>2.28）= 减 NOx 有效；左侧 = 减 NOx 可能反升；
FNR 与屋脊重合 = 指示因子在本机制上可信；显著偏离 → 需按机制重标定（见实验 C）。

### 6.7 变异玩法

| 想研究 | 怎么做 |
|---|---|
| 温度维 | 扩展 fac_batch_omp 的 dt 列表（见 .f90 头注） |
| 光解维 | 扩展 jfac 列表 |
| 更细网格 | `./fac_batch_omp 21 …`（上限 32×32） |
| 只看中间区 | 改网格生成式 `fs(i)=0.3*1.5**(i-1)` 的起点/比值 |

### 6.8 常见坑

| 坑 | 解 |
|---|---|
| 用等距轴画图导致形状失真 | 网格是**几何**的，用 log 轴 |
| 以为 batch 也扫了 fco/fch4/jfac/dt | 默认固定 1/1/1/0，只扫 fnox×fvoc |
| 哨兵超差却继续用 | 先查机制/数值是否被动过（±10–15% 金值） |
| 网格 >32 | 数组上限，会越界 |

---

### 6.9 响应面能力全景：能跑多少种"类 EKMA"曲面？

> **⚠ 先澄清术语**：这里的"四维"**不是空间维**。ozone_lab 与 AtChem2 一样是
> **0-D 箱体模型**（状态向量只有 `ns` 个物种，**无 x/y/z、无平流扩散**）。
> 所谓"四维"指**参数扫描空间**：fnox × fvoc × jfac × ΔT 四个**无量纲缩放因子**
> 的组合。每组组合跑**一个独立的 0-D 箱**——即"**262,144 个互不相干的箱子**"，
> 不是"一个 4D 场"。
>
> **结论先行**：`fac_batch_omp` 是**四维参数扫描驱动器**（fnox × fvoc × jfac × ΔT），
> 容量 **262,144 个独立 0-D 案例**；配 **19 个可作图响应量**，理论上可产出
> **4,980,736 个独立曲面**。日常使用 9×9 = 81 案例 → 1,539 个曲面。

#### 6.9.1 四维参数空间（源码 `fac_batch_omp.f90:415-480`）

> **每个案例都是完整的 0-D 箱**：循环体内 `c_cvcreate` 新建一个 CVODE 实例、
> 状态向量 `N_VNew_Serial(ns)` 只有 74 个物种分量、初值重置为 `y0`
> （见 `:481-492`）。**案例之间零耦合**——这既是"独立箱子"的证明，也是它天然
> 可 OpenMP 并行的原因（每线程一套 `mem/N_Vector/user_data`）。

| 参数轴 | 含义 | 取值范围 | 容量 |
|---|---|---|---|
| **fnox** | NOx 排放缩放因子（× `emissions.config` 基线） | 几何序列 0.3×1.5^(i−1)；`ng` 可调 | **32** |
| **fvoc** | VOC 排放缩放因子（乘 `E_VOC` 组全部物种） | 同上 | **32** |
| **jfac** | 光解缩放（云/气溶胶订正） | CLI 逗号列表，如 `0.6,0.8,1.0,1.2` | **16** |
| **dt** | 温度偏移 [K] | CLI 逗号列表，如 `0,5,10` | **16** |

**注意**：四个轴全是**无量纲缩放因子/偏移量**，不是空间坐标。`fnox=2.0` 意味着
"把 NO 排放设为基线的 2 倍"，不是"向东推进 2 个格"。

```bash
# 二维（默认）：只扫 NOx × VOC
fac_batch_omp 9 out.csv rt.dat

# 三维：加光解维（4 档）
fac_batch_omp 9 out.csv rt.dat "0.6,0.8,1.0,1.2"

# 四维：再加温度维（3 档）→ 9×9×4×3 = 972 案例
fac_batch_omp 9 out.csv rt.dat "0.6,0.8,1.0,1.2" "0,5,10"
```

**已实测验证四维**（2×2×2×2 = 16 案例，12.18 s，零失败）：
```
fnox,fvoc,jfac,dt
0.300,0.450,0.800,0.000
0.450,0.300,0.800,0.000
...
```

#### 6.9.2 19 个可作图响应量（24 列中的非坐标列）

| 类别 | 数量 | 变量 |
|---|---|---|
| **浓度量** | 4 | `O3max` · `O3mean` · `HCHOmax` · `NO2mean` |
| **指示因子（全窗）** | 9 | `FNR_mean` · `H2O2HNO3` · `O3NOy` · `HCHONOy` · `O3NOz` · `NOy_ppb` · `LnOxQ` · `PH2O2PHNO3` · `OPE` |
| **指示因子（正午窗 12–18h）** | 6 | `FNR_pm` · `H2O2HNO3_pm` · `O3NOy_pm` · `HCHONOy_pm` · `O3NOz_pm` · `NOy_pm` |
| 诊断（非曲面） | 1 | `nst`（步数，数值健康检查） |

> **每个响应量 × 每个网格切面 = 一张等值线图**。例如：
> - `O3max` vs (fnox, fvoc) → **经典 O₃ 响应面（屋脊图）**
> - `FNR_mean` vs (fnox, fvoc) → **体制判别图**
> - `OPE` vs (fnox, fvoc) → **NOx 饱和度图**
> - 加 jfac 维后：`O3max` vs (fvoc, jfac) @ fnox=1.5 → **云量敏感性面**
> - 加 ΔT 维后：`O3max` vs (fnox, ΔT) @ fvoc=1 → **温度敏感性面**

#### 6.9.3 组合数统计

| 配置 | 案例数 | × 19 响应量 = 曲面数 |
|---|---|---|
| **9×9（默认）** | 81 | **1,539** |
| 13×13 | 169 | 3,211 |
| 9×9×4 jfac | 324 | 6,156 |
| 9×9×3 ΔT | 243 | 4,617 |
| 9×9×4×3（四维） | 972 | 18,468 |
| 16×16×4×4 | 4,096 | 77,824 |
| **32×32×16×16（满容量）** | **262,144** | **4,980,736** |

#### 6.9.4 其它实验载体（非 batch 网格，但也是"响应面"性质）

| 实验 | 扫描维度 | 案例数 | 产物 |
|---|---|---|---|
| **实验 C `ridge5.py`** | fnox（15 点 logspace）× **多路线**（batch/mcm） | ~30 | `figures/ridge5_routes.png` |
| **实验 F `obs_osse.py`** | **τ_mix 4 档**（2/6/24/72 h）× 站点物种 × 卫星反演 | ~20 | RMSE 技能表 + 过冲指标 |
| **实验 H `sens_adj.py`** | OAT 扰动 ±10% × 各排放组 + 截断 SVD | ~10 | 灵敏度斜率 + 反演对比 |
| **实验 I `scen5.py`** | 6 情景 × 3 路线 | 18 | 方向一致性矩阵 |
| `scan_suite.py` | 9×9 × 6 变量（旧 fac_run 载体，现走 batch） | 81 | `ozone_lab_response_suite.png` |

#### 6.9.5 与"经典 EKMA"的对应关系

| | 经典 EKMA | 本实现（排放空间） |
|---|---|---|
| 模型类型 | **0-D 箱**（Dodge 1977 也是箱模型） | **0-D 箱**（同） |
| 双轴 | [NOx]₀ × [VOC]₀（**初始浓度**） | fnox × fvoc（**排放因子**，乘基线速率） |
| 等值线量 | O₃ 峰值 | O₃max + **另外 18 个响应量** |
| 维度 | 2D | **4D**（可加光解、温度） |
| 物理过程 | 通常无汇 | 含沉降/BC/BLH 膨胀（屋脊存在的前提） |

**所以"类 EKMA 曲面"的数量口径**：

| 口径 | 数量 |
|---|---|
| **双轴 O₃ 等值线（最严格 EKMA 类比）** | **1 种**（fnox×fvoc 的 O3max） |
| **双轴 × 任选响应量** | **19 种**（4 浓度 + 15 指标） |
| **三轴/四轴切片曲面** | 4D 网格的任意 2 维切片 → **6 种轴组合** × 19 响应量 = **114 种曲面族** |
| **理论上限（所有格点组合）** | 262,144 案例（每案例 19 个数） |

#### 6.9.6 经典响应面清单：能实现 4 款，缺 3 款

> 按"是否属于**经典文献响应面类型**"逐条核对（源码实读 + 实测）。

##### ✅ 已实现的 4 类（其中 FNR 单列为 2a，因其地位特殊）

| # | 经典响应面 | 轴 | 本实现方式 | 文献出处 |
|---|---|---|---|---|
| **1** | **排放空间 O₃–NOx–VOC 等值线图**（"EKMA-like"，主用） | fnox × fvoc | `make batch` → `response_surface_v3.png` | Sillman 1995；Dodge 1977 的排放空间变体 |
| **2a** | **★ FNR 响应面**（HCHO/NO₂ 体制判别面） | fnox × fvoc | `make batch && make figures` → `response_surface_v3.png` 的 (b) 面板 + `indicators_vs_ridge.png` 的 9 面板之 (b) | Sillman 1995；Tonnesen-Dennis 2000 |
| **2b** | 其它指示因子等值线（H₂O₂/HNO₃、O₃/NOy、OPE、L_NOx/Q…） | 同上 | 同一 CSV 的其余 14 列 | Kleinman 2005 |
| **3** | **光解（云量）敏感性面** | fnox × jfac | `fac_batch_omp 9 out rt "0.5,0.7,1.0,1.3"` | 常规敏感性实验 |
| **4** | **温度敏感性面** | fnox × ΔT | `fac_batch_omp 9 out rt "1.0" "0,3,6,9"` | 常规敏感性实验（PAN 分解主导） |

**4 轴可组成的二维切面**：C(4,2) = **6 种**（fnox×fvoc / fnox×jfac / fnox×ΔT / fvoc×jfac / fvoc×ΔT / jfac×ΔT），
× 19 响应量 = **114 个曲面**。

##### ❌ 缺失的 3 款（有明确原因）

| # | 经典响应面 | 为什么没有 | 变通办法 |
|---|---|---|---|
| **5** | **真 EKMA（Dodge 1977 原义：初始浓度空间）** | `y0` 直接取自 `runtime.dat` 的 IC 段，**无任何标度乘子**（`:485 yp = y0`）；原载体 `kpp_run` 已随路线退役删除 | 手工改 `initialConcentrations.config` + 外部 shell 循环（无内置网格）；或恢复 IC 标度轴（约 +15 行代码） |
| **6** | **气象多变量面**（BLH / RH / RAIN / 风速） | 只有 TEMP 可经 `dts` 偏移（`eval_met_dt`）；RH/BLH/RAIN 是 config 文件里的常数/序列，无乘子轴 | 改 `meteo.config` + 外部循环；或扩展 `runp_t` 加 `fblh/frh` 字段 |
| **7** | **CO / CH₄ 独立轴** | `scen_fco`/`scen_fch4` 在情景系统里存在，但 **batch 网格中硬编码为 1.0**（`:479 rp%fco = 1.0d0; rp%fch4 = 1.0d0`） | `fac_v3` 走 `SCEN_*` 情景表可设；批量需改这一行 + 加轴（约 +20 行） |

##### 一句话回答

> **经典响应面里，ozone_lab 能实现 4 款**：排放空间 O₃ 等值线（主用）、
> 指示因子等值线、云量（光解）敏感性面、温度敏感性面。
> **缺 3 款**：真 EKMA（IC 空间）、气象多变量面、CO/CH₄ 独立轴——
> 原因分别是"无 IC 标度轴""气象无乘子轴""batch 硬编码固定"，
> 都不是物理限制，而是驱动层未暴露参数（各需 15–20 行代码扩展）。

##### 验证：这 4 款都实测过吗？

| 款 | 实测状态 |
|---|---|
| 1 排放空间 O₃ 等值线 | ✅ 81 案例实测，屋脊 fnox=2.28 / O₃max 125.5 ppb |
| 2 指示因子等值线 | ✅ 同一 CSV 的 19 列，`indicators_vs_ridge.png` 已出图 |
| 3 光解敏感性 | ✅ 四维实测（2×2×2×2=16 案例含 jfac 轴） |
| 4 温度敏感性 | ✅ 同上（含 dt 轴）；SCEN_HOT(+5K) 也有实测 112.6 ppb |
| 5 真 EKMA | ❌ 不可跑（无轴）；`ridge5.py` 里 kpp 的 IC 空间代码已退役 |
| 6 气象多变量 | ❌ 不可跑（需改 config） |
| 7 CO/CH₄ 轴 | ❌ batch 固定 1.0；`fac_v3` + 自定义 `SCEN_*` 可单点试 |

#### 6.9.7 FNR 响应面专项（最核心的一张图）

> **FNR = [HCHO]/[NO₂]** 是 ozone_lab 里唯一**同时具备"物理意义 + 卫星可观测 + 阈值判据"**
> 的响应量，所以它的响应面值得单列。

##### 实测 9×9 响应面（`output/batch_v3_indicators.csv` 的 `FNR_mean` 列）

行 = fnox↑（NOx 排放因子），列 = fvoc→（VOC 排放因子）：

| fnox \ fvoc | 0.30 | 0.45 | 0.68 | 1.01 | 1.52 | 2.28 | 3.42 | 5.13 | 7.69 |
|---|---|---|---|---|---|---|---|---|---|
| **0.300** | 3.06 | 4.40 | 6.35 | **9.03** | 12.72 | 18.39 | 27.59 | 41.15 | 60.26 |
| **0.450** | 2.14 | 3.11 | 4.61 | 6.83 | 9.95 | 14.53 | 21.99 | 34.12 | 50.22 |
| **0.675** | 1.42 | 2.14 | 3.23 | 4.94 | 7.51 | 11.26 | 17.09 | 26.96 | 42.44 |
| **1.012** | 0.68 | 1.41 | 2.21 | **3.45** | 5.42 | 8.46 | 13.09 | 20.58 | 33.36 |
| **1.519** | 0.25 | 0.64 | 1.43 | 2.33 | 3.77 | 6.10 | 9.77 | 15.55 | 25.05 |
| **2.278** | 0.11 | 0.22 | 0.63 | **1.49** | 2.53 | 4.21 | 7.01 | 11.50 | 18.66 |
| **3.417** | 0.07 | 0.10 | 0.18 | **0.64** | 1.60 | 2.81 | 4.83 | 8.21 | 13.68 |
| **5.126** | 0.05 | 0.07 | 0.10 | 0.16 | 0.67 | 1.76 | 3.20 | 5.63 | 9.72 |
| **7.689** | 0.03 | 0.06 | 0.08 | **0.10** | 0.15 | 0.72 | 1.98 | 3.71 | 6.64 |

##### 三条判读规律（都能从上表直接读出）

**① FNR 随 NOx 严格单调下降**（固定任一 fvoc 列，从上往下数值递减）
```
fvoc=1.01 切面: 9.03 → 6.83 → 4.94 → 3.45 → 2.33 → 1.49 → 0.64 → 0.16 → 0.10
                └────────── 单调下降 90 倍 ──────────┘
```
→ **FNR 是可靠的 NOx 暴露度量**（这正是它作为体制判据的基础）。

**② FNR 随 VOC 严格单调上升**（固定任一 fnox 行，从左往右递增）
```
fnox=1.012 行: 0.68 → 1.41 → 2.21 → 3.45 → 5.42 → 8.46 → 13.09 → 20.58 → 33.36
```
→ FNR 也是可靠的 VOC 暴露度量。

**③ 阈值带与屋脊重合**（关键验证）
```
fvoc≈1 切面 FNR 穿越:
  fnox=1.519  FNR=2.33   （仍在 NOx-limited）
  fnox=2.278  FNR=1.49   ← 落入 1–2 过渡带
  fnox=3.417  FNR=0.64   ← 已穿入 <1（VOC-limited）
  而 O₃max 屋脊恰在 fnox=2.28（125.5 ppb）→ 阈值带正中屋脊 ✓
```
→ 这**验证了 Sillman 判据在本机制上的成立**：FNR 穿 1–2 带的位置 = O₃ 由升转跌的位置。

##### 两个图里的 FNR 位置

| 图 | 面板 | 呈现方式 |
|---|---|---|
| `figures/response_surface_v3.png` | **(b) 面板** | FNR 响应面，**log10 色标** + O₃ 屋脊等值线叠加（看两者是否重合） |
| `figures/indicators_vs_ridge.png` | **9 面板之 (b)** | FNR 响应面，等值线取 **[1.0, 2.0, 4.0]**（阈值带显式画出） |

重生成命令：
```bash
make batch && make figures
# → figures/response_surface_v3.png  (3 面板: O₃max / FNR / 屋脊切片)
# → figures/indicators_vs_ridge.png  (3×3 九面板, FNR 在 b 位)
```

##### 与卫星 FNR 的衔接

| 地面 FNR（本响应面） | 卫星 FNR（柱比） |
|---|---|
| HCHO/NO₂ **浓度比**（末 24h 均值） | HCHO柱/NO₂柱 |
| 来自 `FNR_mean`（全窗）/ `FNR_pm`（正午 12–18h） | 需先经 `sat2surf.py` 转表面等效 |
| 阈值 1–2（本机制标定） | 同阈值，但受 BLH 假设影响 |

**用途**：把卫星观测的 FNR 落到这张响应面上 → **反推当前气团落在哪个排放体制** →
据此判断"减 NOx 还是减 VOC"（政策决策链）。

---

#### 6.9.8 实用建议

> 所有下面这些命令跑的都是**一批独立的 0-D 箱**，产出的是**参数-响应表**（CSV），
> 再由 Python 画成等值线图。**不是空间分布图**（那是 CMAQ/WRF-Chem 的事）。

```bash
# ① 要经典 EKMA 式屋脊图（最常用）
make batch && make figures          # → response_surface_v3.png

# ② 要云量敏感性面（第二类常见需求）
cd facsimile_impl
./fac_batch_omp 9 ../output/batch_jf.csv ../runtime_v3.dat "0.5,0.7,1.0,1.3"

# ③ 要温度敏感性面
./fac_batch_omp 9 ../output/batch_dt.csv ../runtime_v3.dat "1.0" "0,3,6,9"

# ④ 要完整四维（注意耗时：972 案例 ≈ 12 分钟）
./fac_batch_omp 9 ../output/batch_4d.csv ../runtime_v3.dat "0.6,0.8,1.0,1.2" "0,5,10"
```

> ⚠ **耗时线性于案例数**：81 案例 ≈ 47–70 s（4 线程），故 972 案例 ≈ 10–14 min，
> 262,144 案例 ≈ 47–70 小时。满容量跑之前先想清楚是否需要。

---

### 6.10 四大分析能力状态核实（观测约束 / 敏感性 / 伴随 / 响应面）

> 逐项对照源码与实跑，含**能力边界**与**不是真伴随**这类必须说清的点。

#### 6.10.-1 六项能力速查（可复制命令）

> 每项能力的最短可执行路径。详细八段式讲解见对应章节。

```bash
# ═══ 0. 前置（一次性）═══
cd ozone_lab
make all && make test                    # 构建 + 自检（24→25 断言）

# ═══ 1. 背景边界（A 用法：定环境）═══
grep -v '^#' input_v3/bc.config          # 看当前背景值(O3 25ppb/CO 61ppb/CH4 1829ppb)
vim input_v3/bc.config                   # 改成你要的背景气候值（常数形式）
make runtime_v3.dat && cd facsimile_impl
./fac_v3 SCEN_BAU ../output/bg_demo ../runtime_v3.dat

# ═══ 2. 观测约束（B 用法：定真值）═══
# ① 准备约束文件（两列 + 单位头）
printf '# unit: ppb\n      0 30.0\n   3600 35.2\n   7200 42.8\n' \
    > ../input_v3/constraints/o3_demo.txt
# ② 登记（把常数行改成 file: 行）
sed -i '' 's|^O3 .*|O3 = file:constraints/o3_demo.txt|' ../input_v3/bc.config
# ③ 重编译 + 运行
cd .. && make runtime_v3.dat && cd facsimile_impl
./fac_v3 SCEN_BAU ../output/assim_demo ../runtime_v3.dat
# ④ 验证同化生效（IPR 的 B 账应显著非零）
cd .. && python3 -c "
import sys; sys.path.insert(0,'diagnostics'); import diagnostics as dg
sp,ipr=dg.ipr_from_fac('output/assim_demo'); print('O3_B =', f\"{ipr['O3']['B']:.3e}\")"

# ═══ 3. 响应面（4D 参数扫描）═══
make batch                               # 9×9 = 81 案例（含 FNR 等 19 响应量）
make figures                             # → response_surface_v3.png + indicators_vs_ridge.png
# 四维：加光解 + 温度
cd facsimile_impl && ./fac_batch_omp 9 ../output/b4d.csv ../runtime_v3.dat \
    "0.6,0.8,1.0,1.2" "0,5,10"

# ═══ 4. FNR 响应面（最核心）═══
make batch && make figures
# → response_surface_v3.png 的 (b) 面板：FNR + O₃ 屋脊叠加
# → indicators_vs_ridge.png 的 (b) 面板：FNR 等值线 [1,2,4]
python3 -c "
import csv; rows=list(csv.DictReader(open('output/batch_v3_indicators.csv')))
for r in rows:
    if abs(float(r['fvoc'])-1.012)<0.01: print(r['fnox'], r['FNR_mean'], r['O3max'])"

# ═══ 5. 情景分析 ═══
grep '^SCEN' input_v3/scenarios.config   # 看 6 个预置情景
cd facsimile_impl
for S in BAU LNNOX HVOC HOT CLOUD ISOPX; do
    ./fac_v3 SCEN_$S ../output/expA_$S ../runtime_v3.dat
done
cd .. && python3 experiments/scen5.py    # → 情景×路线方向一致性矩阵

# ═══ 6. 过程分析（IPR 六账 + IRR 131 反应）═══
cd facsimile_impl && ./fac_v3 SCEN_BAU ../output/pa_demo ../runtime_v3.dat && cd ..
python3 experiments/hourly_pa.py output/pa_demo    # → 逐小时六账 + 反应级 top4 + 一致性
python3 experiments/ipr_closure_check.py          # → 六账闭合校验（1.004–1.019）
python3 -c "
import sys; sys.path.insert(0,'diagnostics'); import diagnostics as dg
sp,ipr=dg.ipr_from_fac('output/pa_demo')
print('O3 六账:', {k: f'{v:+.2e}' for k,v in ipr['O3'].items()})"

# ═══ 7. 敏感性 OAT（局部梯度）═══
python3 experiments/sens_adj.py          # → jacobian.csv (24×9) + SVD 反演

# ═══ 8. 观测系统模拟实验 OSSE（评价同化）═══
python3 experiments/obs_osse.py          # → τ 扫描 + RMSE 技能 + 过冲指标

# ═══ 9. 卫星柱转换（配合观测约束）═══
python3 src/sat2surf.py --selftest       # 先自检
python3 src/sat2surf.py 柱数据.raw ../input_v3/constraints/sat_NO2_surf.txt \
    --blh-file ../input_v3/constraints/met_BLH.txt
```

> **每项对应的详细章节**：① 背景边界 → §3.2⑤ / ② 观测约束 → §14.1 + §10
> ③ 响应面 → §6 / ④ FNR → §6.9.7 / ⑤ 情景 → §5 + §13 / ⑥ 过程 → §8 + §11
> ⑦ OAT → §12 / ⑧ OSSE → §10 / ⑨ 卫星 → §14.4

#### 6.10.0 ⚠ 先分清：敏感性 OAT ≠ 响应面（这是两类不同的分析）

这两个常被混为一谈，但在数学与工程上完全不同：

| | **OAT 敏感性**（`sens_adj.py`） | **响应面**（`fac_batch_omp`） |
|---|---|---|
| **自变量** | **1 个**（逐个参数单独扰动） | **2–4 个**（同时扫描） |
| **扫描方式** | 中心差分 ±10%（**局部线性化**） | 全网格枚举（**全局形状**） |
| **输出形态** | **雅可比矩阵** J（24 观测维 × 9 参数） | **2D 等值线图**（9×9 → 19 个响应量） |
| **回答的问题** | "参数动 1%，结果动多少%"（**梯度**） | "整个空间什么形状？转折点在哪？"（**拓扑**） |
| **成本** | 18 次运行（9 参数 × 2 方向） | 81 次运行（9×9） |
| **结果用途** | 反演、不确定性传播、参数可辨识性 | 体制判别、政策扫描、屋脊定位 |
| **典型产物** | `jacobian.csv` | `response_surface_v3.png` |
| **能否互相替代** | ❌ 不能 | ❌ 不能 |

**为什么不能互替**：
- OAT 只在**标称点附近**有效（±10%），无法告诉你"屋脊在哪"——需要全局扫描
- 响应面给的是**离散格点上的值**，无法直接给梯度——需要差分/拟合
- 实测印证：OAT 斜率 **42.6** ppb/unit vs 网格局部斜率 **52.3**（19% 偏差），
  因为二者在不同尺度上测量（小扰动 vs 大步长跨非线性区）

**正确的关系**：响应面是**全局地图**，OAT 是**局部坡度仪**；地图告诉你山在哪，
坡度仪告诉你脚下多陡。实践中先用响应面定位关键区，再在该点做 OAT 求精确梯度。

#### 6.10.1 观测约束（nudging 同化）—— ✅ 有，且是生产级

##### ⚠ 术语澄清：三个易混概念

| 概念 | 英文 / 代码 | 是什么 | 不是什么 |
|---|---|---|---|
| **① 观测约束** | Observation Constraint<br>`chem_phys_mod.f90` 物理账③ | **一个能力/机制**：把观测序列按 `1/τ` 速率拉进模型（**Newtonian Relaxation / Nudging，牛顿松弛**），属**资料同化（Data Assimilation）的简化形式** | 不是"仿真分析"——它是**动作**，不是**分析** |
| **② OSSE** | Observing System Simulation Experiment<br>`experiments/obs_osse.py` | **一个评价实验**：双胞胎法（真值→采样加噪→同化→RMSE 技能），检验"同化方法有没有用" | 不是同化方法本身——它是**检验方法的实验** |
| **③ 基于观测的模拟仿真分析** | — | **泛称**，可指 OSSE、可指观测约束运行、也可指反演 | **不是严格术语**，不宜当"观测约束"的全称 |

**三者关系**：
```
观测约束 (机制)  ← 被 OSSE (实验) 检验 → 得出技能指标 (RMSE/τ 曲线/B 账)
     ↑
  资料同化 (学科归类: nudging 是其中最简形式)
```

**所以准确说法**：
- 说"**观测约束**" = 指 **nudging 同化机制**（`bc.config` + `TAU_MIX`）
- 说"**OSSE**" = 指 **观测系统模拟实验**（`obs_osse.py`，用于评价）
- 说"**基于观测的模拟仿真分析**" = 口语泛称，**建议避免**在技术文档中使用

> 另外：AtChem2 里对应概念叫 **constrained species（约束物种）**，
> 与本项目的"观测约束"是同一机制的不同叫法。

##### ⚠ 关键：`bc.config` 一个文件承担**两种完全不同的物理含义**

这是最容易混淆的地方——**同一段代码（`chem_phys_mod` 物理账③）、同一个文件（`bc.config`）、
同一个参数（`TAU_MIX`），却在做两件性质不同的事**：

| 用法 | 填什么 | 典型数值 | 物理含义 | 学科归属 |
|---|---|---|---|---|
| **A. 背景边界（默认）** | **常数** | `O3 = 6.2e11`（25 ppb）<br>`CO = 1.5e12`（61 ppb）<br>`CH4 = 4.5e13`（1829 ppb） | **"盒子外面的大气长什么样"**——0-D 箱体与外界交换的空气成分。没有它，箱内 CO/CH₄ 会被耗尽、O₃ 会无限累积 | **边界条件（Boundary Condition）** |
| **B. 观测约束（同化）** | **`file:` 时间序列** | `O3 = file:constraints/o3_station_BJ.txt`<br>（30 → 35 → 43 ppb 逐时） | **"把真实观测塞进模型"**——强制模型跟随观测 | **资料同化（Data Assimilation / Nudging）** |

> **判断方法**（一眼分辨）：看 `bc.config` 里写的是**常数还是 `file:`**。
> 常数 → 背景边界；`file:` → 观测约束。

##### 为什么 0-D 箱**两者都需要**？

0-D 箱是一个"理想化的封闭气块"，但要让它模拟真实大气，必须解决两个问题：

```
问题 1: "盒子外的大气是什么？"     → 需要 背景边界 (A)   否则 CO/CH4 会耗尽, O3 会跑飞
问题 2: "盒子里的结果对不对？"     → 需要 观测约束 (B)   否则只能看模型自己的猜测
```

**举例（当前默认配置的 A 用法，含实测反证）**：

| 物种 | 当前背景值 | 作用 | **若删掉会怎样（实测）** |
|---|---|---|---|
| `CH4` | 4.5e13 = **1829 ppb** | 大气背景甲烷 | **崩到 4.1 ppb（降 446 倍）**——甲烷链整条断掉 |
| `CO` | 1.5e12 = 61 ppb | OH→CO→HO₂ 载体循环 | 靠排放勉强维持 88 ppb（形状失真） |
| `O3` | 6.2e11 = **25 ppb** | 自由对流层背景 | 峰值虚高 110 ppb（无背景约束，O₃ 自涨） |

**实测验证**（删空 bc.config 的背景物种后跑 SCEN_BAU）：
```
有背景边界: CH4 维持 ~1829 ppb (不随时间衰减)
无背景边界: CH4 6.6 → 4.1 ppb  ← 被 OH 消耗殆尽, 甲烷链断裂
             O3max = 109.95 ppb  ← 无背景约束, 数值虚高
```
**结论**：背景边界（A 用法）不是可选项——它是 0-D 箱能代表"真实大气中的一团空气"的**前提**。

**举例（B 用法，观测约束）**：
- 站点逐时 O₃ 观测 → 模型被"拉向"观测 → 用于检验机制或做 OSSE

##### 三者的角色分配（完整图）

| 通道 | 作用 | 类比 |
|---|---|---|
| **emissions.config** | 往箱子里**注入**污染物 | 进水管 |
| **bc.config 常数**（A） | 设定**盒子外**的背景大气 | 房间外的空气 |
| **bc.config `file:`**（B） | 用观测**校准**箱子里的值 | 恒温器（拉向目标） |
| **drydep/wetdep** | 从箱子里**移除** | 排水管 |
| **dilution** | 整体一阶稀释 | 换气扇 |
| **BLH 膨胀** | 混合层抬升导致浓度稀释 | 房间变大 |

##### 所以"观测约束"和"bc 序列张弛"是什么关系？

**是同一件事的两种叫法**（都指 B 用法）：

| 说法 | 指向 |
|---|---|
| 观测约束 | 强调**目的**（让模型跟随观测） |
| bc 序列张弛 | 强调**手段**（用 bc.config 的 `file:` 序列 + TAU_MIX 张弛） |

但**"bc" 本身不等于"观测约束"**——`bc.config` 还可用于**纯背景边界（A 用法，默认就是这种）**，
那时它跟同化无关，只是给 0-D 箱一个合理的环境。

> **一句话**：`bc.config` 是"盒子边界"的通用接口；填常数 = 背景气候值（物理边界），
> 填观测序列 = 数据同化（观测约束）。两者都服务于 0-D 模拟，但**一个定环境、一个定真值**。

| 项 | 实况 |
|---|---|
| 实现位置 | `modules/chem_phys_mod.f90:42-52`（物理账③） |
| 数学形式 | `dy/dt += (观测插值 − y) / TAU_MIX` —— 一阶松弛（Newtonian nudging） |
| 数据通路 | `bc.config` 的 `物种 = file:constraints/xxx.txt` → cfgbuild SERIES → 运行时每内步 `ser_interp(t)` 插值 |
| 强度控制 | `TAU_MIX`（86400=弱 / 21600=中 / 7200=强） |
| 分账追踪 | IPR 第三账 `B`（BC 张弛 + 稀释），可定量看同化贡献 |
| 验证实验 | `experiments/obs_osse.py`（OSSE：真值→采样→同化→RMSE 技能） |
| 实测结论 | 地面 HCHO τ=2h → RMSE **+27%**；卫星反演 E_NOx×1.34 → NO₂ RMSE +74%（观测算子误差被放大） |
| 可约束物种 | **32 个推荐**（详见 `templates/SPECIES_INPUT_MATRIX.md` §4.2）；自由基类**不建议**（寿命≪τ） |
| 关键限制 | 强 nudging 会**破坏光稳态**（NO₂ 尤其敏感，建议弱 τ）；1/τ 线性松弛不是最优估计理论 |

#### 6.10.2 敏感性扰动分析（OAT）—— ✅ 有

| 项 | 实况 |
|---|---|
| 实现位置 | `experiments/sens_adj.py`（主入口 `main()`） |
| 方法 | **OAT 中心差分**：9 个排放物种各 ±10%，重跑 `fac_v3`，构造雅可比 |
| 输出 | `jacobian.csv`（24 观测维 × 9 排放列）+ 各列范数 |
| 实测 | fac_v3 OAT **42.6 ppb/unit fnox** vs batch 网格局部斜率 52.3 → **19% 偏差**（小扰动 vs 大步长的非线性差异，非错误） |
| 成本 | 每个物种 2 次全跑 → 9 物种 = **18 次 fac_v3**（约 2 min） |
| 边界 | 局部线性；扰动 >50% 会跨越体制边界；只覆盖 `EMIT` 列表里的 9 个物种 |

#### 6.10.3 伴随分析 —— ⚠ **不是真伴随，是"有限差分雅可比 + 截断 SVD 反演"**

这是必须纠正的常见误解：

| | 真伴随（adjoint, 如 CMAQ-DDM / GEOS-Chem adjoint） | ozone_lab 实际做的 |
|---|---|---|
| 实现方式 | 反向求解伴随方程（需编译器 AD 或手写伴随代码） | **前向有限差分**：逐参数 ±10% 重跑，拼雅可比 |
| 成本 | **与参数个数无关**（1 次前向 + 1 次反向） | **线性于参数个数**（9 参数 = 18 次前向） |
| 精度 | 解析梯度（机器精度） | 数值梯度（受步长/非线性影响，实测 19% 偏差） |
| 代码状态 | — | **无伴随代码**（`grep adjoint` 零命中） |

**实际提供的是**（`sens_adj.py` ②段）：
```
J = ∂y/∂(E/E0)   (中心差分, 24×9)
SVD 分解 → 截断小奇异值 (S > S0×0.1) → 伪逆 → 反演 ΔE/E
```

**⚠ 实测能力很有限（2026-10-26 实跑）**：
```
奇异值条件数 = 3,153,320        ← 极高病态
保留 1/9 个模态                  ← 9 个排放参数里只有 1 个方向可辨识
```
原因：非线性光化学使 9 个排放物种对 24 维观测的响应**近乎共线**（尤其是
NO/NO₂/CO/ALK/OLE 同属"促进 O₃"方向，ISOP 也高度相关）。所以：

| 能做 | 不能做 |
|---|---|
| 算**单参数**灵敏度（OAT 斜率，实测 42.6 ppb/unit fnox） | 同时反演 9 个排放的**唯一解** |
| 给出"最可辨识的 1 个组合方向" | 逐物种精确反演（条件数 3e6，噪声会放大百万倍） |
| 演示 SVD 正则化方法本身 | 作为实用源解析工具 |

**改进方向**（若要做实用反演）：减少参数（先按排放组聚合到 3–4 个）、
增加观测维度（多时刻 + 多物种）、引入先验协方差（贝叶斯反演）。

> **准确的表述**：ozone_lab 提供的是**线性化反演（有限差分雅可比 + 截断 SVD 正则化）**，
> 不是伴随模式。若需真伴随，需引入 AD 工具（Tapenade/ADIFOR）或手写伴随方程——
> 本项目没有。

#### 6.10.4 响应面分析 —— ✅ 有（4 维参数扫描）

详见 §6.9。一句话：`fac_batch_omp` 扫 fnox×fvoc×jfac×ΔT（最多 262,144 个独立 0-D 案例），
输出 19 个可作图响应量 → 114 个曲面族（6 种二维切面 × 19 响应量）。
**已实现 4 款经典响应面，缺 3 款**（真 EKMA/IC 空间、气象多变量面、CO/CH₄ 轴）。

#### 6.10.5 四能力总表

| 能力 | 状态 | 核心实现 | 关键限制 |
|---|---|---|---|
| **观测约束** | ✅ 生产级 | `chem_phys_mod` 物理账③ + F3 序列 BC | 线性松弛非最优估计；强 τ 破坏光稳态 |
| **敏感性扰动** | ✅ 有 | `sens_adj.py` OAT 中心差分 | 局部线性；成本线性于参数数 |
| **伴随分析** | ⚠ **非真伴随** | 有限差分雅可比 + 截断 SVD 反演 | 数值梯度（19% 偏差）；**实测条件数 3.2e6、仅 1/9 模态可辨识** |
| **响应面** | ✅ 有（4D 参数扫描） | `fac_batch_omp` | 非空间网格；缺 3 款经典面 |

### 6.11 情景分析与过程分析（两项都是核心能力）

#### 6.11.1 情景分析（Scenario Analysis）—— ✅ 有，且是政策沙盘主力

**两层实现**：

| 层 | 机制 | 位置 |
|---|---|---|
| **情景表驱动**（主用） | `scenarios.config` 定义命名情景 → `fac_v3 SCEN_NAME` 直接跑 | `input_v3/scenarios.config` + `runtime_dat_mod.scen_build` |
| **情景矩阵对比** | 6+ 情景 × 多路线 → 方向一致性检验 | `experiments/scen5.py` |

**支持的 6 个情景因子**（源码 `scen_build`）：

| 因子 | 作用对象 | 语义 |
|---|---|---|
| `E_NOx*n` | NO、NO₂ | 排放组乘子 |
| `E_VOC*n` | HCHO/CH₃CHO/ISOP/DMS/MVK/MACR/MSA/SO₂/ALK/OLE/ARO1/ARO2/KET | 排放组乘子 |
| `E_CO*n` | CO | 排放组乘子 |
| `E_CH4*n` | CH₄ | 排放组乘子 |
| `JFAC*n` | 全部光解通道 | 光解缩放 |
| `DT±n` | TEMP | 温度偏移 [K] |

**当前 6 个预置情景**：BAU / LNNOX(NOx×0.3) / HVOC(VOC×3) / HOT(+5K) / CLOUD(JFAC×0.6) / ISOPX(NOx×0.5+VOC×2)

**实测结果**（`scen5.py` 实跑，O₃max [ppb]）：

| 路线 | BAU | LNNOX | HVOC | ISOPX | HOT | CLOUD | IC×2 | BG+50% |
|---|---|---|---|---|---|---|---|---|
| **fac_v3 74sp** | 101.61 | **53.94** | 111.26 | 73.13 | 112.56 | 84.30 | 101.61 | 103.61 |
| mcm_run 29sp | 16.22 | 14.66 | 16.37 | 17.86 | — | — | — | — |
| batch 74sp | 102.48 | 53.95 | 110.64 | 67.57 | — | — | — | — |

**方向一致性**（%Δ vs 各路线 BAU）：
```
fac_v3:  LNNOX −46.9✓  HVOC +9.5✓  ISOPX −28.0✓  HOT +10.8✓  CLOUD −17.0✓  IC×2 +0.0✓  BG+50% +2.0✓
mcm_run: LNNOX −9.6✓   HVOC +0.9✓  ISOPX +10.1✗  ← 跨机制真实分歧（MCM 谱下 ISOP 增 O₃）
batch:   LNNOX −47.4✓  HVOC +8.0✓  ISOPX −34.1✓
```

**能力边界**：
| ✅ 能 | ❌ 不能 |
|---|---|
| 命名情景一键切换（改 config 不动代码） | 分时段情景（如"只限行早高峰"）——因子是全程常数 |
| 4 排放组 + 光解 + 温度 组合 | 单独调某个 VOC（如"只减芳香烃"）——只能整组缩放 |
| 多路线方向一致性检验（mcm 的 ✗ 是真实分歧） | 自发突变/反馈情景（需外部控制回路） |

**扩展方式**：`scenarios.config` 加一行即可，如
`SCEN_MINE = E_NOx*0.7 E_CO*0.5 E_VOC*1 E_CH4*1 JFAC*0.9 DT+2`

#### 6.11.2 过程分析（Process Analysis）—— ✅ 有，是 fac_v3 的核心差异化能力

**两套账并行（独立累积、可互证）**：

| 账本 | 名称 | 粒度 | 内容 | 输出 |
|---|---|---|---|---|
| **IPR** | Integrated Process Rate（过程速率积分） | **过程级**（6 类） | ①机制P ②机制L ③排放 ④干+湿沉 ⑤BC+稀释 ⑥BLH膨胀 | `<pref>_ipr.csv`（74 物种 × 6 账） |
| **IRR** | Integrated Reaction Rate（反应速率积分） | **反应级**（131 条） | 逐反应速率的小时窗精确积分 | `<pref>_irr.csv`（r1..r131） |

**关键设计：IPR 是"小时窗精确积分"，不是快照**
- `fac_v3` 在内步（CVODE 内部子步）累积，输出每小时窗口的积分量（molec/cm³）
- 这是为什么它能做到 **IPR↔IRR 双路径 0.000% 一致**（`hourly_pa.py` 实测）
- 对比：mcm_run 是小时快照，大项对消（P/L 各 ~1e13 → 净值 ~1e10）时有混叠

**实测闭合质量**（`ipr_closure_check.py`）：

| 物种 | 中位比值（六账和 / ΔC） | 判定 |
|---|---|---|
| H₂O₂ | 1.0190 | ✓ |
| CO | 1.0165 | ✓ |
| PAN | 1.0164 | ✓ |
| ISOP | 1.0148 | ✓ |
| CH₃O₂ | 1.0129 | ✓ |
| OH | 1.0155 | ✓ |

> 比值 1.01–1.02 而非严格 1.000 → 梯形积分对小时窗内非线性变化的残差，约 1–2%，可接受。

**三个过程分析工具**：

| 工具 | 输出 | 用途 |
|---|---|---|
| **`hourly_pa.py`** | ①逐小时六账占比（柱状/堆叠）②O₃ 逐小时反应级 top-4 预算 ③IPR↔IRR 一致性检验 | 看"什么过程在什么时候主导" |
| **`ipr_irr5.py`** | 多路线（fac_v3/mcm/batch）IPR 六账占比横向对比 | 看"不同数值路线是否给同一过程结论" |
| **`diagnostics.ipr_from_fac()`** | 任意前缀的 74 物种 × 6 账积分 | 编程式自定义分析 |

**典型产出示例**（O₃ 的六账）：
```
O3 六账: P=+4.63e13  L=−4.16e13  E=0.00  D=−2.48e12  B=−2.26e12  X=0.00
       → 净变化 = +4.5e10 molec/cm³ (≈1.8 ppb)
       → 化学净产占主导，沉降消耗次之，无排放（O₃ 是二次产物）
```

**能力边界**：
| ✅ 能 | ❌ 不能 |
|---|---|
| 每物种每过程的贡献分解（六账） | 空间过程的分解（无空间维） |
| 逐小时演变（何时主导） | 单网格点的过程分析（0-D 箱整体） |
| 反应级 top 预算（131 反应排名） | 单分子/单事件的轨迹追踪 |
| 跨路线过程结论一致性 | 非均相/液相过程的独立账（目前并入机制 P/L） |

#### 6.11.3 四项分析与情景/过程的关系

```
响应面 (Where?)      fac_batch_omp    全局地图: 屋脊在哪 / 体制在哪转折
    │
    ├→ 情景分析 (What if?)   scenarios.config   在关键区试政策组合
    │       └→ 过程分析 (Why?)   IPR + IRR        为什么变 (六账/131 反应)
    │              └→ 敏感性 OAT (How much?)  sens_adj.py   精确梯度 ∂y/∂p
    │
    └→ 反演/伴随 (推断)   FD+SVD   从观测反推排放
```

**四者的依赖链**：**响应面定地图 → 情景试政策 → 过程查原因 → 敏感性求梯度**。

> ⚠ 注意 **响应面与 OAT 不是同一层的替代品**：前者是全局离散地图，后者是
> 单点局部梯度（详见 §6.10.0 对照表）。

#### 6.11.4 六项能力总表（含此前的四项）

| # | 能力 | 类型 | 状态 | 入口 | 核心产物 |
|---|---|---|---|---|---|
| 1 | **观测约束** | 数据同化 | ✅ 生产级 | `bc.config` + `TAU_MIX` | IPR 的 B 账 |
| 2 | **响应面** | **全局扫描** | ✅ 4D 参数空间 | `fac_batch_omp` | 24 列 CSV → 114 曲面族 |
| 3 | **情景分析** | **离散场景对比** | ✅ 政策沙盘 | `scenarios.config` + `scen5.py` | 情景矩阵 + 方向一致性 |
| 4 | **过程分析** | **过程归因** | ✅ **核心差异化** | `hourly_pa.py` + IPR/IRR | 六账 + 131 反应优先级 |
| 5 | **敏感性 OAT** | **局部梯度** | ✅ | `sens_adj.py` | `jacobian.csv`（24×9） |
| 6 | **伴随/反演** | 推断（非真伴随） | ⚠ 仅 1/9 模态可辨识 | `sens_adj.py` ② | FD 雅可比 + 截断 SVD |

**类型分组的含义**：
- **全局扫描**（响应面）：铺满参数空间，看形状/转折
- **离散场景**（情景分析）：几个特定点，看政策响应方向
- **过程归因**（过程分析）：给定一个点，看内部机制
- **局部梯度**（OAT）：给定一个点，看微小扰动响应
- **推断**（反演）：由观测反推参数
- **数据同化**（观测约束）：把观测强行拉进模型

这六类**互为补充而非替代**——各自回答不同层次的问题。

## 7. 实验 C：屋脊的跨路线对比

### 7.1 科学问题

"屋脊—指示因子对应关系"是普适的，还是**依赖机制与排放谱**？

### 7.2 前置条件

- 实验 B 已产出 `output/batch_v3_indicators.csv`
- `mcm_city/mcm_run` 已构建

### 7.3 逐步命令

```bash
python3 experiments/ridge5.py            # fnox ∈ logspace(-1,1,13)，现役路线
python3 experiments/ridge5.py --fast     # 快速版
python3 experiments/fig_ridge5.py        # → figures/ridge5_routes.png
```

### 7.4 V2.3 实测与判读

| 路线 | 屋脊 fnox | 指示因子一致性 | 说明 |
|---|---|---|---|
| batch（74sp 排放空间） | **2.28** | FNR 过渡 0.67 **贴屋脊** | 主结论 |
| mcm（29sp 排放空间） | 0.14 | OPE 谷精确命中（0.00） | 其排放谱下屋脊位置不同 |
| kpp（纯机制，已退役） | — | IC 空间无真屋脊 | 无汇体系单调，路线定位而非缺陷 |

**教学点**：对应关系**依赖 VOC 反应性谱**——这正是 FNR 阈值需按机制/区域
重标定的原因（Sillman 原值 1–2 是特定机制下的标定）。

### 7.5 变异玩法

- 换 fnox 范围/点数：改 `FN = np.logspace(-1,1,13)`
- 只跑两路线：注释掉其它 `run_*`
- 加自己的机制：`route_indicators(pref, 你的.fac)`

### 7.6 常见坑

- 别拿两路线的屋脊 fnox 直接比大小——**排放谱不同，位置本应不同**；
  有意义的是"各自屋脊 vs 各自指标过渡点"的一致性
- 历史"五路线"口径中 fac_run/kpp 已退役，脚本内有守卫与注记

---

## 8. 实验 D：IPR/IRR 过程诊断

### 8.1 科学问题

O₃ 的变化**由谁主导**？拆到逐反应（IRR）与逐过程（IPR）。

### 8.2 前置条件

- 实验 A 已有 `output/expA_BAU_{conc,irr,ipr}.csv`

### 8.3 逐步命令

```bash
# D-1 闭合校验（六账积分 vs ΔC，关键物种）
cd experiments && python3 ipr_closure_check.py
# 期望: ISOP 中位比值 1.014 ✓ / CH3O2 1.013 ✓ / OH 1.016 ✓

# D-2 多路线过程对比
python3 ipr_irr5.py

# D-3 单案例手工读账（Python API）
cd ..
python3 - <<'PY'
import sys; sys.path.insert(0, "diagnostics")
import diagnostics as dg

h, C, ih, I = dg.load("output/expA_BAU")       # (25,),(25,74),(25,),(25,131)
print("窗口行数", len(h), "物种数", C.shape[1], "反应数", I.shape[1])

sp, ipr = dg.ipr_from_fac("output/expA_BAU")   # 任意路线通用
print("O3 六账:", {k: f"{v:.2e}" for k, v in ipr["O3"].items()})
PY
```

### 8.4 预期输出（实测）

```
窗口行数 25 物种数 74 反应数 131
O3 六账: {'P': '4.63e+13', 'L': '-4.16e+13', 'E': '0.00e+00',
          'D': '-2.48e+12', 'B': '-2.26e+12', 'X': '0.00e+00'}
```
**怎么读**：净变化 = P+L+D+B+X = 4.63−4.16−0.248−0.226 ≈ **+4.5e10 molec/cm³**
（≈1.8 ppb 净增），与末 24h O₃ 变化一致 → 账闭合。
`E=0`（O₃ 无排放）与 `X=0`（该窗口 BLH 未增长）都是合理结果。

### 8.5 IRR 关键通道（自己算）

```python
import numpy as np
from diagnostics import load
h, C, ih, I = load("output/expA_BAU")
I_int = np.trapz(I, ih, axis=0)          # 逐反应窗口积分

# 找 HO2+NO 的反应号（.fac 的 % 行序）
mech = [l for l in open("mechanism/mechanism_3_full_chain.fac") if l.startswith("%")]
for k, l in enumerate(mech, 1):
    if "HO2 + NO = NO2 + OH" in l:
        print(f"HO2+NO 是第 {k} 个反应, 窗口积分 = {I_int[k-1]:.3e}")
```

### 8.6 判读要点（含必读边界）

| 现象 | 解读 |
|---|---|
| P 与 L 都很大、净值很小 | **大项对消**——白天 P/L 各 ~1e13 对消到净值 ~1e10；快照采样对此有混叠 |
| fac_v3 的闭合比值 1.004–1.014 | 内步精确累积，可信 |
| mcm 的比值略差 | 小时快照采样，混叠更大（路线特性，非 bug） |
| OH 的六账占比 | OH 寿命 ~1 s ≪ 1 h 窗口，**占比仅指示性** |

### 8.7 变异玩法

- 换物种：`ipr_from_fac` 返回 74 物种全账，任取
- 换窗口：改 `model.parameters` 的天数/步长
- 逐时看账：直接读 `_ipr.csv` 的原始行（实验 G）

### 8.8 常见坑

- ❌ 用 `diagnostics.report()` —— **该函数已移除**（旧文档遗留），
  现在用 `dg.load` + `dg.ipr_from_fac` + `indicators.route_indicators`
- ❌ 把 P/L 绝对值当生成速率——窗口积分单位是 molec/cm³（累积量）
- ❌ 用 OH 的账占比下强结论

---

## 9. 实验 E：指示因子工具箱

### 9.1 科学问题

任意一条路线的输出，**12 个指示因子一键算出**（路线无关、机制无关）。

### 9.2 前置条件

- 任意 `_conc.csv` + `_irr.csv` 对（fac_v3 / batch / mcm 皆可）
- 机制 .fac 带 FORMULA 段（NOy 判定用）

### 9.3 逐步命令

```python
import sys; sys.path.insert(0, "diagnostics")
from indicators import route_indicators, route_timeseries, Mechanism

# ① 末 24h 窗口均值 → 12 指标
d = route_indicators("output/expA_BAU", "mechanism/mechanism_3_full_chain.fac")
for k, v in d.items(): print(f"{k:14s} {v:12.4g}")

# ② 正午窗（12–18h，Sillman 口径）
d_pm = route_indicators("output/expA_BAU", "mechanism/mechanism_3_full_chain.fac",
                        pm_only=True)

# ③ 逐时序列（画日循环）
ts = route_timeseries("output/expA_BAU", "mechanism/mechanism_3_full_chain.fac")
# 键: FNR, H2O2_HNO3, O3_NOy, HCHO_NOy, O3_NOz, NOy_ppb,
#      Q, L_NOx, L_NOx_Q, P_H2O2_HNO3, OPE, O3_ppb

# ④ 复用机制对象（批量算多个案例时省解析）
m = Mechanism("mechanism/mechanism_3_full_chain.fac")
d2 = route_indicators("output/expA_LNNOX", "mechanism/mechanism_3_full_chain.fac", mech=m)
```

### 9.4 预期输出（BAU 实测）

| 指标 | 实测值 | 含义 |
|---|---|---|
| `FNR` | **3.427** | HCHO/NO₂ >2 → NOx-limited |
| `H2O2_HNO3` | 0.227 | >0.2 → NOx 控制区 |
| `O3_NOy` | 5.211 | O₃ 生成效率 |
| `HCHO_NOy` | 0.662 | — |
| `O3_NOz` | 6.600 | — |
| `NOy_ppb` | 10.692 | NOx 存量（17 物种） |
| `Q` | 3.015e10 | 自由基源总速率 |
| `L_NOx` | 9.879e9 | NOx→储库速率 |
| `L_NOx_Q` | **2.011** | >0.5 → 见 §16.3 冲突讨论 |
| `P_H2O2_HNO3` | 0.115 | — |
| `OPE` | 5.319 | 每分子 NOx 造 5.3 个 O₃ |
| `O3_ppb` | 58.429 | 末 24h 均值 |

### 9.5 判读要点

- **FNR 与 L_NOx/Q 在本配置下方向不一致**（3.43 说 NOx-limited，
  2.01 说 NOx 饱和）——真实案例，见 §16.3（HONO 计入 L_NOx 的定义敏感性）。
  **多指标交叉 + 屋脊几何**才是可靠做法。
- 12 指标口径差异：`FNR` 是浓度比，`L_NOx_Q` 是速率比，`OPE` 是过程量比——
  **不可混用同一阈值**。

### 9.6 变异玩法

- 输出到 CSV 做多情景对比表
- `ts` 画日循环：`matplotlib` 画 `FNR` vs `h`
- 用自己的机制：`.fac` 带 FORMULA 段即可（引擎零改动）

### 9.7 常见坑

- ❌ 用旧 API `dg.indicators(C)` / `dg.report()` —— 已移除
- ❌ 把 `pm_only=True` 当"更适合"——它是 Sillman 原始口径，两者可对照看
- ❌ NOy 未含全部 17 物种 → 自动从 FORMULA 提取，别手工加

---

## 10. 实验 F：观测约束 OSSE

### 10.1 科学问题

给定观测，**同化能把模型拉多准？观测算子误差会怎样放大？**

### 10.2 前置条件

- fac_v3 + fac_batch_omp（F3 序列 BC 同化通路）
- 理解 §3.2 ⑤ bc.config 的 TAU_MIX

### 10.3 逐步命令

```bash
python3 experiments/obs_osse.py        # 全自动: 真值→采样→同化→技能评估
```

### 10.4 V2.3 实测与判读

| 实验设置 | 结果 |
|---|---|
| 地面同化 HCHO（τ=2h） | RMSE **+27%** |
| 卫星反演 E_NOx×1.34 | NO₂ RMSE +74%、O₃ +58% |
| B 账同化倾向 | O₃ −1.9e7（占过程总账 ~2%） |

**判读**：
1. 张弛 τ 决定技能上限——τ 太大（弱约束）几乎无效，τ 太小（强约束）压制模型化学
2. **观测算子误差会被放大**（NO₂ +74% >> E_NOx 误差 34%）——为什么卫星反演要谨慎
3. 同化作用可在 `_ipr.csv` 的 `B` 账中定量核算

### 10.5 变异玩法

| 想研究 | 怎么做 |
|---|---|
| τ 敏感性 | 改 bc.config 的 `TAU_MIX`（7200/21600/86400）重跑 |
| 采样频率 | 改约束序列点数（≤128） |
| 观测组合 | bc.config 加多物种 file:（注意总池 ≤8） |
| 自己造真值 | 用另一情景的输出当"观测"做 OSSE |

### 10.6 常见坑

- 物种数 × 序列数 > 8 → 编译拒绝（气象占 2 条时只剩 6 条）
- 单位忘记 `# unit: ppb`
- τ 设成 0 → 除零

---

## 11. 实验 G：逐小时过程分析

### 11.1 科学问题

fac_v3 的 IPR/IRR 是**小时窗精确积分**——逐时看每个过程何时主导，
并自证"IPR↔IRR 双累积路径一致"。

### 11.2 逐步命令

```bash
python3 experiments/hourly_pa.py                        # 缺省现跑 SCEN_BAU
python3 experiments/hourly_pa.py output/expA_LNNOX      # 用已有输出
```

### 11.3 V2.3 实测

```
最大相对差 = 0.000%  (✓ 两累积路径一致)
图: figures/hourly_pa.png
```

### 11.4 判读要点

逐时看三件事：
1. **O₃ 净生成何时转正/转负**：日出后 1–2 h 转正，傍晚转负
2. **BC 同化倾向何时主导**：观测偏离模型时（`B` 账突增）
3. **BLH 膨胀稀释何时抬升**：上午 8–10 点（`X` 账突增）

### 11.5 常见坑

- ❌ 用 mcm 的输出做此实验——快照采样有混叠；过程分析固定用 fac_v3
- ❌ 以为 0.000% 是巧合——它是 IPR 与 IRR 两条独立积分路径的相互验证

---

## 12. 实验 H：灵敏度/伴随分析

### 12.1 科学问题

**∂O₃/∂(排放)** 有多大？哪些反应速率最影响 O₃？

### 12.2 逐步命令

```bash
python3 experiments/sens_adj.py        # OAT 雅可比 + 截断 SVD (~2 min)
```

### 12.3 V2.3 实测与判读

```
fac_v3 OAT (NO+NO₂ ±10%): 42.6 ppb / unit fnox
batch 网格局部斜率:        52.3 ppb / unit fnox
一致性: 19% 偏差
```
- **19% 偏差是方法差异**（OAT 小扰动 vs 网格大步长），非错误
- OAT 给出**局部**斜率；跨越非线性区（如屋脊附近）时必须用网格
- ④ 真速率灵敏度（∂O₃/∂kᵢ 逐反应）需 fac_v3 重编译载体，随 fac_run 退役待重建

### 12.4 变异玩法

- 改扰动幅度：`sens_adj.py [exe_path] [delta_pct]`（如 5、20）
- 加扰动变量：改脚本中的 OAT 表
- 用 SVD 结果做参数辨识（脚本已含截断 SVD 反演框架）

### 12.5 常见坑

- 把 OAT 斜率当全局——非线性区会失真
- 扰动过大（>50%）跨越体制边界 → 结果无意义

---

## 13. 实验 I：情景矩阵

### 13.1 科学问题

同一政策情景，**不同路线的响应方向是否一致**？（跨机制稳健性）

### 13.2 逐步命令

```bash
python3 experiments/scen5.py       # 路线×情景 O₃max 响应矩阵 (%Δ) + 方向一致性
```

### 13.3 V2.3 实测与判读

- fac_v3 / batch：七情景方向**全 ✓**
- **mcm ISOPX 方向相反（✗）**：MCM 速率谱下 ISOP×2 增 O₃，ozone_lab 谱下减

**这个分歧是本次实验最有价值的结果**——它证明：
1. 跨机制结论不能照搬
2. 指示因子阈值需按机制重标定
3. "ISOP 增加对 O₃ 的影响"依赖 NOx 水平与氧化路径分配（MCM 显式异戊二烯
   化学 vs 集总模型）

### 13.4 常见坑

- 把方向不一致当 bug——它是真实的机制差异
- 用 mcm 的绝对 O₃ 与 fac_v3 比大小——机制不同，只能比方向

---

## 14. 六类可选输入数据文件接入操作（templates/ 逐一详解）

> **核心原则（源码实证）**：`templates/*.template` **不被任何程序读取**——它是
> **格式合同**：告诉你原始数据长什么样、怎么换算、落到哪个路径、去哪个 config 登记。
> 真正被程序消费的只有 `input_v3/constraints/*.txt` 与 12 类 config。
>
> **复现验证**（两条命令，输出均为空）：
> ```bash
> grep -rn "templates" --include="*.py" --include="*.f90" --include="*.sh" \
>      src/ modules/ facsimile_impl/ tests/ experiments/
> grep -rn "\.template" src/ modules/ facsimile_impl/ tests/ experiments/
> ```
> 逐 `open()` 点的完整清单见 [§3.1 源码级验证](#31-五段管线程序只认-runtime_v3dat)。
>
> **通用四步**（五类模板全部适用）：
> ```
> ① 按模板格式准备数据  →  ② 落位到 input_v3/constraints/（或直接编辑 emissions.config）
> ③ 在 bc.config / meteo.config 登记一行 file:  →  ④ make runtime_v3.dat 重编译
> ```
> **单位总参考**：[UNITS.md](UNITS.md)（449 行逐子量表 + 换算公式 + 硬编码登记）。

### 14.0 六类模板总览

| # | 模板 | 你提供什么 | 落位到 | 登记到 | 引擎直读? |
|---|---|---|---|---|---|
| ① | `surface_o3_station.template` | 站点 O₃ 序列（ppb） | `constraints/o3_station_<站名>.txt` | `bc.config` | ✓ |
| ② | `surface_aqi_station.template` | AQI 宽表（μg/m³） | **先离线拆分换算**成 ① 格式 | `bc.config` | ✗（需拆分） |
| ③ | `surface_voc_station.template` | VOC 物种序列（一站一物种） | `constraints/voc_<物种>_<站名>.txt` | `bc.config` | ✓ |
| ④ | `satellite_column.template` | 卫星柱浓度（molec/cm²） | **先离线转表面等效 ppb** | `bc.config` | ✗（0-D 无气柱） |
| ⑤ | `emissions_inventory.template` | 排放清单 | **直接编辑** `emissions.config` | —（本身就是配置） | ✓ |
| ⑥ | **`met_series.template`** | **气象逐时序列（K/hPa/%/m/mm·h⁻¹/1）** | `constraints/met_<量>_<站名>.txt` | **`meteo.config`** | ✓ |

**序列文件的通用约束**（违反即 cfgbuild 编译期拒绝）：

| 约束 | 值 | 错误信息关键字 |
|---|---|---|
| 每序列点数 | ≤ 128 | `序列 N 点 > Fortran 端上限 128` |
| 全池序列数 | ≤ 8（气象+BC 共享） | `序列总数 9 > Fortran 端上限 8` |
| 时间轴 | 自模拟起点秒、严格递增 | `时间非单调 @ t=…` |
| 单位指令 | 浓度序列: `# unit: cm3\|ppb\|ppm`（缺省 cm3）<br>**气象序列: 禁止出现** | `未知单位 'x'` / `气象序列不接受单位指令行` |
| 物种名 | ∈ 机制 VARIABLE 段 | `未知物种` |

---

### 14.0a 完整数据链路：模板 → 配置 → runtime.dat（四层，不是两层）

⚠ **常见误解**：以为"模板只是语法参考，只有 runtime.dat 才起效"——这**漏掉了中间层**。
准确说是**四层**，每一层都有明确角色：

```
第 1 层  templates/*.template          ← 格式合同（人读，程序零读）
           │ 你照它的格式，写自己的数据
           ▼
第 2 层  input_v3/*.config             ← ★ 真正的数据源（你编辑这里）
         input_v3/constraints/*.txt        排放/IC/气象/BC/情景 + 序列数据
           │ src/cfgbuild.py 编译（校验+单位换算+归一）
           ▼
第 3 层  runtime_v3.dat                ← 编译产物（自动生成，勿手改）
           │ 唯一被 Fortran 读的文件
           ▼
第 4 层  fac_v3 / fac_batch_omp        ← 运行时（零文本解析）
```

#### 为什么不是"两层"？

| 你的说法 | 实际情况 |
|---|---|
| "模板只是学习语法" | ✅ 部分对：程序**确实不读** templates（`grep` 零引用）；但它除了语法还带**可直接复制的示例数据** |
| "只有 runtime.dat 才起效" | ✅ 对：Fortran **只读** `runtime_v3.dat`（`runtime_dat_mod.f90:72`） |
| **漏掉的** | ❌ **`input_v3/*.config` 才是真正的数据源**——`runtime.dat` 是它的**编译产物**，你不直接编辑 runtime.dat |

#### 实测验证（数据源确认）

```bash
# 只改 input_v3/meteo.config 一个数（TEMP 298.15 → 300.15）
md5 runtime_v3.dat     # 52d75f91... （改前）
# 重新编译
make runtime_v3.dat
md5 runtime_v3.dat     # 673f9065... （改后，已变化 ✓）
```
→ 证明 `runtime.dat` 完全由 `input_v3/` 派生；**改 `input_v3/` 才是"改数据"**。

#### 三层的准确定位

| 层 | 你会不会改 | 作用 |
|---|---|---|
| `templates/*.template` | **只读**（照它抄格式） | 告诉你"文件该长什么样" |
| `input_v3/*.config` + `constraints/*.txt` | **✓ 改这里** | 真正的数据（排放/气象/观测/情景） |
| `runtime_v3.dat` | **✗ 不手改** | 自动生成的中间产物，Fortran 唯一输入 |

#### 一个容易忽略的细节：`constraints/*.txt` 是**间接**读取的

你并不直接"登记"约束文件——而是在 config 里写 `file:` 引用：

```
constraints/bc_O3_demo.txt          ← 序列数据本体
        ▲
        │ 被引用
meteo.config:  BLH = file:constraints/met_BLH.txt     ← 气象序列
bc.config:     O3  = file:constraints/bc_O3_demo.txt  ← 观测约束
        │
        │ cfgbuild 顺着 file: 去读（read_series）
        ▼
runtime_v3.dat 的 SERIES 段（序列被内联成表）
```

代码位置：`cfgbuild.py:178`（meteo 的 file:）、`:197`（bc 的 file:）。

**所以完整表述应该是**：
> **模板告诉你格式 → 你在 `input_v3/` 按格式写数据（config 里用 `file:` 引用 constraints 序列）
> → cfgbuild 编译成 `runtime_v3.dat` → Fortran 只读 runtime.dat 运行。**
>
> 改数据 = 改 `input_v3/` + `make runtime_v3.dat`；
> 改 `runtime.dat` 本身没有意义（会被下次编译覆盖）。

#### ⚠ mtime 陷阱（本项目真实踩过，已修复）

**症状**：改了 config 又改回，`make runtime_v3.dat` 说 `up to date`，但
`runtime_v3.dat` 里其实是**中间那次**的值 → 模拟结果**无声出错**。

**原因**：`make` 用文件 mtime 判断依赖。若两次编辑落在**同一秒**（`sed -i` 连续操作很常见），
mtime 相同 → make 认为不需要重建。

**实测复现**：
```bash
# 原值 TEMP=298.15 → 改成 299.15 → 再改回 298.15，三次都在同一秒内
make runtime_v3.dat       # → "make: 'runtime_v3.dat' is up to date."
md5 runtime_v3.dat        # → 仍是 299.15 那次的版本（错！）
```

**已修复**（Makefile）：
```make
runtime_v3.dat: FORCE $(wildcard input_v3/* ...)
	@$(PY) src/cfgbuild.py input_v3 runtime_v3.dat --mech ...
FORCE:
```
`FORCE` 让 `runtime_v3.dat` **每次 make 都重建**——cfgbuild 只需 0.37 s 且幂等，
强制重建是更安全的默认。

**已加防线**（`tests/functional_test.py` 的 **FT-CFG**）：
```
用 cfgbuild 重新生成一份到临时文件，与仓库内 runtime_v3.dat 逐字节比对
→ 不一致则 FAIL，提示 "⚠ 仓库内 runtime_v3.dat 陈旧, 请 make runtime_v3.dat"
```

**给用户的纪律**：
- 改完任何 `input_v3/` 内容 → **`make runtime_v3.dat`**（现在总会真的重建）
- 手动编译过 Fortran 后 → 确认 `runtime_v3.dat` 是最新的
- 怀疑数据没生效 → 跑 `python tests/functional_test.py`，看 FT-CFG 是否 PASS

---

#### 📋 实验操作 1：逐层验证四层链路（完整可复现）

**目的**：亲手证明"改 `input_v3/` → `runtime_v3.dat` 变化 → 运行结果变化"这条链路，
以及"templates 不参与"。

```bash
cd ozone_lab

# ═══ L1 验证：templates 程序零读 ═══
grep -rn "templates" --include="*.py" --include="*.f90" --include="*.sh" \
     --include="Makefile" src/ modules/ facsimile_impl/ tests/ experiments/
# 期望输出: 空（或仅 sat2surf.py 的注释文字）
# → 证明 templates/*.template 从不被程序打开

# ═══ L2 → L3 验证：改数据源，产物随之变化 ═══
md5 -q runtime_v3.dat                    # 记下改前 md5: 52d75f91...
sed -i '' 's/^TEMP   = 298.15/TEMP   = 300.15/' input_v3/meteo.config
make runtime_v3.dat                      # 重新编译（0.37 s）
md5 -q runtime_v3.dat                    # 改后 md5: 673f9065... ← 已变化 ✓
grep -A1 "^MET" runtime_v3.dat           # 直接看产物里的值
# 期望: MET / C 300.15        ← TEMP 已写入 runtime.dat

# ═══ L3 → L4 验证：运行时确实用新值 ═══
cd facsimile_impl
./fac_v3 SCEN_BAU ../output/l4_test ../runtime_v3.dat
# 期望: O3max ≈ 106.2 ppb（基线 TEMP=298.15 时是 101.61 → 温度升高使 O₃ 上升）
cd ..

# ═══ 还原 ═══
sed -i '' 's/^TEMP   = 300.15/TEMP   = 298.15/' input_v3/meteo.config
make runtime_v3.dat
md5 -q runtime_v3.dat                    # 回到 52d75f91... ✓
```

**实测记录（2026-10-26）**：

| 步骤 | 观察值 | 说明 |
|---|---|---|
| L1 grep | 空输出 | templates 确无程序读取 |
| L2 改前 md5 | `52d75f91...` | 基线 |
| L2 改后 md5 | `673f9065...` | 数据源变更传导到产物 |
| L3 产物内容 | `MET` / `C 300.15` | TEMP 已编译进 runtime.dat |
| L4 运行 | `O3max = 106.24 ppb` | 比基线 101.61 高（温度效应） |
| 还原后 md5 | `52d75f91...` | 完全复原 |

**结论**：`runtime.dat` 是 `input_v3/` 的纯函数；`templates` 只在人脑里起作用。

---

#### 📋 实验操作 2：复现并验证 mtime 陷阱（修复前/后对比）

**目的**：亲手看到"同秒改又改回 → make 不重建 → 无声出错"，再验证 FORCE 修复。

**（A）在有 FORCE 修复的当前版本上跑——应正确重建**

```bash
cd ozone_lab
md5_ok=$(md5 -q runtime_v3.dat)          # 记录正确值
echo "正确值: $md5_ok"

# 同秒内：改 → make → 改回 → make
sed -i '' 's/^TEMP   = 298.15/TEMP   = 299.15/' input_v3/meteo.config
make runtime_v3.dat                       # 跑一次（值 → 299.15）
sed -i '' 's/^TEMP   = 299.15/TEMP   = 298.15/' input_v3/meteo.config
make runtime_v3.dat                       # 再次跑（值应回到 298.15）

md5_now=$(md5 -q runtime_v3.dat)
[ "$md5_now" = "$md5_ok" ] \
  && echo "✓ FORCE 生效：正确重建回原值" \
  || echo "✗ 仍不一致（FORCE 未生效？）"
```

**实测**：`✓ 正确重建为原值`（两次都运行了 cfgbuild）

**（B）如果要看未修复时的表现**（了解原理，勿在正式环境做）

把 Makefile 的 `runtime_v3.dat: FORCE ...` 临时改回 `runtime_v3.dat: ...`（去掉 FORCE），
重复上面流程，会得到：

```
make: `runtime_v3.dat' is up to date.        ← 第二次没重建！
md5 → 14f7a885...（仍是 299.15 那版）        ← 无声出错
```

**为什么同秒会发生**：`sed -i` 两次写文件通常在同一秒内；`make` 只比较
"源 mtime > 目标 mtime"，同秒时不大于 → 判为最新。

**（C）防线验证：让 FT-CFG 检测陈旧 runtime.dat**

```bash
# 故意制造陈旧：手工把 runtime_v3.dat 改成过期内容
cp runtime_v3.dat /tmp/rt_backup.dat
sed -i '' 's/^  C 298.15/  C 999.99/' runtime_v3.dat    # 篡改产物
python3 tests/functional_test.py 2>&1 | grep -A1 "FT-CFG"
# 期望: [FAIL] runtime_v3.dat 与 input_v3 当前内容一致 (非陈旧)
#        ⚠ 仓库内 runtime_v3.dat 陈旧, 请 make runtime_v3.dat
# 恢复:
cp /tmp/rt_backup.dat runtime_v3.dat && rm /tmp/rt_backup.dat
```

---

#### 📋 实验操作 3：完整工作流（改数据 → 看结果）

**这是日常最高频的操作**，三步走：

```bash
cd ozone_lab

# ═══ 步骤 1：改数据（改 input_v3/，不是 templates/ 也不是 runtime.dat）═══
# 例：把 NO 排放基线调成 2 倍
sed -i '' 's/^NO       2.157e+07/NO       4.314e+07/' input_v3/emissions.config
#    或改背景 O₃（bc 常数用法）
sed -i '' 's/^O3  = 6.2e11/O3  = 9.84e11/' input_v3/bc.config   # 25→40 ppb
#    或加观测约束（bc 序列用法）
echo "O3 = file:constraints/bc_O3_demo.txt" >> input_v3/bc.config

# ═══ 步骤 2：重编译数据（现在总会真的重建）═══
make runtime_v3.dat
# 期望: "情景: ['SCEN_BAU', ...]" → 证明本次真的运行了

# ═══ 步骤 3：运行 + 看结果 ═══
cd facsimile_impl && ./fac_v3 SCEN_BAU ../output/mytest ../runtime_v3.dat && cd ..
# 期望: "FAC-V3 ok (…s, nst=…, O3max=… ppb)"

# ═══ 可选：验证结果合理性 ═══
python3 -c "
import sys; sys.path.insert(0,'diagnostics')
from indicators import route_indicators
d = route_indicators('output/mytest', 'mechanism/mechanism_3_full_chain.fac')
print(f\"O3={d['O3_ppb']:.1f} ppb  FNR={d['FNR']:.2f}\")"
```

**常见错误对照**：

| 你改的是 | 结果 | 正确做法 |
|---|---|---|
| `templates/*.template` | **完全没效果**（程序不读） | 改成 `input_v3/` 对应文件 |
| `runtime_v3.dat` | 下次编译被覆盖 | 改 `input_v3/` |
| `input_v3/` 但不 make | 跑的还是旧数据 | **`make runtime_v3.dat`** |
| `input_v3/` + make | ✅ 正确 | — |

---

### 14.0b `input_v3/constraints/` 里的现成文件是什么（重要区分）

⚠ **该目录有两种性质完全不同的文件**，格式相同但用途相反：

| 文件 | 性质 | 被谁引用 | 物理含义 |
|---|---|---|---|
| `met_BLH.txt` / `met_RAIN.txt` | **气象序列**（默认启用，7 点） | `meteo.config` 的 `BLH/RAIN = file:...` | 边界层高度 / 降水（物理环境） |
| **`met_{TEMP,PRES,RH,BLH,RAIN,JFAC}_demo.txt`** | **气象序列样例**（各 121 点逐时） | 同上（6 变量全序列） | 实验 M 的驱动数据 |
| **`bc_O3_demo.txt`** | **观测约束样例** | `bc.config` 的 `O3 = file:...` | 站点 O₃ 观测（真值校准） |
| **`bc_HCHO_demo.txt`** | **观测约束样例** | `bc.config` 的 `HCHO = file:...` | HCHO 观测（FNR 分子） |
| **`bc_NO2_sat_demo.txt`** | **观测约束样例** | `bc.config` 的 `NO2 = file:...` | 卫星反演 NO₂（弱 τ 专用） |

**判断口诀**：
- 前缀 **`met_`** → 气象，登记到 `meteo.config`
- 前缀 **`bc_`** → 观测约束，登记到 `bc.config`

> `met_BLH.txt` **不能**当成 `bc.config` 的参考样式——它是边界层高度（米），
> 不是浓度（ppb）。虽然格式都是"两列 + 数值"，但单位和物理含义完全不同。

**三个 `bc_*_demo.txt` 已实测可用**（2026-10-26）：

| 样例 | 登记行 | 推荐 τ | 实测 B 账 |
|---|---|---|---|
| `bc_O3_demo.txt` | `O3 = file:constraints/bc_O3_demo.txt` | 21600 | **−3.13e12** ✓ |
| `bc_HCHO_demo.txt` | `HCHO = file:constraints/bc_HCHO_demo.txt` | 21600 | **−4.03e11** ✓ |
| `bc_NO2_sat_demo.txt` | `NO2 = file:constraints/bc_NO2_sat_demo.txt` | **43200**（必须弱） | **+3.02e11** ✓ |

**气象序列的最小试跑**（对照实验 M 会自动做这件事）：
```bash
# meteo.config 里把 6 行改成 file:（文件末尾已备好注释块）
# 然后直接跑：常数基准 / 全序列 / 六个单变量，一次出 10 个 run
python experiments/run/run_M_met_series.py
# 看归因：output/run_M/attribution.csv
```

**观测约束的最小试跑**（以 O₃ 为例，3 条命令）：
```bash
# ① bc.config 加一行（保留原有背景常数行）
echo "O3 = file:constraints/bc_O3_demo.txt" >> input_v3/bc.config
# ② 重编译
make runtime_v3.dat
# ③ 跑 + 查同化账
cd facsimile_impl && ./fac_v3 SCEN_BAU ../output/assim_demo ../runtime_v3.dat && cd ..
python3 -c "
import sys; sys.path.insert(0,'diagnostics'); import diagnostics as dg
sp, ipr = dg.ipr_from_fac('output/assim_demo')
print('O3 的 B账(同化) =', f\"{ipr['O3']['B']:+.3e}\", '→ 非零即生效')"
```

---

### 14.1 模板①：地面 O₃ 站点序列（最常用，引擎直读）

**模板**：`templates/surface_o3_station.template`
**现成样例**：`input_v3/constraints/bc_O3_demo.txt`（24 点日循环，实测可用）
**物理语义**：0-D 箱体以 `1/TAU_MIX` 速率向观测序列松弛（nudging 同化，
`chem_phys_mod` 物理账③）。

**原始数据 → 目标格式**：

```bash
# 假设你的原始数据: my_o3.csv, 两列 "小时序号, O3_ppb"
head -3 my_o3.csv
# 0,30.0
# 1,35.2
# 2,42.8

# 转换为模板格式 (t[s] + # unit: 头)
python3 - <<'PY'
import csv
rows = list(csv.reader(open("my_o3.csv")))[1:] if open("my_o3.csv").readline().startswith("h") else list(csv.reader(open("my_o3.csv")))
out = ["# unit: ppb",
       "# station: MYSITE",
       "# lat: 39.99  lon: 116.31",
       "# t[s]  O3"]
for r in rows:
    hh, v = int(r[0]), float(r[1])
    out.append(f"{hh*3600:8d} {v:8.2f}")
open("input_v3/constraints/o3_station_MYSITE.txt", "w").write("\n".join(out) + "\n")
print("写入", len(out)-4, "点")
PY

# 检查点数 (≤128)
grep -vc "^#" input_v3/constraints/o3_station_MYSITE.txt
```

**生产的目标文件**（照模板）：

```
# unit: ppb
# station: MYSITE
# lat: 39.99  lon: 116.31
# t[s]  O3
      0     30.0
   3600     35.2
   7200     42.8
  10800     55.1
```

**登记 + 编译 + 运行**：

```bash
# ① 登记到 bc.config（TAU_MIX 决定张弛强度）
printf '\n# 站点观测同化\nO3 = file:constraints/o3_station_MYSITE.txt\n' >> input_v3/bc.config

# ② 重编译数据（cfgbuild 校验单调/点数/单位/物种）
make runtime_v3.dat

# ③ 运行
cd facsimile_impl && ./fac_v3 SCEN_BAU ../output/assim_mysite ../runtime_v3.dat
```

**验证同化生效**（三种方式）：

```bash
# ① stdout 对比：O3max 应被拉向观测峰值
# ② IPR 的 B 账（BC+稀释）应出现明显非零
python3 - <<'PY'
import sys; sys.path.insert(0, "diagnostics")
import diagnostics as dg
sp, ipr = dg.ipr_from_fac("output/assim_mysite")
print("O3 的 B 账(同化+稀释) =", f"{ipr['O3']['B']:.3e}")
PY
# ③ 与无约束基线对比末 24h O3 均值
```

**TAU_MIX 调参**：`86400`（τ≈1 天=弱，模型主导）→ `21600`（6h=中）→
`7200`（2h=强，观测主导）。τ 太强会压制模型自身化学（OSSE 实验 F 有实证）。

**坑**：
- 时间轴必须是**自模拟起点秒**（本地时/UTC 差异自己换算）
- 缺测时段**直接省略行**（自动线性插值，端点钳制），不要填 0
- 超过 128 点必须先降采样（如 51 点/天 → 每小时 2 点）

---

### 14.2 模板②：AQI 多污染物宽表（需离线拆分换算）

**模板**：`templates/surface_aqi_station.template`
**现成样例**：无（AQI 需先拆分换算；拆分后格式与 `bc_O3_demo.txt` 相同）
**与模板①的区别**：AQI 站是**多污染物同报的宽表**（中国惯例 HJ 193，μg/m³），
需逐物种拆列 + 单位换算后才能接入。

**原始宽表格式**（模板内样例）：

```
# station: BJ-PKU (国控点)  method: HJ 193 连续自动监测
# 时间列: YYYYMMDDHH (北京时 BJST)   质量标记: HJ 193 附录A
yyyyMMddHH  SO2   NO    NO2   O3    CO    PM10  PM2.5  标记
2024060113    8     12    45    102   0.8   55    28     0
2024060114    7      9    38    118   0.7   51    25     0
2024060115    7     15    52    126   0.9   58    30     1
```

**单位换算表**（25 °C、1 atm：`ppb = μg/m³ × 24.45 / M`，M 为摩尔质量）：

| 物种 | M (g/mol) | 换算因子 | 例 |
|---|---|---|---|
| O₃ | 48.00 | **×0.509** | 102 μg/m³ = 51.9 ppb |
| NO₂ | 46.01 | **×0.531** | 45 μg/m³ = 23.9 ppb |
| SO₂ | 64.07 | **×0.382** | 8 μg/m³ = 3.06 ppb |
| CO | 28.01 | mg/m³ **×0.873** → ×1000 | 0.8 mg/m³ = 698 ppb |
| NO | 30.01 | ×0.815（**建议不接入**，见下） | — |
| PM₁₀/PM₂.₅ | — | **不接入**（气相机制无颗粒物种） | — |

**拆分换算脚本**（一次拆出 O₃/NO₂/SO₂/CO 四个文件）：

```bash
python3 - <<'PY'
import csv
FACTOR = {"O3": 0.509, "NO2": 0.531, "SO2": 0.382}
rows = [r for r in csv.reader(open("aqi_raw.txt"), delimiter="\t")
        if r and not r[0].startswith("#")]
hdr, data = rows[0], rows[1:]

for sp, f in list(FACTOR.items()) + [("CO", 0.873)]:
    col = hdr.index(sp)
    flag = hdr.index("标记")
    lines = [f"# unit: ppb", f"# station: BJ-PKU  species: {sp}  source: AQI ×{f}",
             f"# t[s]  {sp}"]
    t = 0
    for r in data:
        if int(r[flag]) != 0:          # 标记≠0 剔除 (HJ 193 质控)
            continue
        v = float(r[col])
        val = v * f * (1000 if sp == "CO" else 1)   # CO: mg/m³ → ppb
        lines.append(f"{t:8d} {val:8.2f}")
        t += 3600
    open(f"input_v3/constraints/aqi_{sp}_BJ-PKU.txt", "w").write("\n".join(lines)+"\n")
    print(f"{sp}: {len(lines)-3} 点")
PY
```

**机制物种映射（哪些 AQI 列能接入）**：

| AQI 列 | 机制物种 | 处理 |
|---|---|---|
| O₃ | O3 | ✓ 直接 |
| NO₂ | NO2 | ✓ 直接，但**建议弱 τ**（NO₂ 受光稳态支配，强 nudging 会劣化，OSSE 实证） |
| SO₂ | SO2 | ✓ 直接 |
| CO | CO | ✓ 直接 |
| NO | — | **不接入**（AQI 的 NO 是 NO−NO₂ 差折算值；机制内 NO 由化学诊断） |
| PM₂.₅/PM₁₀ | — | **不接入**（气相机制无颗粒物种） |

**登记**（可多物种，注意序列池 ≤8）：

```
TAU_MIX = 21600
O3  = file:constraints/aqi_O3_BJ-PKU.txt
SO2 = file:constraints/aqi_SO2_BJ-PKU.txt
CO  = file:constraints/aqi_CO_BJ-PKU.txt
```

**坑**：
- **质控标记必须剔除**（标记≠0 的行），否则异常值会污染插值
- 时间列 `YYYYMMDDHH` 是**北京时**，必须换算为"自模拟起点秒"
- NO 列是差值折算，别当实测值用

---

### 14.3 模板③：VOCs 组份监测（一站一物种）

**模板**：`templates/surface_voc_station.template`
**现成样例**：`input_v3/constraints/bc_HCHO_demo.txt`（HCHO 24 点，实测可用）
**格式**：与模板①**完全相同**（两列 + `# unit: ppb`），只是**一个物种一个文件**。

```
# unit: ppb
# station: BJ-PKU  species: HCHO  method: DNPH-HPLC
# t[s]  HCHO
      0    2.40
   3600    2.85
   7200    3.31
  10800    3.92
  14400    4.55
  18000    5.10
```

**可约束物种**（mechanism_3 内，模板列出的）：

| 分类 | 物种 |
|---|---|
| 直接接入 | HCHO、ISOP、MVK、MACR、CH₃CHO |
| 可接（需谨慎） | CH₄、CO、CH₃O₂、CH₃OOH、PAN、MGLY |
| **不接入** | 苯系物 BTEX（机制无显式芳烃化学，仅 ARO1/ARO2 集总）、C₃-C₈ 烷/烯细分类（机制仅 ALK/OLE 集总）、乙炔/丙烷（可作 CH₄ 通道的弱替代约束） |

**注意**：若观测是 **PAMS-57/TO-15 全谱**，机制里的 `ALK`/`OLE`/`ARO1`/`ARO2`
是**集总物种**——需要你按源谱加权求和成单一等效值，或改为约束其**降解产物**
（如用 HCHO 观测间接约束 VOC 活性）。

**登记**（每物种一行）：

```
TAU_MIX = 21600
HCHO = file:constraints/voc_HCHO_BJ-PKU.txt
ISOP = file:constraints/voc_ISOP_BJ-PKU.txt
```

**物理适用性**：
- **短寿命物种**（HCHO τ~1h、ISOP τ~分钟）nudging 有意义（观测逼近真实）
- **长寿命物种**（CH₄）建议仍用 `emissions.config` 排放驱动，而非观测约束

**坑**：CH₃O₂ 等自由基物种观测极少，若强行约束会破坏自由基平衡（谨慎）。

---

### 14.4 模板④：卫星柱浓度（唯一"引擎不直读"，用 sat2surf.py 转换）

**模板**：`templates/satellite_column.template`
**现成样例**：`input_v3/constraints/bc_NO2_sat_demo.txt`（逐日 6 点，实测可用，需弱 τ）
**关键限制**：ozone_lab 是 **0-D 箱体，无气柱概念**——引擎**不直接读**卫星柱浓度。
必须**离线转换**为地表等效体积混合比后，走模板①同款 BC 通路。

> ✅ **别再手算了**：转换已工具化 —— **`src/sat2surf.py`**（内置两法交叉核对 +
> QA 过滤 + 点数限制 + 量级自检）。手算极易踩下面那个 ×1e9 的坑。

#### 14.4.1 转换工具用法

```bash
# 最小用法（BLH 取常数 1000 m）
python3 src/sat2surf.py 原始柱.raw 输出约束.txt

# 完整选项
python3 src/sat2surf.py <raw.raw> <out.txt> \
    --blh 1000 \                    # 混合层高度常数 [m]
    --blh-file input_v3/constraints/met_BLH.txt \   # 或 BLH 序列, 过境时刻线性插值
    --qa-min 0.75 \                 # QA 阈值 (也可由文件头 qa_filter 解析)
    --species NO2 \                 # 覆盖头部 species
    --t0 2024-06-01 \               # 模拟起点日期 (缺省=首条观测当日 00:00)
    --max-points 128                # 超限自动等间隔抽稀

# 自检（3 组已知值 + 两法互证 + 旧式反证）
python3 src/sat2surf.py --selftest
```

#### 14.4.2 输入：原始柱浓度（照模板）

```
# satellite: TROPOMI/S5P  product: OFFL  L2  v: 2.4.0
# species: NO2_tropo
# overpass: 13:30 local
# qa_filter: qa_value > 0.75
# t[YYYY-MM-DDTHH:MM]  column[molec/cm2]  qa
2024-06-01T13:30   4.2e15   0.92
2024-06-02T13:30   3.8e16   0.88
2024-06-03T13:30   5.1e16   0.95
2024-06-05T13:30   3.2e16   0.91
2024-06-06T13:30   2.8e16   0.89
```
支持时间格式：`YYYY-MM-DDTHH:MM` / `YYYY-MM-DD HH:MM` / `YYYYMMDDHH`。
QA 低于阈值（默认 0.75）的过境自动剔除并打印。

> 📋 **全物种输入矩阵**：74 物种 × 5 通道（IC/BC/EMIS/DRY/WET）的完整可输入清单
> 与物理合理性判断见 [templates/SPECIES_INPUT_MATRIX.md](templates/SPECIES_INPUT_MATRIX.md)。
> 每物种每个通道标注 **✓ 推荐 / △ 可行有前提 / ✗ 无意义**。

#### 14.4.2b 支持哪些卫星产品？（无白名单，凡机制物种皆可）

`sat2surf.py` **没有物种硬编码**——`species` 从文件头解析（或 `--species` 覆盖），
只要能对上机制 `VARIABLE` 段即可。S5P TROPOMI 常用产品到 mechanism_3 的映射：

| TROPOMI 产品 | 机制物种 | 可接入 | 实测/注意 |
|---|---|---|---|
| **HCHO_tropo** | HCHO | **✓ 已验证** | 实测 8.5e15–1.8e16 柱 → **3.46–7.32 ppb**（实测值，合理） |
| **NO2_tropo** | NO2 | **✓ 已验证** | 实测 4.2e15–5.1e16 柱 → 1.71–20.73 ppb |
| SO2_tropo | SO2 | ✓ | 机制有硫链；柱量级小（1e15–1e16） |
| CO | CO | ✓ | TROPOMI CO 为**总柱**，BLH 假设误差更大 |
| CH4 | CH4 | ✓ | 总柱为干空气混合比，换算方式不同（不能直接除 BLH），需先转体积混合比 |
| O3 | O3 | ✓ | TROPOMI O3 是廓线/L2 产品，非常规柱，需专门处理 |
| 气溶胶/云 | — | ✗ | 气相机制无对应物种 |

**HCHO 专项注意**（与 NO₂ 的差异）：

| 项 | HCHO | NO₂ |
|---|---|---|
| 柱量级 | 2e15（清洁）– 3e16（生物源/火） | 1e15（清洁）– 6e16（城市） |
| 转换后地面值 | **1–10 ppb** | **5–20 ppb** |
| 机制内角色 | 有排放（5.8e4）+ 光解 + OH 氧化，τ~1–2 h | 光稳态支配，τ~分钟级 |
| 建议 τ_mix | **中等（21600 s）**——它有真实局地化学，强约束会压制 | **弱（43200+）**——强 nudge 会破坏光稳态（OSSE 实证） |
| 垂直分布 | 主要在对流层低层 → 除 BLH 合理 | 混合较均匀，除 BLH 亦可 |

**HCHO + NO₂ 同时约束 = 卫星 FNR 场景**（推荐组合）：

```bash
# 两个文件分别转换
python3 src/sat2surf.py sat_HCHO_202406.raw input_v3/constraints/sat_HCHO_surf.txt \
        --blh-file input_v3/constraints/met_BLH.txt
python3 src/sat2surf.py sat_NO2_202406.raw  input_v3/constraints/sat_NO2_surf.txt \
        --blh-file input_v3/constraints/met_BLH.txt

# bc.config 同时登记
HCHO = file:constraints/sat_HCHO_surf.txt
NO2  = file:constraints/sat_NO2_surf.txt
```
**实测结果**（本次，6 天 5 次有效过境）：
```
cfgbuild: 序列 ['BLH','RAIN','bc_HCHO','bc_NO2'] ✓
SERIES:   bc_HCHO 8.4993e+10 … / bc_NO2 4.19922e+10 …   (= ppb×2.46e10 ✓)
fac_v3:   O3max = 101.78 ppb, nst = 5750
IPR 分账: HCHO_B = −2.517e+11 ✓ | NO2_B = +5.490e+11 ✓  (双物种同化均生效)
卫星柱比 FNR 均值 0.71  vs  地面 FNR 1.23
  → 差异来自 BLH 假设 + HCHO 局地化学生成（P=1.86e12），合理
```

#### 14.4.3 输出的转换报告（实测）

```
$ python3 src/sat2surf.py sat_NO2_demo.raw sat_NO2_surf.txt
[QA] 阈值 0.75: 剔除 1 条低质量过境 (保留 5)

== 转换结果 → sat_NO2_surf.txt ==
  过境时刻                   柱[molec/cm2]   BLH[m]    表面[ppb]
  2024-06-01T13:30        4.200e+15     1000       1.71
  2024-06-02T13:30        3.800e+16     1000      15.45
  2024-06-03T13:30        5.100e+16     1000      20.73
  2024-06-05T13:30        3.200e+16     1000      13.01
  2024-06-06T13:30        2.800e+16     1000      11.38

  范围 1.71–20.73 ppb | 点数 5 | 两法最大相对差 2.7e-16
  下一步: bc.config 加行  NO2 = file:sat_NO2_surf.txt
```

产出文件（可直接登记）：
```
# unit: ppb
# source: sat2surf.py ... | species=NO2 | QA>=0.75 | BLH=1000m
# method: 柱/BLH/2.46e10 (直观式) 与 混合比法 交叉核对 max_rel=2.7e-16
# t[s]  NO2
   48600     1.707
  135000    15.447
  221400    20.732
  394200    13.008
  480600    11.382
```

#### 14.4.4 转换公式（工具内部实现，两法强制互证）

```
② 直观法（主用，免"柱高"假设）:
   表面等效浓度 [molec/cm³] = 柱浓度 [molec/cm²] / BLH[cm]
   表面等效混合比 [ppb]     = 上式 / 2.46e10

① 混合比法（核对）:
   混合比 = 柱/(2.15e25) × 1e9 / (BLH/8.74e5)
```

**算例**：NO₂ 柱 3.8e16 molec/cm²、BLH = 1 km →
`3.8e16/1e5/2.46e10 = 15.45 ppb`（两法一致到 1e-16，合理量级）。

#### 主准则与参考法的分工（重要）

**直观法 = 主准则**，理由：它只依赖 **2.46e10**（1 ppb 数密度，由 P/kT 直接得出）
与**实测 BLH**，物理上就是"柱浓度摊到混合层厚度 = 层平均浓度"的**定义式**，假设最少。
混合比法用 `2.15e25`（空气柱）与 `8.74 km`（等效柱高）——**这两者本是同一
"标准大气柱"假设的两种表述**，比直观法多一层假设，故仅作参考核对。

**两法分歧时的行为（v2.4 起，按用户要求修改）**：

| 情形 | 旧行为（v2.3） | **新行为（v2.4）** |
|---|---|---|
| 两法相对差 > 1e-9 | `[错误]` + **exit 2 中断** | **`[警告]` 到 stderr，不中断**；一律**以直观法结果写出文件** |
| 退出码 | 2 | **0**（有产物即成功） |

警告示例：
```
[警告] 两法分歧: 最大相对差 3.37e+00 > 1e-9 — 以直观法 (柱/BLH/2.46e10,
假设最少) 为准继续; 混合比法仅供参考, 不中断
```
产物注释同时留痕：`... max_rel=3.4e+00  [警告: 两法分歧, 已按直观法输出]`。

**如何主动制造分歧来验证该路径**：`--scale-height-km 2.0`
（让参考法假设柱混合到 2 km 而非 8.74 km）→ 两法差 2.19 倍，工具仍 exit 0 并
按直观法输出 15.45 ppb（已固化为功能测试断言）。

**默认情况下两法代数等价**（`AIR_SCALE` 由 `AIR_COL/(PPB×1e9)` 派生），
相对差 ~1e-16——所以"分歧"只在显式给出不同柱假设时出现。

> #### ⚠ 已修复的量纲 bug（历史，工具已内置防线）
>
> **原模板公式**：`表面等效浓度 [ppb] ≈ 柱/(空气柱 × 混合层分数)`
>
> **三条独立证据证明它错**：
>
> | # | 证据 | 数值 |
> |---|---|---|
> | 1 | 量纲分析：结果是摩尔混合比（mol/mol），ppb 需 ×1e9 | `3.8e16/(2.15e25×0.1144) = 1.54e-8`，×1e9 才是 **15.45 ppb** |
> | 2 | 与直观式独立对照 | 两法一致到 **0.0%**；原式差 **9 个量级** |
> | 3 | 模板自相矛盾 | 模板"转换后示例"写 4.85 ppb，但原式算 4.2e15 柱只得 1.7e-9 |
>
> **工具防线**（`sat2surf.py`）：
> - 两法交叉核对：分歧**仅警告不中断**（exit 0），一律以直观法（主准则）输出
> - `AIR_SCALE` 由 `AIR_COL/(PPB×1e9)` **解析派生**（手写 8.74e5 会引入 1.9e-5 偏差破坏互证）
> - 输出量级自检（>500 或 <0.001 ppb 报警）
> - 已固化为功能测试 **FT-SAT**（3 组断言，见 §17）

#### 14.4.5 登记与运行（全链路）

```bash
# ① 转换
python3 src/sat2surf.py sat_NO2_demo.raw input_v3/constraints/sat_NO2_surf.txt

# ② 登记
printf '\nNO2 = file:constraints/sat_NO2_surf.txt\n' >> input_v3/bc.config

# ③ 编译（cfgbuild 会打印 "序列: [..., 'bc_NO2']"）
make runtime_v3.dat

# ④ 运行
cd facsimile_impl && ./fac_v3 SCEN_BAU ../output/sat_run ../runtime_v3.dat
```

**全链路实测结果**（本次验证）：
```
cfgbuild: 序列 ['BLH','RAIN','bc_NO2'] ✓
SERIES 值: 4.19922e+10 3.79996e+11 ...  (= ppb × 2.46e10 ✓)
fac_v3:   O3max = 102.31 ppb, nst = 5566
IPR 分账: NO2_B = +5.523e+11   (无约束基线 -4.975e+10 → 同化增量 +5.0e11 ✓)
```

#### 14.4.6 坑与注意

| 项 | 说明 |
|---|---|
| **每日一次过境** | 序列采样 86400 s，cfgbuild 线性插值到逐时窗；τ 建议取较长（如 43200） |
| **BLH 敏感性极强** | 同样柱浓度，BLH 400 m→38.6 ppb，BLH 1375 m→1.24 ppb（实测差 31 倍）——务必用真实 BLH 序列 |
| **qa 掩膜必须做** | 工具默认 0.75 阈值；云污染数据不剔除会严重偏置 |
| 柱≠表面浓度 | 转换必然引入 BLH 假设，这是 0-D 模型的固有边界（不是 bug） |
| FNR 柱比 | 卫星柱比 HCHO/NO₂ 与地面 FNR 的换算关系见 ridge5 实验 |

### 14.5 模板⑤：排放清单（直接对应 emissions.config）

**模板**：`templates/emissions_inventory.template`
**这是唯一"直接就是配置"的模板**——正文语法即 `emissions.config` 语法。

#### 格式与两种单位（V2.3 起支持 ppb/h 指令）

```
<物种名>  <基线速率>  <日型: BOX|FLAT>  <振幅 amp>
```

**单位二选一**（在文件头用指令行声明）：

| 指令 | 含义 | 换算 |
|---|---|---|
| （缺省） | molec/cm³/s | 无 |
| `# unit: ppb/h` | ppb 每小时 | cfgbuild 自动 **×2.46e10/3600 = ×6.833e6** |

**实测验证**（本次）：

```bash
cat > /tmp/t/emissions.config <<'EOF'
# unit: ppb/h
NO      20.0    BOX  0.9
NO2     1.0     BOX  0.9
EOF
python3 src/cfgbuild.py /tmp/t /tmp/rt.dat --mech mechanism/mechanism_3_full_chain.fac
# EMIS 段输出: E 1.366667e+08 BOX 0.9   (20×6.833e6 ✓)
#             E 6.833333e+06 BOX 0.9   (1×6.833e6 ✓)
```

#### 三种来源的换算（若你不用 ppb/h 指令）

**(a) ppb/h → molec/cm³/s**（可直接用 `# unit: ppb/h`）
```
E = C_ppb_per_hour × 2.46e10 / 3600 = C × 6.833e6
```

**(b) t/yr 质量清单 → molec/cm³/s**
```
E = Q[t/yr] × 1e6 / M × 6.022e23 / (3.154e7 × A[cm²] × H[cm])
```
**例**：NO₂ 1000 t/yr、100 km²、BLH 1 km → 4.15e9 molec/cm³/s

**(c) mol/m²/s 面通量**
```
E = F × 6.022e23 / (1e4 × H[cm])
```

**合理性自检**：稳态浓度 ≈ E×τ。NO 的 τ≈79 s（O₃ 滴定）；
若 E=2.157e7 → SS≈1.7e9 ≈ 70 ppb（城市量级 ✓）。

#### 日型与振幅

| 日型 | 公式 | 适用 | 注意 |
|---|---|---|---|
| `BOX` | 06–18 时 `E0×(0.1+0.9·amp·sin²(π(h−6)/12))`；夜间 `0.1×E0` | 交通、工业、光照驱动生物源 | **amp=0 仍排 10%**；amp∈[0,1] |
| `FLAT` | 全天 `E0×1.0` | CH₄、DMS 等恒定源 | 振幅无效 |

**日均与峰值的换算**（设置基线前必须想清）：
```
24h 日均 = E0 × (0.1 + 0.225×amp)    (sin² 只在 6-18 时生效, 权重 12×0.5/24=0.25)
```
- 若清单给的是**日均**，且 amp=1.0 → `E0 = 日均 / 0.55`
- amp=0.9（当前配置）→ 日均 = 0.505×E0

#### 完整操作

```bash
# ① 编辑清单（模板正文可直接用作起点）
cp templates/emissions_inventory.template input_v3/emissions.config
# 或追加自己的行
cat >> input_v3/emissions.config <<'EOF'
# --- 我的清单 (ppb/h 直给, V2.3 特性) ---
# unit: ppb/h
NO      20.0    BOX   0.9     # 早高峰 20 ppb/h
CO      5.0     BOX   0.85
ISOP    1.5     BOX   0.9
EOF
# ② 重编译（cfgbuild 校验物种/单位/日型/振幅）
make runtime_v3.dat
# ③ 运行
cd facsimile_impl && ./fac_v3 SCEN_BAU ../output/myemis ../runtime_v3.dat
```

**选择性输入（重要）**：本文件是**清单**不是**表格**——你**只列需要模拟的物种**即可。

| 你写几行 | 程序行为 |
|---|---|
| **全空文件** | ✅ 合法：无人为排放，纯背景化学（实测可跑，O₃max=0） |
| **只写 2 行**（如 NO + ISOP） | ✅ 合法：其余 72 物种零排放（实测可跑，O₃max=11.6 ppb） |
| **写 15 行**（当前示例） | ✅ 示例配置的完整度，**不是强制要求** |

**未列出的物种 = 零排放**（不是"报错"、不是"用默认值"）。

当前示例配置覆盖 15 个物种（NO/NO₂/CO/HCHO/CH₃CHO/CH₄/ISOP/SO₂/DMS/ALK/OLE/ARO1/ARO2/KET/TERP）。

**情景组映射**（影响 `scenarios.config` 的 E_* 因子）：
```
E_NOx → NO, NO2
E_VOC → HCHO CH3CHO ISOP DMS MVK MACR MSA SO2 ALK OLE ARO1 ARO2 KET
E_CO  → CO
E_CH4 → CH4
```
新排放物种要参与缩放 → 在 `modules/runtime_dat_mod.f90` 的 `build_groups()` 加一行。

**坑**：
- `box` 小写 → 静默按 FLAT 处理（cfgbuild 现已拒绝，但注意错误提示）
- 负振幅 → 白天变负排放（已拒绝）
- 忘记重编译 → 跑的还是旧排放

---

### 14.6 模板⑥：气象序列（`met_series.template`）

**模板**：`templates/met_series.template`（6 条现成样例：`input_v3/constraints/met_*_demo.txt`）

`meteo.config` 的 7 个量都支持 `file:`（两列 + 可选注释）：

```
TEMP  = file:constraints/met_TEMP_mysite.txt      # K
PRES  = file:constraints/met_PRES_mysite.txt      # hPa
RH    = file:constraints/met_RH_mysite.txt        # %
BLH   = file:constraints/met_BLH_mysite.txt       # m，必须 > 0（编译期强制）
RAIN  = file:constraints/met_RAIN_mysite.txt      # mm/h
JFAC  = file:constraints/met_JFAC_mysite.txt      # 0–1（云订正）
RADIAT = file:constraints/met_RADIAT_mysite.txt   # W/m²（保留，当前不参与 J）
```

**⚠ 与模板①的关键差别：气象序列里不能写 `# unit:` 行。**
那条指令是浓度换算因子（ppb→2.46e10），写进气象文件会把 298.15 K
静默变成 7.33e12。`cfgbuild` 现对气象序列调用 `read_series(allow_unit=False)`，
遇到直接报错（`气象序列不接受单位指令行`），并有回归用例守卫。

**注意**：
- BLH 用于物理账②（`vd/BLH`）与账④（`−y·dB/BLH`）——**序列最小值必须 > 0**
- 序列池与 BC 共享 ≤8 条。**6 变量全序列时 BC 只剩 2 条**
- `JFAC` 序列 = 逐时云订正（比 `SCEN_CLOUD` 的全程常数更真实）
- 序列不循环：超出首/末时刻钳制端点值，想覆盖 5 天就写满 120 h

**完整示例（实验 M，10 个 run 实测）**：

```bash
python experiments/run/run_M_met_series.py
# M0_const  六变量全取均值   O3max=111.15
# M1_all    六变量全序列     O3max= 78.05
# 单变量归因: BLH −45.4 | TEMP +16.8 | JFAC +2.6 | RH −1.2 | RAIN +0.2 | PRES 0.00
```

单变量效应**不可加和**（和 −27.0 vs 实际 −33.1），`attribution.csv` 末列的
−137%/+50.6% 是这一非线性的定量标志，不是错误。

---

### 14.7 新增政策情景

```bash
cat >> input_v3/scenarios.config <<'EOF'
SCEN_MINE = E_NOx*0.7 E_CO*0.5 E_CH4*1 E_VOC*1 JFAC*1
SCEN_HOT2 = E_NOx*1 E_VOC*1 E_CO*1 E_CH4*1 JFAC*0.8 DT+3
EOF
make runtime_v3.dat
cd facsimile_impl && ./fac_v3 SCEN_MINE ../output/mine ../runtime_v3.dat
```

**语法**：只支持 `E_*乘法因子` 与 `DT±n`（加法，温度偏移）；
不支持分时段因子。cfgbuild 静态解析，解析失败即编译期拒绝。

---

### 14.8 扩机制（最深的改动）

```bash
# ① .fac 四处改动
vim mechanism/mechanism_3_full_chain.fac
#   (a) VARIABLE 段加物种名（分号结尾）
#   (b) 速率名定义（如需要）: KXXX = A*EXP(-Ea/TEMP)
#   (c) 反应行: % <速率> : <反应物> = <产物> ;
#   (d) FORMULA 段加元素组成: 新物种 C1 H2 O1

# ② 重编译
make v3

# ③ 四道质检关卡（§17）
python audit/element_audit.py mechanism/mechanism_3_full_chain.fac
make test

# ④ 三张表同步（关键教训——OLEOOH 零汇 bug 的根因）
#   (a) emissions.config  若新物种有排放
#   (b) drydep.config     若新物种无化学汇（否则无界累积）
#   (c) build_groups()    若需参与情景缩放
#   (d) NOy 权重表（若含氮）+ ATOL 调优表（若高浓度）

# ⑤ 重编译数据
make runtime_v3.dat
```

**自查脚本**（加物种后必跑，防 OLEOOH 类问题）：

```bash
python3 - <<'PY'
# 化学死端检查: 只产不耗的物种是否有物理汇
ns = 74  # ← 改成你的新物种数
sp = []
for ln in open("runtime_v3.dat"):
    if ln.startswith("SPECIES"):
        sp = ln.split()[1:1+ns]; break
prod, loss = [0]*ns, [0]*ns
for ln in open("facsimile_impl/gen_v3/oz.reac"):
    p = ln.split(); n = int(p[0])
    for k in p[1:1+n]: loss[int(k)-1] += 1
for ln in open("facsimile_impl/gen_v3/oz.prod"):
    p = ln.split(); n = int(p[0]); tk = p[1:]
    for i in range(n): prod[int(tk[2*i])-1] += 1
dry = {l.split()[0] for l in open("input_v3/drydep.config")
       if l.strip() and not l.startswith("#")}
dead = [sp[i] for i in range(ns) if prod[i] > 0 and loss[i] == 0 and sp[i] not in dry]
print("零汇物种 (必须修!):", dead if dead else "无")
PY
```

---

### 14.9 多物种联合约束

```
TAU_MIX = 21600

# 气象占 2 条 (BLH/RAIN) → BC 最多 6 条
O3   = file:constraints/o3_station_A.txt
NO2  = file:constraints/no2_station_A.txt
CO   = file:constraints/co_station_A.txt
HCHO = file:constraints/voc_HCHO_A.txt
ISOP = file:constraints/voc_ISOP_A.txt
# 剩 1 条余量
```

**超池的处置**：合并为更长序列不行（不同物种）；应减少约束物种或
把气象改回常数（释放 2 条池位）。

---

### 14.10 接入后的验证方法（四查）

| 查什么 | 怎么查 | 期望 |
|---|---|---|
| 1. 编译期校验过 | `make runtime_v3.dat` 无报错，末尾打印序列/情景列表 | `序列: ['BLH','RAIN','O3']` |
| 2. 数据真的进去了 | `grep -A3 "^SERIES" runtime_v3.dat` | 能看到你的序列名与点数 |
| 3. 物理账生效 | 读 `_ipr.csv` 的 `B` 列（BC+稀释） | 同化时明显非零 |
| 4. 结果合理 | stdout 的 `O3max` + `route_indicators` 的 `O3_ppb` | 量级符合预期（30–150 ppb） |

```bash
# 一条命令做 2+3
grep -A2 "^SERIES" runtime_v3.dat | head -6
python3 -c "
import sys; sys.path.insert(0,'diagnostics'); import diagnostics as dg
sp, ipr = dg.ipr_from_fac('output/assim_mysite')
print('O3 B账:', f\"{ipr['O3']['B']:.3e}\")"
```

---

### 14.11 常见接入错误对照表

| 症状 | 根因 | 修复 |
|---|---|---|
| `序列 129 点 > 上限 128` | 逐时数据太长 | 降采样（2 点/小时或 1 点/2 小时） |
| `序列总数 9 > 上限 8` | 约束太多 | 减物种或气象改常数 |
| `时间非单调 @ t=…` | 时间轴乱序/时区混用 | 统一换算为自模拟起点秒并排序 |
| `未知单位 'ug/m3'` | 用了未支持的单位 | AQI 必须**离线**换成 ppb（模板②） |
| 结果与无约束完全一致 | 忘记 `make runtime_v3.dat` | 重编译 |
| 结果与观测差很远 | TAU_MIX 太大（弱约束） | 改 21600 或 7200 |
| NO₂ 同化后反而变差 | NO₂ 受光稳态支配，强 nudging 有害 | 用弱 τ 或改 NOx 排放反演 |
| 卫星数据量级离谱 | 混合层分数假设错 | 按 §14.4 公式重算 |
| 排放改了但结果没变 | emission 未重编译 / 物种名拼错 | `make runtime_v3.dat` + 核对 VARIABLE |

---

## 15. 输出契约（含真实样例行）

### 15.1 fac_v3 `_conc.csv`

```
h,NO,NO2,O3,O1D,O,OH,HO2,H2O2,CO,HCHO,HNO3,HONO,NO3,N2O5,CH4,CH3O2,...
   96,  5.83904E+07,  4.92624E+10,  9.75556E+11,  0.00000E+00, -3.10947E-43, ...
```
- 74 列；`h` 从 96 起（末 24 h）；单位 molec/cm³
- 负的 1e-43 量级 = 数值零；表头顺序 = 机制 VARIABLE 顺序（非字母序）

### 15.2 fac_v3 `_irr.csv`
```
h,r1,r2,...,r131
```
每行 = 该小时的各反应速率（molec/cm³/s）。

### 15.3 fac_v3 `_ipr.csv`
```
h,NO_P,NO_L,NO_E,NO_D,NO_B,NO_X,NO2_P,NO2_L,...
```
每物种 6 列；**值是 molec/cm³**（不是速率——小时窗精确积分）。
六账：P 机制产 / L 机制耗 / E 排放 / D 干湿沉 / B BC+稀释 / X BLH 膨胀。

### 15.4 fac_batch_omp `_indicators.csv`（24 列）
见 §6.5 完整列表。

### 15.5 mcm_run `_conc.csv`
```
h,O3,NO,NO2,OH,HO2,HCHO,CO,H2O2,HNO3,NO3,N2O5,HONO,P_O3,L_O3
```
15 列（含 P_O3/L_O3 即时速率）；文件名 = CLI 的 out_pref（v2.3.1 修复后）。

### 15.6 反应号 ↔ 反应式映射

```bash
# 第 k 个反应的反应物索引
sed -n "${k}p" facsimile_impl/gen_v3/oz.reac
# 第 k 个反应的产物与系数
sed -n "${k}p" facsimile_impl/gen_v3/oz.prod
# 或直接看机制文件第 k 个 % 行
grep "^%" mechanism/mechanism_3_full_chain.fac | sed -n "${k}p"
```

---

## 16. 指标解读手册（推导 + 实测 + 冲突案例）

### 16.1 FNR = [HCHO]/[NO₂]（Sillman 1995 / Tonnesen-Dennis 2000）

**推导**：HCHO ≈ ∫[VOC×k_OH]dt（VOC 氧化累积量）；NO₂ ≈ 瞬时 NOx 催化剂存量。
比值反映"燃料消耗 vs 催化剂供给"的平衡。

| FNR | 体制 | 机理解释 | 策略 |
|---|---|---|---|
| **< 1** | VOC-limited（NOx 饱和） | NO₂ 与 OH 反应（→HNO₃）速率高，OH 被压制，自由基池受限 | 减 NOx 释放 OH |
| 1–2 | 过渡带 | 双敏感；阈值随纬度/季节漂移（高纬 ~2.7, Jin 2017） | 需模拟判别 |
| **> 2** | NOx-limited | HO₂/RO₂ 丰富，NOx 是瓶颈 | 减 NOx 线性降 O₃ |

**记忆口诀**：分子是"烧掉的柴"，分母是"手里的火钳"——火钳（NO₂）太多会把
消防员（OH）缠住，所以 FNR 小反而该减 NOx。

### 16.2 屋脊（O₃–NOx 曲线转折）

**定义**：固定 VOC，O₃ 随 NOx 先升后跌的峰值点。
**物理**：左侧 = 加 NOx 促进 HO₂+NO/RO₂+NO（NOx-limited）；右侧 = NOx 过量导致
OH+NO₂→HNO₃ 主导（NOx 饱和）。**V2.3 实测**：fnox=2.28（O₃max 125.5 ppb）。

### 16.3 L_NOx/Q 与 FNR 的冲突案例（重要！）

**实测冲突**：BAU 配置下 FNR=3.427（NOx-limited）而 L_NOx/Q=2.011（>0.5，
Kleinman 判据指向 NOx 饱和）。

**原因分析**：
- `L_NOx` 在 `indicators.py` 定义为 **NO₂+OH→HNO₃ 与 NO+OH→HONO 两条反应之和**
- HONO 是**可逆储库**（HONO+OH→NO₂ 又放回 NOx），把它计入"NOx 损失"会高估 L_NOx
- 本配置 NOy 仅 10.7 ppb（稀缺），NOx 被 OH 快速转化为储库很正常，
  但这不等于"NOx 饱和"

**教训**：
1. **任何单一指示因子都不足以下结论**——必须 FNR + 屋脊几何 + 情景实验三方交叉
2. 指标定义（尤其"损失"与"终止"的边界）决定其数值，**跨库比较前先看定义**
3. 本项目最终以**情景实验（A）+ 屋脊（B）**为准，指标是辅助

### 16.4 OPE = P(O₃)/P(NOz)

每消耗一个 NOx 分子生成多少 O₃。NOx 饱和区触底（V2.3 谷 6.881）——
NOx 过剩时被"浪费"在 HNO₃ 上。

### 16.5 NOy 家族（17 物种，氮守恒的基础）

```
NO + NO₂ + HNO₃ + HONO + NO₃ + N₂O₅(×2) + PAN + MPAN + PNA
  + RONO₂ + OLENIT + ISOPNIT + ISOPDNT + TERPNIT
```
运行期不变量 `Σ w_N·ydot_chem ≈ 0` 覆盖全部 17 物种——V2.3 漏计
PNA/RONO₂/OLENIT/TERPNIT 时，该不变量抓住了 4 个格点的 FAIL。

---

## 17. 质检四关卡（机制/代码变更后必过）

| # | 关卡 | 工具 | 抓什么 | 通过标准 |
|---|---|---|---|---|
| 1 | 元素守恒 | `audit/element_audit.py` | 每反应 C/H/N/O/S/Cl/Br 配平 | 131/131 PASS |
| 2 | 氮守恒不变量 | batch `diag_selfcheck`（运行期） | NOy 权重表漏物种 | FAIL = 0 |
| 3 | 金值哨兵 | `functional_test` SENT 表 | 机制/求解改动的数值漂移 | 4 点位 ±10–15% |
| 4 | mcm_run 对照 | `mcm_city/mcm_run` | 无机块跨机制行为 | 定性一致 |

**金值哨兵 4 点位**（V2.2 标定，V2.3 实测全在容差内）：

| 点位 | 金值 | 列索引 | V2.3 实测 |
|---|---|---|---|
| 屋脊峰 O₃max | 132.10 | 4 | 125.5 ✓ |
| 饱和区 O₃ | 93.68 | 4 | 85.75 ✓ |
| FNR 转折 | 3.539 | 8 | 3.447 ✓ |
| OPE 谷 | 6.487 | 16 | 5.933 ✓ |

---

## 18. 故障排查速查表

### 18.1 构建与运行

| 症状 | 根因 | 处置 |
|---|---|---|
| `Cannot open file 'gen_v3/oz.reac'` | cwd 不对 | `cd facsimile_impl` 再跑 |
| `make` 链接失败 | CVODE 路径/库名 | 检查 `~/soft/atmenv/cvode/lib` |
| `make v3` 说 `Nothing to be done` 但结果没变 | 数据未重编译 | `make runtime_v3.dat` |
| 跑得极慢（>10 s） | 机器负载高 | 看 `uptime`；负载 >核数时性能无意义 |

### 18.2 配置与数据

| 症状 | 根因 | 处置 |
|---|---|---|
| `排放未知物种 'XXX'` | 物种 ∉ 机制 VARIABLE | 核对拼写/大小写 |
| `序列 129 点 > 上限 128` | 数据太密 | 降采样 |
| `序列总数 9 > 上限 8` | 气象+BC 超池 | 合并或减少 |
| `时间非单调 @ t=…` | 时间轴乱序/重复 | 排序去重 |
| `未知单位 'x'` | 单位指令拼错 | `cm3\|ppb\|ppm` |
| `BLH 序列含非正值` | 数据有 0/负 | 修正气象数据 |
| `BOX 振幅 ∉ [0,1]` | 振幅非法 | 改到 [0,1] |
| `日型 'box' 非法` | 大小写错 | 改 `BOX` |

### 18.3 结果异常

| 症状 | 根因 | 处置 |
|---|---|---|
| O₃ 单调爆炸到几百 ppb | 排放单位错（差 1000×） | 按 §3.3 重算 |
| O₃ 恒为 0 | 排放全零 / 光解全零 | 查 emissions.config / JFAC |
| 某物种无界累积 | 无化学汇且未配沉降 | drydep.config 加该物种（OLEOOH 教训） |
| 第 1 天与后 4 天差异大 | 滴定建立期（正常） | 只信末 24h |
| CSV 值 `8.27280-108` | gfortran 3 位指数省略 E | 诊断库已健壮解析；自己读注意 |
| 黎明步长崩溃 | 无预条件 SPGMR（历史） | fac_v3 用 DENSE，勿切 SPGMR |
| 哨兵超差 | 机制/数值被改 | 查四关卡，必要时有意识重标定 |

### 18.4 实验脚本

| 症状 | 处置 |
|---|---|
| `AttributeError: module 'diagnostics' has no attribute 'report'` | 用 `dg.load`/`dg.ipr_from_fac` + `indicators.*`（见 §8.3） |
| ridge5 报 `fac_run 已退役` | 预期行为；现役路线会自动跑 |
| scen5 的 mcm 方向 ✗ | 真实的跨机制分歧，非 bug |

---

## 19. 性能与耗时预算（4 核 MacBook 实测）

| 操作 | 耗时 | 备注 |
|---|---|---|
| cfgbuild 数据重编译 | **0.37 s** | 纯 Python |
| fac_v3 单情景 5 天 | **0.19–0.87 s** | 负载相关 |
| fac_batch_omp 81 案例（4 线程） | **~47 s** | 零失败 |
| mcm_run 5 天（29sp） | 0.09 s CPU | — |
| 机制变更后重编译 | ~8 s | mechc + gfortran |
| `make test` 全链 | <2 min | — |
| 实验 B→D→G 一轮完整诊断 | ~5 min | — |
| `sens_adj.py` | ~2 min | OAT+SVD |
| `obs_osse.py` | ~2–5 min | 多次重跑调库 |

**负载警告**：OpenMP 加速比在系统 load > 核数时测量无意义；
`scan_parallel.sh` 内置负载闸门（SKIP 而非 FAIL）。

---

## 20. Python API 速查

### 20.1 `diagnostics.py`（3 个函数）

```python
import sys; sys.path.insert(0, "diagnostics")
import diagnostics as dg

h, C, irr_h, I = dg.load("output/expA_BAU")
# h: (T,) 小时   C: (T, ns) 浓度   irr_h: (T,)   I: (T, nr) 反应速率

species, ipr = dg.ipr_from_fac("output/expA_BAU")
# species: 74 个物种名列表
# ipr[sp]['P'|'L'|'E'|'D'|'B'|'X'] → 六账积分 [molec/cm3]
```

### 20.2 `indicators.py`（路线无关引擎）

```python
from indicators import Mechanism, route_indicators, route_timeseries

m = Mechanism("mechanism/mechanism_3_full_chain.fac")   # 解析一次可复用

d = route_indicators("output/expA_BAU", "mechanism/mechanism_3_full_chain.fac",
                     mech=m, pm_only=False)
# → 12 指标 dict（见 §9.4）

ts = route_timeseries("output/expA_BAU", "mechanism/mechanism_3_full_chain.fac",
                      mech=m)
# → 12 条逐时序列 dict（见 §9.3 键名）
```

### 20.3 常见二次开发片段

```python
# ① 画 FNR 日循环
import matplotlib.pyplot as plt
ts = route_timeseries("output/expA_BAU", "mechanism/mechanism_3_full_chain.fac")
plt.plot(range(96, 96+len(ts["FNR"])), ts["FNR"])
plt.xlabel("hour"); plt.ylabel("FNR"); plt.axhline(2, ls="--", c="r")

# ② 转 ppb
PPB = 2.46e10
h, C, ih, I = dg.load("output/expA_BAU")
o3_ppb = C[:, 2] / PPB          # 第 3 列是 O3

# ③ 六账守恒自检
species, ipr = dg.ipr_from_fac("output/expA_BAU")
for sp in ("O3", "HNO3", "H2O2"):
    net = sum(ipr[sp].values())
    print(f"{sp:6s} 净变化 {net:.3e}")
```

---

## 附录：最小示例数据集

### A.1 最小可跑配置（3 天城市模拟）

```bash
# 1) 3 天模拟，1 小时输出
cat > input_v3/model.parameters <<'EOF'
3.0        # 天数
3600.0     # 输出步长
1.0e-6     # RTOL
1.0e3      # ATOL
EOF

# 2) 最小排放（只有 NO/CO/ISOP）
cat > input_v3/emissions.config <<'EOF'
NO      2.0e7   BOX  0.9
CO      1.0e7   BOX  0.9
ISOP    2.0e6   BOX  0.9
EOF

# 3) 无序列依赖（常数气象）
cat > input_v3/meteo.config <<'EOF'
TEMP   = 298.15
PRES   = 1013.25
RH     = 60
BLH    = 1000.0
RAIN   = 0.0
JFAC   = 1.0
EOF

# 4) 跑
make runtime_v3.dat
cd facsimile_impl && ./fac_v3 SCEN_BAU ../output/minimal ../runtime_v3.dat
```

### A.2 最小观测约束（24 点序列）

```bash
cat > input_v3/constraints/o3_demo.txt <<'EOF'
# unit: ppb
      0    30.0
   3600    35.0
   7200    42.0
  10800    55.0
  14400    68.0
  18000    82.0
  21600    95.0
  25200   104.0
  28800   108.0
  32400   105.0
  36000    98.0
  39600    88.0
  43200    75.0
  46800    60.0
  50400    45.0
  54000    35.0
  57600    28.0
  61200    22.0
  64800    18.0
  68400    15.0
  72000    13.0
  75600    12.0
  79200    14.0
  82800    20.0
EOF
printf 'TAU_MIX = 21600\nO3 = file:constraints/o3_demo.txt\n' > input_v3/bc.config
make runtime_v3.dat
```

### A.3 验证清单（跑完最小集后）

- [ ] `O3max` 在 50–150 ppb（城市量级）
- [ ] `O3_ppb` 末 24h 均值 30–80 ppb
- [ ] `FNR` 在 0.5–10
- [ ] `ipr_from_fac` 的 O3 六账净变化与浓度变化同量级
- [ ] `_ipr.csv` 中 `O3_B` 列非全零（同化生效）

---

## 与其它文档的关系

| 文档 | 回答 |
|---|---|
| [README.md](README.md) | 项目是什么、路线清单、机制对比、限制与勘误 |
| [OPERATIONS_MANUAL.md](OPERATIONS_MANUAL.md) | 架构地图、路线速查、编译边界、实验 A–I 速览、自纠自查记录 |
| **本文** | **实验怎么做（手把手八段式）+ 数据接入 + 指标解读 + 排错** |
| [figures/README.md](figures/README.md) | 9 张机制图逐图讲解 + 七机制链对比 |
