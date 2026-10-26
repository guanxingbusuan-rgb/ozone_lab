# 输入单位与换算参考手册（INPUT_UNITS）

> **适用范围**：ozone_lab 中**所有**可输入 / 可更改的数值——14 个数据配置文件（逐子量的
> 单位、缺省、换算、校验）、5 类外部数据模板、机制 `.fac` 中的速率与计量、
> 以及 Fortran 侧的硬编码 / 半硬编码系数与参数。
>
> **基线版本**：V2.3（74sp/131rx · 三路线 · 2026-10-26）
> **文档性质**：纯说明（不改代码）；每一条都对应 `src/cfgbuild.py`、
> `modules/runtime_dat_mod.f90`、`modules/chem_phys_mod.f90` 的实际行为。
> **变更单位约定时请同步更新本文档 + `OPERATIONS_MANUAL.md §2b`。**

---

## 目录

- [0. 完备性结论（先看这里）](#0-完备性结论)
- [1. 浓度单位体系（三态 + 换算因子）](#1-浓度单位体系)
- [2. 14 个数据文件逐子量单位表](#2-14-个数据文件逐子量单位表)
- [3. 外部数据模板的单位要求](#3-外部数据模板的单位要求)
- [4. 机制 .fac 的单位约定](#4-机制-fac-的单位约定)
- [5. 硬编码 / 半硬编码参数登记表](#5-硬编码--半硬编码参数登记表)
- [6. 内部换算公式汇编](#6-内部换算公式汇编)
- [7. 校验与拒绝行为清单](#7-校验与拒绝行为清单)
- [8. 已知限制与保留量](#8-已知限制与保留量)
- [9. 换算速查卡](#9-换算速查卡)

---

## 0. 完备性结论

| 类别 | 单位能力 | 结论 |
|---|---|---|
| **浓度类输入**（IC / BC / 排放 / 序列） | `cm3`（molec/cm³）· `ppb` · `ppm` 三态全通；排放额外支持 `ppb/h` | **完备**：用户可用原始单位直接落位，无需手工换算 |
| **气象 / 沉降 / 稀释 / 情景** | 物理单位固定（K/hPa/%/m/mm·h⁻¹/cm·s⁻¹） | **完备**（单位唯一，无需换算） |
| **光解** | `l`（s⁻¹）+ `m`/`n` 无量纲；JMAP 可覆盖任意槽位 | **完备** |
| **机制 `.fac`** | 速率 molec/cm³/s、计量无量纲（分数）、J 槽位 | **完备**（改动 = 重编译，见 §4） |
| **硬编码系数** | mcm_run 罐头参数 / 网格 / ATOL / 组映射 | **有限**：可改但需重编译，清单见 §5 |
| **单位换算的单一事实源** | `src/cfgbuild.py: unit_factor()` | ✓ 三处调用共用（IC / BC / 序列） |

**一句话**：**数值类输入全部支持单位声明（浓度三态 + 排放 ppb/h）；物理参数是固定约定；
只有 mcm_run 罐头案例的排放/初值/沉降是硬编码**——但它是刻意冻结的对照路线，不承担用户数据。

---

## 1. 浓度单位体系

### 1.1 三种可用单位（`unit_factor` 单源函数）

| 单位写法 | 含义 | 换算到 molec/cm³ | 使用位置 |
|---|---|---|---|
| `cm3` / `molec/cm3`（**缺省**） | 数密度 | ×1 | IC 行后缀、BC 内联后缀、序列头部指令 |
| `ppb` | 十亿分之一体积混合比 | **×2.46×10¹⁰** | 同上 |
| `ppm` | 百万分之一 | ×2.46×10¹³（=1e3·PPB） | 同上 |

写法规则：**大小写不敏感**（`.lower()` 归一）。未知单位 → **编译期拒绝**
（errs 门 / `ValueError`），不会静默按 cm³ 处理。

### 1.2 PPB 换算因子的物理含义与口径差 ⚠

```python
PPB = 2.46e10   # src/cfgbuild.py:13  (固定值)
```

该因子对应 **P = 1013.25 hPa, T = 298.15 K** 的理想气体数密度：
`M = 101325/(1.380649e-23×298.15) ≈ 2.46e25 molec/m³ = 2.46e19 molec/cm³`，
1 ppb = 1e-9 × 2.46e19 = **2.46e10**。

**运行时 M 是按实际 PRES/TEMP 计算的**（`eval_env`），所以：

| 情形 | 影响 |
|---|---|
| T = 298.15 K, P = 1013.25 hPa | 无偏差（因子精确） |
| T = 303 K（HOT 情景 +5K） | 以 ppb 写入的 IC/BC 实际混合比偏低 ~1.1% |
| 高原（P = 700 hPa） | 偏低 ~31%（因子仍按 1 atm） |

**这是刻意沿用的 MCM / AtChem2 同款约定**：化学一致性不受影响（速率的 M 依赖用真实
数密度），只有"ppb 标签"是近似。若要严格按实际 T/P 换算，请以 `cm3` 为单位写入并自行换算。

### 1.3 排放的第四种单位（V2.3 新增）

```
# unit: ppb/h          ← 文件头指令（整份文件生效）
NO   10.0   BOX  0.9   ← 10 ppb/h
```
换算：`E[molec/cm³/s] = X[ppb/h] × 2.46e10 / 3600 = X × 6.8333×10⁶`。
未写指令时缺省 `molec/cm³/s`（不变）。

---

## 2. 14 个数据文件逐子量单位表

### 2.1 `initialConcentrations.config` —— 初值（IC）

语法：`物种  值  [单位]`（单位**每行独立**，可混用）

| 子量 | 单位 | 接受值 | 缺省 | 校验 |
|---|---|---|---|---|
| `值` | 由第 3 列决定 | 任意实数 | 未列物种 = **0** | 未知物种 → errs 拒绝；未知单位 → 拒绝 |

```text
NO      5.0e8    cm3      # 数密度
O3      30.0     ppb      # 等价 7.38e11
CH4     1.8      ppm      # 等价 4.428e13
PAN     5.0e9             # 省略单位 → 按 cm3
```

> 注意：`input_v3/initialConcentrations.config` 现有 15 行只是"非零初值"清单，
> 其余 59 个物种按 0 起步（spin-up 由 5 天模拟完成）。

### 2.2 `meteo.config` —— 气象 7 量

语法：`NAME = 常数 | file:constraints/xx.txt`；**序列与常数可逐量混用**。

| 序号 | 键 | 单位 | 缺省 | 换算 / 消费 |
|---|---|---|---|---|
| 1 | `TEMP` | **K** | 298.15 | 直接用于 Arrhenius；情景 `DT` 加上（K） |
| 2 | `PRES` | **hPa** | 1013.25 | ×100 → Pa → M |
| 3 | `RH` | **%**（0–100） | 50 | ÷100 → H₂O 分压（Magnus） |
| 4 | `RADIAT` | W·m⁻² | 0 | **保留量**（见 §8，当前无消费者） |
| 5 | `BLH` | **m** | 1000 | ×100 → cm（干沉 k=vd/BLH） |
| 6 | `RAIN` | **mm/h** | 0 | 负值截零（`max(RAIN,0)`） |
| 7 | `JFAC` | 无量纲 | 1.0 | 云/气溶胶对光解的缩放（乘入 J） |

### 2.3 `model.parameters` —— 模拟控制（4 行，位置固定）

| 行 | 含义 | 单位 | 示例 | 说明 |
|---|---|---|---|---|
| 1 | 模拟时长 | **天** | `5.0` | `nsteps = days×86400/dtout` |
| 2 | 输出步长 | **秒** | `3600.0` | 只写末 24 h（`step ≥ nsteps−24`） |
| 3 | RTOL | 无量纲 | `1.0e-6` | CVODE 相对容差 |
| 4 | ATOL | **molec/cm³** | `1.0e3` | 基线绝对容差；**随后按物种名覆盖**（见 §5.4） |

### 2.4 `bc.config` —— 边界条件（观测张弛）

语法：`TAU_MIX = 秒`；`物种 = 常数 [单位] | file:路径`

| 子量 | 单位 | 缺省 | 校验 |
|---|---|---|---|
| `TAU_MIX` | **s** | 86400（1 天） | — |
| 内联常数 | cm³（默认）/ **ppb / ppm 后缀**（V2.3） | 0（=不约束该物种） | 未知单位/物种 → errs 拒绝 |
| `file:` 序列 | 头部 `# unit:` 指令 | cm³ | 非单调 / >128 点 → 拒绝 |

```text
TAU_MIX = 86400
O3  = 6.2e11            # molec/cm3
O3  = 40 ppb            # 等价写法（V2.3）
O3  = file:constraints/o3_station_BJ-PKU.txt
```

**物理语义**：`dy/dt += (C_obs(t) − y)/TAU_MIX`（nudging），τ 越小约束越强。
`TAU_MIX` 与 IPR 第 3 账的强度直接相关——OSSE 实验（`experiments/obs_osse.py`）扫 2/6/24/72 h。

### 2.5 `dilution.config` —— 稀释

| 子量 | 单位 | 缺省 | 消费 |
|---|---|---|---|
| `DIL` | **s⁻¹**（一阶） | 0 | `dy/dt −= DIL·y`；默认 1.16e-5 ≈ 24 h 通风 |

### 2.6 `drydep.config` —— 干沉降

语法：`物种  vd  on|off`

| 子量 | 单位 | 缺省 | 换算 |
|---|---|---|---|
| `vd` | **cm/s** | 不列 = 0 | `k = vd/BLH[cm]` → **s⁻¹** |
| `on/off` | 布尔 | off（省略第三列） | off 时 vd 不生效 |

参考量级：O₃ 0.8 / HNO₃ 2.0 / H₂O₂ 1.0 / HCHO 0.3 / 颗粒性前体 0.2–0.5 cm/s。

### 2.7 `wetdep.config` —— 湿沉降

语法：`物种  a  b` → `Λ = a · RAIN^b`（RAIN 单位 mm/h）

| 子量 | 单位 | 缺省 |
|---|---|---|
| `a` | **s⁻¹·(mm/h)⁻ᵇ** | 不列 = 0（无湿沉） |
| `b` | 无量纲指数 | — （典型 0.8） |

### 2.8 `emissions.config` —— 排放（V2.3 支持单位指令）

语法：`物种  基线  日型[BOX|FLAT]  振幅`

| 子量 | 单位 | 缺省 | 换算 / 校验 |
|---|---|---|---|
| `基线` | **molec/cm³/s** 或（头部 `# unit: ppb/h`）**ppb/h** | 未列物种 = 0 | ppb/h ×6.8333e6；负值 → 拒绝 |
| `日型` | 枚举 `BOX` / `FLAT` | FLAT | 其它拼写（含小写 `box`）→ **拒绝** |
| `振幅 amp` | 无量纲 ∈ [0,1] | 0 | 越界 → 拒绝 |

**日型公式**（`box_shape`，源码）：
```
BOX : 06–18 时  E₀·(0.1 + 0.9·amp·sin²(π(h−6)/12)) ；18–06 时 E₀·0.1
FLAT: E₀（全天恒定）
```

**日均推导（重要，易错）**：
```
∫₀²⁴ f dh = 12×0.1 [夜间] + ∫₆¹⁸[0.1 + 0.9·amp·sin²] dh
          = 1.2 + 1.2 + 0.9·amp×(12/2)
          = 2.4 + 5.4·amp
日均 = E₀·(2.4 + 5.4·amp)/24 = E₀·(0.1 + 0.225·amp)
```
⚠ 注意系数是 **0.225**，不是 0.45——因为 sin² 只在 12 小时窗口内生效，
其权重是 12×(1/2)/24 = 0.25（而非 0.5）。

**逐时实测表**（源码 `box_shape` 数值验证，24 个半小时中点）：

| amp | 峰值 | 谷值 | **日均** | 逐时分布（0→23 h） |
|---|---|---|---|---|
| 0.0 | 0.100 | 0.100 | **0.1000** | 全天 0.10（⚠ 不是零排放） |
| 0.5 | 0.542 | 0.100 | **0.2125** | 0.10×6 → 0.11 0.17 0.27 0.38 0.48 **0.54 0.54** 0.48 0.38 0.27 0.17 0.11 → 0.10×6 |
| 0.9 | 0.896 | 0.100 | **0.3025** | 0.10×6 → 0.11 0.22 0.40 0.61 0.79 **0.90 0.90** 0.79 0.61 0.40 0.22 0.11 → 0.10×6 |
| 1.0 | 0.985 | 0.100 | **0.3250** | 0.10×6 → 0.12 0.23 0.43 0.67 0.87 **0.98 0.98** 0.87 0.67 0.43 0.23 0.12 → 0.10×6 |

**反算基线**：若清单给的是**日均排放**，则
```
E₀ = 日均 / (0.1 + 0.225·amp)
例：日均 1.0e7 molec/cm³/s、amp=0.9 → E₀ = 1.0e7/0.3025 = 3.31e7
```
若清单给的是**峰值**，则 `E₀ = 峰值/(0.1 + 0.9·amp)`。**务必先确认清单口径**。

#### 2.8b 排放量级参考（判断自己的数是否合理）

**当前配置的实际量级**（`input_v3/emissions.config`，V2.2 校准档）：

| 物种 | E₀ [molec/cm³/s] | 日型/amp | 峰值 [molec/cm³/s] | 折算日均 | 折算峰值 |
|---|---|---|---|---|---|
| NO | 2.157e7 | BOX/0.9 | 1.96e7 | 0.95 ppb/h | 2.87 ppb/h |
| OLE | 2.134e7 | BOX/0.9 | 1.94e7 | 0.94 ppb/h | 2.84 ppb/h |
| ALK | 2.622e7 | BOX/0.9 | 2.39e7 | 1.16 ppb/h | 3.48 ppb/h |
| ARO1 | 1.199e7 | BOX/0.9 | 1.09e7 | 0.53 ppb/h | 1.59 ppb/h |
| ISOP | 2.037e6 | BOX/0.9 | 1.85e6 | 0.09 ppb/h | 0.27 ppb/h |

**典型城市量级参考**（用于判断你的清单是否离谱）：

| 物种 | 城市日均排放 | 峰值（早高峰） | 说明 |
|---|---|---|---|
| NOx（NO+NO₂） | 1–10 ppb/h | 5–50 ppb/h | 交通主导；中国特大城市可到 ~50 |
| CO | 10–50 ppb/h | 50–200 ppb/h | 不完全燃烧 |
| 烷烃 ALK | 1–10 ppb/h | 3–30 ppb/h | LPG/天然气泄漏 + 尾气 |
| 烯烃 OLE | 1–5 ppb/h | 3–15 ppb/h | 尾气（烯烃反应性最强） |
| 芳烃 ARO | 0.5–5 ppb/h | 2–15 ppb/h | 溶剂 + 尾气 |
| 异戊二烯 ISOP | 0.1–2 ppb/h（乡村 5–20） | 0.3–6 / 15–60 | 温度/光照驱动，正午峰 |
| SO₂ | 0.5–10 ppb/h | 2–30 ppb/h | 燃煤/工业 |
| DMS | <0.1 ppb/h（海洋 0.05–0.5） | 恒定 | 浮游植物 |

**自检方法**（稳态法）：
```
SS 浓度 ≈ E / (损失率)
NO:  E=5e6 molec/cm3/s, 对 O₃ 滴定寿命 ~55 s → SS ≈ 2.8e8 ≈ 11 ppb  ✓ 合理
SO2: E=4e5, 干沉 k=0.3cm/s /1e5cm = 3e-6 s^-1 → SS ≈ 1.3e11 ≈ 5.4 ppb ✓
若算出 SS > 1e12（40 ppb）且无强汇 → 排放偏高或漏了汇
```

### 2.9 `scenarios.config` —— 情景因子

语法：`NAME = E_NOx*<f> E_VOC*<f> E_CO*<f> E_CH4*<f> JFAC*<f> [DT+<ΔK>]`

| 子量 | 单位 | 取值 | 语义 |
|---|---|---|---|
| `E_NOx / E_VOC / E_CO / E_CH4` | 无量纲**乘子** | 任意正数 | 乘到物种组基线排放（组映射见 §5.5） |
| `JFAC` | 无量纲乘子 | 任意正数 | 乘入全部光解 J |
| `DT` | **K（绝对温差）** | 如 `+5` | 加到 TEMP（`eval_met` 内 `+dt_scen`） |

### 2.10 `photolysis.config` —— 光解方案

| 子量 | 单位 | 说明 |
|---|---|---|
| `MODE` | 枚举 `MCM` / `DIRECT` | **R16 起两模式语义不同**：`MCM` = `J=jf·l·cosᵐ·exp(−n secχ)`（天顶角调制）；`DIRECT` = `J=jf·l`（不调制，白天恒定取 l、夜间 0）。修前两者逐位相同（死开关） |
| `LAT` | **度**（deg，北纬为正） | SZA 计算 |
| `DOY` | 一年中的第几天（1–365） | 赤纬计算 |
| `JMAP` | 映射表 | `位置=槽位` 或 `位置=D<正午值>`；D 值的单位是 **s⁻¹** |

J 公式（源码 `eval_J`）：
```
decl = 23.45°·sin(2π(DOY−81)/365)
cosχ = sin(LAT)sin(decl) + cos(LAT)cos(decl)cos(时角)
J = JFAC · l · cosχ^m · exp(−n/cosχ)      (cosχ ≤ 1e-3 → J = 0)
```

### 2.11 `photolysis-rates_MCM` —— MCM 光解参数表

列：`j  l  m  n  name  tau`

| 列 | 单位 | 说明 |
|---|---|---|
| `j` | 槽位号（无量纲整数） | 机制 `J<j>` 引用；1–61，本机制用到 1–41 |
| `l` | **s⁻¹** | 天顶角 0 时的光解速率（D 记法，cfgbuild 自动转 E） |
| `m` | 无量纲 | cosχ 幂指数 |
| `n` | 无量纲 | 光程（secχ）衰减指数 |
| `name` | 字符串 | 仅注释 |
| `tau` | 整数标志 | MCM 表保留列，本仓库不消费 |

### 2.12 `constraints/*.txt` —— 时间序列（气象 / BC 共用）

格式：**两列** `t[s]  值`。头部 `# unit: cm3|ppb|ppm` **只能用于 BC 浓度序列**；
气象序列各量用原生单位（见 §2.2），**禁止**出现该指令行。

> ⚠ **这是硬拒绝，不是警告**：`# unit:` 的唯一作用是"乘以浓度换算因子"。
> 气象文件误写 `# unit: ppb`，298.15 K 会被静默变成 7.33e12（×2.46e10）。
> `cfgbuild` 现对气象序列调用 `read_series(allow_unit=False)`，遇到即报错；
> `tests/functional_test.py` 有对应负向用例守卫。模板见
> [`templates/met_series.template`](templates/met_series.template)。

| 子量 | 单位 | 校验 |
|---|---|---|
| `t` | **秒，自模拟起点**，严格递增 | 非严格递增 → 拒绝；外延取端点值 |
| `值` | 该量原生单位（TEMP: K / PRES: hPa / RH: % / BLH: m / RAIN: mm·h⁻¹ / JFAC: 1 / **RADIAT: W·m⁻²，但当前无消费者** / BC: 由 `# unit:` 决定） | ≤128 点（Fortran `ser_t(8,128)`）；气象+BC 共享 8 条池 |

**现有序列文件**：`met_BLH.txt`/`met_RAIN.txt`（默认 7 点）、
`met_{TEMP,PRES,RH,BLH,RAIN,JFAC}_demo.txt`（121 点逐时样例，实验 M 用）、
`bc_{O3,HCHO,NO2_sat}_demo.txt`（24 点 BC 约束样例）。

### 2.13 `mechanism_3_full_chain.fac` —— 见 §4

### 2.14 `emissions`/`drydep`/… 之外的：`input_v3/constraints/` 目录属 2.12

---

## 3. 外部数据模板的单位要求

| 模板 | 原始单位 | 需做的事 | 换算公式 |
|---|---|---|---|
| `surface_o3_station` | ppb | 头部 `# unit: ppb` + 两列 t[s]/值 | 无（cfgbuild 换算） |
| `met_series` | **原生单位**（K / hPa / % / m / mm·h⁻¹ / 1） | **不写** `# unit:` 行 | 无（原样进方程） |
| `surface_voc_station` | ppb（PAMS-57 / TO-15） | 同上，一站一物种 | 物种名须映射到机制物种 |
| `surface_aqi_station` | **μg/m³**（CO 为 mg/m³） | **离线拆分+换算**（模板 §单位换算表） | `ppb = μg/m³ × 24.45 / M`（25 °C,1 atm）：O₃ ×0.509（M=48）、NO₂ ×0.531（M=46）… |
| `satellite_column` | **molec/cm²**（TROPOMI 柱） | **用 `src/sat2surf.py` 转换**（引擎不读柱；工具含两法互证+QA+抽稀） | `ppb = 柱[molec/cm²] / BLH[cm] / 2.46e10`；⚠ 旧式 `柱/(2.15e25×混合层分数)` 缺 ×1e9（差 9 量级），已修正并被 FT-SAT 测试固化 |
| `emissions_inventory` | 任意 → molec/cm³/s 或 ppb/h | 直接对应 emissions.config | §2.8 |

> 模板本身**不被程序读取**（`templates/` 零代码引用）——它是格式合同 + 换算表，
> 数据落位到 `input_v3/constraints/` 或 `emissions.config` 后才进入管线。

---

## 4. 机制 `.fac` 的单位约定

| 成分 | 单位 | 说明 / 校验 |
|---|---|---|
| 命名速率 `K... = A·EXP(-E/T)` | **molec/cm³/s**（双分子）/ s⁻¹（单分子/光解） | Arrhenius 温度用 K；`@n` 为幂律 |
| 三体/Troe | 同上；`M` = **molec/cm³** | 运行时 `eval_env` 由 PRES/TEMP 算 M |
| `H2O` / `O2` | **molec/cm³** | O₂ = 0.2095·M；H₂O 由 Magnus 公式 |
| 光解 `J<n>` | s⁻¹（由 §2.10/2.11 提供） | 槽位必须被 JMAP/表覆盖，否则编译期拒绝 |
| 化学计量系数 | **无量纲，允许分数** | 元素审计强制配平（131/131） |
| `VARIABLE` 物种名 | 无单位 | 重复 → mechc 拒绝 |
| `FORMULA` 元素组成 | 整数元素计数 | 供审计 + NOy 自动判定 |

**修改方式**：改 `.fac` → `python3 audit/element_audit.py <fac>` → `make v3`（重编译）。
**注意**：`O1D+H2O`、`O1D+M`、`RO+O2` 等在 .fac 中已把 H2O/M/O2 **显式乘进速率表达式**
（`KO1DH2O = 2.2D-10*H2O`），不要重复乘。

---

## 5. 硬编码 / 半硬编码参数登记表

### 5.1 mcm_run 罐头参数（`mcm_city/mcm_driver.f90`，改需重编译）

| 行 | 参数 | 值 | 单位 |
|---|---|---|---|
| 35 | `TAU` | 8.64e4 | **s**（24 h 背景张弛） |
| 36 | `O3_BG` | 6.2e11 | molec/cm³（≈25 ppb） |
| 36 | `CO_BG` | 1.5e12 | molec/cm³（≈61 ppb） |
| 36 | `CH4_BG` | 4.5e13 | molec/cm³（≈1.8 ppm） |
| 56–60 | `E(iNO)` 等 5 物种 | 2e7 / 1e6 / 5e6 / 5e3 / 2e5 | molec/cm³/s（NO 处注释"20 ppb/h"） |
| 54 | `vdk` | O₃ 8e-6 / HNO₃ 2e-5 / HCHO 3e-6 | **s⁻¹（已除 BLH 的一阶速率）** |
| 179–181 | `Y(iO3)` 等 8 物种初值 | 3e11 / 5e8 / 2e8 / 2e12 / 4.5e13 / 2e9 / 1e8 / 1e6 | molec/cm³ |
| 184 | `RTOL` / `ATOL` | 1e-5 / 1e4（O₃ 1e9, CH₄ 1e11, CO 1e9） | 无量纲 / molec/cm³ |

> ⚠ `vdk` 的 8e-6 s⁻¹ ≡ 0.8 cm/s ÷ 1000 m（1e5 cm）——**与 fac_v3 的 `vd/BLH` 语义等价**，
> 但存的是除完的结果，不是 cm/s。

### 5.2 fac_batch_omp 网格与容差

| 参数 | 值 | 单位 / 说明 |
|---|---|---|
| 网格 `fs(i)` | `0.3 · 1.5^(i−1)`，i=1..ng | 无量纲排放因子；ng≤32 |
| jfac 列表 `fjs` | CLI 逗号分隔 | 无量纲（≤16 项） |
| ΔT 列表 `dts` | CLI 逗号分隔 | **K**（≤16 项） |
| RTOL/ATOL | 来自 runtime.dat，ATOL 经 `tune_atol_by_name` | §5.4 |

### 5.3 fac_v3 数值控制

| 参数 | 值 | 单位 |
|---|---|---|
| `MAX_NSTEPS` | 200000 | 步 |
| `MAX_STEP` | 100（重启后 10） | **s** |
| 输出窗口 | `step ≥ nsteps−24` | 末 24 个输出步 |

### 5.4 ATOL 按名调优表（`chem_phys_mod.tune_atol_by_name`，单源）

| 物种 | ATOL [molec/cm³] |
|---|---|
| CH₄ | 1e10 |
| CO | 1e9 |
| O₃ | 1e8 |
| HCHO / PAN / MPAN / ISOP / SO₂ / CH₃CHO / MVK / MACR / MGLY / RONO₂ / OLENIT / TERPNIT / TERPOOH | 1e7 |
| DMS / MSA / TERP | 1e6 |
| TERPCHO | 1e8 |
| 其余 | 沿用 `model.parameters` 第 4 行（默认 1e3） |

### 5.5 物种→情景组映射（`runtime_dat_mod.build_groups`，单源）

| 组 | 物种清单 |
|---|---|
| `E_NOx` | NO, NO₂ |
| `E_VOC` | HCHO, CH3CHO, ISOP, DMS, MVK, MACR, MSA, SO₂, ALK, OLE, ARO1, ARO2, KET, **TERP** |
| `E_CO` | CO |
| `E_CH4` | CH₄ |

> 新增排放物种须在此表加一行，否则情景缩放对它不生效（**唯一半硬编码的"数据映射"**）。

### 5.6 其它常量

| 位置 | 常量 | 单位 / 说明 |
|---|---|---|
| `src/cfgbuild.py` | **`PPB = 2.46e10`（唯一定义点）** | molec/cm³ per ppb（§1.2） |
| 其余 Python 引用 | `from cfgbuild import PPB`（R13 起） | sat2surf / 7 实验脚本 / diagnostics 均引用此单源 |
| Fortran 侧 | 无 PPB 常量（输入已由 cfgbuild 换算） | — |
| `cfgbuild` | `PPB/3600 = 6.8333e6` | ppb/h → molec/cm³/s |
| `eval_met` | `dBLH` ±600 s 中心差分 ÷1200 | **cm/s** |
| `eval_env` | `es = 6.112·exp(17.67(T−273.15)/(T−29.65))` | hPa（Magnus） |
| `eval_met` | BLH 输入 m → 内部 ×100 | cm |
| AQI 模板 | `24.45/M` | μg/m³ → ppb（25 °C, 1 atm） |
| 卫星模板 | `2.15e25 molec/cm²` | 标准大气柱 |

---

## 6. 内部换算公式汇编

```text
# 数密度（eval_env, PRES 由 hPa×100 → Pa；1e-6 由 m^-3 → cm^-3）
M    = 6.022e23 · P[Pa] / (8.314 · T[K]) · 1e-6          [molec/cm3]
O2   = 0.2095 · M
es   = 6.112 · exp(17.67(T−273.15)/(T−29.65))             [hPa]   (Magnus)
e    = RH/100 · es
H2O  = 6.022e23 · e[hPa]·100 / (8.314·T) · 1e-6           [molec/cm3]

# 单位换算（cfgbuild.unit_factor）
cm3 → ×1 ；ppb → ×2.46e10 ；ppm → ×2.46e13
ppb/h → ×(2.46e10/3600) = ×6.8333e6                       [molec/cm3/s]

# 光解（eval_J）
decl = 23.45°·sin(2π(DOY−81)/365)
cosχ = sin(lat)sin(decl) + cos(lat)cos(decl)cos(2π(h−12)/24)
J    = JFAC · l · cosχ^m · exp(−n/cosχ)   [s^-1]   (cosχ≤1e-3 → 0)

# 排放日型（box_shape）
BOX  = 0.1 + 0.9·amp·sin²(π(h−6)/12)   (6≤h≤18)，否则 0.1
FLAT = 1.0
日均因子 = 0.1 + 0.225·amp   (见 §2.8 推导: sin² 权重 12×0.5/24)

# 干沉降（phys_tend）
kdry = vd[cm/s] / BLH[cm]                                  [s^-1]
kwet = a · RAIN[mm/h]^b                                    [s^-1]

# BC 张弛
dy/dt += (C_obs(t) − y) / TAU_MIX − DIL·y

# BLH 膨胀稀释（仅 dB>0）
dy/dt −= y · (dBLH/BLH)     dBLH[cm/s], BLH[cm]
```

---

## 7. 校验与拒绝行为清单

| 检查 | 位置 | 行为 |
|---|---|---|
| 未知物种（IC/BC/排放/干沉/湿沉） | cfgbuild | errs → `exit(1)` |
| 未知单位（任何通道） | `unit_factor` | errs / ValueError → 拒绝 |
| 序列非严格递增 / >128 点 | read_series | 拒绝 |
| 排放负基线 / 非法日型 / amp ∉ [0,1] | cfgbuild | 拒绝 |
| 光解槽位无参数（防静默零光解） | cfgbuild 光解完备性校验 | 拒绝 |
| BC 行 token 数 ≠ 1+2·ns（防行截断） | runtime_dat_mod | **FATAL error stop** |
| 机制符号重复（物种/速率名） | mechc.py（R9） | SystemExit |
| 生成码与机制物种/反应数不一致 | FT0（functional_test） | assert 失败（提示 make v3） |
| 元素不守恒 | element_audit | FAIL（非零退出） |

---

## 8. 已知限制与保留量

| # | 项 | 说明 |
|---|---|---|
| 1 | **`RADIAT` 无消费者** | meteo.config 接受并读入第 4 量，但当前无代码使用（保留给 B2 在线光解 `modules/rad_photol.f90` 的集成位）。写入不影响模拟。**R16 起：`RADIAT ≠ 0` 时 cfgbuild 打印明确警告**，不再静默。 |
| 2 | **PPB 因子固定 1 atm/298.15 K** | 非标准 T/P 下 ppb 标签有系统偏差（§1.2）。MCM/AtChem2 同款约定。 |
| 3 | **排放无逐时序列** | 只有 BOX 解析式 / FLAT 两种时间形（唯一"数据进不来"的口）；需要逐时排放请改 §5.6 拓展或离线拟合 amp。 |
| 4 | **排放无源/部门维度** | 数据模型是"物种→(基线,日型,振幅)"一维；多部门清单需离线按物种加和。 |
| 5 | **mcm_run 全部排放/IC/沉降硬编码** | 对照路线刻意冻结（§5.1）；用户数据请走 fac_v3。 |
| 6 | **`photolysis-rates` 的 tau 列不消费** | MCM 表保留列。 |
| 7 | **情景组映射需手改代码** | 新排放物种须在 `build_groups` 加行（§5.5）。 |
| 8 | **序列总数 ≤8（气象+BC 共享池）** | `NSER_MAX=8`；超限编译期拒绝。 |
| 9 | **极痕量物种可出现负值** | O/O₁D 实测 min = −4e−4 / −1.5e−71（相对峰值 <1e-7）。BDF 不保正性是正常行为（10b 坑 #5），**不施加 `max(y,0)`**（会引入 RHS 折点，曾致 15× 步数惩罚）。 |
| 10 | **CSV 三位指数（历史产物）** | R13 前输出用 `ES13.5`，\|指数\|≥100 时静默丢 E（如 `3.07932-174`）。现已改 `ES15.5E3`；旧文件用 `diagnostics/csvfloat.fix_fortran_exp()` 解析。 |
| 11 | **无风速/风向/平流/对流** | 0-D 箱体**结构上没有**这些通道：模型无 x/y/z 坐标、无输运算子，风速风向无处可用（不是"没实现"，是"维度不存在"）。水平输送由 **`dilution.config` 的 `DIL`**（一阶通风率 [s⁻¹]）集总替代，垂直混合由 **BLH** 表达，湍流扩散并入 BLH 膨胀项。实测 `DIL` 是强控制量：0 → O3max 157.5 ppb，1.16e-5 → 101.6，5e-5 → 55.6。需要真实风场驱动的平流/对流 → 必须上三维模式（CMAQ/GEOS-Chem/WRF-Chem），0-D 箱体无法承载。 |

---

## 9. 换算速查卡

```text
1 ppb        = 2.46e10 molec/cm3          (1 atm, 298.15 K)
1 ppm        = 2.46e13 molec/cm3
1 ppb/h      = 6.833e6 molec/cm3/s
1 molec/cm3  ≈ 4.07e-11 ppb
25 ppb O3    = 6.15e11 molec/cm3
1.8 ppm CH4  = 4.43e13 molec/cm3
1 μg/m3 O3   = 0.509 ppb      (M=48)
1 μg/m3 NO2  = 0.531 ppb      (M=46)
1 μg/m3 SO2  = 0.382 ppb      (M=64)
1 μg/m3 CO   = 0.873 ppb      (M=28, 且原始单位 mg/m3)
1 m BLH      = 100 cm
vd 0.8 cm/s @ BLH 1000 m = 8.0e-6 s^-1
DIL 1.16e-5 s^-1 ≈ 24 h 通风一次
```

---

## 相关文档

- [`OPERATIONS_MANUAL.md §2b`](OPERATIONS_MANUAL.md) —— 四类输入/参数/硬编码总表（本手册的摘要版）
- [`mechanism/README.md`](mechanism/README.md) —— 机制逐反应速率与配平方法论
- [`templates/`](templates/) —— 五类外部数据格式合同（各含单位换算表）
- [`README.md §4`](README.md) —— 输入系统 12 类格式详表 + 处理链
