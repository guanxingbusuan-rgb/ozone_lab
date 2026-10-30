# EKMA-FNR/ —— S2 初始浓度正交扰动响应面（HCHO × NOx）

> **一句话**：固定机制/气象，**关闭排放**，把**初始 HCHO 浓度**与**初始 NOx（NO+NO₂ 同因子，
> 保持 NO₂:NO = 3:1）** 按 **0%→1000%（100% 步长，默认）** 正交扰动，跑 11×11=121 个独立 0-D 箱
> （每箱 1 天纯初值化学），得到 O₃max 响应面 + FNR 等值线 + 屋脊。
> 经典 Dodge 1977"初始浓度空间 EKMA"口径；`--pmax/--step` 可改范围与步长。
>
> **现存归档**（以 `output/run_meta.txt` 为准）：**2026-10-30 完整轮** —— **0–1000%/100%（11×11=121）**，
> 1 天、排放关闭；轴基线 NO 1 / NO₂ 3 / HCHO 6 ppb + **事件情景**汇与碳（`DIL=0`、BLH×2、
> CO 500、CH₃CHO 20、ISOP 20 ppb）；**失败 0 / 降级 0 / 氮守恒违例 0**，
> O₃max **34.14–119.21 ppb**，基线点 (0,0) = 34.14 ppb，墙钟 **64.7 s**（4 线程；随负载波动，非回归判据）；
> 屋脊 NO₂ +500%→+1000%（**末列触顶**，见 §8.6）。
> 本文其余章节的数值叙述与这一轮同源。

---

## 1. 实验定义（用户确认口径）

| 项 | 内容 |
|---|---|
| 类型 | **初始浓度扰动**（IC-based EKMA，Dodge 1977） |
| 执行器 | `facsimile_impl/fac_batch_omp`（IC 模式，第 8 参 `IC`、第 9 参 **`NOX:HCHO`**、第 10 参 `LIN`） |
| 固定条件 | ① 机制 V2.5+（冻结）② 本目录气象 `config_override/meteo.config`（5 天晴天序列的第 1 天；**BLH 已 ×2** → 夜间 ~346 m、午后峰 ~3.77 km）③ 其余配置同在 `config_override/` ④ **事件情景假设**：`DIL = 0`（无水平通风）+ 固定物种 CO 500 / CH₃CHO 20 / ISOP 20 ppb |
| **排放** | **关闭**（emissions.config 清空 → 纯初值化学） |
| **扰动轴 1（x，分母）** | **初始 NOx** = `ic1` 走 **NOX 组掩膜**：**NO 与 NO₂ 同因子缩放 → NO₂:NO 恒为 3:1**（基线 NO₂ 3 / NO 1 ppb） |
| **扰动轴 2（y，分子）** | **初始 HCHO**（基线 6 ppb）——只缩 HCHO |
| 不扰动 | O₃（因变量）、H₂O₂/CO/CH₄/CH₃CHO/PAN/ISOP/MVK/MACR/MGLY/SO₂/DMS（保持基线） |
| 扰动范围 | **0%（无扰动）→ 1000%（+1000%）**，步长 **100%**（默认）⇒ 11 点/轴；NO₂ 3→33 ppb、NO 1→11 ppb、HCHO 6→66 ppb |
| 案例数 | 11 × 11 = **121**（默认），每案例 1 天、spinup=0 |
| 输出量 | O₃max（全天最大）+ FNR_mean（HCHO/NO₂ 窗均比）+ 24 列 |

**为什么 1 天 + 关排放**：HCHO/NO₂ 是短寿命物种，5 天排放驱动会完全洗脱初值影响
（S2 退化为 S1 的翻版）；1 天纯初值化学才能让"初始浓度 → O₃"的因果链成立。

**为什么 NOx 轴（纵轴）是 NOx 组而不是只缩 NO₂**：只缩 NO₂ 而 NO 固定在 5 ppb 时，总 NOx 由 NO 底数主导，
低 NO₂ 基线（1 ppb）会让整张面钉在 NOx-limited 一侧、窗内不出屋脊（2026-10-29 实测，§8.6）。
改为 NO+NO₂ 同因子缩放后，NOx 总量真正随轴变化，屋脊、滴定支与 FNR=1/2/4 分界全部落进窗内。

---

## 2. 文件与用法

```
EKMA-FNR/
├── run.py      # ① 读 config_override → FNR 规则 → cfgbuild → n×n 案例 → **全部 CSV 导出**
├── plot.py     # ② 只出图：响应面 + FNR 等值线 + 屋脊标记 + 切片 PNG
├── README.md
├── config_override/                 # 已建立的持久编辑目录（11 个顶层文件 + 17 个 constraints 文件）
│   ├── initialConcentrations.config # 初值基线：NO₂/HCHO 两轴及其余固定物种
│   ├── meteo.config                 # 气象配置，引用 constraints/ 序列
│   ├── photolysis.config / photolysis-rates_MCM
│   ├── model.parameters             # 步长/容差等；天数和 spinup 由实验强制
│   ├── drydep.config / wetdep.config / dilution.config
│   ├── bc.config / emissions.config / scenarios.config
│   └── constraints/                 # 17 个数据文件，含晴天气象 6 条序列
└── output/                          （逐文件数据字典见 §4；**无配置快照目录**）
    ├── runtime_fnr.dat              ① 本轮有效输入记录（cfgbuild 解析出的 IC/EMIS/MODEL/MET 值）
    ├── fnr_ic_base.csv              ② 本次基线初值存档（含缩放归属）
    ├── ekma_fnr_grid_raw.csv        ③ fac_batch_omp 原始 24 列输出
    ├── ekma_fnr_grid_raw.csv_ipr.csv ③a 侧车: 111 物种 × 6 过程账 窗均 [molec/cm3/s]
    ├── ekma_fnr_grid_raw.csv_irr.csv ③b 侧车: r1..r203 逐反应 IRR 窗均 [molec/cm3/s]
    │        （⚠ 2026-10-29 已删除：旧版 batch 坐标列写错；重跑 run.py 自动重建，格式正确）
    ├── ekma_fnr_grid_raw.csv_state.csv ③c 侧车: 实验J 状态量均值对（本实验未消费，备查）
    ├── ekma_fnr_grid.csv            ④ 规整主表（绘图唯一输入）
    ├── ekma_fnr_ridge.csv           ⑤ 屋脊轨迹（run.py 产出）
    ├── ekma_fnr_derived.csv         ⑥ 派生表 PR/FNR_init（run.py 产出）
    ├── ekma_fnr_S2.png / ekma_fnr_slices.png  ⑦ 图件
    └── run_meta.txt                 ⑧ 运行元数据（含"有效输入记录"指向 runtime_fnr.dat）
```

`config_override/` 是**从历史运行输入建立、长期保留的持久输入目录**：`run.py` 只读它，
从不改写它、也不会删除它。
现有基线可直接运行，**不必先准备真实数据**；日常只编辑这里的配置及 `constraints/` 数据，
不改共享层，也不把 `output/` 当作输入编辑区。替换数据见 §8–§9。

**初次准备（与日常运行分开）**：在完整仓库及运行环境中预先构建共享内核；若已构建可跳过。
以下路径以当前仓库位置为例，迁移完整仓库后相应调整：

```bash
cd /Users/wangling/CodeProject/LearnFortran_Projects/ozone_lab
make all                                 # 初次共享构建，不是日常步骤
```

**日常运行**（在 `EKMA-FNR/` 目录）：

```bash
cd /Users/wangling/CodeProject/LearnFortran_Projects/ozone_lab/EKMA-FNR
python run.py                            # ① 默认 0–1000% / 100%，11×11=121 案例
python plot.py                           # ② 单独出图
```

每次完整运行都会自动执行 `cfgbuild` 生成本实验的配置运行时，**不等于重新 `make` 编译共享内核**。
它在系统临时目录里装配配置、再由 `cfgbuild` 编译成运行时（见 §3），**不在 `output/` 留配置快照**；
每轮会更新运行时、基线存档和 CSV；图件需另跑 `plot.py` 才更新。
**重跑（包括冒烟）前，先另存需要保留的整个 `output/`**，避免旧图与新 CSV 混用。
运行仍依赖整个仓库的共享组件，不能把本子目录单独搬走运行（见 §6）。

可选命令（同样在 `EKMA-FNR/` 目录；不是必须接着运行的步骤）：

```bash
python run.py --pmax 500 --step 50         # 自定义网格（例：11×11，0–500%）
python run.py --from /absolute/path/to/ekma_fnr_grid_raw.csv  # 复用 raw 导出，不跑模型
```

`--from` 示例路径须替换为实际 raw 文件的**绝对路径**；实现中的相对路径按仓库根解析，
不是按当前 `EKMA-FNR/` 目录解析。它不装配配置、不执行 cfgbuild、不跑模拟；出图仍需另跑 `python plot.py`。

---

## 3. 输入装配与"有效输入记录"（为什么 output 里没有配置快照）

每次**完整运行**（非 `--from`）的配置链为：

```
本目录 config_override/（唯一配置输入）
    → 临时装配目录（mkdtemp，加 FNR 规则）→ cfgbuild → output/runtime_fnr.dat → fac_batch_omp
```

**装配目录是必需的**：`cfgbuild` 需要一个"一整套配置"的目录作输入，而实验口径（排放清空、
1 天/0 spinup）**不能写回你正在编辑的 `config_override/`**。因此它在系统临时目录里完成，
编译成功后立即删除（编译失败会回显 cfgbuild 的错误输出）。`output/` 不再出现 `in_fnr/`，
避免"看起来像输入入口、实际每轮被重建"的误导。

**本轮实际用了什么输入，看两处：**

| 记录 | 内容 |
|---|---|
| `output/runtime_fnr.dat` | cfgbuild 的解析结果——`IC` 向量、`EMIS` 行（全 0 = 排放关闭）、`MODEL 1 3600 …0`、`MET/PHOT` 与逐时序列的**数值本体**都在里面；这是唯一完整的"有效输入"事实源 |
| `config_override/` + `run_meta.txt` | 你编辑的输入原文，以及本轮的网格、耗时、机制 md5、`输入目录`、`FNR 基线`、`有效输入记录` 路径 |

| 配置项 | 装配时的规则 | 生效位置 |
|---|---|---|
| `initialConcentrations.config` | **取自 `config_override/`（必须存在，无内置回退）** | `runtime_fnr.dat` 的 `IC` 行 + 存档 `fnr_ic_base.csv` |
| `emissions.config` | **强制清空**，输入目录里的排放值不会生效 | `runtime_fnr.dat` 的 `EMIS` 全 0 |
| `model.parameters` | 取自输入目录，仅强制天数 **1**、SPINUP **0** | `runtime_fnr.dat` 的 `MODEL 1 3600 1e-06 1000 0` |
| 其余配置与 `constraints/`（17 文件） | 取自 `config_override/`，按相对路径整体装配 | `runtime_fnr.dat` 的 `MET/SERIES/PHOT/DRY/WET/DIL/BC` 段 |

复现历史轮次：把归档的 `config_override/`（或当初那套配置文件）放回本目录 `config_override/`，
用同一机制、内核与网格重跑。
`runtime_fnr.dat` 可用来**核对**当时生效的数值，但它不是 cfgbuild 的输入，不能直接"重放"。

> `runtime_fnr.dat` 每次完整运行由本次装配重新编译，是本轮有效输入的唯一记录。
> EKMA-emission 同构（`output/runtime_emission.dat`）；两个实验各自操作、各自输出，
> 但仍共享仓库组件。FNR 另施加初值选择、排放关闭和 1 天/0 spinup 规则。

---

## 4. `output/` 数据字典（逐 CSV 逐列）

### 4.1 `fnr_ic_base.csv` —— 基线初值存档（随基线来源更新）

| 列 | 含义 |
|---|---|
| `species` | 物种名（IC 里非零的物种） |
| `value_cm3` | 基线初值 [molec/cm³]（ppb 写入时自动 ×2.46e10 换算） |
| `unit_as_written` | 用户文件里的原始单位（cm3 / ppb / ppm） |
| `scaling` | 缩放归属：`perturbed by ic1 (NOX: NO+NO2 same factor)`（NO 与 NO₂ 同因子）/ `perturbed by ic2 (HCHO)`（×HCHO 因子）/ `fixed baseline`（不缩放，含 O₃） |

表头注释含 `基线来源:` 行，当前归档写 `config_override/initialConcentrations.config`
（`run.py` 的内置 `IC_BASE` 常量已移除，初值只可能来自本地输入目录）。
用途：审计"哪些物种被扰动、哪些固定"；与 `config_override/initialConcentrations.config`
及 `runtime_fnr.dat` 的 `IC` 行互校。

### 4.2 `ekma_fnr_grid_raw.csv` —— fac_batch_omp 原始输出（1 表头 + n² 行；当前归档 11×11=121）

> ⚠ **列名陷阱**：前两列沿用了排放扫描的历史列名 `fnox/fvoc`——**在 IC 模式下实义是
> IC 因子**：`fnox` = **NOx 组缩放因子**（ic1：NO 与 NO₂ 同乘）、`fvoc` = HCHO 缩放因子（ic2）。> 读表时勿按字面当成排放。

| 列 | 类型 | 含义（IC 模式语义） |
|---|---|---|
| `fnox` | 输入因子 | **初始 NOx 缩放因子**（NO 与 NO₂ 同乘），当前归档 1.00–11.00（= 1 + pert_no2_pct/100） |
| `fvoc` | 输入因子 | **初始 HCHO 缩放因子**，当前归档 1.00–11.00（= 1 + pert_hcho_pct/100） |
| `jfac` / `dt` | 输入因子 | 恒 1.0 / 0.0（本实验未启用光解衰减/温度维） |
| `O3max` [ppb] | 响应 | 全 24 h 模拟窗内 O₃ 最大值（**主图色标 / 本实验主响应量**） |
| `O3mean` [ppb] | 响应 | 全窗 O₃ 均值 |
| `HCHOmax` [ppb] | 响应 | 窗内 HCHO 峰值（**演化后浓度**，非初始值） |
| `NO2mean` [ppb] | 响应 | 窗内 NO₂ 均值（演化后；与初始值不同——光稳态重分配） |
| `FNR_mean` | 判据 | 窗均 HCHO/NO₂ 比（Sillman 1995；<1 VOC-limited / 1–2 过渡 / >2 NOx-limited） |
| `H2O2HNO3` | 判据 | 窗均 H₂O₂/HNO₃（>0.2 → NOx 控制区） |
| `O3NOy` / `HCHONOy` / `O3NOz` | 备用判据 | O₃/NOy、HCHO/NOy、O₃/NOz（占列，不进主判据） |
| `NOy_ppb` | 状态量 | 总活性氮（25 物种加权，N₂O₅ 计双 N） |
| `LnOxQ` | 过程量 | L_NOx/Q（Kleinman 2005；>0.5 → NOx 饱和） |
| `PH2O2PHNO3` | 过程量 | P(H₂O₂)/P(HNO₃) 产率比 |
| `OPE` | 过程量 | O₃ 生成效率 = P(O₃)/P(NOz) |
| `*_pm`（6 列） | 判据(午后) | 上列对应量的 **12–18 h 午后口径**（Sillman 原文窗口） |
| `nst` | 诊断 | CVODE 内步数（跨韧性重启累计；异常小 = 该案例收敛失败已降级） |

### 4.3 `ekma_fnr_grid.csv` —— 规整主表（注释头 + 表头 + n² 行；**plot.py 唯一输入**）

| 列 | 含义 |
|---|---|
| `pert_hcho_pct` | 初始 HCHO 扰动百分比 = (fvoc−1)×100，**x 轴（分子）**（范围随网格） |
| `pert_no2_pct` | 初始 NO₂ 扰动百分比 = (fnox−1)×100，**y 轴（分母）**（范围随网格） |
| `O3max_ppb` / `O3mean_ppb` | 由 raw 的 `O3max/O3mean` 改名（语义同上） |
| 其余 18 列 | 与 raw 同名同义（HCHOmax … nst） |

注释头（`#` 行）自记录：网格口径 / 两轴定义 / NOx 组缩放（NO+NO₂ 同因子）与 O₃ 不缩放声明 / 排放关闭 / FNR 定义。

### 4.4 `ekma_fnr_ridge.csv` —— 屋脊轨迹（每 HCHO 列一行；当前归档 11 行）

| 列 | 含义 |
|---|---|
| `pert_hcho_pct` | HCHO 扰动列（步长随网格；当前归档 0–1000，100%） |
| `ridge_pert_no2_pct` | **该列 O₃max 峰值所在的 NO₂ 扰动 %**（体制分界：屋脊**下方**（低 NOx）= NOx-limited，**上方**（高 NOx）= 滴定区） |
| `ridge_O3max_ppb` | 该峰值的 O₃max |

### 4.5 三个侧车（fac_batch_omp 附带产物，本实验未消费，备查/跨实验分析用）

| 文件 | 结构 | 说明 |
|---|---|---|
| `_ipr.csv` | fnox,fvoc,jfac,dt + 111 物种×6 账（`<sp>_P/L/E/D/B/X`） | 末 24 h 过程账均值；**IC 模式下 E 账（排放）恒 0**，可作 IC 生效的自检 |
| `_irr.csv` | fnox,fvoc,jfac,dt + `r1..r203` | 逐反应 IRR 窗均（机制反应行序） |
| `_state.csv` | fnox,fvoc,jfac,dt + NO2/HCHO/HNO3/H2O2 均值对 + Q/L_NOx/P 比 + O₃ + nst | 实验 J 的状态侧车；EKMA-FNR 的坐标已直接用因子，故未消费 |

### 4.7 `ekma_fnr_derived.csv` —— 派生表：PR 峰值率与初始组成比（表头注释 + 11 行）

> **🔒 固化声明**：本导出**内置于 `run.py`**（唯一 CSV 导出者；`plot.py` 只画图）；
> `run.py` 同时在 grid 主表附加 **`FNR_init` 列**（每案例初始 HCHO/NO₂ 比）——
> 均非临时脚本。概念定义见下表与 §4.7 末段。

| 列 | 含义 |
|---|---|
| `pert_hcho_pct` | HCHO 扰动列（步长随网格；当前归档 0–1000，100%） |
| `O3peak_ppb` | 该列沿 NO₂ 轴的 O₃max 峰值（屋脊点） |
| `PR_bau` | 峰值率① = O₃peak / O₃max(基线初值 0%,0%) |
| `PR_vs_baseNO2` | 峰值率② = O₃peak / O₃max(同列 NO₂=0% 端)。⚠ **0% 扰动 = 基线初值，不是零浓度**——是「从基线 NO₂ 爬到屋脊」的放大 |
| `O3_at_baseNO2_ppb` | 同列 NO₂=0% 端（=基线初值）的 O₃max |
| `ridge_pert_no2_pct` | 屋脊 NO₂ 扰动 |
| `FNR_init_ridge` | **屋脊处初始 HCHO/NO₂ 比** = (1+HCHO%)/(1+NO₂%)×H0/N0 —— VNr（排放组成比）在 IC 实验中的对应物。**H0/N0 从 `fnr_ic_base.csv` 读实际基线**（当前基线 6/3 ⇒ 系数 2.0；换基线自动跟随） |

> PR/VNr 概念并非 EKMA-emission 独有：任何规则网格响应面都有"峰值率（屋脊/基准）"
> 与"组成比（两轴因子比）"。EKMA-FNR 中组成比 = **初始 FNR**（这正是平面名称的由来）；
> 实测屋脊处 FNR_init 随 HCHO 上升（**0.33→2.00，非常数**——与排放实验的近常数临界比不同，
> 因为 1 天纯初值化学是瞬态响应而非稳态标度）。

### 4.6 非 CSV 产物

| 文件 | 说明 |
|---|---|
| `runtime_fnr.dat` | cfgbuild 编译产物（Fortran 顺序表控；由临时装配目录生成，**本轮有效输入的完整记录**） |
| `run_meta.txt` | 时间 / 模式 / 网格 / 耗时 / 机制 md5 / IC 缩放归属 / IC 基线与来源 / 排放开关 / 失败·降级 / **氮守恒违例数** / **侧车单位口径** |
| `ekma_fnr_S2.png` / `ekma_fnr_slices.png` | 响应面主图 / NO₂ 切片线图（§5） |

---

## 5. 图怎么读

| 元素 | 含义 |
|---|---|
| 横轴 | 初始 HCHO 扰动 %（0 = 基线 6 ppb；+1000% = 66 ppb）——VOC 轴 |
| 纵轴 | 初始 NO₂ 扰动 %（0 = 基线 3 ppb；+1000% = 33 ppb）——**NO 同因子（NO₂:NO=3:1），即 NOx 轴** |
| 色标 / 黑实线 | O₃max [ppb] 及其等值线 |
| **白虚线 FNR=1/2/4** | HCHO/NO₂ 比值等值线（轴是控制变量，比值不再是过原点射线而是曲线族）—— |
| **白虚线 + 其上黑实线** | **屋脊**：每个 HCHO 列上 O₃max 峰值对应的 NO₂ 扰动 |
| 黑十字 | 基线（0%, 0%）= `config_override/initialConcentrations.config` 的 NO 1 / NO₂ 3 / HCHO 6 ppb |

> **轴序约定**：与 S1（EKMA-emission）统一为 **VOC 在横轴、NOx 在纵轴**，两张图可并排读；
> `run.py` 的 CSV 列名与数据不受影响，只有 `plot.py` 的取列方向变了。

**判读要点**：
- FNR > 4 区（右下，HCHO 高 / NO₂ 低）= NOx-limited：加 NO₂ 增产 O₃；
- FNR < 1 区（左上，NO₂ 高 / HCHO 低）= VOC-limited：加 NO₂ 滴定 O₃（NO + O₃ → NO₂）；
- 屋脊 = 体制分界，其位置随 HCHO 增加而**上移**（NOx 在纵轴）。

**⚠ 本实验的口径注意**：纵轴是 **NOx 组**——`ic1` 同时缩放 NO 与 NO₂（保持 NO₂:NO = 3:1），
所以"NO₂ 扰动 %"读作"NOx 扰动 %"，NO 也按同比例变化；这与"只缩 NO₂、NO 固定"的旧口径不同，
后者在低 NO₂ 基线下会让整张面失去屋脊（见 §8.6）。

---

## 6. 与 EKMA-emission（S1）的区别与共享依赖

| | EKMA-emission（S1） | EKMA-FNR（S2） |
|---|---|---|
| 扰动对象 | 本实验 `config_override/emissions.config`（NOx/NMVOC 排放组） | 本实验 `config_override/initialConcentrations.config`（HCHO / NOx 初值） |
| 排放 | 基线排放 × 因子 | **关闭** |
| 时长 | 4 天（统计窗 = 第 4 天） | 1 天（统计窗 = 唯一那天） |
| 轴语义 | 采样自变量 = 控制量 | 采样自变量 = 控制量（**均已铺满**） |
| runtime | 本目录自编译 `runtime_emission.dat`（← 各自的临时装配目录） | 本目录自编译 `runtime_fnr.dat`（同） |

### 6.1 两个 O₃max 为什么差 2.4 倍——开放箱 vs 封闭瓶

> 数字取自两边**各自默认网格的完整运行**（S1 = 0–625%/25%、4 天/3 spinup；
> S2 = 0–1000%/100%、1 天，均 2026-10-29）。换基线或换汇项后请回读各自的
> `output/run_meta.txt` 与 `ekma_*_grid.csv`，**不要照抄数值，照抄口径**。

**最干净的证据是零扰动点本身**：S1 的 (0%,0%) = **77.12 ppb**，S2 的 (0%,0%) = **34.14 ppb**
（2.3×）；全网格峰值 S1 **281.41** vs S2 **119.21 ppb**（2.4×）。同一机制、同一晴天气象序列，
差异全部来自**实验定义（体系开放度）**，不来自扰动幅度。

| 边界条件 | EKMA-emission（S1，峰 281.4 ppb） | EKMA-FNR（S2，峰 119.2 ppb） |
|---|---|---|
| 排放 | **开**：NOx 基线 1.003 ppb/h（24h 均值，BOX 日型午后峰 ≈3×）、NMVOC 3.716 ppb/h | **关**：`emissions.config` 清空 + Fortran `esc=0` |
| 背景张弛 | **开**：`TAU_MIX = 86400 s`、`O3 = 6.2e11` = **25.2 ppb** 背景，一天换一箱 | **关**：IC 模式强制 `tau_mix = 1e30`（`facsimile_impl/fac_batch_omp.f90:566`） |
| 水平通风 | `dilution.config` 生效（参考配置 1.16e-5 s⁻¹ ≈ 24 h 一箱） | **`DIL = 0`**——本目录默认已改为密闭箱 |
| 混合层 | 晴天日循环：夜间 180 m → 午后 ~1.8 km | **同序列 ×2**：夜间 ~346 m → 午后 ~3.77 km |
| 时长与统计窗 | 模拟 4 天 / spinup 3 天 → 窗 = 第 4 天末 24 步 | 1 天 / spinup 0 → 窗 = 唯一那天的 24 步 |
| 体系性质 | **开放流动反应器**：前体物持续注入 + 背景 O₃ 渗入 → 第 4 天已达准周期稳态 | **封闭瓶**：只有 t=0 那一团空气（NOx 4 ppb、HCHO 6 ppb、O₃ 12 ppb，另有事件情景的 CO 500 / CH₃CHO 20 / ISOP 20 ppb），当天只出不进 |
| 轴的物理量 | 排放**速率**倍数（×7.25） | 初始**浓度**倍数（×11） |

**为什么封闭瓶仍然更低（三条，按量级排序）**

1. **NOx 一次性消耗**：S1 每个响应日注入约 24 ppb 的 NOx（1 ppb/h × 24 h），S2 全天总库存只有
   4 ppb，且不断被 OH 终清除成 HNO₃/PAN——**没有补充**。O₃ 是 NOx 的催化产物，NOx 的
   "周转次数"被自由基寿命限死。
2. **背景 O₃ 的渗入**：S1 有 25.2 ppb 的外部背景以 τ = 1 天持续拉入，S2 完全切断。
3. **起点与时间常数**：S2 从"晨间"起步（O₃ 12 ppb、NOx 4 ppb），一天的时间常数短于自由基与
   NOx 的清除时间，来不及爬到高位；S1 已经跑完 3 天建立期。

**两个数字各自在回答什么（对外引用时的口径）**

- **281.4 ppb** = "把这座城市的 NOx 与 NMVOC 排放同时放大 7.25 倍，第 4 天午后 1 小时 O₃ 能冲到
  多少"——**排放政策空间**（准稳态、有源有汇）。
- **119.2 ppb** = "一团**事件级**晨间空气（NOx 4 ppb／HCHO 6 ppb，CO 500、CH₃CHO 与 ISOP 各
  20 ppb，静稳无通风、混合层深 3.8 km）在完全封闭、无排放、无外部输入下，自身一天最多把
  O₃ 抬到多少"——**观测状态空间的敏感性**（瞬态、封闭）。

**两者不可比大小，只可比形态**（屋脊位置、体制分界、梯度方向）。两图轴序已统一为
**VOC 在横轴、NOx 在纵轴**，所以形态可直接并排对读；但仍只能比无量纲量：
FNR / VNr 等值线位置、屋脊在各自轴百分比上的形状、峰值率 PR。同一坐标点也不能跨图读——
S1 的 (0,0) 是"参考清单"，S2 的 (0,0) 是"那团空气的初值"，两者毫无共同之处。

> **上限由什么决定（实测）**：末 24 步窗均收支里化学消耗占 88%（含 NO+O₃ 滴定——它随 O₃
> 升高自动增大，是**结果**不是**原因**），物理汇（沉降 + 稀释 + 边界层膨胀）只占 12%，但正是
> 这 12% 决定稳态位置。本目录默认已把其中两项设成事件情景（`DIL = 0`、BLH ×2，实测 +86%），
> 再叠加富反应性碳（ISOP/CH₃CHO 5→20、CO 200→500 ppb，+34%）→ 峰值 42.8 → **119.2 ppb**。
> 三条**无效**的直觉：光照已顶格（JFAC=1.0 晴空 + 夏至 + LAT 40°，改低纬反而 −1.2%）；
> **初始 O₃ 不是上限旋钮**（起点 30–60 ppb 会被拉回当日准稳态 ~41 ppb，且 t=0 不进统计窗，
> 它只改 O₃mean 9.9→16 ppb）；**延长时长是反向的**（2 天 → −97.5%，封闭箱被沉降洗空）。

## 7. 运行注意事项

- 极端角点（低 NO₂ × 高 HCHO 等）可能偶发 CVODE 夜间收敛失败——已实现**韧性重启**
  （同容差 + MAX_STEP 10 s）+ **统计保留降级**（已积累 >0 统计的案例警告而非 fatal，
  行仍写出，run_meta/stdout 记录）；零统计失败才计入 `nfail`。
- **`nfail>0` → 执行器末尾 `[FATAL]` + `stop 1`，`run.py` 一律拒出表**（与 EKMA-emission 同口径；
  截断网格不会被静默规整成"完整"响应面）。降级（`nwarn`）不致命，但会写进 `run_meta.txt`。
- **氮守恒不变量已升级为致命门**：每个案例的首个 RHS 快照检一次，违例计数 `M` 进收尾行
  `FAC-BATCH-OMP ok (… 失败 N, 降级 K, 氮守恒违例 M, …)`；`M>0` 与 `N>0` 同样 `stop 1`，
  本 `run.py` 因此**拒出表**。当前归档实测 **M=0**（121 例）。违例原始行
  `[氮守恒 FAIL] |Σ w_N·dy/dt| = …` 也被捕获，汇总成 `run_meta.txt` 的 `氮守恒不变量:` 一行。
- `run_meta.txt` 另含 `侧车单位` 一行：`_ipr.csv/_irr.csv` 是末 24 步**平均倾向** [molec/cm³/s]，
  与 `fac_v3`（CASE 用的三件套）的**窗口积分** [molec/cm³] 相差 dtout 倍，跨程序对账前必须换算。
- `fac_batch_omp` 的 IC 模式第 9 参支持"物种名:物种名"（本实验用 `NOX:HCHO`）、
  `NOX`（NO+NO₂ 组）、`VOC`（活性 VOC 组）三种 spec 写法；EKMA-emission 用 `EMIS` 模式。


---

## 8. 换成你自己的真实数据：能换什么、不能换什么

> 本节回答：如何在已有 `config_override/` 中替换实验数据，**不改代码**。
> 当前默认基线已可直接运行；真实数据替换是可选操作。
> **默认网格 0–1000% / 100%（11×11=121）**；`--pmax/--step` 可改，但改动即改变轮次间可比性。

### 8.0 速览（一句话版）

| | 项 | 改哪里 |
|---|---|---|
| ✅ | **IC 基线 15 物种**（NO₂/HCHO 两轴的 0% 点 + NO/CO/CH₄ 等固定物种） | `config_override/initialConcentrations.config` |
| ✅ | 气象 6 序列 / 光解 LAT·DOY / 干湿沉降 / 稀释 / 步长与容差 | `config_override/` 及其 `constraints/`（唯一输入目录） |
| ❌ | 机制（封版）· 轴定义（ic1=NOX 组 / ic2=HCHO）· O₃ 不缩放 · 排放关闭 · 1天/0spinup · 生成物（runtime/CSV/图）· bc 张弛（无效） | —（实验定义/封版，见 §8.3） |

### 8.1 替换机制总览（三层，各管一类数据）

| 层 | 改哪里 | 管什么 | 生效方式 |
|---|---|---|---|
| **A. 初值基线文件** | 已有 `config_override/initialConcentrations.config`（格式同 configs 模板，单位可 ppb） | **全部初始浓度基线**（两轴的 0% 点 + 全部固定物种） | 本实验唯一 IC 来源；**缺失即报错退出**（`run.py` 自动归档来源） |
| **B. 其余本地配置与序列** | 编辑已有 `config_override/` 文件及 `constraints/` 数据 | **气象 / 光解地点 / 沉降 / 稀释 / 运行参数**等 | 整目录装配进临时目录后施加 FNR 规则（排放关、1 天/0 spinup），再交 cfgbuild（见 §3） |
| **C. CLI `--pmax/--step`** | 命令行 | 扰动幅度/步长 | 默认 1000/100；换网格即可变性由你自己掌控 |

`run.py` 内置的 `IC_BASE` 常量**已移除**：初值只可能来自 `config_override/`，
不存在任何内置缺省初值。

完整目录已经建立（见 §2），无需新建。若选择导入真实站点文件，相关落位示例如下
（仅列相关项，不表示删掉其他已有文件）：

```
EKMA-FNR/config_override/
├── initialConcentrations.config          # 编辑已有初值基线
├── meteo.config                          # 指向你的实测序列
├── photolysis.config                     # LAT/DOY 改为观测地点/日期
├── model.parameters                      # 可调整步长/容差
├── drydep.config / wetdep.config / dilution.config
└── constraints/
    ├── met_TEMP_real.txt                 # 你的实测逐时序列（t[s] K）
    ├── met_PRES_real.txt / met_RH_real.txt
    ├── met_BLH_real.txt / met_RAIN_real.txt / met_JFAC_real.txt
    └── ...（文件名自定，meteo.config 里对应即可）
```

⚠ **口径**：复制 `config_override/` 之后施加 FNR 规则，三类文件的行为各不相同：
- 初值**取输入目录中的 IC 文件**（必须存在，没有内置缺省值）；
- 排放始终强制关闭，`emissions.config` 中的排放值不会生效；
- `model.parameters` **只强制天数 1 和 spinup 0**，其余步长、RTOL/ATOL 等仍可在本地配置。

### 8.2 ✅ 可替换（真实数据域）——逐项

| 数据 | 你的真实数据 | 在哪改 | 格式红线 | 替换后的效果 |
|---|---|---|---|---|
| **IC 基线：轴物种**（NO、NO₂、HCHO） | 观测浓度（如实测 HCHO 6 ppb、NO₂ 3 ppb、NO 1 ppb） | `config_override/initialConcentrations.config` | 单位 cm3/ppb/ppm；物种 ∈ 机制 VARIABLE | **0% 点 = 你写的值**（默认 +1000% = ×11）；整个平面随之平移缩放 |
| **IC 基线：固定物种**（NO/CO/CH₄/CH₃CHO/ISOP/MVK/MACR/MGLY/SO₂/DMS/H₂O₂/PAN） | 同上（观测或清单 morning 值） | 同上（同一文件） | 同上 | 改变化学环境（如 CO/CH₄ 影响 OH 收支）；扰动轴本身不变 |
| **气温/气压/湿度** | 站点实测逐时序列 | `config_override/constraints/met_{TEMP,PRES,RH}_*.txt` + `config_override/meteo.config` 指向 | **TEMP 单位 K**（踩过 °C 的坑：会 t=0 崩溃）；t[s] 严格递增；≤128 点 | 温度影响速率常数/PAN 分解；RH 影响 OH·HO₂ 分配 |
| **边界层高度** | 实测 BLH 逐时序列 | `config_override/constraints/met_BLH_*.txt` | 同上；**不得有非正值**（除零） | 影响稀释账④与干沉速率 |
| **降水** | 实测 RAIN | `config_override/constraints/met_RAIN_*.txt` | mm/h | 湿沉账（晴天=0） |
| **光解缩放 JFAC** | 实测/辐射反演 | `config_override/constraints/met_JFAC_*.txt` | 无量纲 | 直接乘入全部 J |
| **光解地点/季节** | 站点纬度、观测日期 | `config_override/photolysis.config` 的 `LAT`/`DOY` | 度 / 1–365 | SZA 与日长 |
| **干沉降 vd** | 实测或参数化 | `config_override/drydep.config` | cm/s + on/off | 沉降汇强度 |
| **湿沉降 a,b** | 冲刷系数 | `config_override/wetdep.config` | Λ=a·R^b | 湿沉汇强度 |
| **稀释 DIL** | 实测通风率 | `config_override/dilution.config` | s⁻¹（**当前默认 0 = 密闭箱**；城市典型 1.16e-5 ≈ 24 h，静稳 2–4e-6） | 水平交换集总；对 O₃ 上限最敏感的物理汇 |

### 8.3 ❌ 不能替换（实验定义 / 封版约束）——逐项与原因

| 项 | 为什么不能换 |
|---|---|
| **机制 `mechanism_3_full_chain.fac`** | V2.5+ 封版冻结。换机制 = 开 V2.6（全部 EKMA 结果失去可比基准；需走评估流程） |
| **扰动轴物种：ic1=NO₂、ic2=HCHO** | 这就是"S2 = FNR 平面"的定义（FNR=HCHO/NO₂）。`fac_batch_omp` 第 9 参技术上支持任意"物种名:物种名"（如 `NO2:CO`），但那是**另一张平面**，不是本实验 |
| **NO 与 NO₂ 同因子缩放（ic1 = NOX 组）** | 纵轴语义是 **NOx 总量**；NO₂:NO 比例由 IC 基线决定（当前 3:1），改比例只改 IC 文件、不动代码 |
| **O₃ 不缩放** | 因变量原则——扰动因变量 = 自我实现的"响应面" |
| **排放关闭** | 你的决定②——开排放则 IC 在 5 天内被洗脱，S2 退化为 S1 翻版 |
| **1 天 + spinup=0** | 同上（纯初值化学的时间口径） |
| **默认幅度 1000% / 步长 100%** | 屋脊完整落在窗内且降级率最低（0.8%）；250/25、500/20、300/10 均已实测对照，见 §8.6 |
| **`runtime_fnr.dat`、CSV、图件** | 生成物——run.py 每次完整运行重建，手改会被覆盖；`runtime_fnr.dat` 是本轮**有效输入的只读记录**（真实数据请走 A/B 层） |
| **bc.config 的张弛约束** | **替换无效**：IC 模式下 Fortran 强制 `tau_mix=1e30`（关闭 BC），这是"纯初值化学"口径的一部分——观测数据请改走 IC 基线或未来扩展 |

### 8.4 真实数据格式红线（cfgbuild 编译期强制，违反即拒绝）

1. 序列：两列 `t[s] 值`、时间**秒**、严格递增、单序列 **≤128 点**（超限先抽稀）
2. 气象序列**禁止** `# unit:` 指令行（浓度换算因子会把 K 乘成 10¹²——真实踩坑）
3. TEMP 必须开尔文；RH 0–100；BLH > 0（非正值拒绝）
4. IC 物种必须 ∈ 机制 `VARIABLE` 段（大小写敏感）；单位仅 cm3/ppb/ppm
5. 序列池 `NSER_MAX=40` 条（气象 6 条已占 → BC `file:` 最多还能登记 34 条；且见 8.3——**BC 在 IC 模式无效**）

### 8.5 替换后的验证清单（建议顺序）

以下均在 `EKMA-FNR/` 目录执行；冒烟也会覆盖输出，先另存需要保留的 `output/`。

1. `python run.py --pmax 40 --step 10`（5×5 冒烟，~20 s）→ 看 O₃ 量级合理（城市环境几十 ppb）
2. 查 `output/ekma_fnr_grid_raw.csv` 的 `nst` 列：无异常小值（小值 = 收敛失败已降级的案例）
3. 检查响应方向物理合理：低 HCHO 列的 O₃–NO₂ 曲线应有滴定下降支；高 HCHO 列应单调升
4. `output/run_meta.txt` 核对机制 md5、`输入目录` 与参数行；确认 `[配置] 输入目录` 与 `[IC 基线]` 指向 `config_override/`（见 §9.5）
5. 通过后再跑全网格（默认 11×11；`python run.py`），再单独执行 `python plot.py`


---

## 8.6 决议与状态记录（2026-10-28）

### 三项口径决议（用户确认，本实验的"宪法"）

| # | 决议 | 落实 |
|---|---|---|
| ① | **纵轴 = 初始 NOx 组**（NO 与 NO₂ 同因子缩放，保持 IC 里的 NO₂:NO 比例） | `fac_batch_omp` IC 掩膜 `ic1 → NOX 组`（`a9="NOX:HCHO"`） |
| ② | **排放关闭**，纯初值化学 | 装配时清空 `emissions.config`（`runtime_fnr.dat` 的 `EMIS` 全 0）+ Fortran `esc=0` + BC 张弛 `tau_mix=1e30` |
| ③ | **试验基线**（NO 1 / NO₂ 3 / HCHO 6 ppb；O₃ 12 ppb 不缩放） | 唯一来源 `config_override/initialConcentrations.config`；每轮实际生效值存档于 `output/fnr_ic_base.csv` 与 `runtime_fnr.dat` 的 `IC` 行 |

### 运行状态（2026-10-29）

| 项 | 状态 |
|---|---|
| 归档（B+C 事件情景） | 11×11=121 全部完成，**失败 0 / 降级 0 / 氮守恒违例 0**（旧情景下的 1 例窗末降级消失）；墙钟 64.7 s（早期记录的 ~28 s 未复现，耗时随负载波动）；O₃max **34.14–119.21 ppb**（跨度 85.1）；基线点 (0,0) = 34.14 ppb |
| 屋脊 | NO₂ **+500%**（HCHO+0%）→ **+1000%**（HCHO+1000%），单调上移；**末列触顶**（峰落在采样上界，真实峰可能更高）；FNR_mean 覆盖 <1 仅 **1 格**／1–2 有 37 格／>2 有 83 格 → 富反应性碳把体系整体推向 NOx-limited 侧 |
| 网格选择依据 | 旧情景（DIL=1.16e-5、BLH 原值、碳 5/5/200）下 1000%/100% 屋脊 +200→+500%、无触顶；**换 B+C 后屋脊上移至 +500→+1000%，末列触顶** → 要看全峰需 `--pmax 1500`（17×17=289 例）。其他实测对照：500%/20% → 降级 3.1%；300%/10%（低 NO₂ 基线）→ 降级 6.6% 且**窗内无屋脊**；250%/25%（旧基线）→ 降级 2.5% |
| 情景假设（B+C） | `DIL=0` + BLH×2（实测 +86%）与 ISOP/CH₃CHO 5→20、CO 200→500 ppb（+34%）叠加 → 峰值 42.8 → **119.2 ppb**。代价：空气团不再是"普通晨间"，而是**静稳 + 深混合 + 富反应性碳**的事件情景；跨轮次比较必须连 `config_override/` 一起归档 |
| 轴语义依据 | 纵轴取 NOx 组（NO 与 NO₂ 同因子，NO₂:NO 恒 3:1）后，NOx 总量真正随轴变化；"只缩 NO₂ + NO 固定"在低 NO₂ 基线下会让整张面钉在 NOx-limited、窗内不出峰 |
| 输入与产物 | `run.py` 只读 `config_override/`（唯一配置输入，无内置回退）；配置在系统临时目录装配后交 cfgbuild，编译完即清理，`output/` 不留配置快照；`run.py` 是**唯一 CSV 导出者**，`plot.py` 只出两张 PNG |
| `--from` 复用语义 | 读取指定 raw（绝对路径），不装配、不编译、不跑模拟；元数据记录来源与"未计时"；指向本目录 raw 时不自拷贝 IC 存档 |
| 图件文案 | 轴/标题的基线 ppb 与 NO₂:NO 比例从 `fnr_ic_base.csv` 动态读取；图内文案一律英文（matplotlib 默认字体无 CJK，中文会成豆腐块） |
| git | 本目录与共享层改动**均未提交**，等待授权 |

---

## 9. 真实数据实验完整操作手册（手把手，未来你自己的数据）

> **场景**：你拿到了某城市某日的**观测初值**（清晨 06:00 站点浓度）与**实测气象**，
> 要回答"如果那天早晨的 HCHO/NO₂ 初始状态不同，O₃ 会怎样"——即用真实数据重做本实验。
> 全程**不改任何一行代码**、不碰机制（封版）、不碰共享层与其它程序。
> 共享内核按 §2 预先构建后，日常不再 make；没有真实数据时可直接用已有基线运行。

**本节所有命令均在 `EKMA-FNR/` 目录执行，`config_override/` 和 `output/` 都指本目录下的路径：**

```bash
cd /Users/wangling/CodeProject/LearnFortran_Projects/ozone_lab/EKMA-FNR
```

### 9.0 总流程

```
① 备数据（初值 1 个文件 + 气象 6 个文件；使用现有基线可跳过替换）
      ↓
② 编辑已有 config_override/ 及 constraints/  ← 持久用户数据入口
      ↓
③ 冒烟 5×5（--pmax 40 --step 10，~20 s）→ 判读
      ↓
④ 全网格（默认 11×11 = 121 案例）→ ⑤ 出图 → ⑥ 审读
```

全程只有两个入口文件：`run.py`（跑）与 `plot.py`（图）。冒烟及完整运行都会重新装配配置并
更新 CSV，绘图需另跑；每次重跑前先另存需保留的整个 `output/`。

### 9.1 Step 1 —— 准备初始浓度基线文件（你的"0% 点"）

**需要什么数据**：观测morning（约 06:00）浓度的 NO₂、HCHO（两轴）+ NO、O₃、CO、CH₄、
CH₃CHO、ISOP、MVK、MACR、MGLY、SO₂、H₂O₂、PAN（固定环境物种）。

**文件格式**（同 `config_override/initialConcentrations.config`）：
```
物种    值    单位(cm3|ppb|ppm)
```

**单位速查**（ppb → molec/cm³ 乘 2.46e10；直接写 ppb 让程序换算）：

| 物种 | 角色 | 当前试验基线 | 你的实测 |
|---|---|---|---|
| NO₂ | **ic1 扰动轴（NOX 组）** | 3 ppb | ← 填你的值 |
| HCHO | **ic2 扰动轴** | 6 ppb | ← 填你的值 |
| NO | **随 ic1 同因子缩放**（保持 NO₂:NO） | 1 ppb | ← 填你的值 |
| O₃ | 不缩放（因变量起点） | 12 ppb | ← **必须给**（给 0 则光化学起点失真） |
| CO | 固定 | 200 ppb | |
| CH₄ | 固定 | 1.8 ppm | |
| CH₃CHO / ISOP / MVK / MACR / MGLY | 固定 | 5/5/2/2/2 ppb | |
| H₂O₂ / PAN / SO₂ / DMS | 固定 | 0.04/0.2/0.4/0.02 ppb | |

**⚠ 三条口径红线**（写文件前重读）：
1. `NO2` 行 = ic1 轴的 **0% 参考点**（默认 +1000% = 你的值 ×11）；
2. `HCHO` 行 = ic2 轴的 **0% 参考点**；
3. `NO` 行与 `NO2` 行**同因子缩放**（决议①：ic1 = NOX 组），两者比值 = 你写的基线比
   （当前 NO₂ 3 / NO 1 = 3:1）；改比值只改 IC 文件。

### 9.2 Step 2 —— 准备气象序列（你的"环境日"）

**需要 6 个逐时序列文件**（观测日 0→24 h，121 点 hourly）：

| 文件 | 单位 | 来源 |
|---|---|---|
| `met_TEMP_real.txt` | **K**（⚠ 摄氏度 +273.15——°C 会让求解器 t=0 崩溃） | 气象站 |
| `met_PRES_real.txt` | hPa | 气象站 |
| `met_RH_real.txt` | %（0–100） | 气象站 |
| `met_BLH_real.txt` | m（>0，夜间低值） | 探空/再分析 |
| `met_RAIN_real.txt` | mm/h | 气象站 |
| `met_JFAC_real.txt` | 无量纲 0–1（晴天=1；有云按云量衰减） | 辐射/云观测 |

**格式**（每个文件，两列、无 `# unit:` 行）：
```
# 说明注释行（可有）
t[s]  值
0     295.15
3600  294.9
...
86400 293.2
```
红线：时间**秒**、严格递增、单文件 **≤128 点**（121 点 hourly 合规；5 分钟数据需先抽稀）。

### 9.3 Step 3 ——（可选）光解地点/日期与汇

| 想改什么 | 覆盖文件 | 改哪一行 |
|---|---|---|
| 观测地点/日期 | `photolysis.config` | `LAT = 你的纬度`、`DOY = 观测日序数` |
| 干沉降 | `drydep.config` | `物种 vd on` |
| 湿沉降 | `wetdep.config` | `物种 a b` |
| 通风稀释 | `dilution.config` | `DIL = s⁻¹` |

**不要覆盖** `bc.config`（IC 模式下张弛被强制关闭，写了也不生效——§8.3）。

### 9.4 Step 4 —— 编辑已有数据（或可选导入）

`config_override/` 及 `constraints/` **已存在，不需 mkdir**。通常直接编辑
`config_override/initialConcentrations.config`、气象序列及相关配置即可；不要编辑共享层，
也不必去动临时装配目录（每次完整运行都会重建并在编译后清理）。保留 `constraints/` 子目录，`meteo.config` 中的引用应形如
`constraints/met_TEMP_real.txt`，不要把序列摊平到顶层。

若数据在别处，下面仅为**可选导入示例**（当前目录仍是 `EKMA-FNR/`；将 `/path/to/` 换成真实绝对路径）：

```bash
# ① 初值基线（替换已有文件，文件名必须精确）
cp /path/to/initialConcentrations.config config_override/initialConcentrations.config

# ② 气象序列（文件名自定，保留 constraints/ 相对路径）
cp /path/to/met_TEMP_real.txt config_override/constraints/
cp /path/to/met_PRES_real.txt config_override/constraints/
cp /path/to/met_RH_real.txt   config_override/constraints/
cp /path/to/met_BLH_real.txt  config_override/constraints/
cp /path/to/met_RAIN_real.txt config_override/constraints/
cp /path/to/met_JFAC_real.txt config_override/constraints/

# ③ meteo.config（引用 constraints/ 下对应的 6 个文件）
cp /path/to/meteo.config config_override/meteo.config

# ④（按需选择）光解地点/沉降/稀释
cp /path/to/photolysis.config config_override/photolysis.config
cp /path/to/drydep.config config_override/drydep.config
cp /path/to/wetdep.config config_override/wetdep.config
cp /path/to/dilution.config config_override/dilution.config
```

### 9.5 Step 5 —— 冒烟验证（5×5，~20 s）

```bash
python run.py --pmax 40 --step 10          # 当前目录：EKMA-FNR/
```

**预期输出（逐行判读）**：
```
  [配置] 输入目录: <仓库绝对路径>/EKMA-FNR/config_override     ← 唯一配置输入
  [IC 基线] config_override/initialConcentrations.config       ← IC 取自本地文件（无缺省回退）
═══ EKMA-FNR: 5×5 = 25 案例  (NO₂/HCHO 初值扰动 0–40% / 10% ⇒ 因子 1.00–1.40, 线性网格, 排放关闭, 1 天) ═══
  FAC-BATCH-OMP ok (25 案例, 失败 0, 降级 0, xx.xx CPU s, 线程 4)
  网格: 5×5 = 25 案例 (NO₂ 0–40%/10%, HCHO 0–40%/10%)
  产物: EKMA-FNR/output/ekma_fnr_grid.csv
  O3max 范围 …
```

输入目录当前共 **28 个文件（11 顶层 + 17 constraints）**；增删文件不会改变运行方式，
但 `run.py` 只要求 `initialConcentrations.config` 与 `model.parameters` 必须存在，
其余缺哪些由 cfgbuild 在编译期报错。执行器内部诊断（如物种匹配警告）只在出现
`FAC-BATCH`/`FATAL`/`警告` 时回显，需要更多细节时直接手动运行 `facsimile_impl/fac_batch_omp`。
`[IC 基线]` 行只说明采用了本地 IC 文件，**不证明数据来自观测**——即使文件仍是复制来的默认基线也会这样标记。

**判读清单**（任何一条不过 → 先修数据再继续）：
1. 核对 `[配置] 输入目录` 路径，并用 `grep -A1 "^IC" output/runtime_fnr.dat`、`grep "^MODEL" output/runtime_fnr.dat` 核对**实际生效**的初值与时长；
2. 核对 `[IC 基线]` 行指向 `config_override/initialConcentrations.config`（唯一 IC 输入，无缺省回退，也不使用根 IC）；
3. `output/fnr_ic_base.csv` 的来源注释及数值与本地 IC、`runtime_fnr.dat` 的 `IC` 行一致；ppb 值乘 2.46e10 得 `value_cm3`（ppm 乘 2.46e13，cm3 直接使用）；
4. O₃ 量级合理（城市清晨初值 + 1 天 → 20–80 ppb 区间）；
5. `nst` 列无异常小值（有 = 个别案例夜间收敛失败已降级，见 §7）。

### 9.6 Step 6 —— 全网格运行 + 出图

在 `EKMA-FNR/` 目录运行；需要保留冒烟结果时先另存 `output/`：

```bash
python run.py                  # 默认 11×11 = 121 案例（1 天案例），自动 cfgbuild，不重编内核
python plot.py                 # → output/ 仅两张 PNG（屋脊表与派生表已由 run.py 写出）
```

### 9.7 Step 7 —— 结果审读（物理检查表）

- 屋脊是否**单调上移**（HCHO 越多能容忍越多 NO₂）？——应该；
- 低 HCHO 列是否有**滴定下降支**（O₃ 随 NO₂ 先升后降）？——应该；
- FNR=1/2/4 等值线是否落在物理合理的 HCHO/NO₂ 区间？——对照你的实测 FNR；
- O₃max 绝对量级 vs 观测当日 O₃：纯初值化学、无排放，O₃ 只由初值化学演化——
  **不承诺复现观测 O₃ 绝对值**（那是排放驱动的 S1 的职责），本实验回答的是
  "给定初始状态，O₃ 对 HCHO/NO₂ 初值的敏感性"。

### 9.8 常见错误速查（真实数据特有）

| 症状 | 根因 | 处置 |
|---|---|---|
| CVODE 在 `t=0` 崩溃 | TEMP 写成了 **°C**（26 K = −247 °C） | +273.15（本会话踩过的真实坑） |
| `序列 N 点 > 128` | 分钟级观测未抽稀 | 降采样到 hourly |
| `序列时间非单调` | 时间戳含重复/乱序 | 排序去重 |
| `IC 基线未知单位 'ug/m3'` | IC 只认 cm3/ppb/ppm | 先自行换算成 ppb（×2.46e10 得 cm3） |
| `缺少完整输入目录` | `config_override/` 被移动或改名 | 该目录是本实验唯一配置输入，恢复路径即可 |
| `缺必要输入文件: initialConcentrations.config …` | IC / model.parameters 不在 `config_override/` | 补齐同名文件（无内置缺省值） |
| IC 改了但 S2 几乎不动 | 基线里 NO 底数过大或 HCHO 过低 | 用 §9.1 表的基线量级（NO 1 / NO₂ 3 / HCHO 6 ppb） |
| 个别案例 `[警告] ... 收敛失败` | 极端角点（已保留统计） | 看 `nst` 列；网格太极端时收窄 `--pmax` |
| BLH 序列被拒 | 含非正值 | 白天/夜间都给正值（夜间 ≥50 m） |
| 序列池超限 | 气象 6 条 + BC `file:` > 40（`NSER_MAX`） | IC 模式下 BC 本就无效 → 不加 BC 序列 |
| 想换轴物种 | 那是另一张平面 | `fac_batch_omp` 第 9 参 `NO2:CO` 等；但不再是 FNR 实验 |

### 9.9 与"改代码"路线的对比

| 做法 | 适用 | 代价 |
|---|---|---|
| **本手册（编辑 config_override/）** | 换数据（初值/气象/汇/地点） | **零代码改动**；override 目录即实验记录的一部分 |
| 改 `fac_batch_omp` / `run.py` | 换轴物种、改实验口径或新增维度 | 那是新实验（进 V2.6 评估流程）；数据本身不需要改代码 |

---

## 10. 最小化运行 vs AtChem2，与"选择性输入"的精确语义

| | AtChem2 | ozone_lab |
|---|---|---|
| 机制 | MCM Fortran 源码（先编译） | `.fac` → mechc.py → gen_v3（`make all`） |
| 初始浓度 | `initialConcentrations.config` | `config_override/initialConcentrations.config`（未列=0；**本实验的基线就是它**） |
| 气象/光解 | `environmentVariables.config` + 光解内置 | `config_override/` 的 `meteo.config` + `photolysis.config` + `photolysis-rates_MCM`（光解显式必需） |
| 运行控制 | `model.parameters` | `config_override/model.parameters`（天数/spinup 由实验强制，其余取本地值） |

概念三要素一致（机制 + 初始场 + 气象/光解）。文件层面：ozone_lab 的光解是显式必需输入；
`cfgbuild` 打开 10 份配置，其中 6 份可全空。**选择性输入语义 = 未列出 = 零/无**：
排放不输入=纯初值化学（本实验正是如此）；BC 不输入=无约束自由运行——
不存在"内建经验值顶上"，经验基线在随仓库分发的 configs 数据文件里。

## 11. 站点在 S2 响应面上的落点（清单/观测只算坐标，不必灌进模型）

> **核心结论（与 EKMA-emission §9.3 完全对称）**：**初始浓度基线不必改动，观测值也不需要作为
> 模型输入**——S5P/站点实测的 HCHO、NO₂ 只用来算两个百分比，在图上定位。
> 这与 S1 用"清单速率 / 参考清单"定位是同一个逻辑：
>
> | | S1（EKMA-emission） | S2（本实验） |
> |---|---|---|
> | 图的基线 | 参考排放清单 | 参考晨间初值（`fnr_ic_base.csv`） |
> | 你的数据 | NOx/NMVOC **清单速率** | HCHO/NO₂ **观测浓度** |
> | 用途 | 只算坐标 | 只算坐标 |
> | 近似来源 | 你的源谱 ≠ 参考源谱 | 你的**非轴物种** ≠ 参考基线（见 §11.3） |
>
> **S2 比 S1 更直接**：S2 的两根轴**就是**被测的量本身，换算只是"除以基线减一"。

### 11.0 ⚠ 先看网格下界：本图只覆盖 0% → +1000%（**不含低于基线的值**）

`run.py` 写死 `fmin = 1.0`（= 0%）——**没有负扰动**。因此：

| 你的观测 | 与基线之比 | 落点 | 结果 |
|---|---|---|---|
| HCHO 8.1 ppb（基线 6） | 1.35 | **+35%** | ✅ 在图上 |
| HCHO 4.5 ppb（基线 6） | 0.75 | **−25%** | ❌ **落在 0% 下边界外，图上没有** |

当前基线是 NO 1 / NO₂ 3 / HCHO 6 ppb；实际站点若低于此（如清洁郊区 HCHO 3–5 ppb），
模式 ① 无法落点 → 必须用模式 ②（把基线换成你的站点值，站点即原点）。详见 §11.4。

### 11.1 参考基线逐物种明细（落点分母）

来源 `output/fnr_ic_base.csv`（= `config_override/initialConcentrations.config`）；`ppb = value_cm3 / 2.46e10`：

| 物种 | 中文名 | 基线值 [cm³] | 基线 [ppb] | 缩放归属 |
|---|---|---|---|---|
| **HCHO** | **甲醛** | 1.48e11 | **6.0** | **ic2 扰动轴（y）** |
| **NO2** | **二氧化氮** | 7.38e10 | **3.0** | **ic1 扰动轴（x，NOX 组）** |
| NO | 一氧化氮 | 2.46e10 | **1.0** | **随 ic1 同因子缩放**（保持 NO₂:NO = 3:1） |
| O3 | 臭氧 | 3.00e11 | 12.2 | 固定（因变量） |
| CO | 一氧化碳 | 1.23e13 | **500** | 固定（事件情景抬高） |
| CH4 | 甲烷 | 4.50e13 | 1829（1.83 ppm） | 固定 |
| CH3CHO | 乙醛 | 4.92e11 | **20.0** | 固定（事件情景抬高） |
| ISOP | 异戊二烯 | 4.92e11 | **20.0** | 固定（事件情景抬高） |
| MVK | 甲基乙烯基酮 | 4.92e10 | 2.0 | 固定 |
| MACR | 甲基丙烯醛 | 4.92e10 | 2.0 | 固定 |
| MGLY | 甲基乙二醛 | 4.92e10 | 2.0 | 固定 |
| PAN | 过氧乙酰硝酸酯 | 5.00e9 | 0.203 | 固定 |
| H2O2 | 过氧化氢 | 1.00e9 | 0.041 | 固定 |
| SO2 | 二氧化硫 | 1.00e10 | 0.407 | 固定 |
| DMS | 二甲基硫 | 5.00e8 | 0.020 | 固定 |

> **原点 (0%,0%) = 上表的 HCHO 6 / NO₂ 3 / NO 1 ppb 基线，不是你的站点。**

### 11.2 落点公式与实算例子

```
pert_HCHO[%] = (HCHO_obs[ppb] / 6 − 1) × 100     # x 轴（VOC 轴）
pert_NO2[%]  = (NO2_obs [ppb] / 3 − 1) × 100     # y 轴（NO 同因子，NOx 组语义）
```

**例 1（图内，完整分步）**：清晨 HCHO **12** ppb、NO₂ **6** ppb

| 步骤 | HCHO（x） | NO₂（y） |
|---|---|---|
| ① 观测值 | 12 ppb | 6 ppb |
| ② 除以基线 | 12/6 = 2.00 | 6/3 = 2.00 |
| ③ 减 1 × 100 | **+100%** | **+100%** |
| ④ 落点 | **x = +100%** | **y = +100%** |

图上读数（`ekma_fnr_grid.csv` 实测，恰在格点上）：

| 项 | 值 | 来源 |
|---|---|---|
| 该点 O₃max 预测 | **51.27 ppb** | (100,100) 格点；FNR_mean = 3.30 |
| 同列（HCHO+100%）屋脊 NO₂ | **+500%** | `ekma_fnr_ridge.csv`（+0%→500、+100%→500、+200%→600、…、+1000%→1000） |
| 该列峰值 O₃ | 83.07 ppb | 同上 |

**判读**：点在屋脊（NO₂+500%）**下方** → **NOx-limited**（再加 NOx 仍增 O₃，减之则降；
从 51.27 → 83.07 ppb，还有 ~62% 空间）。

**例 2（低于下界，必须换模式）**：清晨 HCHO 4.5 ppb、NO₂ 3.6 ppb
→ HCHO: 4.5/6 = 0.75 → **−25%**（**图外**）；NO₂: 3.6/3 = 1.20 → +20%
→ HCHO 轴越界 → 用模式 ②（§11.4），把 `config_override/initialConcentrations.config` 的 HCHO 改成 4.5、NO₂ 改成 3.6（NO 同步给 1.2 以维持 3:1）。

**例 3（恰在原点）**：HCHO 6 / NO₂ 3 ppb → (0%, 0%)，O₃max = 34.14 ppb。

**⚠ S2 的 O₃ 量级（当前归档 34.1–119.2 ppb）仍低于 S1（15.3–281.4 ppb）**：本实验关闭排放、只跑 1 天，
读数是"这团晨间空气自身一天能生成多少 O₃"，**不是城市全天产量**，两者不可直接比较。

### 11.3 三个必须对齐的口径（否则落点失真）

| 口径 | 基线 | 你要做的 |
|---|---|---|
| **时窗** | 轴是**初始**浓度 | 观测取**日出前（05–06 时）**值作气团初值；S5P 13:30 过境值需按典型日变化折算（如 ×0.6）或在图中如实标注"基线取自 13:30 值" |
| **单位** | ppb | S5P 给的是**柱浓度**（molec/cm²）→ 必须先 `python inputproc/sat2surf.py <raw> <out> --blh-file <blh>` 转表面等效 ppb（内置两法互证 + QA 过滤 + 抽稀） |
| **非轴物种** | 见 §11.1 的 13 个固定物种 | 它们**不能**按你的观测替换（模式①下）。若你站点的 NO/NO₂ 比、CO、生物源 VOC 与基线差别大 → 等值线有系统偏差 → 用模式②复核 |

> 口径提示：纵轴是 **NOx 组**（NO 与 NO₂ 同因子，NO₂:NO 恒 3:1），所以"NO₂ 扰动 %"即"NOx 扰动 %"。
> 若想让 NO 固定而只动 NO₂，需把 `run.py` 的 `IC1_SPEC` 改回 `"NO2"`——但低 NO₂ 基线下窗内会失去屋脊。

### 11.4 两种模式：什么时候必须动基线

| 模式 | 做法 | 适用 | 代价 |
|---|---|---|---|
| **① 参考面 + 落点**（默认） | 基线不动，算两个百分比 | 观测 **≥** 基线（HCHO ≥6、NO₂ ≥3 ppb） | 非轴物种差异带来等值线偏差 |
| **② 基线 = 你的站点** | 编辑 `config_override/initialConcentrations.config` 的 HCHO/NO₂/NO 三行 → 重跑全网格（默认 121 案例） | 观测 **<** 基线（**唯一可行**）；或要求落点精确 | 只能与自己的历史比；多跑一轮（121 例墙钟 ~65 s / 4 核） |

**建议策略**：
1. 先算两个百分比。**任一为负 → 直接走模式 ②**（不改基线的图没有负半轴）；
2. 两轴都在网格内且**远离屋脊**（如与屋脊 NO₂ 差 >50%）→ 模式 ① 即可，结论稳健；
3. 两轴在内但**紧贴屋脊**、或非轴物种与基线差 >2 倍 → 模式 ② 复核（结论一致即可放心）；
4. 任一超出网格上界 → `--pmax 400 --step 10`（或更大）扩窗，或模式 ②。

### 11.5 四种落点操作方式

| 方式 | 操作 |
|---|---|
| **手算** | §11.2 两行公式 → 在 `ekma_fnr_S2.png` 上目测（坐标轴是 %，等值线可读出 O₃ 预测） |
| **CSV 标注** | 自建 `output/obs_points.csv`（列：`pert_no2_pct, pert_hcho_pct, note`），供后续叠加 |
| **✅ 自动化·单站** | `plot.py --obs-single "HCHO=8.0,NO2=4.5,O3=30"` → 红星 + 候选框 + 体制三重校验（§11.8） |
| **✅ 自动化·多站** | `plot.py --obs-more stations.csv`（列 station,HCHO,NO2,O3）→ 编号红星 + 终端汇总表 + 每站候选表（§11.8） |

> 两个自动化参数**一次调用二选一**（同给报错）；观测口径：O₃ = 日最大，HCHO/NO₂ = 日均。

### 11.6 与 S1 落点的分工

| 平面 | 站点落点需要 | 回答的问题 |
|---|---|---|
| **S2（本目录）** | 清晨 HCHO/NO₂ 观测（可测/卫星可比） | "这个观测状态出发，O₃ 会怎样？现在处于什么体制？" |
| **S1（EKMA-emission）** | 排放清单比值（浓度不能直接落） | "如果排放这样改，O₃ 会怎样？"（政策空间） |

### 11.7 额外定位法：仅用 **O₃ + NO₂ 观测** 反查（HCHO 缺测时的替代）

> **动机**：§11.2 的方法需要 **HCHO 与 NO₂ 两个观测**。若手上只有常规监测的 **O₃ + NO₂**
> （无卫星 HCHO），仍可用本实验 `_state.csv` 侧车记录的每个案例 **(NO₂ 均值, O₃max)**
> 做查表反演，反推出 (pert_NO2, pert_HCHO) 坐标。
> **本节为额外方法，不替代 §11.2**（后者两轴本身可观测，更准）。

#### 11.7.1 可直接运行的查表脚本

```python
import sys; sys.path.insert(0, "analysis_lib")
import numpy as np
from csvio import read_table, fcol

hg, G = read_table("EKMA-FNR/output/ekma_fnr_grid.csv")
hs, S = read_table("EKMA-FNR/output/ekma_fnr_grid_raw.csv_state.csv")
ph, pn = fcol(hg, G, "pert_hcho_pct"), fcol(hg, G, "pert_no2_pct")
no2, hcho, o3 = (fcol(hs, S, "NO2mean_ppb"), fcol(hs, S, "HCHOmean_ppb"),
                 fcol(hs, S, "O3max_ppb"))

NO2_obs, O3_obs = 6.5, 29.0          # 你的观测（同窗口口径）
HCHO_obs = 13.5                       # 仅在 use_hcho=True 时使用
use_hcho = False                      # True = 有 HCHO 观测时加入（更准）

cols = np.c_[no2, hcho, o3] if use_hcho else np.c_[no2, o3]
M, SD = cols.mean(0), cols.std(0)
P = (cols - M) / SD
obs = [NO2_obs, HCHO_obs, O3_obs] if use_hcho else [NO2_obs, O3_obs]
q = (np.array(obs) - M) / SD

d = np.linalg.norm(P - q, axis=1)
idx = np.argsort(d)[:10]
print(f"最近距离 {d[idx[0]]:.3f}σ")
print(f"NO₂ 轴候选: {pn[idx].min():.0f}–{pn[idx].max():.0f}%   HCHO 轴候选: "
      f"{ph[idx].min():.0f}–{ph[idx].max():.0f}%")
```

**实测输出**（当前归档 121 例，`NO₂=8.0 ppb / O₃=90 ppb`）：

```text
最近距离 0.090σ
NO₂ 轴候选: 500–700%   HCHO 轴候选: 100–500%
```
> 旧示例 `NO₂=6.5 / O₃=29 ppb` 现在**低于模型最小 O₃（34.14 ppb）**，最近邻距离 1.43σ
> 且候选退化到 NO₂+0% 一列——属于**外推无效**，只能作为"超出覆盖范围"的反例。

候选带比 S1 宽得多（NO₂ 轴 300 个百分点）——印证 §11.7.2 的"区分能力更弱"：
S2 的 O₃ 动态范围只有 26 ppb，不足以精细分辨 HCHO 水平。**用途应限于"判断位于屋脊哪一侧"
这类定性定位。**

#### 11.7.2 实测精度

（下表格式：NO₂ 轴 / HCHO 轴；10 近邻，`seed=0`，200 次抽样）

| 检验 | 仅 O₃+NO₂ | 加 HCHO |
|---|---|---|
| 留一验证轴偏差中位 | 10% / 110% | **20% / 20%** |
| 10% 噪声下轴散布中位 | 100% / 500% | 200% / **300%** |

**解读**：加 HCHO 把 HCHO 轴的歧义显著压下来（留一 110%→20%、噪声散布 500%→300%），代价是 NO₂ 轴
略变差（10%→20%）。S2 的 O₃ 动态范围只有 34.1–119.2 ppb（1 天纯初值化学），比 S1
（15.3–281.4 ppb）窄得多，所以**反演区分能力弱**：O₃ 差 1 ppb 就可能跨好几个 HCHO 水平。
仍可用作**粗定位**（判断落在屋脊哪一侧），不宜作为定量坐标。

#### 11.7.3 三个硬约束

| 约束 | 说明 |
|---|---|
| **统计口径** | 模型端 `O3max_ppb` = 1 天窗口最大、`NO2mean_ppb` = 窗口均值；观测端用对应口径（O₃ 日最大、NO₂ 日均） |
| **多解性** | 同一 (NO₂, O₃) 常对应多个 HCHO 水平 → **必须报告候选范围**，不报单点 |
| **覆盖范围** | 当前归档模型 NO₂ 均值 0.98–13.68 ppb、O₃max 34.14–119.21 ppb；观测超出即外推，无意义 |

#### 11.7.4 四种定位法对比（本实验并存）

| 方法 | 需要的数据 | 精度 | 章节 |
|---|---|---|---|
| **A. HCHO+NO₂ 观测直接落点（首选）** | HCHO（卫星/站点）+ NO₂ | 最高（轴即观测） | §11.2 |
| **B. O₃+NO₂ 反查（本节额外法）** | O₃ + NO₂（常规） | NO₂ 轴 ±10%、HCHO 轴 ±20% | §11.7 |
| **C. 改基线（模式②）** | 同上，但观测低于基线时必须 | 站点即原点 | §11.4 |
| **D. 落 S1** | 排放清单 | 回答排放政策问题 | `EKMA-emission/README` §9.3 |

### 11.8 状态匹配定位法（`--obs-single` 单站 / `--obs-more` 多站——已实现）

> **方法学确认**（用户口径，2026-10-28）：EKMA-FNR 通过**理想晴天静稳态天气**的光化学
> IC 扰动参数空间扫描得到 S2 响应面；**站点定位 = 在参数空间的各模拟情景中，找出与现实
> 监测最相似的情境，把该情境的坐标定义为站点的响应面坐标**。——理解正确，术语叫
> **状态匹配（最近邻反查）**：以"可观测状态量"为匹配空间、"控制因子坐标"为答案。

#### 两种参数（并存，一次调用二选一；混用报错）

| 参数 | 输入 | 图面（三档，与 EKMA-emission §9.4.1 同规则） |
|---|---|---|
| `--obs-single "HCHO=a,NO2=b,O3=c"` | 单站三观测（ppb） | 红星 + 8 近邻**包络框 + 叉点** + 标注框（坐标 / 体制 / `FNR_init`） |
| `--obs-more stations.csv` | CSV 列 `station,HCHO,NO2,O3`（一行一站） | **2–5 站**：红星 + 候选叉点 + `编号: 站名 (体制)`；**>5 站**：红星 + 只标编号 |

> 三档的**完整候选明细一律打印在终端**（rank/距离σ/两轴/`FNR_init`/`FNR_mean`/O₃max/体制）；
> `plot.py` 只写 PNG，不产 CSV。标记会覆盖 `ekma_fnr_S2.png`，要无标记版就重跑不带 `--obs-*`。
> **站名请用 ASCII**：图件不配 CJK 字体，中文站名会渲染成缺字方块（脚本会给出告警；中文明细仍在终端）。
> 标注位置自动避让：右半区向左标、贴底时抬高，避免压住右下角图例。

观测口径：**O₃ = 日最大**，HCHO/NO₂ = **日均**（与模拟窗口均值同口径）。
`--obs-single` 与 `--obs-more` 同时给出 → 报错退出（一次一种输入）。

#### 匹配空间与答案

| | 量 | 来源 | 观测口径 |
|---|---|---|---|
| 匹配维度 1 | `NO2mean_ppb`（窗均 NO₂） | `_state.csv` 侧车 | 日均 ppb |
| 匹配维度 2 | `HCHOmean_ppb`（窗均 HCHO） | 同上 | 日均 ppb |
| 匹配维度 3 | `O3max_ppb`（窗内 O₃ 峰） | 同上 | 日最大 ppb |
| **答案 x** | `pert_hcho_pct`（初始 HCHO 扰动，横轴） | `grid.csv` | — |
| **答案 y** | `pert_no2_pct`（初始 NO₂ 扰动，纵轴） | `grid.csv` | — |
| **体制** | `FNR_init`（初始组成比，<1 / 1–2 / >2 三档） | `grid.csv` | 与图上 FNR 等值线同族 |

#### 用法与输出

```bash
python EKMA-FNR/plot.py --obs-single "HCHO=8.0,NO2=4.5,O3=30"     # 单站
python EKMA-FNR/plot.py --obs-more stations.csv                   # 多站（列 station,HCHO,NO2,O3）
```

**图面标记三档**（统一走 `draw_station`，单站与多站同一入口；与 EKMA-emission §9.4.1 同规则）：

| 站点数 | 红星 | 8 近邻候选（终端表始终给全明细） | 文字标注 |
|---|---|---|---|
| **1 站** | ✅ s=320 + 图例 `station (matched)` / `8-NN candidates` | ✅ **包络框**（红点矩形）+ 叉点 | 左下角框：坐标 + `_regime(FNR_init)` + `FNR_init=值` |
| **2–5 站** | ✅ s=260 | ✅ 仅**叉点**（多框会互相遮挡） | 每站 `编号: 站名 (体制)` |
| **>5 站** | ✅ s=260 | ❌ 不画（避免糊图） | 只标**编号**，明细全在终端候选表 |

```text
终端: 多站汇总表 + 每站 8 近邻候选表（rank/距离σ/两轴/FNR_init/FNR_mean/O3max/体制）
      plot.py 只写 PNG，候选明细不落 CSV（全仓 CSV 由 run.py 独占）
```

⚠ **标记会画进 `ekma_fnr_S2.png` 本身**（覆盖同名文件，不是另存新图）。
存档要的是无标记响应面时，重跑 `python plot.py`（不带 `--obs-*`）即可恢复。
前置条件：`output/ekma_fnr_grid_raw.csv_state.csv` 与 `ekma_fnr_grid.csv` 存在（先跑 `run.py`）。

**实测样例**（当前 11×11 归档，3 站）：

```text
  站点           距离σ    HCHO%    NO₂%  FNR_init  体制
  北大站          0.92      200        0     6.00  NOx-limited  ⚠ 候选跨多体制
  化工园          1.69      200     1000     0.55  VOC-limited
  背景点          0.38        0        0     2.00  transition  ⚠ 候选跨多体制
  体制分布: NOx-limited×1, VOC-limited×1, transition×1
```

**多站点的科学价值**（共面对照）：多站各自独立匹配后，**体制分布表**直接给出"这一天
各站点所处体制的空间分异"——若各站落点聚在响应面同一侧且距离 σ 小 = 密闭瓶假设对这天
整体成立；若各站散落且 σ 大 = 存在输送/日变排放，密闭瓶假设对部分站点失效（本身就是诊断信息）。

#### 四条使用约束（已内置为自动检查/警告）

| 约束 | 实现 |
|---|---|
| **包络检查** | 观测超出模拟状态范围（NO₂mean 0.98–13.68 / HCHOmean 4.44–23.56 / O₃max 34.14–119.21 ppb 之外）→ 逐维 ⚠ 警告"外推不可靠" |
| **多解性** | 同一状态可对应多个 (NO₂, HCHO) 组合（EKMA 双支歧义）→ 强制输出候选**带**（8 近邻）与体制集合，跨体制必警告 |
| **口径对齐** | O₃_obs 必须是**日最大**（对 O3max）、HCHO/NO₂_obs 必须是**日均**（对窗口均值）——用峰值 NO₂ 比均值会系统性错位 |
| **封闭瓶假设** | 匹配隐含"站点当天 ≈ 从某初值出发的密闭瓶演化"。真实当天有排放/输送 → 偏差体现在匹配距离 σ 上：σ 大（>1.5）= 该实验解释不了你的观测日，慎用结论 |
