# EKMA-emission/ —— S1 排放正交扰动响应面

> **一句话**：固定机制/气象/初值，只把 **NOx 排放**与 **E_VOC 排放**按 **0%→625%（25% 步长，默认）**
> 做正交扰动，跑 26×26=676 个独立 0-D 箱（4 天模拟 / 3 天 spinup），得到 O₃max 响应面与屋脊。
>
> ⚠ **现存归档**：`output/` 的 CSV/图 = **2026-10-29 完整运行**，默认网格 **0–625%/25%（26×26=676，4 天/3 天）**，
> 失败 0 / 降级 0 / **氮守恒违例 0**，O₃max 15.26–281.41 ppb，基线点 77.12 ppb（见 `run_meta.txt`）。
> 同参数两轮的**数据行逐位相同**（仅 OpenMP 案例完成序、以及 runtime 首行的临时装配目录名不同），
> 故 `耗时` 随负载波动（实测 899 s 与 1458 s 两次），**不作为回归判据**；回归看数值与守恒计数。
> §9.3/§9.4 的落点与反查算例仍引用被本次覆盖的 41×41 旧归档数值，引用前需按新 CSV 重算。

---

## 1. 实验定义

| 项 | 内容 |
|---|---|
| 类型 | 排放扰动（emission-based EKMA） |
| 执行器 | `facsimile_impl/fac_batch_omp`（C-ABI + OpenMP 批量，一个进程跑完网格） |
| 固定条件 | ① 机制 `mechanism/mechanism_3_full_chain.fac`（V2.5+ 冻结）<br>② 本目录气象 `config_override/meteo.config`（当前为 5 天晴天逐时序列）<br>③ 本目录初值 `config_override/initialConcentrations.config` |
| 扰动轴 | x = **E_VOC**（NMVOC 排放组）、y = **E_NOx**（NO,NO₂）<br>两者均为相对本目录 `config_override/emissions.config` 基线的百分比 |
| 扰动范围 | **0%（无扰动）→ 625%（+625%）**，步长 **25%**（默认）⇒ 26 点/轴，因子 1.00–7.25 |
| 案例数 | 26 × 26 = **676**（默认） |
| 时长 | **模拟 4 天 / SPINUP 3 天**（响应统计 = 末 24 步 = 第 4 天） |
| 输出量 | O₃max（**末 24 个输出步窗口最大**，不是 8h 日最大）+ 24 列指示因子/过程量 |

每轮实际生效的输入以本目录 `output/runtime_emission.dat` 为准（cfgbuild 解析出的 IC/EMIS/MODEL/MET
数值，含被实验改写为 4 天/3 天的 `model.parameters`）；输入原文就是 `config_override/` 本身。
本程序**没有任何全局配置层**（仓库里不存在根 `configs/`），`output/` 里也不留配置快照目录。
机制与执行器路径相对仓库根，其余输入路径相对本目录。

> **E_VOC 组口径**（`modules/runtime_dat_mod.f90::build_groups`，24 物种）：
> HCHO CH3CHO ISOP DMS MVK MACR MSA SO2 ALK OLE ARO1 ARO2 KET TERP
> + V2.4 新增 CSL LALK MEOH ETOH MTBE ESTER SESQ ETH IOLE ALD2
> —— 注意该组含 SO₂/DMS/MSA（硫）三个非 VOC 物种，是本工程既有的集总约定。

---

## 2. 文件与用法

```
EKMA-emission/
├── run.py      # ① 读本目录 config_override → 临时装配 → cfgbuild → n×n 案例 → CSV
├── plot.py     # ② 读 grid.csv → **只出图**（响应面 PNG + NOx 切片 PNG）
├── README.md
├── config_override/   # 已建立的持久输入：11 个顶层文件 + 17 个 constraints 文件
│   ├── meteo.config / initialConcentrations.config / emissions.config
│   ├── bc.config / drydep.config / wetdep.config / dilution.config
│   ├── model.parameters / scenarios.config
│   ├── photolysis.config / photolysis-rates_MCM
│   └── constraints/   # met_sunny5d_*.txt、met_*_demo.txt、bc_*_demo.txt 等
└── output/                          （逐文件数据字典见 §7；结果产物，不是编辑入口）
    ├── runtime_emission.dat          ① 本轮有效输入记录（cfgbuild 解析出的 IC/EMIS/MODEL/MET 值）
    ├── ekma_emission_grid_raw.csv    ② fac_batch_omp 原始 24 列输出
    ├── ekma_emission_grid_raw.csv_ipr.csv   ②a 侧车: 111 物种 × 6 过程账 窗均 [molec/cm3/s]
    ├── ekma_emission_grid_raw.csv_irr.csv   ②b 侧车: r1..r203 逐反应 IRR 窗均 [molec/cm3/s]
    │        （⚠ 2026-10-29 已删除：旧版 batch 坐标列写错；重跑 run.py 自动重建，格式正确）
    ├── ekma_emission_grid_raw.csv_state.csv ②c 侧车: 实验 J 状态量均值对（本实验未消费，备查）
    ├── ekma_emission_grid.csv        ④ 规整主表（plot.py 唯一输入）
    ├── ekma_emission_ridge.csv       ⑤ 屋脊轨迹
    ├── ekma_emission_S1.png / ekma_emission_slices.png  ⑥ 图件
    └── run_meta.txt                  ⑦ 运行元数据
```

`config_override/` 是从历史运行输入建立、**长期保留的持久输入目录**：`run.py` 只读它，
从不改写或删除它。保留现有基线即可直接运行，换数据就直接编辑其中的配置及 `constraints/` 文件，
不必每轮创建目录——仓库里没有需要先行编辑的全局配置层。

**初次准备（与日常流程分开）**：共享内核尚未构建时，在仓库根目录 `ozone_lab/` 执行：

```bash
make all                                 # 准备共享执行器及构建产物；已构建则跳过
```

**日常流程**：共享内核已构建后，在 `EKMA-emission/` 内执行：

```bash
python run.py                            # ① 跑实验（默认 0–625% / 25%，676 案例 / 4 天，约 8–12 min / 4 核）
python plot.py                           # ② 出图

# 自定义网格（例：0–100%，10% 步长 → 11×11）
python run.py --pmax 100 --step 10
python plot.py
```

只改配置/数据**无需重新 make 或重编内核**；`run.py` 会自动调用 `cfgbuild` 重新编译配置。
完整运行会重新装配配置、生成 runtime 和**全部 CSV**（含 grid/ridge/derived；`output/` 不再有配置快照目录），
两张 PNG 需另跑 `plot.py`。
**需要保留历史 `output/` 时，先另存再运行**。本流程仍依赖完整仓库与运行环境，不能只搬走实验子目录运行（见 §6）。

**注意**：`fac_batch_omp` 的第 10 参 `LIN` 选**线性网格**（缺省为等比 `fmin×gr^i`）。
线性网格是"等扰动百分比"轴语义的必要条件——等比网格的步长不是等差。

---

## 3. 图怎么读

| 元素 | 含义 |
|---|---|
| 横轴 | NMVOC 排放扰动 %（0 = 本目录 `config_override/emissions.config` 基线；实际生效值见 `output/runtime_emission.dat` 的 `EMIS` 行） |
| 纵轴 | NOx 排放扰动 % |
| 色标 / 黑实线 | O₃max [ppb] 及其等值线 |
| 白虚线 | **屋脊**：每个 NMVOC 扰动列上 O₃max 达到峰值所对应的 NOx 扰动 |
| 黑十字 | 基线（0%, 0%）= 本轮 `config_override/emissions.config` 的未扰动排放 |

> **轴序约定**：**VOC 在横轴、NOx 在纵轴**，与 EKMA-FNR（S2）一致，两张图可并排读。

**判读要点**：
- 屋脊把平面分成两半——**屋脊下方（NOx 少）= NOx-limited**：加 NOx 增产 O₃；
  **上方（NOx 多）= NOx-saturated/VOC-limited**：再加 NOx 反而抑制 O₃。
- 沿纵轴向上越过屋脊后 O₃ 下降 = 经典"NOx 滴定"。
- 屋脊随 NMVOC 增加而**上移**（VOC 越多能容纳越多 NOx）——这正是 EKMA 的典型形态。

### 3.1 与 EKMA-FNR（S2）的 O₃max 为什么差 2.4 倍——开放箱 vs 封闭瓶

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

**两者不可比大小，只可比形态**（屋脊位置、体制分界、梯度方向）。要跨图比较，只能比无量纲量：
FNR / VNr 等值线位置、屋脊在各自轴百分比上的形状、峰值率 PR。同一坐标点也不能跨图读——
S1 的 (0,0) 是"参考清单"，S2 的 (0,0) 是"那团空气的初值"，两者毫无共同之处。

> S2 侧的对应分析（含"上限由物理汇而非化学决定"的收支实测与三条无效直觉）见
> `EKMA-FNR/README.md` §6.1。

## 4. 网格与轴语义

| 项 | 本目录（EKMA-emission） |
|---|---|
| 网格 | **26×26 线性**（因子 1.00–7.25，默认 `--pmax 625 --step 25`） |
| 语义 | **扰动百分比**扫描（0% = 无扰动 = 基线） |
| 轴 | x = E_VOC 组乘子，y = E_NOx 组乘子（`build_groups` 的分组定义，含 SO₂/DMS/MSA 在 VOC 组） |
| 产物 | 本目录 `output/`（一张图一轮，CSV 全由 `run.py` 写出） |

执行器 `fac_batch_omp` 的缺省（不传 `LIN`）是 **0.3×1.5^i 等比网格**；本程序显式传
`fmin/fmax + LIN`，两种几何不要混用。

---

## 5. 换上你自己的真实数据：完整操作手册

> **核心原则**：扰动是"相对基线的百分比"（`f = 1 + pct/100`），**基线内容永远来自配置文件**，
> `run.py`/`plot.py` 对基线内容零感知。所以换真实数据 = **换基线文件**，然后重跑。
> **默认幅度 625%、步长 25%（26×26=676）**；历史归档轮次用 400%/10%，改用 `--pmax/--step` 即改变可比锚。

### 5.0 日常输入：编辑 override，网格用 CLI

| 调整项 | 做法 | 影响 |
|---|---|---|
| **基线数据** | 直接编辑本目录已有的 `config_override/` 配置及 `constraints/` 数据 | 只改本程序输入，不影响其它程序 |
| **网格** | `python run.py --pmax ... --step ...` | 只改采样网格；默认 625/25 |

`config_override/` 是已建立、长期保留的**持久输入目录**，不需每轮新建。
它按**文件及相对路径整体使用**（含 `constraints/`），不是按配置行合并，也没有任何根配置模板
可回落。当前已有全套 11 个顶层文件和 17 个 constraints 文件。
`run.py` 把它装配到系统临时目录后，只在那份副本上把 `model.parameters` 的天数/SPINUP 改为 4/3，
再交给 `cfgbuild`；装配副本编译完即清理，`config_override/` 本身一个字节都不动。
日志打印 `[配置] 输入目录: …`；**本轮实际生效的输入**查 `output/runtime_emission.dat`（见 §6）。

**换基线只改本目录 `config_override/`**：三个程序各持一套输入，互不共享配置，也不需重编内核。

### 5.1 九类真实数据 → 配置文件的逐类菜谱

#### ① 你的排放清单（NOx × NMVOC 响应面的"基线"，最常换）

**目标格式** `config_override/emissions.config`（每行：`物种  基线  日型  振幅`）：

```text
# unit: ppb/h                      ← 整份按 ppb/h 读入（推荐，清单口径直观）
NO      4.0    BOX 0.9             # NOx：峰型交通源（BOX 06–18h，amp=0.9）
NO2     0.2    BOX 0.9
CO      10.0   BOX 0.7
HCHO    0.35   BOX 0.9             # ← NMVOC 组 24 物种逐个给
ALD2    0.20   BOX 0.9
ISOP    0.30   BOX 0.9             # 生物源（也可 FLAT）
ALK     1.10   BOX 0.9
...                                  # OLE/ARO1/ARO2/KET/TERP/CSL/LALK/... 全部 24 个
SO2    0.50   BOX 0.6
DMS    0.005  FLAT 0.0
```

**从原始清单换算的三条公式**：

| 你手上的量 | 换算 |
|---|---|
| 排放率 `Q [ppb/h]`（小时均值） | 直接写；BOX 日型的 24h 均值 = `Q×(0.1+0.225×amp)`，峰值 = `Q×(0.1+0.9×amp)`。反推基线：`E0 = 日均/(0.1+0.225×amp)` |
| 年清单 `Tt/yr` 某物种 | `E0[ppb/h] = Tt×1e9/(分子量 g/mol) /8760 /2.46e10 ×面积修正`——0-D 箱没有面积，**需要你先定"箱代表的城市体积/混合层"**，把清单摊成体积源强（这是 0-D 的固有步骤，无默认值可代劳） |
| `μ g/m³` 报告值 | 不适用于排放（那是浓度）；浓度走 ③④ |

**放置**：编辑本目录已有的 `config_override/emissions.config`。
**注意**：24 物种名单必须 ⊆ 机制 VARIABLE（大小写敏感）；E_VOC 组覆盖哪些物种见 §1 的
`build_groups` 口径（SO₂/DMS/MSA 也在组内——本工程既有约定）。

#### ② 你的气象观测（逐时 T/P/RH/BLH/RAIN/JFAC）

**目标格式** `config_override/constraints/met_*.txt`（两列，t 单位秒）：

```text
# 2023-07 观测，站点 XX
0        299.82
3600     298.95
7200     298.10
...                       # ≤128 点！5 天逐时 = 121 点正好
```

**从常见原始数据转换**：

- 气象站 CSV（`date,time,temp_C,...`）→ python 转换示例（**此独立数据准备脚本以仓库根目录为 cwd**，与日常子目录 run/plot 命令不同）：

```python
import pandas as pd
df = pd.read_csv("your_station.csv", parse_dates=[["date","time"]])
df = df[(df.time - df.time.iloc[0]).dt.total_seconds() % 3600 == 0]   # 整点
t = ((df.time - df.time.iloc[0]).dt.total_seconds()).astype(int)
for col, out in [("temp_C","TEMP"), ("press_hPa","PRES"), ("rh_pct","RH"), ("blh_m","BLH")]:
    with open(f"EKMA-emission/config_override/constraints/met_{out}_yours.txt","w") as f:
        f.write(f"# your station {out}
# t[s]  {out}\n")
        f.writelines(f"{a:9d}  {b:9.2f}\n" for a, b in zip(t, df[col]))
```

- 单位红线（cfgbuild 编译期强制，违反即拒绝）：**TEMP 必须开尔文**（°C + 273.15；写过 °C
  会 t=0 崩溃——真实踩坑）、RH 0–100、BLH>0、RAIN≥0；**气象序列禁止 `# unit:` 行**。
- 缺测：cfgbuild 只要求时间严格递增，**值的缺口靠线性插值**（端点钳制外延）——先填好时标，缺值段剔除或内插。
- 时窗对齐：模拟从 t=0 开始跑 5 天 → 你的序列必须覆盖 [0, 432000 s]，起点对齐你想要的"第 1 天 0 点"。

**放置**：在本目录 `config_override/constraints/` 编辑现有序列或放入 `met_XXX_yours.txt`（6 个），
直接编辑已有的 `config_override/meteo.config`，将对应引用改到你的文件名；无需从根模板复制。
引用保持 `file:constraints/文件名` 的相对路径形式，不要指向 `output/`。

#### ③ 你的初始浓度

`config_override/initialConcentrations.config`，每行 `物种 值 [cm3|ppb|ppm]`：

```text
NO      5.0   ppb
NO2     5.0   ppb
O3     30.0   ppb
HCHO   10.0   ppb
CO    200.0   ppb
CH4     1.8   ppm
ISOP    0.5   ppb
...                     # 未列物种 = 0
```

对**排放面**影响小（5 天内 IC 被洗脱）；对 EKMA-FNR（IC 面）是基线本身。

#### ④–⑧ 背景/沉降/稀释/时长/光解

| 数据 | override 文件 | 格式 | 效果 |
|---|---|---|---|
| 背景/边界浓度 | `bc.config`（+序列文件放 `config_override/constraints/`） | `TAU_MIX = 秒`；`O3 = 40 ppb` 或 `O3 = file:...`（序列头部可 `# unit: ppb`） | 背景底色；τ 越小拉得越紧 |
| 干沉降 | `drydep.config` | `O3 0.8 on` | 沉降汇 |
| 湿沉降 | `wetdep.config` | `HNO3 2.0e-5 0.8`（Λ=a·R^b） | 湿沉汇 |
| 稀释通风 | `dilution.config` | `DIL = 1.16e-5`（≈24 h 通风） | 强控制量：实测 0→157 ppb、1.16e-5→99 |
| 时长/步长/容差/SPINUP | `model.parameters` | 5 行：天/步长[s]/RTOL/ATOL/SPINUP；SPINUP<天数 | 建议天数≥3、SPINUP≥1（丢弃建立期） |
| 纬度/日期 | `photolysis.config` | `LAT=40.0`、`DOY=172` | SZA→光解强度 |

#### ⑨ 扰动网格本身（可选）

```bash
# 在 EKMA-emission/ 内执行
python run.py --pmax 300 --step 15    # 21×21，0–300%/15%
python plot.py
```
默认 **625%/25%**（26×26=676）；`400/10` 是历史归档轮次的锚，改动即改变可比性。

### 5.2 完整操作流程（在实验子目录内）

**先决条件**：共享内核已按 §2 在仓库根完成构建，仓库依赖与环境齐备。
下列步骤全部在 `EKMA-emission/` 内执行；日常不再要求 make，配置编译由 `run.py` 自动完成。
任何运行（含小网格）前，若需保留已有结果，先将 `output/` 另存。

```bash
# ── 步 1：编辑已有 override 层（也可保持现有基线直接运行）──
# 编辑: config_override/emissions.config（①）
#       config_override/meteo.config + config_override/constraints/ 中的序列（②）
#       config_override/initialConcentrations.config（③）
#       config_override/ 中的 bc/drydep/wetdep/dilution/photolysis 等（④–⑧，按需）
# 目录已完整建立；本程序的输入只有 config_override/，output/ 不是编辑入口

# ── 步 2：小网格冒烟（5×5，~1 min）——换数据后建议先验证──
python run.py --pmax 40 --step 10
#   验证点:
#   a) 日志出现 "[配置] 输入目录: …/EKMA-emission/config_override"（本实验唯一输入源）
#   b) "FAC-BATCH-OMP ok (N 案例, 失败 0, 降级 0, 氮守恒违例 0, …)"
#   c) O3max 量级合理（城市环境几十~200 ppb；若 O₃≈0 查 IC/排放；若爆表查单位）

# ── 步 3：核对本轮"实际生效"的输入 ──
grep "^MODEL" output/runtime_emission.dat          # 期望 MODEL 4 3600 …（天数/spinup 由实验强制）
grep -A2 "^EMIS" output/runtime_emission.dat       # 前几行 = NO/NO2/CO 的 E0 —— 应等于你写的基线
grep -n "file:constraints" config_override/meteo.config   # 气象引用指向你的序列文件？
#   装配副本编译后即清理，output/ 下不再有 in_emission/：核验输入看 config_override/ 与 runtime 两处
#   排放基线若与你写的不同 → 检查文件名/相对路径是否真被读到（VNr 也由这份 runtime 反算）

# ── 步 4：物理合理性抽查（读 output/ekma_emission_grid.csv）──
#   · (0%,0%) 行 = 你基线下的 O₃ —— 与你的经验/观测对比
#   · FNR_mean 范围应在 0.1–10 量级（脱离物理 = 某轴基线错了）
#   · NOy_ppb 应随 NOx 扰动单调

# ── 步 5：全网格正式跑（默认 676 案例 / 4 天，~8–12 min / 4 核）──
python run.py
#   验证: output/run_meta.txt 的 网格/输入目录/失败·降级·氮守恒 三行；任一非 0 时 run.py 已拒出表，
#        按 stdout 的 [FATAL] 与 [氮守恒 FAIL] 原始行定位案例
#   本步重新装配配置，生成 runtime 和 CSV，覆盖冒烟结果

# ── 步 6：出图 ──
python plot.py
#   → output/ekma_emission_S1.png（响应面+屋脊）
#     output/ekma_emission_slices.png（NOx 切片，看拐点在哪）
#   屋脊数值 ekma_emission_ridge.csv 与派生表 ekma_emission_derived.csv 已由步 5 的 run.py 写出

# ── 步 7：科学判读 ──
#   · 屋脊起点 (NMVOC+0%, NOx+X%)：X<200 = 你城市在采样窗内有 NOx 饱和区；
#     X≥400 = 全窗 NOx-limited（减 NOx 一定降 O₃）
#   · 屋脊随 NMVOC 右移的斜率 = VOC 控制的弹性
#   · FNR 等值线（output/ekma_emission_grid.csv 的 FNR_mean 列）可叠加到 S1 上与观测 FNR 对点

# ── 步 8：与另一基线对比（可选）──
#   先另存本轮 output/ 和对应 override 数据，再更换 config_override/ 中的数据
#   保持网格/机制一致，重新执行 python run.py → python plot.py，再并排比较两轮图件

# ── 步 9：归档 ──
#   另存整个 output/（runtime_emission.dat、CSV、图和 run_meta.txt）
#   并连同一份 config_override/（或用 runtime 记录反查生效值）一起归档
#   同时保存对应输入、网格参数、机制与运行环境信息；复跑仍需完整仓库及共享依赖（§6）
```

### 5.3 常见错误对照

| 症状 | 根因 | 处置 |
|---|---|---|
| CVODE 在 t=0 崩溃 | TEMP 写成了 °C（必须 K） | +273.15 重写序列 |
| cfgbuild 拒绝"气象序列含 unit 行" | 气象文件头写了 `# unit: ppb` | 删掉该行（气象用原生单位） |
| `序列 N 点 > 128` | 观测过密 | 降采样到 ≤128 点（1 h 分辨率 5 天=121 正好） |
| O₃max ≈ 0 | 可能是排放/IC 未按预期落位，或输入目录里的文件名拼写错 | 看 `[配置] 输入目录` 日志；用 `grep -A2 "^EMIS" output/runtime_emission.dat` 核对生效的排放基线，再检查数据本身 |
| O₃max 爆表（>500 ppb） | 排放单位错（ppb/h 值写进了 molec/cm³/s 通道） | 检查 `# unit:` 指令与换算 6.8333e6 |
| 响应面平坦无屋脊 | 基线 NOx 太低（背景级）或太高 | 按 §5.1③/①把基线调到城市量级 |
| 其它程序的输出意外变化 | 三个程序共用机制/执行器但各持 `config_override/`——多半是共用的 `.fac` 或 `modules/*.f90` 被动过 | 核对本目录 `output/run_meta.txt` 的机制 md5；本程序数据改在 `config_override/` |
| `字段超限` 类 csv 报错 | 宽表行 >128 KB | 已在 csvio 防御（16 MB）；自定义脚本需同样设置 |

### 5.4 不能替换 / 不能随便动的（动了等于换实验）

| 项 | 为什么不能动 |
|---|---|
| **机制文件** `mechanism/mechanism_3_full_chain.fac` | **V2.5+ 已封版**。换机制 = 换化学 = 另一个实验体系（走 V2.6 流程） |
| **`run.py` 的缩放逻辑**（`f = 1 + pct/100`） | "百分比扰动"的定义本身；改了它 x/y 轴语义即失效 |
| **`build_groups` 分组映射**（Fortran 硬编码） | E_NOx = NO,NO₂；E_VOC = 24 物种（§1）。改分组（如把 SO₂/DMS 踢出 E_VOC）要改 `runtime_dat_mod.f90` 并重编译 |
| **执行器的缺省等比网格** | 不传 `LIN` 时 `fac_batch_omp` 走 `0.3×1.5^i`；本程序显式传 `fmin/fmax + LIN`，别把两种几何混在一轮里 |
| **默认幅度 625% / 步长 25%** | 本轮默认。改网格走 `--pmax/--step` 并在 `run_meta.txt` 留记录，跨网格不可直接比 |

## 6. 输入装配与"有效输入记录"（为什么 output 里没有配置快照）

> 与姊妹实验 EKMA-FNR 同构（详见 `EKMA-FNR/README.md` §3）。

每次**完整运行**（非 `--from`）的配置链为：

```
本目录 config_override/（唯一输入源）
    → 临时装配目录（副本上把天数/SPINUP 改为 4/3）→ cfgbuild → output/runtime_emission.dat → fac_batch_omp
```

**装配目录是必需的**：`cfgbuild` 需要"一整套配置"的目录作输入，而实验口径（4 天/3 天）
**不能写回你正在编辑的 `config_override/`**。它在系统临时目录完成、编译后立即清理，
所以 `output/` 不再有 `in_emission/`——避免出现"看着像输入入口、其实每轮被重建"的误导。

| 记录本轮实际用了什么 | 内容 |
|---|---|
| `output/runtime_emission.dat` | cfgbuild 解析结果：`IC` 向量、`EMIS` 各行 E0/日型/振幅、`MODEL 4 3600 …3`、`MET/SERIES/PHOT/DRY/WET/DIL/BC` 的数值本体——**唯一完整的"有效输入"事实源**；`VNr_group/VNr_nmv` 的 N0/V0 也是从它反算 |
| `config_override/` | 你编辑的输入原文（除天数/SPINUP 由实验强制，其余逐字生效） |
| `output/run_meta.txt` | 本轮网格、时长、耗时、机制 md5、`输入目录`、`基线总量 runtime`、`有效输入记录`、失败/降级数、**氮守恒违例数**、**侧车单位口径** |

> **有效输入记录**：每轮都由本目录装配 → `cfgbuild` → `output/runtime_emission.dat`，
> 该文件就是本轮真正喂进执行器的内容（IC/EMIS/MODEL/MET 的解析值）。
> 仓库里不存在全局 runtime，也没有需要回退的共享配置。

### 操作与输出独立，不等于脱离仓库

1. 不 import 任何编排层代码（配置准备、cfgbuild 调用、batch 调用都在本目录 `run.py` 内实现）；
2. 产物只落本目录 `output/`，不读不写其它程序的 `output/`；
3. 日常在本目录编辑 `config_override/`，执行 `python run.py` → `python plot.py`，
   不重新 `make` 共享内核。

**共享依赖 = "模型本体"四件套**（必要耦合，非杂物）：

| 共享件 | 角色 | 为何必要 |
|---|---|---|
| `mechanism/mechanism_3_full_chain.fac` | 化学机制（V2.5+ 冻结） | 实验对象 |
| `facsimile_impl/fac_batch_omp` | 预先构建的批量积分执行器（EMIS 模式/LIN 网格） | 唯一数值引擎 |
| `inputproc/cfgbuild.py` | 配置 → runtime 编译器 | 每次完整 run 从临时装配目录生成配置运行时，并非重编内核 |
| `analysis_lib/csvio.py` | CSV 读写（Fortran 指数兜底单源） | 数据契约读写 |

外加 Python 与依赖包、执行器所需动态库（CVODE/OpenMP/gfortran 运行时）与兼容环境。
**配置链已完全本地化**：本程序只读 `config_override/`，只写 `output/runtime_emission.dat`。

**归档与复现边界**：
- **一轮实验的完整留档 = `config_override/`（输入原文）+ `output/`（runtime 有效输入记录 + 数据 + 图 + run_meta）**；
  两者一起另存，才能既核验"生效了什么"又重放"用什么跑的"；
- **复现/重跑仍须在完整仓库内**：把要复现的那套输入放回 `config_override/`，配同一机制、
  内核与网格再跑；`output/` 或 `config_override/` 单独拷走都不能脱离仓库运行。

**若将来要把本实验拆出 ozone_lab**：本目录 + 上表四件套 + `modules/`（Fortran 物理模块）+
`facsimile_impl/gen_v3/` + `facsimile_impl/mechc.py`（执行器二进制的编译依赖）。
---

## 7. `output/` 数据字典（逐 CSV 逐列）

> 通用约定：浓度列单位 **ppb**；过程量（Q/L_NOx/P_*/IRR）单位 **molec/cm³/s**；
> "窗口" = 每案例模拟的**末 24 步**；每行 = 一个正交扰动案例（pert_nox × pert_voc 的一个组合）。

### 7.1 `ekma_emission_grid_raw.csv` —— fac_batch_omp 原始输出（表头 + n² 行；当前归档 26×26=676）

| 列 | 类型 | 含义 |
|---|---|---|
| `fnox` | 输入因子 | NOx 排放组缩放因子 = 1 + pert_nox_pct/100（本表 1.00–5.00） |
| `fvoc` | 输入因子 | NMVOC 排放组缩放因子 = 1 + pert_voc_pct/100 |
| `jfac` / `dt` | 输入因子 | 光解/温度维（本实验恒 1.0 / 0.0，未启用） |
| `O3max` [ppb] | **主响应** | 窗口内 O₃ 最大值（主图色标；`grid.csv` 中改名 `O3max_ppb`） |
| `O3mean` [ppb] | 响应 | 窗口 O₃ 均值 |
| `HCHOmax` [ppb] | 响应 | 窗口 HCHO 峰值 |
| `NO2mean` [ppb] | 状态量 | 窗口 NO₂ 均值 |
| `FNR_mean` | 判据 | HCHO/NO₂ 窗均比（Sillman 1995：<1 VOC-limited，1–2 过渡，>2 NOx-limited） |
| `H2O2HNO3` | 判据 | H₂O₂/HNO₃ 窗均比（>0.2 → NOx 控制区） |
| `O3NOy` / `HCHONOy` / `O3NOz` | 备用判据 | O₃/NOy、HCHO/NOy、O₃/NOz（占列，不进主判据） |
| `NOy_ppb` | 状态量 | 总活性氮窗均（25 物种加权，N₂O₅ 计双 N） |
| `LnOxQ` | 过程量 | L_NOx/Q（Kleinman 2005：>0.5 → NOx 饱和） |
| `PH2O2PHNO3` | 过程量 | P(H₂O₂)/P(HNO₃) 产率比 |
| `OPE` | 过程量 | O₃ 生成效率 = P(O₃)/P(NOz) |
| `FNR_pm` 等 6 个 `_pm` 列 | 判据(午后) | 上列对应量的 **12–18 h 午后口径**（Sillman 原文窗口） |
| `nst` | 诊断 | CVODE 内步数（跨韧性重启累计；异常小 = 该案例收敛失败已降级，见 stdout 警告） |

### 7.2 `ekma_emission_grid_raw.csv_ipr.csv` —— 过程账侧车（4 网格列 + 111 物种 × 6 账 = 670 列；单位 = 末 24 步**平均倾向** [molec/cm³/s]，**不是** CASE 那套窗口积分）

| 列模式 | 含义（窗口均值，molec/cm³/s） |
|---|---|
| `<物种>_P` | 机制化学生**成**速率（正） |
| `<物种>_L` | 机制化学**消耗**速率（**输出为负值**） |
| `<物种>_E` | **排放注入**（EMIS 模式下非零且随因子缩放——与 IC 模式恒 0 区分，可作模式自检） |
| `<物种>_D` | 干沉降 + 湿沉降（负） |
| `<物种>_B` | BC 背景张弛 + 稀释（可正可负） |
| `<物种>_X` | 边界层膨胀稀释（负，仅增长段生效） |

用途：把 O₃ 响应差分解到"哪个**过程**造成"（如屋脊右支 O₃ 下降 = `_L` 滴定项涨幅压过 `_P`）。
物种清单 = 机制 VARIABLE 全集（111）；列序 ↔ `runtime_emission.dat` 的 SPECIES 行。

### 7.3 `ekma_emission_grid_raw.csv_irr.csv` —— 反应级侧车（4 网格列 + `r1..r203` = 207 列）

`r<k>` = 机制 `.fac` **第 k 条反应**的 IRR 窗口均值（molec/cm³/s）。反应行序与
`facsimile_impl/gen_v3/oz.reac` 一致（查 `r53` 是哪条反应：按行序数 `.fac` 的 `%` 行）。
用途：把 O₃ 响应差归因到**具体反应**（比 IPR 细一档；如高 NOx 角 NO₂+OH→HNO₃ 通道）。

### 7.4 `ekma_emission_grid_raw.csv_state.csv` —— 状态侧车（27 列，本实验未消费，备查）

| 列 | 含义 |
|---|---|
| 网格 4 列 | fnox, fvoc, jfac, dt |
| `NO2mean_ppb` `HCHOmean_ppb` `HNO3mean_ppb` `H2O2mean_ppb`（+`_pm` 版） | 四状态量**同窗一致均值对** |
| `Q` `L_NOx` `P_H2O2` `P_HNO3` `P_O3` `P_NOz`（+`_pm` 版） | 过程量（自由基周转 / NOx 损耗 / 产率） |
| `O3max_ppb` `O3mean_ppb` | 与主 CSV 同 |
| `nst` | 内步数 |

用途：实验 J"状态空间投影"（轴=状态量）的原始数据。EKMA-emission 的轴是因子本身不需要它；
留作"排放扰动 → 状态空间落点"的后续分析。

### 7.5 `ekma_emission_grid.csv` —— 规整主表（**plot.py 唯一输入**）

| 列 | 含义 |
|---|---|
| `pert_nox_pct` | NOx 排放扰动 % = (fnox−1)×100，0–400 步长 10，**y 轴** |
| `pert_voc_pct` | NMVOC 排放扰动 % = (fvoc−1)×100，0–400 步长 10，**x 轴** |
| `O3max_ppb` / `O3mean_ppb` | raw 的 `O3max/O3mean` 改名（语义同上） |
| 其余 17 列 | 与 raw 同名同义（HCHOmax … nst） |
| `VNr_group` | **该案例的排放 VOC/NOx 比（摩尔比）** = (f_VOC·V₀,group)/(f_NOx·N₀)，VOC 组含 SO₂/DMS/MSA<br>（V₀,group/N₀=3.707；**不可当质量比**——集总类无单一摩尔质量） |
| `VNr_nmv` | 同上但 V₀ 剔除 SO₂/DMS/MSA（更接近 NMVOC） |
| 注释头（`#` 行） | 网格口径 / 两轴定义 / 固定条件——自记录，读表先看注释 |

### 7.6 `ekma_emission_ridge.csv` —— 屋脊轨迹（表头注释 + 26 行）

| 列 | 含义 |
|---|---|
| `pert_voc_pct` | NMVOC 扰动列（当前归档 0–625，25% 步长） |
| `ridge_pert_nox_pct` | **该列 O₃max 峰值所在的 NOx 扰动 %**——体制分界：峰左侧 NOx-limited，右侧滴定区 |
| `ridge_O3max_ppb` | 该峰值的 O₃max |

### 7.7 `run_meta.txt`

模式（完整运行/复用）/ 网格（由数据推断）/ 时长（days/spinup 取自 runtime）/ 耗时 / 机制路径+md5 / 基线总量 runtime /
输入目录 / 有效输入记录（`output/runtime_emission.dat`）/ 失败与降级计数 /
**氮守恒违例数**（逐案例首快照检，违例即 `stop 1`；原始 `[氮守恒 FAIL]` 行由 `run.py` 捕获入库）/ **侧车单位口径**。

> `nfail>0`（零统计求解失败）**或** `氮守恒违例>0` 时执行器末尾 `stop 1`，`run.py` **一律拒绝出表**；
> 不变量是逐案例的首快照检查（每案例一次，确定性的），当前归档实测违例 0。
> `nwarn`（降级）不致命，但必须看 `run_meta.txt` 判断是否可用。

### 7.8 `ekma_emission_derived.csv` —— 派生表：PR 峰值率与 VNr 临界比（表头注释 + 26 行）

> **🔒 固化声明**：本导出**内置于 `run.py` 源码**（`emis_totals()` + 规整段 + 派生段）——
> 每次运行（含 `--from` 复用模式与 `config_override/` 真实数据模式）**自动产出**
> `grid.csv`（含 `VNr_group` / `VNr_nmv` 两列）与 `ekma_emission_derived.csv`，不是临时脚本。
> 概念定义见下表与 §5（PR/VNr 派生量详解）。

| 列 | 含义 |
|---|---|
| `pert_voc_pct` | NMVOC 扰动列（当前归档 0–625，25% 步长） |
| `O3peak_ppb` | 该列沿 NOx 轴的 O₃max 峰值（屋脊点；即 ridge_O3max_ppb 同值） |
| `PR_bau` | **峰值率①** = O₃peak / O₃max(基线 0%,0%)——排放同比例放大的"峰值收益" |
| `PR_vs_baseNOx` | **峰值率②** = O₃peak / O₃max(同列 NOx=0% 端)。⚠ **0% 扰动 = 基线排放，不是零排放**——它是「从基线 NOx 爬到屋脊」的放大，不是「从无 NOx 起」 |
| `ridge_pert_nox_pct` | 屋脊 NOx 扰动（与 `ekma_emission_ridge.csv` 的同名列同值；**越界列 = 采样上界（当前 625）**，即真实峰值在边界之外、未达到） |
| `VNr_ridge_group` | **屋脊处排放 VOC/NOx 摩尔比** = (1+VOC%)/(1+ridgeNOx%) × V₀,group/N₀ |
| `VNr_ridge_nmv` | 同上但 V₀ 剔除 SO₂/DMS/MSA |
| `ridge_at_edge` | 1 = 该列峰值落在 NOx 采样上界（真实峰可能在界外，外推无效） |
| `O3_at_baseNOx_ppb` | 同列 NOx=0% 端（=基线排放）的 O₃max |

**实测与解读**（2026-10-29 归档，26×26 / 0–625%／25%）：

| NMVOC + | O₃峰 | PR_bau | 屋脊 NOx+ | VNr_ridge_group |
|---|---|---|---|---|
| 0% | 94.8 | 1.23 | 125% | **1.65** |
| 50% | 127.6 | 1.65 | 225% | **1.71** |
| 100% | 156.4 | 2.03 | 325% | **1.74** |
| 200% | 206.1 | 2.67 | 525% | **1.78** |
| ≥300% | 245.4→281.4 | 3.18→3.65 | 越界（=625%） | (截断值) |

**两个科学结论**：
1. **临界排放比近常数**：屋脊未越界的列（NMVOC 0–225%，10/26 列）VNr_ridge_group ≈ **1.65–1.79**——经典 EKMA
   理论的"临界 VOC/NOx 排放比"在本机制（含硫集总口径、摩尔比）下 ≈ **1.7**；
   排放组成低于 1.7 → VOC-limited，高于 → NOx-limited。
2. **PR_bau 单调上升**（1.23→3.65）：VOC 越丰富屋脊峰值放大越多；≥300% 列屋脊越界，
   PR 为下界（真实峰值更高）。

---

## 8. 可复现性

- 机制 md5、runtime 路径、网格参数写入 `output/run_meta.txt`；
- CSV 表头注释记录扰动口径与固定条件；
- 每次重跑幂等覆盖本目录 `output/`（不触碰仓库其它产物）。

---

## 9. 最小化运行 vs AtChem2，与"站点在 S1 上的落点"

### 9.1 最小化运行的三要素（与 AtChem2 概念一致）

| | AtChem2 | ozone_lab |
|---|---|---|
| 机制 | MCM Fortran 源码（先编译） | `.fac` → mechc.py → gen_v3（`make all` 先编译） |
| 初始浓度 | `initialConcentrations.config` | `config_override/initialConcentrations.config`（未列=0） |
| 气象/光解 | `environmentVariables.config`（常数或时序）+ 光解内置于发行包 | `config_override/meteo.config` + `photolysis.config` + `photolysis-rates_MCM`（**光解是显式必需输入**） |
| 运行控制 | `model.parameters` | `config_override/model.parameters`（天数/步长/RTOL/ATOL/SPINUP；本程序强制 4 天/3 天） |

**判定**：概念三要素一致（机制 + 初始场 + 气象/光解环境）。文件层面差异：ozone_lab 把光解拆为
**显式必需输入**（语义最小集 = 5 类）；且 `cfgbuild` 无条件打开 10 份配置——**文件层面 10 份必须
存在，其中 6 份（排放/BC/沉降×2/稀释/scenarios）可以全空**（已实测全空照常出结果）。

### 9.2 选择性输入的精确语义：**未列出 = 零/无**（不是"内建经验值顶上"）

| 通道 | 不输入时 | 语义 |
|---|---|---|
| `emissions.config` | 排放 = 0 → **纯初值化学** | "经验值"在随仓库分发的参考清单里（本图 0% 点的基线），是**数据文件**而非代码兜底 |
| `bc.config` | 无约束 → 自由运行 | 一通道两用：常数=背景场；`file:`=观测约束（nudging，τ 控制） |
| `drydep`/`wetdep` | 沉降 = 0 | — |
| `dilution.config` | DIL = 0（密闭盒子） | 参考配置给 1.16e-5（24h 通风） |
| `scenarios.config` | 仅 DEFAULT 情景（全因子=1） | — |

### 9.3 站点在 S1 上的落点：清单只用来算坐标，不必灌进模型

> **先给结论**：**排放基线不用动，清单也不必作为模型输入**——它只做两件事：
> ① 算站点的落点坐标；②（可选）替换基线，让响应面围着你的城市展开。

#### 9.3.1 参考基线的两个总量（落点分母，逐物种实测）

**折算口径**：`emissions.config` 给的是**基准值 E0**（molec/cm³/s），实际排放还乘日型：

```
BOX  : 06–18h  E = E0·(0.1 + 0.9·amp·sin²(π(h−6)/12)) ；夜间 E0·0.1
FLAT : E = E0（全天恒定）
24h 均值 = E0·(0.1 + 0.225·amp)      ← amp=0.9 → 0.3025·E0；amp=0.75 → 0.26875·E0
峰值     = E0·(0.1 + 0.9·amp)        ← amp=0.9 → 0.91·E0
ppb/h 换算: × 2.46e10 / 3600 = ×6.8333e6
```

**① NOx 组（2 物种，全部 BOX/amp=0.9）**

| 物种 | E0 [molec/cm³/s] | 24h 均值 [ppb/h] | 占 E_NOx |
|---|---|---|---|
| NO（一氧化氮） | 2.157e7 | 0.9549 | 95.2% |
| NO₂（二氧化氮） | 1.078e6 | 0.0477 | 4.8% |
| **E_NOx 合计** | | **1.0026 ≈ 1.003** | 100% |

**② E_VOC 组（有排放行的 21 物种，按贡献降序）**

| 物种 | 中文名/含义 | E0 [molec/cm³/s] | 日型/amp | 24h 均值 [ppb/h] | 占比 |
|---|---|---|---|---|---|
| ALK | 烷烃类（集总，C₃代表） | 1.573e7 | BOX/0.90 | 0.6963 | 18.7% |
| OLE | 烯烃类（集总，丙烯代表） | 1.494e7 | BOX/0.90 | 0.6614 | 17.8% |
| SO₂ | 二氧化硫（⚠ 组内但非 VOC） | 1.128e7 | BOX/0.90 | 0.4991 | 13.4% |
| LALK | 长链烷烃（C₆代表） | 1.049e7 | BOX/0.90 | 0.4644 | 12.5% |
| ARO1 | 甲苯类芳烃（集总） | 1.019e7 | BOX/0.90 | 0.4511 | 12.1% |
| ARO2 | 二甲苯类芳烃（集总） | 7.992e6 | BOX/0.90 | 0.3538 | 9.5% |
| ETH | 乙烯 | 3.201e6 | BOX/0.90 | 0.1417 | 3.8% |
| IOLE | 内烯烃（集总） | 3.201e6 | BOX/0.90 | 0.1417 | 3.8% |
| ISOP | 异戊二烯（生物源） | 2.037e6 | BOX/0.90 | 0.0902 | 2.4% |
| CSL | 甲酚类（芳香开环前体） | 1.798e6 | BOX/0.90 | 0.0796 | 2.1% |
| KET | 酮类（丙酮代表） | 1.500e6 | BOX/0.90 | 0.0664 | 1.8% |
| MEOH | 甲醇 | 5.000e5 | BOX/0.90 | 0.0221 | 0.6% |
| ALD2 | 高级醛（丙醛代表） | 3.000e5 | BOX/0.90 | 0.0133 | 0.4% |
| ETOH | 乙醇 | 3.000e5 | BOX/0.90 | 0.0133 | 0.4% |
| ESTER | 酯类（乙酸乙酯代表） | 2.000e5 | BOX/0.90 | 0.0089 | 0.2% |
| MTBE | 甲基叔丁基醚 | 1.000e5 | BOX/0.90 | 0.0044 | 0.1% |
| TERP | 单萜烯（α-蒎烯代表） | 7.500e4 | BOX/0.75 | 0.0029 | 0.1% |
| HCHO | 甲醛 | 5.815e4 | BOX/0.90 | 0.0026 | 0.1% |
| DMS | 二甲基硫（海洋源） | 1.000e4 | FLAT/— | 0.0015 | 0.0% |
| SESQ | 倍半萜烯（β-石竹烯代表） | 2.500e4 | BOX/0.75 | 0.0010 | 0.0% |
| CH3CHO | 乙醛 | 1.744e4 | BOX/0.90 | 0.0008 | 0.0% |
| **E_NMVOC 合计** | | | | **3.7164 ≈ 3.716** | 100% |

> **组内无排放行的 3 物种**：**MVK**（甲基乙烯基酮）、**MACR**（甲基丙烯醛）、**MSA**（甲磺酸）
> ——它们是二次产物，`emissions.config` 里没有对应行（恒 0）。比对时**既不计入分母，也不应
> 计入你的清单**。
>
> **⚠ 口径提示**：E_VOC 组名义 24 物种但含 **SO₂（0.4991 ppb/h，占 13.4%）**——硫排放也在
> E_VOC 组内，这是本工程既有的集总约定（`build_groups`，Fortran 硬编码）。你若只统计
> "真正的 NMVOC"，需自行决定是否剔除 SO₂/DMS/MSA：**用同一口径比两边即可**，但在论文中要
> 写明"E_VOC 含硫组分"。

**③ 关键比值**：**E_NMVOC / E_NOx = 3.716 / 1.003 = 3.707**（后文 VNr 计算用 V₀/N₀=3.707）。

**④ 量级校验**：NOx 1.003 ppb/h @ BLH 1000 m
= 6.85e11 molec/cm²/s = **16.5 t/km²/yr**（以 NO₂ 计）——中国东部城市量级 ✓（对照 UNITS §2.8b）。

**⑤ 复算脚本**（基线若改过，用它刷新数字）：

```bash
python - <<'EOF'
GRP = "HCHO CH3CHO ISOP DMS MVK MACR MSA SO2 ALK OLE ARO1 ARO2 KET TERP "\
      "CSL LALK MEOH ETOH MTBE ESTER SESQ ETH IOLE ALD2".split()
PPB = 2.46e10/3600
nox = voc = 0.0
for ln in open("EKMA-emission/config_override/emissions.config"):
    s = ln.split("#")[0].split()
    if len(s) < 2: continue
    try: e0 = float(s[1])
    except ValueError: continue
    amp = float(s[3]) if len(s) > 3 else 0.0
    day = e0*(0.1+0.225*amp) if (len(s) > 2 and s[2] == "BOX") else e0
    if s[0] in ("NO", "NO2"): nox += day
    if s[0] in GRP: voc += day
print(f"E_NOx = {nox/PPB:.3f} ppb/h ; E_NMVOC = {voc/PPB:.3f} ppb/h ; "
      f"比 = {voc/nox:.3f}")
EOF
```

#### 9.3.2 落点公式与实算例子

**公式**（两行）：

```
pert_NOx[%] = (E_NOx,站点 / 1.003 − 1) × 100
pert_VOC[%] = (E_NMVOC,站点 / 3.716 − 1) × 100
```

**例 1（完整分步）**：你的清单 NOx = 2.5 ppb/h、NMVOC = 8.0 ppb/h

| 步骤 | NOx 轴 | NMVOC 轴 |
|---|---|---|
| ① 你的清单值 | 2.5 ppb/h | 8.0 ppb/h |
| ② 除以参考基线 | 2.5 / 1.003 = **2.492** | 8.0 / 3.716 = **2.153** |
| ③ 减 1，乘 100 | (2.492 − 1) × 100 = **+149%** | (2.153 − 1) × 100 = **+115%** |
| ④ 落点坐标 | **y = +149%** | **x = +115%** |

在 `ekma_emission_S1.png` 上落点 **(x=+115%, y=+149%)**：

| 查表/读图 | 值 | 来源 |
|---|---|---|
| 该点的 O₃max 预测 | **≈ 142 ppb** | `grid.csv` 的 (110,150)=141.2 / (120,150)=143.0 插值 |
| 同列（NMVOC+115%）的屋脊 NOx | **≈ +355%** | `ekma_emission_ridge.csv`：100%→320%、110%→340%、120%→370% 插值 |
| 该列峰值 O₃ | ≈ 166 ppb | 同上 |

**判读**：你的点（NOx +149%）在屋脊（NOx +355%）的**左侧下方** → 该城市处于
**NOx-limited** 区：继续加 NOx 还能增产（从 142 → 166 ppb），**减 NOx 则会降 O₃**。
距屋脊还很远（差 206 个百分点）→ 结论稳健，模式 ① 的近似足够（无需换基线）。

**例 2（超窗情形）**：清单 NOx = 6.0 ppb/h、NMVOC = 20 ppb/h
→ NOx: 6.0/1.003 = 5.98 → **+498%**；NMVOC: 20/3.716 = 5.38 → **+438%**
→ **超出 0–400% 窗**（右上方外）。处置：`--pmax 600 --step 15` 扩窗重跑，或先用模式 ②
把基线换成该清单。

**例 3（恰好等于参考基线）**：NOx 1.003 / NMVOC 3.716 ppb/h → **(0%, 0%)**，落在原点。

**三个必须记住的点**：
- **原点 (0%,0%) = 参考基线**（本轮 `config_override/emissions.config`；生效值见 `runtime_emission.dat`），**不是你的站点**；
- 站点**观测浓度**（O₃/HCHO/NO₂）**不能直接落 S1**——浓度 ≠ 排放；这正是 S2（FNR 平面，
  轴=可观测量）存在的理由：S1 回答"排放政策空间"，S2 回答"观测状态空间"；
- 落点的**体制判读**读的是"点在屋脊的左/右"：左侧 = NOx-limited（加 NOx 增 O₃），
  右侧 = VOC-limited/滴定区（加 NOx 降 O₃）。

#### 9.3.3 三个口径必须对齐（否则比值失真——唯一的技术要害）

| 口径 | 参考基线 | 你要做的 |
|---|---|---|
| **单位** | 体积源强 molec/cm³/s → 折 ppb/h | 两边同单位即可，**换算因子自动抵消**。若你的清单是 `kg/m²/s` 或 `Tt/yr`：把参考基线换成面通量 `E_area = E_vol × BLH[cm]`（默认 BLH 1000 m），或把你的清单换成体积源强——**关键是两边同一箱高假设**。量级参考：NOx 1.003 ppb/h ≈ **16.5 t/km²/yr**（以 NO₂ 计、BLH 1000 m） |
| **物种范围** | E_VOC 组 **21 物种**（MVK/MACR/MSA 无排放行） | 只加这 21 类；你的清单里 PAH、有机酸、CO₂、NH₃ 等**不属 E_VOC，不能计入**。若你只有"NMVoC 总量"，需按源谱拆分并扣除未覆盖物种，否则系统性高估 |
| **时间** | **24h 均值**（非峰值） | 日均比日均、峰比峰。清单若是年均，先除到日均（或两边都用年均——比值不变，只要一致） |

#### 9.3.4 两种模式：什么时候才需要动基线

| 模式 | 做法 | 优点 | 代价 |
|---|---|---|---|
| **① 参考面 + 落点**（默认，推荐先做） | 基线不动，清单只算坐标 | 所有城市落在**同一张图**上可横向比较；零额外计算 | 你的城市源谱/VOC:NOx 若与参考差异大，**屋脊位置有偏差**（落点是近似） |
| **② 基线换你的清单** | 在 `config_override/emissions.config` 写入你的清单 → 重跑全网格（默认 676 案例） | 0% = 你的城市，整张面按你的化学算，落点精确 | 只能与自己的历史比，不再跨城市可比；多跑一轮 |

**建议策略**：
1. 先用模式 ① 算落点。若落点**远离屋脊**（如与屋脊的 NOx 差 >50%）→ 结论稳健，不必换基线；
2. 若落点**紧贴屋脊**、或你的清单与参考差 **>2 倍**、或源谱明显不同（如芳香烃占比高很多）
   → 切模式 ② 复核（图会变；两模式结论一致即可放心）；
3. **网格下界是 0%（`fmin=1.0`），没有负半轴**：两轴比值 **<1 → 负百分比 → 落在图外**；
   比值 **>5 → 超出 +400%**。这两种情况你的城市与参考差太多 → 上模式 ②（或 `--pmax` 扩窗，
   但扩的只是上界，负半轴补不回来）。

#### 9.3.5 没有清单时的替代路径

| 路径 | 做法 | 局限 |
|---|---|---|
| 卫星柱 → 排放反演 | `inputproc/sat2surf.py` 转表面 ppb + 自上而下比值法 | 反演为一步比值法，非严格源解析（见 OPERATIONS §10） |
| 改落到 S2 | 观测 HCHO/NO₂ 直接落 FNR 平面（`EKMA-FNR/README` §11） | 回答的是"观测状态"而非"排放政策"问题 |
| 用文献清单 | 借用同类型城市清单（如 MIX/MEIC 截面） | 需说明来源与代表年 |

### 9.4 额外定位法（✅ 已实现）：仅用 **O₃ + NO₂ 观测** 反查排放坐标（无需 VOC 清单）

> **动机**：§9.3 的方法需要 NOx **与** NMVOC 两份清单。但 O₃ 与 NO₂ 是**常规监测**（站点小时值、
> 卫星柱），若手上没有 VOC 清单，仍可用"观测反查"在 S1 上定位——利用本实验 `_state.csv`
> 侧车已记录的每个案例的 **(NO₂ 窗均, O₃max)** 做查表反演。
> **本节为额外方法，不替代 §9.3 的清单定位**（后者无歧义；前者数据更易得）。

#### 9.4.0 与 EKMA-FNR §11.7/§11.8 的关系

| | EKMA-FNR §11.7/§11.8 | 本节（EKMA-emission） |
|---|---|---|
| 反查语义 | 初值状态 → "这团空气会怎样" | **排放情景 → 政策空间坐标**（与 §9.3 清单法同轴可比） |
| 匹配空间 | NO₂/HCHO 均值 + O₃max（1 天封闭瓶） | NO₂ 均值 + O₃max **(+HCHO 均值可选)**（4 天开放箱稳态） |
| 工具 | `plot.py --obs-single/--obs-more` | `plot.py --obs-single/--obs-more [--with-hcho]`（**同名同构**：同一套三档图面标记 + 同一套终端候选表，`draw_station`） |

#### 9.4.1 用法

```bash
python EKMA-emission/plot.py --obs-single "O3=142,NO2=6.5"                    # 两观测反查
python EKMA-emission/plot.py --obs-single "O3=142,NO2=6.5,HCHO=8.0" --with-hcho   # 三观测（VOC 轴更准）
python EKMA-emission/plot.py --obs-more stations.csv --with-hcho              # 多站点
#   stations.csv 列: station,O3,NO2[,HCHO]（一行一站, ppb）
```

观测口径（**必须与模拟统计窗一致**）：`O₃` = 日最大；`NO₂` = 窗口均值（日均）；`HCHO` = 窗口均值。

**图面标记三档**（`draw_station`，与 EKMA-FNR §11.8 同规则）：

| 站点数 | 红星 | 8 近邻候选（终端表始终给全明细） | 文字标注 |
|---|---|---|---|
| **1 站** | ✅ s=280 + 图例 `station (matched)` / `8-NN candidates` | ✅ **包络框**（候选 VOC%/NOx% 的红点矩形）+ 叉点 | 坐标 + 屋脊侧短标签（`ridge-left: NOx-limited` / `ridge-right: NOx-saturated`）+ `match σ`；两观测时附 `2-obs: VOC axis ambiguous ±40–50%` |
| **2–5 站** | ✅ | ✅ 仅**叉点**（多框会互相遮挡） | `编号: 站名 (VOC% / NOx%)` |
| **>5 站** | ✅ | ❌ 不画（避免糊图） | 只标**编号**，明细在终端表 |

```text
终端: 多站汇总表 + 每站 8 近邻候选表（rank/距离σ/VOC%/NOx%）
      plot.py 只写 PNG，候选明细不落 CSV（全仓 CSV 由 run.py 独占）
```

⚠ **标记会画进 `ekma_emission_S1.png` 本身**（覆盖同名文件）。存档需无标记响应面时，
重跑 `python plot.py`（不带 `--obs-*`）即可恢复。
前置条件：`output/ekma_emission_grid.csv` 与 `ekma_emission_grid_raw.csv_state.csv` 存在（先跑 `run.py`）。
**站名请用 ASCII**：图件不配 CJK 字体，中文站名会渲染成缺字方块（脚本会告警；中文明细仍在终端）；
标注位置自动避让（右半区向左标、贴底抬高），不会压住右下角图例。

#### 9.4.2 反查原理与脚本（--obs-single 实测输出对照）

原理：S1 的 676 案例，每个都有两组坐标——**控制坐标** (pert_voc_pct, pert_nox_pct)（图的轴）与
**状态坐标** (NO2mean_ppb, O3max_ppb)（`_state.csv`，皆可观测）。给定观测 (NO₂, O₃) → 标准化
(可加 HCHO 第三维) → 找最近案例 → 读其控制坐标。**NO₂ 随 NOx 单调 → NOx 轴定位准；VOC 轴
靠 O₃ 形状分辨 → 有歧义（±40–50%）**，加 HCHO 第三维可压缩。

```python
import sys; sys.path.insert(0, "analysis_lib")
import numpy as np
from csvio import read_table, fcol

hg, G = read_table("EKMA-emission/output/ekma_emission_grid.csv")
hs, S = read_table("EKMA-emission/output/ekma_emission_grid_raw.csv_state.csv")
pv, pn = fcol(hg, G, "pert_voc_pct"), fcol(hg, G, "pert_nox_pct")
no2, o3 = fcol(hs, S, "NO2mean_ppb"), fcol(hs, S, "O3max_ppb")

NO2_obs, O3_obs = 6.5, 142.0                       # 你的观测（同窗口口径）
cols = np.c_[no2, o3]
M, SD = cols.mean(0), cols.std(0)
q = (np.array([NO2_obs, O3_obs]) - M) / SD
P = (cols - M) / SD
d = np.linalg.norm(P - q, axis=1)
idx = np.argsort(d)[:10]
print(f"最近距离 {d[idx[0]]:.3f}σ")
print(f"VOC 轴候选: {pv[idx].min():.0f}–{pv[idx].max():.0f}%  (中位 {np.median(pv[idx]):.0f}%)")
print(f"NOx 轴候选: {pn[idx].min():.0f}–{pn[idx].max():.0f}%  (中位 {np.median(pn[idx]):.0f}%)")
```

#### 9.4.3 四条硬约束（自动检查/警告已内置）

| 约束 | 说明 |
|---|---|
| **统计口径** | 模型端 `O3max_ppb` = 末 24 输出步窗口最大、`NO2mean_ppb` = 窗口均值；观测端用对应口径（O₃ 日最大、NO₂ 日均），不可用峰值 NO₂ 比均值 |
| **模型依赖性** | 反查基于"本机制 + 参考源谱"算出的 NO₂–O₃ 关系；你城市源谱不同 → 状态关系偏移 → 坐标有系统误差（与 §9.3 清单法的近似同源，量级可能更大） |
| **多解性（EKMA 歧义）** | 同一 (NO₂, O₃) 可能对应"低 VOC 高 NOx"与"高 VOC 低 NOx"——**必须报告候选范围**，不报单点。加 HCHO 第三维可压缩（见 §9.4.1 `--with-hcho`） |
| **覆盖范围** | 观测超出模拟状态包络（NO₂mean / O3max 的 min–max，脚本会逐维检查并警告）→ 最近点仍是边界值，**外推无意义** |

#### 9.4.4 三种定位法对比（本实验并存，按数据可得性选）

| 方法 | 需要的数据 | 精度 | 章节 |
|---|---|---|---|
| **A. 清单比值（首选）** | NOx + NMVOC 清单 | 高（轴即清单，无歧义） | §9.3 |
| **B. O₃+NO₂ 反查（本节）** | O₃ + NO₂ 观测（常规可得） | NOx 轴 ±10%、VOC 轴 ±40–50% | §9.4 |
| **C. O₃+NO₂+HCHO 反查** | 再加 HCHO（卫星/站点） | 接近 A（VOC 轴歧义被压缩） | §9.4（`--with-hcho`） |

> **推荐**：有清单走 A；无清单但有 HCHO 走 C；只有 O₃/NO₂ 走 B 并**如实报告 VOC 轴的不确定带**。
> 与 EKMA-FNR 的分工见 §9.3.2 末与 `EKMA-FNR/README` §11.6。
