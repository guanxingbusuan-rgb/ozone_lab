# ozone_lab —— 0-D 箱式大气化学模型 · 三个自完备程序

> 机制 V2.5+ 冻结（111 物种 / 203 反应 / 127 命名速率 / 光解槽号至 45，实占 25 槽·26 条光解反应）。
> 仓库由 **3 个程序** + **5 个共享层** 组成；每个程序自带输入快照与结果存档，
> 只读自己的 `config_override/`，产物只落自己的 `output/`。

---

## 1. 目录

```
ozone_lab/
├── CASE-analysis/      个例分析（六账 + 反应级归因 + 前体物 RIR）
├── EKMA-FNR/           S2 状态空间 EKMA（初始 NOx × 初始 HCHO）
├── EKMA-emission/      S1 排放 EKMA（E_NOx × E_VOC）
├── mechanism/          机制 .fac（冻结）+ reference/（MCM 基准）+ audit/（静态审计）
├── facsimile_impl/     执行器：fac_v3 / fac_batch_omp + 生成码 gen_v3/
├── modules/            共享物理层（编译期）：runtime_dat_mod / chem_phys_mod
├── inputproc/          cfgbuild.py（配置→runtime 编译器）、sat2surf.py（卫星柱→表面）
├── analysis_lib/       csvio.py（CSV 契约读写）、indicators.py（指标 / O₃ 收支）
├── Makefile            只建两个二进制；mcmcheck / jcheck / clean-build
├── ENVIRONMENT.md      环境、依赖矩阵、工具链版本、输入硬约束
├── check_env.sh        13 项环境自检（只读）
└── README.md
```

---

## 2. 三个程序

| 程序 | 科学问题 | 扰动轴 / 运行 | 默认网格 | 主要产物 |
|---|---|---|---|---|
| **EKMA-emission** | 排放变化如何移动 O₃ 屋脊与体制？ | `E_NOx` × `E_VOC` 组（排放开） | 0–625% / 25% ⇒ 26×26=676；**4 天 / spinup 3** | `ekma_emission_grid.csv`（含 `VNr_group/VNr_nmv`）、`_derived.csv`（PR、`ridge_at_edge`）、`_ridge.csv`、`ekma_emission_S1.png`；**站点反查**：`plot.py --obs-single`（单站）/ `--obs-more`（多站 CSV）`[--with-hcho]` → 政策空间坐标 + 屋脊侧（图面三档标记与 EKMA-FNR 同构） |
| **EKMA-FNR** | 初始 NOx / HCHO 状态如何决定 O₃ 与 FNR 体制？ | `初始 NOx 组`（NO 与 NO₂ 同因子）× `初始 HCHO`；排放关闭 | 0–1000% / 100% ⇒ 11×11（`--pmax/--step` 可调） | `ekma_fnr_grid.csv`（含 `FNR_init`）、`_derived.csv`（PR、`FNR_init_ridge`）、`_ridge.csv`、`ekma_fnr_S2.png`；**站点定位**：`plot.py --obs-single`（单站）/ `--obs-more`（多站 CSV）→ 三观测状态匹配 + 体制判定 |
| **CASE-analysis** | 这一个个例的 O₃ 由哪些过程 / 反应 / 前体物控制？ | 不扫网格：1 个基线 + RIR 扰动 | RIR 默认 −10%、`--window all`（全逐时；夜间 P≤0 自动剔除） | `ipr_ledger.csv`（六账长表）、`irr_top.csv` / `hourly_rxn_top.csv`、`rir_mean.csv` / `rir_diurnal.csv`（`--per-species` 另出 `rir_species_{mean,diurnal}.csv`）、`p_o3no_baseline.csv`、`case_{ipr,irr,rir,o3_budget,rir_fig5}.png` |

每个程序目录结构一致：

```
<PROGRAM>/
├── run.py            # 唯一入口：配置快照 → cfgbuild → 模拟 → 导出全部 CSV（+run.log/meta）
├── plot.py           # 只读 output/ 的 CSV，出 PNG（不导出数据）
├── README.md
├── config_override/  # 完整配置集（11 个顶层文件 + constraints/；CASE 另有 obs/）——唯一输入
└── output/           # 存档：原始输出 + 分析 CSV + 日志/元数据 + 图
```

---

## 3. 快速开始

```bash
bash check_env.sh                 # 环境自检（13 项，当前 13 OK / 0 MISS）
make all                          # 需要时构建两个二进制（源码未变则跳过）

python CASE-analysis/run.py       # 个例：六账 + IRR + RIR；--no-rir 可跳过 RIR
python CASE-analysis/plot.py

python EKMA-FNR/run.py            # 121 格；--from <raw.csv> 可只重导出
python EKMA-FNR/plot.py

python EKMA-emission/run.py       # 676 格 / 4 天；同样支持 --from
python EKMA-emission/plot.py
```

换个例 / 换输入：直接编辑对应程序的 `config_override/`（气象序列、初值、排放、沉降、光解、情景），
再跑 `run.py`；**不回退任何全局配置**。

---

## 4. 共享层

| 目录 | 作用 | 读取时机 |
|---|---|---|
| `mechanism/mechanism_3_full_chain.fac` | 机制（冻结）；CASE 还用它读反应式 | 运行期 |
| `mechanism/reference/mcm_v3.3.1_export.fac` | MCM 审计基准 | 审计 |
| `facsimile_impl/fac_v3` | 单箱深度执行器（写 post-spinup 全部步的 conc/ipr/irr） | 运行期（CASE） |
| `facsimile_impl/fac_batch_omp` | 批量执行器（OpenMP；末 24 步聚合） | 运行期（两个 EKMA） |
| `facsimile_impl/gen_v3/` | 生成码 + `oz.reac/oz.prod` 计量表 | 运行期（CASE 读计量表）/ 编译期 |
| `inputproc/cfgbuild.py` | 配置 → `runtime.dat` 编译器（唯一） | 运行期 |
| `analysis_lib/csvio.py` | CSV 读写 + gfortran 三位指数兜底 | 运行期 |
| `analysis_lib/indicators.py` | 机制感知指标；`p_o3_gd`（O₃ 收支 Eq4–6） | 运行期（CASE） |
| `modules/*.f90` | 物理层与 runtime 消费端（**静态编入二进制**） | **仅编译期** |

---

## 5. 构建与机制审计

```bash
make all                          # fac_v3 + fac_batch_omp
make clean-build                  # 删二进制/.mod/生成码（不动各程序 output/）
make mcmcheck                     # 48 条速率 vs MCM：最大偏差 8.9%
make jcheck                       # 光解槽位语义：0 错位
python mechanism/audit/element_audit.py mechanism/mechanism_3_full_chain.fac   # 203/203
```

工具链版本、运行期/编译期依赖矩阵、输入格式硬约束（每序列 ≤128 点、序列池 `NSER_MAX=40`、
气象禁 `# unit:`、BLH 必须 >0、SPINUP 校验等）见 [ENVIRONMENT.md](ENVIRONMENT.md)。

---

## 6. 数据契约（三个程序共用）

- CSV 一律由各程序 `run.py` 写出，`plot.py` **只读 CSV、只写 PNG**（站点定位明细也只打印终端）；
  无 JSON / manifest。
- **统计窗**：`O3max/O3mean/HCHOmax/NO2mean` = batch 的**末 24 个输出步**（dtout=3600 时＝最后 1 天），
  不是 8h 日最大。
- **`h` 列语义**：`fac_v3` 三件套（`_conc/_irr/_ipr`）的首列 `h` 是**输出步序号**，不是小时；
  当地小时 = `h × dtout/3600`。CASE 派生表因此带 `hlocal`（本窗起点当地小时）。
- **两种驱动器的账本单位不同**：`fac_v3` 的 `_ipr/_irr` 是**窗口积分** [molec/cm³]；
  `fac_batch_omp` 的 `_ipr/_irr` 侧车是**末 24 步平均倾向** [molec/cm³/s]——相差 dtout 倍，
  跨程序直接对比前必须换算（两侧 `run_meta.txt` 都写了这条口径）。
- `fac_v3` 的 `conc/irr/ipr` 保留**全部 post-spinup 步**；`fac_batch_omp` 为**每案例一行**的聚合。
- `PR_*` 与 `FNR_init` 中的 “0%” 一律指**基线**（不是零排放 / 零浓度）；`FNR_init` 的
  H₀/N₀ 从该程序 `output/fnr_ic_base.csv` 读实际基线。

---

## 7. 已退役（一句话）

通用实验层（`experiments/ configs/ output/ runtime_v3.dat tests/ mcm_city/ docs/ environment/`）
于 2026-10-29 退役，历史内容已删除；需要时从 git HEAD 取回。
