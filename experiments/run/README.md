# experiments/run —— 实验运行脚本（类 CMAQ run.cctm.csh）

> **定位**：**只负责"跑模拟 + 落 CSV 输出"**，不含任何绘图、不含伴随反演。
> 十个实验：A 情景沙盘 · B 响应面 · C 屋脊对比 · D 过程诊断 · E 指示因子 ·
> F OSSE · G 逐小时 · H 敏感性 · I 情景矩阵 · **M 气象序列驱动**。
> 类比 CMAQ：`run.cctm.csh` 只跑化学输运并写输出文件，绘图是另外的 NCL/Python 脚本。
> 本目录对应前者；绘图在 `figures/`。

---

## 1. 与 output/ 和 figures/ 的关系（重要）

```
configs/ ──cfgbuild──> runtime_v3.dat ──fac_v3/fac_batch_omp/mcm_run──> output/run_X/*.csv
                                                                              │
                                                                              ▼
                                                                    figures/*.py （绘图，读 CSV）
```

| 目录 | 角色 | 内容 |
|---|---|---|
| **`output/`** | **实验产物（CSV 数据）** | `run_A/` … `run_M/`，每个实验一个子目录；全部是数值表，可直接 pandas/Excel 打开 |
| **`figures/`** | **可视化** | `fig_response_v3.py` 等**读 output/ 的 CSV** 才是工程图；`gen_all.py`/`gen_cartoon.py` 生成的是**机制示意图**（不依赖输出数据） |

**关键澄清**：
- `figures/gen_all.py`、`gen_cartoon.py` → **机制示意图**（教学图，纯手绘式，不读任何模拟结果）
- `figures/fig_response_v3.py`、`fig_indicators_ridge.py` → **工程图**（读 `output/batch_v3_indicators.csv`）

---

## 2. 输出文件会不会被自动删除？（不会）

| 操作 | 对 `output/` 的影响 |
|---|---|
| 跑任何 `run_*.py` | **只写自己那个子目录**（`output/run_A/` 等），幂等覆盖同名文件 |
| 跑 `make test` | 写 `output/` 根下的 `cmp_*`/`smoke_*`/`ft_*`，**不动 `run_*/`** |
| **`make clean`** | ⚠ **会 `rm -rf output`** —— 这是唯一会删你结果的命令 |
| 跑 `run_all.py` | 各实验覆盖各自目录，不动别人的 |

**建议**：
```bash
# 想保留结果时，别用 make clean，改用：
make clean-build      # 只删二进制与生成码（若不存在该目标见下）
# 或把结果放到仓库外：
OZONE_OUT=$HOME/ozone_results python experiments/run/run_all.py
```
> `make clean` 是"恢复到干净仓库"用的，**做实验期间不要跑它**。

---

## 3. 十个实验脚本（A–I + M）

| 脚本 | 实验 | 输入（假设已备好） | 输出（CSV） | 耗时 |
|---|---|---|---|---|
| `run_A_scenarios.py` | 六情景单箱政策沙盘 | runtime_v3.dat | `<SCEN>_{conc,irr,ipr}.csv` + `summary.csv` | ~3 s |
| `run_B_surface.py` | O₃–NOx–VOC 响应面（4D） | runtime_v3.dat | `batch.csv`（24 列）+ **17 个 `indicator_matrix_*.csv`** + `slices.csv` | ~15 s |
| `run_C_ridge.py` | 屋脊跨路线/跨机制对比 | runtime + mcm_run | `ridge_batch.csv` `ridge_mcm.csv` `ridge_compare.csv` | ~16 s |
| `run_D_ipr.py` | **IPR 六账过程账本** | runtime | `ipr_ledger.csv`（长表：逐小时 × 111 物种 × 6 账，单文件） | ~1 s |
| `run_E_indicators.py` | 指示因子（12 指标） | 任意已有前缀 | `<tag>_indicators.csv` `_timeseries.csv` `compare.csv` | <1 s |
| `run_F_osse.py` | 观测系统模拟实验 | runtime | `obs_samples.csv` `skill_rmse.csv` `skill_delta.csv` + τ 扫描输出 | ~4 s |
| `run_D_irr.py` | **IRR 反应级归因**（物种参数化，默认 O3） | runtime + 计量表 | `irr_top.csv` `hourly_rxn_top.csv` | ~1 s |
| `run_G_consist.py` | **IPR↔IRR 双路径互证**（物种参数化） | runtime | `consistency.csv` | ~1 s |
| `run_H_sensitivity.py` | OAT 敏感性（±10%，18 次） | runtime | `jacobian.csv` `sensitivity.csv` `oat_runs.csv` | ~11 s |
| `run_I_scenario_matrix.py` | 情景 × 路线矩阵 | runtime + mcm_run | `matrix_o3max.csv` `matrix_delta.csv` `direction.csv` | ~4 s |
| `run_M_met_series.py` | **气象观测序列驱动 + 单变量归因** | runtime + `met_*_demo.txt` | `summary.csv` `attribution.csv` `o3_timeseries.csv` | ~6 s |

**总驱动**：`run_all.py`（全跑 ~78 s，顺序 A→B→D→R→G→E→F→H→M→I→C）

**D/R/G 的拆分逻辑**（2026-10 起按依赖解耦）：
`run_D_ipr` 只读 `_ipr.csv` 一个文件（过程账本，长表逐小时×全物种；占比%与全程总量可由表内重算，闭合校验方法见旧版 `ipr_closure_check.py`）；
`run_D_irr` 只读 `_irr.csv`+`gen_v3/oz.reac|oz.prod` 计量表（反应归因，物种可为机制中任意一个）；
`run_G_consist` 是唯一同时读两份 CSV 的脚本——"IPR vs IRR 一致性"必须两边都在，这是设计使然。
`run_D_ipr`/`run_D_irr` 均**不写 manifest.json**，自跑 BAU 的原始模拟输出落 /tmp 读完即删——产物目录只留分析 CSV（前者 1 个，后者 2 个）。
`run_G_consist` 是唯一同时读两份 CSV 的脚本——"IPR vs IRR 一致性"必须两边都在，这是设计使然。

---

## 3.5 气象观测序列（实验 M）

**能像站点约束那样喂气象变量吗？能。** 6 个气象量都支持 `file:` 序列：

| 变量 | 序列样例 | 原生单位 | 物理入口 |
|---|---|---|---|
| `TEMP` | `met_TEMP_demo.txt` | K | Arrhenius 速率 + 数密度 M |
| `PRES` | `met_PRES_demo.txt` | hPa | 数密度 M（三体反应） |
| `RH` | `met_RH_demo.txt` | % | H₂O → O(¹D)+H₂O 供 OH |
| `BLH` | `met_BLH_demo.txt` | m | 体积/稀释/沉降/膨胀 |
| `RAIN` | `met_RAIN_demo.txt` | mm/h | 湿沉降 k_wet |
| `JFAC` | `met_JFAC_demo.txt` | 1 | 全部 J 的乘子 |

每条 121 点（0–120 h 逐时），样例含一段冷槽阵雨 + 一天云盖，方便看到效果。

```bash
# 改 configs/meteo.config（末尾已备好注释块）
# 然后：
make runtime_v3.dat
python experiments/run/run_M_met_series.py
```

**对照设计**：`M0_const`（六变量全取均值）→ `M1_all`（全序列）→ `M2..M7_*_only`
（单变量）。实测主导项是 **BLH（−45 ppb）** 与 **TEMP（+16.8 ppb）**，且二者
**不可加和**（单变量和 −27 vs 实际 −33）——`attribution.csv` 末列出现 −137%/+50.6%
正是"效应不可加"的标志，不是 bug。

**硬规则**：气象序列里**不写 `# unit:` 行**（那是浓度换算指令，写了会被 cfgbuild
直接拒绝；旧版会静默 ×2.46e10）。6 个气象序列会占掉序列池 8 个名额中的 6 个。

---

## 4. 用法

```bash
cd ozone_lab

# ── 前置（一次性）──
make all && make runtime_v3.dat      # 构建二进制 + 编译配置
# mcm 相关实验还需: 见 §5

# ── 全部实验（A–I + M，除反演）──
python experiments/run/run_all.py

# ── 单个实验 ──
python experiments/run/run_A_scenarios.py              # 六情景全跑
python experiments/run/run_A_scenarios.py BAU LNNOX    # 指定情景
python experiments/run/run_B_surface.py 13             # 13×13 网格
python experiments/run/run_B_surface.py 9 0.6,0.8,1.0  # 加光解维
python experiments/run/run_D_ipr.py output/run_A/BAU       # IPR 账本，复用已有输出
python experiments/run/run_D_irr.py HCHO                   # IRR 归因，指定物种
python experiments/run/run_G_consist.py PAN output/run_A/BAU  # 互证，指定物种
python experiments/run/run_H_sensitivity.py 5          # ±5% 扰动
python experiments/run/run_M_met_series.py             # 气象序列对照（10 个 run）
python experiments/run/run_M_met_series.py --const-only  # 只跑常数基准

# ── 换输出目录（不污染仓库）──
OZONE_OUT=$HOME/ozone_results python experiments/run/run_all.py

# ── 列出实验清单 ──
python experiments/run/run_all.py --list
```

---

## 5. mcm_run 的前置构建（C / I 需要）

`make clean` 会删除 mcm 的生成码，需重建：

```bash
# ① 重新生成 MCM 机制码（mechc.py）
cd facsimile_impl && $PY mechc.py ../mcm_city/input/mechanism_mcm_ch4.fac ../mcm_city/gen

# ② 编译 mcm_run
cd ../mcm_city
$FC -O2 -ffree-line-length-none -J gen -I gen gen/oz_mech.f90 mcm_driver.f90 -o mcm_run \
    -L$CVODELIB -Wl,-rpath,$CVODELIB \
    -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial
```
（`$PY`/`$FC`/`$CVODELIB` 见 Makefile 顶部；`make all` 不含 mcm_run）

---

## 6. 输出契约（每个 CSV 可直接读）

| 文件类型 | 格式 | 典型用法 |
|---|---|---|
| `*_conc.csv` | 首列 `h`（步序号=小时@dtout=3600）+ 111 物种列（molec/cm³，SPINUP 后全部逐时） | `pandas.read_csv` → 除以 2.46e10 得 ppb |
| `*_irr.csv` | `h` + `r1..r131`（小时窗**积分**，非速率） | 反应级归因 |
| `*_ipr.csv` | `h` + 74 物种 × 6 账（`P L E D B X`，小时窗积分） | 过程级归因 |
| `indicator_matrix_*.csv` | 行=fnox、列=fvoc 的二维矩阵 | 直接 `plt.contourf` |
| `summary.csv` | 汇总表（每行一个情景/参数） | 直接看/画柱状图 |
| `manifest.json` | 运行元数据（时间/参数/产物清单） | 溯源 |

**单位提醒**：
- `conc` 是 molec/cm³ → 除 **2.46e10** 得 ppb
- `irr`/`ipr` 是**窗口积分**（molec/cm³，不是速率）——不要再对它们做 `trapz`

---

## 7. 与旧脚本的关系（2026-10-27 决议：旧入口已删除）

旧入口 `ridge5/scen5/ipr_irr5/hourly_pa/sens_adj/obs_osse` 已全部移除——
经逐项核查，其功能由 `run/` 全覆盖（跨求解器/五路线历史口径随路线退役一并废止；
原始脚本可查 git 历史）。诊断职能在 `experiments/post/`，绘图在 `experiments/viz/`。
**新脚本（run/）的设计目标**：统一输出目录、统一 CSV 契约、无绘图、幂等、可批量。
