# ozone_lab 三路线实验操作手册

> 手把手级实验详解（每实验 目的/命令/产物/解读/坑）见 [EXPERIMENTS_GUIDE.md](EXPERIMENTS_GUIDE.md)。

> **适用版本**: 2026-10-27（mechanism_3 **111sp/203rx V2.5+ 冻结版** · **3 可执行体** · IPR 精确累积版；
>  历史口径"五路线"中 fac_run/kpp_run 已退役——见路线表与实验退役注记）
> **前置**: `~/soft/atmenv/cvode`（CVODE 2.7.0）、gfortran 13.2、
> Python venv（numpy/matplotlib/netCDF4）。全部为第三方工具链，无 Day xx 依赖
> （KPP/bison/flex 已退役不属前置；版本明细见 [ENVIRONMENT.md](../environment/ENVIRONMENT.md)）。
> **路线一览**:
>
> | 可执行体 | 机制 | 求解内核 | 情景空间 |
> |---|---|---|---|
> | `fac_v3` | 111sp/203rx 工程 (V2.5+) | CVODE-BDF+DENSE (fcmix) | runtime.dat 情景表 |
> | `fac_batch_omp` | 111sp/203rx 工程 (V2.5+) | CVODE-BDF+DENSE（C API 直连，OpenMP） | fnox×fvoc×jfac×ΔT 网格 |
> | `mcm_run` | 29sp 真实 MCM CH₄ | CVODE-BDF+DENSE | 排放缩放 fnox×fvoc |
>
> **实验总目录**（详细操作见 §3–§13）:
> | 实验 | 内容 | 主工具 | 路线覆盖 |
> |---|---|---|---|
> | A | 单箱基线与政策情景 | fac_v3 | v3 |
> | B | 响应面与指示因子屋脊 | fac_batch_omp | v3 双胞胎 |
> | C | 屋脊路线对比（历史五路线口径） | ridge5.py | mcm/batch（fac_run/kpp 已退役） |
> | D | IPR/IRR 过程诊断 | ipr_irr5.py + closure | 现役三路线 |
> | E | 指示因子工具箱 | indicators.py | 任意 (CSV+FORMULA) |
> | F | 观测约束 OSSE | run_F_osse.py | v3 双胞胎 |
> | G | 逐小时过程分析 | hourly_pa.py | v3 (精确积分) |
> | H | 排放敏感性（OAT） | run_H_sensitivity.py | v3（差分法；~~伴随/反演~~ 已移除：0-D 无真伴随，SVD 不可辨识） |
> | I | 情景矩阵（历史五路线口径） | scen5.py | 现役三路线 |
>
> **脚本位置**: `experiments/run/`（10 个实跑脚本）；`experiments/{post,viz}/`（诊断/绘图）
> 运行器）；`analysis_lib/`（诊断库+indicators 引擎）；`audit/`（元素审计）。
> 全部用 venv 解释器手跑或 `make` 目标（test/batch/figures/nc）。
>
> **工作目录约定**: 除非注明，命令均在 `ozone_lab/` 根目录执行；Fortran 可执行体的
> 相对路径参数（`gen_v3`、`../runtime_v3.dat`）要求 cwd 在其实现子目录。
> **Python 解释器**: `experiments/` 脚本依赖 numpy/matplotlib，须用 venv 解释器
> （`$(CONDA_ENV)/bin/python`，即 Makefile 的 `$(PY)`，默认 `$HOME/miniconda3/envs/venv`）；
> 系统 `python3` 缺科学栈——`make` 目标已用 `$(PY)`，直接手跑时请注意。

---

## 架构总览：目录 × 模块地图（全目录详解）

### 七层数据流（每个文件在管线中的位置）

```
┌──────── 表达层 ──── figures/ (机制图集+响应面图)  OPERATIONS_MANUAL/README (文档)
├──────── 实验层 ──── experiments/run/ (10 个实跑脚本, 含体制诊断)
├──────── 诊断层 ──── analysis_lib/ (indicators.py 引擎 + diagnostics.py 账本库)
├──────── 求解器层 ── facsimile_impl/ (fac_v3 + fac_batch_omp)   mcm_city/ (mcm_run)
├──────── 运行时层 ── modules/ (runtime_dat_mod + chem_phys_mod)
├──────── 编译层 ──── inputproc/cfgbuild.py  facsimile_impl/mechc.py
├──────── 合同层 ──── docs/templates/ (6 类外数据模板)   configs/ (12 类生效配置 + constraints/)
└──────── 事实源 ──── mechanism/ (mechanism_3_full_chain.fac) + audit/ (元素守恒审计)
```

### 逐目录详解（16 目录 + 根）

| 目录 | 角色 | 关键文件 → 功能 |
|---|---|---|
| **`mechanism/`** | **机制唯一事实源**（1 个文件） | `mechanism_3_full_chain.fac`：**V2.5+ 冻结版 111sp/203rx**，十段链条骨架（①无机核心+HONO 晨爆 ②NO₃/N₂O₅ 夜化学 Troe ③CH₄+PAN 储库 ④ISOP 二代三分支 ⑤DMS 硫链 ⑥卤素闭环 ⑦集总类 ⑧支撑 ⑨硝酸酯 ⑩单萜+V2.4 扩展）+ V2.5 补正（光解槽位对齐 MCM 官方、GLY 三分支、O₃+HOx、IOLE/ETH-NO₃ 专属速率）；自带 FORMULA 段；**规格书见 mechanism/README.md**（逐反应速率/配平方法论/对比定位）|
| **`src/`** | **两个编译器**（Python，全部输入唯一入口） | `cfgbuild.py`：12 类 config+constraints → runtime_v3.dat，编译期校验（物种存在性/单位指令/单调性/≤128 点/情景表达式）；`mechc.py`：机制编译器（.fac → gen_v3 生成码）；（`mechc_c.py` 已于 2026-10 删除：零消费者） |
| **`modules/`** | **Fortran 共享模块层**（三路线共用） | `runtime_dat_mod.f90`：runtime.dat 唯一消费端（12 段解析+气象插值+SZA→J+BOX 日型+ser_interp+build_groups 组映射+scen_build）；`chem_phys_mod.f90`：**物理四账唯一实现**（排放①/干+湿沉②/BC 序列观测同化③/BLH 膨胀④）+ `tune_atol_by_name`（ATOL 单源）+ `eval_met_dt`/`eval_J_jf`；`rad_photol.f90`：B2 在线光解 demo |
| **`facsimile_impl/`** | **工程主线**（两条活路线） | `gen_v3/`：生成码（oz_mech.f90 + oz.reac/prod 计量表，fac_v3/batch 共享）；`oz_fac_v3.f90`→`fac_v3`：单箱深度（情景名/观测同化/六账 IPR 逐时导出）；`fac_batch_omp.f90`+`cvode_c_mod.f90`→`fac_batch_omp`：批量广度（C-ABI 直连+OpenMP、NOy 17 物种权重表+运行期氮守恒不变量、24 列指示因子）；（原列于此的 `oz_audit.c` C 证人已删除：零调用方，守恒校验由 element_audit/氮守恒不变量/run_D 闭合承担） |
| **`mcm_city/`** | **真实机制对照**（刻意冻结，游离 make 主链外，**零代码耦合**） | `mcm_driver.f90`→`mcm_run`：MCM v3.3.1 CH₄ 子集 29sp/71rx 城市模拟（CLI: days/pref/fnox/fvoc，排放/IC/沉降全硬编码）；`gen/` 独立生成码。**详见 [mcm_city/README.md](../mcm_city/README.md)** |
| **`analysis_lib/`** | **诊断层**（路线无关） | `indicators.py`：route_indicators/route_timeseries——从任意 conc/irr CSV+机制 FORMULA 计算 12 指标（NOy 自动判定，新 N 物种零改动接入）；`diagnostics.py`：经典库（load/indicators/irr_integrate/families/chem_ledger/report） |
| **`experiments/`** | **实验层** | `run/`：**10 个实跑脚本**（A 情景/B 响应面/C 屋脊/D 过程/E 指标/F OSSE/G 逐时/H 敏感性/I 情景矩阵/M 气象序列），只落 CSV 不绘图；`post/`（诊断）`viz/`（绘图） |
| **`docs/templates/`** | **六类外数据格式合同**（零代码直读，人工落位规范） | 地面 O₃/AQI/VOC 站点 ×3 → bc.config 张弛；排放清单 → emissions.config；卫星柱（**语义例外**：离线转表面 ppb）；**气象序列 → meteo.config** |
| **`configs/`** | **工程生效配置**（12 类，唯一输入入口；原 input_v3，2026-10-27 更名） | 12 config + constraints/(气象/BC 序列) → cfgbuild → runtime_v3.dat。旧 `input/` 与 15sp 教学档均已彻底删除 |
| **`audit/`** | **机制质检** | `element_audit.py`：逐反应 C/H/N/O/S/Cl/Br 配平（203/203；储库方向限定防"非物理闭合掩盖真错"）；`mcm_compare.py`：MCM v3.3.1 系数比对 |
| **`tests/`** | **回归防线** | functional_test.py：**62 断言（12 族 FT-\*）**，含 FT0 生成码一致性（物种+反应双守卫）与 4 金值哨兵；smoke_test.sh |
| **`figures/`** | **表达层** | SCI 版 9 图（A-I，Times+宋体 330 DPI PNG+PDF）+ comic 漫画版 9 图（场景化图标）+ 工程图 2 张（response_surface_v3/indicators_vs_ridge）+ 2 绘图脚本 |
| **`output/`** | 全部 CSV 产物 | batch 网格×2、expA 六情景、ridge5/scen5/osse/sens 子目录 |
| **`environment/`** | **环境资料** | `ENVIRONMENT.md`（依赖清单/混合语言分工/从零复现）+ `check_env.sh`（14 项只读体检，自动锚定工程根） |
| **`docs/`** | **操作文档** | `OPERATIONS_MANUAL.md`（本手册）+ `EXPERIMENTS_GUIDE.md`（手把手实验指南） |
| **根目录** | 编排+文档 | `Makefile`（all/v3/batch/test/figures/nc）、`README.md`（项目说明）、`UNITS.md`（单位参考）、`runtime_v3.dat`（编译产物）；文档与环境资料分别在 `docs/` 与 `environment/` |

### 三路线 × 目录映射

| 路线 | 机制 | 生成码 | 物理层 | 输入 | 排放可配 |
|---|---|---|---|---|---|
| ① fac_v3 | gen_v3 (111sp/203rx V2.5+) | mechc.py | chem_phys_mod | runtime_v3.dat 全通路 | ✓ 25 物种已配（任一可加） |
| ② fac_batch_omp | gen_v3（**与①同源**） | mechc.py | chem_phys_mod | runtime_v3.dat（共享） | ✓ 同①基线 × 网格因子 |
| ③ mcm_run | mcm_city/gen (29sp MCM 子集) | mechc.py（独立子树） | mcm_driver 内嵌（冻结） | CLI 四参 | ✗ 刻意冻结 |

### 墓碑（六轮演化的已移除项，均有守卫或文档标注）

fac_run（15sp 教学，ridge5 守卫）· kpp_run/kpp_impl（EKMA IC 空间）· EBI ·
fac_expanded · fac_spgmr · kpp_batch · `rate_sensitivity.py`/`species_constraint.py`
（B3/B4 载体随 15sp 层退役，灵敏度对照由 sens_adj.py 承担——该脚本已于 2026-10-27 删除，OAT 由 run_H 承担）· 15sp 核心档
（15sp 教学机制已删除；`mechanism/` 现仅含 V2.5+ 冻结版 + `reference/` 下的 MCM 官方导出）。

---

## 0. 构建与一键验证

> 环境依赖与复现步骤见 [ENVIRONMENT.md](../environment/ENVIRONMENT.md)；环境体检 `sh environment/check_env.sh`。

```bash
make all          # 两个工程可执行体 (fac_v3/fac_batch_omp)
make mcm          # MCM 对照路线 (mcm_run, 29sp 独立子树) — R15 补建
make test         # 冒烟 1/1 → 功能 67/67 (13 族 FT-*, 含金值哨兵/CFG/DAYS/SPINUP/PHOT/MECHK) → 元素审计 203/203
make run DAYS=3 SCEN=SCEN_BAU PREF=../output/d3   # 单情景一步跑（天数控）
make clean-build  # ⚠ 安全清理: 只删二进制+生成码, output/ 实验结果保留
make clean-output # ⚠ 安全清理: 只删实验结果(交互确认), 二进制保留
make clean        # 🚫 连 output/ 实验结果一并删除 —— 只用于"恢复干净仓库", 做实验期间禁用
```

> **🚨 实验纪律（机制 V2.5+ 封版时代）——三条，违者丢数据/烧时间**：
>
> 1. **机制只读**：`mechanism/` 已封版，实验永远不改 `.fac`。日常循环 =
>    改 `configs/` 输入 → `make runtime_v3.dat`（0.37 s，幂等）→ 直接跑已编译的
>    `fac_v3` / `fac_batch_omp`。只有改 `.fac`/Fortran 才 `make v3`（~8 s）——
>    那是开 V2.6 的动作，不是实验动作。
> 2. **做实验期间严禁 `make clean`**：它会 `rm -rf output/`（**你的全部实验结果**），
>    并删除二进制与 `gen_v3/` 生成码（fac_v3 运行时必需的 oz.reac/oz.prod）。
>    需要清理只允许两个安全变体：`make clean-build`（删二进制+生成码，**output/ 保留**；
>    重建 = `make all && make runtime_v3.dat`，mcm 实验再加 `make mcm`）、
>    `make clean-output`（只删实验结果，交互确认，保留二进制）。
> 3. **runtime_v3.dat 只经 make 生成**：不手改、不用旧文件冒充
>    （FT-CFG 逐字节比对 + Makefile FORCE 强制重建双防线，防同秒 mtime 陷阱）；
>    输入变体走 `experiments/run/_lib.py` 的 `clone_input`（复制到独立目录编译，不碰主 configs/）。

**模拟天数**：`DAYS=` / `--days` **临时覆盖，不改 `model.parameters` 文件**。
实测收敛（末 24h O₃）：1 天 42.69 / 3 天 58.59 / 5 天 59.21 / 7 天 59.23 ppb
→ 建议 ≥3 天。详见 [GUIDE §3⑫-1](EXPERIMENTS_GUIDE.md)。

`make test` 全绿是任何实验的**前置条件**。审计 203/203 由 `audit/element_audit.py`
执行（机制自带 FORMULA 段为配平事实源）。

## 0.5 路线分工与编译边界（工程实践三铁律）

**双胞胎 vs 体外参照**：fac_v3（单箱·深度）与 fac_batch_omp（批量·广度）是一个
产品的两种跑法——机制生成码（gen_v3）、物理层（chem_phys_mod）、输入管线
（runtime.dat）三样全共享；改一次机制/物理，两条路线自动同步，且互为回归
（fac_v3 单点应能被 batch 网格对应格点复现，functional_test 内置该互证）。

**显微镜 vs 望远镜：同一内核、两个窗口**——两者跑的是完全同一个模拟，差别全部在
接口层，不在物理或数值：

| | ① fac_v3 = 显微镜 | ② fac_batch_omp = 望远镜 |
|---|---|---|
| 内核四样 | 全同：机制(gen_v3 111sp/203rx) · 物理(chem_phys_mod 四账) · 数值(CVODE-BDF+DENSE) · 输入基线(runtime_v3.dat) | 同左 |
| 一次运行 | 1 个情景 | 81+ 个格点 |
| 输出深度 | 全深度：逐时 111 物种 + 203 IRR + 六本账逐时 IPR | 压缩：24 列指示因子的末 24h 均值 |
| 回答的问题 | **"为什么"**——打开机器看每个时刻每个过程 | **"在哪里转折"**——曲面上屋脊与体制翻转 |
| 做不到的事 | 扫不出曲面（一次一情景，扫 81 个=调 81 次） | 过程分析（无逐时物种/反应账，只有均值） |

**接力工作流（标准姿势）**：广度找现象 → 深度查原因 → 广度验稳健。
⚠ 精确性注记：fac_v3 与 batch 在同格点**不是逐位相同**——batch 有 jfac/ΔT 覆写路径
与 IPR 均值窗差异，故互证是 **±10~15% 容差的一致性核对**（functional_test 内置）；
"逐位一致"仅指 batch 自身的并行/串行对照。
mcm_run 用完全独立来源的化学（MCM v3.3.1 完整无机组）跑同类情景——若核心
循环写错，它会系统性偏离，因此定位为**机制正确性体外参照**，不参与生产。
确认：它**不在 V2.3 上、也不该在**——29sp 罐头的价值恰在"独立来源+刻意冻结"
（独立 gen、独立子树、游离 make 主链外）；与 ①② 的对照是**跨机制定性对照**
（无机块行为是否系统性一致），非同机制定量复现；下一次被使用的时机 =
mechanism_3 扩链之后拉一下，确认无机块未被新链条带偏（§0.6 第四道关卡）。

**每次实验要重新编译吗？——不用，除非改了机制本身**（实测耗时）：

| 你改了什么 | 要做什么 | 实测耗时 |
|---|---|---|
| 排放/情景/BC/气象/约束 txt/IC/参数（日常 95% 改动） | `make runtime_v3.dat` | **0.37 s** |
| 什么都不改，只换情景名 | 直接 `./fac_v3 SCEN_X …` | **0.19 s** |
| 机制 .fac（加删物种/反应/改速率常数） | `make v3`（mechc 重新生成+gfortran） | ~8 s |
| Fortran 逻辑（物理层/组映射/诊断） | `make v3` | ~8 s |

**三铁律**：
1. **数据 ≠ 代码**——改"数值/形状/用哪个文件"只走 runtime.dat 数据通路；
   改"有哪些反应/物理怎么算"才碰编译器，边界画在 .fac 文件上。
2. **make 替你判断**——依赖图按时间戳：配置改了代码没变时 `make v3` 输出
   `Nothing to be done`（连编译都不触发）。
3. **batch 一次编译 81 案例**——二进制内部循环跑完网格，不存在逐格点编译。

## 0.6 路线选择速查（日常模拟 = ①②；③ 是保险丝）

| 你要做的事 | 用哪条 | 需要编译吗 |
|---|---|---|
| 跑一个情景、改排放/气象、观测同化、看六本账 | ① fac_v3 | 否（只 `make runtime_v3.dat`） |
| 扫描"排放变 ×0.3~×7 会怎样"、定位屋脊与体制 | ② fac_batch_omp | 否 |
| 机制升级后回归对照 / 结果背书 / 反常排查 | ③ mcm_run | 否（罐头案例，CLI 四参） |

①② 是同一产品的两种模式（单箱深度 vs 网格广度），全部生产性工作由它们承担；
只跑单情景时 ① 独立即可，② 用于需要"一整张响应面"时。

**mcm_run 的三个"想起它"时刻**：
1. **机制升级后的回归对照（最重要）**——扩链（如 V2.2/V2.3）之后跑一次：不比 O₃
   数值（机制不同无法直接比），而比无机组行为（NO₃/N₂O₅/HONO 来自 MCM 逐字拷贝的
   独立来源；无机核心被改坏时它会系统性偏离而哨兵可能仍在容差内）。列入机制变更
   checklist 第四道关卡：元素审计 → 氮守恒不变量 → 金值哨兵 → **mcm_run 对照**。
2. **结果可信性背书**——"本机制无机部分与 MCM v3.3.1 独立实现行为一致"这一声明的证据。
3. **反常二分排查**——fac_v3 行为反常时切 mcm 同情景：mcm 也怪 → 共享层
   （物理/输入/求解器）；mcm 正常 → mechanism_3 新增链条。

维护成本≈零（罐头：CLI 四参、无配置、make 主链不碰）——平时零代价，需要时一把就响。

## 0.7 六大分析能力索引（速查 → 详解）

| # | 能力 | 类型 | 一句话 | 命令 | 详解 |
|---|---|---|---|---|---|
| 1 | **观测约束** | 数据同化 | 把观测按 1/τ 拉进模型（nudging） | `bc.config` 填 `file:` + `TAU_MIX` | GUIDE §6.10.1 / §14.1 |
| 2 | **响应面** | 全局扫描 | fnox×fvoc×jfac×ΔT 四维参数空间 | `make batch && make figures` | GUIDE §6 / §6.9 |
| 3 | **情景分析** | 离散场景 | 6 个政策情景一键切换 + 跨路线一致性 | `fac_v3 SCEN_X` + `scen5.py` | GUIDE §5 / §13 / §6.11.1 |
| 4 | **过程分析** | 过程归因 | IPR 六账 + IRR 203 反应（小时窗精确积分） | `hourly_pa.py` + `ipr_closure_check.py` | GUIDE §8 / §11 / §6.11.2 |
| 5 | **敏感性 OAT** | 局部梯度 | 9 排放物种 ±10% → 雅可比（24×9） | `run_H_sensitivity.py` | GUIDE §12 / §6.10.2 |
| ~~6~~ | ~~伴随/反演~~ | 推断 | **已移除（2026-10-27）**：0-D 无真伴随，FD+SVD 实测仅 1/9 模态可辨识） | ~~sens_adj.py~~（已删） | — |
| 7 | **气象序列驱动** | 外部时序 | 6 个气象量（TEMP/PRES/RH/BLH/RAIN/JFAC）逐时观测 → `meteo.config` 的 `file:` | `run_M_met_series.py` | GUIDE §0.3b |

> **三个必须分清的概念对**（详见 GUIDE）：
> - **观测约束**（nudging 机制） ≠ **OSSE**（评价实验）——§6.10.1
> - **OAT**（单参数局部梯度） ≠ **响应面**（多参数全局地图）——§6.10.0
> - **`bc.config` 常数**（背景边界/定环境） ≠ **`bc.config` `file:`**（观测约束/定真值）——§6.10.1
>
> **气象序列 vs 观测约束**：两者**格式相同**（`t[s]` + 值、两列、≤128 点、共享 8 条序列池），
> 但**登记位置和物理身份不同**——`meteo.config` 的 `file:` 是"物理环境"（进方程，
> 不是真值），`bc.config` 的 `file:` 是"观测真值"（经 1/τ 张弛拉向）。且气象文件
> **不允许**写 `# unit:` 浓度换算行（`cfgbuild` 会拒绝），否则 298.15 K 会被静默
> 乘成 7.33e12。

> **模型维度澄清**：本程序是 **0-D 箱体**（状态向量 = 111 物种，无 x/y/z、无平流扩散）。
> "四维"指**参数扫描空间**（排放/光解/温度缩放因子），不是空间网格；`fac_batch_omp`
> 跑的是 **262,144 个互相独立的 0-D 箱**，不是"一个 4D 场"。详见 GUIDE §6.9 与
> [ENVIRONMENT.md §6b.6b](../environment/ENVIRONMENT.md)。

## 1. 机制与输入两条事实源

> **机制规格书**：[`mechanism/README.md`](../mechanism/README.md) —— 111 物种分组表、
> 203 反应逐条速率、45 光解槽位方案、配平方法论、
> 已知省略与 CB6R3/fullchem/MCM 对比定位。

```
mechanism/*.fac (唯一机制事实源)
  ├─ facsimile_impl/mechc.py   → gen*/oz_mech.f90 + oz.reac/prod   (FAC 系求解侧)
  ├─ src/fac2kpp.py            → (KPP 路线 2026-10 退役, kpp_impl/ 已删; 脚本保留)
  └─ inputproc/cfgbuild.py --mech …  → runtime*.dat SPECIES 行            (配置侧)

configs/ (12 类配置) → cfgbuild → runtime_v3.dat → fac_v3 / fac_batch_omp
```

**改机制五步闭环**: 改 `.fac` → `python3 audit/element_audit.py <fac>`（配平审计）
→ `make v3 runtime_v3.dat` → 配置补新物种（IC/排放/沉降）→ `make test`（哨兵拦截漂移）。

## 2. 输入系统与五类外数据模板

12 类配置在 `configs/`；**外部数据**走 `docs/templates/` 六类模板
（详细格式说明内嵌于每个模板文件头）：

> 📋 **全物种输入矩阵**：111 物种 × 5 通道完整清单 →
> [templates/SPECIES_INPUT_MATRIX.md](templates/SPECIES_INPUT_MATRIX.md)

| 模板 | 入口 | 物种/单位 | 物理语义 |
|---|---|---|---|
| `docs/templates/emissions_inventory.template` | `configs/emissions.config` | 机制物种；基线 molec/cm³/s；BOX/FLAT 日型 | 排放账（IPR 第 3 列 E） |
| `docs/templates/surface_o3_station.template` | `bc.config: O3 = file:…` | `# unit: ppb`（cm3/ppb/ppm）+ 两列 t[s]/值 | 向观测张弛（nudging，τ=TAU_MIX） |
| `docs/templates/surface_aqi_station.template` | 宽表拆分转换 → 同上 | AQI 多污染物（HJ 193, μg/m³）→ ppb 换算表 | 同上；PM 不接入；NO₂ 建议弱 τ |
| `docs/templates/surface_voc_station.template` | `bc.config: <物种> = file:…` | 一站一物种一时间序列 + PAMS-57 映射注记 | 同上（短寿命物种为宜） |
| `docs/templates/satellite_column.template` | **`inputproc/sat2surf.py`** 转换 → 站点通路 | TROPOMI 柱浓度 molec/cm² → 表面等效 ppb | 直观式 `柱/BLH/2.46e10`；工具内置两法互证+QA 过滤+抽稀，**引擎不直读柱** |
| **`docs/templates/met_series.template`** | **`meteo.config: TEMP = file:…`** | **原生单位 K/hPa/%/m/mm·h⁻¹/1；禁止 `# unit:` 行** | **物理环境驱动（进方程，非真值）；6 变量可全序列** |

序列文件通用约定：两列 `t[s] 值`、时间严格递增（cfgbuild 编译期校验）、
浓度序列可用 `# unit:` 指令行（气象序列禁止）、缺失时段线性插值端点钳制、
单序列 ≤128 点（Fortran 端 ser_t 上限，超限编译期拒绝）、
**气象+BC 共享 8 条序列池**。

> 📖 **模板的逐步接入操作（原始格式 → 换算 → 落位 → 登记 → 编译 → 验证 → 坑）
> 详见 [EXPERIMENTS_GUIDE.md §14](EXPERIMENTS_GUIDE.md#14-五类可选输入数据文件接入操作templates-逐一详解)**；
> 气象序列的完整对照实验见 [GUIDE §0.3b](EXPERIMENTS_GUIDE.md)；
> 单位换算完整参考见 [UNITS.md](UNITS.md)。

### 2.1 端到端示例：AQI 站点 O₃ → 约束运行（已实测）

```bash
# ① 从 AQI 宽表拆分转换 (μ g/m³ → ppb, ×0.509), 存为序列文件
#    configs/constraints/aqi_O3_BJ-PKU.txt   (照 surface_aqi_station.template 格式)
# ② bc.config 登记约束 (TAU_MIX 决定张弛强度):
#    TAU_MIX = 86400
#    O3 = file:constraints/aqi_O3_BJ-PKU.txt
# ③ 编译配置 (机制唯一事实源直引):
python3 inputproc/cfgbuild.py configs runtime_v3.dat --mech mechanism/mechanism_3_full_chain.fac
#    cfgbuild 编译期校验: 物种存在性/时间单调/单位指令/序列点数≤128/光解槽位完备
# ④ 运行:
cd facsimile_impl && ./fac_v3 SCEN_BAU ../output/aqi_run ../runtime_v3.dat
```

**实测结果**（102→126 μg/m³ 站点序列）：换算 51.9→64.1 ppb → SERIES 编译
→ 运行 O₃ 被张弛拉至 30.9 ppb（无约束 17.9）；IPR 的 B 列即同化倾向分账。
VOCs/卫星模板同路径（卫星先离线转换）；排放清单模板直接对应 emissions.config。

### 2b. 单位、参数与硬编码约定总表（四类全量，cfgbuild 单源 `unit_factor`）

> **完整参考手册（逐子量 + 公式 + 硬编码登记）：[`UNITS.md`](UNITS.md)**。
>
> 用户可输入/可更改的数值分四类。**浓度类单位已全通**（ppb/ppm/cm³ 三态,
> `unit_factor` 单源函数）; **机制速率与物理参数是固定物理约定**（修改 = 改代码重编译）。

#### 表 A · 配置通道（用户数据文件, 单位可声明）

| 通道 | 接受单位 | 换算 | 缺省 | 校验 |
|---|---|---|---|---|
| IC | cm³ / ppb / ppm（每行后缀） | ×2.46e10 / ×2.46e13 | 0 | 未知单位拒绝 |
| BC 序列 file: | 头部 `# unit:` 指令 | 同上 | cm3 | 非单调/超 128 点拒绝 |
| BC 内联常数 | cm³ / ppb / ppm（后缀, V2.3） | 同上 | cm3 | errs 门拒绝 |
| 排放 | molec/cm³/s / **ppb/h**（头部指令, V2.3） | ppb/h ×2.46e10/3600 | molec/cm³/s | 负值/非法日型/振幅拒绝 |
| 气象 7 量 | K / hPa / % / W·m⁻² / m / mm·h⁻¹ | 无需换算 | 各有缺省 | — |
| 干沉降 vd | cm/s | k = vd/BLH | 不列=0 | — |
| 湿沉降 | a: s⁻¹/(mm·h⁻¹)^b | — | 不列=0 | — |
| 光解 JMAP D 值 | s⁻¹（正午） | — | — | 槽位缺失拒绝 |
| model.parameters | 天 / 输出步长[s] / rtol / atol | — | days/dtout≤0 或行数不足 → 拒绝 | 必需；天数用 `--days`/`DAYS=` 覆盖（不改文件），详见 GUIDE §3⑫ |

#### 表 B · .fac 机制参数（半硬编码: 用户可改文本, 改 = 重编译）

| 参数类 | 单位约定 | 修改方式 |
|---|---|---|
| 双分子速率 k | **molec/cm³/s**（Arrhenius: A·exp(−Ea/T)） | 改 .fac 速率表达式 → `make v3` |
| 三体/Troe 速率 | 同上, M = **molec/cm³**（运行时由 PRES/TEMP 计算） | 同上 |
| 光解 l/m/n | l: s⁻¹; J = l·cosᵐ·exp(−n secχ) | photolysis.config（数据层, 免编译） |
| 化学计量系数 | 无量纲（**允许分数**, 如 0.85 CH₃CHO） | 改 .fac → 重编译 |

#### 表 C · Fortran 硬编码参数（改 = 改源码重编译, 已逐项登记）

| 位置 | 内容 | 单位 | 状态 |
|---|---|---|---|
| `mcm_driver.f90:56-60` | E(iNO)=2.0e7 等 5 物种排放 | molec/cm³/s | 冻结（罐头案例） |
| `mcm_driver.f90:179-181` | Y(iO₃)=3.0e11 等 IC | molec/cm³ | 冻结 |
| `fac_batch_omp.f90:69+` | 网格 fs(i)=0.3×1.5^(i-1) | 无量纲因子 | 可改源码 |
| `chem_phys_mod.tune_atol_by_name` | ATOL 按名表（111 物种，未列名者用 model.parameters 第 4 行） | molec/cm³ | **单源**（R9）；RTOL/ATOL 详解见 [GUIDE §3⑬](EXPERIMENTS_GUIDE.md) |
| `runtime_dat_mod.build_groups` | 物种→情景组映射 | — | 单源（含 TERP） |

#### 表 D · 内部计算量（自动, 无用户单位接口）

M/O₂（理想气体, 实际 PRES/TEMP）· H₂O（Magnus 公式 × RH）· SZA（纬度+DOY+时角）·
dBLH（±600 s 中心差分）· 光解 J（l·cosᵐ·exp(−n secχ) × JFAC）。

#### ⚠ PPB 换算因子约定（唯一需要知晓的口径差）

`PPB = 2.46e10` 对应 **1013 hPa / 298.15 K** 的数密度; 运行时 M 按实际 T/P 计算。
非标准 T/P 下, "以 ppb 写入的 IC/BC" 在模型内的实际混合比有 **~1–2%/5 K 偏差**
（HOT 情景 +5 K ≈ +1.1%）——MCM/AtChem2 同款约定, 化学一致性不受影响
（速率的 M 依赖用真实密度）, 仅"ppb 标签"为近似。`.fac` 速率同理。

## 3. 实验 A：单箱基线与政策情景（fac_v3）

```bash
cd facsimile_impl
./fac_v3 SCEN_BAU   ../output/a1_bau  ../runtime_v3.dat    # 基线
./fac_v3 SCEN_LNNOX ../output/a2_ln   ../runtime_v3.dat    # NOx 减 70%
./fac_v3 SCEN_ISOPX ../output/a2_isop ../runtime_v3.dat    # NOx×0.5 VOC×2
# 全情景清单: grep "^SCEN" ../runtime_v3.dat
```

**输出**：`<pref>_conc.csv`（h + 111 物种全列，**spin-up 后全部逐时**；
SPINUP 见 `configs/model.parameters` 第 5 行，0=全程）、`<pref>_irr.csv`
（h + **203 反应** IRR 快照 r1..r203）、`<pref>_ipr.csv`（h + **111×6** 过程账，**ONE_STEP
内步梯形精确积分**，行 h=k 覆盖窗口 (k−1,k]，单位 molec/cm³）。

**实测基线（V2.5+，2026-10-27 run_A 实测）**：O₃max BAU 99.3 / LNNOX 49.4 /
HVOC 95.0 / HOT 113.4 / CLOUD 80.7 / ISOPX 64.5 ppb（末 24h 均值
58.4/30.3/55.4/65.3/47.3/38.4）。相对变化：LNNOX −50%、HVOC（VOC×3）−4.4%、
HOT +14%、CLOUD −19%、ISOPX −35%。BAU 减 NOx 显著降 O₃（NOx-limited 特征保持）；
VOC×3 已近中性——V2.5 光解修复后 BAU 更贴屋脊，VOC 侧弹性减小（V2.3 口径曾 +9.5%）。
> 历史口径（V2.3）：BAU 101.6 / LNNOX 53.9 / HVOC 111.3 / HOT 112.6 / CLOUD 84.3 /
> ISOPX 73.1——V2.4→V2.5+ 的光解槽位对齐与 GLY 三分支等修复使绝对值整体下移 2–15%，
> 体制方向结论不变。

## 4. 实验 B：响应面与指示因子屋脊（fac_batch_omp）

```bash
make batch                                    # 9×9 (fnox×fvoc), 线程级并行
# 四维扩展 (jfac 列表 × ΔT 列表 — SCEN_CLOUD/HOT 的响应面版):
cd facsimile_impl && ./fac_batch_omp 5 ../output/b4.csv ../runtime_v3.dat 0.6,1.0 0,5
```

> ⚠ **jfac 语义（R11 修复后）**：CLI 第 4 参数 `jfac` 是**绝对光解缩放因子**，
> 与 `meteo.config` 的 JFAC 是"同一旋钮的两种设法"——`cfg JFAC=0.8 + CLI 1.0`
> 与 `cfg 1.0 + CLI 0.8` 结果**逐位相同**（旧实现二次相乘，差 8%；已修 + 已加等价断言）。
> 网格维扫描时通常固定 cfg JFAC=1.0、只变 CLI。
>
> ⚠ **末窗样本数**：batch 统计 **24** 个样本（k>nsteps−24），fac_v3 输出 **25** 个（含端点）——
> 均值口径差 ~4%，属双胞胎容差的已知来源。

**输出**：`output/batch_v3_indicators.csv` **24 列** =
网格 4 (fnox,fvoc,jfac,dt) + 地面量 4 (O₃max/O₃mean/HCHOmax/NO₂mean) +
第一类指示因子 6 (FNR/H₂O₂:HNO₃/O₃:NOy/HCHO:NOy/O₃:NOz/NOy) +
第二类过程量 3 (L_NOx:Q/P(H₂O₂):P(HNO₃)/OPE) + 午后口径 6 (*_pm) + nst。
附 `<outf>.csv_ipr.csv`（每案例末 24h IPR 均值）与 `<outf>.csv_irr.csv`
（每案例末 24h IRR 均值，203 反应）。

**判定**：O₃–NOx 屋脊（fvoc=1 列 O₃max 峰位）；FNR 过渡带 1–2（Sillman）；
L_NOx/Q 穿 0.5（Kleinman）；OPE 谷在饱和区。体系诊断表：`python3 experiments/run/run_B_surface.py`。

## 5. 实验 C：五路线屋脊对比（ridge5）

```bash
python3 experiments/run/run_C_ridge.py          # 屋脊跨路线/跨机制对比（CSV）
python3 experiments/viz/fig_ridge5.py        # → figures/ridge5_routes.png
```

逐路线调用 `mcm_run`（fvoc=1），指示因子由
`analysis_lib/indicators.py`（路线无关引擎，CSV+FORMULA 后处理）计算。
**已知结论（V2.3 重跑）**：batch(74sp) 的 FNR 过渡一致性 0.67（贴屋脊）；mcm 的 OPE 谷精确命中（0.00）；
15sp 层退役后两路线对比。

## 6. 实验 D：IPR/IRR 过程诊断（五路线）

```bash
python3 experiments/run/run_D_irr.py            # IRR 反应级归因（top-N + 逐时）
python3 experiments/post/ipr_closure_check.py /tmp/v3w    # fac_v3 分解完备性 (比值≈1)
```

- **IPR 六账语义**（列序 `物种_P/L/E/D/B/X`）：机制产/耗、排放、干+湿沉降、
  BC 张弛+稀释、BLH 膨胀。fac_v3 行值为**内步梯形精确积分**（每小时窗）；
  其余路线（batch）为末 24h 窗口均值（快物种与大项对消预算有混叠局限，见 §9）。
- **完备性判据**：`ipr_closure_check` 各物种比值 ≈1（当前 1.01–1.02，均匀 +2%
  为内步梯形与积分器内部求积的口径差）。
- **跨求解器 IRR**：fac_v3 (fcmix) vs batch (C API+OMP) 同情景 203 反应
  末 12h 均值逐项对照——双胞胎验收。

## 7. 实验 E：指示因子工具箱（indicators.py）

```python
import sys; sys.path.insert(0, "analysis_lib")
from indicators import route_indicators, route_timeseries, Mechanism
d = route_indicators("output/expA_BAU", "mechanism/mechanism_3_full_chain.fac")   # 窗口均值
ts = route_timeseries("output/expA_BAU", "mechanism/mechanism_3_full_chain.fac")  # 逐时序列
```

NOy 物种由 .fac FORMULA 段自动判定（含 N 物种，N₂O₅ 计双 N）。
新机制接入 = 机制文件带 FORMULA 段即可，引擎零改动。

## 8. 辅助工具

> **R16 光解模式**：`photolysis.config` 的 `MODE` 现真实生效——`MCM`（默认）
> 施加天顶角调制 `l·cosᵐ·exp(−n secχ)`；`DIRECT` 不调制（白天恒定 `l`，夜间 0）。
> 修复前两者结果逐位相同（死开关），实测差异 1.3%。
>
> **R13 数值格式**：CSV 输出用 `ES15.5E3`（强制 3 位指数带 E）。历史文件若含
> `3.07932-174` 形态，用 `analysis_lib/csvfloat.py` 的 `fix_fortran_exp()` 解析
> （全部实验脚本已接入）。

```bash
make nc          # CSV → netCDF4 (默认 output/expA_BAU → .nc; 可 NC_SRC=/NC_OUT= 覆盖)
                 # ⚠ 需 <prefix>_conc.csv + _irr.csv 配对 (fac_v3 格式); batch 网格汇总不可用
make figures     # 两张响应面图 (从 output/batch_v3_indicators.csv 重生成)
                 # R12: batch 现直接写此 canonical 名 → make batch && make figures 数据链已闭合
python3 audit/element_audit.py mechanism/mechanism_3_full_chain.fac   # 203/203
python3 experiments/post/ipr_closure_check.py <prefix>                     # 任意 fac_v3 输出
```

## 9. 已知边界（实验解读必读）

1. **快物种 IPR 频闪**：OH 寿命 ~1s ≪ 小时快照——其六账占比仅指示性；
   fac_v3 内步累积后已大幅缓解（完备比值 1.02），快照路线（mcm；fac_run/kpp 已退役）
   仍受限于小时采样。
2. **大项对消预算**：白天 O₃ 的 P/L 各 ~1e7 对消到净值 ~1e5——快照路线的
   小时采样对该净值有混叠；fac_v3 内步累积无此问题。
3. **kpp 纯机制无汇**（已退役，历史口径）：IC 空间 O₃ 单调 → 无真屋脊（路线定位而非缺陷）。
4. **batch 网格上限 32×32**；jfac/ΔT 列表 ≤16 项。
5. **mcm 分支写法**：MCM 的 `CH3O2 = CH3OH` 型支链行非独立反应（元素审计 3 FAIL
   为表征风格，非机制错误）。
6. **哨兵基线**：`tests/functional_test.py` SENT 表 4 金值（屋脊上升段 O₃max 139.04 /
   峰值 O₃max 144.58 / FNR 转折带 5.239 / OPE 饱和区谷 3.832，V2.5+ 标定）——
   机制或求解改动后若超 ±10–15% 需有意识重标定。FT0 防线：生成码 ns/nr 必须等于
   机制 VARIABLE/反应行计数。


## 10. 实验 F：观测约束 OSSE（run_F_osse.py）

```bash
python3 experiments/run/run_F_osse.py        # 双胞胎 OSSE 全自动 (真值→采样→约束→技能)
```

**两条同化通路（按物理语义选择，勿混用）**：

| 观测类型 | 通路 | 配置 | 预期 |
|---|---|---|---|
| 地面 O₃ 站点 | 张弛 nudging | bc.config `O3 = file:…` + TAU_MIX | ⚠ 对化学活跃物种劣化技能 |
| VOCs 组份站点 | 同上 | `<物种> = file:…` | 短寿命物种 (HCHO/ISOP) 相对安全 |

**关键实证**：NO₂ 受光稳态支配——直接 nudge 卫星转化 NO₂ (τ=2h) 会破坏 PSS
（历史版本实测 h→1e-79 崩溃）；站点 nudge 全档劣于自由运行（O₃ 1.25→2.10+）；
**结论：约束排放入口（因果），
勿强约束被动物种（症状）**。IPR B 账 = 同化倾向的过程级度量（τ=6h: O₃ 占 2%,
HCHO 13%, ISOP 1%）。

**五路线观测约束能力矩阵**：fac_v3/fac_batch_omp ✓（chem_phys_mod bc_s 全功能）；
fac_run/mcm_run ✗（常数背景，机制验证定位）；kpp_run ✗（纯机制，IC 即约束）。


## 11. 实验 G：逐小时过程分析（hourly_pa.py）

```bash
python3 experiments/run/run_G_consist.py <fac_v3输出前缀>   # 缺省现跑 SCEN_BAU 基线
```

产出：① 关键物种 (O3/NO2/HCHO/HNO3/ISOP/PAN) 逐小时六账占比表；
② O₃ 逐小时反应级 top 预算（IRR 积分重构）；③ IPR↔IRR 两独立累积路径
互证（当前 0.000%）；④ 日循环图 `figures/hourly_pa.png`。

**契约注意**：fac_v3 的 `<pref>_ipr.csv`/`_irr.csv` 行 = 小时窗精确积分
（molec/cm³），非快照；快照路线（fac_run/kpp/mcm）逐小时分析有混叠局限
（大项对消预算尤甚）——其 IPR 用于量级与占比，不作逐时净值解读。


## 12. 实验 H：排放敏感性 OAT（run_H_sensitivity.py）

> 🚫 **伴随/反演已移除（2026-10-27）**：0-D 箱体不可实现真伴随；原"FD 雅可比 + 截断
> SVD 反演"实测仅 1/9 模态可辨识——功能与载体 `sens_adj.py` 一并删除，
> 

```bash
python3 experiments/run/run_H_sensitivity.py [扰动%]    # OAT 排放扰动 → 雅可比 CSV
```

产出：① 排放 OAT 雅可比（`output/sens/jacobian.csv`，观测 24 维 × 排放 9 物种）；
② 跨路线斜率一致性（OAT 局部 vs batch 网格弦）。
真速率灵敏度谱已随 fac_run 退役（未重建）。


## 13. 实验 I：五路线情景矩阵（scen5.py）

```bash
python3 experiments/run/run_I_scenario_matrix.py    # 多情景 × 多路线矩阵
```

情景集与方向一致性矩阵见 README §10i。**解读**：
- 同向 = 体制一致；反向 = 路线基态处于不同 NOx 体制（查 FNR/屋脊位置再解读）；
- IC×2 ≈ 0 是正确物理（初值 5 天洗脱）；BG+50% 是强杠杆（气候情景）；
- 排放缩放三路径一致性（情景表/config/网格）已验证 ≤0.1%。

## 14. 故障排查手册（按症状索引）

### 14.1 构建与运行时

| 症状 | 根因 | 处置 |
|---|---|---|
| `Cannot open file 'gen_v3/oz.reac'` | cwd 不在 `facsimile_impl/` | `cd facsimile_impl` 再执行 |
| `make: Nothing to be done` 但结果异常 | 时间戳同秒（make 误判） | **已修**：`runtime_v3.dat` 依赖 `FORCE` 每次强制重建；Fortran 侧仍可 `touch <源文件> && make v3`；FT-CFG 测试会检出陈旧 runtime |
| `[read_runtime FATAL] BC 行 token 数 N != 149` | `runtime_v3.dat` 与机制物种数不同步（行被 1024 字符截断是旧 bug，已修） | `make runtime_v3.dat` 刷新；若仍报，检查机制 VARIABLE 行数 |
| `Bad integer for item 1 in list input` @runtime_dat_mod:208 | 同上：runtime.dat 段结构与 ns 不匹配 | 重跑 cfgbuild；核对 `RUN <ns>` 与 `SPECIES` 行一致 |
| 运行瞬间结束、O3max 极小 | 初始 IC 全零 → 光化学无法启动 | 检查 `initialConcentrations.config` 是否有 O₃/NO |
| `fcvode` 返回负值 / `FAIL @step` | 步长崩溃（刚性/非光滑 RHS） | 见 §14.4 数值类 |

### 14.2 配置编译（cfgbuild）

| 症状（errs 输出） | 根因 | 处置 |
|---|---|---|
| `IC 未知物种 'X'` / `BC 未知物种` / `排放未知物种` | 拼写与机制 VARIABLE 不一致（**大小写敏感**） | 对 `runtime_v3.dat` 的 SPECIES 行核对 |
| `IC 未知单位 'ug/m3'` | 只接受 cm3/ppb/ppm | 用模板换算表（AQI: μg/m³→ppb 见 `docs/templates/surface_aqi_station.template`） |
| `排放 'X': 日型 'box' 非法` | 日型大小写敏感 | 用大写 `BOX` / `FLAT` |
| `排放 'X': BOX 振幅 1.2 ∉ [0,1]` | 振幅越界 | 峰值超基线无依据；改用 amp≤1 + 提高基线 |
| `序列时间非单调 @t=…` | constraints 文件 t 列非严格递增 | 排序去重；确认单位为**秒**（非小时/分钟） |
| `序列 N 点 > Fortran 端上限 128` | 单序列超 128 点 | 降采样（如 5 分钟→1 小时） |
| `机制使用光解槽位无参数 (J=[…])` | JMAP/光解表未覆盖机制的 `J<n>` | 在 `photolysis.config` 的 JMAP 加 `槽位=D值` |
| `photorolysis-rates_* 存在但未解析出任何槽位` | 表格式错（D 记法/列序） | 对照表头 `j l m n name tau` |

### 14.3 机制修改（mechc + 审计）

| 症状 | 根因 | 处置 |
|---|---|---|
| `[mechc] VARIABLE 重复物种: [...]` | 物种写了两遍（同行或跨行） | 删除重复（R9 新增守卫） |
| `[mechc] 重复速率定义: KX` | 同名 K 定义两次（后者静默覆盖前者） | 重命名或删除（R9 新增守卫） |
| 审计 `未知物种 '0.42 O2 ;  ! 注释'` | **反应行尾写了 `!` 注释** | 注释移到独立 `*` 行（`.fac` 解析器不剥行内注释） |
| 审计 `净差 C:+0.09 H:+0.21 O:-1.82` | 产物系数不配平 | 按 `mechanism/README §8.3` 的配平案例迭代；注意分数精度（1/4 可，0.145 不可） |
| 审计通过但化学崩（负浓度） | 速率量级错挂（如 KET+OH 挂成 KRO2NO，快 40×） | 用"零引用速率扫描"检查（见 §14.5） |
| FT0 `gen_v3 陈旧! ns=… != …` | 改了 .fac 没重编译 | `make v3` |
| FT0 `nr=… != …` | 加了反应没重编译 | `make v3` |
| 哨兵超容差 | 机制/求解改动导致数值漂移 | 确认物理合理后重标定 SENT 表（±10–15%） |

### 14.4 数值类（步长/收敛）

| 症状 | 根因 | 处置 |
|---|---|---|
| 步长崩到 1e-6，steps 数万 | **RHS 有折点**（`max(y,0)` 类） | 用原始 y（10b 坑 #11：移除后 15× 加速） |
| 黎明时段崩溃/mxstep 耗尽 | 无预条件 SPGMR + 光解跃变 | 单箱用 DENSE（`FCVDENSE`）；批量才考虑 Krylov |
| 重启后结果放飞（MDA8 数值爆） | 韧性重启用了松容差 | `FCVREINIT` 必须同容差 + 小步长（10b 坑 #6） |
| 浓度出现负值且被"钳到 0" | BDF 不保正性 + 钳零破坏守恒 | 不钳零（KPP/CMAQ 惯例）；接受小幅负值 |
| 氮守恒不变量 FAIL（|Σw_N·dy/dt| ~1e3） | 新含 N 物种未进 `noy_add` 权重表 | `fac_batch_omp.f90` 加一行（权重=N 原子数） |

### 14.5 诊断与后处理

| 症状 | 根因 | 处置 |
|---|---|---|
| CSV 值 `8.27280-108` 解析失败 | gfortran 三位指数省略 `E`（ES12.3 时偶发） | `diagnostics._fix_fortran_exp` 已兜底；自定义解析需同款处理 |
| `ipr_closure_check` 比值远离 1 | IPR 与 conc 窗口错位 / 该物种无过程 | 核对 `h` 列对齐；检查物种是否真的被反应触碰 |
| `scipy.io.netcdf_file` 读 netCDF 输出报错 | 文件是 **NETCDF4**（HDF5），`scipy.io.netcdf_file` 只认 NETCDF3 | 用 `netCDF4.Dataset` 或 `xarray.open_dataset(engine="netcdf4")` |
| nc 里找不到 `variables['O3']` | 物种名在**属性** `species_names`，变量是 `conc(species,time)` | 用索引 + 属性名列表 |
| 并行/串行 `diff` 不通过 | **CSV 行序非确定**（OpenMP `schedule(dynamic)` 导致写盘顺序随机） | 先按键 `(fnox,fvoc)` 归位再逐行比，或按行排序后比 |
| batch 耗时突然翻倍 | 系统负载 > ncpu | 看 `uptime`；负载高时的墙钟加速比无意义，应看 CPU 时间 |
| FNR 值 > 10 或 < 0.01 | 边缘格点（极低/极高 NOx） | 正常；判据只在 0.5–5 区间有物理意义 |

### 14.6 环境与工具链

| 症状 | 根因 | 处置 |
|---|---|---|
| `ModuleNotFoundError: numpy` | 用了系统 python3 而非 venv | 用 `$(PY)`（Makefile 变量，默认 `$HOME/miniconda3/envs/venv/bin/python`；可用 `make PY=...` 覆盖） |
| `ld: warning: built for newer macOS` | CVODE dylib 目标版本 12.7 > 本机 12.0 | 无害警告，可忽略 |
| `make batch` 报 `No rule to make target` | `runtime_v3.dat` 缺失 | `make v3 && make batch`（batch 依赖二者） |
| 改了 config 但结果不变 | 没刷新 `runtime_v3.dat` | `make runtime_v3.dat`（0.37 s） |
| 改了 drydep/排放但浓度不变 | 同上是 runtime 未刷新 | 同上（**最高频的"假 bug"**） |

---
