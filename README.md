# ozone_lab —— 对流层臭氧光化学程序族

> **定位**：独立程序系列目录（非 Day xx）——多求解器 × 全谱诊断 × 完备输入系统。
> 两档在役机制（**工程版 mechanism_3 V2.5+ 冻结版：111 物种 203 反应**，十段链条骨架
> + V2.2 集总类 + V2.3 硝酸酯/单萜补链 + V2.4 芳香分级/长链烷烃/醇醚酯/倍半萜/乙烯/
> 内烯烃/高级醛扩展 + V2.5 光解槽位对齐/GLY 三分支/O₃+HOx/烯烃-NO₃ 补正；
> 真实 MCM CH₄ 子集 29 物种 71 反应；15sp 教学档已删除），
> **CVODE 双内核**（fcmix / C-ABI+OpenMP）路线，**双库诊断**（indicators.py 指标引擎 +
> diagnostics.py 六账积分），**12 类输入文件**完备处理（runtime.dat 通路全接线）。
> **路线现状（2026-10-27）**：现役 **fac_v3 / fac_batch_omp / mcm_run** 三条；
> fac_run、kpp_run、kpp_batch、EBI、fac_expanded、fac_spgmr 已退役（详见 §1.1）。
>
> **七大分析能力**：①观测约束（nudging 同化）· ②响应面（4D 参数扫描，最高 262,144 案例）
> · ③情景分析（6 政策情景）· ④过程分析（IPR 六账 + IRR 203 反应）
> · ⑤敏感性 OAT · ⑥反演（⚠ 非真伴随，FD+SVD）
> · ⑦气象序列驱动（6 变量逐时观测）——索引见
> [OPERATIONS_MANUAL §0.7](docs/OPERATIONS_MANUAL.md)，详解见
> [EXPERIMENTS_GUIDE §6.9–6.11](docs/EXPERIMENTS_GUIDE.md) 与 [§0.3b](docs/EXPERIMENTS_GUIDE.md)。
>
> **环境要求**：gfortran 13.2 + SUNDIALS CVODE 2.7.0（含 `--enable-fcmix`）+
> Python 3.9+（numpy/matplotlib）。**路径可通过环境变量覆盖**，无需改文件：
> ```bash
> make all CONDA_ENV=$HOME/miniconda3/envs/ozone \
>          GCC_ROOT=$HOME/gcc-13.2.0 \
>          CVODELIB=$HOME/soft/atmenv/cvode/lib
> ```
> 详见 [ENVIRONMENT.md](environment/ENVIRONMENT.md)；一键体检 `sh environment/check_env.sh`。
>
> **一键复现**（计数为 2026-10-27 实跑值）：
> ```bash
> make test   # 构建 → 冒烟 1/1 → 功能 62/62 → 元素审计 203/203 → 81 案例响应面 + 批量自检
> make v3 && make batch   # 工程版: fac_v3 情景 + fac_batch_omp OpenMP 响应面
> make mcm                # MCM 对照路线 (29sp 独立子树, 可 clean 后重建)
> make figures            # response_surface_v3.png + indicators_vs_ridge.png
> ```
> 本 README 功能描述**逐文件实读源码核对**；与旧文档不符处见 §11 勘误表。

---

## 目录

- [0. 快速开始](#0-快速开始)
- [1. 程序功能总表（根据源码）](#1-程序功能总表根据源码)
  - [1.1 可执行程序清单（现役 3 个）× 机制版本映射](#11-可执行程序清单现役-3-个-机制版本映射)
  - [1.2 Python 工具链（8 个）](#12-python-工具链8-个)
  - [1.3 审计与测试（5 个）](#13-审计与测试5-个)
  - [1.4 机制文件（3 份真源 + 退役档）](#14-机制文件3-份真源--退役档)
  - [1.5 输入系统（12 类配置 → runtime.dat）](#15-输入系统12-类配置--runtimedat)
  - [1.6 输出契约（CSV 精确列序）](#16-输出契约csv-精确列序)
- [2. 可以做的分析与实验（8 类 18 项）](#2-可以做的分析与实验)
  - [A 类：单箱基线与情景实验](#a-类单箱基线与情景实验)
  - [B 类：响应面与 O₃ 体制判别](#b-类响应面与-o₃-体制判别)
  - [C 类：求解器方法学实验](#c-类求解器方法学实验)
  - [D 类：诊断分析（双库）](#d-类诊断分析双库)
  - [E 类：机制工程实验](#e-类机制工程实验)
  - [F 类：数值方法教学实验](#f-类数值方法教学实验)
  - [G 类：性能优化实验](#g-类性能优化实验)
  - [H 类：机制图集与教学](#h-类机制图集与教学)
- [3. 架构总览](#3-架构总览)
- [4. 输入系统](#4-输入系统)（12 类文件格式详表 + 处理链）
- [5. 双求解器实现](#5-双求解器实现原三路线ebi-已退役)
- [6. 诊断库](#6-诊断库)（源码函数级）
- [7. 测试体系](#7-测试体系)
- [8. 运行结果汇总](#8-运行结果汇总)
- [9. 性能分析与优化](#9-性能分析与优化)
- [10. 功能拓展方案](#10-功能拓展方案)
- [10a. 工程落地版 V3（源码级）](#10a-工程落地版-v3源码级)
- [10b. 数值攻坚史（14 坑实录）](#10b-数值攻坚史14-坑实录)
- [11. 已知限制与勘误](#11-已知限制与勘误)
- [12. 参考链接](#12-参考链接)

---

<a name="0-快速开始"></a>
## 0. 快速开始

```bash
# ① 构建 + 全测试链
make test

# ② 工程版单箱情景模拟 (111 物种 203 反应 V2.5, 六情景任选)
make v3
cd facsimile_impl && ./fac_v3 SCEN_BAU ../output/base_v3 ../runtime_v3.dat

# ③ 81 案例响应面 (OpenMP, 输出 24 列指示因子)
make batch            # = v3 + runtime_v3.dat + fac_batch_omp 9

# ④ 诊断一键报告 (指示因子/IRR/family/账本)
python -c "import sys; sys.path.insert(0,'diagnostics')
from indicators import route_indicators
print(route_indicators('output/batch_v3_indicators', 'mechanism/mechanism_3_full_chain.fac'))"

# ⑤ 机制图集 + 工程图
python figures/gen_all.py
make figures

# 路线选择速查: ①② 同内核双窗口 — ①显微镜(单点问为什么,逐时全账) ②望远镜(曲面问在哪转折);
#   标准流: 广度找现象 → 深度查原因 → 广度验稳健 (互证容差 ±10~15%);
#   日常模拟只用 ① fac_v3(单情景/同化) + ② fac_batch_omp(批量响应面);
# ③ mcm_run = 机制正确性体外参照 (MCM 独立来源), 机制扩链后跑一次作第四道关卡。
# 何时编译: 改配置→零编译(make runtime_v3.dat, 0.37s); 改 .fac/Fortran→make v3(~8s)。
```

---

<a name="1-程序功能总表根据源码"></a>
# 1. 程序功能总表（根据源码）

> 每个条目给出源码位置、CLI、**实读确认**的功能与输出（本轮逐行核对
> `oz_fac_v3.f90`、`fac_batch_omp.f90`、`runtime_dat_mod.f90`、`chem_phys_mod.f90`、
> `diagnostics.py`、`Makefile`、`cfgbuild.py` 解析端、全部 tests/audit）。

## 1.1 可执行程序清单（**现役 3 个**；fac_run/kpp_run 已退役见 §1.1b）× 机制版本映射

> **2026-10-27 版本现状**：工程主线合并后仅保留 **fac_v3 / fac_batch_omp / mcm_run**
> 三个可执行体（物理层统一收敛于 `modules/chem_phys_mod.f90`，情景分组单源化于
> `build_groups()`）；fac_run（15sp 教学档）与整个 kpp_impl（KPP 路线）已退役删除，
> 15sp 核心教学档降级为 .fac 参考文档。下文 ①④ 等退役条目保留作历史档案。

| 机制版本 | 唯一事实源 (.fac) | 消费路线 |
|---|---|---|
| ~~15sp / 25rx~~ 教学核心（**已删除**） | —（`mechanism/teaching/` 已于 V2.4+ 冻结时移除） | 无（fac_run/kpp_run 已退役） |
| **111sp / 203rx** 工程 V2.5+ 冻结版 | `mechanism/mechanism_3_full_chain.fac` | `fac_v3` + `fac_batch_omp`（共享 `gen_v3/`） |
| **29sp / 71rx** 真实 MCM CH₄（**独立子树**） | `mcm_city/input/mechanism_mcm_ch4.fac` | `mcm_run`（`make mcm`，见 [mcm_city/README.md](mcm_city/README.md)） |

> **机制规格书（逐反应速率/配平方法论/对比定位）**：见 [`mechanism/README.md`](mechanism/README.md)

30sp/51rx 扩展版已随版本精简退役删除（7 处配平错误，内容被后续机制吸收，
见 §11.1 与版本精简记录）。

### [退役] ~~`facsimile_impl/fac_run`~~ —— 旧教学主力（2026-10 退役；单箱/全物理能力由 fac_v3 全面接管）

- **源码**：`oz_fac_core.f90` + `gen/oz_mech.f90`（mechc.py 编译自 `mechanism/ozone.fac`，
  内容已随 `mechanism/teaching/` 一并删除）
- **CLI**：`fac_run f_NOx f_VOC out_pref`（三参必选，`fac_run 1 1 ../output/base` 为基线）
- **机制**：核心 15sp/25rx（NO NO2 O3 O1D O OH HO2 H2O2 CO HCHO CH4 CH3O2 CH3OOH CH3O HNO3）
- **实读确认**：
  - CVODE fcmix：BDF+Newton，rtol=1e-6；**ATOL 按物种调优**（O₃ 1e8 / CH₄ 1e10 / CO 1e9 / HCHO 1e7）
  - 线性解 DENSE；**SPGMR 切换 = 源码内注释行** `! call FCVSPGMR(0,1,15,0.05d0,IER)`
  - `MAX_NSTEPS=200000`，`MAX_STEP=100 s`；失败 → `FCVREINIT` 同容差 + `MAX_STEP=20 s` 重试一次
  - 物理层内嵌 resid：机制 P/L + 排放 + 沉降 k_dep + **背景松弛仅 O₃/HCHO 两物种**（`(BG−y)/TAU_MIX`）
  - IRR 钩子：每步末 `peff_out(nr)` 快照
- **输出**：
  - `<pref>_conc.csv`：**固定 10 列** `h,O3,NO,NO2,OH,HO2,HCHO,H2O2,HNO3,CH3O2`
    （HNO₃ 取 `Y(15)` 硬编码），**只写末 24 h**（step≥96）
  - `<pref>_irr.csv`：`h,r1..r25`（末 24 h 快照，5 天积分由 Python 梯形完成）

### ① `facsimile_impl/fac_v3` —— 工程版主力（情景名驱动）

- **源码**：`oz_fac_v3.f90`（256 行）+ `modules/runtime_dat_mod.f90`（318 行）+ `gen_v3/`
- **CLI**：`fac_v3 [SCEN_NAME|序号] [out_pref] [runtime.dat]`（缺省 `output/facv3`、`runtime_v3.dat`）
- **机制**：`mechanism/mechanism_3_full_chain.fac` → **111sp/203rx V2.5+ 冻结版**（十段链条骨架 + V2.2 集总类；V2.3 烷基硝酸盐+单萜；V2.4 芳香分级/长链烷烃/醇醚酯/倍半萜/乙烯/内烯烃/高级醛；V2.5 光解槽位对齐/GLY 三分支/O₃+HOx/烯烃-NO₃ 补正）
- **实读确认的物理链**（`facv3_mod.resid_ipr`）：
  ```
  ydot = Σ(机制反应计量 × peff)                     ← 机制 P/L（peff 用原始 y, 不夹 max(y,0)）
       + 排放: emis0·scale·box_shape(BOX) | emis0·scale(FLAT)
       − 干沉降: vd/BLH·y   − 湿沉降: wet_a·max(RAIN,0)^wet_b·y
       + BC 松弛: (bc_c−y)/τ_mix (bc>0 物种) + 全物种稀释 dil_k
       − BLH 膨胀: y·dB/BLH (仅 dB>0 增长段)
  ```
  **六本账 IPR 列**（源码注释逐字）：`1 机制P | 2 机制L | 3 排放 | 4 干+湿沉 | 5 BC松弛+稀释 | 6 BLH膨胀`
- **数值细节**：
  - ATOL 按名调优表：CH₄ 1e10 / CO 1e9 / O₃ 1e8 / HCHO·PAN·MPAN·ISOP·SO₂·CH₃CHO·MVK·MACR·MGLY 1e7 / DMS·MSA 1e6
    （**RTOL/ATOL 含义、实测敏感性与"何时该动"见 [GUIDE §3⑬](docs/EXPERIMENTS_GUIDE.md)**）
  - DENSE（111×111；无预条件 SPGMR 黎明崩塌——10b #12）；`MAX_STEP=100 s`，重启后 10 s
  - 输出目录自动 `mkdir -p`；情景名匹配 → 序号解析 → 失败回落情景 1（`ios` 显式置 0）
- **输出**：`<pref>_conc.csv`（表头 `h,<111 物种名动态生成>`，**末 24 h**）+
  `<pref>_irr.csv`（`h,r1..r203`）+ `<pref>_ipr.csv`（`h` + 111×6 账）
  + stdout（`FAC-V3 ok (…s, nst=…, O3max=… ppb)`）

### ② `facsimile_impl/fac_batch_omp` —— OpenMP 批量响应面

- **源码**：`fac_batch_omp.f90` + `cvode_c_mod.f90`（ISO_C_BINDING 直连 CVODE C API）
- **CLI**：`fac_batch_omp [ng] [outf] [rtfile]`（默认 9、`../output/batch_v3.csv`）
- **实读确认**：
  - 网格 `fs(i)=0.3×1.5^(i-1)`（ng=9 → f∈[0.3,7.7]；上限 32）；**每案例固定 fco=fch4=jfac=1, dts=0**——只扫 f_NOx×f_VOC
  - 每线程独立 mem/N_Vector/user_data；`schedule(dynamic)`；情景生效态 `threadprivate`
  - 启动自检 `diag_selfcheck()`（IRR 提取目标匹配 + 氮守恒不变量首快照）
- **输出**：CSV 表头（源码逐字）：
  ```
  fnox,fvoc,O3max,O3mean,HCHOmax,NO2mean,FNR_mean,
  H2O2HNO3,O3NOy,HCHONOy,O3NOz,NOy_ppb,LnOxQ,PH2O2PHNO3,OPE,nst
  ```
  = 网格坐标 2 + 地面量 4 + **指示因子 9**（FNR/H₂O₂:HNO₃/O₃:NOy/HCHO:NOy/O₃:NOz/NOy/L_NOx·Q⁻¹/P(H₂O₂):P(HNO₃)/OPE）+ 步数。
  另有 `batch_v3_indicators.csv` 同构（诊断列扩展）。

### [退役] ~~`kpp_impl/kpp_run`~~ —— KPP 机制基准（2026-10 退役，kpp_impl/ 整目录删除）

`driver_kpp.f90` + KPP 3.2.0 生成码（`build_kpp.sh` 五补丁链）。**纯化学无物理层**
（算子分裂回写破坏 Rosenbrock 步长历史，实测 OH→1e14）——定位机制精度参照。

### ③ `mcm_city/mcm_run` —— MCM 真实机制城市案例

`mcm_run [days] [out_pref]`。MCM v3.3.1 CH₄ 子集 29sp/71rx/139 命名速率/7 光解；
FAC-SPGMR + 城市排放 + 韧性重启。输出 `mcm_city_conc.csv` + `mcm_city_irr.csv`。

**共享生成码**：`gen*/oz_mech.f90`（ns/nr/nq/mre/mpe + `update_rates` + 计量表）
+ `oz.reac/oz.prod`——机制改一处（`.fac` 为唯一事实源），全家自动同步：
mechc.py → gen*/、fac2kpp.py → KPP 三件套、cfgbuild --mech → runtime_v3.dat。

## 1.2 Python 工具链（8 个）

| # | 脚本 | 实读确认的功能 | 用法 |
|---|---|---|---|
| 1 | `src/cfgbuild.py` | 12 类配置 → Fortran 顺序表控 runtime.dat；段结构（与解析端 `runtime_dat_mod.read_runtime` 逐段对上）：`RUN ns nj nser nscen / SPECIES / MODEL days dtout rtol atol / IC / MET×7 / BC / PHOT lat doy + J idx l m n ×nj / EMIS×ns / DRY×ns / WET×ns / DIL / SCEN name fnox fvoc fco fch4 jfac dt / SERIES name n + t行 + v行`；编译期校验全过才写出；`--days/--dtout` 临时覆盖 | `python src/cfgbuild.py input_v3 runtime_v3.dat` |
| 2 | `facsimile_impl/mechc.py` | **FACSIMILE→Fortran 机制编译器**：`.fac` → `gen_v3/oz_mech.f90`（ns/nr/nq/nj + update_rates）+ `oz.reac/oz.prod` 计量表；V2.5 硬化（行内 `!` 注释剥离 + 未知 token 硬失败） | `make v3`（自动） |
| 3b | `src/sat2surf.py` | **卫星柱 → 地表等效 ppb 转换器**（两法互证/QA 过滤/BLH 序列插值/点数抽稀/量级自检；`--selftest` 内置校验） | `python src/sat2surf.py <raw> <out> [--blh-file ...]` |
| 4 | `diagnostics/diagnostics.py` | CSV 读取与六账积分纯函数（§6.2） | `import diagnostics as dg` |
| 4b | `diagnostics/indicators.py` | **路线无关 12 指标引擎**（NOy 从 FORMULA 段自动判定，§6.1） | `from indicators import route_indicators` |
| 4c | `diagnostics/csvfloat.py` | gfortran 三位指数省略 E 的健壮解析共享件 | 被 diagnostics/nc_output 引用 |
| 5 | `src/fac2kpp.py` | ozone.fac → KPP 三件套转换脚本——**KPP 路线已退役，脚本保留作历史工具（无 make 目标）** | 手动运行 |
| 6 | `figures/fig_indicators_ridge.py` + `fig_response_v3.py` | 消费 `output/batch_v3_indicators.csv` → 两张工程图（`make figures` 调用，自检轮次 3 落盘可复现） | `make figures` |

（另：`modules/nc_output.py` 为 CSV→netCDF 转换，见 `make nc`。）

## 1.3 审计与测试（5 个）

| # | 文件 | 实读确认的功能 | 用法 / 实测 |
|---|---|---|---|
| 1 | `audit/element_audit.py` | 逐反应 C/H/N/O/S/Cl/Br 平衡；储库（H₂O/O₂/CO₂/H₂）自动推断但**限定物理方向**（CO₂/H₂ 只许产物侧、O₂ 反应物侧≤2.5）；111 物种元素表内嵌（机制自带 FORMULA 块覆盖） | `python audit/element_audit.py mechanism/mechanism_3_full_chain.fac` → **203/203** |
| 2 | `audit/mcm_compare.py` | 与 MCM v3.3.1 官方导出逐反应系数比对（反应物集合匹配 / 多分支求和 / 集总映射 / 文档漂移检查） | `python audit/mcm_compare.py` → 最大偏差 8.9% |
| 3 | `tests/smoke_test.sh` | FAC-V3 冒烟（log 含 ok + O₃ 末列非零） | `make smoke` → 1/1 |
| 4 | `tests/functional_test.py` | **62 断言**（FT-CFG/CFGSPEC/DAYS/FREEZE/GIT/MECHK/MP/PHOT/SAT/SYNC/TPL/UNIT 12 族）：编译期负测试 / batch 24 列契约 / 屋脊 / FNR 单调 / IPR 六账闭合 / L_NOx-Q 值域 / 金值哨兵 / 冻结版规格三向一致 / mechc 计量表防御 / output 不入库 | `make functional` |

> ⚠ 原 `oz_audit.c` / `cfg_audit.c`（C 证人）已于 2026-10 删除：从未被编译或调用
> （无 Makefile 规则、无 ctypes 消费者）。同款校验由 `element_audit.py`（元素守恒）、
> `fac_batch_omp` 氮守恒在线不变量、`run_D` + `FT-IPR`（六账闭合）承担。

## 1.4 机制文件（3 份真源 + KPP 生成物）

| 文件 | 规模 | 编译去向 | 实读要点 |
|---|---|---|---|
| ~~`mechanism/ozone.fac`~~（已随 fac_run 退役删除） | 15sp/25rx | — | 15sp 核心教学机制已删除（原 `mechanism/teaching/`） |
| `mechanism/mechanism_3_full_chain.fac` | **111sp/203rx V2.5+** | fac_v3 / fac_batch_omp `gen_v3/` | 十段链条骨架（①无机核心+HONO 晨爆 ②NO₃-N₂O₅ 夜化学 MCM Troe 逐字 ③CH₄+PAN 储库 ④**isoprene 完整二代链** ⑤DMS→MSA ⑥卤素+储库 ⑦集总类 ⑧支撑化学 ⑨硝酸酯家族 ⑩单萜+V2.4 扩展）；GLY 光解三分支；每反应 C/H/N/O/S/Cl/Br 精确配平（203/203） |
| `mcm_city/input/mechanism_mcm_ch4.fac` | 29sp/71rx | mcm_city `gen/` | AtChem2 自带 MCM v3.3.1 CH₄ 子集；完整无机组（NO₃/N₂O₅/HO₂NO₂/HONO）+ Troe（KMT01-03/KBPAN/KFPAN）+ HSO₃/SO₃ 硫链 + CL |
| ~~`kpp_impl/ozone_kpp.{spc,eqn,def}`~~（已随 KPP 路线退役删除） | — | — | `src/fac2kpp.py` 转换脚本保留 |
| （已删）`input_v3/mechanism_3_full_chain.fac` | — | — | 旧副本已删除；cfgbuild 经 `--mech mechanism/mechanism_3_full_chain.fac` 直引真源 |

## 1.5 输入系统（12 类配置 → runtime.dat）

**唯一工程输入目录是 `input_v3/`**（12 类配置，文件名沿用 AtChem2 约定）。
旧的 `input/`（15sp 教学同构份）已于 2026-10 彻底删除——它零消费者且易被误
当工程入口；15sp 教学机制档亦已删除（不再保留死代码）。
解析端（`runtime_dat_mod.read_runtime`）逐段实读确认的 12 段：

```
RUN(ns nj nser nscen) SPECIES MODEL(days dtout rtol atol) IC MET(7行)
BC(τ+2ns列) PHOT(lat doy + J idx l m n ×nj) EMIS(×ns) DRY(×ns)
WET(×ns) DIL SCEN(name fnox fvoc fco fch4 jfac dt ×nscen)
SERIES(name n / t… / v… ×nser)
```

**气象与光解求值**（`eval_met/eval_env/eval_J` 实读）：
- 7 量 TEMP/PRES/RH/RADIAT/BLH/RAIN/JFAC——常数或序列线性插值（**外延取端值**）
- `dBLH`：±600 s 中心差分（cm/s）；BLH 膨胀稀释**仅增长段生效**
- `M=6.022e23·P/(R·T)`，`O₂=0.2095M`，H₂O 由 Magnus 公式
- SZA→J：`decl=23.45°·sin(2π(doy−81)/365)`，`cosz=sinφ·sinδ+cosφ·cosδ·cos(时角)`，
  `J=JFAC·l·cosᶻᵐ·exp(−n/cosχ)`；`cosz≤1e-3 → J=0`
- 排放日型 BOX：`0.1+0.9·amp·sin²(π(h−6)/12)`（06–18 时，夜间 0.1 基底）；FLAT 全天常数

**情景组映射**（`build_groups` 实读——现役 4 组 28 物种，含 V2.4 排放类）：
```
f_NOx → NO, NO2
f_VOC → HCHO, CH3CHO, ISOP, DMS, MVK, MACR, MSA, SO2, ALK, OLE, ARO1, ARO2,
        KET, TERP, CSL, LALK, MEOH, ETOH, MTBE, ESTER, SESQ, ETH, IOLE, ALD2
f_CO  → CO             f_CH4 → CH4        另: jfac 乘 J 通道, dt 加 TEMP
```

## 1.6 输出契约（CSV 精确列序）

| 文件 | 精确列结构（源码逐字） | 消费者 |
|---|---|---|
| ~~fac_run conc~~（退役） | 曾为 `h,O3,...,CH3O2` 10 列；现由 fac_v3 67 列契约取代 | dg.load 需按列名适配 |
| fac_run irr | `h,r1..r25`（末 24 h 快照） | `dg.irr_integrate` |
| fac_v3 conc | `h` + **111 物种** = **112 列**（末 24 h 逐时，molec/cm³） | `diagnostics.load` / 自定义分析 |
| fac_v3 irr | `h` + **203 反应** = **204 列**（**小时窗积分**，非速率） | `dg.irr_integrate` |
| fac_v3 ipr | `h` + **111×6 账** = **667 列**（`P L E D B X`，小时窗积分） | `run_D_process` / 过程归因 |
| batch CSV | `fnox,fvoc,O3max,O3mean,HCHOmax,NO2mean,FNR_mean,H2O2HNO3,O3NOy,HCHONOy,O3NOz,NOy_ppb,LnOxQ,PH2O2PHNO3,OPE,nst`（24 列） | figures 两脚本 |
| mcm_city_* | 29 物种 + 71 反应 | 自定义分析 |

**注意**：conc/irr 只写**末 24 h**（前 96 h 建立准稳态，全程 O₃max 走 stdout）；
要全程序列改驱动的 `step >= nsteps-24` 条件。

---

<a name="2-可以做的分析与实验"></a>
# 2. 可以做的分析与实验

> 8 类 18 项，按"单箱 → 响应面 → 方法学 → 诊断 → 工程 → 教学"递进。
> 每项：**目的 / 命令 / 产物 / 解读**，命令在仓库根目录可复现。

## A 类：单箱基线与情景实验

### A1. 工程版基线 5 天模拟（BAU）

- **目的**：111 物种全链条基线——滴定建立期 → 准稳态、OH/HONO 晨爆、夜化学。
- **命令**：`make v3 && cd facsimile_impl && ./fac_v3 SCEN_BAU ../output/base_v3 ../runtime_v3.dat`
- **产物**：`base_v3_conc.csv`（112 列 = h+111 物种）+ `_irr.csv`（203 反应）+ stdout O3max
- **解读**：实测（V2.5+，2026-10-27）≈2.6 s / 4802 步 / 零重启；末 24h O₃ **58.4 ppb**
  （O₃max 99.3 ppb 走 stdout）。CSV 只含末 24 h（§1.6）——要看建立期改
  `oz_fac_v3.f90` 的 `step >= nsteps-24`。

### A2. 六情景政策沙盘

> **体制解读（V2.5+ 实测，2026-10-27 run_A）**：BAU 位于屋脊左侧 NOx-limited 侧——
> **减 NOx（LNNOX）O₃max −50%**；**VOC×3（HVOC）−4.4%（近中性）**——V2.5 光解修复后
> BAU 更贴屋脊，VOC 侧弹性减小（V2.3 口径曾 +9.5%）。HOT +14%（速率/PAN 分解/排放
> 三温度通路）、CLOUD −19%（光解第一驱动力）、ISOPX −35%（NOx 减半主导）。
>
> ⚠ **R10 勘误（历史记录）**：此前曾以"HVOC 几乎无效（−0.2%）"作为 NOx-limited 的佐证。
> 该数字实为 `scen_build` **漏缩 6 个排放物种**（SO2/ALK/OLE/ARO1/ARO2/KET/TERP）的
> 假象——VOC 根本没加。该 bug 已修复并单源化（`build_groups` + `scen_build` 只读
> grp* 数组，两条路线同源）；历次审计原始报告已随 AUDIT_COMPLETE.md 于 2026-10-26 删除，
> 防线固化为 functional_test 断言。

- **目的**：回答"减谁、减多少、O₃ 怎么变"（情景=排放组缩放×JFAC×ΔT 组合）。
- **命令**：循环 `./fac_v3 SCEN_X ../output/x_v3 ../runtime_v3.dat`
- **实测**：

| 情景 | 因子 | O₃max [ppb] | 解读 |
|---|---|---|---|
| BAU | 全 1 | 99.3 | 基线（真实清单, 屋脊左侧 NOx-limited 侧） |
| LNNOX | NOx×0.3 | 49.4 | **−50%**（NOx-limited 线性响应） |
| HVOC | VOC×3 | 95.0 | −4.4%（近中性；V2.5 后 BAU 更贴屋脊） |
| HOT | T+5K | 113.4 | +14%（速率+PAN 分解+排放三温度通路合力） |
| CLOUD | JFAC×0.6 | 80.7 | −19%（光解第一驱动力） |
| ISOPX | NOx×0.5+VOC×2 | 64.5 | −35%（NOx 减半主导） |

> V2.3 历史口径：BAU 101.6 / LNNOX 53.9 / HVOC 111.3 / HOT 112.6 / CLOUD 84.3 / ISOPX 73.1——
> V2.4→V2.5+ 光解槽位对齐与 GLY 三分支等修复使绝对值整体下移 2–15%，方向结论不变。

### A3. 排放缩放单点（自定义情景）

- **命令**：`input_v3/scenarios.config` 加 `SCEN_MYCASE = E_NOx*1.5 E_VOC*1.0 JFAC*1` →
  `make runtime_v3.dat` → `facsimile_impl/fac_v3 SCEN_MYCASE ../output/mycase ../runtime_v3.dat`
- **解读**：滴定区/稀放区的单点响应由 fac_v3 承载（历史命令 fac_run 已退役；
  FT10 的滴定区判据 1.5 vs 10.4 ppb 已固化为功能测试锚点）。

### A4. 温度敏感性链（三通路分解）

- **目的**：分离速率常数（EXP(−Ea/RT)）、PAN/MPAN 储库（`KMPANTH = 1.6e16·EXP(−13500/T)`，
  Ea≈112 kJ/mol 全机制最强；PAN 热解离走 `KBPAN` Troe 形式）、排放三通路。
- **命令**：SCEN_HOT vs BAU；对比 irr 中 PAN/MPAN 热分解列（r 序号查 `gen_v3/oz.reac`）。

### A5. 光强/云量敏感性

- **命令**：SCEN_CLOUD（JFAC×0.6）；或 meteo.config 的 JFAC 设 `file:` 序列做逐时云过程
- **解读**：CLOUD vs BAU O₃max −18.6 ppb（V2.5+ 实测 80.7 vs 99.3）；O₃ 峰时后移。

### A6. 夜间化学 + HONO 晨爆可视化

- **目的**：六链条②——NO₃/N₂O₅ 夜化学、HONO 表面积累 → 日出爆发（教学机制无此过程）。
- **命令**：A1 产物取 NO3/N2O5/HONO/OH 四列画 5 天日循环。
- **解读**：HONO 夜积累、日出 ~1 h 耗尽（J(HONO)≈1e-3 s⁻¹，τ≈14 min）；N₂O₅ 峰在午夜后。

## B 类：响应面与 O₃ 体制判别

### B1. 9×9 O₃–NOx–VOC 响应面（O₃ 屋脊）

- **命令**：`make batch` → `output/batch_v3.csv`（81×16）→ `make figures`
- **解读**（V2.5+ 哨兵锚定）：fvoc=1 切面 O₃ 随 NOx 升至 **139 ppb（fnox≈2.28，屋脊
  上升段）**、峰值 **144.6（fnox≈3.4）**后转入 NOx-saturated 回撤分支——屋脊即
  NOx-limited→saturated 转折；几何网格 f=0.3×1.5^(i-1)（0.3–7.7）。

### B2. 指示因子全家族 × 屋脊交叉验证

- **目的**：**方法学核心实验**——Sillman（1995）浓度比判据与 Kleinman（2005）过程量
  判据是否与屋脊几何一致。
- **命令**：B1 产物第 8–15 列即 9 个指示因子；`make figures` 出 `indicators_vs_ridge.png`
- **解读**（V2.5+ 哨兵）：指示因子与屋脊的对应已固化为 **4 金值哨兵**——屋脊上升段
  O₃max 139.04 @ fnox 2.28 / 峰值 144.58 @ 3.42 / FNR 5.239 @ fnox 1.01 / OPE 饱和区谷
  3.832 @ fnox 2.28（`make batch` 现跑即校验）。逐指标二维矩阵与体制诊断表见
  `python experiments/run/run_B_surface.py`。
- **判据方向（以本 README 为准）**：FNR<1 → VOC-limited；FNR>2 → NOx-limited；1–2 过渡。

### B3. 响应面体制诊断（run_B）

- **命令**：`make batch`（81 案例）→ `python experiments/run/run_B_surface.py`（响应面 + 体制诊断）
- **产物**：`output/run_B/batch.csv` + 17 个 `indicator_matrix_*.csv` + `slices.csv`
- **解读**：行=fnox、列=fvoc 的逐指标二维矩阵可直接等值线出图；屋脊/FNR 过渡带/
  L_NOx-Q 穿越/OPE 谷在体制诊断表中汇总输出。

### B4. 自定义情景扩展

- **命令**：`input_v3/scenarios.config` 追加 `SCEN_MINE = E_NOx*0.7 E_CO*0.5 JFAC*1` →
  重编 runtime_v3 → `./fac_v3 SCEN_MINE …`
- **注意**：表达式仅乘法/加法；**物种→组映射硬编码于 `scen_build`**（新排放物种需改码）。

## C 类：求解器方法学实验

### C1. 双实现同题互证

- **命令**：`make functional`（FT1/FT3 自动）或手动双实现 diff CSV
- **解读**：KPP 纯机制 10.2 vs FAC 全物理 7.5 ppb——**差值本身是物理层贡献的测量**
  （原第三票 EBI 已退役，见版本精简记录；FAC vs EBI 差 18% 为历史实测）。

### C2. 线性解器对比（DENSE vs SPGMR）[历史——fac_spgmr 已退役]

- **结论保留**：实测 2.9× 加速、O₃ 差 1.3%；但无预条件 SPGMR 黎明步长崩溃（10b #12）——
  小机制 DENSE 反而稳，fac_v3 因此选 DENSE（教训固化为设计决策，标本已删）。

### C3. EBI 收敛阶 [历史——ebi_run 已退役]

- 实测 nit=2 vs nit=4 差 3.0%（<5%）→ CMAQ 惯例成立；FT4 检验随路线退役移除。

### C4. C-ABI 多实例线程安全验证

- **命令**：`fac_batch_omp 9 out_a.csv ../runtime_v3.dat` vs `OMP_NUM_THREADS=1` 复跑 diff
- **解读**：**并行/串行逐位一致**；四坑实录见 `cvode_c_mod.f90` 注释与 10b #13。

## D 类：诊断分析（双库函数级）

| # | 分析 | 函数 | 内容 |
|---|---|---|---|
| D1 | **指标全家族** | `indicators.route_indicators(prefix, fac)` | 12 指标窗口均值（FNR/H₂O₂:HNO₃/O₃:NOy/HCHO:NOy/O₃:NOz/NOy/L_NOx·Q⁻¹/P(H₂O₂):P(HNO₃)/OPE + 午后） |
| D2 | 指标时序 | `indicators.route_timeseries(prefix, fac)` | 同 12 指标的逐时序列（画日变化） |
| D3 | **IPR 六账积分** | `diagnostics.ipr_from_fac(prefix)` | 逐物种 P/L/E/D/B/X 梯形积分（任意 ns×6，路线无关） |
| D4 | CSV 读取 | `diagnostics.load(prefix)` | (h, C, irr_h, I)，自动适配任意物种/反应数 |
| D5 | Fortran 指数修复 | `diagnostics._fix_fortran_exp(tok)` | `8.27280-108` → `8.27280E-108`（gfortran 三位指数遗留格式） |
| D6 | 批量指标（C 内建） | `fac_batch_omp` 24 列 CSV | 见 §6.4；含运行期氮守恒不变量 |

> ⚠ **R12 勘误**：旧版此表列出的 `indicators(C)` / `families` / `chem_ledger` / `report`
> 已随诊断层重构删除。现行 API 见上表。

## E 类：机制工程实验

### E1. 元素守恒审计（机制质检）

- **命令**：`python audit/element_audit.py mechanism/mechanism_3_full_chain.fac` → **203/203**
- **解读**：旧机制 7 处历史 bug 标本（10a.1 全表；标本文件 `mechanism_0_expanded.fac`
  已删除，案例保留于 10a.1 表）——**机制开发审计流程活教材**。

### E2. ~~C 证人双审计~~（已删除）

原 `oz_audit.c`（NOx 望远镜闭合 + O₃ 收支）/ `cfg_audit.c`（序列单调）从未被
编译进任何构建。同款校验现由：`element_audit.py`（元素级）· 氮守恒在线不变量
（`fac_batch_omp` 运行期）· 六账闭合（`run_D` / `FT-IPR`）承担。

### E3. 机制三档切换（复杂度-成本标度）

- **命令**：`mcm_run`（29sp MCM 子集）→ `fac_v3`（111sp V2.5+）
- **解读**：fac_v3 111 物种 0.9 s——RHS 光滑性修复的实证（10b #11）；
  15sp 教学档与 30sp 旧扩展档均已删除。

### E4. MCM 真实机制城市案例

- **命令**：`cd mcm_city && ./mcm_run 5 ../output/mcm_city`
- **解读**：MCM 子集即便"只有甲烷"，无机组完整（NO₃/N₂O₅/HO₂NO₂/HONO）+ Troe——
  MCM"无机块完整+VOC 全显式"哲学实证。

### E5. 机制→C 跨语言编译

- **状态**：~~已退役~~（2026-10 清理）。`src/mechc_c.py`（.fac → C 速率函数）零消费者、
  不在构建链、无测试覆盖，已删除。当前唯一机制编译路径是 `facsimile_impl/mechc.py`。
- **历史**：如需查阅，见 git 历史（commit `1cc2777` 之前）。

## F 类：数值方法教学实验（10b 坑实录精选 8 个）

| # | 坑 | 教学主题 | 复现思路 |
|---|---|---|---|
| 4 | max(cosχ,0) kink 死锁 | BDF 与非光滑 RHS | 日循环换 sin² 对比步数 |
| 5 | 滴定死锁螺旋 | 光化学正反馈失稳 | 极高 NOx 初值跑 fac_v3（历史载体 fac_run 已退役） |
| 6 | 松容差重启放飞（MDA8 604 ppb） | 重启容差纪律 | 韧性重启段改容差 |
| 9 | EBI 日出滴定冲击 | 定点迭代收敛域 | nit=1 跑日出时段 |
| 11 | max(y,0) RHS 折点（**15× 步数差**） | RHS 光滑性决定步长 | fac_v3 源码对照 |
| 12 | 无预条件 SPGMR 黎明崩塌 | Krylov 预条件必要性 | fac_spgmr 黎明时段（该可执行体已退役，教训存续） |
| 13 | C-ABI 四连坑 | Fortran-C 互操作 | cvode_c_mod.f90 注释 |
| 14 | 未初始化 ios 静默回落 | Fortran 局部变量无默认零 | fac_v3 情景匹配段 |

## G 类：性能优化实验（6 项）

| # | 优化 | 基线 | 优化后 | 加速比 | 主题 |
|---|---|---|---|---|---|
| 1 | FAC→SPGMR | 0.60 s | 0.21 s | **2.9×** | 迭代 vs 稠密线性解 |
| 2 | Jacobian 延迟重组 | 0.60 s | 0.60 s | 1.0× | **负结果**（小矩阵无益） |
| 3 | EBI 自适应子步 | 0.35 s | 0.26 s | 1.4× | 自适应 vs 固定 |
| 4 | 进程并行扫描 | 29.0 s | 12.7 s | **2.3×** | embarrassingly parallel |
| 5 | C 代码生成 | — | — | ~2× | 编译后端选择 |
| 6 | KPP 批量 | — | — | 消 fork 开销 | 进程内循环 |

负载敏感项（OpenMP）load>ncpu 时测量无意义——本表数字均在空载下测得；
负载高时墙钟翻倍是假象，应看 CPU 时间。

## H 类：机制图集与教学

- **命令**：`python figures/gen_all.py`
- **产物**：9 版机制图（A 城市全天 / B 日间 / C 夜间 / D 乡村 / E 海洋 / F 城乡 /
  G 城海 / H 机制链条 / I FNR），Times New Roman Bold，330 DPI，PNG+PDF
- **配套**：`figures/README.md`（本地文件，逐图逐箭头教学手册 +
  **七机制链对比**：AtChem2/GEOS-Chem/CB6R3/MOZART/RADM2/FNR/ozone_lab）
  —— 注意 `figures/` 整个目录不入库（23 MB 二进制图），需在本地保留或按上述命令重生成
- **工程图**：`make figures` → response_surface_v3.png + indicators_vs_ridge.png

---

## 3. 架构总览

### 3.1 分层架构图

```
┌────────────────── 用户输入层 ──────────────────────────────────────────┐
│  input_v3/          (12 类配置, 人类可读, # 注释)                      │
│    initialConcentrations / meteo / constraints/                       │
│    bc / photolysis / photolysis-rates_MCM / emissions / dilution      │
│    drydep / wetdep / scenarios / model.parameters                     │
│    （机制 .fac 在 mechanism/，经 --mech 直引唯一真源）                │
└────────────────────────┬───────────────────────────────────────────────┘
                         │  src/cfgbuild.py (校验+归一, 全过才写出)
                         ▼
              runtime.dat / runtime_v3.dat (Fortran 顺序表控, 零文本解析)
                         │  modules/runtime_dat_mod.f90 (纯数据消费端)
    ┌────────────────────┼──────────────────────────┐
    ▼                    ▼                          ▼
 facsimile_impl/                 mcm_city/
 mechc.py 编译器                 MCM v3.3.1 子集 .fac
 CVODE BDF+DENSE                 CVODE-BDF+DENSE
 fac_v3 (情景+观测同化)          mcm_run (29sp 对照)
 fac_batch_omp (C-ABI+OpenMP)
 [已退役: fac_run/kpp_run/EBI/fac_expanded/fac_spgmr/kpp_batch]
    │
    ▼ 统一输出: <pref>_conc.csv + <pref>_irr.csv
    │
 diagnostics/  [纯函数诊断库]
   indicators.py (12 指标引擎) + diagnostics.py (load/ipr_from_fac) + csvfloat.py
    │
 experiments/run/*.py + figures/*.py  [业务层/可视化]
```

### 3.2 多语言角色矩阵

| 语言 | 模块 | 职责 | 为什么是它 |
|---|---|---|---|
| **Python** | `src/cfgbuild.py` | 配置编译（12 类→平面文件） | 文本校验/正则/单位换算/跨文件一致性 |
| **Python** | `diagnostics/` `experiments/run/` `tests/` | 诊断/实验/测试 | 数据分析生态 |
| **Fortran** | 三求解核 + runtime_dat_mod | 刚性 ODE 积分 + 物理层 | 编译性能 + CVODE/KPP 生态 |
| **C** | CVODE 库 | BDF/Newton 数值后端 | SUNDIALS 是 C 库 |
| **Shell** | Makefile + tests/smoke | 构建编排 | 通用构建语言 |

### 3.3 解耦原则

1. **文件即接口**：实现间零代码依赖（.fac + CSV 契约）
2. **生成码共享**：`oz_mech.f90` + `oz.reac/oz.prod` 计量表，现役三路线共用（fac_v3/batch 共享 gen_v3）
3. **诊断库纯函数式**：只消费 CSV
4. **KPP 已退役**（2026-10；`src/fac2kpp.py` 转换脚本保留作历史工具，无 make 目标）
5. **机制可插拔**：换 .fac 即换机制，`spidx` 按名索引
6. **与 Day xx 零代码耦合**（已全量审计）

### 3.4 外部依赖

| 依赖 | 位置 | 性质 |
|---|---|---|
| gfortran/gcc 13.2.0 | `~/gcc-13.2.0` | 编译器 |
| Python 3 + numpy/matplotlib | miniconda venv | 脚本运行时 |
| SUNDIALS CVODE 2.7.0 | `~/soft/atmenv/cvode` | 第三方 ODE 库 |
| ~~KPP 3.2.0~~ | `~/soft/src/KPP-3.2.0` | **已退役**（磁盘残留，构建不用） |

`~/soft` 部署与编译器同性质——本地工具链，不构成对 Day xx 工程的依赖。

---

## 4. 输入系统

### 4.1 文件格式总表

| # | 文件 | 语法 | 单位 | 缺省 | 示例 |
|---|---|---|---|---|---|
| 1 | `mechanism_*.fac` | FACSIMILE 子集（§4.3） | molec/cm³ / s⁻¹ | 必需 | `% KHO2NO : HO2 + NO = NO2 + OH ;` |
| 2 | `initialConcentrations.config` | `物种 值 [单位]` | molec/cm³ / ppb / ppm | 未列=0 | `O3 3.0e11` / `CO 80.0 ppb` |
| 3 | `meteo.config` | `NAME = 常数 \| file:`（7 量） | K / hPa / % / W·m⁻² / m / mm·h⁻¹ / - | 各有默认 | `BLH = file:constraints/met_BLH.txt` |
| 4 | `constraints/*.txt` | 两列 `t[s] 值`，严格递增 | 原生单位 | —— | `3600  1200` |
| 5 | `bc.config` | `TAU_MIX=s` + `物种=常数\|file:` | molec/cm³ | 不列=无 | `O3 = file:constraints/bc_O3.txt` |
| 6 | `photolysis.config` | `MODE=MCM\|DIRECT` + `LAT` + `MAPPING` | - | MCM | `MODE = MCM` |
| 6b | `photolysis-rates_*` | `j l m n name tau`（D 记法） | l: s⁻¹ | 必需(MCM) | `1 6.073D-05 1.743 0.474 J1 1` |
| 7 | `emissions.config` | `物种 基线 日型 振幅` | molec·cm⁻³·s⁻¹ | 不列=0 | `NO 3.0e6 BOX 0.9` |
| 8 | `dilution.config` | `DIL = 值` | s⁻¹ | 0 | `DIL = 0.0` |
| 9 | `drydep.config` | `物种 vd BLHFOLD` | vd: cm/s | 不列=不沉 | `O3 0.80 on` |
| 10 | `wetdep.config` | `物种 A B`（Λ=A·R^B） | s⁻¹/(mm/h)^B | 不列=无 | `HNO3 5.0e-5 0.8` |
| 11 | `scenarios.config` | `NAME = E_*乘 JFAC乘 DT加` | - | 可选 | `SCEN_HVOC = E_NOx*1.0 E_VOC*3.0` |
| 12 | `model.parameters` | 4 行：天/输出步/rtol/atol（**可改**；天数建议用 `--days` 覆盖） | - | 必需 | `5.0` |

### 4.2 配置编译器处理链

```
input_v3/ (12 类文本)
  │ ① .fac → 物种表 + 命名速率 + 反应行      ⑥ BC → 常数/序列
  │ ② 光解表 → l/m/n (D→E 记法)             ⑦ 排放 → 日型
  │ ③ photolysis.config → J 映射            ⑧ 沉降/稀释 → 速率参数
  │ ④ IC → 单位换算 + 存在性校验            ⑨ 情景 → 表达式静态解析
  │ ⑤ 气象 → 常数/序列(校验单调)            ⑩ model.parameters → 控制
  │ ⑪ 错误汇总 → 全部通过才写出
  ▼
runtime.dat (Fortran 顺序表控, 运行时零文本解析)
```

**编译期校验（拒绝能力实测）**：引用不存在物种 → 拒绝；序列非严格递增 → 拒绝；
情景表达式解析失败 → 拒绝；IC 单位不识别 → 拒绝。

### 4.3 机制文件格式（FACSIMILE 子集）

```
* 头部注释 (* 开头)
VARIABLE 物种1 物种2 … ;
KNAME = 表达式 ;                    ← 命名速率 (可引用前面的)
RO2 = CH3O2 ;                       ← 过氧池 (可选)
% 速率表达式 : 反应物 = 产物 ;       ← 反应行
```

| 元素 | 写法 | 生成 | 示例 |
|---|---|---|---|
| 科学计数 | `2.45D-12` | 保留 | Arrhenius A |
| 幂 | `x@n` | `x**n` | `(TEMP/300)@-2.6` |
| 命名速率 | `KNAME` | `q(n)` | `KHO2NO` |
| 光解 | `J<n>` | `JARR(n)` | `J<1>` |
| 环境 | `TEMP M O2 H2O` | 直通哑元 | `EXP(-1775/TEMP)` |
| 三体 | `O + O2 + M = O3` | M/O2 不进计量 | |

### 4.4 单位与换算完备性矩阵（九通道逐项，cfgbuild 单源 unit_factor）

| 输入通道 | 接受单位 | 缺省 | 校验行为 |
|---|---|---|---|
| **IC** | cm³ / **ppb / ppm**（每行后缀） | cm³ | 未知单位 → 编译期拒绝 |
| **BC 序列** file: | 头部 `# unit: cm3/ppb/ppm` | cm3 | 非单调/超 128 点 → 拒绝 |
| **BC 内联常数** | cm³ / **ppb / ppm**（后缀，V2.3 新增） | cm3 | errs 门拒绝 |
| **排放** | molec/cm³/s / **ppb/h**（头部 `# unit:` 指令，V2.3 新增） | molec/cm³/s | 负值/非法日型/振幅越界 → 拒绝 |
| **气象 7 量** | K / hPa / % / W·m⁻² / m / mm·h⁻¹ / 无量纲 | 各有缺省 | — |
| **干沉降 vd** | cm/s | 不列=不沉 | — |
| **湿沉降** | a: s⁻¹/(mm·h⁻¹)^b | 不列=无 | — |
| **光解 JMAP** | D 值: s⁻¹（正午） | — | 槽位缺失 → 拒绝 |
| **model.parameters** | 天 / 输出步长[s] / rtol / atol | 非法值编译期拒绝 | 必需；天数可用 `--days`/`DAYS=` 临时覆盖 |

`unit_factor` 单源函数（read_series/read_ic/BC 三处共用）——新增单位只改一处。
**完备性结论**：V2.3 补齐 BC 内联后缀与排放 ppb/h 指令后，用户数据可用原始单位
（ppb / ppb/h）直接落位，"必须预先手工换算"的场景清零。

**四类全量约定**（机制速率 molec/cm³/s / Fortran 硬编码清单 / 内部计算量 /
⚠ PPB 因子固定 1atm·298.15K 的 ~1–2%/5K 口径差）：
- 摘要版：[`OPERATIONS_MANUAL.md §2b`](docs/OPERATIONS_MANUAL.md)
- **完整版（逐子量 + 换算公式 + 硬编码登记）：[`UNITS.md`](UNITS.md)**

---

## 5. 求解器实现（CVODE 双驱动现役；KPP/EBI 已退役）

### 5.1 对比总表（KPP/EBI 两列为历史档案）

> 现役 = FACSIMILE 机制的两种 CVODE 驱动（① fcmix 单箱 `fac_v3` / ② C-ABI+OpenMP 批量
> `fac_batch_omp`），机制/物理/输入三样全共享。性能列为教学机制时代实测，仅供量级参考。

| | FACSIMILE（现役） | KPP [已退役] | EBI [已退役] |
|---|---|---|---|
| 编译器 | mechc.py | 真 KPP 3.2.0 | 无（手写核） |
| 算法 | CVODE BDF+Newton+DENSE | Rosenbrock+稀疏 LU | 分组定点迭代 |
| 物理层 | 全套内嵌 resid | **纯机制** | 排放+干沉+背景 |
| 自适应步长 | ✓ | ✓ | ✗（12 子步） |
| Jacobian | 解析（FCVDENSESETJAC） | 符号微分自动 | 无 |
| 性能（历史实测） | 0.93 s/run | **0.21 s/run** ★ | 0.87 s/run |
| 精度 | 高（rtol 1e-6） | 高 | 中（nit=2 收敛 3%） |
| 工业对应 | AtChem2 | GEOS-Chem | CMAQ gas |
| 定位 | 业务主力 | 机制基准 | 对照路线 |

### 5.2 FACSIMILE 路线

```
facsimile_impl/
├── mechc.py            FACSIMILE→Fortran 编译器
├── oz_fac_v3.f90       工程版驱动 (111sp + runtime_v3 全通路 + 六本账)
├── fac_batch_omp.f90 + cvode_c_mod.f90   OpenMP 批量 (C-ABI)
├── gen_v3/             生成的 oz_mech.f90 + oz.reac/oz.prod
└── fac_v3 fac_batch_omp
```

核心流程：mechc.py 解析 .fac → `update_rates` + 计量表 → resid 通用稀疏累积 →
CVODE fcmix BDF+Newton → 解析 Jacobian → 拐点韧性重启（FCVREINIT 同容差+小步长）。

### 5.3 KPP 路线（build_kpp.sh 五补丁链）

1. `#DEFVAR` 物种须给元素组成（`O3 = O + O + O;`）2. eqn `^` → `**`
3. `SPARSE_LU_ROW` 生成码必须 `-cpp` 4. 光解 J → Rates 模块全局变量
5. `#INITVALUES` 保 FIX 槽位

### 5.4 EBI 路线 [已退役 2026-10-25]

原为第三对照路线（O/O¹D 闭式稳态 + 快慢两组定点迭代 nit=2 + 12 子步防护），
因无输送耦合计划且互证价值由 KPP 票承担而退役删除——本节保留其设计要点供考古。

---

## 6. 诊断库（双库架构）

诊断层由**两个库**组成，职责分离：

### 6.1 `diagnostics/indicators.py` —— 路线无关指标引擎（主力）

```python
route_indicators(prefix, fac_path)      # → dict: 12 个指标 (窗口均值)
route_timeseries(prefix, fac_path)      # → dict: 12 个指标 (逐时序列)
Mechanism(fac_path)                     # 机制解析器 (.fac → 物种/NOy权/反应表)
```

**特点**：
- **路线无关**：只吃 `<prefix>_conc.csv` + `<prefix>_irr.csv` + 机制 `.fac`
- **NOy 自动判定**：从 `.fac` 的 `FORMULA` 段取含 N 物种（N₂O₅ 计双 N）——新 N 物种零改动接入
- **指标全集**：FNR / H₂O₂:HNO₃ / O₃:NOy / HCHO:NOy / O₃:NOz / NOy / L_NOx·Q⁻¹ / P(H₂O₂):P(HNO₃) / OPE + 午后口径

### 6.2 `diagnostics/diagnostics.py` —— CSV 读取与六账积分（底层）

```python
load(prefix)                # → (h, C(T,ns), irr_h, I(T,nr))  自动适配任意 ns/nr
ipr_from_fac(prefix)        # → (species, {sp: {P/L/E/D/B/X: 积分}})  六账梯形积分
IPR_PROC                    # ["P","L","E","D","B","X"] 六账列序
_fix_fortran_exp(tok)       # gfortran 三位指数(8.27-108) → 标准记法
```

**注意**：`load` 返回的 `C` 已剥去 `h` 列，因此按表头名索引时需 `chdr.index(sp) - 1`
（`experiments/ipr_closure_check.py` 与 `ipr_irr5.py` 均按此约定）。

> ⚠ **R12 勘误**：本节旧版曾列出 `indicators(C)` / `irr_integrate` / `irr_key_processes` /
> `families` / `chem_ledger` / `report` 六个函数——**它们已在诊断层重构中删除**，
> 照抄会报 `AttributeError`。现行 API 以上两节为准。

### 6.3 指示因子阈值（以本 README 为准）

| 因子 | 判据 | 出处 |
|---|---|---|
| FNR=HCHO/NO₂ | **<1 → VOC-limited；>2 → NOx-limited；1–2 过渡** | Sillman 1995 / Tonnesen-Dennis 2000 |
| H₂O₂/HNO₃ | >0.2 → NOx 控制区 | Sillman 1995 |
| L_NOx/Q | 0.5 | Kleinman 2005 |
| OPE=P(O₃)/P(NOz) | NOx 饱和区触底 | Kleinman 2005 |

### 6.4 批量侧（C 代码内建，非 Python）

`fac_batch_omp` 的 24 列指示因子由 Fortran 内建预算表计算（启动期索引化，
`batch_budget_init`）：Q×6 + L_NOx×2 + 产率 41 项，NOy **17 物种**。
提取用防脆断原语 `prod_rate` / `rxn_rate`（按产物系数 / 反应物名多重集匹配，机制改行不改码）。

---

## 7. 测试体系（计数为实跑值）

### 7.1 冒烟 —— 1/1 PASS

FAC-V3 单项（判据：log 含 ok + O₃ 末列非零）。KPP 冒烟已随 KPP 路线退役移除。

### 7.2 功能 —— 62/62 PASS（12 族 FT-\*）

FT-CFG（配置编译 + runtime_v3.dat 逐字节比对）/ CFGSPEC（配置规格负测试）/
DAYS（`--days` 覆盖不改文件）/ **FREEZE（冻结版规格：.fac↔README↔评估报告三向一致，111sp/203rx）**/
GIT（output/ 不入库、.fac 全跟踪）/ **MECHK（mechc 计量表防御：行内 `!` 注释剥离 + 未知 token 硬失败）**/
MP（多前缀）/ PHOT（PHOT 头行 LAT/DOY 生效）/ SAT（sat2surf 五断言）/
SYNC（物种同步：死端/NOy 表/排放分组/ATOL 覆盖）/ TPL（模板完整）/ UNIT（单位换算）；
另含**批量自检段**（batch CSV 契约 / 屋脊存在 / FNR 单调 / L_NOx-Q 值域）与
**回归哨兵**（4 金值 ±10–15%：屋脊上升段 O₃max 139.04 / 峰值 144.58 / FNR 5.239 /
OPE 谷 3.832——机制/求解器改动造成的数值漂移被 make test 直接拦截）。

### 7.3 `make test` 全链

`functional(62/62, 含批量契约与体制诊断断言) → element_audit(203/203) → mcm_compare(8.9%)`

---

## 8. 运行结果汇总

| 实验 | 结果 | 出处 |
|---|---|---|
| FAC 基线 O₃（教学机制） | 7.5 ppb | FT2b |
| KPP 纯机制 | 10.2 ppb | FT3 |
| EBI | 6.2 ppb（vs FAC 18%） | FT3 |
| 教学机制五情景 MDA8 | BAU 8.0 / LNNOX 43.7 / HVOC 15.9 / COLD 8.7 / CLOUD 7.7 | §2-B3 同型 |
| **V3 六情景 O₃max（V2.5+，2026-10-27 run_A）** | **BAU 99.3 / LNNOX 49.4 / HVOC 95.0 / HOT 113.4 / CLOUD 80.7 / ISOPX 64.5**（末 24h 均值 58.4/30.3/55.4/65.3/47.3/38.4；V2.3 口径为 101.6/53.9/111.3/112.6/84.3/73.1，方向结论不变） | 实验 A |
| **屋脊（V2.5+ 哨兵）** | 屋脊上升段 O₃max 139.0 @ fnox 2.28；峰值 144.6 @ 3.42（高 NOx 端回撤） | batch 金值哨兵 |
| **屋脊（V2.3 历史口径）** | fnox=2.28 处 O₃max 125.5 ppb（回撤 90%） | batch_v3_indicators |
| **指示因子-屋脊（V2.4+）** | 屋脊 fnox=2.28；L_NOx/Q 穿 0.5 于 1.52→2.28；OPE 谷 2.28 | run_B / make batch |
| **IPR 闭合** | 各物种中位 1.004–1.007（内步梯形） | ipr_closure_check |
| **IPR↔IRR 一致** | O₃ 逐小时最大相对差 0.000% | hourly_pa |
| **OSSE 技能（实验 F）** | 地面同化 HCHO RMSE +27%；卫星反演 E_NOx×1.34 → NO₂ RMSE +74% | obs_osse |
| **跨求解器双胞胎** | fac_v3 vs batch 六账占比逐项一致；IPR↔IRR 0.000% | ipr_irr5 |
| fac_v3 性能 | ≈2.6 s / 4802 步 / 零重启（5 天 111 物种 V2.5+，run_A 实测） | 实验 A |
| 屋脊（V2.1 教学机制历史值） | fnox≈1.5 处 29→12.6 ppb | batch_v3 |
| 指示因子-屋脊一致（V2.1 历史值） | FNR 穿 1–2 带≈屋脊位；L_NOx/Q 穿 0.5≈fnox 1.0；OPE 谷 13.4 | indicators_vs_ridge |
| ~~fac_v3 性能~~（V2.1 历史值） | 0.19 s / 4965 步 / 零重启（5 天 **53 物种**，已被 111 物种现版取代） | 10a.3 |

---

## 9. 性能分析与优化

### 9.1 瓶颈（教学机制 fac_run）

| 瓶颈 | 占比 |
|---|---|
| 稠密 LU (15×15) | ~40% |
| Jacobian 组装 | ~30% |
| 拐点韧性重启 | ~20% |

### 9.2 优化状态

| # | 方向 | 状态 |
|---|---|---|
| 1 SPGMR | **已验证并退役**（实测 2.9×；fac_v3 单箱选 DENSE——黎明崩塌教训固化为决策） |
| 2 Jacobian 延迟 | **已验证无益**（1.0×，负结果保留） |
| 3 EBI 自适应子步 | **已实现**（1.4×） |
| 4 并行扫描 | **已实现**（scan_parallel 2.3× + fac_batch_omp 2.05×） |
| 5 C 代码生成 | ~~已实现~~（mechc_c.py，2026-10 删除：零消费者） |
| 6 KPP 批量 | 已由 fac_batch_omp 线程级取代（kpp_batch 退役） |

---

## 10. 功能拓展方案

| 拓展 | 状态 | 说明 |
|---|---|---|
| A1 PAN / A2 ISOP / A3 SO₄ / A4 卤素 | **已落地** | 扩展机制 30sp/51rx + mechanism_3 111sp/203rx（isoprene 完整二代链 + V2.3/V2.4 补链） |
| B2 在线光解 | **已落地**（v3 SZA 参数化） | Fast-J 级需接 Day 26 辐射柱 |
| B5 响应面并行 | **已落地** | fac_batch_omp |
| B1 垂直多层 | 未做 | 3-5 层箱 + 层间扩散 |
| B3 浓度约束 | 未做 | 任意物种观测序列约束（机制检验） |
| B4 伴随灵敏度 ∂O₃/∂kᵢ | 未做 | 不确定度传播 |
| B6 netCDF I/O | **已落地** | `make nc`（`modules/nc_output.py`，CSV→NETCDF4；需 conc+irr 配对） |

---

## 10a. 工程落地版 V3（源码级）

### 10a.1 mechanism_3 六链条与元素审计

111 物种 203 反应（C/H/N/O/S/Cl/Br 精确配平）。十段链条见 §1.4。旧机制 7 处科学错误
由元素审计器抓出并修正：

| 错误 | 影响 | 修正 |
|---|---|---|
| O1D+H2O / O1D+M / RO+O2 速率未乘 H2O/M/O2 | O1D 损耗慢 9 个量级，双通道权重错 ~8 倍，OH 源高估 | 速率表达式并入 *H2O/*M/*O2 |
| CH3CHO 光解产物 CH3CO3+HO2+CO | C 超配（2C→3C） | CH3O2+CO+HO2 |
| CH3CO3+NO / +HO2 / +CH3O2 产物缺 CO2 | C 不守恒 | 补 CO2（酰基分解） |
| O3+NO2→NO | 漏 NO3 生成主通道 | →NO3(+O2) |
| Cl+HCHO 5e-13 | 慢 150 倍 | 8.1e-11·exp(−34/T) |
| ClO+NO2→Cl+NO3 / Br+O3 重复行 / SO2 自光解 | Cl 循环失真等 | ClONO2 储库；删除 |

### 10a.2 runtime.dat 全通路

`input_v3/` → cfgbuild → `runtime_v3.dat` → `runtime_dat_mod.f90`（318 行纯数据模块）
全消费：真光解（SZA+MCM l/m/n+JFAC）、时变气象（插值外延端值）、**BLH 膨胀稀释**
（中心差分 dlnH/dt）、排放日型（BOX sin²/FLAT）、干沉（vd/BLH）湿沉（a·R^b）、
BC 松弛、情景因子（排放组×JFAC×ΔT）。顺带修复 cfgbuild 两处存量 bug
（多行 VARIABLE 截断、NO₂ 光解映射错到 O₃→O¹D 参数差 ~200×）。

### 10a.3 fac_v3 数值要点

RHS 去掉 `max(y,0)` 折点（39727→4965 步，**15×**）；单箱 DENSE（无预条件 SPGMR
黎明崩塌）；ATOL 按名调优表（§1.1④）；韧性重启同容差 + MAX_STEP 100→10 s。

### 10a.4 OpenMP 批量（C-ABI）

fcmix COMMON 块限制单进程 1 问题 → `cvode_c_mod.f90` 直连 C API，每线程独立
mem/N_Vector/user_data；81 案例并行/串行**逐位一致**，4 线程 2.05×。
四坑（RHS 返回 int 函数 / N_Vector 二次解引 / CVodeFree(void**) c_loc(mem) /
SetArrayPointer 不清 own_data）见 `cvode_c_mod.f90` 注释。

### 10a.5 自检轮次 3/4（指示因子升级后的结构性预防）

**轮次 3 抓出 4 项**：

| 问题 | 定位方式 | 修复 |
|---|---|---|
| **OPE 的 P(NOz) 双计 HNO₃**（真实数值 bug：旧 OPE 系统性低估 ~27%） | 预算表重构后"数值不变性"对比暴露 → 双算法同跑分解 | 预算表正确枚举（14 项）；OPE 修正后 fvoc=1 切面 148.9→52.1→**13.4（饱和区谷）**→回升 |
| 诊断热路径每次调用做名称多重集匹配 | 编译/剖析 | `batch_budget_init` 启动期解析为索引表；自检改查表并 stop 1 硬失败 |
| 两张关键图内联生成不可复现 | 自检 | 落盘 `figures/fig_indicators_ridge.py` / `fig_response_v3.py` |
| functional_test 汇总口径 | 跑测 | 汇总移末尾；**新增批量自检段**（6×6=36 案例） |

**轮次 4（"能优化解决吗"——结构性预防 + 微优化）**：

1. **运行期氮守恒不变量**（`Σ w_N·ydot_chem ≈ 0`，每线程首快照，阈值 1e-4）——
   element_audit 的动态版，上线即抓出**两个真实 NOy 定义缺陷**：ClONO₂/BrONO₂ 与
   ISOPNO₃ 未计入 NOy（后者曾把 O₃/NOy 与 OPE 悄悄带偏）。NOy 现为 **13 物种**清单。
2. **诊断热路径索引化**（实施中被自己的全局替换咬了两口——自赋值与守恒块误插声明区，
   均被编译警告/越界垃圾当场暴露；教训：字符串替换后必查自赋值）。
3. **回归哨兵**：4 个金值格点 ±10–15%（屋脊峰 29.0 / 饱和区 12.55 / FNR 转折 2.109 /
   OPE 谷 13.35——双计 bug 教训位）。
4. `make figures` 可复现闭环。诚实边界：SPGMR+预条件对 111 物种仍不值得（DENSE 每步
   LU ~4e5 flops 可忽略，实测 SPGMR 黎明崩塌）；全文件金值 CSV 刻意不用（脆断）；
   负载下 wall-time 不作性能结论。

### 10a.6 指示因子全家族 × 屋脊交叉验证

24 列 CSV（§1.1⑤）。NOy=**13 含氮物种**加权（N₂O₅ 计双 N；含 ClONO₂/BrONO₂/ISOPNO₃——
氮守恒不变量抓出补全）。fvoc=1 切面：FNR 5.2→0.2 穿阈值带（1–2）与屋脊（fnox≈1.5）重合；
L_NOx/Q 穿 0.5 于 fnox≈1.0；OPE 谷 13.4——与 Sillman 1995 / Kleinman 2005 文献结论一致。
图由 `figures/fig_indicators_ridge.py` 可复现。

### 10a.7 响应面与体制诊断

fvoc=1 切面经典 **O₃–NOx 屋脊**：29 ppb（fnox≈1.5）→12.6；FNR 指示场与屋脊精确对应
（转折带 FNR≈1.5–2）。低 NOx 角 FNR 4–37（NOx-limited），高 NOx 角 0.01（VOC-limited）。
`output/batch_v3.csv` 供诊断库直接消费；图由 `figures/fig_response_v3.py` 可复现。


## 10c. 全维度优化（P0–P2 批次执行记录，2026-10-25）

基于五路线源码诊断报告执行的全部优化项（不加路线、不动定位）：

### 功能
- **F1 IPR 五路线全通**：`<pref>_ipr.csv` 六过程账（`SPECIES_P/L/E/D/B/X`：机制产/耗/
  排放/干湿沉/BC+稀释/BLH 膨胀）——fac_v3 319 列、fac_run/kpp 91 列、mcm 175 列、
  batch 323 列（末 24h 窗口均值，`<outf>.csv_ipr.csv`）。夜间 O₃ 收支自洽抽查通过。
- **F2 情景批量维度**：batch CLI 第 4/5 参 = jfac/ΔT 列表（如 `0.6,1.0` × `0,5`），
  网格 ng²×nj×nd——SCEN_CLOUD/HOT 的响应面版。
- **F3 序列 BC 接线**：bc.config `file:` 序列经 chem_phys_mod 生效（bc_s 优先于常量）。
- **F5 午后口径**：batch 新增 6 个 `*_pm` 列（12–18h 窗口，Sillman 1995 原文口径）。
- **F4**：`make nc`（CSV→netCDF4 打包）；**P1**：响应面由 `make batch` + `run_B` 承担（scan_suite 已删除）
  （111sp+指示因子，输出体制诊断表：屋脊/FNR 过渡/L_NOx-Q 穿越/OPE 谷）。

### 机制
- **M1 15sp 速率承载修复**（与 mechanism_3 同源三处：O1D+H2O/M、CH3O+O2 并入
  H2O/M/O2）：新基线 O₃ 29.1 ppb（旧 7.5 为 OH 过产 bug 伪象）；KPP 票仅微动
  （10.24 vs 10.2——KPP 的 H₂O/O₂ 一直在计量式里，病只在 FAC 隐式承载）——
  双路线交叉验证修复正确性。FT10 断言随物理翻转（修正后 CH₄-only 玩具为
  NOx-limited，旧"高NOx低O₃"是 bug 伪象）；fac_run 同步移除 max(y,0) 折点。
- **M2 ISOP+O₃**（旧第 97 反应，V2.3 已改为 Criegee OH 产率版 0.26OH+CH₃CHO）：历史基线 0.5MVK+0.5MACR+HCHO+0.5O₂，
  元素精确配平；OH 产率省略已文档化。哨兵漂移 <0.3%。
- **M3 分子式入机制**：.fac 新增 FORMULA 段（机制自带配平说明书，单一事实源），
  mechc/fac2kpp/audit 三解析器支持，审计器双字符元素（Cl/Br）正则解析。

### 模块
- **A1+A2 物理层/情景统一**：`modules/chem_phys_mod.f90` 物理五项唯一实现
  （v3 双胞胎共享），情景语义统一为 esc+dt+jfac 传参——OPE 类双实现 bug 温床根除；
  重构后 O₃max/FNR/哨兵逐位不变（NOy 类比值第 3 位小数微差 = FP 求和顺序，良性）。

### 性能
- **mcm_run 复活**：重生成 gen（旧 gen 用固定 nj=7 的 mechc，J<41>/<51> 越界读垃圾
  → 黎明 mxstep 崩溃——存量缺陷非本轮引入）；DENSE 替换无预条件 SPGMR +
  max(y,0) 移除：5 天模拟 79s 失败 → **0.11 s 完成**（~700×）；光解换 MCM l/m/n
  表+SZA（M4，替 dayf×100 hack）。
- **mechc 能力扩展**：RO2 池单成员展开（→Y(idx)）、N2 常数（N2C）、nj 自动、
  MCM 真实机制全支持（越界修复的根因）。

### 载体修复
- gfortran 3 位指数省略 E 的遗留格式（`8.27280-108`）→ diagnostics/nc_output
  健壮解析。


## 10d. 指示因子屋脊实验（ridge5，2026-10-25；**历史记录**：当时五路线口径，现役三路线）

`experiments/ridge5.py`（图 `experiments/fig_ridge5.py` → `figures/ridge5_routes.png`）：
f_NOx ∈ [0.03, 31.6] 对数 15 点 × fvoc=1，四条可实验路线（fac_v3 与 batch 同机制不重复），
指示因子由**路线无关引擎** `diagnostics/indicators.py` 从 CSV 契约 + FORMULA 段（M3 收益：
NOy 物种自动判定）后处理——batch 之外的路线零 Fortran 改动即可实验。

| 路线 | 情景空间 | O₃屋脊 fnox* | 最一致指示因子 (偏差) |
|---|---|---|---|
| batch 排放 | 排放 | 3.42→2.28 | **FNR 过渡贴屋脊（V2.4 口径）** ✓ |
| mcm_run 29sp MCM 排放 | 排放 | 0.14 | **L_NOx/Q (0%)** — Kleinman 判据精确命中 |
| fac_run 15sp 排放 | 排放 | 1.64 | 无 (<0.5 穿越不存在；浓度比过渡在 0.07，偏离 96%) |
| kpp_run 15sp IC (EKMA) | 初始浓度 | >31.6 (边界*) | 无 (纯化学无汇 → 单调) |

**实验结论**（与文献对照）：
1. Sillman 指示因子-屋脊对应关系**依赖机制与配置**：含完整 VOC 链的机制（V2.3 74sp）成立
   （O₃/NOy 族 12%），真实 MCM 速率下 Kleinman 自由基判据最准（0%），而 CH₄-only
   玩具体系全面失效——呼应 Souri et al. 2020 "单一阈值不可移植"。
2. **纯化学无汇体系没有屋脊**（kpp IC 空间 fnox 到 31.6 仍单调）——EKMA 屋脊的
   存在依赖滴定/沉降等汇过程，这是 kpp 路线重定义为 IC 空间 EKMA 后的方法学产出。
3. nan = 该指示因子在网格内未穿越判据（如 15sp 的 L_NOx/Q 始终 <0.5——CH₄-only
   体系 NOx 永不饱和自由基），非错误；带 `*` = 网格边界伪迹已剔除。

**实验过程中修复的错误/缺失**（记录时口径为五路线）：
- `indicators.py` 引擎两处：`%` 分支误用带 `;` 原文（rhs 尾分号毁全部名字匹配）；
  Q 匹配误用速率表达式 token（FACSIMILE 反应物在计量侧）→ 改 lhs/rhs 多重集
  （与 Fortran 预算表同语义）。
- `mcm_run`：**无情景参数**（缺失）→ +fnox/fvoc CLI；缩放先于 read_stoich 初始化
  （顺序 bug，静默丢失）→ 移后；conc **缺 H₂O₂/HNO₃/NO₃/N₂O₅/HONO 列**（指示因子
  无米下锅）→ 扩至 15 列；irr 格式 off-by-one（71 反应写 70 描述符 → 空尾字段）。
- `kpp_run`：snox/svoc 为死参数（物理层移除后无效应）→ 重定义为 **IC 标度 EKMA**
  （Dodge 1977 浓度空间正统）；IPR 的 E/D 列曾显示从未施加的排放（假账）→ 纯机制
  语义清零；bcv/bcr 死变量清除。


## 10e. 输入模板体系 + 操作手册（2026-10-26）

- **`docs/OPERATIONS_MANUAL.md`**：架构与速查——构建/一键验证、路线分工与编译边界
  （§0.5 显微镜vs望远镜/三铁律）、路线选择速查（§0.6）、两条事实源与改机制闭环、
  五类外数据模板接入、实验 A–I 速查、已知边界、故障排查、自纠自查轮记录附录。
- **`docs/EXPERIMENTS_GUIDE.md`**（手把手级实验指南）：
  环境与首次构建（4 勾选项清单）→ **化学基础（三循环+十段链条+化学时钟+物种分类）**
  → **数据如何进程序（12 类配置逐一详解 + 排放单位换算三种推导 + 落位模板 +
  9 条硬约束）** → 路线决策树与编译边界 → **实验 A–I 逐个八段式**（科学问题/
  前置/逐步命令/实测输出/逐列解读/判读要点/变异玩法/常见坑）→
  **自定义数据接入 Cookbook（6 配方含完整 bash）** → 输出契约（真实样例行）→
  **指标解读手册（含 L_NOx/Q 与 FNR 冲突案例的归因分析）** → 质检四关卡 →
  故障排查（4 类 30 条）→ 性能预算 → **Python API 速查 + 二次开发片段** →
  **最小示例数据集（3 天配置 + 24 点观测 + 验证清单）**。
- **四份文档分工**：README（项目是什么）→ docs/OPERATIONS_MANUAL（架构+速查）→
  docs/EXPERIMENTS_GUIDE（实验怎么做，手把手）→ figures/README（机制图集+七机制对比）。
  环境资料（ENVIRONMENT + check_env.sh）在 `environment/`。
- **`templates/` 五类输入模板**（格式说明内嵌文件头）：
  ① `emissions_inventory.template`（排放清单：物种/molec·cm⁻³·s⁻¹/BOX·FLAT 日型）；
  ② `surface_o3_station.template`（地面 O₃ 站点：ppb 序列 → bc 张弛 nudging）；
  ③ `surface_aqi_station.template`（**常规 AQI 监测站点**：HJ 193 多污染物宽表
  μg/m³ → 拆分转换表 → ppb；O₃/NO₂/SO₂/CO 可接入，PM 不接入，NO 为差折算不接入，
  NO₂ 建议弱 τ 或走排放反演）；
  ④ `surface_voc_station.template`（地面 VOCs 组份：一站一物种 + PAMS-57 映射注记——
  HCHO/ISOP/MVK/MACR/CH₃CHO 直接入，苯系物无芳烃化学不接入）；
  ⑤ `satellite_column.template`（卫星柱浓度：TROPOMI 存档格式 + 离线转换约定——
  柱/空气柱×混合层分数 → 表面等效 ppb；0-D 箱体无气柱概念，引擎不直读柱）。
- **cfgbuild 序列单位指令**：约束序列文件 `# unit: cm3|ppb|ppm`（首 token 容忍行内
  注释）→ 编译期换算 molec/cm³；全链验证：模板站点序列（30→108 ppb）经 bc.config
  → runtime.dat → fac_v3，O₃ 被张弛拉至 39 ppb（无约束 17.9）✓。
- **IPR 精确累积收尾**：fac_v3 的 `<pref>_ipr.csv` 升级为 **ONE_STEP 内步梯形
  精确时间积分**（每行 = 小时窗积分，molec/cm³），消除小时快照对大项对消预算的
  混叠；分解完备性检验 `experiments/ipr_closure_check.py` 全物种比值 **1.01–1.02**
  （均匀 +2% 为内步梯形与积分器内部求积口径差）。快照路线（fac_run/kpp/mcm）
  保留小时快照口径（混叠局限已文档化）。


## 10f. 观测约束 OSSE 实验（2026-10-26）

`experiments/obs_osse.py`——双胞胎观测系统模拟实验：真值运行（NOx/ISOP 排放 ×1.4）
采样加噪 → 三类观测文件（地面 O₃ 站点 / VOCs 组份站点 / 卫星转化 NO₂）→ 基础配置
约束运行 → 技能指标。三条同化通路按物理语义分立：

- **站点 nudge（O₃/HCHO/ISOP）**：τ ∈ {2,6,24,72h} 扫描——**全部劣于自由运行**
  （O₃ RMSE 1.25 → 最优 2.10；HCHO 0.08 → 0.30），且强度越大越差。机理：O₃/NO₂
  受光稳态支配，直接 nudge 被动物种使其与前体物脱耦——约束对抗化学。旧设计
  （卫星 NO₂ 直接强 nudge）曾致 h→1e-79 崩溃（PSS 破坏），已按正确通路重构。
- **卫星 NO₂ → 自上而下排放反演**：E_NOx ×= mean(obs)/mean(free)——NO₂ RMSE
  **0.20 → 0.04 (+82%)**，且 O₃ 连带改善（1.25 → 1.07）：反演修正成因（排放），
  效益外溢到未同化物种。
- **IPR B 账同化倾向**：τ=6h 约束运行 O₃ B=+4.4e6/s（占过程总账 2%）、HCHO 13%、
  ISOP 1%——同化强度的过程级度量。首窗垃圾值缺陷（ipr_w 未初始化）已修复。

**方法论结论**：约束必须走模型的因果结构——排放入口（自上而下）改善全局技能；
被动物种的直接 nudge 是症状级约束，劣化技能。这与同化文献（nudging vs flux
adjustment）一致，在 ozone_lab 得到定量实证。

**cfgbuild 光解完备性校验**（同轮加固）：机制 J<> 槽位缺参数表时编译期拒绝
（此前静默零光解——测试目录漏拷表时曾实际发生）。


## 10g. 逐小时过程分析（hourly PA，2026-10-26）

`experiments/hourly_pa.py`（图 `figures/hourly_pa.png`）——fac_v3 的 CMAQ-PA 式
逐小时过程分析，数据为 ONE_STEP 内步梯形**精确小时窗积分**（IPR 74×6 + IRR 131 反应，
两独立累积路径）。**本轮修复两项真实缺陷**：

1. **IPR 首诊断行含 96h 累计**（非小时窗）——`ipr_w` 窗口累计原仅在写出步更新，
   首诊断行 (h=96) 混入 0→96h 全部积分（O₃ +1.3e11，比正常窗大 6–60 倍；完备性
   检验曾因 six_a 切片恰好跳过该行而漏检）。修复：窗口累计每步推进。
2. **IRR 逐小时积分缺失**——IRR 原为小时边界快照（J 值在边界取值有偏），现与
   IPR 同法内步梯形积分，反应级逐小时预算无偏。

**完备性结论**：修复后 IPR↔IRR 两独立累积路径逐小时互证 **最大相对差 0.000%**
（O₃ 的 P 与 L 全 24 窗）；快照式路线（fac_run/kpp/mcm）的混叠局限文档化。

**过程分析产出**（物理图景自洽）：夜间 O₃ 预算沉降占 68–90%（滴定仅 3%）；
白天光稳态 P≈L≈50%（净产为小差值）；晨昏过渡（h05/h20）清晰可辨；夜间反应级
预算可见卤素通道（r53 Cl+O₃）。


## 10h. 排放+观测灵敏性/伴随分析（sens_adj，2026-10-26）

`experiments/sens_adj.py`——OAT 排放扰动（9 物种 ±10%，config 级）→ 观测雅可比
（4 物种×6 时刻 24 维）→ 截断 SVD 线性化反演 vs OSSE 真值；跨路线斜率一致性；
真速率灵敏度（.fac 拷贝→mechc→重建→±5% 中心差分，**取代 rate_sensitivity.py
的硬编码假实现**——env 变量从未消费的占位输出）。

**发现**：
1. **雅可比结构**：NO 列范数 8.15 支配（NOx 主导 O₃），SO₂/DMS 0.15/0.19（硫化学
   不触 O₃）✓；
2. **可辨识性受限**：观测雅可比条件数 1311，截断后仅 1/9 模态可反演——NO 恢复
   57% (0.23/0.4)，ISOP/NO₂ 与 NO 响应近共线不可分。**排放反演需要更丰富的观测
   向量**（物种×时刻多样性）或先验正则化；
3. **跨路线斜率**：fac_v3 OAT 局部 (±10%) 6.4 vs batch 网格弦 (0.675→1.0125) 8.7
   ppb/unit——27% 偏差符合 O₃–NOx 曲线凸性（弦>局部），非错误；
4. **真速率灵敏度谱**：J<1> 光解 +0.283 / KCH4OH +0.206 / KHO2NO +0.119 /
   KO3NO −0.138（唯一负值，滴定）/ KCOOH +0.028——符号量级全物理合理；
5. **存量确认**：SO₂ 排放不在任何情景缩放组（分组遗漏，schema 级修复待定）；
   排放组映射双实现（runtime_dat_mod/batch_mod）当前一致（一致性风险已记录）。


## 10i. 情景矩阵实验（scen5，2026-10-26；**历史记录**：标题沿用旧口径）

`experiments/scen5.py`——8 情景 × 5 路线的 O₃max 响应矩阵（`output/scen5/`）。
情景集：BAU/LNNOX(0.3)/HVOC(×3)/ISOPX(0.5+2.0)/HOT(+5K)/CLOUD(0.6)/IC×2/BG+50%。

**方向一致性**：LNNOX/HVOC 全 5 路线同向 ✓；ISOPX 4/5——**mcm 反向 (+10.1%)**：
MCM 城市配置 NOx 排放高 7 倍 → NOx 饱和区，减氮+增挥解除滴定 → O₃ 上升。
**同一政策情景在不同 NOx 体制下方向相反**——体制判别（FNR/屋脊）是政策分析
的前置步骤，这正是 ridge5/OSSE 实验链的实践价值。

**情景杠杆排序**（fac_v3）：CLOUD −34% > LNNOX −31% > BG+50% +13% > ISOPX −19%
> HOT +4% > HVOC +4% > IC×2 +0.1%。两点政策启示：① IC×2 五天后 +0.1%——初值
5 天内完全洗脱（确定性模拟的 IC 敏感性集中在 day 1）；② HVOC (醛+异戊二烯
排放×3) 仅 +4%——CH₄ 主导的 VOC 池下，羰基/生物源排放控制的 O₃ 杠杆极弱
(CH₄ 有独立因子未动)。BG（背景 O₃）+13% 是被低估的气候杠杆。

**模拟缺失修复**：fac_v3 情景表无 IC/背景因子 → 实验以 config 级缩放演示
（initialConcentrations.config / bc.config），情景表 schema 扩展待定；
三路径排放缩放一致性验证通过（情景表 19.83 / config 19.83 / batch 网格 19.81）。


## 10j. 机制 V2.2：VOC 降解谱 + 支撑化学补全（2026-10-26）

针对外审第二层（VOC 降解谱）与第三层（支撑化学）缺口：**+14 物种/+22 反应 →
74sp/131rx**（元素审计 131/131 全绿）。VOC 降解谱按 CB6/RADM2 集总思路补全：
OLE（烯烃，含 O₃ 臭氧化通道）、ALK（C2⁺ 烷烃）、ARO1/ARO2（甲苯/二甲苯，开环
产物回灌既有 MGLY/CH₃CO₃→PAN 链）、KET（酮）、GLY（乙二醛）；支撑化学补
PNA（HOx/NOx 温度计）、N₂O₅ 白天光解、OH+HNO₃ 暗再循环、ISOPBO₂+CH₃O₂ 交叉。
ALDX 取消（CH₃CHO 升格为 ALD2 集总目标，排放映射覆盖）；萜烯维持单点省略。

**OH 反应性覆盖校准**（核心验收）：扩链前三大类全缺（0%）→ 校准后
**烯烃+芳香族+烷烃合计 ~26%**（燃烧/溶剂/交通燃料端实质入池）；
残余与文献带的偏差（OVOC 34% 偏高、CO 12% 偏低）为 5 天箱体二次 HCHO
累积的结构性特征（稀释 24h 已缓解、无法完全消除——如实文档化）。
**工程修复**：情景分组单源化（grpvoc 双实现收敛 + SO₂ 顺带归组）、
NO₂→PSS 破坏崩溃、dilution 24h 通风（二次累积上限机制）。

**影响**：O₃max 基线 28.7 → 132 ppb（@fnox 2.28 屋脊）——中国东部城市夏季
污染episode 量级（此前 28.7 偏 clean）；哨兵 4 金值重标定（132.10/93.68/3.539/
6.487）；L_NOx/Q 带宽扩至 [0,1.5]（饱和区合法 >1）。六实验重跑中（本节更新时
batch/FT 已重跑，ridge5/scen5/OSSE/sens_adj 按手册命令可复现，数值结论待复核）。


## 10k. 15sp 层整体退役（2026-10-26）

fac_run + kpp_run + ozone.fac + fac2kpp + KPP 工具链成对退役（无归档，知情决策）。
**路线 5→3**（fac_v3/fac_batch_omp/mcm_run）；KPP 外部依赖移除（依赖缩至
gfortran/Python/CVODE）；FT 瘦身 22→10 项并新增 IPR 分解完备性检查
（比值 1.021 ✓，补偿退役的物理层逐过程检验）。**代价（知情）**：异求解器互证
（CVODE vs Rosenbrock）永久消失——67sp 的正确性主张现为文献速率+元素审计+
预算闭合+体制合理性四层，化学产率独立验证（MCM 显式对照）为最深缺口。
连带清理：diagnostics.py 六个 15sp 死函数、species_constraint.py 孤儿、
cfg_audit.c 孤儿。**实验脚本全部改腿并实跑验证**（ridge5 2 腿/scen5 3 腿/
ipr_irr5 3 路线/sens_adj 3 节）。

## 10b. 数值攻坚史（14 坑实录）

1. 元组解包名错位 → 产物全零
2. 手扫 token 行错位 → 两列纯数值表控根除
3. 物种索引魔数三连坑 → 全部常量化
4. max(cosχ,0) kink 死锁 → sin² 日循环（C¹ 光滑）
5. 滴定死锁螺旋 → 背景松弛
6. 松容差重启放飞（MDA8 604 ppb）→ 同档容差+小步长
7. KPP 算子分裂回写（OH→1e14）→ 定位纯机制
8. EBI L 定义错（Inf 爆）→ L=loss_rate/c
9. EBI 日出滴定冲击 → 12 子步+分母下限
10. Fortran D 记法（float() 崩）→ D→E 转换
11. **max(y,0) RHS 折点（39727→4965 步，15×）** → 原始 y（KPP/CMAQ 惯例）
12. **无预条件 SPGMR 黎明崩塌** → 单箱 DENSE；批量才上预条件 Krylov
13. **C-ABI 四连坑** → cvode_c_mod.f90 注释可查
14. **未初始化 ios 静默回落**（一律跑成情景 1）→ 显式置 0

---

## 11. 已知限制与勘误

### 11.0 自纠自查轮记录（2026-10-26，详见 docs/OPERATIONS_MANUAL 附录 A）

- **第一轮·口径同步**（11 项）：退役口径统一至三路线、死常量/死函数清理、孤儿文件删除。
- **第二轮·外部 4 表项核实**：NOy 权重表 13→17 物种（PNA/RONO₂/OLENIT/TERPNIT 补齐，
  氮守恒 FAIL=0，OPE 分母连带修正）；ATOL 调优表单源化至 chem_phys_mod（双胞胎一改双生效）；
  dg.load 列名索引与 IPR 导出链经查排除。
- **第三轮·汇完备性**（新工具：物种死端检查器）：9 个"只耗不产"=排放初级物种（正常）；
  7 个"只产不耗"中 6 个由 drydep 沉降正常闭合；**抓到真 bug OLEOOH 零汇**（化学零损耗
  + drydep 漏配 → 5 天单调涨 66%）→ 补 `OLEOOH 0.5 on`，IPR 沉降账 −3.7e10 验证生效。
- **第五轮·数据链路与 mtime 陷阱**：确认四层链路（templates 格式合同 → **input_v3/ 真数据源**
  → runtime_v3.dat 编译产物 → Fortran 唯一读）；补 3 个可跑的观测约束样例
  （`bc_O3_demo` / `bc_HCHO_demo` / `bc_NO2_sat_demo`）；修 **mtime 同秒陷阱**（改 config 又改回
  → make 误判 up-to-date → runtime.dat 保留中间值 → 无声出错）：Makefile 的
  `runtime_v3.dat` 加 `FORCE` 每次强制重建，并新增 **FT-CFG** 测试逐字节比对防线
  （功能测试 25→26）。
- **第四轮·mcm_run CLI**：`[days] [out_pref]` 双死参修复（接线+mkdir -p），
  ridge5 改名壳删除；速率名零重复定义核实（KRO2NO×5/KRO2HO2×9 为刻意集总）；
  连带修复状态行 print 格式失配（5 描述符对 6 项 → netf 恒空白，现正常输出 netf=11）。
  ⚠ 方法论：加物种必须同步三张表（NOy 权重/ATOL 调优/drydep 沉降）——建议汇完备性
  检查固化进 element_audit 作第三道关卡。

### 11.1 本轮源码实读勘误（相对旧版 README）

| # | 旧说法 | 源码/实跑实况 |
|---|---|---|
| 1 | 冒烟"三实现各跑基线"4/4 | **实为 3 项：FAC/KPP/V3，EBI 不在冒烟**（smoke_test.sh 两处 `smoke` 调用 + V3 段；实跑 3/3） |
| 2 | — | **`input_v3/mechanism_3_full_chain.fac` 是过时副本**（52sp 旧版，缺 O1D 修复）；真编译源=`mechanism/mechanism_3_full_chain.fac`（Makefile gen_v3）。建议删除副本或加同步检查 |
| 3 | 功能测试计数口径多变（16/16 → 22/22 → 10/10 → 该轮 24/24） | 该轮实跑 24/24（含 R10-R13 七条防线断言 + sat2surf 5 项）；**现版 62/62（12 族 FT-\*），见 §7.2** |
| 4 | conc/irr"5 天全程" | **只写末 24 h**（`step>=nsteps-24`）；全程 O₃max 走 stdout |
| 5 | `ipr_from_fac` 读 FAC ipr 文件 | **返回 None 占位**（六本账实现在 fac_v3 resid 内嵌；教学版无独立 ipr 文件） |
| 6 | FNR 阈值方向 | 标准 Sillman/Tonnesen-Dennis 惯例：**FNR<1 → VOC-limited；>2 → NOx-limited**（源码 diagnostics.py:37 行内注释与版本 I 图曾写反，以本表为准） |

### 11.2 机制/物理限制

1. KPP 无物理层（算子分裂回写破坏 Rosenbrock 步长历史，OH→1e14 实测）
2. 指示因子阈值为 Sillman 原文经验值（需按机制/区域重标定）
3. EBI 12 子步兜底（自适应已 1.4×）
4. 纯箱模型（无垂直/输送）
5. HCHO 代理 VOC
6. v3 机制省略：Criegee 显式物种（ISOP/OLE/TERP+O₃ 用产率法，与 CB6R3/MOZART-4 同级）、N₂O₅ 非均相水解（0-D 无气溶胶面）、DMS→SO₂ 分支；RO₂ 交叉缔合仅 5 对代表（实测占 RO₂ 命运 0.34%）
7. ~~序列 BC 未接线~~ **已修复（F3 接线）**：`chem_phys_mod` 物理账③ 序列优先/常量兜底，fac_v3 与 batch 共享
8. fac_batch_omp：网格上限 32×32；只扫 f_NOx×f_VOC（fco/fch4/jfac/dts 固定）；4 核 2.05×
9. 光解参数：18 个在用槽走 MCM 官方 l/m/n（V2.5 槽位对齐）；仅卤素与 N₂O₅/ClO 共 7 个自定义槽为文献正午值（JMAP），photolysis.config 可替换
10. sp_name 兜底：无 SPECIES 行时退化 SPi，情景组缩放失效（重跑 cfgbuild 恢复）
11. 情景组→物种映射硬编码于 `scen_build`（新排放物种需改码）
12. fac_run CSV 的 HNO₃ 列取 `Y(15)` 硬编码（机制物种序变更需同步）

## 12. 参考链接

- **求解器**: [AtChem2](https://github.com/AtChem/AtChem2) · [KPP](https://github.com/KineticPreProcessor/KPP) · [CMAQ](https://github.com/USEPA/CMAQ)
- **指示因子**: [Sillman 1995](https://doi.org/10.1029/94JD03125) · [Kleinman 2005](https://doi.org/10.1029/2005JD005902)
- **MCM v3.3.1**: http://mcm.leeds.ac.uk/MCM/
- **IPR/IRR**: https://cmaq.readthedocs.io
- **课程内**: Day 20/24/26/28/34/35/36
- **操作手册（架构与速查）**: [docs/OPERATIONS_MANUAL.md](docs/OPERATIONS_MANUAL.md)
- **实验操作指南（手把手级）**: [docs/EXPERIMENTS_GUIDE.md](docs/EXPERIMENTS_GUIDE.md)
- **实验运行脚本（类 run.cctm.csh，只跑模拟+CSV）**: [experiments/run/README.md](experiments/run/README.md) —— 九个实验 A–I，`run_all.py` 一键全跑
- **可输入物种全表（111 物种 × 5 通道）**: [templates/SPECIES_INPUT_MATRIX.md](templates/SPECIES_INPUT_MATRIX.md)
- **运行环境与依赖清单 + 混合语言分工**: [ENVIRONMENT.md](environment/ENVIRONMENT.md)（Python/conda、编译器、CVODE、依赖分层、从零复现、`environment/check_env.sh` 体检；**§6b 每语言实现哪些模块/文件/职责**）
- **机制图集与七机制对比**: `figures/README.md`（**本地文件，不入库**；图由 `figures/gen_all.py` 与 `gen_cartoon.py` 重生成）
