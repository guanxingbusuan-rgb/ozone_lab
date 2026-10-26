# ozone_lab 运行环境与依赖清单（ENVIRONMENT）

> **本文回答**：跑通 ozone_lab 到底需要装什么？版本是多少？装在哪？
> 哪些是**硬依赖**、哪些只在特定路线/工具下才需要？
>
> **实测时间**：2026-10-26 · 全部版本号与路径均为本机实测值。
> **平台**：macOS 12.7.6（Darwin 21.6.0）· x86_64 · 4 逻辑核

---

## 目录

- [1. 依赖总览（一页表）](#1-依赖总览一页表)
- [2. Python / conda 隔离环境](#2-python--conda-隔离环境)
- [3. 编译器](#3-编译器)
- [4. 第三方数值库（CVODE/SUNDIALS）](#4-第三方数值库cvodesundials)
- [5. KPP 与 bison/flex（历史路线）](#5-kpp-与-bisonflex历史路线)
- [6. 各二进制的实际链接依赖](#6-各二进制的实际链接依赖)
- [6b. 混合语言分工：每种语言实现什么、在哪、为什么](#6b-混合语言分工每种语言实现什么在哪为什么)
- [7. 依赖分层：硬依赖 vs 可选依赖](#7-依赖分层硬依赖-vs-可选依赖)
- [8. 从零复现环境（完整步骤）](#8-从零复现环境完整步骤)
- [9. 环境自检脚本](#9-环境自检脚本)
- [10. 常见环境问题](#10-常见环境问题)

---

## 1. 依赖总览（一页表）

| # | 依赖 | 实测版本 | 路径 | 层级 | 缺失后果 |
|---|---|---|---|---|---|
| 1 | **Python** | 3.12.9 | `~/miniconda3/envs/venv` | **硬** | 配置编译/诊断/实验全不可用 |
| 2 | numpy | 1.26.4 | 同上 | **硬** | cfgbuild/诊断失败 |
| 3 | matplotlib | 3.10.0 | 同上 | 半硬 | 图与实验脚本失败；编译/运行不受影响 |
| 4 | scipy | 1.11.4 | 同上 | 可选 | 部分实验脚本（插值/统计） |
| 5 | netCDF4 | 1.7.2 | 同上 | 可选 | `nc_output.py`（CSV→nc） |
| 6 | Pillow (PIL) | 11.1.0 | 同上 | 可选 | 图片后处理/验证 |
| 7 | fontTools | 4.55.3 | 同上 | 可选 | 字体覆盖自检（排版排查用） |
| 8 | **gfortran** | 13.2.0 | `~/gcc-13.2.0/bin` | **硬** | 无法编译求解核 |
| 9 | **gcc** | 13.2.0 | `~/gcc-13.2.0/bin` | 半硬 | 仅 C 审计器/`mechc_c.py` 后端需要 |
| 10 | **SUNDIALS CVODE** | 2.7.0（库 2.9.0/2.7.0 版本号见 §4） | `~/soft/atmenv/cvode` | **硬** | 所有 Fortran 路线无法链接 |
| 11 | libgomp（OpenMP） | 随 gcc 13.2.0 | `~/gcc-13.2.0/lib` | 半硬 | `fac_batch_omp` 无并行（仍可串行） |
| 12 | KPP | 3.2.0 | `~/soft/src/KPP-3.2.0` | **已退役** | 无影响（KPP 路线已删） |
| 13 | bison / flex | 系统 2.3 / 2.6.4 | `/usr/bin` | **已退役** | 无影响（仅 KPP 生成时需要） |

**核心结论**：跑通现役三路线只需 **3 样**——Python 环境（numpy+matplotlib）+
gfortran 13.2.0 + CVODE 2.7.0。KPP/bison/flex **已不在构建链中**。

---

## 2. Python / conda 隔离环境

### 2.1 当前环境实况

```
解释器      : $HOME/miniconda3/envs/venv/bin/python
Python 版本 : 3.12.9
conda 环境名: venv
conda 基座  : $HOME/miniconda3
```

**注意**：该环境名为 `venv` 但由 miniconda 创建（不是 `python -m venv`）。
`conda` 命令本身未装进该 env（只在 base），故 `import conda` 会失败——
这不影响任何功能，项目脚本只用纯 Python 包。

### 2.2 项目实际使用的 Python 包（按用途）

| 包 | 用在哪 | 必需性 |
|---|---|---|
| **numpy** | `src/cfgbuild.py`（校验/换算）、`diagnostics/*`、全部实验脚本 | 硬 |
| **matplotlib** | `figures/*.py`（9 张机制图）、`experiments/fig_*.py`、`make figures` | 半硬 |
| scipy | 部分实验（插值/统计/优化） | 可选 |
| netCDF4 | `modules/nc_output.py`（CSV→netCDF 转换） | 可选 |
| Pillow | 图片检查/裁剪（本次 QA 用） | 可选 |
| fontTools | 字体覆盖自检（图集排版排查） | 可选 |
| csv/json/re/argparse | 全部脚本 | Python 标准库 |

**没有任何 `conda install` 专属包**——所有依赖都是 PyPI 通用包，理论上
`python -m venv` + `pip install numpy matplotlib` 也能跑。

### 2.3 版本要求（最低可行）

| 包 | 最低版本 | 原因 |
|---|---|---|
| Python | **3.9+** | 用了 `match`? 否；但 `dict \|` 合并、f-string 需 3.6+；实测 3.12 |
| numpy | **1.20+** | `np.trapz`/`np.nanmean`/`np.random.default_rng`（1.17+） |
| matplotlib | **3.5+** | 本项目实测 3.10.0；**3.10 起 PDF 后端支持逐字形字体回退**（图集排版依赖此特性，旧版 PDF 中文会变豆腐块） |

> ⚠ **matplotlib 版本不是小事**：`figures/gen_all.py` 与 `gen_cartoon.py` 依赖
> "显式字体列表 + 逐字形回退"，该行为在 3.10 实测可用。若降到 3.4 以下，
> PDF 导出的中文会整串变方块（PNG 仍正常）。

### 2.4 环境变量依赖

```bash
# Makefile 已固化，无需手动设（但知道它做什么有助排错）
export PATH := ~/miniconda3/envs/venv/bin:$(PATH)   # 保证 make 里的 python 是 venv 的
export KPP_HOME := ~/soft/src/KPP-3.2.0             # 历史遗留，已不在构建链
```

---

## 3. 编译器

### 3.1 实测版本

```
GNU Fortran (GCC) 13.2.0   →  ~/gcc-13.2.0/bin/gfortran
gcc (GCC) 13.2.0           →  ~/gcc-13.2.0/bin/gcc
```

**为什么不用系统编译器**：macOS 自带 `clang` 无 Fortran；Apple 的 gfortran
（若有）版本通常过老，且**不提供 `libgomp`/`libquadmath` 的完整布局**。
本项目所有 Fortran 代码为 free-form Fortran 90/2003，实测 13.2.0 通过。

### 3.2 编译选项（Makefile 固化）

```make
FFLAGS := -Wall -O2 -ffree-line-length-none
```

| 选项 | 作用 | 为什么必须 |
|---|---|---|
| `-Wall` | 全警告 | 抓未初始化变量（曾抓出 `ios` 静默回落 bug） |
| `-O2` | 优化 | 运行时间从 ~2 s 降到 ~0.2 s |
| **`-ffree-line-length-none`** | 取消 132 列限制 | 机制生成码（`oz_mech.f90`）单行极长，不加必编译失败 |
| `-fopenmp`（仅 batch） | OpenMP | `fac_batch_omp` 的线程并行 |

### 3.3 编译器版本的下限

| 特性 | 最低 gfortran | 本机 |
|---|---|---|
| `iso_c_binding` + `c_ptr`/`c_f_pointer`（CVODE C API 直连） | 4.6+ | ✓ 13.2 |
| `execute_command_line`（`mkdir -p`） | 4.9+ | ✓ 13.2 |
| Fortran 2008 `block`/`newunit=` | 4.9+ | ✓ 13.2 |
| OpenMP 4.5 | 6+ | ✓ 13.2 |

**建议 9.0+**；实测最稳的是 13.2.0（KPP 时代的 5 补丁链就是针对它调的）。

---

## 4. 第三方数值库（CVODE/SUNDIALS）

### 4.1 版本与路径

```
安装根: ~/soft/atmenv/cvode
头文件: ~/soft/atmenv/cvode/include/{cvode,nvector,sundials}
库    : ~/soft/atmenv/cvode/lib/
```

| 库文件 | 版本号 | 谁在用 |
|---|---|---|
| `libsundials_cvode.2.dylib` | **2.9.0** | 全部 Fortran 路线（核心积分器） |
| `libsundials_nvecserial.2.dylib` | 2.7.0 | 全部（串行向量） |
| `libsundials_fnvecserial.2.dylib` | 2.7.0 | fac_v3 / mcm_run（**fcmix Fortran 接口**） |
| `libsundials_fcvode.a` | 2.7.0 | 静态库（fcmix 符号） |

> **版本号说明**：SUNDIALS 2.7.0 发行版中，`cvode` 模块的 dylib 兼容版本标记为
> 2.9.0，而 `nvecserial` 为 2.7.0——**同一发行版内部模块版本不同是正常的**，
> 不是混装。判断依据看安装目录整体版本（2.7.0）。

### 4.2 两个接口并存（关键）

| 接口 | 符号 | 用在哪 | 特点 |
|---|---|---|---|
| **fcmix（Fortran 封装）** | `FCVMALLOC`/`FCVODE`/`FCVDENSE` | `fac_v3`、`mcm_run` | 简单，但 **COMMON 块限制单进程只能有 1 个 CVODE 问题** |
| **C API（ISO_C_BINDING 直连）** | `CVodeCreate`/`CVode`/`N_VNew_Serial` | `fac_batch_omp` | 每线程独立实例，可 OpenMP 并行 |

**这是为什么 batch 必须用 C API**：fcmix 的全局 COMMON 块让多线程共享同一个
积分器状态，并行必然出错（详见 `facsimile_impl/cvode_c_mod.f90` 注释）。

### 4.3 链接方式

```bash
-L$(CVODELIB) -Wl,-rpath,$(CVODELIB) \
  -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial
```
- `-Wl,-rpath` 把库路径写进二进制 → **运行时不需要设 `DYLD_LIBRARY_PATH`**
- 实测三个二进制均无 `not found` 依赖

### 4.4 若需自行编译 CVODE

```bash
# SUNDIALS 2.7.0 源码
./configure --prefix=$HOME/soft/atmenv/cvode \
            --enable-fcmix --disable-shared --enable-static
# 注意: --enable-fcmix 才能生成 libsundials_fcvode（Fortran 路线必需）
make -j4 && make install
```

**常见坑**：不启用 `--enable-fcmix` 时，`fac_v3` 链接会报
`undefined symbol: fcvmalloc_`——这是 Fortran 接口没编进来。

---

## 5. KPP 与 bison/flex（历史路线）

### 5.1 当前状态：**已退役，不影响任何功能**

| 项 | 状态 |
|---|---|
| `kpp_impl/` 目录 | **已删除**（2026-10 轮） |
| Makefile 中的 KPP 目标 | **已移除**（`grep kpp Makefile` 无输出） |
| `src/fac2kpp.py` | 脚本保留（历史工具，不被 make 调用） |
| `~/soft/src/KPP-3.2.0` | 磁盘上仍在，但**构建不用** |
| `bison`/`flex` | 仅 KPP 代码生成时需要；现未使用 |

**若未来要恢复 KPP 路线**，需要：

| 依赖 | 版本 | 说明 |
|---|---|---|
| KPP | 3.2.0 | `~/soft/src/KPP-3.2.0` |
| **bison** | **3.0+** | ⚠ **系统 `/usr/bin/bison` 只有 2.3，太老**，需 `brew install bison` 并在 PATH 前置 |
| flex | 2.6.4+ | 系统版本够用 |
| gfortran | 13.2.0 | KPP 生成码有 5 处需补丁适配，见 `src/fac2kpp.py` 注释 |

---

## 6. 各二进制的实际链接依赖

实测 `otool -L`（macOS）：

### `facsimile_impl/fac_v3`（工程主线）
```
libsundials_cvode.2.dylib       (~/soft/atmenv/cvode/lib)   ← 2.9.0
libsundials_fnvecserial.2.dylib (~/soft/atmenv/cvode/lib)   ← 2.7.0
libsundials_nvecserial.2.dylib  (~/soft/atmenv/cvode/lib)   ← 2.7.0
libgfortran.5.dylib             (~/gcc-13.2.0/lib)
libquadmath.0.dylib             (~/gcc-13.2.0/lib)
```

### `facsimile_impl/fac_batch_omp`（批量，含 OpenMP）
```
libsundials_cvode.2.dylib       ← 2.9.0
libsundials_nvecserial.2.dylib  ← 2.7.0
libgfortran.5.dylib
libgomp.1.dylib                 ← OpenMP 运行时
libquadmath.0.dylib
# 注意: 无 fnvecserial —— 它用 C API, 不需要 fcmix 封装层
```

### `mcm_city/mcm_run`（MCM 对照）
```
libsundials_cvode.2.dylib
libsundials_fnvecserial.2.dylib
libsundials_nvecserial.2.dylib
libgfortran.5.dylib
libquadmath.0.dylib
```

**共性**：三个二进制都只依赖 **SUNDIALS + gfortran 运行时**（+ batch 的 libgomp），
全部路径已通过 `rpath` 固化，无 `DYLD_LIBRARY_PATH` 需求。

---

## 6b. 混合语言分工：每种语言实现什么、在哪、为什么

> 本项目是**四语言协作**：Python（前端/诊断/可视化）+ Fortran（数值求解核）+
> C（第三方求解库 + 审计证人 + 代码生成后端）+ Shell（构建编排）。
> 下面按语言给出**文件级职责表**——每个文件在系统里扮演什么角色。

### 6b.1 语言分工总图

```
┌─────────────────────────────────────────────────────────────────┐
│ Python (23 个 .py)  ——  前端 + 后端分析 + 可视化                 │
│   src/          配置编译 / 机制编译 / C 后端 / 卫星转换          │
│   diagnostics/  指示因子与过程诊断（纯函数，只吃 CSV）           │
│   experiments/  8 个实验脚本                                     │
│   tests/        功能测试（24 断言）                              │
│   figures/      9+9 张机制图生成                                 │
│   modules/      nc 输出（唯一的 Python-in-modules）              │
└────────────────────────┬────────────────────────────────────────┘
                         │ 生成 runtime_v3.dat（纯文本平面文件）
                         │ 生成 gen_v3/oz_mech.f90（Fortran 源码）
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│ Fortran (9 个 .f90)  ——  数值求解与物理层                        │
│   modules/runtime_dat_mod.f90   输入解析 + 气象/光解求值         │
│   modules/chem_phys_mod.f90     物理四账（双胞胎共享）           │
│   modules/rad_photol.f90        在线辐射光解（独立 demo）        │
│   facsimile_impl/oz_fac_v3.f90  工程主线驱动                     │
│   facsimile_impl/fac_batch_omp.f90  批量并行驱动                 │
│   facsimile_impl/cvode_c_mod.f90    CVODE C API 的 Fortran 绑定  │
│   mcm_city/mcm_driver.f90       MCM 对照案例                     │
│   */gen*/oz_mech.f90            mechc.py 自动生成（勿手改）      │
└────────────────────────┬────────────────────────────────────────┘
                         │ 链接
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│ C  ——  三类角色                                                  │
│   libsundials_*  (第三方)       ODE 求解引擎本体（BDF/Newton）   │
│   audit/oz_audit.c            独立"证人"：NOx 守恒 + O₃ 收支      │
│   (src/mechc_c.py 按需生成)    机制 → C 速率函数（第二执行后端）  │
└─────────────────────────────────────────────────────────────────┘
                         ▼
┌─────────────────────────────────────────────────────────────────┐
│ Shell / Makefile  ——  构建编排与冒烟                              │
│   Makefile          all/v3/mcm/batch/nc/figures/test/clean…      │
│   check_env.sh      环境体检（14 项只读检查）                     │
│   tests/smoke_test.sh  fac_v3 冒烟                                │
└─────────────────────────────────────────────────────────────────┘
```

### 6b.2 Python 层（23 个文件）—— 前端 + 分析 + 可视化

#### A. 输入前端（`src/`，3 个）

| 文件 | 行数 | 角色 | 做什么 |
|---|---|---|---|
| **`src/cfgbuild.py`** | 395 | **配置编译器（唯一前端入口）** | 读 12 类人类可读配置 → 校验（物种/单位/单调/点数/槽位）→ 写 `runtime_v3.dat`（Fortran 唯一输入）。**运行时零文本解析的秘密全在这里** |
| `src/mechc_c.py` | 112 | 机制→C 代码生成 | `.fac` → C 速率函数 → `gcc -shared` → ctypes 直调。绕过 f2py 构建链 |
| **`src/sat2surf.py`** | 282 | 卫星数据转换器 | TROPOMI 柱浓度 → 地表等效 ppb 序列（直观法主准则 + QA 过滤 + BLH 插值 + 点数抽稀） |

#### B. 机制编译器（`facsimile_impl/mechc.py`，1 个）

| 文件 | 行数 | 角色 | 做什么 |
|---|---|---|---|
| **`facsimile_impl/mechc.py`** | 145 | **FACSIMILE→Fortran 编译器** | 解析 `.fac` 的物种表/速率名/反应行 → 生成 `gen_v3/oz_mech.f90`（含 `ns/nr/nq` 常量、`update_rates` 速率函数）+ `oz.reac`/`oz.prod` 计量表。**改机制 → 它重跑 → gfortran 重编**，这是"机制即代码"的实现 |

#### C. 诊断与分析（`diagnostics/`，3 个）

| 文件 | 行数 | 角色 | 做什么 |
|---|---|---|---|
| **`diagnostics/indicators.py`** | 196 | **路线无关指示因子引擎** | 12 指标（FNR/L_NOx·Q⁻¹/OPE…）；从 `.fac` 的 FORMULA 段自动判定 NOy；两套预算表（Q×6、产率）在初始化期解析为索引 |
| `diagnostics/diagnostics.py` | 70 | CSV 读入与 IPR 解析 | `load()`（读 conc/irr）、`ipr_from_fac()`（六账梯形积分） |
| `diagnostics/csvfloat.py` | 50 | 数值解析共享件 | 处理 gfortran 的 `8.27280-108` 三位指数省略 E 的坑（R13） |

#### D. 实验脚本（`experiments/`，8 个）

| 文件 | 行数 | 角色 |
|---|---|---|
| `experiments/obs_osse.py` | 212 | F：观测约束 OSSE（真值→采样→同化→技能） |
| `experiments/scen5.py` | 187 | I：情景 × 路线响应矩阵 |
| `experiments/sens_adj.py` | 176 | H：OAT 灵敏度 + 截断 SVD |
| `experiments/ridge5.py` | 143 | C：屋脊跨路线对比 |
| `experiments/ipr_irr5.py` | 140 | D：多路线 IPR/IRR 对照 |
| `experiments/hourly_pa.py` | 129 | G：逐小时过程分析 |
| `ipr_closure_check.py` | — | D-1：六账闭合校验 |
| `experiments/fig_ridge5.py` | 54 | C 的绘图 |

#### E. 测试与可视化（`tests/` + `figures/` + `scan_suite.py`）

| 文件 | 行数 | 角色 |
|---|---|---|
| **`tests/functional_test.py`** | 369 | **24 条功能断言**（含 FT-SAT 5 条 + 选择性输入 2 条 + 金值哨兵） |
| `figures/gen_all.py` | 915 | 9 张 SCI 版机制图（A–I） |
| `figures/gen_cartoon.py` | 557 | 9 张漫画场景版（comic_A–I） |
| `figures/fig_indicators_ridge.py` / `fig_response_v3.py` | 71 / 58 | 响应面与指示因子工程图 |
| `scan_suite.py` | 61 | 响应面套件运行器 |
| `modules/nc_output.py` | 66 | CSV → netCDF（**唯一的 Python 版 modules/ 文件**） |

### 6b.3 Fortran 层（9 个 .f90）—— 数值求解与物理

| 文件 | 行数 | 类型 | 角色 |
|---|---|---|---|
| **`modules/runtime_dat_mod.f90`** | 364 | module | **输入消费端**：解析 `runtime_v3.dat`（SPECIES/IC/MET/BC/PHOT/EMIS/DRY/WET/DIL/SCEN/SERIES 12 段）+ `eval_met`（气象插值）+ `eval_env`（M/O₂/H₂O）+ `eval_J`（SZA→光解）+ `build_groups`（情景分组）+ `scen_build`。纯数据模块，无求解器依赖，双胞胎共享 |
| **`modules/chem_phys_mod.f90`** | 107 | module | **物理四账唯一实现**（双胞胎共享）：排放 / 干湿沉降 / BC 松弛+稀释 / BLH 膨胀；`tune_atol_by_name`（ATOL 调优单源） |
| `modules/rad_photol.f90` | 82 | module+demo | **在线辐射光解**（Beer-Lambert + 太阳天顶角），独立可执行 demo，未接入主链 |
| **`facsimile_impl/oz_fac_v3.f90`** | 281 | program | **工程主线驱动**：读 runtime → CVODE fcmix BDF+DENSE → 111 物种全列 CSV + 198 IRR + 六账 IPR |
| **`facsimile_impl/fac_batch_omp.f90`** | 622 | program | **批量并行驱动**：网格扫描 + 预算表（NOy 17 物种/Q/L_NOx/OPE）+ 24 列指示因子 + OpenMP |
| `facsimile_impl/cvode_c_mod.f90` | 108 | module | **CVODE C API 的 Fortran 绑定**：`CVodeCreate`/`CVode`/`N_VNew` 等 `bind(C)`，绕开 fcmix 的 COMMON 块限制（多线程必需） |
| `mcm_city/mcm_driver.f90` | 263 | program | **MCM 对照案例**：29sp 子集 + 城市排放 + SPGMR |
| `facsimile_impl/gen_v3/oz_mech.f90` | 247 | **生成码** | mechc.py 产出（勿手改）：`ns=74, nr=131, nq=105, nj=41` + `update_rates` + 计量表 |
| `mcm_city/gen/oz_mech.f90` | 221 | **生成码** | 同上（MCM 29sp/71rx 版） |

**Fortran 层的设计原则**：
- **`modules/` 是共享层**（runtime_dat_mod + chem_phys_mod），fac_v3 与 batch 编译时都包含
- **驱动层各自独立**（facv3_mod / batch_mod / mcm_mod 互不 `use`），只通过 CSV 契约对外
- **生成码不进版本管理逻辑**：改 `.fac` → mechc.py 重生成，绝不手改

### 6b.4 C 层（1 个源文件 + 2 类外部角色）

| 项 | 角色 | 说明 |
|---|---|---|
| **`libsundials_cvode` 等（第三方）** | **ODE 求解引擎本体** | BDF/Newton/SPGMR 的实现是 C；Fortran 通过 fcmix（`FCVODE`）或 C API（`CVode`）调用。**整个项目的数值能力建立在它之上** |
| **`audit/oz_audit.c`**（40 行） | **独立"证人"审计器（C 语言旁证）** | `nox_budget`（NOx 家族望远镜闭合，Kahan 求和）+ `o3_budget`（ΔO₃=∫(P−L)dt）。**用与 Fortran 完全无关的实现算同一笔账**——若两边不一致，说明有一边错了。需手动 `gcc -shared` 编译为 dylib 后由 Python 调 |
| `src/mechc_c.py`（Python，非 C 文件） | **机制→C 代码生成后端** | 按需把 `.fac` 编译成 C 速率函数（`gcc -shared` → ctypes），提供 Fortran 之外的第二条机制执行路径 |

> **C 层的设计意图**：不是"用 C 重写求解器"，而是提供**独立的第二实现**做交叉验证
> + 借助第三方 C 库获得成熟数值算法。当前唯一的 C 源文件 `oz_audit.c` 是历史
> "C 证人"角色，未接入自动构建（需手动编译）。

### 6b.5 Shell / Makefile（2 个 .sh + 1 个 Makefile）

| 文件 | 角色 | 做什么 |
|---|---|---|
| **`Makefile`**（13 个目标） | **构建编排与依赖管理的唯一入口** | `all`（两二进制）/ `v3` / `mcm` / `batch` / `nc` / `figures` / `smoke` / `functional` / `test` / `clean` / `clean-build` / `clean-output` / `FORCE`。**按时间戳判断依赖**：改 config 不触发编译，改 `.fac`/`.f90` 才重编；`runtime_v3.dat` 挂 `FORCE` 强制重建（防同秒编辑漏判） |
| **`tests/smoke_test.sh`** | 冒烟 | fac_v3 跑基线，查输出非空 + ≥24 行 + O₃ 列健全 |
| **`check_env.sh`** | 环境体检 | 14 项只读检查（Python/编译器/CVODE/二进制/依赖解析） |

### 6b.6 数据/文本资产（非代码，但属于接口合同）

| 类型 | 文件 | 角色 |
|---|---|---|
| `.fac`（1 个在役） | `mechanism/mechanism_3_full_chain.fac` | **机制唯一事实源**（111sp/198rx，V2.4） |
| `.fac`（1 个存档） | `mechanism/teaching/mechanism_1_ozone_tutorial.fac` | 15sp 教学档（无路线承载，存档参考） |
| `.config`（9 个） | `input_v3/*.config` | 12 类配置中的生效部分（人类可读） |
| `.template`（6 个） | `templates/*.template` | **格式合同**（程序零直读） |
| `.txt`（序列） | `input_v3/constraints/*.txt` | 观测/气象时间序列（`met_*_demo.txt` 6 条气象样例 + `bc_*_demo.txt` 3 条约束样例） |
| `.dat`、`.csv` | `runtime_v3.dat`、`output/*` | 编译产物 / 实验产物（`output/` 已 gitignore） |
| `.md`（6 个） | README / EXPERIMENTS_GUIDE / OPERATIONS_MANUAL / **ENVIRONMENT** / UNITS / SPECIES_INPUT_MATRIX | 文档体系 |

### 6b.6b 模型维度澄清：0-D 箱 + 参数扫描 ≠ 空间网格

| 概念 | 是什么 | 不是什么 |
|---|---|---|
| **ozone_lab 本体** | **0-D 箱体**：状态向量 = `ns` 个物种（74 个 ODE），无 x/y/z | 不是 3D/4D 空间模型 |
| **`fac_v3`** | 跑**一个** 0-D 箱（一个情景） | — |
| **`fac_batch_omp`** | 跑**一批**独立 0-D 箱（fnox×fvoc×jfac×ΔT 参数组合，最多 262,144 个） | **不是** "一个 4 维空间场" |
| **"四维网格"** | 四维**参数空间**（排放/光解/温度缩放因子） | **不是** 空间坐标网格 |
| 空间输运 | **完全没有**（`grep advect/diffus/transport` 在 modules/ 与驱动层均 0 命中） | 需要空间分布请用 CMAQ / WRF-Chem / GEOS-Chem |

**为什么容易误会**：`fac_batch_omp` 的输出是二维表（fnox × fvoc），画出来像
"等值线图"，视觉上与空间分布图相似；但那两个轴是**缩放因子**不是经纬度。
本质是"**跑 81 个盒子，把结果排成表**"，不是"模拟一个 81 格的空间域"。

与 AtChem2 的关系：**两者都是 0-D 箱**。差别在于 ozone_lab 多了一层
**批量参数扫描驱动器**（`fac_batch_omp`）——这层 AtChem2 没有（它一次一个盒子，
批量靠外部脚本循环）。

---

### 6b.7 语言边界与契约（为什么这样分工）

| 边界 | 契约 | 为什么用这种语言 |
|---|---|---|
| 12 类配置 → `runtime_v3.dat` | 纯文本平面文件（Fortran 顺序表控读） | Python 擅长文本/校验/单位换算；Fortran 不该做正则和容错 |
| `.fac` → `oz_mech.f90` | 生成 Fortran 源码（含速率函数 + 计量表） | 机制是"代码"：把化学式编译成可内联的 Fortran，比解释执行快 1–2 个量级 |
| Fortran → CSV | `<pref>_conc.csv` + `_irr.csv` + `_ipr.csv` | 无二进制依赖，Python 可直接分析；也便于人工检查 |
| CSV → 指标/图 | `diagnostics/` + `figures/` | 数据分析与绘图是 Python 生态的强项 |
| Fortran → CVODE | fcmix（`FCVODE`）或 C API（`CVode`） | 求解器本体是 C；两条路径分别服务"易用"与"并行" |

**一句话总括**：**Python 管"数据怎么说"，Fortran 管"方程怎么算"，C 管"底层怎么解 + 独立怎么验"，Shell 管"什么时候谁该重跑"。**

---

## 7. 依赖分层：硬依赖 vs 可选依赖

### 7.1 硬依赖（缺了跑不了）

| 依赖 | 缺失时的表现 |
|---|---|
| Python 3.9+ | `make` 立即失败（cfgbuild 用 python3） |
| **numpy** | `cfgbuild.py` 导入失败 → 无法生成 runtime.dat |
| **gfortran 13.2.0** | `make all` 编译失败 |
| **CVODE 2.7.0（含 fcmix）** | 链接失败：`undefined symbol: fcvmalloc_` |

### 7.2 半硬依赖（缺了功能受限但核心可跑）

| 依赖 | 缺了会怎样 | 变通 |
|---|---|---|
| matplotlib | 图生成失败，但 cfgbuild/fac_v3 正常 | 只看 CSV |
| libgomp | `fac_batch_omp` 编译失败或串行 | 改用 fac_v3 逐点跑 |
| gcc | C 审计器（`liboz_audit.dylib`）不可编 | 跳过 C 证人，用 Python 审计 |
| 系统 `mkdir`/`sh` | fac_v3 的 `ensure_dir` 与 smoke 脚本失败 | 手动建目录 |

### 7.3 纯可选（特定工具/实验才用）

| 依赖 | 用于 |
|---|---|
| scipy | `experiments/` 部分脚本的插值/优化 |
| netCDF4 | `modules/nc_output.py`（CSV→nc） |
| Pillow / fontTools | 图片检查、字体覆盖自检 |
| KPP / bison / flex | **已退役**，仅未来恢复 KPP 路线时用 |
| FRUIT / openlibm / numdiff | 早期测试框架，现不在链 |

---

## 8. 从零复现环境（完整步骤）

### 8.1 最小可行环境（跑通现役三路线）

```bash
# ── ① Python 环境 ──
conda create -n venv python=3.12 -y
conda activate venv
pip install numpy matplotlib            # 硬依赖
pip install scipy netCDF4 pillow fonttools   # 可选

# ── ② 编译器 ──
# 方式 A: Homebrew
brew install gcc                        # 提供 gfortran-13
# 方式 B: 官方 binary（本项目实测方式）
#   https://gcc.gnu.org/wiki/GFortranBinaries 下载, 解包到 ~/gcc-13.2.0

# ── ③ CVODE 2.7.0 (含 Fortran 接口) ──
mkdir -p ~/soft/src && cd ~/soft/src
curl -LO https://computing.llnl.gov/projects/sundials/download/sundials-2.7.0.tar.gz
tar xzf sundials-2.7.0.tar.gz && cd sundials-2.7.0
./configure --prefix=$HOME/soft/atmenv/cvode \
            --enable-fcmix --disable-shared --enable-static
make -j4 && make install

# ── ④ 改 Makefile 路径（若与默认不同）──
cd /path/to/ozone_lab
# 编辑 Makefile 前 7 行的 PY / FC / CC / CVODELIB

# ── ⑤ 构建 + 自检 ──
make all
make test        # 冒烟 1/1 → 功能 35/35 → 审计 198/198
```

### 8.2 完整环境（含图集与全部实验）

```bash
pip install scipy netCDF4 pillow fonttools
python figures/gen_all.py       # 9 张 SCI 版机制图
python figures/gen_cartoon.py   # 9 张漫画版
python experiments/ridge5.py    # 屋脊跨路线对比
```

---

## 9. 环境自检脚本

把下面保存为 `check_env.sh` 一键体检（**只读，不改环境**）：

```bash
#!/bin/sh
# ozone_lab 环境自检 —— 只读检查, 不修改任何东西
PY=${PY:-python3}
ok=0; bad=0
chk() {  # chk <描述> <命令>
  if eval "$2" >/dev/null 2>&1; then echo "  [OK]   $1"; ok=$((ok+1));
  else echo "  [MISS] $1"; bad=$((bad+1)); fi
}

echo "== Python =="
chk "python 可用 ($($PY --version 2>&1 | cut -d' ' -f2))" "$PY --version"
chk "numpy"      "$PY -c 'import numpy'"
chk "matplotlib" "$PY -c 'import matplotlib'"
chk "scipy (可选)"      "$PY -c 'import scipy'"
chk "netCDF4 (可选)"    "$PY -c 'import netCDF4'"

echo "== 编译器 =="
GF=${FC:-gfortran}
chk "gfortran ($($GF --version 2>/dev/null | head -1))" "$GF --version"
chk "gcc" "gcc --version"

echo "== CVODE =="
CV=${CVODE_ROOT:-$HOME/soft/atmenv/cvode}
chk "libsundials_cvode"     "test -e $CV/lib/libsundials_cvode.dylib -o -e $CV/lib/libsundials_cvode.so"
chk "libsundials_fcvode (fcmix Fortran 接口)" \
    "test -e $CV/lib/libsundials_fcvode.a -o -e $CV/lib/libsundials_fcvode.so"
chk "cvode.h 头文件"        "test -e $CV/include/cvode/cvode.h"

echo "== 项目构建 =="
chk "fac_v3 二进制存在"       "test -x facsimile_impl/fac_v3"
chk "fac_batch_omp 二进制存在" "test -x facsimile_impl/fac_batch_omp"
chk "runtime_v3.dat 存在"     "test -f runtime_v3.dat"

echo "== 二进制依赖可解析 =="
if command -v otool >/dev/null; then
  n=$(otool -L facsimile_impl/fac_v3 2>/dev/null | grep -c "not found")
  [ "$n" = "0" ] && { echo "  [OK]   fac_v3 无未解析依赖"; ok=$((ok+1)); } \
                 || { echo "  [MISS] fac_v3 有 $n 个未解析依赖"; bad=$((bad+1)); }
fi

echo
echo "== 汇总: $ok OK, $bad MISS =="
[ "$bad" = "0" ] && echo "环境就绪 ✓" || echo "请先补齐 MISS 项"
```

运行：

```bash
sh check_env.sh
```

---

## 10. 常见环境问题

| 症状 | 根因 | 处置 |
|---|---|---|
| `undefined symbol: fcvmalloc_` | CVODE 编译时未开 `--enable-fcmix` | 重装：`./configure --enable-fcmix …` |
| `dyld: Library not loaded: libsundials_cvode.2.dylib` | rpath 未写入（手工编译时漏 `-Wl,-rpath`） | 重编或设 `DYLD_LIBRARY_PATH=$HOME/soft/atmenv/cvode/lib` |
| `gfortran: command not found` | PATH 无 `~/gcc-13.2.0/bin` | 用 Makefile 里的绝对路径，或 `export PATH=$HOME/gcc-13.2.0/bin:$PATH` |
| `Free form line length exceeds 132` | 缺 `-ffree-line-length-none` | Makefile 已有；手工编译时别漏 |
| matplotlib 图里中文全是方块 | matplotlib < 3.10 且 PDF 后端无逐字形回退 | 升级到 3.10+，或只输出 PNG |
| `ModuleNotFoundError: cfgbuild`（sat2surf） | 从非 `src/` 目录直接 import | 已修（工具内自动加 `sys.path`）；或 `cd src` 运行 |
| `import conda` 失败 | conda 只在 base 环境 | 无害，项目不需要 conda Python API |
| `bison 2.3` 报版本太老 | 系统 bison 太旧（仅 KPP 需要） | `brew install bison` 并 PATH 前置；**现役路线不用** |
| `make` 用错 Python | PATH 里 python 不是 venv 的 | Makefile 已用绝对路径 `PY := …/venv/bin/python` |
| OpenMP 加速比测不出 | 系统负载 > 核数（实测 load 94 时批处理会掉案例） | 看 `uptime`，负载正常时再测量 |

---

## 附录：本机环境快照（2026-10-26 实测）

```
平台      macOS 12.7.6 (Darwin 21.6.0) x86_64 · 4 逻辑核
Python    3.12.9 @ ~/miniconda3/envs/venv
numpy     1.26.4
matplotlib 3.10.0        ← PDF 中文回退依赖此版本
scipy     1.11.4
netCDF4   1.7.2
Pillow    11.1.0
fontTools 4.55.3
gfortran  13.2.0 @ ~/gcc-13.2.0/bin
gcc       13.2.0 @ ~/gcc-13.2.0/bin
CVODE     2.7.0 @ ~/soft/atmenv/cvode
          （libsundials_cvode 2.9.0 / nvecserial 2.7.0 / fnvecserial 2.7.0）
KPP       3.2.0 @ ~/soft/src/KPP-3.2.0   [已退役, 不参与构建]
bison     2.3 (/usr/bin)                 [仅 KPP 需要; 若恢复需 3.0+]
flex      2.6.4 (/usr/bin)               [同上]
```

### 相关文档

| 文档 | 内容 |
|---|---|
| [README.md](README.md) | 项目总览、路线清单 |
| [EXPERIMENTS_GUIDE.md](EXPERIMENTS_GUIDE.md) | 实验操作（环境章节见 §1） |
| [OPERATIONS_MANUAL.md](OPERATIONS_MANUAL.md) | 架构、速查、编译边界 |
| [UNITS.md](UNITS.md) | 单位与换算参考 |
| [templates/SPECIES_INPUT_MATRIX.md](templates/SPECIES_INPUT_MATRIX.md) | 111 物种 × 5 通道可输入矩阵 |
