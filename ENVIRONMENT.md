# ENVIRONMENT —— ozone_lab 运行与编译环境（2026-10-29 程序级精简版）

> 本文件是当前仓库**唯一**的环境/依赖说明。适用对象：三个自完备程序
> `CASE-analysis/`、`EKMA-FNR/`、`EKMA-emission/` 及其共享层。
> 通用实验层（`experiments/ configs/ output/ runtime_v3.dat tests/ mcm_city/ docs/`）已退役。

---

## 1. 目录与角色

| 目录/文件 | 角色 | 何时被读取 |
|---|---|---|
| `CASE-analysis/` | 个例分析（六账 + IRR + RIR），入口 `run.py` / `plot.py` | 用户运行 |
| `EKMA-FNR/` | S2 初始浓度 EKMA（HCHO × NO₂），入口 `run.py` / `plot.py` | 用户运行 |
| `EKMA-emission/` | S1 排放 EKMA（VOC × NOx），入口 `run.py` / `plot.py` | 用户运行 |
| `mechanism/` | 机制 `.fac`（V2.5+ 冻结）+ `audit/` 三个静态审计 | 运行期（.fac）/ 审计 |
| `facsimile_impl/` | 执行器源码与二进制；`gen_v3/`（生成码 + 计量表 + .mod） | 运行期（二进制 / CASE 读计量表） |
| `modules/` | `runtime_dat_mod.f90`、`chem_phys_mod.f90`（物理层唯一实现） | **仅编译期**（已编入二进制） |
| `inputproc/` | `cfgbuild.py`（配置→runtime 编译器）、`sat2surf.py`（卫星柱→表面） | 运行期（cfgbuild） |
| `analysis_lib/` | `csvio`（CSV 读写）、`indicators`（机制感知指标 / O₃ 收支）等 | 运行期（import） |
| `Makefile` | 只建两个二进制；`mcmcheck` / `jcheck` / `clean-build` | 构建期 |
| `check_env.sh` | 13 项环境自检（只读） | 任意 |
| `README.md` | 仓库总览 | — |

**输入/输出边界**：每个程序只读自己的 `config_override/`（完整配置集：12 个顶层配置 + 17 个 constraints），
产物只落自己的 `output/`；**不回退任何全局配置**，程序之间互不引用。

---

## 2. 运行期依赖矩阵

| 依赖 | CASE-analysis | EKMA-FNR | EKMA-emission |
|---|---|---|---|
| `facsimile_impl/fac_v3` | ✅ | — | — |
| `facsimile_impl/fac_batch_omp` | — | ✅ | ✅ |
| `facsimile_impl/gen_v3/oz.reac|oz.prod` | ✅（IRR 归因计量表） | — | — |
| `inputproc/cfgbuild.py` | ✅ | ✅ | ✅ |
| `mechanism/mechanism_3_full_chain.fac` | ✅ | ✅ | ✅ |
| `analysis_lib/csvio.py` | ✅ | ✅ | ✅ |
| `analysis_lib/indicators.py` | ✅（P(O₃)=G−D） | — | — |
| `modules/*.f90` | ❌（已在二进制内） | ❌ | ❌ |
| Python 包 | `numpy`；`plot.py` 另需 `matplotlib` | 同左 | 同左 |

运行时**不需要**：根 `configs/`、`runtime_v3.dat`、`experiments/`、`tests/`、`docs/`。

---

## 3. 工具链（本机实测）

| 项 | 版本 / 路径 |
|---|---|
| OS / 架构 | macOS x86_64，21.6.0 内核 |
| CPU | 4 逻辑核 / 2 物理核（OpenMP 默认线程 = 逻辑核，可用 `OMP_NUM_THREADS` 覆盖） |
| Shell | `/bin/zsh` |
| Python | 3.12.9（`$HOME/miniconda3/envs/venv/bin/python`） |
| numpy / matplotlib / scipy | 1.26.4 / 3.10.0 / 1.11.4 |
| Fortran | GNU Fortran 13.2.0（`$HOME/gcc-13.2.0/bin/gfortran`） |
| C | gcc 13.2.0（同前缀） |
| SUNDIALS | CVODE 2.9.0（`$HOME/soft/atmenv/cvode/lib`） |
| OpenMP | gfortran `-fopenmp`（`collapse(2)` + `schedule(dynamic)`，案例间并行） |

`Makefile` 支持覆盖：`CONDA_ENV` / `GCC_ROOT` / `PY` / `FC` / `CC` / `CVODELIB`。
`check_env.sh` 的 gfortran 回退顺序：`$FC` → PATH 里的 `gfortran` → `$HOME/gcc-13.2.0/bin/gfortran`。

---

## 4. 编译

```bash
make all            # 只在二进制缺失/源码变更时重建（当前已是最新）
make clean-build    # 删二进制/.mod/生成码（不动各程序 output/）
make mcmcheck       # 机制速率 vs MCM 官方导出
make jcheck         # 光解槽位语义
```

- 两个目标：`facsimile_impl/fac_v3`、`facsimile_impl/fac_batch_omp`；
- 编译命令显式包含 `modules/runtime_dat_mod.f90` 与 `modules/chem_phys_mod.f90`，模块**静态并入**二进制；
- `.mod` 接口文件生成在 `facsimile_impl/gen_v3/`（仅编译期需要）；
- 运行时动态依赖只有 SUNDIALS + libgfortran/libquadmath，无仓库内路径。

**何时必须重新编译**：改 `modules/*.f90`、`mechanism/*.fac`、`facsimile_impl/*.f90`，或换机器/新克隆。
只改各程序的 `config_override/` 或 Python 代码时，**不需要编译**。

---

## 5. 环境自检

```bash
bash check_env.sh     # 13 项：Python/包、编译器、二进制与动态链接、目录结构
```

当前实测：**13 OK / 0 MISS —— 环境就绪 ✓**

---

## 6. 机制静态审计（不跑仿真）

| 脚本 | 检查 | 实测 |
|---|---|---|
| `mechanism/audit/element_audit.py` | 逐反应 C/H/N/O/S/Cl/Br 配平 | 203/203 PASS |
| `mechanism/audit/mcm_compare.py` | 48 条速率 vs MCM，最大偏差 | 8.9%（阈值 10%） |
| `mechanism/audit/jslot_audit.py` | 光解槽位语义 / JMAP / README 一致性 | 0 错位 |

`jslot_audit.py` 的 `photolysis.config` 来源顺序：`configs/`（已退役）→
`CASE-analysis/config_override/` → `EKMA-emission/config_override/`。

**运行期在线门（非静态审计）**：`fac_batch_omp` 每案例首个 RHS 快照检一次
`|Σ w_N·dy/dt| ≤ 1e-4·max(L1,1)`，违例计数与求解失败计数一样 → **`stop 1`**，
两个 EKMA 的 `run.py` 据此拒出表；当前三个存档实测均为 0 违例。

---

## 7. 输入数据硬约束（`cfgbuild.py` 编译期强制）

1. 序列格式：`t[s] 值` 两列、时间严格递增、单序列 **≤128 点**；气象 + BC 共享 **40 条**序列池（NSER_MAX，容纳 26 物种同化约束）；
2. 气象序列**禁止** `# unit:` 指令行（浓度换算因子会放大 2.46e10 倍）；
3. 温度用 K；RH 0–100；**BLH 必须 >0（常数与 `file:` 序列同口径校验）**；IC 单位仅 `cm3|ppb|ppm`；
4. `model.parameters`：天数 >0、输出步长 >0、RTOL/ATOL >0、`0 ≤ SPINUP < 天数`；
5. 排放单位：`molec/cm3/s` 或 `ppb/h`（`# unit:` 行声明，缺省 molec/cm3/s）。

输出侧的两条口径（`h` 列 = 输出步序号；`fac_v3` 账本=窗口积分 / `fac_batch_omp` 侧车=平均倾向，
差 dtout 倍）写在根 [README.md](README.md) §6 数据契约里，跨程序比对前先看它。

---

## 8. 参考基线（归档实测，用于判断环境是否正常）

| 程序 | 口径 | 参考值 |
|---|---|---|
| CASE-analysis | 4 天 / spinup 3 / SCEN_BAU / 26 物种同化 / `--per-species` | `nst≈4454`，O₃max ≈ 149.5 ppb，六账 111 物种 × 25 步，RIR 有效窗 **4–18 时**，36 次模拟 / 212.7 s（不含 `--per-species` 的一轮 ≈ 45 s） |
| EKMA-FNR | 1 天 / 排放关闭 / SPINUP 0；**存档 = 默认 11×11 = 121 案例**（0–1000%/100%） | O₃max 34.14–119.21 ppb，失败/降级/氮守恒违例 0-0-0，墙钟 64.7 s |
| EKMA-emission | 4 天 / SPINUP 3；**存档 = 默认 26×26 = 676 案例**（0–625%/25%，线性网格） | O₃max 15.26–281.41 ppb，基线点 77.12，失败/降级/氮守恒违例 0-0-0；耗时随负载波动（899 s / 1458 s 两轮同参数） |

> 各轮的时间/网格/失败·降级·氮守恒计数以 `output/run_meta.txt` 为准。
> 2026-10-30 重跑对照（三程序全量）：CASE 的 `case_{conc,ipr,irr}.csv` 与 `runtime_case.dat`
> **与原存档逐字节相同**（`nst=4454 / O3max=149.49`）；FNR 归档回到 121 例、基线点 34.14；
> emission 676 例数据行逐位不变（`ridge/derived` 连字节都一致）→ 二进制打印串改动、
> `run.py` 捕获逻辑与逐案例守恒检查均为**数值中性**；守恒门另在 `/tmp` 沙箱用故意改错的
> `noy_add` 权重做过负向验证（9/9 案例触发 → `[FATAL]` + `stop 1` + `run.py` rc=1 拒出表）。

差异超过 ~10% 时先查机制/工具链版本与 `config_override/` 是否被改动。

---

## 9. 已退役（一句话）

通用实验层（`experiments/ configs/ output/ runtime_v3.dat tests/ mcm_city/ docs/`）于 2026-10-29 退役，
其历史内容已删除；需要时从 git HEAD 取回。
