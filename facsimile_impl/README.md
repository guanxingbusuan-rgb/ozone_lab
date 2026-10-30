# facsimile_impl/ —— 仿真引擎层（本项目的"发动机舱"）

> **一句话**：化学积分**真正发生的地方**。两个可执行体（`fac_v3` 单箱深度、
> `fac_batch_omp` 批量广度）+ 一个机制编译器（`mechc.py`）+ 生成码（`gen_v3/`）。
> **本目录全部产物都不入库**（二进制与生成码均可重建）。

---

## 1. 文件清单与角色

| 文件 | 类型 | 入库? | 角色 |
|---|---|---|---|
| `mechc.py` | Python | ✅ | **机制编译器**：`.fac` → `gen_v3/oz_mech.f90` + 计量表 |
| `oz_fac_v3.f90` | Fortran | ✅ | **fac_v3 源码**：单箱全功能驱动 |
| `fac_batch_omp.f90` | Fortran | ✅ | **fac_batch_omp 源码**：批量响应面 + OpenMP |
| `cvode_c_mod.f90` | Fortran | ✅ | CVODE **C API 的 ISO_C_BINDING 包装**（batch 专用） |
| `fac_v3` | 二进制 | ❌ | 可执行体（`make all` 重建） |
| `fac_batch_omp` | 二进制 | ❌ | 可执行体（`make all` 重建） |
| `gen_v3/` | 生成码 | ❌ | `mechc.py` 的产物（`make all` 自动再生） |

## 2. `mechc.py` —— 机制编译器（`.fac` 的唯一消费者）

```
mechanism/mechanism_3_full_chain.fac
    ↓  python mechc.py ../mechanism/... gen_v3
gen_v3/oz_mech.f90   ← update_rates(TEMP,M,O2,H2O,JARR,Y,q,p)：127 个命名速率 + 203 条反应速率
gen_v3/oz.reac       ← 每反应的反应物索引表（供 IPR/IRR 归因）
gen_v3/oz.prod       ← 每反应的产物索引+计量系数表
```

支持的 FACSIMILE 子集：命名速率互引用、`J<n>` 光解索引、`TEMP/M/O2/H2O` 环境量、
`@` 幂运算、Troe 三参数、分数计量、行内 `!` 注释（自动剥离）。
**防呆**（三道，全部编译期响错误）：
1. 重复速率定义 → 报错（历史上 Br+O₃ 重复行静默覆盖过）
2. 行内 `!` 注释 → 解析前剥离（V2.5：曾致反应行尾注释并入产物 token，
   `oz.prod` 静默丢产物 → ODE 与 IPR 双错）
3. 计量表生成时**未知物种 token 硬失败**（V2.5）：白名单 = 机制 VARIABLE 物种 ∪
   储库 {O₂/H₂O/M/N₂/CO₂/H₂}；物种名手误（如 ISOPB0）从此编译期拦截，
   不再静默丢项。

## 3. `fac_v3` —— 单箱深度路线

```bash
cd facsimile_impl && ./fac_v3 <情景名|序号> <输出前缀> <runtime.dat>
# 例: ./fac_v3 SCEN_BAU ../CASE-analysis/output/case ../CASE-analysis/output/runtime_case.dat
```

日常不手敲——`CASE-analysis/run.py` 每轮自动 `cfgbuild` 出 runtime 再 spawn 本二进制。

| 特性 | 说明 |
|---|---|
| 积分 | CVODE BDF + DENSE（fcmix 绑定），内步推进 + 过冲插值回整点 |
| 求解失败韧性 | `FCVREINIT` + 缩 `MAX_STEP` 重试一次 |
| 输出三件套 | `_conc.csv`（h + 111 物种快照）/ `_irr.csv`（h + 203 反应**本窗积分**）/ `_ipr.csv`（h + 111×6 账**本窗积分**）；**保留全部 post-spinup 步**，不是末 24h。首列 `h` 是**输出步序号**（非小时），换算小时 = `h×dtout/3600` |
| 情景 | 读 `runtime.dat` 的 SCEN 段，按名或序号选 |
| 物理层 | `modules/chem_phys_mod`（与 batch 共享，改一处双生效） |

## 4. `fac_batch_omp` —— 批量广度路线

```bash
cd facsimile_impl && ./fac_batch_omp <ng≤64> <out.csv> <runtime.dat> \
    [jfac列表] [dt列表] [fmin] [fmax] [模式] [第9参] [LIN]
# EKMA-FNR 默认（1000%/100% ⇒ 11 点、因子 1–11、IC 模式、排放关闭）:
#   ./fac_batch_omp 11 raw.csv runtime_fnr.dat 1.0 0.0 1 11 IC NOX:HCHO LIN
# EKMA-emission 默认（625%/25% ⇒ 26 点、因子 1–7.25、排放扰动、第二轴=VOC 组）:
#   ./fac_batch_omp 26 raw.csv runtime_emission.dat 1.0 0.0 1 7.25 EMIS 1 LIN
```

| 参数位 | 含义 |
|---|---|
| 1–3 | `ng`（单边格点数，数组上限 `fs(64)`）、输出 CSV 前缀、runtime |
| 4–5 | jfac 列表、ΔT 列表（缺省 `1.0` / `0.0`；**省略时补占位值，勿错位**） |
| 6–7 | `fmin` / `fmax`（缺省 `0.3` 与 `0.3×1.5^(ng-1)`，即等比网格） |
| 8 | 模式：`EMIS`（缺省，排放乘子扫描）/ `IC`（初值扫描，**排放与 BC 张弛同时关闭**） |
| 9 | `IC` 模式 = `ic1spec:ic2spec`（组名 `NOX`/`VOC`/`CO`/`CH4` 或物种名）；`EMIS` 模式 = 第二轴组（1=VOC 2=ISOP 3=CO 4=CH4） |
| 10 | `LIN` = 线性网格（EKMA 的 % 扰动轴用这个），否则等比 |

| 特性 | 说明 |
|---|---|
| 网格 | ng×ng 两轴 × jfac × ΔT（各自上限 64 / 16 / 16） |
| 并行 | OpenMP `collapse(2)` + `schedule(dynamic)`；**案例间并行，单箱内仍单线程** |
| C API | `cvode_c_mod.f90`：fcmix 用 COMMON 块存指针 → 单进程只能 1 个 CVODE 问题；C API 每线程独立 mem 才能并行 |
| 输出 | 24 列指示因子（含 `_pm` 午后口径）+ `_ipr.csv` + `_irr.csv` + `_state.csv`（每案例一行末 24 步均值；**不落全时间序列**）。统计窗 = **末 24 个输出步**（`O3max` 等不是 8h 日最大）。判据只用**官方 5 族 + NOy**（见 [`../analysis_lib/README.md`](../analysis_lib/README.md)），`O3NOy`/`HCHONOy`/`O3NOz` 为备用扩展。结束时打印 `失败 N / 降级 K`：**降级 = 已有统计但窗口不完整**（调用方应记录） |
| 侧车坐标 | `<outf>_irr.csv` / `_ipr.csv` 的前 4 列 = 该案例的 fnox/fvoc/jfac/dt，行序与主 CSV 一致；`<outf>_state.csv` 给同一窗口内的**一致均值对**（`O3max/O3mean` 与 `HCHOmax/NO2mean` 一个取极值一个取均值，不能直接当状态空间两轴用） |
| spin-up | 读取 `spinup_cfg`；要求 post-spinup ≥24 步，否则 FATAL；聚合窗口不混入 spin-up |
| 在线守卫 | 两条都致命：① `nfail>0`（有案例零统计求解失败）；② **氮守恒不变量**——每个案例的**首个 RHS 快照**检一次 `|Σ w_N·dy/dt| ≤ 1e-4·max(L1,1)`，违例计入模块共享 `nNviol`（`!$omp critical` 内累加；`ndbg_once` 是 `threadprivate`，每案例开头重置 ⇒ **逐案例确定性质检**，不是"每线程一次"的随机抽查）。收尾打印 `失败 N / 降级 K / 氮守恒违例 M`，**任一非 0 → `[FATAL]` + `stop 1`**，两个 EKMA 的 `run.py` 据此拒出表并把计数写进 `run_meta.txt` |
| 降级（非致命） | `nwarn`：已有 >0 统计但窗口不完整的案例，行仍写出、不 `stop`；调用方必须在 `run_meta.txt` 里看到它 |
| 账本单位 | `_ipr.csv` / `_irr.csv` 侧车 = 末 24 步的**平均倾向** [molec/cm³/s]（`ipracc/cnt`）；而 `fac_v3` 的同名文件 = **窗口积分** [molec/cm³]。两者相差 dtout 倍，跨驱动对比必须先换算 |
| 预算表 | `batch_budget_init` 启动期解析 Q/L_NOx/产率索引，失败硬停 |

**与 fac_v3 的关系**：共享同一生成码、同一物理层、同一 runtime 通路——
改一次机制/物理，两条路线自动同步，且互为回归对照。

## 5. `gen_v3/` 里的文件

| 文件 | 内容 |
|---|---|
| `oz_mech.f90` | 生成的速率代码（勿手改） |
| `oz.reac` / `oz.prod` | 反应物/产物计量表（IPR/IRR 归因的索引源） |
| `*.mod` | gfortran 模块文件（编译中间产物） |

## 6. 重建与耗时

```bash
make all           # mechc.py 生成 gen_v3 + gfortran 两个二进制，~30 s
make clean-build   # 删二进制/.mod/生成码（不动各程序 output/）
```

改 `.fac` / `modules/*.f90` / 本目录 `.f90` 任一，`make` 会按依赖自动重编；
**只改各程序 `config_override/` 数据不需要重编**（runtime 是纯数据，
`run.py` 每轮自动跑 cfgbuild，0.4 s 重出）。

## 7. 相关目录

| 目录 | 关系 |
|---|---|
| `modules/` | 共享物理/数据层（两二进制都编译它） |
| `mechanism/audit/` | 机制静态审计（element_audit / mcm_compare / jslot_audit） |
| `inputproc/cfgbuild.py` | 上游：12 类配置 → `runtime_*.dat` |
| `CASE-analysis/`、`EKMA-FNR/`、`EKMA-emission/` | 下游：三个自完备程序 |
