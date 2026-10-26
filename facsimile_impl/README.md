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
| `fac_v3` | 二进制 | ❌ | 可执行体（`make v3` 重建） |
| `fac_batch_omp` | 二进制 | ❌ | 可执行体（`make v3` 重建） |
| `gen_v3/` | 生成码 | ❌ | `mechc.py` 的产物（`make v3` 自动再生） |

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
   不再静默丢项。配套 `FT-MECHK` 三断言（注释剥离/硬失败/全量计量表核对）。

## 3. `fac_v3` —— 单箱深度路线

```bash
./fac_v3 <情景名|序号> <输出前缀> <runtime.dat>
# 例: ./fac_v3 SCEN_BAU ../output/run1 ../runtime_v3.dat
```

| 特性 | 说明 |
|---|---|
| 积分 | CVODE BDF + DENSE（fcmix 绑定），内步推进 + 过冲插值回整点 |
| 求解失败韧性 | `FCVREINIT` + 缩 `MAX_STEP` 重试一次 |
| 输出三件套 | `_conc.csv`（h+111 物种，末 24h）/ `_irr.csv`（h+203 反应窗口积分）/ `_ipr.csv`（h+111×6 账） |
| 情景 | 读 `runtime.dat` 的 SCEN 段，按名或序号选 |
| 物理层 | `modules/chem_phys_mod`（与 batch 共享，改一处双生效） |

## 4. `fac_batch_omp` —— 批量广度路线

```bash
./fac_batch_omp <ng≤32> <out.csv> <runtime.dat> [jfac列表] [dt列表]
# 例: ./fac_batch_omp 9 ../output/batch_v3_indicators.csv ../runtime_v3.dat 0.6,1.0 0,5
```

| 特性 | 说明 |
|---|---|
| 网格 | fnox×fvoc 各 ng 点（0.3×1.5^i）× jfac × ΔT，最多 262,144 箱 |
| 并行 | OpenMP `collapse(2)` + `schedule(dynamic)`；**案例间并行，单箱内仍单线程** |
| C API | `cvode_c_mod.f90`：fcmix 用 COMMON 块存指针 → 单进程只能 1 个 CVODE 问题；C API 每线程独立 mem 才能并行 |
| 输出 | 24 列指示因子（含 `_pm` 午后口径）+ `_ipr.csv` + `_irr.csv`（每案例一行窗口均值） |
| 在线守卫 | 氮守恒不变量（每线程首快照自检）+ `nfail>0` → FATAL 退出 |
| 预算表 | `batch_budget_init` 启动期解析 Q/L_NOx/产率索引，失败硬停 |

**与 fac_v3 的关系**：共享同一生成码、同一物理层、同一 runtime 通路——
改一次机制/物理两条路线自动同步，互为回归（FT 契约断言）。

## 5. `gen_v3/` 里的文件

| 文件 | 内容 |
|---|---|
| `oz_mech.f90` | 生成的速率代码（勿手改） |
| `oz.reac` / `oz.prod` | 反应物/产物计量表（IPR/IRR 归因的索引源） |
| `*.mod` | gfortran 模块文件（编译中间产物） |

## 6. 重建与耗时

```bash
make v3        # mechc.py + gfortran 两个二进制，~30 s
make mcm       # mcm_city 对照路线（不在此目录，独立子树）
make clean-build   # 清二进制+生成码，保留 output/
```

改 `.fac` / `modules/*.f90` / 本目录 `.f90` 任一，`make` 会按依赖自动重编；
**只改 `input_v3/` 数据不需要重编**（`runtime_v3.dat` 是纯数据，cfgbuild 0.4 s 重出）。

## 7. 相关目录

| 目录 | 关系 |
|---|---|
| `modules/` | 共享物理/数据层（两二进制都编译它） |
| `mcm_city/` | 独立对照子树（MCM 官方 29sp，零代码耦合） |
| `src/cfgbuild.py` | 上游：12 类配置 → `runtime_v3.dat` |
| `tests/` | 下游：两个二进制的契约与金值哨兵 |
