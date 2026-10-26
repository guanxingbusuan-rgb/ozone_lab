# modules/ —— Fortran 共享模块层（**核心，不可删**）

> **一句话**：这里是 `fac_v3` 与 `fac_batch_omp` 的**共享物理与数据层**，
> 两个二进制都编译它。**不是预留目录，是主链的一部分。**

---

## 1. 四个文件的定位（2 核心 + 1 工具 + 1 预留）

| 文件 | 类型 | 是否进入二进制 | 作用 |
|---|---|---|---|
| `runtime_dat_mod.f90` | **核心** | ✅ 两个都进 | `runtime_v3.dat` 唯一消费端 |
| `chem_phys_mod.f90` | **核心** | ✅ 两个都进 | 物理四账唯一实现 |
| `nc_output.py` | 工具 | ❌（`make nc` 调用） | CSV → netCDF4 打包 |
| `rad_photol.f90` | **预留钩子** | ❌（独立 demo） | B2 在线光解原型，未接主链 |

Makefile 里的编译依赖（可验证）：

```
facsimile_impl/fac_v3:        ... modules/runtime_dat_mod.f90 modules/chem_phys_mod.f90
facsimile_impl/fac_batch_omp: ... modules/runtime_dat_mod.f90 modules/chem_phys_mod.f90
```

---

## 2. 核心模块 ①：`runtime_dat_mod.f90`（15.8 KB）

**职责**：把 `runtime_v3.dat` 的 12 段内容读进内存，并提供求值函数。
**纯数据模块，不依赖 CVODE**，所以任何驱动都能复用。

### 2.1 读入的 12 段

| 段 | 内容 |
|---|---|
| `RUN` | `ns nj nser nscen` |
| `SPECIES` | 物种名表（111 个） |
| `MODEL` | 模拟天数 / 输出步长 / RTOL / ATOL |
| `IC` | 初值向量 |
| `MET` | 7 个气象量（常数 `C` 或序列 `S`） |
| `BC` | `TAU_MIX` + 111 个常数/序列索引 |
| `PHOT` | `MODE lat doy` + `nj` 行 `J idx l m n` |
| `EMIS` | 111 行 `E 基线 日型 振幅` |
| `DRY` | 111 行 `D vd on/off` |
| `WET` | 111 行 `W a b` |
| `DIL` | 稀释率标量 |
| `SCEN` | 情景名 + 6 个因子 |
| `SERIES` | 序列本体（name n / t 行 / v 行） |

### 2.2 对外提供的函数

| 函数 | 作用 | 谁用 |
|---|---|---|
| `eval_met(t, TEMP, PRES, RH, BLH, RAIN, JFAC, dBLH)` | 气象插值（BLH 含导数） | `chem_phys_mod` |
| `eval_env(TEMP, PRES, RH, M, O2, H2O)` | 理想气体 → 数密度 | `chem_phys_mod` |
| `eval_J(t, JARR)` | SZA → 光解速率 `J = jf·l·cosᵐχ·exp(−n/cosχ)` | `chem_phys_mod` / `fac_v3` |
| `ser_interp(k, t)` | 序列线性插值（端点钳制） | BC / 气象 |
| `box_shape(t, amp)` | 排放日型（06–18 时 sin²） | `chem_phys_mod` |
| `build_groups()` | 物种 → 情景组映射（NOx/VOC/CO/CH₄） | `scen_build` |
| `scen_build(isc)` | 应用情景缩放 | 两个驱动 |

### 2.3 为什么它在 modules/ 而不是 facsimile_impl/

因为它被**两条路线共享**。放 `facsimile_impl/` 会让人误以为只有 `fac_v3` 用。
`mcm_city` 故意**不用**它（那棵子树刻意冻结、自带 `mcm_mod`）——这是设计上的隔离。

---

## 3. 核心模块 ②：`chem_phys_mod.f90`（6.0 KB）

**职责**：**物理过程的唯一实现**（两条路线共享，改一处双生效）。

### 3.1 物理四账（IPR 的 E/D/B/X）

| 账 | 过程 | 公式 |
|---|---|---|
| ① 排放 | `phys_tend` ipr(:,1) | `E0 × emis_scale × 日型` |
| ② 干+湿沉降 | ipr(:,2) | `−(vd/BLH + a·RAIN^b)·y` |
| ③ BC 松弛 + 稀释 | ipr(:,3) | `(bc − y)/TAU_MIX − DIL·y` |
| ④ BLH 膨胀 | ipr(:,4) | `−y·dB/BLH`（仅 dB>0） |

> 机制 P/L 两账不在这里（由调用方从 `oz_mech` 计量组装），所以这里是"物理五项"
> 而非六账 —— 六账 = 机制 P/L + 物理四项。

### 3.2 三个辅助子程序

| 子程序 | 作用 | 为什么单源 |
|---|---|---|
| `tune_atol_by_name(a)` | **ATOL 按物种分档**（62 个物种） | 双胞胎必须同容差，否则结果不可比 |
| `eval_met_dt(t, dt, ...)` | 气象 + 情景升温 | `fac_batch_omp` 的 ΔT 轴 |
| `eval_J_jf(t, jf_abs, J)` | 光解 + 绝对 jfac | `fac_batch_omp` 的 jfac 轴（R11-P0 防二次相乘） |

### 3.3 ATOL 表单源的意义（真实教训）

R9 之前 `tune_atol_by_name` 在 `fac_v3` 和 `fac_batch_omp` 里**各存一份**，
一次改动只改一边 → 双胞胎漂移。现在只有这一份，两边 `use` 同一个。
**这是"为什么它必须在 modules/"的直接理由。**

---

## 4. 工具：`nc_output.py`（2.3 KB，`make nc` 调用）

**作用**：把 `<prefix>_conc.csv` + `<prefix>_irr.csv` 打包成 netCDF4 文件。

```bash
make nc NC_SRC=output/run_A/BAU NC_OUT=output/BAU.nc
```

| 项 | 说明 |
|---|---|
| 输入要求 | 必须**成对**存在（`_conc.csv` + `_irr.csv`） |
| 不支持 | `batch` 的 `x_indicators.csv`（那是网格汇总，非时间序列） |
| 输出结构 | `conc(species, time)`，物种名在**属性** `species_names` 里 |

> 坑：用 `scipy.io.netcdf_file` 读会失败——文件是 **NETCDF4（HDF5）**，
> 需用 `netCDF4.Dataset` 或 `xarray.open_dataset(engine="netcdf4")`。

---

## 5. 预留钩子：`rad_photol.f90`（3.6 KB，**未接主链**）

**状态**：独立可执行 demo，**不参与任何构建**（`make all` 不含，`Makefile` 无编译规则）。

| 项 | 说明 |
|---|---|
| 定位 | B2「在线光解」原型的实现（Beer-Lambert + 太阳天顶角 + 云光学厚度） |
| CLI | `rad_photol lat doy cloud_tau` → 输出 7 个 J 值逐时 |
| 与主链关系 | **无**。主链用 `photolysis-rates_MCM` 的 `l/m/n` 参数化（`eval_J`） |
| 保留理由 | `meteo.config` 的 `RADIAT` 字段是**为它预留的**（见 [UNITS.md §8](../UNITS.md)）；将来若要接真实辐射传输，这是起点 |

⚠ **`RADIAT` 当前无消费者**：`meteo.config` 接受并读入，但没有代码用它。
`cfgbuild` 在你填非 0 值时打印警告。这是已知的、有意保留的空位。

> 如果确定不做在线光解，可删 `rad_photol.f90` 并把 `RADIAT` 从 12 类配置里移除
> ——但目前保留成本极低（3.6 KB、不编译、不运行），且它是 `RADIAT` 语义的唯一说明。

---

## 6. 什么情况下会改这里的文件

| 场景 | 改哪个 | 风险 |
|---|---|---|
| 加新物理过程（如非均相） | `chem_phys_mod.f90` | 高——两条路线同时受影响，必须重跑全部哨兵 |
| 改输入段格式 | `runtime_dat_mod.f90` + `src/cfgbuild.py`（**必须成对改**） | 高——段不匹配会 FATAL |
| 新物种参与情景缩放 | `runtime_dat_mod.f90` 的 `build_groups()` | 中——`FT-SYNC` 守卫会抓漏配 |
| 新物种 ATOL 分档 | `chem_phys_mod.f90` 的 `tune_atol_by_name()` | 低——只影响步数，不影响解 |
| CSV 转 netCDF 格式 | `nc_output.py` | 低 |
| 接在线光解 | `rad_photol.f90` → 需新建编译规则 + 接线 | 高 |

**改完必须**：`make clean-build && make v3 && make test`（两条路线都要重建，
`FT-SYNC`/`FT-FREEZE` 守卫会检查同步性）。

---

## 7. 相关文档

| 文档 | 内容 |
|---|---|
| [ops manual §0.4](../docs/OPERATIONS_MANUAL.md) | 目录职责与编译边界 |
| [ENVIRONMENT.md §6b](../environment/ENVIRONMENT.md) | 逐文件行数 / 角色 / 混合语言分工 |
| [UNITS.md §8](../UNITS.md) | `RADIAT` 无消费者说明 |
| [EXPERIMENTS_GUIDE §3](../docs/EXPERIMENTS_GUIDE.md) | 12 类配置的逐段解读（对应 `runtime_dat_mod` 的读取逻辑） |
