# mcm_city —— MCM v3.3.1 CH₄ 子集对照路线（独立子树，非核心程序）

> **一句话**：这是**第三条路线**，用真实 MCM 官方机制跑 29 物种城市光化学，
> 作为 `fac_v3` / `fac_batch_omp` 的**跨机制参照基准**。
> 它与两条核心工程程序**零代码耦合**，不在 `make all` 里，需显式 `make mcm`。

---

## 1. 与核心两条程序的关系（重要）

### 1.1 不相关的地方（无耦合）

| 项 | `fac_v3` / `fac_batch_omp` | **`mcm_run`** |
|---|---|---|
| 机制 | `mechanism/mechanism_3_full_chain.fac`（111sp/198rx） | `mcm_city/input/mechanism_mcm_ch4.fac`（29sp/71rx） |
| 生成码 | `facsimile_impl/gen_v3/` | **`mcm_city/gen/`（自己的）** |
| 输入 | `runtime_v3.dat`（由 `input_v3/` 编译） | **不读任何 runtime 文件**，全部硬编码在 `mcm_driver.f90` |
| 共享模块 | `modules/runtime_dat_mod.f90` + `chem_phys_mod.f90` | **完全不用**（自带 `mcm_mod`） |
| 构建 | `make all` | **`make mcm`（独立目标）** |
| 气象/排放/沉降 | `input_v3/*.config` 可配 | **源码里写死** |

代码级验证（零引用）：

```bash
grep -rn "mcm_city\|mcm_run" facsimile_impl/*.f90 modules/*.f90 src/*.py
# → 无输出
```

唯一共享的是**第三方库 CVODE**（SUNDIALS），以及同一个 `facsimile_impl/mechc.py`
编译器（把 `.fac` 编成 `oz_mech.f90`）。

### 1.2 相关的地方（被实验层引用）

虽然核心程序不依赖它，但**实验层用它做跨机制对比**：

| 实验 | 用法 | 脚本 |
|---|---|---|
| **实验 C**（屋脊跨路线对比） | `run_mcm(5, pref, fnox=..., fvoc=...)` 沿 fnox 扫描，与 batch 74sp/111sp 比屋脊位置 | `experiments/run/run_C_ridge.py` |
| **实验 I**（情景 × 路线矩阵） | 同一情景在 mcm_run 上重跑，检查**方向一致性** | `experiments/run/run_I_scenario_matrix.py` |

结论意义（已在实验输出里观察到）：

- `run_I` 中 mcm 29sp 对 **LNNOX / ISOPX** 出现方向相反（+4.70 vs −47.67 等）
- 这**不是 bug**，而是证明"**简化机制会改变结论方向**"的关键证据——
  29sp CH₄ 子集没有异戊二烯化学，HS 情景在它上面只能表达为"排放加倍但无降解通道"。

> 因此本目录的定位是：**可选的对照参照物**，不是运行 `fac_v3` 的前置依赖。
> 跑 A/B/D/E/F/G/H/M 八个实验**不需要** `mcm_run`；只有 C 和 I 需要。

---

## 2. 机制规模与来源

| 项 | 值 |
|---|---|
| 物种 / 反应 | **29 / 71** |
| 命名速率 | 139 |
| 光解槽位 | 12 |
| 来源 | **MCM v3.3.1 官方 CH₄ 子集**（AtChem2 自带同一份） |
| 机制文件 | `input/mechanism_mcm_ch4.fac`（FACSIMILE 语法，与主机制同构） |
| 光解参数 | `input/photolysis-rates_MCM_v331` |

**为什么选 CH₄ 子集而不是全谱**：这是"真实 MCM 官方机制"里**最小可运行且有完整无机块**
的档（NO₃/N₂O₅/HO₂NO₂/HONO 齐全 + Troe KMT01-03/KBPAN/KFPAN + HSO₃/SO₃ 硫链 + Cl）。
全谱 5832 物种的官方导出在 `mechanism/reference/mcm_v3.3.1_export.fac`，
仅作**系数比对基准**，不可运行（无配套光解表与初始场）。

---

## 3. CLI 与硬编码

```bash
cd mcm_city && ./mcm_run <days> <out_pref> [fnox] [fvoc]

# 例：5 天、基线、无缩放
./mcm_run 5 ../output/mcm_base 1.0 1.0
```

| 参数 | 默认 | 说明 |
|---|---|---|
| `days` | 5 | 模拟天数 |
| `out_pref` | `mcm_city` | 输出前缀（可含路径，已支持 `mkdir -p`；buffer 512 支持长绝对路径） |
| `fnox` | 1.0 | NO/NO2 排放乘子（屋脊实验用） |
| `fvoc` | 1.0 | HCHO 排放乘子 |

### 3.1 全部硬编码（这是"刻意冻结"的直接后果）

| 量 | 值 | 位置 |
|---|---|---|
| 排放 NO | `2.0e7` molec/cm³/s（20 ppb/h 城市交通） | `mcm_driver.f90:57` |
| 排放 NO₂ / CO / CH₄ / HCHO | `1.0e6` / `5.0e6` / `5.0e3` / `2.0e5` | 同上 (58-61) |
| 干沉降 vd | O₃ `8.0e-6` / HNO₃ `2.0e-5` / HCHO `3.0e-6` | `mcm_driver.f90:55` |
| 初值 | 全零 | — |
| 背景松弛 | `TAU=86400`，O₃ `6.2e11` / CO `1.5e12` / CH₄ `4.5e13` | `mcm_driver.f90:36-37` |
| 排放日型 | `sin²(π(t−6h)/12h)`，06–18 时 | `dayf()` |
| TH 气象 | TEMP 288.15 / PRESS 1013.25（见 `input/envvars_reference.config`） | `mcm_driver.f90` |
| 光解 | SZA(40°N, DOY 172)，12 槽位 | `eval_J_mcm()` |

⚠ **想用自己的数据 → 请用 `fac_v3`**。本路线刻意不支持配置化，
目的是保持"官方机制 + 固定场景"这个**不变的参照系**。

---

## 4. 输出

| 文件 | 内容 |
|---|---|
| `<pref>_conc.csv` | **15 列精选**：`h`, O3, NO, NO2, OH, HO2, HCHO, CO, H2O2, HNO3, NO3, N2O5, HONO, **P_O3, L_O3**（O₃ 产/耗通量），末 24h |
| `<pref>_irr.csv` | `h` + **71 反应** = 72 列（小时窗积分） |
| `<pref>_ipr.csv` | `h` + **29 物种 × 6 账** = 175 列（P/L/E/D/B/X，小时窗积分） |

---

## 5. 构建

```bash
make mcm          # = 生成 gen/oz_mech.f90 + 编译 mcm_run
```

`make all` **不包含**它（只构建 `fac_v3` + `fac_batch_omp`）。
`make clean-build` / `make clean` 会删掉 `gen/`，需重新 `make mcm`。

求解器：**CVODE BDF + DENSE（29×29）**。
曾试无预条件 SPGMR，在黎明光解过渡段崩塌（`10b#12`），已回退 DENSE。

---

## 6. 已知限制（刻意保留）

| 限制 | 原因 |
|---|---|
| 排放/IC/沉降/气象全部写死 | 保持官方机制的"纯参照"定位 |
| 只有 29 物种 | MCM CH₄ 子集范围；无 ISOP 化学 |
| 无 runtime 通路 | 不与 `input_v3/` 耦合，避免配置漂移影响参照系 |
| 光解仅 12 槽位 | 子集机制的槽位需求 |
| 情景只能表达 fnox/fvoc | 无 CO/CH₄/JFAC/ΔT 轴（对应 `run_I` 里标 `—` 的格点） |

---

## 7. 相关文档

| 文档 | 内容 |
|---|---|
| [mechanism/README.md](../mechanism/README.md) §10 | 七机制对比（本路线 = "真实 MCM 子集"档） |
| [mechanism/V25_ASSESSMENT.md](../mechanism/V25_ASSESSMENT.md) | 机制升级评估（含"为何 mcm 与主机制方向会相反"） |
| [experiments/run/README.md](../experiments/run/README.md) §5 | C/I 实验对 mcm_run 的前置构建要求 |
| [ENVIRONMENT.md](../ENVIRONMENT.md) §6b | 混合语言分工与子目录职责 |
