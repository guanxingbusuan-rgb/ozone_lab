# analysis_lib/ —— 诊断与指标库（纯 Python，仿真之后的数值消费）

> **一句话**：积分是 `fac_v3` / `fac_batch_omp` 的事；**这里只管"跑完之后"**——
> 读 `_conc.csv` / `_irr.csv` / `_ipr.csv`，算指示因子，兜底解析 Fortran 的怪数值格式。
> 两个文件都是**被 import 的库**，没有独立入口，零第三方依赖（只用 numpy）。

---

## 1. 两个文件

| 文件 | 行数 | 定位 | 谁在用 |
|---|---|---|---|
| `csvio.py` | 115 | **CSV 契约读写的唯一实现**（读/写 + gfortran 指数兜底） | 三个程序的 `run.py` 与 `plot.py`（全部） |
| `indicators.py` | 302 | **机制感知指标**：`Mechanism`（NOy/碳数/计量）+ `p_o3_gd`（O₃ 收支）；`route_*` 为 12 量工具箱 | `Mechanism`/`p_o3_gd` 由 `CASE-analysis/run.py` 在用；**`route_indicators`/`route_timeseries` 当前无在役调用者**，保留作 batch 24 列的手工互证入口 |

依赖方向恒为单向：各程序 `run.py` / `plot.py` → `analysis_lib` →（只读）CSV。
`analysis_lib` 不 import 任何仿真或编排代码，也不写任何目录约定。

---

## 2. `csvio.py` —— CSV 契约

### 为什么需要兜底解析

gfortran 的 `ES13.5` 等格式在 `|指数|≥100` 时**静默丢掉 `E`**，输出形如 `3.07932-174`
的数——Python `float()` 直接抛异常。现行写出格式已改 `ES15.5E3`（恒三位指数），
但解析一律走本模块，不要裸 `float()`。

### API

| 导出 | 作用 |
|---|---|
| `PPB = 2.46e10` | 1 ppb 的数密度（1013.25 hPa / 298.15 K），**全仓唯一定义点** |
| `fix_fortran_exp(tok)` | `3.07932-174` → `3.07932E-174`；正常格式原样返回 |
| `fnum(tok)` | 修复后转 float；空串/不可解析 → `nan`（不抛异常，便于整表解析） |
| `read_table(path)` | → `(header, 二维 float)`；跳过 `write_csv` 写出的 `# ` 注释行 |
| `read_rows(path)` + `fcol(hdr, rows, name)` | 含**字符串列**（体制标签、站点名等）的表用这套；`read_rows` 保原字符串，`fcol` 按列名取数值数组 |
| `write_csv(path, header, rows, note=)` | 写 CSV，`note` 逐行成 `# ` 头；自动建父目录 |

⚠ 数值表用 `read_table`，含标签列的表用 `read_rows + fcol` —— 用错会让文本列变一片 nan。

**导出即用量**：本模块只保留三个程序实际引用的符号。历史上另有
`col` / `scol` / `load_csv_rows` / `last24_mean` / `series_mean` / `series_max` 六个导出，
因**全仓零引用**（响应面聚合在 Fortran 侧完成、逐时统计由 `CASE-analysis/run.py` 自己做）已删除；
需要末 24 行均值时直接 `read_table(...)[-24:]`，不要为此重新加回库函数。

---

## 3. `indicators.py` —— 12 量指示因子引擎

机制感知：物种表、NOy 权重、碳数全部从 `.fac` 的 `VARIABLE`/`FORMULA` 段推导，
所以也是 `fac_batch_omp` 里那 24 列指示因子的 **Python 参照实现**（两边互为核对）。

| API | 作用 |
|---|---|
| `_parse_fac(fac_path)` | 解析 `.fac`：物种表 / `FORMULA` 段（N→NOy 权重、C→碳数）/ 反应表 |
| `class Mechanism` | 物种索引、`noy_w`、`carbon_w`、反应计量；供指标计算用 |
| `Mechanism.p_o3_gd(peff_row)` | 论文 Eq(4)–(6) 的 O₃ 收支：G（HO₂/RO₂+NO）、D（五类汇）、P=G−D |
| `route_indicators(pref, fac, mech=None, pm_only=False, dtout=3600)` | post-spinup 全窗均值指标 → dict（12 个）；**无在役调用者** |
| `route_timeseries(pref, fac, mech=None)` | post-spinup 全部输出行的逐时指标序列；**无在役调用者** |

⚠ `pref` 的首列 `h` 是**输出步序号**（fac_v3 写 step）；`pm_only=True` 的午后窗按
`h × dtout/3600 % 24` 折算，故 `dtout` 必须由调用方从 runtime 的 MODEL 行传入。

**NOz 只有一个定义**：`Mechanism.nox_res = ("NO","NO2","N2O5")`，
`P(NOz)` 的反应白名单与逐时的 `NOy−NOx−N2O5` 都由它导出，并与 `fac_batch_omp` 的 NOz 口径一致。
| `PEROXY_SPECIES` | Eq(5) 的 Σ RO₂+NO 物种白名单（**机制改版新增 RO₂ 必须同步更新**；`p_o3_gd` 启动时校验物种存在） |

### 12 个量分层（官方 5 族是唯一的判据集）

| 类别 | 指标 | 判读 |
|---|---|---|
| **官方 5 族**（进判据/图件/矩阵） | `FNR` = HCHO/NO₂ | <1 VOC-limited，>2 NOx-limited（Sillman 1995） |
| | `H2O2_HNO3` | >0.2 → NOx-limited |
| | `L_NOx_Q` | >0.5 → NOx 饱和（Kleinman 2005） |
| | `P_H2O2_HNO3` | 0.15–0.5（HOx 汇分配） |
| | `OPE` | 低 → 饱和 |
| 状态量 | `NOy_ppb`（含氮物种加和） | 跨情景可比 / 与观测 NOy 对照 |
| 备用扩展（仍计算，不进判据与图件） | `O3_NOy` / `HCHO_NOy` / `O3_NOz` | 与官方 5 族信息重叠 |
| 绝对量 | `Q`（HOx 总循环速率）/ `L_NOx` / `O3_ppb` | 自由基活性 / NOx 汇 / O₃ |

**机制感知的好处**：NOy 定义从 `FORMULA` 段自动推导（N₂O₅ 计双 N），Python 侧加物种零改动。

### 使用示例

```python
import sys; sys.path.insert(0, 'analysis_lib')
from indicators import route_indicators
d = route_indicators('CASE-analysis/output/case', 'mechanism/mechanism_3_full_chain.fac')
print(f"O3={d['O3_ppb']:.1f} ppb  FNR={d['FNR']:.2f}  L_NOx/Q={d['L_NOx_Q']:.2f}")
```

（前缀 = 各程序 `output/` 里 `<name>_conc.csv` 的 `<name>`。）

---

## 4. 不放在这里的东西

| 不在 | 在哪 |
|---|---|
| 绘图 | 各程序目录的 `plot.py`；产物写各程序自己的 `output/`（一图一程序，无共享图库） |
| 编排 / 编译 / 跑模拟 | 各程序目录的 `run.py` |
| 机制与系数审计 | `mechanism/audit/`（element_audit / mcm_compare / jslot_audit） |

---

## 5. 改动风险

- **NOy 权重有两份**：`indicators.py` 从 `FORMULA` 段推导，`fac_batch_omp.f90` 的
  `noy_add()` 硬编码表在 [fac_batch_omp.f90:85](../facsimile_impl/fac_batch_omp.f90) 起。
  机制加含氮物种**两边都要动**，改后人工比对两侧的物种清单与权重（无自动守卫断言）。
- `PEROXY_SPECIES`、`ALKENE_SPECIES`、`nox_res`（NOz 排除集）同理是 Python 侧的硬编码白名单；
  漏更新会让 RIR / O₃ 收支的 G、D 项偏低（`PEROXY_SPECIES` 在 `Mechanism.__init__` 有存在性校验，
  另两个没有）。
- 两个文件必须保持**零第三方依赖**（只用 numpy），保证任何裸 Python 环境可跑。
