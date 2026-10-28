# analysis_lib/ —— 诊断与指标库（Python 共用函数库）（纯 Python，读 CSV 输出，不参与仿真）

> **一句话**：仿真是 `fac_v3` / `fac_batch_omp` 的事；**这里负责"仿真跑完之后"的一切数值消费**——
> 读 `_conc.csv` / `_irr.csv` / `_ipr.csv`，算指示因子、整理 IPR 六账、兜底解析 Fortran 的怪数值格式。
> 三个文件都是**被 import 的库**，没有独立入口脚本。

---

## 1. 三个文件

| 文件 | 行数 | 定位 | 谁在用 |
|---|---|---|---|
| `indicators.py` | ~250 | **12 量指示因子引擎**（面向对象，机制感知；官方 5 族分层） | `run_E_indicators.py`、`CASE-analysis/run.py`/`run_G_consist.py`、`run_B` 的机制解析 |

> `experiments/run/_runner.py` 是**编排层**（跑仿真、spawn 二进制），不属于本目录。
> 依赖方向恒为 `run_* → _runner → (csvio) → 二进制`，analysis_lib 不反向依赖任何编排代码。

---


解决一个真实踩过的坑：gfortran 的 `ES13.5` 等格式在 `|指数|≥100` 时**静默丢掉 E**，
输出形如 `3.07932-174` 的数——Python `float()` 直接抛异常。

| 函数 | 作用 |
|---|---|
| `fix_fortran_exp(tok)` | 把 `3.07932-174` 修成 `3.07932E-174`；正常格式原样返回 |
| `fnum(tok)` | = `float(fix_fortran_exp(tok))`，一行解析 |
| `load_csv_rows(path)` | 读 CSV → 二维 list（已过 `fnum`） |

> 现行输出格式已改 `ES15.5E3`（恒三位指数），新文件不会再出现此问题；
> 但**历史产物**与某些边界值仍可能触发，所以解析一律走这里，不要裸 `float()`。

---


| 函数 | 作用 | 输入 |
|---|---|---|
| `load(prefix)` | 读 `<pref>_conc.csv` + `<pref>_irr.csv` → `(h, C, irr_h, I)` 四个 numpy 数组 | fac_v3 输出前缀 |
| `ipr_from_fac(prefix)` | 读 `<pref>_ipr.csv` → `{物种: {P,L,E,D,B,X}}` 字典（**窗口积分值，勿再 trapz**） | fac_v3 输出前缀 |

设计要点：**薄而稳**。只做"读 + 转 numpy"，所有科学计算在上层（`run_D_*` / `run_G_consist` / 测试）。
`ipr_from_fac` 曾有"对已积分值再做 trapz"的双计 bug（R10），现在函数内注释了警告。

---

## 4. `indicators.py` —— 12 量指示因子引擎（机制感知 + 官方 5 族分层）

这是 `run_E`（指示因子工具箱）和 `CASE-analysis`/`run_G_consist`（反应级归因/互证）的计算核心，
也是 `fac_batch_omp` 里 24 列指示因子的**Python 参照实现**（两边互为校验）。

| API | 作用 |
|---|---|
| `_parse_fac(fac_path)` | 解析 `.fac`：物种表 / **FORMULA 段**（N→NOy 权重，C→碳数）/ 反应表 |
| `class Mechanism` | 封装机制：物种索引、NOy 权重、碳数、反应计量；供指标计算用 |
| `Mechanism.p_o3_gd(peff_row)` | 论文 Eq(4)–(6) 的 O₃ 收支：G（HO₂/RO₂+NO）、D（五类汇）、P=G−D；供 D2 预算与 D3 RIR 共用 |
| `route_indicators(pref, fac, mech=None, pm_only=False)` | post-spinup 全窗均值指标 → dict（12 个） |
| `route_timeseries(pref, fac, mech=None)` | post-spinup 全部输出行的逐时指标序列 |

### 12 个量分层（官方 5 族是唯一的判据集）

| 类别 | 指标 | 判读 |
|---|---|---|
| **官方 5 族**（进判据/图件/矩阵） | `FNR` = HCHO/NO₂ | <1 VOC-limited，>2 NOx-limited（Sillman 1995） |
| | `H2O2_HNO3` | >0.2 → NOx-limited |
| | `L_NOx_Q` | >0.5 → NOx 饱和（Kleinman 2005） |
| | `P_H2O2_HNO3` | 0.15–0.5（HOx 汇分配） |
| | `OPE` | 低 → 饱和 |
| 状态量 | `NOy_ppb`（25 含氮物种加和） | 跨情景可比 / 与观测 NOy 对照 |
| 备用扩展（仍计算，不进判据与图件） | `O3_NOy` / `HCHO_NOy` / `O3_NOz` | 与官方 5 族信息重叠，历史包袱 |
| 绝对量 | `Q`（HOx 总循环速率）/ `L_NOx` / `O3_ppb` | 自由基活性 / NOx 汇 / O₃ |

**机制感知的好处**：NOy 定义从 FORMULA 段自动推导（N₂O₅ 计双 N），
机制加含氮物种**零改动**接入——`FT-SYNC` 守卫断言它与 `fac_batch_omp` 硬编码表逐项一致。

### 使用示例

```python
import sys; sys.path.insert(0, 'analysis_lib')
from indicators import route_indicators
d = route_indicators('output/run_A/BAU', 'mechanism/mechanism_3_full_chain.fac')
print(f"O3={d['O3_ppb']:.1f} ppb  FNR={d['FNR']:.2f}  L_NOx/Q={d['L_NOx_Q']:.2f}")
```

---

## 5. 不放在这里的东西

| 不在 | 原因 |
|---|---|
| 绘图 | 工程图脚本在 `experiments/viz/`，产物在 `output/`；机制图静态归档在 `docs/figures/` |
| 实验编排 | 在 `experiments/run/` |
| 机制/系数审计 | 在 `mechanism/audit/`（element_audit / mcm_compare / jslot_audit） |


## 6. 改动风险

- `indicators.py` 的 NOy 自动判定：改了会破坏与 `fac_batch_omp` 权重表的一致性（守卫会红）。
- 三个文件都必须**零第三方依赖**（只用 numpy），保证任何裸 Python 环境可跑。
