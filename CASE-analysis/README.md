# CASE-analysis/ —— 个例案例分析（六账 + 反应级归因 + 前体物 RIR）

> **一句话**：一个「个例」= 本目录 `config_override/` 的一份完整配置快照；
> `run.py` 跑基线模拟并导出全部分析 CSV，`plot.py` 只画图。
> 取代原 `experiments/run/` 下的 `run_D1_ipr.py`（六账）、`run_D2_irr.py`（反应级 IRR）、`run_D3_rir.py`（RIR）三个入口。

---

## 1. 架构（与 EKMA-emission / EKMA-FNR 同构，自完备）

```
CASE-analysis/
├── run.py              # 唯一入口：配置快照 → cfgbuild → fac_v3 基线 → 导出全部 CSV
├── plot.py             # 只读 output/ 的 CSV，出 3 张图（无数据导出）
├── README.md
├── config_override/    # ★ 输入：完整配置集（气象/初值/排放/沉降/情景/约束…）
└── output/             # 存档：数据 + 日志 + 元数据 + 图
    ├── in_case/                 配置快照（cfgbuild 的真实输入，可追溯）
    ├── runtime_case.dat         本个例专用 runtime
    ├── case_{conc,ipr,irr}.csv  基线 fac_v3 原始输出（post-spinup 全部步）
    ├── ipr_ledger.csv           六账长表 h×物种×(P L E D B X)
    ├── irr_top.csv              反应级净贡献 top-12（默认 O3）
    ├── hourly_rxn_top.csv       逐小时净贡献 top-4
    ├── rir_mean.csv             平均 RIR（Eq8，P(O3) 加权）
    ├── rir_diurnal.csv          逐时 RIR_avg
    ├── p_o3no_baseline.csv      基线 O3 收支（Eq4–6）
    ├── rir_runs/                RIR 扰动的配置快照 + runtime + 原始输出
    ├── runs.csv                 基线与扰动运行记录
    ├── run.log / run_meta.txt   日志与元数据
    └── case_{ipr,irr,rir}.png   图件
```

**自包含边界**：不 import `experiments/run/`；共享仅"模型本体"——
`facsimile_impl/fac_v3`（执行器）、`inputproc/cfgbuild.py`（编译器）、
`mechanism/*.fac`（机制）、`analysis_lib/{csvio,indicators}`（读写与 O3 收支）。

## 2. 用法

```bash
python CASE-analysis/run.py                        # 基线 + 六账 + IRR + RIR(-10%, 8–18h)
python CASE-analysis/run.py --no-rir               # 只做六账/IRR（快）
python CASE-analysis/run.py --pct 20 --groups NOx,CO --window 8-18
python CASE-analysis/run.py --species HCHO --per-species
python CASE-analysis/plot.py [--species O3]
```

换一个例：直接改 `config_override/`（完整配置集，**不回退根 `configs/`**），再跑 `run.py`。
RIR 的扰动通过复制快照并缩放 `emissions.config` 基线实现，产物与运行记录都在 `output/rir_runs/`。

## 3. 口径与判读

| 量 | 口径 |
|---|---|
| 六账 | P 机制产 / L 机制耗 / E 排放 / D 干湿沉 / B 同化+稀释 / X BLH 膨胀；单位 molec/cm³，**每行=本小时窗积分**（勿再 trapz） |
| IRR 净贡献 | Σ_k coef_k(SP)·I_k，产物 +系数、反应物 −1；窗口积分直接求和 |
| RIR | Eq(7) 为内部量不导出；**只导出 Eq(8) 平均 RIR**（P(O3) 加权） |
| P(O3) | Eq(4–6)：G−D（`indicators.Mechanism.p_o3_gd`，与 `experiments/post/o3_budget.py` 同源） |
| S（前体物） | 模拟浓度：NOx=[NO]+[NO₂]、CO、VOC 组按碳数（ppbC）；ΔS=基线−扰动，减排为正 |

RIR>0 表示「减排该前体物会降低 P(O3)」；窗口内 P≤0 或 S 非正的样本自动剔除。
