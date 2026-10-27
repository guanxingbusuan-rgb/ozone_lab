# experiments/ —— 实验层（run 入口 / post 诊断 / viz 绘图）

> **配套关系**：本目录是 `facsimile_impl/` 两个积分执行器（fac_v3 / fac_batch_omp）
> 的实验入口层——`run/` 调驱动 CLI、只跑模拟+落 CSV 到 `output/`；
> `post/` 吃 CSV 出诊断；`viz/` 出图。旧入口（ridge5/scen5/ipr_irr5/hourly_pa）
> 已于 2026-10-27 **整目录删除**——功能经逐项核查由 `run/` 全覆盖（原始脚本可查 git 历史）。

| 子目录 | 职责 | 内容 |
|---|---|---|
| **`run/`** | **实验入口（唯一现役）** | `run_A..run_M` + `run_all.py` + `_lib.py`（驱动封装/clone_input 输入变体/manifest 溯源）；详见 [run/README.md](run/README.md) |
| **`post/`** | 诊断输出（吃 output/ CSV） | `ipr_closure_check.py`（IPR 六账完备性校验） |
| **`viz/`** | 绘图（**已入库**） | `fig_response_v3.py` / `fig_indicators_ridge.py`（`make figures` 调用）、`fig_ridge5.py`；图输出到本地 `figures/`（不入库） |

> 🚫 **伴随/反演已移除（2026-10-27）**：0-D 箱体无编译器 AD、无切线性模式——真伴随
> 不可实现；FD 雅可比 + 截断 SVD 实测仅 1/9 模态可辨识。功能与载体 `sens_adj.py`
> 一并删除；OAT 敏感性（不做反演）由 `run/run_H_sensitivity.py` 承担，
