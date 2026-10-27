# mechanism/audit/ —— 机制质检工具（与仿真解耦，可独立运行；2026-10-27 自仓库根移入 mechanism/，归属机制设计层）

> **一句话**：这里回答"**机制本身对不对**"——元素守恒、系数与 MCM 官方基准的一致性。
> 两个脚本都**不跑仿真**（秒级），都已接入 `make test`，是机制变更后的第一道闸门。

---

## 1. 两个脚本

| 脚本 | 查什么 | 基准 | 耗时 |
|---|---|---|---|
| `element_audit.py` | **逐反应元素守恒**：C/H/N/O/S/Cl/Br 两侧配平 | 机制自带 FORMULA 段（单一事实源） | 0.07 s |
| `mcm_compare.py` | **速率系数 vs MCM v3.3.1 官方导出** + 文档漂移 | `mechanism/reference/mcm_v3.3.1_export.fac` | ~3 s |

## 2. `element_audit.py` —— 元素守恒审计

```bash
python audit/element_audit.py mechanism/mechanism_3_full_chain.fac   # → 198/198 PASS
```

### 设计要点（防"假绿"的三道内建规则）

| 规则 | 防的错 |
|---|---|
| 储库物种（H₂O/O₂/CO₂/H₂/M/N₂）**自动推断**闭合 | 免去手写 100+ 条储库系数 |
| 但**限定物理方向**：CO₂/H₂ 只许产物侧、O₂ 反应物侧 ≤2.5 | 防"随便补个 CO₂ 把真错掩盖" |
| 分子式优先读机制 `FORMULA` 段（`[M3]` 标记），内嵌表仅兜底 | 单一事实源；改机制必须同步 FORMULA，否则报"未知物种" |

**边界**（诚实声明）：它验证**数学配平**，不验证**化学合理性**——
若 FORMULA 段里把某物种分子式写错（如 MACR 写成 C4H7O），审计会拿错式子去配平并 PASS。
分子式本身靠 `mechanism/README.md` 的物种表人工核对。

**退出码**：0=全过，1=有 FAIL（`make test` 依赖此码中断）。

## 3. `mcm_compare.py` —— MCM 系数比对（三道关口的合体）

```bash
python audit/mcm_compare.py              # 详细报告
python audit/mcm_compare.py --quiet      # 仅结论（make test 用）
python audit/mcm_compare.py --threshold 0.05
make mcmcheck                            # 等价于详细模式
```

### 三道关口

| 关口 | 方法 | 当前状态 |
|---|---|---|
| **① 无机反应数值** | 按"反应物集合"匹配 MCM → 同条件（T=298.15K, M=2.46e19）求值比对 | 32 条，最大偏差 8.9% |
| **② 集总物种映射** | 集总类映射到 MCM 代表分子（ISOP→C5H8、TERP→APINENE、SESQ→BCARY…共 18 个）再比对 | 16 条，最大 8.9% |
| **③ 文档漂移** | `mechanism/README.md` 反应表里的速率值 == `.fac` 实际值 | 89 条全一致 |

### 比对口径（防假阳性的四个设计）

1. **按反应物集合匹配，不按产物**——集总化导致产物必然与 MCM 不同
2. MCM 一反应多产率分支 → **求和后比较**（如 MACR+OH 的 0.45/0.47/0.08 三支）
3. 本机制同集合多条（HO₂+HO₂ 的 A/B 双通道）→ **先聚合再比**
4. 输出分两类：**exact**（产物一致，判数值）vs **结构性差异**（集总化取舍，单列不判 FAIL）

### 已知盲区（诚实声明）

- **光解**：`J<n>` 依赖光解表，无法在两机制间直接求值比对 → 本工具不覆盖。
  光解槽位语义需对照 MCM `J<n>` 反应式人工核查
  （2026-10 审计发现 6 处语义错位，已修复并记录于 `mechanism/V25_ASSESSMENT.md`）。
- **产物分支结构**：只比总速率，不比分支拆分是否合理。

**退出码**：0=全部通过，1=超阈值或有文档漂移（接入 `make test`）。

## 4. 与其他目录的分工

| 问题 | 去哪 |
|---|---|
| "反应式配平吗？" | **本目录** `element_audit.py` |
| "速率和 MCM 一致吗？" | **本目录** `mcm_compare.py` |
| "简化掉了什么、影响多大？" | `mechanism/V25_ASSESSMENT.md`（审计+评估总报告） |
| "输出 CSV 数值对吗？" | `tests/functional_test.py`（金值哨兵+契约断言） |
| "六账闭合吗？" | `experiments/run/run_D_ipr.py` + FT-IPR |

## 5. 何时跑

| 场景 | 命令 |
|---|---|
| 改了 `.fac`（任意原因） | `make test`（两个都会自动跑） |
| 只想快速验配平 | `python audit/element_audit.py mechanism/mechanism_3_full_chain.fac` |
| 只想查系数一致性 | `make mcmcheck` |
| 换了 MCM 官方导出（新版本） | 替换 `mechanism/reference/` 下的文件后 `make mcmcheck` |
