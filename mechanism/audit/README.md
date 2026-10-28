# mechanism/audit/ —— 机制质检工具（与仿真解耦，可独立运行；2026-10-27 自仓库根移入 mechanism/，归属机制设计层）

> **一句话**：这里回答"**机制本身对不对**"——元素守恒、系数与 MCM 官方基准的一致性、光解槽位语义。
> 三个脚本都**不跑仿真**（秒级），都已接入 `make test`，是机制变更后的第一道闸门。

---

## 1. 三个脚本

| 脚本 | 查什么 | 基准 | 耗时 |
|---|---|---|---|
| `element_audit.py` | **逐反应元素守恒**：C/H/N/O/S/Cl/Br 两侧配平 | 机制自带 FORMULA 段（单一事实源） | 0.07 s |
| `mcm_compare.py` | **速率系数 vs MCM v3.3.1 官方导出** + 文档漂移 | `mechanism/reference/mcm_v3.3.1_export.fac` | ~3 s |
| `jslot_audit.py` | **光解槽位语义**：机制 `J<n>` 的分子归属 vs MCM 官方该槽归属 + JMAP/README 一致性 | 同上 + `configs/photolysis.config` | ~3 s |

## 2. `element_audit.py` —— 元素守恒审计

```bash
python mechanism/audit/element_audit.py mechanism/mechanism_3_full_chain.fac   # → 203/203 PASS
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
python mechanism/audit/mcm_compare.py              # 详细报告
python mechanism/audit/mcm_compare.py --quiet      # 仅结论（make test 用）
python mechanism/audit/mcm_compare.py --threshold 0.05
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

- **光解**：`J<n>` 依赖光解表，无法在两机制间直接求值比对 → 本工具不覆盖，
  由 `jslot_audit.py` 按**槽位语义**单独把关（见 §4）。
- **产物分支结构**：只比总速率，不比分支拆分是否合理。

**退出码**：0=全部通过，1=超阈值或有文档漂移（接入 `make test`）。

## 4. `jslot_audit.py` —— 光解槽位语义审计（第三维度）

```bash
python mechanism/audit/jslot_audit.py              # 详细报告（逐槽判定表）
python mechanism/audit/jslot_audit.py --quiet      # 仅结论（make test 用）
make jcheck                                        # 等价于详细模式
```

**为什么单独立一维**：光解的 `l/m/n` 完全由**槽号**决定，而槽号是 MCM 官方命名空间里的
"外键"——借错槽 = 借了别的分子的 `l/m/n`。`mcm_compare.py` 比的是速率数值，`J` 值
根本不参与求值，因此对错位**完全无感**。V2.5 曾因此迁移 7 处（2026-10 审计发现 6 处
语义错位，记录于 `mechanism/V25_ASSESSMENT.md` §2.2）；本脚本把那次人工逐槽核对固化成断言。

### 四条不变量

| # | 不变量 | 抓什么错 |
|---|---|---|
| ① | **语义包含**：机制占用的每个 MCM 编号槽，其反应物（剔储库 + 集总映射后）⊆ MCM 官方该槽的分子集合 | 借错槽（如把 N₂O₅ 挂到 J14 = 丙醛槽） |
| ② | **自定义槽封闭**：MCM 无对应分子的槽必须显式落在白名单 `CUSTOM`（25/26/27/28/29/30/45），且 MCM 该槽确实为空，且 JMAP 给了 `D<正午值>` | 悄悄借用 MCM 已占槽；自定义槽漏参数 |
| ③ | **JMAP 卫生**：每个 JMAP 槽都被机制引用（无死参数），且只允许给白名单里的自定义槽 | 陈旧 JMAP 条目；手写值悄悄覆盖 MCM 官方 `l/m/n`（V2.5 前 J13 就是这样偏快 4.4×） |
| ④ | **文档一致**：`mechanism/README.md` 每条 `J<n>` 行（§5 反应表 + §6 光解方案表）的槽号/反应物/产物/参数来源 == `.fac` 实际 | 迁槽后忘改文档（实测抓出 2 处 V2.5 遗留） |

### 判定口径（防假阳/假阴的四个设计）

1. 只取 MCM 中速率 token **恰为** `J<n>`（无产率乘子）的反应作为该槽官方归属——
   MCM 大量用 `J<4>*0.14*0.4` 把 J4 的**日变化形状**借给别的分子，算进来①会退化成永真
2. 反应物侧的 O₂/H₂O/CO₂/H₂/N₂/M 是本机制的记账储库（MCM 不写在反应物侧），
   比对前剔除；储库表与 `element_audit.py` 共用同一份
3. 集总名 → MCM 代表分子复用 `mcm_compare.LUMP_MAP`（单一事实源），另加光解专属代理
   `ISOPOOH → ISOPBOOH`（MCM 六个异戊二烯-OOH 异构体全走 J41）
4. J34/J41 是 MCM 的"通用槽"（分别归属 105 / 1201 个分子），①在那里天然宽松；
   详细报告打印该槽分子数，让读者知道这条断言的实际强度

**当前结果**：26 条光解 / 25 槽 —— 18 个 MCM 官方槽逐槽核对 **0 错位**，
7 个自定义槽白名单封闭，README 52 行一致。

**负向验证**（沙箱改 `.fac`/`photolysis.config` 实测，仓库文件未动）：N₂O₅ 挪回 J14 →
5 处 FAIL 并指名"N₂O₅ 在 MCM 无光解"；ClO 挪到 MCM 已占的 J24 → FAIL；JMAP 加死槽 J22 →
FAIL；JMAP 覆盖官方槽 `13=D4.0e-5`（V2.5 前的真实错法）→ ③ 与 ④ 双双 FAIL；
README §6 槽号改错 → FAIL 并回指 `.fac` 实际槽号。

**退出码**：0=全过，1=有 FAIL，2=缺输入文件。

## 5. 与其他目录的分工

| 问题 | 去哪 |
|---|---|
| "反应式配平吗？" | **本目录** `element_audit.py` |
| "速率和 MCM 一致吗？" | **本目录** `mcm_compare.py` |
| "光解槽号借对了吗？参数有死项吗？" | **本目录** `jslot_audit.py` |
| "简化掉了什么、影响多大？" | `mechanism/V25_ASSESSMENT.md`（审计+评估总报告） |
| "输出 CSV 数值对吗？" | 各程序 `output/run_meta.txt` + 实测复跑（原 `tests/` 已随通用层退役） |
| "六账闭合吗？" | `CASE-analysis/run.py` + FT-IPR |

## 6. 何时跑

| 场景 | 命令 |
|---|---|
| 改了 `.fac`（任意原因） | `make test`（三个都会自动跑） |
| 只想快速验配平 | `python mechanism/audit/element_audit.py mechanism/mechanism_3_full_chain.fac` |
| 只想查系数一致性 | `make mcmcheck` |
| 加了/迁了光解反应，或改了 `photolysis.config` JMAP | `make jcheck` |
| 换了 MCM 官方导出（新版本） | 替换 `mechanism/reference/` 下的文件后 `make mcmcheck && make jcheck` |
