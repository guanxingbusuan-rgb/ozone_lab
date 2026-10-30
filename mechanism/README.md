# mechanism_3_full_chain.fac — V2.5+ 机制规格书（冻结版）

> **ozone_lab 工程机制的唯一事实源**。本文档是该 `.fac` 文件的**逐项完整规格说明**：
> 设计定位、冻结基线、111 物种全表（含分子式与角色）、203 条反应全表（含速率常数与
> 参数来源）、光解方案全表（26 条光解反应 / 槽号至 45）、动力学约定、元素配平方法论、
> 已知省略、对比定位、单位约定、扩展指南与参考文献。
>
> | 项 | 值 |
> |---|---|
> | 版本 | **V2.5+ 冻结版**（2026-10-27，不再迭代） |
> | 规模 | **111 物种 / 203 反应 / 45 个光解槽位 / 127 个命名速率** |
> | 元素审计 | **203/203 PASS**（`mechanism/audit/element_audit.py`）；**全部反应精确平衡** |
> | MCM 一致性 | **48 条可比对，最大偏差 8.9%**（32 无机 + 16 集总映射；34 条产物一致者偏差 0.0%） |
> | 光解槽语义 | **0 错位**（`mechanism/audit/jslot_audit.py`：18 个 MCM 官方槽逐槽核对 + 7 个自定义槽白名单封闭） |
> | 生成码 | `facsimile_impl/gen_v3/oz_mech.f90`（`mechc.py` 编译） |
> | 承载路线 | `fac_v3`（单箱深度）+ `fac_batch_omp`（批量广度），共享同一生成码 |
> | 文件本体 | [`mechanism_3_full_chain.fac`](mechanism_3_full_chain.fac)（FACSIMILE 语法） |
> | 冻结状态 | ✅ **V2.5+ 已冻结**（2026-10-27）——P0 级正确性问题已全部修复；
> 评估与 MCM 审计结论见 [V25_ASSESSMENT.md](V25_ASSESSMENT.md) |

---

## 目录

- [1. 设计定位与适用范围](#1-设计定位与适用范围)
- [2. V2.5+ 冻结基线](#2-v25-冻结基线)
- [3. 反应段总览](#3-反应段总览)
- [4. 物种全表（111，含分子式）](#4-物种全表)
- [5. 反应全表（203 条，含速率）](#5-反应全表)
- [6. 光解方案（26 条光解 / 槽号至 45）](#6-光解方案)
- [7. 动力学约定](#7-动力学约定)
- [8. 元素配平方法论](#8-元素配平方法论)
- [9. 已知省略（工程决策）](#9-已知省略)
- [10. 对比定位](#10-对比定位)
- [11. 单位约定](#11-单位约定)
- [12. 使用与扩展指南](#12-使用与扩展指南)
- [13. 参考文献](#13-参考文献)

---

## 1. 设计定位与适用范围

### 1.1 它是什么

一个**中等规模的对流层臭氧光化学机制**：0-D 箱体模型用的反应集合，覆盖
NOx–O₃–HOx–ROx 三大循环 + 生物源/人为 VOC 降解 + 夜间化学 + 硫链 + 卤素。
规模定位在 **RADM2（63sp/157rx）与 CB6R3（82sp/221rx）之间**，VOC 侧采用
"显式裁剪（CH₄/异戊二烯/单萜）+ 集总类（烯烃/烷烃/芳烃/酮）"的混合表示。

### 1.2 它能回答什么 / 不能回答什么

| 能（设计目标） | 不能（超出定位） |
|---|---|
| 三大循环结构与 O₃ 净生成通道的定量分解（P(O₃) 两通道） | 城市 O₃ 浓度预报（无输运/网格/在线气象） |
| 昼夜化学的完整日循环（含 HONO 晨爆、NO₃/N₂O₅ 夜化学） | 全谱 VOC 响应（无高级醛/甲酚/倍半萜/OVOC 长尾） |
| 排放敏感性（E_KMA 屋脊）与体制判别（FNR/L_NOx·Q⁻¹/OPE） | 气溶胶耦合（SOA 仅以 TERPCHO 代理 + 沉降汇） |
| 观测约束下的机制检验（BC 张弛同化 / OSSE） | 异相化学细节（N₂O₅ 水解用一阶参数化） |
| 情景沙盘（NOx/VOC/温度/云量/IC 因子） | 卤素-气溶胶界面化学（无海盐气溶胶显式） |
| 过程归因（IPR 六本账，逐物种逐过程） | 平流层化学、汞化学、IVOC/SVOC |

### 1.3 设计原则

1. **骨架优先**：每条 VOC 链必须有完整的"引发 → 过氧 → 分支 → 羰基 → 自由基回收/储库"脊柱，不允许断链。
2. **元素严格配平**：每反应 C/H/N/O/S/Cl/Br 闭合（审计器强制），分数化学计量系数允许。
3. **速率可溯源**：JPL/IUPAC 推荐值优先；MCM v3.3.1 逐字拷贝的部分在 `.fac` 注释中标注 `[MCM]`。
4. **复杂度花在 VOC 侧**：无机组从 RADM2 到 MCM 几乎逐字同源，因此本机制的差异化投入全部在 VOC 降解谱。
5. **工程可裁剪**：新链以独立段追加，不改动既有段；撤回只需删段。

---

## 2. V2.5+ 冻结基线

> 本文档只描述 **V2.5+ 冻结版**；历史设计段落不再保留。
> 需要追溯版本沿革时查 git 历史与 [V25_ASSESSMENT.md](V25_ASSESSMENT.md)。

| 项 | 值 |
|---|---|
| 机制规模 | **111 物种 / 203 反应 / 127 命名速率 / 45 光解槽位** |
| 元素审计 | **203/203 PASS**（C/H/N/O/S/Cl/Br 逐反应配平） |
| MCM 对照 | **48 条可比对，最大偏差 8.9%**（32 无机 + 16 集总映射；34 条产物一致者 0.0%） |
| 生成码 | `facsimile_impl/gen_v3/oz_mech.f90` + `oz.reac/oz.prod` 计量表 |
| 承载路线 | `fac_v3`（单箱深度）+ `fac_batch_omp`（批量广度），共享同一生成码 |
| 冻结规则 | 任何机制设计变更 = 开启 V2.6；日常实验只改各程序的 `config_override/`，不改 `.fac` |

### 2.1 扩物种必须同步的 8 处

> 加一个物种不是“只改 .fac”。以下 8 处任一漏掉，都会产生**静默错误**
> （不报错、但结果错）。扩展时逐项验证。

| # | 位置 | 做什么 | 漏掉的后果 | 当前状态 |
|---|---|---|---|---|
| 1 | `mechanism/*.fac` | `VARIABLE` 加物种、`FORMULA` 加分子式、写反应 | 编译失败或 ns 错位 | ✓ 111/203 |
| 2 | `<程序>/config_override/emissions.config` | 若有**一次排放**，加基线+日型 | **新物种恒为 0**，链条空转 | ✓ 25 物种 |
| 3 | `<程序>/config_override/drydep.config` | 若无**气相化学汇**，必须给 vd | **无界累积**（死端） | ✓ 已配 |
| 4 | `modules/runtime_dat_mod.f90` `build_groups()` | 若属 NOx/VOC/CO/CH₄ 排放组 | 情景缩放不生效（HVOC 白做） | ✓ 已配 |
| 5 | `facsimile_impl/fac_batch_omp.f90` `noy_add()` | 若含氮（NOy 家族） | 批量路线**逐案例守恒门直接 `stop 1`**；单箱路线 NOy 指标偏低 | ✓ 25 物种 |
| 6 | `modules/chem_phys_mod.f90` `tune_atol_by_name()` | 若浓度 >1e9（约 40 ppb） | 步数爆炸（精度过剩） | ✓ 已配 |
| 7 | 参考基线（`ENVIRONMENT.md` §8） | 机制改动必然移动输出量级 | **无数值哨兵**（`tests/` 已退役）——漂移无人告知 | ✓ 3 程序实测 |
| 8 | 文档计数 | 本 README + `ENVIRONMENT.md` 的物种/反应计数 | 文档与实际不符 | ✓ 111/203 |

**补正记录**：8 处清单之外还有一类易漏项——“稳定产物缺物理汇”。
`drydep.config` 曾漏配 GLY/CH3OOH/ALKOOH/PNA 四个稳定二次产物，
实测导致 GLY 高估 32%、ALKOOH 49%。已补齐；死端物种由下面的脚本手工筛查
（无数值守卫，扩物种后建议跑一次）。

**两个不需要改的地方**：
- `analysis_lib/indicators.py` 的 NOy 定义**自动**从 FORMULA 段读含 N 物种
  （当前识别 25 个，与 batch 硬编码表逐项一致）——单一事实源的收益。
- `config_override/initialConcentrations.config` 与 `bc.config` **不需要**加新物种：
  二次产物从 0 起步是正确的；背景空气里不该有这些一次 VOC 的边界浓度。

**自动检测死端的脚本**（扩物种后建议跑一次）：

```bash
python3 - <<'EOF'
import re
t = open("mechanism/mechanism_3_full_chain.fac").read()
sp = re.search(r"(?m)^VARIABLE\s*\n(.*?);", t, re.S).group(1).replace(";", " ").split()
gen = "facsimile_impl/gen_v3"
lost = set()
for row in open(gen + "/oz.reac"):
    lost |= {int(x) for x in row.split()[1:]}
em = {sp.index(l.split()[0]) + 1 for l in open("CASE-analysis/config_override/emissions.config")
      if l.split("#")[0].split() and l.split()[0] in sp}
dry = {sp.index(l.split()[0]) + 1 for l in open("CASE-analysis/config_override/drydep.config")
       if len(l.split("#")[0].split()) >= 3 and l.split()[2] == "on"
       and l.split()[0] in sp}
dead = [sp[i-1] for i in range(1, len(sp)+1) if i not in lost and i not in em and i not in dry]
print("零汇且无排放无干沉:", dead or "无")
EOF
```

---

## 3. 反应段总览

V2.5+ 共 **203 条反应**，按化学链分为 10 条主链 + 7 个扩展段：

| 段 | 名称 | 反应数 |
|---|---|---|
| ① | 无机核心：NOx 光稳态 + HOx + HONO 晨爆 | 26 |
| ② | NO₃/N₂O₅ 夜化学（MCM Troe） | 10 |
| ③ | CH₄ 氧化链 + PAN/CH₃CO₃ 储库循环 | 16 |
| ④ | ISOP 完整二代氧化链 | 24 |
| ⑤ | 硫链 DMS→SO₂→SO₃→H₂SO₄ | 4 |
| ⑥ | 卤素催化 O₃ 损耗 + 储库闭环 | 19 |
| ⑦ | VOC 降解谱集总类 | 21 |
| ⑧ | 支撑化学（PNA / N₂O₅ 光解 / OH+HNO₃） | 3 |
| ⑨ | 烷基硝酸盐家族 | 3 |
| ⑩ | 单萜迷你链 | 9 |
| ⑪ | 芳香烃分级 + RO₂+HO₂ 闭环 | 10 |
| ⑫ | 长链烷烃 LALK | 6 |
| ⑬ | 醇/醚/酯类含氧 VOC | 14 |
| ⑭ | 倍半萜 SESQ | 8 |
| ⑮ | 乙烯/内烯烃/高级醛 | 24 |
| ⑯ | MACR 夜化学 | 1 |
| ⑰ | RO₂ 交叉缔合 | 5 |
| | **合计** | **203** |

---

## 4. 物种全表

**111 个演化物种**（`VARIABLE` 段，顺序即 Fortran 索引 1–111；`runtime_*.dat` 的 SPECIES 行同序）。
**过氧自由基池**（`.fac` 的 `RO2 =` 声明）：`CH3O2 CH3CO3 ISOPBO2 ISOPNO3 MVKO2 MACRO2`。

> 分子式来自 `.fac` 的 `FORMULA` 段（元素审计与 NOy 自动判定的单一事实源）。
> "角色"列：R=自由基，P=过氧自由基，S=储库，C=羰基/产物，E=排放初级物种，X=卤素。


#### 无机核心 (14)

| 序号 | 物种 | 分子式 | 角色 |
|---|---|---|---|
| 1 | `NO` | NO |  |
| 2 | `NO2` | NO2 |  |
| 3 | `O3` | O3 |  |
| 4 | `O1D` | O |  |
| 5 | `O` | O |  |
| 6 | `OH` | OH |  |
| 7 | `HO2` | HO2 |  |
| 8 | `H2O2` | H2O2 |  |
| 9 | `CO` | CO |  |
| 10 | `HCHO` | CH2O |  |
| 11 | `HNO3` | HNO3 |  |
| 12 | `HONO` | HNO2 |  |
| 13 | `NO3` | NO3 |  |
| 14 | `N2O5` | N2O5 |  |

#### CH₄ 链 (4)

| 序号 | 物种 | 分子式 | 角色 |
|---|---|---|---|
| 15 | `CH4` | CH4 |  |
| 16 | `CH3O2` | CH3O2 |  |
| 17 | `CH3OOH` | CH4O2 |  |
| 18 | `CH3O` | CH3O |  |

#### PAN / 醛 (3)

| 序号 | 物种 | 分子式 | 角色 |
|---|---|---|---|
| 19 | `CH3CHO` | C2H4O |  |
| 20 | `CH3CO3` | C2H3O3 |  |
| 21 | `PAN` | C2H3NO5 |  |

#### ISOP 一代 (7)

| 序号 | 物种 | 分子式 | 角色 |
|---|---|---|---|
| 22 | `ISOP` | C5H8 |  |
| 23 | `ISOPBO2` | C5H9O3 |  |
| 24 | `ISOPBO` | C5H9O2 |  |
| 25 | `ISOPNO3` | C5H8NO5 |  |
| 26 | `ISOPNIT` | C5H9NO4 |  |
| 27 | `ISOPDNT` | C5H9NO5 |  |
| 28 | `ISOPOOH` | C5H0O3 |  |

#### ISOP 二代 MVK/MACR (10)

| 序号 | 物种 | 分子式 | 角色 |
|---|---|---|---|
| 29 | `MVK` | C4H6O |  |
| 30 | `MVKO2` | C4H7O4 |  |
| 31 | `MVKO` | C4H7O3 |  |
| 32 | `MVKOOH` | C4H8O4 |  |
| 33 | `MACR` | C4H6O |  |
| 34 | `MACRO2` | C4H5O3 |  |
| 35 | `MACROO` | C4H5O2 |  |
| 36 | `MACROOH` | C4H6O3 |  |
| 37 | `MPAN` | C4H5NO5 |  |
| 38 | `MGLY` | C3H4O2 |  |

#### 硫链 (5)

| 序号 | 物种 | 分子式 | 角色 |
|---|---|---|---|
| 39 | `SO2` | SO2 |  |
| 40 | `SO3` | SO3 |  |
| 41 | `H2SO4` | H2SO4 |  |
| 42 | `DMS` | C2H6S |  |
| 43 | `MSA` | CH4O3S |  |

#### 卤素 (10)

| 序号 | 物种 | 分子式 | 角色 |
|---|---|---|---|
| 44 | `Cl` | Cl |  |
| 45 | `ClO` | ClO |  |
| 46 | `ClONO2` | ClNO3 |  |
| 47 | `Br` | Br |  |
| 48 | `BrO` | BrO |  |
| 49 | `BrONO2` | BrNO3 |  |
| 50 | `HCl` | HCl |  |
| 51 | `HBr` | HBr |  |
| 52 | `HOCl` | HOCl |  |
| 53 | `HOBr` | HOBr |  |

#### 集总类 (14)

| 序号 | 物种 | 分子式 | 角色 |
|---|---|---|---|
| 54 | `OLE` | C3H6 |  |
| 55 | `OLERO2` | C3H7O3 |  |
| 56 | `OLEOOH` | C3H8O3 |  |
| 57 | `ALK` | C3H8 |  |
| 58 | `ALKRO2` | C3H7O2 |  |
| 59 | `ALKOOH` | C3H8O2 |  |
| 60 | `ARO1` | C7H8 |  |
| 61 | `ARO1O2` | C7H9O3 |  |
| 62 | `ARO2` | C8H0 |  |
| 63 | `ARO2O2` | C8HO3 |  |
| 64 | `KET` | C3H6O |  |
| 65 | `KETRO2` | C3H5O3 |  |
| 66 | `GLY` | C2H2O2 |  |
| 67 | `PNA` | HNO4 |  |

#### 硝酸酯 + 单萜 (7)

| 序号 | 物种 | 分子式 | 角色 |
|---|---|---|---|
| 68 | `RONO2` | C3H7NO3 |  |
| 69 | `OLENIT` | C3H7NO4 |  |
| 70 | `TERP` | C0H6 |  |
| 71 | `TERPRO2` | C0H7O3 |  |
| 72 | `TERPOOH` | C0H8O5 |  |
| 73 | `TERPNIT` | C0H6NO3 |  |
| 74 | `TERPCHO` | C0H4O2 |  |

#### 扩展 VOC 谱 (37)

| 序号 | 物种 | 分子式 | 角色 |
|---|---|---|---|
| 75 | `CSL` | C7H8O1 | 甲酚（芳香开环关键中间体） |
| 76 | `CSLRO2` | C7H7O3 | 甲酚过氧自由基 |
| 77 | `CSLOOH` | C7H8O5 | 甲酚氢过氧化物 |
| 78 | `CSLNIT` | C7H7NO4 | 甲酚硝酸酯（NOy） |
| 79 | `TOLOOH` | C7H10O5 | 甲苯 RO₂+HO₂ 产物（原缺失） |
| 80 | `XYLOOH` | C8H12O5 | 二甲苯 RO₂+HO₂ 产物（原缺失） |
| 81 | `LALK` | C6H14 | C6 长链烷烃 |
| 82 | `LALKRO2` | C6H13O2 | 长链烷烃 RO₂ |
| 83 | `LALKOOH` | C6H14O4 | 长链烷烃氢过氧化物 |
| 84 | `LALKNIT` | C6H13NO3 | 长链烷烃硝酸酯（NOy） |
| 85 | `MEOH` | CH4O | 甲醇（二次 HCHO 上游） |
| 86 | `ETOH` | C2H6O | 乙醇 |
| 87 | `MTBE` | C5H12O | MTBE 醚 |
| 88 | `MTBERO2` | C5H11O3 | MTBE RO₂ |
| 89 | `MTBEOOH` | C5H12O5 | MTBE 氢过氧化物 |
| 90 | `MTBENIT` | C5H11NO4 | MTBE 硝酸酯（NOy） |
| 91 | `ESTER` | C4H8O2 | 酯类（乙酸乙酯代理） |
| 92 | `ESTERO2` | C4H7O4 | 酯 RO₂ |
| 93 | `ESTEROOH` | C4H8O6 | 酯氢过氧化物 |
| 94 | `ESTERNIT` | C4H7NO5 | 酯硝酸酯（NOy） |
| 95 | `SESQ` | C15H24 | 倍半萜（SOA 主力） |
| 96 | `SESQRO2` | C15H23O3 | 倍半萜 RO₂ |
| 97 | `SESQOOH` | C15H24O5 | 倍半萜氢过氧化物 |
| 98 | `SESQNIT` | C15H23NO4 | 倍半萜硝酸酯（NOy） |
| 99 | `SESQCHO` | C15H22O2 | 倍半萜 SOA 代理终产物 |
| 100 | `ETH` | C2H4 | 乙烯（全球最大排放烯烃） |
| 101 | `ETHO2` | C2H5O3 | 乙烯 RO₂（羟乙基过氧） |
| 102 | `ETHOOH` | C2H6O5 | 乙烯氢过氧化物 |
| 103 | `ETHNIT` | C2H5NO4 | 乙烯硝酸酯（NOy） |
| 104 | `ALD2` | C3H6O | 高级醛（丙醛代理，★HCHO 上游） |
| 105 | `ALD2O2` | C3H5O3 | 丙酰过氧 |
| 106 | `ALD2OOH` | C3H6O5 | 高级醛氢过氧化物 |
| 107 | `PPAN` | C3H5NO5 | 过氧丙酰硝酸酯（NOy） |
| 108 | `IOLE` | C4H8 | 内烯烃（丁烯代理） |
| 109 | `IOLERO2` | C4H9O3 | 内烯烃 RO₂ |
| 110 | `IOLEOOH` | C4H10O5 | 内烯烃氢过氧化物 |
| 111 | `IOLENIT` | C4H9NO4 | 内烯烃硝酸酯（NOy） |

未分组: 无（111 全覆盖）

物种总数: 111


### 4.1 家族速查

| 家族 | 成员 |
|---|---|
| **NOy**（25，batch 权重表 + indicators.py 自动判定） | NO(1) NO₂(1) HNO₃(1) HONO(1) NO₃(1) N₂O₅(2) PAN(1) MPAN(1) ISOPNIT(1) ISOPDNT(1) ClONO₂(1) BrONO₂(1) ISOPNO₃(1) PNA(1) RONO₂(1) OLENIT(1) TERPNIT(1) **CSLNIT(1) LALKNIT(1) MTBENIT(1) ESTERNIT(1) SESQNIT(1) ETHNIT(1) IOLENIT(1) PPAN(1)**（扩展 VOC 谱） |
| **HOx** | OH, HO₂（+ H₂O₂ 为储库） |
| **ROx** | HO₂ + RO₂ 池（6 个显式过氧） |
| **卤素** | Cl, ClO, ClONO₂, HCl, HOCl / Br, BrO, BrONO₂, HBr, HOBr |
| **硫** | SO₂, SO₃, H₂SO₄, DMS, MSA |
| **初级排放（只耗不产）** | CH₄, ISOP, OLE, ALK, ARO1, ARO2, KET, SO₂, DMS, TERP |

---

## 5. 反应全表

> 速率列：`名称 = 表达式`（表达式为 `.fac` 源码逐字）；`@n` = 幂运算，
> `M` = 第三体数密度（molec/cm³），`TEMP` = 温度（K），`EXP` = 自然指数。
> 光解行标 `J<slot>`，参数见 §6。分数系数（如 0.85）为化学计量产率。
> 参数来源：`[MCM]` 标注者取自 MCM v3.3.1 逐字拷贝；其余为 JPL/IUPAC 推荐值或文献代理。


#### ① 无机核心：NOx 光稳态 + HOx + HONO 晨爆（26 条）

| # | 速率 | 反应 |
|---|---|---|
| 1 | `J<4>` (光解) | `NO2 = NO + O` |
| 2 | `J<1>` (光解) | `O3 = O1D` |
| 3 | `J<2>` (光解) | `O3 = O` |
| 4 | `J<3>` (光解) | `H2O2 = OH + OH` |
| 5 | `J<8>` (光解) | `HNO3 = NO2 + OH` |
| 6 | `J<7>` (光解) | `HONO = NO + OH` |
| 7 | 6.0D-34*(TEMP/300)@-2.6*M*O2 | `O + O2 + M = O3` |
| 8 | `KO3NO` = 1.4D-12*EXP(-1310/TEMP) | `O3 + NO = NO2` |
| 9 | `KO3NO2NO3` = 1.4D-13*EXP(-2470/TEMP) | `O3 + NO2 = NO3` |
| 10 | `KONO2` = 5.5D-12*EXP(188/TEMP) | `O + NO2 = NO` |
| 11 | `KMT08` = (K80*K8I)*F8/(K80+K8I) | `NO2 + OH = HNO3` |
| 12 | `KMT07` = (K70*K7I)*F7/(K70+K7I) | `NO + OH = HONO` |
| 12b | `KO3HO2` = 2.03D-16*(TEMP/300)@(4.57)*EXP(693/TEMP) | `HO2 + O3 = OH + 2 O2`（**V2.5 补**） |
| 12c | `KO3OH` = 1.70D-12*EXP(-940/TEMP) | `O3 + OH = HO2 + O2`（**V2.5 补**） |
| 13 | `KHONOOH` = 2.5D-12*EXP(260/TEMP) | `HONO + OH = NO2` |
| 14 | `KO1DH2O` = 2.2D-10*H2O | `O1D + H2O = OH + OH` |
| 15 | `KO1DM` = 2.9D-11*M | `O1D + M = O` |
| 16 | `KCOOH` = 1.44D-13*(1+(M/4.2D19)) | `CO + OH = HO2` |
| 17 | `KHCHOOH` = 5.4D-12*EXP(135/TEMP) | `HCHO + OH = HO2 + CO` |
| 18 | `J<11>` (光解) | `HCHO = HO2 + HO2 + CO` |
| 19 | `J<12>` (光解) | `HCHO = CO` |
| 20 | `KHO2NO` = 3.45D-12*EXP(270/TEMP)*EXP(270/TEMP) | `HO2 + NO = NO2 + OH` |
| 21 | `KHO2HO2A` = 2.20D-13*KMT06*EXP(600/TEMP) | `HO2 + HO2 = H2O2` |
| 22 | `KHO2HO2B` = 1.90D-33*M*KMT06*EXP(980/TEMP) | `HO2 + HO2 = H2O2` |
| 23 | `KH2O2OH` = 2.9D-12*EXP(-160/TEMP) | `H2O2 + OH = HO2` |
| 24 | 4.8D-11*EXP(250/TEMP) | `OH + HO2 = O2` |

#### ② NO₃/N₂O₅ 夜化学（MCM Troe 逐字）（10 条）

| # | 速率 | 反应 |
|---|---|---|
| 25 | `J<5>` (光解) | `NO3 = NO` |
| 26 | `J<6>` (光解) | `NO3 = NO2 + O` |
| 27 | `KMT03` = (K30*K3I)*F3/(K30+K3I) | `NO2 + NO3 = N2O5` |
| 28 | `KMT04` = (K40*K4I)*F4/(K40+K4I) | `N2O5 = NO2 + NO3` |
| 29 | 1.8D-11*EXP(110/TEMP) | `NO + NO3 = NO2 + NO2` |
| 30 | 4.50D-14*EXP(-1260/TEMP) | `NO2 + NO3 = NO + NO2` |
| 31 | 2.0D-11 | `OH + NO3 = HO2 + NO2` |
| 32 | 4.0D-12 | `HO2 + NO3 = OH + NO2` |
| 33 | 5.5D-16 | `NO3 + HCHO = HNO3 + CO + HO2` |
| 34 | 4.0D-04 | `N2O5 = HNO3 + HNO3` |

#### ③ CH₄ 氧化链 + PAN/CH₃CO₃ 储库循环（16 条）

| # | 速率 | 反应 |
|---|---|---|
| 35 | `KCH4OH` = 1.85D-12*EXP(-1690/TEMP) | `CH4 + OH = CH3O2` |
| 36 | `KRO2NO` = 2.7D-12*EXP(360/TEMP) | `CH3O2 + NO = CH3O + NO2` |
| 37 | `KRO2HO2` = 2.91D-13*EXP(1300/TEMP) | `CH3O2 + HO2 = CH3OOH` |
| 38 | 1.03D-13*EXP(365/TEMP) | `CH3O2 + CH3O2 = HCHO + HCHO + HO2 + HO2` |
| 39 | `J<41>` (光解) | `CH3OOH = CH3O + OH` |
| 40 | 5.3D-12 | `CH3OOH + OH = HCHO + OH` |
| 41 | 7.3D-13*EXP(580/TEMP)*O2 | `CH3O + O2 = HCHO + HO2` |
| 42 | `KRO2NO3` = 2.3D-12 | `CH3O2 + NO3 = CH3O + NO2` |
| 43 | `J<13>` (光解) | `CH3CHO = CH3O2 + CO + HO2` |
| 44 | 5.6D-12 | `CH3CHO + OH = CH3CO3` |
| 45 | `KNO3AL` = 1.4D-12*EXP(-1860/TEMP) | `CH3CHO + NO3 = HNO3 + CH3CO3` |
| 46 | `KFPAN` = (KC0*KCI)*FC/(KC0+KCI) | `CH3CO3 + NO2 = PAN` |
| 47 | `KBPAN` = (KD0*KDI)*FD/(KD0+KDI) | `PAN = CH3CO3 + NO2` |
| 48 | 2.0D-11 | `CH3CO3 + NO = CH3O2 + CO2 + NO2` |
| 49 | 4.2D-13 | `CH3CO3 + HO2 = CH3O2 + CO2 + OH` |
| 50 | 2.0D-12 | `CH3CO3 + CH3O2 = HCHO + HO2 + CO2 + CH3O2` |

#### ④ ISOP 完整二代氧化链（24 条）

| # | 速率 | 反应 |
|---|---|---|
| 51 | `KISOPOH` = 2.7D-11*EXP(390/TEMP) | `ISOP + OH = ISOPBO2` |
| 52 | `KISOPNO3` = 3.03D-12*EXP(-446/TEMP) | `ISOP + NO3 = ISOPNO3` |
| 53 | `KISOPO3` = 1.28D-17 | `ISOP + O3 = 0.50 MVK + 0.50 MACR + 0.57 HCHO + 0.15 CH3CHO + 0.13 CO2 + 0.26 OH + 0.38 O2` |
| 54 | `KISOPBO2NO` = 4.0D-12 | `ISOPBO2 + NO = 0.90 ISOPBO + 0.90 NO2 + 0.10 ISOPNIT` |
| 55 | `KRO2HO2` = 2.91D-13*EXP(1300/TEMP) | `ISOPBO2 + HO2 = ISOPOOH` |
| 56 | `KRO2NO3` = 2.3D-12 | `ISOPBO2 + NO3 = ISOPBO + NO2` |
| 57 | 7.3D-13*EXP(580/TEMP)*O2 | `ISOPBO + O2 = 0.50 MVK + 0.50 MACR + HCHO + HO2` |
| 58 | `J<41>` (光解) | `ISOPOOH = ISOPBO + OH`（**V2.5 移至 J41 通用 ROOH 槽**，与 CH₃OOH 共用） |
| 59 | 1.9D-11 | `ISOPOOH + OH = ISOPBO2` |
| 60 | 4.0D-12 | `ISOPNO3 + NO = MACR + HCHO + 2 NO2` |
| 61 | `KRO2HO2` = 2.91D-13*EXP(1300/TEMP) | `ISOPNO3 + HO2 = ISOPDNT` |
| 62 | `KMVKOH` = 2.6D-12*EXP(610/TEMP)*EXP(452/TEMP) | `MVK + OH = MVKO2` |
| 63 | `KRO2NO` = 2.7D-12*EXP(360/TEMP) | `MVKO2 + NO = MVKO + NO2` |
| 64 | `KRO2HO2` = 2.91D-13*EXP(1300/TEMP) | `MVKO2 + HO2 = MVKOOH` |
| 65 | 1.9D-11 | `MVKOOH + OH = MVKO2` |
| 66 | 7.3D-13*EXP(580/TEMP)*O2 | `MVKO + O2 = CH3CO3 + HCHO + HCHO` |
| 67 | `KMACROH` = 8.0D-12*EXP(380/TEMP)*EXP(540/TEMP) | `MACR + OH = MACRO2` |
| 68 | `KFPAN` = (KC0*KCI)*FC/(KC0+KCI) | `MACRO2 + NO2 = MPAN` |
| 69 | `KMPANTH` = 1.6D+16*EXP(-13500/TEMP) | `MPAN = MACRO2 + NO2` |
| 70 | 4.0D-12 | `MACRO2 + NO = MACROO + NO2` |
| 71 | 7.3D-13*EXP(580/TEMP)*O2 | `MACROO + O2 = MGLY + CO + OH` |
| 72 | `KRO2HO2` = 2.91D-13*EXP(1300/TEMP) | `MACRO2 + HO2 = MACROOH` |
| 73 | `J<34>` (光解) | `MGLY + 2 O2 = CH3CO3 + CO + HO2`（**V2.5 移至 J34 甲基乙二醛官方槽**，原 J22 慢 31×） |
| 74 | `KMGLYOH` = 1.9D-12*EXP(575/TEMP) | `MGLY + OH = CH3CO3 + CO` |

#### ⑤ 硫链 DMS→SO₂→SO₃→H₂SO₄（4 条）

| # | 速率 | 反应 |
|---|---|---|
| 75 | `KSO2OH` = 9.0D-13 | `SO2 + OH = SO3 + HO2` |
| 76 | `KSO3HYD` = 1.0D-2 | `SO3 = H2SO4` |
| 77 | `KDMSOH` = 1.12D-11*EXP(-250/TEMP) | `DMS + OH = MSA + CH3O2` |
| 78 | `KDMSNO3` = 1.9D-13*EXP(520/TEMP) | `DMS + NO3 = MSA + CO + NO2` |

#### ⑥ 卤素催化 O₃ 损耗 + 储库闭环（19 条）

| # | 速率 | 反应 |
|---|---|---|
| 79 | `KCLO3` = 2.3D-11*EXP(-200/TEMP) | `Cl + O3 = ClO` |
| 80 | `KCLONO` = 6.4D-12*EXP(290/TEMP) | `ClO + NO = Cl + NO2` |
| 81 | `KCLONO2` = 1.6D-12 | `ClO + NO2 = ClONO2` |
| 82 | `J<26>` (光解) | `ClONO2 = Cl + NO3` |
| 83 | `KCLOHO2` = 2.2D-12*EXP(300/TEMP) | `ClO + HO2 = HOCl` |
| 84 | `J<28>` (光解) | `HOCl = Cl + OH` |
| 85 | `KCLCH4` = 9.6D-12*EXP(-1140/TEMP) | `Cl + CH4 = HCl + CH3O2` |
| 86 | `KCLHCHO` = 8.1D-11*EXP(-34/TEMP) | `Cl + HCHO = HCl + HO2 + CO` |
| 87 | `KHCLOH` = 7.8D-13*EXP(-290/TEMP) | `HCl + OH = Cl` |
| 88 | `KBRO3` = 1.6D-12*EXP(-780/TEMP) | `Br + O3 = BrO` |
| 89 | `KBRONO` = 8.7D-12*EXP(260/TEMP) | `BrO + NO = Br + NO2` |
| 90 | `KBRONO2` = 1.5D-12 | `BrO + NO2 = BrONO2` |
| 91 | `J<27>` (光解) | `BrONO2 = Br + NO3` |
| 92 | `KBROHO2` = 2.5D-12*EXP(520/TEMP) | `BrO + HO2 = HOBr` |
| 93 | `J<29>` (光解) | `HOBr = Br + OH` |
| 94 | `KBRHCHO` = 3.4D-12*EXP(-640/TEMP) | `Br + HCHO = HBr + HO2 + CO` |
| 95 | `KHBROH` = 1.1D-11 | `HBr + OH = Br` |
| 96 | `J<45>` (光解) | `ClO = Cl + O`（**V2.5 移至 J45 自定义槽**，JMAP D2.0e-5） |
| 97 | `J<25>` (光解) | `BrO = Br + O` |

#### ⑦ VOC 降解谱集总类（21 条）

| # | 速率 | 反应 |
|---|---|---|
| 98 | `KOLEOH` = 2.6D-11 | `OLE + OH = OLERO2 + O2` |
| 99 | `KOLEO3` = 1.0D-17 | `OLE + O3 = 1.35 HCHO + 0.75 CH3CHO + 0.15 CO2 + 0.30 OH + 0.15 O2` |
| 100 | `KOLENO3` = 9.4D-15 | `OLE + NO3 + H2O = OLERO2 + HNO3` |
| 101 | 0.97×`KRO2NO` = 2.7D-12*EXP(360/TEMP) | `OLERO2 + NO = 0.50 HCHO + CH3CHO + 0.50 CO + HO2 + NO2 + 0.50 H2` |
| 102 | 0.03×`KRO2NO` = 2.7D-12*EXP(360/TEMP) | `OLERO2 + NO = OLENIT` |
| 103 | `KRO2HO2` = 2.91D-13*EXP(1300/TEMP) | `OLERO2 + HO2 = OLEOOH` |
| 104 | `KALKOH` = 1.1D-12 | `ALK + OH = ALKRO2 + H2O` |
| 105 | 0.85×`KRO2NO` = 2.7D-12*EXP(360/TEMP) | `ALKRO2 + NO = CH3CHO + HCHO + HO2 + NO2` |
| 106 | 0.15×`KRO2NO` = 2.7D-12*EXP(360/TEMP) | `ALKRO2 + NO = RONO2` |
| 107 | `KRO2HO2` = 2.91D-13*EXP(1300/TEMP) | `ALKRO2 + HO2 = ALKOOH + O2` |
| 108 | `KHCLOH` = 7.8D-13*EXP(-290/TEMP) | `ALKOOH + OH = ALKRO2 + H2O` |
| 109 | `KARO1OH` = 5.6D-12 | `ARO1 + OH = ARO1O2 + O2` |
| 110 | `KARO2OH` = 2.31D-11 | `ARO2 + OH = ARO2O2 + O2` |
| 111 | `KRO2NO` = 2.7D-12*EXP(360/TEMP) | `ARO1O2 + NO + O2 = GLY + MGLY + CH3CO3 + NO2` |
| 112 | `KRO2NO` = 2.7D-12*EXP(360/TEMP) | `ARO2O2 + NO + 2.50 O2 = 0.50 GLY + MGLY + 2.00 CH3CO3 + NO2` |
| 113 | `KKETOH` = 8.8D-12*EXP(-1320/TEMP)+1.7D-14*EXP(423/TEMP) | `KET + OH = KETRO2 + H2O` |
| 114 | `KRO2NO` = 2.7D-12*EXP(360/TEMP) | `KETRO2 + NO + O2 = CH3CO3 + HCHO + NO2` |
| 115 | `J<21>` (光解) | `KET = CH3CO3 + CH3O2 + 0.50 O2`（**V2.5 移至 J21 丙酮官方槽**） |
| 116 | `J<31>` (光解) | `GLY = CO + CO + H2`（**V2.5 补分支**，占 60.7%） |
| 116b | `J<32>` (光解) | `GLY = HCHO + CO`（原唯一分支，现占 9.1%） |
| 116c | `J<33>` (光解) | `GLY + 2 O2 = CO + CO + HO2 + HO2`（**V2.5 补分支**，占 30.2%，HO₂ 源） |

#### ⑧ 支撑化学（PNA / N₂O₅ 光解 / OH+HNO₃）（3 条）

| # | 速率 | 反应 |
|---|---|---|
| 117 | `KPNAF` = 1.4D-12 | `HO2 + NO2 = PNA` |
| 118 | `KPNAT` = 1.3D+16*EXP(-10400/TEMP) | `PNA = HO2 + NO2` |
| 119 | `KOHPNA` = 4.5D-12 | `PNA + OH = NO2 + O2 + H2O` |

#### ⑨ 烷基硝酸盐家族（3 条）

| # | 速率 | 反应 |
|---|---|---|
| 120 | `KRONO2OH` = 1.2D-12 | `RONO2 + OH = HNO3 + KET + 0.5 H2O + 0.25 O2` |
| 121 | `KOLENITOH` = 1.2D-12 | `OLENIT + OH = HNO3 + OLERO2` |
| 122 | `KRONO2OH` = 1.2D-12 | `ISOPNIT + OH = HNO3 + ISOPBO` |

#### ⑩ 单萜迷你链（9 条）

| # | 速率 | 反应 |
|---|---|---|
| 123 | `KTERPOH` = 1.2D-11*EXP(440/TEMP) | `TERP + OH = TERPRO2 + H2O` |
| 124 | `KTERPO3` = 8.05D-16*EXP(-640/TEMP) | `TERP + O3 = TERPCHO + 0.80 OH + 0.10 O2 + 0.60 H2` |
| 125 | `KTERPNO3` = 1.2D-12*EXP(490/TEMP) | `TERP + NO3 = TERPNIT` |
| 126 | `KTERPOH` = 1.2D-11*EXP(440/TEMP) | `TERPRO2 + NO = TERPCHO + NO2 + 1.50 H2` |
| 127 | `KRO2HO2` = 2.91D-13*EXP(1300/TEMP) | `TERPRO2 + HO2 = TERPOOH` |
| 128 | `KRO2HO2` = 2.91D-13*EXP(1300/TEMP) | `TERPNIT + HO2 = TERPCHO + HNO3 + H2` |
| 129 | `J<30>` (光解) | `N2O5 = NO2 + NO3`（**V2.5 移至 J30 自定义槽**，JMAP D1.0e-5） |
| 130 | `KMT11` = K1H + K2H | `HNO3 + OH = NO3 + H2O` |
| 131 | `KRO2NO3` = 2.3D-12 | `ISOPBO2 + CH3O2 = ISOPBO + HCHO + HO2` |

#### ⑪ 芳香烃分级 + RO₂+HO₂ 闭环（10 条）

| # | 速率 | 反应 |
|---|---|---|
| 132 | `KCSLOH` | `CSL + OH + 1.5 O2 = 0.5 GLY + 2 MGLY` |
| 133 | `KCSLNO3` | `CSL + NO3 + O2 = CSLRO2 + HNO3` |
| 134 | `KRO2NO` | `CSLRO2 + NO + 3 O2 = GLY + MGLY + 2 CO + NO2 + HO2` |
| 135 | `0.15*KRO2NO` | `CSLRO2 + NO = CSLNIT` |
| 136 | `KRO2HO2` | `CSLRO2 + HO2 = CSLOOH` |
| 137 | `KRO2HO2` | `CSLOOH + OH = CSLRO2 + H2O + O2` |
| 138 | `KRO2HO2` | `ARO1O2 + HO2 = TOLOOH` |
| 139 | `KRO2HO2` | `ARO2O2 + HO2 = XYLOOH` |
| 140 | `KRO2HO2` | `TOLOOH + OH = ARO1O2 + H2O + O2` |
| 141 | `KRO2HO2` | `XYLOOH + OH = ARO2O2 + H2O + O2` |

#### ⑫ 长链烷烃 LALK（6 条）

| # | 速率 | 反应 |
|---|---|---|
| 142 | `KLALKOH` | `LALK + OH + O2 = LALKRO2 + H2O` |
| 143 | `0.85*KRO2NO` | `LALKRO2 + NO + 2 O2 = 3 CH3CHO + NO2 + HO2` |
| 144 | `0.15*KRO2NO` | `LALKRO2 + NO = LALKNIT` |
| 145 | `KRO2HO2` | `LALKRO2 + HO2 = LALKOOH` |
| 146 | `KRO2HO2` | `LALKOOH + OH = LALKRO2 + H2O + O2` |
| 147 | `KRONO2OH` | `LALKNIT + OH + 1.25 O2 = HNO3 + 3 CH3CHO + 0.5 H2O` |

#### ⑬ 醇/醚/酯类含氧 VOC（14 条）

| # | 速率 | 反应 |
|---|---|---|
| 148 | `KMEOHOH` | `MEOH + OH + O2 = HO2 + HCHO + H2O` |
| 149 | `KETOHOH` | `ETOH + OH + O2 = CH3CHO + HO2 + H2O` |
| 150 | `KMTBEOH` | `MTBE + OH + O2 = MTBERO2 + H2O` |
| 151 | `0.85*KRO2NO` | `MTBERO2 + NO + 1.5 O2 = 2 CH3CHO + HCHO + NO2 + HO2` |
| 152 | `0.15*KRO2NO` | `MTBERO2 + NO = MTBENIT` |
| 153 | `KRO2HO2` | `MTBERO2 + HO2 = MTBEOOH` |
| 154 | `KRO2HO2` | `MTBEOOH + OH = MTBERO2 + H2O + O2` |
| 155 | `KRONO2OH` | `MTBENIT + OH + 0.75 O2 = HNO3 + 2 CH3CHO + HCHO + 0.5 H2O` |
| 156 | `KESTEROH` | `ESTER + OH + 1.5 O2 = CH3CHO + CO2 + CH3O2 + H2O` |
| 157 | `0.85*KRO2NO` | `ESTERO2 + NO + O2 = CH3CHO + CO2 + CH3O2 + NO2` |
| 158 | `0.15*KRO2NO` | `ESTERO2 + NO = ESTERNIT` |
| 159 | `KRO2HO2` | `ESTERO2 + HO2 = ESTEROOH` |
| 160 | `KRO2HO2` | `ESTEROOH + OH = ESTERO2 + H2O + O2` |
| 161 | `KRONO2OH` | `ESTERNIT + OH + O2 = HNO3 + CH3CHO + CO2 + CH3O2` |

#### ⑭ 倍半萜 SESQ（8 条）

| # | 速率 | 反应 |
|---|---|---|
| 162 | `KSESQOH` | `SESQ + OH + 1.5 O2 = SESQRO2 + H2O` |
| 163 | `KSESQO3` | `SESQ + O3 = SESQCHO + 0.80 OH + 0.10 O2 + 0.60 H2` |
| 164 | `KSESQNO3` | `SESQ + NO3 + 0.75 O2 = SESQNIT + 0.5 H2O` |
| 165 | `0.85*KRO2NO` | `SESQRO2 + NO + H2O + O2 = SESQCHO + NO2 + HO2 + H2O` |
| 166 | `0.15*KRO2NO` | `SESQRO2 + NO = SESQNIT` |
| 167 | `KRO2HO2` | `SESQRO2 + HO2 = SESQOOH` |
| 168 | `KRO2HO2` | `SESQOOH + OH = SESQRO2 + H2O + O2` |
| 169 | `KRONO2OH` | `SESQNIT + OH + 0.25 O2 = HNO3 + SESQCHO + 0.5 H2O` |

#### ⑮ 乙烯/内烯烃/高级醛（24 条）

| # | 速率 | 反应 |
|---|---|---|
| 170 | `KETHOH` | `ETH + OH + O2 = ETHO2` |
| 171 | `KETHO3` | `ETH + O3 + 0.2 H2O = 2 HCHO + 0.4 OH + 0.4 O2` |
| 172 | `0.85*KRO2NO` | `ETHO2 + NO + O2 = 2 HCHO + NO2 + HO2` |
| 173 | `0.15*KRO2NO` | `ETHO2 + NO = ETHNIT` |
| 174 | `KRO2HO2` | `ETHO2 + HO2 = ETHOOH` |
| 175 | `KRO2HO2` | `ETHOOH + OH = ETHO2 + H2O + O2` |
| 176 | `KRONO2OH` | `ETHNIT + OH + 0.25 O2 = HNO3 + 2 HCHO + 0.5 H2O` |
| 177 | `KALD2OH` | `ALD2 + OH + O2 = ALD2O2 + H2O` |
| 178 | `KNO3AL` | `ALD2 + NO3 + 2.5 O2 = HNO3 + CH3CO3 + CO2 + H2O` |
| 179 | `J<14>` | `ALD2 + 0.5 O2 = CH3CHO + CO + H2`（**V2.5 移至 J14 丙醛官方槽**，J13 归 CH₃CHO） |
| 180 | `KRO2NO` | `ALD2O2 + NO + 2 O2 = CH3CO3 + CO2 + NO2 + H2O` |
| 181 | `KFPAN` | `ALD2O2 + NO2 = PPAN` |
| 182 | `KBPAN` | `PPAN = ALD2O2 + NO2` |
| 183 | `KRO2HO2` | `ALD2O2 + HO2 = ALD2OOH` |
| 184 | `KRO2HO2` | `ALD2OOH + OH = ALD2O2 + H2O + O2` |
| 185 | `KIOLEOH` | `IOLE + OH + O2 = IOLERO2` |
| 186 | `KIOLEO3` | `IOLE + O3 = 2 HCHO + CH3CHO` |
| 187 | `KIOLENO3` = 3.5D-13 | `IOLE + NO3 + 0.5 H2O + 0.25 O2 = IOLENIT`（**V2.5: 内烯专属**，原挪用丙烯 9.4e-15 偏慢 37×） |
| 187b | `3.3D-12*EXP(-2880/TEMP)`（MCM 逐字） | `ETH + NO3 + 0.5 H2O + 0.25 O2 = ETHNIT`（**V2.5+ 补**，k=2.11e-16；ETH 原无直接 NO₃ 汇） |
| 188 | `0.85*KRO2NO` | `IOLERO2 + NO + O2 = 2 CH3CHO + HO2 + NO2` |
| 189 | `0.15*KRO2NO` | `IOLERO2 + NO = IOLENIT` |
| 190 | `KRO2HO2` | `IOLERO2 + HO2 = IOLEOOH` |
| 191 | `KRO2HO2` | `IOLEOOH + OH = IOLERO2 + H2O + O2` |
| 192 | `KOLENITOH` | `IOLENIT + OH + 0.25 O2 = HNO3 + 2 CH3CHO + 0.5 H2O` |

#### ⑯ MACR 夜化学（1 条）

| # | 速率 | 反应 |
|---|---|---|
| 193 | `KNO3AL` | `MACR + NO3 + O2 = HNO3 + MACRO2` |

#### ⑰ RO₂ 交叉缔合（5 条）

| # | 速率 | 反应 |
|---|---|---|
| 194 | `KRO2NO3` | `CH3CO3 + ISOPBO2 = CH3O2 + ISOPBO + CO2` |
| 195 | `KRO2NO3` | `CH3O2 + MACRO2 = CH3O + MACROO + O2` |
| 196 | `KRO2NO3` | `CH3O2 + MVKO2 = CH3O + MVKO + O2` |
| 197 | `KRO2NO3` | `ISOPBO2 + MACRO2 = ISOPBO + MACROO + O2` |
| 198 | `KRO2NO3` | `ETHO2 + CH3O2 = CH3O + 2 HCHO + HO2` |

> **全部 203 条反应精确平衡**（`audit` 无储库依赖，实测确认）。
> 反应行内的 `O2` 出现在反应物侧是物理上必需的（RO₂ 形成需 O₂ 加成），
> 不是配平凑数——`element_audit` 的 O 记账会精确核到每个 O 原子。


### 5.1 关键反应解读（教学要点）

**两条 ★ O₃ 净生成通道**（机制的心脏）：
```
HO₂ + NO → NO₂ + OH      k = 8.6e-12·exp(270/T) ≈ 1.7e-11 @298K
RO₂ + NO → RO + NO₂      k = 2.7e-12·exp(360/T) ≈ 8.7e-12 @298K（集总 RO₂）
```
二者都是"把 NO 氧化成 NO₂ 但不消耗 O₃"——随后 NO₂+hν→NO+O、O+O₂+M→O₃
闭环输出一个 O₃。**P(O₃) = k(HO₂+NO)[HO₂][NO] + k(RO₂+NO)[RO₂][NO]**。

**臭氧化 OH 产率**（夜间/黎明自由基源）：
```
OLE + O₃ → 1.35 HCHO + 0.75 CH3CHO + 0.15 CO2 + 0.30 OH + 0.15 O2   (丙烯, OH 产率 0.30)
ISOP + O₃ → 0.5 MVK + 0.5 MACR + 0.57 HCHO + 0.15 CH3CHO + 0.13 CO2 + 0.26 OH + 0.38 O2   (OH 产率 0.26)
TERP + O₃ → TERPCHO + 0.80 OH + 0.10 O2 + 0.60 H2   (单萜, OH 产率 0.80 — 文献特征值)
```
真实过程经 Criegee 中间体（本机制不显式建模，用产率法等效）：无光解时段
这是 OH 的唯一来源之一，缺了它夜间 OH 会系统性偏低。

**烷基硝酸酯分支**（NOx 远程输送）：
```
0.85×KRO2NO : ALKRO2 + NO → CH3CHO + HCHO + HO2 + NO2     (主通道)
0.15×KRO2NO : ALKRO2 + NO → RONO2                          (硝酸酯, 寿命长 → 可输送)
0.97×KRO2NO : OLERO2 + NO → 0.5HCHO + CH3CHO + 0.5CO + HO2 + NO2 + 0.5H2
0.03×KRO2NO : OLERO2 + NO → OLENIT
0.90/0.10   : ISOPBO2 + NO → ISOPBO + NO2 / ISOPNIT
```
分支必须**速率缩放**（`0.85*K` 主通道 + `0.15*K` 硝酸酯支），不能产物系数缩放——
否则两通道速率和超出总 k（历史首版的 4 处审计 FAIL 即此）。

**PAN 储库的温度开关**：
```
生成 KFPAN = (KC0·KCI)·FC/(KC0+KCI)          Troe, MCM 逐字
热解 KBPAN = (KD0·KDI)·FD/(KD0+KDI)          Ea≈84 kJ/mol 表观
```
低温稳定（对流层上部运输）、高温快解（暖海面/夏季午后释放 NOx）——Ea 是全机制最大之一。

**卤素闭环**（每个 X 双储库）：
```
Cl + O₃ → ClO ;  ClO + NO → Cl + NO₂          催化循环
ClO + NO₂ → ClONO₂ (储库) ;  ClONO₂ + hν → Cl + NO₃   再生
ClO + HO₂ → HOCl (储库) ;  HOCl + hν → Cl + OH        再生
Cl + CH₄/HCHO → HCl + … ;  HCl + OH → Cl              再生
```
净效果 `O₃ + NO → NO₂ + O₂`；一个 Cl 可催化销毁 ~10⁵ 个 O₃。

---

## 6. 光解方案

**机制共 26 条光解反应、占用 25 个不同槽位（最高槽号 J45）**，参数来自
`<程序>/config_override/photolysis-rates_MCM`（MCM 官方表，18 槽）或 `photolysis.config` 的 JMAP 直达值（7 槽）。

> ⚠ **`.fac` 头部第 15 行的注释**（"JMAP 仅用于表内缺失的卤素槽 25-29"）与实际不符：
> JMAP 实为 **7 槽 = 25/26/27/28/29/30/45**（含非卤素的 N₂O₅@30 与 ClO@45）。
> 该行属机制文件内的注释，**不改**——`.fac` 的 md5 被写进每个程序的 `run_meta.txt`
> 作为跨轮次同一性凭据，动一个注释字符就会让现存存档的 md5 全部对不上。
> 本表与 `jslot_audit.py` 的不变量④才是权威口径（它逐行核对 README 与 `.fac` 的光解表）。

| 槽位 | 反应 | l [s⁻¹] | m | n | 参数来源 |
|---|---|---|---|---|---|
| 1 | `O3 = O1D` | 6.073e-05 | 1.743 | 0.474 | MCM 表 |
| 2 | `O3 = O` | 0.0004775 | 0.298 | 0.08 | MCM 表 |
| 3 | `H2O2 = OH + OH` | 1.041e-05 | 0.723 | 0.279 | MCM 表 |
| 4 | `NO2 = NO + O` | 0.01165 | 0.244 | 0.267 | MCM 表 |
| 5 | `NO3 = NO` | 0.02485 | 0.168 | 0.108 | MCM 表 |
| 6 | `NO3 = NO2 + O` | 0.1747 | 0.155 | 0.125 | MCM 表 |
| 7 | `HONO = NO + OH` | 0.002644 | 0.261 | 0.288 | MCM 表 |
| 8 | `HNO3 = NO2 + OH` | 9.312e-07 | 1.23 | 0.307 | MCM 表 |
| 11 | `HCHO = HO2 + HO2 + CO` | 4.642e-05 | 0.762 | 0.353 | MCM 表 |
| 12 | `HCHO = CO` | 6.853e-05 | 0.477 | 0.323 | MCM 表 |
| 13 | `CH3CHO = CH3O2 + CO + HO2` | 7.344e-06 | 1.202 | 0.417 | MCM 表值（**V2.5 恢复官方表值**，原 JMAP 4.0e-5 偏快 4.4×） |
| 14 | `ALD2 + 0.5 O2 = CH3CHO + CO + H2` | 2.879e-05 | 1.067 | 0.358 | MCM 表值（**V2.5 移至 J14 丙醛官方槽**；原此槽被 N₂O₅ 占用） |
| 21 | `KET = CH3CO3 + CH3O2 + 0.50 O2` | 7.992e-07 | 1.578 | 0.271 | MCM 表值（**V2.5 移至 J21 丙酮官方槽**，原 J15 快 33×） |
| 30 | `N2O5 = NO2 + NO3` | 1e-05 | 0.79 | 0.33 | JMAP 直达 D1.0e-5（**V2.5 自 J14 迁出**；MCM 无 N₂O₅ 光解槽） |
| 31 | `GLY = CO + CO + H2` | 6.845e-05 | 0.13 | 0.201 | MCM 表值（**V2.5 补分支**，占 60.7%） |
| 32 | `GLY = HCHO + CO` | 1.032e-05 | 0.13 | 0.201 | MCM 表值（**V2.5 移至 J32 官方槽**） |
| 33 | `GLY + 2 O2 = CO + CO + HO2 + HO2` | 3.802e-05 | 0.644 | 0.312 | MCM 表值（**V2.5 补分支**，占 30.2%，HO₂ 源） |
| 34 | `MGLY + 2 O2 = CH3CO3 + CO + HO2` | 1.537e-04 | 0.17 | 0.208 | MCM 表值（**V2.5 移至 J34 官方槽**，原 J22 慢 31×） |
| 41 | `CH3OOH = CH3O + OH` | 7.649e-06 | 0.682 | 0.279 | MCM 表（通用 ROOH 槽） |
| 41 | `ISOPOOH = ISOPBO + OH` | 7.649e-06 | 0.682 | 0.279 | MCM 表（**V2.5 移至 J41**，与 CH₃OOH 共用，原 J21 慢 12×） |
| 25 | `BrO = Br + O` | 1e-06 | 0.79 | 0.33 | JMAP 直达 D1.0e-6 (m/n 取有机典型 0.79/0.33) |
| 26 | `ClONO2 = Cl + NO3` | 0.0001 | 0.79 | 0.33 | JMAP 直达 D1.0e-4 (m/n 取有机典型 0.79/0.33) |
| 27 | `BrONO2 = Br + NO3` | 0.001 | 0.79 | 0.33 | JMAP 直达 D1.0e-3 (m/n 取有机典型 0.79/0.33) |
| 28 | `HOCl = Cl + OH` | 0.0002 | 0.79 | 0.33 | JMAP 直达 D2.0e-4 (m/n 取有机典型 0.79/0.33) |
| 29 | `HOBr = Br + OH` | 0.0005 | 0.79 | 0.33 | JMAP 直达 D5.0e-4 (m/n 取有机典型 0.79/0.33) |
| 45 | `ClO = Cl + O` | 2e-05 | 0.79 | 0.33 | JMAP 直达 D2.0e-5（**V2.5 自 J24 迁出**；MCM 无 ClO 物种） |


**参数化式**（`eval_J`，源码）：
```
decl  = 23.45° · sin(2π(DOY−81)/365)                    赤纬
cosχ  = sin(LAT)sin(decl) + cos(LAT)cos(decl)cos(2π(h−12)/24)   天顶角余弦
J     = JFAC · l · cosχ^m · exp(−n/cosχ)   [s⁻¹]         (cosχ ≤ 1e-3 → J = 0)
```

**槽位语义**（MCM 官方编号；**V2.5 起逐槽一致**——机制引用的每个 MCM 编号槽位，其分子
归属与官方导出一致，含代理物种映射）：
`1=O₃→O¹D, 2=O₃→O(³P), 3=H₂O₂, 4=NO₂, 5=NO₃→NO, 6=NO₃→NO₂+O, 7=HONO, 8=HNO₃,
11=HCHO(radical), 12=HCHO(molecular), 13=CH₃CHO, 14=C₂H₅CHO(=ALD2),
21=CH₃COCH₃(=KET 丙酮), 31/32/33=GLYOX 三分支(=GLY), 34=MGLYOX(=MGLY),
41=通用 ROOH（CH₃OOH 与 ISOPOOH 共用，MCM 惯例）`

**自定义槽位**（MCM 无对应物种 → 文献正午值，走 `photolysis.config` JMAP）：
`30=N₂O₅(D1.0e-5), 45=ClO(D2.0e-5), 25=BrO, 26=ClONO₂, 27=BrONO₂, 28=HOCl, 29=HOBr`

> ⚠ **新增光解反应前必须查 `mechanism/reference/` 官方导出中该槽位的分子归属**——
> 槽号是 MCM 官方命名空间里的"外键"，借用错槽 = 借用了别的分子的 l/m/n
> （V2.5 曾因此迁移 7 处）。该维度已由 `mechanism/audit/jslot_audit.py` 自动断言
> （`make jcheck`）：机制每个 `J<n>` 的反应物经集总映射后必须
> ⊆ MCM 官方该槽分子集合，MCM 无基准的槽必须落在脚本的 `CUSTOM` 白名单且配 JMAP
> `D<正午值>`。注意 `mcm_compare.py` 是热反应专用，对错位无感——别指望它兜底。

**新增光解反应的规则**：
1. `.fac` 加 `% J<新槽位> : 反应物 = 产物 ;`
2. `photolysis.config` 的 JMAP 加 `新槽位=D<值>`（或 MCM 表加行并恒等映射）
3. MCM 无对应分子的槽 → 同时加进 `jslot_audit.py` 的 `CUSTOM` 白名单（未声明即 FAIL）
4. `cfgbuild` 的光解完备性校验会自动检查"每个槽位有非零 l"，否则**拒绝编译**（防静默零光解）
5. `make jcheck` 复核槽位语义 + JMAP 死参数 + 本文档 §5/§6 两张表的一致性

---

## 7. 动力学约定

### 7.1 速率常数形式

| 形式 | 写法 | 例 |
|---|---|---|
| Arrhenius | `A*EXP(-Ea/TEMP)` | `KHO2NO = 8.6D-12*EXP(270/TEMP)` |
| 幂律修正 | `A*(TEMP/300)@n` | `KCLO3 = 2.3D-11*EXP(-200/TEMP)` |
| 三体（低压极限） | `k0*M*(TEMP/300)@n` | `K30 = 3.6D-30*M*(TEMP/300)@-4.1` |
| Troe | `(k0·k∞)·F/(k0+k∞)` | 见 7.2 |
| 光解 | `J<n>` | 见 §6 |
| 常数 | 直接值 | `KRO2HO2 = 2.91D-13*EXP(1300/TEMP)` |

### 7.2 Troe 参数化（MCM 逐字，以 N₂O₅ 为例）

```
K30 = 3.6D-30*M*(TEMP/300)@-4.1        低压极限 k0
K3I = 1.9D-12*(TEMP/300)@0.2           高压极限 k∞
KR3 = K30/K3I                           约化压力
FC3 = 0.35                              展宽因子
NC3 = 0.75-1.27*LOG10(FC3)              Troe 指数
F3  = 10@(LOG10(FC3)/(1+(LOG10(KR3)/NC3)@2))   弯折函数
KMT03 = (K30*K3I)*F3/(K30+K3I)          最终速率
```

**机制内 Troe 反应**：KMT03/04（NO₂+NO₃⇌N₂O₅）、KMT07/08（OH+NO→HONO、OH+NO₂→HNO₃）、
KMT11（OH+HNO₃，K1+K2 双通道）、KFPAN/KBPAN（PAN 生成/热解）。

### 7.3 环境量约定

| 量 | 定义 | 位置 |
|---|---|---|
| `M`（第三体） | 6.022e23·P[Pa]/(8.314·T)·1e-6 [molec/cm³] | `eval_env` |
| `O2` | 0.2095·M | 同上 |
| `H2O` | RH/100·es(Magnus)·(数密度换算) | 同上 |
| `TEMP` | meteo 的 TEMP + 情景 DT | `eval_met` |

⚠ **注意**：`O1D+H2O`、`O1D+M`、`RO+O2` 等在 `.fac` 中**已把 H2O/M/O2 显式乘进速率**
（如 `KO1DH2O = 2.2D-10*H2O`），不要在机制层再乘一次。

### 7.4 代码内部再乘反应物浓度

`mechc.py` 生成的 `update_rates` 只算"速率常数 × 环境项"；**反应物浓度由驱动侧
再乘**（`peff(k) = p(k)·∏y[reac_id]`，源码 `rates_env`）。因此 `.fac` 的速率表达式
里**不应再写反应物物种名**（除设计上作为环境项处理的 H2O/O2/M）。

---

## 8. 元素配平方法论

### 8.1 审计器设计（`mechanism/audit/element_audit.py`）

- 逐反应检验 **C/H/N/O/S/Cl/Br** 七元素平衡
- 储库物种（H₂O / O₂ / CO₂ / H₂ / M / N₂）**自动推断**两侧需求，但**限定物理方向**：
  - CO₂、H₂ **只许出现在产物侧**（不允许"用 CO₂ 当反应物配平"）
  - O₂ 反应物侧 ≤ 2.5 分子（防止用过量 O₂ 掩盖真实错误）
- 需求出现负值或非整数/半整数 → FAIL
- 退出码 0=全过 / 1=有失败

### 8.2 为什么方向限定这么重要

若允许任意补充储库，任何反应都能"配平"——审计失去意义。方向限定把审计从
"元素算术"升级为"物理可信性检查"。收紧后旧 expanded 机制立刻现出 **7 处 FAIL**
（见 §2 历史修正表）。

### 8.3 配平实战案例（方法示例）

**案例 1：OLE+O₃ 加 OH 产率**
```
初版：OLE + O3 = 2.00 HCHO + 0.50 O2 + 0.30 OH        → 净差 O:−0.08
二版：OLE + O3 = 1.35 HCHO + 0.75 CH3CHO + 0.145 CO2 + 0.30 OH + 0.34 O2
                                                      → 净差 C:+0.005 O:−0.37
终版：OLE + O3 = 1.35 HCHO + 0.75 CH3CHO + 0.15 CO2 + 0.30 OH + 0.15 O2   ✓
```

**案例 2：ISOP+O₃**
```
初版 O2 0.42 → 净差 O:−0.08 → 0.38 ✓
终版：ISOP + O3 = 0.50 MVK + 0.50 MACR + 0.57 HCHO + 0.15 CH3CHO
                  + 0.13 CO2 + 0.26 OH + 0.38 O2
```

**案例 3：TERP+O₃ 的 O/C 不变式困境**（最具方法论价值）
```
约束：C10 有机物 + O₃ 只有 3 个 O → 要容纳 0.80 OH，任何 O/C≥1 的碳产物都不可行
解法：TERPCHO 定为 O/C=0.2 的汇聚型载体（C10H14O2，吸收全部 10 个 C），
      O₃ 的 3 个 O 全交给 0.80 OH + 0.10 O2；H 由 0.60 H2 汇出
终版：TERP + O3 = TERPCHO + 0.80 OH + 0.10 O2 + 0.60 H2   ✓
教训：真实机制靠 O-free 碳键载体（CB6 的 PAR 键）平衡；本机制无此物种，
      TERPCHO 的 O/C 设计即为此让路，因此它只能是 SOA 代理终产物（走沉降汇）。
```

**案例 4：NO₃ 加成不转移 H**
```
TERPNIT = C10H16NO3（加成自由基保持烯烃氢数 16，不是 17）
H17 版审计立即报 H:−1。同理 NO₃ 加成产物一律不加 H。
```

**案例 5：分支拆分的正确写法**
```
错误：主通道产物乘 0.85（OLE + O3 = 0.85 HCHO + …）
       → 两通道总速率仍为 2×k，超出物理 k
正确：速率缩放  0.85*KRO2NO : ALKRO2+NO = 全量产物
                0.15*KRO2NO : ALKRO2+NO = RONO2
```

**案例 6：行内 `!` 注释会破坏审计**
`.fac` 反应行尾的 `!` 注释不会被解析器剥离 → 被当作物种名 → 审计 FAIL。
注释请放在速率定义行或独立的 `*` 行。

---

## 9. 已知省略（工程决策）

> 以下是 V2.5+ 冻结版**有意保留**的简化项；它们不是遗漏，而是与当前 0-D 定位绑定的取舍。

| 省略项 | 理由 | 影响 | 补齐成本 |
|---|---|---|---|
| **Criegee 显式物种** | 用产率法等效（OH 0.26/0.30/0.80 已显式） | 夜间 OH 的 100% 依赖此产率法；中间体谱缺失 | 高（+4sp/+8rx） |
| **烷基硝酸盐水解** | 无气溶胶耦合 | RONO₂→HNO₃ 非均相通道缺失 | 中 |
| **甲酸/乙酸等有机酸** | 城市贡献小 | 未覆盖 | 低-中 |
| **多环芳烃 (PAH)** | 需专门 SOA 模块 | 未覆盖 | 高 |
| **异戊二烯 4 个以上分支** | 现有三分支为教学精简 | 高 NOx 下硝酸酯分配略偏 | 中 |
| **光解槽位 13–29 的 m/n** | 文献典型形状替代 MCM 官方值 | 光解日变化形状近似 | 低（换 JMAP 值） |
| **RO₂ 交叉缔合的温度依赖** | 统一用 `KRO2NO3 = 2.3e-12` | 低温下略偏 | 低 |

---

## 10. 对比定位

| 维度 | MCM v3.3.1 | GEOS-Chem fullchem | CB6R3 (CMAQ) | RADM2 | MOZART-4 | **本机制 V2.5+** |
|---|---|---|---|---|---|---|
| 物种 / 反应 | ~6700 / ~17000 | ~280 / ~850 | 82 / 221 | 63 / 157 | 85 / 157+39 | **111 / 203** |
| 相对 MCM | 100% | 4.2% / 5.0% | 1.2% / 1.3% | 0.9% / 0.9% | 1.3% / 1.2% | **1.7% / 1.2%** |
| VOC 表示 | 全显式 | 混合（显式芳烃四兄弟 + 集总烷烃） | 碳键官能团 | 代表分子 | 大类集总 | 显式裁剪 + 集总类 |
| 芳香烃 | 全显式 | 显式 | TOL/XYL 分级 | 代表分子 | 大类 | **ARO1/ARO2 + CSL 开环链** |
| 烷烃 | 全显式 | 集总 | 碳键 | 代表分子 | 大类 | **ALK(C3) + LALK(C6)** |
| 烯烃 | 全显式 | 集总 | 碳键 | 代表分子 | 大类 | **OLE + ETH + IOLE** |
| 含氧 VOC | 全显式 | 混合 | ✓ | 部分 | 部分 | **MEOH/ETOH/MTBE/ESTER/ALD2** |
| 异戊二烯 | 全显式 | MARV 完整 | 集总 ISPD | 单物种 | 单物种 | **显式二代链（7 物种）** |
| 单萜 / 倍半萜 | 全显式 | 多种 | TERP | ✗ | ✗ | **TERP(5) + SESQ(5)** |
| 卤素 | 基础 CL | Cl/Br/I 完整 | r3m 变体 | ✗ | ✗（TS2 有） | **闭环 + 4 储库** |
| 夜间化学（NO₃/N₂O₅） | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ Troe 逐字 |
| RO₂ 交叉缔合 | 全谱 | ✓ | ✓ | 部分 | 部分 | **5 对** |
| HONO 晨爆 | ✓ | ✓ | 部分 | ✗ | ✗ | ✓ |
| 气溶胶耦合 | ✗ | ✓（SOA/无机） | ✓ | 部分 | ✓ | ✗（TERPCHO/SESQCHO 代理） |
| 元素审计 | 人工 | 人工 | 预平衡 | 人工 | 人工 | **自动 203/203** |
| 城市 OH 反应性覆盖 | 100% | ~95% | ~95% | ~85% | ~85% | **~95%** |
| 定位 | 全谱参照 | 全球定量 | 区域达标 | 历史 | 全球气候 | **过程诊断/教学/情景** |

**复杂度阶梯（本机制所在位置）**：
```
RADM2 (63sp) ── CB6R3 (82sp) ── MOZART-4 (85sp) ── 本机制 (111sp) ── fullchem (~280sp) ── MCM (~6700sp)
```
即：**已超越 CB6R3/MOZART-4 的物种数**（111 vs 82/85），反应数也已接近
（203 vs 221/196）。但要注意这是**数量上的接近，不是结构上的等价**：
- CB6R3 用**碳键官能团**覆盖全谱 VOC（82 物种浓缩了数千个 MCM 物种）；
  本机制用**显式裁剪**，111 个物种里每一个都是具体分子/集总类。
- 因此本机制在**已覆盖的路径上比 CB6 更细**（如 ISOP 二代三分支、TERP/SESQ 三入口），
  但在**未覆盖的 VOC 谱上仍有空白**（甲酸/乙酸、PAH、更多含氧中间体）。
- 对**过程诊断**（IPR/IRR 逐反应归因）本机制更合适；对**区域达标模拟**
  CB6R3 的集总结构在 3D 网格上更经济。

---

## 11. 单位约定

| 量 | 单位 | 说明 |
|---|---|---|
| 速率常数（双分子） | **molec/cm³/s** | Arrhenius `A·exp(−Ea/T)` |
| 速率常数（单分子/光解） | **s⁻¹** | |
| 三体 M / O₂ / H₂O | **molec/cm³** | 运行时由 PRES/TEMP/RH 计算 |
| 化学计量系数 | 无量纲（**允许分数** 0.03–2.5） | 审计强制配平 |
| 光解 l | **s⁻¹**（天顶角 0） | m/n 无量纲 |
| 浓度（机制内部） | **molec/cm³** | IC/BC 可用 ppb/ppm，cfgbuild 换算 |

⚠ **ppb 换算因子固定 1 atm / 298.15 K（2.46e10）**——非标准 T/P 下模型内实际
混合比有 ~1.1%/5 K 偏差（MCM/AtChem2 同款约定，化学一致性不受影响，仅"ppb 标签"为近似）。
完整单位与格式硬约束见 [`../ENVIRONMENT.md`](../ENVIRONMENT.md) §7。

---

## 12. 使用与扩展指南

### 12.1 修改机制的完整流程（五步闭环）

```bash
# ① 改 .fac（物种/反应/速率）
vim mechanism/mechanism_3_full_chain.fac

# ② 元素审计（必须 PASS 才能继续）
python3 mechanism/audit/element_audit.py mechanism/mechanism_3_full_chain.fac
#    → 期望 "全部 N 反应元素守恒 PASS"

# ③ 若改了物种：补配置（IC/排放/沉降/组映射）
#    - <程序>/config_override/initialConcentrations.config（可选）
#    - <程序>/config_override/emissions.config（若该物种有排放）
#    - <程序>/config_override/drydep.config（若需要沉降汇，尤其：无气相汇的产物）
#    - modules/runtime_dat_mod.f90 的 build_groups（若参与情景缩放）
#    - 若带 N：facsimile_impl/fac_batch_omp.f90 的 noy_add（氮守恒不变量会抓）

# ④ 重编译机制 + 刷新配置
make all                     # mechc 重生成 gen_v3 + gfortran 两个二进制
                             # runtime 由各程序 run.py 每轮自动 cfgbuild 刷新

# ⑤ 回归（当前无数值哨兵：tests/ 已退役，靠人工比对）
make mcmcheck && make jcheck # 两道静态审计
python3 mechanism/audit/element_audit.py mechanism/mechanism_3_full_chain.fac
#    机制改动必然移动输出量级 → 与 ENVIRONMENT.md §8 参考基线逐项人工比对
```

### 12.2 新增一条反应的最小清单

| 步骤 | 检查 |
|---|---|
| 在正确的段（①-⑩）加 `% <速率> : 反应物 = 产物 ;` | 段内顺序即 Fortran 反应号 |
| 速率用命名速率（`KNAME = ...` 在速率区） | 不要写行内 `!` 注释 |
| 元素审计 PASS | 分数系数须精确到 1/4 精度（0.145 会差 0.005） |
| `make all` 重生成 gen_v3 | 生成码里的 `nr/ns` 必须与新 `.fac` 一致（无自动断言，人工 `grep ":: ns"` 核对） |
| 跑一轮 CASE-analysis | 输出量级与 ENVIRONMENT.md §8 参考基线比对（**机制改动必然移动基线**） |
| 更新本文档 §3/§5 | 文档-代码一致性 |

### 12.3 常见陷阱（全部实踩）

| 陷阱 | 症状 | 解法 |
|---|---|---|
| 反应行尾写 `!` 注释 | 审计报"未知物种" | 注释移到独立 `*` 行 |
| 分支用产物系数缩放 | 两通道总速率超物理 k | 用 `0.85*K` 速率缩放 |
| 储库产物无沉降汇 | 5 天无界累积（OLEOOH 曾 +66%） | drydep.config 加行 + 刷新 runtime |
| 新排放物种未进 build_groups | 情景缩放对它无效 | build_groups 加一行；**它是唯一组映射源，禁止在 scen_build/scen5 等处复制第二份清单**（R10-P0 教训：复制导致 6 物种漏缩） |
| 新含 N 物种未进 noy_add | 批量路线**直接 `stop 1`**（氮守恒不变量已升级为致命门，逐案例首快照检），单箱路线 NOy 指标偏低 | 加一行（权重按 N 原子数） |
| 改了 .fac 未重编 | 跑的是旧生成码 | `make all`（生成码与二进制一起重建） |
| 改了 drydep/排放未生效 | 汇/排放不生效 | 重跑 `run.py`（每轮自动 cfgbuild） |
| 反应物物种名写进速率表达式 | 浓度被乘两次 | 只写 H2O/O2/M 等环境项 |
| **配置项存在但无消费者** | 写了却静默无效果（R16: MODE=DIRECT 曾是死开关、RADIAT 无消费者） | 每个 runtime.dat 字段必须有读取点；否则实现它/删除它/显式警告，**不允许静默存在** |
| **乘子型参数在两条路径重复施加** | 默认值=1 时完全无害, 非 1 时静默偏离 (R11: batch jfac 被 `eval_J_jf` 与 `eval_J` 各乘一次, 差 8%) | 用**等价路径断言**守着: 两种设法的结果必须相同 |

### 12.4 快速自检命令

```bash
# 物种/反应计数
grep -c "^%" mechanism/mechanism_3_full_chain.fac                 # → 203
awk '/^VARIABLE/{f=1;next} /^[[:space:]]*;/{f=0} f' mechanism/mechanism_3_full_chain.fac \
  | tr -s ' \t' '\n' | grep -c '^[A-Za-z][A-Za-z0-9_]*$'            # → 111

# 生成码与机制一致性
grep ":: ns" facsimile_impl/gen_v3/oz_mech.f90                    # → ns=111, nr=203

# 死端扫描（只产不耗 / 无排放无干沉）: 跑本文 §2.1 末尾的 python3 脚本
```

---

## 13. 参考文献

### 机制方案族
- **MCM v3.3.1**：Jenkin, M.E. et al. (1997) *Atmos. Environ.* 31, 81；Saunders, S.M. et al. (2003) *ACP* 3, 161；Jenkin, M.E. et al. (2003) *ACP* 3, 181（芳香族）；Bloss, C. et al. (2005) *ACP* 5, 641；Jenkin, M.E. et al. (2012) *ACP* 12, 5275（β-石竹烯）；Jenkin, M.E. et al. (2015) *ACP* 15, 11433（异戊二烯）
- **RADM2**：Stockwell, W.R. et al. (1990) *JGR* 95, 16343
- **CB6R3**：Yarwood, G. et al. (2005)；Cao, L. et al. (2021) *ACP* 21, 12687
- **MOZART-4**：Emmons, L.K. et al. (2010) *GMD* 3, 43
- **GEOS-Chem fullchem**：Sherwen, T. et al. (2016) *ACP*（统一卤素）；Fisher, J.A. et al. (2016) *ACP*（MARV 异戊二烯）

### 速率与光解数据
- **JPL 推荐值**：Burkholder, J.B. et al. (2019) *Chemical Kinetics and Photochemical Data for Use in Atmospheric Studies*, Evaluation No. 19, JPL 19-5
- **IUPAC 推荐值**：Atkinson, R. et al. (2004) *Atmos. Chem. Phys.* 4, 1461（气相有机反应）
- **异戊二烯臭氧化 OH 产率**：Atkinson, R. & Arey, J. (2003) *Chem. Rev.* 103, 4605
- **单萜臭氧化**：Lee, A. et al. (2006) *J. Phys. Chem. A*（α-pinene OH 产率 ~0.7–0.9）
- **光解参数化**：MCM 的 `l·cos^m·exp(−n secχ)` 形式；Saunders et al. (2003) *ACP*

### 指示因子与诊断
- **Sillman, S.** (1995) *JGR* 100, 14175（FNR 指示剂思想）
- **Tonnesen, G.S. & Dennis, R.L.** (2000)（HCHO/NO₂ 比值形式）
- **Kleinman, L.I.** (2005) *JGR* 110, D23302（L_NOx/Q 判据）
- **Duncan, B.N. et al.** (2010) *ACP*（OMI 卫星化）
- **Jin, X. et al.** (2017) *GRL*（阈值纬度/季节漂移）

### 工具与模型
- **CVODE/SUNDIALS**：Hindmarsh, A.C. et al. (2005) *ACM TOMS* 31, 363
- **KPP**：Damian, V. et al. (2002) *ACM TOMS* 28, 265
- **AtChem2**：Sommariva, R. et al. (2020) *GMD* 13, 1693
- **CMAQ**：Byun, D. & Schere, K.L. (2006) *Appl. Mech. Rev.* 59, 51
