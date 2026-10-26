# ozone_lab 自纠自查与修复工程全记录（R1–R16）

> **文档定位**：本项目 14 轮系统性自纠自查的**完整技术档案**——每一个被发现的错误、
> 其根因、影响评估、修复方案、验证证据，以及为防止复发而建立的永久防线。
>
> **适用版本**：**V2.3（74 物种 / 131 反应）** · 三路线（fac_v3 / fac_batch_omp / mcm_run）
> **最终状态**：functional_test **24/24 PASS** · 元素审计 **131/131 PASS** ·
> 冒烟 **1/1 PASS** · 从 `make clean` 可完全复现
> **编写时点**：2026-10-26（R16 收官验证通过，验收成功）
>
> ⚠ **本文档是 V2.3 时代的历史档案**。机制已于同日升级到 **V2.4（111 物种 / 198 反应）**
> （见 [mechanism/README.md](mechanism/README.md) §2.1）。本文中所有 "74sp/131rx"
> 计数均为**当时事实**，保留不篡改；当前规模请以 README/mechanism 文档为准。

---

## 目录

- [0. 总览：14 轮战果一览](#0-总览14-轮战果一览)
- [1. 审计方法论与六类探测器](#1-审计方法论与六类探测器)
- [2. R1–R2：防线的建立（哨兵 / 氮守恒 / FT 体系）](#2-r1r2防线的建立)
- [3. R3：数据通路截断与假绿（4 个真 bug）](#3-r3数据通路截断与假绿)
- [4. R4：光解映射回归与文档口径](#4-r4光解映射回归与文档口径)
- [5. R5：氮守恒漏权与双胞胎漂移](#5-r5氮守恒漏权与双胞胎漂移)
- [6. R6：陈旧计数与垃圾副本](#6-r6陈旧计数与垃圾副本)
- [7. R7：防线缺口与失效文档示例](#7-r7防线缺口与失效文档示例)
- [8. R8：生成码保真性终验（最深接口）](#8-r8生成码保真性终验)
- [9. R9：编译器输入守卫与单位缺口](#9-r9编译器输入守卫与单位缺口)
- [10. R10：情景缩放失效（P0）与 IPR 双计](#10-r10情景缩放失效与-ipr-双计)
- [11. R11：jfac 双乘（默认值掩盖 13 轮）](#11-r11jfac-双乘)
- [12. R12：构建系统与文档 API 的断裂](#12-r12构建系统与文档-api-的断裂)
- [13. R13：三位指数静默丢 E](#13-r13三位指数静默丢-e)
- [14. R14：终检（16 项验收全过）](#14-r14终检)
- [14b. R15：收官审计（3 问题 + 1 连锁发现）](#14b-r15收官审计3-问题--1-连锁发现)
- [14c. R16：收官验证（死开关 + 九项健康）](#14c-r16收官验证检验是否真的足够收官)
- [15. 永久防线总表（12 条）](#15-永久防线总表12-条)
- [16. 修复统计与方法学沉淀](#16-修复统计与方法学沉淀)

---

<a name="0-总览14-轮战果一览"></a>
## 0. 总览：14 轮战果一览

| 轮次 | 类型 | 抓到的核心问题 | 级别 |
|---|---|---|---|
| **R1–R2**（前序） | 建设 | 建立**哨兵体系 / 氮守恒不变量 / FT 框架**（哨兵首战即抓出 OPE 双计；氮守恒首战抓出 3 物种漏权） | 建设 |
| **R3** | 自查+修 | **BC 行 1024 字符截断**（致命）、**53-spec 假绿**、**KET 速率错挂 40×**、**TERPCHO 无汇** | P0×2, P1×2 |
| **R4** | 自查+修 | 光解映射回归验证（通过）、comic_H 口径、根目录死 `runtime.dat`、手册标题矛盾 | P2×2, P3×2 |
| **R5** | 自查+修 | **氮守恒 4 物种漏权**（PNA/RONO₂/OLENIT/TERPNIT）、**ATOL 双胞胎漂移** | P1×2 |
| **R6** | 自查+修 | `dsum(15)` 定长疑点（**排除**）、indicators NOy 自动覆盖（**确认**）、**14 处陈旧计数**、2 个垃圾副本 | P3×16 |
| **R7** | 自查+修 | **FT0 只查 ns 不查 nr**（防线缺口）、手册 §7 示例路径失效、架构总览缺指针、attic bak | P2×2, P3×2 |
| **R8** | 自查 | 生成码保真性终验（分数计量 / 速率乘数 / JARR 映射）——**零新 bug** | 验证 |
| **R9** | 自查+修 | **mechc 无重复物种/速率守卫**、**cfgbuild 审计通过**（两道 errs 门超预期）、**单位缺口 2 项**（BC 内联后缀 / 排放 ppb-h）、单位矩阵文档 | P2×2, P3×2 |
| **R10** | 审计+修 | **scen_build 漏缩 6 物种**（政策沙盘失效，P0）、**IPR 12–18h 双计**、**IPR 尾逗号空列**、P2 四处边界、P3 三项 | **P0×1**, P1×1, P2×4, P3×3 |
| **R11** | 审计+修 | **jfac 双乘**（默认值 1.0 掩盖 13 轮，P0）、**scan_parallel 死代码**、`nfail` 不退出、P3 三项 | **P0×1**, P1×1, P2×2, P3×3 |
| **R12** | 审计+修 | **batch→figures 数据链断裂**、**README 引用 6 个已删 API**、`make nc` 必崩、`cfgbuild` 默认指向遗产、30 个陈旧产物、`clean` 漏清 | P1×2, P2×5, P3×1 |
| **R13** | 审计+修 | **三位指数静默丢 E**（7 脚本无保护）、sat2surf 硬编码、`fig_ridge5` NameError、`sys.path` 算错 | P1×1, P2×2, P3×1 |
| **R14** | 终检 | **零新问题**——16 项验收全过 | 验证 |
| **R15** | 收官审计 | README 断链、`clean` 漏 `mcm_city`、**`mcm_run` 无构建目标（可复现性）**、**mechc 不支持行内段头** | P2×2, P1×2 |
| **R16** | 收官验证 | **`MODE=DIRECT` 死开关**、`RADIAT` 死配置；速率/O₃收支/ng边界/网格维全部健康 | P2×1, P3×1 |

**累计**：真 bug **34 个**（P0×4 / P1×8 / P2×15 / P3×7）、排除误报 3 项、
新防线 **12 条 + 7 条语义断言**、文档漂移修正 **40+ 处**、垃圾清理 **6 类**。

---

<a name="1-审计方法论与六类探测器"></a>
## 1. 审计方法论与六类探测器

十四轮下来，真正有效的 bug 探测器只有六类。每一类都至少抓到过一个真 bug：

### 1.1 ① 金值哨兵（golden-value sentinels）

**做法**：在响应面特定格点记录"金值"，测试时比对（±10–15%）。

```python
SENT = [
    ("2.28", 4, 132.10, 0.10),   # O3max 屋脊峰
    ("3.42", 4,  93.68, 0.10),   # O3max 饱和区
    ("1.01", 8,   3.539, 0.15),  # FNR 转折带
    ("2.28", 16,  6.487, 0.15),  # OPE 饱和区谷 ← "双计 bug 的教训位"
]
```

**战果**：首战抓出 **OPE 的 P(NOz) 双计 HNO₃**（详见 §2.1）。
**局限**：只能抓"数值漂移"，抓不住"跑的其实是旧代码"（R3 假绿即此）。

### 1.2 ② 守恒不变量（conservation invariants）

**做法**：在运行时检查物理守恒量，违反即报警。

```fortran
! 氮守恒: Σ w_N · dy/dt(化学) ≈ 0  (每线程首快照, 阈值 1e-4)
nsum = dot_product(wn, ych)
if (abs(nsum) > 1.0d-4*max(l1, 1.0d0)) print *, "[氮守恒 FAIL]", abs(nsum)
```

**战果**：
- 首战抓出 **ClONO₂ / BrONO₂ / ISOPNO₃ 未计入 NOy**（曾把 O₃/NOy 与 OPE 悄悄带偏）
- 再战抓出 **R5 的 V2.3 四个新 N 物种漏权**（PNA/RONO₂/OLENIT/TERPNIT）

**这是本仓库最高价值的探测器**——它抓的是"机制科学性"层面的错误，而非数值噪声。

### 1.3 ③ 等价路径测试（equivalence testing）

**做法**：同一条物理语义的两种实现路径，结果必须相同。**无需金值**。

**战果**：抓出 **R11 的 jfac 双乘**——`cfg JFAC=0.8 + CLI 1.0` 与
`cfg 1.0 + CLI 0.8` 本应完全等价，实测差 8.07%（`J ∝ CLI × JFAC²`）。

**为什么强**：不依赖"正确值是多少"，只依赖"两条路径必须一致"，因此**能发现任何路径实现错误**。

### 1.4 ④ 生成码一致性断言（code-gen parity）

**做法**：断言"生成码里的常量"与"源文件里的声明"一致。

```python
# FT0: gen_v3 的 ns/nr 必须等于机制 VARIABLE 物种数 / 反应行数
assert _ns == _n_sp, "gen_v3 陈旧! 先 make v3"
assert _nr_gen == _n_rx, "gen_v3 陈旧!"
```

**战果**：R3 的 **53-spec 假绿**（测试全绿的其实是陈旧二进制）；R7 补 nr 守卫。

### 1.5 ⑤ 数值路径全扫（exhaustive value sweep）

**做法**：不只看主物种，遍历所有物种/所有时刻，找极端值。

**战果**：
- R13：全扫 74 物种 × 25 时刻 → 发现 O₁D 出现 `3.07932-174`（三位指数丢 E）
- 负浓度检查：全扫发现仅 O/O₁D 两个痕量物种有负值（相对 <1e-7，正常）

### 1.6 ⑥ 静态契约检查（static contract check）

**做法**：解析构建文件 / 文档 / 配置，比对"声称的"与"实际的"。

**战果**：
- R12：解析 Makefile，发现 `batch` 写 A 而 `figures` 读 B（数据链断裂）
- R6/R9/R12：grep 出 40+ 处陈旧计数（53sp/67sp/96rx/97 反应）
- R9：`mechc.py` 无重复物种/速率守卫 → 加守卫

### 1.7 探测器选择原则

| 症状类型 | 该用哪个探测器 |
|---|---|
| 数值不对但不知道差多少 | ① 哨兵 |
| 物理量不守恒 | ② 守恒不变量 |
| 参数只在特定取值下出错 | ③ 等价路径（**首选**） |
| 改了源码但行为没变 | ④ 生成码一致性 |
| 偶发崩溃 / 极端值 | ⑤ 全扫 |
| 文档与代码不符 / 流程断裂 | ⑥ 静态契约 |

---

<a name="2-r1r2防线的建立"></a>
## 2. R1–R2：防线的建立（哨兵 / 氮守恒 / FT 体系）

这两轮不是"修 bug"，而是**建立发现 bug 的能力**。它们本身已各抓出一个真 bug。

### 2.1 哨兵首战：OPE 的 P(NOz) 双计 HNO₃

**症状**：`OPE = P(O₃)/P(NOz)` 系统性偏低 ~27%，但没有任何测试报警。

**根因**（`fac_batch_omp.f90` 的 `batch_budget_init`）：

```fortran
! 旧实现（错误）
pnoz = phno3                              ← ① HNO₃ 产率先算进 P(NOz)
do j = 1, nnoy
    select case (trim(sp_name(noy_id(j))))
    case ("NO", "NO2", "N2O5")
        cycle
    case default
        ! 遍历所有含 N 物种的生成反应
        ... pnoz = pnoz + ...             ← ② 循环里 HNO₃ 又加一次
    end select
end do
```

**量化**：`ΣP(NOz)` 旧 2.65e7 vs 正确 2.09e7，**差值恰为 24h ΣP(HNO₃)**。

**修复**：预算表重构为"正确枚举 14 项"（NOy 物种的生成反应各计一次，不预置 phno3）。

**修复后**：fvoc=1 切面 OPE 从 `148.9 → 52.1 → 13.4（饱和区谷）→ 高 NOx 回升`
——曲线形状终于物理自洽（NOx 饱和区 P(NOz) 占比最大 → OPE 触底）。

**留下的防线**：4 个金值哨兵，其中 `("2.28", 16, 6.487)` 就是**这个 bug 的教训位**。

### 2.2 氮守恒不变量首战：3 物种漏权

**症状**：加入运行期氮守恒自检后，立即在多个格点打印 `[氮守恒 FAIL] |Σ w_N·dy/dt| ~1e3`。

**根因**：`noy_add` 的权重表（当时 13 物种）漏了 **ClONO₂ / BrONO₂ / ISOPNO₃**——
卤素硝酸储库与异戊二烯硝酸过氧都含 N，但没被计入 NOy。

**影响**：
- 运行期守恒违反（N 原子凭空消失）
- **O₃/NOy 与 OPE 被悄悄带偏**（分母少了这三个物种的贡献）
- 且是"静默"的——没有任何报错，只是指示因子数值不对

**修复**：`noy_add` 补 3 行 → NOy 13 物种 → 闭合告警消失。

**留下的防线**：运行期氮守恒不变量（每线程首快照检查）。

### 2.3 自检轮次 3 的另外三项

| 项 | 问题 | 修复 |
|---|---|---|
| 诊断热路径 | 每次调用做名称多重集匹配（O(反应数×物种数)） | `batch_budget_init` 启动期一次性解析为索引表；热路径 O(项数) 求和 |
| 图不可复现 | 两张关键图内联生成（代码在临时脚本里） | 落盘 `figures/fig_indicators_ridge.py` / `fig_response_v3.py` |
| 测试计数口径 | `n_pass` 提前求和导致 16/20 的假数 | 汇总移到末尾重算 + 新增 4 项批量自检（6×6=36 案例） |

### 2.4 自检轮次 4：结构性预防 + 回归哨兵

| 项 | 内容 |
|---|---|
| **运行期氮守恒不变量** | 上线（见 §2.2），上线即抓出 3 物种漏权 |
| **诊断热路径索引化** | 实施时被自己的全局替换咬了两口（自赋值 `spO3=spO3` 与守恒块误插声明区）——均被编译警告/越界垃圾当场暴露；**教训：字符串替换后必查自赋值** |
| **回归哨兵** | 4 个金值格点（屋脊峰 29.0 / 饱和区 12.55 / FNR 转折 2.109 / OPE 谷 13.35——后者正是双计 bug 的教训位） |
| **`make figures`** | 两图从 CSV 一键重生成，可复现闭环 |

---

<a name="3-r3数据通路截断与假绿"></a>
## 3. R3：数据通路截断与假绿（4 个真 bug）

**背景**：V2.3 机制升级（67→74 物种、120→131 反应）刚完成。
本轮第一次系统审计"机制升级后的全链一致性"。

### 3.1 🔴 P0-A：BC 行 1024 字符截断（致命）

**现象**：`fac_v3` / `fac_batch_omp` 启动即崩：

```
At line 208 of file ../modules/runtime_dat_mod.f90
Fortran runtime error: Bad integer for item 1 in list input
```

**诊断过程**（教科书级的"逐层缩小范围"）：

1. 先怀疑 `runtime_v3.dat` 内容——用 Python 复刻 `read_runtime` 的解析逻辑，逐 token 校验 `int_val` → **全部合法**
2. 再编一个只调 `read_runtime` 的探针程序 → 崩溃复现，但 Python 复刻不崩
3. 给 `int_val` 加插桩（打印失败 token）→ 抓到**垃圾 token**（未初始化内存）
4. 给主循环加全轨迹打印 → 定位到 `BC` 段的第一行，`nt=102`（应为 149）
5. 量出 **BC 行实际长度 1117 字符 > `character(len=1024) :: line`** → **行被截断**

**根因**：

```fortran
character(len=1024) :: line       ! 固定 1024 字符行缓冲
```

V2.3 把物种数从 67 增到 74，BC 行 = `tau_mix + 74×(bc_c) + 74×(bc_s)` ≈ **1117 字符**
→ 超出 1024 → 尾部 47 个物种的 `bc_s` 读到**未初始化内存**（垃圾 token）
→ `int_val(垃圾)` 崩溃。

**为什么之前 67 物种时没暴露**：67 物种的 BC 行 ≈ 1140 字符（同样超！）——但垃圾值
恰好 >`nser`，走了"常量兜底"分支，**静默降级为"无序列约束"**。所以它潜伏了很久。

**影响**：

| 影响面 | 说明 |
|---|---|
| **序列 BC 观测同化** | 只有前 27 个物种能正确拿到序列索引，后 47 个静默失效 |
| fac_v3 / fac_batch_omp | 启动崩溃（74 物种时垃圾值恰好落入非法整数） |
| OSSE 实验 | 部分物种的观测约束根本没生效 |

**修复**（双保险）：

```fortran
! ① 行缓冲 1024 → 4096
character(len=4096) :: line   ! V2.3: 74 物种 BC 行 ~1117 字符 > 旧 1024 截断

! ② BC 段硬断言: token 数必须 = 1 + 2*ns
case ("BC")
    if (nt /= 1 + 2*ns) then
        print *, "[read_runtime FATAL] BC 行 token 数 ", nt, " != ", 1 + 2*ns, " (行截断?)"
        error stop
    end if
```

**验证**：探针复跑 → `read OK. ns=74 sp_name(74)=TERPCHO`、零 FAIL。

**留下的防线**：BC 行 token 数断言（将来扩到 100 物种也不会再静默截断）。

### 3.2 🔴 P0-B：53-spec 假绿

**现象**：修完 BC 截断后，`make test` 报 **10/10 全绿**，但 `fac_v3` 仍崩溃。

**根因**：`gen_v3/oz_mech.f90` 里的 `ns=53`（比 V2.2 的 67 还老）。
即**二进制是 53 物种时代的陈旧产物**，而 10/10"全绿"是跑在这个陈旧二进制上的。

**为什么哨兵抓不住**：哨兵读的是 O₃/HCHO 等主物种列——这些物种在 53/67/74 三个版本里
**都存在且行为相近**，所以数值在容差内。测试全绿但跑的是旧代码。

**修复**：手动重跑 `mechc.py` → `make v3` 强制重建 → 验证 CSV 表头 **75 列 = h + 74 物种**
（TERP/RONO₂/OLENIT/TERPCHO/TERPNIT/TERPOOH 全部到位）。

**留下的防线**：**FT0 生成码一致性断言**（后续 R7 又补了 nr 守卫）：

```python
_ns = int(re.search(r"::\s*ns=(\d+)", gen_src).group(1))
_nr = int(re.search(r"\bnr=(\d+)", gen_src).group(1))
assert _ns == 机制物种数, "gen_v3 陈旧! 先 make v3"
assert _nr == 机制反应数, "gen_v3 陈旧!"
```

### 3.3 🟠 P1-A：KET+OH 速率错挂 40 倍

**现象**：写"零引用命名速率扫描"时发现 `KKETOH = 2.2D-13` 定义了却零引用。

**根因**（V2.2 批量写集总类反应时的复制粘贴错误）：

```fortran
% KRO2NO : KET + OH = KETRO2 + H2O  ;    ← 用的是 RO2+NO 的速率(≈8.7e-12)
```

`KKETOH = 2.2e-13`（丙酮+OH，MCM 值）vs `KRO2NO ≈ 8.7e-12` → **快 39.5 倍**。

**影响**：丙酮（KET）的 OH 消耗速率虚高 40 倍 → 白天 HOx 汇被高估、
KET 在大气中寿命被低估 40 倍。

**修复**：反应行改挂 `KKETOH`。

**留下的防线**：零引用命名速率扫描（递归可达性分析，排除 Troe 展开参数）。

### 3.4 🟡 P2：TERPCHO 无汇

**现象**：V2.3 注释里承诺"TERPCHO 汇走 drydep.config"，但**没写**。

**影响**：SOA 代理终产物在气相机制中**无任何去除通道** → 5 天无界累积。

**修复**：
```text
input_v3/drydep.config 补 5 行:
  TERPCHO  0.10  on
  TERPNIT  0.30  on
  TERPOOH  0.30  on
  RONO2    0.25  on
  OLENIT   0.25  on
input_v3/emissions.config 补 TERP 行
build_groups 把 TERP 纳入 E_VOC 组
```

**验证**：IPR 沉降账 5 天积分 = **−4.128e8**（负值 = 汇生效）、
OLEOOH 类同族物种亦补齐。

**教训**：**"窗口内浓度平稳"不是"有汇"的证据**（慢产生的物种短窗内同样近似平稳）；
必须看 **IPR 沉降账是否为负**。

### 3.5 R3 的连带发现：OLEOOH 零汇（引用中已修，但漏了一步）

**背景**：用户引用的报告称 OLEOOH 零汇已修（drydep 加行），但实测 **IPR D 账 = 0.000**。

**核查**：`input_v3/drydep.config:32` 确实有 `OLEOOH 0.5 on`——但
**`runtime_v3.dat` 的 DRY 段没刷新**（改了配置没跑 cfgbuild）。

**修复**：`make runtime_v3.dat` → IPR D 账 = **−3.695e10** ✓（与引用声称的 −3.7e10 吻合）

**留下的教训**：**"零编译"≠"零重编配置"**——drydep/排放/BC 改动后必须刷新 runtime。

---

<a name="4-r4光解映射回归与文档口径"></a>
## 4. R4：光解映射回归与文档口径

### 4.1 ✅ 光解映射交叉验证（通过）

**背景**：历史上曾有"photolysis MAPPING 把 NO₂ 光解映射到 O₃→O¹D 参数"的 bug（差 ~200 倍），
已修复。本轮做**回归验证**。

**方法**：机制 J 槽位 ↔ `photolysis.config` 的 JMAP ↔ `photolysis-rates_MCM` 表，逐槽交叉核对。

| 槽位 | 反应 | 参数来源 | 判定 |
|---|---|---|---|
| 1 / 2 | O₃→O¹D / O₃→O(³P) | MCM 表 | ✓ |
| 4 | NO₂→NO+O | MCM 表 | ✓ |
| 5 / 6 | NO₃→NO / NO₃→NO₂+O | MCM 表 | ✓ |
| 7 / 8 | HONO / HNO₃ | MCM 表 | ✓ |
| 11 / 12 | HCHO radical / molecular | MCM 表 | ✓ |
| 13–29 | CH₃CHO/N₂O₅/KET/GLY/ISOPOOH/MGLY/ClO/BrO/ClONO₂/BrONO₂/HOCl/HOBr | JMAP `D` 直达值 | ✓ |
| 41 | CH₃OOH | MCM 表 | ✓ |

**结论**：恒等映射保持正确，历史修复未回退。

### 4.2 文档口径修正（4 项）

| 项 | 问题 | 修复 |
|---|---|---|
| comic_H 图 | 仍画"五条求解路线"（fac_run/kpp 已退役） | 并行会话已改为三路线；重生成验证 |
| 根目录 `runtime.dat` | fac_run 时代遗物，零代码引用 | 删除 |
| `OPERATIONS_MANUAL` 标题 | "五路线手册·5 可执行体"与 3 行路线表自相矛盾 | 改"三路线手册·3 可执行体" |
| `scan_suite.py` 注释 | "fac_run 保留为冒烟锚点"（已删除） | 改"已退役，由 batch+fac_v3 接管" |

---

<a name="5-r5氮守恒漏权与双胞胎漂移"></a>
## 5. R5：氮守恒漏权与双胞胎漂移

### 5.1 🟠 P1：氮守恒 4 物种漏权（R1-R2 同类问题复发）

**现象**：`scan_suite.py` 运行时打印 4 次 `[氮守恒 FAIL] |Σ w_N·dy/dt| = 800–3000`。

**根因**：V2.3 新增的 4 个含氮物种 **PNA / RONO₂ / OLENIT / TERPNIT** 没进
`fac_batch_omp.f90` 的 `noy_add` 权重表（仍是 V2.2 的 13 物种）。

**这是 §2.2 同一类 bug 的复发**——说明"新加含 N 物种就要补 noy_add"这条规则
在机制升级时被遗漏。

**修复**：
```fortran
call noy_add("PNA", 1)      ! HO2NO2 温度计
call noy_add("RONO2", 1)    ! ALK 烷基硝酸酯 (15% 分支)
call noy_add("OLENIT", 1)   ! OLE 硝酸酯 (3% 分支)
call noy_add("TERPNIT", 1)  ! 单萜硝酸酯 (NO3 加成)
```
→ NOy **13 → 17 物种**，氮守恒 FAIL = 0。

**连带修复**：OPE 分母（`pnoz` 产率枚举按 NOy 物种名自动扫描）随之自动修正 ——
新的硝酸酯生成反应被正确计入 P(NOz)。哨兵 6.881 → 5.933（∈ 6.487±15% ✓）。

**留下的防线**：
- 运行期氮守恒不变量（再次验证其价值）
- `mechanism/README §12.1` 的"五步闭环"清单明确列出"若带 N → 补 noy_add"

### 5.2 🟡 P2：ATOL 双胞胎漂移

**现象**：`fac_v3` 的 `tune_atol_by_name` 已补 V2.3 新物种（R3 时加的），
但 `fac_batch_omp` 的 `tune_atol_` 是**独立的副本**，没有——违反"A1/A2 统一"哲学。

**影响**：同一机制在两个驱动里用不同的 ATOL 调优表 → 数值路径不完全一致
（虽然影响小，但破坏了"双胞胎"概念）。

**修复**：**单源化**——`tune_atol_by_name` 移入 `chem_phys_mod`（共享物理层），
`fac_v3` 删本地副本、`batch` 转调共享版。

```fortran
! chem_phys_mod.f90 (新增)
subroutine tune_atol_by_name(a)
    select case (trim(sp_name(i)))
    case ("CH4"); a(i) = 1.0d10
    ...
    case ("RONO2","OLENIT","TERPNIT","TERPOOH"); a(i) = 1.0d7
    end select
end subroutine
```

**留下的教训**：**同一张表出现在两个文件里 = 漂移风险**。凡是"双胞胎共享语义"的表，
都应单源化。

---

<a name="6-r6陈旧计数与垃圾副本"></a>
## 6. R6：陈旧计数与垃圾副本

### 6.1 排除的两项（查实无 bug）

| 嫌疑 | 查证方法 | 结论 |
|---|---|---|
| `dsum(15)/dg(15)` 定长数组在 ns=74 下越界 | 逐行追写端：`d(1)=H₂O₂/HNO₃ … d(15)=P(HNO₃)`——15 个**固定语义诊断槽位**，不随物种数缩放 | ✅ 无越界 |
| `indicators.py` 的 NOy 判定漏新物种 | 读源：从 `.fac` FORMULA 段**自动**取含 N 物种（N 计数正则）——新 N 物种零改动自动入表 | ✅ 自动覆盖 |

### 6.2 修复：14 处陈旧计数

`ipr_irr5.py`×5、`hourly_pa.py`×1、`fig_ridge5.py`×1、`OPERATIONS_MANUAL.md`×2、
`README.md`×5 —— 全部 `53sp/67sp/96rx/97 反应` → `74sp/131rx`。

### 6.3 清理：2 个垃圾副本

| 文件 | 问题 |
|---|---|
| `mcm_city/input/mechanism_mcm_complete.fac` | 同机制的重排序副本（255 行 diff 仅物种排列序）→ 漂移风险 |
| `input_v3/mechanism_3_full_chain.fac` | 52sp 过时副本（缺 O1D 修复、头注 52 物种）→ 无人引用 |

**教训**：**"副本"是 bug 的温床**——机制文件、ATOL 表、组映射表都应单源。

---

<a name="7-r7防线缺口与失效文档示例"></a>
## 7. R7：防线缺口与失效文档示例

### 7.1 🟡 P2：FT0 只查物种数不查反应数

**问题**：`FT0` 只断言 `ns == 机制物种数`。若**往机制加反应而物种不变**（很常见），
陈旧二进制不会被拦截。

**修复**：补 nr 守卫（首版正则还写错——`::\s*nr=` 匹配不到，因为 nr 与 ns 同行且前无 `::`；
自查修正为 `\bnr=(\d+)`）。

```python
assert _ns == _n_sp, "gen_v3 陈旧! ns=..."
assert _nr_gen == _n_rx, "gen_v3 陈旧! nr=..."
```
实测：`gen ns=74/nr=131 vs 机制物种 74/反应 131` ✓

### 7.2 🟡 P2：手册 §7 示例路径失效

**问题**：`OPERATIONS_MANUAL §7` 的指示因子示例引用 `mechanism/ozone.fac`——
该文件已随 15sp 层退役删除，**照抄会报错**。

**修复**：改用 `mechanism_3_full_chain.fac` + `output/expA_BAU`，并**真跑验证**
（12 keys，FNR=3.427）。

**教训**：**文档里的命令必须被实际执行过一次**——否则就是"看起来能跑"的假文档。

### 7.3 🟢 P3：架构总览缺指针 + attic bak

`OPERATIONS_MANUAL` 架构总览章的 mechanism 行补规格书指针；
根目录 `attic_mcm_ch4_stale.fac.bak`（9 KB 陈旧备份）删除。

---

<a name="8-r8生成码保真性终验"></a>
## 8. R8：生成码保真性终验（最深接口）

**背景**：审计器读 `.fac`（分数计量 ✓），但运行时读的是 **`gen_v3/` 的生成码表格**。
若 `mechc.py` 把 `0.85`/`0.26` 截断成整数，运行时化学就偏离了被审计的机制。
这是"审计体系"与"运行时"之间的**最后一道未验证接口**。

**方法**：从 `runtime_v3.dat` 的 SPECIES 行 + `gen_v3/oz.reac`/`oz.prod` 反解每个反应的计量。

### 8.1 ✅ 分数化学计量保真

反应 53（ISOP+O₃）的产物：

```
29(MVK×0.5)  33(MACR×0.5)  10(HCHO×0.57)  19(CH3CHO×0.15)  6(OH×0.26)
```

五个分数**全部保真**（与被审计的 .fac 逐项一致）；CO₂/O₂ 作为储库被正确剔除
（非状态变量，不写进计量表）。

### 8.2 ✅ 分支速率乘数保真

```fortran
p(101) = 0.97*q(51)   ! OLERO2+NO = HCHO+CH3CHO+CO+HO2+NO2+H2
p(102) = 0.03*q(51)   ! OLERO2+NO = OLENIT
p(105) = 0.85*q(51)   ! ALKRO2+NO = CH3CHO+HCHO+HO2+NO2
p(106) = 0.15*q(51)   ! ALKRO2+NO = RONO2
```

四个乘数全保真（分支用**速率缩放**而非产物缩放，这是 V2.3 配平的正确写法）。

> 附带说明：审计时我最初的 grep 模式没匹配到生成码的 `q(51)` 引用形式，
> 一度以为乘数丢了——**虚惊一场，但这次虚惊促成了更严格的检查方式**。

### 8.3 ✅ JARR 槽位恒等映射

```fortran
p(1)   = JARR(4)    ! NO2 = NO+O
p(2)   = JARR(1)    ! O3 = O1D
p(129) = JARR(14)   ! N2O5 = NO2+NO3
```
与 `photolysis.config` 的"槽位恒等"语义严格一致。

### 8.4 ✅ 其它

| 项 | 结论 |
|---|---|
| batch J 通路 | `eval_J_jf(t, rp%jfac, JARR)` 每格点独立（当时判断正确，但实际有双乘 bug——R11 才发现） |
| fac_v3 IPR 导出列数 | **445 列 = 1 + 74×6** ✓ |

**结论**：R8 零新 bug，但建立了"生成码保真性"这一验证维度（后来 R13 用同样思路
发现了三位指数问题）。

---

<a name="9-r9编译器输入守卫与单位缺口"></a>
## 9. R9：编译器输入守卫与单位缺口

**背景**：`mechc.py` 与 `cfgbuild.py` 是九轮以来**唯一没通读过的两个文件**。

### 9.1 🟡 P2-A：mechc.py 无重复物种守卫

**问题**：`VARIABLE` 段的解析**没有去重**——同一物种写两次会导致
`ns` 虚增 + 索引错位 + 死变量。这正是历史上 `Br+O₃` 重复行的同类问题。

```fortran
! 旧实现
m.species += clean.rstrip(";").split()    ← 直接追加, 不管重复
```

**修复**（双重判定：跨行交集 + 同行自重复）：

```python
new_sp = clean.rstrip(";").split()
dup = sorted(set(new_sp) & set(m.species)) \
      + sorted({x for x in new_sp if new_sp.count(x) > 1})
if dup:
    raise SystemExit(f"[mechc] VARIABLE 重复物种: {sorted(set(dup))} "
                     "(ns 虚增 + 索引错位 — 修复机制文件)")
```

**自测**：构造 `NO NO2 NO` 的 .fac → 正确拒绝 ✓

### 9.2 🟡 P2-B：mechc.py 无速率重名守卫

**问题**：同名 `K` 定义两次 → **后者静默覆盖前者**（也是 `Br+O₃` 同类）。

**修复**：
```python
if any(nm == name for nm, _ in m.rates):
    raise SystemExit(f"[mechc] 重复速率定义: {name} "
                     "(后者静默覆盖前者 — 旧 Br+O3 同类 bug)")
```

**自测**：`KX` 双定义 → 正确拒绝 ✓

### 9.3 ✅ cfgbuild.py 审计通过（超预期）

| 此前担忧 | 实查 |
|---|---|
| `errs` 收集后有没有拦截门？ | **有两道门**（第一阶段 :275、光解 JMAP 后 :307）——"校验不过拒绝写出"承诺成立 |
| 排放日型大小写静默 FLAT 坑？ | **编译期 :212-216 已显式拦截**（非 BOX/FLAT 直接拒绝，且拒绝文案本身就警告了静默行为）；负排放、振幅 ∉[0,1] 同样拦截 |

### 9.4 🟡 P2-C/D：单位能力缺口 2 项（最重要）

**背景**：核查"程序能接受哪些单位、能否换算"。

**修复前的单位能力矩阵**：

| 通道 | 状态 |
|---|---|
| IC | ✓ cm³/ppb/ppm 每行后缀 |
| BC 序列 file: | ✓ 头部 `# unit:` 指令 |
| **BC 内联常数** | ✗ `float(v)` 直读——写 `40 ppb` 会崩或差 2.46e10 倍 |
| **排放** | ✗ 只认 molec/cm³/s（模板自述"待扩展"） |

**修复**：

```python
# ① BC 内联常数支持单位后缀
parts = v.split()
bc_c[sid[k]] = float(parts[0]) * unit_factor(parts[1] if len(parts) > 1 else "cm3")

# ② 排放支持文件头单位指令
# unit: ppb/h          ← 整份文件生效
e0 = float(p[1]) * (PPB / 3600.0 if emit_unit == "ppb/h" else 1.0)

# ③ unit_factor 单源化（read_series/read_ic/BC 三处共用）
def unit_factor(u):
    f = {"cm3": 1.0, "ppb": PPB, "ppm": 1e3*PPB}.get(str(u).lower())
    if f is None: raise ValueError(f"未知单位 '{u}' (cm3|ppb|ppm)")
    return f
```

**端到端自测**：
- `NO 10.0 BOX 0.9` + `# unit: ppb/h` → 基线 **6.8333e7** ✓（= 10 × 2.46e10/3600）
- `O3 = 40 ppb` → BC 常数 **9.8400e11** ✓

**留下的文档**：`UNITS.md`（449 行，后扩到 506 行）——14 个数据文件逐子量单位表 +
硬编码/半硬编码参数登记 + 换算公式汇编 + 校验行为清单 + 换算速查卡。

### 9.5 垃圾大扫除

`__pycache__`×3、`.DS_Store`×2、根目录 stale PNG、`output/` smoke 残留。

---

<a name="10-r10情景缩放失效与-ipr-双计"></a>
## 10. R10：情景缩放失效（P0）与 IPR 双计

**背景**：本轮换扫描面到"未深审的区域"——`runtime_dat_mod` 的情景生效路径、
`fac_batch_omp` 的诊断累积路径。

### 10.1 🔴 P0：scen_build 漏缩 6 个排放物种（政策沙盘失效）

**发现过程**（这是本仓库最戏剧性的一次）：

审计时读到 `scen_build` 调用 `build_groups()` 算出完整组数组，**却又用了一份独立的
`select case` 清单**：

```fortran
subroutine scen_build(isc)
    call build_groups()                       ← 计算了完整的 grp* 组数组
    do i = 1, ns
        select case (trim(sp_name(i)))        ← 但这里又用了一份「独立」清单
        case ("NO", "NO2")                        → E_NOx
        case ("HCHO","CH3CHO","ISOP","DMS","MVK","MACR","MSA")  → E_VOC
        case ("CO")                               → E_CO
        case ("CH4")                              → E_CH4
        end select        ← SO2/ALK/OLE/ARO1/ARO2/KET/TERP 不在任何 case!
    end do
end subroutine
```

而 batch 侧用 `esc_of(rp)` → 直接读 `grpvoc(k)` 等数组（**包含**这 7 个物种）。

**实测铁证**（BAU vs HVOC 的 IPR 排放账 5 天积分）：

| 物种 | BAU | HVOC（VOC×3） | 比值 | 判定 |
|---|---|---|---|---|
| **OLE** | 5.565e11 | 5.543e11 | **0.996** | ✗ 完全没缩放 |
| **ALK / ARO1 / ARO2 / KET / TERP / SO2** | — | — | **0.996** | ✗ 同上 |
| ISOP / HCHO | — | — | 2.988 | ✓ ×3 正常 |

**即：`emissions.config` 里排放量最大的四类 VOC（ALK 2.6e7、OLE 2.1e7、
ARO1 1.2e7、KET 1.5e6）在 fac_v3 的 HVOC/ISOPX 情景中完全不变。**

**最严重的影响在证据链之外**：`README §8` 的"体制解读翻转"论证用
**"HVOC 几乎无效（−0.2%）"** 作为 NOx-limited 的佐证——**该证据完全无效**（VOC 压根没加）。
batch 同情景是 **+8.0%**。此前我把 −0.2% vs +8.0% 解释为"网格最近点差异"，
**实为本 bug 的直接后果**。

**修复**（单源化）：

```fortran
subroutine scen_build(isc)
    ! R10-P0: 单源组映射 —— 只读 build_groups 的 grp* 数组
    call build_groups()
    do i = 1, ns
        if (grpnox(i) == 1) emis_scale(i) = scen_fnox(isc)
        if (grpvoc(i) == 1) emis_scale(i) = scen_fvoc(isc)
        if (grpco(i)  == 1) emis_scale(i) = scen_fco(isc)
        if (grpch4(i) == 1) emis_scale(i) = scen_fch4(isc)
    end do
end subroutine
```

**验证**：
- 全部 VOC 族 ×**3.038** ✓、NO/CO ×1.013（本就不该缩）
- **HVOC O₃max 101.4 → 111.26 ppb（+9.5%）**，与 batch 的 +8.0% 同向一致

**连带修复的文档**：README §8 六情景表、A2 表、"体制解读"论证重写
（明确标注"原 −0.2% 论据失效，结论不变但论据已换"）。

**新增防线**：`R10-P0 HVOC 缩放全部 E_VOC 组物种 (×3)`

### 10.2 🟠 P1：IPR 在 12–18h 窗口双计

**问题**：主循环在 12–18h 窗口内**第二次调用 `diag_indicators`**，
而该子程序对 `iprc` 的累积**不受 `pm` 门控** → 白天 6 个输出步被计两次。

```fortran
call diag_indicators(..., dg, ipracc, .false.)        ← ① 每步累积
if (hh >= 12 .and. hh <= 18) then
    call diag_indicators(..., dgpm, ipracc, .false.)  ← ② 窗内再调，同一 ipracc
```
`diag_indicators` 内 `iprc(...) = iprc(...) + ...` 无条件执行 → 双计。

**影响**：IPR 侧文件 P/L 绝对值偏高 **~25–30%**（不影响主 CSV 的 FNR/OPE 等）。

**修复**：
```fortran
subroutine diag_indicators(t, y, rp, d, iprc, pm, acc_ipr)
    logical, intent(in), optional :: acc_ipr
    logical :: do_ipr
    do_ipr = .true.
    if (present(acc_ipr)) do_ipr = acc_ipr
    ...
    if (do_ipr) then
        ! IPR 累积 (机制 P/L + 物理四账)
    end if
```
第二次调用传 `.false.`。

**验证**：O₃_P 恢复非零（4.42e8）、六账独立、**稳态残差 0.5%**。

### 10.3 🟠 P1（修复中发现）：IPR 尾逗号空列

修复 P1-1 后验证时发现 **IPR 数据行 449 列 vs 表头 448**：

```fortran
write (uipr, '(3(F6.3, ","), F6.1, ",")', advance="no") ...   ← 尾逗号
    + write (uipr, '(",", es13.5)')                            ← 再补前导逗号 → ",," 空字段
```

**后果**：**所有按列名索引的读取器整体错位一位**（审计时正是因此误判 `O3_P=0`）。

**修复**：去掉尾逗号（ipr + irr 两处），列数复验 **448/448、135/135**。

**新增防线**：`R10 batch_ipr/irr 表头/数据列数一致`

### 10.4 🟢 P2 四处边界保护

| # | 位置 | 缺口 | 修复 |
|---|---|---|---|
| P2-1 | `runtime_dat_mod:177` SERIES | `iser` 无上限检查（`ser_n(8)`） | `iser > NSER_MAX` → error stop |
| P2-2 | `fac_batch_omp:25` `noy_id(20)` | 17 个 NOy 已用，仅剩 3 空位且无检查 | 数组 **20→32** + `nnoy >= size(noy_id)` 检查 |
| P2-3 | `runtime_dat_mod:185` | 序列点只由 cfgbuild 保证 ≤128 | `ser_n > 128` → error stop |
| P2-4 | `fac_batch_omp:135` | `spO3==0` 时 `yp(0)` 越界读 | `min(spO3,…,spNO) > 0` 断言 |

### 10.5 🟢 P3 三项 + ⚪ P4

| 项 | 修复 |
|---|---|
| P3-1 旧计数 5 处（scen5/ridge5） | sed → 74sp |
| P3-2 `RUN` 行 `nser` 字段未校验 | 加 `RUN.ns /= ns` → error stop |
| P3-3 `sp_name` 兜底为 SPi（静默致组映射全失效） | 改为**硬失败** |
| P4 排放账 0.2% 波动 | 自检加"同情景两次跑逐位一致" → **确认确定性成立**，0.19% 差异是跨情景采样口径噪声，无需修复 |

---

<a name="11-r11jfac-双乘"></a>
## 11. R11：jfac 双乘（默认值掩盖 13 轮）

**背景**：R10 修完 scen_build 后，审计转向"参数缩放路径"与"实验脚本行为"。

### 11.1 🔴 P0：batch 的 jfac 被乘两次

**完整源码路径**：

```fortran
! fac_batch_omp.f90  batch_rates()
call eval_met_dt(t, rp%dts, TEMP, PRES, RH, BLH, RAIN, JF, dB)
                                       ↑ JF = eval_met 返回的 JFAC (= met_c(7) 或序列)
JF = rp%jfac*JF                        ← ① 第一次乘: CLI jfac × config JFAC
call eval_J_jf(t, JF, JARR)

! chem_phys_mod.f90  eval_J_jf()
jfac_scen = jf                         ← 把 JF 写入 jfac_scen
call eval_J(t, JARR)

! runtime_dat_mod.f90  eval_J()
jf = jfac_scen*met_c(7)                ← ② 第二次乘: 再次乘 config JFAC
JARR(i) = jf*j_l(i)*cosz**j_m(i)*exp(-j_n(i)/cosz)
```

**净效果**：`J ∝ CLI_jfac × (config_JFAC)²`

**决定性实验**（config JFAC=0.8 副本，网格 fnox=fvoc=1.0125）：

| 组合 | O₃max | 相对基线 | 应然 | 判定 |
|---|---|---|---|---|
| 基线 cfg1.0 × CLI1.0 | 102.48 | 1.000 | 1.000 | ✓ |
| A：cfg1.0 × CLI0.8 | **94.91** | 0.9261 | 0.926 | ✓（config=1 时双乘退化为单乘） |
| B：cfg0.8 × CLI1.0 | **87.25** | 0.8514 | 0.926 | ✗ **多降 8.1%** |

**影响**：

| 影响面 | 说明 |
|---|---|
| **jfac 网格维**（CLI 第 4 参） | 传 0.6（模拟 SCEN_CLOUD）时实际生效 **0.36** → O₃ 被压低远超预期 |
| **任何 `JFAC ≠ 1.0` 的 meteo.config** | batch 结果系统性偏低，且与 fac_v3 分叉 |
| **默认配置（JFAC=1.0）** | 双乘退化为单乘，**结果正确** → 这解释了为何 13 轮测试、40+ 次实验、4 个哨兵全部通过 |

**修复（方案 B：绝对语义）**：

```fortran
subroutine eval_J_jf(t, jf_abs, JARR)
    ! jf_abs 是「最终光解缩放因子」(已含 config JFAC × CLI 网格因子)
    ! 先算 config 侧基础值, 归一化掉 eval_J 内的 *cfg_jf
    if (met_s(7) > 0 .and. met_s(7) <= nser) then
        cfg_jf = ser_interp(met_s(7), t)
    else
        cfg_jf = met_c(7)
    end if
    jfac_scen = jf_abs/max(cfg_jf, 1.0d-30)   ! 抵消 eval_J 内的 *cfg_jf
    call eval_J(t, JARR)
    jfac_scen = 1.0d0                          ! 复位, 避免污染后续调用
end subroutine
```

**验证**：A 与 B 现在**逐位相同**（94.91 / **差 0.000%**，修前 8.07%）

**新增防线**：`R11-P0 jfac 等价路径 (cfg0.8×CLI1.0 == cfg1.0×CLI0.8)`
—— 无需金值，只断言两条路径等价。**这是本仓库最强的语义断言。**

### 11.2 🟠 P1：`perf/scan_parallel.sh` 死代码

**问题**：第 18 行调用 **已删除的 `fac_run`** → 脚本 100% 失败；而
`PERF_REPORT §II.5` 仍称它"现行（#4 进程并行扫描）"。

**修复**：删除脚本 + `PERF_REPORT §II.5` 改标"**载体已退役（R11 勘误）**"，
保留负载敏感性的方法学结论并指向 §III.1。

**后续（2026-10-26）**：整个 `perf/` 目录已删除（4 份报告，无代码）。
其中仍在用的结论已并入正文——G 类 6 项优化表在 `README.md`，
负载敏感性/CSV 行序/`netCDF4` 三个坑在 `OPERATIONS_MANUAL.md §14.5`。

### 11.3 🟡 P2-1：双胞胎末窗样本数不一致

| 路线 | 条件 | 样本数 |
|---|---|---|
| fac_v3 | `step >= nsteps-24` | **25**（step 96…120，含端点） |
| fac_batch | `k > nsteps-24` | **24**（k 97…120） |

**影响**：均值口径差 ~4%（1/24），是"双胞胎 ±10–15% 容差"的已知来源。

**决定**：**只改文档，不改代码**（改代码会触发全部哨兵重标定，收益不抵成本）。
写入 `PERF_REPORT §I.3` 与 `OPERATIONS_MANUAL §4`。

### 11.4 🟡 P2-2：batch 失败案例不影响退出码

**问题**：`nfail` 只 print，不 `stop 1` → 81 案例中部分失败时 CI 仍绿。

**修复**：
```fortran
if (nfail > 0) then
    print '(A,I0,A)', "[FATAL] ", nfail, " 个案例求解失败"
    stop 1
end if
```
+ `scan_suite.py` 错误提示点名"含部分格点求解失败 — 见 stdout 的 [FATAL]"。

### 11.5 🟢 P3 三项

旧计数 8 处（scen5 4 + ridge5 4）、`RUN.ns` 校验、`sp_name` 硬失败。

---

<a name="12-r12构建系统与文档-api-的断裂"></a>
## 12. R12：构建系统与文档 API 的断裂

**背景**：前 11 轮审计了机制/求解器/输入/诊断，但**从未审计"审计器、构建系统、文档 API"**
这一层。本轮 2 个 P1 全落在这里。

### 12.1 🔴 P1-A：batch 与 figures 的数据链断裂

**Makefile 事实**：

```make
batch:   ./fac_batch_omp 9 ../output/batch_v3.csv          ← 写 A
figures: output/batch_v3_indicators.csv                     ← 读 B
```

**时间戳实证**：
```
output/batch_v3.csv            Oct 26 05:30   ← make batch 的产物（陈旧）
output/batch_v3_indicators.csv Oct 26 07:44   ← scan_suite.py 的产物（新）
```

**结论**：`make batch && make figures` 会**用上一次 `scan_suite.py` 的陈旧 CSV 重画图**，
新跑出的 81 案例结果**完全不影响图**。而 README/OPERATIONS/Makefile 注释都写着
"`make figures` 从 batch 结果重生成图"——文档与行为不符。

**修复**：`batch` 输出统一为 canonical 名 `batch_v3_indicators.csv`。

**实测验证**：batch 07:54 写 → figures 07:55 读新数据出图 ✓

**新增防线**：`R12 batch 输出名 == figures 依赖名 (数据链不断)`
—— 静态解析 Makefile 比对两个名字。

### 12.2 🔴 P1-B：README 引用 6 个已删除的诊断 API

**事实**：`diagnostics/diagnostics.py` 现仅 **3 个公开符号**（并行会话重构后）：

```python
_fix_fortran_exp(tok)   # Fortran 指数修复
load(prefix)            # → (h, C, irr_h, I)
ipr_from_fac(prefix)    # 六账积分
IPR_PROC                # ["P","L","E","D","B","X"]
```

但 README §6 与 §D 类表仍宣称 **8 个函数**，包括已删除的
`indicators(C)` / `irr_integrate` / `irr_key_processes` / `families` / `chem_ledger` / `report`。

**实测**：`dg.report(...)` → **AttributeError: module 'diagnostics' has no attribute 'report'**

**修复**：README §6 与 §D 表重写为**双库架构**：

- `diagnostics/indicators.py`：`route_indicators(prefix, fac)` / `route_timeseries(...)`
  —— 路线无关 12 指标引擎（**主力**）
- `diagnostics/diagnostics.py`：`load` / `ipr_from_fac` / `IPR_PROC`
  —— CSV 读取与六账积分（**底层**）

+ 两处 R12 勘误说明。

### 12.3 🟡 P2 五项

| # | 问题 | 修复 |
|---|---|---|
| P2-1 | `make nc` 裸跑必崩（`NC_SRC` 未定义 → 指向已退役的 `output/cmp_FACSIMILE`） | 默认改 `output/expA_BAU`；实测生成 25t×74sp×131rxn |
| P2-2 | `cfgbuild` 默认指向 15sp 遗产 `input/` | 默认改 `input_v3` / `runtime_v3.dat` |
| P2-3 | `output/` 30 个陈旧测试产物（fac_run/EBI 时代） | 删除，仅留现行产物 + 4 个实验子目录 |
| P2-4 | `make clean` 漏清 6 类（`oz.reac`/`oz.prod`/`gen_v3/*.mod`/`rad_photol_demo`/`*.dylib`/`*.nc`） | 补齐 |
| P2-5 | `input/` 零消费者 | ~~**不删**（15sp 教学参照），新建 `input/README.md` 标注"非工程入口"~~ → **2026-10 决定改判：已彻底删除**。当时"保留+标注"的策略被证明无效——目录仍在，误用入口依旧存在；15sp 教学机制档移交 `mechanism/teaching/` 存档，其余 12 类配置为 `input_v3/` 的死副本，无保留价值 |

### 12.4 🟢 P3-1：孤儿二进制

`audit/liboz_audit.dylib`（零引用）→ 删除编译产物，保留 `.c` 源。
`modules/rad_photol_demo` 同理（`rad_photol.f90` 可随时重编）。

### 12.5 本轮确认无问题的 4 项

| 检查 | 结论 |
|---|---|
| **审计器自检**（注入 3 类错误到 /tmp 副本） | ✓ 正确抓到"未知物种 XYZ"与 H/O 不配平 |
| Makefile 源依赖完整性 | ✓ **含 `chem_phys_mod.f90`**（R5 改动会自动重编） |
| PPB 常量三处一致性 | ✓ 全部 `2.46e10` |
| `fig_*.py` 列索引 | ✓ 与 24 列表头逐项对应 |

---

<a name="13-r13三位指数静默丢-e"></a>
## 13. R13：三位指数静默丢 E

**背景**：审计转向**新代码**（`src/sat2surf.py`）与**跨工具数值一致性**。

### 13.1 🔴 P1：Fortran 三位指数导致 7 个下游脚本崩溃

**现象**：`fac_v3` 的 `<pref>_conc.csv` 中出现 **`3.07932-174`**（O₁D 列）。
`float('3.07932-174')` → ValueError。

**根因（Fortran 格式语义，实测拍板）**：

```
ES13.5   →  3.07932-174      ✗ 丢失 E
ES14.5   →   3.07932-174     ✗
ES15.5   →    3.07932-174    ✗
ES16.5   →     3.07932-174   ✗
ES17.6   →      3.079320-174 ✗
E15.5    →     0.30793-173   ✗（且尾数归一化方式不同）
ES15.5E3 →   3.07932E-174    ✅ 唯一可用
```

gfortran 的 `ES` 编辑符在 **|指数| ≥ 100** 时，因字段宽度不足以同时容纳 `E` 与 3 位指数，
**静默省略 `E`**。必须显式写 `ESw.dE3` 强制 3 位指数带 E。

**影响面**：

| 文件 | 出现次数 |
|---|---|
| `<pref>_conc.csv` | 2（O₁D 列） |
| `<pref>_ipr.csv` | 0 |
| `<pref>_irr.csv` | 0 |

**下游解析器保护矩阵**（14 个解析器）：

| 有保护（6） | 无保护（7） |
|---|---|
| diagnostics.py / indicators.py / nc_output.py / functional_test.py / ipr_closure_check.py / ipr_irr5.py | **scan_suite / scen5 / ridge5 / fig_ridge5 / hourly_pa / obs_osse / sens_adj** |

**为什么 24/24 测试全绿仍未发现**：所有测试读的是 O₃/HCHO 等主物种列（值域 1e8–1e12，
2 位指数），O₁D 只在极少数场景才跌破 1e-100。

**修复（双管齐下）**：

| 层 | 改动 |
|---|---|
| **源头** | 5 处 `es13.5` → **`es15.5e3`**（oz_fac_v3 三处 + fac_batch_omp 两处） |
| **下游** | 新建 **`diagnostics/csvfloat.py`**（`fix_fortran_exp`/`fnum`/`load_csv_rows`），7 个脚本全部接入 |

**源头验证**：重跑后 O₁D 列 = `3.07932E-174`、`-9.75647E-278`，三位指数 **0 次** ✓

**新增防线**：
```
R13 输出格式无三位无E指数 (ES15.5E3)      发现 0 处 ✓
R13 csvfloat 四种指数形态解析              4 组 ✓
```

### 13.2 🟡 P2-A：sat2surf 的 `2.46e10` 硬编码

7 处直写常量 → 改为 `from cfgbuild import PPB`（单源）。

### 13.3 🟡 P2-B：`oz.prod` 计量表格式（核查后安全）

`mechc.py` 用 `{c:g}`（Python 格式，自动带 E）——无需改。

### 13.4 🟢 P3：极痕量物种负值（文档化）

| 物种 | 最小值 | 相对自身峰值 |
|---|---|---|
| O | −4.016e−04 | −2.15e−08 |
| O₁D | −1.485e−71 | −1.71e−70 |

**判定**：BDF 不保正性是正常行为（10b 坑 #5），**不施加 `max(y,0)`**
（会引入 RHS 非光滑折点，曾致 15× 步数惩罚，10b 坑 #11）。已写入 `UNITS.md §8`。

### 13.5 修复过程中连带发现并修好的 3 个问题

| # | 问题 | 修复 |
|---|---|---|
| 1 | 7 个脚本的 `sys.path` 算错（指向 `experiments/diagnostics`） | 按 `basename(parent)=="experiments"` 判定上溯一级 |
| 2 | `fig_ridge5.py` 的 `NameError: axes`（`zip(PANELS, [axes])` 笔误 + 循环变量遮蔽） | 改 `axp`，**现已跑通出图** |
| 3 | 测试用 sat2surf demo 数据 `qa=0` → 被 QA>0.75 过滤剔光 | 改 0.9（产品正确，测试数据构造错误） |

---

<a name="14-r14终检"></a>
## 14. R14：终检（16 项验收全过）

| # | 验收项 | 结果 |
|---|---|---|
| ① | 全部可执行体在位 | **3/3** ✓ |
| ② | **`make clean` 后从零重建** | ✓ 在独立 `/tmp` 副本上验证：重建后 O₃max **101.61 ppb 与主线完全一致** |
| ③ | 副本上完整测试链 | **24/24 PASS** |
| ④ | 副本上元素审计 | **131/131 PASS** |
| ⑤ | 文档指针全存活 | **10/10** ✓ |
| ⑥ | 三路线数值独立 | fac_v3 101.61 / batch 同格点一致 |
| ⑦ | 旧计数残留（53sp/67sp/96rx/97） | **零残留** ✓ |
| ⑧ | TODO/FIXME/XXX | 唯一命中是测试**故意注入**的 `XXX_UNKNOWN` |
| ⑨ | 编译告警 | 仅 `Unused dummy argument 'y'`——**良性** |
| ⑩ | fac_v3 运行时告警 | **零** ✓ |
| ⑪ | batch 运行时告警 | **零** ✓ |
| ⑫ | **文档数值 vs 实测** | BAU 101.6=101.6 / HVOC 111.3=111.3 / ISOPX 73.1=73.1 全 ✓ |
| ⑬ | output 目录结构 | 18 项，仅现行产物 + 4 个实验子目录 |
| ⑭ | 垃圾清查（`.bak`/`.orig`/`*~`） | **零残留** ✓ |
| ⑮ | 缓存与临时文件 | 已清 |
| ⑯ | **三方计数一致** | `.fac` **74sp/131rx** = `gen_v3 ns=74,nr=131` = `runtime RUN 74` = `SPECIES 行 74` ✓ |

**结论：零新问题，验收成功。**

---

## 14b. R15：收官审计（3 个问题 + 1 个连锁发现）

本轮重点核验"刚写的文档自身准确性"与"全库最终可复现性"。

### 14b.1 🟡 P2：README 断链

`README.md:284` 仍指向已归档的 `AUDIT_R10.md`（实际在 `audit_reports/`）。
**修复**：改指 `AUDIT_COMPLETE.md §10` + `audit_reports/AUDIT_R10.md`。
**复验**：全库 46 处相对链接，**断链 0** ✓

### 14b.2 🟡 P2：`make clean` 未覆盖 `mcm_city/`

`clean` 目标漏清 `mcm_city/`（`.mod`/`mcm_run`/gen 产物）。
**修复**：补一行 `cd mcm_city && rm -f ...`。实测 clean 后 **零残留** ✓

### 14b.3 🔴 连锁发现：`mcm_run` 无重建目标（可复现性破坏）

补 clean 后暴露：**Makefile 从无 `mcm_run` 的构建规则**——二进制是人工编译的，
`make clean` 后**第三条路线永久丢失**。

**修复**：新增两个 Makefile 目标：
```make
mcm_city/gen/oz_mech.f90: mcm_city/input/mechanism_mcm_ch4.fac facsimile_impl/mechc.py
	cd mcm_city && $(PY) ../facsimile_impl/mechc.py input/mechanism_mcm_ch4.fac gen

mcm_city/mcm_run: mcm_city/mcm_driver.f90 mcm_city/gen/oz_mech.f90
	cd mcm_city && $(FC) $(FFLAGS) -J gen -I gen gen/oz_mech.f90 mcm_driver.f90 ... -o mcm_run

mcm: mcm_city/mcm_run
```

### 14b.4 🔴 连锁发现：mechc.py 不支持"行内段头"（R9 守卫误报的真相）

**现象**：`make mcm` 首次执行触发 R9 的重复物种守卫：
```
[mechc] VARIABLE 重复物种: ['H2', 'HCHO'] (ns 虚增 + 索引错位)
```

**误判过程**：先怀疑 MCM 机制文件有重复 → 核实 **29 物种全唯一** →
确认为 **mechc.py 解析器 bug**。

**根因**：MCM 文件把段头写在同一行：

```
HSO3 CO CL O HNO3 ... SA CH3O2 ;FORMULA      ← 第 26 行: 物种表 + 分号 + FORMULA 连写
FORMULA                                       ← 第 56 行: 又一个独立 FORMULA
HCHO       C1 H2 O1                           ← FORMULA 段数据
```

`mechc.py` 的段判定是**整行精确匹配**（`s == "VARIABLE"` / `s.startswith("FORMULA")`）：
- `... CH3O2 ;FORMULA` 既不是 `VARIABLE` 也不以 `FORMULA` 开头 → 落入 `in_vars` 分支，
  因含 `;` 而**提前结束 VARIABLE 段**（物种收对了）
- 但 **FORMULA 模式没开启** → 第 27 行起的元素组成行（`HCHO C1 H2 O1`）被当物种收集
  → 触发重复守卫误报

**影响**：`mcm_run` 从 clean 后**无法重建**（过去一直是手工编译，掩盖了这个 bug）。

**修复**（行内段头归一化）：

```python
# 段标记可能行内连写 ("... CH3O2 ;FORMULA"), 先归一化再判定
_seg = s
if ";" in s:
    _tail = s.split(";", 1)[1].strip()
    if _tail.startswith("FORMULA"):
        _seg = "FORMULA"          # 行内紧跟 FORMULA 段头
if _seg == "FORMULA":
    if in_vars:                   # 先收本行的 VARIABLE 物种
        _heads = s.split(";", 1)[0].split()
        if _heads: m.species += _heads
        in_vars = False
    in_form = True; continue
```

**验证**：
- MCM：`[mechc] 29 物种, 71 反应, 139 命名速率, nj=51` ✓（与原 gen 完全一致）
- V2.3 回归：`74 物种, 131 反应, 105 命名速率, nj=41` ✓（无回归）
- `clean → make v3 + make mcm` → 三二进制全重建、`mcm_run` 跑通（0.100s, nst=7350）

### 14b.5 R15 验收（收官）

| # | 项 | 结果 |
|---|---|---|
| 1 | 全库 46 处相对链接 | **断链 0** ✓ |
| 2 | `make clean` 后残留 | **0 个** `.mod`/`.o`/`.dylib` ✓ |
| 3 | **clean → make v3 + make mcm 全路线重建** | **3/3 二进制重建成功** ✓ |
| 4 | MCM 机制解析 | 29sp/71rx ✓（与历史 gen 一致） |
| 5 | V2.3 机制解析（回归） | 74sp/131rx ✓ |
| 6 | functional_test | **24/24 PASS** |
| 7 | 元素审计 | **131/131 PASS** |
| 8 | 冒烟 | **1/1 PASS** |
| 9 | 三层计数一致（.fac / gen_v3 / runtime / SPECIES） | 全部 74/131 ✓ |
| 10 | AUDIT_COMPLETE 统计自洽（4+6+12+6=28） | ✓ |
| 11 | 文档数值 vs 实测（BAU/HVOC/ISOPX） | 101.6/111.3/73.1 全 ✓ |

**R15 结论：3 个问题（1 断链 + 1 clean 遗漏 + 1 可复现性缺陷）全部修复，
并连带修复了 mechc.py 的行内段头解析缺陷。零遗留问题。**

---

## 14c. R16：收官验证审计（"检验是否真的足够收官"）

**姿态**：不重复已查项，专攻**从未被执行过的代码路径**、物理量级校验、
全配置项消费扫描。

### 14c.1 🟡 P2：`MODE = DIRECT` 是死开关（真问题）

**决定性验证**：`MODE=MCM` 与 `MODE=DIRECT` 的 O₃max **逐位相同（101.61=101.61）**。

**源码证据**：`runtime_dat_mod` 只在解析时"识别" DIRECT（`tok(1) == "MCM" .or. tok(1) == "DIRECT"`），
`eval_J` **从不读 mode**，两种配置走同一条公式。

**影响**：设 `DIRECT` 期望"不随天顶角变化"，实际得到 MCM 参数化——静默无效果。

**修复**（方案 A）：新增 `phot_mode` 变量，`eval_J` 分支：

```fortran
if (phot_mode == 1) then
    JARR(i) = jf*j_l(i)*cosz**j_m(i)*exp(-j_n(i)/cosz)   ! MCM
else
    JARR(i) = jf*j_l(i)                                   ! DIRECT: 不调制
end if
```

**验证**：MCM 101.61 / DIRECT **102.92（差 1.29%）** ✓ 开关生效；
MCM 结果**零回归**（101.61 与修前逐位相同）。

**新增防线**：`R16 MODE=DIRECT 与 MCM 产生不同结果 (非死开关)`

### 14c.2 🟢 P3：`RADIAT` 死配置

`meteo.config` 第 4 量 `RADIAT` 被解析但无任何消费者（已在 `UNITS.md §8` 登记）。
**修复**：`cfgbuild` 在 `RADIAT ≠ 0` 时打印明确警告。**双向验证**：默认（0）无警告、
设 500 时警告出现。

### 14c.3 ✅ 九项健康确认

| # | 项 | 结果 |
|---|---|---|
| 1 | 完整 `make test` 端到端 | 冒烟 1/1 → 24/24 → 131/131 → 81 案例（失败 0）→ 诊断表，**一次通过** |
| 2 | **速率常数 vs JPL/IUPAC** | HO₂+NO **1.01×**、O₃+NO **0.96×**、CH₄+OH **0.99×**、CO+OH **0.99×**、HCHO+OH 1.18× — **四条核心全在 ±5%** |
| 3 | **O₃ 收支闭合** | 末 24h 六账净残差 **0.34%**（P=4.63e13, L=−4.16e13）；化学净占 P 的 10.2% |
| 4 | **`ng` 边界 1 与 32** | 1 案例 0.35s；**1024 案例 865s，输出 1025 行无越界** |
| 5 | **jfac/ΔT 网格维** | jfac 0.6→1.0 升 19.9%、ΔT +5K 升 9.4% — **两维均真实生效、方向正确** |
| 6 | **mcm_run 缩放** | BAU 16.22 / VOC×3 16.37 / NOx×0.3 14.66 — 生效 |
| 7 | **SCEN_HOT 温度通路** | 101.61 → **112.56（+10.8%）** — `dt_scen` 生效 |
| 8 | **全配置项消费扫描（14 项）** | 12 项真实消费，2 项死配置（见上） |
| 9 | 链接与产物 | 46 处链接断链 0；clean 后产物 0 |

### 14c.4 R16 收官判定

**"提前收官"不完全成立**——存在 1 个静默无效果的配置开关（MODE=DIRECT）。
性质温和：默认值 `MCM` 是正确路径，`DIRECT` 属未实现的预留能力，**不影响任何历史结论**。

**修复后判定：收官成立。**

### 14c.5 第八条方法学教训

**⑧ "配置项存在 ≠ 配置项生效"**
前 15 轮验证了"输入能进管线"，但没验证"每个配置项真的被消费"。
本轮的系统扫描（14 项逐项查消费点）才暴露 `MODE` 与 `RADIAT`。

**对策**：凡 `runtime.dat` 里写入的字段，必须在 Fortran 侧有对应读取点，否则
**(a) 实现它 / (b) 删除它 / (c) 显式警告——三者必居其一，不允许静默存在。**

---

<a name="15-永久防线总表12-条"></a>
## 15. 永久防线总表（12 条）

每一轮的修复都留下至少一条机器检查。当前 `make test` 覆盖：

| # | 防线 | 位置 | 抓什么问题 |
|---|---|---|---|
| 1 | **FT0 生成码一致性（ns + nr 双守卫）** | functional_test 开头 | 改了机制/生成码未重编（假绿） |
| 2 | **BC 行 token 数断言** | runtime_dat_mod | 行缓冲截断（静默读垃圾） |
| 3 | **SERIES 上限断言**（序列数 ≤8 / 点数 ≤128） | runtime_dat_mod | 手改 runtime.dat 越界 |
| 4 | **RUN.ns 一致性断言** | runtime_dat_mod | 机制/生成码/runtime 三者不同步 |
| 5 | **sp_name 硬失败**（缺 SPECIES 行） | runtime_dat_mod | 物种名退化为 SPi → 组映射全失效 |
| 6 | **运行期氮守恒不变量**（Σ w_N·dy/dt ≈ 0） | fac_batch_omp | NOy 物种漏权（新 N 物种） |
| 7 | **NOy 数组容量断言**（noy_id(32)） | fac_batch_omp | NOy 物种数超限越界 |
| 8 | **关键物种索引断言**（O3/HCHO/NO2/H2O2/HNO3/NO） | fac_batch_omp | `yp(0)` 越界读 |
| 9 | **batch nfail → stop 1** | fac_batch_omp | 部分格点失败仍报绿 |
| 10 | **元素审计（方向限定）** | audit/element_audit.py | 非物理配平（用储库掩盖真错） |
| 11 | **4 金值哨兵**（屋脊/FNR/OPE 等） | functional_test | 机制/求解改动导致数值漂移 |
| 12 | **IPR 闭合 + IPR↔IRR 双路一致** | functional_test + hourly_pa | 过程账错位/双计 |

**R10–R13 新增的语义断言**（无需金值的"路径正确性"检查）：

| # | 断言 | 抓什么问题 |
|---|---|---|
| A | **R10-P0 HVOC 缩放全部 E_VOC 组物种 (×3)** | 情景组映射漏物种（scen_build 双份清单） |
| B | **R11-P0 jfac 等价路径**（cfg0.8×CLI1.0 == cfg1.0×CLI0.8） | 乘子型参数重复施加 |
| C | **R10 batch_ipr/irr 表头 == 数据列数** | 尾逗号空列（读取器整体错位） |
| D | **R12 batch 输出名 == figures 依赖名** | 构建系统数据链断裂 |
| E | **R13 输出格式无三位无E指数** | Fortran 格式静默丢 E |
| F | **R13 csvfloat 四种指数形态解析** | 解析器不健壮 |
| G | **sat2surf 两法互证 + 旧式反证** | 单位换算差 9 量级 |

---

<a name="16-修复统计与方法学沉淀"></a>
## 16. 修复统计与方法学沉淀

### 16.1 修复统计

| 类别 | 数量 |
|---|---|
| 真 bug 修复 | **32 个**（P0×4 / P1×8 / P2×14 / P3×6，含 R15 的 4 个） |
| 排除的误报 | **3 项**（`dsum(15)` 定长、`indicators` NOy 自动覆盖、`mechc_c` 重复） |
| 新增机器防线 | **12 条常设 + 7 条语义断言** |
| 文档漂移修正 | **40+ 处**（计数口径/失效 API/失效路径） |
| 垃圾清理 | **6 类**（`__pycache__`/`.DS_Store`/`.bak`/陈旧二进制/陈旧副本/陈旧产物） |
| 关键数值验证 | 101.61 / 111.26 / 125.54 / 94.91 等 10+ 项与文档逐项对齐 |

### 16.2 七条核心教训

**① 默认值掩盖 bug**
`JFAC=1.0` 让 jfac 双乘退化为单乘——13 轮测试、40+ 次实验、4 个哨兵全部无法发现。
**对策**：对每个"乘子型参数"都测非 1 值的等价路径。

**② 双胞胎结论矛盾 = 实现差异的强信号**
P0-1（scen_build）之所以被发现，是因为 fac_v3 与 batch 对同一 HVOC 情景给出
**−0.2% vs +8.0%**。我当时用"网格最近点差异"解释了它——**这个解释掩盖了真 bug 两个版本**。
**对策**：当两条共享机制的路线给出相反方向时，**先 diff 实现，不要 diff 参数**。

**③ 越基础的层，出错越隐蔽**
`batch_v3.csv` vs `batch_v3_indicators.csv` 两个近义名让数据链断了多轮无人发现
（因为没人验证"新跑的数据是否真的进了图"）。
**对策**：文件名即契约，跨工具数据流只用一个 canonical 名。

**④ 格式宽度是静默契约**
`ES13.5` 对常规值完美，只在 |指数|≥100 时静默改语义。
**对策**：这类 bug **无法通过正常工况测试发现，必须专门测极端值**。

**⑤ 重复实现 = 假保护**
`_fix_fortran_exp` 复制在 6 个文件里，让人误以为"已经保护了"，实际 7 个文件裸奔。
**对策**：同一张表/同一段逻辑出现在两个文件里 = 漂移与遗漏风险，应单源化。

**⑥ 文档里的命令必须被执行过**
手册 §7 示例引用已删除的 `mechanism/ozone.fac`，照抄必崩；README 宣称的 6 个诊断函数
已删除数月无测试发现。
**对策**：文档声称的"现行"必须有 CI 覆盖，否则标"历史"。

**⑦ 守恒不变量是最高价值的探测器**
氮守恒抓出的都是"机制科学性"层面的错误（NOy 漏权），而非数值噪声——
它能抓出"静默破坏物理守恒"这类最隐蔽的 bug。

### 16.3 最终状态

```
3 条实路线     fac_v3（单箱深度）· fac_batch_omp（批量广度）· mcm_run（MCM 对照）
74 物种        131 反应 · 41 光解槽位 · 105 命名速率
元素审计       131/131 PASS（储库方向限定）
测试链         24/24 PASS（12 条常设防线 + 7 条语义断言 + 4 金值哨兵）
文档           17 份 / 8723 行（含 4 份分轮审计报告 + 本合并记录）
可复现性       从 make clean 可完全重建，结果与主线一致
```

**十六轮自纠自查到此收官。验收成功。**

---

## 附录：分轮审计报告索引

| 文件 | 内容 |
|---|---|
| `audit_reports/AUDIT_R10.md` (319 行) | scen_build 漏缩 / IPR 双计 / 尾逗号 / 边界保护 |
| `audit_reports/AUDIT_R11.md` (271 行) | jfac 双乘 / scan_parallel 死代码 / nfail / 窗口样本数 |
| `audit_reports/AUDIT_R12.md` (222 行) | 数据链断裂 / 文档 API / make nc / clean / input 定位 |
| `audit_reports/AUDIT_R13.md` (180 行) | 三位指数 / sat2surf / fig_ridge5 / sys.path |
| **本文件** `AUDIT_COMPLETE.md` | **R1–R14 完整合并记录（含 R1–R9 的重建）；4 份分轮报告在 `audit_reports/`** |
