# tests/ —— 回归防线（本项目的"质检科"）

> **一句话**：两个文件构成三层防线——**冒烟**（能不能跑）、**契约**（输出格式对不对）、
> **金值哨兵**（数值漂没漂）。`make test` 全绿是任何实验的前置条件。

---

## 1. 两个文件

| 文件 | 断言数 | 耗时 | 入口 |
|---|---|---|---|
| `smoke_test.sh` | 1 | ~1 s | `make smoke` |
| `functional_test.py` | **59 项断言**（源码 52 处 `check`，循环展开后实测 59/59） | ~3 min | `make functional` / `make test` |

---

## 2. `smoke_test.sh` —— 冒烟（1 分钟内知道"死了没"）

只做一件事：跑 `fac_v3 SCEN_BAU`，断言 stdout 含 `ok` 且 O₃ 末列非零。
**Makefile 依赖**：`smoke: all runtime_v3.dat`（runtime 是生成物已 gitignore，
不显式依赖会让全新克隆首次 `make smoke` 直接报"打不开 runtime_v3.dat"）。

## 3. `functional_test.py` —— 功能测试（11 族断言，实测 59 项）

每族针对一类**真实踩过的 bug**，断言数随防线增加：

| FT 族 | 断言 | 防的历史 bug |
|---|---|---|
| `FT-CFG` (2) | 编译期负向：未知物种拒绝 / 气象序列拒绝浓度单位指令 | 未知物种曾运行期才炸；`# unit: ppb` 写进气象文件会静默 ×2.46e10 |
| `FT-MP` (2) | `model.parameters` 拒绝 dtout=0 / days=0 / 缺行 / 非数值 | 旧版零校验，除零后照样打印 ok 并出空 CSV |
| `FT-DAYS` (2) | `--days` 覆盖生效且**不改** model.parameters 文件 | 天数覆盖与文件配置混淆 |
| `FT-PHOT` (2) | 光解 LAT/DOY 真正生效（夏至 vs 冬至 O₃ 差 >20 ppb） | DOY 解析判据 `nt>=4` 永假，换季节静默无效 |
| `FT-SAT` (2) | sat2surf 卫星转换端到端 + 两法互证 | 旧公式漏 ×1e9，差 9 个量级 |
| `FT-SYNC` (4) | 扩机制 8 处同步：无死端 / NOy 表与 indicators 一致 / 排放物种已分组 / ATOL 覆盖度 | 扩物种漏配 → 链条空转或情景失效 |
| `FT-CFGSPEC` (2) | input_v3 无失效物种名 / 无"只产不耗无物理汇"死端 | GLY 等稳定产物漏干沉（实测高估 32–49%） |
| `FT-FREEZE` (8) | 冻结版规格：.fac 计数=声明、README 同步、V25_PLAN/teaching 已删、mcm_city 有定位文档 | 冻结后文档漂移 |
| `FT-GIT` (1) | **output/ 无任何文件被 git 跟踪**且已被 ignore | 防实验产物误入库 |
| `FT-TPL` (1) | templates 全部"原样可用"（无缺 `#` 的说明行） | 4 个模板第 2 行缺 # → 用户照抄即崩 |
| `FT-UNIT` (1) | IC 三态单位（cm3/ppb/ppm）自动换算精确 | 单位换算回归 |
| （核心契约） | batch 24 列×81 案例 / IPR-IRR 表头列数一致 / 屋脊回撤≥20% / FNR 单调 / IPR 六账闭合 / L_NOx-Q 值域 / **4 金值哨兵** | OPE 双计、尾逗号空列、jfac 双乘、scen_build 漏缩 |

### 金值哨兵（抓数值漂移的锚）

```python
SENT = [  # (fnox, CSV列, 期望值, 容差)
    ("2.28", 4, 138.39, 0.10),   # O3max 屋脊上升段
    ("3.42", 4, 123.91, 0.10),   # O3max 峰值
    ("1.01", 8, 4.988, 0.15),    # FNR 转折带
    ("2.28", 16, 4.234, 0.15),   # OPE 饱和区谷
]
```
机制/物理/光解任一改动都会动这些数 → **必须重标**并注明沿革
（V2.3 → V2.4+MCM 对齐 → V2.5 三次沿革都记录在 SENT 注释里）。

### 运行特性

- 自带 `/tmp/ft_*` 中间产物**自清理**（跑完不留垃圾）
- 全程只用已编译二进制 + Python 标准库 + numpy，**不触发编译**
- 退出码 0/1，`make test` 依赖

## 4. 什么时候改这个目录

| 场景 | 改哪 |
|---|---|
| 机制/物理/光解改动 → 数值变 | **重标 SENT 四值**（在注释里记沿革） |
| 新增防线（踩了新坑） | 加一个 `check(...)`，归入相应 FT 族或开新族 |
| 新增输出契约（新 CSV 结构） | 加列数/表头一致性断言 |
| ⚠ 红线 | **不许为让测试变绿而放宽断言**——先查根因，确认是"预期变化"才改期望值并留注释 |

## 5. 与质检体系的关系

```
audit/element_audit.py   元素守恒（机制静态）      ┐
audit/mcm_compare.py     系数 vs MCM（机制静态）    ├─ make test
tests/smoke_test.sh      能跑（运行时冒烟）         │
tests/functional_test.py 契约+哨兵（运行时回归）    ┘
experiments/run/run_D    六账闭合（科学性，实验层）
```
