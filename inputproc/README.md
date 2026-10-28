# inputproc/ —— 数据编译前端（Python 侧的"海关"）

> **一句话**：所有**外部数据**进入引擎前的最后一站——校验、换算、归一，
> 然后压成 Fortran 能直接读的平面文件。**两个工具都编译无关**（纯 Python）。

---

## 1. 两个文件

| 文件 | 行数~ | 定位 | 调用方 |
|---|---|---|---|
| `cfgbuild.py` | 500 | **配置编译器**：12 类 config + constraints → `runtime_v3.dat` | `make runtime_v3.dat`、`experiments/run/_runner.py`（F/H/M 实验变体）、`make run DAYS=` |
| `sat2surf.py` | 300 | **卫星柱 → 地表 ppb 转换器**（离线工具） | 手动运行；`FT-SAT` 测试 |

## 2. `cfgbuild.py` —— 配置编译器（数据链的核心枢纽）

```
configs/ (12 类人类可读配置 + constraints/ 序列)
    ↓  python cfgbuild.py configs runtime_v3.dat [--mech ...] [--days N] [--dtout S]
runtime_v3.dat   ← Fortran 顺序表控平面文件（引擎唯一输入，勿手改）
```

### 编译期完成的事（全过才写出，任何一项失败拒绝编译）

| 校验/换算 | 细节 |
|---|---|
| 物种存在性 | IC/BC/EMIS/DRY/WET 里出现的物种必须在机制 VARIABLE 段 |
| 浓度单位三态 | `cm3`(×1) / `ppb`(×2.46e10) / `ppm`(×2.46e13)，`unit_factor()` 单源；未知单位拒绝 |
| 排放第四单位 | 文件头 `# unit: ppb/h` → ×6.833e6 |
| 序列质量 | 时间严格递增、≤128 点、气象+BC 共享 8 条池、**气象序列禁止 `# unit:` 行**（防静默 ×2.46e10） |
| 气象 | BLH 序列必须 >0（防除零）；RADIAT 非 0 告警（无消费者） |
| 光解完备性 | 机制用到的每个 `J<n>` 槽位必须有非零 l（表值或 JMAP 直给），否则拒绝——防静默零光解 |
| model.parameters | 4 行必填 + 第 5 行可选 SPINUP（0≤SPINUP<days）；行数/数值/days>0/dtout>0 全查；dtout≠3600 与 nsteps<24 警告 |
| 天数覆盖 | `--days N` 只改本次编译产物，**不碰源文件**（`DAYS=` 同源） |
| JMAP | 只用于光解表缺失的槽位（语法 `槽=D正午值`） |

**设计原则**：编译期把一切人能犯的错拦下，Fortran 端**零文本解析、零校验**——
所以"只改数据不重编二进制"（`runtime_v3.dat` 是纯数据）。

### FORCE 机制

`make runtime_v3.dat` 挂 `FORCE` 每次强制重建（仅 0.4 s）——
防"同一秒内改配置又改回，mtime 未变 → make 误判 up-to-date → 静默跑旧数据"（实测踩过）。

## 3. `sat2surf.py` —— 卫星柱 → 地表约束转换器

**为什么存在**：ozone_lab 是 0-D 箱体，没有气柱概念，引擎**不直读卫星柱浓度**
（`docs/templates/satellite_column.template` 明确此边界）。TROPOMI 等卫星数据的正确用法：

```
原始柱浓度 [molec/cm²] ──sat2surf──> 地表等效 ppb 序列 ──> bc.config 张弛约束
```

```bash
python inputproc/sat2surf.py 原始文件 输出文件 --blh 1000        # 单一 BLH
python inputproc/sat2surf.py 原始文件 输出文件 --blh-file xx.txt # BLH 序列插值
python inputproc/sat2surf.py --selftest                         # 内置自检
```

| 特性 | 说明 |
|---|---|
| 主准则 | **直观法**：`柱 / BLH[cm] / 2.46e10`（假设最少） |
| 参考法 | 混合比法（空气柱总量口径），两法强制互证，分歧 >阈值则警告但按直观法继续 |
| QA 过滤 | 卫星 qa 值 <0.75 的过境剔除 |
| 抽稀 | >128 点自动降采样（Fortran 序列上限） |
| 物种无关 | HCHO/NO₂/任一机制物种皆可（无白名单） |

⚠ 历史教训：旧式公式 `柱/(2.15e25×混合层分数)` **漏 ×1e9**（差 9 个量级），
已修正并由 `FT-SAT` 端到端用例固化。

## 4. 不在这里的 Python

| 工具 | 在哪 | 为什么 |
|---|---|---|
| 机制编译器 `mechc.py` | `facsimile_impl/` | 它的产物（gen_v3）与二进制同目录，构建链内聚 |
| 诊断/指标 | `analysis_lib/` | 仿真**之后**的数值消费（inputproc 是仿真**之前**的数据入口） |
| 机制审计 | `mechanism/audit/` | 质检职责 |
| 实验编排 | `experiments/run/` | 业务层 |

**分工口诀**：`inputproc/` 管"数据怎么说"（进引擎前），`facsimile_impl+modules` 管"方程怎么算"，
`analysis_lib/` 管"结果怎么读"（出引擎后）。

## 5. 改动风险

- `unit_factor()` 是三处共用的**换算单源**——改它等于改所有浓度输入的语义（`FT-UNIT` 会红）。
- `read_series(allow_unit=False)` 的气象保护是防 10 个量级静默错误的最后闸门（`FT-CFG` 会红）。
- `cfgbuild` 的所有拒绝/警告文案都有测试断言关键字——改文案先查 `functional_test.py`。
