# input/ —— 15sp 教学遗产（**非工程入口**）

> ⚠ **本目录不是工程主线**。工程配置在 [`../input_v3/`](../input_v3/)、机制在
> [`../mechanism/`](../mechanism/)、数据契约见 [`../UNITS.md`](../UNITS.md)。

## 这是什么

Day 35 时期的**教学版配置集**：15 物种 / 25 反应的核心机制 + 对应的 12 类配置。

| 文件 | 说明 |
|---|---|
| `mechanism_1_ozone_tutorial.fac` | **15sp/25rx 教学机制**（三大循环骨架，无 NO₃/PAN/ISOP） |
| `mechanism_2_mcm_ch4.fac` | MCM v3.3.1 CH₄ 子集（29sp/71rx）——已由 `mcm_city/input/` 接管 |
| 其余 10 个 `.config` | 15sp 对应的配置（排放/气象/BC/沉降/情景等） |

## 为什么保留

1. **知识谱系**：[`mechanism/README.md`](../mechanism/README.md) 的版本沿革（V2.1 起点）
   引用它作为"最小可运行骨架"的参照
2. **教学对照**：讲三大循环时，15sp 版本比 74sp 更适合逐反应板书

## 为什么不要用

`cfgbuild` 的默认输入目录是 **`input_v3`**（R12 修正）；`Makefile` 全部目标也只走
`input_v3`。用本目录编译会得到 15sp 的 `runtime.dat`，与 74sp 生成码不匹配
（`RUN.ns` 断言会在运行时报错）。

```bash
# ✅ 正确
make runtime_v3.dat                    # → input_v3 → runtime_v3.dat

# ❌ 错误（会产出与生成码不匹配的 runtime）
python src/cfgbuild.py input runtime.dat --mech input/mechanism_1_ozone_tutorial.fac
```
