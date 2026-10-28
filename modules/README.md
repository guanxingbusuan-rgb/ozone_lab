# modules/ —— Fortran 共享模块层（编译期依赖，不可删）

> **一句话**：两个 `.f90` 被**静态编入** `fac_v3` 与 `fac_batch_omp`；
> 运行期不读取本目录（二进制里已有）。`modules/` 只在 `make all` 时被编译。

---

## 1. 两个文件（均在役）

| 文件 | 进入二进制 | 作用 |
|---|---|---|
| `runtime_dat_mod.f90` | ✅ fac_v3 + fac_batch_omp | `runtime.dat` 的唯一消费端：SPECIES/IC/MET/BC/PHOT/EMIS/DRY/WET/DIL/SCEN/SERIES 解析；气象序列插值（≤128 点/序列，池上限 8）；SZA→J 参数化；情景组映射 `build_groups` |
| `chem_phys_mod.f90` | ✅ fac_v3 + fac_batch_omp | 物理四账的唯一实现：① 排放（BOX/FLAT 日型 × 情景缩放）② 干+湿沉降 ③ BC 张弛 + 稀释 ④ BLH 膨胀稀释；另含按物种名的 ATOL 调优表与 `eval_J_jf`（CLI jfac × 配置 JFAC 的归一化） |

Makefile 依赖（可验证）：

```
facsimile_impl/fac_v3:        ... modules/runtime_dat_mod.f90 modules/chem_phys_mod.f90
facsimile_impl/fac_batch_omp: ... modules/runtime_dat_mod.f90 modules/chem_phys_mod.f90
```

---

## 2. 关键不变量

1. **单一物理实现**：`fac_v3`（单箱）与 `fac_batch_omp`（批量）共用同一份物理层，
   两条路线的差异只在驱动方式，不在物理。
2. **单源组映射**：`build_groups` 是排放组（NOx / VOC / CO / CH4）的唯一事实源；
   VOC 组含 SO₂/DMS/MSA（工程既有约定，两个 EKMA 的 VNr 因此区分 group / nmv 两列）。
3. **序列硬上限**：`NSER_MAX=8`、每序列 128 点；`cfgbuild` 编译期已拦截，运行时再校验一次。
4. **光解口径**：`eval_J_jf(t, jf_abs, …)` 以“绝对因子”为输入，内部按配置侧 JFAC
   归一化，避免 CLI 与 config 双重相乘。

---

## 3. 修改与重编译

改任意 `.f90` 后必须重新构建，才会进入二进制：

```bash
make all            # 或 make clean-build && make all
bash check_env.sh   # 13 项自检（二进制存在/动态链接/工具链）
```

- 只改各程序的 `config_override/` 或 Python 代码 → **不需要编译**；
- `.mod` 接口文件生成在 `facsimile_impl/gen_v3/`，属编译期产物，可随 `make clean-build` 删除重建。
