# ============================================================
# ozone_lab 根 Makefile —— 三实现统一构建 + 测试一键盘
# ============================================================
# ── 工具链路径 (可用环境变量覆盖, 无需改文件) ──
# 用法示例:
#   make all PY=/usr/bin/python3
#   make all FC=gfortran-13 CC=gcc-13 CVODELIB=/usr/local/lib
#   make all CONDA_ENV=$HOME/miniconda3/envs/ozone
CONDA_ENV ?= $(HOME)/miniconda3/envs/venv
GCC_ROOT  ?= $(HOME)/gcc-13.2.0
PY ?= $(CONDA_ENV)/bin/python
# 注意: FC/CC 是 make 内置变量 (默认 f77/cc), 用 ?= 不会生效,
#       必须用 origin 判断: 命令行/环境覆盖 > 我们的默认 > make 内置
ifeq ($(origin FC), default)
  FC := $(GCC_ROOT)/bin/gfortran
endif
ifeq ($(origin CC), default)
  CC := $(GCC_ROOT)/bin/gcc
endif
CVODELIB ?= $(HOME)/soft/atmenv/cvode/lib
FFLAGS := -Wall -O2 -ffree-line-length-none
export PATH := $(CONDA_ENV)/bin:$(PATH)

.PHONY: all mcm mcmcheck smoke functional test v3 batch run figures nc clean clean-build clean-output FORCE

all: facsimile_impl/fac_v3 facsimile_impl/fac_batch_omp

# ---------- 工程版 V3 (mechanism_3 全链条 + runtime.dat 通路) ----------
FACV3 := facsimile_impl/fac_v3 facsimile_impl/fac_batch_omp

v3: $(FACV3)

facsimile_impl/gen_v3/oz_mech.f90: mechanism/mechanism_3_full_chain.fac facsimile_impl/mechc.py
	cd facsimile_impl && $(PY) mechc.py ../mechanism/mechanism_3_full_chain.fac gen_v3

# FORCE: 每次 make 都重建 runtime_v3.dat (cfgbuild 仅 0.37 s)。
# 原因: 仅靠 mtime 会漏掉"同一秒内改配置又改回"的情形 —— make 误判 up-to-date,
# 导致 runtime.dat 保留中间值, 模拟结果无声出错(实测复现)。
# cfgbuild 幂等 + 极快, 强制重建是更安全的默认。
#
# 模拟天数: 用 DAYS= 临时覆盖, 不改 input_v3/model.parameters 文件
#   make runtime_v3.dat DAYS=3      # 本次只跑 3 天
#   make run DAYS=1 SCEN=SCEN_BAU   # 编译+跑一步到位
DAYS ?=
_CFGBUILD = $(PY) src/cfgbuild.py input_v3 runtime_v3.dat \
            --mech mechanism/mechanism_3_full_chain.fac \
            $(if $(DAYS),--days $(DAYS),)

runtime_v3.dat: FORCE $(wildcard input_v3/* input_v3/constraints/*) src/cfgbuild.py
	@$(_CFGBUILD)

FORCE:

facsimile_impl/fac_v3: facsimile_impl/oz_fac_v3.f90 facsimile_impl/gen_v3/oz_mech.f90 \
    modules/runtime_dat_mod.f90 modules/chem_phys_mod.f90
	cd facsimile_impl && $(FC) $(FFLAGS) -J gen_v3 -I gen_v3 gen_v3/oz_mech.f90 \
	  ../modules/runtime_dat_mod.f90 ../modules/chem_phys_mod.f90 oz_fac_v3.f90 -o fac_v3 \
	  -L$(CVODELIB) -Wl,-rpath,$(CVODELIB) \
	  -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial

facsimile_impl/fac_batch_omp: facsimile_impl/fac_batch_omp.f90 facsimile_impl/cvode_c_mod.f90 \
    facsimile_impl/gen_v3/oz_mech.f90 modules/runtime_dat_mod.f90 modules/chem_phys_mod.f90
	cd facsimile_impl && $(FC) $(FFLAGS) -fopenmp -J gen_v3 -I gen_v3 gen_v3/oz_mech.f90 \
	  ../modules/runtime_dat_mod.f90 ../modules/chem_phys_mod.f90 cvode_c_mod.f90 fac_batch_omp.f90 -o fac_batch_omp \
	  -L$(CVODELIB) -Wl,-rpath,$(CVODELIB) -lsundials_cvode -lsundials_nvecserial

# nc 需要 <prefix>_conc.csv + <prefix>_irr.csv 配对 (fac_v3 输出格式);
# batch 的 x_indicators.csv 是网格汇总, 不是时间序列 → 不可用于 nc
NC_SRC ?= output/expA_BAU
NC_OUT ?= output/expA_BAU.nc

nc:
	$(PY) modules/nc_output.py $(NC_SRC) $(NC_OUT)

# 注: figures/ 整个目录不入库（.gitignore）。新克隆后该目标会提示自备脚本。
figures: output/batch_v3_indicators.csv
	@test -f figures/fig_indicators_ridge.py || { \
	  echo "[figures] 缺少 figures/ 目录（该目录不入库，仅本地保留）"; \
	  echo "          需要本地自备: fig_indicators_ridge.py / fig_response_v3.py /"; \
	  echo "                        gen_all.py / gen_cartoon.py"; \
	  exit 1; }
	cd figures && $(PY) fig_indicators_ridge.py ../output/batch_v3_indicators.csv indicators_vs_ridge.png
	cd figures && $(PY) fig_response_v3.py ../output/batch_v3_indicators.csv response_surface_v3.png

# ---------- MCM 对照路线 (29sp CH4 子集, 独立子树) ----------
mcm_city/gen/oz_mech.f90: mcm_city/input/mechanism_mcm_ch4.fac facsimile_impl/mechc.py
	cd mcm_city && $(PY) ../facsimile_impl/mechc.py input/mechanism_mcm_ch4.fac gen

mcm_city/mcm_run: mcm_city/mcm_driver.f90 mcm_city/gen/oz_mech.f90
	cd mcm_city && $(FC) $(FFLAGS) -J gen -I gen gen/oz_mech.f90 mcm_driver.f90 \
	  -o mcm_run -L$(CVODELIB) -Wl,-rpath,$(CVODELIB) \
	  -lsundials_fcvode -lsundials_cvode -lsundials_fnvecserial -lsundials_nvecserial

mcm: mcm_city/mcm_run

batch: v3 runtime_v3.dat
	cd facsimile_impl && ./fac_batch_omp 9 ../output/batch_v3_indicators.csv ../runtime_v3.dat

# 单情景一步跑: make run [DAYS=n] [SCEN=SCEN_X] [PREF=前缀]
SCEN ?= SCEN_BAU
PREF ?= ../output/run1
run: v3 runtime_v3.dat
	@echo "[run] $(SCEN) · $(if $(DAYS),$(DAYS) 天,默认天数) → $(PREF)"
	cd facsimile_impl && ./fac_v3 $(SCEN) $(PREF) ../runtime_v3.dat

# 注: make all 不含 mcm_run（它是独立子树），需要时显式 make mcm



# ---------- 测试 ----------
# runtime_v3.dat 是生成物（已 gitignore），全新克隆里不存在 →
# smoke 必须显式依赖它，否则 fac_v3 报 "Cannot open ../runtime_v3.dat"
smoke: all runtime_v3.dat
	sh tests/smoke_test.sh

functional: all smoke
	$(PY) tests/functional_test.py

test: functional
	$(PY) audit/element_audit.py mechanism/mechanism_3_full_chain.fac
	$(PY) audit/mcm_compare.py --quiet
	$(PY) scan_suite.py

# 仅比对 MCM 速率系数（快速，不跑仿真）
mcmcheck:
	$(PY) audit/mcm_compare.py

# ── 安全清理变体（做实验期间用这两个，别用 make clean）──
# clean-build: 只删编译产物与生成码，保留 output/ 里所有实验结果
clean-build:
	cd facsimile_impl && rm -f *.o *.mod fac_v3 fac_batch_omp \
	  gen/oz_mech.* gen_v3/oz_mech.* gen_v3/oz.reac gen_v3/oz.prod gen_v3/*.mod
	rm -f modules/rad_photol_demo modules/*.mod audit/*.dylib
	cd mcm_city && rm -f *.o *.mod mcm_run gen/*.mod gen/oz_mech.f90 \
	  gen/oz.reac gen/oz.prod
	find . -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
	@echo "[clean-build] 已删二进制与生成码；output/ 中的实验结果保留"

# clean-output: 只删实验结果（保留二进制），需要显式确认
clean-output:
	@printf "将删除 output/ 下所有实验产物，确认? [y/N] " && read ans && [ "$$ans" = "y" ] && \
	  rm -rf output && echo "[clean-output] 已删除 output/" || echo "[clean-output] 已取消"

clean:
	cd facsimile_impl && rm -f *.o *.mod fac_v3 fac_batch_omp \
	  gen/oz_mech.* gen_v3/oz_mech.* gen_v3/oz.reac gen_v3/oz.prod gen_v3/*.mod
	rm -f modules/rad_photol_demo modules/*.mod audit/*.dylib
	cd mcm_city && rm -f *.o *.mod mcm_run gen/*.mod gen/oz_mech.f90 \
	  gen/oz.reac gen/oz.prod mcm_city_*.csv
	rm -f *.nc output/*.nc
	rm -rf output __pycache__ *.png
	find . -name __pycache__ -type d -not -path "./output/*" -exec rm -rf {} + 2>/dev/null || true
