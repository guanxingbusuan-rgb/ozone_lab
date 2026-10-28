# ============================================================
# ozone_lab —— 程序级精简化后的构建入口（2026-10-29）
# 三个自完备程序: CASE-analysis / EKMA-FNR / EKMA-emission
# 共享层: mechanism/ facsimile_impl/ modules/ inputproc/ analysis_lib/
# 已退役: 通用实验层 (experiments/ configs/ output/ runtime_v3.dat)
# 各程序自带 config_override/ 与 output/，入口是各自的 run.py / plot.py。
# ============================================================
CONDA_ENV ?= $(HOME)/miniconda3/envs/venv
GCC_ROOT  ?= $(HOME)/gcc-13.2.0
PY ?= $(CONDA_ENV)/bin/python
ifeq ($(origin FC), default)
  FC := $(GCC_ROOT)/bin/gfortran
endif
ifeq ($(origin CC), default)
  CC := $(GCC_ROOT)/bin/gcc
endif
CVODELIB ?= $(HOME)/soft/atmenv/cvode/lib
FFLAGS := -Wall -O2 -ffree-line-length-none
export PATH := $(CONDA_ENV)/bin:$(PATH)

.PHONY: all mcmcheck jcheck clean-build

all: facsimile_impl/fac_v3 facsimile_impl/fac_batch_omp

facsimile_impl/gen_v3/oz_mech.f90: mechanism/mechanism_3_full_chain.fac facsimile_impl/mechc.py
	cd facsimile_impl && $(PY) mechc.py ../mechanism/mechanism_3_full_chain.fac gen_v3

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

# 机制静态审计（不跑仿真；不依赖已退役的 configs/）
mcmcheck:
	$(PY) mechanism/audit/mcm_compare.py
jcheck:
	$(PY) mechanism/audit/jslot_audit.py

clean-build:
	cd facsimile_impl && $(RM) *.o *.mod fac_v3 fac_batch_omp \
	  gen/oz_mech.* gen_v3/oz_mech.* gen_v3/oz.reac gen_v3/oz.prod gen_v3/*.mod
	$(RM) modules/*.mod
