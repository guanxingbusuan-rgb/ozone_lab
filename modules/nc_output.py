#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
B6: netCDF 输出 — CSV → netCDF4 (Day 20 栈消费)
  将三实现输出的 <prefix>_conc.csv + <prefix>_irr.csv 打包为 .nc
用法: python nc_output.py <prefix> <out.nc>
"""
import sys, os, csv
import numpy as np
from netCDF4 import Dataset


def _fix_fortran_exp(tok):
    """gfortran 3 位指数省略 E 的遗留格式 → 标准科学计数"""
    import re as _re
    if "E" not in tok and "e" not in tok:
        m = _re.search(r"(\d)-(\d{3})$", tok)
        if m:
            return tok[:m.start(1) + 1] + "E" + tok[m.start(1) + 1:]
    return tok

def csv_to_nc(conc_csv, irr_csv, out_nc, impl_name="FACSIMILE"):
    """CSV 时间序列 → netCDF4"""
    # 读浓度
    rows = list(csv.reader(open(conc_csv)))
    hdr = [h.strip() for h in rows[0]]
    data = np.array([[float(_fix_fortran_exp(x)) for x in r] for r in rows[1:]])
    t_h = data[:, 0]  # 小时
    conc = data[:, 1:]  # 浓度 (T, nsp)

    # 读 IRR
    irr_rows = list(csv.reader(open(irr_csv)))
    irr_data = np.array([[float(_fix_fortran_exp(x)) for x in r] for r in irr_rows[1:]])
    irr = irr_data[:, 1:]  # (T, nr)

    # 写 nc (E14 约定: 声明 (spid, timid) + 转置写入)
    nc = Dataset(out_nc, "w", format="NETCDF4")
    nc.setncattr("implementation", impl_name)
    nc.setncattr("species_names", ",".join(hdr[1:]))
    nc.createDimension("time", len(t_h))
    nc.createDimension("species", conc.shape[1])
    nc.createDimension("reaction", irr.shape[1])

    v = nc.createVariable("time", "f8", ("time",))
    v.units = "hours"
    v[:] = t_h

    # 浓度: 声明 (species, time) + 转置写入 → Python 读 (time, species)
    v = nc.createVariable("conc", "f8", ("species", "time"))
    v.long_name = "concentration"
    v.units = "molec/cm3"
    v[:] = conc.T

    # IRR
    v = nc.createVariable("irr", "f8", ("reaction", "time"))
    v.long_name = "instantaneous reaction rates"
    v.units = "molec/cm3/s"
    v[:] = irr.T

    nc.close()
    print(f"[nc_output] {out_nc}: {conc.shape[0]} t × {conc.shape[1]} sp × {irr.shape[1]} rxn")

if __name__ == "__main__":
    pref = sys.argv[1] if len(sys.argv) > 1 else "output/cmp_FACSIMILE"
    out = sys.argv[2] if len(sys.argv) > 2 else pref + ".nc"
    csv_to_nc(pref + "_conc.csv", pref + "_irr.csv", out)
