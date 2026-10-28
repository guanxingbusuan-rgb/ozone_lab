#!/bin/sh
# ozone_lab 环境自检 —— 只读检查, 不修改任何东西
# 锚定到工程根（本脚本已并入根目录），从任意目录运行均可
cd "$(dirname "$0")" || exit 1
# 优先用项目 venv 的 Python（与 Makefile 一致），否则回退 PATH 里的 python3
if [ -x "$HOME/miniconda3/envs/venv/bin/python" ]; then
  PY=${PY:-$HOME/miniconda3/envs/venv/bin/python}
else
  PY=${PY:-python3}
fi
ok=0; bad=0
chk() {  # chk <描述> <命令>
  if eval "$2" >/dev/null 2>&1; then echo "  [OK]   $1"; ok=$((ok+1));
  else echo "  [MISS] $1"; bad=$((bad+1)); fi
}

echo "== Python =="
chk "python 可用 ($($PY --version 2>&1 | cut -d' ' -f2))" "$PY --version"
chk "numpy"      "$PY -c 'import numpy'"
chk "matplotlib" "$PY -c 'import matplotlib'"
chk "scipy (可选)"      "$PY -c 'import scipy'"
chk "netCDF4 (可选)"    "$PY -c 'import netCDF4'"

echo "== 编译器 =="
GF=${FC:-$(command -v gfortran || echo "$HOME/gcc-13.2.0/bin/gfortran")}
chk "gfortran ($($GF --version 2>/dev/null | head -1))" "$GF --version"
chk "gcc" "gcc --version"

echo "== CVODE =="
CV=${CVODE_ROOT:-$HOME/soft/atmenv/cvode}
chk "libsundials_cvode"     "test -e $CV/lib/libsundials_cvode.dylib -o -e $CV/lib/libsundials_cvode.so"
chk "libsundials_fcvode (fcmix Fortran 接口)" \
    "test -e $CV/lib/libsundials_fcvode.a -o -e $CV/lib/libsundials_fcvode.so"
chk "cvode.h 头文件"        "test -e $CV/include/cvode/cvode.h"

echo "== 项目构建 =="
chk "fac_v3 二进制存在"       "test -x facsimile_impl/fac_v3"
chk "fac_batch_omp 二进制存在" "test -x facsimile_impl/fac_batch_omp"

echo "== 二进制依赖可解析 =="
if command -v otool >/dev/null; then
  n=$(otool -L facsimile_impl/fac_v3 2>/dev/null | grep -c "not found")
  [ "$n" = "0" ] && { echo "  [OK]   fac_v3 无未解析依赖"; ok=$((ok+1)); } \
                 || { echo "  [MISS] fac_v3 有 $n 个未解析依赖"; bad=$((bad+1)); }
fi

echo
echo "== 汇总: $ok OK, $bad MISS =="
[ "$bad" = "0" ] && echo "环境就绪 ✓" || echo "请先补齐 MISS 项"
