#!/bin/sh
# ============================================================
# smoke_test.sh —— 工程版 V3 冒烟 (fac_v3 单路线; fac_run/kpp 已退役)
# ============================================================
set -e
HERE=$(cd "$(dirname "$0")/.." && pwd)
cd $HERE
mkdir -p output
PASS=0; FAIL=0; TOTAL=0

echo "===== 冒烟测试 (工程版 V3) ====="


# ---------- 工程版 V3 冒烟 (mechanism_3 + runtime.dat 通路) ----------
cd "$HERE/facsimile_impl"
./fac_v3 SCEN_BAU /tmp/smoke_v3 ../runtime_v3.dat > /tmp/smoke_v3.log 2>&1
if [ $? -eq 0 ] && [ $(tail -1 /tmp/smoke_v3.log | grep -c "ok") -eq 1 ]; then
  NR=$(wc -l < /tmp/smoke_v3_conc.csv)
  O3=$(awk -F, 'END{print $4}' /tmp/smoke_v3_conc.csv)
  echo "  [PASS] FAC-V3 冒烟 (conc $NR 行, O3末=$O3)"; PASS=$((PASS+1)); TOTAL=$((TOTAL+1))
else
  echo "  [FAIL] FAC-V3 冒烟"; tail -3 /tmp/smoke_v3.log; exit 1
fi
echo "===== 冒烟: $PASS/$TOTAL PASS ====="
[ $FAIL -eq 0 ] || exit 1
