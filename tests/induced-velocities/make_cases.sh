#!/bin/bash
# Scratch OLAF cases for unit induced-velocities, outside the repo: copies of the r-test case
# ad_BAR_OLAF with only TMax, nNWPanels, VelocityMethod (and WrVTK off; for the R8 cases the
# output format: text at full precision) changed.
#   make_cases.sh [dest]   (default $OPT/build/iv-cases)
set -eu
REPO=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
SRC=$REPO/deps/openfast/reg_tests/r-test/modules/aerodyn
DST=${1:-/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/iv-cases}
mkdir -p $DST
ln -sfn $SRC/BAR_Baseline $DST/BAR_Baseline     # the case reads ../BAR_Baseline/...
mk() {  # name TMax nNWPanels VelocityMethod [fullprec]
  local d=$DST/$1
  rm -rf $d; cp -r $SRC/ad_BAR_OLAF $d; rm -f $d/*.outb $d/*.out $d/*.ech
  sed -i -E "s/^( *)6\.0( +TMax)/\1$2\2/" $d/ad_driver.dvr
  sed -i -E "s/^35( +nNWPanels)/$3\1/; s/^2( +VelocityMethod)/$4\1/; s/^2( +WrVTk)/0\1/" $d/OLAF.dat
  if [ "${5:-}" = fullprec ]; then  # R8: text output at full double precision (.outb stores 16-bit scaled channels)
    sed -i -E 's/^"ES15\.8E2"( +OutFmt)/"ES24.16E3"\1/; s/^2( +OutFileFmt)/1\1/' $d/ad_driver.dvr
  fi
  grep -E "TMax|OutFmt|OutFileFmt" $d/ad_driver.dvr | head -3; grep -E "nNWPanels|VelocityMethod|WrVTk" $d/OLAF.dat | head -4
}
mk bar-r8-tree    30.0 120 2 fullprec   # R8: default wake method (tree, Fortran), LL direct (port)
mk bar-r8-direct  30.0 120 1 fullprec   # R8: direct wake (port carries the wake)
mk bar-perf-direct 30.0 300 1  # P1/P3: direct wake, larger
mk bar-perf-part  30.0 300 3   # P1/P3: direct particles
