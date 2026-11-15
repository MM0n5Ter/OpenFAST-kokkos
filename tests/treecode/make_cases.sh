#!/bin/bash
# Scratch OLAF cases for unit treecode (outside the repo), made like tests/induced-velocities/
# make_cases.sh: copies of the r-test case ad_BAR_OLAF with only TMax, nNWPanels, VelocityMethod
# (and WrVTK off) changed. The R8 tree case is unit 2's opt/build/iv-cases/bar-r8-tree (reused).
#   make_cases.sh [dest]   (default $OPT/build/u3-cases)
set -eu
REPO=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
SRC=$REPO/deps/openfast/reg_tests/r-test/modules/aerodyn
DST=${1:-/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/u3-cases}
mkdir -p $DST
ln -sfn $SRC/BAR_Baseline $DST/BAR_Baseline
mk() {  # name TMax nNWPanels VelocityMethod
  local d=$DST/$1
  rm -rf $d; cp -r $SRC/ad_BAR_OLAF $d; rm -f $d/*.outb $d/*.out $d/*.ech
  sed -i -E "s/^( *)6\.0( +TMax)/\1$2\2/" $d/ad_driver.dvr
  sed -i -E "s/^35( +nNWPanels)/$3\1/; s/^2( +VelocityMethod)/$4\1/; s/^2( +WrVTk)/0\1/" $d/OLAF.dat
  grep -E "TMax" $d/ad_driver.dvr | head -1; grep -E "nNWPanels|VelocityMethod|WrVTk" $d/OLAF.dat | head -3
}
mk bar-perf-tree 30.0 300 2   # P2/P3: tree wake (the default method), larger wake than the R8 case
