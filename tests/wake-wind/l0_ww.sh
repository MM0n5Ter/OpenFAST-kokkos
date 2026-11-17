#!/bin/bash
# Unit 4b (wake-wind) L0: aerodyn_driver on variants of r-test ad_BAR_OLAF (and ad_QuadRotor_OLAF), each run three times
# in ONE binary: OFK_WW=fortran (the pin's code path) twice (reproducibility) and the default (H2-H4 on).
# Every output file (.outb, all VTK incl. grids, OFK_WW_DUMP = m%Inflow(1:2) wake velocities and the read V_wind slices
# at every AD_UpdateStates) must be byte-identical. One verdict line per case; details in <work>/l0_ww.log.
#   tests/wake-wind/l0_ww.sh <aerodyn_driver> <work dir> [case ...]
set -u
EXE=$1; W=$2; shift 2
R=$(cd "$(dirname "$0")/../.." && pwd); RT=$R/deps/openfast/reg_tests/r-test/modules/aerodyn
LOG=$W/l0_ww.log; mkdir -p $W/cases $W/runs; : > $LOG
for x in cases runs; do ln -sfn $RT/BAR_Baseline $W/$x/BAR_Baseline; done  # the BAR cases read ../BAR_Baseline
module load gcc/12.3.0 openblas/0.3.28 >/dev/null 2>&1 || true
export OFK_MODE=strict OMP_NUM_THREADS=1
# ---- case builders (BAR base: TMax 2 s, DT = DTfvw = 0.1 s, uniform steady wind, TwrPotent 1, TwrShadowOnWake)
grids() {  # grids <OLAF.dat> "<row>" ...  (row: name type tS tE dt x0 x1 nx y0 y1 ny z0 z1 nz)
  python3 -I - "$@" <<'PY'
import sys; f=sys.argv[1]; rows=sys.argv[2:]; L=open(f).read().split('\n')
i=[k for k,l in enumerate(L) if 'nGridOut' in l][0]; j=i+3
while not L[j].startswith('==='): j+=1
L[i]='%-8d nGridOut           - Number of grid outputs'%len(rows); L[i+3:j]=rows
open(f,'w').write('\n'.join(L))
PY
}
mk() {  # mk <name> <base dir>
  d=$W/cases/$1; rm -rf $d; mkdir -p $W/cases; cp -r $2 $d; echo $d; }
build_case() {
  local c=$1 d
  case $c in
    quad) d=$(mk $c $RT/ad_QuadRotor_OLAF); sed -i -E '0,/TMax/s/^( *)[0-9.]+( +TMax)/\12.0\2/' $d/*.dvr; return;;
  esac
  d=$(mk $c $RT/ad_BAR_OLAF); sed -i -E 's/^( *)6\.0( +TMax)/\12.0\2/' $d/ad_driver.dvr
  case $c in
    base) ;;
    grid3d) python3 -I $R/tests/wake-wind/mkbts.py $d/grid.bts
            sed -i -E 's/^( *)0( +CompInflow)/\11\2/' $d/ad_driver.dvr
            sed -i -E 's/^1( +WindType)/3\1/; s/^"NA"( +Filename_BTS)/"grid.bts"\1/' $d/OpenFAST_BAR_00_InflowFile.dat;;
    g0)     grids $d/OLAF.dat;;
    g2)     grids $d/OLAF.dat '"hA" 2 1.2 1.4 default -80 170. 30 -120. 120. 20 210. 210. 1' \
                              '"vB" 2 0.5 1.0 0.2 -20 60. 9 -60. 60. 7 20. 200. 5';;
    fwlive) # far wake alive before TMax (with 35 NW panels it starts only at 3.5 s): 8 NW panels (needs FWShedVorticity)
            sed -i -E 's/^35( +nNWPanels)/8 \1/; s/^default( +FWShedVorticity)/True\1/' $d/OLAF.dat; grep -q '^8  *nNWPanels' $d/OLAF.dat || exit 3;;
    fw0)    sed -i -E 's/^10( +nFWPanels)/0 \1/; s/^10( +nFWPanelsFree)/0 \1/' $d/OLAF.dat;;
    dtfvw)  sed -i -E 's/^0\.10( +DTfvw)/0.20\1/' $d/OLAF.dat;;
    srcpnl) python3 -I - $RT/ad_Sphere_OLAF/Sphere_17_15.vtk $d/Sphere.vtk <<'PY'
import sys; L=open(sys.argv[1]).read().split('\n'); i=[k for k,l in enumerate(L) if l.startswith('POINTS')][0]
n=int(L[i].split()[1]); v=[]; j=i+1
while len(v)<3*n: v+=[float(x) for x in L[j].split()]; j+=1
assert len(v)==3*n                                     # whole triplets on the consumed lines
for k in range(n): v[3*k]+=-10.0; v[3*k+2]+=95.0      # sphere centre (-10, 0, 100): upstream of the tower, in its potential flow
L[i+1:j]=['%.9g %.9g %.9g'%tuple(v[3*k:3*k+3]) for k in range(n)]
open(sys.argv[2],'w').write('\n'.join(L))
PY
            sed -i -E '/^-+ ADVANCED OPTIONS/a "Sphere.vtk"    SrcPnlFile     - source panels\n1    nSrcPnlUpdate  - panel update interval' $d/OLAF.dat;;
    gedge)  # W4: windows on, just inside and just outside step boundaries (DT 0.1; bDoGrid tolerance DT/2)
            r='-20 60. 6 -40. 40. 5 40. 160. 2'
            grids $d/OLAF.dat "\"e1\" 2 1.2 1.2 default $r" "\"e2\" 2 1.15 1.15 default $r" "\"e3\" 2 1.1500001 1.2499999 default $r" \
                              "\"e4\" 2 0.6499999 0.6500001 default $r" "\"e5\" 2 0.0 0.0 default $r" "\"e6\" 2 1.9 2.5 default $r";;
    gend)   # FVW_FinalWrite reads u(1)'s grid slice at (n+1) DT = TMax: grids starting just before/after TMax + DT/2
            r='-20 60. 6 -40. 40. 5 40. 160. 2'
            grids $d/OLAF.dat "\"f1\" 2 2.04 3.0 default $r" "\"f2\" 2 2.0499999 3.0 default $r" "\"f3\" 2 2.06 3.0 default $r" "\"f4\" 2 1.95 1.95 default $r";;
    *) echo "unknown case $c"; exit 2;;
  esac
}
run() {  # run <case dir> <tag> [env...]
  local d=$1 t=$2; shift 2; local o=$W/runs/$(basename $d).$t
  rm -rf $o; mkdir -p $W/runs; cp -r $d $o
  ( cd $o && env "$@" OFK_WW_DUMP=$o/ww.dump "$EXE" $(ls *.dvr | head -1) > stdout.log 2>&1 ); echo $? > $o/rc
}
same() {  # same <dirA> <dirB>: every output file byte-identical (inputs are copies)
  local a=$1 b=$2 n=0 bad=0
  for f in $(cd $a && find . -type f ! -name stdout.log | sort); do
    n=$((n+1)); cmp -s $a/$f $b/$f || { bad=$((bad+1)); echo "    DIFF $f" >> $LOG; }
  done
  [ "$(cd $a && find . -type f | sort)" = "$(cd $b && find . -type f | sort)" ] || { bad=$((bad+1)); echo "    file lists differ" >> $LOG; }
  echo "$n $bad"
}
declare -A ID=([base]=l0_ww_wake [grid3d]=l0_ww_wake [quad]=l0_ww_wake [g0]=l0_ww_slices [g2]=l0_ww_slices [fw0]=l0_ww_slices [fwlive]=l0_ww_slices \
               [srcpnl]=l0_ww_slices [dtfvw]=l0_ww_slices [gedge]=l0_ww_grid [gend]=l0_ww_grid)
CASES=${*:-base grid3d quad g0 g2 fw0 fwlive srcpnl dtfvw gedge gend}
rc=0
for c in $CASES; do
  echo "=== $c" >> $LOG
  d=$(build_case $c | tail -1); [ -d "$d" ] || d=$W/cases/$c
  run $d ref1 OFK_WW=fortran; run $d ref2 OFK_WW=fortran; run $d ofk
  A=$W/runs/$c.ref1; B=$W/runs/$c.ref2; C=$W/runs/$c.ofk
  r="rc $(cat $A/rc)/$(cat $B/rc)/$(cat $C/rc)"
  read n1 b1 < <(same $A $B); read n2 b2 < <(same $A $C)
  nv=$(cd $A && find . -name '*.vtk' -newer ad_driver.dvr | wc -l); ng=$(cd $A && find . -name '*Grid*.vtk' | wc -l)
  db=$(stat -c %s $A/ww.dump 2>/dev/null || echo 0)
  if [ "$(cat $A/rc)$(cat $B/rc)$(cat $C/rc)" = "000" ] && [ $b1 = 0 ] && [ $b2 = 0 ] && [ $db -gt 0 ]; then v=PASS; else v=FAIL; rc=1; fi
  echo "L0 ${ID[$c]} $c: $v ($r; ref/ref $b1 diff, pin-path/4b $b2 diff of $n2 files; $nv vtk of which $ng grid; dump $db B)"
done
exit $rc
