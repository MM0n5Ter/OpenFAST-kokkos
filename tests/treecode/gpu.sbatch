#!/bin/bash
# GPU jobs of unit treecode. One of:
#   sbatch -J ofk-treecode-mig  -p interactive_gpu --gres=gpu:1 -c 4 gpu.sbatch mig a100
#       correctness on an A100 MIG slice: l0_tc (claim 2: strict CUDA vs Fortran; R2 fast; claim 3;
#       P6), R7 across processes, unit 2's l0_iv with the tree ported (u3-iv-cuda-a100), and
#       (if ov3-cuda-a100 exists) the BAR tree-wake case on the GPU (R8, R7 two processes).
#   sbatch -J ofk-treecode-perf -p h100,h200 --gres=gpu:1 -c 8 gpu.sbatch perf h100
#   sbatch -J ofk-treecode-perf -p a100      --gres=gpu:1 -c 8 gpu.sbatch perf a100
#       P1/P5/P3 timing sweep N = 1e3..1e6 (perf_tc, port only, fast and strict), with and without
#       the per-phase profile (build / tree H2D / CP H2D / kernel / D2H / host add). Fast runs with the
#       device build (default) and with the host build (OFK_TC_DEVBUILD=0) for the A/B.
# Trees (environment, default = the original unit-3 trees): U3=<prefix of the test trees, u3 or u3d>,
# OV=<prefix of the full OpenFAST trees, ov3>. Unit 3's device build used U3=u3d.
#SBATCH -t 02:00:00
#SBATCH -o /weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/u3-gpu.%x.%j.log
set -u
what=$1; arch=$2
source /etc/profile.d/zz-dsai_lmod.sh
module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28 >/dev/null 2>&1
# GPU nodes are AMD EPYC: the module's OpenBLAS (Sapphire Rapids) dies with SIGILL in dgemm_
export LD_LIBRARY_PATH=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/openblas-0.3.28-portable/lib:$LD_LIBRARY_PATH
B=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build
R=/weka/scratch/jhu/zli426/zzhan381/Transformation/OpenFAST-kokkos
T=$R/tests/treecode; TI=$R/tests/induced-velocities
OUT=$B/u3-$what-$arch-$SLURM_JOB_ID; mkdir -p $OUT
PY=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/venv-rtest/bin/python; L=$R/deps/openfast/reg_tests/lib; export PYTHONDONTWRITEBYTECODE=1
export OMP_PROC_BIND=false OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
hostname; nvidia-smi -L; lscpu | grep "Model name"; date
U3=${U3:-u3}; OV=${OV:-ov3}
TC=$B/$U3-cuda-$arch; IVT=$B/$U3-iv-cuda-$arch
echo "trees: $TC $IVT $B/$OV-cuda-$arch"
case $what in
mig)
  echo "=== l0_tc on GPU ($TC)"
  ( cd $TC && ctest -E "ffma" --timeout 1800 > $OUT/tc.ctest.log 2>&1 ); grep -E "Test +#|tests passed" $OUT/tc.ctest.log
  grep -hE "^l0_tc|claim|R2:|P6|FAIL|R7" $TC/Testing/Temporary/LastTest.log | head -40
  if [ -d $IVT ]; then
    echo "=== unit 2 l0_iv on GPU with the tree ported ($IVT)"
    ( cd $IVT && ctest -E "kkoff|ffma|refref" --timeout 1800 > $OUT/iv.ctest.log 2>&1 ); grep -E "Test +#|tests passed" $OUT/iv.ctest.log
    grep -hE "^l0_iv: .*passed|^  (TreePart|Basic|Part|A7) |P6 \(L0\)|FAIL" $IVT/Testing/Temporary/LastTest.log | head -40
    echo "=== unit 2 l0_iv on GPU, device tree build forced (OFK_TC_DEVBUILD=1)"
    ( cd $IVT && OFK_TC_DEVBUILD=1 OFK_TC_DEVBUILD_CHECK=1 ctest -E "kkoff|ffma|refref" --timeout 1800 > $OUT/iv-db.ctest.log 2>&1 ); grep -E "Test +#|tests passed" $OUT/iv-db.ctest.log
    grep -hE "^l0_iv: .*passed|^  TreePart |P6 \(L0\)|FAIL" $IVT/Testing/Temporary/LastTest.log | head -20
  fi
  D=$B/$OV-cuda-$arch/modules/aerodyn/aerodyn_driver
  if [ -x $D ]; then
    C=$B/iv-cases; MT=$IVT/libofk_memtool.so
    for m in strict fast; do bash $TI/run_case.sh $D $C/bar-r8-tree $OUT/bar-r8-tree.gpu-$m OFK_MODE=$m OFK_IV_PROFILE=1; done
    bash $TI/run_case.sh $D $C/bar-r8-tree $OUT/bar-r8-tree.gpu-fast-b OFK_MODE=fast
    bash $TI/run_case.sh $D $C/bar-r8-tree $OUT/bar-r8-tree.gpu-fastdb OFK_MODE=fast OFK_TC_DEVBUILD=1
    bash $TI/run_case.sh $D $C/bar-r8-tree $OUT/bar-r8-tree.gpu-fastdb-b OFK_MODE=fast OFK_TC_DEVBUILD=1
    bash $TI/run_case.sh $D $C/bar-r8-tree $OUT/bar-r8-tree.gpu-fastdb-mem OFK_MODE=fast OFK_TC_DEVBUILD=1 OFK_IV_PROFILE=1 KOKKOS_TOOLS_LIBS=$MT OFK_MEMTOOL_OUT=$OUT/memtool.tree.txt
    tail -1 $OUT/memtool.tree.txt
    echo "--- P6 (GPU, full run, device build): allocations with bytes > 0 after evaluation 3: $(awk '$1=="alloc" && $4>0 && $5>3' $OUT/memtool.tree.txt | wc -l)  (0-byte events: $(awk '$1=="alloc" && $4==0' $OUT/memtool.tree.txt | wc -l))"
    awk '$1=="alloc" && $4>0 && $5>3' $OUT/memtool.tree.txt | head -5
    grep -h "ofk_iv\]" $OUT/bar-r8-tree.gpu-fast/stdout.log $OUT/bar-r8-tree.gpu-fastdb-mem/stdout.log | head -4
    same() { cmp -s <(tail -n +8 $1) <(tail -n +8 $2) && echo identical || echo DIFFER; }
    echo "--- R7 (GPU fast, two processes, from line 8): auto $(same $OUT/bar-r8-tree.gpu-fast/ad_driver.out $OUT/bar-r8-tree.gpu-fast-b/ad_driver.out); device build $(same $OUT/bar-r8-tree.gpu-fastdb/ad_driver.out $OUT/bar-r8-tree.gpu-fastdb-b/ad_driver.out)"
    echo "--- R8 bar-r8-tree on GPU: base1 = Fortran 1 thread (job 1860700), envelope = Fortran built with FMA (job 1860993)"
    $PY $TI/r8.py $L $B/iv-r8-1860700/bar-r8-tree.fomp1/ad_driver.out $B/iv-r8ffma-1860993/bar-r8-tree.ffma/ad_driver.out \
        gpu-strict=$OUT/bar-r8-tree.gpu-strict/ad_driver.out gpu-fast=$OUT/bar-r8-tree.gpu-fast/ad_driver.out \
        gpu-fastdb=$OUT/bar-r8-tree.gpu-fastdb/ad_driver.out 2>&1 | grep -E "max drift|envelope max"
    for r in gpu-strict gpu-fast gpu-fastdb; do echo "    wall $r: $(cat $OUT/bar-r8-tree.$r/time.txt) s"; done
  fi
  ;;
perf)
  cp $TC/perf_tc $OUT/   # the binary this job measured
  for n in 1000 3000 10000 30000 100000 300000 1000000; do
    for m in fast strict; do
      OFK_MODE=$m $OUT/perf_tc $n 3 1 1.5 port 2>&1 | grep csv | sed "s/^/$arch,/" | tee -a $OUT/sweep.csv
      OFK_MODE=$m OFK_IV_PROFILE=1 $OUT/perf_tc $n 3 1 1.5 port > $OUT/prof.$m.$n.log 2>&1
      grep "report" $OUT/prof.$m.$n.log | sed "s/^/  $m n=$n /;s/H2D=.*lattice.*| t\[s\]/| t[s]/"
    done
    # A/B: fast with the host build (uploaded) vs the device build (forced, also below OFK_TC_DEVBUILD_MIN)
    for db in 0 1; do
      OFK_MODE=fast OFK_TC_DEVBUILD=$db $OUT/perf_tc $n 3 1 1.5 port 2>&1 | grep csv | sed "s/^/$arch,devbuild$db,/" | tee -a $OUT/sweep.csv
      OFK_MODE=fast OFK_TC_DEVBUILD=$db OFK_IV_PROFILE=1 $OUT/perf_tc $n 3 1 1.5 port > $OUT/prof.fastdb$db.$n.log 2>&1
      grep "report" $OUT/prof.fastdb$db.$n.log | sed "s/^/  fast devbuild=$db n=$n /;s/H2D=.*lattice.*| t\[s\]/| t[s]/"
    done
  done
  ;;
esac
date
