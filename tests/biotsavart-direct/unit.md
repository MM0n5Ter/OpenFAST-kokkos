# biotsavart-direct

The OLAF direct Biot–Savart kernels `ui_seg` (vortex segments) and `ui_part_nograd` (vortex particles), ported to C++/Kokkos behind their frozen Fortran interfaces. Unit 1 of 3 (architecture.md §7): it establishes the shim/oracle/L0 pattern, the one-body rule, the strict and fast modes, the dispatcher, the FMA switch and the fork build integration.

- **Status (v1, 2026-11-14)**: in the fork, built in all three configurations. Chain claims 1 and 2 hold; claim 3 holds for `ui_part_nograd` and **fails at the declared R1 bound for `ui_seg` because of FMA** (needs a human decision, see "Claim 3 and FMA"). P2 passes in fast mode. R6, R7, A1/A2/A3/A5 pass. P6 fails (known; unit 2).
- **Owner**: agent session `199ad883` (subagent, unit 1). Date: 2026-11-14.
- **Toolchain**: Kokkos 5.3.0 (C++20), gcc/gfortran 12.3.0, cmake 3.31.12, CUDA 12.3, OpenBLAS 0.3.28.

## Files

| Where | What |
|---|---|
| `fork/modules/aerodyn/src/FVW_BiotSavart.f90` | The original file, edited in place: only the bodies of `ui_seg` and `ui_part_nograd` gain `#ifdef OPENFAST_KOKKOS` → `call ofk_…` / `#else` → the original body (+78 lines). Interface blocks and `use iso_c_binding` are local to each routine, so the module's public names are unchanged. |
| `fork/modules/kokkos/src/biotsavart_direct_body.hpp` | **The physics, once**: `seg_term<Reg,T>` and `part_term<Reg,T>`, one target–source contribution each, templated on the regularization branch. Every shape and backend calls these. |
| `fork/modules/kokkos/src/biotsavart_direct_impl.hpp` | Loop shapes (square with `kLanes` partial sums, skinny team reduction), SoA packing, per-mode dispatcher. |
| `fork/modules/kokkos/src/biotsavart_direct_{strict,fast}.cpp` | The two translation units: 7 lines each, `#define OFK_BSD_MODE_NS strict|fast`, `#define OFK_BSD_STRICT 1|0`, `#include "biotsavart_direct_impl.hpp"`. Only the macros and the compile flags differ. |
| `fork/modules/kokkos/src/biotsavart_direct.cpp` | `extern "C"` entry points, argument checks, host↔device views, `OFK_MODE` switch, profiling. No arithmetic on results. |
| `fork/modules/kokkos/src/biotsavart_direct{,_common}.hpp` | C API; shared data types (constants, view types, argument structs). |
| `fork/modules/kokkos/src/ofk_runtime.{hpp,cpp}` | Kokkos lifetime; `ofk::mode()` (reads `OFK_MODE` once, prints `[ofk] OFK_MODE=…` once); `ofk_set_mode` test hook. |
| `fork/modules/kokkos/src/ofk_vector_libm.hpp` | Optional SIMD-clone `exp` (glibc libmvec) for fast host TUs. Off by default (measured slower). |
| `fork/modules/kokkos/CMakeLists.txt` | `ofk_kernels`; per-TU FMA flags; `OFK_FAST_FMA`, `OFK_FAST_HOST_FLAGS`, `OFK_FAST_VECTOR_LIBM`; `kokkos_compilation(GLOBAL)` for CUDA; install/export. |
| `fork/CMakeLists.txt`, `fork/modules/aerodyn/CMakeLists.txt` | `option(OPENFAST_KOKKOS OFF)`, `find_package(Kokkos 5.3...<6)`, `add_subdirectory(modules/kokkos)`; `aerodynlib` gets `OPENFAST_KOKKOS` and links `ofk_kernels`. |
| `tests/biotsavart-direct/` | `CMakeLists.txt` (oracle from the pin, port from the fork), `l0_util.f90`, `l0_ui_seg.f90`, `l0_ui_part_nograd.f90`, `api_fvw_runtests.f90`, `perf_direct.f90`, `r7_threads.sh`, `cmp_outb.py`, `cpu.sbatch`, `gpu.sbatch`, `gpu2.sbatch`, `r6.sbatch`. |

## Contract (frozen, pin `12fa3fc`)

```fortran
subroutine ui_seg(iCPStart, iCPEnd, CPs, iSegStart, iSegEnd, SegPoints, SegConnct, SegGamma, RegFunction, RegParam, Uind_out)
   real(ReKi),     dimension(:,:), intent(in)    :: CPs        ! 3 x nCPs++        (n++: only iCPStart:iCPEnd used)
   integer(IntKi),                 intent(in)    :: iCPStart, iCPEnd, iSegStart, iSegEnd, RegFunction
   real(ReKi),     dimension(:,:), intent(in)    :: SegPoints  ! 3(+) x nSegPTot   (only columns named in SegConnct)
   integer(IntKi), dimension(:,:), intent(in)    :: SegConnct  ! 4 x nSeg++: rows 1-2 = 1-BASED point indices, rows 3-4 = depth/span (unused here)
   real(ReKi),     dimension(:),   intent(in)    :: SegGamma, RegParam   ! nSeg++  (RegParam not read for RegFunction=0)
   real(ReKi),     dimension(:,:), intent(inout) :: Uind_out   ! 3(+) x nCPs++: ACCUMULATED into, columns iCPStart:iCPEnd, rows 1:3 only
! FVW_BiotSavart.f90:135-357.  RegFunction 0..4 (None, Rankine, LambOseen, Vatistas, Offset); other -> print + STOP.

subroutine ui_part_nograd(nCPs, CPs, nPart, Part, Alpha, RegFunction, RegParam, UIout)
   integer(IntKi),             intent(in)    :: nCPs, nPart, RegFunction
   real(ReKi), dimension(:,:), intent(in)    :: CPs, Part, Alpha  ! 3(+) x n++ (only 1:nCPs / 1:nPart used)
   real(ReKi), dimension(:),   intent(in)    :: RegParam          ! nPart++ (not read for RegFunction=0)
   real(ReKi), dimension(:,:), intent(inout) :: UIout             ! 3(+) x nCPs++: ACCUMULATED into, columns 1:nCPs, rows 1:3
! FVW_BiotSavart.f90:360-455.  RegFunction 0..2 (None, Exp, Compact); other -> print + STOP.
```

**C boundary**: `ofk_biotsavart_direct_{ui_seg,ui_part_nograd}_{r8,r4}`. Each 2-D array is passed with its leading dimension `size(A,1)`; ranges and `SegConnct` values stay 1-based (offset at use); `ui_part_nograd` receives `MINNORM` as compiled (audit #1); a Fortran generic interface picks `_r8`/`_r4` from the kind of `ReKi`.

**Callers**: `FVW_Subs.f90:1303,1311` (`InducedVelocitiesAll_Calc`), `:1564,1568` (`LiftingLineInducedVelocities`), `FVW_Tests.f90`. No use outside AeroDyn/OLAF.

## How the architecture is realized

| Concern | Realization |
|---|---|
| One body | `seg_term` / `part_term` (body.hpp). The guards are the Fortran's comparisons, evaluated branch-free (`ok = (den > P) & (n2 > P) & (r0 > P)`; Vatistas `ok = !(den <= P) & …`), then the caller adds `ok ? o : 0`. That is bitwise the Fortran's skip (the accumulator starts at +0 and can never be −0), keeps the NaN difference between the two guard styles, and lets the lane loop if-convert. Fortran `MAX`/`MIN` are compare-selects (`fmax_f`), not `Kokkos::fmax`: same value for non-NaN, and GCC cannot vectorize IEEE `fmax`. |
| Loop shapes | square: `RangePolicy` over targets, `kLanes` partial sums over sources; skinny (fast only): `TeamPolicy`, one team per (target, chunk of 1024 sources), fixed team size 128 on device, then an ordered sum over chunks. Dispatcher: skinny if device, nCP ≤ 1024 and nSeg ≥ 4096 (`OFK_BSD_SHAPE`, `OFK_BSD_SKINNY_MAX_CP`, `OFK_BSD_SKINNY_MIN_SEG`, `OFK_BSD_CHUNK`, `OFK_BSD_TEAM`). Strict always runs square. |
| Data layout | Both modes pack the sources into SoA columns once per call (a copy: values unchanged), so the strict mode is still bitwise and the fast lane loop reads unit-stride data. |
| strict | `kLanes = 1`, `-ffp-contract=off`, nvcc `--fmad=false`. |
| fast (default) | `kLanes = 8` on host backends (`omp simd` over the lanes; a compile-time constant, so the result is independent of thread count and SIMD width), `kLanes = 1` on device; `-ffp-contract=fast -march=x86-64-v3` (FMA needs an ISA that has it), nvcc `--fmad=true`. |
| ODR | Each mode TU wraps body + shapes in `namespace ofk_bsd::strict` / `::fast`. Without that, the linker keeps one COMDAT copy of `seg_term<0,double>` and one mode silently runs the other's arithmetic. |
| Switches | `OFK_MODE=strict|fast` (runtime, printed once); CMake `OFK_FAST_FMA` (ON), `OFK_FAST_HOST_FLAGS` (`-march=x86-64-v3`), `OFK_FAST_VECTOR_LIBM` (OFF). |

Vectorization evidence (GCC 12, `-fopt-info-vec`, fast TU): `impl.hpp:124` (seg lane loop) and `:284` (part lane loop) "loop vectorized using 32 byte vectors" for every branch except the exp-based ones (LambOseen, Exp), whose `exp` call stays scalar. Three changes were needed to get there: `Kokkos::fmax` → compare-select; branch-free guard evaluation (otherwise the loads of Γ/RegParam sink into the branch and need masked loads); column pointers copied into the lambda body.

## Equivalence classes

| Output | Mode / backend | Class | Measured (worst over seeds 1 and 7, every branch) |
|---|---|---|---|
| both | strict, Serial/OpenMP | **B** vs Fortran OPENMP=OFF | 0, no bit difference (incl. −0/NaN pattern) |
| both | strict, CUDA | claim 2: few ulp of Σ\|c\| | 0 for None/Rankine/Vatistas/Offset/part None/Compact; 3.1e-16 (LambOseen), 8.6e-16 (Exp): only the device `exp` differs |
| `UIout` | fast vs strict, OpenMP and CUDA | **R1**, 1e-12 of Σ\|c\| | 1.05e-14: PASS |
| `Uind_out` | fast vs strict, OpenMP and CUDA | **R1**, 1e-12 of Σ\|c\| | **1.1e-11 (None, seed 7), 2.65e-12 (LambOseen), 2.3e-12 (Rankine): FAIL** |
| `Uind_out` | fast with `OFK_FAST_FMA=OFF` vs strict, OpenMP | R1, 1e-12 | 4.6e-15: PASS |

Error metric (architecture.md §4): per target, max over components of |A−B| divided by max over components of Σ_j |c_j| (the oracle called once per source), worst over targets and cases. The old max-normalized metric is printed alongside.

### Claim 3 and FMA (decision needed)

- The `ui_seg` failure is **FMA alone**: with the same blocked order but no contraction, the worst error is 4.6e-15 (job 1859825).
- The pin's own Fortran, compiled with FMA (`-march=x86-64-v3`; gfortran contracts by default), moves by **exactly the same amounts** against the default build: 2.992e-12 / 2.337e-12 (seed 1), 1.105e-11 / 2.65e-12 (seed 7), to 3–4 digits (`l0_ui_seg_ffma`). So it is a property of the segment formula, not of the port: the cross product `ya*zb - za*yb` cancels for targets near the segment's line, and FMA changes that rounding. Upstream builds on aarch64 (FMA in the baseline ISA) already differ from x86-64 builds by this much.
- With FMA, the max-normalized difference to the Fortran grows with N on random clouds: 7.6e-10 at 10⁴×10⁴, 3.0e-8 at 10⁵×10⁵ (None).
- What FMA buys (P2 table, fast vs fast-no-FMA): 37% on seg None (3.73× vs 2.73×), 1–9% on the other seg branches, nothing measurable on the particle kernels (within run-to-run noise).
- Options: (a) `OFK_FAST_FMA=OFF` by default (claim 3 PASS, P2 still PASS except part Exp at 0.966 in one run); (b) keep FMA and declare `ui_seg` fast as "R1-FMA: ≤ the Fortran's own FMA-vs-no-FMA envelope × 1.1 on the same inputs" (the R8-style envelope idea), which needs a sign-off here; (c) per-expression contraction control in the cross product (not expressible in GCC without killing vectorization). **The build default is still FMA ON**, as architecture.md specifies; ctest `l0_ui_seg_strictfast` and `l0_ui_seg_fast` fail until this is decided.

## Audit list → L0 case

| # | Trap | Pin file:line | What the port does | L0 case |
|---|---|---|---|---|
| 1 | `MINNORM = 1e-4` is a default-kind literal; the pin's double build has `-fdefault-real-8`, so it is exactly 1e-4 | `FVW_BiotSavart.f90:20`, `:377` | Passed from the shim; `MINNORM2` computed in C++. Hard-coding `1e-4f` fails L0 (mutation M2b, v0). | part `straddle thresholds` |
| 2 | `PRECISION_UI = epsilon/100` per precision | `:12` | `numeric_limits<T>::epsilon()/T(100)` | seg `straddle thresholds`; SP build (v0) |
| 3 | `MINDENOM = 0` is added | `:18` | Kept literally (equivalent mutant) | — |
| 4 | `SegConnct` 1-based; rows 3/4 depth/span; mirror adds nSegP | `FVW_Subs.f90:980` | Offset at use in the SoA pack; out-of-range indices abort | seg `lattice`, `lattice + ground mirror`, `sub-range` |
| 5 | `Uind` accumulated; sentinel tails | `:148` | Copy-in, one add per target, copy-out of target columns only | every case |
| 6 | Vatistas `<= … cycle` vs nested `>`: they differ for NaN | `:299-306` vs `:184-189` | Same comparisons, `!(x <= P)` vs `(x > P)`, evaluated branch-free | seg `NaN point in one segment` (Vatistas: NaN out; others finite) |
| 7 | Exp skips `exp` for r > 2rc; Compact `merge` evaluates both arms; Exp `cycle` on r2 < MINNORM2 | `:409`, `:415`, `:437-444` | Exp keeps the branch around `exp` (unless `OFK_FAST_VECTOR_LIBM`); `merge` → select; `MAX/MIN` → compare-select; `cycle` → `ok = !(r2 < MINNORM2)` | part `straddle thresholds`, `tiny, huge and zero core radii`, `NaN particle` |
| 8 | "3 × n++": only the used columns may be touched | `:138-148`, `:363-368` | Views cover exactly the used ranges | every case (NaN-poisoned tails) |
| 9 | Every `select case` value | `:168`, `:379` | Template per value; unknown → same print + STOP in the shim | every fixed case × every branch |
| 10 | Non-contiguous actuals | — | gfortran copy-in/out | `non-contiguous …` cases |
| 11 | Leading dimension > 3 | `:176`, `:198` | `ld` per array | `padded …` cases |
| 12 | `-fdefault-real-8` changes test literals | — | explicit `1.0e-4_4` in the test | part straddle |
| 13 | Empty ranges | `:172`, `:382` | nCP ≤ 0 returns before Kokkos; nSeg = 0 adds +0 | `nCP = 0`, `nSeg = 0`, `nPart = 0` |
| 14 | (new) FMA in the segment cross product is ill-conditioned near the segment line | `:186`, `:221`, … | strict: no FMA; fast: FMA (see above) | `l0_ui_seg_strictfast`, `l0_ui_seg_ffma` |

## Gate results (v1)

Build trees in `/weka/scratch/jhu/zli426/zzhan381/Transformation/opt/build/`: `fork-off`, `fork-omp`, `fork-cuda` (full fork), `openfast-pin` (pin, oracle libraries), `biotsavart-direct-v1-omp`, `biotsavart-direct-v1-cuda` (this test tree), `biotsavart-direct-v1-omp-nofma` (`OFK_FAST_FMA=OFF`).

```bash
source /etc/profile.d/zz-dsai_lmod.sh; module load gcc/12.3.0 cmake/3.31.12 cuda/12.3.0 openblas/0.3.28
OPT=/weka/scratch/jhu/zli426/zzhan381/Transformation/opt; export NVCC_WRAPPER_DEFAULT_COMPILER=g++
F="-DCMAKE_BUILD_TYPE=Release -DDOUBLE_PRECISION=ON -DOPENMP=OFF -DBUILD_TESTING=ON -DBLA_VENDOR=OpenBLAS"
cmake -S fork -B $OPT/build/fork-off  $F
cmake -S fork -B $OPT/build/fork-omp  $F -DOPENFAST_KOKKOS=ON -DKokkos_ROOT=$OPT/kokkos-5.3.0-omp -DCMAKE_CXX_COMPILER=g++
cmake -S fork -B $OPT/build/fork-cuda $F -DOPENFAST_KOKKOS=ON -DKokkos_ROOT=$OPT/kokkos-5.3.0-cuda-h100 -DCMAKE_CXX_COMPILER=g++
cmake -S tests/biotsavart-direct -B $OPT/build/biotsavart-direct-v1-omp -DKokkos_ROOT=$OPT/kokkos-5.3.0-omp \
      -DCMAKE_CXX_COMPILER=g++ -DOPENFAST_PIN_BUILD=$OPT/build/openfast-pin -DBLA_VENDOR=OpenBLAS
cmake -S tests/biotsavart-direct -B $OPT/build/biotsavart-direct-v1-cuda -DKokkos_ROOT=$OPT/kokkos-5.3.0-cuda-h100 \
      -DCMAKE_CXX_COMPILER=$OPT/kokkos-5.3.0-cuda-h100/bin/nvcc_wrapper -DOPENFAST_PIN_BUILD=$OPT/build/openfast-pin \
      -DBLA_VENDOR=OpenBLAS -DOFK_L0_RTOL=1e-14
cmake --build <tree> -j 4; cd <tree> && ctest          # or: sbatch cpu.sbatch / gpu.sbatch from the tree
```

| ID | Result | Evidence |
|---|---|---|
| Chain 1 (Fortran vs strict, CPU) | **PASS** | Job 1859729 (Xeon 8480+, 8 threads): `OFK_MODE=strict OFK_L0_RTOL=0`, seeds 1 and 7, ui_seg 290/290 and ui_part 230/230, every branch worst 0, 0 cases with any bit difference. Log `biotsavart-direct-v1-omp/cpu.1859729.log`. |
| Chain 2 (strict CUDA vs strict CPU) | **PASS** | Job 1859659 (H100 NVL): Fortran vs strict on CUDA (= strict CPU by claim 1), rtol 1e-14 of Σ\|c\|: 0 for every branch without `exp`; LambOseen 3.15e-16, Exp 8.6e-16. Log `biotsavart-direct-v1-cuda/gpu.1859659.log`. |
| Chain 3 (strict vs fast), part | **PASS (R1)** | worst 1.05e-14 (OpenMP and CUDA, seeds 1, 7). |
| Chain 3 (strict vs fast), seg | **FAIL at R1 1e-12** | worst 1.1e-11 (OpenMP and CUDA identical to 4 digits). Cause FMA; see "Claim 3 and FMA". With `OFK_FAST_FMA=OFF`: 4.6e-15 PASS (job 1859825). |
| Chain 4 | not in scope | unit 2 |
| R1 coverage of shapes | PASS | Claim 3 on CUDA with `OFK_BSD_SHAPE=auto|skinny|square`: identical results per shape (the skinny team shape also differs from strict only by FMA on these inputs). |
| R2 coverage | PASS (by construction) | Every `RegFunction` value and every guard has an L0 case. |
| R3 constants | PASS | as v0 (MINNORM passed from the shim). |
| R4 sentinels / accumulate | PASS | every L0 case. |
| R6 | **PASS** | Job 1859860: 16/17 pass on fork-off, fork-omp fast and fork-omp strict alike (the 2 `*_OLAF_py` pass on OFF after a rerun; `py_ad_B1n2_OLAF` fails on all three, an upstream path assumption). Strict channel data = OFF in all 15 outputs. See "R6". |
| R7, CPU | **PASS** | In-process port-vs-port (rtol 0) both modes; across processes and `OMP_NUM_THREADS=1/4/8` (`r7_threads.sh`, md5 of the dumped outputs, 1.19 MB seg / 0.99 MB part): identical per mode, strict ≠ fast. Login node rerun after the dump fix (below). |
| R7, GPU | PASS in-process; cross-process PENDING (job 1860007) | Job 1859659: port vs port in one process, rtol 0, both modes. Its cross-process lines were vacuous (empty dumps, bug fixed afterwards). |
| P2 (fast vs Fortran -fopenmp, 8 threads) | **PASS** | 1.07–3.95× at 10⁴ and 10⁵, every branch (table below). |
| P2 strict (reported, not gated) | 0.88–1.30× | seg 1.05–1.30× (the SoA pack helps), part None 0.88–0.97×, Exp 1.0–1.1×, Compact 0.94–1.0× (the Fortran-OpenMP SIMD loop vs our serial order). |
| A1 | PASS | Fork-omp `FVW_BiotSavart.f90.o` vs pin: the 11 defined `_MOD_` symbols identical; new undefined symbols only `ofk_biotsavart_direct_*_r8`. `.mod`: only the `ARRAY_OUTER_DEPENDENCY` attribute differs (as v0). |
| A2 | PASS | All callers compile unchanged: the fork's other AeroDyn files are untouched (`git -C fork diff --stat`: 3 files, 100 insertions) and the full fork builds in three configurations; `api_fvw_runtests` (pin callers vs fork module, FVW_RunTests) passes on OpenMP and CUDA. |
| A3 | PASS | New C symbols only `ofk_*`. |
| A5 | PASS | `OPENFAST_KOKKOS=OFF`: fork-off vs pin, aerodynlib `.text` identical for 18/19 objects; `FVW_BiotSavart.f90.o` differs only in the 3 line-number immediates of the `print`/`STOP` statements (+78 lines); `fvw_biotsavart.mod` md5 identical. `l0_*_kkoff` bitwise. |
| A6 | PASS (no per-instance state) | Stateless per call; the process-wide statics are the Kokkos init flag, the mode and the env config. |
| P6 | FAIL (known) | Views and the SoA buffer are allocated per call. Unit 2 adds the handle. |

### Builds of the full fork

| Config | Result | Notes |
|---|---|---|
| fork-off | builds; `aerodyn_driver`, `openfast`, `libaerodyn_inflow_c_binding.so` | identical to the pin (A5) |
| fork-omp | builds and links all three | CMake picks the C++ linker for everything that links `aerodynlib` (static C++ dependency), adds `-lgfortran` itself; `.so` depends on `libkokkoscore.so.5.3`, `libgomp`, `libstdc++`. No manual `LINKER_LANGUAGE`. |
| fork-cuda | see "fork-cuda" below | `CMAKE_CXX_COMPILER=nvcc_wrapper` for the whole fork **fails**: nvcc cannot compile OpenFAST's own C++ (`modules/awae/src/vtk.cpp` → `nwtc-library/include/fast_float.h(3581): error: incomplete type is not allowed`). Fixed by keeping g++ and `kokkos_compilation(GLOBAL)` in `modules/kokkos/CMakeLists.txt`. |

Two more integration points that will recur for every module: OpenFAST exports its targets (`install(EXPORT OpenFASTLibraries)`), so `ofk_kernels` must be installed into that export set and its include directory must be `$<BUILD_INTERFACE:…>`.

## P2 timings (job 1859729, Xeon 8480+, 8 threads, `OMP_PROC_BIND=spread`)

Fortran = the pin's file built with `-fopenmp` (the OPENMP=ON build, `!$OMP SIMD` reductions); min over reps (5 at 10⁴, 1 at 10⁵). Speedup = Fortran time / port time.

| Routine | T × S | Reg | Fortran [s] | strict [s] | fast [s] | strict × | **fast ×** | fast no-FMA × (1859825) |
|---|---|---|---|---|---|---|---|---|
| ui_seg | 10⁴×10⁴ | None | 0.1320 | 0.1027 | 0.0354 | 1.29 | **3.73** | 2.73 |
| ui_seg | 10⁴×10⁴ | Rankine | 0.1482 | 0.1249 | 0.0510 | 1.19 | **2.90** | 2.87 |
| ui_seg | 10⁴×10⁴ | LambOseen | 0.1615 | 0.1332 | 0.1225 | 1.21 | **1.32** | 1.21 |
| ui_seg | 10⁴×10⁴ | Vatistas | 0.1710 | 0.1625 | 0.0719 | 1.05 | **2.38** | 2.35 |
| ui_seg | 10⁴×10⁴ | Offset | 0.1401 | 0.1079 | 0.0355 | 1.30 | **3.95** | 3.74 |
| ui_seg | 150×3·10⁵ (LL) | Vatistas | 0.0784 | 0.0766 | 0.0358 | 1.02 | **2.19** | 2.15 |
| ui_seg | 10⁵×10⁵ | None | 13.40 | 10.50 | 3.64 | 1.28 | **3.68** | |
| ui_seg | 10⁵×10⁵ | Rankine | 15.12 | 12.70 | 5.26 | 1.19 | **2.87** | |
| ui_seg | 10⁵×10⁵ | LambOseen | 16.41 | 13.55 | 12.55 | 1.21 | **1.31** | |
| ui_seg | 10⁵×10⁵ | Vatistas | 17.38 | 16.45 | 7.36 | 1.06 | **2.36** | |
| ui_seg | 10⁵×10⁵ | Offset | 14.09 | 10.87 | 3.62 | 1.30 | **3.89** | |
| ui_part_nograd | 10⁴×10⁴ | None | 0.0341 | 0.0352 | 0.0213 | 0.97 | **1.60** | 1.75 |
| ui_part_nograd | 10⁴×10⁴ | Exp | 0.0592 | 0.0539 | 0.0554 | 1.10 | **1.07** | 0.97 |
| ui_part_nograd | 10⁴×10⁴ | Compact | 0.0615 | 0.0613 | 0.0445 | 1.00 | **1.38** | 1.28 |
| ui_part_nograd | 10⁵×10⁵ | None | 3.062 | 3.471 | 2.239 | 0.88 | **1.37** | 1.32 |
| ui_part_nograd | 10⁵×10⁵ | Exp | 5.449 | 5.478 | 4.923 | 0.99 | **1.11** | |
| ui_part_nograd | 10⁵×10⁵ | Compact | 5.754 | 6.099 | 4.012 | 0.94 | **1.43** | |

- The exp-based branches (LambOseen, Exp) do not vectorize: their gain is the SoA layout only. Vectorizing them with glibc's libmvec (`OFK_FAST_VECTOR_LIBM=ON`) means evaluating `exp` for every pair, because the Fortran's far-field skip becomes a select; it was **slower** (1 thread, 4000×4000: Exp 0.99× vs 1.12×, LambOseen 0.74× vs 1.31×), even with the argument clamped to the used range (without the clamp, libmvec's underflow path made Exp 0.24×).
- One run per point; on this node run-to-run noise was a few percent.

GPU (job 1859659, H100 NVL vs the Fortran -fopenmp on the job's 8 EPYC cores; per-call allocation and copies included):

| Routine | T × S | Reg | strict × | fast × | Notes |
|---|---|---|---|---|---|
| ui_seg | 10⁴×10⁴ | None / LambOseen / Vatistas | 26.9 / 22.8 / 26.1 | 29.4 / 23.9 / 26.8 | |
| ui_part_nograd | 10⁴×10⁴ | None / Exp / Compact | 10.1 / 15.5 / 10.6 | 10.2 / 16.3 / 11.2 | |
| ui_seg (LL) | 150×3·10⁵ | Vatistas | **0.40** | **31.9** | strict is square only (150 threads); fast picks skinny: kernel 1.0 ms, alloc+H2D 1.2 ms |
| ui_seg | 10⁵×10⁵ | LambOseen | 114 | 121 | |
| ui_part_nograd | 10⁵×10⁵ | Exp | 98 | 101 | kernel 56 ms, alloc+H2D 0.9 ms |

## R6

Job 1859860 (`med`, 8 cores, `OMP_NUM_THREADS=1`, `ctest -R 'aerodyn_utest|OLAF' -j 8`), log `$OPT/build/r6.1859860.log`, outputs `$OPT/build/r6-1859860/{off,omp-fast,omp-strict}`. Script `r6.sbatch`. Python: a venv with numpy and bokeh (`$OPT/venv-rtest`; the r-test scripts need both, otherwise every case "fails" after the simulation completes), plus `pyOpenFAST` installed from a copy of `fork/glue-codes/python`.

| Test (17) | fork-off | fork-omp fast | fork-omp strict |
|---|---|---|---|
| `aerodyn_utest`, `aerodyn_utest_trapping_flags` | Passed | Passed | Passed |
| 10 `ad_*_OLAF` (aerodyn_driver) | Passed | Passed | Passed |
| `HelicalWake_OLAF`, `EllipticalWing_OLAF`, `5MW_Land_QuadRotor_DLL_WSt_OLAF` (openfast) | Passed | Passed | Passed |
| `HelicalWake_OLAF_py`, `EllipticalWing_OLAF_py` (libopenfastlib via Python: Kokkos inside a dlopened library) | Passed (rerun after pyOpenFAST was installed; the job's OFF run predates it) | Passed | Passed |
| `py_ad_B1n2_OLAF` | Failed | Failed | Failed |

**R6: PASS**: no case passes on the baseline and fails on the port. `py_ad_B1n2_OLAF` fails everywhere: its driver looks for `driver_utilities` at `<build>/../reg_tests/r-test`, i.e. it assumes the build tree is `<source>/build`.

Channel data of the 15 `.outb` files compared with `fast_io.load_output` (the files themselves always differ: their header has a timestamp):
- **strict vs OFF: identical in all 15** (chain claim 1 at simulation level). Caveat: `.outb` stores channels quantized (16-bit scaled), so this is identity at output resolution, not a proof of bitwise state.
- fast vs OFF: identical in the 12 `aerodyn_driver` outputs; `EllipticalWing_OLAF` 2.4e-18, `HelicalWake_OLAF` 8.5e-15, `5MW_Land_QuadRotor_DLL_WSt_OLAF` 1.6e-5 (max relative to the channel's max; 30 s of coupled simulation, inside the r-test tolerance). That last one is R8 territory (unit 2).

## GPU job 2 (pending)

Job **1860007** `ofk-biotsavart-direct-gpu2` (`gpu2.sbatch`), submitted 17:53, pending (Priority). Log `$OPT/build/gpu2.1860007.log`, outputs `$OPT/build/gpu2-1860007/`. It runs:
1. R7 across processes on CUDA with the fixed dump (`r7_threads.sh`, both modes, both routines);
2. the full **fork-cuda** build end to end on the GPU: `aerodyn_utest` and the OLAF r-tests (without the 5MW openfast case and the Python cases), `OFK_MODE=fast` then `strict`, `.outb` copied per mode for comparison with `r6-1859860/off` (`cmp_outb`-style channel comparison).

## Lessons (v1)

See the report of this session; the short version:
1. FMA, not reordering, is what separates fast from strict for `ui_seg`, and the Fortran itself shows the same separation under FMA. A tolerance class should be measured against the baseline's own sensitivity to the same switch (FMA here, OpenMP for R8).
2. GCC vectorizes the fixed-lane loop only if the body is branch-free and uses compare-select `max`; `Kokkos::fmax`, early `return`/`cycle`, and conditional loads each block it.
3. The FMA switch must be per TU **and** the instantiations must live in mode-specific namespaces (ODR/COMDAT).
4. Exp-type kernels: keep the far-field branch; vector `exp` everywhere is slower.
5. `nvcc_wrapper` as the global C++ compiler does not work for OpenFAST; `kokkos_compilation(GLOBAL)` with g++ does.
6. Every test that "passes" by comparing files must also check the files are non-empty (the R7 dumps were empty: `newunit=` returns negative unit numbers).

## History

- **v0** (2026-11-14, same session id, staged outside the fork in `modules/` and `port/`): one `.cpp` with the physics inlined in each shape; strict order only; FMA off on host and on (nvcc default) on device.
  - L0 class B on Kokkos-OpenMP; A1/A2/A5 PASS; `FVW_RunTests` passes through the port.
  - P2 FAIL for `ui_part_nograd` (0.48–0.91×; the Fortran vectorizes with `!$OMP SIMD`), PASS for `ui_seg` (1.02×).
  - GPU job 1858779 (H100 NVL, ran 17:32): `ui_seg` failed its 1e-12 bound on CUDA (1.35e-11 max-normalized, seed 1; 1.27e-11 seed 7; every shape), `ui_part_nograd` 7.95e-16. That was the same FMA effect, then unexplained. Determinism port-vs-port: 0. Speedups vs Fortran-OMP-8 on EPYC: seg LambOseen 17× (10⁴), 113× (10⁵). Log `biotsavart-direct-cuda/gpu_l0.1858779.log`.
- **v1** (this card): moved into the fork, fork CMake integration, one body, strict/fast TUs.

## Cost (E1–E6)

- **E1**: C++ 1177 lines in `fork/modules/kokkos` (body 179, impl 312, entry 292, runtime 111, CMake 109); shim +78 lines; fork CMake +22 lines; tests 1112 lines (Fortran, scripts, CMake).
- **E2/E5**: written by the parent into `tools/cost/ledger.csv`.
- **E6**: CPU jobs 1859729 (8 cores, ~7 min), 1859825 (8 cores), 1859817 (R6, 8 cores); GPU job 1859659 (1 GPU, 2 min 15 s) and GPU job 2. Login-node builds: full fork ~10 min per configuration on 3–4 cores.

## GPU results on A100 (job 1862373, 2026-11-14, sm_80 build on a MIG slice; parent session)
- **Claim 2 (strict GPU vs Fortran):**
  - `ui_seg` worst 3.15e-16 (LambOseen only; the other branches are bitwise);
  - `ui_part_nograd` 8.0e-16 (Exp only).
  - The same as on the H100: only the device `exp` differs.
- **Claim 3 (fast vs strict on GPU):**
  - `ui_part_nograd` 9.4e-15, PASS;
  - `ui_seg` 2.99e-12, FAIL against 1e-12. This is the known FMA effect (decisions.md "FMA in fast mode"), identical to CPU.
- **R7 on GPU across processes: PASS.** Bitwise for strict and fast, both kernels.
- `api_fvw_runtests` was not run on the GPU in this job.
