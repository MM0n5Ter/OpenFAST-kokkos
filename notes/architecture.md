# Architecture of the port (v0, 2026-11-14; layout updated to the overlay 2026-11-15)

This note describes the shape every ported unit takes. Three units establish it:
1. the direct kernels;
2. `InducedVelocitiesAll_*` plus lifting-line induction;
3. the particle treecode.

It does not have to be perfect after those three, but it must already carry every high-level idea: API preservation, one source for CPU and GPU, an explicit reproducibility trade-off, the gates, persistent device state and cost tracking.

## 1. Where the code lives
| Path | Role |
|---|---|
| `deps/openfast` | Upstream pin `12fa3fc`, the frozen baseline and oracle. **Never edited.** |
| `src/modules/<mod>/src/*.f90` | Full copies of only the OpenFAST files we change, at their own paths. The module name and public interface stay. Ported bodies become `#ifdef OPENFAST_KOKKOS` shims, and the pin's body stays in `#else`. |
| `src/modules/kokkos/` | `ofk_kernels`: the one C++/Kokkos library for all ported units, so the process has one Kokkos runtime. |
| `src/cmake/ofk_overlay.cmake` | Swaps `src/<path>` for `deps/openfast/<path>` in every target of the pin's unmodified CMake (`-DCMAKE_PROJECT_INCLUDE`). With `-DOPENFAST_KOKKOS=ON` it also adds `ofk_kernels`. See `notes/overlay.md`. The fork it replaced (`fork/`) is retired. |
| `tests/<unit>/` (this repo) | Unit card and L0 tests. The oracle is built from `deps/openfast` with `_ref` renaming, and the port from `src/` where overlaid. |

Build facts from unit 1, which still hold:
- **CUDA:** keep `CMAKE_CXX_COMPILER=g++`; `ofk_kernels` calls `kokkos_compilation(GLOBAL)`. Setting `nvcc_wrapper` for the whole build fails on OpenFAST's own `awae/src/vtk.cpp`.
- **Export set:** `ofk_kernels` is in `install(EXPORT OpenFASTLibraries)`, with a `$<BUILD_INTERFACE:…>` include directory.
- **Linker:** CMake picks the C++ linker automatically because of the static C++ dependency.


## 2. One source, with every difference named
- **One interaction body per physics kernel.** Each `KOKKOS_INLINE_FUNCTION` computes one target–source contribution, templated on the regularization branch. Every backend and every loop shape calls the same body. No backend gets its own copy of the physics.
- **Loop shapes are the only per-backend choice.** They select a policy over the same body:
  - square: `RangePolicy` over targets, serial loop over sources;
  - skinny: `TeamPolicy`, with a team reduction over sources;
  - tree: flat node arrays and an explicit stack.
  - A dispatcher picks the shape from (nTargets, nSources, backend), with thresholds that can be overridden by environment variables.
- **So only three things can make two runs differ, and each is a named switch, never a compiler default:**
  1. summation order;
  2. FMA contraction;
  3. the math library (`exp`, `pow` on the device vs glibc).

## 3. Two modes: strict and fast
| | **strict** | **fast** (default) |
|---|---|---|
| Purpose | Prove the translation, and debug | Production speed |
| Order | Per target, the Fortran's serial order | Fixed blocked order: a constant number of partial sums (SIMD lanes or team reduction), independent of thread count |
| FMA | Off (`-ffp-contract=off`, nvcc `--fmad=false`) | On |
| Data layout | Whatever the strict order needs | Free to change (SoA, sorted, packed) |
| Expected agreement | Bitwise with the Fortran (OPENMP=OFF build) on Serial/OpenMP | ≤ declared tolerance vs strict on the same backend |

- **Both modes are deterministic.** Same input + same backend gives the same bits on every run, so there are no floating-point atomics in sums (criteria R7).
- **Selection:**
  - at runtime, with `OFK_MODE=strict|fast`;
  - no API change;
  - the mode is printed once in the log.
- **Build:** FMA is a per-translation-unit flag, so each unit compiles its body header twice: `<unit>_strict.cpp` (contraction off) and `<unit>_fast.cpp` (contraction on). Both include the same `<unit>_body.hpp`. It is one source, instantiated twice.
- **ODR rule:** each mode's translation unit wraps the body and the loop shapes in a mode-specific namespace (`ofk_<unit>::strict` / `::fast`; the template is `src/modules/kokkos/src/biotsavart_direct_body.hpp`).
  - Without it, the linker merges the inline instantiations, and one mode silently runs the other's arithmetic.
  - A per-TU `-march` has the same hazard. *(Unit 1 lesson.)*
- **FMA can move ill-conditioned formulas far more than reordering does.** The segment cross product `ya*zb - za*yb` cancels for targets near the segment line. There, FMA alone moves `ui_seg` by up to 1.1e-11 relative to Σ|c|, and the pin's own Fortran built with FMA moves by the same amount. Fast mode keeps FMA (user decision, 2026-11-14), and the bound for such formulas is set by the measured `l0_*_ffma` spread (see §4). *(Unit 1 lesson.)*
- **CPU vectorization recipe (GCC 12)** for the fast mode:
  - one fixed set of 8 lanes with a fixed reduction tree, which is thread-count independent for free;
  - Fortran `MAX`/`MIN` written as compare-select, not `Kokkos::fmax`;
  - guards evaluated without branching, accumulating `ok ? o : 0`;
  - column pointers copied into the lambda;
  - `omp simd` on the lane loop;
  - the Fortran's far-field `exp` skip kept (a vector `exp` on every pair was slower);
  - checked with `-fopt-info-vec`.
- **CPU↔GPU bitwise is not a goal** (user decision, 2026-11-14). Getting it would cost GPU speed, through portable math and a fixed order. Strict mode on the GPU is still useful: it isolates math-library differences from ordering differences.

## 4. Equivalence: a chain of claims, each with its own test
| # | Claim | Compared | Test | Bound |
|---|---|---|---|---|
| 1 | The translation is faithful | Fortran (OPENMP=OFF, pin flags) vs strict, Serial/OpenMP | L0 | **0** (class B) |
| 2 | The backend is correct | strict CPU vs strict GPU | L0 | a few ulp, relative to Σ\|contributions\| |
| 3 | The optimizations are safe | strict vs fast, same backend | L0 | ≤ 1e-12 relative to Σ\|contributions\| (class R1), or the declared R2 bound. **FMA kept (user decision):** where FMA moves the formula itself, the bound is 1.1 × the Fortran's own FMA-vs-no-FMA spread on the same inputs |
| 4 | The science is unchanged | fast GPU/CPU whole simulation vs Fortran | r-test, R8 | r-test tolerances, and drift ≤ 2× the Fortran's own envelope. Units 1–2 showed the OMP_1-vs-OMP_N envelope can be exactly 0, so the proposal (pending a user decision) is to use the pin built with FMA in that case. With that envelope, unit 2 passes: 0/542 channels over 2× |

- Bitwise agreement with the Fortran shows only that the C++ does what the Fortran does, on the inputs tested. It is not correctness: the Fortran's own results already change with compiler, flags and OpenMP.
- Claims 3 and 4 are what make a faster, reorganized GPU kernel acceptable.

## 5. Device state: a handle per instance
- **Ownership:** each unit with persistent device data owns an opaque C++ handle. The registry type is not edited.
- **Where the handle lives is decided by the Fortran signatures** (unit 2 lesson):
  - First list, for each entry point, which instance objects arrive and with which intent.
  - Key the handle by the address of an object that **every** entry point receives. In unit 2 that is `m%Sgmt`, and the handle sits in a module-level map guarded by a mutex, capped by least-recently-used eviction.
  - A registry field was rejected for three reasons:
    - `*_CopyMisc` would make two instances share one handle;
    - checkpoints would restore handle ids from another process;
    - the test oracle could no longer share `*_Types` with the port.
- **Cached device data is re-checked by content, not identity:** unit 2 compares the lattice byte for byte against a host shadow on every call.
  - This makes glue copies, rollback, several instances, Python loading and a missing `*_End` cost only speed, never correctness (A7).
  - It needs no version counters threaded through callers.
- **Lifetime:**
  - created in `*_Init`, or lazily on first call;
  - **sized on the first call from the Fortran's maximum allocated extents** (the misc arrays). Doubling is only a fallback: in a growing wake, doubling alone kept reallocating until evaluation 72 of 93 (P6);
  - freed in `*_End`, with a `Kokkos::push_finalize_hook` fallback.
- **Kokkos runtime:** initialized lazily (`ofk_runtime`). If the host program already initialized Kokkos, it is not finalized by us.
- **Residency, phase 1:** Fortran owns the wake state, and the device holds mirrors that are refreshed per call.
- **P6 measurement:** count only allocations with bytes > 0. Kokkos 5.3.0 CUDA emits a 0-byte `TeamScratchMemory` allocation event on every team launch that uses no level-1 scratch. Compare text outputs from line 8 on: the header carries a timestamp.
- **A device algorithm that replaces a host one gets two switches** (unit 3): one that forces it on host backends, so the CPU L0 tests it (`OFK_TC_DEVBUILD=1`), and one that runs both and compares the exact parts bit for bit (`OFK_TC_DEVBUILD_CHECK=1`).
- **Outputs callers can see limit fusion:** the host pack (`Sgmt`, `Part`) still runs because callers can read those arrays. Removing it needs those arrays filled lazily, only when read. *(Corrected 2026-11-17: this note called the pack the GPU path's main Amdahl term. Measured, it is 0.3–0.75 s of a 34–46 s GPU run. The main terms are the host circulation solve (about 20 s) and AeroDyn's wind and tower evaluation at the wind points (about 12 s); see [olaf-unit4.md](olaf-unit4.md).)*
- **Residency, phase 2** (after the three units): the device owns the wake. It needs rollback and copy semantics (criteria A7).

## 6. Small cases and fallback
Below a measured `N_min`, a unit runs in the host execution space (P4). Same source, runtime switch. The r-test cases are small.

## 7. The three units, and what each proves
| Unit | Proves | Done when |
|---|---|---|
| 1. `biotsavart-direct` (`ui_seg`, `ui_part_nograd`). **v1 done 2026-11-14**; FMA policy open. | Shim/oracle/L0 pattern; one body; strict and fast; dispatcher; FMA switch; build integration (`OPENFAST_KOKKOS`) | Claims 1–3 on OpenMP and CUDA; P2 for the fast mode on CPU |
| 2. `induced-velocities` (`InducedVelocitiesAll_Init/Calc/End`, `LiftingLineInducedVelocities`, packing). **CPU gates pass 2026-11-14; GPU pending.** | Handle with persistent buffers (P6); pack fused into the kernel (F1); an OLAF r-test case running end to end on CPU and GPU | Claims 1–4 (R8 on `ad_BAR_OLAF` or the scaled case); A1–A6 on the overlay build |
| 3. `treecode` (`grow_tree_part`, `ui_tree_part`; segment variant deferred). **Done 2026-11-15**; GPU timing sweeps pending. | The irregular pattern: a flat preorder tree with skip pointers (stackless, visits in the Fortran's order); a level-synchronous device build with the Fortran's topology (integer partitions, fixed-tree sums, no sort library, no atomics); the speedup on OLAF's default path | Class R2 against the direct sum (ratio 0.909 = same tree); claims 1–3, R7, R8, P2, P6, A1/A2/A5/A7 pass on CPU and A100 (`tests/treecode/unit.md`); P1/P5 sweep up to 10⁶ particles (CPU done, GPU queued) |

Deferred until after the three: phase-2 residency (wake time-stepping, F3/F4/F9), the host-side circulation solve (F10), and VTK/IO.

**Unit 4 (analysis 2026-11-17, [olaf-unit4.md](olaf-unit4.md)):** the profile picks the host circulation solve (F10-B, class B, Kokkos host space) and the wind-point host work in AeroDyn/InflowWind. Phase-2 residency stays deferred: wake time-stepping is ≤ 0.25% of the runs. Draft card: `tests/circulation-solve/unit.md`.
