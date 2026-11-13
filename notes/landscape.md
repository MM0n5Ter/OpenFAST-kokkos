# Landscape

The general literature on agentic Fortran → C++/Kokkos translation is the same as for our NWChem work: Gupta et al. (kernel level); NICAM, FESOM2 and ORB5X (application ports, standalone); Lahlou et al. (dependency-ordered, Fortran → JAX); Fortran2CPP; ParEval-Repo; RustRepoTrans. The table and the lessons we reuse are in `NWChem-kokkos/notes/landscape.md`. This note covers what is specific to OpenFAST.

Sources: `../investigation/openfast.md` in the parent workspace (2026-10-28) and the upstream issues and PRs it cites.

## GPU and performance work on OpenFAST / OLAF

| Work | What it is | What it means for us |
|---|---|---|
| Upstream issue [#2852](https://github.com/OpenFAST/openfast/issues/2852) "GPU Acceleration for OLAF's Biot-Savart Law?" (2025-06-09, open) | Maintainer: "GPU acceleration of OLAF is not currently supported… NREL is not currently funded to work on that." | Demand acknowledged, nobody upstream working on it |
| Upstream CPU optimization of OLAF (NLR): [#3430](https://github.com/OpenFAST/openfast/pull/3430), [#3457](https://github.com/OpenFAST/openfast/pull/3457), [#3464](https://github.com/OpenFAST/openfast/pull/3464) (merged Sep 2026); #3483, #3499 open | Treecode regularization fix, compact particle core (`RegFunctionPart`), SIMD near-field, leaf-multipole MAC. "Just over 3×" overall. Vectorizing near-field P2P alone cut wall-clock by 41%. | The kernels are a **moving target** and the CPU baseline is strong. We pin `dev` ([decisions.md](decisions.md) #1) and compare against it. Near-field P2P is ≥41% of runtime. |
| Open correctness issues [#3401](https://github.com/OpenFAST/openfast/issues/3401), [#3402](https://github.com/OpenFAST/openfast/issues/3402) | OLAF results differ with and without OpenMP, and are not deterministic | Baselines use `OMP_NUM_THREADS=1` |
| [Sunkai233/openfast-cuda](https://github.com/Sunkai233/openfast-cuda) (2026-07, third party, unreviewed) | Standalone CUDA re-implementations of OpenFAST kernels, including `ui_seg` (`adk_fvw_gpu.cu`). Claims bit-identical results and 17–319× speedups. No treecode, and no integration glue in the repo. | Porting only the direct `ui_seg` kernel would be largely redundant. **The treecode and the in-place, API-preserving method are where the novelty is.** |
| QBlade (TU Berlin) | Lifting-line free-vortex wake on GPU with OpenCL, "2–3 orders of magnitude faster than single-core CPU" ([WES 10, 3045, 2025](https://wes.copernicus.org/articles/10/3045/2025/)) | GPU free-vortex wakes are known in wind energy. Our contribution is the verified migration of a production code, not the algorithm. |
| Older literature | GPU free-wake panel code for NREL Phase VI (J. Phys. Conf. Ser. 524 012100, 2014). GPU treecodes and FMM for vortex particles (Yokota & Barba 2011; Hu, Gumerov & Duraiswami 2013; Stock & Gharakhani). | Design references for the Kokkos treecode: linear octree, stackless traversal |
| ExaWind (AMR-Wind / Nalu-Wind) | GPU CFD coupled to OpenFAST through the C++ API (actuator line). OpenFAST itself stays on CPU. | A downstream consumer, not overlap. See [downstream.md](downstream.md). |

## Where we fit
1. Kernel-level automation works (Gupta), and a direct `ui_seg` CUDA port already exists (openfast-cuda).
2. Nobody has ported OLAF's **treecode**, or done so **in place** behind OpenFAST's interfaces with its own regression suite as the judge.
3. **Our gap:** an *API-preserving, incremental, repository-level* Fortran → Kokkos migration of a code used as a library. It starts at the algorithm that dominates runtime, and every step is gated by API, Regression and Performance checks.
4. Scope claims: OLAF is a mid-fidelity option, and most OpenFAST and FAST.Farm runs use BEM. We say "OLAF free-vortex wake on GPU", not "OpenFAST on GPU".
