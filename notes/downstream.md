# Downstream: who uses OpenFAST as a library

The API gate means **zero source edits for every consumer below**. Each one freezes a part of the API.

| Consumer | How it couples | Frozen surface | Touches the OLAF unit? |
|---|---|---|---|
| **OpenFAST glue** (`openfast`, `openfast-library`) | Fortran calls through the framework API | `AD_*` → `FVW_*` signatures, registry types | Yes: r-test `of_regression` OLAF cases |
| **AeroDyn driver** (`aerodyn_driver`) | Fortran | same | Yes: 9 `ad_*_OLAF` r-test cases |
| **pyOpenFAST** (`glue-codes/python`, `openfast_io`) | `ctypes` on the C bindings (`libaerodyn_inflow_c_binding.so`, IfW, HD, SeaState, MoorDyn) and on `FAST_Library` | Exported C symbols and their argument layouts; shared-library loading | Yes: `py_ad_B1n2_OLAF`, `EllipticalWing_OLAF_py`, `HelicalWake_OLAF_py` |
| **AMR-Wind / Nalu-Wind (ExaWind)** | `glue-codes/openfast-cpp` → `FAST_Library` C API (actuator line, CFD coupling) | `FAST_Library` symbols; restart files | Only if a coupled turbine runs OLAF (unusual: the CFD resolves the wake) |
| **Simulink**, **LabVIEW wavetank** | `FAST_Library` / `WaveTank.f90` C API | Same symbols. Simulink needs a static build (`BUILD_SHARED_LIBS=OFF`). | Indirectly |
| **ROSCO** and other controllers | Bladed-style controller DLL loaded by ServoDyn | `avrSWAP` array layout, DLL calling convention | No |
| **WEIS / WISDEM**, users' scripts | Write input decks and read `.out`/`.outb` through `openfast_io` | Input file formats, output channel names and units, `.outb` binary format | Yes: OLAF output channels |
| **FAST.Farm** | Runs one OpenFAST instance per turbine, plus AWAE/WakeDynamics | Framework API; several instances in one process | *(verify)* whether OLAF per turbine is supported. If it is, the port must allow several instances. |
| **Visualization** (ParaView users, our animation) | VTK files from `WrVTK`, including OLAF's wake (`FVW_IO`) | VTK file names, arrays, structure | Yes |

**Consequences for the port**
1. The ported code may sit **inside a shared library loaded from Python** (`libaerodyn_inflow_c_binding.so`). So Kokkos must work from a `dlopen`ed library and initialize itself; see [build-kokkos-fortran.md](build-kokkos-fortran.md).
2. **Several instances per process** (FAST.Farm, multi-rotor `ad_QuadRotor_OLAF`) mean no per-instance state in C++ globals.
3. A static-library build must keep working (Simulink).
4. Output headers (channel names, units) and VTK structure are checked byte for byte. Values are checked at r-test tolerances.
