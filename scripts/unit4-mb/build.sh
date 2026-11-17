#!/bin/bash
# Unit-4 analysis microbenchmarks (notes/olaf-unit4.md §2). Builds into <outdir>; run there.
#   f10mb [nSpan] [RegFunction]  : one circulation iteration, pin ui_quad_n1/ui_seg_11 vs F10-B replay vs F10-R1 mat-vec
#   ieeemb                       : UniformField_GetVel with and without 'use, intrinsic :: ieee_arithmetic'
# bs.f90 is extracted verbatim from the pin (FVW_BiotSavart.f90 lines 12-33, 54-128, 512-550) at build time.
set -euo pipefail
out=${1:?outdir}; mkdir -p "$out"; H=$(cd "$(dirname "$0")" && pwd)
P=$H/../../deps/openfast/modules/aerodyn/src/FVW_BiotSavart.f90
module load gcc/12.3.0 >/dev/null 2>&1 || true
{ echo "module bs"; echo "   implicit none"; echo "   integer, parameter :: ReKi = 8, IntKi = 4"; sed -n 12,33p $P; echo "contains"
  sed -n 54,128p $P; sed -n 512,550p $P; cat $H/seg_geom.f90; } > $out/bs.f90
F="-O3 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -cpp -fstack-reuse=none"   # the pin's Release flags
cd $out
gfortran $F -c bs.f90 && gfortran $F -c $H/f10mb_main.f90 -o f10mb_main.o && gfortran -o f10mb f10mb_main.o bs.o
gfortran $F -fopenmp -c $H/ieeemb.f90 && gfortran $F -fopenmp -c $H/ieeemain.f90 && gfortran -fopenmp -o ieeemb ieeemain.o ieeemb.o
echo "built $out/f10mb $out/ieeemb"
