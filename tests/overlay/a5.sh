#!/bin/bash
# A5 for the overlay: an overlay build with OPENFAST_KOKKOS=OFF must be the pin's build.
#   tests/overlay/a5.sh <pin build> <overlay OFF build> [target dir, default modules/aerodyn]
# Compares, for every object of aerodynlib: the compile command (flags.make + per-file flags),
# the full disassembly (all code sections, relocations, addresses stripped), and every .mod both trees produced (md5).
set -u
P=$1; O=$2; T=${3:-modules/aerodyn}
module load gcc/12.3.0 >/dev/null 2>&1 || true
tmp=$(mktemp -d); trap 'rm -rf $tmp' EXIT
echo "=== .mod md5 (every .mod present in both ftnmods/)"
same=0; diff=0
for m in $P/ftnmods/*.mod; do
  n=$(basename $m); [ -f $O/ftnmods/$n ] || continue
  if cmp -s $m $O/ftnmods/$n; then same=$((same+1)); else diff=$((diff+1)); echo "  DIFF $n"; fi
done
echo ".mod: $same identical, $diff different"

echo "=== aerodynlib objects: full disassembly (every code section, COMDAT included, with relocations)"
pd=$P/$T/CMakeFiles/aerodynlib.dir; od=$O/$T/CMakeFiles/aerodynlib.dir
dis() { objdump -d -r --no-show-raw-insn "$1" | tail -n +4 | sed -E 's/^ *[0-9a-f]+:\t//; s/<[^>]*\+0x[0-9a-f]+>//g'; }
for po in $(cd $pd && find . -name '*.o' | sort); do
  b=$(basename $po)
  oo=$(find $od -name "$b" | head -1)
  if [ -z "$oo" ]; then echo "  MISSING in overlay build: $b"; continue; fi
  dis $pd/$po > $tmp/a.s; dis $oo > $tmp/b.s
  if cmp -s $tmp/a.s $tmp/b.s; then
    echo "  $b: identical ($(wc -l < $tmp/a.s) lines)"
  else
    nd=$(diff $tmp/a.s $tmp/b.s | grep -c '^<')
    echo "  $b: differs in $nd of $(wc -l < $tmp/a.s) lines; every differing line:"
    diff $tmp/a.s $tmp/b.s | grep '^[<>]' | sed 's/^/      /'
  fi
done
echo "=== aerodynlib compile flags"
diff <(grep -v '^#' $pd/flags.make | sed "s#$P#<build>#g") <(grep -v '^#' $od/flags.make | sed "s#$O#<build>#g") && echo "  flags.make identical"
echo "  per-file flags (build.make: what follows \$(Fortran_FLAGS) in each compile command):"
pf() { grep -E "gfortran .*-c [^ ]*/$2\.f90 " $1/build.make | sed -E 's/.*\$\(Fortran_FLAGS\)(.*) -c .*/[\1 ]/'; }
for f in $(cd $pd && find . -name '*.f90.o' -printf '%f\n' | sed 's/\.f90\.o$//' | sort); do
  a=$(pf $pd $f); b=$(pf $od $f)
  [ "$a" = "$b" ] && echo "    $f: identical $a" || echo "    $f: DIFFER pin $a overlay $b"
done
