// Unit biotsavart-direct: the module constants of FVW_BiotSavart.f90 (pin :12-33) as a Consts<T>
// value. Included only by entry TUs (biotsavart_direct.cpp, induced_velocities.cpp), which are
// compiled with contraction off; the arithmetic here is single operations anyway.
#pragma once

#include "biotsavart_direct_common.hpp"

#include <limits>
#include <type_traits>

namespace ofk_bsd {

// A ReKi literal written as `x_ReKi` in Fortran: rounded once from decimal to T.
template <class T>
constexpr T lit(double d, float f) {
  if constexpr (std::is_same_v<T, float>) return f; else return static_cast<T>(d);
}

// minnorm: FVW_BiotSavart::MINNORM as the Fortran compiler evaluated it (a default-kind literal,
// so it is passed from the Fortran side; audit item 1).
template <class T>
Consts<T> make_consts(T minnorm) {
  Consts<T> k;
  k.precision_ui = std::numeric_limits<T>::epsilon() / T(100);
  k.min_exp_value = T(-10);
  k.part_reg_cut3 = T(8);
  k.part_reg_c2 = lit<T>(1.6, 1.6f);
  k.mindenom = T(0);
  k.minnorm = minnorm;
  k.minnorm2 = minnorm * minnorm;
  k.fourpi_inv = T(0.25) / Kokkos::numbers::pi_v<T>;
  k.lamb_oseen = lit<T>(1.25643, 1.25643f);
  return k;
}

}  // namespace ofk_bsd
