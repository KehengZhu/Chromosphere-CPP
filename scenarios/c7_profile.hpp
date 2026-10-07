/*!
 * @file scenarios/c7_profile.hpp
 * @brief The tabulated Model C7 atmosphere (Avrett & Loeser 2008, Table 26) and
 *        its two interpolants.
 * @ingroup scenarios
 *
 * The C7 profile library shared by the release `model_column` IC (which samples
 * its temperature and ionization fraction here) and the legacy model_c7
 * scenario family (which reads the MODEL_C7 rows directly). It depends on the
 * Grid scalar types only — no solver, no physics-rate library.
 */

#pragma once

#include "core/chromosphere.hpp"

namespace chromosphere {

/** @addtogroup scenarios
 *  @{
 */

/// Number of MODEL_C7 rows forming the standard model_c7 column: chromosphere
/// plus lower transition region, from h = 1003 km to h ≈ 2153 km (T ≈ 23 kK).
constexpr int N_C7_BASE = 32;

/// Total number of MODEL_C7 rows: the standard column extended through the upper
/// transition region and the corona to the top of Table 26 (h ≈ 68 Mm,
/// T ≈ 1.59 MK).
constexpr int N_C7 = 70;

/// Model C7 chromosphere → corona rows (Avrett & Loeser 2008, Table 26), from
/// h = 1003 km upward. Columns: height [km], electron density n_e [m^-3],
/// neutral-hydrogen density n_HI [m^-3], temperature T [K]. Rows
/// 0..N_C7_BASE-1 are the standard model_c7 column; rows N_C7_BASE..N_C7-1 are
/// the upper-TR and coronal extension. Height and temperature increase
/// monotonically down the table.
extern const Real MODEL_C7[N_C7][4];

/// Sample the complete Model C7 atmosphere (Avrett & Loeser 2008 Table 26,
/// photosphere h ≈ −10 km → corona h ≈ 68 Mm) at height `h_km` [km]:
/// temperature `T` [K] linearly, electron density `n_e` [m^-3] and
/// neutral-hydrogen density `n_HI` [m^-3] log-linearly. Concatenates the
/// photosphere rows (h < 1003 km, absent from MODEL_C7) with MODEL_C7; clamps to
/// the table ends. model_column samples its ionization fraction n_e/(n_e+n_HI)
/// here.
void c7_full_profile(Real h_km, Real& T, Real& n_e, Real& n_HI);

/// Shape-preserving C1 (monotone PCHIP) interpolation of the complete Model C7
/// temperature [K] at height `h_km` [km]. Unlike c7_full_profile's
/// piecewise-linear temperature, this has a continuous dT/dh and therefore does
/// not seed knot-scale jumps in conductive heat flux. Used only by the
/// Gamma1-table (release and two-temperature) model_column initialization; the
/// fixed-gamma scenarios keep their existing interpolation exactly.
double c7_full_temperature_pchip(double h_km);

/** @} */

} // namespace chromosphere
