/***************************************************************************
 *            numeric/float_literals.hpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Numeric.
 *
 *  Ariadne Numeric is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Numeric is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric.  If not, see <https://www.gnu.org/licenses/>.
 */

/** @file numeric/float_literals.hpp
 *  @brief Inclusion header for floating-point extended literals.
 */

#ifndef ARIADNE_FLOAT_LITERALS_HPP
#define ARIADNE_FLOAT_LITERALS_HPP

#include "float.decl.hpp"

namespace Ariadne {

/** @brief Construct the typed value from the literal. */
FloatDPError operator""_error(long double lx);
/** @brief Construct the typed value from the literal. */
FloatDP operator""_exact(long double lx);
/** @brief Construct the typed value from the literal. */
FloatDPBall operator""_near(long double lx);
/** @brief Construct the typed value from the literal. */
FloatDPUpperBound operator""_upper(long double lx);
/** @brief Construct the typed value from the literal. */
FloatDPLowerBound operator""_lower(long double lx);
/** @brief Construct the typed value from the literal. */
FloatDPApproximation operator""_approx(long double lx);

} // namespace Ariadne

#endif
