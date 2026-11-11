/***************************************************************************
 *            numeric/double.hpp
 *
 *  Copyright  2008-22  Pieter Collins
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

/** @file numeric/double.hpp
 *  @brief Correctly rounded operations on double-precision numbers.
 */

#ifndef ARIADNE_DOUBLE_HPP
#define ARIADNE_DOUBLE_HPP

namespace Ariadne {

// Correctly rounded functions
/** @brief Return the additive zero with the operand characteristics with the specified rounding mode. */
double nul_rnd(double x);
/** @brief Return the unary positive value with the specified rounding mode. */
double pos_rnd(double x);
/** @brief Return the additive inverse with the specified rounding mode. */
double neg_rnd(double x);
/** @brief Divide the value by two with the specified rounding mode. */
double hlf_rnd(double x);
/** @brief Return the square with the specified rounding mode. */
double sqr_rnd(double x);
/** @brief Return the multiplicative inverse with the specified rounding mode. */
double rec_rnd(double x);
/** @brief Return the sum of the operands with the specified rounding mode. */
double add_rnd(double x1, double x2);
/** @brief Return the difference of the operands with the specified rounding mode. */
double sub_rnd(double x1, double x2);
/** @brief Return the product of the operands with the specified rounding mode. */
double mul_rnd(double x1, double x2);
/** @brief Return the quotient of the operands with the specified rounding mode. */
double div_rnd(double x1, double x2);
/** @brief Compute the product of the first two operands plus the third with the specified rounding mode. */
double fma_rnd(double x1, double x2, double x3);
/** @brief Raise the value to the specified power with the specified rounding mode. */
double pow_rnd(double x, Nat n);
/** @brief Raise the value to the specified power with the specified rounding mode. */
double pow_rnd(double x, int n);
/** @brief Return the square root with the specified rounding mode. */
double sqrt_rnd(double x);
/** @brief Return the exponential with the specified rounding mode. */
double exp_rnd(double x);
/** @brief Return the natural logarithm with the specified rounding mode. */
double log_rnd(double x);
/** @brief Return the sine with the specified rounding mode. */
double sin_rnd(double x);
/** @brief Return the cosine with the specified rounding mode. */
double cos_rnd(double x);
/** @brief Return the tangent with the specified rounding mode. */
double tan_rnd(double x);
/** @brief Return the hyperbolic tangent with the specified rounding mode. */
double tanh_rnd(double x);
/** @brief Return the inverse sine with the specified rounding mode. */
double asin_rnd(double x);
/** @brief Return the inverse cosine with the specified rounding mode. */
double acos_rnd(double x);
/** @brief Return the inverse tangent with the specified rounding mode. */
double atan_rnd(double x);
/** @brief Compute the negated reciprocal with the specified rounding mode. */
double neg_rec_rnd(double x);
/** @brief Evaluate the inverse tangent series with the requested rounding. */
double atan_rnd_series(double x);
/** @brief Return pi at the requested precision with the specified rounding mode. */
double pi_rnd();

/** @brief Return the absolute value with the specified rounding mode. */
double abs_rnd(double x);
/** @brief Return the maximum of the operands with the specified rounding mode. */
double max_rnd(double x1, double x2);
/** @brief Return the minimum of the operands with the specified rounding mode. */
double min_rnd(double x1, double x2);

/** @brief Return the exponential using the local rounded-arithmetic implementation. */
double texp(double x);

/** @brief Return pi at the requested precision with opposite directed rounding. */
double pi_opp();
/** @brief Return the sum of the operands with opposite directed rounding. */
double add_opp(double x, double y);
/** @brief Return the difference of the operands with opposite directed rounding. */
double sub_opp(double x, double y);
/** @brief Return the product of the operands with opposite directed rounding. */
double mul_opp(double x, double y);
/** @brief Return the quotient of the operands with opposite directed rounding. */
double div_opp(double x, double y);
/** @brief Compute the negated reciprocal with opposite directed rounding. */
double neg_rec_opp(double x);

} // namespace Ariadne

#endif
