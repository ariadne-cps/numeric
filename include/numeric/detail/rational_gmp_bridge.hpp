/***************************************************************************
 *            numeric/detail/rational_gmp_bridge.hpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file numeric/detail/rational_gmp_bridge.hpp
 *  \brief Temporary exact conversions at still-GMP-backed numeric boundaries.
 */

#ifndef ARIADNE_RATIONAL_GMP_BRIDGE_HPP
#define ARIADNE_RATIONAL_GMP_BRIDGE_HPP

#include "numeric/rational.hpp"
#include "numeric/detail/integer_gmp_bridge.hpp"

#include <stdexcept>

namespace Ariadne {

//! Transitional adapter, never stored inside Rational or used by its arithmetic.
inline Rational rational_from_gmp(mpq_srcptr source) {
    return Rational(integer_from_gmp(mpq_numref(source)),integer_from_gmp(mpq_denref(source)));
}

//! GMP arithmetic accepts finite fractions only; never export our special tags.
inline void rational_to_gmp(mpq_ptr destination,Rational const& source) {
    if (!is_finite(source)) {
        throw std::domain_error("Cannot export a nonfinite Rational to GMP");
    }
    integer_to_gmp(mpq_numref(destination),source.numerator());
    integer_to_gmp(mpq_denref(destination),source.denominator());
}

//! Independently owned temporary for legacy GMP/MPFR conversion only.
class GmpRationalBoundary
{
    mpq_t _value;

  public:
    GmpRationalBoundary() { mpq_init(_value); }
    explicit GmpRationalBoundary(Rational const& source) : GmpRationalBoundary() {
        rational_to_gmp(_value,source);
    }
    GmpRationalBoundary(GmpRationalBoundary const&) = delete;
    GmpRationalBoundary& operator=(GmpRationalBoundary const&) = delete;
    ~GmpRationalBoundary() { mpq_clear(_value); }
    mpq_ptr mutable_value() noexcept { return _value; }
    mpq_srcptr value() const noexcept { return _value; }
    Rational to_rational() const { return rational_from_gmp(_value); }
};

} // namespace Ariadne

#endif
