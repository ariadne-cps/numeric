/***************************************************************************
 *            numeric/number.cpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file numeric/number.cpp
 *  \brief
 */

#include "utility/module.hpp"
#include "foundation/paradigm.hpp"

#include "casts.hpp"
#include "number.hpp"
#include "lower_number.hpp"
#include "upper_number.hpp"
#include "foundation/logical.hpp"
#include "integer.hpp"
#include "decimal.hpp"
#include "dyadic.hpp"
#include "rational.hpp"
#include "real.hpp"
#include "floatdp.hpp"
#include "floatap.hpp"

#include "number_interface.hpp"
#include "number_wrapper.hpp"

#include "float_ball.hpp"
#include "float_bounds.hpp"
#include "float_upper_bound.hpp"
#include "float_lower_bound.hpp"
#include "float_approximation.hpp"

/************ Number *********************************************************/

namespace Ariadne {

// Required by graded dispatch on directed floating-point bounds.
Approximation<FloatDP> pow(Approximation<FloatDP> const& x, Integer const& z) {
    Int n=z.get_si();
    ARIADNE_PRECONDITION(n==z);
    return Approximation<FloatDP>(pow(approx,x._a,n));
}
Approximation<FloatAP> pow(Approximation<FloatAP> const& x, Integer const& z) {
    Int n=z.get_si();
    ARIADNE_PRECONDITION(n==z);
    return Approximation<FloatAP>(pow(approx,x._a,n));
}

ExactDouble::operator ExactNumber() const { return Dyadic(*this).operator ExactNumber(); }
Integer::operator ExactNumber() const { return ExactNumber(new NumberWrapper<Integer>(*this)); }
Dyadic::operator ExactNumber() const { return ExactNumber(new NumberWrapper<Dyadic>(*this)); }
Rational::operator ExactNumber() const { return ExactNumber(new NumberWrapper<Rational>(*this)); }
Real::operator EffectiveNumber() const { return EffectiveNumber(new NumberWrapper<Real>(*this)); }

FloatDPBounds NumberInterface::_get(ValidatedTag, DoublePrecision pr) const { return this->_get(OrderTag(),pr); }
FloatAPBounds NumberInterface::_get(ValidatedTag, ArbitraryPrecision pr) const { return this->_get(OrderTag(),pr); }

FloatDPBall NumberInterface::_get(MetricTag p, DoublePrecision pr) const { return this->_get(p,pr,pr); }
FloatAPBall NumberInterface::_get(MetricTag p, ArbitraryPrecision pr) const { return this->_get(p,pr,pr); }



DyadicBounds::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<DyadicBounds>(*this)); }
DecimalBounds::operator ValidatedNumber() const { return RationalBounds(*this).operator ValidatedNumber(); }
RationalBounds::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<RationalBounds>(*this)); }

template<> FloatDPApproximation::operator ApproximateNumber() const { return ApproximateNumber(new NumberWrapper<FloatDPApproximation>(*this)); }
//template<> FloatDPLowerBound::operator ValidatedLowerNumber() const { return ValidatedLowerNumber(new NumberWrapper<FloatDPLowerBound>(*this)); }
//template<> FloatDPUpperBound::operator ValidatedUpperNumber() const { return ValidatedUpperNumber(new NumberWrapper<FloatDPUpperBound>(*this)); }
template<> FloatDPBounds::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<FloatDPBounds>(*this)); }
//template<> FloatDPBall::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<FloatDPBall>(*this)); }
FloatDP::operator ExactNumber() const { return ExactNumber(new NumberWrapper<FloatDP>(*this)); }

//template<> FloatDPError::operator ValidatedErrorNumber() const { return ValidatedErrorNumber(new NumberWrapper<FloatDPError>(*this)); }
template<> FloatDPError::operator ValidatedErrorNumber() const { return ValidatedErrorNumber(new NumberWrapper<FloatDP>(cast_exact(*this))); }
// Keep the native error conversion equivalent to the existing DP semantics.
template<> FloatAPError::operator ValidatedErrorNumber() const { return ValidatedErrorNumber(new NumberWrapper<FloatAP>(cast_exact(*this))); }


// Keep the native generic-number conversion alongside the other concrete
// specializations: MSVC must emit it in ariadne-numeric.lib for D100.
template<> FloatAPApproximation::operator ApproximateNumber() const { return ApproximateNumber(new NumberWrapper<FloatAPApproximation>(*this)); }
template<> FloatAPBounds::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<FloatAPBounds>(*this)); }
FloatAP::operator ExactNumber() const { return ExactNumber(new NumberWrapper<FloatAP>(*this)); }

ExactNumber cast_exact(ValidatedUpperNumber const& y) { return ExactNumber(y.handle()); }
ExactNumber cast_exact(ValidatedLowerNumber const& y) { return ExactNumber(y.handle()); }




template<> String class_name<ApproximateNumber>() { return "ApproximateNumber"; }
template<> String class_name<ValidatedLowerNumber>() { return "ValidatedLowerNumber"; }
template<> String class_name<ValidatedUpperNumber>() { return "ValidatedUpperNumber"; }
template<> String class_name<ValidatedNumber>() { return "ValidatedNumber"; }
template<> String class_name<EffectiveNumber>() { return "EffectiveNumber"; }
template<> String class_name<ExactNumber>() { return "ExactNumber"; }

} // namespace Ariadne
