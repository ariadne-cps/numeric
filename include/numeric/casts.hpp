/***************************************************************************
 *            numeric/casts.hpp
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

/** @file numeric/casts.hpp
 *  @brief Inclusion header for casting between different kinds of information.
 */

#ifndef ARIADNE_CASTS_HPP
#define ARIADNE_CASTS_HPP

#include "float.decl.hpp"
#include "float_approximation.hpp"
#include "float_ball.hpp"
#include "float_bounds.hpp"
#include "float_error.hpp"
#include "float_lower_bound.hpp"
#include "float_upper_bound.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"
#include "number.decl.hpp"
#include "utility/macros.hpp"

namespace Ariadne {

// Approximation, directed bounds, balls and errors own their raw float as a data member.
// Reinterpreting the wrapper object as the contained F would rely on object-layout equivalence
// that is not guaranteed by C++. Use the wrapper accessors instead.
template<class F>
    requires DerivedFrom<F, FloatDP> or DerivedFrom<F, FloatMP>
inline const F& cast_exact(const F& x) { return x; }
template<class F>
inline const F& cast_exact(const Approximation<F>& x) { return x.raw(); }
template<class F>
inline const F& cast_exact(const LowerBound<F>& x) { return x.raw(); }
template<class F>
inline const F& cast_exact(const UpperBound<F>& x) { return x.raw(); }
template<class F>
inline const F cast_exact(const Bounds<F>& x) { return cast_exact(Approximation<F>(x)); }
template<class F, class FE>
inline const F& cast_exact(const Ball<F, FE>& x) { return x.value_raw(); }
template<class F>
inline const F& cast_exact(const Error<F>& x) { return x.raw(); }
template<class F>
inline const F& cast_exact(const Rounded<F>& x) { return x.raw(); }

template<class F>
inline const Positive<F> cast_exact(const Positive<Bounds<F>>& t) {
    return Positive<F>(cast_exact(static_cast<const Bounds<F>&>(t)));
}
template<class F>
inline const Positive<F> cast_exact(const Positive<UpperBound<F>>& x) {
    return cast_positive(cast_exact(cast_unsigned(x)));
}
template<class F>
inline const Positive<F> cast_exact(const Positive<LowerBound<F>>& x) {
    return cast_positive(cast_exact(cast_unsigned(x)));
}
template<class F>
inline const Positive<F> cast_exact(const Positive<Approximation<F>>& x) {
    return cast_positive(cast_exact(cast_unsigned(x)));
}

inline const RawFloatDP& cast_raw(const RawFloatDP& x) { return x; }
inline const RawFloatDP& cast_raw(const FloatDPApproximation& x) { return x.raw(); }

template<template<class> class T>
inline const T<RawFloatDP>& cast_raw(const T<RawFloatDP>& t) {
    return reinterpret_cast<const T<RawFloatDP>&>(t);
}
template<template<class> class T>
inline const T<RawFloatDP>& cast_raw(const T<FloatDPApproximation>& t) {
    return reinterpret_cast<const T<RawFloatDP>&>(t);
}

// RawFloatDP and FloatDPApproximation are distinct object types, with no inheritance or
// layout guarantee permitting a reference reinterpretation. Construct the wrapper value explicitly.
inline FloatDPApproximation cast_approximate(const RawFloatDP& x) { return FloatDPApproximation(x); }
inline const FloatDPApproximation& cast_approximate(const FloatDPApproximation& x) { return x; }

template<template<class> class T>
inline const T<FloatDPApproximation>& cast_approximate(const T<RawFloatDP>& t) {
    return reinterpret_cast<const T<FloatDPApproximation>&>(t);
}
template<template<class> class T>
inline const T<FloatDPApproximation>& cast_approximate(const T<FloatDPApproximation>& t) {
    return reinterpret_cast<const T<FloatDPApproximation>&>(t);
}

} // namespace Ariadne

#endif
