/***************************************************************************
 *            numeric/float_approximation.tpl.hpp
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

/** @file numeric/float_approximation.tpl.hpp
 *  @brief Template implementations for float approximation in Ariadne Numeric.
 */

#ifndef ARIADNE_FLOAT_APPROXIMATION_TPL_HPP
#define ARIADNE_FLOAT_APPROXIMATION_TPL_HPP

#include "float_approximation.hpp"
#include "float_bounds.hpp"
#include "float_error.hpp"
#include "float_lower_bound.hpp"
#include "float_upper_bound.hpp"
#include "number_wrapper.hpp"

namespace Ariadne {

template<class F>
Approximation<F>::Approximation(const LowerBound<F>& x) : Approximation<F>(x.raw()) {}
template<class F>
Approximation<F>::Approximation(const UpperBound<F>& x) : Approximation<F>(x.raw()) {}
template<class F>
Approximation<F>::Approximation(const Bounds<F>& x) : Approximation<F>(x.value_raw()) {}
template<class F>
Approximation<F>::Approximation(const Error<F>& x) : Approximation<F>(x.raw()) {}

template<class F>
Approximation<F>::Approximation(const Real& r, PR pr) : Approximation<F>(r.get(pr)) {}
template<class F>
Approximation<F>::Approximation(const ApproximateNumber& y, PR pr) : Approximation<F>(y.get(pr)) {}
template<class F>
Approximation<F>::operator ApproximateNumber() const { return ApproximateNumber(new NumberWrapper<Approximation<F>>(*this)); }

} // namespace Ariadne

#endif
