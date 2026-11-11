/***************************************************************************
 *            numeric/float_upper_bound.tpl.hpp
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

/** @file numeric/float_upper_bound.tpl.hpp
 *  @brief Template implementations for float upper bound in Ariadne Numeric.
 */

#ifndef ARIADNE_FLOAT_UPPER_BOUND_TPL_HPP
#define ARIADNE_FLOAT_UPPER_BOUND_TPL_HPP

#include "float_bounds.hpp"
#include "float_error.hpp"
#include "float_lower_bound.hpp"
#include "float_upper_bound.hpp"
#include "number_wrapper.hpp"
#include "upper_number.hpp"

namespace Ariadne {

template<class F>
UpperBound<F>::UpperBound(const F& x, PR pr) : _u(x, upward, pr) {}
template<class F>
UpperBound<F>::UpperBound(const Bounds<F>& x, PR pr) : _u(x._u, upward, pr) {}
template<class F>
UpperBound<F>::UpperBound(const UpperBound<F>& x, PR pr) : _u(x._u, upward, pr) {}

template<class F>
UpperBound<F>::UpperBound(const Bounds<F>& x) : UpperBound<F>(x.upper_raw()) {}
template<class F>
UpperBound<F>::UpperBound(const Error<F>& x) : UpperBound<F>(x.raw()) {}

template<class F>
UpperBound<F>::UpperBound(const Real& r, PR pr) : UpperBound(r.get(pr)) {}
template<class F>
UpperBound<F>::UpperBound(const ValidatedUpperNumber& y, PR pr) : UpperBound(y.get(pr)) {}
template<class F>
UpperBound<F>::operator ValidatedUpperNumber() const {
    return ValidatedUpperNumber(Handle<const NumberInterface>(new NumberWrapper<UpperBound<F>>(*this)));
}
template<class F>
auto UpperBound<F>::generic() const -> GenericType { return this->operator GenericType(); }

} // namespace Ariadne

#endif
