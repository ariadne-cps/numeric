/***************************************************************************
 *            numeric/float_lower_bound.cpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#include "float_lower_bound.hpp"

#include "float_lower_bound.tpl.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"

namespace Ariadne {

LowerBound<Dyadic>::LowerBound(const LowerBound<FloatDP>& x) : LowerBound(Dyadic(x.raw())) {}
LowerBound<Dyadic>::LowerBound(const LowerBound<FloatMP>& x) : LowerBound(Dyadic(x.raw())) {}
LowerBound<FloatDP> LowerBound<Dyadic>::get(DoublePrecision pr) const { return LowerBound<FloatDP>(FloatDP(this->raw(), down, pr)); }

template<>
LowerBound<FloatDP>::LowerBound(const Real&, DoublePrecision);

template class LowerBound<FloatDP>;
template<>
LowerBound<FloatMP>::LowerBound(const ValidatedLowerNumber& y, MultiplePrecision pr)
    : LowerBound<FloatMP>(y.get(pr)) {}

// Real is first enclosed as exact dyadic bounds by Real::get(pr).
// The lower endpoint is then rounded downward by the native FloatMP.
// This template constructor is defined in float_lower_bound.tpl.hpp;
// generic LowerReal computation creates a FloatMP lower endpoint from
// an already rigorously enclosed FloatMP interval.
// Complete the same outward-rounding and precision-conversion constructors

template class LowerBound<FloatMP>;

template<>
String class_name<LowerBound<FloatDP>>() { return "FloatDPLowerBound"; }
template<>
String class_name<LowerBound<FloatMP>>() { return "FloatMPLowerBound"; }

} // namespace Ariadne
