/***************************************************************************
 *            numeric/float_upper_bound.cpp
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

#include "float_upper_bound.hpp"
#include "float_upper_bound.tpl.hpp"

#include "floatdp.hpp"
#include "floatmp.hpp"
#include "floatap.hpp"

namespace Ariadne {

UpperBound<Dyadic>::UpperBound(UpperBound<FloatDP> const& x) : UpperBound(Dyadic(x.raw())) { }
UpperBound<Dyadic>::UpperBound(UpperBound<FloatMP> const& x) : UpperBound(Dyadic(Rational(x.raw()))) { }
UpperBound<Dyadic>::UpperBound(UpperBound<FloatAP> const& x) : UpperBound(Dyadic(x.raw())) { }
UpperBound<FloatDP> UpperBound<Dyadic>::get(DoublePrecision pr) const { return UpperBound<FloatDP>(FloatDP(this->raw(),up,pr)); }
UpperBound<FloatMP> UpperBound<Dyadic>::get(MultiplePrecision pr) const { return UpperBound<FloatMP>(FloatMP(this->raw(),up,pr)); }

template<> UpperBound<FloatDP>::UpperBound(Real const&, DoublePrecision);

template class UpperBound<FloatDP>;
template<> UpperBound<FloatAP>::UpperBound(ValidatedUpperNumber const& y, ArbitraryPrecision pr)
    : UpperBound<FloatAP>(y.get(pr)) { }

// Upper endpoint construction mirrors LowerBound but rounds upward.
template UpperBound<FloatAP>::UpperBound(Real const&, ArbitraryPrecision);
// Mirror the lower endpoint constructor. The interval itself already
// encloses the real value and upper_raw() is its upper endpoint.
template UpperBound<FloatAP>::UpperBound(Bounds<FloatAP> const&);
 // Match FloatMP's native precision-conversion surface.
template UpperBound<FloatAP>::UpperBound(FloatAP const&, ArbitraryPrecision);
template UpperBound<FloatAP>::UpperBound(Bounds<FloatAP> const&, ArbitraryPrecision);
template UpperBound<FloatAP>::UpperBound(UpperBound<FloatAP> const&, ArbitraryPrecision);


template class UpperBound<FloatMP>;
// Export the complete native directed-bound API just as for FloatMP.
template class UpperBound<FloatAP>;

template<> String class_name<UpperBound<FloatDP>>() { return "FloatDPUpperBound"; }
template<> String class_name<UpperBound<FloatMP>>() { return "FloatMPUpperBound"; }
template<> String class_name<UpperBound<FloatAP>>() { return "FloatAPUpperBound"; }

} // namespace Ariadne
