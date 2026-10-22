/***************************************************************************
 *            numeric/float_approximation.cpp
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

#include "float_approximation.hpp"
#include "float_approximation.tpl.hpp"

#include "floatdp.hpp"
#include "floatmp.hpp"
#include "floatap.hpp"

namespace Ariadne {

template<> Approximation<FloatDP>::Approximation(Real const&, DoublePrecision);
template<> Approximation<FloatDP>::operator ApproximateNumber() const;
template<> Approximation<FloatMP>::operator ApproximateNumber() const;

template class Approximation<FloatDP>;
template class Approximation<FloatMP>;

// Native constructors are implemented in float_approximation.tpl.hpp.
// Instantiate every constructor used by the complete FloatAP approximation
// test suite, rather than suppressing its semantic regression coverage.
template Approximation<FloatAP>::Approximation(Real const&, ArbitraryPrecision);
template Approximation<FloatAP>::Approximation(ApproximateNumber const&, ArbitraryPrecision);
template Approximation<FloatAP>::Approximation(LowerBound<FloatAP> const&);
template Approximation<FloatAP>::Approximation(UpperBound<FloatAP> const&);
template Approximation<FloatAP>::Approximation(Bounds<FloatAP> const&);
template Approximation<FloatAP>::Approximation(Error<FloatAP> const&);

template<> String class_name<Approximation<FloatDP>>() { return "FloatDPApproximation"; }
template<> String class_name<Approximation<FloatMP>>() { return "FloatMPApproximation"; }
template<> String class_name<Approximation<FloatAP>>() { return "FloatAPApproximation"; }

// Generate the native generic-number conversion after class_name is defined:
// NumberWrapper<Approximation<FloatAP>> needs the specialization at link time.
template Approximation<FloatAP>::operator ApproximateNumber() const;

} // namespace Ariadne
