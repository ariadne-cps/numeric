/***************************************************************************
 *            numeric/float_bounds.cpp
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

#include "float_bounds.hpp"
#include "float_bounds.inl.hpp"
#include "float_bounds.tpl.hpp"

#include "floatdp.hpp"
#include "floatap.hpp"

#include "concepts.hpp"

namespace Ariadne {

FloatDP midpoint(Bounds<FloatDP> const& x) { return x.value(); } // DEPRECATED
FloatAP midpoint(Bounds<FloatAP> const& x) { return x.value(); } // DEPRECATED


template<> Bounds<FloatDP>::Bounds(Real const&, DoublePrecision);
template<> Bounds<FloatAP>::Bounds(Real const& value, ArbitraryPrecision pr)
    : Bounds<FloatAP>(value.get(pr)) { }

// NumberWrapper calls this when converting generic validated values
// to native outward-rounded intervals.
template<> Bounds<FloatDP>::operator ValidatedNumber() const;
// This specialization is defined in number.cpp. Declare it before the
// explicit class instantiation below, preventing a second definition on MSVC.
template<> Bounds<FloatAP>::operator ValidatedNumber() const;

template class Bounds<FloatDP>;
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_pi(DoublePrecision);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_sin(Bounds<FloatDP> const&);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_cos(Bounds<FloatDP> const&);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_tan(Bounds<FloatDP> const&);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_trunc(Bounds<FloatDP> const&);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_trunc(Bounds<FloatDP> const&, Nat);
template Integer Operations<Bounds<FloatDP>>::_cast_integer(Bounds<FloatDP> const&);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_add(Bounds<FloatDP> const&, FloatDP const&);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_add(FloatDP const&, Bounds<FloatDP> const&);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_sub(Bounds<FloatDP> const&, FloatDP const&);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_sub(FloatDP const&, Bounds<FloatDP> const&);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_mul(Bounds<FloatDP> const&, FloatDP const&);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_mul(FloatDP const&, Bounds<FloatDP> const&);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_div(Bounds<FloatDP> const&, FloatDP const&);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_div(FloatDP const&, Bounds<FloatDP> const&);
template InputStream& Operations<Bounds<FloatDP>>::_read(InputStream&, Bounds<FloatDP>&);
template Bounds<FloatAP> Operations<Bounds<FloatAP>>::_sin(Bounds<FloatAP> const&);
template Bounds<FloatAP> Operations<Bounds<FloatAP>>::_cos(Bounds<FloatAP> const&);
template Bounds<FloatAP> Operations<Bounds<FloatAP>>::_tan(Bounds<FloatAP> const&);
template OutputStream& Operations<Bounds<FloatAP>>::_write(OutputStream&,Bounds<FloatAP> const&);

template class Bounds<FloatAP>;
template Bounds<FloatAP> Operations<Bounds<FloatAP>>::_pi(ArbitraryPrecision);
template Bounds<FloatAP> Operations<Bounds<FloatAP>>::_trunc(Bounds<FloatAP> const&);
template Bounds<FloatAP> Operations<Bounds<FloatAP>>::_trunc(Bounds<FloatAP> const&, Nat);
template Integer Operations<Bounds<FloatAP>>::_cast_integer(Bounds<FloatAP> const&);
template Bounds<FloatAP> Operations<Bounds<FloatAP>>::_add(Bounds<FloatAP> const&, FloatAP const&);
template Bounds<FloatAP> Operations<Bounds<FloatAP>>::_add(FloatAP const&, Bounds<FloatAP> const&);
template Bounds<FloatAP> Operations<Bounds<FloatAP>>::_sub(Bounds<FloatAP> const&, FloatAP const&);
template Bounds<FloatAP> Operations<Bounds<FloatAP>>::_sub(FloatAP const&, Bounds<FloatAP> const&);
template Bounds<FloatAP> Operations<Bounds<FloatAP>>::_mul(Bounds<FloatAP> const&, FloatAP const&);
template Bounds<FloatAP> Operations<Bounds<FloatAP>>::_mul(FloatAP const&, Bounds<FloatAP> const&);
template Bounds<FloatAP> Operations<Bounds<FloatAP>>::_div(Bounds<FloatAP> const&, FloatAP const&);
template Bounds<FloatAP> Operations<Bounds<FloatAP>>::_div(FloatAP const&, Bounds<FloatAP> const&);
template InputStream& Operations<Bounds<FloatAP>>::_read(InputStream&, Bounds<FloatAP>&);


template<> String class_name<Bounds<FloatDP>>() { return "FloatDPBounds"; }
template<> String class_name<Bounds<FloatAP>>() { return "FloatAPBounds"; }

static_assert(TranscendentalField<Bounds<FloatDP>>);
static_assert(TranscendentalField<Bounds<FloatAP>>);

static_assert(OrderedLattice<Bounds<FloatDP>>);
static_assert(OrderedLattice<Bounds<FloatAP>>);

static_assert(WidenTranscendentalField<FloatDP>);
static_assert(WidenTranscendentalField<FloatAP>);

} // namespace Ariadne
