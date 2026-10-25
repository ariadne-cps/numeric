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
#include "floatmp.hpp"

#include "concepts.hpp"

namespace Ariadne {

FloatDP midpoint(Bounds<FloatDP> const& x) { return x.value(); } // DEPRECATED
FloatMP midpoint(Bounds<FloatMP> const& x) { return x.value(); } // DEPRECATED


template<> Bounds<FloatDP>::Bounds(Real const&, DoublePrecision);
template<> Bounds<FloatMP>::Bounds(Real const& value, MultiplePrecision pr)
    : Bounds<FloatMP>(value.get(pr)) { }

// NumberWrapper calls this when converting generic validated values
// to native outward-rounded intervals.
template<> Bounds<FloatDP>::operator ValidatedNumber() const;
// This specialization is defined in number.cpp. Declare it before the
// explicit class instantiation below, preventing a second definition on MSVC.
template<> Bounds<FloatMP>::operator ValidatedNumber() const;

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
template OutputStream& Operations<Bounds<FloatDP>>::_write(OutputStream&,Bounds<FloatDP> const&);
template InputStream& Operations<Bounds<FloatDP>>::_read(InputStream&, Bounds<FloatDP>&);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_sin(Bounds<FloatMP> const&);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_cos(Bounds<FloatMP> const&);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_tan(Bounds<FloatMP> const&);
template OutputStream& Operations<Bounds<FloatMP>>::_write(OutputStream&,Bounds<FloatMP> const&);

template class Bounds<FloatMP>;
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_pi(MultiplePrecision);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_trunc(Bounds<FloatMP> const&);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_trunc(Bounds<FloatMP> const&, Nat);
template Integer Operations<Bounds<FloatMP>>::_cast_integer(Bounds<FloatMP> const&);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_add(Bounds<FloatMP> const&, FloatMP const&);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_add(FloatMP const&, Bounds<FloatMP> const&);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_sub(Bounds<FloatMP> const&, FloatMP const&);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_sub(FloatMP const&, Bounds<FloatMP> const&);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_mul(Bounds<FloatMP> const&, FloatMP const&);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_mul(FloatMP const&, Bounds<FloatMP> const&);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_div(Bounds<FloatMP> const&, FloatMP const&);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_div(FloatMP const&, Bounds<FloatMP> const&);
template InputStream& Operations<Bounds<FloatMP>>::_read(InputStream&, Bounds<FloatMP>&);


template<> String class_name<Bounds<FloatDP>>() { return "FloatDPBounds"; }
template<> String class_name<Bounds<FloatMP>>() { return "FloatMPBounds"; }

static_assert(TranscendentalField<Bounds<FloatDP>>);
static_assert(TranscendentalField<Bounds<FloatMP>>);

static_assert(OrderedLattice<Bounds<FloatDP>>);
static_assert(OrderedLattice<Bounds<FloatMP>>);

static_assert(WidenTranscendentalField<FloatDP>);
static_assert(WidenTranscendentalField<FloatMP>);

} // namespace Ariadne
