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

#include "concepts.hpp"
#include "float_bounds.inl.hpp"
#include "float_bounds.tpl.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"

namespace Ariadne {

FloatDP midpoint(const Bounds<FloatDP> &x) { return x.value(); } // DEPRECATED
FloatMP midpoint(const Bounds<FloatMP> &x) { return x.value(); } // DEPRECATED

template<>
Bounds<FloatDP>::Bounds(const Real &, DoublePrecision);
template<>
Bounds<FloatMP>::Bounds(const Real &value, MultiplePrecision pr)
    : Bounds<FloatMP>(value.get(pr)) {}

// NumberWrapper calls this when converting generic validated values
// to native outward-rounded intervals.
template<>
Bounds<FloatDP>::operator ValidatedNumber() const;
// This specialization is defined in number.cpp. Declare it before the
// explicit class instantiation below, preventing a second definition on MSVC.
template<>
Bounds<FloatMP>::operator ValidatedNumber() const;

template class Bounds<FloatDP>;
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_pi(DoublePrecision);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_sin(const Bounds<FloatDP> &);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_cos(const Bounds<FloatDP> &);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_tan(const Bounds<FloatDP> &);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_trunc(const Bounds<FloatDP> &);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_trunc(const Bounds<FloatDP> &, Nat);
template Integer Operations<Bounds<FloatDP>>::_cast_integer(const Bounds<FloatDP> &);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_add(const Bounds<FloatDP> &, const FloatDP &);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_add(const FloatDP &, const Bounds<FloatDP> &);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_sub(const Bounds<FloatDP> &, const FloatDP &);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_sub(const FloatDP &, const Bounds<FloatDP> &);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_mul(const Bounds<FloatDP> &, const FloatDP &);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_mul(const FloatDP &, const Bounds<FloatDP> &);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_div(const Bounds<FloatDP> &, const FloatDP &);
template Bounds<FloatDP> Operations<Bounds<FloatDP>>::_div(const FloatDP &, const Bounds<FloatDP> &);
template OutputStream &Operations<Bounds<FloatDP>>::_write(OutputStream &, const Bounds<FloatDP> &);
template InputStream &Operations<Bounds<FloatDP>>::_read(InputStream &, Bounds<FloatDP> &);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_sin(const Bounds<FloatMP> &);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_cos(const Bounds<FloatMP> &);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_tan(const Bounds<FloatMP> &);
template OutputStream &Operations<Bounds<FloatMP>>::_write(OutputStream &, const Bounds<FloatMP> &);

template class Bounds<FloatMP>;
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_pi(MultiplePrecision);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_trunc(const Bounds<FloatMP> &);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_trunc(const Bounds<FloatMP> &, Nat);
template Integer Operations<Bounds<FloatMP>>::_cast_integer(const Bounds<FloatMP> &);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_add(const Bounds<FloatMP> &, const FloatMP &);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_add(const FloatMP &, const Bounds<FloatMP> &);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_sub(const Bounds<FloatMP> &, const FloatMP &);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_sub(const FloatMP &, const Bounds<FloatMP> &);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_mul(const Bounds<FloatMP> &, const FloatMP &);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_mul(const FloatMP &, const Bounds<FloatMP> &);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_div(const Bounds<FloatMP> &, const FloatMP &);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_div(const FloatMP &, const Bounds<FloatMP> &);
template InputStream &Operations<Bounds<FloatMP>>::_read(InputStream &, Bounds<FloatMP> &);

template<>
String class_name<Bounds<FloatDP>>() { return "FloatDPBounds"; }
template<>
String class_name<Bounds<FloatMP>>() { return "FloatMPBounds"; }

static_assert(TranscendentalField<Bounds<FloatDP>>);
static_assert(TranscendentalField<Bounds<FloatMP>>);

static_assert(OrderedLattice<Bounds<FloatDP>>);
static_assert(OrderedLattice<Bounds<FloatMP>>);

static_assert(WidenTranscendentalField<FloatDP>);
static_assert(WidenTranscendentalField<FloatMP>);

} // namespace Ariadne
