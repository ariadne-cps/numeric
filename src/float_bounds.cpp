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

template<> auto Operations<FloatBounds<DoublePrecision>>::_write(OutputStream& os, const FloatBounds<DoublePrecision>& x) -> OutputStream&;
template<> auto Operations<FloatBounds<MultiplePrecision>>::_write(OutputStream& os, const FloatBounds<MultiplePrecision>& x) -> OutputStream&;

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
template class Bounds<FloatMP>;
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_pi(MultiplePrecision);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_sin(Bounds<FloatMP> const&);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_cos(Bounds<FloatMP> const&);
template Bounds<FloatMP> Operations<Bounds<FloatMP>>::_tan(Bounds<FloatMP> const&);
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

Int abslog10floor(FloatMP const&);

namespace {

OutputStream& write_bounds_with_error_places(OutputStream& os, const FloatMP& l, const FloatMP& u, Nat errplc)
{
    static const double log2ten = 3.3219280948873621817;
    using std::max; using std::min;
    if(l==0.0_x && u==0.0_x) { return os << "0.0[:]"; }

    int precplc=static_cast<int>(min(l.precision(),u.precision())/log2ten);
    FloatMP wdth=sub(up,u,l);
    FloatMP mag=max(neg(l),u);
    int log10wdth=max(abslog10floor(wdth),std::numeric_limits<int>::min()+static_cast<int>(errplc));
    int log10mag=abslog10floor(mag);
    int dgtswdth=static_cast<int>(errplc)-(log10wdth+1); // Digits appropriate given width of interval
    //int dgtsbnd=bndplc-(log10mag+1); // Digits appropriate given asked-for precision of bounded objects
    int dgtsprec=precplc-(log10mag+1); // Digits appropriate given precision of objects
    uint dgts=static_cast<uint>(max(min(dgtswdth,dgtsprec),1));
    DecimalPlaces plcs{dgts};

    String lstr=print(l,plcs,MPFR_RNDD);
    String ustr=print(u,plcs,MPFR_RNDU);
    auto lcstr=lstr.c_str();
    size_t cpl=0;
    if((l>=0)==(u>=0) && abslog10floor(l)==abslog10floor(u)) {
        while(lcstr[cpl]!='\0' && lcstr[cpl]==ustr[cpl]) { ++cpl; }

    }

    String ostr=lstr.substr(0,cpl);
    ostr += "[";
    ostr += lstr.substr(cpl);
    ostr += ":";
    ostr += ustr.substr(cpl);
    ostr += "]";
    return os << ostr;
}

} //namespace

template<> OutputStream& Operations<FloatBounds<MultiplePrecision>>::_write(OutputStream& os, const FloatBounds<MultiplePrecision>& x)
{
    return write_bounds_with_error_places(os,x.lower_raw(),x.upper_raw(),FloatError<MultiplePrecision>::output_places);
}

template<> OutputStream& Operations<FloatBounds<DoublePrecision>>::_write(OutputStream& os, const FloatBounds<DoublePrecision>& x)
{
    MultiplePrecision prec(53);
    return write_bounds_with_error_places(os,FloatMP(x.lower_raw(),prec),FloatMP(x.upper_raw(),prec),FloatError<DoublePrecision>::output_places);
}

static_assert(TranscendentalField<Bounds<FloatDP>>);
static_assert(TranscendentalField<Bounds<FloatMP>>);

static_assert(OrderedLattice<Bounds<FloatDP>>);
static_assert(OrderedLattice<Bounds<FloatMP>>);

static_assert(WidenTranscendentalField<FloatDP>);
static_assert(WidenTranscendentalField<FloatMP>);

} // namespace Ariadne
