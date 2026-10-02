/***************************************************************************
 *            numeric/float_ball.cpp
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

#include "float_ball.hpp"
#include "float_ball.tpl.hpp"

namespace Ariadne {

int abslog10floor(double);

template<> OutputStream& Operations<FloatBall<DoublePrecision>>::_write(OutputStream& os, FloatBall<DoublePrecision> const& x);
template<> OutputStream& Operations<FloatBall<MultiplePrecision,DoublePrecision>>::_write(OutputStream& os, FloatBall<MultiplePrecision,DoublePrecision> const& x);
template<> OutputStream& Operations<FloatBall<MultiplePrecision>>::_write(OutputStream& os, FloatBall<MultiplePrecision> const& x);

template<> Ball<FloatDP,FloatDP>::Ball(Real const&, DoublePrecision);

template class Ball<FloatDP,FloatDP>;
template Ball<FloatDP,FloatDP> Operations<Ball<FloatDP,FloatDP>>::_trunc(Ball<FloatDP,FloatDP> const&);
template Ball<FloatDP,FloatDP> Operations<Ball<FloatDP,FloatDP>>::_trunc(Ball<FloatDP,FloatDP> const&, Nat);
template Integer Operations<Ball<FloatDP,FloatDP>>::_cast_integer(Ball<FloatDP,FloatDP> const&);
template InputStream& Operations<Ball<FloatDP,FloatDP>>::_read(InputStream&, Ball<FloatDP,FloatDP>&);
template class Ball<FloatMP,FloatDP>;
template Ball<FloatMP,FloatDP> Operations<Ball<FloatMP,FloatDP>>::_trunc(Ball<FloatMP,FloatDP> const&);
template Ball<FloatMP,FloatDP> Operations<Ball<FloatMP,FloatDP>>::_trunc(Ball<FloatMP,FloatDP> const&, Nat);
template Integer Operations<Ball<FloatMP,FloatDP>>::_cast_integer(Ball<FloatMP,FloatDP> const&);
template InputStream& Operations<Ball<FloatMP,FloatDP>>::_read(InputStream&, Ball<FloatMP,FloatDP>&);
template class Ball<FloatMP,FloatMP>;
template Ball<FloatMP,FloatMP> Operations<Ball<FloatMP,FloatMP>>::_trunc(Ball<FloatMP,FloatMP> const&);
template Ball<FloatMP,FloatMP> Operations<Ball<FloatMP,FloatMP>>::_trunc(Ball<FloatMP,FloatMP> const&, Nat);
template Integer Operations<Ball<FloatMP,FloatMP>>::_cast_integer(Ball<FloatMP,FloatMP> const&);
template InputStream& Operations<Ball<FloatMP,FloatMP>>::_read(InputStream&, Ball<FloatMP,FloatMP>&);

template<> String class_name<Ball<FloatDP>>() { return "FloatDPBall"; }
template<> String class_name<Ball<FloatMP,FloatDP>>() { return "FloatMDPBall"; }
template<> String class_name<Ball<FloatMP>>() { return "FloatMPBall"; }

Ball<FloatDP> FloatDP::pm(Error<FloatDP> const& e) const { return Ball<FloatDP>(*this,e); }
Ball<FloatMP> FloatMP::pm(Error<FloatMP> const& e) const { return Ball<FloatMP>(*this,e); }
Ball<FloatMP,FloatDP> FloatMP::pm(Error<FloatDP> const& e) const { return Ball<FloatMP,FloatDP>(*this,e); }


template<> OutputStream& Operations<FloatBall<MultiplePrecision>>::_write(OutputStream& os, FloatBall<MultiplePrecision> const& x) {
    // Write based on number of correct digits
    static const double log2ten = 3.3219280948873621817;
    static const char pmstr[] = "\xC2\xB1";
    static const char hlfstr[] = "\xC2\xBD";
    FloatMP const& v=x.value_raw();
    FloatMP const& e=x.error_raw();
    double edbl=e.get_d();
    // Compute the number of decimal places to be displayed
    Nat errplc = static_cast<Nat>(FloatError<MultiplePrecision>::output_places);
    Nat log10err = static_cast<Nat>(abslog10floor(edbl));
    Nat dgtserr = errplc-(log10err+1);
    Nat dgtsval = (x.value()==0) ? dgtserr : static_cast<Nat>(std::floor((x.value().precision()+1-x.value().exponent())/log2ten));
    Nat dgts = std::max(std::min(dgtsval,dgtserr),errplc);
    if(edbl==0.0) { dgts = dgtsval; }

    DecimalPlaces plcs{dgts};

    // Get string version of mpfr values
    String vstr=print(v,plcs,MPFR_RNDN);
    String estr=print(e,plcs,MPFR_RNDU);

    // Find position of first significan digit of error
    auto vcstr=vstr.c_str(); auto ecstr=estr.c_str();
    size_t cpl=0;
    if(edbl==0.0) {
        cpl=std::strlen(vcstr);
    } else if(edbl<1.0) {
        const char* vptr = std::strchr(vcstr,'.');
        const char* eptr = std::strchr(ecstr,'.');
        ++vptr; ++eptr;
        while((*eptr)=='0') { ++eptr; ++vptr; }

        cpl = static_cast<size_t>(vptr-vcstr);
    }


    // Chop and concatenate strings
    String ostr=vstr.substr(0,cpl);
    ostr += "[";
    ostr += vstr.substr(cpl);
    ostr += pmstr;
    const size_t ecpl=std::min(cpl,estr.size());
    ostr += estr.substr(ecpl);
    ostr += hlfstr;
    ostr += "]";
    return os << ostr;
}


template<> OutputStream& Operations<FloatBall<MultiplePrecision,DoublePrecision>>::_write(OutputStream& os, FloatBall<MultiplePrecision,DoublePrecision> const& x) {
    MultiplePrecision prec(64);
    return os << FloatBall<MultiplePrecision,MultiplePrecision>(x.value_raw(),FloatMP(x.error_raw(),prec));
}


template<> OutputStream& Operations<FloatBall<DoublePrecision>>::_write(OutputStream& os, FloatBall<DoublePrecision> const& x) {
    MultiplePrecision prec(64);
    return os << FloatBall<MultiplePrecision>(FloatMP(x.value_raw(),prec),FloatMP(x.error_raw(),prec));
}

} // namespace Ariadne
