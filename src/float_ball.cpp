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
#include "floatap.hpp"

namespace Ariadne {

int abslog10floor(double);

template<> OutputStream& Operations<FloatBall<DoublePrecision>>::_write(OutputStream& os, FloatBall<DoublePrecision> const& x);

template<> Ball<FloatDP,FloatDP>::Ball(Real const&, DoublePrecision);

template class Ball<FloatDP,FloatDP>;
template Ball<FloatDP,FloatDP> Operations<Ball<FloatDP,FloatDP>>::_trunc(Ball<FloatDP,FloatDP> const&);
template Ball<FloatDP,FloatDP> Operations<Ball<FloatDP,FloatDP>>::_trunc(Ball<FloatDP,FloatDP> const&, Nat);
template Integer Operations<Ball<FloatDP,FloatDP>>::_cast_integer(Ball<FloatDP,FloatDP> const&);
template InputStream& Operations<Ball<FloatDP,FloatDP>>::_read(InputStream&, Ball<FloatDP,FloatDP>&);
// Explicit instantiations for the FloatAP centre with double-precision error.
template OutputStream& Operations<Ball<FloatAP,FloatDP>>::_write(
    OutputStream&, Ball<FloatAP,FloatDP> const&);


// precisions. The generic template definitions use exact FloatAP values
// and round errors upward in the error precision.
template OutputStream& Operations<Ball<FloatAP,FloatAP>>::_write(
    OutputStream&, Ball<FloatAP,FloatAP> const&);

// NumberGetterMixin and the native interval algorithms require these
// constructors at link time. Use their shared, outward-rounded definitions.


// Export both native centre/error precision combinations with the same
template class Ball<FloatAP,FloatDP>;
template class Ball<FloatAP,FloatAP>;
template Ball<FloatAP,FloatDP> Operations<Ball<FloatAP,FloatDP>>::_trunc(Ball<FloatAP,FloatDP> const&);
template Ball<FloatAP,FloatDP> Operations<Ball<FloatAP,FloatDP>>::_trunc(Ball<FloatAP,FloatDP> const&, Nat);
template Integer Operations<Ball<FloatAP,FloatDP>>::_cast_integer(Ball<FloatAP,FloatDP> const&);
template InputStream& Operations<Ball<FloatAP,FloatDP>>::_read(InputStream&, Ball<FloatAP,FloatDP>&);
template Ball<FloatAP,FloatAP> Operations<Ball<FloatAP,FloatAP>>::_trunc(Ball<FloatAP,FloatAP> const&);
template Ball<FloatAP,FloatAP> Operations<Ball<FloatAP,FloatAP>>::_trunc(Ball<FloatAP,FloatAP> const&, Nat);
template Integer Operations<Ball<FloatAP,FloatAP>>::_cast_integer(Ball<FloatAP,FloatAP> const&);
template InputStream& Operations<Ball<FloatAP,FloatAP>>::_read(InputStream&, Ball<FloatAP,FloatAP>&);


template<> String class_name<Ball<FloatDP>>() { return "FloatDPBall"; }
template<> String class_name<Ball<FloatAP,FloatDP>>() { return "FloatAPDPBall"; }
template<> String class_name<Ball<FloatAP,FloatAP>>() { return "FloatAPBall"; }

Ball<FloatDP> FloatDP::pm(Error<FloatDP> const& e) const { return Ball<FloatDP>(*this,e); }
Ball<FloatAP,FloatAP> FloatAP::pm(Error<FloatAP> const& e) const { return Ball<FloatAP,FloatAP>(*this,e); }
Ball<FloatAP,FloatDP> FloatAP::pm(Error<FloatDP> const& e) const { return Ball<FloatAP,FloatDP>(*this,e); }


template<> OutputStream& Operations<FloatBall<DoublePrecision>>::_write(OutputStream& os, FloatBall<DoublePrecision> const& x) {
    const ArbitraryPrecision prec(std::size_t{64});
    return os << FloatBall<ArbitraryPrecision>(
        FloatAP(ExactDouble(x.value_raw().get_d()),prec),
        FloatAP(ExactDouble(x.error_raw().get_d()),prec));
}

} // namespace Ariadne
