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
#include "floatmp.hpp"

namespace Ariadne {

int abslog10floor(double);

template<>
OutputStream& Operations<FloatBall<DoublePrecision>>::_write(OutputStream& os, const FloatBall<DoublePrecision>& x);

template<>
Ball<FloatDP, FloatDP>::Ball(const Real&, DoublePrecision);

template class Ball<FloatDP, FloatDP>;
template Ball<FloatDP, FloatDP> Operations<Ball<FloatDP, FloatDP>>::_trunc(const Ball<FloatDP, FloatDP>&);
template Ball<FloatDP, FloatDP> Operations<Ball<FloatDP, FloatDP>>::_trunc(const Ball<FloatDP, FloatDP>&, Nat);
template Integer Operations<Ball<FloatDP, FloatDP>>::_cast_integer(const Ball<FloatDP, FloatDP>&);
template InputStream& Operations<Ball<FloatDP, FloatDP>>::_read(InputStream&, Ball<FloatDP, FloatDP>&);
// Explicit instantiations for the FloatMP centre with double-precision error.
template OutputStream& Operations<Ball<FloatMP, FloatDP>>::_write(
    OutputStream&, const Ball<FloatMP, FloatDP>&);

// precisions. The generic template definitions use exact FloatMP values
// and round errors upward in the error precision.
template OutputStream& Operations<Ball<FloatMP, FloatMP>>::_write(
    OutputStream&, const Ball<FloatMP, FloatMP>&);

// NumberGetterMixin and the native interval algorithms require these
// constructors at link time. Use their shared, outward-rounded definitions.

// Export both native centre/error precision combinations with the same
template class Ball<FloatMP, FloatDP>;
template class Ball<FloatMP, FloatMP>;
template Ball<FloatMP, FloatDP> Operations<Ball<FloatMP, FloatDP>>::_trunc(const Ball<FloatMP, FloatDP>&);
template Ball<FloatMP, FloatDP> Operations<Ball<FloatMP, FloatDP>>::_trunc(const Ball<FloatMP, FloatDP>&, Nat);
template Integer Operations<Ball<FloatMP, FloatDP>>::_cast_integer(const Ball<FloatMP, FloatDP>&);
template InputStream& Operations<Ball<FloatMP, FloatDP>>::_read(InputStream&, Ball<FloatMP, FloatDP>&);
template Ball<FloatMP, FloatMP> Operations<Ball<FloatMP, FloatMP>>::_trunc(const Ball<FloatMP, FloatMP>&);
template Ball<FloatMP, FloatMP> Operations<Ball<FloatMP, FloatMP>>::_trunc(const Ball<FloatMP, FloatMP>&, Nat);
template Integer Operations<Ball<FloatMP, FloatMP>>::_cast_integer(const Ball<FloatMP, FloatMP>&);
template InputStream& Operations<Ball<FloatMP, FloatMP>>::_read(InputStream&, Ball<FloatMP, FloatMP>&);

template<>
String class_name<Ball<FloatDP>>() { return "FloatDPBall"; }
template<>
String class_name<Ball<FloatMP, FloatDP>>() { return "FloatMPDPBall"; }
template<>
String class_name<Ball<FloatMP, FloatMP>>() { return "FloatMPBall"; }

Ball<FloatDP> FloatDP::pm(const Error<FloatDP>& e) const { return Ball<FloatDP>(*this, e); }
Ball<FloatMP, FloatMP> FloatMP::pm(const Error<FloatMP>& e) const { return Ball<FloatMP, FloatMP>(*this, e); }
Ball<FloatMP, FloatDP> FloatMP::pm(const Error<FloatDP>& e) const { return Ball<FloatMP, FloatDP>(*this, e); }

template<>
OutputStream& Operations<FloatBall<DoublePrecision>>::_write(OutputStream& os, const FloatBall<DoublePrecision>& x) {
    const MultiplePrecision prec(std::size_t{64});
    return os << FloatBall<MultiplePrecision>(
               FloatMP(ExactDouble(x.value_raw().get_d()), prec),
               FloatMP(ExactDouble(x.error_raw().get_d()), prec));
}

} // namespace Ariadne
