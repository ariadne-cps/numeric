/***************************************************************************
 *            numeric/float_literals.cpp
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

#include <cmath>
#include <limits>

#include "utility/macros.hpp"

#include "float_error.hpp"
#include "float_ball.hpp"
#include "float_lower_bound.hpp"
#include "float_upper_bound.hpp"
#include "float_approximation.hpp"

namespace Ariadne {

FloatError<DoublePrecision> operator""_error(long double lx) {
    double x=static_cast<double>(lx);
    if constexpr (std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits) {
        ARIADNE_PRECONDITION(static_cast<long double>(x)==lx);
    }
    return FloatError<DoublePrecision>(FloatDP(cast_exact(x),dp));
}


Float<DoublePrecision> operator""_exact(long double lx) {
    double x=static_cast<double>(lx);
    if constexpr (std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits) {
        ARIADNE_PRECONDITION(static_cast<long double>(x)==lx);
    }
    return Float<DoublePrecision>(ExactDouble(x),dp);
}

FloatBall<DoublePrecision> operator""_near(long double lx) {
    double x=static_cast<double>(lx);
    double e=0.0;
    if constexpr (std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits) {
        long double le=std::abs(static_cast<long double>(x)-lx);
        e=static_cast<double>(le);
        if(static_cast<long double>(e)<le) {
            e=std::nextafter(e,std::numeric_limits<double>::infinity());
        }
    }
    return FloatBall<DoublePrecision>(FloatDP(cast_exact(x),dp),FloatDP(cast_exact(e),dp));
}


FloatUpperBound<DoublePrecision> operator""_upper(long double lx) {
    double x=static_cast<double>(lx);
    if constexpr (std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits) {
        if(static_cast<long double>(x)<lx) {
            x=std::nextafter(x,std::numeric_limits<double>::infinity());
        }
    }
    return FloatUpperBound<DoublePrecision>(FloatDP(cast_exact(x),dp));
}


FloatLowerBound<DoublePrecision> operator""_lower(long double lx) {
    double x=static_cast<double>(lx);
    if constexpr (std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits) {
        if(static_cast<long double>(x)>lx) {
            x=std::nextafter(x,-std::numeric_limits<double>::infinity());
        }
    }
    return FloatLowerBound<DoublePrecision>(FloatDP(cast_exact(x),dp));
}


FloatApproximation<DoublePrecision> operator""_approx(long double lx) {
    double x=static_cast<double>(lx);
    return FloatApproximation<DoublePrecision>(FloatDP(cast_exact(x),dp));
}

} // namespace Ariadne
