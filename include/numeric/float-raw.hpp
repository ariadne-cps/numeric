/***************************************************************************
 *            float-raw.hpp
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

/*! \file float-raw.hpp
 *  \brief Temporary header for raw floating-point type.
 */

#ifndef ARIADNE_FLOAT_RAW_HPP
#define ARIADNE_FLOAT_RAW_HPP

#include "float.decl.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"
#include "numeric/operators.hpp"

namespace Ariadne {

/*
inline FloatDP nul(exact,FloatDP x) { return nul(x); }
inline FloatDP pos(exact,FloatDP x) { return pos(x); }
inline FloatDP neg(exact,FloatDP x) { return neg(x); }
inline FloatDP hlf(exact,FloatDP x) { return hlf(x); }

inline FloatDP max(exact,FloatDP x1, FloatDP x2) { return max(x1,x2); }
inline FloatDP min(exact,FloatDP x1, FloatDP x2) { return min(x1,x2); }
inline FloatDP abs(exact,FloatDP x) { return abs(x); }


inline FloatDP add(approx,FloatDP x1, FloatDP x2) { return add(near,x1,x2); }
inline FloatDP add(near,FloatDP x1, FloatDP x2) { return add(near,x1,x2); }
inline FloatDP add(up,FloatDP x1, FloatDP x2) { return add(up,x1,x2); }
inline FloatDP add(down,FloatDP x1, FloatDP x2) { return add(down,x1,x2); }

inline FloatDP sub(approx,FloatDP x1, FloatDP x2) { return sub(near,x1,x2); }
inline FloatDP sub(near,FloatDP x1, FloatDP x2) { return sub(near,x1,x2); }
inline FloatDP sub(up,FloatDP x1, FloatDP x2) { return sub(up,x1,x2); }
inline FloatDP sub(down,FloatDP x1, FloatDP x2) { return sub(down,x1,x2); }

inline FloatDP mul(approx,FloatDP x1, FloatDP x2) { return mul(near,x1,x2); }
inline FloatDP mul(near,FloatDP x1, FloatDP x2) { return mul(near,x1,x2); }
inline FloatDP mul(up,FloatDP x1, FloatDP x2) { return mul(up,x1,x2); }
inline FloatDP mul(down,FloatDP x1, FloatDP x2) { return mul(down,x1,x2); }

inline FloatDP div(approx,FloatDP x1, FloatDP x2) { return div(near,x1,x2); }
inline FloatDP div(near,FloatDP x1, FloatDP x2) { return div(near,x1,x2); }
inline FloatDP div(up,FloatDP x1, FloatDP x2) { return div(up,x1,x2); }
inline FloatDP div(down,FloatDP x1, FloatDP x2) { return div(down,x1,x2); }

inline FloatDP pow(approx,FloatDP x, Int n) { return pow(near,x,n); }
inline FloatDP pow(up,FloatDP x, Int n) { return pow(up,x,n); }
inline FloatDP pow(down,FloatDP x, Int n) { return pow(down,x,n); }

inline FloatDP next(down,FloatDP x) { return sub(down,x,FloatDP::min(dp)); }
inline FloatDP next(up,FloatDP x) { return add(up,x,FloatDP::min(dp)); }

inline FloatDP rec(approx,FloatDP x) { return rec(near,x); }
inline FloatDP rec(near,FloatDP x) { return rec(near,x); }
inline FloatDP rec(up,FloatDP x) { return rec(up,x); }
inline FloatDP rec(down,FloatDP x) { return rec(down,x); }

inline FloatDP sqrt(up,FloatDP const& x) { return sqrt(up,x); }
inline FloatDP sqrt(down,FloatDP const& x) { return sqrt(down,x); }

inline FloatDP exp(up,FloatDP const& x) { return exp(up,x); }
inline FloatDP exp(down,FloatDP const& x) { return exp(down,x); }

inline FloatDP log(up,FloatDP const& x) { return log(up,x); }
inline FloatDP log(down,FloatDP const& x) { return log(down,x); }

inline FloatDP sin(up,FloatDP const& x) { return sin(up,x); }
inline FloatDP sin(down,FloatDP const& x) { return sin(down,x); }

inline FloatDP cos(up,FloatDP const& x) { return cos(up,x); }
inline FloatDP cos(down,FloatDP const& x) { return cos(down,x); }

inline FloatDP atan(up,FloatDP const& x) { return atan(up,x); }
inline FloatDP atan(down,FloatDP const& x) { return atan(down,x); }

inline FloatDP pi_up(DoublePrecision const& pr) { return FloatDP::pi_up(pr); }
inline FloatDP pi_approx(DoublePrecision const& pr) { return FloatDP::pi_near(pr); }
inline FloatDP pi_down(DoublePrecision const& pr) { return FloatDP::pi_down(pr); }




















*/
/*
inline FloatDP max(exact,FloatDP const& x1, FloatDP const& x2) { return max(x1,x2); }
inline FloatDP min(exact,FloatDP const& x1, FloatDP const& x2) { return min(x1,x2); }
inline FloatDP abs(exact,FloatDP const& x) { return abs(x); }

inline FloatDP add(FloatDP const& x1, FloatDP const& x2, FloatDP::RoundingModeType rnd) { return apply(Add(),x1,x2,rnd); }
inline FloatDP sub(FloatDP const& x1, FloatDP const& x2, FloatDP::RoundingModeType rnd) { return apply(Sub(),x1,x2,rnd); }
inline FloatDP mul(FloatDP const& x1, FloatDP const& x2, FloatDP::RoundingModeType rnd) { return apply(Mul(),x1,x2,rnd); }
inline FloatDP div(FloatDP const& x1, FloatDP const& x2, FloatDP::RoundingModeType rnd) { return apply(Div(),x1,x2,rnd); }
inline FloatDP pow(FloatDP const& x, Int n, FloatDP::RoundingModeType rnd) { return apply(Pow(),x,n,rnd); }

inline FloatDP nul(exact,FloatDP const& x) { return nul(x); }
inline FloatDP pos(exact,FloatDP const& x) { return pos(x); }
inline FloatDP neg(exact,FloatDP const& x) { return neg(x); }
inline FloatDP hlf(exact,FloatDP const& x) { return hlf(x); }

inline FloatDP add(up,FloatDP const& x1, FloatDP const& x2) { return add(x1,x2,FloatDP::upward); }
inline FloatDP add(down,FloatDP const& x1, FloatDP const& x2) { return add(x1,x2,FloatDP::downward); }
inline FloatDP add(near,FloatDP const& x1, FloatDP const& x2) { return add(x1,x2,FloatDP::to_nearest); }
inline FloatDP add(approx,FloatDP const& x1, FloatDP const& x2) { return add(x1,x2,FloatDP::to_nearest); }

inline FloatDP sub(up,FloatDP const& x1, FloatDP const& x2) { return sub(x1,x2,FloatDP::upward); }
inline FloatDP sub(down,FloatDP const& x1, FloatDP const& x2) { return sub(x1,x2,FloatDP::downward); }
inline FloatDP sub(near,FloatDP const& x1, FloatDP const& x2) { return sub(x1,x2,FloatDP::to_nearest); }
inline FloatDP sub(approx,FloatDP const& x1, FloatDP const& x2) { return sub(x1,x2,FloatDP::to_nearest); }

inline FloatDP mul(up,FloatDP const& x1, FloatDP const& x2) { return mul(x1,x2,FloatDP::upward); }
inline FloatDP mul(down,FloatDP const& x1, FloatDP const& x2) { return mul(x1,x2,FloatDP::downward); }
inline FloatDP mul(near,FloatDP const& x1, FloatDP const& x2) { return mul(x1,x2,FloatDP::to_nearest); }
inline FloatDP mul(approx,FloatDP const& x1, FloatDP const& x2) { return mul(x1,x2,FloatDP::to_nearest); }

inline FloatDP div(up,FloatDP const& x1, FloatDP const& x2) { return div(x1,x2,FloatDP::upward); }
inline FloatDP div(down,FloatDP const& x1, FloatDP const& x2) { return div(x1,x2,FloatDP::downward); }
inline FloatDP div(near,FloatDP const& x1, FloatDP const& x2) { return div(x1,x2,FloatDP::to_nearest); }
inline FloatDP div(approx,FloatDP const& x1, FloatDP const& x2) { return div(x1,x2,FloatDP::to_nearest); }

inline FloatDP fma(up,FloatDP const& x1, FloatDP const& x2, FloatDP const& y) { return add(mul(x1,x2,FloatDP::upward),y,upward); }
inline FloatDP fma(down,FloatDP const& x1, FloatDP const& x2, FloatDP const& y) { return add(mul(x1,x2,FloatDP::downward),y,downward); }
inline FloatDP fma(approx,FloatDP const& x1, FloatDP const& x2, FloatDP const& y) { return add(mul(x1,x2,FloatDP::to_nearest),y,to_nearest); }

inline FloatDP next(up,FloatDP const& x) { return add(up,x,FloatDP::min()); }
inline FloatDP next(down,FloatDP const& x) { return sub(down,x,FloatDP::min()); }

inline FloatDP rad(up,FloatDP const& x1, FloatDP const& x2) { return hlf(sub(up,x2,x1)); }
inline FloatDP med(near,FloatDP const& x1, FloatDP const& x2) { return hlf(add(near,x2,x1)); }

inline FloatDP pow(up,FloatDP const& x, Int n) { return pow(x,n,FloatDP::upward); }
inline FloatDP pow(down,FloatDP const& x, Int n) { return pow(x,n,FloatDP::downward); }
inline FloatDP pow(approx,FloatDP const& x, Int n) { return pow(x,n,FloatDP::to_nearest); }

FloatDP pi_up();
FloatDP pi_down();
FloatDP pi_approx();

inline FloatDP sqrt(approx,FloatDP x) { return std::sqrt(x.dbl); }
inline FloatDP exp(approx,FloatDP x) { return std::exp(x.dbl); }
inline FloatDP log(approx,FloatDP x) { return std::log(x.dbl); }
inline FloatDP sin(approx,FloatDP x) { return std::sin(x.dbl); }
inline FloatDP cos(approx,FloatDP x) { return std::cos(x.dbl); }
inline FloatDP tan(approx,FloatDP x) { return std::tan(x.dbl); }
inline FloatDP asin(approx,FloatDP x) { return std::asin(x.dbl); }
inline FloatDP acos(approx,FloatDP x) { return std::acos(x.dbl); }
inline FloatDP atan(approx,FloatDP x) { return std::atan(x.dbl); }

*/

} // namespace Ariadne

#endif /* ARIADNE_FLOAT_RAW_HPP */
