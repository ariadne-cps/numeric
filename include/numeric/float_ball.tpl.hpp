/***************************************************************************
 *            numeric/float_ball.tpl.hpp
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

/** @file numeric/float_ball.tpl.hpp
 *  @brief Template implementations for float ball in Ariadne Numeric.
 */

#ifndef ARIADNE_FLOAT_BALL_TPL_HPP
#define ARIADNE_FLOAT_BALL_TPL_HPP

#include "decimal.hpp"
#include "dyadic.hpp"
#include "float_ball.hpp"
#include "float_bounds.hpp"
#include "float_error.hpp"
#include "foundation/logical.hpp"
#include "integer.hpp"
#include "number_wrapper.hpp"
#include "rational.hpp"

namespace Ariadne {

namespace {
template<class FE, class FLT>
    requires(not SameAs<FE, FLT>)
inline FE _make_error(FLT const& x) {
    typename FE::PrecisionType pre;
    return FE(Dyadic(x), upward, pre);
}
template<class FE, class FLT>
    requires SameAs<FE, FLT>
inline FE _make_error(FLT const& x) {
    return x;
}
template<class FE, class FLT, class PRE>
inline FE _make_error(FLT const& x, PRE pre) {
    return FE(Dyadic(x), upward, pre);
}
} // namespace

template<class F, class FE>
Ball<F, FE>::Ball(const F& v, const Error<FE>& e)
    : _v(v), _e(e.raw()) {}

template<class F, class FE>
Ball<F, FE>::Ball(const Bounds<F>& x)
    : _v(x.value_raw()), _e(_make_error<FE>(x.error_raw())) {}
template<class F, class FE>
Ball<F, FE>::Ball(const Bounds<F>& x, PRE pre)
    : _v(x.value_raw()), _e(_make_error<FE>(x.error_raw(), pre)) {}

template<class F, class FE>
Ball<F, FE>::Ball(const Ball<F, FE>& x, PR pr)
    : _v(x._v, near, pr), _e(x._e, up, _error_precision<PRE>(pr)) {
    if constexpr (not SameAs<PR, DoublePrecision>) {
        F d = (this->_v >= x._v) ? sub(up, this->_v, x._v) : sub(up, x._v, this->_v);
        _e = add(up, _e, _make_error<FE>(d));
    }
}

template<class F, class FE>
Ball<F, FE>::Ball(const Ball<F, FE>& x, PR pr, PRE pre)
    : _v(x._v, near, pr), _e(x._e, up, pre) {
    F d = (this->_v >= x._v) ? sub(up, this->_v, x._v) : sub(up, x._v, this->_v);
    _e = add(up, _e, _make_error<FE>(d, pre));
}

template<class F, class FE>
Ball<F, FE>::Ball(const ExactDouble& d, PR pr) : _v(d, pr), _e(0, _error_precision<PRE>(pr)) {}
template<class F, class FE>
Ball<F, FE>::Ball(const TwoExp& t, PR pr) : _v(t, pr), _e(0u, _error_precision<PRE>(pr)) {}
template<class F, class FE>
Ball<F, FE>::Ball(const Integer& z, PR pr) : _v(z, near, pr), _e(abs(Dyadic(_v) - z), up, _error_precision<PRE>(pr)) {}
template<class F, class FE>
Ball<F, FE>::Ball(const Dyadic& w, PR pr) : _v(w, near, pr), _e(abs(Dyadic(_v) - w), up, _error_precision<PRE>(pr)) {}
template<class F, class FE>
Ball<F, FE>::Ball(const Decimal& d, PR pr) : Ball(Rational(d), pr) {}
template<class F, class FE>
Ball<F, FE>::Ball(const Rational& q, PR pr) : _v(q, near, pr), _e(abs(Rational(_v) - q), up, _error_precision<PRE>(pr)) {}
template<class F, class FE>
Ball<F, FE>::Ball(const Integer& z, PR pr, PRE pre) : _v(F(z, near, pr)), _e(abs(Dyadic(_v) - z), up, pre) {}
template<class F, class FE>
Ball<F, FE>::Ball(const Dyadic& w, PR pr, PRE pre) : _v(F(w, near, pr)), _e(abs(Dyadic(_v) - w), up, pre) {}
template<class F, class FE>
Ball<F, FE>::Ball(const Rational& q, PR pr, PRE pre) : _v(F(q, near, pr)), _e(abs(Rational(_v) - q), up, pre) {}

template<class F, class FE>
Ball<F, FE>::Ball(const Decimal& dv, const Decimal& de, PrecisionType pr)
    : _v(dv, near, pr), _e(de, up, _error_precision<PRE>(pr)) {}

template<class F, class FE>
Ball<F, FE>::Ball(const Real& r, PR pr) : Ball(r.get(pr)) {}
template<class F, class FE>
Ball<F, FE>::Ball(const ValidatedNumber& y, PR pr) : Ball(y.get(pr, _error_precision<PRE>(pr))) {}
template<class F, class FE>
Ball<F, FE>::Ball(const Real& r, PR pr, PRE pre) : Ball(r.get(pr), pre) {}
template<class F, class FE>
Ball<F, FE>::Ball(const ValidatedNumber& y, PR pr, PRE pre) : Ball(y.get(pr, pre)) {}
template<class F, class FE>
Ball<F, FE>::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<Ball<F, FE>>(*this)); }

template<class F, class FE>
const LowerBound<F> Ball<F, FE>::lower() const { return LowerBound<F>(sub(down, this->_v, this->_e)); }
template<class F, class FE>
const UpperBound<F> Ball<F, FE>::upper() const { return UpperBound<F>(add(up, this->_v, this->_e)); }

template<class F, class FE>
const F Ball<F, FE>::value() const { return F(this->_v); }
template<class F, class FE>
const Error<FE> Ball<F, FE>::error() const { return Error<FE>(this->_e); }

template<class F, class FE>
Ball<F, FE> Ball<F, FE>::pm(const Error<FE>& e) const {
    return Ball<F, FE>(this->_v, add(up, this->_e, e.raw()));
}

template<class F, class FE>
Ball<F, FE> Operations<Ball<F, FE>>::_trunc(const Ball<F, FE>& x) {
    return Ball<F, FE>(trunc(Bounds<F>(x.lower_raw(), x.upper_raw())));
}

template<class F, class FE>
Ball<F, FE> Operations<Ball<F, FE>>::_trunc(const Ball<F, FE>& x, Nat n) {
    return Ball<F, FE>(trunc(Bounds<F>(x.lower_raw(), x.upper_raw()), n));
}

template<class F, class FE>
Integer Operations<Ball<F, FE>>::_cast_integer(const Ball<F, FE>& x) {
    Dyadic w = static_cast<Dyadic>(x.value_raw());
    Integer r = round(w);
    ARIADNE_ASSERT_MSG(abs(w - r) <= x.error_raw(), "cast_integer(Ball<" << class_name<F>() << "," << class_name<FE>() << "> const& x): "
                                                                                                                          "x="
                                                                         << x << "does not model an integer.")
    return r;
}

template<class F, class FE>
InputStream& Operations<Ball<F, FE>>::_read(InputStream& is, Ball<F, FE>& x) {
    static const char pmstr[] = "\xC2\xB1";
    char cpm[3] = {0, 0, 0};
    F _v(x.precision());
    FE _e(x.error_precision());
    auto rnd = F::get_rounding_mode();
    F::set_rounding_to_nearest();
    is >> _v;
    is >> cpm[0] >> cpm[1];
    F::set_rounding_mode(rnd);
    auto rnde = FE::get_rounding_mode();
    FE::set_rounding_upward();
    is >> _e;
    FE::set_rounding_mode(rnde);
    ARIADNE_ASSERT(not is.fail());
    ARIADNE_ASSERT(std::strcmp(cpm, pmstr) == 0);
    x._v = _v;
    x._e = _e;
    return is;
}

template<class F, class FE>
OutputStream& Operations<Ball<F, FE>>::_write(OutputStream& os, const Ball<F, FE>& x) {
    return os << x.value() << "\xC2\xB1" << x.error();
}

} // namespace Ariadne

#endif
