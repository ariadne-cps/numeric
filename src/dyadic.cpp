/***************************************************************************
 *            dyadic.cpp
 *
 *  Copyright  2026
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Numeric.
 *
 *  Ariadne Numeric is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Numeric is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric. If not, see <https://www.gnu.org/licenses/>.
 */

#include "numeric/dyadic.hpp"

#include <limits>

#include "numeric/concepts.hpp"
#include "numeric/decimal.hpp"
#include "numeric/double.hpp"
#include "numeric/float_approximation.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_lower_bound.hpp"
#include "numeric/float_upper_bound.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/rational.hpp"
#include "numeric/twoexp.hpp"
#include "utility/string.hpp"
#include "utility/writable.hpp"

namespace Ariadne {

const Dyadic infty = Dyadic::inf();
Rational::Rational(const Dyadic& x) { *this = x.to_rational(); }
Dyadic::Dyadic(const FloatDP& x) : Dyadic(Rational(x)) {}
Dyadic::Dyadic(const FloatMP& x) : Dyadic(Rational(x)) {}

Dyadic operator+(TwoExp y) { return Dyadic(y); }
Dyadic operator-(TwoExp y) { return -Dyadic(y); }
Dyadic operator*(Dyadic x, TwoExp y) { return x * Dyadic(y); }
Dyadic operator/(Dyadic x, TwoExp y) { return x * Dyadic(rec(y)); }
Dyadic operator/(Int n, Two) { return Dyadic(n, 1u); }
Dyadic operator/(Integer n, Two) { return Dyadic(std::move(n), 1u); }
Dyadic operator/(ExactDouble d, Two) { return hlf(Dyadic(d)); }
Dyadic operator/(Dyadic x, Two) { return hlf(x); }
OutputStream& operator<<(OutputStream& os, TwoExp x) { return os << "2^" << x.exponent(); }
Writer<Dyadic> Dyadic::_default_writer(new DecimalWriter());
Void Dyadic::set_default_writer(Writer<Dyadic> writer) { _default_writer = writer; }
Writer<Dyadic> Dyadic::default_writer() { return _default_writer; }
OutputStream& operator<<(OutputStream& os, const Dyadic& x) {
    return os << Dyadic::default_writer()(x);
}
template<class X>
inline OutputStream& write_infinite(OutputStream& os, const X& x) {
    if (is_nan(x)) {
        os << "NaN";
    } else {
        os << (sgn(x) == Sign::POSITIVE ? "" : "-") << "inf";
    }
    return os;
}

auto DecimalWriter::_write(OutputStream& os, const Dyadic& x) const -> OutputStream& {
    if (is_finite(x)) {
        Dyadic w = x;
        if (w < 0) {
            os << "-";
            w = -w;
        }
        Integer z = floor(w);
        os << z << ".";
        w -= z;
        while (w != 0) {
            w *= 10;
            z = floor(w);
            w -= z;
            os << z;
        }
    } else {
        write_infinite(os, x);
    }
    return os;
}

auto ScientificWriter::_write(OutputStream& os, const Dyadic& x) const -> OutputStream& {
    if (is_finite(x)) {
        Dyadic w = x;
        if (w < 0) {
            os << "-";
            w = -w;
        }
        int e = 0;
        while (w >= 2) {
            w = hlf(w);
            e = e + 1;
        }
        while (w < 1) {
            w = w * 2;
            e = e - 1;
        }
        os << DecimalWriter()(w);
        os << "*2^" << e;
    } else {
        write_infinite(os, x);
    }
    return os;
}

auto FractionWriter::_write(OutputStream& os, const Dyadic& x) const -> OutputStream& {
    if (is_finite(x)) {
        const Rational q(x);
        os << q.numerator();
        Int exp = log2floor(q.denominator());
        if (exp != 0) {
            if (exp == 1) {
                os << "/2";
            } else {
                os << "/2^" << exp;
            }
        }
    } else {
        write_infinite(os, x);
    }
    return os;
}

auto RepresentationWriter<Dyadic>::_write(OutputStream& os, const Dyadic& x) const -> OutputStream& {
    return os << "Dyadic(" << x.mantissa() << "," << x.exponent() << "u)";
}

Dyadic make_dyadic(unsigned long long int n) { return Dyadic(Integer(n)); }
Dyadic hlf(const Integer& n) { return hlf(Dyadic(n)); }
// Exact dyadic endpoint values with directed rounding to each target.
Bounds<FloatMP> Bounds<Dyadic>::get(MultiplePrecision pr) const {
    return Bounds<FloatMP>(FloatMP(lower_raw(), FloatMP::ROUND_DOWNWARD, pr),
                           FloatMP(upper_raw(), FloatMP::ROUND_UPWARD, pr));
}
LowerBound<FloatMP> LowerBound<Dyadic>::get(MultiplePrecision pr) const {
    return LowerBound<FloatMP>(FloatMP(raw(), FloatMP::ROUND_DOWNWARD, pr));
}
UpperBound<FloatMP> UpperBound<Dyadic>::get(MultiplePrecision pr) const {
    return UpperBound<FloatMP>(FloatMP(raw(), FloatMP::ROUND_UPWARD, pr));
}
Approximation<FloatDP> Approximation<Dyadic>::get(DoublePrecision pr) const {
    return Approximation<FloatDP>(FloatDP(raw(), FloatDP::ROUND_TO_NEAREST, pr));
}
Approximation<FloatMP> Approximation<Dyadic>::get(MultiplePrecision pr) const {
    return Approximation<FloatMP>(FloatMP(raw(), FloatMP::ROUND_TO_NEAREST, pr));
}
template<>
String class_name<Dyadic>() { return "Dyadic"; }
template<>
String class_name<DyadicBounds>() { return "DyadicBounds"; }

static_assert(DyadicRing<Dyadic>);
static_assert(OrderedLattice<Dyadic>);
} // namespace Ariadne
