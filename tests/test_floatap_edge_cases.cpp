/***************************************************************************
 *            test_floatap_edge_cases.cpp
 *
 *  Copyright  2026  Luca Geretti
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

#include "numeric/floatap.hpp"
#include "numeric/integer.hpp"
#include "numeric/rational.hpp"
#include <cstdint>
#include <gmp.h>
#include <iostream>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
using K = Ariadne::FloatAP;
using R = K::Rounding;
namespace {
struct Z {
    mpz_t z;
    explicit Z(long n) {
        mpz_init_set_si(z, n);
    }
    ~Z() {
        mpz_clear(z);
    }
};
struct Q {
    mpq_t q;
    Q() {
        mpq_init(q);
    }
    ~Q() {
        mpq_clear(q);
    }
};
void check(bool x, char const* msg) {
    if (!x) throw std::runtime_error(msg);
}
template< class E, class F > void raises(F f, char const* m) {
    bool ok = false;
    try {
        f();
    }
    catch (E const& ) {
        ok = true;
    }
    check(ok, m);
}
K v(long a, long b = 1, unsigned long p = 16, R r = R::nearest) {
    Z n(a), d(b);
    return K::from_ratio(n.z, d.z, 0, p, r);
}
void same(K const& a, K const& b, char const* m) {
    check(a.kind() == b.kind(), m);
    if (a.is_nan()) return;
    check(a.negative() == b.negative(), m);
    if (!a.is_finite()) return;
    Q x, y;
    a.to_rational(x.q);
    b.to_rational(y.q);
    check(mpq_cmp(x.q, y.q) == 0, m);
}
void exponent_bit_length_boundaries() {
    using Ariadne::floatap_exponent_from_bit_length;
    const auto maximum = std::numeric_limits<std::int64_t>::max();

    check(floatap_exponent_from_bit_length(0u) == 0, "zero bit length");
    check(floatap_exponent_from_bit_length(64u) == 64, "ordinary bit length");
    check(floatap_exponent_from_bit_length(static_cast<std::size_t>(maximum)) == maximum,
          "maximum representable bit length");

    if constexpr (std::numeric_limits<std::size_t>::digits >
                  std::numeric_limits<std::int64_t>::digits) {
        raises<std::overflow_error>([&] {
            (void)floatap_exponent_from_bit_length(static_cast<std::size_t>(maximum) + 1u);
        }, "bit length above signed exponent range");
    }
}

void precision_tag() {
    using namespace Ariadne;
    const ArbitraryPrecision via_bits(53_bits);
    const ArbitraryPrecision via_count(std::size_t { 53 });
    const ArbitraryPrecision distinct(113_bits);
    check(via_bits.bits() == 53u, "ArbitraryPrecision Bits constructor");
    check(via_count.bits() == 53u, "ArbitraryPrecision size_t constructor");
    check(distinct.bits() == 113u, "ArbitraryPrecision arbitrary bit count");
    check(via_bits == via_count, "ArbitraryPrecision equality");
    check(!(via_bits == distinct), "ArbitraryPrecision inequality");
}
void int64_boundaries() {
    using namespace Ariadne;
    const ArbitraryPrecision precision(64_bits);
    const K zero(std::int64_t { 0 }, precision);
    const K positive(std::numeric_limits<std::int64_t>::max(), precision);
    const K negative(std::numeric_limits<std::int64_t>::min(), precision);
    const Rational z = static_cast<Rational>(zero);
    const Rational p = static_cast<Rational>(positive);
    const Rational n = static_cast<Rational>(negative);
    check(mpq_cmp_si(z.get_mpq(), 0, 1) == 0, "int64 zero");
    // mpq_cmp_si takes a long, which is 32 bits on 64-bit Windows.
    // Compare against independently parsed, exact 64-bit values instead.
    const Rational expected_positive(String("9223372036854775807"));
    const Rational expected_negative(String("-9223372036854775808"));
    check(mpq_cmp(p.get_mpq(), expected_positive.get_mpq()) == 0,
          "int64 maximum");
    check(mpq_cmp(n.get_mpq(), expected_negative.get_mpq()) == 0,
          "int64 minimum");
    check(is_zero(z) && !positive.negative() && negative.negative(),
          "int64 signs");
}
void rational_integration() {
    using namespace Ariadne;
    const ArbitraryPrecision precision(16_bits);
    const Rational finite(Integer(3), Integer(7));
    const K nearest(finite, R::nearest, precision);
    const K downward(finite, R::downward, precision);
    same(nearest, v(3, 7, 16, R::nearest), "Rational nearest constructor");
    same(downward, v(3, 7, 16, R::downward), "Rational downward constructor");
    same(K(finite, precision), nearest, "Rational default rounding constructor");
    check(nearest.get_d() == static_cast<Rational>(nearest).get_d(),
    "finite FloatAP double conversion");

    same(K(Integer(-19), R::nearest, precision), v(-19),
    "Integer constructor");
    same(K(std::int64_t { - 42 }, precision), v(-42),
    "int64 constructor");

    const K positive_inf(Rational::inf(Sign::POSITIVE), R::nearest, precision);
    const K negative_inf(Rational::inf(Sign::NEGATIVE), R::nearest, precision);
    const K not_a_number(Rational::nan(), R::nearest, precision);
    check(positive_inf.is_inf() && !positive_inf.negative(),
    "positive Rational infinity constructor");
    check(negative_inf.is_inf() && negative_inf.negative(),
    "negative Rational infinity constructor");
    check(not_a_number.is_nan(), "Rational NaN constructor");

    check(mpz_sgn(mpq_denref(static_cast<Rational>(positive_inf).get_mpq())) == 0,
    "positive infinity to Rational");
    check(mpz_sgn(mpq_denref(static_cast<Rational>(not_a_number).get_mpq())) == 0,
    "NaN to Rational");
}

void public_comparisons() {
    using namespace Ariadne;
    using Ariadne::Comparison;
    using Ariadne::Sign;
    const ArbitraryPrecision precision(16_bits);
    const K low = v(-3);
    const K mid = v(0);
    const K high = v(5);
    const K positive_inf = K::inf(Sign::POSITIVE, precision);
    const K negative_inf = K::inf(Sign::NEGATIVE, precision);
    const K not_a_number = K::nan(precision);

    // All three finite comparison outcomes, and both unordered orders.
    check(cmp(low, mid) == Comparison::LESS, "finite comparison less");
    check(cmp(mid, mid) == Comparison::EQUAL, "finite comparison equal");
    check(cmp(high, mid) == Comparison::GREATER, "finite comparison greater");
    check(cmp(not_a_number, mid) == Comparison::INCOMPARABLE, "NaN left comparison");
    check(cmp(mid, not_a_number) == Comparison::INCOMPARABLE, "NaN right comparison");

    // Both infinities, equal infinities, finite/infinite comparisons.
    check(cmp(positive_inf, positive_inf) == Comparison::EQUAL, "+inf equals +inf");
    check(cmp(negative_inf, negative_inf) == Comparison::EQUAL, "-inf equals -inf");
    check(cmp(negative_inf, positive_inf) == Comparison::LESS, "-inf below +inf");
    check(cmp(positive_inf, negative_inf) == Comparison::GREATER, "+inf above -inf");
    check(cmp(positive_inf, high) == Comparison::GREATER, "+inf above finite");
    check(cmp(negative_inf, high) == Comparison::LESS, "-inf below finite");
    check(cmp(high, positive_inf) == Comparison::LESS, "finite below +inf");
    check(cmp(high, negative_inf) == Comparison::GREATER, "finite above -inf");

    // Operators must keep NaN unordered, and agree on all ordered cases.
    check(mid == mid && !(mid != mid), "equal operators");
    check(low != high && !(low == high), "unequal operators");
    check(low < mid && !(high < mid), "less operators");
    check(high > mid && !(low > mid), "greater operators");
    check(low <= mid && mid <= mid && !(high <= mid), "less-or-equal operators");
    check(high >= mid && mid >= mid && !(low >= mid), "greater-or-equal operators");
    check(!(not_a_number == mid) && not_a_number != mid, "NaN equality operators");
    check(!(not_a_number < mid) && !(not_a_number > mid), "NaN strict ordering");
    check(!(not_a_number <= mid) && !(not_a_number >= mid), "NaN inclusive ordering");
    check(!(mid <= not_a_number) && !(mid >= not_a_number), "NaN right operand");

    // ADL-visible free classifiers and precision-tag constructors.
    check(is_nan(not_a_number) && !is_nan(mid), "public NaN classification");
    check(is_inf(positive_inf) && !is_inf(mid), "public infinity classification");
    check(is_finite(mid) && !is_finite(positive_inf), "public finite classification");
    check(is_zero(mid) && !is_zero(low), "public zero classification");
    check(positive_inf.precision().bits() == 16u && not_a_number.precision().bits() == 16u,
    "precision-tag special constructors");
    check(static_cast<Rational>(negative_inf).get_d() < 0.0,
    "negative infinity rational conversion");
}

void public_rounding() {
    const K seven_bits = v(1, 3, 7, R::nearest);
    const K thirteen_bits = v(2, 3, 13, R::nearest);
    const K eleven_bits = v(1, 5, 11, R::nearest);
    const R mode = R::downward;

    // Same kernels, different precision order: the result must use the
    // lower requested precision without dynamic backend dispatch.
    const auto check_pair =[& ](K const& a, K const& b, unsigned long bits) {
        const auto verify =[& ](K const& actual, K const& expected, char const* name) {
            same(actual, expected, name);
            check(actual.precision().bits() == bits, "public operation result precision");
        };
        verify(add(mode, a, b), K::add(a, b, bits, mode), "public add");
        verify(sub(mode, a, b), K::sub(a, b, bits, mode), "public sub");
        verify(mul(mode, a, b), K::mul(a, b, bits, mode), "public mul");
        verify(div(mode, a, b), K::div(a, b, bits, mode), "public div");
    };
    check_pair(seven_bits, thirteen_bits, 7);
    check_pair(thirteen_bits, seven_bits, 7);

    // The ternary precision is the minimum of all three input precisions.
    const auto verify_fma =[& ](K const& a, K const& b, K const& c, unsigned long bits) {
        const K actual = fma(mode, a, b, c);
        same(actual, K::fma(a, b, c, bits, mode), "public fused multiply-add");
        check(actual.precision().bits() == bits, "public FMA result precision");
    };
    verify_fma(seven_bits, thirteen_bits, eleven_bits, 7);
    verify_fma(thirteen_bits, eleven_bits, seven_bits, 7);
    verify_fma(thirteen_bits, seven_bits, eleven_bits, 7);
    verify_fma(eleven_bits, thirteen_bits, thirteen_bits, 11);
}

void conversions() {
    raises < std::invalid_argument >([] {
        K x(1);
    }, "precision");
    raises < std::domain_error >([] {
        v(1, 0);
    }, "denominator");
    same(v(2,-3), v(-2, 3), "signed denominator");
    same(v(-2,-3), v(2, 3), "both negative");
    Z i(-129);
    same(K::from_integer(i.z, 7, R::downward), v(-129, 1, 7, R::downward), "integer");
    same(v(5, 4, 2), v(1, 1, 2), "nearest even low");
    same(v(7, 4, 2), v(2, 1, 2), "nearest even high");
    same(v(5, 4, 2, R::upward), v(3, 2, 2), "round up");
    same(v(-5, 4, 2, R::downward), v(-3, 2, 2), "negative downward");
    same(v(-5, 4, 2, R::upward), v(-1, 1, 2), "negative upward");
    same(v(-5, 4, 2, R::toward_zero), v(-1, 1, 2), "negative toward zero");
    raises < std::invalid_argument >([] {
        v(5, 3, 7, static_cast<R>(99));
    }, "invalid rounding");
    raises < std::overflow_error >([] {
        Z n(7), d(1); K::from_ratio(n.z, d.z, std::numeric_limits<K::Exponent>::max(), 2, R::nearest);
    }, "overflow");
    auto x = v(-11, 8, 8);
    check(x.precision().bits() == 8 && x.exponent() < 0 && mpz_sgn(x.significand()) > 0, "accessors");
    K y = x;
    same(y, x, "copy");
    K z = K::zero(32);
    z = x;
    same(z, x, "assignment");
    K moved(std::move(y));
    same(moved, x, "move");
    K other = v(3);
    z.swap(other);
    same(other, x, "swap");
    K assigned = K::zero(32);
    assigned = std::move(moved);
    same(assigned, x, "move assign");
    check(K::negate(K::nan(8)).is_nan(), "negate nan");
    check(K::negate(K::zero(8)).negative(), "negate zero");
    raises < std::domain_error >([] {
        Q q; K::infinity(16).to_rational(q.q);
    }, "infinite rational");
    raises < std::domain_error >([] {
        Q q; K::nan(16).to_rational(q.q);
    }, "nan rational");
}
// Coverage-driven, but architecturally meaningful boundary behavior:
// no enormous integer allocation is needed to exercise checked exponent arithmetic.
void boundary_cases() {
    const K::Exponent minimum = std::numeric_limits<K::Exponent>::min();
    const K::Exponent maximum = std::numeric_limits<K::Exponent>::max();
    Z seven(7), one(1), eight(8);
    raises < std::overflow_error >([& ] {
        (void) K::from_ratio(seven.z, one.z, maximum, 2, R::nearest);
    }, "positive exponent overflow");
    raises < std::overflow_error >([& ] {
        (void) K::from_ratio(one.z, eight.z, minimum, 7, R::nearest);
    }, "negative exponent overflow");
    // A finite p=2 number with exponent INT64_MIN can be represented
    // without ever allocating 2^63 bits. It exercises checked_sub's
    // INT64_MIN branch through real arithmetic, not a coverage hook.
    Z two(2);
    const K tiny = K::from_ratio(two.z, one.z, minimum, 2, R::nearest);
    check(tiny.exponent() == minimum && mpz_cmp_ui(tiny.significand(), 2) == 0
    && !tiny.negative(), "smallest exponent representation");
    raises < std::overflow_error >([& ] {
        Q exact;
        tiny.to_rational(exact.q);
    }, "exact rational conversion requires unrepresentable shift");

    // 1/(2 * 2^INT64_MIN) = 2^INT64_MAX, represented without huge storage.
    const K unit = K::from_integer(one.z, 2, R::nearest);
    const K huge = K::div(unit, tiny, 2, R::nearest);
    check(huge.exponent() == maximum - 1
    && mpz_cmp_ui(huge.significand(), 2) == 0 && !huge.negative(),
    "division by minimum-exponent number");

    raises < std::overflow_error >([& ] {
        const K twice = K::from_integer(two.z, 2, R::nearest);
        (void) K::div(twice, tiny, 2, R::nearest);
    }, "division exponent overflow");

    if constexpr(std::numeric_limits<unsigned long>::max() >
    static_cast<unsigned long>(std::numeric_limits<K::Exponent>::max())) {
        raises < std::invalid_argument >([] {
            (void) K(std::numeric_limits<unsigned long>::max());
        }, "precision too large for exponent");
    }

    const K nonzero = v(3), positive_zero = K::zero(16);
    check(!nonzero.is_zero(), "nonzero finite is not zero");
    check(positive_zero.is_zero(), "positive zero classification");
    check(!K::nan(16).is_zero(), "NaN is not zero");

    const K positive_infinity = K::infinity(16);
    const K negative_infinity = K::infinity(16, true);
    same(K::fma(positive_infinity, nonzero, positive_infinity, 16, R::nearest),
    positive_infinity, "fma positive infinity plus positive infinity");
    same(K::fma(negative_infinity, nonzero, negative_infinity, 16, R::nearest),
    negative_infinity, "fma negative infinity plus negative infinity");
}

void specials() {
    auto p = K::zero(16), n = K::zero(16, true), x = v(3), neg = v(-3);
    auto inf = K::infinity(16), minus_inf = K::infinity(16, true), nan = K::nan(16);
    check(p.is_finite() && !p.is_inf() && !p.is_nan(), "finite classification");
    check(inf.is_inf() && !inf.is_zero(), "inf classification");
    check(nan.is_nan() && !nan.is_finite(), "nan classification");
    check(K::add(nan, x, 16, R::nearest).is_nan(), "add lhs nan");
    check(K::add(x, nan, 16, R::nearest).is_nan(), "add rhs nan");
    same(K::add(inf, x, 16, R::nearest), inf, "add inf lhs");
    same(K::add(x, minus_inf, 16, R::nearest), minus_inf, "add inf rhs");
    same(K::add(inf, inf, 16, R::nearest), inf, "add inf inf");
    check(K::add(inf, minus_inf, 16, R::nearest).is_nan(), "infinity cancellation");
    same(K::add(n, n, 16, R::nearest), n, "negative zero");
    same(K::add(p, n, 16, R::downward), n, "downward zero");
    same(K::add(x, neg, 16, R::downward), n, "negative cancellation zero");
    same(K::add(x, neg, 16, R::upward), p, "positive cancellation zero");
    check(K::mul(nan, x, 16, R::nearest).is_nan(), "mul lhs nan");
    check(K::mul(x, nan, 16, R::nearest).is_nan(), "mul rhs nan");
    check(K::mul(inf, p, 16, R::nearest).is_nan(), "infinity times zero");
    check(K::mul(p, inf, 16, R::nearest).is_nan(), "zero times infinity");
    same(K::mul(inf, neg, 16, R::nearest), minus_inf, "mul inf");
    same(K::mul(neg, inf, 16, R::nearest), minus_inf, "mul inf rhs");
    same(K::mul(n, x, 16, R::nearest), n, "mul negative zero");
    same(K::mul(x, n, 16, R::nearest), n, "mul rhs zero");
    check(K::div(nan, x, 16, R::nearest).is_nan(), "div lhs nan");
    check(K::div(x, nan, 16, R::nearest).is_nan(), "div rhs nan");
    check(K::div(p, n, 16, R::nearest).is_nan(), "0/0");
    check(K::div(inf, minus_inf, 16, R::nearest).is_nan(), "inf/inf");
    same(K::div(x, n, 16, R::nearest), minus_inf, "finite / negative zero");
    same(K::div(minus_inf, x, 16, R::nearest), minus_inf, "infinity/finite");
    same(K::div(n, x, 16, R::nearest), n, "zero/finite");
    same(K::div(x, minus_inf, 16, R::nearest), n, "finite/infinity");
    check(K::fma(nan, x, x, 16, R::nearest).is_nan(), "fma lhs nan");
    check(K::fma(x, nan, x, 16, R::nearest).is_nan(), "fma rhs nan");
    check(K::fma(x, x, nan, 16, R::nearest).is_nan(), "fma addend nan");
    check(K::fma(p, inf, x, 16, R::nearest).is_nan(), "fma zero infinity");
    check(K::fma(inf, p, x, 16, R::nearest).is_nan(), "fma infinity zero");
    check(K::fma(inf, x, minus_inf, 16, R::nearest).is_nan(), "fma inf cancellation");
    same(K::fma(inf, x, x, 16, R::nearest), inf, "fma inf lhs");
    same(K::fma(x, inf, x, 16, R::nearest), inf, "fma inf rhs");
    same(K::fma(x, x, minus_inf, 16, R::nearest), minus_inf, "fma infinite addend");
    same(K::fma(p, x, x, 16, R::nearest), x, "fma zero lhs");
    same(K::fma(x, p, x, 16, R::nearest), x, "fma zero rhs");
    same(K::fma(n, p, n, 16, R::nearest), n, "fma signed zero");
}
}
int main() {
    try {
        precision_tag();
        exponent_bit_length_boundaries();
        conversions();
        int64_boundaries();
        rational_integration();
        specials();
        boundary_cases();
        public_comparisons();
        public_rounding();
        std::cout << "FloatAP edge cases OK\n";
        return 0;
    }
    catch (std::exception const& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}

