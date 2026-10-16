/***************************************************************************
 *            test_floatap.cpp
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

// GMP-only FloatAP tests: exact rational oracle, no MPFR required.
#include "numeric/floatap.hpp"
#include <gmp.h>
#include <iostream>
#include <stdexcept>
#include <random>
using Ariadne::FloatAP;
namespace {
using R = FloatAP::Rounding;
struct Q {
    mpq_t value;
    Q() {
        mpq_init(value);
    }
    ~Q() {
        mpq_clear(value);
    }
    Q(Q const& ) = delete;
    Q& operator=(Q const& ) = delete;
};
void expect(bool test, char const* msg) {
    if (!test) throw std::runtime_error(msg);
}
void equal(FloatAP const& x, FloatAP const& y, char const* msg) {
    expect(x.kind() == y.kind(), msg);
    if (!x.is_finite()) {
        expect(x.negative() == y.negative(), msg);
        return;
    }
    Q a, b;
    x.to_rational(a.value);
    y.to_rational(b.value);
    expect(mpq_cmp(a.value, b.value) == 0, msg);
    if (x.is_zero()) expect(x.negative() == y.negative(), msg);
}
void test_rounding() {
    const R modes[] = { R::nearest, R::downward, R::upward, R::toward_zero };
    std::mt19937_64 engine(0x5A17ACEu);
    for (unsigned p : { 7u, 53u, 113u, 128u, 256u, 1024u }) {
        for (unsigned i = 0; i < 500;++i) {
            Q a, b, c, x, y, z, oracle;
            long aa = static_cast<long>(engine() % 200000) - 100000;
            long bb = static_cast<long>(engine() % 200000) - 100000;
            if (bb == 0) bb = 1;
            mpq_set_si(a.value, aa, static_cast<unsigned long>(engine() % 300 + 1));
            mpq_set_si(b.value, bb, static_cast<unsigned long>(engine() % 300 + 1));
            mpq_set_si(c.value, static_cast<long>(engine() % 200000) - 100000,
            static_cast<unsigned long>(engine() % 300 + 1));
            mpq_canonicalize(a.value);
            mpq_canonicalize(b.value);
            mpq_canonicalize(c.value);
            for (R mode : modes) {
                auto fa = FloatAP::from_rational(a.value, p, mode);
                auto fb = FloatAP::from_rational(b.value, p, mode);
                auto fc = FloatAP::from_rational(c.value, p, mode);
                fa.to_rational(x.value);
                fb.to_rational(y.value);
                fc.to_rational(z.value);
                auto check =[& ](FloatAP const& computed, mpq_srcptr exact, char const* message) {
                    auto reference = FloatAP::from_rational(exact, p, mode);
                    equal(computed, reference, message);
                };
                mpq_add(oracle.value, x.value, y.value);
                check(FloatAP::add(fa, fb, p, mode), oracle.value, "add");
                mpq_sub(oracle.value, x.value, y.value);
                check(FloatAP::sub(fa, fb, p, mode), oracle.value, "sub");
                mpq_mul(oracle.value, x.value, y.value);
                check(FloatAP::mul(fa, fb, p, mode), oracle.value, "mul");
                if (mpq_sgn(y.value) != 0) {
                    mpq_div(oracle.value, x.value, y.value);
                    check(FloatAP::div(fa, fb, p, mode), oracle.value, "div");
                }
                mpq_mul(oracle.value, x.value, y.value);
                mpq_add(oracle.value, oracle.value, z.value);
                check(FloatAP::fma(fa, fb, fc, p, mode), oracle.value, "fma");
            }
        }
    }
}
void test_comparison_against_exact_rationals() {
    std::mt19937_64 engine(0xC0A4E5u);
    for (unsigned precision : {7u, 53u, 128u, 256u, 1024u}) {
        for (unsigned i = 0; i < 600; ++i) {
            const long numerator_a = static_cast<long>(engine() % 200001u) - 100000L;
            const long numerator_b = static_cast<long>(engine() % 200001u) - 100000L;
            Q rational_a, rational_b, exact_a, exact_b;
            mpq_set_si(rational_a.value, numerator_a, 1u);
            mpq_set_si(rational_b.value, numerator_b, 1u);
            mpq_canonicalize(rational_a.value);
            mpq_canonicalize(rational_b.value);

            // Different scales exercise leading-bit comparisons, discarded
            // low bits, signed zeros, and equal values with different
            // significand/exponent representations.
            const unsigned shift_a = static_cast<unsigned>(engine() % 200u);
            const unsigned shift_b = static_cast<unsigned>(engine() % 200u);
            mpz_mul_2exp(mpq_denref(rational_a.value),
                         mpq_denref(rational_a.value), shift_a);
            mpz_mul_2exp(mpq_denref(rational_b.value),
                         mpq_denref(rational_b.value), shift_b);
            mpq_canonicalize(rational_a.value);
            mpq_canonicalize(rational_b.value);
            const FloatAP a = FloatAP::from_rational(rational_a.value, precision, R::nearest);
            const FloatAP b = FloatAP::from_rational(rational_b.value, precision, R::nearest);
            a.to_rational(exact_a.value);
            b.to_rational(exact_b.value);
            const int oracle = mpq_cmp(exact_a.value, exact_b.value);
            const auto result = cmp(a, b);
            expect((oracle < 0 && result == Ariadne::Comparison::LESS) ||
                   (oracle == 0 && result == Ariadne::Comparison::EQUAL) ||
                   (oracle > 0 && result == Ariadne::Comparison::GREATER),
                   "direct comparison differs from exact rational oracle");
        }
    }
    const FloatAP zero = FloatAP::zero(53);
    const FloatAP negative_zero = FloatAP::zero(53, true);
    expect(cmp(zero, negative_zero) == Ariadne::Comparison::EQUAL,
           "signed zeros compare equal");
    Q exact_one;
    mpq_set_ui(exact_one.value, 1, 1);
    const FloatAP one = FloatAP::from_rational(exact_one.value, 53, R::nearest);
    expect(cmp(one, one) == Ariadne::Comparison::EQUAL, "equal operands");
}
void test_special() {
    using K = FloatAP;
    auto plus = K::zero(64), minus = K::zero(64, true);
    auto inf = K::infinity(64), ninf = K::infinity(64, true), nan = K::nan(64);
    expect(plus.is_zero() && !plus.negative(), "positive zero");
    expect(minus.is_zero() && minus.negative(), "negative zero");
    expect(K::add(minus, minus, 64, R::nearest).negative(), "signed zero addition");
    expect(K::add(plus, minus, 64, R::downward).negative(), "downward zero");
    expect(K::mul(minus, inf, 64, R::nearest).is_nan(), "0 times infinity");
    expect(K::add(inf, ninf, 64, R::nearest).is_nan(), "opposite infinities");
    expect(K::div(plus, plus, 64, R::nearest).is_nan(), "0/0");
    expect(K::div(inf, inf, 64, R::nearest).is_nan(), "inf/inf");
    expect(K::add(nan, plus, 64, R::nearest).is_nan(), "NaN propagation");
    expect(K::div(inf, minus, 64, R::nearest).negative(), "signed infinity");
    bool rejected = false;
    try {
        K(1);
    }
    catch (std::invalid_argument const& ) {
        rejected = true;
    }
    expect(rejected, "precision lower bound");
}
}
int main() {
    try {
        test_special();
        test_rounding();
        test_comparison_against_exact_rationals();
        std::cout << "FloatAP: 3000 randomized inputs x four rounding modes OK\n";
    }
    catch (std::exception const& e) {
        std::cerr << "FloatAP test FAILED: " << e.what() << '\n';
        return 1;
    }
}

