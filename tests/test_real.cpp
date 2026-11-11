/***************************************************************************
 *            test_real.cpp
 *
 *  Copyright  2013-20  Pieter Collins
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

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "foundation/logical.hpp"
#include "numeric/accuracy.hpp"
#include "numeric/builtin.hpp"
#include "numeric/decimal.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_literals.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/integer.hpp"
#include "numeric/rational.hpp"
#include "numeric/reals.hpp"
#include "utility/module.hpp"

namespace Ariadne {

Rational to_rational(String x);
Rational operator""_q(const char* str, std::size_t);

Rational to_rational(String x) {
    Rational q = 0;
    SizeType i = 0;
    while (x[i] != char(0) && x[i] != '.') {
        short d = x[i] - '0';
        q = 10 * q + d;
        ++i;
    }
    if (x[i] == '.') {
        ++i;
    }
    Rational r = 1;
    while (x[i] != char(0)) {
        short d = x[i] - '0';
        r /= 10;
        q = q + d * r;
        ++i;
    }
    return q;
}

Rational operator""_q(const char* str, std::size_t) { return Rational(Decimal(String(str))); }

// Boolean nondeterministic_greater(Real const& x, Rational const& a, Rational const& b);
// template<class Q> Boolean nondeterministic_greater(Real const& x, Q a, Q b) { return nondeterministic_greater(x,Rational(a),Rational(b)); }

} // namespace Ariadne

#include "numeric/floats.hpp"

#include "ariadne-testing.hpp"

using namespace std;
using namespace Ariadne;

namespace {

constexpr bool check_real_concept() {
    return requires(Real x, Real y, Nat m, Int n) {
        +x;
        -x;
        add(x, y);
        sub(x, y);
        mul(x, y);
        div(x, y);
        pow(x, m);
        pow(x, n);
        pos(x);
        neg(x);
        hlf(x);
        sqr(x);
        rec(x);
        sqrt(x);
        exp(x);
        log(x);
        sin(x);
        cos(x);
        tan(x);
        tanh(x);
        asin(x);
        acos(x);
        atan(x);
        abs(x);
        max(x, y);
        min(x, y);
    };
}

static_assert(check_real_concept());
#ifndef ARIADNE_TEST_DATA_DIR
#error "Missing test data directory"
#endif
} // namespace

class TestReal : public Ariadne::Testing::ContextFixture {
  private:
    struct CertifiedPi {
        Rational lower, upper, near;
    };
    static CertifiedPi certified_pi() {
        std::ifstream file(ARIADNE_TEST_DATA_DIR "/floatmp_pi.tsv");
        if (!file)
            throw std::runtime_error("Missing certified pi oracle");
        std::string line;
        CertifiedPi result;
        bool found = false;
        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '#')
                continue;
            std::istringstream row(line);
            std::string b, p, ln, ld, un, ud, nn, nd, extra;
            if (!(row >> b >> p >> ln >> ld >> un >> ud >> nn >> nd) || (row >> extra) || b != "384" || p != "320" || found)
                throw std::runtime_error("Malformed pi oracle");
            result = {Rational(Integer(String(ln)), Integer(String(ld))),
                      Rational(Integer(String(un)), Integer(String(ud))),
                      Rational(Integer(String(nn)), Integer(String(nd)))};
            if (!(result.lower < result.upper && result.near - Dyadic(1, 319u) < result.lower && result.upper < result.near + Dyadic(1, 319u)))
                throw std::runtime_error("Invalid pi interval");
            found = true;
        }
        if (!found)
            throw std::runtime_error("Empty certified pi oracle");
        return result;
    }

  public:
    using ContextFixture::ContextFixture;
    void test_precision_preflight();
    void test_expression_evaluation();
    void test_constructors();
    void test_conversions();
    void test_arithmetic();
    void test_directed();
    void test_transcendental();
    void test_comparison();
    void test_rounding();
    void test_accuracy();
    void test_sequence();

  private:
    // Equality of computable reals is not semidecidably provable. Instead
    // verify an independently certified enclosure at a stated resolution.
    void expect_certified_value(const Real& actual, const Dyadic& expected,
                                std::string description) {
        [[maybe_unused]] auto value_scope = this->ctx().scope(description);
        const DyadicBounds bounds = actual.compute(Accuracy(64_bits)).get();
        ARIADNE_TESTING_PROVE(bounds.lower_raw() <= expected);
        ARIADNE_TESTING_PROVE(expected <= bounds.upper_raw());
        ARIADNE_TESTING_PROVE(bounds.upper_raw() - bounds.lower_raw() <= Dyadic(1, 40u));
    }

    // Real-valued distances cannot generally be decided. Instead certify
    // both endpoints against the original exact dyadic tolerance.
    void expect_certified_within(const Real& actual,
                                 const FloatDPApproximation& reference,
                                 const Dyadic& tolerance,
                                 std::string description) {
        [[maybe_unused]] auto value_scope = this->ctx().scope(description);
        const Dyadic exact_reference(reference.raw());
        const DyadicBounds enclosure =
            actual.compute(Accuracy(128_bits)).get();
        ARIADNE_TESTING_PROVE(enclosure.lower_raw() >= exact_reference - tolerance);
        ARIADNE_TESTING_PROVE(enclosure.upper_raw() <= exact_reference + tolerance);
    }
};

void TestReal::test_precision_preflight() {
    // The precision is validated before any significand storage is allocated.
    if constexpr (std::numeric_limits<std::size_t>::max() > std::numeric_limits<Nat>::max()) {
        const MultiplePrecision excessive(static_cast<std::size_t>(std::numeric_limits<Nat>::max()) + 1u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(Real(1).get(excessive));
    }
}

void TestReal::test_expression_evaluation() {
    FloatDPApproximation::set_output_places(18);
    Effort eff(8u);
    MultiplePrecision mp(128_bits);
    Real x = Real(1) / 2;
    Real y = Real(2);

    std::vector<Real> expressions{
        +x, -x,
        add(x, y), sub(x, y), mul(x, y), div(x, y),
        pow(x, 2u), pow(x, -2),
        pos(x), neg(x), hlf(x), sqr(x), rec(x),
        sqrt(x), exp(x), log(x),
        sin(x), cos(x), tan(x), tanh(x),
        asin(x), acos(x), atan(x),
        abs(x), max(x, y), min(x, y)};

    for (const Real& expression: expressions) {
        FloatDPBounds dp_bounds = expression.compute_get(eff, dp);
        FloatMPBounds ap_bounds = expression.compute_get(eff, mp);
        ARIADNE_TESTING_PROVE(dp_bounds.lower_raw() <= dp_bounds.upper_raw());
        ARIADNE_TESTING_PROVE(ap_bounds.lower_raw() <= ap_bounds.upper_raw());

        std::ostringstream expression_stream;
        expression_stream << expression;
        ARIADNE_TESTING_PROVE(not expression_stream.str().empty());
    }
}

void TestReal::test_conversions() {
    Real one = 1;
    Real pi_ = 4 * atan(one);
    MultiplePrecision mp(128_bits);
    auto pi_dp = pi_.get(dp);
    ARIADNE_TESTING_INSPECT("pi_dp", pi_dp);

    ARIADNE_TESTING_EXPECT_EQUAL(one.get_d(), 1.0);
    ARIADNE_TESTING_EXPECT_EQUAL(one.compute_using(dp).lower_raw(), FloatDP(1, dp));
    ARIADNE_TESTING_EXPECT_EQUAL(one.compute_using(mp).lower_raw(), FloatMP(1, mp));
    const FloatMPBounds native = one.get(mp);
    ARIADNE_TESTING_EXPECT_EQUAL(native.lower_raw(), FloatMP(1, mp));
    ARIADNE_TESTING_EXPECT_EQUAL(native.upper_raw(), FloatMP(1, mp));
    ARIADNE_TESTING_EXPECT_EQUAL(FloatMPLowerBound(native).raw(), FloatMP(1, mp));
    ARIADNE_TESTING_EXPECT_EQUAL(FloatMPUpperBound(native).raw(), FloatMP(1, mp));
    ARIADNE_TESTING_EXPECT_EQUAL(one.compute_using(mp).lower_raw(), FloatMP(1, mp));
    ARIADNE_TESTING_EXPECT_EQUAL(one.compute_get(Effort(8u), mp).upper_raw(), FloatMP(1, mp));
    ARIADNE_TESTING_EXPECT_EQUAL(ValidatedReal(DyadicBounds(Dyadic(1), Dyadic(2)))
                                     .get(mp)
                                     .lower_raw(),
                                 FloatMP(1, mp));

    ARIADNE_TESTING_EXPECT_NO_THROW(FloatDPBall(one, dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(FloatDPBounds(one, dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(FloatDPUpperBound(one, dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(FloatDPLowerBound(one, dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(FloatDPApproximation(one, dp));

    DyadicBounds dyadic_bounds(Dyadic(1), Dyadic(2));
    ARIADNE_TESTING_EXPECT_EQUAL(dyadic_bounds.get(dp).lower_raw(), FloatDP(1, dp));
    ARIADNE_TESTING_EXPECT_EQUAL(dyadic_bounds.get(dp).upper_raw(), FloatDP(2, dp));
    ARIADNE_TESTING_EXPECT_EQUAL(dyadic_bounds.get(mp).lower_raw(), FloatMP(1, mp));
    ARIADNE_TESTING_EXPECT_EQUAL(dyadic_bounds.get(mp).upper_raw(), FloatMP(2, mp));

    std::ostringstream accuracy_stream;
    accuracy_stream << Accuracy(8_bits);
    ARIADNE_TESTING_PROVE(not accuracy_stream.str().empty());

    ARIADNE_TESTING_EXPECT_EQUAL(class_name<Real>(), String("Real"));
    ARIADNE_TESTING_EXPECT_EQUAL(class_name<PositiveReal>(), String("PositiveReal"));

    ApproximateDouble approximate_one(one);
    ARIADNE_TESTING_EXPECT_EQUAL(approximate_one.get_d(), 1.0);
}

void TestReal::test_constructors() {
    Effort eff(2);
    DoublePrecision pr;
    Real xv;
    ARIADNE_TESTING_INSPECT("constructed xv", xv);
    ARIADNE_TESTING_EXPECT_EQUAL(xv.get(pr), 0);
    ARIADNE_TESTING_EXPECT_EQUAL(xv.lower().compute(eff).get().raw(), 0);
    ARIADNE_TESTING_EXPECT_EQUAL(xv.upper().compute(eff).get().raw(), 0);

    NaiveReal naive_zero;
    Real naive_zero_as_real(naive_zero.managed_pointer());
    ARIADNE_TESTING_EXPECT_EQUAL(naive_zero_as_real.get(pr), 0);

    PositiveNaiveReal positive_naive_zero;
    Real positive_naive_zero_as_real(positive_naive_zero.managed_pointer());
    ARIADNE_TESTING_EXPECT_EQUAL(positive_naive_zero_as_real.get(pr), 0);

    PositiveNaiveReal positive_naive_copy(naive_zero);
    ARIADNE_TESTING_EXPECT_EQUAL(Real(positive_naive_copy.managed_pointer()).get(pr), 0);

    PositiveReal positive_real_zero;
    PositiveReal positive_real_one(Real(1));
    ARIADNE_TESTING_EXPECT_EQUAL(positive_real_zero.get(pr), 0);
    ARIADNE_TESTING_EXPECT_EQUAL(positive_real_one.get(pr), 1);

    ARIADNE_TESTING_EXPECT_NO_THROW(PositiveLowerReal{positive_real_one});
    ARIADNE_TESTING_EXPECT_NO_THROW(PositiveLowerReal{Real(1).lower()});
    ARIADNE_TESTING_EXPECT_NO_THROW(PositiveUpperReal{positive_real_one});
    ARIADNE_TESTING_EXPECT_NO_THROW(PositiveUpperReal{Real(1).upper()});
    Real xz(1);
    ARIADNE_TESTING_INSPECT("constructed xz", xz);
    ARIADNE_TESTING_EXPECT_EQUAL(xz.compute(eff).get(), 1);
    ARIADNE_TESTING_EXPECT_NO_THROW(Real(std::uint64_t(2)));
    ARIADNE_TESTING_EXPECT_NO_THROW(Real(std::int64_t(-2)));
    ARIADNE_TESTING_EXPECT_NO_THROW(Real(Dyadic(3, 1u)));
    ARIADNE_TESTING_EXPECT_NO_THROW(Real(Decimal("1.25")));
    ARIADNE_TESTING_EXPECT_NO_THROW(Real(Rational(5, 4)));
    ARIADNE_TESTING_EXPECT_NO_THROW(Real(EffectiveNumber(2)));
    const FloatDP binary_float(1, pr);
    ARIADNE_TESTING_EXPECT_NO_THROW(Real(Dyadic(binary_float)));
    ARIADNE_TESTING_EXPECT_NO_THROW(Real(ExactDouble(binary_float.get_d())));
    Real xe{Dyadic(1.5_exact)};
    ARIADNE_TESTING_INSPECT("constructed xe", xe);
    ARIADNE_TESTING_NOTE("before xe.compute(eff).get()");
    auto xe_bounds = xe.compute(eff).get();
    ARIADNE_TESTING_NOTE("after xe.compute(eff).get()");
    ARIADNE_TESTING_EXPECT_EQUAL(xe_bounds, 1.5_dy);
    ARIADNE_TESTING_NOTE("after xe equality");
    ARIADNE_TESTING_NOTE("before rational literal");
    Rational qn = 1.1_q;
    ARIADNE_TESTING_NOTE("after rational literal: " << qn);
    ARIADNE_TESTING_NOTE("before Real xn(qn)");
    Real xn(qn);
    ARIADNE_TESTING_NOTE("after Real xn(qn): " << xn);
    ARIADNE_TESTING_EXPECT_LESS(Rational(xn.lower().compute(eff).get().raw()), Rational(11, 10));
    ARIADNE_TESTING_EXPECT_GREATER(Rational(xn.upper().compute(eff).get().raw()), Rational(11, 10));
    Real xq(Rational(11, 10));
    ARIADNE_TESTING_INSPECT("constructed xq", xq);
    ARIADNE_TESTING_EXPECT_LESS(Rational(xq.lower().compute(eff).get().raw()), Rational(11, 10));
    ARIADNE_TESTING_EXPECT_GREATER(Rational(xq.upper().compute_get(eff).raw()), Rational(11, 10));

    Real xy(EffectiveNumber(5));
    ARIADNE_TESTING_INSPECT("constructed xy", xy);
    ARIADNE_TESTING_EXPECT_EQUAL(xy.get(double_precision), 5);
    ARIADNE_TESTING_EXPECT_EQUAL(xy.get(MultiplePrecision(128_bits)), 5);
}

void TestReal::test_arithmetic() {
    FloatDPApproximation::set_output_places(18);
    Real x(2.5_dyadic);
    Real y(4.0_dyadic);
    expect_certified_value(x, Dyadic(2.5_dy), "x certified within 2^-40 of 2.5_dy");
    expect_certified_value(y, Dyadic(4.0_x), "y certified within 2^-40 of 4.0_x");
    expect_certified_value(+x, Dyadic(2.5_x), "+x certified within 2^-40 of 2.5_x");
    expect_certified_value(-x, Dyadic(-2.5_x), "-x certified within 2^-40 of -2.5_x");
    expect_certified_value(x + y, Dyadic(6.5_x), "x+y certified within 2^-40 of 6.5_x");
    expect_certified_value(x - y, Dyadic(-1.5_x), "x-y certified within 2^-40 of -1.5_x");
    expect_certified_value(x * y, Dyadic(10.0_x), "x*y certified within 2^-40 of 10.0_x");
    expect_certified_value(x / y, Dyadic(0.625_dy), "x/y certified within 2^-40 of 0.625_dy");
    expect_certified_value(add(x, y), Dyadic(6.5_dy), "add(x,y) certified within 2^-40 of 6.5_dy");
    expect_certified_value(sub(x, y), Dyadic(-1.5_dy), "sub(x,y) certified within 2^-40 of -1.5_dy");
    expect_certified_value(mul(x, y), Dyadic(10.0_dy), "mul(x,y) certified within 2^-40 of 10.0_dy");
    expect_certified_value(div(x, y), Dyadic(0.625_dy), "div(x,y) certified within 2^-40 of 0.625_dy");
    expect_certified_value(pow(x, 3u), Dyadic(15.625_dy), "pow(x,3u) certified within 2^-40 of 15.625_dy");
    expect_certified_value(pow(x, 3), Dyadic(15.625_dy), "pow(x,3) certified within 2^-40 of 15.625_dy");
    expect_certified_value(pos(x), Dyadic(+2.5_dy), "pos(x) certified within 2^-40 of +2.5_dy");
    expect_certified_value(neg(x), Dyadic(-2.5_dy), "neg(x) certified within 2^-40 of -2.5_dy");
    expect_certified_value(hlf(x), Dyadic(1.25_dy), "hlf(x) certified within 2^-40 of 1.25_dy");
    expect_certified_value(sqr(x), Dyadic(6.25_dy), "sqr(x) certified within 2^-40 of 6.25_dy");
    expect_certified_value(rec(y), Dyadic(0.25_dy), "rec(y) certified within 2^-40 of 0.25_dy");
    expect_certified_value(nul(x), Dyadic(0.0_dy), "nul(x) certified within 2^-40 of 0.0_dy");
    ARIADNE_TESTING_PROVE(same(Real(2), Real(2)));
    ARIADNE_TESTING_PROVE(not same(Real(2), Real(3)));
    ARIADNE_TESTING_EXPECT_NO_THROW(mag(Real(-2)).compute_get(Effort(4u)));
    ARIADNE_TESTING_EXPECT_NO_THROW(mag(Real(-2), dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(dist(Real(2), Real(3)).compute_get(Effort(4u)));
}

void TestReal::test_directed() {
    Effort eff(4u);
    MultiplePrecision mp(128_bits);
    Real two_real = 2;
    Real three_real = 3;

    LowerReal ltwo(two_real);
    LowerReal lthree(three_real);
    UpperReal utwo(two_real);
    UpperReal uthree(three_real);

    ARIADNE_TESTING_EXPECT_NO_THROW(ltwo.compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(ltwo.compute_get(eff, dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(ltwo.compute_get(eff, mp));
    ARIADNE_TESTING_PROVE(ltwo.compute_get(eff, mp).raw() <= FloatMP(2, mp));
    ARIADNE_TESTING_EXPECT_NO_THROW(utwo.compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(utwo.compute_get(eff, dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(utwo.compute_get(eff, mp));
    ARIADNE_TESTING_PROVE(utwo.compute_get(eff, mp).raw() >= FloatMP(2, mp));

    std::ostringstream directed_stream;
    directed_stream << ltwo << " " << utwo;
    ARIADNE_TESTING_PROVE(not directed_stream.str().empty());

    ARIADNE_TESTING_EXPECT_NO_THROW(max(ltwo, lthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(min(ltwo, lthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(min(ltwo, three_real).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(min(two_real, lthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(max(utwo, uthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(max(utwo, three_real).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(max(two_real, uthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(min(utwo, uthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(neg(utwo).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(neg(ltwo).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(add(ltwo, lthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(add(utwo, uthree).compute_get(eff));

    PositiveReal ptwo(two_real);
    PositiveReal pthree(three_real);
    ARIADNE_TESTING_EXPECT_NO_THROW(ptwo.compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(max(ptwo, pthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(min(ptwo, pthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(add(ptwo, pthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(mul(ptwo, pthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(div(ptwo, pthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(rec(ptwo).compute_get(eff));

    PositiveLowerReal pltwo(ptwo);
    PositiveLowerReal plthree(pthree);
    PositiveUpperReal putwo(ptwo);
    PositiveUpperReal puthree(pthree);
    ARIADNE_TESTING_EXPECT_NO_THROW(pltwo.compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(putwo.compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(rec(pltwo).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(rec(putwo).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(add(pltwo, plthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(add(putwo, puthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(mul(pltwo, plthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(mul(putwo, puthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(div(pltwo, puthree).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(div(putwo, plthree).compute_get(eff));

    ARIADNE_TESTING_EXPECT_NO_THROW(mul(ltwo, ptwo).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(mul(utwo, ptwo).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(mul(ptwo, ltwo).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(mul(ptwo, utwo).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(div(ltwo, ptwo).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(div(utwo, ptwo).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(div(ptwo, utwo).compute_get(eff));
    ARIADNE_TESTING_EXPECT_NO_THROW(div(ptwo, ltwo).compute_get(eff));

    ValidatedLowerReal validated_lower{DyadicLowerBound{Dyadic(2)}};
    ARIADNE_TESTING_EXPECT_NO_THROW(validated_lower.get());
    ARIADNE_TESTING_EXPECT_NO_THROW(validated_lower.get(dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(validated_lower.get(mp));
    ValidatedUpperReal validated_upper{DyadicUpperBound{Dyadic(2)}};
    ARIADNE_TESTING_EXPECT_NO_THROW(validated_upper.get());
    ARIADNE_TESTING_EXPECT_NO_THROW(validated_upper.get(dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(validated_upper.get(mp));
    std::ostringstream validated_stream;
    validated_stream << validated_lower << " " << validated_upper;
    ARIADNE_TESTING_PROVE(not validated_stream.str().empty());
}

void TestReal::test_transcendental() {
    Dyadic eps{FloatDP::eps(dp)};
    Real x(2.5_dyadic);
    FloatDPApproximation ax = x.get(dp);
    expect_certified_value(sqrt(Real(4)), Dyadic(2.0_dy), "sqrt(Real(4)) certified within 2^-40 of 2.0_dy");
    expect_certified_value(exp(Real(0)), Dyadic(1.0_dy), "exp(Real(0)) certified within 2^-40 of 1.0_dy");
    expect_certified_value(log(Real(1)), Dyadic(0.0_dy), "log(Real(1)) certified within 2^-40 of 0.0_dy");

    expect_certified_within(sqrt(x), sqrt(ax), eps,
                            "sqrt(x) within eps of sqrt(ax)");
    expect_certified_within(exp(x), exp(ax), 8 * eps,
                            "exp(x) within 8*eps of exp(ax)");
    expect_certified_within(log(x), log(ax), eps,
                            "log(x) within eps of log(ax)");
    expect_certified_within(sin(x), sin(ax), 2 * eps,
                            "sin(x) within 2*eps of sin(ax)");
    expect_certified_within(cos(x), cos(ax), 2 * eps,
                            "cos(x) within 2*eps of cos(ax)");
    expect_certified_within(tan(x), tan(ax), 3 * eps,
                            "tan(x) within 3*eps of tan(ax)");
    expect_certified_within(atan(x), atan(ax), eps,
                            "atan(x) within eps of atan(ax)");
}

void TestReal::test_comparison() {
    Effort effort(0);
    Real one = 1;
    Real e = exp(one);
    Real elogpi = e * log(pi);
    //    ARIADNE_TEST_SAME(e,e);
    //    ARIADNE_TESTING_EXPECT_EQUAL(log(e),one);
    ARIADNE_TESTING_PROVE(possibly(check(log(e) == one, effort)));
    ARIADNE_TESTING_PROVE(possibly((log(e) == one).check(effort)));

    ARIADNE_TESTING_INSPECT("pi.get(dp)", pi.get(dp));
    ARIADNE_TESTING_INSPECT("e*log(pi).get(dp)", e * log(pi).get(dp));

    ARIADNE_TESTING_PROVE(not check(pi == elogpi, effort));
    ARIADNE_TESTING_PROVE(check(pi != elogpi, effort));
    ARIADNE_TESTING_PROVE(check(elogpi < pi, effort));
    ARIADNE_TESTING_PROVE(check(elogpi <= pi, effort));
    ARIADNE_TESTING_PROVE(not check(elogpi > pi, effort));
    ARIADNE_TESTING_PROVE(not check(elogpi >= pi, effort));

    Rational zero = 0;
    ARIADNE_TESTING_PROVE(check(e * log(pi) - pi < zero, effort));

    ARIADNE_TESTING_EXPECT_EQUAL(check_sgn(Real(1), effort), true);
    ARIADNE_TESTING_EXPECT_EQUAL(check_sgn(Real(-1), effort), false);
    ARIADNE_TESTING_EXPECT_UNARY_PREDICATE(is_indeterminate, check_sgn(Real(0), effort));

    auto named_equal = eq(Real(1), Real(1));
    auto named_less = lt(Real(1), Real(2));
    auto named_sign = sgn(Real(-1));
    ARIADNE_TESTING_PROVE(possibly(named_equal.check(effort)));
    ARIADNE_TESTING_PROVE(definitely(named_less.check(effort)));
    ARIADNE_TESTING_PROVE(definitely(not named_sign.check(effort)));

    auto* sign_copy = named_sign.repr().pointer()->_copy();
    ARIADNE_TESTING_EXPECT_EQUAL(
        static_cast<int>(sign_copy->_check(effort)),
        static_cast<int>(named_sign.repr().check(effort)));
    delete sign_copy;

    std::ostringstream logical_stream;
    logical_stream << named_equal << " " << named_less << " " << named_sign;
    ARIADNE_TESTING_PROVE(not logical_stream.str().empty());
}

void TestReal::test_rounding() {
    ARIADNE_TESTING_EXPECT_EQUAL(round(Real(0.1_dec) * 10), 1);
    ARIADNE_TESTING_EXPECT_EQUAL(round(Real(0.9_dec) * 10), 9);
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(round(Real(-0.7_dec) * 5), -4);
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(round(Real(-0.7_dec) * 5), -3);
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(round(Real(0.1_dec) * 5), 0);
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(round(Real(0.1_dec) * 5), 1);
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(round(Real(0.9_dec) * 5), 4);
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(round(Real(0.9_dec) * 5), 5);
}

void TestReal::test_accuracy() {
    FloatDPBounds::set_output_places(18);
    Real one = 1;
    Real pi_ = 4 * atan(one);

    const CertifiedPi oracle = certified_pi();
    FloatDPBounds pi_dp(pi_.get(dp));
    ARIADNE_TESTING_INSPECT("constructed pi_dp", pi_dp);
    const MultiplePrecision mp(128_bits);
    FloatMPBounds pi_ap(pi_.get(mp));
    ARIADNE_TESTING_INSPECT("constructed pi_ap", pi_ap);
    ARIADNE_TESTING_PROVE(Rational(pi_dp.lower_raw()) <= oracle.lower);
    ARIADNE_TESTING_PROVE(Rational(pi_dp.upper_raw()) >= oracle.upper);
    ARIADNE_TESTING_PROVE(Rational(pi_ap.lower_raw()) <= oracle.lower);
    ARIADNE_TESTING_PROVE(Rational(pi_ap.upper_raw()) >= oracle.upper);

    const MultiplePrecision ap_high(320_bits);
    Effort effort{256};
    ValidatedReal pi_ord(pi_.compute(effort));
    ARIADNE_TESTING_INSPECT("constructed pi_ord", pi_ord);
    Accuracy accuracy{256_bits};
    FloatMP error(accuracy.error(), ap_high);
    ValidatedReal pi_met(pi_.compute(accuracy));
    ARIADNE_TESTING_INSPECT("constructed pi_met", pi_met);
    FloatMPBounds pi_met_ap(pi_met.get(ap_high));
    ARIADNE_TESTING_INSPECT("constructed pi_met_ap", pi_met_ap);
    ARIADNE_TESTING_INSPECT("pi_met_ap.error()", pi_met_ap.error());
    const FloatMP pi_near(oracle.near, FloatMP::ROUND_TO_NEAREST, ap_high);
    ARIADNE_TESTING_EXPECT_EQUAL(Rational(pi_near), oracle.near);
    ARIADNE_TESTING_INSPECT("abs(sub(up,pi_met_ap.value_raw(),pi_near))", abs(sub(up, pi_met_ap.value_raw(), pi_near)));
    ARIADNE_TESTING_PROVE(pi_met_ap.error() <= error);
    ARIADNE_TESTING_PROVE(rad(up, pi_met_ap.value(), pi_near) <= error.raw());

    Dyadic eps(1, 1024u);
    Real sin_pi(sin(pi_));
    ARIADNE_TESTING_INSPECT("constructed sin_pi", sin_pi);
    ARIADNE_TESTING_INSPECT("choose(sin_pi>-eps,sin_pi<eps)", choose(sin_pi > -eps, sin_pi < eps));
    /*
        ARIADNE_TESTING_INSPECT("nondeterministic_greater(sin_pi,-eps,eps)", nondeterministic_greater(sin_pi,-eps,eps));
        ARIADNE_TESTING_INSPECT("nondeterministic_greater(sin_pi,-1,eps)", nondeterministic_greater(sin_pi,-1,eps));
        ARIADNE_TESTING_INSPECT("nondeterministic_greater(sin_pi,-eps,1)", nondeterministic_greater(sin_pi,-eps,1));
        ARIADNE_TESTING_PROVE(not nondeterministic_greater(sin_pi,eps,2*eps));
        ARIADNE_TESTING_PROVE(nondeterministic_greater(sin_pi,-3*eps,-2*eps));
    */
}

void TestReal::test_sequence() {
    std::function<Dyadic(Natural)> wfn([&](Natural n) { return 1 - Dyadic(1, n); });
    FastCauchySequence<Dyadic> wseq(wfn);
    Real wlim = limit(wseq);
    Real wdirect(wseq);
    MultiplePrecision mp(128_bits);
    ARIADNE_TESTING_EXPECT_NO_THROW(wlim.compute_get(Effort(4u), dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(wlim.compute_get(Effort(4u), mp));
    ARIADNE_TESTING_EXPECT_NO_THROW(wdirect.compute_get(Effort(4u), dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(wdirect.compute_get(Effort(4u), mp));
    std::ostringstream dyadic_limit_stream;
    dyadic_limit_stream << wdirect;
    ARIADNE_TESTING_PROVE(not dyadic_limit_stream.str().empty());

    Sequence<DyadicBounds> bounds_sequence(std::function<DyadicBounds(Natural)>(
        [](Natural) { return DyadicBounds(Dyadic(1)); }));
    auto first_bounds = bounds_sequence[Natural(0u)];
    ARIADNE_TESTING_EXPECT_EQUAL(first_bounds.lower_raw(), Dyadic(1));
    ARIADNE_TESTING_EXPECT_EQUAL(first_bounds.upper_raw(), Dyadic(1));
    SequenceWriter writer(2u);
    std::ostringstream sequence_stream;
    sequence_stream << writer(bounds_sequence);
    ARIADNE_TESTING_PROVE(not sequence_stream.str().empty());

    ConvergentSequence<DyadicBounds> convergent_bounds(bounds_sequence);
    Real bounds_limit(convergent_bounds);
    ARIADNE_TESTING_EXPECT_NO_THROW(bounds_limit.compute_get(Effort(2u)));
    ARIADNE_TESTING_EXPECT_NO_THROW(bounds_limit.compute_get(Effort(2u), dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(bounds_limit.compute_get(Effort(2u), mp));
    std::ostringstream bounds_limit_stream;
    bounds_limit_stream << bounds_limit;
    ARIADNE_TESTING_PROVE(not bounds_limit_stream.str().empty());

    ARIADNE_TESTING_INSPECT("wlim.compute(Accuracy(256_bits))", wlim.compute(Accuracy(256_bits)));

    std::function<Real(Natural)> rfn([&](Natural n) { return exp(Real(-(n + 1u))); });
    FastCauchySequence<Real> rseq(rfn);
    Real rlim = limit(rseq);
    ARIADNE_TESTING_EXPECT_NO_THROW(rlim.compute_get(Effort(4u), dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(rlim.compute_get(Effort(4u), mp));
    std::ostringstream real_limit_stream;
    real_limit_stream << rlim;
    ARIADNE_TESTING_PROVE(not real_limit_stream.str().empty());

    ARIADNE_TESTING_PROVE(abs(rlim.compute(Accuracy(256_bits)).get()) < Dyadic(1, 256u));

    Real pi_ = 4 * atan(Real(1));
    Decimal pi_lower = "3.1415926535897932"_dec;
    Decimal pi_upper = "3.1415926535897933"_dec;
    ARIADNE_TESTING_INSPECT("nondeterministic_greater(pi_,\"3.1415926535897932\"_dec,\"3.1415926535897933\"_dec)", nondeterministic_greater(pi_, "3.1415926535897932"_dec, "3.1415926535897933"_dec));
    ARIADNE_TESTING_PROVE(abs(rlim.compute(Accuracy(256_bits)).get()) < Dyadic(1, 256u));
    ARIADNE_TESTING_PROVE(nondeterministic_greater(pi_, 3.0_dec, 3.1_dec));
    ARIADNE_TESTING_PROVE(not nondeterministic_greater(pi_, 3.2_dec, 3.3_dec));
    ARIADNE_TESTING_EXPECT_EXCEPTION(nondeterministic_greater(pi_, Rational(1), Rational(1)));
    ARIADNE_TESTING_INSPECT("nondeterministic_greater(pi_,\"3.1415926535897932\"_dec,\"3.1415926535897937\"_dec)", nondeterministic_greater(pi_, "3.1415926535897932"_dec, "3.1415926535897937"_dec));
    ARIADNE_TESTING_INSPECT("nondeterministic_greater(pi_,\"3.1415926535897932\"_dec,\"3.1415926535897933\"_dec)", nondeterministic_greater(pi_, "3.1415926535897932"_dec, "3.1415926535897933"_dec));

    ARIADNE_TESTING_INSPECT("pi_.compute(Accuracy(256_bits)).get(MultiplePrecision(256_bits))", pi_.compute(Accuracy(256_bits)).get(MultiplePrecision(256_bits)));
    std::cout << std::boolalpha;
    ARIADNE_TESTING_INSPECT("choose(pi_>pi_lower,pi_<pi_upper)", choose(pi_ > pi_lower, pi_ < pi_upper));
    ARIADNE_TESTING_INSPECT("choose(pi_>pi_lower,pi_<3.2_dec)", choose(pi_ > pi_lower, pi_ < 3.2_dec));

    ARIADNE_TESTING_INSPECT("choose({pi_>pi_lower,Real(4)},{pi_<3.2_dec,Real(3)})", choose({pi_ > pi_lower, Real(4)}, {pi_ < 3.2_dec, Real(3)}));
    ARIADNE_TESTING_INSPECT("choose({pi_>pi_lower,Real(4)},{pi_<pi_upper,Real(3)})", choose({pi_ > pi_lower, Real(4)}, {pi_ < pi_upper, Real(3)}));

    Dyadic pi_upper_bound = pi_.compute(Accuracy(256_bits)).get().upper_raw();
    Real x(pi_ - pi_upper_bound);
    ARIADNE_TESTING_INSPECT("constructed x", x);
    ARIADNE_TESTING_INSPECT("x.compute(Effort(12u)).get(double_precision)", x.compute(Effort(12u)).get(double_precision));
    ARIADNE_TESTING_EXPECT_UNARY_PREDICATE(is_indeterminate, (x >= 0).check(Effort(12u)));
    ARIADNE_TESTING_INSPECT("(x>=0).check(Effort(12u))", (x >= 0).check(Effort(12u)));
    ARIADNE_TESTING_EXPECT_EQUAL((x >= 0).check(Effort(20u)), false);
    //    ARIADNE_TEST_PRINT(unify(x>=0,+x,x<=0,-x));
    ARIADNE_TESTING_INSPECT("when({x>=0,+x},{x<=0,-x}).compute(Effort(12u))", when({x >= 0, +x}, {x <= 0, -x}).compute(Effort(12u)));
    ARIADNE_TESTING_INSPECT("when({x>=0,+x},{x<=0,-x}).compute(Effort(20u))", when({x >= 0, +x}, {x <= 0, -x}).compute(Effort(20u)));
    ARIADNE_TESTING_INSPECT("when({x<=0,-x},{x>=0,+x}).compute(Effort(20u))", when({x <= 0, -x}, {x >= 0, +x}).compute(Effort(20u)));
    ARIADNE_TESTING_INSPECT("when({x>=0,+x},{x<=0,-x}).get(dp)", when({x >= 0, +x}, {x <= 0, -x}).get(dp));
    ARIADNE_TESTING_INSPECT("when({x>=0,+x},{x<=0,-x}).get(MultiplePrecision(192_bits))", when({x >= 0, +x}, {x <= 0, -x}).get(MultiplePrecision(192_bits)));
    ARIADNE_TESTING_INSPECT("when({x>=0,+x},{x<=0,-x}).get(MultiplePrecision(320_bits))", when({x >= 0, +x}, {x <= 0, -x}).get(MultiplePrecision(320_bits)));
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(when({x >= 0, +x}, {x <= 0, 1 + x}).compute(Effort(20u)), 0.5_x);

    ARIADNE_TESTING_EXPECT_EXCEPTION(when({sin(pi) >= 1, +1}, {sin(pi) <= -1, -1}).compute(Effort(12u)));

    ARIADNE_TESTING_INSPECT("(x>=0).check(Effort(0u))", (x >= 0).check(Effort(0u)));
    ARIADNE_TESTING_INSPECT("(x<=0).check(Effort(0u))", (x <= 0).check(Effort(0u)));
    ARIADNE_TESTING_INSPECT("(x>=0).check(Effort(12u))", (x >= 0).check(Effort(12u)));
    ARIADNE_TESTING_INSPECT("(x<=0).check(Effort(12u))", (x <= 0).check(Effort(12u)));
    ARIADNE_TESTING_INSPECT("(x>=0).check(Effort(20u))", (x >= 0).check(Effort(20u)));
    ARIADNE_TESTING_INSPECT("(x<=0).check(Effort(20u))", (x <= 0).check(Effort(20u)));
    ARIADNE_TESTING_INSPECT("x.compute(Effort(12))", x.compute(Effort(12)));
    ARIADNE_TESTING_INSPECT("x.compute(Effort(12))<DyadicBounds(0)", x.compute(Effort(12)) < DyadicBounds(0));
    ARIADNE_TESTING_INSPECT("x.compute(Effort(20u))", x.compute(Effort(20u)));
    ARIADNE_TESTING_INSPECT("x.compute(Effort(20u)).get(MultiplePrecision(320_bits))", x.compute(Effort(20u)).get(MultiplePrecision(320_bits)));

    Real y = 2 * sin(pi / 6) - 1;
    ARIADNE_TESTING_EXPECT_UNARY_PREDICATE(is_indeterminate, (y >= 0).check(Effort(0u)));
    ARIADNE_TESTING_EXPECT_UNARY_PREDICATE(is_indeterminate, (y >= 0).check(Effort(16u)));
}

int main(int argc, char** argv) {
    std::cout << std::setprecision(20);
    std::cerr << std::setprecision(20);
    Ariadne::Testing::Suite suite("real");
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestReal, test_precision_preflight);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestReal, test_expression_evaluation);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestReal, test_constructors);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestReal, test_conversions);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestReal, test_arithmetic);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestReal, test_directed);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestReal, test_transcendental);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestReal, test_comparison);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestReal, test_rounding);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestReal, test_accuracy);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestReal, test_sequence);
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
