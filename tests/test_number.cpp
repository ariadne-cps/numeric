/***************************************************************************
 *            test_number.cpp
 *
 *  Copyright  2016-20  Alberto Casagrande, Pieter Collins
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
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "numeric/builtin.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/float_approximation.hpp"
#include "numeric/float_ball.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_error.hpp"
#include "numeric/float_lower_bound.hpp"
#include "numeric/float_upper_bound.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/lower_number.hpp"
#include "numeric/number.hpp"
#include "numeric/number_wrapper.hpp"
#include "numeric/rational.hpp"
#include "numeric/upper_number.hpp"

#include "ariadne-testing.hpp"

using namespace std;
using namespace Ariadne;

namespace {

template<class Y>
constexpr bool check_number_concept() {
    return requires(Int n, Y y, Y y2) {
        y = Y(2);
        y = 3;
        Y(1);
        y = y + y;
        y = y - y;
        y = y * y;
        y = y / y;
        y += y2;
        y -= y2;
        y *= y2;
        y /= y2;
        y = abs(y);
        y = max(y, y);
        y = min(y, y);
        y = add(y, y);
        y = sub(y, y);
        y = mul(y, y);
        y = div(y, y);
        y = sqrt(y);
        y = exp(y);
        y = log(y);
        y = atan(y);
        y = pos(y);
        y = neg(y);
        y = sqr(y);
        y = rec(y);
        y = sin(y);
        y = cos(y);
        y = tan(y);
        y == y;
        y != y;
        y <= y;
        y >= y;
        y < y;
        y > y;
        y == n;
        y != n;
        y <= n;
        y >= n;
        y < n;
        y > n;
    };
}

template<class L, class U>
constexpr bool check_directed_number_concept() {
    return requires(L l, L l2, U u, U u2) {
        +l;
        -l;
        +u;
        -u;
        l + l2;
        l - u;
        u + u2;
        u - l;
        l += l2;
        l -= u;
        u += u2;
        u -= l;
        add(l, l2);
        sub(l, u);
        add(u, u2);
        sub(u, l);
        sqrt(l);
        exp(l);
        log(l);
        atan(l);
        max(l, l2);
        min(l, l2);
        sqrt(u);
        exp(u);
        log(u);
        atan(u);
        max(u, u2);
        min(u, u2);
        l == u;
        l != u;
        l < u;
        l > u;
        u == l;
        u != l;
        u < l;
        u > l;
    };
}

static_assert(check_number_concept<ApproximateNumber>());
static_assert(check_number_concept<ValidatedNumber>());
static_assert(check_number_concept<EffectiveNumber>());
static_assert(check_number_concept<ExactNumber>());
static_assert(!std::is_abstract_v<NumberWrapper<Integer>>,
              "Every NumberWrapper must override native precision getters");
static_assert(!std::is_abstract_v<NumberWrapper<Dyadic>>);
static_assert(!std::is_abstract_v<NumberWrapper<Rational>>);
static_assert(!std::is_abstract_v<NumberWrapper<FloatDP>>);
static_assert(!std::is_abstract_v<NumberWrapper<FloatMP>>);
static_assert(std::is_convertible_v<FloatMP, ExactNumber>);
// Preserve the FloatMP generic-division call surface directly rather than
// comparing it with an independent exact arithmetic-backed type in this compilation unit.
template<class A, class B>
concept SupportsDirectDivision =
    requires(const A& a, const B& b) { div(a, b); };
static_assert(SupportsDirectDivision<Real, UpperBound<FloatMP>>);
static_assert(SupportsDirectDivision<Rational, LowerBound<FloatMP>>);
static_assert(SupportsDirectDivision<UpperBound<FloatMP>, Integer>);
static_assert(SupportsDirectDivision<LowerBound<FloatMP>, Dyadic>);
static_assert(requires(const FloatMP& a) { { fma(a,a,a) } -> Same<Bounds<FloatMP>>; });
static_assert(requires(const FloatMP& x, const FloatMP& y, const FloatDP& dp,
                       const ExactDouble& d, MultiplePrecision pr,
                       CurrentRoundingMode mode) {
    { sqr(mode, x) } -> Same<FloatMP>;
    { rec(mode, x) } -> Same<FloatMP>;
    { add(mode, x, y) } -> Same<FloatMP>;
    { sub(mode, x, y) } -> Same<FloatMP>;
    { mul(mode, x, y) } -> Same<FloatMP>;
    { div(mode, x, y) } -> Same<FloatMP>;
    { fma(mode, x, y, x) } -> Same<FloatMP>;
    { sqrt(mode, x) } -> Same<FloatMP>;
    { exp(mode, x) } -> Same<FloatMP>;
    { log(mode, x) } -> Same<FloatMP>;
    { sin(mode, x) } -> Same<FloatMP>;
    { cos(mode, x) } -> Same<FloatMP>;
    { tan(mode, x) } -> Same<FloatMP>;
    { asin(mode, x) } -> Same<FloatMP>;
    { acos(mode, x) } -> Same<FloatMP>;
    { atan(mode, x) } -> Same<FloatMP>;
    { add(mode, x, dp) } -> Same<FloatMP>;
    { add(mode, dp, x) } -> Same<FloatMP>;
    { sub(mode, x, d) } -> Same<FloatMP>;
    { sub(mode, d, x) } -> Same<FloatMP>;
    { FloatMP::pi(mode, pr) } -> Same<FloatMP>;
});

static_assert(check_directed_number_concept<ValidatedLowerNumber, ValidatedUpperNumber>());
static_assert(check_directed_number_concept<EffectiveLowerNumber, EffectiveUpperNumber>());

} // namespace

class TestNumbers : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    Void test_operations();
    Void test_dyadic_behaviour();
    Void test_float_value_behaviour();
    Void test_misc();
};

Void TestNumbers::test_dyadic_behaviour() {
    // D59: Exercise NumberInterface's native-precision virtual dispatch,
    // including both-sided outward rounding of exact rational inputs.
    const MultiplePrecision native_precision(128_bits);
    const ExactNumber exact(Dyadic(13, 2u));
    const FloatMPBounds native = exact.get(native_precision);
    ARIADNE_TESTING_PROVE(static_cast<Rational>(native.lower_raw()) <= Rational(13, 4));
    ARIADNE_TESTING_PROVE(static_cast<Rational>(native.upper_raw()) >= Rational(13, 4));
    const ValidatedNumber validated(DyadicBounds(Dyadic(1, 3u), Dyadic(3, 3u)));
    const FloatMPBounds enclosure = validated.get(native_precision);
    ARIADNE_TESTING_PROVE(static_cast<Rational>(enclosure.lower_raw()) <= Rational(1, 8));
    ARIADNE_TESTING_PROVE(static_cast<Rational>(enclosure.upper_raw()) >= Rational(3, 8));
    ExactNumber y(Dyadic(1));
    try {
        ExactNumber r = y + y;
        ARIADNE_TESTING_NOTE("Binary addition on Dyadic within ExactNumber yields " << r.class_name());
        r = y / y;
        ARIADNE_TESTING_NOTE("Binary division on Dyadic within ExactNumber yields " << r.class_name());
    } catch (const DispatchException& e) {
        ARIADNE_TESTING_NOTE("Binary operations on Dyadic give error:\n    " << e.what());
    }
}

Void TestNumbers::test_float_value_behaviour() {
    // Check rigorous enclosure and single-rounding endpoints of native FMA.
    const MultiplePrecision ap_precision(128_bits);
    const FloatMP fa(Rational(1, 3), FloatMP::ROUND_TO_NEAREST, ap_precision);
    const FloatMP fb(Rational(2, 7), FloatMP::ROUND_TO_NEAREST, ap_precision);
    const FloatMP fc(Rational(-1, 5), FloatMP::ROUND_TO_NEAREST, ap_precision);
    const Bounds<FloatMP> enclosed = fma(fa, fb, fc);
    const Rational expected = static_cast<Rational>(fa) * static_cast<Rational>(fb) + static_cast<Rational>(fc);
    ARIADNE_TESTING_PROVE(static_cast<Rational>(enclosed.lower_raw()) <= expected);
    ARIADNE_TESTING_PROVE(expected <= static_cast<Rational>(enclosed.upper_raw()));
    ARIADNE_TESTING_EXPECT_EQUAL(enclosed.lower_raw(), fma(FloatMP::ROUND_DOWNWARD, fa, fb, fc));
    ARIADNE_TESTING_EXPECT_EQUAL(enclosed.upper_raw(), fma(FloatMP::ROUND_UPWARD, fa, fb, fc));

    ExactNumber n(Int(1));
    ExactNumber z(Integer(1));
    ExactNumber w(Dyadic(1));
    ExactNumber q(Rational(1));
    ExactNumber y(FloatDP(1, dp));
    ExactNumber r;
    try {
        ARIADNE_TESTING_NOTE("Binary addition on Float<DP> within ExactNumber yields " << (y + y).class_name());
    } catch (const DispatchException& e) {
        ARIADNE_TESTING_NOTE("Binary addition on Float<DP> within ExactNumber gives error:\n    " << e.what());
    }
    try {
        ARIADNE_TESTING_NOTE("Binary division on Float<DP> within ExactNumber yields " << (y / y).class_name());
    } catch (const DispatchException& e) {
        ARIADNE_TESTING_NOTE("Binary division on Float<DP> within ExactNumber gives error:\n    " << e.what());
    }
    try {
        ARIADNE_TESTING_NOTE("Binary maximum on Float<DP> within ExactNumber yields " << max(y, y).class_name());
    } catch (const DispatchException& e) {
        ARIADNE_TESTING_NOTE("Binary maximum on Float<DP> within ExactNumber gives error:\n    " << e.what());
    }
    try {
        ARIADNE_TESTING_NOTE("Binary addition on Float<DP> and Int within ExactNumber yields " << (y + n).class_name());
    } catch (const DispatchException& e) {
        ARIADNE_TESTING_NOTE("Binary addition on Float<DP> and Int within ExactNumber gives error:\n    " << e.what());
    }
    try {
        ARIADNE_TESTING_NOTE("Binary addition on Float<DP> and Integer within ExactNumber yields " << (y + z).class_name());
    } catch (const DispatchException& e) {
        ARIADNE_TESTING_NOTE("Binary addition on Float<DP> and Integer within ExactNumber gives error:\n    " << e.what());
    }
    try {
        ARIADNE_TESTING_NOTE("Binary addition on Float<DP> and Dyadic within ExactNumber yields " << (y + w).class_name());
    } catch (const DispatchException& e) {
        ARIADNE_TESTING_NOTE("Binary addition on Float<DP> and Dyadic within ExactNumber gives error:\n    " << e.what());
    }
    try {
        ARIADNE_TESTING_NOTE("Binary addition on Float<DP> and Rational within ExactNumber yields " << (y + q).class_name());
    } catch (const DispatchException& e) {
        ARIADNE_TESTING_NOTE("Binary addition on Float<DP> and Rational within ExactNumber gives error:\n    " << e.what());
    }
    try {
        ARIADNE_TESTING_NOTE("Binary addition on Dyadic and Float<DP> within ExactNumber yields " << (w + y).class_name());
    } catch (const DispatchException& e) {
        ARIADNE_TESTING_NOTE("Binary addition on Dyadic and Float<DP> within ExactNumber gives error:\n    " << e.what());
    }
    try {
        ARIADNE_TESTING_NOTE("Binary addition on Rational and Float<DP> within ExactNumber yields " << (q + y).class_name());
    } catch (const DispatchException& e) {
        ARIADNE_TESTING_NOTE("Binary addition on Rational and Float<DP> within ExactNumber gives error:\n    " << e.what());
    }
    y = ExactNumber(FloatMP(1, MultiplePrecision(64_bits)));
    try {
        ARIADNE_TESTING_NOTE("Binary addition on Float<MP> within ExactNumber yields " << (y + y).class_name());
    } catch (const DispatchException& e) {
        ARIADNE_TESTING_NOTE("Binary addition on Float<MP> within ExactNumber gives error:\n    " << e.what());
    }
    try {
        ARIADNE_TESTING_NOTE("Binary division on Float<MP> within ExactNumber yields " << (y / y).class_name());
    } catch (const DispatchException& e) {
        ARIADNE_TESTING_NOTE("Binary division on Float<MP> within ExactNumber gives error:\n    " << e.what());
    }
    try {
        ARIADNE_TESTING_NOTE("Binary maximum on Float<MP> within ExactNumber yields " << max(y, y).class_name());
    } catch (const DispatchException& e) {
        ARIADNE_TESTING_NOTE("Binary maximum on Float<MP> within ExactNumber gives error:\n    " << e.what());
    }
}

Void TestNumbers::test_operations() {
    Nat m = 1u;
    Int n = 1;
    Integer z = 1;
    FloatDP v(3, dp);
    FloatDPBounds b(3, dp);
    ExactNumber yn = n;
    ExactNumber yz = z;
    ExactNumber yv = v;
    ValidatedNumber yb = b;
    ValidatedErrorNumber en = m;
    ValidatedErrorNumber ev = v;

    Dyadic w2(2);
    Dyadic w3(3);
    Rational q2(2);
    Rational q3(3);
    FloatDP x2(2, dp);
    FloatDP x3(3, dp);
    ExactNumber y2(x2);
    ExactNumber y3(x3);

    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, y2 / y3);
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, q2 / y3);
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, y2 / q3);

    ARIADNE_TESTING_EXPECT_EXCEPTION(add(yv, yv));
    ARIADNE_TESTING_INSPECT("add(yb,yb)", (add(yb, yb)));
    ARIADNE_TESTING_INSPECT("add(yb,yv)", (add(yb, yv)));
    ARIADNE_TESTING_INSPECT("add(b,yb)", (add(b, yb)));
    ARIADNE_TESTING_INSPECT("add(yb,yn)", (add(yb, yn)));

    ValidatedErrorNumber ym = max(en, ev);
    ARIADNE_TESTING_INSPECT("max(en,ev)", (max(en, ev)));
    ARIADNE_TESTING_INSPECT("max(en,ev).handle().pointer()", (max(en, ev).handle().pointer()));
    ARIADNE_TESTING_INSPECT("max(en,ev).class_name()", (max(en, ev).class_name()));
    ARIADNE_TESTING_INSPECT("max(en,ev)", (max(en, ev)));

    ARIADNE_TESTING_EXPECT_EXCEPTION(add(ExactNumber(1), ExactNumber(FloatDP(2, dp))));
    ARIADNE_TESTING_EXPECT_EXCEPTION(add(ExactNumber(FloatDP(1, dp)), ExactNumber(FloatDP(2, dp))));

    ARIADNE_TESTING_INSPECT("max(ExactNumber(1),ExactNumber(2))", (max(ExactNumber(1), ExactNumber(2))));
    ARIADNE_TESTING_INSPECT("max(ExactNumber(FloatDP(1,dp)),ExactNumber(FloatDP(2,dp)))", (max(ExactNumber(FloatDP(1, dp)), ExactNumber(FloatDP(2, dp)))));
    ARIADNE_TESTING_INSPECT("max(ExactNumber(1),ExactNumber(FloatDP(2,dp)))", (max(ExactNumber(1), ExactNumber(FloatDP(2, dp)))));

    ARIADNE_TESTING_INSPECT("max(ExactNumber(1),ExactNumber(FloatDP(2,dp)))", (max(ExactNumber(1), ExactNumber(FloatDP(2, dp)))));

    ARIADNE_TESTING_INSPECT("max(ExactNumber(Dyadic(1)),ExactNumber(FloatDP(2,dp)))", (max(ExactNumber(Dyadic(1)), ExactNumber(FloatDP(2, dp)))));
    ARIADNE_TESTING_INSPECT("max(ExactNumber(FloatDP(1,dp)),ExactNumber(Dyadic(2)))", (max(ExactNumber(FloatDP(1, dp)), ExactNumber(Dyadic(2)))));

    ARIADNE_TESTING_INSPECT("min(ExactNumber(1),ExactNumber(2))", (min(ExactNumber(1), ExactNumber(2))));
    ARIADNE_TESTING_INSPECT("min(ExactNumber(FloatDP(1,dp)),ExactNumber(FloatDP(2,dp)))", (min(ExactNumber(FloatDP(1, dp)), ExactNumber(FloatDP(2, dp)))));
    ARIADNE_TESTING_INSPECT("min(ExactNumber(1),ExactNumber(FloatDP(2,dp)))", (min(ExactNumber(1), ExactNumber(FloatDP(2, dp)))));
    ARIADNE_TESTING_INSPECT("min(ExactNumber(FloatDP(1,dp)),ExactNumber(2))", (min(ExactNumber(FloatDP(1, dp)), ExactNumber(2))));
    ARIADNE_TESTING_INSPECT("min(ExactNumber(Dyadic(1)),ExactNumber(FloatDP(2,dp)))", (min(ExactNumber(Dyadic(1)), ExactNumber(FloatDP(2, dp)))));
    ARIADNE_TESTING_INSPECT("min(ExactNumber(FloatDP(1,dp)),ExactNumber(Dyadic(2)))", (min(ExactNumber(FloatDP(1, dp)), ExactNumber(Dyadic(2)))));
}

Void TestNumbers::test_misc() {
    Bounds<FloatDP> x(2, 3, dp);
    ARIADNE_TESTING_INSPECT("x", x);
    //    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedNumber,y,x);
    ValidatedNumber y = x.operator ValidatedNumber();
    ARIADNE_TESTING_INSPECT("y", y);
    ValidatedLowerNumber yyl(y);
    ARIADNE_TESTING_INSPECT("yyl", yyl);

    LowerBound<FloatDP> xl(2u, dp);
    ARIADNE_TESTING_INSPECT("xl", xl);
    //    ARIADNE_TEST_ASSIGN_CONSTRUCT(ValidatedLowerNumber,yl,xl);
    ValidatedLowerNumber yl = xl.operator ValidatedLowerNumber();
    ARIADNE_TESTING_INSPECT("yl", yl);
    ARIADNE_TESTING_EXPECT_NO_THROW(yl = (xl + xl));
    ValidatedLowerNumber z = yl;
    ARIADNE_TESTING_INSPECT("z", z);
    ARIADNE_TESTING_EXPECT_NO_THROW(z = (yl + yl));

    MultiplePrecision mp(128_bits);

    DyadicBounds dyadic_bounds(Dyadic(1), Dyadic(2));
    ValidatedNumber dyadic_number = dyadic_bounds.operator ValidatedNumber();
    ARIADNE_TESTING_EXPECT_EQUAL(dyadic_number.class_name(), String("DyadicBounds"));

    FloatDPError dp_error(1u, dp);
    ValidatedErrorNumber dp_error_number = dp_error.operator ValidatedErrorNumber();
    ARIADNE_TESTING_EXPECT_EQUAL(dp_error_number.class_name(), String("FloatDP"));

    // Verify the native backend's concrete type is preserved through generic
    // validated, approximate, and error wrappers without independent exact arithmetic intermediates.
    FloatMPError ap_error(1u, mp);
    ValidatedErrorNumber ap_error_number = ap_error.operator ValidatedErrorNumber();
    ARIADNE_TESTING_EXPECT_EQUAL(ap_error_number.class_name(), String("FloatMP"));

    FloatMPApproximation ap_approximation(2u, mp);
    ApproximateNumber approximation_number = ap_approximation.operator ApproximateNumber();
    ARIADNE_TESTING_EXPECT_EQUAL(approximation_number.class_name(), String("FloatMPApproximation"));

    FloatMPBounds ap_bounds(2u, mp);
    ValidatedNumber ap_bounds_number = ap_bounds.operator ValidatedNumber();
    ARIADNE_TESTING_EXPECT_EQUAL(ap_bounds_number.class_name(), String("FloatMPBounds"));

    ValidatedUpperNumber upper_dyadic(Dyadic(2));
    ValidatedLowerNumber lower_dyadic(Dyadic(2));
    ARIADNE_TESTING_EXPECT_EQUAL(cast_exact(upper_dyadic).class_name(), String("Dyadic"));
    ARIADNE_TESTING_EXPECT_EQUAL(cast_exact(lower_dyadic).class_name(), String("Dyadic"));

    ValidatedNumber metric_number(Real(Rational(3, 5)));
    ARIADNE_TESTING_EXPECT_NO_THROW(metric_number.ref()._get(MetricTag(), dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(metric_number.ref()._get(MetricTag(), mp));

    ExactNumber integer_one(Integer(1));
    ExactNumber integer_two(Integer(2));
    ExactNumber float_one(FloatDP(1, dp));
    ExactNumber float_two(FloatDP(2, dp));

    ARIADNE_TESTING_PROVE(Ariadne::extract<Integer>(&integer_one.ref()) != nullptr);
    ARIADNE_TESTING_PROVE(Ariadne::extract<FloatDP>(&integer_one.ref()) == nullptr);
    ARIADNE_TESTING_EXPECT_EQUAL(integer_one.ref()._get_q(), Rational(1));

    ARIADNE_TESTING_EXPECT_EQUAL(cmp(float_one.ref(), float_two.ref()), Comparison::LESS);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(float_one.ref(), integer_two.ref()), Comparison::LESS);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(integer_one.ref(), float_two.ref()), Comparison::LESS);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(integer_one.ref(), integer_two.ref()), Comparison::LESS);
    ARIADNE_TESTING_EXPECT_EQUAL(integer_one.ref()._paradigm(), ParadigmCode::EXACT);

    ExactNumber mixed_dp(FloatDP(1, dp));
    ExactNumber mixed_ap(FloatMP(2, mp));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, add(mixed_dp, mixed_ap));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, (void)(mixed_dp < mixed_ap));

    ApproximateNumber unsupported_conversion(FloatDPApproximation(1u, dp));
    ARIADNE_TESTING_EXPECT_THROWS(ParadigmError, unsupported_conversion.ref()._get_q());
    ARIADNE_TESTING_EXPECT_THROWS(ParadigmError, unsupported_conversion.ref()._get(OrderTag(), dp));
    NumberWrapper<FloatDPApproximation> unsupported_wrapper(FloatDPApproximation(1u, dp));
    ARIADNE_TESTING_EXPECT_THROWS(ParadigmError, unsupported_wrapper._get(MetricTag(), dp, dp));

    std::ostringstream interface_stream;
    interface_stream << integer_one.ref();
    ARIADNE_TESTING_PROVE(not interface_stream.str().empty());

    std::ostringstream paradigm_stream;
    paradigm_stream << ParadigmCode::EXACT << ParadigmCode::VALIDATED
                    << ParadigmCode::APPROXIMATE << ParadigmCode::EFFECTIVE;
    ARIADNE_TESTING_EXPECT_EQUAL(paradigm_stream.str(), String("EXACTVALIDATEDAPPROXIMATEUNKNOWN"));

    NumberWrapper<Integer> integer_wrapper(Integer(3));
    NumberInterface* copied = integer_wrapper._copy();
    ARIADNE_TESTING_EXPECT_EQUAL(copied->_class_name(), String("Integer"));
    delete copied;
    NumberInterface* moved = integer_wrapper._move();
    ARIADNE_TESTING_EXPECT_EQUAL(moved->_class_name(), String("Integer"));
    delete moved;

    ARIADNE_TESTING_EXPECT_NO_THROW(pow(integer_one, Int(3)));
    ARIADNE_TESTING_EXPECT_NO_THROW(abs(ExactNumber(Integer(-2))));
    ARIADNE_TESTING_EXPECT_NO_THROW(abs(ExactNumber(Dyadic(-2))));
    ARIADNE_TESTING_EXPECT_NO_THROW(abs(EffectiveNumber(Real(-2))));

    ValidatedLowerNumber positive_dyadic(Dyadic(1));
    ARIADNE_TESTING_PROVE(definitely(positive_dyadic > ValidatedLowerNumber::Zero{}));
    ValidatedLowerNumber positive_float(FloatDPLowerBound(1u, dp));
    ARIADNE_TESTING_PROVE(definitely(positive_float > ValidatedLowerNumber::Zero{}));
    EffectiveLowerNumber positive_real(Real(1));
    ARIADNE_TESTING_PROVE(definitely(positive_real > EffectiveLowerNumber::Zero{}));

    // Comparison results must preserve their logical value and own the returned object.
    ValidatedLowerNumber negative_dyadic(Dyadic(-1));
    ARIADNE_TESTING_PROVE(not definitely(negative_dyadic > ValidatedLowerNumber::Zero{}));
    ARIADNE_TESTING_PROVE(possibly(negative_dyadic < ValidatedLowerNumber::Zero{}));
    ValidatedLowerNumber zero_dyadic(Dyadic(0));
    ARIADNE_TESTING_PROVE(not definitely(zero_dyadic > ValidatedLowerNumber::Zero{}));

    FloatDPApproximation dp_approximation(2u, dp);
    ARIADNE_TESTING_EXPECT_NO_THROW(pow(dp_approximation, Integer(3)));
    FloatMPApproximation ap_power_approximation(2u, mp);
    ARIADNE_TESTING_EXPECT_NO_THROW(pow(ap_power_approximation, Integer(3)));
    // The wrapper's exponent type is Int=int, which can be narrower than
    // the long returned by get_si(). Reuse the independent frozen limits.
    std::ifstream exponent_data(ARIADNE_TEST_DATA_DIR "/integer_conversion.tsv");
    if (!exponent_data)
        throw std::runtime_error("Missing integer_conversion.tsv");
    std::string exponent_line;
    std::set<std::string> exponent_sources;
    while (std::getline(exponent_data, exponent_line)) {
        if (exponent_line.empty() || exponent_line[0] == '#')
            continue;
        std::istringstream fields(exponent_line);
        std::string source, extra;
        int signed64_ok, unsigned64_ok, signed32_ok, unsigned32_ok;
        if (!(fields >> source >> signed64_ok >> unsigned64_ok >> signed32_ok >> unsigned32_ok) || (fields >> extra) || !exponent_sources.insert(source).second)
            throw std::runtime_error("Malformed exponent oracle");
        const Integer exponent{String(source)};
        const bool expected_ok = std::numeric_limits<Int>::digits == 31 ? signed32_ok != 0 : signed64_ok != 0;
        // Binary64 base one verifies every valid boundary, including INT_MIN.
        const Approximation<FloatDP> dp_one(FloatDP(1, dp));
        const Approximation<FloatMP> mp_one(FloatMP(1, MultiplePrecision(64_bits)));
        bool dp_threw = false;
        try {
            ARIADNE_TESTING_EXPECT_EQUAL(pow(dp_one, exponent).raw().get_d(), 1.0);
        } catch (const std::out_of_range&) {
            dp_threw = true;
        }
        ARIADNE_TESTING_EXPECT_EQUAL(dp_threw, !expected_ok);
        // FloatMP forms an exact significand power before final rounding;
        // validate all rejection boundaries before that work, plus small
        // valid powers. Integer's own corpus checks every valid conversion.
        if (!expected_ok || source == "0" || source == "1" || source == "-1") {
            bool mp_threw = false;
            try {
                ARIADNE_TESTING_EXPECT_EQUAL(pow(mp_one, exponent).raw().get_d(), 1.0);
            } catch (const std::out_of_range&) {
                mp_threw = true;
            }
            ARIADNE_TESTING_EXPECT_EQUAL(mp_threw, !expected_ok);
        }
    }
    if (exponent_sources.size() != 27u)
        throw std::runtime_error("Incomplete exponent oracle");

    Integer huge_power{String("9223372036854775808")};
    ARIADNE_TESTING_EXPECT_EXCEPTION(pow(dp_approximation, huge_power));
    ARIADNE_TESTING_EXPECT_EXCEPTION(pow(ap_power_approximation, huge_power));

    ExactNumber unary(Integer(1));
    ARIADNE_TESTING_EXPECT_NO_THROW(nul(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(pos(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(neg(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(sqr(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(hlf(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(rec(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(sqrt(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(exp(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(log(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(sin(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(cos(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(tan(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(tanh(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(asin(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(acos(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(atan(unary));
    ARIADNE_TESTING_EXPECT_NO_THROW(abs(unary));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, unary.ref()._apply(UnaryElementaryOperator(OperatorCode::ADD)));

    DyadicBounds algebraic_dyadic(Dyadic(1), Dyadic(2));
    RationalBounds algebraic_rational(Rational(1), Rational(2));
    ARIADNE_TESTING_PROVE(definitely(algebraic_dyadic > 0));
    ARIADNE_TESTING_PROVE(definitely(algebraic_rational > 0));
    ValidatedNumber dyadic_dispatch = algebraic_dyadic.operator ValidatedNumber();
    ValidatedNumber rational_dispatch = algebraic_rational.operator ValidatedNumber();
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, sqrt(dyadic_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, exp(dyadic_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, log(dyadic_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, sin(dyadic_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, cos(dyadic_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, tan(dyadic_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, tanh(dyadic_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, asin(dyadic_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, acos(dyadic_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, atan(dyadic_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, sqrt(rational_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, exp(rational_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, log(rational_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, sin(rational_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, cos(rational_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, tan(rational_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, tanh(rational_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, asin(rational_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, acos(rational_dispatch));
    ARIADNE_TESTING_EXPECT_THROWS(DispatchException, atan(rational_dispatch));
}

template<class Y>
class TestNumber : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    static void add_cases(Ariadne::Testing::Suite& suite, const std::string& variant);

  private:
    Void test_number_get();
    // A Float approximation's raw value is an exact binary rational. Check
    // the original tolerance without an undecidable approximate comparison.
    template<class F>
    Void expect_exact_approximation_within(
        const Approximation<F>& actual, const Rational& expected,
        const Dyadic& tolerance, const char* expression) {
        auto diagnostic_scope = this->ctx().scope("approximation", expression);
        const Rational exact_value{Dyadic(actual.raw())};
        const Rational exact_error = abs(exact_value - expected);
        ARIADNE_TESTING_EXPECT_LESS_EQUAL(exact_error, Rational(tolerance));
    }
    // Local predicates for the Number fixture's numeric model checks.
    template<class F, class FE>
    Bool models(const Ball<F, FE>& x, const Rational& q) {
        return abs(Dyadic(x.value().raw()) - q) <= Dyadic(x.error().raw());
    }
    template<class F>
    Bool models(const Bounds<F>& x, const Rational& q) {
        return Dyadic(x.lower().raw()) <= q && q <= Dyadic(x.upper().raw());
    }
    template<class F>
    Bool models(const LowerBound<F>& x, const Rational& q) {
        return Dyadic(x.raw()) <= q;
    }
    template<class F>
    Bool models(const UpperBound<F>& x, const Rational& q) {
        return q <= Dyadic(x.raw());
    }

    Void test_class();
    Void test_get();
    Void test_operations();
    Void test_comparisons();
};

template<class Y>
void TestNumber<Y>::add_cases(Ariadne::Testing::Suite& suite, const std::string& variant) {
    using Fixture = TestNumber<Y>;
    auto add = [&](const std::string& name, auto method) {
        suite.add(variant + " / " + name, [method](Ariadne::Testing::Context& context) {
            Fixture fixture(context);
            (fixture.*method)();
        });
    };
    add("class", &Fixture::test_class);
    add("get", &Fixture::test_get);
    add("comparisons", &Fixture::test_comparisons);
}

template<class Y>
void TestNumber<Y>::test_number_get() {
    Rational q = 3 / 5_q;
    Real r(q);
    Y y(r);
    MultiplePrecision mp(192_bits);
    MultiplePrecision ape(128_bits);

    expect_exact_approximation_within(
        ApproximateNumber(y).get(dp), q, Dyadic(TwoExp(-52)),
        "ApproximateNumber(y).get(dp)");
    expect_exact_approximation_within(
        ApproximateNumber(y).get(mp), q, Dyadic(TwoExp(-128)),
        "ApproximateNumber(y).get(mp)");

    ARIADNE_TESTING_PROVE(models(ValidatedLowerNumber(y).get(dp), q));
    ARIADNE_TESTING_PROVE(models(ValidatedLowerNumber(y).get(mp), q));

    ARIADNE_TESTING_PROVE(models(ValidatedUpperNumber(y).get(dp), q));
    ARIADNE_TESTING_PROVE(models(ValidatedUpperNumber(y).get(mp), q));

    ARIADNE_TESTING_PROVE(models(y.get(dp), q));
    ARIADNE_TESTING_PROVE(models(y.get(mp), q));

    ARIADNE_TESTING_PROVE(models(y.get(dp, dp), q));
    ARIADNE_TESTING_PROVE(models(y.get(mp, dp), q));
    ARIADNE_TESTING_PROVE(models(y.get(mp, ape), q));

    ARIADNE_TESTING_COMPILE_SAME_TYPE(decltype(y.get(dp, dp)), FloatDPBall);
    ARIADNE_TESTING_COMPILE_SAME_TYPE(decltype(y.get(mp, dp)), FloatMPDPBall);
    ARIADNE_TESTING_COMPILE_SAME_TYPE(decltype(y.get(mp, ape)), FloatMPBall);
}

template<>
void TestNumber<ExactNumber>::test_number_get() {}

template<>
void TestNumber<ApproximateNumber>::test_number_get() {
    Rational q = 3 / 5_q;
    ApproximateNumber y(q);
    MultiplePrecision mp(128_bits);
    // Dyadic(3e-16) preserves the exact binary64 value of the old threshold.
    expect_exact_approximation_within(y.get(dp), q, Dyadic(3e-16), "y.get(dp)");
    expect_exact_approximation_within(y.get(mp), q, Dyadic(3e-16), "y.get(mp)");

    // Regression tests for native FloatDP/FloatMP conversion.
    Dbl d = -0.6;
    ARIADNE_TESTING_EXPECT_NO_THROW(y = (d));
    ARIADNE_TESTING_EXPECT_NO_THROW(y.get(dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(y.get(mp));

    FloatDPApproximation x(-2.5, double_precision);
    ARIADNE_TESTING_EXPECT_NO_THROW(y = (x));
    ARIADNE_TESTING_EXPECT_NO_THROW(y.get(dp));
    ARIADNE_TESTING_EXPECT_NO_THROW(y.get(mp));
}

template<class Y>
Void TestNumber<Y>::test_get() {
    test_number_get();
}

template<class Y>
Void TestNumber<Y>::test_class() {
    Int n = 1;
    Integer z = 1;
    FloatDP v(3, dp);
    FloatDPBounds b(3, dp);
    ExactNumber yn = n;
    ExactNumber yz = z;
    ExactNumber yv = v;
    ValidatedNumber yb = b;

    ARIADNE_TESTING_EXPECT_EQUAL(yz.class_name(), "Integer");
    ARIADNE_TESTING_EXPECT_EQUAL(yv.class_name(), "FloatDP");
    ARIADNE_TESTING_EXPECT_EQUAL(yb.class_name(), "FloatDPBounds");
}

template<class Y>
Void TestNumber<Y>::test_comparisons() {
}

template<>
Void TestNumber<ExactNumber>::test_comparisons() {
    ExactNumber y1(Rational(2, 3));
    ARIADNE_TESTING_INSPECT("y1", y1);
    ExactNumber y2(Rational(683, 1024));
    ARIADNE_TESTING_INSPECT("y2", y2);
    ExactNumber pinf(+ExactDouble::inf());
    ARIADNE_TESTING_INSPECT("pinf", pinf);
    ExactNumber ninf(-ExactDouble::inf());
    ARIADNE_TESTING_INSPECT("ninf", ninf);

    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator==, y1, y1);
    ARIADNE_TESTING_EXPECT_NOT_BINARY_PREDICATE(operator!=, y1, y1);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<=, y1, y1);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator>=, y1, y1);
    ARIADNE_TESTING_EXPECT_NOT_BINARY_PREDICATE(operator<, y1, y1);
    ARIADNE_TESTING_EXPECT_NOT_BINARY_PREDICATE(operator>, y1, y1);

    ARIADNE_TESTING_EXPECT_NOT_BINARY_PREDICATE(operator==, y1, y2);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator!=, y1, y2);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<=, y1, y2);
    ARIADNE_TESTING_EXPECT_NOT_BINARY_PREDICATE(operator>=, y1, y2);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<, y1, y2);
    ARIADNE_TESTING_EXPECT_NOT_BINARY_PREDICATE(operator>, y1, y2);

    ARIADNE_TESTING_EXPECT_NOT_BINARY_PREDICATE(operator==, y2, y1);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator!=, y2, y1);
    ARIADNE_TESTING_EXPECT_NOT_BINARY_PREDICATE(operator<=, y2, y1);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator>=, y2, y1);
    ARIADNE_TESTING_EXPECT_NOT_BINARY_PREDICATE(operator<, y2, y1);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator>, y2, y1);

    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator==, y2, y2);
    ARIADNE_TESTING_EXPECT_NOT_BINARY_PREDICATE(operator<, y2, y2);

    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator==, y1, y1);

    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator==, pinf, pinf);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator!=, ninf, pinf);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<, ninf, pinf);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator!=, y1, pinf);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator!=, pinf, y1);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<, y1, pinf);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<, ninf, y1);

    ARIADNE_TESTING_PROVE(definitely(ValidatedUpperNumber(3) < ValidatedLowerNumber(4)));
    ARIADNE_TESTING_PROVE(possibly(ValidatedUpperNumber(4) > ValidatedLowerNumber(3)));
    ARIADNE_TESTING_PROVE(not definitely(ValidatedUpperNumber(4) > ValidatedLowerNumber(3)));

    Dyadic w1(1.25_x);
    ARIADNE_TESTING_INSPECT("w1", w1);
    Dyadic w2(0.875_x);
    ARIADNE_TESTING_INSPECT("w2", w2);
    FloatDP d1(w1, dp);
    ARIADNE_TESTING_INSPECT("d1", d1);
    FloatDP d2(w2, dp);
    ARIADNE_TESTING_INSPECT("d2", d2);
    ExactNumber yw1(w1);
    ARIADNE_TESTING_INSPECT("yw1", yw1);
    ExactNumber yw2(w2);
    ARIADNE_TESTING_INSPECT("yw2", yw2);
    ExactNumber yd1(d1);
    ARIADNE_TESTING_INSPECT("yd1", yd1);
    ExactNumber yd2(d2);
    ARIADNE_TESTING_INSPECT("yd2", yd2);

    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator==, yw1, yd1);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator!=, yw1, yd2);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator>, yw1, yd2);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<, yw2, yd1);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator>=, yw1, yd2);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<=, yw2, yd1);

    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator==, yd1, yd1);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator!=, yd1, yd2);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator>, yd1, yd2);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<, yd2, yd1);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator>=, yd1, yd2);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(operator<=, yd2, yd1);
}

template<>
Void TestNumber<EffectiveNumber>::test_comparisons() {
    ARIADNE_TESTING_PROVE(definitely(EffectiveNumber(Real(-inf)) < EffectiveNumber(Real(0))));
    ARIADNE_TESTING_PROVE(definitely(EffectiveNumber(Real(0)) > EffectiveNumber(Real(-inf))));
    ARIADNE_TESTING_PROVE(not possibly(EffectiveNumber(Real(-inf)) > EffectiveNumber(Real(0))));
    ARIADNE_TESTING_PROVE(not possibly(EffectiveNumber(Real(0)) < EffectiveNumber(Real(-inf))));
    ARIADNE_TESTING_PROVE(definitely(not(EffectiveNumber(Real(-inf)) > EffectiveNumber(Real(0)))));
    ARIADNE_TESTING_PROVE(definitely(not(EffectiveNumber(Real(0)) < EffectiveNumber(Real(-inf)))));
    ARIADNE_TESTING_PROVE(definitely(EffectiveNumber(Real(-inf)) <= EffectiveNumber(Real(0))));

    Real r1(1 / 3_q);
    ARIADNE_TESTING_INSPECT("r1", r1);
    Real r2(ExactDouble(0.333333333333333333));
    ARIADNE_TESTING_INSPECT("r2", r2);
    ARIADNE_TESTING_PROVE(not definitely((r1 > r2).check(Effort(0))));
    ARIADNE_TESTING_PROVE(possibly(not(r1 > r2).check(Effort(0))));
    ARIADNE_TESTING_PROVE(definitely((r1 > r2).check(Effort(64))));

    EffectiveNumber y1(r1);
    ARIADNE_TESTING_INSPECT("y1", y1);
    EffectiveNumber y2(r2);
    ARIADNE_TESTING_INSPECT("y2", y2);
    ARIADNE_TESTING_PROVE(not definitely((y1 > y2).check(Effort(0))));
    ARIADNE_TESTING_PROVE(possibly(not(y1 > y2).check(Effort(0))));
    ARIADNE_TESTING_PROVE(definitely((y1 > y2).check(Effort(64))));
}

template<class Y>
class TestDirectedNumber : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    static void add_cases(Ariadne::Testing::Suite& suite, const std::string& variant);

  private:
    Void test_concept();
    Void test_operations();
};

template<class Y>
void TestDirectedNumber<Y>::add_cases(Ariadne::Testing::Suite& suite, const std::string& variant) {
    suite.add(variant + " / operations", [](Ariadne::Testing::Context& context) {
        TestDirectedNumber<Y> fixture(context);
        fixture.test_operations();
    });
}

template<class Y>
Void TestDirectedNumber<Y>::test_operations() {
    if constexpr (Same<Y, ValidatedLowerNumber>) {
        typedef DoublePrecision PR;
        PR pr;
        Rational q(1, 3);
        FloatBounds<PR> x(q, pr);
        FloatLowerBound<PR> xl = x;
        FloatUpperBound<PR> xu = x * 7 / 8;

        ValidatedNumber y(x);
        ValidatedLowerNumber yl(xl);
        ValidatedUpperNumber yu(xu);

        // Tests
        ARIADNE_TESTING_PROVE(not possibly(yu == yl));
        ARIADNE_TESTING_PROVE(definitely(yu != yl));
        ARIADNE_TESTING_PROVE(not definitely(yl > q));
        ARIADNE_TESTING_PROVE(not definitely(-yl < -q));
        ARIADNE_TESTING_PROVE(not definitely(yl + yl + yl > 1));
        ARIADNE_TESTING_PROVE(not definitely(yl + yl + yl > 1));
        ARIADNE_TESTING_PROVE(definitely(yu < yl));
        ARIADNE_TESTING_PROVE(definitely(yl > yu));

        ARIADNE_TESTING_PROVE(not definitely(yl > q));

        static_assert(Same<DifferenceType<ExactNumber, ValidatedUpperNumber>, ValidatedLowerNumber>);
        static_assert(Same<DifferenceType<ValidatedLowerNumber, ValidatedUpperNumber>, ValidatedLowerNumber>);
        static_assert(Same<DifferenceType<ValidatedLowerNumber, ValidatedNumber>, ValidatedLowerNumber>);

        static_assert(Same<DifferenceType<Integer, FloatDPUpperBound>, FloatDPLowerBound>);
        static_assert(Same<DifferenceType<Dyadic, FloatDPUpperBound>, FloatDPLowerBound>);
        static_assert(Same<DifferenceType<FloatDPUpperBound, Integer>, FloatDPUpperBound>);
        static_assert(Same<DifferenceType<FloatDPUpperBound, Dyadic>, FloatDPUpperBound>);

        static_assert(Same<decltype(Add().operator()(xu, xu)), FloatDPUpperBound>);
        static_assert(Same<decltype(Sub().operator()(xu, xl)), FloatDPUpperBound>);
        static_assert(Same<decltype(BinaryElementaryOperator(Add()).operator()(xu, xu)), FloatDPApproximation>);
        static_assert(Same<decltype(BinaryElementaryOperator(Sub()).operator()(xu, xl)), FloatDPApproximation>);

        ARIADNE_TESTING_NOTE("BinaryElementaryOperator returns common arithmetical type, which may differ from type of actual operation.");

        Integer z(5);
        Dyadic w(1.375_x);
        ARIADNE_TESTING_INSPECT("w<yu", (w < yu));
        ARIADNE_TESTING_INSPECT("ExactNumber(w).class_name()", (ExactNumber(w).class_name()));
        ARIADNE_TESTING_INSPECT("ValidatedLowerNumber(w).class_name()", (ValidatedLowerNumber(w).class_name()));

        ARIADNE_TESTING_EXPECT_EQUAL(ExactNumber(w).class_name(), "Dyadic");
        ARIADNE_TESTING_EXPECT_EQUAL((yl + yl).class_name(), "FloatDPLowerBound");
        ARIADNE_TESTING_EXPECT_EQUAL((yl - yu).class_name(), "FloatDPLowerBound");
        ARIADNE_TESTING_EXPECT_EQUAL((yu + yu).class_name(), "FloatDPUpperBound");
        ARIADNE_TESTING_EXPECT_EQUAL((yu - yl).class_name(), "FloatDPUpperBound");
        ARIADNE_TESTING_EXPECT_EQUAL((ExactNumber(z) - yu).class_name(), "FloatDPLowerBound");

        ARIADNE_TESTING_PROVE(not possibly(w < yu));
    }
}

int main(int argc, char** argv) {
    std::cout << std::setprecision(20);
    std::cerr << std::setprecision(20);
    Ariadne::Testing::Suite suite("number");
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestNumbers, test_operations);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestNumbers, test_dyadic_behaviour);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestNumbers, test_float_value_behaviour);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestNumbers, test_misc);
    TestNumber<ApproximateNumber>::add_cases(suite, "ApproximateNumber");
    TestNumber<ValidatedNumber>::add_cases(suite, "ValidatedNumber");
    TestNumber<EffectiveNumber>::add_cases(suite, "EffectiveNumber");
    TestNumber<ExactNumber>::add_cases(suite, "ExactNumber");
    TestDirectedNumber<ValidatedLowerNumber>::add_cases(suite, "ValidatedLowerNumber");
    TestDirectedNumber<ValidatedUpperNumber>::add_cases(suite, "ValidatedUpperNumber");
    TestDirectedNumber<EffectiveLowerNumber>::add_cases(suite, "EffectiveLowerNumber");
    TestDirectedNumber<EffectiveUpperNumber>::add_cases(suite, "EffectiveUpperNumber");
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
