/***************************************************************************
 *            test_float_ball.cpp
 *
 *  Copyright  2006-20  Alberto Casagrande, Pieter Collins
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "numeric/builtin.hpp"
#include "numeric/decimal.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/float_ball.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_error.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "numeric/real.hpp"
#include "test_floats.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;
using namespace std;

template<class PR, class PRE = PR>
class TestFloatBall
    : public TestFloats<PR>,
      public Ariadne::Testing::ContextFixture {
    typedef RawFloat<PR> RawFloatType;
    typedef FloatType<ApproximateTag, PR> FloatApproximationType;
    typedef FloatType<LowerTag, PR> FloatLowerBoundType;
    typedef FloatType<UpperTag, PR> FloatUpperBoundType;
    typedef FloatType<BoundedTag, PR> FloatBoundsType;
    typedef FloatBall<PR, PRE> FloatBallType;
    typedef FloatType<ExactTag, PR> FloatValueType;

  private:
    PR precision;
    PRE error_precision;
    static PRE _make_error_precision(PR const& pr) {
        if constexpr (Same<PR, PRE>) {
            return pr;
        } else {
            return PRE();
        }
    }

  public:
    TestFloatBall(Ariadne::Testing::Context& context, PR prec)
        : ContextFixture(context), precision(prec), error_precision(_make_error_precision(prec)) {};
    TestFloatBall(Ariadne::Testing::Context& context, PR prec, PRE err_prec)
        : ContextFixture(context), precision(prec), error_precision(err_prec) {};

  public:
    using TestFloats<PR>::to_rational;
    Void test_header_api();
    Void test_precision();
    Void test_conversions();
    Void test_rounded_arithmetic();
    static void add_cases(Ariadne::Testing::Suite& suite, const std::string& variant, PR precision);
};

template<class PR, class PRE>
Void TestFloatBall<PR, PRE>::test_header_api() {
    RawFloatType raw_one(1, precision);
    RawFloatType raw_two(2, precision);
    RawFloatType raw_half(Dyadic(1, 1u), precision);
    RawFloat<PRE> raw_error_one(1, error_precision);

    FloatBallType by_precision(precision);
    FloatBallType by_precisions(precision, error_precision);
    typename FloatBallType::CharacteristicsType characteristics_pair(precision, error_precision);
    FloatBallType by_characteristics(characteristics_pair);
    FloatBallType by_raw(raw_one);

    FloatBallType from_twoexp(TwoExp(1), precision);
    FloatBallType from_integer(Integer(2), precision);
    FloatBallType from_dyadic(Dyadic(3, 1u), precision);
    FloatBallType from_decimal(Decimal(String("1.5")), precision);
    Real real_value(Rational(3, 2));
    FloatBallType from_real(real_value, precision);

    ARIADNE_TESTING_EXPECT_EQUAL(from_twoexp.value_raw(), RawFloatType(2, precision));
    ARIADNE_TESTING_EXPECT_EQUAL(from_twoexp.error_raw(), RawFloat<PRE>(0, error_precision));
    ARIADNE_TESTING_EXPECT_EQUAL(from_integer.value_raw(), RawFloatType(2, precision));
    ARIADNE_TESTING_EXPECT_EQUAL(from_dyadic.value_raw(), RawFloatType(Dyadic(3, 1u), precision));
    ARIADNE_TESTING_EXPECT_EQUAL(from_decimal.value_raw(), RawFloatType(Dyadic(3, 1u), precision));
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, from_real, Rational(3, 2));

    ARIADNE_TESTING_EXPECT_EQUAL(by_precision.value_raw(), RawFloatType(0, precision));
    ARIADNE_TESTING_EXPECT_EQUAL(by_precisions.error_raw(), RawFloat<PRE>(0, error_precision));
    ARIADNE_TESTING_EXPECT_EQUAL(by_characteristics.value_raw(), RawFloatType(0, precision));
    ARIADNE_TESTING_EXPECT_EQUAL(by_raw.value_raw(), raw_one);

    ValidatedNumber validated = by_raw.operator ValidatedNumber();
    FloatBallType from_validated(validated, precision);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, from_validated, Rational(1));

    FloatBallType assigned(precision, error_precision);
    assigned = validated;
    ARIADNE_TESTING_EXPECT_NO_THROW(by_raw.create(validated));
    ARIADNE_TESTING_EXPECT_NO_THROW(by_raw.generic());

    ARIADNE_TESTING_EXPECT_EQUAL(by_raw.get_d(), 1.0);

    FloatBallType interval(raw_two, raw_error_one);
    ARIADNE_TESTING_EXPECT_EQUAL(interval.lower().raw(), RawFloatType(1, precision));
    ARIADNE_TESTING_EXPECT_EQUAL(interval.upper().raw(), RawFloatType(3, precision));

    FloatBallType truncated = trunc(from_integer);
    FloatBallType truncated_digits = trunc(from_integer, 8u);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, truncated, Integer(2));
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, truncated_digits, Integer(2));

    StringStream input;
    input << "1" << "\xC2\xB1" << "0";
    FloatBallType parsed(precision, error_precision);
    input >> parsed;
    ARIADNE_TESTING_EXPECT_EQUAL(parsed.value_raw(), RawFloatType(1, precision));
    ARIADNE_TESTING_EXPECT_EQUAL(parsed.error_raw(), RawFloat<PRE>(0, error_precision));

    StringStream incomplete_input;
    incomplete_input << "1" << "\xC2\xB1";
    ARIADNE_TESTING_EXPECT_EXCEPTION(incomplete_input >> parsed);

    StringStream invalid_separator("1xx0");
    ARIADNE_TESTING_EXPECT_EXCEPTION(invalid_separator >> parsed);

    if constexpr (Same<PR, MultiplePrecision>) {
        MultiplePrecision lower_precision(53_bits);

        FloatBallType source_down(Rational(1, 3), precision);
        FloatBallType narrowed_down(source_down, lower_precision);
        FloatBallType narrowed_explicit_down(source_down, lower_precision, error_precision);
        ARIADNE_TESTING_PROVE(narrowed_down.value_raw() < source_down.value_raw());
        ARIADNE_TESTING_PROVE(narrowed_explicit_down.value_raw() < source_down.value_raw());
        ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, narrowed_down, Rational(1, 3));
        ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, narrowed_explicit_down, Rational(1, 3));

        FloatBallType source_up(Rational(1, 5), precision);
        FloatBallType narrowed_up(source_up, lower_precision);
        FloatBallType narrowed_explicit_up(source_up, lower_precision, error_precision);
        ARIADNE_TESTING_PROVE(narrowed_up.value_raw() > source_up.value_raw());
        ARIADNE_TESTING_PROVE(narrowed_explicit_up.value_raw() > source_up.value_raw());
        ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, narrowed_up, Rational(1, 5));
        ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, narrowed_explicit_up, Rational(1, 5));
    }

    auto characteristics = by_raw.characteristics();
    ARIADNE_TESTING_EXPECT_EQUAL(std::get<0>(characteristics), precision);
    ARIADNE_TESTING_EXPECT_EQUAL(std::get<1>(characteristics), error_precision);

    ARIADNE_TESTING_EXPECT_NO_THROW(nul(by_raw));
    ARIADNE_TESTING_EXPECT_NO_THROW(pos(by_raw));
    ARIADNE_TESTING_EXPECT_NO_THROW(neg(by_raw));
    ARIADNE_TESTING_EXPECT_NO_THROW(tanh(by_raw));

    FloatBallType trigonometric(raw_half);
    ARIADNE_TESTING_EXPECT_NO_THROW(asin(trigonometric));
    ARIADNE_TESTING_EXPECT_NO_THROW(acos(trigonometric));

    FloatBallType other(raw_two);

    FloatBoundsType bounds_one(1, precision);
    FloatBallType ball_two(raw_two);
    ARIADNE_TESTING_EXPECT_SAME(add(bounds_one, ball_two), FloatBoundsType(3, precision));
    ARIADNE_TESTING_EXPECT_SAME(sub(bounds_one, ball_two), FloatBoundsType(-1, precision));
    ARIADNE_TESTING_EXPECT_SAME(mul(bounds_one, ball_two), FloatBoundsType(2, precision));
    ARIADNE_TESTING_EXPECT_SAME(div(bounds_one, ball_two), FloatBoundsType(Dyadic(1, 1u), precision));
    ARIADNE_TESTING_EXPECT_SAME(add(ball_two, bounds_one), FloatBoundsType(3, precision));
    ARIADNE_TESTING_EXPECT_SAME(sub(ball_two, bounds_one), FloatBoundsType(1, precision));
    ARIADNE_TESTING_EXPECT_SAME(mul(ball_two, bounds_one), FloatBoundsType(2, precision));
    ARIADNE_TESTING_EXPECT_SAME(div(ball_two, bounds_one), FloatBoundsType(2, precision));

    ARIADNE_TESTING_EXPECT_NO_THROW(eq(by_raw, other));
    ARIADNE_TESTING_EXPECT_NO_THROW(lt(by_raw, other));

    ARIADNE_TESTING_EXPECT_NO_THROW(add(by_raw, raw_two));
    ARIADNE_TESTING_EXPECT_NO_THROW(sub(by_raw, raw_two));
    ARIADNE_TESTING_EXPECT_NO_THROW(mul(by_raw, raw_two));
    ARIADNE_TESTING_EXPECT_NO_THROW(div(by_raw, raw_two));
    ARIADNE_TESTING_EXPECT_NO_THROW(add(raw_two, by_raw));
    ARIADNE_TESTING_EXPECT_NO_THROW(sub(raw_two, by_raw));
    ARIADNE_TESTING_EXPECT_NO_THROW(mul(raw_two, by_raw));
    ARIADNE_TESTING_EXPECT_NO_THROW(div(raw_two, by_raw));

    ARIADNE_TESTING_EXPECT_NO_THROW(max(by_raw, raw_two));
    ARIADNE_TESTING_EXPECT_NO_THROW(min(by_raw, raw_two));
    ARIADNE_TESTING_EXPECT_NO_THROW(max(raw_two, by_raw));
    ARIADNE_TESTING_EXPECT_NO_THROW(min(raw_two, by_raw));

    ARIADNE_TESTING_EXPECT_NO_THROW((void)(by_raw == raw_one));
    ARIADNE_TESTING_EXPECT_NO_THROW((void)(by_raw != raw_two));
    ARIADNE_TESTING_EXPECT_NO_THROW((void)(by_raw < raw_two));
    ARIADNE_TESTING_EXPECT_NO_THROW((void)(by_raw > RawFloatType(0, precision)));
    ARIADNE_TESTING_EXPECT_NO_THROW((void)(by_raw <= raw_one));
    ARIADNE_TESTING_EXPECT_NO_THROW((void)(by_raw >= raw_one));

    ARIADNE_TESTING_EXPECT_NO_THROW((void)(raw_one == by_raw));
    ARIADNE_TESTING_EXPECT_NO_THROW((void)(raw_two != by_raw));
    ARIADNE_TESTING_EXPECT_NO_THROW((void)(RawFloatType(0, precision) < by_raw));
    ARIADNE_TESTING_EXPECT_NO_THROW((void)(raw_two > by_raw));
    ARIADNE_TESTING_EXPECT_NO_THROW((void)(raw_one <= by_raw));
    ARIADNE_TESTING_EXPECT_NO_THROW((void)(raw_one >= by_raw));

    auto fac = factory(by_raw);
    ARIADNE_TESTING_EXPECT_NO_THROW(fac.create(Real(Rational(1, 2))));
    ARIADNE_TESTING_EXPECT_NO_THROW(fac.create(Rational(1, 2)));
    ARIADNE_TESTING_EXPECT_NO_THROW(fac.create(Dyadic(1, 1u)));
    ARIADNE_TESTING_EXPECT_NO_THROW(fac.create(Integer(1)));

    FloatBallType wide(RawFloatType(0, precision), raw_error_one);
    auto squared = sqr(wide);
    ARIADNE_TESTING_PROVE(squared.error_raw() >= RawFloat<PRE>(0, error_precision));
}

template<class PR, class PRE>
Void TestFloatBall<PR, PRE>::test_conversions() {
    Rational one = 1;
    Rational third = one / 3;

    ARIADNE_TESTING_EXPECT_EQUAL(cast_integer(FloatBall<PR, PRE>(FloatBounds<PR>(2, 3, precision))), Integer(3));
    ARIADNE_TESTING_EXPECT_EQUAL(cast_integer(FloatBall<PR, PRE>(FloatBounds<PR>(2.25_dy, 3.25_dy, precision))), Integer(3));
    ARIADNE_TESTING_EXPECT_EXCEPTION(cast_integer(FloatBall<PR, PRE>(FloatBounds<PR>(2.625_dy, 2.875_dy, precision))));

    RawFloatType value(2, precision);
    FloatError<PRE> error(1u, error_precision);
    FloatBallType centred = value.pm(error);
    ARIADNE_TESTING_EXPECT_EQUAL(centred.value_raw(), value);
    ARIADNE_TESTING_EXPECT_EQUAL(centred.error_raw(), RawFloat<PRE>(1, error_precision));

    FloatBallType widened = centred.pm(error);
    ARIADNE_TESTING_EXPECT_EQUAL(widened.value_raw(), value);
    ARIADNE_TESTING_EXPECT_EQUAL(widened.error_raw(), RawFloat<PRE>(2, error_precision));

    if constexpr (Same<PR, DoublePrecision>) {
        ARIADNE_TESTING_EXPECT_EQUAL(class_name<FloatBallType>(), String("FloatDPBall"));
    } else if constexpr (Same<PRE, DoublePrecision>) {
        ARIADNE_TESTING_EXPECT_EQUAL(class_name<FloatBallType>(), String("FloatMPDPBall"));
    } else {
        ARIADNE_TESTING_EXPECT_EQUAL(class_name<FloatBallType>(), String("FloatMPBall"));
    }

    if constexpr (Same<PR, MultiplePrecision> and Same<PRE, MultiplePrecision>) {
        FloatBallType zero_with_unit_error(RawFloatType(0, precision), RawFloat<PRE>(1, error_precision));
        StringStream stream;
        stream << zero_with_unit_error;
        ARIADNE_TESTING_PROVE(not stream.str().empty());
    }
}

template<class PR, class PRE>
Void TestFloatBall<PR, PRE>::test_precision() {
    FloatBallType mx(Rational(1), precision);
    ARIADNE_TESTING_EXPECT_EQUAL(mx.value().precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL(mx.error().precision(), error_precision);
    ARIADNE_TESTING_EXPECT_EQUAL(max(mx, mx).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL(min(mx, mx).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL(abs(mx).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL((mx + mx).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL((mx - mx).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL((mx * mx).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL((mx / mx).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL(sqrt(mx).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL(exp(mx).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL(log(mx).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL(sin(mx).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL(cos(mx).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL(tan(mx).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL(atan(mx).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL((mx + 2u).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL((mx - 2u).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL((mx * 2u).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL((mx / 2u).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL((mx + 2).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL((mx - 2).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL((mx * 2).precision(), precision);
    ARIADNE_TESTING_EXPECT_EQUAL((mx / 2).precision(), precision);
}

template<class PR, class PRE>
Void TestFloatBall<PR, PRE>::test_rounded_arithmetic() {
    Rational one = 1;
    Rational three = 3;
    Rational five = 5;
    Rational six = 6;
    Rational third = one / three;
    Rational fifth = one / five;

    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, FloatBallType(third, precision) + FloatBallType(fifth, precision), third + fifth);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, FloatBallType(third, precision) - FloatBallType(fifth, precision), third - fifth);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, FloatBallType(third, precision) * FloatBallType(fifth, precision), third * fifth);
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, FloatBallType(third, precision) / FloatBallType(fifth, precision), third / fifth);
    //    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models,1/FloatBallType(three,precision)/FloatBallType(fifth,precision),1/three);

    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, pow(FloatBallType(three / five, precision), 4u), pow(three / five, 4u));
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, pow(FloatBallType(three / five, precision), 4), pow(three / five, 4));
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, pow(FloatBallType(three / five, precision), -4), pow(three / five, -4));
    ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(models, pow(FloatBallType(three / five, precision), -7), pow(three / five, -7));
}

template<class PR, class PRE>
void TestFloatBall<PR, PRE>::add_cases(Ariadne::Testing::Suite& suite, const std::string& variant, PR precision) {
    using Fixture = TestFloatBall<PR, PRE>;
    auto add = [&](const std::string& name, auto method) {
        suite.add(variant + " / " + name, [precision, method](Ariadne::Testing::Context& context) {
            Fixture fixture(context, precision);
            (fixture.*method)();
        });
    };
    add("header API", &Fixture::test_header_api);
    add("precision", &Fixture::test_precision);
    add("conversions", &Fixture::test_conversions);
    add("rounded arithmetic", &Fixture::test_rounded_arithmetic);
}

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("float ball");
    TestFloatBall<DoublePrecision, DoublePrecision>::add_cases(suite, "FloatDPBall", dp);
    TestFloatBall<MultiplePrecision>::add_cases(suite, "FloatMPBall", MultiplePrecision(128_bits));
    TestFloatBall<MultiplePrecision, DoublePrecision>::add_cases(suite, "FloatMPDPBall", MultiplePrecision(128_bits));
    return Ariadne::Testing::run_suite(suite, argc, argv);
}

