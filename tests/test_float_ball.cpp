/***************************************************************************
 *            test_float_ball.cpp
 *
 *  Copyright  2006-20  Alberto Casagrande, Pieter Collins
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "numeric/builtin.hpp"
#include "numeric/decimal.hpp"
#include "numeric/rational.hpp"
#include "numeric/number.hpp"
#include "numeric/real.hpp"
#include "numeric/float.decl.hpp"

#include "numeric/float_ball.hpp"

#include "numeric/float_error.hpp"
#include "numeric/float_bounds.hpp"

#include "utility/test.hpp"
#include "test_floats.hpp"

using namespace Ariadne;
using namespace std;

template<class PR, class PRE=PR>
class TestFloatBall
    : public TestFloats<PR>
{
    typedef RawFloat<PR> RawFloatType;
    typedef FloatType<ApproximateTag,PR> FloatApproximationType;
    typedef FloatType<LowerTag,PR> FloatLowerBoundType;
    typedef FloatType<UpperTag,PR> FloatUpperBoundType;
    typedef FloatType<BoundedTag,PR> FloatBoundsType;
    typedef FloatBall<PR,PRE> FloatBallType;
    typedef FloatType<ExactTag,PR> FloatValueType;
  private:
    PR precision; PRE error_precision;
    static PRE _make_error_precision(PR const& pr) { if constexpr (Same<PR,PRE>) { return pr; } else { return PRE(); } }
  public:
    TestFloatBall(PR prec) : precision(prec), error_precision(_make_error_precision(prec)) { };
    TestFloatBall(PR prec, PRE err_prec) : precision(prec), error_precision(err_prec) { };
    Void test();
  private:
    using TestFloats<PR>::to_rational;
    Void test_header_api();
    Void test_precision();
    Void test_conversions();
    Void test_rounded_arithmetic();
};

template<class PR, class PRE> Void
TestFloatBall<PR,PRE>::test()
{
    ARIADNE_TEST_CALL(test_header_api());
    ARIADNE_TEST_CALL(test_precision());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_rounded_arithmetic());
}

template<class PR, class PRE> Void
TestFloatBall<PR,PRE>::test_header_api()
{
    RawFloatType raw_one(1,precision);
    RawFloatType raw_two(2,precision);
    RawFloatType raw_half(Dyadic(1,1u),precision);
    RawFloat<PRE> raw_error_one(1,error_precision);

    FloatBallType by_precision(precision);
    FloatBallType by_precisions(precision,error_precision);
    typename FloatBallType::CharacteristicsType characteristics_pair(precision,error_precision);
    FloatBallType by_characteristics(characteristics_pair);
    FloatBallType by_raw(raw_one);

    FloatBallType from_twoexp(TwoExp(1),precision);
    FloatBallType from_integer(Integer(2),precision);
    FloatBallType from_dyadic(Dyadic(3,1u),precision);
    FloatBallType from_decimal(Decimal(String("1.5")),precision);
    Real real_value(Rational(3,2));
    FloatBallType from_real(real_value,precision);

    ARIADNE_TEST_EQUALS(from_twoexp.value_raw(),RawFloatType(2,precision));
    ARIADNE_TEST_EQUALS(from_twoexp.error_raw(),RawFloat<PRE>(0,error_precision));
    ARIADNE_TEST_EQUALS(from_integer.value_raw(),RawFloatType(2,precision));
    ARIADNE_TEST_EQUALS(from_dyadic.value_raw(),RawFloatType(Dyadic(3,1u),precision));
    ARIADNE_TEST_EQUALS(from_decimal.value_raw(),RawFloatType(Dyadic(3,1u),precision));
    ARIADNE_TEST_BINARY_PREDICATE(models,from_real,Rational(3,2));

    ARIADNE_TEST_EQUALS(by_precision.value_raw(),RawFloatType(0,precision));
    ARIADNE_TEST_EQUALS(by_precisions.error_raw(),RawFloat<PRE>(0,error_precision));
    ARIADNE_TEST_EQUALS(by_characteristics.value_raw(),RawFloatType(0,precision));
    ARIADNE_TEST_EQUALS(by_raw.value_raw(),raw_one);

    ValidatedNumber validated=by_raw.operator ValidatedNumber();
    FloatBallType from_validated(validated,precision);
    ARIADNE_TEST_BINARY_PREDICATE(models,from_validated,Rational(1));

    FloatBallType assigned(precision,error_precision);
    assigned=validated;
    ARIADNE_TEST_EXECUTE(by_raw.create(validated));
    ARIADNE_TEST_EXECUTE(by_raw.generic());

    ARIADNE_TEST_EQUALS(by_raw.get_d(),1.0);

    FloatBallType interval(raw_two,raw_error_one);
    ARIADNE_TEST_EQUALS(interval.lower().raw(),RawFloatType(1,precision));
    ARIADNE_TEST_EQUALS(interval.upper().raw(),RawFloatType(3,precision));

    FloatBallType truncated=trunc(from_integer);
    FloatBallType truncated_digits=trunc(from_integer,8u);
    ARIADNE_TEST_BINARY_PREDICATE(models,truncated,Integer(2));
    ARIADNE_TEST_BINARY_PREDICATE(models,truncated_digits,Integer(2));

    StringStream input;
    input << "1" << "\xC2\xB1" << "0";
    FloatBallType parsed(precision,error_precision);
    input >> parsed;
    ARIADNE_TEST_EQUALS(parsed.value_raw(),RawFloatType(1,precision));
    ARIADNE_TEST_EQUALS(parsed.error_raw(),RawFloat<PRE>(0,error_precision));

    StringStream incomplete_input;
    incomplete_input << "1" << "\xC2\xB1";
    ARIADNE_TEST_FAIL(incomplete_input >> parsed);

    StringStream invalid_separator("1xx0");
    ARIADNE_TEST_FAIL(invalid_separator >> parsed);

    if constexpr (Same<PR,MultiplePrecision>) {
        MultiplePrecision lower_precision(53_bits);
        FloatBallType source(Rational(1,3),precision);
        FloatBallType narrowed(source,lower_precision);
        FloatBallType narrowed_explicit(source,lower_precision,error_precision);
        ARIADNE_TEST_BINARY_PREDICATE(models,narrowed,Rational(1,3));
        ARIADNE_TEST_BINARY_PREDICATE(models,narrowed_explicit,Rational(1,3));
    }

    auto characteristics=by_raw.characteristics();
    ARIADNE_TEST_EQUALS(std::get<0>(characteristics),precision);
    ARIADNE_TEST_EQUALS(std::get<1>(characteristics),error_precision);

    ARIADNE_TEST_EXECUTE(nul(by_raw));
    ARIADNE_TEST_EXECUTE(pos(by_raw));
    ARIADNE_TEST_EXECUTE(neg(by_raw));
    ARIADNE_TEST_EXECUTE(tanh(by_raw));

    FloatBallType trigonometric(raw_half);
    ARIADNE_TEST_EXECUTE(asin(trigonometric));
    ARIADNE_TEST_EXECUTE(acos(trigonometric));

    FloatBallType other(raw_two);

    FloatBoundsType bounds_one(1,precision);
    FloatBallType ball_two(raw_two);
    ARIADNE_TEST_SAME(add(bounds_one,ball_two),FloatBoundsType(3,precision));
    ARIADNE_TEST_SAME(sub(bounds_one,ball_two),FloatBoundsType(-1,precision));
    ARIADNE_TEST_SAME(mul(bounds_one,ball_two),FloatBoundsType(2,precision));
    ARIADNE_TEST_SAME(div(bounds_one,ball_two),FloatBoundsType(Dyadic(1,1u),precision));
    ARIADNE_TEST_SAME(add(ball_two,bounds_one),FloatBoundsType(3,precision));
    ARIADNE_TEST_SAME(sub(ball_two,bounds_one),FloatBoundsType(1,precision));
    ARIADNE_TEST_SAME(mul(ball_two,bounds_one),FloatBoundsType(2,precision));
    ARIADNE_TEST_SAME(div(ball_two,bounds_one),FloatBoundsType(2,precision));

    ARIADNE_TEST_EXECUTE(eq(by_raw,other));
    ARIADNE_TEST_EXECUTE(lt(by_raw,other));

    ARIADNE_TEST_EXECUTE(add(by_raw,raw_two));
    ARIADNE_TEST_EXECUTE(sub(by_raw,raw_two));
    ARIADNE_TEST_EXECUTE(mul(by_raw,raw_two));
    ARIADNE_TEST_EXECUTE(div(by_raw,raw_two));
    ARIADNE_TEST_EXECUTE(add(raw_two,by_raw));
    ARIADNE_TEST_EXECUTE(sub(raw_two,by_raw));
    ARIADNE_TEST_EXECUTE(mul(raw_two,by_raw));
    ARIADNE_TEST_EXECUTE(div(raw_two,by_raw));

    ARIADNE_TEST_EXECUTE(max(by_raw,raw_two));
    ARIADNE_TEST_EXECUTE(min(by_raw,raw_two));
    ARIADNE_TEST_EXECUTE(max(raw_two,by_raw));
    ARIADNE_TEST_EXECUTE(min(raw_two,by_raw));

    ARIADNE_TEST_EXECUTE((void)(by_raw==raw_one));
    ARIADNE_TEST_EXECUTE((void)(by_raw!=raw_two));
    ARIADNE_TEST_EXECUTE((void)(by_raw<raw_two));
    ARIADNE_TEST_EXECUTE((void)(by_raw>RawFloatType(0,precision)));
    ARIADNE_TEST_EXECUTE((void)(by_raw<=raw_one));
    ARIADNE_TEST_EXECUTE((void)(by_raw>=raw_one));

    ARIADNE_TEST_EXECUTE((void)(raw_one==by_raw));
    ARIADNE_TEST_EXECUTE((void)(raw_two!=by_raw));
    ARIADNE_TEST_EXECUTE((void)(RawFloatType(0,precision)<by_raw));
    ARIADNE_TEST_EXECUTE((void)(raw_two>by_raw));
    ARIADNE_TEST_EXECUTE((void)(raw_one<=by_raw));
    ARIADNE_TEST_EXECUTE((void)(raw_one>=by_raw));

    auto fac=factory(by_raw);
    ARIADNE_TEST_EXECUTE(fac.create(Real(Rational(1,2))));
    ARIADNE_TEST_EXECUTE(fac.create(Rational(1,2)));
    ARIADNE_TEST_EXECUTE(fac.create(Dyadic(1,1u)));
    ARIADNE_TEST_EXECUTE(fac.create(Integer(1)));

    FloatBallType wide(RawFloatType(0,precision),raw_error_one);
    auto squared=sqr(wide);
    ARIADNE_TEST_ASSERT(squared.error_raw()>=RawFloat<PRE>(0,error_precision));
}

template<class PR, class PRE> Void
TestFloatBall<PR,PRE>::test_conversions()
{
    Rational one=1;
    Rational third=one/3;

    ARIADNE_TEST_EQUALS(cast_integer(FloatBall<PR,PRE>(FloatBounds<PR>(2,3,precision))),Integer(3));
    ARIADNE_TEST_EQUALS(cast_integer(FloatBall<PR,PRE>(FloatBounds<PR>(2.25_dy,3.25_dy,precision))),Integer(3));
    ARIADNE_TEST_FAIL(cast_integer(FloatBall<PR,PRE>(FloatBounds<PR>(2.625_dy,2.875_dy,precision))));

    RawFloatType value(2,precision);
    FloatError<PRE> error(1u,error_precision);
    FloatBallType centred=value.pm(error);
    ARIADNE_TEST_EQUALS(centred.value_raw(),value);
    ARIADNE_TEST_EQUALS(centred.error_raw(),RawFloat<PRE>(1,error_precision));

    FloatBallType widened=centred.pm(error);
    ARIADNE_TEST_EQUALS(widened.value_raw(),value);
    ARIADNE_TEST_EQUALS(widened.error_raw(),RawFloat<PRE>(2,error_precision));

    if constexpr (Same<PR,DoublePrecision>) {
        ARIADNE_TEST_EQUALS(class_name<FloatBallType>(),String("FloatDPBall"));
    } else if constexpr (Same<PRE,DoublePrecision>) {
        ARIADNE_TEST_EQUALS(class_name<FloatBallType>(),String("FloatMDPBall"));
    } else {
        ARIADNE_TEST_EQUALS(class_name<FloatBallType>(),String("FloatMPBall"));
    }

    if constexpr (Same<PR,MultiplePrecision> and Same<PRE,MultiplePrecision>) {
        FloatBallType zero_with_unit_error(RawFloatType(0,precision),RawFloat<PRE>(1,error_precision));
        StringStream stream;
        stream << zero_with_unit_error;
        ARIADNE_TEST_ASSERT(not stream.str().empty());
    }
}

template<class PR, class PRE> Void
TestFloatBall<PR,PRE>::test_precision()
{
    FloatBallType mx(Rational(1),precision);
    ARIADNE_TEST_EQUALS(mx.value().precision(),precision);
    ARIADNE_TEST_EQUALS(mx.error().precision(),error_precision);
    ARIADNE_TEST_EQUALS(max(mx,mx).precision(),precision);
    ARIADNE_TEST_EQUALS(min(mx,mx).precision(),precision);
    ARIADNE_TEST_EQUALS(abs(mx).precision(),precision);
    ARIADNE_TEST_EQUALS((mx+mx).precision(),precision);
    ARIADNE_TEST_EQUALS((mx-mx).precision(),precision);
    ARIADNE_TEST_EQUALS((mx*mx).precision(),precision);
    ARIADNE_TEST_EQUALS((mx/mx).precision(),precision);
    ARIADNE_TEST_EQUALS(sqrt(mx).precision(),precision);
    ARIADNE_TEST_EQUALS(exp(mx).precision(),precision);
    ARIADNE_TEST_EQUALS(log(mx).precision(),precision);
    ARIADNE_TEST_EQUALS(sin(mx).precision(),precision);
    ARIADNE_TEST_EQUALS(cos(mx).precision(),precision);
    ARIADNE_TEST_EQUALS(tan(mx).precision(),precision);
    ARIADNE_TEST_EQUALS(atan(mx).precision(),precision);
    ARIADNE_TEST_EQUALS((mx+2u).precision(),precision);
    ARIADNE_TEST_EQUALS((mx-2u).precision(),precision);
    ARIADNE_TEST_EQUALS((mx*2u).precision(),precision);
    ARIADNE_TEST_EQUALS((mx/2u).precision(),precision);
    ARIADNE_TEST_EQUALS((mx+2).precision(),precision);
    ARIADNE_TEST_EQUALS((mx-2).precision(),precision);
    ARIADNE_TEST_EQUALS((mx*2).precision(),precision);
    ARIADNE_TEST_EQUALS((mx/2).precision(),precision);
}

template<class PR, class PRE> Void
TestFloatBall<PR,PRE>::test_rounded_arithmetic()
{
    Rational one=1;
    Rational three=3;
    Rational five=5;
    Rational six=6;
    Rational third=one/three;
    Rational fifth=one/five;

    ARIADNE_TEST_BINARY_PREDICATE(models,FloatBallType(third,precision)+FloatBallType(fifth,precision),third+fifth);
    ARIADNE_TEST_BINARY_PREDICATE(models,FloatBallType(third,precision)-FloatBallType(fifth,precision),third-fifth);
    ARIADNE_TEST_BINARY_PREDICATE(models,FloatBallType(third,precision)*FloatBallType(fifth,precision),third*fifth);
    ARIADNE_TEST_BINARY_PREDICATE(models,FloatBallType(third,precision)/FloatBallType(fifth,precision),third/fifth);
//    ARIADNE_TEST_BINARY_PREDICATE(models,1/FloatBallType(three,precision)/FloatBallType(fifth,precision),1/three);

    ARIADNE_TEST_BINARY_PREDICATE(models,pow(FloatBallType(three/five,precision),4u),pow(three/five,4u));
    ARIADNE_TEST_BINARY_PREDICATE(models,pow(FloatBallType(three/five,precision),4),pow(three/five,4));
    ARIADNE_TEST_BINARY_PREDICATE(models,pow(FloatBallType(three/five,precision),-4),pow(three/five,-4));
    ARIADNE_TEST_BINARY_PREDICATE(models,pow(FloatBallType(three/five,precision),-7),pow(three/five,-7));
}





Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestFloatBall<DoublePrecision,DoublePrecision>(dp).test();
    TestFloatBall<MultiplePrecision>(MultiplePrecision(128_bits)).test();
    TestFloatBall<MultiplePrecision,DoublePrecision>(MultiplePrecision(128_bits)).test();

    return ARIADNE_TEST_FAILURES;
}

