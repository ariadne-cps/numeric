/***************************************************************************
 *            test_float.cpp
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
#include <cmath>
#include <limits>

#include "numeric/builtin.hpp"
#include "numeric/double.hpp"
#include "numeric/decimal.hpp"
#include "numeric/rational.hpp"
#include "numeric/number.hpp"
#include "numeric/upper_number.hpp"
#include "numeric/lower_number.hpp"
#include "numeric/float.decl.hpp"


#include "numeric/float_ball.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_lower_bound.hpp"
#include "numeric/float_upper_bound.hpp"
#include "numeric/float_error.hpp"
#include "numeric/float_literals.hpp"
#include "numeric/positive.hpp"
#include "numeric/rounded_float.hpp"
#include "numeric/casts.hpp"

#include "utility/test.hpp"
#include "test_floats.hpp"

using namespace Ariadne;
using namespace std;

Void test_float_literals()
{
    FloatDP half(ExactDouble(0.5),dp);

    auto error=0.5_error;
    auto exact=0.5_exact;
    auto near_value=0.5_near;
    auto upper=0.5_upper;
    auto lower=0.5_lower;
    auto approximation=0.5_approx;

    ARIADNE_TEST_EQUALS(error.raw(),half);
    ARIADNE_TEST_EQUALS(exact,half);
    ARIADNE_TEST_EQUALS(near_value.value_raw(),half);
    ARIADNE_TEST_EQUALS(near_value.error_raw(),FloatDP(0,dp));
    ARIADNE_TEST_EQUALS(upper.raw(),half);
    ARIADNE_TEST_EQUALS(lower.raw(),half);
    ARIADNE_TEST_EQUALS(approximation.raw(),half);

    DecimalPlaces places(3u);
    DecimalPrecision precision(5u);
    ARIADNE_TEST_EQUALS(static_cast<unsigned int>(places),3u);
    ARIADNE_TEST_EQUALS(static_cast<unsigned int>(precision),5u);

    if constexpr (std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits) {
        ARIADNE_TEST_FAIL(operator""_error(0.1L));
        ARIADNE_TEST_FAIL(operator""_exact(0.1L));

        long double near_input=0.1L;
        auto near_inexact=operator""_near(near_input);
        long double near_value_raw=static_cast<long double>(near_inexact.value_raw().get_d());
        long double near_error_raw=static_cast<long double>(near_inexact.error_raw().get_d());
        ARIADNE_TEST_ASSERT(std::abs(near_value_raw-near_input)<=near_error_raw);

        long double upper_input=0.3L;
        auto upper_inexact=operator""_upper(upper_input);
        ARIADNE_TEST_ASSERT(static_cast<long double>(upper_inexact.raw().get_d())>=upper_input);

        long double lower_input=0.1L;
        auto lower_inexact=operator""_lower(lower_input);
        ARIADNE_TEST_ASSERT(static_cast<long double>(lower_inexact.raw().get_d())<=lower_input);
    }
}


template<class PR>
class TestFloat
{
    using PRE=DoublePrecision;
    typedef RawFloat<PR> RawFloatType;
    typedef FloatBounds<PR> FloatBoundsType;
    typedef FloatBall<PR,PRE> FloatBallType;
    typedef Float<PR> FloatType;
  private:
    PR precision;
  public:
    TestFloat(PR prec) : precision(prec) { }
    Void test();
  private:
    Void test_concept();
    Void test_conversions();
    Void test_operations();
    Void test_predicates();
};

template<class PR> Void
TestFloat<PR>::test()
{
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_operations());
    ARIADNE_TEST_CALL(test_predicates());
}

template<class PR> Void
TestFloat<PR>::test_concept()
{
    Float<DoublePrecision>::set_output_places(17);

    PR pr=precision;
    PRE pre;

    Boolean b;
    Nat m=1u;
    Int n=1;
    Integer z=1;
    Dyadic w=1;
    ExactDouble d(1.0);
    TwoExp t(0);
    RawFloatType f(pr);
    FloatType vx(pr);
    FloatType rx(pr);
    FloatBoundsType rbx(pr);

    // Constructors
    rx=FloatType(m,pr); rx=FloatType(n,pr); rx=FloatType(z,pr); rx=FloatType(w,pr);
    rx=FloatType(d,pr); rx=FloatType(t,pr);
    rx=FloatType(pr); rx=FloatType(f);

    // Assignment
    rx=m; rx=n; rx=z; rx=w; rx=d; rx=t;

    // Arithmetic operators
    rx=operator+(vx); rx=operator-(vx); rx=t*vx; rx=vx*t; rx=vx/t;
    rbx=operator+(vx,vx); rbx=operator-(vx,vx); rbx=operator*(vx,vx); rbx=operator/(vx,vx);

    // Exact operations
    rx=nul(vx); rx=pos(vx); rx=neg(vx); rx=hlf(vx);
    rx=mul(vx,t); rx=div(vx,t);

    FloatBall<PR,PRE> rmx=add(vx,vx,pre); rmx=sub(vx,vx,pre); rmx=mul(vx,vx,pre); rmx=div(vx,vx,pre);

    // Arithmetic
    rbx=add(vx,vx); rbx=sub(vx,vx); rbx=mul(vx,vx); rbx=div(vx,vx);
    rbx=sqr(vx); rbx=rec(vx); rbx=pow(vx,m); rbx=pow(vx,n);

    // Order
    rx=max(vx,vx); rx=min(vx,vx); rx=abs(vx);

    // Comparisons
    b=(vx==vx); b=(vx!=vx); b=(vx<=vx); b=(vx>=vx); b=(vx< vx); b=(vx> vx);
}

template<class PR> Void
TestFloat<PR>::test_conversions()
{
    ARIADNE_TEST_EQUALS(cast_integer(Float<PR>(Dyadic(2),precision)),Integer(2));
    ARIADNE_TEST_FAIL(cast_integer(Float<PR>(Dyadic(7,2u),precision)));

    Positive<Float<PR>> positive_value(2u,precision);
    ARIADNE_TEST_EQUALS(cast_unsigned(positive_value),Float<PR>(2u,precision));

    FloatError<PR> error(1u,precision);
    ARIADNE_TEST_EQUALS(cast_exact(error),error.raw());

    FloatApproximation<PR> approximation(2u,precision);
    ARIADNE_TEST_EQUALS(cast_exact(approximation),approximation.raw());

    if constexpr (Same<PR,DoublePrecision>) {
        RawFloatDP raw(2u,dp);
        FloatDPApproximation approximate(2u,dp);
        ARIADNE_TEST_EQUALS(cast_raw(raw),raw);
        ARIADNE_TEST_EQUALS(cast_raw(approximate),approximate.raw());
        ARIADNE_TEST_EQUALS(cast_approximate(raw).raw(),raw);
        ARIADNE_TEST_EQUALS(cast_approximate(approximate).raw(),approximate.raw());

        ARIADNE_TEST_EXECUTE((RoundedFloatDP{approximate}));
        ARIADNE_TEST_EQUALS(FloatDP(String("1.25"),dp),Dyadic(5,2u));
        ARIADNE_TEST_EQUALS(FloatDP(TwoExp(10),dp),Dyadic(TwoExp(10)));
        ARIADNE_TEST_FAIL((FloatDP{TwoExp(1024),dp}));

        Dyadic too_precise{Integer(String("9007199254740993"))};
        ARIADNE_TEST_FAIL((FloatDP{too_precise,dp}));

        FloatDP exact_assignment(0,dp);
        ARIADNE_TEST_EQUALS((exact_assignment=Dyadic(2)),Dyadic(2));
        ARIADNE_TEST_FAIL(exact_assignment=too_precise);
        ARIADNE_TEST_EXECUTE(exact_assignment=Dyadic::nan());
        ARIADNE_TEST_ASSERT(is_nan(exact_assignment));
        ARIADNE_TEST_EXECUTE(exact_assignment=Dyadic::inf(Sign::POSITIVE));
        ARIADNE_TEST_ASSERT(is_inf(exact_assignment));
        ARIADNE_TEST_EXECUTE(exact_assignment=Dyadic::inf(Sign::NEGATIVE));
        ARIADNE_TEST_ASSERT(is_inf(exact_assignment));

        ARIADNE_TEST_EQUALS(integer_cast<Nat>(FloatDP(3u,dp)),Nat(3u));
        ARIADNE_TEST_EQUALS(integer_cast<Int>(FloatDP(-3,dp)),Int(-3));

        FloatDP resized(1,dp);
        resized.set_precision(dp);
        ARIADNE_TEST_EQUALS(resized.precision(),dp);

        Float32 f32(FloatDP(1.25_x,dp),FloatDP::ROUND_TO_NEAREST);
        ARIADNE_TEST_EQUALS(static_cast<FloatDP>(f32),FloatDP(1.25_x,dp));

        ARIADNE_TEST_EQUALS(class_name<double>(),String("double"));
        ARIADNE_TEST_EQUALS(class_name<ApproximateDouble>(),String("ApproximateDouble"));
        ARIADNE_TEST_EQUALS(class_name<ExactDouble>(),String("ExactDouble"));
        ARIADNE_TEST_EQUALS(class_name<Rounded<FloatDP>>(),String("Rounded<FloatDP>"));
    } else if constexpr (Same<PR,MultiplePrecision>) {
        MultiplePrecision pr=precision;
        FloatMP from_string(String("1.25"),pr);
        ARIADNE_TEST_EQUALS(from_string,Dyadic(5,2u));

        FloatMP source(Dyadic(3,2u),pr);
        FloatMP from_raw(source.get_mpfr(),RawPtr());
        ARIADNE_TEST_EQUALS(from_raw,source);

        Dyadic too_wide{Integer(String("18446744073709551617"))};
        ARIADNE_TEST_FAIL(FloatMP(too_wide,MultiplePrecision(64_bits)));
        ARIADNE_TEST_ASSERT(is_nan(FloatMP(Dyadic::nan(),pr)));

        FloatMP exact_assignment(0,MultiplePrecision(64_bits));
        ARIADNE_TEST_EQUALS((exact_assignment=Dyadic(2)),Dyadic(2));
        ARIADNE_TEST_FAIL(exact_assignment=too_wide);

        FloatMP zero(pr);
        ARIADNE_TEST_ASSERT(is_zero(zero));
        ARIADNE_TEST_ASSERT(is_inf(FloatMP::inf(pr)));
        ARIADNE_TEST_ASSERT(is_nan(FloatMP::nan(pr)));

        MultiplePrecision old_default=FloatMP::get_default_precision();
        FloatMP::set_default_precision(MultiplePrecision(96_bits));
        ARIADNE_TEST_EQUALS(FloatMP::get_default_precision(),MultiplePrecision(96_bits));
        FloatMP::set_default_precision(old_default);

        FloatMP resized(Dyadic(3,2u),pr);
        resized.set_precision(MultiplePrecision(96_bits));
        ARIADNE_TEST_EQUALS(resized.precision(),MultiplePrecision(96_bits));
        ARIADNE_TEST_EXECUTE(mpfr_set_si(resized.get_mpfr(),2,MPFR_RNDN));
        ARIADNE_TEST_EQUALS(resized,Dyadic(2));
    }
}

template<class PR> Void
TestFloat<PR>::test_operations()
{
    PR pr=precision;
    DoublePrecision pre;
    FloatType vr(pr);

    Nat m(5);
    Int n(-3);
    TwoExp t(-12);
    Dyadic w(-3,1u);
    FloatType vx(w,pr);
    Dyadic w1(-3,1u), w2(5,2u);
    FloatType vx1(w1,pr), vx2(w2,pr);

    ARIADNE_TEST_EQUALS(vx,w);

    ARIADNE_TEST_EQUALS(FloatType(RawFloatType(1.25_x,pr)),1.25_dy);

    ARIADNE_TEST_EQUALS(FloatType(3u,pr),3.0_dy);
    ARIADNE_TEST_EQUALS(FloatType(-5,pr),-5.0_dy);
    ARIADNE_TEST_EQUALS(FloatType(ExactDouble(1.25),pr),1.25_dy);
    ARIADNE_TEST_EQUALS(FloatType(TwoExp(-3),pr),0.125_dy);
    ARIADNE_TEST_EQUALS(FloatType(Integer(-23),pr),-23.0_dy);
    ARIADNE_TEST_EQUALS(FloatType(Dyadic(-23,3u),pr),-2.875_dy);


    ARIADNE_TEST_EQUALS((vr=3u),3.0_dy);
    ARIADNE_TEST_EQUALS((vr=-5),-5.0_dy);
    ARIADNE_TEST_EQUALS((vr=Integer(-23)),-23.0_dy);
    ARIADNE_TEST_EQUALS((vr=TwoExp(-3)),0.125_dy);
    ARIADNE_TEST_EQUALS((vr=Dyadic(-23,3u)),-2.875_dy);

    ARIADNE_TEST_EQUALS(FloatType(w,pr).operator Dyadic(),w);
    ARIADNE_TEST_EQUALS(FloatType(w,pr).operator Rational(),w);

    ARIADNE_TEST_EQUALS(nul(vx),nul(w));
    ARIADNE_TEST_EQUALS(pos(vx),pos(w));
    ARIADNE_TEST_EQUALS(neg(vx),neg(w));
    ARIADNE_TEST_EQUALS(hlf(vx),hlf(w));

    ARIADNE_TEST_EQUALS(mul(vx,t),mul(w,t));
    ARIADNE_TEST_EQUALS(div(vx,t),w/t);

    ARIADNE_TEST_BINARY_PREDICATE(models,sqr(vx),sqr(w));
    ARIADNE_TEST_BINARY_PREDICATE(models,rec(vx),rec(w));

//    ARIADNE_TEST_BINARY_PREDICATE(models,sqr(vx,pre),sqr(w));
//    ARIADNE_TEST_BINARY_PREDICATE(models,rec(vx,pre),rec(w));

    ARIADNE_TEST_EQUALS(max(vx1,vx2),max(w1,w2));
    ARIADNE_TEST_EQUALS(min(vx1,vx2),min(w1,w2));
    ARIADNE_TEST_EQUALS(abs(vx),abs(w));
    ARIADNE_TEST_EQUALS(mig(vx).raw(),abs(w));
    ARIADNE_TEST_EQUALS(mag(vx).raw(),abs(w));
    ARIADNE_TEST_SAME(mag(vx),PositiveFloatUpperBound<PR>(abs(w),pr));
    ARIADNE_TEST_SAME(mig(vx),PositiveFloatLowerBound<PR>(abs(w),pr));

    ARIADNE_TEST_EQUALS(t*vx,t*w);
    ARIADNE_TEST_EQUALS(vx*t,w*t);
    ARIADNE_TEST_EQUALS(vx/t,w/t);

    ARIADNE_TEST_BINARY_PREDICATE(models,vx1+vx2,w1+w2);
    ARIADNE_TEST_BINARY_PREDICATE(models,vx1-vx2,w1-w2);
    ARIADNE_TEST_BINARY_PREDICATE(models,vx1*vx2,w1*w2);
    ARIADNE_TEST_BINARY_PREDICATE(models,vx1/vx2,w1/w2);

    ARIADNE_TEST_BINARY_PREDICATE(models,add(vx1,vx2),add(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,sub(vx1,vx2),sub(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,mul(vx1,vx2),mul(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,div(vx1,vx2),div(w1,w2));

    ARIADNE_TEST_BINARY_PREDICATE(models,add(vx1,vx2,pre),add(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,sub(vx1,vx2,pre),sub(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,mul(vx1,vx2,pre),mul(w1,w2));
    ARIADNE_TEST_BINARY_PREDICATE(models,div(vx1,vx2,pre),div(w1,w2));


    ARIADNE_TEST_SAME(add(vx1,vx2,pre),FloatBallType(add(w1,w2),pr,pre));

    ARIADNE_TEST_BINARY_PREDICATE(models,pow(vx,m),pow(w,m));
    ARIADNE_TEST_BINARY_PREDICATE(models,pow(vx,n),pow(Rational(w),n));

//    friend Bounds<F> med(F const& x1, F const& x2);
//    friend Bounds<F> rad(F const& x1, F const& x2);

    ARIADNE_TEST_BINARY_PREDICATE(models,sqr(sqrt(abs(vx))),abs(w));
    ARIADNE_TEST_BINARY_PREDICATE(models,log(exp(vx)),w);
    ARIADNE_TEST_BINARY_PREDICATE(models,exp(log(abs(vx))),abs(w));
    ARIADNE_TEST_BINARY_PREDICATE(models,atan(tan(vx)),w);
    ARIADNE_TEST_BINARY_PREDICATE(models,tan(atan(vx)),w);

    ARIADNE_TEST_SAME(sin(vx),sin(FloatBounds<PR>(vx)));
    ARIADNE_TEST_SAME(cos(vx),cos(FloatBounds<PR>(vx)));
    ARIADNE_TEST_SAME(tan(vx),tan(FloatBounds<PR>(vx)));
    ARIADNE_TEST_SAME(tanh(vx),tanh(FloatBounds<PR>(vx)));
    ARIADNE_TEST_SAME(atan(vx),atan(FloatBounds<PR>(vx)));

    FloatBoundsType tanh_wide(-100,100,pr);
    FloatBoundsType tanh_positive(1,2,pr);
    FloatBoundsType tanh_negative(-2,-1,pr);
    auto wide_image=tanh(tanh_wide);
    auto positive_image=tanh(tanh_positive);
    auto negative_image=tanh(tanh_negative);
    ARIADNE_TEST_ASSERT(wide_image.lower_raw()>=FloatType(-1,pr));
    ARIADNE_TEST_ASSERT(wide_image.upper_raw()<=FloatType(1,pr));
    ARIADNE_TEST_ASSERT(wide_image.lower_raw()<FloatType(0,pr));
    ARIADNE_TEST_ASSERT(wide_image.upper_raw()>FloatType(0,pr));
    ARIADNE_TEST_ASSERT(positive_image.lower_raw()>FloatType(0,pr));
    ARIADNE_TEST_ASSERT(positive_image.upper_raw()<FloatType(1,pr));
    ARIADNE_TEST_ASSERT(negative_image.lower_raw()>FloatType(-1,pr));
    ARIADNE_TEST_ASSERT(negative_image.upper_raw()<FloatType(0,pr));

    FloatBallType tanh_ball(FloatBoundsType(-1,1,pr),pre);
    auto tanh_ball_image=tanh(tanh_ball);
    ARIADNE_TEST_ASSERT(tanh_ball_image.lower_raw()<=FloatType(0,pr));
    ARIADNE_TEST_ASSERT(tanh_ball_image.upper_raw()>=FloatType(0,pr));

    FloatLowerBound<PR> tanh_lower(-2,pr);
    FloatUpperBound<PR> tanh_upper(2,pr);
    ARIADNE_TEST_ASSERT(tanh(tanh_lower).raw()>FloatType(-1,pr));
    ARIADNE_TEST_ASSERT(tanh(tanh_upper).raw()<FloatType(1,pr));

//    ARIADNE_TEST_SAME(shft(vx,n),shft(w,n));

    if constexpr (Same<PR,DoublePrecision>) {
        FloatDP one(1,dp);
        FloatDP half(Dyadic(1,1u),dp);
        CurrentRoundingMode current;
        FloatDP::set_rounding_to_nearest();

        ARIADNE_TEST_EQUALS(pow_rnd(FloatDP(2,dp),Int(-2)),FloatDP(Dyadic(1,2u),dp));
        ARIADNE_TEST_EXECUTE(sqrt_rnd(one));
        ARIADNE_TEST_EXECUTE(exp_rnd(one));
        ARIADNE_TEST_EXECUTE(log_rnd(one));
        ARIADNE_TEST_EXECUTE(sin_rnd(one));
        ARIADNE_TEST_EXECUTE(cos_rnd(one));
        ARIADNE_TEST_EXECUTE(tan_rnd(one));
        ARIADNE_TEST_EXECUTE(atan_rnd(one));

        ARIADNE_TEST_EXECUTE(sqr(current,one));
        ARIADNE_TEST_EXECUTE(rec(current,one));
        ARIADNE_TEST_EXECUTE(add(current,one,half));
        ARIADNE_TEST_EXECUTE(sub(current,one,half));
        ARIADNE_TEST_EXECUTE(mul(current,one,half));
        ARIADNE_TEST_EXECUTE(div(current,one,half));
        ARIADNE_TEST_EXECUTE(fma(current,one,half,one));
        ARIADNE_TEST_EXECUTE(pow(current,one,Nat(2u)));
        ARIADNE_TEST_EXECUTE(pow(current,one,Int(-2)));
        ARIADNE_TEST_EXECUTE(sqrt(current,one));
        ARIADNE_TEST_EXECUTE(exp(current,one));
        ARIADNE_TEST_EXECUTE(log(current,one));
        ARIADNE_TEST_EXECUTE(sin(current,one));
        ARIADNE_TEST_EXECUTE(cos(current,one));
        ARIADNE_TEST_EXECUTE(tan(current,one));
        ARIADNE_TEST_EXECUTE(asin(current,half));
        ARIADNE_TEST_EXECUTE(acos(current,half));
        ARIADNE_TEST_EXECUTE(atan(current,one));
        ARIADNE_TEST_EXECUTE(FloatDP::pi(current,dp));

        ARIADNE_TEST_EQUALS(abs(FloatDP::ROUND_TO_NEAREST,FloatDP(-2,dp)),FloatDP(2,dp));
        ARIADNE_TEST_EQUALS(mag(FloatDP::ROUND_TO_NEAREST,FloatDP(-2,dp)),FloatDP(2,dp));
        ARIADNE_TEST_ASSERT(same(FloatDP(2,dp),FloatDP(2,dp)));
        ARIADNE_TEST_ASSERT(not same(FloatDP(2,dp),FloatDP(3,dp)));

        FloatDP scaled(2,dp);
        scaled*=TwoExp(2);
        ARIADNE_TEST_EQUALS(scaled,FloatDP(8,dp));
        scaled/=TwoExp(3);
        ARIADNE_TEST_EQUALS(scaled,FloatDP(1,dp));

        ARIADNE_TEST_EXECUTE(pos_opp(one));
        ARIADNE_TEST_EXECUTE(neg_opp(one));
        ARIADNE_TEST_EXECUTE(sqr_opp(one));
        ARIADNE_TEST_EXECUTE(rec_opp(one));
        ARIADNE_TEST_EXECUTE(add_opp(one,half));
        ARIADNE_TEST_EXECUTE(sub_opp(one,half));
        ARIADNE_TEST_EXECUTE(mul_opp(one,half));
        ARIADNE_TEST_EXECUTE(div_opp(one,half));

        Dyadic downward_target=Dyadic(1)+Dyadic(3,54u);
        FloatDP downward_value(downward_target,FloatDP::ROUND_DOWNWARD,dp);
        ARIADNE_TEST_ASSERT(Dyadic(downward_value)<=downward_target);

        Dyadic negative_downward_target=-Dyadic(1)-Dyadic(3,54u);
        FloatDP negative_downward_value(negative_downward_target,FloatDP::ROUND_DOWNWARD,dp);
        ARIADNE_TEST_ASSERT(Dyadic(negative_downward_value)<=negative_downward_target);

        ARIADNE_TEST_EXECUTE(FloatDP::pi(FloatDP::ROUND_UPWARD,dp));
        ARIADNE_TEST_EXECUTE(FloatDP::pi(FloatDP::ROUND_DOWNWARD,dp));
        ARIADNE_TEST_EXECUTE(FloatDP::pi(FloatDP::ROUND_TO_NEAREST,dp));
        ARIADNE_TEST_FAIL(FloatDP::pi(static_cast<FloatDP::RoundingModeType>(0xffffu),dp));

        ARIADNE_TEST_ASSERT(not one.literal().empty());
        ARIADNE_TEST_ASSERT(not one.literal(FloatDP::ROUND_UPWARD).empty());
        std::ostringstream dp_places_stream;
        write(dp_places_stream,one,DecimalPlaces(3u),FloatDP::ROUND_DOWNWARD);
        ARIADNE_TEST_ASSERT(not dp_places_stream.str().empty());
        std::ostringstream dp_repr_stream;
        repr(dp_repr_stream,one,FloatDP::ROUND_UPWARD);
        ARIADNE_TEST_ASSERT(not dp_repr_stream.str().empty());
        FloatDP::set_output_places(7u);
        ARIADNE_TEST_EQUALS(FloatDP::output_places,Nat(7u));
        ARIADNE_TEST_FAIL(one.literal(FloatDP::ROUND_TOWARD_ZERO));

        DoublePrecision p1;
        DoublePrecision p2;
        ARIADNE_TEST_EXECUTE(max(p1,p2));
        ARIADNE_TEST_ASSERT(p1<=p2);
        std::ostringstream precision_repr;
        repr(precision_repr,p1);
        ARIADNE_TEST_EQUALS(precision_repr.str(),std::string("DoublePrecision()"));

        FloatDP noinit(NoInit{});
        noinit.raw()=FloatDP(2,dp);
        ARIADNE_TEST_EQUALS(noinit,FloatDP(2,dp));
        ARIADNE_TEST_EXECUTE((FloatDP{one,dp}));
        ARIADNE_TEST_ASSERT(is_inf(FloatDP::inf(dp)));
        ARIADNE_TEST_ASSERT(is_finite(one));
        ARIADNE_TEST_ASSERT(is_zero(FloatDP(0,dp)));

        std::istringstream dp_input("+2.5:");
        FloatDP parsed(dp);
        dp_input >> parsed;
        ARIADNE_TEST_EQUALS(parsed,FloatDP(Dyadic(5,1u),dp));
        ARIADNE_TEST_EQUALS(dp_input.peek(),static_cast<int>(':'));

        std::istringstream dp_negative("-.5");
        dp_negative >> parsed;
        ARIADNE_TEST_EQUALS(parsed,FloatDP(Dyadic(-1,1u),dp));

        std::istringstream dp_exp("2E+3");
        dp_exp >> parsed;
        ARIADNE_TEST_EQUALS(parsed,FloatDP(2000,dp));

        std::istringstream dp_exp_delimited("2e3:");
        dp_exp_delimited >> parsed;
        ARIADNE_TEST_EQUALS(parsed,FloatDP(2000,dp));
        ARIADNE_TEST_EQUALS(dp_exp_delimited.peek(),static_cast<int>(':'));

        std::istringstream dp_exp_negative("2e-3");
        ARIADNE_TEST_EXECUTE(dp_exp_negative >> parsed);

        std::istringstream dp_overflow("1e9999");
        ARIADNE_TEST_EXECUTE(dp_overflow >> parsed);
        ARIADNE_TEST_ASSERT(dp_overflow.fail());

        std::istringstream dp_no_digits(".");
        ARIADNE_TEST_EXECUTE(dp_no_digits >> parsed);
        ARIADNE_TEST_ASSERT(dp_no_digits.fail());

        std::istringstream dp_bad_exp("2e+");
        ARIADNE_TEST_EXECUTE(dp_bad_exp >> parsed);
        ARIADNE_TEST_ASSERT(dp_bad_exp.fail());

        std::istringstream dp_prefailed("1");
        dp_prefailed.setstate(std::ios::failbit);
        ARIADNE_TEST_EXECUTE(dp_prefailed >> parsed);
    }

    if constexpr (Same<PR,MultiplePrecision>) {
        MultiplePrecision mpr=precision;
        FloatMP a(2,mpr);
        FloatMP b(1,mpr);
        FloatDP dpa(1,dp);
        ExactDouble ed(0.5);
        CurrentRoundingMode current;

        FloatMP::set_rounding_toward_zero();
        ARIADNE_TEST_EQUALS(fma(current,a,b,b),FloatMP(3,mpr));
        ARIADNE_TEST_EQUALS(pow(current,a,Nat(2u)),FloatMP(4,mpr));
        ARIADNE_TEST_EXECUTE(tan(current,b));
        ARIADNE_TEST_EXECUTE(asin(current,b));
        ARIADNE_TEST_EXECUTE(acos(current,b));
        ARIADNE_TEST_EXECUTE(FloatMP::pi(current,mpr));

        ARIADNE_TEST_EQUALS(add(current,a,dpa),FloatMP(3,mpr));
        ARIADNE_TEST_EQUALS(sub(current,a,dpa),FloatMP(1,mpr));
        ARIADNE_TEST_EQUALS(mul(current,a,dpa),FloatMP(2,mpr));
        ARIADNE_TEST_EQUALS(div(current,a,dpa),FloatMP(2,mpr));
        ARIADNE_TEST_EQUALS(add(current,dpa,a),FloatMP(3,mpr));
        ARIADNE_TEST_EQUALS(sub(current,dpa,a),FloatMP(-1,mpr));
        ARIADNE_TEST_EQUALS(mul(current,dpa,a),FloatMP(2,mpr));
        ARIADNE_TEST_EQUALS(div(current,dpa,a),FloatMP(Dyadic(1,1u),mpr));

        ARIADNE_TEST_EQUALS(add(current,a,ed),FloatMP(Dyadic(5,1u),mpr));
        ARIADNE_TEST_EQUALS(sub(current,a,ed),FloatMP(Dyadic(3,1u),mpr));
        ARIADNE_TEST_EQUALS(mul(current,a,ed),FloatMP(1,mpr));
        ARIADNE_TEST_EQUALS(div(current,a,ed),FloatMP(4,mpr));
        ARIADNE_TEST_EQUALS(add(current,ed,a),FloatMP(Dyadic(5,1u),mpr));
        ARIADNE_TEST_EQUALS(sub(current,ed,a),FloatMP(Dyadic(-3,1u),mpr));
        ARIADNE_TEST_EQUALS(mul(current,ed,a),FloatMP(1,mpr));
        ARIADNE_TEST_EQUALS(div(current,ed,a),FloatMP(Dyadic(1,2u),mpr));

        ARIADNE_TEST_EQUALS(abs(MPFR_RNDN,FloatMP(-2,mpr)),FloatMP(2,mpr));
        ARIADNE_TEST_EQUALS(mag(MPFR_RNDN,FloatMP(-2,mpr)),FloatMP(2,mpr));

        FloatMP inplace(2,mpr);
        ARIADNE_TEST_EQUALS(iadd(MPFR_RNDN,inplace,FloatMP(1,mpr)),FloatMP(3,mpr));
        ARIADNE_TEST_EQUALS(isub(MPFR_RNDN,inplace,FloatMP(1,mpr)),FloatMP(2,mpr));
        ARIADNE_TEST_EQUALS(imul(MPFR_RNDN,inplace,FloatMP(3,mpr)),FloatMP(6,mpr));
        ARIADNE_TEST_EQUALS(idiv(MPFR_RNDN,inplace,FloatMP(2,mpr)),FloatMP(3,mpr));

        ARIADNE_TEST_ASSERT(not a.literal().empty());
        ARIADNE_TEST_ASSERT(not a.literal(MPFR_RNDN).empty());
        std::ostringstream places_stream;
        write(places_stream,a,DecimalPlaces(3u),MPFR_RNDN);
        ARIADNE_TEST_ASSERT(not places_stream.str().empty());
        std::ostringstream repr_stream;
        repr(repr_stream,a,MPFR_RNDN);
        ARIADNE_TEST_ASSERT(not repr_stream.str().empty());

        ARIADNE_TEST_EQUALS(print(a,DecimalPlaces(0u),MPFR_RNDN),String("2."));
        ARIADNE_TEST_EQUALS(print(FloatMP::nan(mpr),DecimalPrecision(3u),MPFR_RNDN),String("nan"));

        std::istringstream float_input("\t\n+2.5:");
        FloatMP parsed(mpr);
        float_input >> parsed;
        ARIADNE_TEST_EQUALS(parsed,FloatMP(Dyadic(5,1u),mpr));
        ARIADNE_TEST_EQUALS(float_input.peek(),static_cast<int>(':'));

        FloatMP::set_output_places(7u);
        ARIADNE_TEST_EQUALS(FloatMP::output_places,Nat(7u));
        ARIADNE_TEST_EQUALS(class_name<Rounded<FloatMP>>(),String("Rounded<FloatMP>"));

        ARIADNE_TEST_EXECUTE(FloatDP(a,FloatDP::ROUND_TO_NEAREST,dp));
        ARIADNE_TEST_EXECUTE(FloatDP(a,FloatDP::ROUND_TOWARD_ZERO,dp));
        ARIADNE_TEST_FAIL(FloatDP(a,static_cast<FloatDP::RoundingModeType>(0xffffu),dp));

        FloatMP::set_rounding_to_nearest();
    }
}

template<class PR> Void
TestFloat<PR>::test_predicates()
{
    PR pr=precision;
    ARIADNE_TEST_BINARY_PREDICATE(eq,Float<PR>(-1,pr),Float<PR>(-1,pr));
    ARIADNE_TEST_BINARY_PREDICATE(not eq,Float<PR>(-1,pr),Float<PR>(-2,pr));
    ARIADNE_TEST_BINARY_PREDICATE(lt,Float<PR>(-2,pr),Float<PR>(-1,pr));
    ARIADNE_TEST_BINARY_PREDICATE(not lt,Float<PR>(-2,pr),Float<PR>(-2,pr));
    ARIADNE_TEST_BINARY_PREDICATE(not lt,Float<PR>(-1,pr),Float<PR>(-2,pr));

    Dyadic w(3,1u);
    Integer zl(1), zu(2);
    Dyadic wl(5,2u), wu(7,2u);
    Rational ql(4,3), qu(5,3);

    ARIADNE_TEST_EQUALS(cmp(Float<PR>(w,pr),w),cmp(w,w));
    ARIADNE_TEST_EQUALS(cmp(Float<PR>(w,pr),zl),cmp(w,zl));
    ARIADNE_TEST_EQUALS(cmp(Float<PR>(w,pr),zu),cmp(w,zu));
    ARIADNE_TEST_EQUALS(cmp(Float<PR>(w,pr),wl),cmp(w,wl));
    ARIADNE_TEST_EQUALS(cmp(Float<PR>(w,pr),wu),cmp(w,wu));
    ARIADNE_TEST_EQUALS(cmp(Float<PR>(w,pr),ql),cmp(w,ql));
    ARIADNE_TEST_EQUALS(cmp(Float<PR>(w,pr),qu),cmp(w,qu));

    if constexpr (Same<PR,DoublePrecision>) {
        Rational zero(0);
        ARIADNE_TEST_EQUALS(cmp(FloatDP::inf(Sign::POSITIVE,dp),zero),Comparison::GREATER);
        ARIADNE_TEST_EQUALS(cmp(FloatDP::inf(Sign::NEGATIVE,dp),zero),Comparison::LESS);
        ARIADNE_TEST_EQUALS(cmp(zero,FloatDP::inf(Sign::POSITIVE,dp)),Comparison::LESS);
        ARIADNE_TEST_EQUALS(cmp(zero,FloatDP::inf(Sign::NEGATIVE,dp)),Comparison::GREATER);
    }

    if constexpr (Same<PR,MultiplePrecision>) {
        FloatMP x(Dyadic(3,1u),pr);
        FloatMP one(1,pr);
        FloatMP two_mp(2,pr);

        ARIADNE_TEST_EQUALS(cmp(two_mp,one),Comparison::GREATER);
        ARIADNE_TEST_EQUALS(cmp(two_mp,Nat(1u)),Comparison::GREATER);
        ARIADNE_TEST_EQUALS(cmp(one,Nat(2u)),Comparison::LESS);

        ARIADNE_TEST_EQUALS(cmp(Int(1),two_mp),Comparison::LESS);
        ARIADNE_TEST_EQUALS(cmp(Int(2),two_mp),Comparison::EQUAL);
        ARIADNE_TEST_EQUALS(cmp(Int(3),two_mp),Comparison::GREATER);

        ARIADNE_TEST_EQUALS(cmp(ExactDouble(1.0),two_mp),Comparison::LESS);
        ARIADNE_TEST_EQUALS(cmp(ExactDouble(2.0),two_mp),Comparison::EQUAL);
        ARIADNE_TEST_EQUALS(cmp(ExactDouble(3.0),two_mp),Comparison::GREATER);

        ARIADNE_TEST_EQUALS(cmp(Nat(1u),x),Comparison::LESS);
        ARIADNE_TEST_EQUALS(cmp(Nat(2u),x),Comparison::GREATER);
        ARIADNE_TEST_EQUALS(cmp(Nat(3u),FloatMP(3,pr)),Comparison::EQUAL);

        FloatDP dp_low(Dyadic(1),dp);
        FloatDP dp_equal(Dyadic(3,1u),dp);
        FloatDP dp_high(Dyadic(2),dp);
        ARIADNE_TEST_EQUALS(cmp(x,dp_low),Comparison::GREATER);
        ARIADNE_TEST_EQUALS(cmp(x,dp_equal),Comparison::EQUAL);
        ARIADNE_TEST_EQUALS(cmp(x,dp_high),Comparison::LESS);
    }
}


Void test_double_runtime()
{
    set_default_builtin_rounding();
    ARIADNE_TEST_EQUALS(get_builtin_rounding_mode(),ROUND_UPWARD);
    set_builtin_rounding_to_nearest();

    ARIADNE_TEST_ASSERT(texp(1.0)>2.7 && texp(1.0)<2.8);
    ARIADNE_TEST_EQUALS(add_opp(1.0,2.0),3.0);
    ARIADNE_TEST_EQUALS(neg_rec_rnd(2.0),-0.5);
    ARIADNE_TEST_EQUALS(neg_rec_opp(2.0),-0.5);

    ARIADNE_TEST_EQUALS(pow_rnd(-2.0,Int(-2)),0.25);
    ARIADNE_TEST_EQUALS(pow_rnd(-2.0,Nat(2u)),4.0);
    ARIADNE_TEST_FAIL(pow_rnd(0.0,Int(-1)));

    const double inf=std::numeric_limits<double>::infinity();
    const double nan=std::numeric_limits<double>::quiet_NaN();
    ARIADNE_TEST_ASSERT(std::isnan(sqrt_rnd(nan)));
    ARIADNE_TEST_ASSERT(std::isinf(sqrt_rnd(inf)));
    ARIADNE_TEST_EQUALS(sqrt_rnd(0.0),0.0);
    ARIADNE_TEST_FAIL(sqrt_rnd(-1.0));
    ARIADNE_TEST_ASSERT(std::isnan(exp_rnd(nan)));
    ARIADNE_TEST_ASSERT(std::isinf(exp_rnd(inf)));
    ARIADNE_TEST_EQUALS(exp_rnd(-inf),0.0);
    ARIADNE_TEST_ASSERT(std::isinf(log_rnd(0.0)));
    ARIADNE_TEST_ASSERT(std::signbit(log_rnd(0.0)));
    ARIADNE_TEST_FAIL(log_rnd(-1.0));
    ARIADNE_TEST_FAIL(atan_rnd_series(0.5));

    set_builtin_rounding_toward_zero();
    ARIADNE_TEST_EXECUTE(pi_rnd());
    ARIADNE_TEST_EXECUTE(pi_opp());
    set_builtin_rounding_to_nearest();

    const double sin_inputs[] = {-4.0,-2.8,-2.0,-1.2,-0.5,0.1,1.0,1.7,2.5,4.0};
    for(double x : sin_inputs) {
        ARIADNE_TEST_ASSERT(std::abs(sin_rnd(x)-std::sin(x))<1e-12);
    }
    ARIADNE_TEST_ASSERT(std::isnan(sin_rnd(inf)));

    ARIADNE_TEST_ASSERT(std::abs(cos_rnd(0.1)-std::cos(0.1))<1e-12);
    ARIADNE_TEST_ASSERT(std::isnan(cos_rnd(inf)));

    ARIADNE_TEST_ASSERT(std::abs(tan_rnd(-4.0)-std::tan(-4.0))<1e-11);
    ARIADNE_TEST_ASSERT(std::isnan(tan_rnd(inf)));

    set_builtin_rounding_toward_zero();
    ARIADNE_TEST_EXECUTE(atan_rnd(-0.5));
    set_builtin_rounding_to_nearest();
}

Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);
    
    test_float_literals();
    test_double_runtime();
    TestFloat<DoublePrecision>(dp).test();
    TestFloat<MultiplePrecision>(MultiplePrecision(128_bits)).test();

    return ARIADNE_TEST_FAILURES;
}

