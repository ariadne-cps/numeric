/***************************************************************************
 *            test_float_bounds.cpp
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
#include "numeric/float.decl.hpp"
#include "numeric/floatmp.hpp"

#include "numeric/float_bounds.hpp"

#include "numeric/float_lower_bound.hpp"
#include "numeric/float_upper_bound.hpp"
#include "numeric/float_error.hpp"
#include "numeric/lower_number.hpp"
#include "numeric/upper_number.hpp"

#include "ariadne-testing.hpp"
#include "test_floats.hpp"

using namespace Ariadne;
using namespace std;

namespace {

template<class PR>
constexpr bool check_float_bounds_concept()
{
    using R=RawFloat<PR>; using B=FloatBounds<PR>;
    return requires(PR pr, Nat m, Int n, R a, R b, B x, B y) {
        B(pr); B(n,pr); B(m,pr); B(a); B(a,b); B(n,m,pr);
        x=a; x=y;
        nul(x); pos(x); neg(x); hlf(x); sqr(x); rec(x);
        +x; -x; x+y; x-y; x*y; x/y;
        add(x,y); sub(x,y); mul(x,y); div(x,y); pow(x,m); pow(x,n);
        max(x,y); min(x,y); abs(x);
        sqrt(x); exp(x); log(x); sin(x); cos(x); tan(x); atan(x);
        x.precision(); x.lower_raw(); x.upper_raw();
    };
}

static_assert(check_float_bounds_concept<DoublePrecision>());
static_assert(check_float_bounds_concept<MultiplePrecision>());

} // namespace


template<class PR>
class TestFloatBounds : public Ariadne::Testing::ContextFixture
{
    typedef RawFloat<PR> RawFloatType;
    typedef FloatType<BoundedTag,PR> FloatBoundsType;
    typedef FloatType<ExactTag,PR> FloatValueType;
  private:
    PR pr;
  public:
    TestFloatBounds(Ariadne::Testing::Context& context, PR prec)
        : ContextFixture(context), pr(prec) { }
  private:
    inline FloatBoundsType make_float_bounds(ExactDouble l, ExactDouble u) { return FloatBoundsType(l,u,pr); }

  public:
    Void test_constructors();
    Void test_header_api();
    Void test_conversions();
    Void test_input();
    Void test_class();
    Void test_comparison();
    Void test_precision();
    Void test_inl_mixed_operations();
    Void test_tpl_api_and_branches();
    Void test_correct_rounded_arithmetic();
    Void test_accurate_rounded_arithmetic();
    Void test_exact_rounded_arithmetic();
    Void test_aliasing();
    Void test_monotone_functions();
    Void test_trigonometric_functions();
    Void regression_tests();
};



template<class PR> Void
TestFloatBounds<PR>::test_header_api()
{
    Decimal decimal_one(String("1.0"));
    Decimal decimal_two(String("2.0"));
    FloatBoundsType decimal_bounds(decimal_one,decimal_two,pr);
    ARIADNE_TESTING_EXPECT_SAME(decimal_bounds,FloatBoundsType(1,2,pr));

    FloatBoundsType twoexp_bounds(TwoExp(1),pr);
    ARIADNE_TESTING_EXPECT_SAME(twoexp_bounds,FloatBoundsType(2,pr));

    RawFloatType raw_one(1,pr);
    FloatBoundsType assigned(pr);
    assigned=raw_one;
    ARIADNE_TESTING_EXPECT_SAME(assigned,FloatBoundsType(1,pr));

    FloatBoundsType source(1,2,pr);

    if constexpr (Same<PR,MultiplePrecision>) {
        RawFloatType lower_mismatched(1,MultiplePrecision(64_bits));
        RawFloatType upper_mismatched(2,MultiplePrecision(128_bits));
        ARIADNE_TESTING_EXPECT_EXCEPTION(FloatBoundsType(lower_mismatched,upper_mismatched));
    }

    ValidatedNumber generic=source.generic();
    assigned=generic;
    ARIADNE_TESTING_EXPECT_SAME(assigned,source);

    FloatBoundsType created=source.create(generic);
    ARIADNE_TESTING_EXPECT_SAME(created,source);

    ARIADNE_TESTING_EXPECT_EQUAL(source.characteristics(),pr);
    ARIADNE_TESTING_EXPECT_NO_THROW(source.generic());

    ARIADNE_TESTING_EXPECT_SAME(sqr(FloatBoundsType(-2,-1,pr)),FloatBoundsType(1,4,pr));
    ARIADNE_TESTING_EXPECT_SAME(sqr(FloatBoundsType(-2,1,pr)),FloatBoundsType(0,4,pr));

    ARIADNE_TESTING_PROVE(is_indeterminate(eq(FloatBoundsType(0,2,pr),FloatBoundsType(1,3,pr))));

    RawFloatType raw_zero(0,pr);
    ARIADNE_TESTING_EXPECT_SAME(asin(raw_zero),FloatBoundsType(0,pr));
    ARIADNE_TESTING_EXPECT_SAME(acos(raw_one),FloatBoundsType(0,pr));

    FloatBoundsType point(raw_one);
    ARIADNE_TESTING_PROVE(not possibly(point!=raw_one));
    ARIADNE_TESTING_PROVE(definitely(point<=raw_one));
    ARIADNE_TESTING_PROVE(definitely(point>=raw_one));

    ARIADNE_TESTING_PROVE(definitely(raw_one==point));
    ARIADNE_TESTING_PROVE(not possibly(raw_one!=point));
    ARIADNE_TESTING_PROVE(not possibly(raw_one<point));
    ARIADNE_TESTING_PROVE(not possibly(raw_one>point));
    ARIADNE_TESTING_PROVE(definitely(raw_one<=point));
    ARIADNE_TESTING_PROVE(definitely(raw_one>=point));
}

template<class PR> Void
TestFloatBounds<PR>::test_precision()
{

    FloatBoundsType bx(Rational(1),pr);
    ARIADNE_TESTING_EXPECT_EQUAL(bx.lower().precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL(bx.upper().precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL(max(bx,bx).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL(min(bx,bx).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL(abs(bx).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL((bx+bx).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL((bx-bx).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL((bx*bx).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL((bx/bx).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL(sqrt(bx).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL(exp(bx).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL(log(bx).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL(sin(bx).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL(cos(bx).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL(tan(bx).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL(atan(bx).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL((bx+2u).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL((bx-2u).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL((bx*2u).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL((bx/2u).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL((bx+2).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL((bx-2).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL((bx*2).precision(),pr);
    ARIADNE_TESTING_EXPECT_EQUAL((bx/2).precision(),pr);

    ARIADNE_TESTING_EXPECT_EQUAL(abs(FloatBoundsType(-1,2,pr)).lower().precision(),pr);

}

template<class PR> Void
TestFloatBounds<PR>::test_inl_mixed_operations()
{
    FloatBoundsType x(-2,3,pr);
    RawFloatType positive(2,pr);
    RawFloatType negative(-2,pr);
    RawFloatType zero(0,pr);

    ARIADNE_TESTING_EXPECT_SAME(
        Operations<FloatBoundsType>::_add(x,positive),
        FloatBoundsType(0,5,pr));
    ARIADNE_TESTING_EXPECT_SAME(
        Operations<FloatBoundsType>::_add(positive,x),
        FloatBoundsType(0,5,pr));

    ARIADNE_TESTING_EXPECT_SAME(
        Operations<FloatBoundsType>::_sub(x,positive),
        FloatBoundsType(-4,1,pr));
    ARIADNE_TESTING_EXPECT_SAME(
        Operations<FloatBoundsType>::_sub(positive,x),
        FloatBoundsType(-1,4,pr));

    ARIADNE_TESTING_EXPECT_SAME(
        Operations<FloatBoundsType>::_mul(x,positive),
        FloatBoundsType(-4,6,pr));
    ARIADNE_TESTING_EXPECT_SAME(
        Operations<FloatBoundsType>::_mul(x,negative),
        FloatBoundsType(-6,4,pr));
    ARIADNE_TESTING_EXPECT_SAME(
        Operations<FloatBoundsType>::_mul(positive,x),
        FloatBoundsType(-4,6,pr));
    ARIADNE_TESTING_EXPECT_SAME(
        Operations<FloatBoundsType>::_mul(negative,x),
        FloatBoundsType(-6,4,pr));

    ARIADNE_TESTING_EXPECT_SAME(
        Operations<FloatBoundsType>::_div(x,positive),
        FloatBoundsType(-1,1.5_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(
        Operations<FloatBoundsType>::_div(x,negative),
        FloatBoundsType(-1.5_x,1,pr));

    auto divided_by_zero=Operations<FloatBoundsType>::_div(x,zero);
    ARIADNE_TESTING_PROVE(is_inf(divided_by_zero.lower_raw()));
    ARIADNE_TESTING_PROVE(is_inf(divided_by_zero.upper_raw()));
    ARIADNE_TESTING_PROVE(divided_by_zero.lower_raw()<zero);
    ARIADNE_TESTING_PROVE(divided_by_zero.upper_raw()>zero);

    FloatBoundsType positive_denominator(1,4,pr);
    FloatBoundsType negative_denominator(-4,-1,pr);
    FloatBoundsType crossing_zero(-1,1,pr);

    ARIADNE_TESTING_EXPECT_SAME(
        Operations<FloatBoundsType>::_div(positive,positive_denominator),
        FloatBoundsType(0.5_x,2,pr));
    ARIADNE_TESTING_EXPECT_SAME(
        Operations<FloatBoundsType>::_div(negative,positive_denominator),
        FloatBoundsType(-2,-0.5_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(
        Operations<FloatBoundsType>::_div(positive,negative_denominator),
        FloatBoundsType(-2,-0.5_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(
        Operations<FloatBoundsType>::_div(negative,negative_denominator),
        FloatBoundsType(0.5_x,2,pr));

    auto divided_by_crossing=Operations<FloatBoundsType>::_div(positive,crossing_zero);
    ARIADNE_TESTING_PROVE(is_inf(divided_by_crossing.lower_raw()));
    ARIADNE_TESTING_PROVE(is_inf(divided_by_crossing.upper_raw()));
    ARIADNE_TESTING_PROVE(divided_by_crossing.lower_raw()<zero);
    ARIADNE_TESTING_PROVE(divided_by_crossing.upper_raw()>zero);
}

template<class PR> Void
TestFloatBounds<PR>::test_tpl_api_and_branches()
{
    FloatLowerBound<PR> concrete_lower(-1,pr);
    FloatUpperBound<PR> concrete_upper(2,pr);
    ValidatedLowerNumber generic_lower(concrete_lower);
    ValidatedUpperNumber generic_upper(concrete_upper);

    FloatBoundsType mixed_upper(concrete_lower,generic_upper);
    FloatBoundsType mixed_lower(generic_lower,concrete_upper);
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(mixed_upper.lower_raw(),RawFloatType(-1,pr));
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(mixed_upper.upper_raw(),RawFloatType(2,pr));
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(mixed_lower.lower_raw(),RawFloatType(-1,pr));
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(mixed_lower.upper_raw(),RawFloatType(2,pr));

    FloatBoundsType base(-1,2,pr);
    FloatError<PR> error(1u,pr);
    FloatBoundsType widened=base.pm(error);
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(widened.lower_raw(),RawFloatType(-2,pr));
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(widened.upper_raw(),RawFloatType(3,pr));

    ARIADNE_TESTING_EXPECT_NO_THROW(trunc(FloatBoundsType(1,1,pr)));
    Dyadic delta(1,30u);
    FloatBoundsType inexact(Dyadic(1)-delta,Dyadic(1)+delta,pr);
    FloatBoundsType truncated=trunc(inexact);
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(truncated.lower_raw(),inexact.lower_raw());
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(truncated.upper_raw(),inexact.upper_raw());
    ARIADNE_TESTING_EXPECT_NO_THROW(trunc(inexact,Nat(40u)));

    // All supported floating-point bounds use the same public
    // zero-interval representation, including native FloatMP.
    {
        std::ostringstream zero_output;
        zero_output << FloatBoundsType(0,pr);
        ARIADNE_TESTING_EXPECT_EQUAL(zero_output.str(),String("0.0[:]"));
    }

    // The endpoints meet, but the two intervals are not equal singleton
    // intervals: the second equality in eq() must evaluate to false.
    {
        FloatBoundsType left(1,2,pr);
        FloatBoundsType right(0,1,pr);
        ARIADNE_TESTING_PROVE(possibly(eq(left,right)));
        ARIADNE_TESTING_PROVE(not definitely(eq(left,right)));
    }

    ARIADNE_TESTING_EXPECT_SAME(cos(FloatBoundsType(-10,10,pr)),FloatBoundsType(-1,1,pr));
    ARIADNE_TESTING_EXPECT_NO_THROW(cos(FloatBoundsType(-3.2_q,0.2_q,pr)));
    ARIADNE_TESTING_EXPECT_NO_THROW(cos(FloatBoundsType(-0.2_q,3.2_q,pr)));
    ARIADNE_TESTING_EXPECT_NO_THROW(cos(FloatBoundsType(1.0_q,4.0_q,pr)));

    FloatBoundsType empty(2,1,pr);
    ARIADNE_TESTING_EXPECT_EXCEPTION(cos(empty));
    ARIADNE_TESTING_EXPECT_EXCEPTION(tan(empty));
}

template<class PR> Void
TestFloatBounds<PR>::test_conversions()
{
    ARIADNE_TESTING_EXPECT_EQUAL(cast_integer(FloatBounds<PR>(Dyadic(2),pr)),Integer(2));
    ARIADNE_TESTING_EXPECT_EQUAL(cast_integer(FloatBounds<PR>(Dyadic(1),Dyadic(2),pr)),Integer(2));
    ARIADNE_TESTING_EXPECT_EQUAL(cast_integer(FloatBounds<PR>(Dyadic(5,2u),Dyadic(2),pr)),Integer(2));
    ARIADNE_TESTING_EXPECT_EQUAL(cast_integer(FloatBounds<PR>(Dyadic(7,2u),Dyadic(9,2u),pr)),Integer(2));
    ARIADNE_TESTING_EXPECT_EXCEPTION(cast_integer(FloatBounds<PR>(Dyadic(9,2u),Dyadic(11,2u),pr)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(cast_integer(FloatBoundsType(-0.75_x,-0.625_x,pr)));
}



// Test that interval arithmetic is rounded correctly,
// without paying attention to accuracy issues.
template<class PR> Void
TestFloatBounds<PR>::test_correct_rounded_arithmetic()
{
    FloatBoundsType onethird=FloatBoundsType(1,pr)/FloatBoundsType(3,pr);
    ARIADNE_TESTING_EXPECT_LESS(onethird.lower_raw(),onethird.upper_raw());
    FloatBoundsType one_approx=onethird*FloatBoundsType(3,pr);
    ARIADNE_TESTING_EXPECT_LESS(one_approx.lower_raw(),1.0_x);
    ARIADNE_TESTING_EXPECT_GREATER(one_approx.upper_raw(),1.0_x);
}


// Test that interval arithmetic gives the most accurate rounded values
template<class PR> Void
TestFloatBounds<PR>::test_accurate_rounded_arithmetic()
{
    if constexpr (Same<PR,DoublePrecision>) {
        ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(1,pr)/FloatBoundsType(3,pr),FloatBoundsType(0.33333333333333331_pr,0.33333333333333337_pr,pr));
        ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(2,pr)/FloatBoundsType(5,pr),FloatBoundsType(0.39999999999999997_pr,0.40000000000000002_pr,pr));

        ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(3,5,pr),-1),FloatBoundsType(0.19999999999999998_pr,0.33333333333333337_pr,pr));
        ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(3,5,pr),-2),FloatBoundsType(0.039999999999999986955_pr,0.11111111111111114658_pr,pr));
    }
}


// Test that interval arithmetic gives exact values if possible
template<class PR> Void
TestFloatBounds<PR>::test_exact_rounded_arithmetic()
{
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(5,7,pr)+FloatBoundsType(2,4,pr),FloatBoundsType(7,11,pr));
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(5,7,pr)-FloatBoundsType(2,6,pr),FloatBoundsType(-1,5,pr));

    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(5,7,pr)*FloatBoundsType(2,4,pr),FloatBoundsType(10,28,pr));
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(5,7,pr)*FloatBoundsType(-2,4,pr),FloatBoundsType(-14,28,pr));
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(5,7,pr)*FloatBoundsType(-4,-2,pr),FloatBoundsType(-28,-10,pr));
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(-7,5,pr)*FloatBoundsType(2,4,pr),FloatBoundsType(-28,20,pr));
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(-7,5,pr)*FloatBoundsType(-2,4,pr),FloatBoundsType(-28,20,pr));
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(-7,5,pr)*FloatBoundsType(-4,-2,pr),FloatBoundsType(-20,28,pr));
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(-7,-5,pr)*FloatBoundsType(2,4,pr),FloatBoundsType(-28,-10,pr));
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(-7,-5,pr)*FloatBoundsType(-2,4,pr),FloatBoundsType(-28,14,pr));
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(-7,-5,pr)*FloatBoundsType(-4,-2,pr),FloatBoundsType(10,28,pr));

    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(5,7,pr)/FloatBoundsType(2,4,pr),FloatBoundsType(1.25_x,3.50_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(5,7,pr)/FloatBoundsType(-4,-2,pr),FloatBoundsType(-3.50_x,-1.25_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(-7,5,pr)/FloatBoundsType(2,4,pr),FloatBoundsType(-3.50_x,2.50_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(-7,5,pr)/FloatBoundsType(-4,-2,pr),FloatBoundsType(-2.50_x,3.5_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(-7,-5,pr)/FloatBoundsType(2,4,pr),FloatBoundsType(-3.50_x,-1.25_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(FloatBoundsType(-7,-5,pr)/FloatBoundsType(-4,-2,pr),FloatBoundsType(1.25_x,3.50_x,pr));

    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(5,7,pr),0u),FloatBoundsType(1,1,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-5,7,pr),0u),FloatBoundsType(1,1,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,5,pr),0u),FloatBoundsType(1,1,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,-5,pr),0u),FloatBoundsType(1,1,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(5,7,pr),1u),FloatBoundsType(5,7,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-5,7,pr),1u),FloatBoundsType(-5,7,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,5,pr),1u),FloatBoundsType(-7,5,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,-5,pr),1u),FloatBoundsType(-7,-5,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(5,7,pr),2u),FloatBoundsType(25,49,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-5,7,pr),2u),FloatBoundsType(0,49,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,5,pr),2u),FloatBoundsType(0,49,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,-5,pr),2u),FloatBoundsType(25,49,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(5,7,pr),3u),FloatBoundsType(125,343,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-5,7,pr),3u),FloatBoundsType(-125,343,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,5,pr),3u),FloatBoundsType(-343,125,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,-5,pr),3u),FloatBoundsType(-343,-125,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(5,7,pr),4u),FloatBoundsType(625,2401,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-5,7,pr),4u),FloatBoundsType(0,2401,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,5,pr),4u),FloatBoundsType(0,2401,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,-5,pr),4u),FloatBoundsType(625,2401,pr));

    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(5,7,pr),0),FloatBoundsType(1,1,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,5,pr),0),FloatBoundsType(1,1,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,-5,pr),0),FloatBoundsType(1,1,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(5,7,pr),1),FloatBoundsType(5,7,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-5,7,pr),1),FloatBoundsType(-5,7,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,5,pr),1),FloatBoundsType(-7,5,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,-5,pr),1),FloatBoundsType(-7,-5,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(5,7,pr),2),FloatBoundsType(25,49,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-5,7,pr),2),FloatBoundsType(0,49,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,5,pr),2),FloatBoundsType(0,49,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,-5,pr),2),FloatBoundsType(25,49,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(5,7,pr),3),FloatBoundsType(125,343,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-5,7,pr),3),FloatBoundsType(-125,343,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,5,pr),3),FloatBoundsType(-343,125,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,-5,pr),3),FloatBoundsType(-343,-125,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(5,7,pr),4),FloatBoundsType(625,2401,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-5,7,pr),4),FloatBoundsType(0,2401,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,5,pr),4),FloatBoundsType(0,2401,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,-5,pr),4),FloatBoundsType(625,2401,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(5,7,pr),5),FloatBoundsType(3125,16807,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-5,7,pr),5),FloatBoundsType(-3125,16807,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,5,pr),5),FloatBoundsType(-16807,3125,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,-5,pr),5),FloatBoundsType(-16807,-3125,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(5,7,pr),7),FloatBoundsType(78125,823543,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-5,7,pr),7),FloatBoundsType(-78125,823543,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,5,pr),7),FloatBoundsType(-823543,78125,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-7,-5,pr),7),FloatBoundsType(-823543,-78125,pr));

    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(2,4,pr),-1),FloatBoundsType(0.25_x,0.5_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-4,-2,pr),-1),FloatBoundsType(-0.5_x,-0.25_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(2,4,pr),-2),FloatBoundsType(0.0625_x,0.25_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-4,-2,pr),-2),FloatBoundsType(0.0625_x,0.25_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(2,4,pr),-3),FloatBoundsType(0.015625_x,0.125_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(pow(FloatBoundsType(-4,-2,pr),-3),FloatBoundsType(-0.125_x,-0.015625_x,pr));

    ARIADNE_TESTING_EXPECT_SAME(rec(FloatBoundsType(2,4,pr)),FloatBoundsType(0.25_x,0.50_x,pr));
    ARIADNE_TESTING_EXPECT_SAME(rec(FloatBoundsType(-4,-2,pr)),FloatBoundsType(-0.50_x,-0.25_x,pr));
}



template<> Void
TestFloatBounds<DoublePrecision>::test_constructors()
{
    FloatDP zero(0,pr);

    // Construct from pair
    FloatBoundsType xd1(FloatDP(1.125_x,pr),FloatDP(2.25_x,pr));
    ARIADNE_TESTING_PROVE(xd1.lower_raw()==1.125_x); ARIADNE_TESTING_PROVE(xd1.upper_raw()==2.25_x);

    // Default constructor
    FloatBoundsType xd2(pr);
    if(xd2.lower_raw()>xd2.upper_raw()) {
        { std::ostringstream warning_stream; warning_stream << "FloatBoundsType default constructor returns an empty set."; this->ctx().warning(warning_stream.str()); }
    } else {
        ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(same,xd2,FloatBoundsType(zero,zero));
    }

    // Constructor from Rational with approximations
    FloatBoundsType xd3(Rational(21,10),pr);
    ARIADNE_TESTING_EXPECT_LESS(Rational(xd3.lower_raw()),Rational(21,10));
    ARIADNE_TESTING_EXPECT_GREATER(Rational(xd3.upper_raw()),Rational(21,10));

    // Constructor from rational bounds
    FloatBoundsType xd4(2.1_q,3.2_q,pr);
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(xd4.lower_raw(),2.1_pr);
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(xd4.upper_raw(),3.2_pr);

    // Approximate constructor from a single value
    FloatBoundsType xd5(Rational(1,3),pr);
    ARIADNE_TESTING_EXPECT_LESS(Rational(xd5.lower_raw()),Rational(1,3));
    ARIADNE_TESTING_EXPECT_GREATER(Rational(xd5.upper_raw()),Rational(1,3));

    // ExactTag constructor from a single value
    FloatBoundsType xd6(FloatDP(1.25_x,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(xd6.lower_raw(),FloatDP(1.25_x,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(xd6.upper_raw(),FloatDP(1.25_x,pr));
}

template<class PR> Void
TestFloatBounds<PR>::test_constructors()
{
    RawFloatType zero(0,pr);

    // Construct from pair
    FloatBoundsType xd1(RawFloatType(1.125_x,pr),RawFloatType(2.25_x,pr));
    ARIADNE_TESTING_PROVE(xd1.lower_raw()==1.125_x); ARIADNE_TESTING_PROVE(xd1.upper_raw()==2.25_x);

    // Default value constructor
    FloatBoundsType xd2(pr);
    if(xd2.lower_raw()>xd2.upper_raw()) {
        { std::ostringstream warning_stream; warning_stream << "FloatBoundsType default constructor returns an empty set."; this->ctx().warning(warning_stream.str()); }
    } else {
        ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(same,xd2,FloatBoundsType(zero,zero));
    }

    // Constructor with approximations
    FloatBoundsType xd3(Rational(21,10),pr);
    ARIADNE_TESTING_EXPECT_LESS(Rational(xd3.lower_raw()),Rational(21,10));
    ARIADNE_TESTING_EXPECT_GREATER(Rational(xd3.upper_raw()),Rational(21,10));

    // Constructor from approximate values
    FloatBoundsType xd4(2.1_q,3.2_q,pr);
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(xd4.lower_raw(),Rational(21,10));
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(xd4.upper_raw(),Rational(32,10));

    // ApproximateTag constructor from a single value
    FloatBoundsType xd5(Rational(1,3),pr);
    ARIADNE_TESTING_EXPECT_LESS(xd5.lower_raw(),Rational(1,3));
    ARIADNE_TESTING_EXPECT_GREATER(xd5.upper_raw(),Rational(1,3));

    // ExactTag constructor from a single value
    FloatBoundsType xd6(RawFloatType(1.25_x,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(xd6.lower_raw(),RawFloatType(1.25_x,pr));
    ARIADNE_TESTING_EXPECT_EQUAL(xd6.upper_raw(),RawFloatType(1.25_x,pr));
}

template<class PR> Void TestFloatBounds<PR>::test_class()
{
    // Test lower, upper, midpoint, radius, width
    RawFloatType one_(1.0_x,pr);
    RawFloatType two_(2.0_x,pr);
    RawFloatType three_(3.0_x,pr);

    // Tests for exact operations
    ARIADNE_TESTING_EXPECT_EQUAL(FloatBoundsType(-0.25_x,0.50_x,pr).lower().raw(),-0.25_x);
    ARIADNE_TESTING_EXPECT_EQUAL(FloatBoundsType(-0.25_x,0.50_x,pr).upper().raw(),0.5_x);
    ARIADNE_TESTING_EXPECT_EQUAL(FloatBoundsType(-0.25_x,0.50_x,pr).value(),0.125_x);
    ARIADNE_TESTING_EXPECT_EQUAL(FloatBoundsType(-0.25_x,0.50_x,pr).error().raw(),0.375_x);

    FloatBoundsType midpoint_bounds(one_,three_);
    ARIADNE_TESTING_EXPECT_EQUAL(midpoint(midpoint_bounds),two_);
    if constexpr (Same<PR,DoublePrecision>) {
        ARIADNE_TESTING_EXPECT_EQUAL(class_name<FloatBoundsType>(),String("FloatDPBounds"));
    } else {
        ARIADNE_TESTING_EXPECT_EQUAL(class_name<FloatBoundsType>(),String("FloatMPBounds"));
    }

    // Tests for inexact operations
    ARIADNE_TESTING_EXPECT_EQUAL((FloatBoundsType(-1,2,pr)/3).lower().raw(),div(down,-one_,three_));
    ARIADNE_TESTING_EXPECT_EQUAL((FloatBoundsType(-1,2,pr)/3).upper().raw(),div(up,two_,three_));
    ARIADNE_TESTING_EXPECT_EQUAL((FloatBoundsType(-1,2,pr)/3).value(),hlf(add(up,div(down,-one_,three_),div(up,two_,three_))));
    ARIADNE_TESTING_EXPECT_EQUAL((FloatBoundsType(-1,2,pr)/3).error().raw(),hlf(sub(up,div(up,two_,three_),div(down,-one_,three_))));
}

template<class PR> Void TestFloatBounds<PR>::test_input()
{
    FloatBoundsType x(pr);

    string input("[1.125,2.75] [0.4,0.6]");
    stringstream iss(input);

    iss >> x;
    ARIADNE_TESTING_EXPECT_EQUAL(x.lower_raw(),Rational(9,8));
    ARIADNE_TESTING_EXPECT_EQUAL(x.upper_raw(),Rational(11,4));

    iss >> x;
    // TODO: Make stream input produce rounded values
    //ARIADNE_TEST_COMPARE(x.lower_raw(),<,Rational(2,5))
    // ARIADNE_TEST_COMPARE(x.upper_raw(),>,Rational(3,5))
    if(not(x.lower_raw()<=Rational(2,5) and x.upper_raw()>=Rational(3,5))) {
        { std::ostringstream warning_stream; warning_stream << "FloatType<BoundedTag,"<<class_name<PR>()<<"> string constructor returns an approximate interval, not an outwardly rounded interval."; this->ctx().warning(warning_stream.str()); }
    }

    {
        stringstream alternate("(1.0:2.0) [1.0;2.0]");
        alternate >> x;
        ARIADNE_TESTING_PROVE(not alternate.fail());
        alternate >> x;
        ARIADNE_TESTING_PROVE(not alternate.fail());
    }

    {
        stringstream invalid_left("{1.0:2.0]");
        ARIADNE_TESTING_EXPECT_EXCEPTION(invalid_left >> x);
    }
    {
        stringstream invalid_middle("[1.0|2.0]");
        ARIADNE_TESTING_EXPECT_EXCEPTION(invalid_middle >> x);
    }
    {
        stringstream invalid_right("[1.0:2.0}");
        ARIADNE_TESTING_EXPECT_EXCEPTION(invalid_right >> x);
    }
    {
        stringstream invalid_number("[x:2.0]");
        ARIADNE_TESTING_EXPECT_EXCEPTION(invalid_number >> x);
    }
}

template<class PR> Void TestFloatBounds<PR>::test_comparison() {
    // FIXME: If using Boost style interval tests, uncomment the line below
    // and comment out the line after
    //ARIADNE_TESTING_PROVE(indeterminate(ivld1==ivld2));
    FloatBoundsType ivl1(1.125_x,2.25_x,pr);
    FloatBoundsType ivl2=ivl1;

    ARIADNE_TESTING_PROVE(!definitely(ivl1==ivl2));
    ARIADNE_TESTING_PROVE(possibly(ivl1==ivl2));
    FloatBoundsType& ivl1ref=ivl1;
    ivl1ref=FloatBoundsType(5.25_x,7.375_x,pr);
    cout << "ivl1ref=" << ivl1ref << endl;
    ARIADNE_TESTING_PROVE(ivl1ref.lower_raw()==5.25_x);
}

template<class PR> Void TestFloatBounds<PR>::test_aliasing() {
    FloatValueType ex2(1.5_x,pr);
    FloatValueType ex3(2.25_x,pr);

    FloatBoundsType vx1(pr);
    FloatBoundsType vx2(1.5_x,2.25_x,pr);
    FloatBoundsType vx3(3.125_x,4.0625_x,pr);

    // Check to make sure aliases are handled correctly
    vx1=vx3; vx1=vx2-vx1; ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(same,vx1,FloatBoundsType(vx2-vx3));
    vx1=vx3; vx1=vx2*vx1; ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(same,vx1,FloatBoundsType(vx2*vx3));
    vx1=vx2; vx1=vx1*vx3; ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(same,vx1,FloatBoundsType(vx2*vx3));
    vx1=vx2; vx1=vx1*ex3; ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(same,vx1,FloatBoundsType(vx2*ex3));
    vx1=vx3; vx1=ex2*vx1; ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(same,vx1,FloatBoundsType(ex2*vx3));
    vx1=vx2; vx1=vx1/vx3; ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(same,vx1,FloatBoundsType(vx2/vx3));
    vx1=vx2; vx1=vx1/ex3; ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(same,vx1,FloatBoundsType(vx2/ex3));
    vx1=vx3; vx1=ex2/vx1; ARIADNE_TESTING_EXPECT_BINARY_PREDICATE(same,vx1,FloatBoundsType(ex2/vx3));
}

template<class PR> Void TestFloatBounds<PR>::test_monotone_functions()
{
    FloatBoundsType two_(2.0_x,pr);
    FloatBoundsType sqrt_two=sqrt(two_);
    FloatBoundsType sqr_sqrt_two=sqr(sqrt_two);
    ARIADNE_TESTING_INSPECT("sqrt_two",sqrt_two);
    ARIADNE_TESTING_INSPECT("sqr_sqrt_two",sqr_sqrt_two);
    if (Same<PR,DoublePrecision>) {
        ARIADNE_TESTING_EXPECT_EQUAL(sqrt_two.lower_raw(),1.4142135623730949_pr);
        ARIADNE_TESTING_EXPECT_EQUAL(sqrt_two.upper_raw(),1.4142135623730951_pr);
        ARIADNE_TESTING_EXPECT_LESS_EQUAL(sqrt_two.lower_raw(),1.4142135623730949_pr);
        ARIADNE_TESTING_EXPECT_GREATER(sqrt_two.lower_raw(),1.4142135623730947_pr);
        ARIADNE_TESTING_EXPECT_GREATER_EQUAL(sqrt_two.upper_raw(),1.4142135623730951_pr);
        ARIADNE_TESTING_EXPECT_LESS(sqrt_two.upper_raw(),1.4142135623730954_pr);
    }
    ARIADNE_TESTING_EXPECT_EQUAL(sqrt_two.lower_raw(),sqrt(down,RawFloatType(2.0_x,pr)));
    ARIADNE_TESTING_EXPECT_EQUAL(sqrt_two.upper_raw(),sqrt(up,RawFloatType(2.0_x,pr)));
    ARIADNE_TESTING_PROVE(sqrt_two.lower_raw()<sqrt_two.upper_raw());
    ARIADNE_TESTING_EXPECT_LESS(sqr_sqrt_two.lower_raw(),2.0_x);
    ARIADNE_TESTING_EXPECT_GREATER(sqr_sqrt_two.upper_raw(),2.0_x);

    FloatBoundsType one(1.0_x,pr);
    FloatBoundsType exp_one=exp(one);
    FloatBoundsType log_exp_one=log(exp_one);
    ARIADNE_TESTING_INSPECT("exp_one",exp_one);
    ARIADNE_TESTING_INSPECT("log_exp_one",log_exp_one);
    if (Same<PR,DoublePrecision>) {
        ARIADNE_TESTING_EXPECT_LESS(exp_one.lower_raw(),2.71828182845905_pr);
        ARIADNE_TESTING_EXPECT_GREATER(exp_one.lower_raw(),2.71828182845903_pr);
        ARIADNE_TESTING_EXPECT_GREATER(exp_one.upper_raw(),2.71828182845904_pr);
        ARIADNE_TESTING_EXPECT_LESS(exp_one.upper_raw(),2.71828182845906_pr);
    }
    ARIADNE_TESTING_EXPECT_EQUAL(exp_one.lower_raw(),exp(down,RawFloatType(1.0_x,pr)));
    ARIADNE_TESTING_EXPECT_EQUAL(exp_one.upper_raw(),exp(up,RawFloatType(1.0_x,pr)));
    ARIADNE_TESTING_PROVE(exp_one.lower_raw()<exp_one.upper_raw());
    ARIADNE_TESTING_EXPECT_LESS(log_exp_one.lower_raw(),1.0_x);
    ARIADNE_TESTING_EXPECT_LESS(log_exp_one.lower_raw(),1.0_x);

    if (Same<PR,DoublePrecision>) {
        FloatBoundsType e(2.7182818284590451_pr,2.7182818284590455_pr,pr);
        FloatBoundsType loge=log(e);
        ARIADNE_TESTING_INSPECT("e",e);
        ARIADNE_TESTING_EXPECT_LESS(loge.lower_raw(),1.0_x);
        ARIADNE_TESTING_EXPECT_GREATER(loge.lower_raw(),0.9999999999998_pr);
        ARIADNE_TESTING_EXPECT_GREATER(loge.upper_raw(),1.0_x);
        ARIADNE_TESTING_EXPECT_LESS(loge.upper_raw(),1.000000000002_pr);
    }

}

template<class PR> Void TestFloatBounds<PR>::test_trigonometric_functions()
{
    try {
        FloatBoundsType x(6.283185307179586_pr,6.283185307179587_pr,pr);
        FloatBoundsType sinx=sin(x);
        ARIADNE_TESTING_INSPECT("x",x);
        ARIADNE_TESTING_EXPECT_LESS(sinx.lower_raw(),0.0_x);
        ARIADNE_TESTING_EXPECT_GREATER(sinx.lower_raw(),-1e-14_pr);
        ARIADNE_TESTING_EXPECT_GREATER(sinx.upper_raw(),0.0_x);
        ARIADNE_TESTING_EXPECT_LESS(sinx.upper_raw(),+1e-14_pr);
        ARIADNE_TESTING_PROVE(sinx.lower_raw()<sinx.upper_raw());
    }
    catch(...) { }

    try {
        FloatBoundsType x(7.0685834705770345_pr,pr);
        FloatBoundsType sinx=sin(x);
        ARIADNE_TESTING_INSPECT("x",x);
        ARIADNE_TESTING_EXPECT_LESS(sinx.lower_raw(),0.7071067811866_pr);
        ARIADNE_TESTING_EXPECT_GREATER(sinx.upper_raw(),0.7071067811865_pr);
        ARIADNE_TESTING_PROVE(sinx.lower_raw()<sinx.upper_raw());
    }
    catch(...) { }

}

template<class PR> Void TestFloatBounds<PR>::regression_tests() {
    RawFloatType inf_=RawFloatType::inf(pr);

    // Regression test; fails dramatically on certain types of rounding
    {
        FloatBoundsType x(1.5707963267948966_pr,1.5707963267948968_pr,pr);
        FloatBoundsType cosx=cos(x);
        ARIADNE_TESTING_INSPECT("x",x);
        ARIADNE_TESTING_EXPECT_LESS(cosx.lower_raw(),0.0_x);
        ARIADNE_TESTING_EXPECT_GREATER(cosx.lower_raw(),-1e-14_pr);
        ARIADNE_TESTING_EXPECT_GREATER(cosx.upper_raw(),0.0_x);
        ARIADNE_TESTING_EXPECT_LESS(cosx.upper_raw(),+1e-14_pr);
        ARIADNE_TESTING_PROVE(cosx.lower_raw()<cosx.upper_raw());
    }

    // Regression test for inverse trigonometric functions at the domain endpoints.
    // These used to reach sqrt_rnd(-inf) through the asin transformation.
    {
        FloatBoundsType minus_one(-1,pr);
        FloatBoundsType plus_one(1,pr);

        FloatBoundsType asin_minus_one=asin(minus_one);
        FloatBoundsType asin_plus_one=asin(plus_one);
        FloatBoundsType acos_minus_one=acos(minus_one);
        FloatBoundsType acos_plus_one=acos(plus_one);

        ARIADNE_TESTING_EXPECT_LESS(asin_minus_one.lower_raw(),-1.5_pr);
        ARIADNE_TESTING_EXPECT_GREATER(asin_minus_one.upper_raw(),-1.6_pr);
        ARIADNE_TESTING_EXPECT_LESS(asin_plus_one.lower_raw(),1.6_pr);
        ARIADNE_TESTING_EXPECT_GREATER(asin_plus_one.upper_raw(),1.5_pr);
        ARIADNE_TESTING_EXPECT_LESS(acos_minus_one.lower_raw(),3.2_pr);
        ARIADNE_TESTING_EXPECT_GREATER(acos_minus_one.upper_raw(),3.1_pr);
        ARIADNE_TESTING_EXPECT_EQUAL(acos_plus_one.lower_raw(),RawFloatType(0,pr));
        ARIADNE_TESTING_EXPECT_EQUAL(acos_plus_one.upper_raw(),RawFloatType(0,pr));
    }

    // Inverse trigonometric interval evaluation intersects with the real domain [-1,1].
    {
        FloatBoundsType below(-2,-2,pr);
        FloatBoundsType above(2,2,pr);
        FloatBoundsType crossing_low(-2,0,pr);
        FloatBoundsType crossing_high(0,2,pr);

        FloatBoundsType asin_below=asin(below);
        FloatBoundsType asin_above=asin(above);
        FloatBoundsType acos_below=acos(below);
        FloatBoundsType acos_above=acos(above);

        ARIADNE_TESTING_PROVE(asin_below.lower_raw()>asin_below.upper_raw());
        ARIADNE_TESTING_PROVE(asin_above.lower_raw()>asin_above.upper_raw());
        ARIADNE_TESTING_PROVE(acos_below.lower_raw()>acos_below.upper_raw());
        ARIADNE_TESTING_PROVE(acos_above.lower_raw()>acos_above.upper_raw());

        FloatBoundsType asin_crossing_low=asin(crossing_low);
        FloatBoundsType asin_crossing_high=asin(crossing_high);
        ARIADNE_TESTING_PROVE(asin_crossing_low.lower_raw()<=asin_crossing_low.upper_raw());
        ARIADNE_TESTING_PROVE(asin_crossing_high.lower_raw()<=asin_crossing_high.upper_raw());
        ARIADNE_TESTING_EXPECT_GREATER_EQUAL(asin_crossing_low.upper_raw(),RawFloatType(0,pr));
        ARIADNE_TESTING_EXPECT_LESS_EQUAL(asin_crossing_high.lower_raw(),RawFloatType(0,pr));
    }

    // Tangent intervals that may contain a pole must return the whole-real hull.
    {
        FloatBoundsType around_half_pi(1.5707_pr,1.5709_pr,pr);
        FloatBoundsType tan_around_half_pi=tan(around_half_pi);
        ARIADNE_TESTING_EXPECT_EQUAL(tan_around_half_pi.lower_raw(),-inf_);
        ARIADNE_TESTING_EXPECT_EQUAL(tan_around_half_pi.upper_raw(),+inf_);
    }

    // Tangent range reduction must remain valid at both pole enclosures,
    // including after translation by one period.
    {
        FloatBoundsType around_positive_half_pi(1.5707_pr,1.5709_pr,pr);
        FloatBoundsType around_negative_half_pi(-1.5709_pr,-1.5707_pr,pr);
        FloatBoundsType around_three_half_pi(4.7123_pr,4.7125_pr,pr);

        FloatBoundsType positive=tan(around_positive_half_pi);
        FloatBoundsType negative=tan(around_negative_half_pi);
        FloatBoundsType translated=tan(around_three_half_pi);

        ARIADNE_TESTING_EXPECT_EQUAL(positive.lower_raw(),-inf_);
        ARIADNE_TESTING_EXPECT_EQUAL(positive.upper_raw(),+inf_);
        ARIADNE_TESTING_EXPECT_EQUAL(negative.lower_raw(),-inf_);
        ARIADNE_TESTING_EXPECT_EQUAL(negative.upper_raw(),+inf_);
        ARIADNE_TESTING_EXPECT_EQUAL(translated.lower_raw(),-inf_);
        ARIADNE_TESTING_EXPECT_EQUAL(translated.upper_raw(),+inf_);
    }

    // Tangent away from poles remains finite and encloses the endpoint values.
    {
        FloatBoundsType x(0.25_x,0.5_x,pr);
        FloatBoundsType tx=tan(x);
        ARIADNE_TESTING_PROVE(tx.lower_raw()<=tx.upper_raw());
        ARIADNE_TESTING_EXPECT_LESS_EQUAL(tx.lower_raw(),tan(down,RawFloatType(0.25_x,pr)));
        ARIADNE_TESTING_EXPECT_GREATER_EQUAL(tx.upper_raw(),tan(up,RawFloatType(0.5_x,pr)));
    }

    // Regression test for dividing by interval with lower endpoint -0.0 or upper endpoint +0.0

    ARIADNE_TESTING_EXPECT_EQUAL((FloatBoundsType(1.0_x,2.0_x,pr)/FloatBoundsType(-0.0_x,1.0_x,pr)).upper_raw(),+inf_);
    ARIADNE_TESTING_EXPECT_EQUAL((FloatBoundsType(1.0_x,2.0_x,pr)/FloatBoundsType(-1.0_x,+0.0_x,pr)).lower_raw(),-inf_);

    ARIADNE_TESTING_EXPECT_EQUAL(rec(FloatBoundsType(-0.0_x,+1.0_x,pr)).upper_raw(),+inf_);
    ARIADNE_TESTING_EXPECT_EQUAL(rec(FloatBoundsType(-1.0_x,+0.0_x,pr)).lower_raw(),-inf_);

    // Reciprocal is undefined on the singleton zero, but intervals crossing
    // zero retain the whole-real hull of their two defined branches.
    {
        FloatBoundsType zero(0,pr);
        FloatBoundsType crossing(-1,1,pr);
        FloatBoundsType rec_zero=rec(zero);
        FloatBoundsType rec_crossing=rec(crossing);
        ARIADNE_TESTING_PROVE(rec_zero.lower_raw()>rec_zero.upper_raw());
        ARIADNE_TESTING_EXPECT_EQUAL(rec_crossing.lower_raw(),-inf_);
        ARIADNE_TESTING_EXPECT_EQUAL(rec_crossing.upper_raw(),+inf_);
    }
}


template<class PR>
void add_float_bounds_cases(Ariadne::Testing::Suite& suite,
                            std::string const& variant, PR precision) {
    using Fixture = TestFloatBounds<PR>;
    auto add = [&](std::string const& name, auto method) {
        suite.add(variant + " / " + name,
            [precision, method](Ariadne::Testing::Context& context) {
                Fixture fixture(context, precision);
                (fixture.*method)();
            });
    };
    add("constructors", &Fixture::test_constructors);
    add("header api", &Fixture::test_header_api);
    add("input", &Fixture::test_input);
    add("class", &Fixture::test_class);
    add("conversions", &Fixture::test_conversions);
    add("comparison", &Fixture::test_comparison);
    add("aliasing", &Fixture::test_aliasing);
    add("precision", &Fixture::test_precision);
    add("inl mixed operations", &Fixture::test_inl_mixed_operations);
    add("tpl api and branches", &Fixture::test_tpl_api_and_branches);
    add("correct rounded arithmetic", &Fixture::test_correct_rounded_arithmetic);
    add("accurate rounded arithmetic", &Fixture::test_accurate_rounded_arithmetic);
    add("exact rounded arithmetic", &Fixture::test_exact_rounded_arithmetic);
    add("monotone functions", &Fixture::test_monotone_functions);
    add("trigonometric functions", &Fixture::test_trigonometric_functions);
    add("regression tests", &Fixture::regression_tests);
}

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("float bounds");
    add_float_bounds_cases<DoublePrecision>(suite, "FloatDP", dp);
    add_float_bounds_cases<MultiplePrecision>(suite, "FloatMP 53", MultiplePrecision(53_bits));
    add_float_bounds_cases<MultiplePrecision>(suite, "FloatMP 128", MultiplePrecision(128_bits));
    return Ariadne::Testing::run_suite(suite, argc, argv);
}

