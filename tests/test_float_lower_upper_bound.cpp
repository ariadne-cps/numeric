/***************************************************************************
 *            test_float_lower_upper_bound.cpp
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
#include "numeric/upper_number.hpp"
#include "numeric/lower_number.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/floats.hpp"

#include "ariadne-testing.hpp"
#include "test_floats.hpp"

using namespace Ariadne;
using namespace std;

namespace {

template<class PR>
constexpr bool check_directed_float_concept()
{
    using L=FloatLowerBound<PR>;
    using U=FloatUpperBound<PR>;
    using PL=PositiveFloatLowerBound<PR>;
    using PU=PositiveFloatUpperBound<PR>;
    return requires(PR pr, Nat m, L l, L l2, U u, U u2, PL pl, PU pu) {
        L(pr); U(pr); L(0,pr); U(0,pr);
        +l; -l; +u; -u;
        l+l2; l-u; u+u2; u-l;
        l+=l2; l-=u; u+=u2; u-=l;
        pos(l); neg(u); hlf(l); add(l,l2); sub(l,u);
        pos(u); neg(l); hlf(u); add(u,u2); sub(u,l);
        sqrt(l); exp(l); log(l); atan(l); max(l,l2); min(l,l2);
        sqrt(u); exp(u); log(u); atan(u); max(u,u2); min(u,u2);
        pl+pl; pl*pl; pl/pu; pl/m;
        pu+pu; pu*pu; pu/pl; pu/m;
        l==u; l!=u; l<u; l>u; u==l; u!=l; u<l; u>l;
        l.precision(); u.precision(); l.raw(); u.raw(); l.generic(); u.generic();
    };
}

static_assert(check_directed_float_concept<DoublePrecision>());
static_assert(check_directed_float_concept<MultiplePrecision>());

} // namespace


template<class PR>
class TestDirectedFloats
    : public TestFloats<PR>, public Ariadne::Testing::ContextFixture
{
    typedef FloatType<ApproximateTag,PR> FloatApproximationType;
    typedef FloatType<LowerTag,PR> FloatLowerBoundType;
    typedef FloatType<UpperTag,PR> FloatUpperBoundType;
    typedef FloatType<BoundedTag,PR> FloatBoundsType;
    typedef FloatType<MetricTag,PR> FloatBallType;
    typedef FloatType<ExactTag,PR> FloatValueType;

    typedef PositiveFloatLowerBound<PR> PositiveFloatLowerBoundType;
    typedef PositiveFloatUpperBound<PR> PositiveFloatUpperBoundType;

  private:
    PR precision;
  public:
    TestDirectedFloats(Ariadne::Testing::Context& context, PR prec)
        : ContextFixture(context), precision(prec) { };
  public:
    Void test_precision();
    Void test_header_api();
    Void test_conversions();
    Void test_validation();
    Void test_rounded_arithmetic();
    Void test_comparison();
};

template<class PR> Void
TestDirectedFloats<PR>::test_precision()
{
    FloatLowerBoundType lx(Rational(1),precision);
    ARIADNE_TESTING_EXPECT_EQUAL(max(lx,lx).precision(),precision);
    ARIADNE_TESTING_EXPECT_EQUAL(min(lx,lx).precision(),precision);
    ARIADNE_TESTING_EXPECT_EQUAL((lx+lx).precision(),precision);
    ARIADNE_TESTING_EXPECT_EQUAL(exp(lx).precision(),precision);
    ARIADNE_TESTING_EXPECT_EQUAL((lx+2u).precision(),precision);
    ARIADNE_TESTING_EXPECT_EQUAL((lx*2u).precision(),precision);
    ARIADNE_TESTING_EXPECT_EQUAL((lx/2u).precision(),precision);
    ARIADNE_TESTING_EXPECT_EQUAL((lx+2).precision(),precision);
    ARIADNE_TESTING_EXPECT_EQUAL((lx-2).precision(),precision);
}

template<class PR> Void
TestDirectedFloats<PR>::test_header_api()
{
    using F=Float<PR>;
    using L=FloatLowerBoundType;
    using U=FloatUpperBoundType;

    L lower_zero(precision);
    U upper_zero(precision);
    ARIADNE_TESTING_EXPECT_EQUAL(lower_zero.raw(),F(0,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(upper_zero.raw(),F(0,precision));

    L lower_twoexp(TwoExp(1),precision);
    U upper_twoexp(TwoExp(1),precision);
    ARIADNE_TESTING_EXPECT_EQUAL(lower_twoexp.raw(),F(2,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(upper_twoexp.raw(),F(2,precision));

    L lower_integer(Integer(2),precision);
    U upper_integer(Integer(2),precision);
    ARIADNE_TESTING_EXPECT_EQUAL(lower_integer.raw(),F(2,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(upper_integer.raw(),F(2,precision));

    Decimal decimal(String("1.5"));
    L lower_decimal(decimal,precision);
    U upper_decimal(decimal,precision);
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(lower_decimal.raw(),Rational(3,2));
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(upper_decimal.raw(),Rational(3,2));

    F raw_two(2,precision);
    L lower_from_raw(raw_two,precision);
    U upper_from_raw(raw_two,precision);
    ARIADNE_TESTING_EXPECT_EQUAL(lower_from_raw.raw(),raw_two);
    ARIADNE_TESTING_EXPECT_EQUAL(upper_from_raw.raw(),raw_two);

    FloatBoundsType bounds(1,2,precision);
    L lower_from_bounds(bounds,precision);
    U upper_from_bounds(bounds,precision);
    ARIADNE_TESTING_EXPECT_EQUAL(lower_from_bounds.raw(),bounds.lower_raw());
    ARIADNE_TESTING_EXPECT_EQUAL(upper_from_bounds.raw(),bounds.upper_raw());

    L lower_assigned(precision);
    lower_assigned=raw_two;
    ARIADNE_TESTING_EXPECT_EQUAL(lower_assigned.raw(),raw_two);
    lower_assigned=bounds;
    ARIADNE_TESTING_EXPECT_EQUAL(lower_assigned.raw(),bounds.lower_raw());

    U upper_assigned(precision);
    upper_assigned=raw_two;
    ARIADNE_TESTING_EXPECT_EQUAL(upper_assigned.raw(),raw_two);
    upper_assigned=bounds;
    ARIADNE_TESTING_EXPECT_EQUAL(upper_assigned.raw(),bounds.upper_raw());

    ValidatedLowerNumber generic_lower=L(Rational(4,3),precision);
    ValidatedUpperNumber generic_upper=U(Rational(5,3),precision);

    lower_assigned=generic_lower;
    upper_assigned=generic_upper;
    ARIADNE_TESTING_EXPECT_EQUAL(lower_assigned.precision(),precision);
    ARIADNE_TESTING_EXPECT_EQUAL(upper_assigned.precision(),precision);

    L lower_created=lower_assigned.create(generic_lower);
    U upper_created=upper_assigned.create(generic_upper);
    U upper_from_lower=lower_assigned.create(generic_upper);
    L lower_from_upper=upper_assigned.create(generic_lower);
    ARIADNE_TESTING_EXPECT_EQUAL(lower_created.precision(),precision);
    ARIADNE_TESTING_EXPECT_EQUAL(upper_created.precision(),precision);
    ARIADNE_TESTING_EXPECT_EQUAL(upper_from_lower.precision(),precision);
    ARIADNE_TESTING_EXPECT_EQUAL(lower_from_upper.precision(),precision);

    ARIADNE_TESTING_EXPECT_NO_THROW(lower_assigned.generic());
    ARIADNE_TESTING_EXPECT_NO_THROW(upper_assigned.generic());
    ARIADNE_TESTING_EXPECT_EQUAL(lower_assigned.characteristics(),precision);
    ARIADNE_TESTING_EXPECT_EQUAL(upper_assigned.characteristics(),precision);
    ARIADNE_TESTING_EXPECT_NO_THROW(lower_assigned.get_d());
    ARIADNE_TESTING_EXPECT_NO_THROW(upper_assigned.get_d());

    FloatBallType ball(F(2,precision),F(1,precision));
    L lower_from_ball(ball);
    U upper_from_ball(ball);
    ARIADNE_TESTING_EXPECT_EQUAL(lower_from_ball.raw(),F(1,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(upper_from_ball.raw(),F(3,precision));

    FloatError<PR> error(1u,precision);
    U upper_from_error(error);
    ARIADNE_TESTING_EXPECT_EQUAL(upper_from_error.raw(),F(1,precision));

    U negative_upper(-1,precision);
    ARIADNE_TESTING_EXPECT_EXCEPTION(PositiveFloatUpperBoundType{negative_upper});

    ValidatedUpperNumber generic_negative=negative_upper;
    PositiveValidatedUpperNumber positive_generic_negative(generic_negative);
    ARIADNE_TESTING_EXPECT_EXCEPTION(PositiveFloatUpperBoundType(positive_generic_negative,precision));

    if constexpr (Same<PR,MultiplePrecision>) {
        FloatDPLowerBound lower_dp(1,dp);
        FloatDPUpperBound upper_dp(1,dp);
        L lower_cross(lower_dp,precision);
        U upper_cross(upper_dp,precision);
        ARIADNE_TESTING_EXPECT_EQUAL(lower_cross.raw(),F(1,precision));
        ARIADNE_TESTING_EXPECT_EQUAL(upper_cross.raw(),F(1,precision));
    }

    L lower_positive(2,precision);
    U upper_positive(2,precision);
    ARIADNE_TESTING_EXPECT_EQUAL(nul(lower_positive).raw(),F(0,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(pos(lower_positive).raw(),F(2,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(hlf(lower_positive).raw(),F(1,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(nul(upper_positive).raw(),F(0,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(pos(upper_positive).raw(),F(2,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(hlf(upper_positive).raw(),F(1,precision));

    ARIADNE_TESTING_EXPECT_GREATER(log(lower_positive).raw(),F(0,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(exp(upper_zero).raw(),F(1,precision));
    ARIADNE_TESTING_EXPECT_NO_THROW(abs(L(-2,precision)));
    ARIADNE_TESTING_EXPECT_NO_THROW(abs(U(-2,precision)));

    L lower_above(2,precision);
    U upper_below(1,precision);
    L lower_below(1,precision);
    U upper_above(2,precision);
    ARIADNE_TESTING_EXPECT_NO_THROW(eq(lower_above,upper_below));
    ARIADNE_TESTING_EXPECT_NO_THROW(eq(lower_below,upper_above));
    ARIADNE_TESTING_EXPECT_NO_THROW(eq(upper_below,lower_above));
    ARIADNE_TESTING_EXPECT_NO_THROW(eq(upper_above,lower_below));
}

template<class PR> Void
TestDirectedFloats<PR>::test_conversions()
{
    Rational one=1;
    Rational four_thirds=4*one/3;
    Rational five_thirds=5*one/3;
    Rational neg_five_thirds=-5*one/3;

    ValidatedLowerNumber l(four_thirds);
    ValidatedUpperNumber u(five_thirds);

    ARIADNE_TESTING_EXPECT_LESS_EQUAL(FloatBoundsType(l,u,precision).lower().raw(),four_thirds);
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(FloatBoundsType(l,u,precision).upper().raw(),five_thirds);

    ARIADNE_TESTING_EXPECT_LESS_EQUAL(FloatLowerBoundType(five_thirds,precision).raw(),five_thirds);
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(FloatLowerBoundType(neg_five_thirds,precision).raw(),neg_five_thirds);
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(FloatUpperBoundType(five_thirds,precision).raw(),five_thirds);
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(FloatUpperBoundType(neg_five_thirds,precision).raw(),neg_five_thirds);

    ARIADNE_TESTING_EXPECT_EQUAL(cast_integer(FloatUpperBound<PR>(Dyadic(5,2u),precision)),Integer(2));
    ARIADNE_TESTING_EXPECT_EQUAL(cast_integer(FloatLowerBound<PR>(Dyadic(11,2u),precision)),Integer(2));
    
    // NaN must not be interpreted as a negative error on either backend.
    // The mixed integer-comparison contract below is tested for FloatMP:
    // FloatDP currently treats NaN differently for <=, >= and !=, so
    // imposing MP/MP semantics on it would change unrelated behavior.
    const Float<PR> nan_value=Float<PR>::nan(precision);
    ARIADNE_TESTING_PROVE(not (nan_value<0));
    ARIADNE_TESTING_PROVE(not (nan_value>0));
    if constexpr (Same<PR,MultiplePrecision>) {
        ARIADNE_TESTING_PROVE(not (nan_value<=0));
        ARIADNE_TESTING_PROVE(not (nan_value>=0));
        ARIADNE_TESTING_PROVE(nan_value!=0);
    }
    ARIADNE_TESTING_EXPECT_NO_THROW(static_cast<void>(FloatError<PR>{nan_value}));

    FloatLowerBoundType lower_bound(five_thirds,precision);
    FloatUpperBoundType upper_bound(five_thirds,precision);
    DyadicLowerBound dyadic_lower_bound(lower_bound);
    DyadicUpperBound dyadic_upper_bound(upper_bound);

    ARIADNE_TESTING_EXPECT_EQUAL(dyadic_lower_bound.raw(),Dyadic(lower_bound.raw()));
    ARIADNE_TESTING_EXPECT_EQUAL(dyadic_upper_bound.raw(),Dyadic(upper_bound.raw()));
    ARIADNE_TESTING_EXPECT_EQUAL(dyadic_lower_bound.get(precision).raw(),lower_bound.raw());
    ARIADNE_TESTING_EXPECT_EQUAL(dyadic_upper_bound.get(precision).raw(),upper_bound.raw());

    if constexpr (Same<PR,DoublePrecision>) {
        ARIADNE_TESTING_EXPECT_EQUAL(class_name<FloatLowerBoundType>(),String("FloatDPLowerBound"));
        ARIADNE_TESTING_EXPECT_EQUAL(class_name<FloatUpperBoundType>(),String("FloatDPUpperBound"));
        ARIADNE_TESTING_EXPECT_EQUAL(class_name<FloatError<PR>>(),String("FloatDPError"));
    } else {
        ARIADNE_TESTING_EXPECT_EQUAL(class_name<FloatLowerBoundType>(),String("FloatMPLowerBound"));
        ARIADNE_TESTING_EXPECT_EQUAL(class_name<FloatUpperBoundType>(),String("FloatMPUpperBound"));
        ARIADNE_TESTING_EXPECT_EQUAL(class_name<FloatError<PR>>(),String("FloatMPError"));
    }
}

template<class PR> Void
TestDirectedFloats<PR>::test_validation() {
    Rational one=1;
    Rational two_=2;
    ARIADNE_TESTING_PROVE(refines(FloatLowerBoundType(one,precision),FloatLowerBoundType(-one,precision)));
    ARIADNE_TESTING_PROVE(refines(FloatUpperBoundType(-one,precision),FloatUpperBoundType(+one,precision)));
    ARIADNE_TESTING_PROVE(refines(FloatUpperBoundType(-two_,precision),FloatUpperBoundType(-one,precision)));
//    ARIADNE_TESTING_PROVE(refines(rec(FloatUpperBoundType(-two,precision)),rec(FloatUpperBoundType(-one,precision))));
}

template<class PR> Void
TestDirectedFloats<PR>::test_rounded_arithmetic() {
    Rational one=1;
    Rational third=one/3;
    Rational fifth=one/3;
    ARIADNE_TESTING_EXPECT_LESS_EQUAL((FloatLowerBoundType(third,precision)+FloatLowerBoundType(fifth,precision)).raw(),third+fifth);
    ARIADNE_TESTING_EXPECT_LESS_EQUAL((FloatLowerBoundType(third,precision)-FloatUpperBoundType(fifth,precision)).raw(),third-fifth);
    ARIADNE_TESTING_PROVE(refines(FloatUpperBoundType(third+fifth,precision),FloatUpperBoundType(third,precision)+FloatUpperBoundType(fifth,precision)));
    ARIADNE_TESTING_PROVE(refines(FloatLowerBoundType(third+fifth,precision),FloatLowerBoundType(third,precision)+FloatLowerBoundType(fifth,precision)));
    ARIADNE_TESTING_EXPECT_LESS_EQUAL((PositiveFloatLowerBoundType(third,precision)*PositiveFloatLowerBoundType(fifth,precision)).raw(),third*fifth);
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL((PositiveFloatUpperBoundType(third,precision)*PositiveFloatUpperBoundType(fifth,precision)).raw(),third*fifth);
    ARIADNE_TESTING_EXPECT_LESS_EQUAL((PositiveFloatLowerBoundType(third,precision)/PositiveFloatUpperBoundType(fifth,precision)).raw(),third/fifth);
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL((PositiveFloatUpperBoundType(third,precision)/PositiveFloatLowerBoundType(fifth,precision)).raw(),third/fifth);

    FloatLowerBoundType lower_one(1,precision);
    FloatUpperBoundType upper_one(1,precision);
    FloatLowerBoundType lower_two(2,precision);
    FloatUpperBoundType upper_two(2,precision);
    ARIADNE_TESTING_EXPECT_EQUAL(nul(lower_two).raw(),FloatValueType(0,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(nul(upper_two).raw(),FloatValueType(0,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(pos(lower_two).raw(),FloatValueType(2,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(pos(upper_two).raw(),FloatValueType(2,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(hlf(lower_two).raw(),FloatValueType(1,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(hlf(upper_two).raw(),FloatValueType(1,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(log(lower_one).raw(),FloatValueType(0,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(exp(FloatUpperBoundType(0,precision)).raw(),FloatValueType(1,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(abs(FloatLowerBoundType(-2,precision)).raw(),FloatValueType(2,precision));
    ARIADNE_TESTING_EXPECT_EQUAL(abs(FloatUpperBoundType(-2,precision)).raw(),FloatValueType(2,precision));
}

template<class PR> Void
TestDirectedFloats<PR>::test_comparison() {
    PR pr=precision;

    {
        PositiveFloatUpperBoundType one(1u,pr);
        PositiveFloatUpperBoundType third=one/3u;

        ARIADNE_TESTING_PROVE(definitely(third+third < 1u));
        ARIADNE_TESTING_PROVE(definitely(third+third <= 1u));
        ARIADNE_TESTING_PROVE(definitely(5u*third < 2u));
        ARIADNE_TESTING_PROVE(definitely(5u*third <= 2u));
        ARIADNE_TESTING_PROVE(possibly(3u*third > 1u));
        ARIADNE_TESTING_PROVE(possibly(3u*third >= 1u));
    }

    {
        FloatLowerBoundType lower_one(1,pr);
        FloatLowerBoundType lower_two(2,pr);
        FloatUpperBoundType upper_one(1,pr);
        FloatUpperBoundType upper_two(2,pr);

        ARIADNE_TESTING_PROVE(not possibly(eq(lower_two,upper_one)));
        ARIADNE_TESTING_PROVE(definitely(lower_two!=upper_one));
        ARIADNE_TESTING_PROVE(is_indeterminate(eq(lower_one,upper_two)));
        ARIADNE_TESTING_PROVE(not possibly(eq(upper_one,lower_two)));
        ARIADNE_TESTING_PROVE(is_indeterminate(eq(upper_two,lower_one)));
    }

    {
        PositiveFloatLowerBoundType one(1u,pr);
        PositiveFloatLowerBoundType two_thirds=one*2u/3u;
        ARIADNE_TESTING_PROVE(possibly(two_thirds+two_thirds < 2u));
        ARIADNE_TESTING_PROVE(possibly(two_thirds+two_thirds <= 2u));
        ARIADNE_TESTING_PROVE(definitely(2u*two_thirds > 1u));
        ARIADNE_TESTING_PROVE(definitely(2u*two_thirds >= 1u));
        ARIADNE_TESTING_PROVE(possibly(3u*two_thirds < 2u));
        ARIADNE_TESTING_PROVE(possibly(3u*two_thirds <= 2u));
    }

}

template<class PR>
void add_directed_float_cases(Ariadne::Testing::Suite& suite, std::string const& variant, PR precision) {
    using Fixture=TestDirectedFloats<PR>;
    auto add=[&](std::string const& name, auto method) {
        suite.add(variant+" / "+name, [precision,method](Ariadne::Testing::Context& context) {
            Fixture fixture(context,precision);
            (fixture.*method)();
        });
    };
    add("precision", &Fixture::test_precision);
    add("header API", &Fixture::test_header_api);
    add("conversions", &Fixture::test_conversions);
    add("validation", &Fixture::test_validation);
    add("rounded arithmetic", &Fixture::test_rounded_arithmetic);
    add("comparison", &Fixture::test_comparison);
}

int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("directed floats");
    add_directed_float_cases<DoublePrecision>(suite,"FloatDP",dp);
    add_directed_float_cases<MultiplePrecision>(suite,"FloatMP",MultiplePrecision(128_bits));
    return Ariadne::Testing::run_suite(suite, argc, argv);
}

