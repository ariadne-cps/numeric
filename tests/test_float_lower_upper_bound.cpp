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
#include "numeric/floats.hpp"

#include "utility/test.hpp"
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
    : public TestFloats<PR>
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
    TestDirectedFloats(PR prec) : precision(prec) { };
    Void test();
  private:
    Void test_precision();
    Void test_conversions();
    Void test_validation();
    Void test_rounded_arithmetic();
    Void test_comparison();
};

template<class PR> Void
TestDirectedFloats<PR>::test()
{
    ARIADNE_TEST_CALL(test_precision());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_validation());
    ARIADNE_TEST_CALL(test_rounded_arithmetic());
    ARIADNE_TEST_CALL(test_comparison());
}

template<class PR> Void
TestDirectedFloats<PR>::test_precision()
{
    FloatLowerBoundType lx(Rational(1),precision);
    ARIADNE_TEST_EQUALS(max(lx,lx).precision(),precision);
    ARIADNE_TEST_EQUALS(min(lx,lx).precision(),precision);
    ARIADNE_TEST_EQUALS((lx+lx).precision(),precision);
    ARIADNE_TEST_EQUALS(exp(lx).precision(),precision);
    ARIADNE_TEST_EQUALS((lx+2u).precision(),precision);
    ARIADNE_TEST_EQUALS((lx*2u).precision(),precision);
    ARIADNE_TEST_EQUALS((lx/2u).precision(),precision);
    ARIADNE_TEST_EQUALS((lx+2).precision(),precision);
    ARIADNE_TEST_EQUALS((lx-2).precision(),precision);
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

    ARIADNE_TEST_COMPARE(FloatBoundsType(l,u,precision).lower().raw(),<=,four_thirds);
    ARIADNE_TEST_COMPARE(FloatBoundsType(l,u,precision).upper().raw(),>=,five_thirds);

    ARIADNE_TEST_COMPARE(FloatLowerBoundType(five_thirds,precision).raw(),<=,five_thirds);
    ARIADNE_TEST_COMPARE(FloatLowerBoundType(neg_five_thirds,precision).raw(),<=,neg_five_thirds);
    ARIADNE_TEST_COMPARE(FloatUpperBoundType(five_thirds,precision).raw(),>=,five_thirds);
    ARIADNE_TEST_COMPARE(FloatUpperBoundType(neg_five_thirds,precision).raw(),>=,neg_five_thirds);

    ARIADNE_TEST_EQUALS(cast_integer(FloatUpperBound<PR>(Dyadic(5,2u),precision)),Integer(2));
    ARIADNE_TEST_EQUALS(cast_integer(FloatLowerBound<PR>(Dyadic(11,2u),precision)),Integer(2));
    
    // Test that FloatError can be constructed from NaN
    ARIADNE_TEST_EXECUTE(FloatError<PR>(Float<PR>::nan(precision)));

    FloatLowerBoundType lower_bound(five_thirds,precision);
    FloatUpperBoundType upper_bound(five_thirds,precision);
    DyadicLowerBound dyadic_lower_bound(lower_bound);
    DyadicUpperBound dyadic_upper_bound(upper_bound);

    ARIADNE_TEST_EQUALS(dyadic_lower_bound.raw(),Dyadic(lower_bound.raw()));
    ARIADNE_TEST_EQUALS(dyadic_upper_bound.raw(),Dyadic(upper_bound.raw()));
    ARIADNE_TEST_EQUALS(dyadic_lower_bound.get(precision).raw(),lower_bound.raw());
    ARIADNE_TEST_EQUALS(dyadic_upper_bound.get(precision).raw(),upper_bound.raw());

    if constexpr (Same<PR,DoublePrecision>) {
        ARIADNE_TEST_EQUALS(class_name<FloatLowerBoundType>(),String("FloatDPLowerBound"));
        ARIADNE_TEST_EQUALS(class_name<FloatUpperBoundType>(),String("FloatDPUpperBound"));
        ARIADNE_TEST_EQUALS(class_name<FloatError<PR>>(),String("FloatDPError"));
    } else {
        ARIADNE_TEST_EQUALS(class_name<FloatLowerBoundType>(),String("FloatMPLowerBound"));
        ARIADNE_TEST_EQUALS(class_name<FloatUpperBoundType>(),String("FloatMPUpperBound"));
        ARIADNE_TEST_EQUALS(class_name<FloatError<PR>>(),String("FloatMPError"));
    }
}

template<class PR> Void
TestDirectedFloats<PR>::test_validation() {
    Rational one=1;
    Rational two_=2;
    ARIADNE_TEST_ASSERT(refines(FloatLowerBoundType(one,precision),FloatLowerBoundType(-one,precision)));
    ARIADNE_TEST_ASSERT(refines(FloatUpperBoundType(-one,precision),FloatUpperBoundType(+one,precision)));
    ARIADNE_TEST_ASSERT(refines(FloatUpperBoundType(-two_,precision),FloatUpperBoundType(-one,precision)));
//    ARIADNE_TEST_ASSERT(refines(rec(FloatUpperBoundType(-two,precision)),rec(FloatUpperBoundType(-one,precision))));
}

template<class PR> Void
TestDirectedFloats<PR>::test_rounded_arithmetic() {
    Rational one=1;
    Rational third=one/3;
    Rational fifth=one/3;
    ARIADNE_TEST_COMPARE((FloatLowerBoundType(third,precision)+FloatLowerBoundType(fifth,precision)).raw(),<=,third+fifth);
    ARIADNE_TEST_COMPARE((FloatLowerBoundType(third,precision)-FloatUpperBoundType(fifth,precision)).raw(),<=,third-fifth);
    ARIADNE_TEST_ASSERT(refines(FloatUpperBoundType(third+fifth,precision),FloatUpperBoundType(third,precision)+FloatUpperBoundType(fifth,precision)));
    ARIADNE_TEST_ASSERT(refines(FloatLowerBoundType(third+fifth,precision),FloatLowerBoundType(third,precision)+FloatLowerBoundType(fifth,precision)));
    ARIADNE_TEST_COMPARE((PositiveFloatLowerBoundType(third,precision)*PositiveFloatLowerBoundType(fifth,precision)).raw(),<=,third*fifth);
    ARIADNE_TEST_COMPARE((PositiveFloatUpperBoundType(third,precision)*PositiveFloatUpperBoundType(fifth,precision)).raw(),>=,third*fifth);
    ARIADNE_TEST_COMPARE((PositiveFloatLowerBoundType(third,precision)/PositiveFloatUpperBoundType(fifth,precision)).raw(),<=,third/fifth);
    ARIADNE_TEST_COMPARE((PositiveFloatUpperBoundType(third,precision)/PositiveFloatLowerBoundType(fifth,precision)).raw(),>=,third/fifth);
}

template<class PR> Void
TestDirectedFloats<PR>::test_comparison() {
    PR pr=precision;

    {
        PositiveFloatUpperBoundType one(1u,pr);
        PositiveFloatUpperBoundType third=one/3u;

        ARIADNE_TEST_ASSERT(definitely(third+third < 1u));
        ARIADNE_TEST_ASSERT(definitely(third+third <= 1u));
        ARIADNE_TEST_ASSERT(definitely(5u*third < 2u));
        ARIADNE_TEST_ASSERT(definitely(5u*third <= 2u));
        ARIADNE_TEST_ASSERT(possibly(3u*third > 1u));
        ARIADNE_TEST_ASSERT(possibly(3u*third >= 1u));
    }

    {
        PositiveFloatLowerBoundType one(1u,pr);
        PositiveFloatLowerBoundType two_thirds=one*2u/3u;
        ARIADNE_TEST_ASSERT(possibly(two_thirds+two_thirds < 2u));
        ARIADNE_TEST_ASSERT(possibly(two_thirds+two_thirds <= 2u));
        ARIADNE_TEST_ASSERT(definitely(2u*two_thirds > 1u));
        ARIADNE_TEST_ASSERT(definitely(2u*two_thirds >= 1u));
        ARIADNE_TEST_ASSERT(possibly(3u*two_thirds < 2u));
        ARIADNE_TEST_ASSERT(possibly(3u*two_thirds <= 2u));
    }

}

Int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    TestDirectedFloats<DoublePrecision>(dp).test();
    TestDirectedFloats<MultiplePrecision>(MultiplePrecision(128_bits)).test();

    return ARIADNE_TEST_FAILURES;
}

