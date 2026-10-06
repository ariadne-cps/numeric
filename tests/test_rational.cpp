/***************************************************************************
 *            test_rational.cpp
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


#include "numeric/rational.hpp"
#include "numeric/builtin.hpp"
#include "numeric/integer.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/decimal.hpp"
#include "numeric/number.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"
#include "foundation/logical.hpp"

#include <iomanip>
#include <sstream>

#include "utility/test.hpp"

using namespace std;
using namespace Ariadne;


class TestRational
{
  public:
    void test();
  private:
    void test_concept();
    void test_literal();
    void test_conversions();
    void test_arithmetic();
    void test_rounding();
    void test_comparisons();
    void test_infinity();
    void test_bounds();

    void test_decimal();

};

void TestRational::test()
{
    ARIADNE_TEST_CALL(test_literal());
    ARIADNE_TEST_CALL(test_conversions());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_rounding());
    ARIADNE_TEST_CALL(test_comparisons());
    ARIADNE_TEST_CALL(test_infinity());
    ARIADNE_TEST_CALL(test_bounds());

    ARIADNE_TEST_CALL(test_decimal());
}

void TestRational::test_concept() {
    unsigned int m=1; unsigned long int lm=1; int n=-2; long int ln=-2; Integer z=-5; Dyadic w=z;
    Rational q, q2; Boolean b;

    q=Rational(); q=Rational(m); q=Rational(lm); q=Rational(n); q=Rational(ln); q=Rational(z); q=Rational(z);
    q2=Rational();
    q=m; q=lm; q=n; q=ln; q=z; q=q2;

    q=+q; q=-q;
    q=q+q; q=q-q; q=q*q; q=q/q;

    q=q+n; q=q-n; q=q*n; q=q/n;
    q=n+q; q=n-q; q=n*q; q=n/q;
    q=q+z; q=q-z; q=q*z; q=q/z;
    q=z+q; q=z-q; q=z*q; q=z/q;
    q=q+w; q=q-w; q=q*w; q=q/w;
    q=w+q; q=w-q; q=w*q; q=w/q;

    q=max(q,q); q=min(q,q); q=abs(q);
    q=pos(q); q=neg(q); q=sqr(q); q=rec(q);

    q=1.5_q; q=3/2_q; q=-1.3_q;

    b=(q==q); b=(q!=q); b=(q<=q); b=(q>=q); b=(q<q); b=(q>q);
    b=(q==n); b=(q!=n); b=(q<=n); b=(q>=n); b=(q<n); b=(q>n);
    b=(n==q); b=(n!=q); b=(n<=q); b=(n>=q); b=(n<q); b=(n>q);
    b=(q==z); b=(q!=z); b=(q<=z); b=(q>=z); b=(q<z); b=(q>z);
    b=(z==q); b=(z!=q); b=(z<=q); b=(z>=q); b=(z<q); b=(z>q);
    b=(q==w); b=(q!=w); b=(q<=w); b=(q>=w); b=(q<w); b=(q>w);
    b=(w==q); b=(w!=q); b=(w<=q); b=(w>=q); b=(w<q); b=(w>q);
}

void TestRational::test_literal() {
    ARIADNE_TEST_CONSTRUCT(Rational,q,(3.25_q));
    ARIADNE_TEST_EQUALS(q,Rational(13,4));
    ARIADNE_TEST_EQUALS(3.25_q,Rational(13,4));
    ARIADNE_TEST_EQUALS(-11.375_q,Rational(-91,8));
    ARIADNE_TEST_EQUALS(10.3_q,Rational(103,10));
    ARIADNE_TEST_EQUALS(0.333333333333333333_q,Rational(1,3));
    ARIADNE_TEST_EQUALS(0.2857142857142857_q,Rational(2,7));
    ARIADNE_TEST_FAIL(0.453591850358036834_q);
    ARIADNE_TEST_FAIL(3.1415926535897931_q);
}

void TestRational::test_conversions() {
    ARIADNE_TEST_EQUAL(Rational(Integer(-3)),Rational(-3,1));
    ARIADNE_TEST_EQUAL(Rational(Dyadic(-13)),Rational(-13));
    ARIADNE_TEST_EQUAL(Rational(Dyadic(-13,3u)),Rational(-13,8));
    ARIADNE_TEST_EQUAL(Rational(Int64(int32_t(-7))),Rational(-7));

    FloatDP fdp(Dyadic(5,2u),dp);
    ARIADNE_TEST_EQUAL(Rational(fdp),Rational(5,4));
    MultiplePrecision mp(128);
    FloatMP fmp(Dyadic(5,2u),mp);
    ARIADNE_TEST_EQUAL(Rational(fmp),Rational(5,4));

    ARIADNE_TEST_EQUAL(Rational(String("7/9")),Rational(7,9));
    ARIADNE_TEST_FAIL(Rational(String("not-a-rational")));

    mpq_t raw;
    mpq_init(raw);
    mpq_set_si(raw,5,6);
    Rational from_raw(raw);
    mpq_clear(raw);
    ARIADNE_TEST_EQUAL(from_raw,Rational(5,6));

    Rational assigned;
    Rational source(11,13);
    assigned=source;
    ARIADNE_TEST_EQUAL(assigned,source);

    std::istringstream input("17/19");
    Rational streamed;
    input >> streamed;
    ARIADNE_TEST_EQUAL(streamed,Rational(17,19));
}

void TestRational::test_arithmetic() {
    ARIADNE_TEST_EQUAL(Rational(4,5)+Rational(-2,7),Rational(18,35));
    ARIADNE_TEST_EQUAL(Rational(-4,5)-Rational(-2,7),Rational(-18,35));
    ARIADNE_TEST_EQUAL(Rational(4,5)*Rational(-2,7),Rational(-8,35));
    ARIADNE_TEST_EQUAL(Rational(4,5)/Rational(-2,7),Rational(-14,5));
    ARIADNE_TEST_EQUAL(div(Integer(3),Integer(4)),Rational(3,4));

    PositiveRational p2=cast_positive(Rational(2));
    PositiveRational p3=cast_positive(Rational(3));
    ARIADNE_TEST_EQUAL(max(Rational(-1),p2),Rational(2));
    ARIADNE_TEST_EQUAL(max(p2,Rational(4)),Rational(4));
    ARIADNE_TEST_EQUAL(max(p2,p3),Rational(3));
    ARIADNE_TEST_EQUAL(min(p2,p3),Rational(2));
    ARIADNE_TEST_EQUAL(mag(Rational(-5,3)),Rational(5,3));
    ARIADNE_TEST_EQUAL(mig(Rational(-5,3)),Rational(5,3));
}

void TestRational::test_rounding() {
    ARIADNE_TEST_EQUALS(round(Rational(0)),Integer(0));
    ARIADNE_TEST_EQUALS(round(Rational(3)),Integer(3));
    ARIADNE_TEST_EQUALS(round(Rational(-11,4)),Integer(-3));
    ARIADNE_TEST_EQUALS(round(Rational(-10,4)),Integer(-3));
    ARIADNE_TEST_EQUALS(round(Rational(-9,4)),Integer(-2));
    ARIADNE_TEST_EQUALS(round(Rational(9,4)),Integer(2));
    ARIADNE_TEST_EQUALS(round(Rational(10,4)),Integer(3));
    ARIADNE_TEST_EQUALS(round(Rational(11,4)),Integer(3));

    ARIADNE_TEST_EQUALS(ceil(Rational(-13,4)),Integer(-3));
    ARIADNE_TEST_EQUALS(ceil(Rational(-12,4)),Integer(-3));
    ARIADNE_TEST_EQUALS(ceil(Rational(0,4)),Integer(0));
    ARIADNE_TEST_EQUALS(ceil(Rational(12,4)),Integer(3));
    ARIADNE_TEST_EQUALS(ceil(Rational(13,4)),Integer(4));

    ARIADNE_TEST_EQUALS(floor(Rational(-13,4)),Integer(-4));
    ARIADNE_TEST_EQUALS(floor(Rational(-12,4)),Integer(-3));
    ARIADNE_TEST_EQUALS(floor(Rational(0,4)),Integer(0));
    ARIADNE_TEST_EQUALS(floor(Rational(12,4)),Integer(3));
    ARIADNE_TEST_EQUALS(floor(Rational(13,4)),Integer(3));
}

void TestRational::test_comparisons() {
    ExactDouble infinity_value=ExactDouble::inf();
    ExactDouble max=ExactDouble(std::numeric_limits<double>::max());
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(-max),Rational(+max));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(-max),Rational(-4,5));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(-max),Rational(2,3));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(-4,5),Rational(-2,7));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,-infinity_value,Rational(18,35));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(18,35),+infinity_value);

    ARIADNE_TEST_EQUALS(cmp(Rational(1),ExactDouble(2)),Comparison::LESS);
    ARIADNE_TEST_EQUALS(cmp(Rational(2),ExactDouble(2)),Comparison::EQUAL);
    ARIADNE_TEST_EQUALS(cmp(Rational(3),ExactDouble(2)),Comparison::GREATER);
    ARIADNE_TEST_EQUALS(cmp(Rational(1),ExactDouble::inf()),Comparison::LESS);
    ARIADNE_TEST_EQUALS(cmp(Rational(1),-ExactDouble::inf()),Comparison::GREATER);
    ARIADNE_TEST_EQUALS(cmp(ExactDouble(2),Rational(1)),Comparison::GREATER);
}

void TestRational::test_infinity() {
    Rational qinf=Rational::inf();
    Rational qninf=Rational::inf(Sign::NEGATIVE);
    Rational qnan=Rational::nan();

    ARIADNE_TEST_ASSERT(is_nan(Rational::nan()));
    ARIADNE_TEST_ASSERT(is_inf(Rational::inf()));
    ARIADNE_TEST_ASSERT(is_inf(Rational::inf(Sign(+1))));
    ARIADNE_TEST_ASSERT(is_inf(Rational::inf(Sign(-1))));
    ARIADNE_TEST_ASSERT(is_finite(Rational(0)));
    ARIADNE_TEST_ASSERT(is_zero(Rational(0)));
    ARIADNE_TEST_EQUALS(Rational::inf(Sign(+1)),Rational::inf());
    ARIADNE_TEST_ASSERT(Rational::inf(Sign(+1))>Rational(0));
    ARIADNE_TEST_ASSERT(Rational::inf(Sign(-1))<Rational(0));

    ARIADNE_TEST_EQUALS(sgn(Rational::inf()), Sign::POSITIVE);
    ARIADNE_TEST_EQUALS(sgn(Rational::nan()), Sign::ZERO);
    ARIADNE_TEST_EQUALS(sgn(-Rational::inf()), Sign::NEGATIVE);

    ARIADNE_TEST_ASSERT(std::isnan(Rational::nan().get_d()));
    ARIADNE_TEST_EQUAL(Rational::inf(Sign::POSITIVE).get_d(),std::numeric_limits<double>::infinity());
    ARIADNE_TEST_EQUAL(Rational::inf(Sign::NEGATIVE).get_d(),-std::numeric_limits<double>::infinity());

    ARIADNE_TEST_BINARY_PREDICATE(operator==,Rational(2,0),Rational(1,0));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(-1,0),Rational(2,0));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(-1,2),Rational(1,0));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(0,1),Rational(1,0));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(3,2),Rational(1,0));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(-1,0),Rational(1,2));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(-1,0),Rational(-3,2));

    ARIADNE_TEST_BINARY_PREDICATE(operator==,Rational::inf(),Rational::inf());
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(-1,4),Rational::inf());
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(0),Rational::inf());
    ARIADNE_TEST_BINARY_PREDICATE(operator<,Rational(3,2),Rational::inf());
    ARIADNE_TEST_BINARY_PREDICATE(operator<,-Rational::inf(),Rational::inf());
    ARIADNE_TEST_BINARY_PREDICATE(operator<,-Rational::inf(),Rational(1,4));
    ARIADNE_TEST_BINARY_PREDICATE(operator<,-Rational::inf(),Rational(-3,2));

    ARIADNE_TEST_BINARY_PREDICATE(operator!=,Rational::inf(),Rational(1,2));
    ARIADNE_TEST_BINARY_PREDICATE(operator>=,Rational::inf(),Rational(1,2));
    ARIADNE_TEST_BINARY_PREDICATE(operator> ,Rational::inf(),Rational(1,2));
    ARIADNE_TEST_BINARY_PREDICATE(operator!=,Rational::inf(Sign::NEGATIVE),Rational(-1,2));
    ARIADNE_TEST_BINARY_PREDICATE(operator<=,Rational::inf(Sign::NEGATIVE),Rational(-1,2));
    ARIADNE_TEST_BINARY_PREDICATE(operator< ,Rational::inf(Sign::NEGATIVE),Rational(-1,2));

    ARIADNE_TEST_ASSERT(is_nan(-Rational::nan()));
    ARIADNE_TEST_EQUAL(-Rational::inf(Sign::POSITIVE),Rational::inf(Sign::NEGATIVE));
    ARIADNE_TEST_EQUAL(-Rational::inf(Sign::NEGATIVE),Rational::inf(Sign::POSITIVE));

    ARIADNE_TEST_ASSERT(is_nan(Rational::inf()+(-Rational::inf())));
    ARIADNE_TEST_EQUALS(Rational::inf()+Rational::inf(),Rational::inf());
    ARIADNE_TEST_EQUALS(Rational::inf()+Rational(-2),Rational::inf());

    ARIADNE_TEST_ASSERT(is_nan(Rational::inf()-Rational::inf()));
    ARIADNE_TEST_EQUALS(Rational(2)-Rational::inf(),-Rational::inf())

    ARIADNE_TEST_ASSERT(is_nan(Rational::inf(Sign(+1))*Rational(0)));
    ARIADNE_TEST_ASSERT(is_nan(Rational::inf(Sign(-1))*Rational(0)));
    ARIADNE_TEST_ASSERT(is_nan(Rational(0)*Rational::inf(Sign(+1))));
    ARIADNE_TEST_ASSERT(is_nan(Rational(0)*Rational::inf(Sign(-1))));
    ARIADNE_TEST_EQUALS(Rational::inf(Sign(+1))*Rational::inf(Sign(+1)),Rational::inf(Sign(+1)));
    ARIADNE_TEST_EQUALS(Rational::inf(Sign(+1))*Rational::inf(Sign(-1)),Rational::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(Rational::inf(Sign(-1))*Rational::inf(Sign(+1)),Rational::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(Rational::inf(Sign(-1))*Rational::inf(Sign(-1)),Rational::inf(Sign(+1)));
    ARIADNE_TEST_EQUALS(Rational(+2)*Rational::inf(Sign(+1)),Rational::inf(Sign(+1)));
    ARIADNE_TEST_EQUALS(Rational(+2)*Rational::inf(Sign(-1)),Rational::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(Rational(-2)*Rational::inf(Sign(+1)),Rational::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(Rational(-2)*Rational::inf(Sign(-1)),Rational::inf(Sign(+1)));
    ARIADNE_TEST_EQUALS(Rational::inf(Sign(+1))*Rational(+2),Rational::inf(Sign(+1)));
    ARIADNE_TEST_EQUALS(Rational::inf(Sign(+1))*Rational(-2),Rational::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(Rational::inf(Sign(-1))*Rational(+2),Rational::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(Rational::inf(Sign(-1))*Rational(-2),Rational::inf(Sign(+1)));

    ARIADNE_TEST_EQUALS(Rational::inf(Sign(+1))/Rational(+2),Rational::inf(Sign(+1)));
    ARIADNE_TEST_EQUALS(Rational::inf(Sign(+1))/Rational(-2),Rational::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(Rational::inf(Sign(-1))/Rational(+2),Rational::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(Rational::inf(Sign(-1))/Rational(-2),Rational::inf(Sign(+1)));
    ARIADNE_TEST_EQUALS(Rational(+2)/Rational::inf(Sign(+1)),Rational(0));
    ARIADNE_TEST_EQUALS(Rational(-2)/Rational::inf(Sign(+1)),Rational(0));
    ARIADNE_TEST_EQUALS(Rational(+2)/Rational::inf(Sign(-1)),Rational(0));
    ARIADNE_TEST_EQUALS(Rational(-2)/Rational::inf(Sign(-1)),Rational(0));

    ARIADNE_TEST_ASSERT(is_nan(abs(Rational::nan())));
    ARIADNE_TEST_EQUALS(abs(Rational::inf()),Rational::inf());
    ARIADNE_TEST_EQUALS(abs(-Rational::inf()),Rational::inf());

    ARIADNE_TEST_ASSERT(is_nan(max(Rational::nan(),Rational(0))));
    ARIADNE_TEST_ASSERT(is_nan(max(Rational::nan(),Rational::inf())));
    ARIADNE_TEST_EQUALS(max(Rational::inf(Sign(-1)),Rational::inf(Sign(-1))),Rational::inf(Sign(-1)));
    ARIADNE_TEST_EQUALS(max(Rational::inf(Sign(-1)),Rational(-2)),Rational(-2));
    ARIADNE_TEST_EQUALS(max(Rational(-2),Rational::inf(Sign(+1))),Rational::inf(Sign(+1)));

    ARIADNE_TEST_EQUALS(Rational(Dyadic::inf(Sign::POSITIVE)),qinf);
    ARIADNE_TEST_EQUALS(Rational(Dyadic::inf(Sign::NEGATIVE)),qninf);
    ARIADNE_TEST_ASSERT(is_nan(Rational(Dyadic::inf(Sign::ZERO))));
}

void TestRational::test_bounds() {
    auto check_bounds=[](RationalBounds const& actual, Rational const& lower, Rational const& upper) {
        ARIADNE_TEST_EQUALS(actual.lower_raw(),lower);
        ARIADNE_TEST_EQUALS(actual.upper_raw(),upper);
    };

    RationalBounds p(2,3);
    RationalBounds n(-3,-2);
    RationalBounds m(-2,3);

    check_bounds(p*p,4,9);
    check_bounds(p*n,-9,-4);
    check_bounds(p*m,-6,9);
    check_bounds(n*p,-9,-4);
    check_bounds(n*n,4,9);
    check_bounds(n*m,-9,6);
    check_bounds(m*p,-6,9);
    check_bounds(m*n,-9,6);
    check_bounds(m*m,-6,9);

    check_bounds(p/p,Rational(2,3),Rational(3,2));
    check_bounds(n/p,Rational(-3,2),Rational(-2,3));
    check_bounds(m/p,Rational(-1),Rational(3,2));
    check_bounds(p/n,Rational(-3,2),Rational(-2,3));
    check_bounds(n/n,Rational(2,3),Rational(3,2));
    check_bounds(m/n,Rational(-3,2),Rational(1));
    check_bounds(p/m,-Rational::inf(),Rational::inf());

    check_bounds(nul(m),0,0);
    check_bounds(pos(m),-2,3);
    check_bounds(neg(m),-3,2);
    check_bounds(hlf(RationalBounds(2,4)),1,2);
    check_bounds(sqr(p),4,9);
    check_bounds(sqr(n),4,9);
    check_bounds(sqr(m),0,9);
    check_bounds(rec(p),Rational(1,3),Rational(1,2));
    check_bounds(rec(n),Rational(-1,2),Rational(-1,3));
    check_bounds(rec(m),-Rational::inf(),Rational::inf());

    check_bounds(add(p,n),-1,1);
    check_bounds(sub(p,n),4,6);
    check_bounds(mul(p,n),-9,-4);
    check_bounds(div(p,n),Rational(-3,2),Rational(-2,3));

    check_bounds(abs(p),2,3);
    check_bounds(abs(n),2,3);
    check_bounds(abs(m),0,3);
    check_bounds(pow(m,Nat(2u)),0,9);
    check_bounds(pow(m,Nat(3u)),-8,27);
    check_bounds(pow(p,Int(-1)),Rational(1,3),Rational(1,2));
    check_bounds(max(p,n),2,3);
    check_bounds(min(p,n),-3,-2);

    RationalBounds a(1,2);
    RationalBounds b(3,4);
    RationalBounds overlap(1,3);
    RationalBounds overlap2(2,4);
    RationalBounds point(2,2);

    ARIADNE_TEST_ASSERT(definitely(point==point));
    ARIADNE_TEST_ASSERT(not possibly(a==b));
    ARIADNE_TEST_ASSERT(is_indeterminate(overlap==overlap2));

    ARIADNE_TEST_ASSERT(definitely(a!=b));
    ARIADNE_TEST_ASSERT(not possibly(point!=point));
    ARIADNE_TEST_ASSERT(is_indeterminate(overlap!=overlap2));

    ARIADNE_TEST_ASSERT(definitely(a<=b));
    ARIADNE_TEST_ASSERT(not possibly(b<=a));
    ARIADNE_TEST_ASSERT(is_indeterminate(overlap<=overlap2));

    ARIADNE_TEST_ASSERT(definitely(b>=a));
    ARIADNE_TEST_ASSERT(not possibly(a>=b));
    ARIADNE_TEST_ASSERT(is_indeterminate(overlap>=overlap2));

    ARIADNE_TEST_ASSERT(definitely(a<b));
    ARIADNE_TEST_ASSERT(not possibly(b<a));
    ARIADNE_TEST_ASSERT(is_indeterminate(RationalBounds(1,2)<RationalBounds(2,3)));

    ARIADNE_TEST_ASSERT(definitely(b>a));
    ARIADNE_TEST_ASSERT(not possibly(a>b));
    ARIADNE_TEST_ASSERT(is_indeterminate(RationalBounds(2,3)>RationalBounds(1,2)));

    ARIADNE_TEST_EQUALS(class_name<RationalBounds>(),String("RationalBounds"));

    auto check_dyadic=[](DyadicBounds const& actual, Dyadic const& lower, Dyadic const& upper) {
        ARIADNE_TEST_EQUALS(actual.lower_raw(),lower);
        ARIADNE_TEST_EQUALS(actual.upper_raw(),upper);
    };
    DyadicBounds dpb(Dyadic(2),Dyadic(3));
    DyadicBounds dnb(Dyadic(-3),Dyadic(-2));
    DyadicBounds dmb(Dyadic(-2),Dyadic(3));
    check_dyadic(dpb*dpb,Dyadic(4),Dyadic(9));
    check_dyadic(dpb*dnb,Dyadic(-9),Dyadic(-4));
    check_dyadic(dpb*dmb,Dyadic(-6),Dyadic(9));
    check_dyadic(dnb*dpb,Dyadic(-9),Dyadic(-4));
    check_dyadic(dnb*dnb,Dyadic(4),Dyadic(9));
    check_dyadic(dnb*dmb,Dyadic(-9),Dyadic(6));
    check_dyadic(dmb*dpb,Dyadic(-6),Dyadic(9));
    check_dyadic(dmb*dnb,Dyadic(-9),Dyadic(6));
    check_dyadic(dmb*dmb,Dyadic(-6),Dyadic(9));
    check_dyadic(mul(dpb,dnb),Dyadic(-9),Dyadic(-4));
    check_dyadic(pow(dmb,Int(3)),Dyadic(-8),Dyadic(27));
}

void TestRational::test_decimal() {
    Decimal d0;
    ARIADNE_TEST_EQUALS(d0,Decimal(0,0u));
    ARIADNE_TEST_CONSTRUCT(Decimal,d1,(23,1u));
    ARIADNE_TEST_CONSTRUCT(Decimal,d2,(-42,2u));
    ARIADNE_TEST_EQUALS(7_decimal,Decimal(7,0u));
    ARIADNE_TEST_EQUALS(8_dec,Decimal(8,0u));
    ARIADNE_TEST_EQUALS(Decimal(0.0),Decimal(0,0u));
    ARIADNE_TEST_EQUALS(Decimal(-3.14),Decimal(-314,2u));
    ARIADNE_TEST_EQUALS(Decimal(10000000000.0),Decimal(10000000000_z,0u));
    ARIADNE_TEST_EQUALS(Decimal("-3.14"),Decimal(-314,2u));
    ARIADNE_TEST_EQUALS(Decimal("+3.14"),Decimal(314,2u));
    ARIADNE_TEST_FAIL(Decimal("1.2.3"));
    ARIADNE_TEST_FAIL(Decimal("1a"));
    ARIADNE_TEST_FAIL(Decimal("/1"));
    ARIADNE_TEST_EQUALS(Decimal("3.141592653589793238462643383279"),Decimal(3141592653589793_z,15u)+Decimal(238462643383279_z,30u));
    ARIADNE_TEST_EQUALS(Decimal("3.141592653589793238462643383279"),Decimal("3141592653589793238462643383279"_z,30u));
    ARIADNE_TEST_FAIL(Decimal d(0.33333333333));
    ARIADNE_TEST_EQUALS(Decimal(Dyadic(7,3u)),Decimal(875,3u));
    ARIADNE_TEST_EQUALS(Decimal(Dyadic(140,1u)),Decimal(70,0u));
    ARIADNE_TEST_FAIL(static_cast<void>(Decimal(Dyadic::inf())));
    ARIADNE_TEST_EQUALS(Decimal("-3.14"),Decimal(-314,2u));
    ARIADNE_TEST_EQUALS(Decimal("-3.1400"),Decimal(-314,2u));
    ARIADNE_TEST_EQUALS(Decimal("-0.0031400"),Decimal(-314,5u));
    ARIADNE_TEST_EQUALS(Decimal("3140"),Decimal(3140,0u));
    ARIADNE_TEST_EQUALS(Decimal("3.14159265")*Decimal("3.14159265"),Decimal("9.8696043785340225"));
    ARIADNE_TEST_PRINT(Decimal("3.14159265")*Decimal("3.14159265"));
    ARIADNE_TEST_EQUALS(Rational(d1),Rational(23,10));
    ARIADNE_TEST_EQUALS(Rational(d2),Rational(-21,50));
    ARIADNE_TEST_EQUAL(Rational(+d2),+Rational(d2));
    ARIADNE_TEST_EQUAL(Rational(-d2),-Rational(d2));
    ARIADNE_TEST_EQUAL(Rational(d1+d2),Rational(d1)+Rational(d2));
    ARIADNE_TEST_EQUAL(Rational(d1-d2),Rational(d1)-Rational(d2));
    ARIADNE_TEST_EQUAL(Rational(d1*d2),Rational(d1)*Rational(d2));
    ARIADNE_TEST_EQUAL(d1/d2,Rational(d1)/Rational(d2));
    ARIADNE_TEST_EQUAL(Decimal(1,2u)/Decimal(1,1u),Rational(1,10));

    Decimal accumulated(1);
    accumulated+=Decimal("2.5");
    ARIADNE_TEST_EQUALS(accumulated,Decimal("3.5"));
    ARIADNE_TEST_EQUALS(nul(d1),Decimal(0));
    ARIADNE_TEST_EQUALS(sqr(Decimal("1.5")),Decimal("2.25"));
    ARIADNE_TEST_EQUALS(hlf(Decimal(24,1u)),Decimal(12,1u));
    ARIADNE_TEST_EQUALS(hlf(Decimal(3)),Decimal(15,1u));
    ARIADNE_TEST_EQUALS(abs(Decimal("-2.5")),Decimal("2.5"));
    ARIADNE_TEST_EQUALS(max(Decimal("-2.5"),Decimal("1.5")),Decimal("1.5"));
    ARIADNE_TEST_EQUALS(min(Decimal("-2.5"),Decimal("1.5")),Decimal("-2.5"));

    ARIADNE_TEST_EQUALS(Decimal("3.014").literal(),String("3.014"));
    ARIADNE_TEST_EQUALS(Decimal(3).literal(),String("3."));
    ARIADNE_TEST_EQUALS(class_name<DecimalBounds>(),String("DecimalBounds"));

    ExactNumber exact_decimal=static_cast<ExactNumber>(Decimal("1.25"));
    ARIADNE_TEST_EQUALS(exact_decimal.class_name(),String("Rational"));
}



int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    ARIADNE_TEST_CLASS(Rational,TestRational());

    return ARIADNE_TEST_FAILURES;
}
