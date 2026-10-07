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
#include "numeric/number.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"
#include "foundation/logical.hpp"

#include <iomanip>
#include <sstream>

#include "utility/test.hpp"

using namespace std;
using namespace Ariadne;

namespace {

constexpr bool check_concept()
{
    return requires(unsigned int m, unsigned long int lm, int n, long int ln, Integer z, Dyadic w, Rational q, Rational q2, Boolean b) {
        q=Rational(); q=Rational(m); q=Rational(lm); q=Rational(n); q=Rational(ln); q=Rational(z); q=Rational(q);
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
    };
}

static_assert(check_concept());

} // namespace



class TestRational
{
  public:
    void test();
  private:
    void test_literal();
    void test_conversions();
    void test_arithmetic();
    void test_rounding();
    void test_comparisons();
    void test_infinity();
    void test_bounds();

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
}

void TestRational::test_literal() {
    ARIADNE_TEST_CONSTRUCT(Rational,q,(3.25_q));
    ARIADNE_TEST_EQUALS(q,Rational(13,4));
    ARIADNE_TEST_EQUALS(3.25_q,Rational(13,4));
    ARIADNE_TEST_EQUALS(-11.375_q,Rational(-91,8));
    ARIADNE_TEST_EQUALS(10.3_q,Rational(103,10));
    ARIADNE_TEST_EQUALS(0.333333333333333333_q,Rational(1,3));
    ARIADNE_TEST_EQUALS(0.2857142857142857_q,Rational(2,7));
    ARIADNE_TEST_EQUALS(operator""_q(-0.333333333333333333L),Rational(-1,3));
    ARIADNE_TEST_FAIL(operator""_q(1.0000000002L));
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
    ARIADNE_TEST_ASSERT(is_nan(Rational(0,0)));

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

    std::istringstream empty_input;
    Rational unchanged(5,7);
    empty_input >> unchanged;
    ARIADNE_TEST_EQUAL(unchanged,Rational(5,7));

    String big_digits(600,'9');
    mpz_t big_integer;
    mpz_init(big_integer);
    ARIADNE_TEST_EQUALS(mpz_set_str(big_integer,big_digits.c_str(),10),0);
    std::ostringstream big_output;
    write(big_output,big_integer);
    mpz_clear(big_integer);
    ARIADNE_TEST_EQUALS(big_output.str(),std::string(600,'9'));
}

void TestRational::test_arithmetic() {
    ARIADNE_TEST_EQUAL(Rational(4,5)+Rational(-2,7),Rational(18,35));
    ARIADNE_TEST_EQUAL(Rational(-4,5)-Rational(-2,7),Rational(-18,35));
    ARIADNE_TEST_EQUAL(Rational(4,5)*Rational(-2,7),Rational(-8,35));
    ARIADNE_TEST_EQUAL(Rational(4,5)/Rational(-2,7),Rational(-14,5));
    ARIADNE_TEST_EQUAL(div(Integer(3),Integer(4)),Rational(3,4));
    ARIADNE_TEST_EQUAL(pow(Integer(2),Int(-3)),Rational(1,8));

    Rational accumulated(3,4);
    ARIADNE_TEST_EQUAL((accumulated+=Rational(1,4)),Rational(1));
    ARIADNE_TEST_EQUAL((accumulated-=Rational(1,2)),Rational(1,2));
    ARIADNE_TEST_EQUAL((accumulated*=Rational(4)),Rational(2));
    ARIADNE_TEST_EQUAL((accumulated/=Rational(8)),Rational(1,4));

    PositiveRational p2=cast_positive(Rational(2));
    PositiveRational p3=cast_positive(Rational(3));
    ARIADNE_TEST_FAIL(cast_positive(Rational(-1)));
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

    ARIADNE_TEST_FAIL(round(Rational::inf()));
    ARIADNE_TEST_FAIL(ceil(Rational::inf()));
    ARIADNE_TEST_FAIL(floor(Rational::inf()));
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
    ARIADNE_TEST_FAIL(cmp(Rational(1),ExactDouble(std::numeric_limits<double>::quiet_NaN())));
}

void TestRational::test_infinity() {
    Rational qinf=Rational::inf();
    Rational qninf=Rational::inf(Sign::NEGATIVE);
    Rational qnan=Rational::nan();

    ARIADNE_TEST_ASSERT(is_inf(Rational(ExactDouble(std::numeric_limits<double>::infinity()))));
    ARIADNE_TEST_ASSERT(is_inf(Rational(ExactDouble(-std::numeric_limits<double>::infinity()))));
    ARIADNE_TEST_ASSERT(is_nan(Rational(ExactDouble(std::numeric_limits<double>::quiet_NaN()))));

    ARIADNE_TEST_ASSERT(is_nan(Rational::nan()));
    ARIADNE_TEST_ASSERT(not is_inf(Rational::nan()));
    ARIADNE_TEST_ASSERT(not is_zero(Rational(1)));
    ARIADNE_TEST_ASSERT(not is_zero(Rational::nan()));
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

    check_bounds(p.pm(1),1,4);
    ARIADNE_TEST_EQUALS(p.lower(),Rational(2));
    ARIADNE_TEST_EQUALS(p.upper(),Rational(3));
    check_bounds(+p,2,3);
    check_bounds(-p,-3,-2);
    check_bounds(p+n,-1,1);
    check_bounds(p-n,4,6);

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
    check_bounds(pow(p,Int(2)),4,9);
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

    RationalBounds left(0,1);
    RationalBounds right(2,3);
    RationalBounds containing(-1,4);
    RationalBounds crossing(1,2);
    ARIADNE_TEST_ASSERT(inconsistent(right,left));
    ARIADNE_TEST_ASSERT(inconsistent(left,right));
    ARIADNE_TEST_ASSERT(not inconsistent(left,crossing));

    ARIADNE_TEST_ASSERT(refines(left,containing));
    ARIADNE_TEST_ASSERT(not refines(containing,left));
    ARIADNE_TEST_ASSERT(not refines(RationalBounds(-2,2),RationalBounds(-2,1)));

    check_bounds(refinement(RationalBounds(0,3),RationalBounds(1,4)),1,3);
    check_bounds(coarsening(RationalBounds(0,3),RationalBounds(1,4)),0,4);

    std::ostringstream bounds_stream;
    bounds_stream << RationalBounds(1,2);
    ARIADNE_TEST_EQUALS(bounds_stream.str(),std::string("[1:2]"));

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
    ARIADNE_TEST_FAIL(pow(dmb,Int(-1)));
}





int main() {
    std::cout<<std::setprecision(20);
    std::cerr<<std::setprecision(20);

    ARIADNE_TEST_CLASS(Rational,TestRational());

    return ARIADNE_TEST_FAILURES;
}
