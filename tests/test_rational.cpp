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
#include "numeric/detail/rational_gmp_bridge.hpp"
#include "numeric/builtin.hpp"
#include "numeric/integer.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/number.hpp"
#include "numeric/floatdp.hpp"
#include "numeric/floatmp.hpp"
#include "foundation/logical.hpp"

#include <array>
#include <algorithm>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <random>
#include <utility>
#include <vector>
#include <sstream>

#include "macro/testing.hpp"

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
    void test_native_contract();
    void test_native_double();
    void test_gmp_differential();

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
    ARIADNE_TEST_CALL(test_native_contract());
    ARIADNE_TEST_CALL(test_native_double());
    ARIADNE_TEST_CALL(test_gmp_differential());
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
    Rational from_raw=rational_from_gmp(raw);
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
    big_output << integer_from_gmp(big_integer);
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

void TestRational::test_native_contract() {
    ARIADNE_TEST_EQUALS(Rational(6,-8).numerator(),Integer(-3));
    ARIADNE_TEST_EQUALS(Rational(6,-8).denominator(),Natural(4u));
    ARIADNE_TEST_EQUALS(Rational(0,-8).denominator(),Natural(1u));
    ARIADNE_TEST_EQUALS(Rational(-6,-8),Rational(3,4));
    ARIADNE_TEST_EQUALS(Rational(std::numeric_limits<std::uint64_t>::max()).numerator(),
                        Integer(String("18446744073709551615")));
    ARIADNE_TEST_EQUALS(Rational(String("  -12 / -18  ")),Rational(2,3));
    ARIADNE_TEST_EQUALS(Rational(String("+0042")),Rational(42));
    for (const char* text : {"", "/2", "2/", "1/2/3", "--1", "1.5"}) {
        ARIADNE_TEST_FAIL(Rational(String(text)));
    }
    ARIADNE_TEST_ASSERT(is_nan(Rational(String("0/0"))));
    ARIADNE_TEST_EQUALS(Rational(String("-27/0")),Rational::inf(Sign::NEGATIVE));
    ARIADNE_TEST_EQUALS(Rational(5)/Rational(0),Rational::inf());
    ARIADNE_TEST_EQUALS(Rational(-5)/Rational(0),Rational::inf(Sign::NEGATIVE));
    ARIADNE_TEST_ASSERT(is_nan(Rational(0)/Rational(0)));
    ARIADNE_TEST_EQUALS(rec(Rational(0)),Rational::inf());
    ARIADNE_TEST_EQUALS(rec(Rational::inf()),Rational(0));
    ARIADNE_TEST_ASSERT(is_nan(Rational::inf()/Rational::inf()));
    ARIADNE_TEST_EQUALS(cmp(Rational::inf(),ExactDouble::inf()),Comparison::EQUAL);
    ARIADNE_TEST_EQUALS(cmp(-Rational::inf(),-ExactDouble::inf()),Comparison::EQUAL);
    ARIADNE_TEST_EQUALS(cmp(Rational::nan(),ExactDouble::inf()),Comparison::INCOMPARABLE);
    const MultiplePrecision precision(128u);
    const FloatMP infinite_float(ExactDouble::inf(),precision);
    const FloatMP finite(ExactDouble(1.0),precision);
    ARIADNE_TEST_EQUALS(cmp(infinite_float,Rational::inf()),Comparison::EQUAL);
    ARIADNE_TEST_EQUALS(cmp(Rational::inf(),infinite_float),Comparison::EQUAL);
    ARIADNE_TEST_EQUALS(cmp(finite,Rational::inf()),Comparison::LESS);
    ARIADNE_TEST_EQUALS(cmp(Rational::nan(),finite),Comparison::INCOMPARABLE);
    ARIADNE_TEST_FAIL(static_cast<void>(GmpRationalBoundary(Rational::nan())));

    // R8: exact cross-products of one-limb magnitudes, including
    // products which exceed 64 bits. Check the independent GMP order.
    {
        const std::uint64_t limit=std::numeric_limits<std::uint64_t>::max();
        const std::array<Rational,4> samples={
            Rational(Integer(limit-17u),Integer(limit-7u)),
            Rational(Integer(limit-19u),Integer(limit-5u)),
            Rational(Integer(limit-13u),Integer(limit-3u)),
            Rational(Integer(limit-11u),Integer(limit-2u))
        };
        for (Rational const& x:samples) {
            for (Rational const& y:samples) {
                GmpRationalBoundary gx(x),gy(y);
                const int raw=mpq_cmp(gx.value(),gy.value());
                const Comparison expected=raw<0 ? Comparison::LESS
                    : raw>0 ? Comparison::GREATER : Comparison::EQUAL;
                ARIADNE_TEST_EQUALS(cmp(x,y),expected);
                ARIADNE_TEST_EQUALS(cmp(-x,-y),
                    raw>0 ? Comparison::LESS
                    : raw<0 ? Comparison::GREATER : Comparison::EQUAL);
            }
        }
    }

    // R28: independent GMP add/sub checks across the 256-bit heap
    // boundary, four effective sign combinations and 8/9-limb cases.
    for (unsigned bits:{127u,128u,255u,256u,257u,319u,320u,448u,512u}) {
        const Integer base=pow(Integer(2),Nat(bits));
        for (unsigned mode=0u;mode<4u;++mode) {
            const Integer a_num=(mode&1u)!=0u ? -(base+13) : base+13;
            const Integer b_num=(mode&2u)!=0u ? -(base+31) : base+31;
            const Rational a(a_num,base*2+17);
            const Rational b(b_num,base*3+19);
            const GmpRationalBoundary oracle_a(a),oracle_b(b);
            GmpRationalBoundary oracle_sum,oracle_diff;
            mpq_add(oracle_sum.mutable_value(),oracle_a.value(),oracle_b.value());
            mpq_sub(oracle_diff.mutable_value(),oracle_a.value(),oracle_b.value());
            const GmpRationalBoundary actual_sum(a+b),actual_diff(a-b);
            ARIADNE_TEST_ASSERT(mpq_equal(actual_sum.value(),oracle_sum.value())!=0);
            ARIADNE_TEST_ASSERT(mpq_equal(actual_diff.value(),oracle_diff.value())!=0);
            ARIADNE_TEST_EQUALS(a-a,Rational(0));
            ARIADNE_TEST_EQUALS(b-b,Rational(0));
        }
    }
    // R30: sparse radix construction with independent GMP oracle.
    // Test high-ratio reduction, exact multiples, nonzero/coprime
    // differences, overflow fallbacks and 4/5 and 128/129 limbs.
    for (unsigned bits:{64u,128u,255u,256u,320u,512u,1024u,2048u,4096u,8192u}) {
        const Integer base=pow(Integer(2),Nat(bits));
        const std::array<std::pair<unsigned,unsigned>,7u> high_pairs={{
            {1u,2u},{1u,3u},{2u,4u},{3u,9u},
            {3u,3u},{2u,3u},{1u,1u}
        }};
        const std::array<std::pair<std::uint64_t,std::uint64_t>,6u> low_pairs={{
            {13u,17u},{63u,126u},{1u,0u},{7u,14u},
            {0u,0u},{std::numeric_limits<std::uint64_t>::max(),
                      std::numeric_limits<std::uint64_t>::max()}
        }};
        for (const auto& [high_n,high_d]:high_pairs) {
            for (const auto& [low_n,low_d]:low_pairs) {
                for (unsigned signs=0u;signs<4u;++signs) {
                    const Integer positive_n=base*Integer(high_n)+Integer(low_n);
                    const Integer positive_d=base*Integer(high_d)+Integer(low_d);
                    const Integer raw_n=(signs&1u)!=0u ? -positive_n : positive_n;
                    const Integer raw_d=(signs&2u)!=0u ? -positive_d : positive_d;
                    const Rational fraction_result(raw_n,raw_d);
                    GmpRationalBoundary expected;
                    integer_to_gmp(mpq_numref(expected.mutable_value()),raw_n);
                    integer_to_gmp(mpq_denref(expected.mutable_value()),raw_d);
                    mpq_canonicalize(expected.mutable_value());
                    const GmpRationalBoundary actual(fraction_result);
                    ARIADNE_TEST_ASSERT(mpq_equal(actual.value(),expected.value())!=0);
                    GmpIntegerBoundary common;
                    mpz_gcd(common.mutable_value(),
                            mpq_numref(actual.value()),mpq_denref(actual.value()));
                    ARIADNE_TEST_ASSERT(mpz_cmp_ui(common.value(),1u)==0);
                    ARIADNE_TEST_ASSERT(mpz_sgn(mpq_denref(actual.value()))>0);
                }
            }
        }
        if (bits>=128u) {
            // A nonzero intermediate word disables the sparse
            // reduction, preserving the original Euclidean path.
            const Integer dense_n=base+pow(Integer(2),Nat(64))+Integer(11);
            const Integer dense_d=base*2+pow(Integer(2),Nat(64))+Integer(17);
            const Rational dense_result(dense_n,dense_d);
            GmpRationalBoundary dense_expected;
            integer_to_gmp(mpq_numref(dense_expected.mutable_value()),dense_n);
            integer_to_gmp(mpq_denref(dense_expected.mutable_value()),dense_d);
            mpq_canonicalize(dense_expected.mutable_value());
            const GmpRationalBoundary dense_actual(dense_result);
            ARIADNE_TEST_ASSERT(mpq_equal(dense_actual.value(),dense_expected.value())!=0);
        }
    }
    // R31: verify GCD shortcuts through real add, subtract and
    // multiply operations against an independent GMP rational oracle.
    // The 2:3 denominator high-word ratio is deliberately *not*
    // an exact integer ratio. Exercise determinant 0, +/-1 and
    // composite values, both operand signs, both orders, and
    // the 4/5/8/9-limb normalization boundaries.
    for (unsigned bits:{128u,192u,255u,256u,320u,448u,512u,
                        1024u,2048u,4096u}) {
        const Integer base=pow(Integer(2),Nat(bits));
        for (unsigned i=0u;i<24u;++i) {
            const Integer an=base+Integer(13u+2u*i);
            const Integer ad=base*2+Integer(17u+3u*i);
            const Integer bn=base+Integer(31u+5u*i);
            const Integer bd=base*3+Integer(19u+7u*i);
            for (unsigned signs=0u;signs<4u;++signs) {
                const Rational left((signs&1u)!=0u ? -an : an,ad);
                const Rational right((signs&2u)!=0u ? -bn : bn,bd);
                const GmpRationalBoundary raw_left(left),raw_right(right);
                GmpRationalBoundary expected_sum,expected_sub,expected_mul,expected_div;
                mpq_add(expected_sum.mutable_value(),
                        raw_left.value(),raw_right.value());
                mpq_sub(expected_sub.mutable_value(),
                        raw_left.value(),raw_right.value());
                mpq_mul(expected_mul.mutable_value(),
                        raw_left.value(),raw_right.value());
                mpq_div(expected_div.mutable_value(),
                        raw_left.value(),raw_right.value());
                const GmpRationalBoundary actual_sum(left+right);
                const GmpRationalBoundary actual_sub(left-right);
                const GmpRationalBoundary actual_mul(left*right);
                const GmpRationalBoundary actual_div(left/right);
                ARIADNE_TEST_ASSERT(mpq_equal(
                    actual_sum.value(),expected_sum.value())!=0);
                ARIADNE_TEST_ASSERT(mpq_equal(
                    actual_sub.value(),expected_sub.value())!=0);
                ARIADNE_TEST_ASSERT(mpq_equal(
                    actual_mul.value(),expected_mul.value())!=0);
                ARIADNE_TEST_ASSERT(mpq_equal(
                    actual_div.value(),expected_div.value())!=0);
                Rational inplace(left);
                inplace+=right;
                const GmpRationalBoundary actual_inplace_add(inplace);
                ARIADNE_TEST_ASSERT(mpq_equal(
                    actual_inplace_add.value(),expected_sum.value())!=0);
                inplace=left;
                inplace*=right;
                const GmpRationalBoundary actual_inplace_mul(inplace);
                ARIADNE_TEST_ASSERT(mpq_equal(
                    actual_inplace_mul.value(),expected_mul.value())!=0);
                inplace=left;
                inplace/=right;
                const GmpRationalBoundary actual_inplace_div(inplace);
                ARIADNE_TEST_ASSERT(mpq_equal(
                    actual_inplace_div.value(),expected_div.value())!=0);
            }
        }
    }
    // R17: high-prefix product bounds must agree with exact GMP
    // comparisons. Cover readily separable and nearly equal ratios,
    // reverse ordering, negatives, and the 8/9-limb stack boundary.
    for (unsigned bits:{128u,256u,448u,512u,1024u,2048u,4096u}) {
        const Integer base=pow(Integer(2),Nat(bits));
        const std::array<Rational,6u> probes={
            Rational(base+13,base*2+17),
            Rational(base+31,base*3+19),
            Rational(base+13,base*2+19),
            Rational(base+1,base+2),
            Rational(base+2,base+3),
            Rational(-(base+29),base*3+7)
        };
        for (Rational const& a:probes) {
            for (Rational const& b:probes) {
                GmpRationalBoundary ga(a),gb(b);
                const int raw=mpq_cmp(ga.value(),gb.value());
                const Comparison expected=raw<0 ? Comparison::LESS
                    : raw>0 ? Comparison::GREATER : Comparison::EQUAL;
                const Comparison reverse=raw>0 ? Comparison::LESS
                    : raw<0 ? Comparison::GREATER : Comparison::EQUAL;
                ARIADNE_TEST_EQUALS(cmp(a,b),expected);
                ARIADNE_TEST_EQUALS(cmp(b,a),reverse);
                ARIADNE_TEST_EQUALS(cmp(-a,-b),reverse);
            }
        }
    }

    // R4: check proper-fraction floor/ceil shortcuts, integral and
    // improper fractions against independent GMP quotients. Include
    // size boundaries around inline storage and long magnitudes.
    for (unsigned bits : {53u,128u,256u,1024u,2048u,4096u}) {
        const Integer base=pow(Integer(2),Nat(bits));
        const Integer below=base-1;
        const Integer above=base+1;
        const std::array<Rational,9> cases={
            Rational(0),Rational(1,3),Rational(-1,3),
            Rational(below,base),Rational(-below,base),
            Rational(base,base),Rational(-base,base),
            Rational(above,base),Rational(-above,base)
        };
        for (Rational const& value : cases) {
            GmpRationalBoundary raw(value);
            GmpIntegerBoundary expected;
            mpz_fdiv_q(expected.mutable_value(),
                       mpq_numref(raw.value()),mpq_denref(raw.value()));
            ARIADNE_TEST_EQUALS(floor(value),expected.to_integer());
            mpz_cdiv_q(expected.mutable_value(),
                       mpq_numref(raw.value()),mpq_denref(raw.value()));
            ARIADNE_TEST_EQUALS(ceil(value),expected.to_integer());
        }
    }

    // R3: exercise stack products (up to 8 limbs each), their exact
    // 8/9-limb boundary, the fallback, signs and opposite orderings.
    for (unsigned bits : {53u,128u,256u,448u,512u,1024u}) {
        const Integer base=pow(Integer(2),Nat(bits));
        const Rational a(base+13,base*3+19);
        const Rational b(base+31,base*2+17);
        GmpRationalBoundary ga(a),gb(b);
        const int raw=mpq_cmp(ga.value(),gb.value());
        const Comparison expected=raw<0 ? Comparison::LESS
            : raw>0 ? Comparison::GREATER : Comparison::EQUAL;
        const Comparison reversed=raw>0 ? Comparison::LESS
            : raw<0 ? Comparison::GREATER : Comparison::EQUAL;
        ARIADNE_TEST_EQUALS(cmp(a,b),expected);
        ARIADNE_TEST_EQUALS(cmp(b,a),reversed);
        ARIADNE_TEST_EQUALS(cmp(-a,-b),reversed);
        ARIADNE_TEST_EQUALS(cmp(-b,-a),expected);
        ARIADNE_TEST_EQUALS(cmp(Rational(0),a),Comparison::LESS);
        ARIADNE_TEST_EQUALS(cmp(a,Rational(0)),Comparison::GREATER);
    }
    const Integer factor=pow(Integer(2),Nat(320))+1;
    const Rational a(factor*7,11),b(11,factor*13);
    ARIADNE_TEST_EQUALS(a*b,Rational(7,13));
    ARIADNE_TEST_EQUALS(a/a,Rational(1));
    ARIADNE_TEST_EQUALS(a-a,Rational(0));
    ARIADNE_TEST_EQUALS((a-a).denominator(),Natural(1u));
    ARIADNE_TEST_EQUALS(Rational(1,6)+Rational(1,10),Rational(4,15));
    ARIADNE_TEST_EQUALS(hlf(Rational(-6,7)),Rational(-3,7));
    ARIADNE_TEST_EQUALS(hlf(Rational(-3,7)),Rational(-3,14));
    ARIADNE_TEST_EQUALS(sqr(Rational(-3,7)),Rational(9,49));
    ARIADNE_TEST_EQUALS(fma(Rational(2,3),Rational(3,7),Rational(1,7)),Rational(3,7));
    ARIADNE_TEST_EQUALS(pow(Rational(-2,3),Int(-3)),Rational(-27,8));
    ARIADNE_TEST_EQUALS(pow(Rational(-1),std::numeric_limits<Int>::min()),Rational(1));
    ARIADNE_TEST_EQUALS(pow(Integer(-1),std::numeric_limits<Int>::min()),Rational(1));
    ARIADNE_TEST_EQUALS(pow(Rational::nan(),Nat(0)),Rational(1));
    ARIADNE_TEST_ASSERT(is_nan(pow(Rational::nan(),Nat(3))));

    Rational source(a),moved(std::move(source));
    ARIADNE_TEST_EQUALS(moved,a);
    ARIADNE_TEST_EQUALS(source,Rational(0));
    Rational destination(2,3);
    destination=std::move(moved);
    ARIADNE_TEST_EQUALS(destination,a);
    ARIADNE_TEST_EQUALS(moved,Rational(2,3));
    Rational& alias=destination;
    destination=std::move(alias);
    ARIADNE_TEST_EQUALS(destination,a);
    destination*=destination;
    ARIADNE_TEST_EQUALS(destination,sqr(a));
    destination/=alias;
    ARIADNE_TEST_EQUALS(destination,Rational(1));
    Rational special=Rational::inf(Sign::NEGATIVE);
    Rational moved_special(std::move(special));
    ARIADNE_TEST_EQUALS(moved_special,Rational::inf(Sign::NEGATIVE));
    ARIADNE_TEST_EQUALS(special,Rational(0));
}

void TestRational::test_native_double() {
    const Integer two53=pow(Integer(2),Nat(53));
    const Integer two1074=pow(Integer(2),Nat(1074));
    ARIADNE_TEST_EQUALS(Rational(1,3).get_d(),0x1.5555555555555p-2);
    ARIADNE_TEST_EQUALS(Rational(-1,3).get_d(),-0x1.5555555555555p-2);
    ARIADNE_TEST_EQUALS(Rational(two53+1,two53).get_d(),1.0);
    ARIADNE_TEST_EQUALS(Rational(1,two1074).get_d(),std::numeric_limits<double>::denorm_min());
    ARIADNE_TEST_EQUALS(Rational(3,two1074*2).get_d(),std::numeric_limits<double>::denorm_min());
    const double underflow=Rational(-1,two1074*2).get_d();
    ARIADNE_TEST_EQUALS(underflow,0.0);
    ARIADNE_TEST_ASSERT(!std::signbit(underflow));
    ARIADNE_TEST_ASSERT(std::isinf(Rational(pow(Integer(2),Nat(1024))).get_d()));
    const std::array<double,9> fixed={0.0,-0.0,1.0,-1.0,0.1,
        std::numeric_limits<double>::denorm_min(),std::numeric_limits<double>::min(),
        std::numeric_limits<double>::max(),-std::numeric_limits<double>::max()};
    for (double value : fixed) {
        ARIADNE_TEST_EQUALS(Rational(ExactDouble(value)).get_d(),value);
    }
    std::mt19937_64 random(0x726174696f6e616cULL);
    bool roundtrips=true;
    for (unsigned int i=0u; i!=2000u; ++i) {
        const double value=std::bit_cast<double>(random());
        if (!std::isfinite(value)) { continue; }
        roundtrips=roundtrips && Rational(ExactDouble(value)).get_d()==value;
    }
    ARIADNE_TEST_ASSERT(roundtrips);
}

void TestRational::test_gmp_differential() {
    std::mt19937_64 random(0x6672616374696f6eULL);
    GmpRationalBoundary ga,gb,expected;
    bool operations_match=true,canonical=true,conversion_matches=true;
    auto check=[&](Rational const& value) {
        GmpRationalBoundary actual(value);
        operations_match=operations_match && mpq_equal(actual.value(),expected.value());
        GmpIntegerBoundary common;
        mpz_gcd(common.mutable_value(),mpq_numref(actual.value()),mpq_denref(actual.value()));
        canonical=canonical && mpz_cmp_ui(common.value(),1u)==0 && mpz_sgn(mpq_denref(actual.value()))>0;
    };
    for (unsigned int i=0u; i!=1200u; ++i) {
        const std::size_t count=1u+(i%17u);
        std::vector<std::uint64_t> words(count);
        auto integer=[&](bool nonzero) {
            for (auto& word : words) { word=random(); }
            if (nonzero) { words.back()|=1u; }
            if (!nonzero && i%13u==0u) { std::fill(words.begin(),words.end(),0u); }
            return Integer::from_magnitude(BigUnsigned(words),(random()&1u)!=0u);
        };
        const Integer an=integer(false),ad=integer(true),bn=integer(false),bd=integer(true);
        const Rational a(an,ad),b(bn,bd);
        integer_to_gmp(mpq_numref(ga.mutable_value()),an);
        integer_to_gmp(mpq_denref(ga.mutable_value()),ad);
        mpq_canonicalize(ga.mutable_value());
        integer_to_gmp(mpq_numref(gb.mutable_value()),bn);
        integer_to_gmp(mpq_denref(gb.mutable_value()),bd);
        mpq_canonicalize(gb.mutable_value());
        mpq_set(expected.mutable_value(),ga.value()); check(a);
        GmpIntegerBoundary integral;
        mpz_fdiv_q(integral.mutable_value(),mpq_numref(ga.value()),mpq_denref(ga.value()));
        operations_match=operations_match && floor(a)==integral.to_integer();
        mpz_cdiv_q(integral.mutable_value(),mpq_numref(ga.value()),mpq_denref(ga.value()));
        operations_match=operations_match && ceil(a)==integral.to_integer();
        mpq_add(expected.mutable_value(),ga.value(),gb.value()); check(a+b);
        mpq_sub(expected.mutable_value(),ga.value(),gb.value()); check(a-b);
        mpq_mul(expected.mutable_value(),ga.value(),gb.value()); check(a*b);
        mpq_mul(expected.mutable_value(),ga.value(),ga.value()); check(sqr(a));
        mpq_div_2exp(expected.mutable_value(),ga.value(),1u); check(hlf(a));
        if (!is_zero(b)) { mpq_div(expected.mutable_value(),ga.value(),gb.value()); check(a/b); }
        if (!is_zero(a)) { mpq_inv(expected.mutable_value(),ga.value()); check(rec(a)); }
        const int order=mpq_cmp(ga.value(),gb.value());
        const Comparison want=order<0 ? Comparison::LESS : order>0 ? Comparison::GREATER : Comparison::EQUAL;
        operations_match=operations_match && cmp(a,b)==want && bool(a==b)==bool(mpq_equal(ga.value(),gb.value()));
        conversion_matches=conversion_matches &&
            std::bit_cast<std::uint64_t>(a.get_d())==std::bit_cast<std::uint64_t>(mpq_get_d(ga.value()));
    }
    // R16: explicit coprime-denominator (g=1), shared factor with
    // h=1 and shared factor with h>1. GMP checks the exact result and
    // canonical positive denominator for both signs and both orders.
    for (unsigned bits : {53u,128u,256u,1024u,2048u,4096u}) {
        const Integer base=pow(Integer(2),Nat(bits));
        const std::array<std::pair<Rational,Rational>,7u> inputs={{
            {Rational(1,base),Rational(1,base+1)},
            {Rational(7,base),Rational(-11,base+1)},
            {Rational(-1,base),Rational(1,base+1)},
            {Rational(1,base),Rational(-1,base)},
            {Rational(1,2*base),Rational(1,2*(base+1))},
            {Rational(-1,2*base),Rational(1,2*(base+1))},
            {Rational(1,6),Rational(1,10)}
        }};
        for (auto const& pair:inputs) {
            for (unsigned order=0u;order<2u;++order) {
                const Rational& a=order==0u ? pair.first : pair.second;
                const Rational& b=order==0u ? pair.second : pair.first;
                GmpRationalBoundary native_a(a),native_b(b);
                mpq_add(expected.mutable_value(),native_a.value(),native_b.value());
                check(a+b);
                // R23: replacement of the left/right operand and of a
                // preexisting special destination preserves the result.
                Rational direct=a;
                direct+=b;
                check(direct);
                direct=b;
                direct=a+direct;
                check(direct);
                direct=Rational::nan();
                direct=a+b;
                check(direct);
                mpq_sub(expected.mutable_value(),native_a.value(),native_b.value());
                check(a-b);
                direct=a;
                direct-=b;
                check(direct);
                direct=b;
                direct=a-direct;
                check(direct);
                direct=Rational::inf(Sign::NEGATIVE);
                direct=a-b;
                check(direct);
                Rational aliased=a;
                aliased+=a;
                mpq_add(expected.mutable_value(),native_a.value(),native_a.value());
                check(aliased);
                direct=a;
                Rational& direct_alias=direct;
                direct+=direct_alias;
                check(direct);
                aliased=a;
                aliased-=a;
                mpq_sub(expected.mutable_value(),native_a.value(),native_a.value());
                check(aliased);
                direct-=direct_alias;
                check(direct);
            }
        }
    }
    // R20: exact GMP cross-sum checks near both sides of the 8-limb
    // stack/fallback boundary, with signs, carries and cancellation.
    for (unsigned bits : {53u,128u,255u,256u,320u,447u,
                          448u,511u,512u,1024u}) {
        const Integer base=pow(Integer(2),Nat(bits));
        const std::array<std::pair<Rational,Rational>,4u> pairs={{
            {Rational(base+13,base*2+17),Rational(base+31,base*3+19)},
            {Rational(base-13,base*3+23),Rational(base+19,base*2+31)},
            {Rational(base+1,base+2),Rational(base+2,base+3)},
            {Rational(base+13,base*2+17),Rational(base+13,base*2+17)}
        }};
        for (auto const& pair:pairs) {
            for (unsigned mask=0u;mask<4u;++mask) {
                const Rational lhs=(mask&1u)!=0u ? -pair.first : pair.first;
                const Rational rhs=(mask&2u)!=0u ? -pair.second : pair.second;
                GmpRationalBoundary oracle_lhs(lhs),oracle_rhs(rhs);
                mpq_add(expected.mutable_value(),oracle_lhs.value(),oracle_rhs.value());
                check(lhs+rhs);
                Rational running=lhs;
                running+=rhs;
                check(running);
                mpq_sub(expected.mutable_value(),oracle_lhs.value(),oracle_rhs.value());
                check(lhs-rhs);
                running=lhs;
                running-=rhs;
                check(running);
            }
        }
    }
    // R19: all four combinations of identity/nonidentity cross-GCD
    // divisors for both multiply and divide, across the inline boundary
    // and long operands. Check exact GMP arithmetic, canonicality,
    // in-place operations, swapped order and independent sign choices.
    for (unsigned bits : {53u,128u,256u,1024u,2048u,4096u}) {
        const Integer base=pow(Integer(2),Nat(bits));
        const Rational a(base+3,base+5);
        const std::array<Rational,4u> factors={
            Rational(base+9,base+7), // both quotient divisors are one
            Rational(base+9,base+3), // shared numerator/denominator
            Rational(base+5,base+9), // opposite shared cross factor
            Rational(base+5,base+3)  // two nontrivial cross factors
        };
        for (Rational const& factor : factors) {
            for (unsigned signs=0u;signs<4u;++signs) {
                const Rational lhs=(signs&1u)!=0u ? -a : a;
                const Rational rhs=(signs&2u)!=0u ? -factor : factor;
                GmpRationalBoundary oracle_lhs(lhs),oracle_rhs(rhs);
                mpq_mul(expected.mutable_value(),oracle_lhs.value(),oracle_rhs.value());
                check(lhs*rhs);
                Rational result=lhs;
                result*=rhs;
                check(result);
                mpq_mul(expected.mutable_value(),oracle_rhs.value(),oracle_lhs.value());
                check(rhs*lhs);
                mpq_div(expected.mutable_value(),oracle_lhs.value(),oracle_rhs.value());
                check(lhs/rhs);
                result=lhs;
                result/=rhs;
                check(result);
                mpq_div(expected.mutable_value(),oracle_rhs.value(),oracle_lhs.value());
                check(rhs/lhs);
            }
        }
        // Preserve the existing alias-aware square and self-division.
        GmpRationalBoundary oracle_self(a);
        mpq_mul(expected.mutable_value(),oracle_self.value(),oracle_self.value());
        Rational self=a;
        self*=self;
        check(self);
        mpq_div(expected.mutable_value(),oracle_self.value(),oracle_self.value());
        self=a;
        Rational& self_alias=self;
        self/=self_alias;
        check(self);
    }
    // R22: scalar exact cancellation, including the even non-power-of-two
    // divisor 14 and both sides of the eight-limb factor-buffer bound.
    for (unsigned bits : {53u,127u,128u,255u,256u,447u,448u,511u,512u,1024u}) {
        const Integer base=pow(Integer(2),Nat(bits));
        const Rational a(15*(base+1),14*base);
        const std::array<Rational,4u> factors={Rational(1),Rational(1,15),
                                             Rational(14),Rational(14,15)};
        ARIADNE_TEST_ASSERT(a*factors[3u]==Rational(base+1,base));
        ARIADNE_TEST_ASSERT(a/Rational(15,14)==Rational(base+1,base));
        for (auto const& factor:factors) {
            for (unsigned mask=0u;mask<4u;++mask) {
                const Rational lhs=(mask&1u)!=0u ? -a : a;
                const Rational rhs=(mask&2u)!=0u ? -factor : factor;
                GmpRationalBoundary oracle_lhs(lhs),oracle_rhs(rhs);
                mpq_mul(expected.mutable_value(),oracle_lhs.value(),oracle_rhs.value());
                check(lhs*rhs);
                Rational running=lhs;
                running*=rhs;
                check(running);
                mpq_div(expected.mutable_value(),oracle_lhs.value(),oracle_rhs.value());
                check(lhs/rhs);
                running=lhs;
                running/=rhs;
                check(running);
            }
        }
    }
    // Independently check truncation at normal/subnormal and overflow edges,
    // under all four host rounding modes. Operands are not formed as doubles.
    const int saved_rounding=std::fegetround();
    const std::array<int,4> modes={FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO};
    // R5: exact high-prefix conversion, 4-to-5-limb boundary, and
    // near-half/near-one quotients with nonzero low tails. Compare raw
    // binary64 patterns with the independent GMP implementation under
    // every rounding mode, including cases requiring the full fallback.
    const std::array<unsigned int,12> shifts={0u,1u,52u,53u,1022u,1023u,1024u,1073u,1074u,1075u,1076u,2048u};
    for (int mode : modes) {
        ARIADNE_TEST_EQUALS(std::fesetround(mode),0);
        // R12: exercise exact one-limb conversion across wide unsigned
        // ratios, repeating fractions, both signs and four host modes.
        // Each result is compared bit-for-bit against GMP.
        std::mt19937_64 small_random(0x1200530012812026ULL);
        for (unsigned sample=0u;sample<512u;++sample) {
            std::uint64_t numerator=small_random()|std::uint64_t{1};
            std::uint64_t denominator=small_random()|std::uint64_t{1};
            if (sample%4u==0u) {
                numerator=(std::uint64_t{1}<<53u)+small_random()%100000u;
                denominator=2u*numerator+small_random()%100000u;
            } else if (sample%4u==1u) {
                numerator=std::uint64_t{1}<<(small_random()%63u);
                denominator=std::uint64_t{1}<<(small_random()%63u);
            } else if (sample%4u==2u) {
                numerator=1u+small_random()%1024u;
            }
            const Rational positive{Integer(numerator),Integer(denominator)};
            const Rational negative=-positive;
            for (const Rational& value:{positive,negative}) {
                const GmpRationalBoundary oracle(value);
                conversion_matches=conversion_matches &&
                    std::bit_cast<std::uint64_t>(value.get_d())==
                    std::bit_cast<std::uint64_t>(mpq_get_d(oracle.value()));
            }
        }
        // R29: sparse two-word prefixes with arbitrary nonzero
        // low tails. GMP is the bit-for-bit binary64 oracle.
        // High-word coefficients cover exact and nonexact leading
        // divisions under all four active host rounding modes.
        for (unsigned bits:{128u,256u,320u,512u,1024u,2048u}) {
            const Integer base=pow(Integer(2),Nat(bits));
            for (std::uint64_t high_n:{1u,3u,7u}) {
                for (std::uint64_t high_d:{2u,3u,5u}) {
                    const Integer n=base*Integer(high_n)
                        +Integer((random()&65535u)+1u);
                    const Integer d=base*Integer(high_d)
                        +Integer((random()&65535u)+1u);
                    const Rational positive(n,d);
                    for (Rational const& value:
                            {positive,-positive}) {
                        const GmpRationalBoundary oracle(value);
                        conversion_matches=conversion_matches &&
                            std::bit_cast<std::uint64_t>(value.get_d())==
                            std::bit_cast<std::uint64_t>(
                                mpq_get_d(oracle.value()));
                    }
                }
            }
        }
        // R13: independent random dense GMP ratios spanning 5..8
        // limbs, both signs, and fractions close to 1/2 and 1/3.
        // The exact certificate must work independently of the active
        // floating-point rounding mode or fall back without error.
        for (unsigned limbs=5u;limbs<=8u;++limbs) {
            std::vector<std::uint64_t> numerator_limbs(limbs);
            std::vector<std::uint64_t> denominator_limbs(limbs);
            for (unsigned sample=0u;sample<96u;++sample) {
                for (unsigned i=0u;i<limbs;++i) {
                    numerator_limbs[i]=random();
                    denominator_limbs[i]=random();
                }
                numerator_limbs.back()=(random()&0x3u)+1u;
                denominator_limbs.back()=(random()&0x3u)+2u;
                const Integer num=Integer::from_magnitude(
                    BigUnsigned(numerator_limbs),sample%2u!=0u);
                const Integer den=Integer::from_magnitude(
                    BigUnsigned(denominator_limbs),false);
                const Rational value(num,den);
                const GmpRationalBoundary oracle(value);
                conversion_matches=conversion_matches &&
                    std::bit_cast<std::uint64_t>(value.get_d())==
                    std::bit_cast<std::uint64_t>(mpq_get_d(oracle.value()));
            }
        }
        for (unsigned bits : {53u,128u,256u,448u,512u,1024u,2048u,4096u}) {
            const Integer base=pow(Integer(2),Nat(bits));
            const Integer twice=base*2;
            // R7 probes both a near-binary-boundary quotient and
            // repeating quotients (not only ratios near unity).
            const std::array<Rational,12> probes={
                Rational(base+13,twice+17),
                Rational(base-13,twice+17),
                Rational(base+31,twice+19),
                Rational(base+1,base-1),
                Rational(base-1,base+1),
                Rational(twice+1,base+1),
                Rational(-(base+13),twice+17),
                Rational(base+1,Integer(1)),
                Rational(base+31,base*3+19),
                Rational(base-1,base*3+17),
                Rational(-(base+13),base*3+19),
                Rational(base*3+1,twice+17)
            };
            for (Rational const& value : probes) {
                const GmpRationalBoundary oracle(value);
                conversion_matches=conversion_matches &&
                    std::bit_cast<std::uint64_t>(value.get_d())==
                    std::bit_cast<std::uint64_t>(mpq_get_d(oracle.value()));
            }
        }
        for (unsigned int shift : shifts) {
            const Integer factor=pow(Integer(2),Nat(shift));
            const Integer n=Integer(random()|1u),d=Integer(random()|1u);
            const std::array<Rational,4> values={Rational(n,d*factor),Rational(-n,d*factor),
                                                Rational(n*factor,d),Rational(-n*factor,d)};
            for (Rational const& value : values) {
                GmpRationalBoundary raw(value);
                const auto actual=std::bit_cast<std::uint64_t>(value.get_d());
                const auto expected_bits=std::bit_cast<std::uint64_t>(mpq_get_d(raw.value()));
                if (actual!=expected_bits) {
                    std::cerr << "Rational double mismatch: mode=" << mode << " shift=" << shift
                              << " actual=" << actual << " expected=" << expected_bits << '\n';
                    conversion_matches=false;
                }
            }
        }
    }
    ARIADNE_TEST_EQUALS(std::fesetround(saved_rounding),0);
    ARIADNE_TEST_ASSERT(operations_match);
    ARIADNE_TEST_ASSERT(canonical);
    ARIADNE_TEST_ASSERT(conversion_matches);
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
