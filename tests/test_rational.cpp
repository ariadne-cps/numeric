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

#include <array>
#include <algorithm>
#include <bit>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <string>
#include <limits>
#include <random>
#include <utility>
#include <vector>
#include <sstream>

#include "macro/testing.hpp"

using namespace std;
using namespace Ariadne;

namespace {

#ifndef ARIADNE_TEST_DATA_DIR
#error "Provide the test data directory from CMake"
#endif
inline constexpr char R8_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_r8.tsv";
inline constexpr char R28_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_r28.tsv";
inline constexpr char R4_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_r4.tsv";
inline constexpr char R16_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_r16.tsv";
inline constexpr char R20_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_r20.tsv";
inline constexpr char R19_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_r19.tsv";
inline constexpr char R22_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_r22.tsv";
inline constexpr char R0_RANDOM_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_random.tsv";
inline constexpr char R30_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_r30.tsv";
inline constexpr char R17_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_r17.tsv";
inline constexpr char R3_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_r3.tsv";
inline constexpr char R7_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_r7.tsv";
inline constexpr char R12_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_r12.tsv";
inline constexpr char R29_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_r29.tsv";
inline constexpr char R5_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/rational_r5.tsv";

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
    void test_differential();

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
    ARIADNE_TEST_CALL(test_differential());
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
    // Exact FloatMP-to-Rational conversion must not depend on another floating-point backend.
    const MultiplePrecision ap_precision(256_bits);
    const FloatMP mp(Dyadic(5,2u),ap_precision);
    ARIADNE_TEST_EQUAL(Rational(mp),Rational(5,4));
    const Integer wide_numerator(String("123456789012345678901234567890123456789"));
    const Dyadic wide_dyadic(wide_numerator,73u);
    const FloatMP wide_ap(wide_dyadic,ap_precision);
    ARIADNE_TEST_EQUAL(Rational(wide_ap),Rational(wide_dyadic));
    ARIADNE_TEST_EQUAL(static_cast<Dyadic>(wide_ap),wide_dyadic);
    ARIADNE_TEST_EQUAL(cast_integer(FloatMP(Integer(-17),ap_precision)),Integer(-17));
    ARIADNE_TEST_FAIL(cast_integer(mp));
    ARIADNE_TEST_ASSERT(is_nan(Rational(FloatMP::nan(ap_precision))));
    ARIADNE_TEST_EQUAL(Rational(FloatMP::inf(Sign::POSITIVE,ap_precision)),
                       Rational::inf(Sign::POSITIVE));
    ARIADNE_TEST_EQUAL(Rational(FloatMP::inf(Sign::NEGATIVE,ap_precision)),
                       Rational::inf(Sign::NEGATIVE));

    ARIADNE_TEST_EQUAL(Rational(String("7/9")),Rational(7,9));
    ARIADNE_TEST_FAIL(Rational(String("not-a-rational")));
    ARIADNE_TEST_ASSERT(is_nan(Rational(0,0)));

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
    std::ostringstream big_output;
    big_output << Integer(big_digits);
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
    const MultiplePrecision precision(128_bits);
    const FloatMP infinite_float(ExactDouble::inf(),precision);
    const FloatMP finite(ExactDouble(1.0),precision);
    ARIADNE_TEST_EQUALS(cmp(infinite_float,Rational::inf()),Comparison::EQUAL);
    ARIADNE_TEST_EQUALS(cmp(Rational::inf(),infinite_float),Comparison::EQUAL);
    ARIADNE_TEST_EQUALS(cmp(finite,Rational::inf()),Comparison::LESS);
    ARIADNE_TEST_EQUALS(cmp(Rational::nan(),finite),Comparison::INCOMPARABLE);

    // R8: original 4x4 near-2^64 exact cross-product ordering.
    // Freeze independent expected signs; test reversed signed ordering too.
    {
        const std::uint64_t limit=std::numeric_limits<std::uint64_t>::max();
        const std::array<Rational,4u> samples={
            Rational(Integer(limit-17u),Integer(limit-7u)),
            Rational(Integer(limit-19u),Integer(limit-5u)),
            Rational(Integer(limit-13u),Integer(limit-3u)),
            Rational(Integer(limit-11u),Integer(limit-2u))
        };
        std::ifstream reference(R8_ORACLE_PATH);
        if (!reference) throw std::runtime_error("Missing R8 oracle TSV");
        std::string header;
        while (reference.peek()=='#') std::getline(reference,header);
        for (unsigned i=0u;i<4u;++i) {
            for (unsigned j=0u;j<4u;++j) {
                unsigned ref_i=0u,ref_j=0u;
                int order=2;
                if (!(reference >> ref_i >> ref_j >> order) ||
                    ref_i!=i || ref_j!=j || order < -1 || order > 1) {
                    throw std::runtime_error("Malformed or unordered R8 TSV row");
                }
                const Comparison expected=order<0 ? Comparison::LESS
                    : order>0 ? Comparison::GREATER : Comparison::EQUAL;
                ARIADNE_TEST_EQUALS(cmp(samples[i],samples[j]),expected);
                const Comparison reversed=order>0 ? Comparison::LESS
                    : order<0 ? Comparison::GREATER : Comparison::EQUAL;
                ARIADNE_TEST_EQUALS(cmp(-samples[i],-samples[j]),reversed);
            }
        }
        std::string extra;
        if (reference >> extra)
            throw std::runtime_error("Extra R8 oracle TSV rows");
    }

    // R28: all original heap thresholds and four signed combinations.
    // Source rational operands and both canonical results are frozen
    // independently, not computed with GMP during the test.
    {
        std::ifstream reference(R28_ORACLE_PATH);
        if (!reference) throw std::runtime_error("Missing R28 oracle TSV");
        std::string header;
        while (reference.peek()=='#') std::getline(reference,header);
        for (unsigned bits:{127u,128u,255u,256u,257u,319u,320u,448u,512u}) {
            for (unsigned mode=0u;mode<4u;++mode) {
                unsigned ref_bits=0u,ref_mode=0u;
                std::string an,ad,bn,bd,sn,sd,dn,dd;
                if (!(reference >> ref_bits >> ref_mode >> an >> ad
                      >> bn >> bd >> sn >> sd >> dn >> dd) ||
                    ref_bits!=bits || ref_mode!=mode) {
                    throw std::runtime_error("Malformed or unordered R28 TSV row");
                }
                const Rational a{Integer(String(an)),Integer(String(ad))};
                const Rational b{Integer(String(bn)),Integer(String(bd))};
                const Rational sum{Integer(String(sn)),Integer(String(sd))};
                const Rational difference{Integer(String(dn)),Integer(String(dd))};
                ARIADNE_TEST_EQUALS(a+b,sum);
                ARIADNE_TEST_EQUALS(a-b,difference);
                ARIADNE_TEST_EQUALS(a-a,Rational(0));
                ARIADNE_TEST_EQUALS(b-b,Rational(0));
            }
        }
        std::string extra;
        if (reference >> extra)
            throw std::runtime_error("Extra R28 oracle TSV rows");
    }

    // R30: original sparse and dense construction cases, including the
    // 128/129-limb boundary. Expected canonical components are frozen
    // independently, without GMP or Rational arithmetic at test time.
    std::ifstream r30_golden(R30_ORACLE_PATH);
    if (!r30_golden) throw std::runtime_error("Missing R30 oracle TSV");
    std::string r30_header;
    while (r30_golden.peek()=='#') std::getline(r30_golden,r30_header);
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
        for (unsigned high_index=0u;high_index<high_pairs.size();++high_index) {
            const auto& [high_n,high_d]=high_pairs[high_index];
            for (unsigned low_index=0u;low_index<low_pairs.size();++low_index) {
                const auto& [low_n,low_d]=low_pairs[low_index];
                for (unsigned signs=0u;signs<4u;++signs) {
                    const Integer positive_n=base*Integer(high_n)+Integer(low_n);
                    const Integer positive_d=base*Integer(high_d)+Integer(low_d);
                    const Integer raw_n=(signs&1u)!=0u ? -positive_n : positive_n;
                    const Integer raw_d=(signs&2u)!=0u ? -positive_d : positive_d;
                    const Rational fraction_result(raw_n,raw_d);
                    unsigned ref_bits=0u,ref_high=0u,ref_low=0u,ref_signs=0u;
                    std::string ref_num,ref_den;
                    if (!(r30_golden >> ref_bits >> ref_high >> ref_low >> ref_signs
                          >> ref_num >> ref_den) ||
                        ref_bits!=bits || ref_high!=high_index ||
                        ref_low!=low_index || ref_signs!=signs) {
                        throw std::runtime_error("Missing or unordered R30 TSV row");
                    }
                    ARIADNE_TEST_EQUALS(fraction_result.numerator(),Integer(String(ref_num)));
                    ARIADNE_TEST_EQUALS(fraction_result.denominator(),Integer(String(ref_den)));
                }
            }
        }
        if (bits>=128u) {
            const Integer dense_n=base+pow(Integer(2),Nat(64))+Integer(11);
            const Integer dense_d=base*2+pow(Integer(2),Nat(64))+Integer(17);
            const Rational dense_result(dense_n,dense_d);
            unsigned ref_bits=0u,ref_high=0u,ref_low=0u,ref_signs=0u;
            std::string ref_num,ref_den;
            if (!(r30_golden >> ref_bits >> ref_high >> ref_low >> ref_signs
                  >> ref_num >> ref_den) ||
                ref_bits!=bits || ref_high!=7u || ref_low!=0u || ref_signs!=0u) {
                throw std::runtime_error("Missing or unordered R30 dense TSV row");
            }
            ARIADNE_TEST_EQUALS(dense_result.numerator(),Integer(String(ref_num)));
            ARIADNE_TEST_EQUALS(dense_result.denominator(),Integer(String(ref_den)));
        }
    }
    std::string r30_extra;
    if (r30_golden >> r30_extra)
        throw std::runtime_error("Extra R30 oracle TSV rows");
    // R31/R32: preserve all 10 * 24 * 4 signed determinant/GCD
    // configurations and every direct/in-place arithmetic check.
    // Both canonical input operands and all four canonical results
    // come from source-keyed exact integer TSV oracles, not GMP.
    for (unsigned bits:{128u,192u,255u,256u,320u,448u,512u,
                        1024u,2048u,4096u}) {
        const Integer base=pow(Integer(2),Nat(bits));
        std::ifstream reference;
        unsigned current_chunk=std::numeric_limits<unsigned>::max();
        for (unsigned i=0u;i<24u;++i) {
            // Chunk only the large TSV files for reliable transfer.
            // This changes no source cases and introduces no sampling.
            const unsigned chunk=bits==4096u ? i/6u
                : bits==2048u ? i/12u : 0u;
            if (chunk!=current_chunk) {
                if (reference.is_open()) {
                    std::string excess;
                    if (reference >> excess) {
                        throw std::runtime_error("Extra R31 TSV rows");
                    }
                    reference.close();
                    reference.clear();
                }
                const std::string file=std::string(ARIADNE_TEST_DATA_DIR)
                    +"/rational_r31_"+std::to_string(bits)
                    +"_"+std::to_string(chunk)+".tsv";
                reference.open(file);
                if (!reference) {
                    throw std::runtime_error("Missing R31 oracle TSV: "+file);
                }
                std::string header;
                while (reference.peek()=='#') std::getline(reference,header);
                current_chunk=chunk;
            }
            const Integer an=base+Integer(13u+2u*i);
            const Integer ad=base*2+Integer(17u+3u*i);
            const Integer bn=base+Integer(31u+5u*i);
            const Integer bd=base*3+Integer(19u+7u*i);
            for (unsigned signs=0u;signs<4u;++signs) {
                const Rational left((signs&1u)!=0u ? -an : an,ad);
                const Rational right((signs&2u)!=0u ? -bn : bn,bd);
                unsigned ref_bits=0u,ref_i=0u,ref_signs=0u;
                std::string left_num,left_den,right_num,right_den;
                std::string sum_num,sum_den,diff_num,diff_den;
                std::string mul_num,mul_den,div_num,div_den;
                if (!(reference >> ref_bits >> ref_i >> ref_signs
                      >> left_num >> left_den >> right_num >> right_den
                      >> sum_num >> sum_den >> diff_num >> diff_den
                      >> mul_num >> mul_den >> div_num >> div_den)
                    || ref_bits!=bits || ref_i!=i || ref_signs!=signs) {
                    throw std::runtime_error("Missing or unordered R31 TSV row");
                }
                // Check canonicality, not merely fraction equivalence.
                const auto check_result=[](Rational const& actual,
                                           std::string const& num,
                                           std::string const& den) {
                    ARIADNE_TEST_EQUALS(actual.numerator(),Integer(String(num)));
                    ARIADNE_TEST_EQUALS(actual.denominator(),Integer(String(den)));
                };
                check_result(left,left_num,left_den);
                check_result(right,right_num,right_den);
                check_result(left+right,sum_num,sum_den);
                check_result(left-right,diff_num,diff_den);
                check_result(left*right,mul_num,mul_den);
                check_result(left/right,div_num,div_den);
                Rational inplace(left);
                inplace+=right;
                check_result(inplace,sum_num,sum_den);
                inplace=left;
                inplace*=right;
                check_result(inplace,mul_num,mul_den);
                inplace=left;
                inplace/=right;
                check_result(inplace,div_num,div_den);
            }
        }
        std::string excess;
        if (reference >> excess) {
            throw std::runtime_error("Extra R31 TSV rows");
        }
    }
    // R17: all original 7*6*6 cross-product comparisons, including
    // reverse/negative order and the 8/9-limb stack threshold.
    std::ifstream r17_golden(R17_ORACLE_PATH);
    if (!r17_golden) throw std::runtime_error("Missing R17 oracle TSV");
    std::string r17_header;
    while (r17_golden.peek()=='#') std::getline(r17_golden,r17_header);
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
        for (unsigned i=0u;i<probes.size();++i) {
            for (unsigned j=0u;j<probes.size();++j) {
                unsigned ref_bits=0u,ref_i=0u,ref_j=0u;
                int raw=2;
                if (!(r17_golden >> ref_bits >> ref_i >> ref_j >> raw) ||
                    ref_bits!=bits || ref_i!=i || ref_j!=j ||
                    raw< -1 || raw>1) {
                    throw std::runtime_error("Missing or unordered R17 TSV row");
                }
                const Comparison expected=raw<0 ? Comparison::LESS
                    : raw>0 ? Comparison::GREATER : Comparison::EQUAL;
                const Comparison reverse=raw>0 ? Comparison::LESS
                    : raw<0 ? Comparison::GREATER : Comparison::EQUAL;
                const Rational& a=probes[i];
                const Rational& b=probes[j];
                ARIADNE_TEST_EQUALS(cmp(a,b),expected);
                ARIADNE_TEST_EQUALS(cmp(b,a),reverse);
                ARIADNE_TEST_EQUALS(cmp(-a,-b),reverse);
            }
        }
    }
    std::string r17_extra;
    if (r17_golden >> r17_extra)
        throw std::runtime_error("Extra R17 oracle TSV rows");

    // R4: exact independently frozen floor/ceil results for all nine
    // fraction classes across all six original bit-width boundaries.
    std::ifstream r4_golden(R4_ORACLE_PATH);
    if (!r4_golden) throw std::runtime_error("Missing R4 floor/ceil TSV");
    std::string r4_header;
    while (r4_golden.peek()=='#') std::getline(r4_golden,r4_header);
    for (unsigned bits : {53u,128u,256u,1024u,2048u,4096u}) {
        const Integer base=pow(Integer(2),Nat(bits));
        const Integer below=base-1;
        const Integer above=base+1;
        const std::array<Rational,9u> cases={
            Rational(0),Rational(1,3),Rational(-1,3),
            Rational(below,base),Rational(-below,base),
            Rational(base,base),Rational(-base,base),
            Rational(above,base),Rational(-above,base)
        };
        for (std::size_t index=0u;index<cases.size();++index) {
            unsigned ref_bits=0u,ref_index=0u;
            std::string numerator,denominator,ref_floor,ref_ceil;
            if (!(r4_golden >> ref_bits >> ref_index >> numerator
                   >> denominator >> ref_floor >> ref_ceil)
                || ref_bits!=bits || ref_index!=index)
                throw std::runtime_error("Missing or unordered R4 TSV row");
            const Rational& value=cases[index];
            ARIADNE_TEST_EQUALS(value.numerator(),Integer(String(numerator)));
            ARIADNE_TEST_EQUALS(value.denominator(),Integer(String(denominator)));
            ARIADNE_TEST_EQUALS(floor(value),Integer(String(ref_floor)));
            ARIADNE_TEST_EQUALS(ceil(value),Integer(String(ref_ceil)));
        }
    }
    std::string extra_r4;
    if (r4_golden >> extra_r4)
        throw std::runtime_error("Extra R4 oracle TSV rows");

    // R3: independently frozen exact ordering at the 8/9-limb boundary.
    std::ifstream r3_golden(R3_ORACLE_PATH);
    if (!r3_golden) throw std::runtime_error("Missing R3 oracle TSV");
    std::string r3_header;
    while (r3_golden.peek()=='#') std::getline(r3_golden,r3_header);
    for (unsigned bits : {53u,128u,256u,448u,512u,1024u}) {
        const Integer base=pow(Integer(2),Nat(bits));
        const Rational a(base+13,base*3+19);
        const Rational b(base+31,base*2+17);
        unsigned ref_bits=0u;
        int raw=2;
        if (!(r3_golden >> ref_bits >> raw) ||
            ref_bits!=bits || raw< -1 || raw>1) {
            throw std::runtime_error("Missing or unordered R3 TSV row");
        }
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
    std::string r3_extra;
    if (r3_golden >> r3_extra)
        throw std::runtime_error("Extra R3 oracle TSV rows");
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

void TestRational::test_differential() {
    std::mt19937_64 random(0x6672616374696f6eULL);
    bool conversion_matches=true;
    // R0: 120 independently frozen representatives of the original
    // 1200 seeded random visits. One representative per ten visits;
    // all 17 input widths and source-sign classes remain represented.
    std::ifstream random_fixture(R0_RANDOM_ORACLE_PATH);
    if (!random_fixture)
        throw std::runtime_error("Missing Rational random TSV");
    std::string random_header;
    while (random_fixture.peek()=='#')
        std::getline(random_fixture,random_header);
    const auto verify_fraction=[](Rational const& actual,
                                  std::string const& expected_text) {
        const std::size_t slash=expected_text.find('/');
        if (slash==std::string::npos)
            throw std::runtime_error("Invalid Rational exact TSV fraction");
        const Integer numerator{String(expected_text.substr(0u,slash))};
        const Integer denominator{String(expected_text.substr(slash+1u))};
        ARIADNE_TEST_EQUALS(actual.numerator(),numerator);
        ARIADNE_TEST_EQUALS(actual.denominator(),denominator);
    };
    for (unsigned sample=0u;sample<120u;++sample) {
        unsigned original=0u;
        std::size_t width=0u;
        std::string an,ad,bn,bd,source_a,source_b;
        std::string sum,difference,product,squared,halved;
        std::string quotient,reciprocal,floored,ceiled,binary64;
        int order=2;
        if (!(random_fixture >> original >> width >> an >> ad >> bn >> bd
              >> source_a >> source_b >> sum >> difference >> product
              >> squared >> halved >> quotient >> reciprocal
              >> floored >> ceiled >> order >> binary64)
            || original/10u!=sample || width!=1u+original%17u
            || order < -1 || order > 1 || binary64.size()!=16u) {
            throw std::runtime_error("Missing or unordered Rational random TSV row");
        }
        const Rational a{Integer(String(an)),Integer(String(ad))};
        const Rational b{Integer(String(bn)),Integer(String(bd))};
        verify_fraction(a,source_a);
        verify_fraction(b,source_b);
        verify_fraction(a+b,sum);
        verify_fraction(a-b,difference);
        verify_fraction(a*b,product);
        verify_fraction(sqr(a),squared);
        verify_fraction(hlf(a),halved);
        if (quotient=="-") {
            ARIADNE_TEST_ASSERT(is_zero(b));
        } else {
            ARIADNE_TEST_ASSERT(!is_zero(b));
            verify_fraction(a/b,quotient);
        }
        if (reciprocal=="-") {
            ARIADNE_TEST_ASSERT(is_zero(a));
        } else {
            ARIADNE_TEST_ASSERT(!is_zero(a));
            verify_fraction(rec(a),reciprocal);
        }
        ARIADNE_TEST_EQUALS(floor(a),Integer(String(floored)));
        ARIADNE_TEST_EQUALS(ceil(a),Integer(String(ceiled)));
        const Comparison want=order<0 ? Comparison::LESS
            : order>0 ? Comparison::GREATER : Comparison::EQUAL;
        ARIADNE_TEST_EQUALS(cmp(a,b),want);
        ARIADNE_TEST_EQUALS(bool(a==b),order==0);
        const std::uint64_t expected_bits=static_cast<std::uint64_t>(
            std::stoull(binary64,nullptr,16));
        ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>(a.get_d()),
                            expected_bits);
    }
    std::string random_extra;
    if (random_fixture >> random_extra)
        throw std::runtime_error("Extra Rational random TSV rows");

    // R29/R13 reuse this engine. Preserve precisely the 1200 original
    // mt19937_64 streams, but perform no discarded GMP calculations.
    // Each original iteration read four word vectors and four sign bits.
    for (unsigned iteration=0u;iteration<1200u;++iteration) {
        const std::size_t count=1u+iteration%17u;
        for (unsigned operand=0u;operand<4u;++operand) {
            for (std::size_t word=0u;word<count;++word) (void)random();
            (void)random();
        }
    }
    // R16: all original operand widths, GCD classes and both orders.
    // Independent exact expectations are frozen with canonical operands.
    std::ifstream r16_file(R16_ORACLE_PATH);
    if (!r16_file) throw std::runtime_error("Missing Rational R16 TSV");
    std::string r16_header;
    while (r16_file.peek()=='#') std::getline(r16_file,r16_header);
    const auto parse_r16=[](std::string const& text) {
        const auto pos=text.find('/');
        if (pos==std::string::npos)
            throw std::runtime_error("Invalid R16 Rational fixture");
        return Rational{Integer(String(text.substr(0u,pos))),
                        Integer(String(text.substr(pos+1u)))};
    };
    for (unsigned bits:{53u,128u,256u,1024u,2048u,4096u}) {
        for (unsigned pair=0u;pair<7u;++pair) {
            for (unsigned order=0u;order<2u;++order) {
                unsigned got_bits=0u,got_pair=0u,got_order=0u;
                std::string raw_a,raw_b,sum,difference,doubled,zero;
                if (!(r16_file >> got_bits >> got_pair >> got_order
                      >> raw_a >> raw_b >> sum >> difference >> doubled >> zero)
                    || got_bits!=bits || got_pair!=pair || got_order!=order)
                    throw std::runtime_error("Missing or unordered R16 TSV row");
                const Rational a=parse_r16(raw_a),b=parse_r16(raw_b);
                const Rational expected_sum=parse_r16(sum);
                const Rational expected_difference=parse_r16(difference);
                const Rational expected_twice=parse_r16(doubled);
                const Rational expected_zero=parse_r16(zero);
                ARIADNE_TEST_EQUALS(a+b,expected_sum);
                ARIADNE_TEST_EQUALS(a-b,expected_difference);
                Rational direct=a;
                direct+=b;
                ARIADNE_TEST_EQUALS(direct,expected_sum);
                direct=b;
                direct=a+direct;
                ARIADNE_TEST_EQUALS(direct,expected_sum);
                direct=Rational::nan();
                direct=a+b;
                ARIADNE_TEST_EQUALS(direct,expected_sum);
                direct=a;
                direct-=b;
                ARIADNE_TEST_EQUALS(direct,expected_difference);
                direct=b;
                direct=a-direct;
                ARIADNE_TEST_EQUALS(direct,expected_difference);
                direct=Rational::inf(Sign::NEGATIVE);
                direct=a-b;
                ARIADNE_TEST_EQUALS(direct,expected_difference);
                Rational aliased=a;
                aliased+=a;
                ARIADNE_TEST_EQUALS(aliased,expected_twice);
                direct=a;
                Rational& direct_alias=direct;
                direct+=direct_alias;
                ARIADNE_TEST_EQUALS(direct,expected_twice);
                aliased=a;
                aliased-=a;
                ARIADNE_TEST_EQUALS(aliased,expected_zero);
                direct-=direct_alias;
                ARIADNE_TEST_EQUALS(direct,expected_zero);
            }
        }
    }
    std::string extra_r16;
    if (r16_file >> extra_r16)
        throw std::runtime_error("Extra Rational R16 TSV rows");
    // R20: all ten limb-storage transition widths, four cancellation
    // patterns and all four sign combinations remain in a frozen oracle.
    std::ifstream r20_fixture(R20_ORACLE_PATH);
    if (!r20_fixture) throw std::runtime_error("Missing Rational R20 TSV");
    std::string r20_header;
    while (r20_fixture.peek()=='#') std::getline(r20_fixture,r20_header);
    for (unsigned bits:{53u,128u,255u,256u,320u,447u,
                         448u,511u,512u,1024u}) {
        for (unsigned pair=0u;pair<4u;++pair) {
            for (unsigned signs=0u;signs<4u;++signs) {
                unsigned ref_bits=0u,ref_pair=0u,ref_signs=0u;
                std::string lhs_text,rhs_text,sum_text,sub_text;
                if (!(r20_fixture >> ref_bits >> ref_pair >> ref_signs
                      >> lhs_text >> rhs_text >> sum_text >> sub_text)
                    || ref_bits!=bits || ref_pair!=pair || ref_signs!=signs) {
                    throw std::runtime_error("Missing or unordered R20 TSV row");
                }
                const Rational lhs=parse_r16(lhs_text),rhs=parse_r16(rhs_text);
                verify_fraction(lhs,lhs_text);
                verify_fraction(rhs,rhs_text);
                verify_fraction(lhs+rhs,sum_text);
                Rational running=lhs;
                running+=rhs;
                verify_fraction(running,sum_text);
                verify_fraction(lhs-rhs,sub_text);
                running=lhs;
                running-=rhs;
                verify_fraction(running,sub_text);
            }
        }
    }
    std::string r20_extra;
    if (r20_fixture >> r20_extra)
        throw std::runtime_error("Extra R20 TSV rows");
    // R19: exact source-keyed product and division oracle across all
    // four cross-GCD classes, both sign masks and both operand orders.
    std::ifstream r19_file(R19_ORACLE_PATH);
    if (!r19_file) throw std::runtime_error("Missing Rational R19 TSV");
    std::string r19_header;
    while (r19_file.peek()=='#') std::getline(r19_file,r19_header);
    for (unsigned bits:{53u,128u,256u,1024u,2048u,4096u}) {
        for (unsigned factor=0u;factor<5u;++factor) {
            for (unsigned signs=0u;signs<(factor==4u?1u:4u);++signs) {
                unsigned ref_bits=0u,ref_factor=0u,ref_signs=0u;
                std::string lhs_text,rhs_text,mul_text,div_text,reverse_text;
                if (!(r19_file >> ref_bits >> ref_factor >> ref_signs
                      >> lhs_text >> rhs_text >> mul_text >> div_text >> reverse_text)
                    || ref_bits!=bits || ref_factor!=factor || ref_signs!=signs)
                    throw std::runtime_error("Missing or unordered R19 TSV row");
                const Rational lhs=parse_r16(lhs_text),rhs=parse_r16(rhs_text);
                verify_fraction(lhs,lhs_text);
                verify_fraction(rhs,rhs_text);
                verify_fraction(lhs*rhs,mul_text);
                verify_fraction(rhs*lhs,mul_text);
                Rational assigned=lhs;
                assigned*=rhs;
                verify_fraction(assigned,mul_text);
                verify_fraction(lhs/rhs,div_text);
                assigned=lhs;
                assigned/=rhs;
                verify_fraction(assigned,div_text);
                verify_fraction(rhs/lhs,reverse_text);
                if (factor==4u) {
                    Rational self=lhs;
                    self*=self;
                    verify_fraction(self,mul_text);
                    self=lhs;
                    Rational& alias=self;
                    self/=alias;
                    verify_fraction(self,div_text);
                }
            }
        }
    }
    std::string r19_extra;
    if (r19_file >> r19_extra)
        throw std::runtime_error("Extra R19 TSV rows");
    // R22: source-keyed exact cancellation across all original
    // factor/sign combinations and ten eight-limb buffer boundaries.
    // Keep every structurally distinct source case, not just one/ten.
    std::ifstream r22_fixture(R22_ORACLE_PATH);
    if (!r22_fixture) throw std::runtime_error("Missing Rational R22 TSV");
    std::string r22_header;
    while (r22_fixture.peek()=='#') std::getline(r22_fixture,r22_header);
    for (unsigned bits:{53u,127u,128u,255u,256u,447u,448u,
                        511u,512u,1024u}) {
        const Integer base=pow(Integer(2),Nat(bits));
        const Rational a(15*(base+1),14*base);
        const std::array<Rational,4u> factors={
            Rational(1),Rational(1,15),Rational(14),Rational(14,15)};
        ARIADNE_TEST_ASSERT(a*factors[3u]==Rational(base+1,base));
        ARIADNE_TEST_ASSERT(a/Rational(15,14)==Rational(base+1,base));
        for (unsigned factor=0u;factor<4u;++factor) {
            for (unsigned mask=0u;mask<4u;++mask) {
                unsigned ref_bits=0u,ref_factor=0u,ref_mask=0u;
                std::string lhs_text,rhs_text,product_text,quotient_text;
                if (!(r22_fixture >> ref_bits >> ref_factor >> ref_mask
                      >> lhs_text >> rhs_text >> product_text >> quotient_text)
                    || ref_bits!=bits || ref_factor!=factor || ref_mask!=mask)
                    throw std::runtime_error("Missing or unordered R22 TSV row");
                const Rational lhs=(mask&1u)!=0u ? -a : a;
                const Rational rhs=(mask&2u)!=0u ? -factors[factor] : factors[factor];
                verify_fraction(lhs,lhs_text);
                verify_fraction(rhs,rhs_text);
                verify_fraction(lhs*rhs,product_text);
                Rational assigned=lhs;
                assigned*=rhs;
                verify_fraction(assigned,product_text);
                verify_fraction(lhs/rhs,quotient_text);
                assigned=lhs;
                assigned/=rhs;
                verify_fraction(assigned,quotient_text);
            }
        }
    }
    std::string r22_extra;
    if (r22_fixture >> r22_extra)
        throw std::runtime_error("Extra R22 TSV rows");
    // Independently check truncation at normal/subnormal and overflow edges,
    // under all four host rounding modes. Operands are not formed as doubles.
    const int saved_rounding=std::fegetround();
    const std::array<int,4> modes={FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO};
    // R5: exact high-prefix conversion, 4-to-5-limb boundary, and
    // near-half/near-one quotients with nonzero low tails. Compare raw
    // binary64 patterns with the independent GMP implementation under
    // every rounding mode, including cases requiring the full fallback.
    const std::array<unsigned int,12> shifts={0u,1u,52u,53u,1022u,1023u,1024u,1073u,1074u,1075u,1076u,2048u};
    std::ifstream r7_golden(R7_ORACLE_PATH);
    if (!r7_golden) throw std::runtime_error("Missing R7 oracle TSV");
    std::string r7_header;
    while (r7_golden.peek()=='#') std::getline(r7_golden,r7_header);
    std::ifstream r12_golden(R12_ORACLE_PATH);
    if (!r12_golden) throw std::runtime_error("Missing R12 oracle TSV");
    std::string r12_header;
    while (r12_golden.peek()=='#') std::getline(r12_golden,r12_header);
    std::ifstream r29_golden(R29_ORACLE_PATH);
    if (!r29_golden) throw std::runtime_error("Missing R29 oracle TSV");
    std::string r29_header;
    while (r29_golden.peek()=='#') std::getline(r29_golden,r29_header);
    std::ifstream r5_golden(R5_ORACLE_PATH);
    if (!r5_golden) throw std::runtime_error("Missing R5 oracle TSV");
    std::string r5_header;
    while (r5_golden.peek()=='#') std::getline(r5_golden,r5_header);
    for (unsigned mode_index=0u;mode_index<modes.size();++mode_index) {
        const int mode=modes[mode_index];
        ARIADNE_TEST_EQUALS(std::fesetround(mode),0);
        // R12: all original one-limb sources and both signs, with
        // independently frozen bit-exact binary64 results. Validate each
        // source against the original mt19937_64 sequence as well.
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
            unsigned ref_mode=0u,ref_sample=0u;
            std::uint64_t ref_num=0u,ref_den=0u;
            std::uint64_t expected_positive=0u,expected_negative=0u;
            if (!(r12_golden >> ref_mode >> ref_sample >> ref_num >> ref_den
                  >> expected_positive >> expected_negative)
                || ref_mode!=mode_index || ref_sample!=sample
                || ref_num!=numerator || ref_den!=denominator) {
                throw std::runtime_error("Missing, misordered or source-mismatched R12 TSV row");
            }
            const Rational positive{Integer(numerator),Integer(denominator)};
            const Rational negative=-positive;
            ARIADNE_TEST_EQUALS(
                std::bit_cast<std::uint64_t>(positive.get_d()),expected_positive);
            ARIADNE_TEST_EQUALS(
                std::bit_cast<std::uint64_t>(negative.get_d()),expected_negative);
        }
        // R29: preserve every original sparse-prefix source and both signs.
        // The fixture checks the exact PRNG tails and canonical binary64 bits.
        for (unsigned bits:{128u,256u,320u,512u,1024u,2048u}) {
            const Integer base=pow(Integer(2),Nat(bits));
            for (std::uint64_t high_n:{1u,3u,7u}) {
                for (std::uint64_t high_d:{2u,3u,5u}) {
                    const Integer n=base*Integer(high_n)
                        +Integer((random()&65535u)+1u);
                    const Integer d=base*Integer(high_d)
                        +Integer((random()&65535u)+1u);
                    unsigned ref_mode=0u,ref_bits=0u;
                    std::uint64_t ref_high_n=0u,ref_high_d=0u;
                    std::string ref_n,ref_d,positive_token,negative_token;
                    // Parse before validating: Integer comparison returns an
                    // Ariadne Boolean, whose overloaded || cannot short-circuit.
                    if (!(r29_golden >> ref_mode >> ref_bits >> ref_high_n
                          >> ref_high_d >> ref_n >> ref_d
                          >> positive_token >> negative_token)) {
                        throw std::runtime_error("Malformed or missing R29 TSV row");
                    }
                    if (ref_mode!=mode_index || ref_bits!=bits
                        || ref_high_n!=high_n || ref_high_d!=high_d
                        || ref_n.empty() || ref_d.empty()
                        || std::string(n.literal().c_str())!=ref_n
                        || std::string(d.literal().c_str())!=ref_d) {
                        throw std::runtime_error("Out-of-order or source-mismatched R29 TSV row");
                    }
                    const std::uint64_t positive_bits=std::stoull(positive_token);
                    const std::uint64_t negative_bits=std::stoull(negative_token);
                    const Rational positive(n,d);
                    const Rational negative=-positive;
                    ARIADNE_TEST_EQUALS(
                        std::bit_cast<std::uint64_t>(positive.get_d()),positive_bits);
                    ARIADNE_TEST_EQUALS(
                        std::bit_cast<std::uint64_t>(negative.get_d()),negative_bits);
                }
            }
        }
        // R13: all original 4 widths * 96 seeded dense cases, per host
        // mode. The source-keyed fixture retains both full signed operands.
        for (unsigned limbs=5u;limbs<=8u;++limbs) {
            std::vector<std::uint64_t> numerator_limbs(limbs);
            std::vector<std::uint64_t> denominator_limbs(limbs);
            const std::string fixture_path=std::string(ARIADNE_TEST_DATA_DIR)
                +"/rational_r13_mode"+std::to_string(mode_index)
                +"_limbs"+std::to_string(limbs)+".tsv";
            std::ifstream fixture(fixture_path);
            if (!fixture) throw std::runtime_error("Missing R13 oracle TSV: "+fixture_path);
            std::string header;
            while (fixture.peek()=='#') std::getline(fixture,header);
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
                unsigned ref_mode=0u,ref_limbs=0u,ref_sample=0u;
                std::string ref_n,ref_d,bits_token;
                if (!(fixture >> ref_mode >> ref_limbs >> ref_sample
                      >> ref_n >> ref_d >> bits_token)) {
                    throw std::runtime_error("Malformed or missing R13 TSV row");
                }
                if (ref_mode!=mode_index || ref_limbs!=limbs
                    || ref_sample!=sample || ref_n.empty() || ref_d.empty()
                    || std::string(num.literal().c_str())!=ref_n
                    || std::string(den.literal().c_str())!=ref_d) {
                    throw std::runtime_error("Out-of-order or source-mismatched R13 TSV row");
                }
                const std::uint64_t ref_bits=std::stoull(bits_token);
                const Rational value(num,den);
                ARIADNE_TEST_EQUALS(std::bit_cast<std::uint64_t>(value.get_d()),ref_bits);
            }
            std::string extra;
            if (fixture >> extra) throw std::runtime_error("Extra R13 TSV rows");
        }
        // R7: preserve all twelve rational source cases at eight widths
        // and all four host modes; compare exact independently frozen
        // binary64 bit patterns, not a rounded/toleranced double.
        for (unsigned bits : {53u,128u,256u,448u,512u,1024u,2048u,4096u}) {
            const Integer base=pow(Integer(2),Nat(bits));
            const Integer twice=base*2;
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
            for (unsigned probe=0u;probe<probes.size();++probe) {
                unsigned ref_mode=0u,ref_bits=0u,ref_probe=0u;
                std::uint64_t expected_bits=0u;
                if (!(r7_golden >> ref_mode >> ref_bits >> ref_probe >> expected_bits)
                    || ref_mode!=mode_index || ref_bits!=bits
                    || ref_probe!=probe) {
                    throw std::runtime_error("Missing or unordered R7 TSV row");
                }
                conversion_matches=conversion_matches &&
                    std::bit_cast<std::uint64_t>(probes[probe].get_d())==
                    expected_bits;
            }
        }
        // R5: all original shifts and four signed, scaled operands.
        // Frozen raw binary64 patterns are independent of Numeric and GMP.
        for (unsigned int shift : shifts) {
            const Integer factor=pow(Integer(2),Nat(shift));
            const Integer n=Integer(random()|1u),d=Integer(random()|1u);
            unsigned ref_mode=0u,ref_shift=0u;
            std::string ref_n,ref_d;
            std::array<std::string,4> bit_tokens;
            if (!(r5_golden >> ref_mode >> ref_shift >> ref_n >> ref_d
                  >> bit_tokens[0] >> bit_tokens[1]
                  >> bit_tokens[2] >> bit_tokens[3])) {
                throw std::runtime_error("Malformed or missing R5 TSV row");
            }
            if (ref_mode!=mode_index || ref_shift!=shift
                || ref_n.empty() || ref_d.empty()
                || std::string(n.literal().c_str())!=ref_n
                || std::string(d.literal().c_str())!=ref_d) {
                throw std::runtime_error("Out-of-order or source-mismatched R5 TSV row");
            }
            std::array<std::uint64_t,4> expected_bits{};
            for (unsigned index=0u;index<4u;++index) {
                expected_bits[index]=std::stoull(bit_tokens[index]);
            }
            const std::array<Rational,4> values={Rational(n,d*factor),Rational(-n,d*factor),
                                                Rational(n*factor,d),Rational(-n*factor,d)};
            for (unsigned index=0u;index<values.size();++index) {
                const auto actual=std::bit_cast<std::uint64_t>(values[index].get_d());
                if (actual!=expected_bits[index]) {
                    std::cerr << "Rational double mismatch: mode=" << mode << " shift=" << shift
                              << " case=" << index << " actual=" << actual
                              << " expected=" << expected_bits[index] << '\n';
                    conversion_matches=false;
                }
            }
        }
    }
    ARIADNE_TEST_EQUALS(std::fesetround(saved_rounding),0);
    std::string r7_extra;
    if (r7_golden >> r7_extra)
        throw std::runtime_error("Extra R7 oracle TSV rows");
    std::string r12_extra;
    if (r12_golden >> r12_extra)
        throw std::runtime_error("Extra R12 oracle TSV rows");
    std::string r29_extra;
    if (r29_golden >> r29_extra)
        throw std::runtime_error("Extra R29 oracle TSV rows");
    std::string r5_extra;
    if (r5_golden >> r5_extra)
        throw std::runtime_error("Extra R5 oracle TSV rows");
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
