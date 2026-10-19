/***************************************************************************
 *            test_integer.cpp
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

#include "numeric/integer.hpp"
#include "numeric/detail/integer_gmp_bridge.hpp"
#include "numeric/int.hpp"
#include "numeric/twoexp.hpp"
#include "numeric/accuracy.hpp"
#include "numeric/sign.hpp"
#include "foundation/logical.hpp"
#include "utility/string.hpp"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include <gmp.h>

#include "macro/testing.hpp"

using namespace std;
using namespace Ariadne;

namespace {

constexpr bool check_concept()
{
    return requires(unsigned int m, unsigned long int lm, int n, long int ln, Integer z, Integer z2, Boolean b) {
        z=Integer(); z=Integer(m); z=Integer(lm); z=Integer(n); z=Integer(ln); z=Integer(z);
        z2=Integer();
        z=m; z=lm; z=n; z=ln; z=z2;
        z=+z; z=-z;
        z=z+z; z=z-z; z=z*z;
        z=z+n; z=z-n; z=z*n;
        z=n+z; z=n-z; z=n*z;
        z=1_z; z=-1_z;
        b=(z==z); b=(z!=z); b=(z<=z); b=(z>=z); b=(z<z); b=(z>z);
        b=(z==n); b=(z!=n); b=(z<=n); b=(z>=n); b=(z<n); b=(z>n);
        b=(n==z); b=(n!=z); b=(n<=z); b=(n>=z); b=(n<z); b=(n>z);
    };
}

static_assert(check_concept());

} // namespace



class TestInteger
{
  public:
    void test();
  private:
    void test_literal();
    void test_constructors();
    void test_arithmetic();
    void test_comparisons();
    void test_small_types();
    void test_combinatorics();
    void test_native();
    void test_signed_regression();
};

void TestInteger::test()
{
    ARIADNE_TEST_CALL(test_constructors());
    ARIADNE_TEST_CALL(test_comparisons());
    ARIADNE_TEST_CALL(test_literal());
    ARIADNE_TEST_CALL(test_arithmetic());
    ARIADNE_TEST_CALL(test_small_types());
    ARIADNE_TEST_CALL(test_combinatorics());
    ARIADNE_TEST_CALL(test_native());
    ARIADNE_TEST_CALL(test_signed_regression());
}

void TestInteger::test_literal() {
    ARIADNE_TEST_CONSTRUCT(Integer,z,(3_z));
    ARIADNE_TEST_EQUALS(z,3);
    ARIADNE_TEST_EQUALS(z,Integer(3));
    ARIADNE_TEST_EQUALS(100000_z,Integer(100000));
    ARIADNE_TEST_EQUAL(1000000000000_z,sqr(Integer(1000000)));
    ARIADNE_TEST_EQUALS(1000000000000_z-sqr(Integer(1000000)),0);
    ARIADNE_TEST_EQUALS(4294967295_z,4294967295u);
    ARIADNE_TEST_EQUALS(-2147483647_z,-2147483647);
    ARIADNE_TEST_EQUALS(4294967295_z,Integer(4294967295u));
    ARIADNE_TEST_EQUALS(4611686016279904256_z,Integer(2147483647)*2147483647+2147483647);
}

void TestInteger::test_constructors() {
    int m=2147483647;
    unsigned int um=2147483647u;
    unsigned long int ulm=um;
    unsigned long long int ullm=um;
    int n=-2147483647;
    long int ln=n;
    long long int lln=n;
    String sn="-2147483647";
    String sz="314159265358979323846264";

    ARIADNE_TEST_ASSERT((not Constructible<Integer,float>));
    ARIADNE_TEST_ASSERT((not Constructible<Integer,double>));

    ARIADNE_TEST_CONSTRUCT(Integer,zum,(um));
    ARIADNE_TEST_EQUALS(zum.get_si(),m);
    ARIADNE_TEST_CONSTRUCT(Integer,zulm,(ulm));
    ARIADNE_TEST_EQUALS(zulm.get_si(),m);
    ARIADNE_TEST_CONSTRUCT(Integer,zullm,(ullm));
    ARIADNE_TEST_EQUALS(zullm.get_si(),m);
    ARIADNE_TEST_CONSTRUCT(Integer,zn,(n));
    ARIADNE_TEST_EQUALS(zn.get_si(),n);
    ARIADNE_TEST_CONSTRUCT(Integer,zln,(ln));
    ARIADNE_TEST_EQUALS(zln.get_si(),n);
    ARIADNE_TEST_CONSTRUCT(Integer,zlln,(lln));
    ARIADNE_TEST_EQUALS(zlln.get_si(),n);
    ARIADNE_TEST_CONSTRUCT(Integer,zsn,(sn));
    ARIADNE_TEST_EQUALS(zsn.get_si(),n);
    ARIADNE_TEST_CONSTRUCT(Integer,zsz,(sz));
    ARIADNE_TEST_EQUALS(zsz,Integer(314159265358)*Integer(1000000000000)+Integer(979323846264));

    ARIADNE_TEST_CONSTRUCT(Integer,z1,(0));
    ARIADNE_TEST_EQUALS(z1.get_si(),0);
    ARIADNE_TEST_CONSTRUCT(Integer,z2,(-3));
    ARIADNE_TEST_EQUALS(z2.get_si(),-3);
    ARIADNE_TEST_CONSTRUCT(Integer,z3,(ullm*ullm+ullm));
    ARIADNE_TEST_EQUALS(z3,zum*zum+zum);
    ARIADNE_TEST_CONSTRUCT(Integer,z4,(lln*lln+lln));
    ARIADNE_TEST_EQUALS(z4,zn*zn+zn);

    ARIADNE_TEST_EQUALS(Integer(Nat32(7)),Integer(7));
    mpz_t raw;
    mpz_init_set_si(raw,-17);
    ARIADNE_TEST_EQUALS(integer_from_gmp(raw),Integer(-17));
    mpz_clear(raw);

    Integer assigned(1);
    Integer source(9);
    assigned=source;
    ARIADNE_TEST_EQUALS(assigned,Integer(9));
}

void TestInteger::test_arithmetic() {
    ARIADNE_TEST_EQUALS(+Integer(-5),-5);
    ARIADNE_TEST_EQUALS(-Integer(-5), 5);
    ARIADNE_TEST_EQUALS(Integer(-5)+Integer(2),-3);
    ARIADNE_TEST_EQUALS(Integer(-5)-Integer(2),-7);
    ARIADNE_TEST_EQUALS(Integer(-5)*Integer(2),-10);

    ARIADNE_TEST_EQUALS(pos(Integer(-5)),-5);
    ARIADNE_TEST_EQUALS(neg(Integer(-5)), 5);
    ARIADNE_TEST_EQUALS(sqr(Integer(-5)),25);
    ARIADNE_TEST_EQUALS(pow(Integer(-5),3u),-125);


    ARIADNE_TEST_EQUALS((Integer)max(Integer(5),Integer(3)),5);
    ARIADNE_TEST_EQUALS((Integer)max(Integer(-5),Integer(-3)),-3);
    ARIADNE_TEST_EQUALS((Integer)min(Integer(5),Integer(3)),3);
    ARIADNE_TEST_EQUALS((Integer)min(Integer(-5),Integer(-3)),-5);
    ARIADNE_TEST_EQUALS(abs(Integer(-5)),5);
    ARIADNE_TEST_EQUALS(abs(Integer( 0)),0);
    ARIADNE_TEST_EQUALS(abs(Integer(+5)),5);

    Integer mutable_value(5);
    --mutable_value;
    ARIADNE_TEST_EQUALS(mutable_value,Integer(4));
    mutable_value-=Integer(3);
    ARIADNE_TEST_EQUALS(mutable_value,Integer(1));

    ARIADNE_TEST_ASSERT(not is_nan(Integer(1)));
    ARIADNE_TEST_ASSERT(not is_inf(Integer(1)));
    ARIADNE_TEST_ASSERT(is_finite(Integer(1)));
    ARIADNE_TEST_ASSERT(is_zero(Integer(0)));
    ARIADNE_TEST_ASSERT(not is_zero(Integer(1)));
    ARIADNE_TEST_ASSERT(not is_zero(Integer(-1)));
    ARIADNE_TEST_EQUALS(sgn(Integer(-1)),Sign::NEGATIVE);
    ARIADNE_TEST_EQUALS(sgn(Integer(0)),Sign::ZERO);
    ARIADNE_TEST_EQUALS(sgn(Integer(1)),Sign::POSITIVE);

    ARIADNE_TEST_EQUALS(cmp(Integer(2),Integer(3)),Comparison::LESS);
    ARIADNE_TEST_EQUALS(cmp(Integer(3),Integer(3)),Comparison::EQUAL);
    ARIADNE_TEST_EQUALS(cmp(Integer(4),Integer(3)),Comparison::GREATER);
}

void TestInteger::test_small_types() {
    Nat32 n32;
    Nat64 n64;
    Int32 i32;
    Int64 i64;
    ARIADNE_TEST_EQUALS(n32.get_ui(),uint32_t(0));
    ARIADNE_TEST_EQUALS(n64.get_ui(),uint64_t(0));
    ARIADNE_TEST_EQUALS(i32.get_si(),int32_t(0));
    ARIADNE_TEST_EQUALS(i64.get_si(),int64_t(0));

    ARIADNE_TEST_EQUALS(Nat32(uint16_t(7)).get_ui(),uint32_t(7));
    ARIADNE_TEST_EQUALS(Nat32(int16_t(7)).get_ui(),uint32_t(7));
    ARIADNE_TEST_EQUALS(Nat64(uint32_t(8)).get_ui(),uint64_t(8));
    ARIADNE_TEST_EQUALS(Nat64(int32_t(8)).get_ui(),uint64_t(8));
    ARIADNE_TEST_EQUALS(Int32(uint16_t(9)).get_si(),int32_t(9));
    ARIADNE_TEST_EQUALS(Int32(int16_t(-9)).get_si(),int32_t(-9));
    ARIADNE_TEST_EQUALS(Int64(uint32_t(10)).get_si(),int64_t(10));
    ARIADNE_TEST_EQUALS(Int64(uint64_t(11)).get_si(),int64_t(11));
    ARIADNE_TEST_EQUALS(Int64(int16_t(-10)).get_si(),int64_t(-10));

    ARIADNE_TEST_FAIL(Nat32(uint64_t(4294967296ULL)));
    ARIADNE_TEST_FAIL(Nat32(int64_t(-1)));
    ARIADNE_TEST_FAIL(Nat32(int64_t(4294967296LL)));
    ARIADNE_TEST_FAIL(Nat64(int64_t(-1)));
    ARIADNE_TEST_FAIL(Int32(uint64_t(2147483648ULL)));
    ARIADNE_TEST_FAIL(Int32(int64_t(-2147483649LL)));
    ARIADNE_TEST_FAIL(Int32(int64_t(2147483648LL)));
    ARIADNE_TEST_FAIL(Int64(uint64_t(9223372036854775808ULL)));

    Bits bits=8_bits;
    ARIADNE_TEST_EQUALS(static_cast<unsigned long int>(bits),8ul);
    std::ostringstream bits_stream;
    bits_stream << bits;
    ARIADNE_TEST_EQUALS(bits_stream.str(),std::string("8_bits"));

    Accuracy direct_accuracy(Dyadic(1,8u));
    ARIADNE_TEST_EQUALS(direct_accuracy.error(),Dyadic(1,8u));
    ARIADNE_TEST_EQUALS(accuracy(8_bits).error(),Dyadic(1,8u));

    TwoExp p(3);
    TwoExp from_two(two);
    ARIADNE_TEST_EQUALS(p.exponent(),3);
    ARIADNE_TEST_EQUALS(from_two.exponent(),1);
    ARIADNE_TEST_EQUALS(p.get_d(),8.0);
    ARIADNE_TEST_EQUALS(rec(p).exponent(),-3);
    ARIADNE_TEST_EQUALS(Ariadne::exp2(-2).exponent(),-2);
    ARIADNE_TEST_EQUALS((two^static_cast<Nat>(4)).exponent(),4);
    ARIADNE_TEST_EQUALS((two^static_cast<Int>(-4)).exponent(),-4);
    ARIADNE_TEST_EQUALS(pow(two,-5).exponent(),-5);

    ARIADNE_TEST_EQUALS(-Sign::NEGATIVE,Sign::POSITIVE);
    ARIADNE_TEST_EQUALS(Sign::NEGATIVE*Sign::NEGATIVE,Sign::POSITIVE);
    ARIADNE_TEST_EQUALS(Sign::NEGATIVE*Sign::POSITIVE,Sign::NEGATIVE);

    ARIADNE_TEST_EQUALS(class_name<uint>(),String("uint"));
    ARIADNE_TEST_EQUALS(class_name<int>(),String("int"));
    ARIADNE_TEST_EQUALS(class_name<Integer>(),String("Integer"));
    ARIADNE_TEST_EQUALS(class_name<Natural>(),String("Natural"));

    std::ostringstream sign_stream;
    sign_stream << Sign::ZERO << " " << Sign::NEGATIVE << " " << Sign::POSITIVE;
    ARIADNE_TEST_EQUALS(sign_stream.str(),std::string("ZERO NEGATIVE POSITIVE"));

    std::ostringstream comparison_stream;
    comparison_stream << LESS << " " << EQUAL << " " << GREATER << " " << INCOMPARABLE;
    ARIADNE_TEST_EQUALS(comparison_stream.str(),std::string("LESS EQUAL GREATER INCOMPARABLE"));
    std::ostringstream invalid_comparison_stream;
    ARIADNE_TEST_FAIL(invalid_comparison_stream << static_cast<Comparison>(42));

    Natural natural(2u);
    Natural other(3u);
    ++natural;
    natural+=other;
    ARIADNE_TEST_EQUALS(natural,Integer(6));
    ARIADNE_TEST_EQUALS(Natural(2u)+Natural(3u),Integer(5));
    ARIADNE_TEST_EQUALS(Natural(2u)*Natural(3u),Integer(6));
    ARIADNE_TEST_EQUALS(log2floor(Natural(0u)),-1);
    ARIADNE_TEST_EQUALS(log2floor(Natural(1u)),0);
    ARIADNE_TEST_EQUALS(log2floor(Natural(8u)),3);
    ARIADNE_TEST_EQUALS(cast_positive(Integer(4)),Integer(4));
    ARIADNE_TEST_FAIL(Positive<Integer>(Integer(-1)));
    ARIADNE_TEST_FAIL(Natural(Integer(-1)));
}

void TestInteger::test_combinatorics() {
    ARIADNE_TEST_EQUALS(fac(uint8_t(0)),uint32_t(1));
    ARIADNE_TEST_EQUALS(fac(uint8_t(5)),uint32_t(120));
    ARIADNE_TEST_FAIL(fac(uint8_t(13)));
    ARIADNE_TEST_EQUALS(bin(uint8_t(5),uint8_t(0)),uint32_t(1));
    ARIADNE_TEST_EQUALS(bin(uint8_t(5),uint8_t(2)),uint32_t(10));
    ARIADNE_TEST_EQUALS(bin(uint8_t(5),uint8_t(6)),uint32_t(0));
    ARIADNE_TEST_FAIL(bin(uint8_t(5),uint8_t(7)));
    ARIADNE_TEST_FAIL(bin(uint8_t(32),uint8_t(1)));

    ARIADNE_TEST_EQUALS(fac(uint16_t(0)),uint16_t(1));
    ARIADNE_TEST_EQUALS(fac(uint16_t(5)),uint16_t(120));
    ARIADNE_TEST_FAIL(fac(uint16_t(9)));
    ARIADNE_TEST_EQUALS(bin(uint16_t(5),uint16_t(0)),uint16_t(1));
    ARIADNE_TEST_EQUALS(bin(uint16_t(5),uint16_t(2)),uint16_t(10));
    ARIADNE_TEST_EQUALS(bin(uint16_t(5),uint16_t(6)),uint16_t(0));
    ARIADNE_TEST_FAIL(bin(uint16_t(5),uint16_t(7)));
    ARIADNE_TEST_FAIL(bin(uint16_t(16),uint16_t(1)));

    ARIADNE_TEST_EQUALS(fac(uint32_t(0)),uint32_t(1));
    ARIADNE_TEST_EQUALS(fac(uint32_t(5)),uint32_t(120));
    ARIADNE_TEST_FAIL(fac(uint32_t(13)));
    ARIADNE_TEST_EQUALS(bin(uint32_t(5),uint32_t(0)),uint32_t(1));
    ARIADNE_TEST_EQUALS(bin(uint32_t(5),uint32_t(2)),uint32_t(10));
    ARIADNE_TEST_EQUALS(bin(uint32_t(5),uint32_t(6)),uint32_t(0));
    ARIADNE_TEST_FAIL(bin(uint32_t(5),uint32_t(7)));
    ARIADNE_TEST_FAIL(bin(uint32_t(31),uint32_t(1)));

    ARIADNE_TEST_EQUALS(fac(uint64_t(0)),uint64_t(1));
    ARIADNE_TEST_EQUALS(fac(uint64_t(5)),uint64_t(120));
    ARIADNE_TEST_FAIL(fac(uint64_t(21)));
    ARIADNE_TEST_EQUALS(bin(uint64_t(5),uint64_t(0)),uint64_t(1));
    ARIADNE_TEST_EQUALS(bin(uint64_t(5),uint64_t(2)),uint64_t(10));
    ARIADNE_TEST_EQUALS(bin(uint64_t(5),uint64_t(6)),uint64_t(0));
    ARIADNE_TEST_FAIL(bin(uint64_t(5),uint64_t(7)));
    ARIADNE_TEST_FAIL(bin(uint64_t(63),uint64_t(1)));
}

void TestInteger::test_comparisons() {
    ARIADNE_TEST_COMPARE(Integer(3),==,3);
    ARIADNE_TEST_COMPARE(3,==,Integer(3));

    ARIADNE_TEST_COMPARE(Integer(2),==,Integer(2));
    ARIADNE_TEST_COMPARE(Integer(0),==,Integer(-0));
    ARIADNE_TEST_COMPARE(Integer(2),!=,Integer(-2));
    ARIADNE_TEST_COMPARE(Integer(2),!=,Integer(-3));
    ARIADNE_TEST_COMPARE(Integer(2),<=,Integer(23));
    ARIADNE_TEST_COMPARE(Integer(2),<=,Integer(3));
    ARIADNE_TEST_COMPARE(Integer(2),>=,Integer(2));
    ARIADNE_TEST_COMPARE(Integer(2),>=,Integer(-3));
    ARIADNE_TEST_COMPARE(Integer(2),< ,Integer(3));
    ARIADNE_TEST_COMPARE(Integer(2),> ,Integer(-3));
}


void TestInteger::test_native() {
    const Integer zero(0);
    const Integer one(1);
    const Integer minus_one(-1);
    ARIADNE_TEST_ASSERT(Integer(String("-0")).literal()=="0");
    ARIADNE_TEST_ASSERT(Integer(String("+00000042")).literal()=="42");
    ARIADNE_TEST_ASSERT(Integer(String("00000")).literal()=="0");
    ARIADNE_TEST_ASSERT(Integer(String("18446744073709551616")).literal()
                        =="18446744073709551616");
    ARIADNE_TEST_ASSERT(Integer(std::numeric_limits<std::int64_t>::min())
        .get<std::int64_t>()==std::numeric_limits<std::int64_t>::min());
    ARIADNE_TEST_ASSERT(Integer(std::numeric_limits<std::uint64_t>::max())
        .get<std::uint64_t>()==std::numeric_limits<std::uint64_t>::max());
    ARIADNE_TEST_ASSERT(pow(Integer(0),Nat(0u))==1);
    ARIADNE_TEST_ASSERT(pow(Integer(-7),Nat(3u))==-343);
    ARIADNE_TEST_ASSERT(pow(Integer(-7),Nat(4u))==2401);
    ARIADNE_TEST_ASSERT(sqr(Integer(-13))==169);
    ARIADNE_TEST_ASSERT(abs(Integer(-13))==13);
    ARIADNE_TEST_ASSERT(log2floor(Natural(0u))==-1);
    ARIADNE_TEST_ASSERT(log2floor(Natural(8u))==3);
    ARIADNE_TEST_ASSERT(quot(Integer(-11),Integer(3))==-3);
    ARIADNE_TEST_ASSERT(rem(Integer(-11),Integer(3))==-2);
    ARIADNE_TEST_ASSERT(quot(Integer(11),Integer(-3))==-3);
    ARIADNE_TEST_ASSERT(rem(Integer(11),Integer(-3))==2);

    for (const char* invalid: {"","-","+","-+1"," 1","1 ","1a","1.0"}) {
        bool threw=false;
        try { (void)Integer(String(invalid)); }
        catch (std::invalid_argument const&) { threw=true; }
        ARIADNE_TEST_ASSERT(threw);
    }
    bool negative_unsigned=false,positive_overflow=false,division_zero=false;
    try { (void)minus_one.get<std::uint64_t>(); }
    catch (std::out_of_range const&) { negative_unsigned=true; }
    try { (void)Integer(String("9223372036854775808")).get<std::int64_t>(); }
    catch (std::out_of_range const&) { positive_overflow=true; }
    try { (void)quot(one,zero); }
    catch (std::domain_error const&) { division_zero=true; }
    ARIADNE_TEST_ASSERT(negative_unsigned && positive_overflow && division_zero);

    for (const char* a: {"-11","11","-12","12","0"}) {
        for (const char* b: {"-3","3","-5","5"}) {
            Integer dividend{String(a)};
            Integer divisor{String(b)};
            const Integer quotient=quot(dividend,divisor);
            const Integer remainder=rem(dividend,divisor);
            ARIADNE_TEST_ASSERT(dividend==quotient*divisor+remainder);
            ARIADNE_TEST_ASSERT(abs(remainder)<abs(divisor));
            ARIADNE_TEST_ASSERT(remainder==0 || sgn(remainder)==sgn(dividend));
        }
    }

    struct GMP {
        mpz_t value;
        GMP() { mpz_init(value); }
        ~GMP() { mpz_clear(value); }
        GMP(GMP const&)=delete;
        GMP& operator=(GMP const&)=delete;
    };
    GMP a,b,expected,quotient,remainder;

    // B41: compare direct signed subtraction with independent GMP
    // across all sign pairs, equal/unequal widths and long common
    // high prefixes. This also checks x-=y and the x-=x alias case.
    std::mt19937_64 sub_engine(0xB41D1FFULL);
    for (const std::size_t count:
            {1u,2u,7u,8u,16u,32u,64u,128u}) {
        for (unsigned repetition=0u;repetition<24u;++repetition) {
            std::vector<std::uint64_t> left(count),right(count);
            for (auto& word:left) word=sub_engine();
            for (auto& word:right) word=sub_engine();
            switch (repetition%6u) {
            case 0u: // Exact equal values.
                right=left;
                break;
            case 1u: // Lower limb differs, high prefix cancels.
                right=left;
                left[0u]=std::numeric_limits<std::uint64_t>::max();
                right[0u]=0u;
                break;
            case 2u: // Reversed ordering of low word.
                right=left;
                left[0u]=0u;
                right[0u]=std::numeric_limits<std::uint64_t>::max();
                break;
            case 3u: // Same sign and opposite sign over zero.
                std::fill(left.begin(),left.end(),0u);
                break;
            case 4u: // Different canonical widths.
                right.back()=0u;
                break;
            default: // Independent dense values.
                break;
            }
            mpz_import(a.value,left.size(),-1,sizeof(std::uint64_t),
                       0,0,left.data());
            mpz_import(b.value,right.size(),-1,sizeof(std::uint64_t),
                       0,0,right.data());
            for (unsigned sign_mode=0u;sign_mode<4u;++sign_mode) {
                if ((sign_mode&1u)!=0u) mpz_neg(a.value,a.value);
                if ((sign_mode&2u)!=0u) mpz_neg(b.value,b.value);
                const Integer lhs=integer_from_gmp(a.value);
                const Integer rhs=integer_from_gmp(b.value);
                mpz_sub(expected.value,a.value,b.value);
                const GmpIntegerBoundary native_sub(lhs-rhs);
                ARIADNE_TEST_ASSERT(
                    mpz_cmp(native_sub.value(),expected.value)==0);
                Integer assigned=lhs;
                assigned-=rhs;
                const GmpIntegerBoundary assigned_sub(assigned);
                ARIADNE_TEST_ASSERT(
                    mpz_cmp(assigned_sub.value(),expected.value)==0);
                Integer aliased=lhs;
                aliased-=aliased;
                ARIADNE_TEST_ASSERT(aliased==0);
                ARIADNE_TEST_ASSERT(!aliased.negative());
                if ((sign_mode&1u)!=0u) mpz_neg(a.value,a.value);
                if ((sign_mode&2u)!=0u) mpz_neg(b.value,b.value);
            }
        }
    }

    // B32: GMP is the independent oracle for the first-set-bit
    // multiplication elision, including 0^0, signs and exponents
    // with leading zero bits or several nonzero bits.
    for (const char* literal:
            {"0","1","-1","2","-2","17","-17",
             "18446744073709551616","-18446744073709551616",
             "340282366920938463463374607431768211455"}) {
        const Integer input{String(literal)};
        ARIADNE_TEST_ASSERT(mpz_set_str(a.value,literal,10)==0);
        for (unsigned exponent=0u;exponent<=12u;++exponent) {
            mpz_pow_ui(expected.value,a.value,exponent);
            const GmpIntegerBoundary native_pow(
                pow(input,Nat(exponent)));
            ARIADNE_TEST_ASSERT(
                mpz_cmp(native_pow.value(),expected.value)==0);
        }
    }

    // B33: alias-aware multiplication is selected only when two
    // references denote the exact same Integer. A distinct copy with
    // identical value must still produce the same GMP-oracle square.
    // Include the self-assignment x*=x contract, negative signs,
    // zero, limb saturation and the 7/8-limb dispatch boundary.
    std::mt19937_64 alias_square_engine(0xB33A11A5ULL);
    for (const std::size_t count:
            {1u,7u,8u,9u,16u,17u,32u,64u,128u}) {
        for (unsigned mode=0u;mode<16u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=alias_square_engine();
            switch (mode%8u) {
            case 0u:
                std::fill(words.begin(),words.end(),0u);
                break;
            case 1u:
                std::fill(words.begin(),words.end(),
                          std::numeric_limits<std::uint64_t>::max());
                break;
            case 2u:
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
                break;
            case 3u:
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()|=1u;
                break;
            case 4u:
                for (std::size_t i=1u;i+1u<count;++i) words[i]=0u;
                break;
            default:
                words.back()|=std::uint64_t{1}<<63u;
                break;
            }
            mpz_import(a.value,words.size(),-1,sizeof(std::uint64_t),
                       0,0,words.data());
            if (mode%2u!=0u) mpz_neg(a.value,a.value);
            const Integer x=integer_from_gmp(a.value);
            const Integer distinct_copy=x;
            mpz_mul(expected.value,a.value,a.value);

            const GmpIntegerBoundary aliased_product(x*x);
            const GmpIntegerBoundary distinct_product(x*distinct_copy);
            ARIADNE_TEST_ASSERT(
                mpz_cmp(aliased_product.value(),expected.value)==0);
            ARIADNE_TEST_ASSERT(
                mpz_cmp(distinct_product.value(),expected.value)==0);

            Integer self_assigned=x;
            self_assigned*=self_assigned;
            const GmpIntegerBoundary assigned_product(self_assigned);
            ARIADNE_TEST_ASSERT(
                mpz_cmp(assigned_product.value(),expected.value)==0);
            ARIADNE_TEST_ASSERT(!self_assigned.negative());
        }
    }

    // B37: large signed powers exercise the same ARM64 specialized
    // square used by pow3 without changing the exponentiation schedule.
    std::mt19937_64 arm_pow_engine(0xB37C03ULL);
    for (const std::size_t count: {7u,8u,9u,16u,32u,64u}) {
        for (unsigned repetition=0u;repetition<12u;++repetition) {
            std::vector<std::uint64_t> digits(count);
            for (auto& word:digits) word=arm_pow_engine();
            if (repetition%5u==0u) {
                std::fill(digits.begin(),digits.end(),
                          std::numeric_limits<std::uint64_t>::max());
            } else if (repetition%5u==1u) {
                std::fill(digits.begin(),digits.end(),0u);
                digits.back()=std::uint64_t{1}<<63u;
            } else {
                digits.back()|=std::uint64_t{1}<<63u;
            }
            mpz_import(a.value,digits.size(),-1,sizeof(std::uint64_t),
                       0,0,digits.data());
            if (repetition%2u!=0u) mpz_neg(a.value,a.value);
            const Integer input=integer_from_gmp(a.value);
            for (const unsigned exponent: {2u,3u,5u}) {
                mpz_pow_ui(expected.value,a.value,exponent);
                const GmpIntegerBoundary native_pow(
                    pow(input,Nat(exponent)));
                ARIADNE_TEST_ASSERT(
                    mpz_cmp(native_pow.value(),expected.value)==0);
            }
        }
    }

    // B38: exercise pow3 on general zero-row patterns and both signs.
    // These are genuine arbitrary-precision powers, not a shortcut
    // based on the benchmark's particular operand construction.
    std::mt19937_64 sparse_power_engine(0xB38C03ULL);
    for (const std::size_t count: {8u,16u,32u,64u,128u}) {
        for (unsigned repetition=0u;repetition<24u;++repetition) {
            std::vector<std::uint64_t> digits(count);
            for (auto& word:digits) word=sparse_power_engine();
            if (repetition%6u==0u) {
                std::fill(digits.begin(),digits.end(),0u);
                digits.back()=std::uint64_t{1}<<63u;
            } else if (repetition%6u==1u) {
                std::fill(digits.begin(),digits.end(),0u);
                digits.front()=13u;
                digits.back()=std::uint64_t{1}<<63u;
            } else if (repetition%6u==2u) {
                for (std::size_t i=0u;i<count;i+=2u) digits[i]=0u;
            } else if (repetition%6u==3u) {
                std::fill(digits.begin(),digits.end(),
                          std::numeric_limits<std::uint64_t>::max());
            }
            digits.back()|=std::uint64_t{1}<<63u;
            mpz_import(a.value,digits.size(),-1,sizeof(std::uint64_t),
                       0,0,digits.data());
            if (repetition%2u!=0u) mpz_neg(a.value,a.value);
            const Integer input=integer_from_gmp(a.value);
            mpz_pow_ui(expected.value,a.value,3u);
            const GmpIntegerBoundary native_cube(pow(input,Nat(3u)));
            ARIADNE_TEST_ASSERT(
                mpz_cmp(native_cube.value(),expected.value)==0);
        }
    }

    std::mt19937_64 engine(0x1A7E6EEDULL);
    for (unsigned repetition=0u;repetition<1700u;++repetition) {
        std::string lhs=std::to_string(engine());
        std::string rhs=std::to_string(engine()|std::uint64_t{1});
        for (unsigned j=0u;j<repetition%7u;++j) {
            lhs+=std::to_string(engine());
            rhs+=std::to_string(engine());
        }
        if (repetition%2u!=0u) lhs.insert(lhs.begin(),'-');
        if (repetition%3u!=0u) rhs.insert(rhs.begin(),'-');
        if (repetition%19u==0u) lhs="0";
        const Integer x{String(lhs)};
        const Integer y{String(rhs)};
        ARIADNE_TEST_ASSERT(mpz_set_str(a.value,lhs.c_str(),10)==0);
        ARIADNE_TEST_ASSERT(mpz_set_str(b.value,rhs.c_str(),10)==0);
        mpz_add(expected.value,a.value,b.value);
        GmpIntegerBoundary sum(x+y);
        ARIADNE_TEST_ASSERT(mpz_cmp(sum.value(),expected.value)==0);
        mpz_sub(expected.value,a.value,b.value);
        GmpIntegerBoundary difference(x-y);
        ARIADNE_TEST_ASSERT(mpz_cmp(difference.value(),expected.value)==0);
        mpz_mul(expected.value,a.value,b.value);
        GmpIntegerBoundary product(x*y);
        ARIADNE_TEST_ASSERT(mpz_cmp(product.value(),expected.value)==0);
        mpz_tdiv_qr(quotient.value,remainder.value,a.value,b.value);
        const Integer q=quot(x,y);
        const Integer r=rem(x,y);
        GmpIntegerBoundary native_q(q),native_r(r);
        ARIADNE_TEST_ASSERT(mpz_cmp(native_q.value(),quotient.value)==0);
        ARIADNE_TEST_ASSERT(mpz_cmp(native_r.value(),remainder.value)==0);
        ARIADNE_TEST_ASSERT(q*y+r==x);
        ARIADNE_TEST_ASSERT(abs(r)<abs(y));
        const auto converted=integer_from_gmp(a.value);
        ARIADNE_TEST_ASSERT(converted==x);
        ARIADNE_TEST_ASSERT(converted.literal()==x.literal());
    }
}

// Preserve all B17 signed-magnitude regression scenarios on the
// production Integer type. The earlier B18 differential corpus (1700
// independent cases) is kept in test_native() above.
void TestInteger::test_signed_regression() {
    const Integer zero(0);
    const Integer one(1);
    const Integer minus_one(-1);
    ARIADNE_TEST_ASSERT(sgn(zero)==Sign::ZERO);
    ARIADNE_TEST_ASSERT(!Integer(String("-0")).negative());
    ARIADNE_TEST_ASSERT(Integer(String("-00000")).literal()=="0");
    ARIADNE_TEST_ASSERT(Integer(String("+000123")).literal()=="123");
    ARIADNE_TEST_ASSERT(Integer::from_magnitude(BigUnsigned{},true)==zero);
    ARIADNE_TEST_ASSERT(!Integer::from_magnitude(BigUnsigned{},true).negative());
    ARIADNE_TEST_ASSERT(cmp(one,minus_one)==Comparison::GREATER);
    ARIADNE_TEST_ASSERT(cmp(minus_one,zero)==Comparison::LESS);
    ARIADNE_TEST_ASSERT(cmp(zero,zero)==Comparison::EQUAL);
    ARIADNE_TEST_ASSERT(Integer(std::numeric_limits<std::int64_t>::min())
        .get<std::int64_t>()==std::numeric_limits<std::int64_t>::min());
    ARIADNE_TEST_ASSERT(Integer(std::numeric_limits<std::uint64_t>::max())
        .get<std::uint64_t>()==std::numeric_limits<std::uint64_t>::max());
    ARIADNE_TEST_ASSERT(pow(zero,Nat(0u)).literal()=="1");
    ARIADNE_TEST_ASSERT(sqr(Integer(-13))==169);
    ARIADNE_TEST_ASSERT(quot(Integer(-7),Integer(2))==-3);

    for (const char* invalid: {"","-","+"," 1","1 ","1_0","--1"}) {
        bool threw=false;
        try { (void)Integer(String(invalid)); }
        catch (std::invalid_argument const&) { threw=true; }
        ARIADNE_TEST_ASSERT(threw);
    }
    bool negative_to_unsigned=false,positive_overflow=false;
    try { (void)minus_one.get<std::uint64_t>(); }
    catch (std::out_of_range const&) { negative_to_unsigned=true; }
    try { (void)Integer(String("9223372036854775808"))
                  .get<std::int64_t>(); }
    catch (std::out_of_range const&) { positive_overflow=true; }
    ARIADNE_TEST_ASSERT(negative_to_unsigned && positive_overflow);

    for (const char* left: {"-11","11","-12","12","0"}) {
        for (const char* right: {"-3","3","-5","5"}) {
            const Integer a{String(left)},b{String(right)};
            const Integer q=quot(a,b),r=rem(a,b);
            ARIADNE_TEST_ASSERT(cmp(q*b+r,a)==Comparison::EQUAL);
            ARIADNE_TEST_ASSERT(abs(r)<abs(b));
            ARIADNE_TEST_ASSERT(is_zero(r) || sgn(r)==sgn(a));
        }
    }

    struct GMP {
        mpz_t value;
        GMP() { mpz_init(value); }
        ~GMP() { mpz_clear(value); }
        GMP(GMP const&)=delete;
        GMP& operator=(GMP const&)=delete;
    };
    GMP a,b,expected,quotient,remainder;
    const auto check=[&](Integer const& actual,mpz_srcptr oracle) {
        GmpIntegerBoundary converted(actual);
        ARIADNE_TEST_ASSERT(mpz_cmp(converted.value(),oracle)==0);
        ARIADNE_TEST_ASSERT(!is_zero(actual) || !actual.negative());
        ARIADNE_TEST_ASSERT((mpz_sgn(oracle)==0)==is_zero(actual));
        ARIADNE_TEST_ASSERT((mpz_sgn(oracle)<0)==actual.negative());
    };

    // Retain the full old 2500-case corpus (seed, digit patterns and
    // signed divisors) in addition to B18's existing 1700-case suite.
    std::mt19937_64 engine(0x51C0FFEEULL);
    for (unsigned rep=0u;rep<2500u;++rep) {
        std::string ad=std::to_string(engine());
        std::string bd=std::to_string(engine()|std::uint64_t{1});
        for (unsigned j=0u;j<rep%8u;++j) {
            ad+=std::to_string(engine());
            bd+=std::to_string(engine());
        }
        if (rep%2u!=0u) ad.insert(ad.begin(),'-');
        if (rep%3u!=0u) bd.insert(bd.begin(),'-');
        if (rep%17u==0u) ad="0";

        const Integer x{String(ad)},y{String(bd)};
        ARIADNE_TEST_ASSERT(mpz_set_str(a.value,ad.c_str(),10)==0);
        ARIADNE_TEST_ASSERT(mpz_set_str(b.value,bd.c_str(),10)==0);
        check(x,a.value);
        check(y,b.value);
        mpz_add(expected.value,a.value,b.value);
        check(x+y,expected.value);
        mpz_sub(expected.value,a.value,b.value);
        check(x-y,expected.value);
        mpz_mul(expected.value,a.value,b.value);
        check(x*y,expected.value);

        mpz_tdiv_qr(quotient.value,remainder.value,a.value,b.value);
        const Integer q=quot(x,y),r=rem(x,y);
        check(q,quotient.value);
        check(r,remainder.value);
        ARIADNE_TEST_ASSERT(q*y+r==x);
        ARIADNE_TEST_ASSERT(abs(r)<abs(y));
        ARIADNE_TEST_ASSERT(is_zero(r) || sgn(r)==sgn(x));

        const unsigned bits=rep%130u;
        const BigUnsigned scale=shift_left(
            BigUnsigned(std::uint64_t{1}),bits);
        const Integer divisor=Integer::from_magnitude(scale,false);

        // The former temporary signed right shift truncated toward zero,
        // unlike an arithmetic/floor shift for negative inputs.
        // Integer has no signed-shift API: verify the same contract
        // through its native magnitude shift and signed quotient.
        mpz_tdiv_q_2exp(expected.value,a.value,
                        static_cast<mp_bitcnt_t>(bits));
        const Integer shifted_right=Integer::from_magnitude(
            shift_right(x.magnitude(),bits),x.negative());
        check(shifted_right,expected.value);
        check(quot(x,divisor),expected.value);

        // The migrated signed left shift keeps the source sign.
        mpz_mul_2exp(expected.value,a.value,
                     static_cast<mp_bitcnt_t>(bits));
        const Integer shifted_left=Integer::from_magnitude(
            shift_left(x.magnitude(),bits),x.negative());
        check(shifted_left,expected.value);
    }

    bool quotient_zero=false,remainder_zero=false;
    try { (void)quot(one,zero); }
    catch (std::domain_error const&) { quotient_zero=true; }
    try { (void)rem(one,zero); }
    catch (std::domain_error const&) { remainder_zero=true; }
    ARIADNE_TEST_ASSERT(quotient_zero && remainder_zero);
}

int main() {
    ARIADNE_TEST_CLASS(Integer,TestInteger());

    return ARIADNE_TEST_FAILURES;
}
