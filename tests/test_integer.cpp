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
#include "numeric/int.hpp"
#include "numeric/twoexp.hpp"
#include "numeric/accuracy.hpp"
#include "numeric/sign.hpp"
#include "foundation/logical.hpp"
#include "utility/string.hpp"

#include <algorithm>
#include <array>
#include <span>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>


inline constexpr char ARIADNE_INTEGER_SUB_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_sub.tsv";

inline constexpr char ARIADNE_INTEGER_B17_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_b17.tsv";
inline constexpr char ARIADNE_INTEGER_B18_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_b18.tsv";

inline constexpr char ARIADNE_INTEGER_B33_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_b33.tsv";
inline constexpr char ARIADNE_INTEGER_B37_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_b37.tsv";
inline constexpr char ARIADNE_INTEGER_B38_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_b38.tsv";

inline constexpr char ARIADNE_INTEGER_B66_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_b66.tsv";

#ifndef ARIADNE_TEST_DATA_DIR
#error "Provide the absolute test data directory from CMake"
#endif

inline constexpr char ARIADNE_INTEGER_POWER_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_power.tsv";
inline constexpr char ARIADNE_INTEGER_SQRT_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_sqrt.tsv";

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

// Parse the independently frozen hexadecimal oracle values without GMP.
BigUnsigned sqrt_oracle_hex(std::string const& hex)
{
    if (hex.empty() || (hex.size()>1u && hex.front()=='0')) {
        throw std::invalid_argument("Invalid Integer square-root oracle hex");
    }
    std::vector<std::uint64_t> words((hex.size()+15u)/16u,0u);
    for (std::size_t i=0u;i<hex.size();++i) {
        const char digit=hex[hex.size()-1u-i];
        const unsigned nibble=digit>='0' && digit<='9'
            ? static_cast<unsigned>(digit-'0')
            : digit>='a' && digit<='f'
                ? static_cast<unsigned>(digit-'a'+10)
                : 16u;
        if (nibble==16u) {
            throw std::invalid_argument("Malformed Integer square-root oracle hex");
        }
        words[i/16u]|=static_cast<std::uint64_t>(nibble)
                       <<static_cast<unsigned>(4u*(i%16u));
    }
    return BigUnsigned(words);
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
    void test_owned_results();
    void test_sqrt();
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
    ARIADNE_TEST_CALL(test_owned_results());
    ARIADNE_TEST_CALL(test_sqrt());
}

void TestInteger::test_sqrt()
{
    // Reuse the 161 independent restoring-digit sqrt/root/remainder
    // triples, not GMP or another call to the function under test.
    std::ifstream golden(ARIADNE_INTEGER_SQRT_ORACLE_PATH);
    if (!golden) throw std::runtime_error("Missing Integer sqrt oracle TSV");
    std::size_t cases=0u;
    std::string line;
    while (std::getline(golden,line)) {
        if (line.empty() || line.front()=='#') continue;
        std::istringstream fields(line);
        std::string input_hex,root_hex,remainder_hex,extra;
        if (!(fields>>input_hex>>root_hex>>remainder_hex) || (fields>>extra)) {
            throw std::runtime_error("Malformed Integer sqrt oracle row");
        }
        const BigUnsigned input_magnitude=sqrt_oracle_hex(input_hex);
        const Integer input=Integer::from_magnitude(input_magnitude,false);
        const Integer preserved=input;
        const Integer expected_root=Integer::from_magnitude(
            sqrt_oracle_hex(root_hex),false);
        const Integer expected_remainder=Integer::from_magnitude(
            sqrt_oracle_hex(remainder_hex),false);
        const Integer root=isqrt(input);
        const auto [paired_root,remainder]=sqrtrem(input);
        ARIADNE_TEST_ASSERT(root==expected_root);
        ARIADNE_TEST_ASSERT(paired_root==expected_root);
        ARIADNE_TEST_ASSERT(remainder==expected_remainder);
        ARIADNE_TEST_ASSERT(!root.negative() && !remainder.negative());
        ARIADNE_TEST_ASSERT(!root.negative() && !paired_root.negative());
        ARIADNE_TEST_ASSERT(paired_root*paired_root+remainder==input);
        ARIADNE_TEST_ASSERT(remainder<=paired_root+paired_root);
        ARIADNE_TEST_ASSERT(input==preserved);

        // Only zero accepts a canonicalized negative-zero source.
        const Integer negative_input=Integer::from_magnitude(
            input_magnitude,true);
        if (input_magnitude.is_zero()) {
            ARIADNE_TEST_ASSERT(!negative_input.negative());
            ARIADNE_TEST_ASSERT(isqrt(negative_input)==Integer(0));
            const auto [zero_root,zero_remainder]=sqrtrem(negative_input);
            ARIADNE_TEST_ASSERT(zero_root==Integer(0));
            ARIADNE_TEST_ASSERT(zero_remainder==Integer(0));
        } else {
            bool root_threw=false,remainder_threw=false;
            try { (void)isqrt(negative_input); }
            catch (std::domain_error const&) { root_threw=true; }
            try { (void)sqrtrem(negative_input); }
            catch (std::domain_error const&) { remainder_threw=true; }
            ARIADNE_TEST_ASSERT(root_threw && remainder_threw);
            ARIADNE_TEST_ASSERT(negative_input.negative());
            ARIADNE_TEST_ASSERT(compare(
                negative_input.magnitude(),input_magnitude)==0);
        }
        ++cases;
    }
    if (golden.bad() || cases!=161u) {
        throw std::runtime_error("Incomplete Integer sqrt oracle fixture");
    }
}

void TestInteger::test_owned_results() {
    // Exercise inline and heap magnitude ownership independently of
    // decimal parsing. Moving a result must preserve its full value,
    // leave the caller's moved-from magnitude zero and normalize -0.
    for (const std::size_t count:{0u,1u,4u,5u,8u,9u,64u}) {
        std::vector<std::uint64_t> words(count,
            std::numeric_limits<std::uint64_t>::max());
        const BigUnsigned source(words);
        for (const bool negative:{false,true}) {
            const Integer copied=Integer::from_magnitude(source,negative);
            ARIADNE_TEST_ASSERT(compare(copied.magnitude(),source)==0);
            ARIADNE_TEST_ASSERT(copied.negative()==(negative && count!=0u));
            BigUnsigned magnitude=source;
            Integer owned=Integer::from_magnitude(std::move(magnitude),negative);
            ARIADNE_TEST_ASSERT(magnitude.is_zero());
            ARIADNE_TEST_ASSERT(owned==copied);
            const String text=owned.literal();
            ARIADNE_TEST_ASSERT(Integer(text)==copied);
            ARIADNE_TEST_ASSERT(compare(source,BigUnsigned(words))==0);
            std::ostringstream stream;
            stream << owned;
            ARIADNE_TEST_ASSERT(stream.str()==std::string(text));
            owned=Integer(7);
            ARIADNE_TEST_ASSERT(Integer(text)==copied);
        }
    }
    ARIADNE_TEST_ASSERT(Integer::from_magnitude(BigUnsigned{},true).literal()=="0");
    ARIADNE_TEST_ASSERT(Integer(-1).literal()=="-1");
    ARIADNE_TEST_ASSERT(Integer(1).literal()=="1");
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

    // B41: six structural classes x four sign pairs x eight widths.
    // Source and expected values are independently frozen, with no GMP
    // calls in this differential test. Self-aliasing remains covered.
    std::ifstream sub_fixture(ARIADNE_INTEGER_SUB_ORACLE_PATH);
    if (!sub_fixture) throw std::runtime_error("Missing Integer subtraction TSV");
    std::string sub_header;
    while (sub_fixture.peek()=='#') std::getline(sub_fixture,sub_header);
    for (const std::size_t width:{1u,2u,7u,8u,16u,32u,64u,128u}) {
        for (unsigned mode=0u;mode<6u;++mode) {
            for (unsigned signs=0u;signs<4u;++signs) {
                std::size_t ref_width=0u;
                unsigned ref_mode=0u,ref_signs=0u;
                std::string raw_lhs,raw_rhs,raw_expected;
                if (!(sub_fixture >> ref_width >> ref_mode >> ref_signs >>
                                  raw_lhs >> raw_rhs >> raw_expected) ||
                    ref_width!=width || ref_mode!=mode || ref_signs!=signs) {
                    throw std::runtime_error("Missing or unordered Integer subtraction TSV row");
                }
                const Integer lhs{String(raw_lhs)},rhs{String(raw_rhs)};
                const Integer expected_sub{String(raw_expected)};
                ARIADNE_TEST_EQUALS(lhs-rhs,expected_sub);
                Integer assigned=lhs;
                assigned-=rhs;
                ARIADNE_TEST_EQUALS(assigned,expected_sub);
                Integer aliased=lhs;
                Integer const& same_object=aliased;
                aliased-=same_object;
                ARIADNE_TEST_ASSERT(aliased==0);
                ARIADNE_TEST_ASSERT(!aliased.negative());
            }
        }
    }
    std::string sub_extra;
    if (sub_fixture >> sub_extra)
        throw std::runtime_error("Extra Integer subtraction TSV rows");

    // B32: independent frozen exact signed-power reference; retain all
    // sign/width cases and sample exponent bit patterns at boundaries.
    std::ifstream pow_fixture(ARIADNE_INTEGER_POWER_ORACLE_PATH);
    if (!pow_fixture) throw std::runtime_error("Missing Integer power TSV");
    std::string pow_header;
    while (pow_fixture.peek()=='#') std::getline(pow_fixture,pow_header);
    for (const char* literal:
            {"0","1","-1","2","-2","17","-17",
             "18446744073709551616","-18446744073709551616",
             "340282366920938463463374607431768211455"}) {
        const Integer input{String(literal)};
        for (unsigned exponent: {0u,1u,2u,3u,5u,7u,8u,9u,11u,12u}) {
            std::string ref_base,ref_result;
            unsigned ref_exponent=0u;
            if (!(pow_fixture >> ref_base >> ref_exponent >> ref_result) ||
                ref_base!=literal || ref_exponent!=exponent) {
                throw std::runtime_error("Missing or unordered Integer power TSV row");
            }
            ARIADNE_TEST_EQUALS(pow(input,Nat(exponent)),
                                Integer(String(ref_result)));
        }
    }
    std::string pow_extra;
    if (pow_fixture >> pow_extra) {
        throw std::runtime_error("Extra Integer power TSV rows");
    }

    // B66: original 18 exponent boundaries, across each original width.
    // Four distinct input-pattern classes are retained from the source
    // mt19937_64 sequence. Results were frozen using independent exact
    // integer arithmetic (no GMP or Integer-under-test at runtime).
    std::ifstream exponent_oracle(ARIADNE_INTEGER_B66_ORACLE_PATH);
    if (!exponent_oracle) throw std::runtime_error("Missing Integer B66 TSV");
    std::string exponent_header;
    while (exponent_oracle.peek()=='#')
        std::getline(exponent_oracle,exponent_header);
    for (const std::size_t width:{0u,1u,2u,3u,4u}) {
        for (unsigned mode=0u;mode<4u;++mode) {
            for (unsigned exponent:
                    {0u,1u,2u,3u,5u,7u,8u,9u,12u,15u,16u,17u,
                     31u,32u,33u,63u,64u,65u}) {
                std::size_t ref_width=0u;
                unsigned ref_mode=0u,ref_exponent=0u;
                std::string raw_base,raw_result;
                if (!(exponent_oracle >> ref_width >> ref_mode >> ref_exponent
                        >> raw_base >> raw_result)
                    || ref_width!=width || ref_mode!=mode
                    || ref_exponent!=exponent) {
                    throw std::runtime_error("Missing or unordered B66 TSV row");
                }
                const Integer base{String(raw_base)};
                ARIADNE_TEST_EQUALS(pow(base,Nat(exponent)),
                                    Integer(String(raw_result)));
            }
        }
    }
    std::string extra_exponent;
    if (exponent_oracle >> extra_exponent)
        throw std::runtime_error("Extra B66 TSV rows");

    // B33: independent exact source-keyed squares for all eight distinct
    // structural classes at every original width; includes x*=x aliasing.
    std::ifstream square_oracle(ARIADNE_INTEGER_B33_ORACLE_PATH);
    if (!square_oracle) throw std::runtime_error("Missing Integer B33 TSV");
    std::string square_comment;
    while (square_oracle.peek()=='#') std::getline(square_oracle,square_comment);
    for (const std::size_t width:{1u,7u,8u,9u,16u,17u,32u,64u,128u}) {
        for (unsigned mode=0u;mode<8u;++mode) {
            std::size_t ref_width=0u;
            unsigned ref_mode=0u;
            std::string source,expected;
            if (!(square_oracle >> ref_width >> ref_mode >> source >> expected)
                || ref_width!=width || ref_mode!=mode) {
                throw std::runtime_error("Missing or unordered B33 TSV row");
            }
            const Integer x{String(source)},copy=x;
            const Integer square{String(expected)};
            ARIADNE_TEST_EQUALS(x*x,square);
            ARIADNE_TEST_EQUALS(x*copy,square);
            Integer assigned=x;
            assigned*=assigned;
            ARIADNE_TEST_EQUALS(assigned,square);
            ARIADNE_TEST_ASSERT(!assigned.negative());
        }
    }
    std::string extra_square;
    if (square_oracle >> extra_square)
        throw std::runtime_error("Extra B33 TSV rows");

    // B37: keep the saturation/high-bit/dense sign classes across
    // the 7/8/9-limb square dispatch and all original wider widths.
    std::ifstream large_power_oracle(ARIADNE_INTEGER_B37_ORACLE_PATH);
    if (!large_power_oracle) throw std::runtime_error("Missing Integer B37 TSV");
    std::string large_power_comment;
    while (large_power_oracle.peek()=='#')
        std::getline(large_power_oracle,large_power_comment);
    for (const std::size_t width:{7u,8u,9u,16u,32u,64u}) {
        for (unsigned repetition=0u;repetition<4u;++repetition) {
            for (unsigned exponent:{2u,3u,5u}) {
                std::size_t ref_width=0u;
                unsigned ref_repetition=0u,ref_exponent=0u;
                std::string source,expected;
                if (!(large_power_oracle >> ref_width >> ref_repetition
                        >> ref_exponent >> source >> expected)
                    || ref_width!=width || ref_repetition!=repetition
                    || ref_exponent!=exponent) {
                    throw std::runtime_error("Missing or unordered B37 TSV row");
                }
                const Integer x{String(source)};
                ARIADNE_TEST_EQUALS(pow(x,Nat(exponent)),
                                    Integer(String(expected)));
            }
        }
    }
    std::string extra_large_power;
    if (large_power_oracle >> extra_large_power)
        throw std::runtime_error("Extra B37 TSV rows");

    // B38: all six sparse/dense structural cases at 8..128 limbs,
    // including the alternating-zero, saturated and both sign classes.
    std::ifstream cubic_oracle(ARIADNE_INTEGER_B38_ORACLE_PATH);
    if (!cubic_oracle) throw std::runtime_error("Missing Integer B38 TSV");
    std::string cubic_comment;
    while (cubic_oracle.peek()=='#') std::getline(cubic_oracle,cubic_comment);
    for (const std::size_t width:{8u,16u,32u,64u,128u}) {
        for (unsigned repetition=0u;repetition<6u;++repetition) {
            std::size_t ref_width=0u;
            unsigned ref_repetition=0u;
            std::string source,expected;
            if (!(cubic_oracle >> ref_width >> ref_repetition
                 >> source >> expected)
                || ref_width!=width || ref_repetition!=repetition) {
                throw std::runtime_error("Missing or unordered B38 TSV row");
            }
            ARIADNE_TEST_EQUALS(pow(Integer(String(source)),Nat(3u)),
                                Integer(String(expected)));
        }
    }
    std::string extra_cubic;
    if (cubic_oracle >> extra_cubic)
        throw std::runtime_error("Extra B38 TSV rows");

    // B18: frozen exact signed arithmetic on 170 strategically sampled
    // originals out of 1700; expected results never use GMP at runtime.
    std::ifstream golden(ARIADNE_INTEGER_B18_ORACLE_PATH);
    if (!golden) throw std::runtime_error("Missing Integer B18 oracle TSV");
    std::string comment;
    while (golden.peek()=='#') std::getline(golden,comment);
    for (unsigned sample=0u;sample<170u;++sample) {
        unsigned rep=0u;
        std::string lhs,rhs,sum,sub,mul,q_expected,r_expected;
        if (!(golden >> rep >> lhs >> rhs >> sum >> sub >> mul
                     >> q_expected >> r_expected) ||
            rep!=10u*sample+(7u*sample+3u)%10u) {
            throw std::runtime_error("Missing or out-of-order Integer B18 oracle");
        }
        const Integer x{String(lhs)},y{String(rhs)};
        ARIADNE_TEST_EQUALS(x+y,Integer(String(sum)));
        ARIADNE_TEST_EQUALS(x-y,Integer(String(sub)));
        ARIADNE_TEST_EQUALS(x*y,Integer(String(mul)));
        const Integer q=quot(x,y),r=rem(x,y);
        ARIADNE_TEST_EQUALS(q,Integer(String(q_expected)));
        ARIADNE_TEST_EQUALS(r,Integer(String(r_expected)));
        ARIADNE_TEST_ASSERT(q*y+r==x);
        ARIADNE_TEST_ASSERT(abs(r)<abs(y));
        ARIADNE_TEST_ASSERT(r==0 || sgn(r)==sgn(x));
        ARIADNE_TEST_ASSERT(x.literal()==Integer(String(lhs)).literal());
    }
    std::string extra;
    if (golden >> extra)
        throw std::runtime_error("Extra Integer B18 oracle rows");
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

    // B17: 250 source-keyed references, one per ten original visits.
    // All 130 shift counts survive along with signed and zero cases.
    std::ifstream golden(ARIADNE_INTEGER_B17_ORACLE_PATH);
    if (!golden) throw std::runtime_error("Missing Integer B17 oracle TSV");
    std::string comment;
    while (golden.peek()=='#') std::getline(golden,comment);
    const auto check=[&](Integer const& actual,std::string const& reference) {
        ARIADNE_TEST_EQUALS(actual,Integer(String(reference)));
        ARIADNE_TEST_ASSERT(!is_zero(actual) || !actual.negative());
        ARIADNE_TEST_ASSERT((reference=="0")==is_zero(actual));
        ARIADNE_TEST_ASSERT((reference.front()=='-')==actual.negative());
    };
    for (unsigned sample=0u;sample<250u;++sample) {
        unsigned rep=0u;
        std::string ad,bd,sum,sub,mul,q_expected,r_expected,shr,shl;
        if (!(golden >> rep >> ad >> bd >> sum >> sub >> mul >>
                     q_expected >> r_expected >> shr >> shl) ||
            rep!=10u*sample+(7u*sample+3u)%10u) {
            throw std::runtime_error("Missing or unordered Integer B17 oracle");
        }
        const Integer x{String(ad)},y{String(bd)};
        check(x,ad);
        check(y,bd);
        check(x+y,sum);
        check(x-y,sub);
        check(x*y,mul);
        const Integer q=quot(x,y),r=rem(x,y);
        check(q,q_expected);
        check(r,r_expected);
        ARIADNE_TEST_ASSERT(q*y+r==x);
        ARIADNE_TEST_ASSERT(abs(r)<abs(y));
        ARIADNE_TEST_ASSERT(is_zero(r) || sgn(r)==sgn(x));

        const unsigned bits=rep%130u;
        const BigUnsigned scale=shift_left(
            BigUnsigned(std::uint64_t{1}),bits);
        const Integer divisor=Integer::from_magnitude(scale,false);
        const Integer shifted_right=Integer::from_magnitude(
            shift_right(x.magnitude(),bits),x.negative());
        check(shifted_right,shr);
        check(quot(x,divisor),shr);
        const Integer shifted_left=Integer::from_magnitude(
            shift_left(x.magnitude(),bits),x.negative());
        check(shifted_left,shl);
    }
    std::string extra;
    if (golden >> extra)
        throw std::runtime_error("Extra Integer B17 oracle rows");

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
