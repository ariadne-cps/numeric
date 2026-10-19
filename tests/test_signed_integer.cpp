/***************************************************************************
 *            test_signed_integer.cpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "numeric/signed_integer.hpp"

#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

#include <gmp.h>
#include "macro/testing.hpp"

using namespace Ariadne;

namespace {

struct GmpInteger {
    mpz_t value;
    GmpInteger() { mpz_init(value); }
    ~GmpInteger() { mpz_clear(value); }
    GmpInteger(GmpInteger const&) = delete;
    GmpInteger& operator=(GmpInteger const&) = delete;
};

void check(SignedInteger const& actual,mpz_srcptr expected)
{
    const auto decimal=actual.decimal();
    GmpInteger native;
    ARIADNE_TEST_ASSERT(mpz_set_str(native.value,decimal.c_str(),10)==0);
    ARIADNE_TEST_ASSERT(mpz_cmp(native.value,expected)==0);
    ARIADNE_TEST_ASSERT(actual.sign()==mpz_sgn(expected));
    ARIADNE_TEST_ASSERT(!actual.is_zero() || !actual.is_negative());
}

class TestSignedInteger
{
  public:
    void test();
};

void TestSignedInteger::test()
{
    const auto zero=SignedInteger(std::int64_t{0});
    const auto one=SignedInteger(std::int64_t{1});
    const auto minus_one=SignedInteger(std::int64_t{-1});
    ARIADNE_TEST_ASSERT(zero.sign()==0);
    ARIADNE_TEST_ASSERT(!SignedInteger("-0").is_negative());
    ARIADNE_TEST_ASSERT(SignedInteger("-00000").decimal()=="0");
    ARIADNE_TEST_ASSERT(SignedInteger("+000123").decimal()=="123");
    ARIADNE_TEST_ASSERT(compare(one,minus_one)>0);
    ARIADNE_TEST_ASSERT(compare(minus_one,zero)<0);
    ARIADNE_TEST_ASSERT(compare(zero,zero)==0);
    ARIADNE_TEST_ASSERT(SignedInteger(std::numeric_limits<std::int64_t>::min())
        .get_i64()==std::numeric_limits<std::int64_t>::min());
    ARIADNE_TEST_ASSERT(SignedInteger(std::numeric_limits<std::uint64_t>::max())
        .get_u64()==std::numeric_limits<std::uint64_t>::max());
    ARIADNE_TEST_ASSERT(pow(zero,0u).decimal()=="1");
    ARIADNE_TEST_ASSERT(square(SignedInteger("-13")).decimal()=="169");
    ARIADNE_TEST_ASSERT(shift_right(SignedInteger("-7"),1u).decimal()=="-3");
    for (const auto* invalid: {"","-","+"," 1","1 ","1_0","--1"}) {
        bool threw=false;
        try { (void)SignedInteger(invalid); }
        catch (std::invalid_argument const&) { threw=true; }
        ARIADNE_TEST_ASSERT(threw);
    }
    bool negative_to_unsigned=false,positive_overflow=false;
    try { (void)minus_one.get_u64(); }
    catch (std::out_of_range const&) { negative_to_unsigned=true; }
    try { (void)SignedInteger("9223372036854775808").get_i64(); }
    catch (std::out_of_range const&) { positive_overflow=true; }
    ARIADNE_TEST_ASSERT(negative_to_unsigned && positive_overflow);

    for (const auto* left: {"-11","11","-12","12","0"}) {
        for (const auto* right: {"-3","3","-5","5"}) {
            const SignedInteger a(left),b(right);
            const auto [q,r]=divrem(a,b);
            ARIADNE_TEST_ASSERT(compare(add(mul(q,b),r),a)==0);
            ARIADNE_TEST_ASSERT(compare(r.magnitude(),b.magnitude())<0);
            ARIADNE_TEST_ASSERT(r.is_zero() || r.sign()==a.sign());
        }
    }
    std::mt19937_64 engine(0x51C0FFEEULL);
    GmpInteger a,b,expected,quotient,remainder;
    for (unsigned rep=0u;rep<2500u;++rep) {
        std::string ad=std::to_string(engine()),bd=std::to_string(engine()|1u);
        for (unsigned j=0u;j<rep%8u;++j) {
            ad+=std::to_string(engine());
            bd+=std::to_string(engine());
        }
        if (rep%2u!=0u) ad.insert(ad.begin(),'-');
        if (rep%3u!=0u) bd.insert(bd.begin(),'-');
        if (rep%17u==0u) ad="0";
        const SignedInteger x(ad),y(bd);
        ARIADNE_TEST_ASSERT(mpz_set_str(a.value,ad.c_str(),10)==0);
        ARIADNE_TEST_ASSERT(mpz_set_str(b.value,bd.c_str(),10)==0);
        check(x,a.value);
        mpz_add(expected.value,a.value,b.value);
        check(add(x,y),expected.value);
        mpz_sub(expected.value,a.value,b.value);
        check(sub(x,y),expected.value);
        mpz_mul(expected.value,a.value,b.value);
        check(mul(x,y),expected.value);
        mpz_tdiv_qr(quotient.value,remainder.value,a.value,b.value);
        auto [q,r]=divrem(x,y);
        check(q,quotient.value);
        check(r,remainder.value);
        ARIADNE_TEST_ASSERT(compare(add(mul(q,y),r),x)==0);
        const unsigned bits=rep%130u;
        const auto shifted=shift_right(x,bits);
        mpz_tdiv_q_2exp(expected.value,a.value,
                        static_cast<mp_bitcnt_t>(bits));
        check(shifted,expected.value);
    }
    bool zero_division=false;
    try { (void)divrem(one,zero); }
    catch (std::domain_error const&) { zero_division=true; }
    ARIADNE_TEST_ASSERT(zero_division);
}

} // namespace

int main()
{
    ARIADNE_TEST_CLASS(SignedInteger,TestSignedInteger());
    return ARIADNE_TEST_FAILURES;
}
