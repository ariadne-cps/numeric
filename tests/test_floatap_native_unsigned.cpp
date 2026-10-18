/***************************************************************************
 *            test_floatap_native_unsigned.cpp
 *
 *  Copyright  2026  Ariadne contributors
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

/*! \file test_floatap_native_unsigned.cpp
 *  \brief Exact unsigned limb invariants and secondary GMP differential tests.
 */

#include "floatap_native_unsigned.hpp"
#include "floatap_native_ratio.hpp"
#include "macro/testing.hpp"

#include <gmp.h>

#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace Ariadne;
using namespace Ariadne::FloatAPNative;

namespace {

struct GmpInteger {
    mpz_t value;
    GmpInteger() { mpz_init(value); }
    ~GmpInteger() { mpz_clear(value); }
    GmpInteger(GmpInteger const&)=delete;
    GmpInteger& operator=(GmpInteger const&)=delete;
};

void import_exact(mpz_ptr output,Unsigned const& input) {
    mpz_import(output,input.size(),-1,sizeof(std::uint32_t),0,0,input.data());
}

void assert_exact(Unsigned const& input,mpz_srcptr expected) {
    GmpInteger actual;
    import_exact(actual.value,input);
    ARIADNE_TEST_ASSERT(mpz_cmp(actual.value,expected)==0);
    ARIADNE_TEST_ASSERT(input.is_zero() == (mpz_sgn(expected)==0));
    if (!input.is_zero()) {
        ARIADNE_TEST_ASSERT(input.word(input.size()-1u)!=0u);
        ARIADNE_TEST_ASSERT(input.bit(input.bit_length()-1u));
    }
}

Unsigned make_value(unsigned bits,std::uint64_t seed) {
    Unsigned result;
    std::uint64_t state=seed;
    for (unsigned i=0u; i<bits; i+=32u) {
        state=state*6364136223846793005ull+1442695040888963407ull;
        const std::uint32_t word=static_cast<std::uint32_t>(state>>32);
        const unsigned count=std::min(32u,bits-i);
        const std::uint32_t mask=count==32u ? 0xffffffffu
            : (std::uint32_t{1}<<count)-1u;
        result.set_word(i/32u,word&mask);
    }
    if (bits!=0u) { result.set_bit(bits-1u); }
    return result;
}

void test_small_oracles() {
    ARIADNE_TEST_ASSERT(Unsigned{}.is_zero());
    ARIADNE_TEST_ASSERT(Unsigned(1u).bit_length()==1u);
    ARIADNE_TEST_ASSERT(Unsigned(0xffffffffu).bit_length()==32u);
    ARIADNE_TEST_ASSERT(Unsigned(0x100000000ull).bit_length()==33u);
    ARIADNE_TEST_ASSERT(Unsigned(0x100000000ull).word(1u)==1u);
    ARIADNE_TEST_ASSERT(compare(Unsigned(7u),Unsigned(8u))<0);
    ARIADNE_TEST_ASSERT(add(Unsigned(0xffffffffu),Unsigned(1u)) ==
                        Unsigned(0x100000000ull));
    ARIADNE_TEST_ASSERT(sub(Unsigned(0x100000000ull),Unsigned(1u)) ==
                        Unsigned(0xffffffffu));
    ARIADNE_TEST_ASSERT(mul(Unsigned(0xffffffffu),Unsigned(0xffffffffu)) ==
                        Unsigned(0xfffffffe00000001ull));
    ARIADNE_TEST_ASSERT(shift_left(Unsigned(3u),32u) ==
                        Unsigned(0x300000000ull));
    ARIADNE_TEST_ASSERT(shift_right(Unsigned(0x300000000ull),32u) ==
                        Unsigned(3u));
    ARIADNE_TEST_ASSERT(divrem(Unsigned(17u),Unsigned(5u)).first ==
                        Unsigned(3u));
    ARIADNE_TEST_ASSERT(divrem(Unsigned(17u),Unsigned(5u)).second ==
                        Unsigned(2u));
    ARIADNE_TEST_ASSERT(divrem(Unsigned(3u),Unsigned(7u)).second ==
                        Unsigned(3u));
    ARIADNE_TEST_ASSERT(divrem(Unsigned(7u),Unsigned(7u)).first ==
                        Unsigned(1u));
    ARIADNE_TEST_ASSERT(shift_left(Unsigned{},10000u).is_zero());
    ARIADNE_TEST_ASSERT(shift_right(Unsigned(3u),10000u).is_zero());

    Unsigned boundary;
    boundary.set_bit(1024u);
    ARIADNE_TEST_ASSERT(boundary.bit_length()==1025u);
    ARIADNE_TEST_ASSERT(boundary.size()==33u);
    ARIADNE_TEST_ASSERT(sub(boundary,Unsigned(1u)).bit_length()==1024u);

    bool underflow=false;
    try { (void)sub(Unsigned(1u),Unsigned(2u)); }
    catch (std::domain_error const&) { underflow=true; }
    ARIADNE_TEST_ASSERT(underflow);
    bool zero_division=false;
    try { (void)divrem(Unsigned(1u),Unsigned{}); }
    catch (std::domain_error const&) { zero_division=true; }
    ARIADNE_TEST_ASSERT(zero_division);
}

void test_arithmetic_differential() {
    GmpInteger a,b,expected,quotient,remainder;
    const unsigned sizes[]={1u,2u,31u,32u,33u,53u,64u,65u,
                            127u,128u,129u,256u,512u,1024u,2048u};
    for (unsigned bits : sizes) {
        for (unsigned index=0u; index<8u; ++index) {
            const Unsigned x=make_value(bits,17u+index);
            const Unsigned y=make_value(bits/2u+1u,211u+index);
            import_exact(a.value,x);
            import_exact(b.value,y);

            mpz_add(expected.value,a.value,b.value);
            assert_exact(add(x,y),expected.value);
            mpz_mul(expected.value,a.value,b.value);
            assert_exact(mul(x,y),expected.value);
            mpz_mul(expected.value,a.value,a.value);
            assert_exact(mul(x,x),expected.value);

            mpz_sub(expected.value,a.value,b.value);
            assert_exact(sub(x,y),expected.value);
            ARIADNE_TEST_ASSERT(compare(x,y)>=0);

            for (unsigned shift : {0u,1u,31u,32u,33u,65u,257u}) {
                mpz_mul_2exp(expected.value,a.value,shift);
                assert_exact(shift_left(x,shift),expected.value);
                mpz_fdiv_q_2exp(expected.value,a.value,shift);
                assert_exact(shift_right(x,shift),expected.value);
            }
            mpz_tdiv_qr(quotient.value,remainder.value,a.value,b.value);
            const auto answer=divrem(x,y);
            assert_exact(answer.first,quotient.value);
            assert_exact(answer.second,remainder.value);
            ARIADNE_TEST_ASSERT(compare(answer.second,y)<0);
            ARIADNE_TEST_ASSERT(add(mul(answer.first,y),answer.second)==x);
        }
    }
}

} // namespace


void test_ratio_oracles() {
    for (unsigned long bits : {2ul,53ul,128ul,1024ul,1100ul,1536ul}) {
        const Unsigned numerator(1u),denominator(3u);
        for (bool negative : {false,true}) {
            for (auto mode : {RatioRounding::NearestEven,
                              RatioRounding::Downward,
                              RatioRounding::Upward,
                              RatioRounding::TowardZero}) {
                const auto result=round_ratio(numerator,denominator,bits,negative,mode);
                const Unsigned lower=divrem(
                    shift_left(numerator,static_cast<std::size_t>(bits+1ul)),
                    denominator).first;
                // 1/3 lies in [1/4,1/2): k=-2, exponent=-p-1.
                // The quotient is floor(2^(p+1)/3); its discarded
                // fraction is 1/3 or 2/3, so ties never occur.
                const auto quotient=divrem(
                    shift_left(Unsigned(1u),static_cast<std::size_t>(bits+1ul)),
                    Unsigned(3u));
                const bool round_up=mode==RatioRounding::NearestEven
                    ? compare(shift_left(quotient.second,1u),Unsigned(3u))>0
                    : mode==RatioRounding::Downward ? negative
                    : mode==RatioRounding::Upward ? !negative : false;
                const Unsigned expected=round_up
                    ? add(lower,Unsigned(1u)) : lower;
                ARIADNE_TEST_ASSERT(result.significand==expected);
                ARIADNE_TEST_ASSERT(result.exponent==
                                    -static_cast<std::int64_t>(bits)-1);
                ARIADNE_TEST_ASSERT(!result.carry);
            }
        }
    }
    // Exact half-integer quotients: ties must round to even.
    const auto even=round_ratio(Unsigned(5u),Unsigned(4u),2u,false,
                                RatioRounding::NearestEven);
    ARIADNE_TEST_ASSERT(even.significand==Unsigned(2u));
    ARIADNE_TEST_ASSERT(even.exponent==-1);
    const auto odd=round_ratio(Unsigned(7u),Unsigned(4u),2u,false,
                               RatioRounding::NearestEven);
    ARIADNE_TEST_ASSERT(odd.significand==Unsigned(2u));
    ARIADNE_TEST_ASSERT(odd.exponent==0);
    ARIADNE_TEST_ASSERT(odd.carry);
}

int main() {
    ARIADNE_TEST_CALL(test_small_oracles());
    ARIADNE_TEST_CALL(test_ratio_oracles());
    ARIADNE_TEST_CALL(test_arithmetic_differential());
    return ARIADNE_TEST_FAILURES;
}
