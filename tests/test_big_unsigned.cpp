/***************************************************************************
 *            test_big_unsigned.cpp
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
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "numeric/big_unsigned.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gmp.h>
#include "macro/testing.hpp"

using namespace Ariadne;

namespace {

void check_against_gmp(BigUnsigned const& x, mpz_srcptr expected)
{
    mpz_t native;
    mpz_init(native);
    const auto digits = x.limbs();
    mpz_import(native,digits.size(),-1,sizeof(BigUnsigned::Limb),0,0,digits.data());
    ARIADNE_TEST_ASSERT(mpz_cmp(native,expected)==0);
    ARIADNE_TEST_ASSERT(x.is_zero() == (mpz_sgn(native)==0));
    ARIADNE_TEST_ASSERT(x.bit_length() ==
        (mpz_sgn(native)==0 ? 0u : mpz_sizeinbase(native,2)));
    ARIADNE_TEST_ASSERT(digits.empty() || digits.back()!=0u);
    mpz_clear(native);
}

struct GmpInteger {
    mpz_t value;
    GmpInteger() { mpz_init(value); }
    ~GmpInteger() { mpz_clear(value); }
    GmpInteger(GmpInteger const&) = delete;
    GmpInteger& operator=(GmpInteger const&) = delete;
};

class TestBigUnsigned
{
  public:
    void test() {
        ARIADNE_TEST_CALL(test_storage());
        ARIADNE_TEST_CALL(test_boundaries());
        ARIADNE_TEST_CALL(test_randomized());
        ARIADNE_TEST_CALL(test_multiplication());
        ARIADNE_TEST_CALL(test_division_limb());
        ARIADNE_TEST_CALL(test_division());
        ARIADNE_TEST_CALL(test_gcd_decimal());
        ARIADNE_TEST_CALL(test_views_and_scratch());
        ARIADNE_TEST_CALL(test_multiply_kernel_paths());
    }

  private:
    void test_storage();
    void test_boundaries();
    void test_randomized();
    void test_multiplication();
    void test_division_limb();
    void test_division();
    void test_gcd_decimal();
    void test_views_and_scratch();
    void test_multiply_kernel_paths();
};

void TestBigUnsigned::test_storage()
{
    BigUnsigned zero;
    ARIADNE_TEST_ASSERT(zero.is_zero());
    ARIADNE_TEST_ASSERT(zero.limbs().empty());
    ARIADNE_TEST_ASSERT(zero.bit_length()==0u);

    BigUnsigned one(std::uint64_t{1});
    ARIADNE_TEST_ASSERT(one.limb(0)==1u);
    ARIADNE_TEST_ASSERT(one.limb(500)==0u);
    ARIADNE_TEST_ASSERT(one.bit_length()==1u);

    std::array<std::uint64_t,7> raw{1,2,3,4,5,0,0};
    BigUnsigned long_value(raw);
    ARIADNE_TEST_ASSERT(long_value.size()==5u);
    ARIADNE_TEST_ASSERT(long_value.limb(4)==5u);

    BigUnsigned copy(long_value);
    ARIADNE_TEST_ASSERT(compare(copy,long_value)==0);

    BigUnsigned moved(std::move(copy));
    ARIADNE_TEST_ASSERT(copy.is_zero());
    ARIADNE_TEST_ASSERT(compare(moved,long_value)==0);

    copy=one;
    copy=std::move(moved);
    ARIADNE_TEST_ASSERT(moved.is_zero());
    ARIADNE_TEST_ASSERT(compare(copy,long_value)==0);

    zero.swap(copy);
    ARIADNE_TEST_ASSERT(zero.size()==5u);
    ARIADNE_TEST_ASSERT(copy.is_zero());
    const std::array<std::uint64_t,5> zeros{};
    ARIADNE_TEST_ASSERT(BigUnsigned(zeros).is_zero());
}

void TestBigUnsigned::test_boundaries()
{
    const auto max=std::numeric_limits<std::uint64_t>::max();
    const BigUnsigned all(max);
    const BigUnsigned one(std::uint64_t{1});
    ARIADNE_TEST_ASSERT(add(all,one).size()==2u);
    ARIADNE_TEST_ASSERT(add(all,one).limb(0)==0u);
    ARIADNE_TEST_ASSERT(add(all,one).limb(1)==1u);
    ARIADNE_TEST_ASSERT(compare(sub(all,all),BigUnsigned{})==0);
    ARIADNE_TEST_ASSERT(compare(shift_right(shift_left(all,257u),257u),all)==0);
    ARIADNE_TEST_ASSERT(compare(shift_right(all,64u),BigUnsigned{})==0);
    ARIADNE_TEST_ASSERT(compare(shift_left(one,256u),BigUnsigned(
        std::array<std::uint64_t,5>{0,0,0,0,1}))==0);
    ARIADNE_TEST_ASSERT(compare(shift_right(shift_left(one,63u),63u),one)==0);

    bool threw = false;
    try {
        (void)sub(one,all);
    } catch (std::underflow_error const&) {
        threw = true;
    }
    ARIADNE_TEST_ASSERT(threw);
    ARIADNE_TEST_ASSERT(compare(add(all,BigUnsigned{}),all)==0);

    // Scan across equal high limbs in heap storage. The first
    // differing limb is near the bottom, or does not exist at all.
    std::array<std::uint64_t,64> equal_high{};
    equal_high.fill(max);
    std::array<std::uint64_t,64> low_different=equal_high;
    low_different[0u]=max-1u;
    const BigUnsigned high(equal_high);
    const BigUnsigned low(low_different);
    ARIADNE_TEST_ASSERT(compare(low,high)<0);
    ARIADNE_TEST_ASSERT(compare(high,low)>0);
    ARIADNE_TEST_ASSERT(compare(high,high)==0);
    ARIADNE_TEST_ASSERT(compare(high,one)>0);
    low_different[63u]=max-1u;
    ARIADNE_TEST_ASSERT(compare(BigUnsigned(low_different),high)<0);
}

void TestBigUnsigned::test_multiplication()
{
    const auto max = std::numeric_limits<std::uint64_t>::max();
    const BigUnsigned one(std::uint64_t{1});
    const std::array<std::uint64_t,4> four_max{max,max,max,max};
    const BigUnsigned four(four_max);

    ARIADNE_TEST_ASSERT(mul_low(four,one,0u).is_zero());
    ARIADNE_TEST_ASSERT(compare(mul_low(four,one,1u),BigUnsigned(max)) == 0);
    ARIADNE_TEST_ASSERT(compare(mul(four,one),four) == 0);
    ARIADNE_TEST_ASSERT(compare(square(one),one) == 0);

    std::mt19937_64 engine(0x1282561024ULL);
    GmpInteger a,b,c,t;

    static_assert(std::numeric_limits<mp_bitcnt_t>::max() >= 64u*80u);
    for (unsigned repetition=0; repetition<6000u; ++repetition) {
        const std::size_t n = engine()%34u;
        const std::size_t m = engine()%34u;
        std::array<std::uint64_t,34> da{}, db{};

        for (std::size_t j=0u; j<n; ++j) da[j] = engine();
        for (std::size_t j=0u; j<m; ++j) db[j] = engine();
        if (repetition%7u == 0u) {
            for (std::size_t j=0u; j<n; ++j) da[j] = max;
        }
        if (repetition%11u == 0u) {
            for (std::size_t j=0u; j<m; ++j) db[j] = max;
        }
        if (repetition%13u == 0u) {
            for (std::size_t j=0u; j<n; ++j) da[j] = 0u;
        }
        if (repetition%17u == 0u) {
            for (std::size_t j=0u; j<m; ++j) db[j] = 0u;
        }

        const BigUnsigned x(std::span(da.data(),n));
        const BigUnsigned y(std::span(db.data(),m));
        mpz_import(a.value,n,-1,8,0,0,da.data());
        mpz_import(b.value,m,-1,8,0,0,db.data());

        mpz_mul(c.value,a.value,b.value);
        check_against_gmp(mul(x,y),c.value);
        ARIADNE_TEST_ASSERT(compare(mul(x,y),mul(y,x)) == 0);

        mpz_mul(t.value,a.value,a.value);
        check_against_gmp(square(x),t.value);
        check_against_gmp(mul(x,x),t.value);

        const std::size_t low_limbs = repetition%5u == 0u ? 0u :
            repetition%5u == 1u ? 1u : repetition%5u == 2u ? 4u :
            repetition%5u == 3u ? 5u : engine()%80u;
        mpz_fdiv_r_2exp(t.value,c.value,
                        static_cast<mp_bitcnt_t>(64u*low_limbs));
        check_against_gmp(mul_low(x,y,low_limbs),t.value);
    }
}

void TestBigUnsigned::test_views_and_scratch()
{
    const auto max = std::numeric_limits<std::uint64_t>::max();
    BigUnsigned::Scratch scratch;
    ARIADNE_TEST_ASSERT(scratch.capacity() >= 16u);

    auto short_view = scratch.acquire(4u);
    for (auto limb : short_view) {
        ARIADNE_TEST_ASSERT(limb == 0u);
    }
    short_view[0] = 11u;
    BigUnsigned x;
    x.assign(short_view);
    ARIADNE_TEST_ASSERT(x.size() == 1u);
    ARIADNE_TEST_ASSERT(x.limb(0) == 11u);

    auto heap_view = scratch.acquire(35u);
    ARIADNE_TEST_ASSERT(scratch.capacity() >= 35u);
    for (std::size_t i=0u; i<heap_view.size(); ++i) {
        heap_view[i] = max;
    }
    x.assign(heap_view);
    ARIADNE_TEST_ASSERT(x.size() == 35u);
    const auto* old_heap = heap_view.data();
    heap_view = scratch.acquire(25u);
    ARIADNE_TEST_ASSERT(heap_view.data() == old_heap);
    for (auto limb : heap_view) {
        ARIADNE_TEST_ASSERT(limb == 0u);
    }

    x.assign(x.limbs());
    ARIADNE_TEST_ASSERT(x.size() == 35u);
    x.assign(x.limbs().subspan(3u,2u));
    ARIADNE_TEST_ASSERT(x.size() == 2u);
    x.assign(x.limbs().first(1u));
    ARIADNE_TEST_ASSERT(x.size() == 1u);

    GmpInteger aa,bb,expected;
    BigUnsigned::Scratch output;
    std::mt19937_64 engine(0xB22B22ULL);
    static_assert(std::numeric_limits<mp_bitcnt_t>::max() >= 65u*64u);
    for (unsigned repetition=0u; repetition<5000u; ++repetition) {
        const std::size_t n = engine()%28u;
        const std::size_t m = engine()%28u;
        std::array<std::uint64_t,28> da{},db{};
        for (std::size_t i=0u; i<n; ++i) da[i] = engine();
        for (std::size_t i=0u; i<m; ++i) db[i] = engine();
        if (repetition%7u == 0u) {
            for (std::size_t i=0u; i<n; ++i) da[i] = max;
        }
        if (repetition%11u == 0u) {
            for (std::size_t i=0u; i<m; ++i) db[i] = max;
        }
        const BigUnsigned lhs(std::span(da.data(),n));
        const BigUnsigned rhs(std::span(db.data(),m));
        mpz_import(aa.value,n,-1,8,0,0,da.data());
        mpz_import(bb.value,m,-1,8,0,0,db.data());

        const std::size_t count = std::max(lhs.size(),rhs.size());
        auto sum_view = output.acquire(count+2u);
        const auto carry = BigUnsigned::add_into(
            lhs.limbs(),rhs.limbs(),sum_view.first(count));
        sum_view[count] = carry;
        ARIADNE_TEST_ASSERT(sum_view[count+1u] == 0u);
        mpz_add(expected.value,aa.value,bb.value);
        check_against_gmp(BigUnsigned(sum_view),expected.value);

        auto difference_view = output.acquire(count+1u);
        const auto borrow = BigUnsigned::sub_into(
            lhs.limbs(),rhs.limbs(),difference_view);
        ARIADNE_TEST_ASSERT(difference_view[count] == 0u);
        if (mpz_cmp(aa.value,bb.value) >= 0) {
            ARIADNE_TEST_ASSERT(borrow == 0u);
            mpz_sub(expected.value,aa.value,bb.value);
            check_against_gmp(BigUnsigned(difference_view),expected.value);
        } else if (count != 0u) {
            ARIADNE_TEST_ASSERT(borrow == 1u);
        }

        const std::size_t low = repetition%5u==0u ? 0u :
            repetition%5u==1u ? 1u : repetition%5u==2u ? 4u :
            repetition%5u==3u ? 15u : engine()%65u;
        auto product_view = output.acquire(low);
        BigUnsigned::mul_low_into(lhs.limbs(),rhs.limbs(),product_view);
        mpz_mul(expected.value,aa.value,bb.value);
        mpz_fdiv_r_2exp(expected.value,expected.value,
                        static_cast<mp_bitcnt_t>(64u*low));
        check_against_gmp(BigUnsigned(product_view),expected.value);
    }
}

void TestBigUnsigned::test_multiply_kernel_paths()
{
    const auto max = std::numeric_limits<std::uint64_t>::max();
    std::mt19937_64 engine(0xB3B3B3ULL);
    BigUnsigned::Scratch scratch;
    BigUnsigned::Scratch karatsuba_workspace;
    GmpInteger a,b,expected,rounded;
    static_assert(std::numeric_limits<mp_bitcnt_t>::max() >= 64u*80u);

    // Multiple writes reuse the same result space; both implementations
    // must overwrite all digits, including untouched high destinations.
    for (unsigned repetition=0u; repetition<6000u; ++repetition) {
        const std::size_t n = engine()%38u;
        const std::size_t m = engine()%38u;
        std::array<std::uint64_t,38> da{},db{};
        for (std::size_t i=0u; i<n; ++i) da[i] = engine();
        for (std::size_t i=0u; i<m; ++i) db[i] = engine();
        if (repetition%7u==0u) std::fill_n(da.begin(),n,max);
        if (repetition%11u==0u) std::fill_n(db.begin(),m,max);
        if (repetition%13u==0u) std::fill_n(da.begin(),n,0u);
        if (repetition%17u==0u) std::fill_n(db.begin(),m,0u);
        const BigUnsigned lhs(std::span(da.data(),n));
        const BigUnsigned rhs(std::span(db.data(),m));
        mpz_import(a.value,n,-1,8,0,0,da.data());
        mpz_import(b.value,m,-1,8,0,0,db.data());
        mpz_mul(expected.value,a.value,b.value);

        const std::size_t length = lhs.size()+rhs.size()+2u;
        auto selected = scratch.acquire(length);
        BigUnsigned::mul_into(lhs.limbs(),rhs.limbs(),selected);
        check_against_gmp(BigUnsigned(selected),expected.value);
        const auto selected_copy = std::vector<std::uint64_t>(
            selected.begin(),selected.end());

        selected = scratch.acquire(length);
        BigUnsigned::mul_into_portable(lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),selected_copy.begin()));
        check_against_gmp(BigUnsigned(selected),expected.value);

        // The experimental Comba kernel must match the B3 portable
        // reference word-for-word, also for oversized full destinations.
        selected = scratch.acquire(length);
        BigUnsigned::mul_into_comba(lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),selected_copy.begin()));
        check_against_gmp(BigUnsigned(selected),expected.value);

        selected = scratch.acquire(length);
        BigUnsigned::mul_into_fused(lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),selected_copy.begin()));
        check_against_gmp(BigUnsigned(selected),expected.value);

        selected = scratch.acquire(length);
        BigUnsigned::mul_into_arm64_carry(lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),selected_copy.begin()));
        check_against_gmp(BigUnsigned(selected),expected.value);

        // The explicit x86-64 carry-chain uses the same GMP oracle,
        // including all-ones, empty, asymmetric and reused output cases.
        selected=scratch.acquire(length);
        BigUnsigned::mul_into_x86_64_carry(
            lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),selected_copy.begin()));
        check_against_gmp(BigUnsigned(selected),expected.value);

        // Explicit B3 schoolbook remains available even if the
        // selected kernel changes on GCC x86-64.
        selected = scratch.acquire(length);
        BigUnsigned::mul_into_schoolbook(lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),selected_copy.begin()));
        check_against_gmp(BigUnsigned(selected),expected.value);

        const auto needed=BigUnsigned::karatsuba_scratch_limbs(
            lhs.size(),rhs.size());
        auto work=karatsuba_workspace.acquire(needed);
        selected=scratch.acquire(length);
        BigUnsigned::mul_into_karatsuba(
            lhs.limbs(),rhs.limbs(),selected,work);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),selected_copy.begin()));
        check_against_gmp(BigUnsigned(selected),expected.value);

        const std::size_t low = repetition%5u==0u ? 0u :
            repetition%5u==1u ? 1u : repetition%5u==2u ? 4u :
            repetition%5u==3u ? 5u : engine()%80u;
        selected = scratch.acquire(low);
        BigUnsigned::mul_low_into(lhs.limbs(),rhs.limbs(),selected);
        const auto low_copy = std::vector<std::uint64_t>(
            selected.begin(),selected.end());
        mpz_fdiv_r_2exp(rounded.value,expected.value,
                        static_cast<mp_bitcnt_t>(64u*low));
        check_against_gmp(BigUnsigned(selected),rounded.value);
        selected = scratch.acquire(low);
        BigUnsigned::mul_low_into_portable(
            lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),low_copy.begin()));

        selected = scratch.acquire(low);
        BigUnsigned::mul_low_into_comba(lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),low_copy.begin()));
        check_against_gmp(BigUnsigned(selected),rounded.value);

        selected=scratch.acquire(low);
        BigUnsigned::mul_low_into_fused(lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),low_copy.begin()));
        check_against_gmp(BigUnsigned(selected),rounded.value);

        selected=scratch.acquire(low);
        BigUnsigned::mul_low_into_arm64_carry(
            lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),low_copy.begin()));
        check_against_gmp(BigUnsigned(selected),rounded.value);

        selected=scratch.acquire(low);
        BigUnsigned::mul_low_into_x86_64_carry(
            lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),low_copy.begin()));
        check_against_gmp(BigUnsigned(selected),rounded.value);

        selected=scratch.acquire(low);
        BigUnsigned::mul_low_into_schoolbook(
            lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),low_copy.begin()));
        check_against_gmp(BigUnsigned(selected),rounded.value);
    }

    const std::array<std::uint64_t,1> one{1u};
    std::array<std::uint64_t,1> too_short{};
    bool selected_threw = false;
    bool portable_threw = false;
    bool comba_threw = false;
    bool fused_threw = false;
    bool schoolbook_threw = false;
    bool carry_threw = false;
    bool x86_carry_threw = false;
    bool karatsuba_threw = false;
    ARIADNE_TEST_ASSERT(
        !BigUnsigned::selected_multiply_uses_fused()
        || BigUnsigned::has_native_multiply());
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
    ARIADNE_TEST_ASSERT(BigUnsigned::has_x86_64_carry_multiply());
#endif
#if defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
    // Native ARM64 uses the specialized carry-chain, not common fused.
    ARIADNE_TEST_ASSERT(BigUnsigned::has_arm64_multiply());
    ARIADNE_TEST_ASSERT(!BigUnsigned::selected_multiply_uses_fused());
#endif
    try {
        BigUnsigned::mul_into(one,one,too_short);
    } catch (std::length_error const&) {
        selected_threw = true;
    }
    try {
        BigUnsigned::mul_into_portable(one,one,too_short);
    } catch (std::length_error const&) {
        portable_threw = true;
    }
    try {
        BigUnsigned::mul_into_comba(one,one,too_short);
    } catch (std::length_error const&) {
        comba_threw = true;
    }
    try {
        BigUnsigned::mul_into_fused(one,one,too_short);
    } catch (std::length_error const&) {
        fused_threw = true;
    }
    try {
        BigUnsigned::mul_into_arm64_carry(one,one,too_short);
    } catch (std::length_error const&) {
        carry_threw = true;
    }
    try {
        BigUnsigned::mul_into_x86_64_carry(one,one,too_short);
    } catch (std::length_error const&) {
        x86_carry_threw = true;
    }
    try {
        BigUnsigned::mul_into_schoolbook(one,one,too_short);
    } catch (std::length_error const&) {
        schoolbook_threw = true;
    }
    try {
        BigUnsigned::mul_into_karatsuba(one,one,too_short,{});
    } catch (std::length_error const&) {
        karatsuba_threw = true;
    }
    ARIADNE_TEST_ASSERT(selected_threw && portable_threw &&
                        comba_threw && fused_threw && carry_threw &&
                        x86_carry_threw && schoolbook_threw &&
                        karatsuba_threw);

    // Missing recursive scratch must fail before altering output.
    const std::array<std::uint64_t,17> long_ones{
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
    };
    std::array<std::uint64_t,34> guarded;
    guarded.fill(max);
    ARIADNE_TEST_ASSERT(
        BigUnsigned::karatsuba_scratch_limbs(17u,17u)>0u);
    ARIADNE_TEST_ASSERT(
        BigUnsigned::karatsuba_scratch_limbs(17u,17u,8u) ==
        BigUnsigned::karatsuba_scratch_limbs(17u,17u));
    ARIADNE_TEST_ASSERT(
        BigUnsigned::karatsuba_scratch_limbs(17u,17u,64u)==0u);
    bool invalid_cutoff_threw=false;
    try {
        (void)BigUnsigned::karatsuba_scratch_limbs(17u,17u,7u);
    } catch (std::invalid_argument const&) {
        invalid_cutoff_threw=true;
    }
    ARIADNE_TEST_ASSERT(invalid_cutoff_threw);
    bool scratch_threw=false;
    try {
        BigUnsigned::mul_into_karatsuba(
            long_ones,long_ones,guarded,{});
    } catch (std::length_error const&) {
        scratch_threw=true;
    }
    ARIADNE_TEST_ASSERT(scratch_threw);
    for (auto value: guarded) {
        ARIADNE_TEST_ASSERT(value==max);
    }

    // B8 needs coverage beyond the B3/B6 38-limb randomized range.
    std::mt19937_64 wider_engine(0xB8B8B8ULL);
    for (unsigned repetition=0u; repetition<700u; ++repetition) {
        const std::size_t n=repetition%7u==0u ? 128u :
                            9u+wider_engine()%95u;
        const std::size_t m=repetition%11u==0u ? 1u :
                            9u+wider_engine()%95u;
        std::array<std::uint64_t,128> da{},db{};
        for (std::size_t i=0u;i<n;++i) da[i]=wider_engine();
        for (std::size_t i=0u;i<m;++i) db[i]=wider_engine();
        if (repetition%13u==0u) std::fill_n(da.begin(),n,max);
        if (repetition%17u==0u) std::fill_n(db.begin(),m,max);
        if (repetition%19u==0u) std::fill_n(da.begin(),n,0u);
        if (repetition%23u==0u) std::fill_n(db.begin(),m,0u);
        const BigUnsigned lhs(std::span(da.data(),n));
        const BigUnsigned rhs(std::span(db.data(),m));
        mpz_import(a.value,n,-1,8,0,0,da.data());
        mpz_import(b.value,m,-1,8,0,0,db.data());
        mpz_mul(expected.value,a.value,b.value);
        for (std::size_t cutoff: {8u,16u,32u,64u}) {
            auto product_view=scratch.acquire(lhs.size()+rhs.size()+2u);
            auto work=karatsuba_workspace.acquire(
                BigUnsigned::karatsuba_scratch_limbs(
                    lhs.size(),rhs.size(),cutoff));
            BigUnsigned::mul_into_karatsuba(
                lhs.limbs(),rhs.limbs(),product_view,work,cutoff);
            check_against_gmp(BigUnsigned(product_view),expected.value);
            for (std::size_t i=lhs.size()+rhs.size();
                    i<product_view.size();++i) {
                ARIADNE_TEST_ASSERT(product_view[i]==0u);
            }
        }
    }

    // Karatsuba output is prefilled with nonzero data to prove that the
    // z0/z2 recursion fully overwrites even unbalanced and odd splits.
    // The scratch guard limbs expose any off-by-one workspace overruns.
    // Reuse both buffers across cutoffs, including a cutoff with no split.
    for (auto [n,m]: {std::pair<std::size_t,std::size_t>{17u,17u},
                       {31u,33u}, {63u,67u}, {79u,80u}, {127u,129u},
                       {255u,257u}, {1u,180u}, {8u,8u}}) {
        std::vector<std::uint64_t> da(n,max),db(m,max);
        for (std::size_t i=0u;i<n;++i) {
            if (i%5u==1u) da[i]=wider_engine();
        }
        for (std::size_t i=0u;i<m;++i) {
            if (i%7u==2u) db[i]=wider_engine();
        }
        mpz_import(a.value,n,-1,8,0,0,da.data());
        mpz_import(b.value,m,-1,8,0,0,db.data());
        mpz_mul(expected.value,a.value,b.value);
        std::vector<std::uint64_t> product(n+m+3u,max);
        for (std::size_t cutoff: {8u,16u,32u,64u,128u}) {
            const auto needed=BigUnsigned::karatsuba_scratch_limbs(n,m,cutoff);
            std::vector<std::uint64_t> guards(needed+2u,max);
            const auto work=BigUnsigned::MutableLimbView(
                guards.data()+1u,needed);
            std::fill(product.begin(),product.end(),max);
            BigUnsigned::mul_into_karatsuba(da,db,product,work,cutoff);
            check_against_gmp(BigUnsigned(product),expected.value);
            for (std::size_t i=n+m;i<product.size();++i) {
                ARIADNE_TEST_ASSERT(product[i]==0u);
            }
            ARIADNE_TEST_ASSERT(guards.front()==max);
            ARIADNE_TEST_ASSERT(guards.back()==max);
        }
    }

}

void TestBigUnsigned::test_division_limb()
{
    const auto max=std::numeric_limits<std::uint64_t>::max();
    // Fixed independent cases that exercise the 64-bit remainder carry.
    const std::array<std::uint64_t,2> all_ones{max,max};
    const std::array<std::uint64_t,2> power_of_base{0u,1u};
    const auto [whole,zero_remainder]=divrem_limb(BigUnsigned(all_ones),max);
    ARIADNE_TEST_ASSERT(compare(whole,BigUnsigned(
        std::array<std::uint64_t,2>{1u,1u}))==0);
    ARIADNE_TEST_ASSERT(zero_remainder==0u);
    const auto [one,remainder]=divrem_limb(BigUnsigned(power_of_base),max);
    ARIADNE_TEST_ASSERT(compare(one,BigUnsigned(std::uint64_t{1}))==0);
    ARIADNE_TEST_ASSERT(remainder==1u);
    // A nonzero high half of the remainder must survive on LLP64 hosts.
    const std::uint64_t wide_remainder=(std::uint64_t{1}<<40u)+17u;
    const auto [zero_quotient,wide_rem]=divrem_limb(
        BigUnsigned(wide_remainder),max);
    ARIADNE_TEST_ASSERT(zero_quotient.is_zero());
    ARIADNE_TEST_ASSERT(wide_rem==wide_remainder);
    ARIADNE_TEST_ASSERT(divrem_limb(BigUnsigned{},1u).first.is_zero());

    GmpInteger n,d,q,r,check,native_remainder;
    std::mt19937_64 engine(0xD1D1D1ULL);
    BigUnsigned::Scratch scratch;
    // Cover zero, all-ones, exact powers of radix, the inline
    // 4/5-limb boundary, and large dividends with varied divisors.
    for (unsigned repetition=0u;repetition<6000u;++repetition) {
        const std::size_t count=repetition%23u==0u ? 0u :
                                repetition%29u==0u ? 64u :
                                repetition%7u+1u;
        std::vector<std::uint64_t> digits(count);
        for (auto& limb:digits) limb=engine();
        if (repetition%13u==0u) {
            std::fill(digits.begin(),digits.end(),max);
        }
        if (repetition%17u==0u) {
            std::fill(digits.begin(),digits.end(),0u);
        }
        if (repetition%19u==0u && !digits.empty()) {
            std::fill(digits.begin(),digits.end(),0u);
            digits.back()=1u;
        }
        const BigUnsigned value(digits);
        const std::uint64_t divisor=repetition%11u==0u ? 1u :
                                    repetition%11u==1u ? max :
                                    repetition%11u==2u ? std::uint64_t{1}<<63u :
                                    repetition%11u==3u ? 2u :
                                    engine()|std::uint64_t{1};
        mpz_import(n.value,digits.size(),-1,sizeof(std::uint64_t),
                   0,0,digits.data());
        mpz_import(d.value,1u,-1,sizeof(std::uint64_t),0,0,&divisor);
        mpz_tdiv_qr(q.value,r.value,n.value,d.value);
        const auto [actual,rem]=divrem_limb(value,divisor);
        check_against_gmp(actual,q.value);
        // mpz_get_ui returns unsigned long (32 bits on Windows x64).
        // Import all 64 bits of the native remainder for an exact oracle.
        mpz_import(native_remainder.value,1u,-1,sizeof(rem),0,0,&rem);
        ARIADNE_TEST_ASSERT(mpz_cmp(native_remainder.value,r.value)==0);

        auto output=scratch.acquire(value.size()+2u);
        std::fill(output.begin(),output.end(),max);
        const auto low=BigUnsigned::divrem_limb_into(
            value.limbs(),divisor,output);
        ARIADNE_TEST_ASSERT(low==rem);
        check_against_gmp(BigUnsigned(output),q.value);
        ARIADNE_TEST_ASSERT(output[value.size()]==0u);
        ARIADNE_TEST_ASSERT(output[value.size()+1u]==0u);
        mpz_mul(check.value,q.value,d.value);
        mpz_add(check.value,check.value,r.value);
        ARIADNE_TEST_ASSERT(mpz_cmp(check.value,n.value)==0);
    }

    // Target 128/64 single-digit division including the maximum
    // possible quotient (dividend high word == divisor-1). GMP is the
    // independent oracle and the remainder comparison preserves all
    // 64 bits on Windows LLP64.
    const std::array<std::uint64_t,8> edge_divisors{
        1u,2u,3u,std::uint64_t{0xffffffffu},
        std::uint64_t{1}<<32u,std::uint64_t{1}<<63u,
        max-1u,max};
    const std::array<std::uint64_t,6> edge_lows{
        0u,1u,std::uint64_t{0xffffffffu},
        std::uint64_t{1}<<63u,max-1u,max};
    for (const std::uint64_t divisor:edge_divisors) {
        const std::array<std::uint64_t,4> highs{
            0u,divisor>1u?1u:0u,divisor/2u,divisor-1u};
        mpz_import(d.value,1u,-1,sizeof(divisor),0,0,&divisor);
        for (const std::uint64_t high:highs) {
            for (const std::uint64_t low:edge_lows) {
                const std::array<std::uint64_t,2> dividend{low,high};
                const BigUnsigned value(dividend);
                mpz_import(n.value,dividend.size(),-1,
                           sizeof(std::uint64_t),0,0,dividend.data());
                mpz_tdiv_qr(q.value,r.value,n.value,d.value);
                const auto [native_q,native_r]=divrem_limb(value,divisor);
                check_against_gmp(native_q,q.value);
                mpz_import(native_remainder.value,1u,-1,
                           sizeof(native_r),0,0,&native_r);
                ARIADNE_TEST_ASSERT(mpz_cmp(native_remainder.value,r.value)==0);
                ARIADNE_TEST_ASSERT(native_r<divisor);
                if (high==divisor-1u && low==max) {
                    ARIADNE_TEST_ASSERT(native_q.limb(0u)==max);
                    ARIADNE_TEST_ASSERT(native_r==divisor-1u);
                }
            }
        }
    }

    const std::array<std::uint64_t,2> numerator{max,max};
    std::array<std::uint64_t,3> output{max,max,max};
    bool zero_threw=false,owning_zero_threw=false,short_threw=false;
    try {
        (void)BigUnsigned::divrem_limb_into(numerator,0u,output);
    } catch (std::domain_error const&) {
        zero_threw=true;
    }
    try {
        (void)divrem_limb(BigUnsigned(numerator),0u);
    } catch (std::domain_error const&) {
        owning_zero_threw=true;
    }
    try {
        (void)BigUnsigned::divrem_limb_into(
            numerator,3u,std::span(output).first(1u));
    } catch (std::length_error const&) {
        short_threw=true;
    }
    ARIADNE_TEST_ASSERT(zero_threw && owning_zero_threw && short_threw);
    for (auto limb:output) {
        ARIADNE_TEST_ASSERT(limb==max);
    }
}


void TestBigUnsigned::test_division()
{
    const auto max=std::numeric_limits<std::uint64_t>::max();
    const BigUnsigned zero;
    const BigUnsigned one(std::uint64_t{1});
    const BigUnsigned radix(std::array<std::uint64_t,2>{0u,1u});
    const BigUnsigned radix_squared(
        std::array<std::uint64_t,3>{0u,0u,1u});
    const auto [q0,r0]=divrem(zero,radix);
    ARIADNE_TEST_ASSERT(q0.is_zero() && r0.is_zero());
    const auto [q1,r1]=divrem(one,radix);
    ARIADNE_TEST_ASSERT(q1.is_zero() && compare(r1,one)==0);
    const auto [q2,r2]=divrem(radix,radix);
    ARIADNE_TEST_ASSERT(compare(q2,one)==0 && r2.is_zero());
    const auto [q3,r3]=divrem(radix_squared,radix);
    ARIADNE_TEST_ASSERT(compare(q3,radix)==0 && r3.is_zero());

    GmpInteger a,b,q,r,check;
    std::mt19937_64 engine(0xB15D1A1DULL);
    for (unsigned repetition=0u;repetition<6000u;++repetition) {
        const std::size_t n=repetition%31u==0u?96u:
                             repetition%7u==0u?32u:engine()%17u;
        const std::size_t m=repetition%37u==0u?1u:
                             repetition%11u==0u?96u:
                             1u+engine()%20u;
        std::vector<std::uint64_t> da(n),db(m);
        for (auto& digit:da) digit=engine();
        for (auto& digit:db) digit=engine();
        if (repetition%13u==0u) std::fill(da.begin(),da.end(),max);
        if (repetition%17u==0u) std::fill(db.begin(),db.end(),max);
        if (repetition%23u==0u) std::fill(da.begin(),da.end(),0u);
        if (repetition%29u==0u) {
            std::fill(db.begin(),db.end(),0u);
            db.back()=1u;
        }
        if (repetition%43u==0u) {
            std::fill(db.begin(),db.end(),0u);
            db.back()=std::uint64_t{1}<<63u;
        }
        if (repetition%47u==0u) db.back()=1u;
        if (repetition%53u==0u && n==m) da=db;
        const BigUnsigned dividend(da),divisor(db);
        mpz_import(a.value,dividend.size(),-1,sizeof(std::uint64_t),
                   0,0,dividend.limbs().data());
        mpz_import(b.value,divisor.size(),-1,sizeof(std::uint64_t),
                   0,0,divisor.limbs().data());
        mpz_tdiv_qr(q.value,r.value,a.value,b.value);
        const auto [actual_q,actual_r]=divrem(dividend,divisor);
        check_against_gmp(actual_q,q.value);
        check_against_gmp(actual_r,r.value);
        // Independent single-output paths must match both GMP and the
        // complete divrem result without using unrequested output objects.
        check_against_gmp(quot(dividend,divisor),q.value);
        check_against_gmp(rem(dividend,divisor),r.value);
        ARIADNE_TEST_ASSERT(compare(quot(dividend,divisor),actual_q)==0);
        ARIADNE_TEST_ASSERT(compare(rem(dividend,divisor),actual_r)==0);
        ARIADNE_TEST_ASSERT(compare(actual_r,divisor)<0);
        ARIADNE_TEST_ASSERT(
            compare(add(mul(actual_q,divisor),actual_r),dividend)==0);

        const auto needed=BigUnsigned::divrem_scratch_limbs(
            dividend.size(),divisor.size());
        std::vector<std::uint64_t> qout(dividend.size()+2u,max);
        std::vector<std::uint64_t> rout(
            std::min(dividend.size(),divisor.size())+2u,max);
        std::vector<std::uint64_t> guards(needed+2u,max);
        const auto scratch=BigUnsigned::MutableLimbView(
            guards.data()+1u,needed);
        BigUnsigned::divrem_into(
            dividend.limbs(),divisor.limbs(),qout,rout,scratch);
        check_against_gmp(BigUnsigned(qout),q.value);
        check_against_gmp(BigUnsigned(rout),r.value);
        for (std::size_t i=dividend.size();i<qout.size();++i)
            ARIADNE_TEST_ASSERT(qout[i]==0u);
        for (std::size_t i=std::min(dividend.size(),divisor.size());
                i<rout.size();++i)
            ARIADNE_TEST_ASSERT(rout[i]==0u);
        ARIADNE_TEST_ASSERT(guards.front()==max);
        ARIADNE_TEST_ASSERT(guards.back()==max);
        mpz_mul(check.value,q.value,b.value);
        mpz_add(check.value,check.value,r.value);
        ARIADNE_TEST_ASSERT(mpz_cmp(check.value,a.value)==0);
    }

    // Force all 64 divisor-normalization shifts. Vary the dividend's
    // top words and saturated bottom words to exercise quotient
    // estimates, downward corrections and full-width remainders.
    for (unsigned shift=0u;shift<64u;++shift) {
        const std::uint64_t top=max>>shift;
        for (unsigned mode=0u;mode<10u;++mode) {
            const std::array<std::uint64_t,2> divisor_words{
                (mode%2u)==0u ? max : std::uint64_t{1},top};
            const std::array<std::uint64_t,4> dividend_words{
                (mode%3u)==0u ? max : std::uint64_t{mode+1u},
                (mode%3u)==1u ? 0u : max,
                (mode%2u)==0u ? max : std::uint64_t{1},
                (mode%3u)==2u ? top-1u : top};
            const BigUnsigned dividend_value(dividend_words);
            const BigUnsigned divisor_value(divisor_words);
            mpz_import(a.value,dividend_words.size(),-1,
                       sizeof(std::uint64_t),0,0,dividend_words.data());
            mpz_import(b.value,divisor_words.size(),-1,
                       sizeof(std::uint64_t),0,0,divisor_words.data());
            mpz_tdiv_qr(q.value,r.value,a.value,b.value);
            const auto [native_q,native_r]=divrem(
                dividend_value,divisor_value);
            check_against_gmp(native_q,q.value);
            check_against_gmp(native_r,r.value);
            check_against_gmp(quot(dividend_value,divisor_value),q.value);
            check_against_gmp(rem(dividend_value,divisor_value),r.value);
            ARIADNE_TEST_ASSERT(compare(native_r,divisor_value)<0);
            ARIADNE_TEST_ASSERT(compare(
                add(mul(native_q,divisor_value),native_r),
                dividend_value)==0);
        }
    }

    // Focused zero, short dividend, exact quotient and single-limb
    // divisor cases; keep both divided-output implementations checked.
    for (const auto& numerator: {
            BigUnsigned{},BigUnsigned(std::uint64_t{1}),
            BigUnsigned(std::array<std::uint64_t,2>{0u,1u}),
            BigUnsigned(std::array<std::uint64_t,4>{max,max,max,max})}) {
        for (const auto& denominator: {
                BigUnsigned(std::uint64_t{1}),BigUnsigned(max),
                BigUnsigned(std::array<std::uint64_t,2>{max,max}),
                BigUnsigned(std::array<std::uint64_t,5>{0u,0u,0u,0u,1u})}) {
            const auto both=divrem(numerator,denominator);
            ARIADNE_TEST_ASSERT(compare(
                quot(numerator,denominator),both.first)==0);
            ARIADNE_TEST_ASSERT(compare(
                rem(numerator,denominator),both.second)==0);
        }
    }

    const std::array<std::uint64_t,3> dividend{max,max,max};
    const std::array<std::uint64_t,2> divisor{max,max};
    std::array<std::uint64_t,4> qout{max,max,max,max};
    std::array<std::uint64_t,3> rout{max,max,max};
    const auto needed=BigUnsigned::divrem_scratch_limbs(3u,2u);
    std::vector<std::uint64_t> scratch(needed,max);
    bool zero_threw=false,owning_zero_threw=false;
    bool short_quotient=false,short_remainder=false,short_scratch=false;
    try {
        BigUnsigned::divrem_into(
            dividend,BigUnsigned{}.limbs(),qout,rout,scratch);
    } catch (std::domain_error const&) { zero_threw=true; }
    try {
        (void)divrem(BigUnsigned(dividend),zero);
    } catch (std::domain_error const&) { owning_zero_threw=true; }
    bool quotient_zero_threw=false,remainder_zero_threw=false;
    try { (void)quot(BigUnsigned(dividend),zero); }
    catch (std::domain_error const&) { quotient_zero_threw=true; }
    try { (void)rem(BigUnsigned(dividend),zero); }
    catch (std::domain_error const&) { remainder_zero_threw=true; }
    ARIADNE_TEST_ASSERT(quotient_zero_threw && remainder_zero_threw);
    try {
        BigUnsigned::divrem_into(
            dividend,divisor,std::span(qout).first(2u),rout,scratch);
    } catch (std::length_error const&) { short_quotient=true; }
    try {
        BigUnsigned::divrem_into(
            dividend,divisor,qout,std::span(rout).first(1u),scratch);
    } catch (std::length_error const&) { short_remainder=true; }
    try {
        BigUnsigned::divrem_into(
            dividend,divisor,qout,rout,std::span(scratch).first(needed-1u));
    } catch (std::length_error const&) { short_scratch=true; }
    ARIADNE_TEST_ASSERT(zero_threw && owning_zero_threw &&
                        short_quotient && short_remainder && short_scratch);
    for (auto digit:qout) ARIADNE_TEST_ASSERT(digit==max);
    for (auto digit:rout) ARIADNE_TEST_ASSERT(digit==max);
    for (auto digit:scratch) ARIADNE_TEST_ASSERT(digit==max);
    bool overflow_threw=false;
    try {
        (void)BigUnsigned::divrem_scratch_limbs(
            std::numeric_limits<std::size_t>::max(),2u);
    } catch (std::length_error const&) { overflow_threw=true; }
    ARIADNE_TEST_ASSERT(overflow_threw);
}


void TestBigUnsigned::test_gcd_decimal()
{
    const auto zero=BigUnsigned{};
    const auto one=BigUnsigned(std::uint64_t{1});
    ARIADNE_TEST_ASSERT(gcd(zero,zero).is_zero());
    ARIADNE_TEST_ASSERT(compare(gcd(zero,one),one)==0);
    ARIADNE_TEST_ASSERT(compare(gcd(one,zero),one)==0);
    ARIADNE_TEST_ASSERT(compare(
        gcd(BigUnsigned(std::uint64_t{48}),BigUnsigned(std::uint64_t{18})),
        BigUnsigned(std::uint64_t{6}))==0);
    ARIADNE_TEST_ASSERT(to_decimal(zero)=="0");
    ARIADNE_TEST_ASSERT(to_decimal(one)=="1");
    ARIADNE_TEST_ASSERT(to_decimal(from_decimal("000"))=="0");
    ARIADNE_TEST_ASSERT(to_decimal(from_decimal("00018446744073709551616"))
                        =="18446744073709551616");
    ARIADNE_TEST_ASSERT(to_decimal(from_decimal(
        "340282366920938463463374607431768211455"))
        =="340282366920938463463374607431768211455");

    for (const auto* bad: {"","-1","+1"," 1","1 ","1.2","0x1","12a"}) {
        bool threw=false;
        try { (void)from_decimal(bad); }
        catch (std::invalid_argument const&) { threw=true; }
        ARIADNE_TEST_ASSERT(threw);
    }

    std::mt19937_64 engine(0xB16B00B5ULL);
    GmpInteger a,b,expected;
    for (unsigned repetition=0u;repetition<2800u;++repetition) {
        const std::size_t n=repetition%113u==0u?96u:engine()%20u;
        const std::size_t m=repetition%71u==0u?32u:engine()%12u;
        std::vector<std::uint64_t> left(n),right(m);
        for (auto& digit:left) digit=engine();
        for (auto& digit:right) digit=engine();
        if (repetition%17u==0u) std::fill(left.begin(),left.end(),0u);
        if (repetition%13u==0u) std::fill(left.begin(),left.end(),
                                           std::numeric_limits<std::uint64_t>::max());
        if (repetition%11u==0u) std::fill(right.begin(),right.end(),0u);
        const BigUnsigned x(left),y(right);
        mpz_import(a.value,x.size(),-1,sizeof(std::uint64_t),
                   0,0,x.limbs().data());
        mpz_import(b.value,y.size(),-1,sizeof(std::uint64_t),
                   0,0,y.limbs().data());
        mpz_gcd(expected.value,a.value,b.value);
        check_against_gmp(gcd(x,y),expected.value);
        const auto decimal=to_decimal(x);
        ARIADNE_TEST_ASSERT(compare(from_decimal(decimal),x)==0);
        ARIADNE_TEST_ASSERT(decimal.front()>='0' && decimal.front()<='9');
        if (!x.is_zero()) ARIADNE_TEST_ASSERT(decimal.front()!='0');
        const auto decimal_gmp=std::string(mpz_sizeinbase(a.value,10u)+2u,'\0');
        std::string oracle=decimal_gmp;
        mpz_get_str(oracle.data(),10,a.value);
        oracle.resize(std::char_traits<char>::length(oracle.c_str()));
        ARIADNE_TEST_ASSERT(decimal==oracle);
    }
    for (unsigned repetition=0u;repetition<1800u;++repetition) {
        const std::size_t length=repetition%149u==0u?1900u:
                                 1u+engine()%370u;
        std::string digits(length,'0');
        for (auto& digit:digits) {
            digit=static_cast<char>('0'+static_cast<int>(engine()%10u));
        }
        const auto parsed=from_decimal(digits);
        ARIADNE_TEST_ASSERT(to_decimal(parsed)==
            to_decimal(from_decimal(to_decimal(parsed))));
        ARIADNE_TEST_ASSERT(mpz_set_str(a.value,digits.c_str(),10)==0);
        check_against_gmp(parsed,a.value);
    }
    // B24: long-input radix-10^18 parsing. Probe each side of the
    // 309/310-digit selector, 18-digit group carry boundaries, and
    // lengths beyond the 2048/4096-bit benchmark workloads. GMP
    // parses the ORIGINAL input, not a native formatted roundtrip.
    const auto check_long_parse=[&](std::string const& input) {
        ARIADNE_TEST_ASSERT(
            mpz_set_str(a.value,input.c_str(),10)==0);
        const BigUnsigned parsed=from_decimal(input);
        check_against_gmp(parsed,a.value);
        std::string oracle(mpz_sizeinbase(a.value,10u)+2u,'\0');
        mpz_get_str(oracle.data(),10,a.value);
        oracle.resize(std::char_traits<char>::length(oracle.c_str()));
        ARIADNE_TEST_ASSERT(to_decimal(parsed)==oracle);
    };
    for (const std::size_t length:
            {308u,309u,310u,311u,312u,616u,617u,1234u,1900u}) {
        check_long_parse(std::string(length,'9'));
        check_long_parse("1"+std::string(length-2u,'0')+"1");
        // The long path is selected by raw input size, even if
        // leading zeros mean that the stored magnitude is small.
        check_long_parse(std::string(length-1u,'0')+"1");
        std::string mixed(length,'0');
        mixed.front()='8';
        mixed[18u]='7';
        mixed.back()='6';
        check_long_parse(mixed);
        std::string alternating(length,'9');
        for (std::size_t i=1u;i<length;i+=2u) {
            alternating[i]='0';
        }
        check_long_parse(alternating);
    }
    // The long-string fast path must never bypass the shared
    // validation of every ASCII digit (including the final byte).
    for (const std::size_t invalid_position:
            {0u,155u,309u,310u}) {
        std::string invalid(311u,'7');
        invalid[invalid_position]='x';
        bool threw=false;
        try { (void)from_decimal(invalid); }
        catch (std::invalid_argument const&) { threw=true; }
        ARIADNE_TEST_ASSERT(threw);
    }

    // Balanced decimal formatting B22: exact radix-10^(9*2^k)
    // boundaries, long zero runs and carry chains, plus signed-free
    // decimal formatting beyond the standard 4096-bit benchmark.
    // The independent GMP oracle uses binary export, not from_decimal.
    const auto check_decimal_boundary=[&](mpz_srcptr reference) {
        const std::size_t words_needed=
            (mpz_sizeinbase(reference,2u)+63u)/64u+1u;
        std::vector<std::uint64_t> words(words_needed,0u);
        std::size_t used=0u;
        mpz_export(words.data(),&used,-1,sizeof(std::uint64_t),
                   0,0,reference);
        const BigUnsigned native(std::span(words.data(),used));
        std::string expected_decimal(
            mpz_sizeinbase(reference,10u)+2u,'\0');
        mpz_get_str(expected_decimal.data(),10,reference);
        expected_decimal.resize(
            std::char_traits<char>::length(expected_decimal.c_str()));

        const std::string decimal=to_decimal(native);
        ARIADNE_TEST_ASSERT(decimal==expected_decimal);
        ARIADNE_TEST_ASSERT(
            compare(from_decimal(decimal),native)==0);
    };
    // B25: independent GMP-oracle regression around the 8/9-limb
    // radix-10^9 / radix-10^18 selector and at recursive leaf sizes.
    // Include all-ones carries, powers of the binary limb radix,
    // long internal zeros and random 64-bit limb patterns.
    std::mt19937_64 leaf_engine(0xB250C0DEULL);
    for (const std::size_t count:
            {7u,8u,9u,10u,15u,16u,17u,32u,64u,128u}) {
        for (unsigned mode=0u;mode<5u;++mode) {
            std::vector<std::uint64_t> words(count,0u);
            if (mode==0u) {
                words.back()=1u;
            } else if (mode==1u) {
                std::fill(words.begin(),words.end(),
                          std::numeric_limits<std::uint64_t>::max());
            } else if (mode==2u) {
                words.front()=1u;
                words.back()=1u;
            } else {
                for (auto& word:words) word=leaf_engine();
                if (mode==4u) {
                    for (std::size_t i=1u;i+1u<count;++i) words[i]=0u;
                }
                words.back()|=std::uint64_t{1}<<63u;
            }
            mpz_import(a.value,words.size(),-1,sizeof(std::uint64_t),
                       0,0,words.data());
            check_decimal_boundary(a.value);
        }
    }

    // B25.1: on supported ARM64 these exercise the 18-digit leaf
    // formatter, while x86-64/portable builds must match GMP via
    // the B19 9-digit fallback. Probe carries and zero-padding at
    // decimal block boundaries near the 8/9/16-limb cutoffs.
    for (const unsigned long groups: {8ul,9ul,10ul,16ul,17ul,18ul}) {
        const unsigned long digits=18ul*groups;
        mpz_ui_pow_ui(a.value,10ul,digits);
        check_decimal_boundary(a.value); // exact power of 10^18
        mpz_sub_ui(a.value,a.value,1ul);
        check_decimal_boundary(a.value); // every low block is 9...
        mpz_add_ui(a.value,a.value,2ul);
        check_decimal_boundary(a.value); // one plus power; zero runs
    }

    for (unsigned k=0u;k<=9u;++k) {
        const unsigned long digits=9ul*(1ul<<k);
        mpz_ui_pow_ui(a.value,10ul,digits);
        check_decimal_boundary(a.value);             // 10^digits
        mpz_sub_ui(a.value,a.value,1ul);
        check_decimal_boundary(a.value);             // 99...99
        mpz_add_ui(a.value,a.value,2ul);
        check_decimal_boundary(a.value);             // 10^digits+1
        mpz_add_ui(a.value,a.value,8ul);
        check_decimal_boundary(a.value);             // 10^digits+9

        // A nonzero suffix after many complete zero-filled radix
        // groups; tests width padding in recursive low partitions.
        const std::string blocks=
            std::string("123456789")+std::string(digits,'0')+
            std::string("000000001");
        ARIADNE_TEST_ASSERT(
            mpz_set_str(a.value,blocks.c_str(),10)==0);
        check_decimal_boundary(a.value);
    }
}

void TestBigUnsigned::test_randomized()
{
    std::mt19937_64 engine(0xBA17B00DULL);
    GmpInteger a,b,c;

    // The tested shift is at most 1029 bits, even when GMP uses 32-bit
    // mp_bitcnt_t (notably with 64-bit size_t on Windows/MSVC).
    static_assert(std::numeric_limits<mp_bitcnt_t>::max() >= 1029u);

    for (unsigned repetition=0; repetition<10000u; ++repetition) {
        const std::size_t n = engine()%20u;
        const std::size_t m = engine()%20u;
        std::array<std::uint64_t,20> na{}, nb{};
        for (std::size_t j=0; j<n; ++j) na[j] = engine();
        for (std::size_t j=0; j<m; ++j) nb[j] = engine();
        if (repetition%7u == 0u) na.fill(std::numeric_limits<std::uint64_t>::max());
        if (repetition%11u == 0u) nb.fill(0u);
        const BigUnsigned x(std::span(na.data(),n));
        const BigUnsigned y(std::span(nb.data(),m));
        mpz_import(a.value,n,-1,8,0,0,na.data());
        mpz_import(b.value,m,-1,8,0,0,nb.data());
        ARIADNE_TEST_ASSERT(compare(x,y)==mpz_cmp(a.value,b.value) ||
            (compare(x,y)<0 && mpz_cmp(a.value,b.value)<0) ||
            (compare(x,y)>0 && mpz_cmp(a.value,b.value)>0));
        mpz_add(c.value,a.value,b.value);
        check_against_gmp(add(x,y),c.value);
        if (mpz_cmp(a.value,b.value) >= 0) {
            mpz_sub(c.value,a.value,b.value);
            check_against_gmp(sub(x,y),c.value);
        }
        const std::size_t bits = engine()%1030u;
        const auto gmp_bits = static_cast<mp_bitcnt_t>(bits);

        mpz_mul_2exp(c.value,a.value,gmp_bits);
        check_against_gmp(shift_left(x,bits),c.value);

        mpz_fdiv_q_2exp(c.value,a.value,gmp_bits);
        check_against_gmp(shift_right(x,bits),c.value);
    }
}

} // namespace

int main()
{
    ARIADNE_TEST_CLASS(BigUnsigned,TestBigUnsigned());
    return ARIADNE_TEST_FAILURES;
}
