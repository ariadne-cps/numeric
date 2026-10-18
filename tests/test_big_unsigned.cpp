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
#include <utility>

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
        ARIADNE_TEST_CALL(test_views_and_scratch());
        ARIADNE_TEST_CALL(test_multiply_kernel_paths());
    }

  private:
    void test_storage();
    void test_boundaries();
    void test_randomized();
    void test_multiplication();
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
    }

    const std::array<std::uint64_t,1> one{1u};
    std::array<std::uint64_t,1> too_short{};
    bool selected_threw = false;
    bool portable_threw = false;
    bool comba_threw = false;
    bool fused_threw = false;
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
    ARIADNE_TEST_ASSERT(
        selected_threw && portable_threw && comba_threw && fused_threw);
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
