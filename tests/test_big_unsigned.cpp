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
        ARIADNE_TEST_CALL(test_subtraction_common_prefix());
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
    void test_subtraction_common_prefix();
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

// B41: the highest differing limb determines both unsigned order
// and the maximum required result width. Exercise common-prefix
// elimination against an independent GMP oracle in both directions,
// including underflow, canonical zero, dense inputs and long borrows.
void TestBigUnsigned::test_subtraction_common_prefix()
{
    const auto max=std::numeric_limits<std::uint64_t>::max();
    std::mt19937_64 engine(0xB41CA11CEULL);
    GmpInteger a,b,expected;
    for (const std::size_t count:
            {0u,1u,2u,7u,8u,16u,31u,32u,63u,64u,65u,128u}) {
        for (unsigned repetition=0u;repetition<64u;++repetition) {
            std::vector<std::uint64_t> left(count),right(count);
            for (auto& word:left) word=engine();
            for (auto& word:right) word=engine();
            switch (repetition%9u) {
            case 0u: // Exact equality: result width zero.
                right=left;
                break;
            case 1u: // Only the lowest word differs.
                right=left;
                if (count!=0u) left[0u]=max;
                break;
            case 2u: // Opposite ordering in the low word.
                right=left;
                if (count!=0u) right[0u]=max;
                break;
            case 3u: // Very long equal high prefix.
                right=left;
                if (count!=0u) {
                    left[0u]=max;
                    right[0u]=0u;
                }
                break;
            case 4u: // Dense, saturated, long carry/borrow chain.
                std::fill(left.begin(),left.end(),max);
                std::fill(right.begin(),right.end(),max);
                if (count!=0u) right.front()=0u;
                break;
            case 5u: // Different canonical widths, left shorter.
                if (count!=0u) left.back()=0u;
                break;
            case 6u: // Different canonical widths, right shorter.
                if (count!=0u) right.back()=0u;
                break;
            case 7u: // Two remote nonzero limbs.
                std::fill(left.begin(),left.end(),0u);
                std::fill(right.begin(),right.end(),0u);
                if (count!=0u) {
                    left.front()=max;
                    left.back()|=std::uint64_t{1}<<63u;
                    right.front()=1u;
                    right.back()|=std::uint64_t{1}<<63u;
                }
                break;
            default: // Independent dense random values.
                break;
            }
            const BigUnsigned lhs(left),rhs(right);
            mpz_import(a.value,lhs.size(),-1,sizeof(std::uint64_t),
                       0,0,lhs.limbs().data());
            mpz_import(b.value,rhs.size(),-1,sizeof(std::uint64_t),
                       0,0,rhs.limbs().data());
            const auto verify=[&](BigUnsigned const& first,
                                  BigUnsigned const& second,
                                  mpz_srcptr x,mpz_srcptr y) {
                mpz_sub(expected.value,x,y);
                if (mpz_sgn(expected.value)<0) {
                    bool threw=false;
                    try { (void)sub(first,second); }
                    catch (std::underflow_error const&) { threw=true; }
                    ARIADNE_TEST_ASSERT(threw);
                } else {
                    check_against_gmp(sub(first,second),expected.value);
                }
            };
            verify(lhs,rhs,a.value,b.value);
            verify(rhs,lhs,b.value,a.value);
        }
    }
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

    // B35: the owning product must agree with the independent GMP
    // oracle across the selective 63/64-limb Karatsuba boundary,
    // odd lengths, near-balanced and strongly unbalanced operands.
    // These tests run on every architecture, including the old
    // unfused/reference multiplication paths.
    std::mt19937_64 karatsuba_selected_engine(0xB35CA7A5ULL);
    for (const auto& [n,m]:
            {std::pair<std::size_t,std::size_t>{63u,63u},
             {63u,64u}, {64u,63u}, {64u,64u}, {64u,65u},
             {65u,64u}, {65u,65u}, {80u,64u}, {96u,64u},
             {96u,96u}, {128u,64u}, {128u,96u},
             {96u,128u}, {128u,128u}}) {
        for (unsigned repetition=0u;repetition<32u;++repetition) {
            std::vector<std::uint64_t> left(n),right(m);
            for (auto& word:left) word=karatsuba_selected_engine();
            for (auto& word:right) word=karatsuba_selected_engine();
            if (repetition%7u==0u) {
                std::fill(left.begin(),left.end(),max);
                std::fill(right.begin(),right.end(),max);
            } else if (repetition%7u==1u) {
                std::fill(left.begin(),left.end(),0u);
                std::fill(right.begin(),right.end(),0u);
                left.front()=13u;
                right.front()=31u;
            } else if (repetition%7u==2u) {
                for (std::size_t i=1u;i+1u<n;++i) left[i]=0u;
                for (std::size_t i=1u;i+1u<m;++i) right[i]=0u;
            } else if (repetition%7u==3u) {
                for (std::size_t i=0u;i<n;i+=2u) left[i]=0u;
                for (std::size_t i=0u;i<m;i+=2u) right[i]=0u;
            }
            left.back()|=std::uint64_t{1}<<63u;
            right.back()|=std::uint64_t{1}<<63u;

            const BigUnsigned lhs(left),rhs(right);
            mpz_import(a.value,n,-1,sizeof(std::uint64_t),
                       0,0,left.data());
            mpz_import(b.value,m,-1,sizeof(std::uint64_t),
                       0,0,right.data());
            mpz_mul(c.value,a.value,b.value);
            check_against_gmp(mul(lhs,rhs),c.value);
            check_against_gmp(mul(rhs,lhs),c.value);
        }
    }

    // B40: regression corpus for the 63/64 active-row dispatch
    // crossover and for 31/32/33-limb Karatsuba leaf products.
    // GMP is the independent oracle in both operand orders, with
    // sparse/dense boundaries and low-product truncations. Other
    // architectures execute the identical cases as a control.
    std::mt19937_64 arm_small_dispatch_engine(0xB4004040ULL);
    for (const auto& [n,m]:
            {std::pair<std::size_t,std::size_t>{7u,64u},
             {8u,64u}, {16u,16u}, {31u,32u}, {32u,31u},
             {32u,32u}, {32u,33u}, {33u,32u}, {33u,33u},
             {63u,63u}, {63u,64u}, {64u,63u},
             {64u,64u}, {64u,65u}, {65u,64u},
             {64u,128u}, {128u,64u}, {128u,128u}}) {
        for (unsigned repetition=0u;repetition<16u;++repetition) {
            std::vector<std::uint64_t> left(n),right(m);
            for (auto& word:left) word=arm_small_dispatch_engine();
            for (auto& word:right) word=arm_small_dispatch_engine();
            const auto zeros_for=[](std::size_t count,
                                    unsigned mode) {
                switch (mode%8u) {
                case 0u: return std::size_t{0u};
                case 1u: return std::size_t{1u};
                case 2u: return count/8u;
                case 3u: return count/8u+1u;
                case 4u: return count/2u;
                case 5u: return count-1u;
                case 6u: return std::size_t{2u};
                default: return count/4u;
                }
            };
            const std::size_t zeros_left=zeros_for(n,repetition);
            const std::size_t zeros_right=zeros_for(m,repetition+3u);
            for (std::size_t i=0u;i<zeros_left && i+1u<n;++i) {
                left[i]=0u;
            }
            for (std::size_t i=0u;i<zeros_right && i+1u<m;++i) {
                right[i]=0u;
            }
            if (repetition==15u) {
                std::fill(left.begin(),left.end(),max);
                std::fill(right.begin(),right.end(),max);
            }
            left.back()|=std::uint64_t{1}<<63u;
            right.back()|=std::uint64_t{1}<<63u;
            const BigUnsigned lhs(left),rhs(right);
            mpz_import(a.value,lhs.size(),-1,sizeof(std::uint64_t),
                       0,0,lhs.limbs().data());
            mpz_import(b.value,rhs.size(),-1,sizeof(std::uint64_t),
                       0,0,rhs.limbs().data());
            mpz_mul(c.value,a.value,b.value);
            check_against_gmp(mul(lhs,rhs),c.value);
            check_against_gmp(mul(rhs,lhs),c.value);
            const std::size_t low=(repetition%2u==0u)
                ? std::size_t{1u} : n+m-1u;
            mpz_fdiv_r_2exp(t.value,c.value,
                            static_cast<mp_bitcnt_t>(64u*low));
            check_against_gmp(mul_low(lhs,rhs,low),t.value);
            check_against_gmp(mul_low(rhs,lhs,low),t.value);
            mpz_mul(t.value,a.value,a.value);
            check_against_gmp(square(lhs),t.value);

            // Exercise explicit cutoff-32 Karatsuba recursion and
            // its small leaves, including sparse and dense operands.
            if (n>=64u && m>=64u) {
                const std::size_t required=
                    BigUnsigned::karatsuba_scratch_limbs(
                        n,m,32u);
                std::vector<std::uint64_t> scratch(required);
                std::vector<std::uint64_t> output(n+m);
                BigUnsigned::mul_into_karatsuba(
                    lhs.limbs(),rhs.limbs(),output,scratch,32u);
                check_against_gmp(BigUnsigned(output),c.value);
            }
        }
    }

    // B39: cover both ARM64 dense/sparse template dispatch thresholds
    // and the owning balanced-Karatsuba 63/64-limb boundary. Test data
    // are independent from the benchmark pool; expected values come
    // solely from GMP. Also run on non-ARM64 architectures as controls.
    std::mt19937_64 adaptive_engine(0xB39DA7AULL);
    for (const auto& [n,m]:
            {std::pair<std::size_t,std::size_t>{7u,8u},
             {8u,8u}, {15u,16u}, {16u,16u},
             {32u,32u}, {63u,63u}, {63u,64u},
             {64u,63u}, {64u,64u}, {64u,65u},
             {65u,64u}, {64u,85u}, {64u,86u},
             {85u,64u}, {86u,64u}, {96u,64u},
             {128u,64u}, {64u,128u}, {128u,128u}}) {
        for (unsigned repetition=0u;repetition<20u;++repetition) {
            std::vector<std::uint64_t> left(n),right(m);
            for (auto& word:left) word=adaptive_engine();
            for (auto& word:right) word=adaptive_engine();
            const auto zero_rows=[&](std::size_t count,
                                     unsigned variant) {
                switch (variant%7u) {
                case 0u: return std::size_t{0u}; // dense
                case 1u: return std::size_t{1u};
                case 2u: return count/8u;       // density boundary
                case 3u: return count/8u+1u;
                case 4u: return count/2u;
                case 5u: return count-2u;       // strongly sparse
                default: return count/4u;
                }
            };
            for (std::size_t i=0u;
                    i<zero_rows(n,repetition) && i+1u<n;++i) {
                left[i]=0u;
            }
            const unsigned right_mode=
                repetition%3u==0u ? repetition : repetition+3u;
            for (std::size_t i=0u;
                    i<zero_rows(m,right_mode) && i+1u<m;++i) {
                right[i]=0u;
            }
            left.back()|=std::uint64_t{1}<<63u;
            right.back()|=std::uint64_t{1}<<63u;
            const BigUnsigned lhs(left),rhs(right);
            mpz_import(a.value,lhs.size(),-1,sizeof(std::uint64_t),
                       0,0,lhs.limbs().data());
            mpz_import(b.value,rhs.size(),-1,sizeof(std::uint64_t),
                       0,0,rhs.limbs().data());
            mpz_mul(c.value,a.value,b.value);
            check_against_gmp(mul(lhs,rhs),c.value);
            check_against_gmp(mul(rhs,lhs),c.value);
            mpz_mul(t.value,a.value,a.value);
            check_against_gmp(square(lhs),t.value);
            const std::size_t low=(repetition%2u==0u) ?
                (n+m)/2u : (n+m-1u);
            mpz_set(t.value,c.value);
            mpz_fdiv_r_2exp(
                t.value,t.value,static_cast<mp_bitcnt_t>(64u*low));
            check_against_gmp(mul_low(lhs,rhs,low),t.value);
        }
    }

    // B54: GCC x86 B38/B40 symmetric-square sparse/dense dispatch
    // and B39's owning Karatsuba density gate. Probe the exact
    // 63/64 size boundary, the 1/8 zero cutoff, sparse powers,
    // alternating limbs, saturated carries and random dense values.
    // All expected full values originate independently from GMP.
    std::mt19937_64 b54_engine(0xB54A6C4ULL);
    for (const auto& [n,m]:
            {std::pair<std::size_t,std::size_t>{32u,32u},
             {33u,33u}, {63u,63u}, {63u,64u},
             {64u,63u}, {64u,64u}, {64u,65u},
             {65u,64u}, {65u,65u}, {64u,85u},
             {64u,86u}, {85u,64u}, {86u,64u},
             {96u,96u}, {127u,128u}, {128u,127u},
             {128u,128u}}) {
        for (unsigned mode=0u;mode<32u;++mode) {
            std::vector<std::uint64_t> left(n),right(m);
            for (auto& word:left) word=b54_engine();
            for (auto& word:right) word=b54_engine();
            const auto apply_pattern=[&](std::vector<std::uint64_t>& words,
                                         unsigned pattern) {
                const std::size_t count=words.size();
                switch (pattern%8u) {
                case 0u: break; // dense random
                case 1u:
                    std::fill(words.begin(),words.end(),max);
                    break;
                case 2u:
                    std::fill(words.begin(),words.end(),0u);
                    words.front()=13u;
                    break;
                case 3u:
                    std::fill_n(words.begin(),count/8u,0u);
                    break;
                case 4u:
                    std::fill_n(words.begin(),count/8u+1u,0u);
                    break;
                case 5u:
                    for (std::size_t i=0u;i<count;i+=2u)
                        words[i]=0u;
                    break;
                case 6u:
                    for (std::size_t i=1u;i+1u<count;++i)
                        words[i]=0u;
                    break;
                default:
                    for (std::size_t i=0u;i<count;i+=8u)
                        words[i]=0u;
                    break;
                }
                words.back()|=std::uint64_t{1}<<63u;
            };
            apply_pattern(left,mode);
            apply_pattern(right,mode+3u);
            const BigUnsigned lhs(left),rhs(right);
            mpz_import(a.value,n,-1,sizeof(std::uint64_t),
                       0,0,left.data());
            mpz_import(b.value,m,-1,sizeof(std::uint64_t),
                       0,0,right.data());
            mpz_mul(c.value,a.value,b.value);
            check_against_gmp(mul(lhs,rhs),c.value);
            check_against_gmp(mul(rhs,lhs),c.value);

            std::vector<std::uint64_t> fused(n+m);
            BigUnsigned::mul_into_fused(
                lhs.limbs(),rhs.limbs(),fused);
            check_against_gmp(BigUnsigned(fused),c.value);
        }
    }

    for (const std::size_t n:
            {7u,8u,9u,16u,31u,32u,33u,
             63u,64u,65u,66u,96u,128u}) {
        for (unsigned mode=0u;mode<64u;++mode) {
            std::vector<std::uint64_t> words(n);
            for (auto& word:words) word=b54_engine();
            if (mode%8u==0u) {
                std::fill(words.begin(),words.end(),max);
            } else if (mode%8u==1u) {
                std::fill(words.begin(),words.end(),0u);
                words.front()=13u;
            } else if (mode%8u==2u) {
                std::fill_n(words.begin(),n/8u,0u);
            } else if (mode%8u==3u) {
                std::fill_n(words.begin(),n/8u+1u,0u);
            } else if (mode%8u==4u) {
                for (std::size_t i=0u;i<n;i+=2u)
                    words[i]=0u;
            } else if (mode%8u==5u) {
                for (std::size_t i=1u;i+1u<n;++i)
                    words[i]=0u;
            } else if (mode%8u==6u) {
                for (std::size_t i=0u;i<n;i+=8u)
                    words[i]=0u;
            }
            words.back()|=std::uint64_t{1}<<63u;
            const BigUnsigned input(words);
            mpz_import(a.value,n,-1,sizeof(std::uint64_t),
                       0,0,words.data());
            mpz_mul(t.value,a.value,a.value);
            check_against_gmp(square(input),t.value);
            check_against_gmp(mul(input,input),t.value);
        }
    }

    // B37: stress the 7/8-limb ARM64 symmetric-square crossover,
    // long carry chains and diagonal/off-diagonal propagation.
    // Also exercise the same public test on x86-64 and fallback builds.
    std::mt19937_64 arm_square_engine(0xB37A11ULL);
    for (const std::size_t count:
            {1u,4u,7u,8u,9u,15u,16u,17u,32u,33u,64u,65u,128u}) {
        for (unsigned repetition=0u;repetition<96u;++repetition) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=arm_square_engine();
            switch (repetition%7u) {
            case 0u:
                std::fill(words.begin(),words.end(),max);
                break;
            case 1u:
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
                break;
            case 2u:
                for (std::size_t i=1u;i+1u<count;++i) words[i]=0u;
                break;
            case 3u:
                for (std::size_t i=0u;i<count;i+=2u) words[i]=0u;
                break;
            default:
                break;
            }
            const BigUnsigned input(words);
            mpz_import(a.value,input.size(),-1,sizeof(std::uint64_t),
                       0,0,input.limbs().data());
            mpz_mul(t.value,a.value,a.value);
            check_against_gmp(square(input),t.value);
            check_against_gmp(mul(input,input),t.value);
        }
    }

    // B38: cover sparse zero-row elimination in full and truncated
    // ARM64 products, including both asymmetric operand orders. The
    // GMP oracle remains independent of all native multiplication
    // kernels, and these cases also execute on other architectures.
    std::mt19937_64 arm_sparse_engine(0xB38C0FFEEULL);
    for (const auto& [n,m]:
            {std::pair<std::size_t,std::size_t>{1u,3u},
             {4u,7u}, {7u,8u}, {8u,8u}, {9u,16u},
             {16u,32u}, {32u,64u}, {64u,32u},
             {64u,128u}, {128u,64u}, {64u,129u},
             {129u,64u}, {96u,128u}, {128u,96u}}) {
        for (unsigned repetition=0u;repetition<24u;++repetition) {
            std::vector<std::uint64_t> left(n),right(m);
            for (auto& word:left) word=arm_sparse_engine();
            for (auto& word:right) word=arm_sparse_engine();
            switch (repetition%8u) {
            case 0u: // Sparse first operand: one high nonzero limb.
                std::fill(left.begin(),left.end(),0u);
                left.back()=std::uint64_t{1}<<63u;
                break;
            case 1u: // Sparse first operand with two remote limbs.
                std::fill(left.begin(),left.end(),0u);
                left.front()=13u;
                left.back()=std::uint64_t{1}<<63u;
                break;
            case 2u: // Sparse second operand (reversed orientation).
                std::fill(right.begin(),right.end(),0u);
                right.front()=17u;
                right.back()=std::uint64_t{1}<<63u;
                break;
            case 3u: // Sparse both, nonzero carry-prone endpoints.
                std::fill(left.begin(),left.end(),0u);
                std::fill(right.begin(),right.end(),0u);
                left.front()=max;
                right.front()=max;
                left.back()=max;
                right.back()=max;
                break;
            case 4u: // Every second row contributes no product.
                for (std::size_t i=0u;i<n;i+=2u) left[i]=0u;
                break;
            case 5u: // Dense saturated products and long carries.
                std::fill(left.begin(),left.end(),max);
                std::fill(right.begin(),right.end(),max);
                break;
            default: // Dense random and mixed patterns.
                break;
            }
            left.back()|=std::uint64_t{1}<<63u;
            right.back()|=std::uint64_t{1}<<63u;
            const BigUnsigned lhs(left),rhs(right);
            mpz_import(a.value,lhs.size(),-1,sizeof(std::uint64_t),
                       0,0,lhs.limbs().data());
            mpz_import(b.value,rhs.size(),-1,sizeof(std::uint64_t),
                       0,0,rhs.limbs().data());
            mpz_mul(c.value,a.value,b.value);
            check_against_gmp(mul(lhs,rhs),c.value);
            check_against_gmp(mul(rhs,lhs),c.value);
            const std::array<std::size_t,6> low_lengths{
                0u,1u,n/2u,n,n+m-1u,n+m+2u};
            for (const std::size_t low_limbs:low_lengths) {
                mpz_set(t.value,c.value);
                mpz_fdiv_r_2exp(t.value,t.value,
                    static_cast<mp_bitcnt_t>(64u*low_limbs));
                check_against_gmp(mul_low(lhs,rhs,low_limbs),t.value);
                check_against_gmp(mul_low(rhs,lhs,low_limbs),t.value);
            }
        }
    }

    // B32: validate the symmetric squaring threshold and large carry
    // chains against GMP; verify mul(x,x) independently remains exact.
    // Include dense/sparse, saturated and exact high-limb boundaries.
    std::mt19937_64 square_engine(0xB3205A7EULL);
    for (const std::size_t count:
            {1u,2u,7u,8u,9u,12u,16u,32u,64u,96u,128u}) {
        for (unsigned repetition=0u;repetition<192u;++repetition) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=square_engine();
            if (repetition%13u==0u) {
                std::fill(words.begin(),words.end(),max);
            } else if (repetition%17u==0u) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
            } else if (repetition%19u==0u) {
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()=max;
            } else if (repetition%23u==0u) {
                std::fill(words.begin(),words.end(),0u);
            }
            const BigUnsigned value(words);
            mpz_import(a.value,value.size(),-1,sizeof(std::uint64_t),
                       0,0,value.limbs().data());
            mpz_mul(t.value,a.value,a.value);
            check_against_gmp(square(value),t.value);
            check_against_gmp(mul(value,value),t.value);
        }
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

    // B52: data-driven 32x32 GCC x86-64 full-product dispatch
    // to the baseline MULQ/ADCQ kernel. GMP is the independent
    // product oracle. Exercise 31/32/33-limb boundaries, including
    // noncanonical zero data, sparse rows, all-ones, extra output
    // capacity and low products. No workload fixtures are changed.
    std::mt19937_64 b52_engine(0xB5203232ULL);
    for (const std::size_t n:{31u,32u,33u}) {
        for (const std::size_t m:{31u,32u,33u}) {
            for (unsigned mode=0u;mode<32u;++mode) {
                std::vector<std::uint64_t> left(n),right(m);
                for (auto& digit:left) digit=b52_engine();
                for (auto& digit:right) digit=b52_engine();
                if (mode%8u==0u) {
                    std::fill(left.begin(),left.end(),max);
                    std::fill(right.begin(),right.end(),max);
                } else if (mode%8u==1u) {
                    std::fill(left.begin(),left.end(),0u);
                    std::fill(right.begin(),right.end(),0u);
                    left.back()=1u;
                    right.back()=1u;
                } else if (mode%8u==2u) {
                    for (std::size_t i=0u;i<n-1u;i+=2u)
                        left[i]=0u;
                    for (std::size_t i=0u;i<m-1u;i+=2u)
                        right[i]=0u;
                } else if (mode%8u==3u) {
                    std::fill(left.begin(),left.end(),0u);
                    std::fill(right.begin(),right.end(),0u);
                    left.front()=max;
                    right.front()=max;
                } else if (mode%8u==4u) {
                    for (std::size_t i=1u;i+1u<n;++i)
                        left[i]=0u;
                    for (std::size_t i=1u;i+1u<m;++i)
                        right[i]=0u;
                }
                const BigUnsigned lhs(left),rhs(right);
                mpz_import(a.value,lhs.size(),-1,
                           sizeof(std::uint64_t),0,0,lhs.limbs().data());
                mpz_import(b.value,rhs.size(),-1,
                           sizeof(std::uint64_t),0,0,rhs.limbs().data());
                mpz_mul(expected.value,a.value,b.value);
                const std::size_t full=lhs.size()+rhs.size();
                std::vector<std::uint64_t> selected(full,max);
                BigUnsigned::mul_into(
                    lhs.limbs(),rhs.limbs(),selected);
                check_against_gmp(BigUnsigned(selected),expected.value);
                check_against_gmp(mul(lhs,rhs),expected.value);

                std::vector<std::uint64_t> fused(full,max);
                BigUnsigned::mul_into_fused(
                    lhs.limbs(),rhs.limbs(),fused);
                ARIADNE_TEST_ASSERT(selected==fused);

                std::vector<std::uint64_t> oversized(full+2u,max);
                BigUnsigned::mul_into(
                    lhs.limbs(),rhs.limbs(),oversized);
                check_against_gmp(BigUnsigned(oversized),expected.value);
                ARIADNE_TEST_ASSERT(oversized[full]==0u);
                ARIADNE_TEST_ASSERT(oversized[full+1u]==0u);

                const std::size_t low=full==0u?0u:full-1u;
                std::vector<std::uint64_t> truncated(low,max);
                BigUnsigned::mul_low_into(
                    lhs.limbs(),rhs.limbs(),truncated);
                mpz_fdiv_r_2exp(rounded.value,expected.value,
                                 static_cast<mp_bitcnt_t>(64u*low));
                check_against_gmp(BigUnsigned(truncated),rounded.value);
            }
        }
    }

    // B53: independent GMP oracle for x86 GCC fused and GNU/Clang
    // assembly zero-row skipping. Include 31/32/33 B52 boundaries,
    // 63/64/65 adaptive-density thresholds, unbalanced rows and
    // truncated/oversized outputs. Raw views include noncanonical
    // zero-high limbs, so skipping must not depend on normalization.
    std::mt19937_64 b53_engine(0xB530B040ULL);
    for (const auto& [n,m]:
            {std::pair<std::size_t,std::size_t>{1u,1u},
             {7u,64u}, {16u,16u}, {31u,32u}, {32u,31u},
             {32u,32u}, {32u,33u}, {33u,32u}, {63u,64u},
             {64u,63u}, {64u,64u}, {64u,65u}, {65u,64u},
             {64u,128u}, {128u,64u}, {128u,128u}}) {
        for (unsigned mode=0u;mode<24u;++mode) {
            std::vector<std::uint64_t> left(n),right(m);
            for (auto& word:left) word=b53_engine();
            for (auto& word:right) word=b53_engine();
            switch (mode%8u) {
            case 0u:
                std::fill(left.begin(),left.end(),max);
                std::fill(right.begin(),right.end(),max);
                break;
            case 1u:
                std::fill(left.begin(),left.end(),0u);
                std::fill(right.begin(),right.end(),0u);
                break;
            case 2u:
                for (std::size_t i=0u;i<n;i+=2u) left[i]=0u;
                for (std::size_t i=0u;i<m;i+=2u) right[i]=0u;
                break;
            case 3u:
                std::fill(left.begin()+1u,left.end(),0u);
                std::fill(right.begin()+1u,right.end(),0u);
                break;
            case 4u:
                std::fill_n(left.begin(),n/2u,0u);
                std::fill_n(right.begin(),m/2u,0u);
                break;
            case 5u:
                std::fill(left.begin(),left.end(),0u);
                left.back()=max;
                std::fill(right.begin(),right.end(),0u);
                right.back()=max;
                break;
            case 6u:
                for (std::size_t i=0u;i<n;i+=8u) left[i]=0u;
                for (std::size_t i=0u;i<m;i+=8u) right[i]=0u;
                break;
            default: break;
            }
            mpz_import(a.value,n,-1,sizeof(std::uint64_t),
                       0,0,left.data());
            mpz_import(b.value,m,-1,sizeof(std::uint64_t),
                       0,0,right.data());
            mpz_mul(expected.value,a.value,b.value);

            for (const std::size_t length:
                    {std::size_t{0u},std::size_t{1u},
                     (n+m)/2u,n+m-1u,n+m,n+m+2u}) {
                std::vector<std::uint64_t> selected(length,max);
                std::vector<std::uint64_t> fused(length,max);
                std::vector<std::uint64_t> carry(length,max);
                BigUnsigned::mul_low_into(left,right,selected);
                BigUnsigned::mul_low_into_fused(left,right,fused);
                BigUnsigned::mul_low_into_x86_64_carry(
                    left,right,carry);
                ARIADNE_TEST_ASSERT(selected==fused);
                ARIADNE_TEST_ASSERT(selected==carry);
                mpz_fdiv_r_2exp(
                    rounded.value,expected.value,
                    static_cast<mp_bitcnt_t>(64u*length));
                check_against_gmp(BigUnsigned(selected),rounded.value);
            }
        }
    }

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

    // B26: stress exact 128/64 DIVQ preconditions. Each numerator
    // has high < divisor (the hardware quotient must fit in a limb).
    // Exercise saturated quotient, high=0, high=divisor-1,
    // random carry, and divisors with and without their top bit set.
    // All expected results come from GMP, never the native Wide path.
    std::mt19937_64 wide_engine(0xB260D1B1ULL);
    for (unsigned repetition=0u;repetition<2048u;++repetition) {
        const std::uint64_t divisor=repetition%11u==0u ? max :
            repetition%11u==1u ? std::uint64_t{1}<<63u :
            repetition%11u==2u ? 3u :
            repetition%11u==3u ? std::uint64_t{1}<<32u :
            wide_engine()|std::uint64_t{1};
        const std::uint64_t high=repetition%5u==0u ?
            divisor-1u : repetition%5u==1u ? 0u :
            wide_engine()%divisor;
        const std::uint64_t low=repetition%7u==0u ?
            max : repetition%7u==1u ? 0u : wide_engine();
        const std::array<std::uint64_t,2> dividend{low,high};
        const BigUnsigned value(dividend);
        mpz_import(n.value,dividend.size(),-1,sizeof(std::uint64_t),
                   0,0,dividend.data());
        mpz_import(d.value,1u,-1,sizeof(divisor),0,0,&divisor);
        mpz_tdiv_qr(q.value,r.value,n.value,d.value);
        const auto [native_q,native_r]=divrem_limb(value,divisor);
        check_against_gmp(native_q,q.value);
        mpz_import(native_remainder.value,1u,-1,sizeof(native_r),
                   0,0,&native_r);
        ARIADNE_TEST_ASSERT(mpz_cmp(native_remainder.value,r.value)==0);
        ARIADNE_TEST_ASSERT(native_r<divisor);
    }

    // B48: the ARM64 32-bit-divisor fast path must return the exact
    // full quotient and 64-bit remainder, independently of the
    // production divide_wide implementation. Include saturated and
    // sparse operands across inline/heap storage boundaries and
    // the 2^32 divisor selector boundary.
    std::mt19937_64 narrow_engine(0xB48032D0ULL);
    for (const std::size_t count:
            {1u,2u,3u,4u,5u,8u,16u,32u,64u,96u}) {
        for (const std::uint64_t divisor:
                {1ULL,2ULL,3ULL,7ULL,9ULL,10ULL,31ULL,63ULL,
                 255ULL,65535ULL,2147483648ULL,4294967295ULL,
                 4294967296ULL,18446744073709551615ULL}) {
            for (unsigned mode=0u;mode<4u;++mode) {
                std::vector<std::uint64_t> words(count);
                for (auto& limb:words) limb=narrow_engine();
                if (mode==0u) {
                    std::fill(words.begin(),words.end(),max);
                } else if (mode==1u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=1u;
                } else if (mode==2u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.front()=max;
                    words.back()|=1u;
                }
                const BigUnsigned value(words);
                mpz_import(n.value,value.size(),-1,sizeof(std::uint64_t),
                           0,0,value.limbs().data());
                mpz_import(d.value,1u,-1,sizeof(divisor),0,0,&divisor);
                mpz_tdiv_qr(q.value,r.value,n.value,d.value);
                const auto [got_quotient,got_remainder]=
                    divrem_limb(value,divisor);
                check_against_gmp(got_quotient,q.value);
                mpz_import(native_remainder.value,1u,-1,
                           sizeof(got_remainder),0,0,&got_remainder);
                ARIADNE_TEST_ASSERT(
                    mpz_cmp(native_remainder.value,r.value)==0);
                auto out=scratch.acquire(value.size()+2u);
                std::fill(out.begin(),out.end(),max);
                const auto low=BigUnsigned::divrem_limb_into(
                    value.limbs(),divisor,out);
                ARIADNE_TEST_ASSERT(low==got_remainder);
                check_against_gmp(BigUnsigned(out),q.value);
                ARIADNE_TEST_ASSERT(out[value.size()]==0u);
                ARIADNE_TEST_ASSERT(out[value.size()+1u]==0u);
            }
        }
    }

    // B49: exercise the pre-inverted 128/32 quotient path on ARM64
    // across 7/8/9 limb dispatch boundaries and longer dividends.
    // GMP is the independent oracle for full quotient, remainder and
    // all output limbs. Test divisor 1, 2^32-1 and 2^32 fallback too.
    std::mt19937_64 reciprocal_engine(0xB490C0DEULL);
    for (const std::size_t count:
            {7u,8u,9u,16u,32u,64u,96u,128u}) {
        for (const std::uint64_t divisor:
                {1ULL,2ULL,3ULL,5ULL,7ULL,9ULL,31ULL,
                 63ULL,255ULL,65535ULL,2147483648ULL,
                 4294967295ULL,4294967296ULL,
                 18446744073709551615ULL}) {
            for (unsigned mode=0u;mode<7u;++mode) {
                std::vector<std::uint64_t> words(count);
                for (auto& word:words) word=reciprocal_engine();
                if (mode==0u) {
                    std::fill(words.begin(),words.end(),max);
                } else if (mode==1u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=1u;
                } else if (mode==2u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.front()=max;
                    words.back()=1u;
                } else if (mode==3u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=std::uint64_t{1}<<63u;
                } else if (mode==4u) {
                    std::fill(words.begin(),words.end(),0u);
                } else if (mode==5u) {
                    for (std::size_t i=0u;i<count;i+=2u) {
                        words[i]=0u;
                    }
                }
                const BigUnsigned value(words);
                mpz_import(n.value,value.size(),-1,sizeof(std::uint64_t),
                           0,0,value.limbs().data());
                mpz_import(d.value,1u,-1,sizeof(divisor),0,0,&divisor);
                mpz_tdiv_qr(q.value,r.value,n.value,d.value);
                const auto [actual_q,actual_r]=
                    divrem_limb(value,divisor);
                check_against_gmp(actual_q,q.value);
                mpz_import(native_remainder.value,1u,-1,
                           sizeof(actual_r),0,0,&actual_r);
                ARIADNE_TEST_ASSERT(
                    mpz_cmp(native_remainder.value,r.value)==0);

                auto output=scratch.acquire(value.size()+2u);
                std::fill(output.begin(),output.end(),max);
                const auto rem_out=BigUnsigned::divrem_limb_into(
                    value.limbs(),divisor,output);
                ARIADNE_TEST_ASSERT(rem_out==actual_r);
                check_against_gmp(BigUnsigned(output),q.value);
                ARIADNE_TEST_ASSERT(output[value.size()]==0u);
                ARIADNE_TEST_ASSERT(output[value.size()+1u]==0u);
            }
        }
    }

    // B57: normalized pi1 selection for >=64 raw quotient limbs
    // on GCC x86-64, and universal one-limb powers-of-two.
    // GMP is an independent oracle for saturated, random, sparse,
    // denormalized and 63/64/65-limb boundary views.
    std::mt19937_64 b57_engine(0xB570641ULL);
    for (const std::size_t count:
            {0u,1u,7u,8u,16u,31u,32u,63u,64u,65u,96u,128u}) {
        for (const std::uint64_t divisor:
                {1ULL,2ULL,3ULL,5ULL,7ULL,35ULL,255ULL,
                 65535ULL,4294967295ULL,4294967296ULL,
                 4294967297ULL,9223372036854775807ULL,
                 9223372036854775808ULL,9223372036854775809ULL,
                 18446744073709551615ULL}) {
            for (unsigned mode=0u;mode<8u;++mode) {
                std::vector<std::uint64_t> words(count);
                for (auto& word:words) word=b57_engine();
                if (mode==0u) {
                    std::fill(words.begin(),words.end(),max);
                } else if (mode==1u) {
                    std::fill(words.begin(),words.end(),0u);
                } else if (mode==2u && !words.empty()) {
                    std::fill(words.begin(),words.end(),0u);
                    words.front()=1u;
                } else if (mode==3u && !words.empty()) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=1u;
                } else if (mode==4u && !words.empty()) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=max;
                } else if (mode==5u) {
                    for (std::size_t i=0u;i<count;i+=2u)
                        words[i]=0u;
                } else if (mode==6u) {
                    for (std::size_t i=0u;i<count;i+=3u)
                        words[i]=max;
                }
                mpz_import(n.value,count,-1,sizeof(std::uint64_t),
                           0,0,words.data());
                mpz_import(d.value,1u,-1,sizeof(std::uint64_t),
                           0,0,&divisor);
                mpz_tdiv_qr(q.value,r.value,n.value,d.value);
                auto out=scratch.acquire(count+2u);
                std::fill(out.begin(),out.end(),max);
                const auto limb_remainder=BigUnsigned::divrem_limb_into(
                    words,divisor,out);
                check_against_gmp(BigUnsigned(out),q.value);
                mpz_import(native_remainder.value,1u,-1,
                           sizeof(limb_remainder),0,0,&limb_remainder);
                ARIADNE_TEST_ASSERT(
                    mpz_cmp(native_remainder.value,r.value)==0);
                ARIADNE_TEST_ASSERT(out[count]==0u);
                ARIADNE_TEST_ASSERT(out[count+1u]==0u);
                const BigUnsigned value(words);
                const auto [owned_q,owned_r]=
                    divrem_limb(value,divisor);
                check_against_gmp(owned_q,q.value);
                ARIADNE_TEST_ASSERT(owned_r==limb_remainder);
            }
        }
    }

    // B61: exact single-limb remainder-only pi1 with a GCC x86-64
    // dispatch boundary at 7/8/9 limbs; compare against independent
    // GMP for zero, saturated, sparse and random dividends.
    // Also guard the existing 63/64/65-limb quotient pi1 boundary.
    std::mt19937_64 b61_engine(0xB610D1A6ULL);
    for (const std::size_t count:
            {0u,1u,7u,8u,9u,15u,16u,31u,32u,63u,64u,65u,128u}) {
        for (const std::uint64_t divisor:
                {1ULL,2ULL,3ULL,5ULL,7ULL,35ULL,255ULL,
                 65535ULL,4294967295ULL,4294967296ULL,
                 4294967297ULL,9223372036854775807ULL,
                 9223372036854775808ULL,9223372036854775809ULL,
                 18446744073709551615ULL}) {
            for (unsigned mode=0u;mode<8u;++mode) {
                std::vector<std::uint64_t> words(count);
                for (auto& word:words) word=b61_engine();
                if (mode==0u) {
                    std::fill(words.begin(),words.end(),max);
                } else if (mode==1u) {
                    std::fill(words.begin(),words.end(),0u);
                } else if (mode==2u && !words.empty()) {
                    std::fill(words.begin(),words.end(),0u);
                    words.front()=1u;
                } else if (mode==3u && !words.empty()) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=1u;
                } else if (mode==4u && !words.empty()) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=max;
                } else if (mode==5u) {
                    for (std::size_t i=0u;i<count;i+=2u)
                        words[i]=0u;
                } else if (mode==6u) {
                    for (std::size_t i=0u;i<count;i+=3u)
                        words[i]=max;
                }
                const BigUnsigned value(words);
                const BigUnsigned by(divisor);
                const auto limbs=value.limbs();
                mpz_import(n.value,limbs.size(),-1,sizeof(std::uint64_t),
                           0,0,limbs.data());
                mpz_import(d.value,1u,-1,sizeof(divisor),0,0,&divisor);
                mpz_tdiv_qr(q.value,r.value,n.value,d.value);
                check_against_gmp(rem(value,by),r.value);
                const auto [full_q,full_r]=divrem(value,by);
                check_against_gmp(full_q,q.value);
                check_against_gmp(full_r,r.value);
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

    // B57: all 3 output modes of normalized multiprecision
    // division by 2^(64*k+s), including high/low limb boundaries
    // and dividend shorter than divisor. Verify against GMP both
    // as owning objects and with guarded caller-owned storage.
    std::mt19937_64 b57_power_engine(0xB570B172ULL);
    for (const std::size_t nlimbs:
            {1u,2u,8u,16u,32u,64u,65u,96u,128u}) {
        for (const std::size_t mlimbs:
                {1u,2u,3u,8u,16u,32u,64u,65u}) {
            for (const unsigned bit:{0u,1u,17u,32u,63u}) {
                for (unsigned mode=0u;mode<6u;++mode) {
                    std::vector<std::uint64_t> left(nlimbs);
                    std::vector<std::uint64_t> right(mlimbs,0u);
                    right.back()=std::uint64_t{1}<<bit;
                    for (auto& limb:left) limb=b57_power_engine();
                    if (mode==0u) {
                        std::fill(left.begin(),left.end(),max);
                    } else if (mode==1u) {
                        std::fill(left.begin(),left.end(),0u);
                    } else if (mode==2u) {
                        std::fill(left.begin(),left.end(),0u);
                        left.front()=1u;
                    } else if (mode==3u) {
                        std::fill(left.begin(),left.end(),0u);
                        left.back()=1u;
                    } else if (mode==4u) {
                        for (std::size_t i=0u;i<nlimbs;i+=2u)
                            left[i]=0u;
                    }
                    const BigUnsigned dividend(left),divisor(right);
                    mpz_import(a.value,dividend.size(),-1,
                               sizeof(std::uint64_t),0,0,
                               dividend.limbs().data());
                    mpz_import(b.value,divisor.size(),-1,
                               sizeof(std::uint64_t),0,0,
                               divisor.limbs().data());
                    mpz_tdiv_qr(q.value,r.value,a.value,b.value);
                    const auto [both_q,both_r]=
                        divrem(dividend,divisor);
                    check_against_gmp(both_q,q.value);
                    check_against_gmp(both_r,r.value);
                    check_against_gmp(quot(dividend,divisor),q.value);
                    check_against_gmp(rem(dividend,divisor),r.value);
                    const std::size_t needed=
                        BigUnsigned::divrem_scratch_limbs(
                            dividend.size(),divisor.size());
                    std::vector<std::uint64_t> qwords(
                        dividend.size()+2u,max);
                    std::vector<std::uint64_t> rwords(
                        std::min(dividend.size(),divisor.size())+2u,max);
                    std::vector<std::uint64_t> storage(
                        needed+2u,max);
                    BigUnsigned::divrem_into(
                        dividend.limbs(),divisor.limbs(),
                        qwords,rwords,BigUnsigned::MutableLimbView(
                            storage.data()+1u,needed));
                    check_against_gmp(BigUnsigned(qwords),q.value);
                    check_against_gmp(BigUnsigned(rwords),r.value);
                    ARIADNE_TEST_ASSERT(
                        qwords[dividend.size()]==0u &&
                        qwords[dividend.size()+1u]==0u);
                    const std::size_t m=
                        std::min(dividend.size(),divisor.size());
                    ARIADNE_TEST_ASSERT(rwords[m]==0u &&
                                        rwords[m+1u]==0u);
                    ARIADNE_TEST_ASSERT(storage.front()==max);
                    ARIADNE_TEST_ASSERT(storage.back()==max);
                }
            }
        }
    }

    // B55: GMP differential coverage of the native x86-64
    // MULQ/ADCQ digit-product subtraction selected at divisor
    // length 8. Preserve the 7-limb C++ control, 8/9 boundary,
    // 16/32/64-limb denominators, and both divisor operand orders.
    // Full/joint and single-output modes must independently agree
    // with GMP even after one- or two-digit corrections.
    std::mt19937_64 b55_engine(0xB550DEC1ULL);
    for (const std::size_t m:{7u,8u,9u,16u,32u,64u}) {
        for (const std::size_t extra:{2u,9u,32u}) {
            const std::size_t n=m+extra;
            for (unsigned mode=0u;mode<24u;++mode) {
                std::vector<std::uint64_t> left(n),right(m);
                for (auto& word:left) word=b55_engine();
                for (auto& word:right) word=b55_engine();
                switch (mode%8u) {
                case 0u:
                    std::fill(left.begin(),left.end(),max);
                    std::fill(right.begin(),right.end(),max);
                    break;
                case 1u:
                    std::fill(left.begin(),left.end(),0u);
                    left.back()=max;
                    std::fill(right.begin(),right.end(),0u);
                    right.back()=max;
                    break;
                case 2u:
                    std::fill(left.begin(),left.end(),0u);
                    left.front()=max;
                    left.back()=max;
                    std::fill(right.begin(),right.end(),0u);
                    right.front()=3u;
                    right.back()=max;
                    break;
                case 3u:
                    for (std::size_t i=0u;i<n;i+=2u) left[i]=0u;
                    for (std::size_t i=0u;i<m;i+=2u) right[i]=0u;
                    break;
                case 4u:
                    for (std::size_t i=1u;i+1u<n;++i) left[i]=0u;
                    for (std::size_t i=1u;i+1u<m;++i) right[i]=0u;
                    break;
                case 5u:
                    std::fill_n(left.begin(),n/2u,0u);
                    std::fill_n(right.begin(),m/2u,0u);
                    break;
                case 6u:
                    for (std::size_t i=1u;i<n;i+=3u) left[i]=max;
                    for (std::size_t i=1u;i<m;i+=3u) right[i]=max;
                    break;
                default: break;
                }
                left.back()|=std::uint64_t{1}<<63u;
                right.back()|=std::uint64_t{1}<<63u;
                const BigUnsigned dividend(left),divisor(right);
                mpz_import(a.value,dividend.size(),-1,
                           sizeof(std::uint64_t),0,0,
                           dividend.limbs().data());
                mpz_import(b.value,divisor.size(),-1,
                           sizeof(std::uint64_t),0,0,
                           divisor.limbs().data());
                mpz_tdiv_qr(q.value,r.value,a.value,b.value);
                const auto [actual_q,actual_r]=
                    divrem(dividend,divisor);
                check_against_gmp(actual_q,q.value);
                check_against_gmp(actual_r,r.value);
                check_against_gmp(quot(dividend,divisor),q.value);
                check_against_gmp(rem(dividend,divisor),r.value);
                const auto needed=BigUnsigned::divrem_scratch_limbs(
                    dividend.size(),divisor.size());
                std::vector<std::uint64_t> qwords(
                    dividend.size()+2u,max);
                std::vector<std::uint64_t> rwords(
                    divisor.size()+2u,max);
                std::vector<std::uint64_t> workspace(
                    needed+2u,max);
                BigUnsigned::divrem_into(
                    dividend.limbs(),divisor.limbs(),qwords,rwords,
                    BigUnsigned::MutableLimbView(
                        workspace.data()+1u,needed));
                check_against_gmp(BigUnsigned(qwords),q.value);
                check_against_gmp(BigUnsigned(rwords),r.value);
                ARIADNE_TEST_ASSERT(qwords[dividend.size()]==0u);
                ARIADNE_TEST_ASSERT(qwords[dividend.size()+1u]==0u);
                ARIADNE_TEST_ASSERT(rwords[divisor.size()]==0u);
                ARIADNE_TEST_ASSERT(rwords[divisor.size()+1u]==0u);
                ARIADNE_TEST_ASSERT(workspace.front()==max);
                ARIADNE_TEST_ASSERT(workspace.back()==max);
            }
        }
    }

    // B27: exercise the fused subtract/add-back implementation,
    // including quotient estimates that need TWO corrections.
    // These fixtures were found by the independent GMP stress probe,
    // not by adjusting the expected result to the native algorithm.
    const auto check_fused_division=[&](
            BigUnsigned const& numerator,BigUnsigned const& denominator) {
        mpz_import(a.value,numerator.size(),-1,sizeof(std::uint64_t),
                   0,0,numerator.limbs().data());
        mpz_import(b.value,denominator.size(),-1,sizeof(std::uint64_t),
                   0,0,denominator.limbs().data());
        mpz_tdiv_qr(q.value,r.value,a.value,b.value);
        const auto [native_q,native_r]=divrem(numerator,denominator);
        check_against_gmp(native_q,q.value);
        check_against_gmp(native_r,r.value);
        check_against_gmp(quot(numerator,denominator),q.value);
        check_against_gmp(rem(numerator,denominator),r.value);
        const auto needed=BigUnsigned::divrem_scratch_limbs(
            numerator.size(),denominator.size());
        if (denominator.size()>=2u &&
                numerator.size()>=denominator.size()) {
            ARIADNE_TEST_ASSERT(needed==
                numerator.size()+denominator.size()+1u);
        }
        ARIADNE_TEST_ASSERT(compare(native_r,denominator)<0);
    };
    check_fused_division(
        BigUnsigned(std::array<std::uint64_t,4>{
            max,max,max,max}),
        BigUnsigned(std::array<std::uint64_t,2>{
            0xe04e3de411f558a7ULL,0xba43f747b13ac11aULL}));
    check_fused_division(
        BigUnsigned(std::array<std::uint64_t,5>{
            0x65a601688b16b11bULL,0x5b7d0bd5828e0a0cULL,
            0x7ad9be6c73e93fdaULL,0xf370dcf1b5187984ULL,
            0x48293917315898ccULL}),
        BigUnsigned(std::array<std::uint64_t,2>{
            0xc28ba8299ce45291ULL,0x91b35c306fef072dULL}));

    std::mt19937_64 fused_engine(0xB27071ULL);
    for (unsigned repetition=0u;repetition<2048u;++repetition) {
        const std::size_t n=repetition%31u==0u?96u:
                            2u+fused_engine()%32u;
        const std::size_t m=repetition%17u==0u?2u:
                            2u+fused_engine()%n;
        std::vector<std::uint64_t> left(n),right(m);
        for (auto& word:left) word=fused_engine();
        for (auto& word:right) word=fused_engine();
        if (repetition%13u==0u) std::fill(left.begin(),left.end(),max);
        if (repetition%19u==0u) {
            std::fill(right.begin(),right.end(),max);
        }
        if (repetition%23u==0u) {
            std::fill(right.begin(),right.end(),0u);
            right.back()=1u;
        }
        if (repetition%29u==0u) {
            std::fill(right.begin(),right.end(),0u);
            right.back()=std::uint64_t{1}<<63u;
        }
        check_fused_division(BigUnsigned(left),BigUnsigned(right));
    }

    // B48: preserve quotient-only, remainder-only and joint-output
    // semantics for the small divisors used in the unchanged
    // 32-case official benchmark; GMP remains the independent oracle.
    std::mt19937_64 small_divisor_engine(0xB48A11D1ULL);
    for (const std::size_t count:
            {2u,4u,5u,8u,16u,32u,64u,96u}) {
        for (const std::uint64_t divisor:
                {1ULL,3ULL,5ULL,7ULL,9ULL,16ULL,31ULL,
                 34ULL,63ULL,255ULL,65535ULL,4294967295ULL,
                 4294967296ULL,18446744073709551615ULL}) {
            for (unsigned mode=0u;mode<3u;++mode) {
                std::vector<std::uint64_t> words(count);
                for (auto& limb:words) limb=small_divisor_engine();
                if (mode==0u) {
                    std::fill(words.begin(),words.end(),max);
                } else if (mode==1u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=1u;
                }
                check_fused_division(
                    BigUnsigned(words),BigUnsigned(divisor));
            }
        }
    }

    // B28: equal highest limbs guarantee a quotient of either 0 or 1.
    // Exercise both orderings, equality, borrow propagation and long
    // exact remainders, independently against GMP in all output modes.
    std::mt19937_64 equal_high_engine(0xB280C27ULL);
    for (unsigned repetition=0u;repetition<2048u;++repetition) {
        const std::size_t count=repetition%37u==0u?128u:
                                2u+equal_high_engine()%64u;
        std::vector<std::uint64_t> left(count),right(count);
        for (auto& word:left) word=equal_high_engine();
        for (auto& word:right) word=equal_high_engine();
        right.back()|=1u;
        left.back()=right.back();
        if (repetition%4u==0u) {
            left=right;
        } else if (repetition%4u==1u) {
            left=right;
            ++left.front();
        } else if (repetition%4u==2u) {
            right=left;
            ++right.front();
        }
        const BigUnsigned numerator(left),denominator(right);
        check_fused_division(numerator,denominator);

        // The public joint-output contract still requires the same
        // scratch size and clears surplus output limbs. This fast path
        // must not write into that scratch (including its guard limbs).
        const auto needed=BigUnsigned::divrem_scratch_limbs(
            numerator.size(),denominator.size());
        std::vector<std::uint64_t> qout(numerator.size()+2u,max);
        std::vector<std::uint64_t> rout(denominator.size()+2u,max);
        std::vector<std::uint64_t> guards(needed+2u,max);
        BigUnsigned::divrem_into(
            numerator.limbs(),denominator.limbs(),qout,rout,
            BigUnsigned::MutableLimbView(guards.data()+1u,needed));
        check_against_gmp(BigUnsigned(qout),q.value);
        check_against_gmp(BigUnsigned(rout),r.value);
        for (const auto word:guards) {
            ARIADNE_TEST_ASSERT(word==max);
        }
        for (std::size_t i=numerator.size();i<qout.size();++i) {
            ARIADNE_TEST_ASSERT(qout[i]==0u);
        }
        for (std::size_t i=denominator.size();i<rout.size();++i) {
            ARIADNE_TEST_ASSERT(rout[i]==0u);
        }
    }

    // B36: equal-size quotients that can be classified from the
    // most-significant limbs without Knuth normalization. GMP is the
    // independent oracle for q=0, q=1, q>=2, and boundary ratios.
    // Exercise all three public division modes via check_fused_division.
    std::mt19937_64 high_ratio_engine(0xB360A11ULL);
    for (const std::size_t count:
            {2u,3u,4u,7u,8u,16u,32u,64u,96u,128u}) {
        for (unsigned repetition=0u;repetition<128u;++repetition) {
            std::vector<std::uint64_t> left(count),right(count);
            for (auto& limb:left) limb=high_ratio_engine();
            for (auto& limb:right) limb=high_ratio_engine();
            // These pairs are canonical; low words can be saturated,
            // sparse or identical regardless of the chosen high ratio.
            if (repetition%11u==0u) {
                std::fill(left.begin(),left.end(),max);
                std::fill(right.begin(),right.end(),max);
            } else if (repetition%11u==1u) {
                std::fill(left.begin(),left.end(),0u);
                std::fill(right.begin(),right.end(),0u);
            }
            const std::uint64_t high=std::uint64_t{1}<<63u;
            switch (repetition%8u) {
            case 0u: // q=0: strictly smaller top limb.
                left.back()=high-1u; right.back()=high;
                break;
            case 1u: // q=1: strictly greater top, below 2*bh.
                left.back()=high+123u; right.back()=high;
                break;
            case 2u: // q=0 or 1: equal top, B28 path.
                left.back()=high; right.back()=high;
                break;
            case 3u: // q>=2: outside the new fast path.
                left.back()=high; right.back()=2u;
                break;
            case 4u: // q=1 at the largest possible high limb.
                left.back()=max; right.back()=high;
                break;
            case 5u: // boundary ah==2*bh, must use long division.
                left.back()=high; right.back()=high/2u;
                break;
            case 6u: // q=0: top of divisor is saturated.
                left.back()=high; right.back()=max;
                break;
            default: // arbitrary canonical high limbs.
                left.back()|=1u; right.back()|=1u;
                break;
            }
            check_fused_division(BigUnsigned(left),BigUnsigned(right));
        }
    }

    // B29: deterministic sparse normalized divisors with long runs
    // of zero limbs. Verify both the sparse and dense selection shapes
    // against GMP through joint, quotient-only and remainder-only modes.
    std::mt19937_64 sparse_engine(0xB290D1F5ULL);
    for (unsigned repetition=0u;repetition<1024u;++repetition) {
        const std::size_t n=16u+sparse_engine()%113u;
        const std::size_t m=8u+sparse_engine()%(n-7u);
        std::vector<std::uint64_t> left(n),right(m,0u);
        for (auto& word:left) word=sparse_engine();
        right.front()=sparse_engine();
        right.back()=std::uint64_t{1}<<(repetition%64u);
        if (repetition%3u==0u) {
            right[m/2u]=sparse_engine();
        }
        if (repetition%5u==0u) {
            // Partially dense divisors must retain the original kernel.
            for (std::size_t i=1u;i+1u<m;++i) {
                if (i%3u==0u) right[i]=sparse_engine();
            }
        }
        if (repetition%7u==0u) {
            std::fill(left.begin(),left.end(),max);
        }
        if (repetition%11u==0u) {
            left.back()=1u;
        }
        check_fused_division(BigUnsigned(left),BigUnsigned(right));
    }

    // B45: independent GMP checks of whole-zero low divisor limbs.
    // Include inline/heap boundaries, equal/small/large quotients,
    // arbitrary low remainders, all output modes and guard contracts.
    std::mt19937_64 shifted_engine(0xB450D1D1ULL);
    for (const std::size_t m:{2u,3u,5u,8u,16u,32u,64u}) {
        for (const std::size_t zeros:
                {std::size_t{1},m/2u,m-1u}) {
            for (unsigned repetition=0u;repetition<24u;++repetition) {
                const std::size_t n=m+
                    (repetition%4u==0u?0u:
                     repetition%4u==1u?1u:
                     repetition%4u==2u?m:2u*m);
                std::vector<std::uint64_t> left(n),right(m);
                for (auto& word:left) word=shifted_engine();
                for (auto& word:right) word=shifted_engine();
                std::fill_n(right.begin(),zeros,0u);
                right.back()|=std::uint64_t{1}<<63u;
                left.back()|=std::uint64_t{1}<<63u;
                if (repetition%6u==0u && n==m) {
                    left=right;
                } else if (repetition%6u==1u && n==m) {
                    left.back()=right.back()-1u;
                } else if (repetition%6u==2u) {
                    std::fill_n(left.begin(),zeros,0u);
                } else if (repetition%6u==3u) {
                    std::fill(left.begin(),left.end(),max);
                }
                const BigUnsigned numerator(left),denominator(right);
                check_fused_division(numerator,denominator);
                const auto needed=BigUnsigned::divrem_scratch_limbs(
                    numerator.size(),denominator.size());
                std::vector<std::uint64_t> qout(
                    numerator.size()+2u,max);
                std::vector<std::uint64_t> rout(
                    std::min(numerator.size(),denominator.size())+2u,max);
                std::vector<std::uint64_t> guards(needed+2u,max);
                BigUnsigned::divrem_into(
                    numerator.limbs(),denominator.limbs(),qout,rout,
                    BigUnsigned::MutableLimbView(
                        guards.data()+1u,needed));
                check_against_gmp(BigUnsigned(qout),q.value);
                check_against_gmp(BigUnsigned(rout),r.value);
                for (std::size_t i=numerator.size();i<qout.size();++i) {
                    ARIADNE_TEST_ASSERT(qout[i]==0u);
                }
                const std::size_t rem_count=
                    std::min(numerator.size(),denominator.size());
                for (std::size_t i=rem_count;i<rout.size();++i) {
                    ARIADNE_TEST_ASSERT(rout[i]==0u);
                }
                ARIADNE_TEST_ASSERT(guards.front()==max);
                ARIADNE_TEST_ASSERT(guards.back()==max);
            }
        }
    }

    // B45: the decimal divisors of B22 recursion contain large
    // powers of two; exercise the exact official sparse input shape.
    for (const unsigned digits:{72u,144u,288u,576u}) {
        const BigUnsigned decimal_power=from_decimal(
            std::string("1")+std::string(digits,'0'));
        for (const std::size_t bits:{2048u,4096u}) {
            const BigUnsigned binary_power=shift_left(
                BigUnsigned(std::uint64_t{1}),bits-1u);
            for (const std::uint64_t offset:
                    {0u,1u,13u,127u,1023u}) {
                const BigUnsigned numerator=add(
                    binary_power,BigUnsigned(offset));
                check_fused_division(numerator,decimal_power);
            }
        }
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
    // B24: long-input radix-10^18 parsing. Probe the historical
    // 309/310-digit boundary (both now wide after B47), 18-digit
    // group carry boundaries, and
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
    // B47: the wide input path now starts at 32 digits, including
    // 128-, 256- and 1024-bit benchmark magnitudes. Stress first
    // short groups, full 18-digit boundaries and both sides of
    // the old 309/310 threshold against independently parsed GMP.
    for (const std::size_t length:
            {16u,17u,18u,19u,31u,32u,33u,35u,36u,37u,38u,
             39u,40u,53u,54u,55u,76u,77u,78u,79u,
             307u,308u,309u,310u,311u,617u}) {
        check_long_parse(std::string(length,'9'));
        check_long_parse(std::string(length,'0'));
        check_long_parse("1"+std::string(length-2u,'0')+"1");
        check_long_parse(std::string(length-1u,'0')+"7");
        std::string alternating(length,'0');
        for (std::size_t i=0u;i<length;i+=2u) {
            alternating[i]='9';
        }
        check_long_parse(alternating);
        std::string random_digits(length,'0');
        for (auto& digit:random_digits) {
            digit=static_cast<char>('0'+engine()%10u);
        }
        check_long_parse(random_digits);
    }
    // B50: radix-10^19 parsing must handle all boundaries of
    // complete 19-digit groups, the 31/32-digit selector, the
    // 308-digit 1024-bit input, and 2048/4096-bit inputs. Check
    // exact parsed magnitudes against GMP, including leading zeros,
    // all nines, sparse digits, and independently randomized text.
    for (const std::size_t length:
            {18u,19u,20u,31u,32u,33u,36u,37u,38u,39u,40u,
             56u,57u,58u,75u,76u,77u,78u,94u,95u,96u,
             113u,114u,115u,151u,152u,153u,307u,308u,
             309u,310u,311u,616u,617u,618u,619u,
             1233u,1234u,1235u}) {
        check_long_parse(std::string(length,'9'));
        check_long_parse(std::string(length,'0'));
        check_long_parse("1"+std::string(length-2u,'0')+"1");
        check_long_parse(std::string(length-1u,'0')+"7");
        std::string alternating(length,'0');
        for (std::size_t i=0u;i<length;i+=2u) alternating[i]='9';
        check_long_parse(alternating);
        std::string random_digits(length,'0');
        for (auto& digit:random_digits) {
            digit=static_cast<char>('0'+engine()%10u);
        }
        check_long_parse(random_digits);
    }

    // B47: exact powers 2^(bits-1) plus offsets, as used by
    // construct_decimal, checked with GMP independent of the
    // benchmark; include both signs via the established Integer
    // benchmark verification.
    for (const unsigned bits:{128u,256u,1024u,2048u,4096u}) {
        mpz_set_ui(a.value,1u);
        mpz_mul_2exp(a.value,a.value,bits-1u);
        for (const unsigned long offset:{0ul,13ul,31ul,75ul,187ul}) {
            mpz_add_ui(b.value,a.value,offset);
            const std::size_t buffer_size=mpz_sizeinbase(b.value,10u)+2u;
            std::string decimal(buffer_size,'\0');
            mpz_get_str(decimal.data(),10,b.value);
            decimal.resize(std::char_traits<char>::length(decimal.c_str()));
            check_long_parse(decimal);
        }
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

    // B30/B30.1: preserve the 8/9-limb boundary and 10^18
    // decimal carry fixtures against an independent GMP oracle.
    // ARM64 selects radix 10^18; x86-64 and other targets now
    // exercise the restored B29 radix-10^9 implementation.
    std::mt19937_64 divq_decimal_engine(0xB300DEC1ULL);
    for (const std::size_t count:
            {8u,9u,10u,12u,16u,17u,32u,64u}) {
        for (unsigned mode=0u;mode<16u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=divq_decimal_engine();
            if (mode%8u==0u) {
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()=1u;
            } else if (mode%8u==1u) {
                std::fill(words.begin(),words.end(),
                          std::numeric_limits<std::uint64_t>::max());
            } else if (mode%8u==2u) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=1u;
            } else if (mode%8u==3u) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::numeric_limits<std::uint64_t>::max();
            } else if (mode%8u==4u) {
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()=std::uint64_t{1}<<63u;
            } else if (mode%8u==5u) {
                for (std::size_t i=1u;i+1u<count;++i) {
                    words[i]=0u;
                }
            }
            words.back()|=1u;
            mpz_import(a.value,words.size(),-1,sizeof(std::uint64_t),
                       0,0,words.data());
            check_decimal_boundary(a.value);
        }
    }

    // B42: GMP-oracle regression at the 8/9/16/17-limb wide
    // decimal formatter boundaries, plus large dense/sparse inputs.
    std::mt19937_64 b42_engine(0xB42DEC1A7ULL);
    for (const std::size_t count:
            {8u,9u,10u,15u,16u,17u,18u,32u,33u,64u,65u,128u,256u}) {
        for (unsigned mode=0u;mode<16u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& w:words) w=b42_engine();
            switch (mode%8u) {
            case 0u: std::fill(words.begin(),words.end(),0u); break;
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
                words.back()=std::uint64_t{1}<<63u;
                break;
            case 4u:
                for (std::size_t i=0u;i<count;i+=2u) words[i]=0u;
                break;
            case 5u:
                std::fill(words.begin(),words.end(),
                          std::uint64_t{999999999999999999ULL});
                break;
            case 6u: words.back()=1u; break;
            default: break;
            }
            mpz_import(a.value,words.size(),-1,sizeof(std::uint64_t),
                       0,0,words.data());
            check_decimal_boundary(a.value);
        }
    }

    // B46: reciprocal-based ARM64 radix-10^18 leaf divrem replaces
    // the unsigned 128/64 runtime divide. Exercise random dense and
    // sparse 9..16-limb leaves, their 8/9 dispatch edge, the 16/17
    // recursion boundary, carries and long internal zero runs. On
    // every host compare against GMP, not a native roundtrip alone.
    std::mt19937_64 b46_engine(0xB460DEC1A7ULL);
    for (const std::size_t count:
            {8u,9u,10u,11u,12u,13u,14u,15u,16u,17u,32u}) {
        for (unsigned mode=0u;mode<32u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=b46_engine();
            if (mode%8u==0u) {
                std::fill(words.begin(),words.end(),
                          std::numeric_limits<std::uint64_t>::max());
            } else if (mode%8u==1u) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=1u;
            } else if (mode%8u==2u) {
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()=1u;
            } else if (mode%8u==3u) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
            } else if (mode%8u==4u) {
                std::fill_n(words.begin(),count/2u,0u);
            } else if (mode%8u==5u) {
                for (std::size_t j=0u;j<count;j+=2u) {
                    words[j]=0u;
                }
            } else if (mode%8u==6u) {
                std::fill(words.begin(),words.end(),
                          std::uint64_t{999999999999999999ULL});
            }
            words.back()|=1u;
            mpz_import(a.value,words.size(),-1,sizeof(std::uint64_t),
                       0,0,words.data());
            check_decimal_boundary(a.value);
        }
    }

    // B51: the multiply-only B46 wide leaf is now also selected
    // for x86-64 GCC/Clang at 9+ limbs. Independent GMP oracles
    // cover the 8/9-limb selector, the 16/17-limb recursive cutoff,
    // dense/sparse operands, carry boundaries and zero runs. Keep
    // the portable/MSVC radix-10^9 regression path equally checked.
    std::mt19937_64 b51_engine(0xB510D3C1ULL);
    for (const std::size_t count:
            {7u,8u,9u,10u,12u,15u,16u,17u,18u,32u,33u,64u}) {
        for (unsigned mode=0u;mode<24u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=b51_engine();
            switch (mode%8u) {
            case 0u:
                std::fill(words.begin(),words.end(),
                          std::numeric_limits<std::uint64_t>::max());
                break;
            case 1u:
                std::fill(words.begin(),words.end(),0u);
                words.back()=1u;
                break;
            case 2u:
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()=1u;
                break;
            case 3u:
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
                break;
            case 4u:
                std::fill_n(words.begin(),count/2u,0u);
                break;
            case 5u:
                for (std::size_t i=0u;i<count;i+=2u) {
                    words[i]=0u;
                }
                break;
            case 6u:
                std::fill(words.begin(),words.end(),
                          std::uint64_t{999999999999999999ULL});
                break;
            default: break;
            }
            words.back()|=1u;
            mpz_import(a.value,words.size(),-1,sizeof(std::uint64_t),
                       0,0,words.data());
            check_decimal_boundary(a.value);
        }
    }

    // B59: oracle-checked x86 GCC 16 -> 32-limb recursion cutoff,
    // especially 31/32/33-limb transition and 63/64/65 inputs.
    // Other hosts use the unmodified 16-limb selector. Reuse the
    // independent GMP binary-export + mpz_get_str oracle.
    std::mt19937_64 b59_engine(0xB590DEC1ULL);
    for (const std::size_t count:
            {15u,16u,17u,31u,32u,33u,63u,64u,65u}) {
        for (unsigned mode=0u;mode<24u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=b59_engine();
            switch (mode%8u) {
            case 0u:
                std::fill(words.begin(),words.end(),0u);
                words.back()=1u;
                break;
            case 1u:
                std::fill(words.begin(),words.end(),
                          std::numeric_limits<std::uint64_t>::max());
                break;
            case 2u:
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()=std::numeric_limits<std::uint64_t>::max();
                break;
            case 3u:
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
                break;
            case 4u:
                for (std::size_t i=0u;i<count;i+=2u) words[i]=0u;
                break;
            case 5u:
                std::fill(words.begin(),words.end(),
                          std::uint64_t{999999999999999999ULL});
                break;
            case 6u:
                for (std::size_t i=1u;i+1u<count;++i) words[i]=0u;
                break;
            default: break;
            }
            // Keep the raw magnitude at the intended selector size.
            words.back()|=std::uint64_t{1}<<63u;
            mpz_import(a.value,words.size(),-1,sizeof(std::uint64_t),
                       0,0,words.data());
            check_decimal_boundary(a.value);
        }
    }

    // B43: append-only recursive formatting must preserve the
    // *entire* fixed-width lower field even when the remainder
    // itself needs further recursive divisions. GMP independently
    // supplies the reference, never a native decimal roundtrip.
    for (const unsigned long digits:
            {144ul,288ul,576ul,1152ul,2304ul}) {
        mpz_ui_pow_ui(a.value,10ul,digits);
        for (const unsigned long factor:
                {1ul,7ul,1234567ul,9999999ul}) {
            mpz_mul_ui(b.value,a.value,factor);
            for (const unsigned long offset:
                    {0ul,1ul,11ul,99ul}) {
                mpz_add_ui(b.value,b.value,offset);
                check_decimal_boundary(b.value);
                mpz_sub_ui(b.value,b.value,offset);
            }
        }
    }

    // B31: the new pair-wise base-10^9 renderer must preserve every
    // zero-padded two-digit group (00, 01, 09, 10, 99), carry across
    // decimal blocks, and the variable-width most-significant block.
    // GMP supplies the independent exact decimal oracle on every host.
    for (const unsigned long digits:
            {9ul,18ul,27ul,72ul,144ul,288ul,576ul}) {
        mpz_ui_pow_ui(a.value,10ul,digits);
        for (const unsigned long offset:
                {0ul,1ul,9ul,10ul,11ul,99ul,100ul,101ul,999ul}) {
            mpz_set(b.value,a.value);
            mpz_add_ui(b.value,b.value,offset);
            check_decimal_boundary(b.value);
            if (offset!=0ul) {
                mpz_sub_ui(b.value,a.value,offset);
                check_decimal_boundary(b.value);
            }
        }
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
