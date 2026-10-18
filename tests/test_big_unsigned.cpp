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
    const auto digits=x.limbs();
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
    }
  private:
    void test_storage();
    void test_boundaries();
    void test_randomized();
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
    bool threw=false;
    try { (void)sub(one,all); } catch (std::underflow_error const&) { threw=true; }
    ARIADNE_TEST_ASSERT(threw);
    ARIADNE_TEST_ASSERT(compare(add(all,BigUnsigned{}),all)==0);
}

void TestBigUnsigned::test_randomized()
{
    std::mt19937_64 engine(0xBA17B00DULL);
    GmpInteger a,b,c;
    for (unsigned repetition=0; repetition<10000u; ++repetition) {
        const std::size_t n=engine()%20u, m=engine()%20u;
        std::array<std::uint64_t,20> na{},nb{};
        for(std::size_t j=0;j<n;++j) na[j]=engine();
        for(std::size_t j=0;j<m;++j) nb[j]=engine();
        if (repetition%7u==0u) na.fill(std::numeric_limits<std::uint64_t>::max());
        if (repetition%11u==0u) nb.fill(0u);
        const BigUnsigned x(std::span(na.data(),n));
        const BigUnsigned y(std::span(nb.data(),m));
        mpz_import(a.value,n,-1,8,0,0,na.data());
        mpz_import(b.value,m,-1,8,0,0,nb.data());
        ARIADNE_TEST_ASSERT(compare(x,y)==mpz_cmp(a.value,b.value) ||
            (compare(x,y)<0 && mpz_cmp(a.value,b.value)<0) ||
            (compare(x,y)>0 && mpz_cmp(a.value,b.value)>0));
        mpz_add(c.value,a.value,b.value);
        check_against_gmp(add(x,y),c.value);
        if(mpz_cmp(a.value,b.value)>=0) {
            mpz_sub(c.value,a.value,b.value);
            check_against_gmp(sub(x,y),c.value);
        }
        const std::size_t bits=engine()%1030u;
        mpz_mul_2exp(c.value,a.value,bits);
        check_against_gmp(shift_left(x,bits),c.value);
        mpz_fdiv_q_2exp(c.value,a.value,bits);
        check_against_gmp(shift_right(x,bits),c.value);
    }
}

} // namespace

int main()
{
    ARIADNE_TEST_CLASS(BigUnsigned,TestBigUnsigned());
    return ARIADNE_TEST_FAILURES;
}
