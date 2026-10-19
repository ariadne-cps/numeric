/***************************************************************************
 *            numeric/big_unsigned.cpp
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

/*! \file numeric/big_unsigned.cpp
 *  \brief Shared unsigned arithmetic, fused multiplication, and Karatsuba.
 */

#include "numeric/big_unsigned.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <charconv>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#if defined(_MSC_VER) && defined(_M_X64) && !defined(__clang__)
#include <immintrin.h>
#endif

namespace Ariadne {

#if defined(ARIADNE_BIG_UNSIGNED_X86_64)
namespace BigUnsignedDetail {
void mul_low_into_x86_64(BigUnsigned::LimbView a,
                         BigUnsigned::LimbView b,
                         BigUnsigned::MutableLimbView output) noexcept;
void mul_low_into_comba_x86_64(BigUnsigned::LimbView a,
                               BigUnsigned::LimbView b,
                               BigUnsigned::MutableLimbView output) noexcept;
void mul_low_into_carry_x86_64(BigUnsigned::LimbView a,
                               BigUnsigned::LimbView b,
                               BigUnsigned::MutableLimbView output) noexcept;
} // namespace BigUnsignedDetail
#endif

#if defined(ARIADNE_BIG_UNSIGNED_ARM64)
namespace BigUnsignedDetail {
void mul_low_into_carry_arm64(BigUnsigned::LimbView a,
                              BigUnsigned::LimbView b,
                              BigUnsigned::MutableLimbView output) noexcept;
} // namespace BigUnsignedDetail
#endif

BigUnsigned::MutableLimbView BigUnsigned::Scratch::acquire(std::size_t count) {
    if (count <= inline_capacity) {
        std::fill_n(_inline_limbs.data(),count,Limb{0});
        return {_inline_limbs.data(),count};
    }

    // A new request reuses existing vector capacity whenever possible.
    _heap.resize(count);
    std::fill(_heap.begin(),_heap.end(),Limb{0});
    return {_heap.data(),count};
}

std::size_t BigUnsigned::Scratch::capacity() const noexcept {
    return std::max(inline_capacity,_heap.capacity());
}

BigUnsigned::BigUnsigned(Limb value) noexcept {
    _inline_limbs[0] = value;
    _size = value != 0u ? 1u : 0u;
}

BigUnsigned::BigUnsigned(std::span<Limb const> digits) {
    set_limbs(digits);
}

BigUnsigned::BigUnsigned(BigUnsigned&& other) noexcept
    : _inline_limbs(other._inline_limbs),
      _heap(std::move(other._heap)),
      _size(other._size)
{
    other._size = 0u;
}

BigUnsigned& BigUnsigned::operator=(BigUnsigned&& other) noexcept {
    if (this != &other) {
        _inline_limbs = other._inline_limbs;
        _heap = std::move(other._heap);
        _size = other._size;

        other._size = 0u;
    }

    return *this;
}

BigUnsigned::Limb* BigUnsigned::mutable_data() noexcept {
    return _size <= inline_capacity ? _inline_limbs.data() : _heap.data();
}

void BigUnsigned::resize_for_result(std::size_t limbs) {
    if (limbs > inline_capacity) {
        _heap.resize(limbs,0u);
    } else {
        _inline_limbs.fill(0u);
    }

    _size = limbs;
}

void BigUnsigned::normalize() noexcept {
    auto* digits = mutable_data();

    while (_size != 0u && digits[_size-1u] == 0u) {
        --_size;
    }

    if (_size <= inline_capacity && !_heap.empty()) {
        // Read the significant digits before releasing heap storage.
        std::copy_n(_heap.data(),_size,_inline_limbs.data());
        _heap.clear();
    }
}

BigUnsigned::LimbView BigUnsigned::limbs() const noexcept {
    return _size <= inline_capacity
        ? std::span<Limb const>(_inline_limbs.data(),_size)
        : std::span<Limb const>(_heap.data(),_size);
}

BigUnsigned::Limb BigUnsigned::limb(std::size_t index) const noexcept {
    return index < _size ? limbs()[index] : 0u;
}

void BigUnsigned::assign(LimbView digits) {
    set_limbs(digits);
}

void BigUnsigned::set_limbs(std::span<Limb const> digits) {
    // The empty sequence is the unique representation of zero.
    while (!digits.empty() && digits.back() == 0u) {
        digits = digits.first(digits.size()-1u);
    }

    if (digits.size() <= inline_capacity) {
        // Copy before mutating state: digits may refer to our own storage.
        std::array<Limb,inline_capacity> replacement{};
        std::copy(digits.begin(),digits.end(),replacement.begin());

        _inline_limbs = replacement;
        _size = digits.size();
        _heap.clear();
    } else {
        // Construct first so allocation failures leave this value intact.
        std::vector<Limb> replacement(digits.begin(),digits.end());

        _heap.swap(replacement);
        _size = digits.size();
    }
}

std::size_t BigUnsigned::bit_length() const noexcept {
    if (_size == 0u) {
        return 0u;
    }

    return 64u * (_size-1u) +
        64u - static_cast<unsigned>(std::countl_zero(limb(_size-1u)));
}

void BigUnsigned::swap(BigUnsigned& other) noexcept {
    using std::swap;

    swap(_inline_limbs,other._inline_limbs);
    _heap.swap(other._heap);
    swap(_size,other._size);
}

int compare(BigUnsigned const& a, BigUnsigned const& b) noexcept {
    if (a._size != b._size) {
        return a._size < b._size ? -1 : 1;
    }

    // Choose inline or heap storage just once per operand. In
    // particular, avoid repeated limb() bounds/storage checks while
    // scanning long equal prefixes from the most significant limb.
    const auto left=a.limbs();
    const auto right=b.limbs();
    for (std::size_t i=left.size(); i!=0u; --i) {
        const auto lhs=left[i-1u];
        const auto rhs=right[i-1u];
        if (lhs!=rhs) {
            return lhs<rhs ? -1 : 1;
        }
    }

    return 0;
}

BigUnsigned::Limb BigUnsigned::add_into(
        LimbView a,LimbView b,MutableLimbView output)
{
    const std::size_t count = std::max(a.size(),b.size());
    if (output.size() < count) {
        throw std::length_error("BigUnsigned add_into output too small");
    }

    Limb carry = 0u;
    for (std::size_t i=0u; i<count; ++i) {
        const Limb left = i<a.size() ? a[i] : 0u;
        const Limb right = i<b.size() ? b[i] : 0u;
        const Limb sum = left+right;
        const Limb value = sum+carry;

        carry = static_cast<Limb>(sum<left || value<sum);
        output[i] = value;
    }

    std::fill(output.subspan(count).begin(),output.end(),Limb{0});
    return carry;
}

BigUnsigned::Limb BigUnsigned::sub_into(
        LimbView a,LimbView b,MutableLimbView output)
{
    const std::size_t count = std::max(a.size(),b.size());
    if (output.size() < count) {
        throw std::length_error("BigUnsigned sub_into output too small");
    }

    Limb borrow = 0u;
    for (std::size_t i=0u; i<count; ++i) {
        const Limb left = i<a.size() ? a[i] : 0u;
        const Limb right = i<b.size() ? b[i] : 0u;
        const Limb difference = left-right;
        const Limb value = difference-borrow;

        borrow = static_cast<Limb>(left<right || difference<borrow);
        output[i] = value;
    }

    std::fill(output.subspan(count).begin(),output.end(),Limb{0});
    return borrow;
}

BigUnsigned add(BigUnsigned const& a, BigUnsigned const& b) {
    const std::size_t count = std::max(a._size,b._size);
    if (count == std::numeric_limits<std::size_t>::max()) {
        throw std::length_error("BigUnsigned addition size overflow");
    }

    BigUnsigned result;
    result.resize_for_result(count+1u);
    auto out = BigUnsigned::MutableLimbView(result.mutable_data(),count+1u);
    out[count] = BigUnsigned::add_into(a.limbs(),b.limbs(),out.first(count));
    result.normalize();
    return result;
}

BigUnsigned sub(BigUnsigned const& a, BigUnsigned const& b) {
    if (compare(a,b)<0) {
        throw std::underflow_error("BigUnsigned negative result");
    }

    BigUnsigned result;
    result.resize_for_result(a._size);
    const auto out = BigUnsigned::MutableLimbView(result.mutable_data(),a._size);
    (void)BigUnsigned::sub_into(a.limbs(),b.limbs(),out);
    result.normalize();
    return result;
}

BigUnsigned shift_left(BigUnsigned const& a, std::size_t bits) {
    if (a.is_zero()) {
        return BigUnsigned{};
    }

    const std::size_t words = bits/64u;
    const unsigned offset = static_cast<unsigned>(bits%64u);

    if (words >= std::numeric_limits<std::size_t>::max()-a._size-1u) {
        throw std::length_error("BigUnsigned left shift size overflow");
    }

    BigUnsigned result;
    result.resize_for_result(a._size+words+1u);
    auto* digits = result.mutable_data();

    BigUnsigned::Limb carry = 0u;
    for (std::size_t i=0u; i<a._size; ++i) {
        const auto value = a.limb(i);

        digits[words+i] = (value<<offset)|carry;
        carry = offset == 0u ? 0u : (value>>(64u-offset));
    }

    digits[words+a._size] = carry;
    result.normalize();

    return result;
}

BigUnsigned shift_right(BigUnsigned const& a, std::size_t bits) {
    const std::size_t words = bits/64u;

    if (words >= a._size) {
        return BigUnsigned{};
    }

    const unsigned offset = static_cast<unsigned>(bits%64u);
    BigUnsigned result;
    result.resize_for_result(a._size-words);
    auto* digits = result.mutable_data();

    for (std::size_t i=words; i<a._size; ++i) {
        const auto low = a.limb(i);
        const auto high = i+1u < a._size ? a.limb(i+1u) : 0u;

        digits[i-words] = offset == 0u ? low :
            ((low>>offset)|(high<<(64u-offset)));
    }

    result.normalize();

    return result;
}

// Divide a 128-bit (high:low) numerator by a nonzero 64-bit divisor.
// Caller invariant: high < divisor, so the quotient fits one limb.
// Used for both single-limb division and multiprecision quotient estimates.
namespace {

// Retain the correctness-first, ISA-independent C++20 fallback.
[[maybe_unused]] std::pair<std::uint64_t,std::uint64_t> divide_wide_portable(
        std::uint64_t high,std::uint64_t low,std::uint64_t divisor) noexcept
{
    std::uint64_t remainder=high;
    std::uint64_t quotient=0u;
    for (unsigned bit=64u;bit!=0u;--bit) {
        const std::uint64_t overflow=remainder>>63u;
        remainder=(remainder<<1u)|((low>>(bit-1u))&std::uint64_t{1});
        if (overflow!=0u || remainder>=divisor) {
            remainder-=divisor;
            quotient|=std::uint64_t{1}<<(bit-1u);
        }
    }
    return {quotient,remainder};
}

std::pair<std::uint64_t,std::uint64_t> divide_wide(
        std::uint64_t high,std::uint64_t low,std::uint64_t divisor) noexcept
{
    // Dividing a single 64-bit limb does not require wide arithmetic.
    if (high==0u) {
        return {low/divisor,low%divisor};
    }
#if defined(_MSC_VER) && defined(_M_X64) && !defined(__clang__)
    // MSVC x64 intrinsic available since Visual Studio 2019.
    unsigned __int64 remainder=0u;
    const unsigned __int64 quotient=_udiv128(
        high,low,divisor,&remainder);
    return {quotient,remainder};
#elif defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && \
        !defined(_MSC_VER)
    // B26: DIVQ divides RDX:RAX by a 64-bit operand, returning
    // quotient in RAX and remainder in RDX. The caller guarantees
    // 0 <= high < divisor, so the quotient fits in 64 bits and the
    // hardware divide cannot overflow. Both registers are read/write
    // operands; the divisor is distinct from them. No BMI2/ADX,
    // extra build flags or CPU dispatch are required.
    __asm__ volatile(
        "divq %[divisor]"
        : "+a" (low), "+d" (high)
        : [divisor] "r" (divisor)
        : "cc");
    return {low,high};
#elif defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    // The extension is supported here; suppress only -Wpedantic.
    __extension__ using Wide=unsigned __int128;
    const Wide dividend=(static_cast<Wide>(high)<<64u)|
                        static_cast<Wide>(low);
    return {static_cast<std::uint64_t>(dividend/divisor),
            static_cast<std::uint64_t>(dividend%divisor)};
#else
    return divide_wide_portable(high,low,divisor);
#endif
}

} // namespace

BigUnsigned::Limb BigUnsigned::divrem_limb_into(
        LimbView a,Limb divisor,MutableLimbView quotient)
{
    if (divisor==0u) {
        throw std::domain_error("BigUnsigned division by zero");
    }
    if (quotient.size()<a.size()) {
        throw std::length_error("BigUnsigned quotient output too small");
    }

    Limb remainder=0u;
    for (std::size_t i=a.size();i!=0u;--i) {
        const auto [digit,next]=divide_wide(
            remainder,a[i-1u],divisor);
        quotient[i-1u]=digit;
        remainder=next;
    }
    std::fill(quotient.subspan(a.size()).begin(),quotient.end(),Limb{0});
    return remainder;
}

std::pair<BigUnsigned,BigUnsigned::Limb> divrem_limb(
        BigUnsigned const& a,BigUnsigned::Limb divisor)
{
    BigUnsigned::Scratch scratch;
    auto digits=scratch.acquire(a.size());
    const auto remainder=BigUnsigned::divrem_limb_into(
        a.limbs(),divisor,digits);
    return {BigUnsigned(digits),remainder};
}


std::size_t BigUnsigned::divrem_scratch_limbs(
        std::size_t dividend_limbs,std::size_t divisor_limbs)
{
    if (divisor_limbs<2u || dividend_limbs<divisor_limbs) {
        return 0u;
    }
    const auto max=std::numeric_limits<std::size_t>::max();
    if (divisor_limbs>max-1u ||
            dividend_limbs>max-1u-divisor_limbs) {
        throw std::length_error("BigUnsigned division scratch size overflow");
    }
    // B27: the fused subtraction eliminates the digit-product buffer.
    // Only the normalized dividend (n+1) and divisor (m) remain.
    return dividend_limbs+divisor_limbs+1u;
}

namespace {

// Portable 64x64 -> 128 fallback for the fused division step. The
// same algebra as the reference multiply kernel: no wide integer or
// optional CPU instruction is required.
[[maybe_unused]] std::pair<BigUnsigned::Limb,BigUnsigned::Limb>
multiply_digit_portable(BigUnsigned::Limb x,BigUnsigned::Limb y) noexcept
{
    using Limb=BigUnsigned::Limb;
    constexpr Limb mask=0xffffffffULL;
    const Limb x0=x&mask,x1=x>>32u,y0=y&mask,y1=y>>32u;
    const Limb p00=x0*y0,p01=x0*y1,p10=x1*y0,p11=x1*y1;
    const Limb middle=(p00>>32u)+(p01&mask)+(p10&mask);
    return {(p00&mask)|(middle<<32u),
            p11+(p01>>32u)+(p10>>32u)+(middle>>32u)};
}

// Subtract digit*divisor from the (m+1)-limb partial dividend in
// place, fusing multiplication and subtraction into one pass.
// Each carry is the high product limb plus the subtraction borrow.
// The wide product+carry is < 2^128, so carry always fits 64 bits.
// Return whether the signed mathematical difference is negative.
template<bool SparseZero>
bool subtract_digit_product(BigUnsigned::MutableLimbView window,
                            BigUnsigned::LimbView divisor,
                            BigUnsigned::Limb digit) noexcept
{
    using Limb=BigUnsigned::Limb;
    Limb carry=0u;
    for (std::size_t i=0u;i<divisor.size();++i) {
        if constexpr (SparseZero) {
            // If this normalized divisor limb is zero, the high
            // product vanishes. Only the incoming carry/borrow must
            // be subtracted; the dense kernel stays unchanged.
            if (divisor[i]==0u) {
                const Limb old=window[i];
                window[i]=old-carry;
                carry=static_cast<Limb>(old<carry);
                continue;
            }
        }
        Limb low,high;
#if defined(_MSC_VER) && defined(_M_X64) && !defined(__clang__)
        high=0u;
        low=_umul128(divisor[i],digit,&high);
        const Limb previous=low;
        low+=carry;
        high+=static_cast<Limb>(low<previous);
#elif defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
        __extension__ using Wide=unsigned __int128;
        const Wide product=static_cast<Wide>(divisor[i])*digit+carry;
        low=static_cast<Limb>(product);
        high=static_cast<Limb>(product>>64u);
#else
        const auto [product_low,product_high]=
            multiply_digit_portable(divisor[i],digit);
        low=product_low+carry;
        high=product_high+static_cast<Limb>(low<product_low);
#endif
        const Limb old=window[i];
        window[i]=old-low;
        carry=high+static_cast<Limb>(old<low);
    }
    const std::size_t m=divisor.size();
    const Limb old=window[m];
    window[m]=old-carry;
    return old<carry;
}

// Undo one unit of an overestimated quotient digit. When the
// subtraction was negative, the working window holds its two's
// complement modulo B^(m+1). After adding the divisor back, a
// nonzero high limb means the estimate is still too large.
// Knuth's normalized estimate needs at most two such add-backs.
void add_back_divisor(BigUnsigned::MutableLimbView window,
                      BigUnsigned::LimbView divisor) noexcept
{
    using Limb=BigUnsigned::Limb;
    Limb carry=0u;
    for (std::size_t i=0u;i<divisor.size();++i) {
        const Limb old=window[i];
        const Limb sum=old+divisor[i];
        const Limb result=sum+carry;
        carry=static_cast<Limb>((sum<old)||(result<sum));
        window[i]=result;
    }
    window[divisor.size()]+=carry;
}

} // namespace

template<bool WantQuotient,bool WantRemainder>
void BigUnsigned::divrem_selected_into(
        LimbView a,LimbView b,MutableLimbView quotient,
        MutableLimbView remainder,MutableLimbView scratch)
{
    static_assert(WantQuotient || WantRemainder);
    if (b.empty()) {
        throw std::domain_error("BigUnsigned division by zero");
    }
    const std::size_t quotient_limbs=
        a.size()<b.size() ? 0u : a.size()-b.size()+1u;
    const std::size_t remainder_limbs=std::min(a.size(),b.size());
    if constexpr (WantQuotient) {
        if (quotient.size()<quotient_limbs) {
            throw std::length_error("BigUnsigned quotient output too small");
        }
    }
    if constexpr (WantRemainder) {
        if (remainder.size()<remainder_limbs) {
            throw std::length_error("BigUnsigned remainder output too small");
        }
    }
    const auto required=divrem_scratch_limbs(a.size(),b.size());
    if (scratch.size()<required) {
        throw std::length_error("BigUnsigned division scratch too small");
    }

    if constexpr (WantQuotient) {
        std::fill(quotient.begin(),quotient.end(),Limb{0});
    }
    if constexpr (WantRemainder) {
        std::fill(remainder.begin(),remainder.end(),Limb{0});
    }
    if (a.size()<b.size()) {
        if constexpr (WantRemainder) {
            std::copy(a.begin(),a.end(),remainder.begin());
        }
        return;
    }
    if (b.size()==1u) {
        if constexpr (WantQuotient) {
            const Limb low=divrem_limb_into(a,b[0],quotient);
            if constexpr (WantRemainder) {
                remainder[0u]=low;
            }
        } else {
            // No quotient buffer or object is allocated for rem-only.
            Limb low=0u;
            for (std::size_t i=a.size();i!=0u;--i) {
                low=divide_wide(low,a[i-1u],b[0u]).second;
            }
            remainder[0u]=low;
        }
        return;
    }

    const std::size_t n=a.size();
    const std::size_t m=b.size();

    // B28: equal leading limbs imply 0 <= a/b < 2 for equal-length
    // normalized magnitudes. Determine the exact quotient (0 or 1)
    // without normalizing or estimating a quotient digit.
    // Preserve the checked public scratch requirement even on this path.
    if (n==m && a.back()==b.back()) {
        if constexpr (WantRemainder) {
            // The final borrow determines whether a<b. When a>=b,
            // this one subtraction already yields the exact remainder.
            Limb borrow=0u;
            for (std::size_t i=0u;i<m;++i) {
                const Limb difference=a[i]-b[i];
                remainder[i]=difference-borrow;
                borrow=static_cast<Limb>((a[i]<b[i]) ||
                                         (difference<borrow));
            }
            if (borrow!=0u) {
                // The quotient is zero, so the remainder is a.
                std::copy(a.begin(),a.end(),remainder.begin());
            }
            if constexpr (WantQuotient) {
                quotient[0u]=static_cast<Limb>(borrow==0u);
            }
        } else if constexpr (WantQuotient) {
            // Quotient-only does not need to construct a remainder.
            for (std::size_t i=m;i!=0u;--i) {
                if (a[i-1u]!=b[i-1u]) {
                    quotient[0u]=static_cast<Limb>(a[i-1u]>b[i-1u]);
                    return;
                }
            }
            quotient[0u]=1u;
        }
        return;
    }

    auto u=scratch.first(n+1u);
    auto v=scratch.subspan(n+1u,m);

    // All three compile-time output modes share the same normalized
    // quotient estimation, in-place subtraction and correction.
    const unsigned shift=static_cast<unsigned>(std::countl_zero(b.back()));
    Limb carry=0u;
    for (std::size_t i=0u;i<n;++i) {
        const auto digit=a[i];
        u[i]=(digit<<shift)|carry;
        carry=shift==0u?0u:digit>>(64u-shift);
    }
    u[n]=carry;
    carry=0u;
    for (std::size_t i=0u;i<m;++i) {
        const auto digit=b[i];
        v[i]=(digit<<shift)|carry;
        carry=shift==0u?0u:digit>>(64u-shift);
    }

    // B29: only sufficiently large, sparse normalized divisors
    // justify testing individual zero limbs. The selection is made
    // once per division, and the dense digit loop has no extra check.
    const std::size_t quotient_digits=n-m+1u;
    const bool use_sparse= m>=8u && quotient_digits>=8u &&
        static_cast<std::size_t>(
            std::count(v.begin(),v.end(),Limb{0})) >= m/2u+m%2u;

    const auto divide_digits=[&]<bool SparseZero>() {
        for (std::size_t counter=quotient_digits;counter!=0u;--counter) {
            const std::size_t j=counter-1u;
            const Limb high=u[j+m];
            const Limb low=u[j+m-1u];
            Limb digit=high>=v[m-1u] ?
                std::numeric_limits<Limb>::max() :
                divide_wide(high,low,v[m-1u]).first;
            const auto window=u.subspan(j,m+1u);
            bool overestimated=
                subtract_digit_product<SparseZero>(window,v,digit);
            while (overestimated) {
                --digit;
                add_back_divisor(window,v);
                overestimated=window.back()!=0u;
            }
            if constexpr (WantQuotient) {
                quotient[j]=digit;
            }
        }
    };
    if (use_sparse) {
        divide_digits.template operator()<true>();
    } else {
        divide_digits.template operator()<false>();
    }

    if constexpr (WantRemainder) {
        // The final partial dividend already holds the remainder.
        for (std::size_t i=0u;i<m;++i) {
            const Limb high=i+1u<m?u[i+1u]:Limb{0};
            remainder[i]=shift==0u?u[i]:
                (u[i]>>shift)|(high<<(64u-shift));
        }
    }
}

void BigUnsigned::divrem_into(
        LimbView a,LimbView b,MutableLimbView quotient,
        MutableLimbView remainder,MutableLimbView scratch)
{
    // Keep the historical public output sizes and fail-before-mutation
    // checks, even though specialized internal outputs can be smaller.
    if (b.empty()) {
        throw std::domain_error("BigUnsigned division by zero");
    }
    if (quotient.size()<a.size() ||
            remainder.size()<std::min(a.size(),b.size())) {
        throw std::length_error("BigUnsigned division output too small");
    }
    divrem_selected_into<true,true>(a,b,quotient,remainder,scratch);
}

std::pair<BigUnsigned,BigUnsigned> divrem(
        BigUnsigned const& a,BigUnsigned const& b)
{
    if (b.is_zero()) {
        throw std::domain_error("BigUnsigned division by zero");
    }
    BigUnsigned quotient;
    BigUnsigned remainder;
    quotient.resize_for_result(a.size());
    const auto rem_size=std::min(a.size(),b.size());
    remainder.resize_for_result(rem_size);
    BigUnsigned::Scratch workspace;
    const auto scratch=workspace.acquire(
        BigUnsigned::divrem_scratch_limbs(a.size(),b.size()));
    BigUnsigned::divrem_into(
        a.limbs(),b.limbs(),
        BigUnsigned::MutableLimbView(quotient.mutable_data(),a.size()),
        BigUnsigned::MutableLimbView(remainder.mutable_data(),rem_size),
        scratch);
    quotient.normalize();
    remainder.normalize();
    return {std::move(quotient),std::move(remainder)};
}

BigUnsigned quot(BigUnsigned const& a,BigUnsigned const& b)
{
    if (b.is_zero()) {
        throw std::domain_error("BigUnsigned division by zero");
    }
    const auto count=a.size()<b.size()?0u:a.size()-b.size()+1u;
    BigUnsigned quotient;
    quotient.resize_for_result(count);
    BigUnsigned::Scratch workspace;
    const auto scratch=workspace.acquire(
        BigUnsigned::divrem_scratch_limbs(a.size(),b.size()));
    BigUnsigned::divrem_selected_into<true,false>(
        a.limbs(),b.limbs(),
        BigUnsigned::MutableLimbView(quotient.mutable_data(),count),
        {},scratch);
    quotient.normalize();
    return quotient;
}

BigUnsigned rem(BigUnsigned const& a,BigUnsigned const& b)
{
    if (b.is_zero()) {
        throw std::domain_error("BigUnsigned division by zero");
    }
    const auto count=std::min(a.size(),b.size());
    BigUnsigned remainder;
    remainder.resize_for_result(count);
    BigUnsigned::Scratch workspace;
    const auto scratch=workspace.acquire(
        BigUnsigned::divrem_scratch_limbs(a.size(),b.size()));
    BigUnsigned::divrem_selected_into<false,true>(
        a.limbs(),b.limbs(),{},
        BigUnsigned::MutableLimbView(remainder.mutable_data(),count),
        scratch);
    remainder.normalize();
    return remainder;
}

BigUnsigned gcd(BigUnsigned a,BigUnsigned b)
{
    while (!b.is_zero()) {
        auto division=divrem(a,b);
        a=std::move(b);
        b=std::move(division.second);
    }
    return a;
}

BigUnsigned from_decimal(std::string_view input)
{
    if (input.empty()) {
        throw std::invalid_argument("BigUnsigned empty decimal input");
    }
    for (const char digit:input) {
        if (digit<'0' || digit>'9') {
            throw std::invalid_argument("BigUnsigned invalid decimal digit");
        }
    }

    // B24: for long strings, group 18 decimal digits per pass and
    // accumulate an exact 64x60-bit product in an unsigned 128-bit
    // intermediate. This halves the number of full-limb scans. The
    // carry is < 10^18, so (2^64-1)*10^18+carry fits in 128 bits.
    // Keep B19's portable radix-10^9 implementation below unchanged
    // for <=309 digits and for toolchains without native 128-bit ints.
#if defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    if (input.size()>309u) {
        __extension__ using Wide=unsigned __int128;
        constexpr std::size_t group=18u;
        std::vector<BigUnsigned::Limb> digits;
        std::size_t pos=0u;
        while (pos<input.size()) {
            const std::size_t count=std::min(group,input.size()-pos);
            BigUnsigned::Limb radix=1u;
            BigUnsigned::Limb chunk=0u;
            for (std::size_t i=0u;i<count;++i) {
                radix*=10u;
                chunk=chunk*10u+
                    static_cast<BigUnsigned::Limb>(input[pos+i]-'0');
            }
            BigUnsigned::Limb carry=chunk;
            for (auto& digit:digits) {
                const Wide product=static_cast<Wide>(digit)*radix+carry;
                digit=static_cast<BigUnsigned::Limb>(product);
                carry=static_cast<BigUnsigned::Limb>(product>>64u);
            }
            if (carry!=0u) {
                digits.push_back(carry);
            }
            pos+=count;
        }
        return BigUnsigned(digits);
    }
#endif

    // Multiply in base-10^9 blocks using 32-bit half-products, never
    // assuming compiler support for a 128-bit integer.
    constexpr std::size_t group=9u;
    constexpr BigUnsigned::Limb low_mask=0xffffffffULL;
    std::vector<BigUnsigned::Limb> digits;
    std::size_t pos=0u;
    while (pos<input.size()) {
        const std::size_t count=std::min(group,input.size()-pos);
        BigUnsigned::Limb radix=1u;
        BigUnsigned::Limb chunk=0u;
        for (std::size_t i=0u;i<count;++i) {
            radix*=10u;
            chunk=chunk*10u+
                static_cast<BigUnsigned::Limb>(input[pos+i]-'0');
        }
        BigUnsigned::Limb carry=chunk;
        for (auto& digit:digits) {
            const BigUnsigned::Limb low=(digit&low_mask)*radix+carry;
            const BigUnsigned::Limb high=
                (digit>>32u)*radix+(low>>32u);
            digit=(high<<32u)|(low&low_mask);
            carry=high>>32u;
        }
        if (carry!=0u) {
            digits.push_back(carry);
        }
        pos+=count;
    }
    return BigUnsigned(digits);
}

namespace {

// B25.1: enable wide-radix decimal output only on native ARM64.
// Its performance is better on ARM64 but substantially worse on x86-64.
// Select using the compiler-provided target macro, not a CMake option
// or a build-specific definition; leave B19's radix-10^9 fallback
// unchanged for x86-64, MSVC and other targets. The 128-bit division
// processes 18 decimal digits per quotient/remainder step on ARM64.
#if defined(__aarch64__) && defined(__SIZEOF_INT128__) && \
    (defined(__GNUC__) || defined(__clang__))
std::string to_decimal_blocks_wide(BigUnsigned const& value)
{
    __extension__ using Wide=unsigned __int128;
    constexpr BigUnsigned::Limb radix=1000000000000000000ULL;
    constexpr std::size_t block_width=18u;

    const auto limbs=value.limbs();
    std::vector<BigUnsigned::Limb> work(limbs.begin(),limbs.end());
    std::vector<BigUnsigned::Limb> blocks;
    while (!work.empty()) {
        BigUnsigned::Limb remainder=0u;
        for (std::size_t i=work.size();i!=0u;--i) {
            // remainder < 10^18 and digit < 2^64, so the
            // numerator is < 10^18 * 2^64 and fits in Wide.
            const Wide numerator=(static_cast<Wide>(remainder)<<64u)|
                                  static_cast<Wide>(work[i-1u]);
            work[i-1u]=static_cast<BigUnsigned::Limb>(numerator/radix);
            remainder=static_cast<BigUnsigned::Limb>(numerator%radix);
        }
        blocks.push_back(remainder);
        while (!work.empty() && work.back()==0u) {
            work.pop_back();
        }
    }

    std::string result=std::to_string(blocks.back());
    for (std::size_t i=blocks.size()-1u;i!=0u;--i) {
        const std::string block=std::to_string(blocks[i-1u]);
        result.append(block_width-block.size(),'0');
        result+=block;
    }
    return result;
}
#endif

// B19's quick base-10^9 division remains the small-value and
// portable fallback leaf algorithm.
std::string to_decimal_blocks(BigUnsigned const& value)
{
    if (value.is_zero()) {
        return "0";
    }
#if defined(__aarch64__) && defined(__SIZEOF_INT128__) && \
    (defined(__GNUC__) || defined(__clang__))
    // Retain the original path through the 8-limb boundary;
    // the 18-digit divider is used for >=9 limbs (>=513 bits).
    if (value.size()>8u) {
        return to_decimal_blocks_wide(value);
    }
#endif

    constexpr BigUnsigned::Limb radix=1000000000u;
    constexpr std::size_t block_width=9u;
    const auto limbs=value.limbs();
    std::vector<BigUnsigned::Limb> work(limbs.begin(),limbs.end());
    std::vector<unsigned> blocks;
    // B31: avoid repeated growth of the decimal-block container.
    // A 64-bit limb contributes at most three base-10^9 blocks.
    if (limbs.size()<=blocks.max_size()/3u) {
        blocks.reserve(limbs.size()*3u);
    }
    while (!work.empty()) {
        // radix < 2^32. Two 64/32 divisions avoid a general 128/64
        // division, and the per-limb output update is in-place.
        std::uint64_t remainder=0u;
        for (std::size_t i=work.size();i!=0u;--i) {
            const BigUnsigned::Limb digit=work[i-1u];
            const std::uint64_t upper=(remainder<<32u)|(digit>>32u);
            const std::uint64_t quotient_high=upper/radix;
            remainder=upper%radix;
            const std::uint64_t lower=(remainder<<32u)|
                                      (digit&BigUnsigned::Limb{0xffffffffu});
            const std::uint64_t quotient_low=lower/radix;
            remainder=lower%radix;
            work[i-1u]=(quotient_high<<32u)|quotient_low;
        }
        blocks.push_back(static_cast<unsigned>(remainder));
        while (!work.empty() && work.back()==0u) {
            work.pop_back();
        }
    }

    // B31: format the base-10^9 blocks into one exactly sized
    // string. Writing from the end preserves leading zeros in each
    // nine-digit block without temporary strings or repeated appends.
    static constexpr auto digit_pairs=[] {
        std::array<char,200> table{};
        for (unsigned pair_value=0u;pair_value<100u;++pair_value) {
            table[2u*pair_value]=static_cast<char>('0'+pair_value/10u);
            table[2u*pair_value+1u]=static_cast<char>('0'+pair_value%10u);
        }
        return table;
    }();

    std::array<char,10> top_buffer{};
    const auto [top_end,error]=std::to_chars(
        top_buffer.data(),top_buffer.data()+top_buffer.size(),
        blocks.back());
    if (error!=std::errc{}) {
        throw std::length_error("BigUnsigned decimal head overflow");
    }
    const std::size_t top_width=static_cast<std::size_t>(
        top_end-top_buffer.data());
    std::string result;
    if (blocks.size()-1u>(result.max_size()-top_width)/block_width) {
        throw std::length_error("BigUnsigned decimal output overflow");
    }
    result.resize(top_width+(blocks.size()-1u)*block_width);
    std::copy_n(top_buffer.data(),top_width,result.data());

    char* output=result.data()+result.size();
    for (std::size_t i=0u;i+1u<blocks.size();++i) {
        unsigned block=blocks[i]; // Each block is strictly below 10^9.
        for (unsigned pair_index=0u;pair_index<4u;++pair_index) {
            const unsigned pair=block%100u;
            block/=100u;
            output-=2u;
            output[0u]=digit_pairs[2u*pair];
            output[1u]=digit_pairs[2u*pair+1u];
        }
        *--output=static_cast<char>('0'+block);
    }
    return result;
}

// Each power equals 10^(9*2^i); width is its exact number of
// decimal digits after the split. Computed once per conversion.
struct DecimalSplitPower {
    BigUnsigned value;
    std::size_t width;
};

constexpr std::size_t decimal_recursive_cutoff_limbs=16u;

std::string to_decimal_recursive(
        BigUnsigned const& value,
        std::vector<DecimalSplitPower> const& powers)
{
    if (value.size()<=decimal_recursive_cutoff_limbs) {
        return to_decimal_blocks(value);
    }

    // Take the largest 10^(9*2^i) no wider than half the input's
    // bit length; this yields a near-balanced split without a
    // floating-point estimate of the number of decimal digits.
    const std::size_t midpoint=value.bit_length()/2u;
    std::size_t count=powers.size();
    while (count!=0u && powers[count-1u].value.bit_length()>midpoint) {
        --count;
    }
    if (count==0u) {
        return to_decimal_blocks(value);
    }

    auto const& divisor=powers[count-1u];
    auto [high,low]=divrem(value,divisor.value);
    std::string result=to_decimal_recursive(high,powers);
    const std::string suffix=to_decimal_recursive(low,powers);
    // Since 0 <= low < 10^width, suffix.size() <= width.
    result.append(divisor.width-suffix.size(),'0');
    result+=suffix;
    return result;
}

} // namespace

std::string to_decimal(BigUnsigned const& value)
{
    if (value.is_zero()) {
        return "0";
    }
    if (value.size()<=decimal_recursive_cutoff_limbs) {
        return to_decimal_blocks(value);
    }

    // Squaring successive powers of 10^9 costs logarithmically many
    // products. Each recursive step uses the shared exact divrem;
    // no GMP boundary conversions, global caches, or ISA assumptions.
    std::vector<DecimalSplitPower> powers;
    powers.push_back({BigUnsigned(std::uint64_t{1000000000u}),9u});
    const std::size_t midpoint=value.bit_length()/2u;
    while (true) {
        const auto& last=powers.back();
        const std::size_t bits=last.value.bit_length();
        // The square has at least 2*bits-1 bits; avoid the large
        // unnecessary final square when it cannot fit the midpoint.
        if (bits>(midpoint+1u)/2u) {
            break;
        }
        if (last.width>std::numeric_limits<std::size_t>::max()/2u) {
            throw std::length_error("BigUnsigned decimal width overflow");
        }
        BigUnsigned squared=square(last.value);
        if (squared.bit_length()>midpoint) {
            break;
        }
        const std::size_t width=last.width*2u;
        powers.push_back({std::move(squared),width});
    }
    return to_decimal_recursive(value,powers);
}

namespace {

// Compiler-supported wide accumulation, shared by GCC/Clang x86-64 and ARM64.
// The intermediate is bounded by (B-1)^2+2*(B-1)=B^2-1, B=2^64.
// Preserve the historical MSVC x64 native B3 fallback. On other toolchains
// without unsigned __int128, retain portable multiplication as a fallback.
// Inputs and output must be disjoint; all output limbs are overwritten.
void multiply_limbs_fused(BigUnsigned::LimbView a,
                          BigUnsigned::LimbView b,
                          BigUnsigned::MutableLimbView output) noexcept
{
#if defined(_MSC_VER) && defined(_M_X64) && defined(ARIADNE_BIG_UNSIGNED_X86_64)
    BigUnsignedDetail::mul_low_into_x86_64(a,b,output);
#elif defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    __extension__ using Wide = unsigned __int128;
    std::fill(output.begin(),output.end(),BigUnsigned::Limb{0});
    for (std::size_t i=0u; i<a.size() && i<output.size(); ++i) {
        BigUnsigned::Limb carry=0u;
        const std::size_t terms=std::min(b.size(),output.size()-i);
        for (std::size_t j=0u; j<terms; ++j) {
            const Wide total=static_cast<Wide>(a[i])*b[j]
                             +output[i+j]+carry;
            output[i+j]=static_cast<BigUnsigned::Limb>(total);
            carry=static_cast<BigUnsigned::Limb>(total>>64u);
        }
        if (i+terms<output.size()) {
            output[i+terms]=carry;
        }
    }
#else
    BigUnsigned::mul_low_into_portable(a,b,output);
#endif
}

} // namespace

bool BigUnsigned::has_native_multiply() noexcept {
#if defined(ARIADNE_BIG_UNSIGNED_X86_64) || defined(ARIADNE_BIG_UNSIGNED_ARM64)
    return true;
#else
    return false;
#endif
}

bool BigUnsigned::has_arm64_multiply() noexcept {
#if defined(ARIADNE_BIG_UNSIGNED_ARM64)
    return true;
#else
    return false;
#endif
}

bool BigUnsigned::has_x86_64_carry_multiply() noexcept {
#if defined(ARIADNE_BIG_UNSIGNED_X86_64) && defined(__x86_64__) \
    && (defined(__GNUC__) || defined(__clang__)) && !defined(_MSC_VER)
    return true;
#else
    return false;
#endif
}

bool BigUnsigned::has_x86_64_multiply() noexcept {
#if defined(ARIADNE_BIG_UNSIGNED_X86_64)
    return true;
#else
    return false;
#endif
}

bool BigUnsigned::selected_multiply_uses_fused() noexcept {
#if defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED)
    return true;
#else
    return false;
#endif
}

void BigUnsigned::mul_low_into_comba(
        LimbView a,LimbView b,MutableLimbView output)
{
#if defined(ARIADNE_BIG_UNSIGNED_X86_64)
    BigUnsignedDetail::mul_low_into_comba_x86_64(a,b,output);
#else
    mul_low_into_portable(a,b,output);
#endif
}

void BigUnsigned::mul_low_into_fused(
        LimbView a,LimbView b,MutableLimbView output)
{
#if defined(ARIADNE_BIG_UNSIGNED_X86_64) || defined(ARIADNE_BIG_UNSIGNED_ARM64)
    multiply_limbs_fused(a,b,output);
#else
    mul_low_into_portable(a,b,output);
#endif
}

void BigUnsigned::mul_low_into_arm64_carry(
        LimbView a,LimbView b,MutableLimbView output)
{
#if defined(ARIADNE_BIG_UNSIGNED_ARM64)
    BigUnsignedDetail::mul_low_into_carry_arm64(a,b,output);
#else
    mul_low_into_portable(a,b,output);
#endif
}

void BigUnsigned::mul_low_into_x86_64_carry(
        LimbView a,LimbView b,MutableLimbView output)
{
#if defined(ARIADNE_BIG_UNSIGNED_X86_64)
    BigUnsignedDetail::mul_low_into_carry_x86_64(a,b,output);
#else
    mul_low_into_portable(a,b,output);
#endif
}

void BigUnsigned::mul_low_into_schoolbook(
        LimbView a,LimbView b,MutableLimbView output)
{
#if defined(ARIADNE_BIG_UNSIGNED_X86_64)
    // Reproduce the B3 selected behavior, including its one-limb
    // portable fallback, to preserve a strict A/B control.
    if (a.size() >= 2u && b.size() >= 2u) {
        BigUnsignedDetail::mul_low_into_x86_64(a,b,output);
        return;
    }
#endif
    mul_low_into_portable(a,b,output);
}

void BigUnsigned::mul_low_into(LimbView a,LimbView b,MutableLimbView output) {
#if defined(ARIADNE_BIG_UNSIGNED_ARM64)
    // Apple Silicon measurements and native GMP differential tests validate
    // the explicit carry-chain path; fused remains available for comparison.
    BigUnsignedDetail::mul_low_into_carry_arm64(a,b,output);
#elif defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED)
    // GCC x86-64 continues to select the shared unsigned __int128 kernel.
    multiply_limbs_fused(a,b,output);
#else
    mul_low_into_schoolbook(a,b,output);
#endif
}

namespace {

void check_full_product_size(
        BigUnsigned::LimbView a,BigUnsigned::LimbView b,
        BigUnsigned::MutableLimbView output)
{
    if (a.empty() || b.empty()) {
        return;
    }
    if (a.size() > std::numeric_limits<std::size_t>::max()-b.size()
            || output.size() < a.size()+b.size()) {
        throw std::length_error("BigUnsigned mul_into output too small");
    }
}

} // namespace

void BigUnsigned::mul_into(LimbView a,LimbView b,MutableLimbView output) {
    check_full_product_size(a,b,output);
    mul_low_into(a,b,output);
}

void BigUnsigned::mul_into_arm64_carry(
        LimbView a,LimbView b,MutableLimbView output)
{
    check_full_product_size(a,b,output);
    mul_low_into_arm64_carry(a,b,output);
}

void BigUnsigned::mul_into_x86_64_carry(
        LimbView a,LimbView b,MutableLimbView output)
{
    check_full_product_size(a,b,output);
    mul_low_into_x86_64_carry(a,b,output);
}

void BigUnsigned::mul_into_portable(
        LimbView a,LimbView b,MutableLimbView output)
{
    check_full_product_size(a,b,output);
    mul_low_into_portable(a,b,output);
}

void BigUnsigned::mul_into_comba(
        LimbView a,LimbView b,MutableLimbView output)
{
    check_full_product_size(a,b,output);
    mul_low_into_comba(a,b,output);
}

void BigUnsigned::mul_into_fused(
        LimbView a,LimbView b,MutableLimbView output)
{
    check_full_product_size(a,b,output);
    mul_low_into_fused(a,b,output);
}

void BigUnsigned::mul_into_schoolbook(
        LimbView a,LimbView b,MutableLimbView output)
{
    check_full_product_size(a,b,output);
    mul_low_into_schoolbook(a,b,output);
}

BigUnsigned mul_low(BigUnsigned const& a,BigUnsigned const& b,
                    std::size_t limb_count)
{
    if (a.is_zero() || b.is_zero() || limb_count==0u) {
        return BigUnsigned{};
    }

    if (a._size > std::numeric_limits<std::size_t>::max()-b._size) {
        throw std::length_error("BigUnsigned product size overflow");
    }

    const std::size_t result_size = std::min(limb_count,a._size+b._size);
    BigUnsigned result;
    result.resize_for_result(result_size);
    BigUnsigned::mul_low_into(
        a.limbs(),b.limbs(),
        BigUnsigned::MutableLimbView(result.mutable_data(),result_size));
    result.normalize();
    return result;
}

BigUnsigned mul(BigUnsigned const& a,BigUnsigned const& b)
{
    if (a.is_zero() || b.is_zero()) {
        return BigUnsigned{};
    }

    if (a._size > std::numeric_limits<std::size_t>::max()-b._size) {
        throw std::length_error("BigUnsigned product size overflow");
    }

    BigUnsigned result;
    result.resize_for_result(a._size+b._size);
    BigUnsigned::mul_into(
        a.limbs(),b.limbs(),
        BigUnsigned::MutableLimbView(result.mutable_data(),a._size+b._size));
    result.normalize();
    return result;
}

// B32: x86-64 GCC's selected fused schoolbook kernel evaluates each
// off-diagonal product twice in a*a. For >=8 limbs a symmetric square
// requires each a[i]*a[j] (i<j) only once. The existing output buffer
// first collects these cross-products, doubles them, and finally adds
// the diagonal squares. The two carry passes are exact in radix 2^64.
// Keep the old selected kernel on ARM64, Clang, MSVC and non-fused
// GCC builds until their own native benchmarks justify dispatch.
#if defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED) && \
    defined(__x86_64__) && defined(__SIZEOF_INT128__) && \
    defined(__GNUC__) && !defined(__clang__) && !defined(_MSC_VER)
namespace {
void square_limbs_symmetric(BigUnsigned::LimbView a,
                            BigUnsigned::MutableLimbView output) noexcept
{
    using Limb=BigUnsigned::Limb;
    __extension__ using Wide=unsigned __int128;
    const std::size_t n=a.size();
    std::fill(output.begin(),output.end(),Limb{0});

    // Accumulate sum_{i<j} a[i]*a[j]*B^(i+j) with one wide product
    // per unordered pair. Slot i+n has not been touched by prior rows.
    for (std::size_t i=0u;i+1u<n;++i) {
        Limb carry=0u;
        for (std::size_t j=i+1u;j<n;++j) {
            const std::size_t k=i+j;
            const Wide sum=static_cast<Wide>(a[i])*a[j]+
                           output[k]+carry;
            output[k]=static_cast<Limb>(sum);
            carry=static_cast<Limb>(sum>>64u);
        }
        output[i+n]=carry;
    }

    // Every cross term occurs twice in a^2. The doubled cross sum
    // fits in 2*n limbs because it is <= a^2 < B^(2*n).
    Limb carry=0u;
    for (std::size_t k=0u;k<output.size();++k) {
        const Limb word=output[k];
        output[k]=(word<<1u)|carry;
        carry=word>>63u;
    }

    // Add the diagonal squares in ascending digit order, preserving
    // the carry between consecutive two-limb diagonal products.
    carry=0u;
    for (std::size_t i=0u;i<n;++i) {
        const Wide diagonal=static_cast<Wide>(a[i])*a[i];
        const std::size_t k=2u*i;
        const Wide low=static_cast<Wide>(output[k])+
                       static_cast<Limb>(diagonal)+carry;
        output[k]=static_cast<Limb>(low);
        carry=static_cast<Limb>(low>>64u);
        const Wide high=static_cast<Wide>(output[k+1u])+
                        static_cast<Limb>(diagonal>>64u)+carry;
        output[k+1u]=static_cast<Limb>(high);
        carry=static_cast<Limb>(high>>64u);
    }
}
} // namespace
#endif

BigUnsigned square(BigUnsigned const& a)
{
#if defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED) && \
    defined(__x86_64__) && defined(__SIZEOF_INT128__) && \
    defined(__GNUC__) && !defined(__clang__) && !defined(_MSC_VER)
    if (a._size>=8u) {
        if (a._size>std::numeric_limits<std::size_t>::max()/2u) {
            throw std::length_error("BigUnsigned square size overflow");
        }
        BigUnsigned result;
        result.resize_for_result(2u*a._size);
        square_limbs_symmetric(
            a.limbs(),BigUnsigned::MutableLimbView(
                result.mutable_data(),2u*a._size));
        result.normalize();
        return result;
    }
#endif
    return mul(a,a);
}


// Portable, allocation-free Karatsuba implementation (formerly separate TU).
// Kept in the common arithmetic source alongside selected limb kernels.
namespace {

using Limb = BigUnsigned::Limb;
using Read = BigUnsigned::LimbView;
using Write = BigUnsigned::MutableLimbView;

// This is a recursion cutoff, not a proven performance crossover.
// Small cases use the current architecture-selected schoolbook kernel.
constexpr std::size_t karatsuba_default_cutoff=8u;
constexpr std::size_t karatsuba_minimum_cutoff=8u;

std::size_t checked_sum(std::size_t a,std::size_t b)
{
    if (a>std::numeric_limits<std::size_t>::max()-b) {
        throw std::length_error("BigUnsigned Karatsuba size overflow");
    }
    return a+b;
}

bool should_split(std::size_t a,std::size_t b,
                  std::size_t cutoff) noexcept
{
    const std::size_t half=std::max(a,b)/2u;
    return std::min(a,b)>cutoff
        && a>half && b>half;
}

// Exact recursion-frame bound: sumA, sumB, cross product, plus the
// maximum of the three sequential recursive calls' workspaces.
// No nested frame overlaps the live parent sums or cross product.
std::size_t karatsuba_workspace(
        std::size_t a,std::size_t b,std::size_t cutoff)
{
    if (!should_split(a,b,cutoff)) {
        return 0u;
    }

    const std::size_t half=std::max(a,b)/2u;
    const std::size_t a0=std::min(a,half);
    const std::size_t b0=std::min(b,half);
    const std::size_t a1=a-a0;
    const std::size_t b1=b-b0;
    const std::size_t sa=checked_sum(std::max(a0,a1),1u);
    const std::size_t sb=checked_sum(std::max(b0,b1),1u);
    const std::size_t pair=checked_sum(sa,sb);
    const std::size_t frame=checked_sum(pair,pair);

    const std::size_t nested=std::max({
        karatsuba_workspace(a0,b0,cutoff),
        karatsuba_workspace(a1,b1,cutoff),
        karatsuba_workspace(sa,sb,cutoff)
    });
    return checked_sum(frame,nested);
}

// An extra high zero limb makes addition exact even when both halves
// are saturated. Arguments must not overlap the destination.
void sum_halves(Read a,Read b,Write output) noexcept
{
    Limb carry=0u;
    for (std::size_t i=0u;i<output.size();++i) {
        const Limb left=i<a.size() ? a[i] : Limb{0};
        const Limb right=i<b.size() ? b[i] : Limb{0};
        const Limb sum=left+right;
        const Limb value=sum+carry;
        carry=static_cast<Limb>(sum<left || value<sum);
        output[i]=value;
    }
}

// The caller has proved lhs >= rhs as unsigned integers. Subtracting
// from lhs itself avoids an extra allocated intermediate buffer.
void subtract_in_place(Write lhs,Read rhs)
{
    Limb borrow=0u;
    for (std::size_t i=0u;i<lhs.size();++i) {
        const Limb right=i<rhs.size() ? rhs[i] : Limb{0};
        const Limb left=lhs[i];
        const Limb difference=left-right;
        const Limb value=difference-borrow;
        borrow=static_cast<Limb>(left<right || difference<borrow);
        lhs[i]=value;
    }
    if (borrow!=0u) {
        throw std::logic_error("BigUnsigned Karatsuba negative cross term");
    }
}

// Accumulate the shifted cross term, touching only limbs affected by
// term or a propagated carry. In the no-carry case, the remaining output
// limbs are already correct and need not be read or rewritten.
// A widened sum-of-halves temporary can extend beyond the destination:
// every discarded word must still be checked for zero.
void add_shifted(Write output,Read term)
{
    const std::size_t common=std::min(output.size(),term.size());
    Limb carry=0u;
    std::size_t i=0u;
    for (;i<common;++i) {
        const Limb left=output[i];
        const Limb right=term[i];
        const Limb sum=left+right;
        const Limb value=sum+carry;
        carry=static_cast<Limb>(sum<left || value<sum);
        output[i]=value;
    }
    while (carry!=0u && i<output.size()) {
        const Limb left=output[i];
        output[i]=left+carry;
        carry=static_cast<Limb>(output[i]<left);
        ++i;
    }
    if (carry!=0u) {
        throw std::logic_error("BigUnsigned Karatsuba output overflow");
    }
    for (std::size_t j=output.size();j<term.size();++j) {
        if (term[j]!=0u) {
            throw std::logic_error("BigUnsigned Karatsuba nonzero high term");
        }
    }
}

// output.size()==a.size()+b.size(); scratch size is precomputed.
// z0 and z2 occupy disjoint output slices, so no extra result buffer
// is required for either of them. Only sums and cross term use scratch.
void karatsuba_recursive(
        Read a,Read b,Write output,Write scratch,
        std::size_t cutoff)
{
    if (!should_split(a.size(),b.size(),cutoff)) {
        BigUnsigned::mul_low_into(a,b,output);
        return;
    }

    // should_split() ensures both operands extend beyond half. Thus
    // a0.size()==b0.size()==half, z0 occupies [0,2*half) and z2
    // occupies [2*half,a.size()+b.size()). Their recursive base cases
    // overwrite every limb; clearing the entire parent is redundant.
    const std::size_t half=std::max(a.size(),b.size())/2u;
    const Read a0=a.first(std::min(a.size(),half));
    const Read b0=b.first(std::min(b.size(),half));
    const Read a1=a.subspan(a0.size());
    const Read b1=b.subspan(b0.size());

    karatsuba_recursive(
        a0,b0,output.first(a0.size()+b0.size()),scratch,cutoff);
    karatsuba_recursive(
        a1,b1,output.subspan(2u*half,a1.size()+b1.size()),
        scratch,cutoff);

    const std::size_t sa=std::max(a0.size(),a1.size())+1u;
    const std::size_t sb=std::max(b0.size(),b1.size())+1u;
    const std::size_t pair=sa+sb;

    Write sum_a=scratch.first(sa);
    Write sum_b=scratch.subspan(sa,sb);
    Write cross=scratch.subspan(pair,pair);
    Write nested=scratch.subspan(2u*pair);

    sum_halves(a0,a1,sum_a);
    sum_halves(b0,b1,sum_b);
    karatsuba_recursive(sum_a,sum_b,cross,nested,cutoff);

    // (a0+a1)(b0+b1)-a0*b0-a1*b1 = a0*b1+a1*b0.
    subtract_in_place(cross,output.first(a0.size()+b0.size()));
    subtract_in_place(
        cross,output.subspan(2u*half,a1.size()+b1.size()));
    add_shifted(output.subspan(half),cross);
}

} // namespace

std::size_t BigUnsigned::karatsuba_scratch_limbs(
        std::size_t a_limbs,std::size_t b_limbs)
{
    return karatsuba_scratch_limbs(
        a_limbs,b_limbs,karatsuba_default_cutoff);
}

std::size_t BigUnsigned::karatsuba_scratch_limbs(
        std::size_t a_limbs,std::size_t b_limbs,std::size_t cutoff)
{
    if (cutoff<karatsuba_minimum_cutoff) {
        throw std::invalid_argument("BigUnsigned Karatsuba cutoff must be >= 8");
    }
    (void)checked_sum(a_limbs,b_limbs);
    return karatsuba_workspace(a_limbs,b_limbs,cutoff);
}

void BigUnsigned::mul_into_karatsuba(
        LimbView a,LimbView b,
        MutableLimbView output,MutableLimbView scratch)
{
    mul_into_karatsuba(
        a,b,output,scratch,karatsuba_default_cutoff);
}

void BigUnsigned::mul_into_karatsuba(
        LimbView a,LimbView b,
        MutableLimbView output,MutableLimbView scratch,
        std::size_t cutoff)
{
    if (cutoff<karatsuba_minimum_cutoff) {
        throw std::invalid_argument("BigUnsigned Karatsuba cutoff must be >= 8");
    }
    if (a.empty() || b.empty()) {
        std::fill(output.begin(),output.end(),Limb{0});
        return;
    }

    const std::size_t count=checked_sum(a.size(),b.size());
    if (output.size()<count) {
        throw std::length_error("BigUnsigned Karatsuba output too small");
    }
    const std::size_t required=karatsuba_workspace(a.size(),b.size(),cutoff);
    if (scratch.size()<required) {
        throw std::length_error("BigUnsigned Karatsuba scratch too small");
    }

    // Recursion writes all count result limbs, including its base cases.
    // Only extra destination words need explicit clearing.
    std::fill(output.subspan(count).begin(),output.end(),Limb{0});
    karatsuba_recursive(
        a,b,output.first(count),scratch.first(required),cutoff);
}

} // namespace Ariadne
