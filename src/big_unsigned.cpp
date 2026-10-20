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
// B37: exact symmetric full square; requires 2*a.size() output limbs,
// disjoint from the canonical nonempty input, and at least eight limbs.
void square_into_symmetric_arm64(BigUnsigned::LimbView a,
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
    if (a._size<b._size) {
        throw std::underflow_error("BigUnsigned negative result");
    }

    const auto left=a.limbs();
    const auto right=b.limbs();
    std::size_t active=a._size;

    // B41: if the canonical operands have equal width, identical
    // leading limbs cancel exactly. The first different limb decides
    // the unsigned ordering; the result needs at most 'active' limbs.
    // This avoids allocating/subtracting dozens of high zero-result
    // limbs for nearly equal large integers. A dense random pair
    // generally differs in its first (highest) inspected limb.
    if (a._size==b._size) {
        while (active!=0u && left[active-1u]==right[active-1u]) {
            --active;
        }
        if (active==0u) {
            return BigUnsigned{};
        }
        if (left[active-1u]<right[active-1u]) {
            throw std::underflow_error("BigUnsigned negative result");
        }
    }

    BigUnsigned result;
    result.resize_for_result(active);
    const auto out=BigUnsigned::MutableLimbView(
        result.mutable_data(),active);
    // The high differing limb (or unequal canonical lengths) proves
    // a>=b, so subtraction cannot borrow out of the active range.
    (void)BigUnsigned::sub_into(
        left.first(active),right.first(std::min(active,right.size())),out);
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
#if defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
    // B48: a 32-bit divisor admits an exact 128/32 quotient via
    // two 64/32 hardware divisions, instead of the compiler's
    // 128/64 __udivti3 helper on ARM64. Since high<divisor<2^32,
    // both partial numerators fit 64 bits and the quotient halves
    // fit 32 bits. The same exact result serves joint, quotient-only
    // and remainder-only callers. Leave other targets and all
    // larger divisors on their existing paths.
    if (divisor<=std::uint64_t{0xffffffffu}) {
        constexpr std::uint64_t mask=0xffffffffULL;
        const std::uint64_t upper=(high<<32u)|(low>>32u);
        const std::uint64_t quotient_high=upper/divisor;
        const std::uint64_t upper_remainder=upper%divisor;
        const std::uint64_t lower=
            (upper_remainder<<32u)|(low&mask);
        const std::uint64_t quotient_low=lower/divisor;
        return {(quotient_high<<32u)|quotient_low,lower%divisor};
    }
#endif
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

// B57: exact 64-bit preinverted normalized 2-by-1 division.
// For B=2^64 and B/2 <= d < B, inverse=floor((B^2-1)/d)-B.
// The quotient approximation uses one wide product, two carry
// corrections and no DIVQ per limb; d is fixed for the whole stream.
// B63: GCC x86-64 selects it for >=16-limb quotient-producing
// callers, based on the GMP-verified B62 16/32/64-limb timings,
// and (B61) >=8-limb single-limb remainder-only callers.
// Clang, shorter inputs and other targets retain DIVQ; no BMI2/ADX.
#if defined(__x86_64__) && defined(__GNUC__) && \
    !defined(__clang__) && !defined(_MSC_VER) && \
    defined(__SIZEOF_INT128__)
std::pair<std::uint64_t,std::uint64_t> divide_preinverted_x86_64(
        std::uint64_t high,std::uint64_t low,
        std::uint64_t normalized,std::uint64_t inverse) noexcept
{
    __extension__ using Wide=unsigned __int128;
    const Wide estimate=static_cast<Wide>(high)*inverse;
    const Wide sum=estimate+
        (static_cast<Wide>(high+1u)<<64u)+low;
    std::uint64_t quotient=static_cast<std::uint64_t>(sum>>64u);
    const std::uint64_t estimate_low=static_cast<std::uint64_t>(sum);
    std::uint64_t remainder=low-quotient*normalized;
    const std::uint64_t low_estimate=
        static_cast<std::uint64_t>(remainder>estimate_low);
    quotient-=low_estimate;
    remainder+=low_estimate*normalized;
    if (remainder>=normalized) {
        ++quotient;
        remainder-=normalized;
    }
    return {quotient,remainder};
}

// B61: same exact B57 normalized pi1 used by quotient-producing division.
// Rem-only inputs need the remainder, not an entire quotient array.
std::uint64_t remainder_preinverted_x86_64(
        BigUnsigned::LimbView words,std::uint64_t divisor) noexcept
{
    const unsigned shift=static_cast<unsigned>(std::countl_zero(divisor));
    const std::uint64_t normalized=divisor<<shift;
    const std::uint64_t inverse=divide_wide(
        ~normalized,std::numeric_limits<std::uint64_t>::max(),
        normalized).first;
    std::uint64_t remaining=shift==0u ? std::uint64_t{0} :
        words.back()>>(64u-shift);
    for (std::size_t i=words.size();i!=0u;--i) {
        const std::uint64_t low=shift==0u ? words[i-1u] :
            (words[i-1u]<<shift)|
            (i>1u ? words[i-2u]>>(64u-shift) : std::uint64_t{0});
        remaining=divide_preinverted_x86_64(
            remaining,low,normalized,inverse).second;
    }
    return remaining>>shift;
}
#endif

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

    // B57: division by 2^k is an exact right shift and mask,
    // valid on every target; handle k=0 without an invalid shift.
    if (std::has_single_bit(divisor)) {
        const unsigned bits=static_cast<unsigned>(std::countr_zero(divisor));
        if (bits==0u) {
            std::copy(a.begin(),a.end(),quotient.begin());
        } else {
            for (std::size_t i=0u;i<a.size();++i) {
                const Limb upper=i+1u<a.size()?a[i+1u]:Limb{0};
                quotient[i]=(a[i]>>bits)|(upper<<(64u-bits));
            }
        }
        std::fill(quotient.subspan(a.size()).begin(),quotient.end(),
                  Limb{0});
        return a.empty()?Limb{0}:(a.front()&(divisor-1u));
    }

#if defined(__x86_64__) && defined(__GNUC__) && \
    !defined(__clang__) && !defined(_MSC_VER) && \
    defined(__SIZEOF_INT128__)
    // B63: the exact B62 GMP-verified pi1 kernel beats DIVQ
    // at 16, 32 and 64 limbs on the user's GCC x86-64 host.
    // Preserve DIVQ below 16 limbs and on other compiler/ISA paths.
    // Construct the exact reciprocal with one native 128/64 DIVQ:
    // (B^2-1)/normalized - B = ((B-1-normalized)*B+B-1)/normalized.
    if (a.size()>=16u) {
        const unsigned shift=static_cast<unsigned>(std::countl_zero(divisor));
        const Limb normalized=divisor<<shift;
        const Limb inverse=divide_wide(
            ~normalized,std::numeric_limits<Limb>::max(),
            normalized).first;
        Limb rem=shift==0u ? Limb{0} :
            a.back()>>(64u-shift);
        for (std::size_t i=a.size();i!=0u;--i) {
            const Limb low=shift==0u ? a[i-1u] :
                (a[i-1u]<<shift)|
                (i>1u ? a[i-2u]>>(64u-shift) : Limb{0});
            const auto [digit,next]=divide_preinverted_x86_64(
                rem,low,normalized,inverse);
            quotient[i-1u]=digit;
            rem=next;
        }
        std::fill(quotient.subspan(a.size()).begin(),quotient.end(),
                  Limb{0});
        return rem>>shift;
    }
#endif

#if defined(__aarch64__) && defined(__SIZEOF_INT128__) && \
    (defined(__GNUC__) || defined(__clang__))
    // B49: amortize a small divisor's reciprocal over a whole
    // multi-limb quotient. B48 needed two native UDIV instructions
    // for every limb (two exact 64/32 steps). For >=8 limbs, one
    // precomputed R=floor(2^64/divisor) replaces both divisions
    // with UMULH and at most one correction per 32-bit half.
    // The small-input path stays unchanged to avoid setup overhead.
    if (a.size()>=8u && divisor>1u &&
            divisor<=Limb{0xffffffffu}) {
        __extension__ using Wide=unsigned __int128;
        const Limb maximum=std::numeric_limits<Limb>::max();
        const Limb reciprocal=maximum/divisor+
            static_cast<Limb>(maximum%divisor==divisor-1u);
        const auto divide_half=[divisor,reciprocal](Limb numerator) {
            // Numerator < divisor*2^32 < 2^64.
            // R=floor(2^64/divisor) yields a quotient estimate
            // either exact or one low; residual is <2*divisor.
            Limb digit=static_cast<Limb>(
                (static_cast<Wide>(numerator)*reciprocal)>>64u);
            Limb residual=numerator-digit*divisor;
            if (residual>=divisor) {
                ++digit;
                residual-=divisor;
            }
            return std::pair<Limb,Limb>{digit,residual};
        };
        Limb remainder=0u;
        for (std::size_t i=a.size();i!=0u;--i) {
            const Limb word=a[i-1u];
            const Limb upper=(remainder<<32u)|(word>>32u);
            const auto [high_digit,partial]=divide_half(upper);
            const Limb lower=(partial<<32u)|
                (word&Limb{0xffffffffu});
            const auto [low_digit,next]=divide_half(lower);
            quotient[i-1u]=(high_digit<<32u)|low_digit;
            remainder=next;
        }
        std::fill(quotient.subspan(a.size()).begin(),quotient.end(),
                  Limb{0});
        return remainder;
    }
#endif
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

// B62: caller-owned quotient comparison for the SAME non-power-of-two
// single-limb inputs, without changing production's >=64 limb pi1 gate.
// Both helpers retain the public kernel's exact quotient-output contract.
// Defined only for the GCC x86-64 configuration under test.
#if defined(__x86_64__) && defined(__GNUC__) && \
    !defined(__clang__) && !defined(_MSC_VER) && \
    defined(__SIZEOF_INT128__)
BigUnsigned::Limb quotient_limb_pi1_into_diagnostic(
        BigUnsigned::LimbView a,BigUnsigned::Limb divisor,
        BigUnsigned::MutableLimbView quotient)
{
    using Limb=BigUnsigned::Limb;
    if (divisor==0u) {
        throw std::domain_error("B62 quotient pi1 division by zero");
    }
    if (quotient.size()<a.size()) {
        throw std::length_error("B62 quotient pi1 output too small");
    }
    if (a.empty() || std::has_single_bit(divisor)) {
        return BigUnsigned::divrem_limb_into(a,divisor,quotient);
    }
    const unsigned shift=static_cast<unsigned>(std::countl_zero(divisor));
    const Limb normalized=divisor<<shift;
    const Limb inverse=divide_wide(
        ~normalized,std::numeric_limits<Limb>::max(),normalized).first;
    Limb remaining=shift==0u ? Limb{0} : a.back()>>(64u-shift);
    for (std::size_t i=a.size();i!=0u;--i) {
        const Limb low=shift==0u ? a[i-1u] :
            (a[i-1u]<<shift)|
            (i>1u ? a[i-2u]>>(64u-shift) : Limb{0});
        const auto [digit,next]=divide_preinverted_x86_64(
            remaining,low,normalized,inverse);
        quotient[i-1u]=digit;
        remaining=next;
    }
    std::fill(quotient.begin()+static_cast<std::ptrdiff_t>(a.size()),
              quotient.end(),Limb{0});
    return remaining>>shift;
}

BigUnsigned::Limb quotient_limb_divq_into_diagnostic(
        BigUnsigned::LimbView a,BigUnsigned::Limb divisor,
        BigUnsigned::MutableLimbView quotient)
{
    using Limb=BigUnsigned::Limb;
    if (divisor==0u) {
        throw std::domain_error("B62 quotient DIVQ division by zero");
    }
    if (quotient.size()<a.size()) {
        throw std::length_error("B62 quotient DIVQ output too small");
    }
    Limb remaining=0u;
    for (std::size_t i=a.size();i!=0u;--i) {
        const auto [digit,next]=divide_wide(remaining,a[i-1u],divisor);
        quotient[i-1u]=digit;
        remaining=next;
    }
    std::fill(quotient.begin()+static_cast<std::ptrdiff_t>(a.size()),
              quotient.end(),Limb{0});
    return remaining;
}
#endif

// B61: independent diagnostic reference to the *old* owning single-limb
// DIVQ remainder. This preserves the B60 baseline after pi1 promotion.
BigUnsigned rem_limb_divq_diagnostic(
        BigUnsigned const& a,BigUnsigned const& b)
{
    if (b.size()!=1u) {
        return rem(a,b);
    }
    const auto words=a.limbs();
    const BigUnsigned::Limb divisor=b.limb(0u);
    BigUnsigned::Limb remaining=0u;
    for (std::size_t i=words.size();i!=0u;--i) {
        remaining=divide_wide(remaining,words[i-1u],divisor).second;
    }
    return BigUnsigned(remaining);
}

// Preserve the B60 diagnostic pi1 comparison on the exact same operands.
// Production now selects the identical preinverse kernel only on GCC x86-64
// for non-power-of-two single-limb divisors and >=8-limb dividends.
BigUnsigned rem_limb_preinverse_diagnostic(
        BigUnsigned const& a,BigUnsigned const& b)
{
#if defined(__x86_64__) && defined(__GNUC__) && \
    !defined(__clang__) && !defined(_MSC_VER) && \
    defined(__SIZEOF_INT128__)
    if (b.size()==1u && a.size()>=8u &&
            !std::has_single_bit(b.limb(0u))) {
        return BigUnsigned(remainder_preinverted_x86_64(
            a.limbs(),b.limb(0u)));
    }
#endif
    return rem(a,b);
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

// B55: x86-64 GNU/Clang baseline MULQ/ADDQ/ADCQ/SUBQ/ADCQ
// specialization for the existing fused digit-product subtraction.
// It computes high:low = divisor[i]*digit + carry, subtracts
// low from window[i], then folds the subtraction borrow into high.
// The exact recurrence and (m+1)-limb two's-complement result
// match subtract_digit_product<SparseZero>. Skipping a zero divisor
// limb still subtracts the incoming carry/borrow. No BMI2/ADX.
// Production selects this kernel for divisor lengths >=8 limbs;
// smaller divisors retain the C++ wide implementation.
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && \
    !defined(_MSC_VER)
template<bool SparseZero>
bool subtract_digit_product_x86_64(
        BigUnsigned::MutableLimbView window,
        BigUnsigned::LimbView divisor,
        BigUnsigned::Limb digit) noexcept
{
    using Limb=BigUnsigned::Limb;
    Limb carry=0u;
    for (std::size_t i=0u;i<divisor.size();++i) {
        if constexpr (SparseZero) {
            if (divisor[i]==0u) {
                const Limb old=window[i];
                window[i]=old-carry;
                carry=static_cast<Limb>(old<carry);
                continue;
            }
        }
        Limb low=divisor[i],high;
        Limb old=window[i];
        __asm__ volatile(
            "mulq %[digit]\n\t"
            "addq %[carry], %%rax\n\t"
            "adcq $0, %%rdx\n\t"
            "subq %%rax, %[old]\n\t"
            "adcq $0, %%rdx\n\t"
            : "+&a" (low), "=&d" (high), [old] "+&r" (old)
            : [digit] "r" (digit), [carry] "r" (carry)
            : "cc");
        window[i]=old;
        carry=high;
    }
    const std::size_t m=divisor.size();
    const Limb old=window[m];
    window[m]=old-carry;
    return old<carry;
}
#endif

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

    // B57: detect an exact power of B=2^64 times a power of two.
    // For b=B^k*2^s, q=a>>(64*k+s), r=a mod b. This bypasses
    // all quotient estimates and product-subtractions, for every
    // output mode, including remainder-only. No temporary storage.
    const bool power_of_two=std::has_single_bit(b.back()) &&
        (b.size()==1u ||
         (b.front()==0u &&
          std::all_of(b.begin()+1u,b.end()-1u,
                      [](Limb word) { return word==0u; })));
    if (power_of_two) {
        const std::size_t skip=b.size()-1u;
        const unsigned bits=static_cast<unsigned>(std::countr_zero(b.back()));
        if constexpr (WantQuotient) {
            for (std::size_t i=0u;i<quotient_limbs;++i) {
                const std::size_t index=i+skip;
                const Limb upper=index+1u<a.size() ?
                    a[index+1u] : Limb{0};
                quotient[i]=bits==0u ? a[index] :
                    (a[index]>>bits)|(upper<<(64u-bits));
            }
        }
        if constexpr (WantRemainder) {
            std::copy_n(a.begin(),skip,remainder.begin());
            remainder[skip]=a[skip]&(b.back()-1u);
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
            // B61: exact normalized pi1 speeds up 8+ limb rem-only
            // division on GCC x86-64, without quotient allocation.
            // Preserve short DIVQ, power-of-two, ARM64 and Clang paths.
#if defined(__x86_64__) && defined(__GNUC__) && \
    !defined(__clang__) && !defined(_MSC_VER) && \
    defined(__SIZEOF_INT128__)
            if (a.size()>=8u) {
                remainder[0u]=remainder_preinverted_x86_64(a,b[0u]);
            } else
#endif
            {
                Limb low=0u;
                for (std::size_t i=a.size();i!=0u;--i) {
                    low=divide_wide(low,a[i-1u],b[0u]).second;
                }
                remainder[0u]=low;
            }
        }
        return;
    }

    // B45: remove an exact B^k factor (B=2^64) from the divisor.
    // a = a_hi*B^k+a_lo and b = b_hi*B^k, 0 <= a_lo < B^k:
    // q = a_hi/b_hi and r = (a_hi%b_hi)*B^k+a_lo.
    // Reuse the already checked and zero-initialized output/scratch;
    // the high-span division needs strictly less scratch storage.
    if (b.front()==0u) {
        std::size_t zero_limbs=1u;
        while (zero_limbs+1u<b.size() && b[zero_limbs]==0u) {
            ++zero_limbs;
        }
        divrem_selected_into<WantQuotient,WantRemainder>(
            a.subspan(zero_limbs),b.subspan(zero_limbs),quotient,
            WantRemainder ? remainder.subspan(zero_limbs)
                          : MutableLimbView{},
            scratch);
        if constexpr (WantRemainder) {
            std::copy_n(a.begin(),zero_limbs,remainder.begin());
        }
        return;
    }

    const std::size_t n=a.size();
    const std::size_t m=b.size();

    // B36: for equal-length canonical operands the high limb can
    // establish an exact quotient of zero or one without normalization.
    // If ah<bh, then a<b. If bh<ah<2*bh, then b<=a<2*b.
    // Test ah-bh<bh without doubling bh (which could overflow).
    // Equal high limbs retain B28's exact borrow/compare path.
    // All output and scratch requirements were checked above.
    if (n==m) {
        const Limb ah=a.back(),bh=b.back();
        if (ah<bh) {
            // Quotient is already zero-initialized.
            if constexpr (WantRemainder) {
                std::copy(a.begin(),a.end(),remainder.begin());
            }
            return;
        }
        if (ah==bh || ah-bh<bh) {
            if constexpr (WantRemainder) {
                Limb borrow=0u;
                for (std::size_t i=0u;i<m;++i) {
                    const Limb difference=a[i]-b[i];
                    remainder[i]=difference-borrow;
                    borrow=static_cast<Limb>((a[i]<b[i]) ||
                                             (difference<borrow));
                }
                if (borrow!=0u) {
                    // Possible only if the high limbs are equal.
                    std::copy(a.begin(),a.end(),remainder.begin());
                }
                if constexpr (WantQuotient) {
                    quotient[0u]=static_cast<Limb>(borrow==0u);
                }
            } else if constexpr (WantQuotient) {
                if (ah>bh) {
                    // Strict high-limb dominance proves quotient 1.
                    quotient[0u]=1u;
                } else {
                    // B28: equal high limbs require a lower-limb scan.
                    for (std::size_t i=m;i!=0u;--i) {
                        if (a[i-1u]!=b[i-1u]) {
                            quotient[0u]=static_cast<Limb>(a[i-1u]>b[i-1u]);
                            return;
                        }
                    }
                    quotient[0u]=1u;
                }
            }
            return;
        }
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
            // B55: a single selector outside the hot inner limb loop.
            // >=8-limb x86 divisors use an exact native MULQ/ADCQ
            // multiply/subtract recurrence; all other targets and
            // smaller denominators keep the original wide C++ code.
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__)) && \
    !defined(_MSC_VER)
            bool overestimated=m>=8u ?
                subtract_digit_product_x86_64<SparseZero>(
                    window,v,digit) :
                subtract_digit_product<SparseZero>(window,v,digit);
#else
            bool overestimated=
                subtract_digit_product<SparseZero>(window,v,digit);
#endif
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

    // B47/B50: use native wide multiplication for medium/long
    // decimal inputs; B47 lowered the threshold to 32 digits.
    // B50 packs 19 digits per 64-bit chunk (10^19 < 2^64), one
    // extra digit per full limb scan versus radix 10^18. Keep the
    // original B19 radix-10^9 portable path for smaller inputs
    // and toolchains without unsigned 128-bit arithmetic.
#if defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    if (input.size()>=32u) {
        __extension__ using Wide=unsigned __int128;
        constexpr std::size_t group=19u;
        constexpr BigUnsigned::Limb full_radix=10000000000000000000ULL;
        std::vector<BigUnsigned::Limb> digits;
        // ceil(decimal_digits/19) bounds the necessary 64-bit
        // limbs because 10^19 < 2^64, even for leading zeros.
        digits.reserve(input.size()/group+
                       static_cast<std::size_t>(input.size()%group!=0u));
        // Consume a possibly short first group, then multiply by
        // the compile-time constant 10^19 on all full groups.
        // Every multiply+carry remains strictly below 2^128.
        const std::size_t first=input.size()%group==0u ?
            group : input.size()%group;
        std::size_t pos=0u;
        while (pos<input.size()) {
            const std::size_t count=pos==0u ? first : group;
            BigUnsigned::Limb radix=full_radix;
            if (count!=group) {
                radix=1u;
                for (std::size_t i=0u;i<count;++i) radix*=10u;
            }
            BigUnsigned::Limb chunk=0u;
            for (std::size_t i=0u;i<count;++i) {
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

// B51: B46's exact reciprocal removed the runtime/inline 128/64
// division that made earlier B25/B30 radix-10^18 x86-64 attempts
// slower than the base-10^9 implementation. Permit the new
// multiply-only 18-digit formatter on both AArch64 and x86-64
// GCC/Clang with 128-bit integer support. The B19 base-10^9 path
// remains selected for <=8 limbs and unsupported toolchains,
// including MSVC. No CPU extension (BMI2/ADX) is required.
#if (defined(__aarch64__) || defined(__x86_64__)) && \
    defined(__SIZEOF_INT128__) && \
    (defined(__GNUC__) || defined(__clang__))
std::string to_decimal_blocks_wide(BigUnsigned const& value)
{
    __extension__ using Wide=unsigned __int128;
    using Limb=BigUnsigned::Limb;
    constexpr Limb radix=1000000000000000000ULL;
    constexpr std::size_t block_width=18u;

    // B46: replace the compiler's 128/64 runtime division (__udivti3
    // on ARM64 Clang) with an exact constant-divisor reciprocal.
    // R=floor(2^128/radix) (2^128 is not divisible by radix).
    // For N=high*2^64+low and high<radix, floor(N*R/2^128) differs
    // from floor(N/radix) by at most one, always on the low side:
    // 0 <= N/radix-N*R/2^128 < radix/2^64 < 1.
    // The provisional remainder is <2*radix<2^64, so its low
    // 64-bit subtraction is the exact nonnegative remainder.
    constexpr Wide reciprocal=(~Wide{0})/static_cast<Wide>(radix);
    constexpr Limb reciprocal_high=static_cast<Limb>(reciprocal>>64u);
    constexpr Limb reciprocal_low=static_cast<Limb>(reciprocal);
    static_assert(reciprocal_high==18u);
    static_assert(radix<=std::numeric_limits<Limb>::max()/2u);

    const auto limbs=value.limbs();
    std::vector<BigUnsigned::Limb> work(limbs.begin(),limbs.end());
    std::vector<BigUnsigned::Limb> blocks;
    // Two base-10^18 blocks per 64-bit limb are always sufficient.
    if (limbs.size()<=blocks.max_size()/2u) {
        blocks.reserve(limbs.size()*2u);
    }
    while (!work.empty()) {
        BigUnsigned::Limb remainder=0u;
        for (std::size_t i=work.size();i!=0u;--i) {
            // High 128 bits of (remainder*2^64+digit)*R.
            // Compute its quotient with only 64x64->128 products.
            // The cross-term fits in Wide because remainder<radix.
            const Limb digit=work[i-1u];
            const Wide cross=
                static_cast<Wide>(remainder)*reciprocal_low+
                static_cast<Wide>(digit)*reciprocal_high+
                ((static_cast<Wide>(digit)*reciprocal_low)>>64u);
            Limb quotient=remainder*reciprocal_high+
                          static_cast<Limb>(cross>>64u);
            Limb residual=digit-quotient*radix;
            if (residual>=radix) {
                ++quotient;
                residual-=radix;
            }
            work[i-1u]=quotient;
            remainder=residual;
        }
        blocks.push_back(remainder);
        while (!work.empty() && work.back()==0u) {
            work.pop_back();
        }
    }

    // B42: one exactly sized output and no per-block temporary strings.
    static constexpr auto digit_pairs=[] {
        std::array<char,200> table{};
        for (unsigned p=0u;p<100u;++p) {
            table[2u*p]=static_cast<char>('0'+p/10u);
            table[2u*p+1u]=static_cast<char>('0'+p%10u);
        }
        return table;
    }();
    std::array<char,20> head{};
    const auto [end,error]=std::to_chars(
        head.data(),head.data()+head.size(),blocks.back());
    if (error!=std::errc{}) {
        throw std::length_error("BigUnsigned decimal head overflow");
    }
    const std::size_t top_width=static_cast<std::size_t>(end-head.data());
    std::string result;
    if (blocks.size()-1u>(result.max_size()-top_width)/block_width) {
        throw std::length_error("BigUnsigned decimal output overflow");
    }
    result.resize(top_width+(blocks.size()-1u)*block_width);
    std::copy_n(head.data(),top_width,result.data());
    char* output=result.data()+result.size();
    for (std::size_t i=0u;i+1u<blocks.size();++i) {
        BigUnsigned::Limb block=blocks[i];
        for (unsigned p=0u;p<block_width/2u;++p) {
            const unsigned pair=static_cast<unsigned>(block%100u);
            block/=100u;
            output-=2u;
            output[0u]=digit_pairs[2u*pair];
            output[1u]=digit_pairs[2u*pair+1u];
        }
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
#if (defined(__aarch64__) || defined(__x86_64__)) && \
    defined(__SIZEOF_INT128__) && \
    (defined(__GNUC__) || defined(__clang__))
    // B51: reuse B46's exact reciprocal formatter on both targets
    // for >=9 limbs. Keep 8-limb and smaller inputs on radix-10^9
    // so the small-value path and its allocation profile stay fixed.
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

// B59: the user's GMP-verified x86-64 GCC B58 experiment favored 32
// limbs at 2048 and 4096 bits. Preserve the measured 16-limb
// production baseline on all other compiler/architecture combinations.
#if defined(__x86_64__) && defined(__GNUC__) && !defined(__clang__) && \
    !defined(_MSC_VER) && defined(__SIZEOF_INT128__)
constexpr std::size_t decimal_recursive_cutoff_limbs=32u;
#else
constexpr std::size_t decimal_recursive_cutoff_limbs=16u;
#endif

// B43: append into a single owning output string throughout the
// divide-and-conquer recursion. min_digits is the exact field width
// required by an ancestor's low half; only a leaf emits zero padding.
// No leading zeros are introduced for the root (min_digits==0).
template<std::size_t CutoffLimbs>
void append_decimal_recursive_with_cutoff(
        BigUnsigned const& value,
        std::vector<DecimalSplitPower> const& powers,
        std::string& output,std::size_t min_digits)
{
    if (value.size()>CutoffLimbs) {
        // Preserve B22's balanced power selection and exact divrem.
        const std::size_t midpoint=value.bit_length()/2u;
        std::size_t count=powers.size();
        while (count!=0u &&
                powers[count-1u].value.bit_length()>midpoint) {
            --count;
        }
        if (count!=0u) {
            auto const& divisor=powers[count-1u];
            auto [high,low]=divrem(value,divisor.value);

            // The quotient contributes the high digits. If an
            // ancestor requested a fixed-width decimal field, the
            // high part carries only the width not owned by 'low'.
            const std::size_t high_width=min_digits>divisor.width
                ? min_digits-divisor.width : 0u;
            append_decimal_recursive_with_cutoff<CutoffLimbs>(high,powers,output,high_width);
            // 0 <= low < 10^width, so its representation must occupy
            // exactly 'width' characters, including leading zeros.
            append_decimal_recursive_with_cutoff<CutoffLimbs>(low,powers,output,divisor.width);
            return;
        }
    }

    const std::string leaf=to_decimal_blocks(value);
    if (leaf.size()<min_digits) {
        output.append(min_digits-leaf.size(),'0');
    }
    output+=leaf;
}


// One shared exact implementation for the selected production cutoff
// and the B58/B59 diagnostic controls, avoiding algorithm drift.
template<std::size_t CutoffLimbs>
std::string decimal_with_cutoff(BigUnsigned const& value)
{
    if (value.is_zero()) {
        return "0";
    }
    if (value.size()<=CutoffLimbs) {
        return to_decimal_blocks(value);
    }

    // Squaring successive powers of 10^9 costs logarithmically many
    // products. Each recursive step uses the shared exact divrem;
    // no GMP boundary conversions, global caches, or ISA assumptions.
    std::vector<DecimalSplitPower> powers;
    powers.push_back({BigUnsigned(std::uint64_t{1000000000u}),9u});
    const std::size_t bit_length=value.bit_length();
    const std::size_t midpoint=bit_length/2u;
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
    std::string result;
    // floor(bits/3)+1 safely bounds decimal digits since
    // log10(2)<1/3. Reserve once and avoid repeated intermediate
    // string allocation/copying, for both x86-64 and ARM64.
    const std::size_t digit_capacity=bit_length/3u+1u;
    if (digit_capacity<=result.max_size()) {
        result.reserve(digit_capacity);
    }
    append_decimal_recursive_with_cutoff<CutoffLimbs>(value,powers,result,0u);
    return result;
}

} // namespace

std::string to_decimal(BigUnsigned const& value)
{
    return decimal_with_cutoff<decimal_recursive_cutoff_limbs>(value);
}

// B59 diagnostic control: preserve the original production-16 baseline.
std::string to_decimal_cutoff16_diagnostic(BigUnsigned const& value)
{
    return decimal_with_cutoff<16u>(value);
}

// B58 diagnostic control: preserve the experimentally measured 32 variant.
std::string to_decimal_cutoff32_diagnostic(BigUnsigned const& value)
{
    return decimal_with_cutoff<32u>(value);
}

namespace {

// Compiler-supported wide accumulation, shared by GCC/Clang x86-64 and ARM64.
// The intermediate is bounded by (B-1)^2+2*(B-1)=B^2-1, B=2^64.
// Preserve the historical MSVC x64 native B3 fallback. On other toolchains
// without unsigned __int128, retain portable multiplication as a fallback.
// Inputs and output must be disjoint; all output limbs are overwritten.
// B53: x86-64 GNU counterpart of ARM64 B38/B40: two compile-time
// recurrences, one skipping exact-zero outer rows and one without a
// branch in dense rows. A skipped row contributes zero and its final
// carry slot has never been touched by an earlier row. No input aliasing.
#if defined(ARIADNE_BIG_UNSIGNED_X86_64) && \
    defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED) && \
    defined(__x86_64__) && defined(__SIZEOF_INT128__) && \
    defined(__GNUC__) && !defined(__clang__) && !defined(_MSC_VER)
bool use_sparse_fused_rows_x86_64(BigUnsigned::LimbView rows) noexcept
{
    const std::size_t threshold=std::max<std::size_t>(
        2u,rows.size()/8u);
    std::size_t zeros=0u;
    for (const BigUnsigned::Limb digit:rows) {
        if (digit==0u && ++zeros>=threshold) {
            return true;
        }
    }
    return false;
}

template<bool SkipZeroRows>
void multiply_limbs_fused_x86_64_impl(
        BigUnsigned::LimbView a,BigUnsigned::LimbView b,
        BigUnsigned::MutableLimbView output) noexcept
{
    __extension__ using Wide=unsigned __int128;
    using Limb=BigUnsigned::Limb;
    std::fill(output.begin(),output.end(),Limb{0});
    for (std::size_t i=0u;i<a.size() && i<output.size();++i) {
        const Limb ai=a[i];
        if constexpr (SkipZeroRows) {
            if (ai==0u) {
                continue;
            }
        }
        Limb carry=0u;
        const std::size_t terms=std::min(b.size(),output.size()-i);
        for (std::size_t j=0u;j<terms;++j) {
            const Wide total=static_cast<Wide>(ai)*b[j]+
                             output[i+j]+carry;
            output[i+j]=static_cast<Limb>(total);
            carry=static_cast<Limb>(total>>64u);
        }
        if (i+terms<output.size()) {
            output[i+terms]=carry;
        }
    }
}
#endif

void multiply_limbs_fused(BigUnsigned::LimbView a,
                          BigUnsigned::LimbView b,
                          BigUnsigned::MutableLimbView output) noexcept
{
#if defined(_MSC_VER) && defined(_M_X64) && defined(ARIADNE_BIG_UNSIGNED_X86_64)
    BigUnsignedDetail::mul_low_into_x86_64(a,b,output);
#elif defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
#if defined(ARIADNE_BIG_UNSIGNED_X86_64) && \
    defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED) && \
    defined(__x86_64__) && defined(__GNUC__) && \
    !defined(__clang__) && !defined(_MSC_VER)
    // B40's ARM64 threshold: for small rows avoid an extra density
    // prescan, for large rows pick a branch-free dense recurrence
    // unless enough exact-zero rows justify skipping them.
    const std::size_t active=std::min(a.size(),output.size());
    if (active<64u || b.size()<64u ||
            use_sparse_fused_rows_x86_64(a.first(active))) {
        multiply_limbs_fused_x86_64_impl<true>(a,b,output);
    } else {
        multiply_limbs_fused_x86_64_impl<false>(a,b,output);
    }
#else
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
#endif
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
    // B54: on the SAME x86-64 host, B53's updated fused recurrence
    // beat the B52 32x32 MULQ/ADCQ special case on dense inputs
    // (801.142 vs 858.091 ns). Re-evaluate fused for all GCC sizes,
    // including 32-limb Karatsuba leaves. The forced assembly entry
    // point remains available for independent benchmarking.
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

#if defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED) && \
    defined(__x86_64__) && defined(__SIZEOF_INT128__) && \
    defined(__GNUC__) && !defined(__clang__) && !defined(_MSC_VER)
    // B72: B71's independent GMP-verified same-library-TU stack
    // construction was faster for short owning products on GCC x86-64,
    // especially four-limb inputs. Use the SAME selected fused kernel,
    // with bounded stack workspace. The returned BigUnsigned is fully
    // owned and normalized; the source buffers never escape.
    // Above four limbs, retain the previous heap and Karatsuba paths.
    if (a._size<=BigUnsigned::inline_capacity &&
            b._size<=BigUnsigned::inline_capacity) {
        std::array<BigUnsigned::Limb,2u*BigUnsigned::inline_capacity> digits{};
        const std::size_t count=a._size+b._size;
        BigUnsigned::mul_into(
            a.limbs(),b.limbs(),
            BigUnsigned::MutableLimbView(digits.data(),count));
        return BigUnsigned(BigUnsigned::LimbView(digits.data(),count));
    }
#endif

    BigUnsigned result;
    result.resize_for_result(a._size+b._size);

#if defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED) && \
    defined(__x86_64__) && defined(__SIZEOF_INT128__) && \
    defined(__GNUC__) && !defined(__clang__) && !defined(_MSC_VER)
    // B54: B35's balanced, large-size Karatsuba dispatch is
    // worthwhile for dense operands; B53 made sparse schoolbook
    // dramatically cheaper. Port ARM64 B39's exact density gate
    // so sparse powers skip unnecessary recursion and scratch.
    // Avoid multiplication in the size/balance comparison.
    const std::size_t smaller=std::min(a._size,b._size);
    const std::size_t larger=std::max(a._size,b._size);
    if (smaller>=64u && smaller>=larger-larger/4u) {
        const auto sufficiently_dense=[](BigUnsigned::LimbView words) {
            const std::size_t permitted_zeros=words.size()/8u;
            std::size_t zeros=0u;
            for (const BigUnsigned::Limb word:words) {
                if (word==0u && ++zeros>permitted_zeros) {
                    return false;
                }
            }
            return true;
        };
        if (sufficiently_dense(a.limbs()) &&
                sufficiently_dense(b.limbs())) {
            constexpr std::size_t karatsuba_cutoff=32u;
            const std::size_t scratch_limbs=
                BigUnsigned::karatsuba_scratch_limbs(
                    a._size,b._size,karatsuba_cutoff);
            std::vector<BigUnsigned::Limb> scratch(scratch_limbs);
            BigUnsigned::mul_into_karatsuba(
                a.limbs(),b.limbs(),
                BigUnsigned::MutableLimbView(
                    result.mutable_data(),a._size+b._size),
                scratch,karatsuba_cutoff);
            result.normalize();
            return result;
        }
    }
#endif
#if defined(ARIADNE_BIG_UNSIGNED_ARM64) && defined(__aarch64__) && \
    defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    // B39: on ARM64, Karatsuba with cutoff 32 beats the dense B37
    // schoolbook for balanced owning 64-limb products in the user's
    // independent preallocated benchmark. Allocate private scratch
    // only when both magnitudes are large, balanced and substantially
    // dense. Sparse powers/products retain B38's zero-row kernel.
    const std::size_t smaller=std::min(a._size,b._size);
    const std::size_t larger=std::max(a._size,b._size);
    if (smaller>=64u && smaller>=larger-larger/4u) {
        const auto sufficiently_dense=[](BigUnsigned::LimbView words) {
            const std::size_t permitted_zeros=words.size()/8u;
            std::size_t zeros=0u;
            for (const BigUnsigned::Limb word:words) {
                if (word==0u && ++zeros>permitted_zeros) {
                    return false;
                }
            }
            return true;
        };
        if (sufficiently_dense(a.limbs()) &&
                sufficiently_dense(b.limbs())) {
            constexpr std::size_t karatsuba_cutoff=32u;
            const std::size_t scratch_limbs=
                BigUnsigned::karatsuba_scratch_limbs(
                    a._size,b._size,karatsuba_cutoff);
            std::vector<BigUnsigned::Limb> scratch(scratch_limbs);
            BigUnsigned::mul_into_karatsuba(
                a.limbs(),b.limbs(),
                BigUnsigned::MutableLimbView(
                    result.mutable_data(),a._size+b._size),
                scratch,karatsuba_cutoff);
            result.normalize();
            return result;
        }
    }
#endif

    BigUnsigned::mul_into(
        a.limbs(),b.limbs(),
        BigUnsigned::MutableLimbView(result.mutable_data(),a._size+b._size));
    result.normalize();
    return result;
}

// B71: short owning multiplication construction A/B control.
// Keep BOTH alternatives in the same library TU as production mul().
// They use a local 8-limb stack result for <=4-limb operands and then
// construct a canonical BigUnsigned. The only kernel difference between
// the two diagnostic variants is forced fused vs selected dispatch.
// No production selector, headers, ABI contract or arithmetic changes.
BigUnsigned mul_short_stack_selected_diagnostic(
        BigUnsigned const& a,BigUnsigned const& b)
{
#if defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED) && \
    defined(__x86_64__) && defined(__SIZEOF_INT128__) && \
    defined(__GNUC__) && !defined(__clang__) && !defined(_MSC_VER)
    if (a.size()<=4u && b.size()<=4u) {
        if (a.is_zero() || b.is_zero()) return BigUnsigned{};
        std::array<BigUnsigned::Limb,8u> output{};
        const std::size_t count=a.size()+b.size();
        BigUnsigned::mul_into(
            a.limbs(),b.limbs(),
            BigUnsigned::MutableLimbView(output.data(),count));
        return BigUnsigned(BigUnsigned::LimbView(output.data(),count));
    }
#endif
    return mul(a,b);
}

BigUnsigned mul_short_stack_fused_diagnostic(
        BigUnsigned const& a,BigUnsigned const& b)
{
#if defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED) && \
    defined(__x86_64__) && defined(__SIZEOF_INT128__) && \
    defined(__GNUC__) && !defined(__clang__) && !defined(_MSC_VER)
    if (a.size()<=4u && b.size()<=4u) {
        if (a.is_zero() || b.is_zero()) return BigUnsigned{};
        std::array<BigUnsigned::Limb,8u> output{};
        const std::size_t count=a.size()+b.size();
        BigUnsigned::mul_into_fused(
            a.limbs(),b.limbs(),
            BigUnsigned::MutableLimbView(output.data(),count));
        return BigUnsigned(BigUnsigned::LimbView(output.data(),count));
    }
#endif
    return mul(a,b);
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
template<bool SkipZeroRows>
void square_limbs_symmetric_impl(
        BigUnsigned::LimbView a,
        BigUnsigned::MutableLimbView output) noexcept
{
    using Limb=BigUnsigned::Limb;
    __extension__ using Wide=unsigned __int128;
    const std::size_t n=a.size();
    std::fill(output.begin(),output.end(),Limb{0});

    // Accumulate sum_{i<j} a[i]*a[j]*B^(i+j) with one wide product
    // per unordered pair. Slot i+n has not been touched by prior rows.
    for (std::size_t i=0u;i+1u<n;++i) {
        const Limb ai=a[i];
        // B54/B38: skipped rows contribute zero. Earlier rows cannot
        // touch output[i+n], which was cleared at function entry.
        if constexpr (SkipZeroRows) {
            if (ai==0u) {
                continue;
            }
        }
        Limb carry=0u;
        for (std::size_t j=i+1u;j<n;++j) {
            const std::size_t k=i+j;
            const Wide sum=static_cast<Wide>(ai)*a[j]+
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

// B54: mirror ARM64 B40's adaptive square selection. The triangular
// part has n-1 active outer rows. For <=64 limbs skip the prescan;
// otherwise choose a branch-free dense variant unless the existing
// B53 zero-row classifier detects enough exact-zero rows. Doubling
// and all diagonal products are performed in BOTH variants.
void square_limbs_symmetric(BigUnsigned::LimbView a,
                            BigUnsigned::MutableLimbView output) noexcept
{
    if (a.size()<=64u ||
            use_sparse_fused_rows_x86_64(a.first(a.size()-1u))) {
        square_limbs_symmetric_impl<true>(a,output);
    } else {
        square_limbs_symmetric_impl<false>(a,output);
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
#if defined(ARIADNE_BIG_UNSIGNED_ARM64) && defined(__aarch64__) && \
    defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    // B37: native Apple Silicon uses a symmetric square instead of
    // computing a generic product twice for each off-diagonal pair.
    // Leave small inputs on the already-tuned carry-chain schoolbook.
    if (a._size>=8u) {
        if (a._size>std::numeric_limits<std::size_t>::max()/2u) {
            throw std::length_error("BigUnsigned square size overflow");
        }
        BigUnsigned result;
        result.resize_for_result(2u*a._size);
        BigUnsignedDetail::square_into_symmetric_arm64(
            a.limbs(),BigUnsigned::MutableLimbView(
                result.mutable_data(),2u*a._size));
        result.normalize();
        return result;
    }
#endif
    return mul(a,a);
}


// B65: diagnostic-only GCC short symmetric square control.
// Reuse the B32 exact triangular/doubling/diagonal recurrence at 1..7
// limbs, where the production dispatcher still prefers schoolbook.
// A bounded stack output avoids inventing additional heap allocations.
// No production selector or public header is modified.
BigUnsigned square_short_symmetric_diagnostic(BigUnsigned const& input)
{
#if defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED) && \
    defined(__x86_64__) && defined(__SIZEOF_INT128__) && \
    defined(__GNUC__) && !defined(__clang__) && !defined(_MSC_VER)
    if (!input.is_zero() && input.size()<=7u) {
        std::array<BigUnsigned::Limb,14u> digits{};
        const std::size_t count=2u*input.size();
        square_limbs_symmetric(
            input.limbs(),BigUnsigned::MutableLimbView(digits.data(),count));
        return BigUnsigned(BigUnsigned::LimbView(digits.data(),count));
    }
#endif
    return square(input);
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
