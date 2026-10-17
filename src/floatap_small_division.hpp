/***************************************************************************
 *  floatap_small_division.hpp
 *  Unsigned fixed-capacity 32-bit-limb arithmetic for FloatAP.
 *  No GMP/MPFR types or compiler-specific wide integers.
 ***************************************************************************/
#ifndef ARIADNE_NUMERIC_FLOATAP_SMALL_DIVISION_HPP
#define ARIADNE_NUMERIC_FLOATAP_SMALL_DIVISION_HPP

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>

namespace Ariadne::FloatAPSmallDivision {

constexpr std::size_t capacity = 18;
constexpr std::uint64_t base = std::uint64_t{1} << 32;

struct Unsigned {
    std::array<std::uint32_t,capacity> limb{};
    std::size_t size = 0;
    void trim() noexcept {
        while (size != 0 && limb[size-1] == 0) --size;
    }
};

inline std::size_t bit_length(Unsigned const& x) noexcept {
    return x.size == 0 ? 0 :
        32 * (x.size-1) + 32 - std::countl_zero(x.limb[x.size-1]);
}

inline int compare(Unsigned const& a, Unsigned const& b) noexcept {
    if (a.size != b.size) return a.size < b.size ? -1 : 1;
    for (std::size_t i = a.size; i != 0; --i) {
        if (a.limb[i-1] != b.limb[i-1])
            return a.limb[i-1] < b.limb[i-1] ? -1 : 1;
    }
    return 0;
}

inline bool shift_left(Unsigned const& x, std::size_t shift, Unsigned& out) noexcept {
    out = {};
    if (x.size == 0) return true;
    const std::size_t words = shift / 32;
    const unsigned bits = static_cast<unsigned>(shift % 32);
    if (words >= capacity || x.size > capacity-words) return false;
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < x.size; ++i) {
        const std::uint64_t part = (std::uint64_t{x.limb[i]} << bits) | carry;
        out.limb[i+words] = static_cast<std::uint32_t>(part);
        carry = part >> 32;
    }
    out.size = words+x.size;
    if (carry != 0) {
        if (out.size == capacity) return false;
        out.limb[out.size++] = static_cast<std::uint32_t>(carry);
    }
    return true;
}

inline void shift_right(Unsigned const& x, unsigned bits, Unsigned& out) noexcept {
    out = {};
    if (x.size == 0) return;
    if (bits == 0) { out = x; return; }
    std::uint32_t carry = 0;
    for (std::size_t i = x.size; i != 0; --i) {
        const std::uint32_t here = x.limb[i-1];
        out.limb[i-1] = static_cast<std::uint32_t>(
            (here >> bits) | (std::uint64_t{carry} << (32-bits)));
        carry = here;
    }
    out.size = x.size;
    out.trim();
}

inline bool increment(Unsigned& x) noexcept {
    for (std::size_t i = 0; i < x.size; ++i) {
        if (++x.limb[i] != 0) return true;
    }
    if (x.size == capacity) return false;
    x.limb[x.size++] = 1;
    return true;
}

// Normalized base-2^32 long division (Knuth algorithm D).
// A 64/32 quotient handles the single-limb denominator case.
// Uses only standard uint32_t/uint64_t and stack storage.
inline bool divrem(Unsigned const& dividend, Unsigned const& divisor,
                   Unsigned& quotient, Unsigned& remainder) noexcept {
    quotient = {};
    remainder = {};
    if (divisor.size == 0) return false;
    if (compare(dividend, divisor) < 0) {
        remainder = dividend;
        return true;
    }
    if (divisor.size == 1) {
        quotient.size = dividend.size;
        std::uint64_t r = 0;
        for (std::size_t i = dividend.size; i != 0; --i) {
            const std::uint64_t partial = (r << 32) | dividend.limb[i-1];
            quotient.limb[i-1] =
                static_cast<std::uint32_t>(partial / divisor.limb[0]);
            r = partial % divisor.limb[0];
        }
        quotient.trim();
        if (r != 0) {
            remainder.limb[0] = static_cast<std::uint32_t>(r);
            remainder.size = 1;
        }
        return true;
    }
    const unsigned normalization =
        std::countl_zero(divisor.limb[divisor.size-1]);
    Unsigned v, u;
    if (!shift_left(divisor, normalization, v) ||
        !shift_left(dividend, normalization, u)) return false;
    const std::size_t n = v.size;
    const std::size_t m = u.size-n;
    if (u.size == capacity) return false; // sentinel top word
    u.limb[u.size] = 0;
    quotient.size = m+1;

    for (std::size_t offset = m+1; offset != 0; --offset) {
        const std::size_t j = offset-1;
        const std::uint64_t top = u.limb[j+n];
        const std::uint64_t below = u.limb[j+n-1];
        std::uint64_t estimate, estimated_remainder;
        if (top >= v.limb[n-1]) {
            estimate = base-1;
            estimated_remainder = below + v.limb[n-1];
        } else {
            const std::uint64_t numerator = (top << 32) | below;
            estimate = numerator / v.limb[n-1];
            estimated_remainder = numerator % v.limb[n-1];
        }
        while (estimated_remainder < base &&
               estimate * v.limb[n-2] >
               estimated_remainder * base + u.limb[j+n-2]) {
            --estimate;
            estimated_remainder += v.limb[n-1];
        }
        std::uint64_t carry = 0, borrow = 0;
        for (std::size_t i = 0; i < n; ++i) {
            const std::uint64_t product = estimate * v.limb[i] + carry;
            carry = product >> 32;
            const std::uint64_t subtrahend = (product & (base-1)) + borrow;
            borrow = u.limb[j+i] < subtrahend;
            u.limb[j+i] = static_cast<std::uint32_t>(
                std::uint64_t{u.limb[j+i]} - subtrahend);
        }
        const std::uint64_t final_subtrahend = carry+borrow;
        const bool overestimate = u.limb[j+n] < final_subtrahend;
        u.limb[j+n] = static_cast<std::uint32_t>(
            std::uint64_t{u.limb[j+n]} - final_subtrahend);
        if (overestimate) {
            --estimate;
            std::uint64_t addition_carry = 0;
            for (std::size_t i = 0; i < n; ++i) {
                const std::uint64_t sum = std::uint64_t{u.limb[j+i]} +
                                          v.limb[i] + addition_carry;
                u.limb[j+i] = static_cast<std::uint32_t>(sum);
                addition_carry = sum >> 32;
            }
            u.limb[j+n] += static_cast<std::uint32_t>(addition_carry);
        }
        quotient.limb[j] = static_cast<std::uint32_t>(estimate);
    }
    quotient.trim();
    Unsigned normalized_remainder;
    normalized_remainder.size = n;
    for (std::size_t i = 0; i < n; ++i)
        normalized_remainder.limb[i] = u.limb[i];
    normalized_remainder.trim();
    shift_right(normalized_remainder, normalization, remainder);
    return true;
}

} // namespace Ariadne::FloatAPSmallDivision
#endif
