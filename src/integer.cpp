/***************************************************************************
 *            numeric/integer.cpp
 *
 *  Copyright  2013-22  Pieter Collins
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

#include "utility/stdlib.hpp"

#include "integer.hpp"

#include "macro/general.hpp"
#include "utility/string.hpp"
#include "foundation/logical.hpp"

#include "numeric/concepts.hpp"

#include <bit>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace Ariadne {

Comparison cmp(Integer const& z1, Integer const& z2);

Integer::~Integer() = default;

Integer::Integer(BigUnsigned&& magnitude,bool negative) noexcept
    : _magnitude(std::move(magnitude)),
      _negative(negative && !_magnitude.is_zero())
{
}

Integer Integer::from_magnitude(BigUnsigned magnitude,bool negative) noexcept {
    return Integer(std::move(magnitude),negative);
}

Integer::Integer() = default;

Integer::Integer(Nat32 m) : _magnitude(static_cast<std::uint64_t>(m.get_ui())) {
}

Integer::Integer(Int32 n) : Integer(Int64(n.get_si())) {
}

Integer::Integer(Nat64 m) : _magnitude(m.get_ui()) {
}

Integer::Integer(Int64 n) {
    const std::int64_t value=n.get_si();
    _negative=value<0;
    const std::uint64_t magnitude=_negative
        ? static_cast<std::uint64_t>(-(value+1))+std::uint64_t{1}
        : static_cast<std::uint64_t>(value);
    _magnitude=BigUnsigned(magnitude);
}

Integer::Integer(String const& str) {
    std::string_view digits(str.c_str(),str.size());
    if (digits.empty()) {
        throw std::invalid_argument("Integer empty decimal literal");
    }
    const bool negative=digits.front()=='-';
    if (negative || digits.front()=='+') {
        digits.remove_prefix(1u);
    }
    _magnitude=from_decimal(digits);
    _negative=negative && !_magnitude.is_zero();
}

Integer::Integer(const Integer& z) = default;
Integer::Integer(Integer&& z) noexcept
    : _magnitude(std::move(z._magnitude)),_negative(z._negative)
{
    z._negative=false;
}
Integer& Integer::operator=(const Integer& z) = default;
Integer& Integer::operator=(Integer&& z) noexcept {
    if (this!=&z) {
        _magnitude=std::move(z._magnitude);
        _negative=z._negative;
        z._negative=false;
    }
    return *this;
}

long int Integer::get_si() const {
    return get<long int>();
}

Integer& operator++(Integer& z) {
    z=add(z,Integer(1));
    return z;
}

Integer& operator--(Integer& z) {
    z=sub(z,Integer(1));
    return z;
}

Integer& operator+=(Integer& z1, Integer const& z2) {
    z1=add(z1,z2);
    return z1;
}

Integer& operator-=(Integer& z1, Integer const& z2) {
    z1=sub(z1,z2);
    return z1;
}

Integer& operator*=(Integer& z1, Integer const& z2) {
    z1=mul(z1,z2);
    return z1;
}

Integer nul(Integer const&) {
    return Integer{};
}

Integer pos(Integer const& z) {
    return z;
}

Integer neg(Integer const& z) {
    return Integer::from_magnitude(z._magnitude,!z._negative);
}

Natural sqr(Integer const& z) {
    return Natural(Integer::from_magnitude(
        square(z._magnitude),false));
}

Integer add(Integer const& z1, Integer const& z2) {
    if (z1._negative==z2._negative) {
        return Integer::from_magnitude(
            add(z1._magnitude,z2._magnitude),z1._negative);
    }
    const int comparison=compare(z1._magnitude,z2._magnitude);
    if (comparison==0) {
        return Integer{};
    }
    if (comparison>0) {
        return Integer::from_magnitude(
            sub(z1._magnitude,z2._magnitude),z1._negative);
    }
    return Integer::from_magnitude(
        sub(z2._magnitude,z1._magnitude),z2._negative);
}

Integer sub(Integer const& z1, Integer const& z2) {
    // B41: a-b is performed directly from the signed magnitudes.
    // The former add(a,neg(b)) materialized a full copy of b before
    // determining the arithmetic path. Opposite signs add magnitudes,
    // equal signs subtract their absolute values in the larger-first
    // order; no temporary negated operand is needed.
    if (z1._negative!=z2._negative) {
        return Integer::from_magnitude(
            add(z1._magnitude,z2._magnitude),z1._negative);
    }
    const int comparison=compare(z1._magnitude,z2._magnitude);
    if (comparison==0) {
        return Integer{};
    }
    if (comparison>0) {
        return Integer::from_magnitude(
            sub(z1._magnitude,z2._magnitude),z1._negative);
    }
    return Integer::from_magnitude(
        sub(z2._magnitude,z1._magnitude),!z1._negative);
}

Integer mul(Integer const& z1, Integer const& z2) {
#if defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED) && \
    defined(__x86_64__) && defined(__SIZEOF_INT128__) && \
    defined(__GNUC__) && !defined(__clang__) && !defined(_MSC_VER)
    // B33: the exact same Integer object is multiplied by itself.
    // Dispatch to B32's symmetric kernel without inspecting operand
    // values or changing the non-aliased multiplication path. Keep
    // other architectures and smaller magnitudes unchanged.
    if (&z1==&z2 && z1._magnitude.size()>=8u) {
        return Integer::from_magnitude(square(z1._magnitude),false);
    }
#endif
#if defined(ARIADNE_BIG_UNSIGNED_ARM64) && defined(__aarch64__) && \
    defined(__SIZEOF_INT128__) && (defined(__GNUC__) || defined(__clang__))
    // B37: on native ARM64, only actual self-aliasing uses the new
    // symmetric square; equal-valued distinct objects still use mul.
    // Size-first keeps the already-fast 1..7-limb path unchanged.
    if (z1._magnitude.size()>=8u && &z1==&z2) {
        return Integer::from_magnitude(square(z1._magnitude),false);
    }
#endif
    return Integer::from_magnitude(
        mul(z1._magnitude,z2._magnitude),z1._negative!=z2._negative);
}

// Signed integer division truncates toward zero; the remainder
// has the dividend's sign. Both agree with GMP's mpz_tdiv_* family.
Integer quot(Integer const& z1, Integer const& z2) {
    return Integer::from_magnitude(quot(z1._magnitude,z2._magnitude),
                                   z1._negative!=z2._negative);
}

Integer rem(Integer const& z1, Integer const& z2) {
    return Integer::from_magnitude(rem(z1._magnitude,z2._magnitude),
                                   z1._negative);
}

Integer pow(Integer const& z, Nat exponent) {
#if defined(ARIADNE_BIG_UNSIGNED_GCC_FUSED) && \
    defined(__x86_64__) && defined(__SIZEOF_INT128__) && \
    defined(__GNUC__) && !defined(__clang__) && !defined(_MSC_VER)
    // B66: exact, generic MSB-first exponentiation by squaring.
    // B65 GMP-verified this schedule for powers other than three and
    // found lower owning-result latency on the user's GCC x86-64 host
    // at all six official sizes. Borrow the original base until a
    // square is necessary: do not copy it merely to start the loop.
    const auto n=static_cast<std::uint64_t>(exponent);
    const bool negative=z._negative && (n&std::uint64_t{1})!=0u;
    if (n==0u) {
        return Integer::from_magnitude(
            BigUnsigned(std::uint64_t{1}),false);
    }
    if (n==1u) {
        return Integer::from_magnitude(z._magnitude,negative);
    }
    const BigUnsigned& base=z._magnitude;
    const BigUnsigned* current=&base;
    BigUnsigned value;
    for (std::uint64_t bit=std::bit_floor(n)>>1u;
            bit!=0u;bit>>=1u) {
        // current aliases the input initially and value thereafter;
        // square/mul fully materialize their results before assignment.
        value=square(*current);
        current=&value;
        if ((n&bit)!=0u) {
            value=mul(value,base);
        }
    }
    return Integer::from_magnitude(std::move(value),negative);
#else
    // Preserve B32's right-to-left algorithm for ARM64, Clang, MSVC
    // and other builds until native performance can be validated.
    BigUnsigned value(std::uint64_t{1});
    BigUnsigned base=z._magnitude;
    auto n=static_cast<std::uint64_t>(exponent);
    const bool negative=z._negative && (n&std::uint64_t{1})!=0u;
    bool has_factor=false;
    while (n!=0u) {
        if ((n&std::uint64_t{1})!=0u) {
            // B32: multiplying the initial accumulator 1 by the first
            // selected power is unnecessary. Copy that power instead;
            // keep the same squaring and multiplication schedule.
            if (!has_factor) {
                value=base;
                has_factor=true;
            } else {
                value=mul(value,base);
            }
        }
        n>>=1u;
        if (n!=0u) {
            base=square(base);
        }
    }
    return Integer::from_magnitude(std::move(value),negative);
#endif
}

Integer min(Integer const& z1,Integer const& z2) {
    return (z1<z2)?z1:z2;
}

Integer max(Integer const& z1,Integer const& z2) {
    return (z1>z2)?z1:z2;
}

Natural abs(Integer const& z) {
    return Natural(Integer::from_magnitude(z._magnitude,false));
}

Bool is_nan(Integer const&) { return false; }
Bool is_inf(Integer const&) { return false; }
Bool is_finite(Integer const&) { return true; }
Bool is_zero(Integer const& z) { return z._magnitude.is_zero(); }

Sign sgn(Integer const& z) {
    return z._magnitude.is_zero() ? Sign::ZERO
        : z._negative ? Sign::NEGATIVE : Sign::POSITIVE;
}

Comparison cmp(Integer const& z1, Integer const& z2) {
    if (z1._negative!=z2._negative) {
        return z1._negative ? Comparison::LESS : Comparison::GREATER;
    }
    int c=compare(z1._magnitude,z2._magnitude);
    if (z1._negative) {
        c=-c;
    }
    return c==0 ? Comparison::EQUAL
        : (c>0 ? Comparison::GREATER : Comparison::LESS);
}

Boolean eq(Integer const& z1, Integer const& z2) {
    return cmp(z1,z2)==Comparison::EQUAL;
}

Boolean lt(Integer const& z1, Integer const& z2) {
    return cmp(z1,z2)==Comparison::LESS;
}

String Integer::literal() const {
    String result;
    // String's std::string constructor copies even an rvalue. Move
    // the formatter's owned buffer into its standard-string base.
    static_cast<std::string&>(result)=to_decimal(_magnitude);
    if (_negative) {
        result.insert(result.begin(),'-');
    }
    return result;
}

OutputStream& operator<<(OutputStream& os, Integer const& z) {
    return os << z.literal();
}

Integer operator""_z(unsigned long long int n) {
    return Integer(Nat64(n));
}

Integer operator""_z(const char* str, std::size_t) {
    return Integer(String(str));
}


template<> String class_name<uint>() { return "uint"; }

template<> String class_name<int>() { return "int"; }

template<> String class_name<Integer>() { return "Integer"; }

template<> String class_name<Natural>() { return "Natural"; }

Int log2floor(Natural const& z) {
    if (z==0) { return -1; }
    return static_cast<Int>(z.magnitude().bit_length()-1u);
}

OutputStream& operator<<(OutputStream& os, Sign s) {
    return os << ( (s==Sign::ZERO) ? "ZERO" : (s==Sign::NEGATIVE) ? "NEGATIVE" : "POSITIVE" );
}

static_assert(Ring<Integer>);
static_assert(OrderedLattice<Integer>);


} // namespace Ariadne
