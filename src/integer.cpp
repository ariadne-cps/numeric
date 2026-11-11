/***************************************************************************
 *            integer.cpp
 *
 *  Copyright  2013-22  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Numeric.
 *
 *  Ariadne Numeric is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Numeric is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "integer.hpp"

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "foundation/logical.hpp"
#include "numeric/concepts.hpp"
#include "utility/macros.hpp"
#include "utility/stdlib.hpp"
#include "utility/string.hpp"

namespace Ariadne {

Comparison cmp(const Integer& z1, const Integer& z2);

Integer::Integer() = default;

Integer::Integer(Nat32 m) : _magnitude(static_cast<std::uint64_t>(m.get_ui())) {
}

Integer::Integer(Int32 n) : Integer(Int64(n.get_si())) {
}

Integer::Integer(Nat64 m) : _magnitude(m.get_ui()) {
}

Integer::Integer(Int64 n) {
    const std::int64_t value = n.get_si();
    _negative = value < 0;
    const std::uint64_t magnitude = _negative
                                      ? static_cast<std::uint64_t>(-(value + 1)) + std::uint64_t{1}
                                      : static_cast<std::uint64_t>(value);
    _magnitude = BigUnsigned(magnitude);
}

Integer::Integer(const String& str) {
    std::string_view digits(str.c_str(), str.size());
    if (digits.empty()) {
        throw std::invalid_argument("Integer empty decimal literal");
    }
    const bool negative = digits.front() == '-';
    if (negative || digits.front() == '+') {
        digits.remove_prefix(1u);
    }
    _magnitude = from_decimal(digits);
    _negative = negative && !_magnitude.is_zero();
}

long int Integer::get_si() const {
    return get<long int>();
}

Integer& operator++(Integer& z) {
    z = add(z, Integer(1));
    return z;
}

Integer& operator--(Integer& z) {
    z = sub(z, Integer(1));
    return z;
}

Integer& operator+=(Integer& z1, const Integer& z2) {
    z1 = add(z1, z2);
    return z1;
}

Integer& operator-=(Integer& z1, const Integer& z2) {
    z1 = sub(z1, z2);
    return z1;
}

Integer& operator*=(Integer& z1, const Integer& z2) {
    z1 = mul(z1, z2);
    return z1;
}

Integer nul(const Integer&) {
    return Integer{};
}

Integer pos(const Integer& z) {
    return z;
}

Integer neg(const Integer& z) {
    return Integer::from_magnitude(z._magnitude, !z._negative);
}

Natural sqr(const Integer& z) {
    return Natural(Integer::from_magnitude(
        square(z._magnitude), false));
}

Integer add(const Integer& z1, const Integer& z2) {
    if (z1._negative == z2._negative) {
        return Integer::from_magnitude(
            add(z1._magnitude, z2._magnitude), z1._negative);
    }
    const int comparison = compare(z1._magnitude, z2._magnitude);
    if (comparison == 0) {
        return Integer{};
    }
    if (comparison > 0) {
        return Integer::from_magnitude(
            sub(z1._magnitude, z2._magnitude), z1._negative);
    }
    return Integer::from_magnitude(
        sub(z2._magnitude, z1._magnitude), z2._negative);
}

Integer sub(const Integer& z1, const Integer& z2) {
    // B41: a-b is performed directly from the signed magnitudes.
    // The former add(a,neg(b)) materialized a full copy of b before
    // determining the arithmetic path. Opposite signs add magnitudes,
    // equal signs subtract their absolute values in the larger-first
    // order; no temporary negated operand is needed.
    if (z1._negative != z2._negative) {
        return Integer::from_magnitude(
            add(z1._magnitude, z2._magnitude), z1._negative);
    }
    const int comparison = compare(z1._magnitude, z2._magnitude);
    if (comparison == 0) {
        return Integer{};
    }
    if (comparison > 0) {
        return Integer::from_magnitude(
            sub(z1._magnitude, z2._magnitude), z1._negative);
    }
    return Integer::from_magnitude(
        sub(z2._magnitude, z1._magnitude), !z1._negative);
}

Integer mul(const Integer& z1, const Integer& z2) {
    // The unsigned layer owns alias-aware square/kernel selection.
    // Only the product sign belongs to the signed front end.
    return Integer::from_magnitude(
        mul(z1._magnitude, z2._magnitude), z1._negative != z2._negative);
}

// Signed integer division truncates toward zero; the remainder
// has the dividend's sign. Both agree with GMP's mpz_tdiv_* family.
Integer quot(const Integer& z1, const Integer& z2) {
    return Integer::from_magnitude(quot(z1._magnitude, z2._magnitude),
                                   z1._negative != z2._negative);
}

Integer rem(const Integer& z1, const Integer& z2) {
    return Integer::from_magnitude(rem(z1._magnitude, z2._magnitude),
                                   z1._negative);
}

Integer isqrt(const Integer& z) {
    if (z._negative) {
        throw std::domain_error("Integer square root of negative value");
    }
    return Integer::from_magnitude(isqrt(z._magnitude), false);
}

std::pair<Integer, Integer> sqrtrem(const Integer& z) {
    if (z._negative) {
        throw std::domain_error("Integer square root of negative value");
    }
    auto [root, remainder] = sqrtrem(z._magnitude);
    return {
        Integer::from_magnitude(std::move(root), false),
        Integer::from_magnitude(std::move(remainder), false)};
}

Integer pow(const Integer& z, Nat exponent) {
    // Only sign parity and negative-zero normalization belong here.
    // The unsigned engine selects the existing per-platform schedule.
    const auto n = static_cast<std::uint64_t>(exponent);
    const bool negative = z._negative && (n & std::uint64_t{1}) != 0u;
    return Integer::from_magnitude(pow(z._magnitude, n), negative);
}

Integer min(const Integer& z1, const Integer& z2) {
    return (z1 < z2) ? z1 : z2;
}

Integer max(const Integer& z1, const Integer& z2) {
    return (z1 > z2) ? z1 : z2;
}

Natural abs(const Integer& z) {
    return Natural(Integer::from_magnitude(z._magnitude, false));
}

Bool is_nan(const Integer&) { return false; }
Bool is_inf(const Integer&) { return false; }
Bool is_finite(const Integer&) { return true; }
Bool is_zero(const Integer& z) { return z._magnitude.is_zero(); }

Sign sgn(const Integer& z) {
    return z._magnitude.is_zero() ? Sign::ZERO
         : z._negative            ? Sign::NEGATIVE
                                  : Sign::POSITIVE;
}

Comparison cmp(const Integer& z1, const Integer& z2) {
    if (z1._negative != z2._negative) {
        return z1._negative ? Comparison::LESS : Comparison::GREATER;
    }
    int c = compare(z1._magnitude, z2._magnitude);
    if (z1._negative) {
        c = -c;
    }
    return c == 0 ? Comparison::EQUAL
                  : (c > 0 ? Comparison::GREATER : Comparison::LESS);
}

Boolean eq(const Integer& z1, const Integer& z2) {
    return cmp(z1, z2) == Comparison::EQUAL;
}

Boolean lt(const Integer& z1, const Integer& z2) {
    return cmp(z1, z2) == Comparison::LESS;
}

String Integer::literal() const {
    String result;
    // String's std::string constructor copies even an rvalue. Move
    // the formatter's owned buffer into its standard-string base.
    static_cast<std::string&>(result) = to_decimal(_magnitude);
    if (_negative) {
        result.insert(result.begin(), '-');
    }
    return result;
}

OutputStream& operator<<(OutputStream& os, const Integer& z) {
    return os << z.literal();
}

Integer operator""_z(unsigned long long int n) {
    return Integer(Nat64(n));
}

Integer operator""_z(const char* str, std::size_t) {
    return Integer(String(str));
}

template<>
String class_name<uint>() { return "uint"; }

template<>
String class_name<int>() { return "int"; }

template<>
String class_name<Integer>() { return "Integer"; }

template<>
String class_name<Natural>() { return "Natural"; }

Int log2floor(const Natural& z) {
    if (z == 0) {
        return -1;
    }
    return static_cast<Int>(z.magnitude().bit_length() - 1u);
}

OutputStream& operator<<(OutputStream& os, Sign s) {
    return os << ((s == Sign::ZERO) ? "ZERO" : (s == Sign::NEGATIVE) ? "NEGATIVE"
                                                                     : "POSITIVE");
}

static_assert(Ring<Integer>);
static_assert(OrderedLattice<Integer>);

} // namespace Ariadne
