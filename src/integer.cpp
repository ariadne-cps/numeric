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

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace Ariadne {

Comparison cmp(Integer const& z1, Integer const& z2);

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
    // The unsigned layer owns alias-aware square/kernel selection.
    // Only the product sign belongs to the signed front end.
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
    // Only sign parity and negative-zero normalization belong here.
    // The unsigned engine selects the existing per-platform schedule.
    const auto n=static_cast<std::uint64_t>(exponent);
    const bool negative=z._negative && (n&std::uint64_t{1})!=0u;
    return Integer::from_magnitude(pow(z._magnitude,n),negative);
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
