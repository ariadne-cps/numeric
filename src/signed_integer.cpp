/***************************************************************************
 *            numeric/signed_integer.cpp
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

/*! \file numeric/signed_integer.cpp
 *  \brief GMP-independent signed magnitude arithmetic.
 */

#include "numeric/signed_integer.hpp"

#include <limits>
#include <stdexcept>
#include <utility>

namespace Ariadne {

SignedInteger::SignedInteger(BigUnsigned magnitude,bool negative) noexcept
    : _magnitude(std::move(magnitude)),
      _negative(negative && !_magnitude.is_zero())
{
}

SignedInteger::SignedInteger(std::int64_t value) noexcept
    : _magnitude(value<0
        ? static_cast<std::uint64_t>(-(value+1))+std::uint64_t{1}
        : static_cast<std::uint64_t>(value)),
      _negative(value<0)
{
}

SignedInteger::SignedInteger(std::uint64_t value) noexcept
    : _magnitude(value)
{
}

SignedInteger::SignedInteger(std::string_view decimal)
{
    if (decimal.empty()) {
        throw std::invalid_argument("SignedInteger empty decimal input");
    }
    const bool negative=decimal.front()=='-';
    if (negative || decimal.front()=='+') {
        decimal.remove_prefix(1u);
    }
    _magnitude=from_decimal(decimal);
    _negative=negative && !_magnitude.is_zero();
}

SignedInteger SignedInteger::from_magnitude(
        BigUnsigned magnitude,bool negative) noexcept
{
    return SignedInteger(std::move(magnitude),negative);
}

int SignedInteger::sign() const noexcept
{
    return is_zero()?0:(_negative?-1:1);
}

std::string SignedInteger::decimal() const
{
    const std::string value=to_decimal(_magnitude);
    return _negative?"-"+value:value;
}

std::uint64_t SignedInteger::get_u64() const
{
    if (_negative || _magnitude.size()>1u) {
        throw std::out_of_range("SignedInteger does not fit uint64_t");
    }
    return _magnitude.limb(0u);
}

std::int64_t SignedInteger::get_i64() const
{
    if (_magnitude.size()>1u) {
        throw std::out_of_range("SignedInteger does not fit int64_t");
    }
    const std::uint64_t value=_magnitude.limb(0u);
    const auto signed_max=static_cast<std::uint64_t>(
        std::numeric_limits<std::int64_t>::max());
    if (_negative) {
        if (value>signed_max+std::uint64_t{1}) {
            throw std::out_of_range("SignedInteger does not fit int64_t");
        }
        if (value==signed_max+std::uint64_t{1}) {
            return std::numeric_limits<std::int64_t>::min();
        }
        return -static_cast<std::int64_t>(value);
    }
    if (value>signed_max) {
        throw std::out_of_range("SignedInteger does not fit int64_t");
    }
    return static_cast<std::int64_t>(value);
}

int compare(SignedInteger const& a,SignedInteger const& b) noexcept
{
    if (a._negative!=b._negative) {
        return a._negative?-1:1;
    }
    const int magnitude_order=compare(a._magnitude,b._magnitude);
    return a._negative?-magnitude_order:magnitude_order;
}

SignedInteger neg(SignedInteger const& a)
{
    return SignedInteger::from_magnitude(a._magnitude,!a._negative);
}

SignedInteger add(SignedInteger const& a,SignedInteger const& b)
{
    if (a._negative==b._negative) {
        return SignedInteger::from_magnitude(
            add(a._magnitude,b._magnitude),a._negative);
    }
    const int magnitudes=compare(a._magnitude,b._magnitude);
    if (magnitudes==0) {
        return SignedInteger{};
    }
    if (magnitudes>0) {
        return SignedInteger::from_magnitude(
            sub(a._magnitude,b._magnitude),a._negative);
    }
    return SignedInteger::from_magnitude(
        sub(b._magnitude,a._magnitude),b._negative);
}

SignedInteger sub(SignedInteger const& a,SignedInteger const& b)
{
    return add(a,neg(b));
}

SignedInteger mul(SignedInteger const& a,SignedInteger const& b)
{
    return SignedInteger::from_magnitude(
        mul(a._magnitude,b._magnitude),a._negative!=b._negative);
}

SignedInteger square(SignedInteger const& a)
{
    return SignedInteger::from_magnitude(square(a._magnitude),false);
}

SignedInteger pow(SignedInteger const& a,std::uint64_t exponent)
{
    BigUnsigned result(std::uint64_t{1});
    BigUnsigned base=a._magnitude;
    const bool negative=a._negative && (exponent&std::uint64_t{1})!=0u;
    while (exponent!=0u) {
        if ((exponent&std::uint64_t{1})!=0u) {
            result=mul(result,base);
        }
        exponent>>=1u;
        if (exponent!=0u) {
            base=square(base);
        }
    }
    return SignedInteger::from_magnitude(std::move(result),negative);
}

std::pair<SignedInteger,SignedInteger> divrem(
        SignedInteger const& a,SignedInteger const& b)
{
    auto [quotient,remainder]=divrem(a._magnitude,b._magnitude);
    return {
        SignedInteger::from_magnitude(
            std::move(quotient),a._negative!=b._negative),
        SignedInteger::from_magnitude(std::move(remainder),a._negative)
    };
}

SignedInteger shift_left(SignedInteger const& a,std::size_t bits)
{
    return SignedInteger::from_magnitude(
        shift_left(a._magnitude,bits),a._negative);
}

SignedInteger shift_right(SignedInteger const& a,std::size_t bits)
{
    return SignedInteger::from_magnitude(
        shift_right(a._magnitude,bits),a._negative);
}

} // namespace Ariadne
