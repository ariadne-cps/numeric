/***************************************************************************
 *            numeric/detail/checked_arithmetic.hpp
 *
 *  Checked scalar arithmetic for working precisions and binary shifts.
 *  These operations validate metadata before any storage is allocated.
 ****************************************************************************/
#ifndef ARIADNE_NUMERIC_DETAIL_CHECKED_ARITHMETIC_HPP
#define ARIADNE_NUMERIC_DETAIL_CHECKED_ARITHMETIC_HPP

#include <concepts>
#include <limits>
#include <stdexcept>

namespace Ariadne::Detail {

template<std::unsigned_integral T>
constexpr T checked_unsigned_add(T left,T right,char const* message) {
    if(left>std::numeric_limits<T>::max()-right) throw std::overflow_error(message);
    return left+right;
}

template<std::unsigned_integral To,std::unsigned_integral From>
constexpr To checked_unsigned_cast(From value,char const* message) {
    if constexpr(std::numeric_limits<To>::digits<std::numeric_limits<From>::digits) {
        if(value>static_cast<From>(std::numeric_limits<To>::max()))
            throw std::overflow_error(message);
    }
    return static_cast<To>(value);
}

} // namespace Ariadne::Detail
#endif
