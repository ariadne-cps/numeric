/***************************************************************************
 *            numeric/float_operations.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/** @file numeric/float_operations.hpp
 *  @brief Common header for floating-point number operations.
 */

#ifndef ARIADNE_FLOAT_OPERATIONS_HPP
#define ARIADNE_FLOAT_OPERATIONS_HPP

#include "arithmetic.hpp"
#include "dyadic.hpp"
#include "float.decl.hpp"
#include "number.decl.hpp"
#include "operators.hpp"
#include "rational.hpp"
#include "utility/macros.hpp"

namespace Ariadne {

/** @brief Declare float operations for the supported operand types. */
template<class X, class R = X>
struct DeclareFloatOperations
    : DeclareRealOperations<X, R>,
      DeclareComparisonOperations<X, LessTrait<X>, EqualsTrait<X>> {
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream&, const X&);
    /** @brief Read the value from the input stream. */
    friend InputStream& operator>>(InputStream&, X&);
};

/** @brief Declare directed float operations for the supported operand types. */
template<class X, class NX = OppositeTrait<X>>
struct DeclareDirectedFloatOperations
    : DeclareDirectedNumericOperations<X, NX>,
      DeclareMixedDirectedGroupOperators<X, NX, GenericType<X>, GenericType<NX>> {
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream&, const X&);
    /** @brief Read the value from the input stream. */
    friend InputStream& operator>>(InputStream&, X&);
};

/** @brief Declare positive directed float operations for the supported operand types. */
template<class PX, class QPX = OppositeTrait<PX>>
struct DeclarePositiveDirectedFloatOperations
    : DeclarePositiveDirectedNumericOperations<PX, QPX> {
};

/** @brief Declare positive float operations for the supported operand types. */
template<class PX>
struct DeclarePositiveFloatOperations
    : DeclarePositiveDirectedFloatOperations<PX, PX> {
};

/** @brief Implement float operations for the supported operand types. */
template<class X, class R = X>
struct DefineFloatOperations
    : public DefineFieldOperators<X, R>,
      public DefineComparisonOperators<X, LessTrait<X>, EqualsTrait<X>>,
      public DefineConcreteGenericOperators<X>,
      public DeclareFloatOperations<X> {
};

/** @brief Implement directed float operations for the supported operand types. */
template<class X, class NX>
struct DefineDirectedFloatOperations
    : public DefineDirectedGroupOperators<X, NX>,
      public DefineDirectedGroupOperators<NX, X>,
      public DefineDirectedComparisonOperators<X, NX, LessTrait<X>, EqualsTrait<X>>,
      public DefineDirectedComparisonOperators<NX, X, LessTrait<NX>, EqualsTrait<NX>>,
      public DefineConcreteGenericOperators<X> {
};

/** @brief Dispatch float operations for the supported operand types. */
template<class X, class R = X>
struct DispatchFloatOperations
    : DispatchNumericOperations<X, R>,
      DispatchComparisonOperations<X, LessTrait<X>, EqualsTrait<X>>,
      DefineConcreteGenericOperators<X> {
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const X& x) { return Operations<X>::_write(os, x); }
    /** @brief Read the value from the input stream. */
    friend InputStream& operator>>(InputStream& is, X& x) { return Operations<X>::_read(is, x); }
};

/** @brief Dispatch directed float operations for the supported operand types. */
template<class X>
struct DispatchDirectedFloatOperations
    : DispatchDirectedNumericOperations<X, OppositeTrait<X>>,
      DispatchDirectedNumericOperations<OppositeTrait<X>, X>,
      DispatchDirectedComparisonOperations<X, OppositeTrait<X>, LessTrait<X>, EqualsTrait<X>>,
      DispatchDirectedComparisonOperations<OppositeTrait<X>, X, LessTrait<OppositeTrait<X>>, EqualsTrait<OppositeTrait<X>>>,
      DefineConcreteGenericOperators<X> {
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const X& x) { return Operations<X>::_write(os, x); }
    /** @brief Read the value from the input stream. */
    friend InputStream& operator>>(InputStream& is, X& x) { return Operations<X>::_read(is, x); }
};

/** @brief Dispatch positive directed float operations for the supported operand types. */
template<class PX, class QPX = OppositeTrait<PX>>
struct DispatchPositiveDirectedFloatOperations
    : DispatchPositiveDirectedNumericOperations<PX, QPX>,
      DispatchPositiveDirectedNumericOperations<QPX, PX>,
      DefineConcreteGenericOperators<PX> {
};

/** @brief Dispatch positive float operations for the supported operand types. */
template<class PX>
struct DispatchPositiveFloatOperations
    : DispatchPositiveDirectedNumericOperations<PX, PX>,
      DefineConcreteGenericOperators<PX> {
};

template<class Y1, class Y2>
concept HaveCmp = requires(Y1 y1, Y2 y2) {
    { cmp(y1, y2) } -> SameAs<Comparison>;
};
/** @brief Compare the operands using the comparison result type. */
template<class OP, class Y1, class Y2>
    requires HaveCmp<Y1, Y2>
Boolean _cmp(OP op, Y1 const& y1, Y2 const& y2) {
    const Comparison comparison = cmp(y1, y2);
    // INCOMPARABLE is encoded as -128, not as an ordered value.
    // Comparing that sentinel with EQUAL (=0) makes NaN spuriously less
    // than ordinary numbers. Preserve the same unordered behavior as independent exact arithmetic.
    if (comparison == Comparison::INCOMPARABLE) {
        return SameAs<OP, Ne>;
    }
    return op(comparison, Comparison::EQUAL);
}

/** @brief Test equality of the operands. */
template<class Y1, class Y2>
    requires HaveCmp<Y1, Y2>
Boolean operator==(Y1 const& y1, Y2 const& y2) { return _cmp(Eq(), y1, y2); }
/** @brief Test inequality of the operands. */
template<class Y1, class Y2>
    requires HaveCmp<Y1, Y2>
Boolean operator!=(Y1 const& y1, Y2 const& y2) { return _cmp(Ne(), y1, y2); }
/** @brief Test whether the first operand is less than the second. */
template<class Y1, class Y2>
    requires HaveCmp<Y1, Y2>
Boolean operator<(Y1 const& y1, Y2 const& y2) { return _cmp(Lt(), y1, y2); }
/** @brief Test whether the first operand is greater than the second. */
template<class Y1, class Y2>
    requires HaveCmp<Y1, Y2>
Boolean operator>(Y1 const& y1, Y2 const& y2) { return _cmp(Gt(), y1, y2); }
/** @brief Test whether the first operand is less than or equal to the second. */
template<class Y1, class Y2>
    requires HaveCmp<Y1, Y2>
Boolean operator<=(Y1 const& y1, Y2 const& y2) { return _cmp(Le(), y1, y2); }
/** @brief Test whether the first operand is greater than or equal to the second. */
template<class Y1, class Y2>
    requires HaveCmp<Y1, Y2>
Boolean operator>=(Y1 const& y1, Y2 const& y2) { return _cmp(Ge(), y1, y2); }

// Mixed Float-ValidatedNumber comparisons
/** @brief Test equality of the operands. */
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator==(const Y& y1, const F& x2) { return Bounds<F>(y1, x2.precision()) == x2; }
/** @brief Test inequality of the operands. */
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator!=(const Y& y1, const F& x2) { return Bounds<F>(y1, x2.precision()) != x2; }
/** @brief Test whether the first operand is less than the second. */
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator<(const Y& y1, const F& x2) { return Bounds<F>(y1, x2.precision()) < x2; }
/** @brief Test whether the first operand is greater than the second. */
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator>(const Y& y1, const F& x2) { return Bounds<F>(y1, x2.precision()) > x2; }
/** @brief Test whether the first operand is less than or equal to the second. */
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator<=(const Y& y1, const F& x2) { return Bounds<F>(y1, x2.precision()) <= x2; }
/** @brief Test whether the first operand is greater than or equal to the second. */
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator>=(const Y& y1, const F& x2) { return Bounds<F>(y1, x2.precision()) >= x2; }
/** @brief Test equality of the operands. */
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator==(const F& x1, const Y& y2) { return x1 == Bounds<F>(y2, x1.precision()); }
/** @brief Test inequality of the operands. */
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator!=(const F& x1, const Y& y2) { return x1 != Bounds<F>(y2, x1.precision()); }
/** @brief Test whether the first operand is less than the second. */
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator<(const F& x1, const Y& y2) { return x1 < Bounds<F>(y2, x1.precision()); }
/** @brief Test whether the first operand is greater than the second. */
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator>(const F& x1, const Y& y2) { return x1 > Bounds<F>(y2, x1.precision()); }
/** @brief Test whether the first operand is less than or equal to the second. */
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator<=(const F& x1, const Y& y2) { return x1 <= Bounds<F>(y2, x1.precision()); }
/** @brief Test whether the first operand is greater than or equal to the second. */
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator>=(const F& x1, const Y& y2) { return x1 >= Bounds<F>(y2, x1.precision()); }

// Mixed Float-ValidatedNumber operations
/** @brief Return the sum or unary positive value. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator+(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator+(x1, F(y2, x1.precision()));
    } else {
        return operator+(x1, Bounds<F>(y2, x1.precision()));
    }
}
/** @brief Return the difference or additive inverse. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator-(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator-(x1, F(y2, x1.precision()));
    } else {
        return operator-(x1, Bounds<F>(y2, x1.precision()));
    }
}
/** @brief Return the product. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator*(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator*(x1, F(y2, x1.precision()));
    } else {
        return operator*(x1, Bounds<F>(y2, x1.precision()));
    }
}
/** @brief Return the quotient. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator/(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator/(x1, F(y2, x1.precision()));
    } else {
        return operator/(x1, Bounds<F>(y2, x1.precision()));
    }
}
/** @brief Return the sum or unary positive value. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator+(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator+(F(y1, x2.precision()), x2);
    } else {
        return operator+(Bounds<F>(y1, x2.precision()), x2);
    }
}
/** @brief Return the difference or additive inverse. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator-(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator-(F(y1, x2.precision()), x2);
    } else {
        return operator-(Bounds<F>(y1, x2.precision()), x2);
    }
}
/** @brief Return the product. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator*(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator*(F(y1, x2.precision()), x2);
    } else {
        return operator*(Bounds<F>(y1, x2.precision()), x2);
    }
}
/** @brief Return the quotient. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator/(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator/(F(y1, x2.precision()), x2);
    } else {
        return operator/(Bounds<F>(y1, x2.precision()), x2);
    }
}

/** @brief Return the sum of the operands. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> add(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return add(x1, F(y2, x1.precision()));
    } else {
        return add(x1, Bounds<F>(y2, x1.precision()));
    }
}
/** @brief Return the difference of the operands. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> sub(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return sub(x1, F(y2, x1.precision()));
    } else {
        return sub(x1, Bounds<F>(y2, x1.precision()));
    }
}
/** @brief Return the product of the operands. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> mul(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return mul(x1, F(y2, x1.precision()));
    } else {
        return mul(x1, Bounds<F>(y2, x1.precision()));
    }
}
/** @brief Return the quotient of the operands. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> div(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return div(x1, F(y2, x1.precision()));
    } else {
        return div(x1, Bounds<F>(y2, x1.precision()));
    }
}
/** @brief Return the sum of the operands. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> add(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return add(F(y1, x2.precision()), x2);
    } else {
        return add(Bounds<F>(y1, x2.precision()), x2);
    }
}
/** @brief Return the difference of the operands. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> sub(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return sub(F(y1, x2.precision()), x2);
    } else {
        return sub(Bounds<F>(y1, x2.precision()), x2);
    }
}
/** @brief Return the product of the operands. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> mul(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return mul(F(y1, x2.precision()), x2);
    } else {
        return mul(Bounds<F>(y1, x2.precision()), x2);
    }
}
/** @brief Return the quotient of the operands. */
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> div(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return div(F(y1, x2.precision()), x2);
    } else {
        return div(Bounds<F>(y1, x2.precision()), x2);
    }
}

/** @brief Return the maximum of the operands. */
template<ARawFloat F, AValidatedNumber Y>
decltype(auto) max(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return max(x1, F(y2, x1.precision()));
    } else {
        return max(x1, Bounds<F>(y2, x1.precision()));
    }
}
/** @brief Return the minimum of the operands. */
template<ARawFloat F, AValidatedNumber Y>
decltype(auto) min(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return min(x1, F(y2, x1.precision()));
    } else {
        return min(x1, Bounds<F>(y2, x1.precision()));
    }
}
/** @brief Return the maximum of the operands. */
template<ARawFloat F, AValidatedNumber Y>
decltype(auto) max(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return max(F(y1, x2.precision()), x2);
    } else {
        return max(Bounds<F>(y1, x2.precision()), x2);
    }
}
/** @brief Return the minimum of the operands. */
template<ARawFloat F, AValidatedNumber Y>
decltype(auto) min(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return min(F(y1, x2.precision()), x2);
    } else {
        return min(Bounds<F>(y1, x2.precision()), x2);
    }
}

template<class X>
class Bounds;

/** @brief Declare rounded operations for the supported operand types. */
template<class FLT>
class DeclareRoundedOperations {
    /** @brief Return the sum or unary positive value. */
    friend Bounds<FloatDP> operator+(const FloatDP& x1, const FloatDP& x2);
    /** @brief Return the difference or additive inverse. */
    friend Bounds<FloatDP> operator-(const FloatDP& x1, const FloatDP& x2);
    /** @brief Return the product. */
    friend Bounds<FloatDP> operator*(const FloatDP& x1, const FloatDP& x2);
    /** @brief Return the quotient. */
    friend Bounds<FloatDP> operator/(const FloatDP& x1, const FloatDP& x2);

    /** @brief Return the multiplicative inverse. */
    friend Bounds<FloatDP> rec(const FloatDP& x);
    /** @brief Return the sum of the operands. */
    friend Bounds<FloatDP> add(const FloatDP& x1, const FloatDP& x2);
    /** @brief Return the difference of the operands. */
    friend Bounds<FloatDP> sub(const FloatDP& x1, const FloatDP& x2);
    /** @brief Return the product of the operands. */
    friend Bounds<FloatDP> mul(const FloatDP& x1, const FloatDP& x2);
    /** @brief Return the quotient of the operands. */
    friend Bounds<FloatDP> div(const FloatDP& x1, const FloatDP& x2);
    /** @brief Compute the product of the first two operands plus the third. */
    friend Bounds<FloatDP> fma(const FloatDP& x1, const FloatDP& x2, const FloatDP& x3);
    /** @brief Raise the value to the specified power. */
    friend Bounds<FloatDP> pow(const FloatDP& x, Nat m);
    /** @brief Raise the value to the specified power. */
    friend Bounds<FloatDP> pow(const FloatDP& x, Int n);
    /** @brief Return the square root. */
    friend Bounds<FloatDP> sqrt(const FloatDP& x);
    /** @brief Return the exponential. */
    friend Bounds<FloatDP> exp(const FloatDP& x);
    /** @brief Return the natural logarithm. */
    friend Bounds<FloatDP> log(const FloatDP& x);
    /** @brief Return the sine. */
    friend Bounds<FloatDP> sin(const FloatDP& x);
    /** @brief Return the cosine. */
    friend Bounds<FloatDP> cos(const FloatDP& x);
    /** @brief Return the tangent. */
    friend Bounds<FloatDP> tan(const FloatDP& x);
    /** @brief Return the inverse sine. */
    friend Bounds<FloatDP> asin(const FloatDP& x);
    /** @brief Return the inverse cosine. */
    friend Bounds<FloatDP> acos(const FloatDP& x);
    /** @brief Return the inverse tangent. */
    friend Bounds<FloatDP> atan(const FloatDP& x);
};

} // namespace Ariadne

#include "float_factory.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"

#endif
