/***************************************************************************
 *            numeric/float_operations.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

template<class X, class R = X>
struct DeclareFloatOperations
    : DeclareRealOperations<X, R>,
      DeclareComparisonOperations<X, LessTrait<X>, EqualsTrait<X>> {
    friend OutputStream& operator<<(OutputStream&, const X&);
    friend InputStream& operator>>(InputStream&, X&);
};

template<class X, class NX = OppositeTrait<X>>
struct DeclareDirectedFloatOperations
    : DeclareDirectedNumericOperations<X, NX>,
      DeclareMixedDirectedGroupOperators<X, NX, GenericType<X>, GenericType<NX>> {
    friend OutputStream& operator<<(OutputStream&, const X&);
    friend InputStream& operator>>(InputStream&, X&);
};

template<class PX, class QPX = OppositeTrait<PX>>
struct DeclarePositiveDirectedFloatOperations
    : DeclarePositiveDirectedNumericOperations<PX, QPX> {
};

template<class PX>
struct DeclarePositiveFloatOperations
    : DeclarePositiveDirectedFloatOperations<PX, PX> {
};

template<class X, class R = X>
struct DefineFloatOperations
    : public DefineFieldOperators<X, R>,
      public DefineComparisonOperators<X, LessTrait<X>, EqualsTrait<X>>,
      public DefineConcreteGenericOperators<X>,
      public DeclareFloatOperations<X> {
};

template<class X, class NX>
struct DefineDirectedFloatOperations
    : public DefineDirectedGroupOperators<X, NX>,
      public DefineDirectedGroupOperators<NX, X>,
      public DefineDirectedComparisonOperators<X, NX, LessTrait<X>, EqualsTrait<X>>,
      public DefineDirectedComparisonOperators<NX, X, LessTrait<NX>, EqualsTrait<NX>>,
      public DefineConcreteGenericOperators<X> {
};

template<class X, class R = X>
struct DispatchFloatOperations
    : DispatchNumericOperations<X, R>,
      DispatchComparisonOperations<X, LessTrait<X>, EqualsTrait<X>>,
      DefineConcreteGenericOperators<X> {
    friend OutputStream& operator<<(OutputStream& os, const X& x) { return Operations<X>::_write(os, x); }
    friend InputStream& operator>>(InputStream& is, X& x) { return Operations<X>::_read(is, x); }
};

template<class X>
struct DispatchDirectedFloatOperations
    : DispatchDirectedNumericOperations<X, OppositeTrait<X>>,
      DispatchDirectedNumericOperations<OppositeTrait<X>, X>,
      DispatchDirectedComparisonOperations<X, OppositeTrait<X>, LessTrait<X>, EqualsTrait<X>>,
      DispatchDirectedComparisonOperations<OppositeTrait<X>, X, LessTrait<OppositeTrait<X>>, EqualsTrait<OppositeTrait<X>>>,
      DefineConcreteGenericOperators<X> {
    friend OutputStream& operator<<(OutputStream& os, const X& x) { return Operations<X>::_write(os, x); }
    friend InputStream& operator>>(InputStream& is, X& x) { return Operations<X>::_read(is, x); }
};

template<class PX, class QPX = OppositeTrait<PX>>
struct DispatchPositiveDirectedFloatOperations
    : DispatchPositiveDirectedNumericOperations<PX, QPX>,
      DispatchPositiveDirectedNumericOperations<QPX, PX>,
      DefineConcreteGenericOperators<PX> {
};

template<class PX>
struct DispatchPositiveFloatOperations
    : DispatchPositiveDirectedNumericOperations<PX, PX>,
      DefineConcreteGenericOperators<PX> {
};

template<class Y1, class Y2>
concept HaveCmp = requires(Y1 y1, Y2 y2) {
    { cmp(y1, y2) } -> SameAs<Comparison>;
};
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

template<class Y1, class Y2>
    requires HaveCmp<Y1, Y2>
Boolean operator==(Y1 const& y1, Y2 const& y2) { return _cmp(Eq(), y1, y2); }
template<class Y1, class Y2>
    requires HaveCmp<Y1, Y2>
Boolean operator!=(Y1 const& y1, Y2 const& y2) { return _cmp(Ne(), y1, y2); }
template<class Y1, class Y2>
    requires HaveCmp<Y1, Y2>
Boolean operator<(Y1 const& y1, Y2 const& y2) { return _cmp(Lt(), y1, y2); }
template<class Y1, class Y2>
    requires HaveCmp<Y1, Y2>
Boolean operator>(Y1 const& y1, Y2 const& y2) { return _cmp(Gt(), y1, y2); }
template<class Y1, class Y2>
    requires HaveCmp<Y1, Y2>
Boolean operator<=(Y1 const& y1, Y2 const& y2) { return _cmp(Le(), y1, y2); }
template<class Y1, class Y2>
    requires HaveCmp<Y1, Y2>
Boolean operator>=(Y1 const& y1, Y2 const& y2) { return _cmp(Ge(), y1, y2); }

// Mixed Float-ValidatedNumber comparisons
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator==(const Y& y1, const F& x2) { return Bounds<F>(y1, x2.precision()) == x2; }
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator!=(const Y& y1, const F& x2) { return Bounds<F>(y1, x2.precision()) != x2; }
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator<(const Y& y1, const F& x2) { return Bounds<F>(y1, x2.precision()) < x2; }
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator>(const Y& y1, const F& x2) { return Bounds<F>(y1, x2.precision()) > x2; }
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator<=(const Y& y1, const F& x2) { return Bounds<F>(y1, x2.precision()) <= x2; }
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator>=(const Y& y1, const F& x2) { return Bounds<F>(y1, x2.precision()) >= x2; }
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator==(const F& x1, const Y& y2) { return x1 == Bounds<F>(y2, x1.precision()); }
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator!=(const F& x1, const Y& y2) { return x1 != Bounds<F>(y2, x1.precision()); }
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator<(const F& x1, const Y& y2) { return x1 < Bounds<F>(y2, x1.precision()); }
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator>(const F& x1, const Y& y2) { return x1 > Bounds<F>(y2, x1.precision()); }
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator<=(const F& x1, const Y& y2) { return x1 <= Bounds<F>(y2, x1.precision()); }
template<ARawFloat F, ANonExactValidatedNumber Y>
ValidatedKleenean operator>=(const F& x1, const Y& y2) { return x1 >= Bounds<F>(y2, x1.precision()); }

// Mixed Float-ValidatedNumber operations
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator+(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator+(x1, F(y2, x1.precision()));
    } else {
        return operator+(x1, Bounds<F>(y2, x1.precision()));
    }
}
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator-(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator-(x1, F(y2, x1.precision()));
    } else {
        return operator-(x1, Bounds<F>(y2, x1.precision()));
    }
}
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator*(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator*(x1, F(y2, x1.precision()));
    } else {
        return operator*(x1, Bounds<F>(y2, x1.precision()));
    }
}
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator/(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator/(x1, F(y2, x1.precision()));
    } else {
        return operator/(x1, Bounds<F>(y2, x1.precision()));
    }
}
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator+(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator+(F(y1, x2.precision()), x2);
    } else {
        return operator+(Bounds<F>(y1, x2.precision()), x2);
    }
}
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator-(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator-(F(y1, x2.precision()), x2);
    } else {
        return operator-(Bounds<F>(y1, x2.precision()), x2);
    }
}
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator*(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator*(F(y1, x2.precision()), x2);
    } else {
        return operator*(Bounds<F>(y1, x2.precision()), x2);
    }
}
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> operator/(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return operator/(F(y1, x2.precision()), x2);
    } else {
        return operator/(Bounds<F>(y1, x2.precision()), x2);
    }
}

template<ARawFloat F, AValidatedNumber Y>
Bounds<F> add(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return add(x1, F(y2, x1.precision()));
    } else {
        return add(x1, Bounds<F>(y2, x1.precision()));
    }
}
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> sub(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return sub(x1, F(y2, x1.precision()));
    } else {
        return sub(x1, Bounds<F>(y2, x1.precision()));
    }
}
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> mul(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return mul(x1, F(y2, x1.precision()));
    } else {
        return mul(x1, Bounds<F>(y2, x1.precision()));
    }
}
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> div(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return div(x1, F(y2, x1.precision()));
    } else {
        return div(x1, Bounds<F>(y2, x1.precision()));
    }
}
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> add(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return add(F(y1, x2.precision()), x2);
    } else {
        return add(Bounds<F>(y1, x2.precision()), x2);
    }
}
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> sub(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return sub(F(y1, x2.precision()), x2);
    } else {
        return sub(Bounds<F>(y1, x2.precision()), x2);
    }
}
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> mul(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return mul(F(y1, x2.precision()), x2);
    } else {
        return mul(Bounds<F>(y1, x2.precision()), x2);
    }
}
template<ARawFloat F, AValidatedNumber Y>
Bounds<F> div(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return div(F(y1, x2.precision()), x2);
    } else {
        return div(Bounds<F>(y1, x2.precision()), x2);
    }
}

template<ARawFloat F, AValidatedNumber Y>
decltype(auto) max(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return max(x1, F(y2, x1.precision()));
    } else {
        return max(x1, Bounds<F>(y2, x1.precision()));
    }
}
template<ARawFloat F, AValidatedNumber Y>
decltype(auto) min(const F& x1, const Y& y2) {
    if constexpr (AnExactDyadic<Y>) {
        return min(x1, F(y2, x1.precision()));
    } else {
        return min(x1, Bounds<F>(y2, x1.precision()));
    }
}
template<ARawFloat F, AValidatedNumber Y>
decltype(auto) max(const Y& y1, const F& x2) {
    if constexpr (AnExactDyadic<Y>) {
        return max(F(y1, x2.precision()), x2);
    } else {
        return max(Bounds<F>(y1, x2.precision()), x2);
    }
}
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

template<class FLT>
class DeclareRoundedOperations {
    friend Bounds<FloatDP> operator+(const FloatDP& x1, const FloatDP& x2);
    friend Bounds<FloatDP> operator-(const FloatDP& x1, const FloatDP& x2);
    friend Bounds<FloatDP> operator*(const FloatDP& x1, const FloatDP& x2);
    friend Bounds<FloatDP> operator/(const FloatDP& x1, const FloatDP& x2);

    friend Bounds<FloatDP> rec(const FloatDP& x);
    friend Bounds<FloatDP> add(const FloatDP& x1, const FloatDP& x2);
    friend Bounds<FloatDP> sub(const FloatDP& x1, const FloatDP& x2);
    friend Bounds<FloatDP> mul(const FloatDP& x1, const FloatDP& x2);
    friend Bounds<FloatDP> div(const FloatDP& x1, const FloatDP& x2);
    friend Bounds<FloatDP> fma(const FloatDP& x1, const FloatDP& x2, const FloatDP& x3);
    friend Bounds<FloatDP> pow(const FloatDP& x, Nat m);
    friend Bounds<FloatDP> pow(const FloatDP& x, Int n);
    friend Bounds<FloatDP> sqrt(const FloatDP& x);
    friend Bounds<FloatDP> exp(const FloatDP& x);
    friend Bounds<FloatDP> log(const FloatDP& x);
    friend Bounds<FloatDP> sin(const FloatDP& x);
    friend Bounds<FloatDP> cos(const FloatDP& x);
    friend Bounds<FloatDP> tan(const FloatDP& x);
    friend Bounds<FloatDP> asin(const FloatDP& x);
    friend Bounds<FloatDP> acos(const FloatDP& x);
    friend Bounds<FloatDP> atan(const FloatDP& x);
};

} // namespace Ariadne

#include "float_factory.hpp"
#include "floatdp.hpp"
#include "floatmp.hpp"

#endif
