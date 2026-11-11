/***************************************************************************
 *            numeric/positive.hpp
 *
 *  Copyright  2008-22  Pieter Collins
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

/** @file numeric/positive.hpp
 *  @brief Positive numbers.
 */

#ifndef ARIADNE_POSITIVE_HPP
#define ARIADNE_POSITIVE_HPP

#include "float.decl.hpp"
#include "foundation/logical.decl.hpp"
#include "number.decl.hpp"

namespace Ariadne {

/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<ARawFloat F>
class Positive<F> : public F {
    using PR = typename F::PrecisionType;

  public:
    /** @brief Construct a \c Positive with the specified characteristics. */
    Positive() : F() {}
    /** @brief Construct a \c Positive from pr. */
    explicit Positive(PR const& pr) : F(pr) {}
    /** @brief Construct a \c Positive from m, pr. */
    template<BuiltinUnsignedIntegral M>
    Positive(M m, PR pr) : F(m, pr) {}
    /** @brief Construct a \c Positive from ex, pr. */
    Positive(const TwoExp& ex, PR pr) : F(ex, pr) {}
    /** @brief Construct a \c Positive from w, pr. */
    explicit Positive(const Dyadic& w, PR pr) : F(w, pr) {}
    /** @brief Construct a \c Positive from x. */
    explicit Positive(const F& x) : F(x) {}

  public:
    /** @brief Return the sum or unary positive value. */
    friend Positive<Bounds<F>> operator+(const Positive<F>& x1, const Positive<F>& x2) {
        return cast_positive(cast_unsigned(x1) + cast_unsigned(x2));
    }
    /** @brief Return the product. */
    friend Positive<Bounds<F>> operator*(const Positive<F>& x1, const Positive<F>& x2) {
        return cast_positive(cast_unsigned(x1) * cast_unsigned(x2));
    }
    /** @brief Return the quotient. */
    friend Positive<Bounds<F>> operator/(const Positive<F>& x1, const Positive<F>& x2) {
        return cast_positive(cast_unsigned(x1) / cast_unsigned(x2));
    }
    /** @brief Return the maximum of the operands. */
    friend Positive<F> max(const Positive<F>& x1, const Positive<F>& x2) {
        return cast_positive(max(cast_unsigned(x1), cast_unsigned(x2)));
    }
    /** @brief Return the maximum of the operands. */
    friend Positive<F> max(const Positive<F>& x1, const F& x2) {
        return cast_positive(max(cast_unsigned(x1), x2));
    }
    /** @brief Return the maximum of the operands. */
    friend Positive<F> max(const F& x1, const Positive<F>& x2) {
        return cast_positive(max(x1, cast_unsigned(x2)));
    }
    /** @brief Return the minimum of the operands. */
    friend Positive<F> min(const Positive<F>& x1, const Positive<F>& x2) {
        return cast_positive(min(cast_unsigned(x1), cast_unsigned(x2)));
    }

    /** @brief Return the additive zero with the operand characteristics. */
    friend Positive<F> nul(const Positive<F>& x) { return Positive<F>(nul(static_cast<const F&>(x))); }
    /** @brief Return the unary positive value. */
    friend Positive<F> pos(const Positive<F>& x) { return Positive<F>(pos(static_cast<const F&>(x))); }
    /** @brief Divide the value by two. */
    friend Positive<F> hlf(const Positive<F>& x) { return Positive<F>(hlf(static_cast<const F&>(x))); }
    /** @brief Raise the value to the specified power. */
    friend Positive<Bounds<F>> pow(const Positive<F>& x, Nat m) {
        return pow(Positive<Bounds<F>>(x), m);
    }
    /** @brief Raise the value to the specified power. */
    friend Positive<Bounds<F>> pow(const Positive<F>& x, Int n) {
        return pow(Positive<Bounds<F>>(x), n);
    }

    /** @brief Return the sum or unary positive value. */
    friend Error<F> operator+(const Positive<F>& v1, const Error<F>& e2);
};

/** @brief Convert to the positive representation. */
template<class F>
inline Positive<F> cast_positive(const F& x) {
    return Positive<F>(x);
}
/** @brief Convert to the positive representation. */
template<class F>
inline Positive<F> cast_positive(const Positive<F>& x) {
    return x;
}
/** @brief Convert to the unsigned representation. */
template<class F>
inline const F& cast_unsigned(const Positive<F>& x) {
    return static_cast<const F&>(x);
}

} // namespace Ariadne

#endif
