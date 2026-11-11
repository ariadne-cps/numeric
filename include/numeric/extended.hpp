/***************************************************************************
 *            numeric/extended.hpp
 *
 *  Copyright  2017-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

/** @file numeric/extended.hpp
 *  @brief Extended real numbers.
 */

#ifndef ARIADNE_EXTENDED_HPP
#define ARIADNE_EXTENDED_HPP

#include "numeric/dyadic.hpp"
#include "numeric/integer.hpp"
#include "numeric/rational.hpp"
#include "numeric/sign.hpp"

namespace Ariadne {

/** @brief Compare the operands using the comparison result type. */
inline Comparison cmp(int n1, int n2) {
    return n1 == n2 ? Comparison::EQUAL : n1 < n2 ? Comparison::LESS
                                                  : Comparison::GREATER;
}

/** @brief Operations supported by finite objects. */
template<class X>
class FiniteOperations;
/** @brief Operations supported by extension objects. */
template<class X>
class ExtensionOperations;

//! \ingroup NumericModule
//! \brief A number of type \a X, extended with infinite and not-a-number values.
template<class X>
class ExtendedOperations {
    typedef FiniteOperations<X> Finite;
    typedef ExtensionOperations<X> Extension;

  public:
    /** @brief Test whether the value is not a number. */
    static Bool is_nan(const X& x) { return Extension::is_nan(x); }
    /** @brief Test whether the value is infinite. */
    static Bool is_inf(const X& x) { return Extension::is_inf(x); }
    /** @brief Test whether the value is finite. */
    static Bool is_finite(const X& x) { return Extension::is_finite(x); }
    /** @brief Test whether the value is zero. */
    static Bool is_zero(const X& x) { return Extension::is_zero(x); }

    /** @brief Set the nan. */
    static Void set_nan(X& r) { Extension::set_nan(r); }
    /** @brief Set the inf. */
    static Void set_inf(X& r, Sign sgn) { Extension::set_inf(r, sgn); }
    /** @brief Set the zero. */
    static Void set_zero(X& r) { Extension::set_zero(r); }

    /** @brief Return the sign of the value. */
    static Sign sgn(const X& x) { return Extension::sgn(x); }

    /** @brief Return the absolute value. */
    static Sign abs(Sign s) { return (s == Sign::NEGATIVE) ? -s : s; }
    /** @brief Raise the value to the specified power. */
    static Sign pow(Sign s, Nat m) { return (s == Sign::NEGATIVE && (m % 2 == 0)) ? Sign::POSITIVE : s; }

    /** @brief Return the additive inverse. */
    static Void neg(X& r, const X& x) {
        if (is_finite(x)) {
            Finite::neg(r, x);
        } else {
            set_inf(r, -sgn(x));
        }
    }

    /** @brief Return the sum of the operands. */
    static Void add(X& r, const X& x1, const X& x2) {
        if (is_finite(x1) and is_finite(x2)) {
            Finite::add(r, x1, x2);
        } else if (is_finite(x1)) {
            set_inf(r, sgn(x2));
        } else if (is_finite(x2)) {
            set_inf(r, sgn(x1));
        } else if (sgn(x1) == sgn(x2)) {
            set_inf(r, sgn(x1));
        } else {
            set_nan(r);
        }
    }

    /** @brief Return the difference of the operands. */
    static Void sub(X& r, const X& x1, const X& x2) {
        if (is_finite(x1) and is_finite(x2)) {
            Finite::sub(r, x1, x2);
        } else if (is_finite(x1)) {
            set_inf(r, -sgn(x2));
        } else if (is_finite(x2)) {
            set_inf(r, sgn(x1));
        } else if (sgn(x1) == (-sgn(x2))) {
            set_inf(r, sgn(x1));
        } else {
            set_nan(r);
        }
    }
    /** @brief Return the product of the operands. */
    static Void mul(X& r, const X& x1, const X& x2) {
        if (is_finite(x1) and is_finite(x2)) {
            Finite::mul(r, x1, x2);
        } else {
            set_inf(r, sgn(x1) * sgn(x2));
        }
    }
    /** @brief Return the quotient of the operands. */
    static Void div(X& r, const X& x1, const X& x2) {
        if (is_finite(x1) and is_finite(x2)) {
            Finite::div(r, x1, x2);
        } else {
            rec(r, x2);
            mul(r, x1, r);
        }
    }

    /** @brief Return the square. */
    static Void sqr(X& r, const X& x) {
        if (is_finite(x)) {
            Finite::sqr(r, x);
        } else {
            set_inf(r, abs(sgn(x)));
        }
    }
    /** @brief Raise the value to the specified power. */
    static Void pow(X& r, const X& x, Nat m) {
        if (is_finite(x)) {
            Finite::pow(r, x, m);
        } else if (is_nan(x)) {
            set_nan(r);
        } else {
            set_inf(r, (sgn(x) == Sign::POSITIVE || (m % 2 == 0)) ? Sign::POSITIVE : Sign::NEGATIVE);
        }
    }
    /** @brief Raise the value to the specified power. */
    static Void pow(X& r, const X& x, Int n) {
        if (is_finite(x)) {
            Finite::pow(r, x, n);
        } else {
            assert(false);
        }
    }

    /** @brief Divide the value by two. */
    static Void hlf(X& r, const X& x) {
        if (is_finite(x)) {
            Finite::hlf(r, x);
        } else {
            set_inf(r, sgn(x));
        }
    }
    /** @brief Return the multiplicative inverse. */
    static Void rec(X& r, const X& x) {
        if (is_finite(x)) {
            Finite::rec(r, x);
        } else if (is_inf(x)) {
            set_zero(r);
        } else {
            set_nan(r);
        }
    }

    /** @brief Return the absolute value. */
    static Void abs(X& r, const X& x) {
        if (is_finite(x)) {
            Finite::abs(r, x);
        } else {
            set_inf(r, abs(sgn(x)));
        }
    }
    /** @brief Return the maximum of the operands. */
    static Void max(X& r, const X& x1, const X& x2) {
        if (is_finite(x1) and is_finite(x2)) {
            Finite::max(r, x1, x2);
        } else if (is_nan(x1) or is_nan(x2)) {
            set_nan(r);
        } else if ((is_inf(x1) and sgn(x1) == Sign::POSITIVE) or (is_inf(x2) and sgn(x2) == Sign::POSITIVE)) {
            set_inf(r, Sign::POSITIVE);
        } else if (is_inf(x1)) {
            Finite::set(r, x2);
        } else {
            Finite::set(r, x1);
        }
    } // As x1==-inf or x2==-inf
    /** @brief Return the minimum of the operands. */
    static Void min(X& r, const X& x1, const X& x2) {
        if (is_finite(x1) and is_finite(x2)) {
            Finite::min(r, x1, x2);
        } else if (is_nan(x1) or is_nan(x2)) {
            set_nan(r);
        } else if ((is_inf(x1) and sgn(x1) == Sign::NEGATIVE) or (is_inf(x2) and sgn(x2) == Sign::NEGATIVE)) {
            set_inf(r, Sign::NEGATIVE);
        } else if (is_inf(x1)) {
            Finite::set(r, x2);
        } else {
            Finite::set(r, x1);
        }
    } // As x1==+inf or x2==+inf

    /** @brief Test equality of the operands. */
    static Boolean eq(const X& x1, const X& x2) {
        if (is_finite(x1) and is_finite(x2)) {
            return Finite::eq(x1, x2);
        } else if (is_finite(x1) or is_finite(x2)) {
            return false;
        } else if (is_nan(x1) or is_nan(x2)) {
            return false;
        } else {
            return sgn(x1) == sgn(x2);
        }
    }

    /** @brief Compare the operands using the comparison result type. */
    static Comparison cmp(const X& x1, const X& x2) {
        if (is_finite(x1) and is_finite(x2)) {
            return Finite::cmp(x1, x2);
        } else if (is_nan(x1) or is_nan(x2)) {
            return Comparison::INCOMPARABLE;
        } else if (is_finite(x1)) {
            return 0 < (int)sgn(x2) ? Comparison::LESS : Comparison::GREATER;
        } else if (is_finite(x2)) {
            return (int)sgn(x1) > 0 ? Comparison::GREATER : Comparison::LESS;
        } else {
            return Ariadne::cmp((int)sgn(x1), (int)sgn(x2));
        }
    }

    /** @brief Write the value to the output stream. */
    static OutputStream& write(OutputStream& os, const X& x) {
        if (is_finite(x)) {
            return os << x;
        } else if (is_nan(x)) {
            return os << "nan";
        } else {
            if (sgn(x) == Sign::NEGATIVE) {
                os << '-';
            }
            return os << "inf";
        }
    }
};

} // namespace Ariadne

#endif
