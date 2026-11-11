/***************************************************************************
 *            numeric/builtin.hpp
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

/** @file numeric/builtin.hpp
 *  @brief Inclusion header for wrapper for builtin numbers.
 */

#ifndef ARIADNE_BUILTIN_HPP
#define ARIADNE_BUILTIN_HPP

#include <cassert>
#include <cfenv>
#include <cfloat>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

#include "foundation/logical.hpp"
#include "numeric/number.decl.hpp"
#include "twoexp.hpp"
#include "utility/macros.hpp"
#include "utility/metaprogramming.hpp"

namespace Ariadne {

struct ExactTag;
enum class Comparison : ComparableEnumerationType;

/** @brief Compile-time predicate for conversion to a built-in double. */
template<class D>
struct HasGetD {
    /** @brief Evaluate the compile-time or runtime predicate for the supplied arguments. */
    template<class DD, class = decltype(std::declval<DD>().get_d())>
    static std::true_type test(int);
    /** @brief Evaluate the compile-time or runtime predicate for the supplied arguments. */
    template<class DD>
    static std::false_type test(...);
    static const bool value = decltype(test<D>(1))::value;
};

//! \ingroup NumericModule
//! \brief A wrapper around a builtin double-precision floating-point number,
//! indicating that the stored value is an \em approximation to the value of a real quantity.
class ApproximateDouble {
    double _d;

  public:
    typedef ApproximateTag Paradigm;
    typedef ApproximateDouble NumericType;
    /** @brief Construct a \c ApproximateDouble with the specified characteristics. */
    ApproximateDouble() : _d() {}
    /** @brief Construct a \c ApproximateDouble from d. */
    ApproximateDouble(double d) : _d(d) {}
    /** @brief Construct a \c ApproximateDouble from x. */
    ApproximateDouble(const ExactDouble& x);
    /** @brief Construct a \c ApproximateDouble from r0. */
    ApproximateDouble(const Real& r0);
    /** @brief Convert to double. */
    explicit operator double() const { return _d; }
    /** @brief Convert to a built-in double. */
    double get_d() const { return this->_d; }
    friend ApproximateDouble operator""_a(long double lx) {
        double x = static_cast<double>(lx);
        return ApproximateDouble(x);
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, ApproximateDouble x) { return os << x._d; }

    /** @brief Return the additive zero with the operand characteristics. */
    friend ApproximateDouble nul(ApproximateDouble) { return ApproximateDouble(0.0); }
    /** @brief Return the unary positive value. */
    friend ApproximateDouble pos(ApproximateDouble x) { return ApproximateDouble(+x._d); }
    /** @brief Return the additive inverse. */
    friend ApproximateDouble neg(ApproximateDouble x) { return ApproximateDouble(-x._d); }
    /** @brief Return the sum of the operands. */
    friend ApproximateDouble add(ApproximateDouble x1, ApproximateDouble x2) { return ApproximateDouble(x1._d + x2._d); }
    /** @brief Return the difference of the operands. */
    friend ApproximateDouble sub(ApproximateDouble x1, ApproximateDouble x2) { return ApproximateDouble(x1._d - x2._d); }
    /** @brief Return the product of the operands. */
    friend ApproximateDouble mul(ApproximateDouble x1, ApproximateDouble x2) { return ApproximateDouble(x1._d * x2._d); }
    /** @brief Return the quotient of the operands. */
    friend ApproximateDouble div(ApproximateDouble x1, ApproximateDouble x2) { return ApproximateDouble(x1._d / x2._d); }
    /** @brief Return the absolute value. */
    friend ApproximateDouble abs(ApproximateDouble x) { return ApproximateDouble(x._d >= 0 ? x._d : -x._d); }
    /** @brief Return the maximum of the operands. */
    friend ApproximateDouble max(ApproximateDouble x1, ApproximateDouble x2) { return ApproximateDouble(x1._d >= x2._d ? x1._d : x2._d); }
    /** @brief Return the minimum of the operands. */
    friend ApproximateDouble min(ApproximateDouble x1, ApproximateDouble x2) { return ApproximateDouble(x1._d <= x2._d ? x1._d : x2._d); }
    /** @brief Return an upper bound on the absolute value. */
    friend ApproximateDouble mag(ApproximateDouble x) { return abs(x); }
    /** @brief Return a lower bound on the absolute value. */
    friend ApproximateDouble mig(ApproximateDouble x) { return abs(x); }

    /** @brief Return the sum or unary positive value. */
    friend ApproximateDouble operator+(ApproximateDouble x) { return ApproximateDouble(+x._d); }
    /** @brief Return the difference or additive inverse. */
    friend ApproximateDouble operator-(ApproximateDouble x) { return ApproximateDouble(-x._d); }
    /** @brief Return the sum or unary positive value. */
    friend ApproximateDouble operator+(ApproximateDouble x1, ApproximateDouble x2) { return ApproximateDouble(x1._d + x2._d); }
    /** @brief Return the difference or additive inverse. */
    friend ApproximateDouble operator-(ApproximateDouble x1, ApproximateDouble x2) { return ApproximateDouble(x1._d - x2._d); }
    /** @brief Return the product. */
    friend ApproximateDouble operator*(ApproximateDouble x1, ApproximateDouble x2) { return ApproximateDouble(x1._d * x2._d); }
    /** @brief Return the quotient. */
    friend ApproximateDouble operator/(ApproximateDouble x1, ApproximateDouble x2) { return ApproximateDouble(x1._d / x2._d); }
    /** @brief Return the product. */
    friend ApproximateDouble operator*(ApproximateDouble x, TwoExp p) { return ApproximateDouble(x.get_d() * p.get_d()); }
    /** @brief Return the quotient. */
    friend ApproximateDouble operator/(ApproximateDouble x, TwoExp p) { return ApproximateDouble(x.get_d() / p.get_d()); }
    /** @brief Add the operand in place. */
    friend ApproximateDouble& operator+=(ApproximateDouble& x1, ApproximateDouble x2) {
        x1._d += x2._d;
        return x1;
    }
    /** @brief Subtract the operand in place. */
    friend ApproximateDouble& operator-=(ApproximateDouble& x1, ApproximateDouble x2) {
        x1._d -= x2._d;
        return x1;
    }
    /** @brief Multiply by the operand in place. */
    friend ApproximateDouble& operator*=(ApproximateDouble& x1, ApproximateDouble x2) {
        x1._d *= x2._d;
        return x1;
    }
    /** @brief Divide by the operand in place. */
    friend ApproximateDouble& operator/=(ApproximateDouble& x1, ApproximateDouble x2) {
        x1._d /= x2._d;
        return x1;
    }

    /** @brief Test whether the operands have the same representation. */
    friend Bool same(ApproximateDouble x1, ApproximateDouble x2) { return x1._d == x2._d; }
    /** @brief Test equality of the operands. */
    friend ApproximateKleenean operator==(ApproximateDouble x1, ApproximateDouble x2) { return x1._d == x2._d; }
    /** @brief Test inequality of the operands. */
    friend ApproximateKleenean operator!=(ApproximateDouble x1, ApproximateDouble x2) { return x1._d != x2._d; }
    /** @brief Test whether the first operand is greater than the second. */
    friend ApproximateKleenean operator>(ApproximateDouble x1, ApproximateDouble x2) { return x1._d > x2._d; }
    /** @brief Test whether the first operand is less than the second. */
    friend ApproximateKleenean operator<(ApproximateDouble x1, ApproximateDouble x2) { return x1._d < x2._d; }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ApproximateKleenean operator>=(ApproximateDouble x1, ApproximateDouble x2) { return x1._d >= x2._d; }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ApproximateKleenean operator<=(ApproximateDouble x1, ApproximateDouble x2) { return x1._d <= x2._d; }
};

//! \ingroup NumericModule
//! \brief Indicate that a floating-point literal only represents an double-precision approximation to a number.
inline ApproximateDouble operator""_a(long double lx);

//! \ingroup NumericModule
//! \brief A wrapper around a builtin double-precision floating-point number,
//! indicating that the stored value is the \em exact value of a real quantity.
class ExactDouble {
    double _d;
    /** @brief Convert to a signed integer after checking the representable range. */
    static double _checked_signed(long long n) {
        double d = static_cast<double>(n);
        ARIADNE_PRECONDITION(d < std::ldexp(1.0, 63));
        ARIADNE_PRECONDITION(static_cast<long long>(d) == n);
        return d;
    }
    /** @brief Convert to an unsigned integer after checking the representable range. */
    static double _checked_unsigned(unsigned long long n) {
        double d = static_cast<double>(n);
        ARIADNE_PRECONDITION(d < std::ldexp(1.0, 64));
        ARIADNE_PRECONDITION(static_cast<unsigned long long>(d) == n);
        return d;
    }
#if LDBL_MANT_DIG > DBL_MANT_DIG
    /** @brief Convert a long double after checking exact representability. */
    static double _checked_long_double(long double x) {
        double d = static_cast<double>(x);
        ARIADNE_PRECONDITION(std::isnan(d) || static_cast<long double>(d) == x);
        return d;
    }
#endif

  public:
    typedef ExactTag Paradigm;
    /** @brief Convert to a built-in double. */
    double get_d() const { return this->_d; }
    /** @brief Construct a \c ExactDouble with the specified characteristics. */
    ExactDouble() : _d() {}
    /** @brief Construct a \c ExactDouble from n. */
    template<BuiltinIntegral N>
        requires std::is_signed_v<N>
    ExactDouble(N n) : _d(_checked_signed(static_cast<long long>(n))) {}
    /** @brief Construct a \c ExactDouble from n. */
    template<BuiltinIntegral N>
        requires(!std::is_signed_v<N>)
    ExactDouble(N n) : _d(_checked_unsigned(static_cast<unsigned long long>(n))) {}
    /** @brief Construct a \c ExactDouble from x. */
    explicit ExactDouble(float x) : _d(static_cast<double>(x)) {}
    /** @brief Construct a \c ExactDouble from x. */
    explicit ExactDouble(double x) : _d(x) {}
#if LDBL_MANT_DIG > DBL_MANT_DIG
    /** @brief Construct a \c ExactDouble from x. */
    explicit ExactDouble(long double x) : _d(_checked_long_double(x)) {}
#else
    /** @brief Construct a \c ExactDouble from x. */
    explicit ExactDouble(long double x) : _d(static_cast<double>(x)) {}
#endif
    /** @brief Return an infinite value. */
    static ExactDouble infinity() { return ExactDouble(std::numeric_limits<double>::infinity()); }
    /** @brief Convert to exact number. */
    operator ExactNumber() const;
    /** @brief Return the additive zero with the operand characteristics. */
    friend ExactDouble nul(ExactDouble) { return ExactDouble(0.0); }
    /** @brief Return the absolute value. */
    friend ExactDouble abs(ExactDouble x) { return ExactDouble(std::abs(x._d)); }
    /** @brief Return an infinite value. */
    static ExactDouble inf() { return ExactDouble(std::numeric_limits<double>::infinity()); }
    /** @brief Return a not-a-number value. */
    static ExactDouble nan() { return ExactDouble(std::numeric_limits<double>::quiet_NaN()); }
    /** @brief Test whether the value is finite. */
    friend Bool is_finite(ExactDouble x) { return std::isfinite(x._d); }
    /** @brief Return the sum or unary positive value. */
    friend ExactDouble operator+(ExactDouble x) { return ExactDouble(+x._d); }
    /** @brief Return the difference or additive inverse. */
    friend ExactDouble operator-(ExactDouble x) { return ExactDouble(-x._d); }
    /** @brief Return the sum or unary positive value. */
    friend ApproximateDouble operator+(ApproximateDouble x1, ApproximateDouble x2);
    /** @brief Return the difference or additive inverse. */
    friend ApproximateDouble operator-(ApproximateDouble x1, ApproximateDouble x2);
    /** @brief Return the product. */
    friend ApproximateDouble operator*(ApproximateDouble x1, ApproximateDouble x2);
    /** @brief Return the quotient. */
    friend ApproximateDouble operator/(ApproximateDouble x1, ApproximateDouble x2);
    /** @brief Return the product. */
    friend ExactDouble operator*(ExactDouble x, TwoExp p) { return ExactDouble(x.get_d() * p.get_d()); }
    /** @brief Return the quotient. */
    friend ExactDouble operator/(ExactDouble x, TwoExp p) { return ExactDouble(x.get_d() / p.get_d()); }
    /** @brief Compare the operands using the comparison result type. */
    friend Comparison cmp(const ExactDouble& x1, const ExactDouble& x2) {
        return Comparison((x1._d == x2._d) ? 0 : (x1._d < x2._d) ? -1
                                                                 : +1);
    }
    /** @brief Test equality of the operands. */
    friend Boolean operator==(const ExactDouble& x1, const ExactDouble& x2) { return x1._d == x2._d; }
    /** @brief Test inequality of the operands. */
    friend Boolean operator!=(const ExactDouble& x1, const ExactDouble& x2) { return x1._d != x2._d; }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend Boolean operator>=(const ExactDouble& x1, const ExactDouble& x2) { return x1._d >= x2._d; }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend Boolean operator<=(const ExactDouble& x1, const ExactDouble& x2) { return x1._d <= x2._d; }
    /** @brief Test whether the first operand is greater than the second. */
    friend Boolean operator>(const ExactDouble& x1, const ExactDouble& x2) { return x1._d > x2._d; }
    /** @brief Test whether the first operand is less than the second. */
    friend Boolean operator<(const ExactDouble& x1, const ExactDouble& x2) { return x1._d < x2._d; }
    friend ExactDouble operator""_x(long double lx) { return ExactDouble(lx); }
    friend ExactDouble operator""_pr(long double lx) {
#if LDBL_MANT_DIG > DBL_MANT_DIG
        // A directed rounding mode must not change a precision-rounded literal.
        const int saved_rounding = std::fegetround();
        ARIADNE_PRECONDITION(saved_rounding != -1);
        if (saved_rounding != FE_TONEAREST) {
            ARIADNE_PRECONDITION(std::fesetround(FE_TONEAREST) == 0);
        }
        // Force the conversion under the selected mode, not at compile time.
        volatile long double input = lx;
        const double x = static_cast<double>(input);
        if (saved_rounding != FE_TONEAREST) {
            ARIADNE_PRECONDITION(std::fesetround(saved_rounding) == 0);
        }
#else
        const double x = static_cast<double>(lx);
#endif
        return ExactDouble(x);
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, ExactDouble x) { return os << std::setprecision(18) << x.get_d(); }
};

//! \ingroup ExtendedLiteralsSubModule
//! \brief Indicate that a floating-point literal is the exact value of a number in double-precision.
//! \details The input \a lx is converted first converted to a double-precision value \a x.
//! If \a lx and \a x differ, then the decimal literal input to \a lx almost certainly does not represent a double-precision number exactly, and an assertion fails. <p/>
//! For example, <c>0.625_x</c> yields the exact double-precision value \f$5/2^3\f$, but <c>0.6_x</c> fails, since \f$3/5\f$ is not exactly representable as a double-precision number.
ExactDouble operator""_x(long double lx);

//! \ingroup ExtendedLiteralsSubModule
//! \brief Indicate that a floating-point literal is sufficiently precise that it can safely be taken to represent the closest double-precision floating-point number.
inline ExactDouble exact(long double lx) { return operator""_x(lx); }

//! \ingroup ExtendedLiteralsSubModule
//! \brief Indicate that a floating-point literal is sufficiently precise that it can safely be taken to represent the closest double-precision floating-point number.
//! \details For example, writing <c>0.14285714285714285_pr</c> claims that the exact decimal "0.14285714285714285" should be interpreted as the closest double-precision floating-point number, which has exact value \f$0.142857142857142849212692681248881854116916656494140625\f$.
inline ExactDouble operator""_pr(long double lx);

#ifdef DOXYGRN
//! \ingroup ExtendedLiteralsSubModule
//! \brief Indicate that a floating-point literal is the exact value of a number in double-precision. (Python-only syntax)
ExactDouble x_(long double lx);

//! \ingroup ExtendedLiteralsSubModule
//! \brief Indicate that a floating-point literal is sufficiently precise that it can safely be taken to represent the closest double-precision floating-point number. (Python-only syntax)
ExactDouble pr_(long double lx);
#endif // DOXYGEN

inline ExactDouble cast_exact(double d) { return ExactDouble(d); }
/** @brief Convert to the exact representation. */
inline ExactDouble cast_exact(ApproximateDouble ax) { return ExactDouble(ax.get_d()); }

static const ExactDouble inf = ExactDouble(std::numeric_limits<double>::infinity());

template<class X>
class Positive;
/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<>
class Positive<ExactDouble> : public ExactDouble {
  public:
    /** @brief Construct a \c Positive from x. */
    explicit Positive(ExactDouble x) : ExactDouble(x) {}
};
/** @brief Convert to the positive representation. */
inline Positive<ExactDouble> cast_positive(ExactDouble x) { return Positive<ExactDouble>(x); }

inline ApproximateDouble::ApproximateDouble(const ExactDouble& x) : _d(x.get_d()) {}

} // namespace Ariadne

#endif
