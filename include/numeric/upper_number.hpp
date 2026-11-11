/***************************************************************************
 *            numeric/upper_number.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/** @file numeric/upper_number.hpp
 *  @brief Generic upper-numbers
 */

#ifndef ARIADNE_UPPER_NUMBER_HPP
#define ARIADNE_UPPER_NUMBER_HPP

#include "arithmetic.hpp"
#include "dyadic.hpp"
#include "float.decl.hpp"
#include "foundation/logical.decl.hpp"
#include "foundation/logical.hpp"
#include "foundation/paradigm.hpp"
#include "integer.hpp"
#include "number.decl.hpp"
#include "number_interface.hpp"
#include "rational.hpp"
#include "real.hpp"
#include "utility/handle.hpp"
#include "utility/prototype.hpp"

namespace Ariadne {

template<class X>
struct IsNumber;
template<class P>
class UpperNumber;
/** @brief Compile-time predicate testing whether a type is number. */
template<class P>
struct IsNumber<UpperNumber<P>> : True {};

template<class P>
Positive<UpperNumber<P>> cast_positive(UpperNumber<P> y);

//! \ingroup NumericModule
//! \brief Generic upper (real) numbers with computational paradigm \a P, which may be %EffectiveTag or %ValidatedTag.
template<class P>
class UpperNumber {
    static_assert(Same<P, EffectiveTag> or Same<P, ValidatedTag>, "P must be a paradigm");
    friend class LowerNumber<P>;

  private:

  public:
    Handle<const NumberInterface> _handle;
    /** @brief Construct a \c UpperNumber from p. */
    explicit UpperNumber(NumberInterface* p) : _handle(p) {}

  private:

  public:
    /** @brief Construct a \c UpperNumber from h. */
    explicit UpperNumber(Handle<const NumberInterface> h) : _handle(h) {}
    /** @brief Return the handle to the represented object. */
    Handle<const NumberInterface> handle() const { return this->_handle; }

  private:
    /** @brief Return a reference to the stored object. */
    const NumberInterface& ref() const { return this->_handle.reference(); }
    /** @brief Return a pointer to the stored object. */
    const NumberInterface* ptr() const { return this->_handle.pointer(); }

  public:
    typedef P Paradigm;
    typedef UpperNumber<P> NumericType;

    /** @brief Construct a \c UpperNumber with the specified characteristics. */
    UpperNumber() : UpperNumber(Integer(0)) {}

    //! \brief Construct from a UpperNumber of a stronger paradigm
    template<StrongerThan<P> SP>
    UpperNumber(const UpperNumber<SP>& y) : UpperNumber(y.handle()) {}
    //! \brief Construct from a type convertible to a Number.
    template<ConvertibleTo<Number<P>> X>
    UpperNumber(const X& x) : UpperNumber(Number<P>(x).handle()) {}

    /** @brief Return the requested element or stored value. */
    template<class PR>
    FloatUpperBound<PR> get(PR pr) const { return this->ref()._get(UpperTag(), pr); }

    /** @brief Extract the requested concrete representation. */
    template<class X>
    X extract() const;

    /** @brief Return the sum or unary positive value. */
    friend UpperNumber<P> operator+(const UpperNumber<P>& y) { return pos(y); }
    /** @brief Return the difference or additive inverse. */
    friend LowerNumber<P> operator-(const UpperNumber<P>& y) { return neg(y); }
    /** @brief Return the sum or unary positive value. */
    friend UpperNumber<P> operator+(const UpperNumber<P>& y1, const UpperNumber<P>& y2) { return add(y1, y2); }
    /** @brief Return the difference or additive inverse. */
    friend UpperNumber<P> operator-(const UpperNumber<P>& y1, const LowerNumber<P>& y2) { return sub(y1, y2); }
    /** @brief Add the operand in place. */
    friend UpperNumber<P>& operator+=(UpperNumber<P>& y1, const UpperNumber<P>& y2) { return y1 = y1 + y2; }
    /** @brief Subtract the operand in place. */
    friend UpperNumber<P>& operator-=(UpperNumber<P>& y1, const LowerNumber<P>& y2) { return y1 = y1 - y2; }

    /** @brief Return the unary positive value. */
    friend UpperNumber<P> pos(const UpperNumber<P>& y) { return UpperNumber<P>(y.ref()._apply(Pos())); }
    /** @brief Return the additive inverse. */
    friend LowerNumber<P> neg(const UpperNumber<P>& y) { return LowerNumber<P>(y.ref()._apply(Neg())); }
    /** @brief Return the sum of the operands. */
    friend UpperNumber<P> add(const UpperNumber<P>& y1, const UpperNumber<P>& y2) { return UpperNumber<P>(y1.ref()._apply(Add(), &y2.ref())); }
    /** @brief Return the difference of the operands. */
    friend UpperNumber<P> sub(const UpperNumber<P>& y1, const LowerNumber<P>& y2) { return UpperNumber<P>(y1.ref()._apply(Sub(), &y2.handle().reference())); }

    /** @brief Return the square root. */
    friend UpperNumber<P> sqrt(const UpperNumber<P>& y) { return UpperNumber<P>(y.ref()._apply(Sqrt())); }
    /** @brief Return the exponential. */
    friend UpperNumber<P> exp(const UpperNumber<P>& y) { return UpperNumber<P>(y.ref()._apply(Exp())); }
    /** @brief Return the natural logarithm. */
    friend UpperNumber<P> log(const UpperNumber<P>& y) { return UpperNumber<P>(y.ref()._apply(Log())); }
    /** @brief Return the inverse tangent. */
    friend UpperNumber<P> atan(const UpperNumber<P>& y) { return UpperNumber<P>(y.ref()._apply(Atan())); }
    /** @brief Return the hyperbolic tangent. */
    friend UpperNumber<P> tanh(const UpperNumber<P>& y) { return UpperNumber<P>(y.ref()._apply(Tanh())); }

    /** @brief Return the maximum of the operands. */
    friend UpperNumber<P> max(const UpperNumber<P>& y1, const UpperNumber<P>& y2) { return UpperNumber<P>(y1.ref()._apply(Max(), &y2.ref())); }
    /** @brief Return the minimum of the operands. */
    friend UpperNumber<P> min(const UpperNumber<P>& y1, const UpperNumber<P>& y2) { return UpperNumber<P>(y1.ref()._apply(Min(), &y2.ref())); }

    /** @brief Test equality of the operands. */
    friend UpperLogicalType<P> operator==(const UpperNumber<P>& y1, const LowerNumber<P>& y2) {
        return Detail::logical_type_from_owned_pointer<EqualityInformation<P>>(y1.ref()._apply(BinaryComparisonOperator(Equal()), y2.handle().pointer()));
    }
    /** @brief Test inequality of the operands. */
    friend LowerLogicalType<P> operator!=(const UpperNumber<P>& y1, const LowerNumber<P>& y2) {
        return not(y1 == y2);
    }
    /** @brief Test whether the first operand is less than the second. */
    friend LowerLogicalType<P> operator<(const UpperNumber<P>& y1, const LowerNumber<P>& y2) {
        return Detail::logical_type_from_owned_pointer<LessThanInformation<P>>(y1.ref()._apply(BinaryComparisonOperator(Less()), y2.handle().pointer()));
    }
    /** @brief Test whether the first operand is greater than the second. */
    friend UpperLogicalType<P> operator>(const UpperNumber<P>& y1, const LowerNumber<P>& y2) {
        return y2 < y1;
    }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend LowerLogicalType<P> operator<=(const UpperNumber<P>& y1, const LowerNumber<P>& y2) {
        return not(y2 < y1);
    }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend UpperLogicalType<P> operator>=(const UpperNumber<P>& y1, const LowerNumber<P>& y2) {
        return not(y1 < y2);
    }
    /** @brief Test whether the first operand is less than the second. */
    friend UpperLogicalType<P> operator<(const LowerNumber<P>& y1, const UpperNumber<P>& y2);
    /** @brief Test whether the first operand is greater than the second. */
    friend LowerLogicalType<P> operator>(const LowerNumber<P>& y1, const UpperNumber<P>& y2);
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend UpperLogicalType<P> operator<=(const LowerNumber<P>& y1, const UpperNumber<P>& y2);
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend LowerLogicalType<P> operator>=(const LowerNumber<P>& y1, const UpperNumber<P>& y2);

    /** @brief Return the display name of the type. */
    String class_name() const { return this->ref()._class_name(); }

    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const UpperNumber<P>& y) { return y.ref()._write(os); }
};

/** @brief Convert to the exact representation. */
ExactNumber cast_exact(const ValidatedUpperNumber& y);

/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<class P>
class Positive<UpperNumber<P>> : public UpperNumber<P> {
    /** @brief Return the unsigned magnitude of the integer. */
    friend const UpperNumber<P>& unsign(const Positive<UpperNumber<P>>& y) { return y; }

  public:
    /** @brief Construct a \c Positive with the specified characteristics. */
    Positive() : UpperNumber<P>() {}
    /** @brief Construct a \c Positive from y. */
    explicit Positive(const UpperNumber<P>& y) : UpperNumber<P>(y) {}
    /** @brief Construct a \c Positive from n. */
    template<BuiltinUnsignedIntegral N>
    Positive(N n) : UpperNumber<P>(n) {}
    /** @brief Copy-construct the object from the supplied instance. */
    template<class N>
        requires Constructible<ExactNumber, N>
    Positive(const Positive<N>& n) : UpperNumber<P>(ExactNumber(static_cast<const N&>(n))) {}
    /** @brief Construct a \c Positive from n. */
    template<class N>
        requires Constructible<ExactNumber, N> and (not BuiltinIntegral<N>)
    Positive(const N& n) : UpperNumber<P>(ExactNumber(n)) {}
    /** @brief Convert to upper number< p >. */
    explicit operator UpperNumber<P>() const { return *this; }

    /** @brief Return the product of the operands. */
    friend UpperNumber<P> mul(const UpperNumber<P>& y1, const UpperNumber<P>& y2);

    /** @brief Return the sum or unary positive value. */
    friend Positive<UpperNumber<P>> operator+(const Positive<UpperNumber<P>>& y1, const Positive<UpperNumber<P>>& y2) {
        return cast_positive(add(unsign(y1), unsign(y2)));
    }
    /** @brief Return the product. */
    friend Positive<UpperNumber<P>> operator*(const Positive<UpperNumber<P>>& y1, const Positive<UpperNumber<P>>& y2) {
        return cast_positive(mul(unsign(y1), unsign(y2)));
    }
    /** @brief Return the sum of the operands. */
    friend Positive<UpperNumber<P>> add(const Positive<UpperNumber<P>>& y1, const Positive<UpperNumber<P>>& y2) {
        return cast_positive(add(unsign(y1), unsign(y2)));
    }
    /** @brief Return the product of the operands. */
    friend Positive<UpperNumber<P>> mul(const Positive<UpperNumber<P>>& y1, const Positive<UpperNumber<P>>& y2) {
        return cast_positive(Number<P>(y1.ref()._apply(Mul(), y2.ref())));
    }
    /** @brief Return the maximum of the operands. */
    friend Positive<UpperNumber<P>> max(const Positive<UpperNumber<P>>& y1, const Positive<UpperNumber<P>>& y2) {
        return cast_positive(max(unsign(y1), unsign(y2)));
    }
};

/** @brief Convert to the positive representation. */
template<class P>
Positive<UpperNumber<P>> cast_positive(UpperNumber<P> y) {
    return Positive<UpperNumber<P>>(y);
}
template<class P>
Positive<ExactNumber> cast_exact(const Positive<UpperNumber<P>>& y);

} // namespace Ariadne

#endif
