/***************************************************************************
 *            numeric/lower_number.hpp
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

/** @file numeric/lower_number.hpp
 *  @brief Generic lower-numbers
 */

#ifndef ARIADNE_LOWER_NUMBER_HPP
#define ARIADNE_LOWER_NUMBER_HPP

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
class LowerNumber;
/** @brief Compile-time predicate testing whether a type is number. */
template<class P>
struct IsNumber<LowerNumber<P>> : True {};

template<class P>
Positive<LowerNumber<P>> cast_positive(LowerNumber<P> y);

//! \ingroup NumericModule
//! \brief Generic lower (real) numbers with computational paradigm \a P, which may be %EffectiveTag or %ValidatedTag.
template<class P>
class LowerNumber {
    static_assert(Same<P, EffectiveTag> or Same<P, ValidatedTag>, "P must be a paradigm");
    friend class UpperNumber<P>;

  private:

  public:
    Handle<const NumberInterface> _handle;
    /** @brief Construct a \c LowerNumber from p. */
    explicit LowerNumber(NumberInterface* p) : _handle(p) {}

  private:

  public:
    /** @brief Construct a \c LowerNumber from h. */
    explicit LowerNumber(Handle<const NumberInterface> h) : _handle(h) {}
    /** @brief Return the handle to the represented object. */
    Handle<const NumberInterface> handle() const { return this->_handle; }

  private:
    /** @brief Return a reference to the stored object. */
    const NumberInterface& ref() const { return this->_handle.reference(); }
    /** @brief Return a pointer to the stored object. */
    const NumberInterface* ptr() const { return this->_handle.pointer(); }

  public:
    typedef P Paradigm;
    typedef LowerNumber<P> NumericType;

    /** @brief Construct a \c LowerNumber with the specified characteristics. */
    LowerNumber() : LowerNumber(Integer(0)) {}

    //! \brief Construct from a LowerNumber of a stronger paradigm
    template<StrongerThan<P> SP>
    LowerNumber(const LowerNumber<SP>& y) : LowerNumber(y.handle()) {}
    //! \brief Construct from a type convertible to a Number.
    template<ConvertibleTo<Number<P>> X>
    LowerNumber(const X& x) : LowerNumber(Number<P>(x).handle()) {}

    /** @brief Return the requested element or stored value. */
    template<class PR>
    FloatLowerBound<PR> get(PR pr) const { return this->ref()._get(LowerTag(), pr); }

    /** @brief Extract the requested concrete representation. */
    template<class X>
    X extract() const;

    /** @brief Return the sum or unary positive value. */
    friend LowerNumber<P> operator+(const LowerNumber<P>& y) { return pos(y); }
    /** @brief Return the difference or additive inverse. */
    friend UpperNumber<P> operator-(const LowerNumber<P>& y) { return neg(y); }
    /** @brief Return the sum or unary positive value. */
    friend LowerNumber<P> operator+(const LowerNumber<P>& y1, const LowerNumber<P>& y2) { return add(y1, y2); }
    /** @brief Return the difference or additive inverse. */
    friend LowerNumber<P> operator-(const LowerNumber<P>& y1, const UpperNumber<P>& y2) { return sub(y1, y2); }
    /** @brief Add the operand in place. */
    friend LowerNumber<P>& operator+=(LowerNumber<P>& y1, const LowerNumber<P>& y2) { return y1 = y1 + y2; }
    /** @brief Subtract the operand in place. */
    friend LowerNumber<P>& operator-=(LowerNumber<P>& y1, const UpperNumber<P>& y2) { return y1 = y1 - y2; }

    /** @brief Return the unary positive value. */
    friend LowerNumber<P> pos(const LowerNumber<P>& y) { return LowerNumber<P>(y.ref()._apply(Pos())); }
    /** @brief Return the additive inverse. */
    friend UpperNumber<P> neg(const LowerNumber<P>& y) { return UpperNumber<P>(y.ref()._apply(Neg())); }
    /** @brief Return the sum of the operands. */
    friend LowerNumber<P> add(const LowerNumber<P>& y1, const LowerNumber<P>& y2) { return LowerNumber<P>(y1.ref()._apply(Add(), &y2.ref())); }
    /** @brief Return the difference of the operands. */
    friend LowerNumber<P> sub(const LowerNumber<P>& y1, const UpperNumber<P>& y2) { return LowerNumber<P>(y1.ref()._apply(Sub(), &y2.handle().reference())); }

    /** @brief Return the square root. */
    friend LowerNumber<P> sqrt(const LowerNumber<P>& y) { return LowerNumber<P>(y.ref()._apply(Sqrt())); }
    /** @brief Return the exponential. */
    friend LowerNumber<P> exp(const LowerNumber<P>& y) { return LowerNumber<P>(y.ref()._apply(Exp())); }
    /** @brief Return the natural logarithm. */
    friend LowerNumber<P> log(const LowerNumber<P>& y) { return LowerNumber<P>(y.ref()._apply(Log())); }
    /** @brief Return the inverse tangent. */
    friend LowerNumber<P> atan(const LowerNumber<P>& y) { return LowerNumber<P>(y.ref()._apply(Atan())); }
    /** @brief Return the hyperbolic tangent. */
    friend LowerNumber<P> tanh(const LowerNumber<P>& y) { return LowerNumber<P>(y.ref()._apply(Tanh())); }

    /** @brief Return the maximum of the operands. */
    friend LowerNumber<P> max(const LowerNumber<P>& y1, const LowerNumber<P>& y2) { return LowerNumber<P>(y1.ref()._apply(Max(), &y2.ref())); }
    /** @brief Return the minimum of the operands. */
    friend LowerNumber<P> min(const LowerNumber<P>& y1, const LowerNumber<P>& y2) { return LowerNumber<P>(y1.ref()._apply(Min(), &y2.ref())); }

    /** @brief Test equality of the operands. */
    friend UpperLogicalType<P> operator==(const LowerNumber<P>& y1, const UpperNumber<P>& y2) {
        return Detail::logical_type_from_owned_pointer<EqualityInformation<P>>(y1.ref()._apply(BinaryComparisonOperator(Equal()), y2.handle().pointer()));
    }
    /** @brief Test inequality of the operands. */
    friend LowerLogicalType<P> operator!=(const LowerNumber<P>& y1, const UpperNumber<P>& y2) {
        return not(y1 == y2);
    }
    /** @brief Test whether the first operand is less than the second. */
    friend UpperLogicalType<P> operator<(const LowerNumber<P>& y1, const UpperNumber<P>& y2) {
        return Detail::logical_type_from_owned_pointer<LessThanInformation<P>>(y1.ref()._apply(BinaryComparisonOperator(Less()), y2.handle().pointer()));
    }
    /** @brief Test whether the first operand is greater than the second. */
    friend LowerLogicalType<P> operator>(const LowerNumber<P>& y1, const UpperNumber<P>& y2) {
        return y2 < y1;
    }
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend UpperLogicalType<P> operator<=(const LowerNumber<P>& y1, const UpperNumber<P>& y2) {
        return not(y2 < y1);
    }
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend LowerLogicalType<P> operator>=(const LowerNumber<P>& y1, const UpperNumber<P>& y2) {
        return not(y1 < y2);
    }
    /** @brief Test whether the first operand is less than the second. */
    friend LowerLogicalType<P> operator<(const UpperNumber<P>& y1, const LowerNumber<P>& y2);
    /** @brief Test whether the first operand is greater than the second. */
    friend UpperLogicalType<P> operator>(const UpperNumber<P>& y1, const LowerNumber<P>& y2);
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend LowerLogicalType<P> operator<=(const UpperNumber<P>& y1, const LowerNumber<P>& y2);
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend UpperLogicalType<P> operator>=(const UpperNumber<P>& y1, const LowerNumber<P>& y2);

    /** @brief Return the display name of the type. */
    String class_name() const { return this->ref()._class_name(); }

    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const LowerNumber<P>& y) { return y.ref()._write(os); }

    /** @brief Additive-zero operation token. */
    struct Zero {};
    /** @brief Test whether the first operand is greater than the second. */
    friend LowerLogicalType<P> operator>(const LowerNumber<P>& y, const Zero&) {
        return LowerLogicalType<P>(Detail::logical_type_from_owned_pointer<P>(y.ref()._is_pos()));
    }
    /** @brief Test whether the first operand is less than the second. */
    friend UpperLogicalType<P> operator<(const LowerNumber<P>& y, const Zero&) {
        return not(y > Zero{});
    }
};

/** @brief Convert to the exact representation. */
ExactNumber cast_exact(const ValidatedLowerNumber& y);

/** @brief Nonnegative numeric representation supporting sign-restricted operations. */
template<class P>
class Positive<LowerNumber<P>> : public LowerNumber<P> {
    /** @brief Return the unsigned magnitude of the integer. */
    friend const LowerNumber<P>& unsign(const Positive<LowerNumber<P>>& y) { return y; }

  public:
    /** @brief Construct a \c Positive with the specified characteristics. */
    Positive() : LowerNumber<P>() {}
    /** @brief Construct a \c Positive from y. */
    explicit Positive(const LowerNumber<P>& y) : LowerNumber<P>(y) {}
    /** @brief Construct a \c Positive from n. */
    template<BuiltinUnsignedIntegral N>
    Positive(N n) : LowerNumber<P>(n) {}
    /** @brief Copy-construct the object from the supplied instance. */
    template<class N>
        requires Constructible<ExactNumber, N>
    Positive(const Positive<N>& n) : LowerNumber<P>(ExactNumber(static_cast<const N&>(n))) {}
    /** @brief Construct a \c Positive from n. */
    template<class N>
        requires Constructible<ExactNumber, N> and (not BuiltinIntegral<N>)
    Positive(const N& n) : LowerNumber<P>(ExactNumber(n)) {}
    /** @brief Convert to lower number< p >. */
    explicit operator LowerNumber<P>() const { return *this; }

    /** @brief Return the product of the operands. */
    friend LowerNumber<P> mul(const LowerNumber<P>& y1, const LowerNumber<P>& y2);

    /** @brief Return the sum or unary positive value. */
    friend Positive<LowerNumber<P>> operator+(const Positive<LowerNumber<P>>& y1, const Positive<LowerNumber<P>>& y2) {
        return cast_positive(add(unsign(y1), unsign(y2)));
    }
    /** @brief Return the product. */
    friend Positive<LowerNumber<P>> operator*(const Positive<LowerNumber<P>>& y1, const Positive<LowerNumber<P>>& y2) {
        return cast_positive(mul(unsign(y1), unsign(y2)));
    }
    /** @brief Return the sum of the operands. */
    friend Positive<LowerNumber<P>> add(const Positive<LowerNumber<P>>& y1, const Positive<LowerNumber<P>>& y2) {
        return cast_positive(add(unsign(y1), unsign(y2)));
    }
    /** @brief Return the product of the operands. */
    friend Positive<LowerNumber<P>> mul(const Positive<LowerNumber<P>>& y1, const Positive<LowerNumber<P>>& y2) {
        return cast_positive(Number<P>(y1.ref()._apply(Mul(), y2.ref())));
    }
    /** @brief Return the maximum of the operands. */
    friend Positive<LowerNumber<P>> max(const Positive<LowerNumber<P>>& y1, const Positive<LowerNumber<P>>& y2) {
        return cast_positive(max(unsign(y1), unsign(y2)));
    }
};

/** @brief Convert to the positive representation. */
template<class P>
Positive<LowerNumber<P>> cast_positive(LowerNumber<P> y) {
    return Positive<LowerNumber<P>>(y);
}
template<class P>
Positive<ExactNumber> cast_exact(const Positive<LowerNumber<P>>& y);

} // namespace Ariadne

#endif
