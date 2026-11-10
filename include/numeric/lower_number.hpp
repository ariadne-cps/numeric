/***************************************************************************
 *            numeric/lower_number.hpp
 *
 *  Copyright  2013-20  Pieter Collins
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

/*! \file numeric/lower_number.hpp
 *  \brief Generic lower-numbers
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
    explicit LowerNumber(NumberInterface* p) : _handle(p) {}

  private:

  public:
    explicit LowerNumber(Handle<const NumberInterface> h) : _handle(h) {}
    Handle<const NumberInterface> handle() const { return this->_handle; }

  private:
    const NumberInterface& ref() const { return this->_handle.reference(); }
    const NumberInterface* ptr() const { return this->_handle.pointer(); }

  public:
    typedef P Paradigm;
    typedef LowerNumber<P> NumericType;

    LowerNumber() : LowerNumber(Integer(0)) {}

    //! \brief Construct from a LowerNumber of a stronger paradigm
    template<StrongerThan<P> SP>
    LowerNumber(const LowerNumber<SP>& y) : LowerNumber(y.handle()) {}
    //! \brief Construct from a type convertible to a Number.
    template<ConvertibleTo<Number<P>> X>
    LowerNumber(const X& x) : LowerNumber(Number<P>(x).handle()) {}

    template<class PR>
    FloatLowerBound<PR> get(PR pr) const { return this->ref()._get(LowerTag(), pr); }

    template<class X>
    X extract() const;

    friend LowerNumber<P> operator+(const LowerNumber<P>& y) { return pos(y); }
    friend UpperNumber<P> operator-(const LowerNumber<P>& y) { return neg(y); }
    friend LowerNumber<P> operator+(const LowerNumber<P>& y1, const LowerNumber<P>& y2) { return add(y1, y2); }
    friend LowerNumber<P> operator-(const LowerNumber<P>& y1, const UpperNumber<P>& y2) { return sub(y1, y2); }
    friend LowerNumber<P>& operator+=(LowerNumber<P>& y1, const LowerNumber<P>& y2) { return y1 = y1 + y2; }
    friend LowerNumber<P>& operator-=(LowerNumber<P>& y1, const UpperNumber<P>& y2) { return y1 = y1 - y2; }

    friend LowerNumber<P> pos(const LowerNumber<P>& y) { return LowerNumber<P>(y.ref()._apply(Pos())); }
    friend UpperNumber<P> neg(const LowerNumber<P>& y) { return UpperNumber<P>(y.ref()._apply(Neg())); }
    friend LowerNumber<P> add(const LowerNumber<P>& y1, const LowerNumber<P>& y2) { return LowerNumber<P>(y1.ref()._apply(Add(), &y2.ref())); }
    friend LowerNumber<P> sub(const LowerNumber<P>& y1, const UpperNumber<P>& y2) { return LowerNumber<P>(y1.ref()._apply(Sub(), &y2.handle().reference())); }

    friend LowerNumber<P> sqrt(const LowerNumber<P>& y) { return LowerNumber<P>(y.ref()._apply(Sqrt())); }
    friend LowerNumber<P> exp(const LowerNumber<P>& y) { return LowerNumber<P>(y.ref()._apply(Exp())); }
    friend LowerNumber<P> log(const LowerNumber<P>& y) { return LowerNumber<P>(y.ref()._apply(Log())); }
    friend LowerNumber<P> atan(const LowerNumber<P>& y) { return LowerNumber<P>(y.ref()._apply(Atan())); }
    friend LowerNumber<P> tanh(const LowerNumber<P>& y) { return LowerNumber<P>(y.ref()._apply(Tanh())); }

    friend LowerNumber<P> max(const LowerNumber<P>& y1, const LowerNumber<P>& y2) { return LowerNumber<P>(y1.ref()._apply(Max(), &y2.ref())); }
    friend LowerNumber<P> min(const LowerNumber<P>& y1, const LowerNumber<P>& y2) { return LowerNumber<P>(y1.ref()._apply(Min(), &y2.ref())); }

    friend UpperLogicalType<P> operator==(const LowerNumber<P>& y1, const UpperNumber<P>& y2) {
        return Detail::logical_type_from_owned_pointer<EqualityInformation<P>>(y1.ref()._apply(BinaryComparisonOperator(Equal()), y2.handle().pointer()));
    }
    friend LowerLogicalType<P> operator!=(const LowerNumber<P>& y1, const UpperNumber<P>& y2) {
        return not(y1 == y2);
    }
    friend UpperLogicalType<P> operator<(const LowerNumber<P>& y1, const UpperNumber<P>& y2) {
        return Detail::logical_type_from_owned_pointer<LessThanInformation<P>>(y1.ref()._apply(BinaryComparisonOperator(Less()), y2.handle().pointer()));
    }
    friend LowerLogicalType<P> operator>(const LowerNumber<P>& y1, const UpperNumber<P>& y2) {
        return y2 < y1;
    }
    friend UpperLogicalType<P> operator<=(const LowerNumber<P>& y1, const UpperNumber<P>& y2) {
        return not(y2 < y1);
    }
    friend LowerLogicalType<P> operator>=(const LowerNumber<P>& y1, const UpperNumber<P>& y2) {
        return not(y1 < y2);
    }
    friend LowerLogicalType<P> operator<(const UpperNumber<P>& y1, const LowerNumber<P>& y2);
    friend UpperLogicalType<P> operator>(const UpperNumber<P>& y1, const LowerNumber<P>& y2);
    friend LowerLogicalType<P> operator<=(const UpperNumber<P>& y1, const LowerNumber<P>& y2);
    friend UpperLogicalType<P> operator>=(const UpperNumber<P>& y1, const LowerNumber<P>& y2);

    String class_name() const { return this->ref()._class_name(); }

    friend OutputStream& operator<<(OutputStream& os, const LowerNumber<P>& y) { return y.ref()._write(os); }

    struct Zero {};
    friend LowerLogicalType<P> operator>(const LowerNumber<P>& y, const Zero&) {
        return LowerLogicalType<P>(Detail::logical_type_from_owned_pointer<P>(y.ref()._is_pos()));
    }
    friend UpperLogicalType<P> operator<(const LowerNumber<P>& y, const Zero&) {
        return not(y > Zero{});
    }
};

ExactNumber cast_exact(const ValidatedLowerNumber& y);

template<class P>
class Positive<LowerNumber<P>> : public LowerNumber<P> {
    friend const LowerNumber<P>& unsign(const Positive<LowerNumber<P>>& y) { return y; }

  public:
    Positive() : LowerNumber<P>() {}
    explicit Positive(const LowerNumber<P>& y) : LowerNumber<P>(y) {}
    template<BuiltinUnsignedIntegral N>
    Positive(N n) : LowerNumber<P>(n) {}
    template<class N>
        requires Constructible<ExactNumber, N>
    Positive(const Positive<N>& n) : LowerNumber<P>(ExactNumber(static_cast<const N&>(n))) {}
    template<class N>
        requires Constructible<ExactNumber, N> and (not BuiltinIntegral<N>)
    Positive(const N& n) : LowerNumber<P>(ExactNumber(n)) {}
    explicit operator LowerNumber<P>() const { return *this; }

    friend LowerNumber<P> mul(const LowerNumber<P>& y1, const LowerNumber<P>& y2);

    friend Positive<LowerNumber<P>> operator+(const Positive<LowerNumber<P>>& y1, const Positive<LowerNumber<P>>& y2) {
        return cast_positive(add(unsign(y1), unsign(y2)));
    }
    friend Positive<LowerNumber<P>> operator*(const Positive<LowerNumber<P>>& y1, const Positive<LowerNumber<P>>& y2) {
        return cast_positive(mul(unsign(y1), unsign(y2)));
    }
    friend Positive<LowerNumber<P>> add(const Positive<LowerNumber<P>>& y1, const Positive<LowerNumber<P>>& y2) {
        return cast_positive(add(unsign(y1), unsign(y2)));
    }
    friend Positive<LowerNumber<P>> mul(const Positive<LowerNumber<P>>& y1, const Positive<LowerNumber<P>>& y2) {
        return cast_positive(Number<P>(y1.ref()._apply(Mul(), y2.ref())));
    }
    friend Positive<LowerNumber<P>> max(const Positive<LowerNumber<P>>& y1, const Positive<LowerNumber<P>>& y2) {
        return cast_positive(max(unsign(y1), unsign(y2)));
    }
};

template<class P>
Positive<LowerNumber<P>> cast_positive(LowerNumber<P> y) {
    return Positive<LowerNumber<P>>(y);
}
template<class P>
Positive<ExactNumber> cast_exact(const Positive<LowerNumber<P>>& y);

} // namespace Ariadne

#endif
