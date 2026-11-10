/***************************************************************************
 *            numeric/upper_number.hpp
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

/*! \file numeric/upper_number.hpp
 *  \brief Generic upper-numbers
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
    explicit UpperNumber(NumberInterface* p) : _handle(p) {}

  private:

  public:
    explicit UpperNumber(Handle<const NumberInterface> h) : _handle(h) {}
    Handle<const NumberInterface> handle() const { return this->_handle; }

  private:
    const NumberInterface& ref() const { return this->_handle.reference(); }
    const NumberInterface* ptr() const { return this->_handle.pointer(); }

  public:
    typedef P Paradigm;
    typedef UpperNumber<P> NumericType;

    UpperNumber() : UpperNumber(Integer(0)) {}

    //! \brief Construct from a UpperNumber of a stronger paradigm
    template<StrongerThan<P> SP>
    UpperNumber(const UpperNumber<SP>& y) : UpperNumber(y.handle()) {}
    //! \brief Construct from a type convertible to a Number.
    template<ConvertibleTo<Number<P>> X>
    UpperNumber(const X& x) : UpperNumber(Number<P>(x).handle()) {}

    template<class PR>
    FloatUpperBound<PR> get(PR pr) const { return this->ref()._get(UpperTag(), pr); }

    template<class X>
    X extract() const;

    friend UpperNumber<P> operator+(const UpperNumber<P>& y) { return pos(y); }
    friend LowerNumber<P> operator-(const UpperNumber<P>& y) { return neg(y); }
    friend UpperNumber<P> operator+(const UpperNumber<P>& y1, const UpperNumber<P>& y2) { return add(y1, y2); }
    friend UpperNumber<P> operator-(const UpperNumber<P>& y1, const LowerNumber<P>& y2) { return sub(y1, y2); }
    friend UpperNumber<P>& operator+=(UpperNumber<P>& y1, const UpperNumber<P>& y2) { return y1 = y1 + y2; }
    friend UpperNumber<P>& operator-=(UpperNumber<P>& y1, const LowerNumber<P>& y2) { return y1 = y1 - y2; }

    friend UpperNumber<P> pos(const UpperNumber<P>& y) { return UpperNumber<P>(y.ref()._apply(Pos())); }
    friend LowerNumber<P> neg(const UpperNumber<P>& y) { return LowerNumber<P>(y.ref()._apply(Neg())); }
    friend UpperNumber<P> add(const UpperNumber<P>& y1, const UpperNumber<P>& y2) { return UpperNumber<P>(y1.ref()._apply(Add(), &y2.ref())); }
    friend UpperNumber<P> sub(const UpperNumber<P>& y1, const LowerNumber<P>& y2) { return UpperNumber<P>(y1.ref()._apply(Sub(), &y2.handle().reference())); }

    friend UpperNumber<P> sqrt(const UpperNumber<P>& y) { return UpperNumber<P>(y.ref()._apply(Sqrt())); }
    friend UpperNumber<P> exp(const UpperNumber<P>& y) { return UpperNumber<P>(y.ref()._apply(Exp())); }
    friend UpperNumber<P> log(const UpperNumber<P>& y) { return UpperNumber<P>(y.ref()._apply(Log())); }
    friend UpperNumber<P> atan(const UpperNumber<P>& y) { return UpperNumber<P>(y.ref()._apply(Atan())); }
    friend UpperNumber<P> tanh(const UpperNumber<P>& y) { return UpperNumber<P>(y.ref()._apply(Tanh())); }

    friend UpperNumber<P> max(const UpperNumber<P>& y1, const UpperNumber<P>& y2) { return UpperNumber<P>(y1.ref()._apply(Max(), &y2.ref())); }
    friend UpperNumber<P> min(const UpperNumber<P>& y1, const UpperNumber<P>& y2) { return UpperNumber<P>(y1.ref()._apply(Min(), &y2.ref())); }

    friend UpperLogicalType<P> operator==(const UpperNumber<P>& y1, const LowerNumber<P>& y2) {
        return Detail::logical_type_from_owned_pointer<EqualityInformation<P>>(y1.ref()._apply(BinaryComparisonOperator(Equal()), y2.handle().pointer()));
    }
    friend LowerLogicalType<P> operator!=(const UpperNumber<P>& y1, const LowerNumber<P>& y2) {
        return not(y1 == y2);
    }
    friend LowerLogicalType<P> operator<(const UpperNumber<P>& y1, const LowerNumber<P>& y2) {
        return Detail::logical_type_from_owned_pointer<LessThanInformation<P>>(y1.ref()._apply(BinaryComparisonOperator(Less()), y2.handle().pointer()));
    }
    friend UpperLogicalType<P> operator>(const UpperNumber<P>& y1, const LowerNumber<P>& y2) {
        return y2 < y1;
    }
    friend LowerLogicalType<P> operator<=(const UpperNumber<P>& y1, const LowerNumber<P>& y2) {
        return not(y2 < y1);
    }
    friend UpperLogicalType<P> operator>=(const UpperNumber<P>& y1, const LowerNumber<P>& y2) {
        return not(y1 < y2);
    }
    friend UpperLogicalType<P> operator<(const LowerNumber<P>& y1, const UpperNumber<P>& y2);
    friend LowerLogicalType<P> operator>(const LowerNumber<P>& y1, const UpperNumber<P>& y2);
    friend UpperLogicalType<P> operator<=(const LowerNumber<P>& y1, const UpperNumber<P>& y2);
    friend LowerLogicalType<P> operator>=(const LowerNumber<P>& y1, const UpperNumber<P>& y2);

    String class_name() const { return this->ref()._class_name(); }

    friend OutputStream& operator<<(OutputStream& os, const UpperNumber<P>& y) { return y.ref()._write(os); }
};

ExactNumber cast_exact(const ValidatedUpperNumber& y);

template<class P>
class Positive<UpperNumber<P>> : public UpperNumber<P> {
    friend const UpperNumber<P>& unsign(const Positive<UpperNumber<P>>& y) { return y; }

  public:
    Positive() : UpperNumber<P>() {}
    explicit Positive(const UpperNumber<P>& y) : UpperNumber<P>(y) {}
    template<BuiltinUnsignedIntegral N>
    Positive(N n) : UpperNumber<P>(n) {}
    template<class N>
        requires Constructible<ExactNumber, N>
    Positive(const Positive<N>& n) : UpperNumber<P>(ExactNumber(static_cast<const N&>(n))) {}
    template<class N>
        requires Constructible<ExactNumber, N> and (not BuiltinIntegral<N>)
    Positive(const N& n) : UpperNumber<P>(ExactNumber(n)) {}
    explicit operator UpperNumber<P>() const { return *this; }

    friend UpperNumber<P> mul(const UpperNumber<P>& y1, const UpperNumber<P>& y2);

    friend Positive<UpperNumber<P>> operator+(const Positive<UpperNumber<P>>& y1, const Positive<UpperNumber<P>>& y2) {
        return cast_positive(add(unsign(y1), unsign(y2)));
    }
    friend Positive<UpperNumber<P>> operator*(const Positive<UpperNumber<P>>& y1, const Positive<UpperNumber<P>>& y2) {
        return cast_positive(mul(unsign(y1), unsign(y2)));
    }
    friend Positive<UpperNumber<P>> add(const Positive<UpperNumber<P>>& y1, const Positive<UpperNumber<P>>& y2) {
        return cast_positive(add(unsign(y1), unsign(y2)));
    }
    friend Positive<UpperNumber<P>> mul(const Positive<UpperNumber<P>>& y1, const Positive<UpperNumber<P>>& y2) {
        return cast_positive(Number<P>(y1.ref()._apply(Mul(), y2.ref())));
    }
    friend Positive<UpperNumber<P>> max(const Positive<UpperNumber<P>>& y1, const Positive<UpperNumber<P>>& y2) {
        return cast_positive(max(unsign(y1), unsign(y2)));
    }
};

template<class P>
Positive<UpperNumber<P>> cast_positive(UpperNumber<P> y) {
    return Positive<UpperNumber<P>>(y);
}
template<class P>
Positive<ExactNumber> cast_exact(const Positive<UpperNumber<P>>& y);

} // namespace Ariadne

#endif
