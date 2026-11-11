/***************************************************************************
 *            numeric/validated_real.hpp
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

/** @file numeric/validated_real.hpp
 *  @brief Validated real numbers
 */

#ifndef ARIADNE_VALIDATED_REAL_HPP
#define ARIADNE_VALIDATED_REAL_HPP

#include "foundation/logical.decl.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/number.decl.hpp"
#include "utility/handle.hpp"
#include "utility/pointer.hpp"
#include "utility/typedefs.hpp"

namespace Ariadne {

class ValidatedRealInterface;

//! \ingroup NumericModule
//! \brief A generic class representing rigorous bounds on a real number.
//! \see Real
class ValidatedReal
    : public Handle<const ValidatedRealInterface> {
  public:
    typedef ValidatedRealInterface Interface;

    /** @brief Construct a \c ValidatedReal with the specified characteristics. */
    ValidatedReal(const DyadicBounds&);
    /** @brief Convert to dyadic bounds. */
    operator DyadicBounds() const;

    /** @brief Return the stored value. */
    Dyadic value();
    /** @brief Return the error bound. */
    Dyadic error();
    /** @brief Return the lower endpoint. */
    Dyadic lower();
    /** @brief Return the upper endpoint. */
    Dyadic upper();

    //! \brief Get dyadic bounds for the number.
    //! \details It is not always clear how this function can be implemented,
    //!   and it may be removed or modified in future versions.
    DyadicBounds get() const;
    //! \brief Get the bounds for the number, representing in double precision.
    FloatDPBounds get(DoublePrecision) const;
    //! \brief Get the bounds for the number, representing to precision \a pr.
    //! Note that increasing the accuracy typically does not yield arbitrarily
    //! tight bounds, as the object is already an approximation to finite accuracy.
    FloatMPBounds get(MultiplePrecision pr) const;
    //! \brief Write to an output stream.
    friend OutputStream& operator<<(OutputStream&, const ValidatedReal&);

    /** @brief Return the sign of the value. */
    friend ValidatedKleenean sgn(const DyadicBounds& r);
    /** @brief Test equality of the operands. */
    friend ValidatedKleenean operator==(const DyadicBounds& vr1, const DyadicBounds& vr2);
    /** @brief Test inequality of the operands. */
    friend ValidatedKleenean operator!=(const DyadicBounds& vr1, const DyadicBounds& vr2);
    /** @brief Test whether the first operand is less than or equal to the second. */
    friend ValidatedKleenean operator<=(const DyadicBounds& vr1, const DyadicBounds& vr2);
    /** @brief Test whether the first operand is greater than or equal to the second. */
    friend ValidatedKleenean operator>=(const DyadicBounds& vr1, const DyadicBounds& vr2);
    /** @brief Test whether the first operand is less than the second. */
    friend ValidatedKleenean operator<(const DyadicBounds& vr1, const DyadicBounds& vr2);
    /** @brief Test whether the first operand is greater than the second. */
    friend ValidatedKleenean operator>(const DyadicBounds& vr1, const DyadicBounds& vr2);
};

} // namespace Ariadne

#endif
