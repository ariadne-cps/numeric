/***************************************************************************
 *            numeric/field.hpp
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

/** @file numeric/field.hpp
 *  @brief Abstractions for ring and field (numeric) types.
 */

#ifndef ARIADNE_FIELD_HPP
#define ARIADNE_FIELD_HPP

namespace Ariadne {

//! \brief Interface for a ring.
template<class T>
class RingInterface
    : public virtual WritableInterface {
  public:
    //! \brief Virtual destructor.
    virtual ~RingInterface<T>() = default;
    //! \brief Create a dynamically-allocated copy.
    inline T* _create_copy() const { return this->_pos(); }
    //! \brief Create the zero element in the same algebra as the current object.
    inline T* _create_zero() const { return this->_nul(); }

    //! \brief The square of the element.
    virtual T* _sqr() const = 0;

    /** @brief Return a copy of the object. */
    inline T* _copy() const { return this->_create_copy(); }
};

//! \brief Interface for a field.
template<class T>
class FieldInterface
    : public virtual RingInterface<T> {
  public:
    /** @brief Return the quotient of the operands. */
    virtual T* _div(const T&) const = 0;
    /** @brief Return the multiplicative inverse. */
    virtual T* _rec() const = 0;
    /** @brief Raise the value to the specified power. */
    virtual T* _pow(Int n) const = 0;
};

//! \brief Interface for a field supporting transcendental operations.
template<class T>
class TranscendentalInterface
    : public virtual FieldInterface<T> {
  public:
    /** @brief Return the square root. */
    virtual T* _sqrt() const = 0;
    /** @brief Return the exponential. */
    virtual T* _exp() const = 0;
    /** @brief Return the natural logarithm. */
    virtual T* _log() const = 0;
    /** @brief Return the sine. */
    virtual T* _sin() const = 0;
    /** @brief Return the cosine. */
    virtual T* _cos() const = 0;
    /** @brief Return the tangent. */
    virtual T* _tan() const = 0;
    /** @brief Return the inverse tangent. */
    virtual T* _atan() const = 0;
};

} // namespace Ariadne

#endif
