/***************************************************************************
 *            numeric/number_interface.hpp
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

/** @file numeric/number_interface.hpp
 *  @brief Declarations and utilities for number interface in Ariadne Numeric.
 */

#ifndef ARIADNE_NUMBER_INTERFACE
#define ARIADNE_NUMBER_INTERFACE

#include "float.decl.hpp"
#include "number.decl.hpp"
#include "operators.hpp"
#include "utility/clonable.hpp"
#include "utility/writable.hpp"

namespace Ariadne {

/************ Number *********************************************************/

class NumberInterface;
class UnaryOperatorInterface;
template<class... YS>
struct Aware;
template<class X>
class NumberWrapper;

namespace Detail {
class LogicalInterface;
}
using Detail::LogicalInterface;

/** @brief Polymorphic interface for number. */
class NumberInterface
    : public virtual WritableInterface,
      public virtual ClonableInterface
// Remove 'enable_shared_from_this' due to "relocation against undefined hidden symbol when making a shared object" error in g++
//    , public std::enable_shared_from_this<NumberInterface>
{
    template<class X>
    friend class NumberWrapper;
    friend class Handle<NumberInterface>;

  public:
    /** @brief Destroy the object and release its owned resources. */
    virtual ~NumberInterface() = default;

  public:
    /** @brief Return a copy of the object. */
    virtual NumberInterface* _copy() const = 0;
    /** @brief Return a dynamically allocated object by moving the stored value. */
    virtual NumberInterface* _move() = 0;

    /** @brief Apply the supplied operation to the operands. */
    virtual NumberInterface* _apply(UnaryElementaryOperator op) const = 0;
    /** @brief Apply the supplied operation to the operands. */
    virtual NumberInterface* _apply(BinaryElementaryOperator op, const NumberInterface* y) const = 0;
    /** @brief Apply the operation with reversed operand order. */
    virtual NumberInterface* _rapply(BinaryElementaryOperator op, const NumberInterface* y) const = 0;
    /** @brief Apply the supplied operation to the operands. */
    virtual NumberInterface* _apply(GradedElementaryOperator op, Int n) const = 0;

    /** @brief Apply the supplied operation to the operands. */
    virtual LogicalInterface* _apply(BinaryComparisonOperator op, const NumberInterface* y) const = 0;
    /** @brief Apply the operation with reversed operand order. */
    virtual LogicalInterface* _rapply(BinaryComparisonOperator op, const NumberInterface* y) const = 0;

    /** @brief Return the q. */
    virtual Rational _get_q() const = 0;

    /** @brief Return the requested element or stored value. */
    virtual FloatDPBall _get(MetricTag, DoublePrecision, DoublePrecision) const = 0;
    /** @brief Return the requested element or stored value. */
    virtual FloatDPBounds _get(OrderTag, DoublePrecision) const = 0;
    /** @brief Return the requested element or stored value. */
    virtual FloatDPUpperBound _get(UpperTag, DoublePrecision) const = 0;
    /** @brief Return the requested element or stored value. */
    virtual FloatDPLowerBound _get(LowerTag, DoublePrecision) const = 0;
    /** @brief Return the requested element or stored value. */
    virtual FloatDPApproximation _get(ApproximateTag, DoublePrecision) const = 0;

    /** @brief Return the requested element or stored value. */
    virtual FloatMPDPBall _get(MetricTag, MultiplePrecision, DoublePrecision) const = 0;
    /** @brief Return the requested element or stored value. */
    virtual FloatMPBall _get(MetricTag, MultiplePrecision, MultiplePrecision) const = 0;
    /** @brief Return the requested element or stored value. */
    virtual FloatMPBounds _get(OrderTag, MultiplePrecision) const = 0;
    /** @brief Return the requested element or stored value. */
    virtual FloatMPUpperBound _get(UpperTag, MultiplePrecision) const = 0;
    /** @brief Return the requested element or stored value. */
    virtual FloatMPLowerBound _get(LowerTag, MultiplePrecision) const = 0;
    /** @brief Return the requested element or stored value. */
    virtual FloatMPApproximation _get(ApproximateTag, MultiplePrecision) const = 0;

    // FIXME: Only used in templates; should not really be needed
    /** @brief Return the requested element or stored value. */
    FloatDPBall _get(MetricTag p, DoublePrecision pr) const;
    /** @brief Return the requested element or stored value. */
    FloatMPBall _get(MetricTag p, MultiplePrecision pr) const;

    /** @brief Return the requested element or stored value. */
    FloatDPBounds _get(ValidatedTag, DoublePrecision pr) const;
    /** @brief Return the requested element or stored value. */
    FloatMPBounds _get(ValidatedTag, MultiplePrecision pr) const;

    /** @brief Test whether the represented number is positive. */
    virtual LogicalInterface* _is_pos() const = 0;

    /** @brief Return the computational paradigm of the represented value. */
    virtual ParadigmCode _paradigm() const = 0;
    /** @brief Return the display name of the type. */
    virtual String _class_name() const = 0;
};

} // namespace Ariadne

#endif /* ARIADNE_NUMBER_INTERFACE */
