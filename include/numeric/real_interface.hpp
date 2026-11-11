/***************************************************************************
 *            numeric/real_interface.hpp
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

/** @file numeric/real_interface.hpp
 *  @brief Declarations and utilities for real interface in Ariadne Numeric.
 */

#ifndef ARIADNE_REAL_INTERFACE_HPP
#define ARIADNE_REAL_INTERFACE_HPP

#include <iosfwd>

#include "utility/typedefs.hpp"

namespace Ariadne {

class Real;
class ValidatedReal;

class Dyadic;

class DoublePrecision;
class MultiplePrecision;
template<class... PRS>
class Float;
using FloatDP = Float<DoublePrecision>;
using FloatMP = Float<MultiplePrecision>;

template<class F>
class Bounds;
using DyadicBounds = Bounds<Dyadic>;
using FloatDPBounds = Bounds<FloatDP>;
using FloatMPBounds = Bounds<FloatMP>;

class Real;

/** @brief Polymorphic interface for real. */
class RealInterface {
  public:
    /** @brief Destroy the object and release its owned resources. */
    virtual ~RealInterface() = default;
    /** @brief Compute the value at the requested accuracy or effort. */
    virtual ValidatedReal _compute(Accuracy) const = 0;
    /** @brief Compute the value at the requested accuracy or effort. */
    virtual ValidatedReal _compute(Effort) const = 0;
    /** @brief Compute the value and extract the requested numeric representation. */
    virtual DyadicBounds _compute_get(Effort) const = 0;
    /** @brief Compute the value and extract the requested numeric representation. */
    virtual FloatDPBounds _compute_get(Effort, DoublePrecision) const = 0;
    /** @brief Compute the value and extract the requested numeric representation. */
    virtual FloatMPBounds _compute_get(Effort, MultiplePrecision) const = 0;

  public:
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

/** @brief Polymorphic interface for validated real. */
class ValidatedRealInterface {
  public:
    /** @brief Destroy the object and release its owned resources. */
    virtual ~ValidatedRealInterface() = default;
    /** @brief Return the requested element or stored value. */
    virtual DyadicBounds _get() const = 0;
    /** @brief Return the requested element or stored value. */
    virtual FloatDPBounds _get(DoublePrecision) const = 0;
    /** @brief Return the requested element or stored value. */
    virtual FloatMPBounds _get(MultiplePrecision) const = 0;

  public:
    /** @brief Write the value to the output stream. */
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

} // namespace Ariadne

#endif /* ARIADNE_REAL_INTERFACE_HPP */
