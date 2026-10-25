/***************************************************************************
 *            numeric/real_interface.hpp
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

#ifndef ARIADNE_REAL_INTERFACE_HPP
#define ARIADNE_REAL_INTERFACE_HPP

#include <iosfwd>
#include "utility/typedefs.hpp"

namespace Ariadne {

class Real;
class ValidatedReal;

class Dyadic;

class DoublePrecision;
class ArbitraryPrecision;
template<class... PRS> class Float;
using FloatDP=Float<DoublePrecision>;
using FloatAP=Float<ArbitraryPrecision>;

template<class F> class Bounds;
using DyadicBounds = Bounds<Dyadic>;
using FloatDPBounds = Bounds<FloatDP>;
using FloatAPBounds = Bounds<FloatAP>;


class Real;

class RealInterface {
  public:
    virtual ~RealInterface() = default;
    virtual ValidatedReal _compute(Accuracy) const = 0;
    virtual ValidatedReal _compute(Effort) const = 0;
    virtual DyadicBounds _compute_get(Effort) const = 0;
    virtual FloatDPBounds _compute_get(Effort, DoublePrecision) const = 0;
    virtual FloatAPBounds _compute_get(Effort, ArbitraryPrecision) const = 0;
  public:
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

class ValidatedRealInterface {
  public:
    virtual ~ValidatedRealInterface() = default;
    virtual DyadicBounds _get() const = 0;
    virtual FloatDPBounds _get(DoublePrecision) const = 0;
    virtual FloatAPBounds _get(ArbitraryPrecision) const = 0;
  public:
    virtual OutputStream& _write(OutputStream& os) const = 0;
};

} // namespace Ariadne

#endif /* ARIADNE_REAL_INTERFACE_HPP */
