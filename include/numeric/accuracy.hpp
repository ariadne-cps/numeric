/***************************************************************************
 *            numeric/accuracy.hpp
 *
 *  Copyright  2019-20  Pieter Collins
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

/** @file numeric/accuracy.hpp
 *  @brief A class representing a bound on the accuracy of a computation in a metric space.
 */

#ifndef ARIADNE_ACCURACY_HPP
#define ARIADNE_ACCURACY_HPP

#include <iosfwd>

#include "bits.hpp"
#include "dyadic.hpp"
#include "utility/typedefs.hpp"

namespace Ariadne {

class Bits;

//! \ingroup NumericModule
//! \brief The accuracy of a computation of a value in a metric space.
class Accuracy {
    Dyadic _error;

  public:
    //! <p/>
    /** @brief Construct a \c Accuracy from precision. */
    Accuracy(Bits precision);
    //! <p/>
    /** @brief Construct a \c Accuracy from error. */
    Accuracy(Dyadic error) : _error(error) {}
    //! <p/>
    /** @brief Return the error bound. */
    const Dyadic& error() const { return this->_error; }
    //! <p/>
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const Accuracy& acc);
};
//! \relates Accuracy \brief <p/>
/** @brief Return the requested accuracy. */
inline Accuracy accuracy(Bits acc) { return Accuracy(acc); }

} // namespace Ariadne

#endif // ARIADNE_ACCURACY_HPP
