/***************************************************************************
 *            numeric/float_error.tpl.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/** @file numeric/float_error.tpl.hpp
 *  @brief Template implementations for float error in Ariadne Numeric.
 */

#ifndef ARIADNE_FLOAT_ERROR_TPL_HPP
#define ARIADNE_FLOAT_ERROR_TPL_HPP

#include "float_bounds.hpp"
#include "float_error.hpp"
#include "upper_number.hpp"

namespace Ariadne {

template<class F>
Error<F>::Error(const PositiveBounds<F>& x)
    : _e(x._u) {}

template<class F>
Error<F>::Error(const Positive<F>& x)
    : _e(cast_unsigned(x)) {}

template<class F>
Error<F>& Error<F>::operator=(ValidatedErrorNumber y) {
    return *this = cast_positive(y.get(this->precision()));
}

template<class F>
auto Error<F>::generic() const -> ValidatedErrorNumber {
    return this->operator ValidatedErrorNumber();
}

template<class F>
struct Operations<Error<F>> {
    static OutputStream& _write(OutputStream& os, const Error<F>& x) {
        return write(os, x.raw(), DecimalPrecision{Error<F>::output_places}, upward);
    }
    static InputStream& _read(InputStream& is, Error<F>& x) {
        F e(x.precision());
        is >> e;
        if (!is) {
            return is;
        }
        if (e < 0) {
            ARIADNE_THROW(std::runtime_error, "operator>>(InputStream&,Error<F>&)", "negative error bound " << e);
        }
        x.raw() = e;
        return is;
    }
};

} // namespace Ariadne

#endif
