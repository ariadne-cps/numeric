/***************************************************************************
 *            numeric/declarations.hpp
 *
 *  Copyright  2011-26  Pieter Collins
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

/** @file numeric/declarations.hpp
 *  @brief Declarations and utilities for declarations in Ariadne Numeric.
 */

#ifndef ARIADNE_NUMERIC_DECLARATIONS_HPP
#define ARIADNE_NUMERIC_DECLARATIONS_HPP

#include <iosfwd>

#include "foundation/logical.decl.hpp"
#include "foundation/paradigm.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/number.decl.hpp"
#include "utility/metaprogramming.hpp"
#include "utility/typedefs.hpp"

namespace Ariadne {

template<class X>
struct InformationTypedef;
template<>
struct InformationTypedef<Real> {
    typedef EffectiveTag Type;
};
template<class P>
struct InformationTypedef<Number<P>> {
    typedef P Type;
};
template<class X>
using InformationTag = typename InformationTypedef<X>::Type;

} // namespace Ariadne

#endif /* ARIADNE_NUMERIC_DECLARATIONS_HPP */
