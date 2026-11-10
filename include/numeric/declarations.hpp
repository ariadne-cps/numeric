/***************************************************************************
 *            numeric/declarations.hpp
 *
 *  Copyright  2011-26  Pieter Collins
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


#ifndef ARIADNE_NUMERIC_DECLARATIONS_HPP
#define ARIADNE_NUMERIC_DECLARATIONS_HPP

#include <iosfwd>

#include "utility/metaprogramming.hpp"
#include "utility/typedefs.hpp"

#include "foundation/paradigm.hpp"
#include "foundation/logical.decl.hpp"

#include "numeric/number.decl.hpp"
#include "numeric/float.decl.hpp"

namespace Ariadne {


template<class X> struct InformationTypedef;
template<> struct InformationTypedef<Real> { typedef EffectiveTag Type; };
template<class P> struct InformationTypedef<Number<P>> { typedef P Type; };
template<class X> using InformationTag = typename InformationTypedef<X>::Type;

} // namespace Ariadne

#endif /* ARIADNE_NUMERIC_DECLARATIONS_HPP */
