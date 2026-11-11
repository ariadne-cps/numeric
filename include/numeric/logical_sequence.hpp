/***************************************************************************
 *            numeric/logical_sequence.hpp
 *
 *  Copyright  2026  Luca Geretti
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

/** @file numeric/logical_sequence.hpp
 *  @brief Declarations and utilities for logical sequence in Ariadne Numeric.
 */

#ifndef ARIADNE_NUMERIC_LOGICAL_SEQUENCE_HPP
#define ARIADNE_NUMERIC_LOGICAL_SEQUENCE_HPP
#include "foundation/logical.hpp"
namespace Ariadne {
template<class X>
class Sequence;
/** @brief Return the logical disjunction. */
LowerKleenean disjunction(const Sequence<LowerKleenean>&);
/** @brief Return the logical conjunction. */
UpperKleenean conjunction(const Sequence<UpperKleenean>&);
} // namespace Ariadne
#endif
