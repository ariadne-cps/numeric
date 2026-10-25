/***************************************************************************
 *            numeric/concepts.cpp
 *
 *  Copyright  2020  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
1 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "concepts.hpp"
#include "archetypes.hpp"

#include "builtin.hpp"
#include "dyadic.hpp"
#include "decimal.hpp"
#include "rational.hpp"
#include "real.hpp"

#include "floatdp.hpp"
#include "floatap.hpp"
#include "float_ball.hpp"
#include "float_bounds.hpp"
#include "float_upper_bound.hpp"
#include "float_lower_bound.hpp"
#include "float_approximation.hpp"



namespace Ariadne {

static_assert(Ring<Integer> && OrderedLattice<Integer>);
static_assert(DyadicRing<Dyadic> && OrderedLattice<Dyadic>);
static_assert(Field<Rational> && OrderedLattice<Rational>);
static_assert(TranscendentalField<Real> && OrderedLattice<Real>);


static_assert(IsRoundedLatticeField<FloatDP>);
static_assert(IsRoundedLatticeField<FloatAP>);
static_assert(IsRoundedLatticeField<RoundedArchetype>);

static_assert(OrderedLatticeTranscendentalField<Ball<FloatDP,FloatDP>>);
static_assert(OrderedLatticeTranscendentalField<Ball<FloatAP,FloatDP>>);
static_assert(OrderedLatticeTranscendentalField<Ball<FloatAP,FloatAP>>);

static_assert(OrderedLatticeTranscendentalField<Bounds<FloatDP>>);
static_assert(OrderedLatticeTranscendentalField<Bounds<FloatAP>>);

static_assert(OrderedLatticeTranscendentalField<Approximation<FloatDP>>);
static_assert(OrderedLatticeTranscendentalField<Approximation<FloatAP>>);

static_assert(DirectedGroups<LowerBound<FloatDP>,UpperBound<FloatDP>>);
static_assert(DirectedGroup<LowerBound<FloatDP>>);

static_assert(Lattice<LowerBound<FloatDP>>);
static_assert(DirectedGroup<LowerBound<FloatDP>>);
static_assert(Monotone<LowerBound<FloatDP>>);

static_assert(DirectedSemiField<Positive<UpperBound<FloatDP>>>);
static_assert(DirectedSemiField<Positive<UpperBound<FloatAP>>>);
static_assert(DirectedSemiField<Positive<LowerBound<FloatDP>>>);
static_assert(DirectedSemiField<Positive<LowerBound<FloatAP>>>);

} // namespace Ariadne
