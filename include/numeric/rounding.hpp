/***************************************************************************
 *            numeric/rounding.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! \file numeric/rounding.hpp
 *  \brief Functions to set and retrieve the processor rounding mode.
 *  May be platform-dependent.
 */

#ifndef ARIADNE_ROUNDING_HPP
#define ARIADNE_ROUNDING_HPP

#include <iosfwd>
#include "utility/typedefs.hpp"

#include <cfenv>

namespace Ariadne {

typedef int rounding_mode_t;

const rounding_mode_t ROUND_TO_NEAREST  = FE_TONEAREST;
const rounding_mode_t ROUND_DOWNWARD    = FE_DOWNWARD;
const rounding_mode_t ROUND_UPWARD      = FE_UPWARD;
const rounding_mode_t ROUND_TOWARD_ZERO = FE_TOWARDZERO;

inline void set_builtin_rounding_to_nearest() { std::fesetround(FE_TONEAREST); }
inline void set_builtin_rounding_downward() { std::fesetround(FE_DOWNWARD); }
inline void set_builtin_rounding_upward() { std::fesetround(FE_UPWARD); }
inline void set_builtin_rounding_toward_zero() { std::fesetround(FE_TOWARDZERO); }

inline void set_builtin_rounding_mode(rounding_mode_t rnd) { std::fesetround(rnd); }
inline rounding_mode_t get_builtin_rounding_mode() { return std::fegetround(); }

} // namespace Ariadne

/************  Import MPFR rounding mode controls *******************/

#include <mpfr.h>

/************  Publicly-accessible rounding-mode changing *******************/

namespace Ariadne {


//!@{
//! \ingroup NumericModule
//! \name Rounding mode control

//! \brief The rounding mode type used for builtin floating-point objects such as FloatDP. \ingroup NumericModule
typedef rounding_mode_t BuiltinRoundingModeType;
//! \brief The rounding mode type used for multiple-precision floating-point objects such as FloatMP. \ingroup NumericModule
typedef mpfr_rnd_t MPFRRoundingModeType;

//! \brief The floating-point environment value for rounding arithmetic to the nearest exactly-representable value.
extern const BuiltinRoundingModeType ROUND_TO_NEAREST;
//! \brief The floating-point environment value for upwards-rounded arithmetic.
extern const BuiltinRoundingModeType ROUND_DOWNWARD;
//! \brief The floating-point environment value for downwards-rounded arithmetic.
extern const BuiltinRoundingModeType ROUND_UPWARD;
//! \brief The floating-point environment value for rounding arithmetic to zero.
extern const BuiltinRoundingModeType ROUND_TOWARD_ZERO;

#ifdef DOXYGEN
//! \brief Set the builtin rounding mode. \ingroup NumericModule
void set_builtin_rounding_mode(BuiltinRoundingModeType rnd);
//! \brief Get the current rounding mode. \ingroup NumericModule
BuiltinRoundingModeType get_builtin_rounding_mode();

//! \brief Set the rounding mode to nearest. \ingroup NumericModule
void set_builtin_rounding_to_nearest();
//! \brief Set the rounding mode to downwards rounding. \ingroup NumericModule
void set_builtin_rounding_downward();
//! \brief Set the rounding mode to upwards rounding. \ingroup NumericModule
void set_builtin_rounding_upward();
//! \brief Set the rounding mode to towards-zero rounding. \ingroup NumericModule
void set_builtin_rounding_toward_zero();

//! \brief Set the rounding mode to the expected default rounding mode.
void set_default_builtin_rounding();
#endif

//! \brief The rounding mode type used for multiple-precision floating-point objects such as FloatMP. \ingroup NumericModule
typedef mpfr_rnd_t MPFRRoundingModeType;

//! \brief Tag class for downward rounding. Constants \ref downward, \ref down. \ingroup NumericModule
struct RoundDownward {
    constexpr operator BuiltinRoundingModeType() const { return ROUND_DOWNWARD; }
    constexpr operator MPFRRoundingModeType() const { return MPFR_RNDD; }
};
//! \brief Tag class for rounding to nearest. Constants \ref to_nearest, \ref near. \ingroup NumericModule
struct RoundToNearest {
    constexpr operator BuiltinRoundingModeType() const { return ROUND_TO_NEAREST; }
    constexpr operator MPFRRoundingModeType() const { return MPFR_RNDN; }
};
//! \brief Tag class for upward rounding. Constants \ref upward, \ref up. \ingroup NumericModule
struct RoundUpward {
    constexpr operator BuiltinRoundingModeType() const { return ROUND_UPWARD; }
    constexpr operator MPFRRoundingModeType() const { return MPFR_RNDU; }
};
//! \brief Tag class for rounding towards zero. %Constant \ref toward_zero. \ingroup NumericModule
struct RoundTowardZero {
    constexpr operator BuiltinRoundingModeType() const { return ROUND_TOWARD_ZERO; }
    constexpr operator MPFRRoundingModeType() const { return MPFR_RNDZ; }
};
//! \brief Tag class for approximate rounding. %Constant \ref approx. \ingroup NumericModule
struct RoundApproximately {
    constexpr operator BuiltinRoundingModeType() const { return ROUND_TO_NEAREST; }
    constexpr operator MPFRRoundingModeType() const { return MPFR_RNDN; }
};

struct CurrentRoundingMode { };
static const CurrentRoundingMode rounded = CurrentRoundingMode();

//! \brief General rounding mode class. \ingroup NumericModule
class Rounding {
    BuiltinRoundingModeType _rbp; MPFRRoundingModeType _rmp;
  public:
    Rounding(BuiltinRoundingModeType rbp, MPFRRoundingModeType rmp) : _rbp(rbp), _rmp(rmp) { } //!< <p/>
    Rounding(RoundDownward) : Rounding(ROUND_DOWNWARD,MPFR_RNDD) { } //!< <p/>
    Rounding(RoundToNearest) : Rounding(ROUND_TO_NEAREST,MPFR_RNDN) { } //!< <p/>
    Rounding(RoundUpward) :  Rounding(ROUND_UPWARD,MPFR_RNDU) { } //!< <p/>
    operator BuiltinRoundingModeType() const { return _rbp; } //!< <p/>
    operator MPFRRoundingModeType() const { return _rmp; } //!< <p/>
    friend OutputStream& operator<<(OutputStream& os, Rounding const& rnd) {
        return os << ( rnd._rbp == ROUND_TO_NEAREST ? "near" : (rnd._rbp == ROUND_DOWNWARD ? "down" : "up") ); } //!< <p/>
};


const RoundDownward downward = RoundDownward(); //!< Round exact answer downward to a representable value. Synonymous with \ref down. \ingroup NumericModule
const RoundToNearest to_nearest = RoundToNearest(); //!< Round exact answer to a nearest representable value. Synonymous with \ref near. \ingroup NumericModule
const RoundUpward upward = RoundUpward(); //!< Round exact answer upward to a representable value. Synonymous with \ref up. \ingroup NumericModule
const RoundTowardZero toward_zero = RoundTowardZero(); //!< Round exact answer to a representable value at least as close to zero. \ingroup NumericModule
const RoundApproximately approximately = RoundApproximately(); //!< Round exact answer to some close representable value, which need not be the nearest. Synonymous with \ref approx. \ingroup NumericModule
using RoundApprox = RoundApproximately; //!< . \ingroup NumericModule

const RoundDownward down = downward; //!< Round exact answer downward to a representable value. Synonymous with \ref downward. \ingroup NumericModule
const RoundToNearest near = to_nearest; //!< Round exact answer to a nearest representable value. Synonymous with \ref to_nearest. \ingroup NumericModule
const RoundUpward up = upward; //!< Round exact answer upward to a representable value. Synonymous with \ref upward. \ingroup NumericModule
const RoundApproximately approx = approximately; //!< Round exact answer to some close representable value, which need not be the nearest. Synonymous with \ref approximately. \ingroup NumericModule

//!@}

} // namespace Ariadne

#endif // ARIADNE_ROUNDING_HPP

