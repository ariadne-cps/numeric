/***************************************************************************
 *            test_rounding_mode.cpp
 *
 *  Copyright  2008-20  Pieter Collins, Luca Geretti
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

#include <iostream>
#include <sstream>

#include "numeric/rounding.hpp"
#include "numeric/floatap.hpp"
#include "numeric/rational.hpp"

#include "macro/testing.hpp"

namespace Ariadne { }

using namespace Ariadne;

class TestRoundingMode
{
  public:
    TestRoundingMode() { }
    void test() const;
};

void TestRoundingMode::test() const {
    #if defined ARIADNE_C99_ROUNDING
        std::cout << "Using standard fenv.hpp C header file for setting the rounding mode." << std::endl;
    #elif defined ARIADNE_BOOST_ROUNDING
        #if defined BOOST_NUMERIC_INTERVAL_DETAIL_C99_ROUNDING_CONTROL_HPP
            std::cout << "Using Boost interval library standard fenv.hpp C header for setting the rounding mode." << std::endl;
        #else
            std::cout << "Using Boost interval library hardware rounding for setting the rounding mode." << std::endl;
        #endif
    #elif defined ARIADNE_GCC_ROUNDING
        std::cout << "Using ordinary GCC inline assembler for setting the rounding mode." << std::endl;
    #elif defined ARIADNE_EGCC_ROUNDING
        std::cout << "Using extended GCC inline assembler for setting the rounding mode." << std::endl;
    #elif defined ARIADNE_SSE_ROUNDING
        std::cout << "Using SSE <xmmintrin.h> header file for setting the rounding mode." << std::endl;
    #elif defined ARIADNE_MSVC_ROUNDING
        std::cout << "Using Microsoft Visual Studio inline assembler for setting the rounding mode." << std::endl;
    #else
        std::cout << "Error: no rounding mode available." << std::endl;
        ++ARIADNE_TEST_FAILURES;
    #endif

    auto original_rounding=get_builtin_rounding_mode();

    set_builtin_rounding_to_nearest();
    ARIADNE_TEST_EQUALS(get_builtin_rounding_mode(),ROUND_TO_NEAREST);
    set_builtin_rounding_downward();
    ARIADNE_TEST_EQUALS(get_builtin_rounding_mode(),ROUND_DOWNWARD);
    set_builtin_rounding_upward();
    ARIADNE_TEST_EQUALS(get_builtin_rounding_mode(),ROUND_UPWARD);
    set_builtin_rounding_toward_zero();
    ARIADNE_TEST_EQUALS(get_builtin_rounding_mode(),ROUND_TOWARD_ZERO);
    set_builtin_rounding_mode(original_rounding);

    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(downward),ROUND_DOWNWARD);
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(to_nearest),ROUND_TO_NEAREST);
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(upward),ROUND_UPWARD);
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(toward_zero),ROUND_TOWARD_ZERO);
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(approximately),ROUND_TO_NEAREST);

    // Verify the same rounding tags map to the native arbitrary-precision
    // backend, using only the native rounding representation.
    ARIADNE_TEST_ASSERT(FloatAP::Rounding(downward)==FloatAP::ROUND_DOWNWARD);
    ARIADNE_TEST_ASSERT(FloatAP::Rounding(to_nearest)==FloatAP::ROUND_TO_NEAREST);
    ARIADNE_TEST_ASSERT(FloatAP::Rounding(upward)==FloatAP::ROUND_UPWARD);
    ARIADNE_TEST_ASSERT(FloatAP::Rounding(toward_zero)==FloatAP::ROUND_TOWARD_ZERO);
    ARIADNE_TEST_ASSERT(FloatAP::Rounding(approximately)==FloatAP::ROUND_TO_NEAREST);

    Rounding rnd_near(to_nearest);
    Rounding rnd_down(downward);
    Rounding rnd_up(upward);
    Rounding rnd_zero(toward_zero);
    Rounding rnd_approx(approximately);
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(rnd_near),ROUND_TO_NEAREST);
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(rnd_down),ROUND_DOWNWARD);
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(rnd_up),ROUND_UPWARD);
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(rnd_zero),ROUND_TOWARD_ZERO);
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(rnd_approx),ROUND_TO_NEAREST);

    const auto saved_native=FloatAP::get_rounding_mode();
    FloatAP::set_rounding_to_nearest();
    ARIADNE_TEST_ASSERT(FloatAP::get_rounding_mode()==FloatAP::ROUND_TO_NEAREST);
    FloatAP::set_rounding_downward();
    ARIADNE_TEST_ASSERT(FloatAP::get_rounding_mode()==FloatAP::ROUND_DOWNWARD);
    FloatAP::set_rounding_upward();
    ARIADNE_TEST_ASSERT(FloatAP::get_rounding_mode()==FloatAP::ROUND_UPWARD);
    FloatAP::set_rounding_toward_zero();
    ARIADNE_TEST_ASSERT(FloatAP::get_rounding_mode()==FloatAP::ROUND_TOWARD_ZERO);
    FloatAP::set_rounding_mode(saved_native);

    const ArbitraryPrecision native_precision(128_bits);
    const Rational third(1,3);
    const FloatAP lower(third,FloatAP::ROUND_DOWNWARD,native_precision);
    const FloatAP upper(third,FloatAP::ROUND_UPWARD,native_precision);
    const FloatAP toward_zero_positive(third,FloatAP::ROUND_TOWARD_ZERO,native_precision);
    const FloatAP toward_zero_negative(-third,FloatAP::ROUND_TOWARD_ZERO,native_precision);
    ARIADNE_TEST_ASSERT(Rational(lower)<third);
    ARIADNE_TEST_ASSERT(third<Rational(upper));
    ARIADNE_TEST_EQUALS(toward_zero_positive,lower);
    ARIADNE_TEST_EQUALS(toward_zero_negative,-lower);

    std::ostringstream stream;
    stream << rnd_near << " " << rnd_down << " " << rnd_up << " " << rnd_zero;
    ARIADNE_TEST_EQUALS(stream.str(),std::string("near down up zero"));
}


int main() {
   TestRoundingMode().test();
   return ARIADNE_TEST_FAILURES;
}

