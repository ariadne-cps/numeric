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
#include "numeric/floatmp.hpp"
#include "numeric/rational.hpp"

#include "ariadne-testing.hpp"

namespace Ariadne { }

using namespace Ariadne;

class TestRoundingMode : public Ariadne::Testing::ContextFixture
{
  public:
    using ContextFixture::ContextFixture;
    void test();
};

void TestRoundingMode::test() {
    #if defined ARIADNE_C99_ROUNDING
        ARIADNE_TESTING_NOTE("Using standard fenv.hpp C header file for setting the rounding mode.");
    #elif defined ARIADNE_BOOST_ROUNDING
        #if defined BOOST_NUMERIC_INTERVAL_DETAIL_C99_ROUNDING_CONTROL_HPP
            ARIADNE_TESTING_NOTE("Using Boost interval library standard fenv.hpp C header for setting the rounding mode.");
        #else
            ARIADNE_TESTING_NOTE("Using Boost interval library hardware rounding for setting the rounding mode.");
        #endif
    #elif defined ARIADNE_GCC_ROUNDING
        ARIADNE_TESTING_NOTE("Using ordinary GCC inline assembler for setting the rounding mode.");
    #elif defined ARIADNE_EGCC_ROUNDING
        ARIADNE_TESTING_NOTE("Using extended GCC inline assembler for setting the rounding mode.");
    #elif defined ARIADNE_SSE_ROUNDING
        ARIADNE_TESTING_NOTE("Using SSE <xmmintrin.h> header file for setting the rounding mode.");
    #elif defined ARIADNE_MSVC_ROUNDING
        ARIADNE_TESTING_NOTE("Using Microsoft Visual Studio inline assembler for setting the rounding mode.");
    #else
        ARIADNE_TESTING_NOTE("Error: no rounding mode available.");
        ARIADNE_TESTING_PROVE( false);
    #endif

    auto original_rounding=get_builtin_rounding_mode();

    set_builtin_rounding_to_nearest();
    ARIADNE_TESTING_EXPECT_EQUAL( get_builtin_rounding_mode(),ROUND_TO_NEAREST);
    set_builtin_rounding_downward();
    ARIADNE_TESTING_EXPECT_EQUAL( get_builtin_rounding_mode(),ROUND_DOWNWARD);
    set_builtin_rounding_upward();
    ARIADNE_TESTING_EXPECT_EQUAL( get_builtin_rounding_mode(),ROUND_UPWARD);
    set_builtin_rounding_toward_zero();
    ARIADNE_TESTING_EXPECT_EQUAL( get_builtin_rounding_mode(),ROUND_TOWARD_ZERO);

    // The generic setter must also apply each mode, not silently ignore it
    // on architectures without the x86/SSE control register.
    for (auto mode : {ROUND_TO_NEAREST, ROUND_DOWNWARD,
                      ROUND_UPWARD, ROUND_TOWARD_ZERO}) {
        set_builtin_rounding_mode(mode);
        ARIADNE_TESTING_EXPECT_EQUAL( get_builtin_rounding_mode(),mode);
    }
    set_builtin_rounding_mode(original_rounding);

    ARIADNE_TESTING_EXPECT_EQUAL( static_cast<BuiltinRoundingModeType>(downward),ROUND_DOWNWARD);
    ARIADNE_TESTING_EXPECT_EQUAL( static_cast<BuiltinRoundingModeType>(to_nearest),ROUND_TO_NEAREST);
    ARIADNE_TESTING_EXPECT_EQUAL( static_cast<BuiltinRoundingModeType>(upward),ROUND_UPWARD);
    ARIADNE_TESTING_EXPECT_EQUAL( static_cast<BuiltinRoundingModeType>(toward_zero),ROUND_TOWARD_ZERO);
    ARIADNE_TESTING_EXPECT_EQUAL( static_cast<BuiltinRoundingModeType>(approximately),ROUND_TO_NEAREST);

    // Verify the same rounding tags map to the native multiple-precision
    // backend, using only the native rounding representation.
    ARIADNE_TESTING_PROVE( FloatMP::Rounding(downward)==FloatMP::ROUND_DOWNWARD);
    ARIADNE_TESTING_PROVE( FloatMP::Rounding(to_nearest)==FloatMP::ROUND_TO_NEAREST);
    ARIADNE_TESTING_PROVE( FloatMP::Rounding(upward)==FloatMP::ROUND_UPWARD);
    ARIADNE_TESTING_PROVE( FloatMP::Rounding(toward_zero)==FloatMP::ROUND_TOWARD_ZERO);
    ARIADNE_TESTING_PROVE( FloatMP::Rounding(approximately)==FloatMP::ROUND_TO_NEAREST);

    Rounding rnd_builtin(ROUND_DOWNWARD);
    ARIADNE_TESTING_EXPECT_EQUAL( static_cast<BuiltinRoundingModeType>(rnd_builtin),ROUND_DOWNWARD);
    Rounding rnd_near(to_nearest);
    Rounding rnd_down(downward);
    Rounding rnd_up(upward);
    Rounding rnd_zero(toward_zero);
    Rounding rnd_approx(approximately);
    ARIADNE_TESTING_EXPECT_EQUAL( static_cast<BuiltinRoundingModeType>(rnd_near),ROUND_TO_NEAREST);
    ARIADNE_TESTING_EXPECT_EQUAL( static_cast<BuiltinRoundingModeType>(rnd_down),ROUND_DOWNWARD);
    ARIADNE_TESTING_EXPECT_EQUAL( static_cast<BuiltinRoundingModeType>(rnd_up),ROUND_UPWARD);
    ARIADNE_TESTING_EXPECT_EQUAL( static_cast<BuiltinRoundingModeType>(rnd_zero),ROUND_TOWARD_ZERO);
    ARIADNE_TESTING_EXPECT_EQUAL( static_cast<BuiltinRoundingModeType>(rnd_approx),ROUND_TO_NEAREST);

    const auto saved_native=FloatMP::get_rounding_mode();
    FloatMP::set_rounding_to_nearest();
    ARIADNE_TESTING_PROVE( FloatMP::get_rounding_mode()==FloatMP::ROUND_TO_NEAREST);
    FloatMP::set_rounding_downward();
    ARIADNE_TESTING_PROVE( FloatMP::get_rounding_mode()==FloatMP::ROUND_DOWNWARD);
    FloatMP::set_rounding_upward();
    ARIADNE_TESTING_PROVE( FloatMP::get_rounding_mode()==FloatMP::ROUND_UPWARD);
    FloatMP::set_rounding_toward_zero();
    ARIADNE_TESTING_PROVE( FloatMP::get_rounding_mode()==FloatMP::ROUND_TOWARD_ZERO);
    FloatMP::set_rounding_mode(saved_native);

    const MultiplePrecision native_precision(128_bits);
    const Rational third(1,3);
    const FloatMP lower(third,FloatMP::ROUND_DOWNWARD,native_precision);
    const FloatMP upper(third,FloatMP::ROUND_UPWARD,native_precision);
    const FloatMP toward_zero_positive(third,FloatMP::ROUND_TOWARD_ZERO,native_precision);
    const FloatMP toward_zero_negative(-third,FloatMP::ROUND_TOWARD_ZERO,native_precision);
    ARIADNE_TESTING_PROVE( Rational(lower)<third);
    ARIADNE_TESTING_PROVE( third<Rational(upper));
    ARIADNE_TESTING_EXPECT_EQUAL( toward_zero_positive,lower);
    ARIADNE_TESTING_EXPECT_EQUAL( toward_zero_negative,-lower);

    std::ostringstream stream;
    stream << rnd_near << " " << rnd_down << " " << rnd_up << " " << rnd_zero;
    ARIADNE_TESTING_EXPECT_EQUAL( stream.str(),std::string("near down up zero"));
}


int main(int argc, char** argv) {
    Ariadne::Testing::Suite suite("rounding_mode");
    suite.add("rounding mode contracts", [](Ariadne::Testing::Context& context) {
        TestRoundingMode(context).test();
    });
    return Ariadne::Testing::run_suite(suite, argc, argv);
}

