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

#include "utility/test.hpp"

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
    ARIADNE_TEST_EQUALS(static_cast<int>(static_cast<MPFRRoundingModeType>(downward)),static_cast<int>(MPFR_RNDD));
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(to_nearest),ROUND_TO_NEAREST);
    ARIADNE_TEST_EQUALS(static_cast<int>(static_cast<MPFRRoundingModeType>(to_nearest)),static_cast<int>(MPFR_RNDN));
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(upward),ROUND_UPWARD);
    ARIADNE_TEST_EQUALS(static_cast<int>(static_cast<MPFRRoundingModeType>(upward)),static_cast<int>(MPFR_RNDU));
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(toward_zero),ROUND_TOWARD_ZERO);
    ARIADNE_TEST_EQUALS(static_cast<int>(static_cast<MPFRRoundingModeType>(toward_zero)),static_cast<int>(MPFR_RNDZ));
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(approximately),ROUND_TO_NEAREST);
    ARIADNE_TEST_EQUALS(static_cast<int>(static_cast<MPFRRoundingModeType>(approximately)),static_cast<int>(MPFR_RNDN));

    Rounding rnd_near(to_nearest);
    Rounding rnd_down(downward);
    Rounding rnd_up(upward);
    Rounding rnd_zero(ROUND_TOWARD_ZERO,MPFR_RNDZ);
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(rnd_near),ROUND_TO_NEAREST);
    ARIADNE_TEST_EQUALS(static_cast<int>(static_cast<MPFRRoundingModeType>(rnd_down)),static_cast<int>(MPFR_RNDD));
    ARIADNE_TEST_EQUALS(static_cast<BuiltinRoundingModeType>(rnd_up),ROUND_UPWARD);
    ARIADNE_TEST_EQUALS(static_cast<int>(static_cast<MPFRRoundingModeType>(rnd_zero)),static_cast<int>(MPFR_RNDZ));

    std::ostringstream stream;
    stream << rnd_near << " " << rnd_down << " " << rnd_up;
    ARIADNE_TEST_EQUALS(stream.str(),std::string("near down up"));
}


int main() {
   TestRoundingMode().test();
   return ARIADNE_TEST_FAILURES;
}

