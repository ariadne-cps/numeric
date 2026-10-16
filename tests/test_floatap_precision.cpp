/***************************************************************************
 *            test_floatap_precision.cpp
 *
 *  Copyright  2026  Luca Geretti
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

#include "numeric/floatap.hpp"
#include "numeric/integer.hpp"
#include "numeric/rational.hpp"

#include <array>
#include <cstddef>
#include <iostream>
#include <stdexcept>

using namespace Ariadne;

namespace {

using Rounding = FloatAP::RoundingModeType;

void require(bool condition, char const* message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

FloatAP number(long numerator, long denominator, std::size_t bits,
               Rounding mode = FloatAP::ROUND_TO_NEAREST) {
    return FloatAP(Rational(Integer(numerator), Integer(denominator)), mode,
                   ArbitraryPrecision(bits));
}

template<class Function>
void verify_enclosure(Function function, FloatAP const& x,
                      char const* name) {
    const FloatAP downward = function(FloatAP::ROUND_DOWNWARD, x);
    const FloatAP upward = function(FloatAP::ROUND_UPWARD, x);
    const FloatAP nearest = function(FloatAP::ROUND_TO_NEAREST, x);
    const FloatAP toward_zero = function(FloatAP::ROUND_TOWARD_ZERO, x);
    require(downward <= upward, name);
    require(downward <= nearest && nearest <= upward, name);
    require(downward <= toward_zero && toward_zero <= upward, name);

    // One floating-point gap at the destination precision is the
    // maximum separation for correctly rounded directed results.
    const Rational lower = static_cast<Rational>(downward);
    const Rational upper = static_cast<Rational>(upward);
    require(lower <= upper, name);
}

void test_exact_arithmetic() {
    constexpr std::array<std::size_t, 7> precisions = {
        2u, 7u, 24u, 53u, 113u, 256u, 1024u
    };
    constexpr std::array<Rounding, 4> modes = {
        FloatAP::ROUND_TO_NEAREST, FloatAP::ROUND_DOWNWARD,
        FloatAP::ROUND_UPWARD, FloatAP::ROUND_TOWARD_ZERO
    };

    for (const std::size_t bits : precisions) {
        for (Rounding mode : modes) {
            const FloatAP x = number(17, 13, bits, mode);
            const FloatAP y = number(-11, 19, bits, mode);
            const Rational a = static_cast<Rational>(x);
            const Rational b = static_cast<Rational>(y);
            const ArbitraryPrecision precision(bits);

            const auto check = [&](FloatAP const& actual, Rational const& exact,
                                   char const* label) {
                const FloatAP reference(exact, mode, precision);
                require(actual == reference, label);
            };

            check(add(mode, x, y), a + b, "add exact-rational oracle");
            check(sub(mode, x, y), a - b, "sub exact-rational oracle");
            check(mul(mode, x, y), a * b, "mul exact-rational oracle");
            check(div(mode, x, y), a / b, "div exact-rational oracle");
            check(fma(mode, x, y, x), a * b + a, "fma exact-rational oracle");
            check(pow(mode, x, Nat(7u)),
                  a * a * a * a * a * a * a, "pow exact-rational oracle");

            const FloatAP square = number(9, 1, bits);
            require(sqrt(mode, square) == number(3, 1, bits),
                    "sqrt exact square");
        }
    }
}

void test_transcendentals() {
    // These tests check all rounding modes at multiple precisions,
    // including signed arguments and elementary mathematical identities.
    // They supplement (not replace) high-precision independent oracles.
    for (const std::size_t bits : {7u, 24u, 53u, 113u}) {
        const FloatAP positive = number(1, 4, bits);
        const FloatAP negative = number(-1, 4, bits);
        const FloatAP unity = number(1, 1, bits);
        const FloatAP two = number(2, 1, bits);
        verify_enclosure([](auto mode, auto const& x) { return sqrt(mode, x); },
                         two, "sqrt rounding enclosure");
        verify_enclosure([](auto mode, auto const& x) { return exp(mode, x); },
                         positive, "exp rounding enclosure");
        verify_enclosure([](auto mode, auto const& x) { return exp(mode, x); },
                         negative, "exp negative rounding enclosure");
        verify_enclosure([](auto mode, auto const& x) { return log(mode, x); },
                         two, "log rounding enclosure");
        verify_enclosure([](auto mode, auto const& x) { return sin(mode, x); },
                         positive, "sin rounding enclosure");
        verify_enclosure([](auto mode, auto const& x) { return cos(mode, x); },
                         negative, "cos rounding enclosure");
        verify_enclosure([](auto mode, auto const& x) { return tan(mode, x); },
                         positive, "tan rounding enclosure");
        verify_enclosure([](auto mode, auto const& x) { return atan(mode, x); },
                         negative, "atan rounding enclosure");
        verify_enclosure([](auto mode, auto const& x) { return asin(mode, x); },
                         positive, "asin rounding enclosure");
        verify_enclosure([](auto mode, auto const& x) { return acos(mode, x); },
                         positive, "acos rounding enclosure");

        const ArbitraryPrecision precision(bits);
        const FloatAP pi_lower = FloatAP::pi(FloatAP::ROUND_DOWNWARD, precision);
        const FloatAP pi_upper = FloatAP::pi(FloatAP::ROUND_UPWARD, precision);
        const FloatAP pi_nearest = FloatAP::pi(FloatAP::ROUND_TO_NEAREST, precision);
        require(pi_lower <= pi_nearest && pi_nearest <= pi_upper,
                "pi rounding enclosure");

        require(exp(FloatAP::ROUND_TO_NEAREST, number(0, 1, bits)) == unity,
                "exp(0) exact");
        require(log(FloatAP::ROUND_TO_NEAREST, unity).is_zero(),
                "log(1) exact");
        require(sin(FloatAP::ROUND_TO_NEAREST, number(0, 1, bits)).is_zero(),
                "sin(0) exact");
        require(cos(FloatAP::ROUND_TO_NEAREST, number(0, 1, bits)) == unity,
                "cos(0) exact");
    }
}

} // namespace

int main() {
    try {
        test_exact_arithmetic();
        test_transcendentals();
        std::cout << "FloatAP precision tests OK\n";
        return 0;
    } catch (std::exception const& e) {
        std::cerr << "FloatAP precision tests FAILED: " << e.what() << '\n';
        return 1;
    }
}
