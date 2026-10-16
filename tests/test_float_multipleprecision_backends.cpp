/***************************************************************************
 *            test_float_multipleprecision_backends.cpp
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

// Shared differential conformance test; MPFR is the oracle, not a runtime backend.
#include "numeric/floatmp.hpp"
#include "numeric/rational.hpp"
#ifdef ARIADNE_NUMERIC_ENABLE_FLOATAP
#include "numeric/floatap.hpp"
#endif
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace Ariadne;

template< class F >
void verify(char const* backend, unsigned bits, char const* op,
F const& actual, FloatMP const& expected) {
    // Convert exactly, not via double: the latter loses differences above 53 bits.
    const Rational obtained = static_cast<Rational>(actual);
    const Rational oracle = static_cast<Rational>(expected);
    if (!(obtained == oracle)) {
        std::ostringstream out;
        out << backend << " " << bits << " bits " << op
        << ": obtained " << obtained << ", MPFR " << oracle;
        throw std::runtime_error(out.str());
    }
}

template< class F >
void check_mode(char const* backend, unsigned bits,
typename F::RoundingModeType mode, FloatMP::RoundingModeType reference_mode,
typename F::PrecisionType precision,
Rational const& lhs, Rational const& rhs) {
    F x(lhs, mode, precision), y(rhs, mode, precision);
    const MultiplePrecision oracle_precision = bits == 53 ? MultiplePrecision(53_bits)
    : bits == 128 ? MultiplePrecision(128_bits)
    : bits == 256 ? MultiplePrecision(256_bits) : MultiplePrecision(1024_bits);
    FloatMP reference_x(lhs, reference_mode, oracle_precision);
    FloatMP reference_y(rhs, reference_mode, oracle_precision);
    verify(backend, bits, "convert lhs", x, reference_x);
    verify(backend, bits, "convert rhs", y, reference_y);
    verify(backend, bits, "add", add(mode, x, y), add(reference_mode, reference_x, reference_y));
    verify(backend, bits, "sub", sub(mode, x, y), sub(reference_mode, reference_x, reference_y));
    verify(backend, bits, "mul", mul(mode, x, y), mul(reference_mode, reference_x, reference_y));
    verify(backend, bits, "div", div(mode, x, y), div(reference_mode, reference_x, reference_y));
    verify(backend, bits, "fma", fma(mode, x, y, x), fma(reference_mode, reference_x, reference_y, reference_x));
#ifndef ARIADNE_NUMERIC_ENABLE_FLOATAP
    verify(backend, bits, "sqrt", sqrt(mode, x), sqrt(reference_mode, reference_x));
    verify(backend, bits, "exp", exp(mode, x), exp(reference_mode, reference_x));
    verify(backend, bits, "log", log(mode, x), log(reference_mode, reference_x));
    verify(backend, bits, "sin", sin(mode, x), sin(reference_mode, reference_x));
    verify(backend, bits, "cos", cos(mode, x), cos(reference_mode, reference_x));
    verify(backend, bits, "tan", tan(mode, x), tan(reference_mode, reference_x));
    verify(backend, bits, "asin", asin(mode, x), asin(reference_mode, reference_x));
    verify(backend, bits, "acos", acos(mode, x), acos(reference_mode, reference_x));
    verify(backend, bits, "atan", atan(mode, x), atan(reference_mode, reference_x));
    verify(backend, bits, "pow_7", pow(mode, x, Nat(7u)), pow(reference_mode, reference_x, Nat(7u)));
    verify(backend, bits, "pi", F::pi(mode, precision), FloatMP::pi(reference_mode, oracle_precision));
#endif
}

template< class F >
void run_backend(char const* backend) {
    for (unsigned bits : { 53u, 128u, 256u, 1024u }) {
        using Precision = typename F::PrecisionType;
        Precision precision = bits == 53 ? Precision(53_bits) : bits == 128 ? Precision(128_bits)
        : bits == 256 ? Precision(256_bits) : Precision(1024_bits);
        const Integer scale = pow(Integer(2), Nat(bits));
        // Both operands are in (0,1); transcendental functions remain in-domain.
        const Rational lhs(scale + Integer(13), 2* scale + Integer(17));
        const Rational rhs(scale + Integer(31), 3* scale + Integer(19));
        check_mode < F >(backend, bits, F::ROUND_TO_NEAREST, FloatMP::ROUND_TO_NEAREST, precision, lhs, rhs);
        check_mode < F >(backend, bits, F::ROUND_DOWNWARD, FloatMP::ROUND_DOWNWARD, precision, lhs, rhs);
        check_mode < F >(backend, bits, F::ROUND_UPWARD, FloatMP::ROUND_UPWARD, precision, lhs, rhs);
        check_mode < F >(backend, bits, F::ROUND_TOWARD_ZERO, FloatMP::ROUND_TOWARD_ZERO, precision, lhs, rhs);

        // Independent property: directed conversion encloses the exact rational.
        F low(lhs, F::ROUND_DOWNWARD, precision), high(lhs, F::ROUND_UPWARD, precision);
        if (!(static_cast<Rational>(low) <= lhs && lhs <= static_cast<Rational>(high))) {
            throw std::runtime_error(std::string(backend) + ": directed conversion failed to enclose input");
        }
        F positive_infinity = F::inf(Sign::POSITIVE, precision);
        F not_a_number = F::nan(precision);
        if (!is_inf(positive_infinity) || !is_nan(not_a_number)) {
            throw std::runtime_error(std::string(backend) + ": special-value classification failed");
        }
    }
}

int main() {
    try {
        run_backend < FloatMP >("mpfr");
#ifdef ARIADNE_NUMERIC_ENABLE_FLOATAP
        run_backend < FloatAP >("floatap");
#endif
        std::cout << "float multiprecision backend conformance OK\n";
        return 0;
    }
    catch (std::exception const& e) {
        std::cerr << "float multiprecision backend conformance FAILED: " << e.what() << '\n';
        return 1;
    }
}

