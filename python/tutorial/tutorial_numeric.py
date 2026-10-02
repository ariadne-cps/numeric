#!/usr/bin/python3

##############################################################################
#            tutorial_numeric.py
#
#  Copyright  2026  Luca Geretti
##############################################################################

# This file is part of Ariadne.
#
# Ariadne is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Ariadne is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Ariadne. If not, see <https://www.gnu.org/licenses/>.

from pyariadne_numeric import *


def show(label, value):
    print(f"{label}: {value}")


def section(title):
    print(f"\n=== {title} ===")


def tutorial_numeric():
    # 1. Exact numbers
    section("Exact numbers")

    integer = Integer(5)
    show("Integer 5", integer)

    dyadic = Dyadic(11, 3)
    show("Dyadic 11/8", dyadic)

    # Python uses ** for exponentiation.
    same_dyadic = 11 / (two ** 3)
    show("11/(two**3)", same_dyadic)

    decimal = Decimal("9.81")
    show("Decimal 9.81", decimal)

    rational = Rational(11, 8)
    show("Rational 11/8", rational)

    dyadic_from_integer = Dyadic(integer)
    rational_from_decimal = Rational(decimal)
    show("Dyadic converted from Integer", dyadic_from_integer)
    show("Rational converted from Decimal", rational_from_decimal)

    # 2. Exact binary floating-point values
    section("Exact binary floating-point values")

    exact_double = ExactDouble(1.375)
    show("ExactDouble 1.375", exact_double)

    # Python has no C++ user-defined literals; x_(...) is the binding of _x.
    exact_double = x_(1.375)
    show("x_(1.375)", exact_double)

    dyadic_from_double = Dyadic(exact_double)
    show("Dyadic converted from x_(1.375)", dyadic_from_double)

    # 3. Computable real numbers
    section("Computable real numbers")

    one = Real(1)
    two_real = Real(2)
    one_third = Real(Rational(1, 3))
    show("Real 1/3", one_third)

    arithmetic_expression = sqr(one_third) + one_third - one_third / two_real
    show("Arithmetic expression", arithmetic_expression)

    sqrt_two = sqrt(two_real)
    exponential = exp(one)
    logarithm = log(two_real)
    sine = sin(one_third)
    cosine = cos(one_third)
    arctangent = atan(one)

    show("sqrt(2)", sqrt_two)
    show("exp(1)", exponential)
    show("log(2)", logarithm)
    show("sin(1/3)", sine)
    show("cos(1/3)", cosine)
    show("atan(1)", arctangent)

    pi_value = 4 * atan(Real(1))
    show("pi as a Real", pi_value)

    # 4. Precision and concrete representations
    section("Floating-point precision")

    dp = DoublePrecision()
    mp = MultiplePrecision(128)

    pi_dp = pi_value.get(dp)
    pi_mp = pi_value.get(mp)

    show("pi with double precision", pi_dp)
    show("pi with 128-bit multiple precision", pi_mp)
    show("double-precision lower bound", pi_dp.lower())
    show("double-precision upper bound", pi_dp.upper())
    show("double-precision error", pi_dp.error())
    show("multiple-precision error", pi_mp.error())

    # 5. Accuracy is different from precision
    section("Requested accuracy")

    requested_accuracy = Accuracy(bits=128)
    accurate_pi = pi_value.compute(requested_accuracy)
    show("pi computed to 128 bits of accuracy", accurate_pi)

    accurate_pi_mp = accurate_pi.get(mp)
    show("validated pi at 128-bit precision", accurate_pi_mp)
    show("validated pi error", accurate_pi_mp.error())

    # 6. Effort-controlled computation
    section("Effort-controlled computation")

    pi_with_small_effort = pi_value.compute(Effort(8))
    pi_with_more_effort = pi_value.compute(Effort(64))

    show("pi with Effort(8)", pi_with_small_effort)
    show("pi with Effort(64)", pi_with_more_effort)

    # 7. Raw floating-point numbers
    section("Raw floating-point numbers")

    raw_dp = FloatDP(x_(1.75), dp)
    raw_mp = FloatMP(x_(1.75), mp)

    show("FloatDP 1.75", raw_dp)
    show("FloatMP 1.75", raw_mp)

    # 8. Validated floating-point bounds
    section("Validated floating-point bounds")

    decimal_bounds = FloatDPBounds(Decimal("1.2"), dp)
    show("Bounds enclosing decimal 1.2", decimal_bounds)

    interval = FloatDPBounds(Rational(11, 10), Rational(14, 10), dp)
    show("Bounds enclosing [11/10,14/10]", interval)

    interval_squared = sqr(interval)
    interval_transformed = exp(interval) + interval_squared

    show("sqr(interval)", interval_squared)
    show("exp(interval) + sqr(interval)", interval_transformed)
    show("interval value", interval.value())
    show("interval error", interval.error())

    # 9. Multiple-precision validated bounds
    section("Multiple-precision bounds")

    mp_decimal_bounds = FloatMPBounds(Decimal("1.2"), mp)
    mp_interval = FloatMPBounds(Rational(11, 10), Rational(14, 10), mp)

    show("128-bit bounds for decimal 1.2", mp_decimal_bounds)
    show("128-bit bounds for [11/10,14/10]", mp_interval)
    show("128-bit decimal bound error", mp_decimal_bounds.error())

    # 10. Approximate floating-point numbers
    section("Floating-point approximations")

    approximation_dp = FloatDPApproximation(1.23, dp)
    approximation_mp = FloatMPApproximation(1.23, mp)

    show("double-precision approximation to 1.23", approximation_dp)
    show("multiple-precision approximation to 1.23", approximation_mp)

    validated_123 = FloatDPBounds(Decimal("1.23"), dp)
    show("validated enclosure of exact decimal 1.23", validated_123)

    # 11. Ball representation
    section("Validated balls")

    pi_ball = FloatMPBall(pi_value.compute(Accuracy(bits=128)).get(mp))
    show("pi represented as a multiple-precision ball", pi_ball)
    show("ball error", pi_ball.error())

    # 12. Generic numeric paradigms
    section("Generic numeric paradigms")

    exact_number = ExactNumber(Rational(1, 3))
    effective_number = EffectiveNumber(pi_value)
    validated_number = ValidatedNumber(pi_dp)
    approximate_number = ApproximateNumber(approximation_dp)

    show("ExactNumber", exact_number)
    show("EffectiveNumber", effective_number)
    show("ValidatedNumber", validated_number)
    show("ApproximateNumber", approximate_number)

    show("ExactNumber dynamic representation", exact_number.class_name())
    show("ValidatedNumber dynamic representation", validated_number.class_name())

    # 13. Comparisons of computable real numbers
    section("Comparisons with finite effort")

    exact_third = Real(Rational(1, 3))
    binary_third = Real(ExactDouble(0.333333333333333333))
    comparison = exact_third > binary_third

    show("1/3 > exact binary double, Effort(0)", comparison.check(Effort(0)))
    show("1/3 > exact binary double, Effort(64)", comparison.check(Effort(64)))

    # 14. A complete rigorous computation
    section("Complete rigorous computation")

    x = sqrt(Real(2))
    result = exp(-x) * sin(pi_value / 4)

    validated_result = result.compute(Accuracy(bits=128))
    result_bounds = validated_result.get(mp)

    show("result as a Real", result)
    show("result with requested accuracy", validated_result)
    show("128-bit rigorous bounds", result_bounds)
    show("rigorous error bound", result_bounds.error())


if __name__ == "__main__":
    tutorial_numeric()
