/***************************************************************************
 *            tutorial_numeric.cpp
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

#include <iostream>

#include "ariadne-numeric.hpp"

using namespace Ariadne;


/*
 * Print a labelled value.
 *
 * Ariadne numeric types provide stream output operators, so the same helper
 * can be used for exact numbers, floating-point values, validated bounds and
 * computable real numbers.
 */
template<class T> void print(const char* label, T const& value) {
    std::cout << label << ": " << value << '\n';
}


/*
 * Print a visual separator between the different parts of the tutorial.
 */
void section(const char* title) {
    std::cout << "\n=== " << title << " ===\n";
}


void tutorial_numeric() {
    //! [Numeric tutorial]

    /*
     * ----------------------------------------------------------------------
     * 1. Exact numbers
     * ----------------------------------------------------------------------
     *
     * Ariadne distinguishes several exact representations of numbers.
     *
     * Integer  represents arbitrary-precision integers.
     * Dyadic   represents rational numbers whose denominator is a power of 2.
     * Decimal  represents numbers given by a finite decimal expansion.
     * Rational represents arbitrary rational numbers.
     *
     * These classes are exact: constructing and manipulating them does not
     * introduce floating-point rounding errors.
     */

    section("Exact numbers");

    Integer integer(5);
    print("Integer 5", integer);

    /*
     * A dyadic number has the form p / 2^q.
     *
     * Here 11 / 2^3 = 11 / 8 exactly. The second argument is unsigned in
     * order to make clear that it is the exponent of the denominator.
     */
    Dyadic dyadic(11,3u);
    print("Dyadic 11/8", dyadic);

    /*
     * The same value can be written using Ariadne's 'two' object.
     *
     * Notice that ^ is not exponentiation in ordinary C++; Ariadne overloads
     * it for this numeric object.
     */
    Dyadic same_dyadic = 11/(two^3u);
    print("11/(two^3)", same_dyadic);

    /*
     * Decimal should be preferred when the source value is conceptually a
     * decimal quantity. The string is parsed exactly as the decimal number
     * 9.81; it is not first converted to a binary floating-point value.
     */
    Decimal decimal("9.81");
    print("Decimal 9.81", decimal);

    /*
     * Rational stores an exact quotient of arbitrary-precision integers.
     */
    Rational rational(11,8);
    print("Rational 11/8", rational);

    /*
     * Exact numeric types can be converted into one another whenever the
     * destination representation can represent the value.
     */
    Dyadic dyadic_from_integer(integer);
    Rational rational_from_decimal(decimal);
    print("Dyadic converted from Integer", dyadic_from_integer);
    print("Rational converted from Decimal", rational_from_decimal);


    /*
     * ----------------------------------------------------------------------
     * 2. Exact binary floating-point values
     * ----------------------------------------------------------------------
     *
     * A C++ double is normally understood as an approximation to some intended
     * mathematical value. Ariadne also provides ExactDouble for the different
     * situation in which the binary floating-point value itself is the exact
     * mathematical value.
     *
     * The _x literal is a convenient way to express this intention.
     */

    section("Exact binary floating-point values");

    ExactDouble exact_double(1.375);
    print("ExactDouble 1.375", exact_double);

    exact_double = 1.375_x;
    print("1.375_x", exact_double);

    /*
     * 1.375 is exactly representable in binary, so it can be converted exactly
     * to a Dyadic number.
     */
    Dyadic dyadic_from_double(exact_double);
    print("Dyadic converted from 1.375_x", dyadic_from_double);


    /*
     * ----------------------------------------------------------------------
     * 3. Computable real numbers
     * ----------------------------------------------------------------------
     *
     * Real represents an effectively computable real number rather than a
     * fixed floating-point approximation.
     *
     * A Real can therefore represent values such as sqrt(2), exp(1) or pi,
     * even though those values cannot be stored exactly in a finite ordinary
     * floating-point representation.
     */

    section("Computable real numbers");

    Real one(1);
    Real two_real(2);
    Real one_third(Rational(1,3));
    print("Real 1/3", one_third);

    /*
     * Arithmetic on Real values constructs another computable real number.
     */
    Real arithmetic_expression = sqr(one_third) + one_third - one_third/two_real;
    print("Arithmetic expression", arithmetic_expression);

    /*
     * Standard algebraic and transcendental operations are also available.
     */
    Real sqrt_two = sqrt(two_real);
    Real exponential = exp(one);
    Real logarithm = log(two_real);
    Real sine = sin(one_third);
    Real cosine = cos(one_third);
    Real arctangent = atan(one);

    print("sqrt(2)", sqrt_two);
    print("exp(1)", exponential);
    print("log(2)", logarithm);
    print("sin(1/3)", sine);
    print("cos(1/3)", cosine);
    print("atan(1)", arctangent);

    /*
     * A familiar example is pi = 4 atan(1).
     *
     * At this stage pi_value is still a Real. We have described the number,
     * but have not yet chosen a concrete floating-point representation for it.
     */
    Real pi_value = 4*atan(Real(1));
    print("pi as a Real", pi_value);


    /*
     * ----------------------------------------------------------------------
     * 4. Precision and concrete representations
     * ----------------------------------------------------------------------
     *
     * Precision describes the floating-point format used for a concrete
     * representation.
     *
     * DoublePrecision selects the ordinary double-precision backend.
     * MultiplePrecision allows the number of binary precision bits to be
     * chosen explicitly.
     */

    section("Floating-point precision");

    DoublePrecision dp;
    MultiplePrecision mp(128);

    /*
     * Asking a Real for a concrete representation produces validated bounds
     * enclosing the exact real value.
     *
     * The exact same Real can therefore be evaluated using different
     * floating-point precisions.
     */
    FloatDPBounds pi_dp = pi_value.get(dp);
    FloatMPBounds pi_mp = pi_value.get(mp);

    print("pi with double precision", pi_dp);
    print("pi with 128-bit multiple precision", pi_mp);

    /*
     * Bounds expose their lower and upper endpoints explicitly.
     */
    print("double-precision lower bound", pi_dp.lower());
    print("double-precision upper bound", pi_dp.upper());

    /*
     * The error gives a rigorous bound around the representative value.
     */
    print("double-precision error", pi_dp.error());
    print("multiple-precision error", pi_mp.error());


    /*
     * ----------------------------------------------------------------------
     * 5. Accuracy is different from precision
     * ----------------------------------------------------------------------
     *
     * Precision concerns the floating-point representation being used.
     *
     * Accuracy instead specifies how close the computed result must be to the
     * exact mathematical value.
     *
     * These are related concepts, but they are not interchangeable.
     */

    section("Requested accuracy");

    Accuracy requested_accuracy(128_bits);

    /*
     * compute(Accuracy(...)) asks Ariadne to refine the computation until the
     * requested error requirement is achieved.
     */
    ValidatedReal accurate_pi = pi_value.compute(requested_accuracy);
    print("pi computed to 128 bits of accuracy", accurate_pi);

    /*
     * Once computed, the validated result can again be materialised using a
     * chosen floating-point precision.
     */
    FloatMPBounds accurate_pi_mp = accurate_pi.get(mp);
    print("validated pi at 128-bit precision", accurate_pi_mp);
    print("validated pi error", accurate_pi_mp.error());


    /*
     * ----------------------------------------------------------------------
     * 6. Effort-controlled computation
     * ----------------------------------------------------------------------
     *
     * Sometimes the desired stopping criterion is computational effort rather
     * than a prescribed metric accuracy.
     *
     * Increasing Effort asks the implementation to work harder to improve the
     * result. Unlike Accuracy, Effort is not itself an error tolerance.
     */

    section("Effort-controlled computation");

    ValidatedReal pi_with_small_effort = pi_value.compute(Effort(8));
    ValidatedReal pi_with_more_effort = pi_value.compute(Effort(64));

    print("pi with Effort(8)", pi_with_small_effort);
    print("pi with Effort(64)", pi_with_more_effort);


    /*
     * ----------------------------------------------------------------------
     * 7. Raw floating-point numbers
     * ----------------------------------------------------------------------
     *
     * FloatDP and FloatMP are concrete floating-point values.
     *
     * Unlike Bounds, they do not by themselves describe an interval known to
     * contain an exact mathematical value. They represent a floating-point
     * value at the selected precision.
     */

    section("Raw floating-point numbers");

    FloatDP raw_dp(1.75_x,dp);
    FloatMP raw_mp(1.75_x,mp);

    print("FloatDP 1.75", raw_dp);
    print("FloatMP 1.75", raw_mp);


    /*
     * ----------------------------------------------------------------------
     * 8. Validated floating-point bounds
     * ----------------------------------------------------------------------
     *
     * Bounds are one of the central tools for rigorous numerics.
     *
     * FloatDPBounds and FloatMPBounds represent an unknown mathematical value
     * together with guaranteed lower and upper floating-point bounds.
     *
     * Arithmetic is outwardly rounded so that the resulting interval still
     * encloses every possible exact result.
     */

    section("Validated floating-point bounds");

    /*
     * Decimal 1.2 is not exactly representable in binary floating point.
     *
     * Constructing bounds from Decimal therefore produces two adjacent binary
     * values enclosing the exact decimal value 1.2.
     */
    FloatDPBounds decimal_bounds(Decimal("1.2"),dp);
    print("Bounds enclosing decimal 1.2", decimal_bounds);

    /*
     * Bounds can also represent an explicit interval of possible values.
     */
    FloatDPBounds interval(Rational(11,10),Rational(14,10),dp);
    print("Bounds enclosing [11/10,14/10]", interval);

    /*
     * The enclosure property is preserved by arithmetic.
     */
    FloatDPBounds interval_squared = sqr(interval);
    FloatDPBounds interval_transformed = exp(interval) + interval_squared;

    print("sqr(interval)", interval_squared);
    print("exp(interval) + sqr(interval)", interval_transformed);

    /*
     * The midpoint-like representative value and its rigorous error estimate
     * are both available.
     */
    print("interval value", interval.value());
    print("interval error", interval.error());


    /*
     * ----------------------------------------------------------------------
     * 9. Multiple-precision validated bounds
     * ----------------------------------------------------------------------
     *
     * The same validated semantics are available with arbitrary floating-point
     * precision.
     */

    section("Multiple-precision bounds");

    FloatMPBounds mp_decimal_bounds(Decimal("1.2"),mp);
    FloatMPBounds mp_interval(Rational(11,10),Rational(14,10),mp);

    print("128-bit bounds for decimal 1.2", mp_decimal_bounds);
    print("128-bit bounds for [11/10,14/10]", mp_interval);

    /*
     * Increasing floating-point precision normally produces a much tighter
     * enclosure of an exactly specified value.
     */
    print("128-bit decimal bound error", mp_decimal_bounds.error());


    /*
     * ----------------------------------------------------------------------
     * 10. Approximate floating-point numbers
     * ----------------------------------------------------------------------
     *
     * Approximation expresses different semantics from Bounds.
     *
     * A FloatDPApproximation is intended for ordinary approximate numerical
     * computation. It stores an approximation, but does not provide the same
     * enclosure guarantee as FloatDPBounds.
     */

    section("Floating-point approximations");

    FloatDPApproximation approximation_dp(1.23,dp);
    FloatMPApproximation approximation_mp(1.23,mp);

    print("double-precision approximation to 1.23", approximation_dp);
    print("multiple-precision approximation to 1.23", approximation_mp);

    /*
     * Here the C++ literal 1.23 is first interpreted as a built-in floating
     * point value. This is appropriate for approximate computation, but it is
     * deliberately different from constructing validated bounds from
     * Decimal("1.23").
     */
    FloatDPBounds validated_123(Decimal("1.23"),dp);
    print("validated enclosure of exact decimal 1.23", validated_123);


    /*
     * ----------------------------------------------------------------------
     * 11. Ball representation
     * ----------------------------------------------------------------------
     *
     * A validated number can also be represented by a centre value together
     * with an explicit error radius.
     *
     * Bounds and balls encode rigorous numerical information in different
     * forms and are useful in different algorithms.
     */

    section("Validated balls");

    FloatMPBall pi_ball(pi_value.compute(Accuracy(128_bits)).get(mp));
    print("pi represented as a multiple-precision ball", pi_ball);
    print("ball error", pi_ball.error());


    /*
     * ----------------------------------------------------------------------
     * 12. Generic numeric paradigms
     * ----------------------------------------------------------------------
     *
     * Ariadne also provides generic Number wrappers describing the kind of
     * information available about a value.
     *
     * ExactNumber       : an exact representation is available.
     * EffectiveNumber   : the value can be computed to increasing accuracy.
     * ValidatedNumber   : rigorous finite information enclosing the value.
     * ApproximateNumber : approximate numerical information.
     *
     * These generic types are useful when an algorithm cares about the
     * semantics of the information rather than the concrete numeric class.
     */

    section("Generic numeric paradigms");

    ExactNumber exact_number(Rational(1,3));
    EffectiveNumber effective_number(pi_value);
    ValidatedNumber validated_number(pi_dp);
    ApproximateNumber approximate_number(approximation_dp);

    print("ExactNumber", exact_number);
    print("EffectiveNumber", effective_number);
    print("ValidatedNumber", validated_number);
    print("ApproximateNumber", approximate_number);

    print("ExactNumber dynamic representation", exact_number.class_name());
    print("ValidatedNumber dynamic representation", validated_number.class_name());


    /*
     * ----------------------------------------------------------------------
     * 13. Comparisons of computable real numbers
     * ----------------------------------------------------------------------
     *
     * Comparison deserves special attention.
     *
     * For general computable real numbers, finite computation may not always
     * be sufficient to decide whether two values are ordered. Ariadne
     * therefore does not reduce every comparison immediately to an ordinary
     * C++ bool.
     *
     * A comparison can be checked using a finite amount of Effort. More
     * effort may refine an initially indeterminate result.
     */

    section("Comparisons with finite effort");

    Real exact_third(Rational(1,3));

    /*
     * The C++ literal below is first rounded to the nearest binary double, and
     * ExactDouble then treats that binary value itself as exact.
     *
     * Consequently exact_third and binary_third represent two different exact
     * real numbers, although they are extremely close.
     */
    Real binary_third(ExactDouble(0.333333333333333333));
    auto comparison = exact_third > binary_third;

    print("1/3 > exact binary double, Effort(0)", comparison.check(Effort(0)));
    print("1/3 > exact binary double, Effort(64)", comparison.check(Effort(64)));


    /*
     * ----------------------------------------------------------------------
     * 14. A complete rigorous computation
     * ----------------------------------------------------------------------
     *
     * We finish by putting the main ideas together.
     *
     *  - Define the mathematical quantity using Real.
     *  - Perform algebraic and transcendental operations symbolically at the
     *    level of computable real numbers.
     *  - Request a rigorous accuracy.
     *  - Materialise the result using a chosen floating-point precision.
     */

    section("Complete rigorous computation");

    Real x = sqrt(Real(2));
    Real result = exp(-x) * sin(pi_value/4);

    ValidatedReal validated_result = result.compute(Accuracy(128_bits));
    FloatMPBounds result_bounds = validated_result.get(mp);

    print("result as a Real", result);
    print("result with requested accuracy", validated_result);
    print("128-bit rigorous bounds", result_bounds);
    print("rigorous error bound", result_bounds.error());

    //! [Numeric tutorial]
}


int main() {
    tutorial_numeric();
    return 0;
}
