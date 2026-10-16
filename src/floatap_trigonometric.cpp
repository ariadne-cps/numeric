/***************************************************************************
 *            numeric/floatap_trigonometric.cpp
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

/*! \file numeric/floatap_trigonometric.cpp
 *  \brief Correctly rounded trigonometry using exact rational bounds.
 */

#include "numeric/floatap.hpp"
#include "numeric/rational.hpp"

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Ariadne {
namespace {

// Arithmetic on these fractions is exact. The approximation error is
// represented only by the endpoints of Interval, never by floating point.
class Fraction {
    mpq_t _value;

public:
    Fraction() {
        mpq_init(_value);
    }

    explicit Fraction(long integer) : Fraction() {
        mpq_set_si(_value, integer, 1);
    }

    explicit Fraction(mpq_srcptr rational) : Fraction() {
        mpq_set(_value, rational);
    }

    Fraction(mpz_srcptr numerator, unsigned long denominator_bits) : Fraction() {
        mpz_set(mpq_numref(_value), numerator);
        mpz_set_ui(mpq_denref(_value), 1);
        mpz_mul_2exp(mpq_denref(_value), mpq_denref(_value), denominator_bits);
        mpq_canonicalize(_value);
    }

    Fraction(Fraction const& other) : Fraction(other._value) { }

    Fraction(Fraction&& other) noexcept : Fraction() {
        mpq_swap(_value, other._value);
    }

    ~Fraction() {
        mpq_clear(_value);
    }

    Fraction& operator=(Fraction other) noexcept {
        mpq_swap(_value, other._value);
        return *this;
    }

    mpq_srcptr raw() const noexcept {
        return _value;
    }

    int sign() const noexcept {
        return mpq_sgn(_value);
    }

    friend int compare(Fraction const& a, Fraction const& b) {
        return mpq_cmp(a._value, b._value);
    }

    friend Fraction operator+(Fraction const& a, Fraction const& b) {
        Fraction result;
        mpq_add(result._value, a._value, b._value);
        return result;
    }

    friend Fraction operator-(Fraction const& a, Fraction const& b) {
        Fraction result;
        mpq_sub(result._value, a._value, b._value);
        return result;
    }

    friend Fraction operator*(Fraction const& a, Fraction const& b) {
        Fraction result;
        mpq_mul(result._value, a._value, b._value);
        return result;
    }

    friend Fraction operator/(Fraction const& a, Fraction const& b) {
        if (b.sign() == 0) {
            throw std::domain_error("FloatAP zero rational denominator");
        }
        Fraction result;
        mpq_div(result._value, a._value, b._value);
        return result;
    }
};

Fraction absolute(Fraction const& x) {
    return x.sign() < 0 ? Fraction{} - x : x;
}

struct Interval {
    Fraction lower;
    Fraction upper;
};

Interval point(Fraction const& value) {
    return {value, value};
}

Interval enclosure(Fraction const& value, Fraction const& error) {
    return {value - error, value + error};
}

Interval add(Interval const& a, Interval const& b) {
    return {a.lower + b.lower, a.upper + b.upper};
}

Interval subtract(Interval const& a, Interval const& b) {
    return {a.lower - b.upper, a.upper - b.lower};
}

Interval negative(Interval const& value) {
    return {Fraction{} - value.upper, Fraction{} - value.lower};
}

Interval multiply(Interval const& a, Interval const& b) {
    const Fraction candidates[4] = {
        a.lower * b.lower, a.lower * b.upper,
        a.upper * b.lower, a.upper * b.upper
    };
    Fraction lower(candidates[0]);
    Fraction upper(candidates[0]);
    for (unsigned i = 1; i < 4; ++i) {
        if (compare(candidates[i], lower) < 0) {
            lower = candidates[i];
        }
        if (compare(candidates[i], upper) > 0) {
            upper = candidates[i];
        }
    }
    return {std::move(lower), std::move(upper)};
}

Interval divide(Interval const& numerator, Interval const& denominator) {
    if (denominator.lower.sign() <= 0 && denominator.upper.sign() >= 0) {
        throw std::domain_error("FloatAP trigonometric interval crosses a pole");
    }
    const Interval reciprocal = {
        Fraction(1L) / denominator.upper,
        Fraction(1L) / denominator.lower
    };
    return multiply(numerator, reciprocal);
}

Fraction tolerance(unsigned long bits) {
    mpz_t one;
    mpz_init_set_ui(one, 1);
    const Fraction result(one, bits);
    mpz_clear(one);
    return result;
}

// For |x| <= 1, both sin and cos Taylor terms decrease in magnitude.
// The absolute value of the first omitted term bounds the entire
// alternating remainder; using a symmetric enclosure also handles x<0.
std::pair<Interval, Interval> small_sine_cosine(Fraction const& x,
                                               unsigned long work) {
    const Fraction square = x * x;
    const Fraction limit = tolerance(work);
    Fraction sine_term(x);
    Fraction cosine_term(1L);
    Fraction sine_sum;
    Fraction cosine_sum;

    for (unsigned long n = 0; ; ++n) {
        if (n > static_cast<unsigned long>(std::numeric_limits<long>::max() / 4)) {
            throw std::overflow_error("FloatAP trigonometric series index");
        }

        if (n % 2 == 0) {
            sine_sum = sine_sum + sine_term;
            cosine_sum = cosine_sum + cosine_term;
        } else {
            sine_sum = sine_sum - sine_term;
            cosine_sum = cosine_sum - cosine_term;
        }

        const long index = static_cast<long>(n);
        sine_term = sine_term * square /
                    Fraction((2 * index + 2) * (2 * index + 3));
        cosine_term = cosine_term * square /
                      Fraction((2 * index + 1) * (2 * index + 2));

        if (compare(absolute(sine_term), limit) < 0 &&
            compare(absolute(cosine_term), limit) < 0) {
            return {
                enclosure(sine_sum, absolute(sine_term)),
                enclosure(cosine_sum, absolute(cosine_term))
            };
        }
    }
}

// For |x| <= 1/2, the remainder of the alternating atan series is
// bounded by the first omitted term. Its sign can be ignored because
// we enclose the answer symmetrically.
Interval small_atan(Fraction const& x, unsigned long work) {
    const Fraction square = x * x;
    const Fraction limit = tolerance(work);
    Fraction power(x);
    Fraction sum;

    for (unsigned long n = 0; ; ++n) {
        if (n > static_cast<unsigned long>(std::numeric_limits<long>::max() / 4)) {
            throw std::overflow_error("FloatAP arctangent series index");
        }

        const long index = static_cast<long>(n);
        const Fraction term = power / Fraction(2 * index + 1);
        sum = n % 2 == 0 ? sum + term : sum - term;
        power = power * square;

        const Fraction error = absolute(power / Fraction(2 * index + 3));
        if (compare(error, limit) < 0) {
            return enclosure(sum, error);
        }
    }
}

Interval pi_interval(unsigned long work) {
    // Machin: pi = 16 atan(1/5) - 4 atan(1/239).
    const Interval a = small_atan(Fraction(1L) / Fraction(5L), work);
    const Interval b = small_atan(Fraction(1L) / Fraction(239L), work);
    return subtract(multiply(point(Fraction(16L)), a),
                    multiply(point(Fraction(4L)), b));
}

Interval atan_interval(Fraction const& x, unsigned long work) {
    const Fraction half = Fraction(1L) / Fraction(2L);
    if (compare(absolute(x), half) <= 0) {
        return small_atan(x, work);
    }
    if (x.sign() < 0) {
        return negative(atan_interval(Fraction{} - x, work));
    }

    const Interval quarter_pi = divide(pi_interval(work), point(Fraction(4L)));
    if (compare(x, Fraction(2L)) <= 0) {
        // atan(x) = pi/4 + atan((x-1)/(x+1)), |argument| <= 1/3.
        const Fraction reduced = (x - Fraction(1L)) / (x + Fraction(1L));
        return add(quarter_pi, small_atan(reduced, work));
    }
    // atan(x) = pi/2 - atan(1/x) for x > 2.
    return subtract(multiply(point(Fraction(2L)), quarter_pi),
                    small_atan(Fraction(1L) / x, work));
}

std::pair<Interval, Interval> sine_cosine_intervals(Fraction x,
                                                   unsigned long work) {
    unsigned halvings = 0;
    while (compare(absolute(x), Fraction(1L)) > 0) {
        if (++halvings > 4096u) {
            throw std::overflow_error("FloatAP trigonometric argument too large");
        }
        x = x / Fraction(2L);
    }

    auto results = small_sine_cosine(x, work + halvings + 16u);
    for (unsigned n = 0; n < halvings; ++n) {
        const Interval sine = multiply(point(Fraction(2L)),
                                       multiply(results.first, results.second));
        const Interval cosine = subtract(multiply(results.second, results.second),
                                          multiply(results.first, results.first));
        results = {sine, cosine};
    }
    return results;
}

Interval square_root_interval(Fraction const& x, unsigned long work) {
    if (x.sign() < 0) {
        throw std::domain_error("FloatAP inverse trigonometric domain");
    }
    if (x.sign() == 0) {
        return point(Fraction{});
    }
    if (work > std::numeric_limits<mp_bitcnt_t>::max() / 2u) {
        throw std::overflow_error("FloatAP radical precision too large");
    }

    mpz_t scaled, quotient, root, squared;
    mpz_inits(scaled, quotient, root, squared, nullptr);
    struct Cleanup {
        mpz_ptr a, b, c, d;
        ~Cleanup() { mpz_clears(a, b, c, d, nullptr); }
    } cleanup{scaled, quotient, root, squared};

    mpz_mul_2exp(scaled, mpq_numref(x.raw()), 2u * work);
    mpz_fdiv_q(quotient, scaled, mpq_denref(x.raw()));
    mpz_sqrt(root, quotient);
    const Fraction lower(root, work);

    mpz_mul(squared, root, root);
    mpz_mul(squared, squared, mpq_denref(x.raw()));
    if (mpz_cmp(squared, scaled) == 0) {
        return point(lower);
    }

    mpz_add_ui(root, root, 1u);
    return {lower, Fraction(root, work)};
}

Interval asin_interval(Fraction const& x, unsigned long work) {
    // asin(x) = 2 atan(x / (1 + sqrt(1-x*x))).
    // Unlike atan(x/sqrt(1-x*x)), the denominator remains >= 1.
    const Interval radical = square_root_interval(
        Fraction(1L) - x * x, work);
    const Interval denominator = add(point(Fraction(1L)), radical);
    const Interval argument = divide(point(x), denominator);
    const Interval left = atan_interval(argument.lower, work);
    const Interval right = atan_interval(argument.upper, work);
    return {Fraction(2L) * left.lower, Fraction(2L) * right.upper};
}

enum class Function {
    sine, cosine, tangent, arctangent, arcsine, arccosine,
    pi, half_pi, negative_half_pi
};

Interval evaluate(Function function, Fraction const& x, unsigned long work) {
    switch (function) {
        case Function::sine:
            return sine_cosine_intervals(x, work).first;
        case Function::cosine:
            return sine_cosine_intervals(x, work).second;
        case Function::tangent: {
            const auto [sine, cosine] = sine_cosine_intervals(x, work);
            return divide(sine, cosine);
        }
        case Function::arctangent:
            return atan_interval(x, work);
        case Function::arcsine:
            return asin_interval(x, work);
        case Function::arccosine:
            return subtract(divide(pi_interval(work), point(Fraction(2L))),
                            asin_interval(x, work));
        case Function::pi:
            return pi_interval(work);
        case Function::half_pi:
            return divide(pi_interval(work), point(Fraction(2L)));
        case Function::negative_half_pi:
            return negative(divide(pi_interval(work), point(Fraction(2L))));
        default:
            throw std::invalid_argument("FloatAP unknown trigonometric function");
    }
}

// A result is returned only when all real values in its certified
// enclosure round to exactly the same FloatAP. This includes directed
// rounding and correctly handles the nearest-even halfway boundaries.
FloatAP certified(Function function, Fraction const& x,
                  FloatAP::RoundingModeType mode, FloatAP::PrecisionType precision) {
    const unsigned long bits = static_cast<unsigned long>(precision.bits());
    unsigned long extra = 32u;
    for (;;) {
        if (extra > std::numeric_limits<unsigned long>::max() - bits - 32u) {
            throw std::overflow_error("FloatAP trigonometric precision too large");
        }

        const Interval bounds = evaluate(function, x, bits + extra);
        const FloatAP low = FloatAP::from_rational(bounds.lower.raw(), bits, mode);
        const FloatAP high = FloatAP::from_rational(bounds.upper.raw(), bits, mode);
        if (low == high) {
            return low;
        }

        if (extra > std::numeric_limits<unsigned long>::max() / 2u) {
            throw std::overflow_error("FloatAP trigonometric refinement overflow");
        }
        extra *= 2u;
    }
}

Fraction rational_argument(FloatAP const& x) {
    const Rational rational = static_cast<Rational>(x);
    return Fraction(rational.get_mpq());
}

FloatAP unary(Function function, FloatAP::RoundingModeType mode,
              FloatAP const& x) {
    const unsigned long bits = static_cast<unsigned long>(x.precision().bits());
    if (x.is_nan()) {
        return FloatAP::nan(bits);
    }

    if (x.is_inf()) {
        if (function == Function::arctangent) {
            return certified(x.negative() ? Function::negative_half_pi :
                           Function::half_pi, Fraction{}, mode, x.precision());
        }
        return FloatAP::nan(bits);
    }

    if (x.is_zero()) {
        switch (function) {
            case Function::sine:
            case Function::tangent:
            case Function::arctangent:
            case Function::arcsine:
                return FloatAP::zero(bits, x.negative());
            case Function::cosine:
                return FloatAP(std::int64_t{1}, x.precision());
            case Function::arccosine:
                return certified(Function::arccosine, Fraction{}, mode, x.precision());
            case Function::pi:
            case Function::half_pi:
            case Function::negative_half_pi:
            default:
                break;
        }
    }

    const Fraction argument = rational_argument(x);
    if (function == Function::arcsine || function == Function::arccosine) {
        if (compare(absolute(argument), Fraction(1L)) > 0) {
            return FloatAP::nan(bits);
        }
        if (compare(argument, Fraction(1L)) == 0) {
            if (function == Function::arccosine) {
                return FloatAP::zero(bits);
            }
            return certified(Function::half_pi, Fraction{}, mode, x.precision());
        }
        if (compare(argument, Fraction(-1L)) == 0) {
            return certified(function == Function::arcsine ?
                             Function::negative_half_pi : Function::pi,
                             Fraction{}, mode, x.precision());
        }
    }
    return certified(function, argument, mode, x.precision());
}

} // namespace

FloatAP sin(FloatAP::RoundingModeType mode, FloatAP const& x) {
    return unary(Function::sine, mode, x);
}

FloatAP cos(FloatAP::RoundingModeType mode, FloatAP const& x) {
    return unary(Function::cosine, mode, x);
}

FloatAP tan(FloatAP::RoundingModeType mode, FloatAP const& x) {
    return unary(Function::tangent, mode, x);
}

FloatAP atan(FloatAP::RoundingModeType mode, FloatAP const& x) {
    return unary(Function::arctangent, mode, x);
}

FloatAP asin(FloatAP::RoundingModeType mode, FloatAP const& x) {
    return unary(Function::arcsine, mode, x);
}

FloatAP acos(FloatAP::RoundingModeType mode, FloatAP const& x) {
    return unary(Function::arccosine, mode, x);
}

FloatAP Float<ArbitraryPrecision>::pi(RoundingModeType mode, PrecisionType precision) {
    return certified(Function::pi, Fraction{}, mode, precision);
}

} // namespace Ariadne
