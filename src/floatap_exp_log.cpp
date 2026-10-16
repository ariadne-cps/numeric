/***************************************************************************
 *            numeric/floatap_exp_log.cpp
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
#include "numeric/rational.hpp"

#include <limits>
#include <stdexcept>
#include <utility>

namespace Ariadne {
namespace {

// All approximations are exact GMP rationals; MPFR is never called.
class Fraction {
    mpq_t _q;

public:
    Fraction() { mpq_init(_q); }
    explicit Fraction(long n) : Fraction() { mpq_set_si(_q, n, 1); }
    explicit Fraction(mpq_srcptr q) : Fraction() { mpq_set(_q, q); }
    Fraction(Fraction const& other) : Fraction(other._q) { }
    Fraction(Fraction&& other) noexcept : Fraction() { mpq_swap(_q, other._q); }
    ~Fraction() { mpq_clear(_q); }

    Fraction& operator=(Fraction other) noexcept {
        mpq_swap(_q, other._q);
        return *this;
    }

    mpq_srcptr raw() const noexcept { return _q; }
    int sign() const noexcept { return mpq_sgn(_q); }

    friend int compare(Fraction const& a, Fraction const& b) {
        return mpq_cmp(a._q, b._q);
    }
    friend Fraction operator+(Fraction const& a, Fraction const& b) {
        Fraction result;
        mpq_add(result._q, a._q, b._q);
        return result;
    }
    friend Fraction operator-(Fraction const& a, Fraction const& b) {
        Fraction result;
        mpq_sub(result._q, a._q, b._q);
        return result;
    }
    friend Fraction operator*(Fraction const& a, Fraction const& b) {
        Fraction result;
        mpq_mul(result._q, a._q, b._q);
        return result;
    }
    friend Fraction operator/(Fraction const& a, Fraction const& b) {
        Fraction result;
        mpq_div(result._q, a._q, b._q);
        return result;
    }
};

struct Interval {
    Fraction low;
    Fraction high;
};

Fraction abs_fraction(Fraction const& x) {
    return x.sign() < 0 ? Fraction{} - x : x;
}

Fraction epsilon(unsigned long bits) {
    mpq_t q;
    mpq_init(q);
    mpq_set_ui(q, 1, 1);
    mpz_mul_2exp(mpq_denref(q), mpq_denref(q), bits);
    Fraction result(q);
    mpq_clear(q);
    return result;
}

Interval centered(Fraction const& value, Fraction const& error) {
    return {value - error, value + error};
}

// For |x| <= 1/2 the exponential tail after term n is at most
// twice the absolute value of the next term.
Interval exp_series(Fraction const& x, unsigned long work) {
    Fraction term(1);
    Fraction sum(1);
    const Fraction limit = epsilon(work);
    for (unsigned long n = 1; ; ++n) {
        if (n == std::numeric_limits<long>::max()) {
            throw std::overflow_error("FloatAP exponential series index");
        }
        term = term * x / Fraction(static_cast<long>(n));
        sum = sum + term;
        const Fraction next = term * x / Fraction(static_cast<long>(n + 1));
        const Fraction error = Fraction(2) * abs_fraction(next);
        if (compare(error, limit) < 0) {
            return centered(sum, error);
        }
    }
}

Interval exp_interval(Fraction x, unsigned long work) {
    unsigned halvings = 0;
    while (compare(abs_fraction(x), Fraction(1) / Fraction(2)) > 0) {
        if (++halvings > 2048) {
            throw std::overflow_error("FloatAP exp argument range");
        }
        x = x / Fraction(2);
    }
    Interval interval = exp_series(x, work + halvings + 8);
    for (unsigned i = 0; i < halvings; ++i) {
        // exp(x) is positive, so endpoint multiplication is monotone.
        interval = {interval.low * interval.low,
                    interval.high * interval.high};
    }
    return interval;
}

// For |z| <= 1/3, log((1+z)/(1-z)) = 2 sum z^(2n+1)/(2n+1).
// The tail is bounded by 4*|z^(2n+3)|/(2n+3).
Interval log_series(Fraction const& z, unsigned long work) {
    Fraction power(z);
    Fraction sum;
    const Fraction square = z * z;
    const Fraction limit = epsilon(work);
    for (unsigned long n = 0; ; ++n) {
        if (n >= (static_cast<unsigned long>(std::numeric_limits<long>::max()) - 3) / 2) {
            throw std::overflow_error("FloatAP logarithm series index");
        }
        const long denominator = static_cast<long>(2 * n + 1);
        sum = sum + power / Fraction(denominator);
        power = power * square;
        const Fraction error = Fraction(4) * abs_fraction(power) /
                               Fraction(denominator + 2);
        if (compare(error, limit) < 0) {
            return centered(Fraction(2) * sum, error);
        }
    }
}

Interval log_interval(Fraction x, unsigned long work) {
    long power = 0;
    while (compare(x, Fraction(2)) > 0) {
        if (++power > 2048) {
            throw std::overflow_error("FloatAP log argument range");
        }
        x = x / Fraction(2);
    }
    while (compare(x, Fraction(1) / Fraction(2)) < 0) {
        if (--power < -2048) {
            throw std::overflow_error("FloatAP log argument range");
        }
        x = x * Fraction(2);
    }
    const Fraction z = (x - Fraction(1)) / (x + Fraction(1));
    Interval result = log_series(z, work);
    if (power != 0) {
        const Interval ln2 = log_series(Fraction(1) / Fraction(3), work);
        const Fraction count(power < 0 ? -power : power);
        if (power > 0) {
            result = {result.low + count * ln2.low,
                      result.high + count * ln2.high};
        } else {
            result = {result.low - count * ln2.high,
                      result.high - count * ln2.low};
        }
    }
    return result;
}

template<class Evaluate>
FloatAP certified(Evaluate evaluate, FloatAP::RoundingModeType mode,
                  FloatAP::PrecisionType precision) {
    const unsigned long bits = static_cast<unsigned long>(precision.bits());
    unsigned long extra = 32;
    for (;;) {
        if (extra > std::numeric_limits<unsigned long>::max() - bits - 16) {
            throw std::overflow_error("FloatAP certified precision overflow");
        }
        const Interval interval = evaluate(bits + extra);
        const FloatAP lower = FloatAP::from_rational(interval.low.raw(), bits, mode);
        const FloatAP upper = FloatAP::from_rational(interval.high.raw(), bits, mode);
        if (lower == upper) {
            return lower;
        }
        if (extra > std::numeric_limits<unsigned long>::max() / 2) {
            throw std::overflow_error("FloatAP certified refinement overflow");
        }
        extra *= 2;
    }
}

Fraction exact_argument(FloatAP const& x) {
    const Rational q = static_cast<Rational>(x);
    return Fraction(q.get_mpq());
}

} // namespace

FloatAP exp(FloatAP::RoundingModeType mode, FloatAP const& x) {
    const unsigned long bits = static_cast<unsigned long>(x.precision().bits());
    if (x.is_nan()) {
        return FloatAP::nan(bits);
    }
    if (x.is_inf()) {
        return x.negative() ? FloatAP::zero(bits) : FloatAP::infinity(bits);
    }
    if (x.is_zero()) {
        return FloatAP(std::int64_t{1}, x.precision());
    }
    const Fraction input = exact_argument(x);
    return certified([&](unsigned long work) {
        return exp_interval(input, work);
    }, mode, x.precision());
}

FloatAP log(FloatAP::RoundingModeType mode, FloatAP const& x) {
    const unsigned long bits = static_cast<unsigned long>(x.precision().bits());
    if (x.is_nan() || (x.negative() && !x.is_zero())) {
        return FloatAP::nan(bits);
    }
    if (x.is_zero()) {
        return FloatAP::infinity(bits, true);
    }
    if (x.is_inf()) {
        return FloatAP::infinity(bits);
    }
    const Fraction input = exact_argument(x);
    if (compare(input, Fraction(1)) == 0) {
        return FloatAP::zero(bits);
    }
    return certified([&](unsigned long work) {
        return log_interval(input, work);
    }, mode, x.precision());
}

} // namespace Ariadne
