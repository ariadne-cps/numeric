// Part of Ariadne Numeric. GPL-3.0-or-later.
#ifndef ARIADNE_NUMERIC_DETAIL_FLOATAP_KERNEL_HPP
#define ARIADNE_NUMERIC_DETAIL_FLOATAP_KERNEL_HPP

#include <gmp.h>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Ariadne::detail {

// Implementation kernel for Float<ArbitraryPrecision>, independent of MPFR.
// A finite nonzero value is (-1)^negative * significand * 2^exponent,
// with precisely precision significant bits. GMP is temporary limb storage.
class FloatAPKernel {
public:
    enum class Rounding { nearest, downward, upward, toward_zero };
    enum class Kind { finite, infinity, nan };
    using Exponent = std::int64_t;

private:
    mpz_t _significand;
    unsigned long _precision;
    Exponent _exponent = 0;
    Kind _kind = Kind::finite;
    bool _negative = false;

    static Exponent checked_add(Exponent a, Exponent b) {
        constexpr auto max=std::numeric_limits<Exponent>::max();
        constexpr auto min=std::numeric_limits<Exponent>::min();
        if ((b>0 && a>max-b) || (b<0 && a<min-b)) {
            throw std::overflow_error("FloatAP exponent overflow");
        }
        return a+b;
    }
    static Exponent checked_sub(Exponent a, Exponent b) {
        if (b==std::numeric_limits<Exponent>::min()) {
            return checked_add(checked_add(a,std::numeric_limits<Exponent>::max()),1);
        }
        return checked_add(a,-b);
    }
    static mp_bitcnt_t shift_amount(Exponent n) {
        if (n<0 || static_cast<std::uint64_t>(n)>std::numeric_limits<mp_bitcnt_t>::max()) {
            throw std::overflow_error("FloatAP shift out of range");
        }
        return static_cast<mp_bitcnt_t>(n);
    }
    static Exponent bit_length(mpz_srcptr x) {
        const auto length=mpz_sizeinbase(x,2);
        if (length>static_cast<std::uint64_t>(std::numeric_limits<Exponent>::max())) {
            throw std::overflow_error("FloatAP operand too large");
        }
        return static_cast<Exponent>(length);
    }
    void set_zero(bool negative) {
        _kind=Kind::finite;
        _negative=negative;
        _exponent=0;
        mpz_set_ui(_significand,0);
    }

public:
    explicit FloatAPKernel(unsigned long precision) : _precision(precision) {
        if (precision<2 || precision>static_cast<std::uint64_t>(std::numeric_limits<Exponent>::max())) {
            throw std::invalid_argument("FloatAP precision must be >=2 and fit in exponent range");
        }
        mpz_init(_significand);
    }
    ~FloatAPKernel() { mpz_clear(_significand); }
    FloatAPKernel(FloatAPKernel const& x)
      : _precision(x._precision), _exponent(x._exponent), _kind(x._kind), _negative(x._negative) {
        mpz_init_set(_significand,x._significand);
    }
    FloatAPKernel(FloatAPKernel&& x) noexcept
      : _precision(x._precision), _exponent(x._exponent), _kind(x._kind), _negative(x._negative) {
        mpz_init(_significand);
        mpz_swap(_significand,x._significand);
    }
    FloatAPKernel& operator=(FloatAPKernel x) noexcept { swap(x); return *this; }
    void swap(FloatAPKernel& x) noexcept {
        mpz_swap(_significand,x._significand);
        std::swap(_precision,x._precision);
        std::swap(_exponent,x._exponent);
        std::swap(_kind,x._kind);
        std::swap(_negative,x._negative);
    }

    unsigned long precision() const noexcept { return _precision; }
    Exponent exponent() const noexcept { return _exponent; }
    mpz_srcptr significand() const noexcept { return _significand; }
    Kind kind() const noexcept { return _kind; }
    bool negative() const noexcept { return _negative; }
    bool is_zero() const noexcept { return _kind==Kind::finite && mpz_sgn(_significand)==0; }
    bool is_finite() const noexcept { return _kind==Kind::finite; }
    bool is_inf() const noexcept { return _kind==Kind::infinity; }
    bool is_nan() const noexcept { return _kind==Kind::nan; }

    static FloatAPKernel zero(unsigned long p, bool negative=false) {
        FloatAPKernel z(p); z.set_zero(negative); return z;
    }
    static FloatAPKernel infinity(unsigned long p, bool negative=false) {
        FloatAPKernel z(p); z._kind=Kind::infinity; z._negative=negative; return z;
    }
    static FloatAPKernel nan(unsigned long p) {
        FloatAPKernel z(p); z._kind=Kind::nan; return z;
    }

    // Round signed numerator/positive denominator * 2^scale with one rounding.
    static FloatAPKernel from_ratio(mpz_srcptr numerator, mpz_srcptr denominator,
                                    Exponent scale, unsigned long precision,
                                    Rounding mode, bool negative_zero=false);
    static FloatAPKernel from_integer(mpz_srcptr n, unsigned long precision, Rounding mode);
    static FloatAPKernel from_rational(mpq_srcptr value, unsigned long precision, Rounding mode);
    // Exact conversion, not an approximation to double.
    void to_rational(mpq_ptr out) const;
    static FloatAPKernel negate(FloatAPKernel x);
    static FloatAPKernel add(FloatAPKernel const& a, FloatAPKernel const& b,
                             unsigned long precision, Rounding mode);
    static FloatAPKernel sub(FloatAPKernel const& a, FloatAPKernel const& b,
                             unsigned long precision, Rounding mode);
    static FloatAPKernel mul(FloatAPKernel const& a, FloatAPKernel const& b,
                             unsigned long precision, Rounding mode);
    static FloatAPKernel div(FloatAPKernel const& a, FloatAPKernel const& b,
                             unsigned long precision, Rounding mode);
    // Exact product followed by a single rounding after addition.
    static FloatAPKernel fma(FloatAPKernel const& a, FloatAPKernel const& b,
                             FloatAPKernel const& c, unsigned long precision, Rounding mode);
};

} // namespace Ariadne::detail
#endif
