/***************************************************************************
 *            numeric/floatap.cpp
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

/*! \file numeric/floatap.cpp
 *  \brief Native arbitrary-precision floating-point operations.
 */

#include "numeric/floatap.hpp"
#include "numeric/integer.hpp"
#include "numeric/rational.hpp"
namespace Ariadne {
FloatAP FloatAP::from_ratio(mpz_srcptr numerator, mpz_srcptr denominator,
    Exponent scale, unsigned long precision, Rounding mode, bool negative_zero) {
    FloatAP result(precision);
    if (mpz_sgn(denominator)==0) throw std::domain_error("FloatAP zero denominator");
    if (mpz_sgn(numerator)==0) { result.set_zero(negative_zero); return result; }
    result._negative=(mpz_sgn(numerator)<0)!=(mpz_sgn(denominator)<0);
    mpz_t n,d,q,r,w; mpz_inits(n,d,q,r,w,nullptr);
    struct Cleanup { mpz_ptr a,b,c,d,e; ~Cleanup(){mpz_clears(a,b,c,d,e,nullptr);} } cleanup{n,d,q,r,w};
    mpz_abs(n,numerator); mpz_abs(d,denominator);
    Exponent k=checked_sub(bit_length(n),bit_length(d));
    if (k>=0) { mpz_mul_2exp(w,d,shift_amount(k)); if(mpz_cmp(n,w)<0) k=checked_sub(k,1); }
    else { mpz_mul_2exp(w,n,shift_amount(checked_sub(0,k))); if(mpz_cmp(w,d)<0) k=checked_sub(k,1); }
    const Exponent shift=checked_sub(static_cast<Exponent>(precision-1),k);
    if (shift>=0) {
        mpz_mul_2exp(w,n,shift_amount(shift)); mpz_tdiv_qr(q,r,w,d); mpz_set(w,d);
    } else {
        mpz_mul_2exp(w,d,shift_amount(checked_sub(0,shift))); mpz_tdiv_qr(q,r,n,w);
    }
    bool increment=false;
    if (mpz_sgn(r)) {
        switch(mode) {
        case Rounding::nearest:
            mpz_mul_2exp(r,r,1);
            increment=mpz_cmp(r,w)>0 || (mpz_cmp(r,w)==0 && mpz_odd_p(q));
            break;
        case Rounding::downward: increment=result._negative; break;
        case Rounding::upward: increment=!result._negative; break;
        case Rounding::toward_zero: break;
        default:
            throw std::invalid_argument("FloatAP unknown rounding mode");
        }
    }
    if (increment) mpz_add_ui(q,q,1);
    result._exponent=checked_add(scale,checked_sub(k,static_cast<Exponent>(precision-1)));
    if (bit_length(q)>static_cast<Exponent>(precision)) {
        mpz_tdiv_q_2exp(q,q,1); result._exponent=checked_add(result._exponent,1);
    }
    mpz_set(result._significand,q); return result;
}
FloatAP FloatAP::from_integer(mpz_srcptr n,unsigned long p,Rounding mode) {
    mpz_t one;mpz_init_set_ui(one,1);
    struct Cleanup{mpz_ptr z;~Cleanup(){mpz_clear(z);}} cleanup{one};
    return from_ratio(n,one,0,p,mode);
}
FloatAP FloatAP::from_rational(mpq_srcptr r,unsigned long p,Rounding mode) {
    return from_ratio(mpq_numref(r),mpq_denref(r),0,p,mode);
}
void FloatAP::to_rational(mpq_ptr out) const {
    if(!is_finite()) throw std::domain_error("FloatAP nonfinite to Rational");
    mpz_set(mpq_numref(out),_significand);
    if(_negative) mpz_neg(mpq_numref(out),mpq_numref(out));
    mpz_set_ui(mpq_denref(out),1);
    if(_exponent>0) mpz_mul_2exp(mpq_numref(out),mpq_numref(out),shift_amount(_exponent));
    if(_exponent<0) mpz_mul_2exp(mpq_denref(out),mpq_denref(out),shift_amount(checked_sub(0,_exponent)));
    mpq_canonicalize(out);
}
FloatAP FloatAP::negate(FloatAP x) {
    if(!x.is_nan()) x._negative=!x._negative;
    return x;
}
FloatAP FloatAP::add(FloatAP const& a,FloatAP const& b,
    unsigned long p,Rounding mode) {
    if(a.is_nan()||b.is_nan()) return nan(p);
    if(a.is_inf()||b.is_inf()) {
        if(a.is_inf()&&b.is_inf()&&a._negative!=b._negative) return nan(p);
        return infinity(p,a.is_inf()?a._negative:b._negative);
    }
    if(a.is_zero()&&b.is_zero()) return zero(p,a._negative==b._negative?a._negative:mode==Rounding::downward);
    const Exponent low=a._exponent<b._exponent?a._exponent:b._exponent;
    mpz_t x,y,sum,one;mpz_inits(x,y,sum,one,nullptr);
    struct Cleanup{mpz_ptr a,b,c,d;~Cleanup(){mpz_clears(a,b,c,d,nullptr);}} cleanup{x,y,sum,one};
    mpz_mul_2exp(x,a._significand,shift_amount(checked_sub(a._exponent,low)));
    mpz_mul_2exp(y,b._significand,shift_amount(checked_sub(b._exponent,low)));
    if(a._negative) mpz_neg(x,x);
    if(b._negative) mpz_neg(y,y);
    mpz_add(sum,x,y);mpz_set_ui(one,1);
    return from_ratio(sum,one,low,p,mode,mode==Rounding::downward);
}
FloatAP FloatAP::sub(FloatAP const& a,FloatAP const& b,
    unsigned long p,Rounding mode) {return add(a,negate(b),p,mode);}
FloatAP FloatAP::mul(FloatAP const& a,FloatAP const& b,
    unsigned long p,Rounding mode) {
    if(a.is_nan()||b.is_nan()||(a.is_inf()&&b.is_zero())||(a.is_zero()&&b.is_inf())) return nan(p);
    const bool neg=a._negative!=b._negative;
    if(a.is_inf()||b.is_inf()) return infinity(p,neg);
    if(a.is_zero()||b.is_zero()) return zero(p,neg);
    mpz_t n,one;mpz_inits(n,one,nullptr);
    struct Cleanup{mpz_ptr a,b;~Cleanup(){mpz_clears(a,b,nullptr);}} cleanup{n,one};
    mpz_mul(n,a._significand,b._significand);
    if(neg) mpz_neg(n,n);
    mpz_set_ui(one,1);
    return from_ratio(n,one,checked_add(a._exponent,b._exponent),p,mode);
}
FloatAP FloatAP::div(FloatAP const& a,FloatAP const& b,
    unsigned long p,Rounding mode) {
    if(a.is_nan()||b.is_nan()||(a.is_zero()&&b.is_zero())||(a.is_inf()&&b.is_inf())) return nan(p);
    const bool neg=a._negative!=b._negative;
    if(a.is_inf()||b.is_zero()) return infinity(p,neg);
    if(a.is_zero()||b.is_inf()) return zero(p,neg);
    mpz_t n;mpz_init_set(n,a._significand);
    struct Cleanup{mpz_ptr z;~Cleanup(){mpz_clear(z);}} cleanup{n};
    if(neg) mpz_neg(n,n);
    return from_ratio(n,b._significand,checked_sub(a._exponent,b._exponent),p,mode);
}
FloatAP FloatAP::fma(FloatAP const& a,FloatAP const& b,
    FloatAP const& c,unsigned long p,Rounding mode) {
    if(a.is_nan()||b.is_nan()||c.is_nan()||(a.is_zero()&&b.is_inf())||(a.is_inf()&&b.is_zero())) return nan(p);
    const bool neg=a._negative!=b._negative;
    if(a.is_inf()||b.is_inf()) {
        if(c.is_inf()&&c._negative!=neg) return nan(p);
        return infinity(p,neg);
    }
    if(c.is_inf()) return infinity(p,c._negative);
    if(a.is_zero()||b.is_zero()) return add(zero(p,neg),c,p,mode);
    mpz_t product;mpz_init(product);
    struct Cleanup{mpz_ptr z;~Cleanup(){mpz_clear(z);}} cleanup{product};
    mpz_mul(product,a._significand,b._significand);
    FloatAP exact(static_cast<unsigned long>(mpz_sizeinbase(product,2)));
    mpz_set(exact._significand,product);
    exact._negative=neg;
    exact._exponent=checked_add(a._exponent,b._exponent);
    return add(exact,c,p,mode);
}


namespace {
// A signed 64-bit integer is imported directly, without a decimal
// intermediate. In particular, INT64_MIN is handled without negating
// a signed value outside its range.
FloatAP float_from_int64(std::int64_t n, ArbitraryPrecision precision) {
    mpz_t value;
    mpz_init(value);
    struct Cleanup {
        mpz_ptr value;
        ~Cleanup() { mpz_clear(value); }
    } cleanup{value};

    std::uint64_t magnitude=static_cast<std::uint64_t>(n);
    if (n<0) {
        magnitude=std::uint64_t{0}-magnitude;
    }
    mpz_import(value,1,1,sizeof(magnitude),0,0,&magnitude);
    if (n<0) {
        mpz_neg(value,value);
    }
    return FloatAP::from_integer(value,static_cast<unsigned long>(precision.bits()),
                                 FloatAP::ROUND_TO_NEAREST);
}
}

FloatAP::Float(std::int64_t n, PrecisionType pr)
    : Float(float_from_int64(n,pr)) { }

FloatAP::Float(Integer const& n, RoundingModeType mode, PrecisionType pr)
    : Float(from_integer(n.get_mpz(),static_cast<unsigned long>(pr.bits()),mode)) { }

namespace {
// Rational declares is_nan/is_inf as friend functions: ADL finds them,
// whereas qualified lookup (Ariadne::is_nan) does not.
bool rational_is_nan(Rational const& q) { return is_nan(q); }
bool rational_is_inf(Rational const& q) { return is_inf(q); }
}

FloatAP::Float(Rational const& q, RoundingModeType mode, PrecisionType pr)
    : Float(pr) {
    if (rational_is_nan(q)) {
        *this=nan(_precision);
    } else if (rational_is_inf(q)) {
        *this=infinity(_precision,mpz_sgn(mpq_numref(q.get_mpq()))<0);
    } else {
        *this=from_rational(q.get_mpq(),_precision,mode);
    }
}

FloatAP::operator Rational() const {
    if (is_nan()) return Rational::nan();
    if (is_inf()) return Rational::inf(negative()?Sign::NEGATIVE:Sign::POSITIVE);
    Rational q;
    to_rational(q._mpq);
    return q;
}

double FloatAP::get_d() const {
    return static_cast<Rational>(*this).get_d();
}

Comparison cmp(FloatAP const& a, FloatAP const& b) {
    if (a.is_nan()||b.is_nan()) return Comparison::INCOMPARABLE;
    if (a.is_inf()||b.is_inf()) {
        if (a.is_inf()&&b.is_inf()) {
            return a.negative()==b.negative()?Comparison::EQUAL:
                (a.negative()?Comparison::LESS:Comparison::GREATER);
        }
        return a.is_inf()?(a.negative()?Comparison::LESS:Comparison::GREATER):
                          (b.negative()?Comparison::GREATER:Comparison::LESS);
    }
    mpq_t x,y;
    mpq_inits(x,y,nullptr);
    struct Cleanup { mpq_ptr x,y; ~Cleanup(){mpq_clears(x,y,nullptr);} } cleanup{x,y};
    a.to_rational(x);
    b.to_rational(y);
    const int order=mpq_cmp(x,y);
    return order<0?Comparison::LESS:(order>0?Comparison::GREATER:Comparison::EQUAL);
}

} // namespace Ariadne
