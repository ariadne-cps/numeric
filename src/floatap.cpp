// Part of Ariadne Numeric. GPL-3.0-or-later.
#include "numeric/detail/floatap_kernel.hpp"
namespace Ariadne::detail {
FloatAPKernel FloatAPKernel::from_ratio(mpz_srcptr numerator, mpz_srcptr denominator,
    Exponent scale, unsigned long precision, Rounding mode, bool negative_zero) {
    FloatAPKernel result(precision);
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
FloatAPKernel FloatAPKernel::from_integer(mpz_srcptr n,unsigned long p,Rounding mode) {
    mpz_t one;mpz_init_set_ui(one,1);
    struct Cleanup{mpz_ptr z;~Cleanup(){mpz_clear(z);}} cleanup{one};
    return from_ratio(n,one,0,p,mode);
}
FloatAPKernel FloatAPKernel::from_rational(mpq_srcptr r,unsigned long p,Rounding mode) {
    return from_ratio(mpq_numref(r),mpq_denref(r),0,p,mode);
}
void FloatAPKernel::to_rational(mpq_ptr out) const {
    if(!is_finite()) throw std::domain_error("FloatAP nonfinite to Rational");
    mpz_set(mpq_numref(out),_significand);
    if(_negative) mpz_neg(mpq_numref(out),mpq_numref(out));
    mpz_set_ui(mpq_denref(out),1);
    if(_exponent>0) mpz_mul_2exp(mpq_numref(out),mpq_numref(out),shift_amount(_exponent));
    if(_exponent<0) mpz_mul_2exp(mpq_denref(out),mpq_denref(out),shift_amount(checked_sub(0,_exponent)));
    mpq_canonicalize(out);
}
FloatAPKernel FloatAPKernel::negate(FloatAPKernel x) {
    if(!x.is_nan()) x._negative=!x._negative;
    return x;
}
FloatAPKernel FloatAPKernel::add(FloatAPKernel const& a,FloatAPKernel const& b,
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
FloatAPKernel FloatAPKernel::sub(FloatAPKernel const& a,FloatAPKernel const& b,
    unsigned long p,Rounding mode) {return add(a,negate(b),p,mode);}
FloatAPKernel FloatAPKernel::mul(FloatAPKernel const& a,FloatAPKernel const& b,
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
FloatAPKernel FloatAPKernel::div(FloatAPKernel const& a,FloatAPKernel const& b,
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
FloatAPKernel FloatAPKernel::fma(FloatAPKernel const& a,FloatAPKernel const& b,
    FloatAPKernel const& c,unsigned long p,Rounding mode) {
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
    FloatAPKernel exact(static_cast<unsigned long>(mpz_sizeinbase(product,2)));
    mpz_set(exact._significand,product);
    exact._negative=neg;
    exact._exponent=checked_add(a._exponent,b._exponent);
    return add(exact,c,p,mode);
}
} // namespace Ariadne::detail
