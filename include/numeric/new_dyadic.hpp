/***************************************************************************
 * new_dyadic.hpp
 * Experimental native arbitrary-precision dyadic; Dyadic remains unchanged.
 ****************************************************************************/
#ifndef ARIADNE_NEW_DYADIC_HPP
#define ARIADNE_NEW_DYADIC_HPP

#include "numeric/integer.hpp"
#include "numeric/rational.hpp"
#include "numeric/decimal.hpp"
#include "numeric/twoexp.hpp"
#include "numeric/builtin.hpp"
#include <type_traits>
#include <array>
#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace Ariadne {

// Exact finite value: coefficient / 2^denominator_exponent.
// Canonical form has no factor of two in the coefficient if exponent>0.
// GMP is deliberately absent from representation and arithmetic.
class NewDyadic {
public:
    enum class Kind : unsigned char { Finite, PositiveInfinity, NegativeInfinity, NaN };
private:
    Integer _coefficient;
    std::size_t _denominator_exponent=0u;
    Kind _kind=Kind::Finite;

    static Integer shifted(Integer const& x,std::size_t bits) {
        if (bits==0u || x.magnitude().is_zero()) { return x; }
        return Integer::from_magnitude(shift_left(x.magnitude(),bits),x.negative());
    }

    void normalize() {
        if (_kind!=Kind::Finite || _coefficient.magnitude().is_zero()) {
            _denominator_exponent=0u;
            return;
        }
        auto const& magnitude=_coefficient.magnitude();
        std::size_t trailing=0u;
        for (std::size_t i=0u;i<magnitude.size();++i) {
            const std::uint64_t limb=magnitude.limb(i);
            if (limb!=0u) {
                trailing+=static_cast<std::size_t>(std::countr_zero(limb));
                break;
            }
            trailing+=64u;
        }
        const std::size_t reduce=std::min(trailing,_denominator_exponent);
        if (reduce!=0u) {
            _coefficient=Integer::from_magnitude(
                shift_right(magnitude,reduce),_coefficient.negative());
            _denominator_exponent-=reduce;
        }
    }

    static bool has_fraction(BigUnsigned const& n,std::size_t exponent) {
        if (exponent==0u || n.is_zero()) { return false; }
        if (exponent>=n.bit_length()) { return true; }
        const auto truncated=shift_right(n,exponent);
        return compare(shift_left(truncated,exponent),n)!=0;
    }

    static Integer truncated(NewDyadic const& x) {
        return Integer::from_magnitude(
            shift_right(x._coefficient.magnitude(),x._denominator_exponent),
            x._coefficient.negative());
    }

    // Exact comparison of positive n/2^en and m/2^em without
    // constructing shifted operands.
    static int compare_scaled(BigUnsigned const& n,std::size_t en,
                              BigUnsigned const& m,std::size_t em) noexcept {
        if (n.is_zero() || m.is_zero()) {
            return n.is_zero() ? (m.is_zero() ? 0 : -1) : 1;
        }
        const std::size_t nb=n.bit_length(),mb=m.bit_length();
        const bool ng=nb>=en,mg=mb>=em;
        if (ng!=mg) { return ng ? 1 : -1; }
        if (ng) {
            if (nb-en!=mb-em) { return nb-en>mb-em ? 1 : -1; }
        } else {
            if (en-nb!=em-mb) { return en-nb<em-mb ? 1 : -1; }
        }
        return nb>=mb ? compare_aligned(n,m) : -compare_aligned(m,n);
    }

    // Align shorter positive magnitude virtually to the longer
    // magnitude's most significant bit. No owning BigUnsigned.
    static int compare_aligned(BigUnsigned const& a,
                               BigUnsigned const& b) noexcept {
        const std::size_t shift=a.bit_length()-b.bit_length();
        const std::size_t whole=shift/64u;
        const unsigned part=static_cast<unsigned>(shift%64u);
        for (std::size_t i=a.size();i!=0u;--i) {
            const std::size_t index=i-1u;
            std::uint64_t word=0u;
            if (index>=whole) {
                const std::size_t j=index-whole;
                if (j<b.size()) {
                    word=part==0u ? b.limb(j) : b.limb(j)<<part;
                }
                if (part!=0u && j!=0u && j-1u<b.size()) {
                    word|=b.limb(j-1u)>>(64u-part);
                }
            }
            const std::uint64_t actual=a.limb(index);
            if (actual!=word) { return actual<word ? -1 : 1; }
        }
        return 0;
    }

    static void require_finite(NewDyadic const& x) {
        if (x._kind!=Kind::Finite) {
            throw std::domain_error("non-finite NewDyadic cannot be rounded");
        }
    }

    static int signum(NewDyadic const& x) noexcept {
        if (x._kind==Kind::PositiveInfinity) { return 1; }
        if (x._kind==Kind::NegativeInfinity) { return -1; }
        if (x._kind==Kind::NaN || x._coefficient.magnitude().is_zero()) { return 0; }
        return x._coefficient.negative() ? -1 : 1;
    }

public:
    NewDyadic()=default;
    explicit NewDyadic(Integer numerator,std::size_t exponent=0u)
        : _coefficient(std::move(numerator)),_denominator_exponent(exponent) {
        normalize();
    }
    NewDyadic(Integer const& numerator,Natural exponent)
        : NewDyadic(numerator,static_cast<std::size_t>(exponent.get_si())) { }
    template<BuiltinSignedIntegral N>
    NewDyadic(Integer const&,N) = delete;
    template<BuiltinIntegral N>
    NewDyadic(N n) : NewDyadic(Integer(n)) { }

    // Same value-conversion entry points as Dyadic; all arithmetic
    // remains native, and only exact dyadic Rational values are accepted.
    explicit NewDyadic(Rational const& q) {
        if (is_nan(q)) { *this=nan(); return; }
        if (is_inf(q)) {
            *this=inf(sgn(q)==Sign::POSITIVE ? 1 : -1);
            return;
        }
        auto const denominator=q.denominator();
        auto const& magnitude=denominator.magnitude();
        if (magnitude.is_zero() || magnitude.bit_length()==0u) {
            throw std::invalid_argument("invalid Dyadic denominator");
        }
        const std::size_t exp=magnitude.bit_length()-1u;
        if (compare(magnitude,shift_left(BigUnsigned(1u),exp))!=0) {
            throw std::invalid_argument("Rational is not a Dyadic");
        }
        _coefficient=q.numerator();
        _denominator_exponent=exp;
        normalize();
    }
    explicit NewDyadic(ExactDouble const& d) : NewDyadic(Rational(d)) { }
    explicit NewDyadic(double d) : NewDyadic(Rational(ExactDouble(d))) { }
    explicit NewDyadic(String const& s) : NewDyadic(Rational(Decimal(s))) { }
    NewDyadic(TwoExp const& power) {
        const Int e=power.exponent();
        if (e>=0) {
            _coefficient=shifted(Integer(1),static_cast<std::size_t>(e));
        } else {
            // Avoid signed overflow for INT_MIN.
            _coefficient=Integer(1);
            _denominator_exponent=static_cast<std::size_t>(-(e+1))+1u;
        }
    }
    Integer mantissa() const { return _coefficient; }
    Int exponent() const {
        if (_denominator_exponent>
              static_cast<std::size_t>(std::numeric_limits<Int>::max())) {
            throw std::overflow_error("NewDyadic exponent out of range");
        }
        return static_cast<Int>(_denominator_exponent);
    }
    String literal() const {
        if (_kind==Kind::NaN) { return String("NaN"); }
        if (_kind==Kind::PositiveInfinity) { return String("inf"); }
        if (_kind==Kind::NegativeInfinity) { return String("-inf"); }
        if (_denominator_exponent==0u) { return _coefficient.literal(); }
        if (_denominator_exponent>std::numeric_limits<Nat>::max()) {
            throw std::overflow_error("NewDyadic literal exponent out of range");
        }
        return Decimal(_coefficient*pow(Integer(5),
            static_cast<Nat>(_denominator_exponent)),
            static_cast<Nat>(_denominator_exponent)).literal();
    }

    static NewDyadic inf(int sign=1) {
        NewDyadic result;
        result._kind=sign==0 ? Kind::NaN :
                     sign>0 ? Kind::PositiveInfinity : Kind::NegativeInfinity;
        return result;
    }
    static NewDyadic nan() { return inf(0); }
    static NewDyadic inf(Sign sign) {
        return inf(sign==Sign::ZERO ? 0 :
                   sign==Sign::POSITIVE ? 1 : -1);
    }

    Kind kind() const noexcept { return _kind; }
    Integer const& coefficient() const noexcept { return _coefficient; }
    std::size_t denominator_exponent() const noexcept {
        return _denominator_exponent;
    }
    bool is_finite() const noexcept { return _kind==Kind::Finite; }
    bool is_zero() const noexcept {
        return is_finite() && _coefficient.magnitude().is_zero();
    }

    Rational to_rational() const {
        if (_kind==Kind::NaN) { return Rational::nan(); }
        if (_kind==Kind::PositiveInfinity) { return Rational::inf(Sign::POSITIVE); }
        if (_kind==Kind::NegativeInfinity) { return Rational::inf(Sign::NEGATIVE); }
        if (_denominator_exponent==0u) { return Rational(_coefficient); }
        const Integer denominator=Integer::from_magnitude(
            shift_left(BigUnsigned(1u),_denominator_exponent),false);
        return Rational(_coefficient,denominator);
    }
    // Truncate to IEEE binary64 from the native magnitude.
    // GMP mpf_get_d returns positive zero on underflow even
    // for negative operands. No Rational allocation.
    double get_d() const {
        if (_kind==Kind::NaN) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (_kind==Kind::PositiveInfinity) {
            return std::numeric_limits<double>::infinity();
        }
        if (_kind==Kind::NegativeInfinity) {
            return -std::numeric_limits<double>::infinity();
        }
        const BigUnsigned& n=_coefficient.magnitude();
        if (n.is_zero()) { return 0.0; }
        const bool negative=_coefficient.negative();
        const std::size_t bits=n.bit_length();
        int exponent;
        if (bits>_denominator_exponent) {
            const std::size_t positive=bits-1u-_denominator_exponent;
            if (positive>1023u) {
                return negative ? -std::numeric_limits<double>::infinity()
                                : std::numeric_limits<double>::infinity();
            }
            exponent=static_cast<int>(positive);
        } else {
            const std::size_t negative_power=_denominator_exponent-bits+1u;
            if (negative_power>1074u) { return 0.0; }
            exponent=-static_cast<int>(negative_power);
        }
        std::uint64_t significand;
        if (bits<=53u) {
            significand=n.limb(0u)<<(53u-bits);
        } else {
            const std::size_t first=bits-53u;
            const std::size_t index=first/64u;
            const unsigned shift=static_cast<unsigned>(first%64u);
            significand=n.limb(index)>>shift;
            if (shift!=0u && index+1u<n.size()) {
                significand|=n.limb(index+1u)<<(64u-shift);
            }
            significand&=(std::uint64_t{1}<<53u)-1u;
        }
        double value;
        if (exponent>=-1022) {
            value=std::ldexp(static_cast<double>(significand),exponent-52);
        } else {
            const unsigned drop=static_cast<unsigned>(-1022-exponent);
            value=std::ldexp(static_cast<double>(significand>>drop),-1074);
        }
        return value==0.0 ? 0.0 : (negative ? -value : value);
    }

    friend NewDyadic operator-(NewDyadic const& x) {
        if (x._kind==Kind::NaN) { return nan(); }
        if (x._kind==Kind::PositiveInfinity) { return inf(-1); }
        if (x._kind==Kind::NegativeInfinity) { return inf(1); }
        return NewDyadic(-x._coefficient,x._denominator_exponent);
    }
    friend NewDyadic operator+(NewDyadic const& x) { return x; }
    friend NewDyadic nul(NewDyadic const&) { return NewDyadic(0); }
    friend NewDyadic pos(NewDyadic const& x) { return x; }
    friend NewDyadic neg(NewDyadic const& x) { return -x; }
    friend NewDyadic add(NewDyadic const& a,NewDyadic const& b) {
        return a+b;
    }
    friend NewDyadic sub(NewDyadic const& a,NewDyadic const& b) {
        return a-b;
    }
    friend NewDyadic mul(NewDyadic const& a,NewDyadic const& b) {
        return a*b;
    }
    friend NewDyadic fma(NewDyadic const& a,NewDyadic const& b,
                          NewDyadic const& c) { return a*b+c; }
    NewDyadic& operator+=(NewDyadic const& rhs) {
        *this=*this+rhs; return *this;
    }
    NewDyadic& operator-=(NewDyadic const& rhs) {
        *this=*this-rhs; return *this;
    }
    NewDyadic& operator*=(NewDyadic const& rhs) {
        *this=*this*rhs; return *this;
    }

  private:
    // R36: on ARM64 the single-buffer accumulator improves add/sub
    // across 5..64-limb magnitudes; x86-64 benefits through 16 limbs
    // but regresses at 32..64 limbs. Keep a single exact algorithm,
    // changing only the architecture-specific crossover to fallback.
    // All other architectures conservatively use the shorter crossover.
    static constexpr std::size_t aligned_add_sub_max_limbs() noexcept {
#if defined(__aarch64__) || defined(_M_ARM64)
        return 64u;
#else
        return 16u;
#endif
    }

    static bool prefer_aligned_add_sub(NewDyadic const& a,
                                       NewDyadic const& b) noexcept {
        const std::size_t as=a._coefficient.magnitude().size();
        const std::size_t bs=b._coefficient.magnitude().size();
        const std::size_t limit=aligned_add_sub_max_limbs();
        const std::size_t max_exponent=std::max(
            a._denominator_exponent,b._denominator_exponent);
        const std::size_t min_exponent=std::min(
            a._denominator_exponent,b._denominator_exponent);
        return (as>=5u || bs>=5u) && as<=limit && bs<=limit &&
               max_exponent-min_exponent<=63u;
    }

    // R35: for aligned finite signed arithmetic on 5..64-limb
    // magnitudes with an exponent gap <=63, scan borrowed virtual
    // shifted limbs and produce one canonical owned magnitude.
    // No shifted Integer operands or temporary negation.
    static NewDyadic aligned_add_sub(NewDyadic const& a,NewDyadic const& b,
                                      bool subtract_second) {
        using Limb=BigUnsigned::Limb;
        const BigUnsigned& an=a._coefficient.magnitude();
        const BigUnsigned& bn=b._coefficient.magnitude();
        if (an.is_zero()) { return subtract_second ? -b : b; }
        if (bn.is_zero()) { return a; }

        const std::size_t exponent=std::max(a._denominator_exponent,
                                            b._denominator_exponent);
        const unsigned shift_a=static_cast<unsigned>(
            exponent-a._denominator_exponent);
        const unsigned shift_b=static_cast<unsigned>(
            exponent-b._denominator_exponent);
        const std::size_t a_bits=an.bit_length()+shift_a;
        const std::size_t b_bits=bn.bit_length()+shift_b;
        const std::size_t width=(std::max(a_bits,b_bits)+63u)/64u;
        // Up to 65 shifted limbs, plus one final addition carry.
        std::array<Limb,66u> output{};

        const auto digit=[](BigUnsigned const& n,std::size_t i,
                            unsigned shift)->Limb {
            if (shift==0u) { return n.limb(i); }
            const Limb low=n.limb(i)<<shift;
            return i==0u ? low : low|(n.limb(i-1u)>>(64u-shift));
        };

        const bool a_negative=a._coefficient.negative();
        const bool b_negative=b._coefficient.negative()!=subtract_second;
        const bool adding=a_negative==b_negative;
        bool negative=a_negative,reverse=false;
        if (!adding) {
            int order=0;
            if (a_bits!=b_bits) {
                order=a_bits>b_bits ? 1 : -1;
            } else {
                for (std::size_t i=width;i!=0u;--i) {
                    const Limb av=digit(an,i-1u,shift_a);
                    const Limb bv=digit(bn,i-1u,shift_b);
                    if (av!=bv) {
                        order=av>bv ? 1 : -1;
                        break;
                    }
                }
            }
            if (order==0) { return NewDyadic(); }
            reverse=order<0;
            negative=reverse ? b_negative : a_negative;
        }

        Limb carry_or_borrow=0u;
        for (std::size_t i=0u;i<width;++i) {
            Limb av=digit(an,i,shift_a),bv=digit(bn,i,shift_b);
            if (reverse) { std::swap(av,bv); }
            if (adding) {
                const Limb partial=av+bv;
                const Limb result=partial+carry_or_borrow;
                carry_or_borrow=static_cast<Limb>(
                    partial<av || result<partial);
                output[i]=result;
            } else {
                const Limb partial=av-bv;
                const Limb result=partial-carry_or_borrow;
                carry_or_borrow=static_cast<Limb>(
                    av<bv || partial<carry_or_borrow);
                output[i]=result;
            }
        }
        std::size_t used=width;
        if (adding && carry_or_borrow!=0u) {
            output[used++]=carry_or_borrow;
        }
        return NewDyadic(Integer::from_magnitude(
            BigUnsigned(BigUnsigned::LimbView(output.data(),used)),negative),
            exponent);
    }

  public:
    friend NewDyadic operator+(NewDyadic const& a,NewDyadic const& b) {
        if (!a.is_finite() || !b.is_finite()) {
            if (a._kind==Kind::NaN || b._kind==Kind::NaN) { return nan(); }
            if (!a.is_finite() && !b.is_finite()) {
                return a._kind==b._kind ? a : nan();
            }
            return a.is_finite() ? b : a;
        }
        const std::size_t e=std::max(a._denominator_exponent,
                                      b._denominator_exponent);
        if (prefer_aligned_add_sub(a,b)) {
            return aligned_add_sub(a,b,false);
        }
        return NewDyadic(shifted(a._coefficient,e-a._denominator_exponent)
                       +shifted(b._coefficient,e-b._denominator_exponent),e);
    }
    friend NewDyadic operator-(NewDyadic const& a,NewDyadic const& b) {
        if (a.is_finite() && b.is_finite()) {
            if (prefer_aligned_add_sub(a,b)) {
                return aligned_add_sub(a,b,true);
            }
        }
        return a+(-b);
    }
    friend NewDyadic operator*(NewDyadic const& a,NewDyadic const& b) {
        if (!a.is_finite() || !b.is_finite()) {
            if (a._kind==Kind::NaN || b._kind==Kind::NaN ||
                    signum(a)==0 || signum(b)==0) { return nan(); }
            return inf(signum(a)==signum(b) ? 1 : -1);
        }
        if (a._denominator_exponent>
                std::numeric_limits<std::size_t>::max()-b._denominator_exponent) {
            throw std::overflow_error("NewDyadic binary exponent overflow");
        }
        return NewDyadic(a._coefficient*b._coefficient,
                         a._denominator_exponent+b._denominator_exponent);
    }
    friend NewDyadic operator/(NewDyadic x,int denominator) {
        if (denominator!=2) {
            throw std::domain_error("only divide by two is a Dyadic operation");
        }
        if (x.is_finite() && !x.is_zero()) {
            if (x._denominator_exponent==std::numeric_limits<std::size_t>::max()) {
                throw std::overflow_error("NewDyadic binary exponent overflow");
            }
            ++x._denominator_exponent;
            x.normalize();
        }
        return x;
    }
    friend NewDyadic hlf(NewDyadic const& x) { return x/2; }

    friend NewDyadic pow(NewDyadic const& x,Nat exponent) {
        NewDyadic result(1);
        NewDyadic factor(x);
        for (unsigned n=exponent;n!=0u;n>>=1u) {
            if ((n&1u)!=0u) { result=result*factor; }
            if (n>1u) { factor=factor*factor; }
        }
        return result;
    }
    friend Boolean eq(NewDyadic const& a,NewDyadic const& b) {
        if (a._kind==Kind::NaN || b._kind==Kind::NaN) { return false; }
        if (a._kind!=b._kind) { return false; }
        if (!a.is_finite()) { return true; }
        return a._denominator_exponent==b._denominator_exponent &&
               static_cast<bool>(a._coefficient==b._coefficient);
    }
    friend Comparison cmp(NewDyadic const& a,NewDyadic const& b) {
        if (a._kind==Kind::NaN || b._kind==Kind::NaN) {
            // cmp does not model unordered; Boolean comparisons
            // explicitly reject NaN below, as with the legacy type.
            return Comparison::EQUAL;
        }
        if (a._kind==b._kind && !a.is_finite()) {
            return Comparison::EQUAL;
        }
        if (a._kind==Kind::NegativeInfinity ||
                b._kind==Kind::PositiveInfinity) {
            return Comparison::LESS;
        }
        if (a._kind==Kind::PositiveInfinity ||
                b._kind==Kind::NegativeInfinity) {
            return Comparison::GREATER;
        }
        const int sa=signum(a),sb=signum(b);
        if (sa!=sb) { return sa<sb ? Comparison::LESS : Comparison::GREATER; }
        if (sa==0) { return Comparison::EQUAL; }
        const int magnitude=compare_scaled(
            a._coefficient.magnitude(),a._denominator_exponent,
            b._coefficient.magnitude(),b._denominator_exponent);
        const int order=sa<0 ? -magnitude : magnitude;
        return order<0 ? Comparison::LESS :
               order>0 ? Comparison::GREATER : Comparison::EQUAL;
    }
    friend Boolean lt(NewDyadic const& a,NewDyadic const& b) {
        if (a._kind==Kind::NaN || b._kind==Kind::NaN) { return false; }
        return cmp(a,b)==Comparison::LESS;
    }
    friend Sign sgn(NewDyadic const& x) {
        return signum(x)>0 ? Sign::POSITIVE :
               signum(x)<0 ? Sign::NEGATIVE : Sign::ZERO;
    }
    friend Bool is_nan(NewDyadic const& x) {
        return x._kind==Kind::NaN;
    }
    friend Bool is_inf(NewDyadic const& x) {
        return x._kind==Kind::PositiveInfinity ||
               x._kind==Kind::NegativeInfinity;
    }
    friend Bool is_finite(NewDyadic const& x) { return x.is_finite(); }
    friend Bool is_zero(NewDyadic const& x) { return x.is_zero(); }
    friend Boolean operator==(NewDyadic const& a,NewDyadic const& b) { return eq(a,b); }
    friend Boolean operator!=(NewDyadic const& a,NewDyadic const& b) { return !eq(a,b); }

    friend Boolean operator<(NewDyadic const& a,NewDyadic const& b) { return lt(a,b); }
    friend Boolean operator>(NewDyadic const& a,NewDyadic const& b) { return lt(b,a); }
    friend Boolean operator<=(NewDyadic const& a,NewDyadic const& b) { return !lt(b,a); }
    friend Boolean operator>=(NewDyadic const& a,NewDyadic const& b) { return !lt(a,b); }
    friend Rational pow(NewDyadic const& x,Int exponent) {
        if (exponent>=0) {
            return pow(x,static_cast<Nat>(exponent)).to_rational();
        }
        const Nat positive=static_cast<Nat>(-(exponent+1))+Nat(1);
        return Rational(1)/pow(x,positive).to_rational();
    }
    friend NewDyadic min(NewDyadic const& a,NewDyadic const& b) {
        if (a._kind==Kind::NaN || b._kind==Kind::NaN) { return nan(); }
        return b<a ? b : a;
    }
    friend NewDyadic max(NewDyadic const& a,NewDyadic const& b) {
        if (a._kind==Kind::NaN || b._kind==Kind::NaN) { return nan(); }
        return a<b ? b : a;
    }
    friend Positive<NewDyadic> abs(NewDyadic const& x);
    friend Positive<NewDyadic> sqr(NewDyadic const& x);
    friend Integer floor(NewDyadic const& x) {
        require_finite(x);
        Integer z=truncated(x);
        if (x._coefficient.negative() &&
                has_fraction(x._coefficient.magnitude(),x._denominator_exponent)) {
            z=z-Integer(1);
        }
        return z;
    }
    friend Integer ceil(NewDyadic const& x) {
        require_finite(x);
        Integer z=truncated(x);
        if (!x._coefficient.negative() &&
                has_fraction(x._coefficient.magnitude(),x._denominator_exponent)) {
            z=z+Integer(1);
        }
        return z;
    }
    friend Integer round(NewDyadic const& x) {
        require_finite(x);
        const auto& magnitude=x._coefficient.magnitude();
        BigUnsigned integer=shift_right(magnitude,x._denominator_exponent);
        // Round ties away from zero, matching current Dyadic.
        if (x._denominator_exponent!=0u) {
            const auto half=shift_right(magnitude,x._denominator_exponent-1u);
            if ((half.limb(0u)&1u)!=0u) {
                integer=add(integer,BigUnsigned(1u));
            }
        }
        return Integer::from_magnitude(std::move(integer),x._coefficient.negative());
    }
};

template<> class Positive<NewDyadic> : public NewDyadic {
public:
    Positive() : NewDyadic() { }
    template<BuiltinUnsignedIntegral M>
    Positive(M n) : NewDyadic(n) { }
    Positive(int) = delete;
    explicit Positive(NewDyadic const& x) : NewDyadic(x) {
        if (sgn(x)==Sign::NEGATIVE) {
            throw std::domain_error("negative Positive<NewDyadic>");
        }
    }
};
inline Positive<NewDyadic> abs(NewDyadic const& x) {
    return Positive<NewDyadic>(sgn(x)==Sign::NEGATIVE ? -x : x);
}
inline Positive<NewDyadic> sqr(NewDyadic const& x) {
    return Positive<NewDyadic>(x*x);
}

} // namespace Ariadne
#endif // ARIADNE_NEW_DYADIC_HPP
