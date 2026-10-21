/***************************************************************************
 *            numeric/dyadic.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/*! \file numeric/dyadic.hpp
 *  \brief Dyadic numbers.
 */

#ifndef ARIADNE_DYADIC_HPP
#define ARIADNE_DYADIC_HPP

#include "numeric/integer.hpp"
#include "numeric/rational.hpp"
#include "numeric/decimal.hpp"
#include "numeric/twoexp.hpp"
#include "numeric/builtin.hpp"
#include "numeric/number.decl.hpp"
#include "utility/writable.hpp"
#include "utility/handle.hpp"
#include "foundation/logical.hpp"
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
class Dyadic;
extern const Dyadic infty;
class InfinityException : public std::runtime_error {
  public: using std::runtime_error::runtime_error;
};
// Exact finite value: coefficient / 2^denominator_exponent.
// Canonical form has no factor of two in the coefficient if exponent>0.
// GMP is deliberately absent from representation and arithmetic.
class Dyadic {
    static Writer<Dyadic> _default_writer;
public:
    typedef ExactTag Paradigm;
    static Void set_default_writer(Writer<Dyadic> writer);
    static Writer<Dyadic> default_writer();
    operator ExactNumber() const;
    friend OutputStream& operator<<(OutputStream&,Dyadic const&);

    enum class Kind : unsigned char { Finite, PositiveInfinity, NegativeInfinity, NaN };
private:
    Integer _coefficient;
    std::size_t _denominator_exponent=0u;
    Kind _kind=Kind::Finite;

    static std::size_t checked_natural_exponent(Natural const& exp) {
        const BigUnsigned& magnitude=exp.magnitude();
        if (magnitude.bit_length()>std::numeric_limits<std::size_t>::digits) {
            throw std::overflow_error("Dyadic denominator exponent too large");
        }
        const std::size_t value=exp.template get<std::size_t>();
        // Match the original Natural overload's signed Int exponent range.
        if (!std::in_range<Int>(value)) {
            throw std::overflow_error("Dyadic denominator exponent too large");
        }
        return value;
    }
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

    static Integer truncated(Dyadic const& x) {
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

    static void require_finite(Dyadic const& x) {
        if (x._kind!=Kind::Finite) {
            throw std::domain_error("non-finite Dyadic cannot be rounded");
        }
    }

    static int signum(Dyadic const& x) noexcept {
        if (x._kind==Kind::PositiveInfinity) { return 1; }
        if (x._kind==Kind::NegativeInfinity) { return -1; }
        if (x._kind==Kind::NaN || x._coefficient.magnitude().is_zero()) { return 0; }
        return x._coefficient.negative() ? -1 : 1;
    }

public:
    Dyadic()=default;
    Dyadic(Integer numerator,std::size_t exponent=0u)
        : _coefficient(std::move(numerator)),_denominator_exponent(exponent) {
        normalize();
    }
    Dyadic(Integer const& numerator,Natural exponent)
        : Dyadic(numerator,checked_natural_exponent(exponent)) { }
    template<BuiltinSignedIntegral N>
    Dyadic(Integer const&,N) = delete;
    template<BuiltinIntegral N>
    Dyadic(N n) : Dyadic(Integer(n)) { }

    // Same value-conversion entry points as Dyadic; all arithmetic
    // remains native, and only exact dyadic Rational values are accepted.
    explicit Dyadic(Rational const& q) {
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
    Dyadic(ExactDouble const& d) : Dyadic(Rational(d)) { }
    explicit Dyadic(double d) : Dyadic(Rational(ExactDouble(d))) { }
    explicit Dyadic(String const& s) : Dyadic(Rational(Decimal(s))) { }
    explicit Dyadic(FloatDP const& d);
    explicit Dyadic(FloatAP const& d);
    Dyadic(TwoExp const& power) {
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
            throw std::overflow_error("Dyadic exponent out of range");
        }
        return static_cast<Int>(_denominator_exponent);
    }
    String literal() const {
        if (_kind==Kind::NaN) { return String("NaN"); }
        if (_kind==Kind::PositiveInfinity) { return String("inf"); }
        if (_kind==Kind::NegativeInfinity) { return String("-inf"); }
        if (_denominator_exponent==0u) { return _coefficient.literal(); }
        if (_denominator_exponent>std::numeric_limits<Nat>::max()) {
            throw std::overflow_error("Dyadic literal exponent out of range");
        }
        return Decimal(_coefficient*pow(Integer(5),
            static_cast<Nat>(_denominator_exponent)),
            static_cast<Nat>(_denominator_exponent)).literal();
    }

    static Dyadic inf() { return inf(1); }
    static Dyadic inf(int sign) {
        Dyadic result;
        result._kind=sign==0 ? Kind::NaN :
                     sign>0 ? Kind::PositiveInfinity : Kind::NegativeInfinity;
        return result;
    }
    static Dyadic nan() { return inf(0); }
    static Dyadic inf(Sign sign) {
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

    friend Dyadic operator-(Dyadic const& x) {
        if (x._kind==Kind::NaN) { return nan(); }
        if (x._kind==Kind::PositiveInfinity) { return inf(-1); }
        if (x._kind==Kind::NegativeInfinity) { return inf(1); }
        return Dyadic(-x._coefficient,x._denominator_exponent);
    }
    friend Dyadic operator+(Dyadic const& x) { return x; }
    friend Dyadic nul(Dyadic const&) { return Dyadic(0); }
    friend Dyadic pos(Dyadic const& x) { return x; }
    friend Dyadic neg(Dyadic const& x) { return -x; }
    friend Dyadic add(Dyadic const& a,Dyadic const& b) {
        return a+b;
    }
    friend Dyadic sub(Dyadic const& a,Dyadic const& b) {
        return a-b;
    }
    friend Dyadic mul(Dyadic const& a,Dyadic const& b) {
        return a*b;
    }
    friend Dyadic fma(Dyadic const& a,Dyadic const& b,
                          Dyadic const& c) { return a*b+c; }
    Dyadic& operator+=(Dyadic const& rhs) {
        *this=*this+rhs; return *this;
    }
    Dyadic& operator-=(Dyadic const& rhs) {
        *this=*this-rhs; return *this;
    }
    Dyadic& operator*=(Dyadic const& rhs) {
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

    static bool prefer_aligned_add_sub(Dyadic const& a,
                                       Dyadic const& b) noexcept {
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

    // R38: use one owned magnitude allocation directly in BigUnsigned.
    // Large inputs bypass both the two-shift fallback and the R35
    // stack-buffer/owning-buffer double materialization. The equality
    // and sign decision use virtual magnitude comparison only.
    static bool prefer_owned_shifted_add_sub(Dyadic const& a,
                                              Dyadic const& b) noexcept {
#if defined(__aarch64__) || defined(_M_ARM64)
        constexpr std::size_t min_limbs=32u;
#else
        constexpr std::size_t min_limbs=17u;
#endif
        const std::size_t largest=std::max(
            a._coefficient.magnitude().size(),
            b._coefficient.magnitude().size());
        const std::size_t high=std::max(
            a._denominator_exponent,b._denominator_exponent);
        const std::size_t low=std::min(
            a._denominator_exponent,b._denominator_exponent);
        return largest>=min_limbs && high-low<=63u;
    }

    static Dyadic owned_shifted_add_sub(
        Dyadic const& a,Dyadic const& b,bool subtract_second) {
        const BigUnsigned& an=a._coefficient.magnitude();
        const BigUnsigned& bn=b._coefficient.magnitude();
        const std::size_t exponent=std::max(
            a._denominator_exponent,b._denominator_exponent);
        const unsigned shift_a=static_cast<unsigned>(
            exponent-a._denominator_exponent);
        const unsigned shift_b=static_cast<unsigned>(
            exponent-b._denominator_exponent);
        const bool negative_a=a._coefficient.negative();
        const bool negative_b=b._coefficient.negative()!=subtract_second;
        const bool same_sign=negative_a==negative_b;

        BigUnsigned magnitude;
        bool result_negative=negative_a;
        if (same_sign) {
            magnitude=shift_a==0u && shift_b==0u
                ? add(an,bn)
                : BigUnsigned::add_shifted(
                    an.limbs(),shift_a,bn.limbs(),shift_b);
        } else {
            const int order=compare_scaled(
                an,a._denominator_exponent,
                bn,b._denominator_exponent);
            if (order==0) { return Dyadic(); }
            if (order>0) {
                magnitude=shift_a==0u && shift_b==0u
                    ? sub(an,bn)
                    : BigUnsigned::sub_shifted(
                        an.limbs(),shift_a,bn.limbs(),shift_b);
            } else {
                magnitude=shift_a==0u && shift_b==0u
                    ? sub(bn,an)
                    : BigUnsigned::sub_shifted(
                        bn.limbs(),shift_b,an.limbs(),shift_a);
                result_negative=negative_b;
            }
        }
        return Dyadic(Integer::from_magnitude(
            std::move(magnitude),result_negative),exponent);
    }

    // R35: for aligned finite signed arithmetic on 5..64-limb
    // magnitudes with an exponent gap <=63, scan borrowed virtual
    // shifted limbs and produce one canonical owned magnitude.
    // No shifted Integer operands or temporary negation.
    static Dyadic aligned_add_sub(Dyadic const& a,Dyadic const& b,
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
            if (order==0) { return Dyadic(); }
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
        return Dyadic(Integer::from_magnitude(
            BigUnsigned(BigUnsigned::LimbView(output.data(),used)),negative),
            exponent);
    }

  public:
    friend Dyadic operator+(Dyadic const& a,Dyadic const& b) {
        if (!a.is_finite() || !b.is_finite()) {
            if (a._kind==Kind::NaN || b._kind==Kind::NaN) { return nan(); }
            if (!a.is_finite() && !b.is_finite()) {
                return a._kind==b._kind ? a : nan();
            }
            return a.is_finite() ? b : a;
        }
        const std::size_t e=std::max(a._denominator_exponent,
                                      b._denominator_exponent);
        if (prefer_owned_shifted_add_sub(a,b)) {
            return owned_shifted_add_sub(a,b,false);
        }
        if (prefer_aligned_add_sub(a,b)) {
            return aligned_add_sub(a,b,false);
        }
        return Dyadic(shifted(a._coefficient,e-a._denominator_exponent)
                       +shifted(b._coefficient,e-b._denominator_exponent),e);
    }
    friend Dyadic operator-(Dyadic const& a,Dyadic const& b) {
        if (a.is_finite() && b.is_finite()) {
            if (prefer_owned_shifted_add_sub(a,b)) {
                return owned_shifted_add_sub(a,b,true);
            }
            if (prefer_aligned_add_sub(a,b)) {
                return aligned_add_sub(a,b,true);
            }
        }
        return a+(-b);
    }
    friend Dyadic operator*(Dyadic const& a,Dyadic const& b) {
        if (!a.is_finite() || !b.is_finite()) {
            if (a._kind==Kind::NaN || b._kind==Kind::NaN ||
                    signum(a)==0 || signum(b)==0) { return nan(); }
            return inf(signum(a)==signum(b) ? 1 : -1);
        }
        if (a._denominator_exponent>
                std::numeric_limits<std::size_t>::max()-b._denominator_exponent) {
            throw std::overflow_error("Dyadic binary exponent overflow");
        }
        return Dyadic(a._coefficient*b._coefficient,
                         a._denominator_exponent+b._denominator_exponent);
    }
    friend Dyadic operator/(Dyadic x,int denominator) {
        if (denominator!=2) {
            throw std::domain_error("only divide by two is a Dyadic operation");
        }
        if (x.is_finite() && !x.is_zero()) {
            if (x._denominator_exponent==std::numeric_limits<std::size_t>::max()) {
                throw std::overflow_error("Dyadic binary exponent overflow");
            }
            ++x._denominator_exponent;
            x.normalize();
        }
        return x;
    }
    friend Dyadic hlf(Dyadic const& x) { return x/2; }
    friend Rational operator/(Dyadic const& a,Dyadic const& b) {
        return div(a.to_rational(),b.to_rational());
    }
    // Quotients and reciprocals of dyadics are exact rational numbers.
    friend Rational rec(Dyadic const& x) {
        return rec(x.to_rational());
    }
    friend Rational div(Dyadic const& a,Dyadic const& b) {
        return div(a.to_rational(),b.to_rational());
    }
    friend Rational div(Dyadic const& a,Integer const& b) {
        return div(a.to_rational(),Rational(b));
    }
    friend Rational div(Integer const& a,Dyadic const& b) {
        return div(Rational(a),b.to_rational());
    }

    friend Dyadic pow(Dyadic const& x,Nat exponent) {
        Dyadic result(1);
        Dyadic factor(x);
        for (unsigned n=exponent;n!=0u;n>>=1u) {
            if ((n&1u)!=0u) { result=result*factor; }
            if (n>1u) { factor=factor*factor; }
        }
        return result;
    }
    friend Boolean eq(Dyadic const& a,Dyadic const& b) {
        if (a._kind==Kind::NaN || b._kind==Kind::NaN) { return false; }
        if (a._kind!=b._kind) { return false; }
        if (!a.is_finite()) { return true; }
        return a._denominator_exponent==b._denominator_exponent &&
               static_cast<bool>(a._coefficient==b._coefficient);
    }
    friend Comparison cmp(Dyadic const& a,Dyadic const& b) {
        if (a._kind==Kind::NaN || b._kind==Kind::NaN) {
            return Comparison::INCOMPARABLE;
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
    friend Comparison cmp(Dyadic const& a,Int b) { return cmp(a,Dyadic(b)); }
    friend Comparison cmp(Int a,Dyadic const& b) { return cmp(Dyadic(a),b); }
    friend Boolean lt(Dyadic const& a,Dyadic const& b) {
        if (a._kind==Kind::NaN || b._kind==Kind::NaN) { return false; }
        return cmp(a,b)==Comparison::LESS;
    }
    friend Sign sgn(Dyadic const& x) {
        return signum(x)>0 ? Sign::POSITIVE :
               signum(x)<0 ? Sign::NEGATIVE : Sign::ZERO;
    }
    friend Bool is_nan(Dyadic const& x) {
        return x._kind==Kind::NaN;
    }
    friend Bool is_inf(Dyadic const& x) {
        return x._kind==Kind::PositiveInfinity ||
               x._kind==Kind::NegativeInfinity;
    }
    friend Bool is_finite(Dyadic const& x) { return x.is_finite(); }
    friend Bool is_zero(Dyadic const& x) { return x.is_zero(); }
    friend Boolean operator==(Dyadic const& a,Dyadic const& b) { return eq(a,b); }
    friend Boolean operator!=(Dyadic const& a,Dyadic const& b) { return !eq(a,b); }

    friend Boolean operator<(Dyadic const& a,Dyadic const& b) { return lt(a,b); }
    friend Boolean operator>(Dyadic const& a,Dyadic const& b) { return lt(b,a); }
    friend Boolean operator<=(Dyadic const& a,Dyadic const& b) { return !lt(b,a); }
    friend Boolean operator>=(Dyadic const& a,Dyadic const& b) { return !lt(a,b); }
    friend Rational pow(Dyadic const& x,Int exponent) {
        if (exponent>=0) {
            return pow(x,static_cast<Nat>(exponent)).to_rational();
        }
        const Nat positive=static_cast<Nat>(-(exponent+1))+Nat(1);
        return Rational(1)/pow(x,positive).to_rational();
    }
    friend Dyadic min(Dyadic const& a,Dyadic const& b) {
        if (a._kind==Kind::NaN || b._kind==Kind::NaN) { return nan(); }
        return b<a ? b : a;
    }
    friend Dyadic max(Dyadic const& a,Dyadic const& b) {
        if (a._kind==Kind::NaN || b._kind==Kind::NaN) { return nan(); }
        return a<b ? b : a;
    }
    friend Positive<Dyadic> abs(Dyadic const& x);
    friend Positive<Dyadic> sqr(Dyadic const& x);
    friend Integer floor(Dyadic const& x) {
        require_finite(x);
        Integer z=truncated(x);
        if (x._coefficient.negative() &&
                x._denominator_exponent!=0u) {
            z=z-Integer(1);
        }
        return z;
    }
    friend Integer ceil(Dyadic const& x) {
        require_finite(x);
        Integer z=truncated(x);
        if (!x._coefficient.negative() &&
                x._denominator_exponent!=0u) {
            z=z+Integer(1);
        }
        return z;
    }
    friend Integer round(Dyadic const& x) {
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

class DecimalWriter : public WriterInterface<Dyadic> {
    virtual OutputStream& _write(OutputStream& os,Dyadic const& x) const final override;
};
class ScientificWriter : public WriterInterface<Dyadic> {
    virtual OutputStream& _write(OutputStream& os,Dyadic const& x) const final override;
};
class FractionWriter : public WriterInterface<Dyadic> {
    virtual OutputStream& _write(OutputStream& os,Dyadic const& x) const final override;
};
template<> class RepresentationWriter<Dyadic> : public WriterInterface<Dyadic> {
    virtual OutputStream& _write(OutputStream& os,Dyadic const& x) const final override;
};
template<> class Positive<Dyadic> : public Dyadic {
public:
    Positive() : Dyadic() { }
    template<BuiltinUnsignedIntegral M>
    Positive(M n) : Dyadic(n) { }
    Positive(int) = delete;
    explicit Positive(Dyadic const& x) : Dyadic(x) {
        if (sgn(x)==Sign::NEGATIVE) {
            throw std::domain_error("negative Positive<Dyadic>");
        }
    }
    explicit Positive(Dyadic&& x) : Dyadic(std::move(x)) {
        if (sgn(*this)==Sign::NEGATIVE) {
            throw std::domain_error("negative Positive<Dyadic>");
        }
    }
};
inline Positive<Dyadic> cast_positive(Dyadic const& x) { return Positive<Dyadic>(x); }
inline Positive<Dyadic> abs(Dyadic const& x) {
    // D42: copy the canonical magnitude exactly once, without forming
    // a negated Dyadic or renormalizing an already canonical value.
    Dyadic result;
    result._coefficient=Integer::from_magnitude(
        x._coefficient.magnitude(),false);
    result._denominator_exponent=x._denominator_exponent;
    result._kind=x._kind==Dyadic::Kind::NegativeInfinity
        ? Dyadic::Kind::PositiveInfinity : x._kind;
    return Positive<Dyadic>(std::move(result));
}
inline Positive<Dyadic> sqr(Dyadic const& x) {
    return Positive<Dyadic>(x*x);
}

// Existing native core implementation ends here. Wrapper specializations
// live in this same header, matching the final numeric/dyadic.hpp layout.
} // namespace Ariadne

#include "numeric/arithmetic.hpp"
#include "utility/handle.hpp"
#include "utility/writable.hpp"
#include "foundation/logical.hpp"

namespace Ariadne {

using PositiveDyadic = Positive<Dyadic>;

//! \relates Dyadic
//! \name Type synonyms
//!@{
//! \ingroup NumericModule
using DyadicBounds = Bounds<Dyadic>; //!< Alias for dyadic bounds on a number.
using DyadicBall = Ball<Dyadic,Dyadic>; //!< Alias for ball about a number with dyadic value and error.
using DyadicUpperBound = UpperBound<Dyadic>; //!< Alias for dyadic upper bound for a number.
using DyadicLowerBound = LowerBound<Dyadic>; //!< Alias for dyadic lower bound for a number.
using DyadicApproximation = Approximation<Dyadic>; //!< Alias for dyadic approximateion to a number.

using PositiveDyadicApproximation = Positive<Approximation<Dyadic>>; //!< <p/>
using PositiveDyadicLowerBound = Positive<LowerBound<Dyadic>>; //!< <p/>
using PositiveDyadicUpperBound = Positive<UpperBound<Dyadic>>; //!< <p/>
using PositiveDyadicBounds = Positive<Bounds<Dyadic>>; //!< <p/>
//!@}

template<> class Bounds<Dyadic> {
    Dyadic _l, _u;
  public:
    typedef ValidatedTag Paradigm;
    Bounds(Dyadic w) : _l(w), _u(w) { }
    template<ConvertibleTo<Dyadic> X> Bounds(X const& x)
        : Bounds(Dyadic(x)) { }
    Bounds(Dyadic l, Dyadic u) : _l(l), _u(u) { }
    template<class X> requires Constructible<Dyadic,X> Bounds(Bounds<X> const& x)
        : Bounds(Dyadic(x.lower_raw()),Dyadic(x.upper_raw())) { }
    operator ValidatedNumber() const;
    Bounds<Dyadic> pm(Dyadic e) { return Bounds<Dyadic>(_l-e,_u+e); }
    Dyadic lower() const { return _l; }
    Dyadic upper() const { return _u; }
    Dyadic lower_raw() const { return _l; }
    Dyadic upper_raw() const { return _u; }
    Bounds<FloatDP> get(DoublePrecision pr) const;
    Bounds<FloatAP> get(ArbitraryPrecision pr) const;
    Bounds<FloatMP> get(MultiplePrecision pr) const;
    friend Bounds<Dyadic> operator+(Bounds<Dyadic> const& w) { return Bounds<Dyadic>(+w._l,w._u); }
    friend Bounds<Dyadic> operator-(Bounds<Dyadic> const& w) { return Bounds<Dyadic>(-w._u,-w._l); }
    friend Bounds<Dyadic> operator+(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        return Bounds<Dyadic>(w1._l+w2._l,w1._u+w2._u); }
    friend Bounds<Dyadic> operator-(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        return Bounds<Dyadic>(w1._l-w2._u,w1._u-w2._l); }
    friend Bounds<Dyadic> operator*(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        // Exact multiplication enclosure: all four vertices, no rounding.
        Dyadic p=w1._l*w2._l, q=w1._l*w2._u;
        Dyadic r=w1._u*w2._l, s=w1._u*w2._u;
        return Bounds<Dyadic>(min(min(p,q),min(r,s)),
                                 max(max(p,q),max(r,s)));
    }
    friend Bounds<Rational> operator/(Bounds<Rational> const& w1, Bounds<Rational> const& w2);
    friend Bounds<Dyadic> add(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        return Bounds<Dyadic>(w1._l+w2._l,w1._u+w2._u); }
    friend Bounds<Dyadic> sub(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        return Bounds<Dyadic>(w1._l-w2._u,w1._u-w2._l); }
    friend Bounds<Dyadic> mul(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) { return w1*w2; }
    friend Bounds<Rational> div(Bounds<Rational> const& w1, Bounds<Rational> const& w2);
    friend Bounds<Dyadic> pow(Bounds<Dyadic> const& w, Nat m) {
        if (m==0u) { return Bounds<Dyadic>(Dyadic(1)); }
        if (m%2u!=0u) {
            return Bounds<Dyadic>(pow(w._l,m),pow(w._u,m));
        }
        if (w._l>=0) {
            return Bounds<Dyadic>(pow(w._l,m),pow(w._u,m));
        }
        if (w._u<=0) {
            return Bounds<Dyadic>(pow(w._u,m),pow(w._l,m));
        }
        return Bounds<Dyadic>(Dyadic(0),
                                 max(pow(w._l,m),pow(w._u,m)));
    }
    friend Bounds<Dyadic> pow(Bounds<Dyadic> const& w, Int m) {
        if (m<0) { throw std::domain_error("negative Dyadic interval power"); }
        return pow(w,static_cast<Nat>(m));
    }
    friend Bounds<Dyadic> nul(Bounds<Dyadic> const& w) {
        return Bounds<Dyadic>(nul(w._l),nul(w._u)); }
    friend Bounds<Dyadic> pos(Bounds<Dyadic> const& w) {
        return Bounds<Dyadic>(pos(w._l),pos(w._u)); }
    friend Bounds<Dyadic> neg(Bounds<Dyadic> const& w) {
        return Bounds<Dyadic>(neg(w._u),neg(w._l)); }
    friend Bounds<Dyadic> sqr(Bounds<Dyadic> const& w) {
        if(w._l>0) { return Bounds<Dyadic>(sqr(w._l),sqr(w._u)); }
        else if(w._u<0) { return Bounds<Dyadic>(sqr(w._u),sqr(w._l)); }
        else { return Bounds<Dyadic>(nul(w._l),max(sqr(w._l),sqr(w._u))); } }
    friend Bounds<Dyadic> hlf(Bounds<Dyadic> const& w) { return Bounds<Dyadic>(hlf(w._l),hlf(w._u)); }
    friend Bounds<Rational> rec(Bounds<Rational> const& w);
    friend Bounds<Dyadic> abs(Bounds<Dyadic> const& w) { return Bounds<Dyadic>(max(min(w._l,-w._u),0),max(-w._l,w._u)); }
    friend Bounds<Dyadic> max(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) { return Bounds<Dyadic>(max(w1._l,w2._l),max(w1._u,w2._u)); }
    friend Bounds<Dyadic> min(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) { return Bounds<Dyadic>(min(w1._l,w2._l),min(w1._u,w2._u)); }

    friend ValidatedKleenean sgn(Bounds<Dyadic> const& w) {
        if (w._l>0) { return true; } else if (w._u<0) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator==(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        if (w1._l>=w2._u && w1._u<=w2._l) { return true; } else if (w1._u< w2._l || w1._l> w2._u) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator!=(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        if (w1._u< w2._l || w1._l> w2._u) { return true; } else if (w1._l>=w2._u && w1._u<=w2._l) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator<=(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        if (w1._u<=w2._l) { return true; } else if (w1._l> w2._u) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator>=(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        if (w1._l>=w2._u) { return true; } else if (w1._u< w2._l) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator< (Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        if (w1._u< w2._l) { return true; } else if (w1._l>=w2._u) { return false; } else { return indeterminate; } }
    friend ValidatedKleenean operator> (Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        if (w1._l> w2._u) { return true; } else if (w1._u<=w2._l) { return false; } else { return indeterminate; } }
    friend Bool consistent(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        return w1._l<=w2._u && w1._u >= w2._l; }
    friend Boolean inconsistent(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        return w1._l>w2._u || w1._u < w2._l; }
    friend Boolean refines(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        return w1._l>=w2._l and w1._u<=w2._u; }
    friend Bounds<Dyadic> refinement(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        return Bounds<Dyadic>(max(w1._l,w2._l),min(w1._u,w2._u)); }
    friend Bounds<Dyadic> coarsening(Bounds<Dyadic> const& w1, Bounds<Dyadic> const& w2) {
        return Bounds<Dyadic>(min(w1._l,w2._l),max(w1._u,w2._u)); }
    friend OutputStream& operator<<(OutputStream& os, Bounds<Dyadic> y) { return os << "[" << y._l << ":" << y._u << "]"; }
};


template<> class Ball<Dyadic,Dyadic> {
    Dyadic _v, _e;
  public:
    Ball(Dyadic w) : _v(w), _e(0) { }
    Ball(Dyadic v, Dyadic e) : _v(v), _e(e) { }
    template<class X, class XE> requires Constructible<Dyadic,X> and Constructible<Dyadic,XE>
    Ball(Ball<X,XE> const& x)
        : Ball(Dyadic(x.value_raw()),Dyadic(x.error_raw())) { }
    explicit Ball(Bounds<Dyadic> const& w) : _v(hlf(w.lower()+w.upper())), _e(hlf(w.upper()-w.lower())) { }
    explicit operator Bounds<Dyadic>() const { return Bounds<Dyadic>(_v-_e,_v+_e); }
    Dyadic value() const { return _v; }
    Dyadic error() const { return _e; }
    Dyadic value_raw() const { return _v; }
    Dyadic error_raw() const { return _e; }
    friend DyadicBall operator+(DyadicBall const& w1, DyadicBall const& w2) { return DyadicBall(w1._v+w2._v,w1._e+w2._e); }
    friend DyadicBall operator-(DyadicBall const& w1, DyadicBall const& w2) { return DyadicBall(w1._v-w2._v,w1._e+w2._e); }
    friend DyadicBall operator*(DyadicBall const& w1, DyadicBall const& w2) { return DyadicBall(w1._v*w2._v,abs(w1._v)*w2._e+w1._e*abs(w2._v)+w1._e*w2._e); }
    friend DyadicBall abs(DyadicBall const& w) {
        if (abs(w._v)>=w._e) { return DyadicBall(abs(w._v),w._e); } else { Dyadic av=hlf(max(w._e-w._v,w._v+w._e)); return DyadicBall(av,av); } }
    friend ValidatedKleenean operator<(DyadicBall const& w1, DyadicBall const& w2) {
        if (w1._v+w1._e<w2._v-w2._e) { return true; } else if (w1._v-w1._e >= w2._v+w2._e) { return false; } else { return indeterminate; } }
    friend Boolean inconsistent(DyadicBall const& w1, DyadicBall const& w2) { return abs(w1._v-w2._v) > w1._e+w2._e; }
    friend Boolean refines(DyadicBall const& w1, DyadicBall const& w2) { return abs(w1._v-w2._v)+w1._e <= w2._e; }
    friend DyadicBall refinement(DyadicBall const& w1, DyadicBall const& w2) { return DyadicBall(refinement(DyadicBounds(w1),DyadicBounds(w2))); }
    friend DyadicBall coarsening(DyadicBall const& w1, DyadicBall const& w2) { return DyadicBall(coarsening(DyadicBounds(w1),DyadicBounds(w2))); }
    friend OutputStream& operator<<(OutputStream& os, DyadicBall y) { return os << "[" << y._v << ":" << y._e << "]"; }
};

template<> class LowerBound<Dyadic> {
    Dyadic _l;
  public:
    LowerBound(Dyadic l) : _l(l) { }
    LowerBound(Bounds<Dyadic> lu) : _l(lu.lower_raw()) { }
    LowerBound(LowerBound<FloatDP> const& x);
    LowerBound(LowerBound<FloatMP> const& x);
    Dyadic raw() const { return _l; }
    LowerBound<FloatDP> get(DoublePrecision pr) const;
    LowerBound<FloatAP> get(ArbitraryPrecision pr) const;
    LowerBound<FloatMP> get(MultiplePrecision pr) const;
    friend OutputStream& operator<<(OutputStream& os, LowerBound<Dyadic> const& y) { return os << y._l << ":"; }
};

template<> class UpperBound<Dyadic> {
    Dyadic _u;
  public:
    explicit UpperBound(Dyadic u) : _u(u) { }
    UpperBound(Bounds<Dyadic> lu) : _u(lu.upper_raw()) { }
    UpperBound(UpperBound<FloatDP> const& x);
    UpperBound(UpperBound<FloatMP> const& x);
    Dyadic raw() const { return _u; }
    UpperBound<FloatDP> get(DoublePrecision pr) const;
    UpperBound<FloatAP> get(ArbitraryPrecision pr) const;
    UpperBound<FloatMP> get(MultiplePrecision pr) const;
    friend OutputStream& operator<<(OutputStream& os, UpperBound<Dyadic> const& y) { return os << ":" << y._u; }
};

template<> class Approximation<Dyadic> {
    Dyadic _a;
  public:
    Approximation(Dyadic a) : _a(a) { }
    Approximation(Approximation<FloatDP> const& x);
    Dyadic raw() const { return _a; }
    Approximation<FloatDP> get(DoublePrecision pr) const;
    Approximation<FloatAP> get(ArbitraryPrecision pr) const;
    friend OutputStream& operator<<(OutputStream& os, Approximation<Dyadic> const& y) {
        return os << y._a;
    }
};

template<> class Positive<Bounds<Dyadic>> : public Bounds<Dyadic> { public: Positive(Bounds<Dyadic> w) : Bounds<Dyadic>(w) { }; };
template<> class Positive<LowerBound<Dyadic>> : public LowerBound<Dyadic> { public: Positive(LowerBound<Dyadic> w) : LowerBound<Dyadic>(w) { }; };
template<> class Positive<UpperBound<Dyadic>> : public UpperBound<Dyadic> { public: Positive(UpperBound<Dyadic> w) : UpperBound<Dyadic>(w) { } ; };
template<> class Positive<Approximation<Dyadic>> : public Approximation<Dyadic> { public: Positive(Approximation<Dyadic> w) : Approximation<Dyadic>(w) { } };



inline Dyadic operator""_dyadic(long double x) { return Dyadic(static_cast<double>(x)); }
inline Dyadic operator""_dy(long double x) { return operator""_dyadic(x); }
inline Dyadic operator""_q2(long double x) { return operator""_dyadic(x); }
inline Dyadic operator""_bin(long double x) { return operator""_dyadic(x); }
inline Dyadic dy_(long double x) { return operator""_dyadic(x); }
Dyadic make_dyadic(unsigned long long int);
Dyadic hlf(Integer const& n);
// Exact mixed dyadic/rational enclosure operations. Overloads prevent
// ambiguous implicit promotions between distinct interval specializations.
inline RationalBounds rational_enclosure(DyadicBounds const& x) {
    return RationalBounds(Rational(x.lower_raw()),Rational(x.upper_raw()));
}
inline RationalBounds div(DyadicBounds const& a,DyadicBounds const& b) {
    return div(rational_enclosure(a),rational_enclosure(b));
}
#define ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY(NAME) \
inline RationalBounds NAME(RationalBounds const& q,DyadicBounds const& d) { \
    return NAME(q,rational_enclosure(d)); \
} \
inline RationalBounds NAME(DyadicBounds const& d,RationalBounds const& q) { \
    return NAME(rational_enclosure(d),q); \
}
ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY(add)
ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY(sub)
ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY(mul)
ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY(div)
ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY(max)
ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY(min)
#undef ARIADNE_DYADIC_RATIONAL_BOUNDS_BINARY

#define ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE(OP) \
inline ValidatedKleenean operator OP(RationalBounds const& q,DyadicBounds const& d) { \
    return q OP rational_enclosure(d); \
} \
inline ValidatedKleenean operator OP(DyadicBounds const& d,RationalBounds const& q) { \
    return rational_enclosure(d) OP q; \
}
ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE(==)
ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE(!=)
ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE(<)
ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE(>)
ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE(<=)
ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE(>=)
#undef ARIADNE_DYADIC_RATIONAL_BOUNDS_COMPARE

} // namespace Ariadne
#endif // ARIADNE_DYADIC_HPP
