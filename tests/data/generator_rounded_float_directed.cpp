// SPDX-License-Identifier: MIT
// Reproduce the 128-bit directed exp/log/atan oracle using rigorous integer intervals.
// Standalone standard C++20; never uses Numeric, GMP, MPFR, FloatMP or host floats.
#include "oracle_exact_integer.hpp"
#include <bit>
#include <cstdint>
#include <limits>
#include <tuple>

namespace {
using independent_oracle::Nat;
using independent_oracle::row;

std::size_t bit_length(Nat const& x) {
    return x.zero()?0u:(x.w.size()-1u)*64u+64u-
        static_cast<std::size_t>(std::countl_zero(x.w.back()));
}
Nat multiply(Nat const& a,Nat const& b) {
    if(a.zero()||b.zero()) return {};
    // Portable base-2^32 grade-school multiplication: each accumulator
    // fits exactly in uint64_t, unlike a direct 64-bit limb product.
    const auto unpack=[](Nat const& n) {
        std::vector<std::uint32_t> out;
        out.reserve(2u*n.w.size());
        for(auto word:n.w) {
            out.push_back(static_cast<std::uint32_t>(word));
            out.push_back(static_cast<std::uint32_t>(word>>32u));
        }
        while(!out.empty()&&out.back()==0u) out.pop_back();
        return out;
    };
    const auto x=unpack(a),y=unpack(b);
    std::vector<std::uint32_t> p(x.size()+y.size(),0u);
    for(std::size_t i=0;i<x.size();++i) {
        std::uint64_t carry=0u;
        for(std::size_t j=0;j<y.size();++j) {
            const std::uint64_t value=std::uint64_t{x[i]}*y[j]+p[i+j]+carry;
            p[i+j]=static_cast<std::uint32_t>(value);
            carry=value>>32u;
        }
        std::size_t k=i+y.size();
        while(carry) {
            if(k==p.size()) p.push_back(0u);
            const std::uint64_t value=std::uint64_t{p[k]}+carry;
            p[k]=static_cast<std::uint32_t>(value);
            carry=value>>32u;
            ++k;
        }
    }
    std::vector<std::uint64_t> limbs((p.size()+1u)/2u,0u);
    for(std::size_t i=0;i<p.size();++i)
        limbs[i/2u] |= std::uint64_t{p[i]}<<((i%2u)*32u);
    return Nat(std::move(limbs));
}
std::pair<Nat,Nat> divide(Nat remainder,Nat const& divisor) {
    if(divisor.zero()) throw std::invalid_argument("divide by zero");
    Nat quotient;
    while(remainder.compare(divisor)>=0) {
        std::size_t shift=bit_length(remainder)-bit_length(divisor);
        Nat trial=divisor.shl(shift);
        if(trial.compare(remainder)>0) trial=divisor.shl(--shift);
        remainder=Nat::sub(remainder,trial);
        quotient=Nat::add(quotient,Nat::power2(shift));
    }
    return {std::move(quotient),std::move(remainder)};
}
struct Signed {
    Nat magnitude;
    bool negative=false;
    Signed()=default;
    Signed(Nat n,bool sign=false):magnitude(std::move(n)),negative(sign&&!magnitude.zero()){}
};
Signed neg(Signed x) {if(!x.magnitude.zero())x.negative=!x.negative;return x;}
int compare(Signed const& a,Signed const& b) {
    if(a.negative!=b.negative)return a.negative?-1:1;
    const int c=a.magnitude.compare(b.magnitude);
    return a.negative?-c:c;
}
Signed add(Signed const& a,Signed const& b) {
    if(a.negative==b.negative)return {Nat::add(a.magnitude,b.magnitude),a.negative};
    const int c=a.magnitude.compare(b.magnitude);
    if(c==0)return {};
    return c>0?Signed{Nat::sub(a.magnitude,b.magnitude),a.negative}
              :Signed{Nat::sub(b.magnitude,a.magnitude),b.negative};
}
Signed sub(Signed const& a,Signed const& b) {return add(a,neg(b));}
Signed mul(Signed const& a,Signed const& b) {
    return {multiply(a.magnitude,b.magnitude),a.negative!=b.negative};
}
Signed times(Signed const& a,std::uint64_t factor) {
    return {a.magnitude.mul_small(factor),a.negative};
}
Signed shifted(Signed const& a,std::size_t bits) {
    return {a.magnitude.shl(bits),a.negative};
}
Signed floor_divide(Signed const& a,Nat const& positive) {
    auto [q,r]=divide(a.magnitude,positive);
    if(a.negative&&!r.zero())q=Nat::add(q,Nat(1u));
    return {std::move(q),a.negative};
}
Signed ceil_divide(Signed const& a,Nat const& positive) {
    auto [q,r]=divide(a.magnitude,positive);
    if(!a.negative&&!r.zero())q=Nat::add(q,Nat(1u));
    return {std::move(q),a.negative};
}
struct Interval {Signed lo,hi;};
constexpr std::size_t scale_bits=1400u;
Nat const& scale() {static const Nat b=Nat::power2(scale_bits);return b;}
Signed exact(std::uint64_t value) {return {Nat(value)};}
Interval number(std::uint64_t n) {
    const auto fixed=Signed{Nat(n).shl(scale_bits)};
    return {fixed,fixed};
}
Interval fraction(Signed const& numerator,Nat const& denominator) {
    const Signed scaled=shifted(numerator,scale_bits);
    return {floor_divide(scaled,denominator),ceil_divide(scaled,denominator)};
}
Interval ratio(std::uint64_t num,std::uint64_t den) {
    return fraction(exact(num),Nat(den));
}
Interval add(Interval const& a,Interval const& b) {
    return {add(a.lo,b.lo),add(a.hi,b.hi)};
}
Interval neg(Interval const& a) {return {neg(a.hi),neg(a.lo)};}
Interval mul(Interval const& a,Interval const& b) {
    const std::array<Signed,4u> products={mul(a.lo,b.lo),mul(a.lo,b.hi),
        mul(a.hi,b.lo),mul(a.hi,b.hi)};
    Signed low=products[0],high=products[0];
    for(auto const& p:products) {
        if(compare(p,low)<0)low=p;
        if(compare(p,high)>0)high=p;
    }
    return {floor_divide(low,scale()),ceil_divide(high,scale())};
}
Interval square(Interval const& a) {
    if(a.lo.negative)return mul(a,a);
    return {floor_divide(mul(a.lo,a.lo),scale()),
            ceil_divide(mul(a.hi,a.hi),scale())};
}
Interval divide_positive(Interval const& a,Interval const& b) {
    if(b.lo.negative||b.lo.magnitude.zero())
        throw std::invalid_argument("nonpositive interval denominator");
    const std::array<Signed,4u> n={shifted(a.lo,scale_bits),shifted(a.lo,scale_bits),
                                  shifted(a.hi,scale_bits),shifted(a.hi,scale_bits)};
    const std::array<Nat,4u> d={b.lo.magnitude,b.hi.magnitude,
                                b.lo.magnitude,b.hi.magnitude};
    Signed low=floor_divide(n[0],d[0]),high=ceil_divide(n[0],d[0]);
    for(std::size_t i=1;i<n.size();++i) {
        const auto l=floor_divide(n[i],d[i]);
        const auto h=ceil_divide(n[i],d[i]);
        if(compare(l,low)<0)low=l;
        if(compare(h,high)>0)high=h;
    }
    return {std::move(low),std::move(high)};
}
Interval scale(Interval const& a,std::uint64_t factor) {
    return {times(a.lo,factor),times(a.hi,factor)};
}
Interval power(Interval base,unsigned exponent) {
    Interval result=number(1u);
    for(;exponent!=0u;exponent>>=1u) {
        if(exponent&1u)result=mul(result,base);
        if(exponent>1u)base=square(base);
    }
    return result;
}
Interval exponential_abs(std::uint64_t n,std::uint64_t d) {
    const Interval z=ratio(n,16u*d);
    Interval term=number(1u),sum=term;
    for(unsigned k=1u;k<=241u;++k) {
        term=divide_positive(mul(term,z),number(k));
        sum=add(sum,term);
    }
    const Interval next=divide_positive(mul(term,z),number(242u));
    sum.hi=add(sum.hi,times(next.hi,2u));
    for(unsigned j=0u;j<4u;++j)sum=square(sum);
    return sum;
}
Interval logarithm_positive(std::uint64_t n,std::uint64_t d) {
    const std::uint64_t diff=n>d?n-d:d-n;
    const Interval y=ratio(diff,n+d),y2=square(y);
    Interval term=y,sum=number(0u);
    for(unsigned k=0u;k<271u;++k) {
        sum=add(sum,divide_positive(term,number(2u*k+1u)));
        term=mul(term,y2);
    }
    sum=scale(sum,2u);
    const Interval tail=power(y,543u);
    sum.hi=add(sum.hi,times(tail.hi,4u));
    return n>d?sum:neg(sum);
}
Interval arctangent_unit(std::uint64_t n,std::uint64_t d) {
    const Interval z=ratio(n,d),z2=square(z);
    Interval term=z,sum=number(0u);
    for(unsigned k=0u;k<271u;++k) {
        const auto part=divide_positive(term,number(2u*k+1u));
        sum=add(sum,(k%2u)==0u?part:neg(part));
        term=mul(term,z2);
    }
    // Alternating-series remainder is enclosed conservatively by
    // +/- the first omitted unsigned power (without its divisor).
    return {sub(sum.lo,term.hi),add(sum.hi,term.hi)};
}
Interval arctangent_signed(std::int64_t n,std::uint64_t d) {
    const auto abs=static_cast<std::uint64_t>(n<0?-n:n);
    const Interval pi4=add(scale(arctangent_unit(1u,5u),4u),
                           neg(arctangent_unit(1u,239u)));
    const std::uint64_t diff=abs>=d?abs-d:d-abs;
    const Interval change=arctangent_unit(diff,abs+d);
    const Interval result=add(pi4,abs>=d?change:neg(change));
    return n<0?neg(result):result;
}
std::pair<std::string,std::string> dyadic_fraction(Signed value,std::int64_t exponent) {
    if(value.magnitude.zero())return {"0","1"};
    if(exponent<0) {
        const auto canceled=std::min(static_cast<std::size_t>(-exponent),
                                     value.magnitude.trailing_zero_bits());
        value.magnitude=value.magnitude.shr(canceled);
        exponent+=static_cast<std::int64_t>(canceled);
    }
    Nat numerator=value.magnitude;
    if(exponent>=0)numerator=numerator.shl(static_cast<std::size_t>(exponent));
    const std::string n=(value.negative?"-":"")+numerator.dec();
    const std::string d=exponent<0?Nat::power2(static_cast<std::size_t>(-exponent)).dec():"1";
    return {n,d};
}
std::pair<std::string,std::string> round_binary(Signed const& input,bool upward) {
    if(input.magnitude.zero())return {"0","1"};
    const auto bits=bit_length(input.magnitude);
    const std::int64_t exponent=static_cast<std::int64_t>(bits)-
                                 static_cast<std::int64_t>(scale_bits)-128;
    Signed quantized;
    if(bits>=128u) {
        const Nat divisor=Nat::power2(bits-128u);
        quantized=upward?ceil_divide(input,divisor):floor_divide(input,divisor);
    } else quantized=shifted(input,128u-bits);
    return dyadic_fraction(std::move(quantized),exponent);
}
std::array<std::string,4u> certify(Interval const& interval) {
    auto lower=round_binary(interval.lo,false);
    auto upper=round_binary(interval.hi,true);
    if(lower!=round_binary(interval.hi,false)||upper!=round_binary(interval.lo,true))
        throw std::runtime_error("128-bit directed roundings not certified");
    return {lower.first,lower.second,upper.first,upper.second};
}
std::string generate() {
    std::string out=
        "# Certified 128-bit directed exp/log/atan oracle, replacing runtime MPFR references.\n"
        "# kind original_num original_den down_num down_den up_num up_den\n"
        "# Preserves the 24 exact original inputs and 48 strict comparisons of test_rounded_float.\n"
        "# Generated by integer-only interval arithmetic using B=2^1400 fixed point.\n"
        "# exp: 241 Taylor terms on |x|/16 + 2*next term; four outward squarings;\n"
        "# negative inputs use outward reciprocal. Bound exp(1/8)<2 establishes tail.\n"
        "# log: 271 positive atanh terms in y=|x-1|/(x+1), plus 4*y^543 tail.\n"
        "# atan: 271 alternating Taylor terms + next term bound; reduce near +/-1\n"
        "# via pi/4 + atan((abs(x)-1)/(abs(x)+1)), where pi/4 is certified\n"
        "# from 4*atan(1/5)-atan(1/239) using the same integer intervals.\n"
        "# The 128-bit DOWN and UP rounding of both interval endpoints must coincide.\n"
        "# No FloatAP, FloatMP, MPFR, GMP or floating-point transcendentals as oracle.\n"
        "# FloatDP output is checked strictly > 128-bit upper or < 128-bit lower,\n"
        "# matching the original tests; NOT for equality at 53-bit precision.\n";
    const auto append=[&](char const* kind,std::int64_t n,std::uint64_t d,Interval const& value) {
        const auto expected=certify(value);
        out+=row({kind,std::to_string(n),std::to_string(d),expected[0],expected[1],
                  expected[2],expected[3]});
    };
    const std::array<std::pair<std::int64_t,std::uint64_t>,12u> exp_sources={{
        {-2,1},{-3,2},{-1,1},{-3,4},{-1,2},{-3,8},
        {3,8},{1,2},{3,4},{1,1},{3,2},{2,1}}};
    for(auto [n,d]:exp_sources) {
        auto value=exponential_abs(static_cast<std::uint64_t>(n<0?-n:n),d);
        if(n<0)value=divide_positive(number(1u),value);
        append("exp",n,d,value);
    }
    const std::array<std::pair<std::int64_t,std::uint64_t>,6u> log_sources={{
        {3,1},{2,1},{3,4},{5,8},{1,2},{3,8}}};
    for(auto [n,d]:log_sources)
        append("log",n,d,logarithm_positive(static_cast<std::uint64_t>(n),d));
    constexpr std::int64_t around=281474976710656LL;
    const std::array<std::pair<std::int64_t,std::uint64_t>,6u> atan_sources={{
        {around+1,static_cast<std::uint64_t>(around)},
        {1,1},{around-1,static_cast<std::uint64_t>(around)},
        {1-around,static_cast<std::uint64_t>(around)},
        {-1,1},{-around-1,static_cast<std::uint64_t>(around)}}};
    for(auto [n,d]:atan_sources)append("atan",n,d,arctangent_signed(n,d));
    return out;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"rounded_float_directed.tsv",generate(),write);
        return 0;
    } catch(std::exception const& e) {
        std::cerr<<"RoundedFloat directed generator: "<<e.what()<<'\n';
        return 1;
    }
}
