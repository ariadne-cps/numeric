// SPDX-License-Identifier: MIT
// Independently certified FloatMP arctangent fixtures, standard C++20 integers only.
// Standalone standard C++20; never uses Numeric, GMP, independent exact arithmetic, independent exact arithmetic or host floats.
#include "oracle_exact_integer.hpp"
#include <bit>
#include <cstdint>
#include <limits>
#include <tuple>
#include <set>
#include <map>

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
constexpr std::size_t scale_bits=1800u;
Nat const& scale() {static const Nat b=Nat::power2(scale_bits);return b;}
Interval number(std::uint64_t n) {
    const auto fixed=Signed{Nat(n).shl(scale_bits)};
    return {fixed,fixed};
}
Interval fraction(Signed const& numerator,Nat const& denominator) {
    const Signed scaled=shifted(numerator,scale_bits);
    return {floor_divide(scaled,denominator),ceil_divide(scaled,denominator)};
}
Interval add(Interval const& a,Interval const& b) {
    return {add(a.lo,b.lo),add(a.hi,b.hi)};
}
Interval neg(Interval const& a) {return {neg(a.hi),neg(a.lo)};}
Interval scale(Interval const& a,std::uint64_t factor) {return {times(a.lo,factor),times(a.hi,factor)};}
[[maybe_unused]] Interval mul(Interval const& a,Interval const& b) {
    const std::array<Signed,4u> products={mul(a.lo,b.lo),mul(a.lo,b.hi),
        mul(a.hi,b.lo),mul(a.hi,b.hi)};
    Signed low=products[0],high=products[0];
    for(auto const& p:products) {
        if(compare(p,low)<0)low=p;
        if(compare(p,high)>0)high=p;
    }
    return {floor_divide(low,scale()),ceil_divide(high,scale())};
}
// Positive fixed-point operations avoid general interval division.
struct Positive {Nat lo,hi;};
std::pair<Nat,std::uint32_t> divide_small(Nat n,std::uint32_t d) {
    if(d==0u)throw std::invalid_argument("zero small divisor");
    std::uint64_t rem=0;
    for(std::size_t i=n.w.size();i>0;--i) {
        const std::uint64_t word=n.w[i-1u];
        const std::uint64_t top=(rem<<32u)|(word>>32u);
        const std::uint64_t qt=top/d;
        rem=top%d;
        const std::uint64_t bot=(rem<<32u)|(word&0xffffffffULL);
        const std::uint64_t qb=bot/d;
        rem=bot%d;
        n.w[i-1u]=(qt<<32u)|qb;
    }
    n.trim();
    return {std::move(n),static_cast<std::uint32_t>(rem)};
}
Positive mul_positive(Positive const& a,Positive const& b) {
    const Nat low=multiply(a.lo,b.lo);
    const Nat high=multiply(a.hi,b.hi);
    Nat lower=low.shr(scale_bits),upper=high.shr(scale_bits);
    if(high.compare(upper.shl(scale_bits))>0)upper=Nat::add(upper,Nat(1u));
    return {std::move(lower),std::move(upper)};
}
Positive divide_positive_small(Positive const& p,std::uint32_t d) {
    auto [lo,rlo]=divide_small(p.lo,d);
    auto [hi,rhi]=divide_small(p.hi,d);
    (void)rlo;
    if(rhi)hi=Nat::add(hi,Nat(1u));
    return {std::move(lo),std::move(hi)};
}
// A rigorous alternating Taylor enclosure for 0<=x<=1/2.
Interval arctangent_unit(Nat const& n,Nat const& d) {
    if(n.zero())return number(0u);
    if(n.compare(d)>0)throw std::invalid_argument("atan Taylor argument exceeds one");
    const auto z=fraction(Signed(n),d);
    const Positive input{z.lo.magnitude,z.hi.magnitude};
    const Positive z2=mul_positive(input,input);
    Positive term=input;
    Interval sum=number(0u);
    const Nat threshold=Nat::power2(180u);
    for(unsigned k=0u;k<1150u;++k) {
        const Positive part=divide_positive_small(term,2u*k+1u);
        const Interval signed_part{Signed(part.lo),Signed(part.hi)};
        sum=add(sum,(k%2u)==0u?signed_part:neg(signed_part));
        term=mul_positive(term,z2);
        if(term.hi.compare(threshold)<=0)break;
    }
    // Alternating-series bound: the first omitted term has absolute
    // value <= term.hi * 2^-1800; a wider symmetric interval is safe.
    const Signed tail(term.hi);
    return {sub(sum.lo,tail),add(sum.hi,tail)};
}
Interval arctangent_positive(Nat const& n,Nat const& d) {
    if(n.zero())return number(0u);
    static const Interval pi4=add(scale(arctangent_unit(Nat(1u),Nat(5u)),4u),
                                 neg(arctangent_unit(Nat(1u),Nat(239u))));
    const int comparison=n.compare(d);
    if(comparison>=0) {
        if(comparison==0)return pi4;
        // atan(x)=pi/2-atan(1/x), x>1 (the original source x>=2).
        return add(scale(pi4,2u),neg(arctangent_unit(d,n)));
    }
    if(n.shl(1u).compare(d)<=0)return arctangent_unit(n,d);
    // atan(x)=pi/4-atan((1-x)/(1+x)), 1/2<x<1.
    return add(pi4,neg(arctangent_unit(Nat::sub(d,n),Nat::add(n,d))));
}
std::pair<std::string,std::string> dyadic_fraction(Signed value,std::int64_t exponent) {
    if(value.magnitude.zero())return {"0","1"};
    if(exponent<0) {
        const auto cut=std::min(static_cast<std::size_t>(-exponent),
                                value.magnitude.trailing_zero_bits());
        value.magnitude=value.magnitude.shr(cut);
        exponent+=static_cast<std::int64_t>(cut);
    }
    Nat numerator=value.magnitude;
    if(exponent>=0)numerator=numerator.shl(static_cast<std::size_t>(exponent));
    return {(value.negative?"-":"")+numerator.dec(),
            exponent<0?Nat::power2(static_cast<std::size_t>(-exponent)).dec():"1"};
}
std::pair<std::string,std::string> round_binary(Signed const& input,unsigned p,char mode) {
    if(input.magnitude.zero())return {"0","1"};
    const auto nbits=bit_length(input.magnitude);
    const auto grid=static_cast<std::int64_t>(nbits)-
        static_cast<std::int64_t>(scale_bits)-static_cast<std::int64_t>(p);
    Signed rounded;
    if(nbits>p) {
        const auto shift=nbits-p;
        const Nat denominator=Nat::power2(shift);
        auto [q,r]=divide(input.magnitude,denominator);
        bool inc=false;
        if(!r.zero()) {
            if(mode=='N') {
                const int c=r.shl(1u).compare(denominator);
                inc=c>0||(c==0&&!q.zero()&&(q.w.front()&1u)!=0u);
            } else if(mode=='D')inc=input.negative;
            else if(mode=='U')inc=!input.negative;
            else if(mode!='Z')throw std::invalid_argument("round mode");
        }
        if(inc)q=Nat::add(q,Nat(1u));
        rounded=Signed(std::move(q),input.negative);
    } else rounded=shifted(input,p-nbits);
    return dyadic_fraction(std::move(rounded),grid);
}
std::pair<std::string,std::string> certify(Interval const& interval,unsigned p,char mode) {
    const auto first=round_binary(interval.lo,p,mode);
    const auto second=round_binary(interval.hi,p,mode);
    if(first!=second)throw std::runtime_error("atan interval not certified at bits="+
                                               std::to_string(p)+" mode="+mode);
    return first;
}
struct Rational {Nat n,d;};
Rational reduce(Rational a) {
    if(a.d.zero())throw std::invalid_argument("rational denominator zero");
    if(a.n.zero())return {Nat(),Nat(1u)};
    auto gcd_unsigned=[](Nat x,Nat y) {
        if(x.zero())return y;
        if(y.zero())return x;
        const auto common=std::min(x.trailing_zero_bits(),y.trailing_zero_bits());
        x=x.shr(x.trailing_zero_bits());
        y=y.shr(y.trailing_zero_bits());
        while(x.compare(y)!=0) {
            if(x.compare(y)>0) {x=Nat::sub(x,y);x=x.shr(x.trailing_zero_bits());}
            else {y=Nat::sub(y,x);y=y.shr(y.trailing_zero_bits());}
        }
        return x.shl(common);
    };
    const Nat gcd=gcd_unsigned(a.n,a.d);
    a.n=divide(a.n,gcd).first;
    a.d=divide(a.d,gcd).first;
    return a;
}
Rational rational(std::uint64_t n,std::uint64_t d=1u) {
    return reduce({Nat(n),Nat(d)});
}
Rational power_ratio(std::uint64_t n,unsigned shift) {
    return reduce({Nat(n),Nat::power2(shift)});
}
Rational plus(Rational const& a,Rational const& b) {
    return reduce({Nat::add(multiply(a.n,b.d),multiply(b.n,a.d)),
                   multiply(a.d,b.d)});
}
Rational minus(Rational const& a,Rational const& b) {
    const Nat first=multiply(a.n,b.d),second=multiply(b.n,a.d);
    if(first.compare(second)<0)throw std::invalid_argument("negative magnitude source");
    return reduce({Nat::sub(first,second),multiply(a.d,b.d)});
}
std::int64_t binary_exponent(Nat const& n,Nat const& d) {
    std::int64_t e=static_cast<std::int64_t>(bit_length(n))-
                   static_cast<std::int64_t>(bit_length(d));
    if(e>=0) {
        if(n.compare(d.shl(static_cast<std::size_t>(e)))<0)--e;
    } else if(n.shl(static_cast<std::size_t>(-e)).compare(d)<0)--e;
    return e;
}
Rational nearest(Rational const& a,unsigned p) {
    if(a.n.zero())return rational(0u);
    const auto grid=binary_exponent(a.n,a.d)-static_cast<std::int64_t>(p)+1;
    const Nat numerator=grid<0?a.n.shl(static_cast<std::size_t>(-grid)):a.n;
    const Nat denominator=grid>=0?a.d.shl(static_cast<std::size_t>(grid)):a.d;
    auto [q,r]=divide(numerator,denominator);
    const int c=r.shl(1u).compare(denominator);
    if(c>0||(c==0&&!q.zero()&&(q.w.front()&1u)!=0u))q=Nat::add(q,Nat(1u));
    if(grid<0)return reduce({std::move(q),Nat::power2(static_cast<std::size_t>(-grid))});
    return reduce({q.shl(static_cast<std::size_t>(grid)),Nat(1u)});
}
struct Generator {
    std::string out=
        "# Certified correctly rounded FloatMP atan references (D87).\n"
        "# bits mode original_input_num original_input_den expected_num expected_den\n"
        "# Original Rational is the lookup key; nearest/even input rounding is independently exact.\n"
        "# B=2^1800; outward alternating atan Taylor bounds k=0..1149 and next-term tail.\n"
        "# For 0<=x<=1/2: atan(x) direct; for 1/2<x<=1: pi/4-atan((1-x)/(1+x)).\n"
        "# For x>1: pi/2-atan(1/x). Handle negative inputs by oddness.\n"
        "# pi/4=4*atan(1/5)-atan(1/239) certified with the same integer intervals.\n"
        "# Emit rows ONLY when both exact interval endpoints round to same binary rational.\n"
        "# No GMP, independent exact arithmetic, independent exact arithmetic, FloatMP or floating-point transcendental calculations.\n";
    std::set<std::tuple<unsigned,char,std::string,std::string>> seen;
    std::map<std::pair<std::string,std::string>,Interval> cache;
    std::size_t count=0;
    void emit(unsigned p,Rational const& value,bool negative) {
        const std::string num=(negative&&!value.n.zero()?"-":"")+value.n.dec();
        const std::string den=value.d.dec();
        bool needed=false;
        for(char mode:{'N','D','U','Z'})
            if(!seen.contains({p,mode,num,den}))needed=true;
        if(!needed)return;
        const auto rounded=nearest(value,p);
        const auto key=std::make_pair(rounded.n.dec(),rounded.d.dec());
        auto it=cache.find(key);
        if(it==cache.end())it=cache.emplace(key,arctangent_positive(rounded.n,rounded.d)).first;
        const auto bounds=negative?neg(it->second):it->second;
        for(char mode:{'N','D','U','Z'}) {
            if(!seen.emplace(p,mode,num,den).second)continue;
            const auto output=certify(bounds,p,mode);
            out+=row({std::to_string(p),std::string(1u,mode),num,den,
                      output.first,output.second});
            ++count;
        }
    }
    void both(unsigned p,Rational const& value) {
        emit(p,value,true);
        emit(p,value,false);
    }
    void generate() {
        // D87: original centered atan sources, revisiting precision 53.
        for(unsigned p:{53u,7u,128u,256u,1024u,53u}) {
            const auto half=rational(1u,2u);
            const auto ulp=power_ratio(1u,p);
            const auto cutoff=power_ratio(1u,(p+33u)/2u);
            const std::array<Rational,16u> magnitude={
                half,minus(half,ulp),plus(half,ulp),
                minus(half,rational_shift(ulp,4u)),plus(half,rational_shift(ulp,4u)),
                minus(half,cutoff),plus(half,cutoff),
                minus(half,rational_shift(cutoff,2u)),plus(half,rational_shift(cutoff,2u)),
                rational(1u,4u),rational(3u,4u),rational(5u,7u),
                rational(1u,3u),rational(2u,3u),rational(1u),rational(2u)};
            for(auto const& v:magnitude)both(p,v);
        }
        if(count!=640u)throw std::runtime_error("D87 count "+std::to_string(count));
        out+="# D88: additional independently certified single_recurrence_atan probes.\n";
        // D88: exactly the original Dyadic sources; no oracle data literals.
        const std::array<Rational,8u> d88={
            power_ratio(1u,64u),power_ratio(255u,9u),power_ratio(1u,1u),
            power_ratio(257u,9u),rational(1u),power_ratio(1023u,9u),
            rational(2u),power_ratio(1025u,9u)};
        for(unsigned p:{7u,53u,128u,256u,1024u})
            for(auto const& v:d88)both(p,v);
        if(count!=840u)throw std::runtime_error("D88 count "+std::to_string(count));
        out+="# D89: certified atan cache and centered argument-reduction regression keys.\n";
        const auto eps=power_ratio(1u,40u),tiny=power_ratio(1u,100u);
        const std::array<Rational,20u> d89={
            rational(0u),rational(1u,8u),
            minus(rational(1u,4u),eps),rational(1u,4u),plus(rational(1u,4u),eps),
            minus(rational(1u,3u),eps),rational(1u,3u),plus(rational(1u,3u),eps),
            minus(rational(1u,2u),eps),minus(rational(1u,2u),tiny),rational(1u,2u),
            plus(rational(1u,2u),tiny),plus(rational(1u,2u),eps),
            minus(rational(5u,7u),eps),rational(5u,7u),plus(rational(5u,7u),eps),
            rational(3u,4u),rational(1u),rational(2u),rational(3u)};
        for(unsigned p:{53u,7u,128u,1024u,256u,53u})
            for(auto const& v:d89)both(p,v);
        if(count!=1340u)throw std::runtime_error("D89 count "+std::to_string(count));
        out+="# Exact dyadic centered slope correction: numerator divisible by five.\n";
        both(128u,minus(rational(1u,2u),power_ratio(5u,128u)));
        if(count!=1348u)throw std::runtime_error("exact correction count "+std::to_string(count));
    }
    static Rational rational_shift(Rational const& a,std::uint64_t factor) {
        return reduce({a.n.mul_small(factor),a.d});
    }
};
}
int main(int argc,char**argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        Generator generator;
        generator.generate();
        independent_oracle::output(dir,"floatmp_atan.tsv",generator.out,write);
        return 0;
    }catch(std::exception const& e) {
        std::cerr<<"FloatMP atan generator: "<<e.what()<<'\n';
        return 1;
    }
}
