// SPDX-License-Identifier: MIT
// Independent correctly rounded FloatAP asin/acos references, C++20 only.
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
    if(first!=second)throw std::runtime_error("inverse-trigonometric interval not certified at bits="+
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

// Digit-pair restoring square root of a nonnegative arbitrary-width integer.
// Invariant: input = root^2 + remainder when all input bit pairs consumed.
std::pair<Nat,Nat> isqrt_remainder(Nat input) {
    if(input.zero())return {Nat(),Nat()};
    const std::size_t top=2u*((bit_length(input)-1u)/2u);
    Nat bit=Nat::power2(top),root;
    while(!bit.zero()) {
        const Nat trial=Nat::add(root,bit);
        if(input.compare(trial)>=0) {
            input=Nat::sub(input,trial);
            root=Nat::add(root.shr(1u),bit);
        } else root=root.shr(1u);
        bit=bit.shr(2u);
    }
    return {std::move(root),std::move(input)};
}
struct SqrtInterval { Nat lo,hi; };
SqrtInterval sqrt_one_minus_square(Rational const& input) {
    const Nat d2=multiply(input.d,input.d);
    const Nat n2=multiply(input.n,input.n);
    if(n2.compare(d2)>0)throw std::invalid_argument("asin/acos input outside [-1,1]");
    const Nat exact_num=Nat::sub(d2,n2).shl(2u*scale_bits);
    const auto [q,quotient_remainder]=divide(exact_num,d2);
    const auto [root,sqrt_remainder]=isqrt_remainder(q);
    const bool exact=quotient_remainder.zero()&&sqrt_remainder.zero();
    return {root,exact?root:Nat::add(root,Nat(1u))};
}
Interval arcsine_magnitude(Rational const& input) {
    if(input.n.zero())return number(0u);
    if(input.n.compare(input.d)==0)return scale(arctangent_positive(Nat(1u),Nat(1u)),2u);
    const auto s=sqrt_one_minus_square(input);
    // atan(x/(1+sqrt(1-x*x)))*2, using explicit rational bounds.
    // At fixed point B, the two rational bounds for the atan argument are:
    //     n*B / (d*(B+upper_sqrt))
    //     n*B / (d*(B+lower_sqrt))
    const Nat numerator=input.n.shl(scale_bits);
    const Nat lower_d=multiply(input.d,Nat::add(scale(),s.hi));
    const Nat upper_d=multiply(input.d,Nat::add(scale(),s.lo));
    const Interval lower=arctangent_positive(numerator,lower_d);
    const Interval upper=arctangent_positive(numerator,upper_d);
    return scale({lower.lo,upper.hi},2u);
}
struct Generator {
    std::string asin_out=
        "# Correctly rounded asin certified golden references (D86).\n"
        "# bits mode source_numerator source_denominator result_numerator result_denominator\n"
        "# Original exact test-source input is the key; input conversion is independently\n"
        "# nearest/even rounded before exact integer-rational interval certification.\n"
        "# B=2^1800; outward sqrt (integer isqrt), atan alternating Taylor k=0..1149;\n"
        "# pi/2 = 8 atan(1/5)-2 atan(1/239); asin=2 atan(x/(1+sqrt(1-x^2))).\n"
        "# Accept a value only when BOTH enclosure endpoints round to the same\n"
        "# exact binary Rational at the requested bits and direction.\n"
        "# No GMP, independent exact arithmetic, independent exact arithmetic, FloatAP, or external transcendental oracle.\n";
    std::string acos_out=asin_out;
    std::set<std::tuple<unsigned,char,std::string,std::string>> seen;
    std::map<std::pair<std::string,std::string>,Interval> cache;
    std::size_t count=0;
    Generator() {
        constexpr const char* asin_title="Correctly rounded asin certified golden references";
        constexpr const char* acos_title="Correctly rounded acos certified golden references";
        const auto index=acos_out.find(asin_title);
        if(index==std::string::npos)throw std::runtime_error("header string error");
        acos_out.replace(index,std::char_traits<char>::length(asin_title),acos_title);
    }
    void emit(unsigned p,Rational const& val,bool negative) {
        // The golden key is the original exact source, not the rounded input.
        const auto num=(negative&&!val.n.zero()?"-":"")+val.n.dec();
        const auto den=val.d.dec();
        bool needed=false;
        for(char m:{'N','D','U','Z'})
            if(!seen.contains({p,m,num,den}))needed=true;
        if(!needed)return;
        const auto rounded=nearest(val,p);
        const auto key=std::make_pair(rounded.n.dec(),rounded.d.dec());
        auto it=cache.find(key);
        if(it==cache.end())it=cache.emplace(key,arcsine_magnitude(rounded)).first;
        const Interval signed_asin=negative?neg(it->second):it->second;
        const Interval pi2=scale(arctangent_positive(Nat(1u),Nat(1u)),2u);
        const Interval signed_acos=rounded.n.compare(rounded.d)==0&&!negative?
            number(0u):add(pi2,neg(signed_asin));
        for(char mode:{'N','D','U','Z'}) {
            if(!seen.emplace(p,mode,num,den).second)continue;
            const auto a=certify(signed_asin,p,mode);
            const auto c=certify(signed_acos,p,mode);
            const std::vector<std::string> keyfields={std::to_string(p),std::string(1u,mode),num,den};
            auto result=keyfields;
            result.push_back(a.first);result.push_back(a.second);
            asin_out+=row(result);
            result=keyfields;
            result.push_back(c.first);result.push_back(c.second);
            acos_out+=row(result);
            ++count;
        }
    }
    void emit_both(unsigned p,Rational const& val) {
        emit(p,val,true);
        emit(p,val,false);
    }
    void generate() {
        const auto half=rational(1u,2u);
        for(unsigned p:{53u,7u,128u,256u,1024u,53u}) {
            const auto ulp=power_ratio(1u,p);
            const auto threshold=power_ratio(1u,(p+33u)/2u);
            const std::array<Rational,13u> sources={
                half,minus(half,ulp),plus(half,ulp),
                minus(half,scale_ratio(ulp,4u)),plus(half,scale_ratio(ulp,4u)),
                minus(half,threshold),plus(half,threshold),
                minus(half,scale_ratio(threshold,2u)),plus(half,scale_ratio(threshold,2u)),
                rational(1u,4u),rational(3u,4u),rational(1u,3u),rational(2u,3u)
            };
            for(auto const& source:sources)emit_both(p,source);
        }
        if(count!=520u)throw std::runtime_error("D86 row count "+std::to_string(count));
        asin_out+="# D90: certified one-evaluation asin enclosure and endpoint probes.\n";
        acos_out+="# D90: certified shared acos enclosure and endpoint probes.\n";
        const auto eps=power_ratio(1u,40u);
        for(unsigned p:{53u,128u,256u,1024u}) {
            const auto signed_source=[&](Rational const& v,bool negative){emit(p,v,negative);};
            signed_source(rational(0u),false);
            signed_source(rational(1u,64u),false);
            signed_source(rational(1u,64u),true);
            signed_source(half,false);
            signed_source(half,true);
            signed_source(minus(half,eps),false);
            signed_source(plus(half,eps),false);
            signed_source(minus(half,eps),true);
            signed_source(plus(half,eps),true);
            signed_source(rational(3u,4u),false);
            signed_source(rational(3u,4u),true);
            signed_source(rational(7u,8u),false);
            signed_source(rational(7u,8u),true);
            signed_source(minus(rational(1u),eps),false);
            signed_source(minus(rational(1u),eps),true);
            signed_source(rational(1u),false);
            signed_source(rational(1u),true);
        }
        if(count!=728u)throw std::runtime_error("D90 row count "+std::to_string(count));
    }
    static Rational scale_ratio(Rational const& x,std::uint64_t k) {
        return reduce({x.n.mul_small(k),x.d});
    }
};
}
int main(int argc,char**argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        Generator generator;
        generator.generate();
        independent_oracle::output(dir,"floatap_asin.tsv",generator.asin_out,write);
        independent_oracle::output(dir,"floatap_acos.tsv",generator.acos_out,write);
        return 0;
    }catch(std::exception const& e) {
        std::cerr<<"FloatAP asin/acos generator: "<<e.what()<<'\n';
        return 1;
    }
}
