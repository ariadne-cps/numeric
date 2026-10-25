// SPDX-License-Identifier: MIT
// Independently certified FloatAP logarithm fixture, standard C++20 integers only.
#include "oracle_exact_integer.hpp"
#include <set>
#include <tuple>

namespace {
using independent_oracle::Nat;
using independent_oracle::row;
constexpr std::size_t work_bits=1800u;
std::size_t bits(Nat const& n) {
    return n.zero()?0u:(n.w.size()-1u)*64u+64u-
        static_cast<std::size_t>(std::countl_zero(n.w.back()));
}
Nat multiply(Nat const& a,Nat const& b) {
    if(a.zero()||b.zero())return {};
    const auto unpack=[](Nat const& n){
        std::vector<std::uint32_t> d;
        d.reserve(2u*n.w.size());
        for(auto word:n.w){d.push_back(static_cast<std::uint32_t>(word));
            d.push_back(static_cast<std::uint32_t>(word>>32u));}
        while(!d.empty()&&d.back()==0u)d.pop_back();
        return d;
    };
    const auto x=unpack(a),y=unpack(b);
    std::vector<std::uint32_t> p(x.size()+y.size(),0u);
    for(std::size_t i=0;i<x.size();++i){
        std::uint64_t carry=0;
        for(std::size_t j=0;j<y.size();++j){
            const std::uint64_t v=std::uint64_t{x[i]}*y[j]+p[i+j]+carry;
            p[i+j]=static_cast<std::uint32_t>(v);
            carry=v>>32u;
        }
        std::size_t k=i+y.size();
        while(carry){
            if(k==p.size())p.push_back(0u);
            const auto v=std::uint64_t{p[k]}+carry;
            p[k]=static_cast<std::uint32_t>(v);
            carry=v>>32u;
            ++k;
        }
    }
    std::vector<std::uint64_t> limbs((p.size()+1u)/2u,0u);
    for(std::size_t i=0;i<p.size();++i)
        limbs[i/2u]|=std::uint64_t{p[i]}<<((i%2u)*32u);
    return Nat(std::move(limbs));
}
std::pair<Nat,Nat> divrem(Nat numerator,Nat const& denominator){
    if(denominator.zero())throw std::invalid_argument("zero denominator");
    Nat quotient;
    while(numerator.compare(denominator)>=0){
        std::size_t shift=bits(numerator)-bits(denominator);
        Nat subtrahend=denominator.shl(shift);
        if(subtrahend.compare(numerator)>0)subtrahend=denominator.shl(--shift);
        numerator=Nat::sub(numerator,subtrahend);
        quotient=Nat::add(quotient,Nat::power2(shift));
    }
    return {std::move(quotient),std::move(numerator)};
}
Nat ceildiv(Nat const& n,Nat const& d){
    auto [q,r]=divrem(n,d);
    return r.zero()?q:Nat::add(q,Nat(1u));
}
struct Interval {Nat lo,hi;};
Nat const& B(){static const Nat b=Nat::power2(work_bits);return b;}
Interval fixed_rational(Nat const& n,Nat const& d){
    if(d.zero())throw std::invalid_argument("zero divisor");
    auto [q,r]=divrem(n.shl(work_bits),d);
    return {q,r.zero()?q:Nat::add(q,Nat(1u))};
}
Interval add(Interval const& a,Interval const& b){
    return {Nat::add(a.lo,b.lo),Nat::add(a.hi,b.hi)};
}
Interval mul(Interval const& a,Interval const& b){
    return {multiply(a.lo,b.lo).shr(work_bits),
        ceildiv(multiply(a.hi,b.hi),B())};
}
Interval divby(Interval const& a,std::uint64_t denominator){
    const Nat d(denominator);
    return {divrem(a.lo,d).first,ceildiv(a.hi,d)};
}
Interval log_abs(Nat const& n,Nat const& d){
    const int order=n.compare(d);
    if(order==0)return {{},{}};
    const Nat diff=order>0?Nat::sub(n,d):Nat::sub(d,n);
    const Nat sum=Nat::add(n,d);
    const Interval y=fixed_rational(diff,sum),y2=mul(y,y);
    Interval term=y,total{{},{}};
    for(unsigned k=0;k<=1100u;++k){
        total=add(total,divby(term,2u*k+1u));
        term=mul(term,y2);
    }
    // Omitted positive sum bounded by 2*y^2205 in atanh scale;
    // after outer factor 2, the final error is <= 4*y^2205.
    // term here encloses y^2203, one more multiplication encloses y^2205.
    const Interval tail=mul(term,y2);
    return {total.lo.shl(1u),Nat::add(total.hi.shl(1u),tail.hi.mul_small(4u))};
}
struct Fraction{Nat n,d;bool negative=false;};
Fraction reduce_dyadic(Nat magnitude,std::int64_t e,bool negative=false){
    if(magnitude.zero())return {Nat(),Nat(1u),false};
    if(e<0){
        const auto cut=std::min(static_cast<std::size_t>(-e),magnitude.trailing_zero_bits());
        magnitude=magnitude.shr(cut);
        e+=static_cast<std::int64_t>(cut);
    }
    if(e>=0)magnitude=magnitude.shl(static_cast<std::size_t>(e));
    return {std::move(magnitude),e<0?Nat::power2(static_cast<std::size_t>(-e)):Nat(1u),negative};
}
std::string number(Fraction const& f){return (f.negative?"-":"")+f.n.dec();}
Nat quantize(Nat q,Nat const& r,Nat const& d,char mode,bool negative){
    if(r.zero())return q;
    bool increment=false;
    if(mode=='N'){
        const int c=r.shl(1u).compare(d);
        increment=c>0||(c==0&&!q.zero()&&(q.w.front()&1u)!=0u);
    }else if(mode=='D')increment=negative;
    else if(mode=='U')increment=!negative;
    else if(mode!='Z')throw std::invalid_argument("rounding mode");
    return increment?Nat::add(q,Nat(1u)):q;
}
std::int64_t exponent(Nat const& n,Nat const& d){
    std::int64_t e=static_cast<std::int64_t>(bits(n))-static_cast<std::int64_t>(bits(d));
    if(e>=0){if(n.compare(d.shl(static_cast<std::size_t>(e)))<0)--e;}
    else if(n.shl(static_cast<std::size_t>(-e)).compare(d)<0)--e;
    return e;
}
Fraction rounded_rational(Nat const& n,Nat const& d,unsigned precision,char mode,bool negative=false){
    if(n.zero())return {Nat(),Nat(1u),false};
    const std::int64_t grid=exponent(n,d)-static_cast<std::int64_t>(precision)+1;
    const Nat numerator=grid<0?n.shl(static_cast<std::size_t>(-grid)):n;
    const Nat denominator=grid>=0?d.shl(static_cast<std::size_t>(grid)):d;
    auto [q,r]=divrem(numerator,denominator);
    return reduce_dyadic(quantize(std::move(q),r,denominator,mode,negative),grid,negative);
}
Fraction rounded_fixed(Nat const& fixed,unsigned precision,char mode,bool negative){
    if(fixed.zero())return {Nat(),Nat(1u),false};
    const auto grid=static_cast<std::int64_t>(bits(fixed))-
        static_cast<std::int64_t>(work_bits)-static_cast<std::int64_t>(precision);
    const auto shift=static_cast<std::int64_t>(bits(fixed))-
        static_cast<std::int64_t>(precision);
    if(shift<=0)return reduce_dyadic(fixed.shl(static_cast<std::size_t>(-shift)),grid,negative);
    const auto k=static_cast<std::size_t>(shift);
    Nat q=fixed.shr(k);
    const Nat r=Nat::sub(fixed,q.shl(k));
    return reduce_dyadic(quantize(std::move(q),r,Nat::power2(k),mode,negative),grid,negative);
}
Fraction certify(Interval const& magnitude,unsigned precision,char mode,bool negative){
    // For a negative result, the lower numeric bound is -magnitude.hi.
    const auto answer=rounded_fixed(magnitude.lo,precision,mode,negative);
    const auto other=rounded_fixed(magnitude.hi,precision,mode,negative);
    if(number(answer)!=number(other)||answer.d.compare(other.d)!=0)
        throw std::runtime_error("logarithm rounding not certified at precision "+
                                 std::to_string(precision)+" mode "+mode);
    return answer;
}
Fraction half_offset(unsigned power,std::int64_t offset){
    // 1/2 + offset*2^-power for the original centered logarithm cases.
    const Nat d=Nat::power2(power);
    Nat n=Nat::power2(power-1u);
    if(offset<0)n=Nat::sub(n,Nat(static_cast<std::uint64_t>(-offset)));
    else n=Nat::add(n,Nat(static_cast<std::uint64_t>(offset)));
    return {std::move(n),std::move(d),false};
}
std::string generate(){
    std::string out=
        "# Certified correctly rounded FloatAP logarithm values, D84.\n"
        "# precision_bits mode input_num input_den expected_num expected_den\n"
        "# Source: centered_log existing 15 Rational inputs, 6 precision visits, 4 modes (360 checks).\n"
        "# Each input value is independently rounded nearest-even to the original precision.\n"
        "# For x>0 use y=abs(x-1)/(x+1)<=3/5; abs(log(x))=2*atanh(y).\n"
        "# Set B=2^1800. All integer fixed-point products and quotients are rounded outward.\n"
        "# Sum y^(2k+1)/(2k+1) for k=0..1100. Upper remainder <= 4*y^(2205).\n"
        "# The sign of log(x) is recovered by reversing bounds when x<1.\n"
        "# Accept ONLY when lower and upper endpoints yield the SAME signed binary rounding.\n"
        "# Modes: N nearest-even, D downward, U upward, Z toward zero.\n"
        "# No GMP, independent exact arithmetic, FloatAP or floating-point transcendental oracle used in generation.\n"
        "# 284 unique data rows for 360 original checks.\n";
    std::set<std::tuple<unsigned,char,std::string,std::string>> keys;
    std::size_t count=0;
    const auto append=[&](unsigned p,Fraction const& arg){
        const std::string n=number(arg),d=arg.d.dec();
        const bool negative=arg.n.compare(arg.d)<0;
        const auto bound=(arg.n.compare(Nat(8u))==0 && arg.d.compare(Nat(1u))==0)||
            (arg.n.compare(Nat(1u))==0 && arg.d.compare(Nat(8u))==0)
            ? [&](){const auto v=log_abs(Nat(2u),Nat(1u));return Interval{v.lo.mul_small(3u),v.hi.mul_small(3u)};}()
            : log_abs(arg.n,arg.d);
        for(char mode:{'N','D','U','Z'}){
            if(!keys.emplace(p,mode,n,d).second)continue;
            Fraction ans;
            try{ans=certify(bound,p,mode,negative);}
            catch(std::exception const& e){throw std::runtime_error(std::string(e.what())+" input "+n+"/"+d);}
            out+=row({std::to_string(p),std::string(1u,mode),n,d,
                      number(ans),ans.d.dec()});
            ++count;
        }
    };
    const auto make=[](std::uint64_t n,std::uint64_t d){return Fraction{Nat(n),Nat(d),false};};
    for(unsigned p:{53u,7u,128u,256u,1024u,53u}){
        const unsigned cutoff=(p+33u)/2u;
        const std::array<Fraction,15u> inputs={
            make(1u,2u),half_offset(p,-1),half_offset(p,1),
            half_offset(p,-4),half_offset(p,4),
            half_offset(cutoff,-1),half_offset(cutoff,1),
            half_offset(cutoff,-2),half_offset(cutoff,2),
            make(1u,4u),make(3u,4u),make(1u,3u),make(2u,3u),
            make(1u,1u),make(2u,1u)
        };
        for(auto const& input:inputs){
            const auto rounded=rounded_rational(input.n,input.d,p,'N');
            append(p,rounded);
        }
    }
    // D96 original (not pre-quantized!) input keys, and their original order.
    out+="# D96: log range reduction, cached log(2), boundary and cancellation probes.\n";
    const Nat unit=Nat::power2(40u);
    const auto epsilon=[&](std::uint64_t twice,std::int64_t adjust){
        Nat n=Nat(twice).shl(39u);
        if(adjust<0)n=Nat::sub(n,Nat(static_cast<std::uint64_t>(-adjust)));
        else n=Nat::add(n,Nat(static_cast<std::uint64_t>(adjust)));
        return Fraction{std::move(n),unit,false};
    };
    const std::array<Fraction,13u> d96={
        make(1,2),make(3,4),make(3,2),make(2,1),make(1,8),make(8,1),
        epsilon(2,-1),epsilon(2,1),
        Fraction{Nat::sub(Nat(3u).shl(38u),Nat(1u)),unit,false},
        Fraction{Nat::add(Nat(3u).shl(38u),Nat(1u)),unit,false},
        epsilon(3,-1),epsilon(3,1),epsilon(1,1)
    };
    for(unsigned p:{53u,128u,256u,1024u})
        for(const auto& source:d96){
            // Canonicalize the exact dyadic source (not via FloatAP input rounding).
            Nat n=source.n,d=source.d;
            const auto cancel=std::min(n.trailing_zero_bits(),d.trailing_zero_bits());
            n=n.shr(cancel);d=d.shr(cancel);
            append(p,{std::move(n),std::move(d),false});
        }
    if(count!=444u)throw std::runtime_error("unexpected log oracle row count: "+std::to_string(count));
    return out;
}
}
int main(int argc,char**argv){
    try{
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"floatap_log.tsv",generate(),write);
        return 0;
    }catch(std::exception const& e){
        std::cerr<<"FloatAP log generator: "<<e.what()<<'\n';
        return 1;
    }
}
