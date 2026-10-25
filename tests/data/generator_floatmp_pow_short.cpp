// SPDX-License-Identifier: MIT
// Independent exact rational-to-dyadic integer power oracle, standard C++20.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;
struct Dyadic {
    Nat mantissa;
    std::int64_t exponent=0;
    bool negative=false;
};
std::size_t bit_length(Nat const& n) {
    return n.zero()?0u:(n.w.size()-1u)*64u+
        64u-static_cast<std::size_t>(std::countl_zero(n.w.back()));
}
std::pair<Nat,Nat> divrem(Nat remainder,Nat const& denominator) {
    if(denominator.zero())throw std::invalid_argument("zero divisor");
    Nat quotient;
    while(remainder.compare(denominator)>=0) {
        std::size_t shift=bit_length(remainder)-bit_length(denominator);
        Nat aligned=denominator.shl(shift);
        if(aligned.compare(remainder)>0)aligned=denominator.shl(--shift);
        remainder=Nat::sub(remainder,aligned);
        quotient=Nat::add(quotient,Nat::power2(shift));
    }
    return {std::move(quotient),std::move(remainder)};
}
std::int64_t binary_exponent(Nat const& n,Nat const& d) {
    std::int64_t e=static_cast<std::int64_t>(bit_length(n))-
        static_cast<std::int64_t>(bit_length(d));
    if(e>=0) {
        if(n.compare(d.shl(static_cast<std::size_t>(e)))<0)--e;
    } else if(n.shl(static_cast<std::size_t>(-e)).compare(d)<0)--e;
    return e;
}
Nat rounded(Nat q,Nat const& r,Nat const& d,char mode,bool negative) {
    if(r.zero())return q;
    bool increment=false;
    if(mode=='N') {
        const int cmp=r.shl(1u).compare(d);
        increment=cmp>0||(cmp==0&&!q.zero()&&(q.w.front()&1u)!=0u);
    } else if(mode=='D')increment=negative;
    else if(mode=='U')increment=!negative;
    else if(mode!='Z')throw std::invalid_argument("rounding mode");
    return increment?Nat::add(q,Nat(1u)):q;
}
Dyadic round_source(Nat const& n,Nat const& d,bool negative,unsigned precision) {
    if(n.zero())return {};
    const std::int64_t grid=binary_exponent(n,d)-static_cast<std::int64_t>(precision)+1;
    const Nat numerator=grid<0?n.shl(static_cast<std::size_t>(-grid)):n;
    const Nat denominator=grid>=0?d.shl(static_cast<std::size_t>(grid)):d;
    auto [q,r]=divrem(numerator,denominator);
    return {rounded(std::move(q),r,denominator,'N',negative),grid,negative};
}
Nat multiply(Nat lhs,Nat rhs) {
    Nat result;
    while(!rhs.zero()) {
        if((rhs.w.front()&1u)!=0u)result=Nat::add(result,lhs);
        rhs=rhs.shr(1u);
        if(!rhs.zero())lhs=lhs.shl(1u);
    }
    return result;
}
Nat power(Nat base,unsigned exponent) {
    Nat result(1u);
    while(exponent!=0u) {
        if((exponent&1u)!=0u)result=multiply(result,base);
        exponent>>=1u;
        if(exponent!=0u)base=multiply(base,base);
    }
    return result;
}
Dyadic exact_power(Dyadic const& value,unsigned degree) {
    if(degree==0u)return {Nat(1u),0,false};
    if(value.mantissa.zero())return {};
    return {power(value.mantissa,degree),value.exponent*static_cast<std::int64_t>(degree),
        value.negative&&(degree&1u)!=0u};
}
Dyadic round_dyadic(Dyadic const& value,unsigned precision,char mode) {
    if(value.mantissa.zero())return {};
    const std::int64_t e=static_cast<std::int64_t>(bit_length(value.mantissa))-1+value.exponent;
    const std::int64_t grid=e-static_cast<std::int64_t>(precision)+1;
    const std::int64_t difference=value.exponent-grid;
    if(difference>=0)return {value.mantissa.shl(static_cast<std::size_t>(difference)),grid,value.negative};
    const std::size_t shift=static_cast<std::size_t>(-difference);
    Nat q=value.mantissa.shr(shift);
    const Nat r=Nat::sub(value.mantissa,q.shl(shift));
    const Nat d=Nat::power2(shift);
    return {rounded(std::move(q),r,d,mode,value.negative),grid,value.negative};
}
std::pair<std::string,std::string> fraction(Dyadic const& value) {
    if(value.mantissa.zero())return {"0","1"};
    Nat magnitude=value.mantissa;
    std::int64_t e=value.exponent;
    if(e<0) {
        const std::size_t cancel=std::min(static_cast<std::size_t>(-e),magnitude.trailing_zero_bits());
        magnitude=magnitude.shr(cancel);
        e+=static_cast<std::int64_t>(cancel);
    }
    if(e>=0)magnitude=magnitude.shl(static_cast<std::size_t>(e));
    const Nat denominator=e<0?Nat::power2(static_cast<std::size_t>(-e)):Nat(1u);
    return {(value.negative?"-":"")+magnitude.dec(),denominator.dec()};
}
struct Case {
    Nat scaled;
    Nat numerator;
    Nat denominator;
    bool negative;
    unsigned degree;
};
std::vector<Case> source_cases(unsigned precision) {
    std::vector<Case> cases;
    Nat radix=Nat::power2(precision);
    const Nat half=Nat::power2(precision-1u);
    const Nat low=Nat::power2(precision+3u);
    const Nat high=low.shl(1u);
    const auto add=[&](Nat n,Nat d,bool negative,unsigned degree) {
        Nat scaled=d.compare(Nat(1u))==0?n.shl(precision):n;
        for(const auto& x:cases)if(x.negative==negative&&x.degree==degree&&x.scaled.compare(scaled)==0)return;
        cases.push_back({std::move(scaled),std::move(n),std::move(d),negative,degree});
    };
    for(unsigned degree:{2u,3u,7u,8u,13u,16u,17u}) {
        Nat threshold=Nat::power2(precision>2u?(precision-3u)/2u:0u);
        const std::uint64_t t=divrem(threshold,Nat(degree)).first.zero()?0u:
            divrem(threshold,Nat(degree)).first.w.front();
        std::vector<std::uint64_t> around{0u,1u};
        if(t>0u) around.push_back(t-1u);
        around.push_back(t); around.push_back(t+1u);
        for(bool negative:{true,false}) {
            for(auto o:around) {
                Nat increment(o);
                Nat new_low=Nat::add(half,increment);
                if(new_low.compare(radix)<0) add(new_low,radix,negative,degree);
                Nat new_large=Nat::add(low,increment.shl(4u));
                if(new_large.compare(high)<0)add(new_large,Nat(1u),negative,degree);
            }
            const Nat upper_threshold=divrem(Nat::power2((precision-2u)/2u),Nat(degree)).first;
            const std::uint64_t upper_t=upper_threshold.zero()?0u:upper_threshold.w.front();
            std::vector<std::uint64_t> upper{1u};
            if(upper_t>0u)upper.push_back(upper_t-1u);
            upper.push_back(upper_t);upper.push_back(upper_t+1u);
            for(auto o:upper) {
                if(o==0u)continue;
                Nat delta(o);
                if(delta.compare(radix)<0)add(Nat::sub(radix,delta),radix,negative,degree);
                Nat delta_large=delta.shl(4u);
                if(delta_large.compare(low)<0)add(Nat::sub(high,delta_large),Nat(1u),negative,degree);
            }
        }
    }
    std::sort(cases.begin(),cases.end(),[](Case const& a,Case const& b) {
        if(a.negative!=b.negative)return a.negative;
        const int cmp=a.scaled.compare(b.scaled);
        if(cmp!=0)return a.negative?cmp>0:cmp<0;
        return a.degree<b.degree;
    });
    return cases;
}
std::string generate() {
    std::string out=
        "# Independent exact short-power boundary oracle.\n"
        "# bits mode input_num input_den degree expected_num expected_den\n"
        "# Generated offline with Python signed integers and exact Fraction arithmetic.\n"
        "# Covers both binade edges, both signs, exponents, p=2/7/53/64/65, degrees\n"
        "# 2/3/7/8/13/16/17, and the last accepted/first rejected shortcut distances.\n"
        "# Input is first rounded nearest-even, exact power is then rounded ONCE.\n";
    for(unsigned precision:{2u,7u,53u,64u,65u}) {
        const auto cases=source_cases(precision);
        for(char mode:{'D','N','U','Z'}) {
            for(const auto& input:cases) {
                Nat n(input.numerator.w),d(input.denominator.w);
                const Dyadic initial=round_source(n,d,input.negative,precision);
                const Dyadic exact=exact_power(initial,input.degree);
                const auto [num,den]=fraction(round_dyadic(exact,precision,mode));
                const auto [src_num,src_den]=fraction(Dyadic{input.numerator,
                    input.denominator.compare(Nat(1u))==0?0:-static_cast<std::int64_t>(precision),input.negative});
                out+=row({std::to_string(precision),std::string(1u,mode),src_num,src_den,
                    std::to_string(input.degree),num,den});
            }
        }
    }
    return out;
}
}
int main(int argc,char**argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"floatmp_pow_short.tsv",generate(),write);
        return 0;
    } catch(std::exception const& ex) {
        std::cerr<<"FloatMP pow generator: "<<ex.what()<<'\n';
        return 1;
    }
}
