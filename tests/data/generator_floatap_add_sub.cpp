// SPDX-License-Identifier: MIT
// Reproducible exact FloatAP add/sub fixture, using only independent integer limbs.
// Standalone standard C++20: no Numeric, GMP, MPFR, FloatMP, Boost, or host floats.
#include "oracle_exact_integer.hpp"
#include <set>
#include <tuple>

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
Dyadic nearest_ratio(Nat const& n,Nat const& d,unsigned precision) {
    if(n.zero())return {};
    const std::int64_t grid=binary_exponent(n,d)-static_cast<std::int64_t>(precision)+1;
    const Nat numerator=grid<0?n.shl(static_cast<std::size_t>(-grid)):n;
    const Nat denominator=grid>=0?d.shl(static_cast<std::size_t>(grid)):d;
    auto [q,r]=divrem(numerator,denominator);
    return {rounded(std::move(q),r,denominator,'N',false),grid,false};
}
Dyadic exact_sum(Dyadic const& a,Dyadic const& b,bool subtract) {
    if(a.mantissa.zero()){
        Dyadic answer=b;
        if(subtract && !answer.mantissa.zero())answer.negative=!answer.negative;
        return answer;
    }
    if(b.mantissa.zero())return a;
    const auto exponent=std::min(a.exponent,b.exponent);
    const Nat left=a.mantissa.shl(static_cast<std::size_t>(a.exponent-exponent));
    const Nat right=b.mantissa.shl(static_cast<std::size_t>(b.exponent-exponent));
    const bool right_negative=b.negative!=subtract;
    if(a.negative==right_negative)
        return {Nat::add(left,right),exponent,a.negative};
    const int order=left.compare(right);
    if(order==0)return {};
    return order>0?Dyadic{Nat::sub(left,right),exponent,a.negative}
                  :Dyadic{Nat::sub(right,left),exponent,right_negative};
}
Dyadic round_dyadic(Dyadic const& source,unsigned precision,char mode) {
    if(source.mantissa.zero())return {};
    const std::int64_t e=static_cast<std::int64_t>(bit_length(source.mantissa))-1+source.exponent;
    const std::int64_t grid=e-static_cast<std::int64_t>(precision)+1;
    const std::int64_t difference=source.exponent-grid;
    if(difference>=0)return {source.mantissa.shl(static_cast<std::size_t>(difference)),grid,source.negative};
    const auto shift=static_cast<std::size_t>(-difference);
    Nat q=source.mantissa.shr(shift);
    const Nat r=Nat::sub(source.mantissa,q.shl(shift));
    const Nat d=Nat::power2(shift);
    return {rounded(std::move(q),r,d,mode,source.negative),grid,source.negative};
}
std::pair<std::string,std::string> fraction(Dyadic const& value) {
    if(value.mantissa.zero())return {"0","1"};
    Nat magnitude=value.mantissa;
    std::int64_t e=value.exponent;
    if(e<0) {
        const auto cancel=std::min(static_cast<std::size_t>(-e),magnitude.trailing_zero_bits());
        magnitude=magnitude.shr(cancel);
        e+=static_cast<std::int64_t>(cancel);
    }
    if(e>=0)magnitude=magnitude.shl(static_cast<std::size_t>(e));
    const Nat d=e<0?Nat::power2(static_cast<std::size_t>(-e)):Nat(1u);
    return {(value.negative?"-":"")+magnitude.dec(),d.dec()};
}
std::string generate() {
    std::string out=
        "# FloatAP exact add/sub oracle v1: 10x10 ordered pairs x 2 operations x 4 rounding modes x 5 precisions.\n"
        "# precision_bits mode input_num input_den expected_num expected_den\n"
        "# Modes N nearest/ties-even, D downward, U upward, Z toward zero.\n"
        "# Operands are exact FloatAP-representable dyadics: 1, -1, rounded(7/11), its negation,\n"
        "# +/-2^(-bits), +/- powers of two with exponent +/-200 and +/-1400.\n"
        "# Input is the exact sum or difference BEFORE final rounding.\n"
        "# Values computed using independent arbitrary-precision integer/rational algorithms.\n"
        "# D80; no FloatAP/FloatMP/GMP/MPFR involved in generation.\n";
    std::set<std::tuple<unsigned,char,std::string,std::string>> seen;
    for(unsigned precision:{7u,53u,128u,256u,1024u}) {
        const Dyadic one{Nat(1u),0,false};
        const Dyadic minus_one{Nat(1u),0,true};
        const Dyadic ratio=nearest_ratio(Nat(7u),Nat(11u),precision);
        const Dyadic minus_ratio{ratio.mantissa,ratio.exponent,true};
        const Dyadic tiny{Nat(1u),-static_cast<std::int64_t>(precision),false};
        const Dyadic minus_tiny{Nat(1u),-static_cast<std::int64_t>(precision),true};
        const std::array<Dyadic,10u> values={one,minus_one,ratio,minus_ratio,
            tiny,minus_tiny,{Nat(1u),-200,false},{Nat(1u),200,false},
            {Nat(1u),-1400,false},{Nat(1u),1400,false}};
        for(auto const& a:values)for(auto const& b:values) {
            const Dyadic sum=exact_sum(a,b,false),difference=exact_sum(a,b,true);
            for(char mode:{'N','D','U','Z'}) {
                for(auto const& exact:{sum,difference}) {
                    const auto source=fraction(exact);
                    if(!seen.emplace(precision,mode,source.first,source.second).second)continue;
                    const auto output=fraction(round_dyadic(exact,precision,mode));
                    out+=row({std::to_string(precision),std::string(1u,mode),source.first,
                              source.second,output.first,output.second});
                }
            }
        }
    }
    return out;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"floatap_add_sub.tsv",generate(),write);
        return 0;
    } catch(std::exception const& e) {
        std::cerr<<"FloatAP add/sub generator: "<<e.what()<<'\n';
        return 1;
    }
}
