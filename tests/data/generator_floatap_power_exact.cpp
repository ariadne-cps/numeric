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
std::string generate() {
    std::string out=
        "# FloatAP integer pow golden reference v1\n"
        "# precision_bits mode input_num input_den exponent expected_num expected_den\n"
        "# Source is the exact signed Rational supplied to the existing FloatAP constructor.\n"
        "# Offline oracle FIRST rounds signed Rational input to nearest-even at precision_bits.\n"
        "# It then computes the exact integer power of the rounded dyadic,\n"
        "# and finally rounds ONCE at the same precision in the requested mode.\n"
        "# This preserves the original input list and its 0^0 = +1 behavior.\n"
        "# Signed-zero behavior is checked separately using dedicated assertions.\n"
        "# Derived solely by signed arbitrary-precision integer division and exponentiation;\n"
        "# no FloatAP, FloatMP, GMP, MPFR, or floating-point rounding invoked.\n"
        "# Four modes: N nearest even, D downward, U upward, Z toward zero.\n"
        "# 5 bit widths x 8 magnitudes x 2 signs x 8 powers x 4 modes = 2560\n"
        "# exercised cases, 2400 distinct reference keys (negative and positive zero merge).\n";
    const std::array<std::pair<std::uint64_t,std::uint64_t>,8u> inputs={{{0u,1u},{1u,1u},{1u,2u},{3u,4u},
        {17u,16u},{5u,3u},{3u,1u},{131u,128u}}};
    for(unsigned precision:{7u,53u,128u,256u,1024u}) {
        for(const auto& [n,d]:inputs) {
            for(bool negative:{true,false}) {
                if(n==0u&&negative)continue;
                const Dyadic base=round_source(Nat(n),Nat(d),negative,precision);
                for(unsigned degree:{0u,1u,2u,3u,7u,8u,13u,31u}) {
                    const Dyadic exact=exact_power(base,degree);
                    for(char mode:{'N','D','U','Z'}) {
                        const Dyadic final=round_dyadic(exact,precision,mode);
                        const auto [num,den]=fraction(final);
                        out+=row({std::to_string(precision),std::string(1u,mode),
                            (negative?"-":"")+std::to_string(n),std::to_string(d),
                            std::to_string(degree),num,den});
                    }
                }
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
        independent_oracle::output(dir,"floatap_pow.tsv",generate(),write);
        return 0;
    } catch(std::exception const& ex) {
        std::cerr<<"FloatAP pow generator: "<<ex.what()<<'\n';
        return 1;
    }
}
