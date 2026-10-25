// SPDX-License-Identifier: MIT
// Reproduce FloatAP directed rational bounds without external arithmetic libraries.
// Standalone C++20; no Numeric, GMP, independent exact arithmetic, Boost or FloatAP.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;

std::size_t bit_length(Nat const& n) {
    return n.zero() ? 0u : (n.w.size()-1u)*64u +
        64u-static_cast<std::size_t>(std::countl_zero(n.w.back()));
}
std::pair<Nat,Nat> divrem(Nat remainder,Nat const& denominator) {
    if(denominator.zero()) throw std::invalid_argument("zero denominator");
    Nat quotient;
    while(remainder.compare(denominator)>=0) {
        std::size_t shift=bit_length(remainder)-bit_length(denominator);
        Nat scaled=denominator.shl(shift);
        if(scaled.compare(remainder)>0)scaled=denominator.shl(--shift);
        remainder=Nat::sub(remainder,scaled);
        quotient=Nat::add(quotient,Nat::power2(shift));
    }
    return {quotient,remainder};
}
Nat exact_divide(Nat const& n,Nat const& d) {
    const auto [q,r]=divrem(n,d);
    if(!r.zero())throw std::runtime_error("non-exact canonical fraction");
    return q;
}
struct Fraction {
    Nat numerator;
    Nat denominator;
    bool negative;
    Fraction(Nat n,Nat d,bool sign=false):numerator(std::move(n)),denominator(std::move(d)),negative(sign&&!numerator.zero()) {
        const Nat gcd=Nat::gcd(numerator,denominator);
        numerator=exact_divide(numerator,gcd);
        denominator=exact_divide(denominator,gcd);
    }
    std::string num() const {return (negative?"-":"")+numerator.dec();}
    std::string den() const {return denominator.dec();}
};
std::pair<Fraction,Fraction> bounds(std::uint64_t raw_n,std::uint64_t raw_d,unsigned precision,bool negative) {
    Nat n(raw_n),d(raw_d);
    if(n.zero()) throw std::invalid_argument("zero has a unique exact bound");
    std::int64_t exponent=static_cast<std::int64_t>(bit_length(n))-static_cast<std::int64_t>(bit_length(d));
    if((exponent>=0 ? n.compare(d.shl(static_cast<std::size_t>(exponent)))<0
                    : n.shl(static_cast<std::size_t>(-exponent)).compare(d)<0))--exponent;
    const std::int64_t unit_exp=exponent-static_cast<std::int64_t>(precision)+1;
    Nat scaled_n(raw_n),scaled_d(raw_d);
    if(unit_exp<0)scaled_n=n.shl(static_cast<std::size_t>(-unit_exp));
    else scaled_d=d.shl(static_cast<std::size_t>(unit_exp));
    auto [q,r]=divrem(scaled_n,scaled_d);
    Nat hi=r.zero()?q:Nat::add(q,Nat(1u));
    Nat denominator(1u);
    if(unit_exp<0)denominator=denominator.shl(static_cast<std::size_t>(-unit_exp));
    else {q=q.shl(static_cast<std::size_t>(unit_exp));hi=hi.shl(static_cast<std::size_t>(unit_exp));}
    if(negative)return {Fraction(hi,denominator,true),Fraction(q,denominator,true)};
    return {Fraction(q,denominator),Fraction(hi,denominator)};
}
std::string generate() {
    std::string result=
        "# Independent FloatAP directed bound oracle, integer arithmetic only.\n"
        "# bits numerator denominator lower_num lower_den upper_num upper_den\n"
        "# For each exact rational q and precision p, find k=floor(log2(abs(q)))\n"
        "# by integer bit-length and shift comparisons; let e=k-p+1.\n"
        "# lower=floor(q/2^e)*2^e; upper=ceil(q/2^e)*2^e.\n"
        "# Frozen as reduced rationals. No FloatAP, independent exact arithmetic, independent exact arithmetic or GMP involved.\n"
        "# Six precisions x fifteen signed, zero and dyadic inputs = 90 rows.\n";
    constexpr std::array<std::pair<std::int64_t,std::uint64_t>,15u> inputs={{{0,1},
        {1,3},{-1,3},{5,3},{-5,3},{7,11},{-7,11},{1,2},{-1,2},
        {9,8},{-9,8},{123456789,1000000007},{-123456789,1000000007},{17,5},{-17,5}}};
    for(unsigned precision:{7u,32u,53u,64u,128u,192u}) {
        for(const auto& [num,den]:inputs) {
            const std::uint64_t magnitude=static_cast<std::uint64_t>(num<0?-num:num);
            if(magnitude==0u) {
                result+=row({std::to_string(precision),"0","1","0","1","0","1"});
                continue;
            }
            const auto [lower,upper]=bounds(magnitude,den,precision,num<0);
            result+=row({std::to_string(precision),std::to_string(num),std::to_string(den),
                lower.num(),lower.den(),upper.num(),upper.den()});
        }
    }
    return result;
}
}
int main(int argc,char**argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"floatap_directed_bounds.tsv",generate(),write);
        return 0;
    }catch(std::exception const& e){
        std::cerr<<"FloatAP directed bounds generator: "<<e.what()<<'\n';
        return 1;
    }
}
