// SPDX-License-Identifier: MIT
// Certified exact FloatMP square-root fixture, standard C++20 only.
#include "oracle_exact_integer.hpp"
namespace {
using independent_oracle::Nat;
using independent_oracle::row;
struct Dyadic { Nat mantissa; std::int64_t exponent=0; };
std::size_t bits(Nat const& a) {return a.zero()?0u:(a.w.size()-1u)*64u+64u-static_cast<std::size_t>(std::countl_zero(a.w.back()));}
std::pair<Nat,Nat> divrem(Nat remain,Nat const& den) {
    if(den.zero())throw std::invalid_argument("divide by zero");
    Nat quotient;
    while(remain.compare(den)>=0) {
        std::size_t shift=bits(remain)-bits(den);
        Nat shifted=den.shl(shift);
        if(shifted.compare(remain)>0)shifted=den.shl(--shift);
        remain=Nat::sub(remain,shifted);
        quotient=Nat::add(quotient,Nat::power2(shift));
    }
    return {std::move(quotient),std::move(remain)};
}
std::int64_t exponent(Nat const& n,Nat const& d) {
    std::int64_t value=static_cast<std::int64_t>(bits(n))-static_cast<std::int64_t>(bits(d));
    if(value>=0) {if(n.compare(d.shl(static_cast<std::size_t>(value)))<0)--value;}
    else if(n.shl(static_cast<std::size_t>(-value)).compare(d)<0)--value;
    return value;
}
Dyadic round_source(Nat n,Nat d,unsigned precision) {
    if(n.zero())throw std::invalid_argument("nonpositive sqrt input");
    const std::int64_t e=exponent(n,d)-static_cast<std::int64_t>(precision)+1;
    if(e>=0)d=d.shl(static_cast<std::size_t>(e));
    else n=n.shl(static_cast<std::size_t>(-e));
    auto [q,r]=divrem(n,d);
    const int cmp=r.shl(1u).compare(d);
    if(cmp>0||(cmp==0 && (q.w.front()&1u)!=0u))q=Nat::add(q,Nat(1u));
    return {std::move(q),e};
}
std::pair<std::string,std::string> fraction(Dyadic const& d) {
    if(d.mantissa.zero())return {"0","1"};
    Nat n=d.mantissa;
    std::int64_t e=d.exponent;
    if(e<0) {
        const auto cancel=std::min(static_cast<std::size_t>(-e),n.trailing_zero_bits());
        n=n.shr(cancel);
        e+=static_cast<std::int64_t>(cancel);
    }
    if(e>=0)n=n.shl(static_cast<std::size_t>(e));
    return {n.dec(), e<0?Nat::power2(static_cast<std::size_t>(-e)).dec():"1"};
}
Nat multiply(Nat a,Nat b) {
    Nat result;
    while(!b.zero()) {
        if((b.w.front()&1u)!=0u)result=Nat::add(result,a);
        b=b.shr(1u);
        if(!b.zero())a=a.shl(1u);
    }
    return result;
}
Nat isqrt(Nat const& n) {
    if(n.zero())return {};
    Nat root,remainder;
    const std::size_t pairs=(bits(n)+1u)/2u;
    for(std::size_t pos=pairs;pos>0;--pos) {
        const std::size_t bit=2u*(pos-1u);
        const std::uint64_t pair=(n.w[bit/64u]>>(bit%64u))&3u;
        remainder=Nat::add(remainder.shl(2u),Nat(pair));
        const Nat candidate=Nat::add(root.shl(2u),Nat(1u));
        if(remainder.compare(candidate)>=0) {
            remainder=Nat::sub(remainder,candidate);
            root=Nat::add(root.shl(1u),Nat(1u));
        } else root=root.shl(1u);
    }
    return root;
}
Dyadic square_root(Dyadic const& input,unsigned precision,char mode) {
    const std::int64_t x_log=static_cast<std::int64_t>(bits(input.mantissa))-1+input.exponent;
    const std::int64_t floor_half=x_log>=0 ? x_log/2 : -((-x_log+1)/2);
    const std::int64_t e=floor_half-static_cast<std::int64_t>(precision)+1;
    const std::int64_t shift=input.exponent-2*e;
    const Nat numerator=shift>=0?input.mantissa.shl(static_cast<std::size_t>(shift)):input.mantissa;
    const Nat denominator=shift<0?Nat::power2(static_cast<std::size_t>(-shift)):Nat(1u);
    const auto [floor_integer,_]=divrem(numerator,denominator);
    (void)_;
    Nat root=isqrt(floor_integer);
    const Nat square=multiply(root,root);
    if(multiply(square,denominator).compare(numerator)>0)throw std::runtime_error("sqrt lower invariant");
    const Nat next=Nat::add(root,Nat(1u));
    if(multiply(multiply(next,next),denominator).compare(numerator)<=0)throw std::runtime_error("sqrt upper invariant");
    const bool exact=multiply(square,denominator).compare(numerator)==0;
    bool increment=false;
    if(mode=='N') {
        const Nat mid=Nat::add(root.shl(1u),Nat(1u));
        const int cmp=numerator.shl(2u).compare(multiply(multiply(mid,mid),denominator));
        increment=cmp>0||(cmp==0 && !root.zero() && (root.w.front()&1u)!=0u);
    } else if(mode=='U')increment=!exact;
    else if(mode!='D'&&mode!='Z')throw std::invalid_argument("mode");
    if(increment)root=Nat::add(root,Nat(1u));
    return {std::move(root),e};
}
std::vector<Dyadic> original_inputs(unsigned precision) {
    std::vector<Dyadic> result;
    const auto append=[&](Nat n,Nat d) {result.push_back(round_source(std::move(n),std::move(d),precision));};
    for(auto [n,d]:std::array<std::pair<std::uint64_t,std::uint64_t>,10u>{{{1u,1u},{2u,1u},{3u,1u},{4u,1u},{5u,1u},
        {9u,1u},{15u,16u},{17u,16u},{1u,4u},{9u,16u}}})append(Nat(n),Nat(d));
    const Nat radix=Nat::power2(precision);
    const Nat one(1u);
    append(Nat::sub(radix,one),radix);
    append(Nat::add(radix,one),radix);
    append(Nat::sub(radix.mul_small(4u),one),radix);
    append(Nat::add(radix.mul_small(4u),one),radix);
    const Dyadic scaled=round_source(Nat(5u),Nat(4u),precision);
    for(std::int64_t shift:{-2000,-127,127,2000})result.push_back({scaled.mantissa,scaled.exponent+shift});
    return result;
}
std::string generate() {
    std::string out=
        "# Certified FloatMP sqrt oracle (D83), exact arithmetic only.\n"
        "# precision_bits mode input_num input_den expected_num expected_den\n"
        "# Input is the correctly nearest-even-rounded positive Rational used by the original test.\n"
        "# 5 precisions x (14 input rationals plus 4 shifted operands) x 4 modes: 360 original checks.\n"
        "# Duplicate rounded inputs are deduplicated; all original checks still execute.\n"
        "# For x=n/d>0, set e=floor(floor(log2(n/d))/2)-(p-1).\n"
        "# A/B = (n/d)*2^(-2e). Calculate q=isqrt(floor(A/B)).\n"
        "# Verify q^2*B <= A < (q+1)^2*B (integer square root floor).\n"
        "# Nearest: compare 4*A to (2*q+1)^2*B; ties-to-even. Directed: exact iff A=q^2*B.\n"
        "# Result q*2^e, canonical rational. No FloatMP, independent exact arithmetic, GMP or independent exact arithmetic required.\n";
    for(unsigned p:{7u,53u,128u,256u,1024u}) {
        std::vector<std::pair<std::string,std::string>> seen;
        for(const auto& input:original_inputs(p)) {
            const auto key=fraction(input);
            if(std::find(seen.begin(),seen.end(),key)!=seen.end())continue;
            seen.push_back(key);
            for(char mode:{'N','D','U','Z'}) {
                const auto ans=fraction(square_root(input,p,mode));
                out+=row({std::to_string(p),std::string(1u,mode),key.first,key.second,ans.first,ans.second});
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
        independent_oracle::output(dir,"floatmp_sqrt.tsv",generate(),write);
        return 0;
    } catch(std::exception const& e) {
        std::cerr<<"FloatMP sqrt generator: "<<e.what()<<'\n';
        return 1;
    }
}
