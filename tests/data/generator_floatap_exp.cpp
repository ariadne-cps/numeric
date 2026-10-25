// SPDX-License-Identifier: MIT
// Rigorous integer-only exponential oracle for the original FloatAP fixtures.
// Standalone standard C++20; no GMP, independent exact arithmetic, Numeric, Boost, or host floats.
#include "oracle_exact_integer.hpp"
#include <array>
#include <bit>
#include <cstdint>
#include <stdexcept>

namespace {
using independent_oracle::Nat;
using independent_oracle::row;
constexpr std::size_t work_bits=1400u;
std::size_t bits(Nat const& n) {
    return n.zero()?0u:(n.w.size()-1u)*64u+64u-
        static_cast<std::size_t>(std::countl_zero(n.w.back()));
}
Nat multiply(Nat const& a,Nat const& b) {
    if(a.zero()||b.zero())return {};
    const auto unpack=[](Nat const& n) {
        std::vector<std::uint32_t> out;
        out.reserve(n.w.size()*2u);
        for(auto v:n.w) {
            out.push_back(static_cast<std::uint32_t>(v));
            out.push_back(static_cast<std::uint32_t>(v>>32u));
        }
        while(!out.empty()&&out.back()==0u)out.pop_back();
        return out;
    };
    const auto x=unpack(a),y=unpack(b);
    std::vector<std::uint32_t> p(x.size()+y.size(),0u);
    for(std::size_t i=0;i<x.size();++i){
        std::uint64_t carry=0u;
        for(std::size_t j=0;j<y.size();++j){
            const auto v=std::uint64_t{x[i]}*y[j]+p[i+j]+carry;
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
std::pair<Nat,Nat> divrem(Nat n,Nat const& d){
    if(d.zero())throw std::invalid_argument("zero denominator");
    Nat q;
    while(n.compare(d)>=0){
        std::size_t shift=bits(n)-bits(d);
        Nat trial=d.shl(shift);
        if(trial.compare(n)>0)trial=d.shl(--shift);
        n=Nat::sub(n,trial);
        q=Nat::add(q,Nat::power2(shift));
    }
    return {std::move(q),std::move(n)};
}
Nat ceildiv(Nat const& n,Nat const& d){
    auto [q,r]=divrem(n,d);
    return r.zero()?q:Nat::add(q,Nat(1u));
}
struct Interval{Nat lo,hi;};
Nat const& B(){static Nat b=Nat::power2(work_bits);return b;}
Interval integer(std::uint64_t n){const Nat fixed=Nat(n).shl(work_bits);return {fixed,fixed};}
Interval ratio(std::uint64_t n,Nat const& d){
    const Nat scaled=Nat(n).shl(work_bits);
    auto [q,r]=divrem(scaled,d);
    return {q,r.zero()?q:Nat::add(q,Nat(1u))};
}
Interval add(Interval const& a,Interval const& b){
    return {Nat::add(a.lo,b.lo),Nat::add(a.hi,b.hi)};
}
Interval mul(Interval const& a,Interval const& b){
    return {multiply(a.lo,b.lo).shr(work_bits),
        ceildiv(multiply(a.hi,b.hi),B())};
}
Interval sqr(Interval const& a){
    return {multiply(a.lo,a.lo).shr(work_bits),
        ceildiv(multiply(a.hi,a.hi),B())};
}
Interval divpositive(Interval const& a,Interval const& b){
    if(b.lo.zero())throw std::invalid_argument("division interval encloses zero");
    return {divrem(a.lo.shl(work_bits),b.hi).first,
        ceildiv(a.hi.shl(work_bits),b.lo)};
}
Interval exp_abs(std::uint64_t n,Nat const& d){
    const auto z=ratio(n,d.shl(4u));
    Interval term=integer(1u),sum=term;
    for(unsigned k=1u;k<=165u;++k){
        term=divpositive(mul(term,z),integer(k));
        sum=add(sum,term);
    }
    const auto next=divpositive(mul(term,z),integer(166u));
    sum.hi=Nat::add(sum.hi,next.hi.mul_small(2u));
    for(unsigned j=0u;j<4u;++j)sum=sqr(sum);
    return sum;
}
std::pair<std::string,std::string> fraction(Nat magnitude,std::int64_t exponent){
    if(magnitude.zero())return {"0","1"};
    if(exponent<0){
        const auto cut=std::min(static_cast<std::size_t>(-exponent),magnitude.trailing_zero_bits());
        magnitude=magnitude.shr(cut);
        exponent+=static_cast<std::int64_t>(cut);
    }
    if(exponent>=0)magnitude=magnitude.shl(static_cast<std::size_t>(exponent));
    return {magnitude.dec(),exponent<0?Nat::power2(static_cast<std::size_t>(-exponent)).dec():"1"};
}
std::pair<std::string,std::string> round(Nat const& fixed,unsigned precision,char mode){
    if(fixed.zero())return {"0","1"};
    const std::size_t nbits=bits(fixed);
    const auto grid=static_cast<std::int64_t>(nbits)-
        static_cast<std::int64_t>(work_bits)-static_cast<std::int64_t>(precision);
    const auto shift=static_cast<std::int64_t>(nbits)-static_cast<std::int64_t>(precision);
    if(shift<=0)return fraction(fixed.shl(static_cast<std::size_t>(-shift)),grid);
    const auto amount=static_cast<std::size_t>(shift);
    Nat q=fixed.shr(amount);
    const Nat r=Nat::sub(fixed,q.shl(amount));
    bool increment=false;
    if(mode=='N'){
        const int c=r.shl(1u).compare(Nat::power2(amount));
        increment=c>0||(c==0&&!q.zero()&&(q.w.front()&1u)!=0u);
    }else if(mode=='U')increment=!r.zero();
    else if(mode!='D'&&mode!='Z')throw std::invalid_argument("mode");
    if(increment)q=Nat::add(q,Nat(1u));
    return fraction(std::move(q),grid);
}
std::pair<std::string,std::string> certify(Interval const& interval,
                                          unsigned p,char mode){
    const auto low=round(interval.lo,p,mode);
    if(low!=round(interval.hi,p,mode))
        throw std::runtime_error("exp interval not certified at precision " +
                                 std::to_string(p) + " mode " + mode);
    return low;
}
std::string generate(){
    std::string out=
        "# Certified exp(x) correctly rounded results; independent of GMP/independent exact arithmetic/FloatAP.\n"
        "# precision_bits mode input_num input_den expected_num expected_den\n"
        "# W=1400, B=2^W, t=abs(x)/16 <= 7/16, Taylor N=165.\n"
        "# Every term t^k/k! bounded with outward floor/ceil fixed-point division.\n"
        "# Omitted positive tail bounded by twice the upper estimate of term N+1.\n"
        "# Both bounds squared four times, outward rounded; for negative x invert.\n"
        "# Every stored result is the identical rational binary rounding of both endpoints.\n"
        "# N nearest-even, D downward, U upward, Z toward zero.\n"
        "# 18 inputs, 4 bit widths, 4 rounding modes: 288 certified rows.\n";
    const std::array<std::pair<std::uint64_t,Nat>,9u> inputs={{{1,Nat(4u)},{3,Nat(4u)},{1,Nat(2u)},{1,Nat(1u)},{3,Nat(1u)},{7,Nat(1u)},{1,Nat::power2(16u)},{1,Nat::power2(64u)},{1023,Nat(1024u)}}};
    for(const auto& [n,d]:inputs){
        for(bool negative:{false,true}){
            Interval bounds=exp_abs(n,d);
            if(negative)bounds=divpositive(integer(1u),bounds);
            for(unsigned p:{53u,128u,256u,1024u})
                for(char mode:{'N','D','U','Z'}){
                    const auto expected=certify(bounds,p,mode);
                    out+=row({std::to_string(p),std::string(1u,mode),
                         (negative?"-":"")+std::to_string(n),d.dec(),
                         expected.first,expected.second});
                }
        }
    }
    return out;
}
}
int main(int argc,char**argv){
    try{
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"floatap_exp.tsv",generate(),write);
        return 0;
    }catch(std::exception const& e){
        std::cerr<<"FloatAP exp generator: "<<e.what()<<'\n';
        return 1;
    }
}
