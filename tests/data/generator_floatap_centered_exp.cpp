// SPDX-License-Identifier: MIT
// Rigorously certified original centered FloatAP exponential fixture.
// Standalone standard C++20; no GMP, independent exact arithmetic, Numeric, Boost, or host floats.
#include "oracle_exact_integer.hpp"
#include <array>
#include <bit>
#include <cstdint>
#include <stdexcept>
#include <map>
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
Interval ratio(Nat const& n,Nat const& d){
    const Nat scaled=n.shl(work_bits);
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
Interval exp_abs(Nat const& n,Nat const& d){
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
std::int64_t binary_exponent(Nat const& n,Nat const& d) {
    std::int64_t e=static_cast<std::int64_t>(bits(n))-
        static_cast<std::int64_t>(bits(d));
    if(e>=0) {
        if(n.compare(d.shl(static_cast<std::size_t>(e)))<0)--e;
    } else if(n.shl(static_cast<std::size_t>(-e)).compare(d)<0)--e;
    return e;
}
std::pair<Nat,Nat> nearest_binary(Nat const& n,Nat const& d,unsigned precision) {
    if(n.zero())return {Nat(),Nat(1u)};
    const std::int64_t grid=binary_exponent(n,d)-static_cast<std::int64_t>(precision)+1;
    const Nat numerator=grid<0?n.shl(static_cast<std::size_t>(-grid)):n;
    const Nat denominator=grid>=0?d.shl(static_cast<std::size_t>(grid)):d;
    auto [q,r]=divrem(numerator,denominator);
    const int c=r.shl(1u).compare(denominator);
    if(c>0||(c==0&&!q.zero()&&(q.w.front()&1u)!=0u))q=Nat::add(q,Nat(1u));
    Nat output_num=q,output_den(1u);
    if(grid<0)output_den=Nat::power2(static_cast<std::size_t>(-grid));
    else output_num=output_num.shl(static_cast<std::size_t>(grid));
    if(!output_num.zero()) {
        const auto cancel=std::min(output_num.trailing_zero_bits(),output_den.trailing_zero_bits());
        output_num=output_num.shr(cancel);
        output_den=output_den.shr(cancel);
    }
    return {std::move(output_num),std::move(output_den)};
}
std::pair<Nat,Nat> half_offset(unsigned shift,std::uint64_t distance,bool above) {
    // 1/2 +/- distance / 2^shift, with shift >= 2.
    Nat n=Nat::power2(shift-1u);
    const Nat offset(distance);
    n=above?Nat::add(n,offset):Nat::sub(n,offset);
    Nat d=Nat::power2(shift);
    if(!n.zero()) {
        const auto cancel=std::min(n.trailing_zero_bits(),d.trailing_zero_bits());
        n=n.shr(cancel);
        d=d.shr(cancel);
    }
    return {std::move(n),std::move(d)};
}
std::string generate(){
    std::string out=
        "# Certified correctly rounded FloatAP centered-exponential oracle, D85.\n"
        "# precision_bits mode rounded_input_num rounded_input_den expected_num expected_den\n"
        "# Original tests: 6 visits (53,7,128,256,1024,53) x 15 magnitudes x 2 signs x 4 modes = 720 checks.\n"
        "# 568 unique references: 142 distinct rounded signed inputs x 4 rounding modes.\n"
        "# Original source Rational is rounded independently to nearest even at each precision.\n"
        "# Exact certificate: t=abs(x)/16 <= 1/8, B=2^1800, outward fixed-point Taylor k=0..165.\n"
        "# The positive remainder is at most twice the upper bound for term k=166.\n"
        "# Four outward-rounded squarings enclose exp(abs(x)); reciprocal bounds for x<0.\n"
        "# Both rational interval endpoints must round identically in the indicated mode.\n"
        "# N: nearest even; D: downward; U: upward; Z: toward zero.\n"
        "# No GMP, independent exact arithmetic, FloatAP, independent exact arithmetic or floating-point exponential oracle used.\n";
    std::set<std::tuple<unsigned,char,std::string,std::string>> seen;
    std::map<std::pair<std::string,std::string>,Interval> intervals;
    std::size_t count=0;
    for(unsigned p:{53u,7u,128u,256u,1024u,53u}) {
        const unsigned cutoff=(p+33u)/2u;
        const std::array<std::pair<Nat,Nat>,15u> sources={{
            {Nat(1u),Nat(2u)},
            half_offset(p,1u,false),half_offset(p,1u,true),
            half_offset(p,4u,false),half_offset(p,4u,true),
            half_offset(cutoff,1u,false),half_offset(cutoff,1u,true),
            half_offset(cutoff,2u,false),half_offset(cutoff,2u,true),
            {Nat(1u),Nat(4u)},{Nat(3u),Nat(4u)},
            {Nat(1u),Nat(3u)},{Nat(2u),Nat(3u)},
            {Nat(1u),Nat(1u)},{Nat(2u),Nat(1u)}
        }};
        for(auto const& [n,d]:sources) {
            const auto [rounded_n,rounded_d]=nearest_binary(n,d,p);
            const std::string n_string=rounded_n.dec(),d_string=rounded_d.dec();
            const auto cache_key=std::make_pair(n_string,d_string);
            for(bool negative:{true,false}) {
                const std::string signed_n=(negative?"-":"")+n_string;
                bool needed=false;
                for(char mode:{'N','D','U','Z'})
                    if(!seen.contains({p,mode,signed_n,d_string}))needed=true;
                if(!needed)continue;
                auto found=intervals.find(cache_key);
                if(found==intervals.end())
                    found=intervals.emplace(cache_key,exp_abs(rounded_n,rounded_d)).first;
                Interval bounds=found->second;
                if(negative)bounds=divpositive(integer(1u),bounds);
                for(char mode:{'N','D','U','Z'}) {
                    if(!seen.emplace(p,mode,signed_n,d_string).second)continue;
                    const auto result=certify(bounds,p,mode);
                    out+=row({std::to_string(p),std::string(1u,mode),
                              signed_n,d_string,result.first,result.second});
                    ++count;
                }
            }
        }
    }
    if(count!=568u)throw std::runtime_error("centered exp row count "+std::to_string(count));
    return out;
}
}
int main(int argc,char**argv){
    try{
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"floatap_centered_exp.tsv",generate(),write);
        return 0;
    }catch(std::exception const& e){
        std::cerr<<"FloatAP centered exp generator: "<<e.what()<<'\n';
        return 1;
    }
}
