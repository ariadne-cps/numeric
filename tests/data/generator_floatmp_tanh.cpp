// SPDX-License-Identifier: MIT
// Recreate certified tanh fixture using independent standard C++20 integer arithmetic.
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
constexpr std::size_t work_bits=2304u;
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
    const auto z=ratio(n,d.shl(5u));
    Interval term=integer(1u),sum=term;
    for(unsigned k=1u;k<=409u;++k){
        term=divpositive(mul(term,z),integer(k));
        sum=add(sum,term);
    }
    const auto next=divpositive(mul(term,z),integer(410u));
    sum.hi=Nat::add(sum.hi,next.hi.mul_small(2u));
    for(unsigned j=0u;j<6u;++j)sum=sqr(sum);
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
[[maybe_unused]] std::pair<std::string,std::string> certify(Interval const& interval,
                                          unsigned p,char mode){
    const auto low=round(interval.lo,p,mode);
    if(low!=round(interval.hi,p,mode))
        throw std::runtime_error("tanh interval not certified at precision " +
                                 std::to_string(p) + " mode " + mode);
    return low;
}

std::int64_t binary_exponent(Nat const& n,Nat const& d) {
    std::int64_t e=static_cast<std::int64_t>(bits(n))-static_cast<std::int64_t>(bits(d));
    if(e>=0) {if(n.compare(d.shl(static_cast<std::size_t>(e)))<0)--e;}
    else if(n.shl(static_cast<std::size_t>(-e)).compare(d)<0)--e;
    return e;
}
struct Rational {Nat n,d;bool negative=false;};
Rational dyadic(Nat n,unsigned denominator_bits,bool negative=false) {
    if(n.zero())return {{},Nat(1u),false};
    const auto cancel=std::min(n.trailing_zero_bits(),static_cast<std::size_t>(denominator_bits));
    n=n.shr(cancel);
    return {std::move(n),Nat::power2(denominator_bits-cancel),negative};
}
Rational fraction_of(std::uint64_t n,unsigned denom_bits=0u,bool neg=false) {
    return dyadic(Nat(n),denom_bits,neg);
}
Rational centered(unsigned p,unsigned shift,int distance,bool negative) {
    // 1/2 + distance * 2^-shift where distance is -1,0,+1.
    const unsigned power=std::max(shift,1u);
    Nat n=Nat::power2(power-1u);
    if(distance<0)n=Nat::sub(n,Nat(static_cast<std::uint64_t>(-distance)));
    else n=Nat::add(n,Nat(static_cast<std::uint64_t>(distance)));
    (void)p;
    return dyadic(std::move(n),power,negative);
}
Rational near_cutoff(unsigned p,unsigned cutoff,int side,int tiny,bool negative) {
    // 1/2 + side*2^-cutoff + tiny*2^-(p+1).
    const unsigned power=p+1u;
    Nat n=Nat::power2(power-1u);
    Nat base=Nat::power2(power-cutoff);
    if(side<0)n=Nat::sub(n,base);
    else n=Nat::add(n,base);
    if(tiny<0)n=Nat::sub(n,Nat(1u));
    else n=Nat::add(n,Nat(1u));
    return dyadic(std::move(n),power,negative);
}
int compare_rational(Rational const& a,Rational const& b) {
    if(a.negative!=b.negative)return a.negative?-1:1;
    const auto an=multiply(a.n,b.d),bn=multiply(b.n,a.d);
    const int c=an.compare(bn);
    return a.negative?-c:c;
}
Rational nearest(Rational const& source,unsigned p) {
    if(source.n.zero())return fraction_of(0u);
    const auto exponent=binary_exponent(source.n,source.d)-static_cast<std::int64_t>(p)+1;
    const Nat n=exponent<0?source.n.shl(static_cast<std::size_t>(-exponent)):source.n;
    const Nat d=exponent>=0?source.d.shl(static_cast<std::size_t>(exponent)):source.d;
    auto [q,r]=divrem(n,d);
    const int cmp=r.shl(1u).compare(d);
    if(cmp>0||(cmp==0&&!q.zero()&&(q.w.front()&1u)!=0u))q=Nat::add(q,Nat(1u));
    if(exponent<0)return dyadic(std::move(q),static_cast<unsigned>(-exponent),source.negative);
    return dyadic(q.shl(static_cast<std::size_t>(exponent)),0u,source.negative);
}
Interval tanh_abs(Rational const& x) {
    if(x.n.zero())return {Nat(),Nat()};
    if(x.n.compare(x.d.mul_small(900u))>=0)
        throw std::runtime_error("saturation must use an open upper bound");
    // E=exp(2*x)>1. tanh(x)=(E-1)/(E+1) is monotone in E.
    const auto exponential=exp_abs(x.n,x.d); // exp(|x|/32), six squarings -> exp(2|x|).
    if(exponential.lo.compare(B())<0)throw std::runtime_error("E below 1");
    const auto lo=Nat::sub(exponential.lo,B());
    const auto hi=Nat::sub(exponential.hi,B());
    return {divrem(lo.shl(work_bits),Nat::add(exponential.lo,B())).first,
            ceildiv(hi.shl(work_bits),Nat::add(exponential.hi,B()))};
}
std::pair<std::string,std::string> round_signed(Nat const& fixed,unsigned p,char mode,bool negative) {
    if(fixed.zero())return {"0","1"};
    const auto magnitude_mode=mode=='N'?'N':
        mode=='Z'?'D':(negative?(mode=='D'?'U':'D'):mode);
    auto result=round(fixed,p,magnitude_mode);
    if(negative&&result.first!="0")result.first="-"+result.first;
    return result;
}
std::pair<std::string,std::string> certify_signed(Interval const& a,unsigned p,char mode,bool negative) {
    const auto lo=round_signed(a.lo,p,mode,negative);
    const auto hi=round_signed(a.hi,p,mode,negative);
    if(lo!=hi)throw std::runtime_error("tanh interval not certified p="+std::to_string(p)+
                                      " mode="+mode);
    return lo;
}
std::pair<std::string,std::string> saturate(unsigned p,char mode,bool negative) {
    // For x >= 900: exp(2*x)>2^1800 since e>2, so
    // 0 < 1-tanh(x) = 2/(exp(2*x)+1) < 2^-1799.
    // Every source in this case has p<=1024. Strict tanh(x)<1 means
    // downward rounding is PREDECESSOR(1), not 1; nearest rounds to 1.
    const bool to_one=mode=='N'||(negative?mode=='D':mode=='U');
    const auto previous=fraction(Nat::sub(Nat::power2(p),Nat(1u)),-static_cast<std::int64_t>(p));
    const auto result=to_one?std::make_pair(std::string("1"),std::string("1")):previous;
    return {negative?"-"+result.first:result.first,result.second};
}
std::string generate() {
    std::string out=
        "# Independent correctly rounded tanh oracle.\n"
        "# bits mode input_num input_den expected_num expected_den\n"
        "# Exact input is first rounded nearest-even to p bits.\n"
        "# Generated offline using only standard C++20 unsigned integer arithmetic:\n"
        "# outward fixed-point exp Taylor recurrence on |t|<=1/4, tail <= last term,\n"
        "# outward repeated squaring, then (exp(2|x|)-1)/(exp(2|x|)+1).\n"
        "# Both rigorous endpoint bounds independently round to the same literal.\n"
        "# No Numeric, libm, GMP, independent exact arithmetic or floating-point expectation generation.\n"
        "# Four modes, both signs, near-half cutoff edges, tiny inputs and saturation.\n";
    constexpr std::array<unsigned,13u> precisions={2u,7u,53u,64u,65u,95u,96u,127u,128u,129u,256u,1024u,2048u};
    const auto append_signed=[](std::vector<Rational>& values,Rational a){
        values.push_back(a);
        if(!a.n.zero()) {a.negative=true;values.push_back(std::move(a));}
    };
    std::map<std::tuple<unsigned,std::string,std::string>,Interval> intervals;
    std::size_t rows=0;
    for(unsigned p:precisions) {
        std::vector<Rational> sources;
        const bool original=(p==2u||p==7u||p==53u||p==64u||p==128u||p==256u||p==1024u);
        if(original) {
            for(auto [n,e]:std::array<std::pair<std::uint64_t,unsigned>,8u>{{
                {1000u,0u},{8u,0u},{2u,0u},{1u,0u},
                {3u,2u},{1u,1u},{1u,3u},{1u,80u}}})
                append_signed(sources,fraction_of(n,e));
            sources.push_back(fraction_of(0u));
        }
        if(p==53u) {
            // Binary64 fixed-point rounding across a 64-bit discarded-word
            // boundary, including an even retained significand with sticky bits.
            append_signed(sources,fraction_of(1u,11u));
            append_signed(sources,fraction_of(9007199254740991ULL,64u));
        }
        if(p!=2u&&p!=7u) {
            const unsigned cutoff=original?(p+31u)/2u:(p+33u)/2u;
            append_signed(sources,centered(p,1u,0,false));
            for(int side:{-1,1}) {
                append_signed(sources,centered(p,p,side,false));
                if(original) {
                    append_signed(sources,centered(p,cutoff,side,false));
                    append_signed(sources,centered(p,cutoff+1u,side,false));
                } else {
                    append_signed(sources,centered(p,cutoff,side,false));
                    for(int d:{-1,1})
                        append_signed(sources,near_cutoff(p,cutoff,side,d,false));
                }
            }
        }
        std::sort(sources.begin(),sources.end(),[](auto const& a,auto const& b){return compare_rational(a,b)<0;});
        sources.erase(std::unique(sources.begin(),sources.end(),[](auto const& a,auto const& b){return compare_rational(a,b)==0;}),sources.end());
        const auto expected=(p==2u||p==7u)?17u:p==53u?33u:original?29u:18u;
        if(sources.size()!=expected)throw std::runtime_error("incorrect source count p="+std::to_string(p));
        for(char mode:{'D','N','U','Z'})for(auto const& source:sources) {
            const std::string input_n=(source.negative?"-":"")+source.n.dec();
            const std::string input_d=source.d.dec();
            const Rational rounded=nearest(source,p);
            auto [it,inserted]=intervals.try_emplace({p,rounded.n.dec(),rounded.d.dec()});
            const bool saturated=rounded.n.compare(rounded.d.mul_small(900u))>=0;
            if(inserted&&!saturated)it->second=tanh_abs(rounded);
            const auto result=saturated?saturate(p,mode,source.negative):
                certify_signed(it->second,p,mode,source.negative);
            out+=row({std::to_string(p),std::string(1u,mode),input_n,input_d,result.first,result.second});
            ++rows;
        }
    }
    if(rows!=1164u)throw std::runtime_error("tanh row count: "+std::to_string(rows));
    return out;
}
}
int main(int argc,char**argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"floatmp_tanh.tsv",generate(),write);
        return 0;
    }catch(std::exception const& e) {
        std::cerr<<"FloatMP tanh generator: "<<e.what()<<'\n';
        return 1;
    }
}
