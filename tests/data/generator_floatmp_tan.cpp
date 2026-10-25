// SPDX-License-Identifier: MIT
// Independently certified tangent fixture, standard C++20 only.
// No Numeric, GMP, independent exact arithmetic, independent exact arithmetic, Boost or host floating-point arithmetic.
#include "oracle_exact_integer.hpp"
#include <array>
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
    // 32-bit schoolbook products avoid non-standard wide integer types.
    const auto unpack=[](Nat const& n){
        std::vector<std::uint32_t> d;
        d.reserve(n.w.size()*2u);
        for(auto word:n.w) {
            d.push_back(static_cast<std::uint32_t>(word));
            d.push_back(static_cast<std::uint32_t>(word>>32u));
        }
        while(!d.empty()&&d.back()==0u)d.pop_back();
        return d;
    };
    const auto x=unpack(a),y=unpack(b);
    std::vector<std::uint32_t> p(x.size()+y.size(),0u);
    for(std::size_t i=0;i<x.size();++i) {
        std::uint64_t carry=0u;
        for(std::size_t j=0;j<y.size();++j) {
            const std::uint64_t v=std::uint64_t{x[i]}*y[j]+p[i+j]+carry;
            p[i+j]=static_cast<std::uint32_t>(v);
            carry=v>>32u;
        }
        std::size_t k=i+y.size();
        while(carry) {
            if(k==p.size())p.push_back(0u);
            const std::uint64_t v=std::uint64_t{p[k]}+carry;
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
std::pair<Nat,Nat> divrem(Nat numerator,Nat const& d) {
    if(d.zero())throw std::invalid_argument("zero divisor");
    Nat quotient;
    while(numerator.compare(d)>=0) {
        std::size_t sh=bits(numerator)-bits(d);
        Nat trial=d.shl(sh);
        if(trial.compare(numerator)>0)trial=d.shl(--sh);
        numerator=Nat::sub(numerator,trial);
        quotient=Nat::add(quotient,Nat::power2(sh));
    }
    return {std::move(quotient),std::move(numerator)};
}
Nat ceildiv(Nat const& n,Nat const& d) {
    auto [q,r]=divrem(n,d);
    return r.zero()?q:Nat::add(q,Nat(1u));
}
struct Interval{Nat lo,hi;};
Nat const& B() {static const Nat b=Nat::power2(work_bits);return b;}
Interval zero() {return {{},{}};}
Interval one() {return {B(),B()};}
Interval ratio(Nat const& n,Nat const& d) {
    auto [q,r]=divrem(n.shl(work_bits),d);
    return {q,r.zero()?q:Nat::add(q,Nat(1u))};
}
Interval add(Interval const& a,Interval const& b) {
    return {Nat::add(a.lo,b.lo),Nat::add(a.hi,b.hi)};
}
Interval mul(Interval const& a,Interval const& b) {
    return {multiply(a.lo,b.lo).shr(work_bits),
            ceildiv(multiply(a.hi,b.hi),B())};
}
Interval divby(Interval const& a,std::uint64_t d) {
    const Nat den(d);
    return {divrem(a.lo,den).first,ceildiv(a.hi,den)};
}
struct Signed {Nat mag;bool negative=false;};
Signed signed_sub(Nat const& a,Nat const& b) {
    const int c=a.compare(b);
    if(c==0)return {};
    return c>0?Signed{Nat::sub(a,b),false}:Signed{Nat::sub(b,a),true};
}
Signed signed_add(Signed a,Nat const& positive) {
    if(!a.negative)return {Nat::add(a.mag,positive),false};
    return signed_sub(positive,a.mag);
}
Signed signed_minus(Signed a,Nat const& positive) {
    if(a.negative)return {Nat::add(a.mag,positive),true};
    return signed_sub(a.mag,positive);
}
struct SignedInterval {Signed lo,hi;};
SignedInterval alternating_series(Interval const& x,bool cosine) {
    // At |x|<=7/2, the alternating Taylor tails after 240 terms are
    // bounded by the first omitted term (terms decrease beyond k=240).
    Interval term=cosine?one():x;
    const Interval square=mul(x,x);
    Interval positive=zero(),negative=zero();
    for(unsigned k=0u;k<240u;++k) {
        if(k%2u==0u)positive=add(positive,term);
        else negative=add(negative,term);
        const std::uint64_t first=cosine?2u*k+1u:2u*k+2u;
        term=divby(mul(term,square),first*(first+1u));
    }
    // After k=239, the k=240 term is positive. Include a symmetric
    // remainder for robust enclosure, even at hard-to-round values.
    Signed lo=signed_sub(positive.lo,negative.hi);
    Signed hi=signed_sub(positive.hi,negative.lo);
    lo=signed_minus(std::move(lo),term.hi);
    hi=signed_add(std::move(hi),term.hi);
    return {std::move(lo),std::move(hi)};
}
SignedInterval negate(SignedInterval value) {
    std::swap(value.lo,value.hi);
    if(!value.lo.mag.zero())value.lo.negative=!value.lo.negative;
    if(!value.hi.mag.zero())value.hi.negative=!value.hi.negative;
    return value;
}
struct Dyadic {Nat mantissa;std::int64_t exponent=0;bool negative=false;};
std::int64_t binary_exponent(Nat const& n,Nat const& d) {
    std::int64_t e=static_cast<std::int64_t>(bits(n))-
        static_cast<std::int64_t>(bits(d));
    if(e>=0) {if(n.compare(d.shl(static_cast<std::size_t>(e)))<0)--e;}
    else if(n.shl(static_cast<std::size_t>(-e)).compare(d)<0)--e;
    return e;
}
Nat round_mag(Nat q,Nat const& rem,Nat const& d,char mode,bool negative) {
    if(rem.zero())return q;
    bool inc=false;
    if(mode=='N') {
        const int cmp=rem.shl(1u).compare(d);
        inc=cmp>0||(cmp==0&&!q.zero()&&(q.w.front()&1u)!=0u);
    }else if(mode=='D')inc=negative;
    else if(mode=='U')inc=!negative;
    else if(mode!='Z')throw std::invalid_argument("bad rounding mode");
    return inc?Nat::add(q,Nat(1u)):q;
}
Dyadic nearest(Nat const& n,Nat const& d,unsigned precision) {
    if(n.zero())return {};
    const std::int64_t grid=binary_exponent(n,d)-static_cast<std::int64_t>(precision)+1;
    const Nat a=grid<0?n.shl(static_cast<std::size_t>(-grid)):n;
    const Nat b=grid>=0?d.shl(static_cast<std::size_t>(grid)):d;
    auto [q,r]=divrem(a,b);
    return {round_mag(std::move(q),r,b,'N',false),grid,false};
}
std::pair<std::string,std::string> fraction(Nat mag,std::int64_t e,bool negative=false) {
    if(mag.zero())return {"0","1"};
    if(e<0) {
        const auto cut=std::min(static_cast<std::size_t>(-e),mag.trailing_zero_bits());
        mag=mag.shr(cut);
        e+=static_cast<std::int64_t>(cut);
    }
    if(e>=0)mag=mag.shl(static_cast<std::size_t>(e));
    return {(negative?"-":"")+mag.dec(),
            e<0?Nat::power2(static_cast<std::size_t>(-e)).dec():"1"};
}
std::pair<std::string,std::string> exact_key(Nat const& n,Nat const& d,bool negative) {
    // All input sources are binary rationals; cancel common powers of two.
    if(n.zero())return {"0","1"};
    if(d.zero()||d.w.size()>30u)throw std::invalid_argument("invalid source den");
    const auto shift=d.trailing_zero_bits();
    if(Nat::power2(shift).compare(d)!=0)throw std::invalid_argument("non-dyadic input");
    return fraction(n,-static_cast<std::int64_t>(shift),negative);
}
std::pair<std::string,std::string> round_fixed(Signed const& value,unsigned p,char mode) {
    if(value.mag.zero())return {"0","1"};
    const std::int64_t grid=static_cast<std::int64_t>(bits(value.mag))-
        static_cast<std::int64_t>(work_bits)-static_cast<std::int64_t>(p);
    const std::int64_t shift=static_cast<std::int64_t>(bits(value.mag))-
        static_cast<std::int64_t>(p);
    if(shift<=0)
        return fraction(value.mag.shl(static_cast<std::size_t>(-shift)),grid,value.negative);
    const auto amount=static_cast<std::size_t>(shift);
    Nat q=value.mag.shr(amount);
    const Nat rem=Nat::sub(value.mag,q.shl(amount));
    q=round_mag(std::move(q),rem,Nat::power2(amount),mode,value.negative);
    return fraction(std::move(q),grid,value.negative);
}
std::pair<std::string,std::string> certify(SignedInterval const& value,unsigned p,char mode) {
    const auto l=round_fixed(value.lo,p,mode);
    const auto h=round_fixed(value.hi,p,mode);
    if(l!=h)throw std::runtime_error("unresolved sin/cos interval at p="+
                                      std::to_string(p)+" mode="+mode);
    return l;
}


int compare_signed(Signed const& a,Signed const& b) {
    if(a.negative!=b.negative)return a.negative?-1:1;
    const int c=a.mag.compare(b.mag);
    return a.negative?-c:c;
}
Signed floor_fraction(Signed const& a,Nat const& d) {
    auto [q,r]=divrem(a.mag.shl(work_bits),d);
    if(a.negative&&!r.zero())q=Nat::add(q,Nat(1u));
    return {std::move(q),a.negative};
}
Signed ceil_fraction(Signed const& a,Nat const& d) {
    auto [q,r]=divrem(a.mag.shl(work_bits),d);
    if(!a.negative&&!r.zero())q=Nat::add(q,Nat(1u));
    return {std::move(q),a.negative};
}
SignedInterval divide_intervals(SignedInterval numerator,SignedInterval denominator) {
    // cos is strictly nonzero throughout each enclosure, even for |x|>pi/2.
    if(denominator.lo.mag.zero()||denominator.hi.mag.zero() ||
       denominator.lo.negative!=denominator.hi.negative)
        throw std::runtime_error("cos interval contains zero");
    if(denominator.lo.negative) {
        numerator=negate(numerator);
        denominator=negate(denominator);
    }
    const std::array<Signed,4u> n={numerator.lo,numerator.lo,
                                  numerator.hi,numerator.hi};
    const std::array<Nat,4u> d={denominator.lo.mag,denominator.hi.mag,
                                denominator.lo.mag,denominator.hi.mag};
    Signed low=floor_fraction(n[0],d[0]);
    Signed high=ceil_fraction(n[0],d[0]);
    for(std::size_t i=1;i<n.size();++i) {
        const auto lower=floor_fraction(n[i],d[i]);
        const auto upper=ceil_fraction(n[i],d[i]);
        if(compare_signed(lower,low)<0)low=lower;
        if(compare_signed(upper,high)>0)high=upper;
    }
    return {std::move(low),std::move(high)};
}
struct Generator {
    std::string out=
        "# D92: independent correctly rounded tangent golden references.\n"
        "# bits mode original_input_num original_input_den expected_num expected_den\n"
        "# Original rational argument is the key, independently converted nearest/even.\n"
        "# B=2^1800, outward fixed-point alternating sine/cosine Taylor bounds\n"
        "# k=0..239 plus next positive omitted term; cos lower bound strictly >0.\n"
        "# Divide sin enclosure by cos enclosure with outward floor/ceil, apply oddness.\n"
        "# Emit only when the bounds round to identical exact binary rationals.\n"
        "# No GMP, independent exact arithmetic, independent exact arithmetic, FloatMP or floating transcendental oracle.\n";
    std::set<std::tuple<unsigned,char,std::string,std::string>> seen;
    std::map<std::pair<std::string,std::string>,SignedInterval> cache;
    std::size_t count=0;
    void emit(unsigned p,Nat const& n,Nat const& d,bool negative) {
        const auto [input_num,input_den]=exact_key(n,d,negative);
        bool needed=false;
        for(char mode:{'N','D','U','Z'})
            if(!seen.contains({p,mode,input_num,input_den}))needed=true;
        if(!needed)return;
        const Dyadic rounded=nearest(n,d,p);
        const auto [rounded_num,rounded_den]=fraction(rounded.mantissa,rounded.exponent);
        const auto key=std::make_pair(rounded_num,rounded_den);
        auto it=cache.find(key);
        if(it==cache.end()) {
            Nat actual_n=rounded.mantissa,actual_d(1u);
            if(rounded.exponent>=0)actual_n=actual_n.shl(static_cast<std::size_t>(rounded.exponent));
            else actual_d=actual_d.shl(static_cast<std::size_t>(-rounded.exponent));
            const Interval x=ratio(actual_n,actual_d);
            const auto sine=alternating_series(x,false);
            const auto cosine=alternating_series(x,true);
            it=cache.emplace(key,divide_intervals(sine,cosine)).first;
        }
        const auto tangent=negative?negate(it->second):it->second;
        for(char mode:{'N','D','U','Z'}) {
            if(!seen.emplace(p,mode,input_num,input_den).second)continue;
            const auto expected=certify(tangent,p,mode);
            out+=row({std::to_string(p),std::string(1u,mode),input_num,input_den,
                     expected.first,expected.second});
            ++count;
        }
    }
    void half_offset(unsigned p,unsigned power,std::uint64_t distance,bool above,bool negative) {
        // 1/2 +/- distance * 2^-power
        Nat n=Nat::power2(power-1u);
        const Nat delta(distance);
        n=above?Nat::add(n,delta):Nat::sub(n,delta);
        emit(p,n,Nat::power2(power),negative);
    }
    void generate() {
        // D92: local sin/cos/tan original visitor order.
        for(unsigned p:{53u,7u,128u,192u,224u,256u,384u,1024u,53u})
            for(unsigned i=0;i<9u;++i)
                for(bool negative:{true,false}) {
                    if(i<=4u)half_offset(p,p,i==0u?0u:i<=2u?1u:8u,i==0u||i==2u||i==4u,negative);
                    if(i==5u)emit(p,Nat(3u),Nat(8u),negative);
                    if(i==6u)emit(p,Nat(5u),Nat(8u),negative);
                    if(i==7u)emit(p,Nat(1u),Nat(4u),negative);
                    if(i==8u)emit(p,Nat(3u),Nat(4u),negative);
                }
        if(count!=576u)throw std::runtime_error("D92 row count " + std::to_string(count));
        out+="# D93: centered_tan exact cutoff and fallback input probes.\n";
        // D93: centered tan at the p/2+16 cutoff.
        for(unsigned p:{53u,7u,256u,128u,1024u,53u})
            for(unsigned i=0;i<9u;++i)
                for(bool negative:{true,false}) {
                    const unsigned shift=p/2u+16u;
                    if(i<=6u)half_offset(p,shift,i==0u?0u:i<=2u?1u:i<=4u?2u:8u,
                        i==0u||i==2u||i==4u||i==6u,negative);
                    if(i==7u)emit(p,Nat(1u),Nat(4u),negative);
                    if(i==8u)emit(p,Nat(3u),Nat(4u),negative);
                }
        if(count!=816u)throw std::runtime_error("D93 row count " + std::to_string(count));
        out+="# D94: single-sqrt cosine tangent regression (including x>1 and binary-rounded p=7 sources).\n";
        constexpr std::array<std::pair<std::uint64_t,unsigned>,7u> d94={{{1u,64u},{1u,1u},
            {3u,2u},{1023u,10u},{1u,0u},{1025u,10u},{3u,1u}}};
        for(unsigned p:{7u,53u,128u,256u,1024u})
            for(const auto& [numerator,shift]:d94)
                for(bool negative:{false,true})
                    emit(p,Nat(numerator),Nat::power2(shift),negative);
        if(count!=1016u)throw std::runtime_error("D94 row count " + std::to_string(count));
        out+="# D95: exact dyadic trigonometric recurrence and reduction probes.\n";
        constexpr std::array<std::pair<std::uint64_t,unsigned>,11u> d95={{{1u,2u},{3u,2u},
            {5u,1u},{7u,1u},{2u,1u},{7u,3u},{1023u,10u},{31u,5u},
            {17u,4u},{1u,64u},{1u,128u}}};
        for(unsigned p:{53u,128u,256u,1024u})
            for(const auto& [numerator,shift]:d95)
                for(bool negative:{false,true})
                    emit(p,Nat(numerator),Nat::power2(shift),negative);
        if(count!=1208u)throw std::runtime_error("D95 row count " + std::to_string(count));
    }
};
}
int main(int argc,char**argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        Generator generator;
        generator.generate();
        independent_oracle::output(dir,"floatmp_tan.tsv",generator.out,write);
        return 0;
    }catch(std::exception const& e) {
        std::cerr<<"FloatMP tan generator: "<<e.what()<<'\n';
        return 1;
    }
}
