// SPDX-License-Identifier: MIT
// Independently certified sine and cosine fixtures, standard C++20 only.
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

std::string header(bool cosine) {
    std::string s=cosine?
        "# D91: independently certified correctly rounded cos values.\n":
        "# D91: independently certified correctly rounded sin values.\n";
    return s+
        "# bits mode original_input_num original_input_den expected_num expected_den\n"
        "# 7 precision visits, 7 magnitudes, both signs and 4 modes; revisit 128 bits.\n"
        "# The ORIGINAL source Rational is the key, with an independent exact nearest/even\n"
        "# binary input conversion before the trigonometric evaluation.\n"
        "# Fixed-point B=2^1800, x squared enclosed outward; 240 alternating Taylor\n"
        "# sin/cos terms with outward floor/ceil for EVERY product and division.\n"
        "# First omitted positive term bounds remainder; negative sine uses oddness.\n"
        "# Both exact interval endpoints round to the IDENTICAL p-bit binary Rational.\n"
        "# No GMP, independent exact arithmetic, independent exact arithmetic, FloatMP or floating point transcendental oracle used.\n";
}
struct Oracle {
    std::string sin_data=header(false),cos_data=header(true);
    std::set<std::tuple<unsigned,char,std::string,std::string>> seen;
    // Reuse source interval for both functions and all modes.
    std::map<std::pair<std::string,std::string>,
             std::pair<SignedInterval,SignedInterval>> cache;
    std::size_t count=0u;
    void emit(unsigned p,Nat const& numerator,Nat const& denominator,bool negative) {
        const auto [input_num,input_den]=exact_key(numerator,denominator,negative);
        bool needed=false;
        for(char mode:{'N','D','U','Z'})
            if(!seen.contains({p,mode,input_num,input_den}))needed=true;
        if(!needed)return;
        Dyadic rounded=nearest(numerator,denominator,p);
        const auto [rounded_num,rounded_den]=
            fraction(rounded.mantissa,rounded.exponent);
        const auto cache_key=std::make_pair(rounded_num,rounded_den);
        auto it=cache.find(cache_key);
        if(it==cache.end()) {
            Nat n=rounded.mantissa,d(1u);
            if(rounded.exponent>=0)n=n.shl(static_cast<std::size_t>(rounded.exponent));
            else d=d.shl(static_cast<std::size_t>(-rounded.exponent));
            const Interval x=ratio(n,d);
            const SignedInterval sine=alternating_series(x,false);
            const SignedInterval cosine=alternating_series(x,true);
            it=cache.emplace(cache_key,std::make_pair(sine,cosine)).first;
        }
        const auto sine=negative?negate(it->second.first):it->second.first;
        const auto& cosine=it->second.second;
        for(char mode:{'N','D','U','Z'}) {
            if(!seen.emplace(p,mode,input_num,input_den).second)continue;
            const auto sr=certify(sine,p,mode);
            const auto cr=certify(cosine,p,mode);
            sin_data+=row({std::to_string(p),std::string(1u,mode),input_num,input_den,
                           sr.first,sr.second});
            cos_data+=row({std::to_string(p),std::string(1u,mode),input_num,input_den,
                           cr.first,cr.second});
            ++count;
        }
    }
    void near_half(unsigned p,unsigned shift,std::int64_t offset,bool sign) {
        // Rational 1/2 + offset * 2^-shift (D91 and D92).
        const unsigned denominator_bits=std::max(shift,1u);
        Nat n=Nat::power2(denominator_bits-1u);
        const Nat adjustment=Nat(static_cast<std::uint64_t>(offset<0?-offset:offset))
                .shl(static_cast<std::size_t>(denominator_bits-shift));
        n=offset<0?Nat::sub(n,adjustment):Nat::add(n,adjustment);
        emit(p,n,Nat::power2(denominator_bits),sign);
    }
    void generate() {
        for(unsigned p:{96u,128u,192u,224u,256u,1024u,128u}) {
            const unsigned shift=p/2u+16u;
            for(unsigned i=0;i<7u;++i)
                for(bool negative:{true,false}) {
                    if(i==0u)emit(p,Nat(1u),Nat(2u),negative);
                    if(i>=1u&&i<=4u)near_half(p,shift,
                        i==1u?-1:i==2u?1:i==3u?-2:2,negative);
                    if(i==5u)emit(p,Nat(1u),Nat(4u),negative);
                    if(i==6u)emit(p,Nat(3u),Nat(4u),negative);
                }
        }
        if(count!=336u)throw std::runtime_error("D91 row count: "+std::to_string(count));
        sin_data+="# D92: locally centered sin/cos/tan regression probes.\n";
        cos_data+="# D92: locally centered sin/cos/tan regression probes.\n";
        for(unsigned p:{53u,7u,128u,192u,224u,256u,384u,1024u,53u}) {
            for(unsigned i=0;i<9u;++i)
                for(bool negative:{true,false}) {
                    if(i<=4u)near_half(p,p,i==0u?0:i==1u?-1:i==2u?1:i==3u?-8:8,negative);
                    if(i==5u)emit(p,Nat(3u),Nat(8u),negative);
                    if(i==6u)emit(p,Nat(5u),Nat(8u),negative);
                    if(i==7u)emit(p,Nat(1u),Nat(4u),negative);
                    if(i==8u)emit(p,Nat(3u),Nat(4u),negative);
                }
        }
        if(count!=792u)throw std::runtime_error("D92 row count: "+std::to_string(count));
        sin_data+="# D95: exact dyadic trigonometric recurrence and reduction probes.\n";
        cos_data+="# D95: exact dyadic trigonometric recurrence and reduction probes.\n";
        constexpr std::array<std::pair<std::uint64_t,unsigned>,11u> d95={{{1,2},{3,2},{5,1},{7,1},{2,1},{7,3},
            {1023,10},{31,5},{17,4},{1,64},{1,128}}};
        for(unsigned p:{53u,128u,256u,1024u})
            for(const auto& [numerator,shift]:d95)
                for(bool negative:{false,true})
                    emit(p,Nat(numerator),Nat::power2(shift),negative);
        if(count!=1080u)throw std::runtime_error("D95 row count: "+std::to_string(count));
    }
};
}
int main(int argc,char**argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        Oracle oracle;
        oracle.generate();
        independent_oracle::output(dir,"floatmp_sin.tsv",oracle.sin_data,write);
        independent_oracle::output(dir,"floatmp_cos.tsv",oracle.cos_data,write);
        return 0;
    }catch(std::exception const& e) {
        std::cerr<<"FloatMP sin/cos generator: "<<e.what()<<'\n';
        return 1;
    }
}
