// SPDX-License-Identifier: MIT
// Reproduce all exact rational FloatAP rounding fixtures from source test loops.
// Standard C++20 and independent integer limbs only, no Numeric/GMP/MPFR/Boost.
#include "oracle_exact_integer.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <set>
#include <tuple>
#include <utility>

namespace {
using independent_oracle::Nat;
using independent_oracle::row;

std::size_t bit_length(Nat const& a) {
    return a.zero()?0:(a.w.size()-1)*64+64-static_cast<std::size_t>(std::countl_zero(a.w.back()));
}
Nat product(Nat const& a,Nat const& b) {
    if(a.zero()||b.zero())return {};
    const auto unpack=[](Nat const& n) {
        std::vector<std::uint32_t> x;
        for(auto w:n.w) {
            x.push_back(static_cast<std::uint32_t>(w));
            x.push_back(static_cast<std::uint32_t>(w>>32u));
        }
        while(!x.empty()&&x.back()==0)x.pop_back();
        return x;
    };
    const auto x=unpack(a),y=unpack(b);
    std::vector<std::uint32_t> z(x.size()+y.size(),0);
    for(std::size_t i=0;i<x.size();++i) {
        std::uint64_t carry=0;
        for(std::size_t j=0;j<y.size();++j) {
            const std::uint64_t v=std::uint64_t{x[i]}*y[j]+z[i+j]+carry;
            z[i+j]=static_cast<std::uint32_t>(v);
            carry=v>>32u;
        }
        for(std::size_t k=i+y.size();carry;++k) {
            if(k>=z.size())z.push_back(0u);
            const std::uint64_t v=std::uint64_t{z[k]}+carry;
            z[k]=static_cast<std::uint32_t>(v);
            carry=v>>32u;
        }
    }
    std::vector<std::uint64_t> limbs((z.size()+1)/2,0);
    for(std::size_t k=0;k<z.size();++k)limbs[k/2]|=std::uint64_t{z[k]}<<((k%2)*32);
    return Nat(std::move(limbs));
}
std::pair<Nat,Nat> divrem(Nat a,Nat const& b) {
    if(b.zero())throw std::invalid_argument("division by zero");
    Nat q;
    while(a.compare(b)>=0) {
        std::size_t shift=bit_length(a)-bit_length(b);
        Nat trial=b.shl(shift);
        if(trial.compare(a)>0)trial=b.shl(--shift);
        a=Nat::sub(a,trial);
        q=Nat::add(q,Nat::power2(shift));
    }
    return {std::move(q),std::move(a)};
}
struct Rat { Nat n,d; bool neg=false; };
Rat reduce(Nat n,Nat d,bool negative=false) {
    if(d.zero())throw std::runtime_error("zero denominator");
    if(n.zero())return {Nat(),Nat(1),false};
    const auto binary_gcd=[](Nat x,Nat y) {
        if(x.zero())return y;
        if(y.zero())return x;
        const std::size_t common=std::min(x.trailing_zero_bits(),y.trailing_zero_bits());
        x=x.shr(x.trailing_zero_bits());
        y=y.shr(y.trailing_zero_bits());
        while(x.compare(y)!=0) {
            if(x.compare(y)>0){x=Nat::sub(x,y);x=x.shr(x.trailing_zero_bits());}
            else {y=Nat::sub(y,x);y=y.shr(y.trailing_zero_bits());}
        }
        return x.shl(common);
    };
    const Nat g=binary_gcd(n,d);
    return {divrem(n,g).first,divrem(d,g).first,negative};
}
Rat scalar(std::uint64_t n,std::uint64_t d=1,bool negative=false) {
    return reduce(Nat(n),Nat(d),negative);
}
Rat integer(std::int64_t n) {
    return scalar(static_cast<std::uint64_t>(n<0?-n:n),1,n<0);
}
Rat neg(Rat a) {
    if(!a.n.zero())a.neg=!a.neg;
    return a;
}
Rat add(Rat const& a,Rat const& b) {
    const Nat x=product(a.n,b.d),y=product(b.n,a.d),d=product(a.d,b.d);
    if(a.neg==b.neg)return reduce(Nat::add(x,y),d,a.neg);
    const int c=x.compare(y);
    return c>=0?reduce(Nat::sub(x,y),d,a.neg):reduce(Nat::sub(y,x),d,b.neg);
}
Rat sub(Rat const& a,Rat const& b) {return add(a,neg(b));}
Rat mul(Rat const& a,Rat const& b) {
    return reduce(product(a.n,b.n),product(a.d,b.d),a.neg!=b.neg);
}
Rat divide(Rat const& a,Rat const& b) {
    return reduce(product(a.n,b.d),product(a.d,b.n),a.neg!=b.neg);
}
std::int64_t binary_exponent(Nat const& n,Nat const& d) {
    std::int64_t e=static_cast<std::int64_t>(bit_length(n))-static_cast<std::int64_t>(bit_length(d));
    if(e>=0) {if(n.compare(d.shl(static_cast<std::size_t>(e)))<0)--e;}
    else if(n.shl(static_cast<std::size_t>(-e)).compare(d)<0)--e;
    return e;
}
Rat rounded(Rat const& a,unsigned precision,char mode) {
    if(a.n.zero())return scalar(0);
    const auto e=binary_exponent(a.n,a.d)-static_cast<std::int64_t>(precision)+1;
    const Nat n=e<0?a.n.shl(static_cast<std::size_t>(-e)):a.n;
    const Nat d=e>=0?a.d.shl(static_cast<std::size_t>(e)):a.d;
    auto [q,r]=divrem(n,d);
    bool increment=false;
    if(!r.zero()) {
        if(mode=='N') {
            const int cmp=r.shl(1).compare(d);
            increment=cmp>0||(cmp==0&&(q.w.front()&1u)!=0);
        } else if(mode=='D')increment=a.neg;
        else if(mode=='U')increment=!a.neg;
        else if(mode!='Z')throw std::invalid_argument("invalid rounding mode");
    }
    if(increment)q=Nat::add(q,Nat(1));
    if(e<0)return reduce(q,Nat::power2(static_cast<std::size_t>(-e)),a.neg);
    return reduce(q.shl(static_cast<std::size_t>(e)),Nat(1),a.neg);
}
std::string decimal(Rat const& a) {return (a.neg?"-":"")+a.n.dec();}
struct Generator {
    std::string data=
        "# FloatAP exact rational rounding oracle v1\n"
        "# precision_bits mode input_num input_den expected_num expected_den\n"
        "# N nearest ties-to-even; D down; U up; Z toward zero\n"
        "# Base-10 integers, positive denominators; exact rational values.\n"
        "# Derived from independent integer divrem and ties-to-even (no FloatAP, FloatMP, GMP, MPFR).\n"
        "# This is NOT a transcendental oracle. Cases 4794.\n";
    std::set<std::tuple<unsigned,char,std::string,std::string>> seen;
    std::size_t count=0;
    void emit(unsigned p,char mode,Rat const& source) {
        const auto input_n=decimal(source),input_d=source.d.dec();
        if(!seen.emplace(p,mode,input_n,input_d).second)return;
        const Rat target=rounded(source,p,mode);
        data+=row({std::to_string(p),std::string(1,mode),input_n,input_d,
                   decimal(target),target.d.dec()});
        ++count;
    }
    void modes(unsigned p,Rat const& source) {
        for(char mode:{'N','D','U','Z'})emit(p,mode,source);
    }
    void generate() {
        // Original exact-ratio constructor and 1/3 direct conversion probes.
        for(unsigned p:{2u,7u,31u,32u,33u,53u,64u,65u,127u,128u,129u,
                        255u,256u,257u,1023u,1024u,1100u,1536u}) {
            modes(p,scalar(1,3));modes(p,scalar(1,3,true));
            modes(p,reduce(Nat::add(Nat::power2(1800),Nat(7)),
                           Nat::add(Nat::power2(1500),Nat(11))));
        }
        for(Rat const& v:{scalar(5,4),scalar(7,4),scalar(5,4,true),scalar(7,4,true)})
            modes(2u,v);
        if(count!=232u)throw std::runtime_error("first section count "+std::to_string(count));
        data+="# Additional broad quotient, sparse denominator and Knuth correction boundary references\n";
        for(unsigned p:{53u,128u,256u,1024u,1536u})for(unsigned off:{1u,7u,19u,257u}) {
            const Nat power=Nat::power2(p);
            const Rat source=reduce(Nat::add(power,Nat(off)),
                Nat::sub(power.mul_small(2),Nat(off+2u)));
            modes(p,source);modes(p,neg(source));
        }
        for(unsigned p:{128u,256u,1024u,1536u})for(unsigned off:{1u,17u,257u}) {
            const Nat power=Nat::power2(p);
            const Rat source=reduce(Nat::add(power,Nat(off)),
                Nat::add(power.mul_small(2),Nat(off+2u)));
            modes(p,source);modes(p,neg(source));
        }
        for(unsigned p:{65u,128u,129u,255u,256u,257u,511u,512u,
                        1023u,1024u,1025u,1536u})for(unsigned off:{1u,3u,257u}) {
            const Nat power=Nat::power2(p);
            const Nat x=Nat::sub(power,Nat(off)),y=Nat::add(power,Nat(off+2u));
            for(bool reverse:{false,true}) {
                const Rat source=reverse?reduce(y,x):reduce(x,y);
                modes(p,source);modes(p,neg(source));
            }
        }
        if(count!=1064u)throw std::runtime_error("second section count "+std::to_string(count));
        data+="# D79: exact rational, dyadic, small-limb division and rounded integer arithmetic oracles\n";
        for(unsigned p:{7u,53u,128u,256u,1024u})
            for(int num:{-17,-3,1,17})
                for(int den:{-16,-7,-3,3,7,16}) {
                    const bool negative=(num<0)!=(den<0);
                    const Rat r=scalar(static_cast<unsigned>(num<0?-num:num),
                                       static_cast<unsigned>(den<0?-den:den),negative);
                    modes(p,r);
                }
        if(count!=1384u)throw std::runtime_error("D79 ratios count "+std::to_string(count));
        for(unsigned p:{7u,53u,128u,256u,1024u}) {
            const Nat power=Nat::power2(p);
            for(unsigned sh:{1u,7u,31u,130u}) {
                const Nat den=Nat::power2(sh);
                for(const auto& [n,negative]:std::array<std::pair<Nat,bool>,4u>{{
                    {Nat::sub(power,Nat(1)),false},{Nat::add(power,Nat(1)),false},
                    {Nat::sub(power,Nat(1)),true},{Nat::add(power,Nat(1)),true}}}) {
                    modes(p,reduce(n,den,negative));
                }
            }
        }
        if(count!=1704u)throw std::runtime_error("D79 dyadic count "+std::to_string(count));
        for(unsigned p:{7u,53u,64u,128u,256u})for(unsigned i=0;i<32u;++i) {
            const Nat base=Nat::power2(p);
            const Rat exact=reduce(Nat::add(base,Nat(13u+2u*i)),
                                   Nat::add(base.mul_small(2u),Nat(17u+3u*i)));
            const Rat divisor=reduce(Nat::add(base,Nat(31u+5u*i)),
                                     Nat::add(base.mul_small(3u),Nat(19u+7u*i)));
            for(bool negative:{true,false}) {
                const Rat source=negative?neg(exact):exact;
                const Rat input=rounded(source,p,'N');
                const Rat other=rounded(divisor,p,'N');
                const Rat quotient=divide(input,other);
                for(char mode:{'N','D','U','Z'}) {
                    emit(p,mode,source);
                    emit(p,mode,quotient);
                }
            }
        }
        if(count!=4224u)throw std::runtime_error("D79 quotient count "+std::to_string(count));
        for(unsigned p:{7u,53u,128u,256u,1024u}) {
            const Nat power=Nat::power2(p);
            for(const auto& val:std::array<Rat,6u>{{
                reduce(Nat::add(power,Nat(1)),Nat(1)),
                reduce(Nat::sub(power,Nat(1)),Nat(1)),
                reduce(Nat::add(power,Nat(1)),Nat(1),true),
                reduce(Nat::sub(power,Nat(1)),Nat(1),true),
                integer(3),integer(-3)}}) {
                for(char mode:{'N','D','U','Z'}) {
                    const Rat rounded_a=rounded(val,p,mode);
                    const Rat minus_three=integer(-3),one=integer(1);
                    emit(p,mode,val);
                    emit(p,mode,add(rounded_a,minus_three));
                    emit(p,mode,sub(rounded_a,minus_three));
                    emit(p,mode,mul(rounded_a,minus_three));
                    emit(p,mode,add(mul(rounded_a,minus_three),one));
                }
            }
        }
        if(count!=4794u)throw std::runtime_error("D79 total count "+std::to_string(count));
    }
};
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        Generator g;g.generate();
        independent_oracle::output(dir,"floatap_rounding.tsv",g.data,write);
    }catch(std::exception const& e){
        std::cerr<<"FloatAP rounding generator: "<<e.what()<<'\n';
        return 1;
    }
}
