// SPDX-License-Identifier: MIT
// Standalone exact R0 random Rational fixture generator. Standard C++20 only.
// No Numeric/GMP/MPFR/FloatMP/Boost, and no fixture input is read.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;

constexpr std::array<unsigned,120> selected_offsets={
    0,3,6,9,0,3,5,8,0,6,6,7,0,9,3,2,1,2,2,5,
    5,9,2,1,1,1,4,7,8,6,3,0,9,4,0,0,3,6,5,4,
    4,5,9,7,1,6,7,4,5,8,4,0,5,5,8,0,2,4,1,4,
    7,4,1,6,3,6,8,5,0,0,6,9,7,4,8,4,7,3,0,6,
    6,5,0,4,0,2,1,8,2,7,4,1,4,0,4,7,6,7,6,0,
    4,9,4,1,9,6,3,6,4,5,8,5,8,6,7,7,1,1,9,0
};
Nat multiply(Nat a,Nat b) {
    Nat result;
    while(!b.zero()) {
        if((b.w.front()&1u)!=0u)result=Nat::add(result,a);
        b=b.shr(1u);
        if(!b.zero())a=a.shl(1u);
    }
    return result;
}
std::size_t bit_length(Nat const& n) {
    return n.zero()?0u:(n.w.size()-1u)*64u+
        64u-static_cast<std::size_t>(std::countl_zero(n.w.back()));
}
std::pair<Nat,Nat> quotient_remainder(Nat remainder,Nat const& divisor) {
    if(divisor.zero()) throw std::invalid_argument("division by zero");
    Nat quotient;
    while(remainder.compare(divisor)>=0) {
        std::size_t shift=bit_length(remainder)-bit_length(divisor);
        Nat aligned=divisor.shl(shift);
        if(aligned.compare(remainder)>0)aligned=divisor.shl(--shift);
        remainder=Nat::sub(remainder,aligned);
        quotient=Nat::add(quotient,Nat::power2(shift));
    }
    return {quotient,remainder};
}
Nat exact_divide(Nat const& a,Nat const& b) {
    auto [q,r]=quotient_remainder(a,b);
    if(!r.zero())throw std::runtime_error("non-exact rational reduction");
    return q;
}
struct Signed {
    Nat mag;
    bool neg=false;
    Signed(Nat n,bool sign=false):mag(std::move(n)),neg(sign&&!mag.zero()){}
    Signed opposite() const {return Signed(mag,!neg);}
    std::string dec() const {return (neg?"-":"")+mag.dec();}
};
Signed sum(Signed const& a,Signed const& b) {
    if(a.neg==b.neg)return Signed(Nat::add(a.mag,b.mag),a.neg);
    const int cmp=a.mag.compare(b.mag);
    if(cmp>=0)return Signed(Nat::sub(a.mag,b.mag),a.neg);
    return Signed(Nat::sub(b.mag,a.mag),b.neg);
}
Signed product(Signed const& a,Signed const& b) {
    return Signed(multiply(a.mag,b.mag),a.neg!=b.neg);
}
struct Fraction {
    Signed num;
    Nat den;
    Fraction(Signed n,Signed d):num(std::move(n)),den(std::move(d.mag)) {
        if(den.zero())throw std::invalid_argument("zero rational denominator");
        num.neg=(num.neg!=d.neg)&&!num.mag.zero();
        reduce();
    }
    Fraction(Signed n,Nat d):num(std::move(n)),den(std::move(d)) {
        if(den.zero())throw std::invalid_argument("zero rational denominator");
        reduce();
    }
    void reduce() {
        const Nat g=Nat::gcd(num.mag,den);
        if(g.compare(Nat(1u))!=0) {
            num.mag=exact_divide(num.mag,g);
            den=exact_divide(den,g);
        }
    }
    std::string str() const {return num.dec()+"/"+den.dec();}
};
Fraction add(Fraction const& a,Fraction const& b,bool subtract=false) {
    const Signed x=product(a.num,Signed(b.den));
    Signed y=product(b.num,Signed(a.den));
    if(subtract)y=y.opposite();
    return Fraction(sum(x,y),multiply(a.den,b.den));
}
Fraction mul(Fraction const& a,Fraction const& b) {
    return Fraction(product(a.num,b.num),multiply(a.den,b.den));
}
Fraction div(Fraction const& a,Fraction const& b) {
    if(b.num.mag.zero())throw std::invalid_argument("fraction divided by zero");
    return Fraction(Signed(multiply(a.num.mag,b.den),a.num.neg!=b.num.neg),
        multiply(a.den,b.num.mag));
}
std::pair<Signed,Signed> floor_ceil(Fraction const& a) {
    auto [q,r]=quotient_remainder(a.num.mag,a.den);
    const Nat rounded=r.zero()?q:Nat::add(q,Nat(1u));
    return a.num.neg?std::pair<Signed,Signed>{Signed(rounded,true),Signed(q,true)}
                    :std::pair<Signed,Signed>{Signed(q),Signed(rounded)};
}
std::uint64_t truncated_binary64(Fraction const& a) {
    const Nat& n=a.num.mag;
    const Nat& d=a.den;
    if(n.zero())return 0u;
    const std::uint64_t sign=a.num.neg?(std::uint64_t{1}<<63u):0u;
    std::int64_t exponent=static_cast<std::int64_t>(bit_length(n))-
        static_cast<std::int64_t>(bit_length(d));
    if((exponent>=0?n.compare(d.shl(static_cast<std::size_t>(exponent)))<0:
         n.shl(static_cast<std::size_t>(-exponent)).compare(d)<0))--exponent;
    if(exponent>1023)return sign|UINT64_C(0x7ff0000000000000);
    if(exponent< -1074)return 0u;
    if(exponent< -1022) {
        auto [mantissa,r]=quotient_remainder(n.shl(1074u),d);
        (void)r;
        if(mantissa.zero())return 0u;
        if(mantissa.w.size()!=1u)throw std::runtime_error("subnormal overflow");
        return sign|mantissa.w.front();
    }
    const std::int64_t shift=52-exponent;
    const Nat top=shift>=0?n.shl(static_cast<std::size_t>(shift)):n;
    const Nat bottom=shift<0?d.shl(static_cast<std::size_t>(-shift)):d;
    auto [mantissa,r]=quotient_remainder(top,bottom);
    (void)r;
    if(mantissa.w.size()!=1u || mantissa.w.front()<(std::uint64_t{1}<<52u)
       ||mantissa.w.front()>=(std::uint64_t{1}<<53u))throw std::runtime_error("bad mantissa");
    return sign|(static_cast<std::uint64_t>(exponent+1023)<<52u)
        |(mantissa.w.front()-(std::uint64_t{1}<<52u));
}
std::string padded_hex(std::uint64_t bits) {
    std::ostringstream stream;
    stream<<std::hex<<std::setfill('0')<<std::setw(16)<<bits;
    return stream.str();
}
int order(Fraction const& a,Fraction const& b) {
    const Signed x=product(a.num,Signed(b.den));
    const Signed y=product(b.num,Signed(a.den));
    if(x.neg!=y.neg)return x.neg?-1:1;
    int cmp=x.mag.compare(y.mag);
    return x.neg?-cmp:cmp;
}
std::string generate() {
    std::string out=
        "# R0 frozen original-source signed Rational reference; no GMP/MPFR runtime arithmetic.\n"
        "# original_iteration limb_count raw_an raw_ad raw_bn raw_bd canon_a canon_b sum diff product square half quotient inverse floor ceil order mpq_get_d_bits_hex\n"
        "# 120 cases selected adaptively, exactly one from each ten original mt19937_64(0x6672616374696f6e) visits.\n"
        "# All 17 limb widths and original sign, zero and order classes preserved; GMP binary64 truncation cross-checked with independent compiled C++ values (visits 0..9).\n";
    std::mt19937_64 rng(0x6672616374696f6eULL);
    for(unsigned iteration=0u;iteration<1200u;++iteration) {
        const std::size_t count=1u+iteration%17u;
        std::array<Signed,4u> values={Signed(Nat()),Signed(Nat()),Signed(Nat()),Signed(Nat())};
        for(std::size_t operand=0u;operand<4u;++operand) {
            std::vector<std::uint64_t> words(count);
            for(auto& word:words)word=rng();
            if(operand%2u!=0u)words.back()|=1u;
            else if(iteration%13u==0u)std::fill(words.begin(),words.end(),0u);
            const Nat value(std::move(words));
            const bool negative=(rng()&1u)!=0u;
            values[operand]=Signed(value,negative);
        }
        if(iteration%10u!=selected_offsets[iteration/10u])continue;
        const Fraction a(values[0],values[1]);
        const Fraction b(values[2],values[3]);
        const Fraction total=add(a,b),difference=add(a,b,true);
        const Fraction result_product=mul(a,b),square=mul(a,a);
        const Fraction half(a.num,a.den.shl(1u));
        const auto [floored,ceiled]=floor_ceil(a);
        out+=row({std::to_string(iteration),std::to_string(count),values[0].dec(),
            values[1].dec(),values[2].dec(),values[3].dec(),a.str(),b.str(),
            total.str(),difference.str(),result_product.str(),square.str(),half.str(),
            b.num.mag.zero()?"-":div(a,b).str(),
            a.num.mag.zero()?"-":Fraction(Signed(a.den,a.num.neg),a.num.mag).str(),
            floored.dec(),ceiled.dec(),std::to_string(order(a,b)),
            padded_hex(truncated_binary64(a))});
    }
    return out;
}
}
int main(int argc,char**argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"rational_random.tsv",generate(),write);
        return 0;
    }catch(std::exception const& error) {
        std::cerr<<"R0 source generator: "<<error.what()<<'\n';
        return 1;
    }
}
