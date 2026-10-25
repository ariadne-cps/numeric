// SPDX-License-Identifier: MIT
// Exact two-stage (FloatAP p-bit -> IEEE binary64) rounding oracle, C++20.
// No GMP, independent exact arithmetic, FloatAP, Numeric, Boost, or host floating-point conversions.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;

constexpr std::uint64_t negative_bit=std::uint64_t{1}<<63u;
constexpr std::uint64_t infinity_bits=UINT64_C(0x7ff0000000000000);
constexpr std::uint64_t max_finite=UINT64_C(0x7fefffffffffffff);

std::size_t bit_length(Nat const& value) {
    return value.zero()?0u:(value.w.size()-1u)*64u+
        64u-static_cast<std::size_t>(std::countl_zero(value.w.back()));
}
std::pair<Nat,Nat> divrem(Nat remainder,Nat const& denominator) {
    if(denominator.zero())throw std::invalid_argument("zero divisor");
    Nat result;
    while(remainder.compare(denominator)>=0) {
        std::size_t shift=bit_length(remainder)-bit_length(denominator);
        Nat aligned=denominator.shl(shift);
        if(aligned.compare(remainder)>0)aligned=denominator.shl(--shift);
        remainder=Nat::sub(remainder,aligned);
        result=Nat::add(result,Nat::power2(shift));
    }
    return {result,remainder};
}
std::int64_t exponent(Nat const& numerator,Nat const& denominator) {
    if(numerator.zero()||denominator.zero())throw std::invalid_argument("zero exponent source");
    std::int64_t exponent=static_cast<std::int64_t>(bit_length(numerator))-
        static_cast<std::int64_t>(bit_length(denominator));
    if(exponent>=0) {
        if(numerator.compare(denominator.shl(static_cast<std::size_t>(exponent)))<0)--exponent;
    } else {
        if(numerator.shl(static_cast<std::size_t>(-exponent)).compare(denominator)<0)--exponent;
    }
    return exponent;
}
struct Rational {
    Nat magnitude;
    Nat denominator;
    bool negative=false;
    Rational(Nat numerator,Nat den,bool sign=false):magnitude(std::move(numerator)),denominator(std::move(den)),negative(sign&&!magnitude.zero()) {
        if(denominator.zero())throw std::invalid_argument("zero rational denominator");
    }
    std::string num() const {return (negative?"-":"")+magnitude.dec();}
    std::string den() const {return denominator.dec();}
};
Nat quantize(Nat const& numerator,Nat const& denominator,std::int64_t grid_exponent,char mode,bool negative) {
    Nat scaled_n,scaled_d;
    scaled_n.w.assign(numerator.w.begin(),numerator.w.end());
    scaled_d.w.assign(denominator.w.begin(),denominator.w.end());
    if(grid_exponent>=0)scaled_d=denominator.shl(static_cast<std::size_t>(grid_exponent));
    else scaled_n=numerator.shl(static_cast<std::size_t>(-grid_exponent));
    auto [q,r]=divrem(scaled_n,scaled_d);
    if(r.zero())return q;
    bool increment=false;
    if(mode=='N') {
        const int relation=r.shl(1u).compare(scaled_d);
        increment=relation>0 || (relation==0 && !q.zero() && (q.w.front()&1u)!=0u);
    } else if(mode=='D')increment=negative;
    else if(mode=='U')increment=!negative;
    else if(mode!='Z')throw std::invalid_argument("unknown rounding mode");
    return increment?Nat::add(q,Nat(1u)):q;
}
Rational rounded_to_precision(Rational const& source,unsigned precision) {
    if(source.magnitude.zero())throw std::invalid_argument("input list contains no zero");
    const std::int64_t grid=exponent(source.magnitude,source.denominator)-
        static_cast<std::int64_t>(precision)+1;
    Nat significand=quantize(source.magnitude,source.denominator,grid,'N',source.negative);
    if(grid>=0)return Rational(significand.shl(static_cast<std::size_t>(grid)),Nat(1u),source.negative);
    return Rational(std::move(significand),Nat::power2(static_cast<std::size_t>(-grid)),source.negative);
}
std::uint64_t binary64(Rational const& value,char mode) {
    const bool negative=value.negative;
    const std::uint64_t sign=negative?negative_bit:0u;
    if(value.magnitude.zero())return sign;
    std::int64_t power=exponent(value.magnitude,value.denominator);
    if(power>1023) {
        const bool away=mode=='N'||(mode=='D'&&negative)||(mode=='U'&&!negative);
        return sign|(away?infinity_bits:max_finite);
    }
    const std::int64_t grid=std::max(std::int64_t{-1074},power-52);
    Nat q=quantize(value.magnitude,value.denominator,grid,mode,negative);
    if(q.zero())return sign;
    if(grid==-1074) {
        if(q.compare(Nat::power2(52u))>=0) return sign|(std::uint64_t{1}<<52u);
        if(q.w.size()!=1u)throw std::runtime_error("subnormal mantissa width");
        return sign|q.w.front();
    }
    if(q.compare(Nat::power2(53u))==0){q=q.shr(1u);++power;}
    if(power>1023)return sign|infinity_bits;
    if(q.w.size()!=1u || q.w.front()<(std::uint64_t{1}<<52u))
        throw std::runtime_error("normal mantissa width");
    return sign|(static_cast<std::uint64_t>(power+1023)<<52u)
        |(q.w.front()-(std::uint64_t{1}<<52u));
}
std::vector<Rational> original_inputs() {
    const Nat one(1u);
    const Nat two53=Nat::power2(53u);
    return {
        Rational(one,one),Rational(one,one,true),Rational(one,Nat(3u)),
        Rational(Nat(17u),Nat(13u),true),
        Rational(Nat::add(two53,one),one),Rational(Nat::add(two53,Nat(3u)),one),
        Rational(one,Nat::power2(1200u)),Rational(one,Nat::power2(1200u),true),
        Rational(one,Nat::power2(1076u)),Rational(one,Nat::power2(1075u)),
        Rational(one,Nat::power2(1075u),true),Rational(one,Nat::power2(1074u)),
        Rational(one,Nat::power2(1022u)),Rational(Nat::power2(1023u),one),
        Rational(Nat::power2(1024u),one),Rational(Nat::power2(1024u),one,true),
        Rational(Nat::sub(Nat::power2(54u),one).shl(970u),one)
    };
}
std::string generate() {
    std::string out=
        "# D98: independent IEEE 754 binary64 bit-pattern oracle for FloatAP::get_d.\n"
        "# bits mode original_num original_den expected_unsigned_64_bits 1\n"
        "# All expected outputs are EXACT unsigned integer bit-patterns, NOT real-valued doubles.\n"
        "# A source Rational is independently rounded to p-bit binary nearest/even BEFORE\n"
        "# IEEE 754 binary64 conversion, including p=53 input ties and inexact rationals.\n"
        "# Normal 53-bit mantissa, subnormal 2^-1074 grid, zero and overflow handled by\n"
        "# exact integer quotient/remainder and midpoint-even tests in N/D/U/Z modes.\n"
        "# Sign bit is retained for zeros; infinities and maximum finite values encoded\n"
        "# as unsigned 64-bit integers, avoiding any host floating point dependencies.\n"
        "# No independent exact arithmetic, GMP, independent exact arithmetic, FloatAP or host double used as a reference.\n";
    for(unsigned precision:{53u,128u,256u,1024u}) {
        for(const auto& source:original_inputs()) {
            const Rational rounded=rounded_to_precision(source,precision);
            for(char mode:{'N','D','U','Z'}) {
                out+=row({std::to_string(precision),std::string(1u,mode),
                    source.num(),source.den(),std::to_string(binary64(rounded,mode)),"1"});
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
        independent_oracle::output(dir,"floatap_get_d.tsv",generate(),write);
        return 0;
    } catch(std::exception const& e) {
        std::cerr<<"FloatAP get_d generator: "<<e.what()<<'\n';
        return 1;
    }
}
