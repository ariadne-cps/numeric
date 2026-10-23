// SPDX-License-Identifier: MIT
// Reproduce existing R3/R4/R8, B32 and R14 exact golden fixtures.
// Standalone C++20; no Numeric, GMP, MPFR, Boost or other oracle dependency.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;

// Binary shift-and-add multiplication, independent of Numeric's multiplier.
Nat multiply(Nat lhs, Nat rhs) {
    Nat result;
    while (!rhs.zero()) {
        if ((rhs.w.front() & 1u) != 0u) result = Nat::add(result,lhs);
        rhs = rhs.shr(1);
        if (!rhs.zero()) lhs = lhs.shl(1);
    }
    return result;
}
Nat power(Nat base, unsigned exponent) {
    Nat result(1);
    while (exponent != 0u) {
        if ((exponent & 1u) != 0u) result=multiply(result,base);
        exponent >>= 1u;
        if (exponent != 0u) base=multiply(base,base);
    }
    return result;
}
std::size_t bit_length(Nat const& x) {
    if (x.zero()) return 0;
    return (x.w.size()-1u)*64u + 64u-static_cast<std::size_t>(std::countl_zero(x.w.back()));
}
std::pair<Nat,Nat> divide(Nat remainder,Nat const& denominator) {
    if (denominator.zero()) throw std::invalid_argument("division by zero");
    Nat quotient;
    while (remainder.compare(denominator)>=0) {
        std::size_t shift=bit_length(remainder)-bit_length(denominator);
        Nat shifted=denominator.shl(shift);
        if (shifted.compare(remainder)>0) shifted=denominator.shl(--shift);
        remainder=Nat::sub(remainder,shifted);
        quotient=Nat::add(quotient,Nat::power2(shift));
    }
    return {quotient,remainder};
}
std::string signed_decimal(Nat const& value,bool negative) {
    return (negative && !value.zero() ? "-" : "")+value.dec();
}
std::string r3() {
    std::string data="# R3 exact signed rational ordering, independent integer cross-product.\n"
                     "# bits expected_sign (-1/0/1).\n";
    for (unsigned bits:{53u,128u,256u,448u,512u,1024u}) {
        const Nat base=Nat::power2(bits);
        const Nat an=Nat::add(base,Nat(13));
        const Nat ad=Nat::add(base.mul_small(3),Nat(19));
        const Nat bn=Nat::add(base,Nat(31));
        const Nat bd=Nat::add(base.mul_small(2),Nat(17));
        const int order=multiply(an,bd).compare(multiply(bn,ad));
        data+=row({std::to_string(bits),std::to_string(order)});
    }
    return data;
}
std::string r4() {
    std::string data="# R4 frozen exact rational floor/ceil expected values, independently computed.\n"
                     "# bits sample_index numerator denominator floor ceil\n"
                     "# All six bit-width boundaries x nine original fraction classes retained.\n";
    for (unsigned bits:{53u,128u,256u,1024u,2048u,4096u}) {
        const Nat base=Nat::power2(bits);
        const Nat below=Nat::sub(base,Nat(1));
        const Nat above=Nat::add(base,Nat(1));
        auto add=[&](unsigned index,Nat const& n,bool neg,Nat const& d) {
            const auto [quotient,remainder]=divide(n,d);
            const Nat rounded=Nat::add(quotient,Nat(remainder.zero()?0u:1u));
            const std::string floor=signed_decimal(neg?rounded:quotient,neg);
            const std::string ceil=signed_decimal(neg?quotient:rounded,neg);
            data+=row({std::to_string(bits),std::to_string(index),signed_decimal(n,neg),
                       d.dec(),floor,ceil});
        };
        add(0,Nat(),false,Nat(1));
        add(1,Nat(1),false,Nat(3));
        add(2,Nat(1),true,Nat(3));
        add(3,below,false,base);
        add(4,below,true,base);
        add(5,Nat(1),false,Nat(1));
        add(6,Nat(1),true,Nat(1));
        add(7,above,false,base);
        add(8,above,true,base);
    }
    return data;
}
std::string r8() {
    std::string data="# R8 original near-2^64 one-limb fractions; exact independent integer cross-multiplication.\n"
                     "# lhs_index rhs_index cmp_result (-1/0/1); 4x4 comparisons retained.\n";
    const std::uint64_t max=UINT64_MAX;
    const std::array<std::pair<Nat,Nat>,4> v={{{Nat(max-17u),Nat(max-7u)},
        {Nat(max-19u),Nat(max-5u)},{Nat(max-13u),Nat(max-3u)},
        {Nat(max-11u),Nat(max-2u)}}};
    for (unsigned i=0;i<v.size();++i)for(unsigned j=0;j<v.size();++j) {
        const int order=multiply(v[i].first,v[j].second).compare(
                        multiply(v[j].first,v[i].second));
        data+=row({std::to_string(i),std::to_string(j),std::to_string(order)});
    }
    return data;
}
std::string integer_power() {
    std::string data="# B32 signed integer powers; independent exact decimal integer arithmetic.\n"
                     "# original_base_decimal exponent result_decimal; zero power includes 0^0 = 1.\n"
                     "# Selects 10 exponent boundaries (0/1/2, sparse and adjacent bits) from original 13.\n";
    const Nat two64=Nat::power2(64);
    const Nat max128=Nat::sub(Nat::power2(128),Nat(1));
    const std::array<std::pair<Nat,bool>,10> values={{{Nat(),false},{Nat(1),false},
        {Nat(1),true},{Nat(2),false},{Nat(2),true},{Nat(17),false},
        {Nat(17),true},{two64,false},{two64,true},{max128,false}}};
    for (auto const& [base,negative]:values) {
        for (unsigned exponent:{0u,1u,2u,3u,5u,7u,8u,9u,11u,12u}) {
            const Nat result=power(base,exponent);
            data+=row({signed_decimal(base,negative),std::to_string(exponent),
                       signed_decimal(result,negative && (exponent%2u)!=0u)});
        }
    }
    return data;
}
std::string r14() {
    std::string data="# R14: original mt19937_64(0x1282561024) inputs, independent exact integer square oracles.\n"
                     "# count mode canonical_input_hex canonical_square_hex\n"
                     "# Each width keeps four distinct patterns: zero, all-ones, top bit, dense random.\n"
                     "# All original 12 RNG visits per width are still consumed at runtime to preserve downstream seeded regression.\n";
    std::mt19937_64 engine(0x1282561024ULL);
    for (std::size_t count:{1u,4u,7u,8u,9u,16u,64u}) {
        for (unsigned mode=0u;mode<12u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=engine();
            if (mode%4u==0u) std::fill(words.begin(),words.end(),0u);
            else if (mode%4u==1u) std::fill(words.begin(),words.end(),UINT64_MAX);
            else if (mode%4u==2u) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
            } else words.back()|=std::uint64_t{1}<<63u;
            if (mode>=4u) continue;
            const Nat input(std::move(words));
            data+=row({std::to_string(count),std::to_string(mode),input.hex(),
                       multiply(input,input).hex()});
        }
    }
    return data;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"rational_r3.tsv",r3(),write);
        independent_oracle::output(dir,"rational_r4.tsv",r4(),write);
        independent_oracle::output(dir,"rational_r8.tsv",r8(),write);
        independent_oracle::output(dir,"integer_power.tsv",integer_power(),write);
        independent_oracle::output(dir,"big_unsigned_r14.tsv",r14(),write);
        return 0;
    } catch(std::exception const& e) {
        std::cerr<<"Exact small-fixture generator: "<<e.what()<<'\n';
        return 1;
    }
}
