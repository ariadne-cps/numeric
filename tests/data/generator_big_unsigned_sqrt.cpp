// SPDX-License-Identifier: MIT
// Independent integer restoring square root, standard C++20 only.
// Deterministic SplitMix64 replaces the unrecoverable original entropy.
#include "oracle_exact_integer.hpp"
#include <bit>
#include <cstdint>
#include <iomanip>
#include <sstream>
namespace {
using independent_oracle::Nat;
// Portable, fully specified SplitMix64 random source: unsigned 64-bit wrap.
// Each draw fills one limb, low limb first; set highest bit to exact width.
struct SplitMix64 {
    std::uint64_t state=0x526f6f74534d3230ULL;
    std::uint64_t next() {
        state+=0x9e3779b97f4a7c15ULL;
        auto x=state;
        x=(x^(x>>30u))*0xbf58476d1ce4e5b9ULL;
        x=(x^(x>>27u))*0x94d049bb133111ebULL;
        return x^(x>>31u);
    }
};
Nat random_exact_width(SplitMix64& prng,unsigned width,bool zero_low_limb=false) {
    if(width==0u)throw std::invalid_argument("zero random width");
    std::vector<std::uint64_t> words((width+63u)/64u);
    for(auto& w:words)w=prng.next();
    const unsigned valid_top=width%64u;
    if(valid_top!=0u)words.back()&=(std::uint64_t{1}<<valid_top)-1u;
    words.back()|=std::uint64_t{1}<<((width-1u)%64u);
    if(zero_low_limb) {
        if(width<=64u)throw std::invalid_argument("not enough limbs for sparse sample");
        words.front()=0u;
    }
    return Nat(std::move(words));
}
std::string hex(Nat const& n) {
    if(n.zero())return "0";
    std::ostringstream o;
    o<<std::hex<<n.w.back();
    for(std::size_t i=n.w.size()-1u;i>0u;--i)
        o<<std::setw(16)<<std::setfill('0')<<n.w[i-1u];
    return o.str();
}
std::size_t bits(Nat const& n) {
    return n.zero()?0u:(n.w.size()-1u)*64u+64u-
        static_cast<std::size_t>(std::countl_zero(n.w.back()));
}
std::pair<Nat,Nat> restoring_sqrtrem(Nat rem) {
    if(rem.zero())return {Nat(),Nat()};
    Nat bit=Nat::power2(2u*((bits(rem)-1u)/2u)),root;
    while(!bit.zero()) {
        const Nat trial=Nat::add(root,bit);
        if(rem.compare(trial)>=0) {
            rem=Nat::sub(rem,trial);
            root=Nat::add(root.shr(1u),bit);
        } else root=root.shr(1u);
        bit=bit.shr(2u);
    }
    return {std::move(root),std::move(rem)};
}
Nat independent_square(Nat const& n) {
    Nat result;
    for(std::size_t word=0;word<n.w.size();++word) {
        std::uint64_t v=n.w[word];
        for(unsigned b=0;b<64u;++b) {
            if(v&1u)result=Nat::add(result,n.shl(word*64u+b));
            v>>=1u;
            if(v==0u)break;
        }
    }
    return result;
}
std::string generate() {
    std::string out=
        "# BigUnsigned integer sqrtrem: canonical input_hex expected_floor_root_hex expected_remainder_hex\n"
        "# Source: SplitMix64 seed 0x526f6f74534d3230, modulo-2^64 words, low limb first.\n"
        "# Source widths: 53/64/65 (3 each), 127/128/256/512/1024/2048/4096 (4 each).\n"
        "# Each four-sample group includes a forced-zero low 64-bit limb at index 1.\n"
        "# Eleven more generated roots give square-1, square, square+1 cases.\n"
        "# Oracle: independent radix-4 restoring sqrtrem and shift/add square bounds.\n"
        "# No GMP/independent exact arithmetic or BigUnsigned output used to compute expected values.\n";
    std::size_t count=0u;
    const auto append=[&](Nat const& n) {
        const auto [root,rem]=restoring_sqrtrem(n);
        const Nat sq=independent_square(root);
        if(sq.compare(n)>0||Nat::sub(n,sq).compare(rem)!=0||
           independent_square(Nat::add(root,Nat(1u))).compare(n)<=0)
            throw std::runtime_error("independent square root validation failed");
        out+=hex(n)+" "+hex(root)+" "+hex(rem)+"\n";
        ++count;
    };
    for(std::uint64_t n:{0u,1u,2u,3u,4u,7u,8u,9u,10u,15u,16u,
                         17u,24u,25u,26u,63u,64u,65u,255u,256u,257u})
        append(Nat(n));
    append(Nat(5u));
    for(unsigned p:{7u,31u,32u,53u,63u,64u,65u,66u,127u,128u,129u,
                    255u,256u,257u,511u,512u,513u,1023u,1024u,1025u,
                    2047u,2048u,4096u}) {
        const Nat power=Nat::power2(p);
        append(Nat::sub(power,Nat(1u)));
        append(power);
        append(Nat::add(power,Nat(1u)));
    }
    SplitMix64 prng;
    for(unsigned width:{53u,64u,65u})
        for(unsigned i=0u;i<3u;++i)append(random_exact_width(prng,width));
    for(unsigned width:{127u,128u,256u,512u,1024u,2048u,4096u})
        for(unsigned i=0u;i<4u;++i)
            append(random_exact_width(prng,width,i==1u));
    for(unsigned width:{31u,32u,63u,64u,65u,127u,128u,256u,512u,1024u,2048u}) {
        const auto root=random_exact_width(prng,width);
        const auto sq=independent_square(root);
        append(Nat::sub(sq,Nat(1u)));append(sq);append(Nat::add(sq,Nat(1u)));
    }
    if(count!=161u)throw std::runtime_error("wrong sqrt row count");
    return out;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"big_unsigned_sqrt.tsv",generate(),write);
        return 0;
    } catch(std::exception const& e) {
        std::cerr<<"BigUnsigned sqrt generator: "<<e.what()<<'\n';
        return 1;
    }
}
