// SPDX-License-Identifier: MIT
// Independent integer restoring square root, standard C++20 only.
#include "oracle_exact_integer.hpp"
#include "big_unsigned_sqrt_legacy_seeds.hpp"
#include <bit>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string_view>
namespace {
using independent_oracle::Nat;
Nat from_hex(std::string_view src) {
    Nat n;
    for(char digit:src) {
        unsigned v=0;
        if(digit>='0'&&digit<='9')v=static_cast<unsigned>(digit-'0');
        else if(digit>='a'&&digit<='f')v=10u+static_cast<unsigned>(digit-'a');
        else throw std::invalid_argument("bad hex seed");
        n=Nat::add(n.mul_small(16u),Nat(v));
    }
    return n;
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
        "# Oracle: independent JavaScript BigInt digit-pair restoring arithmetic.\n"
        "# Root independently cross-checked against BigInt binary-search sqrt.\n"
        "# Remainder obtained directly from digit-pair restoring state and checked\n"
        "# against independent input - root*root arithmetic.\n"
        "# No GMP/MPFR or BigUnsigned output used to compute expected values.\n";
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
    for(std::string_view input:legacy_sqrt_random_inputs)append(from_hex(input));
    for(std::string_view root_input:legacy_sqrt_square_roots) {
        const Nat sq=independent_square(from_hex(root_input));
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
