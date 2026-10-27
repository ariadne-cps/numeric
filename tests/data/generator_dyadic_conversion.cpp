// SPDX-License-Identifier: MIT
// Independent truncating dyadic -> binary64, assembled as IEEE integer bits.
#include "oracle_exact_integer.hpp"
namespace {
using independent_oracle::Nat;
std::uint64_t binary64(Nat const& n,unsigned denominator_exponent,bool negative) {
    if (n.zero()) return 0u;
    const auto width=static_cast<unsigned>((n.w.size()-1u)*64u+64u-
                                          static_cast<unsigned>(std::countl_zero(n.w.back())));
    const int exponent=static_cast<int>(width)-1-static_cast<int>(denominator_exponent);
    const std::uint64_t sign=negative?std::uint64_t{1}<<63u:0u;
    if (exponent>1023) return sign|UINT64_C(0x7ff0000000000000);
    if (exponent< -1074) return 0u; // Dyadic's documented canonical underflow zero.
    const int grid=std::max(-1074,exponent-52);
    const int shift=static_cast<int>(denominator_exponent)+grid;
    const auto quotient=shift>=0?n.shr(static_cast<unsigned>(shift)):n.shl(static_cast<unsigned>(-shift));
    if (quotient.zero()) return 0u;
    if (quotient.w.size()!=1u) throw std::logic_error("binary64 significand overflow");
    if (exponent< -1022) return sign|quotient.w.front();
    return sign|(static_cast<std::uint64_t>(exponent+1023)<<52u)|
        (quotient.w.front()-(std::uint64_t{1}<<52u));
}
std::string fixture() {
    const std::array<Nat,11> values={Nat(),Nat(1),Nat(3),
        Nat::sub(Nat::power2(53),Nat(1)),Nat::add(Nat::power2(53),Nat(1)),
        Nat::add(Nat::power2(64),Nat(3)),Nat::add(Nat::power2(127),Nat(1)),Nat::power2(116),
        Nat::power2(1023),Nat::power2(1024),Nat::add(Nat::power2(1024),Nat(1))};
    std::string result="# Exact truncating dyadic binary64 oracle, integer bit assembly.\n"
        "# id coefficient denominator_exponent binary64_bits exact_binary64 denominator floor ceil\n";
    unsigned id=0u;
    for (auto const& n:values) for (unsigned exponent:{0u,1u,53u,64u,1022u,1074u,1075u,1200u}) {
        for (bool negative:{false,true}) {
            if (negative && n.zero()) continue;
            const auto bits=binary64(n,exponent,negative);
            const auto width=n.zero()?0u:(n.w.size()-1u)*64u+64u-
                static_cast<unsigned>(std::countl_zero(n.w.back()));
            const auto trailing=n.zero()?0u:n.trailing_zero_bits();
            const auto grid=std::max(-1074,static_cast<int>(width)-1-static_cast<int>(exponent)-52);
            const int discarded=static_cast<int>(exponent)+grid;
            const bool exact=n.zero() || ((bits&UINT64_C(0x7ff0000000000000))!=UINT64_C(0x7ff0000000000000)
                && (discarded<=0 || static_cast<std::size_t>(discarded)<=trailing));
            const auto integer=n.shr(exponent);
            const bool fraction=integer.shl(exponent).compare(n)!=0;
            const auto away=Nat::add(integer,Nat(fraction?1u:0u));
            const auto signed_integer=[&](Nat const& value) {
                return (negative && !value.zero()?"-":"")+value.dec();
            };
            result+=independent_oracle::row({std::to_string(id++),(negative?"-":"")+n.dec(),
                std::to_string(exponent),std::to_string(bits),exact?"1":"0",
                Nat::power2(exponent).dec(),signed_integer(negative?away:integer),
                signed_integer(negative?integer:away)});
        }
    }
    return result;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"dyadic_conversion.tsv",fixture(),write);
    } catch (std::exception const& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
