// SPDX-License-Identifier: MIT
// Independent exponent-range arithmetic: arbitrary-size signed exponents,
// bounded unsigned significand rounding. No Numeric/GMP/MPFR/Boost.
#include "oracle_exact_integer.hpp"
namespace {
using independent_oracle::Nat;
using I=std::int64_t;
struct Signed {
    Nat magnitude;
    bool negative;
};
Signed number(I x) { return {Nat(x<0?static_cast<std::uint64_t>(-(x+1))+1u:static_cast<std::uint64_t>(x)),x<0}; }
Signed negate(Signed x) { if (!x.magnitude.zero()) x.negative=!x.negative; return x; }
Signed add(Signed a,Signed b) {
    if (a.negative==b.negative) return {Nat::add(a.magnitude,b.magnitude),a.negative};
    const int order=a.magnitude.compare(b.magnitude);
    Nat n=order>=0?Nat::sub(a.magnitude,b.magnitude):Nat::sub(b.magnitude,a.magnitude);
    return {n,!n.zero()&&(order>=0?a.negative:b.negative)};
}
bool fits(Signed const& x) {
    const int comparison=x.magnitude.compare(Nat::power2(63u));
    return x.negative?comparison<=0:comparison<0;
}
std::string text(Signed const& x) { return (x.negative?"-":"")+x.magnitude.dec(); }
I exponent(std::string const& token) {
    return token=="min"?std::numeric_limits<I>::min():std::numeric_limits<I>::max();
}
std::uint64_t significand(std::string const& token,unsigned p) {
    return token=="min"?std::uint64_t{1}<<(p-1u):(std::uint64_t{1}<<p)-1u;
}
std::uint64_t rounded(std::uint64_t n,std::uint64_t d,char mode) {
    auto q=n/d; const auto r=n%d;
    if (mode=='U' && r!=0u) ++q;
    else if (mode=='N' && (2u*r>d || (2u*r==d && q%2u!=0u))) ++q;
    return q;
}
std::string fixture() {
    std::string out="# Independent native FloatMP exponent-range oracle; no huge powers of two allocated.\n"
        "# id operation mode precision lhs rhs_or_shift status expected_significand expected_exponent\n";
    unsigned id=0u;
    for (unsigned p:{2u,7u}) for (std::string const op:{"mul","div","shift"})
    for (std::string const a:{"min","max"}) for (std::string const b:{"min","max","-1","0","1"}) {
        if ((op=="shift")!=(b!="min"&&b!="max")) continue;
        for (char mode:{'N','D','U','Z'}) {
            if (op=="shift" && mode!='N') continue;
            Signed scale=number(exponent(a));
            auto coefficient=significand(a,p);
            if (op=="shift") scale=add(scale,number(std::stoll(b)));
            else if (op=="mul") {
                const auto product=coefficient*significand(b,p);
                const auto drop=static_cast<unsigned>(std::bit_width(product))-p;
                coefficient=rounded(product,std::uint64_t{1}<<drop,mode);
                scale=add(add(scale,number(exponent(b))),number(static_cast<I>(drop)));
            } else {
                const auto denominator=significand(b,p);
                const int ratio_exponent=coefficient<denominator?-1:0;
                const unsigned shift=p-1u+static_cast<unsigned>(-ratio_exponent);
                coefficient=rounded(coefficient<<shift,denominator,mode);
                scale=add(add(scale,negate(number(exponent(b)))),number(-static_cast<I>(shift)));
            }
            if (coefficient==(std::uint64_t{1}<<p)) {
                coefficient>>=1u;scale=add(scale,number(1));
            }
            const bool valid=fits(scale);
            out+=independent_oracle::row({std::to_string(id++),op,std::string(1,mode),
                std::to_string(p),a,b,valid?"finite":"overflow",valid?std::to_string(coefficient):"0",
                valid?text(scale):"0"});
        }
    }
    return out;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"floatmp_exponent_limits.tsv",fixture(),write);
    } catch (std::exception const& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
