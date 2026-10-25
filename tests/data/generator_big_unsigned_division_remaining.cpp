// SPDX-License-Identifier: MIT
// Reproduce remaining BigUnsigned division oracle fixtures.
// Standalone C++20; no Numeric, GMP, Boost or external runtime.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;

Nat multiply(Nat left,Nat right) {
    Nat result;
    while (!right.zero()) {
        if ((right.w.front()&1u)!=0u) result=Nat::add(result,left);
        right=right.shr(1u);
        if (!right.zero()) left=left.shl(1u);
    }
    return result;
}
std::size_t bit_length(Nat const& x) {
    return x.zero()?0u:(x.w.size()-1u)*64u+
        64u-static_cast<std::size_t>(std::countl_zero(x.w.back()));
}
std::pair<Nat,Nat> divide(Nat remainder,Nat const& divisor) {
    if (divisor.zero()) throw std::invalid_argument("zero divisor");
    Nat quotient;
    while (remainder.compare(divisor)>=0) {
        std::size_t shift=bit_length(remainder)-bit_length(divisor);
        Nat term=divisor.shl(shift);
        if (term.compare(remainder)>0) term=divisor.shl(--shift);
        remainder=Nat::sub(remainder,term);
        quotient=Nat::add(quotient,Nat::power2(shift));
    }
    if (remainder.compare(divisor)>=0 ||
        Nat::add(multiply(quotient,divisor),remainder).compare(
            Nat::add(multiply(quotient,divisor),remainder))!=0)
        throw std::runtime_error("invalid division result");
    return {quotient,remainder};
}

std::array<std::string,6u> a29_general() {
    std::array<std::string,6u> out;
    for (std::size_t part=0u;part<out.size();++part) {
        out[part]="# A29 general division oracle part "+std::to_string(part)+
            ": visits "+std::to_string(part*1000u)+".."+
            std::to_string(part*1000u+999u)+".\n"
            "# repetition n m dividend_hex divisor_hex quotient_hex remainder_hex\n"
            "# Exact mt19937_64(0xB15D1A1D) source stream and source transforms preserved.\n"
            "# No Numeric/GMP/Boost dependency.\n";
    }
    std::mt19937_64 engine(0xB15D1A1DULL);
    const auto max=UINT64_MAX;
    for (unsigned repetition=0u;repetition<6000u;++repetition) {
        const std::size_t n=repetition%31u==0u?96u:
            repetition%7u==0u?32u:static_cast<std::size_t>(engine()%17u);
        const std::size_t m=repetition%37u==0u?1u:
            repetition%11u==0u?96u:1u+static_cast<std::size_t>(engine()%20u);
        std::vector<std::uint64_t> da(n),db(m);
        for (auto& digit:da) digit=engine();
        for (auto& digit:db) digit=engine();
        if (repetition%13u==0u) std::fill(da.begin(),da.end(),max);
        if (repetition%17u==0u) std::fill(db.begin(),db.end(),max);
        if (repetition%23u==0u) std::fill(da.begin(),da.end(),0u);
        if (repetition%29u==0u) {
            std::fill(db.begin(),db.end(),0u);
            db.back()=1u;
        }
        if (repetition%43u==0u) {
            std::fill(db.begin(),db.end(),0u);
            db.back()=std::uint64_t{1}<<63u;
        }
        if (repetition%47u==0u) db.back()=1u;
        if (repetition%53u==0u && n==m) da=db;
        const Nat dividend(std::move(da)),divisor(std::move(db));
        const auto [quotient,remainder]=divide(dividend,divisor);
        if (Nat::add(multiply(quotient,divisor),remainder).compare(dividend)!=0 ||
            remainder.compare(divisor)>=0)
            throw std::runtime_error("invalid A29 Euclidean identity");
        out[repetition/1000u]+=row({
            std::to_string(repetition),std::to_string(n),std::to_string(m),
            dividend.hex(),divisor.hex(),quotient.hex(),remainder.hex()});
    }
    return out;
}
}

int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        const auto a29=a29_general();
        for (std::size_t part=0u;part<a29.size();++part)
            independent_oracle::output(
                dir,"big_unsigned_a29_division_"+std::to_string(part)+".tsv",
                a29[part],write);
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<"BigUnsigned remaining division generator: "
                 <<error.what()<<'\n';
        return 1;
    }
}
