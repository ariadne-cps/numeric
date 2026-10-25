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

std::array<std::string,2u> a30_b28() {
    std::array<std::string,2u> out;
    for (std::size_t part=0u;part<out.size();++part) {
        out[part]="# A30 B28 equal-high division oracle part "+
            std::to_string(part)+": visits "+std::to_string(part*1024u)+
            ".."+std::to_string(part*1024u+1023u)+".\n"
            "# repetition count numerator_hex denominator_hex quotient_hex remainder_hex\n"
            "# Exact mt19937_64(0xB280C27) source stream and equal-high transforms preserved.\n"
            "# No Numeric/GMP/Boost dependency.\n";
    }
    std::mt19937_64 engine(0xB280C27ULL);
    for (unsigned repetition=0u;repetition<2048u;++repetition) {
        const std::size_t count=repetition%37u==0u?128u:
            2u+static_cast<std::size_t>(engine()%64u);
        std::vector<std::uint64_t> left(count),right(count);
        for (auto& word:left) word=engine();
        for (auto& word:right) word=engine();
        right.back()|=1u;
        left.back()=right.back();
        if (repetition%4u==0u) {
            left=right;
        } else if (repetition%4u==1u) {
            left=right;
            ++left.front();
        } else if (repetition%4u==2u) {
            right=left;
            ++right.front();
        }
        const Nat numerator(std::move(left)),denominator(std::move(right));
        const auto [quotient,remainder]=divide(numerator,denominator);
        if (Nat::add(multiply(quotient,denominator),remainder)
                .compare(numerator)!=0 ||
            remainder.compare(denominator)>=0)
            throw std::runtime_error("invalid A30 Euclidean identity");
        out[repetition/1024u]+=row({
            std::to_string(repetition),std::to_string(count),
            numerator.hex(),denominator.hex(),quotient.hex(),remainder.hex()});
    }
    return out;
}


std::array<std::string,2u> a31_b36() {
    std::array<std::string,2u> out{
        "# A31 B36 high-ratio division oracle part 0.\n"
        "# ordinal count repetition numerator_hex denominator_hex quotient_hex remainder_hex\n"
        "# Exact mt19937_64(0xB360A11) stream and high-ratio classifications preserved.\n"
        "# No Numeric/GMP/Boost dependency.\n",
        "# A31 B36 high-ratio division oracle part 1.\n"
        "# ordinal count repetition numerator_hex denominator_hex quotient_hex remainder_hex\n"
        "# Exact mt19937_64(0xB360A11) stream and high-ratio classifications preserved.\n"
        "# No Numeric/GMP/Boost dependency.\n"};
    std::mt19937_64 engine(0xB360A11ULL);
    const auto max=UINT64_MAX;
    std::size_t ordinal=0u;
    for (const std::size_t count:
            {2u,3u,4u,7u,8u,16u,32u,64u,96u,128u}) {
        for (unsigned repetition=0u;repetition<128u;++repetition) {
            std::vector<std::uint64_t> left(count),right(count);
            for (auto& limb:left) limb=engine();
            for (auto& limb:right) limb=engine();
            if (repetition%11u==0u) {
                std::fill(left.begin(),left.end(),max);
                std::fill(right.begin(),right.end(),max);
            } else if (repetition%11u==1u) {
                std::fill(left.begin(),left.end(),0u);
                std::fill(right.begin(),right.end(),0u);
            }
            const std::uint64_t high=std::uint64_t{1}<<63u;
            switch (repetition%8u) {
            case 0u: left.back()=high-1u; right.back()=high; break;
            case 1u: left.back()=high+123u; right.back()=high; break;
            case 2u: left.back()=high; right.back()=high; break;
            case 3u: left.back()=high; right.back()=2u; break;
            case 4u: left.back()=max; right.back()=high; break;
            case 5u: left.back()=high; right.back()=high/2u; break;
            case 6u: left.back()=high; right.back()=max; break;
            default: left.back()|=1u; right.back()|=1u; break;
            }
            const Nat numerator(std::move(left)),denominator(std::move(right));
            const auto [quotient,remainder]=divide(numerator,denominator);
            out[ordinal/640u]+=row({
                std::to_string(ordinal),std::to_string(count),
                std::to_string(repetition),numerator.hex(),denominator.hex(),
                quotient.hex(),remainder.hex()});
            ++ordinal;
        }
    }
    if (ordinal!=1280u) throw std::runtime_error("wrong A31 row count");
    return out;
}


std::array<std::string,2u> a32_b29() {
    std::array<std::string,2u> out{
        "# A32 B29 sparse division oracle part 0: visits 0..511.\n"
        "# repetition n m numerator_hex denominator_hex quotient_hex remainder_hex\n"
        "# Exact mt19937_64(0xB290D1F5) conditional source stream preserved.\n"
        "# No Numeric/GMP/Boost dependency.\n",
        "# A32 B29 sparse division oracle part 1: visits 512..1023.\n"
        "# repetition n m numerator_hex denominator_hex quotient_hex remainder_hex\n"
        "# Exact mt19937_64(0xB290D1F5) conditional source stream preserved.\n"
        "# No Numeric/GMP/Boost dependency.\n"};
    std::mt19937_64 engine(0xB290D1F5ULL);
    const auto max=UINT64_MAX;
    for (unsigned repetition=0u;repetition<1024u;++repetition) {
        const std::size_t n=16u+static_cast<std::size_t>(engine()%113u);
        const std::size_t m=8u+static_cast<std::size_t>(engine()%(n-7u));
        std::vector<std::uint64_t> left(n),right(m,0u);
        for (auto& word:left) word=engine();
        right.front()=engine();
        right.back()=std::uint64_t{1}<<(repetition%64u);
        if (repetition%3u==0u) right[m/2u]=engine();
        if (repetition%5u==0u)
            for (std::size_t i=1u;i+1u<m;++i)
                if (i%3u==0u) right[i]=engine();
        if (repetition%7u==0u) std::fill(left.begin(),left.end(),max);
        if (repetition%11u==0u) left.back()=1u;
        const Nat numerator(std::move(left)),denominator(std::move(right));
        const auto [quotient,remainder]=divide(numerator,denominator);
        out[repetition/512u]+=row({
            std::to_string(repetition),std::to_string(n),std::to_string(m),
            numerator.hex(),denominator.hex(),quotient.hex(),remainder.hex()});
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
        const auto a30=a30_b28();
        for (std::size_t part=0u;part<a30.size();++part)
            independent_oracle::output(
                dir,"big_unsigned_a30_b28_"+std::to_string(part)+".tsv",
                a30[part],write);
        const auto a31=a31_b36();
        for (std::size_t part=0u;part<a31.size();++part)
            independent_oracle::output(
                dir,"big_unsigned_a31_b36_"+std::to_string(part)+".tsv",
                a31[part],write);
        const auto a32=a32_b29();
        for (std::size_t part=0u;part<a32.size();++part)
            independent_oracle::output(
                dir,"big_unsigned_a32_b29_"+std::to_string(part)+".tsv",
                a32[part],write);
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<"BigUnsigned remaining division generator: "
                 <<error.what()<<'\n';
        return 1;
    }
}
