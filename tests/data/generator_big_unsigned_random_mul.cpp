// SPDX-License-Identifier: MIT
// Original-source exact BigUnsigned randomized multiplication fixture.
// Standalone C++20; no Numeric, GMP, MPFR, Boost or FloatMP.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;
constexpr std::string_view offsets =
    "042396695972444984963916999444443075045497348994999952969944989839924949497910619699994444444494049444999646099962494499949994941996925149469494148949699494299444459481644549474449499445944928997430354545494444746424548449437395949484349895849924409899749495374042044849994431493444929443941662191994384549530594496749704894499699491643944999474441994409949449494444990494934139441904644444940095749944394844149996899201924599339990740805960945456099490430445915438909673159747479917441254949849990961960744029444284689419171408925344449439429944789596062569944449944948448464492944944744044896897459";
static_assert(offsets.size()==600u);

Nat multiply(Nat a,Nat b) {
    Nat answer;
    while(!b.zero()) {
        if((b.w.front()&1u)!=0u)answer=Nat::add(answer,a);
        b=b.shr(1u);
        if(!b.zero())a=a.shl(1u);
    }
    return answer;
}
Nat low(Nat p,std::size_t count) {
    if(p.w.size()>count)p.w.resize(count);
    p.trim();
    return p;
}
std::string generate() {
    std::string result=
        "# Frozen original-source BigUnsigned randomized multiplication; independent arbitrary-precision exact integer oracle.\n"
        "# rep source_n source_m low_limbs source_x_hex source_y_hex full_product_hex square_x_hex low_product_hex\n"
        "# std::mt19937_64(0x1282561024) after R14 draws, first 10 original source metadata independently compared with compiled C++.\n"
        "# One repetition per consecutive block of ten; coverage-aware greedy choice, not naive fixed modulus.\n"
        "# The original 34 width values are all represented for both operands; all low-limit classes survive.\n";
    std::mt19937_64 rng(0x1282561024ULL);
    // R14 precedes the random multiplication loop and uses the same RNG.
    for(std::size_t count:{1u,4u,7u,8u,9u,16u,64u})
        for(unsigned mode=0u;mode<12u;++mode)
            for(std::size_t j=0u;j<count;++j)(void)rng();
    for(unsigned rep=0u;rep<6000u;++rep) {
        const std::size_t n=rng()%34u;
        const std::size_t m=rng()%34u;
        std::vector<std::uint64_t> left(n),right(m);
        for(auto& v:left)v=rng();
        for(auto& v:right)v=rng();
        if(rep%7u==0u)std::fill(left.begin(),left.end(),UINT64_MAX);
        if(rep%11u==0u)std::fill(right.begin(),right.end(),UINT64_MAX);
        if(rep%13u==0u)std::fill(left.begin(),left.end(),0u);
        if(rep%17u==0u)std::fill(right.begin(),right.end(),0u);
        const std::size_t trunc=rep%5u==0u?0u:rep%5u==1u?1u:
            rep%5u==2u?4u:rep%5u==3u?5u:rng()%80u;
        if(rep%10u != static_cast<unsigned>(offsets[rep/10u]-'0'))continue;
        const Nat a(std::move(left)),b(std::move(right));
        const Nat p=multiply(a,b);
        result+=row({std::to_string(rep),std::to_string(n),std::to_string(m),
                     std::to_string(trunc),a.hex(),b.hex(),p.hex(),
                     multiply(a,a).hex(),low(p,trunc).hex()});
    }
    return result;
}
}
int main(int argc,char**argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"big_unsigned_random_mul.tsv",generate(),write);
        return 0;
    } catch(std::exception const& e) {
        std::cerr<<"Random multiplication generator: "<<e.what()<<'\n';
        return 1;
    }
}
