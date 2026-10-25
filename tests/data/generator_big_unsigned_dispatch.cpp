// SPDX-License-Identifier: MIT
// Reproduce BigUnsigned B39/B40/B54 exact fixture source streams.
// Independent standard C++20; no Numeric, GMP, independent exact arithmetic, Boost or independent exact arithmetic.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;
using Words=std::vector<std::uint64_t>;
using Pair=std::pair<std::size_t,std::size_t>;

Nat multiply(Nat lhs,Nat rhs) {
    Nat answer;
    while (!rhs.zero()) {
        if ((rhs.w.front() & 1u)!=0u) answer=Nat::add(answer,lhs);
        rhs=rhs.shr(1u);
        if(!rhs.zero()) lhs=lhs.shl(1u);
    }
    return answer;
}
Nat truncated(Nat product,std::size_t count) {
    if (product.w.size()>count) product.w.resize(count);
    product.trim();
    return product;
}
constexpr std::array<Pair,19> b39_widths={{{7u,8u},{8u,8u},{15u,16u},{16u,16u},
    {32u,32u},{63u,63u},{63u,64u},{64u,63u},{64u,64u},{64u,65u},
    {65u,64u},{64u,85u},{64u,86u},{85u,64u},{86u,64u},{96u,64u},
    {128u,64u},{64u,128u},{128u,128u}}};
constexpr std::array<Pair,18> b40_widths={{{7u,64u},{8u,64u},{16u,16u},{31u,32u},
    {32u,31u},{32u,32u},{32u,33u},{33u,32u},{33u,33u},{63u,63u},
    {63u,64u},{64u,63u},{64u,64u},{64u,65u},{65u,64u},
    {64u,128u},{128u,64u},{128u,128u}}};
constexpr std::array<Pair,17> b54_widths={{{32u,32u},{33u,33u},{63u,63u},
    {63u,64u},{64u,63u},{64u,64u},{64u,65u},{65u,64u},{65u,65u},
    {64u,85u},{64u,86u},{85u,64u},{86u,64u},{96u,96u},
    {127u,128u},{128u,127u},{128u,128u}}};

std::size_t b39_zero_rows(std::size_t count,unsigned mode) {
    switch(mode%7u) {
        case 0u:return 0u;
        case 1u:return 1u;
        case 2u:return count/8u;
        case 3u:return count/8u+1u;
        case 4u:return count/2u;
        case 5u:return count-2u;
        default:return count/4u;
    }
}
std::string make_b39() {
    std::string out=
        "# B39 independent product/square/truncated product reference from original mt19937_64(0xB39DA7A).\n"
        "# n m repetition lhs_hex rhs_hex product_hex lhs_square_hex truncated_hex\n"
        "# 19 original operand width pairs. Repetitions 0..6 and 9 cover all 7 left and right zero-row density modes plus both low-product lengths.\n"
        "# Original 20 RNG visits per width pair consumed offline; references from exact independent integers.\n";
    std::mt19937_64 rng(0xB39DA7AULL);
    for(const auto& [n,m]:b39_widths) {
        for(unsigned rep=0u;rep<20u;++rep) {
            Words left(n),right(m);
            for(auto& x:left)x=rng();
            for(auto& x:right)x=rng();
            const auto left_zeros=b39_zero_rows(n,rep);
            const unsigned rhs_mode=rep%3u==0u ? rep : rep+3u;
            const auto right_zeros=b39_zero_rows(m,rhs_mode);
            for(std::size_t i=0u;i<left_zeros && i+1u<n;++i) left[i]=0u;
            for(std::size_t i=0u;i<right_zeros && i+1u<m;++i) right[i]=0u;
            left.back()|=std::uint64_t{1}<<63u;
            right.back()|=std::uint64_t{1}<<63u;
            if(rep>6u && rep!=9u) continue;
            const Nat lhs(std::move(left)),rhs(std::move(right));
            const Nat result=multiply(lhs,rhs);
            const auto low=(rep%2u)==0u ? (n+m)/2u : n+m-1u;
            out+=row({std::to_string(n),std::to_string(m),std::to_string(rep),
                lhs.hex(),rhs.hex(),result.hex(),multiply(lhs,lhs).hex(),
                truncated(result,low).hex()});
        }
    }
    return out;
}
std::size_t b40_zero_rows(std::size_t count,unsigned mode) {
    switch(mode%8u) {
        case 0u:return 0u;
        case 1u:return 1u;
        case 2u:return count/8u;
        case 3u:return count/8u+1u;
        case 4u:return count/2u;
        case 5u:return count-1u;
        case 6u:return 2u;
        default:return count/4u;
    }
}
std::string make_b40() {
    std::string out=
        "# B40 independent exact unsigned multiplication/short-product/square oracle.\n"
        "# n m original_rep original_lhs_hex original_rhs_hex product_hex truncated_hex lhs_square_hex\n"
        "# All 18 source width pairs and 8 zero-row density classes plus original saturated rep 15.\n"
        "# Original C++ mt19937_64(0xB4004040) full source stream; expected results calculated offline with exact integers.\n";
    std::mt19937_64 rng(0xB4004040ULL);
    for(const auto& [n,m]:b40_widths) {
        for(unsigned rep=0u;rep<16u;++rep) {
            Words left(n),right(m);
            for(auto& x:left)x=rng();
            for(auto& x:right)x=rng();
            const auto zero_a=b40_zero_rows(n,rep),zero_b=b40_zero_rows(m,rep+3u);
            for(std::size_t i=0u;i<zero_a && i+1u<n;++i)left[i]=0u;
            for(std::size_t i=0u;i<zero_b && i+1u<m;++i)right[i]=0u;
            if(rep==15u) {
                std::fill(left.begin(),left.end(),UINT64_MAX);
                std::fill(right.begin(),right.end(),UINT64_MAX);
            }
            left.back()|=std::uint64_t{1}<<63u;
            right.back()|=std::uint64_t{1}<<63u;
            if(rep>7u && rep!=15u)continue;
            const Nat lhs(std::move(left)),rhs(std::move(right));
            const Nat product=multiply(lhs,rhs);
            const auto low=rep%2u==0u ? 1u : n+m-1u;
            out+=row({std::to_string(n),std::to_string(m),std::to_string(rep),
                lhs.hex(),rhs.hex(),product.hex(),truncated(product,low).hex(),
                multiply(lhs,lhs).hex()});
        }
    }
    return out;
}
void b54_pattern(Words& words,unsigned mode) {
    const std::size_t count=words.size();
    switch(mode%8u) {
        case 0u:break;
        case 1u:std::fill(words.begin(),words.end(),UINT64_MAX);break;
        case 2u:std::fill(words.begin(),words.end(),0u);words.front()=13u;break;
        case 3u:std::fill_n(words.begin(),count/8u,0u);break;
        case 4u:std::fill_n(words.begin(),count/8u+1u,0u);break;
        case 5u:for(std::size_t i=0u;i<count;i+=2u)words[i]=0u;break;
        case 6u:for(std::size_t i=1u;i+1u<count;++i)words[i]=0u;break;
        default:for(std::size_t i=0u;i<count;i+=8u)words[i]=0u;break;
    }
    words.back()|=std::uint64_t{1}<<63u;
}
std::mt19937_64 init_b54_product_stream() {
    return std::mt19937_64(0xB54A6C4ULL);
}
std::string make_b54_product(std::mt19937_64& rng) {
    std::string out=
        "# B54 frozen exact full-product oracle from mt19937_64(0xB54A6C4).\n"
        "# n m mode lhs_hex rhs_hex exact_product_hex\n"
        "# 17 width pairs, all 8 distinct pattern%8 modes; original RNG retained for downstream square suite.\n";
    for(const auto& [n,m]:b54_widths) {
        for(unsigned mode=0u;mode<32u;++mode) {
            Words left(n),right(m);
            for(auto& x:left)x=rng();
            for(auto& x:right)x=rng();
            b54_pattern(left,mode);
            b54_pattern(right,mode+3u);
            if(mode>=8u)continue;
            const Nat lhs(std::move(left)),rhs(std::move(right));
            out+=row({std::to_string(n),std::to_string(m),std::to_string(mode),
                lhs.hex(),rhs.hex(),multiply(lhs,rhs).hex()});
        }
    }
    return out;
}
std::string make_b54_square(std::mt19937_64& rng) {
    std::string out=
        "# B54 signed-independent unsigned square oracle from exact integer arithmetic.\n"
        "# original_count original_mode source_hex expected_square_hex\n"
        "# All 13 widths and eight mode%8 classes, 832 original checks compressed to 104.\n"
        "# Original mt19937_64(0xB54A6C4) advanced through 17*32 product corpus before these 64 visits per width.\n"
        "# Only original-source selected values are frozen; expected square never computed using BigUnsigned.\n";
    for(std::size_t count:{7u,8u,9u,16u,31u,32u,33u,63u,64u,65u,66u,96u,128u}) {
        for(unsigned mode=0u;mode<64u;++mode) {
            Words words(count);
            for(auto& x:words)x=rng();
            if(mode%8u==0u)std::fill(words.begin(),words.end(),UINT64_MAX);
            else if(mode%8u==1u){std::fill(words.begin(),words.end(),0u);words.front()=13u;}
            else if(mode%8u==2u)std::fill_n(words.begin(),count/8u,0u);
            else if(mode%8u==3u)std::fill_n(words.begin(),count/8u+1u,0u);
            else if(mode%8u==4u){for(std::size_t i=0u;i<count;i+=2u)words[i]=0u;}
            else if(mode%8u==5u){for(std::size_t i=1u;i+1u<count;++i)words[i]=0u;}
            else if(mode%8u==6u){for(std::size_t i=0u;i<count;i+=8u)words[i]=0u;}
            words.back()|=std::uint64_t{1}<<63u;
            if(mode>=8u)continue;
            const Nat lhs(std::move(words));
            out+=row({std::to_string(count),std::to_string(mode),lhs.hex(),multiply(lhs,lhs).hex()});
        }
    }
    return out;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"big_unsigned_b39.tsv",make_b39(),write);
        independent_oracle::output(dir,"big_unsigned_b40.tsv",make_b40(),write);
        auto rng=init_b54_product_stream();
        independent_oracle::output(dir,"big_unsigned_b54_product.tsv",make_b54_product(rng),write);
        independent_oracle::output(dir,"big_unsigned_b54_square.tsv",make_b54_square(rng),write);
        return 0;
    } catch(std::exception const& e) {
        std::cerr<<"B39/B40/B54 generator: "<<e.what()<<'\n';
        return 1;
    }
}
