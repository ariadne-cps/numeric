// SPDX-License-Identifier: MIT
// Source-reproducible Integer and BigUnsigned structural exact oracles.
// Standalone C++20. No Numeric, GMP, MPFR, FloatMP or Boost.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;
using Words = std::vector<std::uint64_t>;

Nat multiply(Nat lhs, Nat rhs) {
    Nat result;
    while (!rhs.zero()) {
        if ((rhs.w.front() & 1u)!=0u) result=Nat::add(result,lhs);
        rhs=rhs.shr(1u);
        if (!rhs.zero()) lhs=lhs.shl(1u);
    }
    return result;
}
Nat power(Nat base,unsigned exponent) {
    Nat result(1u);
    while (exponent!=0u) {
        if ((exponent & 1u)!=0u) result=multiply(result,base);
        exponent>>=1u;
        if (exponent!=0u) base=multiply(base,base);
    }
    return result;
}
struct Signed {
    Nat magnitude;
    bool negative;
    Signed(Nat n,bool s=false):magnitude(std::move(n)),negative(s&&!magnitude.zero()) {}
    std::string dec() const {return std::string(negative?"-":"")+magnitude.dec();}
};
Signed subtract(Signed const& a,Signed const& b) {
    if (a.negative!=b.negative)
        return Signed(Nat::add(a.magnitude,b.magnitude),a.negative);
    const int cmp=a.magnitude.compare(b.magnitude);
    if (cmp>=0) return Signed(Nat::sub(a.magnitude,b.magnitude),a.negative);
    return Signed(Nat::sub(b.magnitude,a.magnitude),!a.negative);
}
std::string unsigned_sub(Nat const& a,Nat const& b) {
    return a.compare(b)<0 ? "-" : Nat::sub(a,b).hex();
}
std::string integer_b66() {
    std::string out=
        "# B66: independent exact integer exponent oracle, no runtime GMP and no implementation-under-test.\n"
        "# original_word_count mode exponent signed_source_decimal signed_result_decimal\n"
        "# Original std::mt19937_64(0xB660F00D) stream consumed for 12 modes at each width.\n"
        "# Modes 0,1,2,3 retain zero, saturated negative, single-top-bit positive and dense negative.\n"
        "# All 18 original exponents kept, including 0^0, 31/32/33, 63/64/65 boundaries.\n"
        "# Offline expected result generated with exact arbitrary-precision signed exponentiation.\n";
    std::mt19937_64 rng(0xB660F00DULL);
    for (std::size_t count:{0u,1u,2u,3u,4u}) {
        for (unsigned mode=0u;mode<12u;++mode) {
            Words words(count);
            for (auto& v:words) v=rng();
            if (mode%4u==0u) std::fill(words.begin(),words.end(),0u);
            else if (mode%4u==1u) std::fill(words.begin(),words.end(),UINT64_MAX);
            else if (mode%4u==2u && !words.empty()) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
            }
            if (mode>=4u) continue;
            const Nat base(std::move(words));
            const Signed original(base,(mode&1u)!=0u);
            for (unsigned exponent:{0u,1u,2u,3u,5u,7u,8u,9u,12u,15u,16u,17u,
                                    31u,32u,33u,63u,64u,65u}) {
                const Signed result(power(base,exponent),original.negative && (exponent&1u)!=0u);
                out+=row({std::to_string(count),std::to_string(mode),std::to_string(exponent),
                          original.dec(),result.dec()});
            }
        }
    }
    return out;
}
std::string integer_sub() {
    std::string out=
        "# Independent B41 signed subtraction fixtures (C++ mt19937_64 seed 0xB41D1FF).\n"
        "# width structural_mode sign_mask original_lhs_decimal original_rhs_decimal expected_difference_decimal\n"
        "# Six distinct structural classes, all sign pairs, all eight operand widths.\n";
    std::mt19937_64 rng(0xB41D1FFULL);
    for (std::size_t count:{1u,2u,7u,8u,16u,32u,64u,128u}) {
        for (unsigned mode=0u;mode<24u;++mode) {
            Words left(count),right(count);
            for (auto& x:left) x=rng();
            for (auto& x:right) x=rng();
            switch (mode%6u) {
                case 0u: right=left;break;
                case 1u: right=left;left[0u]=UINT64_MAX;right[0u]=0u;break;
                case 2u: right=left;left[0u]=0u;right[0u]=UINT64_MAX;break;
                case 3u: std::fill(left.begin(),left.end(),0u);break;
                case 4u: right.back()=0u;break;
                default: break;
            }
            if (mode>=6u) continue;
            const Nat a(left),b(right);
            for (unsigned signs=0u;signs<4u;++signs) {
                const Signed lhs(a,(signs&1u)!=0u),rhs(b,(signs&2u)!=0u);
                out+=row({std::to_string(count),std::to_string(mode),std::to_string(signs),
                    lhs.dec(),rhs.dec(),subtract(lhs,rhs).dec()});
            }
        }
    }
    return out;
}
std::string big_unsigned_sub() {
    std::string out=
        "# B41: source-keyed unsigned subtraction oracle; exact independent integer arithmetic.\n"
        "# width class original_left_hex original_right_hex left_minus_right_hex right_minus_left_hex\n"
        "# '-' is expected underflow. 0-width case and nine classes at each other width.\n"
        "# Samples extracted from the original mt19937_64(0xB41CA11CE) stream; no GMP/BigUnsigned for expected values.\n";
    std::mt19937_64 rng(0xB41CA11CEULL);
    for (std::size_t count:{0u,1u,2u,7u,8u,16u,31u,32u,63u,64u,65u,128u}) {
        for (unsigned mode=0u;mode<64u;++mode) {
            Words left(count),right(count);
            for (auto& x:left) x=rng();
            for (auto& x:right) x=rng();
            switch(mode%9u) {
                case 0u: right=left;break;
                case 1u: right=left;if(count!=0u)left[0u]=UINT64_MAX;break;
                case 2u: right=left;if(count!=0u)right[0u]=UINT64_MAX;break;
                case 3u: right=left;if(count!=0u){left[0u]=UINT64_MAX;right[0u]=0u;}break;
                case 4u:
                    std::fill(left.begin(),left.end(),UINT64_MAX);
                    std::fill(right.begin(),right.end(),UINT64_MAX);
                    if(count!=0u) right.front()=0u;
                    break;
                case 5u: if(count!=0u)left.back()=0u;break;
                case 6u: if(count!=0u)right.back()=0u;break;
                case 7u:
                    std::fill(left.begin(),left.end(),0u);
                    std::fill(right.begin(),right.end(),0u);
                    if(count!=0u){
                        left.front()=UINT64_MAX;
                        left.back()|=std::uint64_t{1}<<63u;
                        right.front()=1u;
                        right.back()|=std::uint64_t{1}<<63u;
                    }
                    break;
                default:break;
            }
            if(mode>=(count==0u?1u:9u)) continue;
            const Nat a(std::move(left)),b(std::move(right));
            out+=row({std::to_string(count),std::to_string(mode),a.hex(),b.hex(),
                unsigned_sub(a,b),unsigned_sub(b,a)});
        }
    }
    return out;
}
std::string integer_b33() {
    std::string out=
        "# Independent exact signed square oracle B33, original mt19937_64 seed 0xB33A11A5.\n"
        "# width mode source_signed_decimal expected_square_decimal.\n"
        "# Eight distinct structural modes per width (0..7). Identity, distinct-copy and *= aliasing are all checked.\n";
    std::mt19937_64 rng(0xB33A11A5ULL);
    for (std::size_t count:{1u,7u,8u,9u,16u,17u,32u,64u,128u}) {
        for(unsigned mode=0u;mode<16u;++mode) {
            Words words(count);
            for(auto& x:words)x=rng();
            switch(mode%8u) {
                case 0u: std::fill(words.begin(),words.end(),0u);break;
                case 1u: std::fill(words.begin(),words.end(),UINT64_MAX);break;
                case 2u: std::fill(words.begin(),words.end(),0u);words.back()=std::uint64_t{1}<<63u;break;
                case 3u: std::fill(words.begin(),words.end(),0u);words.front()=1u;words.back()|=1u;break;
                case 4u: for(std::size_t i=1u;i+1u<count;++i)words[i]=0u;break;
                default:words.back()|=std::uint64_t{1}<<63u;break;
            }
            if(mode>=8u)continue;
            const Nat value(std::move(words));
            const Signed original(value,(mode%2u)!=0u);
            out+=row({std::to_string(count),std::to_string(mode),original.dec(),multiply(value,value).dec()});
        }
    }
    return out;
}
std::string integer_b37() {
    std::string out=
        "# Independent exact signed power oracle B37, original mt19937_64 seed 0xB37C03.\n"
        "# width repetition exponent source_signed_decimal result_signed_decimal.\n"
        "# Four structural modes per width: saturated+, high-bit-only-, dense+, dense-; exponents 2,3,5.\n";
    std::mt19937_64 rng(0xB37C03ULL);
    for (std::size_t count:{7u,8u,9u,16u,32u,64u}) {
        for(unsigned rep=0u;rep<12u;++rep) {
            Words words(count);
            for(auto& x:words)x=rng();
            if(rep%5u==0u) std::fill(words.begin(),words.end(),UINT64_MAX);
            else if(rep%5u==1u){std::fill(words.begin(),words.end(),0u);words.back()=std::uint64_t{1}<<63u;}
            else words.back()|=std::uint64_t{1}<<63u;
            if(rep>=4u)continue;
            const Nat value(std::move(words));
            const Signed original(value,(rep%2u)!=0u);
            for(unsigned exponent:{2u,3u,5u}) {
                out+=row({std::to_string(count),std::to_string(rep),std::to_string(exponent),
                    original.dec(),Signed(power(value,exponent),original.negative && (exponent%2u)!=0u).dec()});
            }
        }
    }
    return out;
}
std::string integer_b38() {
    std::string out=
        "# Independent exact signed pow3 oracle B38, original mt19937_64 seed 0xB38C03.\n"
        "# width repetition source_signed_decimal result_signed_decimal.\n"
        "# Six structural modes per width: top-bit, top+13, alternating zero, saturated, two dense; both signs.\n";
    std::mt19937_64 rng(0xB38C03ULL);
    for (std::size_t count:{8u,16u,32u,64u,128u}) {
        for(unsigned rep=0u;rep<24u;++rep) {
            Words words(count);
            for(auto& x:words)x=rng();
            if(rep%6u==0u){std::fill(words.begin(),words.end(),0u);words.back()=std::uint64_t{1}<<63u;}
            else if(rep%6u==1u){std::fill(words.begin(),words.end(),0u);words.front()=13u;words.back()=std::uint64_t{1}<<63u;}
            else if(rep%6u==2u){for(std::size_t i=0u;i<count;i+=2u)words[i]=0u;}
            else if(rep%6u==3u)std::fill(words.begin(),words.end(),UINT64_MAX);
            words.back()|=std::uint64_t{1}<<63u;
            if(rep>=6u)continue;
            const Nat value(std::move(words));
            const Signed original(value,(rep%2u)!=0u);
            out+=row({std::to_string(count),std::to_string(rep),original.dec(),
                Signed(power(value,3u),original.negative).dec()});
        }
    }
    return out;
}
}
int main(int argc,char**argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"integer_b66.tsv",integer_b66(),write);
        independent_oracle::output(dir,"integer_sub.tsv",integer_sub(),write);
        independent_oracle::output(dir,"big_unsigned_sub.tsv",big_unsigned_sub(),write);
        independent_oracle::output(dir,"integer_b33.tsv",integer_b33(),write);
        independent_oracle::output(dir,"integer_b37.tsv",integer_b37(),write);
        independent_oracle::output(dir,"integer_b38.tsv",integer_b38(),write);
        return 0;
    } catch (std::exception const& e) {
        std::cerr<<"Structural generator: "<<e.what()<<'\n';
        return 1;
    }
}
