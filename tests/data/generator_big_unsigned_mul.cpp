// SPDX-License-Identifier: MIT
// Reproduce B35/B37/B38 multiplication and A1 modular accumulation cases.
// Standalone C++20; no Numeric, GMP, independent exact arithmetic, Boost or independent exact arithmetic.
#include "oracle_exact_integer.hpp"

namespace {
using independent_oracle::Nat;
using independent_oracle::row;
using Words=std::vector<std::uint64_t>;
using Pair=std::pair<std::size_t,std::size_t>;

// Binary shift/add avoids Numeric's multiplication, squaring and low-mul kernels.
Nat multiply(Nat left, Nat right) {
    Nat result;
    while (!right.zero()) {
        if ((right.w.front()&1u)!=0u) result=Nat::add(result,left);
        right=right.shr(1u);
        if (!right.zero()) left=left.shl(1u);
    }
    return result;
}
Nat low_product(Nat value,std::size_t limbs) {
    if(value.w.size()>limbs) value.w.resize(limbs);
    value.trim();
    return value;
}
// A2: reproduce all 11*192 original B32 square source visits exactly.
// Each mt19937_64 draw is consumed before applying the original transforms;
// expected squares use the independent shift/add exact multiplication above.
std::string b32() {
    std::string out=
        "# B32 independent exact square oracle; all 2,112 original source visits.\n"
        "# count repetition input_hex square_hex\n"
        "# Original mt19937_64(0xB3205A7E) stream; all draws preserved, even when overwritten.\n"
        "# Expected squares computed with independent binary shift/add multiplication.\n"
        "# No Numeric, GMP, independent exact arithmetic, independent exact arithmetic or Boost.\n";
    std::mt19937_64 rng(0xB3205A7EULL);
    std::size_t rows=0u;
    for(std::size_t count:{1u,2u,7u,8u,9u,12u,16u,32u,64u,96u,128u}) {
        for(unsigned repetition=0u;repetition<192u;++repetition) {
            Words words(count);
            for(auto& word:words)word=rng();
            if(repetition%13u==0u)std::fill(words.begin(),words.end(),UINT64_MAX);
            else if(repetition%17u==0u) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
            } else if(repetition%19u==0u) {
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()=UINT64_MAX;
            } else if(repetition%23u==0u) {
                std::fill(words.begin(),words.end(),0u);
            }
            const Nat input(std::move(words));
            out+=row({std::to_string(count),std::to_string(repetition),
                      input.hex(),multiply(input,input).hex()});
            ++rows;
        }
    }
    if(rows!=2112u)throw std::runtime_error("wrong B32 row count");
    return out;
}

std::string b35() {
    std::string out=
      "# B35 frozen independent exact Karatsuba product oracle.\n"
      "# n m repetition lhs_hex rhs_hex product_hex\n"
      "# Sources from std::mt19937_64(0xB35CA7A5); cross-checked against separately compiled C++ first five draws.\n"
      "# All 14 original width pairs retained. Original 32 random draws per width pair consumed before selecting original modes 0..4.\n"
      "# Structural classes: all-ones/carry, sparse low-only, sparse endpoints, alternate-zero rows, dense random.\n";
    std::mt19937_64 rng(0xB35CA7A5ULL);
    for(const auto& [n,m]:std::initializer_list<Pair>{{63u,63u},{63u,64u},{64u,63u},{64u,64u},
        {64u,65u},{65u,64u},{65u,65u},{80u,64u},{96u,64u},{96u,96u},
        {128u,64u},{128u,96u},{96u,128u},{128u,128u}}) {
        for(unsigned mode=0u;mode<32u;++mode) {
            Words left(n),right(m);
            for(auto& x:left) x=rng();
            for(auto& x:right) x=rng();
            if(mode%7u==0u) {
                std::fill(left.begin(),left.end(),UINT64_MAX);
                std::fill(right.begin(),right.end(),UINT64_MAX);
            } else if(mode%7u==1u) {
                std::fill(left.begin(),left.end(),0u);
                std::fill(right.begin(),right.end(),0u);
                left.front()=13u;
                right.front()=31u;
            } else if(mode%7u==2u) {
                for(std::size_t i=1u;i+1u<n;++i) left[i]=0u;
                for(std::size_t i=1u;i+1u<m;++i) right[i]=0u;
            } else if(mode%7u==3u) {
                for(std::size_t i=0u;i<n;i+=2u) left[i]=0u;
                for(std::size_t i=0u;i<m;i+=2u) right[i]=0u;
            }
            left.back()|=std::uint64_t{1}<<63u;
            right.back()|=std::uint64_t{1}<<63u;
            if(mode>=5u) continue;
            const Nat lhs(std::move(left)),rhs(std::move(right));
            out+=row({std::to_string(n),std::to_string(m),std::to_string(mode),
                lhs.hex(),rhs.hex(),multiply(lhs,rhs).hex()});
        }
    }
    return out;
}
std::string b37() {
    std::string out=
      "# B37 exact independent unsigned squares, original mt19937_64(0xB37A11) input stream.\n"
      "# width original_repetition input_hex square_hex\n"
      "# All 13 limb widths and 7 distinct repetition%7 classes; 1248 -> 91.\n";
    std::mt19937_64 rng(0xB37A11ULL);
    for(std::size_t count:{1u,4u,7u,8u,9u,15u,16u,17u,32u,33u,64u,65u,128u}) {
        for(unsigned mode=0u;mode<96u;++mode) {
            Words words(count);
            for(auto& x:words) x=rng();
            switch(mode%7u) {
                case 0u: std::fill(words.begin(),words.end(),UINT64_MAX);break;
                case 1u: std::fill(words.begin(),words.end(),0u);words.back()=std::uint64_t{1}<<63u;break;
                case 2u: for(std::size_t i=1u;i+1u<count;++i)words[i]=0u;break;
                case 3u: for(std::size_t i=0u;i<count;i+=2u)words[i]=0u;break;
                default: break;
            }
            if(mode>=7u) continue;
            const Nat value(std::move(words));
            out+=row({std::to_string(count),std::to_string(mode),value.hex(),multiply(value,value).hex()});
        }
    }
    return out;
}
std::string b38() {
    std::string out=
      "# B38 independently frozen full and six truncated product references; original mt19937_64(0xB38C0FFEE) source.\n"
      "# n m repetition lhs_hex rhs_hex full_hex low0 low1 low_n_div_2 low_n low_n_plus_m_minus_1 low_n_plus_m_plus_2\n"
      "# All 14 original width pairs and eight distinct repetition%8 classes (336 -> 112).\n";
    std::mt19937_64 rng(0xB38C0FFEEULL);
    for(const auto& [n,m]:std::initializer_list<Pair>{{1u,3u},{4u,7u},{7u,8u},{8u,8u},
        {9u,16u},{16u,32u},{32u,64u},{64u,32u},{64u,128u},{128u,64u},
        {64u,129u},{129u,64u},{96u,128u},{128u,96u}}) {
        for(unsigned mode=0u;mode<24u;++mode) {
            Words left(n),right(m);
            for(auto& x:left) x=rng();
            for(auto& x:right) x=rng();
            switch(mode%8u) {
                case 0u:
                    std::fill(left.begin(),left.end(),0u);
                    left.back()=std::uint64_t{1}<<63u;
                    break;
                case 1u:
                    std::fill(left.begin(),left.end(),0u);
                    left.front()=13u;
                    left.back()=std::uint64_t{1}<<63u;
                    break;
                case 2u:
                    std::fill(right.begin(),right.end(),0u);
                    right.front()=17u;
                    right.back()=std::uint64_t{1}<<63u;
                    break;
                case 3u:
                    std::fill(left.begin(),left.end(),0u);
                    std::fill(right.begin(),right.end(),0u);
                    left.front()=UINT64_MAX;
                    right.front()=UINT64_MAX;
                    left.back()=UINT64_MAX;
                    right.back()=UINT64_MAX;
                    break;
                case 4u: for(std::size_t i=0u;i<n;i+=2u)left[i]=0u;break;
                case 5u:
                    std::fill(left.begin(),left.end(),UINT64_MAX);
                    std::fill(right.begin(),right.end(),UINT64_MAX);
                    break;
                default:break;
            }
            left.back()|=std::uint64_t{1}<<63u;
            right.back()|=std::uint64_t{1}<<63u;
            if(mode>=8u) continue;
            const Nat lhs(std::move(left)),rhs(std::move(right));
            const Nat product=multiply(lhs,rhs);
            const std::array<std::size_t,6u> lengths={0u,1u,n/2u,n,n+m-1u,n+m+2u};
            std::vector<std::string> columns={std::to_string(n),std::to_string(m),
                std::to_string(mode),lhs.hex(),rhs.hex(),product.hex()};
            for(std::size_t k:lengths) columns.push_back(low_product(product,k).hex());
            out+=row(columns);
        }
    }
    return out;
}


std::string a25_b52() {
    std::string out=
        "# A25 B52 multiply-kernel oracle: all 288 original n/m/mode visits.\n"
        "# n m mode lhs_hex rhs_hex product_hex low_limbs low_product_hex\n"
        "# Exact mt19937_64(0xB5203232) source stream; independent exact integer multiplication.\n"
        "# No Numeric/GMP/Boost dependency.\n";
    std::mt19937_64 rng(0xB5203232ULL);
    std::size_t rows=0u;
    for (const std::size_t n:{31u,32u,33u}) {
        for (const std::size_t m:{31u,32u,33u}) {
            for (unsigned mode=0u;mode<32u;++mode) {
                Words left(n),right(m);
                for (auto& digit:left) digit=rng();
                for (auto& digit:right) digit=rng();
                if (mode%8u==0u) {
                    std::fill(left.begin(),left.end(),UINT64_MAX);
                    std::fill(right.begin(),right.end(),UINT64_MAX);
                } else if (mode%8u==1u) {
                    std::fill(left.begin(),left.end(),0u);
                    std::fill(right.begin(),right.end(),0u);
                    left.back()=1u;
                    right.back()=1u;
                } else if (mode%8u==2u) {
                    for (std::size_t i=0u;i<n-1u;i+=2u) left[i]=0u;
                    for (std::size_t i=0u;i<m-1u;i+=2u) right[i]=0u;
                } else if (mode%8u==3u) {
                    std::fill(left.begin(),left.end(),0u);
                    std::fill(right.begin(),right.end(),0u);
                    left.front()=UINT64_MAX;
                    right.front()=UINT64_MAX;
                } else if (mode%8u==4u) {
                    for (std::size_t i=1u;i+1u<n;++i) left[i]=0u;
                    for (std::size_t i=1u;i+1u<m;++i) right[i]=0u;
                }
                const Nat lhs(std::move(left)),rhs(std::move(right));
                const Nat product=multiply(lhs,rhs);
                const std::size_t full=lhs.w.size()+rhs.w.size();
                const std::size_t low=full==0u?0u:full-1u;
                out+=row({std::to_string(n),std::to_string(m),
                          std::to_string(mode),lhs.hex(),rhs.hex(),
                          product.hex(),std::to_string(low),
                          low_product(product,low).hex()});
                ++rows;
            }
        }
    }
    if (rows!=288u) throw std::runtime_error("wrong A25 B52 row count");
    return out;
}

// A1: preserve every original 9*32*3*2 seeded product-accumulation visit.
// Expected values use exact signed arithmetic modulo the entire destination
// width, not any BigUnsigned accumulation or multiplication implementation.
std::string accumulation() {
    std::string out=
        "# Independent BigUnsigned accumulate_product_into exact modular oracle, A1.\n"
        "# length trial extra subtract left_hex right_hex initial_hex result_hex escaped\n"
        "# Inputs replay all original mt19937_64(0x23ACC00) draws, including overwritten ones.\n"
        "# Product uses independent radix-2^32 arithmetic; escape and low result modulo 2^(64*count).\n"
        "# Each 9 x 32 x 3 x 2 case preserves its original requested buffer count and subtraction direction.\n"
        "# No Numeric, GMP, independent exact arithmetic, independent exact arithmetic or Boost linked to this standalone C++20 generator.\n";
    std::mt19937_64 random(0x23ACC00ULL);
    std::size_t count_rows=0u;
    for(std::size_t length:{0u,1u,2u,4u,5u,8u,9u,16u,64u}) {
        for(unsigned trial=0u;trial<32u;++trial) {
            const std::size_t other=trial%3u==0u?1u:length;
            Words x(length),y(other);
            for(auto& word:x)word=random();
            for(auto& word:y)word=random();
            if(trial%4u==0u)std::fill(x.begin(),x.end(),UINT64_MAX);
            if(trial%4u==1u) {
                std::fill(x.begin(),x.end(),0u);
                if(!x.empty()){x.back()=1u;x.front()|=1u;}
            }
            if(trial%4u==2u)std::fill(y.begin(),y.end(),0u);
            const Nat left(x),right(y);
            const Nat full=multiply(left,right);
            for(std::size_t extra:{0u,1u,3u}) {
                const std::size_t count=left.w.size()+right.w.size()+extra;
                Words original(count);
                for(auto& word:original)word=random();
                if(trial%3u==0u)std::fill(original.begin(),original.end(),0u);
                if(trial%3u==1u)std::fill(original.begin(),original.end(),UINT64_MAX);
                const Nat initial(original),modulus=Nat::power2(64u*count);
                for(bool subtract:{false,true}) {
                    Nat expected;
                    bool escaped=false;
                    if(!subtract) {
                        const Nat sum=Nat::add(initial,full);
                        escaped=sum.compare(modulus)>=0;
                        expected=low_product(sum,count);
                    } else if(initial.compare(full)>=0) {
                        expected=low_product(Nat::sub(initial,full),count);
                    } else {
                        escaped=true;
                        const Nat deficit=low_product(Nat::sub(full,initial),count);
                        expected=deficit.zero()?Nat():Nat::sub(modulus,deficit);
                    }
                    out+=row({std::to_string(length),std::to_string(trial),
                              std::to_string(extra),subtract?"1":"0",
                              left.hex(),right.hex(),initial.hex(),
                              expected.hex(),escaped?"1":"0"});
                    ++count_rows;
                }
            }
        }
    }
    if(count_rows!=1728u)throw std::runtime_error("wrong accumulation case count");
    return out;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"big_unsigned_b35.tsv",b35(),write);
        independent_oracle::output(dir,"big_unsigned_b37.tsv",b37(),write);
        independent_oracle::output(dir,"big_unsigned_b38.tsv",b38(),write);
        independent_oracle::output(dir,"big_unsigned_b32.tsv",b32(),write);
        independent_oracle::output(dir,"big_unsigned_a25_b52.tsv",a25_b52(),write);
        independent_oracle::output(dir,"big_unsigned_accumulation.tsv",accumulation(),write);
        return 0;
    } catch(std::exception const& ex) {
        std::cerr<<"BigUnsigned multiplication/accumulation generator: "<<ex.what()<<'\n';
        return 1;
    }
}
