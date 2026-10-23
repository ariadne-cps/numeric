// SPDX-License-Identifier: MIT
// Exact, standalone B72 multiplication/squaring fixture generator (C++20).
#include "oracle_exact_integer.hpp"
using independent_oracle::Nat;
using independent_oracle::row;
namespace {
Nat multiply(Nat a,Nat b) {
    Nat result;
    while(!b.zero()) {
        if((b.w.front()&1u)!=0u) result=Nat::add(result,a);
        b=b.shr(1u);
        if(!b.zero()) a=a.shl(1u);
    }
    return result;
}
std::string generate_b72() {
    std::string out=
      "# B72: independent exact 4/5-limb multiplication and squaring oracles.\n"
      "# n m mode canonical_lhs_hex canonical_rhs_hex canonical_product_hex canonical_lhs_square_hex\n"
      "# Preserve all 36 width pairs (0..5) and all eight mode%8 classes, from the original mt19937_64(0xB7204455) stream.\n"
      "# Nonselected modes still advance the original PRNG, preserving original source-keyed sample identity.\n";
    std::mt19937_64 rng(0xB7204455ULL);
    for(std::size_t n:{0u,1u,2u,3u,4u,5u})
    for(std::size_t m:{0u,1u,2u,3u,4u,5u})
    for(unsigned mode=0u;mode<32u;++mode) {
        std::array<std::uint64_t,5u> a{},b{};
        for(std::size_t i=0u;i<n;++i) a[i]=rng();
        for(std::size_t i=0u;i<m;++i) b[i]=rng();
        if(mode%8u==0u) {
            std::fill(a.begin(),a.end(),UINT64_MAX);
            std::fill(b.begin(),b.end(),UINT64_MAX);
        } else if(mode%8u==1u) {
            std::fill(a.begin(),a.end(),0u);
            std::fill(b.begin(),b.end(),0u);
        } else if(mode%8u==2u) {
            std::fill(a.begin(),a.end(),0u);
            std::fill(b.begin(),b.end(),0u);
            if(n!=0u) a[n-1u]=1u;
            if(m!=0u) b[m-1u]=1u;
        } else if(mode%8u==3u) {
            for(std::size_t i=0u;i<n;i+=2u) a[i]=0u;
            for(std::size_t i=0u;i<m;i+=2u) b[i]=0u;
        } else if(mode%8u==4u) {
            for(std::size_t i=1u;i<n;++i) a[i]=0u;
            for(std::size_t i=1u;i<m;++i) b[i]=0u;
        }
        if(mode>=8u) continue;
        const Nat lhs(std::vector<std::uint64_t>(a.begin(),a.begin()+static_cast<std::ptrdiff_t>(n)));
        const Nat rhs(std::vector<std::uint64_t>(b.begin(),b.begin()+static_cast<std::ptrdiff_t>(m)));
        out+=row({std::to_string(n),std::to_string(m),std::to_string(mode),lhs.hex(),rhs.hex(),
                  multiply(lhs,rhs).hex(),multiply(lhs,lhs).hex()});
    }
    return out;
}
std::string generate_r17() {
    std::string out=
      "# R17 independent signed fraction cross-product comparison oracle.\n"
      "# bits lhs_index rhs_index expected_sign (-1/0/1).\n";
    struct Probe { Nat numerator,denominator; bool negative; };
    for(unsigned bits:{128u,256u,448u,512u,1024u,2048u,4096u}) {
        const Nat base=Nat::power2(bits);
        const std::array<Probe,6u> p={{{Nat::add(base,Nat(13u)),Nat::add(base.mul_small(2u),Nat(17u)),false},
             {Nat::add(base,Nat(31u)),Nat::add(base.mul_small(3u),Nat(19u)),false},
             {Nat::add(base,Nat(13u)),Nat::add(base.mul_small(2u),Nat(19u)),false},
             {Nat::add(base,Nat(1u)),Nat::add(base,Nat(2u)),false},
             {Nat::add(base,Nat(2u)),Nat::add(base,Nat(3u)),false},
             {Nat::add(base,Nat(29u)),Nat::add(base.mul_small(3u),Nat(7u)),true}}};
        for(unsigned i=0u;i<p.size();++i)for(unsigned j=0u;j<p.size();++j) {
            int cmp=0;
            if(p[i].negative!=p[j].negative) cmp=p[i].negative?-1:1;
            else {
                cmp=multiply(p[i].numerator,p[j].denominator).compare(
                    multiply(p[j].numerator,p[i].denominator));
                if(p[i].negative) cmp=-cmp;
            }
            out+=row({std::to_string(bits),std::to_string(i),std::to_string(j),std::to_string(cmp)});
        }
    }
    return out;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"big_unsigned_b72.tsv",generate_b72(),write);
        independent_oracle::output(dir,"rational_r17.tsv",generate_r17(),write);
        return 0;
    } catch(std::exception const& e) {
        std::cerr<<"B72 generator: "<<e.what()<<'\n';
        return 1;
    }
}
