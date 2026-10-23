// SPDX-License-Identifier: MIT
// Reproduce B76's 2,800 exact golden GCD + decimal oracle cases.
// Standalone C++20; no Numeric, GMP, MPFR, Boost or external runtime.
#include "oracle_exact_integer.hpp"
using independent_oracle::Nat;

int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        std::mt19937_64 engine(0xB16B00B5ULL);
        constexpr unsigned size=200;
        for(unsigned part=0;part<14;++part) {
            std::ostringstream golden;
            golden<<"# B76: exact frozen sources and independent gcd/decimal oracles for every original random iteration.\n"
                  <<"# repetition generated_n generated_m x_hex y_hex gcd_hex x_decimal\n"
                  <<"# Part "<<part<<"; "<<part*size<<".."<<(part+1)*size-1
                  <<"; original mt19937_64 seed 0xB16B00B5.\n";
            for(unsigned repetition=part*size;repetition<(part+1)*size;++repetition) {
                const std::size_t n=repetition%113u==0u?96u:engine()%20u;
                const std::size_t m=repetition%71u==0u?32u:engine()%12u;
                std::vector<std::uint64_t> left(n),right(m);
                for(auto& word:left) word=engine();
                for(auto& word:right) word=engine();
                if(repetition%17u==0u) std::fill(left.begin(),left.end(),0u);
                if(repetition%13u==0u) std::fill(left.begin(),left.end(),UINT64_MAX);
                if(repetition%11u==0u) std::fill(right.begin(),right.end(),0u);
                const Nat x(std::move(left)),y(std::move(right));
                golden<<independent_oracle::row({std::to_string(repetition),std::to_string(n),std::to_string(m),
                    x.hex(),y.hex(),Nat::gcd(x,y).hex(),x.dec()});
            }
            independent_oracle::output(dir,"big_unsigned_b76_gcd_decimal_"+std::to_string(part)+".tsv",golden.str(),write);
        }
        return 0;
    } catch(std::exception const& e) {std::cerr<<"B76 generator: "<<e.what()<<'\n';return 1;}
}
