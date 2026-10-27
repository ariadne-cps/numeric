// SPDX-License-Identifier: MIT
// Exact small dyadic constructor values; bounded integer reduction only.
#include "oracle_exact_integer.hpp"
#include <numeric>
namespace {
std::string fixture() {
    std::string result="# Exact FloatMP construction oracle: all sources fit in eight significant bits.\n"
        "# id source_numerator source_denominator precision_bits expected_numerator expected_denominator\n";
    unsigned id=0u;
    for (int n:{-129,-5,-1,0,1,5,129}) for (int d:{1,4}) for (unsigned p:{8u,53u,512u}) {
        const int divisor=std::gcd(n,d);
        result+=independent_oracle::row({std::to_string(id++),std::to_string(n),std::to_string(d),
            std::to_string(p),std::to_string(n/divisor),std::to_string(d/divisor)});
    }
    return result;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"floatmp_construction.tsv",fixture(),write);
    } catch (std::exception const& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
