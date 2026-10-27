// SPDX-License-Identifier: MIT
// Exact unary binary64 wrappers. Integer bit-field transformations only.
#include "oracle_exact_integer.hpp"
namespace {
std::string fixture() {
    constexpr std::uint64_t sign=std::uint64_t{1}<<63u;
    constexpr std::uint64_t unit=std::uint64_t{1}<<52u;
    std::string result="# Exact binary64 unary wrappers; expected IEEE bit patterns.\n# input zero positive negative half\n";
    // Zero, 1/2, 3/2, and 2^60. All halves are exactly representable.
    for (std::uint64_t magnitude:{std::uint64_t{0},std::uint64_t{1022}*unit,
                                  std::uint64_t{1023}*unit+unit/2u,std::uint64_t{1083}*unit}) {
        for (std::uint64_t negative:{std::uint64_t{0},sign}) {
            const auto input=magnitude|negative;
            const auto half=(magnitude==0u?0u:magnitude-unit)|negative;
            const auto hex=[](std::uint64_t bits) {
                std::ostringstream out; out<<std::hex<<bits; return out.str();
            };
            result+=independent_oracle::row({hex(input),"0",hex(input),hex(input^sign),hex(half)});
        }
    }
    return result;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"floatdp_wrappers.tsv",fixture(),write);
    } catch (std::exception const& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
