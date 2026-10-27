// SPDX-License-Identifier: MIT
// Standalone exact decimal oracle: bounded signed integer arithmetic only.
#include "oracle_exact_integer.hpp"

namespace {
struct DecimalInput { std::int64_t coefficient; unsigned scale; };
std::int64_t ten_power(unsigned n) {
    std::int64_t result=1;
    while (n--!=0u) result*=10;
    return result;
}
DecimalInput canonical(std::int64_t coefficient,unsigned scale) {
    while (scale!=0u && coefficient%10==0) { coefficient/=10; --scale; }
    return {coefficient,scale};
}
std::string fixture() {
    // |coefficient| <= 1250, scale <= 3: all products and rescalings fit int64_t.
    const std::array<DecimalInput,8> inputs={{{0,0},{125,2},{1250,3},{-125,2},
                                            {3,0},{-7,1},{1,3},{-1,3}}};
    std::string result="# Independent exact decimal arithmetic.\n"
        "# id lhs_coefficient lhs_scale rhs_coefficient rhs_scale sum_coefficient sum_scale product_coefficient product_scale unequal\n";
    unsigned id=0;
    for (auto a:inputs) for (auto b:inputs) {
        const unsigned scale=std::max(a.scale,b.scale);
        const auto left=a.coefficient*ten_power(scale-a.scale);
        const auto right=b.coefficient*ten_power(scale-b.scale);
        const auto sum=canonical(left+right,scale);
        const auto product=canonical(a.coefficient*b.coefficient,a.scale+b.scale);
        result+=independent_oracle::row({std::to_string(id++),std::to_string(a.coefficient),
            std::to_string(a.scale),std::to_string(b.coefficient),std::to_string(b.scale),
            std::to_string(sum.coefficient),std::to_string(sum.scale),
            std::to_string(product.coefficient),std::to_string(product.scale),left!=right?"1":"0"});
    }
    return result;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"decimal_exact.tsv",fixture(),write);
    } catch (std::exception const& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
