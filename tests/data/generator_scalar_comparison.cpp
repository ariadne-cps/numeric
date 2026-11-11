/***************************************************************************
 *            generator_scalar_comparison.cpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Numeric.
 *
 *  Ariadne Numeric is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Numeric is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric. If not, see <https://www.gnu.org/licenses/>.
 */

#include "oracle_exact_integer.hpp"
namespace {
std::string fixture() {
    std::string result="# Independent Int/Rational comparisons.\n# lhs rhs_numerator rhs_denominator comparison\n";
    for (int lhs:{-1,0,1}) {
        for (auto [n,d]:std::array<std::pair<int,int>,9>{{{-3,2},{0,1},{1,2},{1,1},
                                                       {3,2},{-1,1},{0,0},{1,0},{-1,0}}}) {
            const int order=d==0 ? (n==0?-128:n>0?-1:1) :
                lhs*d<n?-1:lhs*d>n?1:0;
            result+=independent_oracle::row({std::to_string(lhs),std::to_string(n),
                                            std::to_string(d),std::to_string(order)});
        }
    }
    return result;
}
}
int main(int argc,char** argv) {
    try {
        bool write=false;
        const auto dir=independent_oracle::data_directory(argc,argv,write);
        independent_oracle::output(dir,"scalar_comparison.tsv",fixture(),write);
    } catch (std::exception const& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
