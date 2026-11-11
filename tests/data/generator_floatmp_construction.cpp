/***************************************************************************
 *            generator_floatmp_construction.cpp
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

#include <numeric>

#include "oracle_exact_integer.hpp"
namespace {
std::string fixture() {
    std::string result = "# Exact FloatMP construction oracle: all sources fit in eight significant bits.\n"
                         "# id source_numerator source_denominator precision_bits expected_numerator expected_denominator\n";
    unsigned id = 0u;
    for (int n: {-129, -5, -1, 0, 1, 5, 129})
        for (int d: {1, 4})
            for (unsigned p: {8u, 53u, 512u}) {
                const int divisor = std::gcd(n, d);
                result += independent_oracle::row({std::to_string(id++), std::to_string(n), std::to_string(d),
                                                   std::to_string(p), std::to_string(n / divisor), std::to_string(d / divisor)});
            }
    return result;
}
} // namespace
int main(int argc, char** argv) {
    try {
        bool write = false;
        const auto dir = independent_oracle::data_directory(argc, argv, write);
        independent_oracle::output(dir, "floatmp_construction.tsv", fixture(), write);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
