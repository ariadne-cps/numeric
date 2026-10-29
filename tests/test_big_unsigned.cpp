/***************************************************************************
 *            test_big_unsigned.cpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "numeric/big_unsigned.hpp"

#include <algorithm>
#include <bit>
#include <fstream>
#include <initializer_list>
#include <array>
#include <cstdint>
#include <cfenv>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "macro/testing.hpp"

namespace Ariadne {
BigUnsigned rem_limb_divq_diagnostic(BigUnsigned const&,BigUnsigned const&);
BigUnsigned rem_limb_preinverse_diagnostic(BigUnsigned const&,BigUnsigned const&);
BigUnsigned mul_short_stack_selected_diagnostic(BigUnsigned const&,BigUnsigned const&);
BigUnsigned mul_short_stack_fused_diagnostic(BigUnsigned const&,BigUnsigned const&);
BigUnsigned square_short_symmetric_diagnostic(BigUnsigned const&);
std::string to_decimal_cutoff16_diagnostic(BigUnsigned const&);
std::string to_decimal_cutoff32_diagnostic(BigUnsigned const&);
}

using namespace Ariadne;

namespace {

#ifndef ARIADNE_TEST_DATA_DIR
#error "Provide the absolute test data directory from CMake"
#endif

inline constexpr char ARIADNE_BIG_UNSIGNED_SHORT_DIV_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_short_div.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_DIV128_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_div128by64.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_DIV_LIMB_RANDOM_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_div_limb_random.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B48_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b48.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B49_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b49.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B57_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b57.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B61_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b61.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_EXACT_LIMB_SMALL_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_exact_limb_small.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_EXACT_LIMB_65_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_exact_limb_65.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_EXACT_LIMB_129A_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_exact_limb_129a.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_EXACT_LIMB_129B_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_exact_limb_129b.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B57_POWER_SHORT_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b57_power_short.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B57_POWER_MEDIUM_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b57_power_medium.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B57_POWER_LARGE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b57_power_large.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B55_SHORT_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b55_short.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B55_MEDIUM_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b55_medium.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B55_LARGE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b55_large.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B27_FIRST_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b27_first.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B27_SECOND_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b27_second.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B48_MULTIPRECISION_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b48_multiprecision.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B123_FIRST_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b123_randomized_first.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B124_SECOND_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b124_randomized_second.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B125_THIRD_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b125_randomized_third.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B126_FOURTH_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b126_randomized_fourth.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B127_FIFTH_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b127_randomized_fifth.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B128_SIXTH_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b128_randomized_sixth.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B129_SEVENTH_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b129_randomized_seventh.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B130_EIGHTH_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b130_randomized_eighth.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_R6_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_r6.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_R8_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_r8.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_R10_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_r10.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_R11_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_r11.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_R27_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_r27.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_R21_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_r21.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B76_GCD_ORACLE_PREFIX[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b76_gcd_decimal_";
inline constexpr char ARIADNE_BIG_UNSIGNED_SQRT_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_sqrt.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_ACCUMULATION_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_accumulation.tsv";

unsigned golden_hex_nibble(char c)
{
    if (c>='0' && c<='9') return static_cast<unsigned>(c-'0');
    if (c>='a' && c<='f') return static_cast<unsigned>(c-'a'+10);
    throw std::invalid_argument("Invalid hexadecimal BigUnsigned golden digit");
}

BigUnsigned big_unsigned_from_golden_hex(std::string const& hex)
{
    if (hex.empty() || (hex.size()>1u && hex.front()=='0')) {
        throw std::invalid_argument("Noncanonical BigUnsigned golden input");
    }
    std::vector<std::uint64_t> words((hex.size()+15u)/16u,0u);
    for (std::size_t i=0u;i<hex.size();++i) {
        words[i/16u] |= std::uint64_t{golden_hex_nibble(hex[hex.size()-1u-i])}
                        << static_cast<unsigned>(4u*(i%16u));
    }
    return BigUnsigned(words);
}

void check_against_precomputed_hex(BigUnsigned const& value, std::string const& hex)
{
    if (hex.empty() || (hex.size()>1u && hex.front()=='0')) {
        throw std::invalid_argument("Noncanonical BigUnsigned golden input");
    }
    const BigUnsigned expected=big_unsigned_from_golden_hex(hex);
    const std::size_t bits=(hex=="0") ? 0u :
        4u*(hex.size()-1u)+static_cast<std::size_t>(
            std::bit_width(golden_hex_nibble(hex.front())));
    const auto limbs=value.limbs();
    ARIADNE_TEST_ASSERT(compare(value,expected)==0);
    ARIADNE_TEST_ASSERT(value.is_zero()==(bits==0u));
    ARIADNE_TEST_ASSERT(value.bit_length()==bits);
    ARIADNE_TEST_ASSERT(limbs.empty() || limbs.back()!=0u);
}



void check_div128_golden(std::ifstream& golden,char kind,std::size_t ordinal,
                         std::uint64_t divisor,std::uint64_t low,
                         std::uint64_t high,BigUnsigned const& dividend,
                         BigUnsigned const& quotient,std::uint64_t remainder)
{
    std::string ref_kind,divisor_hex,low_hex,high_hex,source_hex;
    std::string quotient_hex,remainder_hex;
    std::size_t ref_ordinal=0u;
    if (!(golden >> ref_kind >> ref_ordinal >> divisor_hex >> low_hex
                 >> high_hex >> source_hex >> quotient_hex >> remainder_hex)
        || ref_kind!=std::string(1,kind) || ref_ordinal!=ordinal)
        throw std::runtime_error("Missing, malformed or unordered B114 div128 TSV row");
    check_against_precomputed_hex(BigUnsigned(divisor),divisor_hex);
    check_against_precomputed_hex(BigUnsigned(low),low_hex);
    check_against_precomputed_hex(BigUnsigned(high),high_hex);
    check_against_precomputed_hex(dividend,source_hex);
    check_against_precomputed_hex(quotient,quotient_hex);
    check_against_precomputed_hex(BigUnsigned(remainder),remainder_hex);
    ARIADNE_TEST_ASSERT(remainder<divisor);
    check_against_precomputed_hex(dividend,source_hex);
}


struct LimbDivisionGolden {
    std::ifstream file;
    explicit LimbDivisionGolden(char const* path):file(path) {
        if (!file) throw std::runtime_error("Missing B116 limb division TSV");
        std::string header;
        while (file.peek()=='#') std::getline(file,header);
    }
    void check(std::size_t count,std::uint64_t divisor,unsigned mode,
               BigUnsigned const& source,BigUnsigned const& quotient,
               std::uint64_t remainder) {
        std::size_t ref_count=0u;
        std::uint64_t ref_divisor=0u;
        unsigned ref_mode=0u;
        std::string source_hex,quotient_hex,remainder_hex;
        if (!(file >> ref_count >> ref_divisor >> ref_mode
                   >> source_hex >> quotient_hex >> remainder_hex) ||
            ref_count!=count || ref_divisor!=divisor || ref_mode!=mode)
            throw std::runtime_error("Missing or unordered B116 limb division row");
        check_against_precomputed_hex(source,source_hex);
        check_against_precomputed_hex(quotient,quotient_hex);
        check_against_precomputed_hex(BigUnsigned(remainder),remainder_hex);
        ARIADNE_TEST_ASSERT(remainder<divisor);
        check_against_precomputed_hex(source,source_hex);
    }
    void finish() {
        std::string extra;
        if (file >> extra) throw std::runtime_error("Extra B116 limb division rows");
    }
};

inline constexpr char ARIADNE_BIG_UNSIGNED_B54_PRODUCT_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b54_product.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B54_SQUARE_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b54_square.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B32_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b32.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B37_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b37.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B38_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b38.tsv";

inline constexpr char ARIADNE_BIG_UNSIGNED_B39_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b39.tsv";

inline constexpr char ARIADNE_BIG_UNSIGNED_B40_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b40.tsv";

inline constexpr char ARIADNE_BIG_UNSIGNED_B35_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b35.tsv";

inline constexpr char ARIADNE_BIG_UNSIGNED_RANDOM_MUL_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_random_mul.tsv";

inline constexpr char ARIADNE_BIG_UNSIGNED_R14_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_r14.tsv";
inline constexpr char ARIADNE_BIG_UNSIGNED_B72_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_b72.tsv";

inline constexpr char ARIADNE_BIG_UNSIGNED_SUB_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_sub.tsv";

inline constexpr char ARIADNE_BIG_UNSIGNED_POWER_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_power.tsv";

void check_against_power_golden(std::ifstream& golden,char kind,
                                std::size_t param1,std::size_t param2,
                                unsigned exponent,
                                BigUnsigned const& input,
                                BigUnsigned const& actual)
{
    char row_kind='?';
    std::size_t row_param1=0u,row_param2=0u;
    unsigned row_exponent=0u;
    std::string original_hex,expected_hex;
    if (!(golden >> row_kind >> row_param1 >> row_param2
                 >> row_exponent >> original_hex >> expected_hex)) {
        throw std::runtime_error("Missing or malformed BigUnsigned power row");
    }
    if (row_kind!=kind || row_param1!=param1 || row_param2!=param2
        || row_exponent!=exponent) {
        throw std::runtime_error("Unordered BigUnsigned power golden row");
    }
    check_against_precomputed_hex(actual,expected_hex);
    check_against_precomputed_hex(input,original_hex);
}

// Source-keyed independent exact GCD, including source immutability,
// canonical representation and every original bit-width boundary.
void require_gcd_golden(std::ifstream& golden,char const* group) {
    if (!golden)
        throw std::runtime_error(std::string("Missing BigUnsigned ")+group+" GCD TSV");
    std::string header;
    while (golden.peek()=='#') std::getline(golden,header);
}
void ensure_gcd_golden_end(std::ifstream& golden,char const* group) {
    std::string extra;
    if (golden >> extra)
        throw std::runtime_error(std::string("Extra BigUnsigned ")+group+" GCD TSV rows");
}
void check_gcd_golden(std::ifstream& golden,
                      std::initializer_list<std::string> keys,
                      BigUnsigned const& left,BigUnsigned const& right) {
    for (auto const& expected_key:keys) {
        std::string key;
        if (!(golden >> key) || key!=expected_key)
            throw std::runtime_error("Missing or unordered BigUnsigned GCD TSV row");
    }
    std::string lhs,rhs,want;
    if (!(golden >> lhs >> rhs >> want))
        throw std::runtime_error("Malformed BigUnsigned GCD TSV row");
    check_against_precomputed_hex(left,lhs);
    check_against_precomputed_hex(right,rhs);
    check_against_precomputed_hex(gcd(left,right),want);
    check_against_precomputed_hex(gcd(right,left),want);
    check_against_precomputed_hex(left,lhs);
    check_against_precomputed_hex(right,rhs);
}


class TestBigUnsigned
{
  public:
    void test() {
        ARIADNE_TEST_CALL(test_storage());
        ARIADNE_TEST_CALL(test_boundaries());
        ARIADNE_TEST_CALL(test_subtraction_common_prefix());
        ARIADNE_TEST_CALL(test_randomized());
        ARIADNE_TEST_CALL(test_multiplication());
        ARIADNE_TEST_CALL(test_division_limb());
        ARIADNE_TEST_CALL(test_exact_division_limb());
        ARIADNE_TEST_CALL(test_product_accumulation());
        ARIADNE_TEST_CALL(test_division());
        ARIADNE_TEST_CALL(test_gcd_decimal());
        ARIADNE_TEST_CALL(test_views_and_scratch());
        ARIADNE_TEST_CALL(test_multiply_kernel_paths());
        ARIADNE_TEST_CALL(test_power());
        ARIADNE_TEST_CALL(test_sqrt());
        ARIADNE_TEST_CALL(test_short_owned_division());
    }

  private:
    void test_storage();
    void test_boundaries();
    void test_subtraction_common_prefix();
    void test_randomized();
    void test_multiplication();
    void test_division_limb();
    void test_exact_division_limb();
    void test_product_accumulation();
    void test_division();
    void test_gcd_decimal();
    void test_views_and_scratch();
    void test_multiply_kernel_paths();
    void test_power();
    void test_sqrt();
    void test_short_owned_division();
};


void TestBigUnsigned::test_sqrt()
{
    // Exact squares whose floating seeds round on either side under
    // directed modes. The expected roots are integer identities.
    const auto saved_mode=std::fegetround();
    const std::array<BigUnsigned::Limb,2u> high_square={1u,0x4000000000000001ULL};
    const std::array<BigUnsigned::Limb,2u> max_square={1u,0xfffffffffffffffeULL};
    for (auto mode : {FE_TONEAREST,FE_DOWNWARD,FE_UPWARD,FE_TOWARDZERO}) {
        std::fesetround(mode);
        for (auto const& [input,expected] : {
                std::pair{BigUnsigned(0xfffffffe00000001ULL),std::uint64_t{0xffffffffULL}},
                std::pair{BigUnsigned(high_square),std::uint64_t{0x8000000000000001ULL}},
                std::pair{BigUnsigned(max_square),std::uint64_t{0xffffffffffffffffULL}}}) {
            const auto [root,remainder]=sqrtrem(input);
            ARIADNE_TEST_EQUALS(root.limb(0u),expected);
            ARIADNE_TEST_ASSERT(remainder.is_zero());
            ARIADNE_TEST_EQUALS(isqrt(input).limb(0u),expected);
        }
    }
    std::fesetround(saved_mode);
    // Independent, frozen exact triples generated with BigInt digit-pair
    // restoring arithmetic (and cross-checked with a binary-search oracle).
    // This block never consults GMP, even though older tests elsewhere do.
    std::ifstream golden(ARIADNE_BIG_UNSIGNED_SQRT_ORACLE_PATH);
    if (!golden) throw std::runtime_error("Missing BigUnsigned sqrt golden TSV");
    std::string line;
    std::size_t rows=0u;
    while (std::getline(golden,line)) {
        if (line.empty() || line.front()=='#') continue;
        std::istringstream fields(line);
        std::string input_hex,root_hex,remainder_hex,extra;
        if (!(fields>>input_hex>>root_hex>>remainder_hex) || (fields>>extra)) {
            throw std::runtime_error("Malformed BigUnsigned sqrt golden row");
        }
        const BigUnsigned input=big_unsigned_from_golden_hex(input_hex);
        const BigUnsigned original=input;
        const BigUnsigned root=isqrt(input);
        check_against_precomputed_hex(root,root_hex);
        const auto [root_with_remainder,remainder]=sqrtrem(input);
        check_against_precomputed_hex(root_with_remainder,root_hex);
        check_against_precomputed_hex(remainder,remainder_hex);
        // Algebraic properties supplement, never replace the frozen oracle.
        ARIADNE_TEST_ASSERT(compare(square(root),input)<=0);
        const BigUnsigned successor=add(root,BigUnsigned(std::uint64_t{1}));
        ARIADNE_TEST_ASSERT(compare(square(successor),input)>0);
        ARIADNE_TEST_ASSERT(compare(
            add(square(root_with_remainder),remainder),input)==0);
        ARIADNE_TEST_ASSERT(compare(remainder,
            add(root_with_remainder,root_with_remainder))<=0);
        ARIADNE_TEST_ASSERT(compare(input,original)==0);
        ++rows;
    }
    if (golden.bad() || rows!=166u) {
        throw std::runtime_error("Incomplete BigUnsigned sqrt golden TSV");
    }

    // Known roots exercise every even normalization shift, odd/even
    // limb counts, recursive splits and both workspace storage paths.
    // The expected answers below follow directly from (s-1)^2 < s^2-1
    // and s^2+2*s = (s+1)^2-1; no square-root oracle is involved.
    const BigUnsigned one(std::uint64_t{1}),two(std::uint64_t{2});
    for (std::size_t count:{2u,3u,4u,5u,7u,8u,9u,15u,16u,17u,
                           31u,32u,33u,64u,101u,102u,103u,128u}) {
        for (unsigned leading=0u;leading<64u;++leading) {
            std::vector<BigUnsigned::Limb> words(count);
            for (std::size_t i=0u;i<count;++i) {
                words[i]=0x9e3779b97f4a7c15ULL*(i+1u);
            }
            words.back()=std::numeric_limits<BigUnsigned::Limb>::max()>>leading;
            const BigUnsigned s(words),squared=square(s),twice=add(s,s);
            const auto verify=[&](BigUnsigned const& input,
                                  BigUnsigned const& expected_root,
                                  BigUnsigned const& expected_remainder) {
                const BigUnsigned original=input;
                ARIADNE_TEST_ASSERT(compare(isqrt(input),expected_root)==0);
                const auto [root,remainder]=sqrtrem(input);
                ARIADNE_TEST_ASSERT(compare(root,expected_root)==0);
                ARIADNE_TEST_ASSERT(compare(remainder,expected_remainder)==0);
                ARIADNE_TEST_ASSERT(compare(input,original)==0);
            };
            verify(squared,s,BigUnsigned{});
            verify(sub(squared,one),sub(s,one),sub(twice,two));
            verify(add(squared,one),s,one);
            verify(add(squared,twice),s,twice);
        }
    }
}

void TestBigUnsigned::test_power()
{
    const BigUnsigned zero,one(std::uint64_t{1});
    ARIADNE_TEST_ASSERT(compare(pow(zero,0u),one)==0);
    ARIADNE_TEST_ASSERT(pow(zero,std::numeric_limits<std::uint64_t>::max()).is_zero());
    ARIADNE_TEST_ASSERT(compare(pow(one,std::numeric_limits<std::uint64_t>::max()),one)==0);
    ARIADNE_TEST_ASSERT(compare(pow(BigUnsigned(std::uint64_t{3}),5u),
                              BigUnsigned(std::uint64_t{243}))==0);
    const std::array<std::uint64_t,2u> two64{0u,1u};
    const std::array<std::uint64_t,3u> two128{0u,0u,1u};
    ARIADNE_TEST_ASSERT(compare(pow(BigUnsigned(two64),2u),BigUnsigned(two128))==0);

    // Deterministic boundary-class sampling; skipped visits still consume RNG.
    // Exact source/result reference pairs were frozen independently
    // using native multiple-precision integer exponentiation.
    std::ifstream golden(ARIADNE_BIG_UNSIGNED_POWER_ORACLE_PATH);
    if (!golden) throw std::runtime_error("Missing BigUnsigned power golden TSV");
    std::string header;
    while (golden.peek()=='#') std::getline(golden,header);
    std::mt19937_64 generator(0xB7400A11ULL);
    for (const std::size_t count:{0u,1u,2u,3u,4u,5u,8u,9u,16u,17u}) {
        for (unsigned mode=0u;mode<12u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=generator();
            if (mode%4u==0u) std::fill(words.begin(),words.end(),0u);
            if (mode%4u==1u) std::fill(words.begin(),words.end(),
                std::numeric_limits<std::uint64_t>::max());
            if (mode%4u==2u && !words.empty()) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
            }
            const BigUnsigned input(words);
            for (const unsigned exponent:{0u,1u,2u,3u,4u,5u,7u,8u,9u,15u,16u,17u}) {
                if (mode!=0u && mode!=1u && mode!=2u && mode!=3u &&
                    mode!=7u) continue;
                if (exponent!=0u && exponent!=1u && exponent!=2u &&
                    exponent!=8u && exponent!=16u && exponent!=17u) continue;
                check_against_power_golden(
                    golden,'R',count,mode,exponent,input,pow(input,exponent));
            }
        }
    }
    // R15: unsigned powers now own both platform-specific binary
    // schedules. Probe irregular bit patterns and sign-independent
    // magnitude behavior against independent exact golden values,
    // with borrowed inputs unchanged.
    const std::array<std::uint64_t,3u> sample_limbs={
        std::uint64_t{0xfffffffffffffffdu},std::uint64_t{17u},
        std::uint64_t{1u}
    };
    const BigUnsigned sample_base(sample_limbs);
    for (const unsigned exponent:{0u,1u,2u,3u,5u,6u,7u,8u,9u,
                                 12u,15u,16u,17u,31u,32u,33u}) {
        check_against_power_golden(
            golden,'S',0u,0u,exponent,sample_base,pow(sample_base,exponent));
    }

    // Probe both sides of the conservative 1024-bit result bound.
    // Dense 65/341-bit bases also need the seventeenth temporary
    // product limb, even though the canonical result fits in sixteen.
    for (const unsigned bits:{1u,2u,53u,64u,65u,127u,128u,129u,
            255u,256u,257u,341u,342u,511u,512u,513u,1023u,1024u,1025u}) {
        // Build exactly 2^bits-1 directly from little-endian words,
        // preserving the original operand and all three exponent cases.
        std::vector<std::uint64_t> words((bits+63u)/64u,
                                         std::numeric_limits<std::uint64_t>::max());
        if (bits%64u != 0u) {
            words.back() = (std::uint64_t{1}<<(bits%64u))-1u;
        }
        const BigUnsigned input(words);
        const unsigned limit=1024u/bits;
        for (const unsigned exponent:{limit==0u ? 0u : limit-1u,limit,limit+1u}) {
            check_against_power_golden(
                golden,'B',bits,0u,exponent,input,pow(input,exponent));
        }
    }
    std::string trailing;
    if (golden >> trailing) {
        throw std::runtime_error("Extra BigUnsigned power golden rows");
    }
}

void TestBigUnsigned::test_short_owned_division()
{
    // Keep all 36 width pairs, sample modes by behavior (zero, maximal,
    // sparse, equal, borrow, dense); RNG consumes every original visit.
    // Exact source-keyed golden references were generated independently
    // of GMP and BigUnsigned, with the standardized mt19937_64 stream.
    std::ifstream golden(ARIADNE_BIG_UNSIGNED_SHORT_DIV_ORACLE_PATH);
    if (!golden) throw std::runtime_error("Missing BigUnsigned division golden TSV");
    std::string header;
    while (golden.peek()=='#') std::getline(golden,header);
    std::mt19937_64 generator(0xB74D1A1DEULL);
    for (std::size_t n=0u;n<=5u;++n) {
        for (std::size_t m=0u;m<=5u;++m) {
            for (unsigned mode=0u;mode<24u;++mode) {
                std::vector<std::uint64_t> left(n),right(m);
                for (auto& word:left) word=generator();
                for (auto& word:right) word=generator();
                if (mode%6u==0u) std::fill(left.begin(),left.end(),0u);
                if (mode%6u==1u) {
                    std::fill(left.begin(),left.end(),
                        std::numeric_limits<std::uint64_t>::max());
                    std::fill(right.begin(),right.end(),1u);
                }
                if (mode%6u==2u && !right.empty()) {
                    std::fill(right.begin(),right.end(),0u);
                    right.back()=std::uint64_t{1}<<(mode%64u);
                }
                if (mode%6u==3u && n==m) right=left;
                if (mode%6u==4u && n==m && !right.empty()) {
                    right=left;
                    right[0u]=0u;
                    left[0u]=std::numeric_limits<std::uint64_t>::max();
                }
                if (mode!=0u && mode!=1u && mode!=2u &&
                    mode!=3u && mode!=4u && mode!=5u &&
                    mode!=11u && mode!=17u && mode!=23u) continue;
                const BigUnsigned a(left),b(right);
                std::size_t ref_n=0u,ref_m=0u;
                unsigned ref_mode=0u;
                std::string ref_a,ref_b,ref_q,ref_r;
                if (!(golden >> ref_n >> ref_m >> ref_mode
                             >> ref_a >> ref_b >> ref_q >> ref_r)) {
                    throw std::runtime_error("Missing short-division golden row");
                }
                if (ref_n!=n || ref_m!=m || ref_mode!=mode) {
                    throw std::runtime_error("Incorrect short-division golden row order");
                }
                // Check the ORIGINAL operands as well as the resulting
                // quotient and remainder, rather than deriving lookup
                // keys from possibly incorrect arithmetic outputs.
                check_against_precomputed_hex(a,ref_a);
                check_against_precomputed_hex(b,ref_b);
                if (b.is_zero()) {
                    if (ref_q!="-" || ref_r!="-") {
                        throw std::runtime_error("Expected zero-divisor golden markers");
                    }
                    bool q_threw=false,r_threw=false;
                    try { (void)quot(a,b); } catch (std::domain_error const&) { q_threw=true; }
                    try { (void)rem(a,b); } catch (std::domain_error const&) { r_threw=true; }
                    ARIADNE_TEST_ASSERT(q_threw && r_threw);
                } else {
                    if (ref_q=="-" || ref_r=="-") {
                        throw std::runtime_error("Missing nonzero-division golden values");
                    }
                    check_against_precomputed_hex(quot(a,b),ref_q);
                    check_against_precomputed_hex(rem(a,b),ref_r);
                }
            }
        }
    }
    std::string extra;
    if (golden >> extra) {
        throw std::runtime_error("Extra short-division golden rows");
    }
    // Moved-from inline digits may still contain their previous bits.
    // Scalar addition must use the canonical size, not stale storage.
    BigUnsigned moved(std::numeric_limits<std::uint64_t>::max());
    const BigUnsigned owned(std::move(moved));
    const BigUnsigned one(std::uint64_t{1});
    ARIADNE_TEST_ASSERT(compare(add(moved,one),one)==0);
    ARIADNE_TEST_ASSERT(compare(add(one,moved),one)==0);
    const std::array<std::uint64_t,2u> carry{0u,1u};
    ARIADNE_TEST_ASSERT(compare(add(owned,one),BigUnsigned(carry))==0);
}

void TestBigUnsigned::test_storage()
{
    BigUnsigned zero;
    ARIADNE_TEST_ASSERT(zero.is_zero());
    ARIADNE_TEST_ASSERT(zero.limbs().empty());
    ARIADNE_TEST_ASSERT(zero.bit_length()==0u);

    BigUnsigned one(std::uint64_t{1});
    ARIADNE_TEST_ASSERT(one.limb(0)==1u);
    ARIADNE_TEST_ASSERT(one.limb(500)==0u);
    ARIADNE_TEST_ASSERT(one.bit_length()==1u);

    std::array<std::uint64_t,7> raw{1,2,3,4,5,0,0};
    BigUnsigned long_value(raw);
    ARIADNE_TEST_ASSERT(long_value.size()==5u);
    ARIADNE_TEST_ASSERT(long_value.limb(4)==5u);

    BigUnsigned copy(long_value);
    ARIADNE_TEST_ASSERT(compare(copy,long_value)==0);

    BigUnsigned moved(std::move(copy));
    ARIADNE_TEST_ASSERT(copy.is_zero());
    ARIADNE_TEST_ASSERT(compare(moved,long_value)==0);

    copy=one;
    copy=std::move(moved);
    ARIADNE_TEST_ASSERT(moved.is_zero());
    ARIADNE_TEST_ASSERT(compare(copy,long_value)==0);

    // Self-move must preserve both inline and heap-owned values.
    // Use an alias to exercise the public assignment without a compiler
    // diagnostic for a syntactically obvious self-move.
    auto* heap_alias=&copy;
    copy=std::move(*heap_alias);
    ARIADNE_TEST_ASSERT(compare(copy,long_value)==0);
    auto* inline_alias=&one;
    one=std::move(*inline_alias);
    ARIADNE_TEST_ASSERT(one.limb(0u)==1u && one.size()==1u);

    zero.swap(copy);
    ARIADNE_TEST_ASSERT(zero.size()==5u);
    ARIADNE_TEST_ASSERT(copy.is_zero());
    const std::array<std::uint64_t,5> zeros{};
    ARIADNE_TEST_ASSERT(BigUnsigned(zeros).is_zero());
}

void TestBigUnsigned::test_boundaries()
{
    const auto max=std::numeric_limits<std::uint64_t>::max();
    struct DiscardedRow { std::array<BigUnsigned::Limb,2u> words; std::size_t limit; bool expected; };
    for(auto const& row:{
        DiscardedRow{{0u,0u},0u,false},DiscardedRow{{1u,0u},0u,false},
        DiscardedRow{{1u,0u},1u,true},DiscardedRow{{2u,0u},1u,false},
        DiscardedRow{{2u,0u},2u,true},DiscardedRow{{0x8000000000000000ULL,0u},63u,false},
        DiscardedRow{{0x8000000000000000ULL,0u},64u,true},
        DiscardedRow{{0u,1u},64u,false},DiscardedRow{{0u,1u},65u,true},
        DiscardedRow{{0u,2u},65u,false},DiscardedRow{{0u,2u},66u,true},
        DiscardedRow{{0u,0u},128u,false},
        DiscardedRow{{0u,0u},std::numeric_limits<std::size_t>::max(),false},
        DiscardedRow{{0u,1u},std::numeric_limits<std::size_t>::max(),true}}) {
        ARIADNE_TEST_EQUALS(BigUnsigned::any_bits_below(row.words,row.limit),row.expected);
    }
    ARIADNE_TEST_ASSERT(!BigUnsigned::any_bits_below({},std::numeric_limits<std::size_t>::max()));
    const BigUnsigned all(max);
    const BigUnsigned one(std::uint64_t{1});

    const std::array<BigUnsigned::Limb,3u> noncanonical_shift{1u,0u,0u};
    std::array<BigUnsigned::Limb,4u> shift_destination{7u,7u,7u,7u};
    ARIADNE_TEST_EQUALS(BigUnsigned::shift_right_into(
        noncanonical_shift,0u,shift_destination),std::size_t{1u});
    ARIADNE_TEST_ASSERT((shift_destination==std::array<BigUnsigned::Limb,4u>{1u,0u,0u,0u}));
    ARIADNE_TEST_EQUALS(BigUnsigned::shift_right_into(
        noncanonical_shift,std::numeric_limits<std::size_t>::max(),{}),std::size_t{0u});
    ARIADNE_TEST_EQUALS(BigUnsigned::shift_right_into({},0u,{}),std::size_t{0u});
    std::array<BigUnsigned::Limb,1u> short_shift_output{7u};
    ARIADNE_TEST_FAIL(BigUnsigned::shift_right_into(noncanonical_shift,0u,short_shift_output));
    ARIADNE_TEST_EQUALS(short_shift_output[0u],BigUnsigned::Limb{7u});
    ARIADNE_TEST_ASSERT(add(all,one).size()==2u);
    ARIADNE_TEST_ASSERT(add(all,one).limb(0)==0u);
    ARIADNE_TEST_ASSERT(add(all,one).limb(1)==1u);
    ARIADNE_TEST_ASSERT(compare(sub(all,all),BigUnsigned{})==0);
    ARIADNE_TEST_ASSERT(compare(shift_right(shift_left(all,257u),257u),all)==0);
    ARIADNE_TEST_ASSERT(compare(shift_right(all,64u),BigUnsigned{})==0);
    ARIADNE_TEST_ASSERT(compare(shift_left(one,256u),BigUnsigned(
        std::array<std::uint64_t,5>{0,0,0,0,1}))==0);
    ARIADNE_TEST_ASSERT(compare(shift_right(shift_left(one,63u),63u),one)==0);

    bool threw = false;
    try {
        (void)sub(one,all);
    } catch (std::underflow_error const&) {
        threw = true;
    }
    ARIADNE_TEST_ASSERT(threw);
    ARIADNE_TEST_ASSERT(compare(add(all,BigUnsigned{}),all)==0);

    // MSVC checked std::span iterators require both algorithm endpoints to
    // originate from the same view. Exercise unequal-width add/sub paths,
    // which reach add_into/sub_into instead of the one-limb fast paths.
    const std::array<std::uint64_t,2u> two64_words{0u,1u};
    const BigUnsigned two64(two64_words);
    const std::array<std::uint64_t,2u> two64_plus_one_words{1u,1u};
    ARIADNE_TEST_ASSERT(compare(
        add(two64,one),BigUnsigned(two64_plus_one_words))==0);
    ARIADNE_TEST_ASSERT(compare(sub(two64,one),all)==0);

    // Scan across equal high limbs in heap storage. The first
    // differing limb is near the bottom, or does not exist at all.
    std::array<std::uint64_t,64> equal_high{};
    equal_high.fill(max);
    std::array<std::uint64_t,64> low_different=equal_high;
    low_different[0u]=max-1u;
    const BigUnsigned high(equal_high);
    const BigUnsigned low(low_different);
    ARIADNE_TEST_ASSERT(compare(low,high)<0);
    ARIADNE_TEST_ASSERT(compare(high,low)>0);
    ARIADNE_TEST_ASSERT(compare(high,high)==0);
    ARIADNE_TEST_ASSERT(compare(high,one)>0);
    low_different[63u]=max-1u;
    ARIADNE_TEST_ASSERT(compare(BigUnsigned(low_different),high)<0);
}

// B41: the highest differing limb determines both unsigned order
// and the maximum required result width. Exercise common-prefix
// elimination against the independent B133 oracle in both directions,
// including underflow, canonical zero, dense inputs and long borrows.
void TestBigUnsigned::test_subtraction_common_prefix()
{
    // One source-keyed representative of each of the nine original
    // structural classes at every nonzero width. Include the zero width.
    // No GMP or test-subject calculation is used to derive expectations.
    std::ifstream golden(ARIADNE_BIG_UNSIGNED_SUB_ORACLE_PATH);
    if (!golden) throw std::runtime_error("Missing BigUnsigned subtraction TSV");
    std::string header;
    while (golden.peek()=='#') std::getline(golden,header);
    for (const std::size_t width:
            {0u,1u,2u,7u,8u,16u,31u,32u,63u,64u,65u,128u}) {
        const unsigned classes=width==0u ? 1u : 9u;
        for (unsigned mode=0u;mode<classes;++mode) {
            std::size_t row_width=0u;
            unsigned row_mode=0u;
            std::string lhs_hex,rhs_hex,forward_hex,reverse_hex;
            if (!(golden >> row_width >> row_mode >> lhs_hex >> rhs_hex
                         >> forward_hex >> reverse_hex)
                || row_width!=width || row_mode!=mode) {
                throw std::runtime_error("Missing or out-of-order subtraction TSV row");
            }
            const BigUnsigned lhs=big_unsigned_from_golden_hex(lhs_hex);
            const BigUnsigned rhs=big_unsigned_from_golden_hex(rhs_hex);
            check_against_precomputed_hex(lhs,lhs_hex);
            check_against_precomputed_hex(rhs,rhs_hex);
            const auto verify=[&](BigUnsigned const& a,BigUnsigned const& b,
                                  std::string const& expected_hex) {
                if (expected_hex=="-") {
                    bool threw=false;
                    try { (void)sub(a,b); }
                    catch (std::underflow_error const&) { threw=true; }
                    ARIADNE_TEST_ASSERT(threw);
                } else {
                    check_against_precomputed_hex(sub(a,b),expected_hex);
                }
            };
            verify(lhs,rhs,forward_hex);
            verify(rhs,lhs,reverse_hex);
        }
    }
    std::string extra;
    if (golden >> extra) throw std::runtime_error("Extra subtraction TSV rows");
}


void TestBigUnsigned::test_multiplication()
{
    const auto max = std::numeric_limits<std::uint64_t>::max();
    const BigUnsigned one(std::uint64_t{1});
    const std::array<std::uint64_t,4> four_max{max,max,max,max};
    const BigUnsigned four(four_max);

    ARIADNE_TEST_ASSERT(mul_low(four,one,0u).is_zero());
    ARIADNE_TEST_ASSERT(compare(mul_low(four,one,1u),BigUnsigned(max)) == 0);
    ARIADNE_TEST_ASSERT(compare(mul(four,one),four) == 0);
    ARIADNE_TEST_ASSERT(compare(square(one),one) == 0);

    std::mt19937_64 engine(0x1282561024ULL);
    // R14: independent frozen self-product references. Preserve all
    // original draws from the shared engine: later randomized tests
    // must continue to use exactly the same inputs.
    std::ifstream square_golden(ARIADNE_BIG_UNSIGNED_R14_ORACLE_PATH);
    if (!square_golden) throw std::runtime_error("Missing R14 square TSV");
    std::string square_header;
    while (square_golden.peek()=='#') std::getline(square_golden,square_header);
    for (const std::size_t count:{1u,4u,7u,8u,9u,16u,64u}) {
        for (unsigned mode=0u;mode<12u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=engine();
            if (mode%4u==0u) {
                std::fill(words.begin(),words.end(),0u);
            } else if (mode%4u==1u) {
                std::fill(words.begin(),words.end(),
                          std::numeric_limits<std::uint64_t>::max());
            } else if (mode%4u==2u) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
            } else {
                words.back()|=std::uint64_t{1}<<63u;
            }
            if (mode>=4u) continue; // Four structurally distinct classes.
            std::size_t ref_count=0u;
            unsigned ref_mode=0u;
            std::string ref_input,ref_square;
            if (!(square_golden >> ref_count >> ref_mode
                               >> ref_input >> ref_square)
                || ref_count!=count || ref_mode!=mode) {
                throw std::runtime_error("Missing or unordered R14 square row");
            }
            const BigUnsigned input(words),distinct(input);
            check_against_precomputed_hex(input,ref_input);
            check_against_precomputed_hex(distinct,ref_input);
            check_against_precomputed_hex(mul(input,input),ref_square);
            check_against_precomputed_hex(mul(input,distinct),ref_square);
            check_against_precomputed_hex(square(input),ref_square);
        }
    }
    std::string extra_square;
    if (square_golden >> extra_square)
        throw std::runtime_error("Extra R14 square oracle rows");

    // B72: independent frozen product and square expectations for
    // every width pair on either side of the owning 4/5-limb threshold.
    // All original random draws are consumed before selecting a class.
    std::ifstream product_golden(ARIADNE_BIG_UNSIGNED_B72_ORACLE_PATH);
    if (!product_golden) throw std::runtime_error("Missing B72 product TSV");
    std::string product_header;
    while (product_golden.peek()=='#')
        std::getline(product_golden,product_header);
    std::mt19937_64 short_product_engine(0xB7204455ULL);
    for (const std::size_t n:{0u,1u,2u,3u,4u,5u}) {
        for (const std::size_t m:{0u,1u,2u,3u,4u,5u}) {
            for (unsigned mode=0u;mode<32u;++mode) {
                std::array<std::uint64_t,5u> left{},right{};
                for (std::size_t j=0u;j<n;++j)
                    left[j]=short_product_engine();
                for (std::size_t j=0u;j<m;++j)
                    right[j]=short_product_engine();
                if (mode%8u==0u) {
                    std::fill(left.begin(),left.end(),max);
                    std::fill(right.begin(),right.end(),max);
                } else if (mode%8u==1u) {
                    std::fill(left.begin(),left.end(),0u);
                    std::fill(right.begin(),right.end(),0u);
                } else if (mode%8u==2u) {
                    std::fill(left.begin(),left.end(),0u);
                    std::fill(right.begin(),right.end(),0u);
                    if (n!=0u) left[n-1u]=1u;
                    if (m!=0u) right[m-1u]=1u;
                } else if (mode%8u==3u) {
                    for (std::size_t j=0u;j<n;j+=2u) left[j]=0u;
                    for (std::size_t j=0u;j<m;j+=2u) right[j]=0u;
                } else if (mode%8u==4u) {
                    for (std::size_t j=1u;j<n;++j) left[j]=0u;
                    for (std::size_t j=1u;j<m;++j) right[j]=0u;
                }
                if (mode>=8u) continue; // All eight behavior classes.
                std::size_t ref_n=0u,ref_m=0u;
                unsigned ref_mode=0u;
                std::string ref_lhs,ref_rhs,ref_product,ref_square;
                if (!(product_golden >> ref_n >> ref_m >> ref_mode >>
                                      ref_lhs >> ref_rhs >>
                                      ref_product >> ref_square)
                    || ref_n!=n || ref_m!=m || ref_mode!=mode) {
                    throw std::runtime_error("Missing or unordered B72 product row");
                }
                const BigUnsigned lhs(BigUnsigned::LimbView(left.data(),n));
                const BigUnsigned rhs(BigUnsigned::LimbView(right.data(),m));
                check_against_precomputed_hex(lhs,ref_lhs);
                check_against_precomputed_hex(rhs,ref_rhs);
                check_against_precomputed_hex(mul(lhs,rhs),ref_product);
                check_against_precomputed_hex(mul(rhs,lhs),ref_product);
                check_against_precomputed_hex(mul(lhs,lhs),ref_square);
                check_against_precomputed_hex(square(lhs),ref_square);
            }
        }
    }
    std::string extra_product;
    if (product_golden >> extra_product)
        throw std::runtime_error("Extra B72 product oracle rows");

    // Independently frozen original-source random multiplication corpus.
    // One representative per ten source visits, stratified over all 34
    // operand widths, truncation lengths, dense, sparse and zero classes.
    std::ifstream random_golden(ARIADNE_BIG_UNSIGNED_RANDOM_MUL_ORACLE_PATH);
    if (!random_golden)
        throw std::runtime_error("Missing BigUnsigned random multiplication TSV");
    std::string random_header;
    while (random_golden.peek()=='#')
        std::getline(random_golden,random_header);
    unsigned previous_repetition=0u;
    for (unsigned sample=0u;sample<600u;++sample) {
        unsigned repetition=0u;
        std::size_t original_n=0u,original_m=0u,low_limbs=0u;
        std::string source_x,source_y,product,squared,truncated;
        if (!(random_golden >> repetition >> original_n >> original_m
               >> low_limbs >> source_x >> source_y >> product
               >> squared >> truncated) ||
            repetition/10u!=sample ||
            (sample!=0u && repetition<=previous_repetition) ||
            original_n>33u || original_m>33u || low_limbs>79u) {
            throw std::runtime_error(
                "Missing or out-of-order random multiplication TSV row");
        }
        previous_repetition=repetition;
        const BigUnsigned x=big_unsigned_from_golden_hex(source_x);
        const BigUnsigned y=big_unsigned_from_golden_hex(source_y);
        ARIADNE_TEST_ASSERT(x.size()<=original_n);
        ARIADNE_TEST_ASSERT(y.size()<=original_m);
        check_against_precomputed_hex(mul(x,y),product);
        check_against_precomputed_hex(mul(y,x),product);
        check_against_precomputed_hex(square(x),squared);
        check_against_precomputed_hex(mul(x,x),squared);
        check_against_precomputed_hex(mul_low(x,y,low_limbs),truncated);
    }
    std::string random_extra;
    if (random_golden >> random_extra)
        throw std::runtime_error("Extra random multiplication TSV rows");

    // B35: freeze five behaviorally distinct modes for every original
    // Karatsuba width pair (including 63/64/65 and unbalanced lengths).
    // Both operand orders compare against one independent exact reference.
    std::ifstream b35_golden(ARIADNE_BIG_UNSIGNED_B35_ORACLE_PATH);
    if (!b35_golden) throw std::runtime_error("Missing BigUnsigned B35 oracle TSV");
    std::string b35_header;
    while (b35_golden.peek()=='#') std::getline(b35_golden,b35_header);
    for (const auto& [n,m]:
            {std::pair<std::size_t,std::size_t>{63u,63u},
             {63u,64u}, {64u,63u}, {64u,64u}, {64u,65u},
             {65u,64u}, {65u,65u}, {80u,64u}, {96u,64u},
             {96u,96u}, {128u,64u}, {128u,96u},
             {96u,128u}, {128u,128u}}) {
        for (unsigned repetition=0u;repetition<5u;++repetition) {
            std::size_t source_n=0u,source_m=0u;
            unsigned source_repetition=0u;
            std::string lhs_hex,rhs_hex,product_hex;
            if (!(b35_golden >> source_n >> source_m >> source_repetition
                        >> lhs_hex >> rhs_hex >> product_hex)
                || source_n!=n || source_m!=m ||
                   source_repetition!=repetition) {
                throw std::runtime_error("Missing or unordered B35 TSV row");
            }
            const BigUnsigned lhs=big_unsigned_from_golden_hex(lhs_hex);
            const BigUnsigned rhs=big_unsigned_from_golden_hex(rhs_hex);
            ARIADNE_TEST_ASSERT(lhs.size()==n);
            ARIADNE_TEST_ASSERT(rhs.size()==m);
            const BigUnsigned expected_product=
                big_unsigned_from_golden_hex(product_hex);
            check_against_precomputed_hex(mul(lhs,rhs),product_hex);
            ARIADNE_TEST_ASSERT(compare(mul(rhs,lhs),expected_product)==0);
        }
    }
    std::string b35_extra;
    if (b35_golden >> b35_extra)
        throw std::runtime_error("Extra B35 TSV oracle rows");

    // B40: all eight zero-row density classes and the all-ones
    // saturation case on every original Karatsuba/dispatch width pair.
    // Independent source-keyed oracle, including truncated products.
    std::ifstream b40_golden(ARIADNE_BIG_UNSIGNED_B40_ORACLE_PATH);
    if (!b40_golden) throw std::runtime_error("Missing B40 TSV");
    std::string b40_header;
    while (b40_golden.peek()=='#') std::getline(b40_golden,b40_header);
    for (const auto& [n,m]:
            {std::pair<std::size_t,std::size_t>{7u,64u},
             {8u,64u}, {16u,16u}, {31u,32u}, {32u,31u},
             {32u,32u}, {32u,33u}, {33u,32u}, {33u,33u},
             {63u,63u}, {63u,64u}, {64u,63u},
             {64u,64u}, {64u,65u}, {65u,64u},
             {64u,128u}, {128u,64u}, {128u,128u}}) {
        for (unsigned ordinal=0u;ordinal<9u;++ordinal) {
            const unsigned repetition=ordinal==8u?15u:ordinal;
            std::size_t source_n=0u,source_m=0u;
            unsigned source_rep=0u;
            std::string lhs_hex,rhs_hex,prod_hex,low_hex,square_hex;
            if (!(b40_golden >> source_n >> source_m >> source_rep
                    >> lhs_hex >> rhs_hex >> prod_hex >> low_hex >> square_hex)
                || source_n!=n || source_m!=m || source_rep!=repetition) {
                throw std::runtime_error("Missing or unordered B40 TSV row");
            }
            const BigUnsigned lhs=big_unsigned_from_golden_hex(lhs_hex);
            const BigUnsigned rhs=big_unsigned_from_golden_hex(rhs_hex);
            ARIADNE_TEST_ASSERT(lhs.size()==n && rhs.size()==m);
            check_against_precomputed_hex(mul(lhs,rhs),prod_hex);
            check_against_precomputed_hex(mul(rhs,lhs),prod_hex);
            const std::size_t low=repetition%2u==0u?1u:n+m-1u;
            check_against_precomputed_hex(mul_low(lhs,rhs,low),low_hex);
            check_against_precomputed_hex(mul_low(rhs,lhs,low),low_hex);
            check_against_precomputed_hex(square(lhs),square_hex);
            if (n>=64u && m>=64u) {
                const std::size_t required=
                    BigUnsigned::karatsuba_scratch_limbs(n,m,32u);
                std::vector<std::uint64_t> scratch(required),output(n+m);
                BigUnsigned::mul_into_karatsuba(
                    lhs.limbs(),rhs.limbs(),output,scratch,32u);
                check_against_precomputed_hex(BigUnsigned(output),prod_hex);
            }
        }
    }
    std::string b40_extra;
    if (b40_golden >> b40_extra)
        throw std::runtime_error("Extra B40 TSV rows");

    // B39: 19 original dispatch-width pairs, all seven density modes
    // on BOTH operands, with both low-product truncation lengths.
    // Original source-keyed values and independent exact results in TSV.
    std::ifstream b39_golden(ARIADNE_BIG_UNSIGNED_B39_ORACLE_PATH);
    if (!b39_golden) throw std::runtime_error("Missing B39 TSV");
    std::string b39_header;
    while (b39_golden.peek()=='#') std::getline(b39_golden,b39_header);
    for (const auto& [n,m]:
            {std::pair<std::size_t,std::size_t>{7u,8u},
             {8u,8u}, {15u,16u}, {16u,16u},
             {32u,32u}, {63u,63u}, {63u,64u},
             {64u,63u}, {64u,64u}, {64u,65u},
             {65u,64u}, {64u,85u}, {64u,86u},
             {85u,64u}, {86u,64u}, {96u,64u},
             {128u,64u}, {64u,128u}, {128u,128u}}) {
        for (const unsigned repetition:{0u,1u,2u,3u,4u,5u,6u,9u}) {
            std::size_t ref_n=0u,ref_m=0u;
            unsigned ref_repetition=0u;
            std::string lhs_hex,rhs_hex,product_hex,square_hex,low_hex;
            if (!(b39_golden >> ref_n >> ref_m >> ref_repetition
                        >> lhs_hex >> rhs_hex >> product_hex
                        >> square_hex >> low_hex)
                || ref_n!=n || ref_m!=m || ref_repetition!=repetition) {
                throw std::runtime_error("Missing or unordered B39 TSV row");
            }
            const BigUnsigned lhs=big_unsigned_from_golden_hex(lhs_hex);
            const BigUnsigned rhs=big_unsigned_from_golden_hex(rhs_hex);
            ARIADNE_TEST_ASSERT(lhs.size()==n && rhs.size()==m);
            const BigUnsigned product=mul(lhs,rhs);
            check_against_precomputed_hex(product,product_hex);
            ARIADNE_TEST_ASSERT(compare(mul(rhs,lhs),product)==0);
            check_against_precomputed_hex(square(lhs),square_hex);
            const std::size_t low=repetition%2u==0u
                ? (n+m)/2u : n+m-1u;
            check_against_precomputed_hex(mul_low(lhs,rhs,low),low_hex);
        }
    }
    std::string b39_extra;
    if (b39_golden >> b39_extra)
        throw std::runtime_error("Extra B39 TSV rows");

    // B54: all 17 width pairs and 8 distinct density patterns.
    std::ifstream b54_golden(ARIADNE_BIG_UNSIGNED_B54_PRODUCT_ORACLE_PATH);
    if (!b54_golden) throw std::runtime_error("Missing B54 TSV");
    std::string b54_header;
    while (b54_golden.peek()=='#') std::getline(b54_golden,b54_header);
    for (const auto& [n,m]:std::initializer_list<std::pair<std::size_t,std::size_t>>{{32u,32u},{33u,33u},{63u,63u},{63u,64u},{64u,63u},{64u,64u},{64u,65u},{65u,64u},{65u,65u},{64u,85u},{64u,86u},{85u,64u},{86u,64u},{96u,96u},{127u,128u},{128u,127u},{128u,128u}}) {
        for (unsigned mode=0u;mode<8u;++mode) {
            std::size_t ref_n=0u,ref_m=0u;
            unsigned ref_mode=0u;
            std::string lhs_hex,rhs_hex,result_hex;
            if (!(b54_golden >> ref_n >> ref_m >> ref_mode
                  >> lhs_hex >> rhs_hex >> result_hex)
                || ref_n!=n || ref_m!=m || ref_mode!=mode)
                throw std::runtime_error("Missing or unordered B54 TSV row");
            const BigUnsigned lhs=big_unsigned_from_golden_hex(lhs_hex);
            const BigUnsigned rhs=big_unsigned_from_golden_hex(rhs_hex);
            ARIADNE_TEST_ASSERT(lhs.size()==n && rhs.size()==m);
            check_against_precomputed_hex(mul(lhs,rhs),result_hex);
            check_against_precomputed_hex(mul(rhs,lhs),result_hex);
            std::vector<std::uint64_t> fused(n+m);
            BigUnsigned::mul_into_fused(lhs.limbs(),rhs.limbs(),fused);
            check_against_precomputed_hex(BigUnsigned(fused),result_hex);
        }
    }
    std::string b54_extra;
    if (b54_golden >> b54_extra) throw std::runtime_error("Extra B54 TSV rows");

    // B54 squares: all original widths and eight distinct mode%8
    // structures, with the original stream frozen independently.
    std::ifstream b54_squares(ARIADNE_BIG_UNSIGNED_B54_SQUARE_ORACLE_PATH);
    if (!b54_squares) throw std::runtime_error("Missing B54 squares TSV");
    std::string b54_square_header;
    while (b54_squares.peek()=='#')
        std::getline(b54_squares,b54_square_header);
    for (const std::size_t n:
            {7u,8u,9u,16u,31u,32u,33u,
             63u,64u,65u,66u,96u,128u}) {
        for (unsigned mode=0u;mode<8u;++mode) {
            std::size_t ref_n=0u;
            unsigned ref_mode=0u;
            std::string input_hex,result_hex;
            if (!(b54_squares >> ref_n >> ref_mode >> input_hex >> result_hex)
                || ref_n!=n || ref_mode!=mode)
                throw std::runtime_error("Missing or unordered B54 square TSV row");
            const BigUnsigned input=big_unsigned_from_golden_hex(input_hex);
            ARIADNE_TEST_ASSERT(input.size()==n);
            check_against_precomputed_hex(square(input),result_hex);
            check_against_precomputed_hex(mul(input,input),result_hex);
        }
    }
    std::string b54_square_extra;
    if (b54_squares >> b54_square_extra)
        throw std::runtime_error("Extra B54 square TSV rows");

    // B37: all seven structural square classes at every original width.
    std::ifstream b37_golden(ARIADNE_BIG_UNSIGNED_B37_ORACLE_PATH);
    if (!b37_golden) throw std::runtime_error("Missing B37 oracle TSV");
    std::string b37_header;
    while (b37_golden.peek()=='#') std::getline(b37_golden,b37_header);
    for (const std::size_t count:
            {1u,4u,7u,8u,9u,15u,16u,17u,32u,33u,64u,65u,128u}) {
        for (unsigned repetition=0u;repetition<7u;++repetition) {
            std::size_t ref_count=0u;
            unsigned ref_repetition=0u;
            std::string raw_input,raw_square;
            if (!(b37_golden >> ref_count >> ref_repetition
                  >> raw_input >> raw_square)
                || ref_count!=count || ref_repetition!=repetition)
                throw std::runtime_error("Missing or unordered B37 oracle row");
            const BigUnsigned input=big_unsigned_from_golden_hex(raw_input);
            ARIADNE_TEST_ASSERT(input.size()<=count);
            check_against_precomputed_hex(square(input),raw_square);
            check_against_precomputed_hex(mul(input,input),raw_square);
        }
    }
    std::string b37_extra;
    if (b37_golden >> b37_extra)
        throw std::runtime_error("Extra B37 oracle rows");

    // B38: retain each sparse/dense mode at every asymmetry/width.
    // Full and six truncated products use independent exact TSV values.
    std::ifstream b38_golden(ARIADNE_BIG_UNSIGNED_B38_ORACLE_PATH);
    if (!b38_golden) throw std::runtime_error("Missing B38 oracle TSV");
    std::string b38_header;
    while (b38_golden.peek()=='#') std::getline(b38_golden,b38_header);
    for (const auto& [n,m]:
            std::initializer_list<std::pair<std::size_t,std::size_t>>{{1u,3u},{4u,7u},{7u,8u},{8u,8u},{9u,16u},{16u,32u},{32u,64u},{64u,32u},{64u,128u},{128u,64u},{64u,129u},{129u,64u},{96u,128u},{128u,96u}}) {
        for (unsigned repetition=0u;repetition<8u;++repetition) {
            std::size_t ref_n=0u,ref_m=0u;
            unsigned ref_rep=0u;
            std::string lhs_hex,rhs_hex,product_hex;
            std::array<std::string,6u> low_hex;
            if (!(b38_golden >> ref_n >> ref_m >> ref_rep
                 >> lhs_hex >> rhs_hex >> product_hex
                 >> low_hex[0] >> low_hex[1] >> low_hex[2]
                 >> low_hex[3] >> low_hex[4] >> low_hex[5])
                || ref_n!=n || ref_m!=m || ref_rep!=repetition)
                throw std::runtime_error("Missing or unordered B38 oracle row");
            const BigUnsigned lhs=big_unsigned_from_golden_hex(lhs_hex);
            const BigUnsigned rhs=big_unsigned_from_golden_hex(rhs_hex);
            ARIADNE_TEST_ASSERT(lhs.size()==n && rhs.size()==m);
            check_against_precomputed_hex(mul(lhs,rhs),product_hex);
            check_against_precomputed_hex(mul(rhs,lhs),product_hex);
            const std::array<std::size_t,6u> low_lengths{
                0u,1u,n/2u,n,n+m-1u,n+m+2u};
            for (std::size_t i=0u;i<low_lengths.size();++i) {
                check_against_precomputed_hex(
                    mul_low(lhs,rhs,low_lengths[i]),low_hex[i]);
                check_against_precomputed_hex(
                    mul_low(rhs,lhs,low_lengths[i]),low_hex[i]);
            }
        }
    }
    std::string b38_extra;
    if (b38_golden >> b38_extra)
        throw std::runtime_error("Extra B38 oracle rows");

    // B113/B32: retain all 2,112 original random-source visits, including
    // overwritten mt19937_64 draws. Exact expectations are independent TSV
    // squares, not GMP computations at runtime.
    std::ifstream b32_golden(ARIADNE_BIG_UNSIGNED_B32_ORACLE_PATH);
    if (!b32_golden) throw std::runtime_error("Missing B32 square TSV");
    std::string b32_header;
    while (b32_golden.peek()=='#') std::getline(b32_golden,b32_header);
    std::mt19937_64 square_engine(0xB3205A7EULL);
    for (const std::size_t count:
            {1u,2u,7u,8u,9u,12u,16u,32u,64u,96u,128u}) {
        for (unsigned repetition=0u;repetition<192u;++repetition) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=square_engine();
            if (repetition%13u==0u) {
                std::fill(words.begin(),words.end(),max);
            } else if (repetition%17u==0u) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
            } else if (repetition%19u==0u) {
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()=max;
            } else if (repetition%23u==0u) {
                std::fill(words.begin(),words.end(),0u);
            }
            std::size_t source_count=0u;
            unsigned source_repetition=0u;
            std::string input_hex,square_hex;
            if (!(b32_golden >> source_count >> source_repetition
                             >> input_hex >> square_hex)
                || source_count!=count || source_repetition!=repetition)
                throw std::runtime_error("Missing or unordered B32 TSV row");
            const BigUnsigned value(words);
            ARIADNE_TEST_ASSERT(value.size()<=count);
            check_against_precomputed_hex(value,input_hex);
            check_against_precomputed_hex(square(value),square_hex);
            check_against_precomputed_hex(mul(value,value),square_hex);
            check_against_precomputed_hex(value,input_hex);
        }
    }
    std::string b32_extra;
    if (b32_golden >> b32_extra)
        throw std::runtime_error("Extra B32 TSV rows");
}

void TestBigUnsigned::test_views_and_scratch()
{
    const auto max = std::numeric_limits<std::uint64_t>::max();
    BigUnsigned::Scratch scratch;
    ARIADNE_TEST_ASSERT(scratch.capacity() >= 16u);

    auto short_view = scratch.acquire(4u);
    for (auto limb : short_view) {
        ARIADNE_TEST_ASSERT(limb == 0u);
    }
    short_view[0] = 11u;
    BigUnsigned x;
    x.assign(short_view);
    ARIADNE_TEST_ASSERT(x.size() == 1u);
    ARIADNE_TEST_ASSERT(x.limb(0) == 11u);

    auto heap_view = scratch.acquire(35u);
    ARIADNE_TEST_ASSERT(scratch.capacity() >= 35u);
    for (std::size_t i=0u; i<heap_view.size(); ++i) {
        heap_view[i] = max;
    }
    x.assign(heap_view);
    ARIADNE_TEST_ASSERT(x.size() == 35u);
    const auto* old_heap = heap_view.data();
    heap_view = scratch.acquire(25u);
    ARIADNE_TEST_ASSERT(heap_view.data() == old_heap);
    for (auto limb : heap_view) {
        ARIADNE_TEST_ASSERT(limb == 0u);
    }

    x.assign(x.limbs());
    ARIADNE_TEST_ASSERT(x.size() == 35u);
    x.assign(x.limbs().subspan(3u,2u));
    ARIADNE_TEST_ASSERT(x.size() == 2u);
    x.assign(x.limbs().first(1u));
    ARIADNE_TEST_ASSERT(x.size() == 1u);

    BigUnsigned::Scratch output;
    std::mt19937_64 engine(0xB22B22ULL);
    std::array<std::ifstream,5u> b134_golden;
    for (std::size_t part=0u;part<b134_golden.size();++part) {
        const std::string path=std::string(ARIADNE_TEST_DATA_DIR)+
            "/big_unsigned_b134_views_"+std::to_string(part)+".tsv";
        b134_golden[part].open(path);
        if (!b134_golden[part])
            throw std::runtime_error("Missing B134 views TSV: "+path);
        std::string header;
        while (b134_golden[part].peek()=='#')
            std::getline(b134_golden[part],header);
    }
    for (unsigned repetition=0u; repetition<5000u; ++repetition) {
        const std::size_t n = engine()%28u;
        const std::size_t m = engine()%28u;
        std::array<std::uint64_t,28> da{},db{};
        for (std::size_t i=0u; i<n; ++i) da[i] = engine();
        for (std::size_t i=0u; i<m; ++i) db[i] = engine();
        if (repetition%7u == 0u) {
            for (std::size_t i=0u; i<n; ++i) da[i] = max;
        }
        if (repetition%11u == 0u) {
            for (std::size_t i=0u; i<m; ++i) db[i] = max;
        }
        const BigUnsigned lhs(std::span(da.data(),n));
        const BigUnsigned rhs(std::span(db.data(),m));

        const std::size_t low = repetition%5u==0u ? 0u :
            repetition%5u==1u ? 1u : repetition%5u==2u ? 4u :
            repetition%5u==3u ? 15u : engine()%65u;
        unsigned source_repetition=0u;
        std::size_t source_n=0u,source_m=0u,source_low=0u;
        std::string lhs_hex,rhs_hex,sum_hex,difference_hex,product_hex;
        std::ifstream& golden=b134_golden[repetition/1000u];
        if (!(golden >> source_repetition >> source_n >> source_m
              >> source_low >> lhs_hex >> rhs_hex >> sum_hex
              >> difference_hex >> product_hex) ||
            source_repetition!=repetition || source_n!=n ||
            source_m!=m || source_low!=low)
            throw std::runtime_error(
                "Missing or unordered B134 views TSV row");
        check_against_precomputed_hex(lhs,lhs_hex);
        check_against_precomputed_hex(rhs,rhs_hex);

        const std::size_t count = std::max(lhs.size(),rhs.size());
        auto sum_view = output.acquire(count+2u);
        const auto carry = BigUnsigned::add_into(
            lhs.limbs(),rhs.limbs(),sum_view.first(count));
        sum_view[count] = carry;
        ARIADNE_TEST_ASSERT(sum_view[count+1u] == 0u);
        check_against_precomputed_hex(BigUnsigned(sum_view),sum_hex);

        auto difference_view = output.acquire(count+1u);
        const auto borrow = BigUnsigned::sub_into(
            lhs.limbs(),rhs.limbs(),difference_view);
        ARIADNE_TEST_ASSERT(difference_view[count] == 0u);
        if (difference_hex!="-") {
            ARIADNE_TEST_ASSERT(compare(lhs,rhs)>=0);
            ARIADNE_TEST_ASSERT(borrow == 0u);
            check_against_precomputed_hex(
                BigUnsigned(difference_view),difference_hex);
        } else {
            ARIADNE_TEST_ASSERT(compare(lhs,rhs)<0);
            if (count != 0u) ARIADNE_TEST_ASSERT(borrow == 1u);
        }

        auto product_view = output.acquire(low);
        BigUnsigned::mul_low_into(lhs.limbs(),rhs.limbs(),product_view);
        check_against_precomputed_hex(BigUnsigned(product_view),product_hex);
        check_against_precomputed_hex(lhs,lhs_hex);
        check_against_precomputed_hex(rhs,rhs_hex);
    }
    for (auto& golden:b134_golden) {
        std::string extra;
        if (golden >> extra)
            throw std::runtime_error("Extra B134 views TSV rows");
    }
    // B135/R38: preserve all shifted add/sub source cases and shift
    // combinations using independent source-keyed TSV expectations.
    std::ifstream b135_sources(
        std::string(ARIADNE_TEST_DATA_DIR)+
        "/big_unsigned_b135_shifted_sources.tsv");
    if (!b135_sources)
        throw std::runtime_error("Missing B135 shifted source TSV");
    std::string b135_header;
    while (b135_sources.peek()=='#') std::getline(b135_sources,b135_header);
    std::array<std::ifstream,10u> b135_golden;
    for (std::size_t part=0u;part<b135_golden.size();++part) {
        const std::string path=std::string(ARIADNE_TEST_DATA_DIR)+
            "/big_unsigned_b135_shifted_"+std::to_string(part)+".tsv";
        b135_golden[part].open(path);
        if (!b135_golden[part])
            throw std::runtime_error("Missing B135 shifted TSV: "+path);
        std::string header;
        while (b135_golden[part].peek()=='#')
            std::getline(b135_golden[part],header);
    }
    const auto b135_part=[](std::size_t ordinal) {
        return ordinal<5775u ? ordinal/825u :
            7u+(ordinal-5775u)/275u;
    };
    std::mt19937_64 shifted_engine(0xB38038ULL);
    std::size_t source_ordinal=0u,ordinal=0u;
    for (std::size_t width:{0u,1u,4u,16u,17u,31u,32u,33u,64u,65u,129u}) {
        for (unsigned case_index=0u;case_index<24u;++case_index) {
            const std::size_t length_b=case_index%5u==0u ? width :
                case_index%5u==1u ? 0u :
                case_index%5u==2u ? 1u :
                case_index%5u==3u ? width/2u : width+1u;
            std::vector<std::uint64_t> da(width),db(length_b);
            for (auto& word:da) { word=shifted_engine(); }
            for (auto& word:db) { word=shifted_engine(); }
            if (case_index%7u==0u) {
                std::fill(da.begin(),da.end(),
                          std::numeric_limits<std::uint64_t>::max());
                std::fill(db.begin(),db.end(),
                          std::numeric_limits<std::uint64_t>::max());
            }
            const BigUnsigned a(da),b(db);
            std::size_t source_ref=0u,source_width=0u,source_length_b=0u;
            unsigned source_case=0u;
            std::string a_hex,b_hex;
            if (!(b135_sources >> source_ref >> source_width >> source_case
                  >> source_length_b >> a_hex >> b_hex) ||
                source_ref!=source_ordinal++ || source_width!=width ||
                source_case!=case_index || source_length_b!=length_b)
                throw std::runtime_error(
                    "Missing or unordered B135 shifted source row");
            check_against_precomputed_hex(a,a_hex);
            check_against_precomputed_hex(b,b_hex);

            for (unsigned sa:{0u,1u,7u,31u,63u}) {
                for (unsigned sb:{0u,1u,7u,31u,63u}) {
                    std::ifstream& golden=b135_golden[b135_part(ordinal)];
                    std::size_t ref_ordinal=0u,ref_width=0u;
                    unsigned ref_case=0u,ref_sa=0u,ref_sb=0u;
                    int ref_comparison=0;
                    std::string sum_hex,difference_hex;
                    if (!(golden >> ref_ordinal >> ref_width >> ref_case
                          >> ref_sa >> ref_sb >> ref_comparison
                          >> sum_hex >> difference_hex) ||
                        ref_ordinal!=ordinal++ || ref_width!=width ||
                        ref_case!=case_index || ref_sa!=sa || ref_sb!=sb ||
                        ref_comparison < -1 || ref_comparison > 1)
                        throw std::runtime_error(
                            "Missing or unordered B135 shifted oracle row");
                    check_against_precomputed_hex(
                        BigUnsigned::add_shifted(
                            a.limbs(),sa,b.limbs(),sb),sum_hex);
                    if (ref_comparison>=0) {
                        check_against_precomputed_hex(
                            BigUnsigned::sub_shifted(
                                a.limbs(),sa,b.limbs(),sb),difference_hex);
                    } else {
                        check_against_precomputed_hex(
                            BigUnsigned::sub_shifted(
                                b.limbs(),sb,a.limbs(),sa),difference_hex);
                    }
                }
            }
            check_against_precomputed_hex(a,a_hex);
            check_against_precomputed_hex(b,b_hex);
        }
    }
    if (source_ordinal!=264u || ordinal!=6600u)
        throw std::runtime_error("Wrong B135 shifted oracle row count");
    std::string b135_source_extra;
    if (b135_sources >> b135_source_extra)
        throw std::runtime_error("Extra B135 shifted source rows");
    for (auto& golden:b135_golden) {
        std::string extra;
        if (golden >> extra)
            throw std::runtime_error("Extra B135 shifted oracle rows");
    }
    // R38: shifted subtraction rejects a negative result.
    ARIADNE_TEST_ASSERT(([] {
        const BigUnsigned one(std::uint64_t{1}),two(std::uint64_t{2});
        try {
            (void)BigUnsigned::sub_shifted(one.limbs(),0u,two.limbs(),0u);
        } catch (std::underflow_error const&) {
            return true;
        }
        return false;
    }()));
}

void TestBigUnsigned::test_multiply_kernel_paths()
{
    const auto max = std::numeric_limits<std::uint64_t>::max();
    std::mt19937_64 engine(0xB3B3B3ULL);
    BigUnsigned::Scratch scratch;
    BigUnsigned::Scratch karatsuba_workspace;

    // B52: data-driven 32x32 GCC x86-64 full-product dispatch
    // to the baseline MULQ/ADCQ kernel. GMP is the independent
    // product oracle. Exercise 31/32/33-limb boundaries, including
    // noncanonical zero data, sparse rows, all-ones, extra output
    // capacity and low products. No workload fixtures are changed.
    std::ifstream b136_b52(
        std::string(ARIADNE_TEST_DATA_DIR)+"/big_unsigned_b136_b52.tsv");
    if (!b136_b52) throw std::runtime_error("Missing B136 B52 TSV");
    std::string b136_b52_header;
    while (b136_b52.peek()=='#') std::getline(b136_b52,b136_b52_header);
    std::mt19937_64 b52_engine(0xB5203232ULL);
    for (const std::size_t n:{31u,32u,33u}) {
        for (const std::size_t m:{31u,32u,33u}) {
            for (unsigned mode=0u;mode<32u;++mode) {
                std::vector<std::uint64_t> left(n),right(m);
                for (auto& digit:left) digit=b52_engine();
                for (auto& digit:right) digit=b52_engine();
                if (mode%8u==0u) {
                    std::fill(left.begin(),left.end(),max);
                    std::fill(right.begin(),right.end(),max);
                } else if (mode%8u==1u) {
                    std::fill(left.begin(),left.end(),0u);
                    std::fill(right.begin(),right.end(),0u);
                    left.back()=1u;
                    right.back()=1u;
                } else if (mode%8u==2u) {
                    for (std::size_t i=0u;i<n-1u;i+=2u)
                        left[i]=0u;
                    for (std::size_t i=0u;i<m-1u;i+=2u)
                        right[i]=0u;
                } else if (mode%8u==3u) {
                    std::fill(left.begin(),left.end(),0u);
                    std::fill(right.begin(),right.end(),0u);
                    left.front()=max;
                    right.front()=max;
                } else if (mode%8u==4u) {
                    for (std::size_t i=1u;i+1u<n;++i)
                        left[i]=0u;
                    for (std::size_t i=1u;i+1u<m;++i)
                        right[i]=0u;
                }
                const BigUnsigned lhs(left),rhs(right);
                std::size_t ref_n=0u,ref_m=0u,ref_low=0u;
                unsigned ref_mode=0u;
                std::string lhs_hex,rhs_hex,product_hex,low_hex;
                if (!(b136_b52 >> ref_n >> ref_m >> ref_mode
                      >> lhs_hex >> rhs_hex >> product_hex
                      >> ref_low >> low_hex) ||
                    ref_n!=n || ref_m!=m || ref_mode!=mode)
                    throw std::runtime_error(
                        "Missing or unordered B136 B52 TSV row");
                check_against_precomputed_hex(lhs,lhs_hex);
                check_against_precomputed_hex(rhs,rhs_hex);
                const std::size_t full=lhs.size()+rhs.size();
                std::vector<std::uint64_t> selected(full,max);
                BigUnsigned::mul_into(
                    lhs.limbs(),rhs.limbs(),selected);
                check_against_precomputed_hex(
                    BigUnsigned(selected),product_hex);
                check_against_precomputed_hex(mul(lhs,rhs),product_hex);

                std::vector<std::uint64_t> fused(full,max);
                BigUnsigned::mul_into_fused(
                    lhs.limbs(),rhs.limbs(),fused);
                ARIADNE_TEST_ASSERT(selected==fused);

                std::vector<std::uint64_t> oversized(full+2u,max);
                BigUnsigned::mul_into(
                    lhs.limbs(),rhs.limbs(),oversized);
                check_against_precomputed_hex(
                    BigUnsigned(oversized),product_hex);
                ARIADNE_TEST_ASSERT(oversized[full]==0u);
                ARIADNE_TEST_ASSERT(oversized[full+1u]==0u);

                const std::size_t low=full==0u?0u:full-1u;
                if (ref_low!=low)
                    throw std::runtime_error("Wrong B136 B52 low length");
                std::vector<std::uint64_t> truncated(low,max);
                BigUnsigned::mul_low_into(
                    lhs.limbs(),rhs.limbs(),truncated);
                check_against_precomputed_hex(
                    BigUnsigned(truncated),low_hex);
                check_against_precomputed_hex(lhs,lhs_hex);
                check_against_precomputed_hex(rhs,rhs_hex);
            }
        }
    }
    std::string b136_b52_extra;
    if (b136_b52 >> b136_b52_extra)
        throw std::runtime_error("Extra B136 B52 TSV rows");

    // B53: independent GMP oracle for x86 GCC fused and GNU/Clang
    // assembly zero-row skipping. Include 31/32/33 B52 boundaries,
    // 63/64/65 adaptive-density thresholds, unbalanced rows and
    // truncated/oversized outputs. Raw views include noncanonical
    // zero-high limbs, so skipping must not depend on normalization.
    std::ifstream b137_b53(
        std::string(ARIADNE_TEST_DATA_DIR)+"/big_unsigned_b137_b53.tsv");
    if (!b137_b53) throw std::runtime_error("Missing B137 B53 TSV");
    std::string b137_b53_header;
    while (b137_b53.peek()=='#') std::getline(b137_b53,b137_b53_header);
    std::mt19937_64 b53_engine(0xB530B040ULL);
    for (const auto& [n,m]:
            {std::pair<std::size_t,std::size_t>{1u,1u},
             {7u,64u}, {16u,16u}, {31u,32u}, {32u,31u},
             {32u,32u}, {32u,33u}, {33u,32u}, {63u,64u},
             {64u,63u}, {64u,64u}, {64u,65u}, {65u,64u},
             {64u,128u}, {128u,64u}, {128u,128u}}) {
        for (unsigned mode=0u;mode<24u;++mode) {
            std::vector<std::uint64_t> left(n),right(m);
            for (auto& word:left) word=b53_engine();
            for (auto& word:right) word=b53_engine();
            switch (mode%8u) {
            case 0u:
                std::fill(left.begin(),left.end(),max);
                std::fill(right.begin(),right.end(),max);
                break;
            case 1u:
                std::fill(left.begin(),left.end(),0u);
                std::fill(right.begin(),right.end(),0u);
                break;
            case 2u:
                for (std::size_t i=0u;i<n;i+=2u) left[i]=0u;
                for (std::size_t i=0u;i<m;i+=2u) right[i]=0u;
                break;
            case 3u:
                std::fill(left.begin()+1u,left.end(),0u);
                std::fill(right.begin()+1u,right.end(),0u);
                break;
            case 4u:
                std::fill_n(left.begin(),n/2u,0u);
                std::fill_n(right.begin(),m/2u,0u);
                break;
            case 5u:
                std::fill(left.begin(),left.end(),0u);
                left.back()=max;
                std::fill(right.begin(),right.end(),0u);
                right.back()=max;
                break;
            case 6u:
                for (std::size_t i=0u;i<n;i+=8u) left[i]=0u;
                for (std::size_t i=0u;i<m;i+=8u) right[i]=0u;
                break;
            default: break;
            }
            const BigUnsigned lhs(left),rhs(right);
            std::size_t ref_n=0u,ref_m=0u;
            unsigned ref_mode=0u;
            std::string lhs_hex,rhs_hex;
            if (!(b137_b53 >> ref_n >> ref_m >> ref_mode
                  >> lhs_hex >> rhs_hex) ||
                ref_n!=n || ref_m!=m || ref_mode!=mode)
                throw std::runtime_error(
                    "Missing or unordered B137 B53 TSV row");
            check_against_precomputed_hex(lhs,lhs_hex);
            check_against_precomputed_hex(rhs,rhs_hex);
            for (const std::size_t length:
                    {std::size_t{0u},std::size_t{1u},
                     (n+m)/2u,n+m-1u,n+m,n+m+2u}) {
                std::size_t ref_length=0u;
                std::string low_hex;
                if (!(b137_b53 >> ref_length >> low_hex) ||
                    ref_length!=length)
                    throw std::runtime_error(
                        "Wrong B137 B53 low-product length");
                std::vector<std::uint64_t> selected(length,max);
                std::vector<std::uint64_t> fused(length,max);
                std::vector<std::uint64_t> carry(length,max);
                BigUnsigned::mul_low_into(left,right,selected);
                BigUnsigned::mul_low_into_fused(left,right,fused);
                BigUnsigned::mul_low_into_x86_64_carry(
                    left,right,carry);
                ARIADNE_TEST_ASSERT(selected==fused);
                ARIADNE_TEST_ASSERT(selected==carry);
                check_against_precomputed_hex(
                    BigUnsigned(selected),low_hex);
            }
        }
    }
    std::string b137_b53_extra;
    if (b137_b53 >> b137_b53_extra)
        throw std::runtime_error("Extra B137 B53 TSV rows");

    // B138: preserve all 6,000 original general-kernel visits using six
    // sequential source-keyed TSV segments, while keeping one RNG stream.
    std::array<std::ifstream,6u> b138_general;
    for (std::size_t part=0u;part<b138_general.size();++part) {
        const std::string path=std::string(ARIADNE_TEST_DATA_DIR)+
            "/big_unsigned_b138_general_"+std::to_string(part)+".tsv";
        b138_general[part].open(path);
        if (!b138_general[part])
            throw std::runtime_error("Missing B138 general multiply TSV: "+path);
        std::string header;
        while (b138_general[part].peek()=='#')
            std::getline(b138_general[part],header);
    }
    // Multiple writes reuse the same result space; both implementations
    // must overwrite all digits, including untouched high destinations.
    for (unsigned repetition=0u; repetition<6000u; ++repetition) {
        const std::size_t n = engine()%38u;
        const std::size_t m = engine()%38u;
        std::array<std::uint64_t,38> da{},db{};
        for (std::size_t i=0u; i<n; ++i) da[i] = engine();
        for (std::size_t i=0u; i<m; ++i) db[i] = engine();
        if (repetition%7u==0u) std::fill_n(da.begin(),n,max);
        if (repetition%11u==0u) std::fill_n(db.begin(),m,max);
        if (repetition%13u==0u) std::fill_n(da.begin(),n,0u);
        if (repetition%17u==0u) std::fill_n(db.begin(),m,0u);
        const BigUnsigned lhs(std::span(da.data(),n));
        const BigUnsigned rhs(std::span(db.data(),m));
        const std::size_t low = repetition%5u==0u ? 0u :
            repetition%5u==1u ? 1u : repetition%5u==2u ? 4u :
            repetition%5u==3u ? 5u : engine()%80u;
        std::ifstream& golden=b138_general[repetition/1000u];
        unsigned ref_repetition=0u;
        std::size_t ref_n=0u,ref_m=0u,ref_low=0u;
        std::string lhs_hex,rhs_hex,product_hex,low_hex;
        if (!(golden >> ref_repetition >> ref_n >> ref_m >> ref_low
              >> lhs_hex >> rhs_hex >> product_hex >> low_hex) ||
            ref_repetition!=repetition || ref_n!=n || ref_m!=m ||
            ref_low!=low)
            throw std::runtime_error(
                "Missing or unordered B138 general multiply TSV row");
        check_against_precomputed_hex(lhs,lhs_hex);
        check_against_precomputed_hex(rhs,rhs_hex);

        const std::size_t length = lhs.size()+rhs.size()+2u;
        auto selected = scratch.acquire(length);
        BigUnsigned::mul_into(lhs.limbs(),rhs.limbs(),selected);
        check_against_precomputed_hex(BigUnsigned(selected),product_hex);
        const auto selected_copy = std::vector<std::uint64_t>(
            selected.begin(),selected.end());

        selected = scratch.acquire(length);
        BigUnsigned::mul_into_portable(lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),selected_copy.begin()));
        check_against_precomputed_hex(BigUnsigned(selected),product_hex);

        // The experimental Comba kernel must match the B3 portable
        // reference word-for-word, also for oversized full destinations.
        selected = scratch.acquire(length);
        BigUnsigned::mul_into_comba(lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),selected_copy.begin()));
        check_against_precomputed_hex(BigUnsigned(selected),product_hex);

        selected = scratch.acquire(length);
        BigUnsigned::mul_into_fused(lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),selected_copy.begin()));
        check_against_precomputed_hex(BigUnsigned(selected),product_hex);

        selected = scratch.acquire(length);
        BigUnsigned::mul_into_arm64_carry(lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),selected_copy.begin()));
        check_against_precomputed_hex(BigUnsigned(selected),product_hex);

        // The explicit x86-64 carry-chain uses the same GMP oracle,
        // including all-ones, empty, asymmetric and reused output cases.
        selected=scratch.acquire(length);
        BigUnsigned::mul_into_x86_64_carry(
            lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),selected_copy.begin()));
        check_against_precomputed_hex(BigUnsigned(selected),product_hex);

        // Explicit B3 schoolbook remains available even if the
        // selected kernel changes on GCC x86-64.
        selected = scratch.acquire(length);
        BigUnsigned::mul_into_schoolbook(lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),selected_copy.begin()));
        check_against_precomputed_hex(BigUnsigned(selected),product_hex);

        const auto needed=BigUnsigned::karatsuba_scratch_limbs(
            lhs.size(),rhs.size());
        auto work=karatsuba_workspace.acquire(needed);
        selected=scratch.acquire(length);
        BigUnsigned::mul_into_karatsuba(
            lhs.limbs(),rhs.limbs(),selected,work);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),selected_copy.begin()));
        check_against_precomputed_hex(BigUnsigned(selected),product_hex);

        selected = scratch.acquire(low);
        BigUnsigned::mul_low_into(lhs.limbs(),rhs.limbs(),selected);
        const auto low_copy = std::vector<std::uint64_t>(
            selected.begin(),selected.end());
        check_against_precomputed_hex(BigUnsigned(selected),low_hex);
        selected = scratch.acquire(low);
        BigUnsigned::mul_low_into_portable(
            lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),low_copy.begin()));

        selected = scratch.acquire(low);
        BigUnsigned::mul_low_into_comba(lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),low_copy.begin()));
        check_against_precomputed_hex(BigUnsigned(selected),low_hex);

        selected=scratch.acquire(low);
        BigUnsigned::mul_low_into_fused(lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),low_copy.begin()));
        check_against_precomputed_hex(BigUnsigned(selected),low_hex);

        selected=scratch.acquire(low);
        BigUnsigned::mul_low_into_arm64_carry(
            lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),low_copy.begin()));
        check_against_precomputed_hex(BigUnsigned(selected),low_hex);

        selected=scratch.acquire(low);
        BigUnsigned::mul_low_into_x86_64_carry(
            lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),low_copy.begin()));
        check_against_precomputed_hex(BigUnsigned(selected),low_hex);

        selected=scratch.acquire(low);
        BigUnsigned::mul_low_into_schoolbook(
            lhs.limbs(),rhs.limbs(),selected);
        ARIADNE_TEST_ASSERT(std::equal(
            selected.begin(),selected.end(),low_copy.begin()));
        check_against_precomputed_hex(BigUnsigned(selected),low_hex);
    }

    for (std::size_t part=0u;part<b138_general.size();++part) {
        std::string extra;
        if (b138_general[part] >> extra)
            throw std::runtime_error("Extra B138 general multiply TSV rows");
    }

    const std::array<std::uint64_t,1> one{1u};
    std::array<std::uint64_t,1> too_short{};
    bool selected_threw = false;
    bool portable_threw = false;
    bool comba_threw = false;
    bool fused_threw = false;
    bool schoolbook_threw = false;
    bool carry_threw = false;
    bool x86_carry_threw = false;
    bool karatsuba_threw = false;
    ARIADNE_TEST_ASSERT(
        !BigUnsigned::selected_multiply_uses_fused()
        || BigUnsigned::has_native_multiply());
#if defined(__x86_64__) && (defined(__GNUC__) || defined(__clang__))
    ARIADNE_TEST_ASSERT(BigUnsigned::has_x86_64_carry_multiply());
#endif
#if defined(__aarch64__) && (defined(__GNUC__) || defined(__clang__))
    // Native ARM64 uses the specialized carry-chain, not common fused.
    ARIADNE_TEST_ASSERT(BigUnsigned::has_arm64_multiply());
    ARIADNE_TEST_ASSERT(!BigUnsigned::selected_multiply_uses_fused());
#endif
    try {
        BigUnsigned::mul_into(one,one,too_short);
    } catch (std::length_error const&) {
        selected_threw = true;
    }
    try {
        BigUnsigned::mul_into_portable(one,one,too_short);
    } catch (std::length_error const&) {
        portable_threw = true;
    }
    try {
        BigUnsigned::mul_into_comba(one,one,too_short);
    } catch (std::length_error const&) {
        comba_threw = true;
    }
    try {
        BigUnsigned::mul_into_fused(one,one,too_short);
    } catch (std::length_error const&) {
        fused_threw = true;
    }
    try {
        BigUnsigned::mul_into_arm64_carry(one,one,too_short);
    } catch (std::length_error const&) {
        carry_threw = true;
    }
    try {
        BigUnsigned::mul_into_x86_64_carry(one,one,too_short);
    } catch (std::length_error const&) {
        x86_carry_threw = true;
    }
    try {
        BigUnsigned::mul_into_schoolbook(one,one,too_short);
    } catch (std::length_error const&) {
        schoolbook_threw = true;
    }
    try {
        BigUnsigned::mul_into_karatsuba(one,one,too_short,{});
    } catch (std::length_error const&) {
        karatsuba_threw = true;
    }
    ARIADNE_TEST_ASSERT(selected_threw && portable_threw &&
                        comba_threw && fused_threw && carry_threw &&
                        x86_carry_threw && schoolbook_threw &&
                        karatsuba_threw);

    // Missing recursive scratch must fail before altering output.
    const std::array<std::uint64_t,17> long_ones{
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
    };
    std::array<std::uint64_t,34> guarded;
    guarded.fill(max);
    ARIADNE_TEST_ASSERT(
        BigUnsigned::karatsuba_scratch_limbs(17u,17u)>0u);
    ARIADNE_TEST_ASSERT(
        BigUnsigned::karatsuba_scratch_limbs(17u,17u,8u) ==
        BigUnsigned::karatsuba_scratch_limbs(17u,17u));
    ARIADNE_TEST_ASSERT(
        BigUnsigned::karatsuba_scratch_limbs(17u,17u,64u)==0u);
    bool invalid_cutoff_threw=false;
    try {
        (void)BigUnsigned::karatsuba_scratch_limbs(17u,17u,7u);
    } catch (std::invalid_argument const&) {
        invalid_cutoff_threw=true;
    }
    ARIADNE_TEST_ASSERT(invalid_cutoff_threw);
    bool scratch_threw=false;
    try {
        BigUnsigned::mul_into_karatsuba(
            long_ones,long_ones,guarded,{});
    } catch (std::length_error const&) {
        scratch_threw=true;
    }
    ARIADNE_TEST_ASSERT(scratch_threw);
    for (auto value: guarded) {
        ARIADNE_TEST_ASSERT(value==max);
    }

    // B139: B8 retains all 700 original sources and four cutoff visits
    // per source, with an independent exact product TSV.
    std::ifstream b139_b8(
        std::string(ARIADNE_TEST_DATA_DIR)+"/big_unsigned_b139_b8.tsv");
    if (!b139_b8) throw std::runtime_error("Missing B139 B8 TSV");
    std::string b139_b8_header;
    while (b139_b8.peek()=='#') std::getline(b139_b8,b139_b8_header);
    std::mt19937_64 wider_engine(0xB8B8B8ULL);
    for (unsigned repetition=0u; repetition<700u; ++repetition) {
        const std::size_t n=repetition%7u==0u ? 128u :
                            9u+wider_engine()%95u;
        const std::size_t m=repetition%11u==0u ? 1u :
                            9u+wider_engine()%95u;
        std::array<std::uint64_t,128> da{},db{};
        for (std::size_t i=0u;i<n;++i) da[i]=wider_engine();
        for (std::size_t i=0u;i<m;++i) db[i]=wider_engine();
        if (repetition%13u==0u) std::fill_n(da.begin(),n,max);
        if (repetition%17u==0u) std::fill_n(db.begin(),m,max);
        if (repetition%19u==0u) std::fill_n(da.begin(),n,0u);
        if (repetition%23u==0u) std::fill_n(db.begin(),m,0u);
        const BigUnsigned lhs(std::span(da.data(),n));
        const BigUnsigned rhs(std::span(db.data(),m));
        unsigned ref_repetition=0u;
        std::size_t ref_n=0u,ref_m=0u;
        std::string lhs_hex,rhs_hex,product_hex;
        if (!(b139_b8 >> ref_repetition >> ref_n >> ref_m
              >> lhs_hex >> rhs_hex >> product_hex) ||
            ref_repetition!=repetition || ref_n!=n || ref_m!=m)
            throw std::runtime_error("Missing or unordered B139 B8 TSV row");
        check_against_precomputed_hex(lhs,lhs_hex);
        check_against_precomputed_hex(rhs,rhs_hex);
        for (std::size_t cutoff: {8u,16u,32u,64u}) {
            auto product_view=scratch.acquire(lhs.size()+rhs.size()+2u);
            auto work=karatsuba_workspace.acquire(
                BigUnsigned::karatsuba_scratch_limbs(
                    lhs.size(),rhs.size(),cutoff));
            BigUnsigned::mul_into_karatsuba(
                lhs.limbs(),rhs.limbs(),product_view,work,cutoff);
            check_against_precomputed_hex(
                BigUnsigned(product_view),product_hex);
            for (std::size_t i=lhs.size()+rhs.size();
                    i<product_view.size();++i) {
                ARIADNE_TEST_ASSERT(product_view[i]==0u);
            }
        }
    }

    std::string b139_b8_extra;
    if (b139_b8 >> b139_b8_extra)
        throw std::runtime_error("Extra B139 B8 TSV rows");

    // B139: fixed Karatsuba sources continue the same wider_engine stream.
    std::ifstream b139_fixed(
        std::string(ARIADNE_TEST_DATA_DIR)+
        "/big_unsigned_b139_karatsuba_fixed.tsv");
    if (!b139_fixed) throw std::runtime_error("Missing B139 fixed Karatsuba TSV");
    std::string b139_fixed_header;
    while (b139_fixed.peek()=='#')
        std::getline(b139_fixed,b139_fixed_header);
    std::size_t b139_fixed_index=0u;
    // Karatsuba output is prefilled with nonzero data to prove that the
    // z0/z2 recursion fully overwrites even unbalanced and odd splits.
    // The scratch guard limbs expose any off-by-one workspace overruns.
    // Reuse both buffers across cutoffs, including a cutoff with no split.
    for (auto [n,m]: {std::pair<std::size_t,std::size_t>{17u,17u},
                       {31u,33u}, {63u,67u}, {79u,80u}, {127u,129u},
                       {255u,257u}, {1u,180u}, {8u,8u}}) {
        std::vector<std::uint64_t> da(n,max),db(m,max);
        for (std::size_t i=0u;i<n;++i) {
            if (i%5u==1u) da[i]=wider_engine();
        }
        for (std::size_t i=0u;i<m;++i) {
            if (i%7u==2u) db[i]=wider_engine();
        }
        const BigUnsigned lhs(da),rhs(db);
        std::size_t ref_index=0u,ref_n=0u,ref_m=0u;
        std::string lhs_hex,rhs_hex,product_hex;
        if (!(b139_fixed >> ref_index >> ref_n >> ref_m
              >> lhs_hex >> rhs_hex >> product_hex) ||
            ref_index!=b139_fixed_index++ || ref_n!=n || ref_m!=m)
            throw std::runtime_error(
                "Missing or unordered B139 fixed Karatsuba TSV row");
        check_against_precomputed_hex(lhs,lhs_hex);
        check_against_precomputed_hex(rhs,rhs_hex);
        std::vector<std::uint64_t> product(n+m+3u,max);
        for (std::size_t cutoff: {8u,16u,32u,64u,128u}) {
            const auto needed=BigUnsigned::karatsuba_scratch_limbs(n,m,cutoff);
            std::vector<std::uint64_t> guards(needed+2u,max);
            const auto work=BigUnsigned::MutableLimbView(
                guards.data()+1u,needed);
            std::fill(product.begin(),product.end(),max);
            BigUnsigned::mul_into_karatsuba(da,db,product,work,cutoff);
            check_against_precomputed_hex(
                BigUnsigned(product),product_hex);
            for (std::size_t i=n+m;i<product.size();++i) {
                ARIADNE_TEST_ASSERT(product[i]==0u);
            }
            ARIADNE_TEST_ASSERT(guards.front()==max);
            ARIADNE_TEST_ASSERT(guards.back()==max);
        }
    }
    if (b139_fixed_index!=8u)
        throw std::runtime_error("Wrong B139 fixed Karatsuba row count");
    std::string b139_fixed_extra;
    if (b139_fixed >> b139_fixed_extra)
        throw std::runtime_error("Extra B139 fixed Karatsuba TSV rows");

}

void TestBigUnsigned::test_division_limb()
{
    const auto max=std::numeric_limits<std::uint64_t>::max();
    // Fixed independent cases that exercise the 64-bit remainder carry.
    const std::array<std::uint64_t,2> all_ones{max,max};
    const std::array<std::uint64_t,2> power_of_base{0u,1u};
    const auto [whole,zero_remainder]=divrem_limb(BigUnsigned(all_ones),max);
    ARIADNE_TEST_ASSERT(compare(whole,BigUnsigned(
        std::array<std::uint64_t,2>{1u,1u}))==0);
    ARIADNE_TEST_ASSERT(zero_remainder==0u);
    const auto [one,remainder]=divrem_limb(BigUnsigned(power_of_base),max);
    ARIADNE_TEST_ASSERT(compare(one,BigUnsigned(std::uint64_t{1}))==0);
    ARIADNE_TEST_ASSERT(remainder==1u);
    // A nonzero high half of the remainder must survive on LLP64 hosts.
    const std::uint64_t wide_remainder=(std::uint64_t{1}<<40u)+17u;
    const auto [zero_quotient,wide_rem]=divrem_limb(
        BigUnsigned(wide_remainder),max);
    ARIADNE_TEST_ASSERT(zero_quotient.is_zero());
    ARIADNE_TEST_ASSERT(wide_rem==wide_remainder);
    ARIADNE_TEST_ASSERT(divrem_limb(BigUnsigned{},1u).first.is_zero());

    std::ifstream random_div_golden(ARIADNE_BIG_UNSIGNED_DIV_LIMB_RANDOM_ORACLE_PATH);
    if (!random_div_golden)
        throw std::runtime_error("Missing B115 random limb division TSV");
    std::string random_div_header;
    while (random_div_golden.peek()=='#')
        std::getline(random_div_golden,random_div_header);
    std::mt19937_64 engine(0xD1D1D1ULL);
    BigUnsigned::Scratch scratch;
    // Retain all 6,000 original RNG visits and size/transform cases.
    for (unsigned repetition=0u;repetition<6000u;++repetition) {
        const std::size_t count=repetition%23u==0u ? 0u :
                                repetition%29u==0u ? 64u :
                                repetition%7u+1u;
        std::vector<std::uint64_t> digits(count);
        for (auto& limb:digits) limb=engine();
        if (repetition%13u==0u) {
            std::fill(digits.begin(),digits.end(),max);
        }
        if (repetition%17u==0u) {
            std::fill(digits.begin(),digits.end(),0u);
        }
        if (repetition%19u==0u && !digits.empty()) {
            std::fill(digits.begin(),digits.end(),0u);
            digits.back()=1u;
        }
        const BigUnsigned value(digits);
        const std::uint64_t divisor=repetition%11u==0u ? 1u :
                                    repetition%11u==1u ? max :
                                    repetition%11u==2u ? std::uint64_t{1}<<63u :
                                    repetition%11u==3u ? 2u :
                                    engine()|std::uint64_t{1};
        unsigned source_repetition=0u;
        std::uint64_t source_divisor=0u;
        std::string source_hex,quotient_hex,remainder_hex;
        if (!(random_div_golden >> source_repetition >> source_divisor
                                >> source_hex >> quotient_hex >> remainder_hex)
            || source_repetition!=repetition || source_divisor!=divisor)
            throw std::runtime_error("Missing or unordered B115 random division row");
        check_against_precomputed_hex(value,source_hex);
        const auto [actual,rem]=divrem_limb(value,divisor);
        check_against_precomputed_hex(actual,quotient_hex);
        check_against_precomputed_hex(BigUnsigned(rem),remainder_hex);
        ARIADNE_TEST_ASSERT(rem<divisor);

        auto output=scratch.acquire(value.size()+2u);
        std::fill(output.begin(),output.end(),max);
        const auto low=BigUnsigned::divrem_limb_into(
            value.limbs(),divisor,output);
        ARIADNE_TEST_ASSERT(low==rem);
        check_against_precomputed_hex(BigUnsigned(output),quotient_hex);
        ARIADNE_TEST_ASSERT(output[value.size()]==0u);
        ARIADNE_TEST_ASSERT(output[value.size()+1u]==0u);
        check_against_precomputed_hex(value,source_hex);
    }
    std::string random_div_extra;
    if (random_div_golden >> random_div_extra)
        throw std::runtime_error("Extra B115 random division TSV rows");

    // B114: independent source-keyed exact quotient and remainder for all
    // 192 original 128/64 edge visits and all 2,048 original B26 visits.
    // Retain the original input order and conditional PRNG consumption.
    std::ifstream div128_golden(ARIADNE_BIG_UNSIGNED_DIV128_ORACLE_PATH);
    if (!div128_golden) throw std::runtime_error("Missing B114 div128 TSV");
    std::string div128_header;
    while (div128_golden.peek()=='#')
        std::getline(div128_golden,div128_header);

    // Target 128/64 single-digit division including saturated quotient,
    // divisor boundary values, high=divisor-1 and wide remainders.
    const std::array<std::uint64_t,8> edge_divisors{
        1u,2u,3u,std::uint64_t{0xffffffffu},
        std::uint64_t{1}<<32u,std::uint64_t{1}<<63u,
        max-1u,max};
    const std::array<std::uint64_t,6> edge_lows{
        0u,1u,std::uint64_t{0xffffffffu},
        std::uint64_t{1}<<63u,max-1u,max};
    std::size_t edge_index=0u;
    for (const std::uint64_t divisor:edge_divisors) {
        const std::array<std::uint64_t,4> highs{
            0u,divisor>1u?1u:0u,divisor/2u,divisor-1u};
        for (const std::uint64_t high:highs) {
            for (const std::uint64_t low:edge_lows) {
                const std::array<std::uint64_t,2> dividend{low,high};
                const BigUnsigned value(dividend);
                const auto [native_q,native_r]=divrem_limb(value,divisor);
                check_div128_golden(div128_golden,'E',edge_index++,
                                    divisor,low,high,value,native_q,native_r);
                ARIADNE_TEST_ASSERT(native_r<divisor);
                if (high==divisor-1u && low==max) {
                    ARIADNE_TEST_ASSERT(native_q.limb(0u)==max);
                    ARIADNE_TEST_ASSERT(native_r==divisor-1u);
                }
            }
        }
    }

    // B26: preserve all 2,048 original 128/64 DIVQ precondition visits.
    // Original seed, conditional draws and high<divisor are unchanged.
    std::mt19937_64 wide_engine(0xB260D1B1ULL);
    for (unsigned repetition=0u;repetition<2048u;++repetition) {
        const std::uint64_t divisor=repetition%11u==0u ? max :
            repetition%11u==1u ? std::uint64_t{1}<<63u :
            repetition%11u==2u ? 3u :
            repetition%11u==3u ? std::uint64_t{1}<<32u :
            wide_engine()|std::uint64_t{1};
        const std::uint64_t high=repetition%5u==0u ?
            divisor-1u : repetition%5u==1u ? 0u :
            wide_engine()%divisor;
        const std::uint64_t low=repetition%7u==0u ?
            max : repetition%7u==1u ? 0u : wide_engine();
        const std::array<std::uint64_t,2> dividend{low,high};
        const BigUnsigned value(dividend);
        const auto [native_q,native_r]=divrem_limb(value,divisor);
        check_div128_golden(div128_golden,'B',repetition,
                            divisor,low,high,value,native_q,native_r);
        ARIADNE_TEST_ASSERT(native_r<divisor);
    }
    std::string div128_extra;
    if (div128_golden >> div128_extra)
        throw std::runtime_error("Extra B114 div128 TSV rows");

    // B48: the ARM64 32-bit-divisor fast path must return the exact
    // full quotient and 64-bit remainder, independently of the
    // production divide_wide implementation. Include saturated and
    // sparse operands across inline/heap storage boundaries and
    // the 2^32 divisor selector boundary.
    LimbDivisionGolden golden_b48(ARIADNE_BIG_UNSIGNED_B48_ORACLE_PATH);
    std::mt19937_64 narrow_engine(0xB48032D0ULL);
    for (const std::size_t count:
            {1u,2u,3u,4u,5u,8u,16u,32u,64u,96u}) {
        for (const std::uint64_t divisor:
                {1ULL,2ULL,3ULL,7ULL,9ULL,10ULL,31ULL,63ULL,
                 255ULL,65535ULL,2147483648ULL,4294967295ULL,
                 4294967296ULL,18446744073709551615ULL}) {
            for (unsigned mode=0u;mode<4u;++mode) {
                std::vector<std::uint64_t> words(count);
                for (auto& limb:words) limb=narrow_engine();
                if (mode==0u) {
                    std::fill(words.begin(),words.end(),max);
                } else if (mode==1u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=1u;
                } else if (mode==2u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.front()=max;
                    words.back()|=1u;
                }
                const BigUnsigned value(words);
                const auto [got_quotient,got_remainder]=
                    divrem_limb(value,divisor);

                golden_b48.check(count,divisor,mode,value,got_quotient,got_remainder);
                auto out=scratch.acquire(value.size()+2u);
                std::fill(out.begin(),out.end(),max);
                const auto low=BigUnsigned::divrem_limb_into(
                    value.limbs(),divisor,out);
                ARIADNE_TEST_ASSERT(low==got_remainder);
                ARIADNE_TEST_ASSERT(compare(BigUnsigned(out),got_quotient)==0);
                ARIADNE_TEST_ASSERT(out[value.size()]==0u);
                ARIADNE_TEST_ASSERT(out[value.size()+1u]==0u);
            }
        }
    }

    golden_b48.finish();

    // B49: exercise the pre-inverted 128/32 quotient path on ARM64
    // across 7/8/9 limb dispatch boundaries and longer dividends.
    // GMP is the independent oracle for full quotient, remainder and
    // all output limbs. Test divisor 1, 2^32-1 and 2^32 fallback too.
    LimbDivisionGolden golden_b49(ARIADNE_BIG_UNSIGNED_B49_ORACLE_PATH);
    std::mt19937_64 reciprocal_engine(0xB490C0DEULL);
    for (const std::size_t count:
            {7u,8u,9u,16u,32u,64u,96u,128u}) {
        for (const std::uint64_t divisor:
                {1ULL,2ULL,3ULL,5ULL,7ULL,9ULL,31ULL,
                 63ULL,255ULL,65535ULL,2147483648ULL,
                 4294967295ULL,4294967296ULL,
                 18446744073709551615ULL}) {
            for (unsigned mode=0u;mode<7u;++mode) {
                std::vector<std::uint64_t> words(count);
                for (auto& word:words) word=reciprocal_engine();
                if (mode==0u) {
                    std::fill(words.begin(),words.end(),max);
                } else if (mode==1u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=1u;
                } else if (mode==2u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.front()=max;
                    words.back()=1u;
                } else if (mode==3u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=std::uint64_t{1}<<63u;
                } else if (mode==4u) {
                    std::fill(words.begin(),words.end(),0u);
                } else if (mode==5u) {
                    for (std::size_t i=0u;i<count;i+=2u) {
                        words[i]=0u;
                    }
                }
                const BigUnsigned value(words);
                const auto [actual_q,actual_r]=
                    divrem_limb(value,divisor);

                golden_b49.check(count,divisor,mode,value,actual_q,actual_r);

                auto output=scratch.acquire(value.size()+2u);
                std::fill(output.begin(),output.end(),max);
                const auto rem_out=BigUnsigned::divrem_limb_into(
                    value.limbs(),divisor,output);
                ARIADNE_TEST_ASSERT(rem_out==actual_r);
                ARIADNE_TEST_ASSERT(compare(BigUnsigned(output),actual_q)==0);
                ARIADNE_TEST_ASSERT(output[value.size()]==0u);
                ARIADNE_TEST_ASSERT(output[value.size()+1u]==0u);
            }
        }
    }

    golden_b49.finish();

    // B57/B63: normalized pi1 selection for >=16 raw quotient
    // limbs on GCC x86-64, and universal one-limb powers-of-two.
    // Independent B57 TSV covers saturated, random, sparse and denormalized
    // operands and the 15/16/17, 31/32/33, 63/64/65 boundaries.
    LimbDivisionGolden golden_b57(ARIADNE_BIG_UNSIGNED_B57_ORACLE_PATH);
    std::mt19937_64 b57_engine(0xB570641ULL);
    for (const std::size_t count:
            {0u,1u,7u,8u,15u,16u,17u,31u,32u,33u,63u,64u,65u,96u,128u}) {
        for (const std::uint64_t divisor:
                {1ULL,2ULL,3ULL,5ULL,7ULL,35ULL,255ULL,
                 65535ULL,4294967295ULL,4294967296ULL,
                 4294967297ULL,9223372036854775807ULL,
                 9223372036854775808ULL,9223372036854775809ULL,
                 18446744073709551615ULL}) {
            for (unsigned mode=0u;mode<8u;++mode) {
                std::vector<std::uint64_t> words(count);
                for (auto& word:words) word=b57_engine();
                if (mode==0u) {
                    std::fill(words.begin(),words.end(),max);
                } else if (mode==1u) {
                    std::fill(words.begin(),words.end(),0u);
                } else if (mode==2u && !words.empty()) {
                    std::fill(words.begin(),words.end(),0u);
                    words.front()=1u;
                } else if (mode==3u && !words.empty()) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=1u;
                } else if (mode==4u && !words.empty()) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=max;
                } else if (mode==5u) {
                    for (std::size_t i=0u;i<count;i+=2u)
                        words[i]=0u;
                } else if (mode==6u) {
                    for (std::size_t i=0u;i<count;i+=3u)
                        words[i]=max;
                }
                auto out=scratch.acquire(count+2u);
                std::fill(out.begin(),out.end(),max);
                const auto limb_remainder=BigUnsigned::divrem_limb_into(
                    words,divisor,out);
                const BigUnsigned value(words);
                golden_b57.check(count,divisor,mode,value,
                                 BigUnsigned(out),limb_remainder);
                ARIADNE_TEST_ASSERT(out[count]==0u);
                ARIADNE_TEST_ASSERT(out[count+1u]==0u);
                const auto [owned_q,owned_r]=
                    divrem_limb(value,divisor);
                ARIADNE_TEST_ASSERT(compare(owned_q,BigUnsigned(out))==0);
                ARIADNE_TEST_ASSERT(owned_r==limb_remainder);
            }
        }
    }

    golden_b57.finish();

    // B61: exact single-limb remainder-only pi1 with a GCC x86-64
    // dispatch boundary at 7/8/9 limbs; use independent B61
    // source-keyed TSV for zero, saturated, sparse and random dividends.
    // Also guard the existing 63/64/65-limb quotient pi1 boundary.
    LimbDivisionGolden golden_b61(ARIADNE_BIG_UNSIGNED_B61_ORACLE_PATH);
    std::mt19937_64 b61_engine(0xB610D1A6ULL);
    for (const std::size_t count:
            {0u,1u,7u,8u,9u,15u,16u,31u,32u,63u,64u,65u,128u}) {
        for (const std::uint64_t divisor:
                {1ULL,2ULL,3ULL,5ULL,7ULL,35ULL,255ULL,
                 65535ULL,4294967295ULL,4294967296ULL,
                 4294967297ULL,9223372036854775807ULL,
                 9223372036854775808ULL,9223372036854775809ULL,
                 18446744073709551615ULL}) {
            for (unsigned mode=0u;mode<8u;++mode) {
                std::vector<std::uint64_t> words(count);
                for (auto& word:words) word=b61_engine();
                if (mode==0u) {
                    std::fill(words.begin(),words.end(),max);
                } else if (mode==1u) {
                    std::fill(words.begin(),words.end(),0u);
                } else if (mode==2u && !words.empty()) {
                    std::fill(words.begin(),words.end(),0u);
                    words.front()=1u;
                } else if (mode==3u && !words.empty()) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=1u;
                } else if (mode==4u && !words.empty()) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=max;
                } else if (mode==5u) {
                    for (std::size_t i=0u;i<count;i+=2u)
                        words[i]=0u;
                } else if (mode==6u) {
                    for (std::size_t i=0u;i<count;i+=3u)
                        words[i]=max;
                }
                const BigUnsigned value(words);
                const BigUnsigned by(divisor);
                const auto [full_q,full_r]=divrem(value,by);
                ARIADNE_TEST_ASSERT(full_r.size()<=1u);
                golden_b61.check(count,divisor,mode,value,full_q,
                                 full_r.limb(0u));
                ARIADNE_TEST_ASSERT(compare(rem(value,by),full_r)==0);
            }
        }
    }

    golden_b61.finish();

    const std::array<std::uint64_t,2> numerator{max,max};
    std::array<std::uint64_t,3> output{max,max,max};
    bool zero_threw=false,owning_zero_threw=false,short_threw=false;
    try {
        (void)BigUnsigned::divrem_limb_into(numerator,0u,output);
    } catch (std::domain_error const&) {
        zero_threw=true;
    }
    try {
        (void)divrem_limb(BigUnsigned(numerator),0u);
    } catch (std::domain_error const&) {
        owning_zero_threw=true;
    }
    try {
        (void)BigUnsigned::divrem_limb_into(
            numerator,3u,std::span(output).first(1u));
    } catch (std::length_error const&) {
        short_threw=true;
    }
    ARIADNE_TEST_ASSERT(zero_threw && owning_zero_threw && short_threw);
    for (auto limb:output) {
        ARIADNE_TEST_ASSERT(limb==max);
    }
}


void TestBigUnsigned::test_exact_division_limb()
{
    using Limb=BigUnsigned::Limb;
    const Limb max=std::numeric_limits<Limb>::max();
    // Fixed values independent of the native product and division paths.
    ARIADNE_TEST_ASSERT(to_decimal(quot_exact(
        from_decimal("340282366920938463463374607431768211455"),BigUnsigned(Limb{3})))
        =="113427455640312821154458202477256070485");
    ARIADNE_TEST_ASSERT(to_decimal(quot_exact(
        from_decimal("36893488147419103230"),BigUnsigned(max)))=="2");
    std::mt19937_64 random(0x22E1AC7ULL);
    std::vector<Limb> divisors={1u,2u,3u,6u,7u,10u,12u,17u,65537u,max,max-1u};
    for (unsigned bit=0u;bit<64u;++bit) {
        divisors.push_back(Limb{1}<<bit);
        divisors.push_back((random()>>(63u-bit))|(Limb{1}<<bit));
    }
    for (Limb odd=1u;odd<32u;odd+=2u) { divisors.push_back(odd); }
    // B118: original 6,200 source visits and RNG draws. Four independently
    // certified fixtures avoid ever computing the exact dividend with GMP.
    std::array<std::ifstream,4> golden{
        std::ifstream(ARIADNE_BIG_UNSIGNED_EXACT_LIMB_SMALL_PATH),
        std::ifstream(ARIADNE_BIG_UNSIGNED_EXACT_LIMB_65_PATH),
        std::ifstream(ARIADNE_BIG_UNSIGNED_EXACT_LIMB_129A_PATH),
        std::ifstream(ARIADNE_BIG_UNSIGNED_EXACT_LIMB_129B_PATH)};
    for (auto& file:golden) {
        if (!file) throw std::runtime_error("Missing B118 exact limb TSV");
        std::string header;
        while (file.peek()=='#') std::getline(file,header);
    }
    for (std::size_t length : {0u,1u,2u,4u,5u,8u,9u,16u,65u,129u}) {
        for (std::size_t divisor_index=0u;
                divisor_index<divisors.size();++divisor_index) {
            const Limb divisor=divisors[divisor_index];
            for (unsigned mode=0u;mode<4u;++mode) {
                std::vector<Limb> words(length);
                for (auto& word:words) { word=random(); }
                if (mode==1u) { std::fill(words.begin(),words.end(),max); }
                if (mode==2u) { std::fill(words.begin(),words.end(),Limb{0}); }
                if (mode==3u) {
                    std::fill(words.begin(),words.end(),Limb{0});
                    if (!words.empty()) { words.back()=1u; words.front()|=1u; }
                }
                const std::size_t group=length<=16u ? 0u :
                    length==65u ? 1u : divisor_index<78u ? 2u : 3u;
                std::size_t ref_length=0u,ref_index=0u;
                unsigned ref_mode=0u;
                Limb ref_divisor=0u;
                std::string quotient_hex,dividend_hex;
                if (!(golden[group] >> ref_length >> ref_index >> ref_mode
                      >> ref_divisor >> quotient_hex >> dividend_hex) ||
                    ref_length!=length || ref_index!=divisor_index ||
                    ref_mode!=mode || ref_divisor!=divisor)
                    throw std::runtime_error("Missing or unordered B118 exact limb row");
                const BigUnsigned expected_quotient(words);
                check_against_precomputed_hex(expected_quotient,quotient_hex);
                const BigUnsigned source=big_unsigned_from_golden_hex(dividend_hex);
                const auto source_before=source;
                const std::size_t count=source.size();
                std::vector<Limb> output(count+2u,max);
                BigUnsigned::divexact_limb_into(source.limbs(),divisor,output);
                check_against_precomputed_hex(BigUnsigned(output),quotient_hex);
                check_against_precomputed_hex(
                    quot_exact(source,BigUnsigned(divisor)),quotient_hex);
                check_against_precomputed_hex(source,dividend_hex);
                ARIADNE_TEST_ASSERT(compare(source,source_before)==0);
                ARIADNE_TEST_ASSERT(output[count]==0u && output[count+1u]==0u);
            }
        }
    }
    for (auto& file:golden) {
        std::string extra;
        if (file >> extra) throw std::runtime_error("Extra B118 exact limb TSV rows");
    }
    const std::array<Limb,2u> input={1u,1u};
    std::array<Limb,2u> output={max,max};
    bool zero_threw=false,short_threw=false;
    try { BigUnsigned::divexact_limb_into(input,0u,output); }
    catch (std::domain_error const&) { zero_threw=true; }
    try { BigUnsigned::divexact_limb_into(input,3u,std::span(output).first(1u)); }
    catch (std::length_error const&) { short_threw=true; }
    ARIADNE_TEST_ASSERT(zero_threw && short_threw && output[0u]==max && output[1u]==max);
    for (Limb divisor:std::array<Limb,4u>{2u,3u,6u,max}) {
        bool inexact_threw=false;
        try { BigUnsigned::divexact_limb_into(input,divisor,output); }
        catch (std::domain_error const&) { inexact_threw=true; }
        ARIADNE_TEST_ASSERT(inexact_threw);
    }
}

void TestBigUnsigned::test_product_accumulation()
{
    using Limb=BigUnsigned::Limb;
    const Limb max=std::numeric_limits<Limb>::max();
    const std::array<Limb,1u> one={1u};
    std::array<Limb,2u> fixed={max,max};
    ARIADNE_TEST_ASSERT(BigUnsigned::accumulate_product_into(one,one,fixed,false)==1u);
    ARIADNE_TEST_ASSERT(fixed[0u]==0u && fixed[1u]==0u);
    ARIADNE_TEST_ASSERT(BigUnsigned::accumulate_product_into(one,one,fixed,true)==1u);
    ARIADNE_TEST_ASSERT(fixed[0u]==max && fixed[1u]==max);

    // Keep all 1,728 original seeded visits and their full RNG consumption.
    // Reference values and source keys are independently regenerated by
    // generator_big_unsigned_accumulation.cpp, not computed with GMP here.
    std::ifstream golden(ARIADNE_BIG_UNSIGNED_ACCUMULATION_ORACLE_PATH);
    if (!golden) throw std::runtime_error("Missing BigUnsigned accumulation golden TSV");
    std::string header;
    while (golden.peek()=='#') std::getline(golden,header);
    std::mt19937_64 random(0x23ACC00ULL);
    for (std::size_t length : {0u,1u,2u,4u,5u,8u,9u,16u,64u}) {
        for (unsigned trial=0u;trial<32u;++trial) {
            const std::size_t other=trial%3u==0u ? 1u : length;
            std::vector<Limb> x(length),y(other);
            for (auto& word:x) { word=random(); }
            for (auto& word:y) { word=random(); }
            if (trial%4u==0u) { std::fill(x.begin(),x.end(),max); }
            if (trial%4u==1u) {
                std::fill(x.begin(),x.end(),Limb{0});
                if (!x.empty()) { x.back()=1u; x.front()|=1u; }
            }
            if (trial%4u==2u) { std::fill(y.begin(),y.end(),Limb{0}); }
            const BigUnsigned left(x),right(y);
            for (std::size_t extra : {0u,1u,3u}) {
                const std::size_t count=left.size()+right.size()+extra;
                std::vector<Limb> original(count);
                for (auto& word:original) { word=random(); }
                if (trial%3u==0u) { std::fill(original.begin(),original.end(),Limb{0}); }
                if (trial%3u==1u) { std::fill(original.begin(),original.end(),max); }
                for (bool subtract : {false,true}) {
                    std::size_t row_length=0u,row_extra=0u;
                    unsigned row_trial=0u,row_subtract=0u,row_escaped=0u;
                    std::string lhs_hex,rhs_hex,initial_hex,result_hex;
                    if (!(golden >> row_length >> row_trial >> row_extra >> row_subtract
                                 >> lhs_hex >> rhs_hex >> initial_hex
                                 >> result_hex >> row_escaped))
                        throw std::runtime_error("Missing or malformed BigUnsigned accumulation TSV row");
                    if (row_length!=length || row_trial!=trial || row_extra!=extra
                        || row_subtract!=static_cast<unsigned>(subtract) || row_escaped>1u)
                        throw std::runtime_error("Out-of-order BigUnsigned accumulation TSV row");
                    check_against_precomputed_hex(left,lhs_hex);
                    check_against_precomputed_hex(right,rhs_hex);
                    check_against_precomputed_hex(BigUnsigned(original),initial_hex);

                    auto output=original;
                    const Limb escaped=BigUnsigned::accumulate_product_into(
                        left.limbs(),right.limbs(),output,subtract);
                    ARIADNE_TEST_ASSERT(escaped==static_cast<Limb>(row_escaped));
                    ARIADNE_TEST_ASSERT(output.size()==count);
                    check_against_precomputed_hex(BigUnsigned(output),result_hex);
                    check_against_precomputed_hex(left,lhs_hex);
                    check_against_precomputed_hex(right,rhs_hex);
                    check_against_precomputed_hex(BigUnsigned(original),initial_hex);
                }
            }
        }
    }
    std::string extra_row;
    if (golden >> extra_row)
        throw std::runtime_error("Extra BigUnsigned accumulation TSV rows");

    // Preserve the original error/unchanged-destination regression.
    bool short_threw=false;
    fixed={7u,11u};
    try { (void)BigUnsigned::accumulate_product_into(one,one,std::span(fixed).first(1u),false); }
    catch (std::length_error const&) { short_threw=true; }
    ARIADNE_TEST_ASSERT(short_threw && fixed[0u]==7u && fixed[1u]==11u);
}

void TestBigUnsigned::test_division()
{
    const auto max=std::numeric_limits<std::uint64_t>::max();
    const BigUnsigned zero;
    const BigUnsigned one(std::uint64_t{1});
    const BigUnsigned radix(std::array<std::uint64_t,2>{0u,1u});
    const BigUnsigned radix_squared(
        std::array<std::uint64_t,3>{0u,0u,1u});
    const auto [q0,r0]=divrem(zero,radix);
    ARIADNE_TEST_ASSERT(q0.is_zero() && r0.is_zero());
    const auto [q1,r1]=divrem(one,radix);
    ARIADNE_TEST_ASSERT(q1.is_zero() && compare(r1,one)==0);
    const auto [q2,r2]=divrem(radix,radix);
    ARIADNE_TEST_ASSERT(compare(q2,one)==0 && r2.is_zero());
    const auto [q3,r3]=divrem(radix_squared,radix);
    ARIADNE_TEST_ASSERT(compare(q3,radix)==0 && r3.is_zero());

    std::array<std::ifstream,6u> b140_general;
    for (std::size_t part=0u;part<b140_general.size();++part) {
        const std::string path=std::string(ARIADNE_TEST_DATA_DIR)+
            "/big_unsigned_b140_division_"+std::to_string(part)+".tsv";
        b140_general[part].open(path);
        if (!b140_general[part])
            throw std::runtime_error("Missing B140 general division TSV: "+path);
        std::string header;
        while (b140_general[part].peek()=='#')
            std::getline(b140_general[part],header);
    }
    std::mt19937_64 engine(0xB15D1A1DULL);
    for (unsigned repetition=0u;repetition<6000u;++repetition) {
        const std::size_t n=repetition%31u==0u?96u:
                             repetition%7u==0u?32u:engine()%17u;
        const std::size_t m=repetition%37u==0u?1u:
                             repetition%11u==0u?96u:
                             1u+engine()%20u;
        std::vector<std::uint64_t> da(n),db(m);
        for (auto& digit:da) digit=engine();
        for (auto& digit:db) digit=engine();
        if (repetition%13u==0u) std::fill(da.begin(),da.end(),max);
        if (repetition%17u==0u) std::fill(db.begin(),db.end(),max);
        if (repetition%23u==0u) std::fill(da.begin(),da.end(),0u);
        if (repetition%29u==0u) {
            std::fill(db.begin(),db.end(),0u);
            db.back()=1u;
        }
        if (repetition%43u==0u) {
            std::fill(db.begin(),db.end(),0u);
            db.back()=std::uint64_t{1}<<63u;
        }
        if (repetition%47u==0u) db.back()=1u;
        if (repetition%53u==0u && n==m) da=db;
        const BigUnsigned dividend(da),divisor(db);

        std::ifstream& golden=b140_general[repetition/1000u];
        unsigned ref_repetition=0u;
        std::size_t ref_n=0u,ref_m=0u;
        std::string dividend_hex,divisor_hex,quotient_hex,remainder_hex;
        if (!(golden >> ref_repetition >> ref_n >> ref_m
              >> dividend_hex >> divisor_hex >> quotient_hex >> remainder_hex) ||
            ref_repetition!=repetition || ref_n!=n || ref_m!=m)
            throw std::runtime_error(
                "Missing or unordered B140 general division TSV row");
        check_against_precomputed_hex(dividend,dividend_hex);
        check_against_precomputed_hex(divisor,divisor_hex);

        const auto [actual_q,actual_r]=divrem(dividend,divisor);
        check_against_precomputed_hex(actual_q,quotient_hex);
        check_against_precomputed_hex(actual_r,remainder_hex);
        check_against_precomputed_hex(quot(dividend,divisor),quotient_hex);
        check_against_precomputed_hex(rem(dividend,divisor),remainder_hex);
        ARIADNE_TEST_ASSERT(compare(quot(dividend,divisor),actual_q)==0);
        ARIADNE_TEST_ASSERT(compare(rem(dividend,divisor),actual_r)==0);
        ARIADNE_TEST_ASSERT(compare(actual_r,divisor)<0);
        ARIADNE_TEST_ASSERT(
            compare(add(mul(actual_q,divisor),actual_r),dividend)==0);

        const auto needed=BigUnsigned::divrem_scratch_limbs(
            dividend.size(),divisor.size());
        std::vector<std::uint64_t> qout(dividend.size()+2u,max);
        std::vector<std::uint64_t> rout(
            std::min(dividend.size(),divisor.size())+2u,max);
        std::vector<std::uint64_t> guards(needed+2u,max);
        const auto scratch=BigUnsigned::MutableLimbView(
            guards.data()+1u,needed);
        BigUnsigned::divrem_into(
            dividend.limbs(),divisor.limbs(),qout,rout,scratch);
        check_against_precomputed_hex(BigUnsigned(qout),quotient_hex);
        check_against_precomputed_hex(BigUnsigned(rout),remainder_hex);
        for (std::size_t i=dividend.size();i<qout.size();++i)
            ARIADNE_TEST_ASSERT(qout[i]==0u);
        for (std::size_t i=std::min(dividend.size(),divisor.size());
                i<rout.size();++i)
            ARIADNE_TEST_ASSERT(rout[i]==0u);
        ARIADNE_TEST_ASSERT(guards.front()==max);
        ARIADNE_TEST_ASSERT(guards.back()==max);
        check_against_precomputed_hex(dividend,dividend_hex);
        check_against_precomputed_hex(divisor,divisor_hex);
    }
    for (std::size_t part=0u;part<b140_general.size();++part) {
        std::string extra;
        if (b140_general[part] >> extra)
            throw std::runtime_error("Extra B140 general division TSV rows");
    }
// B57: all 3 output modes of normalized multiprecision
    // division by 2^(64*k+s), including high/low limb boundaries
    // and dividend shorter than divisor. Verify each output against
    // independently certified exact Euclidean TSV records.
    std::array<std::ifstream,3> b57_power_golden{
        std::ifstream(ARIADNE_BIG_UNSIGNED_B57_POWER_SHORT_PATH),
        std::ifstream(ARIADNE_BIG_UNSIGNED_B57_POWER_MEDIUM_PATH),
        std::ifstream(ARIADNE_BIG_UNSIGNED_B57_POWER_LARGE_PATH)};
    for (auto& file:b57_power_golden) {
        if (!file) throw std::runtime_error("Missing B119 B57 power TSV");
        std::string header;
        while (file.peek()=='#') std::getline(file,header);
    }
    std::mt19937_64 b57_power_engine(0xB570B172ULL);
    for (const std::size_t nlimbs:
            {1u,2u,8u,16u,32u,64u,65u,96u,128u}) {
        for (const std::size_t mlimbs:
                {1u,2u,3u,8u,16u,32u,64u,65u}) {
            for (const unsigned bit:{0u,1u,17u,32u,63u}) {
                for (unsigned mode=0u;mode<6u;++mode) {
                    std::vector<std::uint64_t> left(nlimbs);
                    std::vector<std::uint64_t> right(mlimbs,0u);
                    right.back()=std::uint64_t{1}<<bit;
                    for (auto& limb:left) limb=b57_power_engine();
                    if (mode==0u) {
                        std::fill(left.begin(),left.end(),max);
                    } else if (mode==1u) {
                        std::fill(left.begin(),left.end(),0u);
                    } else if (mode==2u) {
                        std::fill(left.begin(),left.end(),0u);
                        left.front()=1u;
                    } else if (mode==3u) {
                        std::fill(left.begin(),left.end(),0u);
                        left.back()=1u;
                    } else if (mode==4u) {
                        for (std::size_t i=0u;i<nlimbs;i+=2u)
                            left[i]=0u;
                    }
                    const BigUnsigned dividend(left),divisor(right);
                    const std::size_t group=nlimbs<=32u?0u:
                        nlimbs<=65u?1u:2u;
                    std::size_t ref_n=0u,ref_m=0u;
                    unsigned ref_bit=0u,ref_mode=0u;
                    std::string dividend_hex,divisor_hex,quotient_hex,remainder_hex;
                    if (!(b57_power_golden[group] >>
                          ref_n >> ref_m >> ref_bit >> ref_mode >>
                          dividend_hex >> divisor_hex >> quotient_hex >> remainder_hex)
                        || ref_n!=nlimbs || ref_m!=mlimbs ||
                           ref_bit!=bit || ref_mode!=mode)
                        throw std::runtime_error("Missing or unordered B119 B57 power row");
                    check_against_precomputed_hex(dividend,dividend_hex);
                    check_against_precomputed_hex(divisor,divisor_hex);
                    const auto [both_q,both_r]=divrem(dividend,divisor);
                    check_against_precomputed_hex(both_q,quotient_hex);
                    check_against_precomputed_hex(both_r,remainder_hex);
                    check_against_precomputed_hex(quot(dividend,divisor),quotient_hex);
                    check_against_precomputed_hex(rem(dividend,divisor),remainder_hex);
                    ARIADNE_TEST_ASSERT(compare(both_r,divisor)<0);
                    check_against_precomputed_hex(dividend,dividend_hex);
                    const std::size_t needed=
                        BigUnsigned::divrem_scratch_limbs(
                            dividend.size(),divisor.size());
                    std::vector<std::uint64_t> qwords(
                        dividend.size()+2u,max);
                    std::vector<std::uint64_t> rwords(
                        std::min(dividend.size(),divisor.size())+2u,max);
                    std::vector<std::uint64_t> storage(
                        needed+2u,max);
                    BigUnsigned::divrem_into(
                        dividend.limbs(),divisor.limbs(),
                        qwords,rwords,BigUnsigned::MutableLimbView(
                            storage.data()+1u,needed));
                    check_against_precomputed_hex(BigUnsigned(qwords),quotient_hex);
                    check_against_precomputed_hex(BigUnsigned(rwords),remainder_hex);
                    ARIADNE_TEST_ASSERT(
                        qwords[dividend.size()]==0u &&
                        qwords[dividend.size()+1u]==0u);
                    const std::size_t m=
                        std::min(dividend.size(),divisor.size());
                    ARIADNE_TEST_ASSERT(rwords[m]==0u &&
                                        rwords[m+1u]==0u);
                    ARIADNE_TEST_ASSERT(storage.front()==max);
                    ARIADNE_TEST_ASSERT(storage.back()==max);
                }
            }
        }
    }

    for (auto& file:b57_power_golden) {
        std::string extra;
        if (file >> extra) throw std::runtime_error("Extra B119 B57 power rows");
    }

    // B55: independently frozen numerical coverage of the native x86-64
    // MULQ/ADCQ digit-product subtraction selected at divisor
    // length 8. Preserve the 7-limb C++ control, 8/9 boundary,
    // 16/32/64-limb denominators, and both divisor operand orders.
    // Full/joint and single-output modes must match exact TSV results
    // even after one- or two-digit corrections.
    std::array<std::ifstream,3> b55_golden{
        std::ifstream(ARIADNE_BIG_UNSIGNED_B55_SHORT_PATH),
        std::ifstream(ARIADNE_BIG_UNSIGNED_B55_MEDIUM_PATH),
        std::ifstream(ARIADNE_BIG_UNSIGNED_B55_LARGE_PATH)};
    for (auto& file:b55_golden) {
        if (!file) throw std::runtime_error("Missing B55 exact TSV");
        std::string header;
        while (file.peek()=='#') std::getline(file,header);
    }
    std::mt19937_64 b55_engine(0xB550DEC1ULL);
    for (const std::size_t m:{7u,8u,9u,16u,32u,64u}) {
        for (const std::size_t extra:{2u,9u,32u}) {
            const std::size_t n=m+extra;
            for (unsigned mode=0u;mode<24u;++mode) {
                std::vector<std::uint64_t> left(n),right(m);
                for (auto& word:left) word=b55_engine();
                for (auto& word:right) word=b55_engine();
                switch (mode%8u) {
                case 0u:
                    std::fill(left.begin(),left.end(),max);
                    std::fill(right.begin(),right.end(),max);
                    break;
                case 1u:
                    std::fill(left.begin(),left.end(),0u);
                    left.back()=max;
                    std::fill(right.begin(),right.end(),0u);
                    right.back()=max;
                    break;
                case 2u:
                    std::fill(left.begin(),left.end(),0u);
                    left.front()=max;
                    left.back()=max;
                    std::fill(right.begin(),right.end(),0u);
                    right.front()=3u;
                    right.back()=max;
                    break;
                case 3u:
                    for (std::size_t i=0u;i<n;i+=2u) left[i]=0u;
                    for (std::size_t i=0u;i<m;i+=2u) right[i]=0u;
                    break;
                case 4u:
                    for (std::size_t i=1u;i+1u<n;++i) left[i]=0u;
                    for (std::size_t i=1u;i+1u<m;++i) right[i]=0u;
                    break;
                case 5u:
                    std::fill_n(left.begin(),n/2u,0u);
                    std::fill_n(right.begin(),m/2u,0u);
                    break;
                case 6u:
                    for (std::size_t i=1u;i<n;i+=3u) left[i]=max;
                    for (std::size_t i=1u;i<m;i+=3u) right[i]=max;
                    break;
                default: break;
                }
                left.back()|=std::uint64_t{1}<<63u;
                right.back()|=std::uint64_t{1}<<63u;
                const BigUnsigned dividend(left),divisor(right);
                const std::size_t group=m<=9u?0u:m<=32u?1u:2u;
                std::size_t ref_m=0u,ref_extra=0u;
                unsigned ref_mode=0u;
                std::string lhs_hex,rhs_hex,quotient_hex,remainder_hex;
                if (!(b55_golden[group] >> ref_m >> ref_extra >> ref_mode
                       >> lhs_hex >> rhs_hex >> quotient_hex >> remainder_hex)
                    || ref_m!=m || ref_extra!=extra || ref_mode!=mode)
                    throw std::runtime_error("Missing or unordered B55 exact TSV row");
                check_against_precomputed_hex(dividend,lhs_hex);
                check_against_precomputed_hex(divisor,rhs_hex);
                const auto [actual_q,actual_r]=divrem(dividend,divisor);
                check_against_precomputed_hex(actual_q,quotient_hex);
                check_against_precomputed_hex(actual_r,remainder_hex);
                check_against_precomputed_hex(quot(dividend,divisor),quotient_hex);
                check_against_precomputed_hex(rem(dividend,divisor),remainder_hex);
                check_against_precomputed_hex(dividend,lhs_hex);
                check_against_precomputed_hex(divisor,rhs_hex);
                const auto needed=BigUnsigned::divrem_scratch_limbs(
                    dividend.size(),divisor.size());
                std::vector<std::uint64_t> qwords(
                    dividend.size()+2u,max);
                std::vector<std::uint64_t> rwords(
                    divisor.size()+2u,max);
                std::vector<std::uint64_t> workspace(
                    needed+2u,max);
                BigUnsigned::divrem_into(
                    dividend.limbs(),divisor.limbs(),qwords,rwords,
                    BigUnsigned::MutableLimbView(
                        workspace.data()+1u,needed));
                check_against_precomputed_hex(BigUnsigned(qwords),quotient_hex);
                check_against_precomputed_hex(BigUnsigned(rwords),remainder_hex);
                ARIADNE_TEST_ASSERT(qwords[dividend.size()]==0u);
                ARIADNE_TEST_ASSERT(qwords[dividend.size()+1u]==0u);
                ARIADNE_TEST_ASSERT(rwords[divisor.size()]==0u);
                ARIADNE_TEST_ASSERT(rwords[divisor.size()+1u]==0u);
                ARIADNE_TEST_ASSERT(workspace.front()==max);
                ARIADNE_TEST_ASSERT(workspace.back()==max);
            }
        }
    }

    for (auto& file:b55_golden) {
        std::string extra;
        if (file >> extra) throw std::runtime_error("Extra B55 exact TSV rows");
    }

    // B27: exercise the fused subtract/add-back implementation,
    // including quotient estimates that need TWO corrections.
    // These fixtures were found by the independent GMP stress probe,
    // not by adjusting the expected result to the native algorithm.
    std::array<std::ifstream,2> b27_golden{
        std::ifstream(ARIADNE_BIG_UNSIGNED_B27_FIRST_PATH),
        std::ifstream(ARIADNE_BIG_UNSIGNED_B27_SECOND_PATH)};
    for (auto& file:b27_golden) {
        if (!file) throw std::runtime_error("Missing B27 independent TSV");
        std::string header;
        while (file.peek()=='#') std::getline(file,header);
    }
    const auto check_b27=[&](char kind,unsigned index,
            BigUnsigned const& numerator,BigUnsigned const& denominator) {
        const std::size_t group=kind=='F'||index<1024u?0u:1u;
        std::string ref_kind,a_hex,b_hex,q_hex,r_hex;
        unsigned ref_index=0u;
        if (!(b27_golden[group] >> ref_kind >> ref_index >>
              a_hex >> b_hex >> q_hex >> r_hex) ||
            ref_kind!=std::string(1,kind) || ref_index!=index)
            throw std::runtime_error("Missing or unordered B27 oracle row");
        check_against_precomputed_hex(numerator,a_hex);
        check_against_precomputed_hex(denominator,b_hex);
        const auto [qvalue,rvalue]=divrem(numerator,denominator);
        check_against_precomputed_hex(qvalue,q_hex);
        check_against_precomputed_hex(rvalue,r_hex);
        check_against_precomputed_hex(quot(numerator,denominator),q_hex);
        check_against_precomputed_hex(rem(numerator,denominator),r_hex);
        const auto needed=BigUnsigned::divrem_scratch_limbs(
            numerator.size(),denominator.size());
        if (denominator.size()>=2u &&
                numerator.size()>=denominator.size())
            ARIADNE_TEST_ASSERT(needed==
                numerator.size()+denominator.size()+1u);
        ARIADNE_TEST_ASSERT(compare(rvalue,denominator)<0);
        check_against_precomputed_hex(numerator,a_hex);
        check_against_precomputed_hex(denominator,b_hex);
    };
    check_b27('F',0u,
        BigUnsigned(std::array<std::uint64_t,4>{
            max,max,max,max}),
        BigUnsigned(std::array<std::uint64_t,2>{
            0xe04e3de411f558a7ULL,0xba43f747b13ac11aULL}));
    check_b27('F',1u,
        BigUnsigned(std::array<std::uint64_t,5>{
            0x65a601688b16b11bULL,0x5b7d0bd5828e0a0cULL,
            0x7ad9be6c73e93fdaULL,0xf370dcf1b5187984ULL,
            0x48293917315898ccULL}),
        BigUnsigned(std::array<std::uint64_t,2>{
            0xc28ba8299ce45291ULL,0x91b35c306fef072dULL}));

    std::mt19937_64 fused_engine(0xB27071ULL);
    for (unsigned repetition=0u;repetition<2048u;++repetition) {
        const std::size_t n=repetition%31u==0u?96u:
                            2u+fused_engine()%32u;
        const std::size_t m=repetition%17u==0u?2u:
                            2u+fused_engine()%n;
        std::vector<std::uint64_t> left(n),right(m);
        for (auto& word:left) word=fused_engine();
        for (auto& word:right) word=fused_engine();
        if (repetition%13u==0u) std::fill(left.begin(),left.end(),max);
        if (repetition%19u==0u) {
            std::fill(right.begin(),right.end(),max);
        }
        if (repetition%23u==0u) {
            std::fill(right.begin(),right.end(),0u);
            right.back()=1u;
        }
        if (repetition%29u==0u) {
            std::fill(right.begin(),right.end(),0u);
            right.back()=std::uint64_t{1}<<63u;
        }
        check_b27('R',repetition,BigUnsigned(left),BigUnsigned(right));
    }

    for (auto& file:b27_golden) {
        std::string extra;
        if (file >> extra) throw std::runtime_error("Extra B27 oracle rows");
    }

    // B48: preserve quotient-only, remainder-only and joint-output
    // semantics for the small divisors used in the unchanged
    // 32-case official benchmark; the independent TSV is the oracle.
    std::ifstream b48_golden(ARIADNE_BIG_UNSIGNED_B48_MULTIPRECISION_PATH);
    if (!b48_golden) throw std::runtime_error("Missing B122 B48 TSV");
    std::string b48_header;
    while (b48_golden.peek()=='#') std::getline(b48_golden,b48_header);
    std::mt19937_64 small_divisor_engine(0xB48A11D1ULL);
    for (const std::size_t count:
            {2u,4u,5u,8u,16u,32u,64u,96u}) {
        for (const std::uint64_t divisor:
                {1ULL,3ULL,5ULL,7ULL,9ULL,16ULL,31ULL,
                 34ULL,63ULL,255ULL,65535ULL,4294967295ULL,
                 4294967296ULL,18446744073709551615ULL}) {
            for (unsigned mode=0u;mode<3u;++mode) {
                std::vector<std::uint64_t> words(count);
                for (auto& limb:words) limb=small_divisor_engine();
                if (mode==0u) {
                    std::fill(words.begin(),words.end(),max);
                } else if (mode==1u) {
                    std::fill(words.begin(),words.end(),0u);
                    words.back()=1u;
                }
                std::size_t ref_count=0u;
                std::uint64_t ref_divisor=0u;
                unsigned ref_mode=0u;
                std::string input_hex,quotient_hex,remainder_hex;
                if (!(b48_golden >> ref_count >> ref_divisor >> ref_mode
                      >> input_hex >> quotient_hex >> remainder_hex)
                    || ref_count!=count || ref_divisor!=divisor || ref_mode!=mode)
                    throw std::runtime_error("Missing or unordered B122 B48 row");
                const BigUnsigned numerator(words),denominator(divisor);
                check_against_precomputed_hex(numerator,input_hex);
                const auto [qvalue,rvalue]=divrem(numerator,denominator);
                check_against_precomputed_hex(qvalue,quotient_hex);
                check_against_precomputed_hex(rvalue,remainder_hex);
                check_against_precomputed_hex(quot(numerator,denominator),quotient_hex);
                check_against_precomputed_hex(rem(numerator,denominator),remainder_hex);
                const auto needed=BigUnsigned::divrem_scratch_limbs(
                    numerator.size(),denominator.size());
                if (denominator.size()>=2u &&
                    numerator.size()>=denominator.size())
                    ARIADNE_TEST_ASSERT(needed==
                        numerator.size()+denominator.size()+1u);
                ARIADNE_TEST_ASSERT(compare(rvalue,denominator)<0);
                check_against_precomputed_hex(numerator,input_hex);
            }
        }
    }

    std::string b48_extra;
    if (b48_golden >> b48_extra)
        throw std::runtime_error("Extra B122 B48 rows");

    // B141/B28: equal-high fast-path cases use two sequential exact TSVs.
    std::array<std::ifstream,2u> b141_b28;
    for (std::size_t part=0u;part<b141_b28.size();++part) {
        const std::string path=std::string(ARIADNE_TEST_DATA_DIR)+
            "/big_unsigned_b141_b28_"+std::to_string(part)+".tsv";
        b141_b28[part].open(path);
        if (!b141_b28[part])
            throw std::runtime_error("Missing B141 B28 TSV: "+path);
        std::string header;
        while (b141_b28[part].peek()=='#')
            std::getline(b141_b28[part],header);
    }
    std::mt19937_64 equal_high_engine(0xB280C27ULL);
    for (unsigned repetition=0u;repetition<2048u;++repetition) {
        const std::size_t count=repetition%37u==0u?128u:
                                2u+equal_high_engine()%64u;
        std::vector<std::uint64_t> left(count),right(count);
        for (auto& word:left) word=equal_high_engine();
        for (auto& word:right) word=equal_high_engine();
        right.back()|=1u;
        left.back()=right.back();
        if (repetition%4u==0u) {
            left=right;
        } else if (repetition%4u==1u) {
            left=right;
            ++left.front();
        } else if (repetition%4u==2u) {
            right=left;
            ++right.front();
        }
        const BigUnsigned numerator(left),denominator(right);
        std::ifstream& golden=b141_b28[repetition/1024u];
        unsigned ref_repetition=0u;
        std::size_t ref_count=0u;
        std::string numerator_hex,denominator_hex,quotient_hex,remainder_hex;
        if (!(golden >> ref_repetition >> ref_count
              >> numerator_hex >> denominator_hex
              >> quotient_hex >> remainder_hex) ||
            ref_repetition!=repetition || ref_count!=count)
            throw std::runtime_error(
                "Missing or unordered B141 B28 TSV row");
        check_against_precomputed_hex(numerator,numerator_hex);
        check_against_precomputed_hex(denominator,denominator_hex);
        const auto [native_q,native_r]=divrem(numerator,denominator);
        check_against_precomputed_hex(native_q,quotient_hex);
        check_against_precomputed_hex(native_r,remainder_hex);
        check_against_precomputed_hex(
            quot(numerator,denominator),quotient_hex);
        check_against_precomputed_hex(
            rem(numerator,denominator),remainder_hex);
        ARIADNE_TEST_ASSERT(compare(native_r,denominator)<0);

        // The public joint-output contract still requires the same
        // scratch size and clears surplus output limbs. This fast path
        // must not write into that scratch (including its guard limbs).
        const auto needed=BigUnsigned::divrem_scratch_limbs(
            numerator.size(),denominator.size());
        std::vector<std::uint64_t> qout(numerator.size()+2u,max);
        std::vector<std::uint64_t> rout(denominator.size()+2u,max);
        std::vector<std::uint64_t> guards(needed+2u,max);
        BigUnsigned::divrem_into(
            numerator.limbs(),denominator.limbs(),qout,rout,
            BigUnsigned::MutableLimbView(guards.data()+1u,needed));
        check_against_precomputed_hex(BigUnsigned(qout),quotient_hex);
        check_against_precomputed_hex(BigUnsigned(rout),remainder_hex);
        for (const auto word:guards) {
            ARIADNE_TEST_ASSERT(word==max);
        }
        for (std::size_t i=numerator.size();i<qout.size();++i) {
            ARIADNE_TEST_ASSERT(qout[i]==0u);
        }
        for (std::size_t i=denominator.size();i<rout.size();++i) {
            ARIADNE_TEST_ASSERT(rout[i]==0u);
        }
    }

    for (std::size_t part=0u;part<b141_b28.size();++part) {
        std::string extra;
        if (b141_b28[part] >> extra)
            throw std::runtime_error("Extra B141 B28 TSV rows");
    }

    // B142/B36: classify equal-size quotient paths against two exact TSVs.
    std::array<std::ifstream,2u> b142_b36;
    for (std::size_t part=0u;part<b142_b36.size();++part) {
        const std::string path=std::string(ARIADNE_TEST_DATA_DIR)+
            "/big_unsigned_b142_b36_"+std::to_string(part)+".tsv";
        b142_b36[part].open(path);
        if (!b142_b36[part])
            throw std::runtime_error("Missing B142 B36 TSV: "+path);
        std::string header;
        while (b142_b36[part].peek()=='#')
            std::getline(b142_b36[part],header);
    }
    std::mt19937_64 high_ratio_engine(0xB360A11ULL);
    std::size_t b142_ordinal=0u;
    for (const std::size_t count:
            {2u,3u,4u,7u,8u,16u,32u,64u,96u,128u}) {
        for (unsigned repetition=0u;repetition<128u;++repetition) {
            std::vector<std::uint64_t> left(count),right(count);
            for (auto& limb:left) limb=high_ratio_engine();
            for (auto& limb:right) limb=high_ratio_engine();
            // These pairs are canonical; low words can be saturated,
            // sparse or identical regardless of the chosen high ratio.
            if (repetition%11u==0u) {
                std::fill(left.begin(),left.end(),max);
                std::fill(right.begin(),right.end(),max);
            } else if (repetition%11u==1u) {
                std::fill(left.begin(),left.end(),0u);
                std::fill(right.begin(),right.end(),0u);
            }
            const std::uint64_t high=std::uint64_t{1}<<63u;
            switch (repetition%8u) {
            case 0u: // q=0: strictly smaller top limb.
                left.back()=high-1u; right.back()=high;
                break;
            case 1u: // q=1: strictly greater top, below 2*bh.
                left.back()=high+123u; right.back()=high;
                break;
            case 2u: // q=0 or 1: equal top, B28 path.
                left.back()=high; right.back()=high;
                break;
            case 3u: // q>=2: outside the new fast path.
                left.back()=high; right.back()=2u;
                break;
            case 4u: // q=1 at the largest possible high limb.
                left.back()=max; right.back()=high;
                break;
            case 5u: // boundary ah==2*bh, must use long division.
                left.back()=high; right.back()=high/2u;
                break;
            case 6u: // q=0: top of divisor is saturated.
                left.back()=high; right.back()=max;
                break;
            default: // arbitrary canonical high limbs.
                left.back()|=1u; right.back()|=1u;
                break;
            }
            const BigUnsigned numerator(left),denominator(right);
            std::ifstream& golden=b142_b36[b142_ordinal/640u];
            std::size_t ref_ordinal=0u,ref_count=0u;
            unsigned ref_repetition=0u;
            std::string numerator_hex,denominator_hex,quotient_hex,remainder_hex;
            if (!(golden >> ref_ordinal >> ref_count >> ref_repetition
                  >> numerator_hex >> denominator_hex
                  >> quotient_hex >> remainder_hex) ||
                ref_ordinal!=b142_ordinal || ref_count!=count ||
                ref_repetition!=repetition)
                throw std::runtime_error(
                    "Missing or unordered B142 B36 TSV row");
            ++b142_ordinal;
            check_against_precomputed_hex(numerator,numerator_hex);
            check_against_precomputed_hex(denominator,denominator_hex);
            const auto [native_q,native_r]=divrem(numerator,denominator);
            check_against_precomputed_hex(native_q,quotient_hex);
            check_against_precomputed_hex(native_r,remainder_hex);
            check_against_precomputed_hex(
                quot(numerator,denominator),quotient_hex);
            check_against_precomputed_hex(
                rem(numerator,denominator),remainder_hex);
            ARIADNE_TEST_ASSERT(compare(native_r,denominator)<0);
        }
    }

    if (b142_ordinal!=1280u)
        throw std::runtime_error("Wrong B142 B36 row count");
    for (std::size_t part=0u;part<b142_b36.size();++part) {
        std::string extra;
        if (b142_b36[part] >> extra)
            throw std::runtime_error("Extra B142 B36 TSV rows");
    }

    // B143/B29: deterministic sparse/dense divisor shapes use exact TSVs.
    std::array<std::ifstream,2u> b143_b29;
    for (std::size_t part=0u;part<b143_b29.size();++part) {
        const std::string path=std::string(ARIADNE_TEST_DATA_DIR)+
            "/big_unsigned_b143_b29_"+std::to_string(part)+".tsv";
        b143_b29[part].open(path);
        if (!b143_b29[part])
            throw std::runtime_error("Missing B143 B29 TSV: "+path);
        std::string header;
        while (b143_b29[part].peek()=='#')
            std::getline(b143_b29[part],header);
    }
    std::mt19937_64 sparse_engine(0xB290D1F5ULL);
    for (unsigned repetition=0u;repetition<1024u;++repetition) {
        const std::size_t n=16u+sparse_engine()%113u;
        const std::size_t m=8u+sparse_engine()%(n-7u);
        std::vector<std::uint64_t> left(n),right(m,0u);
        for (auto& word:left) word=sparse_engine();
        right.front()=sparse_engine();
        right.back()=std::uint64_t{1}<<(repetition%64u);
        if (repetition%3u==0u) {
            right[m/2u]=sparse_engine();
        }
        if (repetition%5u==0u) {
            // Partially dense divisors must retain the original kernel.
            for (std::size_t i=1u;i+1u<m;++i) {
                if (i%3u==0u) right[i]=sparse_engine();
            }
        }
        if (repetition%7u==0u) {
            std::fill(left.begin(),left.end(),max);
        }
        if (repetition%11u==0u) {
            left.back()=1u;
        }
        const BigUnsigned numerator(left),denominator(right);
        std::ifstream& golden=b143_b29[repetition/512u];
        unsigned ref_repetition=0u;
        std::size_t ref_n=0u,ref_m=0u;
        std::string numerator_hex,denominator_hex,quotient_hex,remainder_hex;
        if (!(golden >> ref_repetition >> ref_n >> ref_m
              >> numerator_hex >> denominator_hex
              >> quotient_hex >> remainder_hex) ||
            ref_repetition!=repetition || ref_n!=n || ref_m!=m)
            throw std::runtime_error(
                "Missing or unordered B143 B29 TSV row");
        check_against_precomputed_hex(numerator,numerator_hex);
        check_against_precomputed_hex(denominator,denominator_hex);
        const auto [native_q,native_r]=divrem(numerator,denominator);
        check_against_precomputed_hex(native_q,quotient_hex);
        check_against_precomputed_hex(native_r,remainder_hex);
        check_against_precomputed_hex(
            quot(numerator,denominator),quotient_hex);
        check_against_precomputed_hex(
            rem(numerator,denominator),remainder_hex);
        ARIADNE_TEST_ASSERT(compare(native_r,denominator)<0);
    }

    for (std::size_t part=0u;part<b143_b29.size();++part) {
        std::string extra;
        if (b143_b29[part] >> extra)
            throw std::runtime_error("Extra B143 B29 TSV rows");
    }

    // B144/B45: preserve all zero-low-limb and decimal-power cases
    // with one sequential source-keyed exact fixture.
    std::ifstream b144_b45(
        std::string(ARIADNE_TEST_DATA_DIR)+"/big_unsigned_b144_b45.tsv");
    if (!b144_b45) throw std::runtime_error("Missing B144 B45 TSV");
    std::string b144_header;
    while (b144_b45.peek()=='#') std::getline(b144_b45,b144_header);
    const auto check_a33=[&](char kind,std::size_t key1,std::size_t key2,
            std::uint64_t key3,BigUnsigned const& numerator,
            BigUnsigned const& denominator) {
        std::string ref_kind,numerator_hex,denominator_hex,quotient_hex,remainder_hex;
        std::size_t ref_key1=0u,ref_key2=0u;
        std::uint64_t ref_key3=0u;
        if (!(b144_b45 >> ref_kind >> ref_key1 >> ref_key2 >> ref_key3
              >> numerator_hex >> denominator_hex >> quotient_hex >> remainder_hex) ||
            ref_kind!=std::string(1,kind) || ref_key1!=key1 ||
            ref_key2!=key2 || ref_key3!=key3)
            throw std::runtime_error("Missing or unordered B144 B45 TSV row");
        check_against_precomputed_hex(numerator,numerator_hex);
        check_against_precomputed_hex(denominator,denominator_hex);
        const auto [native_q,native_r]=divrem(numerator,denominator);
        check_against_precomputed_hex(native_q,quotient_hex);
        check_against_precomputed_hex(native_r,remainder_hex);
        check_against_precomputed_hex(quot(numerator,denominator),quotient_hex);
        check_against_precomputed_hex(rem(numerator,denominator),remainder_hex);
        ARIADNE_TEST_ASSERT(compare(native_r,denominator)<0);
        ARIADNE_TEST_ASSERT(compare(
            add(mul(native_q,denominator),native_r),numerator)==0);
        return std::pair<std::string,std::string>{quotient_hex,remainder_hex};
    };

    std::mt19937_64 shifted_engine(0xB450D1D1ULL);
    for (const std::size_t m:{2u,3u,5u,8u,16u,32u,64u}) {
        for (const std::size_t zeros:
                {std::size_t{1},m/2u,m-1u}) {
            for (unsigned repetition=0u;repetition<24u;++repetition) {
                const std::size_t n=m+
                    (repetition%4u==0u?0u:
                     repetition%4u==1u?1u:
                     repetition%4u==2u?m:2u*m);
                std::vector<std::uint64_t> left(n),right(m);
                for (auto& word:left) word=shifted_engine();
                for (auto& word:right) word=shifted_engine();
                std::fill_n(right.begin(),zeros,0u);
                right.back()|=std::uint64_t{1}<<63u;
                left.back()|=std::uint64_t{1}<<63u;
                if (repetition%6u==0u && n==m) {
                    left=right;
                } else if (repetition%6u==1u && n==m) {
                    left.back()=right.back()-1u;
                } else if (repetition%6u==2u) {
                    std::fill_n(left.begin(),zeros,0u);
                } else if (repetition%6u==3u) {
                    std::fill(left.begin(),left.end(),max);
                }
                const BigUnsigned numerator(left),denominator(right);
                const auto expected=check_a33(
                    'S',m,zeros,repetition,numerator,denominator);
                const auto needed=BigUnsigned::divrem_scratch_limbs(
                    numerator.size(),denominator.size());
                std::vector<std::uint64_t> qout(
                    numerator.size()+2u,max);
                std::vector<std::uint64_t> rout(
                    std::min(numerator.size(),denominator.size())+2u,max);
                std::vector<std::uint64_t> guards(needed+2u,max);
                BigUnsigned::divrem_into(
                    numerator.limbs(),denominator.limbs(),qout,rout,
                    BigUnsigned::MutableLimbView(
                        guards.data()+1u,needed));
                check_against_precomputed_hex(BigUnsigned(qout),expected.first);
                check_against_precomputed_hex(BigUnsigned(rout),expected.second);
                for (std::size_t i=numerator.size();i<qout.size();++i)
                    ARIADNE_TEST_ASSERT(qout[i]==0u);
                const std::size_t rem_count=
                    std::min(numerator.size(),denominator.size());
                for (std::size_t i=rem_count;i<rout.size();++i)
                    ARIADNE_TEST_ASSERT(rout[i]==0u);
                ARIADNE_TEST_ASSERT(guards.front()==max);
                ARIADNE_TEST_ASSERT(guards.back()==max);
            }
        }
    }

    for (const unsigned digits:{72u,144u,288u,576u}) {
        const BigUnsigned decimal_power=from_decimal(
            std::string("1")+std::string(digits,'0'));
        for (const std::size_t bits:{2048u,4096u}) {
            const BigUnsigned binary_power=shift_left(
                BigUnsigned(std::uint64_t{1}),bits-1u);
            for (const std::uint64_t offset:
                    {0u,1u,13u,127u,1023u}) {
                const BigUnsigned numerator=add(
                    binary_power,BigUnsigned(offset));
                (void)check_a33(
                    'D',digits,bits,offset,numerator,decimal_power);
            }
        }
    }
    std::string b144_extra;
    if (b144_b45 >> b144_extra)
        throw std::runtime_error("Extra B144 B45 TSV rows");

    // B145: force all 64 divisor-normalization shifts and preserve all
    // ten original dividend modes with a deterministic source-keyed TSV.
    std::ifstream b145_normalization(
        std::string(ARIADNE_TEST_DATA_DIR)+
        "/big_unsigned_b145_normalization.tsv");
    if (!b145_normalization)
        throw std::runtime_error("Missing B145 normalization TSV");
    std::string b145_header;
    while (b145_normalization.peek()=='#')
        std::getline(b145_normalization,b145_header);
    for (unsigned shift=0u;shift<64u;++shift) {
        const std::uint64_t top=max>>shift;
        for (unsigned mode=0u;mode<10u;++mode) {
            const std::array<std::uint64_t,2> divisor_words{
                (mode%2u)==0u ? max : std::uint64_t{1},top};
            const std::array<std::uint64_t,4> dividend_words{
                (mode%3u)==0u ? max : std::uint64_t{mode+1u},
                (mode%3u)==1u ? 0u : max,
                (mode%2u)==0u ? max : std::uint64_t{1},
                (mode%3u)==2u ? top-1u : top};
            const BigUnsigned dividend_value(dividend_words);
            const BigUnsigned divisor_value(divisor_words);
            unsigned ref_shift=0u,ref_mode=0u;
            std::string dividend_hex,divisor_hex,quotient_hex,remainder_hex;
            if (!(b145_normalization >> ref_shift >> ref_mode
                  >> dividend_hex >> divisor_hex
                  >> quotient_hex >> remainder_hex) ||
                ref_shift!=shift || ref_mode!=mode)
                throw std::runtime_error(
                    "Missing or unordered B145 normalization TSV row");
            check_against_precomputed_hex(dividend_value,dividend_hex);
            check_against_precomputed_hex(divisor_value,divisor_hex);
            const auto [native_q,native_r]=divrem(
                dividend_value,divisor_value);
            check_against_precomputed_hex(native_q,quotient_hex);
            check_against_precomputed_hex(native_r,remainder_hex);
            check_against_precomputed_hex(
                quot(dividend_value,divisor_value),quotient_hex);
            check_against_precomputed_hex(
                rem(dividend_value,divisor_value),remainder_hex);
            ARIADNE_TEST_ASSERT(compare(native_r,divisor_value)<0);
            ARIADNE_TEST_ASSERT(compare(
                add(mul(native_q,divisor_value),native_r),
                dividend_value)==0);
        }
    }
    std::string b145_extra;
    if (b145_normalization >> b145_extra)
        throw std::runtime_error("Extra B145 normalization TSV rows");

    // Focused zero, short dividend, exact quotient and single-limb
    // divisor cases; keep both divided-output implementations checked.
    for (const auto& numerator: {
            BigUnsigned{},BigUnsigned(std::uint64_t{1}),
            BigUnsigned(std::array<std::uint64_t,2>{0u,1u}),
            BigUnsigned(std::array<std::uint64_t,4>{max,max,max,max})}) {
        for (const auto& denominator: {
                BigUnsigned(std::uint64_t{1}),BigUnsigned(max),
                BigUnsigned(std::array<std::uint64_t,2>{max,max}),
                BigUnsigned(std::array<std::uint64_t,5>{0u,0u,0u,0u,1u})}) {
            const auto both=divrem(numerator,denominator);
            ARIADNE_TEST_ASSERT(compare(
                quot(numerator,denominator),both.first)==0);
            ARIADNE_TEST_ASSERT(compare(
                rem(numerator,denominator),both.second)==0);
        }
    }

    const std::array<std::uint64_t,3> dividend{max,max,max};
    const std::array<std::uint64_t,2> divisor{max,max};
    std::array<std::uint64_t,4> qout{max,max,max,max};
    std::array<std::uint64_t,3> rout{max,max,max};
    const auto needed=BigUnsigned::divrem_scratch_limbs(3u,2u);
    std::vector<std::uint64_t> scratch(needed,max);
    bool zero_threw=false,owning_zero_threw=false;
    bool short_quotient=false,short_remainder=false,short_scratch=false;
    try {
        BigUnsigned::divrem_into(
            dividend,BigUnsigned{}.limbs(),qout,rout,scratch);
    } catch (std::domain_error const&) { zero_threw=true; }
    try {
        (void)divrem(BigUnsigned(dividend),zero);
    } catch (std::domain_error const&) { owning_zero_threw=true; }
    bool quotient_zero_threw=false,remainder_zero_threw=false;
    try { (void)quot(BigUnsigned(dividend),zero); }
    catch (std::domain_error const&) { quotient_zero_threw=true; }
    try { (void)rem(BigUnsigned(dividend),zero); }
    catch (std::domain_error const&) { remainder_zero_threw=true; }
    ARIADNE_TEST_ASSERT(quotient_zero_threw && remainder_zero_threw);
    try {
        BigUnsigned::divrem_into(
            dividend,divisor,std::span(qout).first(2u),rout,scratch);
    } catch (std::length_error const&) { short_quotient=true; }
    try {
        BigUnsigned::divrem_into(
            dividend,divisor,qout,std::span(rout).first(1u),scratch);
    } catch (std::length_error const&) { short_remainder=true; }
    try {
        BigUnsigned::divrem_into(
            dividend,divisor,qout,rout,std::span(scratch).first(needed-1u));
    } catch (std::length_error const&) { short_scratch=true; }
    ARIADNE_TEST_ASSERT(zero_threw && owning_zero_threw &&
                        short_quotient && short_remainder && short_scratch);
    for (auto digit:qout) ARIADNE_TEST_ASSERT(digit==max);
    for (auto digit:rout) ARIADNE_TEST_ASSERT(digit==max);
    for (auto digit:scratch) ARIADNE_TEST_ASSERT(digit==max);
    bool overflow_threw=false;
    try {
        (void)BigUnsigned::divrem_scratch_limbs(
            std::numeric_limits<std::size_t>::max(),2u);
    } catch (std::length_error const&) { overflow_threw=true; }
    ARIADNE_TEST_ASSERT(overflow_threw);
}


void TestBigUnsigned::test_gcd_decimal()
{
    const auto zero=BigUnsigned{};
    const auto one=BigUnsigned(std::uint64_t{1});
    ARIADNE_TEST_ASSERT(gcd(zero,zero).is_zero());
    ARIADNE_TEST_ASSERT(compare(gcd(zero,one),one)==0);
    ARIADNE_TEST_ASSERT(compare(gcd(one,zero),one)==0);
    ARIADNE_TEST_ASSERT(compare(
        gcd(BigUnsigned(std::uint64_t{48}),BigUnsigned(std::uint64_t{18})),
        BigUnsigned(std::uint64_t{6}))==0);
    // R6: original 7 widths x 4 offsets x both operand orders.
    {
        std::ifstream golden(ARIADNE_BIG_UNSIGNED_R6_ORACLE_PATH);
        require_gcd_golden(golden,"R6");
        for (const unsigned bits:{128u,192u,256u,257u,320u,512u,1024u}) {
            for (const std::uint64_t offset:{1u,3u,17u,65537u}) {
                const auto high=shift_left(BigUnsigned(std::uint64_t{1}),bits);
                const auto left=add(high,BigUnsigned(offset*3u));
                const auto right=add(high,BigUnsigned(offset*5u));
                for (unsigned order=0u;order<2u;++order) {
                    const auto& a=order==0u ? left : right;
                    const auto& b=order==0u ? right : left;
                    check_gcd_golden(golden,{std::to_string(bits),
                        std::to_string(offset),std::to_string(order)},a,b);
                }
            }
        }
        ensure_gcd_golden_end(golden,"R6");
    }
    // R8: original conditional PRNG draws and all wide/small cases.
    {
        std::ifstream golden(ARIADNE_BIG_UNSIGNED_R8_ORACLE_PATH);
        require_gcd_golden(golden,"R8");
        std::mt19937_64 random(0x8890ABCDEF234ULL);
        for (unsigned i=0u;i<640u;++i) {
            const std::uint64_t x=i%11u==0u ? 0u : random();
            const std::uint64_t y=i%17u==0u ? 0u : random();
            const BigUnsigned a(x),b(y);
            check_gcd_golden(golden,{"S",std::to_string(i),"0"},a,b);
        }
        for (unsigned bits:{128u,255u,256u,257u,1024u,4096u}) {
            const BigUnsigned wide=add(
                shift_left(one,bits),BigUnsigned(std::uint64_t{77}));
            for (std::uint64_t d:std::array<std::uint64_t,6u>{
                    1u,2u,3u,17u,65537u,
                    std::numeric_limits<std::uint64_t>::max()}) {
                const BigUnsigned small(d);
                check_gcd_golden(golden,
                    {"W",std::to_string(bits),std::to_string(d)},wide,small);
            }
        }
        ensure_gcd_golden_end(golden,"R8");
    }
    // R10: borrowed operands, both orders, zero and self-aliasing.
    {
        std::ifstream golden(ARIADNE_BIG_UNSIGNED_R10_ORACLE_PATH);
        require_gcd_golden(golden,"R10");
        for (unsigned bits:{53u,128u,255u,256u,257u,
                           512u,1024u,2048u,4096u}) {
            const BigUnsigned base=shift_left(one,bits);
            const BigUnsigned close=add(base,BigUnsigned(std::uint64_t{41}));
            const BigUnsigned apart=add(
                shift_left(one,bits+3u),BigUnsigned(std::uint64_t{93}));
            const BigUnsigned multiple=mul(close,BigUnsigned(std::uint64_t{7}));
            const std::array<std::pair<BigUnsigned,BigUnsigned>,5u> cases={{
                {close,apart},{apart,close},{close,multiple},
                {apart,BigUnsigned(std::uint64_t{65537})},{multiple,zero}
            }};
            for (unsigned index=0u;index<cases.size();++index) {
                check_gcd_golden(golden,{std::to_string(bits),
                    std::to_string(index)},cases[index].first,cases[index].second);
            }
            check_gcd_golden(golden,{std::to_string(bits),"5"},close,close);
        }
        ensure_gcd_golden_end(golden,"R10");
    }
    // R11: original Euclidean quotient=1,2,3 and exact boundary cases.
    {
        std::ifstream golden(ARIADNE_BIG_UNSIGNED_R11_ORACLE_PATH);
        require_gcd_golden(golden,"R11");
        for (unsigned bits:{53u,128u,192u,255u,256u,257u,
                           320u,512u,1024u,2048u,4096u}) {
            const BigUnsigned unit=shift_left(one,bits);
            const BigUnsigned divisor=add(unit,BigUnsigned(std::uint64_t{5}));
            const BigUnsigned twice=mul(unit,BigUnsigned(std::uint64_t{2}));
            const BigUnsigned thrice=mul(unit,BigUnsigned(std::uint64_t{3}));
            const std::array<BigUnsigned,5u> dividends={
                add(twice,BigUnsigned(std::uint64_t{7})),
                add(twice,BigUnsigned(std::uint64_t{10})),
                add(twice,BigUnsigned(std::uint64_t{17})),
                add(thrice,BigUnsigned(std::uint64_t{23})),
                divisor
            };
            for (unsigned index=0u;index<dividends.size();++index) {
                check_gcd_golden(golden,{std::to_string(bits),
                    std::to_string(index)},dividends[index],divisor);
            }
        }
        ensure_gcd_golden_end(golden,"R11");
    }
    // R26: preserve all 9 widths, 40 trials and 9 quotient variants.
    {
        std::mt19937_64 engine(0x52632633474344ULL);
        for (const std::size_t count:{2u,3u,4u,5u,8u,16u,127u,128u,129u}) {
            const std::string path=std::string(ARIADNE_TEST_DATA_DIR)+
                "/big_unsigned_r26_"+std::to_string(count)+".tsv";
            std::ifstream golden(path);
            require_gcd_golden(golden,"R26");
            const BigUnsigned top=shift_left(one,64u*(count-1u));
            for (unsigned trial=0u;trial<40u;++trial) {
                const BigUnsigned low(engine()|1u);
                const BigUnsigned b=add(top,low);
                std::size_t ref_count=0u;
                unsigned ref_trial=0u;
                std::string low_hex,b_hex;
                std::array<std::string,9u> results;
                if (!(golden >> ref_count >> ref_trial >> low_hex >> b_hex))
                    throw std::runtime_error("Missing R26 TSV row");
                for (auto& value:results)
                    if (!(golden >> value)) throw std::runtime_error("Malformed R26 TSV row");
                if (ref_count!=count || ref_trial!=trial)
                    throw std::runtime_error("Unordered R26 TSV row");
                check_against_precomputed_hex(low,low_hex);
                check_against_precomputed_hex(b,b_hex);
                unsigned index=0u;
                for (std::uint64_t q:{2u,3u,4u}) {
                    for (unsigned variant=0u;variant<3u;++variant) {
                        const BigUnsigned product=mul(b,BigUnsigned(q));
                        const BigUnsigned a=variant==0u ? product :
                            variant==1u ? add(product,one) :
                            sub(add(product,b),one);
                        check_against_precomputed_hex(gcd(a,b),results[index]);
                        check_against_precomputed_hex(gcd(b,a),results[index]);
                        check_against_precomputed_hex(b,b_hex);
                        ++index;
                    }
                }
                check_against_precomputed_hex(low,low_hex);
            }
            ensure_gcd_golden_end(golden,"R26");
        }
    }
    // R27: complete original shared-prefix and scalar-tail random streams.
    {
        std::ifstream golden(ARIADNE_BIG_UNSIGNED_R27_ORACLE_PATH);
        require_gcd_golden(golden,"R27");
        std::mt19937_64 random(0x5277705245464958ULL);
        for (std::size_t count:{2u,3u,4u,5u,7u,8u,9u}) {
            for (unsigned repetition=0u;repetition<96u;++repetition) {
                std::vector<std::uint64_t> left(count),right(count);
                for (auto& word:left) word=random();
                left.back()|=std::uint64_t{1}<<63u;
                right=left;
                const std::uint64_t difference=
                    repetition%5u==0u ? 0u :
                    repetition%5u==1u ? 1u :
                    repetition%5u==2u ? 6u :
                    repetition%5u==3u ? 65537u :
                    std::uint64_t{1}<<63u;
                left.front()=std::numeric_limits<std::uint64_t>::max();
                right.front()=left.front()-difference;
                check_gcd_golden(golden,{"P",std::to_string(count),
                    std::to_string(repetition)},BigUnsigned(left),BigUnsigned(right));
            }
        }
        for (std::size_t count:{1u,2u,3u,4u,5u,7u,8u,9u}) {
            for (unsigned repetition=0u;repetition<96u;++repetition) {
                std::vector<std::uint64_t> limbs(count);
                for (auto& word:limbs) word=random();
                limbs.back()|=std::uint64_t{1}<<63u;
                const std::uint64_t d=repetition%6u==0u ? 1u :
                    repetition%6u==1u ? 2u :
                    repetition%6u==2u ? 65537u :
                    repetition%6u==3u ? std::uint64_t{1}<<63u :
                    repetition%6u==4u ?
                        std::numeric_limits<std::uint64_t>::max() :
                        random()|1u;
                check_gcd_golden(golden,{"S",std::to_string(count),
                    std::to_string(repetition)},BigUnsigned(limbs),BigUnsigned(d));
            }
        }
        ensure_gcd_golden_end(golden,"R27");
    }
    // R21: preserve the fixed-factor identity, all 40 dense samples,
    // source immutability and both borrowed-input orders.
    {
        const BigUnsigned factor=add(shift_left(one,64u),BigUnsigned(std::uint64_t{3}));
        std::mt19937_64 random(0x21BCD128ULL);
        std::ifstream golden(ARIADNE_BIG_UNSIGNED_R21_ORACLE_PATH);
        require_gcd_golden(golden,"R21");
        for (std::size_t count:{2u,4u,5u,8u,16u,64u,127u,128u,129u,160u}) {
            const BigUnsigned power=shift_left(one,64u*count-2u);
            const BigUnsigned left=mul(sub(power,one),factor);
            const BigUnsigned right=mul(add(power,one),factor);
            ARIADNE_TEST_ASSERT(compare(gcd(left,right),factor)==0);
            ARIADNE_TEST_ASSERT(compare(gcd(right,left),factor)==0);
            for (unsigned trial=0u;trial<4u;++trial) {
                std::vector<std::uint64_t> x(count),y(count);
                for (auto& word:x) word=random();
                for (auto& word:y) word=random();
                x.back()|=std::uint64_t{1}<<63u;
                y.back()|=std::uint64_t{1}<<63u;
                check_gcd_golden(golden,{std::to_string(count),
                    std::to_string(trial)},BigUnsigned(x),BigUnsigned(y));
            }
        }
        ensure_gcd_golden_end(golden,"R21");
    }
    // 10^1152 has 3827 bits and midpoint 1913. The next split power
    // 10^576 has 1914 bits, so cutoff-32 rejects the final square and
    // must select a smaller power during its later recursive calls.
    const BigUnsigned cutoff32_boundary=pow(BigUnsigned(std::uint64_t{10}),std::uint64_t{1152});
    const std::string cutoff32_decimal="1"+std::string(1152u,'0');
    ARIADNE_TEST_EQUALS(to_decimal_cutoff32_diagnostic(cutoff32_boundary),cutoff32_decimal);
    ARIADNE_TEST_EQUALS(to_decimal_cutoff16_diagnostic(cutoff32_boundary),cutoff32_decimal);
    // The root split uses 10^1152. Its low block 10^1000 is larger
    // than 32 limbs and splits again while preserving all 1152 digits.
    const BigUnsigned padded_low=add(
        pow(BigUnsigned(std::uint64_t{10}),std::uint64_t{2304}),
        pow(BigUnsigned(std::uint64_t{10}),std::uint64_t{1000}));
    const std::string padded_decimal="1"+std::string(1303u,'0')+"1"+std::string(1000u,'0');
    ARIADNE_TEST_EQUALS(to_decimal_cutoff32_diagnostic(padded_low),padded_decimal);
    ARIADNE_TEST_ASSERT(to_decimal(zero)=="0");
    ARIADNE_TEST_ASSERT(to_decimal(one)=="1");
    ARIADNE_TEST_ASSERT(to_decimal(from_decimal("000"))=="0");
    ARIADNE_TEST_ASSERT(to_decimal(from_decimal("00018446744073709551616"))
                        =="18446744073709551616");
    ARIADNE_TEST_ASSERT(to_decimal(from_decimal(
        "340282366920938463463374607431768211455"))
        =="340282366920938463463374607431768211455");

    for (const auto* bad: {"","-1","+1"," 1","1 ","1.2","0x1","12a"}) {
        bool threw=false;
        try { (void)from_decimal(bad); }
        catch (std::invalid_argument const&) { threw=true; }
        ARIADNE_TEST_ASSERT(threw);
    }

    std::mt19937_64 engine(0xB16B00B5ULL);
    // B76: every original seeded GCD/decimal case has an independent,
    // immutable exact TSV oracle, including its original source operands.
    std::ifstream b76_fixture;
    for (unsigned repetition=0u;repetition<2800u;++repetition) {
        if (repetition%200u==0u) {
            if (b76_fixture.is_open()) {
                std::string extra;
                if (b76_fixture >> extra)
                    throw std::runtime_error("Extra B76 GCD/decimal TSV rows");
                b76_fixture.close();
                b76_fixture.clear();
            }
            const std::string path=std::string(ARIADNE_BIG_UNSIGNED_B76_GCD_ORACLE_PREFIX)
                +std::to_string(repetition/200u)+".tsv";
            b76_fixture.open(path);
            if (!b76_fixture)
                throw std::runtime_error("Missing B76 GCD/decimal TSV: "+path);
            std::string header;
            while (b76_fixture.peek()=='#') std::getline(b76_fixture,header);
        }
        const std::size_t n=repetition%113u==0u?96u:engine()%20u;
        const std::size_t m=repetition%71u==0u?32u:engine()%12u;
        std::vector<std::uint64_t> left(n),right(m);
        for (auto& digit:left) digit=engine();
        for (auto& digit:right) digit=engine();
        if (repetition%17u==0u) std::fill(left.begin(),left.end(),0u);
        if (repetition%13u==0u) std::fill(left.begin(),left.end(),
                                           std::numeric_limits<std::uint64_t>::max());
        if (repetition%11u==0u) std::fill(right.begin(),right.end(),0u);
        const BigUnsigned x(left),y(right);

        unsigned source_iteration=0u;
        std::size_t source_n=0u,source_m=0u;
        std::string source_x,source_y,expected_gcd,expected_decimal;
        if (!(b76_fixture >> source_iteration >> source_n >> source_m
              >> source_x >> source_y >> expected_gcd >> expected_decimal))
            throw std::runtime_error("Missing or malformed B76 GCD/decimal TSV row");
        if (source_iteration!=repetition || source_n!=n || source_m!=m)
            throw std::runtime_error("Wrong source index/width in B76 TSV");

        check_against_precomputed_hex(x,source_x);
        check_against_precomputed_hex(y,source_y);
        check_against_precomputed_hex(gcd(x,y),expected_gcd);
        check_against_precomputed_hex(gcd(y,x),expected_gcd);
        const auto decimal=to_decimal(x);
        ARIADNE_TEST_ASSERT(compare(from_decimal(decimal),x)==0);
        ARIADNE_TEST_ASSERT(decimal.front()>='0' && decimal.front()<='9');
        if (!x.is_zero()) ARIADNE_TEST_ASSERT(decimal.front()!='0');
        ARIADNE_TEST_ASSERT(decimal==expected_decimal);
        check_against_precomputed_hex(x,source_x);
        check_against_precomputed_hex(y,source_y);
    }
    std::string b76_extra;
    if (b76_fixture >> b76_extra)
        throw std::runtime_error("Extra B76 GCD/decimal TSV rows");
    // B131: preserve the original 1,800 randomized decimal strings and
    // exact mt19937_64 stream, but consume source-keyed independent TSV
    // expectations instead of computing the numerical oracle with GMP.
    std::ifstream b131_decimal_random(
        std::string(ARIADNE_TEST_DATA_DIR)+
        "/big_unsigned_b131_decimal_random.tsv");
    if (!b131_decimal_random)
        throw std::runtime_error("Missing B131 randomized decimal TSV");
    std::string b131_header;
    while (b131_decimal_random.peek()=='#')
        std::getline(b131_decimal_random,b131_header);
    for (unsigned repetition=0u;repetition<1800u;++repetition) {
        const std::size_t length=repetition%149u==0u?1900u:
                                 1u+engine()%370u;
        std::string digits(length,'0');
        for (auto& digit:digits) {
            digit=static_cast<char>('0'+static_cast<int>(engine()%10u));
        }
        unsigned source_repetition=0u;
        std::size_t source_length=0u;
        std::string source_digits,expected_hex,expected_decimal;
        if (!(b131_decimal_random >> source_repetition >> source_length
              >> source_digits >> expected_hex >> expected_decimal) ||
            source_repetition!=repetition || source_length!=length ||
            source_digits!=digits)
            throw std::runtime_error(
                "Missing or unordered B131 randomized decimal TSV row");
        const auto parsed=from_decimal(digits);
        check_against_precomputed_hex(parsed,expected_hex);
        ARIADNE_TEST_ASSERT(to_decimal(parsed)==
            to_decimal(from_decimal(to_decimal(parsed))));
        ARIADNE_TEST_ASSERT(to_decimal(parsed)==expected_decimal);
        ARIADNE_TEST_EQUALS(to_decimal_cutoff16_diagnostic(parsed),expected_decimal);
        ARIADNE_TEST_EQUALS(to_decimal_cutoff32_diagnostic(parsed),expected_decimal);
    }
    std::string b131_extra;
    if (b131_decimal_random >> b131_extra)
        throw std::runtime_error("Extra B131 randomized decimal TSV rows");
    // B132: all remaining structured decimal parsing cases use one
    // sequential source-keyed fixture instead of a runtime GMP oracle.
    std::ifstream b132_decimal_structured(
        std::string(ARIADNE_TEST_DATA_DIR)+
        "/big_unsigned_b132_decimal_structured.tsv");
    if (!b132_decimal_structured)
        throw std::runtime_error("Missing B132 structured decimal TSV");
    std::string b132_header;
    while (b132_decimal_structured.peek()=='#')
        std::getline(b132_decimal_structured,b132_header);
    std::size_t b132_ordinal=0u;
    const auto read_a21=[&]() {
        std::size_t source_ordinal=0u;
        std::string source_decimal,expected_hex,expected_decimal;
        if (!(b132_decimal_structured >> source_ordinal >> source_decimal
              >> expected_hex >> expected_decimal) ||
            source_ordinal!=b132_ordinal++)
            throw std::runtime_error(
                "Missing or unordered B132 structured decimal TSV row");
        return std::array<std::string,3u>{
            std::move(source_decimal),std::move(expected_hex),
            std::move(expected_decimal)};
    };
    const auto check_long_parse=[&](std::string const& input) {
        const auto reference=read_a21();
        if (reference[0u]!=input)
            throw std::runtime_error(
                "Wrong source in B132 structured decimal TSV row");
        const BigUnsigned parsed=from_decimal(input);
        check_against_precomputed_hex(parsed,reference[1u]);
        ARIADNE_TEST_ASSERT(to_decimal(parsed)==reference[2u]);
        return parsed;
    };
    const auto check_long_parse_fixture_source=[&]() {
        const auto reference=read_a21();
        const BigUnsigned parsed=from_decimal(reference[0u]);
        check_against_precomputed_hex(parsed,reference[1u]);
        ARIADNE_TEST_ASSERT(to_decimal(parsed)==reference[2u]);
        return parsed;
    };
    // Direct result construction may reserve more limbs than the
    // canonical value needs. Check heap-to-inline normalization and
    // reject malformed bytes at every position around both the
    // 32-digit radix and 77/78-digit storage boundaries.
    for (const std::size_t length:{1u,19u,31u,32u,76u,77u,78u,95u}) {
        check_long_parse(std::string(length,'0'));
        check_long_parse(std::string(length-1u,'0')+"7");
        check_long_parse(std::string(length,'9'));
        for (std::size_t position=0u;position<length;++position) {
            for (const char invalid:{'x','\0'}) {
                std::string input(length,'9');
                input[position]=invalid;
                bool threw=false;
                try { (void)from_decimal(input); }
                catch (std::invalid_argument const&) { threw=true; }
                ARIADNE_TEST_ASSERT(threw);
            }
        }
    }
    for (const std::size_t length:
            {308u,309u,310u,311u,312u,616u,617u,1234u,1900u}) {
        check_long_parse(std::string(length,'9'));
        check_long_parse("1"+std::string(length-2u,'0')+"1");
        // The long path is selected by raw input size, even if
        // leading zeros mean that the stored magnitude is small.
        check_long_parse(std::string(length-1u,'0')+"1");
        std::string mixed(length,'0');
        mixed.front()='8';
        mixed[18u]='7';
        mixed.back()='6';
        check_long_parse(mixed);
        std::string alternating(length,'9');
        for (std::size_t i=1u;i<length;i+=2u) {
            alternating[i]='0';
        }
        check_long_parse(alternating);
    }
    // B47: the wide input path now starts at 32 digits, including
    // 128-, 256- and 1024-bit benchmark magnitudes. Stress first
    // short groups, full 18-digit boundaries and both sides of
    // the old 309/310 threshold against the independent B132 TSV.
    for (const std::size_t length:
            {16u,17u,18u,19u,31u,32u,33u,35u,36u,37u,38u,
             39u,40u,53u,54u,55u,76u,77u,78u,79u,
             307u,308u,309u,310u,311u,617u}) {
        check_long_parse(std::string(length,'9'));
        check_long_parse(std::string(length,'0'));
        check_long_parse("1"+std::string(length-2u,'0')+"1");
        check_long_parse(std::string(length-1u,'0')+"7");
        std::string alternating(length,'0');
        for (std::size_t i=0u;i<length;i+=2u) {
            alternating[i]='9';
        }
        check_long_parse(alternating);
        std::string random_digits(length,'0');
        for (auto& digit:random_digits) {
            digit=static_cast<char>('0'+engine()%10u);
        }
        check_long_parse(random_digits);
    }
    // B50: radix-10^19 parsing must handle all boundaries of
    // complete 19-digit groups, the 31/32-digit selector, the
    // 308-digit 1024-bit input, and 2048/4096-bit inputs. Check
    // exact parsed magnitudes against B132, including leading zeros,
    // all nines, sparse digits, and independently randomized text.
    for (const std::size_t length:
            {18u,19u,20u,31u,32u,33u,36u,37u,38u,39u,40u,
             56u,57u,58u,75u,76u,77u,78u,94u,95u,96u,
             113u,114u,115u,151u,152u,153u,307u,308u,
             309u,310u,311u,616u,617u,618u,619u,
             1233u,1234u,1235u}) {
        check_long_parse(std::string(length,'9'));
        check_long_parse(std::string(length,'0'));
        check_long_parse("1"+std::string(length-2u,'0')+"1");
        check_long_parse(std::string(length-1u,'0')+"7");
        std::string alternating(length,'0');
        for (std::size_t i=0u;i<length;i+=2u) alternating[i]='9';
        check_long_parse(alternating);
        std::string random_digits(length,'0');
        for (auto& digit:random_digits) {
            digit=static_cast<char>('0'+engine()%10u);
        }
        check_long_parse(random_digits);
    }

    // B47: preserve exact powers 2^(bits-1) plus offsets. The source
    // decimal is frozen in B132; native shift/add additionally ties each
    // sequential fixture row to the original bits/offset case.
    for (const unsigned bits:{128u,256u,1024u,2048u,4096u}) {
        for (const unsigned long offset:{0ul,13ul,31ul,75ul,187ul}) {
            const BigUnsigned parsed=check_long_parse_fixture_source();
            const BigUnsigned expected=add(
                shift_left(one,bits-1u),BigUnsigned(offset));
            ARIADNE_TEST_ASSERT(compare(parsed,expected)==0);
        }
    }
    if (b132_ordinal!=484u)
        throw std::runtime_error("Wrong B132 structured decimal row count");
    std::string b132_extra;
    if (b132_decimal_structured >> b132_extra)
        throw std::runtime_error("Extra B132 structured decimal TSV rows");

    // The long-string fast path must never bypass the shared
    // validation of every ASCII digit (including the final byte).
    for (const std::size_t invalid_position:
            {0u,155u,309u,310u}) {
        std::string invalid(311u,'7');
        invalid[invalid_position]='x';
        bool threw=false;
        try { (void)from_decimal(invalid); }
        catch (std::invalid_argument const&) { threw=true; }
        ARIADNE_TEST_ASSERT(threw);
    }

    // B133: all decimal-formatting cases use a source-keyed independent
    // TSV oracle. Native sources keep the original construction paths.
    std::ifstream b133_decimal_format(
        std::string(ARIADNE_TEST_DATA_DIR)+
        "/big_unsigned_b133_decimal_format.tsv");
    if (!b133_decimal_format)
        throw std::runtime_error("Missing B133 decimal formatting TSV");
    std::string b133_header;
    while (b133_decimal_format.peek()=='#')
        std::getline(b133_decimal_format,b133_header);
    std::size_t b133_ordinal=0u;
    const auto check_decimal_boundary=[&](BigUnsigned const& native) {
        std::size_t source_ordinal=0u;
        std::string source_hex,expected_decimal;
        if (!(b133_decimal_format >> source_ordinal >> source_hex
              >> expected_decimal) ||
            source_ordinal!=b133_ordinal++)
            throw std::runtime_error(
                "Missing or unordered B133 decimal formatting TSV row");
        check_against_precomputed_hex(native,source_hex);
        const std::string decimal=to_decimal(native);
        ARIADNE_TEST_ASSERT(decimal==expected_decimal);
        ARIADNE_TEST_ASSERT(
            compare(from_decimal(decimal),native)==0);
    };
    // B25: independent B133 oracle regression around the 8/9-limb
    // radix-10^9 / radix-10^18 selector and at recursive leaf sizes.
    // Include all-ones carries, powers of the binary limb radix,
    // long internal zeros and random 64-bit limb patterns.
    std::mt19937_64 leaf_engine(0xB250C0DEULL);
    for (const std::size_t count:
            {7u,8u,9u,10u,15u,16u,17u,32u,64u,128u}) {
        for (unsigned mode=0u;mode<5u;++mode) {
            std::vector<std::uint64_t> words(count,0u);
            if (mode==0u) {
                words.back()=1u;
            } else if (mode==1u) {
                std::fill(words.begin(),words.end(),
                          std::numeric_limits<std::uint64_t>::max());
            } else if (mode==2u) {
                words.front()=1u;
                words.back()=1u;
            } else {
                for (auto& word:words) word=leaf_engine();
                if (mode==4u) {
                    for (std::size_t i=1u;i+1u<count;++i) words[i]=0u;
                }
                words.back()|=std::uint64_t{1}<<63u;
            }
            check_decimal_boundary(BigUnsigned(words));
        }
    }

    // B25.1: on supported ARM64 these exercise the 18-digit leaf
    // formatter, while x86-64/portable builds must match B133 via
    // the B19 9-digit fallback. Probe carries and zero-padding at
    // decimal block boundaries near the 8/9/16-limb cutoffs.
    for (const unsigned long groups: {8ul,9ul,10ul,16ul,17ul,18ul}) {
        const unsigned long digits=18ul*groups;
        const BigUnsigned power=from_decimal(
            "1"+std::string(digits,'0'));
        check_decimal_boundary(power);
        check_decimal_boundary(sub(power,one));
        check_decimal_boundary(add(power,one));
    }

    // B30/B30.1: preserve the 8/9-limb boundary and 10^18
    // decimal carry fixtures against the independent B133 oracle.
    // ARM64 selects radix 10^18; x86-64 and other targets now
    // exercise the restored B29 radix-10^9 implementation.
    std::mt19937_64 divq_decimal_engine(0xB300DEC1ULL);
    for (const std::size_t count:
            {8u,9u,10u,12u,16u,17u,32u,64u}) {
        for (unsigned mode=0u;mode<16u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=divq_decimal_engine();
            if (mode%8u==0u) {
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()=1u;
            } else if (mode%8u==1u) {
                std::fill(words.begin(),words.end(),
                          std::numeric_limits<std::uint64_t>::max());
            } else if (mode%8u==2u) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=1u;
            } else if (mode%8u==3u) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::numeric_limits<std::uint64_t>::max();
            } else if (mode%8u==4u) {
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()=std::uint64_t{1}<<63u;
            } else if (mode%8u==5u) {
                for (std::size_t i=1u;i+1u<count;++i) {
                    words[i]=0u;
                }
            }
            words.back()|=1u;
            check_decimal_boundary(BigUnsigned(words));
        }
    }

    // B42: B133 oracle regression at the 8/9/16/17-limb wide
    // decimal formatter boundaries, plus large dense/sparse inputs.
    std::mt19937_64 b42_engine(0xB42DEC1A7ULL);
    for (const std::size_t count:
            {8u,9u,10u,15u,16u,17u,18u,32u,33u,64u,65u,128u,256u}) {
        for (unsigned mode=0u;mode<16u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& w:words) w=b42_engine();
            switch (mode%8u) {
            case 0u: std::fill(words.begin(),words.end(),0u); break;
            case 1u:
                std::fill(words.begin(),words.end(),
                          std::numeric_limits<std::uint64_t>::max());
                break;
            case 2u:
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
                break;
            case 3u:
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()=std::uint64_t{1}<<63u;
                break;
            case 4u:
                for (std::size_t i=0u;i<count;i+=2u) words[i]=0u;
                break;
            case 5u:
                std::fill(words.begin(),words.end(),
                          std::uint64_t{999999999999999999ULL});
                break;
            case 6u: words.back()=1u; break;
            default: break;
            }
            check_decimal_boundary(BigUnsigned(words));
        }
    }

    // B46: reciprocal-based ARM64 radix-10^18 leaf divrem replaces
    // the unsigned 128/64 runtime divide. Exercise random dense and
    // sparse 9..16-limb leaves, their 8/9 dispatch edge, the 16/17
    // recursion boundary, carries and long internal zero runs. On
    // every host compare against B133, not a native roundtrip alone.
    std::mt19937_64 b46_engine(0xB460DEC1A7ULL);
    for (const std::size_t count:
            {8u,9u,10u,11u,12u,13u,14u,15u,16u,17u,32u}) {
        for (unsigned mode=0u;mode<32u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=b46_engine();
            if (mode%8u==0u) {
                std::fill(words.begin(),words.end(),
                          std::numeric_limits<std::uint64_t>::max());
            } else if (mode%8u==1u) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=1u;
            } else if (mode%8u==2u) {
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()=1u;
            } else if (mode%8u==3u) {
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
            } else if (mode%8u==4u) {
                std::fill_n(words.begin(),count/2u,0u);
            } else if (mode%8u==5u) {
                for (std::size_t j=0u;j<count;j+=2u) {
                    words[j]=0u;
                }
            } else if (mode%8u==6u) {
                std::fill(words.begin(),words.end(),
                          std::uint64_t{999999999999999999ULL});
            }
            words.back()|=1u;
            check_decimal_boundary(BigUnsigned(words));
        }
    }

    // B51: the multiply-only B46 wide leaf is now also selected
    // for x86-64 GCC/Clang at 9+ limbs. Independent B133 oracles
    // cover the 8/9-limb selector, the 16/17-limb recursive cutoff,
    // dense/sparse operands, carry boundaries and zero runs. Keep
    // the portable/MSVC radix-10^9 regression path equally checked.
    std::mt19937_64 b51_engine(0xB510D3C1ULL);
    for (const std::size_t count:
            {7u,8u,9u,10u,12u,15u,16u,17u,18u,32u,33u,64u}) {
        for (unsigned mode=0u;mode<24u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=b51_engine();
            switch (mode%8u) {
            case 0u:
                std::fill(words.begin(),words.end(),
                          std::numeric_limits<std::uint64_t>::max());
                break;
            case 1u:
                std::fill(words.begin(),words.end(),0u);
                words.back()=1u;
                break;
            case 2u:
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()=1u;
                break;
            case 3u:
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
                break;
            case 4u:
                std::fill_n(words.begin(),count/2u,0u);
                break;
            case 5u:
                for (std::size_t i=0u;i<count;i+=2u) {
                    words[i]=0u;
                }
                break;
            case 6u:
                std::fill(words.begin(),words.end(),
                          std::uint64_t{999999999999999999ULL});
                break;
            default: break;
            }
            words.back()|=1u;
            check_decimal_boundary(BigUnsigned(words));
        }
    }

    // B59: oracle-checked x86 GCC 16 -> 32-limb recursion cutoff,
    // especially 31/32/33-limb transition and 63/64/65 inputs.
    // Other hosts use the unmodified 16-limb selector. Reuse the
    // independent B133 source-keyed oracle.
    std::mt19937_64 b59_engine(0xB590DEC1ULL);
    for (const std::size_t count:
            {15u,16u,17u,31u,32u,33u,63u,64u,65u}) {
        for (unsigned mode=0u;mode<24u;++mode) {
            std::vector<std::uint64_t> words(count);
            for (auto& word:words) word=b59_engine();
            switch (mode%8u) {
            case 0u:
                std::fill(words.begin(),words.end(),0u);
                words.back()=1u;
                break;
            case 1u:
                std::fill(words.begin(),words.end(),
                          std::numeric_limits<std::uint64_t>::max());
                break;
            case 2u:
                std::fill(words.begin(),words.end(),0u);
                words.front()=1u;
                words.back()=std::numeric_limits<std::uint64_t>::max();
                break;
            case 3u:
                std::fill(words.begin(),words.end(),0u);
                words.back()=std::uint64_t{1}<<63u;
                break;
            case 4u:
                for (std::size_t i=0u;i<count;i+=2u) words[i]=0u;
                break;
            case 5u:
                std::fill(words.begin(),words.end(),
                          std::uint64_t{999999999999999999ULL});
                break;
            case 6u:
                for (std::size_t i=1u;i+1u<count;++i) words[i]=0u;
                break;
            default: break;
            }
            // Keep the raw magnitude at the intended selector size.
            words.back()|=std::uint64_t{1}<<63u;
            check_decimal_boundary(BigUnsigned(words));
        }
    }

    // B43: append-only recursive formatting must preserve the
    // *entire* fixed-width lower field even when the remainder
    // itself needs further recursive divisions. B133 independently
    // supplies the reference, never a native decimal roundtrip.
    for (const unsigned long digits:
            {144ul,288ul,576ul,1152ul,2304ul}) {
        const BigUnsigned power=from_decimal(
            "1"+std::string(digits,'0'));
        for (const unsigned long factor:
                {1ul,7ul,1234567ul,9999999ul}) {
            const BigUnsigned base=mul(power,BigUnsigned(factor));
            for (const unsigned long offset:
                    {0ul,1ul,11ul,99ul}) {
                check_decimal_boundary(
                    add(base,BigUnsigned(offset)));
            }
        }
    }

    // B31: the new pair-wise base-10^9 renderer must preserve every
    // zero-padded two-digit group (00, 01, 09, 10, 99), carry across
    // decimal blocks, and the variable-width most-significant block.
    // B133 supplies the independent exact decimal oracle on every host.
    for (const unsigned long digits:
            {9ul,18ul,27ul,72ul,144ul,288ul,576ul}) {
        const BigUnsigned power=from_decimal(
            "1"+std::string(digits,'0'));
        for (const unsigned long offset:
                {0ul,1ul,9ul,10ul,11ul,99ul,100ul,101ul,999ul}) {
            check_decimal_boundary(add(power,BigUnsigned(offset)));
            if (offset!=0ul)
                check_decimal_boundary(sub(power,BigUnsigned(offset)));
        }
    }

    for (unsigned k=0u;k<=9u;++k) {
        const unsigned long digits=9ul*(1ul<<k);
        const BigUnsigned power=from_decimal(
            "1"+std::string(digits,'0'));
        check_decimal_boundary(power);
        check_decimal_boundary(sub(power,one));
        check_decimal_boundary(add(power,one));
        check_decimal_boundary(add(power,BigUnsigned(std::uint64_t{9})));

        // A nonzero suffix after many complete zero-filled radix
        // groups; tests width padding in recursive low partitions.
        const std::string blocks=
            std::string("123456789")+std::string(digits,'0')+
            std::string("000000001");
        check_decimal_boundary(from_decimal(blocks));
    }
    if (b133_ordinal!=1509u)
        throw std::runtime_error("Wrong B133 decimal formatting row count");
    std::string b133_extra;
    if (b133_decimal_format >> b133_extra)
        throw std::runtime_error("Extra B133 decimal formatting TSV rows");
}

void TestBigUnsigned::test_randomized()
{
    std::mt19937_64 engine(0xBA17B00DULL);

    std::ifstream b123_golden(ARIADNE_BIG_UNSIGNED_B123_FIRST_PATH);
    if (!b123_golden) throw std::runtime_error("Missing B123 randomized TSV");
    std::string b123_header;
    while (b123_golden.peek()=='#') std::getline(b123_golden,b123_header);
    std::ifstream b124_golden(ARIADNE_BIG_UNSIGNED_B124_SECOND_PATH);
    if (!b124_golden) throw std::runtime_error("Missing B124 randomized TSV");
    std::string b124_header;
    while (b124_golden.peek()=='#') std::getline(b124_golden,b124_header);
    std::ifstream b125_golden(ARIADNE_BIG_UNSIGNED_B125_THIRD_PATH);
    if (!b125_golden) throw std::runtime_error("Missing B125 randomized TSV");
    std::string b125_header;
    while (b125_golden.peek()=='#') std::getline(b125_golden,b125_header);
    std::ifstream b126_golden(ARIADNE_BIG_UNSIGNED_B126_FOURTH_PATH);
    if (!b126_golden) throw std::runtime_error("Missing B126 randomized TSV");
    std::string b126_header;
    while (b126_golden.peek()=='#') std::getline(b126_golden,b126_header);
    std::ifstream b127_golden(ARIADNE_BIG_UNSIGNED_B127_FIFTH_PATH);
    if (!b127_golden) throw std::runtime_error("Missing B127 randomized TSV");
    std::string b127_header;
    while (b127_golden.peek()=='#') std::getline(b127_golden,b127_header);
    std::ifstream b128_golden(ARIADNE_BIG_UNSIGNED_B128_SIXTH_PATH);
    if (!b128_golden) throw std::runtime_error("Missing B128 randomized TSV");
    std::string b128_header;
    while (b128_golden.peek()=='#') std::getline(b128_golden,b128_header);
    std::ifstream b129_golden(ARIADNE_BIG_UNSIGNED_B129_SEVENTH_PATH);
    if (!b129_golden) throw std::runtime_error("Missing B129 randomized TSV");
    std::string b129_header;
    while (b129_golden.peek()=='#') std::getline(b129_golden,b129_header);
    std::ifstream b130_golden(ARIADNE_BIG_UNSIGNED_B130_EIGHTH_PATH);
    if (!b130_golden) throw std::runtime_error("Missing B130 randomized TSV");
    std::string b130_header;
    while (b130_golden.peek()=='#') std::getline(b130_golden,b130_header);
    for (unsigned repetition=0; repetition<10000u; ++repetition) {
        const std::size_t n = engine()%20u;
        const std::size_t m = engine()%20u;
        std::array<std::uint64_t,20> na{}, nb{};
        for (std::size_t j=0; j<n; ++j) na[j] = engine();
        for (std::size_t j=0; j<m; ++j) nb[j] = engine();
        if (repetition%7u == 0u) na.fill(std::numeric_limits<std::uint64_t>::max());
        if (repetition%11u == 0u) nb.fill(0u);
        const BigUnsigned x(std::span(na.data(),n));
        const BigUnsigned y(std::span(nb.data(),m));
        const std::size_t bits=engine()%1030u;
        std::ifstream& golden=repetition<1250u?b123_golden:
            repetition<2500u?b124_golden:
            repetition<3750u?b125_golden:
            repetition<5000u?b126_golden:
            repetition<6250u?b127_golden:
            repetition<7500u?b128_golden:
            repetition<8750u?b129_golden:b130_golden;
        unsigned ref_rep=0u;
        std::size_t ref_n=0u,ref_m=0u,ref_bits=0u;
        int ref_cmp=0;
        std::string xhex,yhex,sumhex,diffhex,lhex,rhex;
        if (!(golden >> ref_rep >> ref_n >> ref_m >> ref_bits
            >> xhex >> yhex >> ref_cmp >> sumhex >> diffhex
            >> lhex >> rhex) ||
            ref_rep!=repetition || ref_n!=n || ref_m!=m ||
            ref_bits!=bits || ref_cmp < -1 || ref_cmp > 1)
            throw std::runtime_error("Missing or unordered B123-B130 randomized row");
        check_against_precomputed_hex(x,xhex);
        check_against_precomputed_hex(y,yhex);
        ARIADNE_TEST_ASSERT(compare(x,y)==ref_cmp);
        check_against_precomputed_hex(add(x,y),sumhex);
        if (ref_cmp>=0) {
            if (diffhex=="-")
                throw std::runtime_error("Missing B123-B130 subtraction oracle");
            check_against_precomputed_hex(sub(x,y),diffhex);
        } else if (diffhex!="-")
            throw std::runtime_error("Unexpected B123-B130 subtraction oracle");
        check_against_precomputed_hex(shift_left(x,bits),lhex);
        check_against_precomputed_hex(shift_right(x,bits),rhex);
        // The borrowed kernel uses the same independently frozen result,
        // including raw leading zeros and a poisoned oversized destination.
        std::array<BigUnsigned::Limb,22u> shifted_words;
        shifted_words.fill(std::numeric_limits<BigUnsigned::Limb>::max());
        const auto shifted_used=BigUnsigned::shift_right_into(
            BigUnsigned::LimbView(na.data(),n),bits,shifted_words);
        check_against_precomputed_hex(BigUnsigned(BigUnsigned::LimbView(
            shifted_words.data(),shifted_used)),rhex);
        ARIADNE_TEST_ASSERT(std::all_of(shifted_words.begin()+shifted_used,
            shifted_words.end(),[](auto digit){return digit==0u;}));
        check_against_precomputed_hex(x,xhex);
        check_against_precomputed_hex(y,yhex);
    }
    std::string b123_extra;
    if (b123_golden >> b123_extra)
        throw std::runtime_error("Extra B123 randomized TSV rows");
    std::string b124_extra;
    if (b124_golden >> b124_extra)
        throw std::runtime_error("Extra B124 randomized TSV rows");
    std::string b125_extra;
    if (b125_golden >> b125_extra)
        throw std::runtime_error("Extra B125 randomized TSV rows");
    std::string b126_extra;
    if (b126_golden >> b126_extra)
        throw std::runtime_error("Extra B126 randomized TSV rows");
    std::string b127_extra;
    if (b127_golden >> b127_extra)
        throw std::runtime_error("Extra B127 randomized TSV rows");
    std::string b128_extra;
    if (b128_golden >> b128_extra)
        throw std::runtime_error("Extra B128 randomized TSV rows");
    std::string b129_extra;
    if (b129_golden >> b129_extra)
        throw std::runtime_error("Extra B129 randomized TSV rows");
    std::string b130_extra;
    if (b130_golden >> b130_extra)
        throw std::runtime_error("Extra B130 randomized TSV rows");
}

} // namespace

void test_big_unsigned_contracts_and_controls() {
    using Limb=BigUnsigned::Limb;
    const BigUnsigned a(17u), b(5u);
    const std::array<Limb,2u> wide={3u,1u};
    std::array<Limb,2u> output={91u,92u};
    std::array<Limb,2u> remainder={81u,82u};
    const auto empty=BigUnsigned::MutableLimbView{};
    ARIADNE_TEST_FAIL(BigUnsigned::add_into(a.limbs(),b.limbs(),empty));
    ARIADNE_TEST_FAIL(BigUnsigned::sub_into(a.limbs(),b.limbs(),empty));
    ARIADNE_TEST_EQUALS(output[0u],91u);
    ARIADNE_TEST_FAIL(BigUnsigned::add_shifted(a.limbs(),64u,b.limbs(),0u));
    ARIADNE_TEST_FAIL(BigUnsigned::add_shifted(a.limbs(),0u,b.limbs(),64u));
    ARIADNE_TEST_FAIL(BigUnsigned::sub_shifted(a.limbs(),64u,b.limbs(),0u));
    ARIADNE_TEST_FAIL(BigUnsigned::sub_shifted(a.limbs(),0u,b.limbs(),64u));
    ARIADNE_TEST_FAIL(BigUnsigned::divrem_into(a.limbs(),{},output,remainder,{}));
    ARIADNE_TEST_FAIL(BigUnsigned::divrem_into(a.limbs(),b.limbs(),empty,remainder,{}));
    ARIADNE_TEST_FAIL(BigUnsigned::divrem_into(a.limbs(),b.limbs(),output,empty,{}));
    const std::array<Limb,2u> divisor_words={2u,1u};
    ARIADNE_TEST_FAIL(BigUnsigned::divrem_into(wide,divisor_words,output,remainder,{}));
    // Every preflight failure leaves both destination spans untouched.
    ARIADNE_TEST_EQUALS(output[0u],91u);
    ARIADNE_TEST_EQUALS(output[1u],92u);
    ARIADNE_TEST_EQUALS(remainder[0u],81u);
    ARIADNE_TEST_EQUALS(remainder[1u],82u);
    // Scratch contents are unspecified by the public API. Poison every
    // word: (B^2+5)/(B+3) = B-3, with remainder 14, B=2^64.
    const std::array<Limb,3u> dividend_words={5u,0u,1u};
    const std::array<Limb,2u> dirty_divisor={3u,1u};
    std::array<Limb,3u> dirty_quotient;
    std::array<Limb,2u> dirty_remainder;
    std::array<Limb,6u> dirty_scratch;
    dirty_scratch.fill(std::numeric_limits<Limb>::max());
    BigUnsigned::divrem_into(dividend_words,dirty_divisor,dirty_quotient,dirty_remainder,dirty_scratch);
    ARIADNE_TEST_EQUALS(dirty_quotient[0u],std::numeric_limits<Limb>::max()-2u);
    ARIADNE_TEST_EQUALS(dirty_quotient[1u],0u);
    ARIADNE_TEST_EQUALS(dirty_quotient[2u],0u);
    ARIADNE_TEST_EQUALS(dirty_remainder[0u],14u);
    ARIADNE_TEST_EQUALS(dirty_remainder[1u],0u);
    ARIADNE_TEST_FAIL(BigUnsigned::mul_into_karatsuba(a.limbs(),b.limbs(),output,{},7u));
    ARIADNE_TEST_FAIL(BigUnsigned::divrem_scratch_limbs(
        std::numeric_limits<std::size_t>::max(),2u));
    ARIADNE_TEST_FAIL(BigUnsigned::divrem_scratch_limbs(
        std::numeric_limits<std::size_t>::max(),std::numeric_limits<std::size_t>::max()));
    ARIADNE_TEST_FAIL(BigUnsigned::karatsuba_scratch_limbs(10u,10u,7u));
    ARIADNE_TEST_FAIL(BigUnsigned::karatsuba_scratch_limbs(
        std::numeric_limits<std::size_t>::max(),1u));
    ARIADNE_TEST_EQUALS(BigUnsigned::karatsuba_scratch_limbs(
        std::numeric_limits<std::size_t>::max(),0u),0u);
    for (auto const& input : {BigUnsigned{},a,BigUnsigned(wide)}) {
        // 17 mod 5 = 2; (2^64+3) mod 5 = 4.
        const Limb expected=input.is_zero() ? 0u : input.size()==1u ? 2u : 4u;
        ARIADNE_TEST_EQUALS(rem_limb_divq_diagnostic(input,b).limb(0u),expected);
        ARIADNE_TEST_EQUALS(rem_limb_preinverse_diagnostic(input,b).limb(0u),expected);
    }
    ARIADNE_TEST_ASSERT(compare(rem_limb_divq_diagnostic(a,BigUnsigned(wide)),a)==0);
    ARIADNE_TEST_EQUALS(mul_short_stack_selected_diagnostic(a,b).limb(0u),85u);
    ARIADNE_TEST_EQUALS(mul_short_stack_fused_diagnostic(a,b).limb(0u),85u);
    ARIADNE_TEST_ASSERT(mul_short_stack_selected_diagnostic(a,BigUnsigned{}).is_zero());
    ARIADNE_TEST_ASSERT(mul_short_stack_fused_diagnostic(BigUnsigned{},b).is_zero());
    ARIADNE_TEST_EQUALS(square_short_symmetric_diagnostic(a).limb(0u),289u);
    ARIADNE_TEST_ASSERT(square_short_symmetric_diagnostic(BigUnsigned{}).is_zero());
    ARIADNE_TEST_EQUALS(to_decimal_cutoff16_diagnostic(BigUnsigned{}),std::string("0"));
    ARIADNE_TEST_EQUALS(to_decimal_cutoff32_diagnostic(BigUnsigned{}),std::string("0"));
#if defined(__aarch64__) || defined(_M_ARM64)
    // CMake builds the ARM64 assembly kernel only for GNU/Clang.
    // MSVC ARM64 deliberately uses the tested portable C++20 fallback.
#if defined(__GNUC__) || defined(__clang__)
    ARIADNE_TEST_ASSERT(BigUnsigned::has_native_multiply());
    ARIADNE_TEST_ASSERT(BigUnsigned::has_arm64_multiply());
#else
    ARIADNE_TEST_ASSERT(!BigUnsigned::has_native_multiply());
    ARIADNE_TEST_ASSERT(!BigUnsigned::has_arm64_multiply());
#endif
    ARIADNE_TEST_ASSERT(!BigUnsigned::has_x86_64_carry_multiply());
    ARIADNE_TEST_ASSERT(!BigUnsigned::has_x86_64_multiply());
#elif defined(__x86_64__) || defined(_M_X64)
    ARIADNE_TEST_ASSERT(BigUnsigned::has_native_multiply());
    // The MULQ/ADCQ carry-chain kernel uses GNU-style inline assembly.
    // MSVC (including clang-cl) retains the native x86-64 kernel but
    // does not build this GNU/Clang-only carry variant.
#if (defined(__GNUC__) || defined(__clang__)) && !defined(_MSC_VER)
    ARIADNE_TEST_ASSERT(BigUnsigned::has_x86_64_carry_multiply());
#else
    ARIADNE_TEST_ASSERT(!BigUnsigned::has_x86_64_carry_multiply());
#endif
    ARIADNE_TEST_ASSERT(BigUnsigned::has_x86_64_multiply());
#endif
}

int main()
{
    ARIADNE_TEST_CALL(test_big_unsigned_contracts_and_controls());
    ARIADNE_TEST_CLASS(BigUnsigned,TestBigUnsigned());
    return ARIADNE_TEST_FAILURES;
}
