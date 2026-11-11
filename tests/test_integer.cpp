/***************************************************************************
 *            test_integer.cpp
 *
 *  Copyright  2013-20  Pieter Collins
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "foundation/logical.hpp"
#include "numeric/accuracy.hpp"
#include "numeric/int.hpp"
#include "numeric/integer.hpp"
#include "numeric/sign.hpp"
#include "numeric/twoexp.hpp"
#include "utility/string.hpp"

inline constexpr char ARIADNE_INTEGER_SUB_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_sub.tsv";

inline constexpr char ARIADNE_INTEGER_I17_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_i17.tsv";
inline constexpr char ARIADNE_INTEGER_I18_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_i18.tsv";

inline constexpr char ARIADNE_INTEGER_I33_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_i33.tsv";
inline constexpr char ARIADNE_INTEGER_I37_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_i37.tsv";
inline constexpr char ARIADNE_INTEGER_I38_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_i38.tsv";

inline constexpr char ARIADNE_INTEGER_I66_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_i66.tsv";

#ifndef ARIADNE_TEST_DATA_DIR
#error "Provide the absolute test data directory from CMake"
#endif

inline constexpr char ARIADNE_INTEGER_POWER_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/integer_power.tsv";
inline constexpr char ARIADNE_INTEGER_SQRT_ORACLE_PATH[] = ARIADNE_TEST_DATA_DIR "/big_unsigned_sqrt.tsv";

#include "ariadne-testing.hpp"

using namespace std;
using namespace Ariadne;

namespace {

constexpr bool check_concept() {
    return requires(unsigned int m, unsigned long int lm, int n, long int ln, Integer z, Integer z2, Boolean b) {
        z = Integer();
        z = Integer(m);
        z = Integer(lm);
        z = Integer(n);
        z = Integer(ln);
        z = Integer(z);
        z2 = Integer();
        z = m;
        z = lm;
        z = n;
        z = ln;
        z = z2;
        z = +z;
        z = -z;
        z = z + z;
        z = z - z;
        z = z * z;
        z = z + n;
        z = z - n;
        z = z * n;
        z = n + z;
        z = n - z;
        z = n * z;
        z = 1_z;
        z = -1_z;
        b = (z == z);
        b = (z != z);
        b = (z <= z);
        b = (z >= z);
        b = (z < z);
        b = (z > z);
        b = (z == n);
        b = (z != n);
        b = (z <= n);
        b = (z >= n);
        b = (z < n);
        b = (z > n);
        b = (n == z);
        b = (n != z);
        b = (n <= z);
        b = (n >= z);
        b = (n < z);
        b = (n > z);
    };
}

// Parse the independently frozen hexadecimal oracle values without GMP.
BigUnsigned sqrt_oracle_hex(const std::string &hex) {
    if (hex.empty() || (hex.size() > 1u && hex.front() == '0')) {
        throw std::invalid_argument("Invalid Integer square-root oracle hex");
    }
    std::vector<std::uint64_t> words((hex.size() + 15u) / 16u, 0u);
    for (std::size_t i = 0u; i < hex.size(); ++i) {
        const char digit = hex[hex.size() - 1u - i];
        const unsigned nibble = digit >= '0' && digit <= '9'
                                  ? static_cast<unsigned>(digit - '0')
                              : digit >= 'a' && digit <= 'f'
                                  ? static_cast<unsigned>(digit - 'a' + 10)
                                  : 16u;
        if (nibble == 16u) {
            throw std::invalid_argument("Malformed Integer square-root oracle hex");
        }
        words[i / 16u] |= static_cast<std::uint64_t>(nibble)
                       << static_cast<unsigned>(4u * (i % 16u));
    }
    return BigUnsigned(words);
}

static_assert(check_concept());

} // namespace

class TestInteger : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    void test_literal();
    void test_constructors();
    void test_arithmetic();
    void test_comparisons();
    void test_small_types();
    void test_combinatorics();
    void test_native();
    void test_signed_regression();
    void test_owned_results();
    void test_sqrt();
};

void TestInteger::test_sqrt() {
    // Reuse the 161 independent restoring-digit sqrt/root/remainder
    // triples, not GMP or another call to the function under test.
    std::ifstream golden(ARIADNE_INTEGER_SQRT_ORACLE_PATH);
    if (!golden)
        throw std::runtime_error("Missing Integer sqrt oracle TSV");
    std::size_t cases = 0u;
    std::string line;
    while (std::getline(golden, line)) {
        if (line.empty() || line.front() == '#')
            continue;
        std::istringstream fields(line);
        std::string input_hex, root_hex, remainder_hex, extra;
        if (!(fields >> input_hex >> root_hex >> remainder_hex) || (fields >> extra)) {
            throw std::runtime_error("Malformed Integer sqrt oracle row");
        }
        const BigUnsigned input_magnitude = sqrt_oracle_hex(input_hex);
        const Integer input = Integer::from_magnitude(input_magnitude, false);
        const Integer preserved = input;
        const Integer expected_root = Integer::from_magnitude(
            sqrt_oracle_hex(root_hex), false);
        const Integer expected_remainder = Integer::from_magnitude(
            sqrt_oracle_hex(remainder_hex), false);
        const Integer root = isqrt(input);
        const auto [paired_root, remainder] = sqrtrem(input);
        ARIADNE_TESTING_PROVE(root == expected_root);
        ARIADNE_TESTING_PROVE(paired_root == expected_root);
        ARIADNE_TESTING_PROVE(remainder == expected_remainder);
        ARIADNE_TESTING_PROVE(!root.negative() && !remainder.negative());
        ARIADNE_TESTING_PROVE(!root.negative() && !paired_root.negative());
        ARIADNE_TESTING_PROVE(paired_root * paired_root + remainder == input);
        ARIADNE_TESTING_PROVE(remainder <= paired_root + paired_root);
        ARIADNE_TESTING_PROVE(input == preserved);

        // Only zero accepts a canonicalized negative-zero source.
        const Integer negative_input = Integer::from_magnitude(
            input_magnitude, true);
        if (input_magnitude.is_zero()) {
            ARIADNE_TESTING_PROVE(!negative_input.negative());
            ARIADNE_TESTING_PROVE(isqrt(negative_input) == Integer(0));
            const auto [zero_root, zero_remainder] = sqrtrem(negative_input);
            ARIADNE_TESTING_PROVE(zero_root == Integer(0));
            ARIADNE_TESTING_PROVE(zero_remainder == Integer(0));
        } else {
            bool root_threw = false, remainder_threw = false;
            try {
                (void)isqrt(negative_input);
            } catch (const std::domain_error &) {
                root_threw = true;
            }
            try {
                (void)sqrtrem(negative_input);
            } catch (const std::domain_error &) {
                remainder_threw = true;
            }
            ARIADNE_TESTING_PROVE(root_threw && remainder_threw);
            ARIADNE_TESTING_PROVE(negative_input.negative());
            ARIADNE_TESTING_PROVE(compare(
                                      negative_input.magnitude(), input_magnitude) == 0);
        }
        ++cases;
    }
    if (golden.bad() || cases != 166u) {
        throw std::runtime_error("Incomplete Integer sqrt oracle fixture");
    }
}

void TestInteger::test_owned_results() {
    // Exercise inline and heap magnitude ownership independently of
    // decimal parsing. Moving a result must preserve its full value,
    // leave the caller's moved-from magnitude zero and normalize -0.
    for (const std::size_t count: {0u, 1u, 4u, 5u, 8u, 9u, 64u}) {
        std::vector<std::uint64_t> words(count,
                                         std::numeric_limits<std::uint64_t>::max());
        const BigUnsigned source(words);
        for (const bool negative: {false, true}) {
            const Integer copied = Integer::from_magnitude(source, negative);
            ARIADNE_TESTING_PROVE(compare(copied.magnitude(), source) == 0);
            ARIADNE_TESTING_PROVE(copied.negative() == (negative && count != 0u));
            BigUnsigned magnitude = source;
            Integer owned = Integer::from_magnitude(std::move(magnitude), negative);
            ARIADNE_TESTING_PROVE(magnitude.is_zero());
            ARIADNE_TESTING_PROVE(owned == copied);
            const String text = owned.literal();
            ARIADNE_TESTING_PROVE(Integer(text) == copied);
            ARIADNE_TESTING_PROVE(compare(source, BigUnsigned(words)) == 0);
            std::ostringstream stream;
            stream << owned;
            ARIADNE_TESTING_PROVE(stream.str() == std::string(text));
            owned = Integer(7);
            ARIADNE_TESTING_PROVE(Integer(text) == copied);
        }
    }
    ARIADNE_TESTING_PROVE(Integer::from_magnitude(BigUnsigned{}, true).literal() == "0");
    ARIADNE_TESTING_PROVE(Integer(-1).literal() == "-1");
    ARIADNE_TESTING_PROVE(Integer(1).literal() == "1");
}

void TestInteger::test_literal() {
    Integer z(3_z);
    ARIADNE_TESTING_INSPECT("constructed z", z);
    ARIADNE_TESTING_EXPECT_EQUAL(z, 3);
    ARIADNE_TESTING_EXPECT_EQUAL(z, Integer(3));
    ARIADNE_TESTING_EXPECT_EQUAL(100000_z, Integer(100000));
    ARIADNE_TESTING_EXPECT_EQUAL(1000000000000_z, sqr(Integer(1000000)));
    ARIADNE_TESTING_EXPECT_EQUAL(1000000000000_z - sqr(Integer(1000000)), 0);
    ARIADNE_TESTING_EXPECT_EQUAL(4294967295_z, 4294967295u);
    ARIADNE_TESTING_EXPECT_EQUAL(-2147483647_z, -2147483647);
    ARIADNE_TESTING_EXPECT_EQUAL(4294967295_z, Integer(4294967295u));
    ARIADNE_TESTING_EXPECT_EQUAL(4611686016279904256_z, Integer(2147483647) * 2147483647 + 2147483647);
}

void TestInteger::test_constructors() {
    int m = 2147483647;
    unsigned int um = 2147483647u;
    unsigned long int ulm = um;
    unsigned long long int ullm = um;
    int n = -2147483647;
    long int ln = n;
    long long int lln = n;
    String sn = "-2147483647";
    String sz = "314159265358979323846264";

    ARIADNE_TESTING_PROVE((not Constructible<Integer, float>));
    ARIADNE_TESTING_PROVE((not Constructible<Integer, double>));

    Integer zum(um);
    ARIADNE_TESTING_INSPECT("constructed zum", zum);
    ARIADNE_TESTING_EXPECT_EQUAL(zum.get_si(), m);
    Integer zulm(ulm);
    ARIADNE_TESTING_INSPECT("constructed zulm", zulm);
    ARIADNE_TESTING_EXPECT_EQUAL(zulm.get_si(), m);
    Integer zullm(ullm);
    ARIADNE_TESTING_INSPECT("constructed zullm", zullm);
    ARIADNE_TESTING_EXPECT_EQUAL(zullm.get_si(), m);
    Integer zn(n);
    ARIADNE_TESTING_INSPECT("constructed zn", zn);
    ARIADNE_TESTING_EXPECT_EQUAL(zn.get_si(), n);
    Integer zln(ln);
    ARIADNE_TESTING_INSPECT("constructed zln", zln);
    ARIADNE_TESTING_EXPECT_EQUAL(zln.get_si(), n);
    Integer zlln(lln);
    ARIADNE_TESTING_INSPECT("constructed zlln", zlln);
    ARIADNE_TESTING_EXPECT_EQUAL(zlln.get_si(), n);
    Integer zsn(sn);
    ARIADNE_TESTING_INSPECT("constructed zsn", zsn);
    ARIADNE_TESTING_EXPECT_EQUAL(zsn.get_si(), n);
    Integer zsz(sz);
    ARIADNE_TESTING_INSPECT("constructed zsz", zsz);
    ARIADNE_TESTING_EXPECT_EQUAL(zsz, Integer(314159265358) * Integer(1000000000000) + Integer(979323846264));

    Integer z1(0);
    ARIADNE_TESTING_INSPECT("constructed z1", z1);
    ARIADNE_TESTING_EXPECT_EQUAL(z1.get_si(), 0);
    Integer z2(-3);
    ARIADNE_TESTING_INSPECT("constructed z2", z2);
    ARIADNE_TESTING_EXPECT_EQUAL(z2.get_si(), -3);
    Integer z3(ullm * ullm + ullm);
    ARIADNE_TESTING_INSPECT("constructed z3", z3);
    ARIADNE_TESTING_EXPECT_EQUAL(z3, zum * zum + zum);
    Integer z4(lln * lln + lln);
    ARIADNE_TESTING_INSPECT("constructed z4", z4);
    ARIADNE_TESTING_EXPECT_EQUAL(z4, zn * zn + zn);

    ARIADNE_TESTING_EXPECT_EQUAL(Integer(Nat32(7)), Integer(7));

    Integer assigned(1);
    Integer source(9);
    assigned = source;
    ARIADNE_TESTING_EXPECT_EQUAL(assigned, Integer(9));
}

void TestInteger::test_arithmetic() {
    ARIADNE_TESTING_EXPECT_EQUAL(+Integer(-5), -5);
    ARIADNE_TESTING_EXPECT_EQUAL(-Integer(-5), 5);
    ARIADNE_TESTING_EXPECT_EQUAL(Integer(-5) + Integer(2), -3);
    ARIADNE_TESTING_EXPECT_EQUAL(Integer(-5) - Integer(2), -7);
    ARIADNE_TESTING_EXPECT_EQUAL(Integer(-5) * Integer(2), -10);

    ARIADNE_TESTING_EXPECT_EQUAL(pos(Integer(-5)), -5);
    ARIADNE_TESTING_EXPECT_EQUAL(neg(Integer(-5)), 5);
    ARIADNE_TESTING_EXPECT_EQUAL(sqr(Integer(-5)), 25);
    ARIADNE_TESTING_EXPECT_EQUAL(pow(Integer(-5), 3u), -125);

    ARIADNE_TESTING_EXPECT_EQUAL((Integer)max(Integer(5), Integer(3)), 5);
    ARIADNE_TESTING_EXPECT_EQUAL((Integer)max(Integer(-5), Integer(-3)), -3);
    ARIADNE_TESTING_EXPECT_EQUAL((Integer)min(Integer(5), Integer(3)), 3);
    ARIADNE_TESTING_EXPECT_EQUAL((Integer)min(Integer(-5), Integer(-3)), -5);
    ARIADNE_TESTING_EXPECT_EQUAL(abs(Integer(-5)), 5);
    ARIADNE_TESTING_EXPECT_EQUAL(abs(Integer(0)), 0);
    ARIADNE_TESTING_EXPECT_EQUAL(abs(Integer(+5)), 5);

    Integer mutable_value(5);
    --mutable_value;
    ARIADNE_TESTING_EXPECT_EQUAL(mutable_value, Integer(4));
    mutable_value -= Integer(3);
    ARIADNE_TESTING_EXPECT_EQUAL(mutable_value, Integer(1));

    ARIADNE_TESTING_PROVE(not is_nan(Integer(1)));
    ARIADNE_TESTING_PROVE(not is_inf(Integer(1)));
    ARIADNE_TESTING_PROVE(is_finite(Integer(1)));
    ARIADNE_TESTING_PROVE(is_zero(Integer(0)));
    ARIADNE_TESTING_PROVE(not is_zero(Integer(1)));
    ARIADNE_TESTING_PROVE(not is_zero(Integer(-1)));
    ARIADNE_TESTING_EXPECT_EQUAL(sgn(Integer(-1)), Sign::NEGATIVE);
    ARIADNE_TESTING_EXPECT_EQUAL(sgn(Integer(0)), Sign::ZERO);
    ARIADNE_TESTING_EXPECT_EQUAL(sgn(Integer(1)), Sign::POSITIVE);

    ARIADNE_TESTING_EXPECT_EQUAL(cmp(Integer(2), Integer(3)), Comparison::LESS);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(Integer(3), Integer(3)), Comparison::EQUAL);
    ARIADNE_TESTING_EXPECT_EQUAL(cmp(Integer(4), Integer(3)), Comparison::GREATER);
}

void TestInteger::test_small_types() {
    ARIADNE_TESTING_EXPECT_EQUAL(Integer(Int32(-9)), Integer(-9));
    const Positive<Integer> positive_default;
    const Natural natural_default;
    ARIADNE_TESTING_EXPECT_EQUAL(positive_default, Integer(0));
    ARIADNE_TESTING_EXPECT_EQUAL(natural_default, Integer(0));

    // Exercise both concrete get() instantiations, including their public
    // failure contracts for multi-limb and excessive negative magnitudes.
    std::ifstream conversions(ARIADNE_TEST_DATA_DIR "/integer_conversion.tsv");
    if (!conversions)
        throw std::runtime_error("Missing integer_conversion.tsv");
    std::string conversion_line;
    std::set<std::string> conversion_sources;
    while (std::getline(conversions, conversion_line)) {
        if (conversion_line.empty() || conversion_line[0] == '#')
            continue;
        std::istringstream fields(conversion_line);
        std::string source, extra;
        int signed_ok, unsigned_ok, signed32_ok, unsigned32_ok;
        if (!(fields >> source >> signed_ok >> unsigned_ok >> signed32_ok >> unsigned32_ok) || (fields >> extra) || !conversion_sources.insert(source).second)
            throw std::runtime_error("Malformed integer_conversion.tsv");
        const Integer value{String(source)};
        [[maybe_unused]] auto conversion_scope = this->ctx().scope("integer conversion oracle", source);
        auto check_conversion = [&]<class T>() {
            constexpr int digits = std::numeric_limits<T>::digits;
            static_assert(digits == 31 || digits == 32 || digits == 63 || digits == 64);
            const int expected_ok = digits == 63 ? signed_ok : digits == 64 ? unsigned_ok
                                                         : digits == 31     ? signed32_ok
                                                                            : unsigned32_ok;
            bool threw = false;
            try {
                ARIADNE_TESTING_EXPECT_EQUAL(std::to_string(value.get<T>()), source);
            } catch (const std::out_of_range &) {
                threw = true;
            }
            ARIADNE_TESTING_EXPECT_EQUAL(threw, expected_ok == 0);
        };
        // These four public destination types already occur in the report.
        check_conversion.template operator()<int>();
        check_conversion.template operator()<long>();
        check_conversion.template operator()<unsigned long>();
        check_conversion.template operator()<long long>();
        check_conversion.template operator()<unsigned long long>();
    }
    if (conversion_sources.size() != 27u)
        throw std::runtime_error("Incomplete integer_conversion.tsv");

    Nat32 n32;
    Nat64 n64;
    Int32 i32;
    Int64 i64;
    ARIADNE_TESTING_EXPECT_EQUAL(n32.get_ui(), uint32_t(0));
    ARIADNE_TESTING_EXPECT_EQUAL(n64.get_ui(), uint64_t(0));
    ARIADNE_TESTING_EXPECT_EQUAL(i32.get_si(), int32_t(0));
    ARIADNE_TESTING_EXPECT_EQUAL(i64.get_si(), int64_t(0));

    ARIADNE_TESTING_EXPECT_EQUAL(Nat32(uint16_t(7)).get_ui(), uint32_t(7));
    ARIADNE_TESTING_EXPECT_EQUAL(Nat32(int16_t(7)).get_ui(), uint32_t(7));
    ARIADNE_TESTING_EXPECT_EQUAL(Nat64(uint32_t(8)).get_ui(), uint64_t(8));
    ARIADNE_TESTING_EXPECT_EQUAL(Nat64(int32_t(8)).get_ui(), uint64_t(8));
    ARIADNE_TESTING_EXPECT_EQUAL(Int32(uint16_t(9)).get_si(), int32_t(9));
    ARIADNE_TESTING_EXPECT_EQUAL(Int32(int16_t(-9)).get_si(), int32_t(-9));
    ARIADNE_TESTING_EXPECT_EQUAL(Int64(uint32_t(10)).get_si(), int64_t(10));
    ARIADNE_TESTING_EXPECT_EQUAL(Int64(uint64_t(11)).get_si(), int64_t(11));
    ARIADNE_TESTING_EXPECT_EQUAL(Int64(int16_t(-10)).get_si(), int64_t(-10));

    ARIADNE_TESTING_EXPECT_EXCEPTION(Nat32(uint64_t(4294967296ULL)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(Nat32(int64_t(-1)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(Nat32(int64_t(4294967296LL)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(Nat64(int64_t(-1)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(Int32(uint64_t(2147483648ULL)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(Int32(int64_t(-2147483649LL)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(Int32(int64_t(2147483648LL)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(Int64(uint64_t(9223372036854775808ULL)));

    Bits bits = 8_bits;
    ARIADNE_TESTING_EXPECT_EQUAL(static_cast<unsigned long int>(bits), 8ul);
    std::ostringstream bits_stream;
    bits_stream << bits;
    ARIADNE_TESTING_EXPECT_EQUAL(bits_stream.str(), std::string("8_bits"));

    Accuracy direct_accuracy(Dyadic(1, 8u));
    ARIADNE_TESTING_EXPECT_EQUAL(direct_accuracy.error(), Dyadic(1, 8u));
    ARIADNE_TESTING_EXPECT_EQUAL(accuracy(8_bits).error(), Dyadic(1, 8u));

    TwoExp p(3);
    TwoExp from_two(two);
    ARIADNE_TESTING_EXPECT_EQUAL(p.exponent(), 3);
    ARIADNE_TESTING_EXPECT_EQUAL(from_two.exponent(), 1);
    ARIADNE_TESTING_EXPECT_EQUAL(p.get_d(), 8.0);
    ARIADNE_TESTING_EXPECT_EQUAL(rec(p).exponent(), -3);
    ARIADNE_TESTING_EXPECT_EQUAL(Ariadne::exp2(-2).exponent(), -2);
    ARIADNE_TESTING_EXPECT_EQUAL((two ^ static_cast<Nat>(4)).exponent(), 4);
    ARIADNE_TESTING_EXPECT_EQUAL((two ^ static_cast<Int>(-4)).exponent(), -4);
    ARIADNE_TESTING_EXPECT_EQUAL(pow(two, -5).exponent(), -5);

    ARIADNE_TESTING_EXPECT_EQUAL(-Sign::NEGATIVE, Sign::POSITIVE);
    ARIADNE_TESTING_EXPECT_EQUAL(Sign::NEGATIVE * Sign::NEGATIVE, Sign::POSITIVE);
    ARIADNE_TESTING_EXPECT_EQUAL(Sign::NEGATIVE * Sign::POSITIVE, Sign::NEGATIVE);

    ARIADNE_TESTING_EXPECT_EQUAL(class_name<uint>(), String("uint"));
    ARIADNE_TESTING_EXPECT_EQUAL(class_name<int>(), String("int"));
    ARIADNE_TESTING_EXPECT_EQUAL(class_name<Integer>(), String("Integer"));
    ARIADNE_TESTING_EXPECT_EQUAL(class_name<Natural>(), String("Natural"));

    std::ostringstream sign_stream;
    sign_stream << Sign::ZERO << " " << Sign::NEGATIVE << " " << Sign::POSITIVE;
    ARIADNE_TESTING_EXPECT_EQUAL(sign_stream.str(), std::string("ZERO NEGATIVE POSITIVE"));

    std::ostringstream comparison_stream;
    comparison_stream << LESS << " " << EQUAL << " " << GREATER << " " << INCOMPARABLE;
    ARIADNE_TESTING_EXPECT_EQUAL(comparison_stream.str(), std::string("LESS EQUAL GREATER INCOMPARABLE"));
    std::ostringstream invalid_comparison_stream;
    ARIADNE_TESTING_EXPECT_EXCEPTION(invalid_comparison_stream << static_cast<Comparison>(42));

    Natural natural(2u);
    Natural other(3u);
    ++natural;
    natural += other;
    ARIADNE_TESTING_EXPECT_EQUAL(natural, Integer(6));
    ARIADNE_TESTING_EXPECT_EQUAL(Natural(2u) + Natural(3u), Integer(5));
    ARIADNE_TESTING_EXPECT_EQUAL(Natural(2u) * Natural(3u), Integer(6));
    ARIADNE_TESTING_EXPECT_EQUAL(log2floor(Natural(0u)), -1);
    ARIADNE_TESTING_EXPECT_EQUAL(log2floor(Natural(1u)), 0);
    ARIADNE_TESTING_EXPECT_EQUAL(log2floor(Natural(8u)), 3);
    ARIADNE_TESTING_EXPECT_EQUAL(cast_positive(Integer(4)), Integer(4));
    ARIADNE_TESTING_EXPECT_EXCEPTION(Positive<Integer>(Integer(-1)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(Natural(Integer(-1)));
}

void TestInteger::test_combinatorics() {
    ARIADNE_TESTING_EXPECT_EQUAL(fac(uint8_t(0)), uint32_t(1));
    ARIADNE_TESTING_EXPECT_EQUAL(fac(uint8_t(5)), uint32_t(120));
    ARIADNE_TESTING_EXPECT_EXCEPTION(fac(uint8_t(13)));
    ARIADNE_TESTING_EXPECT_EQUAL(bin(uint8_t(5), uint8_t(0)), uint32_t(1));
    ARIADNE_TESTING_EXPECT_EQUAL(bin(uint8_t(5), uint8_t(2)), uint32_t(10));
    ARIADNE_TESTING_EXPECT_EQUAL(bin(uint8_t(5), uint8_t(6)), uint32_t(0));
    ARIADNE_TESTING_EXPECT_EXCEPTION(bin(uint8_t(5), uint8_t(7)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(bin(uint8_t(32), uint8_t(1)));

    ARIADNE_TESTING_EXPECT_EQUAL(fac(uint16_t(0)), uint16_t(1));
    ARIADNE_TESTING_EXPECT_EQUAL(fac(uint16_t(5)), uint16_t(120));
    ARIADNE_TESTING_EXPECT_EXCEPTION(fac(uint16_t(9)));
    ARIADNE_TESTING_EXPECT_EQUAL(bin(uint16_t(5), uint16_t(0)), uint16_t(1));
    ARIADNE_TESTING_EXPECT_EQUAL(bin(uint16_t(5), uint16_t(2)), uint16_t(10));
    ARIADNE_TESTING_EXPECT_EQUAL(bin(uint16_t(5), uint16_t(6)), uint16_t(0));
    ARIADNE_TESTING_EXPECT_EXCEPTION(bin(uint16_t(5), uint16_t(7)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(bin(uint16_t(16), uint16_t(1)));

    ARIADNE_TESTING_EXPECT_EQUAL(fac(uint32_t(0)), uint32_t(1));
    ARIADNE_TESTING_EXPECT_EQUAL(fac(uint32_t(5)), uint32_t(120));
    ARIADNE_TESTING_EXPECT_EXCEPTION(fac(uint32_t(13)));
    ARIADNE_TESTING_EXPECT_EQUAL(bin(uint32_t(5), uint32_t(0)), uint32_t(1));
    ARIADNE_TESTING_EXPECT_EQUAL(bin(uint32_t(5), uint32_t(2)), uint32_t(10));
    ARIADNE_TESTING_EXPECT_EQUAL(bin(uint32_t(5), uint32_t(6)), uint32_t(0));
    ARIADNE_TESTING_EXPECT_EXCEPTION(bin(uint32_t(5), uint32_t(7)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(bin(uint32_t(31), uint32_t(1)));

    ARIADNE_TESTING_EXPECT_EQUAL(fac(uint64_t(0)), uint64_t(1));
    ARIADNE_TESTING_EXPECT_EQUAL(fac(uint64_t(5)), uint64_t(120));
    ARIADNE_TESTING_EXPECT_EXCEPTION(fac(uint64_t(21)));
    ARIADNE_TESTING_EXPECT_EQUAL(bin(uint64_t(5), uint64_t(0)), uint64_t(1));
    ARIADNE_TESTING_EXPECT_EQUAL(bin(uint64_t(5), uint64_t(2)), uint64_t(10));
    ARIADNE_TESTING_EXPECT_EQUAL(bin(uint64_t(5), uint64_t(6)), uint64_t(0));
    ARIADNE_TESTING_EXPECT_EXCEPTION(bin(uint64_t(5), uint64_t(7)));
    ARIADNE_TESTING_EXPECT_EXCEPTION(bin(uint64_t(63), uint64_t(1)));
}

void TestInteger::test_comparisons() {
    ARIADNE_TESTING_EXPECT_EQUAL(Integer(3), 3);
    ARIADNE_TESTING_EXPECT_EQUAL(3, Integer(3));

    ARIADNE_TESTING_EXPECT_EQUAL(Integer(2), Integer(2));
    ARIADNE_TESTING_EXPECT_EQUAL(Integer(0), Integer(-0));
    ARIADNE_TESTING_EXPECT_NOT_EQUAL(Integer(2), Integer(-2));
    ARIADNE_TESTING_EXPECT_NOT_EQUAL(Integer(2), Integer(-3));
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(Integer(2), Integer(23));
    ARIADNE_TESTING_EXPECT_LESS_EQUAL(Integer(2), Integer(3));
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(Integer(2), Integer(2));
    ARIADNE_TESTING_EXPECT_GREATER_EQUAL(Integer(2), Integer(-3));
    ARIADNE_TESTING_EXPECT_LESS(Integer(2), Integer(3));
    ARIADNE_TESTING_EXPECT_GREATER(Integer(2), Integer(-3));
}

void TestInteger::test_native() {
    const Integer zero(0);
    const Integer one(1);
    const Integer minus_one(-1);
    ARIADNE_TESTING_PROVE(Integer(String("-0")).literal() == "0");
    ARIADNE_TESTING_PROVE(Integer(String("+00000042")).literal() == "42");
    ARIADNE_TESTING_PROVE(Integer(String("00000")).literal() == "0");
    ARIADNE_TESTING_PROVE(Integer(String("18446744073709551616")).literal() == "18446744073709551616");
    ARIADNE_TESTING_PROVE(Integer(std::numeric_limits<std::int64_t>::min())
                              .get<std::int64_t>() == std::numeric_limits<std::int64_t>::min());
    ARIADNE_TESTING_PROVE(Integer(std::numeric_limits<std::uint64_t>::max())
                              .get<std::uint64_t>() == std::numeric_limits<std::uint64_t>::max());
    ARIADNE_TESTING_PROVE(pow(Integer(0), Nat(0u)) == 1);
    ARIADNE_TESTING_PROVE(pow(Integer(-7), Nat(3u)) == -343);
    ARIADNE_TESTING_PROVE(pow(Integer(-7), Nat(4u)) == 2401);
    ARIADNE_TESTING_PROVE(sqr(Integer(-13)) == 169);
    ARIADNE_TESTING_PROVE(abs(Integer(-13)) == 13);
    ARIADNE_TESTING_PROVE(log2floor(Natural(0u)) == -1);
    ARIADNE_TESTING_PROVE(log2floor(Natural(8u)) == 3);
    ARIADNE_TESTING_PROVE(quot(Integer(-11), Integer(3)) == -3);
    ARIADNE_TESTING_PROVE(rem(Integer(-11), Integer(3)) == -2);
    ARIADNE_TESTING_PROVE(quot(Integer(11), Integer(-3)) == -3);
    ARIADNE_TESTING_PROVE(rem(Integer(11), Integer(-3)) == 2);

    for (const char *invalid: {"", "-", "+", "-+1", " 1", "1 ", "1a", "1.0"}) {
        bool threw = false;
        try {
            (void)Integer(String(invalid));
        } catch (const std::invalid_argument &) {
            threw = true;
        }
        ARIADNE_TESTING_PROVE(threw);
    }
    bool negative_unsigned = false, positive_overflow = false, division_zero = false;
    try {
        (void)minus_one.get<std::uint64_t>();
    } catch (const std::out_of_range &) {
        negative_unsigned = true;
    }
    try {
        (void)Integer(String("9223372036854775808")).get<std::int64_t>();
    } catch (const std::out_of_range &) {
        positive_overflow = true;
    }
    try {
        (void)quot(one, zero);
    } catch (const std::domain_error &) {
        division_zero = true;
    }
    ARIADNE_TESTING_PROVE(negative_unsigned && positive_overflow && division_zero);

    for (const char *a: {"-11", "11", "-12", "12", "0"}) {
        for (const char *b: {"-3", "3", "-5", "5"}) {
            Integer dividend{String(a)};
            Integer divisor{String(b)};
            const Integer quotient = quot(dividend, divisor);
            const Integer remainder = rem(dividend, divisor);
            ARIADNE_TESTING_PROVE(dividend == quotient * divisor + remainder);
            ARIADNE_TESTING_PROVE(abs(remainder) < abs(divisor));
            ARIADNE_TESTING_PROVE(remainder == 0 || sgn(remainder) == sgn(dividend));
        }
    }

    // I41: six structural classes x four sign pairs x eight widths.
    // Source and expected values are independently frozen, with no GMP
    // calls in this differential test. Self-aliasing remains covered.
    std::ifstream sub_fixture(ARIADNE_INTEGER_SUB_ORACLE_PATH);
    if (!sub_fixture)
        throw std::runtime_error("Missing Integer subtraction TSV");
    std::string sub_header;
    while (sub_fixture.peek() == '#')
        std::getline(sub_fixture, sub_header);
    for (const std::size_t width: {1u, 2u, 7u, 8u, 16u, 32u, 64u, 128u}) {
        for (unsigned mode = 0u; mode < 6u; ++mode) {
            for (unsigned signs = 0u; signs < 4u; ++signs) {
                std::size_t ref_width = 0u;
                unsigned ref_mode = 0u, ref_signs = 0u;
                std::string raw_lhs, raw_rhs, raw_expected;
                if (!(sub_fixture >> ref_width >> ref_mode >> ref_signs >>
                      raw_lhs >> raw_rhs >> raw_expected) ||
                    ref_width != width || ref_mode != mode || ref_signs != signs) {
                    throw std::runtime_error("Missing or unordered Integer subtraction TSV row");
                }
                const Integer lhs{String(raw_lhs)}, rhs{String(raw_rhs)};
                const Integer expected_sub{String(raw_expected)};
                ARIADNE_TESTING_EXPECT_EQUAL(lhs - rhs, expected_sub);
                Integer assigned = lhs;
                assigned -= rhs;
                ARIADNE_TESTING_EXPECT_EQUAL(assigned, expected_sub);
                Integer aliased = lhs;
                const Integer &same_object = aliased;
                aliased -= same_object;
                ARIADNE_TESTING_PROVE(aliased == 0);
                ARIADNE_TESTING_PROVE(!aliased.negative());
            }
        }
    }
    std::string sub_extra;
    if (sub_fixture >> sub_extra)
        throw std::runtime_error("Extra Integer subtraction TSV rows");

    // B32: independent frozen exact signed-power reference; retain all
    // sign/width cases and sample exponent bit patterns at boundaries.
    std::ifstream pow_fixture(ARIADNE_INTEGER_POWER_ORACLE_PATH);
    if (!pow_fixture)
        throw std::runtime_error("Missing Integer power TSV");
    std::string pow_header;
    while (pow_fixture.peek() == '#')
        std::getline(pow_fixture, pow_header);
    for (const char *literal:
         {"0", "1", "-1", "2", "-2", "17", "-17",
          "18446744073709551616", "-18446744073709551616",
          "340282366920938463463374607431768211455"}) {
        const Integer input{String(literal)};
        for (unsigned exponent: {0u, 1u, 2u, 3u, 5u, 7u, 8u, 9u, 11u, 12u}) {
            std::string ref_base, ref_result;
            unsigned ref_exponent = 0u;
            if (!(pow_fixture >> ref_base >> ref_exponent >> ref_result) ||
                ref_base != literal || ref_exponent != exponent) {
                throw std::runtime_error("Missing or unordered Integer power TSV row");
            }
            ARIADNE_TESTING_EXPECT_EQUAL(pow(input, Nat(exponent)),
                                         Integer(String(ref_result)));
        }
    }
    std::string pow_extra;
    if (pow_fixture >> pow_extra) {
        throw std::runtime_error("Extra Integer power TSV rows");
    }

    // I66: original 18 exponent boundaries, across each original width.
    // Four distinct input-pattern classes are retained from the source
    // mt19937_64 sequence. Results were frozen using independent exact
    // integer arithmetic (no GMP or Integer-under-test at runtime).
    std::ifstream exponent_oracle(ARIADNE_INTEGER_I66_ORACLE_PATH);
    if (!exponent_oracle)
        throw std::runtime_error("Missing Integer I66 TSV");
    std::string exponent_header;
    while (exponent_oracle.peek() == '#')
        std::getline(exponent_oracle, exponent_header);
    for (const std::size_t width: {0u, 1u, 2u, 3u, 4u}) {
        for (unsigned mode = 0u; mode < 4u; ++mode) {
            for (unsigned exponent:
                 {0u, 1u, 2u, 3u, 5u, 7u, 8u, 9u, 12u, 15u, 16u, 17u,
                  31u, 32u, 33u, 63u, 64u, 65u}) {
                std::size_t ref_width = 0u;
                unsigned ref_mode = 0u, ref_exponent = 0u;
                std::string raw_base, raw_result;
                if (!(exponent_oracle >> ref_width >> ref_mode >> ref_exponent >> raw_base >> raw_result) || ref_width != width || ref_mode != mode || ref_exponent != exponent) {
                    throw std::runtime_error("Missing or unordered I66 TSV row");
                }
                const Integer base{String(raw_base)};
                ARIADNE_TESTING_EXPECT_EQUAL(pow(base, Nat(exponent)),
                                             Integer(String(raw_result)));
            }
        }
    }
    std::string extra_exponent;
    if (exponent_oracle >> extra_exponent)
        throw std::runtime_error("Extra I66 TSV rows");

    // I33: independent exact source-keyed squares for all eight distinct
    // structural classes at every original width; includes x*=x aliasing.
    std::ifstream square_oracle(ARIADNE_INTEGER_I33_ORACLE_PATH);
    if (!square_oracle)
        throw std::runtime_error("Missing Integer I33 TSV");
    std::string square_comment;
    while (square_oracle.peek() == '#')
        std::getline(square_oracle, square_comment);
    for (const std::size_t width: {1u, 7u, 8u, 9u, 16u, 17u, 32u, 64u, 128u}) {
        for (unsigned mode = 0u; mode < 8u; ++mode) {
            std::size_t ref_width = 0u;
            unsigned ref_mode = 0u;
            std::string source, expected;
            if (!(square_oracle >> ref_width >> ref_mode >> source >> expected) || ref_width != width || ref_mode != mode) {
                throw std::runtime_error("Missing or unordered I33 TSV row");
            }
            const Integer x{String(source)}, copy = x;
            const Integer square{String(expected)};
            ARIADNE_TESTING_EXPECT_EQUAL(x * x, square);
            ARIADNE_TESTING_EXPECT_EQUAL(x * copy, square);
            Integer assigned = x;
            assigned *= assigned;
            ARIADNE_TESTING_EXPECT_EQUAL(assigned, square);
            ARIADNE_TESTING_PROVE(!assigned.negative());
        }
    }
    std::string extra_square;
    if (square_oracle >> extra_square)
        throw std::runtime_error("Extra I33 TSV rows");

    // B37: keep the saturation/high-bit/dense sign classes across
    // the 7/8/9-limb square dispatch and all original wider widths.
    std::ifstream large_power_oracle(ARIADNE_INTEGER_I37_ORACLE_PATH);
    if (!large_power_oracle)
        throw std::runtime_error("Missing Integer B37 TSV");
    std::string large_power_comment;
    while (large_power_oracle.peek() == '#')
        std::getline(large_power_oracle, large_power_comment);
    for (const std::size_t width: {7u, 8u, 9u, 16u, 32u, 64u}) {
        for (unsigned repetition = 0u; repetition < 4u; ++repetition) {
            for (unsigned exponent: {2u, 3u, 5u}) {
                std::size_t ref_width = 0u;
                unsigned ref_repetition = 0u, ref_exponent = 0u;
                std::string source, expected;
                if (!(large_power_oracle >> ref_width >> ref_repetition >> ref_exponent >> source >> expected) || ref_width != width || ref_repetition != repetition || ref_exponent != exponent) {
                    throw std::runtime_error("Missing or unordered B37 TSV row");
                }
                const Integer x{String(source)};
                ARIADNE_TESTING_EXPECT_EQUAL(pow(x, Nat(exponent)),
                                             Integer(String(expected)));
            }
        }
    }
    std::string extra_large_power;
    if (large_power_oracle >> extra_large_power)
        throw std::runtime_error("Extra B37 TSV rows");

    // B38: all six sparse/dense structural cases at 8..128 limbs,
    // including the alternating-zero, saturated and both sign classes.
    std::ifstream cubic_oracle(ARIADNE_INTEGER_I38_ORACLE_PATH);
    if (!cubic_oracle)
        throw std::runtime_error("Missing Integer B38 TSV");
    std::string cubic_comment;
    while (cubic_oracle.peek() == '#')
        std::getline(cubic_oracle, cubic_comment);
    for (const std::size_t width: {8u, 16u, 32u, 64u, 128u}) {
        for (unsigned repetition = 0u; repetition < 6u; ++repetition) {
            std::size_t ref_width = 0u;
            unsigned ref_repetition = 0u;
            std::string source, expected;
            if (!(cubic_oracle >> ref_width >> ref_repetition >> source >> expected) || ref_width != width || ref_repetition != repetition) {
                throw std::runtime_error("Missing or unordered B38 TSV row");
            }
            ARIADNE_TESTING_EXPECT_EQUAL(pow(Integer(String(source)), Nat(3u)),
                                         Integer(String(expected)));
        }
    }
    std::string extra_cubic;
    if (cubic_oracle >> extra_cubic)
        throw std::runtime_error("Extra B38 TSV rows");

    // I18: frozen exact signed arithmetic on 170 strategically sampled
    // originals out of 1700; expected results never use GMP at runtime.
    std::ifstream golden(ARIADNE_INTEGER_I18_ORACLE_PATH);
    if (!golden)
        throw std::runtime_error("Missing Integer I18 oracle TSV");
    std::string comment;
    while (golden.peek() == '#')
        std::getline(golden, comment);
    for (unsigned sample = 0u; sample < 170u; ++sample) {
        unsigned rep = 0u;
        std::string lhs, rhs, sum, sub, mul, q_expected, r_expected;
        if (!(golden >> rep >> lhs >> rhs >> sum >> sub >> mul >> q_expected >> r_expected) ||
            rep != 10u * sample + (7u * sample + 3u) % 10u) {
            throw std::runtime_error("Missing or out-of-order Integer I18 oracle");
        }
        const Integer x{String(lhs)}, y{String(rhs)};
        ARIADNE_TESTING_EXPECT_EQUAL(x + y, Integer(String(sum)));
        ARIADNE_TESTING_EXPECT_EQUAL(x - y, Integer(String(sub)));
        ARIADNE_TESTING_EXPECT_EQUAL(x * y, Integer(String(mul)));
        const Integer q = quot(x, y), r = rem(x, y);
        ARIADNE_TESTING_EXPECT_EQUAL(q, Integer(String(q_expected)));
        ARIADNE_TESTING_EXPECT_EQUAL(r, Integer(String(r_expected)));
        ARIADNE_TESTING_PROVE(q * y + r == x);
        ARIADNE_TESTING_PROVE(abs(r) < abs(y));
        ARIADNE_TESTING_PROVE(r == 0 || sgn(r) == sgn(x));
        ARIADNE_TESTING_PROVE(x.literal() == Integer(String(lhs)).literal());
    }
    std::string extra;
    if (golden >> extra)
        throw std::runtime_error("Extra Integer I18 oracle rows");
}

// Preserve all I17 signed-magnitude regression scenarios on the
// production Integer type. The earlier I18 differential corpus (1700
// independent cases) is kept in test_native() above.
void TestInteger::test_signed_regression() {
    const Integer zero(0);
    const Integer one(1);
    const Integer minus_one(-1);
    ARIADNE_TESTING_PROVE(sgn(zero) == Sign::ZERO);
    ARIADNE_TESTING_PROVE(!Integer(String("-0")).negative());
    ARIADNE_TESTING_PROVE(Integer(String("-00000")).literal() == "0");
    ARIADNE_TESTING_PROVE(Integer(String("+000123")).literal() == "123");
    ARIADNE_TESTING_PROVE(Integer::from_magnitude(BigUnsigned{}, true) == zero);
    ARIADNE_TESTING_PROVE(!Integer::from_magnitude(BigUnsigned{}, true).negative());
    ARIADNE_TESTING_PROVE(cmp(one, minus_one) == Comparison::GREATER);
    ARIADNE_TESTING_PROVE(cmp(minus_one, zero) == Comparison::LESS);
    ARIADNE_TESTING_PROVE(cmp(zero, zero) == Comparison::EQUAL);
    ARIADNE_TESTING_PROVE(Integer(std::numeric_limits<std::int64_t>::min())
                              .get<std::int64_t>() == std::numeric_limits<std::int64_t>::min());
    ARIADNE_TESTING_PROVE(Integer(std::numeric_limits<std::uint64_t>::max())
                              .get<std::uint64_t>() == std::numeric_limits<std::uint64_t>::max());
    ARIADNE_TESTING_PROVE(pow(zero, Nat(0u)).literal() == "1");
    ARIADNE_TESTING_PROVE(sqr(Integer(-13)) == 169);
    ARIADNE_TESTING_PROVE(quot(Integer(-7), Integer(2)) == -3);

    for (const char *invalid: {"", "-", "+", " 1", "1 ", "1_0", "--1"}) {
        bool threw = false;
        try {
            (void)Integer(String(invalid));
        } catch (const std::invalid_argument &) {
            threw = true;
        }
        ARIADNE_TESTING_PROVE(threw);
    }
    bool negative_to_unsigned = false, positive_overflow = false;
    try {
        (void)minus_one.get<std::uint64_t>();
    } catch (const std::out_of_range &) {
        negative_to_unsigned = true;
    }
    try {
        (void)Integer(String("9223372036854775808"))
            .get<std::int64_t>();
    } catch (const std::out_of_range &) {
        positive_overflow = true;
    }
    ARIADNE_TESTING_PROVE(negative_to_unsigned && positive_overflow);

    for (const char *left: {"-11", "11", "-12", "12", "0"}) {
        for (const char *right: {"-3", "3", "-5", "5"}) {
            const Integer a{String(left)}, b{String(right)};
            const Integer q = quot(a, b), r = rem(a, b);
            ARIADNE_TESTING_PROVE(cmp(q * b + r, a) == Comparison::EQUAL);
            ARIADNE_TESTING_PROVE(abs(r) < abs(b));
            ARIADNE_TESTING_PROVE(is_zero(r) || sgn(r) == sgn(a));
        }
    }

    // I17: 250 source-keyed references, one per ten original visits.
    // All 130 shift counts survive along with signed and zero cases.
    std::ifstream golden(ARIADNE_INTEGER_I17_ORACLE_PATH);
    if (!golden)
        throw std::runtime_error("Missing Integer I17 oracle TSV");
    std::string comment;
    while (golden.peek() == '#')
        std::getline(golden, comment);
    const auto check = [&](const Integer &actual, const std::string &reference) {
        ARIADNE_TESTING_EXPECT_EQUAL(actual, Integer(String(reference)));
        ARIADNE_TESTING_PROVE(!is_zero(actual) || !actual.negative());
        ARIADNE_TESTING_PROVE((reference == "0") == is_zero(actual));
        ARIADNE_TESTING_PROVE((reference.front() == '-') == actual.negative());
    };
    for (unsigned sample = 0u; sample < 250u; ++sample) {
        unsigned rep = 0u;
        std::string ad, bd, sum, sub, mul, q_expected, r_expected, shr, shl;
        if (!(golden >> rep >> ad >> bd >> sum >> sub >> mul >>
              q_expected >> r_expected >> shr >> shl) ||
            rep != 10u * sample + (7u * sample + 3u) % 10u) {
            throw std::runtime_error("Missing or unordered Integer I17 oracle");
        }
        const Integer x{String(ad)}, y{String(bd)};
        check(x, ad);
        check(y, bd);
        check(x + y, sum);
        check(x - y, sub);
        check(x * y, mul);
        const Integer q = quot(x, y), r = rem(x, y);
        check(q, q_expected);
        check(r, r_expected);
        ARIADNE_TESTING_PROVE(q * y + r == x);
        ARIADNE_TESTING_PROVE(abs(r) < abs(y));
        ARIADNE_TESTING_PROVE(is_zero(r) || sgn(r) == sgn(x));

        const unsigned bits = rep % 130u;
        const BigUnsigned scale = shift_left(
            BigUnsigned(std::uint64_t{1}), bits);
        const Integer divisor = Integer::from_magnitude(scale, false);
        const Integer shifted_right = Integer::from_magnitude(
            shift_right(x.magnitude(), bits), x.negative());
        check(shifted_right, shr);
        check(quot(x, divisor), shr);
        const Integer shifted_left = Integer::from_magnitude(
            shift_left(x.magnitude(), bits), x.negative());
        check(shifted_left, shl);
    }
    std::string extra;
    if (golden >> extra)
        throw std::runtime_error("Extra Integer I17 oracle rows");

    bool quotient_zero = false, remainder_zero = false;
    try {
        (void)quot(one, zero);
    } catch (const std::domain_error &) {
        quotient_zero = true;
    }
    try {
        (void)rem(one, zero);
    } catch (const std::domain_error &) {
        remainder_zero = true;
    }
    ARIADNE_TESTING_PROVE(quotient_zero && remainder_zero);
}

int main(int argc, char **argv) {
    Ariadne::Testing::Suite suite("integer");
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestInteger, test_constructors);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestInteger, test_comparisons);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestInteger, test_literal);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestInteger, test_arithmetic);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestInteger, test_small_types);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestInteger, test_combinatorics);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestInteger, test_native);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestInteger, test_signed_regression);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestInteger, test_owned_results);
    ARIADNE_TESTING_ADD_TEST_TO(suite, TestInteger, test_sqrt);
    return Ariadne::Testing::run_suite(suite, argc, argv);
}
