/***************************************************************************
 *            benchmark_exact_arithmetic_workload.hpp
 *
 *  Copyright  2026  Ariadne contributors
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

/*! \file benchmark_exact_arithmetic_workload.hpp
 *  \brief Shared direct-GMP and native integer/rational benchmarks.
 */

#ifndef ARIADNE_NUMERIC_BENCHMARK_EXACT_ARITHMETIC_WORKLOAD_HPP
#define ARIADNE_NUMERIC_BENCHMARK_EXACT_ARITHMETIC_WORKLOAD_HPP

#include "numeric/integer.hpp"
#include "numeric/rational.hpp"

#include <gmp.h>

#include <algorithm>
#include <bit>
#include <chrono>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace AriadneBenchmark {

using namespace Ariadne;

// RAII controls only the lifetime of raw GMP objects. Each operation
// below calls mpz_* or mpq_* directly, not the old Ariadne wrappers.
struct RawInteger {
    mpz_t value;

    RawInteger() { mpz_init(value); }

    explicit RawInteger(std::string const& decimal) : RawInteger() {
        if (mpz_set_str(value,decimal.c_str(),10)!=0) {
            throw std::invalid_argument("invalid GMP integer input");
        }
    }

    explicit RawInteger(mpz_srcptr source) : RawInteger() {
        mpz_set(value,source);
    }

    RawInteger(RawInteger const& other) : RawInteger(other.value) { }

    RawInteger(RawInteger&& other) noexcept : RawInteger() {
        mpz_swap(value,other.value);
    }

    RawInteger& operator=(RawInteger const& other) {
        if (this!=&other) { mpz_set(value,other.value); }
        return *this;
    }

    RawInteger& operator=(RawInteger&& other) noexcept {
        if (this!=&other) { mpz_swap(value,other.value); }
        return *this;
    }

    ~RawInteger() { mpz_clear(value); }
};

struct RawRational {
    mpq_t value;

    RawRational() { mpq_init(value); }

    RawRational(RawInteger const& numerator, RawInteger const& denominator)
        : RawRational()
    {
        mpq_set_num(value,numerator.value);
        mpq_set_den(value,denominator.value);
        mpq_canonicalize(value);
    }

    RawRational(RawRational const& other) : RawRational() {
        mpq_set(value,other.value);
    }

    RawRational(RawRational&& other) noexcept : RawRational() {
        mpq_swap(value,other.value);
    }

    RawRational& operator=(RawRational const& other) {
        if (this!=&other) { mpq_set(value,other.value); }
        return *this;
    }

    RawRational& operator=(RawRational&& other) noexcept {
        if (this!=&other) { mpq_swap(value,other.value); }
        return *this;
    }

    ~RawRational() { mpq_clear(value); }
};

inline std::string decimal(mpz_srcptr value) {
    std::vector<char> buffer(mpz_sizeinbase(value,10)+3u);
    mpz_get_str(buffer.data(),10,value);
    return std::string(buffer.data());
}

inline std::string exact_value(Integer const& value) {
    return std::string(value.literal());
}

inline std::string exact_value(RawInteger const& value) {
    return decimal(value.value);
}

inline std::string exact_value(Rational const& value) {
    return exact_value(value.numerator()) + "/" + exact_value(value.denominator());
}

inline std::string exact_value(RawRational const& value) {
    return decimal(mpq_numref(value.value)) + "/"
         + decimal(mpq_denref(value.value));
}

inline std::string exact_value(int value) {
    return std::to_string(value);
}

inline std::string exact_value(double value) {
    // Match binary64 conversions exactly, including the sign of zero.
    const std::uint64_t bits=std::bit_cast<std::uint64_t>(value);
    return std::to_string(bits);
}

template<class Value> void keep_result_alive(Value const& value) {
#if defined(__GNUC__) || defined(__clang__)
    asm volatile("" : : "g"(&value) : "memory");
#else
    static void const* volatile sink=nullptr;
    sink=&value;
#endif
}

template<class Operation>
double exact_arithmetic_median_ns(Operation&& operation, unsigned iterations,
                                   unsigned samples)
{
    std::vector<double> measurements;
    measurements.reserve(samples);
    for (unsigned sample=0; sample<=samples; ++sample) {
        auto const begin=std::chrono::steady_clock::now();
        for (unsigned i=0; i<iterations; ++i) {
            auto result=operation(i);
            keep_result_alive(result);
        }
        auto const end=std::chrono::steady_clock::now();
        if (sample!=0u) {
            measurements.push_back(
                std::chrono::duration<double,std::nano>(end-begin).count()
                / static_cast<double>(iterations));
        }
    }
    std::sort(measurements.begin(),measurements.end());
    return measurements[samples/2u];
}

struct ExactArithmeticDatum {
    std::string backend;
    std::string operation;
    unsigned bits;
    double nanoseconds;
};

// Numeric value specifications are generated once from exact GMP integers.
// Both backends receive the same decimal strings. No GMP generator work
// is performed within a timed operation except explicit construction tests.
struct IntegerSpec {
    std::string a;
    std::string b;
};

inline std::vector<IntegerSpec> integer_specifications(unsigned bits) {
    std::vector<IntegerSpec> entries;
    entries.reserve(32u);
    RawInteger base, left, right;
    mpz_set_ui(base.value,1ul);
    mpz_mul_2exp(base.value,base.value,static_cast<mp_bitcnt_t>(bits-1u));
    for (unsigned i=0; i<32u; ++i) {
        mpz_add_ui(left.value,base.value,13ul+2ul*i);
        mpz_add_ui(right.value,base.value,31ul+5ul*i);
        if (i%4u==0u || i%4u==3u) { mpz_neg(left.value,left.value); }
        if (i%4u==1u || i%4u==3u) { mpz_neg(right.value,right.value); }
        entries.push_back({decimal(left.value),decimal(right.value)});
    }
    return entries;
}

// Additional division-only operands; preserve the original historical pool.
inline std::vector<IntegerSpec> integer_division_specifications(unsigned bits) {
    std::vector<IntegerSpec> entries;
    entries.reserve(32u);
    RawInteger a,b,base;
    mpz_set_ui(base.value,1ul);
    mpz_mul_2exp(base.value,base.value,static_cast<mp_bitcnt_t>(bits-1u));
    for (unsigned i=0u;i<32u;++i) {
        mpz_add_ui(a.value,base.value,101ul+2ul*i);
        if (i%4u==0u) {
            mpz_set_ui(b.value,3ul+i);
        } else if (i%4u==1u) {
            mpz_set_ui(b.value,1ul);
            mpz_mul_2exp(b.value,b.value,static_cast<mp_bitcnt_t>(bits/2u));
            mpz_add_ui(b.value,b.value,1ul+i);
        } else if (i%4u==2u) {
            mpz_set(b.value,a.value);
            mpz_sub_ui(b.value,b.value,1ul+i);
        } else {
            mpz_set_ui(b.value,2ul+i);
        }
        if (i%4u==2u || i%4u==3u) mpz_neg(a.value,a.value);
        if (i%4u==1u || i%4u==3u) mpz_neg(b.value,b.value);
        entries.push_back({decimal(a.value),decimal(b.value)});
    }
    return entries;
}

struct RationalSpec {
    std::string numerator_a;
    std::string denominator_a;
    std::string numerator_b;
    std::string denominator_b;
};

inline std::vector<RationalSpec> rational_specifications(unsigned bits) {
    std::vector<RationalSpec> entries;
    entries.reserve(32u);
    RawInteger base, an, ad, bn, bd;
    mpz_set_ui(base.value,1ul);
    mpz_mul_2exp(base.value,base.value,static_cast<mp_bitcnt_t>(bits));
    for (unsigned i=0; i<32u; ++i) {
        mpz_add_ui(an.value,base.value,13ul+2ul*i);
        mpz_mul_ui(ad.value,base.value,2ul);
        mpz_add_ui(ad.value,ad.value,17ul+3ul*i);
        mpz_add_ui(bn.value,base.value,31ul+5ul*i);
        mpz_mul_ui(bd.value,base.value,3ul);
        mpz_add_ui(bd.value,bd.value,19ul+7ul*i);
        if (i%4u==0u || i%4u==3u) { mpz_neg(an.value,an.value); }
        if (i%4u==1u || i%4u==3u) { mpz_neg(bn.value,bn.value); }
        // Include nontrivial common factors and reduction costs.
        if (i%4u==2u) {
            mpz_mul_ui(an.value,an.value,7ul);
            mpz_mul_ui(ad.value,ad.value,7ul);
            mpz_mul_ui(bn.value,bn.value,11ul);
            mpz_mul_ui(bd.value,bd.value,11ul);
        }
        entries.push_back({decimal(an.value),decimal(ad.value),
                           decimal(bn.value),decimal(bd.value)});
    }
    return entries;
}

struct IntegerOperations {
    using Number=Integer;

    static Number parse(std::string const& value) {
        return Integer(String(value));
    }
    static Number copy(Number const& a) { return Number(a); }
    static int compare(Number const& a,Number const& b) {
        const auto order=cmp(a,b);
        return order==Comparison::LESS ? -1 : order==Comparison::GREATER ? 1 : 0;
    }
    static Number negate(Number const& a) { return -a; }
    static Number add(Number const& a,Number const& b) { return a+b; }
    static Number sub(Number const& a,Number const& b) { return a-b; }
    static Number mul(Number const& a,Number const& b) { return a*b; }
    static Number square(Number const& a) { return a*a; }
    static Number quotient(Number const& a,Number const& b) {
        return quot(a,b);
    }
    static Number remainder(Number const& a,Number const& b) {
        return rem(a,b);
    }
    static Number power3(Number const& a) { return pow(a,Nat(3u)); }
    static std::string to_decimal(Number const& a) {
        return std::string(a.literal());
    }
};

struct GmpIntegerOperations {
    using Number=RawInteger;

    static Number parse(std::string const& value) { return Number(value); }
    static Number copy(Number const& a) { return Number(a); }
    static int compare(Number const& a,Number const& b) {
        const int c=mpz_cmp(a.value,b.value);
        return c<0 ? -1 : c>0 ? 1 : 0;
    }
    static Number negate(Number const& a) {
        Number result;
        mpz_neg(result.value,a.value);
        return result;
    }
    static Number add(Number const& a,Number const& b) {
        Number result;
        mpz_add(result.value,a.value,b.value);
        return result;
    }
    static Number sub(Number const& a,Number const& b) {
        Number result;
        mpz_sub(result.value,a.value,b.value);
        return result;
    }
    static Number mul(Number const& a,Number const& b) {
        Number result;
        mpz_mul(result.value,a.value,b.value);
        return result;
    }
    static Number square(Number const& a) {
        Number result;
        mpz_mul(result.value,a.value,a.value);
        return result;
    }
    static Number quotient(Number const& a,Number const& b) {
        Number result;
        // Native Integer uses truncation toward zero for quotient.
        mpz_tdiv_q(result.value,a.value,b.value);
        return result;
    }
    static Number remainder(Number const& a,Number const& b) {
        Number result;
        // The remainder has the dividend sign, matching Integer.
        mpz_tdiv_r(result.value,a.value,b.value);
        return result;
    }
    static Number power3(Number const& a) {
        Number result;
        mpz_pow_ui(result.value,a.value,3ul);
        return result;
    }
    static std::string to_decimal(Number const& a) {
        return decimal(a.value);
    }
};

template<class Operations>
struct IntegerPool {
    using Number=typename Operations::Number;
    static constexpr unsigned pool_size=32u;

    std::vector<IntegerSpec> spec,division_spec;
    std::vector<Number> a,b,division_a,division_b;

    explicit IntegerPool(unsigned bits)
        : spec(integer_specifications(bits)),
          division_spec(integer_division_specifications(bits)) {
        a.reserve(pool_size);
        b.reserve(pool_size);
        division_a.reserve(pool_size);
        division_b.reserve(pool_size);
        for (auto const& item : spec) {
            a.push_back(Operations::parse(item.a));
            b.push_back(Operations::parse(item.b));
        }
        for (auto const& item : division_spec) {
            division_a.push_back(Operations::parse(item.a));
            division_b.push_back(Operations::parse(item.b));
        }
    }
};

template<class Operations,class Visit>
void integer_operations(IntegerPool<Operations> const& pool,Visit&& report) {
    constexpr unsigned size=IntegerPool<Operations>::pool_size;
    auto index=[](unsigned i) { return i%size; };
    auto const& a=pool.a;
    auto const& b=pool.b;
    report("construct_decimal",[&](unsigned i) {
        return Operations::parse(pool.spec[index(i)].a);
    });
    report("copy",[&](unsigned i) {
        return Operations::copy(a[index(i)]);
    });
    report("compare",[&](unsigned i) {
        return Operations::compare(a[index(i)],b[index(i)]);
    });
    report("negate",[&](unsigned i) {
        return Operations::negate(a[index(i)]);
    });
    report("add",[&](unsigned i) {
        return Operations::add(a[index(i)],b[index(i)]);
    });
    report("sub",[&](unsigned i) {
        return Operations::sub(a[index(i)],b[index(i)]);
    });
    report("mul",[&](unsigned i) {
        return Operations::mul(a[index(i)],b[index(i)]);
    });
    report("square",[&](unsigned i) {
        return Operations::square(a[index(i)]);
    });
    report("quotient",[&](unsigned i) {
        return Operations::quotient(a[index(i)],b[index(i)]);
    });
    report("remainder",[&](unsigned i) {
        return Operations::remainder(a[index(i)],b[index(i)]);
    });
    report("quotient_nontrivial",[&](unsigned i) {
        return Operations::quotient(pool.division_a[index(i)],
                                    pool.division_b[index(i)]);
    });
    report("remainder_nontrivial",[&](unsigned i) {
        return Operations::remainder(pool.division_a[index(i)],
                                     pool.division_b[index(i)]);
    });
    report("pow3",[&](unsigned i) {
        return Operations::power3(a[index(i)]);
    });
    report("to_decimal",[&](unsigned i) {
        return Operations::to_decimal(a[index(i)]);
    });
}

struct RationalOperations {
    using Number=Rational;
    using IntegerType=Integer;

    static IntegerType parse_integer(std::string const& value) {
        return IntegerType(String(value));
    }
    static Number construct(IntegerType const& n,IntegerType const& d) {
        return Rational(n,d);
    }
    static Number copy(Number const& a) { return Number(a); }
    static int compare(Number const& a,Number const& b) {
        return a<b ? -1 : b<a ? 1 : 0;
    }
    static Number negate(Number const& a) { return -a; }
    static Number add(Number const& a,Number const& b) { return a+b; }
    static Number sub(Number const& a,Number const& b) { return a-b; }
    static Number mul(Number const& a,Number const& b) { return a*b; }
    static Number divide(Number const& a,Number const& b) { return a/b; }
    static Number square(Number const& a) { return a*a; }
    static Number reciprocal(Number const& a) { return rec(a); }
    static IntegerType numerator(Number const& a) { return a.get_num(); }
    static IntegerType denominator(Number const& a) { return a.get_den(); }
    static IntegerType floor_value(Number const& a) { return floor(a); }
    static IntegerType ceil_value(Number const& a) { return ceil(a); }
    static double to_double(Number const& a) { return a.get_d(); }
};

struct GmpRationalOperations {
    using Number=RawRational;
    using IntegerType=RawInteger;

    static IntegerType parse_integer(std::string const& value) {
        return IntegerType(value);
    }
    static Number construct(IntegerType const& n,IntegerType const& d) {
        return Number(n,d);
    }
    static Number copy(Number const& a) { return Number(a); }
    static int compare(Number const& a,Number const& b) {
        const int c=mpq_cmp(a.value,b.value);
        return c<0 ? -1 : c>0 ? 1 : 0;
    }
    static Number negate(Number const& a) {
        Number result;
        mpq_neg(result.value,a.value);
        return result;
    }
    static Number add(Number const& a,Number const& b) {
        Number result;
        mpq_add(result.value,a.value,b.value);
        return result;
    }
    static Number sub(Number const& a,Number const& b) {
        Number result;
        mpq_sub(result.value,a.value,b.value);
        return result;
    }
    static Number mul(Number const& a,Number const& b) {
        Number result;
        mpq_mul(result.value,a.value,b.value);
        return result;
    }
    static Number divide(Number const& a,Number const& b) {
        Number result;
        mpq_div(result.value,a.value,b.value);
        return result;
    }
    static Number square(Number const& a) {
        Number result;
        mpq_mul(result.value,a.value,a.value);
        return result;
    }
    static Number reciprocal(Number const& a) {
        Number result;
        mpq_inv(result.value,a.value);
        return result;
    }
    static IntegerType numerator(Number const& a) {
        return IntegerType(mpq_numref(a.value));
    }
    static IntegerType denominator(Number const& a) {
        return IntegerType(mpq_denref(a.value));
    }
    static IntegerType floor_value(Number const& a) {
        IntegerType result;
        mpz_fdiv_q(result.value,mpq_numref(a.value),mpq_denref(a.value));
        return result;
    }
    static IntegerType ceil_value(Number const& a) {
        IntegerType result;
        mpz_cdiv_q(result.value,mpq_numref(a.value),mpq_denref(a.value));
        return result;
    }
    static double to_double(Number const& a) { return mpq_get_d(a.value); }
};

template<class Operations>
struct RationalPool {
    using Number=typename Operations::Number;
    using IntegerType=typename Operations::IntegerType;
    static constexpr unsigned pool_size=32u;

    std::vector<RationalSpec> spec;
    std::vector<IntegerType> numerator_a,denominator_a,numerator_b,denominator_b;
    std::vector<Number> a,b;

    explicit RationalPool(unsigned bits) : spec(rational_specifications(bits)) {
        numerator_a.reserve(pool_size);
        denominator_a.reserve(pool_size);
        numerator_b.reserve(pool_size);
        denominator_b.reserve(pool_size);
        a.reserve(pool_size);
        b.reserve(pool_size);
        for (auto const& item : spec) {
            numerator_a.push_back(Operations::parse_integer(item.numerator_a));
            denominator_a.push_back(Operations::parse_integer(item.denominator_a));
            numerator_b.push_back(Operations::parse_integer(item.numerator_b));
            denominator_b.push_back(Operations::parse_integer(item.denominator_b));
            a.push_back(Operations::construct(numerator_a.back(),denominator_a.back()));
            b.push_back(Operations::construct(numerator_b.back(),denominator_b.back()));
        }
    }
};

template<class Operations,class Visit>
void rational_operations(RationalPool<Operations> const& pool,Visit&& report) {
    constexpr unsigned size=RationalPool<Operations>::pool_size;
    auto index=[](unsigned i) { return i%size; };
    auto const& a=pool.a;
    auto const& b=pool.b;
    report("construct_fraction",[&](unsigned i) {
        const unsigned j=index(i);
        return Operations::construct(pool.numerator_a[j],pool.denominator_a[j]);
    });
    report("copy",[&](unsigned i) {
        return Operations::copy(a[index(i)]);
    });
    report("compare",[&](unsigned i) {
        return Operations::compare(a[index(i)],b[index(i)]);
    });
    report("negate",[&](unsigned i) {
        return Operations::negate(a[index(i)]);
    });
    report("add",[&](unsigned i) {
        return Operations::add(a[index(i)],b[index(i)]);
    });
    report("sub",[&](unsigned i) {
        return Operations::sub(a[index(i)],b[index(i)]);
    });
    report("mul",[&](unsigned i) {
        return Operations::mul(a[index(i)],b[index(i)]);
    });
    report("div",[&](unsigned i) {
        return Operations::divide(a[index(i)],b[index(i)]);
    });
    report("square",[&](unsigned i) {
        return Operations::square(a[index(i)]);
    });
    report("reciprocal",[&](unsigned i) {
        return Operations::reciprocal(a[index(i)]);
    });
    report("numerator",[&](unsigned i) {
        return Operations::numerator(a[index(i)]);
    });
    report("denominator",[&](unsigned i) {
        return Operations::denominator(a[index(i)]);
    });
    report("floor",[&](unsigned i) {
        return Operations::floor_value(a[index(i)]);
    });
    report("ceil",[&](unsigned i) {
        return Operations::ceil_value(a[index(i)]);
    });
    report("to_double",[&](unsigned i) {
        return Operations::to_double(a[index(i)]);
    });
}

enum class ExactArithmeticKind { Integer, Rational };

template<ExactArithmeticKind Kind,class Operations>
auto make_exact_pool(unsigned bits) {
    if constexpr (Kind==ExactArithmeticKind::Integer) {
        return IntegerPool<Operations>(bits);
    } else {
        return RationalPool<Operations>(bits);
    }
}

template<ExactArithmeticKind Kind,class Operations,class Pool,class Visit>
void for_each_exact_operation(Pool const& pool,Visit&& report) {
    if constexpr (Kind==ExactArithmeticKind::Integer) {
        integer_operations(pool,std::forward<Visit>(report));
    } else {
        rational_operations(pool,std::forward<Visit>(report));
    }
}

} // namespace AriadneBenchmark

#endif // ARIADNE_NUMERIC_BENCHMARK_EXACT_ARITHMETIC_WORKLOAD_HPP
