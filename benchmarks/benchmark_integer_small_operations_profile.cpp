/***************************************************************************
 *            benchmark_integer_small_operations_profile.cpp
 *
 *  Copyright 2026 Ariadne contributors
 ***************************************************************************/

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

/*! \file benchmark_integer_small_operations_profile.cpp
 *  \brief B67 owning Integer versus owning magnitude primitive diagnostics.
 *
 *  Diagnostic only; the original 84-operation official benchmark remains
 *  completely unchanged. Signed Integer/GMP rows have identical semantics.
 *  "abs core" rows use the corresponding absolute magnitudes, omit the
 *  signed wrapper, and are NOT an additive decomposition of either row.
 *  Unlike a benchmark-specific shortcut, all paths run production APIs.
 */

#include "benchmark_exact_arithmetic_report.hpp"

#include <cstddef>
#include <array>
#include <cstdint>
#include <limits>
#include <random>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace {

using namespace Ariadne;
using namespace AriadneBenchmark;

void check_magnitude(BigUnsigned const& native,mpz_srcptr expected,
                     std::string const& context)
{
    RawInteger actual;
    const auto limbs=native.limbs();
    mpz_import(actual.value,limbs.size(),-1,
               sizeof(BigUnsigned::Limb),0,0,limbs.data());
    if (mpz_cmp(actual.value,expected)!=0) {
        throw std::runtime_error("B67 GMP unsigned oracle mismatch: "+context);
    }
}

std::string_view abs_digits(std::string const& signed_decimal)
{
    std::string_view view(signed_decimal);
    if (!view.empty() && (view.front()=='-' || view.front()=='+')) {
        view.remove_prefix(1u);
    }
    return view;
}

// B68 candidates preserve the COMPLETE signed value/result contract.
// Only String adapter and sign-handling steps differ from production.
Integer parse_direct_signed_diagnostic(std::string const& source)
{
    std::string_view digits(source);
    if (digits.empty()) {
        throw std::invalid_argument("Integer empty decimal literal");
    }
    const bool negative=digits.front()=='-';
    if (negative || digits.front()=='+') {
        digits.remove_prefix(1u);
    }
    return Integer::from_magnitude(from_decimal(digits),negative);
}

std::string format_direct_signed_diagnostic(Integer const& value)
{
    std::string digits=to_decimal(value.magnitude());
    if (value.negative()) {
        digits.insert(digits.begin(),'-');
    }
    return digits;
}

// This is semantically identical to Integer::mul for nonaliased inputs;
// unlike unsigned B67 core rows, the result retains signed semantics.
Integer multiply_direct_signed_diagnostic(Integer const& first,
                                          Integer const& second)
{
    return Integer::from_magnitude(
        mul(first.magnitude(),second.magnitude()),
        first.negative()!=second.negative());
}

Integer square_direct_signed_diagnostic(Integer const& input)
{
    return Integer::from_magnitude(square(input.magnitude()),false);
}

void verify_signed_diagnostics(unsigned bits,
                               IntegerPool<IntegerOperations> const& native,
                               IntegerPool<GmpIntegerOperations> const& gmp)
{
    RawInteger expected,got;
    auto check=[&](Integer const& value,mpz_srcptr oracle,
                   std::string const& context) {
        const auto words=value.magnitude().limbs();
        mpz_import(got.value,words.size(),-1,
                   sizeof(BigUnsigned::Limb),0,0,words.data());
        if (value.negative()) {
            mpz_neg(got.value,got.value);
        }
        if (mpz_cmp(got.value,oracle)!=0) {
            throw std::runtime_error("B68 signed GMP mismatch: "+context);
        }
    };

    for (unsigned i=0u;i<32u;++i) {
        const std::string context=std::to_string(bits)+" bits #"+
                                  std::to_string(i);
        check(parse_direct_signed_diagnostic(native.spec[i].a),
              gmp.a[i].value,context+" parse a");
        check(parse_direct_signed_diagnostic(native.spec[i].b),
              gmp.b[i].value,context+" parse b");
        if (format_direct_signed_diagnostic(native.a[i])!=
                decimal(gmp.a[i].value) ||
                format_direct_signed_diagnostic(native.b[i])!=
                decimal(gmp.b[i].value)) {
            throw std::runtime_error("B68 signed decimal output differs: "+
                                     context);
        }

        mpz_mul(expected.value,gmp.a[i].value,gmp.b[i].value);
        check(multiply_direct_signed_diagnostic(native.a[i],native.b[i]),
              expected.value,context+" multiply");
        mpz_mul(expected.value,gmp.a[i].value,gmp.a[i].value);
        check(square_direct_signed_diagnostic(native.a[i]),
              expected.value,context+" square");
    }

    // The original pool contains negative inputs but not negative zero.
    // Preserve normalization and handling of leading signs explicitly.
    for (const std::string text:{"0","-0","+0","1","-1",
                                  "18446744073709551616",
                                  "-18446744073709551616"}) {
        const Integer ref{String(text)};
        const Integer candidate=parse_direct_signed_diagnostic(text);
        if (candidate.negative()!=ref.negative() ||
                format_direct_signed_diagnostic(candidate)!=
                std::string(ref.literal())) {
            throw std::runtime_error("B68 sign/zero normalization mismatch");
        }
    }
}

template<class Reference,class Production,class Candidate>
void report_signed_candidate(char const* name,
                             Reference&& reference,Production&& production,
                             Candidate&& candidate,
                             unsigned iterations,unsigned samples)
{
    const double g=exact_arithmetic_median_ns(
        std::forward<Reference>(reference),iterations,samples);
    const double n=exact_arithmetic_median_ns(
        std::forward<Production>(production),iterations,samples);
    const double c=exact_arithmetic_median_ns(
        std::forward<Candidate>(candidate),iterations,samples);
    std::ostringstream line;
    line << std::left << std::setw(22) << name
         << std::fixed << std::setprecision(3)
         << std::right << std::setw(12) << g
         << std::setw(14) << n
         << std::setw(14) << c
         << std::setw(14) << c/n << "x"
         << std::setw(14) << c/g << "x";
    LOGGING_PRINTLN(line.str())
}

// B69: kernel-only exact 1..4 limb comparisons. Both forced algorithms
// overwrite precisely the same caller-owned output span; neither creates
// a BigUnsigned object, performs an allocation or changes production.
enum class B69Kernel { Fused, X86Carry, Selected };

template<B69Kernel Kernel>
std::array<BigUnsigned::Limb,8u> b69_preallocated(
        BigUnsigned const& first,BigUnsigned const& second)
{
    using Limb=BigUnsigned::Limb;
    if (first.size()>4u || second.size()>4u) {
        throw std::length_error("B69 diagnostic supports at most 4 limbs");
    }
    std::array<Limb,8u> output{};
    if (!first.is_zero() && !second.is_zero()) {
        const auto used=first.size()+second.size();
        const auto target=BigUnsigned::MutableLimbView(
            output.data(),used);
        if constexpr (Kernel==B69Kernel::Fused) {
            BigUnsigned::mul_into_fused(first.limbs(),second.limbs(),target);
        } else if constexpr (Kernel==B69Kernel::X86Carry) {
            BigUnsigned::mul_into_x86_64_carry(
                first.limbs(),second.limbs(),target);
        } else {
            // B70: same public kernel dispatch used by owning mul().
            BigUnsigned::mul_into(first.limbs(),second.limbs(),target);
        }
    }
    return output;
}

// Same signed Integer ownership on both forced paths; the only
// difference is the multiplication kernel. This deliberately differs
// from the production mul() allocation strategy, so forced/production
// ratios are whole-result controls, NOT isolated kernel timings.
template<B69Kernel Kernel>
Integer b69_owned(Integer const& first,Integer const& second)
{
    const auto& a=first.magnitude();
    const auto& b=second.magnitude();
    auto digits=b69_preallocated<Kernel>(a,b);
    const std::size_t used=a.size()+b.size();
    return Integer::from_magnitude(
        BigUnsigned(BigUnsigned::LimbView(digits.data(),used)),
        first.negative()!=second.negative());
}

void b69_check_signed(Integer const& actual,mpz_srcptr reference,
                      std::string const& label)
{
    RawInteger check;
    auto const words=actual.magnitude().limbs();
    mpz_import(check.value,words.size(),-1,sizeof(BigUnsigned::Limb),
               0,0,words.data());
    if (actual.negative()) mpz_neg(check.value,check.value);
    if (mpz_cmp(check.value,reference)!=0) {
        throw std::runtime_error("B69 signed GMP mismatch: "+label);
    }
}

void b69_check_preallocated(BigUnsigned const& a,BigUnsigned const& b,
                           mpz_srcptr oracle,std::string const& context)
{
    for (const B69Kernel kernel:
            {B69Kernel::Fused,B69Kernel::X86Carry,B69Kernel::Selected}) {
        const auto digits=kernel==B69Kernel::Fused ?
            b69_preallocated<B69Kernel::Fused>(a,b) :
            kernel==B69Kernel::X86Carry ?
            b69_preallocated<B69Kernel::X86Carry>(a,b) :
            b69_preallocated<B69Kernel::Selected>(a,b);
        RawInteger candidate;
        mpz_import(candidate.value,a.size()+b.size(),-1,
                   sizeof(BigUnsigned::Limb),0,0,digits.data());
        if (mpz_cmp(candidate.value,oracle)!=0) {
            throw std::runtime_error(
                "B69 preallocated GMP mismatch: "+context);
        }
    }
}

void verify_b69(unsigned bits,IntegerPool<IntegerOperations> const& native,
                IntegerPool<GmpIntegerOperations> const& gmp)
{
    if (bits>256u) return;
    RawInteger reference,abs_a,abs_b;
    for (unsigned i=0u;i<32u;++i) {
        const auto& a=native.a[i];
        const auto& b=native.b[i];
        const std::string context=std::to_string(bits)+" bits #"+
                                  std::to_string(i);
        mpz_mul(reference.value,gmp.a[i].value,gmp.b[i].value);
        b69_check_signed(b69_owned<B69Kernel::Fused>(a,b),
                         reference.value,context+" fused multiply");
        b69_check_signed(b69_owned<B69Kernel::X86Carry>(a,b),
                         reference.value,context+" carry multiply");
        b69_check_signed(b69_owned<B69Kernel::Selected>(a,b),
                         reference.value,context+" selected multiply");
        mpz_abs(abs_a.value,gmp.a[i].value);
        mpz_abs(abs_b.value,gmp.b[i].value);
        mpz_mul(reference.value,abs_a.value,abs_b.value);
        b69_check_preallocated(a.magnitude(),b.magnitude(),
                               reference.value,context+" prealloc multiply");

        mpz_mul(reference.value,gmp.a[i].value,gmp.a[i].value);
        b69_check_signed(b69_owned<B69Kernel::Fused>(a,a),
                         reference.value,context+" fused square");
        b69_check_signed(b69_owned<B69Kernel::X86Carry>(a,a),
                         reference.value,context+" carry square");
        b69_check_signed(b69_owned<B69Kernel::Selected>(a,a),
                         reference.value,context+" selected square");
        mpz_mul(reference.value,abs_a.value,abs_a.value);
        b69_check_preallocated(a.magnitude(),a.magnitude(),
                               reference.value,context+" prealloc square");
    }

    // Additional non-benchmark cases exercise canonical zero, each
    // limb length 1..4, saturated/dense/sparse inputs and both signs.
    std::mt19937_64 generator(0xB690ABCDULL);
    for (std::size_t n=0u;n<=4u;++n) {
        for (std::size_t m=0u;m<=4u;++m) {
            for (unsigned mode=0u;mode<12u;++mode) {
                std::array<std::uint64_t,4u> left{},right{};
                for (std::size_t j=0u;j<n;++j) left[j]=generator();
                for (std::size_t j=0u;j<m;++j) right[j]=generator();
                if (mode%4u==0u) {
                    std::fill(left.begin(),left.end(),~std::uint64_t{0});
                    std::fill(right.begin(),right.end(),~std::uint64_t{0});
                } else if (mode%4u==1u) {
                    std::fill(left.begin(),left.end(),0u);
                    std::fill(right.begin(),right.end(),0u);
                    if (n>0u) left[n-1u]=std::uint64_t{1};
                    if (m>0u) right[m-1u]=std::uint64_t{1};
                } else if (mode%4u==2u) {
                    for (std::size_t j=0u;j<n;j+=2u) left[j]=0u;
                    for (std::size_t j=0u;j<m;j+=2u) right[j]=0u;
                }
                const BigUnsigned a(BigUnsigned::LimbView(left.data(),n));
                const BigUnsigned b(BigUnsigned::LimbView(right.data(),m));
                mpz_import(abs_a.value,n,-1,sizeof(std::uint64_t),
                           0,0,left.data());
                mpz_import(abs_b.value,m,-1,sizeof(std::uint64_t),
                           0,0,right.data());
                mpz_mul(reference.value,abs_a.value,abs_b.value);
                const std::string context="short "+std::to_string(n)+
                    "x"+std::to_string(m)+" mode "+std::to_string(mode);
                b69_check_preallocated(a,b,reference.value,context);
                const Integer signed_a=Integer::from_magnitude(
                    a,(mode&1u)!=0u);
                const Integer signed_b=Integer::from_magnitude(
                    b,(mode&2u)!=0u);
                if (((mode&1u)!=0u)!=((mode&2u)!=0u)) {
                    mpz_neg(reference.value,reference.value);
                }
                b69_check_signed(
                    b69_owned<B69Kernel::Fused>(signed_a,signed_b),
                    reference.value,context+" fused");
                b69_check_signed(
                    b69_owned<B69Kernel::X86Carry>(signed_a,signed_b),
                    reference.value,context+" carry");
                b69_check_signed(
                    b69_owned<B69Kernel::Selected>(signed_a,signed_b),
                    reference.value,context+" selected");
            }
        }
    }
}

template<class Gmp,class Production,class Fused,class Carry>
void report_b69_owning(char const* label,
        Gmp&& gmp,Production&& production,Fused&& fused,Carry&& carry,
        unsigned iterations,unsigned samples)
{
    const double reference=exact_arithmetic_median_ns(
        std::forward<Gmp>(gmp),iterations,samples);
    const double selected=exact_arithmetic_median_ns(
        std::forward<Production>(production),iterations,samples);
    const double f=exact_arithmetic_median_ns(
        std::forward<Fused>(fused),iterations,samples);
    const double c=exact_arithmetic_median_ns(
        std::forward<Carry>(carry),iterations,samples);
    LOGGING_PRINTLN(std::left << std::setw(19) << label
                    << std::right << std::fixed << std::setprecision(3)
                    << " | GMP=" << reference << " ns"
                    << " | production=" << selected << " ns"
                    << " | forced_fused=" << f << " ns"
                    << " | forced_carry=" << c << " ns"
                    << " | carry/fused=" << c/f << "x"
                    << " | carry/GMP=" << c/reference << "x")
}

template<class Fused,class Carry>
void report_b69_kernel(char const* label,Fused&& fused,Carry&& carry,
                       unsigned iterations,unsigned samples)
{
    const double f=exact_arithmetic_median_ns(
        std::forward<Fused>(fused),iterations,samples);
    const double c=exact_arithmetic_median_ns(
        std::forward<Carry>(carry),iterations,samples);
    LOGGING_PRINTLN(std::left << std::setw(19) << label
                    << std::right << std::fixed << std::setprecision(3)
                    << " | fused=" << f << " ns"
                    << " | carry=" << c << " ns"
                    << " | carry/fused=" << c/f << "x")
}

// B70: GCC diagnostic hint prevents the compiler from inlining the
// exact same forced fused whole owning operation into the timing
// lambda. It does NOT reproduce every property of another TU.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
Integer b70_fused_noinline(Integer const& a,Integer const& b)
{
    return b69_owned<B69Kernel::Fused>(a,b);
}

// Compare a non-inlined non-aliased signed wrapper calling the real
// owning BigUnsigned multiplication to production Integer::mul.
// Both use the same signed result type and mathematical operation.
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
Integer b70_core_noinline(Integer const& a,Integer const& b)
{
    return Integer::from_magnitude(
        mul(a.magnitude(),b.magnitude()),
        a.negative()!=b.negative());
}

template<class Baseline,class Experimental>
void b70_pair(char const* operation,Baseline&& baseline,
              Experimental&& experimental,
              unsigned iterations,unsigned samples)
{
    const double n=exact_arithmetic_median_ns(
        std::forward<Baseline>(baseline),iterations,samples);
    const double candidate=exact_arithmetic_median_ns(
        std::forward<Experimental>(experimental),iterations,samples);
    LOGGING_PRINTLN(std::left<<std::setw(25)<<operation
                    <<std::right<<std::fixed<<std::setprecision(3)
                    <<" | baseline="<<n<<" ns"
                    <<" | candidate="<<candidate<<" ns"
                    <<" | candidate/baseline="<<candidate/n<<"x")
}

void report_b70(unsigned bits,
                IntegerPool<IntegerOperations> const& native,
                unsigned iterations,unsigned samples)
{
    if (bits>256u) return;
    if (!BigUnsigned::selected_multiply_uses_fused() ||
            !BigUnsigned::has_x86_64_carry_multiply()) return;
    const auto idx=[](unsigned i)->unsigned {return i%32u;};
    LOGGING_PRINTLN(" ")
    LOGGING_PRINTLN("B70 production-dispatch vs ownership/call-boundary "
                    "controls | "<<bits<<" bits")
    b70_pair("mul selected vs fused kernel",
        [&](unsigned i) {
            return b69_preallocated<B69Kernel::Fused>(
                native.a[idx(i)].magnitude(),native.b[idx(i)].magnitude());
        },[&](unsigned i) {
            return b69_preallocated<B69Kernel::Selected>(
                native.a[idx(i)].magnitude(),native.b[idx(i)].magnitude());
        },iterations,samples);
    b70_pair("square selected vs fused",
        [&](unsigned i) {
            auto const& a=native.a[idx(i)].magnitude();
            return b69_preallocated<B69Kernel::Fused>(a,a);
        },[&](unsigned i) {
            auto const& a=native.a[idx(i)].magnitude();
            return b69_preallocated<B69Kernel::Selected>(a,a);
        },iterations,samples);
    b70_pair("mul owned selected/fused",
        [&](unsigned i) {return b69_owned<B69Kernel::Fused>(
            native.a[idx(i)],native.b[idx(i)]);},
        [&](unsigned i) {return b69_owned<B69Kernel::Selected>(
            native.a[idx(i)],native.b[idx(i)]);},
        iterations,samples);
    b70_pair("square owned selected/fused",
        [&](unsigned i) {return b69_owned<B69Kernel::Fused>(
            native.a[idx(i)],native.a[idx(i)]);},
        [&](unsigned i) {return b69_owned<B69Kernel::Selected>(
            native.a[idx(i)],native.a[idx(i)]);},
        iterations,samples);
    b70_pair("mul forced inline/noinline",
        [&](unsigned i) {return b69_owned<B69Kernel::Fused>(
            native.a[idx(i)],native.b[idx(i)]);},
        [&](unsigned i) {return b70_fused_noinline(
            native.a[idx(i)],native.b[idx(i)]);},
        iterations,samples);
    b70_pair("square forced inline/noinline",
        [&](unsigned i) {return b69_owned<B69Kernel::Fused>(
            native.a[idx(i)],native.a[idx(i)]);},
        [&](unsigned i) {return b70_fused_noinline(
            native.a[idx(i)],native.a[idx(i)]);},
        iterations,samples);
    b70_pair("mul production/core noinline",
        [&](unsigned i) {return IntegerOperations::mul(
            native.a[idx(i)],native.b[idx(i)]);},
        [&](unsigned i) {return b70_core_noinline(
            native.a[idx(i)],native.b[idx(i)]);},
        iterations,samples);
    b70_pair("square production/core noinline",
        [&](unsigned i) {return IntegerOperations::square(
            native.a[idx(i)]);},
        [&](unsigned i) {return b70_core_noinline(
            native.a[idx(i)],native.a[idx(i)]);},
        iterations,samples);
    LOGGING_PRINTLN("B70 controls share GMP-verified B69 arithmetic "
                    "and the SAME 32 signed official inputs; "
                    "results are independently timed. "
                    "The noinline hint does not isolate full linking, "
                    "ABI or allocation costs.")
}

void report_b69(unsigned bits,
                IntegerPool<IntegerOperations> const& native,
                IntegerPool<GmpIntegerOperations> const& gmp,
                unsigned iterations,unsigned samples)
{
    if (bits>256u) return;
    if (!BigUnsigned::selected_multiply_uses_fused() ||
            !BigUnsigned::has_x86_64_carry_multiply()) {
        LOGGING_PRINTLN("B69 skipped: requires GCC x86-64 fused selection "
                        "and native x86 carry kernel.")
        return;
    }
    verify_b69(bits,native,gmp);
    const auto ix=[](unsigned i)->unsigned { return i%32u; };
    LOGGING_PRINTLN(" ")
    LOGGING_PRINTLN("B69 matched x86-64 1..4-limb fused/carry diagnostics | "
                    << bits << " bits | GMP oracle PASS")
    report_b69_kernel("mul_prealloc",
        [&](unsigned i) {
            const auto& a=native.a[ix(i)];
            const auto& b=native.b[ix(i)];
            return b69_preallocated<B69Kernel::Fused>(
                a.magnitude(),b.magnitude());
        },
        [&](unsigned i) {
            const auto& a=native.a[ix(i)];
            const auto& b=native.b[ix(i)];
            return b69_preallocated<B69Kernel::X86Carry>(
                a.magnitude(),b.magnitude());
        },iterations,samples);
    report_b69_kernel("square_prealloc",
        [&](unsigned i) {
            auto const& a=native.a[ix(i)].magnitude();
            return b69_preallocated<B69Kernel::Fused>(a,a);
        },
        [&](unsigned i) {
            auto const& a=native.a[ix(i)].magnitude();
            return b69_preallocated<B69Kernel::X86Carry>(a,a);
        },iterations,samples);

    report_b69_owning("mul_signed_owned",
        [&](unsigned i) { return GmpIntegerOperations::mul(
            gmp.a[ix(i)],gmp.b[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::mul(
            native.a[ix(i)],native.b[ix(i)]); },
        [&](unsigned i) { return b69_owned<B69Kernel::Fused>(
            native.a[ix(i)],native.b[ix(i)]); },
        [&](unsigned i) { return b69_owned<B69Kernel::X86Carry>(
            native.a[ix(i)],native.b[ix(i)]); },
        iterations,samples);
    report_b69_owning("square_signed_owned",
        [&](unsigned i) { return GmpIntegerOperations::square(
            gmp.a[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::square(
            native.a[ix(i)]); },
        [&](unsigned i) { return b69_owned<B69Kernel::Fused>(
            native.a[ix(i)],native.a[ix(i)]); },
        [&](unsigned i) { return b69_owned<B69Kernel::X86Carry>(
            native.a[ix(i)],native.a[ix(i)]); },
        iterations,samples);
    verify_b69(bits,native,gmp);
    report_b70(bits,native,iterations,samples);
    LOGGING_PRINTLN("B69/B70 GMP oracle PASS before/after timing; "
                    "prealloc rows are allocation-free kernel calls, "
                    "owning rows have matched signed Integer results. "
                    "No rows are additive; no change to production.")
}

void verify(unsigned bits,IntegerPool<IntegerOperations> const& native,
            IntegerPool<GmpIntegerOperations> const& gmp)
{
    RawInteger first,second,expected;
    for (unsigned i=0u;i<32u;++i) {
        const auto& a=native.a[i].magnitude();
        const auto& b=native.b[i].magnitude();
        mpz_abs(first.value,gmp.a[i].value);
        mpz_abs(second.value,gmp.b[i].value);
        const std::string label=std::to_string(bits)+" bits #"+
                                std::to_string(i);
        check_magnitude(a,first.value,label+" input a");
        check_magnitude(b,second.value,label+" input b");
        check_magnitude(from_decimal(abs_digits(native.spec[i].a)),
                        first.value,label+" decimal construction");
        if (to_decimal(a)!=decimal(first.value)) {
            throw std::runtime_error("B67 decimal magnitude mismatch: "+label);
        }

        mpz_mul(expected.value,first.value,second.value);
        check_magnitude(mul(a,b),expected.value,label+" mul");
        mpz_mul(expected.value,first.value,first.value);
        check_magnitude(square(a),expected.value,label+" square");
        mpz_pow_ui(expected.value,first.value,3u);
        check_magnitude(mul(square(a),a),expected.value,label+" cube");

        const auto& numerator=native.division_a[i].magnitude();
        const auto& divisor=native.division_b[i].magnitude();
        mpz_abs(first.value,gmp.division_a[i].value);
        mpz_abs(second.value,gmp.division_b[i].value);
        check_magnitude(numerator,first.value,label+" division a");
        check_magnitude(divisor,second.value,label+" division b");
        if (mpz_sgn(second.value)==0) {
            throw std::runtime_error("B67 unexpected zero divisor");
        }
        mpz_tdiv_q(expected.value,first.value,second.value);
        check_magnitude(quot(numerator,divisor),
                        expected.value,label+" quotient");
        mpz_tdiv_r(expected.value,first.value,second.value);
        check_magnitude(rem(numerator,divisor),
                        expected.value,label+" remainder");
    }
}

template<class GMP,class Native,class Core>
void report_row(char const* name,GMP&& gmp,Native&& native,
                Core&& core,unsigned iterations,unsigned samples)
{
    const double g=exact_arithmetic_median_ns(
        std::forward<GMP>(gmp),iterations,samples);
    const double n=exact_arithmetic_median_ns(
        std::forward<Native>(native),iterations,samples);
    const double c=exact_arithmetic_median_ns(
        std::forward<Core>(core),iterations,samples);
    std::ostringstream line;
    line << std::left << std::setw(21) << name
         << std::fixed << std::setprecision(3)
         << std::right << std::setw(11) << g
         << std::setw(11) << n
         << std::setw(11) << c
         << std::setw(11) << n/g << "x"
         << std::setw(11) << c/g << "x"
         << std::setw(11) << n/c << "x";
    LOGGING_PRINTLN(line.str())
}

void measure(unsigned bits,unsigned iterations,unsigned samples)
{
    IntegerPool<IntegerOperations> native(bits);
    IntegerPool<GmpIntegerOperations> gmp(bits);
    verify(bits,native,gmp);
    const auto ix=[](unsigned i)->unsigned { return i%32u; };
    LOGGING_PRINTLN(" ")
    LOGGING_PRINTLN(bits<<" bits | 32 GMP-verified OFFICIAL integer/division input pairs"
                   <<" | iterations="<<iterations<<" | samples="<<samples)
    LOGGING_PRINTLN(std::left<<std::setw(21)<<"operation"
                   <<std::right<<std::setw(11)<<"GMP ns"
                   <<std::setw(11)<<"Integer ns"
                   <<std::setw(11)<<"abs core ns"
                   <<std::setw(12)<<"Integer/GMP"
                   <<std::setw(12)<<"core/GMP"
                   <<std::setw(12)<<"Integer/core")

    report_row("construct_decimal",
        [&](unsigned i) {
            return GmpIntegerOperations::parse(gmp.spec[ix(i)].a);
        },[&](unsigned i) {
            return IntegerOperations::parse(native.spec[ix(i)].a);
        },[&](unsigned i) {
            return from_decimal(abs_digits(native.spec[ix(i)].a));
        },iterations,samples);

    report_row("copy",
        [&](unsigned i) { return GmpIntegerOperations::copy(gmp.a[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::copy(native.a[ix(i)]); },
        [&](unsigned i) { return BigUnsigned(native.a[ix(i)].magnitude()); },
        iterations,samples);

    report_row("compare",
        [&](unsigned i) { return GmpIntegerOperations::compare(gmp.a[ix(i)],gmp.b[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::compare(native.a[ix(i)],native.b[ix(i)]); },
        [&](unsigned i) { return compare(native.a[ix(i)].magnitude(),native.b[ix(i)].magnitude()); },
        iterations,samples);

    report_row("mul",
        [&](unsigned i) { return GmpIntegerOperations::mul(gmp.a[ix(i)],gmp.b[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::mul(native.a[ix(i)],native.b[ix(i)]); },
        [&](unsigned i) { return mul(native.a[ix(i)].magnitude(),native.b[ix(i)].magnitude()); },
        iterations,samples);

    report_row("square",
        [&](unsigned i) { return GmpIntegerOperations::square(gmp.a[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::square(native.a[ix(i)]); },
        [&](unsigned i) { return square(native.a[ix(i)].magnitude()); },
        iterations,samples);

    report_row("quotient",
        [&](unsigned i) { return GmpIntegerOperations::quotient(gmp.division_a[ix(i)],gmp.division_b[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::quotient(native.division_a[ix(i)],native.division_b[ix(i)]); },
        [&](unsigned i) { return quot(native.division_a[ix(i)].magnitude(),native.division_b[ix(i)].magnitude()); },
        iterations,samples);

    report_row("remainder",
        [&](unsigned i) { return GmpIntegerOperations::remainder(gmp.division_a[ix(i)],gmp.division_b[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::remainder(native.division_a[ix(i)],native.division_b[ix(i)]); },
        [&](unsigned i) { return rem(native.division_a[ix(i)].magnitude(),native.division_b[ix(i)].magnitude()); },
        iterations,samples);

    report_row("pow3",
        [&](unsigned i) { return GmpIntegerOperations::power3(gmp.a[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::power3(native.a[ix(i)]); },
        [&](unsigned i) {
            auto const& a=native.a[ix(i)].magnitude();
            return mul(square(a),a);
        },iterations,samples);

    report_row("to_decimal",
        [&](unsigned i) { return GmpIntegerOperations::to_decimal(gmp.a[ix(i)]); },
        [&](unsigned i) { return IntegerOperations::to_decimal(native.a[ix(i)]); },
        [&](unsigned i) { return to_decimal(native.a[ix(i)].magnitude()); },
        iterations,samples);

    verify(bits,native,gmp);
    verify_signed_diagnostics(bits,native,gmp);
    LOGGING_PRINTLN(" ")
    LOGGING_PRINTLN("B68 SAME SIGNED SEMANTICS: independent owning A/B"
                    " comparisons | "+std::to_string(bits)+" bits")
    LOGGING_PRINTLN(std::left << std::setw(22) << "operation"
                    << std::right << std::setw(12) << "GMP ns"
                    << std::setw(14) << "production ns"
                    << std::setw(14) << "candidate ns"
                    << std::setw(15) << "candidate/prod"
                    << std::setw(16) << "candidate/GMP")
    report_signed_candidate("construct_signed",
        [&](unsigned i) {
            return GmpIntegerOperations::parse(gmp.spec[ix(i)].a);
        },[&](unsigned i) {
            return IntegerOperations::parse(native.spec[ix(i)].a);
        },[&](unsigned i) {
            return parse_direct_signed_diagnostic(native.spec[ix(i)].a);
        },iterations,samples);
    report_signed_candidate("to_decimal_signed",
        [&](unsigned i) {
            return GmpIntegerOperations::to_decimal(gmp.a[ix(i)]);
        },[&](unsigned i) {
            return IntegerOperations::to_decimal(native.a[ix(i)]);
        },[&](unsigned i) {
            return format_direct_signed_diagnostic(native.a[ix(i)]);
        },iterations,samples);
    report_signed_candidate("mul_signed",
        [&](unsigned i) {
            return GmpIntegerOperations::mul(gmp.a[ix(i)],gmp.b[ix(i)]);
        },[&](unsigned i) {
            return IntegerOperations::mul(native.a[ix(i)],native.b[ix(i)]);
        },[&](unsigned i) {
            return multiply_direct_signed_diagnostic(
                native.a[ix(i)],native.b[ix(i)]);
        },iterations,samples);
    report_signed_candidate("square_signed",
        [&](unsigned i) {
            return GmpIntegerOperations::square(gmp.a[ix(i)]);
        },[&](unsigned i) {
            return IntegerOperations::square(native.a[ix(i)]);
        },[&](unsigned i) {
            return square_direct_signed_diagnostic(native.a[ix(i)]);
        },iterations,samples);
    verify_signed_diagnostics(bits,native,gmp);
    LOGGING_PRINTLN("B68: every candidate preserves the signed GMP result. "
                    "Rows are timed independently, not additive. "
                    "A/B timings can still contain host drift.")
    report_b69(bits,native,gmp,iterations,samples);
    LOGGING_PRINTLN("abs core = unsigned production BigUnsigned operations "
                    "on the SAME magnitudes. Signed overhead cannot "
                    "be inferred as an exact difference because "
                    "signed controls can take different branches.")
    LOGGING_PRINTLN("All rows timed independently; do not sum components.")
}

} // namespace

int main(int argc,const char* argv[])
{
    try {
        BenchmarkCommandLine::initialize_logging();
        ExactArithmeticSelection selection;
        selection.sizes={53u,128u,256u,1024u,2048u,4096u};
        CommandLineInterface command_line;
        command_line.add(logging_command_line());
        command_line.add(CommandLineContribution(
            "numeric-integer-small-operations-profile",
            exact_arithmetic_parsers(selection,false)));
        const auto arguments=BenchmarkCommandLine::normalize_legacy_arguments(
            argc,argv,{"iterations","samples"},
            {"--iterations","-i","--samples","-n","--bits","-b",
             "--scheduler","-s","--theme","-t","--verbosity","-v"});
        if (!command_line.acquire(arguments)) {
            return BenchmarkCommandLine::requests_help(arguments) ? 0 : 1;
        }
        for (const unsigned bits:selection.sizes) {
            measure(bits,selection.iterations,selection.samples);
        }
        LOGGING_PRINTLN(" ")
        LOGGING_PRINTLN("Numerical verification PASS: all magnitude core "
                        "operations, B68 signed candidates, B69 forced kernels "
                        "and B70 matched selected/forced controls "
                        "versus independent GMP.")
        return 0;
    } catch (std::exception const& error) {
        std::cerr<<"B67 owning-vs-core diagnostic FAILED: "
                 <<error.what()<<'\n';
        return 1;
    }
}
