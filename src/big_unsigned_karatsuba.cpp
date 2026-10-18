/***************************************************************************
 *            numeric/big_unsigned_karatsuba.cpp
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

/*! \file numeric/big_unsigned_karatsuba.cpp
 *  \brief Experimental allocation-free Karatsuba full-product kernel.
 */

#include "numeric/big_unsigned.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>

namespace Ariadne {

namespace {

using Limb = BigUnsigned::Limb;
using Read = BigUnsigned::LimbView;
using Write = BigUnsigned::MutableLimbView;

// This is a recursion cutoff, not a proven performance crossover.
// Small cases use the current architecture-selected schoolbook kernel.
constexpr std::size_t karatsuba_base_limbs=8u;

std::size_t checked_sum(std::size_t a,std::size_t b)
{
    if (a>std::numeric_limits<std::size_t>::max()-b) {
        throw std::length_error("BigUnsigned Karatsuba size overflow");
    }
    return a+b;
}

bool should_split(std::size_t a,std::size_t b) noexcept
{
    const std::size_t half=std::max(a,b)/2u;
    return std::min(a,b)>karatsuba_base_limbs
        && a>half && b>half;
}

// Exact recursion-frame bound: sumA, sumB, cross product, plus the
// maximum of the three sequential recursive calls' workspaces.
// No nested frame overlaps the live parent sums or cross product.
std::size_t karatsuba_workspace(std::size_t a,std::size_t b)
{
    if (!should_split(a,b)) {
        return 0u;
    }

    const std::size_t half=std::max(a,b)/2u;
    const std::size_t a0=std::min(a,half);
    const std::size_t b0=std::min(b,half);
    const std::size_t a1=a-a0;
    const std::size_t b1=b-b0;
    const std::size_t sa=checked_sum(std::max(a0,a1),1u);
    const std::size_t sb=checked_sum(std::max(b0,b1),1u);
    const std::size_t pair=checked_sum(sa,sb);
    const std::size_t frame=checked_sum(pair,pair);

    const std::size_t nested=std::max({
        karatsuba_workspace(a0,b0),
        karatsuba_workspace(a1,b1),
        karatsuba_workspace(sa,sb)
    });
    return checked_sum(frame,nested);
}

// An extra high zero limb makes addition exact even when both halves
// are saturated. Arguments must not overlap the destination.
void sum_halves(Read a,Read b,Write output) noexcept
{
    Limb carry=0u;
    for (std::size_t i=0u;i<output.size();++i) {
        const Limb left=i<a.size() ? a[i] : Limb{0};
        const Limb right=i<b.size() ? b[i] : Limb{0};
        const Limb sum=left+right;
        const Limb value=sum+carry;
        carry=static_cast<Limb>(sum<left || value<sum);
        output[i]=value;
    }
}

// The caller has proved lhs >= rhs as unsigned integers. Subtracting
// from lhs itself avoids an extra allocated intermediate buffer.
void subtract_in_place(Write lhs,Read rhs)
{
    Limb borrow=0u;
    for (std::size_t i=0u;i<lhs.size();++i) {
        const Limb right=i<rhs.size() ? rhs[i] : Limb{0};
        const Limb left=lhs[i];
        const Limb difference=left-right;
        const Limb value=difference-borrow;
        borrow=static_cast<Limb>(left<right || difference<borrow);
        lhs[i]=value;
    }
    if (borrow!=0u) {
        throw std::logic_error("BigUnsigned Karatsuba negative cross term");
    }
}

// Accumulate a shifted cross term into the exact full product output.
// A widened sum-of-halves temporary can have more high zero limbs
// than fit in the final output; those discarded words must be zero.
void add_shifted(Write output,Read term)
{
    Limb carry=0u;
    for (std::size_t i=0u;i<output.size();++i) {
        const Limb left=output[i];
        const Limb right=i<term.size() ? term[i] : Limb{0};
        const Limb sum=left+right;
        const Limb value=sum+carry;
        carry=static_cast<Limb>(sum<left || value<sum);
        output[i]=value;
    }
    if (carry!=0u) {
        throw std::logic_error("BigUnsigned Karatsuba output overflow");
    }
    for (std::size_t i=output.size();i<term.size();++i) {
        if (term[i]!=0u) {
            throw std::logic_error("BigUnsigned Karatsuba nonzero high term");
        }
    }
}

// output.size()==a.size()+b.size(); scratch size is precomputed.
// z0 and z2 occupy disjoint output slices, so no extra result buffer
// is required for either of them. Only sums and cross term use scratch.
void karatsuba_recursive(Read a,Read b,Write output,Write scratch)
{
    if (!should_split(a.size(),b.size())) {
        BigUnsigned::mul_low_into(a,b,output);
        return;
    }

    std::fill(output.begin(),output.end(),Limb{0});
    const std::size_t half=std::max(a.size(),b.size())/2u;
    const Read a0=a.first(std::min(a.size(),half));
    const Read b0=b.first(std::min(b.size(),half));
    const Read a1=a.subspan(a0.size());
    const Read b1=b.subspan(b0.size());

    karatsuba_recursive(
        a0,b0,output.first(a0.size()+b0.size()),scratch);
    karatsuba_recursive(
        a1,b1,output.subspan(2u*half,a1.size()+b1.size()),scratch);

    const std::size_t sa=std::max(a0.size(),a1.size())+1u;
    const std::size_t sb=std::max(b0.size(),b1.size())+1u;
    const std::size_t pair=sa+sb;

    Write sum_a=scratch.first(sa);
    Write sum_b=scratch.subspan(sa,sb);
    Write cross=scratch.subspan(pair,pair);
    Write nested=scratch.subspan(2u*pair);

    sum_halves(a0,a1,sum_a);
    sum_halves(b0,b1,sum_b);
    karatsuba_recursive(sum_a,sum_b,cross,nested);

    // (a0+a1)(b0+b1)-a0*b0-a1*b1 = a0*b1+a1*b0.
    subtract_in_place(cross,output.first(a0.size()+b0.size()));
    subtract_in_place(
        cross,output.subspan(2u*half,a1.size()+b1.size()));
    add_shifted(output.subspan(half),cross);
}

} // namespace

std::size_t BigUnsigned::karatsuba_scratch_limbs(
        std::size_t a_limbs,std::size_t b_limbs)
{
    (void)checked_sum(a_limbs,b_limbs);
    return karatsuba_workspace(a_limbs,b_limbs);
}

void BigUnsigned::mul_into_karatsuba(
        LimbView a,LimbView b,
        MutableLimbView output,MutableLimbView scratch)
{
    if (a.empty() || b.empty()) {
        std::fill(output.begin(),output.end(),Limb{0});
        return;
    }

    const std::size_t count=checked_sum(a.size(),b.size());
    if (output.size()<count) {
        throw std::length_error("BigUnsigned Karatsuba output too small");
    }
    const std::size_t required=karatsuba_workspace(a.size(),b.size());
    if (scratch.size()<required) {
        throw std::length_error("BigUnsigned Karatsuba scratch too small");
    }

    std::fill(output.begin(),output.end(),Limb{0});
    karatsuba_recursive(a,b,output.first(count),scratch.first(required));
}

} // namespace Ariadne
