/***************************************************************************
 *            numeric/sequence.hpp
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

/** @file numeric/sequence.hpp
 *  @brief Declarations and utilities for sequence in Ariadne Numeric.
 */

#ifndef ARIADNE_SEQUENCE_HPP
#define ARIADNE_SEQUENCE_HPP

#include <functional>

#include "utility/typedefs.hpp"
#include "utility/writable.hpp"

namespace Ariadne {

class Natural;

class Dyadic;
class Rational;
class Real;
template<class Y>
struct CompletionTypedef;
template<>
struct CompletionTypedef<Dyadic> {
    typedef Real Type;
};
template<>
struct CompletionTypedef<Rational> {
    typedef Real Type;
};
template<>
struct CompletionTypedef<Real> {
    typedef Real Type;
};
template<class Y>
using CompletionType = typename CompletionTypedef<Y>::Type;

//! \brief A function \f$\mathbb{N} \to X\f$.
template<class X>
class Sequence {
    std::function<X(Natural)> _fn;

  public:
    Sequence(const std::function<X(Natural)>& fn) : _fn(fn) {}
    X operator[](const Natural& n) const { return _fn(n); }
};

//! \brief A convergent sequence in \f$X\f$, with no further information about the convergence rate
template<class X>
class ConvergentSequence : public Sequence<X> {
  public:
    ConvergentSequence(const std::function<X(Natural)>& fn) : Sequence<X>(fn) {}
    ConvergentSequence(const Sequence<X>& seq) : Sequence<X>(seq) {}
};

//! An alternating sequence in a partially-ordered space (X,,≤), satisfying \f$x_n \in [x_{n-1},x_{n-2}]\f$.
template<class X>
class AlternatingSequence : public ConvergentSequence<X> {
  public:
    explicit AlternatingSequence(const Sequence<X>& seq) : ConvergentSequence<X>(seq) {}
};

//! \brief An increasing sequence in a partially-ordered space (X,≤).
template<class X>
class IncreasingSequence : public Sequence<X> {
  public:
    IncreasingSequence(std::function<X(Natural)> fn) : Sequence<X>(fn) {}
    IncreasingSequence(const Sequence<X>& seq) : Sequence<X>(seq) {}
};
//! \brief A decreasing sequence in a partially-ordered space (X,≤).
template<class X>
class DecreasingSequence : public Sequence<X> {
  public:
    DecreasingSequence(std::function<X(Natural)> fn) : Sequence<X>(fn) {}
    DecreasingSequence(const Sequence<X>& seq) : Sequence<X>(seq) {}
};

//! \brief A fast-converging Cauchy sequence in a metric space (X,d), satisfying \f$ d(x_m,x_n) \leq 2^{-min(m,n)}\f$.
template<class X>
class FastCauchySequence : public Sequence<X> {
  public:
    FastCauchySequence(const std::function<X(Natural)>& fn) : Sequence<X>(fn) {}
    FastCauchySequence(const Sequence<X>& seq) : Sequence<X>(seq) {}
    friend CompletionType<X> limit(const FastCauchySequence<X>&);
};

class SequenceWriter {
    Nat _num;

  public:
    SequenceWriter(const Nat& n) : _num(n) {}
    template<class T>
    inline WritableTemporary<Sequence<T>, SequenceWriter> operator()(const Sequence<T>& seq) const {
        return WritableTemporary<Sequence<T>, SequenceWriter>(*this, seq);
    }
    template<class T>
    OutputStream& _write(OutputStream& os, const Sequence<T>& seq) const {
        os << "[";
        for (Nat i = 0u; i != _num; ++i) {
            if (i != 0) {
                os << ",";
            }
            os << seq[i];
        }
        os << "]";
        return os;
    }
};

} // namespace Ariadne

#endif
