/***************************************************************************
 *            numeric/operators.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

/** @file numeric/operators.hpp
 *  @brief Numerical operator classes
 */

#ifndef ARIADNE_OPERATORS_HPP
#define ARIADNE_OPERATORS_HPP

#include <cassert>
#include <iostream>

#include "foundation/logical.decl.hpp"
#include "number.decl.hpp"
#include "utility/string.hpp"
#include "utility/variant.hpp"
#include "utility/variant.inl.hpp"

namespace Ariadne {

/** @brief Logical result associated with numeric comparisons. */
template<class X>
struct Logic;
/** @brief Logical result associated with numeric comparisons. */
template<>
struct Logic<String> {
    typedef Boolean Type;
};
/** @brief Logical result associated with numeric comparisons. */
template<>
struct Logic<Integer> {
    typedef Boolean Type;
};
/** @brief Logical result associated with numeric comparisons. */
template<>
struct Logic<Real> {
    typedef Kleenean Type;
};
template<class X>
using LogicType = typename Logic<X>::Type;

template<class X1, class X2>
using RingType = DifferenceType<ProductType<X1, X2>, ProductType<X2, X1>>;
template<class X1, class X2>
using FieldType = DifferenceType<QuotientType<X1, X2>, QuotientType<X2, X1>>;
template<class X>
using TranscendentalType = decltype(sin(declval<X>()));

/** @brief Callable operation object for . */
class Operator {
  public:
    enum class Code : ComparableEnumerationType;
    enum class Kind : ComparableEnumerationType;
    /** @brief Return the category discriminator. */
    friend Kind get_kind(Code);

  private:
    Code _code;
    Kind _kind;

  public:
    /** @brief Construct a \c Operator from c, k. */
    Operator(Code c, Kind k) : _code(c), _kind(k) {}
    /** @brief Construct a \c Operator from c. */
    Operator(Code c) : _code(c), _kind(get_kind(_code)) {}
    /** @brief Construct a \c Operator from op. */
    template<class OP>
    Operator(const OP& op) : _code(op.code()), _kind(op.kind()) {}

    /** @brief Convert to code. */
    operator Code() const { return _code; }
    /** @brief Return the operation or variant discriminator. */
    Code code() const { return _code; }
    /** @brief Return the category discriminator. */
    Kind kind() const { return _kind; }
};

using OperatorCode = Operator::Code;
using OperatorKind = Operator::Kind;

enum class Operator::Kind : ComparableEnumerationType {
    VARIABLE,
    COORDINATE,
    NULLARY,
    UNARY,
    BINARY,
    TERNARY,
    GRADED,
    SCALAR,
    PREDICATE,
    COMPARISON
};

enum class Operator::Code : ComparableEnumerationType {
    CNST, // A constant value
    VAR,  // A named variable
    IND,  // A numbered index
    ADD,  // Addition
    SUB,  // Subtraction
    MUL,  // Multiplication
    DIV,  // Division
    FMA,  // Fused multiply-and-add
    POW,  // Integer power
    ROOT, // Integer root
    NUL,  // Zero
    POS,  // Unary plus
    NEG,  // Unary negation
    HLF,  // Halve
    REC,  // Reciprocal
    SQR,  // Square
    SQRT, // Square root
    EXP,  // Natural exponential
    LOG,  // Natural logarithm
    SIN,  // Sine
    COS,  // Cosine
    TAN,  // Tangent
    ASIN, // ArcSine
    ACOS, // ArcCosine
    ATAN, // ArcTangent
    SADD, // Addition of a scalar constant
    SSUB, // Scalar Subtraction
    SMUL, // Multiplication by a scalar constant
    SDIV, // Scalar Division
    SFMA, // Fused multiply-and-add with scalar constants
    ABS,  // Absolute value
    MAX,  // Maximum
    MIN,  // Minimum
    NOT,  // Logical not
    AND,  // Logical and
    OR,   // Logical or
    XOR,  // Logical exclusive or
    IMPL, // Logical implication
    ITOR, // Conversion of Int to Real
    GET,  // Subscript operator
    VEC,  // Vectorisation operator
    PUSH,
    PULL,
    TANH,      // Hyperbolic tangent
    EQ = -6,   // Equal
    NEQ = -5,  // Not equal
    GEQ = -4,  // Greater or equal
    LEQ = -3,  // Less or equal
    GT = -2,   // Greater than
    LT = -1,   // Less than
    SGN = -7,  // Compare with zero
    SUBS = -8, // Compare as a subset
    DISJ = -9  // Compare disjointness
};

/** @brief Write the value to the output stream. */
OutputStream& operator<<(OutputStream& os, const OperatorKind& knd);
/** @brief Write the value to the output stream. */
OutputStream& operator<<(OutputStream&, const OperatorCode& op);
/** @brief Return the name. */
const char* name(const OperatorCode& op);
/** @brief Return the display symbol. */
const char* symbol(const OperatorCode& op);

/** @brief Stored numeric operation object. */
template<class OBJ>
struct Object {
    /** @brief Convert to the common interface or value type. */
    OBJ const& upcast() const { return static_cast<OBJ const&>(*this); }
};
/** @brief Object storing an operation and its operands. */
template<class OP>
struct OperatorObject : Object<OP> {};
/** @brief Object storing a comparison operation and its operands. */
template<class CMP>
struct ComparisonObject : Object<CMP> {};

/** @brief Write the value to the output stream. */
template<class OP>
inline OutputStream& operator<<(OutputStream& os, const OperatorObject<OP>& op) {
    return os << op.upcast().code();
}

/** @brief Write the value to the output stream. */
template<class CMP>
inline OutputStream& operator<<(OutputStream& os, const ComparisonObject<CMP>& cmp) {
    return os << cmp.upcast().code();
}

/** @brief Positive-sign comparison operation. */
struct GtrZero {};
/** @brief Negative-sign comparison operation. */
struct LessZero {};

/** @brief Greater-than comparison operation. */
struct Gtr : ComparisonObject<Gtr> {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 > a2) { return a1 > a2; }
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::GT; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::COMPARISON; }
};

/** @brief Ordering predicate using  order. */
struct Less : ComparisonObject<Less> {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 < a2) { return a1 < a2; }
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::LT; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::COMPARISON; }
};

/** @brief Greater-than-or-equal comparison operation. */
struct Geq : ComparisonObject<Geq> {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 >= a2) { return a1 >= a2; }
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::GEQ; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::COMPARISON; }
};

/** @brief Less-than-or-equal comparison operation. */
struct Leq : ComparisonObject<Leq> {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 <= a2) { return a1 <= a2; }
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::LEQ; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::COMPARISON; }
};

/** @brief Equality comparison operation. */
struct Equal : ComparisonObject<Equal> {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 == a2) { return a1 == a2; }
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::EQ; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::COMPARISON; }
};

/** @brief Inequality comparison operation. */
struct Unequal : ComparisonObject<Unequal> {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 != a2) { return a1 != a2; }
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::NEQ; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::COMPARISON; }
};

using Eq = Equal;
using Ne = Unequal;
using Lt = Less;
using Gt = Gtr;
using Le = Leq;
using Ge = Geq;

/** @brief Logical exclusive-disjunction operation. */
struct XOrOp : OperatorObject<XOrOp> {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 xor a2) { return a1 xor a2; }
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::XOR; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
};
/** @brief Logical conjunction operation. */
struct AndOp : OperatorObject<AndOp> {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 and a2) { return a1 and a2; }
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::AND; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
};
/** @brief Logical disjunction operation. */
struct OrOp : OperatorObject<OrOp> {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 or a2) { return a1 or a2; }
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::OR; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
};
/** @brief Logical negation operation. */
struct NotOp : OperatorObject<NotOp> {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(not a) { return not a; }
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::NOT; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
};

/** @brief Constant-expression operation token. */
struct Cnst : OperatorObject<Cnst> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::CNST; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::NULLARY; }
};

/** @brief Index-expression operation token. */
struct Ind : OperatorObject<Ind> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::IND; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::COORDINATE; }
};

/** @brief Variable-expression operation token. */
struct Var : OperatorObject<Var> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::VAR; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::VARIABLE; }
};

/** @brief Addition operation token. */
struct Plus {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(+a) { return +a; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 + a2) { return a1 + a2; }
};
/** @brief Subtraction operation token. */
struct Minus {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(-a) { return -a; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 - a2) { return a1 - a2; }
};
/** @brief Multiplication operation token. */
struct Times {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 * a2) { return a1 * a2; }
};
/** @brief Divisibility predicate. */
struct Divides {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(a1 / a2) { return a1 / a2; }
};

/** @brief Addition operation with explicit rounding. */
struct RAdd : OperatorObject<RAdd> {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(add(a2, a1)) { return add(a2, a1); }
};
/** @brief Subtraction operation with explicit rounding. */
struct RSub : OperatorObject<RSub> {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(sub(a2, a1)) { return sub(a2, a1); }
};
/** @brief Multiplication operation with explicit rounding. */
struct RMul : OperatorObject<RMul> {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(mul(a2, a1)) { return mul(a2, a1); }
};
/** @brief Division operation with explicit rounding. */
struct RDiv : OperatorObject<RDiv> {
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(div(a2, a1)) { return div(a2, a1); }
};

/** @brief Tag selecting first argument. */
struct FirstArgumentTag {};
/** @brief Tag selecting second argument. */
struct SecondArgumentTag {};
/** @brief Tag selecting third argument. */
struct ThirdArgumentTag {};

/** @brief Binary addition operation. */
struct Add : OperatorObject<Add> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::ADD; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(add(a1, a2)) { return add(a1, a2); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a1, const X&, FirstArgumentTag) const { return a1.create(1); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X&, const X& a2, SecondArgumentTag) const { return a2.create(1); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X&, const D& d1, const X&, const D& d2) const { return d1 + d2; }
};
/** @brief Binary subtraction operation. */
struct Sub : OperatorObject<Sub> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::SUB; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(sub(a1, a2)) { return sub(a1, a2); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a1, const X&, FirstArgumentTag) const { return a1.create(+1); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X&, const X& a2, SecondArgumentTag) const { return a2.create(-1); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X&, const D& d1, const X&, const D& d2) const { return d1 - d2; }
};
/** @brief Binary multiplication operation. */
struct Mul : OperatorObject<Mul> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::MUL; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(mul(a1, a2)) { return mul(a1, a2); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X&, const X& a2, FirstArgumentTag) const { return a2; }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a1, const X&, SecondArgumentTag) const { return a1; }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const { return a2 * d1 + a1 * d2; }
};
/** @brief Binary division operation. */
struct Div : OperatorObject<Div> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::DIV; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(div(a1, a2)) { return div(a1, a2); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X&, const X& a2, FirstArgumentTag) const { return rec(a2); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a1, const X& a2, SecondArgumentTag) const { return div(neg(a1), sqr(a2)); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const { return (d1 - d2 * (a1 / a2)) / a2; }
};

/** @brief Integer-graded power operation. */
struct Pow : OperatorObject<Pow> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::POW; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::GRADED; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A, class N>
    auto operator()(A&& a, N&& n) const -> decltype(pow(a, n)) { return pow(a, n); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class N>
    X derivative(const X& a, N n) const { return n * pow(a, n - 1); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D, class N>
    D derivative(const X& a, const D& d, N n) const { return derivative(a, n) * d; }
};
/** @brief Integer-graded root operation. */
struct Root : OperatorObject<Root> {
    /** @brief Return the operation or variant discriminator. */
    OperatorCode code() const { return OperatorCode::ROOT; }
    /** @brief Return the category discriminator. */
    OperatorKind kind() const { return OperatorKind::GRADED; }
    /** @brief Return the root of the specified order. */
    template<class A, class N>
    friend auto root(A&& a, N&& n) -> decltype(exp(log(a) / n)) { return exp(log(a) / n); }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A, class N>
    auto operator()(A&& a, N&& n) const -> decltype(root(a, n)) { return root(a, n); }
};
/** @brief Ternary multiply-add operation. */
struct Fma : OperatorObject<Fma> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::FMA; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::TERNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2, class A3>
    auto operator()(A1&& a1, A2&& a2, A3&& a3) const -> decltype(fma(a1, a2, a3)) { return fma(a1, a2, a3); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X&, const X& a2, const X&, FirstArgumentTag) const { return a2; }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a1, const X&, const X&, SecondArgumentTag) const { return a1; }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X&, const X&, const X& a3, ThirdArgumentTag) const { return nul(a3) + 1; }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a1, const D& d1, const X& a2, const D& d2, const X&, const D& d3) const {
        return a2 * d1 + a1 * d2 + d3;
    }
};

/** @brief Unary operation creating an additive zero. */
struct Nul : OperatorObject<Nul> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::NUL; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(nul(a)) { return nul(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return nul(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X&, const D& d) const { return nul(d); }
};
/** @brief Unary positive operation. */
struct Pos : OperatorObject<Pos> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::POS; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(pos(a)) { return pos(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return nul(a) + 1; }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X&, const D& d) const { return d; }
};
/** @brief Unary negation operation. */
struct Neg : OperatorObject<Neg> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::NEG; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(neg(a)) { return neg(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return nul(a) - 1; }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X&, const D& d) const { return -d; }
};
/** @brief Unary division-by-two operation. */
struct Hlf : OperatorObject<Hlf> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::HLF; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(hlf(a)) { return hlf(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return hlf(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X&, const D& d) const { return hlf(d); }
};
/** @brief Unary reciprocal operation. */
struct Rec : OperatorObject<Rec> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::REC; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(rec(a)) { return rec(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return rec(neg(sqr(a))); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a, const D& d) const { return d / (neg(sqr(a))); }
};
/** @brief Unary squaring operation. */
struct Sqr : OperatorObject<Sqr> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::SQR; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(sqr(a)) { return sqr(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return 2 * a; }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a, const D& d) const { return (2 * a) * d; }
};
/** @brief Unary square-root operation. */
struct Sqrt : OperatorObject<Sqrt> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::SQRT; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(sqrt(a)) { return sqrt(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return rec(2 * sqrt(a)); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a, const D& d) const { return d / (2 * sqrt(a)); }
};
/** @brief Unary exponential operation. */
struct Exp : OperatorObject<Exp> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::EXP; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(exp(a)) { return exp(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return exp(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a, const D& d) const { return exp(a) * d; }
};
/** @brief Unary natural-logarithm operation. */
struct Log : OperatorObject<Log> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::LOG; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(log(a)) { return log(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return rec(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a, const D& d) const { return d / a; }
};
/** @brief Unary sine operation. */
struct Sin : OperatorObject<Sin> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::SIN; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(sin(a)) { return sin(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return cos(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a, const D& d) const { return derivative(a) * d; }
};
/** @brief Unary cosine operation. */
struct Cos : OperatorObject<Cos> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::COS; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(cos(a)) { return cos(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return neg(sin(a)); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a, const D& d) const { return derivative(a) * d; }
};
/** @brief Unary tangent operation. */
struct Tan : OperatorObject<Tan> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::TAN; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(tan(a)) { return tan(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return rec(sqr(cos(a))); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a, const D& d) const { return derivative(a) * d; }
};
/** @brief Unary hyperbolic-tangent operation. */
struct Tanh : OperatorObject<Tanh> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::TANH; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(tanh(a)) { return tanh(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return 1 - sqr(tanh(a)); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a, const D& d) const { return derivative(a) * d; }
};
/** @brief Unary inverse-sine operation. */
struct Asin : OperatorObject<Asin> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::ASIN; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(asin(a)) { return asin(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return rec(1 - sqrt(a)); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a, const D& d) const { return derivative(a) * d; }
};
/** @brief Unary inverse-cosine operation. */
struct Acos : OperatorObject<Acos> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::ACOS; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(acos(a)) { return acos(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return rec(1 - sqrt(a)); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a, const D& d) const { return derivative(a) * d; }
};
/** @brief Unary inverse-tangent operation. */
struct Atan : OperatorObject<Atan> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::ATAN; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(atan(a)) { return atan(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& a) const { return rec(1 + sqr(a)); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a, const D& d) const { return derivative(a) * d; }
};

/** @brief Binary maximum operation. */
struct Max : OperatorObject<Max> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::MAX; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(max(a1, a2)) { return max(a1, a2); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const { return a1 >= a2 ? d1 : d2; }
};

/** @brief Binary minimum operation. */
struct Min : OperatorObject<Min> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::MIN; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::BINARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A1, class A2>
    auto operator()(A1&& a1, A2&& a2) const -> decltype(min(a1, a2)) { return min(a1, a2); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a1, const D& d1, const X& a2, const D& d2) const { return a1 <= a2 ? d1 : d2; }
};

/** @brief Unary absolute-value operation. */
struct Abs : OperatorObject<Abs> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::ABS; }
    /** @brief Return the category discriminator. */
    static constexpr OperatorKind kind() { return OperatorKind::UNARY; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(abs(a)) { return abs(a); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X, class D>
    D derivative(const X& a, const D&) const { return a >= 0 ? a : -a; }
};

/** @brief Unary sign operation. */
struct Sgn : ComparisonObject<Sgn> {
    /** @brief Return the operation or variant discriminator. */
    static constexpr OperatorCode code() { return OperatorCode::SGN; }
    /** @brief Return the category discriminator. */
    static OperatorKind kind() { return OperatorKind::PREDICATE; }
    /** @brief Apply the callable to the supplied arguments. */
    template<class A>
    auto operator()(A&& a) const -> decltype(sgn(a)) { return sgn(a); }
    /** @brief Apply the callable to the supplied arguments. */
    Kleenean operator()(const Real& a) const;
};

/*
template<class OP> struct CodedOperator {
    typedef OP CodeType;
    operator Operator() const { return Operator(static_cast<Operator::Code>(_code)); }
    CodeType code() const { return _code; }
  protected:
    CodeType _code;
};
struct SpecialOperator : CodedOperator<SpecialOperator::Code> {
    enum class Code : char;
};
*/

/** @brief Tagged variant of the supported operation types. */
template<class... OPS>
class OperatorVariant
    : public CodedVariant<OperatorCode, OPS...> {
  public:
    /** @brief Construct a \c OperatorVariant from op. */
    template<class OP>
        requires OneOf<OP, OPS...>
    OperatorVariant(OP op) : CodedVariant<OperatorCode, OPS...>(op) {}
    /** @brief Construct a \c OperatorVariant from code. */
    explicit OperatorVariant(OperatorCode code) : CodedVariant<OperatorCode, OPS...>(code) {}
    /** @brief Return the category discriminator. */
    OperatorKind kind() const {
        return this->accept([](auto op) { return op.kind(); });
    }
    /** @brief Apply the callable to the supplied arguments. */
    template<class... AS>
    decltype(auto) operator()(AS&&... as) const {
        return this->accept([&as...](auto op) { return op(std::forward<AS>(as)...); });
    }
    /** @brief Evaluate with the supplied arguments. */
    template<class... AS>
    decltype(auto) call(AS&&... as) const {
        return this->accept([&as...](auto op) { return op(std::forward<AS>(as)...); });
    }
    /** @brief Evaluate using the requested result representation. */
    template<class R, class... AS>
    R call_as(AS&&... as) const {
        return this->accept([&as...](auto op) { return static_cast<R>(op(std::forward<AS>(as)...)); });
    }
    /** @brief Evaluate at the supplied argument. */
    template<class... AS>
    friend decltype(auto) evaluate(const OperatorVariant<OPS...>& ops, AS&&... as) {
        return ops.call(as...);
    }
    /** @brief Evaluate using the requested result representation. */
    template<class R, class... AS>
    friend R evaluate_as(const OperatorVariant<OPS...>& ops, AS&&... as) {
        return ops.call_as(as...);
    }
    /** @brief Write the value to the output stream. */
    friend OutputStream& operator<<(OutputStream& os, const OperatorVariant<OPS...>& op) {
        op.accept([&os](auto op_) { os << op_; });
        return os;
    }
};

/** @brief Callable operation object for unary logical. */
struct UnaryLogicalOperator : OperatorVariant<NotOp> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
};
/** @brief Callable operation object for unary comparison. */
struct UnaryComparisonOperator : OperatorVariant<Sgn> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
};
/** @brief Callable operation object for unary ring. */
struct UnaryRingOperator : OperatorVariant<Neg> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
};
/** @brief Callable operation object for unary arithmetic. */
struct UnaryArithmeticOperator : OperatorVariant<Nul, Pos, Neg, Sqr, Rec> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
};
/** @brief Callable operation object for unary transcendental. */
struct UnaryTranscendentalOperator : OperatorVariant<Pos, Neg, Sqr, Hlf, Rec, Sqrt, Exp, Log, Sin, Cos, Tan, Tanh, Atan> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
};
/** @brief Callable operation object for unary elementary. */
struct UnaryElementaryOperator : OperatorVariant<Nul, Pos, Neg, Sqr, Hlf, Rec, Sqrt, Exp, Log, Sin, Cos, Tan, Tanh, Asin, Acos, Atan, Abs> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    decltype(sin(declval<X>())) operator()(X&& x) const {
        typedef decltype(sin(x)) R;
        return this->accept([&x](auto op) { return static_cast<R>(op(std::forward<X>(x))); });
    }
};
/** @brief Callable operation object for unary lattice. */
struct UnaryLatticeOperator : OperatorVariant<Abs> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
};
/** @brief Callable operation object for binary logical. */
struct BinaryLogicalOperator : OperatorVariant<AndOp, OrOp> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
};
/** @brief Callable operation object for binary comparison. */
struct BinaryComparisonOperator : OperatorVariant<Equal, Unequal, Less, Gtr, Leq, Geq> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
};
/** @brief Callable operation object for binary lattice. */
struct BinaryLatticeOperator : OperatorVariant<Max, Min> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
};
/** @brief Callable operation object for binary ring. */
struct BinaryRingOperator : OperatorVariant<Add, Sub, Mul> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
    /** @brief Apply the callable to the supplied arguments. */
    template<class X1, class X2>
    ArithmeticType<X1, X2> operator()(X1&& x1, X2&& x2) const {
        return this->accept([&x1, &x2](auto op) { return static_cast<ArithmeticType<X1, X2>>(op(std::forward<X1>(x1), std::forward<X2>(x2))); });
    }
};
/** @brief Callable operation object for binary field. */
struct BinaryFieldOperator : OperatorVariant<Add, Sub, Mul, Div> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
    /** @brief Apply the callable to the supplied arguments. */
    template<class X1, class X2>
    QuotientType<X1, X2> operator()(X1&& x1, X2&& x2) const {
        return this->accept([&x1, &x2](auto op) { return static_cast<QuotientType<X1, X2>>(op(std::forward<X1>(x1), std::forward<X2>(x2))); });
    }
};
using BinaryArithmeticOperator = BinaryFieldOperator;
/** @brief Callable operation object for binary elementary. */
struct BinaryElementaryOperator : OperatorVariant<Add, Sub, Mul, Div, Max, Min> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
    /** @brief Apply the callable to the supplied arguments. */
    template<class X1, class X2>
    QuotientType<X1, X2> operator()(X1&& x1, X2&& x2) const {
        return this->accept([&x1, &x2](auto op) { return static_cast<QuotientType<X1, X2>>(op(std::forward<X1>(x1), std::forward<X2>(x2))); });
    }
};
/** @brief Callable operation object for graded ring. */
struct GradedRingOperator : OperatorVariant<Pow> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
    /** @brief Apply the callable to the supplied arguments. */
    template<class X, class N>
    decltype(pow(declval<X>(), declval<N>())) operator()(X&& x, const N& n) const {
        return this->accept([&x, &n](auto op) { return static_cast<decltype(pow(declval<X>(), declval<N>()))>(op(std::forward<X>(x), n)); });
    }
};
/** @brief Callable operation object for graded field. */
struct GradedFieldOperator : OperatorVariant<Pow> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
    /** @brief Apply the callable to the supplied arguments. */
    template<class X, class N>
    decltype(pow(declval<X>(), declval<N>())) operator()(X&& x, const N& n) const {
        return this->accept([&x, &n](auto op) { return static_cast<decltype(pow(declval<X>(), declval<N>()))>(op(std::forward<X>(x), n)); });
    }
};
/** @brief Callable operation object for graded elementary. */
struct GradedElementaryOperator : OperatorVariant<Pow> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
    /** @brief Apply the callable to the supplied arguments. */
    template<class X, class N>
    decltype(pow(declval<X>(), declval<N>())) operator()(X&& x, const N& n) const {
        return this->accept([&x, &n](auto op) { return static_cast<decltype(pow(declval<X>(), declval<N>()))>(op(std::forward<X>(x), n)); });
    }
};
/** @brief Callable operation object for ternary arithmetic. */
struct TernaryArithmeticOperator : OperatorVariant<Fma> {
    /** @brief Construct a \c OperatorVariant from code. */
    using OperatorVariant::OperatorVariant;
};

/** @brief Return the inverse. */
constexpr Pos inverse(Pos) { return Pos(); }
/** @brief Return the inverse. */
constexpr Neg inverse(Neg) { return Neg(); }
/** @brief Return the inverse. */
constexpr Rec inverse(Rec) { return Rec(); }
/** @brief Return the inverse. */
constexpr Sqrt inverse(Sqr) { return Sqrt(); }
/** @brief Return the inverse. */
constexpr Sqr inverse(Sqrt) { return Sqr(); }
/** @brief Return the inverse. */
constexpr Log inverse(Exp) { return Log(); }
/** @brief Return the inverse. */
constexpr Exp inverse(Log) { return Exp(); }
/** @brief Return the inverse. */
constexpr Atan inverse(Tan) { return Atan(); }
/** @brief Return the inverse. */
constexpr Tan inverse(Atan) { return Tan(); }

template<class OP>
using InverseType = decltype(inverse(declval<OP>()));

template<class OP>
concept HasInverse = requires(OP op) { inverse(op); };

namespace {
/** @brief Test whether the two operations are inverses. */
template<class OP, class... OPS>
Bool _are_inverses(OP const& op1, const OperatorVariant<OPS...>& ops2) {
    if constexpr (HasInverse<decltype(op1)>) {
        return inverse(op1).code() == ops2.code();
    } else {
        return false;
    }
}
} // namespace

/** @brief Test whether the two operations are inverses. */
template<class... OPS>
Bool are_inverses(const OperatorVariant<OPS...>& ops1, const OperatorVariant<OPS...>& ops2) {
    return ops1.accept([&ops2](auto op1) { return _are_inverses(op1, ops2); });
}

/** @brief Callable operation object for unary. */
class UnaryOperator {
    OperatorCode _op;

  public:
    /** @brief Construct a \c UnaryOperator from op. */
    template<class OP>
    UnaryOperator(OP op) : _op(op) {}
    /** @brief Apply the callable to the supplied arguments. */
    template<class X>
    X operator()(const X& x) { return compute(_op, x); }
    /** @brief Return the derivative with respect to the selected argument. */
    template<class X>
    X derivative(const X& x) { return derivative(_op, x); }
};

/** @brief Write the value to the output stream. */
inline OutputStream& operator<<(OutputStream& os, const Less&) { return os << "<"; }
/** @brief Write the value to the output stream. */
inline OutputStream& operator<<(OutputStream& os, const Gtr&) { return os << ">"; }
/** @brief Write the value to the output stream. */
inline OutputStream& operator<<(OutputStream& os, const Leq&) { return os << "<="; }
/** @brief Write the value to the output stream. */
inline OutputStream& operator<<(OutputStream& os, const Geq&) { return os << ">="; }
/** @brief Write the value to the output stream. */
inline OutputStream& operator<<(OutputStream& os, const Equal&) { return os << "=="; }
/** @brief Write the value to the output stream. */
inline OutputStream& operator<<(OutputStream& os, const Unequal&) { return os << "!="; }

/** @brief Write the value to the output stream. */
inline OutputStream& operator<<(OutputStream& os, const AndOp&) { return os << "&&"; }
/** @brief Write the value to the output stream. */
inline OutputStream& operator<<(OutputStream& os, const OrOp&) { return os << "||"; }
/** @brief Write the value to the output stream. */
inline OutputStream& operator<<(OutputStream& os, const NotOp&) { return os << "!"; }

/** @brief Write the value to the output stream. */
inline OutputStream& operator<<(OutputStream& os, const Plus&) { return os << "+"; }
/** @brief Write the value to the output stream. */
inline OutputStream& operator<<(OutputStream& os, const Minus&) { return os << "-"; }
/** @brief Write the value to the output stream. */
inline OutputStream& operator<<(OutputStream& os, const Times&) { return os << "*"; }
/** @brief Write the value to the output stream. */
inline OutputStream& operator<<(OutputStream& os, const Divides&) { return os << "/"; }

} // namespace Ariadne

#endif // ARIADNE_OPERATORS_HPP
