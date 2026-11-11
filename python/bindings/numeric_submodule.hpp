/***************************************************************************
 *            python/bindings/numeric_submodule.hpp
 *
 *  Copyright  2021  Pieter Collins
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

/** @file python/bindings/numeric_submodule.hpp
 *  @brief Declarations and utilities for numeric submodule in Ariadne Numeric.
 */

#ifndef ARIADNE_PYTHON_NUMERIC_SUBMODULE_HPP
#define ARIADNE_PYTHON_NUMERIC_SUBMODULE_HPP

#include <pybind11/operators.h>

#include "pybind11.hpp"

#include "numeric-utilities.hpp"

namespace Ariadne {

/** @brief Return the Python class tag for the selected numeric precision. */
template<class T>
std::string numeric_class_tag();
/** @brief Return the Python class tag for the selected numeric precision. */
template<>
inline std::string numeric_class_tag<DoublePrecision>() { return "DP"; }
/** @brief Return the Python class tag for the selected numeric precision. */
template<>
inline std::string numeric_class_tag<MultiplePrecision>() { return "MP"; }

/** @brief Compile-time Python template-name customization. */
template<>
struct PythonTemplateName<Positive> {
    /** @brief Return the requested element or stored value. */
    static std::string get() { return "Positive"; }
};

/** @brief Compile-time Python template-name customization. */
template<>
struct PythonTemplateName<Float> {};
/** @brief Compile-time Python template-name customization. */
template<>
struct PythonTemplateName<Approximation> {
    /** @brief Return the requested element or stored value. */
    static std::string get() { return "Approximation"; }
};
/** @brief Compile-time Python template-name customization. */
template<>
struct PythonTemplateName<LowerBound> {
    /** @brief Return the requested element or stored value. */
    static std::string get() { return "LowerBound"; }
};
/** @brief Compile-time Python template-name customization. */
template<>
struct PythonTemplateName<UpperBound> {
    /** @brief Return the requested element or stored value. */
    static std::string get() { return "UpperBound"; }
};
/** @brief Compile-time Python template-name customization. */
template<>
struct PythonTemplateName<Bounds> {
    /** @brief Return the requested element or stored value. */
    static std::string get() { return "Bounds"; }
};
/** @brief Compile-time Python template-name customization. */
template<>
struct PythonTemplateName<Ball> {
    /** @brief Return the requested element or stored value. */
    static std::string get() { return "Ball"; }
};
/** @brief Compile-time Python template-name customization. */
template<>
struct PythonTemplateName<Error> {
    /** @brief Return the requested element or stored value. */
    static std::string get() { return "Error"; }
};
/** @brief Compile-time Python template-name customization. */
template<>
struct PythonTemplateName<Rounded> {
    /** @brief Return the requested element or stored value. */
    static std::string get() { return "Rounded"; }
};

/** @brief Compile-time Python template-name customization. */
template<>
struct PythonTemplateName<Number> {
    /** @brief Return the requested element or stored value. */
    static std::string get() { return "Number"; }
};

/** @brief Compile-time Python class-name customization. */
template<>
struct PythonClassName<FloatDP> {
    /** @brief Return the requested element or stored value. */
    static std::string get() { return "FloatDP"; }
};
/** @brief Compile-time Python class-name customization. */
template<>
struct PythonClassName<FloatMP> {
    /** @brief Return the requested element or stored value. */
    static std::string get() { return "FloatMP"; }
};

/** @brief Compile-time Python class-name customization. */
template<template<class> class T, class X>
struct PythonClassName<T<X>> {
    /** @brief Return the requested element or stored value. */
    static std::string get() { return python_class_name<X>() + python_template_name<T>(); }
};

/** @brief Compile-time Python class-name customization. */
template<class F, class FE>
struct PythonClassName<Ball<F, FE>> {
    /** @brief Return the requested element or stored value. */
    static std::string get() {
        typedef typename FE::PrecisionType PRE;
        return python_class_name<F>() + numeric_class_tag<PRE>() + "Ball";
    }
};
/** @brief Compile-time Python class-name customization. */
template<class F>
struct PythonClassName<Ball<F>> {
    /** @brief Return the requested element or stored value. */
    static std::string get() { return python_class_name<F>() + "Ball"; }
};

/** @brief Write the value to the output stream. */
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Dyadic>& repr);
/** @brief Write the value to the output stream. */
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Decimal>& repr);
/** @brief Write the value to the output stream. */
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Rational>& repr);
/** @brief Write the value to the output stream. */
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Real>& repr);
/** @brief Write the value to the output stream. */
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatDP>& repr);
/** @brief Write the value to the output stream. */
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatMP>& repr);

/** @brief Write the value to the output stream. */
template<class P>
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Number<P>>& repr);

/** @brief Write the value to the output stream. */
template<class F>
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Approximation<F>>& repr);
/** @brief Write the value to the output stream. */
template<class F>
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<LowerBound<F>>& repr);
/** @brief Write the value to the output stream. */
template<class F>
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<UpperBound<F>>& repr);
/** @brief Write the value to the output stream. */
template<class F>
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Bounds<F>>& repr);
/** @brief Write the value to the output stream. */
template<class F, class FE>
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Ball<F, FE>>& repr);
/** @brief Write the value to the output stream. */
template<class F>
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<F>& repr);
/** @brief Write the value to the output stream. */
template<class FE>
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Error<FE>>& repr);

/** @brief Write the value to the output stream. */
OutputStream& operator<<(OutputStream& os, const PythonLiteral<Dyadic>& lit);
/** @brief Write the value to the output stream. */
OutputStream& operator<<(OutputStream& os, const PythonLiteral<Decimal>& lit);
/** @brief Write the value to the output stream. */
OutputStream& operator<<(OutputStream& os, const PythonLiteral<Rational>& lit);
/** @brief Write the value to the output stream. */
OutputStream& operator<<(OutputStream& os, const PythonLiteral<Real>& lit);

/** @brief Write the value to the output stream. */
OutputStream& operator<<(OutputStream& os, const PythonLiteral<FloatDP>& lit);
/** @brief Write the value to the output stream. */
OutputStream& operator<<(OutputStream& os, const PythonLiteral<FloatMP>& lit);
/** @brief Write the value to the output stream. */
template<class F>
OutputStream& operator<<(OutputStream& os, const PythonLiteral<Approximation<F>>& lit);
/** @brief Write the value to the output stream. */
template<class F>
OutputStream& operator<<(OutputStream& os, const PythonLiteral<LowerBound<F>>& lit);
/** @brief Write the value to the output stream. */
template<class F>
OutputStream& operator<<(OutputStream& os, const PythonLiteral<UpperBound<F>>& lit);
/** @brief Write the value to the output stream. */
template<class F>
OutputStream& operator<<(OutputStream& os, const PythonLiteral<Bounds<F>>& lit);
/** @brief Write the value to the output stream. */
template<class F, class FE>
OutputStream& operator<<(OutputStream& os, const PythonLiteral<Ball<F, FE>>& lit);
/** @brief Write the value to the output stream. */
template<class F>
OutputStream& operator<<(OutputStream& os, const PythonLiteral<F>& lit);
/** @brief Write the value to the output stream. */
template<class FE>
OutputStream& operator<<(OutputStream& os, const PythonLiteral<Error<FE>>& lit);

/** @brief Convert a Python value or literal to the requested numeric type. */
template<class X>
X from_python_object_or_literal(pybind11::handle h) {
    if constexpr (Constructible<X, std::string>) {
        try {
            return X(pybind11::cast<std::string>(h));
        } catch (...) {
        }
    }
    return pybind11::cast<X>(h);
}
/** @brief Convert a Python value or literal to the requested numeric type. */
template<>
ValidatedNumber from_python_object_or_literal<ValidatedNumber>(pybind11::handle h);
/** @brief Convert a Python value or literal to the requested numeric type. */
template<>
ValidatedUpperNumber from_python_object_or_literal<ValidatedUpperNumber>(pybind11::handle h);
/** @brief Convert a Python value or literal to the requested numeric type. */
template<>
ValidatedLowerNumber from_python_object_or_literal<ValidatedLowerNumber>(pybind11::handle h);
/** @brief Convert a Python value or literal to the requested numeric type. */
template<>
ApproximateNumber from_python_object_or_literal<ApproximateNumber>(pybind11::handle h);

} // namespace Ariadne

#endif /* ARIADNE_PYTHON_NUMERIC_SUBMODULE_HPP */
