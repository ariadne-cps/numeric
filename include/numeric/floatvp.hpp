// Part of Ariadne Numeric. GPL-3.0-or-later.
#ifndef ARIADNE_FLOATVP_HPP
#define ARIADNE_FLOATVP_HPP

#include "numeric/bits.hpp"
#include "numeric/float.decl.hpp"
#include <cstddef>
#include <ostream>

namespace Ariadne {

// The precision tag is deliberately independent of MPFR types.
class VariablePrecision {
    std::size_t _bits;
public:
    explicit VariablePrecision(Bits count) : _bits(static_cast<unsigned long>(count)) { }
    explicit VariablePrecision(std::size_t count) : _bits(count) { }
    std::size_t bits() const noexcept { return _bits; }
    friend bool operator==(VariablePrecision a, VariablePrecision b) noexcept {
        return a._bits==b._bits;
    }
};

// Float<VariablePrecision> is not yet specialized. BUILD_FLOATVP_COMPARISON
// remains disabled until src/floatvp.cpp and the specialization are provided.
using FloatVP = Float<VariablePrecision>;
}
#endif
