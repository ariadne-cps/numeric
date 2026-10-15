// Part of Ariadne Numeric. GPL-3.0-or-later.
#ifndef ARIADNE_FLOATNP_HPP
#define ARIADNE_FLOATNP_HPP

#include "numeric/bits.hpp"
#include "numeric/float.decl.hpp"
#include <cstddef>
#include <ostream>

namespace Ariadne {

// The precision tag is deliberately independent of MPFR types.
class NativePrecision {
    std::size_t _bits;
public:
    explicit NativePrecision(Bits count) : _bits(static_cast<unsigned long>(count)) { }
    explicit NativePrecision(std::size_t count) : _bits(count) { }
    std::size_t bits() const noexcept { return _bits; }
    friend bool operator==(NativePrecision a, NativePrecision b) noexcept {
        return a._bits==b._bits;
    }
};

// Float<NativePrecision> is not yet specialized. BUILD_FLOATNP_COMPARISON
// remains disabled until src/floatnp.cpp and the specialization are provided.
using FloatNP = Float<NativePrecision>;
}
#endif
