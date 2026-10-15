// Part of Ariadne Numeric. GPL-3.0-or-later.
#ifndef ARIADNE_FLOATAP_HPP
#define ARIADNE_FLOATAP_HPP

#include "numeric/bits.hpp"
#include "numeric/float.decl.hpp"
#include <cstddef>
#include <ostream>

namespace Ariadne {

// The precision tag is deliberately independent of MPFR types.
class ArbitraryPrecision {
    std::size_t _bits;
public:
    explicit ArbitraryPrecision(Bits count) : _bits(static_cast<unsigned long>(count)) { }
    explicit ArbitraryPrecision(std::size_t count) : _bits(count) { }
    std::size_t bits() const noexcept { return _bits; }
    friend bool operator==(ArbitraryPrecision a, ArbitraryPrecision b) noexcept {
        return a._bits==b._bits;
    }
};

// Float<ArbitraryPrecision> is not yet specialized. BUILD_FLOATAP_COMPARISON
// remains disabled until src/floatap.cpp and the specialization are provided.
using FloatAP = Float<ArbitraryPrecision>;
}
#endif
