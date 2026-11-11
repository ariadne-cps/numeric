/***************************************************************************
 *            logical_sequence.cpp
 *
 *  Copyright  2026  Ariadne contributors
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Numeric. If not, see <https://www.gnu.org/licenses/>.
 */

#include "numeric/logical_sequence.hpp"

#include "numeric/integer.hpp"
#include "numeric/sequence.hpp"
namespace Ariadne {
namespace {
class SequenceDisjunction : public LogicalInterface {
    Sequence<LowerKleenean> _sequence;

  public:
    explicit SequenceDisjunction(Sequence<LowerKleenean> s) : _sequence(s) {}

  private:
    LogicalInterface* _copy() const override { return new SequenceDisjunction(*this); }
    LogicalValue _check(Effort e) const override {
        for (Natural k = 0u; k != e.work(); ++k)
            if (definitely(_sequence[k].check(e)))
                return LogicalValue::TRUE;
        return LogicalValue::INDETERMINATE;
    }
    OutputStream& _write(OutputStream& os) const override { return os << "disjunction(" << _sequence[0u] << "," << _sequence[1u] << "," << _sequence[2u] << ",...)"; }
};
class SequenceConjunction : public LogicalInterface {
    Sequence<UpperKleenean> _sequence;

  public:
    explicit SequenceConjunction(Sequence<UpperKleenean> s) : _sequence(s) {}

  private:
    LogicalInterface* _copy() const override { return new SequenceConjunction(*this); }
    LogicalValue _check(Effort e) const override {
        for (Natural k = 0u; k != e.work(); ++k)
            if (definitely(not _sequence[k].check(e)))
                return LogicalValue::FALSE;
        return LogicalValue::INDETERMINATE;
    }
    OutputStream& _write(OutputStream& os) const override { return os << "conjunction(" << _sequence[0u] << "," << _sequence[1u] << "," << _sequence[2u] << ",...)"; }
};
} // namespace
LowerKleenean disjunction(const Sequence<LowerKleenean>& s) {
    return LowerKleenean(LogicalHandle(make_handle<const SequenceDisjunction>(s)));
}
UpperKleenean conjunction(const Sequence<UpperKleenean>& s) {
    return UpperKleenean(LogicalHandle(make_handle<const SequenceConjunction>(s)));
}
} // namespace Ariadne
