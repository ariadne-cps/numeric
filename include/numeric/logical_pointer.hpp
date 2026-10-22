#ifndef ARIADNE_NUMERIC_LOGICAL_POINTER_HPP
#define ARIADNE_NUMERIC_LOGICAL_POINTER_HPP

#include <memory>

#include "foundation/logical.hpp"

namespace Ariadne::Detail {

// The comparison dispatch API returns an owning pointer. LogicalHandle
// retains its object, while a value-based logical result must release it.
template<class P> LogicalType<P> logical_type_from_owned_pointer(LogicalInterface* ptr) {
    std::unique_ptr<LogicalInterface> owner(ptr);
    if constexpr (ConstructibleFrom<LogicalType<P>,LogicalHandle>) {
        return LogicalType<P>(LogicalHandle(std::shared_ptr<LogicalInterface>(std::move(owner))));
    } else {
        return LogicalType<P>(logical_value_from_pointer(owner.get()));
    }
}

} // namespace Ariadne::Detail

#endif
