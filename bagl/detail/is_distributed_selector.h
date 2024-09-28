// Copyright 2012 The Trustees of Indiana University.
// Authors: Jeremiah Willcock, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

// Selector to determine whether a selector is distributed_s (can only be true
// if "bagl/distributed/selector.h" has been included) so that we can
// disable various sequential-graph-only traits specializations for distributed
// graphs.

#ifndef BAGL_BAGL_DETAIL_IS_DISTRIBUTED_SELECTOR_H_
#define BAGL_BAGL_DETAIL_IS_DISTRIBUTED_SELECTOR_H_

#include <type_traits>

namespace bagl::distributed_detail {
template <typename>
struct is_distributed_selector : std::false_type {};

template <typename T>
constexpr bool is_distributed_selector_v = is_distributed_selector<T>::value;
}  // namespace bagl::distributed_detail

#endif  // BAGL_BAGL_DETAIL_IS_DISTRIBUTED_SELECTOR_H_
