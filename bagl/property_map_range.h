// (C) Copyright Jeremy Siek, 2001.
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PROPERTY_MAP_RANGE_H_
#define BAGL_BAGL_PROPERTY_MAP_RANGE_H_

#include <ranges>
#include <type_traits>

#include "bagl/property_map.h"
#include "bagl/property.h"

namespace bagl {

//======================================================================
// property range, generalized from ideas by Francois Faure

template <typename PropertyMap, std::ranges::input_range KeyRange>
auto make_property_map_range(PropertyMap pmap, KeyRange&& k_range) {
  if constexpr (concepts::LvaluePropertyMap<PropertyMap, std::ranges::range_value_t<KeyRange>>) {
    return std::forward<KeyRange>(k_range) | std::views::transform([pmap](const auto& k) { return pmap[k]; });
  } else {
    return std::forward<KeyRange>(k_range) | std::views::transform([pmap](const auto& k) { return get(pmap, k); });
  }
}

template <typename Graph, typename PropertyTag>
auto make_graph_property_map_range(Graph& graph, PropertyTag tag) {
  if constexpr (std::is_same_v<property_kind_t<PropertyTag>, edge_property_tag>) {
    return make_property_map_range(get(tag, graph), edges(graph));
  } else if constexpr (std::is_same_v<property_kind_t<PropertyTag>, vertex_property_tag>) {
    return make_property_map_range(get(tag, graph), vertices(graph));
  } else {
    return make_property_map_range(get(tag, graph), std::ranges::single_view{std::ref(graph)});
  }
}

}  // namespace bagl

#endif  // BAGL_BAGL_PROPERTY_MAP_RANGE_H_
