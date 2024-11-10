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
  if constexpr (put_get_helper_detail::IndexableBy<PropertyMap, std::ranges::range_reference_t<KeyRange>>) {
    return std::forward<KeyRange>(k_range) |
           std::views::transform([pmap](auto&& k) -> decltype(auto) { return pmap[std::forward<decltype(k)>(k)]; });
  } else {
    return std::forward<KeyRange>(k_range) | std::views::transform([pmap](auto&& k) -> decltype(auto) {
             return get(pmap, std::forward<decltype(k)>(k));
           });
  }
}

template <typename Graph>
auto make_single_graph_range(Graph& g) {
  return std::ranges::single_view{&g} | std::views::transform([](Graph* pg) -> Graph& { return *pg; });
}

template <typename Graph, typename PropertyTag>
auto make_graph_property_map_range(Graph& g, PropertyTag tag) {
  if constexpr (std::is_same_v<property_kind_t<PropertyTag>, edge_property_tag>) {
    return make_property_map_range(get(tag, g), edges(g));
  } else if constexpr (std::is_same_v<property_kind_t<PropertyTag>, vertex_property_tag>) {
    return make_property_map_range(get(tag, g), vertices(g));
  } else {
    return make_property_map_range(get(tag, g), make_single_graph_range(g));
  }
}

}  // namespace bagl

#endif  // BAGL_BAGL_PROPERTY_MAP_RANGE_H_
