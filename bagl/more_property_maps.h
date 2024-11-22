// Copyright 2013-2024 Mikael Persson.

#ifndef BAGL_BAGL_MORE_PROPERTY_MAPS_H_
#define BAGL_BAGL_MORE_PROPERTY_MAPS_H_

#include <utility>

#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/property_map.h"

namespace bagl {

//======== Whole-bundle property-map ==========

// This property-map delivers the entire vertex, edge or graph bundle associated to a
// vertex / edge descriptor of a graph, or the graph itself.
// This property-map relies on the operator[] from the graph.
template <typename Graph, typename BundleTag>
struct whole_bundle_property_map : public put_get_helper<whole_bundle_property_map<Graph, BundleTag>> {
 public:
  static constexpr bool is_vertex_bundle_v = is_vertex_property_kind_v<BundleTag>;
  static constexpr bool is_edge_bundle_v = is_edge_property_kind_v<BundleTag>;
  using non_const_graph = std::remove_cv_t<Graph>;
  using value_type = std::conditional_t<
      is_vertex_bundle_v, vertex_bundle_type<non_const_graph>,
      std::conditional_t<is_edge_bundle_v, edge_bundle_type<non_const_graph>, graph_bundle_type<non_const_graph>>>;

  explicit whole_bundle_property_map(Graph* pg) : pg_(pg) {}
  whole_bundle_property_map() = default;

  template <typename Key>
  decltype(auto) operator[](Key&& k) const {
    if constexpr (is_vertex_bundle_v || is_edge_bundle_v) {
      return (*pg_)[std::forward<Key>(k)];
    } else {
      return std::forward<Key>(k)[graph_bundle];
    }
  }

 private:
  Graph* pg_ = nullptr;
};

//======== Tagged-in-property property-map ==========

// This property-map uses a graph's "get_property" function to get the
// whole property and then extract the value associated to a given tag.
// Unlikely to be useful for direct use, but is generally what is formed
// by calls like get(vertex_distance, g).
template <typename T, typename Graph, typename PropertyMapTag>
struct tagged_in_property_property_map
    : public put_get_helper<tagged_in_property_property_map<T, Graph, PropertyMapTag>> {
 private:
  Graph* pg_ = nullptr;

 public:
  using value_type = T;

  explicit tagged_in_property_property_map(Graph* pg, PropertyMapTag /*tag*/ = {}) : pg_(pg) {}
  tagged_in_property_property_map() = default;

  template <typename Key>
  decltype(auto) operator[](Key&& k) const {
    if constexpr (std::is_same_v<std::remove_cv_t<value_type>, no_property>) {
      return no_property{};
    } else if constexpr (is_vertex_property_kind_v<PropertyMapTag> || is_edge_property_kind_v<PropertyMapTag>) {
      return get_property_value(get_property(*pg_, std::forward<Key>(k)), PropertyMapTag{});
    } else {
      return get_property_value(get_property(std::forward<Key>(k), graph_all), PropertyMapTag{});
    }
  }
};

//======== Bundle-data-member property-map ==========

// This property-map delivers a data-member of the vertex or edge bundle associated to a
// vertex / edge descriptor of a graph. This is similar to the property-map obtained by
// calling get(&SomeBundle::SomeMember, my_graph), and can be used to implement such a
// functionality for custom graph classes.
// This property-map relies of the existence of an operator[] for the graph and
// for the relevant descriptor type. This operator[] should deliver the vertex/edge bundle.
template <typename T, typename Graph, typename BundleTag>
class bundle_member_property_map : public put_get_helper<bundle_member_property_map<T, Graph, BundleTag>> {
 public:
  using bundle_map = whole_bundle_property_map<Graph, BundleTag>;
  using bundle_type = typename bundle_map::value_type;
  using member_ptr_type = T bundle_type::*;

 private:
  bundle_map bmap_;
  member_ptr_type mem_ptr_ = nullptr;

 public:
  using value_type = T;

  bundle_member_property_map(Graph* pg, member_ptr_type mem_ptr) : bmap_(pg), mem_ptr_(mem_ptr) {}
  bundle_member_property_map() = default;

  template <typename Key>
  decltype(auto) operator[](Key&& k) const {
    return bmap_[std::forward<Key>(k)].*mem_ptr_;
  }
};

// This function template can be used to construct a property-map that takes a
// vertex-descriptor into a graph and maps it to a data member of the vertex-bundle.
// This property-map is constructed from data-member map and a graph.
// `member_prop`: The data-member property-map.
// `g`: The graph.
// Returns a property-map that can map a vertex-descriptor to a vertex-bundle data-member.
template <typename MemberMap, typename Graph>
auto bundle_prop_to_vertex_prop(MemberMap member_prop, Graph& g) {
  return composite_property_map(member_prop, whole_bundle_property_map<Graph, vertex_bundle_t>(&g));
}

// This function template can be used to construct a property-map that takes a
// edge-descriptor into a graph and maps it to a data member of the edge-bundle.
// This property-map is constructed from data-member map and a graph.
// `member_prop`: The data-member property-map.
// `g`: The graph.
// Returns a property-map that can map a edge-descriptor to a edge-bundle data-member.
template <typename MemberMap, typename Graph>
auto bundle_prop_to_edge_prop(MemberMap member_prop, Graph& g) {
  return composite_property_map(member_prop, whole_bundle_property_map<Graph, edge_bundle_t>(&g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_MORE_PROPERTY_MAPS_H_
