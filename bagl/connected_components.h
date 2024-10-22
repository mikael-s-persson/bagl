// Copyright 1997-2001 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_CONNECTED_COMPONENTS_H_
#define BAGL_BAGL_CONNECTED_COMPONENTS_H_

#include "bagl/depth_first_search.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/two_bit_color_map.h"

namespace bagl {

namespace connected_components_detail {

// This visitor is used both in the connected_components algorithm
// and in the kosaraju strong components algorithm during the
// second DFS traversal.
template <class ComponentsMap>
class components_recorder : public dfs_visitor<> {
  using CompType = property_traits_value_t<ComponentsMap>;

 public:
  components_recorder(ComponentsMap c, CompType& c_count) : component_(c), count_(c_count) {}

  template <class Vertex, class Graph>
  void start_vertex(Vertex /*unused*/, Graph& /*unused*/) {
    if (count_ == std::numeric_limits<CompType>::max()) {
      count_ = 0;  // start counting components at zero
    } else {
      ++count_;
    }
  }
  template <class Vertex, class Graph>
  void discover_vertex(Vertex u, Graph& /*unused*/) {
    put(component_, u, count_);
  }

 protected:
  ComponentsMap component_;
  CompType& count_;
};

}  // namespace connected_components_detail

// This function computes the connected components of an undirected
// graph using a single application of depth first search.

template <concepts::VertexListGraph G, concepts::WritableVertexPropertyMap<G> ComponentMap,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap>
property_traits_value_t<ComponentMap> connected_components(const G& g, ComponentMap c, ColorMap color) {
  static_assert(is_undirected_graph_v<G>);
  if (num_vertices(g) == 0) {
    return 0;
  }

  using CompType = property_traits_value_t<ComponentMap>;
  // c_count initialized to "nil" (with nil represented by (max)())
  CompType c_count = std::numeric_limits<CompType>::max();
  connected_components_detail::components_recorder vis(c, c_count);
  depth_first_search(g, vis, color);
  return c_count + 1;
}

template <concepts::VertexListGraph G, concepts::WritableVertexPropertyMap<G> ComponentMap>
property_traits_value_t<ComponentMap> connected_components(const G& g, ComponentMap c) {
  return connected_components(g, c, two_bit_color_map(num_vertices(g), get(vertex_index, g)).ref());
}

}  // namespace bagl

#endif  // BAGL_BAGL_CONNECTED_COMPONENTS_H_
