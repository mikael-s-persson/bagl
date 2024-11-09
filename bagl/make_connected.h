// Copyright 2007 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_MAKE_CONNECTED_H_
#define BAGL_BAGL_MAKE_CONNECTED_H_

#include <iterator>
#include <tuple>
#include <vector>

#include "bagl/connected_components.h"
#include "bagl/planar_detail/add_edge_visitors.h"
#include "bagl/planar_detail/bucket_sort.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"

namespace bagl {

template <concepts::VertexListGraph G, concepts::ReadableVertexIndexMap<G> VertexIndexMap, typename AddEdgeVisitor>
void make_connected(G& g, VertexIndexMap vm, AddEdgeVisitor& vis) {
  using vertex_t = graph_vertex_descriptor_t<G>;

  auto component = vector_property_map(num_vertices(g), vm, std::size_t{});
  std::vector<vertex_t> vertices_by_component(num_vertices(g));

  std::size_t num_components = connected_components(g, component.ref());

  if (num_components < 2) {
    return;
  }

  auto v_rg = vertices(g);
  std::copy(v_rg.begin(), v_rg.end(), vertices_by_component.begin());

  bucket_sort(vertices_by_component.begin(), vertices_by_component.end(), component.ref());

  auto ci_end = vertices_by_component.end();
  auto ci_prev = vertices_by_component.begin();
  if (ci_prev == ci_end) {
    return;
  }

  for (auto ci = std::next(ci_prev); ci != ci_end; ci_prev = ci, ++ci) {
    if (component[*ci_prev] != component[*ci]) {
      vis.visit_vertex_pair(*ci_prev, *ci, g);
    }
  }
}

template <concepts::VertexListGraph G, concepts::ReadableVertexIndexMap<G> VertexIndexMap>
void make_connected(G& g, VertexIndexMap vm) {
  default_add_edge_visitor vis;
  make_connected(g, vm, vis);
}

template <concepts::VertexListGraph G>
void make_connected(G& g) {
  make_connected(g, get(vertex_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_MAKE_CONNECTED_H_
