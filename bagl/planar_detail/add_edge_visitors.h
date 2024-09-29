// Copyright 2007 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PLANAR_DETAIL_ADD_EDGE_VISITORS_H_
#define BAGL_BAGL_PLANAR_DETAIL_ADD_EDGE_VISITORS_H_

#include "bagl/property_map.h"

namespace bagl {

struct default_add_edge_visitor {
  template <typename Graph, typename Vertex>
  void visit_vertex_pair(Vertex u, Vertex v, Graph& g) {
    add_edge(u, v, g);
  }
};

template <typename EdgeIndexMap>
struct edge_index_update_visitor {
  using edge_index_value_t = property_traits_value_t<EdgeIndexMap>;

  edge_index_update_visitor(EdgeIndexMap em, edge_index_value_t next_index_available)
      : em_(em), next_index_(next_index_available) {}

  template <typename Graph, typename Vertex>
  void visit_vertex_pair(Vertex u, Vertex v, Graph& g) {
    auto [e, added] = add_edge(u, v, g);
    if (added) {
      put(em_, e, next_index_++);
    }
  }

 private:
  EdgeIndexMap em_;
  edge_index_value_t next_index_;
};

}  // namespace bagl

#endif  // BAGL_BAGL_PLANAR_DETAIL_ADD_EDGE_VISITORS_H_
