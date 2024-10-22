// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_TOPOLOGICAL_SORT_H_
#define BAGL_BAGL_TOPOLOGICAL_SORT_H_

#include "bagl/depth_first_search.h"
#include "bagl/exception.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"
#include "bagl/visitors.h"

namespace bagl {

// Topological sort visitor
//
// This visitor merely writes the linear ordering into an
// OutputIterator. The OutputIterator could be something like an
// ostream_iterator, or it could be a back/front_insert_iterator.
// Note that if it is a back_insert_iterator, the recorded order is
// the reverse topological order. On the other hand, if it is a
// front_insert_iterator, the recorded order is the topological
// order.
//
template <typename OutputIterator>
struct topo_sort_visitor : dfs_visitor<> {
  explicit topo_sort_visitor(OutputIterator iter) : iter_(iter) {}

  template <typename Edge, typename Graph>
  void back_edge(const Edge&, Graph&) {
    throw not_a_dag();
  }

  template <typename Vertex, typename Graph>
  void finish_vertex(const Vertex& u, Graph& /*unused*/) {
    *iter_++ = u;
  }

  OutputIterator iter_;
};

// Topological Sort
//
// The topological sort algorithm creates a linear ordering
// of the vertices such that if edge (u,v) appears in the graph,
// then u comes before v in the ordering. The graph must
// be a directed acyclic graph (DAG). The implementation
// consists mainly of a call to depth-first search.
//
template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator>
void topological_sort(G& g, OutputIterator result) {
  depth_first_search(g, topo_sort_visitor(result), two_bit_color_map(num_vertices(g), get(vertex_index, g)).ref());
}

template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void topological_sort(const G& g, OutputIterator result, ColorMap color) {
  depth_first_search(g, topo_sort_visitor(result), color);
}

template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator,
          concepts::ReadableVertexIndexMap<G> VertexIndexMap>
void topological_sort(const G& g, OutputIterator result, VertexIndexMap v_index) {
  depth_first_search(g, topo_sort_visitor(result), two_bit_color_map(num_vertices(g), v_index).ref());
}

}  // namespace bagl

#endif  // BAGL_BAGL_TOPOLOGICAL_SORT_H_
