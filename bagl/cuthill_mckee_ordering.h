// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Copyright 2004, 2005 Trustees of Indiana University
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek,
//          Doug Gregor, D. Kevin McGrath
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_CUTHILL_MCKEE_ORDERING_H_
#define BAGL_BAGL_CUTHILL_MCKEE_ORDERING_H_

#include <algorithm>

#include "bagl/depth_first_search.h"
#include "bagl/detail/sparse_ordering.h"
#include "bagl/graph_utility.h"
#include "bagl/properties.h"
#include "bagl/vector_property_map.h"

// (Reverse) Cuthill-McKee Algorithm for matrix reordering

namespace bagl {

namespace cuthill_mckee_detail {

template <typename OutputIterator, typename Buffer, typename DegreeMap>
class bfs_rcm_visitor : public default_bfs_visitor {
 public:
  bfs_rcm_visitor(OutputIterator* iter, Buffer* b, DegreeMap deg) : permutation_(iter), q_ptr_(b), degree_(deg) {}
  template <class Vertex, class Graph>
  void examine_vertex(Vertex u, Graph& /*unused*/) {
    *(*permutation_)++ = u;
    index_begin_ = q_ptr_->size();
  }
  template <class Vertex, class Graph>
  void finish_vertex(Vertex /*unused*/, Graph& /*unused*/) {
    std::sort(q_ptr_->begin() + index_begin_, q_ptr_->end(),
              [this](const auto& lhs, const auto& rhs) { return get(degree_, lhs) < get(degree_, rhs); });
  }

 protected:
  OutputIterator* permutation_ = nullptr;
  int index_begin_ = 0;
  Buffer* q_ptr_ = nullptr;
  DegreeMap degree_;
};

}  // namespace cuthill_mckee_detail

// Reverse Cuthill-McKee algorithm with a given starting Vertex.
//
// If user provides a reverse iterator, this will be a reverse-cuthill-mckee
// algorithm, otherwise it will be a standard CM algorithm
template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, concepts::ReadableVertexPropertyMap<G> DegreeMap>
OutputIterator cuthill_mckee_ordering(const G& g, std::deque<graph_vertex_descriptor_t<G>> vertex_queue,
                                      OutputIterator permutation, ColorMap color, DegreeMap degree) {
  using Vertex = graph_vertex_descriptor_t<G>;
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;

  sparse_detail::sparse_ordering_queue<Vertex> q;

  // create a bfs_rcm_visitor as defined above
  cuthill_mckee_detail::bfs_rcm_visitor vis(&permutation, &q, degree);

  // initialize the color map
  for (auto u : vertices(g)) {
    put(color, u, Color::white());
  }

  while (!vertex_queue.empty()) {
    Vertex s = vertex_queue.front();
    vertex_queue.pop_front();

    // call BFS with visitor
    breadth_first_visit(g, std::ranges::single_view{s}, q, vis, color);
  }
  return permutation;
}

// This is the case where only a single starting vertex is supplied.
template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, concepts::ReadableVertexPropertyMap<G> DegreeMap>
OutputIterator cuthill_mckee_ordering(const G& g, graph_vertex_descriptor_t<G> s, OutputIterator permutation,
                                      ColorMap color, DegreeMap degree) {
  std::deque<graph_vertex_descriptor_t<G>> vertex_queue;
  vertex_queue.push_front(s);
  return cuthill_mckee_ordering(g, std::move(vertex_queue), permutation, color, degree);
}

// This is the version of CM which selects its own starting vertex
template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, concepts::ReadableVertexPropertyMap<G> DegreeMap>
OutputIterator cuthill_mckee_ordering(const G& g, OutputIterator permutation, ColorMap color, DegreeMap degree) {
  if (vertices(g).empty()) {
    return permutation;
  }

  using Vertex = graph_vertex_descriptor_t<G>;
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;

  std::deque<Vertex> vertex_queue;

  // Mark everything white
  for (auto v : vertices(g)) {
    put(color, v, Color::white());
  }

  // Find one vertex from each connected component
  for (auto v : vertices(g)) {
    if (get(color, v) == Color::white()) {
      depth_first_visit(g, v, default_dfs_visitor(), color);
      vertex_queue.push_back(v);
    }
  }

  // Find starting nodes for all vertices
  // TBD: How to do this with a directed graph?
  for (auto& v : vertex_queue) {
    v = find_starting_node(g, v, color, degree);
  }

  return cuthill_mckee_ordering(g, std::move(vertex_queue), permutation, color, degree);
}

template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator,
          concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
OutputIterator cuthill_mckee_ordering(const G& g, OutputIterator permutation, VertexIndexMap index_map) {
  if (vertices(g).empty()) {
    return permutation;
  }

  return cuthill_mckee_ordering(g, permutation,
                                vector_property_map(num_vertices(g), index_map, default_color_type::white_color).ref(),
                                out_degree_property_map(g));
}

template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator>
OutputIterator cuthill_mckee_ordering(const G& g, OutputIterator permutation) {
  return cuthill_mckee_ordering(g, permutation, get(vertex_index, g));
}
}  // namespace bagl

#endif  // BAGL_BAGL_CUTHILL_MCKEE_ORDERING_H_
