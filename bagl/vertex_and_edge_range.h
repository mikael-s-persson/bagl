// Copyright 2004 The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_VERTEX_AND_EDGE_RANGE_H_
#define BAGL_BAGL_VERTEX_AND_EDGE_RANGE_H_

#include <ranges>

#include "bagl/graph_traits.h"

namespace bagl {

namespace vertex_and_edge_range_detail {
template <typename Graph, std::ranges::input_range VertexRange, std::ranges::input_range EdgeRange>
class vertex_and_edge_range {
  using traits_type = graph_traits<Graph>;

 public:
  using directed_category = typename traits_type::directed_category;
  using edge_parallel_category = typename traits_type::edge_parallel_category;
  struct traversal_category : public virtual vertex_list_graph_tag, public virtual edge_list_graph_tag {};

  using vertex_descriptor = std::ranges::range_value_t<VertexRange>;
  using edge_descriptor = std::ranges::range_value_t<EdgeRange>;

  static vertex_descriptor null_vertex() { return traits_type::null_vertex(); }

  vertex_and_edge_range(const Graph& g, VertexRange v_range, std::size_t n, EdgeRange e_range, std::size_t m)
      : g_(&g), v_range_(std::move(v_range)), m_num_vertices_(n), e_range_(std::move(e_range)), m_num_edges_(m) {}

  vertex_and_edge_range(const Graph& g, VertexRange v_range, EdgeRange e_range)
      : g_(&g),
        v_range_(std::move(v_range)),
        m_num_vertices_(std::ranges::distance(v_range_)),
        e_range_(std::move(e_range)),
        m_num_edges_(std::ranges::distance(e_range_)) {}

  const Graph* g_ = nullptr;
  VertexRange v_range_;
  std::size_t m_num_vertices_ = 0;
  EdgeRange e_range_;
  std::size_t m_num_edges_ = 0;
};

template <typename Graph, typename VertexRange, typename EdgeRange>
auto vertices(const vertex_and_edge_range<Graph, VertexRange, EdgeRange>& g) {
  return g.v_range_;
}

template <typename Graph, typename VertexRange, typename EdgeRange>
auto num_vertices(const vertex_and_edge_range<Graph, VertexRange, EdgeRange>& g) {
  return g.m_num_vertices_;
}

template <typename Graph, typename VertexRange, typename EdgeRange>
auto edges(const vertex_and_edge_range<Graph, VertexRange, EdgeRange>& g) {
  return g.e_range_;
}

template <typename Graph, typename VertexRange, typename EdgeRange>
auto num_edges(const vertex_and_edge_range<Graph, VertexRange, EdgeRange>& g) {
  return g.m_num_edges_;
}

template <typename Graph, typename VertexRange, typename EdgeRange>
auto source(graph_edge_descriptor_t<Graph> e, const vertex_and_edge_range<Graph, VertexRange, EdgeRange>& g) {
  return source(e, *g.g_);
}

template <typename Graph, typename VertexRange, typename EdgeRange>
auto target(graph_edge_descriptor_t<Graph> e, const vertex_and_edge_range<Graph, VertexRange, EdgeRange>& g) {
  return target(e, *g.g_);
}
}  // namespace vertex_and_edge_range_detail

// Make type alias so that free functions are not in top namespace.

template <typename Graph, std::ranges::input_range VertexRange, std::ranges::input_range EdgeRange>
using vertex_and_edge_range = vertex_and_edge_range_detail::vertex_and_edge_range<Graph, VertexRange, EdgeRange>;

}  // namespace bagl
#endif  // BAGL_BAGL_VERTEX_AND_EDGE_RANGE_H_
