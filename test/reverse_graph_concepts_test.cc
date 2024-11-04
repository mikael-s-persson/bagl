// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <string>

#include "bagl/adjacency_list.h"
#include "bagl/graph_archetypes.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/reverse_graph.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

// Check const reverse_graph
TEST(ReverseGraphTest, ConstTest) {
  using AdjList = adjacency_list<vec_s, vec_s, bidirectional_s, property<vertex_color_t, int>,
                                 property<edge_weight_t, int>, property<graph_name_t, std::string> >;
  using Graph = reverse_graph<AdjList>;
  static_assert(concepts::VertexListGraph<Graph>);
  static_assert(concepts::IncidenceGraph<Graph>);
  static_assert(concepts::BidirectionalGraph<Graph>);
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;
  static_assert(concepts::ReadablePropertyGraph<Graph, Vertex, vertex_color_t>);
  static_assert(concepts::ReadablePropertyGraph<Graph, Edge, edge_weight_t>);
  static_assert(concepts::ReadablePropertyGraph<Graph, Edge, edge_underlying_t>);
  AdjList g;
  Graph gr(g);
  get_property(gr, graph_name);
}

// Check non-const reverse_graph
TEST(ReverseGraphTest, NonConstTest) {
  using AdjList = adjacency_list<vec_s, vec_s, bidirectional_s, property<vertex_color_t, int>,
                                 property<edge_weight_t, int>, property<graph_name_t, std::string> >;
  using Graph = reverse_graph<AdjList, AdjList&>;
  static_assert(concepts::VertexListGraph<Graph>);
  static_assert(concepts::IncidenceGraph<Graph>);
  static_assert(concepts::BidirectionalGraph<Graph>);
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;
  static_assert(concepts::PropertyGraph<Graph, Vertex, vertex_color_t>);
  static_assert(concepts::PropertyGraph<Graph, Edge, edge_weight_t>);
  static_assert(concepts::ReadablePropertyGraph<Graph, Edge, edge_underlying_t>);
  AdjList g;
  Graph gr(g);
  get_property(gr, graph_name_t());
  set_property(gr, graph_name_t(), "foo");
}

}  // namespace
}  // namespace bagl
