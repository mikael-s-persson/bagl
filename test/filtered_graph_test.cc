// Copyright 1997, 1998, 1999, 2000, 2002 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// (c) Copyright Juergen Hunold 2012
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/filtered_graph.h"

#include "bagl/adjacency_list.h"
#include "bagl/dijkstra_shortest_paths.h"
#include "bagl/graph_archetypes.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

// Check filtered_graph
TEST(FilteredGraphTest, BasicGraphConcepts) {
  using Graph = adjacency_list<vec_s, vec_s, directed_s, no_property, property<edge_residual_capacity_t, int>>;
  using ResCapMap = property_map_t<Graph, edge_residual_capacity_t>;
  using ResGraph = filtered_graph<Graph, is_residual_edge<ResCapMap>>;
  using Edge = graph_edge_descriptor_t<ResGraph>;

  static_assert(concepts::VertexListGraph<ResGraph>);
  static_assert(concepts::EdgeListGraph<ResGraph>);
  static_assert(concepts::IncidenceGraph<ResGraph>);
  static_assert(concepts::AdjacencyGraph<ResGraph>);
  static_assert(concepts::PropertyGraph<ResGraph, Edge, edge_residual_capacity_t>);
  static_assert(!is_undirected_graph_v<ResGraph>);
}
// Check filtered_graph with bidirectional adjacency_list
TEST(FilteredGraphTest, BidirectionalGraphConcepts) {
  using Graph = adjacency_list<vec_s, vec_s, bidirectional_s, no_property, property<edge_residual_capacity_t, int>>;
  using ResCapMap = property_map_t<Graph, edge_residual_capacity_t>;
  using ResGraph = filtered_graph<Graph, is_residual_edge<ResCapMap>>;
  static_assert(concepts::BidirectionalGraph<ResGraph>);
  static_assert(!is_undirected_graph_v<ResGraph>);
}
// Check filtered_graph with undirected adjacency_list
TEST(FilteredGraphTest, UndirectedGraphConcepts) {
  using Graph = adjacency_list<vec_s, vec_s, undirected_s, no_property, property<edge_residual_capacity_t, int>>;
  using ResCapMap = property_map_t<Graph, edge_residual_capacity_t>;
  using ResGraph = filtered_graph<Graph, is_residual_edge<ResCapMap>>;
  static_assert(concepts::BidirectionalGraph<ResGraph>);
  static_assert(is_undirected_graph_v<ResGraph>);
}

TEST(FilteredGraphTest, MakeFilteredGraph) {
  adjacency_list<> g;
  auto fg = make_filtered_graph(g, keep_all(), [](auto) { return true; });
  static_assert(concepts::VertexListGraph<decltype(fg)>);
}

TEST(FilteredGraphTest, Dijkstra) {
  using Graph = adjacency_list<set_s, vec_s, undirected_s, no_property, property<edge_weight_t, double>>;

  Graph g(5);
  add_edge(0, 1, g, 0.1);
  add_edge(0, 2, g, 1.0);
  add_edge(1, 2, g, 1.0);
  add_edge(0, 3, g, -0.2);
  add_edge(1, 3, g, 0.2);
  add_edge(3, 2, g, -0.1);
  add_edge(3, 4, g, 0.1);
  add_edge(0, 4, g, -0.5);
  auto node_src = *vertices(g).begin();
  auto pred = vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<Graph>::null_vertex());
  auto dist = vector_property_map(num_vertices(g), get(vertex_index, g), double{0});

  auto fg = make_filtered_graph(g, [&g](auto e) { return get(edge_weight, g, e) > 0.0; });
  dijkstra_shortest_paths(fg, std::ranges::single_view{node_src}, pred.ref(), dist.ref(), get(edge_weight, fg),
                          get(vertex_index, fg), default_dijkstra_visitor());
  EXPECT_THAT(dist[0], ::testing::DoubleEq(0.0));
  EXPECT_EQ(pred[0], 0);
  EXPECT_THAT(dist[1], ::testing::DoubleEq(0.1));
  EXPECT_EQ(pred[1], 0);
  EXPECT_THAT(dist[2], ::testing::DoubleEq(1.0));
  EXPECT_EQ(pred[2], 0);
  EXPECT_THAT(dist[3], ::testing::DoubleEq(0.3));
  EXPECT_EQ(pred[3], 1);
  EXPECT_THAT(dist[4], ::testing::DoubleEq(0.4));
  EXPECT_EQ(pred[4], 3);
}

}  // namespace
}  // namespace bagl
