// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/closeness_centrality.h"

#include <iostream>
#include <vector>

#include "bagl/constant_property_map.h"
#include "bagl/directed_graph.h"
#include "bagl/floyd_warshall_shortest.h"
#include "bagl/graph_traits.h"
#include "bagl/undirected_graph.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

// number of vertices in the graph
constexpr unsigned N = 5;

template <typename Graph>
void build_graph(Graph& g, std::vector<graph_vertex_descriptor_t<Graph>>& v) {
  // add vertices
  for (auto& u : v) {
    u = add_vertex(g);
  }

  // add edges
  ASSERT_TRUE(add_edge(v[0], v[1], g).second);
  ASSERT_TRUE(add_edge(v[1], v[2], g).second);
  ASSERT_TRUE(add_edge(v[2], v[0], g).second);
  ASSERT_TRUE(add_edge(v[3], v[4], g).second);
  ASSERT_TRUE(add_edge(v[4], v[0], g).second);
}

TEST(ClosenessCentrality, UndirectedGraph) {
  using Graph = undirected_graph<>;
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;

  Graph g;
  std::vector<Vertex> v(N);
  build_graph(g, v);

  auto cm = vector_property_map(num_vertices(g), get(vertex_index, g), double{0.0});
  auto dm = vector_property_map(num_vertices(g), get(vertex_index, g),
                                vector_property_map(num_vertices(g), get(vertex_index, g), int{0}));
  for (auto u : vertices(g)) {
    dm[u] = vector_property_map(num_vertices(g), get(vertex_index, g), int{0});
  }

  auto wm = constant_property_map<Edge, int>(1);

  floyd_warshall_all_pairs_shortest_paths(g, dm, wm);
  all_closeness_centralities(g, dm, cm);

  EXPECT_THAT(cm[v[0]], ::testing::DoubleEq(1.0 / 5.0));
  EXPECT_THAT(cm[v[1]], ::testing::DoubleEq(1.0 / 7.0));
  EXPECT_THAT(cm[v[2]], ::testing::DoubleEq(1.0 / 7.0));
  EXPECT_THAT(cm[v[3]], ::testing::DoubleEq(1.0 / 9.0));
  EXPECT_THAT(cm[v[4]], ::testing::DoubleEq(1.0 / 6.0));
}

TEST(ClosenessCentrality, DirectedGraph) {
  using Graph = directed_graph<>;
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;

  Graph g;
  std::vector<Vertex> v(N);
  build_graph(g, v);

  auto cm = vector_property_map(num_vertices(g), get(vertex_index, g), double{0.0});
  auto dm = vector_property_map(num_vertices(g), get(vertex_index, g),
                                vector_property_map(num_vertices(g), get(vertex_index, g), int{0}));
  for (auto u : vertices(g)) {
    dm[u] = vector_property_map(num_vertices(g), get(vertex_index, g), int{0});
  }

  auto wm = constant_property_map<Edge, int>(1);

  floyd_warshall_all_pairs_shortest_paths(g, dm, wm);
  all_closeness_centralities(g, dm, cm);

  EXPECT_THAT(cm[v[0]], ::testing::DoubleEq(0.0));
  EXPECT_THAT(cm[v[1]], ::testing::DoubleEq(0.0));
  EXPECT_THAT(cm[v[2]], ::testing::DoubleEq(0.0));
  EXPECT_THAT(cm[v[3]], ::testing::DoubleEq(1.0 / 10.0));
  EXPECT_THAT(cm[v[4]], ::testing::DoubleEq(0.0));
}

}  // namespace
}  // namespace bagl
