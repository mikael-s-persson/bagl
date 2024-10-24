// Copyright (C) 2002 Trustees of Indiana University
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/adjacency_list.h"
#include "bagl/dag_shortest_paths.h"
#include "bagl/graph_traits.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(DagLongestPaths, Basic) {
  using Graph = adjacency_list<vec_s, vec_s, directed_s, no_property, property<edge_weight_t, int> >;

  Graph g;

  auto v0 = add_vertex(g);
  auto v1 = add_vertex(g);
  auto v2 = add_vertex(g);
  auto v3 = add_vertex(g);

  auto [e, e_added] = add_edge(v0, v1, g);
  EXPECT_TRUE(e_added);
  put(edge_weight, g, e, 1);

  std::tie(e, e_added) = add_edge(v1, v2, g);
  EXPECT_TRUE(e_added);
  put(edge_weight, g, e, 1);

  std::tie(e, e_added) = add_edge(v3, v1, g);
  EXPECT_TRUE(e_added);
  put(edge_weight, g, e, 5);

  auto dist = vector_property_map(num_vertices(g), get(vertex_index, g), int{0});

  dag_shortest_paths(
      g, v0, dist.ref(), get(edge_weight, g), two_bit_color_map(num_vertices(g), get(vertex_index, g)).ref(),
      vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<Graph>::null_vertex()).ref(),
      default_dijkstra_visitor(), std::greater<>(), std::plus<>(), std::numeric_limits<int>::min(), 0);

  EXPECT_EQ(dist[v2], 2);
}

}  // namespace
}  // namespace bagl
