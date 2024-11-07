// Copyright 2013 University of Warsaw.
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/successive_shortest_path_nonnegative_weights.h"

#include "bagl/find_flow_cost.h"
#include "bagl/properties.h"
#include "test/min_cost_max_flow_utils.h"

namespace bagl {
namespace {

TEST(SuccShortestPathAugmentation, Test) {
  using Graph = sample_graph::Graph;
  using Edge = graph_edge_descriptor_t<Graph>;
  Graph g;
  auto rev = vector_property_map(0, get(edge_index, g), Edge{});
  auto [s, t] = sample_graph::get_sample_graph(g, rev.ref());

  int N = num_vertices(g);
  auto idx = get(vertex_index, g);
  auto edge_pred = vector_property_map(N, idx, Edge{});

  successive_shortest_path_nonnegative_weights(g, s, t, get(edge_capacity, g), get(edge_residual_capacity, g),
                                               get(edge_weight, g), rev.ref(), idx, edge_pred.ref());

  int cost = find_flow_cost(g);
  EXPECT_EQ(cost, 29);
}

TEST(SuccShortestPathAugmentation, Test2) {
  using Graph = sample_graph::Graph;
  using Edge = graph_edge_descriptor_t<Graph>;
  Graph g;
  auto rev = vector_property_map(0, get(edge_index, g), Edge{});
  auto [s, t] = sample_graph::get_sample_graph_2(g, rev.ref());

  int N = num_vertices(g);
  auto idx = get(vertex_index, g);
  auto edge_pred = vector_property_map(N, idx, Edge{});

  successive_shortest_path_nonnegative_weights(g, s, t, get(edge_capacity, g), get(edge_residual_capacity, g),
                                               get(edge_weight, g), rev.ref(), idx, edge_pred.ref());

  int cost = find_flow_cost(g);
  EXPECT_EQ(cost, 7);
}

}  // namespace
}  // namespace bagl
