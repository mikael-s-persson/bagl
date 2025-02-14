// Copyright 2013 University of Warsaw.
// Authors: Piotr Wygocki
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/cycle_canceling.h"

#include "bagl/edmonds_karp_max_flow.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/properties.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "test/min_cost_max_flow_utils.h"

namespace bagl {
namespace {

TEST(PathAugmentation, Test) {
  using Graph = sample_graph::Graph;
  using Edge = graph_edge_descriptor_t<Graph>;
  Graph g;
  auto rev = vector_property_map(0, get(edge_index, g), Edge{});
  auto [s, t] = sample_graph::get_sample_graph_2(g, rev.ref());

  int N = num_vertices(g);
  auto idx = get(vertex_index, g);
  auto edge_pred = vector_property_map(N, idx, Edge{});

  edmonds_karp_max_flow(g, s, t, get(edge_capacity, g), get(edge_residual_capacity, g), rev.ref(),
                        two_bit_color_map(N, idx).ref(), edge_pred.ref());
  cycle_canceling(g, get(edge_weight, g), rev.ref(), get(edge_residual_capacity, g), edge_pred.ref());

  int cost = find_flow_cost(g);
  EXPECT_EQ(cost, 7);
}

TEST(CycleCanceling, Test) {
  using Graph = sample_graph::Graph;
  using Edge = graph_edge_descriptor_t<Graph>;
  Graph g;
  auto rev = vector_property_map(0, get(edge_index, g), Edge{});
  auto [s, t] = sample_graph::get_sample_graph(g, rev.ref());

  int N = num_vertices(g);
  auto idx = get(vertex_index, g);
  auto edge_pred = vector_property_map(N, idx, Edge{});

  edmonds_karp_max_flow(g, s, t, get(edge_capacity, g), get(edge_residual_capacity, g), rev.ref(),
                        two_bit_color_map(N, idx).ref(), edge_pred.ref());
  cycle_canceling(g, get(edge_weight, g), rev.ref(), get(edge_residual_capacity, g), edge_pred.ref());

  int cost = find_flow_cost(g);
  EXPECT_EQ(cost, 29);
}

}  // namespace
}  // namespace bagl
