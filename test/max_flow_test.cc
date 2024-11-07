//  Copyright (c) 2006, Stephan Diederich
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <algorithm>
#include <iostream>
#include <iterator>
#include <random>
#include <unordered_map>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/boykov_kolmogorov_max_flow.h"
#include "bagl/directed_graph.h"
#include "bagl/edmonds_karp_max_flow.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"
#include "bagl/push_relabel_max_flow.h"
#include "bagl/random.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "test/min_cost_max_flow_utils.h"

/***************
 * test which compares results of the three different max_flow implementations
 * command line parameters:
 *   number_of_vertices: defaults to 100
 *   number_of_edges:    defaults to 1000
 *   seeed:              defaults to 1
 ***************/

namespace bagl {
namespace {

TEST(MaxFlowTest, RandomGraph) {
  using Graph = adjacency_list<vec_s, vec_s, directed_s,
                               property<vertex_color_t, default_color_type, property<vertex_distance_t, int> >,
                               property<edge_capacity_t, int, property<edge_residual_capacity_t, int> > >;
  using Edge = graph_edge_descriptor_t<Graph>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  std::size_t n_verts = 100;
  std::size_t n_edges = 1000;
  std::size_t seed = 42;

  Graph g;
  const int cap_low = 1;
  const int cap_high = 1000;

  // init random numer generator
  std::mt19937 gen(seed);
  // generate graph
  generate_random_graph(g, n_verts, n_edges, gen);

  std::uniform_int_distribution<> int_dist(cap_low, cap_high);
  // init edge-capacities
  randomize_property<edge_capacity_t>(g, [&]() { return int_dist(gen); });

  // get source and sink node
  Vertex u = random_vertex(g, gen);
  Vertex v = random_vertex(g, gen);
  while (v == u) {
    v = random_vertex(g, gen);
  }

  // add reverse edges (ugly... how to do better?!)
  std::unordered_map<Edge, Edge> rev_cont;
  auto cap = get(edge_capacity, g);
  std::vector<Edge> edges_copy;
  std::ranges::copy(edges(g), std::back_inserter(edges_copy));
  while (!edges_copy.empty()) {
    Edge old_edge = edges_copy.back();
    edges_copy.pop_back();
    auto [new_edge, inserted] = add_edge(target(old_edge, g), source(old_edge, g), g);
    ASSERT_TRUE(inserted);
    rev_cont[old_edge] = new_edge;
    rev_cont[new_edge] = old_edge;
    cap[new_edge] = 0;
  }

  auto pred = vector_property_map(num_vertices(g), get(vertex_index, g), Edge{});
  auto rev = associative_property_map(rev_cont);

  auto bk = boykov_kolmogorov_max_flow(g, get(edge_capacity, g), get(edge_residual_capacity, g), rev, pred.ref(),
                                       get(vertex_color, g), get(vertex_distance, g), get(vertex_index, g), u, v);
  auto push_relabel =
      push_relabel_max_flow(g, u, v, get(edge_capacity, g), get(edge_residual_capacity, g), rev, get(vertex_index, g));
  auto edmonds_karp = edmonds_karp_max_flow(g, u, v, get(edge_capacity, g), get(edge_residual_capacity, g), rev,
                                            get(vertex_color, g), pred.ref());

  EXPECT_EQ(bk, push_relabel);
  EXPECT_EQ(push_relabel, edmonds_karp);
}

struct TestEdgeProp {
  int capacity;
  int residual_capacity;
  int weight;
  std::size_t id;
};
struct TestNodeProp {
  float dist;
  default_color_type color;
};

TEST(MaxFlowTest, BundledGraph) {
  using Graph = adjacency_list<list_s, vec_s, directed_s, TestNodeProp, TestEdgeProp>;
  using Edge = graph_edge_descriptor_t<Graph>;

  Graph g;

  auto capacity = get(&TestEdgeProp::capacity, g);
  auto residual_capacity = get(&TestEdgeProp::residual_capacity, g);
  auto col = get(&TestNodeProp::color, g);

  auto pred = vector_property_map(num_vertices(g), get(vertex_index, g), Edge{});
  std::unordered_map<Edge, Edge> rev_cont;
  auto rev = associative_property_map(rev_cont);

  auto [s, t] = sample_graph::get_sample_graph(g, rev);

  auto flow_value = edmonds_karp_max_flow(g, s, t, capacity, residual_capacity, rev, col, pred.ref());

  EXPECT_EQ(flow_value, 4);
}

}  // namespace
}  // namespace bagl
