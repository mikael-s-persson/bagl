//  (C) Copyright Jeremy Siek 2004
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <cstdlib>
#include <iostream>
#include <random>
#include <sstream>

#include "bagl/adjacency_list.h"
#include "bagl/bellman_ford_shortest_paths.h"
#include "bagl/dijkstra_shortest_paths.h"
#include "bagl/graph_utility.h"
#include "bagl/random.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "test/new_york_state_map.h"

namespace bagl {
namespace {

TEST(BellmanFordSearchTest, Search) {
  // specify some types
  using Graph = test::state_map_type;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  std::mt19937 gen(42);

  Graph g = test::make_new_york_state_map();

  auto name = get(vertex_name, g);
  auto idx = get(vertex_index, g);
  auto weight = get(edge_weight, g);

  auto pm = vector_property_map(num_vertices(g), idx, graph_traits<Graph>::null_vertex());
  auto dm = vector_property_map(num_vertices(g), idx, std::numeric_limits<float>::max());

  for (int i = 0; i < 5; ++i) {
    // pick random start/goal
    Vertex start = gen() % num_vertices(g);

    auto print_problem_def = [&]() {
      std::stringstream out;
      print_graph(g, get(vertex_bundle, g), weight, out);
      out << "\nStart vertex: " << name[start];
      return out.str();
    };
    // Print problem definition if any checks fail.
    SCOPED_TRACE(print_problem_def());

    bellman_ford_shortest_paths(g, default_bellman_visitor(), start, num_vertices(g), weight, pm, dm);
    for (auto v : vertices(g)) {
      if (pm[v] == v) {
        EXPECT_EQ(v, start);
        continue;
      }
      EXPECT_NE(dm[v], std::numeric_limits<float>::max());
      // Recorded predecessor is the best.
      for (auto e : out_edges(v, g)) {
        EXPECT_GE(dm[target(e, g)] + get(weight, e), dm[v] - 1e-5F);
      }
    }
  }
}

TEST(BellmanFordSearchTest, CompareToDijkstra) {
  // specify some types
  using Graph = adjacency_list<list_s, vec_s, directed_s, no_property, property<edge_weight_t, float>>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  std::mt19937 gen(42);

  for (int i = 0; i < 5; ++i) {
    Graph g;
    generate_random_graph(g, 10, 50, gen);

    auto idx = get(vertex_index, g);
    auto weight = get(edge_weight, g);
    for (auto e : edges(g)) {
      std::uniform_real_distribution<float> uni_dist;
      put(weight, e, uni_dist(gen) + 1e-5F);
    }

    auto pm = vector_property_map(num_vertices(g), idx, graph_traits<Graph>::null_vertex());
    auto dm = vector_property_map(num_vertices(g), idx, default_inf_v<float>);
    auto djpm = vector_property_map(num_vertices(g), idx, graph_traits<Graph>::null_vertex());
    auto djdm = vector_property_map(num_vertices(g), idx, default_inf_v<float>);

    // pick random start/goal
    Vertex start = gen() % num_vertices(g);

    auto print_problem_def = [&]() {
      std::stringstream out;
      print_graph(g, get(vertex_index, g), weight, out);
      out << "\nStart vertex: " << start;
      return out.str();
    };
    // Print problem definition if any checks fail.
    SCOPED_TRACE(print_problem_def());

    // call astar named parameter interface

    bellman_ford_shortest_paths(g, default_bellman_visitor(), start, num_vertices(g), weight, pm, dm);
    dijkstra_shortest_paths(g, start, djpm, djdm, weight, idx, default_dijkstra_visitor());

    for (auto u : vertices(g)) {
      EXPECT_EQ(pm[u], djpm[u]);
      EXPECT_THAT(dm[u], ::testing::FloatNear(djdm[u], 1.0e-5F));
    }
  }
}

}  // namespace
}  // namespace bagl
