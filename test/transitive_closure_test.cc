// Copyright (C) 2001 Vladimir Prus <ghost@cs.msu.su>
// Copyright (C) 2001 Jeremy Siek <jsiek@cs.indiana.edu>
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/transitive_closure.h"

#include <cstdlib>
#include <ctime>
#include <iostream>
#include <random>

#include "bagl/adjacency_list.h"
#include "bagl/graph_utility.h"
#include "bagl/properties.h"
#include "bagl/vector_as_graph.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

void generate_graph(std::mt19937& gen, int n, double p, std::vector<std::vector<int>>& r1) {
  std::uniform_real_distribution<> real_dist;
  r1.clear();
  r1.resize(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (real_dist(gen) < p) {
        r1[i].push_back(j);
      }
    }
  }
}

template <typename Graph>
std::size_t num_incident(graph_vertex_descriptor_t<Graph> u, graph_vertex_descriptor_t<Graph> v, const Graph& g) {
  std::size_t d = 0;
  for (auto e : out_edges(u, g)) {
    if (target(e, g) == v) {
      ++d;
    }
  }
  return d;
}

// (i,j) is in E' iff j is reachable from i
// Hmm, is_reachable does not detect when there is a non-trivial path
// from i to i. It always returns true for is_reachable(i,i).
// This needs to be fixed/worked around.
template <typename Graph, typename GraphTC>
bool check_transitive_closure(Graph& g, GraphTC& tc) {
  for (auto u : vertices(g)) {
    for (auto v : vertices(g)) {
      auto [e_g, g_has_edge] = edge(u, v, g);
      auto num_tc = num_incident(u, v, tc);
      if (u == v) {
        if (g_has_edge) {
          EXPECT_EQ(num_tc, 1) << "u: " << u;
        } else {
          bool can_reach = false;
          for (auto k : adjacent_vertices(u, g)) {
            auto color_map_vec =
                vector_property_map(num_vertices(g), get(vertex_index, g), default_color_type::white_color);
            if (is_reachable(k, u, g, color_map_vec.ref())) {
              can_reach = true;
              break;
            }
          }
          if (can_reach) {
            EXPECT_EQ(num_tc, 1) << "u: " << u;
          } else {
            EXPECT_EQ(num_tc, 0) << "u: " << u;
          }
        }
      } else {
        auto color_map_vec =
            vector_property_map(num_vertices(g), get(vertex_index, g), default_color_type::white_color);
        if (is_reachable(u, v, g, color_map_vec.ref())) {
          EXPECT_EQ(num_tc, 1) << "u: " << u << " v: " << v;
        } else {
          EXPECT_EQ(num_tc, 0) << "u: " << u << " v: " << v;
        }
      }
    }
  }
  return true;
}

TEST(TransitiveClosureTest, RandomGraphs) {
  std::mt19937 gen{42};
  std::uniform_int_distribution<> int_dist{1, 20};
  std::uniform_real_distribution<> real_dist;

  for (int i = 0; i < 20; ++i) {
    int n = int_dist(gen);
    double p = real_dist(gen);

    std::vector<std::vector<int>> g1;
    generate_graph(gen, n, p, g1);

    std::vector<std::vector<int>> g1_c(g1);

    vector_as_graph g1_g(g1);
    std::vector<std::vector<int>> g1_tc;
    vector_as_graph g1_tc_g(g1_tc);
    transitive_closure(g1_g, g1_tc_g, get(vertex_index, g1_g));

    EXPECT_TRUE(check_transitive_closure(g1_g, g1_tc_g)) << "Original graph was:\n"
                                                         << print_graph_to_string(g1_g) << "\nResult is:\n"
                                                         << print_graph_to_string(g1_tc_g);
  }
}

TEST(TransitiveClosureTest, Basic) {
  using graph_t = adjacency_list<>;
  std::vector<std::pair<std::size_t, std::size_t>> expected_edges = {{0, 2}, {1, 0}, {1, 2}, {1, 4},
                                                                     {3, 0}, {3, 2}, {4, 2}, {4, 3}};

  graph_t g(5, expected_edges);

  graph_t g_tc;
  transitive_closure(g, g_tc);

  std::vector<std::pair<std::size_t, std::size_t>> resulting_edges;
  auto vindex_tc = get(vertex_index, g_tc);
  for (auto e : edges(g_tc)) {
    resulting_edges.emplace_back(get(vindex_tc, source(e, g_tc)), get(vindex_tc, target(e, g_tc)));
  }
  EXPECT_THAT(resulting_edges, ::testing::IsSupersetOf(expected_edges));

  EXPECT_TRUE(check_transitive_closure(g, g_tc)) << "Original graph was:\n"
                                                 << print_graph_to_string(g) << "\nResult is:\n"
                                                 << print_graph_to_string(g_tc);
}

}  // namespace
}  // namespace bagl
