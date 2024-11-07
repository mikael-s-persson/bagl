// Copyright (c) 2005 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <cstdlib>
#include <iostream>
#include <random>

#include "bagl/adjacency_list.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/max_cardinality_matching.h"
#include "bagl/random.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(MatchingTest, RandomGraph) {
  using Graph = adjacency_list<vec_s, vec_s, undirected_s>;

  int n = 10;
  int m = 30;

  Graph g(n);

  std::mt19937 gen(42);
  while (num_edges(g) < m) {
    auto u = random_vertex(g, gen);
    auto v = random_vertex(g, gen);
    while (u == v || edge(u, v, g).second) {
      v = random_vertex(g, gen);
    }
    auto [e, e_added] = add_edge(u, v, g);
    EXPECT_TRUE(e_added);
  }

  auto mate = vector_property_map(n, get(vertex_index, g), graph_traits<Graph>::null_vertex());
  EXPECT_TRUE(checked_edmonds_maximum_cardinality_matching(g, mate.ref()));

  if (::testing::Test::HasFailure()) {
    std::cerr << "\nMatching: ";
    for (auto v : vertices(g)) {
      if (mate[v] != graph_traits<Graph>::null_vertex() && v < mate[v]) {
        std::cerr << "{" << v << "," << mate[v] << "}, ";
      }
    }
  }

  // Now remove an edge from the random_mate matching.
  for (auto v : vertices(g)) {
    if (mate[v] != graph_traits<Graph>::null_vertex()) {
      mate[mate[v]] = graph_traits<Graph>::null_vertex();
      mate[v] = graph_traits<Graph>::null_vertex();
      break;
    }
  }

  //...and run the matching verifier - it should tell us that the matching
  // isn't a maximum matching.
  EXPECT_FALSE(maximum_cardinality_matching_verify(g, mate.ref(), get(vertex_index, g)));

  // find a greedy matching on the graph
  auto greedy_mate = vector_property_map(n, get(vertex_index, g), graph_traits<Graph>::null_vertex());
  greedy_matching(g, greedy_mate.ref());

  if (matching_size(g, mate.ref()) > matching_size(g, greedy_mate.ref())) {
    EXPECT_FALSE(maximum_cardinality_matching_verify(g, greedy_mate.ref(), get(vertex_index, g)));
  }

  if (::testing::Test::HasFailure()) {
    std::cerr << "Graph:\n";
    print_graph(g, std::cerr);

    std::cerr << "\nGreedy matching: ";
    for (auto v : vertices(g)) {
      if (greedy_mate[v] != graph_traits<Graph>::null_vertex() && v < greedy_mate[v]) {
        std::cerr << "{" << v << "," << greedy_mate[v] << "}, ";
      }
    }
  }
}

}  // namespace
}  // namespace bagl
