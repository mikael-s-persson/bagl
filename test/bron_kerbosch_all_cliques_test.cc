// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/bron_kerbosch_all_cliques.h"

#include <algorithm>
#include <iostream>
#include <iterator>
#include <map>
#include <random>
#include <vector>

#include "bagl/directed_graph.h"
#include "bagl/erdos_renyi_generator.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/undirected_graph.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

// TODO: This is probably not a very good test. We should be able to define
// an identity test - something like find the clique of K5.

struct clique_validator {
  template <typename Clique, typename Graph>
  void clique(const Clique& c, const Graph& g) const {
    // Simply assert that each vertex in the clique is connected
    // to all others in the clique.
    auto c_end = c.end();
    for (auto i = c.begin(); i != c_end; ++i) {
      for (auto j = c.begin(); j != c_end; ++j) {
        if (i != j) {
          EXPECT_TRUE(edge(*i, *j, g).second);
        }
      }
    }
  }
};

template <typename Graph>
void test() {
  // Generate random graph with N vertices and probability P
  // of edge connection.
  constexpr std::size_t n = 20;
  constexpr double p = 0.1;

  std::mt19937 rng(42);
  Graph g(n, erdos_renyi_range<is_undirected_graph_v<Graph>>(rng, n, p));
  renumber_indices(g);
  print_graph(g);

  bron_kerbosch_all_cliques(g, clique_validator{});
}

TEST(BronKerboschAllCliques, UndirectedGraph) {
  using Graph = undirected_graph<>;
  test<Graph>();
}

TEST(BronKerboschAllCliques, DirectedGraph) {
  using Graph = directed_graph<>;
  test<Graph>();
}

}  // namespace
}  // namespace bagl
