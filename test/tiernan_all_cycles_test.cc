// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/tiernan_all_cycles.h"

#include <cstddef>
#include <iostream>
#include <random>

#include "bagl/directed_graph.h"
#include "bagl/erdos_renyi_generator.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/undirected_graph.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

struct cycle_validator {
  cycle_validator(std::size_t& c) : cycles(c) {}

  template <typename Path, typename Graph>
  void cycle(const Path& p, const Graph& g) {
    ++cycles;
    // Check to make sure that each of the vertices in the path
    // is truly connected and that the back is connected to the
    // front - it's not validating that we find all paths, just
    // that the paths are valid.
    auto last = std::prev(p.end());
    for (auto i = p.begin(); i != last; ++i) {
      auto j = std::next(i);
      EXPECT_TRUE(edge(*i, *j, g).second);
    }
    EXPECT_TRUE(edge(p.back(), p.front(), g).second);
  }

  std::size_t& cycles;
};

struct cycle_printer {
  template <typename Path, typename Graph>
  void cycle(const Path& p, const Graph& g) {
    ++cycles;
    std::cout << "Cycle " << cycles << ": ";
    for (auto v : p) {
      std::cout << " " << get(vertex_index, g, v);
    }
    std::cout << std::endl;
  }
  std::size_t cycles = 0;
};

template <typename Graph>
void test() {
  using ERIter = erdos_renyi_iterator<std::mt19937, Graph>;

  // Generate random graph with N vertices and probability P
  // of edge connection.
  constexpr std::size_t n = 20;
  constexpr double p = 0.1;
  std::mt19937 rng(42);

  Graph g(n, std::ranges::subrange<ERIter, ERIter>(ERIter(rng, n, p), ERIter()));
  renumber_indices(g);

  std::size_t cycles = 0;
  cycle_validator vis(cycles);
  tiernan_all_cycles(g, vis);
  if constexpr (is_undirected_graph_v<Graph>) {
    EXPECT_EQ(vis.cycles, 12);
  } else {
    EXPECT_EQ(vis.cycles, 96);
  }
  if (::testing::Test::HasFailure()) {
    print_graph(g);
    cycle_printer print_vis{};
    tiernan_all_cycles(g, print_vis);
  }
}

TEST(TiernanAllCyclesTest, UndirectedGraph) {
  using Graph = undirected_graph<>;
  test<Graph>();
}

TEST(TiernanAllCyclesTest, DirectedGraph) {
  using Graph = directed_graph<>;
  test<Graph>();
}

}  // namespace
}  // namespace bagl
