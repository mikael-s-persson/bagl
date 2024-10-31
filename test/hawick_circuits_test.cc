// (C) Copyright 2013 Louis Dionne
// Modified from `tiernan_all_cycles_test.cc`.
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/hawick_circuits.h"

#include <cstddef>
#include <iostream>
#include <iterator>
#include <random>

#include "bagl/directed_graph.h"
#include "bagl/erdos_renyi_generator.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/undirected_graph.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

struct cycle_validator {
  explicit cycle_validator(std::size_t& number_of_cycles) : cycles(number_of_cycles) {}

  template <typename Path, typename Graph>
  void cycle(Path const& p, Graph const& g) {
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

template <typename Graph, typename Algorithm>
void test_one(Algorithm algorithm, std::size_t num_cycles_expected) {
  // Generate random graph with N vertices and probability P
  // of edge connection.
  static std::size_t const N = 20;
  static double const P = 0.1;

  std::mt19937 rng{42};
  Graph g(N, erdos_renyi_range<is_undirected_graph_v<Graph>>(rng, N, P));
  renumber_indices(g);

  std::size_t cycles = 0;
  cycle_validator vis(cycles);
  algorithm(g, vis);

  EXPECT_EQ(vis.cycles, num_cycles_expected);
}

template <typename Algorithm>
void cycle_test(Algorithm const& algorithm, std::size_t num_cycles_expected_undirected,
                std::size_t num_cycles_expected_directed) {
  test_one<undirected_graph<> >(algorithm, num_cycles_expected_undirected);

  test_one<directed_graph<> >(algorithm, num_cycles_expected_directed);
}

struct call_hawick_circuits {
  std::size_t max_length;
  call_hawick_circuits(std::size_t ml = 0) : max_length(ml) {}

  template <typename Graph, typename Visitor>
  void operator()(Graph const& g, Visitor& v) const {
    hawick_circuits(g, v, get(vertex_index, g), max_length);
  }
};

struct call_hawick_unique_circuits {
  std::size_t max_length;
  call_hawick_unique_circuits(std::size_t ml = 0) : max_length(ml) {}

  template <typename Graph, typename Visitor>
  void operator()(Graph const& g, Visitor& v) const {
    hawick_unique_circuits(g, v, get(vertex_index, g), max_length);
  }
};

TEST(HawickCircuitsTest, Basic) {
  // The last two arguments to cycle_test() are the expected (correct)
  // number of circuits in the undirected and directed test graphs.

  cycle_test(call_hawick_circuits(), 32, 96);
  cycle_test(call_hawick_unique_circuits(), 32, 96);

  // Correct values for max_length = 0 to 10
  // undirected
  std::size_t nc1[] = {32, 0, 20, 22, 24, 28, 28, 30, 32, 32, 32};
  // directed
  std::size_t nc2[] = {96, 0, 3, 5, 10, 17, 28, 36, 49, 62, 69};
  // undirected, unique
  std::size_t nc3[] = {32, 0, 20, 22, 24, 28, 28, 30, 32, 32, 32};
  // directed, unique
  std::size_t nc4[] = {96, 0, 3, 5, 10, 17, 28, 36, 49, 62, 69};
  for (unsigned int ml = 0; ml <= 10; ++ml) {
    cycle_test(call_hawick_circuits(ml), nc1[ml], nc2[ml]);
    cycle_test(call_hawick_unique_circuits(ml), nc3[ml], nc4[ml]);
  }
}

}  // namespace
}  // namespace bagl
