// Copyright (c) 2005 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <cstdint>
#include <iostream>
#include <random>
#include <tuple>
#include <type_traits>
#include <utility>

#include "bagl/adjacency_list.h"
#include "bagl/adjacency_matrix.h"
#include "bagl/max_cardinality_matching.h"
#include "bagl/random.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

using VecGraph = adjacency_list<vec_s, vec_s, undirected_s>;

using ListGraph = adjacency_list<list_s, list_s, undirected_s, property<vertex_index_t, std::size_t>>;

using MatrixGraph = adjacency_matrix<undirected_s>;

template <typename Graph>
void vertex_index_install(Graph&) {}

void vertex_index_install(ListGraph& g) {
  std::size_t i = 0;
  for (auto v : vertices(g)) {
    put(vertex_index, g, v, i++);
  }
}

template <typename T>
class MatchingTest : public ::testing::Test {
 public:
  using Graph = std::tuple_element_t<0, T>;
  static constexpr std::size_t num_v = std::tuple_element_t<1, T>::value;
  static constexpr std::size_t double_num_v = num_v * 2;

  static void complete_graph(Graph& g, int n) {
    // creates the complete graph on n vertices
    g = Graph(n);
    auto v_rg = vertices(g);
    for (auto vi = v_rg.begin(); vi != v_rg.end(); ++vi) {
      for (auto wi = std::next(vi); wi != v_rg.end(); ++wi) {
        EXPECT_TRUE(add_edge(*vi, *wi, g).second);
      }
    }
  }

  static void gabows_graph(Graph& g, int n) {
    // creates a graph with 2n vertices, consisting of the complete graph
    // on n vertices plus n vertices of degree one, each adjacent to one
    // vertex in the complete graph. without any initial matching, this
    // graph forces edmonds' algorithm into worst-case behavior.

    g = Graph(2 * n);

    auto v_rg = vertices(g);

    auto halfway = v_rg.begin();
    for (int i = 0; i < n; ++i) {
      ++halfway;
    }

    for (auto ui = v_rg.begin(); ui != halfway; ++ui) {
      for (auto vi = std::next(ui); vi != halfway; ++vi) {
        EXPECT_TRUE(add_edge(*ui, *vi, g).second);
      }
    }

    for (auto ui = v_rg.begin(); halfway != v_rg.end(); ++ui, ++halfway) {
      EXPECT_TRUE(add_edge(*ui, *halfway, g).second);
    }
  }
};
TYPED_TEST_SUITE_P(MatchingTest);

TYPED_TEST_P(MatchingTest, Edmonds) {
  using Base = MatchingTest<TypeParam>;
  using Graph = Base::Graph;
  constexpr std::size_t num_v = Base::num_v;
  constexpr std::size_t double_num_v = Base::double_num_v;
  Graph g(double_num_v);
  Base::complete_graph(g, double_num_v);

  vertex_index_install(g);
  auto vm = get(vertex_index, g);

  // find a maximum cardinality matching using edmonds' blossom-shrinking
  // algorithm, starting with an empty matching.
  auto edmonds_mate = vector_property_map(double_num_v, vm, graph_traits<Graph>::null_vertex());
  empty_matching(g, edmonds_mate.ref());
  matching_loop(edmonds_mate.ref(), edmonds_augmenting_path_finder(g, edmonds_mate.ref(), vm));
  EXPECT_TRUE(maximum_cardinality_matching_verify(g, edmonds_mate.ref(), vm));
  EXPECT_TRUE(is_a_matching(g, edmonds_mate.ref()));
  EXPECT_EQ(matching_size(g, edmonds_mate.ref()), num_v);

  // Now remove an edge from the edmonds_mate matching.
  for (auto v : vertices(g)) {
    if (edmonds_mate[v] != graph_traits<Graph>::null_vertex()) {
      edmonds_mate[edmonds_mate[v]] = graph_traits<Graph>::null_vertex();
      edmonds_mate[v] = graph_traits<Graph>::null_vertex();
      break;
    }
  }

  //...and run the matching verifier - it should tell us that the matching
  // isn't a maximum matching.
  EXPECT_FALSE(maximum_cardinality_matching_verify(g, edmonds_mate.ref(), vm));
}

TYPED_TEST_P(MatchingTest, Greedy) {
  using Base = MatchingTest<TypeParam>;
  using Graph = Base::Graph;
  constexpr std::size_t num_v = Base::num_v;
  constexpr std::size_t double_num_v = Base::double_num_v;
  Graph g(double_num_v);
  Base::complete_graph(g, double_num_v);

  vertex_index_install(g);
  auto vm = get(vertex_index, g);

  // find a greedy matching
  auto greedy_mate = vector_property_map(double_num_v, vm, graph_traits<Graph>::null_vertex());
  greedy_matching(g, greedy_mate.ref());
  EXPECT_TRUE(maximum_cardinality_matching_verify(g, greedy_mate.ref(), vm));
  EXPECT_TRUE(is_a_matching(g, greedy_mate.ref()));
  EXPECT_EQ(matching_size(g, greedy_mate.ref()), num_v);
}

TYPED_TEST_P(MatchingTest, ExtraGreedy) {
  using Base = MatchingTest<TypeParam>;
  using Graph = Base::Graph;
  constexpr std::size_t num_v = Base::num_v;
  constexpr std::size_t double_num_v = Base::double_num_v;
  Graph g(double_num_v);
  Base::complete_graph(g, double_num_v);

  vertex_index_install(g);
  auto vm = get(vertex_index, g);

  // find an extra greedy matching
  auto extra_greedy_mate = vector_property_map(double_num_v, vm, graph_traits<Graph>::null_vertex());
  extra_greedy_matching(g, extra_greedy_mate.ref());
  EXPECT_TRUE(maximum_cardinality_matching_verify(g, extra_greedy_mate.ref(), vm));
  EXPECT_TRUE(is_a_matching(g, extra_greedy_mate.ref()));
  EXPECT_EQ(matching_size(g, extra_greedy_mate.ref()), num_v);
}

TYPED_TEST_P(MatchingTest, Gabow) {
  using Base = MatchingTest<TypeParam>;
  using Graph = Base::Graph;
  constexpr std::size_t num_v = Base::num_v;
  constexpr std::size_t double_num_v = Base::double_num_v;
  Graph g(double_num_v);
  Base::gabows_graph(g, num_v);

  vertex_index_install(g);
  auto vm = get(vertex_index, g);

  // gabow's graph always has a perfect matching. it's also a good example of
  // why one should initialize with the extra_greedy_matching in most cases.

  auto gabow_mate = vector_property_map(double_num_v, vm, graph_traits<Graph>::null_vertex());
  EXPECT_TRUE(checked_edmonds_maximum_cardinality_matching(g, gabow_mate.ref()));
  EXPECT_EQ(matching_size(g, gabow_mate.ref()), num_v);
}

TYPED_TEST_P(MatchingTest, Random) {
  using Base = MatchingTest<TypeParam>;
  using Graph = Base::Graph;
  constexpr std::size_t double_num_v = Base::double_num_v;
  Graph g(double_num_v);

  vertex_index_install(g);
  auto vm = get(vertex_index, g);

  std::mt19937 gen(42);
  while (num_edges(g) < 4 * double_num_v) {
    auto u = random_vertex(g, gen);
    auto v = random_vertex(g, gen);
    while (u == v) {
      v = random_vertex(g, gen);
    }
    add_edge(u, v, g);
  }

  auto random_mate = vector_property_map(double_num_v, vm, graph_traits<Graph>::null_vertex());
  EXPECT_TRUE(checked_edmonds_maximum_cardinality_matching(g, random_mate.ref()));

  // Now remove an edge from the random_mate matching.
  for (auto v : vertices(g)) {
    if (random_mate[v] != graph_traits<Graph>::null_vertex()) {
      random_mate[random_mate[v]] = graph_traits<Graph>::null_vertex();
      random_mate[v] = graph_traits<Graph>::null_vertex();
      break;
    }
  }

  //...and run the matching verifier - it should tell us that the matching
  // isn't a maximum matching.
  EXPECT_FALSE(maximum_cardinality_matching_verify(g, random_mate.ref(), vm));
}

REGISTER_TYPED_TEST_SUITE_P(MatchingTest, Edmonds, Greedy, ExtraGreedy, Gabow, Random);

using alltests = ::testing::Types<std::tuple<VecGraph, std::integral_constant<std::size_t, 10>>,
                                  std::tuple<ListGraph, std::integral_constant<std::size_t, 10>>,
                                  std::tuple<MatrixGraph, std::integral_constant<std::size_t, 10>>,
                                  std::tuple<VecGraph, std::integral_constant<std::size_t, 20>>,
                                  std::tuple<ListGraph, std::integral_constant<std::size_t, 20>>,
                                  std::tuple<MatrixGraph, std::integral_constant<std::size_t, 20>>,
                                  std::tuple<VecGraph, std::integral_constant<std::size_t, 21>>,
                                  std::tuple<ListGraph, std::integral_constant<std::size_t, 21>>,
                                  std::tuple<MatrixGraph, std::integral_constant<std::size_t, 21>>>;

INSTANTIATE_TYPED_TEST_SUITE_P(AllMatchingTest, MatchingTest, alltests);

}  // namespace
}  // namespace bagl
