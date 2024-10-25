// Copyright (C) 2005 Jong Soo Park <jongsoo.park -at- gmail.com>
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/dominator_tree.h"

#include <algorithm>
#include <array>
#include <iostream>
#include <utility>
#include <vector>

#include "bagl/adjacency_list.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

struct DominatorCorrectnessTestSet {
  using edge = std::pair<std::size_t, std::size_t>;

  std::size_t num_of_vertices = 0;
  std::vector<edge> edges;
  std::vector<int> correct_idoms;
};

using Graph = adjacency_list<list_s, list_s, bidirectional_s, property<vertex_index_t, std::size_t>, no_property>;

TEST(DominatorTreeTest, TestSets) {
  using edge = DominatorCorrectnessTestSet::edge;

  std::array<DominatorCorrectnessTestSet, 7> test_set{};

  // Tarjan's paper
  test_set[0].num_of_vertices = 13;
  test_set[0].edges =
      std::vector<edge>({edge(0, 1), edge(0, 2),  edge(0, 3),  edge(1, 4),  edge(2, 1),  edge(2, 4),  edge(2, 5),
                         edge(3, 6), edge(3, 7),  edge(4, 12), edge(5, 8),  edge(6, 9),  edge(7, 9),  edge(7, 10),
                         edge(8, 5), edge(8, 11), edge(9, 11), edge(10, 9), edge(11, 0), edge(11, 9), edge(12, 8)});
  test_set[0].correct_idoms = std::vector<int>({std::numeric_limits<int>::max(), 0, 0, 0, 0, 0, 3, 3, 0, 0, 7, 0, 4});

  // Appel. p441. figure 19.4
  test_set[1].num_of_vertices = 7;
  test_set[1].edges = std::vector<edge>(
      {edge(0, 1), edge(1, 2), edge(1, 3), edge(2, 4), edge(2, 5), edge(4, 6), edge(5, 6), edge(6, 1)});
  test_set[1].correct_idoms = std::vector<int>({std::numeric_limits<int>::max(), 0, 1, 1, 2, 2, 2});

  // Appel. p449. figure 19.8
  test_set[2].num_of_vertices = 13;
  test_set[2].edges =
      std::vector<edge>({edge(0, 1), edge(0, 2), edge(1, 3), edge(1, 6), edge(2, 4), edge(2, 7), edge(3, 5), edge(3, 6),
                         edge(4, 7), edge(4, 2), edge(5, 8), edge(5, 10), edge(6, 9), edge(7, 12), edge(8, 11),
                         edge(9, 8), edge(10, 11), edge(11, 1), edge(11, 12)});
  test_set[2].correct_idoms = std::vector<int>({std::numeric_limits<int>::max(), 0, 0, 1, 2, 3, 1, 2, 1, 6, 5, 1, 0});

  test_set[3].num_of_vertices = 8;
  test_set[3].edges = std::vector<edge>(
      {edge(0, 1), edge(1, 2), edge(1, 3), edge(2, 7), edge(3, 4), edge(4, 5), edge(4, 6), edge(5, 7), edge(6, 4)});
  test_set[3].correct_idoms = std::vector<int>({std::numeric_limits<int>::max(), 0, 1, 1, 3, 4, 4, 1});

  // Muchnick. p256. figure 8.21
  test_set[4].num_of_vertices = 8;
  test_set[4].edges = std::vector<edge>(
      {edge(0, 1), edge(1, 2), edge(2, 3), edge(2, 4), edge(3, 2), edge(4, 5), edge(4, 6), edge(5, 7), edge(6, 7)});
  test_set[4].correct_idoms = std::vector<int>({std::numeric_limits<int>::max(), 0, 1, 2, 2, 4, 4, 4});

  // Muchnick. p253. figure 8.18
  test_set[5].num_of_vertices = 8;
  test_set[5].edges = std::vector<edge>(
      {edge(0, 1), edge(0, 2), edge(1, 6), edge(2, 3), edge(2, 4), edge(3, 7), edge(5, 7), edge(6, 7)});
  test_set[5].correct_idoms =
      std::vector<int>({std::numeric_limits<int>::max(), 0, 0, 2, 2, std::numeric_limits<int>::max(), 1, 0});

  // Cytron's paper, fig. 9
  test_set[6].num_of_vertices = 14;
  test_set[6].edges =
      std::vector<edge>({edge(0, 1), edge(0, 13), edge(1, 2), edge(2, 3), edge(2, 7), edge(3, 4), edge(3, 5),
                         edge(4, 6), edge(5, 6), edge(6, 8), edge(7, 8), edge(8, 9), edge(9, 10), edge(9, 11),
                         edge(10, 11), edge(11, 9), edge(11, 12), edge(12, 2), edge(12, 13)});
  test_set[6].correct_idoms =
      std::vector<int>({std::numeric_limits<int>::max(), 0, 1, 2, 3, 3, 3, 2, 2, 8, 9, 9, 11, 0});

  for (const auto& test_case : test_set) {
    Graph g(test_case.num_of_vertices, test_case.edges);

    auto vindex = get(vertex_index, g);
    std::size_t j = 0;
    for (auto u : vertices(g)) {
      put(vindex, u, j++);
    }

    // Lengauer-Tarjan dominator tree algorithm
    auto dom_tree_pred_map = vector_property_map(num_vertices(g), vindex, graph_traits<Graph>::null_vertex());

    lengauer_tarjan_dominator_tree(g, *vertices(g).begin(), dom_tree_pred_map.ref());

    std::vector<int> idom(num_vertices(g));
    for (auto u : vertices(g)) {
      if (get(dom_tree_pred_map, u) != graph_traits<Graph>::null_vertex()) {
        idom[get(vindex, u)] = get(vindex, get(dom_tree_pred_map, u));
      } else {
        idom[get(vindex, u)] = std::numeric_limits<int>::max();
      }
    }

    // dominator tree correctness test
    EXPECT_THAT(idom, ::testing::ElementsAreArray(test_case.correct_idoms));

    // compare results of fast version and slow version of dominator tree
    auto dom_tree_pred_map2 = vector_property_map(num_vertices(g), vindex, graph_traits<Graph>::null_vertex());

    iterative_bit_vector_dominator_tree(g, *vertices(g).begin(), dom_tree_pred_map2.ref());

    std::vector<int> idom2(num_vertices(g));
    for (auto u : vertices(g)) {
      if (get(dom_tree_pred_map2, u) != graph_traits<Graph>::null_vertex()) {
        idom2[get(vindex, u)] = get(vindex, get(dom_tree_pred_map2, u));
      } else {
        idom2[get(vindex, u)] = std::numeric_limits<int>::max();
      }
    }

    // dominator tree correctness test
    EXPECT_THAT(idom2, ::testing::ElementsAreArray(test_case.correct_idoms));

    for (auto u : vertices(g)) {
      EXPECT_EQ(dom_tree_pred_map[u], dom_tree_pred_map2[u]) << "Vertex: " << get(vindex, u);
    }
  }
}

}  // namespace
}  // namespace bagl
