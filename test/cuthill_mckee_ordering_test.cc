// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//          Doug Gregor, D. Kevin McGrath
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/cuthill_mckee_ordering.h"

#include <array>
#include <iostream>
#include <utility>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/bandwidth.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

class CuthillMcKeeOrderingTest : public ::testing::Test {
 public:
  using Graph = adjacency_list<vec_s, vec_s, undirected_s, property<vertex_color_t, default_color_type> >;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  CuthillMcKeeOrderingTest() : g(10) {
    using Pair = std::pair<std::size_t, std::size_t>;
    std::array edges = {Pair(0, 3),   // a-d
                        Pair(0, 5),   // a-f
                        Pair(1, 2),   // b-c
                        Pair(1, 4),   // b-e
                        Pair(1, 6),   // b-g
                        Pair(1, 9),   // b-j
                        Pair(2, 3),   // c-d
                        Pair(2, 4),   // c-e
                        Pair(3, 5),   // d-f
                        Pair(3, 8),   // d-i
                        Pair(4, 6),   // e-g
                        Pair(5, 6),   // f-g
                        Pair(5, 7),   // f-h
                        Pair(6, 7)};  // g-h

    for (auto p : edges) {
      add_edge(p.first, p.second, g);
    }

    EXPECT_EQ(bandwidth(g), 8);
  }

  Graph g;
};

TEST_F(CuthillMcKeeOrderingTest, StartingAtSix) {
  auto vindex = get(vertex_index, g);
  std::vector<Vertex> inv_perm(num_vertices(g));
  std::vector<std::size_t> perm(num_vertices(g));
  Vertex s = *std::next(vertices(g).begin(), 6);
  // reverse cuthill_mckee_ordering
  cuthill_mckee_ordering(g, s, inv_perm.rbegin(), get(vertex_color, g), get(vertex_degree, g));
  std::vector<std::size_t> inv_perm_ids(num_vertices(g));
  std::transform(inv_perm.begin(), inv_perm.end(), inv_perm_ids.begin(), [&vindex](auto v) { return vindex[v]; });
  EXPECT_THAT(inv_perm_ids, ::testing::ElementsAreArray({8, 3, 0, 9, 2, 5, 1, 4, 7, 6}));

  for (std::size_t c = 0; c != inv_perm.size(); ++c) {
    perm[vindex[inv_perm[c]]] = c;
  }
  EXPECT_EQ(bandwidth(g, iterator_property_map(&perm[0], vindex)), 4);
}
TEST_F(CuthillMcKeeOrderingTest, StartingAtZero) {
  auto vindex = get(vertex_index, g);
  std::vector<Vertex> inv_perm(num_vertices(g));
  std::vector<std::size_t> perm(num_vertices(g));
  Vertex s = *vertices(g).begin();
  // reverse cuthill_mckee_ordering
  cuthill_mckee_ordering(g, s, inv_perm.rbegin(), get(vertex_color, g), get(vertex_degree, g));
  std::vector<std::size_t> inv_perm_ids(num_vertices(g));
  std::transform(inv_perm.begin(), inv_perm.end(), inv_perm_ids.begin(), [&vindex](auto v) { return vindex[v]; });
  EXPECT_THAT(inv_perm_ids, ::testing::ElementsAreArray({9, 1, 4, 6, 7, 2, 8, 5, 3, 0}));

  for (std::size_t c = 0; c != inv_perm.size(); ++c) {
    perm[vindex[inv_perm[c]]] = c;
  }
  EXPECT_EQ(bandwidth(g, iterator_property_map(&perm[0], vindex)), 4);
}

TEST_F(CuthillMcKeeOrderingTest, Basic) {
  auto vindex = get(vertex_index, g);
  std::vector<Vertex> inv_perm(num_vertices(g));
  std::vector<std::size_t> perm(num_vertices(g));
  // reverse cuthill_mckee_ordering
  cuthill_mckee_ordering(g, inv_perm.rbegin());
  std::vector<std::size_t> inv_perm_ids(num_vertices(g));
  std::transform(inv_perm.begin(), inv_perm.end(), inv_perm_ids.begin(), [&vindex](auto v) { return vindex[v]; });
  EXPECT_THAT(inv_perm_ids, ::testing::ElementsAreArray({0, 8, 5, 7, 3, 6, 4, 2, 1, 9}));

  for (std::size_t c = 0; c != inv_perm.size(); ++c) {
    perm[vindex[inv_perm[c]]] = c;
  }
  EXPECT_EQ(bandwidth(g, iterator_property_map(&perm[0], vindex)), 4);
}

}  // namespace
}  // namespace bagl
