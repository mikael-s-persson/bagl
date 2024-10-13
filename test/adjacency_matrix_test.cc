// Copyright Cromwell D. Enage 2004
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/adjacency_matrix.h"

#include <algorithm>  // For std::sort
#include <iostream>
#include <type_traits>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/zip_range.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

template <typename T>
class AdjacencyMatrixTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(AdjacencyMatrixTest);

TYPED_TEST_P(AdjacencyMatrixTest, Basic) {
  using Graph1 = typename TypeParam::first_type;
  using Graph2 = typename TypeParam::second_type;
  typedef typename property_map<Graph1, vertex_index_t>::type IndexMap1;
  typedef typename property_map<Graph2, vertex_index_t>::type IndexMap2;

  Graph1 g1(24);
  Graph2 g2(24);

  add_edge(1, 2, g1);
  add_edge(1, 2, g2);
  add_edge(2, 10, g1);
  add_edge(2, 10, g2);
  add_edge(2, 5, g1);
  add_edge(2, 5, g2);
  add_edge(3, 10, g1);
  add_edge(3, 10, g2);
  add_edge(3, 0, g1);
  add_edge(3, 0, g2);
  add_edge(4, 5, g1);
  add_edge(4, 5, g2);
  add_edge(4, 0, g1);
  add_edge(4, 0, g2);
  add_edge(5, 14, g1);
  add_edge(5, 14, g2);
  add_edge(6, 3, g1);
  add_edge(6, 3, g2);
  add_edge(7, 17, g1);
  add_edge(7, 17, g2);
  add_edge(7, 11, g1);
  add_edge(7, 11, g2);
  add_edge(8, 17, g1);
  add_edge(8, 17, g2);
  add_edge(8, 1, g1);
  add_edge(8, 1, g2);
  add_edge(9, 11, g1);
  add_edge(9, 11, g2);
  add_edge(9, 1, g1);
  add_edge(9, 1, g2);
  add_edge(10, 19, g1);
  add_edge(10, 19, g2);
  add_edge(10, 15, g1);
  add_edge(10, 15, g2);
  add_edge(10, 8, g1);
  add_edge(10, 8, g2);
  add_edge(11, 19, g1);
  add_edge(11, 19, g2);
  add_edge(11, 15, g1);
  add_edge(11, 15, g2);
  add_edge(11, 4, g1);
  add_edge(11, 4, g2);
  add_edge(12, 19, g1);
  add_edge(12, 19, g2);
  add_edge(12, 8, g1);
  add_edge(12, 8, g2);
  add_edge(12, 4, g1);
  add_edge(12, 4, g2);
  add_edge(13, 15, g1);
  add_edge(13, 15, g2);
  add_edge(13, 8, g1);
  add_edge(13, 8, g2);
  add_edge(13, 4, g1);
  add_edge(13, 4, g2);
  add_edge(14, 22, g1);
  add_edge(14, 22, g2);
  add_edge(14, 12, g1);
  add_edge(14, 12, g2);
  add_edge(15, 22, g1);
  add_edge(15, 22, g2);
  add_edge(15, 6, g1);
  add_edge(15, 6, g2);
  add_edge(16, 12, g1);
  add_edge(16, 12, g2);
  add_edge(16, 6, g1);
  add_edge(16, 6, g2);
  add_edge(17, 20, g1);
  add_edge(17, 20, g2);
  add_edge(18, 9, g1);
  add_edge(18, 9, g2);
  add_edge(19, 23, g1);
  add_edge(19, 23, g2);
  add_edge(19, 18, g1);
  add_edge(19, 18, g2);
  add_edge(20, 23, g1);
  add_edge(20, 23, g2);
  add_edge(20, 13, g1);
  add_edge(20, 13, g2);
  add_edge(21, 18, g1);
  add_edge(21, 18, g2);
  add_edge(21, 13, g1);
  add_edge(21, 13, g2);
  add_edge(22, 21, g1);
  add_edge(22, 21, g2);
  add_edge(23, 16, g1);
  add_edge(23, 16, g2);

  IndexMap1 index_map1 = get(vertex_index, g1);
  IndexMap2 index_map2 = get(vertex_index, g2);

  for (auto [v1, v2] : zip_range(vertices(g1), vertices(g2))) {
    EXPECT_EQ(get(index_map1, v1), get(index_map2, v2));
    std::vector<std::size_t> v1_neighbors;
    for (auto a1 : adjacent_vertices(v1, g1)) {
      v1_neighbors.emplace_back(get(index_map1, a1));
    }
    std::vector<std::size_t> v2_neighbors;
    for (auto a2 : adjacent_vertices(v2, g2)) {
      v2_neighbors.emplace_back(get(index_map2, a2));
    }
    EXPECT_THAT(v1_neighbors, ::testing::UnorderedElementsAreArray(v2_neighbors));
  }

  for (auto [v1, v2] : zip_range(vertices(g1), vertices(g2))) {
    EXPECT_EQ(get(index_map1, v1), get(index_map2, v2));
    std::vector<std::size_t> v1_neighbors;
    for (auto e1 : out_edges(v1, g1)) {
      v1_neighbors.emplace_back(get(index_map1, target(e1, g1)));
    }
    std::vector<std::size_t> v2_neighbors;
    for (auto e2 : out_edges(v2, g2)) {
      v2_neighbors.emplace_back(get(index_map2, target(e2, g2)));
    }
    EXPECT_THAT(v1_neighbors, ::testing::UnorderedElementsAreArray(v2_neighbors));
  }

  for (auto [v1, v2] : zip_range(vertices(g1), vertices(g2))) {
    EXPECT_EQ(get(index_map1, v1), get(index_map2, v2));
    std::vector<std::size_t> v1_neighbors;
    for (auto e1 : in_edges(v1, g1)) {
      v1_neighbors.emplace_back(get(index_map1, source(e1, g1)));
    }
    std::vector<std::size_t> v2_neighbors;
    for (auto e2 : in_edges(v2, g2)) {
      v2_neighbors.emplace_back(get(index_map2, source(e2, g2)));
    }
    EXPECT_THAT(v1_neighbors, ::testing::UnorderedElementsAreArray(v2_neighbors));
  }

  // Test construction from a range of pairs
  std::vector<std::pair<int, int>> edge_pairs_g1;
  for (auto e : edges(g1)) {
    edge_pairs_g1.push_back(std::make_pair(get(index_map1, source(e, g1)), get(index_map1, target(e, g1))));
  }
  Graph2 g3(num_vertices(g1), edge_pairs_g1);
  EXPECT_EQ(num_vertices(g1), num_vertices(g3));
  std::vector<std::pair<int, int>> edge_pairs_g3;
  IndexMap2 index_map3 = get(vertex_index_t(), g3);
  for (auto e : edges(g3)) {
    edge_pairs_g3.push_back(std::make_pair(get(index_map3, source(e, g3)), get(index_map3, target(e, g3))));
  }
  // Normalize the edge pairs for comparison
  if (is_undirected_graph_v<Graph1> || is_undirected_graph_v<Graph2>) {
    for (size_t i = 0; i < edge_pairs_g1.size(); ++i) {
      if (edge_pairs_g1[i].first < edge_pairs_g1[i].second) {
        std::swap(edge_pairs_g1[i].first, edge_pairs_g1[i].second);
      }
    }
    for (size_t i = 0; i < edge_pairs_g3.size(); ++i) {
      if (edge_pairs_g3[i].first < edge_pairs_g3[i].second) {
        std::swap(edge_pairs_g3[i].first, edge_pairs_g3[i].second);
      }
    }
  }
  std::sort(edge_pairs_g1.begin(), edge_pairs_g1.end());
  std::sort(edge_pairs_g3.begin(), edge_pairs_g3.end());
  edge_pairs_g1.erase(std::unique(edge_pairs_g1.begin(), edge_pairs_g1.end()), edge_pairs_g1.end());
  edge_pairs_g3.erase(std::unique(edge_pairs_g3.begin(), edge_pairs_g3.end()), edge_pairs_g3.end());
  EXPECT_EQ(edge_pairs_g1, edge_pairs_g3);
}

TYPED_TEST_P(AdjacencyMatrixTest, RemoveEdge) {
  using Graph2 = typename TypeParam::second_type;
  // Build a 2-vertex graph
  Graph2 g(2);
  add_edge(0, 1, g);
  EXPECT_EQ(num_vertices(g), 2);
  EXPECT_EQ(num_edges(g), 1);
  remove_edge(0, 1, g);
  EXPECT_EQ(num_edges(g), 0);

  // Make sure we don't decrement the edge count if the edge doesn't actually
  // exist.
  remove_edge(0, 1, g);
  EXPECT_EQ(num_edges(g), 0);
}

REGISTER_TYPED_TEST_SUITE_P(AdjacencyMatrixTest, Basic, RemoveEdge);

using AdjacentyMatrixTypes =
    ::testing::Types<std::pair<adjacency_list<set_s, vec_s, undirected_s>, adjacency_matrix<undirected_s>>,
                     std::pair<adjacency_list<set_s, vec_s, bidirectional_s>, adjacency_matrix<directed_s>>>;
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyMatrix, AdjacencyMatrixTest, AdjacentyMatrixTypes);

}  // namespace
}  // namespace bagl
