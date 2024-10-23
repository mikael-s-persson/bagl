// (C) Copyright David Gleich 2007
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/core_numbers.h"

#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(CoreNumbersTest, Test1) {
  // core numbers of sample graph
  using Graph = adjacency_list<vec_s, vec_s, undirected_s>;

  Graph g(21);
  add_edge(0, 1, g);
  add_edge(1, 2, g);
  add_edge(1, 3, g);
  add_edge(2, 3, g);
  add_edge(1, 4, g);
  add_edge(3, 4, g);
  add_edge(4, 5, g);
  add_edge(4, 6, g);
  add_edge(5, 6, g);
  add_edge(4, 7, g);
  add_edge(5, 7, g);
  add_edge(6, 7, g);
  add_edge(7, 8, g);
  add_edge(3, 9, g);
  add_edge(8, 9, g);
  add_edge(8, 10, g);
  add_edge(9, 10, g);
  add_edge(10, 11, g);
  add_edge(10, 12, g);
  add_edge(3, 13, g);
  add_edge(9, 13, g);
  add_edge(3, 14, g);
  add_edge(9, 14, g);
  add_edge(13, 14, g);
  add_edge(16, 17, g);
  add_edge(16, 18, g);
  add_edge(17, 19, g);
  add_edge(18, 19, g);
  add_edge(19, 20, g);

  std::vector<int> core_nums(num_vertices(g));
  core_numbers(g, iterator_property_map(core_nums.begin(), get(vertex_index, g)));

  std::vector<int> correct = {1, 2, 2, 3, 3, 3, 3, 3, 2, 3, 2, 1, 1, 3, 3, 0, 2, 2, 2, 2, 1};
  EXPECT_THAT(core_nums, ::testing::ElementsAreArray(correct));
}

TEST(CoreNumbersTest, Test2) {
  // core numbers of sample graph
  using Graph = adjacency_list<list_s, vec_s, undirected_s, no_property, property<edge_weight_t, int> >;
  std::size_t num_nodes = 3;
  using Edge = std::pair<std::size_t, std::size_t>;

  std::array<Edge, 3> edge_array = {Edge(0, 1), Edge(0, 2), Edge(1, 2)};
  std::array<int, 3> weights = {-1, -2, -2};

  Graph g(num_nodes, std::ranges::ref_view(edge_array), std::ranges::ref_view(weights));

  std::vector<int> core_nums(num_vertices(g));
  weighted_core_numbers(g, iterator_property_map(core_nums.begin(), get(vertex_index, g)));

  EXPECT_THAT(core_nums, ::testing::ElementsAreArray({-1, -1, -4}));
}

TEST(CoreNumbersTest, Test3) {
  // core numbers of a directed graph, the core numbers of a directed
  // cycle are always one
  using Graph = adjacency_list<vec_s, vec_s, directed_s>;
  std::size_t num_nodes = 5;
  using Edge = std::pair<std::size_t, std::size_t>;

  std::array<Edge, 5> edge_array = {Edge(0, 1), Edge(1, 2), Edge(2, 3), Edge(3, 4), Edge(4, 0)};

  Graph g(num_nodes, std::ranges::ref_view(edge_array));

  std::vector<int> core_nums(num_vertices(g));
  core_numbers(g, iterator_property_map(core_nums.begin(), get(vertex_index, g)));

  EXPECT_THAT(core_nums, ::testing::ElementsAreArray({1, 1, 1, 1, 1}));
}

}  // namespace
}  // namespace bagl
