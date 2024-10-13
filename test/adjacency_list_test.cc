// (C) Copyright Andrew Sutton 2009
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/adjacency_list.h"

#include <iostream>

#include "bagl/graph_selectors.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

template <typename T>
class AdjacencyListUniqueEdgeTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(AdjacencyListUniqueEdgeTest);

template <typename T>
class AdjacencyListTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(AdjacencyListTest);

TYPED_TEST_P(AdjacencyListUniqueEdgeTest, NonLoop) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  // Build a graph with 1 edge and turn it into a loop.
  Graph g(5);
  Vertex u = *vertices(g).begin();
  Vertex v = *std::next(vertices(g).begin(), 2);
  add_edge(u, v, g);
  EXPECT_EQ(num_vertices(g), 5);
  EXPECT_EQ(num_edges(g), 1);
  remove_edge(u, v, g);
  EXPECT_EQ(num_edges(g), 0);
}

TYPED_TEST_P(AdjacencyListTest, NonLoop) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  // Build a graph with 1 edge and turn it into a loop.
  Graph g(5);
  Vertex u = *vertices(g).begin();
  Vertex v = *std::next(vertices(g).begin(), 2);
  add_edge(u, v, g);
  add_edge(u, v, g);
  EXPECT_EQ(num_vertices(g), 5);
  EXPECT_EQ(num_edges(g), 2);
  remove_edge(u, v, g);
  EXPECT_EQ(num_edges(g), 0);
}

TYPED_TEST_P(AdjacencyListUniqueEdgeTest, Loop) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g(5);
  Vertex v = *std::next(vertices(g).begin(), 2);
  add_edge(v, v, g);
  EXPECT_EQ(num_vertices(g), 5);
  EXPECT_EQ(num_edges(g), 1);
  remove_edge(v, v, g);
  EXPECT_EQ(num_edges(g), 0);
}

TYPED_TEST_P(AdjacencyListTest, Loop) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g(5);
  Vertex v = *std::next(vertices(g).begin(), 2);
  add_edge(v, v, g);
  add_edge(v, v, g);
  EXPECT_EQ(num_vertices(g), 5);
  EXPECT_EQ(num_edges(g), 2);
  remove_edge(v, v, g);
  EXPECT_EQ(num_edges(g), 0);
}

using na = no_property;

template <typename Kind>
struct test_graphs {
  using VVL = adjacency_list<vec_s, vec_s, Kind, na, na, na>;
  using PVL = adjacency_list<pool_s, vec_s, Kind, na, na, na>;
  using LVL = adjacency_list<list_s, vec_s, Kind, na, na, na>;
  using MVL = adjacency_list<multiset_s, vec_s, Kind, na, na, na>;
  using UMVL = adjacency_list<unordered_multiset_s, vec_s, Kind, na, na, na>;

  using VLL = adjacency_list<vec_s, list_s, Kind, na, na, na>;
  using PLL = adjacency_list<pool_s, list_s, Kind, na, na, na>;
  using LLL = adjacency_list<list_s, list_s, Kind, na, na, na>;
  using MLL = adjacency_list<multiset_s, list_s, Kind, na, na, na>;
  using UMLL = adjacency_list<unordered_multiset_s, list_s, Kind, na, na, na>;

  using VPL = adjacency_list<vec_s, pool_s, Kind, na, na, na>;
  using PPL = adjacency_list<pool_s, pool_s, Kind, na, na, na>;
  using LPL = adjacency_list<list_s, pool_s, Kind, na, na, na>;
  using MPL = adjacency_list<multiset_s, pool_s, Kind, na, na, na>;
  using UMPL = adjacency_list<unordered_multiset_s, pool_s, Kind, na, na, na>;

  using SVL = adjacency_list<set_s, vec_s, Kind, na, na, na>;
  using USVL = adjacency_list<unordered_set_s, vec_s, Kind, na, na, na>;

  using SLL = adjacency_list<set_s, list_s, Kind, na, na, na>;
  using USLL = adjacency_list<unordered_set_s, list_s, Kind, na, na, na>;

  using SPL = adjacency_list<set_s, pool_s, Kind, na, na, na>;
  using USPL = adjacency_list<unordered_set_s, pool_s, Kind, na, na, na>;

  using multigraphs = ::testing::Types<VVL, PVL, LVL, MVL, UMVL, VLL, PLL, LLL, MLL, UMLL, VPL, PPL, LPL, MPL, UMPL>;
  using uniqueedges = ::testing::Types<SVL, USVL, SLL, USLL, SPL, USPL>;
};

REGISTER_TYPED_TEST_SUITE_P(AdjacencyListUniqueEdgeTest, NonLoop, Loop);
REGISTER_TYPED_TEST_SUITE_P(AdjacencyListTest, NonLoop, Loop);

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListUndirUnique, AdjacencyListUniqueEdgeTest,
                               test_graphs<undirected_s>::uniqueedges);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListUndirMulti, AdjacencyListUniqueEdgeTest,
                               test_graphs<undirected_s>::multigraphs);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListUndirMulti, AdjacencyListTest, test_graphs<undirected_s>::multigraphs);

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListDirUnique, AdjacencyListUniqueEdgeTest,
                               test_graphs<directed_s>::uniqueedges);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListDirMulti, AdjacencyListUniqueEdgeTest,
                               test_graphs<directed_s>::multigraphs);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListDirMulti, AdjacencyListTest, test_graphs<directed_s>::multigraphs);

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListBidirUnique, AdjacencyListUniqueEdgeTest,
                               test_graphs<bidirectional_s>::uniqueedges);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListBidirMulti, AdjacencyListUniqueEdgeTest,
                               test_graphs<bidirectional_s>::multigraphs);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListBidirMulti, AdjacencyListTest, test_graphs<bidirectional_s>::multigraphs);

}  // namespace
}  // namespace bagl
