// (C) Copyright Andrew Sutton 2009
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/adjacency_list.h"

#include <cstddef>
#include <iostream>
#include <random>
#include <string_view>

#include "bagl/graph_selectors.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/properties.h"
#include "bagl/random.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

// Create special test id properties.
struct vertex_test_id_t {
  using kind = vertex_property_tag;
  static constexpr std::string_view name = "vertex_test_id";
};
constexpr vertex_test_id_t vertex_test_id = {};
struct edge_test_id_t {
  using kind = edge_property_tag;
  static constexpr std::string_view name = "edge_test_id";
};
constexpr edge_test_id_t edge_test_id = {};

template <typename Graph>
void check_vertex_cleared(const Graph& g, graph_vertex_descriptor_t<Graph> v) {
  const std::size_t v_id = get(vertex_test_id, g, v);
  for (auto u : vertices(g)) {
    auto a_rg = adjacent_vertices(u, g);
    auto a_it = std::find_if(a_rg.begin(), a_rg.end(), [&](auto a) { return get(vertex_test_id, g, a) == v_id; });
    EXPECT_EQ(a_it, a_rg.end()) << "Should not have found vertex " << v_id << ". But found it adjacent to "
                                << get(vertex_test_id, g, u);
  }
}

template <class Graph>
std::size_t count_edges(const Graph& g) {
  return std::ranges::distance(edges(g));
}

template <typename Graph>
void check_edge_added(Graph& g, graph_edge_descriptor_t<Graph> e, graph_vertex_descriptor_t<Graph> a,
                      graph_vertex_descriptor_t<Graph> b, std::size_t correct_id, bool inserted) {
  EXPECT_EQ(source(e, g), a) << "Expected vertex " << get(vertex_test_id, g, a) << " as source of "
                             << get(edge_test_id, g, e) << ", but got " << get(vertex_test_id, g, source(e, g));
  EXPECT_EQ(target(e, g), b) << "Expected vertex " << get(vertex_test_id, g, b) << " as target of "
                             << get(edge_test_id, g, e) << ", but got " << get(vertex_test_id, g, target(e, g));
  EXPECT_TRUE(is_adjacent(g, a, b));
  EXPECT_TRUE(in_edge_set(g, e));
  if (inserted) {
    EXPECT_EQ(get(edge_test_id, g, e), correct_id);
  } else {
    EXPECT_EQ(get(edge_test_id, g, e), get(edge_test_id, g, edge(a, b, g).first));
  }
  EXPECT_EQ(num_edges(g), count_edges(g));
}

template <typename T>
class AdjacencyListUniqueEdgeTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(AdjacencyListUniqueEdgeTest);

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

template <typename T>
class AdjacencyListTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(AdjacencyListTest);

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

TYPED_TEST_P(AdjacencyListTest, RemoveEdge) {
  using Graph = TypeParam;
  Graph g(2);
  auto u = *vertices(g).begin();
  auto v = *std::next(vertices(g).begin(), 1);
  auto [e, b] = add_edge(u, v, 42, g);
  EXPECT_TRUE(b);
  EXPECT_EQ(num_edges(g), 1);
  EXPECT_EQ(get(edge_test_id, g, e), 42);
  remove_edge(e, g);
  EXPECT_EQ(num_edges(g), 0);
}

TYPED_TEST_P(AdjacencyListTest, RemoveMiddleVertex) {
  using Graph = TypeParam;
  Graph g(3);
  auto u = *vertices(g).begin();
  auto m = *std::next(vertices(g).begin(), 1);
  auto v = *std::next(vertices(g).begin(), 2);
  auto [e1, e1_exists] = add_edge(v, u, g);
  EXPECT_TRUE(e1_exists);
  put(edge_test_id, g, e1, 1234);
  remove_vertex(m, g);
  EXPECT_EQ(num_vertices(g), 2);
  EXPECT_EQ(num_edges(g), 1);

  // Find e1 again, everything might have been invalidated.
  u = *vertices(g).begin();
  v = *std::next(vertices(g).begin(), 1);
  std::tie(e1, e1_exists) = edge(v, u, g);
  EXPECT_EQ(get(edge_test_id, g, e1), 1234);
}

TYPED_TEST_P(AdjacencyListTest, RemoveEdgeTwoEdges) {
  using Graph = TypeParam;
  Graph g(2);
  auto u = *vertices(g).begin();
  auto v = *std::next(vertices(g).begin(), 1);

  add_edge(u, v, 42, g);
  add_edge(v, u, 17, g);

  auto [e1, e1_found] = edge(u, v, g);
  EXPECT_TRUE(e1_found);
  auto [e2, e2_found] = edge(v, u, g);
  EXPECT_TRUE(e2_found);
  EXPECT_EQ(num_edges(g), 2);
  if constexpr (is_undirected_graph_v<Graph>) {
    EXPECT_THAT(get(edge_test_id, g, e1), ::testing::AnyOf(42, 17));
    EXPECT_THAT(get(edge_test_id, g, e2), ::testing::AnyOf(42, 17));
  } else {
    EXPECT_EQ(get(edge_test_id, g, e1), 42);
    EXPECT_EQ(get(edge_test_id, g, e2), 17);
  }
  remove_edge(e1, g);
  EXPECT_EQ(num_edges(g), 1);

  // e2 has been invalidated, so grab it again
  std::tie(e2, e2_found) = edge(v, u, g);
  EXPECT_TRUE(e2_found);
  if constexpr (is_undirected_graph_v<Graph>) {
    EXPECT_THAT(get(edge_test_id, g, e2), ::testing::AnyOf(42, 17));
  } else {
    EXPECT_EQ(get(edge_test_id, g, e2), 17);
  }

  remove_edge(e2, g);
}

TYPED_TEST_P(AdjacencyListTest, RandomMutations) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;

  std::size_t N = 5;

  Graph g;

  std::size_t current_vertex_id = 0;
  std::size_t current_edge_id = 0;

  property_map_t<Graph, vertex_test_id_t> vertex_id_map = get(vertex_test_id, g);
  property_map_t<Graph, edge_test_id_t> edge_id_map = get(edge_test_id, g);

  for (std::size_t k = 0; k < N; ++k) {
    add_vertex(current_vertex_id++, g);
  }

  // Generate random values from a fixed seed.
  std::mt19937 gen{42};

  for (std::size_t j = 0; j < 10; ++j) {
    // add_edge
    std::cout << "Adding edges..." << std::endl;
    for (std::size_t i = 0; i < 6; ++i) {
      Vertex b = graph_traits<Graph>::null_vertex();
      Vertex a = random_vertex(g, gen);
      do {
        b = random_vertex(g, gen);
      } while (a == b);  // don't do self edges

      std::cout << "add_edge(" << vertex_id_map[a] << "," << vertex_id_map[b] << ")" << std::endl;
      auto [e, inserted] = add_edge(a, b, current_edge_id++, g);
      check_edge_added(g, e, a, b, current_edge_id - 1, inserted);
    }
    std::cout << "After adding edges:" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>\nEdges:" << std::endl;
    print_edges2(g, vertex_id_map, edge_id_map);
    std::cout << "==========================\nGraph:" << std::endl;
    print_graph(g, vertex_id_map);
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

    // remove_edge(u, v, g)
    std::cout << "Removing (u,v) edges..." << std::endl;
    for (std::size_t i = 0; i < 2; ++i) {
      Edge e = random_edge(g, gen);
      auto [a, b] = incident(e, g);
      std::cout << "remove_edge(" << vertex_id_map[a] << "," << vertex_id_map[b] << ")" << std::endl;
      remove_edge(a, b, g);
      EXPECT_FALSE(is_adjacent(g, a, b));
      EXPECT_FALSE(in_edge_set(g, a, b));
      EXPECT_EQ(num_edges(g), count_edges(g));
    }
    std::cout << "After removing (u,v) edges:" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>\nEdges:" << std::endl;
    print_edges(g, vertex_id_map);
    std::cout << "==========================\nGraph:" << std::endl;
    print_graph(g, vertex_id_map);
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

    // remove_edge(e, g)
    std::cout << "Removing (e) edges..." << std::endl;
    for (std::size_t i = 0; i < 2; ++i) {
      Edge e = random_edge(g, gen);
      auto [a, b] = incident(e, g);
      std::cout << "remove_edge(" << vertex_id_map[a] << "," << vertex_id_map[b] << ")" << std::endl;
      std::size_t old_E = num_edges(g);
      remove_edge(e, g);
      EXPECT_EQ(old_E, num_edges(g) + 1);
      EXPECT_EQ(num_edges(g), count_edges(g));
    }
    std::cout << "After removing (e) edges:" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>\nEdges:" << std::endl;
    print_edges(g, vertex_id_map);
    std::cout << "==========================\nVertices:" << std::endl;
    print_vertices(g, vertex_id_map);
    std::cout << "==========================\nGraph:" << std::endl;
    print_graph(g, vertex_id_map);
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

    // add_vertex
    std::cout << "Adding vertices again..." << std::endl;
    std::size_t old_N = num_vertices(g);
    auto vid = add_vertex(g);
    auto vidp1 = add_vertex(g);
    vertex_id_map[vid] = current_vertex_id++;
    vertex_id_map[vidp1] = current_vertex_id++;

    std::cout << "After adding vertices:" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>\nVertices:" << std::endl;
    print_vertices(g, vertex_id_map);
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

    // make sure the two added vertices are in the graph's vertex set
    EXPECT_TRUE(in_vertex_set(g, vid)) << vertex_id_map[vid] << " not in vertices(g)";
    EXPECT_TRUE(in_vertex_set(g, vidp1)) << vertex_id_map[vidp1] << " not in vertices(g)";

    // make sure the vertices do not have any out edges yet
    EXPECT_EQ(out_degree(vid, g), 0) << vertex_id_map[vid] << " has unexpected out edges";
    EXPECT_TRUE(out_edges(vid, g).empty()) << vertex_id_map[vid] << " has unexpected out edges";
    EXPECT_EQ(out_degree(vidp1, g), 0) << vertex_id_map[vidp1] << " has unexpected out edges";
    EXPECT_TRUE(out_edges(vidp1, g).empty()) << vertex_id_map[vidp1] << " has unexpected out edges";

    // make sure the vertices do not yet appear in any of the edges
    for (auto e : edges(g)) {
      EXPECT_NE(source(e, g), vid) << vertex_id_map[vid] << " is the unexpected source of edge " << edge_id_map[e];
      EXPECT_NE(target(e, g), vid) << vertex_id_map[vid] << " is the unexpected target of edge " << edge_id_map[e];
      EXPECT_NE(source(e, g), vidp1) << vertex_id_map[vidp1] << " is the unexpected source of edge " << edge_id_map[e];
      EXPECT_NE(target(e, g), vidp1) << vertex_id_map[vidp1] << " is the unexpected target of edge " << edge_id_map[e];
    }
    // Make sure num_vertices(g) has been updated
    N = num_vertices(g);
    EXPECT_EQ(N, old_N + 2);

    // add_edge again
    std::cout << "Adding edges again..." << std::endl;
    for (std::size_t i = 0; i < 2; ++i) {
      Vertex a = random_vertex(g, gen);
      Vertex b = random_vertex(g, gen);
      while (a == vid) {
        a = random_vertex(g, gen);
      }
      while (b == vidp1) {
        b = random_vertex(g, gen);
      }
      std::cout << "add_edge(" << vertex_id_map[vid] << "," << vertex_id_map[a] << ")" << std::endl;
      auto [e, inserted] = add_edge(vid, a, current_edge_id++, g);
      check_edge_added(g, e, vid, a, current_edge_id - 1, inserted);

      std::cout << "add_edge(" << vertex_id_map[b] << "," << vertex_id_map[vidp1] << ")" << std::endl;
      // add_edge without property
      std::tie(e, inserted) = add_edge(b, vidp1, g);
      if (inserted) {
        edge_id_map[e] = current_edge_id;
      }
      ++current_edge_id;
      check_edge_added(g, e, b, vidp1, current_edge_id - 1, inserted);
    }
    std::cout << "After adding edges:" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>\nEdges:" << std::endl;
    print_edges2(g, vertex_id_map, edge_id_map);
    std::cout << "==========================\nGraph:" << std::endl;
    print_graph(g, vertex_id_map);
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

    // clear_vertex
    Vertex c = random_vertex(g, gen);
    std::cout << "Clearing vertex " << vertex_id_map[c] << std::endl;
    clear_vertex(c, g);
    check_vertex_cleared(g, c);
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>\nEdges:" << std::endl;
    print_edges2(g, vertex_id_map, edge_id_map);
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

    // remove_vertex
    std::cout << "Removing vertex " << vertex_id_map[c] << std::endl;
    std::size_t old_c_id = vertex_id_map[c];
    old_N = num_vertices(g);
    remove_vertex(c, g);
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>\nEdges:" << std::endl;
    print_graph(g, vertex_id_map);
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

    EXPECT_EQ(num_vertices(g), old_N - 1);
    EXPECT_EQ(std::ranges::distance(vertices(g)), old_N - 1);
    for (auto v : vertices(g)) {
      EXPECT_NE(vertex_id_map[v], old_c_id);
    }
  }
}

template <typename Kind>
struct test_graphs {
  using vprop = property<vertex_test_id_t, std::size_t>;
  using eprop = property<edge_test_id_t, std::size_t>;

  using VVL = adjacency_list<vec_s, vec_s, Kind, vprop, eprop>;
  using PVL = adjacency_list<pool_s, vec_s, Kind, vprop, eprop>;
  using LVL = adjacency_list<list_s, vec_s, Kind, vprop, eprop>;
  using MVL = adjacency_list<multiset_s, vec_s, Kind, vprop, eprop>;
  using UMVL = adjacency_list<unordered_multiset_s, vec_s, Kind, vprop, eprop>;

  using VLL = adjacency_list<vec_s, list_s, Kind, vprop, eprop>;
  using PLL = adjacency_list<pool_s, list_s, Kind, vprop, eprop>;
  using LLL = adjacency_list<list_s, list_s, Kind, vprop, eprop>;
  using MLL = adjacency_list<multiset_s, list_s, Kind, vprop, eprop>;
  using UMLL = adjacency_list<unordered_multiset_s, list_s, Kind, vprop, eprop>;

  using VPL = adjacency_list<vec_s, pool_s, Kind, vprop, eprop>;
  using PPL = adjacency_list<pool_s, pool_s, Kind, vprop, eprop>;
  using LPL = adjacency_list<list_s, pool_s, Kind, vprop, eprop>;
  using MPL = adjacency_list<multiset_s, pool_s, Kind, vprop, eprop>;
  using UMPL = adjacency_list<unordered_multiset_s, pool_s, Kind, vprop, eprop>;

  using SVL = adjacency_list<set_s, vec_s, Kind, vprop, eprop>;
  using USVL = adjacency_list<unordered_set_s, vec_s, Kind, vprop, eprop>;

  using SLL = adjacency_list<set_s, list_s, Kind, vprop, eprop>;
  using USLL = adjacency_list<unordered_set_s, list_s, Kind, vprop, eprop>;

  using SPL = adjacency_list<set_s, pool_s, Kind, vprop, eprop>;
  using USPL = adjacency_list<unordered_set_s, pool_s, Kind, vprop, eprop>;

  using multigraphs = ::testing::Types<VVL, PVL, LVL, MVL, UMVL, VLL, PLL, LLL, MLL, UMLL, VPL, PPL, LPL, MPL, UMPL>;
  using uniqueedges = ::testing::Types<SVL, USVL, SLL, USLL, SPL, USPL>;
};

REGISTER_TYPED_TEST_SUITE_P(AdjacencyListUniqueEdgeTest, NonLoop, Loop);
REGISTER_TYPED_TEST_SUITE_P(AdjacencyListTest, NonLoop, Loop, RandomMutations, RemoveEdge, RemoveMiddleVertex,
                            RemoveEdgeTwoEdges);

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
