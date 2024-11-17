// Copyright 2024 Mikael Persson

#ifndef BAGL_TEST_GRAPH_MUTATION_TEST_SUITE_H_
#define BAGL_TEST_GRAPH_MUTATION_TEST_SUITE_H_

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <random>
#include <string>
#include <string_view>

#include "bagl/filtered_graph.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_selectors.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/random.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {

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

struct vertex_test_bundle_t {
  std::string str_value;
};

struct edge_test_bundle_t {
  std::string str_value;
};

// Use the following vertex and edge properties in the graphs to be tested.
using test_vertex_property = property<vertex_test_id_t, std::size_t, vertex_test_bundle_t>;
using test_edge_property = property<edge_test_id_t, std::size_t, edge_test_bundle_t>;

template <concepts::IncidenceGraph Graph>
requires concepts::VertexListGraph<Graph>
void check_vertex_cleared(const Graph& g, graph_vertex_descriptor_t<Graph> v) {
  const std::size_t v_id = get(vertex_test_id, g, v);
  for (auto u : vertices(g)) {
    auto oe_rg = out_edges(u, g);
    auto oe_it =
        std::find_if(oe_rg.begin(), oe_rg.end(), [&](auto e) { return get(vertex_test_id, g, target(e, g)) == v_id; });
    EXPECT_EQ(oe_it, oe_rg.end()) << "Should not have found vertex " << v_id << ". But found it adjacent to "
                                  << get(vertex_test_id, g, u);
  }
}

template <concepts::Graph Graph>
std::size_t count_edges(const Graph& g) {
  if constexpr (concepts::EdgeListGraph<Graph>) {
    return std::ranges::distance(edges(g));
  } else if constexpr (concepts::VertexListGraph<Graph> && concepts::IncidenceGraph<Graph>) {
    std::size_t e_count = 0;
    for (auto u : vertices(g)) {
      e_count += std::ranges::distance(out_edges(u, g));
    }
    if constexpr (is_undirected_graph_v<Graph>) {
      e_count /= 2;
    }
    return e_count;
  } else {
    // Probably will fail on checks.
    return 0;
  }
}

template <concepts::Graph Graph>
std::size_t num_edges_or_count(const Graph& g) {
  if constexpr (concepts::EdgeListGraph<Graph>) {
    return num_edges(g);
  } else {
    return count_edges(g);
  }
}

template <concepts::IncidenceGraph Graph>
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
  EXPECT_EQ(num_edges_or_count(g), count_edges(g));
}

template <concepts::VertexMutableGraph Graph>
graph_vertex_descriptor_t<Graph> add_vertex_and_check(Graph& g, std::size_t new_id, std::string new_str = "") {
  if constexpr (concepts::VertexMutablePropertyGraph<Graph>) {
    graph_vertex_descriptor_t<Graph> v = graph_traits<Graph>::null_vertex();
    EXPECT_NO_THROW(v = add_vertex(g, new_id, new_str));
    auto vid = get(vertex_test_id, g);
    EXPECT_EQ(vid[v], new_id);
    EXPECT_EQ(g[v].str_value, new_str);
    return v;
  } else {
    graph_vertex_descriptor_t<Graph> v = graph_traits<Graph>::null_vertex();
    EXPECT_NO_THROW(v = add_vertex(g));
    auto vid = get(vertex_test_id, g);
    put(vid, v, new_id);
    g[v].str_value = new_str;
    EXPECT_EQ(vid[v], new_id);
    EXPECT_EQ(g[v].str_value, new_str);
    return v;
  }
}

template <concepts::EdgeMutableGraph Graph>
std::pair<graph_edge_descriptor_t<Graph>, bool> add_edge_and_check(Graph& g, graph_vertex_descriptor_t<Graph> u,
                                                                   graph_vertex_descriptor_t<Graph> v,
                                                                   std::size_t new_id, std::string new_str = "") {
  if constexpr (concepts::EdgeMutablePropertyGraph<Graph>) {
    graph_edge_descriptor_t<Graph> e;
    bool e_added = false;
    EXPECT_NO_THROW(std::tie(e, e_added) = add_edge(u, v, g, new_id, new_str));
    if (e_added) {
      auto eid = get(edge_test_id, g);
      EXPECT_EQ(eid[e], new_id);
      EXPECT_EQ(g[e].str_value, new_str);
    }
    return {e, e_added};
  } else {
    graph_edge_descriptor_t<Graph> e;
    bool e_added = false;
    EXPECT_NO_THROW(std::tie(e, e_added) = add_edge(u, v, g));
    if (e_added) {
      auto eid = get(edge_test_id, g);
      put(eid, e, new_id);
      g[e].str_value = new_str;
      EXPECT_EQ(eid[e], new_id);
      EXPECT_EQ(g[e].str_value, new_str);
    }
    return {e, e_added};
  }
}

template <concepts::VertexListGraph Graph>
std::vector<std::size_t> get_graph_vertex_ids(const Graph& g, std::size_t expected_count) {
  auto vid = get(vertex_test_id, g);
  std::vector<std::size_t> vp_list;
  for (auto v : vertices(g)) {
    vp_list.push_back(vid[v]);
  }
  return vp_list;
}

template <concepts::EdgeListGraph Graph>
std::vector<std::size_t> get_graph_edge_ids(const Graph& g, std::size_t expected_count) {
  auto eid = get(edge_test_id, g);
  std::vector<std::size_t> ep_list;
  for (auto e : edges(g)) {
    ep_list.push_back(eid[e]);
  }
  return ep_list;
}

template <typename Graph>
bool allow_new_edge(Graph& g, graph_vertex_descriptor_t<Graph> a, graph_vertex_descriptor_t<Graph> b) {
  if (a == b) {
    // Don't add self edges.
    return false;
  }
  if constexpr (!allows_parallel_edges_v<Graph>) {
    // Search in out-edges to avoid requiring adjacency-matrix, we don't care about performance.
    auto oe_rg = out_edges(a, g);
    // Don't add parallel edges.
    return std::find_if(oe_rg.begin(), oe_rg.end(), [b, &g](auto e) { return target(e, g) == b; }) == oe_rg.end();
  } else {
    return true;
  }
}

// NOTE: Avoid using 'auto' in tests, if reasonable, to ensure functions match the expected types.

template <typename T>
class GraphMutationTest : public ::testing::Test {
  // These requirements are pretty much the minimal for a basic mutable graph.
  // A special-purpose test suite is more appropriate if these can't be met.
  static_assert(concepts::IncidenceGraph<T>);
  static_assert(concepts::VertexListGraph<T>);
  static_assert(concepts::VertexMutableGraph<T>);
  static_assert(concepts::EdgeMutableGraph<T>);
  // Further capabilities, such as edge-list, adjacency, parallel edges, bidir/undir, etc.,
  // enhance the test (e.g., edge lists are checked if possible, parallel edges are added if possible, etc.).
};
TYPED_TEST_SUITE_P(GraphMutationTest);

TYPED_TEST_P(GraphMutationTest, SimpleEdge) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  // Build a graph with 1 edge, maybe parallel.
  Graph g;

  Vertex u = add_vertex(g);
  add_vertex(g);
  Vertex v = add_vertex(g);
  add_vertex(g);
  add_vertex(g);

  std::size_t expected_edge_num = 0;
  add_edge(u, v, g);
  ++expected_edge_num;
  if constexpr (allows_parallel_edges_v<Graph>) {
    add_edge(u, v, g);
    ++expected_edge_num;
  }
  EXPECT_EQ(num_vertices(g), 5);
  EXPECT_EQ(num_edges_or_count(g), expected_edge_num);
  remove_edge(u, v, g);
  // Expect all parallel edges were removed.
  EXPECT_EQ(num_edges_or_count(g), 0);
}

TYPED_TEST_P(GraphMutationTest, SelfEdge) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;

  add_vertex(g);
  add_vertex(g);
  Vertex v = add_vertex(g);
  add_vertex(g);
  add_vertex(g);

  std::size_t expected_edge_num = 0;
  EXPECT_TRUE(add_edge(v, v, g).second);
  ++expected_edge_num;
  if constexpr (allows_parallel_edges_v<Graph>) {
    EXPECT_TRUE(add_edge(v, v, g).second);
    ++expected_edge_num;
  }
  EXPECT_EQ(num_vertices(g), 5);
  EXPECT_EQ(num_edges_or_count(g), expected_edge_num);
  remove_edge(v, v, g);
  // Expect all parallel edges were removed.
  EXPECT_EQ(num_edges_or_count(g), 0);
}

TYPED_TEST_P(GraphMutationTest, RemoveEdge) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;

  Vertex u = add_vertex(g);
  Vertex v = add_vertex(g);
  EXPECT_EQ(num_vertices(g), 2);

  auto [e, b] = add_edge_and_check(g, u, v, 42);
  EXPECT_TRUE(b);
  EXPECT_EQ(num_edges_or_count(g), 1);
  remove_edge(e, g);
  EXPECT_EQ(num_edges_or_count(g), 0);
}

TYPED_TEST_P(GraphMutationTest, RemoveMiddleVertex) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;

  Vertex u = add_vertex(g);
  Vertex m = add_vertex(g);
  Vertex v = add_vertex(g);
  EXPECT_EQ(num_vertices(g), 3);

  auto [e1, e1_exists] = add_edge_and_check(g, v, u, 1234);
  EXPECT_TRUE(e1_exists);
  remove_vertex(m, g);
  EXPECT_EQ(num_vertices(g), 2);
  EXPECT_EQ(num_edges_or_count(g), 1);

  // Find e1 again, vertices might have been invalidated.
  u = *vertices(g).begin();
  v = *std::next(vertices(g).begin(), 1);
  std::tie(e1, e1_exists) = edge(v, u, g);
  EXPECT_TRUE(e1_exists);
  EXPECT_EQ(get(edge_test_id, g, e1), 1234);
}

TYPED_TEST_P(GraphMutationTest, RemoveEdgeTwoEdges) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;

  property_map_t<Graph, edge_test_id_t> edge_id_map = get(edge_test_id, g);

  Vertex u = add_vertex(g);
  Vertex v = add_vertex(g);

  auto [e1, e1_added] = add_edge_and_check(g, u, v, 42);
  EXPECT_TRUE(e1_added);
  auto [e2, e2_added] = add_edge_and_check(g, v, u, 17);
  EXPECT_TRUE(e2_added);
  EXPECT_EQ(num_edges_or_count(g), 2);

  if constexpr (concepts::AdjacencyMatrix<Graph>) {
    auto [e1_f, e1_found] = edge(u, v, g);
    EXPECT_TRUE(e1_found);
    auto [e2_f, e2_found] = edge(v, u, g);
    EXPECT_TRUE(e2_found);
    if constexpr (!is_undirected_graph_v<Graph>) {
      EXPECT_EQ(e1_f, e1);
      EXPECT_EQ(e2_f, e2);
    }
  }
  if constexpr (is_undirected_graph_v<Graph>) {
    std::array eprops{edge_id_map[e1], edge_id_map[e2]};
    EXPECT_THAT(eprops, ::testing::UnorderedElementsAre(42, 17));
  } else {
    std::array eprops{edge_id_map[e1], edge_id_map[e2]};
    EXPECT_THAT(eprops, ::testing::ElementsAre(42, 17));
  }
  remove_edge(e1, g);
  EXPECT_EQ(num_edges_or_count(g), 1);

  // e2 might be invalidated, so find it again
  std::tie(e2, e2_added) = edge(v, u, g);
  EXPECT_TRUE(e2_added);
  if constexpr (is_undirected_graph_v<Graph>) {
    EXPECT_THAT(edge_id_map[e2], ::testing::AnyOf(42, 17));
  } else {
    EXPECT_EQ(edge_id_map[e2], 17);
  }

  remove_edge(e2, g);
  EXPECT_EQ(num_edges_or_count(g), 0);
}

// This is the real big test.
// It makes multiple passes randomly adding/removing vertices and edges,
// and keeps checking the integrity of the graph.
TYPED_TEST_P(GraphMutationTest, RandomMutations) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;

  // N needs to be reasonable large to be able to add enough non-parallel edges.
  // But, not too big, so we get good degrees.
  std::size_t N = 10;

  Graph g;

  std::size_t current_vertex_id = 0;
  std::size_t current_edge_id = 0;

  property_map_t<Graph, vertex_test_id_t> vertex_id_map = get(vertex_test_id, g);
  property_map_t<Graph, edge_test_id_t> edge_id_map = get(edge_test_id, g);

  for (std::size_t k = 0; k < N; ++k) {
    if constexpr (concepts::VertexMutablePropertyGraph<Graph>) {
      add_vertex(g, current_vertex_id++);
    } else {
      auto v = add_vertex(g);
      vertex_id_map[v] = current_vertex_id++;
    }
  }

  // Generate random values from a fixed seed.
  std::mt19937 gen{42};

  for (std::size_t j = 0; j < 10; ++j) {
    // add_edge
    std::cout << "Adding edges..." << std::endl;
    for (std::size_t i = 0; i < 6; ++i) {
      Vertex a = random_vertex(g, gen);
      Vertex b = random_vertex(g, gen);
      while (!allow_new_edge(g, a, b)) {
        b = random_vertex(g, gen);
      }

      std::cout << "add_edge(" << vertex_id_map[a] << "," << vertex_id_map[b] << ")" << std::endl;
      if constexpr (concepts::EdgeMutablePropertyGraph<Graph>) {
        auto [e, inserted] = add_edge(a, b, g, current_edge_id++);
        check_edge_added(g, e, a, b, current_edge_id - 1, inserted);
      } else {
        auto [e, inserted] = add_edge(a, b, g);
        edge_id_map[e] = current_edge_id++;
        check_edge_added(g, e, a, b, current_edge_id - 1, inserted);
      }
    }
    std::cout << "After adding edges:" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    print_graph(g, vertex_id_map, edge_id_map);
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
      EXPECT_EQ(num_edges_or_count(g), count_edges(g));
    }
    std::cout << "After removing (u,v) edges:" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    print_graph(g, vertex_id_map, edge_id_map);
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

    // remove_edge(e, g)
    std::cout << "Removing (e) edges..." << std::endl;
    for (std::size_t i = 0; i < 2; ++i) {
      Edge e = random_edge(g, gen);
      auto [a, b] = incident(e, g);
      std::cout << "remove_edge(" << vertex_id_map[a] << "," << vertex_id_map[b] << ")" << std::endl;
      std::size_t old_E = num_edges_or_count(g);
      if constexpr (concepts::EdgeMutablePropertyGraph<Graph>) {
        std::size_t expected_ep = edge_id_map[e];
        test_edge_property ep;
        remove_edge(e, g, &ep);
        EXPECT_EQ(get_property_value(ep, edge_test_id), expected_ep);
      } else {
        remove_edge(e, g);
      }
      EXPECT_EQ(old_E, num_edges_or_count(g) + 1);
      EXPECT_EQ(num_edges_or_count(g), count_edges(g));
    }
    std::cout << "After removing (e) edges:" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    print_graph(g, vertex_id_map, edge_id_map);
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

    // add_vertex
    std::cout << "Adding vertices again..." << std::endl;
    std::size_t old_N = num_vertices(g);
    auto vid = add_vertex(g);
    auto vidp1 = add_vertex(g);
    vertex_id_map[vid] = current_vertex_id++;
    vertex_id_map[vidp1] = current_vertex_id++;

    std::cout << "After adding vertices:" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    print_graph(g, vertex_id_map, edge_id_map);
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
    if constexpr (concepts::EdgeListGraph<Graph>) {
      for (auto e : edges(g)) {
        EXPECT_NE(source(e, g), vid) << vertex_id_map[vid] << " is the unexpected source of edge " << edge_id_map[e];
        EXPECT_NE(target(e, g), vid) << vertex_id_map[vid] << " is the unexpected target of edge " << edge_id_map[e];
        EXPECT_NE(source(e, g), vidp1) << vertex_id_map[vidp1] << " is the unexpected source of edge "
                                       << edge_id_map[e];
        EXPECT_NE(target(e, g), vidp1) << vertex_id_map[vidp1] << " is the unexpected target of edge "
                                       << edge_id_map[e];
      }
    }
    // Make sure num_vertices(g) has been updated
    N = num_vertices(g);
    EXPECT_EQ(N, old_N + 2);

    // add_edge again
    std::cout << "Adding edges again..." << std::endl;
    for (std::size_t i = 0; i < 2; ++i) {
      Vertex a = random_vertex(g, gen);
      while (!allow_new_edge(g, vid, a)) {
        a = random_vertex(g, gen);
      }
      Vertex b = random_vertex(g, gen);
      while (!allow_new_edge(g, b, vidp1)) {
        b = random_vertex(g, gen);
      }
      std::cout << "add_edge(" << vertex_id_map[vid] << "," << vertex_id_map[a] << ")" << std::endl;
      auto [e, inserted] = add_edge(vid, a, g);
      if (inserted) {
        edge_id_map[e] = current_edge_id++;
      }
      check_edge_added(g, e, vid, a, current_edge_id - 1, inserted);

      std::cout << "add_edge(" << vertex_id_map[b] << "," << vertex_id_map[vidp1] << ")" << std::endl;
      // add_edge without property
      std::tie(e, inserted) = add_edge(b, vidp1, g);
      if (inserted) {
        edge_id_map[e] = current_edge_id++;
      }
      check_edge_added(g, e, b, vidp1, current_edge_id - 1, inserted);
    }
    std::cout << "After adding edges:" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    print_graph(g, vertex_id_map, edge_id_map);
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

    // clear_vertex
    Vertex c = random_vertex(g, gen);
    std::cout << "Clearing vertex " << vertex_id_map[c] << std::endl;
    clear_vertex(c, g);
    check_vertex_cleared(g, c);
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    print_graph(g, vertex_id_map, edge_id_map);
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

    // remove_vertex
    std::cout << "Removing vertex " << vertex_id_map[c] << std::endl;
    std::size_t old_c_id = vertex_id_map[c];
    old_N = num_vertices(g);
    if constexpr (concepts::VertexMutablePropertyGraph<Graph>) {
      std::size_t expected_vp = vertex_id_map[c];
      test_vertex_property vp;
      remove_vertex(c, g, &vp);
      EXPECT_EQ(get_property_value(vp, vertex_test_id), expected_vp);
    } else {
      remove_vertex(c, g);
    }
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    print_graph(g, vertex_id_map, edge_id_map);
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

    EXPECT_EQ(num_vertices(g), old_N - 1);
    EXPECT_EQ(std::ranges::distance(vertices(g)), old_N - 1);
    for (auto v : vertices(g)) {
      EXPECT_NE(vertex_id_map[v], old_c_id);
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(GraphMutationTest, SimpleEdge, SelfEdge, RandomMutations, RemoveEdge, RemoveMiddleVertex,
                            RemoveEdgeTwoEdges);

template <typename T>
class GraphTreeMutationTest : public ::testing::Test {
 public:
  // These requirements are pretty much the minimal for a basic mutable graph.
  // A special-purpose test suite is more appropriate if these can't be met.
  static_assert(concepts::IncidenceGraph<T>);
  static_assert(concepts::VertexListGraph<T>);
  static_assert(concepts::VertexMutableGraph<T>);
  static_assert(concepts::EdgeMutableGraph<T>);
  // Further capabilities, such as edge-list, adjacency, parallel edges, bidir/undir, etc.,
  // enhance the test (e.g., edge lists are checked if possible, parallel edges are added if possible, etc.).

  template <concepts::BidirectionalGraph Graph>
  static std::pair<std::vector<std::size_t>, std::vector<std::size_t>> get_graph_in_edge_ids(
      const Graph& g, graph_vertex_descriptor_t<Graph> v, std::size_t expected_count) {
    EXPECT_EQ(in_degree(v, g), expected_count);
    auto vid = get(vertex_test_id, g);
    auto eid = get(edge_test_id, g);
    std::vector<std::size_t> ep_list;
    std::vector<std::size_t> vp_list;
    for (graph_edge_descriptor_t<Graph> e : in_edges(v, g)) {
      if constexpr (is_undirected_graph_v<Graph>) {
        EXPECT_THAT(eid[e], ::testing::AnyOf(vid[source(e, g)] * 1000 + vid[target(e, g)],
                                             vid[source(e, g)] + vid[target(e, g)] * 1000));
        if (target(e, g) == v) {
          vp_list.push_back(vid[source(e, g)]);
        } else {
          vp_list.push_back(vid[target(e, g)]);
        }
      } else {
        EXPECT_EQ(eid[e], (vid[source(e, g)] * 1000 + vid[target(e, g)]));
        vp_list.push_back(vid[source(e, g)]);
      }
      ep_list.push_back(eid[e]);
    }
    return {ep_list, vp_list};
  }

  template <concepts::IncidenceGraph Graph>
  static std::pair<std::vector<std::size_t>, std::vector<std::size_t>> get_graph_out_edge_ids(
      const Graph& g, graph_vertex_descriptor_t<Graph> u, std::size_t expected_count) {
    EXPECT_EQ(out_degree(u, g), expected_count);
    auto vid = get(vertex_test_id, g);
    auto eid = get(edge_test_id, g);
    std::vector<std::size_t> ep_list;
    std::vector<std::size_t> vp_list;
    for (graph_edge_descriptor_t<Graph> e : out_edges(u, g)) {
      if constexpr (is_undirected_graph_v<Graph>) {
        EXPECT_THAT(eid[e], ::testing::AnyOf(vid[source(e, g)] * 1000 + vid[target(e, g)],
                                             vid[source(e, g)] + vid[target(e, g)] * 1000));
        if (source(e, g) == u) {
          vp_list.push_back(vid[target(e, g)]);
        } else {
          vp_list.push_back(vid[source(e, g)]);
        }
      } else {
        EXPECT_EQ(eid[e], vid[source(e, g)] * 1000 + vid[target(e, g)]);
        vp_list.push_back(vid[target(e, g)]);
      }
      ep_list.push_back(eid[e]);
    }
    return {ep_list, vp_list};
  }

  template <concepts::IncidenceGraph Graph>
  static void check_fullbranch_integrity(const Graph& g, graph_vertex_descriptor_t<Graph> u) {
    auto vid = get(vertex_test_id, g);
    auto eid = get(edge_test_id, g);
    if constexpr (is_undirected_graph_v<Graph>) {
      for (auto v : vertices(g)) {
        for (auto e : out_edges(v, g)) {
          EXPECT_THAT(eid[e], ::testing::AnyOf(vid[source(e, g)] * 1000 + vid[target(e, g)],
                                               vid[source(e, g)] + vid[target(e, g)] * 1000));
        }
      }
    } else {
      if (out_degree(u, g) == 0) {
        return;
      }
      EXPECT_EQ(out_degree(u, g), 4);
      for (auto e : out_edges(u, g)) {
        EXPECT_EQ(eid[e], vid[source(e, g)] * 1000 + vid[target(e, g)]);
        check_fullbranch_integrity(g, target(e, g));
      }
    }
  }
};
TYPED_TEST_SUITE_P(GraphTreeMutationTest);

TYPED_TEST_P(GraphTreeMutationTest, MutateTree) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;

  Graph g;
  Vertex v_root = graph_traits<Graph>::null_vertex();

  property_map_t<Graph, vertex_test_id_t> vertex_id_map = get(vertex_test_id, g);
  property_map_t<Graph, edge_test_id_t> edge_id_map = get(edge_test_id, g);

  /* MutableGraph */
  EXPECT_NO_THROW(v_root = add_vertex(g));
  EXPECT_EQ(num_vertices(g), 1);
  EXPECT_NO_THROW(remove_vertex(v_root, g));
  EXPECT_EQ(num_vertices(g), 0);

  /* MutablePropertyGraph */
  v_root = add_vertex_and_check(g, 10, "10");
  vertex_id_map[v_root] = 1;
  EXPECT_EQ(vertex_id_map[v_root], 1);
  g[v_root].str_value = "1";
  EXPECT_EQ(g[v_root].str_value, "1");

  /* MutableGraph */
  std::array<int, 4> vp_rc = {2, 3, 4, 5};
  std::array<int, 4> ep_rc = {1002, 1003, 1004, 1005};
  std::array<Vertex, 4> v_rc{};
  std::array<Edge, 4> e_rc{};
  for (int i = 0; i < 4; ++i) {
    v_rc[i] = add_vertex_and_check(g, vp_rc[i], std::to_string(vp_rc[i]));
    bool edge_added_success = false;
    EXPECT_NO_THROW(std::tie(e_rc[i], edge_added_success) = add_edge(v_root, v_rc[i], g));
    EXPECT_TRUE(edge_added_success);
    edge_id_map[e_rc[i]] = ep_rc[i];
    check_edge_added(g, e_rc[i], v_root, v_rc[i], ep_rc[i], edge_added_success);
  }
  EXPECT_EQ(num_vertices(g), 5);

  /* MutablePropertyGraph */
  std::array<int, 4> vp_rc1c = {6, 7, 8, 9};
  std::array<int, 4> ep_rc1c = {2006, 2007, 2008, 2009};
  std::array<Vertex, 4> v_rc1c{};
  std::array<Edge, 4> e_rc1c{};
  for (std::size_t i = 0; i < 4; ++i) {
    v_rc1c[i] = add_vertex_and_check(g, vp_rc1c[i], std::to_string(vp_rc1c[i]));
    bool edge_added_success = false;
    std::tie(e_rc1c[i], edge_added_success) =
        add_edge_and_check(g, v_rc[0], v_rc1c[i], ep_rc1c[i], std::to_string(ep_rc1c[i]));
    EXPECT_TRUE(edge_added_success);
  }
  EXPECT_EQ(num_vertices(g), 9);

  EXPECT_EQ(vertex_id_map[v_root], 1);
  {
    /* IncidenceGraph */
    {
      auto [e_ids, v_ids] = this->get_graph_out_edge_ids(g, v_root, 4);
      EXPECT_THAT(e_ids, ::testing::UnorderedElementsAreArray({1002, 1003, 1004, 1005}));
      EXPECT_THAT(v_ids, ::testing::UnorderedElementsAreArray({2, 3, 4, 5}));
    }

    /* BidirectionalGraph */
    if constexpr (is_bidirectional_graph_v<Graph>) {
      if constexpr (is_undirected_graph_v<Graph>) {
        auto [e_ids, v_ids] = this->get_graph_in_edge_ids(g, v_rc[0], 5);
        EXPECT_THAT(e_ids, ::testing::UnorderedElementsAreArray({1002, 2006, 2007, 2008, 2009}));
        EXPECT_THAT(v_ids, ::testing::UnorderedElementsAreArray({1, 6, 7, 8, 9}));
      } else {
        auto [e_ids, v_ids] = this->get_graph_in_edge_ids(g, v_rc[0], 1);
        EXPECT_THAT(e_ids, ::testing::UnorderedElementsAreArray({1002}));
        EXPECT_THAT(v_ids, ::testing::UnorderedElementsAreArray({1}));
      }
    }

    /* IncidenceGraph */
    if constexpr (is_undirected_graph_v<Graph>) {
      auto [e_ids, v_ids] = this->get_graph_out_edge_ids(g, v_rc[0], 5);
      EXPECT_THAT(e_ids, ::testing::UnorderedElementsAreArray({1002, 2006, 2007, 2008, 2009}));
      EXPECT_THAT(v_ids, ::testing::UnorderedElementsAreArray({1, 6, 7, 8, 9}));
    } else {
      auto [e_ids, v_ids] = this->get_graph_out_edge_ids(g, v_rc[0], 4);
      EXPECT_THAT(e_ids, ::testing::UnorderedElementsAreArray({2006, 2007, 2008, 2009}));
      EXPECT_THAT(v_ids, ::testing::UnorderedElementsAreArray({6, 7, 8, 9}));
    }
  }

  /* MutablePropertyGraph (with rvalue-ref) */
  std::array<int, 4> vp_rc2c = {10, 11, 12, 13};
  std::array<int, 4> ep_rc2c = {3010, 3011, 3012, 3013};
  std::array<Vertex, 4> v_rc2c{};
  std::array<Edge, 4> e_rc2c{};
  for (std::size_t i = 0; i < 4; ++i) {
    v_rc2c[i] = add_vertex_and_check(g, std::move(vp_rc2c[i]));
    bool e_added = false;
    std::tie(e_rc2c[i], e_added) = add_edge_and_check(g, v_rc[1], v_rc2c[i], ep_rc2c[i], std::to_string(ep_rc2c[i]));
    EXPECT_TRUE(e_added);
  }
  EXPECT_EQ(num_vertices(g), 13);

  /* IncidenceGraph */
  if constexpr (is_undirected_graph_v<Graph>) {
    auto [e_ids, v_ids] = this->get_graph_out_edge_ids(g, v_rc[1], 5);
    EXPECT_THAT(e_ids, ::testing::UnorderedElementsAreArray({1003, 3010, 3011, 3012, 3013}));
    EXPECT_THAT(v_ids, ::testing::UnorderedElementsAreArray({1, 10, 11, 12, 13}));
  } else {
    auto [e_ids, v_ids] = this->get_graph_out_edge_ids(g, v_rc[1], 4);
    EXPECT_THAT(e_ids, ::testing::UnorderedElementsAreArray({3010, 3011, 3012, 3013}));
    EXPECT_THAT(v_ids, ::testing::UnorderedElementsAreArray({10, 11, 12, 13}));
  }

  /* Copying function */
  this->check_fullbranch_integrity(g, v_root);
  {
    Graph* p_g_cpy = nullptr;
    EXPECT_NO_THROW(p_g_cpy = new Graph(g));
    this->check_fullbranch_integrity(*p_g_cpy, *(vertices(*p_g_cpy).begin()));
    EXPECT_NO_THROW(delete p_g_cpy);
  }

  {
    Graph g_cpy;
    EXPECT_NO_THROW(g_cpy = g);
    this->check_fullbranch_integrity(g_cpy, *(vertices(g_cpy).begin()));
  }

  {
    Graph* p_g_mv = nullptr;
    EXPECT_NO_THROW(p_g_mv = new Graph(std::move(g)));
    this->check_fullbranch_integrity(*p_g_mv, *(vertices(*p_g_mv).begin()));
    EXPECT_NO_THROW(g = std::move(*p_g_mv));
    this->check_fullbranch_integrity(g, *(vertices(g).begin()));
    EXPECT_NO_THROW(delete p_g_mv);
    v_root = *(vertices(g).begin());
  }

  /* MutableGraph */
  EXPECT_NO_THROW(clear_vertex(v_rc[0], g));
  check_vertex_cleared(g, v_rc[0]);

  if constexpr (is_undirected_graph_v<Graph>) {
    /* IncidenceGraph */
    EXPECT_EQ(out_degree(v_rc[0], g), 0);
    EXPECT_EQ(out_degree(v_root, g), 3);

    /* BidirectionalGraph */
    if constexpr (is_bidirectional_graph_v<Graph>) {
      EXPECT_EQ(in_degree(v_rc[0], g), 0);
      EXPECT_EQ(in_degree(v_rc1c[0], g), 0);
      EXPECT_EQ(in_degree(v_rc1c[1], g), 0);
      EXPECT_EQ(in_degree(v_rc1c[2], g), 0);
      EXPECT_EQ(in_degree(v_rc1c[3], g), 0);
    }
  } else {
    /* IncidenceGraph */
    EXPECT_EQ(out_degree(v_rc[0], g), 0);
    EXPECT_EQ(out_degree(v_root, g), 3);

    /* BidirectionalGraph */
    if constexpr (is_bidirectional_graph_v<Graph>) {
      EXPECT_EQ(in_degree(v_rc[0], g), 0);
      EXPECT_EQ(in_degree(v_rc1c[0], g), 0);
      EXPECT_EQ(in_degree(v_rc1c[1], g), 0);
      EXPECT_EQ(in_degree(v_rc1c[2], g), 0);
      EXPECT_EQ(in_degree(v_rc1c[3], g), 0);
    }
  }

  /* VertexListGraph */
  EXPECT_EQ(num_vertices(g), 13);
  {
    auto v_ids = get_graph_vertex_ids(g, 13);
    EXPECT_THAT(v_ids, ::testing::UnorderedElementsAreArray({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}));
  }

  /* EdgeListGraph */
  EXPECT_EQ(num_edges_or_count(g), 7);
  if constexpr (concepts::EdgeListGraph<Graph>) {
    auto e_ids = get_graph_edge_ids(g, 7);
    EXPECT_THAT(e_ids, ::testing::UnorderedElementsAreArray({1003, 1004, 1005, 3010, 3011, 3012, 3013}));
  }

  /* MutableGraph */
  EXPECT_NO_THROW(remove_edge(v_rc[1], v_rc2c[2], g));

  /* VertexListGraph */
  EXPECT_EQ(num_vertices(g), 13);
  {
    auto v_ids = get_graph_vertex_ids(g, 13);
    EXPECT_THAT(v_ids, ::testing::UnorderedElementsAreArray({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}));
  }

  /* EdgeListGraph */
  EXPECT_EQ(num_edges_or_count(g), 6);
  if constexpr (concepts::EdgeListGraph<Graph>) {
    auto e_ids = get_graph_edge_ids(g, 6);
    EXPECT_THAT(e_ids, ::testing::UnorderedElementsAreArray({1003, 1004, 1005, 3010, 3011, 3013}));
  }

  /* MutableGraph */
  std::pair<Edge, bool> last_e_of_rc2;
  EXPECT_NO_THROW(last_e_of_rc2 = edge(v_rc[1], v_rc2c[3], g));
  EXPECT_TRUE(last_e_of_rc2.second);
  EXPECT_NO_THROW(remove_edge(last_e_of_rc2.first, g));

  /* VertexListGraph */
  EXPECT_EQ(num_vertices(g), 13);
  {
    auto v_ids = get_graph_vertex_ids(g, 13);
    EXPECT_THAT(v_ids, ::testing::UnorderedElementsAreArray({1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13}));
  }

  /* EdgeListGraph */
  EXPECT_EQ(num_edges_or_count(g), 5);
  if constexpr (concepts::EdgeListGraph<Graph>) {
    auto e_ids = get_graph_edge_ids(g, 5);
    EXPECT_THAT(e_ids, ::testing::UnorderedElementsAreArray({1003, 1004, 1005, 3010, 3011}));
  }

  /* MutableGraph */
  EXPECT_NO_THROW(clear_vertex(v_rc2c[0], g));
  EXPECT_NO_THROW(remove_vertex(v_rc2c[0], g));

  /* VertexListGraph */
  EXPECT_EQ(num_vertices(g), 12);
  {
    auto v_ids = get_graph_vertex_ids(g, 12);
    EXPECT_THAT(v_ids, ::testing::UnorderedElementsAreArray({1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13}));
  }

  /* EdgeListGraph */
  EXPECT_EQ(num_edges_or_count(g), 4);
  if constexpr (concepts::EdgeListGraph<Graph>) {
    auto e_ids = get_graph_edge_ids(g, 4);
    EXPECT_THAT(e_ids, ::testing::UnorderedElementsAreArray({1003, 1004, 1005, 3011}));
  }
}

REGISTER_TYPED_TEST_SUITE_P(GraphTreeMutationTest, MutateTree);

}  // namespace bagl

#endif  // BAGL_TEST_GRAPH_MUTATION_TEST_SUITE_H_
