// Copyright 2005 The Trustees of Indiana University.
//  Authors: Jeremiah Willcock
//           Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <algorithm>
#include <ctime>
#include <iostream>
#include <limits>
#include <random>
#include <ranges>
#include <string>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/betweenness_centrality.h"
#include "bagl/compressed_sparse_row_graph.h"
#include "bagl/erdos_renyi_generator.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/kruskal_min_spanning_tree.h"
#include "bagl/properties.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

typedef adjacency_list<> GraphT;

struct VertexData {
  int index;
};

struct EdgeData {
  int index_e;
};

typedef compressed_sparse_row_graph<directed_s, VertexData, EdgeData> CSRGraphT;
typedef compressed_sparse_row_graph<bidirectional_s, VertexData, EdgeData> BidirCSRGraphT;

template <class G1, class G2>
void assert_graphs_equal(const G1& g1, const G2& g2) {
  EXPECT_EQ(num_vertices(g1), num_vertices(g2));
  EXPECT_EQ(num_edges(g1), num_edges(g2));

  for (auto v1 : vertices(g1)) {
    auto v2 = v1;

    EXPECT_EQ(v1, v2);

    EXPECT_EQ(out_degree(v1, g1), out_degree(v2, g2));
    std::vector<std::size_t> edges1;
    for (auto e1 : out_edges(v1, g1)) {
      EXPECT_EQ(source(e1, g1), v1);
      edges1.push_back(target(e1, g1));
    }
    std::vector<std::size_t> edges2;
    for (auto e2 : out_edges(v2, g2)) {
      EXPECT_EQ(source(e2, g2), v2);
      edges2.push_back(target(e2, g2));
    }
    EXPECT_THAT(edges1, ::testing::UnorderedElementsAreArray(edges2));
  }

  std::vector<std::pair<std::size_t, std::size_t>> all_edges1;
  std::vector<std::pair<std::size_t, std::size_t>> all_edges2;
  for (auto e1 : edges(g1)) {
    all_edges1.emplace_back(source(e1, g1), target(e1, g1));
  }
  for (auto e2 : edges(g2)) {
    all_edges2.emplace_back(source(e2, g2), target(e2, g2));
  }
  EXPECT_THAT(all_edges1, ::testing::UnorderedElementsAreArray(all_edges2));
}

template <typename Structure>
void check_consistency_one(const Structure& g) {
  // Do a bunch of tests on the graph internal data
  // Check that rowstart_ entries are valid, and that entries after
  // m_last_source + 1 are all zero
  EXPECT_EQ(g.rowstart_[0], 0);
  for (std::size_t i = 0; i < g.rowstart_.size() - 1; ++i) {
    EXPECT_GE(g.rowstart_[i + 1], g.rowstart_[i]);
    EXPECT_LE(g.rowstart_[i + 1], g.rowstart_.back());
  }
  // Check that column_ entries are within range
  for (std::size_t i = 0; i < g.rowstart_.back(); ++i) {
    EXPECT_LT(g.column_[i], g.rowstart_.size() - 1);
  }
}

template <typename G>
void assert_bidir_equal_in_both_dirs(const G& g) {
  EXPECT_EQ(g.forward_.rowstart_.size(), g.backward_.rowstart_.size());
  EXPECT_EQ(g.forward_.column_.size(), g.backward_.column_.size());
  typedef graph_vertex_descriptor_t<G> Vertex;
  std::vector<std::tuple<std::size_t, Vertex, Vertex>> edges_forward;
  std::vector<std::tuple<std::size_t, Vertex, Vertex>> edges_backward;
  for (Vertex i = 0; i < g.forward_.rowstart_.size() - 1; ++i) {
    for (std::size_t j = g.forward_.rowstart_[i]; j < g.forward_.rowstart_[i + 1]; ++j) {
      edges_forward.emplace_back(j, i, g.forward_.column_[j]);
    }
  }
  for (Vertex i = 0; i < g.backward_.rowstart_.size() - 1; ++i) {
    for (std::size_t j = g.backward_.rowstart_[i]; j < g.backward_.rowstart_[i + 1]; ++j) {
      edges_backward.emplace_back(g.backward_.edge_properties_[j], g.backward_.column_[j], i);
    }
  }
  EXPECT_THAT(edges_forward, ::testing::UnorderedElementsAreArray(edges_backward));
}

template <typename Graph>
void check_consistency(const Graph& g) {
  check_consistency_one(g.forward_);
  if constexpr (is_bidirectional_graph_v<Graph>) {
    check_consistency_one(g.backward_);
    assert_bidir_equal_in_both_dirs(g);
  }
}

template <typename OrigGraph>
void graph_test(const OrigGraph& g) {
  // Check copying of a graph
  CSRGraphT g2(g);
  check_consistency(g2);
  EXPECT_EQ(std::ranges::distance(edges(g2)), num_edges(g2));
  assert_graphs_equal(g, g2);

  // Check constructing a graph from iterators
  CSRGraphT g3(edges_are_sorted,
               edges(g2) | std::views::transform(csr_detail::edge_to_index_pair(g2, get(vertex_index, g2))),
               num_vertices(g));
  check_consistency(g3);
  EXPECT_EQ(std::ranges::distance(edges(g3)), num_edges(g3));
  assert_graphs_equal(g2, g3);

  // Check constructing a graph using in-place modification of vectors
  {
    std::vector<std::size_t> sources(num_edges(g2));
    std::vector<std::size_t> targets(num_edges(g2));
    std::size_t idx = 0;
    // Edges actually sorted
    for (auto e : edges(g2)) {
      sources[idx] = source(e, g2);
      targets[idx] = target(e, g2);
      ++idx;
    }
    CSRGraphT g3a(construct_inplace_from_sources_and_targets, sources, targets, num_vertices(g2));
    check_consistency(g3a);
    assert_graphs_equal(g2, g3a);
  }
  {
    std::vector<std::size_t> sources(num_edges(g2));
    std::vector<std::size_t> targets(num_edges(g2));
    std::size_t idx = 0;
    // Edges reverse-sorted
    for (auto e : edges(g2)) {
      sources[num_edges(g2) - 1 - idx] = source(e, g2);
      targets[num_edges(g2) - 1 - idx] = target(e, g2);
      ++idx;
    }
    CSRGraphT g3a(construct_inplace_from_sources_and_targets, sources, targets, num_vertices(g2));
    check_consistency(g3a);
    assert_graphs_equal(g2, g3a);
  }
  {
    std::vector<std::size_t> sources(num_edges(g2));
    std::vector<std::size_t> targets(num_edges(g2));
    std::size_t idx = 0;
    // Edges scrambled using Fisher-Yates shuffle (Durstenfeld variant) from wikipedia
    for (auto e : edges(g2)) {
      sources[idx] = source(e, g2);
      targets[idx] = target(e, g2);
      ++idx;
    }
    std::mt19937 gen(1);
    if (num_edges(g) != 0) {
      for (std::size_t i = num_edges(g) - 1; i > 0; --i) {
        std::size_t scrambled = std::uniform_int_distribution<>(0, i)(gen);
        if (scrambled == i) {
          continue;
        }
        using std::swap;
        swap(sources[i], sources[scrambled]);
        swap(targets[i], targets[scrambled]);
      }
    }
    CSRGraphT g3a(construct_inplace_from_sources_and_targets, sources, targets, num_vertices(g2));
    check_consistency(g3a);
    assert_graphs_equal(g2, g3a);
  }

  // Check edge_from_index (and implicitly the edge_index property map) for
  // each edge in g2
  std::size_t last_src = 0;
  for (auto e : edges(g2)) {
    EXPECT_EQ(edge_from_index(get(edge_index, g2, e), g2), e);
    std::size_t src = get(vertex_index, g2, source(e, g2));
    (void)get(vertex_index, g2, target(e, g2));
    EXPECT_GE(src, last_src);
    last_src = src;
  }

  // Check out edge iteration and vertex iteration for sortedness
  std::size_t last_vertex = 0;
  bool first_iter = true;
  for (auto v2 : vertices(g2)) {
    std::size_t v = get(vertex_index, g2, v2);
    if (!first_iter) {
      EXPECT_GT(v, last_vertex);
    }
    last_vertex = v;
    first_iter = false;

    for (auto e : out_edges(v2, g2)) {
      EXPECT_EQ(source(e, g2), v2);
    }
  }

  // Run brandes_betweenness_centrality, which touches on a whole lot
  // of things, including VertexListGraph and IncidenceGraph
  auto vertex_centralities = vector_property_map(num_vertices(g3), get(vertex_index, g3), double{0.0});
  auto edge_centralities = vector_property_map(num_edges(g3), get(edge_index, g3), double{0.0});
  brandes_betweenness_centrality(g3, vertex_centralities.ref(), edge_centralities.ref());
  // Extra qualifications for aCC

  // Invert the edge centralities and use these as weights to
  // Kruskal's MST algorithm, which will test the EdgeListGraph
  // capabilities.
  double max_val = std::numeric_limits<double>::max();
  for (auto e : edges(g3)) {
    edge_centralities[e] = (edge_centralities[e] == 0.0 ? max_val : 1.0 / edge_centralities[e]);
  }

  using Edge = graph_edge_descriptor_t<CSRGraphT>;
  std::vector<Edge> mst_edges;
  mst_edges.reserve(num_vertices(g3));
  kruskal_minimum_spanning_tree(g3, std::back_inserter(mst_edges), edge_centralities.ref());
}

void graph_test(int nnodes, double density, int seed) {
  std::mt19937 gen(seed);
  GraphT g(nnodes, erdos_renyi_range<is_undirected_graph_v<GraphT>>(gen, nnodes, density));
  graph_test(g);
}

TEST(CompressedSparseRowGraphTest, GraphProp) {
  using GraphPropCSRGraphT =
      compressed_sparse_row_graph<directed_s, no_property, no_property, property<graph_name_t, std::string>>;

  GraphPropCSRGraphT g;
  EXPECT_EQ(get_property(g, graph_name), "");
  set_property(g, graph_name, "beep");
  EXPECT_EQ(get_property(g, graph_name), "beep");
}

struct Vertex {
  double centrality;
};

struct Edge {
  Edge(double w) : weight(w), centrality(0.0) {}

  double weight;
  double centrality;
};

TEST(CompressedSparseRowGraphTest, VertexAndEdgeProp) {
  using CSRGraphWithPropsT = compressed_sparse_row_graph<directed_s, Vertex, Edge>;

  std::vector<std::pair<int, int>> edges_init = {{0, 1}, {0, 3}, {1, 2}, {3, 1}, {3, 4}, {4, 2}};
  std::vector<double> weights = {1.0, 1.0, 0.5, 1.0, 1.0, 0.5};
  std::vector<double> centrality = {0.0, 1.5, 0.0, 1.0, 0.5};

  CSRGraphWithPropsT g(edges_are_sorted, edges_init, weights, 5, 6);
  brandes_betweenness_centrality_weighted(g, get(&Vertex::centrality, g), get(&Edge::weight, g),
                                          get(&Edge::centrality, g), get(vertex_index, g));

  for (auto v : vertices(g)) {
    EXPECT_EQ(g[v].centrality, centrality[v]) << "At vertex " << get(vertex_index, g, v);
  }
}

TEST(CompressedSparseRowGraphTest, EmptyGraph) {
  CSRGraphT g;
  graph_test(g);
}

// Optionally accept a seed value
constexpr int test_seed = 42;

TEST(CompressedSparseRowGraphTest, Random_1000_001) { graph_test(1000, 0.001, test_seed); }

TEST(CompressedSparseRowGraphTest, Random_1000_0005) { graph_test(1000, 0.0005, test_seed); }

TEST(CompressedSparseRowGraphTest, SmallGraph) {
  std::vector<std::pair<int, int>> unsorted_edges = {{5, 0}, {3, 2}, {4, 1}, {4, 0}, {0, 2}, {5, 2}};
  CSRGraphT g(edges_are_unsorted, unsorted_edges, 6);

  // Test vertex and edge bundle access
  auto v1 = *vertices(g).begin();
  (void)static_cast<VertexData&>(get(get(vertex_bundle, g), v1));
  (void)static_cast<const VertexData&>(get(get(vertex_bundle, std::as_const(g)), v1));
  (void)static_cast<VertexData&>(get(vertex_bundle, g, v1));
  (void)static_cast<const VertexData&>(get(vertex_bundle, std::as_const(g), v1));
  put(vertex_bundle, g, v1, VertexData{});
  auto e1 = *edges(g).begin();
  (void)static_cast<EdgeData&>(get(get(edge_bundle, g), e1));
  (void)static_cast<const EdgeData&>(get(get(edge_bundle, std::as_const(g)), e1));
  (void)static_cast<EdgeData&>(get(edge_bundle, g, e1));
  (void)static_cast<const EdgeData&>(get(edge_bundle, std::as_const(g), e1));
  put(edge_bundle, g, e1, EdgeData{});

  CSRGraphT g2(edges_are_unsorted_multi_pass, unsorted_edges, 6);
  graph_test(g);
  graph_test(g2);
  assert_graphs_equal(g, g2);
  BidirCSRGraphT g2b(edges_are_unsorted_multi_pass, unsorted_edges, 6);
  graph_test(g2b);
  assert_graphs_equal(g, g2b);
  // Check in edge access
  auto ie_rg = in_edges(*vertices(g2b).begin(), g2b);
  static_assert(std::ranges::input_range<decltype(ie_rg)>);
  (void)ie_rg;

  // Test building a graph using add_edges on unsorted lists
  CSRGraphT g3(edges_are_unsorted, std::ranges::empty_view<std::pair<int, int>>(), 6);
  add_edges(std::ranges::ref_view(unsorted_edges) | std::views::take(3), g3);
  std::array<EdgeData, 3> edge_data;
  add_edges(std::ranges::ref_view(unsorted_edges) | std::views::drop(3), edge_data, g3);
  graph_test(g3);
  assert_graphs_equal(g, g3);
}

}  // namespace
}  // namespace bagl
