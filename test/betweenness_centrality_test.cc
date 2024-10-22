// Copyright 2004 The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/betweenness_centrality.h"

#include <queue>
#include <random>
#include <stack>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/compressed_sparse_row_graph.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

const double error_tolerance = 0.0025;

struct weighted_edge {
  int source;
  int target;
  double weight;
};

template <typename T>
class BetweennessCentralityDirectedTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(BetweennessCentralityDirectedTest);

TYPED_TEST_P(BetweennessCentralityDirectedTest, WeightedEdges) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g(5);

  std::vector<Vertex> v_by_index{num_vertices(g)};
  int index = 0;
  for (auto v : vertices(g)) {
    put(vertex_index, g, v, index);
    v_by_index[index] = v;
    ++index;
  }

  weighted_edge edge_init[6] = {{0, 1, 1.0}, {0, 3, 1.0}, {1, 2, 0.5}, {3, 1, 1.0}, {3, 4, 1.0}, {4, 2, 0.5}};
  for (auto [src, tgt, wght] : edge_init) {
    auto [e, e_added] = add_edge(v_by_index[src], v_by_index[tgt], g);
    ASSERT_TRUE(e_added);
    put(edge_weight, g, e, wght);  // 1.0);
  }

  auto centrality = vector_property_map(num_vertices(g), get(vertex_index, g), double{0.0});
  brandes_betweenness_centrality_weighted(g, centrality.ref(), get(edge_weight, g));

  double correct_centrality[5] = {0.0, 1.5, 0.0, 1.0, 0.5};
  for (auto v : vertices(g)) {
    EXPECT_EQ(centrality[v], correct_centrality[get(vertex_index, g, v)]);
  }
}

REGISTER_TYPED_TEST_SUITE_P(BetweennessCentralityDirectedTest, WeightedEdges);

struct unweighted_edge {
  int source;
  int target;
};

template <typename T>
class BetweennessCentralityUndirectedTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(BetweennessCentralityUndirectedTest);

template <typename Graph>
void run_unweighted_test(std::span<unweighted_edge> edge_init, std::span<double> correct_centrality,
                         std::span<double> correct_edge_centrality) {
  ASSERT_TRUE(correct_edge_centrality.empty() || correct_edge_centrality.size() == edge_init.size());
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g(correct_centrality.size());

  std::vector<Vertex> v_by_index{num_vertices(g)};
  int index = 0;
  for (auto v : vertices(g)) {
    put(vertex_index, g, v, index);
    v_by_index[index] = v;
    ++index;
  }

  index = 0;
  for (auto [src, tgt] : edge_init) {
    auto [e, e_added] = add_edge(v_by_index[src], v_by_index[tgt], g);
    ASSERT_TRUE(e_added);
    put(edge_weight, g, e, 1.0);
    put(edge_index, g, e, index);
    ++index;
  }

  auto vindex = get(vertex_index, g);
  auto eindex = get(edge_index, g);

  auto centrality1 = vector_property_map(num_vertices(g), vindex, double{0.0});
  auto edge_centrality1 = vector_property_map(edge_init.size(), eindex, double{0.0});
  brandes_betweenness_centrality(g, centrality1.ref(), edge_centrality1.ref());

  auto centrality2 = vector_property_map(num_vertices(g), vindex, double{0.0});
  auto edge_centrality2 = vector_property_map(edge_init.size(), eindex, double{0.0});
  brandes_betweenness_centrality_weighted(g, centrality2.ref(), edge_centrality2.ref(), get(edge_weight, g), vindex);

  auto edge_centrality3 = vector_property_map(edge_init.size(), eindex, double{0.0});
  brandes_betweenness_centrality(g, dummy_property_map(), edge_centrality3.ref());

  for (auto v : vertices(g)) {
    EXPECT_THAT(centrality1[v], ::testing::DoubleEq(centrality2[v]));
    EXPECT_THAT(centrality1[v], ::testing::DoubleNear(correct_centrality[vindex[v]], error_tolerance))
        << "Vertex: " << vindex[v];
  }

  for (auto e : edges(g)) {
    EXPECT_THAT(edge_centrality1[e], ::testing::DoubleEq(edge_centrality2[e]));
    EXPECT_THAT(edge_centrality1[e], ::testing::DoubleEq(edge_centrality3[e]));
    if (!correct_edge_centrality.empty()) {
      EXPECT_THAT(edge_centrality1[e], ::testing::DoubleNear(correct_edge_centrality[eindex[e]], error_tolerance))
          << "Edge: " << eindex[e];
    }
  }
}

TYPED_TEST_P(BetweennessCentralityUndirectedTest, UnweightedEdges1) {
  unweighted_edge ud_edge_init1[5] = {{0, 1}, {0, 3}, {1, 2}, {3, 2}, {2, 4}};
  double ud_centrality1[5] = {0.5, 1.0, 3.5, 1.0, 0.0};
  run_unweighted_test<TypeParam>(ud_edge_init1, ud_centrality1, {});
}

TYPED_TEST_P(BetweennessCentralityUndirectedTest, UnweightedEdges2) {
  // Example borrowed from the JUNG test suite
  unweighted_edge ud_edge_init2[10] = {
      {0, 1}, {0, 6}, {1, 2}, {1, 3}, {2, 4}, {3, 4}, {4, 5}, {5, 8}, {7, 8}, {6, 7},
  };
  double ud_centrality2[9] = {0.2142 * 28, 0.2797 * 28, 0.0892 * 28, 0.0892 * 28, 0.2797 * 28,
                              0.2142 * 28, 0.1666 * 28, 0.1428 * 28, 0.1666 * 28};
  double ud_edge_centrality2[10] = {10.66666, 9.33333, 6.5, 6.5, 6.5, 6.5, 10.66666, 9.33333, 8.0, 8.0};

  run_unweighted_test<TypeParam>(ud_edge_init2, ud_centrality2, ud_edge_centrality2);
}

TYPED_TEST_P(BetweennessCentralityUndirectedTest, Wheel) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g(15);
  Vertex center = *vertices(g).begin();

  std::vector<Vertex> v_by_index(num_vertices(g));
  int index = 0;
  for (auto v : vertices(g)) {
    put(vertex_index, g, v, index);
    v_by_index[index] = v;
    if (v != center) {
      auto [e, e_added] = add_edge(v, center, g);
      ASSERT_TRUE(e_added);
      put(edge_weight, g, e, 1.0);
    }
    ++index;
  }

  auto vindex = get(vertex_index, g);

  auto centrality1 = vector_property_map(num_vertices(g), vindex, double{0.0});
  brandes_betweenness_centrality(g, centrality1.ref());

  auto centrality2 = vector_property_map(num_vertices(g), vindex, double{0.0});
  brandes_betweenness_centrality_weighted(g, centrality2.ref(), get(edge_weight, g), vindex);

  relative_betweenness_centrality(g, centrality1.ref());
  relative_betweenness_centrality(g, centrality2.ref());

  for (auto v : vertices(g)) {
    EXPECT_THAT(centrality1[v], ::testing::DoubleEq(centrality2[v]));
    if (v == center) {
      EXPECT_THAT(centrality1[v], ::testing::DoubleEq(1.0));
    } else {
      EXPECT_THAT(centrality1[v], ::testing::DoubleEq(0.0));
    }
  }

  double dominance = central_point_dominance(g, centrality2.ref());
  EXPECT_THAT(dominance, ::testing::DoubleEq(1.0));
}

template <typename MutableGraph>
void randomly_add_edges(MutableGraph& g, double edge_probability) {
  std::mt19937 gen(76);
  std::uniform_real_distribution<> uni_dist;

  auto v_rg = vertices(g);
  for (auto vi = v_rg.begin(), vi_end = v_rg.end(); vi != vi_end; ++vi) {
    auto v = *vi;
    auto wi = is_undirected_graph_v<MutableGraph> ? vi : v_rg.begin();
    while (wi != vi_end) {
      auto w = *wi++;
      if (v != w && uni_dist(gen) < edge_probability) {
        add_edge(v, w, g);
      }
    }
  }
}

template <typename Graph, typename VertexIndexMap, typename CentralityMap>
void simple_unweighted_betweenness_centrality(const Graph& g, VertexIndexMap vindex, CentralityMap centrality) {
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using CentralityType = property_traits_value_t<CentralityMap>;

  for (Vertex v : vertices(g)) {
    put(centrality, v, CentralityType{});
  }

  for (Vertex s : vertices(g)) {
    // S <-- empty stack
    std::stack<Vertex> S;

    // P[w] <-- empty list, w \in V
    std::vector<std::vector<Vertex>> predecessors(num_vertices(g));

    // sigma[t] <-- 0, t \in V
    std::vector<std::size_t> sigma(num_vertices(g), 0);

    // sigma[s] <-- 1
    sigma[get(vindex, s)] = 1;

    // d[t] <-- -1, t \in V
    std::vector<int> d(num_vertices(g), -1);

    // d[s] <-- 0
    d[get(vindex, s)] = 0;

    // Q <-- empty queue
    std::queue<Vertex> Q;

    // enqueue s --> Q
    Q.push(s);

    while (!Q.empty()) {
      // dequeue v <-- Q
      auto v = Q.front();
      Q.pop();

      // push v --> S
      S.push(v);

      for (Vertex w : adjacent_vertices(v, g)) {
        // w found for the first time?
        if (d[get(vindex, w)] < 0) {
          // enqueue w --> Q
          Q.push(w);

          // d[w] <-- d[v] + 1
          d[get(vindex, w)] = d[get(vindex, v)] + 1;
        }

        // shortest path to w via v?
        if (d[get(vindex, w)] == d[get(vindex, v)] + 1) {
          // sigma[w] = sigma[w] + sigma[v]
          sigma[get(vindex, w)] += sigma[get(vindex, v)];

          // append v --> P[w]
          predecessors[get(vindex, w)].push_back(v);
        }
      }
    }

    // delta[v] <-- 0, v \in V
    std::vector<CentralityType> delta(num_vertices(g), CentralityType{});

    // S returns vertices in order of non-increasing distance from s
    while (!S.empty()) {
      // pop w <-- S
      Vertex w = S.top();
      S.pop();

      const auto& w_preds = predecessors[get(vindex, w)];
      for (Vertex v : w_preds) {
        // delta[v] <-- delta[v] + (sigma[v]/sigma[w])*(1 + delta[w])
        delta[get(vindex, v)] +=
            (static_cast<CentralityType>(sigma[get(vindex, v)]) / sigma[get(vindex, w)]) * (1 + delta[get(vindex, w)]);
      }

      if (w != s) {
        // C_B[w] <-- C_B[w] + delta[w]
        centrality[w] += delta[get(vindex, w)];
      }
    }
  }

  if constexpr (is_undirected_graph_v<Graph>) {
    for (Vertex v : vertices(g)) {
      put(centrality, v, get(centrality, v) / CentralityType{2});
    }
  }
}

TYPED_TEST_P(BetweennessCentralityUndirectedTest, RandomUnweightedEdges) {
  using Graph = TypeParam;
  int n = 300;
  Graph g(n);

  int index = 0;
  for (auto v : vertices(g)) {
    put(vertex_index, g, v, index++);
  }

  randomly_add_edges(g, 0.20);

  // Direct translation of Brandes' algorithm...
  auto centrality1 = vector_property_map(num_vertices(g), get(vertex_index, g), double{0.0});
  simple_unweighted_betweenness_centrality(g, get(vertex_index, g), centrality1.ref());

  // Real version, unweighted...
  auto centrality2 = vector_property_map(num_vertices(g), get(vertex_index, g), double{0.0});
  brandes_betweenness_centrality(g, centrality2.ref());

  for (auto v : vertices(g)) {
    EXPECT_THAT(centrality2[v], ::testing::DoubleNear(centrality1[v], error_tolerance))
        << "Vertex: " << get(vertex_index, g, v);
  }

  for (auto e : edges(g)) {
    put(edge_weight, g, e, 1.0);
  }

  auto centrality3 = vector_property_map(num_vertices(g), get(vertex_index, g), double{0.0});
  brandes_betweenness_centrality_weighted(g, centrality3.ref(), get(edge_weight, g), get(vertex_index, g));

  for (auto v : vertices(g)) {
    EXPECT_THAT(centrality3[v], ::testing::DoubleNear(centrality1[v], error_tolerance))
        << "Vertex: " << get(vertex_index, g, v);
  }
}

REGISTER_TYPED_TEST_SUITE_P(BetweennessCentralityUndirectedTest, UnweightedEdges1, UnweightedEdges2, Wheel,
                            RandomUnweightedEdges);

using TestEdgeProp = property<edge_weight_t, double, property<edge_index_t, std::size_t>>;
using TestGraph = adjacency_list<list_s, list_s, undirected_s, property<vertex_index_t, int>, TestEdgeProp>;
using TestDigraph = adjacency_list<list_s, list_s, directed_s, property<vertex_index_t, int>, TestEdgeProp>;

INSTANTIATE_TYPED_TEST_SUITE_P(BetweennessCentralityDirected, BetweennessCentralityDirectedTest, TestDigraph);
INSTANTIATE_TYPED_TEST_SUITE_P(BetweennessCentralityUndirected, BetweennessCentralityUndirectedTest, TestGraph);

}  // namespace
}  // namespace bagl
