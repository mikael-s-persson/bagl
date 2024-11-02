// Copyright (C) 2012 Flavio De Lorenzi (fdlorenzi@gmail.com)
// Revision History:
//   8 April 2013: Fixed a typo in random_functor. (Flavio De Lorenzi)
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/vf2_sub_graph_iso.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <random>

#include "bagl/adjacency_list.h"
#include "bagl/bfl_d_ary_tree.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/property_map.h"
#include "bagl/random.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

struct get_hit_and_stop {
  explicit get_hit_and_stop(bool& got_hit, bool stop = true) : got_hit_(&got_hit), stop_(stop) {}

  template <typename Map1To2, typename Map2To1>
  bool operator()(Map1To2 /*unused*/, Map2To1 /*unused*/) {
    *got_hit_ = true;
    return stop_;
  }

  bool* got_hit_;
  bool stop_;
};

struct false_predicate {
  template <typename VertexOrEdge1, typename VertexOrEdge2>
  bool operator()(VertexOrEdge1 /*unused*/, VertexOrEdge2 /*unused*/) const {
    return false;
  }
};

using TestGraph = adjacency_list<vec_s, vec_s, bidirectional_s>;

TEST(GraphIsomorphismTest, Empty) {
  TestGraph g_empty;
  bool got_hit = false;
  EXPECT_TRUE(vf2_graph_iso(g_empty, g_empty, get_hit_and_stop{got_hit}));
  // even empty matches are reported
  EXPECT_TRUE(got_hit);
}

TEST(SubGraphIsomorphismTest, EmptyVsEmpty) {
  TestGraph g_empty;
  bool got_hit = false;
  EXPECT_TRUE(vf2_subgraph_iso(g_empty, g_empty, get_hit_and_stop{got_hit}));
  // even empty matches are reported
  EXPECT_TRUE(got_hit);
}

TEST(SubGraphIsomorphismTest, EmptyVsNonEmpty) {
  TestGraph g_empty;
  TestGraph g_large{1};
  bool got_hit = false;
  EXPECT_TRUE(vf2_subgraph_iso(g_empty, g_large, get_hit_and_stop{got_hit}));
  // even empty matches are reported
  EXPECT_TRUE(got_hit);
}

TEST(SubGraphMonomorphismTest, EmptyVsEmpty) {
  TestGraph g_empty;
  bool got_hit = false;
  EXPECT_TRUE(vf2_subgraph_mono(g_empty, g_empty, get_hit_and_stop{got_hit}));
  // even empty matches are reported
  EXPECT_TRUE(got_hit);
}

TEST(SubGraphMonomorphismTest, EmptyVsNonEmpty) {
  TestGraph g_empty;
  TestGraph g_large{1};
  bool got_hit = false;
  EXPECT_TRUE(vf2_subgraph_mono(g_empty, g_large, get_hit_and_stop{got_hit}));
  // even empty matches are reported
  EXPECT_TRUE(got_hit);
}

TEST(GraphIsomorphismTest, SizeMismatch) {
  TestGraph g_small{1};
  TestGraph g_large{2, std::ranges::single_view{std::pair{0, 1}}};
  bool got_hit = false;
  EXPECT_FALSE(vf2_graph_iso(g_small, g_large, get_hit_and_stop{got_hit}));
  EXPECT_FALSE(got_hit);
}

TEST(GraphIsomorphismTest, VertexMismatch) {
  TestGraph g_large{2, std::ranges::single_view{std::pair{0, 1}}};
  bool got_hit = false;
  EXPECT_FALSE(vf2_graph_iso(g_large, g_large, get_hit_and_stop{got_hit}, get(vertex_index, g_large),
                             get(vertex_index, g_large), vertex_order_by_mult(g_large), false_predicate{},
                             false_predicate{}));
  EXPECT_FALSE(got_hit);
}

TEST(GraphIsomorphismTest, FindAll) {
  TestGraph g_large{2, std::ranges::single_view{std::pair{0, 1}}};
  bool got_hit = false;
  EXPECT_TRUE(vf2_graph_iso(g_large, g_large, get_hit_and_stop{got_hit, false}));
  EXPECT_TRUE(got_hit);
}

TEST(GraphIsomorphismTest, FindFirst) {
  TestGraph g_large{2, std::ranges::single_view{std::pair{0, 1}}};
  bool got_hit = false;
  EXPECT_TRUE(vf2_graph_iso(g_large, g_large, get_hit_and_stop{got_hit}));
  EXPECT_TRUE(got_hit);
}

TEST(SubGraphIsomorphismTest, SizeMismatch) {
  TestGraph g_small{1};
  TestGraph g_large{2, std::ranges::single_view{std::pair{0, 1}}};
  bool got_hit = false;
  EXPECT_FALSE(vf2_subgraph_iso(g_large, g_small, get_hit_and_stop{got_hit}));
  EXPECT_FALSE(got_hit);
}

TEST(SubGraphIsomorphismTest, VertexMismatch) {
  TestGraph g_large{2, std::ranges::single_view{std::pair{0, 1}}};
  bool got_hit = false;
  EXPECT_FALSE(vf2_subgraph_iso(g_large, g_large, get_hit_and_stop{got_hit}, get(vertex_index, g_large),
                                get(vertex_index, g_large), vertex_order_by_mult(g_large), false_predicate{},
                                false_predicate{}));
  EXPECT_FALSE(got_hit);
}

TEST(SubGraphIsomorphismTest, FindAll) {
  TestGraph g_large{2, std::ranges::single_view{std::pair{0, 1}}};
  bool got_hit = false;
  EXPECT_TRUE(vf2_subgraph_iso(g_large, g_large, get_hit_and_stop{got_hit, false}));
  EXPECT_TRUE(got_hit);
}

TEST(SubGraphIsomorphismTest, FindFirst) {
  TestGraph g_large{2, std::ranges::single_view{std::pair{0, 1}}};
  bool got_hit = false;
  EXPECT_TRUE(vf2_subgraph_iso(g_large, g_large, get_hit_and_stop{got_hit}));
  EXPECT_TRUE(got_hit);
}

TEST(SubGraphMonomorphismTest, SizeMismatch) {
  TestGraph g_small{1};
  TestGraph g_large{2, std::ranges::single_view{std::pair{0, 1}}};
  bool got_hit = false;
  EXPECT_FALSE(vf2_subgraph_mono(g_large, g_small, get_hit_and_stop{got_hit}));
  EXPECT_FALSE(got_hit);
}

TEST(SubGraphMonomorphismTest, VertexMismatch) {
  TestGraph g_large{2, std::ranges::single_view{std::pair{0, 1}}};
  bool got_hit = false;
  EXPECT_FALSE(vf2_subgraph_mono(g_large, g_large, get_hit_and_stop{got_hit}, get(vertex_index, g_large),
                                 get(vertex_index, g_large), vertex_order_by_mult(g_large), false_predicate{},
                                 false_predicate{}));
  EXPECT_FALSE(got_hit);
}

TEST(SubGraphMonomorphismTest, FindAll) {
  TestGraph g_large{2, std::ranges::single_view{std::pair{0, 1}}};
  bool got_hit = false;
  EXPECT_TRUE(vf2_subgraph_mono(g_large, g_large, get_hit_and_stop{got_hit, false}));
  EXPECT_TRUE(got_hit);
}

TEST(SubGraphMonomorphismTest, FindFirst) {
  TestGraph g_large{2, std::ranges::single_view{std::pair{0, 1}}};
  bool got_hit = false;
  EXPECT_TRUE(vf2_subgraph_mono(g_large, g_large, get_hit_and_stop{got_hit}));
  EXPECT_TRUE(got_hit);
}

template <typename Graph1, typename Graph2>
void randomly_permute_graph(Graph1& g1, const Graph2& g2) {
  EXPECT_LE(num_vertices(g1), num_vertices(g2));
  EXPECT_EQ(num_edges(g1), 0);

  using Vertex2 = graph_vertex_descriptor_t<Graph2>;

  std::mt19937 gen{42};

  // Decide new order
  std::vector<Vertex2> orig_vertices;
  std::ranges::copy(vertices(g2), std::back_inserter(orig_vertices));
  std::shuffle(orig_vertices.begin(), orig_vertices.end(), gen);
  auto vertex_map = vector_property_map(num_vertices(g2), get(vertex_index, g2), graph_traits<Graph1>::null_vertex());

  std::size_t i = 0;
  for (auto v : vertices(g1)) {
    vertex_map[orig_vertices[i]] = v;
    put(vertex_name, g1, v, get(vertex_name, g2, orig_vertices[i]));
    ++i;
  }

  for (auto e : edges(g2)) {
    auto si = vertex_map[source(e, g2)];
    auto ti = vertex_map[target(e, g2)];
    if ((si != graph_traits<Graph1>::null_vertex()) && (ti != graph_traits<Graph1>::null_vertex())) {
      EXPECT_TRUE(add_edge(si, ti, g1, get(edge_name, g2, e)).second);
    }
  }
}

template <typename Graph>
void generate_random_digraph(Graph& g, double edge_probability, int max_parallel_edges,
                             double parallel_edge_probability, int max_edge_name, int max_vertex_name) {
  EXPECT_GE(edge_probability, 0);
  EXPECT_LE(edge_probability, 1);
  EXPECT_GE(parallel_edge_probability, 0);
  EXPECT_LE(parallel_edge_probability, 1);
  EXPECT_GE(max_parallel_edges, 0);
  EXPECT_GE(max_edge_name, 0);
  EXPECT_GE(max_vertex_name, 0);

  std::mt19937 gen;
  std::uniform_real_distribution<> dist_real;

  for (auto u : vertices(g)) {
    for (auto v : vertices(g)) {
      if (dist_real(gen) <= edge_probability) {
        EXPECT_TRUE(add_edge(u, v, g).second);
        for (int i = 0; i < max_parallel_edges; ++i) {
          if (dist_real(gen) <= parallel_edge_probability) {
            EXPECT_TRUE(add_edge(u, v, g).second);
          }
        }
      }
    }
  }

  {
    std::uniform_int_distribution<> dist_int(0, max_edge_name);
    randomize_property<vertex_name_t>(g, [&gen, &dist_int]() { return dist_int(gen); });
  }

  {
    std::uniform_int_distribution<> dist_int(0, max_vertex_name);
    randomize_property<edge_name_t>(g, [&gen, &dist_int]() { return dist_int(gen); });
  }
}

template <typename Graph1, typename Graph2, typename EdgeEqPredicate, typename VertexEqPredicate>
struct test_callback {
  test_callback(const Graph1& graph1, const Graph2& graph2, EdgeEqPredicate edge_comp, VertexEqPredicate vertex_comp)
      : graph1_(graph1), graph2_(graph2), edge_comp_(edge_comp), vertex_comp_(vertex_comp) {}

  template <typename CorrMap1To2, typename CorrMap2To1>
  bool operator()(CorrMap1To2 f, CorrMap2To1) {
    bool verified = verify_vf2_subgraph_iso(graph1_, graph2_, f, edge_comp_, vertex_comp_);
    EXPECT_TRUE(verified) << "Graph1:\n"
                          << print_graph_to_string(graph1_) << "\nGraph2:\n"
                          << print_graph_to_string(graph2_);
    return true;
  }

 private:
  const Graph1& graph1_;
  const Graph2& graph2_;
  EdgeEqPredicate edge_comp_;
  VertexEqPredicate vertex_comp_;
};

TEST(SubGraphIsomorphismTest, RandomGraph) {
  int n1 = 10;
  int n2 = 20;
  double edge_probability = 0.4;
  int max_parallel_edges = 2;
  double parallel_edge_probability = 0.4;
  int max_edge_name = 5;
  int max_vertex_name = 5;

  using EdgeProp = property<edge_name_t, int>;
  using VertexProp = property<vertex_name_t, int, property<vertex_index_t, int> >;

  using Graph1 = adjacency_list<list_s, list_s, bidirectional_s, VertexProp, EdgeProp>;
  using Graph2 = adjacency_list<vec_s, vec_s, bidirectional_s, VertexProp, EdgeProp>;

  Graph1 g1(n1);
  Graph2 g2(n2);
  generate_random_digraph(g2, edge_probability, max_parallel_edges, parallel_edge_probability, max_edge_name,
                          max_vertex_name);
  randomly_permute_graph(g1, g2);

  int v_idx = 0;
  for (auto v : vertices(g1)) {
    put(vertex_index, g1, v, v_idx++);
  }

  // Create vertex and edge predicates
  auto vertex_comp = property_map_equivalent(get(vertex_name, g1), get(vertex_name, g2));
  auto edge_comp = property_map_equivalent(get(edge_name, g1), get(edge_name, g2));

  // Create callback
  test_callback callback(g1, g2, edge_comp, vertex_comp);

  EXPECT_TRUE(vf2_subgraph_iso(g1, g2, std::ref(callback), get(vertex_index, g1), get(vertex_index, g2),
                               vertex_order_by_mult(g1), edge_comp, vertex_comp));

  if (num_vertices(g1) == num_vertices(g2)) {
    EXPECT_TRUE(vf2_graph_iso(g1, g2, std::ref(callback), get(vertex_index, g1), get(vertex_index, g2),
                              vertex_order_by_mult(g1), edge_comp, vertex_comp));
  }

  if (::testing::Test::HasFailure()) {
    std::cout << "Graph1:\n";
    print_graph(g1, get(vertex_name, g1), get(edge_name, g1));
    std::cout << "\nGraph2:\n";
    print_graph(g2, get(vertex_name, g2), get(edge_name, g2));
  }
}

}  // namespace
}  // namespace bagl
