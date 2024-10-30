// Copyright (C) 2001-20044 Douglas Gregor (dgregor at cs dot indiana dot edu)
// Revision History:
// 29 Nov 2001    Jeremy Siek
//      Changed to use Boost.Random.
// 29 Nov 2001    Doug Gregor
//      Initial checkin.
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/isomorphism.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <unordered_set>

#include "bagl/adjacency_list.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/property_map.h"
#include "bagl/subgraph.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

typedef std::mt19937 random_generator_type;

template <typename Graph1, typename Graph2>
auto randomly_permute_graph(const Graph1& g1, Graph2& g2) {
  // Need a clean graph to start with
  EXPECT_EQ(num_vertices(g2), 0);
  EXPECT_EQ(num_edges(g2), 0);

  using Vertex1 = graph_vertex_descriptor_t<Graph1>;

  random_generator_type gen{42};

  // Decide new order
  std::vector<Vertex1> orig_vertices;
  std::ranges::copy(vertices(g1), std::back_inserter(orig_vertices));
  std::shuffle(orig_vertices.begin(), orig_vertices.end(), gen);
  auto vertex_map = vector_property_map(num_vertices(g1), get(vertex_index, g1), graph_traits<Graph2>::null_vertex());

  for (std::size_t i = 0; i < num_vertices(g1); ++i) {
    vertex_map[orig_vertices[i]] = add_vertex(g2);
  }

  for (auto e : edges(g1)) {
    add_edge(vertex_map[source(e, g1)], vertex_map[target(e, g1)], g2);
  }

  return vertex_map;
}

template <typename Graph>
void generate_random_digraph(Graph& g, double edge_probability) {
  random_generator_type random_gen;
  std::uniform_real_distribution<double> distrib;
  auto v_rg = vertices(g);
  for (auto u = v_rg.begin(); u != v_rg.end(); ++u) {
    auto v = u;
    ++v;
    for (; v != v_rg.end(); ++v) {
      if (distrib(random_gen) <= edge_probability) {
        add_edge(*u, *v, g);
      }
    }
  }
}

TEST(IsomorphismTest, Second) {
  using Graph1 = adjacency_list<vec_s, vec_s, bidirectional_s>;
  using Graph2 = adjacency_list<list_s, list_s, bidirectional_s, property<vertex_index_t, int> >;

  Graph1 g1(2);
  auto v1_rg = vertices(g1);
  add_edge(*v1_rg.begin(), *std::next(v1_rg.begin()), g1);
  add_edge(*std::next(v1_rg.begin()), *std::next(v1_rg.begin()), g1);
  Graph2 g2;
  randomly_permute_graph(g1, g2);

  int v_idx = 0;
  for (auto v : vertices(g2)) {
    put(vertex_index, g2, v, v_idx++);
  }

  auto vindex1 = get(vertex_index, g1);
  auto vindex2 = get(vertex_index, g2);

  auto mapping = vector_property_map(num_vertices(g1), vindex1, graph_traits<Graph2>::null_vertex());

  bool isomorphism_correct = isomorphism(g1, g2, mapping.ref(), make_degree_invariant(g1, vindex1),
                                         make_degree_invariant(g2, vindex2), 0, vindex1, vindex2);
  EXPECT_TRUE(isomorphism_correct);

  bool verify_correct = verify_isomorphism(g1, g2, mapping.ref());
  EXPECT_TRUE(verify_correct);

  if (!isomorphism_correct || !verify_correct) {
    print_graph(g1);
    print_graph(g2);
  }
}

TEST(IsomorphismTest, First) {
  using Graph1 = adjacency_list<vec_s, vec_s, bidirectional_s>;
  using Graph2 = adjacency_list<list_s, list_s, bidirectional_s, property<vertex_index_t, int> >;

  int n = 30;
  double edge_probability = 0.45;

  Graph1 g1(n);
  generate_random_digraph(g1, edge_probability);
  Graph2 g2;
  randomly_permute_graph(g1, g2);

  int v_idx = 0;
  for (auto v : vertices(g2)) {
    put(vertex_index, g2, v, v_idx++);
  }

  auto vindex1 = get(vertex_index, g1);
  auto vindex2 = get(vertex_index, g2);

  auto mapping = vector_property_map(num_vertices(g1), vindex1, graph_traits<Graph2>::null_vertex());

  bool isomorphism_correct = isomorphism(g1, g2, mapping.ref(), make_degree_invariant(g1, vindex1),
                                         make_degree_invariant(g2, vindex2), 0, vindex1, vindex2);
  EXPECT_TRUE(isomorphism_correct);

  bool verify_correct = verify_isomorphism(g1, g2, mapping.ref());
  EXPECT_TRUE(verify_correct);

  if (!isomorphism_correct || !verify_correct) {
    print_graph(g1);
    print_graph(g2);
  }
}

template <typename Graph>
struct ColorFunctor {
  ColorFunctor(const Graph& g) : g_(g) {}

  auto operator()(graph_vertex_descriptor_t<Graph> v) const { return get(vertex_color, g_, v); }

  auto max() const {
    auto max_result = std::numeric_limits<int>::min();
    for (auto v : vertices(g_)) {
      max_result = std::max(max_result, this->operator()(v));
    }
    return max_result + 1;
  }

  const Graph& g_;
};

TEST(IsomorphismTest, Colored) {
  using Graph1 = adjacency_list<vec_s, vec_s, bidirectional_s, property<vertex_color_t, int> >;
  using Graph2 =
      adjacency_list<list_s, list_s, bidirectional_s, property<vertex_index_t, int, property<vertex_color_t, int> > >;

  int n = 30;
  double edge_probability = 0.45;

  Graph1 g1(n);
  generate_random_digraph(g1, edge_probability);

  Graph2 g2;
  auto vertex_map = randomly_permute_graph(g1, g2);

  std::vector<int> colors(n);
  auto midpoint = colors.begin() + n / 2;
  std::fill(colors.begin(), midpoint, 0);
  std::fill(midpoint, colors.end(), 1);

  random_generator_type gen{42};
  std::shuffle(colors.begin(), colors.end(), gen);

  int v_idx = 0;
  for (auto v : vertices(g1)) {
    put(get(vertex_color, g1), v, colors[v_idx]);
    put(get(vertex_color, g2), vertex_map[v], colors[v_idx]);
    v_idx += 1;
  }

  v_idx = 0;
  for (auto v : vertices(g2)) {
    put(get(vertex_index_t(), g2), v, v_idx);
    v_idx += 1;
  }

  auto vindex1 = get(vertex_index, g1);
  auto vindex2 = get(vertex_index, g2);

  auto mapping = vector_property_map(num_vertices(g1), vindex1, graph_traits<Graph2>::null_vertex());

  bool isomorphism_correct =
      isomorphism(g1, g2, mapping.ref(), ColorFunctor<Graph1>(g1), ColorFunctor<Graph2>(g2), 0, vindex1, vindex2);
  EXPECT_TRUE(isomorphism_correct);

  bool verify_correct = verify_isomorphism(g1, g2, mapping.ref());
  EXPECT_TRUE(verify_correct);

  if (!isomorphism_correct || !verify_correct) {
    print_graph(g1);
    print_graph(g2);
  }

  for (auto v : vertices(g1)) {
    EXPECT_NE(mapping[v], Graph2::null_vertex()) << "At " << vindex1[v];
  }

  // Map is bijective if each vertex of the second graph occurs only once
  std::unordered_set<std::size_t> v2_set;
  for (auto v : vertices(g1)) {
    EXPECT_TRUE(v2_set.insert(vindex2[mapping[v]]).second)
        << "At " << vindex1[v] << " mapping to " << vindex2[mapping[v]];
  }
}

}  // namespace
}  // namespace bagl
