// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/degree_centrality.h"

#include <vector>

#include "bagl/directed_graph.h"
#include "bagl/graph_traits.h"
#include "bagl/undirected_graph.h"
#include "bagl/vector_property_map.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

// useful types
// number of vertices in the graph
constexpr std::size_t N = 5;

template <typename Graph>
void build_graph(Graph& g, std::vector<graph_vertex_descriptor_t<Graph> >& v) {
  // add vertices
  for (std::size_t i = 0; i < N; ++i) {
    v[i] = add_vertex(g);
  }

  // add edges
  add_edge(v[0], v[1], g);
  add_edge(v[1], v[2], g);
  add_edge(v[2], v[0], g);
  add_edge(v[3], v[4], g);
  add_edge(v[4], v[0], g);
}

TEST(DegreeCentralityTest, Undirected) {
  using Graph = undirected_graph<>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;
  std::vector<Vertex> v(N);
  build_graph(g, v);

  auto cm = vector_property_map(num_vertices(g), get(vertex_index, g), std::size_t{0});
  all_degree_centralities(g, cm.ref());

  EXPECT_EQ(cm[v[0]], 3);
  EXPECT_EQ(cm[v[1]], 2);
  EXPECT_EQ(cm[v[2]], 2);
  EXPECT_EQ(cm[v[3]], 1);
  EXPECT_EQ(cm[v[4]], 2);
}

TEST(DegreeCentralityTest, Influence) {
  using Graph = directed_graph<>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;
  std::vector<Vertex> v(N);
  build_graph(g, v);

  auto cm = vector_property_map(num_vertices(g), get(vertex_index, g), std::size_t{0});
  all_influence_values(g, cm.ref());

  EXPECT_EQ(cm[v[0]], 1);
  EXPECT_EQ(cm[v[1]], 1);
  EXPECT_EQ(cm[v[2]], 1);
  EXPECT_EQ(cm[v[3]], 1);
  EXPECT_EQ(cm[v[4]], 1);
}

TEST(DegreeCentralityTest, Prestige) {
  using Graph = directed_graph<>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;
  std::vector<Vertex> v(N);
  build_graph(g, v);

  auto cm = vector_property_map(num_vertices(g), get(vertex_index, g), std::size_t{0});
  all_prestige_values(g, cm.ref());

  EXPECT_EQ(cm[v[0]], 2);
  EXPECT_EQ(cm[v[1]], 1);
  EXPECT_EQ(cm[v[2]], 1);
  EXPECT_EQ(cm[v[3]], 0);
  EXPECT_EQ(cm[v[4]], 1);
}

}  // namespace
}  // namespace bagl
