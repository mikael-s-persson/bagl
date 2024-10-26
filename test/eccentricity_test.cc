// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/eccentricity.h"

#include <iostream>

#include "bagl/directed_graph.h"
#include "bagl/floyd_warshall_shortest.h"
#include "bagl/graph_traits.h"
#include "bagl/matrix_property_map.h"
#include "bagl/single_property_map.h"
#include "bagl/undirected_graph.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

template <typename Graph>
void build_graph(Graph& g, std::vector<graph_vertex_descriptor_t<Graph>>& v) {
  // add vertices
  v.resize(5);
  for (auto& u : v) {
    u = add_vertex(g);
  }

  // add edges
  add_edge(v[0], v[1], g);
  add_edge(v[1], v[2], g);
  add_edge(v[2], v[0], g);
  add_edge(v[3], v[4], g);
  add_edge(v[4], v[0], g);
}

TEST(EccentricityTest, Undirected) {
  using Graph = undirected_graph<>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;
  std::vector<Vertex> v;
  build_graph(g, v);

  auto em = vector_property_map(num_vertices(g), get(vertex_index, g), int{0});
  auto dm = matrix_property_store(num_vertices(g), get(vertex_index, g), int{0});
  auto wm = single_property_map(1);

  floyd_warshall_all_pairs_shortest_paths(g, dm.ref(), wm);
  all_eccentricities(g, dm.ref(), em.ref());
  int rad = radius(g, em.ref());
  int dia = diameter(g, em.ref());

  EXPECT_EQ(em[v[0]], 2);
  EXPECT_EQ(em[v[1]], 3);
  EXPECT_EQ(em[v[2]], 3);
  EXPECT_EQ(em[v[3]], 3);
  EXPECT_EQ(em[v[4]], 2);
  EXPECT_EQ(rad, 2);
  EXPECT_EQ(dia, 3);
}

TEST(EccentricityTest, Directed) {
  using Graph = directed_graph<>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;
  std::vector<Vertex> v;
  build_graph(g, v);

  auto em = vector_property_map(num_vertices(g), get(vertex_index, g), int{0});
  auto dm = matrix_property_store(num_vertices(g), get(vertex_index, g), int{0});
  auto wm = single_property_map(1);

  floyd_warshall_all_pairs_shortest_paths(g, dm.ref(), wm);
  all_eccentricities(g, dm.ref(), em.ref());
  int rad = radius(g, em.ref());
  int dia = diameter(g, em.ref());

  int inf = numeric_values<int>::infinity();
  EXPECT_EQ(em[v[0]], inf);
  EXPECT_EQ(em[v[1]], inf);
  EXPECT_EQ(em[v[2]], inf);
  EXPECT_EQ(em[v[3]], 4);
  EXPECT_EQ(em[v[4]], inf);
  EXPECT_EQ(rad, 4);
  EXPECT_EQ(dia, inf);
}

}  // namespace
}  // namespace bagl
