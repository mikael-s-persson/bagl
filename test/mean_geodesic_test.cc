// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <iostream>

#include "bagl/directed_graph.h"
#include "bagl/floyd_warshall_shortest.h"
#include "bagl/geodesic_distance.h"
#include "bagl/graph_traits.h"
#include "bagl/matrix_property_map.h"
#include "bagl/properties.h"
#include "bagl/single_property_map.h"
#include "bagl/undirected_graph.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

template <typename Graph>
void build_graph(Graph& g, std::vector<graph_vertex_descriptor_t<Graph> >& v) {
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

TEST(MeanGeodesicTest, UndirectedGraph) {
  using Graph = undirected_graph<>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;
  std::vector<Vertex> v;
  build_graph(g, v);

  auto centralities = vector_property_map(num_vertices(g), get(vertex_index, g), double{0});
  auto distances = matrix_property_store(num_vertices(g), get(vertex_index, g), int{0});

  auto cm = centralities.ref();
  auto dm = distances.ref();
  auto wm = single_property_map(1);

  floyd_warshall_all_pairs_shortest_paths(g, dm, wm);
  double geo1 = all_mean_geodesics(g, dm, cm);
  double geo2 = small_world_distance(g, cm);

  EXPECT_THAT(cm[v[0]], ::testing::DoubleEq(5.0 / 4.0));
  EXPECT_THAT(cm[v[1]], ::testing::DoubleEq(7.0 / 4.0));
  EXPECT_THAT(cm[v[2]], ::testing::DoubleEq(7.0 / 4.0));
  EXPECT_THAT(cm[v[3]], ::testing::DoubleEq(9.0 / 4.0));
  EXPECT_THAT(cm[v[4]], ::testing::DoubleEq(6.0 / 4.0));
  EXPECT_THAT(geo1, ::testing::DoubleEq(34.0 / 20.0));
  EXPECT_THAT(geo1, ::testing::DoubleEq(geo2));
}

TEST(MeanGeodesicTest, DirectedGraph) {
  using Graph = directed_graph<>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;
  std::vector<Vertex> v;
  build_graph(g, v);

  auto centralities = vector_property_map(num_vertices(g), get(vertex_index, g), double{0});
  auto distances = matrix_property_store(num_vertices(g), get(vertex_index, g), int{0});

  auto cm = centralities.ref();
  auto dm = distances.ref();
  auto wm = single_property_map(1);

  floyd_warshall_all_pairs_shortest_paths(g, dm, wm);
  double geo1 = all_mean_geodesics(g, dm, cm);
  double geo2 = small_world_distance(g, cm);

  double inf = numeric_values<double>::infinity();
  EXPECT_EQ(cm[v[0]], inf);
  EXPECT_EQ(cm[v[1]], inf);
  EXPECT_EQ(cm[v[2]], inf);
  EXPECT_THAT(cm[v[3]], ::testing::DoubleEq(10.0 / 4.0));
  EXPECT_EQ(cm[v[4]], inf);
  EXPECT_EQ(geo1, inf);
  EXPECT_EQ(geo1, geo2);
}

}  // namespace
}  // namespace bagl
