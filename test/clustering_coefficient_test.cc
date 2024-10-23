// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/clustering_coefficient.h"

#include <iostream>

#include "bagl/directed_graph.h"
#include "bagl/graph_traits.h"
#include "bagl/undirected_graph.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

// number of vertices in the graph
constexpr std::size_t N = 5;

template <typename Graph>
void build_graph(Graph& g, std::vector<graph_vertex_descriptor_t<Graph> >& v) {
  // add vertices
  for (auto& u : v) {
    u = add_vertex(g);
  }

  // add edges
  ASSERT_TRUE(add_edge(v[0], v[1], g).second);
  ASSERT_TRUE(add_edge(v[1], v[2], g).second);
  ASSERT_TRUE(add_edge(v[2], v[0], g).second);
  ASSERT_TRUE(add_edge(v[3], v[4], g).second);
  ASSERT_TRUE(add_edge(v[4], v[0], g).second);
}

// TODO: write a test for directed clustering coefficient.

TEST(ClusteringCoefficientTest, UndirectedGraph) {
  using Graph = undirected_graph<>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;
  std::vector<Vertex> v(N);
  build_graph(g, v);

  auto cm = vector_property_map(num_vertices(g), get(vertex_index, g), double{0});

  EXPECT_EQ(num_paths_through_vertex(g, v[0]), 3);
  EXPECT_EQ(num_paths_through_vertex(g, v[1]), 1);
  EXPECT_EQ(num_paths_through_vertex(g, v[2]), 1);
  EXPECT_EQ(num_paths_through_vertex(g, v[3]), 0);
  EXPECT_EQ(num_paths_through_vertex(g, v[4]), 1);

  EXPECT_EQ(num_triangles_on_vertex(g, v[0]), 1);
  EXPECT_EQ(num_triangles_on_vertex(g, v[1]), 1);
  EXPECT_EQ(num_triangles_on_vertex(g, v[2]), 1);
  EXPECT_EQ(num_triangles_on_vertex(g, v[3]), 0);
  EXPECT_EQ(num_triangles_on_vertex(g, v[4]), 0);

  // TODO: Need a FP approximation to assert here.
  EXPECT_THAT(clustering_coefficient(g, v[0]), testing::DoubleEq(double(1) / 3.0));
  EXPECT_THAT(clustering_coefficient(g, v[1]), testing::DoubleEq(1));
  EXPECT_THAT(clustering_coefficient(g, v[2]), testing::DoubleEq(1));
  EXPECT_THAT(clustering_coefficient(g, v[3]), testing::DoubleEq(0));
  EXPECT_THAT(clustering_coefficient(g, v[4]), testing::DoubleEq(0));

  all_clustering_coefficients(g, cm.ref());

  // TODO: Need a FP approximation to assert here.
  EXPECT_THAT(cm[v[0]], testing::DoubleEq(double(1) / 3.0));
  EXPECT_THAT(cm[v[1]], testing::DoubleEq(1));
  EXPECT_THAT(cm[v[2]], testing::DoubleEq(1));
  EXPECT_THAT(cm[v[3]], testing::DoubleEq(0));
  EXPECT_THAT(cm[v[4]], testing::DoubleEq(0));

  // I would have used check_close, but apparently, that requires
  // me to link this against a library - which I don't really want
  // to do. Basically, this makes sure that that coefficient is
  // within some tolerance (like 1/10 million).
  double coef = mean_clustering_coefficient(g, cm.ref());
  EXPECT_THAT(coef, ::testing::DoubleNear((7.0 / 15.0), 1e-7));
}

}  // namespace
}  // namespace bagl
