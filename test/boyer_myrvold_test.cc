// Copyright 2007 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/adjacency_list.h"
#include "bagl/boyer_myrvold_planar_test.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

struct VertexIndexUpdater {
  template <typename Graph>
  static void reset(Graph& g) {
    auto index = get(vertex_index, g);
    std::size_t cnt = 0;
    for (auto v : vertices(g)) {
      put(index, v, cnt++);
    }
  }
};

struct NoVertexIndexUpdater {
  template <typename Graph>
  static void reset(Graph&) {}
};

template <typename T>
class BoyerMyrvoldTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(BoyerMyrvoldTest);

TYPED_TEST_P(BoyerMyrvoldTest, K_5) {
  using Graph = std::tuple_element_t<0, TypeParam>;
  using VertexIndexUpdater = std::tuple_element_t<1, TypeParam>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;
  Vertex v1 = add_vertex(g);
  Vertex v2 = add_vertex(g);
  Vertex v3 = add_vertex(g);
  Vertex v4 = add_vertex(g);
  Vertex v5 = add_vertex(g);
  VertexIndexUpdater::reset(g);

  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v1, v2, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v1, v3, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v1, v4, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v1, v5, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v2, v3, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v2, v4, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v2, v5, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v3, v4, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v3, v5, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));

  // This edge should make the graph non-planar
  add_edge(v4, v5, g);
  EXPECT_FALSE(boyer_myrvold_planarity_test(g));
}

TYPED_TEST_P(BoyerMyrvoldTest, K_3_3) {
  using Graph = std::tuple_element_t<0, TypeParam>;
  using VertexIndexUpdater = std::tuple_element_t<1, TypeParam>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;
  Vertex v1 = add_vertex(g);
  Vertex v2 = add_vertex(g);
  Vertex v3 = add_vertex(g);
  Vertex v4 = add_vertex(g);
  Vertex v5 = add_vertex(g);
  Vertex v6 = add_vertex(g);
  VertexIndexUpdater::reset(g);

  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v1, v4, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v1, v5, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v1, v6, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v2, v4, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v2, v5, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v2, v6, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v3, v4, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  add_edge(v3, v5, g);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));

  // This edge should make the graph non-planar
  add_edge(v3, v6, g);
  EXPECT_FALSE(boyer_myrvold_planarity_test(g));
}

// This test creates a maximal planar graph on num_vertices vertices,
// then, if num_vertices is at least 5, adds an additional edge to
// create a non-planar graph.

TYPED_TEST_P(BoyerMyrvoldTest, MaximalPlanar) {
  using Graph = std::tuple_element_t<0, TypeParam>;
  using VertexIndexUpdater = std::tuple_element_t<1, TypeParam>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  for (std::size_t num_v : {3, 6, 10, 20, 50}) {
    Graph g;
    std::vector<Vertex> vmap;
    for (std::size_t i = 0; i < num_v; ++i) {
      vmap.push_back(add_vertex(g));
    }

    VertexIndexUpdater::reset(g);

    EXPECT_TRUE(boyer_myrvold_planarity_test(g));
    // create a cycle
    for (std::size_t i = 0; i < num_v; ++i) {
      add_edge(vmap[i], vmap[(i + 1) % num_v], g);
      EXPECT_TRUE(boyer_myrvold_planarity_test(g));
    }

    // triangulate the interior of the cycle.
    for (std::size_t i = 2; i < num_v - 1; ++i) {
      add_edge(vmap[0], vmap[i], g);
      EXPECT_TRUE(boyer_myrvold_planarity_test(g));
    }

    // triangulate the exterior of the cycle.
    for (std::size_t i = 3; i < num_v; ++i) {
      add_edge(vmap[1], vmap[i], g);
      EXPECT_TRUE(boyer_myrvold_planarity_test(g));
    }

    // Now add an additional edge, forcing the graph to be non-planar.
    if (num_v > 4) {
      add_edge(vmap[2], vmap[4], g);
      EXPECT_FALSE(boyer_myrvold_planarity_test(g));
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(BoyerMyrvoldTest, K_5, K_3_3, MaximalPlanar);

using VPropIndex = property<vertex_index_t, std::size_t>;
using TestTypes =
    ::testing::Types<std::tuple<adjacency_list<vec_s, vec_s, undirected_s>, NoVertexIndexUpdater>,
                     std::tuple<adjacency_list<vec_s, list_s, undirected_s, VPropIndex>, VertexIndexUpdater>,
                     std::tuple<adjacency_list<list_s, vec_s, undirected_s>, NoVertexIndexUpdater>,
                     std::tuple<adjacency_list<list_s, list_s, undirected_s, VPropIndex>, VertexIndexUpdater> >;
INSTANTIATE_TYPED_TEST_SUITE_P(BoyerMyrvoldOnAdjList, BoyerMyrvoldTest, TestTypes);

}  // namespace
}  // namespace bagl
