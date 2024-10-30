// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <vector>

#include "bagl/directed_graph.h"
#include "bagl/undirected_graph.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

template <typename Graph>
void test() {
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using IndexMap = property_map_t<Graph, vertex_index_t>;

  constexpr std::size_t N = 5;

  Graph g;
  IndexMap vindex = get(vertex_index, g);
  (void)vindex;

  // build up the graph
  std::vector<Vertex> v(N);
  for (auto& u : v) {
    u = add_vertex(g);
  }

  // after the first build, we should have these conditions
  EXPECT_EQ(max_vertex_index(g), N);
  for (std::size_t i = 0; i < N; ++i) {
    EXPECT_EQ(get_vertex_index(v[i], g), i);
  }

  // Remove all vertices and then re-add them.
  for (std::size_t i = 0; i < N; ++i) {
    remove_vertex(v[i], g);
  }
  EXPECT_EQ(num_vertices(g), 0);
  for (auto& u : v) {
    u = add_vertex(g);
  }
  EXPECT_EQ(num_vertices(g), N);

  // Before renumbering, our vertices should be off by N. In other words,
  // we're not in a good state.
  EXPECT_EQ(max_vertex_index(g), 10);
  for (std::size_t i = 0; i < N; ++i) {
    EXPECT_EQ(get_vertex_index(v[i], g), N + i);
  }

  // Renumber vertices
  renumber_vertex_indices(g);

  // And we should be back to the initial condition
  EXPECT_EQ(max_vertex_index(g), N);
  for (std::size_t i = 0; i < N; ++i) {
    EXPECT_EQ(get_vertex_index(v[i], g), i);
  }
}

TEST(IndexGraphTest, UndirectedGraph) { test<undirected_graph<> >(); }

TEST(IndexGraphTest, DirectedGraph) { test<directed_graph<> >(); }

// Make sure that graphs constructed over n vertices will actually number
// their vertices correctly.
template <typename Graph>
void build() {
  using IndexMap = property_map_t<Graph, vertex_index_t>;

  static const size_t N = 5;

  Graph g(N);
  EXPECT_EQ(max_vertex_index(g), N);

  IndexMap vindex = get(vertex_index, g);
  (void)vindex;

  // Each vertex should be numbered correctly.
  std::size_t x = 0;
  for (auto u : vertices(g)) {
    EXPECT_EQ(get_vertex_index(u, g), x);
    x++;
  }
}

TEST(IndexGraphBuild, UndirectedGraph) { build<undirected_graph<> >(); }

}  // namespace
}  // namespace bagl
