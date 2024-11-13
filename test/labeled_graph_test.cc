// (C) Copyright 2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/labeled_graph.h"

#include <ranges>
#include <string>

#include "bagl/directed_graph.h"
#include "bagl/graph_concepts.h"
#include "bagl/undirected_graph.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(LabeledGraphTest, Concepts) {
  // The labeled mutable graph hides the add_ and remove_ vertex functions
  // from the mutable graph concept, so VertexMutableGraph will not be
  // tested here.
  {
    using Graph = labeled_graph<directed_graph<>, unsigned>;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::EdgeMutableGraph<Graph>);
  }
  {
    using Graph = labeled_graph<undirected_graph<>, unsigned>;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::EdgeMutableGraph<Graph>);
  }
}

// Utility Functions and Types

struct Actor {
  Actor() {}
  Actor(std::string s) : name(std::move(s)) {}
  std::string name;
};

struct Movie {
  Movie() {}
  Movie(std::string s) : name(std::move(s)) {}
  std::string name;
};

template <typename Graph>
void init_graph(Graph& g) {
  for (int i = 0; i < 6; ++i) {
    add_vertex(g, i);
  }
}

template <typename Graph>
void label_graph(Graph& g) {
  int x = 0;
  for (auto v : vertices(g)) {
    label_vertex(v, x++, g);
  }
}

template <typename Graph>
void build_graph(Graph& g) {
  // This is the graph shown on the wikipedia page for Graph Theory.
  add_edge_by_label(5, 3, g);
  add_edge_by_label(3, 4, g);
  add_edge_by_label(3, 2, g);
  add_edge_by_label(4, 1, g);
  add_edge_by_label(4, 0, g);
  add_edge_by_label(2, 1, g);
  add_edge_by_label(1, 0, g);
  EXPECT_EQ(num_vertices(g), 6);
  EXPECT_EQ(num_edges(g), 7);
}

TEST(LabeledGraphTest, Norm) {
  {
    using Graph = labeled_graph<undirected_graph<>, int>;
    Graph g;
    init_graph(g);
    build_graph(g);
  }

  {
    using Graph = labeled_graph<directed_graph<>, int>;
    Graph g;
    init_graph(g);
    build_graph(g);
  }
}

TEST(LabeledGraphTest, Temp) {
  using Graph = undirected_graph<>;
  using LabGraph = labeled_graph<Graph*, int>;
  Graph g(6);
  LabGraph lg(&g);
  label_graph(lg);
  build_graph(lg);
}

TEST(LabeledGraphTest, Bacon) {
  std::string bacon("Kevin Bacon");
  std::string slater("Christian Slater");
  Movie murder("Murder in the First");
  {
    using Graph = labeled_graph<undirected_graph<Actor, Movie>, std::string>;
    Graph g;
    add_vertex(g, bacon);
    add_vertex(g, slater);
    add_edge_by_label(bacon, slater, g, murder);
  }

  {
    using Graph = labeled_graph<directed_graph<Actor, Movie>, std::string>;
    Graph g;
    add_vertex(g, bacon);
    add_vertex(g, slater);
    add_edge_by_label(bacon, slater, g, murder);
  }
}

TEST(LabeledGraphTest, RemoveLabeledVertex) {
  using Graph = labeled_graph<directed_graph<>, std::string>;

  Graph g;
  add_vertex(g, "foo");

  auto v = g.vertex("foo");
  EXPECT_NE(v, graph_traits<Graph>::null_vertex());

  remove_vertex("foo", g);

  auto v2 = g.vertex("foo");
  EXPECT_EQ(v2, graph_traits<Graph>::null_vertex());
}

TEST(LabeledGraphTest, MultipleAssocContainer) {
  using Graph = labeled_graph<adjacency_list<list_s, list_s, directed_s>, std::string, multiset_s>;

  Graph g;
  add_vertex(g, "test");
  add_vertex(g, "test");

  EXPECT_EQ(num_vertices(g), 2);

  remove_vertex("test", g);
  EXPECT_EQ(num_vertices(g), 1);

  remove_vertex("test", g);
  EXPECT_EQ(num_vertices(g), 0);
}

}  // namespace
}  // namespace bagl
