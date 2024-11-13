// Copyright 2004-5 Trustees of Indiana University
// Author: Ronald Garcia
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/graphviz.h"

#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <utility>

#include "bagl/adjacency_list.h"
#include "bagl/compressed_sparse_row_graph.h"
#include "bagl/graph_traits.h"
#include "bagl/more_property_maps.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

using namespace std::string_literals;

using NodeId = std::string;
using EdgeId = std::pair<NodeId, NodeId>;

using Mass = float;
using Weight = double;
using Masses = std::map<NodeId, Mass>;
using Weights = std::map<EdgeId, Weight>;

struct Fixture {
  std::string graphviz_text;
  std::size_t correct_num_vertices;
  Masses masses;
  Weights weights;
};

Fixture sample_directed_basic() {
  return {
      "digraph { a  node [mass = 7.7] c e [mass = 6.66] }",
      3,
      Masses({{"a"s, 0.0f}, {"c"s, 7.7f}, {"e"s, 6.66f}}),
      Weights(),
  };
}

Fixture sample_directed_basic_aliased() {
  return {
      "digraph { a  node [mass = 7.7] \"a\" e [mass = 6.66] }",
      2,
      Masses({{"a"s, 0.0f}, {"e"s, 6.66f}}),
      Weights(),
  };
}

Fixture sample_directed_full() {
  return {
      "digraph { a -> b eDge [weight = 7.7] "
      "c -> d e-> f [weight = 6.66] "
      "d ->e->a [weight=.5]}",
      6,
      Masses(),
      Weights({{EdgeId("a", "b"), 0.0},
               {EdgeId("c", "d"), 7.7},
               {EdgeId("e", "f"), 6.66},
               {EdgeId("d", "e"), 0.5},
               {EdgeId("e", "a"), 0.5}}),
  };
}

Fixture sample_undirected_basic() {
  return {
      "graph { a  nodE [mass = 7.7] c e [mass =\\\n6.66] }",
      3,
      Masses({{"a"s, 0.0f}, {"c"s, 7.7f}, {"e"s, 6.66f}}),
      Weights(),
  };
}

Fixture sample_undirected_full() {
  return {
      "graph { a -- b eDge [weight = 7.7] "
      "c -- d e -- f [weight = 6.66] }",
      6,
      Masses(),
      Weights({{EdgeId("a", "b"), 0.0}, {EdgeId("c", "d"), 7.7}, {EdgeId("e", "f"), 6.66}}),
  };
}

template <typename Graph, typename NameMap, typename MassMap, typename WeightMap>
void test_graph(std::string const& text, Graph& g, std::size_t correct_num_vertices, Masses const& expected_masses,
                Weights const& expected_weights, std::string const& node_id, std::string const& g_name,
                NameMap name_map, MassMap mass_map, WeightMap weight_map) {
  // Construct a graph and set up the dynamic_property_maps.
  dynamic_properties dp(ignore_other_properties);
  dp.property<graph_vertex_descriptor_t<Graph>>(node_id, name_map);
  dp.property<graph_vertex_descriptor_t<Graph>>("mass", mass_map);
  dp.property<graph_edge_descriptor_t<Graph>>("weight", weight_map);
  dp.property("name", graph_name, g);

  std::istringstream is(text);
  EXPECT_TRUE(read_graphviz(is, g, dp, node_id));
  // check correct vertex count
  EXPECT_EQ(num_vertices(g), correct_num_vertices);
  // check masses
  if (!expected_masses.empty()) {
    // assume that all the masses have been set
    // for each vertex:
    for (auto i : vertices(g)) {
      std::string node_name = get(name_map, i);
      Mass node_mass = get(mass_map, i);
      EXPECT_NE(expected_masses.find(node_name), expected_masses.end());
      Mass ref_mass = expected_masses.find(node_name)->second;
      //  - compare the mass to the result in the table
      EXPECT_THAT(node_mass, ::testing::FloatNear(ref_mass, 0.0001f));
    }
  }
  // check weights
  if (!expected_weights.empty()) {
    // assume that all weights have been set
    /// for each edge:
    for (auto e : edges(g)) {
      // - get its name
      auto edge_name = std::pair(get(name_map, source(e, g)), get(name_map, target(e, g)));
      // - get its weight
      Weight edge_weight = get(weight_map, e);
      EXPECT_NE(expected_weights.find(edge_name), expected_weights.end());
      Weight ref_weight = expected_weights.find(edge_name)->second;
      // - compare the weight to the result in the table
      EXPECT_THAT(edge_weight, ::testing::DoubleNear(ref_weight, 0.000001))
          << "For edge: " << edge_name.first << " -- " << edge_name.second;
    }
  }
  if (!g_name.empty()) {
    std::string parsed_name = get_property(g, graph_name);
    EXPECT_EQ(parsed_name, g_name);
  }
}

template <typename Graph, typename NameMap, typename MassMap, typename WeightMap>
void test_graph(Fixture const& sample, Graph& g, std::string const& node_id, std::string const& g_name,
                NameMap name_map, MassMap mass_map, WeightMap weight_map) {
  test_graph(sample.graphviz_text, g, sample.correct_num_vertices, sample.masses, sample.weights, node_id, g_name,
             name_map, mass_map, weight_map);
}

template <typename Graph>
void test_graph(std::string const& dottext, std::size_t correct_num_vertices, Masses const& masses,
                Weights const& weights, std::string const& node_id = "node_id", std::string const& g_name = "") {
  Graph g;
  test_graph(dottext, g, correct_num_vertices, masses, weights, node_id, g_name, get(vertex_name, g),
             get(vertex_color, g), get(edge_weight, g));
}

template <typename Graph>
void test_graph(Fixture const& sample, std::string const& node_id = "node_id", std::string const& g_name = "") {
  Graph g;
  test_graph(sample.graphviz_text, g, sample.correct_num_vertices, sample.masses, sample.weights, node_id, g_name,
             get(vertex_name, g), get(vertex_color, g), get(edge_weight, g));
}

namespace Models {
using VertexProp = property<vertex_name_t, std::string, property<vertex_color_t, Mass>>;
using EdgeProp = property<edge_weight_t, Weight>;
using GraphProp = property<graph_name_t, std::string>;

using Graph = adjacency_list<vec_s, vec_s, undirected_s, VertexProp, EdgeProp, GraphProp>;
using DiGraph = adjacency_list<vec_s, vec_s, directed_s, VertexProp, EdgeProp, GraphProp>;
using DiGraphNoParallel = adjacency_list<set_s, vec_s, directed_s, VertexProp, EdgeProp, GraphProp>;

struct VertexBundle {
  std::string name;
  Mass mass;
};
struct EdgeBundle {
  Weight weight;
};
using CSRBundledGraph = compressed_sparse_row_graph<directed_s, VertexBundle, EdgeBundle, GraphProp>;
using CSRGraph = compressed_sparse_row_graph<directed_s, no_property, EdgeBundle, GraphProp>;
}  // namespace Models

// Basic directed graph tests
TEST(GraphVizTest, BasicDirectedGraph1) { test_graph<Models::DiGraph>(sample_directed_basic()); }

TEST(GraphVizTest, BasicDirectedGraph2) { test_graph<Models::DiGraph>(sample_directed_basic_aliased()); }

TEST(GraphVizTest, BasicDirectedGraph3) { test_graph<Models::DiGraph>(sample_directed_full()); }

// undirected graph with alternate node_id property name
TEST(GraphVizTest, UndirectedGraphAlternateNodeId) {
  test_graph<Models::Graph>(sample_undirected_basic(), "nodenames");
}

// Basic undirected graph tests
TEST(GraphVizTest, BasicUndirectedGraph1) { test_graph<Models::Graph>(sample_undirected_basic()); }

TEST(GraphVizTest, BasicUndirectedGraph2) { test_graph<Models::Graph>(sample_undirected_full()); }

// Mismatch directed graph test
TEST(GraphVizTest, MismatchDirectedGraph) {
  EXPECT_THROW(test_graph<Models::DiGraph>(sample_undirected_basic()), undirected_graph_error);
}

// Mismatch undirected graph test
TEST(GraphVizTest, MismatchUndirectedGraph) {
  EXPECT_THROW(test_graph<Models::Graph>(sample_directed_basic()), directed_graph_error);
}

// Complain about parallel edges
TEST(GraphVizTest, ParallelEdges) {
  Fixture parallel{
      "diGraph { a -> b [weight = 7.7]  a -> b [weight = 7.7] }",
      2,
      Masses(),
      Weights({{EdgeId("a", "b"), 7.7}}),
  };
  test_graph<Models::DiGraph>(parallel);
  EXPECT_THROW(test_graph<Models::DiGraphNoParallel>(parallel), bad_parallel_edge);
}

// Graph Property Test 1
TEST(GraphVizTest, GraphProperty1) {
  Fixture named{
      "digraph { graph [name=\"foo \\\"escaped\\\"\"]  a  c e [mass = 6.66] "
      "}",
      3,
      Masses({{"a"s, 0.0f}, {"c"s, 0.0f}, {"e"s, 6.66f}}),
      Weights(),
  };
  test_graph<Models::DiGraph>(named, "", "foo \"escaped\"");
}

// Graph Property Test 2
TEST(GraphVizTest, GraphProperty2) {
  Fixture named{
      "digraph { name=\"fo\"+ \"o\"  a  c e [mass = 6.66] }",
      3,
      Masses({{"a"s, 0.0f}, {"c"s, 0.0f}, {"e"s, 6.66f}}),
      Weights(),
  };
  test_graph<Models::DiGraph>(named, "", "foo");
}

// Graph Property Test 3 (HTML)
TEST(GraphVizTest, GraphProperty3) {
  std::string graph_name =
      "<html title=\"x'\" title2='y\"'>foo<b><![CDATA[><bad "
      "tag&>]]>bar</b>\n<br/>\nbaz</html>";
  Fixture html_named{
      "digraph { name=" + graph_name + "  a  c e [mass = 6.66] }",
      3,
      Masses({{"a"s, 0.0f}, {"c"s, 0.0f}, {"e"s, 6.66f}}),
      Weights(),
  };
  test_graph<Models::DiGraph>(html_named, "", graph_name);
}

// Comments embedded in strings
TEST(GraphVizTest, CommentsEmbeddedInStrings) {
  std::string gv(
      "digraph { "
      "a0 [ label = \"//depot/path/to/file_14#4\" ];"
      "a1 [ label = \"//depot/path/to/file_29#9\" ];"
      "a0 -> a1 [ color=gray ];"
      "}");
  test_graph<Models::DiGraph>(gv, 2, Masses(), Weights());
}

TEST(GraphVizTest, BasicCSRDirectedGraph) {
  auto sample = sample_directed_full();

  typedef Models::CSRBundledGraph Graph;
  Graph g;
  test_graph(sample, g, "node_id", "", get(&Models::VertexBundle::name, g), get(&Models::VertexBundle::mass, g),
             get(&Models::EdgeBundle::weight, g));
}

TEST(GraphVizTest, BasicCSRDirectedGraphExtProps) {
  auto sample = sample_directed_full();
  using Models::CSRGraph;
  CSRGraph g;
  auto vidx = get(vertex_index, g);

  auto vname = vector_property_map(0, vidx, std::string{});
  auto vmass = vector_property_map(0, vidx, Mass{});
  // External edge properties don't work due to invalidating edge index in CSR graph.
  auto eweight = get(&Models::EdgeBundle::weight, g);

  test_graph<CSRGraph>(sample, g, "node_id", "", vname.ref(), vmass.ref(), eweight);
}

TEST(GraphVizTest, Subgraphs) {
  // on the BGL side, the new parser doesn't support subgraphs
  // however, the docs promise to support reading them on the input side as
  // "syntactic sugar".
  for (auto gv : {
           Fixture{"digraph {}"},
           Fixture{"digraph { 1 -> {} }", 1},
           Fixture{"digraph { 1 -> {2} }", 2},
           Fixture{"digraph { 1; { 2; 3; } }", 3},
           Fixture{"digraph { { 2; 3; } 1; }", 3},
           Fixture{"digraph { 1; subgraph { 2; 3; } }", 3},
           Fixture{"digraph { 1 -> subgraph { 2; 3; } }", 3},
           Fixture{"digraph { 1 -> subgraph hello { 2; 3; } }", 3},
           Fixture{"digraph { 1 -> subgraph clust_Hello { 2; 3; } }", 3},
           Fixture{"digraph { 1 -> subgraph \"hello\" { 2; 3; } }", 3},
           Fixture{"digraph { {2} -> subgraph \"hello\" {{{{ 2; 3; }}}} }", 2},
       }) {
    test_graph<Models::DiGraph>(gv);
  }
}

// issue #364
TEST(GraphVizTest, SubgraphNestingLimit) {
  auto independent_nests = [=](unsigned level) {
    auto sg = std::string(level, '{') + " 2; 3; " + std::string(level, '}');
    test_graph<Models::DiGraph>({"digraph{1->" + sg + "}", 3});
    test_graph<Models::DiGraph>({"digraph{1->" + sg + ";4->" + sg + "}", 4});
  };

  constexpr unsigned limit = 255;
  independent_nests(1);
  independent_nests(limit / 2);
  independent_nests(limit - 1);
  independent_nests(limit);  // edge-case
  EXPECT_THROW(independent_nests(limit + 1), bad_graphviz_syntax);
}

}  // namespace
}  // namespace bagl
