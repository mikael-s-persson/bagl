// Copyright 2024 Mikael Persson

#include "bagl/undirected_graph.h"

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "gtest/gtest.h"
#include "test/graph_mutation_test_suite.h"
#include "test/graph_properties_test_suite.h"

namespace bagl {
namespace {

using Graph = undirected_graph<test_vertex_property, test_edge_property>;

TEST(UndirectedGraph, ConceptChecks) {
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;
  static_assert(concepts::VertexListGraph<Graph>);
  static_assert(concepts::EdgeListGraph<Graph>);
  static_assert(concepts::IncidenceGraph<Graph>);
  static_assert(concepts::AdjacencyGraph<Graph>);
  static_assert(concepts::MutableBidirectionalGraph<Graph>);
  static_assert(concepts::MutableEdgeListGraph<Graph>);
  static_assert(concepts::VertexMutablePropertyGraph<Graph>);
  static_assert(concepts::EdgeMutablePropertyGraph<Graph>);
  static_assert(concepts::LvaluePropertyGraph<Graph, Vertex, vertex_test_id_t>);
  static_assert(concepts::LvaluePropertyGraph<Graph, Edge, edge_test_id_t>);
}

INSTANTIATE_TYPED_TEST_SUITE_P(UndirectedGraphMutation, GraphMutationTest, Graph);
INSTANTIATE_TYPED_TEST_SUITE_P(UndirectedGraphTreeMutation, GraphTreeMutationTest, Graph);

using UndirPType = std::tuple<undirected_graph<VertexPropTest, EdgePropTest, GraphPropTest>, PropMapMaps>;
using UndirBType = std::tuple<undirected_graph<VertexBundleTest, EdgeBundleTest, GraphBundleTest>, BundleMaps>;

INSTANTIATE_TYPED_TEST_SUITE_P(UndirectedGraphProps, GraphPropertiesTest, UndirPType);
INSTANTIATE_TYPED_TEST_SUITE_P(UndirectedGraphBundles, GraphPropertiesTest, UndirBType);

}  // namespace
}  // namespace bagl
