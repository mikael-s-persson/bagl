// Copyright 2009 Trustees of Indiana University.
// Authors: Michael Hansen, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <type_traits>

#include "bagl/graph_archetypes.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/grid_graph.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

template <typename T>
class GridGraphConceptsTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(GridGraphConceptsTest);

TYPED_TEST_P(GridGraphConceptsTest, Basic) {
  constexpr int Dims = TypeParam::value;
  using Graph = grid_graph<Dims>;
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;

  static_assert(concepts::BidirectionalGraph<Graph>);
  static_assert(concepts::VertexListGraph<Graph>);
  static_assert(concepts::EdgeListGraph<Graph>);
  static_assert(concepts::IncidenceGraph<Graph>);
  static_assert(concepts::AdjacencyGraph<Graph>);
  static_assert(concepts::AdjacencyMatrix<Graph>);
  static_assert(concepts::ReadablePropertyGraph<Graph, Vertex, vertex_index_t>);
  static_assert(concepts::ReadablePropertyGraph<Graph, Edge, edge_index_t>);
}

REGISTER_TYPED_TEST_SUITE_P(GridGraphConceptsTest, Basic);

using GridGraphDimensions = ::testing::Types<std::integral_constant<int, 1>, std::integral_constant<int, 2>,
                                             std::integral_constant<int, 3>, std::integral_constant<int, 4>>;

INSTANTIATE_TYPED_TEST_SUITE_P(GridGraphConcepts, GridGraphConceptsTest, GridGraphDimensions);

}  // namespace
}  // namespace bagl
