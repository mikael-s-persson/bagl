// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/graph_concepts.h"

#include "bagl/graph_archetypes.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(GraphConceptsTest, Archetypes) {
  // Check graph concepts againt their archetypes
  using Vertex = semiregular_archetype<1>;

  using Graph1 = incidence_graph_archetype<Vertex, directed_tag, allow_parallel_edge_tag>;
  static_assert(concepts::IncidenceGraph<Graph1>);

  using Graph2 = adjacency_graph_archetype<Vertex, directed_tag, allow_parallel_edge_tag>;
  static_assert(concepts::AdjacencyGraph<Graph2>);

  using Graph3 = vertex_list_graph_archetype<Vertex, directed_tag, allow_parallel_edge_tag>;
  static_assert(concepts::VertexListGraph<Graph3>);

  static_assert(concepts::ColorValue<color_value_archetype>);

  using Graph4 = incidence_graph_archetype<Vertex, directed_tag, allow_parallel_edge_tag>;
  using Graph5 = property_graph_archetype<Graph4, vertex_color_t, color_value_archetype>;
  static_assert(concepts::PropertyGraph<Graph5, Vertex, vertex_color_t>);
}

}  // namespace
}  // namespace bagl
