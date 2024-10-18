// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <cstddef>
#include <list>
#include <vector>

#include "bagl/graph_archetypes.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/vector_as_graph.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(VectorAsGraphTest, Concepts) {
  using Vector = std::vector<std::list<int> >;
  static_assert(concepts::VertexListGraph<vector_as_graph<Vector> >);
  static_assert(concepts::IncidenceGraph<vector_as_graph<Vector> >);
  static_assert(concepts::AdjacencyGraph<vector_as_graph<Vector> >);
}

}  // namespace
}  // namespace bagl
