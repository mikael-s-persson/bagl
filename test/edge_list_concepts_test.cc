// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <cstddef>

#include "bagl/edge_list.h"
#include "bagl/graph_archetypes.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(EdgeListTest, Concepts) {
  using E = std::pair<int, int>;

  using EdgeList = edge_list<E*, E, std::ptrdiff_t, std::random_access_iterator_tag>;

  using Edge = graph_edge_descriptor_t<EdgeList>;

  static_assert(concepts::EdgeListGraph<EdgeList>);
  static_assert(concepts::ReadablePropertyGraph<EdgeList, Edge, edge_index_t>);
}

}  // namespace
}  // namespace bagl
