// Copyright 2004 The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/sequential_vertex_coloring.h"

#include <utility>

#include "bagl/adjacency_list.h"
#include "bagl/vector_property_map.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(SequentialVertexColoringTest, Basic) {
  using Graph = adjacency_list<list_s, vec_s, undirected_s>;

  using Edge = std::pair<std::size_t, std::size_t>;
  enum nodes { A, B, C, D, E, n };
  std::array edge_array = {Edge(A, C), Edge(B, B), Edge(B, D), Edge(B, E), Edge(C, B),
                           Edge(C, D), Edge(D, E), Edge(E, A), Edge(E, B)};
  Graph g(n, edge_array);

  // Test with the normal order
  auto color = vector_property_map(num_vertices(g), get(vertex_index, g), std::size_t{0});
  std::size_t num_colors = sequential_vertex_coloring(g, color.ref());
  EXPECT_EQ(num_colors, 3);
  EXPECT_EQ(get(color, A), 0);
  EXPECT_EQ(get(color, B), 0);
  EXPECT_EQ(get(color, C), 1);
  EXPECT_EQ(get(color, D), 2);
  EXPECT_EQ(get(color, E), 1);
}

}  // namespace
}  // namespace bagl
