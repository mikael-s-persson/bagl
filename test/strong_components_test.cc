// Copyright 2014 Alexander Lauser.
// Authors: Alexander Lauser
// Copyright 2024 Mikael Persson - Modernized to C++20

// This test is an adapted version of the MWE for Bug #10231 (showing
// incorrect root_map computation).

/* Output should be:
The example graph:
a --> b
b --> a c
c --> b

Vertex a is in component 0 and has root 0
Vertex b is in component 0 and has root 0
Vertex c is in component 0 and has root 0
*/

#include "bagl/strong_components.h"

#include <iostream>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/vector_property_map.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(StrongComponentsTest, Basic) {
  using Graph = adjacency_list<vec_s, vec_s, directed_s>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;

  Vertex a = add_vertex(g);
  Vertex b = add_vertex(g);
  Vertex c = add_vertex(g);

  add_edge(a, b, g);
  add_edge(b, a, g);

  add_edge(c, b, g);
  add_edge(b, c, g);

  auto vindex = get(vertex_index, g);
  auto component = vector_property_map(num_vertices(g), vindex, int{0});
  auto discover_time = vector_property_map(num_vertices(g), vindex, int{0});
  auto root = vector_property_map(num_vertices(g), vindex, graph_traits<Graph>::null_vertex());
  strong_components(g, component.ref(), root.ref(), discover_time.ref(), vindex);

  for (auto v : vertices(g)) {
    EXPECT_EQ(component[v], 0) << "Vertex: " << vindex[v];
    EXPECT_EQ(root[v], 0) << "Vertex: " << vindex[v];
  }

  if (::testing::Test::HasFailure()) {
    for (auto v : vertices(g)) {
      std::cout << "Vertex " << vindex[v] << " is in component " << component[v] << " and has root " << root[v]
                << std::endl;
    }
  }
}

}  // namespace
}  // namespace bagl
