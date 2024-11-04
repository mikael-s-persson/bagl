// Copyright 2017 Felix Salfelder
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/minimum_degree_ordering.h"

#include <map>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(MinimumDegreeOrdering, Empty) {
  using Graph = adjacency_list<vec_s, vec_s, directed_s>;
  std::size_t n = 10;
  Graph g(n);

  auto id = get(vertex_index, g);
  auto degree = vector_property_map(n, id, 0);
  auto supernode_sizes = vector_property_map(n, id, 1);
  auto io = vector_property_map(n, id, 0);
  auto o = vector_property_map(n, id, 0);

  minimum_degree_ordering(g, degree.ref(), io.ref(), o.ref(), supernode_sizes.ref(), 0, id);

  for (std::size_t k = 0; k < n; ++k) {
    EXPECT_EQ(o[io[k]], k);
  }
}

}  // namespace
}  // namespace bagl
