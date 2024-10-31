// Copyright 2009 The Trustees of Indiana University.
//  Authors: Nicholas Edmonds
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <random>

#include "bagl/adjacency_list.h"
#include "bagl/erdos_renyi_generator.h"
#include "bagl/graph_traits.h"
#include "bagl/mesh_graph_generator.h"
#include "bagl/rmat_graph_generator.h"
#include "bagl/small_world_generator.h"
#include "bagl/ssca_graph_generator.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

// TODO: Currently, this is just a compilation test and basic input range test.

using Graph = adjacency_list<vec_s, vec_s, directed_s>;

constexpr std::size_t N = 100;
constexpr std::size_t M = 1000;
constexpr double p = 0.05;

constexpr double a = 0.57;
constexpr double b = 0.19;
constexpr double c = 0.19;
constexpr double d = 0.05;

// Test Erdos-Renyi generator
TEST(GeneratorTest, ErdosRenyi) {
  std::mt19937 gen{42};
  auto e_rg = erdos_renyi_range<is_undirected_graph_v<Graph>>(gen, N, p);

  auto it = e_rg.begin();
  while (it != e_rg.end()) {
    ++it;
  }
  // TODO: Test some property of the generated graph.

  EXPECT_EQ(it, e_rg.end());
}

TEST(GeneratorTest, SortedErdosRenyi) {
  std::mt19937 gen{42};
  auto e_rg = sorted_erdos_renyi_range<is_undirected_graph_v<Graph>>(gen, N, p);

  auto it = e_rg.begin();
  while (it != e_rg.end()) {
    ++it;
  }
  // TODO: Test some property of the generated graph.

  EXPECT_EQ(it, e_rg.end());
}

// Test Small World generator
TEST(GeneratorTest, SmallWorld) {
  std::mt19937 gen{42};
  auto e_rg = small_world_range(gen, N, M, p);

  auto it = e_rg.begin();
  while (it != e_rg.end()) {
    ++it;
  }
  // TODO: Test some property of the generated graph.

  EXPECT_EQ(it, e_rg.end());
}

// Test SSCA generator
TEST(GeneratorTest, SSCA) {
  std::mt19937 gen{42};
  auto e_rg = ssca_range(gen, N, 5, 0.5, 5, p);

  auto it = e_rg.begin();
  while (it != e_rg.end()) {
    ++it;
  }
  // TODO: Test some property of the generated graph.

  EXPECT_EQ(it, e_rg.end());
}

// Test Mesh generator
TEST(GeneratorTest, Mesh) {
  auto e_rg = mesh_range<is_undirected_graph_v<Graph>>(N, N);

  auto it = e_rg.begin();
  while (it != e_rg.end()) {
    ++it;
  }
  // TODO: Test some property of the generated graph.

  EXPECT_EQ(it, e_rg.end());
}

// Test R-MAT generator
TEST(GeneratorTest, RMat) {
  std::mt19937 gen{42};
  auto e_rg = rmat_range(gen, N, M, a, b, c, d);

  auto it = e_rg.begin();
  while (it != e_rg.end()) {
    ++it;
  }
  // TODO: Test some property of the generated graph.

  EXPECT_EQ(it, e_rg.end());
}

TEST(GeneratorTest, UniqueRMat) {
  std::mt19937 gen{42};
  auto e_rg = unique_rmat_range<is_undirected_graph_v<Graph>>(gen, N, M, a, b, c, d);

  auto it = e_rg.begin();
  while (it != e_rg.end()) {
    ++it;
  }
  // TODO: Test some property of the generated graph.

  EXPECT_EQ(it, e_rg.end());
}

TEST(GeneratorTest, SortedRMat) {
  std::mt19937 gen{42};
  auto e_rg = sorted_unique_rmat_range<is_undirected_graph_v<Graph>>(gen, N, M, a, b, c, d);

  auto it = e_rg.begin();
  while (it != e_rg.end()) {
    ++it;
  }
  // TODO: Test some property of the generated graph.

  EXPECT_EQ(it, e_rg.end());
}

TEST(GeneratorTest, PermutedRMat) {
  std::mt19937 gen{42};
  auto e_rg = sorted_unique_rmat_range<is_undirected_graph_v<Graph>>(gen, N, M, a, b, c, d, true);

  auto it = e_rg.begin();
  while (it != e_rg.end()) {
    ++it;
  }
  // TODO: Test some property of the generated graph.

  EXPECT_EQ(it, e_rg.end());
}

}  // namespace
}  // namespace bagl
