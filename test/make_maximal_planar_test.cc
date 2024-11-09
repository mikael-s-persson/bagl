// Copyright 2007 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/make_maximal_planar.h"

#include "bagl/adjacency_list.h"
#include "bagl/boyer_myrvold_planar_test.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

template <typename Graph>
void update_vindex(Graph& g) {
  if constexpr (!std::is_same_v<typename Graph::vertex_list_selector, vec_s>) {
    auto vindex = get(vertex_index, g);
    std::size_t cnt = 0;
    for (auto v : vertices(g)) {
      put(vindex, v, cnt++);
    }
  }
}

template <typename Graph>
void update_eindex(Graph& g) {
  auto eindex = get(edge_index, g);
  std::size_t cnt = 0;
  for (auto e : edges(g)) {
    put(eindex, e, cnt++);
  }
}

template <typename Graph>
void make_cycle(Graph& g, std::size_t size) {
  typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;

  vertex_t first_vertex = add_vertex(g);
  vertex_t prev_vertex = first_vertex;

  for (std::size_t i = 1; i < size; ++i) {
    vertex_t curr_vertex = add_vertex(g);
    add_edge(curr_vertex, prev_vertex, g);
    prev_vertex = curr_vertex;
  }

  add_edge(first_vertex, prev_vertex, g);
}

template <typename T>
class MakeMaximalPlanarTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(MakeMaximalPlanarTest);

TYPED_TEST_P(MakeMaximalPlanarTest, Basic) {
  using Graph = std::tuple_element_t<0, TypeParam>;
  using Edge = graph_edge_descriptor_t<Graph>;
  static constexpr std::size_t num_v = std::tuple_element_t<1, TypeParam>::value;

  Graph g;
  make_cycle(g, num_v);
  update_vindex(g);
  update_eindex(g);

  auto embedding = vector_property_map(num_vertices(g), get(vertex_index, g), std::vector<Edge>{});
  for (auto v : vertices(g)) {
    std::ranges::copy(out_edges(v, g), std::back_inserter(embedding[v]));
  }

  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  make_maximal_planar(g, embedding.ref());
  update_eindex(g);

  // A graph is maximal planar exactly when it's both
  // planar and has 3 * num_vertices(g) - 6 edges.
  EXPECT_EQ(num_edges(g), 3 * num_vertices(g) - 6);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
}

REGISTER_TYPED_TEST_SUITE_P(MakeMaximalPlanarTest, Basic);

using VVGraph = adjacency_list<vec_s, vec_s, undirected_s, no_property, property<edge_index_t, std::size_t>>;
using VLGraph = adjacency_list<vec_s, list_s, undirected_s, property<vertex_index_t, std::size_t>,
                               property<edge_index_t, std::size_t>>;
using LVGraph = adjacency_list<list_s, vec_s, undirected_s, no_property, property<edge_index_t, std::size_t>>;
using LLGraph = adjacency_list<list_s, list_s, undirected_s, property<vertex_index_t, std::size_t>,
                               property<edge_index_t, std::size_t>>;

using alltests = ::testing::Types<std::tuple<VVGraph, std::integral_constant<std::size_t, 10>>,
                                  std::tuple<VVGraph, std::integral_constant<std::size_t, 50>>,
                                  std::tuple<VLGraph, std::integral_constant<std::size_t, 3>>,
                                  std::tuple<VLGraph, std::integral_constant<std::size_t, 30>>,
                                  std::tuple<LVGraph, std::integral_constant<std::size_t, 15>>,
                                  std::tuple<LVGraph, std::integral_constant<std::size_t, 45>>,
                                  std::tuple<LLGraph, std::integral_constant<std::size_t, 8>>,
                                  std::tuple<LLGraph, std::integral_constant<std::size_t, 19>>>;

INSTANTIATE_TYPED_TEST_SUITE_P(AllMakeMaximalPlanarTest, MakeMaximalPlanarTest, alltests);

}  // namespace
}  // namespace bagl
