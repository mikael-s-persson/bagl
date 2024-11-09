// Copyright 2007 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/make_connected.h"

#include <type_traits>

#include "bagl/adjacency_list.h"
#include "bagl/connected_components.h"
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
void make_disconnected_cycles(Graph& g, int num_cycles, int cycle_size) {
  // This graph will consist of num_cycles cycles, each of which
  // has cycle_size vertices and edges. The entire graph has
  // num_cycles * cycle_size vertices and edges, and requires
  // num_cycles - 1 edges to make it connected
  for (int i = 0; i < num_cycles; ++i) {
    auto first_vertex = add_vertex(g);
    auto curr_vertex = first_vertex;
    for (int j = 1; j < cycle_size; ++j) {
      auto prev_vertex = curr_vertex;
      curr_vertex = add_vertex(g);
      add_edge(prev_vertex, curr_vertex, g);
    }
    add_edge(curr_vertex, first_vertex, g);
  }
}

template <typename T>
class MakeConnectedTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(MakeConnectedTest);

TYPED_TEST_P(MakeConnectedTest, Basic) {
  using Graph = std::tuple_element_t<0, TypeParam>;
  static constexpr std::size_t num_cycles = std::tuple_element_t<1, TypeParam>::value;
  static constexpr std::size_t cycle_size = std::tuple_element_t<2, TypeParam>::value;

  Graph g;
  make_disconnected_cycles(g, num_cycles, cycle_size);
  update_eindex(g);
  update_vindex(g);
  auto components = vector_property_map(num_vertices(g), get(vertex_index, g), 0);
  EXPECT_EQ(connected_components(g, components.ref()), static_cast<int>(num_cycles));
  make_connected(g);
  EXPECT_EQ(connected_components(g, components.ref()), 1);
  EXPECT_EQ(num_edges(g), num_cycles * cycle_size + num_cycles - 1);

  // Now make sure that no edges are added to an already connected graph
  // when you call make_connected again

  auto old_num_edges = num_edges(g);
  make_connected(g);
  EXPECT_EQ(num_edges(g), old_num_edges);
}

REGISTER_TYPED_TEST_SUITE_P(MakeConnectedTest, Basic);

using VVGraph = adjacency_list<vec_s, vec_s, undirected_s, no_property, property<edge_index_t, std::size_t>>;
using VLGraph = adjacency_list<vec_s, list_s, undirected_s, property<vertex_index_t, std::size_t>,
                               property<edge_index_t, std::size_t>>;
using LVGraph = adjacency_list<list_s, vec_s, undirected_s, no_property, property<edge_index_t, std::size_t>>;
using LLGraph = adjacency_list<list_s, list_s, undirected_s, property<vertex_index_t, std::size_t>,
                               property<edge_index_t, std::size_t>>;

template <std::size_t I>
using IntC = std::integral_constant<std::size_t, I>;

using alltests = ::testing::Types<std::tuple<VVGraph, IntC<10>, IntC<10>>, std::tuple<LVGraph, IntC<20>, IntC<20>>,
                                  std::tuple<VLGraph, IntC<30>, IntC<30>>, std::tuple<LLGraph, IntC<40>, IntC<40>>>;

INSTANTIATE_TYPED_TEST_SUITE_P(AllMakeConnectedTest, MakeConnectedTest, alltests);

}  // namespace
}  // namespace bagl
