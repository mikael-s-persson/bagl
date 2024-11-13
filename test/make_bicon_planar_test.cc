// Copyright 2007 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/adjacency_list.h"
#include "bagl/biconnected_components.h"
#include "bagl/boyer_myrvold_planar_test.h"
#include "bagl/make_biconnected_planar.h"
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
void make_line_graph(Graph& g, int size) {
  auto prev_vertex = add_vertex(g);
  for (int i = 1; i < size; ++i) {
    auto curr_vertex = add_vertex(g);
    add_edge(curr_vertex, prev_vertex, g);
    prev_vertex = curr_vertex;
  }
}

template <typename T>
class MakeBiconnPlanarTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(MakeBiconnPlanarTest);

TYPED_TEST_P(MakeBiconnPlanarTest, LineGraph) {
  using Graph = std::tuple_element_t<0, TypeParam>;
  using Edge = graph_edge_descriptor_t<Graph>;
  static constexpr std::size_t num_v = std::tuple_element_t<1, TypeParam>::value;
  Graph g;
  make_line_graph(g, num_v);
  update_vindex(g);
  update_eindex(g);

  auto embedding = vector_property_map(num_vertices(g), get(vertex_index, g), std::vector<Edge>{});
  for (auto v : vertices(g)) {
    std::ranges::copy(out_edges(v, g), std::back_inserter(embedding[v]));
  }

  EXPECT_GT(biconnected_components(g, vector_property_map(num_edges(g), get(edge_index, g), int{0}).ref()), 1);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
  make_biconnected_planar(g, embedding.ref());
  update_eindex(g);
  EXPECT_EQ(biconnected_components(g, vector_property_map(num_edges(g), get(edge_index, g), int{0}).ref()), 1);
  EXPECT_TRUE(boyer_myrvold_planarity_test(g));
}

REGISTER_TYPED_TEST_SUITE_P(MakeBiconnPlanarTest, LineGraph);

using VVGraph = adjacency_list<vec_s, vec_s, undirected_s, no_property, property<edge_index_t, std::size_t>>;
using VLGraph = adjacency_list<vec_s, list_s, undirected_s, property<vertex_index_t, std::size_t>,
                               property<edge_index_t, std::size_t>>;
using LVGraph = adjacency_list<list_s, vec_s, undirected_s, no_property, property<edge_index_t, std::size_t>>;
using LLGraph = adjacency_list<list_s, list_s, undirected_s, property<vertex_index_t, std::size_t>,
                               property<edge_index_t, std::size_t>>;
using SLGraph = adjacency_list<set_s, list_s, undirected_s, property<vertex_index_t, std::size_t>,
                               property<edge_index_t, std::size_t>>;

template <std::size_t I>
using IntC = std::integral_constant<std::size_t, I>;

using alltests =
    ::testing::Types<std::tuple<VVGraph, IntC<10>>, std::tuple<VVGraph, IntC<50>>, std::tuple<VLGraph, IntC<3>>,
                     std::tuple<VLGraph, IntC<30>>, std::tuple<LVGraph, IntC<15>>, std::tuple<LVGraph, IntC<45>>,
                     std::tuple<LLGraph, IntC<8>>, std::tuple<LLGraph, IntC<19>>, std::tuple<SLGraph, IntC<13>>,
                     std::tuple<SLGraph, IntC<20>>>;

INSTANTIATE_TYPED_TEST_SUITE_P(AllMakeBiconnPlanarTest, MakeBiconnPlanarTest, alltests);

}  // namespace
}  // namespace bagl
