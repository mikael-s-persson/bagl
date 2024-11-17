// Copyright 2024 Mikael Persson

#ifndef BAGL_TEST_GRAPH_PROPERTIES_TEST_SUITE_H_
#define BAGL_TEST_GRAPH_PROPERTIES_TEST_SUITE_H_

#include <iostream>
#include <string>

#include "bagl/adjacency_list.h"
#include "bagl/graph_traits.h"
#include "bagl/linked_tree.h"
#include "bagl/tree_adaptor.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {

using VertexPropTest = property<vertex_name_t, std::string, property<vertex_distance_t, int>>;

using EdgePropTest = property<edge_name_t, std::string, property<edge_weight_t, int>>;

using GraphPropTest = property<graph_name_t, std::string>;

struct PropMapMaps {
  using vname_t = vertex_name_t;
  using vdistance_t = vertex_distance_t;
  using ename_t = edge_name_t;
  using eweight_t = edge_weight_t;
  using gname_t = graph_name_t;

  static constexpr vname_t vname = vertex_name;
  static constexpr vdistance_t vdistance = vertex_distance;
  static constexpr ename_t ename = edge_name;
  static constexpr eweight_t eweight = edge_weight;
  static constexpr gname_t gname = graph_name;
};

struct VertexBundleTest {
  std::string name;
  int distance;
};

struct EdgeBundleTest {
  std::string name;
  int weight;
};

struct GraphBundleTest {
  std::string name;
};

struct BundleMaps {
  using vname_t = std::string VertexBundleTest::*;
  using vdistance_t = int VertexBundleTest::*;
  using ename_t = std::string EdgeBundleTest::*;
  using eweight_t = int EdgeBundleTest::*;
  using gname_t = std::string GraphBundleTest::*;

  static constexpr vname_t vname = &VertexBundleTest::name;
  static constexpr vdistance_t vdistance = &VertexBundleTest::distance;
  static constexpr ename_t ename = &EdgeBundleTest::name;
  static constexpr eweight_t eweight = &EdgeBundleTest::weight;
  static constexpr gname_t gname = &GraphBundleTest::name;
};

template <typename T>
class GraphPropertiesTest : public ::testing::Test {
 public:
  using Graph = std::tuple_element_t<0, T>;
  using Maps = std::tuple_element_t<1, T>;
  static_assert(concepts::IncidenceGraph<Graph>);
  static_assert(concepts::VertexListGraph<Graph>);
};
TYPED_TEST_SUITE_P(GraphPropertiesTest);

TYPED_TEST_P(GraphPropertiesTest, AllCases) {
  using Maps = std::tuple_element_t<1, TypeParam>;
  using Graph = std::tuple_element_t<0, TypeParam>;
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;

  using VNameMap = property_map_t<Graph, typename Maps::vname_t>;
  using VDistMap = property_map_t<Graph, typename Maps::vdistance_t>;
  using ENameMap = property_map_t<Graph, typename Maps::ename_t>;
  using EWeightMap = property_map_t<Graph, typename Maps::eweight_t>;

  using VNameCMap = property_map_const_t<Graph, typename Maps::vname_t>;
  using VDistCMap = property_map_const_t<Graph, typename Maps::vdistance_t>;
  using ENameCMap = property_map_const_t<Graph, typename Maps::ename_t>;
  using EWeightCMap = property_map_const_t<Graph, typename Maps::eweight_t>;

  Graph g;

  VNameMap vname_map;
  EXPECT_NO_THROW(vname_map = get(Maps::vname, g));
  VDistMap dist_map;
  EXPECT_NO_THROW(dist_map = get(Maps::vdistance, g));
  ENameMap ename_map;
  EXPECT_NO_THROW(ename_map = get(Maps::ename, g));
  EWeightMap weight_map;
  EXPECT_NO_THROW(weight_map = get(Maps::eweight, g));

  const Graph& cg = g;

  VNameCMap vname_cmap;
  EXPECT_NO_THROW(vname_cmap = get(Maps::vname, cg));
  VDistCMap dist_cmap;
  EXPECT_NO_THROW(dist_cmap = get(Maps::vdistance, cg));
  ENameCMap ename_cmap;
  EXPECT_NO_THROW(ename_cmap = get(Maps::ename, cg));
  EWeightCMap weight_cmap;
  EXPECT_NO_THROW(weight_cmap = get(Maps::eweight, cg));

  Vertex v_root = create_root(g);
  Vertex v1;
  Edge e1;
  bool c1_added;
  std::tie(v1, e1, c1_added) = add_child(v_root, g);
  EXPECT_TRUE(c1_added);

  EXPECT_NO_THROW(put(vname_map, v_root, "root-vertex"));
  EXPECT_EQ(get(vname_map, v_root), "root-vertex");
  EXPECT_EQ(get(vname_cmap, v_root), "root-vertex");

  EXPECT_NO_THROW(put(dist_map, v_root, 42));
  EXPECT_EQ(get(dist_map, v_root), 42);
  EXPECT_EQ(get(dist_cmap, v_root), 42);

  EXPECT_NO_THROW(put(ename_map, e1, "root-edge"));
  EXPECT_EQ(get(ename_map, e1), "root-edge");
  EXPECT_EQ(get(ename_cmap, e1), "root-edge");

  EXPECT_NO_THROW(put(weight_map, e1, 69));
  EXPECT_EQ(get(weight_map, e1), 69);
  EXPECT_EQ(get(weight_cmap, e1), 69);

  EXPECT_NO_THROW(vname_map[v1] = "child-vertex");
  EXPECT_EQ(vname_map[v1], "child-vertex");
  EXPECT_EQ(vname_cmap[v1], "child-vertex");
  EXPECT_EQ(get(vname_map, v_root), "root-vertex");
  EXPECT_EQ(get(vname_cmap, v_root), "root-vertex");

  EXPECT_NO_THROW(dist_map[v1] = 21);
  EXPECT_EQ(dist_map[v1], 21);
  EXPECT_EQ(dist_cmap[v1], 21);
  EXPECT_EQ(get(dist_map, v_root), 42);
  EXPECT_EQ(get(dist_cmap, v_root), 42);

  EXPECT_NO_THROW(ename_map[e1] = "root-edge-again");
  EXPECT_EQ(ename_map[e1], "root-edge-again");
  EXPECT_EQ(ename_cmap[e1], "root-edge-again");

  EXPECT_NO_THROW(weight_map[e1] = 96);
  EXPECT_EQ(weight_map[e1], 96);
  EXPECT_EQ(weight_cmap[e1], 96);
}

REGISTER_TYPED_TEST_SUITE_P(GraphPropertiesTest, AllCases);

}  // namespace bagl

#endif  // BAGL_TEST_GRAPH_PROPERTIES_TEST_SUITE_H_
