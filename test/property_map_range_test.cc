// Copyright 2024 Mikael Persson

#include "bagl/property_map_range.h"

#include <ranges>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

template <typename Rg>
auto as_vector(Rg&& rg) {
  using Value = std::ranges::range_value_t<Rg&&>;
  return std::vector<Value>(rg.begin(), rg.end());
}

TEST(PropertyMapRangeTest, KeyRanges) {
  vector_property_map<int> pmap;
  put(pmap, 0, 0);
  put(pmap, 1, 2);
  put(pmap, 2, 4);
  put(pmap, 3, 6);
  put(pmap, 4, 8);
  EXPECT_THAT(as_vector(make_property_map_range(pmap.ref(), std::ranges::iota_view(0, 5))),
              ::testing::ElementsAreArray({0, 2, 4, 6, 8}));
}

struct TestVertBundle {
  int id;
  std::string name;
};

struct TestEdgeBundle {
  int id;
  std::string name;
};

struct TestGraphBundle {
  std::string name;
};

TEST(PropertyMapRangeTest, GraphRanges) {
  using Graph =
      adjacency_list<vec_s, vec_s, directed_s, property<vertex_distance_t, int, TestVertBundle>,
                     property<edge_weight_t, int, TestEdgeBundle>, property<graph_name_t, int, TestGraphBundle>>;

  Graph g;
  set_property(g, graph_name, 42);
  g[graph_bundle] = TestGraphBundle{"g_bundle"};

  add_vertex(g, 0, 0, "a");
  add_vertex(g, 2, 1, "b");
  add_vertex(g, 4, 2, "c");
  add_vertex(g, 6, 3, "d");
  add_vertex(g, 8, 4, "e");

  add_edge(0, 1, g, 10, 0, "ab");
  add_edge(1, 2, g, 20, 1, "bc");
  add_edge(2, 3, g, 30, 2, "cd");
  add_edge(3, 4, g, 40, 3, "de");

  EXPECT_THAT(as_vector(make_graph_property_map_range(g, vertex_distance)),
              ::testing::ElementsAreArray({0, 2, 4, 6, 8}));
  EXPECT_THAT(as_vector(make_graph_property_map_range(g, edge_weight)),
              ::testing::UnorderedElementsAreArray({10, 20, 30, 40}));
  EXPECT_THAT(as_vector(make_graph_property_map_range(g, graph_name)), ::testing::ElementsAreArray({42}));

  EXPECT_THAT(as_vector(make_property_map_range(get(&TestVertBundle::name, g), vertices(g))),
              ::testing::ElementsAreArray({"a", "b", "c", "d", "e"}));
  EXPECT_THAT(as_vector(make_property_map_range(get(&TestEdgeBundle::name, g), edges(g))),
              ::testing::UnorderedElementsAreArray({"ab", "bc", "cd", "de"}));
  EXPECT_THAT(as_vector(make_property_map_range(get(&TestGraphBundle::name, g), make_single_graph_range(g))),
              ::testing::ElementsAreArray({"g_bundle"}));
}

}  // namespace
}  // namespace bagl
