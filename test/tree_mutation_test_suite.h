// Copyright 2024 Mikael Persson

#ifndef BAGL_TEST_TREE_MUTATION_TEST_SUITE_H_
#define BAGL_TEST_TREE_MUTATION_TEST_SUITE_H_

#include <iostream>
#include <iterator>

#include "bagl/adjacency_list.h"
#include "bagl/graph_concepts.h"
#include "bagl/linked_tree.h"
#include "bagl/property.h"
#include "bagl/tree_adaptor.h"
#include "bagl/tree_traits.h"
#include "gtest/gtest.h"

namespace bagl {

// Create special test id properties.
struct vertex_test_id_t {
  using kind = vertex_property_tag;
  static constexpr std::string_view name = "vertex_test_id";
};
constexpr vertex_test_id_t vertex_test_id = {};
struct edge_test_id_t {
  using kind = edge_property_tag;
  static constexpr std::string_view name = "edge_test_id";
};
constexpr edge_test_id_t edge_test_id = {};

struct vertex_test_bundle_t {
  std::string str_value;
};

struct edge_test_bundle_t {
  std::string str_value;
};

// Use the following vertex and edge properties in the graphs to be tested.
using test_vertex_property = property<vertex_test_id_t, std::size_t, vertex_test_bundle_t>;
using test_edge_property = property<edge_test_id_t, std::size_t, edge_test_bundle_t>;

template <concepts::VertexListGraph TreeType>
void do_final_vertex_check(const TreeType& g) {
  auto vid = get(vertex_test_id, g);
  std::vector<std::size_t> all_vertices;
  for (auto v : vertices(g)) {
    all_vertices.push_back(vid[v]);
  }
  EXPECT_THAT(all_vertices, ::testing::UnorderedElementsAreArray({1, 2, 4, 5, 6, 7, 8, 9}));
}

template <concepts::Graph TreeType>
void do_final_vertex_check(const TreeType& /*unused*/) {}

template <concepts::VertexListGraph TreeType>
void check_tree_vertex_count(const TreeType& g, std::size_t expected_count) {
  EXPECT_EQ(num_vertices(g), expected_count);
}

template <concepts::Graph TreeType>
void check_tree_vertex_count(const TreeType& /*unused*/, std::size_t /*unused*/) {}

template <concepts::BidirectionalGraph TreeType>
void do_in_edge_check(const TreeType& g, graph_vertex_descriptor_t<TreeType> v, std::size_t parent_value,
                      std::size_t cur_value) {
  using Vertex = graph_vertex_descriptor_t<TreeType>;
  auto vid = get(vertex_test_id, g);
  auto eid = get(edge_test_id, g);

  Vertex u = Vertex();
  EXPECT_NO_THROW(u = parent(v, g));
  EXPECT_EQ(vid[u], parent_value);

  for (auto e : in_edges(v, g)) {
    EXPECT_EQ(eid[e], parent_value * 1000 + cur_value);
    EXPECT_EQ(vid[source(e, g)], parent_value);
    EXPECT_EQ(vid[target(e, g)], cur_value);
  }
}

template <concepts::Graph TreeType>
void do_in_edge_check(const TreeType& /*unused*/, graph_vertex_descriptor_t<TreeType> /*unused*/,
                      std::size_t /*unused*/, std::size_t /*unused*/) {}

template <concepts::AdjacencyMatrix TreeType>
void do_edge_check(const TreeType& g, graph_vertex_descriptor_t<TreeType> u, graph_vertex_descriptor_t<TreeType> v,
                   std::size_t u_value, std::size_t v_value) {
  using Edge = graph_edge_descriptor_t<TreeType>;
  auto vid = get(vertex_test_id, g);
  auto eid = get(edge_test_id, g);

  Edge e = Edge();
  bool e_found = false;
  EXPECT_NO_THROW(std::tie(e, e_found) = edge(u, v, g));
  EXPECT_TRUE(e_found);
  EXPECT_EQ(eid[e], u_value * 1000 + v_value);
  EXPECT_EQ(vid[source(e, g)], u_value);
  EXPECT_EQ(vid[target(e, g)], v_value);
}

template <concepts::Graph TreeType>
void do_edge_check(const TreeType& /*unused*/, graph_vertex_descriptor_t<TreeType> /*unused*/,
                   graph_vertex_descriptor_t<TreeType> /*unused*/, std::size_t /*unused*/, std::size_t /*unused*/) {}

template <typename TreeType>
void check_fullbranch_integrity(const TreeType& g, graph_vertex_descriptor_t<TreeType> u) {
  auto vid = get(vertex_test_id, g);
  auto eid = get(edge_test_id, g);

  if (out_degree(u, g) == 0) {
    return;
  }

  EXPECT_EQ(out_degree(u, g), 4);
  for (auto e : out_edges(u, g)) {
    EXPECT_EQ(eid[e], vid[source(e, g)] * 1000 + vid[target(e, g)]);
    check_fullbranch_integrity(g, target(e, g));
  }
}

template <typename T>
class TreeMutationTest : public ::testing::Test {
  static_assert(concepts::IncidenceGraph<T>);
};
TYPED_TEST_SUITE_P(TreeMutationTest);

TYPED_TEST_P(TreeMutationTest, AllCases) {
  using TreeType = TypeParam;
  using Vertex = graph_vertex_descriptor_t<TreeType>;
  using Edge = graph_edge_descriptor_t<TreeType>;

  TreeType g;
  auto vid = get(vertex_test_id, g);
  auto eid = get(edge_test_id, g);

  /* MutableTree, Graph */
  Vertex v_root = graph_traits<TreeType>::null_vertex();
  EXPECT_NO_THROW(v_root = create_root(g));
  check_tree_vertex_count(g, 1);

  /* MutableTree */
  EXPECT_NO_THROW(remove_branch(v_root, g));
  check_tree_vertex_count(g, 0);

  /* MutableTree */
  std::size_t vp_r = 1;
  EXPECT_NO_THROW(v_root = create_root(g, vp_r));
  check_tree_vertex_count(g, 1);
  EXPECT_EQ(vid[v_root], 1);

  /* MutablePropertyTree */
  std::vector<test_vertex_property> props;
  EXPECT_NO_THROW(remove_branch(v_root, g, back_inserter(props)));
  check_tree_vertex_count(g, 0);
  EXPECT_EQ(props.size(), 1);
  EXPECT_EQ(get_property_value(props[0], vertex_test_id), 1);
  props.clear();

  /* MutablePropertyTree */
  EXPECT_NO_THROW(v_root = create_root(g, std::move(vp_r)));
  check_tree_vertex_count(g, 1);
  EXPECT_EQ(vid[v_root], 1);

  /* MutablePropertyTree */
  std::array<int, 4> vp_rc = {2, 3, 4, 5};
  std::array<int, 4> ep_rc = {1002, 1003, 1004, 1005};
  std::array<Vertex, 4> v_rc;
  std::array<Edge, 4> e_rc;
  for (int i = 0; i < 4; ++i) {
    bool c_added = false;
    EXPECT_NO_THROW(std::tie(v_rc[i], e_rc[i], c_added) = add_child(v_root, g, vp_rc[i], ep_rc[i]));
    EXPECT_TRUE(c_added);
  }
  check_tree_vertex_count(g, 5);

  /* Tree(IncidenceGraph) */
  EXPECT_EQ(out_degree(v_root, g), 4);

  /* MutablePropertyTree */
  std::array<std::size_t, 4> vp_rc1c = {6, 7, 8, 9};
  std::array<std::size_t, 4> ep_rc1c = {2006, 2007, 2008, 2009};
  std::array<Vertex, 4> v_rc1c;
  std::array<Edge, 4> e_rc1c;
  for (std::size_t i = 0; i < 4; ++i) {
    bool c_added = false;
    EXPECT_NO_THROW(std::tie(v_rc1c[i], e_rc1c[i], c_added) = add_child(v_rc[0], g, vp_rc1c[i], ep_rc1c[i]));
    EXPECT_TRUE(c_added);
  }
  check_tree_vertex_count(g, 9);

  /* Tree */
  EXPECT_NO_THROW(v_root = tree_root(g));
  EXPECT_EQ(vid[v_root], 1);
  {
    /* Tree(IncidenceGraph) */
    std::vector<std::size_t> e_list;
    for (auto e : out_edges(v_root, g)) {
      EXPECT_EQ(eid[e], (vid[source(e, g)] * 1000 + vid[target(e, g)]));
      e_list.push_back(eid[e]);
    }
    EXPECT_THAT(e_list, ::testing::UnorderedElementsAreArray({1002, 1003, 1004, 1005}));

    /* Tree */
    std::vector<std::size_t> vp_list;
    for (auto v : children(v_root, g)) {
      vp_list.push_back(vid[v]);
    }
    EXPECT_THAT(vp_list, ::testing::UnorderedElementsAreArray({2, 3, 4, 5}));

    /* Tree */
    std::vector<Vertex> v_list;
    for (auto v : children(v_root, g)) {
      if (vid[v] == 2) {
        /* Tree(IncidenceGraph) */
        std::vector<std::size_t> e_list2;
        for (auto e : out_edges(v, g)) {
          EXPECT_EQ(eid[e], (vid[source(e, g)] * 1000 + vid[target(e, g)]));
          e_list2.push_back(eid[e]);
        }
        EXPECT_THAT(e_list2, ::testing::UnorderedElementsAreArray({2006, 2007, 2008, 2009}));
      }
    }
  }

  /* MutablePropertyTree */
  std::array<std::size_t, 4> vp_rc2c = {10, 11, 12, 13};
  std::array<std::size_t, 4> ep_rc2c = {3010, 3011, 3012, 3013};
  std::array<Vertex, 4> v_rc2c;
  std::array<Edge, 4> e_rc2c;
  for (std::size_t i = 0; i < 4; ++i) {
    bool c_added = false;
    EXPECT_NO_THROW(std::tie(v_rc2c[i], e_rc2c[i], c_added) =
                        add_child(v_rc[1], g, std::move(vp_rc2c[i]), std::move(ep_rc2c[i])));
    EXPECT_TRUE(c_added);
  }
  check_tree_vertex_count(g, 13);

  {
    /* Tree(IncidenceGraph) */
    std::vector<std::size_t> e_list;
    for (auto e : out_edges(v_rc[1], g)) {
      EXPECT_EQ(eid[e], (vid[source(e, g)] * 1000 + vid[target(e, g)]));
      e_list.push_back(eid[e]);
    }
    EXPECT_THAT(e_list, ::testing::UnorderedElementsAreArray({3010, 3011, 3012, 3013}));

    /* Tree */
    std::vector<std::size_t> vp_list;
    for (auto v : children(v_rc[1], g)) {
      vp_list.push_back(vid[v]);
    }
    EXPECT_THAT(vp_list, ::testing::UnorderedElementsAreArray({10, 11, 12, 13}));
  }

  /* Copying function */
  check_fullbranch_integrity(g, v_root);
  {
    TreeType* p_g_cpy = nullptr;
    EXPECT_NO_THROW(p_g_cpy = new TreeType(g));
    check_fullbranch_integrity(*p_g_cpy, tree_root(*p_g_cpy));
    EXPECT_NO_THROW(delete p_g_cpy);
  }

  {
    TreeType g_cpy;
    EXPECT_NO_THROW(g_cpy = g);
    check_fullbranch_integrity(g_cpy, tree_root(g_cpy));
  }

  {
    TreeType* p_g_mv = nullptr;
    EXPECT_NO_THROW(p_g_mv = new TreeType(std::move(g)));
    check_fullbranch_integrity(*p_g_mv, tree_root(*p_g_mv));
    EXPECT_NO_THROW(g = std::move(*p_g_mv));
    check_fullbranch_integrity(g, tree_root(g));
    EXPECT_NO_THROW(delete p_g_mv);
    v_root = tree_root(g);
  }

  /* MutableTree */
  EXPECT_NO_THROW(remove_branch(v_rc[0], g));
  check_tree_vertex_count(g, 8);

  EXPECT_EQ(vid[v_root], 1);
  {
    /* Tree(IncidenceGraph) */
    std::vector<std::size_t> e_list;
    for (auto e : out_edges(v_root, g)) {
      EXPECT_EQ(eid[e], (vid[source(e, g)] * 1000 + vid[target(e, g)]));
      e_list.push_back(eid[e]);
    }
    EXPECT_THAT(e_list, ::testing::UnorderedElementsAreArray({1003, 1004, 1005}));

    /* Tree */
    std::vector<std::size_t> vp_list;
    for (auto v : children(v_root, g)) {
      vp_list.push_back(vid[v]);
    }
    EXPECT_THAT(vp_list, ::testing::UnorderedElementsAreArray({3, 4, 5}));

    /* Tree */
    std::vector<Vertex> v_list;
    for (auto v : children(v_root, g)) {
      if (vid[v] == 3) {
        /* Tree(IncidenceGraph) */
        std::vector<std::size_t> e_list2;
        for (auto e : out_edges(v, g)) {
          EXPECT_EQ(eid[e], (vid[source(e, g)] * 1000 + vid[target(e, g)]));
          e_list2.push_back(eid[e]);
        }
        EXPECT_THAT(e_list2, ::testing::UnorderedElementsAreArray({3010, 3011, 3012, 3013}));
      }
    }
  }

  /* MutablePropertyTree */
  EXPECT_NO_THROW(remove_branch(v_rc[1], g, std::back_inserter(props)));
  EXPECT_EQ(props.size(), 5);
  EXPECT_EQ(get_property_value(props[0], vertex_test_id), 3);  // the first vertex should be the root of the branch.
  std::vector<std::size_t> props_id;
  std::transform(props.begin(), props.end(), std::back_inserter(props_id),
                 [](const auto& p) { return get_property_value(p, vertex_test_id); });
  EXPECT_THAT(props_id, ::testing::UnorderedElementsAreArray({3, 10, 11, 12, 13}));
  check_tree_vertex_count(g, 3);

  EXPECT_EQ(vid[v_root], 1);
  {
    /* Tree(IncidenceGraph) */
    std::vector<std::size_t> e_list;
    for (auto e : out_edges(v_root, g)) {
      EXPECT_EQ(eid[e], (vid[source(e, g)] * 1000 + vid[target(e, g)]));
      e_list.push_back(eid[e]);
    }
    EXPECT_THAT(e_list, ::testing::UnorderedElementsAreArray({1004, 1005}));

    /* Tree */
    std::vector<std::size_t> vp_list;
    for (auto v : children(v_root, g)) {
      vp_list.push_back(vid[v]);
    }
    EXPECT_THAT(vp_list, ::testing::UnorderedElementsAreArray({4, 5}));
  }

  /* MutablePropertyTree */
  bool c_added = false;
  EXPECT_NO_THROW(std::tie(v_rc[0], e_rc[0], c_added) = add_child(v_root, g, vp_rc[0], ep_rc[0]));
  EXPECT_TRUE(c_added);
  for (std::size_t i = 0; i < 4; ++i) {
    EXPECT_NO_THROW(std::tie(v_rc1c[i], e_rc1c[i], c_added) = add_child(v_rc[0], g, vp_rc1c[i], ep_rc1c[i]));
    EXPECT_TRUE(c_added);
  }
  check_tree_vertex_count(g, 8);

  EXPECT_EQ(vid[v_root], 1);
  {
    /* Tree(IncidenceGraph) */
    std::vector<std::size_t> e_list;
    for (auto e : out_edges(v_root, g)) {
      EXPECT_EQ(eid[e], (vid[source(e, g)] * 1000 + vid[target(e, g)]));
      e_list.push_back(eid[e]);
    }
    EXPECT_THAT(e_list, ::testing::UnorderedElementsAreArray({1002, 1004, 1005}));

    /* Tree */
    std::vector<std::size_t> vp_list;
    for (auto v : children(v_root, g)) {
      vp_list.push_back(vid[v]);
    }
    EXPECT_THAT(vp_list, ::testing::UnorderedElementsAreArray({2, 4, 5}));

    /* Tree */
    std::vector<Vertex> v_list;
    for (auto v : children(v_root, g)) {
      if (vid[v] == 2) {
        /* Tree(IncidenceGraph) */
        std::vector<std::size_t> e_list2;
        for (auto e : out_edges(v, g)) {
          EXPECT_EQ(eid[e], (vid[source(e, g)] * 1000 + vid[target(e, g)]));
          e_list2.push_back(eid[e]);
        }
        EXPECT_THAT(e_list2, ::testing::UnorderedElementsAreArray({2006, 2007, 2008, 2009}));
      }
    }
  }

  /* VertexListGraph */
  do_final_vertex_check(g);
  /* BidirectionalGraph */
  do_in_edge_check(g, v_rc1c[1], 2, 7);
  /* AdjacencyMatrix */
  do_edge_check(g, v_rc[0], v_rc1c[2], 2, 8);
}

REGISTER_TYPED_TEST_SUITE_P(TreeMutationTest, AllCases);

}  // namespace bagl

#endif  // BAGL_TEST_TREE_MUTATION_TEST_SUITE_H_
