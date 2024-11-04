//  (C) Copyright Jeremy Siek 2004
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/subgraph.h"

#include <random>
#include <set>

#include "bagl/adjacency_list.h"
#include "bagl/filtered_graph.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/random.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

// UNDER CONSTRUCTION

namespace bagl {
namespace {

template <typename Graph>
struct graph_test {
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;

  struct ignore_vertex {
    ignore_vertex() {}
    ignore_vertex(Vertex v) : u(v) {}
    bool operator()(Vertex x) const { return x != u; }
    Vertex u;
  };
  struct ignore_edge {
    ignore_edge() {}
    ignore_edge(Edge f) : e(f) {}
    bool operator()(Edge x) const { return x != e; }
    Edge e;
  };
  struct ignore_edges {
    ignore_edges(Vertex u, Vertex v, const Graph& gr) : s(u), t(v), g(gr) {}
    bool operator()(Edge x) const { return !(source(x, g) == s && target(x, g) == t); }
    Vertex s;
    Vertex t;
    const Graph& g;
  };

  //=========================================================================
  // Traversal Operations

  void test_incidence_graph(const std::vector<Vertex>& vertex_set, const std::vector<Edge>& edge_set, const Graph& g) {
    for (Vertex u : vertex_set) {
      std::vector<Vertex> adj;
      for (auto e : edge_set) {
        if (source(e, g) == u) {
          adj.push_back(target(e, g));
        }
      }

      auto oe_rg = out_edges(u, g);
      EXPECT_EQ(out_degree(u, g), adj.size());
      EXPECT_EQ(std::distance(oe_rg.begin(), oe_rg.end()), out_degree(u, g));
      for (Edge e : oe_rg) {
        EXPECT_EQ(source(e, g), u);
        EXPECT_NE(std::ranges::count(adj, target(e, g)), 0);
      }
    }
  }

  void test_bidirectional_graph(const std::vector<Vertex>& vertex_set, const std::vector<Edge>& edge_set,
                                const Graph& g) {
    for (Vertex v : vertex_set) {
      std::vector<Vertex> inv_adj;
      for (auto e : edge_set) {
        if (target(e, g) == v) {
          inv_adj.push_back(source(e, g));
        }
      }

      auto oe_rg = in_edges(v, g);
      EXPECT_EQ(in_degree(v, g), inv_adj.size());
      EXPECT_EQ(std::distance(oe_rg.begin(), oe_rg.end()), in_degree(v, g));
      for (Edge e : oe_rg) {
        EXPECT_EQ(target(e, g), v);
        EXPECT_NE(std::ranges::count(inv_adj, source(e, g)), 0);
      }
    }
  }

  void test_adjacency_graph(const std::vector<Vertex>& vertex_set, const std::vector<Edge>& edge_set, const Graph& g) {
    for (Vertex u : vertex_set) {
      std::vector<Vertex> adj;
      for (auto e : edge_set) {
        if (source(e, g) == u) {
          adj.push_back(target(e, g));
        }
      }

      auto adj_rg = adjacent_vertices(u, g);
      EXPECT_EQ(std::distance(adj_rg.begin(), adj_rg.end()), adj.size());
      for (Vertex v : adj_rg) {
        EXPECT_NE(std::ranges::count(adj, v), 0);
      }
    }
  }

  void test_vertex_list_graph(const std::vector<Vertex>& vertex_set, const Graph& g) {
    auto v_rg = vertices(g);
    EXPECT_EQ(num_vertices(g), vertex_set.size());
    EXPECT_EQ(std::distance(v_rg.begin(), v_rg.end()), num_vertices(g));
    for (Vertex v : v_rg) {
      EXPECT_NE(std::ranges::count(vertex_set, v), 0);
    }
  }

  void test_edge_list_graph(const std::vector<Vertex>& vertex_set, const std::vector<Edge>& edge_set, const Graph& g) {
    auto e_rg = edges(g);
    EXPECT_EQ(num_edges(g), edge_set.size());
    EXPECT_EQ(std::distance(e_rg.begin(), e_rg.end()), num_edges(g));
    for (Edge e : e_rg) {
      EXPECT_NE(std::find_if(edge_set.begin(), edge_set.end(), connects(source(e, g), target(e, g), g)),
                edge_set.end());
      EXPECT_NE(std::ranges::count(vertex_set, source(e, g)), 0);
      EXPECT_NE(std::ranges::count(vertex_set, target(e, g)), 0);
    }
  }

  void test_adjacency_matrix(const std::vector<Vertex>& vertex_set, const std::vector<Edge>& edge_set, const Graph& g) {
    for (auto u : vertex_set) {
      for (auto v : vertex_set) {
        auto [e, e_found] = edge(u, v, g);
        if (e_found) {
          EXPECT_NE(std::find_if(edge_set.begin(), edge_set.end(), connects(source(e, g), target(e, g), g)),
                    edge_set.end());
        }
      }
    }
  }

  //=========================================================================
  // Mutating Operations

  void test_add_vertex(Graph& g) {
    Graph cpy;
    auto iso_map = vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<Graph>::null_vertex());
    copy_graph(g, cpy, iso_map.ref());

    EXPECT_TRUE(verify_isomorphism(g, cpy, iso_map.ref()));

    Vertex v = add_vertex(g);

    EXPECT_EQ(num_vertices(g), num_vertices(cpy) + 1);
    EXPECT_EQ(out_degree(v, g), 0);

    // Make sure the rest of the graph stayed the same
    EXPECT_TRUE(verify_isomorphism(make_filtered_graph(g, keep_all(), ignore_vertex(v)), cpy, iso_map.ref()));
  }

  void test_add_edge(Vertex u, Vertex v, Graph& g) {
    Graph cpy;
    auto iso_map = vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<Graph>::null_vertex());
    copy_graph(g, cpy, iso_map.ref());

    bool parallel_edge_exists = (std::ranges::count(adjacent_vertices(u, g), v) != 0);

    auto [e, added] = add_edge(u, v, g);

    if (is_undirected(g) && u == v) {
      EXPECT_FALSE(added);
    } else if (parallel_edge_exists) {
      if (allows_parallel_edges(g)) {
        EXPECT_TRUE(added);
      } else {
        EXPECT_FALSE(added);
      }
    } else {
      EXPECT_TRUE(added);
    }

    if (added) {
      EXPECT_EQ(num_edges(g), num_edges(cpy) + 1);
      EXPECT_NE(std::ranges::count(out_edges(u, g), e), 0);
      EXPECT_TRUE(verify_isomorphism(make_filtered_graph(g, ignore_edge(e)), cpy, iso_map.ref()));
    } else {
      if (!(is_undirected(g) && u == v)) {
        // e should be a parallel edge
        EXPECT_EQ(source(e, g), u);
        EXPECT_EQ(target(e, g), v);
      }
      // The graph should not be changed.
      EXPECT_TRUE(verify_isomorphism(g, cpy, iso_map.ref()));
    }
  }  // test_add_edge()

  void test_remove_edge(Vertex u, Vertex v, Graph& g) {
    Graph cpy;
    auto iso_map = vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<Graph>::null_vertex());
    copy_graph(g, cpy, iso_map.ref());

    std::size_t occurances = std::ranges::count(adjacent_vertices(u, g), v);

    remove_edge(u, v, g);

    EXPECT_EQ(num_edges(g) + occurances, num_edges(cpy));
    EXPECT_TRUE(verify_isomorphism(g, make_filtered_graph(cpy, ignore_edges(u, v, cpy)), iso_map.ref()));
  }

  void test_remove_edge(Edge e, Graph& g) {
    Graph cpy;
    auto iso_map = vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<Graph>::null_vertex());
    copy_graph(g, cpy, iso_map.ref());

    Vertex u = source(e, g), v = target(e, g);
    std::size_t occurances = count(adjacent_vertices(u, g), v);

    remove_edge(e, g);

    EXPECT_EQ(num_edges(g) + 1, num_edges(cpy));
    EXPECT_EQ(std::ranges::count(adjacent_vertices(u, g), v) + 1, occurances);
    EXPECT_TRUE(verify_isomorphism(g, make_filtered_graph(cpy, ignore_edge(e)), iso_map.ref()));
  }

  void test_clear_vertex(Vertex v, Graph& g) {
    Graph cpy;
    auto iso_map = vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<Graph>::null_vertex());
    copy_graph(g, cpy, iso_map.ref());

    clear_vertex(v, g);

    EXPECT_EQ(out_degree(v, g), 0);
    EXPECT_EQ(num_vertices(g), num_vertices(cpy));
    EXPECT_TRUE(verify_isomorphism(g, make_filtered_graph(cpy, keep_all(), ignore_vertex(v)), iso_map.ref()));
  }

  //=========================================================================
  // Property Map

  template <typename PropVal, typename PropertyTag>
  void test_readable_vertex_property_graph(const std::vector<PropVal>& vertex_prop, PropertyTag tag, const Graph& g) {
    using PMap = property_map_const_t<Graph, PropertyTag>;
    PMap pmap = get(tag, g);
    auto i = vertex_prop.begin();
    for (auto v : vertices(g)) {
      auto pval1 = get(pmap, v);
      auto pval2 = get(tag, g, v);
      EXPECT_EQ(pval1, pval2);
      EXPECT_EQ(pval1, *i++);
    }
  }

  template <typename PropVal, typename PropertyTag>
  void test_vertex_property_graph(const std::vector<PropVal>& vertex_prop, PropertyTag tag, Graph& g) {
    using PMap = property_map_t<Graph, PropertyTag>;
    PMap pmap = get(tag, g);
    auto i = vertex_prop.begin();
    for (auto v : vertices(g)) {
      put(pmap, v, *i++);
    }

    test_readable_vertex_property_graph(vertex_prop, tag, g);

    for (auto v : vertices(g)) {
      put(pmap, v, vertex_prop[0]);
    }

    auto j = vertex_prop.begin();
    for (auto v : vertices(g)) {
      put(tag, g, v, *j++);
    }

    test_readable_vertex_property_graph(vertex_prop, tag, g);
  }

  template <typename PropVal, typename PropertyTag>
  void test_readable_edge_property_graph(const std::vector<PropVal>& edge_prop, PropertyTag tag, const Graph& g) {
    using PMap = property_map_const_t<Graph, PropertyTag>;
    PMap pmap = get(tag, g);
    auto i = edge_prop.begin();
    for (auto e : edges(g)) {
      auto pval1 = get(pmap, e);
      auto pval2 = get(tag, g, e);
      EXPECT_EQ(pval1, pval2);
      EXPECT_EQ(pval1, *i++);
    }
  }

  template <typename PropVal, typename PropertyTag>
  void test_edge_property_graph(const std::vector<PropVal>& edge_prop, PropertyTag tag, Graph& g) {
    using PMap = property_map_t<Graph, PropertyTag>;
    PMap pmap = get(tag, g);
    auto i = edge_prop.begin();
    for (auto e : edges(g)) {
      put(pmap, e, *i++);
    }

    test_readable_edge_property_graph(edge_prop, tag, g);

    for (auto e : edges(g)) {
      put(pmap, e, edge_prop[0]);
    }

    auto j = edge_prop.begin();
    for (auto e : edges(g)) {
      put(tag, g, e, *j++);
    }

    test_readable_edge_property_graph(edge_prop, tag, g);
  }
};

// This is a helper function to recusively compare two subgraphs,
// including the index for every local edges and their children.
template <typename SubGraph>
void sub_cmp(SubGraph const& g1, SubGraph const& g2) {
  EXPECT_EQ(g1.is_root(), g2.is_root());
  EXPECT_EQ(num_vertices(g1), num_vertices(g2));
  EXPECT_EQ(num_edges(g1), num_edges(g2));
  auto e1_rg = edges(g1);
  auto e2_rg = edges(g2);
  auto e1_i = e1_rg.begin();
  auto e2_i = e2_rg.begin();
  for (; e1_i != e1_rg.end(); ++e1_i, ++e2_i) {
    EXPECT_EQ(get(edge_index, g1, *e1_i), get(edge_index, g2, *e2_i));
  }
  auto g1_rg = g1.children();
  auto g2_rg = g2.children();
  auto g1_i = g1_rg.begin();
  auto g2_i = g2_rg.begin();
  for (; g1_i != g1_rg.end() && g2_i != g2_rg.end(); ++g1_i, ++g2_i) {
    sub_cmp(*g1_i, *g2_i);
  }
  EXPECT_EQ(g1_i, g1_rg.end());
  EXPECT_EQ(g2_i, g2_rg.end());
}

TEST(SubGraphTest, Basic) {
  using Graph = adjacency_list<vec_s, vec_s, bidirectional_s, property<vertex_color_t, int>,
                               property<edge_index_t, std::size_t, property<edge_weight_t, int> > >;
  using SubGraph = subgraph<Graph>;
  using Vertex = graph_vertex_descriptor_t<SubGraph>;
  using Edge = graph_edge_descriptor_t<SubGraph>;

  std::mt19937 gen{42};
  for (int t = 0; t < 100; t += 5) {
    SubGraph g;
    int N = t + 2;
    std::vector<Vertex> vertex_set;
    std::vector<Edge> edge_set;
    generate_random_graph(g, N, N * 2, gen, std::back_inserter(vertex_set), std::back_inserter(edge_set));

    graph_test<SubGraph> gt;

    gt.test_incidence_graph(vertex_set, edge_set, g);
    gt.test_bidirectional_graph(vertex_set, edge_set, g);
    gt.test_adjacency_graph(vertex_set, edge_set, g);
    gt.test_vertex_list_graph(vertex_set, g);
    gt.test_edge_list_graph(vertex_set, edge_set, g);
    gt.test_adjacency_matrix(vertex_set, edge_set, g);

    std::vector<Vertex> sub_vertex_set;
    std::vector<Vertex> sub_global_map;
    std::vector<Vertex> global_sub_map(num_vertices(g), graph_traits<Graph>::null_vertex());
    std::vector<Edge> sub_edge_set;

    SubGraph& g_s = g.create_subgraph();

    const std::size_t Nsub = N / 2;

    // Collect a set of random vertices to put in the subgraph
    std::set<Vertex> verts;
    while (verts.size() < Nsub) {
      verts.insert(random_vertex(g, gen));
    }

    for (Vertex v_global : verts) {
      Vertex v = add_vertex(v_global, g_s);
      sub_vertex_set.push_back(v);
      sub_global_map.push_back(v_global);
      global_sub_map[v_global] = v;
    }

    // compute induced edges
    for (auto e : edges(g)) {
      if ((std::ranges::count(sub_global_map, source(e, g)) != 0) &&
          (std::ranges::count(sub_global_map, target(e, g)) != 0)) {
        auto [e_sub, e_sub_found] = edge(global_sub_map[source(e, g)], global_sub_map[target(e, g)], g_s);
        EXPECT_TRUE(e_sub_found);
        if (e_sub_found) {
          sub_edge_set.push_back(e_sub);
        }
      }
    }

    gt.test_incidence_graph(sub_vertex_set, sub_edge_set, g_s);
    gt.test_bidirectional_graph(sub_vertex_set, sub_edge_set, g_s);
    gt.test_adjacency_graph(sub_vertex_set, sub_edge_set, g_s);
    gt.test_vertex_list_graph(sub_vertex_set, g_s);
    gt.test_edge_list_graph(sub_vertex_set, sub_edge_set, g_s);
    gt.test_adjacency_matrix(sub_vertex_set, sub_edge_set, g_s);

    if (num_vertices(g_s) == 0) {
      continue;
    }

    std::vector<int> weights;
    for (std::size_t i = 0; i < num_vertices(g_s); ++i) {
      weights.push_back(i * 2);
    }
    gt.test_vertex_property_graph(weights, vertex_color, g_s);

    weights.clear();
    for (std::size_t i = 0; i < num_edges(g_s); ++i) {
      weights.push_back(i * 2);
    }
    gt.test_edge_property_graph(weights, edge_weight, g_s);
  }
}

struct test_node_bundle {
  int color;
};

struct test_arc_bundle {
  int weight;
};

TEST(SubGraphTest, Bundled) {
  using ArcProp = property<edge_index_t, std::size_t, test_arc_bundle>;
  using Graph = adjacency_list<vec_s, vec_s, bidirectional_s, test_node_bundle, ArcProp>;
  using SubGraph = subgraph<Graph>;
  using Vertex = graph_vertex_descriptor_t<SubGraph>;
  using Edge = graph_edge_descriptor_t<SubGraph>;

  std::mt19937 gen{42};
  for (int t = 0; t < 100; t += 5) {
    SubGraph g;
    int N = t + 2;
    std::vector<Vertex> vertex_set;
    std::vector<Edge> edge_set;

    generate_random_graph(g, N, N * 2, gen, std::back_inserter(vertex_set), std::back_inserter(edge_set));
    graph_test<SubGraph> gt;

    gt.test_incidence_graph(vertex_set, edge_set, g);
    gt.test_bidirectional_graph(vertex_set, edge_set, g);
    gt.test_adjacency_graph(vertex_set, edge_set, g);
    gt.test_vertex_list_graph(vertex_set, g);
    gt.test_edge_list_graph(vertex_set, edge_set, g);
    gt.test_adjacency_matrix(vertex_set, edge_set, g);

    std::vector<Vertex> sub_vertex_set;
    std::vector<Vertex> sub_global_map;
    std::vector<Vertex> global_sub_map(num_vertices(g), graph_traits<SubGraph>::null_vertex());
    std::vector<Edge> sub_edge_set;

    SubGraph& g_s = g.create_subgraph();

    const std::size_t Nsub = N / 2;

    // Collect a set of random vertices to put in the subgraph
    std::set<Vertex> verts;
    while (verts.size() < Nsub) {
      verts.insert(random_vertex(g, gen));
    }

    for (Vertex v_global : verts) {
      Vertex v = add_vertex(v_global, g_s);
      sub_vertex_set.push_back(v);
      sub_global_map.push_back(v_global);
      global_sub_map[v_global] = v;
    }

    // compute induced edges
    for (auto e : edges(g)) {
      if ((std::ranges::count(sub_global_map, source(e, g)) != 0) &&
          (std::ranges::count(sub_global_map, target(e, g)) != 0)) {
        auto [e_sub, e_sub_found] = edge(global_sub_map[source(e, g)], global_sub_map[target(e, g)], g_s);
        EXPECT_TRUE(e_sub_found);
        if (e_sub_found) {
          sub_edge_set.push_back(e_sub);
        }
      }
    }
    gt.test_incidence_graph(sub_vertex_set, sub_edge_set, g_s);
    gt.test_bidirectional_graph(sub_vertex_set, sub_edge_set, g_s);
    gt.test_adjacency_graph(sub_vertex_set, sub_edge_set, g_s);
    gt.test_vertex_list_graph(sub_vertex_set, g_s);
    gt.test_edge_list_graph(sub_vertex_set, sub_edge_set, g_s);
    gt.test_adjacency_matrix(sub_vertex_set, sub_edge_set, g_s);

    if (num_vertices(g_s) == 0) {
      continue;
    }

    std::vector<int> weights;
    for (std::size_t i = 0; i < num_vertices(g_s); ++i) {
      weights.push_back(i * 2);
    }
    gt.test_vertex_property_graph(weights, &test_node_bundle::color, g_s);

    weights.clear();
    for (std::size_t i = 0; i < num_edges(g_s); ++i) {
      weights.push_back(i * 2);
    }
    gt.test_edge_property_graph(weights, &test_arc_bundle::weight, g_s);
  }
}

TEST(SubGraphTest, LocalGlobalTaggedProperties) {
  using Graph = adjacency_list<vec_s, vec_s, bidirectional_s, property<vertex_color_t, int>,
                               property<edge_index_t, std::size_t, property<edge_weight_t, int> > >;
  using SubGraph = subgraph<Graph>;

  // Create a graph with some vertices.
  SubGraph g(5);
  auto v_rg = vertices(g);
  auto v0 = *v_rg.begin();
  auto v1 = *std::next(v_rg.begin());

  // Create a child subgraph and add some vertices.
  SubGraph& sg = g.create_subgraph();
  auto v = add_vertex(v0, sg);
  add_vertex(v1, sg);
  auto [e, e_added] = add_edge(v0, v1, sg);
  EXPECT_TRUE(e_added);

  using DefaultVertMap = property_map_t<SubGraph, vertex_color_t>;
  DefaultVertMap map = get(vertex_color, g);
  EXPECT_EQ(get(map, v), 0);
  put(map, v, 5);
  EXPECT_EQ(get(map, v), 5);

  using GlobalVertProp = global_property<vertex_color_t>;
  using GlobalVertMap = property_map_t<SubGraph, GlobalVertProp>;
  GlobalVertMap groot = get(global(vertex_color), g);
  GlobalVertMap gsub = get(global(vertex_color), sg);
  EXPECT_EQ(get(groot, v), 5);
  EXPECT_EQ(get(gsub, v), 5);
  put(gsub, v, 10);
  EXPECT_EQ(get(groot, v), 10);
  EXPECT_EQ(get(gsub, v), 10);
  EXPECT_EQ(get(map, v), 10);

  using LocalVertProp = local_property<vertex_color_t>;
  using LocalVertMap = property_map_t<SubGraph, LocalVertProp>;
  LocalVertMap lroot = get(local(vertex_color), g);  // Actually global!
  LocalVertMap lsub = get(local(vertex_color), sg);
  EXPECT_EQ(get(lroot, v), 10);  // Recall it's 10 from above!
  EXPECT_EQ(get(lsub, v), 0);
  put(lsub, v, 5);
  EXPECT_EQ(get(lsub, v), 5);
  EXPECT_EQ(get(lroot, v), 10);  // Don't change the root prop
  EXPECT_EQ(get(map, v), 10);    // Don't change the root prop

  using DefaultEdgeMap = property_map_t<SubGraph, edge_weight_t>;
  DefaultEdgeMap emap = get(edge_weight, g);
  EXPECT_EQ(get(emap, e), 0);
  put(emap, e, 5);
  EXPECT_EQ(get(emap, e), 5);

  using GlobalEdgeProp = global_property<edge_weight_t>;
  using GlobalEdgeMap = property_map_t<SubGraph, GlobalEdgeProp>;
  GlobalEdgeMap geroot = get(global(edge_weight), g);
  GlobalEdgeMap gesub = get(global(edge_weight), sg);
  EXPECT_EQ(get(geroot, e), 5);
  EXPECT_EQ(get(gesub, e), 5);
  put(gesub, e, 10);
  EXPECT_EQ(get(geroot, e), 10);
  EXPECT_EQ(get(gesub, e), 10);
  EXPECT_EQ(get(emap, e), 10);

  using LocalEdgeProp = local_property<edge_weight_t>;
  using LocalEdgeMap = property_map_t<SubGraph, LocalEdgeProp>;
  LocalEdgeMap leroot = get(local(edge_weight), g);  // Actually global!
  LocalEdgeMap lesub = get(local(edge_weight), sg);
  EXPECT_EQ(get(leroot, e), 10);  // Recall it's 10 from above!
  EXPECT_EQ(get(lesub, e), 0);
  put(lesub, e, 5);
  EXPECT_EQ(get(lesub, e), 5);
  EXPECT_EQ(get(leroot, e), 10);  // Don't change the root prop
  EXPECT_EQ(get(emap, e), 10);    // Don't change the root prop
}

TEST(SubGraphTest, LocalGlobalBundledProperties) {
  using ArcProp = property<edge_index_t, std::size_t, test_arc_bundle>;
  using Graph = adjacency_list<vec_s, vec_s, bidirectional_s, test_node_bundle, ArcProp>;
  using SubGraph = subgraph<Graph>;

  // Create a graph with some vertices.
  SubGraph g(5);
  auto v_rg = vertices(g);
  auto v0 = *v_rg.begin();
  auto v1 = *std::next(v_rg.begin());

  // Create a child subgraph and add some vertices.
  SubGraph& sg = g.create_subgraph();
  auto v = add_vertex(v0, sg);
  add_vertex(v1, sg);
  auto [e, e_added] = add_edge(v0, v1, sg);
  EXPECT_TRUE(e_added);

  using DefaultVertMap = property_map_t<SubGraph, int test_node_bundle::*>;
  DefaultVertMap map = get(&test_node_bundle::color, g);
  EXPECT_EQ(get(map, v), 0);
  put(map, v, 5);
  EXPECT_EQ(get(map, v), 5);

  using GlobalVertProp = global_property<int test_node_bundle::*>;
  using GlobalVertMap = property_map_t<SubGraph, GlobalVertProp>;
  GlobalVertMap groot = get(global(&test_node_bundle::color), g);
  GlobalVertMap gsub = get(global(&test_node_bundle::color), sg);
  EXPECT_EQ(get(groot, v), 5);
  EXPECT_EQ(get(gsub, v), 5);
  put(gsub, v, 10);
  EXPECT_EQ(get(groot, v), 10);
  EXPECT_EQ(get(gsub, v), 10);
  EXPECT_EQ(get(map, v), 10);

  using LocalVertProp = local_property<int test_node_bundle::*>;
  using LocalVertMap = property_map_t<SubGraph, LocalVertProp>;
  LocalVertMap lroot = get(local(&test_node_bundle::color), g);  // Actually global!
  LocalVertMap lsub = get(local(&test_node_bundle::color), sg);
  EXPECT_EQ(get(lroot, v), 10);  // Recall it's 10 from above!
  EXPECT_EQ(get(lsub, v), 0);
  put(lsub, v, 5);
  EXPECT_EQ(get(lsub, v), 5);
  EXPECT_EQ(get(lroot, v), 10);  // Don't change the root prop
  EXPECT_EQ(get(map, v), 10);    // Don't change the root prop

  using DefaultEdgeMap = property_map_t<SubGraph, int test_arc_bundle::*>;
  DefaultEdgeMap emap = get(&test_arc_bundle::weight, g);
  EXPECT_EQ(get(emap, e), 0);
  put(emap, e, 5);
  EXPECT_EQ(get(emap, e), 5);

  using GlobalEdgeProp = global_property<int test_arc_bundle::*>;
  using GlobalEdgeMap = property_map_t<SubGraph, GlobalEdgeProp>;
  GlobalEdgeMap geroot = get(global(&test_arc_bundle::weight), g);
  GlobalEdgeMap gesub = get(global(&test_arc_bundle::weight), sg);
  EXPECT_EQ(get(geroot, e), 5);
  EXPECT_EQ(get(gesub, e), 5);
  put(gesub, e, 10);
  EXPECT_EQ(get(geroot, e), 10);
  EXPECT_EQ(get(gesub, e), 10);
  EXPECT_EQ(get(emap, e), 10);
  sg[global(e)].weight = 20;
  EXPECT_EQ(get(geroot, e), 20);
  EXPECT_EQ(get(gesub, e), 20);
  EXPECT_EQ(get(emap, e), 20);

  using LocalEdgeProp = local_property<int test_arc_bundle::*>;
  using LocalEdgeMap = property_map_t<SubGraph, LocalEdgeProp>;
  LocalEdgeMap leroot = get(local(&test_arc_bundle::weight), g);  // Actually global!
  LocalEdgeMap lesub = get(local(&test_arc_bundle::weight), sg);
  EXPECT_EQ(get(leroot, e), 20);  // Recall it's 10 from above!
  EXPECT_EQ(get(lesub, e), 0);
  put(lesub, e, 5);
  EXPECT_EQ(get(lesub, e), 5);
  EXPECT_EQ(get(leroot, e), 20);  // Don't change the root prop
  EXPECT_EQ(get(emap, e), 20);    // Don't change the root prop
  sg[local(e)].weight = 10;
  EXPECT_EQ(get(lesub, e), 10);
  EXPECT_EQ(get(leroot, e), 20);  // Don't change the root prop
  EXPECT_EQ(get(emap, e), 20);    // Don't change the root prop
}

// This block generates a complete graph with 8 vertices,
// and puts the first and last four of the vertices into two children.
// Do these again to the children, so there are 4 grandchildren with 2
// vertices for each. Use the copy constructor to generate a copy and
// compare with the original one.
TEST(SubGraphTest, SimpleCopy) {
  using Graph = adjacency_list<vec_s, vec_s, bidirectional_s, property<vertex_color_t, int>,
                               property<edge_index_t, std::size_t, property<edge_weight_t, int> > >;
  using SubGraph = subgraph<Graph>;

  SubGraph g1;

  for (size_t i = 0; i < 8; i++) {
    add_vertex(g1);
  }
  for (auto u : vertices(g1)) {
    for (auto v : vertices(g1)) {
      if (u != v) {
        add_edge(u, v, g1);
      }
    }
  }
  auto v1_rg = vertices(g1);
  auto vi = std::next(v1_rg.begin(), 4);
  g1.create_subgraph(std::ranges::subrange(v1_rg.begin(), vi));
  g1.create_subgraph(std::ranges::subrange(std::next(vi), v1_rg.end()));
  auto c1_rg = g1.children();
  auto gi1 = c1_rg.begin();
  auto gi2 = std::next(gi1);
  auto cv1_rg = vertices(*gi1);
  auto cv1i = std::next(cv1_rg.begin(), 2);
  auto cv2_rg = vertices(*gi2);
  auto cv2i = std::next(cv2_rg.begin(), 2);
  (*gi1).create_subgraph(std::ranges::subrange(cv1_rg.begin(), cv1i));
  (*gi1).create_subgraph(std::ranges::subrange(std::next(cv1i), cv1_rg.end()));
  (*gi2).create_subgraph(std::ranges::subrange(cv2_rg.begin(), cv2i));
  (*gi2).create_subgraph(std::ranges::subrange(std::next(cv2i), cv2_rg.end()));
  SubGraph g2(g1);
  sub_cmp(g1, g2);
}

// A regression test: the copy constructor of subgraph did not
// copy one of the members, so local_edge->global_edge mapping
// was broken.
TEST(SubGraphTest, RegressionTestEdgeMapping) {
  using Graph = adjacency_list<vec_s, vec_s, bidirectional_s, property<vertex_color_t, int>,
                               property<edge_index_t, std::size_t, property<edge_weight_t, int> > >;
  using SubGraph = subgraph<Graph>;

  SubGraph g;
  auto v1 = add_vertex(g);
  auto v2 = add_vertex(g);
  add_edge(v1, v2, g);

  SubGraph sub = g.create_subgraph(vertices(g));

  for (auto e : edges(sub)) {
    // This used to segfault.
    get(edge_weight, sub, e);
  }
}

// A regression test: the copy constructor of subgraph did not
// copy one of the members, so local_edge->global_edge mapping
// was broken.
TEST(SubGraphTest, RegressionTestEdgeBundled) {
  using ArcProp = property<edge_index_t, std::size_t, test_arc_bundle>;
  using Graph = adjacency_list<vec_s, vec_s, bidirectional_s, test_node_bundle, ArcProp>;
  using SubGraph = subgraph<Graph>;

  SubGraph g;
  auto v1 = add_vertex(g);
  auto v2 = add_vertex(g);
  add_edge(v1, v2, g);

  SubGraph sub = g.create_subgraph(vertices(g));

  for (auto e : edges(sub)) {
    // This used to segfault.
    get(&test_arc_bundle::weight, sub, e);
  }
}

TEST(SubGraphTest, SimpleGraphTest) {
  using Graph = subgraph<adjacency_list<vec_s, vec_s, directed_s, no_property, property<edge_index_t, int> > >;

  const int N = 6;
  Graph g0(N);

  // for conveniently refering to vertices in g0
  enum { A, B, C, D, E, F };

  Graph& g1 = g0.create_subgraph();
  Graph& g2 = g1.create_subgraph();

  EXPECT_EQ(&g1.parent(), &g0);
  EXPECT_EQ(&g2.parent(), &g1);

  // for conveniently refering to vertices in g1
  enum { A1, B1, C1 };
  // for conveniently refering to vertices in g2
  enum { A2, B2, C2 };

  add_vertex(C, g1);  // global vertex C becomes local A1 for g1
  add_vertex(E, g1);  // global vertex E becomes local B1 for g1
  add_vertex(F, g1);  // global vertex F becomes local C1 for g1

  add_vertex(C, g2);  // global vertex C becomes local A2 for g2
  add_vertex(E, g2);  // global vertex E becomes local B2 for g2

  EXPECT_EQ(num_vertices(g0), 6);
  EXPECT_EQ(num_vertices(g1), 3);
  EXPECT_EQ(num_vertices(g2), 2);

  // add edges to root graph
  add_edge(A, B, g0);
  add_edge(B, C, g0);
  add_edge(B, D, g0);
  add_edge(E, F, g0);

  EXPECT_EQ(num_edges(g0), 4);
  EXPECT_EQ(num_edges(g1), 1);
  EXPECT_EQ(num_edges(g2), 0);

  // add edges to g1
  add_edge(A1, B1, g1);
  EXPECT_EQ(num_edges(g0), 5);
  EXPECT_EQ(num_edges(g1), 2);
  EXPECT_EQ(num_edges(g2), 1);
  // num_vertices stays the same
  EXPECT_EQ(num_vertices(g0), 6);
  EXPECT_EQ(num_vertices(g1), 3);
  EXPECT_EQ(num_vertices(g2), 2);
}

TEST(SubGraphTest, AddVerticesTest) {
  using Graph = subgraph<adjacency_list<vec_s, vec_s, directed_s, no_property, property<edge_index_t, int> > >;

  const int N = 3;
  Graph g0(N);
  Graph& g1 = g0.create_subgraph();
  Graph& g2 = g1.create_subgraph();

  EXPECT_EQ(&g1.parent(), &g0);
  EXPECT_EQ(&g2.parent(), &g1);

  // add vertices to g2
  auto n1 = add_vertex(0, g2);
  auto n2 = add_vertex(1, g2);
  // check if the global vertex 2 is equal to the returned local vertex
  EXPECT_TRUE(g2.find_vertex(0).second);
  EXPECT_EQ(g2.find_vertex(0).first, n1);
  EXPECT_TRUE(g2.find_vertex(1).second);
  EXPECT_EQ(g2.find_vertex(1).first, n2);
  // and check if this vertex is also present in g1
  EXPECT_TRUE(g1.find_vertex(0).second);
  EXPECT_EQ(g1.local_to_global(g1.find_vertex(0).first), 0);
  EXPECT_TRUE(g1.find_vertex(0).second);
  EXPECT_EQ(g1.local_to_global(g1.find_vertex(1).first), 1);

  // num_vertices stays the same
  EXPECT_EQ(num_vertices(g0), 3);
  EXPECT_EQ(num_vertices(g1), 2);
  EXPECT_EQ(num_vertices(g2), 2);

  // add vertices to g1
  auto n3 = add_vertex(2, g1);
  // check if the global vertex 2 is equal to the returned local vertex
  EXPECT_TRUE(g1.find_vertex(2).second);
  EXPECT_EQ(g1.find_vertex(2).first, n3);
  // num_vertices stays the same
  EXPECT_EQ(num_vertices(g0), 3);
  EXPECT_EQ(num_vertices(g1), 3);
  EXPECT_EQ(num_vertices(g2), 2);

  // add vertices to g1
  auto n4 = add_vertex(g1);

  // check if the new local vertex is also in the global graph
  EXPECT_TRUE(g0.find_vertex(g1.local_to_global(n4)).second);
  // check if the new local vertex is not in the subgraphs
  EXPECT_FALSE(g2.find_vertex(n4).second);

  // num_vertices stays the same
  EXPECT_EQ(num_vertices(g0), 4);
  EXPECT_EQ(num_vertices(g1), 4);
  EXPECT_EQ(num_vertices(g2), 2);

  // add vertices to g0
  auto n5 = add_vertex(g0);

  // check if the new local vertex is not in the subgraphs
  EXPECT_FALSE(g1.find_vertex(n5).second);
  EXPECT_FALSE(g2.find_vertex(n5).second);

  // num_vertices stays the same
  EXPECT_EQ(num_vertices(g0), 5);
  EXPECT_EQ(num_vertices(g1), 4);
  EXPECT_EQ(num_vertices(g2), 2);

  if (::testing::Test::HasFailure()) {
    std::cerr << "All g0 vertices: " << std::endl;
    for (auto v : g0.m_local_vertex) {
      std::cerr << g0.local_to_global(v.first) << std::endl;
    }
    std::cerr << "All g1 vertices: " << std::endl;
    for (auto v : g1.m_local_vertex) {
      std::cerr << g1.local_to_global(v.first) << std::endl;
    }
    std::cerr << "All g2 vertices: " << std::endl;
    for (auto v : g2.m_local_vertex) {
      std::cerr << g2.local_to_global(v.first) << std::endl;
    }
  }
}

TEST(SubGraphTest, AddEdgesTest) {
  using Graph = subgraph<adjacency_list<vec_s, vec_s, directed_s, no_property, property<edge_index_t, int> > >;

  const int N = 3;
  Graph g0(N);
  Graph& g1 = g0.create_subgraph();
  Graph& g2 = g1.create_subgraph();

  EXPECT_EQ(&g1.parent(), &g0);
  EXPECT_EQ(&g2.parent(), &g1);

  // add vertices
  add_vertex(0, g2);
  add_vertex(1, g2);
  EXPECT_EQ(num_vertices(g1), 2);
  EXPECT_EQ(num_vertices(g2), 2);

  // add edge to g0 which needs propagation
  add_edge(0, 1, g0);

  EXPECT_EQ(num_edges(g0), 1);
  EXPECT_EQ(num_edges(g1), 1);
  EXPECT_EQ(num_edges(g2), 1);
  // num_vertices stays the same
  EXPECT_EQ(num_vertices(g0), 3);
  EXPECT_EQ(num_vertices(g1), 2);
  EXPECT_EQ(num_vertices(g2), 2);

  // add edge to g0 without propagation
  add_edge(1, 2, g0);

  EXPECT_EQ(num_edges(g0), 2);
  EXPECT_EQ(num_edges(g1), 1);
  EXPECT_EQ(num_edges(g2), 1);
  // num_vertices stays the same
  EXPECT_EQ(num_vertices(g0), 3);
  EXPECT_EQ(num_vertices(g1), 2);
  EXPECT_EQ(num_vertices(g2), 2);

  // add vertex 2 to g2/g1 with edge propagation
  auto n = add_vertex(2, g2);
  EXPECT_EQ(g2.local_to_global(n), 2);

  EXPECT_EQ(num_edges(g0), 2);
  EXPECT_EQ(num_edges(g1), 2);
  EXPECT_EQ(num_edges(g2), 2);
  // num_vertices stays the same
  EXPECT_EQ(num_vertices(g0), 3);
  EXPECT_EQ(num_vertices(g1), 3);
  EXPECT_EQ(num_vertices(g2), 3);

  // add edge to g2 with propagation upwards
  add_edge(0, 2, g2);

  EXPECT_EQ(num_edges(g0), 3);
  EXPECT_EQ(num_edges(g1), 3);
  EXPECT_EQ(num_edges(g2), 3);
  // num_vertices stays the same
  EXPECT_EQ(num_vertices(g0), 3);
  EXPECT_EQ(num_vertices(g1), 3);
  EXPECT_EQ(num_vertices(g2), 3);

  if (::testing::Test::HasFailure()) {
    std::cerr << "All g0 vertices: " << std::endl;
    for (auto v : g0.m_local_vertex) {
      std::cerr << g0.local_to_global(v.first) << std::endl;
    }
    std::cerr << "All g1 vertices: " << std::endl;
    for (auto v : g1.m_local_vertex) {
      std::cerr << g1.local_to_global(v.first) << std::endl;
    }
    std::cerr << "All g2 vertices: " << std::endl;
    for (auto v : g2.m_local_vertex) {
      std::cerr << g2.local_to_global(v.first) << std::endl;
    }
    std::cerr << "All g0 edges: " << std::endl;
    for (auto e : edges(g0)) {
      std::cerr << source(e, g0) << "->" << target(e, g0) << std::endl;
    }
    std::cerr << "All g1 edges: " << std::endl;
    for (auto e : edges(g1)) {
      std::cerr << source(e, g1) << "->" << target(e, g1) << std::endl;
    }
    std::cerr << "All g2 edges: " << std::endl;
    for (auto e : edges(g2)) {
      std::cerr << source(e, g2) << "->" << target(e, g2) << std::endl;
    }
  }
}

}  // namespace
}  // namespace bagl
