// Copyright 2008
// Author: Matyas W Egyhazy
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/metric_tsp_approx.h"

#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/adjacency_matrix.h"
#include "bagl/graph_traits.h"
#include "bagl/simple_point.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

struct ComparePointByX {
  template <typename T>
  bool operator()(const simple_point<T>& l, const simple_point<T>& r) const {
    return (l.x > r.x);
  }
};

// add edges to the graph (for each node connect it to all other nodes)
template <typename VertexListGraph, typename PointContainer, typename WeightMap, typename VertexIndexMap>
void connect_all_euclidean(VertexListGraph& g, const PointContainer& points, WeightMap wmap, VertexIndexMap vmap) {
  auto v_rg = vertices(g);
  for (auto src = v_rg.begin(); src != v_rg.end(); ++src) {
    for (auto dest = std::next(src); dest != v_rg.end(); ++dest) {
      const double dx = points[vmap[*src]].x - points[vmap[*dest]].x;
      const double dy = points[vmap[*src]].y - points[vmap[*dest]].y;
      const double w = std::sqrt(dx * dx + dy * dy);

      auto [e, e_added] = add_edge(*src, *dest, g);
      EXPECT_TRUE(e_added);
      wmap[e] = w;
    }
  }
}

// Create a randomly generated point scatter time execution
TEST(MetricTSPApproxTest, Scalability) {
  const std::size_t numpts = 1000;

  using Graph =
      adjacency_matrix<undirected_s, no_property, property<edge_weight_t, double, property<edge_index_t, std::size_t>>>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  std::mt19937 rng(42);
  std::uniform_real_distribution<> uni_dist(0.01, double{numpts * 2});

  std::set<simple_point<double>, ComparePointByX> points;
  simple_point<double> pnt;

  while (points.size() < numpts) {
    pnt.x = uni_dist(rng);
    pnt.y = uni_dist(rng);
    points.insert(pnt);
  }

  Graph g(numpts);
  auto weight_map = get(edge_weight, g);
  std::vector<simple_point<double>> point_vec(points.begin(), points.end());

  connect_all_euclidean(g, point_vec, weight_map, get(vertex_index, g));

  std::vector<Vertex> c;
  double len = 0.0;

  // Run the TSP approx, creating the visitor on the fly.
  metric_tsp_approx(g, tsp_tour_len_visitor(g, std::back_inserter(c), len, weight_map));

  EXPECT_EQ(c.size(), num_vertices(g) + 1);
  EXPECT_EQ(c.front(), c.back());
  // A reasonable len upper-bound for this size of random problems.
  // This isn't really strict, just basic sanity check that results are not wildly off.
  EXPECT_LT(len, 100'000);
}

template <typename PositionVec, typename IdVec>
void check_adj_list(const PositionVec& pos_vec, const IdVec& expected_tour) {
  using Graph = adjacency_list<list_s, list_s, undirected_s>;
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;

  std::vector<Vertex> c;
  std::unordered_map<Vertex, std::size_t, graph_descriptor_hash_t<Vertex>> v_map;
  auto v_pmap = associative_property_map(v_map);
  std::unordered_map<Edge, double, graph_descriptor_hash_t<Edge>> w_map;
  auto w_pmap = associative_property_map(w_map);

  Graph g(pos_vec.size());

  // create vertex index map
  std::size_t idx = 0;
  for (auto v : vertices(g)) {
    v_pmap[v] = idx++;
  }

  connect_all_euclidean(g, pos_vec, w_pmap, v_pmap);

  metric_tsp_approx_from_vertex(g, *vertices(g).begin(), w_pmap, v_pmap, tsp_tour_visitor(std::back_inserter(c)));
  std::vector<std::size_t> c_id;
  for (auto v : c) {
    c_id.push_back(v_pmap[v]);
  }
  EXPECT_EQ(c_id.size(), num_vertices(g) + 1);
  EXPECT_EQ(c_id.front(), c_id.back());
  EXPECT_THAT(c_id, ::testing::ElementsAreArray(expected_tour));
}

TEST(MetricTSPApproxTest, Basic) {
  using Graph = adjacency_matrix<undirected_s, no_property, property<edge_weight_t, double>>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  // Open the graph file, failing if one isn't given on the command line.
  std::ifstream fin("test/testdata/metric_tsp_approx.graph");
  ASSERT_TRUE(fin);

  std::string line;
  std::vector<simple_point<double>> position_vec;

  while (getline(fin, line)) {
    simple_point<double> vertex;

    std::size_t idx = line.find(",");
    std::string x_str = line.substr(0, idx);
    std::string y_str = line.substr(idx + 1, line.size() - idx);

    vertex.x = std::stod(x_str);
    vertex.y = std::stod(y_str);

    position_vec.push_back(vertex);
  }

  fin.close();

  const std::vector<Vertex> expected_tour = {0, 1, 2, 7, 3, 4, 5, 6, 0};

  Graph g(position_vec.size());
  auto w_map = get(edge_weight, g);
  auto v_map = get(vertex_index, g);

  connect_all_euclidean(g, position_vec, w_map, v_map);

  std::vector<Vertex> c;
  metric_tsp_approx_tour(g, std::back_inserter(c));
  EXPECT_EQ(c.size(), num_vertices(g) + 1);
  EXPECT_EQ(c.front(), c.back());
  EXPECT_THAT(c, ::testing::ElementsAreArray(expected_tour));
  c.clear();

  check_adj_list(position_vec, expected_tour);

  metric_tsp_approx_from_vertex(g, *vertices(g).begin(), w_map, v_map, tsp_tour_visitor(std::back_inserter(c)));
  EXPECT_EQ(c.size(), num_vertices(g) + 1);
  EXPECT_EQ(c.front(), c.back());
  EXPECT_THAT(c, ::testing::ElementsAreArray(expected_tour));
  c.clear();

  double len = 0.0;
  EXPECT_NO_THROW(metric_tsp_approx(g, tsp_tour_len_visitor(g, std::back_inserter(c), len, w_map)));
  EXPECT_EQ(c.size(), num_vertices(g) + 1);
  EXPECT_EQ(c.front(), c.back());
  EXPECT_THAT(c, ::testing::ElementsAreArray(expected_tour));

  ASSERT_GT(c.size(), 3);
  std::pair triangle_edge = {c[2], c[3]};
  c.clear();

  // if the graph is not fully connected then some of the
  // assumed triangle-inequality edges may not exist
  auto [tri_e, tri_e_found] = edge(triangle_edge.first, triangle_edge.second, g);
  EXPECT_TRUE(tri_e_found);
  remove_edge(tri_e, g);

  // Make sure that we can actually trap incomplete graphs.
  len = 0.0;
  EXPECT_THROW(metric_tsp_approx(g, tsp_tour_len_visitor(g, std::back_inserter(c), len, w_map)), bad_graph);
}

}  // namespace
}  // namespace bagl
