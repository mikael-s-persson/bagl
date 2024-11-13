// Copyright 2004 The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/circle_layout.h"
#include "bagl/fruchterman_reingold.h"
#include "bagl/graph_traits.h"
#include "bagl/kamada_kawai_spring_layout.h"
#include "bagl/point_traits.h"
#include "bagl/random_layout.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

using TestPoint = square_topology<>::point_type;
using TestGraph = adjacency_list<list_s, list_s, undirected_s,
                                 // Vertex properties
                                 property<vertex_index_t, int, property<vertex_position_t, TestPoint>>,
                                 // Edge properties
                                 property<edge_weight_t, double>>;

template <typename Graph, typename PositionMap, typename Topology>
void print_graph_layout(const Graph& g, PositionMap position, const Topology& topology) {
  // Find min/max ranges
  auto min_point = position[*vertices(g).begin()];
  auto max_point = min_point;
  for (auto v : vertices(g)) {
    min_point = topology.pointwise_min(min_point, position[v]);
    max_point = topology.pointwise_max(max_point, position[v]);
  }

  for (auto y = min_point[1]; y <= max_point[1]; ++y) {
    for (auto x = min_point[0]; x <= max_point[0]; ++x) {
      // Find vertex at this position
      std::size_t index = 0;
      for (auto v : vertices(g)) {
        if (position[v][0] == x && position[v][1] == y) {
          break;
        }
        ++index;
      }

      if (index == num_vertices(g)) {
        std::cout << ' ';
      } else {
        std::cout << static_cast<char>(index + 'A');
      }
    }
    std::cout << std::endl;
  }
}

template <typename Graph, typename PositionMap>
void dump_graph_layout(std::string name, const Graph& g, PositionMap position) {
  std::ofstream out((name + ".dot").c_str());
  out << "graph " << name << " {" << std::endl;

  for (auto v : vertices(g)) {
    out << "  n" << get(vertex_index, g, v) << "[ pos=\"" << (position[v][0] + 25) << ", " << (position[v][1] + 25)
        << "\" ];\n";
  }

  for (auto e : edges(g)) {
    out << "  n" << get(vertex_index, g, source(e, g)) << " -- n" << get(vertex_index, g, target(e, g)) << ";\n";
  }
  out << "}\n";
}

TEST(LayoutTest, Circle) {
  using Graph = TestGraph;
  constexpr std::size_t n = 5;

  Graph g(n);

  // Initialize vertex indices
  std::size_t i = 0;
  for (auto v : vertices(g)) {
    put(vertex_index, g, v, i);
  }

  circle_graph_layout(g, get(vertex_position, g), 10.0);

  std::cout << "Regular polygon layout with " << n << " points.\n";
  square_topology<> topology;
  print_graph_layout(g, get(vertex_position, g), topology);
}

struct kamada_kawai_done {
  kamada_kawai_done() : last_delta() {}

  template <typename Graph>
  bool operator()(double delta_p, graph_vertex_descriptor_t<Graph> /*p*/, const Graph& /*g*/, bool global) {
    if (global) {
      double diff = last_delta - delta_p;
      if (diff < 0) {
        diff = -diff;
      }
      last_delta = delta_p;
      return diff < 0.01;
    } else {
      return delta_p < 0.01;
    }
  }

  double last_delta;
};

TEST(LayoutTest, Triangle) {
  using Graph = TestGraph;

  Graph g;

  auto u = add_vertex(g);
  put(vertex_index, g, u, 0);
  auto v = add_vertex(g);
  put(vertex_index, g, v, 1);
  auto w = add_vertex(g);
  put(vertex_index, g, w, 2);

  auto e1 = add_edge(u, v, g).first;
  put(edge_weight, g, e1, 1.0);
  auto e2 = add_edge(v, w, g).first;
  put(edge_weight, g, e2, 1.0);
  auto e3 = add_edge(w, u, g).first;
  put(edge_weight, g, e3, 1.0);

  circle_graph_layout(g, get(vertex_position, g), 25.0);

  EXPECT_TRUE(kamada_kawai_spring_layout(g, get(vertex_position, g), get(edge_weight, g), square_topology<>(50.0),
                                         side_length(50.0)));

  std::cout << "Triangle layout (Kamada-Kawai).\n";
  print_graph_layout(g, get(vertex_position, g), square_topology<>(50.0));
}

TEST(LayoutTest, Cube) {
  using Graph = TestGraph;
  enum { A, B, C, D, E, F, G, H };
  std::vector<std::pair<std::size_t, std::size_t>> cube_edges = {{A, E}, {A, B}, {A, D}, {B, F}, {B, C}, {C, D},
                                                                 {C, G}, {D, H}, {E, H}, {E, F}, {F, G}, {G, H}};

  Graph g(8, cube_edges);

  std::size_t i = 0;
  for (auto v : vertices(g)) {
    put(vertex_index, g, v, i++);
  }

  for (auto e : edges(g)) {
    put(edge_weight, g, e, 1.0);
    std::cerr << "(" << static_cast<char>(get(vertex_index, g, source(e, g)) + 'A') << ", "
              << static_cast<char>(get(vertex_index, g, target(e, g)) + 'A') << ") ";
  }
  std::cerr << std::endl;

  circle_graph_layout(g, get(vertex_position, g), 25.0);

  EXPECT_TRUE(kamada_kawai_spring_layout(g, get(vertex_position, g), get(edge_weight, g), square_topology<>(50.0),
                                         side_length(50.0), kamada_kawai_done()));

  std::cout << "Cube layout (Kamada-Kawai).\n";
  print_graph_layout(g, get(vertex_position, g), square_topology<>(50.));

  dump_graph_layout("/tmp/cube", g, get(vertex_position, g));

  std::mt19937 gen(42);
  using Topology = square_topology<>;
  Topology topology(gen, 50.0);
  auto displacements = vector_property_map(num_vertices(g), get(vertex_index, g), Topology::point_difference_type());
  rectangle_topology<> rect_top(gen, 0, 0, 50, 50);
  random_graph_layout(g, get(vertex_position, g), rect_top);

  fruchterman_reingold_force_directed_layout(g, get(vertex_position, g), topology, square_distance_attractive_force(),
                                             square_distance_repulsive_force(), all_force_pairs(),
                                             linear_cooling<double>(100), displacements.ref());

  std::cout << "Cube layout (Fruchterman-Reingold).\n";
  print_graph_layout(g, get(vertex_position, g), square_topology<>(50.));

  dump_graph_layout("/tmp/cube-fr", g, get(vertex_position, g));
}

TEST(LayoutTest, Triangular) {
  using Graph = TestGraph;
  enum { A, B, C, D, E, F, G, H, I, J };
  std::vector<std::pair<std::size_t, std::size_t>> triangular_edges = {{A, B}, {A, C}, {B, C}, {B, D}, {B, E}, {C, E},
                                                                       {C, F}, {D, E}, {D, G}, {D, H}, {E, F}, {E, H},
                                                                       {E, I}, {F, I}, {F, J}, {G, H}, {H, I}, {I, J}};

  Graph g(10, triangular_edges);

  std::size_t i = 0;
  for (auto v : vertices(g)) {
    put(vertex_index, g, v, i++);
  }

  for (auto e : edges(g)) {
    put(edge_weight, g, e, 1.0);
    std::cerr << "(" << static_cast<char>(get(vertex_index, g, source(e, g)) + 'A') << ", "
              << static_cast<char>(get(vertex_index, g, target(e, g)) + 'A') << ") ";
  }
  std::cerr << std::endl;

  using Topology = square_topology<>;
  std::mt19937 gen(42);
  Topology topology(gen, 50.0);
  Topology::point_type origin;
  origin[0] = origin[1] = 50.0;
  Topology::point_difference_type extent;
  extent[0] = extent[1] = 50.0;

  circle_graph_layout(g, get(vertex_position, g), 25.0);

  EXPECT_TRUE(kamada_kawai_spring_layout(g, get(vertex_position, g), get(edge_weight, g), topology, side_length(50.0),
                                         kamada_kawai_done()));

  std::cout << "Triangular layout (Kamada-Kawai).\n";
  print_graph_layout(g, get(vertex_position, g), square_topology<>(50.));

  dump_graph_layout("/tmp/triangular-kk", g, get(vertex_position, g));

  rectangle_topology<> rect_top(gen, -25, -25, 25, 25);
  random_graph_layout(g, get(vertex_position, g), rect_top);

  dump_graph_layout("/tmp/random", g, get(vertex_position, g));

  fruchterman_reingold_force_directed_layout(
      g, get(vertex_position, g), topology, square_distance_attractive_force(), square_distance_repulsive_force(),
      make_grid_force_pairs(topology, get(vertex_position, g), g), linear_cooling<double>(100));

  std::cout << "Triangular layout (Fruchterman-Reingold).\n";
  print_graph_layout(g, get(vertex_position, g), square_topology<>(50.));

  dump_graph_layout("/tmp/triangular-fr", g, get(vertex_position, g));
}

TEST(LayoutTest, Disconnected) {
  using Graph = TestGraph;
  enum { A, B, C, D, E, F, G, H };
  std::vector<std::pair<std::size_t, std::size_t>> triangular_edges = {
      {A, B}, {B, C}, {C, A}, {D, E}, {E, F}, {F, G}, {G, H}, {H, D}, {D, F}, {F, H}, {H, E}, {E, G}, {G, D}};

  Graph g(8, triangular_edges);

  std::size_t i = 0;
  for (auto v : vertices(g)) {
    put(vertex_index, g, v, i++);
  }

  for (auto e : edges(g)) {
    put(edge_weight, g, e, 1.0);
    std::cerr << "(" << static_cast<char>(get(vertex_index, g, source(e, g)) + 'A') << ", "
              << static_cast<char>(get(vertex_index, g, target(e, g)) + 'A') << ") ";
  }
  std::cerr << std::endl;

  circle_graph_layout(g, get(vertex_position, g), 25.0);

  EXPECT_TRUE(kamada_kawai_spring_layout(g, get(vertex_position, g), get(edge_weight, g), square_topology<>(50.0),
                                         side_length(50.0), kamada_kawai_done()));

  std::mt19937 gen(42);
  rectangle_topology<> rect_top(gen, -25, -25, 25, 25);
  random_graph_layout(g, get(vertex_position, g), rect_top);

  using Topology = square_topology<>;
  Topology topology(gen, 50.0);
  fruchterman_reingold_force_directed_layout(
      g, get(vertex_position, g), topology, square_distance_attractive_force(), square_distance_repulsive_force(),
      make_grid_force_pairs(topology, get(vertex_position, g), g), linear_cooling<double>(50));

  std::cout << "Disconnected layout (Fruchterman-Reingold).\n";
  print_graph_layout(g, get(vertex_position, g), square_topology<>(50.));

  dump_graph_layout("/tmp/disconnected-fr", g, get(vertex_position, g));
}

}  // namespace
}  // namespace bagl
