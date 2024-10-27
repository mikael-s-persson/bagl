// Copyright 2004 The Trustees of Indiana University.
//  Authors: Jeremiah Willcock
//           Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GURSOY_ATUN_LAYOUT_H_
#define BAGL_BAGL_GURSOY_ATUN_LAYOUT_H_

// Gürsoy-Atun graph layout, based on:
// "Neighbourhood Preserving Load Balancing: A Self-Organizing Approach"
// in 6th International Euro-Par Conference Munich, Germany, August 29 –
// September 1, 2000 Proceedings, pp 234-241
// https://doi.org/10.1007/3-540-44520-X_32

#include <algorithm>
#include <cassert>
#include <cmath>
#include <exception>
#include <memory>
#include <queue>
#include <random>
#include <vector>

#include "bagl/breadth_first_search.h"
#include "bagl/dijkstra_shortest_paths.h"
#include "bagl/filtered_graph.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/single_property_map.h"
#include "bagl/topology.h"
#include "bagl/two_bit_color_map.h"
#include "bagl/vector_property_map.h"
#include "bagl/visitors.h"

namespace bagl {

namespace gursoy_atun_detail {

struct over_distance_limit : public std::exception {};

template <typename PositionMap, typename NodeDistanceMap, typename Topology>
struct update_position_visitor {
  PositionMap position_map;
  NodeDistanceMap node_distance;
  const Topology* space;
  typename Topology::point_type input_point;
  double distance_limit;
  double learning_constant;
  double falloff_ratio;

  template <typename Graph, typename Vertex>
  void examine_vertex(Vertex v, const Graph& /*unused*/) const {
    if (get(node_distance, v) > distance_limit) {
      throw over_distance_limit();
    }
    auto old_position = get(position_map, v);
    double distance = get(node_distance, v);
    double fraction = learning_constant * std::pow(falloff_ratio, distance * distance);
    put(position_map, v, space->move_position_toward(old_position, fraction, input_point));
  }
};

template <typename Graph, typename NodeDistanceMap, typename UpdatePosition, typename EdgeWeightMap>
void gursoy_shortest(const Graph& g, graph_vertex_descriptor_t<Graph> s, NodeDistanceMap node_distance,
                     UpdatePosition update_position, EdgeWeightMap weight) {
  if constexpr (std::is_same_v<EdgeWeightMap, single_property_map<double>>) {
    buffer_queue<graph_vertex_descriptor_t<Graph>> q;
    breadth_first_search(g, s, q, make_bfs_visitor(distance_recorder_on_tree_edge(node_distance), update_position),
                         two_bit_color_map(num_vertices(g), get(vertex_index, g)).ref());
  } else {
    dijkstra_shortest_paths(
        g, s, vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<Graph>::null_vertex()).ref(),
        vector_property_map(num_vertices(g), get(vertex_index, g), property_traits_value_t<NodeDistanceMap>{}).ref(),
        weight, get(vertex_index, g),
        make_dijkstra_visitor(distance_recorder_on_edge_relaxed(node_distance), update_position));
  }
}

}  // namespace gursoy_atun_detail

template <concepts::VertexListGraph G, typename Topology, concepts::ReadWriteVertexPropertyMap<G> PositionMap,
          concepts::ReadableVertexIndexMap<G> VertexIndexMap, concepts::ReadableEdgePropertyMap<G> EdgeWeightMap>
requires concepts::IncidenceGraph<G>
void gursoy_atun_step(const G& g, const Topology& space, PositionMap position, double diameter,
                      double learning_constant, VertexIndexMap vindex, EdgeWeightMap weight) {
  auto distance_from_input = vector_property_map(num_vertices(g), vindex, double{0.0});
  auto node_distance = vector_property_map(num_vertices(g), vindex, double{0.0});
  auto input_point = space.random_point();
  auto min_distance_loc = graph_traits<G>::null_vertex();
  double min_distance = 0.0;
  bool min_distance_unset = true;
  for (auto v : vertices(g)) {
    double this_distance = space.distance(get(position, v), input_point);
    put(distance_from_input, v, this_distance);
    if (min_distance_unset || this_distance < min_distance) {
      min_distance = this_distance;
      min_distance_loc = v;
    }
    min_distance_unset = false;
  }
  assert(!min_distance_unset);  // Graph must have at least one vertex
  try {
    gursoy_atun_detail::gursoy_shortest(
        g, min_distance_loc, node_distance.ref(),
        gursoy_atun_detail::update_position_visitor{position, node_distance.ref(), &space, input_point, diameter,
                                                    learning_constant, std::exp(-1. / (2 * diameter * diameter))},
        weight);
  } catch (const gursoy_atun_detail::over_distance_limit&) {
    // Thrown to break out of BFS or Dijkstra early
  }
}

template <concepts::VertexListGraph G, typename Topology, concepts::ReadWriteVertexPropertyMap<G> PositionMap,
          concepts::ReadableVertexIndexMap<G> VertexIndexMap, concepts::ReadableEdgePropertyMap<G> EdgeWeightMap>
requires concepts::IncidenceGraph<G>
void gursoy_atun_refine(const G& g, const Topology& space, PositionMap position, int nsteps, double diameter_initial,
                        double diameter_final, double learning_constant_initial, double learning_constant_final,
                        VertexIndexMap vindex, EdgeWeightMap weight) {
  double diameter_ratio = diameter_final / diameter_initial;
  double learning_constant_ratio = learning_constant_final / learning_constant_initial;
  for (int round = 0; round < nsteps; ++round) {
    double part_done = static_cast<double>(round) / (nsteps - 1);
    double diameter = diameter_initial * std::pow(diameter_ratio, part_done);
    double learning_constant = learning_constant_initial * std::pow(learning_constant_ratio, part_done);
    gursoy_atun_step(g, space, position, diameter, learning_constant, vindex, weight);
  }
}

template <concepts::VertexListGraph G, typename Topology, concepts::ReadWriteVertexPropertyMap<G> PositionMap,
          concepts::ReadableVertexIndexMap<G> VertexIndexMap, concepts::ReadableEdgePropertyMap<G> EdgeWeightMap>
requires concepts::IncidenceGraph<G>
void gursoy_atun_layout(const G& graph, const Topology& space, PositionMap position, int nsteps,
                        double diameter_initial, double diameter_final, double learning_constant_initial,
                        double learning_constant_final, VertexIndexMap vindex, EdgeWeightMap weight) {
  for (auto v : vertices(graph)) {
    put(position, v, space.random_point());
  }
  gursoy_atun_refine(graph, space, position, nsteps, diameter_initial, diameter_final, learning_constant_initial,
                     learning_constant_final, vindex, weight);
}

template <concepts::VertexListGraph G, typename Topology, concepts::ReadWriteVertexPropertyMap<G> PositionMap,
          concepts::ReadableVertexIndexMap<G> VertexIndexMap>
requires concepts::IncidenceGraph<G>
void gursoy_atun_layout(const G& g, const Topology& space, PositionMap position, int nsteps, double diameter_initial,
                        double diameter_final, double learning_constant_initial, double learning_constant_final,
                        VertexIndexMap vindex) {
  gursoy_atun_layout(g, space, position, nsteps, diameter_initial, diameter_final, learning_constant_initial,
                     learning_constant_final, vindex, single_property_map(double{1.0}));
}

template <concepts::VertexListGraph G, typename Topology, concepts::ReadWriteVertexPropertyMap<G> PositionMap>
requires concepts::IncidenceGraph<G>
void gursoy_atun_layout(const G& g, const Topology& space, PositionMap position, int nsteps, double diameter_initial,
                        double diameter_final = 1.0, double learning_constant_initial = 0.8,
                        double learning_constant_final = 0.2) {
  gursoy_atun_layout(g, space, position, nsteps, diameter_initial, diameter_final, learning_constant_initial,
                     learning_constant_final, get(vertex_index, g));
}

template <concepts::VertexListGraph G, typename Topology, concepts::ReadWriteVertexPropertyMap<G> PositionMap>
requires concepts::IncidenceGraph<G>
void gursoy_atun_layout(const G& g, const Topology& space, PositionMap position, int nsteps) {
  gursoy_atun_layout(g, space, position, nsteps, std::sqrt(static_cast<double>(num_vertices(g))));
}

template <concepts::VertexListGraph G, typename Topology, concepts::ReadWriteVertexPropertyMap<G> PositionMap>
requires concepts::IncidenceGraph<G>
void gursoy_atun_layout(const G& g, const Topology& space, PositionMap position) {
  gursoy_atun_layout(g, space, position, num_vertices(g));
}

template <concepts::VertexListGraph G, typename Topology, concepts::ReadWriteVertexPropertyMap<G> PositionMap,
          concepts::ReadableEdgePropertyMap<G> EdgeWeightMap>
requires concepts::IncidenceGraph<G>
void gursoy_atun_layout(const G& g, const Topology& space, PositionMap position, EdgeWeightMap weight) {
  gursoy_atun_layout(g, space, position, num_vertices(g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_GURSOY_ATUN_LAYOUT_H_
