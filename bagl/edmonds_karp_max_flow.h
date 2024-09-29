// Copyright 2000 University of Notre Dame.
// Authors: Jeremy G. Siek, Andrew Lumsdaine, Lie-Quan Lee
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_EDMONDS_KARP_MAX_FLOW_H_
#define BAGL_BAGL_EDMONDS_KARP_MAX_FLOW_H_

#include <algorithm>  // for std::min and std::max
#include <queue>
#include <vector>

#include "bagl/breadth_first_search.h"
#include "bagl/detail/augment.h"
#include "bagl/filtered_graph.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"

namespace bagl {

// The "labeling" algorithm from "Network Flows" by Ahuja, Magnanti,
// Orlin.  I think this is the same as or very similar to the original
// Edmonds-Karp algorithm.  This solves the maximum flow problem.

template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> CapacityEdgeMap,
          concepts::ReadWriteEdgePropertyMap<G> ResidualCapacityEdgeMap,
          concepts::ReadableEdgePropertyMap<G> ReverseEdgeMap, concepts::ReadWriteVertexPropertyMap<G> ColorMap,
          concepts::ReadWriteVertexPropertyMap<G> PredEdgeMap>
requires std::same_as<graph_edge_descriptor_t<G>, property_traits_value_t<PredEdgeMap>> &&
    std::convertible_to<property_traits_value_t<CapacityEdgeMap>, property_traits_value_t<ResidualCapacityEdgeMap>>
auto edmonds_karp_max_flow(const G& g, graph_vertex_descriptor_t<G> src, graph_vertex_descriptor_t<G> sink,
                           CapacityEdgeMap cap,          // get(edge_capacity, g)
                           ResidualCapacityEdgeMap res,  // get(edge_residual_capacity, g)
                           ReverseEdgeMap rev,           // get(edge_reverse, g)
                           ColorMap color,               // get(vertex_color, g)
                           PredEdgeMap pred) {
  using ColorValue = typename property_traits<ColorMap>::value_type;
  using Color = color_traits<ColorValue>;

  for (auto u : vertices(g)) {
    for (auto e : out_edges(u, g)) {
      put(res, e, get(cap, e));
    }
  }

  put(color, sink, Color::gray());
  while (get(color, sink) != Color::white()) {
    std::queue<graph_vertex_descriptor_t<G>> q;
    breadth_first_search(augment_detail::residual_graph(g, res), src, q,
                         make_bfs_visitor(edge_predecessor_recorder_on_tree_edge(pred)), color);
    if (get(color, sink) != Color::white()) {
      augment_detail::augment(g, src, sink, pred, res, rev);
    }
  }  // while

  property_traits_value_t<CapacityEdgeMap> flow = 0;
  for (auto e : out_edges(src, g)) {
    flow += (get(cap, e) - get(res, e));
  }
  return flow;
}  // edmonds_karp_max_flow()

template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> CapacityEdgeMap,
          concepts::ReadWriteEdgePropertyMap<G> ResidualCapacityEdgeMap,
          concepts::ReadableEdgePropertyMap<G> ReverseEdgeMap, concepts::ReadWriteVertexPropertyMap<G> PredEdgeMap>
requires std::same_as<graph_edge_descriptor_t<G>, property_traits_value_t<PredEdgeMap>> &&
    std::convertible_to<property_traits_value_t<CapacityEdgeMap>, property_traits_value_t<ResidualCapacityEdgeMap>>
auto edmonds_karp_max_flow(const G& g, graph_vertex_descriptor_t<G> src, graph_vertex_descriptor_t<G> sink,
                           CapacityEdgeMap cap,          // get(edge_capacity, g)
                           ResidualCapacityEdgeMap res,  // get(edge_residual_capacity, g)
                           ReverseEdgeMap rev,           // get(edge_reverse, g)
                           PredEdgeMap pred) {
  return edmonds_karp_max_flow(g, src, sink, cap, res, rev,
                               two_bit_color_map(num_vertices_or_zero(g), get(vertex_index, g)), pred);
}

template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> CapacityEdgeMap,
          concepts::ReadableEdgePropertyMap<G> ReverseEdgeMap, concepts::ReadWriteVertexPropertyMap<G> ColorMap,
          concepts::ReadWriteVertexPropertyMap<G> PredEdgeMap>
requires std::same_as<graph_edge_descriptor_t<G>, property_traits_value_t<PredEdgeMap>>
auto edmonds_karp_max_flow(const G& g, graph_vertex_descriptor_t<G> src, graph_vertex_descriptor_t<G> sink,
                           CapacityEdgeMap cap,  // get(edge_capacity, g)
                           ReverseEdgeMap rev,   // get(edge_reverse, g)
                           ColorMap color,       // get(vertex_color, g)
                           PredEdgeMap pred) {
  using FlowValue = property_traits_value_t<CapacityEdgeMap>;
  return edmonds_karp_max_flow(g, src, sink, cap,
                               vector_property_map(num_vertices_or_zero(g), get(vertex_index, g), FlowValue{}), rev,
                               color, pred);
}

template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> CapacityEdgeMap,
          concepts::ReadableEdgePropertyMap<G> ReverseEdgeMap, concepts::ReadWriteVertexPropertyMap<G> PredEdgeMap>
requires std::same_as<graph_edge_descriptor_t<G>, property_traits_value_t<PredEdgeMap>>
auto edmonds_karp_max_flow(const G& g, graph_vertex_descriptor_t<G> src, graph_vertex_descriptor_t<G> sink,
                           CapacityEdgeMap cap,  // get(edge_capacity, g)
                           ReverseEdgeMap rev,   // get(edge_reverse, g)
                           PredEdgeMap pred) {
  using FlowValue = property_traits_value_t<CapacityEdgeMap>;
  return edmonds_karp_max_flow(g, src, sink, cap,
                               vector_property_map(num_vertices_or_zero(g), get(vertex_index, g), FlowValue{}), rev,
                               two_bit_color_map(num_vertices_or_zero(g), get(vertex_index, g)), pred);
}

template <concepts::IncidenceGraph G, concepts::ReadableVertexPropertyMap<G> VertexIndexMap,
          concepts::ReadableEdgePropertyMap<G> CapacityEdgeMap, concepts::ReadableEdgePropertyMap<G> ReverseEdgeMap,
          concepts::ReadWriteVertexPropertyMap<G> PredEdgeMap>
requires std::same_as<graph_edge_descriptor_t<G>, property_traits_value_t<PredEdgeMap>>
auto edmonds_karp_max_flow(const G& g, graph_vertex_descriptor_t<G> src, graph_vertex_descriptor_t<G> sink,
                           VertexIndexMap v_index,  // get(vertex_index, g)
                           CapacityEdgeMap cap,     // get(edge_capacity, g)
                           ReverseEdgeMap rev,      // get(edge_reverse, g)
                           PredEdgeMap pred) {
  using FlowValue = property_traits_value_t<CapacityEdgeMap>;
  return edmonds_karp_max_flow(g, src, sink, cap, vector_property_map(num_vertices_or_zero(g), v_index, FlowValue{}),
                               rev, two_bit_color_map(num_vertices_or_zero(g), v_index), pred);
}

template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> CapacityEdgeMap,
          concepts::ReadableEdgePropertyMap<G> ReverseEdgeMap>
auto edmonds_karp_max_flow(const G& g, graph_vertex_descriptor_t<G> src, graph_vertex_descriptor_t<G> sink,
                           CapacityEdgeMap cap,   // get(edge_capacity, g)
                           ReverseEdgeMap rev) {  // get(edge_reverse, g)
  using FlowValue = property_traits_value_t<CapacityEdgeMap>;
  using Edge = graph_edge_descriptor_t<G>;
  return edmonds_karp_max_flow(g, src, sink, cap,
                               vector_property_map(num_vertices_or_zero(g), get(vertex_index, g), FlowValue{}), rev,
                               two_bit_color_map(num_vertices_or_zero(g), get(vertex_index, g)),
                               vector_property_map(num_vertices_or_zero(g), get(vertex_index, g), Edge{}));
}

template <concepts::IncidenceGraph G, concepts::ReadableVertexPropertyMap<G> VertexIndexMap,
          concepts::ReadableEdgePropertyMap<G> CapacityEdgeMap, concepts::ReadableEdgePropertyMap<G> ReverseEdgeMap>
auto edmonds_karp_max_flow(const G& g, graph_vertex_descriptor_t<G> src, graph_vertex_descriptor_t<G> sink,
                           VertexIndexMap v_index,  // get(vertex_index, g)
                           CapacityEdgeMap cap,     // get(edge_capacity, g)
                           ReverseEdgeMap rev) {    // get(edge_reverse, g)
  using FlowValue = property_traits_value_t<CapacityEdgeMap>;
  using Edge = graph_edge_descriptor_t<G>;
  return edmonds_karp_max_flow(g, src, sink, cap, vector_property_map(num_vertices_or_zero(g), v_index, FlowValue{}),
                               rev, two_bit_color_map(num_vertices_or_zero(g), v_index),
                               vector_property_map(num_vertices_or_zero(g), v_index, Edge{}));
}

}  // namespace bagl

#endif  // BAGL_BAGL_EDMONDS_KARP_MAX_FLOW_H_
