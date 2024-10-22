// Copyright 2002 Indiana University.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DAG_SHORTEST_PATHS_H_
#define BAGL_BAGL_DAG_SHORTEST_PATHS_H_

#include "bagl/dijkstra_shortest_paths.h"
#include "bagl/property_map.h"
#include "bagl/topological_sort.h"
#include "bagl/vector_property_map.h"

// single-source shortest paths for a Directed Acyclic Graph (DAG)

namespace bagl {

// Initalize distances and call depth first search
template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::PropertyComparator<DistanceMap> Compare, concepts::PropertyCombinator<DistanceMap> Combine>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>> &&
    std::convertible_to<property_traits_value_t<WeightMap>, property_traits_value_t<DistanceMap>>
void dag_shortest_paths(const G& g, graph_vertex_descriptor_t<G> start,
                        DistanceMap distance,  // get(vertex_distance, g)
                        WeightMap weight,      // get(edge_weight, g)
                        ColorMap color,        // get(vertex_color, g) or
                                               // two_bit_color_map(num_vertices(g), get(vertex_index, g)).ref()
                        PredecessorMap pred,   // get(vertex_predecessor, g)
                        V vis,                 // default_dfs_visitor()
                        Compare compare, Combine combine, property_traits_value_t<DistanceMap> inf,
                        property_traits_value_t<DistanceMap> zero) {
  std::vector<graph_vertex_descriptor_t<G>> rev_topo_order;
  rev_topo_order.reserve(num_vertices(g));

  // Call 'depth_first_visit', not 'topological_sort', because we don't
  // want to traverse the entire graph, only vertices reachable from 'start',
  // and 'topological_sort' will traverse everything. The logic below
  // is the same as for 'topological_sort', only we call 'depth_first_visit'
  // and 'topological_sort' calls 'depth_first_search'.
  depth_first_visit(g, start, topo_sort_visitor(std::back_inserter(rev_topo_order)), color);

  for (auto u : vertices(g)) {
    put(distance, u, inf);
    put(pred, u, u);
  }

  put(distance, start, zero);
  vis.discover_vertex(start, g);
  for (auto u : rev_topo_order | std::views::reverse) {
    vis.examine_vertex(u, g);
    for (auto e : out_edges(u, g)) {
      vis.discover_vertex(target(e, g), g);
      bool decreased = relax(e, g, weight, pred, distance, combine, compare);
      if (decreased) {
        vis.edge_relaxed(e, g);
      } else {
        vis.edge_not_relaxed(e, g);
      }
    }
    vis.finish_vertex(u, g);
  }
}

template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, concepts::ReadWriteVertexPropertyMap<G> PredecessorMap>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>> &&
    std::convertible_to<property_traits_value_t<WeightMap>, property_traits_value_t<DistanceMap>>
void dag_shortest_paths(const G& g, graph_vertex_descriptor_t<G> start, DistanceMap distance, WeightMap weight,
                        ColorMap color, PredecessorMap pred, V vis) {
  using D = property_traits_value_t<WeightMap>;
  dag_shortest_paths(g, start, distance, weight, color, pred, vis, std::less<>(), closed_plus<>(), default_inf_v<D>,
                     default_zero_v<D>);
}

template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, concepts::ReadWriteVertexPropertyMap<G> PredecessorMap>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>>
void dag_shortest_paths(const G& g, graph_vertex_descriptor_t<G> start, WeightMap weight, ColorMap color,
                        PredecessorMap pred, V vis) {
  using D = property_traits_value_t<WeightMap>;
  dag_shortest_paths(g, start, vector_property_map(num_vertices(g), get(vertex_index, g), default_zero_v<D>).ref(),
                     weight, color, pred, vis);
}

template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>> &&
    std::convertible_to<property_traits_value_t<WeightMap>, property_traits_value_t<DistanceMap>>
void dag_shortest_paths(const G& g, graph_vertex_descriptor_t<G> start, DistanceMap distance, WeightMap weight,
                        PredecessorMap pred, V vis) {
  dag_shortest_paths(g, start, distance, weight, two_bit_color_map(num_vertices(g), get(vertex_index, g)).ref(), pred,
                     vis);
}

template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>>
void dag_shortest_paths(const G& g, graph_vertex_descriptor_t<G> start, WeightMap weight, PredecessorMap pred, V vis) {
  using D = property_traits_value_t<WeightMap>;
  dag_shortest_paths(g, start, vector_property_map(num_vertices(g), get(vertex_index, g), default_zero_v<D>).ref(),
                     weight, two_bit_color_map(num_vertices(g), get(vertex_index, g)).ref(), pred, vis);
}

template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>>
void dag_shortest_paths(const G& g, graph_vertex_descriptor_t<G> start, WeightMap weight, PredecessorMap pred, V vis,
                        VertexIndexMap v_index) {
  using D = property_traits_value_t<WeightMap>;
  dag_shortest_paths(g, start, vector_property_map(num_vertices(g), v_index, default_zero_v<D>).ref(), weight,
                     two_bit_color_map(num_vertices(g), v_index).ref(), pred, vis);
}

}  // namespace bagl

#endif  // BAGL_BAGL_DAG_SHORTEST_PATHS_H_
