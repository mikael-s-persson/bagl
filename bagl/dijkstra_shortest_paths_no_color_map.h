// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Copyright 2009 Trustees of Indiana University.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek, Michael Hansen
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DIJKSTRA_SHORTEST_PATHS_NO_COLOR_MAP_H_
#define BAGL_BAGL_DIJKSTRA_SHORTEST_PATHS_NO_COLOR_MAP_H_

#include "bagl/d_ary_heap.h"
#include "bagl/dijkstra_shortest_paths.h"
#include "bagl/graph_traits.h"
#include "bagl/relax.h"
#include "bagl/vector_property_map.h"

namespace bagl {

// No init version
template <concepts::IncidenceGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, concepts::ReadableVertexPropertyMap<G> VertexIndexMap,
          concepts::PropertyComparator<DistanceMap> Compare, concepts::PropertyCombinator<DistanceMap> Combine>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>> &&
    std::convertible_to<property_traits_value_t<WeightMap>, property_traits_value_t<DistanceMap>>
void dijkstra_shortest_paths_no_color_map_no_init(const G& g, graph_vertex_descriptor_t<G> start,
                                                  PredecessorMap predecessor_map,  // get(vertex_predecessor, g)
                                                  DistanceMap distance_map,        // get(vertex_distance, g)
                                                  WeightMap weight_map,            // get(edge_weight, g)
                                                  VertexIndexMap index_map,        // get(vertex_index, g)
                                                  Compare compare, Combine combine,
                                                  property_traits_value_t<DistanceMap> d_infinity, V visitor) {
  // Default - use d-ary heap (d = 4)
  auto vertex_queue = make_d_ary_heap_indirect<graph_vertex_descriptor_t<G>, 4>(
      distance_map, vector_property_map{num_vertices_or_zero(g), index_map, std::size_t{}}, compare);

  // Add vertex to the queue
  vertex_queue.push(start);

  // Starting vertex will always be the first discovered vertex
  visitor.discover_vertex(start, g);

  while (!vertex_queue.empty()) {
    auto min_vertex = vertex_queue.top();
    vertex_queue.pop();

    visitor.examine_vertex(min_vertex, g);

    // Check if any other vertices can be reached
    auto min_vertex_distance = get(distance_map, min_vertex);

    if (!distance_compare(min_vertex_distance, d_infinity)) {
      // This is the minimum vertex, so all other vertices are unreachable
      return;
    }

    // Examine neighbors of min_vertex
    for (auto current_edge : out_edges(min_vertex, g)) {
      visitor.examine_edge(current_edge, g);

      // Check if the edge has a negative weight
      if (distance_compare(get(weight_map, current_edge), d_infinity)) {
        throw negative_edge();
      }

      // Extract the neighboring vertex and get its distance
      auto neighbor_vertex = target(current_edge, g);
      auto neighbor_vertex_distance = get(distance_map, neighbor_vertex);
      bool is_neighbor_undiscovered = !compare(neighbor_vertex_distance, d_infinity);

      // Attempt to relax the edge
      bool was_edge_relaxed =
          relax_target(current_edge, g, weight_map, predecessor_map, distance_map, combine, compare);

      if (was_edge_relaxed) {
        visitor.edge_relaxed(current_edge, g);
        if (is_neighbor_undiscovered) {
          visitor.discover_vertex(neighbor_vertex, g);
          vertex_queue.push(neighbor_vertex);
        } else {
          vertex_queue.update(neighbor_vertex);
        }
      } else {
        visitor.edge_not_relaxed(current_edge, g);
      }

    }  // end out edge iteration

    visitor.finish_vertex(min_vertex, g);
  }  // end while queue not empty
}

template <concepts::IncidenceGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>> &&
    std::convertible_to<property_traits_value_t<WeightMap>, property_traits_value_t<DistanceMap>>
void dijkstra_shortest_paths_no_color_map_no_init(const G& g, graph_vertex_descriptor_t<G> start,
                                                  PredecessorMap predecessor_map,  // get(vertex_predecessor, g)
                                                  DistanceMap distance_map,        // get(vertex_distance, g)
                                                  WeightMap weight_map,            // get(edge_weight, g)
                                                  VertexIndexMap index_map,        // get(vertex_index, g)
                                                  V vis) {
  using D = property_traits_value_t<WeightMap>;
  dijkstra_shortest_paths_no_color_map_no_init(g, start, predecessor_map, distance_map, weight_map, index_map,
                                               std::less<>(), closed_plus<>(), default_inf_v<D>, vis);
}

template <concepts::IncidenceGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>>
void dijkstra_shortest_paths_no_color_map_no_init(const G& g, graph_vertex_descriptor_t<G> start,
                                                  PredecessorMap predecessor_map,  // get(vertex_predecessor, g)
                                                  WeightMap weight_map,            // get(edge_weight, g)
                                                  VertexIndexMap index_map,        // get(vertex_index, g)
                                                  V vis) {
  using D = property_traits_value_t<WeightMap>;
  dijkstra_shortest_paths_no_color_map_no_init(
      g, start, predecessor_map, make_vector_property_map(num_vertices_or_zero(g), index_map, default_zero_v<D>),
      weight_map, index_map, std::less<>(), closed_plus<>(), default_inf_v<D>, vis);
}

template <concepts::IncidenceGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>> &&
    std::convertible_to<property_traits_value_t<WeightMap>, property_traits_value_t<DistanceMap>>
void dijkstra_shortest_paths_no_color_map_no_init(const G& g, graph_vertex_descriptor_t<G> start,
                                                  PredecessorMap predecessor_map,  // get(vertex_predecessor, g)
                                                  DistanceMap distance_map,        // get(vertex_distance, g)
                                                  WeightMap weight_map,            // get(edge_weight, g)
                                                  V vis) {
  using D = property_traits_value_t<WeightMap>;
  dijkstra_shortest_paths_no_color_map_no_init(g, start, predecessor_map, distance_map, weight_map,
                                               get(vertex_index, g), std::less<>(), closed_plus<>(), default_inf_v<D>,
                                               vis);
}

template <concepts::IncidenceGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadableEdgePropertyMap<G> WeightMap>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>>
void dijkstra_shortest_paths_no_color_map_no_init(const G& g, graph_vertex_descriptor_t<G> start,
                                                  PredecessorMap predecessor_map,  // get(vertex_predecessor, g)
                                                  WeightMap weight_map,            // get(edge_weight, g)
                                                  V vis) {
  using D = property_traits_value_t<WeightMap>;
  dijkstra_shortest_paths_no_color_map_no_init(
      g, start, predecessor_map,
      make_vector_property_map(num_vertices_or_zero(g), get(vertex_index, g), default_zero_v<D>), weight_map,
      get(vertex_index, g), std::less<>(), closed_plus<>(), default_inf_v<D>, vis);
}

// Full init version
template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, concepts::ReadableVertexPropertyMap<G> VertexIndexMap,
          concepts::PropertyComparator<DistanceMap> Compare, concepts::PropertyCombinator<DistanceMap> Combine>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>> &&
    std::convertible_to<property_traits_value_t<WeightMap>, property_traits_value_t<DistanceMap>>
void dijkstra_shortest_paths_no_color_map(const G& g, graph_vertex_descriptor_t<G> start,
                                          PredecessorMap predecessor_map,  // get(vertex_predecessor, g)
                                          DistanceMap distance_map,        // get(vertex_distance, g)
                                          WeightMap weight_map,            // get(edge_weight, g)
                                          VertexIndexMap index_map,        // get(vertex_index, g)
                                          Compare compare, Combine combine,
                                          property_traits_value_t<DistanceMap> d_infinity,
                                          property_traits_value_t<DistanceMap> d_zero, V visitor) {
  // Initialize vertices
  for (auto current_vertex : vertices(g)) {
    visitor.initialize_vertex(current_vertex, g);

    // Default all distances to infinity
    put(distance_map, current_vertex, d_infinity);

    // Default all vertex predecessors to the vertex itself
    put(predecessor_map, current_vertex, current_vertex);
  }

  // Set distance for start to zero
  put(distance_map, start, d_zero);

  // Pass everything on to the no_init version
  dijkstra_shortest_paths_no_color_map_no_init(g, start, predecessor_map, distance_map, weight_map, index_map, combine,
                                               compare, d_infinity, visitor);
}

template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>> &&
    std::convertible_to<property_traits_value_t<WeightMap>, property_traits_value_t<DistanceMap>>
void dijkstra_shortest_paths_no_color_map(const G& g, graph_vertex_descriptor_t<G> start,
                                          PredecessorMap predecessor_map,  // get(vertex_predecessor, g)
                                          DistanceMap distance_map,        // get(vertex_distance, g)
                                          WeightMap weight_map,            // get(edge_weight, g)
                                          VertexIndexMap index_map,        // get(vertex_index, g)
                                          V vis) {
  using D = property_traits_value_t<WeightMap>;
  dijkstra_shortest_paths_no_color_map(g, start, predecessor_map, distance_map, weight_map, index_map, std::less<>(),
                                       closed_plus<>(), default_inf_v<D>, default_zero_v<D>, vis);
}

template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>>
void dijkstra_shortest_paths_no_color_map(const G& g, graph_vertex_descriptor_t<G> start,
                                          PredecessorMap predecessor_map,  // get(vertex_predecessor, g)
                                          WeightMap weight_map,            // get(edge_weight, g)
                                          VertexIndexMap index_map,        // get(vertex_index, g)
                                          V vis) {
  using D = property_traits_value_t<WeightMap>;
  dijkstra_shortest_paths_no_color_map(
      g, start, predecessor_map, make_vector_property_map(num_vertices_or_zero(g), index_map, default_zero_v<D>),
      weight_map, index_map, std::less<>(), closed_plus<>(), default_inf_v<D>, default_zero_v<D>, vis);
}

template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>> &&
    std::convertible_to<property_traits_value_t<WeightMap>, property_traits_value_t<DistanceMap>>
void dijkstra_shortest_paths_no_color_map(const G& g, graph_vertex_descriptor_t<G> start,
                                          PredecessorMap predecessor_map,  // get(vertex_predecessor, g)
                                          DistanceMap distance_map,        // get(vertex_distance, g)
                                          WeightMap weight_map,            // get(edge_weight, g)
                                          V vis) {
  using D = property_traits_value_t<WeightMap>;
  dijkstra_shortest_paths_no_color_map(g, start, predecessor_map, distance_map, weight_map, get(vertex_index, g),
                                       std::less<>(), closed_plus<>(), default_inf_v<D>, default_zero_v<D>, vis);
}

template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadableEdgePropertyMap<G> WeightMap>
requires std::same_as<graph_vertex_descriptor_t<G>, property_traits_value_t<PredecessorMap>>
void dijkstra_shortest_paths_no_color_map(const G& g, graph_vertex_descriptor_t<G> start,
                                          PredecessorMap predecessor_map,  // get(vertex_predecessor, g)
                                          WeightMap weight_map,            // get(edge_weight, g)
                                          V vis) {
  using D = property_traits_value_t<WeightMap>;
  dijkstra_shortest_paths_no_color_map(
      g, start, predecessor_map,
      make_vector_property_map(num_vertices_or_zero(g), get(vertex_index, g), default_zero_v<D>), weight_map,
      get(vertex_index, g), std::less<>(), closed_plus<>(), default_inf_v<D>, default_zero_v<D>, vis);
}

}  // namespace bagl

#endif  // BAGL_BAGL_DIJKSTRA_SHORTEST_PATHS_NO_COLOR_MAP_H_
