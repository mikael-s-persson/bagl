// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//   04 April 2001: Added named parameter variant. (Jeremy Siek)
//   01 April 2001: Modified to use new <boost/limits.hpp> header. (JMaddock)
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DIJKSTRA_SHORTEST_PATHS_H_
#define BAGL_BAGL_DIJKSTRA_SHORTEST_PATHS_H_

#include <concepts>
#include <functional>
#include <limits>
#include <memory>
#include <type_traits>

#include "bagl/breadth_first_search.h"
#include "bagl/d_ary_heap.h"
#include "bagl/exception.h"
#include "bagl/graph_traits.h"
#include "bagl/mutable_queue.h"
#include "bagl/property_map.h"
#include "bagl/relax.h"
#include "bagl/two_bit_color_map.h"
#include "bagl/vector_property_map.h"

namespace bagl {

/**
 * @brief Updates a particular value in a queue used by Dijkstra's
 * algorithm.
 *
 * This routine is called by Dijkstra's algorithm after it has
 * decreased the distance from the seed vertex to the given @p
 * vertex. By default, this routine will just call @c
 * Q.update(vertex). However, other queues may provide more
 * specialized versions of this routine.
 *
 * @param Q             the queue that will be updated.
 * @param vertex        the vertex whose distance has been updated
 * @param old_distance  the previous distance to @p vertex
 */
template <typename Buffer, typename Vertex, typename DistanceType>
void dijkstra_queue_update(Buffer& q, Vertex vertex, DistanceType old_distance) {
  (void)old_distance;
  q.update(vertex);
}

namespace concepts {

template <typename Visitor, typename Graph>
concept DijkstraVisitor = std::copy_constructible<Visitor> &&
    requires(Visitor& vis, Graph& g, graph_vertex_descriptor_t<Graph> u, graph_edge_descriptor_t<Graph> e) {
  vis.initialize_vertex(u, g);
  vis.discover_vertex(u, g);
  vis.examine_vertex(u, g);
  vis.examine_edge(e, g);
  vis.edge_relaxed(e, g);
  vis.edge_not_relaxed(e, g);
  vis.finish_vertex(u, g);
};

}  // namespace concepts

template <class Visitors = null_visitor>
class dijkstra_visitor : public bfs_visitor<Visitors> {
 public:
  dijkstra_visitor() = default;
  explicit dijkstra_visitor(Visitors vis) : bfs_visitor<Visitors>(vis) {}

  template <class Edge, class Graph>
  void edge_relaxed(Edge e, Graph& g) {
    visitors_detail::invoke_edge_relaxed_on_all(this->vis_, e, g);
  }
  template <class Edge, class Graph>
  void edge_not_relaxed(Edge e, Graph& g) {
    visitors_detail::invoke_edge_not_relaxed_on_all(this->vis_, e, g);
  }

  template <class Edge, class Graph>
  void tree_edge(Edge u, Graph& g) = delete;
};

template <typename Visitors>
auto make_dijkstra_visitor(Visitors vis) {
  return dijkstra_visitor(std::move(vis));
}
using default_dijkstra_visitor = dijkstra_visitor<>;

namespace dijkstra_detail {

template <class UniformCostVisitor, class UpdatableQueue, class WeightMap, class PredecessorMap, class DistanceMap,
          class BinaryFunction, class BinaryPredicate>
struct dijkstra_bfs_visitor {
  using D = property_traits_value_t<DistanceMap>;
  using W = property_traits_value_t<WeightMap>;

  dijkstra_bfs_visitor(UniformCostVisitor vis, UpdatableQueue& q, WeightMap w, PredecessorMap p, DistanceMap d,
                       BinaryFunction combine, BinaryPredicate compare, D zero)
      : vis_(vis),
        q_(q),
        weight_(w),
        predecessor_(p),
        distance_(d),
        combine_(combine),
        compare_(compare),
        zero_(zero) {}

  template <class Edge, class Graph>
  void tree_edge(Edge e, Graph& g) {
    bool decreased = relax_target(e, g, weight_, predecessor_, distance_, combine_, compare_);
    if (decreased) {
      vis_.edge_relaxed(e, g);
    } else {
      vis_.edge_not_relaxed(e, g);
    }
  }
  template <class Edge, class Graph>
  void gray_target(Edge e, Graph& g) {
    D old_distance = get(distance_, target(e, g));

    bool decreased = relax_target(e, g, weight_, predecessor_, distance_, combine_, compare_);
    if (decreased) {
      dijkstra_queue_update(q_, target(e, g), old_distance);
      vis_.edge_relaxed(e, g);
    } else {
      vis_.edge_not_relaxed(e, g);
    }
  }

  template <class Vertex, class Graph>
  void initialize_vertex(Vertex u, Graph& g) {
    vis_.initialize_vertex(u, g);
  }
  template <class Edge, class Graph>
  void non_tree_edge(Edge, Graph&) {}
  template <class Vertex, class Graph>
  void discover_vertex(Vertex u, Graph& g) {
    vis_.discover_vertex(u, g);
  }
  template <class Vertex, class Graph>
  void examine_vertex(Vertex u, Graph& g) {
    vis_.examine_vertex(u, g);
  }
  template <class Edge, class Graph>
  void examine_edge(Edge e, Graph& g) {
    // Test for negative-weight edges:
    //
    // Reasons that other comparisons do not work:
    //
    // compare_(e_weight, D(0)):
    //    compare_ only needs to work on distances, not weights, and
    //    those types do not need to be the same (bug 8398,
    //    https://svn.boost.org/trac/boost/ticket/8398).
    // compare_(combine_(source_dist, e_weight), source_dist):
    //    if combine_ is project2nd (as in prim_minimum_spanning_tree),
    //    this test will claim that the edge weight is negative whenever
    //    the edge weight is less than source_dist, even if both of
    //    those are positive (bug 9012,
    //    https://svn.boost.org/trac/boost/ticket/9012).
    // compare_(combine_(e_weight, source_dist), source_dist):
    //    would fix project2nd issue, but documentation only requires
    //    that combine_ be able to take a distance and a weight (in
    //    that order) and return a distance.

    // W e_weight = get(m_weight, e);
    // sd_plus_ew = source_dist + e_weight.
    // D sd_plus_ew = combine_(source_dist, e_weight);
    // sd_plus_2ew = source_dist + 2 * e_weight.
    // D sd_plus_2ew = combine_(sd_plus_ew, e_weight);
    // The test here is equivalent to e_weight < 0 if combine_ has a
    // cancellation law, but always returns false when combine_ is a
    // projection operator.
    if (compare_(combine_(zero_, get(weight_, e)), zero_)) {
      throw negative_edge();
    }
    // End of test for negative-weight edges.

    vis_.examine_edge(e, g);
  }
  template <class Edge, class Graph>
  void black_target(Edge, Graph&) {}
  template <class Vertex, class Graph>
  void finish_vertex(Vertex u, Graph& g) {
    vis_.finish_vertex(u, g);
  }

  UniformCostVisitor vis_;
  UpdatableQueue& q_;
  WeightMap weight_;
  PredecessorMap predecessor_;
  DistanceMap distance_;
  BinaryFunction combine_;
  BinaryPredicate compare_;
  D zero_;
};

template <concepts::IncidenceGraph G, class IndexMap, class Value>
auto make_default_index_map(const G& g, const IndexMap& index) {
  if constexpr (is_vertex_list_graph_v<G>) {
    return make_vector_property_map<Value>(num_vertices(g), index);
  } else {
    return make_vector_property_map<Value>(index);
  }
}

template <concepts::IncidenceGraph G, class IndexMap>
auto make_default_color_map(const G& g, const IndexMap& index) {
  if constexpr (is_vertex_list_graph_v<G>) {
    return two_bit_color_map(num_vertices(g), index);
  } else {
    return make_vector_property_map<two_bit_color_type>(index);
  }
}
}  // namespace dijkstra_detail

// Call breadth first search with default color map.
template <concepts::IncidenceGraph G, std::ranges::input_range Seeds, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, class IndexMap, class Compare,
          class Combine>
void dijkstra_shortest_paths_no_init(const G& g, Seeds seeds, PredecessorMap predecessor, DistanceMap distance,
                                     WeightMap weight, IndexMap index_map, Compare compare, Combine combine,
                                     property_traits_value_t<DistanceMap> zero, V vis) {
  dijkstra_shortest_paths_no_init(g, seeds, predecessor, distance, weight, index_map, compare, combine, zero, vis,
                                  dijkstra_detail::make_default_color_map(g, index_map));
}

// Call breadth first search with default color map.
template <concepts::IncidenceGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, class IndexMap, class Compare,
          class Combine>
void dijkstra_shortest_paths_no_init(const G& g, graph_vertex_descriptor_t<G> s, PredecessorMap predecessor,
                                     DistanceMap distance, WeightMap weight, IndexMap index_map, Compare compare,
                                     Combine combine, property_traits_value_t<DistanceMap> zero, V vis) {
  dijkstra_shortest_paths_no_init(g, std::ranges::single_view(s), predecessor, distance, weight, index_map, compare,
                                  combine, zero, vis);
}

// Call breadth first search
template <concepts::IncidenceGraph G, std::ranges::input_range Seeds, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, class IndexMap, class Compare,
          class Combine, concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void dijkstra_shortest_paths_no_init(const G& g, Seeds seeds, PredecessorMap predecessor, DistanceMap distance,
                                     WeightMap weight, IndexMap index_map, Compare compare, Combine combine,
                                     property_traits_value_t<DistanceMap> zero, V vis, ColorMap color) {
  using Vertex = graph_vertex_descriptor_t<G>;

  // Now the default: use a d-ary heap
  auto index_in_heap = dijkstra_detail::make_default_index_map(g, index_map);
  auto q = make_d_ary_heap_indirect<Vertex, 4>(distance, index_in_heap, compare);

  dijkstra_detail::dijkstra_bfs_visitor bfs_vis{vis, q, weight, predecessor, distance, combine, compare, zero};

  breadth_first_visit(g, seeds, q, bfs_vis, color);
}

// Call breadth first search
template <concepts::IncidenceGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, class IndexMap, class Compare,
          class Combine, concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void dijkstra_shortest_paths_no_init(const G& g, graph_vertex_descriptor_t<G> s, PredecessorMap predecessor,
                                     DistanceMap distance, WeightMap weight, IndexMap index_map, Compare compare,
                                     Combine combine, property_traits_value_t<DistanceMap> zero, V vis,
                                     ColorMap color) {
  dijkstra_shortest_paths_no_init(g, std::ranges::single_view(s), predecessor, distance, weight, index_map, compare,
                                  combine, zero, vis, color);
}

// Initialize distances and call breadth first search with default color map
template <concepts::VertexListGraph G, std::ranges::input_range Seeds, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, class IndexMap, class Compare,
          class Combine>
void dijkstra_shortest_paths(const G& g, Seeds seeds, PredecessorMap predecessor, DistanceMap distance,
                             WeightMap weight, IndexMap index_map, Compare compare, Combine combine,
                             property_traits_value_t<DistanceMap> inf, property_traits_value_t<DistanceMap> zero,
                             V vis) {
  two_bit_color_map<IndexMap> color(num_vertices(g), index_map);
  dijkstra_shortest_paths(g, seeds, predecessor, distance, weight, index_map, compare, combine, inf, zero, vis, color);
}

// Initialize distances and call breadth first search with default color map
template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, class IndexMap, class Compare,
          class Combine>
void dijkstra_shortest_paths(const G& g, graph_vertex_descriptor_t<G> s, PredecessorMap predecessor,
                             DistanceMap distance, WeightMap weight, IndexMap index_map, Compare compare,
                             Combine combine, property_traits_value_t<DistanceMap> inf,
                             property_traits_value_t<DistanceMap> zero, V vis) {
  dijkstra_shortest_paths(g, std::ranges::single_view(s), predecessor, distance, weight, index_map, compare, combine,
                          inf, zero, vis);
}

// Initialize distances and call breadth first search
template <concepts::VertexListGraph G, std::ranges::input_range Seeds, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, class IndexMap, class Compare,
          class Combine, concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void dijkstra_shortest_paths(const G& g, Seeds seeds, PredecessorMap predecessor, DistanceMap distance,
                             WeightMap weight, IndexMap index_map, Compare compare, Combine combine,
                             property_traits_value_t<DistanceMap> inf, property_traits_value_t<DistanceMap> zero, V vis,
                             ColorMap color) {
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;
  for (auto u : vertices(g)) {
    vis.initialize_vertex(u, g);
    put(distance, u, inf);
    put(predecessor, u, u);
    put(color, u, Color::white());
  }
  for (auto u : seeds) {
    put(distance, u, zero);
  }

  dijkstra_shortest_paths_no_init(g, seeds, predecessor, distance, weight, index_map, compare, combine, zero, vis,
                                  color);
}

// Initialize distances and call breadth first search
template <concepts::VertexListGraph G, concepts::DijkstraVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, class IndexMap, class Compare,
          class Combine, concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void dijkstra_shortest_paths(const G& g, graph_vertex_descriptor_t<G> s, PredecessorMap predecessor,
                             DistanceMap distance, WeightMap weight, IndexMap index_map, Compare compare,
                             Combine combine, property_traits_value_t<DistanceMap> inf,
                             property_traits_value_t<DistanceMap> zero, V vis, ColorMap color) {
  dijkstra_shortest_paths(g, std::ranges::single_view(s), predecessor, distance, weight, index_map, compare, combine,
                          inf, zero, vis, color);
}

}  // namespace bagl

#endif  // BAGL_BAGL_DIJKSTRA_SHORTEST_PATHS_H_
