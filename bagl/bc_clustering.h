// Copyright 2004 The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_BC_CLUSTERING_H_
#define BAGL_BAGL_BC_CLUSTERING_H_

#include <algorithm>
#include <vector>

#include "bagl/betweenness_centrality.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/property_map.h"

namespace bagl {

// Threshold termination function for the betweenness centrality
// clustering algorithm.
template <typename T>
struct bc_clustering_threshold {
  using centrality_type = T;

  // Terminate clustering when maximum absolute edge centrality is
  // below the given threshold.
  explicit bc_clustering_threshold(T threshold) : threshold_(threshold), dividend_(1.0) {}

  // Terminate clustering when the maximum edge centrality is below
  // the given threshold.
  // When normalize is true, the threshold is compared against the
  // normalized edge centrality based on the input graph; otherwise,
  // the threshold is compared against the absolute edge centrality.
  template <typename Graph>
  bc_clustering_threshold(T threshold, const Graph& g, bool normalize = true) : threshold_(threshold), dividend_(1.0) {
    if (normalize) {
      auto n = num_vertices(g);
      dividend_ = T((n - 1) * (n - 2)) / T(2);
    }
  }

  // Returns true when the given maximum edge centrality (potentially
  // normalized) falls below the threshold.
  template <typename Graph, typename Edge>
  bool operator()(T max_centrality, Edge /*unused*/, const Graph& /*unused*/) {
    return (max_centrality / dividend_) < threshold_;
  }

 protected:
  T threshold_;
  T dividend_;
};

/** Graph clustering based on edge betweenness centrality.
 *
 * This algorithm implements graph clustering based on edge
 * betweenness centrality. It is an iterative algorithm, where in each
 * step it compute the edge betweenness centrality (via @ref
 * brandes_betweenness_centrality) and removes the edge with the
 * maximum betweenness centrality. The @p done function object
 * determines when the algorithm terminates (the edge found when the
 * algorithm terminates will not be removed).
 *
 * @param g The graph on which clustering will be performed. The type
 * of this parameter (@c MutableGraph) must be a model of the
 * VertexListGraph, IncidenceGraph, EdgeListGraph, and Mutable Graph
 * concepts.
 *
 * @param done The function object that indicates termination of the
 * algorithm. It must be a ternary function object thats accepts the
 * maximum centrality, the descriptor of the edge that will be
 * removed, and the graph @p g.
 *
 * @param e_centrality (UTIL/OUT) The property map that will store
 * the betweenness centrality for each edge. When the algorithm
 * terminates, it will contain the edge centralities for the
 * graph. The type of this property map must model the
 * ReadWritePropertyMap concept. Defaults to an @c
 * iterator_property_map whose value type is
 * @c Done::centrality_type and using @c get(edge_index, g) for the
 * index map.
 *
 * @param v_index (IN) The property map that maps vertices to
 * indices in the range @c [0, num_vertices(g)). This type of this
 * property map must model the ReadablePropertyMap concept and its
 * value type must be an integral type. Defaults to
 * @c get(v_index, g).
 */
template <typename MutableGraph, typename Done, typename EdgeCentralityMap, typename VertexIndexMap>
void betweenness_centrality_clustering(MutableGraph& g, Done done, EdgeCentralityMap e_centrality,
                                       VertexIndexMap v_index) {
  if (has_no_edges(g)) {
    return;
  }

  // Function object that compares the centrality of edges
  auto cmp = [&e_centrality](auto e, auto f) { return get(e_centrality, e) < get(e_centrality, f); };

  bool is_done = false;
  do {
    brandes_betweenness_centrality(g, dummy_property_map(), e_centrality, v_index);
    auto edges_rg = edges(g);
    auto e = *std::max_element(edges_rg.begin(), edges_rg.end(), cmp);
    is_done = done(get(e_centrality, e), e, g);
    if (!is_done) {
      remove_edge(e, g);
    }
  } while (!is_done && !has_no_edges(g));
}

template <typename MutableGraph, typename Done, typename EdgeCentralityMap>
void betweenness_centrality_clustering(MutableGraph& g, Done done, EdgeCentralityMap e_centrality) {
  betweenness_centrality_clustering(g, done, e_centrality, get(vertex_index, g));
}

template <typename MutableGraph, typename Done>
void betweenness_centrality_clustering(MutableGraph& g, Done done) {
  using centrality_type = typename Done::centrality_type;
  betweenness_centrality_clustering(
      g, done, make_vector_property_map<centrality_type>(num_edges(g), get(edge_index, g)), get(vertex_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_BC_CLUSTERING_H_
