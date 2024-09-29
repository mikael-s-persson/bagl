// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DEGREE_CENTRALITY_H_
#define BAGL_BAGL_DEGREE_CENTRALITY_H_

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"

namespace bagl {

template <concepts::IncidenceGraph Graph>
struct influence_measure {
  auto operator()(graph_vertex_descriptor_t<Graph> v, const Graph& g) { return out_degree(v, g); }
};

template <concepts::BidirectionalGraph Graph>
struct prestige_measure {
  auto operator()(graph_vertex_descriptor_t<Graph> v, const Graph& g) { return in_degree(v, g); }
};

template <concepts::Graph Graph, concepts::DegreeMeasure<Graph> Measure>
auto degree_centrality(const Graph& g, graph_vertex_descriptor_t<Graph> v, Measure measure) {
  return measure(v, g);
}

template <concepts::IncidenceGraph Graph>
auto degree_centrality(const Graph& g, graph_vertex_descriptor_t<Graph> v) {
  return degree_centrality(g, v, influence_measure<Graph>());
}

// These are alias functions, intended to provide a more expressive interface.

template <concepts::IncidenceGraph Graph>
auto influence(const Graph& g, graph_vertex_descriptor_t<Graph> v) {
  return degree_centrality(g, v, influence_measure<Graph>());
}

template <concepts::BidirectionalGraph Graph>
auto prestige(const Graph& g, graph_vertex_descriptor_t<Graph> v) {
  return degree_centrality(g, v, prestige_measure<Graph>());
}

template <concepts::VertexListGraph Graph, concepts::WritableVertexPropertyMap<Graph> CentralityMap,
          concepts::DegreeMeasure<Graph> Measure>
void all_degree_centralities(const Graph& g, CentralityMap cent, Measure measure) {
  for (auto v : vertices(g)) {
    put(cent, v, degree_centrality(g, v, measure));
  }
}

template <concepts::VertexListGraph Graph, concepts::WritableVertexPropertyMap<Graph> CentralityMap>
requires concepts::IncidenceGraph<Graph>
void all_degree_centralities(const Graph& g, CentralityMap cent) {
  all_degree_centralities(g, cent, influence_measure<Graph>());
}

// More helper functions for computing influence and prestige.
// I hate the names of these functions, but influence and prestige
// don't pluralize too well.

template <concepts::VertexListGraph Graph, concepts::WritableVertexPropertyMap<Graph> CentralityMap>
requires concepts::IncidenceGraph<Graph>
void all_influence_values(const Graph& g, CentralityMap cent) {
  all_degree_centralities(g, cent, influence_measure<Graph>());
}

template <concepts::VertexListGraph Graph, concepts::WritableVertexPropertyMap<Graph> CentralityMap>
requires concepts::BidirectionalGraph<Graph>
void all_prestige_values(const Graph& g, CentralityMap cent) {
  all_degree_centralities(g, cent, prestige_measure<Graph>());
}

}  // namespace bagl

#endif  // BAGL_BAGL_DEGREE_CENTRALITY_H_
