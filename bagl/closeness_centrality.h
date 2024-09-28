// Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_CLOSENESS_CENTRALITY_H_
#define BAGL_BAGL_CLOSENESS_CENTRALITY_H_

#include "bagl/exterior_property.h"
#include "bagl/geodesic.h"
#include "bagl/graph_concepts.h"
#include "bagl/properties.h"

namespace bagl {

template <concepts::NumericValue DistanceType, concepts::NumericValue ResultType,
          std::invocable<ResultType> Reciprocal = reciprocal>
struct closeness_measure : geodesic_measure<DistanceType, ResultType> {
  using base_type = geodesic_measure<DistanceType, ResultType>;

  template <typename Graph>
  ResultType operator()(DistanceType d, const Graph& /*unused*/) {
    return (d == base_type::infinite_distance()) ? base_type::zero_result() : rec(result_type(d));
  }
  Reciprocal rec;
};

template <typename Graph, typename DistanceMap>
auto measure_closeness(const Graph& /*unused*/, DistanceMap /*unused*/) {
  return closeness_measure<property_traits_value_t<DistanceMap>, double>();
}

template <typename T, typename Graph, typename DistanceMap>
auto measure_closeness(const Graph& /*unused*/, DistanceMap /*unused*/) {
  return closeness_measure<property_traits_value_t<DistanceMap>, T>();
}

template <typename T, typename Graph, typename DistanceMap, typename Reciprocal>
auto measure_closeness(const Graph& /*unused*/, DistanceMap /*unused*/) {
  return closeness_measure<property_traits_value_t<DistanceMap>, T, Reciprocal>();
}

template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> DistanceMap,
          concepts::DistanceMeasure<G> Measure, typename Combinator>
auto closeness_centrality(const G& g, DistanceMap dist, Measure measure, Combinator combine) {
  using Distance = property_traits_value_t<DistanceMap>;
  Distance n = combine_distances(g, dist, combine, Distance(0));
  return measure(n, g);
}

template <concepts::Graph G, concepts::ReadableVertexPropertyMap<G> DistanceMap, concepts::DistanceMeasure<G> Measure>
auto closeness_centrality(const G& g, DistanceMap dist, Measure measure) {
  return closeness_centrality(g, dist, measure, std::plus<>());
}

template <concepts::Graph G, concepts::ReadableVertexPropertyMap<G> DistanceMap>
double closeness_centrality(const G& g, DistanceMap dist) {
  return closeness_centrality(g, dist, measure_closeness(g, dist));
}

template <typename T, concepts::Graph G, concepts::ReadableVertexPropertyMap<G> DistanceMap>
T closeness_centrality(const G& g, DistanceMap dist) {
  return closeness_centrality(g, dist, measure_closeness<T>(g, dist));
}

template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> DistanceMatrixMap,
          concepts::WritableVertexPropertyMap<G> CentralityMap, concepts::DistanceMeasure<G> Measure>
requires concepts::ReadableVertexPropertyMap<property_traits_value_t<DistanceMatrixMap>, G>
void all_closeness_centralities(const G& g, DistanceMatrixMap dist, CentralityMap cent, Measure measure) {
  for (auto v : vertices(g)) {
    put(cent, v, closeness_centrality(g, get(dist, v), measure));
  }
}

template <concepts::Graph G, concepts::ReadableVertexPropertyMap<G> DistanceMatrixMap,
          concepts::WritableVertexPropertyMap<G> CentralityMap>
requires concepts::ReadableVertexPropertyMap<property_traits_value_t<DistanceMatrixMap>, G>
void all_closeness_centralities(const G& g, DistanceMatrixMap dist, CentralityMap cent) {
  using DistanceMap = property_traits_value_t<DistanceMatrixMap>;
  using Result = property_traits_value_t<CentralityMap>;
  all_closeness_centralities(g, dist, cent, measure_closeness<Result>(g, DistanceMap()));
}

}  // namespace bagl

#endif  // BAGL_BAGL_CLOSENESS_CENTRALITY_H_
