// Copyright 2007 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GEODESIC_DISTANCE_H_
#define BAGL_BAGL_GEODESIC_DISTANCE_H_

#include "bagl/geodesic.h"
#include "bagl/graph_concepts.h"
#include "bagl/property_map.h"

namespace bagl {

template <concepts::VertexListGraph G, concepts::NumericValue DistanceType, concepts::NumericValue ResultType,
          typename Divides = std::divides<ResultType> >
struct mean_geodesic_measure : public geodesic_measure<DistanceType, ResultType> {
  using base_type = geodesic_measure<DistanceType, ResultType>;
  using distance_type = typename base_type::distance_type;
  using result_type = typename base_type::result_type;

  result_type operator()(distance_type d, const G& g) {
    return (d == base_type::infinite_distance()) ? base_type::infinite_result()
                                                 : div(result_type(d), result_type(num_vertices(g) - 1));
  }
  Divides div;
};

template <concepts::VertexListGraph G, typename DistanceMap>
auto measure_mean_geodesic(const G&, DistanceMap) {
  return mean_geodesic_measure<G, property_traits_value_t<DistanceMap>, double>();
}

template <typename T, typename Graph, typename DistanceMap>
auto measure_mean_geodesic(const Graph&, DistanceMap) {
  return mean_geodesic_measure<Graph, property_traits_value_t<DistanceMap>, T>();
}

// This is a little different because it's expected that the result type
// should (must?) be the same as the distance type. There's a type of
// transitivity in this thinking... If the average of distances has type
// X then the average of x's should also be type X. Is there a case where this
// is not true?
//
// This type is a little under-genericized... It needs generic parameters
// for addition and division.
template <concepts::VertexListGraph G, concepts::NumericValue DistanceType>
struct mean_graph_distance_measure : public geodesic_measure<DistanceType, DistanceType> {
  using base_type = geodesic_measure<DistanceType, DistanceType>;
  using distance_type = typename base_type::distance_type;
  using result_type = typename base_type::result_type;

  result_type operator()(distance_type d, const G& g) {
    if (d == base_type::infinite_distance()) {
      return base_type::infinite_result();
    }

    return d / result_type(num_vertices(g));
  }
};

template <typename Graph, typename DistanceMap>
auto measure_graph_mean_geodesic(const Graph&, DistanceMap) {
  return mean_graph_distance_measure<Graph, property_traits_value_t<DistanceMap> >();
}

template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> DistanceMap,
          concepts::DistanceMeasure<G, property_traits_value_t<DistanceMap>> Measure,
          concepts::PropertyCombinator<DistanceMap> Combine>
auto mean_geodesic(const G& g, DistanceMap dist, Measure measure, Combine combine) {
  using Distance = property_traits_value_t<DistanceMap>;
  Distance n = combine_distances(g, dist, combine, Distance(0));
  return measure(n, g);
}

template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> DistanceMap,
          concepts::DistanceMeasure<G, property_traits_value_t<DistanceMap>> Measure>
auto mean_geodesic(const G& g, DistanceMap dist, Measure measure) {
  using Distance = typename Measure::distance_type;
  return mean_geodesic(g, dist, measure, std::plus<Distance>());
}

template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> DistanceMap>
double mean_geodesic(const G& g, DistanceMap dist) {
  return mean_geodesic(g, dist, measure_mean_geodesic(g, dist));
}

template <typename T, concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> DistanceMap>
T mean_geodesic(const G& g, DistanceMap dist) {
  return mean_geodesic(g, dist, measure_mean_geodesic<T>(g, dist));
}

template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> DistanceMatrixMap,
          concepts::WritableVertexPropertyMap<G> GeodesicMap,
          concepts::DistanceMeasure<G, property_traits_value_t<property_traits_value_t<DistanceMatrixMap>>> Measure>
property_traits_value_t<GeodesicMap> all_mean_geodesics(const G& g, DistanceMatrixMap dist, GeodesicMap geo,
                                                        Measure measure) {
  using Result = typename Measure::result_type;
  // NOTE: We could compute the mean geodesic here by performing additional
  // computations (i.e., adding and dividing). However, I don't really feel
  // like fully genericizing the entire operation yet so I'm not going to.

  Result inf = numeric_values<Result>::infinity();
  Result sum = numeric_values<Result>::zero();
  for (auto v : vertices(g)) {
    auto dm = get(dist, v);
    Result r = mean_geodesic(g, dm, measure);
    put(geo, v, r);

    // compute the sum along with geodesics
    if (r == inf) {
      sum = inf;
    } else if (sum != inf) {
      sum += r;
    }
  }

  // return the average of averages.
  return sum / Result(num_vertices(g));
}

template <concepts::Graph G, concepts::ReadableVertexPropertyMap<G> DistanceMatrixMap,
          concepts::WritableVertexPropertyMap<G> GeodesicMap>
property_traits_value_t<GeodesicMap> all_mean_geodesics(const G& g, DistanceMatrixMap dist, GeodesicMap geo) {
  using DistanceMap = property_traits_value_t<DistanceMatrixMap>;
  using Result = property_traits_value_t<GeodesicMap>;

  return all_mean_geodesics(g, dist, geo, measure_mean_geodesic<Result>(g, DistanceMap()));
}

template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> GeodesicMap,
          concepts::DistanceMeasure<G, property_traits_value_t<GeodesicMap>> Measure>
auto small_world_distance(const G& g, GeodesicMap geo, Measure measure) {
  using Distance = property_traits_value_t<GeodesicMap>;
  auto sum = combine_distances(g, geo, std::plus<>(), Distance(0));
  return measure(sum, g);
}

template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> GeodesicMap>
property_traits_value_t<GeodesicMap> small_world_distance(const G& g, GeodesicMap geo) {
  return small_world_distance(g, geo, measure_graph_mean_geodesic(g, geo));
}

}  // namespace bagl

#endif  // BAGL_BAGL_GEODESIC_DISTANCE_H_
