// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_ECCENTRICITY_H_
#define BAGL_BAGL_ECCENTRICITY_H_

#include <iterator>

#include "bagl/geodesic.h"
#include "bagl/graph_concepts.h"
#include "bagl/property_map.h"

namespace bagl {

template <concepts::Graph G, concepts::ReadableVertexPropertyMap<G> DistanceMap,
          concepts::PropertyCombinator<DistanceMap> Combine>
auto eccentricity(const G& g, DistanceMap dist, Combine combine) {
  return combine_distances(g, dist, combine);
}

template <concepts::Graph G, concepts::ReadableVertexPropertyMap<G> DistanceMap>
auto eccentricity(const G& g, DistanceMap dist) {
  return eccentricity(g, dist, maximize());
}

template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> DistanceMatrix,
          concepts::WritableVertexPropertyMap<G> EccentricityMap>
auto all_eccentricities(const G& g, const DistanceMatrix& dist, EccentricityMap ecc) {
  using Eccentricity = property_traits_value_t<EccentricityMap>;

  Eccentricity r = numeric_values<Eccentricity>::infinity();
  Eccentricity d = numeric_values<Eccentricity>::zero();
  for (auto v : vertices(g)) {
    auto dm = get(dist, v);
    auto e = eccentricity(g, dm);
    put(ecc, v, e);

    // track the radius and diameter at the same time
    r = std::min(r, e);
    d = std::max(d, e);
  }
  return std::pair{r, d};
}

template <concepts::VertexListGraph G, concepts::WritableVertexPropertyMap<G> EccentricityMap>
auto radius_and_diameter(const G& g, EccentricityMap ecc) {
  using Eccentricity = property_traits_value_t<EccentricityMap>;

  Eccentricity radius = numeric_values<Eccentricity>::infinity();
  Eccentricity diameter = numeric_values<Eccentricity>::zero();
  for (auto v : vertices(g)) {
    Eccentricity cur = get(ecc, v);
    radius = std::min(radius, cur);
    diameter = std::max(diameter, cur);
  }
  return std::pair{radius, diameter};
}

template <concepts::VertexListGraph G, concepts::WritableVertexPropertyMap<G> EccentricityMap>
auto radius(const G& g, EccentricityMap ecc) {
  return radius_and_diameter(g, ecc).first;
}

template <concepts::VertexListGraph G, concepts::WritableVertexPropertyMap<G> EccentricityMap>
auto diameter(const G& g, EccentricityMap ecc) {
  return radius_and_diameter(g, ecc).second;
}

}  // namespace bagl

#endif  // BAGL_BAGL_ECCENTRICITY_H_
