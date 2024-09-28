// (C) Copyright 2007 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GEODESIC_H_
#define BAGL_BAGL_GEODESIC_H_

#include <functional>

#include "bagl/graph_concepts.h"
#include "bagl/numeric_values.h"
#include "bagl/property_map.h"

// TODO: Should this really be in detail?

namespace bagl {
// This is a very good discussion on centrality measures. While I can't
// say that this has been the motivating factor for the design and
// implementation of ths centrality framework, it does provide a single
// point of reference for defining things like degree and closeness
// centrality. Plus, the bibliography seems fairly complete.
//
//     @article{citeulike:1144245,
//         author = {Borgatti, Stephen  P. and Everett, Martin  G.},
//         citeulike-article-id = {1144245},
//         doi = {10.1016/j.socnet.2005.11.005},
//         journal = {Social Networks},
//         month = {October},
//         number = {4},
//         pages = {466--484},
//         priority = {0},
//         title = {A Graph-theoretic perspective on centrality},
//         url = {https://doi.org/10.1016/j.socnet.2005.11.005},
//             volume = {28},
//             year = {2006}
//         }
//     }

// Note that this assumes T == property_traits<DistanceMap>::value_type
// and that the args and return of combine are also T.
template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> DistanceMap,
          std::relation<property_traits_value_t<DistanceMap>, property_traits_value_t<DistanceMap>> Combinator>
          requires concepts::NumericValue<property_traits_value_t<DistanceMap>>
auto combine_distances(const G& g, DistanceMap dist, Combinator combine, property_traits_value_t<DistanceMap> init) {
  using Distance = property_traits_value_t<DistanceMap>;
  using DistanceNumbers = numeric_values<Distance>;

  // If there's ever an infinite distance, then we simply return
  // infinity. Note that this /will/ include the a non-zero
  // distance-to-self in the combined values. However, this is usually
  // zero, so it shouldn't be too problematic.
  Distance ret = init;
  for (auto v : vertices(g)) {
    if (get(dist, v) != DistanceNumbers::infinity()) {
      ret = combine(ret, get(dist, v));
    } else {
      ret = DistanceNumbers::infinity();
      break;
    }
  }
  return ret;
}

// Similar to std::plus<T>, but maximizes parameters
// rather than adding them.
struct maximize {
  template <typename T>
  T operator()(const T& x, const T& y) const {
    return std::max(x, y);
  }
};

// Another helper, like maximize() to help abstract functional
// concepts. This is trivially instantiated for builtin numeric
// types, but should be specialized for those types that have
// discrete notions of reciprocals.
struct reciprocal {
  template <typename T>
  T operator()(const T& t) {
    return T{1} / t;
  }
};

// This type defines the basic facilities used for computing values
// based on the geodesic distances between vertices. Examples include
// closeness centrality and mean geodesic distance.
template <concepts::NumericValue DistanceType, concepts::NumericValue ResultType>
struct geodesic_measure {
  static auto infinite_distance() { return numeric_values<DistanceType>::infinity(); }
  static auto infinite_result() { return numeric_values<ResultType>::infinity(); }
  static auto zero_result() { return numeric_values<ResultType>::zero(); }
};

}  // namespace bagl

#endif  // BAGL_BAGL_GEODESIC_H_
