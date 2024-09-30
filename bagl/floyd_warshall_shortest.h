// Copyright 2002 Rensselaer Polytechnic Institute
//  Authors: Lauren Foutz, Scott Hill
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_FLOYD_WARSHALL_SHORTEST_H_
#define BAGL_BAGL_FLOYD_WARSHALL_SHORTEST_H_

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/numeric_values.h"
#include "bagl/property.h"
#include "bagl/property_map.h"
#include "bagl/relax.h"

namespace bagl {

namespace floyd_warshall_detail {
template <typename T, typename BinaryPredicate>
T min_with_compare(const T& x, const T& y, const BinaryPredicate& compare) {
  if (compare(x, y)) {
    return x;
  } else {
    return y;
  }
}

template <typename M>
using property_matrix_value_t = property_traits_value_t<property_traits_value_t<M>>;

template <concepts::VertexListGraph G, typename DistanceMatrix,
          concepts::PropertyComparator<property_traits_value_t<DistanceMatrix>> Compare,
          concepts::PropertyCombinator<property_traits_value_t<DistanceMatrix>> Combine, typename Infinity,
          typename Zero>
bool floyd_warshall_dispatch(const G& g, DistanceMatrix& d, const Compare& compare, const Combine& combine,
                             property_matrix_value_t<DistanceMatrix> inf,
                             property_matrix_value_t<DistanceMatrix> zero) {
  for (auto k : vertices(g)) {
    for (auto i : vertices(g)) {
      if (d[i][k] != inf) {
        for (auto j : vertices(g)) {
          if (d[k][j] != inf) {
            d[i][j] = min_with_compare(d[i][j], combine(d[i][k], d[k][j]), compare);
          }
        }
      }
    }
  }

  for (auto v : vertices(g)) {
    if (compare(d[v][v], zero)) {
      return false;
    }
  }
  return true;
}
}  // namespace floyd_warshall_detail

template <concepts::VertexListGraph G, typename DistanceMatrix,
          concepts::PropertyComparator<property_traits_value_t<DistanceMatrix>> Compare,
          concepts::PropertyCombinator<property_traits_value_t<DistanceMatrix>> Combine>
bool floyd_warshall_initialized_all_pairs_shortest_paths(
    const G& g, DistanceMatrix& d, Compare compare, Combine combine,
    floyd_warshall_detail::property_matrix_value_t<DistanceMatrix> inf,
    floyd_warshall_detail::property_matrix_value_t<DistanceMatrix> zero) {
  return floyd_warshall_detail::floyd_warshall_dispatch(g, d, compare, combine, inf, zero);
}

template <concepts::VertexListGraph G, typename DistanceMatrix>
bool floyd_warshall_initialized_all_pairs_shortest_paths(const G& g, DistanceMatrix& d) {
  using D = floyd_warshall_detail::property_matrix_value_t<DistanceMatrix>;
  return floyd_warshall_detail::floyd_warshall_dispatch(g, d, std::less<>(), closed_plus<>(),
                                                        numeric_values<D>::infinity(), numeric_values<D>::zero());
}

template <concepts::VertexAndEdgeListGraph G, typename DistanceMatrix, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::PropertyComparator<property_traits_value_t<DistanceMatrix>> Compare,
          concepts::PropertyCombinator<property_traits_value_t<DistanceMatrix>> Combine>
bool floyd_warshall_all_pairs_shortest_paths(const G& g, DistanceMatrix& d, const WeightMap& w, Compare compare,
                                             Combine combine,
                                             floyd_warshall_detail::property_matrix_value_t<DistanceMatrix> inf,
                                             floyd_warshall_detail::property_matrix_value_t<DistanceMatrix> zero) {
  for (auto u : vertices(g)) {
    for (auto v : vertices(g)) {
      d[u][v] = inf;
    }
  }

  for (auto u : vertices(g)) {
    d[u][u] = zero;
  }

  for (auto e : edges(g)) {
    if (d[source(e, g)][target(e, g)] != inf) {
      d[source(e, g)][target(e, g)] =
          floyd_warshall_detail::min_with_compare(get(w, e), d[source(e, g)][target(e, g)], compare);
    } else {
      d[source(e, g)][target(e, g)] = get(w, e);
    }
  }

  if constexpr (is_undirected_graph_v<G>) {
    for (auto e : edges(g)) {
      if (d[target(e, g)][source(e, g)] != inf) {
        d[target(e, g)][source(e, g)] =
            floyd_warshall_detail::min_with_compare(get(w, e), d[target(e, g)][source(e, g)], compare);
      } else {
        d[target(e, g)][source(e, g)] = get(w, e);
      }
    }
  }

  return floyd_warshall_detail::floyd_warshall_dispatch(g, d, compare, combine, inf, zero);
}

template <concepts::VertexAndEdgeListGraph G, typename DistanceMatrix, concepts::ReadableEdgePropertyMap<G> WeightMap>
bool floyd_warshall_all_pairs_shortest_paths(const G& g, DistanceMatrix& d, const WeightMap& w) {
  using D = floyd_warshall_detail::property_matrix_value_t<DistanceMatrix>;
  return floyd_warshall_all_pairs_shortest_paths(g, d, w, std::less<>(), closed_plus<>(), numeric_values<D>::infinity(),
                                                 numeric_values<D>::zero());
}

}  // namespace bagl

#endif  // BAGL_BAGL_FLOYD_WARSHALL_SHORTEST_H_
