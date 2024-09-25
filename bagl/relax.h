// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek,
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_RELAX_H_
#define BAGL_BAGL_RELAX_H_

#include <functional>
#include <limits>  // for numeric limits

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"

namespace bagl {

// The following version of the plus functor prevents
// problems due to overflow at positive infinity.

template <class T>
struct closed_plus {
  const T inf;

  constexpr closed_plus() : inf((std::numeric_limits<T>::max)()) {}
  constexpr explicit closed_plus(T a_inf) : inf(a_inf) {}

  constexpr T operator()(const T& a, const T& b) const {
    if (a == inf || b == inf) {
      return inf;
    }
    return a + b;
  }
};

template <concepts::IncidenceGraph G, concepts::ReadablePropertyMap<graph_edge_descriptor_t<G>> WeightMap,
          concepts::ReadWritePropertyMap<graph_vertex_descriptor_t<G>> PredecessorMap,
          concepts::ReadWritePropertyMap<graph_vertex_descriptor_t<G>> DistanceMap, class BinaryFunction,
          class BinaryPredicate>
bool relax(graph_edge_descriptor_t<G> e, const G& g, const WeightMap& w, PredecessorMap& p, DistanceMap& d,
           const BinaryFunction& combine, const BinaryPredicate& compare) {
  auto u = source(e, g);
  auto v = target(e, g);
  const auto d_u = get(d, u);
  const auto d_v = get(d, v);
  const auto& w_e = get(w, e);

  // The seemingly redundant comparisons after the distance puts are to
  // ensure that extra floating-point precision in x87 registers does not
  // lead to relax() returning true when the distance did not actually
  // change.
  if (compare(combine(d_u, w_e), d_v)) {
    put(d, v, combine(d_u, w_e));
    if (compare(get(d, v), d_v)) {
      put(p, v, u);
      return true;
    }

    return false;
  }
  if (is_undirected_graph_v<G> && compare(combine(d_v, w_e), d_u)) {
    put(d, u, combine(d_v, w_e));
    if (compare(get(d, u), d_u)) {
      put(p, u, v);
      return true;
    }

    return false;
  }
  return false;
}

template <concepts::IncidenceGraph G, concepts::ReadablePropertyMap<graph_edge_descriptor_t<G>> WeightMap,
          concepts::ReadWritePropertyMap<graph_vertex_descriptor_t<G>> PredecessorMap,
          concepts::ReadWritePropertyMap<graph_vertex_descriptor_t<G>> DistanceMap, class BinaryFunction,
          class BinaryPredicate>
bool relax_target(graph_edge_descriptor_t<G> e, const G& g, const WeightMap& w, PredecessorMap& p, DistanceMap& d,
                  const BinaryFunction& combine, const BinaryPredicate& compare) {
  const auto u = source(e, g);
  const auto v = target(e, g);
  const auto d_u = get(d, u);
  const auto d_v = get(d, v);
  const auto& w_e = get(w, e);

  // The seemingly redundant comparisons after the distance puts are to
  // ensure that extra floating-point precision in x87 registers does not
  // lead to relax() returning true when the distance did not actually
  // change.
  if (compare(combine(d_u, w_e), d_v)) {
    put(d, v, combine(d_u, w_e));
    if (compare(get(d, v), d_v)) {
      put(p, v, u);
      return true;
    }
  }
  return false;
}

template <concepts::IncidenceGraph G, concepts::ReadablePropertyMap<graph_edge_descriptor_t<G>> WeightMap,
          concepts::ReadWritePropertyMap<graph_vertex_descriptor_t<G>> PredecessorMap,
          concepts::ReadWritePropertyMap<graph_vertex_descriptor_t<G>> DistanceMap>
bool relax(graph_edge_descriptor_t<G> e, const G& g, WeightMap w, PredecessorMap p, DistanceMap d) {
  using D = property_traits_value_t<DistanceMap>;
  using Combine = closed_plus<D>;
  using Compare = std::less<D>;
  return relax(e, g, w, p, d, Combine(), Compare());
}

}  // namespace bagl

#endif  // BAGL_BAGL_RELAX_H_
