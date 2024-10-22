// Copyright Daniel Trebbien 2010.
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_STOER_WAGNER_MIN_CUT_H_
#define BAGL_BAGL_STOER_WAGNER_MIN_CUT_H_

#include <cassert>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/buffer_concepts.h"
#include "bagl/d_ary_heap.h"
#include "bagl/exception.h"
#include "bagl/graph_traits.h"
#include "bagl/maximum_adjacency_search.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"

namespace bagl {

namespace min_cut_detail {
/**
 * \brief Performs a phase of the Stoer-Wagner min-cut algorithm
 *
 * Performs a phase of the Stoer-Wagner min-cut algorithm.
 *
 * As described by Stoer & Wagner (1997), a phase is simply a maximum
 * adjacency search (also called a maximum cardinality search), which
 * results in the selection of two vertices \em s and \em t, and, as a side
 * product, a minimum <em>s</em>-<em>t</em> cut of the input graph. Here,
 * the input graph is basically \p g, but some vertices are virtually
 * assigned to others as a way of viewing \p g as a graph with some sets of
 * vertices merged together.
 *
 * This implementation is a translation of pseudocode by Professor Uri
 * Zwick, School of Computer Science, Tel Aviv University.
 *
 * \pre \p g is a connected, undirected graph
 * \param[in] g the input graph
 * \param[in] assignments a read/write property map from each vertex to the
 *                        vertex that it is assigned to
 * \param[in] assigned_vertices a list of vertices that are assigned to
 *                             others
 * \param[in] weights a readable property map from each edge to its
 *                    weight (a non-negative value)
 * \param[out] pq a keyed, updatable max-priority queue
 * \returns a tuple (\em s, \em t, \em w) of the "<em>s</em>" and
 *          "<em>t</em>" of the minimum <em>s</em>-<em>t</em> cut and the
 *          cut weight \em w of the minimum <em>s</em>-<em>t</em> cut.
 * \see http://www.cs.tau.ac.il/~zwick/grad-algo-08/gmc.pdf
 *
 * \author Daniel Trebbien
 * \date 2010-09-11
 */
template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> VertexAssignmentMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, concepts::KeyedUpdatableQueue KeyedUpdatablePriorityQueue>
requires concepts::IncidenceGraph<G> &&
    std::convertible_to<graph_vertex_descriptor_t<G>, property_traits_value_t<VertexAssignmentMap>>
        std::tuple<graph_vertex_descriptor_t<G>, graph_vertex_descriptor_t<G>, property_traits_value_t<WeightMap>>
        stoer_wagner_phase(const G& g, VertexAssignmentMap assignments,
                           const std::unordered_set<graph_vertex_descriptor_t<G>>& assigned_vertices, WeightMap weights,
                           KeyedUpdatablePriorityQueue& pq) {
  using weight_type = property_traits_value_t<WeightMap>;

  assert(pq.empty());
  auto keys = pq.keys();

  for (auto v : vertices(g)) {
    if (v == get(assignments, v)) {  // foreach u \in V do
      put(keys, v, weight_type{});

      pq.push(v);
    }
  }

  assert(pq.size() >= 2);

  auto s = graph_traits<G>::null_vertex();
  auto t = graph_traits<G>::null_vertex();
  weight_type w;
  while (!pq.empty()) {       // while PQ \neq {} do
    const auto u = pq.top();  // u = extractmax(PQ)
    w = get(keys, u);
    pq.pop();

    s = t;
    t = u;

    for (auto e : out_edges(u, g)) {
      const auto v = get(assignments, target(e, g));
      if (pq.contains(v)) {
        // if v \in PQ then
        // increasekey(PQ, v, wA(v) + w(u, v))
        put(keys, v, get(keys, v) + get(weights, e));
        pq.update(v);
      }
    }

    for (auto u_prime : assigned_vertices) {
      if (get(assignments, u_prime) == u) {
        for (auto e : out_edges(u_prime, g)) {
          const auto v = get(assignments, target(e, g));
          if (pq.contains(v)) {
            // if v \in PQ then
            // increasekey(PQ, v, wA(v) + w(u, v))
            put(keys, v, get(keys, v) + get(weights, e));
            pq.update(v);
          }
        }
      }
    }
  }

  return {s, t, w};
}

/**
 * \brief Computes a min-cut of the input graph
 *
 * Computes a min-cut of the input graph using the Stoer-Wagner algorithm.
 *
 * \pre \p g is a connected, undirected graph
 * \pre <code>pq.empty()</code>
 * \param[in] g the input graph
 * \param[in] weights a readable property map from each edge to its weight
 * (a non-negative value) \param[out] parities a writable property map from
 * each vertex to a bool type object for distinguishing the two vertex sets
 * of the min-cut \param[out] assignments a read/write property map from
 * each vertex to a \c vertex_descriptor object. This map serves as work
 * space, and no particular meaning should be derived from property values
 *     after completion of the algorithm.
 * \param[out] pq a keyed, updatable max-priority queue
 * \returns the cut weight of the min-cut
 * \see
 * http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.114.6687&rep=rep1&type=pdf
 * \see
 * http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.31.614&rep=rep1&type=pdf
 *
 * \author Daniel Trebbien
 * \date 2010-09-11
 */
template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::WritableVertexPropertyMap<G> ParityMap, concepts::ReadWriteVertexPropertyMap<G> VertexAssignmentMap,
          concepts::KeyedUpdatableQueue KeyedUpdatablePriorityQueue>
requires concepts::IncidenceGraph<G> &&
    std::convertible_to<graph_vertex_descriptor_t<G>, property_traits_value_t<VertexAssignmentMap>>
        property_traits_value_t<WeightMap> stoer_wagner_min_cut(const G& g, WeightMap weights, ParityMap parities,
                                                                VertexAssignmentMap assignments,
                                                                KeyedUpdatablePriorityQueue& pq) {
  using vertex_descriptor = graph_vertex_descriptor_t<G>;
  using weight_type = property_traits_value_t<WeightMap>;
  using parity_type = property_traits_value_t<ParityMap>;

  std::size_t n = num_vertices(g);

  std::unordered_set<vertex_descriptor> assigned_vertices;

  // initialize `assignments` (all vertices are initially assigned to
  // themselves)
  for (auto v : vertices(g)) {
    put(assignments, v, v);
  }

  auto [s, t, bestW] = stoer_wagner_phase(g, assignments, assigned_vertices, weights, pq);
  assert(s != t);
  for (auto v : vertices(g)) {
    put(parities, v, parity_type(v == t ? 1 : 0));
  }
  put(assignments, t, s);
  assigned_vertices.insert(t);
  --n;

  for (; n >= 2; --n) {
    weight_type w;
    std::tie(s, t, w) = stoer_wagner_phase(g, assignments, assigned_vertices, weights, pq);
    assert(s != t);

    if (w < bestW) {
      for (auto v : vertices(g)) {
        put(parities, v, parity_type(get(assignments, v) == t ? 1 : 0));

        if (get(assignments, v) == t) {
          // all vertices that were assigned to t are now
          // assigned to s
          put(assignments, v, s);
        }
      }

      bestW = w;
    } else {
      for (auto v : vertices(g)) {
        if (get(assignments, v) == t) {
          // all vertices that were assigned to t are now
          // assigned to s
          put(assignments, v, s);
        }
      }
    }
    put(assignments, t, s);
    assigned_vertices.insert(t);
  }

  assert(pq.empty());

  return bestW;
}
}  // namespace min_cut_detail

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::WritableVertexPropertyMap<G> ParityMap, concepts::ReadWriteVertexPropertyMap<G> VertexAssignmentMap,
          concepts::KeyedUpdatableQueue KeyedUpdatablePriorityQueue>
requires concepts::IncidenceGraph<G> &&
    std::convertible_to<graph_vertex_descriptor_t<G>, property_traits_value_t<VertexAssignmentMap>>
        property_traits_value_t<WeightMap> stoer_wagner_min_cut(
            const G& g, WeightMap weights,
            ParityMap parities,               // dummy_property_map()
            VertexAssignmentMap assignments,  // vector_property_map(num_vertices(g), get(vertex_index, g),
                                              // graph_traits<G>::null_vertex()).ref()
            KeyedUpdatablePriorityQueue&
                pq) {  // make_d_ary_heap_indirect<graph_vertex_descriptor_t<G>, 4>(..., std::greater<>())
  static_assert(is_undirected_graph_v<G>);

  std::size_t n = num_vertices(g);
  if (n < 2) {
    throw bad_graph("the input graph must have at least two vertices.");
  }
  if (!pq.empty()) {
    throw std::invalid_argument("the max-priority queue must be empty initially.");
  }

  return min_cut_detail::stoer_wagner_min_cut(g, weights, parities, assignments, pq);
}

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::WritableVertexPropertyMap<G> ParityMap, concepts::KeyedUpdatableQueue KeyedUpdatablePriorityQueue>
requires concepts::IncidenceGraph<G> property_traits_value_t<WeightMap> stoer_wagner_min_cut(
    const G& g, WeightMap weights, ParityMap parities, KeyedUpdatablePriorityQueue& pq) {
  return stoer_wagner_min_cut(
      g, weights, parities,
      vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<G>::null_vertex()).ref(), pq);
}

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::KeyedUpdatableQueue KeyedUpdatablePriorityQueue>
requires concepts::IncidenceGraph<G> property_traits_value_t<WeightMap> stoer_wagner_min_cut(
    const G& g, WeightMap weights, KeyedUpdatablePriorityQueue& pq) {
  return stoer_wagner_min_cut(g, weights, dummy_property_map(), pq);
}

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::WritableVertexPropertyMap<G> ParityMap, concepts::KeyedUpdatableQueue KeyedUpdatablePriorityQueue,
          concepts::ReadableVertexIndexMap<G> VertexIndexMap>
requires concepts::IncidenceGraph<G> property_traits_value_t<WeightMap> stoer_wagner_min_cut(
    const G& g, WeightMap weights, ParityMap parities, KeyedUpdatablePriorityQueue& pq, VertexIndexMap v_index) {
  return stoer_wagner_min_cut(g, weights, parities,
                              vector_property_map(num_vertices(g), v_index, graph_traits<G>::null_vertex()).ref(), pq);
}

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::KeyedUpdatableQueue KeyedUpdatablePriorityQueue, concepts::ReadableVertexIndexMap<G> VertexIndexMap>
requires concepts::IncidenceGraph<G> property_traits_value_t<WeightMap> stoer_wagner_min_cut(
    const G& g, WeightMap weights, KeyedUpdatablePriorityQueue& pq, VertexIndexMap v_index) {
  return stoer_wagner_min_cut(g, weights, dummy_property_map(),
                              vector_property_map(num_vertices(g), v_index, graph_traits<G>::null_vertex()).ref(), pq);
}

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::WritableVertexPropertyMap<G> ParityMap, concepts::ReadableVertexIndexMap<G> VertexIndexMap>
requires concepts::IncidenceGraph<G> property_traits_value_t<WeightMap> stoer_wagner_min_cut(const G& g,
                                                                                             WeightMap weights,
                                                                                             ParityMap parities,
                                                                                             VertexIndexMap v_index) {
  using weight_type = property_traits_value_t<WeightMap>;
  auto prio = vector_property_map(num_vertices(g), v_index, weight_type{});
  auto index_in_heap = vector_property_map(num_vertices(g), v_index, std::size_t{0});
  auto pq =
      make_d_ary_heap_indirect<graph_vertex_descriptor_t<G>, 4>(prio.ref(), index_in_heap.ref(), std::greater<>());
  return stoer_wagner_min_cut(g, weights, parities, pq, v_index);
}

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadableVertexIndexMap<G> VertexIndexMap>
requires concepts::IncidenceGraph<G> property_traits_value_t<WeightMap> stoer_wagner_min_cut(const G& g,
                                                                                             WeightMap weights,
                                                                                             VertexIndexMap v_index) {
  return stoer_wagner_min_cut(g, weights, dummy_property_map(), v_index);
}

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap>
requires concepts::IncidenceGraph<G> property_traits_value_t<WeightMap> stoer_wagner_min_cut(const G& g,
                                                                                             WeightMap weights) {
  return stoer_wagner_min_cut(g, weights, get(vertex_index, g));
}

}  // namespace bagl

#endif  // !BAGL_BAGL_STOER_WAGNER_MIN_CUT_H_
