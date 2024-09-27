// Copyright 2012 Fernando Vilas
//           2010 Daniel Trebbien
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_MAXIMUM_ADJACENCY_SEARCH_H_
#define BAGL_BAGL_MAXIMUM_ADJACENCY_SEARCH_H_

// The maximum adjacency search algorithm was originally part of the
// Stoer-Wagner min cut implementation by Daniel Trebbien. It has been
// broken out into its own file to be a public search algorithm, with
// visitor concepts.

// This is an implementation of the maximum adjacency search on an
// undirected graph. It allows a visitor object to perform some
// operation on each vertex as that vertex is visited.
//
// The algorithm runs as follows:
//
// Initialize all nodes to be unvisited (reach count = 0)
//   and call vis.initialize_vertex
// For i = number of nodes in graph downto 1
//   Select the unvisited node with the highest reach count
//     The user provides the starting node to break the first tie,
//     but future ties are broken arbitrarily
//   Visit the node by calling vis.start_vertex
//   Increment the reach count for all unvisited neighbors
//     and call vis.examine_edge for each of these edges
//   Mark the node as visited and call vis.finish_vertex

#include <concepts>
#include <set>
#include <tuple>

#include "bagl/buffer_concepts.h"
#include "bagl/exception.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/vebl_d_ary_tree.h"
#include "bagl/visitors.h"

namespace bagl {

namespace concepts {
template <typename V, typename G>
concept MASVisitor = std::copy_constructible<V> &&
    requires(const V& vis, const G& g, graph_vertex_descriptor_t<G> u, graph_edge_descriptor_t<G> e) {
  vis.initialize_vertex(u, g);
  vis.start_vertex(u, g);
  vis.examine_edge(e, g);
  vis.finish_vertex(u, g);
};

}  // namespace concepts

template <typename Visitors = null_visitor>
class mas_visitor {
 public:
  mas_visitor() = default;
  explicit mas_visitor(Visitors vis) : vis_(vis) {}

  template <typename Vertex, typename Graph>
  void initialize_vertex(Vertex u, Graph& g) {
    visitors_detail::invoke_initialize_vertex_on_all(vis_, u, g);
  }

  template <typename Vertex, typename Graph>
  void start_vertex(Vertex u, Graph& g) {
    visitors_detail::invoke_start_vertex_on_all(vis_, u, g);
  }

  template <typename Edge, typename Graph>
  void examine_edge(Edge e, Graph& g) {
    visitors_detail::invoke_examine_edge_on_all(vis_, e, g);
  }

  template <typename Vertex, typename Graph>
  void finish_vertex(Vertex u, Graph& g) {
    visitors_detail::invoke_finish_vertex_on_all(vis_, u, g);
  }

 protected:
  Visitors vis_;
};
template <typename Visitors>
auto make_mas_visitor(Visitors vis) {
  return mas_visitor<Visitors>(vis);
}
using default_mas_visitor = mas_visitor<>;

namespace mas_detail {
template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::MASVisitor<G> V, concepts::ReadWriteVertexPropertyMap<G> VertexAssignmentMap,
          concepts::KeyedUpdatableQueue PriorityQueue>
requires concepts::IncidenceGraph<G>
void maximum_adjacency_search(const G& g, WeightMap weights, V vis, graph_vertex_descriptor_t<G> start,
                              VertexAssignmentMap assignments, PriorityQueue pq) {
  using Vertex = graph_vertex_descriptor_t<G>;
  using Weight = property_traits_value_t<WeightMap>;

  std::set<Vertex> assigned_verts;

  // initialize `assignments` (all vertices are initially
  // assigned to themselves)
  for (auto v : vertices(g)) {
    put(assignments, v, v);
  }

  auto keys = pq.keys();

  // set number of visited neighbors for all vertices to 0
  // foreach u \in V do
  for (auto v : vertices(g)) {
    if (v == get(assignments, v)) {
      put(keys, v, Weight(0));
      vis.initialize_vertex(v, g);
      pq.push(v);
    }
  }
  assert(pq.size() >= 2);

  // Give the starting vertex high priority
  put(keys, start, get(keys, start) + num_vertices(g) + 1);
  pq.update(start);

  // start traversing the graph
  // Vertex s, t;
  // Weight w;

  // while PQ \neq {} do
  while (!pq.empty()) {
    // u = extractmax(PQ)
    const Vertex u = pq.top();
    /* Weight w = */ (void)get(keys, u);
    vis.start_vertex(u, g);
    pq.pop();

    // foreach (u, v) \in E do
    for (auto e : out_edges(u, g)) {
      vis.examine_edge(e, g);

      const Vertex v = get(assignments, target(e, g));

      // if v \in PQ then
      if (pq.contains(v)) {
        // increasekey(PQ, v, wA(v) + w(u, v))
        put(keys, v, get(keys, v) + get(weights, e));
        pq.update(v);
      }
    }

    // TODO: This loop will never execute ?!?
    for (auto u_prime : assigned_verts) {
      if (get(assignments, u_prime) != u) {
        continue;
      }
      for (auto e : out_edges(u_prime, g)) {  // foreach (u, v) \in E do
        vis.examine_edge(e, g);

        const Vertex v = get(assignments, target(e, g));

        // if v \in PQ then
        if (pq.contains(v)) {
          // increasekey(PQ, v, wA(v) + w(u, v))
          put(keys, v, get(keys, v) + get(weights, e));
          pq.update(v);
        }
      }
    }
    vis.finish_vertex(u, g);
  }
}
}  // namespace mas_detail

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::MASVisitor<G> V, concepts::ReadWriteVertexPropertyMap<G> VertexAssignmentMap,
          concepts::KeyedUpdatableQueue PriorityQueue>
requires concepts::IncidenceGraph<G>
void maximum_adjacency_search(const G& g, WeightMap weights, V vis, graph_vertex_descriptor_t<G> start,
                              VertexAssignmentMap assignments, PriorityQueue pq) {
  static_assert(is_undirected_graph_v<G>);
  static_assert(std::is_convertible_v<graph_vertex_descriptor_t<G>, property_traits_value_t<VertexAssignmentMap>>);

  if (num_vertices(g) < 2) {
    throw bad_graph("the input graph must have at least two vertices.");
  }
  if (!pq.empty()) {
    throw std::invalid_argument("the max-priority queue must be empty initially.");
  }

  mas_detail::maximum_adjacency_search(g, weights, vis, start, assignments, pq);
}

}  // end namespace bagl

#endif  // BAGL_BAGL_MAXIMUM_ADJACENCY_SEARCH_H_
