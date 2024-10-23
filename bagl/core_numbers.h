// Copyright 2007 Stanford University
// Authors: David Gleich
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_CORE_NUMBERS_H_
#define BAGL_BAGL_CORE_NUMBERS_H_

#include <concepts>

#include "bagl/d_ary_heap.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/single_property_map.h"
#include "bagl/visitors.h"

// History
//
// 30 July 2007
// Added visitors to the implementation
//
// 8 February 2008
// Fixed headers and missing typename

namespace bagl {

// A linear time O(m) algorithm to compute the indegree core number
// of a graph for unweighted graphs.
//
// and a O((n+m) log n) algorithm to compute the in-edge-weight core
// numbers of a weighted graph.
//
// The linear algorithm comes from:
// Vladimir Batagelj and Matjaz Zaversnik, "An O(m) Algorithm for Cores
// Decomposition of Networks."  Sept. 1 2002.

namespace concepts {

template <typename V, typename G>
concept CoreNumbersVisitor = std::copy_constructible<V> &&
    requires(V& vis, const G& g, graph_vertex_descriptor_t<G> u, graph_edge_descriptor_t<G> e) {
  vis.examine_edge(e, g);
  vis.examine_vertex(u, g);
  vis.finish_vertex(u, g);
};

}  // namespace concepts

template <typename Visitors = null_visitors>
class core_numbers_visitor {
 public:
  core_numbers_visitor() = default;
  explicit core_numbers_visitor(Visitors vis) : vis_(vis) {}

  template <class Vertex, class Graph>
  void examine_vertex(Vertex u, Graph& g) {
    visitors_detail::invoke_examine_vertex_on_all(vis_, u, g);
  }

  template <class Edge, class Graph>
  void examine_edge(Edge e, Graph& g) {
    visitors_detail::invoke_examine_edge_on_all(vis_, e, g);
  }

  template <class Vertex, class Graph>
  void finish_vertex(Vertex u, Graph& g) {
    visitors_detail::invoke_finish_vertex_on_all(vis_, u, g);
  }

 protected:
  Visitors vis_;
};

template <typename... Visitors>
auto make_core_numbers_visitor(Visitors&&... vis) {
  if constexpr (sizeof...(Visitors) == 0) {
    return core_numbers_visitor<>();
  } else {
    return core_numbers_visitor(std::tuple<std::decay_t<Visitors>...>(std::forward<Visitors>(vis)...));
  }
}

using default_core_numbers_visitor = core_numbers_visitor<>;

namespace core_numbers_detail {

// the core numbers start as the indegree or inweight.  This function
// will initialize these values
template <concepts::IncidenceGraph G, concepts::ReadWriteVertexPropertyMap<G> CoreMap,
          concepts::ReadableEdgePropertyMap<G> EdgeWeightMap>
requires concepts::VertexListGraph<G>
void compute_in_degree_map(const G& g, CoreMap d, EdgeWeightMap wm) {
  for (auto v : vertices(g)) {
    put(d, v, 0);
  }
  for (auto v : vertices(g)) {
    for (auto e : out_edges(v, g)) {
      put(d, target(e, g), get(d, target(e, g)) + get(wm, e));
    }
  }
}

// the version for weighted graphs is a little different
template <concepts::IncidenceGraph G, concepts::ReadWriteVertexPropertyMap<G> CoreMap,
          concepts::ReadableEdgePropertyMap<G> EdgeWeightMap, typename MutableQueue, concepts::CoreNumbersVisitor<G> V>
property_traits_value_t<CoreMap> core_numbers_impl(const G& g, CoreMap c, EdgeWeightMap wm, MutableQueue& q, V vis) {
  property_traits_value_t<CoreMap> v_cn = {};
  while (!q.empty()) {
    // remove v from the q, and then decrease the core numbers
    // of its successors
    auto v = q.top();
    vis.examine_vertex(v, g);
    q.pop();
    v_cn = get(c, v);
    for (auto e : out_edges(v, g)) {
      vis.examine_edge(e, g);
      auto u = target(e, g);
      // if c[u] > c[v], then u is still in the graph,
      if (get(c, u) > v_cn) {
        // remove the edge
        put(c, u, get(c, u) - get(wm, e));
        if (q.contains(u)) {
          q.update(u);
        }
      }
    }
    vis.finish_vertex(v, g);
  }
  return v_cn;
}

template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> CoreMap,
          concepts::ReadableEdgePropertyMap<G> EdgeWeightMap, concepts::ReadWriteVertexPropertyMap<G> IndexMap,
          concepts::CoreNumbersVisitor<G> V>
property_traits_value_t<CoreMap> core_numbers_dispatch(const G& g, CoreMap c, EdgeWeightMap wm, IndexMap im, V vis) {
  // build the mutable queue
  using Vertex = graph_vertex_descriptor_t<G>;
  auto index_in_heap = make_vector_property_map<std::size_t>(im);
  auto q = make_d_ary_heap_indirect<Vertex, 4>(c, index_in_heap.ref());

  for (auto v : vertices(g)) {
    q.push(v);
  }
  return core_numbers_impl(g, c, wm, q, vis);
}

// the version for the unweighted case
// for this functions CoreMap must be initialized
// with the in degree of each vertex
template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> CoreMap,
          concepts::ReadWriteVertexPropertyMap<G> PositionMap, concepts::CoreNumbersVisitor<G> V>
property_traits_value_t<CoreMap> core_numbers_impl(const G& g, CoreMap c, PositionMap pos, V vis) {
  using Vertex = graph_vertex_descriptor_t<G>;
  using CoreValue = property_traits_value_t<CoreMap>;

  // store the vertex core numbers
  CoreValue v_cn = {};

  // compute the maximum degree (degrees are in the coremap)
  CoreValue max_deg = {};
  for (auto v : vertices(g)) {
    max_deg = std::max(max_deg, get(c, v));
  }

  // store the vertices in bins by their degree
  // allocate two extra locations to ease boundary cases
  std::vector<std::size_t> bin(max_deg + 2);
  for (auto v : vertices(g)) {
    ++bin[get(c, v)];
  }

  // this loop sets bin[d] to the starting position of vertices
  // with degree d in the vert array for the bucket sort
  std::size_t cur_pos = 0;
  for (std::size_t cur_deg = 0; cur_deg < max_deg + 2; ++cur_deg) {
    std::size_t tmp = bin[cur_deg];
    bin[cur_deg] = cur_pos;
    cur_pos += tmp;
  }

  // perform the bucket sort with pos and vert so that
  // pos[0] is the vertex of smallest degree
  std::vector<Vertex> vert(num_vertices(g));
  for (auto v : vertices(g)) {
    std::size_t p = bin[get(c, v)];
    put(pos, v, p);
    vert[p] = v;
    ++bin[get(c, v)];
  }
  // we abused bin while placing the vertices, now,
  // we need to restore it
  std::copy_backward(bin.begin(), bin.end() - 2, bin.end() - 1);
  // now simulate removing the vertices
  for (std::size_t i = 0; i < num_vertices(g); ++i) {
    Vertex v = vert[i];
    vis.examine_vertex(v, g);
    v_cn = get(c, v);
    for (auto e : out_edges(v, g)) {
      vis.examine_edge(e, g);
      Vertex u = target(e, g);
      // if c[u] > c[v], then u is still in the graph,
      if (get(c, u) > v_cn) {
        std::size_t deg_u = get(c, u);
        std::size_t pos_u = get(pos, u);
        // w is the first vertex with the same degree as u
        // (this is the resort operation!)
        std::size_t pos_w = bin[deg_u];
        Vertex w = vert[pos_w];
        if (u != v) {
          // swap u and w
          put(pos, u, pos_w);
          put(pos, w, pos_u);
          vert[pos_w] = u;
          vert[pos_u] = w;
        }
        // now, the vertices array is sorted assuming
        // we perform the following step
        // start the set of vertices with degree of u
        // one into the future (this now points at vertex
        // w which we swapped with u).
        ++bin[deg_u];
        // we are removing v from the graph, so u's degree
        // decreases
        put(c, u, get(c, u) - 1);
      }
    }
    vis.finish_vertex(v, g);
  }
  return v_cn;
}

}  // namespace core_numbers_detail

// non-named parameter version for the unweighted case
template <concepts::IncidenceGraph G, concepts::ReadWriteVertexPropertyMap<G> CoreMap,
          concepts::CoreNumbersVisitor<G> V>
property_traits_value_t<CoreMap> core_numbers(const G& g, CoreMap c, V vis) {
  core_numbers_detail::compute_in_degree_map(g, c, single_property_map<property_traits_value_t<CoreMap>>(1));
  return core_numbers_detail::core_numbers_impl(
      g, c, make_vector_property_map<std::size_t>(num_vertices(g), get(vertex_index, g)).ref(), vis);
}

// non-named paramter version for the unweighted case
template <concepts::IncidenceGraph G, concepts::ReadWriteVertexPropertyMap<G> CoreMap>
property_traits_value_t<CoreMap> core_numbers(const G& g, CoreMap c) {
  return core_numbers(g, c, make_core_numbers_visitor(null_visitor()));
}

// non-named parameter version for the weighted case
template <concepts::IncidenceGraph G, concepts::ReadWriteVertexPropertyMap<G> CoreMap,
          concepts::ReadableEdgePropertyMap<G> EdgeWeightMap, concepts::ReadWriteVertexPropertyMap<G> VertexIndexMap,
          concepts::CoreNumbersVisitor<G> V>
requires concepts::VertexListGraph<G> property_traits_value_t<CoreMap> core_numbers(const G& g, CoreMap c,
                                                                                    EdgeWeightMap wm,
                                                                                    VertexIndexMap vim, V vis) {
  core_numbers_detail::compute_in_degree_map(g, c, wm);
  return core_numbers_detail::core_numbers_dispatch(g, c, wm, vim, vis);
}

// non-named parameter version for the weighted case
//    template <typename Graph, typename CoreMap, typename EdgeWeightMap>
//    property_traits_value_t<CoreMap>
//    core_numbers(Graph& g, CoreMap c, EdgeWeightMap wm)
//    {
//        core_numbers_detail::compute_in_degree_map(g,c,wm);
//        return core_numbers_detail::core_numbers_dispatch(g,c,wm,get(vertex_index,g),
//            make_core_numbers_visitor(null_visitor()));
//    }

template <concepts::IncidenceGraph G, concepts::ReadWriteVertexPropertyMap<G> CoreMap,
          concepts::CoreNumbersVisitor<G> V>
property_traits_value_t<CoreMap> weighted_core_numbers(const G& g, CoreMap c, V vis) {
  return core_numbers(g, c, get(edge_weight, g), get(vertex_index, g), vis);
}

template <concepts::IncidenceGraph G, concepts::ReadWriteVertexPropertyMap<G> CoreMap>
property_traits_value_t<CoreMap> weighted_core_numbers(const G& g, CoreMap c) {
  return weighted_core_numbers(g, c, make_core_numbers_visitor(null_visitor()));
}

}  // namespace bagl

#endif  // BAGL_BAGL_CORE_NUMBERS_H_
