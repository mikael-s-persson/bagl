// Copyright (C) 2005-2009 Jongsoo Park <jongsoo.park -at- gmail.com>
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DOMINATOR_H_
#define BAGL_BAGL_DOMINATOR_H_

#include <deque>
#include <set>

#include "bagl/depth_first_search.h"
#include "bagl/graph_concepts.h"
#include "bagl/vector_property_map.h"

// Dominator tree computation

namespace bagl {
namespace detail {
// An extended time_stamper which also records vertices for each dfs number
template <class TimeMap, class VertexVector, class TimeT>
class time_stamper_with_vertex_on_discover_vertex {
 public:
  time_stamper_with_vertex_on_discover_vertex(TimeMap time_map, VertexVector& v, TimeT& t)
      : time_pa_(time_map), time_(t), v_(v) {}

  template <class Graph>
  void discover_vertex(const property_traits_key_t<TimeMap>& v, const Graph& g) {
    put(time_pa_, v, ++time_);
    v_[time_] = v;
  }

 private:
  TimeMap time_pa_;
  TimeT& time_;
  VertexVector& v_;
};

template <concepts::BidirectionalGraph Graph, concepts::ReadableVertexIndexMap<Graph> IndexMap,
          concepts::ReadWriteVertexPropertyMap<Graph> TimeMap, concepts::ReadWriteVertexPropertyMap<Graph> PredMap,
          concepts::WritableVertexPropertyMap<Graph> DomTreePredMap>
class dominator_visitor {
  using Vertex = graph_vertex_descriptor_t<Graph>;

 public:
  // `g` [in] the target graph of the dominator tree
  // `entry` [in] the entry node of g
  // `index_map` [in] the vertex index map for g
  // `dom_tree_pred_map` [out] the immediate dominator map (parent map in dominator tree)
  dominator_visitor(const Graph& g, const Vertex& entry, const IndexMap& index_map, DomTreePredMap dom_tree_pred_map)
      : semi_map_(num_vertices(g), index_map, graph_traits<Graph>::null_vertex()),
        ancestor_map_(num_vertices(g), index_map, graph_traits<Graph>::null_vertex()),
        best_map_(num_vertices(g), index_map, graph_traits<Graph>::null_vertex()),
        bucket_map_(num_vertices(g), index_map),
        entry_(entry),
        dom_tree_pred_map_(dom_tree_pred_map),
        num_of_vertices_(num_vertices(g)),
        samedom_map(num_vertices(g), index_map, graph_traits<Graph>::null_vertex()) {}

  void operator()(const Vertex& n, const TimeMap& df_num_map, const PredMap& parent_map, const Graph& g) {
    if (n == entry_) {
      return;
    }

    const Vertex p(get(parent_map, n));
    Vertex s(p);

    // 1. Calculate the semidominator of n,
    // based on the semidominator thm.
    // * Semidominator thm. : To find the semidominator of a node n,
    //   consider all predecessors v of n in the CFG (Control Flow
    //   Graph).
    //  - If v is a proper ancestor of n in the spanning tree
    //    (so dfnum(v) < dfnum(n)), then v is a candidate for semi(n)
    //  - If v is a non-ancestor of n (so dfnum(v) > dfnum(n))
    //    then for each u that is an ancestor of v (or u = v),
    //    Let semi(u) be a candidate for semi(n)
    //   of all these candidates, the one with lowest dfnum is
    //   the semidominator of n.

    // For each predecessor of n
    for (auto e : in_edges(n, g)) {
      const Vertex v = source(e, g);
      // To deal with unreachable nodes
      if (get(df_num_map, v) < 0 || get(df_num_map, v) >= num_of_vertices_) {
        continue;
      }

      Vertex s2;
      if (get(df_num_map, v) <= get(df_num_map, n)) {
        s2 = v;
      } else {
        s2 = get(semi_map_, ancestor_with_lowest_semi(v, df_num_map));
      }

      if (get(df_num_map, s2) < get(df_num_map, s)) {
        s = s2;
      }
    }
    put(semi_map_, n, s);

    // 2. Calculation of n's dominator is deferred until
    // the path from s to n has been linked into the forest
    get(bucket_map_, s).push_back(n);
    get(ancestor_map_, n) = p;
    get(best_map_, n) = n;

    // 3. Now that the path from p to v has been linked into
    // the spanning forest, these lines calculate the dominator of v,
    // based on the dominator thm., or else defer the calculation
    // until y's dominator is known
    // * Dominator thm. : On the spanning-tree path below semi(n) and
    //   above or including n, let y be the node
    //   with the smallest-numbered semidominator. Then,
    //
    //  idom(n) = semi(n) if semi(y)=semi(n) or
    //            idom(y) if semi(y) != semi(n)
    for (auto v : get(bucket_map_, p)) {
      const Vertex y(ancestor_with_lowest_semi(v, df_num_map));
      if (get(semi_map_, y) == get(semi_map_, v)) {
        put(dom_tree_pred_map_, v, p);
      } else {
        put(samedom_map, v, y);
      }
    }

    get(bucket_map_, p).clear();
  }

 protected:
  // Evaluate function in Tarjan's path compression
  Vertex ancestor_with_lowest_semi(const Vertex& v, const TimeMap& df_num_map) {
    const Vertex a(get(ancestor_map_, v));

    if (get(ancestor_map_, a) != graph_traits<Graph>::null_vertex()) {
      const Vertex b(ancestor_with_lowest_semi(a, df_num_map));

      put(ancestor_map_, v, get(ancestor_map_, a));

      if (get(df_num_map, get(semi_map_, b)) < get(df_num_map, get(semi_map_, get(best_map_, v)))) {
        put(best_map_, v, b);
      }
    }

    return get(best_map_, v);
  }

  vector_property_map<Vertex, IndexMap> semi_map_;
  vector_property_map<Vertex, IndexMap> ancestor_map_;
  vector_property_map<Vertex, IndexMap> best_map_;
  vector_property_map<std::deque<Vertex>, IndexMap> bucket_map_;

  const Vertex& entry_;
  DomTreePredMap dom_tree_pred_map_;
  const std::size_t num_of_vertices_;

 public:
  PredMap samedom_map;
};

}  // namespace detail

// Build dominator tree using Lengauer-Tarjan algorithm.
// It takes O((V+E)log(V+E)) time.
//
// df_num_map, parent_map and vertices_by_df_num have dfs results corresponding index_map.
// If dfs has already run before, this function would be good for saving computations.
//
// Unreachable nodes must be masked as graph_traits<Graph>::null_vertex in parent_map.
//
// Unreachable nodes must be masked as std::numeric_limits<std::size_t>::max() in df_num_map.
//
// `dom_tree_pred_map` [out] : immediate dominator map (parent map in dom. tree)
//
// Reference Appel. p. 452~453. algorithm 19.9, 19.10.
//
// TODO: Optimization in Finding Dominators in Practice, Loukas Georgiadis
template <concepts::BidirectionalGraph Graph, concepts::ReadableVertexIndexMap<Graph> IndexMap,
          concepts::ReadWriteVertexPropertyMap<Graph> TimeMap, concepts::ReadWriteVertexPropertyMap<Graph> PredMap,
          class VertexVector, concepts::WritableVertexPropertyMap<Graph> DomTreePredMap>
void lengauer_tarjan_dominator_tree_without_dfs(const Graph& g, const graph_vertex_descriptor_t<Graph>& entry,
                                                const IndexMap& index_map, TimeMap df_num_map, PredMap parent_map,
                                                VertexVector& vertices_by_df_num, DomTreePredMap dom_tree_pred_map) {
  // Typedefs and concept check
  using Vertex = graph_vertex_descriptor_t<Graph>;

  const std::size_t num_of_vertices = num_vertices(g);
  if (num_of_vertices == 0) {
    return;
  }

  // 1. Visit each vertex in reverse post order and calculate sdom.
  detail::dominator_visitor<Graph, IndexMap, TimeMap, PredMap, DomTreePredMap> visitor(g, entry, index_map,
                                                                                       dom_tree_pred_map);

  for (std::size_t i = 0; i < num_of_vertices; ++i) {
    const Vertex u(vertices_by_df_num[num_of_vertices - 1 - i]);
    if (u != graph_traits<Graph>::null_vertex()) {
      visitor(u, df_num_map, parent_map, g);
    }
  }

  // 2. Now all the deferred dominator calculations,
  // based on the second clause of the dominator thm., are performed
  for (std::size_t i = 0; i < num_of_vertices; ++i) {
    const Vertex n(vertices_by_df_num[i]);

    if (n == entry || n == graph_traits<Graph>::null_vertex()) {
      continue;
    }

    Vertex u = get(visitor.samedom_map, n);
    if (u != graph_traits<Graph>::null_vertex()) {
      put(dom_tree_pred_map, n, get(dom_tree_pred_map, u));
    }
  }
}

// Unlike lengauer_tarjan_dominator_tree_without_dfs, dfs is run in this function and
// the result is written to df_num_map, parent_map, vertices.
//
// If the result of dfs required after this algorithm, this function can eliminate the need of rerunning dfs.
template <concepts::BidirectionalGraph Graph, concepts::ReadableVertexIndexMap<Graph> IndexMap, class TimeMap,
          concepts::ReadWriteVertexPropertyMap<Graph> PredMap, class VertexVector,
          concepts::WritableVertexPropertyMap<Graph> DomTreePredMap>
void lengauer_tarjan_dominator_tree(const Graph& g, const graph_vertex_descriptor_t<Graph>& entry,
                                    const IndexMap& index_map, TimeMap df_num_map, PredMap parent_map,
                                    VertexVector& vertices_by_df_num, DomTreePredMap dom_tree_pred_map) {
  // 1. Depth first visit
  const std::size_t num_of_vertices = num_vertices(g);
  if (num_of_vertices == 0) {
    return;
  }

  auto time = std::numeric_limits<std::size_t>::max();
  depth_first_visit(
      g, entry,
      make_dfs_visitor(predecessor_recorder_on_tree_edge(parent_map),
                       detail::time_stamper_with_vertex_on_discover_vertex(df_num_map, vertices_by_df_num, time)),
      vector_property_map(num_of_vertices, index_map, color_traits<default_color_type>::white()));

  // 2. Run main algorithm.
  lengauer_tarjan_dominator_tree_without_dfs(g, entry, index_map, df_num_map, parent_map, vertices_by_df_num,
                                             dom_tree_pred_map);
}

// Use vertex_index as IndexMap and make df_num_map, parent_map, vertices_by_df_num internally.
// If we don't need the result of dfs (df_num_map, parent_map, vertices_by_df_num),
// this function would be more convenient one.
template <concepts::BidirectionalGraph Graph, concepts::WritableVertexPropertyMap<Graph> DomTreePredMap>
void lengauer_tarjan_dominator_tree(const Graph& g, const graph_vertex_descriptor_t<Graph>& entry,
                                    DomTreePredMap dom_tree_pred_map) {
  // typedefs
  using Vertex = graph_vertex_descriptor_t<Graph>;

  // Make property maps
  const std::size_t num_of_vertices = num_vertices(g);
  if (num_of_vertices == 0) {
    return;
  }

  auto index_map = get(vertex_index, g);
  auto df_num_map = vector_property_map(num_of_vertices, index_map, std::size_t{0});
  auto parent_map = vector_property_map(num_of_vertices, index_map, graph_traits<Graph>::null_vertex());
  std::vector<Vertex> vertices_by_df_num(num_of_vertices, graph_traits<Graph>::null_vertex());

  // Run main algorithm
  lengauer_tarjan_dominator_tree(g, entry, index_map, df_num_map, parent_map, vertices_by_df_num, dom_tree_pred_map);
}

// Muchnick. p. 182, 184
//
// Using iterative bit vector analysis
template <concepts::BidirectionalGraph Graph, concepts::ReadableVertexIndexMap<Graph> IndexMap,
          concepts::WritableVertexPropertyMap<Graph> DomTreePredMap>
void iterative_bit_vector_dominator_tree(const Graph& g, const graph_vertex_descriptor_t<Graph>& entry,
                                         const IndexMap& index_map, DomTreePredMap dom_tree_pred_map) {
  using Vertex = graph_vertex_descriptor_t<Graph>;

  // 1. Finding dominator
  // 1.1. Initialize
  const std::size_t num_of_vertices = num_vertices(g);
  if (num_of_vertices == 0) {
    return;
  }

  auto v_rg = vertices(g);
  const std::set<Vertex> N(v_rg.begin(), v_rg.end());

  bool change = true;

  auto dom_map = vector_property_map(num_of_vertices, index_map, N);
  get(dom_map, entry).clear();
  get(dom_map, entry).insert(entry);

  while (change) {
    change = false;
    for (auto v : vertices(g)) {
      if (v == entry) {
        continue;
      }

      std::set<Vertex> T(N);

      for (auto e : in_edges(v, g)) {
        const Vertex p = source(e, g);

        std::set<Vertex> temp_set;
        std::set_intersection(T.begin(), T.end(), get(dom_map, p).begin(), get(dom_map, p).end(),
                              std::inserter(temp_set, temp_set.begin()));
        T.swap(temp_set);
      }

      T.insert(v);
      if (T != get(dom_map, v)) {
        change = true;
        get(dom_map, v).swap(T);
      }
    }  // end of for (boost::tie(vi, viend) = vertices(g)
  }    // end of while(change)

  // 2. Build dominator tree
  for (auto v : vertices(g)) {
    get(dom_map, v).erase(v);
  }

  Graph dom_tree(num_of_vertices);

  for (auto v : vertices(g)) {
    if (v == entry) {
      continue;
    }

    // We have to iterate through copied dominator set
    auto& dom_v = get(dom_map, v);
    const std::set<Vertex> temp_set = dom_v;
    for (auto s : temp_set) {
      const auto& dom_s = get(dom_map, s);
      for (auto t = dom_v.begin(); t != dom_v.end();) {
        auto old_t = t;
        ++t;  // Done early because t may become invalid
        if (*old_t == s) {
          continue;
        }
        if (dom_s.find(*old_t) != dom_s.end()) {
          dom_v.erase(old_t);
        }
      }
    }
  }

  for (auto v : vertices(g)) {
    if (v != entry && get(dom_map, v).size() == 1) {
      Vertex temp = *get(dom_map, v).begin();
      put(dom_tree_pred_map, v, temp);
    }
  }
}

template <concepts::BidirectionalGraph Graph, concepts::WritableVertexPropertyMap<Graph> DomTreePredMap>
void iterative_bit_vector_dominator_tree(const Graph& g, const graph_vertex_descriptor_t<Graph>& entry,
                                         DomTreePredMap dom_tree_pred_map) {
  iterative_bit_vector_dominator_tree(g, entry, get(vertex_index, g), dom_tree_pred_map);
}
}  // namespace bagl

#endif  // BAGL_BAGL_DOMINATOR_H_
