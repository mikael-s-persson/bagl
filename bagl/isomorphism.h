// Copyright (C) 2001 Jeremy Siek, Douglas Gregor, Brian Osman
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_ISOMORPHISM_H_
#define BAGL_BAGL_ISOMORPHISM_H_

#include <algorithm>
#include <cassert>
#include <iterator>
#include <memory>
#include <ranges>
#include <unordered_map>
#include <utility>
#include <vector>

#include "bagl/depth_first_search.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/vector_property_map.h"

namespace bagl {

namespace isomorphism_detail {
template <typename G, typename Inv>
using default_invariant_count_map =
    std::unordered_map<decltype(std::declval<Inv>()(std::declval<graph_vertex_descriptor_t<G>>())), std::size_t>;

// When you can't use goto, use an algorithmic class and pretend you're not committing a sin.
// If you just want to watch the world burn, use goto's inside an algorithmic class.
template <typename Graph1, typename Graph2, typename IsoMapping, typename Invariant1, typename Invariant2,
          typename IndexMap1, typename IndexMap2,
          typename InvariantCountMap = default_invariant_count_map<Graph1, Invariant1>>
class isomorphism_algo {
  using Vertex1 = graph_vertex_descriptor_t<Graph1>;
  using Vertex2 = graph_vertex_descriptor_t<Graph2>;
  using Edge1 = graph_edge_descriptor_t<Graph1>;
  using InvariantValue = decltype(std::declval<Invariant1>()(std::declval<Vertex1>()));

  const Graph1& g1_;
  const Graph2& g2_;
  IsoMapping f_;
  Invariant1 invariant1_;
  Invariant2 invariant2_;
  IndexMap1 index_map1_;
  IndexMap2 index_map2_;

  std::vector<Vertex1> dfs_vertices_;
  vector_property_map<int, IndexMap1> dfs_num_;
  std::vector<Edge1> ordered_edges_;

  vector_property_map<bool, IndexMap2> in_s_;

  int num_edges_on_k_ = 0;

  struct record_dfs_order : default_dfs_visitor {
    record_dfs_order(std::vector<Vertex1>& v, std::vector<Edge1>& e) : vertices(&v), edges(&e) {}

    void discover_vertex(Vertex1 v, const Graph1& /*unused*/) const { vertices->push_back(v); }
    void examine_edge(Edge1 e, const Graph1& /*unused*/) const { edges->push_back(e); }
    std::vector<Vertex1>* vertices;
    std::vector<Edge1>* edges;
  };

 public:
  isomorphism_algo(const Graph1& g1, const Graph2& g2, IsoMapping f, Invariant1 invariant1, Invariant2 invariant2,
                   std::size_t /* max_invariant */, IndexMap1 index_map1, IndexMap2 index_map2)
      : g1_(g1),
        g2_(g2),
        f_(f),
        invariant1_(std::move(invariant1)),
        invariant2_(std::move(invariant2)),
        index_map1_(index_map1),
        index_map2_(index_map2),
        in_s_(vector_property_map(num_vertices(g1_), index_map2_, false)) {}

  // Generates map of invariant multiplicity from sorted invariants
  template <std::ranges::forward_range ForwardRange>
  InvariantCountMap multiplicities(const ForwardRange& rg) {
    assert(std::ranges::is_sorted(rg));
    InvariantCountMap invar_multiplicity;

    if (rg.empty()) {
      return invar_multiplicity;
    }

    auto first = rg.begin();
    InvariantValue invar = *first;
    auto inserted = invar_multiplicity.emplace(invar, 1).first;
    for (++first; first != rg.end(); ++first) {
      if (*first == invar) {
        inserted->second += 1;
      } else {
        invar = *first;
        inserted = invar_multiplicity.emplace(invar, 1).first;
      }
    }

    return invar_multiplicity;
  }

  bool test_isomorphism() {
    // reset isomapping
    for (auto v : vertices(g1_)) {
      f_[v] = graph_traits<Graph2>::null_vertex();
    }

    // Calculate all invariants of g1_ and g2_, sort and compare
    std::vector<InvariantValue> invar1_array;
    invar1_array.reserve(num_vertices(g1_));
    for (auto v : vertices(g1_)) {
      invar1_array.push_back(invariant1_(v));
    }
    std::ranges::sort(invar1_array);

    std::vector<InvariantValue> invar2_array;
    invar2_array.reserve(num_vertices(g2_));
    for (auto v : vertices(g2_)) {
      invar2_array.push_back(invariant2_(v));
    }
    std::ranges::sort(invar2_array);
    if (invar1_array != invar2_array) {
      return false;
    }

    // Sort vertices by the multiplicity of their invariants
    std::vector<Vertex1> v_mult;
    for (auto v : vertices(g1_)) {
      v_mult.push_back(v);
    }
    std::ranges::sort(v_mult, [this, multiplicity = multiplicities(invar1_array)](const Vertex1& x, const Vertex1& y) {
      auto x_multiplicity_iter = multiplicity.find(invariant1_(x));
      assert(x_multiplicity_iter != multiplicity.end());
      auto y_multiplicity_iter = multiplicity.find(invariant1_(y));
      assert(y_multiplicity_iter != multiplicity.end());
      return *x_multiplicity_iter < *y_multiplicity_iter;
    });

    using Color = color_traits<default_color_type>;
    auto color_map = vector_property_map(num_vertices(g1_), index_map1_, Color::white());
    record_dfs_order dfs_visitor(dfs_vertices_, ordered_edges_);
    for (auto u : v_mult) {
      if (color_map[u] == Color::white()) {
        dfs_visitor.start_vertex(u, g1_);
        depth_first_visit(g1_, u, dfs_visitor, color_map.ref());
      }
    }
    // Create the dfs_num array and dfs_num_map
    dfs_num_ = vector_property_map(num_vertices(g1_), index_map1_, int{0});
    std::size_t n = 0;
    for (auto v : dfs_vertices_) {
      dfs_num_[v] = n++;
    }

    std::ranges::sort(ordered_edges_, [this](const Edge1& e1, const Edge1& e2) {
      int u1 = dfs_num_[source(e1, g1_)];
      int v1 = dfs_num_[target(e1, g1_)];
      int u2 = dfs_num_[source(e2, g1_)];
      int v2 = dfs_num_[target(e2, g1_)];
      int m1 = std::max(u1, v1);
      int m2 = std::max(u2, v2);
      // lexicographical comparison
      return std::tuple(m1, u1, v1) < std::tuple(m2, u2, v2);
    });

    int dfs_num_k = -1;
    return match(dfs_num_k);
  }

 private:
  struct match_continuation {
    enum { pos_g2_vertex_loop, pos_fi_adj_loop, pos_dfs_num } position;
    std::optional<decltype(partial_view(vertices(std::declval<Graph2>())))> g2_verts;
    std::optional<decltype(partial_view(adjacent_vertices(std::declval<Vertex2>(), std::declval<Graph2>())))> fi_adj;
    typename std::vector<Edge1>::iterator iter;
    int dfs_num_k;
  };

  bool match(int dfs_num_k) {
    std::vector<match_continuation> k;
    auto iter = ordered_edges_.begin();
    decltype(match_continuation::g2_verts) g2_verts = partial_view(vertices(g2_));
    decltype(match_continuation::fi_adj) fi_adj;
    Vertex1 i;
    Vertex1 j;

  recur:
    if (iter != ordered_edges_.end()) {
      i = source(*iter, g1_);
      j = target(*iter, g1_);
      if (dfs_num_[i] > dfs_num_k) {
        g2_verts = partial_view(vertices(g2_));
        while (!g2_verts->empty()) {
          {
            Vertex2 u = *g2_verts->begin();
            Vertex1 kp1 = dfs_vertices_[dfs_num_k + 1];
            if (invariant1_(kp1) == invariant2_(u) && !in_s_[u]) {
              f_[kp1] = u;
              in_s_[u] = true;
              num_edges_on_k_ = 0;

              match_continuation new_k;
              new_k.position = match_continuation::pos_g2_vertex_loop;
              new_k.g2_verts = *g2_verts;
              new_k.iter = iter;
              new_k.dfs_num_k = dfs_num_k;
              k.push_back(new_k);
              ++dfs_num_k;
              goto recur;
            }
          }
        g2_loop_k:
          g2_verts->move_to_next();
        }
      } else if (dfs_num_[j] > dfs_num_k) {
        {
          Vertex1 vk = dfs_vertices_[dfs_num_k];
          num_edges_on_k_ -= std::ranges::count_if(adjacent_vertices(f_[vk], g2_), [this](auto u) { return in_s_[u]; });

          for (int jj = 0; jj < dfs_num_k; ++jj) {
            num_edges_on_k_ -= std::ranges::count(adjacent_vertices(f_[dfs_vertices_[jj]], g2_), f_[vk]);
          }
        }

        if (num_edges_on_k_ != 0) {
          goto return_point_false;
        }
        fi_adj = partial_view(adjacent_vertices(f_[i], g2_));
        while (!fi_adj->empty()) {
          {
            Vertex2 v = *fi_adj->begin();
            if (invariant2_(v) == invariant1_(j) && !in_s_[v]) {
              f_[j] = v;
              in_s_[v] = true;
              num_edges_on_k_ = 1;
              int next_k = std::max(dfs_num_k, std::max(dfs_num_[i], dfs_num_[j]));
              match_continuation new_k;
              new_k.position = match_continuation::pos_fi_adj_loop;
              new_k.fi_adj = *fi_adj;
              new_k.iter = iter;
              new_k.dfs_num_k = dfs_num_k;
              ++iter;
              dfs_num_k = next_k;
              k.push_back(new_k);
              goto recur;
            }
          }
        fi_adj_loop_k:
          fi_adj->move_to_next();
        }
      } else {
        if (std::ranges::count(adjacent_vertices(f_[i], g2_), f_[j]) != 0) {
          ++num_edges_on_k_;
          match_continuation new_k;
          new_k.position = match_continuation::pos_dfs_num;
          k.push_back(new_k);
          ++iter;
          goto recur;
        }
      }
    } else {
      // At this point, there may still be null vertices in the
      // mapping for disconnected vertices
      map_disconnected_vertices();
      return true;
    }

    {
    return_point_false:
      if (k.empty()) {
        return false;
      }
      const match_continuation& this_k = k.back();
      switch (this_k.position) {
        case match_continuation::pos_g2_vertex_loop: {
          g2_verts = *this_k.g2_verts;
          iter = this_k.iter;
          dfs_num_k = this_k.dfs_num_k;
          k.pop_back();
          in_s_[*g2_verts->begin()] = false;
          i = source(*iter, g1_);
          j = target(*iter, g1_);
          goto g2_loop_k;
        }
        case match_continuation::pos_fi_adj_loop: {
          fi_adj = *this_k.fi_adj;
          iter = this_k.iter;
          dfs_num_k = this_k.dfs_num_k;
          k.pop_back();
          in_s_[*fi_adj->begin()] = false;
          i = source(*iter, g1_);
          j = target(*iter, g1_);
          goto fi_adj_loop_k;
        }
        case match_continuation::pos_dfs_num: {
          k.pop_back();
          goto return_point_false;
        }
        default: {
          assert(!"Bad position");
          abort();
        }
      }
    }
  }

  void map_disconnected_vertices() {
    std::vector<Vertex1> unmatched_g1_vertices;
    for (auto v : vertices(g1_)) {
      if (f_[v] == graph_traits<Graph2>::null_vertex()) {
        unmatched_g1_vertices.push_back(v);
      }
    }

    if (!unmatched_g1_vertices.empty()) {
      std::unordered_multimap<InvariantValue, Vertex2> unmatched_invariants;
      for (auto v : vertices(g2_)) {
        if (!in_s_[v]) {
          unmatched_invariants.emplace(invariant2_(v), v);
        }
      }

      for (auto v : unmatched_g1_vertices) {
        InvariantValue unmatched_g1_vertex_invariant = invariant1_(v);
        auto matching_invariant = unmatched_invariants.find(unmatched_g1_vertex_invariant);
        assert(matching_invariant != unmatched_invariants.end());
        f_[v] = matching_invariant->second;
        unmatched_invariants.erase(matching_invariant);
      }
    }
  }
};

template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> InDegreeMap>
void compute_in_degree(const G& g, InDegreeMap in_degree_map) {
  for (auto v : vertices(g)) {
    put(in_degree_map, v, 0);
  }

  for (auto u : vertices(g)) {
    for (auto v : adjacent_vertices(u, g)) {
      put(in_degree_map, v, get(in_degree_map, v) + 1);
    }
  }
}

}  // namespace isomorphism_detail

template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> InDegreeMap>
class degree_vertex_invariant {
 public:
  degree_vertex_invariant(const G& g, InDegreeMap&& in_degree_map) : g_(&g), in_degree_map_(std::move(in_degree_map)) {
    for (auto v : vertices(*g_)) {
      max_vertex_in_degree_ = std::max(max_vertex_in_degree_, get(in_degree_map_, v));
      max_vertex_out_degree_ = std::max(max_vertex_out_degree_, out_degree(v, *g_));
    }
  }

  std::size_t operator()(graph_vertex_descriptor_t<G> v) {
    return (max_vertex_in_degree_ + 1) * out_degree(v, *g_) + get(in_degree_map_, v);
  }
  // The largest possible vertex invariant number
  [[nodiscard]] std::size_t max() const { return (max_vertex_in_degree_ + 1) * (max_vertex_out_degree_ + 1); }

 private:
  const G* g_;
  InDegreeMap in_degree_map_;
  std::size_t max_vertex_in_degree_ = 0;
  std::size_t max_vertex_out_degree_ = 0;
};

// Count actual number of vertices, even in filtered graphs.
template <concepts::VertexListGraph Graph>
auto count_vertices(const Graph& g) {
  return std::ranges::distance(vertices(g));
}

template <concepts::VertexListGraph G, concepts::ReadableVertexIndexMap<G> Index>
auto make_degree_invariant(const G& g, const Index& index) {
  auto pm = vector_property_map(num_vertices(g), index, std::size_t{0});
  isomorphism_detail::compute_in_degree(g, pm.ref());
  return degree_vertex_invariant(g, std::move(pm));
}

template <concepts::VertexAndEdgeListGraph Graph1, concepts::VertexListGraph Graph2,
          concepts::ReadWriteVertexPropertyMap<Graph1> IsoMapping,
          std::invocable<graph_vertex_descriptor_t<Graph1>> Invariant1,
          std::invocable<graph_vertex_descriptor_t<Graph2>> Invariant2,
          concepts::ReadableVertexIndexMap<Graph1> IndexMap1, concepts::ReadableVertexIndexMap<Graph2> IndexMap2>
bool isomorphism(
    const Graph1& g1, const Graph2& g2,
    IsoMapping f,           // vector_property_map(num_vertices(g1), index_map1, graph_traits<Graph2>::null_vertex())
    Invariant1 invariant1,  // make_degree_invariant(g1, index_map1)
    Invariant2 invariant2,  // make_degree_invariant(g2, index_map2)
    std::size_t /* max_invariant */,
    IndexMap1 index_map1,  // get(vertex_index, g1)
    IndexMap2 index_map2)  // get(vertex_index, g2)

{
  using Vertex1 = graph_vertex_descriptor_t<Graph1>;
  using Vertex2 = graph_vertex_descriptor_t<Graph2>;

  using Invariant1Value = decltype(invariant1(std::declval<Vertex1>()));
  using Invariant2Value = decltype(invariant2(std::declval<Vertex2>()));
  static_assert(std::is_same_v<Invariant1Value, Invariant2Value>);

  // Property map requirements
  using IsoMappingValue = property_traits_value_t<IsoMapping>;
  static_assert(std::is_convertible_v<IsoMappingValue, Vertex2>);

  using IndexMap1Value = property_traits_value_t<IndexMap1>;
  static_assert(std::is_convertible_v<IndexMap1Value, std::size_t>);

  using IndexMap2Value = property_traits_value_t<IndexMap2>;
  static_assert(std::is_convertible_v<IndexMap2Value, std::size_t>);

  if (count_vertices(g1) != count_vertices(g2)) {
    return false;
  }
  if (count_vertices(g1) == 0 && count_vertices(g2) == 0) {
    return true;
  }

  isomorphism_detail::isomorphism_algo algo(g1, g2, f, std::move(invariant1), std::move(invariant2), 0, index_map1, index_map2);
  return algo.test_isomorphism();
}

// Verify that the given mapping iso_map from the vertices of g1 to the
// vertices of g2 describes an isomorphism.
// Note: this could be made much faster by specializing based on the graph
// concepts modeled, but since we're verifying an O(n^(lg n)) algorithm,
// O(n^4) won't hurt us.
template <typename Graph1, typename Graph2, typename IsoMap>
bool verify_isomorphism(const Graph1& g1, const Graph2& g2, IsoMap iso_map) {
  for (auto e1 : edges(g1)) {
    bool found_edge = false;
    for (auto e2 : edges(g2)) {
      if (source(e2, g2) == get(iso_map, source(e1, g1)) && target(e2, g2) == get(iso_map, target(e1, g1))) {
        found_edge = true;
      }
    }

    if (!found_edge) {
      return false;
    }
  }

  return true;
}

}  // namespace bagl

#endif  // BAGL_BAGL_ISOMORPHISM_H_
