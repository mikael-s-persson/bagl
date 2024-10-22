// Copyright (c) 2004 Kristopher Beevers
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_ASTAR_SEARCH_H_
#define BAGL_BAGL_ASTAR_SEARCH_H_

#include <functional>
#include <limits>
#include <vector>

#include "bagl/breadth_first_search.h"
#include "bagl/d_ary_heap.h"
#include "bagl/exception.h"
#include "bagl/function_property_map.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"
#include "bagl/relax.h"
#include "bagl/vector_property_map.h"
#include "bagl/visitors.h"

namespace bagl {

namespace concepts {

template <typename H, typename G>
concept AStarHeuristic = std::copy_constructible<H> && requires(const H& h, graph_vertex_descriptor_t<G> u) {
  { h(u) } -> std::regular;
};

template <typename V, typename G>
concept AStarVisitor = std::copy_constructible<V> &&
    requires(V& vis, const G& g, graph_vertex_descriptor_t<G> u, graph_edge_descriptor_t<G> e) {
  vis.initialize_vertex(u, g);
  vis.discover_vertex(u, g);
  vis.examine_vertex(u, g);
  vis.examine_edge(e, g);
  vis.edge_relaxed(e, g);
  vis.edge_not_relaxed(e, g);
  vis.black_target(e, g);
  vis.finish_vertex(u, g);
};

}  // namespace concepts

template <typename G, typename CostType>
class astar_heuristic {
 public:
  auto operator()(graph_vertex_descriptor_t<G> /*u*/) { return CostType{}; }
};

template <typename Visitors = null_visitors>
class astar_visitor : public bfs_visitor<Visitors> {
 public:
  astar_visitor() = default;
  explicit astar_visitor(Visitors vis) : bfs_visitor<Visitors>(std::move(vis)) {}

  template <typename Edge, typename Graph>
  void edge_relaxed(Edge e, const Graph& g) {
    visitors_detail::invoke_edge_relaxed_on_all(this->vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void edge_not_relaxed(Edge e, const Graph& g) {
    visitors_detail::invoke_edge_not_relaxed_on_all(this->vis_, e, g);
  }

  template <typename Edge, typename Graph>
  void tree_edge(Edge e, const Graph& g) = delete;
  template <typename Edge, typename Graph>
  void non_tree_edge(Edge e, const Graph& g) = delete;
};
template <typename... Visitors>
auto make_astar_visitor(Visitors&&... vis) {
  if constexpr (sizeof...(Visitors) == 0) {
    return astar_visitor<>();
  } else {
    return astar_visitor<std::tuple<std::decay_t<Visitors>...>>(std::forward<Visitors>(vis)...);
  }
}
using default_astar_visitor = astar_visitor<>;

namespace astar_detail {

template <class AStarHeuristic, class UniformCostVisitor, class UpdatableQueue, class PredecessorMap, class CostMap,
          class DistanceMap, class WeightMap, class ColorMap, class BinaryFunction, class BinaryPredicate>
struct astar_bfs_visitor {
  using C = property_traits_value_t<CostMap>;
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;
  using distance_type = property_traits_value_t<DistanceMap>;

  astar_bfs_visitor(AStarHeuristic h, UniformCostVisitor vis, UpdatableQueue& Q, PredecessorMap p, CostMap c,
                    DistanceMap d, WeightMap w, ColorMap col, BinaryFunction combine, BinaryPredicate compare, C zero)
      : h_(h),
        vis_(vis),
        q_(Q),
        predecessor_(p),
        cost_(c),
        distance_(d),
        weight_(w),
        color_(col),
        combine_(combine),
        compare_(compare),
        zero_(zero) {}

  template <class Vertex, class Graph>
  void initialize_vertex(Vertex u, const Graph& g) {
    vis_.initialize_vertex(u, g);
  }
  template <class Vertex, class Graph>
  void discover_vertex(Vertex u, const Graph& g) {
    vis_.discover_vertex(u, g);
  }
  template <class Vertex, class Graph>
  void examine_vertex(Vertex u, const Graph& g) {
    vis_.examine_vertex(u, g);
  }
  template <class Vertex, class Graph>
  void finish_vertex(Vertex u, const Graph& g) {
    vis_.finish_vertex(u, g);
  }
  template <class Edge, class Graph>
  void examine_edge(Edge e, const Graph& g) {
    if (compare_(get(weight_, e), zero_)) {
      throw negative_edge();
    }
    vis_.examine_edge(e, g);
  }
  template <class Edge, class Graph>
  void non_tree_edge(Edge, const Graph&) {}

  template <class Edge, class Graph>
  void tree_edge(Edge e, const Graph& g) {
    bool decreased = relax(e, g, weight_, predecessor_, distance_, combine_, compare_);

    if (decreased) {
      vis_.edge_relaxed(e, g);
      put(cost_, target(e, g), combine_(get(distance_, target(e, g)), h_(target(e, g))));
    } else {
      vis_.edge_not_relaxed(e, g);
    }
  }

  template <class Edge, class Graph>
  void gray_target(Edge e, const Graph& g) {
    bool decreased = relax(e, g, weight_, predecessor_, distance_, combine_, compare_);

    if (decreased) {
      put(cost_, target(e, g), combine_(get(distance_, target(e, g)), h_(target(e, g))));
      q_.update(target(e, g));
      vis_.edge_relaxed(e, g);
    } else {
      vis_.edge_not_relaxed(e, g);
    }
  }

  template <class Edge, class Graph>
  void black_target(Edge e, const Graph& g) {
    bool decreased = relax(e, g, weight_, predecessor_, distance_, combine_, compare_);

    if (decreased) {
      vis_.edge_relaxed(e, g);
      put(cost_, target(e, g), combine_(get(distance_, target(e, g)), h_(target(e, g))));
      q_.emplace(target(e, g));
      put(color_, target(e, g), Color::gray());
      vis_.black_target(e, g);
    } else {
      vis_.edge_not_relaxed(e, g);
    }
  }

  AStarHeuristic h_;
  UniformCostVisitor vis_;
  UpdatableQueue& q_;
  PredecessorMap predecessor_;
  CostMap cost_;
  DistanceMap distance_;
  WeightMap weight_;
  ColorMap color_;
  BinaryFunction combine_;
  BinaryPredicate compare_;
  C zero_;
};

}  // namespace astar_detail

template <concepts::IncidenceGraph G, concepts::AStarHeuristic<G> H, concepts::AStarVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadWriteVertexPropertyMap<G> CostMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, typename VertexIndexMap, typename CompareFunction,
          typename CombineFunction>
void astar_search_no_init(const G& g, graph_vertex_descriptor_t<G> s, H h, V vis, PredecessorMap predecessor,
                          CostMap cost, DistanceMap distance, WeightMap weight, ColorMap color,
                          VertexIndexMap index_map, CompareFunction compare, CombineFunction combine,
                          property_traits_value_t<CostMap> zero) {
  using Vertex = graph_vertex_descriptor_t<G>;
  auto index_in_heap = make_vector_property_map<std::size_t>(index_map);
  auto q = make_d_ary_heap_indirect<Vertex, 4>(cost, index_in_heap, compare);

  astar_detail::astar_bfs_visitor bfs_vis{h,      vis,   q,       predecessor, cost, distance,
                                          weight, color, combine, compare,     zero};

  breadth_first_visit(g, s, q, bfs_vis, color);
}

template <concepts::IncidenceGraph G, concepts::AStarHeuristic<G> H, concepts::AStarVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadWriteVertexPropertyMap<G> CostMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap, concepts::ReadableEdgePropertyMap<G> WeightMap,
          typename CompareFunction, typename CombineFunction>
void astar_search_no_init_tree(const G& g, graph_vertex_descriptor_t<G> s, H h, V vis, PredecessorMap predecessor,
                               CostMap cost, DistanceMap distance, WeightMap weight, CompareFunction compare,
                               CombineFunction combine, property_traits_value_t<CostMap> zero) {
  using Vertex = graph_vertex_descriptor_t<G>;
  using Distance = property_traits_value_t<DistanceMap>;
  using DVPair = std::pair<Distance, Vertex>;

  auto q =
      make_d_ary_heap_indirect<Vertex, 4>(make_function_property_map<DVPair>([](const auto& p) { return p.first; }),
                                          null_property_map<DVPair, std::size_t>(), compare);

  vis.discover_vertex(s, g);
  q.emplace(get(cost, s), s);
  while (!q.empty()) {
    auto [v_rank, v] = q.top();
    q.pop();
    vis.examine_vertex(v, g);
    for (auto e : out_edges(v, g)) {
      Vertex w = target(e, g);
      vis.examine_edge(e, g);
      Distance e_weight = get(weight, e);
      if (compare(e_weight, zero)) {
        throw negative_edge();
      }
      bool decreased = relax(e, g, weight, predecessor, distance, combine, compare);
      if (decreased) {
        vis.edge_relaxed(e, g);
        Distance w_rank = combine(get(distance, w), h(w));
        put(cost, w, w_rank);
        vis.discover_vertex(w, g);
        q.emplace(std::move(w_rank), w);
      } else {
        vis.edge_not_relaxed(e, g);
      }
    }
    vis.finish_vertex(v, g);
  }
}

// Non-named parameter interface
template <concepts::VertexListGraph G, concepts::AStarHeuristic<G> H, concepts::AStarVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadWriteVertexPropertyMap<G> CostMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, concepts::ReadableVertexIndexMap<G> VertexIndexMap,
          typename CompareFunction, typename CombineFunction>
void astar_search(const G& g, graph_vertex_descriptor_t<G> s, H h, V vis, PredecessorMap predecessor, CostMap cost,
                  DistanceMap distance, WeightMap weight, VertexIndexMap index_map, ColorMap color,
                  CompareFunction compare, CombineFunction combine, property_traits_value_t<CostMap> inf,
                  property_traits_value_t<CostMap> zero) {
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;
  for (auto u : vertices(g)) {
    put(color, u, Color::white());
    put(distance, u, inf);
    put(cost, u, inf);
    put(predecessor, u, u);
    vis.initialize_vertex(u, g);
  }
  put(distance, s, zero);
  put(cost, s, h(s));

  astar_search_no_init(g, s, h, vis, predecessor, cost, distance, weight, color, index_map, compare, combine, zero);
}

// Version with default cost arithmetic (float, double, int..)
template <concepts::VertexListGraph G, concepts::AStarHeuristic<G> H, concepts::AStarVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadWriteVertexPropertyMap<G> CostMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, concepts::ReadableVertexIndexMap<G> VertexIndexMap>
void astar_search(const G& g, graph_vertex_descriptor_t<G> s, H h, V vis, PredecessorMap predecessor, CostMap cost,
                  DistanceMap distance, WeightMap weight, VertexIndexMap index_map, ColorMap color) {
  using CostValue = property_traits_value_t<WeightMap>;
  astar_search(g, s, h, vis, predecessor, cost, distance, weight, index_map, color, std::less<>(), closed_plus<>(),
               std::numeric_limits<CostValue>::max(), CostValue{});
}

// Non-named parameter interface
template <concepts::VertexListGraph G, concepts::AStarHeuristic<G> H, concepts::AStarVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadWriteVertexPropertyMap<G> CostMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap, concepts::ReadableEdgePropertyMap<G> WeightMap,
          typename CompareFunction, typename CombineFunction>
void astar_search_tree(const G& g, graph_vertex_descriptor_t<G> s, H h, V vis, PredecessorMap predecessor, CostMap cost,
                       DistanceMap distance, WeightMap weight, CompareFunction compare, CombineFunction combine,
                       property_traits_value_t<CostMap> inf, property_traits_value_t<CostMap> zero) {
  for (auto u : vertices(g)) {
    put(distance, u, inf);
    put(cost, u, inf);
    put(predecessor, u, u);
    vis.initialize_vertex(u, g);
  }
  put(distance, s, zero);
  put(cost, s, h(s));

  astar_search_no_init_tree(g, s, h, vis, predecessor, cost, distance, weight, compare, combine, zero);
}

template <concepts::VertexListGraph G, concepts::AStarHeuristic<G> H, concepts::AStarVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::ReadWriteVertexPropertyMap<G> CostMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap, concepts::ReadableEdgePropertyMap<G> WeightMap>
void astar_search_tree(const G& g, graph_vertex_descriptor_t<G> s, H h, V vis, PredecessorMap predecessor, CostMap cost,
                       DistanceMap distance, WeightMap weight) {
  using CostValue = property_traits_value_t<WeightMap>;
  astar_search_tree(g, s, h, vis, predecessor, cost, distance, weight, std::less<>(), closed_plus<>(),
                    std::numeric_limits<CostValue>::max(), CostValue{});
}

}  // namespace bagl

#endif  // BAGL_BAGL_ASTAR_SEARCH_H_
