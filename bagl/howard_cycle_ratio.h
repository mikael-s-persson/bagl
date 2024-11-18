// Copyright (C) 2006-2009 Dmitry Bufistov and Andrey Parfenov
// Authors: Dmitry Bufistov, Andrey Parfenov
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_HOWARD_CYCLE_RATIO_H_
#define BAGL_BAGL_HOWARD_CYCLE_RATIO_H_

#include <algorithm>
#include <functional>
#include <limits>
#include <list>
#include <queue>
#include <tuple>
#include <type_traits>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"

// The implementation of the maximum/minimum cycle ratio/mean algorithm.

namespace bagl {

// The mcr_float is like numeric_limits, but only for floating point types
// and only defines infinity() and epsilon(). This class is primarily used
// to encapsulate a less-precise epsilon than natively supported by the
// floating point type.
template <typename Float = double>
struct mcr_float {
  using value_type = Float;
  static constexpr Float infinity() { return std::numeric_limits<Float>::infinity(); }
  static constexpr Float epsilon() { return Float{-0.005}; }
};

namespace cycle_ratio_detail {

struct min_comparator_props {
  using comparator = std::greater<>;
  static constexpr int multiplier = 1;
};

struct max_comparator_props {
  using comparator = std::less<>;
  static constexpr int multiplier = -1;
};

template <typename FloatTraits, typename ComparatorProps>
struct float_wrapper {
  using value_type = typename FloatTraits::value_type;
  using comparator_props_t = ComparatorProps;
  using comparator = typename ComparatorProps::comparator;

  static constexpr auto infinity() { return FloatTraits::infinity() * ComparatorProps::multiplier; }

  static constexpr auto epsilon() { return FloatTraits::epsilon() * ComparatorProps::multiplier; }
};

// Calculates optimum (maximum/minimum) cycle ratio of a directed
// graph. Uses  Howard's iteration policy algorithm. </br>(It is described
// in the paper "Experimental Analysis of the Fastest Optimum Cycle Ratio
// and Mean Algorithm" by Ali Dasdan).
template <typename FloatTraits, typename Graph, typename VertexIndexMap, typename EdgeWeight1, typename EdgeWeight2>
class mcr_howard {
 public:
  using float_t = typename FloatTraits::value_type;
  using cmp_props_t = typename FloatTraits::comparator_props_t;
  using comparator_t = typename FloatTraits::comparator;
  enum my_color_type { my_white = 0, my_black };
  using vertex_t = graph_vertex_descriptor_t<Graph>;
  using edge_t = graph_edge_descriptor_t<Graph>;
  using distance_map_t = vector_property_map<float_t, VertexIndexMap>;  // V -> float_t

  using vcol_t = std::vector<my_color_type>;
  using policy_t = vector_property_map<edge_t, VertexIndexMap>;  // Vertex -> Edge
  using color_map_t = iterator_property_map<typename vcol_t::iterator, VertexIndexMap>;

  using pinel_t = std::list<vertex_t>;  // The in_edges list of the policy graph
  using inedges_t = vector_property_map<pinel_t, VertexIndexMap>;
  using critical_cycle_t = std::vector<edge_t>;

  // Bad  vertex flag. If true, then the vertex is "bad".
  // Vertex is "bad" if its out_degree is equal to zero.
  using badv_t = vector_property_map<bool, VertexIndexMap>;

  /*!
   * Constructor
   * \param g = (V, E) - a directed multigraph.
   * \param vim  Vertex Index Map. Read property Map: V -> [0,
   * num_vertices(g)). \param ewm  edge weight map. Read property map: E
   * -> R \param ew2m  edge weight map. Read property map: E -> R+ \param
   * infty A big enough value to guaranty that there exist a cycle with
   *  better ratio.
   * \param cmp The compare operator for float_ts.
   */
  mcr_howard(const Graph& g, VertexIndexMap vim, EdgeWeight1 ewm, EdgeWeight2 ew2m, bool store_cc = false)
      : g_(g),
        v_index_(vim),
        weight1_(ewm),
        weight2_(ew2m),
        cycle_ratio_bound_(mcr_bound()),
        cycle_ratio_(cycle_ratio_bound_),
        v_count_(num_vertices(g_)),
        distance_(v_count_, v_index_, 0),
        policy_(v_count_, v_index_, edge_t{}),
        in_edge_list_(v_count_, v_index_, {}),
        bad_v_(v_count_, v_index_, false),
        color_cycle_(v_count_),
        color_bfs_(v_count_),
        store_critical_cycle_(store_cc) {}

  /*!
   * \return maximum/minimum_{for all cycles C}
   *         [sum_{e in C} w1(e)] / [sum_{e in C} w2(e)],
   * or FloatTraits::infinity() if graph has no cycles.
   */
  float_t ocr_howard() {
    construct_policy_graph();
    int k = 0;
    float_t mcr = 0;
    do {
      mcr = policy_mcr();
      ++k;
    } while (try_improve_policy(mcr) && k < 100);  // To avoid infinite loop

    constexpr float_t tolerance = -0.00000001 * cmp_props_t::multiplier;
    if (compare_(mcr, cycle_ratio_bound_ + tolerance)) {
      return FloatTraits::infinity();
    }
    return mcr;
  }

  void get_critical_cycle(critical_cycle_t& cc) {
    if (store_critical_cycle_) {
      cc.swap(critical_cycle_);
    }
  }

 private:
  void store_critical_edge(edge_t ed, critical_cycle_t& cc) {
    if (store_critical_cycle_) {
      cc.push_back(ed);
    }
  }

  void store_critical_cycle(critical_cycle_t& cc) {
    if (store_critical_cycle_) {
      critical_cycle_.swap(cc);
    }
  }

  /*!
   * \return lower/upper bound for the maximal/minimal cycle ratio
   */
  float_t mcr_bound() {
    float_t cz = std::numeric_limits<float_t>::max();  // Closest to zero value
    float_t s = 0;
    constexpr float_t tolerance = std::numeric_limits<float_t>::epsilon();
    for (auto v : vertices(g_)) {
      for (auto e : out_edges(v, g_)) {
        s += std::abs(weight1_[e]);
        float_t a = std::abs(weight2_[e]);
        if (a > tolerance && a < cz) {
          cz = a;
        }
      }
    }
    return cmp_props_t::multiplier * (s / cz);
  }

  /*!
   *  Constructs an arbitrary policy graph.
   */
  void construct_policy_graph() {
    sink_ = graph_traits<Graph>().null_vertex();
    for (auto v : vertices(g_)) {
      auto eo_rg = out_edges(v, g_);
      if (eo_rg.empty()) {
        if (sink_ == graph_traits<Graph>().null_vertex()) {
          sink_ = v;
        }
        bad_v_[v] = true;
        in_edge_list_[sink_].push_back(v);
      } else {
        auto ei = eo_rg.begin();
        auto e_max = *ei;
        for (++ei; ei != eo_rg.end(); ++ei) {
          if (compare_(weight1_[e_max], weight1_[*ei])) {
            e_max = *ei;
          }
        }
        in_edge_list_[target(e_max, g_)].push_back(v);
        policy_[v] = e_max;
      }
    }
  }
  /*! Sets the distance value for all vertices "v" such that there is
   * a path from "v" to "sv". It does "inverse" breadth first visit of the
   * policy graph, starting from the vertex "sv".
   */
  void mcr_bfv(vertex_t sv, float_t cr) {
    auto vcm = iterator_property_map(color_bfs_.begin(), v_index_);
    buffer_queue<vertex_t> q;
    vcm[sv] = my_black;
    q.push(sv);
    while (!q.empty()) {
      vertex_t v = q.top();
      q.pop();
      for (auto u : in_edge_list_[v])
      // For all in_edges of the policy graph
      {
        if (u != sv) {
          if (bad_v_[u]) {
            distance_[u] = distance_[v] + cycle_ratio_bound_ - cr;
          } else {
            distance_[u] = distance_[v] + weight1_[policy_[u]] - weight2_[policy_[u]] * cr;
          }
          vcm[u] = my_black;
          q.push(u);
        }
      }
    }
  }

  /*!
   * \param sv an arbitrary (undiscovered) vertex of the policy graph.
   * \return a vertex in the policy graph that belongs to a cycle.
   * Performs a depth first visit until a cycle edge is found.
   */
  vertex_t find_cycle_vertex(vertex_t sv) {
    vertex_t gv = sv;
    std::fill(color_cycle_.begin(), color_cycle_.end(), my_white);
    auto cm = iterator_property_map(color_cycle_.begin(), v_index_);
    do {
      cm[gv] = my_black;
      if (!bad_v_[gv]) {
        gv = target(policy_[gv], g_);
      } else {
        gv = sink_;
      }
    } while (cm[gv] != my_black);
    return gv;
  }

  /*!
   * \param sv - vertex that belongs to a cycle in the policy graph.
   */
  float_t cycle_ratio(vertex_t sv) {
    if (sv == sink_) {
      return cycle_ratio_bound_;
    }
    float_t sums_num{0};
    float_t sums_den{0};
    vertex_t v = sv;
    critical_cycle_t cc;
    do {
      store_critical_edge(policy_[v], cc);
      sums_num += weight1_[policy_[v]];
      sums_den += weight2_[policy_[v]];
      v = target(policy_[v], g_);
    } while (v != sv);
    float_t cr = sums_num / sums_den;
    if (compare_(cycle_ratio_, cr)) {
      cycle_ratio_ = cr;
      store_critical_cycle(cc);
    }
    return cr;
  }

  /*!
   *  Finds the optimal cycle ratio of the policy graph
   */
  float_t policy_mcr() {
    std::fill(color_bfs_.begin(), color_bfs_.end(), my_white);
    auto vcm = iterator_property_map(color_bfs_.begin(), v_index_);
    float_t mcr = cycle_ratio_bound_;
    // For undiscovered vertices
    auto is_undiscovered = [&](vertex_t u) { return my_white == vcm[u]; };
    for (auto v : vertices(g_) | std::views::filter(is_undiscovered)) {
      vertex_t gv = find_cycle_vertex(v);
      float_t cr = cycle_ratio(gv);
      mcr_bfv(gv, cr);
      if (compare_(mcr, cr)) {
        mcr = cr;
      }
    }
    return mcr;
  }

  /*!
   * Changes the edge m_policy[s] to the new_edge.
   */
  void improve_policy(vertex_t s, edge_t new_edge) {
    vertex_t t = target(policy_[s], g_);
    in_edge_list_[t].erase(std::find(in_edge_list_[t].begin(), in_edge_list_[t].end(), s));
    policy_[s] = new_edge;
    t = target(new_edge, g_);
    in_edge_list_[t].push_back(s);  /// Maintain in_edge list
  }

  /*!
   * A negative cycle detector.
   */
  bool try_improve_policy(float_t cr) {
    bool improved = false;
    constexpr auto tolerance = FloatTraits::epsilon();
    for (auto v : vertices(g_)) {
      if (!bad_v_[v]) {
        for (auto e : out_edges(v, g_)) {
          vertex_t t = target(e, g_);
          // Current distance from v to some vertex
          float_t dis = weight1_[e] - weight2_[e] * cr + distance_[t];
          if (compare_(distance_[v] + tolerance, dis)) {
            improve_policy(v, e);
            distance_[v] = dis;
            improved = true;
          }
        }
      } else {
        float_t dis = cycle_ratio_bound_ - cr + distance_[sink_];
        if (compare_(distance_[v] + tolerance, dis)) {
          distance_[v] = dis;
        }
      }
    }
    return improved;
  }

  const Graph& g_;
  VertexIndexMap v_index_;
  EdgeWeight1 weight1_;
  EdgeWeight2 weight2_;
  comparator_t compare_;
  float_t cycle_ratio_bound_;  //> The lower/upper bound to the maximal/minimal cycle ratio
  float_t cycle_ratio_;        //>The best cycle ratio that has been found so far

  std::size_t v_count_;      //>The number of the vertices in the graph
  distance_map_t distance_;  //>Distance map

  policy_t policy_;  //>The interface for the policy graph

  inedges_t in_edge_list_;  //>Policy graph, input edges list

  badv_t bad_v_;  // Marks "bad" vertices

  vcol_t color_cycle_;
  vcol_t color_bfs_;  // Color maps
  vertex_t sink_;     // To convert any graph to "good"

  bool store_critical_cycle_ = false;
  critical_cycle_t critical_cycle_;  // Critical cycle
};

/*!
 * \param g a directed multigraph.
 * \param vim Vertex Index Map. A map V->[0, num_vertices(g))
 * \param ewm Edge weight1 map.
 * \param ew2m Edge weight2 map.
 * \param pcc  pointer to the critical edges list.
 * \return Optimum cycle ratio of g or FloatTraits::infinity() if g has no
 * cycles.
 */
template <typename FT, concepts::VertexListGraph TG, concepts::ReadableVertexPropertyMap<TG> TVIM,
          concepts::ReadableEdgePropertyMap<TG> TEW1, concepts::ReadableEdgePropertyMap<TG> TEW2, typename EV>
requires concepts::IncidenceGraph<TG>
auto optimum_cycle_ratio(const TG& g, TVIM vim, TEW1 ewm, TEW2 ew2m, EV* pcc) {
  static_assert(is_directed_graph_v<TG>);

  mcr_howard<FT, TG, TVIM, TEW1, TEW2> obj(g, vim, ewm, ew2m, (pcc != nullptr));
  const double ocr = obj.ocr_howard();
  if (pcc != nullptr) {
    obj.get_critical_cycle(*pcc);
  }
  return ocr;
}
}  // namespace cycle_ratio_detail

// Algorithms
// Maximum Cycle Ratio

template <typename FloatTraits, typename Graph, typename VertexIndexMap, typename EdgeWeight1Map,
          typename EdgeWeight2Map>
auto maximum_cycle_ratio(const Graph& g, VertexIndexMap vim, EdgeWeight1Map ew1m, EdgeWeight2Map ew2m,
                         std::vector<graph_edge_descriptor_t<Graph>>* pcc = nullptr,
                         FloatTraits /*unused*/ = FloatTraits()) {
  using Traits = cycle_ratio_detail::float_wrapper<FloatTraits, cycle_ratio_detail::max_comparator_props>;
  return cycle_ratio_detail::optimum_cycle_ratio<Traits>(g, vim, ew1m, ew2m, pcc);
}

template <typename Graph, typename VertexIndexMap, typename EdgeWeight1Map, typename EdgeWeight2Map>
double maximum_cycle_ratio(const Graph& g, VertexIndexMap vim, EdgeWeight1Map ew1m, EdgeWeight2Map ew2m,
                           std::vector<graph_edge_descriptor_t<Graph>>* pcc = nullptr) {
  return maximum_cycle_ratio(g, vim, ew1m, ew2m, pcc, mcr_float<>());
}

// Minimum Cycle Ratio

template <typename FloatTraits, typename Graph, typename VertexIndexMap, typename EdgeWeight1Map,
          typename EdgeWeight2Map>
auto minimum_cycle_ratio(const Graph& g, VertexIndexMap vim, EdgeWeight1Map ew1m, EdgeWeight2Map ew2m,
                         std::vector<graph_edge_descriptor_t<Graph>>* pcc = nullptr,
                         FloatTraits /*unused*/ = FloatTraits()) {
  using Traits = cycle_ratio_detail::float_wrapper<FloatTraits, cycle_ratio_detail::min_comparator_props>;
  return cycle_ratio_detail::optimum_cycle_ratio<Traits>(g, vim, ew1m, ew2m, pcc);
}

template <typename Graph, typename VertexIndexMap, typename EdgeWeight1Map, typename EdgeWeight2Map>
double minimum_cycle_ratio(const Graph& g, VertexIndexMap vim, EdgeWeight1Map ew1m, EdgeWeight2Map ew2m,
                           std::vector<graph_edge_descriptor_t<Graph>>* pcc = nullptr) {
  return minimum_cycle_ratio(g, vim, ew1m, ew2m, pcc, mcr_float<>());
}

// Maximum Cycle Mean

template <typename FloatTraits, typename Graph, typename VertexIndexMap, typename EdgeWeightMap>
auto maximum_cycle_mean(const Graph& g, VertexIndexMap vim, EdgeWeightMap ewm,
                        std::vector<graph_edge_descriptor_t<Graph>>* pcc = nullptr, FloatTraits ft = FloatTraits()) {
  using Weight = std::decay_t<property_traits_value_t<EdgeWeightMap>>;
  return maximum_cycle_ratio(g, vim, ewm, single_property_map(Weight{1}), pcc, ft);
}

template <typename Graph, typename VertexIndexMap, typename EdgeWeightMap>
double maximum_cycle_mean(const Graph& g, VertexIndexMap vim, EdgeWeightMap ewm,
                          std::vector<graph_edge_descriptor_t<Graph>>* pcc = nullptr) {
  return maximum_cycle_mean(g, vim, ewm, pcc, mcr_float<>());
}

// Minimum Cycle Mean

template <typename FloatTraits, typename Graph, typename VertexIndexMap, typename EdgeWeightMap>
auto minimum_cycle_mean(const Graph& g, VertexIndexMap vim, EdgeWeightMap ewm,
                        std::vector<graph_edge_descriptor_t<Graph>>* pcc = nullptr, FloatTraits ft = FloatTraits()) {
  using Weight = std::decay_t<property_traits_value_t<EdgeWeightMap>>;
  return minimum_cycle_ratio(g, vim, ewm, single_property_map(Weight{1}), pcc, ft);
}

template <typename Graph, typename VertexIndexMap, typename EdgeWeightMap>
double minimum_cycle_mean(const Graph& g, VertexIndexMap vim, EdgeWeightMap ewm,
                          std::vector<graph_edge_descriptor_t<Graph>>* pcc = nullptr) {
  return minimum_cycle_mean(g, vim, ewm, pcc, mcr_float<>());
}

}  // namespace bagl

#endif  // BAGL_BAGL_HOWARD_CYCLE_RATIO_H_
