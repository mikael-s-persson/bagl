// Copyright 2000 University of Notre Dame.
// Authors: Jeremy G. Siek, Andrew Lumsdaine, Lie-Quan Lee
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PUSH_RELABEL_MAX_FLOW_H_
#define BAGL_BAGL_PUSH_RELABEL_MAX_FLOW_H_

#include <algorithm>
#include <cassert>
#include <iosfwd>
#include <limits>
#include <list>
#include <queue>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/partial_range.h"
#include "bagl/properties.h"

namespace bagl {

namespace max_flow_detail {

// This implementation is based on Goldberg's
// "On Implementing Push-Relabel Method for the Maximum Flow Problem"
// by B.V. Cherkassky and A.V. Goldberg, IPCO '95, pp. 157--171
// and on the h_prf.c and hi_pr.c code written by the above authors.

// This implements the highest-label version of the push-relabel method
// with the global relabeling and gap relabeling heuristics.

// The terms "rank", "distance", "height" are synonyms in
// Goldberg's implementation, paper and in the CLR.  A "layer" is a
// group of vertices with the same distance. The vertices in each
// layer are categorized as active or inactive.  An active vertex
// has positive excess flow and its distance is less than n (it is
// not blocked).

template <class Vertex>
struct preflow_layer {
  std::list<Vertex> active_vertices;
  std::list<Vertex> inactive_vertices;
};

template <concepts::VertexListGraph Graph, concepts::ReadableEdgePropertyMap<Graph> EdgeCapacityMap,
          concepts::ReadWriteEdgePropertyMap<Graph> ResidualCapacityEdgeMap,
          concepts::ReadableEdgePropertyMap<Graph> ReverseEdgeMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
requires concepts::IncidenceGraph<Graph>
class push_relabel {
 public:
  using FlowValue = property_traits_value_t<EdgeCapacityMap>;
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;
  using PartialOutEdges = partial_view<decltype(out_edges(std::declval<Vertex>(), std::declval<Graph>()))>;

  using Layer = preflow_layer<Vertex>;
  using ListIterator = typename std::list<Vertex>::iterator;
  using LayerArray = std::vector<Layer>;

  using ColorTraits = color_traits<default_color_type>;

  //=======================================================================
  // Some helper predicates

  [[nodiscard]] bool is_admissible(Vertex u, Vertex v) { return get(distance, u) == get(distance, v) + 1; }
  [[nodiscard]] bool is_residual_edge(Edge a) { return 0 < get(residual_capacity, a); }
  [[nodiscard]] bool is_saturated(Edge a) { return get(residual_capacity, a) == 0; }

  //=======================================================================
  // Layer List Management Functions

  void add_to_active_list(Vertex u, Layer& layer) {
    layer.active_vertices.push_front(u);
    max_active = std::max(get(distance, u), max_active);
    min_active = std::min(get(distance, u), min_active);
    layer_list_ptr[u] = layer.active_vertices.begin();
  }
  void remove_from_active_list(Vertex u) { layers[get(distance, u)].active_vertices.erase(layer_list_ptr[u]); }

  void add_to_inactive_list(Vertex u, Layer& layer) {
    layer.inactive_vertices.push_front(u);
    layer_list_ptr[u] = layer.inactive_vertices.begin();
  }
  void remove_from_inactive_list(Vertex u) { layers[get(distance, u)].inactive_vertices.erase(layer_list_ptr[u]); }

  //=======================================================================
  // initialization
  push_relabel(Graph& g_, EdgeCapacityMap cap, ResidualCapacityEdgeMap res, ReverseEdgeMap reverse_edge_, Vertex src_,
               Vertex sink_, VertexIndexMap idx)
      : g(g_),
        n(num_vertices(g_)),
        capacity(cap),
        src(src_),
        sink(sink_),
        index(idx),
        excess_flow(num_vertices(g_), idx, FlowValue{}),
        current(num_vertices(g_), idx, PartialOutEdges{}),
        distance(num_vertices(g_), idx, std::size_t{0}),
        color(num_vertices(g_), idx, default_color_type::white_color),
        reverse_edge(reverse_edge_),
        residual_capacity(res),
        layers(num_vertices(g_)),
        layer_list_ptr(num_vertices(g_), idx, layers.front().inactive_vertices.end()) {
    // Don't count the reverse edges
    std::size_t m = num_edges(g) / 2;
    nm = alpha * n + m;

    // Initialize flow to zero which means initializing
    // the residual capacity to equal the capacity.
    for (auto u : vertices(g)) {
      for (auto e : out_edges(u, g)) {
        put(residual_capacity, e, get(capacity, e));
      }
    }

    for (auto u : vertices(g)) {
      put(excess_flow, u, 0);
      current[u] = partial_view(out_edges(u, g));
    }

    bool overflow_detected = false;
    FlowValue test_excess = 0;

    for (auto e : out_edges(src, g)) {
      if (target(e, g) != src) {
        test_excess += get(residual_capacity, e);
      }
    }
    if (test_excess > std::numeric_limits<FlowValue>::max()) {
      overflow_detected = true;
    }

    if (overflow_detected) {
      put(excess_flow, src, std::numeric_limits<FlowValue>::max());
    } else {
      put(excess_flow, src, 0);
      for (auto a : out_edges(src, g)) {
        Vertex tgt = target(a, g);
        if (tgt == src) {
          continue;
        }
        ++push_count;
        FlowValue delta = get(residual_capacity, a);
        put(residual_capacity, a, get(residual_capacity, a) - delta);
        Edge rev = get(reverse_edge, a);
        put(residual_capacity, rev, get(residual_capacity, rev) + delta);
        put(excess_flow, tgt, get(excess_flow, tgt) + delta);
      }
    }
    max_distance = num_vertices(g) - 1;
    max_active = 0;
    min_active = n;

    for (auto u : vertices(g)) {
      if (u == sink) {
        put(distance, u, 0);
        continue;
      }
      if (u == src && !overflow_detected) {
        put(distance, u, n);
      } else {
        put(distance, u, 1);
      }

      if (get(excess_flow, u) > 0) {
        add_to_active_list(u, layers[1]);
      } else if (get(distance, u) < n) {
        add_to_inactive_list(u, layers[1]);
      }
    }

  }  // push_relabel constructor

  //=======================================================================
  // This is a breadth-first search over the residual graph
  // (well, actually the reverse of the residual graph).
  // Would be cool to have a graph view adaptor for hiding certain
  // edges, like the saturated (non-residual) edges in this case.
  // Goldberg's implementation abused "distance" for the coloring.
  void global_distance_update() {
    ++update_count;
    for (auto u : vertices(g)) {
      put(color, u, ColorTraits::white());
      put(distance, u, n);
    }
    put(color, sink, ColorTraits::gray());
    put(distance, sink, 0);

    for (std::size_t l = 0; l <= max_distance; ++l) {
      layers[l].active_vertices.clear();
      layers[l].inactive_vertices.clear();
    }

    max_distance = max_active = 0;
    min_active = n;

    Q.push(sink);
    while (!Q.empty()) {
      Vertex u = Q.top();
      Q.pop();
      std::size_t d_v = get(distance, u) + 1;

      for (auto a : out_edges(u, g)) {
        Vertex v = target(a, g);
        if (get(color, v) != ColorTraits::white() || !is_residual_edge(get(reverse_edge, a))) {
          continue;
        }

        put(distance, v, d_v);
        put(color, v, ColorTraits::gray());
        current[v] = partial_view(out_edges(v, g));
        max_distance = std::max(d_v, max_distance);

        if (get(excess_flow, v) > 0) {
          add_to_active_list(v, layers[d_v]);
        } else {
          add_to_inactive_list(v, layers[d_v]);
        }

        Q.push(v);
      }
    }
  }  // global_distance_update()

  //=======================================================================
  // This function is called "push" in Goldberg's h_prf implementation,
  // but it is called "discharge" in the paper and in hi_pr.c.
  void discharge(Vertex u) {
    assert(get(excess_flow, u) > 0);
    while (true) {
      auto& current_u = current[u];
      auto ai = current_u.begin();
      for (; ai != current_u.end(); ++ai) {
        Edge a = *ai;
        if (!is_residual_edge(a)) {
          continue;
        }
        Vertex v = target(a, g);
        if (!is_admissible(u, v)) {
          continue;
        }
        ++push_count;
        if (v != sink && get(excess_flow, v) == 0) {
          remove_from_inactive_list(v);
          add_to_active_list(v, layers[get(distance, v)]);
        }
        push_flow(a);
        if (get(excess_flow, u) == 0) {
          break;
        }
      }

      Layer& layer = layers[get(distance, u)];
      std::size_t du = get(distance, u);

      if (ai != current_u.end()) {
        // i is no longer active
        current_u.move_begin_to(ai);
        add_to_inactive_list(u, layer);
        break;
      }
      // i must be relabeled
      relabel_distance(u);
      if (layer.active_vertices.empty() && layer.inactive_vertices.empty()) {
        gap(du);
      }
      if (get(distance, u) == n) {
        break;
      }
    }
  }  // discharge()

  //=======================================================================
  // This corresponds to the "push" update operation of the paper,
  // not the "push" function in Goldberg's h_prf.c implementation.
  // The idea is to push the excess flow from from vertex u to v.
  void push_flow(Edge u_v) {
    Vertex u = source(u_v, g);
    Vertex v = target(u_v, g);

    FlowValue flow_delta = std::min(get(excess_flow, u), get(residual_capacity, u_v));

    put(residual_capacity, u_v, get(residual_capacity, u_v) - flow_delta);
    Edge rev = get(reverse_edge, u_v);
    put(residual_capacity, rev, get(residual_capacity, rev) + flow_delta);

    put(excess_flow, u, get(excess_flow, u) - flow_delta);
    put(excess_flow, v, get(excess_flow, v) + flow_delta);
  }  // push_flow()

  //=======================================================================
  // The main purpose of this routine is to set distance[v]
  // to the smallest value allowed by the valid labeling constraints,
  // which are:
  // distance[t] = 0
  // distance[u] <= distance[v] + 1   for every residual edge (u,v)
  //
  std::size_t relabel_distance(Vertex u) {
    ++relabel_count;
    work_since_last_update += beta;

    std::size_t min_distance = num_vertices(g);
    put(distance, u, min_distance);

    // Examine the residual out-edges of vertex i, choosing the
    // edge whose target vertex has the minimal distance.
    auto& current_u = current[u];
    auto min_edge_iter = current_u.begin();
    for (auto ai = current_u.base_begin(); ai != current_u.end(); ++ai) {
      ++work_since_last_update;
      Edge a = *ai;
      Vertex v = target(a, g);
      if (is_residual_edge(a) && get(distance, v) < min_distance) {
        min_distance = get(distance, v);
        min_edge_iter = ai;
      }
    }
    ++min_distance;
    if (min_distance < n) {
      put(distance, u, min_distance);  // this is the main action
      current_u.move_begin_to(min_edge_iter);
      max_distance = std::max(min_distance, max_distance);
    }
    return min_distance;
  }  // relabel_distance()

  //=======================================================================
  // cleanup beyond the gap
  void gap(std::size_t empty_distance) {
    ++gap_count;

    // distance of layer before the current layer
    std::size_t r = empty_distance - 1;

    // Set the distance for the vertices beyond the gap to "infinity".
    for (auto l = std::next(layers.begin(), empty_distance + 1), l_end = std::next(layers.begin(), max_distance);
         l != l_end; ++l) {
      for (auto i : l->inactive_vertices) {
        put(distance, i, n);
        ++gap_node_count;
      }
      l->inactive_vertices.clear();
    }
    max_distance = r;
    max_active = r;
  }

  //=======================================================================
  // This is the core part of the algorithm, "phase one".
  FlowValue maximum_preflow() {
    work_since_last_update = 0;

    while (max_active >= min_active) {
      Layer& layer = layers[max_active];
      auto u_iter = layer.active_vertices.begin();

      if (u_iter == layer.active_vertices.end()) {
        --max_active;
      } else {
        Vertex u = *u_iter;
        remove_from_active_list(u);

        discharge(u);

        if (work_since_last_update * global_update_frequency > nm) {
          global_distance_update();
          work_since_last_update = 0;
        }
      }
    }

    return get(excess_flow, sink);
  }  // maximum_preflow()

  //=======================================================================
  // remove excess flow, the "second phase"
  // This does a DFS on the reverse flow graph of nodes with excess flow.
  // If a cycle is found, cancel it.
  // Return the nodes with excess flow in topological order.
  //
  // Unlike the prefl_to_flow() implementation, we use
  //   "color" instead of "distance" for the DFS labels
  //   "parent" instead of nl_prev for the DFS tree
  //   "topo_next" instead of nl_next for the topological ordering
  void convert_preflow_to_flow() {
    std::vector<Vertex> parent(n);
    std::vector<Vertex> topo_next(n);

    Vertex tos = graph_traits<Graph>::null_vertex();
    Vertex bos = graph_traits<Graph>::null_vertex();
    bool bos_null = true;

    // handle self-loops
    for (auto u : vertices(g)) {
      for (auto a : out_edges(u, g)) {
        if (target(a, g) == u) {
          put(residual_capacity, a, get(capacity, a));
        }
      }
    }

    // initialize
    for (auto u : vertices(g)) {
      put(color, u, ColorTraits::white());
      parent[get(index, u)] = u;
      current[u] = partial_view(out_edges(u, g));
    }
    // eliminate flow cycles and topologically order the vertices
    for (auto u : vertices(g)) {
      if (u == src || u == sink || get(color, u) != ColorTraits::white() || get(excess_flow, u) <= 0) {
        continue;
      }
      auto r = u;
      put(color, r, ColorTraits::gray());
      while (true) {
        for (; !current[u].empty(); current[u].move_to_next()) {
          Edge a = *current[u].begin();
          if (get(capacity, a) != 0 || !is_residual_edge(a)) {
            continue;
          }
          Vertex v = target(a, g);
          if (get(color, v) == ColorTraits::white()) {
            put(color, v, ColorTraits::gray());
            parent[get(index, v)] = u;
            u = v;
            break;
          }

          if (get(color, v) != ColorTraits::gray()) {
            continue;
          }

          // find minimum flow on the cycle
          FlowValue delta = get(residual_capacity, a);
          while (true) {
            Edge b = *current[v].begin();
            delta = std::min(delta, get(residual_capacity, b));
            if (v == u) {
              break;
            }
            v = target(b, g);
          }
          // remove delta flow units
          v = u;
          while (true) {
            a = *current[v].begin();
            put(residual_capacity, a, get(residual_capacity, a) - delta);
            Edge rev = get(reverse_edge, a);
            put(residual_capacity, rev, get(residual_capacity, rev) + delta);
            v = target(a, g);
            if (v == u) {
              break;
            }
          }

          // back-out of DFS to the first saturated
          // edge
          auto restart = u;
          for (v = target(*current[u].begin(), g); v != u; v = target(a, g)) {
            a = *current[v].begin();
            if (get(color, v) != ColorTraits::white() && !is_saturated(a)) {
              continue;
            }
            put(color, target(a, g), ColorTraits::white());
            if (get(color, v) != ColorTraits::white()) {
              restart = v;
            }
          }
          if (restart != u) {
            u = restart;
            current[u].move_to_nex();
            break;
          }
        }

        if (!current[u].empty()) {
          continue;
        }
        // scan of i is complete
        put(color, u, ColorTraits::black());
        if (u != src) {
          if (bos_null) {
            bos = u;
            bos_null = false;
            tos = u;
          } else {
            topo_next[get(index, u)] = tos;
            tos = u;
          }
        }
        if (u != r) {
          u = parent[get(index, u)];
          current[u].move_to_nex();
        } else {
          break;
        }
      }
    }

    if (bos_null) {
      return;
    }
    // return excess flows
    // note that the sink is not on the stack but must be pushed too.
    auto u = tos;
    while (true) {
      auto a_rg = out_edges(u, g);
      for (auto ai = a_rg.begin(); get(excess_flow, u) > 0 && ai != a_rg.end(); ++ai) {
        if (get(capacity, *ai) == 0 && is_residual_edge(*ai)) {
          push_flow(*ai);
        }
      }
      if (u == bos) {
        break;
      }
      u = topo_next[get(index, u)];
    }
  }  // convert_preflow_to_flow()

  //=======================================================================
  bool is_flow() {
    // check edge flow values
    for (auto u : vertices(g)) {
      for (auto a : out_edges(u, g)) {
        if (get(capacity, a) <= 0) {
          continue;
        }
        auto res_a = get(residual_capacity, a);
        auto rev_a = get(reverse_edge, a);
        auto res_rev_a = get(residual_capacity, rev_a);
        if ((res_a < 0) || (res_rev_a < 0) || (res_a + res_rev_a != get(capacity, a) + get(capacity, rev_a))) {
          return false;
        }
      }
    }

    // check conservation
    for (auto u : vertices(g)) {
      if (u != src && u != sink) {
        if (get(excess_flow, u) != 0) {
          return false;
        }
        FlowValue sum = 0;
        for (auto a : out_edges(u, g)) {
          if (get(capacity, a) > 0) {
            sum -= get(capacity, a) - get(residual_capacity, a);
          } else {
            sum += get(residual_capacity, a);
          }
        }
        if (get(excess_flow, u) != sum) {
          return false;
        }
      }
    }

    return true;
  }

  bool is_optimal() {
    // check if mincut is saturated...
    global_distance_update();
    return get(distance, src) >= n;
  }

  void print_statistics(std::ostream& os) const {
    os << "pushes:     " << push_count << std::endl
       << "relabels:   " << relabel_count << std::endl
       << "updates:    " << update_count << std::endl
       << "gaps:       " << gap_count << std::endl
       << "gap nodes:  " << gap_node_count << std::endl
       << std::endl;
  }

  void print_flow_values(std::ostream& os) const {
    os << "flow values" << std::endl;
    for (auto u : vertices(g)) {
      for (auto e : out_edges(u, g)) {
        if (get(capacity, e) > 0) {
          os << u << " " << target(e, g) << " " << (get(capacity, e) - get(residual_capacity, e)) << std::endl;
        }
      }
    }
    os << std::endl;
  }

  //=======================================================================

  Graph& g;
  std::size_t n;
  std::size_t nm;
  EdgeCapacityMap capacity;
  Vertex src;
  Vertex sink;
  VertexIndexMap index;

  // will need to use random_access_property_map with these
  vector_property_map<FlowValue, VertexIndexMap> excess_flow;
  vector_property_map<PartialOutEdges, VertexIndexMap> current;
  vector_property_map<std::size_t, VertexIndexMap> distance;
  vector_property_map<default_color_type, VertexIndexMap> color;

  // Edge Property Maps that must be interior to the graph
  ReverseEdgeMap reverse_edge;
  ResidualCapacityEdgeMap residual_capacity;

  LayerArray layers;
  vector_property_map<ListIterator, VertexIndexMap> layer_list_ptr;
  std::size_t max_distance;  // maximal distance
  std::size_t max_active;    // maximal distance with active node
  std::size_t min_active;    // minimal distance with active node
  std::queue<Vertex> Q;

  // Statistics counters
  std::size_t push_count = 0;
  std::size_t update_count = 0;
  std::size_t relabel_count = 0;
  std::size_t gap_count = 0;
  std::size_t gap_node_count = 0;

  static constexpr double global_update_frequency = 0.5;
  static constexpr std::size_t alpha = 6;
  static constexpr std::size_t beta = 12;

  std::size_t work_since_last_update = 0;
};

}  // namespace max_flow_detail

template <concepts::VertexListGraph Graph, concepts::ReadableEdgePropertyMap<Graph> CapacityEdgeMap,
          concepts::ReadWriteEdgePropertyMap<Graph> ResidualCapacityEdgeMap,
          concepts::ReadableEdgePropertyMap<Graph> ReverseEdgeMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
requires concepts::IncidenceGraph<Graph>
auto push_relabel_max_flow(Graph& g, graph_vertex_descriptor_t<Graph> src, graph_vertex_descriptor_t<Graph> sink,
                           CapacityEdgeMap cap, ResidualCapacityEdgeMap res, ReverseEdgeMap rev,
                           VertexIndexMap index_map) {
  max_flow_detail::push_relabel algo{g, cap, res, rev, src, sink, index_map};
  auto flow = algo.maximum_preflow();

  algo.convert_preflow_to_flow();

  assert(algo.is_flow());
  assert(algo.is_optimal());

  return flow;
}

template <concepts::VertexListGraph Graph>
requires concepts::IncidenceGraph<Graph>
auto push_relabel_max_flow(Graph& g, graph_vertex_descriptor_t<Graph> src, graph_vertex_descriptor_t<Graph> sink) {
  return push_relabel_max_flow(g, src, sink, get(edge_capacity, g), get(edge_residual_capacity, g),
                               get(edge_reverse, g), get(vertex_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_PUSH_RELABEL_MAX_FLOW_H_
