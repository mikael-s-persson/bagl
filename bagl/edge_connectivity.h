// Copyright 2000 University of Notre Dame.
// Authors: Jeremy G. Siek, Andrew Lumsdaine, Lie-Quan Lee
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_EDGE_CONNECTIVITY_H_
#define BAGL_BAGL_EDGE_CONNECTIVITY_H_

// WARNING: not-yet fully tested!

#include <algorithm>
#include <set>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/edmonds_karp_max_flow.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"

namespace bagl {

namespace edge_connectivity_detail {

template <concepts::VertexListGraph G>
auto min_degree_vertex(const G& g) {
  auto p = graph_traits<G>::null_vertex();
  std::size_t delta = std::numeric_limits<std::size_t>::max();

  for (auto v : vertices(g)) {
    if (degree(v, g) < delta) {
      delta = degree(v, g);
      p = v;
    }
  }
  return std::pair{p, delta};
}

template <concepts::AdjacencyGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator>
void neighbors(const G& g, graph_vertex_descriptor_t<G> u, OutputIterator result) {
  for (auto v : adjacent_vertices(u, g)) {
    *result++ = v;
  }
}

template <concepts::AdjacencyGraph G, std::ranges::input_range VertexRange,
          std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator>
void neighbors(const G& g, VertexRange v_rg, OutputIterator result) {
  for (auto v : v_rg) {
    neighbors(g, v, result);
  }
}

}  // namespace edge_connectivity_detail

// O(m n)
template <concepts::VertexListGraph G, class OutputIterator>
std::size_t edge_connectivity(G& g, OutputIterator disconnecting_set) {
  //-------------------------------------------------------------------------
  // Type Definitions
  using Vertex = graph_vertex_descriptor_t<G>;
  using Color = color_traits<default_color_type>;

  using Tr_edge_desc = typename adjacency_list_traits<vec_s, vec_s, directed_s>::edge_descriptor_or_void;
  using FlowGraph =
      adjacency_list<vec_s, vec_s, directed_s, no_property,
                     property<edge_capacity_t, std::size_t,
                              property<edge_residual_capacity_t, std::size_t, property<edge_reverse_t, Tr_edge_desc>>>>;
  using FlowEdge = graph_edge_descriptor_t<FlowGraph>;

  //-------------------------------------------------------------------------
  // Create a network flow graph out of the undirected graph
  FlowGraph flow_g(num_vertices(g));

  auto cap = get(edge_capacity, flow_g);
  auto res_cap = get(edge_residual_capacity, flow_g);
  auto rev_edge = get(edge_reverse, flow_g);

  for (auto e : edges(g)) {
    auto u = source(e, g);
    auto v = target(e, g);
    auto [e1, inserted1] = add_edge(u, v, flow_g);
    put(cap, e1, 1);
    auto [e2, inserted2] = add_edge(v, u, flow_g);
    put(cap, e2, 1);  // not sure about this
    put(rev_edge, e1, e2);
    put(rev_edge, e2, e1);
  }

  //-------------------------------------------------------------------------
  // The Algorithm

  auto [p, delta] = edge_connectivity_detail::min_degree_vertex(g);
  std::vector<Vertex> s_star;
  s_star.push_back(p);
  std::size_t alpha_star = delta;
  std::set<Vertex> s;
  s.insert(p);
  std::set<Vertex> neighbor_s;
  neighbor_s.insert(p);
  edge_connectivity_detail::neighbors(g, s.begin(), s.end(), std::inserter(neighbor_s, neighbor_s.begin()));

  std::vector<Vertex> non_neighbor_s;
  auto v_rg = vertices(g);
  std::set_difference(v_rg.begin(), v_rg.end(), neighbor_s.begin(), neighbor_s.end(),
                      std::back_inserter(non_neighbor_s));

  auto color = vector_property_map(num_vertices_or_zero(g), get(vertex_index, g), Color::white());
  auto pred = vector_property_map(num_vertices_or_zero(g), get(vertex_index, g), FlowEdge{});
  while (!non_neighbor_s.empty()) {  // at most n - 1 times
    auto k = non_neighbor_s.front();

    std::size_t alpha_s_k = edmonds_karp_max_flow(flow_g, p, k, cap, res_cap, rev_edge, color.ref(), pred.ref());

    if (alpha_s_k < alpha_star) {
      alpha_star = alpha_s_k;
      s_star.clear();
      for (auto u : vertices(flow_g)) {
        if (color[u] != Color::white()) {
          s_star.push_back(u);
        }
      }
    }
    s.insert(k);
    neighbor_s.insert(k);
    edge_connectivity_detail::neighbors(g, k, std::inserter(neighbor_s, neighbor_s.begin()));
    non_neighbor_s.clear();
    auto u_rg = vertices(g);
    std::set_difference(u_rg.begin(), u_rg.end(), neighbor_s.begin(), neighbor_s.end(),
                        std::back_inserter(non_neighbor_s));
  }
  //-------------------------------------------------------------------------
  // Compute edges of the cut [s*, ~s*]
  std::vector<bool> in_s_star(num_vertices(g), false);
  for (auto u : s_star) {
    in_s_star[u] = true;
  }

  std::size_t c = 0;
  for (auto u : s_star) {
    for (auto e : out_edges(u, g)) {
      if (!in_s_star[target(e, g)]) {
        *disconnecting_set++ = e;
        ++c;
      }
    }
  }
  return c;
}

}  // namespace bagl

#endif  // BAGL_BAGL_EDGE_CONNECTIVITY_H_
