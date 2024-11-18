// Copyright 2007 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_IS_KURATOWSKI_SUBGRAPH_H_
#define BAGL_BAGL_IS_KURATOWSKI_SUBGRAPH_H_

#include <algorithm>
#include <map>
#include <ranges>
#include <set>
#include <tuple>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/isomorphism.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"

namespace bagl {

namespace subgraph_detail {

template <typename Graph>
Graph make_k_5() {
  Graph k_5(5);
  auto v_rg = vertices(k_5);
  for (auto vi = v_rg.begin(); vi != v_rg.end(); ++vi) {
    for (auto inner_vi = std::next(vi); inner_vi != v_rg.end(); ++inner_vi) {
      add_edge(*vi, *inner_vi, k_5);
    }
  }
  return k_5;
}

template <typename Graph>
Graph make_k_3_3() {
  Graph k_3_3(6);
  auto v_rg = vertices(k_3_3);
  auto bipartition_start = std::next(std::next(std::next(v_rg.begin())));
  for (auto vi = v_rg.begin(); vi != bipartition_start; ++vi) {
    for (auto inner_vi = bipartition_start; inner_vi != v_rg.end(); ++inner_vi) {
      add_edge(*vi, *inner_vi, k_3_3);
    }
  }
  return k_3_3;
}

template <typename AdjacencyList, typename Vertex>
void contract_edge(AdjacencyList& neighbors, Vertex u, Vertex v) {
  // Remove u from v's neighbor list
  neighbors[v].erase(std::remove(neighbors[v].begin(), neighbors[v].end(), u), neighbors[v].end());

  // Replace any references to u with references to v
  for (Vertex u_neighbor : neighbors[u]) {
    std::replace(neighbors[u_neighbor].begin(), neighbors[u_neighbor].end(), u, v);
  }

  // Remove v from u's neighbor list
  neighbors[u].erase(std::remove(neighbors[u].begin(), neighbors[u].end(), v), neighbors[u].end());

  // Add everything in u's neighbor list to v's neighbor list
  std::copy(neighbors[u].begin(), neighbors[u].end(), std::back_inserter(neighbors[v]));

  // Clear u's neighbor list
  neighbors[u].clear();
}

enum target_graph_t { tg_k_3_3, tg_k_5 };

}  // namespace subgraph_detail

template <typename Graph, std::ranges::input_range ERange, typename VertexIndexMap>
bool is_kuratowski_subgraph(const Graph& g, ERange e_rg, VertexIndexMap vm) {
  using Vertex = graph_vertex_descriptor_t<Graph>;

  using small_graph_t = adjacency_list<vec_s, vec_s, undirected_s>;

  subgraph_detail::target_graph_t target_graph = subgraph_detail::tg_k_3_3;  // unless we decide otherwise later

  static auto k_5 = subgraph_detail::make_k_5<small_graph_t>();

  static auto k_3_3 = subgraph_detail::make_k_3_3<small_graph_t>();

  std::size_t n_vertices(num_vertices(g));
  std::size_t max_num_edges(3 * n_vertices - 5);

  auto neighbors = vector_property_map(n_vertices, vm, std::vector<Vertex>{});

  std::size_t count = 0;
  for (auto e : e_rg) {
    if (count++ > max_num_edges) {
      return false;
    }

    Vertex u(source(e, g));
    Vertex v(target(e, g));

    neighbors[u].push_back(v);
    neighbors[v].push_back(u);
  }

  for (std::size_t max_size = 2; max_size < 5; ++max_size) {
    for (Vertex v : vertices(g)) {
      // a hack to make sure we don't contract the middle edge of a path
      // of four degree-3 vertices
      if (max_size == 4 && neighbors[v].size() == 3) {
        if (neighbors[neighbors[v][0]].size() + neighbors[neighbors[v][1]].size() + neighbors[neighbors[v][2]].size() <
            11) {
          // so, it has two degree-3 neighbors
          continue;
        }
      }

      while (neighbors[v].size() > 0 && neighbors[v].size() < max_size) {
        // Find one of v's neighbors u such that v and u
        // have no neighbors in common. We'll look for such a
        // neighbor with a naive cubic-time algorithm since the
        // max size of any of the neighbor sets we'll consider
        // merging is 3

        bool neighbor_sets_intersect = false;

        Vertex min_u = graph_traits<Graph>::null_vertex();
        for (Vertex u : neighbors[v]) {
          neighbor_sets_intersect = false;
          for (auto u_neighbor : neighbors[u]) {
            for (auto inner_v_neighbor : neighbors[v]) {
              if (u_neighbor == inner_v_neighbor) {
                neighbor_sets_intersect = true;
                break;
              }
            }
            if (neighbor_sets_intersect) {
              break;
            }
          }
          if (!neighbor_sets_intersect &&
              (min_u == graph_traits<Graph>::null_vertex() || neighbors[u].size() < neighbors[min_u].size())) {
            min_u = u;
          }
        }

        if (min_u == graph_traits<Graph>::null_vertex()) {
          // Exited the loop without finding an appropriate neighbor
          // of v, so v must be a lost cause. Move on to other
          // vertices.
          break;
        }

        subgraph_detail::contract_edge(neighbors, min_u, v);
      }  // end iteration over v's neighbors

    }  // end iteration through vertices v

    if (max_size == 3) {
      // check to see whether we should go on to find a k_5
      for (auto v : vertices(g)) {
        if (neighbors[v].size() == 4) {
          target_graph = subgraph_detail::tg_k_5;
          break;
        }
      }

      if (target_graph == subgraph_detail::tg_k_3_3) {
        break;
      }
    }

  }  // end iteration through max degree 2,3, and 4

  // Now, there should only be 5 or 6 vertices with any neighbors. Find them.

  std::vector<Vertex> main_vertices;
  for (auto v : vertices(g)) {
    if (!neighbors[v].empty()) {
      main_vertices.push_back(v);
    }
  }

  // create a graph isomorphic to the contracted graph to test
  // against k_5 and k_3_3
  small_graph_t contracted_graph(main_vertices.size());
  std::map<Vertex, graph_vertex_descriptor_t<small_graph_t> > contracted_vertex_map;

  auto sv_rg = vertices(contracted_graph);
  auto si = sv_rg.begin();
  for (auto itr = main_vertices.begin(); itr != main_vertices.end() && si != sv_rg.end(); ++itr, ++si) {
    contracted_vertex_map[*itr] = *si;
  }

  for (auto i : main_vertices) {
    for (auto j : neighbors[i]) {
      if (get(vm, i) < get(vm, j)) {
        add_edge(contracted_vertex_map[i], contracted_vertex_map[j], contracted_graph);
      }
    }
  }

  if (target_graph == subgraph_detail::tg_k_5) {
    auto vm_ident = identity_property_map{};
    return isomorphism(
        k_5, contracted_graph,
        vector_property_map(num_vertices(k_5), vm_ident, graph_traits<small_graph_t>::null_vertex()).ref(),
        make_degree_invariant(k_5, vm_ident), make_degree_invariant(contracted_graph, vm_ident), 0, vm_ident, vm_ident);
  } else  // target_graph == tg_k_3_3
  {
    auto vm_ident = identity_property_map{};
    return isomorphism(
        k_3_3, contracted_graph,
        vector_property_map(num_vertices(k_3_3), vm_ident, graph_traits<small_graph_t>::null_vertex()).ref(),
        make_degree_invariant(k_3_3, vm_ident), make_degree_invariant(contracted_graph, vm_ident), 0, vm_ident,
        vm_ident);
  }
}

template <typename Graph, std::ranges::input_range ERange>
bool is_kuratowski_subgraph(const Graph& g, ERange e_rg) {
  return is_kuratowski_subgraph(g, e_rg, get(vertex_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_IS_KURATOWSKI_SUBGRAPH_H_
