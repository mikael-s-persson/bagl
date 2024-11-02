// Copyright (C) 2001 Vladimir Prus <ghost@cs.msu.su>
// Copyright (C) 2001 Jeremy Siek <jsiek@cs.indiana.edu>
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_TRANSITIVE_CLOSURE_H_
#define BAGL_BAGL_TRANSITIVE_CLOSURE_H_

#include <algorithm>  // for std::min and std::max
#include <functional>
#include <ranges>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/graph_concepts.h"
#include "bagl/strong_components.h"
#include "bagl/topological_sort.h"
#include "bagl/vector_property_map.h"

namespace bagl {

template <concepts::VertexListGraph Graph, concepts::MutableGraph GraphTC, typename G_to_TC_VertexMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
requires concepts::AdjacencyGraph<Graph>
void transitive_closure(const Graph& g, GraphTC& tc, G_to_TC_VertexMap g_to_tc_map, VertexIndexMap index_map) {
  if (vertices(g).empty()) {
    return;
  }
  using vertex = graph_vertex_descriptor_t<Graph>;
  auto component_number = vector_property_map(num_vertices(g), index_map, std::size_t{0});

  const std::size_t num_scc = strong_components(
      g, component_number.ref(),
      vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<Graph>::null_vertex()).ref(),
      vector_property_map(num_vertices(g), get(vertex_index, g), std::size_t{0}).ref(), index_map);

  std::vector<std::vector<vertex> > components;
  build_component_lists(g, num_scc, component_number.ref(), components);

  adjacency_list<vec_s, vec_s, directed_s> cg(num_scc);
  for (std::size_t s = 0; s < components.size(); ++s) {
    std::vector<std::size_t> adj;
    for (std::size_t i = 0; i < components[s].size(); ++i) {
      vertex u = components[s][i];
      for (auto v : adjacent_vertices(u, g)) {
        std::size_t t = component_number[v];
        if (s != t) {
          // Avoid loops in the condensation graph
          adj.push_back(t);
        }
      }
    }
    std::sort(adj.begin(), adj.end());
    const auto di = std::unique(adj.begin(), adj.end());

    for (auto i = adj.begin(); i != di; ++i) {
      add_edge(s, *i, cg);
    }
  }

  std::vector<std::size_t> topo_order;
  std::vector<std::size_t> topo_number(num_vertices(cg));
  topological_sort(cg, std::back_inserter(topo_order));
  std::reverse(topo_order.begin(), topo_order.end());
  std::size_t n = 0;
  for (auto i : topo_order) {
    topo_number[i] = n++;
  }

  std::vector<std::vector<std::size_t> > cg_vec(num_vertices(cg));
  for (std::size_t i = 0; i < num_vertices(cg); ++i) {
    auto pr = adjacent_vertices(i, cg);
    cg_vec[i].assign(pr.begin(), pr.end());
    std::sort(cg_vec[i].begin(), cg_vec[i].end(),
              [&topo_number](std::size_t r, std::size_t s) { return topo_number[r] < topo_number[s]; });
  }

  std::vector<std::vector<std::size_t> > chains;
  {
    std::vector<std::size_t> in_a_chain(cg_vec.size());
    for (std::size_t v : topo_order) {
      if (!in_a_chain[v]) {
        chains.resize(chains.size() + 1);
        auto& chain = chains.back();
        for (;;) {
          chain.push_back(v);
          in_a_chain[v] = true;

          auto next =
              std::find_if(cg_vec[v].begin(), cg_vec[v].end(), [&in_a_chain](std::size_t u) { return !in_a_chain[u]; });

          if (next != cg_vec[v].end()) {
            v = *next;
          } else {
            break;  // end of chain, dead-end
          }
        }
      }
    }
  }
  std::vector<std::size_t> chain_number(cg_vec.size());
  std::vector<std::size_t> pos_in_chain(cg_vec.size());
  for (std::size_t i = 0; i < chains.size(); ++i) {
    for (std::size_t j = 0; j < chains[i].size(); ++j) {
      std::size_t v = chains[i][j];
      chain_number[v] = i;
      pos_in_chain[v] = j;
    }
  }

  std::size_t inf = std::numeric_limits<std::size_t>::max();
  std::vector<std::vector<std::size_t> > successors(cg_vec.size(), std::vector<std::size_t>(chains.size(), inf));
  for (std::size_t u : topo_order | std::views::reverse) {
    for (std::size_t v : cg_vec[u]) {
      if (topo_number[v] < successors[u][chain_number[v]]) {
        // Succ(u) = Succ(u) U Succ(v)
        for (std::size_t k = 0; k < successors[u].size(); ++k) {
          successors[u][k] = std::min(successors[u][k], successors[v][k]);
        }
        // Succ(u) = Succ(u) U {v}
        successors[u][chain_number[v]] = topo_number[v];
      }
    }
  }

  for (auto& i : cg_vec) {
    i.clear();
  }
  for (std::size_t i = 0; i < cg_vec.size(); ++i) {
    for (std::size_t j = 0; j < chains.size(); ++j) {
      std::size_t topo_num = successors[i][j];
      if (topo_num < inf) {
        std::size_t v = topo_order[topo_num];
        for (std::size_t k = pos_in_chain[v]; k < chains[j].size(); ++k) {
          cg_vec[i].push_back(chains[j][k]);
        }
      }
    }
  }

  // Add vertices to the transitive closure graph
  for (auto v : vertices(g)) {
    g_to_tc_map[v] = add_vertex(tc);
  }
  // Add edges between all the vertices in two adjacent SCCs
  for (auto si = cg_vec.begin(), si_end = cg_vec.end(); si != si_end; ++si) {
    std::size_t s = si - cg_vec.begin();
    for (std::size_t t : cg_vec[s]) {
      for (std::size_t k = 0; k < components[s].size(); ++k) {
        for (std::size_t l = 0; l < components[t].size(); ++l) {
          add_edge(g_to_tc_map[components[s][k]], g_to_tc_map[components[t][l]], tc);
        }
      }
    }
  }
  // Add edges connecting all vertices in a SCC
  for (std::size_t i = 0; i < components.size(); ++i) {
    if (components[i].size() > 1) {
      for (std::size_t k = 0; k < components[i].size(); ++k) {
        for (std::size_t l = 0; l < components[i].size(); ++l) {
          vertex u = components[i][k];
          vertex v = components[i][l];
          add_edge(g_to_tc_map[u], g_to_tc_map[v], tc);
        }
      }
    }
  }

  // Find loopbacks in the original graph.
  // Need to add it to transitive closure.
  for (auto v : vertices(g)) {
    for (auto a : adjacent_vertices(v, g)) {
      if (a == v && components[component_number[v]].size() == 1) {
        add_edge(g_to_tc_map[v], g_to_tc_map[v], tc);
      }
    }
  }
}

template <concepts::VertexListGraph Graph, concepts::MutableGraph GraphTC, typename G_to_TC_VertexMap>
requires concepts::AdjacencyGraph<Graph>
void transitive_closure(const Graph& g, GraphTC& tc, G_to_TC_VertexMap g_to_tc_map) {
  transitive_closure(g, tc, g_to_tc_map, get(vertex_index, g));
}

template <typename Graph, typename GraphTC>
void transitive_closure(const Graph& g, GraphTC& tc) {
  auto index_map = get(vertex_index, g);
  auto g_to_tc_map = vector_property_map(num_vertices(g), index_map, graph_traits<GraphTC>::null_vertex());
  transitive_closure(g, tc, g_to_tc_map.ref(), index_map);
}

template <concepts::EdgeMutableGraph G>
requires concepts::AdjacencyMatrix<G> && concepts::VertexListGraph<G>
void warshall_transitive_closure(G& g) {
  // Matrix form:
  // for k
  //  for i
  //    if A[i,k]
  //      for j
  //        A[i,j] = A[i,j] | A[k,j]
  for (auto k : vertices(g)) {
    for (auto i : vertices(g)) {
      if (edge(i, k, g).second) {
        for (auto j : vertices(g)) {
          if (!edge(i, j, g).second && edge(k, j, g).second) {
            add_edge(i, j, g);
          }
        }
      }
    }
  }
}

template <concepts::EdgeMutableGraph G>
requires concepts::AdjacencyMatrix<G> && concepts::VertexListGraph<G>
void warren_transitive_closure(G& g) {
  // Make sure second loop will work
  if (vertices(g).empty()) {
    return;
  }

  // for i = 2 to n
  //    for k = 1 to i - 1
  //      if A[i,k]
  //        for j = 1 to n
  //          A[i,j] = A[i,j] | A[k,j]
  auto v_rg = vertices(g);
  for (auto ic = std::next(v_rg.begin()); ic != v_rg.end(); ++ic) {
    for (auto kc = v_rg.begin(); kc != ic; ++kc) {
      if (edge(*ic, *kc, g).second) {
        for (auto j : vertices(g)) {
          if (!edge(*ic, j, g).second && edge(*kc, j, g).second) {
            add_edge(*ic, j, g);
          }
        }
      }
    }
  }
  //  for i = 1 to n - 1
  //    for k = i + 1 to n
  //      if A[i,k]
  //        for j = 1 to n
  //          A[i,j] = A[i,j] | A[k,j]
  for (auto ic = v_rg.begin(); std::next(ic) != v_rg.end(); ++ic) {
    for (auto kc = std::next(ic); kc != v_rg.end(); ++kc) {
      if (edge(*ic, *kc, g).second) {
        for (auto j : vertices(g)) {
          if (!edge(*ic, j, g).second && edge(*kc, j, g).second) {
            add_edge(*ic, j, g);
          }
        }
      }
    }
  }
}

}  // namespace bagl

#endif  // BAGL_BAGL_TRANSITIVE_CLOSURE_H_
