// (C) Copyright 2009 Eric Bose-Wolf
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_TRANSITIVE_REDUCTION_H_
#define BAGL_BAGL_TRANSITIVE_REDUCTION_H_

#include <algorithm>  //std::find
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/topological_sort.h"
#include "bagl/vector_property_map.h"

// also I didn't got all of the concepts thin. Am I suppose to check
// for all concepts, which are needed for functions I call? (As if I
// wouldn't do that, the users would see the functions called by
// complaining about missings concepts, which would be clearly an error
// message revealing internal implementation and should therefore be avoided?)

// the pseudocode which I followed implementing this algorithmn was taken
// from the german book Algorithmische Graphentheorie by Volker Turau
// it is proposed to be of O(n + nm_red ) where n is the number
// of vertices and m_red is the number of edges in the transitive
// reduction, but I think my implementation spoiled this up at some point
// indicated below.

namespace bagl {

template <concepts::VertexListGraph Graph, concepts::MutableGraph GraphTR,
          concepts::LvalueVertexPropertyMap<Graph> OrigToTrVertexMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
requires concepts::IncidenceGraph<Graph>
void transitive_reduction(const Graph& g, GraphTR& tr, OrigToTrVertexMap g_to_tr_map, VertexIndexMap g_index_map) {
  using Vertex = graph_vertex_descriptor_t<Graph>;

  std::vector<Vertex> topo_order;
  topological_sort(g, std::back_inserter(topo_order));

  auto topo_number = vector_property_map(num_vertices(g), g_index_map, std::size_t{0});
  {
    std::size_t n = 0;
    for (auto i : topo_order | std::views::reverse) {
      topo_number[i] = n++;
    }
  }

  std::vector<std::vector<bool> > edge_in_closure(num_vertices(g), std::vector<bool>(num_vertices(g), false));
  for (auto i : topo_order | std::views::reverse) {
    g_to_tr_map[i] = add_vertex(tr);
  }

  for (auto v : topo_order) {
    std::size_t i = topo_number[v];
    edge_in_closure[i][i] = true;
    std::vector<Vertex> neighbors;

    // I have to collect the successors of v and traverse them in
    // ascending topological order. I didn't know a better way, how to
    // do that. So what I'm doint is, collection the successors of v here
    for (auto e : out_edges(v, g)) {
      neighbors.push_back(target(e, g));
    }

    // and run through all vertices in topological order
    for (auto u : topo_order | std::views::reverse) {
      // looking if they are successors of v
      if (std::find(neighbors.begin(), neighbors.end(), u) == neighbors.end()) {
        continue;
      }
      std::size_t j = topo_number[u];
      if (edge_in_closure[i][j]) {
        continue;
      }
      for (std::size_t k = j; k < num_vertices(g); ++k) {
        if (not edge_in_closure[i][k]) {
          // here we need edge_in_closure to be in
          // topological order,
          edge_in_closure[i][k] = edge_in_closure[j][k];
        }
      }
      // therefore we only access edge_in_closure only through
      // topo_number property_map
      add_edge(g_to_tr_map[v], g_to_tr_map[u], tr);
    }  // for reverse topo_order
  }    // for topo_order
}

}  // namespace bagl

#endif  // BAGL_BAGL_TRANSITIVE_REDUCTION_H_
