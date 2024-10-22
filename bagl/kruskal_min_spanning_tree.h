// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_KRUSKAL_MIN_SPANNING_TREE_H_
#define BAGL_BAGL_KRUSKAL_MIN_SPANNING_TREE_H_

#include <functional>
#include <queue>
#include <vector>

#include "bagl/disjoint_sets.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/indirect_cmp.h"
#include "bagl/property_map.h"

namespace bagl {

// Kruskal's algorithm for Minimum Spanning Tree
//
// This is a greedy algorithm to calculate the Minimum Spanning Tree
// for an undirected graph with weighted edges. The output will be a
// set of edges.
//

template <concepts::VertexAndEdgeListGraph Graph, std::output_iterator<graph_edge_descriptor_t<Graph>> OutputIterator,
          concepts::ReadWriteVertexPropertyMap<Graph> Rank, concepts::ReadWriteVertexPropertyMap<Graph> Parent,
          concepts::ReadableEdgePropertyMap<Graph> Weight>
requires std::totally_ordered<property_traits_value_t<Weight>> &&
    std::convertible_to<property_traits_value_t<Parent>, graph_vertex_descriptor_t<Graph>> &&
    std::integral<property_traits_value_t<Rank>>
void kruskal_minimum_spanning_tree(const Graph& g, OutputIterator spanning_tree_edges, Rank rank, Parent parent,
                                   Weight weight) {
  if (num_vertices(g) == 0) {
    return;  // Nothing to do in this case
  }
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;

  disjoint_sets<Rank, Parent> dset(rank, parent);

  for (auto u : vertices(g)) {
    dset.make_set(u);
  }

  using weight_greater = indirect_cmp<Weight, std::greater<>>;
  std::priority_queue<Edge, std::vector<Edge>, weight_greater> q(weight_greater{weight});
  /*push all edge into Q*/
  for (auto e : edges(g)) {
    q.push(e);
  }

  while (!q.empty()) {
    Edge e = q.top();
    q.pop();
    Vertex u = dset.find_set(source(e, g));
    Vertex v = dset.find_set(target(e, g));
    if (u != v) {
      *spanning_tree_edges++ = e;
      dset.link(u, v);
    }
  }
}

template <concepts::VertexAndEdgeListGraph Graph, std::output_iterator<graph_edge_descriptor_t<Graph>> OutputIterator>
void kruskal_minimum_spanning_tree(const Graph& g, OutputIterator spanning_tree_edges) {
  std::size_t n = num_vertices(g);
  if (n == 0) {
    return;  // Nothing to do in this case
  }
  kruskal_minimum_spanning_tree(
      g, spanning_tree_edges, vector_property_map(n, get(vertex_index, g), std::size_t{0}).ref(),
      vector_property_map(n, get(vertex_index, g), graph_traits<Graph>::null_vertex()).ref(), get(edge_weight, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_KRUSKAL_MIN_SPANNING_TREE_H_
