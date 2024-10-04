// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Copyright 2004 The Trustees of Indiana University
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_SEQUENTIAL_VERTEX_COLORING_H_
#define BAGL_BAGL_SEQUENTIAL_VERTEX_COLORING_H_

#include <limits>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"

/* This algorithm is to find coloring of a graph

   Algorithm:
   Let G = (V,E) be a graph with vertices (somehow) ordered v_1, v_2, ...,
   v_n. For k = 1, 2, ..., n the sequential algorithm assigns v_k to the
   smallest possible color.

   Reference:

   Thomas F. Coleman and Jorge J. More, Estimation of sparse Jacobian
   matrices and graph coloring problems. J. Numer. Anal. V20, P187-209, 1983

   v_k is stored as o[k] here.

   The color of the vertex v will be stored in color[v].
   i.e., vertex v belongs to coloring color[v] */

namespace bagl {

template <concepts::VertexListGraph G, concepts::ReadablePropertyMap<std::size_t> OrderPA,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap>
requires concepts::AdjacencyGraph<G>
auto sequential_vertex_coloring(const G& g, OrderPA order, ColorMap color) {
  using Vertex = graph_vertex_descriptor_t<G>;

  std::size_t max_color = 0;
  const std::size_t n_v = num_vertices(g);

  // We need to keep track of which colors are used by
  // adjacent vertices. We do this by marking the colors
  // that are used. The mark array contains the mark
  // for each color. The length of mark is the
  // number of vertices since the maximum possible number of colors
  // is the number of vertices.
  std::vector<std::size_t> mark(n_v, std::numeric_limits<std::size_t>::max());

  // Initialize colors
  for (auto v : vertices(g)) {
    put(color, v, n_v - 1);
  }

  // Determine the color for every vertex one by one
  for (std::size_t i = 0; i < n_v; i++) {
    Vertex current = get(order, i);

    // Mark the colors of vertices adjacent to current.
    // i can be the value for marking since i increases successively
    for (auto v : adjacent_vertices(current, g)) {
      mark[get(color, *v)] = i;
    }

    // Next step is to assign the smallest un-marked color
    // to the current vertex.
    std::size_t j = 0;

    // Scan through all useable colors, find the smallest possible
    // color that is not used by neighbors.  Note that if mark[j]
    // is equal to i, color j is used by one of the current vertex's
    // neighbors.
    while (j < max_color && mark[j] == i) {
      ++j;
    }

    if (j == max_color) {
      // All colors are used up. Add one more color
      ++max_color;
    }

    // At this point, j is the smallest possible color
    put(color, current, j);  // Save the color of vertex current
  }

  return max_color;
}

template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> ColorMap>
requires concepts::AdjacencyGraph<G>
auto sequential_vertex_coloring(const G& g, ColorMap color) {
  using Vertex = graph_vertex_descriptor_t<G>;
  auto v_rg = vertices(g);
  std::vector<Vertex> order(v_rg.begin(), v_rg.end());
  return sequential_vertex_coloring(g, iterator_property_map(order.begin(), identity_property_map()), color);
}

}  // namespace bagl

#endif  // BAGL_BAGL_SEQUENTIAL_VERTEX_COLORING_H_
