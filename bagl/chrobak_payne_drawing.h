// Copyright 2007 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_CHROBAK_PAYNE_DRAWING_H_
#define BAGL_BAGL_CHROBAK_PAYNE_DRAWING_H_

#include <list>
#include <stack>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"

namespace bagl {

namespace chrobak_payne_detail {

template <typename Graph, typename VertexToVertexMap, typename VertexTo1DCoordMap>
void accumulate_offsets(graph_vertex_descriptor_t<Graph> v, std::size_t offset, const Graph& g, VertexTo1DCoordMap x,
                        VertexTo1DCoordMap delta_x, VertexToVertexMap left, VertexToVertexMap right) {
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using StackEntry = std::pair<Vertex, std::size_t>;
  std::stack<StackEntry> st;
  st.emplace(v, offset);
  while (!st.empty()) {
    std::tie(v, offset) = st.top();
    st.pop();
    if (v != graph_traits<Graph>::null_vertex()) {
      x[v] += delta_x[v] + offset;
      st.emplace(left[v], x[v]);
      st.emplace(right[v], x[v]);
    }
  }
}

}  // namespace chrobak_payne_detail

template <concepts::VertexListGraph G, concepts::LvalueVertexPropertyMap<G> PlanarEmbedding,
          std::ranges::forward_range FwdRange, concepts::MutableLvalueVertexPropertyMap<G> GridPositionMap,
          concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
void chrobak_payne_straight_line_drawing(const G& g, PlanarEmbedding embedding, const FwdRange& ordering,
                                         GridPositionMap drawing, VertexIndexMap vm) {
  using Vertex = graph_vertex_descriptor_t<G>;
  const Vertex null_v = graph_traits<G>::null_vertex();

  const std::size_t n_verts = num_vertices(g);
  auto left = make_vector_property_map(n_verts, vm, null_v);
  auto right = make_vector_property_map(n_verts, vm, null_v);
  auto seen_as_right = make_vector_property_map(n_verts, vm, std::size_t{0});
  auto seen = make_vector_property_map(n_verts, vm, std::size_t{0});
  auto delta_x = make_vector_property_map(n_verts, vm, std::size_t{0});
  auto y = make_vector_property_map(n_verts, vm, std::size_t{0});
  auto x = make_vector_property_map(n_verts, vm, std::size_t{0});
  auto installed = make_vector_property_map(n_verts, vm, false);

  std::size_t timestamp = 1;
  std::vector<Vertex> installed_neighbors;

  auto itr = ordering.begin();
  Vertex v1 = *itr;
  ++itr;
  Vertex v2 = *itr;
  ++itr;
  Vertex v3 = *itr;
  ++itr;

  delta_x[v2] = 1;
  delta_x[v3] = 1;

  y[v1] = 0;
  y[v2] = 0;
  y[v3] = 1;

  right[v1] = v3;
  right[v3] = v2;

  installed[v1] = installed[v2] = installed[v3] = true;

  for (auto itr_end = ordering.end(); itr != itr_end; ++itr) {
    Vertex v = *itr;

    // First, find the leftmost and rightmost neighbor of v on the outer
    // cycle of the embedding.
    // Note: since we're moving clockwise through the edges adjacent to v,
    // we're actually moving from right to left among v's neighbors on the
    // outer face (since v will be installed above them all) looking for
    // the leftmost and rightmost installed neigbhors

    Vertex leftmost = null_v;
    Vertex rightmost = null_v;

    installed_neighbors.clear();

    Vertex prev_vertex = null_v;
    for (auto e : embedding[v]) {
      auto curr_vertex = source(e, g) == v ? target(e, g) : source(e, g);

      // Skip any self-loops or parallel edges
      if (curr_vertex == v || curr_vertex == prev_vertex) {
        continue;
      }

      if (installed[curr_vertex]) {
        seen[curr_vertex] = timestamp;

        if (right[curr_vertex] != null_v) {
          seen_as_right[right[curr_vertex]] = timestamp;
        }
        installed_neighbors.push_back(curr_vertex);
      }

      prev_vertex = curr_vertex;
    }

    for (auto u : installed_neighbors) {
      if (right[u] == null_v || seen[right[u]] != timestamp) {
        rightmost = u;
      }
      if (seen_as_right[u] != timestamp) {
        leftmost = u;
      }
    }

    ++timestamp;

    // stretch gaps
    ++delta_x[right[leftmost]];
    ++delta_x[rightmost];

    // adjust offsets
    std::size_t delta_p_q = 0;
    Vertex stopping_vertex = right[rightmost];
    for (Vertex temp = right[leftmost]; temp != stopping_vertex; temp = right[temp]) {
      delta_p_q += delta_x[temp];
    }

    delta_x[v] = ((y[rightmost] + delta_p_q) - y[leftmost]) / 2;
    y[v] = y[leftmost] + delta_x[v];
    delta_x[rightmost] = delta_p_q - delta_x[v];

    bool leftmost_and_rightmost_adjacent = right[leftmost] == rightmost;
    if (!leftmost_and_rightmost_adjacent) {
      delta_x[right[leftmost]] -= delta_x[v];
    }

    // install v
    if (!leftmost_and_rightmost_adjacent) {
      left[v] = right[leftmost];
      Vertex next_to_rightmost = null_v;
      for (Vertex temp = leftmost; temp != rightmost; temp = right[temp]) {
        next_to_rightmost = temp;
      }

      right[next_to_rightmost] = null_v;
    } else {
      left[v] = null_v;
    }

    right[leftmost] = v;
    right[v] = rightmost;
    installed[v] = true;
  }

  chrobak_payne_detail::accumulate_offsets(*ordering.begin(), 0, g, x, delta_x, left, right);

  for (auto v : vertices(g)) {
    drawing[v].x = x[v];
    drawing[v].y = y[v];
  }
}

template <concepts::VertexListGraph G, concepts::LvalueVertexPropertyMap<G> PlanarEmbedding,
          std::ranges::forward_range FwdRange, concepts::MutableLvalueVertexPropertyMap<G> GridPositionMap>
void chrobak_payne_straight_line_drawing(const G& g, PlanarEmbedding embedding, const FwdRange& ordering,
                                         GridPositionMap drawing) {
  chrobak_payne_straight_line_drawing(g, embedding, ordering, drawing, get(vertex_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_CHROBAK_PAYNE_DRAWING_H_
