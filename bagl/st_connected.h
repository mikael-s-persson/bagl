// Copyright (C) 2006 The Trustees of Indiana University.
//  Authors: Douglas Gregor, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_ST_CONNECTED_H_
#define BAGL_BAGL_ST_CONNECTED_H_

#include <queue>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"
#include "bagl/two_bit_color_map.h"

namespace bagl {

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> ColorMap>
requires concepts::IncidenceGraph<Graph>
bool st_connected(const Graph& g, graph_vertex_descriptor_t<Graph> s, graph_vertex_descriptor_t<Graph> t,
                  ColorMap color) {
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  // Set all vertices to white (unvisited)
  for (auto v : vertices(g)) {
    put(color, v, Color::white());
  }

  // Vertices found from the source are grey
  put(color, s, Color::gray());

  // Vertices found from the target are green
  put(color, t, Color::green());
  std::queue<Vertex> q;
  q.push(s);
  q.push(t);

  while (!q.empty()) {
    Vertex u = q.top();
    q.pop();
    auto u_color = get(color, u);

    for (auto e : out_edges(u, g)) {
      Vertex v = target(e, g);
      auto v_color = get(color, v);
      if (v_color == Color::white()) {
        // We have not seen "v" before; mark it with the same color as u
        put(color, v, u_color);

        // Push it on the queue
        q.push(v);
      } else if (v_color != Color::black() && u_color != v_color) {
        // Colors have collided. We're done!
        return true;
      }
    }
    // u is done, so mark it black
    put(color, u, Color::black());
  }

  return false;
}

template <concepts::VertexListGraph Graph>
requires concepts::IncidenceGraph<Graph>
bool st_connected(const Graph& g, graph_vertex_descriptor_t<Graph> s, graph_vertex_descriptor_t<Graph> t) {
  return st_connected(g, s, t, two_bit_color_map(num_vertices(g), get(vertex_index, g)).ref());
}

template <concepts::VertexListGraph Graph, concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
requires concepts::IncidenceGraph<Graph>
bool st_connected(const Graph& g, graph_vertex_descriptor_t<Graph> s, graph_vertex_descriptor_t<Graph> t,
                  VertexIndexMap v_index) {
  return st_connected(g, s, t, two_bit_color_map(num_vertices(g), v_index).ref());
}

}  // namespace bagl

#endif  // BAGL_BAGL_ST_CONNECTED_H_
