// Copyright 2013 Maciej Piechotka
// Authors: Maciej Piechotka
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_EDGE_COLORING_H_
#define BAGL_BAGL_EDGE_COLORING_H_

#include <algorithm>
#include <limits>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/numeric_values.h"
#include "bagl/properties.h"

// This algorithm is to find coloring of an edges
//
// Reference:
//
// Misra, J., & Gries, D. (1992). A constructive proof of Vizing's
// theorem. In Information Processing Letters.

namespace bagl {
namespace edge_coloring_detail {
template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> ColorMap>
bool is_free(const G& g, ColorMap color, graph_vertex_descriptor_t<G> u, property_traits_value_t<ColorMap> free_color) {
  if (free_color == numeric_values<decltype(free_color)>::infinity()) {
    return false;
  }
  for (auto e : out_edges(u, g)) {
    if (get(color, e) == free_color) {
      return false;
    }
  }
  return true;
}

template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> ColorMap>
std::vector<graph_vertex_descriptor_t<G>> maximal_fan(const G& g, ColorMap color, graph_vertex_descriptor_t<G> x,
                                                      graph_vertex_descriptor_t<G> y) {
  using Vertex = graph_vertex_descriptor_t<G>;
  std::vector<Vertex> fan;
  fan.push_back(y);
  bool extended = false;
  do {
    extended = false;
    for (auto e : out_edges(x, g)) {
      Vertex v = target(e, g);
      if (is_free(g, color, fan.back(), get(color, e)) && std::find(fan.begin(), fan.end(), v) == fan.end()) {
        fan.push_back(v);
        extended = true;
      }
    }
  } while (extended);
  return fan;
}
template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> ColorMap>
auto find_free_color(const G& g, ColorMap color, graph_vertex_descriptor_t<G> u) {
  auto c = numeric_values<property_traits_value_t<ColorMap>>::zero();
  while (!is_free(g, color, u, c)) {
    c++;
  }
  return c;
}

template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> ColorMap>
void invert_cd_path(const G& g, ColorMap color, graph_vertex_descriptor_t<G> x, graph_edge_descriptor_t<G> eold,
                    property_traits_value_t<ColorMap> c, property_traits_value_t<ColorMap> d) {
  put(color, eold, d);
  for (auto e : out_edges(x, g)) {
    if (get(color, e) == d && e != eold) {
      invert_cd_path(g, color, target(e, g), e, d, c);
      return;
    }
  }
}

template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> ColorMap>
void invert_cd_path(const G& g, ColorMap color, graph_vertex_descriptor_t<G> x, property_traits_value_t<ColorMap> c,
                    property_traits_value_t<ColorMap> d) {
  for (auto e : out_edges(x, g)) {
    if (get(color, e) == d) {
      invert_cd_path(g, color, target(e, g), e, d, c);
      return;
    }
  }
}

template <concepts::IncidenceGraph G, concepts::ReadWriteEdgePropertyMap<G> ColorMap,
          std::ranges::input_range VertexRange>
requires concepts::AdjacencyMatrix<G>
void rotate_fan(const G& g, ColorMap color, graph_vertex_descriptor_t<G> x, VertexRange v_rg) {
  if (v_rg.empty()) {
    return;
  }
  auto begin = v_rg.begin();
  auto previous = edge(x, *begin, g).first;
  for (begin++; begin != v_rg.end(); begin++) {
    auto current = edge(x, *begin, g).first;
    put(color, previous, get(color, current));
    previous = current;
  }
}
}  // namespace edge_coloring_detail

template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> ColorMap>
property_traits_value_t<ColorMap> color_edge(const G& g, ColorMap color, graph_edge_descriptor_t<G> e) {
  auto x = source(e, g);
  auto y = target(e, g);
  auto fan = edge_coloring_detail::maximal_fan(g, color, x, y);
  auto c = edge_coloring_detail::find_free_color(g, color, x);
  auto d = edge_coloring_detail::find_free_color(g, color, fan.back());
  edge_coloring_detail::invert_cd_path(g, color, x, c, d);
  auto w = std::find_if(fan.begin(), fan.end(), [g, color, d](auto u) { return is_free(g, color, u, d); });
  edge_coloring_detail::rotate_fan(g, color, x, fan.begin(), w + 1);
  put(color, edge(x, *w, g).first, d);
  return std::max(c, d);
}

template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> ColorMap>
requires concepts::EdgeListGraph<G> property_traits_value_t<ColorMap> edge_coloring(const G& g, ColorMap color) {
  using ColorValue = property_traits_value_t<ColorMap>;
  for (auto e : edges(g)) {
    put(color, e, numeric_values<ColorValue>::infinity());
  }
  auto colors = numeric_values<ColorValue>::zero();
  for (auto e : edges(g)) {
    colors = std::max(colors, color_edge(g, color, e) + 1);
  }
  return colors;
}

}  // namespace bagl

#endif  // BAGL_BAGL_EDGE_COLORING_H_
