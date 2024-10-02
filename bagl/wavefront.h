// Copyright 2002 Marc Wintermantel (wintermantel@even-ag.ch)
// ETH Zurich, Center of Structure Technologies
// (https://web.archive.org/web/20050307090307/http://www.structures.ethz.ch/)
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_WAVEFRONT_H_
#define BAGL_BAGL_WAVEFRONT_H_

#include <algorithm>  // for std::min and std::max
#include <cmath>
#include <vector>

#include "bagl/bandwidth.h"
#include "bagl/graph_traits.h"

namespace bagl {

template <concepts::VertexListGraph G, concepts::ReadableVertexIndexMap<G> VertexIndexMap>
std::size_t ith_wavefront(graph_vertex_descriptor_t<G> i, const G& g, VertexIndexMap index) {
  std::size_t b = 1;
  const std::size_t index_i = index[i];
  std::vector<bool> rows_active(num_vertices(g), false);

  rows_active[index_i] = true;

  for (auto v : vertices(g)) {
    if (index[v] > index_i) {
      continue;
    }
    for (auto e : out_edges(v, g)) {
      auto w = target(e, g);
      if ((index[w] >= index_i) && (!rows_active[index[w]])) {
        b++;
        rows_active[index[w]] = true;
      }
    }
  }

  return b;
}

template <concepts::VertexListGraph G>
std::size_t ith_wavefront(graph_vertex_descriptor_t<G> v, const G& g) {
  return ith_wavefront(v, g, get(vertex_index, g));
}

template <concepts::VertexListGraph G, concepts::ReadableVertexIndexMap<G> VertexIndexMap>
std::size_t max_wavefront(const G& g, VertexIndexMap index) {
  std::size_t b = 0;
  for (auto v : vertices(g)) {
    b = std::max(b, ith_wavefront(v, g, index));
  }
  return b;
}

template <concepts::VertexListGraph G>
std::size_t max_wavefront(const G& g) {
  return max_wavefront(g, get(vertex_index, g));
}

template <concepts::VertexListGraph G, concepts::ReadableVertexIndexMap<G> VertexIndexMap>
double aver_wavefront(const G& g, VertexIndexMap index) {
  double b = 0;
  for (auto v : vertices(g)) {
    b += ith_wavefront(v, g, index);
  }

  b /= num_vertices(g);
  return b;
}

template <concepts::VertexListGraph G>
double aver_wavefront(const G& g) {
  return aver_wavefront(g, get(vertex_index, g));
}

template <concepts::VertexListGraph G, concepts::ReadableVertexIndexMap<G> VertexIndexMap>
double rms_wavefront(const G& g, VertexIndexMap index) {
  double b = 0;
  for (auto v : vertices(g)) {
    b += std::pow(static_cast<double>(ith_wavefront(v, g, index)), 2.0);
  }

  b /= num_vertices(g);

  return std::sqrt(b);
}

template <concepts::VertexListGraph G>
double rms_wavefront(const G& g) {
  return rms_wavefront(g, get(vertex_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_WAVEFRONT_H_
