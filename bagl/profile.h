// Copyright 2002 Marc Wintermantel (wintermantel@even-ag.ch)
// ETH Zurich, Center of Structure Technologies
// (https://web.archive.org/web/20050307090307/http://www.structures.ethz.ch/)
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PROFILE_H_
#define BAGL_BAGL_PROFILE_H_

#include <limits>

#include "bagl/bandwidth.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"

namespace bagl {

template <concepts::VertexListGraph G, concepts::ReadableVertexIndexMap<G> VertexIndex>
auto profile(const G& g, VertexIndex index) {
  std::size_t b = 0;
  for (auto v : vertices(g)) {
    b += ith_bandwidth(v, g, index) + 1;
  }

  return b;
}

template <concepts::VertexListGraph G>
auto profile(const G& g) {
  return profile(g, get(vertex_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_PROFILE_H_
