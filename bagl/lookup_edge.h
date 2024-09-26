// Copyright 2009 Trustees of Indiana University
// Author: Jeremiah Willcock
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_LOOKUP_EDGE_H_
#define BAGL_BAGL_LOOKUP_EDGE_H_

#include <utility>

#include "bagl/graph_concepts.h"

// lookup_edge: a function that acts like edge() but falls back to out_edges()
// and a search when edge() is not provided.

namespace bagl {

template <concepts::Graph G>
std::pair<graph_edge_descriptor_t<G>, bool> lookup_edge(graph_vertex_descriptor_t<G> src,
                                                        graph_vertex_descriptor_t<G> tgt, const G& g) {
  if constexpr (concepts::AdjacencyMatrix<G>) {
    return edge(src, tgt, g);
  } else {
    for (auto e : out_edges(src, g)) {
      if (target(e, g) == tgt) {
        return {e, true};
      }
    }
    return {{}, false};
  }
}

}  // namespace bagl

#endif  // BAGL_BAGL_LOOKUP_EDGE_H_
