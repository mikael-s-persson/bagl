// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_TRANSPOSE_GRAPH_H_
#define BAGL_BAGL_TRANSPOSE_GRAPH_H_

#include "bagl/copy.h"
#include "bagl/graph_traits.h"
#include "bagl/reverse_graph.h"

namespace bagl {

template <concepts::VertexListGraph GIn, concepts::MutableGraph GOut>
void transpose_graph(const GIn& g, GOut& g_t) {
  reverse_graph r(g);
  copy_graph(r, g_t);
}

}  // namespace bagl

#endif  // BAGL_BAGL_TRANSPOSE_GRAPH_H_
