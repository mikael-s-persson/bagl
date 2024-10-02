// Copyright 2002 Indiana University.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GRAPH_SELECTORS_H_
#define BAGL_BAGL_GRAPH_SELECTORS_H_

#include <type_traits>

namespace bagl {

//===========================================================================
// Selectors for the Directed template parameter of adjacency_list
// and adjacency_matrix.

struct directed_s {
  static constexpr bool is_directed = true;
  static constexpr bool is_bidir = false;
  using is_directed_t = std::true_type;
  using is_bidir_t = std::false_type;
};
struct undirected_s {
  static constexpr bool is_directed = false;
  static constexpr bool is_bidir = false;
  using is_directed_t = std::false_type;
  using is_bidir_t = std::false_type;
};
struct bidirectional_s {
  static constexpr bool is_directed = true;
  static constexpr bool is_bidir = true;
  using is_directed_t = std::true_type;
  using is_bidir_t = std::true_type;
};

// A type selector that denotes the use of some default value.
struct default_s {};

}  // namespace bagl

#endif  // BAGL_BAGL_GRAPH_SELECTORS_H_
