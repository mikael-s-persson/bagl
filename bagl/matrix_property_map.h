// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_MATRIX_PROPERTY_MAP_H_
#define BAGL_BAGL_MATRIX_PROPERTY_MAP_H_

#include "bagl/container_property_map.h"

namespace bagl {

// This property map is built specifically for property maps over
// matrices. Like the basic property map over a container, this builds
// the property abstraction over a matrix (usually a vector of vectors)
// and returns property maps over the nested containers.
template <typename Graph, typename Key, typename Matrix>
struct matrix_property_map : put_get_helper<matrix_property_map<Graph, Key, Matrix> > {
  // abstract the indexing keys
  using indexer_type = typename index_detail::choose_indexer<Graph, Key>::indexer_type;

  // aliases for the nested container and its corresponding map
  using container_type = typename Matrix::value_type;
  using map_type = container_property_map<Graph, Key, container_type>;

  using key_type = Key;

  // This property map doesn't really provide access to nested containers,
  // but returns property maps over them. Since property maps are all
  // copy-constructible (or should be anyways), we never return references.
  // As such, this property is only readable, but not writable. Curiously,
  // the inner property map is actually an lvalue pmap.
  using value_type = map_type;
  using reference = map_type;
  using category = readable_property_map_tag;

  matrix_property_map() = default;

  matrix_property_map(Matrix& m, const Graph& g) : matrix_(&m), graph_(&g) {}

  reference operator[](key_type k) const {
    auto x = indexer_type::index(k, *graph_);
    return map_type((*matrix_)[x], *graph_);
  }

 private:
  Matrix* matrix_ = nullptr;
  const Graph* graph_ = nullptr;
};
}  // namespace bagl

#endif
