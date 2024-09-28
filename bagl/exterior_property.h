// Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_EXTERIOR_PROPERTY_H_
#define BAGL_BAGL_EXTERIOR_PROPERTY_H_

#include <vector>

#include "bagl/container_property_map.h"
#include "bagl/graph_traits.h"
#include "bagl/matrix_property_map.h"

namespace bagl {
namespace exterior_property_detail {
// The vector matrix provides a little abstraction over vector
// types that makes matrices easier to work with.
template <typename Value>
struct vector_matrix {
  using container_type = std::vector<Value>;
  using matrix_type = std::vector<container_type>;

  using value_type = container_type;
  using reference = container_type&;
  using const_reference = const container_type&;
  using pointer = container_type*;
  using size_type = typename matrix_type::size_type;

  // Instantiate the matrix over n elements (creates an n by n matrix).
  // The graph has to be passed in order to ensure the index maps
  // are constructed correctly when returning indexible elements.
  explicit vector_matrix(size_type n) : matrix(n, container_type(n)) {}

  reference operator[](size_type n) { return matrix[n]; }
  const_reference operator[](size_type n) const { return matrix[n]; }

  matrix_type matrix;
};
}  // namespace exterior_property_detail

// The exterior_property metafunction defines an appropriate set of types for
// creating an exterior property. An exterior property is comprised of a both
// a container and a property map that acts as its abstraction. An extension
// of this metafunction will select an appropriate "matrix" property that
// records values for pairs of vertices.
//
// @todo This does not currently support the ability to define exterior
// properties for graph types that do not model the IndexGraph concepts. A
// solution should not be especially difficult, but will require an extension
// of type traits to affect the type selection.
template <typename Graph, typename Key, typename Value>
struct exterior_property {
  using key_type = Key;
  using value_type = Value;

  using container_type = std::vector<Value>;
  using map_type = container_property_map<Graph, Key, container_type>;

  using matrix_type = exterior_property_detail::vector_matrix<Value>;
  using matrix_map_type = matrix_property_map<Graph, Key, matrix_type>;

  exterior_property() = delete;
};

// Define a the container and property map types requried to create an exterior
// vertex property for the given value type. The Graph parameter is required to
// model the VertexIndexGraph concept.
template <typename Graph, typename Value>
struct exterior_vertex_property {
  using property_type = exterior_property<Graph, graph_vertex_descriptor_t<Graph>, Value>;
  using key_type = typename property_type::key_type;
  using value_type = typename property_type::value_type;
  using container_type = typename property_type::container_type;
  using map_type = typename property_type::map_type;
  using matrix_type = typename property_type::matrix_type;
  using matrix_map_type = typename property_type::matrix_map_type;
};

// Define a the container and property map types requried to create an exterior
// edge property for the given value type. The Graph parameter is required to
// model the EdgeIndexGraph concept.
template <typename Graph, typename Value>
struct exterior_edge_property {
  using property_type = exterior_property<Graph, graph_edge_descriptor_t<Graph>, Value>;
  using key_type = typename property_type::key_type;
  using value_type = typename property_type::value_type;
  using container_type = typename property_type::container_type;
  using map_type = typename property_type::map_type;
  using matrix_type = typename property_type::matrix_type;
  using matrix_map_type = typename property_type::matrix_map_type;
};

}  // namespace bagl

#endif  // BAGL_BAGL_EXTERIOR_PROPERTY_H_
