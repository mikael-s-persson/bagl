// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

// Revision History:
//   17 March 2006: Fixed a bug: when updating the degree a vertex
//                  could be moved to a wrong bucket. (Roman Dementiev)
//

#ifndef BAGL_BAGL_SMALLEST_LAST_ORDERING_H_
#define BAGL_BAGL_SMALLEST_LAST_ORDERING_H_
/*
   The smallest-last ordering is defined for the loopless graph G with
   vertices a(j), j = 1,2,...,n where a(j) is the j-th column of A and
   with edge (a(i),a(j)) if and only if columns i and j have a
   non-zero in the same row position.  The smallest-last ordering is
   determined recursively by letting list(k), k = n,...,1 be a column
   with least degree in the subgraph spanned by the un-ordered
   columns.
 */
#include <algorithm>
#include <vector>

#include "bagl/bucket_sorter.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"

namespace bagl {

template <class VertexListGraph, class Order, class Degree, class Marker>
void smallest_last_vertex_ordering(const VertexListGraph& G, Order order, Degree degree, Marker marker) {
  using Vertex = graph_vertex_descriptor_t<VertexListGraph>;

  using ID = property_map_t<VertexListGraph, vertex_index_t>;
  using BucketSorter = bucket_sorter<std::size_t, Vertex, Degree, ID>;

  const std::size_t num = num_vertices(G);
  BucketSorter degree_bucket_sorter(num, num, degree, get(vertex_index, G));

  smallest_last_vertex_ordering(G, order, degree, marker, degree_bucket_sorter);
}

template <class VertexListGraph, class Order, class Degree, class Marker, class BucketSorter>
void smallest_last_vertex_ordering(const VertexListGraph& G, Order order, Degree degree, Marker marker,
                                   BucketSorter& degree_buckets) {
  using Vertex = graph_vertex_descriptor_t<VertexListGraph>;

  const std::size_t num = num_vertices(G);
  for (auto v : vertices(G)) {
    put(marker, v, num);
    put(degree, v, out_degree(v, G));
    degree_buckets.push(v);
  }

  std::size_t minimum_degree = 0;
  std::size_t current_order = num - 1;

  while (true) {
    auto minimum_degree_stack = degree_buckets[minimum_degree];
    while (minimum_degree_stack.empty()) {
      minimum_degree_stack = degree_buckets[++minimum_degree];
    }

    Vertex node = minimum_degree_stack.top();
    put(order, current_order, node);

    if (current_order == 0) {
      // find all vertices
      break;
    }

    minimum_degree_stack.pop();
    put(marker, node, 0);  // node has been ordered.

    for (auto v : adjacent_vertices(node, G)) {
      if (get(marker, v) > current_order) {
        // v is unordered vertex
        // mark the columns adjacent to node
        put(marker, v, current_order);

        // delete v from the bucket sorter
        degree_buckets.remove(v);

        // It is possible minimum degree goes down
        // Here we keep tracking it.
        put(degree, v, get(degree, v) - 1);
        minimum_degree = std::min(minimum_degree, get(degree, v));

        // reinsert *v in the bucket sorter with the new degree
        degree_buckets.push(v);
      }
    }

    current_order--;
  }

  // at this point, order[i] = v_i;
}

template <class VertexListGraph, class Order>
void smallest_last_vertex_ordering(const VertexListGraph& G, Order order) {
  smallest_last_vertex_ordering(G, order, vector_property_map(num_vertices(G), get(vertex_index, G), std::size_t{0}),
                                vector_property_map(num_vertices(G), get(vertex_index, G), std::size_t{0}));
}

template <class VertexListGraph>
std::vector<typename graph_traits<VertexListGraph>::vertex_descriptor> smallest_last_vertex_ordering(
    const VertexListGraph& G) {
  std::vector<graph_vertex_descriptor_t<VertexListGraph> > o(num_vertices(G));
  smallest_last_vertex_ordering(G, make_iterator_property_map(o.begin(), typed_identity_property_map<std::size_t>()));
  return o;
}

}  // namespace bagl

#endif  // BAGL_BAGL_SMALLEST_LAST_ORDERING_H_
