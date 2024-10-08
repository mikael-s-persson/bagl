// Copyright (c) 2006, Stephan Diederich
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <ostream>
#include <string>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"

namespace bagl {

// Writes maximal flow problem in extended DIMACS format to an output stream.
// Vertex indices are read from an IndexMap and shiftet by 1 so their new range is [1..num_vertices(g)].

template <concepts::EdgeListGraph Graph, concepts::ReadableEdgePropertyMap<Graph> CapacityMap,
          concepts::ReadableVertexIndexMap<Graph> IndexMap>
void write_dimacs_max_flow(const Graph& g, CapacityMap capacity, IndexMap idx, graph_vertex_descriptor_t<Graph> src,
                           graph_vertex_descriptor_t<Graph> sink, std::ostream& out) {
  out << "c DIMACS max-flow file generated from bagl::write_dimacs_max_flow" << std::endl;
  // print problem description "max" and number of verts and edges
  out << "p max " << num_vertices(g) << " " << num_edges(g) << std::endl;
  // say which one is source
  out << "n " << get(idx, src) + 1 << " s" << std::endl;
  // say which one is sink
  out << "n " << get(idx, sink) + 1 << " t" << std::endl;

  // output the edges
  for (auto e : edges(g)) {
    out << "a " << idx[source(e, g)] + 1 << " " << idx[target(e, g)] + 1 << " " << get(capacity, e) << std::endl;
  }
}

}  // namespace bagl
