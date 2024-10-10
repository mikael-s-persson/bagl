// Copyright (C) 2006  Tiago de Paula Peixoto <tiago@forked.de>
// Copyright (C) 2004  The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
//           Tiago de Paula Peixoto
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GRAPHML_H_
#define BAGL_BAGL_GRAPHML_H_

#include <cstddef>
#include <exception>
#include <sstream>

#include "bagl/adjacency_list.h"
#include "bagl/dynamic_graph.h"
#include "bagl/dynamic_property_map.h"
#include "bagl/exception.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/property_map.h"

namespace bagl {

// Graph reader exceptions
struct parse_error : public graph_exception {
  explicit parse_error(std::string err) : error(std::move(err)), statement("parse error: " + error) {}
  [[nodiscard]] const char* what() const noexcept override { return statement.c_str(); }
  std::string error;
  std::string statement;
};

void read_graphml(std::istream& in, dynamic_graph_mutator& g, std::size_t desired_idx);

template <concepts::MutableGraph Graph>
void read_graphml(std::istream& in, Graph& g, dynamic_properties& dp, std::size_t desired_idx = 0) {
  dynamic_graph_mutator_wrapper<Graph> mg(g, dp);
  read_graphml(in, mg, desired_idx);
}

void write_graphml(std::ostream& out, const dynamic_graph_observer& g, bool ordered_vertices = false);

template <concepts::VertexAndEdgeListGraph Graph, concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
void write_graphml(std::ostream& out, const Graph& g, const dynamic_properties& dp, bool ordered_vertices = false) {
  dynamic_graph_observer_wrapper<Graph> mg(g, dp);
  write_graphml(out, mg, ordered_vertices);
}

}  // namespace bagl

#endif  // BAGL_BAGL_GRAPHML_H_
