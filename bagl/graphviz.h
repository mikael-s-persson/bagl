// Copyright 2001 University of Notre Dame.
// Copyright 2003 Jeremy Siek
// Authors: Lie-Quan Lee, Jeremy Siek, and Douglas Gregor
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GRAPHVIZ_H_
#define BAGL_BAGL_GRAPHVIZ_H_

#include <cstdio>  // for FILE
#include <fstream>
#include <iostream>
#include <map>
#include <ostream>
#include <string>
#include <tuple>
#include <type_traits>

#include "bagl/adjacency_list.h"
#include "bagl/compressed_sparse_row_graph.h"
#include "bagl/dynamic_graph.h"
#include "bagl/dynamic_property_map.h"
#include "bagl/exception.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/subgraph.h"

namespace bagl {

// Default writer of graphviz attributes, writes nothing.
struct default_writer {
  void operator()(std::ostream&, const std::any&) const {}
};

// Writer of graphviz attributes from all properties in a dynamic properties map.
class dynamic_properties_graphviz_writer {
 public:
  explicit dynamic_properties_graphviz_writer(const dynamic_properties& dp, bool with_brackets = true)
      : dp_(&dp), with_brackets_(with_brackets) {}

  void operator()(std::ostream& out, const std::any& key) const;

 private:
  const dynamic_properties* dp_ = nullptr;
  bool with_brackets_ = true;
};

using graphviz_attr_writer = std::function<void(std::ostream&, const std::any&)>;

void write_graphviz(std::ostream& out, const dynamic_graph_observer& g, graphviz_attr_writer vpw,
                    graphviz_attr_writer epw, graphviz_attr_writer gpw, const std::string& node_id_pmap);

template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> NodeIDMap>
void write_graphviz(std::ostream& out, const Graph& g, dynamic_properties& dp, graphviz_attr_writer vpw,
                    graphviz_attr_writer epw, graphviz_attr_writer gpw, const std::string& node_id_pmap_name,
                    NodeIDMap id) {
  dynamic_graph_observer_wrapper<Graph> mg(g, dp.property<graph_vertex_descriptor_t<Graph>>(node_id_pmap_name, id));
  write_graphviz(out, mg, vpw, epw, gpw, node_id_pmap_name);
}

template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> NodeIDMap>
void write_graphviz(std::ostream& out, const Graph& g, dynamic_properties& dp, const std::string& node_id_pmap_name,
                    NodeIDMap id) {
  dynamic_graph_observer_wrapper<Graph> mg(g, dp.property<graph_vertex_descriptor_t<Graph>>(node_id_pmap_name, id));
  write_graphviz(out, mg,
                 /*vertex_writer=*/dynamic_properties_graphviz_writer(dp),
                 /*edge_writer=*/dynamic_properties_graphviz_writer(dp),
                 /*graph_writer=*/dynamic_properties_graphviz_writer(dp, false), node_id_pmap_name);
}

template <concepts::VertexListGraph Graph>
void write_graphviz(std::ostream& out, const Graph& g, dynamic_properties& dp,
                    const std::string& node_id_prefix = "node_") {
  auto node_id_pmap = transform_property_map{
      [node_id_prefix](std::size_t i) { return node_id_prefix + std::to_string(i); }, get(vertex_index, g)};
  write_graphviz(out, g, dp, "node_id_prefixed_vindex", node_id_pmap);
}

template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> NodeIDMap>
void write_graphviz(std::ostream& out, const Graph& g, graphviz_attr_writer vpw,
                    graphviz_attr_writer epw,
                    graphviz_attr_writer gpw, const std::string& node_id_pmap_name,
                    NodeIDMap id) {
  dynamic_properties dp;
  dynamic_graph_observer_wrapper<Graph> mg(g, dp.property<graph_vertex_descriptor_t<Graph>>(node_id_pmap_name, id));
  write_graphviz(out, mg, vpw, epw, gpw);
}

template <concepts::VertexListGraph Graph>
std::enable_if_t<!std::is_convertible_v<const Graph&, const dynamic_graph_observer&>> write_graphviz(
    std::ostream& out, const Graph& g, graphviz_attr_writer vpw = default_writer{},
    graphviz_attr_writer epw = default_writer{}, graphviz_attr_writer gpw = default_writer{},
    const std::string& node_id_prefix = "node_") {
  auto node_id_pmap = transform_property_map{
      [node_id_prefix](std::size_t i) { return node_id_prefix + std::to_string(i); }, get(vertex_index, g)};
  write_graphviz(out, g, vpw, epw, gpw, "node_id_prefixed_vindex", node_id_pmap);
}

/////////////////////////////////////////////////////////////////////////////
// Graph reader exceptions
/////////////////////////////////////////////////////////////////////////////
struct bad_graphviz_syntax : public graph_exception {
  std::string errmsg;
  bad_graphviz_syntax(const std::string& a_errmsg) : errmsg(a_errmsg) {}
  const char* what() const noexcept override { return errmsg.c_str(); }
  ~bad_graphviz_syntax() noexcept override {}
};

bool read_graphviz(const std::string& data, dynamic_graph_mutator& mg, const std::string& node_id = "node_id");

// Parse the passed string as a GraphViz dot file.
template <typename MutableGraph>
bool read_graphviz(const std::string& data, MutableGraph& graph, dynamic_properties& dp,
                   const std::string& node_id = "node_id") {
  dynamic_graph_mutator_wrapper<MutableGraph> mg(graph, dp);
  return read_graphviz(data, mg, node_id);
}

// Parse the passed iterator range as a GraphViz dot file.
template <typename InputIterator, typename MutableGraph>
bool read_graphviz(InputIterator user_first, InputIterator user_last, MutableGraph& graph, dynamic_properties& dp,
                   const std::string& node_id = "node_id") {
  return read_graphviz(std::string(user_first, user_last), graph, dp, node_id);
}

// Parse the passed stream as a GraphViz dot file.
template <typename MutableGraph>
bool read_graphviz(std::istream& in, MutableGraph& graph, dynamic_properties& dp,
                   const std::string& node_id = "node_id") {
  typedef std::istream_iterator<char> is_t;
  in >> std::noskipws;
  return read_graphviz(is_t(in), is_t(), graph, dp, node_id);
}

}  // namespace bagl

#endif  // BAGL_BAGL_GRAPHVIZ_H_
