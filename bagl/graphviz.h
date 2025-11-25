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

//
// This is an example of using the read and write functions to load and store a modified graphviz dot file.
//
// In this example, we just read a graphviz graph with some basic common properties, and then change the
// background color of the graph and write it out.
//
//  struct VertexBundle {
//    std::string name;   // Stores node id's (aka names).
//    std::string shape;  // Stores the 'shape' attribute of each node.
//  };
//  struct EdgeBundle {
//    std::string color;  // Stores the 'color' attribute of each edge.
//    std::string label;  // Stores the 'label' attribute of each edge.
//  };
//  struct GraphBundle {
//    std::string name;    // Stores the name of the graph (e.g., the "some_name" in "digraph some_name {..").
//    std::string rankdir; // Stores the 'rankdir' attribute of the graph.
//  };
//
//  using DiGraph = bagl::adjacency_list<bagl::vec_s, bagl::vec_s, bagl::directed_s,
//                                       VertexBundle, EdgeBundle, GraphBundle>;
//  DiGraph g;
//
//  // Set up the dynamic properties for all attributes we want to read / capture from the dot file.
//  bagl::dynamic_properties dp(bagl::ignore_other_properties);
//  dp.property<bagl::graph_vertex_descriptor_t<DiGraph>>("name", get(&VertexBundle::name, g));
//  dp.property<bagl::graph_vertex_descriptor_t<DiGraph>>("shape", get(&VertexBundle::shape, g));
//  dp.property<bagl::graph_edge_descriptor_t<DiGraph>>("color", get(&EdgeBundle::color, g));
//  dp.property<bagl::graph_edge_descriptor_t<DiGraph>>("label", get(&EdgeBundle::label, g));
//  dp.property<DiGraph*>("graph_name", get(&GraphBundle::name, g));
//  dp.property<DiGraph*>("rankdir", get(&GraphBundle::rankdir, g));
//
//  // Read dot file from std::cin. Storing the node id's in the "name" property.
//  if (!read_graphviz(std::cin, g, dp, "name")) {
//    exit(1);
//  }
//
//  // Add a single-valued property for the graph to set the "bgcolor" attribute to "lightblue".
//  dp.property<DiGraph*>("bgcolor", bagl::single_property_map("lightblue"));
//  // Write out the graph, with the attributes we captured and the new background color.
//  write_graphviz(std::cout, g, dp, "name", bagl::dummy_property_map{});
//

// Default writer of graphviz attributes, writes nothing.
struct default_writer {
  void operator()(std::ostream&, const std::string&, const std::any&) const {}
};

// Writer of graphviz attributes from all properties in a dynamic properties map.
class dynamic_properties_graphviz_writer {
 public:
  explicit dynamic_properties_graphviz_writer(const dynamic_properties& dp, bool with_brackets = true)
      : dp_(&dp), with_brackets_(with_brackets) {}

  void operator()(std::ostream& out, const std::string& node_id_pmap, const std::any& key) const;

 private:
  const dynamic_properties* dp_ = nullptr;
  bool with_brackets_ = true;
};

using graphviz_attr_writer = std::function<void(std::ostream&, const std::string&, const std::any&)>;

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
  if constexpr (!std::is_same_v<NodeIDMap, dummy_property_map>) {
    dp.property<graph_vertex_descriptor_t<Graph>>(node_id_pmap_name, id);
  }
  dynamic_graph_observer_wrapper<const Graph> mg(g, dp);
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
  dynamic_graph_observer_wrapper<const Graph> mg(g,
                                                 dp.property<graph_vertex_descriptor_t<Graph>>(node_id_pmap_name, id));
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
