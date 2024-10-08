// Copyright (C) 2006  Tiago de Paula Peixoto <tiago@forked.de>
// Copyright (C) 2004  The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
//           Tiago de Paula Peixoto
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GRAPHML_H_
#define BAGL_BAGL_GRAPHML_H_

#include <any>
#include <exception>
#include <sstream>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>

#include "bagl/adjacency_list.h"
#include "bagl/dynamic_property_map.h"
#include "bagl/exception.h"
#include "bagl/graph_traits.h"

namespace bagl {

// Graph reader exceptions
struct parse_error : public graph_exception {
  explicit parse_error(std::string err) : error(std::move(err)), statement("parse error: " + error) {}
  [[nodiscard]] const char* what() const noexcept override { return statement.c_str(); }
  std::string error;
  std::string statement;
};

class mutate_graph {
 public:
  mutate_graph(const mutate_graph&) = delete;
  mutate_graph(mutate_graph&&) = delete;
  mutate_graph& operator=(const mutate_graph&) = delete;
  mutate_graph& operator=(mutate_graph&&) = delete;
  virtual ~mutate_graph() = default;

  [[nodiscard]] virtual bool is_directed() const = 0;

  virtual std::any do_add_vertex() = 0;
  virtual std::pair<std::any, bool> do_add_edge(std::any source, std::any target) = 0;

  virtual void set_graph_property(const std::string& name, const std::string& value, const std::string& value_type) = 0;

  virtual void set_vertex_property(const std::string& name, std::any vertex, const std::string& value,
                                   const std::string& value_type) = 0;

  virtual void set_edge_property(const std::string& name, std::any edge, const std::string& value,
                                 const std::string& value_type) = 0;
};

template <concepts::MutableGraph Graph>
class mutate_graph_impl : public mutate_graph {
  using vertex_descriptor = graph_vertex_descriptor_t<Graph>;
  using edge_descriptor = graph_edge_descriptor_t<Graph>;

 public:
  mutate_graph_impl(Graph& g, dynamic_properties& dp) : g_(g), dp_(dp) {}

  [[nodiscard]] bool is_directed() const override { return is_directed_graph_v<Graph>; }

  std::any do_add_vertex() override { return std::any(add_vertex(g_)); }

  std::pair<std::any, bool> do_add_edge(std::any source, std::any target) override {
    auto [e, added] = add_edge(std::any_cast<vertex_descriptor>(source), std::any_cast<vertex_descriptor>(target), g_);
    return std::make_pair(std::any(e), added);
  }

  void set_graph_property(const std::string& name, const std::string& value, const std::string& value_type) override {
    bool type_found = put_property(name, &g_, value, value_type);
    if (!type_found) {
      throw parse_error("unrecognized type \"" + value_type + "\" for key " + name);
    }
  }

  void set_vertex_property(const std::string& name, std::any vertex, const std::string& value,
                           const std::string& value_type) override {
    bool type_found = put_property(name, std::any_cast<vertex_descriptor>(vertex), value, value_type);
    if (!type_found) {
      throw parse_error("unrecognized type \"" + value_type + "\" for key " + name);
    }
  }

  void set_edge_property(const std::string& name, std::any edge, const std::string& value,
                         const std::string& value_type) override {
    bool type_found = put_property(name, std::any_cast<edge_descriptor>(edge), value, value_type);
    if (!type_found) {
      throw parse_error("unrecognized type \"" + value_type + "\" for key " + name);
    }
  }

  template <typename Key>
  bool put_property(const std::string& name, const Key& key, const std::string& value, const std::string& value_type) {
    if (value_type == "boolean") {
      put(name, dp_, key, dynamic_pmap_detail::read_value<bool>(value));
    } else if (value_type == "int") {
      put(name, dp_, key, dynamic_pmap_detail::read_value<std::int32_t>(value));
    } else if (value_type == "long") {
      put(name, dp_, key, dynamic_pmap_detail::read_value<std::int64_t>(value));
    } else if (value_type == "float") {
      put(name, dp_, key, dynamic_pmap_detail::read_value<float>(value));
    } else if (value_type == "double") {
      put(name, dp_, key, dynamic_pmap_detail::read_value<double>(value));
    } else if (value_type == "string") {
      put(name, dp_, key, dynamic_pmap_detail::read_value<std::string>(value));
    } else {
      return false;
    }
    return true;
  }

 protected:
  Graph& g_;
  dynamic_properties& dp_;
};

void read_graphml(std::istream& in, mutate_graph& g, std::size_t desired_idx);

template <concepts::MutableGraph Graph>
void read_graphml(std::istream& in, Graph& g, dynamic_properties& dp, std::size_t desired_idx = 0) {
  mutate_graph_impl<Graph> mg(g, dp);
  read_graphml(in, mg, desired_idx);
}

template <concepts::VertexAndEdgeListGraph Graph, concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
void write_graphml(std::ostream& out, const Graph& g, VertexIndexMap v_index, const dynamic_properties& dp,
                   bool ordered_vertices = false) {
  using edge_descriptor = graph_edge_descriptor_t<Graph>;
  using vertex_descriptor = graph_vertex_descriptor_t<Graph>;

  auto encode_char_entities = [](const std::string& s) {
    // Don't do anything for empty strings.
    if (s.empty()) {
      return s;
    }

    std::string r;
    // To properly round-trip spaces and not uglify the XML beyond
    // recognition, we have to encode them IF the text contains only spaces.
    if (s.find_first_not_of(' ') == std::string::npos) {
      // The first will suffice.
      r = "&#32;" + std::string(s.size() - 1, ' ');
    } else {
      for (char c : s) {
        switch (c) {
          case '<':
            r += "&lt;";
            break;
          case '>':
            r += "&gt;";
            break;
          case '&':
            r += "&amp;";
            break;
          case '"':
            r += "&quot;";
            break;
          case '\'':
            r += "&apos;";
            break;
          default:
            r += c;
            break;
        }
      }
    }
    return r;
  };

  constexpr bool graph_is_directed = is_directed_graph_v<Graph>;

  out << R"""(<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
)""";

  const std::unordered_map<std::type_index, std::string_view> type_names = {
      {std::type_index(typeid(bool)), "boolean"},       {std::type_index(typeid(std::int16_t)), "int"},
      {std::type_index(typeid(std::uint16_t)), "int"},  {std::type_index(typeid(std::int32_t)), "int"},
      {std::type_index(typeid(std::uint32_t)), "int"},  {std::type_index(typeid(std::int64_t)), "long"},
      {std::type_index(typeid(std::uint64_t)), "long"}, {std::type_index(typeid(float)), "float"},
      {std::type_index(typeid(double)), "double"},      {std::type_index(typeid(long double)), "double"},
      {std::type_index(typeid(std::string)), "string"},
  };
  std::unordered_map<std::string, std::string> graph_key_ids;
  std::unordered_map<std::string, std::string> vertex_key_ids;
  std::unordered_map<std::string, std::string> edge_key_ids;
  int key_count = 0;

  // Output keys
  for (auto [name, pmap] : dp) {
    std::string key_id = "key" + std::to_string(key_count++);
    std::string_view classification = "graph";
    if (pmap->key() == typeid(Graph*)) {
      graph_key_ids[name] = key_id;
      classification = "graph";
    } else if (pmap->key() == typeid(vertex_descriptor)) {
      vertex_key_ids[name] = key_id;
      classification = "node";
    } else if (pmap->key() == typeid(edge_descriptor)) {
      edge_key_ids[name] = key_id;
      classification = "edge";
    } else {
      continue;
    }
    std::string type_name = "string";
    auto it = type_names.find(std::type_index(pmap->value()));
    if (it != type_names.end()) {
      type_name = it->second;
    }
    out << "  <key id=\"" << encode_char_entities(key_id) << "\" for=\"" << classification << "\""
        << " attr.name=\"" << name << "\""
        << " attr.type=\"" << type_name << "\" />\n";
  }

  out << "  <graph id=\"G\" edgedefault=\"" << (graph_is_directed ? "directed" : "undirected") << "\""
      << " parse.nodeids=\"" << (ordered_vertices ? "canonical" : "free") << "\""
      << " parse.edgeids=\"canonical\" parse.order=\"nodesfirst\">\n";

  // Output graph data
  for (auto [name, pmap] : dp) {
    if (pmap->key() != typeid(Graph*)) {
      continue;
    }
    // The const_cast here is just to get typeid correct for property
    // map key; the graph should not be mutated using it.
    out << "   <data key=\"" << graph_key_ids[name] << "\">"
        << encode_char_entities(pmap->get_string(const_cast<Graph*>(&g))) << "</data>\n";
  }

  for (auto v : vertices(g)) {
    out << "    <node id=\"n" << get(v_index, v) << "\">\n";
    // Output data
    for (auto [name, pmap] : dp) {
      if (pmap->key() != typeid(vertex_descriptor)) {
        continue;
      }
      out << "      <data key=\"" << vertex_key_ids[name] << "\">" << encode_char_entities(pmap->get_string(v))
          << "</data>\n";
    }
    out << "    </node>\n";
  }

  std::size_t edge_count = 0;
  for (auto e : edges(g)) {
    out << "    <edge id=\"e" << edge_count++ << "\" source=\"n" << get(v_index, source(e, g)) << "\" target=\"n"
        << get(v_index, target(e, g)) << "\">\n";

    // Output data
    for (auto [name, pmap] : dp) {
      if (pmap->key() != typeid(edge_descriptor)) {
        continue;
      }
      out << "      <data key=\"" << edge_key_ids[name] << "\">" << encode_char_entities(pmap->get_string(e))
          << "</data>\n";
    }
    out << "    </edge>\n";
  }

  out << "  </graph>\n"
      << "</graphml>\n";
}

template <concepts::VertexAndEdgeListGraph Graph>
void write_graphml(std::ostream& out, const Graph& g, const dynamic_properties& dp, bool ordered_vertices = false) {
  write_graphml(out, g, get(vertex_index, g), dp, ordered_vertices);
}

}  // namespace bagl

#endif  // BAGL_BAGL_GRAPHML_H_
