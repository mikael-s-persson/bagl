
#include <any>
#include <ranges>
#include <typeindex>

#include "bagl/graphml.h"
#include "bagl/dynamic_graph.h"

namespace bagl {

namespace {
std::string encode_char_entities(const std::string& s) {
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
}
}  // namespace

void write_graphml(std::ostream& out, const dynamic_graph_observer& g, bool ordered_vertices) {
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
  for (auto [name, pmap] : g.get_properties()) {
    std::string key_id = "key" + std::to_string(key_count++);
    std::string_view key_class_str = "graph";
    auto key_class = g.classify_key(pmap->key());
    if (key_class == dynamic_graph_observer::key_type::graph) {
      graph_key_ids[name] = key_id;
      key_class_str = "graph";
    } else if (key_class == dynamic_graph_observer::key_type::vertex) {
      vertex_key_ids[name] = key_id;
      key_class_str = "node";
    } else if (key_class == dynamic_graph_observer::key_type::edge) {
      edge_key_ids[name] = key_id;
      key_class_str = "edge";
    } else {
      continue;
    }
    std::string type_name = "string";
    auto it = type_names.find(std::type_index(pmap->value()));
    if (it != type_names.end()) {
      type_name = it->second;
    }
    out << "  <key id=\"" << encode_char_entities(key_id) << "\" for=\"" << key_class_str << "\""
        << " attr.name=\"" << name << "\""
        << " attr.type=\"" << type_name << "\" />\n";
  }

  out << "  <graph id=\"G\" edgedefault=\"" << (g.is_directed() ? "directed" : "undirected") << "\""
      << " parse.nodeids=\"" << (ordered_vertices ? "canonical" : "free") << "\""
      << " parse.edgeids=\"canonical\" parse.order=\"nodesfirst\">\n";

  // Output graph data
  for (auto [name, pmap] : g.get_properties()) {
    if (!g.is_graph_key(pmap->key())) {
      continue;
    }
    out << "   <data key=\"" << graph_key_ids[name] << "\">"
        << encode_char_entities(pmap->get_string(&g)) << "</data>\n";
  }

  for (auto v : g.get_vertices()) {
    out << "    <node id=\"n" << g.get_index_of(v) << "\">\n";
    // Output data
    for (auto [name, pmap] : g.get_properties()) {
      if (!g.is_vertex_key(pmap->key())) {
        continue;
      }
      out << "      <data key=\"" << vertex_key_ids[name] << "\">" << encode_char_entities(pmap->get_string(v))
          << "</data>\n";
    }
    out << "    </node>\n";
  }

  for (auto e : g.get_edges()) {
    out << "    <edge id=\"e" << g.get_index_of(e) << "\" source=\"n" << g.get_index_of(g.get_source(e)) << "\" target=\"n"
        << g.get_index_of(g.get_target(e)) << "\">\n";

    // Output data
    for (auto [name, pmap] : g.get_properties()) {
      if (!g.is_edge_key(pmap->key())) {
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

} // namespace bagl
