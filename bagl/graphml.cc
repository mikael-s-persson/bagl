
#include "bagl/graphml.h"

#include <any>
#include <iterator>
#include <limits>
#include <ranges>
#include <string>
#include <string_view>
#include <typeindex>
#include <unordered_map>

#include "bagl/dynamic_graph.h"
#include "tinyxml2.h"

namespace bagl {

#if 1

namespace {

// The kinds of keys. Not all of these are supported
enum class key_kind { graph_key, node_key, edge_key, hyperedge_key, port_key, endpoint_key, all_key, graphml_key };

template <typename HandlePropFunc>
void get_graphs(const tinyxml2::XMLElement* top, std::size_t desired_idx /* or max for all */, bool is_root,
                std::vector<const tinyxml2::XMLElement*>& result, const HandlePropFunc& handle_graph_property) {
  std::size_t current_idx = 0;
  bool is_first = is_root;
  for (const tinyxml2::XMLElement* n = top->FirstChildElement("graph"); n != nullptr;
       n = n->NextSiblingElement("graph")) {
    if (current_idx == desired_idx || desired_idx == all_graphs_desired) {
      result.push_back(n);
      if (is_first) {
        is_first = false;
        for (const tinyxml2::XMLElement* data_elem = n->FirstChildElement("data"); data_elem != nullptr;
             data_elem = data_elem->NextSiblingElement("data")) {
          std::string key = data_elem->Attribute("key");
          std::string value = data_elem->GetText();
          handle_graph_property(key, value);
        }
      }

      get_graphs(n, all_graphs_desired, false, result, handle_graph_property);
      if (desired_idx != all_graphs_desired) {
        break;
      }
    }
    ++current_idx;
  }
}

}  // namespace

void read_graphml(std::istream& in, dynamic_graph_mutator& g, std::size_t desired_idx) {
  tinyxml2::XMLDocument doc;

  const std::string in_content{std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
  doc.Parse(in_content.data(), in_content.size());

  const tinyxml2::XMLElement* doc_graphml = doc.FirstChildElement("graphml");

  std::unordered_map<std::string, key_kind> keys;
  std::unordered_map<std::string, std::string> key_name;
  std::unordered_map<std::string, std::string> key_type;
  std::unordered_map<std::string, std::string> key_default;

  for (const tinyxml2::XMLElement* doc_key = doc_graphml->FirstChildElement("key"); doc_key != nullptr;
       doc_key = doc_key->NextSiblingElement("key")) {
    std::string attr_id = doc_key->Attribute("id");
    std::string attr_for = doc_key->Attribute("for");
    std::string attr_name = doc_key->Attribute("attr.name");
    std::string attr_type = doc_key->Attribute("attr.type");
    key_kind kind = key_kind::all_key;
    if (attr_for == "graph") {
      kind = key_kind::graph_key;
    } else if (attr_for == "node") {
      kind = key_kind::node_key;
    } else if (attr_for == "edge") {
      kind = key_kind::edge_key;
    } else if (attr_for == "hyperedge") {
      kind = key_kind::hyperedge_key;
    } else if (attr_for == "port") {
      kind = key_kind::port_key;
    } else if (attr_for == "endpoint") {
      kind = key_kind::endpoint_key;
    } else if (attr_for == "all") {
      kind = key_kind::all_key;
    } else if (attr_for == "graphml") {
      kind = key_kind::graphml_key;
    } else {
      throw parse_error("Attribute for is not valid: " + attr_for);
    }
    keys[attr_id] = kind;
    key_name[attr_id] = attr_name;
    key_type[attr_id] = attr_type;
    if (const tinyxml2::XMLElement* default_xmlval = doc_key->FirstChildElement("default"); default_xmlval != nullptr) {
      key_default[attr_id] = default_xmlval->GetText();
    }
  }

  for (const auto& [key_id, default_val] : key_default) {
    if (keys[key_id] == key_kind::graph_key) {
      g.set_graph_property(key_name[key_id], default_val, key_type[key_id]);
    }
  }
  std::vector<const tinyxml2::XMLElement*> graphs;
  get_graphs(doc_graphml, desired_idx, true, graphs, [&](const std::string& key_id, const std::string& value) {
    g.set_graph_property(key_name[key_id], value, key_type[key_id]);
  });

  // Handle potentially new vertices.
  std::unordered_map<std::string, std::any> vertex;
  const auto handle_vertex = [&](const std::string& v) {
    if (vertex.find(v) == vertex.end()) {
      vertex[v] = g.do_add_vertex();
      for (const auto& [key_id, default_val] : key_default) {
        if (keys[key_id] == key_kind::node_key) {
          g.set_vertex_property(key_name[key_id], vertex[v], default_val, key_type[key_id]);
        }
      }
    }
  };

  for (const tinyxml2::XMLElement* gr : graphs) {
    // Search for nodes / vertices
    for (const tinyxml2::XMLElement* node = gr->FirstChildElement("node"); node != nullptr;
         node = node->NextSiblingElement("node")) {
      std::string id = node->Attribute("id");
      handle_vertex(id);
      for (const tinyxml2::XMLElement* data_elem = node->FirstChildElement("data"); data_elem != nullptr;
           data_elem = data_elem->NextSiblingElement("data")) {
        std::string key_id = data_elem->Attribute("key");
        std::string value = data_elem->GetText();
        g.set_vertex_property(key_name[key_id], vertex[id], value, key_type[key_id]);
      }
    }
  }

  std::vector<std::any> edge;

  for (const tinyxml2::XMLElement* gr : graphs) {
    const char* gr_edgedefault = gr->Attribute("edgedefault");
    bool default_directed = true;
    if (gr_edgedefault != nullptr && gr_edgedefault != std::string_view("directed")) {
      default_directed = false;
    }
    // Search for edges
    for (const tinyxml2::XMLElement* e = gr->FirstChildElement("edge"); e != nullptr;
         e = e->NextSiblingElement("edge")) {
      std::string source = e->Attribute("source");
      std::string target = e->Attribute("target");
      const char* local_directed = e->Attribute("directed");
      bool is_directed = (local_directed == nullptr ? default_directed : local_directed == std::string_view{"true"});
      if (is_directed != g.is_directed()) {
        if (is_directed) {
          throw directed_graph_error();
        } else {
          throw undirected_graph_error();
        }
      }

      handle_vertex(source);
      handle_vertex(target);

      std::any source_desc = vertex[source];
      std::any target_desc = vertex[target];

      auto [e_desc, e_added] = g.do_add_edge(source_desc, target_desc);
      if (!e_added) {
        throw bad_parallel_edge(source, target);
      }

      edge.push_back(e_desc);

      for (const auto& [key_id, default_val] : key_default) {
        if (keys[key_id] == key_kind::edge_key) {
          g.set_edge_property(key_name[key_id], e_desc, default_val, key_type[key_id]);
        }
      }

      for (const tinyxml2::XMLElement* data_elem = e->FirstChildElement("data"); data_elem != nullptr;
           data_elem = data_elem->NextSiblingElement("data")) {
        std::string key_id = data_elem->Attribute("key");
        std::string value = data_elem->GetText();
        g.set_edge_property(key_name[key_id], e_desc, value, key_type[key_id]);
      }
    }
  }
}

void write_graphml(std::ostream& out, const dynamic_graph_observer& g, bool ordered_vertices) {
  tinyxml2::XMLPrinter doc;

  doc.PushDeclaration(R"(version="1.0" encoding="UTF-8")");

  doc.OpenElement("graphml");
  doc.PushAttribute("xmlns", "http://graphml.graphdrawing.org/xmlns");
  doc.PushAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
  doc.PushAttribute("xsi:schemaLocation",
                    "http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd");

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
    const char* key_class_str = "graph";
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
    doc.OpenElement("key");
    doc.PushAttribute("id", key_id.c_str());
    doc.PushAttribute("for", key_class_str);
    doc.PushAttribute("attr.name", name.c_str());
    doc.PushAttribute("attr.type", type_name.c_str());
    doc.CloseElement();  // key
  }

  doc.OpenElement("graph");
  doc.PushAttribute("id", "G");
  doc.PushAttribute("edgedefault", g.is_directed() ? "directed" : "undirected");
  doc.PushAttribute("parse.nodeids", ordered_vertices ? "canonical" : "free");
  doc.PushAttribute("parse.edgeids", "canonical");
  doc.PushAttribute("parse.order", "nodesfirst");

  // Output graph data
  for (auto [name, pmap] : g.get_properties()) {
    if (!g.is_graph_key(pmap->key())) {
      continue;
    }
    doc.OpenElement("data");
    doc.PushAttribute("key", graph_key_ids[name].c_str());
    doc.PushText(pmap->get_string(g.get_graph_key()).c_str());
    doc.CloseElement();  // data
  }

  for (auto v : g.get_vertices()) {
    doc.OpenElement("node");
    doc.PushAttribute("id", ("n" + std::to_string(g.get_index_of(v))).c_str());
    // Output data
    for (auto [name, pmap] : g.get_properties()) {
      if (!g.is_vertex_key(pmap->key())) {
        continue;
      }
      doc.OpenElement("data");
      doc.PushAttribute("key", vertex_key_ids[name].c_str());
      doc.PushText(pmap->get_string(v).c_str());
      doc.CloseElement();  // data
    }
    doc.CloseElement();  // node
  }

  for (auto e : g.get_edges()) {
    doc.OpenElement("edge");
    doc.PushAttribute("id", ("e" + std::to_string(g.get_index_of(e))).c_str());
    doc.PushAttribute("source", ("n" + std::to_string(g.get_index_of(g.get_source(e)))).c_str());
    doc.PushAttribute("target", ("n" + std::to_string(g.get_index_of(g.get_target(e)))).c_str());
    // Output data
    for (auto [name, pmap] : g.get_properties()) {
      if (!g.is_edge_key(pmap->key())) {
        continue;
      }
      doc.OpenElement("data");
      doc.PushAttribute("key", edge_key_ids[name].c_str());
      doc.PushText(pmap->get_string(e).c_str());
      doc.CloseElement();  // data
    }
    doc.CloseElement();  // edge
  }

  doc.CloseElement();  // graph

  doc.CloseElement();  // graphml

  out << doc.CStr();
}

#else  // if 1

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
    out << "   <data key=\"" << graph_key_ids[name] << "\">" << encode_char_entities(pmap->get_string(&g))
        << "</data>\n";
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
    out << "    <edge id=\"e" << g.get_index_of(e) << "\" source=\"n" << g.get_index_of(g.get_source(e))
        << "\" target=\"n" << g.get_index_of(g.get_target(e)) << "\">\n";

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

#endif

}  // namespace bagl
