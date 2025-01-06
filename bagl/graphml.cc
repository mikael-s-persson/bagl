
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

  using any_vertex_desc = dynamic_graph_mutator::vertex_descriptor;
  using any_edge_desc = dynamic_graph_mutator::edge_descriptor;

  // Handle potentially new vertices.
  std::unordered_map<std::string, any_vertex_desc> vertex;
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

  std::vector<any_edge_desc> edge;

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

      any_vertex_desc source_desc = vertex[source];
      any_vertex_desc target_desc = vertex[target];

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
      assert(false);
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
    doc.PushAttribute("id", ("n" + std::to_string(g.get_index_of_vertex(v))).c_str());
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
    doc.PushAttribute("id", ("e" + std::to_string(g.get_index_of_edge(e))).c_str());
    doc.PushAttribute("source", ("n" + std::to_string(g.get_index_of_vertex(g.get_source(e)))).c_str());
    doc.PushAttribute("target", ("n" + std::to_string(g.get_index_of_vertex(g.get_target(e)))).c_str());
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

}  // namespace bagl
