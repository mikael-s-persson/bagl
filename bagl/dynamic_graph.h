// Copyright 2024 Mikael Persson

#ifndef BAGL_BAGL_DYNAMIC_GRAPH_H_
#define BAGL_BAGL_DYNAMIC_GRAPH_H_

#include <any>
#include <cstddef>
#include <ranges>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>

#include "bagl/any_iterator.h"
#include "bagl/dynamic_property_map.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/property_map.h"
#include "bagl/subgraph.h"

namespace bagl {

class dynamic_graph_mutator {
 public:
  dynamic_graph_mutator() = default;
  dynamic_graph_mutator(const dynamic_graph_mutator&) = delete;
  dynamic_graph_mutator(dynamic_graph_mutator&&) = delete;
  dynamic_graph_mutator& operator=(const dynamic_graph_mutator&) = delete;
  dynamic_graph_mutator& operator=(dynamic_graph_mutator&&) = delete;
  virtual ~dynamic_graph_mutator() = default;

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
class dynamic_graph_mutator_wrapper : public dynamic_graph_mutator {
  using vertex_descriptor = graph_vertex_descriptor_t<Graph>;
  using edge_descriptor = graph_edge_descriptor_t<Graph>;

 public:
  dynamic_graph_mutator_wrapper(Graph& g, dynamic_properties& dp) : g_(g), dp_(dp) {}

  [[nodiscard]] bool is_directed() const override { return is_directed_graph_v<Graph>; }

  std::any do_add_vertex() override { return std::any(add_vertex(g_)); }

  std::pair<std::any, bool> do_add_edge(std::any source, std::any target) override {
    auto [e, added] = add_edge(std::any_cast<vertex_descriptor>(source), std::any_cast<vertex_descriptor>(target), g_);
    return std::make_pair(std::any(e), added);
  }

  void set_graph_property(const std::string& name, const std::string& value, const std::string& value_type) override {
    bool type_found = put_property(name, &g_, value, value_type);
    if (!type_found) {
      throw std::bad_any_cast();
    }
  }

  void set_vertex_property(const std::string& name, std::any vertex, const std::string& value,
                           const std::string& value_type) override {
    bool type_found = put_property(name, std::any_cast<vertex_descriptor>(vertex), value, value_type);
    if (!type_found) {
      throw std::bad_any_cast();
    }
  }

  void set_edge_property(const std::string& name, std::any edge, const std::string& value,
                         const std::string& value_type) override {
    bool type_found = put_property(name, std::any_cast<edge_descriptor>(edge), value, value_type);
    if (!type_found) {
      throw std::bad_any_cast();
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

// This is a type-erasure interface for a graph that provides some of the common interfaces as
// virtual member functions (rather than statically dispatched free functions).
// Vertices and edges are represented as std::any.
// Properties are dispatched through a `dynamic_properties` structure.
// NOTE: Member functions have the "get_" prefix to avoid name-clashing with free functions.
// NOTE: It should be obvious that this type-erased interface will be much slower and unsafe
// than a statically-known graph type.
class dynamic_graph_observer {
 public:
  dynamic_graph_observer() = default;
  dynamic_graph_observer(const dynamic_graph_observer&) = delete;
  dynamic_graph_observer(dynamic_graph_observer&&) = delete;
  dynamic_graph_observer& operator=(const dynamic_graph_observer&) = delete;
  dynamic_graph_observer& operator=(dynamic_graph_observer&&) = delete;
  virtual ~dynamic_graph_observer() = default;

  [[nodiscard]] virtual bool is_directed() const = 0;

  enum class key_type {
    graph,
    vertex,
    edge,
    unknown,
  };
  [[nodiscard]] virtual key_type classify_key(const std::type_info& key_info) const = 0;
  [[nodiscard]] bool is_graph_key(const std::type_info& key_info) const {
    return classify_key(key_info) == key_type::graph;
  }
  [[nodiscard]] bool is_vertex_key(const std::type_info& key_info) const {
    return classify_key(key_info) == key_type::vertex;
  }
  [[nodiscard]] bool is_edge_key(const std::type_info& key_info) const {
    return classify_key(key_info) == key_type::edge;
  }

  [[nodiscard]] virtual const dynamic_properties& get_properties() const = 0;

  [[nodiscard]] virtual std::any get_graph_key() const = 0;

  // Subgraph
  [[nodiscard]] virtual bool is_root() const { return true; }
  [[nodiscard]] virtual const dynamic_graph_observer& get_root() const { return *this; }
  [[nodiscard]] virtual std::size_t get_num_children() const { return 0; }
  [[nodiscard]] virtual any_range<const dynamic_graph_observer&> get_children() const {
    return make_any_range_to<const dynamic_graph_observer&>(std::ranges::empty_view<const dynamic_graph_observer>());
  }

  // VertexListGraph
  [[nodiscard]] virtual std::size_t get_num_vertices() const = 0;
  [[nodiscard]] virtual any_range<std::any> get_vertices() const = 0;

  // EdgeListGraph
  [[nodiscard]] virtual std::size_t get_num_edges() const = 0;
  [[nodiscard]] virtual any_range<std::any> get_edges() const = 0;

  // IncidenceGraph
  [[nodiscard]] virtual std::size_t get_out_degree(const std::any& u) const = 0;
  [[nodiscard]] virtual any_range<std::any> get_out_edges(const std::any& u) const = 0;

  [[nodiscard]] virtual std::any get_source(const std::any& e) const = 0;
  [[nodiscard]] virtual std::any get_target(const std::any& e) const = 0;

  // AdjacencyMatrix
  [[nodiscard]] virtual std::pair<std::any, bool> get_edge(const std::any& u, const std::any& v) const = 0;

  // ReadableVertexIndexMap / ReadableEdgeIndexMap
  [[nodiscard]] virtual std::size_t get_index_of(const std::any& u_or_e) const = 0;
};

// This type-erases a given graph type.
template <concepts::Graph Graph>
class dynamic_graph_observer_wrapper : public dynamic_graph_observer {
 public:
  dynamic_graph_observer_wrapper(const Graph& g, const dynamic_properties& dp) : g_(g), dp_(dp) {
    if constexpr (!has_property_map_v<Graph, edge_index_t> && concepts::EdgeListGraph<Graph>) {
      // Build edge index.
      eindex_.reserve(num_edges(g_));
      std::size_t i = 0;
      for (auto e : edges(g_)) {
        eindex_.emplace(e, i++);
      }
    }
    if constexpr (!has_property_map_v<Graph, vertex_index_t> && concepts::VertexListGraph<Graph>) {
      // Build vertex index.
      vindex_.reserve(num_vertices(g_));
      std::size_t i = 0;
      for (auto v : vertices(g_)) {
        vindex_.emplace(v, i++);
      }
    }
  }

  [[nodiscard]] bool is_directed() const override { return is_directed_graph_v<Graph>; }

  using vertex_descriptor = graph_vertex_descriptor_t<Graph>;
  using edge_descriptor = graph_edge_descriptor_t<Graph>;
  using key_type = dynamic_graph_observer::key_type;

  [[nodiscard]] key_type classify_key(const std::type_info& key_info) const override {
    if (key_info == typeid(const Graph*) || key_info == typeid(Graph*)) {
      return key_type::graph;
    }
    if (key_info == typeid(vertex_descriptor)) {
      return key_type::vertex;
    }
    if (key_info == typeid(edge_descriptor)) {
      return key_type::edge;
    }
    return key_type::unknown;
  }

  [[nodiscard]] const dynamic_properties& get_properties() const override { return dp_; }

  [[nodiscard]] std::any get_graph_key() const override { return &g_; }

  // VertexListGraph
  [[nodiscard]] std::size_t get_num_vertices() const override {
    if constexpr (concepts::VertexListGraph<Graph>) {
      return num_vertices(g_);
    } else {
      return 0;
    }
  }
  [[nodiscard]] any_range<std::any> get_vertices() const override {
    if constexpr (concepts::VertexListGraph<Graph>) {
      return make_any_range_to_any(vertices(g_));
    } else {
      return make_any_range_to_any(std::ranges::empty_view<vertex_descriptor>());
    }
  }

  // EdgeListGraph
  [[nodiscard]] std::size_t get_num_edges() const override {
    if constexpr (concepts::EdgeListGraph<Graph>) {
      return num_edges(g_);
    } else {
      return 0;
    }
  }
  [[nodiscard]] any_range<std::any> get_edges() const override {
    if constexpr (concepts::EdgeListGraph<Graph>) {
      return make_any_range_to_any(edges(g_));
    } else {
      return make_any_range_to_any(std::ranges::empty_view<edge_descriptor>());
    }
  }

  // IncidenceGraph
  [[nodiscard]] std::size_t get_out_degree(const std::any& u) const override {
    if constexpr (concepts::IncidenceGraph<Graph>) {
      return out_degree(std::any_cast<vertex_descriptor>(u), g_);
    } else {
      return 0;
    }
  }
  [[nodiscard]] any_range<std::any> get_out_edges(const std::any& u) const override {
    if constexpr (concepts::IncidenceGraph<Graph>) {
      return make_any_range_to_any(out_edges(std::any_cast<vertex_descriptor>(u), g_));
    } else {
      return make_any_range_to_any(std::ranges::empty_view<edge_descriptor>());
    }
  }

  [[nodiscard]] std::any get_source(const std::any& e) const override {
    if constexpr (concepts::IncidenceGraph<Graph>) {
      return source(std::any_cast<edge_descriptor>(e), g_);
    } else {
      return graph_traits<Graph>::null_vertex();
    }
  }
  [[nodiscard]] std::any get_target(const std::any& e) const override {
    if constexpr (concepts::IncidenceGraph<Graph>) {
      return target(std::any_cast<edge_descriptor>(e), g_);
    } else {
      return graph_traits<Graph>::null_vertex();
    }
  }

  // AdjacencyMatrix
  [[nodiscard]] std::pair<std::any, bool> get_edge(const std::any& u, const std::any& v) const override {
    if constexpr (concepts::AdjacencyMatrix<Graph>) {
      return edge(std::any_cast<vertex_descriptor>(u), std::any_cast<vertex_descriptor>(v), g_);
    } else {
      return {edge_descriptor{}, false};
    }
  }

  // ReadableVertexIndexMap / ReadableEdgeIndexMap
  [[nodiscard]] std::size_t get_index_of(const std::any& u_or_e) const override {
    if (u_or_e.type() == typeid(vertex_descriptor)) {
      if constexpr (has_property_map_v<Graph, vertex_index_t>) {
        return get(vertex_index, g_, std::any_cast<vertex_descriptor>(u_or_e));
      } else {
        auto it = vindex_.find(std::any_cast<vertex_descriptor>(u_or_e));
        return it != vindex_.end() ? it->second : 0;
      }
    }
    if (u_or_e.type() == typeid(edge_descriptor)) {
      if constexpr (has_property_map_v<Graph, edge_index_t>) {
        return get(edge_index, g_, std::any_cast<edge_descriptor>(u_or_e));
      } else {
        auto it = eindex_.find(std::any_cast<edge_descriptor>(u_or_e));
        return it != eindex_.end() ? it->second : 0;
      }
    }
    // Maybe graph key.
    return 0;
  }

 protected:
  const Graph& g_;
  const dynamic_properties& dp_;
  std::map<vertex_descriptor, std::size_t> vindex_;
  std::map<edge_descriptor, std::size_t> eindex_;
};

template <concepts::Graph Graph>
class dynamic_graph_observer_wrapper<subgraph<Graph>> : public dynamic_graph_observer_wrapper<Graph> {
  using Base = dynamic_graph_observer_wrapper<Graph>;

 public:
  dynamic_graph_observer_wrapper(const subgraph<Graph>& sub_g, const dynamic_properties& dp)
      : Base(sub_g.root().underlying(), dp), sub_g_(sub_g), parent_(*this) {
    children_.reserve(sub_g.num_children());
    for (const subgraph<Graph>& child : sub_g.children()) {
      children_.emplace_back(std::make_unique<dynamic_graph_observer_wrapper<subgraph<Graph>>>(child, dp));
    }
  }
  dynamic_graph_observer_wrapper(const subgraph<Graph>& sub_g, const dynamic_properties& dp,
                                 const dynamic_graph_observer& parent)
      : Base(sub_g.root().underlying(), dp), sub_g_(sub_g), parent_(parent) {
    children_.reserve(sub_g.num_children());
    for (const subgraph<Graph>& child : sub_g.children()) {
      children_.emplace_back(std::make_unique<dynamic_graph_observer_wrapper<subgraph<Graph>>>(child, dp));
    }
  }

  using vertex_descriptor = graph_vertex_descriptor_t<Graph>;
  using edge_descriptor = graph_edge_descriptor_t<Graph>;
  using key_type = dynamic_graph_observer::key_type;

  // Subgraph
  [[nodiscard]] bool is_root() const override { return sub_g_.is_root(); }
  [[nodiscard]] const dynamic_graph_observer& get_root() const override {
    if (parent_ == *this) {
      return *this;
    }
    return parent_;
  }
  [[nodiscard]] any_range<const dynamic_graph_observer&> get_children() const override {
    return make_any_range_to<const dynamic_graph_observer&>(std::ranges::ref_view(children_) | std::views::transform([](const auto& child_ptr) {
                                   return std::as_const(*child_ptr);
                                 }));
  }

  [[nodiscard]] key_type classify_key(const std::type_info& key_info) const override {
    if (key_info == typeid(const subgraph<Graph>*) || key_info == typeid(subgraph<Graph>*)) {
      return key_type::graph;
    }
    if (key_info == typeid(vertex_descriptor)) {
      return key_type::vertex;
    }
    if (key_info == typeid(edge_descriptor)) {
      return key_type::edge;
    }
    return key_type::unknown;
  }

  [[nodiscard]] std::any get_graph_key() const override { return &sub_g_; }

  // VertexListGraph
  [[nodiscard]] std::size_t get_num_vertices() const override {
    if constexpr (concepts::VertexListGraph<Graph>) {
      return num_vertices(sub_g_);
    } else {
      return 0;
    }
  }
  [[nodiscard]] any_range<std::any> get_vertices() const override {
    if constexpr (concepts::VertexListGraph<Graph>) {
      return make_any_range_to_any(vertices(sub_g_) | std::views::transform([this](const vertex_descriptor& v) {
                                     return sub_g_.local_to_global(v);
                                   }));
    } else {
      return make_any_range_to_any(std::ranges::empty_view<vertex_descriptor>());
    }
  }

  // EdgeListGraph
  [[nodiscard]] std::size_t get_num_edges() const override {
    if constexpr (concepts::EdgeListGraph<Graph>) {
      return num_edges(sub_g_);
    } else {
      return 0;
    }
  }
  [[nodiscard]] any_range<std::any> get_edges() const override {
    if constexpr (concepts::EdgeListGraph<Graph>) {
      return make_any_range_to_any(edges(sub_g_) | std::views::transform([this](const edge_descriptor& e) {
                                     return sub_g_.local_to_global(e);
                                   }));
    } else {
      return make_any_range_to_any(std::ranges::empty_view<edge_descriptor>());
    }
  }

  // IncidenceGraph
  [[nodiscard]] std::size_t get_out_degree(const std::any& u) const override {
    if constexpr (concepts::IncidenceGraph<Graph>) {
      return out_degree(sub_g_.global_to_local(std::any_cast<vertex_descriptor>(u)), sub_g_);
    } else {
      return 0;
    }
  }
  [[nodiscard]] any_range<std::any> get_out_edges(const std::any& u) const override {
    if constexpr (concepts::IncidenceGraph<Graph>) {
      return make_any_range_to_any(
          out_edges(sub_g_.global_to_local(std::any_cast<vertex_descriptor>(u)), sub_g_) |
          std::views::transform([this](const edge_descriptor& e) { return sub_g_.local_to_global(e); }));
    } else {
      return make_any_range_to_any(std::ranges::empty_view<edge_descriptor>());
    }
  }

  // AdjacencyMatrix
  [[nodiscard]] std::pair<std::any, bool> get_edge(const std::any& u, const std::any& v) const override {
    if constexpr (concepts::AdjacencyMatrix<Graph>) {
      auto [e_global, e_found] =
          edge(std::any_cast<vertex_descriptor>(u), std::any_cast<vertex_descriptor>(v), sub_g_.root().underlying());
      if (e_found) {
        return {e_global, sub_g_.find_edge(e_global).second};
      }
    }
    return {edge_descriptor{}, false};
  }

 protected:
  const subgraph<Graph>& sub_g_;
  std::vector<std::unique_ptr<dynamic_graph_observer>> children_;
  const dynamic_graph_observer& parent_;
};

}  // namespace bagl

#endif  // BAGL_BAGL_DYNAMIC_GRAPH_H_
