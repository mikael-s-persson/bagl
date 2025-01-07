// Copyright 2024 Mikael Persson

#ifndef BAGL_BAGL_DYNAMIC_GRAPH_H_
#define BAGL_BAGL_DYNAMIC_GRAPH_H_

#include <any>
#include <cstddef>
#include <iostream>
#include <ranges>
#include <type_traits>
#include <typeinfo>
#include <unordered_map>

#include "bagl/adjacency_range.h"
#include "bagl/any_iterator.h"
#include "bagl/dynamic_property_map.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/property_map.h"
#include "bagl/subgraph.h"

namespace bagl {

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

  using directed_category = bidirectional_tag;
  using edge_parallel_category = allow_parallel_edge_tag;
  struct traversal_category : bidirectional_graph_tag,
                              adjacency_graph_tag,
                              vertex_list_graph_tag,
                              edge_list_graph_tag {};

  struct vertex_descriptor {
    std::any base;
    explicit vertex_descriptor(const std::any& b) : base(b) {}
    vertex_descriptor() : vertex_descriptor(std::any{}) {}
    operator std::any&() { return base; }              // NOLINT
    operator const std::any&() const { return base; }  // NOLINT
  };

  struct edge_descriptor {
    std::any base;
    explicit edge_descriptor(const std::any& b) : base(b) {}
    edge_descriptor() : edge_descriptor(std::any{}) {}
    operator std::any&() { return base; }              // NOLINT
    operator const std::any&() const { return base; }  // NOLINT
  };

  using graph_property_type = std::any;
  using vertex_property_type = std::any;
  using edge_property_type = std::any;
  using graph_bundled = std::any;
  using vertex_bundled = std::any;
  using edge_bundled = std::any;

  using edge_range = any_range<edge_descriptor>;
  using out_edge_range = edge_range;
  using in_edge_range = edge_range;

  using vertex_range = any_range<vertex_descriptor>;

  [[nodiscard]] virtual bool is_directed() const = 0;
  [[nodiscard]] virtual bool allows_parallel_edges() const = 0;
  [[nodiscard]] virtual bool is_incidence_graph() const = 0;
  [[nodiscard]] virtual bool is_bidirectional_graph() const = 0;
  [[nodiscard]] virtual bool is_vertex_list_graph() const = 0;
  [[nodiscard]] virtual bool is_edge_list_graph() const = 0;
  [[nodiscard]] virtual bool is_adjacency_matrix() const = 0;

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

  [[nodiscard]] virtual graph_property_type get_graph_property() const = 0;
  [[nodiscard]] virtual vertex_property_type get_vertex_property(const vertex_descriptor& u) const = 0;
  [[nodiscard]] virtual edge_property_type get_edge_property(const edge_descriptor& e) const = 0;

  [[nodiscard]] virtual graph_bundled get_graph_bundled() const = 0;
  [[nodiscard]] virtual vertex_bundled get_vertex_bundled(const vertex_descriptor& u) const = 0;
  [[nodiscard]] virtual edge_bundled get_edge_bundled(const edge_descriptor& e) const = 0;

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
  [[nodiscard]] virtual vertex_range get_vertices() const = 0;

  // EdgeListGraph
  [[nodiscard]] virtual std::size_t get_num_edges() const = 0;
  [[nodiscard]] virtual edge_range get_edges() const = 0;

  // IncidenceGraph
  [[nodiscard]] virtual std::size_t get_out_degree(const vertex_descriptor& u) const = 0;
  [[nodiscard]] virtual edge_range get_out_edges(const vertex_descriptor& u) const = 0;

  // BidirectionalGraph
  [[nodiscard]] virtual std::size_t get_in_degree(const vertex_descriptor& u) const = 0;
  [[nodiscard]] virtual edge_range get_in_edges(const vertex_descriptor& u) const = 0;
  [[nodiscard]] virtual std::size_t get_degree(const vertex_descriptor& u) const = 0;

  [[nodiscard]] virtual vertex_descriptor get_source(const edge_descriptor& e) const = 0;
  [[nodiscard]] virtual vertex_descriptor get_target(const edge_descriptor& e) const = 0;

  // AdjacencyMatrix
  [[nodiscard]] virtual std::pair<edge_descriptor, bool> get_edge(const vertex_descriptor& u,
                                                                  const vertex_descriptor& v) const = 0;

  // ReadableVertexIndexMap / ReadableEdgeIndexMap
  [[nodiscard]] virtual std::size_t get_index_of_edge(const edge_descriptor& e) const = 0;
  [[nodiscard]] virtual std::size_t get_index_of_vertex(const vertex_descriptor& u) const = 0;

  // Equality of vertices or edges
  [[nodiscard]] virtual bool are_vertices_equal(const vertex_descriptor& u, const vertex_descriptor& v) const = 0;
  [[nodiscard]] virtual bool are_edges_equal(const edge_descriptor& e, const edge_descriptor& f) const = 0;
};

class dynamic_graph_mutator : public dynamic_graph_observer {
 public:
  using vertex_descriptor = dynamic_graph_observer::vertex_descriptor;
  using edge_descriptor = dynamic_graph_observer::edge_descriptor;

  dynamic_graph_mutator() = default;
  ~dynamic_graph_mutator() override = default;

  [[nodiscard]] virtual dynamic_properties& get_mutable_properties() = 0;

  virtual vertex_descriptor do_add_vertex() = 0;
  virtual void do_clear_vertex(const vertex_descriptor& u) = 0;
  virtual void do_remove_vertex(const vertex_descriptor& u) = 0;

  virtual std::pair<edge_descriptor, bool> do_add_edge(const vertex_descriptor& source,
                                                       const vertex_descriptor& target) = 0;
  virtual void do_remove_edge(const edge_descriptor& e) = 0;

  virtual void do_remove_out_edge_if(const vertex_descriptor& u, std::function<bool(const edge_descriptor&)> pred) = 0;
  virtual void do_remove_in_edge_if(const vertex_descriptor& v, std::function<bool(const edge_descriptor&)> pred) = 0;
  virtual void do_remove_edge_if(std::function<bool(const edge_descriptor&)> pred) = 0;

  virtual void set_graph_property(const std::string& name, const std::string& value, const std::string& value_type) = 0;

  virtual void set_vertex_property(const std::string& name, const vertex_descriptor& vertex, const std::string& value,
                                   const std::string& value_type) = 0;

  virtual void set_edge_property(const std::string& name, const edge_descriptor& edge, const std::string& value,
                                 const std::string& value_type) = 0;
};

// This type-erases a given graph type.
template <concepts::Graph Graph, typename Base = dynamic_graph_observer>
class dynamic_graph_observer_wrapper : public Base {
  using DynProps = std::conditional_t<std::is_const_v<Graph>, const dynamic_properties, dynamic_properties>;

 public:
  using vertex_descriptor = typename Base::vertex_descriptor;
  using edge_descriptor = typename Base::edge_descriptor;

  using graph_property_type = typename Base::graph_property_type;
  using vertex_property_type = typename Base::vertex_property_type;
  using edge_property_type = typename Base::edge_property_type;
  using graph_bundled = typename Base::graph_bundled;
  using vertex_bundled = typename Base::vertex_bundled;
  using edge_bundled = typename Base::edge_bundled;

  using edge_range = typename Base::edge_range;
  using out_edge_range = typename Base::out_edge_range;
  using in_edge_range = typename Base::in_edge_range;
  using vertex_range = typename Base::vertex_range;

  dynamic_graph_observer_wrapper(Graph& g, DynProps& dp) : g_(g), dp_(dp) {
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

  [[nodiscard]] std::add_const_t<Graph>& real_graph() const { return g_; }

  [[nodiscard]] bool is_directed() const override { return is_directed_graph_v<Graph>; }
  [[nodiscard]] bool allows_parallel_edges() const override { return allows_parallel_edges_v<Graph>; }
  [[nodiscard]] bool is_incidence_graph() const override { return is_incidence_graph_v<Graph>; }
  [[nodiscard]] bool is_bidirectional_graph() const override { return is_bidirectional_graph_v<Graph>; }
  [[nodiscard]] bool is_vertex_list_graph() const override { return is_vertex_list_graph_v<Graph>; }
  [[nodiscard]] bool is_edge_list_graph() const override { return is_edge_list_graph_v<Graph>; }
  [[nodiscard]] bool is_adjacency_matrix() const override { return is_adjacency_matrix_v<Graph>; }

  using real_vertex_descriptor = graph_vertex_descriptor_t<Graph>;
  using real_edge_descriptor = graph_edge_descriptor_t<Graph>;
  using key_type = dynamic_graph_observer::key_type;

  [[nodiscard]] key_type classify_key(const std::type_info& key_info) const override {
    if (key_info == typeid(Graph*)) {
      return key_type::graph;
    }
    if constexpr (std::is_const_v<Graph>) {
      if (key_info == typeid(std::remove_const_t<Graph>*)) {
        return key_type::graph;
      }
    } else {
      if (key_info == typeid(std::add_const_t<Graph>*)) {
        return key_type::graph;
      }
    }
    if (key_info == typeid(real_vertex_descriptor)) {
      return key_type::vertex;
    }
    if (key_info == typeid(real_edge_descriptor)) {
      return key_type::edge;
    }
    return key_type::unknown;
  }

  [[nodiscard]] const dynamic_properties& get_properties() const override { return dp_; }

  [[nodiscard]] graph_property_type get_graph_property() const override {
    return graph_property_type{get_property(g_, graph_all)};
  }
  [[nodiscard]] vertex_property_type get_vertex_property(const vertex_descriptor& u) const override {
    return vertex_property_type{get_property(g_, std::any_cast<real_vertex_descriptor>(u))};
  }
  [[nodiscard]] edge_property_type get_edge_property(const edge_descriptor& e) const override {
    return edge_property_type{get_property(g_, std::any_cast<real_edge_descriptor>(e))};
  }

  [[nodiscard]] graph_bundled get_graph_bundled() const override { return graph_bundled{get(g_, graph_bundle)}; }
  [[nodiscard]] vertex_bundled get_vertex_bundled(const vertex_descriptor& u) const override {
    return vertex_bundled{get(vertex_bundle, g_, std::any_cast<real_vertex_descriptor>(u))};
  }
  [[nodiscard]] edge_bundled get_edge_bundled(const edge_descriptor& e) const override {
    return edge_bundled{get(edge_bundle, g_, std::any_cast<real_edge_descriptor>(e))};
  }

  [[nodiscard]] std::any get_graph_key() const override { return &g_; }

  // VertexListGraph
  [[nodiscard]] std::size_t get_num_vertices() const override {
    if constexpr (concepts::VertexListGraph<Graph>) {
      return num_vertices(g_);
    } else {
      return 0;
    }
  }
  [[nodiscard]] vertex_range get_vertices() const override {
    if constexpr (concepts::VertexListGraph<Graph>) {
      return make_any_range_to<vertex_descriptor>(vertices(g_) | std::views::transform(real_to_any_vertex_desc()));
    } else {
      return make_any_range_to<vertex_descriptor>(std::ranges::empty_view<vertex_descriptor>());
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
  [[nodiscard]] edge_range get_edges() const override {
    if constexpr (concepts::EdgeListGraph<Graph>) {
      return make_any_range_to<edge_descriptor>(edges(g_) | std::views::transform(real_to_any_edge_desc()));
    } else {
      return make_any_range_to<edge_descriptor>(std::ranges::empty_view<edge_descriptor>());
    }
  }

  // IncidenceGraph
  [[nodiscard]] std::size_t get_out_degree(const vertex_descriptor& u) const override {
    if constexpr (concepts::IncidenceGraph<Graph>) {
      return out_degree(std::any_cast<real_vertex_descriptor>(u), g_);
    } else {
      return 0;
    }
  }
  [[nodiscard]] edge_range get_out_edges(const vertex_descriptor& u) const override {
    if constexpr (concepts::IncidenceGraph<Graph>) {
      return make_any_range_to<edge_descriptor>(out_edges(std::any_cast<real_vertex_descriptor>(u), g_) |
                                                std::views::transform(real_to_any_edge_desc()));
    } else {
      return make_any_range_to<edge_descriptor>(std::ranges::empty_view<edge_descriptor>());
    }
  }

  // BidirectionalGraph
  [[nodiscard]] std::size_t get_in_degree(const vertex_descriptor& u) const override {
    if constexpr (concepts::BidirectionalGraph<Graph>) {
      return in_degree(std::any_cast<real_vertex_descriptor>(u), g_);
    } else {
      return 0;
    }
  }
  [[nodiscard]] edge_range get_in_edges(const vertex_descriptor& u) const override {
    if constexpr (concepts::BidirectionalGraph<Graph>) {
      return make_any_range_to<edge_descriptor>(in_edges(std::any_cast<real_vertex_descriptor>(u), g_) |
                                                std::views::transform(real_to_any_edge_desc()));
    } else {
      return make_any_range_to<edge_descriptor>(std::ranges::empty_view<edge_descriptor>());
    }
  }
  [[nodiscard]] std::size_t get_degree(const vertex_descriptor& u) const override {
    if constexpr (concepts::BidirectionalGraph<Graph>) {
      return degree(std::any_cast<real_vertex_descriptor>(u), g_);
    } else {
      return 0;
    }
  }

  [[nodiscard]] vertex_descriptor get_source(const edge_descriptor& e) const override {
    if constexpr (concepts::IncidenceGraph<Graph>) {
      return vertex_descriptor{source(std::any_cast<real_edge_descriptor>(e), g_)};
    } else {
      return vertex_descriptor{graph_traits<Graph>::null_vertex()};
    }
  }
  [[nodiscard]] vertex_descriptor get_target(const edge_descriptor& e) const override {
    if constexpr (concepts::IncidenceGraph<Graph>) {
      return vertex_descriptor{target(std::any_cast<real_edge_descriptor>(e), g_)};
    } else {
      return vertex_descriptor{graph_traits<Graph>::null_vertex()};
    }
  }

  // AdjacencyMatrix
  [[nodiscard]] std::pair<edge_descriptor, bool> get_edge(const vertex_descriptor& u,
                                                          const vertex_descriptor& v) const override {
    if constexpr (concepts::AdjacencyMatrix<Graph>) {
      auto [e, e_exists] = edge(std::any_cast<real_vertex_descriptor>(u), std::any_cast<real_vertex_descriptor>(v), g_);
      return {edge_descriptor{e}, e_exists};
    } else {
      return {edge_descriptor{real_edge_descriptor{}}, false};
    }
  }

  // ReadableVertexIndexMap / ReadableEdgeIndexMap
  [[nodiscard]] std::size_t get_index_of_vertex(const vertex_descriptor& u) const override {
    if constexpr (has_property_map_v<Graph, vertex_index_t>) {
      return get(vertex_index, g_, std::any_cast<real_vertex_descriptor>(u));
    } else {
      auto it = vindex_.find(std::any_cast<real_vertex_descriptor>(u));
      return it != vindex_.end() ? it->second : 0;
    }
  }
  [[nodiscard]] std::size_t get_index_of_edge(const edge_descriptor& e) const override {
    if constexpr (has_property_map_v<Graph, edge_index_t>) {
      return get(edge_index, g_, std::any_cast<real_edge_descriptor>(e));
    } else {
      auto it = eindex_.find(std::any_cast<real_edge_descriptor>(e));
      return it != eindex_.end() ? it->second : 0;
    }
  }

  // Equality of vertices or edges
  [[nodiscard]] bool are_vertices_equal(const vertex_descriptor& u, const vertex_descriptor& v) const override {
    return std::any_cast<real_vertex_descriptor>(u) == std::any_cast<real_vertex_descriptor>(v);
  }
  [[nodiscard]] bool are_edges_equal(const edge_descriptor& e, const edge_descriptor& f) const override {
    return std::any_cast<real_edge_descriptor>(e) == std::any_cast<real_edge_descriptor>(f);
  }

 protected:
  Graph& g_;
  DynProps& dp_;
  std::unordered_map<real_vertex_descriptor, std::size_t, graph_descriptor_hash_t<real_vertex_descriptor>> vindex_;
  std::unordered_map<real_edge_descriptor, std::size_t, graph_descriptor_hash_t<real_edge_descriptor>> eindex_;

  [[nodiscard]] auto real_to_any_vertex_desc() const {
    return [](const real_vertex_descriptor& v) { return vertex_descriptor{v}; };
  }
  [[nodiscard]] auto real_to_any_edge_desc() const {
    return [](const real_edge_descriptor& e) { return edge_descriptor{e}; };
  }
};

template <typename SubGraph, concepts::Graph Graph, typename DynBase>
class dynamic_subgraph_observer_wrapper : public dynamic_graph_observer_wrapper<Graph, DynBase> {
  using Base = dynamic_graph_observer_wrapper<Graph, DynBase>;
  using DynProps = std::conditional_t<std::is_const_v<Graph>, const dynamic_properties, dynamic_properties>;
  using ParentObs = std::conditional_t<std::is_const_v<Graph>, const DynBase, DynBase>;

 public:
  using vertex_descriptor = typename Base::vertex_descriptor;
  using edge_descriptor = typename Base::edge_descriptor;
  using vertex_range = typename Base::vertex_range;
  using edge_range = typename Base::edge_range;

  dynamic_subgraph_observer_wrapper(SubGraph& sub_g, DynProps& dp)
      : Base(sub_g.root().underlying(), dp), sub_g_(sub_g), parent_(*this) {
    children_.reserve(sub_g.num_children());
    for (SubGraph& child : sub_g.children()) {
      children_.emplace_back(std::make_unique<dynamic_subgraph_observer_wrapper<SubGraph, Graph, DynBase>>(child, dp));
    }
  }
  dynamic_subgraph_observer_wrapper(SubGraph& sub_g, DynProps& dp, ParentObs& parent)
      : Base(sub_g.root().underlying(), dp), sub_g_(sub_g), parent_(parent) {
    children_.reserve(sub_g.num_children());
    for (SubGraph& child : sub_g.children()) {
      children_.emplace_back(std::make_unique<dynamic_subgraph_observer_wrapper<SubGraph, Graph, DynBase>>(child, dp));
    }
  }

  using real_vertex_descriptor = graph_vertex_descriptor_t<Graph>;
  using real_edge_descriptor = graph_edge_descriptor_t<Graph>;
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
    return make_any_range_to<const dynamic_graph_observer&>(
        std::ranges::ref_view(children_) |
        std::views::transform([](const auto& child_ptr) { return std::as_const(*child_ptr); }));
  }

  [[nodiscard]] key_type classify_key(const std::type_info& key_info) const override {
    if (key_info == typeid(SubGraph*)) {
      return key_type::graph;
    }
    if constexpr (std::is_const_v<SubGraph>) {
      if (key_info == typeid(std::remove_const_t<SubGraph>*)) {
        return key_type::graph;
      }
    } else {
      if (key_info == typeid(std::add_const_t<SubGraph>*)) {
        return key_type::graph;
      }
    }
    if (key_info == typeid(real_vertex_descriptor)) {
      return key_type::vertex;
    }
    if (key_info == typeid(real_edge_descriptor)) {
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
  [[nodiscard]] vertex_range get_vertices() const override {
    if constexpr (concepts::VertexListGraph<Graph>) {
      return make_any_range_to<vertex_descriptor>(vertices(sub_g_) | std::views::transform(real_to_any_vertex_desc()));
    } else {
      return make_any_range_to<vertex_descriptor>(std::ranges::empty_view<real_vertex_descriptor>());
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
  [[nodiscard]] edge_range get_edges() const override {
    if constexpr (concepts::EdgeListGraph<Graph>) {
      return make_any_range_to<edge_descriptor>(edges(sub_g_) | std::views::transform(real_to_any_edge_desc()));
    } else {
      return make_any_range_to<edge_descriptor>(std::ranges::empty_view<real_edge_descriptor>());
    }
  }

  // IncidenceGraph
  [[nodiscard]] std::size_t get_out_degree(const vertex_descriptor& u) const override {
    if constexpr (concepts::IncidenceGraph<Graph>) {
      return out_degree(sub_g_.global_to_local(std::any_cast<real_vertex_descriptor>(u)), sub_g_);
    } else {
      return 0;
    }
  }
  [[nodiscard]] edge_range get_out_edges(const vertex_descriptor& u) const override {
    if constexpr (concepts::IncidenceGraph<Graph>) {
      return make_any_range_to<edge_descriptor>(
          out_edges(sub_g_.global_to_local(std::any_cast<real_vertex_descriptor>(u)), sub_g_) |
          std::views::transform(real_to_any_edge_desc()));
    } else {
      return make_any_range_to<edge_descriptor>(std::ranges::empty_view<real_edge_descriptor>());
    }
  }

  // IncidenceGraph
  [[nodiscard]] std::size_t get_in_degree(const vertex_descriptor& u) const override {
    if constexpr (concepts::BidirectionalGraph<Graph>) {
      return in_degree(sub_g_.global_to_local(std::any_cast<real_vertex_descriptor>(u)), sub_g_);
    } else {
      return 0;
    }
  }
  [[nodiscard]] edge_range get_in_edges(const vertex_descriptor& u) const override {
    if constexpr (concepts::BidirectionalGraph<Graph>) {
      return make_any_range_to<edge_descriptor>(
          in_edges(sub_g_.global_to_local(std::any_cast<real_vertex_descriptor>(u)), sub_g_) |
          std::views::transform(real_to_any_edge_desc()));
    } else {
      return make_any_range_to<edge_descriptor>(std::ranges::empty_view<real_edge_descriptor>());
    }
  }
  [[nodiscard]] std::size_t get_degree(const vertex_descriptor& u) const override {
    if constexpr (concepts::BidirectionalGraph<Graph>) {
      return degree(sub_g_.global_to_local(std::any_cast<real_vertex_descriptor>(u)), sub_g_);
    } else {
      return 0;
    }
  }

  // AdjacencyMatrix
  [[nodiscard]] std::pair<edge_descriptor, bool> get_edge(const vertex_descriptor& u,
                                                          const vertex_descriptor& v) const override {
    if constexpr (concepts::AdjacencyMatrix<Graph>) {
      auto [e_global, e_found] = edge(std::any_cast<real_vertex_descriptor>(u),
                                      std::any_cast<real_vertex_descriptor>(v), sub_g_.root().underlying());
      if (e_found) {
        return {edge_descriptor{e_global}, sub_g_.find_edge(e_global).second};
      }
    }
    return {edge_descriptor{real_edge_descriptor{}}, false};
  }

 protected:
  SubGraph& sub_g_;
  std::vector<std::unique_ptr<DynBase>> children_;
  ParentObs& parent_;

  [[nodiscard]] auto real_to_any_vertex_desc() const {
    return [this](const real_vertex_descriptor& v) { return vertex_descriptor{sub_g_.local_to_global(v)}; };
  }
  [[nodiscard]] auto real_to_any_edge_desc() const {
    return [this](const real_edge_descriptor& e) { return edge_descriptor{sub_g_.local_to_global(e)}; };
  }
};

template <concepts::Graph Graph, typename DynBase>
class dynamic_graph_observer_wrapper<const subgraph<Graph>, DynBase>
    : public dynamic_subgraph_observer_wrapper<const subgraph<Graph>, const Graph, DynBase> {
  using Base = dynamic_subgraph_observer_wrapper<const subgraph<Graph>, const Graph, DynBase>;

 public:
  dynamic_graph_observer_wrapper(const subgraph<Graph>& sub_g, const dynamic_properties& dp) : Base(sub_g, dp) {}
  dynamic_graph_observer_wrapper(const subgraph<Graph>& sub_g, const dynamic_properties& dp, const DynBase& parent)
      : Base(sub_g, dp, parent) {}
};

template <concepts::Graph Graph, typename DynBase>
class dynamic_graph_observer_wrapper<subgraph<Graph>, DynBase>
    : public dynamic_subgraph_observer_wrapper<subgraph<Graph>, Graph, DynBase> {
  using Base = dynamic_subgraph_observer_wrapper<subgraph<Graph>, Graph, DynBase>;

 public:
  dynamic_graph_observer_wrapper(subgraph<Graph>& sub_g, dynamic_properties& dp) : Base(sub_g, dp) {}
  dynamic_graph_observer_wrapper(subgraph<Graph>& sub_g, dynamic_properties& dp, DynBase& parent)
      : Base(sub_g, dp, parent) {}
};

template <concepts::MutableGraph Graph>
class dynamic_graph_mutator_wrapper : public dynamic_graph_observer_wrapper<Graph, dynamic_graph_mutator> {
  using Base = dynamic_graph_observer_wrapper<Graph, dynamic_graph_mutator>;
  using real_vertex_descriptor = graph_vertex_descriptor_t<Graph>;
  using real_edge_descriptor = graph_edge_descriptor_t<Graph>;

 public:
  using vertex_descriptor = typename Base::vertex_descriptor;
  using edge_descriptor = typename Base::edge_descriptor;

  dynamic_graph_mutator_wrapper(Graph& g, dynamic_properties& dp) : Base(g, dp) {}

  [[nodiscard]] const Graph& real_graph() const { return this->g_; }
  [[nodiscard]] Graph& real_graph() { return this->g_; }

  [[nodiscard]] dynamic_properties& get_mutable_properties() override {
    return this->dp_;
  }

  vertex_descriptor do_add_vertex() override { return vertex_descriptor{add_vertex(this->g_)}; }

  void do_clear_vertex(const vertex_descriptor& u) override {
    clear_vertex(std::any_cast<real_vertex_descriptor>(u), this->g_);
  }

  void do_remove_vertex(const vertex_descriptor& u) override {
    remove_vertex(std::any_cast<real_vertex_descriptor>(u), this->g_);
  }

  std::pair<edge_descriptor, bool> do_add_edge(const vertex_descriptor& source,
                                               const vertex_descriptor& target) override {
    auto [e, added] = add_edge(std::any_cast<real_vertex_descriptor>(source),
                               std::any_cast<real_vertex_descriptor>(target), this->g_);
    return {edge_descriptor{e}, added};
  }

  void do_remove_edge(const edge_descriptor& e) override {
    remove_edge(std::any_cast<real_edge_descriptor>(e), this->g_);
  }

  void do_remove_out_edge_if(const vertex_descriptor& u, std::function<bool(const edge_descriptor&)> pred) override {
    if constexpr (concepts::MutableIncidenceGraph<Graph>) {
      remove_out_edge_if(
          std::any_cast<real_vertex_descriptor>(u),
          [p = std::move(pred)](const real_edge_descriptor& e) { return p(edge_descriptor{e}); }, this->g_);
    }
  }

  void do_remove_in_edge_if(const vertex_descriptor& u, std::function<bool(const edge_descriptor&)> pred) override {
    if constexpr (concepts::MutableBidirectionalGraph<Graph>) {
      remove_in_edge_if(
          std::any_cast<real_vertex_descriptor>(u),
          [p = std::move(pred)](const real_edge_descriptor& e) { return p(edge_descriptor{e}); }, this->g_);
    }
  }

  void do_remove_edge_if(std::function<bool(const edge_descriptor&)> pred) override {
    if constexpr (concepts::MutableEdgeListGraph<Graph>) {
      remove_edge_if([p = std::move(pred)](const real_edge_descriptor& e) { return p(edge_descriptor{e}); }, this->g_);
    }
  }

  void set_graph_property(const std::string& name, const std::string& value, const std::string& value_type) override {
    bool type_found = put_property(name, &this->g_, value, value_type);
    if (!type_found) {
      throw std::bad_any_cast();
    }
  }

  void set_vertex_property(const std::string& name, const vertex_descriptor& vertex, const std::string& value,
                           const std::string& value_type) override {
    bool type_found = put_property(name, std::any_cast<real_vertex_descriptor>(vertex), value, value_type);
    if (!type_found) {
      throw std::bad_any_cast();
    }
  }

  void set_edge_property(const std::string& name, const edge_descriptor& edge, const std::string& value,
                         const std::string& value_type) override {
    bool type_found = put_property(name, std::any_cast<real_edge_descriptor>(edge), value, value_type);
    if (!type_found) {
      throw std::bad_any_cast();
    }
  }

  template <typename Key>
  bool put_property(const std::string& name, const Key& key, const std::string& value, const std::string& value_type) {
    if (value_type == "boolean") {
      put(name, this->dp_, key, dynamic_pmap_detail::read_value<bool>(value));
    } else if (value_type == "int") {
      put(name, this->dp_, key, dynamic_pmap_detail::read_value<std::int32_t>(value));
    } else if (value_type == "long") {
      put(name, this->dp_, key, dynamic_pmap_detail::read_value<std::int64_t>(value));
    } else if (value_type == "float") {
      put(name, this->dp_, key, dynamic_pmap_detail::read_value<float>(value));
    } else if (value_type == "double") {
      put(name, this->dp_, key, dynamic_pmap_detail::read_value<double>(value));
    } else if (value_type == "string") {
      put(name, this->dp_, key, dynamic_pmap_detail::read_value<std::string>(value));
    } else {
      return false;
    }
    return true;
  }
};

/***************************** IncidenceGraph *************************************************/

inline auto source(const dynamic_graph_observer::edge_descriptor& e, const dynamic_graph_observer& g) {
  return g.get_source(e);
}

inline auto target(const dynamic_graph_observer::edge_descriptor& e, const dynamic_graph_observer& g) {
  return g.get_target(e);
}

inline auto out_edges(const dynamic_graph_observer::vertex_descriptor& v, const dynamic_graph_observer& g) {
  return g.get_out_edges(v);
}

inline std::size_t out_degree(const dynamic_graph_observer::vertex_descriptor& v, const dynamic_graph_observer& g) {
  return g.get_out_degree(v);
}

/***************************** BidirectionalGraph *********************************************/

inline auto in_edges(const dynamic_graph_observer::vertex_descriptor& v, const dynamic_graph_observer& g) {
  return g.get_in_edges(v);
}

inline std::size_t in_degree(const dynamic_graph_observer::vertex_descriptor& v, const dynamic_graph_observer& g) {
  return g.get_in_degree(v);
}

inline std::size_t degree(const dynamic_graph_observer::vertex_descriptor& v, const dynamic_graph_observer& g) {
  return g.get_degree(v);
}

/***************************** VertexListGraph ************************************************/

inline auto vertices(const dynamic_graph_observer& g) { return g.get_vertices(); }

inline auto num_vertices(const dynamic_graph_observer& g) { return g.get_num_vertices(); }

inline auto vertex(std::size_t i, const dynamic_graph_observer& g) { return *std::next(g.get_vertices().begin(), i); }

/***************************** EdgeListGraph **************************************************/

inline auto edges(const dynamic_graph_observer& g) { return g.get_edges(); }

inline auto num_edges(const dynamic_graph_observer& g) { return g.get_num_edges(); }

/***************************** AdjacencyGraph *************************************************/

inline auto adjacent_vertices(const dynamic_graph_observer::vertex_descriptor& v, const dynamic_graph_observer& g) {
  return adjacency_range(g.get_out_edges(v), g);
}

/***************************** InvAdjacencyGraph **********************************************/

inline auto inv_adjacent_vertices(const dynamic_graph_observer::vertex_descriptor& v, const dynamic_graph_observer& g) {
  return inv_adjacency_range(g.get_in_edges(v), g);
}

/***************************** AdjacencyMatrix ************************************************/

inline auto edge(const dynamic_graph_observer::vertex_descriptor& u, const dynamic_graph_observer::vertex_descriptor& v,
                 const dynamic_graph_observer& g) {
  return g.get_edge(u, v);
}

/***************************** MutableGraph ***************************************************/

inline void clear_vertex(const dynamic_graph_mutator::vertex_descriptor& v, dynamic_graph_mutator& g) {
  g.do_clear_vertex(v);
}

inline void remove_vertex(const dynamic_graph_mutator::vertex_descriptor& v, dynamic_graph_mutator& g) {
  g.do_remove_vertex(v);
}

inline void remove_edge(const dynamic_graph_mutator::vertex_descriptor& u,
                        const dynamic_graph_mutator::vertex_descriptor& v, dynamic_graph_mutator& g) {
  auto [e, e_found] = g.get_edge(u, v);
  while (e_found) {
    g.do_remove_edge(e);
    std::tie(e, e_found) = g.get_edge(u, v);
  }
}

inline void remove_edge(const dynamic_graph_mutator::edge_descriptor& e, dynamic_graph_mutator& g) {
  g.do_remove_edge(e);
}

/***************************** MutableIncidenceGraph ******************************************/

template <typename EdgePred>
void remove_out_edge_if(const dynamic_graph_mutator::vertex_descriptor& u, EdgePred pred, dynamic_graph_mutator& g) {
  g.do_remove_out_edge_if(u, pred);
}

/***************************** MutableBidirectionalGraph **************************************/

template <typename EdgePred>
void remove_in_edge_if(const dynamic_graph_mutator::vertex_descriptor& v, EdgePred pred, dynamic_graph_mutator& g) {
  g.do_remove_in_edge_if(v, pred);
}

/***************************** MutableEdgeListGraph *******************************************/

template <typename EdgePred>
void remove_edge_if(EdgePred pred, dynamic_graph_mutator& g) {
  g.do_remove_edge_if(pred);
}

}  // namespace bagl

#endif  // BAGL_BAGL_DYNAMIC_GRAPH_H_
