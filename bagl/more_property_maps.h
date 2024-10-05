// Copyright 2013-2024 Mikael Persson.

#ifndef BAGL_BAGL_MORE_PROPERTY_MAPS_H_
#define BAGL_BAGL_MORE_PROPERTY_MAPS_H_

#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

namespace bagl {

//======== Whole-bundle property-map ==========

// This property-map delivers the entire vertex or edge bundle associated to a
// vertex / edge descriptor of a graph.
// This property-map relies of the existence of an operator[] for the graph and
// for the relevant descriptor type. This operator[] should deliver the vertex/edge bundle.
template <typename Graph, typename PropertyMapTag>
struct whole_bundle_property_map : public put_get_helper<whole_bundle_property_map<Graph, PropertyMapTag>> {
 public:
  static constexpr bool is_vertex_bundle_v = std::is_same_v<PropertyMapTag, vertex_bundle_t>;
  static constexpr bool is_edge_bundle_v = std::is_same_v<PropertyMapTag, edge_bundle_t>;
  static constexpr bool is_const_graph_v = std::is_const_v<Graph>;
  using non_const_graph = std::remove_cv_t<Graph>;
  using value_type = std::conditional_t<is_vertex_bundle_v, vertex_bundle_type<non_const_graph>, std::conditional_t<is_edge_bundle_v, edge_bundle_type<non_const_graph>, graph_bundle_type<non_const_graph>>>;
  using reference = std::conditional_t<is_const_graph_v, const value_type&, value_type&>;
  using key_type = std::conditional_t<is_vertex_bundle_v, graph_vertex_descriptor_t<non_const_graph>,
                                      std::conditional_t<is_edge_bundle_v, graph_edge_descriptor_t<non_const_graph>, Graph&>>;
  using category = std::conditional_t<is_const_graph_v, readable_property_map_tag, lvalue_property_map_tag>;

  explicit whole_bundle_property_map(Graph* pg) : pg_(pg) {}
  whole_bundle_property_map() = default;
  reference operator[](key_type k) const {
    if constexpr (is_vertex_bundle_v || is_edge_bundle_v) {
      return (*pg_)[k];
    } else {
      return k[graph_bundle];
    }
  }

 private:
  Graph* pg_ = nullptr;
};

//======== Tagged-in-property property-map ==========

// This property-map uses a graph's "get_property" function to get the
// whole property and then extract the value associated to a given tag.
template <typename T, typename Graph, typename PropertyMapTag>
struct tagged_in_property_property_map
    : public put_get_helper<tagged_in_property_property_map<T, Graph, PropertyMapTag>> {
 private:
  Graph* pg_ = nullptr;

 public:
  static constexpr bool is_vertex_prop_v = std::is_same_v<property_kind_t<PropertyMapTag>, vertex_property_tag>;
  static constexpr bool is_edge_prop_v = std::is_same_v<property_kind_t<PropertyMapTag>, edge_property_tag>;
  static constexpr bool is_const_graph_v = std::is_const_v<Graph>;
  using non_const_graph = std::remove_cv_t<Graph>;
  using value_type = T;
  using reference = T&;
  using key_type =
      std::conditional_t<is_vertex_prop_v, graph_vertex_descriptor_t<non_const_graph>, std::conditional_t<is_edge_prop_v, graph_edge_descriptor_t<non_const_graph>, Graph&>>;
  using category = std::conditional_t<std::is_const_v<T>, readable_property_map_tag, lvalue_property_map_tag>;

  explicit tagged_in_property_property_map(Graph* pg, PropertyMapTag /*tag*/ = {}) : pg_(pg) {}
  tagged_in_property_property_map() = default;
  reference operator[](key_type k) const {
    if constexpr (is_vertex_prop_v || is_edge_prop_v) {
      return get_property_value(pg_->get_property(k), PropertyMapTag{});
    } else {
      return get_property_value(k.get_property(graph_all), PropertyMapTag{});
    }
  }
};

//======== Property-graph property-map ==========

// This property-map uses a graph's "get" function to obtain the
// property value associated to a given descriptor.
template <typename T, typename Graph, typename PropertyMapTag>
struct propgraph_property_map : public put_get_helper<propgraph_property_map<T, Graph, PropertyMapTag>> {
 private:
  Graph* pg_ = nullptr;

 public:
  static constexpr bool is_vertex_prop_v = std::is_same_v<property_kind_t<PropertyMapTag>, vertex_property_tag>;
  using value_type = T;
  using reference = T&;
  using key_type =
      std::conditional_t<is_vertex_prop_v, graph_vertex_descriptor_t<Graph>, graph_edge_descriptor_t<Graph>>;
  using category = std::conditional_t<std::is_const_v<T>, readable_property_map_tag, lvalue_property_map_tag>;

  explicit propgraph_property_map(Graph* pg, PropertyMapTag /*tag*/ = {}) : pg_(pg) {}
  propgraph_property_map() = default;
  reference operator[](key_type k) const { return get(PropertyMapTag{}, *pg_, k); }
};

//======== Bundle-data-member property-map ==========

// This property-map delivers a data-member of the vertex or edge bundle associated to a
// vertex / edge descriptor of a graph. This is similar to the property-map obtained by
// calling get(&SomeBundle::SomeMember, my_graph), and can be used to implement such a
// functionality for custom graph classes.
// This property-map relies of the existence of an operator[] for the graph and
// for the relevant descriptor type. This operator[] should deliver the vertex/edge bundle.
template <typename T, typename Graph, typename PropertyMapTag>
class bundle_member_property_map : public put_get_helper<bundle_member_property_map<T, Graph, PropertyMapTag>> {
 public:
  using self = bundle_member_property_map<T, Graph, PropertyMapTag>;
  static constexpr bool is_vertex_bundle_v = std::is_same_v<PropertyMapTag, vertex_bundle_t>;
  static constexpr bool is_edge_bundle_v = std::is_same_v<PropertyMapTag, edge_bundle_t>;
  static constexpr bool is_const_graph_v = std::is_const_v<Graph>;
  using non_const_graph = std::remove_cv_t<Graph>;
  using bundle_type =
      std::conditional_t<is_vertex_bundle_v, vertex_bundle_type<non_const_graph>,
                         std::conditional_t<is_edge_bundle_v, edge_bundle_type<non_const_graph>, graph_bundle_type<non_const_graph>>>;
  using member_ptr_type = T bundle_type::*;

 private:
  Graph* pg_ = nullptr;
  member_ptr_type mem_ptr_ = nullptr;

 public:
  using value_type = T;
  using reference = T&;
  using key_type = std::conditional_t<is_vertex_bundle_v, graph_vertex_descriptor_t<non_const_graph>,
                                      std::conditional_t<is_edge_bundle_v, graph_edge_descriptor_t<non_const_graph>, Graph&>>;
  using category = std::conditional_t<is_const_graph_v, readable_property_map_tag, lvalue_property_map_tag>;

  bundle_member_property_map(Graph* pg, member_ptr_type mem_ptr) : pg_(pg), mem_ptr_(mem_ptr) {}
  bundle_member_property_map() = default;
  reference operator[](key_type k) const {
    if constexpr (is_vertex_bundle_v || is_edge_bundle_v) {
      return (*pg_)[k].*mem_ptr_;
    } else {
      return k[graph_bundle].*mem_ptr_;
    }
  }
};

//======== Sub-object put-get helper ==========

// This is a CRTP-style base-class that can be used to imbue a property-map class
// with put and get functions. The property-map class is only required to provide
// an operator[].
// This helper class applies to property-maps whose value is a subobject of the key object.
template <typename LvaluePropertyMap>
struct subobject_put_get_helper {};

template <typename PropertyMap, typename K>
const auto& get(const subobject_put_get_helper<PropertyMap>& pa, const K& k) {
  return static_cast<const PropertyMap&>(pa)[k];
}

template <typename PropertyMap, typename K>
auto& get(const subobject_put_get_helper<PropertyMap>& pa, K& k) {
  return static_cast<const PropertyMap&>(pa)[k];
}

template <typename PropertyMap, typename K, typename V>
void put(const subobject_put_get_helper<PropertyMap>& pa, K& k, V&& v) {
  static_cast<const PropertyMap&>(pa)[k] = std::forward<V>(v);
}

//======== Self property-map ==========

// This property-map is an "identity" map that preserves the reference
// semantics. In other words, it takes a key by reference and delivers
// it back unchanged (by reference).
template <typename T>
struct self_property_map : public subobject_put_get_helper<self_property_map<T>> {
  using self = self_property_map;

  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using key_type = T;
  using category = std::conditional_t<std::is_const_v<T>, readable_property_map_tag, lvalue_property_map_tag>;

  self_property_map() = default;
  reference operator[](reference p) const { return p; }
  const_reference operator[](const_reference p) const { return p; }
};

//======== Data-member property-map ==========

// This property-map class can be used to map an object (by reference) to one
// of its data-members (by reference).
template <typename T, typename PropertyType>
class data_member_property_map : public subobject_put_get_helper<data_member_property_map<T, PropertyType>> {
 public:
  using member_ptr_type = T PropertyType::*;
  using self = data_member_property_map<T, PropertyType>;

 private:
  member_ptr_type mem_ptr_ = nullptr;

 public:
  using value_type = T;
  using reference = T&;
  using const_reference = const T&;
  using key_type = PropertyType;
  using category = std::conditional_t<std::is_const_v<T>, readable_property_map_tag, lvalue_property_map_tag>;

  explicit data_member_property_map(member_ptr_type aMemPtr) : mem_ptr_(aMemPtr) {}
  data_member_property_map() = default;
  reference operator[](key_type& p) const { return p.*mem_ptr_; }
  const_reference operator[](const key_type& p) const { return p.*mem_ptr_; }
};

// This property-map class can be used to map an object (by reference) to one
// of its data-members (by reference).
// This specialization applies to const objects (and const data-member).
template <typename T, typename PropertyType>
class data_member_property_map<const T, const PropertyType>
    : public subobject_put_get_helper<data_member_property_map<const T, const PropertyType>> {
 public:
  using value_type = T;
  using reference = const T&;
  using const_reference = const T&;
  using key_type = const PropertyType;

  using member_ptr_type = T key_type::*;
  using self = data_member_property_map<const T, const PropertyType>;

 private:
  member_ptr_type mem_ptr_ = nullptr;

 public:
  using category = std::conditional_t<std::is_const_v<T>, readable_property_map_tag, lvalue_property_map_tag>;

  explicit data_member_property_map(member_ptr_type aMemPtr) : mem_ptr_(aMemPtr) {}
  data_member_property_map() = default;
  reference operator[](key_type& p) const { return p.*mem_ptr_; }
};

//======== Composite property-map ==========

// This property-map allows for the composition of two property-maps.
// The resulting property-map has the value-type of the output-map and
// the key-type of the input-map.
template <typename OutputMap, typename InputMap>
class composite_property_map {
 public:  // private:  would be private is friends were more portable.
  OutputMap prop_out;
  InputMap prop_in;

  using value_type = property_traits_value_t<OutputMap>;
  using key_type = std::remove_cv_t<property_traits_key_t<InputMap>>;
  using category = property_traits_category_t<OutputMap>;
  using reference = property_traits_reference_t<OutputMap>;
  using const_reference = const reference;

  explicit composite_property_map(OutputMap aPropOut, InputMap aPropIn = {}) : prop_out(aPropOut), prop_in(aPropIn) {}
  composite_property_map() = default;

  reference operator[](const key_type& k) const { return prop_out[prop_in[k]]; }
  reference operator[](key_type& k) const { return prop_out[prop_in[k]]; }
};

template <typename OutputMap, typename InputMap>
auto get(const composite_property_map<OutputMap, InputMap>& m,
         const typename composite_property_map<OutputMap, InputMap>::key_type& p) {
  return m.prop_out[m.prop_in[p]];
}

template <typename OutputMap, typename InputMap, typename V>
void put(const composite_property_map<OutputMap, InputMap>& m,
         const typename composite_property_map<OutputMap, InputMap>::key_type& p, V&& value) {
  put(m.prop_out, m.prop_in[p], std::forward<V>(value));
}

template <typename OutputMap, typename InputMap, typename V>
void put(const composite_property_map<OutputMap, InputMap>& m,
         typename composite_property_map<OutputMap, InputMap>::key_type& p, V&& value) {
  put(m.prop_out, m.prop_in[p], std::forward<V>(value));
}

// This function template can be used to construct a property-map that takes a
// vertex-descriptor into a graph and maps it to a data member of the vertex-bundle.
// This property-map is constructed from bundle-to-data-member map and a graph.
// `bundle_prop`: The bundle-to-data-member property-map.
// `g`: The graph.
// Returns a property-map that can map a vertex-descriptor to a vertex-bundle data-member.
template <typename BundleMemberMap, typename Graph>
auto bundle_prop_to_vertex_prop(BundleMemberMap bundle_prop, Graph& g) {
  return composite_property_map(bundle_prop, whole_bundle_property_map<Graph, vertex_bundle_t>(&g));
}

// This function template can be used to construct a property-map that takes a
// edge-descriptor into a graph and maps it to a data member of the edge-bundle.
// This property-map is constructed from bundle-to-data-member map and a graph.
// `bundle_prop`: The bundle-to-data-member property-map.
// `g`: The graph.
// Returns a property-map that can map a edge-descriptor to a edge-bundle data-member.
template <typename BundleMemberMap, typename Graph>
auto bundle_prop_to_edge_prop(BundleMemberMap bundle_prop, Graph& g) {
  return composite_property_map(bundle_prop, whole_bundle_property_map<Graph, edge_bundle_t>(&g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_MORE_PROPERTY_MAPS_H_
