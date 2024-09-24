// Copyright 2007 Douglas Gregor
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_NAMED_GRAPH_H_
#define BAGL_BAGL_NAMED_GRAPH_H_

#include <cassert>
#include <functional>
#include <optional>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include "bagl/container_traits.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property.h"

namespace bagl {

// Provides support for named vertices in graphs, allowing one to more
// easily associate unique external names (URLs, city names, employee
// ID numbers, etc.) with the vertices of a graph.

/*******************************************************************
 * User-customized traits                                          *
 *******************************************************************/

/**
 * @brief Trait used to extract the internal vertex name from a vertex
 * property.
 *
 * To enable the use of internal vertex names in a graph type,
 * specialize the @c internal_vertex_name trait for your graph
 * property (e.g., @c a City class, which stores information about the
 * vertices in a road map).
 */
template <typename VertexProperty>
struct internal_vertex_name {
  /**
   *  The @c type field provides a function object that extracts a key
   *  from the @c VertexProperty. The function object type must have a
   *  nested @c result_type that provides the type of the key. For
   *  more information, see the @c KeyExtractor concept in the
   *  Boost.MultiIndex documentation: @c type must either be @c void
   *  (if @c VertexProperty does not have an internal vertex name) or
   *  a model of @c KeyExtractor.
   */
  using type = void;
};

/**
 * Extract the internal vertex name from a @c property structure by
 * looking at its base.
 */
template <typename Tag, typename T, typename Base>
struct internal_vertex_name<property<Tag, T, Base>> : internal_vertex_name<Base> {};

template <typename VertexProperty>
using internal_vertex_name_t = typename internal_vertex_name<VertexProperty>::type;

/**
 * Construct an instance of @c VertexProperty directly from its
 * name. This function object should be used within the @c
 * internal_vertex_constructor trait.
 */
template <typename VertexProperty>
struct vertex_from_name {
 private:
  using extract_name_type = internal_vertex_name_t<VertexProperty>;

  using vertex_name_type = std::remove_cv_t<std::remove_reference_t<typename extract_name_type::result_type>>;

 public:
  using argument_type = vertex_name_type;
  using result_type = VertexProperty;

  VertexProperty operator()(const vertex_name_type& name) { return VertexProperty(name); }
};

/**
 * Throw an exception whenever one tries to construct a @c
 * VertexProperty instance from its name.
 */
template <typename VertexProperty>
struct cannot_add_vertex {
 private:
  using extract_name_type = internal_vertex_name_t<VertexProperty>;

  using vertex_name_type = std::remove_cv_t<std::remove_reference_t<typename extract_name_type::result_type>>;

 public:
  using argument_type = vertex_name_type;
  using result_type = VertexProperty;

  VertexProperty operator()(const vertex_name_type&) {
    assert(false && "add_vertex: unable to create a vertex from its name");
  }
};

/**
 * @brief Trait used to construct an instance of a @c VertexProperty,
 * which is a class type that stores the properties associated with a
 * vertex in a graph, from just the name of that vertex property. This
 * operation is used when an operation is required to map from a
 * vertex name to a vertex descriptor (e.g., to add an edge outgoing
 * from that vertex), but no vertex by the name exists. The function
 * object provided by this trait will be used to add new vertices
 * based only on their names. Since this cannot be done for arbitrary
 * types, the default behavior is to throw an exception when this
 * routine is called, which requires that all named vertices be added
 * before adding any edges by name.
 */
template <typename VertexProperty>
struct internal_vertex_constructor {
  /**
   * The @c type field provides a function object that constructs a
   * new instance of @c VertexProperty from the name of the vertex (as
   * determined by @c internal_vertex_name). The function object shall
   * accept a vertex name and return a @c VertexProperty. Predefined
   * options include:
   *
   *   @c vertex_from_name<VertexProperty>: construct an instance of
   *   @c VertexProperty directly from the name.
   *
   *   @c cannot_add_vertex<VertexProperty>: the default value, which
   *   throws an @c std::runtime_error if one attempts to add a vertex
   *   given just the name.
   */
  using type = cannot_add_vertex<VertexProperty>;
};

/**
 * Extract the internal vertex constructor from a @c property structure by
 * looking at its base.
 */
template <typename Tag, typename T, typename Base>
struct internal_vertex_constructor<property<Tag, T, Base>> : internal_vertex_constructor<Base> {};

/*******************************************************************
 * Named graph mixin                                               *
 *******************************************************************/

/**
 * named_graph is a mixin that provides names for the vertices of a
 * graph, including a mapping from names to vertices. Graph types that
 * may or may not be have vertex names (depending on the properties
 * supplied by the user) should use maybe_named_graph.
 *
 * Template parameters:
 *
 *   Graph: the graph type that derives from named_graph
 *
 *   Vertex: the type of a vertex descriptor in Graph. Note: we cannot
 *   use graph_traits here, because the Graph is not yet defined.
 *
 *   VertexProperty: the type stored with each vertex in the Graph.
 */
template <typename Graph, typename Vertex, typename VertexProperty>
class named_graph {
 public:
  /// The type of the function object that extracts names from vertex
  /// properties.
  using extract_name_type = internal_vertex_name_t<VertexProperty>;
  /// The type of the "bundled" property, from which the name can be
  /// extracted.
  using bundled_vertex_property_type = typename lookup_one_property<VertexProperty, vertex_bundle_t>::type;

  /// The type of the function object that generates vertex properties
  /// from names, for the implicit addition of vertices.
  using vertex_constructor_type = typename internal_vertex_constructor<VertexProperty>::type;

  /// The type used to name vertices in the graph
  using vertex_name_type = std::decay_t<typename extract_name_type::result_type>;

  /// The type of vertex descriptors in the graph
  using vertex_descriptor = Vertex;

  /// The type that maps names to vertices
  using named_vertices_type = std::unordered_map<vertex_name_type, Vertex>;

  /// Construct an instance of the named graph mixin, using the given
  /// function object to extract a name from the bundled property
  /// associated with a vertex.
  explicit named_graph(const extract_name_type& a_extract = extract_name_type(),
                       const vertex_constructor_type& a_vertex_constructor = vertex_constructor_type())
      : extract(a_extract), named_vertices(), vertex_constructor(a_vertex_constructor) {}

  /// Notify the named_graph that we have added the given vertex. The
  /// name of the vertex will be added to the mapping.
  void added_vertex(Vertex vertex) { named_vertices.emplace(extract_name(derived()[vertex]), vertex); }

  /// Notify the named_graph that we are removing the given
  /// vertex. The name of the vertex will be removed from the mapping.
  template <typename VertexIterStability>
  void removing_vertex(Vertex vertex, VertexIterStability) {
    static_assert(std::is_base_of_v<container_traits_detail::stable_tag, VertexIterStability>,
                  "Named graphs cannot use vecS as vertex container and remove "
                  "vertices; the lack of vertex descriptor stability (which iterator "
                  "stability is a proxy for) means that the name -> vertex mapping "
                  "would need to be completely rebuilt after each deletion.  See "
                  "https://svn.boost.org/trac/boost/ticket/7863 for more information "
                  "and a test case.");
    const auto& vname = extract_name(derived()[vertex]);
    named_vertices.erase(vname);
  }

  /// Notify the named_graph that we are clearing the graph.
  /// This will clear out all of the name->vertex mappings
  void clearing_graph() { named_vertices.clear(); }

  /// Retrieve the derived instance
  Graph& derived() { return static_cast<Graph&>(*this); }
  const Graph& derived() const { return static_cast<const Graph&>(*this); }

  /// Extract the name from a vertex property instance
  auto extract_name(const bundled_vertex_property_type& property) { return extract(property); }

  /// Search for a vertex that has the given property (based on its
  /// name)
  std::optional<vertex_descriptor> vertex_by_property(const bundled_vertex_property_type& property) {
    return find_vertex(extract_name(property), *this);
  }

  /// The function object that extracts names from vertex properties.
  extract_name_type extract;

  /// Mapping from names to vertices
  named_vertices_type named_vertices;

  /// Constructs a vertex from the name of that vertex
  vertex_constructor_type vertex_constructor;
};

/// Helper macro containing the template parameters of named_graph
#define BGL_NAMED_GRAPH_PARAMS typename Graph, typename Vertex, typename VertexProperty
/// Helper macro containing the named_graph<...> instantiation
#define BGL_NAMED_GRAPH named_graph<Graph, Vertex, VertexProperty>

/// Retrieve the vertex associated with the given name
template <BGL_NAMED_GRAPH_PARAMS>
std::optional<Vertex> find_vertex(typename BGL_NAMED_GRAPH::vertex_name_type const& name, const BGL_NAMED_GRAPH& g) {
  /// Look for a vertex with the given name
  auto iter = g.named_vertices.find(name);
  if (iter == g.named_vertices.end()) {
    return std::nullopt;  // vertex not found
  }
  return *iter;
}

/// Retrieve the vertex associated with the given name, or add a new
/// vertex with that name if no such vertex is available.
/// Note: This is enabled only when the vertex property type is different
///       from the vertex name to avoid ambiguous overload problems with
///       the add_vertex() function that takes a vertex property.
template <BGL_NAMED_GRAPH_PARAMS>
std::enable_if_t<!std::is_same_v<typename BGL_NAMED_GRAPH::vertex_name_type, VertexProperty>, Vertex> add_vertex(
    typename BGL_NAMED_GRAPH::vertex_name_type const& name, BGL_NAMED_GRAPH& g) {
  if (std::optional<Vertex> vertex = find_vertex(name, g); vertex.has_value()) {
    /// We found the vertex, so return it
    return vertex.value();
  }
  /// There is no vertex with the given name, so create one
  return add_vertex(g.vertex_constructor(name), g.derived());
 
}

/// Add an edge using vertex names to refer to the vertices
template <BGL_NAMED_GRAPH_PARAMS>
auto add_edge(typename BGL_NAMED_GRAPH::vertex_name_type const& u_name,
              typename BGL_NAMED_GRAPH::vertex_name_type const& v_name, BGL_NAMED_GRAPH& g) {
  return add_edge(add_vertex(u_name, g.derived()), add_vertex(v_name, g.derived()), g.derived());
}

/// Add an edge using vertex descriptors or names to refer to the vertices
template <BGL_NAMED_GRAPH_PARAMS>
auto add_edge(typename BGL_NAMED_GRAPH::vertex_descriptor const& u,
              typename BGL_NAMED_GRAPH::vertex_name_type const& v_name, BGL_NAMED_GRAPH& g) {
  return add_edge(u, add_vertex(v_name, g.derived()), g.derived());
}

/// Add an edge using vertex descriptors or names to refer to the vertices
template <BGL_NAMED_GRAPH_PARAMS>
auto add_edge(typename BGL_NAMED_GRAPH::vertex_name_type const& u_name,
              typename BGL_NAMED_GRAPH::vertex_descriptor const& v, BGL_NAMED_GRAPH& g) {
  return add_edge(add_vertex(u_name, g.derived()), v, g.derived());
}

// Overloads to support EdgeMutablePropertyGraph graphs
template <BGL_NAMED_GRAPH_PARAMS>
auto add_edge(typename BGL_NAMED_GRAPH::vertex_descriptor const& u,
              typename BGL_NAMED_GRAPH::vertex_name_type const& v_name,
              typename edge_property_type<Graph>::type const& p, BGL_NAMED_GRAPH& g) {
  return add_edge(u, add_vertex(v_name, g.derived()), p, g.derived());
}

template <BGL_NAMED_GRAPH_PARAMS>
auto add_edge(typename BGL_NAMED_GRAPH::vertex_name_type const& u_name,
              typename BGL_NAMED_GRAPH::vertex_descriptor const& v, typename edge_property_type<Graph>::type const& p,
              BGL_NAMED_GRAPH& g) {
  return add_edge(add_vertex(u_name, g.derived()), v, p, g.derived());
}

template <BGL_NAMED_GRAPH_PARAMS>
auto add_edge(typename BGL_NAMED_GRAPH::vertex_name_type const& u_name,
              typename BGL_NAMED_GRAPH::vertex_name_type const& v_name,
              typename edge_property_type<Graph>::type const& p, BGL_NAMED_GRAPH& g) {
  return add_edge(add_vertex(u_name, g.derived()), add_vertex(v_name, g.derived()), p, g.derived());
}

#undef BGL_NAMED_GRAPH
#undef BGL_NAMED_GRAPH_PARAMS

/*******************************************************************
 * Maybe named graph mixin                                         *
 *******************************************************************/

/**
 * A graph mixin that can provide a mapping from names to vertices,
 * and use that mapping to simplify creation and manipulation of
 * graphs.
 */
template <typename Graph, typename Vertex, typename VertexProperty,
          typename ExtractName = internal_vertex_name_t<VertexProperty>>
struct maybe_named_graph : public named_graph<Graph, Vertex, VertexProperty> {};

/**
 * A graph mixin that can provide a mapping from names to vertices,
 * and use that mapping to simplify creation and manipulation of
 * graphs. This partial specialization turns off this functionality
 * when the @c VertexProperty does not have an internal vertex name.
 */
template <typename Graph, typename Vertex, typename VertexProperty>
struct maybe_named_graph<Graph, Vertex, VertexProperty, void> {
  /// The type of the "bundled" property, from which the name can be
  /// extracted.
  using bundled_vertex_property_type = typename lookup_one_property<VertexProperty, vertex_bundle_t>::type;

  /// Notify the named_graph that we have added the given vertex. This
  /// is a no-op.
  void added_vertex(Vertex) {}

  /// Notify the named_graph that we are removing the given
  /// vertex. This is a no-op.
  template <typename VertexIterStability>
  void removing_vertex(Vertex, VertexIterStability) {}

  /// Notify the named_graph that we are clearing the graph. This is a
  /// no-op.
  void clearing_graph() {}

  /// Search for a vertex that has the given property (based on its
  /// name). This always returns nullopt
  std::optional<Vertex> vertex_by_property(const bundled_vertex_property_type&) { return std::nullopt; }
};

}  // namespace bagl

#endif  // BAGL_BAGL_NAMED_GRAPH_H_
