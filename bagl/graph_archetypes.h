// Copyright 2002 Indiana University.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GRAPH_ARCHETYPES_H_
#define BAGL_BAGL_GRAPH_ARCHETYPES_H_

#include "bagl/property_map.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"

namespace bagl {

namespace graph_archetypes_detail {
class dummy_constructor {};
template <class T = int>
class null_archetype {
 private:
  null_archetype() = default;
  null_archetype(const null_archetype&) = default;
  null_archetype& operator=(const null_archetype&) = default;
  null_archetype(null_archetype&&)  noexcept = default;
  null_archetype& operator=(null_archetype&&)  noexcept = default;

 public:
  explicit null_archetype(dummy_constructor /*unused*/) {}
};

struct null_graph_archetype : public null_archetype<> {
  struct traversal_category {};
};

template <class T, int I = 0>
class input_iterator_archetype {
 private:
  using self = input_iterator_archetype;
 public:
  using iterator_category = std::input_iterator_tag;
  using value_type = T;
  struct reference {
    // NOLINTNEXTLINE
    operator const value_type&() const { return std::declval<T>(); }
  };
  using pointer = const T*;
  using difference_type = std::ptrdiff_t;
  bool operator==(const self&) const { return true; }
  bool operator!=(const self&) const { return true; }
  reference operator*() const { return reference(); }
  self& operator++() { return *this; }
  self operator++(int) { return *this; }
};
}  // namespace graph_archetypes_detail

//===========================================================================
template < typename Vertex, typename Directed, typename ParallelCategory,
    typename Base = graph_archetypes_detail::null_graph_archetype >
struct incidence_graph_archetype : public Base
{
    using base_trav_cat = typename Base::traversal_category;
    struct traversal_category : public incidence_graph_tag, public base_trav_cat
    {
    };
    using vertex_descriptor = Vertex;
    using degree_size_type = unsigned int;
    using vertices_size_type = unsigned int;
    using edges_size_type = unsigned int;
    struct edge_descriptor
    {
        edge_descriptor() = default;
        explicit edge_descriptor(const graph_archetypes_detail::dummy_constructor&) {}
        bool operator==(const edge_descriptor&) const { return false; }
        bool operator!=(const edge_descriptor&) const { return false; }
    };
    using out_edge_iterator = graph_archetypes_detail::input_iterator_archetype< edge_descriptor >;

    using directed_category = Directed;
    using edge_parallel_category = ParallelCategory;

    using adjacency_iterator = void;
    using in_edge_iterator = void;
    using vertex_iterator = void;
    using edge_iterator = void;

    static vertex_descriptor null_vertex() { return vertex_descriptor(); }
};
template < typename V, typename D, typename P, typename B >
V source(
    const typename incidence_graph_archetype< V, D, P, B >::edge_descriptor&,
    const incidence_graph_archetype< V, D, P, B >&)
{
    return V(graph_archetypes_detail::dummy_constructor{});
}
template < typename V, typename D, typename P, typename B >
V target(
    const typename incidence_graph_archetype< V, D, P, B >::edge_descriptor&,
    const incidence_graph_archetype< V, D, P, B >&)
{
    return V(graph_archetypes_detail::dummy_constructor{});
}

template < typename V, typename D, typename P, typename B >
auto out_edges(const V&, const incidence_graph_archetype< V, D, P, B >&)
{
    using Iter = typename incidence_graph_archetype< V, D, P, B >::out_edge_iterator;
    return std::pair{Iter{}, Iter{}};
}

template < typename V, typename D, typename P, typename B >
int out_degree(const V&, const incidence_graph_archetype< V, D, P, B >&)
{
    return 0;
}

//===========================================================================
template < typename Vertex, typename Directed, typename ParallelCategory,
    typename Base = graph_archetypes_detail::null_graph_archetype >
struct adjacency_graph_archetype : public Base
{
    using base_trav_cat = typename Base::traversal_category;
    struct traversal_category : public adjacency_graph_tag, public base_trav_cat
    {
    };
    using vertex_descriptor = Vertex;
    using degree_size_type = unsigned int;
    using vertices_size_type = unsigned int;
    using edges_size_type = unsigned int;
    using edge_descriptor = void;
    using adjacency_iterator = graph_archetypes_detail::input_iterator_archetype< Vertex >;

    using directed_category = Directed;
    using edge_parallel_category = ParallelCategory;

    using in_edge_iterator = void;
    using out_edge_iterator = void;
    using vertex_iterator = void;
    using edge_iterator = void;

    static vertex_descriptor null_vertex() { return vertex_descriptor(); }
};

template < typename V, typename D, typename P, typename B >
auto adjacent_vertices(const V&, const adjacency_graph_archetype< V, D, P, B >&)
{
    using Iter = typename adjacency_graph_archetype< V, D, P, B >::adjacency_iterator;
    return std::pair{Iter{}, Iter{}};
}

template < typename V, typename D, typename P, typename B >
int out_degree(const V&, const adjacency_graph_archetype< V, D, P, B >&)
{
    return 0;
}

//===========================================================================
template < typename Vertex, typename Directed, typename ParallelCategory,
    typename Base = graph_archetypes_detail::null_graph_archetype >
struct vertex_list_graph_archetype : public Base
{
    using Incidence = incidence_graph_archetype< Vertex, Directed, ParallelCategory >;
    using Adjacency = adjacency_graph_archetype< Vertex, Directed, ParallelCategory >;

    using base_trav_cat = typename Base::traversal_category;
    struct traversal_category : public vertex_list_graph_tag,
                                public base_trav_cat
    {
    };
    using vertex_descriptor = Vertex;
    using degree_size_type = unsigned int;
    using edge_descriptor = typename Incidence::edge_descriptor;
    using out_edge_iterator = typename Incidence::out_edge_iterator;
    using adjacency_iterator = typename Adjacency::adjacency_iterator;

    using vertex_iterator = graph_archetypes_detail::input_iterator_archetype< Vertex >;
    using vertices_size_type = unsigned int;
    using edges_size_type = unsigned int;

    using directed_category = Directed;
    using edge_parallel_category = ParallelCategory;

    using in_edge_iterator = void;
    using edge_iterator = void;

    static vertex_descriptor null_vertex() { return vertex_descriptor(); }
};

template < typename V, typename D, typename P, typename B >
auto vertices(const vertex_list_graph_archetype< V, D, P, B >&)
{
    using Iter = typename vertex_list_graph_archetype< V, D, P, B >::vertex_iterator;
    return std::pair{Iter{}, Iter{}};
}

template < typename V, typename D, typename P, typename B >
int num_vertices(const vertex_list_graph_archetype< V, D, P, B >&)
{
    return 0;
}

// ambiguously inherited from incidence graph and adjacency graph
template < typename V, typename D, typename P, typename B >
int out_degree(const V&, const vertex_list_graph_archetype< V, D, P, B >&)
{
    return 0;
}

//===========================================================================

struct property_graph_archetype_tag
{
};

template < typename GraphArchetype, typename Property, typename ValueArch >
struct property_graph_archetype : public GraphArchetype
{
    using graph_tag = property_graph_archetype_tag;
    using vertex_property_type = ValueArch;
    using edge_property_type = ValueArch;
};

struct choose_edge_property_map_archetype
{
    template < typename Graph, typename Property, typename Tag > struct bind_
    {
        using type = mutable_lvalue_property_map_archetype<
            typename Graph::edge_descriptor, Property >;
        using const_type = lvalue_property_map_archetype< typename Graph::edge_descriptor,
            Property >;
    };
};
template <> struct edge_property_selector< property_graph_archetype_tag >
{
    using type = choose_edge_property_map_archetype;
};

struct choose_vertex_property_map_archetype
{
    template < typename Graph, typename Property, typename Tag > struct bind_
    {
        using type = mutable_lvalue_property_map_archetype<
            typename Graph::vertex_descriptor, Property >;
        using const_type = lvalue_property_map_archetype<
            typename Graph::vertex_descriptor, Property >;
    };
};

template <> struct vertex_property_selector< property_graph_archetype_tag >
{
    using type = choose_vertex_property_map_archetype;
};

template < typename G, typename P, typename V >
property_map_t< property_graph_archetype< G, P, V >, P > get(
    P, property_graph_archetype< G, P, V >&)
{
    return {};
}

template < typename G, typename P, typename V >
property_map_const_t< property_graph_archetype< G, P, V >, P > get(
    P, const property_graph_archetype< G, P, V >&)
{
    return {};
}

template < typename G, typename P, typename K, typename V >
auto get(P p, const property_graph_archetype< G, P, V >& g, K k)
{
    return get(get(p, g), k);
}

template < typename G, typename P, typename V, typename Key >
void put(
    P p, property_graph_archetype< G, P, V >& g, const Key& key, const V& value)
{
    auto pmap = get(p, g);
    put(pmap, key, value);
}

struct color_value_archetype
{
    color_value_archetype() = default;
    explicit color_value_archetype(graph_archetypes_detail::dummy_constructor) {}
    bool operator==(const color_value_archetype&) const { return true; }
    bool operator!=(const color_value_archetype&) const { return true; }
};
template <> struct color_traits< color_value_archetype >
{
    static color_value_archetype white()
    {
        return color_value_archetype(
            graph_archetypes_detail::dummy_constructor{});
    }
    static color_value_archetype gray()
    {
        return color_value_archetype(
            graph_archetypes_detail::dummy_constructor{});
    }
    static color_value_archetype black()
    {
        return color_value_archetype(
            graph_archetypes_detail::dummy_constructor{});
    }
};

template < typename T > class buffer_archetype
{
public:
    void push(const T&) {}
    void pop() {}
    T& top() { return std::declval< T >(); }
    const T& top() const { return std::declval< T >(); }
    bool empty() const { return true; }
};

} // namespace bagl

#endif // BAGL_BAGL_GRAPH_ARCHETYPES_H_
