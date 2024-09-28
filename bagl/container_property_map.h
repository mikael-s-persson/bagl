// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_CONTAINER_PROPERTY_MAP_H_
#define BAGL_BAGL_CONTAINER_PROPERTY_MAP_H_

#include "bagl/detail/index.h"
#include "bagl/property_map.h"

namespace bagl
{
// This is an adapter built over the iterator property map with
// more useful uniform construction semantics. Specifically, this
// requires the container rather than the iterator and the graph
// rather than the optional index map.
template < typename Graph, typename Key, typename Container >
struct container_property_map
: put_get_helper<container_property_map< Graph, Key, Container > >
{
    using indexer_type = typename index_detail::choose_indexer<Graph, Key>::indexer_type;
    using index_type = typename indexer_type::index_type;
    using map_type = iterator_property_map<typename Container::iterator, property_map_const_t<Graph, index_type>>;
    using key_type = typename map_type::key_type;
    using value_type = typename map_type::value_type;
    using reference = typename map_type::reference;
    using category = typename map_type::category;

    // The default constructor will *probably* crash if its actually
    // used for referencing vertices since the underlying iterator
    // map points past the end of an unknown container.
    container_property_map() = default;

    // This is the preferred constructor. It is invoked over the container
    // and the graph explicitly. This requires that the underlying iterator
    // map use the indices of the vertices in g rather than the default
    // identity map.
    //
    // Note the const-cast this ensures the reference type of the
    // vertex index map is non-const, which happens to be an
    // artifact of passing const graph references.
    container_property_map(Container& c, const Graph& g)
    : m_map(c.begin(), indexer_type::index_map(g))
    {
    }

    // The [] operator delegates to the underlying map/
    reference operator[](const key_type& k) const { return m_map[k]; }

    map_type m_map;
};

} // namespace bagl

#endif // BAGL_BAGL_CONTAINER_PROPERTY_MAP_H_
