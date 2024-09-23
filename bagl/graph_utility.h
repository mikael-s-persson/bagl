// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GRAPH_UTILITY_H_
#define BAGL_BAGL_GRAPH_UTILITY_H_

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <type_traits>

#include "bagl/depth_first_search.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"

namespace bagl {

// Provide an undirected graph interface alternative to the
// the source() and target() edge functions.
template <typename UndirectedGraph>
auto incident(graph_edge_descriptor_t<UndirectedGraph> e, UndirectedGraph& g) {
  return std::make_pair(source(e, g), target(e, g));
}

// Provide an undirected graph interface alternative
// to the out_edges() function.
template <typename Graph>
auto incident_edges(graph_vertex_descriptor_t<Graph> u, Graph& g) {
  return out_edges(u, g);
}

template <typename Graph>
auto opposite(graph_edge_descriptor_t<Graph> e, graph_vertex_descriptor_t<Graph> v, const Graph& g) {
  using vertex_descriptor = graph_vertex_descriptor_t<Graph>;
  if (v == source(e, g)) {
    return target(e, g);
  } else if (v == target(e, g)) {
    return source(e, g);
  } else {
    return vertex_descriptor{};
  }
}

//===========================================================================
// Some handy predicates

template <typename Vertex, typename Graph>
auto incident_from(Vertex u, const Graph& g) {
  return [u, pg = &g](const auto& e) { return source(e, *pg) == u; };
}

template <typename Vertex, typename Graph>
auto incident_to(Vertex u, const Graph& g) {
  return [u, pg = &g](const auto& e) { return target(e, *pg) == u; };
}

template <typename Vertex, typename Graph>
auto incident_on(Vertex u, const Graph& g) {
  return [u, pg = &g](const auto& e) { return source(e, *pg) == u || target(e, *pg) == u; };
}

template <typename Vertex, typename Graph>
auto connects(Vertex u, Vertex v, const Graph& g) {
  return [u, v, pg = &g](const auto& e) {
    if constexpr (is_directed_graph_v<Graph>) {
      return source(e, *pg) == u && target(e, *pg) == v;
    } else {
      return (source(e, *pg) == u && target(e, *pg) == v) || (source(e, *pg) == v && target(e, *pg) == u);
    }
  };
}

template <typename IncidenceGraph, typename Name>
void print_in_edges(const IncidenceGraph& G, Name name, std::ostream& os = std::cout) {
  for (auto u : vertices(G)) {
    os << get(name, u) << " <-- ";
    for (auto e : in_edges(u, G)) {
      os << get(name, source(e, G)) << " ";
    }
    os << '\n';
  }
}

template <typename IncidenceGraph, typename Name>
void print_graph(const IncidenceGraph& G, Name name, std::ostream& os = std::cout) {
  const char* arrow = (is_directed_graph_v<IncidenceGraph> ? " --> " : " <--> ");
  for (auto u : vertices(G)) {
    os << get(name, u) << arrow;
    for (auto e : out_edges(u, G)) {
      os << get(name, target(e, G)) << " ";
    }
    os << '\n';
  }
}
template <typename IncidenceGraph>
void print_graph(const IncidenceGraph& G, std::ostream& os = std::cout) {
  print_graph(G, get(vertex_index, G), os);
}

template <typename EdgeListGraph, typename Name>
void print_edges(const EdgeListGraph& G, Name name, std::ostream& os = std::cout) {
  for (auto e : edges(G)) {
    os << "(" << get(name, source(e, G)) << "," << get(name, target(e, G)) << ") ";
  }
  os << '\n';
}

template <typename EdgeListGraph, typename VertexName, typename EdgeName>
void print_edges2(const EdgeListGraph& G, VertexName vname, EdgeName ename, std::ostream& os = std::cout) {
  for (auto e : edges(G)) {
    os << get(ename, e) << "(" << get(vname, source(e, G)) << "," << get(vname, target(e, G)) << ") ";
  }
  os << '\n';
}

template <typename VertexListGraph, typename Name>
void print_vertices(const VertexListGraph& G, Name name, std::ostream& os = std::cout) {
  for (auto u : vertices(G)) {
    os << get(name, u) << " ";
  }
  os << '\n';
}

template <typename Graph, typename Vertex>
bool is_adjacent(Graph& g, Vertex a, Vertex b) {
  const auto adj_rg = adjacent_vertices(a, g);
  const auto adj_it = std::find(adj_rg.begin(), adj_rg.end(), b);
  if (adj_it == adj_rg.end()) {
    return false;
  }

  const auto oe_rg = out_edges(a, g);
  const auto oe_it = std::find_if(oe_rg.begin(), oe_rg.end(), incident_to(b, g));
  if (oe_it == oe_rg.end()) {
    return false;
  }

  if constexpr (std::is_convertible_v<graph_directed_category_t<Graph>, bidirectional_tag>) {
    const auto ie_rg = in_edges(b, g);
    const auto ie_it = std::find_if(ie_rg.begin(), ie_rg.end(), incident_from(a, g));
    if (ie_it == ie_rg.end()) {
      return false;
    }
  }

  return true;
}

template <typename Graph, typename Edge>
bool in_edge_set(Graph& g, Edge e) {
  const auto e_rg = edges(g);
  return std::find(e_rg.begin(), e_rg.end(), e) != e_rg.end();
}

template <typename Graph, typename Vertex>
bool in_vertex_set(Graph& g, Vertex v) {
  const auto v_rg = vertices(g);
  return std::find(v_rg.begin(), v_rg.end(), v) != v_rg.end();
}

template <typename Graph, typename Vertex>
bool in_edge_set(Graph& g, Vertex u, Vertex v) {
  for (auto e : edges(g)) {
    if (source(e, g) == u && target(e, g) == v) {
      return true;
    }
  }
  return false;
}

// is x a descendant of y?
template <typename ParentMap>
bool is_descendant(property_value_t<ParentMap> x, property_value_t<ParentMap> y, ParentMap parent) {
  // x is the root of the tree
  while (get(parent, x) != x) {
    if (get(parent, x) == y) {
      return true;
    }
    x = get(parent, x);
  }
  return false;
}

// is y reachable from x?
template <typename IncidenceGraph, typename VertexColorMap>
inline bool is_reachable(graph_vertex_descriptor_<IncidenceGraph> x, graph_vertex_descriptor_<IncidenceGraph> y,
                         const IncidenceGraph& g,
                         VertexColorMap color)  // should start out white for every vertex
{
  using ColorValue = property_value_t<VertexColorMap>;
  dfs_visitor<> vis;
  depth_first_visit(g, x, vis, color);
  return get(color, y) != color_traits<ColorValue>::white();
}

// Is the undirected graph connected?
// Is the directed graph strongly connected?
template <typename VertexListGraph, typename VertexColorMap>
inline bool is_connected(const VertexListGraph& g, VertexColorMap color) {
  using ColorValue = property_value_t<VertexColorMap>;
  using Color = color_traits<ColorValue>;
  for (auto u : vertices(g)) {
    for (auto v : vertices(g)) {
      if (u != v) {
        for (auto c : vertices(g)) {
          put(color, c, Color::white());
        }
        if (!is_reachable(u, v, g, color)) {
          return false;
        }
      }
    }
  }
  return true;
}

template <typename Graph>
bool is_self_loop(graph_edge_descriptor_t<Graph> e, const Graph& g) {
  return source(e, g) == target(e, g);
}

// Functor for remove_parallel_edges: edge property of the removed edge is
// added to the remaining
template <typename EdgeProperty>
struct add_removed_edge_property {
  add_removed_edge_property(EdgeProperty ep) : ep(ep) {}

  template <typename Edge>
  void operator()(Edge stay, Edge away) {
    put(ep, stay, get(ep, stay) + get(ep, away));
  }
  EdgeProperty ep;
};

// Same as above: edge property is capacity here
template <typename Graph>
struct add_removed_edge_capacity : add_removed_edge_property<property_map_t<Graph, edge_capacity_t> > {
  using base = add_removed_edge_property<property_map_t<Graph, edge_capacity_t> >;
  add_removed_edge_capacity(Graph& g) : base(get(edge_capacity, g)) {}
};

template <typename Graph>
bool has_no_vertices(const Graph& g) {
  auto [a, b] = vertices(g);
  return (a == b);
}

template <typename Graph>
bool has_no_edges(const Graph& g) {
  auto [a, b] = edges(g);
  return (a == b);
}

template <typename Graph>
bool has_no_out_edges(const graph_vertex_descriptor_t<Graph>& v, const Graph& g) {
  auto [a, b] = out_edges(v, g);
  return (a == b);
}

template <typename PropertyIn, typename PropertyOut, typename Graph>
void copy_vertex_property(PropertyIn p_in, PropertyOut p_out, Graph& g) {
  for (auto u : vertices(g)) {
    put(p_out, u, get(p_in, g));
  }
}

template <typename PropertyIn, typename PropertyOut, typename Graph>
void copy_edge_property(PropertyIn p_in, PropertyOut p_out, Graph& g) {
  for (auto e : edges(g)) {
    put(p_out, e, get(p_in, g));
  }
}

// Return true if property_map1 and property_map2 differ
// for any of the vertices in graph.
template <typename PropertyMapFirst, typename PropertyMapSecond, typename Graph>
bool are_property_maps_different(const PropertyMapFirst property_map1, const PropertyMapSecond property_map2,
                                 const Graph& g) {
  for (auto u : vertices(g)) {
    if (get(property_map1, u) != get(property_map2, u)) {
      return true;
    }
  }

  return false;
}

}  // namespace bagl

#endif  // BAGL_BAGL_GRAPH_UTILITY_H_
