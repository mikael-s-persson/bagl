// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GRAPH_UTILITY_H_
#define BAGL_BAGL_GRAPH_UTILITY_H_

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
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
  }
  if (v == target(e, g)) {
    return source(e, g);
  }
  return vertex_descriptor{};
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

// This property map is a best-effort to "stringify" a vertex or edge descriptor.
// This is mostly meant to be used a default for `print_graph`.
template <typename G, typename Desc>
class printable_descriptor_property_map : public put_get_helper<printable_descriptor_property_map<G, Desc>> {
 public:
  using key_type = Desc;
  using value_type = std::string;
  using reference = value_type;

  explicit printable_descriptor_property_map(const G& g) : g_(&g) {}
  value_type operator[](const key_type& k) const {
    // Try:
    //  1) vertex/edge name maps
    //  2) vertex/edge index maps
    //  3) descriptor itself as an integral type
    //  4) address of property
    //  5) address of bundle
    if constexpr (std::is_same_v<Desc, graph_vertex_descriptor_t<G>> && has_property_map_v<G, vertex_name_t>) {
      const auto& vname = get(vertex_name, *g_, k);
      if constexpr (std::is_convertible_v<decltype(vname), std::string>) {
        return std::string{vname};
      } else {
        return std::to_string(vname);
      }
    } else if constexpr (std::is_same_v<Desc, graph_edge_descriptor_t<G>> && has_property_map_v<G, edge_name_t>) {
      const auto& ename = get(edge_name, *g_, k);
      if constexpr (std::is_convertible_v<decltype(ename), std::string>) {
        return std::string{ename};
      } else {
        return std::to_string(ename);
      }
    } else if constexpr (std::is_same_v<Desc, graph_vertex_descriptor_t<G>> && has_property_map_v<G, vertex_index_t>) {
      return std::to_string(get(vertex_index, *g_, k));
    } else if constexpr (std::is_same_v<Desc, graph_edge_descriptor_t<G>> && has_property_map_v<G, edge_index_t>) {
      return std::to_string(get(edge_index, *g_, k));
    } else if constexpr (std::is_integral_v<Desc>) {
      return std::to_string(k);
    } else if constexpr (std::is_same_v<Desc, graph_vertex_descriptor_t<G>> && has_property_map_v<G, vertex_all_t>) {
      return (std::stringstream{} << std::hex << reinterpret_cast<std::uintptr_t>(&get(vertex_all, *g_, k))).str();
    } else if constexpr (std::is_same_v<Desc, graph_edge_descriptor_t<G>> && has_property_map_v<G, edge_all_t>) {
      return (std::stringstream{} << std::hex << reinterpret_cast<std::uintptr_t>(&get(edge_all, *g_, k))).str();
    } else {
      return (std::stringstream{} << std::hex << reinterpret_cast<std::uintptr_t>(&(*g_)[k])).str();
    }
  }

 private:
  const G* g_;
};

template <concepts::IncidenceGraph G, concepts::ReadableVertexPropertyMap<G> VertexName,
          concepts::ReadableEdgePropertyMap<G> EdgeName>
void print_graph(const G& g, VertexName vname, EdgeName ename, std::ostream& os = std::cout) {
  const char* arrow = (is_directed_graph_v<G> ? " --> " : " <-> ");
  if constexpr (concepts::VertexListGraph<G>) {
    // Print vertices as "u --> v1 v2 v3.." on each line.
    os << "Vertices:\n";
    for (auto u : vertices(g)) {
      os << get(vname, u) << arrow;
      for (auto e : out_edges(u, g)) {
        os << get(vname, target(e, g)) << " ";
      }
      os << '\n';
    }
  }
  if constexpr (concepts::EdgeListGraph<G>) {
    // Print edges as "e ( u --> v )" on each line.
    os << "Edges:\n";
    for (auto e : edges(g)) {
      os << get(ename, e) << " ( " << get(vname, source(e, g)) << arrow << get(vname, target(e, g)) << " )\n";
    }
  }
}
template <concepts::IncidenceGraph G, concepts::ReadableVertexPropertyMap<G> VertexName>
void print_graph(const G& g, VertexName vname, std::ostream& os = std::cout) {
  print_graph(g, vname, printable_descriptor_property_map<G, graph_edge_descriptor_t<G>>(g), os);
}
template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> EdgeName>
void print_graph(const G& g, EdgeName ename, std::ostream& os = std::cout) {
  print_graph(g, printable_descriptor_property_map<G, graph_vertex_descriptor_t<G>>(g), ename, os);
}
template <concepts::IncidenceGraph G>
void print_graph(const G& g, std::ostream& os = std::cout) {
  print_graph(g, printable_descriptor_property_map<G, graph_vertex_descriptor_t<G>>(g),
              printable_descriptor_property_map<G, graph_edge_descriptor_t<G>>(g), os);
}

template <concepts::Graph G, typename Vertex>
bool is_adjacent(G& g, Vertex a, Vertex b) {
  if constexpr (concepts::AdjacencyGraph<G>) {
    auto adj_rg = adjacent_vertices(a, g);
    const auto adj_it = std::find(adj_rg.begin(), adj_rg.end(), b);
    if (adj_it == adj_rg.end()) {
      return false;
    }
  }

  if constexpr (concepts::IncidenceGraph<G>) {
    auto oe_rg = out_edges(a, g);
    const auto oe_it = std::find_if(oe_rg.begin(), oe_rg.end(), incident_to(b, g));
    if (oe_it == oe_rg.end()) {
      return false;
    }
  }

  if constexpr (is_bidirectional_graph_v<G>) {
    auto ie_rg = in_edges(b, g);
    const auto ie_it = std::find_if(ie_rg.begin(), ie_rg.end(), incident_from(a, g));
    if (ie_it == ie_rg.end()) {
      return false;
    }
  }

  return true;
}

template <concepts::Graph G, typename Vertex>
auto num_adjacent_vertices(G& g, Vertex a) {
  if constexpr (concepts::AdjacencyGraph<G>) {
    return std::ranges::distance(adjacent_vertices(a, g));
  } else if constexpr (concepts::IncidenceGraph<G>) {
    return std::ranges::distance(out_edges(a, g));
  }
}

template <concepts::Graph G, typename Edge>
bool in_edge_set(G& g, Edge e) {
  if constexpr (concepts::EdgeListGraph<G>) {
    auto e_rg = edges(g);
    return std::find(e_rg.begin(), e_rg.end(), e) != e_rg.end();
  } else if constexpr (concepts::VertexListGraph<G> && concepts::IncidenceGraph<G>) {
    for (auto u : vertices(g)) {
      auto oe_rg = out_edges(u, g);
      if (std::find(oe_rg.begin(), oe_rg.end(), e) != oe_rg.end()) {
        return true;
      }
    }
    return false;
  } else {
    // Not in edge set, because there is no edge set to inspect.
    return false;
  }
}

template <concepts::Graph G, typename Vertex>
bool in_vertex_set(G& g, Vertex v) {
  if constexpr (concepts::VertexListGraph<G>) {
    auto v_rg = vertices(g);
    return std::find(v_rg.begin(), v_rg.end(), v) != v_rg.end();
  } else if constexpr (concepts::EdgeListGraph<G> && concepts::IncidenceGraph<G>) {
    for (auto e : edges(g)) {
      if (target(e, g) == v || source(e, g) == v) {
        return true;
      }
    }
    return false;
  } else {
    // Not in vertex set, because there is no vertex set to inspect.
    return false;
  }
}

template <concepts::IncidenceGraph G, typename Vertex>
bool in_edge_set(G& g, Vertex u, Vertex v) {
  if constexpr (concepts::EdgeListGraph<G>) {
    for (auto e : edges(g)) {
      if (connects(u, v, g)(e)) {
        return true;
      }
    }
    return false;
  } else if constexpr (concepts::VertexListGraph<G>) {
    // We don't just check out_edges(u,g) because this is used in
    // tests to exhaustively check the graph.
    for (auto a : vertices(g)) {
      auto oe_rg = out_edges(a, g);
      if (std::find_if(oe_rg.begin(), oe_rg.end(), connects(u, v, g)) != oe_rg.end()) {
        return true;
      }
    }
    return false;
  } else {
    // Not in edge set, because there is no edge set to inspect.
    return false;
  }
}

// is x a descendant of y?
template <typename ParentMap>
bool is_descendant(property_traits_value_t<ParentMap> x, property_traits_value_t<ParentMap> y, ParentMap parent) {
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
bool is_reachable(graph_vertex_descriptor_t<IncidenceGraph> x, graph_vertex_descriptor_t<IncidenceGraph> y,
                  const IncidenceGraph& g,
                  VertexColorMap color)  // should start out white for every vertex
{
  using ColorValue = property_traits_value_t<VertexColorMap>;
  dfs_visitor<> vis;
  depth_first_visit(g, x, vis, color);
  return get(color, y) != color_traits<ColorValue>::white();
}

// Is the undirected graph connected?
// Is the directed graph strongly connected?
template <typename VertexListGraph, typename VertexColorMap>
bool is_connected(const VertexListGraph& g, VertexColorMap color) {
  using ColorValue = property_traits_value_t<VertexColorMap>;
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
  explicit add_removed_edge_property(EdgeProperty a_ep) : ep(a_ep) {}

  template <typename Edge>
  void operator()(Edge stay, Edge away) {
    put(ep, stay, get(ep, stay) + get(ep, away));
  }
  EdgeProperty ep;
};

// Same as above: edge property is capacity here
template <typename Graph>
struct add_removed_edge_capacity : add_removed_edge_property<property_map_t<Graph, edge_capacity_t>> {
  using base = add_removed_edge_property<property_map_t<Graph, edge_capacity_t>>;
  explicit add_removed_edge_capacity(Graph& g) : base(get(edge_capacity, g)) {}
};

template <typename G>
auto num_vertices_or_zero(const G& g) {
  if constexpr (concepts::VertexListGraph<G>) {
    return num_vertices(g);
  } else {
    return 0;
  }
}

template <typename G>
auto num_edges_or_zero(const G& g) {
  if constexpr (concepts::EdgeListGraph<G>) {
    return num_edges(g);
  } else {
    return 0;
  }
}

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
