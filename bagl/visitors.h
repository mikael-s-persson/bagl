// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_VISITORS_H_
#define BAGL_BAGL_VISITORS_H_

#include <iosfwd>
#include <tuple>

#include "bagl/graph_traits.h"
#include "bagl/has_trait_member.h"
#include "bagl/property_map.h"

namespace bagl {

//========================================================================
// Event Tags

namespace visitors_detail {

BAGL_GRAPH_HAS_MEMBER_FUNCTION(initialize_vertex)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(start_vertex)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(visit_vertex)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(discover_vertex)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(cycle)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(examine_edge)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(examine_in_edge)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(examine_out_edge)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(examine_vertex)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(tree_edge)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(non_tree_edge)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(tree_in_edge)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(tree_out_edge)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(non_tree_in_edge)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(non_tree_out_edge)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(edge_relaxed)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(edge_not_relaxed)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(edge_minimized)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(edge_not_minimized)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(gray_target)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(black_target)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(gray_source)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(black_source)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(back_edge)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(forward_or_cross_edge)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(finish_edge)
BAGL_GRAPH_HAS_MEMBER_FUNCTION(finish_vertex)

}  // namespace visitors_detail

//========================================================================
// base_visitor and null_visitor

// needed for MSVC workaround
template <class Visitor>
struct base_visitor {
  template <class T, class Graph>
  void operator()(T /*unused*/, Graph& /*unused*/) {}
};

struct null_visitor : public base_visitor<null_visitor> {
  template <class T, class Graph>
  void operator()(T /*unused*/, Graph& /*unused*/) {}
};

using null_visitors = std::tuple<null_visitor>;

//========================================================================
// predecessor_recorder

template <typename PredecessorMap>
struct predecessor_recorder_on_tree_edge {
  explicit predecessor_recorder_on_tree_edge(PredecessorMap pa) : predecessor(std::move(pa)) {}
  template <class Edge, class Graph>
  void tree_edge(Edge e, const Graph& g) {
    put(predecessor, target(e, g), source(e, g));
  }
  PredecessorMap predecessor;
};

//========================================================================
// edge_predecessor_recorder

template <typename PredEdgeMap>
struct edge_predecessor_recorder_on_edge_relaxed {
  explicit edge_predecessor_recorder_on_edge_relaxed(PredEdgeMap pa) : predecessor(pa) {}
  template <class Edge, class Graph>
  void edge_relaxed(Edge e, const Graph& g) {
    put(predecessor, target(e, g), e);
  }
  PredEdgeMap predecessor;
};

template <typename PredEdgeMap>
struct edge_predecessor_recorder_on_tree_edge {
  explicit edge_predecessor_recorder_on_tree_edge(PredEdgeMap pa) : predecessor(pa) {}
  template <class Edge, class Graph>
  void tree_edge(Edge e, const Graph& g) {
    put(predecessor, target(e, g), e);
  }
  PredEdgeMap predecessor;
};

//========================================================================
// distance_recorder

template <class DistanceMap, class Tag>
struct distance_recorder : public base_visitor<distance_recorder<DistanceMap, Tag> > {
  explicit distance_recorder(DistanceMap pa) : m_distance(pa) {}
  template <class Edge, class Graph>
  void operator()(Edge e, const Graph& g) {
    auto u = source(e, g);
    auto v = target(e, g);
    put(m_distance, v, get(m_distance, u) + 1);
  }
  DistanceMap m_distance;
};
template <class DistanceMap, class Tag>
distance_recorder<DistanceMap, Tag> record_distances(DistanceMap pa, Tag) {
  return distance_recorder<DistanceMap, Tag>(pa);
}

//========================================================================
// time_stamper

template <class TimeMap, class TimeT, class Tag>
struct time_stamper : public base_visitor<time_stamper<TimeMap, TimeT, Tag> > {
  time_stamper(TimeMap pa, TimeT& t) : m_time_pa(pa), m_time(t) {}
  template <class Vertex, class Graph>
  void operator()(Vertex u, const Graph&) {
    put(m_time_pa, u, ++m_time);
  }
  TimeMap m_time_pa;
  TimeT& m_time;
};
template <class TimeMap, class TimeT, class Tag>
time_stamper<TimeMap, TimeT, Tag> stamp_times(TimeMap pa, TimeT& time_counter, Tag) {
  return time_stamper<TimeMap, TimeT, Tag>(pa, time_counter);
}

//========================================================================
// property_writer

template <class PA, class OutputIterator, class Tag>
struct property_writer : public base_visitor<property_writer<PA, OutputIterator, Tag> > {
  property_writer(PA pa, OutputIterator out) : m_pa(pa), m_out(out) {}

  template <class T, class Graph>
  void operator()(T x, Graph&) {
    *m_out++ = get(m_pa, x);
  }
  PA m_pa;
  OutputIterator m_out;
};
template <class PA, class OutputIterator, class Tag>
property_writer<PA, OutputIterator, Tag> write_property(PA pa, OutputIterator out, Tag) {
  return property_writer<PA, OutputIterator, Tag>(pa, out);
}

//========================================================================
// property_put

/**
 * Functor which just sets a given value to a vertex or edge in a property map.
 */

template <typename PropertyMap, typename EventTag>
struct property_put {
  property_put(PropertyMap property_map, typename property_traits<PropertyMap>::value_type value)
      : property_map_(property_map), value_(value) {}

  template <typename VertexOrEdge, typename Graph>
  void operator()(VertexOrEdge v, const Graph&) {
    put(property_map_, v, value_);
  }

 private:
  PropertyMap property_map_;
  typename property_traits<PropertyMap>::value_type value_;
};

/**
 * Creates a property_put functor which just sets a given value to a vertex or
 * edge.
 *
 * @param property_map Given writeable property map
 * @param value Fixed value of the map
 * @param tag Event Filter
 * @return The functor.
 */

template <typename PropertyMap, typename EventTag>
inline property_put<PropertyMap, EventTag> put_property(PropertyMap property_map,
                                                        typename property_traits<PropertyMap>::value_type value,
                                                        EventTag /*unused*/) {
  return property_put<PropertyMap, EventTag>(property_map, value);
}

}  // namespace bagl

#endif  // BAGL_BAGL_VISITORS_H_
