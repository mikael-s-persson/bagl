// Copyright 2007 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_IS_STRAIGHT_LINE_DRAWING_H_
#define BAGL_BAGL_IS_STRAIGHT_LINE_DRAWING_H_

#include <algorithm>
#include <map>
#include <set>
#include <tuple>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/planar_detail/bucket_sort.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

namespace bagl {

// Return true exactly when the line segments s1 = ((x1,y1), (x2,y2)) and
// s2 = ((a1,b1), (a2,b2)) intersect in a point other than the endpoints of
// the line segments. The one exception to this rule is when s1 = s2, in
// which case false is returned - this is to accomodate multiple edges
// between the same pair of vertices, which shouldn't invalidate the straight
// line embedding. A tolerance variable epsilon can also be used, which
// defines how far away from the endpoints of s1 and s2 we want to consider
// an intersection.
inline bool intersects(double x1, double y1, double x2, double y2, double a1, double b1, double a2, double b2,
                       double epsilon = 0.000001) {
  if (x1 - x2 == 0) {
    std::swap(x1, a1);
    std::swap(y1, b1);
    std::swap(x2, a2);
    std::swap(y2, b2);
  }

  if (x1 - x2 == 0) {
    // two vertical line segments
    double min_y = std::min(y1, y2);
    double max_y = std::max(y1, y2);
    double min_b = std::min(b1, b2);
    double max_b = std::max(b1, b2);
    return (max_y > max_b && max_b > min_y) || (max_b > max_y && max_y > min_b);
  }

  double x_diff = x1 - x2;
  double y_diff = y1 - y2;
  double a_diff = a2 - a1;
  double b_diff = b2 - b1;

  double beta_denominator = b_diff - (y_diff / x_diff) * a_diff;

  if (beta_denominator == 0) {
    // parallel lines
    return false;
  }

  double beta = (b2 - y2 - (y_diff / x_diff) * (a2 - x2)) / beta_denominator;
  double alpha = (a2 - x2 - beta * (a_diff)) / x_diff;

  double upper_bound = 1 - epsilon;
  double lower_bound = 0 + epsilon;

  return (beta < upper_bound && beta > lower_bound && alpha < upper_bound && alpha > lower_bound);
}

template <concepts::EdgeListGraph Graph, concepts::LvalueVertexPropertyMap<Graph> GridPositionMap>
bool is_straight_line_drawing(const Graph& g, GridPositionMap drawing) {
  using edge_t = graph_edge_descriptor_t<Graph>;

  using x_coord_t = std::size_t;
  using y_coord_t = std::size_t;
  using edge_event_t = std::tuple<edge_t, x_coord_t, y_coord_t>;
  using edge_event_queue_t = std::vector<edge_event_t>;

  using active_map_key_t = std::tuple<y_coord_t, y_coord_t, x_coord_t, x_coord_t>;
  using active_map_value_t = edge_t;
  using active_map_t = std::map<active_map_key_t, active_map_value_t>;

  edge_event_queue_t edge_event_queue;
  active_map_t active_edges;

  for (auto e : edges(g)) {
    auto s = source(e, g);
    auto t = target(e, g);
    edge_event_queue.push_back(
        make_tuple(e, static_cast<std::size_t>(drawing[s].x), static_cast<std::size_t>(drawing[s].y)));
    edge_event_queue.push_back(
        make_tuple(e, static_cast<std::size_t>(drawing[t].x), static_cast<std::size_t>(drawing[t].y)));
  }

  // Order by edge_event_queue by first, then second coordinate
  // (bucket_sort is a stable sort.)
  bucket_sort(edge_event_queue.begin(), edge_event_queue.end(), property_map_tuple_adaptor<edge_event_t, 2>());

  bucket_sort(edge_event_queue.begin(), edge_event_queue.end(), property_map_tuple_adaptor<edge_event_t, 1>());

  for (const auto& [e, x_c, y_c] : edge_event_queue) {
    auto source_v = source(e, g);
    auto target_v = target(e, g);
    if (drawing[source_v].y > drawing[target_v].y) {
      std::swap(source_v, target_v);
    }

    active_map_key_t key(get(drawing, source_v).y, get(drawing, target_v).y, get(drawing, source_v).x,
                         get(drawing, target_v).x);

    auto a_itr = active_edges.find(key);
    if (a_itr == active_edges.end()) {
      active_edges[key] = e;
    } else {
      auto before = active_edges.end();
      if (a_itr != active_edges.begin()) {
        before = std::prev(a_itr);
      }
      auto after = std::next(a_itr);

      if (before != active_edges.end()) {
        auto f = before->second;
        auto e_source = source(e, g);
        auto e_target = target(e, g);
        auto f_source = source(f, g);
        auto f_target = target(f, g);

        if (intersects(drawing[e_source].x, drawing[e_source].y, drawing[e_target].x, drawing[e_target].y,
                       drawing[f_source].x, drawing[f_source].y, drawing[f_target].x, drawing[f_target].y)) {
          return false;
        }
      }

      if (after != active_edges.end()) {
        auto f = after->second;
        auto e_source = source(e, g);
        auto e_target = target(e, g);
        auto f_source = source(f, g);
        auto f_target = target(f, g);

        if (intersects(drawing[e_source].x, drawing[e_source].y, drawing[e_target].x, drawing[e_target].y,
                       drawing[f_source].x, drawing[f_source].y, drawing[f_target].x, drawing[f_target].y)) {
          return false;
        }
      }

      active_edges.erase(a_itr);
    }
  }

  return true;
}

template <concepts::EdgeListGraph Graph, concepts::LvalueVertexPropertyMap<Graph> GridPositionMap,
          typename VertexIndexMap>
bool is_straight_line_drawing(const Graph& g, GridPositionMap drawing, VertexIndexMap /*unused*/) {
  return is_straight_line_drawing(g, drawing);
}

}  // namespace bagl

#endif  // BAGL_BAGL_IS_STRAIGHT_LINE_DRAWING_H_
