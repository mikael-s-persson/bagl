// Copyright 2004, 2005 The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_FRUCHTERMAN_REINGOLD_H_
#define BAGL_BAGL_FRUCHTERMAN_REINGOLD_H_

#include <algorithm>  // for std::min and std::max
#include <cmath>
#include <functional>
#include <list>
#include <numeric>  // for std::accumulate
#include <vector>

#include "bagl/filtered_graph.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"
#include "bagl/topology.h"  // For topology concepts
#include "bagl/vector_property_map.h"

namespace bagl {

struct square_distance_attractive_force {
  template <typename Graph, typename T>
  T operator()(graph_edge_descriptor_t<Graph> /*unused*/, T k, T d, const Graph& /*unused*/) const {
    return d * d / k;
  }
};

struct square_distance_repulsive_force {
  template <typename Graph, typename T>
  T operator()(graph_vertex_descriptor_t<Graph> /*unused*/,
               graph_vertex_descriptor_t<Graph> /*unused*/, T k, T d, const Graph& /*unused*/) const {
    return k * k / d;
  }
};

template <typename T>
struct linear_cooling {
  explicit linear_cooling(std::size_t iterations) : temp_(T(iterations) / T(10)), step_(0.1) {}

  linear_cooling(std::size_t iterations, T temp) : temp_(temp), step_(temp / T(iterations)) {}

  T operator()() {
    T old_temp = temp_;
    temp_ -= step_;
    if (temp_ < T(0)) {
      temp_ = T(0);
    }
    return old_temp;
  }

 private:
  T temp_;
  T step_;
};

struct all_force_pairs {
  template <typename Graph, typename ApplyForce>
  void operator()(const Graph& g, ApplyForce apply_force) {
    for (auto v_rg = partial_view(vertices(g)); !v_rg.empty(); v_rg.move_to_next()) {
      auto u_rg = v_rg;
      for (u_rg.move_to_next(); !u_rg.empty(); u_rg.move_to_next()) {
        apply_force(*u_rg.begin(), *v_rg.begin());
        apply_force(*v_rg.begin(), *u_rg.begin());
      }
    }
  }
};

template <typename Topology, typename PositionMap>
struct grid_force_pairs {
  using Point = typename property_traits<PositionMap>::value_type;
  static_assert(Point::dimensions == 2);

  template <typename Graph>
  explicit grid_force_pairs(const Topology& topology, PositionMap position, const Graph& g)
      : topology_(topology),
        position_(position),
        two_k_(2. * topology_.volume(topology_.extent()) / std::sqrt(static_cast<double>(num_vertices(g)))) {}

  template <typename Graph, typename ApplyForce>
  void operator()(const Graph& g, ApplyForce apply_force) {
    using Vertex = typename graph_traits<Graph>::vertex_descriptor;
    using Bucket = std::list<Vertex>;
    using Buckets = std::vector<Bucket>;

    std::size_t columns = std::size_t(topology_.extent()[0] / two_k_ + 1.);
    std::size_t rows = std::size_t(topology_.extent()[1] / two_k_ + 1.);
    Buckets buckets(rows * columns);
    for (auto v : vertices(g)) {
      std::size_t column = std::size_t((get(position_, v)[0] + topology_.extent()[0] / 2) / two_k_);
      std::size_t row = std::size_t((get(position_, v)[1] + topology_.extent()[1] / 2) / two_k_);

      if (column >= columns) {
        column = columns - 1;
      }
      if (row >= rows) {
        row = rows - 1;
      }
      buckets[row * columns + column].push_back(v);
    }

    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t column = 0; column < columns; ++column) {
        Bucket& bucket = buckets[row * columns + column];
        for (auto u = bucket.begin(); u != bucket.end(); ++u) {
          // Repulse vertices in this bucket
          auto v = u;
          for (++v; v != bucket.end(); ++v) {
            apply_force(*u, *v);
            apply_force(*v, *u);
          }

          std::size_t adj_start_row = row == 0 ? 0 : row - 1;
          std::size_t adj_end_row = row == rows - 1 ? row : row + 1;
          std::size_t adj_start_column = column == 0 ? 0 : column - 1;
          std::size_t adj_end_column = column == columns - 1 ? column : column + 1;
          for (std::size_t other_row = adj_start_row; other_row <= adj_end_row; ++other_row) {
            for (std::size_t other_column = adj_start_column; other_column <= adj_end_column; ++other_column) {
              if (other_row != row || other_column != column) {
                // Repulse vertices in this bucket
                Bucket& other_bucket = buckets[other_row * columns + other_column];
                for (v = other_bucket.begin(); v != other_bucket.end(); ++v) {
                  double dist = topology_.distance(get(position_, *u), get(position_, *v));
                  if (dist < two_k_) {
                    apply_force(*u, *v);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

 private:
  const Topology& topology_;
  PositionMap position_;
  double two_k_;
};

template <typename PositionMap, typename Topology, typename Graph>
auto make_grid_force_pairs(const Topology& topology, const PositionMap& position, const Graph& g) {
  return grid_force_pairs<Topology, PositionMap>(topology, position, g);
}

template <typename Graph, typename PositionMap, typename Topology>
void scale_graph(const Graph& g, PositionMap position, const Topology& topology,
                 typename Topology::point_type upper_left, typename Topology::point_type lower_right) {
  if (num_vertices(g) == 0) {
    return;
  }

  using Point = typename Topology::point_type;
  using PointDiff = typename Topology::point_difference_type;

  // Find min/max ranges
  Point min_point = get(position, *vertices(g).first);
  Point max_point = min_point;
  for (auto v : vertices(g)) {
    min_point = topology.pointwise_min(min_point, get(position, v));
    max_point = topology.pointwise_max(max_point, get(position, v));
  }

  Point old_origin = topology.move_position_toward(min_point, 0.5, max_point);
  Point new_origin = topology.move_position_toward(upper_left, 0.5, lower_right);
  PointDiff old_size = topology.difference(max_point, min_point);
  PointDiff new_size = topology.difference(lower_right, upper_left);

  // Scale to bounding box provided
  for (auto v : vertices(g)) {
    PointDiff relative_loc = topology.difference(get(position, v), old_origin);
    relative_loc = (relative_loc / old_size) * new_size;
    put(position, v, topology.adjust(new_origin, relative_loc));
  }
}

namespace fr_detail {

template <typename Topology, typename PropMap, typename Vertex>
void maybe_jitter_point(const Topology& topology, const PropMap& pm, Vertex v,
                        const typename Topology::point_type& p2) {
  double too_close = topology.norm(topology.extent()) / 10000.;
  if (topology.distance(get(pm, v), p2) < too_close) {
    put(pm, v, topology.move_position_toward(get(pm, v), 1. / 200, topology.random_point()));
  }
}

template <typename Topology, typename PositionMap, typename DisplacementMap, typename RepulsiveForce, typename Graph>
struct fr_apply_force {
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Point = typename Topology::point_type;
  using PointDiff = typename Topology::point_difference_type;

  fr_apply_force(const Topology& topology, const PositionMap& position, const DisplacementMap& displacement,
                 RepulsiveForce repulsive_force, double k, const Graph& g)
      : topology_(topology),
        position_(position),
        displacement_(displacement),
        repulsive_force_(repulsive_force),
        k_(k),
        g_(g) {}

  void operator()(Vertex u, Vertex v) {
    if (u != v) {
      // When the vertices land on top of each other, move the
      // first vertex away from the boundaries.
      maybe_jitter_point(topology_, position_, u, get(position_, v));

      double dist = topology_.distance(get(position_, u), get(position_, v));
      typename Topology::point_difference_type dispv = get(displacement_, v);
      if (dist == 0.) {
        for (std::size_t i = 0; i < Point::dimensions; ++i) {
          dispv[i] += 0.01;
        }
      } else {
        double fr = repulsive_force_(u, v, k_, dist, g_);
        dispv += (fr / dist) * topology_.difference(get(position_, v), get(position_, u));
      }
      put(displacement_, v, dispv);
    }
  }

 private:
  const Topology& topology_;
  PositionMap position_;
  DisplacementMap displacement_;
  RepulsiveForce repulsive_force_;
  double k_;
  const Graph& g_;
};

}  // namespace fr_detail

template <typename Topology, concepts::VertexAndEdgeListGraph G, concepts::ReadWriteVertexPropertyMap<G> PositionMap,
          typename AttractiveForce, typename RepulsiveForce, typename ForcePairs, typename Cooling,
          concepts::ReadWriteVertexPropertyMap<G> DisplacementMap>
void fruchterman_reingold_force_directed_layout(const G& g, PositionMap position, const Topology& topology,
                                                AttractiveForce attractive_force,  // square_distance_attractive_force()
                                                RepulsiveForce repulsive_force,    // square_distance_repulsive_force()
                                                ForcePairs force_pairs,  // make_grid_force_pairs(topology, position, g)
                                                Cooling cool,            // linear_cooling<double>(100)
                                                DisplacementMap displacement) {
  using Vertex = graph_vertex_descriptor_t<G>;
  using PointDiff = typename Topology::point_difference_type;

  double volume = topology.volume(topology.extent());

  // assume positions are initialized randomly
  double k = std::pow(volume / num_vertices(g), 1.0 / static_cast<double>(PointDiff::dimensions));

  fr_detail::fr_apply_force apply_force{topology, position, displacement, repulsive_force, k, g};

  do {
    // Calculate repulsive forces
    for (auto v : vertices(g)) {
      put(displacement, v, PointDiff{});
    }
    force_pairs(g, apply_force);

    // Calculate attractive forces
    for (auto e : edges(g)) {
      Vertex v = source(e, g);
      Vertex u = target(e, g);

      // When the vertices land on top of each other, move the
      // first vertex away from the boundaries.
      fr_detail::maybe_jitter_point(topology, position, u, get(position, v));

      PointDiff delta = topology.difference(get(position, v), get(position, u));
      double dist = topology.distance(get(position, u), get(position, v));
      double fa = attractive_force(e, k, dist, g);

      put(displacement, v, get(displacement, v) - (fa / dist) * delta);
      put(displacement, u, get(displacement, u) + (fa / dist) * delta);
    }

    if (double temp = cool(); temp != 0.0) {
      // Update positions
      for (auto v : vertices(g)) {
        double disp_size = topology.norm(get(displacement, v));
        put(position, v,
            topology.adjust(get(position, v), get(displacement, v) * (std::min(disp_size, temp) / disp_size)));
        put(position, v, topology.bound(get(position, v)));
      }
    } else {
      break;
    }
  } while (true);
}

template <typename Topology, concepts::VertexAndEdgeListGraph G, concepts::ReadWriteVertexPropertyMap<G> PositionMap,
          typename AttractiveForce, typename RepulsiveForce, typename ForcePairs, typename Cooling>
void fruchterman_reingold_force_directed_layout(const G& g, PositionMap position, const Topology& topology,
                                                AttractiveForce attractive_force,  // square_distance_attractive_force()
                                                RepulsiveForce repulsive_force,    // square_distance_repulsive_force()
                                                ForcePairs force_pairs,  // make_grid_force_pairs(topology, position, g)
                                                Cooling cool) {          // linear_cooling<double>(100)
  using PointDiff = typename Topology::point_difference_type;
  fruchterman_reingold_force_directed_layout(g, position, topology, attractive_force, repulsive_force, force_pairs,
                                             cool,
                                             vector_property_map(num_vertices(g), get(vertex_index, g), PointDiff{}));
}

}  // namespace bagl

#endif  // BAGL_BAGL_FRUCHTERMAN_REINGOLD_H_
