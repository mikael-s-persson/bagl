// Copyright 2004 The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_KAMADA_KAWAI_SPRING_LAYOUT_H_
#define BAGL_BAGL_KAMADA_KAWAI_SPRING_LAYOUT_H_

#include <cmath>
#include <iostream>
#include <iterator>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "bagl/graph_traits.h"
#include "bagl/johnson_all_pairs_shortest.h"
#include "bagl/topology.h"

namespace bagl {
namespace kamada_kawai_detail {
/**
 * Denotes an edge or display area side length used to scale a
 * Kamada-Kawai drawing.
 */
template <bool Edge, typename T>
struct edge_or_side {
  T value;
};

/**
 * Compute the edge length from an edge length. This is trivial.
 */
template <typename Graph, typename DistanceMap, typename IndexMap, typename T>
T compute_edge_length(const Graph&, DistanceMap, IndexMap, edge_or_side<true, T> length) {
  return length.value;
}

/**
 * Compute the edge length based on the display area side
   length. We do this by dividing the side length by the largest
   shortest distance between any two vertices in the graph.
 */
template <typename Graph, typename DistanceMap, typename IndexMap, typename T>
T compute_edge_length(const Graph& g, DistanceMap distance, IndexMap index, edge_or_side<false, T> length) {
  T result(0);
  auto v_rg = vertices(g);
  for (auto ui = v_rg.begin(); ui != v_rg.end(); ++ui) {
    auto vi = ui;
    for (++vi; vi != v_rg.end(); ++vi) {
      T dij = distance[get(index, *ui)][get(index, *vi)];
      if (dij > result) {
        result = dij;
      }
    }
  }
  return length.value / result;
}

/**
 * Dense linear solver for fixed-size matrices.
 */
template <std::size_t Size>
struct linear_solver {
  // Indices in mat are (row, column)
  // template <typename Vec>
  // static Vec solve(double mat[Size][Size], Vec rhs);
};

template <>
struct linear_solver<1> {
  template <typename Vec>
  static Vec solve(const std::array<std::array<double, 1>, 1>& mat, Vec rhs) {
    return rhs / mat[0][0];
  }
};

// These are from http://en.wikipedia.org/wiki/Cramer%27s_rule
template <>
struct linear_solver<2> {
  template <typename Vec>
  static Vec solve(const std::array<std::array<double, 2>, 2>& mat, Vec rhs) {
    double denom = mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
    double x_num = rhs[0] * mat[1][1] - rhs[1] * mat[0][1];
    double y_num = mat[0][0] * rhs[1] - mat[1][0] * rhs[0];
    Vec result;
    result[0] = x_num / denom;
    result[1] = y_num / denom;
    return result;
  }
};

template <>
struct linear_solver<3> {
  template <typename Vec>
  static Vec solve(const std::array<std::array<double, 3>, 3>& mat, Vec rhs) {
    double denom = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) -
                   mat[1][0] * (mat[0][1] * mat[2][2] - mat[2][1] * mat[0][2]) +
                   mat[2][0] * (mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2]);
    double x_num = rhs[0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) -
                   rhs[1] * (mat[0][1] * mat[2][2] - mat[2][1] * mat[0][2]) +
                   rhs[2] * (mat[0][1] * mat[1][2] - mat[1][1] * mat[0][2]);
    double y_num = mat[0][0] * (rhs[1] * mat[2][2] - rhs[2] * mat[1][2]) -
                   mat[1][0] * (rhs[0] * mat[2][2] - rhs[2] * mat[0][2]) +
                   mat[2][0] * (rhs[0] * mat[1][2] - rhs[1] * mat[0][2]);
    double z_num = mat[0][0] * (mat[1][1] * rhs[2] - mat[2][1] * rhs[1]) -
                   mat[1][0] * (mat[0][1] * rhs[2] - mat[2][1] * rhs[0]) +
                   mat[2][0] * (mat[0][1] * rhs[1] - mat[1][1] * rhs[0]);
    Vec result;
    result[0] = x_num / denom;
    result[1] = y_num / denom;
    result[2] = z_num / denom;
    return result;
  }
};

/**
 * Implementation of the Kamada-Kawai spring layout algorithm.
 */
template <typename Topology, concepts::VertexAndEdgeListGraph Graph,
          concepts::LvalueVertexPropertyMap<Graph> PositionMap, concepts::ReadableEdgePropertyMap<Graph> WeightMap,
          typename EdgeOrSideLength, typename Done, concepts::ReadableVertexIndexMap<Graph> VertexIndexMap,
          concepts::BasicMatrix<property_traits_value_t<WeightMap>> DistanceMatrix,
          concepts::BasicMatrix<property_traits_value_t<WeightMap>> SpringStrengthMatrix,
          concepts::ReadWriteVertexPropertyMap<Graph> PartialDerivativeMap>
class kamada_kawai_spring_layout_impl {
 public:
  using weight_type = property_traits_value_t<WeightMap>;
  using Point = typename Topology::point_type;
  using point_difference_type = typename Topology::point_difference_type;
  using deriv_type = point_difference_type;
  using vertex_descriptor = graph_vertex_descriptor_t<Graph>;

  kamada_kawai_spring_layout_impl(const Topology& topology, const Graph& g, PositionMap position, WeightMap weight,
                                  EdgeOrSideLength edge_or_side_length, Done done, weight_type spring_constant,
                                  VertexIndexMap index, DistanceMatrix distance, SpringStrengthMatrix spring_strength,
                                  PartialDerivativeMap partial_derivatives)
      : topology_(topology),
        g_(g),
        position_(position),
        weight_(weight),
        edge_or_side_length_(edge_or_side_length),
        done_(done),
        spring_constant_(spring_constant),
        index_(index),
        distance_(distance),
        spring_strength_(spring_strength),
        partial_derivatives_(partial_derivatives) {}

  // Compute contribution of vertex i to the first partial
  // derivatives (dE/dx_m, dE/dy_m) (for vertex m)
  deriv_type compute_partial_derivative(vertex_descriptor m, vertex_descriptor i) {
    using std::sqrt;

    deriv_type result;
    if (i != m) {
      point_difference_type diff = topology_.difference(position_[m], position_[i]);
      weight_type dist = topology_.norm(diff);
      result = spring_strength_[get(index_, m)][get(index_, i)] *
               (diff - distance_[get(index_, m)][get(index_, i)] / dist * diff);
    }

    return result;
  }

  // Compute partial derivatives dE/dx_m and dE/dy_m
  deriv_type compute_partial_derivatives(vertex_descriptor m) {
    using std::sqrt;

    deriv_type result;

    // TBD: looks like an accumulate to me
    for (auto i : vertices(g_)) {
      deriv_type deriv = compute_partial_derivative(m, i);
      result += deriv;
    }

    return result;
  }

  // The actual Kamada-Kawai spring layout algorithm implementation
  bool run() {
    using std::sqrt;

    // Compute d_{ij} and place it in the distance matrix
    if (!johnson_all_pairs_shortest_paths(g_, distance_, index_, weight_)) {
      return false;
    }

    // Compute L based on side length (if needed), or retrieve L
    weight_type edge_length = compute_edge_length(g_, distance_, index_, edge_or_side_length_);

    // Compute l_{ij} and k_{ij}
    const weight_type k = spring_constant_;
    auto v_rg = vertices(g_);
    for (auto ui = v_rg.begin(); ui != v_rg.end(); ++ui) {
      auto vi = ui;
      for (++vi; vi != v_rg.end(); ++vi) {
        weight_type dij = distance_[get(index_, *ui)][get(index_, *vi)];
        if (dij == std::numeric_limits<weight_type>::max()) {
          return false;
        }
        distance_[get(index_, *ui)][get(index_, *vi)] = edge_length * dij;
        distance_[get(index_, *vi)][get(index_, *ui)] = edge_length * dij;
        spring_strength_[get(index_, *ui)][get(index_, *vi)] = k / (dij * dij);
        spring_strength_[get(index_, *vi)][get(index_, *ui)] = k / (dij * dij);
      }
    }

    // Compute Delta_i and find max
    vertex_descriptor p = *vertices(g_).begin();
    weight_type delta_p(0);

    for (auto u : v_rg) {
      deriv_type deriv = compute_partial_derivatives(u);
      put(partial_derivatives_, u, deriv);

      weight_type delta = topology_.norm(deriv);

      if (delta > delta_p) {
        p = u;
        delta_p = delta;
      }
    }

    while (!done_(delta_p, p, g_, true)) {
      // The contribution p makes to the partial derivatives of
      // each vertex. Computing this (at O(n) cost) allows us to
      // update the delta_i values in O(n) time instead of O(n^2)
      // time.
      std::vector<deriv_type> p_partials(num_vertices(g_));
      for (auto u : v_rg) {
        p_partials[get(index_, u)] = compute_partial_derivative(u, p);
      }

      do {
        // For debugging, compute the energy value E
        double e = 0.;
        for (auto ui = v_rg.begin(); ui != v_rg.end(); ++ui) {
          auto vi = ui;
          for (++vi; vi != v_rg.end(); ++vi) {
            double dist = topology_.distance(position_[*ui], position_[*vi]);
            weight_type k_ij = spring_strength_[get(index_, *ui)][get(index_, *vi)];
            weight_type l_ij = distance_[get(index_, *ui)][get(index_, *vi)];
            e += .5 * k_ij * (dist - l_ij) * (dist - l_ij);
          }
        }
        // std::cerr << "E = " << e << std::endl;

        // Compute the elements of the Jacobian
        // From
        // http://www.cs.panam.edu/~rfowler/papers/1994_kumar_fowler_A_Spring_UTPACSTR.pdf
        // with the bugs fixed in the off-diagonal case
        std::array<std::array<weight_type, Point::dimensions>, Point::dimensions> d_e_d_d{};
        for (std::size_t i = 0; i < Point::dimensions; ++i) {
          for (std::size_t j = 0; j < Point::dimensions; ++j) {
            d_e_d_d[i][j] = 0.;
          }
        }
        for (auto ui = v_rg.begin(); ui != v_rg.end(); ++ui) {
          auto u = *ui;
          if (u != p) {
            point_difference_type diff = topology_.difference(position_[p], position_[u]);
            weight_type dist = topology_.norm(diff);
            weight_type dist_squared = dist * dist;
            weight_type inv_dist_cubed = 1. / (dist_squared * dist);
            weight_type k_mi = spring_strength_[get(index_, p)][get(index_, u)];
            weight_type l_mi = distance_[get(index_, p)][get(index_, u)];
            for (std::size_t i = 0; i < Point::dimensions; ++i) {
              for (std::size_t j = 0; j < Point::dimensions; ++j) {
                if (i == j) {
                  d_e_d_d[i][i] += k_mi * (1 + (l_mi * (diff[i] * diff[i] - dist_squared) * inv_dist_cubed));
                } else {
                  d_e_d_d[i][j] += k_mi * l_mi * diff[i] * diff[j] * inv_dist_cubed;
                  // dE_d_d[i][j] += k_mi * l_mi *
                  // sqrt(hypot(diff[i], diff[j])) *
                  // inv_dist_cubed;
                }
              }
            }
          }
        }

        deriv_type d_e_d = get(partial_derivatives_, p);

        // Solve dE_d_d * delta = -dE_d to get delta
        point_difference_type delta = -linear_solver<Point::dimensions>::solve(d_e_d_d, d_e_d);

        // Move p by delta
        position_[p] = topology_.adjust(position_[p], delta);

        // Recompute partial derivatives and delta_p
        deriv_type deriv = compute_partial_derivatives(p);
        put(partial_derivatives_, p, deriv);

        delta_p = topology_.norm(deriv);
      } while (!done_(delta_p, p, g_, false));

      // Select new p by updating each partial derivative and
      // delta
      vertex_descriptor old_p = p;
      for (auto ui = v_rg.begin(); ui != v_rg.end(); ++ui) {
        deriv_type old_deriv_p = p_partials[get(index_, *ui)];
        deriv_type old_p_partial = compute_partial_derivative(*ui, old_p);
        deriv_type deriv = get(partial_derivatives_, *ui);

        deriv += old_p_partial - old_deriv_p;

        put(partial_derivatives_, *ui, deriv);
        weight_type delta = topology_.norm(deriv);

        if (delta > delta_p) {
          p = *ui;
          delta_p = delta;
        }
      }
    }

    return true;
  }

 private:
  const Topology& topology_;
  const Graph& g_;
  PositionMap position_;
  WeightMap weight_;
  EdgeOrSideLength edge_or_side_length_;
  Done done_;
  weight_type spring_constant_;
  VertexIndexMap index_;
  DistanceMatrix distance_;
  SpringStrengthMatrix spring_strength_;
  PartialDerivativeMap partial_derivatives_;
};
}  // namespace kamada_kawai_detail

/// States that the given quantity is an edge length.
template <typename T>
auto edge_length(T x) {
  return kamada_kawai_detail::edge_or_side<true, T>(x);
}

/// States that the given quantity is a display area side length.
template <typename T>
auto side_length(T x) {
  return kamada_kawai_detail::edge_or_side<false, T>(x);
}

/**
 * \brief Determines when to terminate layout of a particular graph based
 * on a given relative tolerance.
 */
template <typename T = double>
struct layout_tolerance {
  explicit layout_tolerance(const T& tolerance = T(0.001)) : tolerance_(tolerance) {}

  template <typename Graph>
  bool operator()(T delta_p, graph_vertex_descriptor_t<Graph> /*p*/, const Graph& /*g*/, bool global) {
    if (global) {
      if (last_energy_ == (std::numeric_limits<T>::max)()) {
        last_energy_ = delta_p;
        return false;
      }

      T diff = last_energy_ - delta_p;
      if (diff < T(0)) {
        diff = -diff;
      }
      bool done = (delta_p == T(0) || diff / last_energy_ < tolerance_);
      last_energy_ = delta_p;
      return done;
    }

    if (last_local_energy_ == (std::numeric_limits<T>::max)()) {
      last_local_energy_ = delta_p;
      return delta_p == T(0);
    }

    T diff = last_local_energy_ - delta_p;
    bool done = (delta_p == T(0) || (diff / last_local_energy_) < tolerance_);
    last_local_energy_ = delta_p;
    return done;
  }

 private:
  T tolerance_;
  T last_energy_ = std::numeric_limits<T>::max();
  T last_local_energy_ = std::numeric_limits<T>::max();
};

/** \brief Kamada-Kawai spring layout for undirected graphs.
 *
 * This algorithm performs graph layout (in two dimensions) for
 * connected, undirected graphs. It operates by relating the layout
 * of graphs to a dynamic spring system and minimizing the energy
 * within that system. The strength of a spring between two vertices
 * is inversely proportional to the square of the shortest distance
 * (in graph terms) between those two vertices. Essentially,
 * vertices that are closer in the graph-theoretic sense (i.e., by
 * following edges) will have stronger springs and will therefore be
 * placed closer together.
 *
 * Prior to invoking this algorithm, it is recommended that the
 * vertices be placed along the vertices of a regular n-sided
 * polygon.
 *
 * \param g (IN) must be a model of Vertex List Graph, Edge List
 * Graph, and Incidence Graph and must be undirected.
 *
 * \param position (OUT) must be a model of Lvalue Property Map,
 * where the value type is a class containing fields @c x and @c y
 * that will be set to the @c x and @c y coordinates of each vertex.
 *
 * \param weight (IN) must be a model of Readable Property Map,
 * which provides the weight of each edge in the graph @p g.
 *
 * \param topology (IN) must be a topology object (see topology.hpp),
 * which provides operations on points and differences between them.
 *
 * \param edge_or_side_length (IN) provides either the unit length
 * @c e of an edge in the layout or the length of a side @c s of the
 * display area, and must be either @c boost::edge_length(e) or @c
 * boost::side_length(s), respectively.
 *
 * \param done (IN) is a 4-argument function object that is passed
 * the current value of delta_p (i.e., the energy of vertex @p p),
 * the vertex @p p, the graph @p g, and a boolean flag indicating
 * whether @p delta_p is the maximum energy in the system (when @c
 * true) or the energy of the vertex being moved. Defaults to @c
 * layout_tolerance instantiated over the value type of the weight
 * map.
 *
 * \param spring_constant (IN) is the constant multiplied by each
 * spring's strength. Larger values create systems with more energy
 * that can take longer to stabilize; smaller values create systems
 * with less energy that stabilize quickly but do not necessarily
 * result in pleasing layouts. The default value is 1.
 *
 * \param index (IN) is a mapping from vertices to index values
 * between 0 and @c num_vertices(g). The default is @c
 * get(vertex_index,g).
 *
 * \param distance (UTIL/OUT) will be used to store the distance
 * from every vertex to every other vertex, which is computed in the
 * first stages of the algorithm. This value's type must be a model
 * of BasicMatrix with value type equal to the value type of the
 * weight map. The default is a vector of vectors.
 *
 * \param spring_strength (UTIL/OUT) will be used to store the
 * strength of the spring between every pair of vertices. This
 * value's type must be a model of BasicMatrix with value type equal
 * to the value type of the weight map. The default is a vector of
 * vectors.
 *
 * \param partial_derivatives (UTIL) will be used to store the
 * partial derivates of each vertex with respect to the @c x and @c
 * y coordinates. This must be a Read/Write Property Map whose value
 * type is a pair with both types equivalent to the value type of
 * the weight map. The default is an iterator property map.
 *
 * \returns @c true if layout was successful or @c false if a
 * negative weight cycle was detected.
 */
template <typename Topology, concepts::VertexAndEdgeListGraph Graph,
          concepts::LvalueVertexPropertyMap<Graph> PositionMap, concepts::ReadableEdgePropertyMap<Graph> WeightMap,
          typename T, bool EdgeOrSideLength, typename Done, concepts::ReadableVertexIndexMap<Graph> VertexIndexMap,
          concepts::BasicMatrix<property_traits_value_t<WeightMap>> DistanceMatrix,
          concepts::BasicMatrix<property_traits_value_t<WeightMap>> SpringStrengthMatrix,
          concepts::ReadWriteVertexPropertyMap<Graph> PartialDerivativeMap>
bool kamada_kawai_spring_layout(const Graph& g, PositionMap position, WeightMap weight, const Topology& topology,
                                kamada_kawai_detail::edge_or_side<EdgeOrSideLength, T> edge_or_side_length, Done done,
                                property_traits_value_t<WeightMap> spring_constant, VertexIndexMap index,
                                DistanceMatrix distance, SpringStrengthMatrix spring_strength,
                                PartialDerivativeMap partial_derivatives) {
  static_assert(is_undirected_graph_v<Graph>);

  kamada_kawai_detail::kamada_kawai_spring_layout_impl<
      Topology, Graph, PositionMap, WeightMap, kamada_kawai_detail::edge_or_side<EdgeOrSideLength, T>, Done,
      VertexIndexMap, DistanceMatrix, SpringStrengthMatrix, PartialDerivativeMap>
      alg(topology, g, position, weight, edge_or_side_length, done, spring_constant, index, distance, spring_strength,
          partial_derivatives);
  return alg.run();
}

/**
 * \overload
 */
template <typename Topology, concepts::VertexAndEdgeListGraph Graph,
          concepts::LvalueVertexPropertyMap<Graph> PositionMap, concepts::ReadableEdgePropertyMap<Graph> WeightMap,
          typename T, bool EdgeOrSideLength, typename Done, concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
bool kamada_kawai_spring_layout(const Graph& g, PositionMap position, WeightMap weight, const Topology& topology,
                                kamada_kawai_detail::edge_or_side<EdgeOrSideLength, T> edge_or_side_length, Done done,
                                property_traits_value_t<WeightMap> spring_constant, VertexIndexMap index) {
  using weight_type = property_traits_value_t<WeightMap>;
  using point_difference_type = typename Topology::point_difference_type;

  auto n = num_vertices(g);
  using weight_vec = std::vector<weight_type>;

  std::vector<weight_vec> distance(n, weight_vec(n));
  std::vector<weight_vec> spring_strength(n, weight_vec(n));

  return kamada_kawai_spring_layout(g, position, weight, topology, edge_or_side_length, done, spring_constant, index,
                                    distance.begin(), spring_strength.begin(),
                                    vector_property_map(n, index, point_difference_type()).ref());
}

/**
 * \overload
 */
template <typename Topology, concepts::VertexAndEdgeListGraph Graph,
          concepts::LvalueVertexPropertyMap<Graph> PositionMap, concepts::ReadableEdgePropertyMap<Graph> WeightMap,
          typename T, bool EdgeOrSideLength, typename Done>
bool kamada_kawai_spring_layout(const Graph& g, PositionMap position, WeightMap weight, const Topology& topology,
                                kamada_kawai_detail::edge_or_side<EdgeOrSideLength, T> edge_or_side_length, Done done,
                                property_traits_value_t<WeightMap> spring_constant) {
  return kamada_kawai_spring_layout(g, position, weight, topology, edge_or_side_length, done, spring_constant,
                                    get(vertex_index, g));
}

/**
 * \overload
 */
template <typename Topology, concepts::VertexAndEdgeListGraph Graph,
          concepts::LvalueVertexPropertyMap<Graph> PositionMap, concepts::ReadableEdgePropertyMap<Graph> WeightMap,
          typename T, bool EdgeOrSideLength, typename Done>
bool kamada_kawai_spring_layout(const Graph& g, PositionMap position, WeightMap weight, const Topology& topology,
                                kamada_kawai_detail::edge_or_side<EdgeOrSideLength, T> edge_or_side_length, Done done) {
  using weight_type = property_traits_value_t<WeightMap>;
  return kamada_kawai_spring_layout(g, position, weight, topology, edge_or_side_length, done, weight_type(1));
}

/**
 * \overload
 */
template <typename Topology, concepts::VertexAndEdgeListGraph Graph,
          concepts::LvalueVertexPropertyMap<Graph> PositionMap, concepts::ReadableEdgePropertyMap<Graph> WeightMap,
          typename T, bool EdgeOrSideLength>
bool kamada_kawai_spring_layout(const Graph& g, PositionMap position, WeightMap weight, const Topology& topology,
                                kamada_kawai_detail::edge_or_side<EdgeOrSideLength, T> edge_or_side_length) {
  using weight_type = property_traits_value_t<WeightMap>;
  return kamada_kawai_spring_layout(g, position, weight, topology, edge_or_side_length, layout_tolerance<weight_type>(),
                                    weight_type(1.0), get(vertex_index, g));
}
}  // namespace bagl

#endif  // BAGL_BAGL_KAMADA_KAWAI_SPRING_LAYOUT_H_
