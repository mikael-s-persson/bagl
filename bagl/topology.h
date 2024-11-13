// Copyright 2009 The Trustees of Indiana University.
// Authors: Jeremiah Willcock, Douglas Gregor, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_TOPOLOGY_H_
#define BAGL_BAGL_TOPOLOGY_H_

#include <array>
#include <cmath>
#include <memory>
#include <numbers>
#include <random>

// Classes and concepts to represent points in a space, with distance and move
// operations (used for Gurson-Atun layout), plus other things like bounding
// boxes used for other layout algorithms.

namespace bagl {

/***********************************************************
 * Topologies                                              *
 ***********************************************************/
template <std::size_t Dims>
class convex_topology {
 public:
  struct point {
    static constexpr std::size_t dimensions = Dims;
    constexpr point() = default;
    constexpr double& operator[](std::size_t i) { return values_[i]; }
    constexpr double operator[](std::size_t i) const { return values_[i]; }

   private:
    std::array<double, Dims> values_{};
  };

  struct point_difference {
    static constexpr std::size_t dimensions = Dims;
    constexpr point_difference() = default;
    constexpr double& operator[](std::size_t i) { return values_[i]; }
    constexpr double operator[](std::size_t i) const { return values_[i]; }

    friend point_difference operator+(const point_difference& a, const point_difference& b) {
      point_difference result;
      for (std::size_t i = 0; i < Dims; ++i) {
        result[i] = a[i] + b[i];
      }
      return result;
    }

    friend point_difference& operator+=(point_difference& a, const point_difference& b) {
      for (std::size_t i = 0; i < Dims; ++i) {
        a[i] += b[i];
      }
      return a;
    }

    friend point_difference operator-(const point_difference& a) {
      point_difference result;
      for (std::size_t i = 0; i < Dims; ++i) {
        result[i] = -a[i];
      }
      return result;
    }

    friend point_difference operator-(const point_difference& a, const point_difference& b) {
      point_difference result;
      for (std::size_t i = 0; i < Dims; ++i) {
        result[i] = a[i] - b[i];
      }
      return result;
    }

    friend point_difference& operator-=(point_difference& a, const point_difference& b) {
      for (std::size_t i = 0; i < Dims; ++i) {
        a[i] -= b[i];
      }
      return a;
    }

    friend point_difference operator*(const point_difference& a, const point_difference& b) {
      point_difference result;
      for (std::size_t i = 0; i < Dims; ++i) {
        result[i] = a[i] * b[i];
      }
      return result;
    }

    friend point_difference operator*(const point_difference& a, double b) {
      point_difference result;
      for (std::size_t i = 0; i < Dims; ++i) {
        result[i] = a[i] * b;
      }
      return result;
    }

    friend point_difference operator*(double a, const point_difference& b) {
      point_difference result;
      for (std::size_t i = 0; i < Dims; ++i) {
        result[i] = a * b[i];
      }
      return result;
    }

    friend point_difference operator/(const point_difference& a, const point_difference& b) {
      point_difference result;
      for (std::size_t i = 0; i < Dims; ++i) {
        result[i] = (b[i] == 0.) ? 0. : a[i] / b[i];
      }
      return result;
    }

    friend double dot(const point_difference& a, const point_difference& b) {
      double result = 0;
      for (std::size_t i = 0; i < Dims; ++i) {
        result += a[i] * b[i];
      }
      return result;
    }

   private:
    std::array<double, Dims> values_{};
  };

  using point_type = point;
  using point_difference_type = point_difference;

  [[nodiscard]] double distance(point a, point b) const {
    double dist = 0.;
    for (std::size_t i = 0; i < Dims; ++i) {
      double diff = b[i] - a[i];
      dist = std::hypot(dist, diff);
    }
    // Exact properties of the distance are not important, as long as
    // < on what this returns matches real distances; l_2 is used because
    // Fruchterman-Reingold also uses this code and it relies on l_2.
    return dist;
  }

  [[nodiscard]] point move_position_toward(point a, double fraction, point b) const {
    point result;
    for (std::size_t i = 0; i < Dims; ++i) {
      result[i] = a[i] + (b[i] - a[i]) * fraction;
    }
    return result;
  }

  [[nodiscard]] point_difference difference(point a, point b) const {
    point_difference result;
    for (std::size_t i = 0; i < Dims; ++i) {
      result[i] = a[i] - b[i];
    }
    return result;
  }

  [[nodiscard]] point adjust(point a, point_difference delta) const {
    point result;
    for (std::size_t i = 0; i < Dims; ++i) {
      result[i] = a[i] + delta[i];
    }
    return result;
  }

  [[nodiscard]] point pointwise_min(point a, point b) const {
    point result;
    for (std::size_t i = 0; i < Dims; ++i) {
      result[i] = std::min(a[i], b[i]);
    }
    return result;
  }

  [[nodiscard]] point pointwise_max(point a, point b) const {
    point result;
    for (std::size_t i = 0; i < Dims; ++i) {
      result[i] = std::max(a[i], b[i]);
    }
    return result;
  }

  [[nodiscard]] double norm(point_difference delta) const {
    double n = 0.;
    for (std::size_t i = 0; i < Dims; ++i) {
      n = std::hypot(n, delta[i]);
    }
    return n;
  }

  [[nodiscard]] double volume(point_difference delta) const {
    double n = 1.;
    for (std::size_t i = 0; i < Dims; ++i) {
      n *= delta[i];
    }
    return n;
  }
};

template <std::size_t Dims, typename RandomNumberGenerator = std::mt19937>
class hypercube_topology : public convex_topology<Dims> {
  using rand_t = std::uniform_real_distribution<double>;

 public:
  using point_type = typename convex_topology<Dims>::point_type;
  using point_difference_type = typename convex_topology<Dims>::point_difference_type;

  explicit hypercube_topology(double scaling = 1.0)
      : gen_ptr_(std::make_shared<RandomNumberGenerator>()), rand_(std::make_shared<rand_t>()), scaling_(scaling) {}

  explicit hypercube_topology(RandomNumberGenerator& gen, double scaling = 1.0)
      : gen_ptr_(std::make_shared<RandomNumberGenerator>(gen)), rand_(std::make_shared<rand_t>()), scaling_(scaling) {}

  [[nodiscard]] point_type random_point() const {
    point_type p;
    for (std::size_t i = 0; i < Dims; ++i) {
      p[i] = (*rand_)(*gen_ptr_) * scaling_;
    }
    return p;
  }

  [[nodiscard]] point_type bound(point_type a) const {
    point_type p;
    for (std::size_t i = 0; i < Dims; ++i) {
      p[i] = std::min(scaling_, std::max(-scaling_, a[i]));
    }
    return p;
  }

  [[nodiscard]] double distance_from_boundary(point_type a) const {
    using std::abs;
    static_assert(Dims >= 1);
    double dist = abs(scaling_ - a[0]);
    for (std::size_t i = 1; i < Dims; ++i) {
      dist = std::min(dist, abs(scaling_ - a[i]));
    }
    return dist;
  }

  [[nodiscard]] point_type center() const {
    point_type result;
    for (std::size_t i = 0; i < Dims; ++i) {
      result[i] = scaling_ * .5;
    }
    return result;
  }

  [[nodiscard]] point_type origin() const {
    point_type result;
    for (std::size_t i = 0; i < Dims; ++i) {
      result[i] = 0;
    }
    return result;
  }

  [[nodiscard]] point_difference_type extent() const {
    point_difference_type result;
    for (std::size_t i = 0; i < Dims; ++i) {
      result[i] = scaling_;
    }
    return result;
  }

 private:
  std::shared_ptr<RandomNumberGenerator> gen_ptr_;
  std::shared_ptr<rand_t> rand_;
  double scaling_;
};

template <typename RandomNumberGenerator = std::mt19937>
class square_topology : public hypercube_topology<2, RandomNumberGenerator> {
  using inherited = hypercube_topology<2, RandomNumberGenerator>;

 public:
  explicit square_topology(double scaling = 1.0) : inherited(scaling) {}

  explicit square_topology(RandomNumberGenerator& gen, double scaling = 1.0) : inherited(gen, scaling) {}
};

template <typename RandomNumberGenerator = std::mt19937>
class rectangle_topology : public convex_topology<2> {
  using rand_t = std::uniform_real_distribution<double>;

 public:
  rectangle_topology(double left, double top, double right, double bottom)
      : gen_ptr_(std::make_shared<RandomNumberGenerator>()),
        rand_(std::make_shared<rand_t>()),
        left_(std::min(left, right)),
        top_(std::min(top, bottom)),
        right_(std::max(left, right)),
        bottom_(std::max(top, bottom)) {}

  rectangle_topology(RandomNumberGenerator& gen, double left, double top, double right, double bottom)
      : gen_ptr_(std::make_shared<RandomNumberGenerator>(gen)),
        rand_(std::make_shared<rand_t>()),
        left_(std::min(left, right)),
        top_(std::min(top, bottom)),
        right_(std::max(left, right)),
        bottom_(std::max(top, bottom)) {}

  using point_type = typename convex_topology<2>::point_type;
  using point_difference_type = typename convex_topology<2>::point_difference_type;

  [[nodiscard]] point_type random_point() const {
    point_type p;
    p[0] = (*rand_)(*gen_ptr_) * (right_ - left_) + left_;
    p[1] = (*rand_)(*gen_ptr_) * (bottom_ - top_) + top_;
    return p;
  }

  [[nodiscard]] point_type bound(point_type a) const {
    point_type p;
    p[0] = std::min(right_, std::max(left_, a[0]));
    p[1] = std::min(bottom_, std::max(top_, a[1]));
    return p;
  }

  [[nodiscard]] double distance_from_boundary(point_type a) const {
    using std::abs;
    double dist = abs(left_ - a[0]);
    dist = std::min(dist, abs(right_ - a[0]));
    dist = std::min(dist, abs(top_ - a[1]));
    dist = std::min(dist, abs(bottom_ - a[1]));
    return dist;
  }

  [[nodiscard]] point_type center() const {
    point_type result;
    result[0] = (left_ + right_) / 2.;
    result[1] = (top_ + bottom_) / 2.;
    return result;
  }

  [[nodiscard]] point_type origin() const {
    point_type result;
    result[0] = left_;
    result[1] = top_;
    return result;
  }

  [[nodiscard]] point_difference_type extent() const {
    point_difference_type result;
    result[0] = right_ - left_;
    result[1] = bottom_ - top_;
    return result;
  }

 private:
  std::shared_ptr<RandomNumberGenerator> gen_ptr_;
  std::shared_ptr<rand_t> rand_;
  double left_;
  double top_;
  double right_;
  double bottom_;
};

template <typename RandomNumberGenerator = std::mt19937>
class cube_topology : public hypercube_topology<3, RandomNumberGenerator> {
  using inherited = hypercube_topology<3, RandomNumberGenerator>;

 public:
  explicit cube_topology(double scaling = 1.0) : inherited(scaling) {}

  explicit cube_topology(RandomNumberGenerator& gen, double scaling = 1.0) : inherited(gen, scaling) {}
};

template <std::size_t Dims, typename RandomNumberGenerator = std::mt19937>
class ball_topology : public convex_topology<Dims> {
  using rand_t = std::uniform_real_distribution<double>;

 public:
  using point_type = typename convex_topology<Dims>::point_type;
  using point_difference_type = typename convex_topology<Dims>::point_difference_type;

  explicit ball_topology(double radius = 1.0)
      : gen_ptr_(std::make_shared<RandomNumberGenerator>()), rand_(std::make_shared<rand_t>()), radius_(radius) {}

  explicit ball_topology(RandomNumberGenerator& gen, double radius = 1.0)
      : gen_ptr_(std::make_shared<RandomNumberGenerator>(gen)), rand_(std::make_shared<rand_t>()), radius_(radius) {}

  [[nodiscard]] point_type random_point() const {
    point_type p;
    double dist_sum = 0.0;
    do {
      dist_sum = 0.0;
      for (std::size_t i = 0; i < Dims; ++i) {
        double x = (*rand_)(*gen_ptr_) * 2 * radius_ - radius_;
        p[i] = x;
        dist_sum += x * x;
      }
    } while (dist_sum > radius_ * radius_);
    return p;
  }

  [[nodiscard]] point_type bound(point_type a) const {
    double r = 0.;
    for (std::size_t i = 0; i < Dims; ++i) {
      r = std::hypot(r, a[i]);
    }
    if (r <= radius_) {
      return a;
    }
    double scaling_factor = radius_ / r;
    point_type p;
    for (std::size_t i = 0; i < Dims; ++i) {
      p[i] = a[i] * scaling_factor;
    }
    return p;
  }

  [[nodiscard]] double distance_from_boundary(point_type a) const {
    double r = 0.;
    for (std::size_t i = 0; i < Dims; ++i) {
      r = std::hypot(r, a[i]);
    }
    return radius_ - r;
  }

  [[nodiscard]] point_type center() const {
    point_type result;
    for (std::size_t i = 0; i < Dims; ++i) {
      result[i] = 0;
    }
    return result;
  }

  [[nodiscard]] point_type origin() const {
    point_type result;
    for (std::size_t i = 0; i < Dims; ++i) {
      result[i] = -radius_;
    }
    return result;
  }

  [[nodiscard]] point_difference_type extent() const {
    point_difference_type result;
    for (std::size_t i = 0; i < Dims; ++i) {
      result[i] = 2. * radius_;
    }
    return result;
  }

 private:
  std::shared_ptr<RandomNumberGenerator> gen_ptr_;
  std::shared_ptr<rand_t> rand_;
  double radius_;
};

template <typename RandomNumberGenerator = std::mt19937>
class circle_topology : public ball_topology<2, RandomNumberGenerator> {
  using inherited = ball_topology<2, RandomNumberGenerator>;

 public:
  explicit circle_topology(double radius = 1.0) : inherited(radius) {}

  explicit circle_topology(RandomNumberGenerator& gen, double radius = 1.0) : inherited(gen, radius) {}
};

template <typename RandomNumberGenerator = std::mt19937>
class sphere_topology : public ball_topology<3, RandomNumberGenerator> {
  using inherited = ball_topology<3, RandomNumberGenerator>;

 public:
  explicit sphere_topology(double radius = 1.0) : inherited(radius) {}

  explicit sphere_topology(RandomNumberGenerator& gen, double radius = 1.0) : inherited(gen, radius) {}
};

template <typename RandomNumberGenerator = std::mt19937>
class heart_topology {
  // Heart is defined as the union of three shapes:
  // Square w/ corners (+-1000, -1000), (0, 0), (0, -2000)
  // Circle centered at (-500, -500) radius 500*sqrt(2)
  // Circle centered at (500, -500) radius 500*sqrt(2)
  // Bounding box (-1000, -2000) - (1000, 500*(sqrt(2) - 1))

  struct point {
    point() {
      values_[0] = 0.0;
      values_[1] = 0.0;
    }
    point(double x, double y) {
      values_[0] = x;
      values_[1] = y;
    }

    double& operator[](std::size_t i) { return values_[i]; }
    double operator[](std::size_t i) const { return values_[i]; }

   private:
    std::array<double, 2> values_{};
  };

  using rand_t = std::uniform_real_distribution<double>;

 public:
  using point_type = point;

  heart_topology() : gen_ptr_(std::make_shared<RandomNumberGenerator>()), rand_(std::make_shared<rand_t>()) {}

  explicit heart_topology(RandomNumberGenerator& gen)
      : gen_ptr_(std::make_shared<RandomNumberGenerator>(gen)), rand_(std::make_shared<rand_t>()) {}

  point random_point() const {
    point result;
    do {
      result[0] = (*rand_)(*gen_ptr_) *
                  (1000 + 1000 * std::numbers::sqrt2_v<double>)-(500 + 500 * std::numbers::sqrt2_v<double>);
      result[1] = (*rand_)(*gen_ptr_) * (2000 + 500 * (std::numbers::sqrt2_v<double> - 1)) - 2000;
    } while (!in_heart(result));
    return result;
  }

  [[nodiscard]] bool in_heart(point p) const {
    using std::abs;

    if (p[1] < abs(p[0]) - 2000) {
      return false;  // Bottom
    }
    if (p[1] <= -1000) {
      return true;  // Diagonal of square
    }
    if (std::hypot(p[0] - -500, p[1] - -500) <= 500. * std::numbers::sqrt2_v<double>) {
      return true;  // Left circle
    }
    if (std::hypot(p[0] - 500, p[1] - -500) <= 500. * std::numbers::sqrt2_v<double>) {
      return true;  // Right circle
    }
    return false;
  }

  [[nodiscard]] bool segment_within_heart(point p1, point p2) const {
    // Assumes that p1 and p2 are within the heart
    if ((p1[0] < 0) == (p2[0] < 0)) {
      return true;  // Same side of symmetry line
    }
    if (p1[0] == p2[0]) {
      return true;  // Vertical
    }
    double slope = (p2[1] - p1[1]) / (p2[0] - p1[0]);
    double intercept = p1[1] - p1[0] * slope;
    return (intercept <= 0);
  }

  // Not going to provide clipping to bounding region or distance from
  // boundary

  double distance(point a, point b) const {
    if (segment_within_heart(a, b)) {
      // Straight line
      return std::hypot(b[0] - a[0], b[1] - a[1]);
    }
    // Straight line bending around (0, 0)
    return std::hypot(a[0], a[1]) + std::hypot(b[0], b[1]);
  }

  point move_position_toward(point a, double fraction, point b) const {
    if (segment_within_heart(a, b)) {
      // Straight line
      return point(a[0] + (b[0] - a[0]) * fraction, a[1] + (b[1] - a[1]) * fraction);
    }
    double distance_to_point_a = std::hypot(a[0], a[1]);
    double distance_to_point_b = std::hypot(b[0], b[1]);
    double location_of_point = distance_to_point_a / (distance_to_point_a + distance_to_point_b);
    if (fraction < location_of_point) {
      return point(a[0] * (1 - fraction / location_of_point), a[1] * (1 - fraction / location_of_point));
    }
    return point(b[0] * ((fraction - location_of_point) / (1 - location_of_point)),
                 b[1] * ((fraction - location_of_point) / (1 - location_of_point)));
  }

 private:
  std::shared_ptr<RandomNumberGenerator> gen_ptr_;
  std::shared_ptr<rand_t> rand_;
};

}  // namespace bagl

#endif  // BAGL_BAGL_TOPOLOGY_H_
