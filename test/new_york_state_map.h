// Copyright 2024 Mikael Persson

#ifndef BAGL_TEST_NEW_YORK_STATE_MAP_H_
#define BAGL_TEST_NEW_YORK_STATE_MAP_H_

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>

#include "bagl/adjacency_list.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"

namespace bagl {
namespace test {

// auxiliary types
struct location {
  float lat = 0.0F;
  float lon = 0.0F;
};
struct travel_time {
  float minutes = 0.0F;
  explicit travel_time(float m = float()) : minutes(m) {}
};
inline std::ostream& operator<<(std::ostream& o, travel_time f) { return o << f.minutes; }
inline travel_time operator+(travel_time a, travel_time b) { return travel_time(a.minutes + b.minutes); }
inline bool operator==(travel_time a, travel_time b) { return a.minutes == b.minutes; }
inline bool operator<(travel_time a, travel_time b) { return a.minutes < b.minutes; }

class city {
 public:
  city() = default;
  city(std::string n, location l) : name(std::move(n)), loc(l) {}

  std::string name;
  location loc{};
};
inline std::ostream& operator<<(std::ostream& out, const city& c) {
  out << "\"" << c.name << "\" (" << c.loc.lat << ", " << c.loc.lon << ")";
  return out;
}

using state_map_type = adjacency_list<list_s, vec_s, undirected_s, city, travel_time>;

// bird-flight travel-time heuristic
class travel_time_heuristic {
 public:
  travel_time_heuristic(const state_map_type& m, graph_vertex_descriptor_t<state_map_type> goal)
      : location_(get(&city::loc, m)), goal_(goal) {}
  float operator()(graph_vertex_descriptor_t<state_map_type> u) const {
    float min_lat_deg = std::min(location_[goal_].lat, location_[u].lat);
    float min_lat_rad = std::numbers::pi_v<float> * min_lat_deg / 180.0F;
    static constexpr float deg_to_miles = 69.1F;
    float tlat = (location_[goal_].lat - location_[u].lat) * deg_to_miles;
    float tlon = (location_[goal_].lon - location_[u].lon) * std::cos(min_lat_rad) * deg_to_miles;
    // One mile a minute.
    return std::sqrt(tlon * tlon + tlat * tlat);
  }

 private:
  decltype(get(&city::loc, std::declval<const state_map_type&>())) location_;
  graph_vertex_descriptor_t<state_map_type> goal_;
};

struct nys {
  enum cities {
    Troy,
    LakePlacid,
    Plattsburgh,
    Massena,
    Watertown,
    Utica,
    Syracuse,
    Rochester,
    Buffalo,
    Ithaca,
    Binghamton,
    Woodstock,
    NewYork,
    N
  };
};

state_map_type make_new_york_state_map();

}  // namespace test

inline auto get(vertex_name_t /*unused*/, const test::state_map_type& m) { return get(&test::city::name, m); }
inline auto get(vertex_name_t /*unused*/, test::state_map_type& m) { return get(&test::city::name, m); }

inline auto get(edge_weight_t /*unused*/, const test::state_map_type& m) { return get(&test::travel_time::minutes, m); }
inline auto get(edge_weight_t /*unused*/, test::state_map_type& m) { return get(&test::travel_time::minutes, m); }

}  // namespace bagl

#endif  // BAGL_TEST_NEW_YORK_STATE_MAP_H_
