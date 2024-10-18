// Copyright (c) 2004 Kristopher Beevers
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/astar_search.h"

#include <cmath>
#include <iostream>
#include <list>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/random.h"
#include "bagl/vector_property_map.h"
#include "bagl/visitors.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

// auxiliary types
struct location {
  float lat;
  float lon;  // lat, long
};
struct cost {
  float v;
  explicit cost(float a_v = float()) : v(a_v) {}
};
std::ostream& operator<<(std::ostream& o, cost f) { return o << f.v; }
cost operator+(cost a, cost b) { return cost(a.v + b.v); }
bool operator==(cost a, cost b) { return a.v == b.v; }
bool operator<(cost a, cost b) { return a.v < b.v; }

class city {
 public:
  city() = default;
  city(std::string n, location l) : name(n), loc(l) {}
  friend std::ostream& operator<<(std::ostream& out, const city& c) {
    out << "\"" << c.name << "\" (" << c.loc.lat << ", " << c.loc.lon << ")";
    return out;
  }

  std::string name;
  location loc{};
};

// euclidean distance heuristic
template <class LocMap, class Vertex>
class distance_heuristic {
 public:
  distance_heuristic(LocMap l, Vertex goal) : location_(l), goal_(goal) {}
  cost operator()(Vertex u) const {
    float dlon = location_[goal_].lon - location_[u].lon;
    float dlat = location_[goal_].lat - location_[u].lat;
    return cost{std::sqrt(dlon * dlon + dlat * dlat)};
  }

 private:
  LocMap location_;
  Vertex goal_;
};

struct found_goal {};  // exception for termination

// visitor that terminates when we find the goal
template <class Vertex>
class astar_goal_visitor {
 public:
  astar_goal_visitor(Vertex goal) : m_goal(goal) {}
  template <class Graph>
  void examine_vertex(Vertex u, const Graph& g) {
    if (u == m_goal) {
      throw found_goal();
    }
  }

 private:
  Vertex m_goal;
};

TEST(AStarSearchTest, Search) {
  // specify some types
  using Graph = adjacency_list<list_s, vec_s, undirected_s, city, property<edge_weight_t, cost> >;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  // specify data
  enum nodes {
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
  std::tuple<int, int, cost> edge_array[] = {
      {Troy, Utica, cost(96)},
      {Troy, LakePlacid, cost(134)},
      {Troy, Plattsburgh, cost(143)},
      {LakePlacid, Plattsburgh, cost(65)},
      {Plattsburgh, Massena, cost(115)},
      {LakePlacid, Massena, cost(133)},
      {Massena, Watertown, cost(117)},
      {Watertown, Utica, cost(116)},
      {Watertown, Syracuse, cost(74)},
      {Utica, Syracuse, cost(56)},
      {Syracuse, Rochester, cost(84)},
      {Rochester, Buffalo, cost(73)},
      {Syracuse, Ithaca, cost(69)},
      {Ithaca, Binghamton, cost(70)},
      {Ithaca, Rochester, cost(116)},
      {Binghamton, Troy, cost(147)},
      {Binghamton, Woodstock, cost(173)},
      {Binghamton, NewYork, cost(183)},
      {Syracuse, Binghamton, cost(74)},
      {Woodstock, Troy, cost(71)},
      {Woodstock, NewYork, cost(124)},
  };

  // create graph
  Graph g(N);
  for (auto [u, v, c] : edge_array) {
    auto [e, inserted] = add_edge(u, v, g, c);
    ASSERT_TRUE(inserted);
  }

  g[Troy] = {"Troy", {42.73, 73.68}};
  g[LakePlacid] = {"Lake Placid", {44.28, 73.99}};
  g[Plattsburgh] = {"Plattsburgh", {44.70, 73.46}};
  g[Massena] = {"Massena", {44.93, 74.89}};
  g[Watertown] = {"Watertown", {43.97, 75.91}};
  g[Utica] = {"Utica", {43.10, 75.23}};
  g[Syracuse] = {"Syracuse", {43.04, 76.14}};
  g[Rochester] = {"Rochester", {43.17, 77.61}};
  g[Buffalo] = {"Buffalo", {42.89, 78.86}};
  g[Ithaca] = {"Ithaca", {42.44, 76.50}};
  g[Binghamton] = {"Binghamton", {42.10, 75.91}};
  g[Woodstock] = {"Woodstock", {42.04, 74.11}};
  g[NewYork] = {"New York", {40.67, 73.94}};

  auto name = get(&city::name, g);
  auto idx = get(vertex_index, g);
  auto weight = get(edge_weight, g);

  auto pm = vector_property_map(num_vertices(g), idx, graph_traits<Graph>::null_vertex());
  auto cm = vector_property_map(num_vertices(g), idx, cost(std::numeric_limits<float>::max()));
  auto dm = vector_property_map(num_vertices(g), idx, cost(std::numeric_limits<float>::max()));
  auto colorm = vector_property_map(num_vertices(g), idx, default_color_type::white_color);

  for (int i = 0; i < 5; ++i) {
    // pick random start/goal
    std::mt19937 gen(42);
    Vertex start = gen() % num_vertices(g);
    Vertex goal = gen() % num_vertices(g);

    auto print_problem_def = [&]() {
      std::stringstream out;
      print_graph(g, get(vertex_bundle, g), weight, out);
      out << "\nStart vertex: " << name[start];
      out << "\nGoal vertex: " << name[goal];
      return out.str();
    };
    // Print problem definition if any checks fail.
    SCOPED_TRACE(print_problem_def());

    // call astar named parameter interface
    EXPECT_THROW(astar_search(g, start, distance_heuristic(get(&city::loc, g), goal),
                              make_astar_visitor(astar_goal_visitor(goal)), pm, cm, dm, weight, idx, colorm,
                              std::less<>(), std::plus<>(), cost(std::numeric_limits<float>::max()), cost(0)),
                 found_goal)
        << "Did not find a path to the goal at all!";
    std::vector<Vertex> shortest_path;
    cost path_weight_sum{0};
    for (Vertex v = goal;; v = pm[v]) {
      shortest_path.push_back(v);
      if (pm[v] == v) {
        // Expect to end up back at start.
        EXPECT_EQ(v, start);
        break;
      }
      auto [e, e_found] = edge(pm[v], v, g);
      ASSERT_TRUE(e_found) << "Edge from " << name[pm[v]] << " to " << name[v] << " does not exist!";
      EXPECT_THAT((dm[pm[v]] + get(weight, e)).v, ::testing::FloatNear(dm[v].v, 1.0e-5F));
      path_weight_sum = path_weight_sum + get(weight, e);
    }
    std::ranges::reverse(shortest_path);
    auto print_solution = [&]() {
      std::stringstream out;
      out << "Shortest path from " << name[start] << " to " << name[goal] << ": ";
      bool first_u = true;
      for (auto spi : shortest_path) {
        if (!first_u) {
          out << " -> ";
        }
        out << name[spi];
        first_u = false;
      }
      out << "\nTotal travel time: " << dm[goal];
      return out.str();
    };
    // Print solution if any checks fail.
    SCOPED_TRACE(print_solution());

    EXPECT_THAT(dm[goal].v, ::testing::FloatNear(path_weight_sum.v, 1.0e-5F));
    for (auto v : vertices(g)) {
      if (pm[v] == v) {
        continue;
      }
      // Recorded predecessor is the best.
      for (auto e : out_edges(v, g)) {
        EXPECT_GE((dm[target(e, g)] + get(weight, e)).v, dm[v].v - 1e-5F);
      }
    }
    cost early_stop_cost = dm[goal];
    // Run astar exhaustively, to be sure.
    astar_search(g, start, distance_heuristic(get(&city::loc, g), goal), default_astar_visitor(), pm, cm, dm, weight,
                 idx, colorm, std::less<>(), std::plus<>(), cost(std::numeric_limits<float>::max()), cost(0));
    EXPECT_THAT(early_stop_cost.v, ::testing::FloatNear(dm[goal].v, 1.0e-5F));
  }
}

}  // namespace
}  // namespace bagl
