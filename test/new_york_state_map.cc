// Copyright 2024 Mikael Persson

#include "test/new_york_state_map.h"

#include "gtest/gtest.h"

namespace bagl::test {

state_map_type make_new_york_state_map() {
  std::tuple<int, int, float> edge_array[] = {
      {nys::Troy, nys::Utica, 96},
      {nys::Troy, nys::LakePlacid, 134},
      {nys::Troy, nys::Plattsburgh, 143},
      {nys::LakePlacid, nys::Plattsburgh, 65},
      {nys::Plattsburgh, nys::Massena, 115},
      {nys::LakePlacid, nys::Massena, 133},
      {nys::Massena, nys::Watertown, 117},
      {nys::Watertown, nys::Utica, 116},
      {nys::Watertown, nys::Syracuse, 74},
      {nys::Utica, nys::Syracuse, 56},
      {nys::Syracuse, nys::Rochester, 84},
      {nys::Rochester, nys::Buffalo, 73},
      {nys::Syracuse, nys::Ithaca, 69},
      {nys::Ithaca, nys::Binghamton, 70},
      {nys::Ithaca, nys::Rochester, 116},
      {nys::Binghamton, nys::Troy, 147},
      {nys::Binghamton, nys::Woodstock, 173},
      {nys::Binghamton, nys::NewYork, 183},
      {nys::Syracuse, nys::Binghamton, 74},
      {nys::Woodstock, nys::Troy, 71},
      {nys::Woodstock, nys::NewYork, 124},
  };

  // create graph
  state_map_type g(nys::N);
  for (auto [u, v, c] : edge_array) {
    auto [e, inserted] = add_edge(u, v, g, travel_time{c});
    EXPECT_TRUE(inserted);
  }

  g[nys::Troy] = {"Troy", {42.73, 73.68}};
  g[nys::LakePlacid] = {"Lake Placid", {44.28, 73.99}};
  g[nys::Plattsburgh] = {"Plattsburgh", {44.70, 73.46}};
  g[nys::Massena] = {"Massena", {44.93, 74.89}};
  g[nys::Watertown] = {"Watertown", {43.97, 75.91}};
  g[nys::Utica] = {"Utica", {43.10, 75.23}};
  g[nys::Syracuse] = {"Syracuse", {43.04, 76.14}};
  g[nys::Rochester] = {"Rochester", {43.17, 77.61}};
  g[nys::Buffalo] = {"Buffalo", {42.89, 78.86}};
  g[nys::Ithaca] = {"Ithaca", {42.44, 76.50}};
  g[nys::Binghamton] = {"Binghamton", {42.10, 75.91}};
  g[nys::Woodstock] = {"Woodstock", {42.04, 74.11}};
  g[nys::NewYork] = {"New York", {40.67, 73.94}};

  return g;
}

}  // namespace bagl::test
