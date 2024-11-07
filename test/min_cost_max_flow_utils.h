// Copyright 2013 University of Warsaw.
// Authors: Piotr Wygocki
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_TEST_MIN_COST_MAX_FLOW_UTILS_H_
#define BAGL_TEST_MIN_COST_MAX_FLOW_UTILS_H_

#include <cstdlib>
#include <iostream>

#include "bagl/adjacency_list.h"
#include "bagl/graph_traits.h"
#include "bagl/vector_property_map.h"
#include "gtest/gtest.h"

namespace bagl {

struct sample_graph {
  using Graph =
      adjacency_list<vec_s, vec_s, directed_s, no_property,
                     property<edge_capacity_t, int,
                              property<edge_residual_capacity_t, int,
                                       property<edge_weight_t, int, property<edge_index_t, std::size_t> > > > >;

  template <typename OGraph, typename OReversed>
  class edge_adder {
   public:
    using OVertex = graph_vertex_descriptor_t<OGraph>;
    edge_adder(OGraph& g, OReversed& rev) : g_(g), rev_(rev) {}
    void add(OVertex v, OVertex w, int weight, int capacity) {
      auto e = add_raw(v, w, weight, capacity);
      auto f = add_raw(w, v, -weight, 0);
      rev_[e] = f;
      rev_[f] = e;
    }

   private:
    auto add_raw(OVertex v, OVertex w, int weight, int capacity) {
      auto [e, b] = add_edge(*std::next(vertices(g_).begin(), v), *std::next(vertices(g_).begin(), w), g_, capacity, 0,
                             weight, num_edges(g_));
      EXPECT_TRUE(b);
      return e;
    }
    OGraph& g_;
    OReversed& rev_;
  };

  template <typename OGraph, typename OReversed>
  static auto get_sample_graph(OGraph& g, OReversed rev) {
    const std::size_t N(6);

    for (std::size_t i = 0; i < N; ++i) {
      add_vertex(g);
    }

    auto s = *vertices(g).begin();
    auto t = *std::next(vertices(g).begin(), 5);

    edge_adder ea(g, rev);

    ea.add(0, 1, 4, 2);
    ea.add(0, 2, 2, 2);

    ea.add(1, 3, 2, 2);
    ea.add(1, 4, 1, 1);
    ea.add(2, 3, 1, 1);
    ea.add(2, 4, 1, 1);

    ea.add(3, 5, 4, 20);
    ea.add(4, 5, 2, 20);

    return std::pair{s, t};
  }

  template <typename OGraph, typename OReversed>
  static auto get_sample_graph_2(OGraph& g, OReversed rev) {
    const std::size_t N(5);

    for (std::size_t i = 0; i < N; ++i) {
      add_vertex(g);
    }

    auto s = *vertices(g).begin();
    auto t = *std::next(vertices(g).begin(), 4);

    edge_adder ea(g, rev);

    ea.add(0, 1, 4, 2);
    ea.add(0, 2, 2, 2);
    ea.add(1, 2, 2, 2);
    ea.add(2, 3, 1, 1);
    ea.add(2, 4, 1, 1);
    ea.add(3, 4, 1, 1);

    ea.add(1, 0, 2, 2);
    ea.add(2, 0, 1, 1);
    ea.add(2, 1, 5, 2);
    ea.add(3, 2, 1, 1);
    ea.add(4, 2, 2, 2);
    ea.add(4, 3, 1, 3);

    return std::pair{s, t};
  }
};
}  // namespace bagl

#endif  // BAGL_TEST_MIN_COST_MAX_FLOW_UTILS_H_
