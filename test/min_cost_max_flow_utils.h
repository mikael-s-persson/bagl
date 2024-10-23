// Copyright 2013 University of Warsaw.
// Authors: Piotr Wygocki
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_TEST_MIN_COST_MAX_FLOW_UTILS_H_
#define BAGL_TEST_MIN_COST_MAX_FLOW_UTILS_H_

#include <cstdlib>
#include <iostream>

#include "bagl/adjacency_list.h"
#include "bagl/vector_property_map.h"
#include "gtest/gtest.h"

namespace bagl {

struct SampleGraph {
  using Graph =
      adjacency_list<vec_s, vec_s, directed_s, no_property,
                     property<edge_capacity_t, int,
                              property<edge_residual_capacity_t, int,
                                       property<edge_weight_t, int, property<edge_index_t, std::size_t> > > > >;
  using EdgeIndex = property_map_t<Graph, edge_index_t>;
  using Edge = graph_edge_descriptor_t<Graph>;
  using Reversed = property_map_ref<vector_property_map<Edge, EdgeIndex> >;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  class EdgeAdder {
   public:
    EdgeAdder(Graph& g, Reversed& rev) : m_g(g), m_rev(rev) {}
    void addEdge(Vertex v, Vertex w, int weight, int capacity) {
      auto e = add(v, w, weight, capacity);
      auto f = add(w, v, -weight, 0);
      m_rev[e] = f;
      m_rev[f] = e;
    }

   private:
    Edge add(Vertex v, Vertex w, int weight, int capacity) {
      auto [e, b] = add_edge(*std::next(vertices(m_g).begin(), v), *std::next(vertices(m_g).begin(), w), m_g, capacity,
                             0, weight, num_edges(m_g));
      EXPECT_TRUE(b);
      return e;
    }
    Graph& m_g;
    Reversed& m_rev;
  };

  static void getSampleGraph(Graph& g, Vertex& s, Vertex& t, Reversed rev) {
    std::size_t N(6);

    for (std::size_t i = 0; i < N; ++i) {
      add_vertex(g);
    }

    s = 0;
    t = 5;

    EdgeAdder ea(g, rev);

    ea.addEdge(0, 1, 4, 2);
    ea.addEdge(0, 2, 2, 2);

    ea.addEdge(1, 3, 2, 2);
    ea.addEdge(1, 4, 1, 1);
    ea.addEdge(2, 3, 1, 1);
    ea.addEdge(2, 4, 1, 1);

    ea.addEdge(3, 5, 4, 20);
    ea.addEdge(4, 5, 2, 20);
  }

  static void getSampleGraph2(Graph& g, Vertex& s, Vertex& t, Reversed rev) {
    const std::size_t N(5);

    for (std::size_t i = 0; i < N; ++i) {
      add_vertex(g);
    }

    s = 0;
    t = 4;

    EdgeAdder ea(g, rev);

    ea.addEdge(0, 1, 4, 2);
    ea.addEdge(0, 2, 2, 2);
    ea.addEdge(1, 2, 2, 2);
    ea.addEdge(2, 3, 1, 1);
    ea.addEdge(2, 4, 1, 1);
    ea.addEdge(3, 4, 1, 1);

    ea.addEdge(1, 0, 2, 2);
    ea.addEdge(2, 0, 1, 1);
    ea.addEdge(2, 1, 5, 2);
    ea.addEdge(3, 2, 1, 1);
    ea.addEdge(4, 2, 2, 2);
    ea.addEdge(4, 3, 1, 3);
  }
};
}  // namespace bagl

#endif  // BAGL_TEST_MIN_COST_MAX_FLOW_UTILS_H_
