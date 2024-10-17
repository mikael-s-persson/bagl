// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Copyright (C) Vladimir Prus 2003
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_RANDOM_H_
#define BAGL_BAGL_RANDOM_H_

#include <cassert>
#include <iostream>
#include <iterator>
#include <random>
#include <type_traits>

#include "bagl/adjacency_list.h"
#include "bagl/copy.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property.h"

namespace bagl {

// grab a random vertex from the graph's vertex set
template <concepts::VertexListGraph Graph, std::uniform_random_bit_generator RandomNumGen>
graph_vertex_descriptor_t<Graph> random_vertex(const Graph& g, RandomNumGen& gen) {
  std::size_t num_v = num_vertices(g);
  if (num_v == 0) {
    return graph_traits<Graph>::null_vertex();
  }
  if (num_v == 1) {
    return *vertices(g).begin();
  }
  std::uniform_int_distribution<std::size_t> distrib(0, num_v - 1);
  std::size_t n = distrib(gen);
  return *(std::next(vertices(g).begin(), n));
}

template <concepts::EdgeListGraph Graph, std::uniform_random_bit_generator RandomNumGen>
graph_edge_descriptor_t<Graph> random_edge(const Graph& g, RandomNumGen& gen) {
  std::size_t num_e = num_edges(g);
  if (num_e == 0) {
    return graph_edge_descriptor_t<Graph>{};
  }
  if (num_e == 1) {
    return *edges(g).begin();
  }
  std::uniform_int_distribution<std::size_t> distrib(0, num_e - 1);
  std::size_t n = distrib(gen);
  return *(std::next(edges(g).begin(), n));
}

template <concepts::IncidenceGraph Graph, std::uniform_random_bit_generator RandomNumGen>
graph_edge_descriptor_t<Graph> random_out_edge(const Graph& g, graph_vertex_descriptor_t<Graph> src,
                                               RandomNumGen& gen) {
  std::size_t num_e = out_degree(src, g);
  if (num_e == 0) {
    return graph_edge_descriptor_t<Graph>{};
  }
  if (num_e == 1) {
    return *out_edges(src, g).begin();
  }
  std::uniform_int_distribution<std::size_t> distrib(0, num_e - 1);
  std::size_t n = distrib(gen);
  return *(std::next(out_edges(src, g).begin(), n));
}

template <concepts::IncidenceGraph Graph, concepts::ReadableEdgePropertyMap<Graph> WeightMap,
          std::uniform_random_bit_generator RandomNumGen>
graph_edge_descriptor_t<Graph> weighted_random_out_edge(const Graph& g, graph_vertex_descriptor_t<Graph> src,
                                                        WeightMap weight, RandomNumGen& gen) {
  using weight_type = property_traits_value_t<WeightMap>;
  weight_type weight_sum(0);
  for (auto e : out_edges(src, g)) {
    weight_sum += get(weight, e);
  }
  std::uniform_real_distribution<weight_type> ur(0, weight_sum);
  weight_type chosen_weight = ur(gen);
  for (auto e : out_edges(src, g)) {
    weight_type w = get(weight, e);
    if (chosen_weight < w) {
      return e;
    }
    chosen_weight -= w;
  }
  assert(false);  // Should not get here
  return graph_edge_descriptor_t<Graph>{};
}

namespace random_detail {
struct dummy_property_copier {
  template <class V1, class V2>
  void operator()(const V1& /*unused*/, const V2& /*unused*/) const {}
};
}  // namespace random_detail

template <concepts::MutableGraph G, std::uniform_random_bit_generator RandNumGen>
void generate_random_graph(G& g, std::size_t v_count, std::size_t e_count, RandNumGen& gen, bool allow_parallel = true,
                            bool self_edges = false) {
  // When parallel edges are not allowed, we create a new graph which
  // does not allow parallel edges, construct it and copy back.
  // This is not efficient if 'g' already disallow parallel edges,
  // but that's task for later.
  if (!allow_parallel) {
    using select = std::conditional_t<is_directed_graph_v<G>, directed_s, undirected_s>;
    adjacency_list<unordered_set_s, vec_s, select> g2;
    generate_random_graph(g2, v_count, e_count, gen, true, self_edges);
    copy_graph(g2, g, random_detail::dummy_property_copier(), random_detail::dummy_property_copier());
  } else {
    for (std::size_t i = 0; i < v_count; ++i) {
      add_vertex(g);
    }

    std::size_t not_inserted_counter = 0;  // Number of edge insertion failures
    std::size_t num_vertices_squared = num_vertices(g) * num_vertices(g);
    for (std::size_t j = 0; j < e_count; /* Increment in body */) {
      auto a = random_vertex(g, gen);
      auto b = a;
      do {
        b = random_vertex(g, gen);
      } while (!self_edges && a == b);
      auto [e, inserted] = add_edge(a, b, g);
      if (inserted) {
        ++j;
      } else {
        ++not_inserted_counter;
      }
      if (not_inserted_counter >= num_vertices_squared) {
        return; /* Rather than looping forever on complete graph */
      }
    }
  }
}

template <concepts::MutableGraph G, std::uniform_random_bit_generator RandNumGen,
          std::output_iterator<graph_vertex_descriptor_t<G>> VertexOutputIterator,
          std::output_iterator<graph_edge_descriptor_t<G>> EdgeOutputIterator>
void generate_random_graph(G& g, std::size_t v_count, std::size_t e_count, RandNumGen& gen,
                           VertexOutputIterator vertex_out, EdgeOutputIterator edge_out, bool self_edges = false) {
  for (std::size_t i = 0; i < v_count; ++i) {
    *vertex_out++ = add_vertex(g);
  }

  std::size_t not_inserted_counter = 0; /* Number of edge insertion failures */
  std::size_t num_vertices_squared = num_vertices(g) * num_vertices(g);
  for (std::size_t j = 0; j < e_count; /* Increment in body */) {
    auto a = random_vertex(g, gen);
    auto b = a;
    do {
      b = random_vertex(g, gen);
    } while (!self_edges && a == b);
    auto [e, inserted] = add_edge(a, b, g);
    if (inserted) {
      *edge_out++ = std::make_pair(source(e, g), target(e, g));
      ++j;
    } else {
      ++not_inserted_counter;
    }
    if (not_inserted_counter >= num_vertices_squared) {
      return; /* Rather than looping forever on complete graph */
    }
  }
}

namespace random_detail {

template <typename Property, concepts::VertexListGraph G, std::invocable<> RandomGenerator>
void randomize_property(G& g, RandomGenerator& rg, Property /*unused*/, vertex_property_tag /*unused*/) {
  auto pm = get(Property(), g);
  for (auto v : vertices(g)) {
    pm[v] = rg();
  }
}

template <typename Property, concepts::EdgeListGraph G, std::invocable<> RandomGenerator>
void randomize_property(G& g, RandomGenerator& rg, Property /*unused*/, edge_property_tag /*unused*/) {
  auto pm = get(Property(), g);
  for (auto e : edges(g)) {
    pm[e] = rg();
  }
}
}  // namespace random_detail

template <typename Property, typename G, std::invocable<> RandomGenerator>
void randomize_property(G& g, RandomGenerator& rg) {
  random_detail::randomize_property(g, rg, Property(), property_kind_t<Property>());
}

}  // namespace bagl

#endif  // BAGL_BAGL_RANDOM_H_
