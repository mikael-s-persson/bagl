// Copyright 2010 The Trustees of Indiana University.
//  Authors: Jeremiah Willcock
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_LOOP_ERASED_RANDOM_WALK_H_
#define BAGL_BAGL_LOOP_ERASED_RANDOM_WALK_H_

#include <cassert>
#include <iterator>
#include <vector>

#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/random.h"

namespace bagl {

struct loop_erased_random_walk_stuck : std::exception {
  [[nodiscard]] const char* what() const noexcept override {
    return "Loop-erased random walk found a vertex with no out-edges";
  }
};

// Do a loop-erased random walk from vertex s to any vertex colored black (or
// actually any color other than white or gray) in the color map.  The color
// white is for vertices that are not part of the path, while gray is for
// those that are on the path (for cycle detection).  The vector path is used
// for temporary storage and as the result of the algorithm; while all
// elements of the path except the last have their colors set to gray upon
// return.  Vertex s must start off colored white.
//
// Useful references:
// http://everything2.com/title/loop-erased+random+walk
// Wikipedia page on "Loop-Erased Random Walk"
//
// For examples of next_edge, see below:
//  `unweighted_random_out_edge_gen`
//  `weighted_random_out_edge_gen`
template <concepts::IncidenceGraph G, concepts::ReadWriteVertexPropertyMap<G> ColorMap, typename NextEdge>
void loop_erased_random_walk(const G& g, graph_vertex_descriptor_t<G> s, NextEdge next_edge, ColorMap color,
                             std::vector<graph_vertex_descriptor_t<G> >& path) {
  using vertex_descriptor = graph_vertex_descriptor_t<G>;
  using edge_descriptor = graph_edge_descriptor_t<G>;
  using color_t = property_traits_value_t<ColorMap>;
  using color_gen = color_traits<color_t>;

  assert(get(color, s) == color_gen::white());
  path.clear();
  path.push_back(s);
  put(color, s, color_gen::gray());
  while (true) {
    edge_descriptor e = next_edge(s, g);
    vertex_descriptor t = target(e, g);
    color_t t_color = get(color, t);
    if (t_color == color_gen::white()) {
      path.push_back(t);
      put(color, t, color_gen::gray());
      s = t;
    } else if (t_color == color_gen::gray()) {
      // Found a loop; delete from path from the first occurrence of t to
      // the end, coloring vertices white.
      auto it = std::find(path.begin(), path.end(), t);
      assert(it != path.end());
      ++it;
      for (auto j = it; j != path.end(); ++j) {
        put(color, *j, color_gen::white());
      }
      path.erase(it, path.end());
      s = t;
    } else {
      // Done
      path.push_back(t);
      break;
    }
  }
}

template <concepts::IncidenceGraph G, std::uniform_random_bit_generator Gen>
class unweighted_random_out_edge_gen {
  Gen& gen_;

 public:
  explicit unweighted_random_out_edge_gen(Gen& gen) : gen_(gen) {}

  auto operator()(graph_vertex_descriptor_t<G> src, const G& g) const {
    if (out_degree(src, g) == 0) {
      throw loop_erased_random_walk_stuck();
    }
    return random_out_edge(g, src, gen_);
  }
};

template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          std::uniform_random_bit_generator Gen>
class weighted_random_out_edge_gen {
  WeightMap weight_;
  Gen& gen_;

 public:
  weighted_random_out_edge_gen(const WeightMap& weight, Gen& gen) : weight_(weight), gen_(gen) {}

  auto operator()(graph_vertex_descriptor_t<G> src, const G& g) const {
    if (out_degree(src, g) == 0) {
      throw loop_erased_random_walk_stuck();
    }
    return weighted_random_out_edge(g, src, weight_, gen_);
  }
};

}  // namespace bagl

#endif  // BAGL_BAGL_LOOP_ERASED_RANDOM_WALK_H_
