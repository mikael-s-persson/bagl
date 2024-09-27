// Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_TIERNAN_ALL_CYCLES_H_
#define BAGL_BAGL_TIERNAN_ALL_CYCLES_H_

#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"

namespace bagl {
namespace concepts {
template <typename V, typename Path, typename G>
concept CycleVisitor = std::copy_constructible<V> && requires(V& vis, const Path& p, const G& g) {
  vis.cycle(p, g);
};
}  // namespace concepts

// The implementation of this algorithm is a reproduction of the Teirnan
// approach for directed graphs: bibtex follows
//
//     @article{362819,
//         author = {James C. Tiernan},
//         title = {An efficient search algorithm to find the elementary
//         circuits of a graph}, journal = {Commun. ACM}, volume = {13}, number
//         = {12}, year = {1970}, issn = {0001-0782}, pages = {722--726}, doi =
//         {http://doi.acm.org/10.1145/362814.362819},
//             publisher = {ACM Press},
//             address = {New York, NY, USA},
//         }
//
// It should be pointed out that the author does not provide a complete analysis
// for either time or space. This is in part, due to the fact that it's a fairly
// input sensitive problem related to the density and construction of the graph,
// not just its size.
//
// I've also taken some liberties with the interpretation of the algorithm -
// I've basically modernized it to use real data structures (no more arrays and
// matrices). Oh... and there's explicit control structures - not just gotos.
//
// The problem is definitely NP-complete, an unbounded implementation of this
// will probably run for quite a while on a large graph. The conclusions
// of this paper also reference a Paton algorithm for undirected graphs as being
// much more efficient (apparently based on spanning trees). Although not
// implemented, it can be found here:
//
//     @article{363232,
//         author = {Keith Paton},
//         title = {An algorithm for finding a fundamental set of cycles of a
//         graph}, journal = {Commun. ACM}, volume = {12}, number = {9}, year =
//         {1969}, issn = {0001-0782}, pages = {514--518}, doi =
//         {http://doi.acm.org/10.1145/363219.363232},
//             publisher = {ACM Press},
//             address = {New York, NY, USA},
//         }

// The default cycle visitor provides an empty visit function for cycle visitors.
struct cycle_visitor {
  template <typename Path, typename Graph>
  void cycle(const Path& p, const Graph& g) const {}
};

// The min_max_cycle_visitor simultaneously records the minimum and maximum cycles in a graph.
struct min_max_cycle_visitor {
  min_max_cycle_visitor(std::size_t& min_, std::size_t& max_) : minimum(min_), maximum(max_) {}

  template <typename Path, typename Graph>
  void cycle(const Path& p, const Graph& /*g*/) {
    std::size_t len = p.size();
    minimum = std::min(minimum, len);
    maximum = std::max(maximum, len);
  }
  std::size_t& minimum;
  std::size_t& maximum;
};

inline min_max_cycle_visitor find_min_max_cycle(std::size_t& min_, std::size_t& max_) { return {min_, max_}; }

namespace cycles_detail {
template <typename G, typename Path>
inline bool is_vertex_in_path(const G& /*unused*/, graph_vertex_descriptor_t<G> v, const Path& p) {
  return (std::find(p.begin(), p.end(), v) != p.end());
}

template <concepts::VertexIndexGraph G, typename ClosedMatrix>
inline bool is_path_closed(const G& g, graph_vertex_descriptor_t<G> u, graph_vertex_descriptor_t<G> v,
                           const ClosedMatrix& closed) {
  // the path from u to v is closed if v can be found in the list
  // of closed vertices associated with u.
  auto r = closed[get(vertex_index, g, u)];
  return static_cast<bool>(std::find(r.begin(), r.end(), v) != r.end());
}

template <concepts::IncidenceGraph G, typename Path, typename ClosedMatrix>
requires concepts::VertexIndexGraph<G>
bool can_extend_path(const G& g, graph_edge_descriptor_t<G> e, const Path& p, const ClosedMatrix& m) {
  // get the vertices in question
  auto u = source(e, g);
  auto v = target(e, g);

  // conditions for allowing a traversal along this edge are:
  // 1. the index of v must be greater than that at which the
  //    path is rooted (p.front()).
  // 2. the vertex v cannot already be in the path
  // 3. the vertex v cannot be closed to the vertex u

  bool indices = get(vertex_index, g, p.front()) < get(vertex_index, g, v);
  bool path = !is_vertex_in_path(g, v, p);
  bool closed = !is_path_closed(g, u, v, m);
  return indices && path && closed;
}

template <concepts::IncidenceGraph G, typename Path>
bool can_wrap_path(const G& g, const Path& p) {
  // iterate over the out-edges of the back, looking for the
  // front of the path. also, we can't travel along the same
  // edge that we did on the way here, but we don't quite have the
  // stringent requirements that we do in can_extend_path().
  auto u = p.back();
  auto v = p.front();
  for (auto e : out_edges(u, g)) {
    if (target(e, g) == v) {
      return true;
    }
  }
  return false;
}

template <concepts::IncidenceGraph G, typename Path, typename ClosedMatrix>
graph_vertex_descriptor_t<G> extend_path(const G& g, Path& p, ClosedMatrix& closed) {
  // get the current vertex
  auto u = p.back();
  auto ret = graph_traits<G>::null_vertex();

  // AdjacencyIterator i, end;
  for (auto e : out_edges(u, g)) {
    auto v = target(e, g);

    // if we can actually extend along this edge,
    // then that's what we want to do
    if (can_extend_path(g, e, p, closed)) {
      p.push_back(v);  // add the vertex to the path
      ret = v;
      break;
    }
  }
  return ret;
}

template <concepts::Graph G, typename Path, typename ClosedMatrix>
bool exhaust_paths(const G& g, Path& p, ClosedMatrix& closed) {
  // if there's more than one vertex in the path, this closes
  // of some possible routes and returns true. otherwise, if there's
  // only one vertex left, the vertex has been used up
  if (p.size() > 1) {
    // get the last and second to last vertices, popping the last
    // vertex off the path
    auto last = p.back();
    p.pop_back();
    auto prev = p.back();

    // reset the closure for the last vertex of the path and
    // indicate that the last vertex in p is now closed to
    // the next-to-last vertex in p
    closed[get(vertex_index, g, last)].clear();
    closed[get(vertex_index, g, prev)].push_back(last);
    return true;
  }
  return false;
}

template <typename G>
using default_path_t = std::vector<graph_vertex_descriptor_t<G>>;

template <concepts::VertexListGraph G, concepts::CycleVisitor<default_path_t<G>, G> V>
void all_cycles_from_vertex(const G& g, graph_vertex_descriptor_t<G> v, V vis, std::size_t minlen, std::size_t maxlen) {
  using Vertex = graph_vertex_descriptor_t<G>;
  using Path = default_path_t<G>;
  using ClosedMatrix = std::vector<Path>;

  Path p;
  ClosedMatrix closed(num_vertices(g), Path{});
  const Vertex null = graph_traits<G>::null_vertex();

  // each path investigation starts at the ith vertex
  p.push_back(v);

  while (1) {
    // extend the path until we've reached the end or the
    // maxlen-sized cycle
    Vertex j = null;
    while (((j = extend_path(g, p, closed)) != null) && (p.size() < maxlen)) {
    }

    // if we're done extending the path and there's an edge
    // connecting the back to the front, then we should have
    // a cycle.
    if (can_wrap_path(g, p) && p.size() >= minlen) {
      vis.cycle(p, g);
    }

    if (!exhaust_paths(g, p, closed)) {
      break;
    }
  }
}

// Select the minimum allowable length of a cycle based on the directedness
// of the graph - 2 for directed, 3 for undirected.
template <typename G>
constexpr std::size_t min_cycles = (is_undirected_graph_v<G> ? 3 : 2);

}  // namespace cycles_detail

template <concepts::VertexListGraph G, concepts::CycleVisitor<cycles_detail::default_path_t<G>, G> V>
void tiernan_all_cycles(const G& g, V vis, std::size_t minlen, std::size_t maxlen) {
  for (auto u : vertices(g)) {
    cycles_detail::all_cycles_from_vertex(g, u, vis, minlen, maxlen);
  }
}

template <concepts::VertexListGraph G, concepts::CycleVisitor<cycles_detail::default_path_t<G>, G> V>
void tiernan_all_cycles(const G& g, V vis, std::size_t maxlen) {
  tiernan_all_cycles(g, vis, cycles_detail::min_cycles<G>, maxlen);
}

template <concepts::VertexListGraph G, concepts::CycleVisitor<cycles_detail::default_path_t<G>, G> V>
void tiernan_all_cycles(const G& g, V vis) {
  tiernan_all_cycles(g, vis, cycles_detail::min_cycles<G>, std::numeric_limits<std::size_t>::max());
}

template <concepts::VertexListGraph G>
std::pair<std::size_t, std::size_t> tiernan_girth_and_circumference(const G& g) {
  std::size_t accum_min = std::numeric_limits<std::size_t>::max();
  std::size_t accum_max = 0;
  tiernan_all_cycles(g, find_min_max_cycle(accum_min, accum_max));

  // if this is the case, the graph is acyclic...
  if (accum_max == 0) {
    accum_max = accum_min;
  }

  return {accum_min, accum_max};
}

template <concepts::VertexListGraph G>
std::size_t tiernan_girth(const G& g) {
  return tiernan_girth_and_circumference(g).first;
}

template <concepts::VertexListGraph G>
std::size_t tiernan_circumference(const G& g) {
  return tiernan_girth_and_circumference(g).second;
}

}  // namespace bagl

#endif  // BAGL_BAGL_TIERNAN_ALL_CYCLES_H_
