// Copyright (c) Aaron Windsor 2007
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PLANAR_FACE_TRAVERSAL_H_
#define BAGL_BAGL_PLANAR_FACE_TRAVERSAL_H_

#include <iterator>
#include <ranges>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/vector_property_map.h"

namespace bagl {

namespace concepts {

template <typename Visitor, typename Graph>
concept PlanarFaceTraversalVisitor = requires(Visitor& vis, Graph& g, graph_vertex_descriptor_t<Graph> u,
                                              graph_edge_descriptor_t<Graph> e) {
  vis.begin_traversal();
  vis.begin_face();
  vis.next_edge(e, g);
  vis.next_vertex(u, g);
  vis.end_face();
  vis.end_traversal();
};

}  // namespace concepts

struct planar_face_traversal_visitor {
  void begin_traversal() {}

  void begin_face() {}

  template <typename Edge, typename Graph>
  void next_edge(Edge /*unused*/, const Graph& /*unused*/) {}

  template <typename Vertex, typename Graph>
  void next_vertex(Vertex /*unused*/, const Graph& /*unused*/) {}

  void end_face() {}

  void end_traversal() {}
};

template <concepts::VertexAndEdgeListGraph G, concepts::LvalueVertexPropertyMap<G> PlanarEmbedding,
          concepts::PlanarFaceTraversalVisitor<G> Visitor, concepts::ReadableEdgeIndexMap<G> EdgeIndexMap>
void planar_face_traversal(const G& g, PlanarEmbedding embedding, Visitor& vis, EdgeIndexMap em) {
  using vertex_t = graph_vertex_descriptor_t<G>;
  using edge_t = graph_edge_descriptor_t<G>;

  auto visited = vector_property_map(num_edges(g), em, std::unordered_set<vertex_t>{});
  auto next_edge = vector_property_map(num_edges(g), em, std::unordered_map<vertex_t, edge_t>{});

  vis.begin_traversal();

  // Initialize the next_edge property map. This map is initialized from the
  // PlanarEmbedding so that get(next_edge, e)[v] is the edge that comes
  // after e in the clockwise embedding around vertex v.

  for (vertex_t v : vertices(g)) {
    auto p_rg = std::ranges::ref_view(embedding[v]);
    for (auto pi = p_rg.begin(); pi != p_rg.end(); ++pi) {
      edge_t e(*pi);
      next_edge[e][v] = std::next(pi) == p_rg.end() ? *p_rg.begin() : *std::next(pi);
    }
  }

  // Take a copy of the edges in the graph here, since we want to accomodate
  // face traversals that add edges to the graph (for triangulation, in
  // particular) and don't want to use invalidated edge iterators.
  // Also, while iterating over all edges in the graph, we single out
  // any self-loops, which need some special treatment in the face traversal.

  std::vector<edge_t> self_loops;
  std::vector<edge_t> edges_cache;
  edges_cache.reserve(num_edges(g));
  for (edge_t e : edges(g)) {
    edges_cache.push_back(e);
    if (source(e, g) == target(e, g)) {
      self_loops.push_back(e);
    }
  }

  // Iterate over all edges in the graph
  for (edge_t e : edges_cache) {
    // Iterate over both vertices in the current edge
    for (vertex_t v : {source(e, g), target(e, g)}) {
      auto* e_visited = &visited[e];
      auto e_visited_found = e_visited->find(v);

      if (e_visited_found == e_visited->end()) {
        vis.begin_face();
      }

      while (e_visited->find(v) == e_visited->end()) {
        vis.next_vertex(v, g);
        vis.next_edge(e, g);
        e_visited->insert(v);
        v = source(e, g) == v ? target(e, g) : source(e, g);
        e = next_edge[e][v];
        e_visited = &visited[e];
      }

      if (e_visited_found == e_visited->end()) {
        vis.end_face();
      }
    }
  }

  // Iterate over all self-loops, visiting them once separately
  // (they've already been visited once, this visitation is for
  // the "inside" of the self-loop)
  for (edge_t e : self_loops) {
    vis.begin_face();
    vis.next_edge(e, g);
    vis.next_vertex(source(e, g), g);
    vis.end_face();
  }

  vis.end_traversal();
}

template <concepts::VertexAndEdgeListGraph G, concepts::LvalueVertexPropertyMap<G> PlanarEmbedding,
          concepts::PlanarFaceTraversalVisitor<G> Visitor>
void planar_face_traversal(const G& g, PlanarEmbedding embedding, Visitor& vis) {
  planar_face_traversal(g, embedding, vis, get(edge_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_PLANAR_FACE_TRAVERSAL_H_
