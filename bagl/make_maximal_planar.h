// Copyright 2007 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_MAKE_MAXIMAL_PLANAR_H_
#define BAGL_BAGL_MAKE_MAXIMAL_PLANAR_H_

#include <algorithm>
#include <iterator>
#include <tuple>
#include <vector>

#include "bagl/biconnected_components.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/planar_detail/add_edge_visitors.h"
#include "bagl/planar_face_traversal.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"

namespace bagl {

template <typename Graph, typename VertexIndexMap, typename AddEdgeVisitor>
struct triangulation_visitor : public planar_face_traversal_visitor {
  using vertex_t = graph_vertex_descriptor_t<Graph>;
  using edge_t = graph_edge_descriptor_t<Graph>;
  using vertex_vector_t = std::vector<vertex_t>;
  using vertex_to_size_map_t = vector_property_map<std::size_t, VertexIndexMap>;

  triangulation_visitor(Graph& arg_g, VertexIndexMap arg_vm, AddEdgeVisitor arg_add_edge_visitor)
      : g_(arg_g),
        vm_(arg_vm),
        add_edge_visitor_(arg_add_edge_visitor),
        marked_(num_vertices(g_), vm_, 0),
        degree_(num_vertices(g_), vm_, 0) {
    for (auto v : vertices(g_)) {
      put(degree_, v, out_degree(v, g_));
    }
  }

  void next_vertex(vertex_t v, const Graph& /*g*/) {
    // Self-loops will appear as consecutive vertices in the list of
    // vertices on a face. We want to skip these.
    if (!vertices_on_face_.empty() && (vertices_on_face_.back() == v || vertices_on_face_.front() == v)) {
      return;
    }
    vertices_on_face_.push_back(v);
  }

  void end_face() {
    ++timestamp_;

    if (vertices_on_face_.size() <= 3) {
      // At most three vertices on this face - don't need to triangulate
      vertices_on_face_.clear();
      return;
    }

    // Find vertex on face of minimum degree
    std::size_t min_degree = num_vertices(g_);
    auto min_degree_vertex_itr = vertices_on_face_.end();
    for (auto fi = vertices_on_face_.begin(); fi != vertices_on_face_.end(); ++fi) {
      std::size_t deg = get(degree_, *fi);
      if (deg < min_degree) {
        min_degree_vertex_itr = fi;
        min_degree = deg;
      }
    }

    // To simplify some of the manipulations, we'll re-arrange
    // vertices_on_face so that it still contains the same
    // (counter-clockwise) order of the vertices on this face, but now the
    // min_degree_vertex is the first element in vertices_on_face.
    std::rotate(vertices_on_face_.begin(), min_degree_vertex_itr, vertices_on_face_.end());

    // Mark all of the min degree vertex's neighbors
    for (auto a : adjacent_vertices(vertices_on_face_.front(), g_)) {
      put(marked_, a, timestamp_);
    }

    // The iterator manipulations on the next two lines are safe because
    // vertices_on_face.size() > 3 (from the first test in this function)
    auto marked_neighbor = vertices_on_face_.end();
    auto fi_end = std::prev(vertices_on_face_.end());
    for (auto fi = std::next(std::next(vertices_on_face_.begin())); fi != fi_end; ++fi) {
      if (get(marked_, *fi) == timestamp_) {
        marked_neighbor = fi;
        break;
      }
    }

    if (marked_neighbor == vertices_on_face_.end()) {
      add_edge_range(vertices_on_face_[0], std::next(std::next(vertices_on_face_.begin())),
                     std::prev(vertices_on_face_.end()));
    } else {
      add_edge_range(vertices_on_face_[1], std::next(marked_neighbor), vertices_on_face_.end());
      add_edge_range(*std::next(marked_neighbor), std::next(std::next(vertices_on_face_.begin())), marked_neighbor);
    }

    // reset for the next face
    vertices_on_face_.clear();
  }

 private:
  template <typename FaceIter>
  void add_edge_range(vertex_t anchor, FaceIter fi, FaceIter fi_end) {
    for (; fi != fi_end; ++fi) {
      vertex_t v(*fi);
      add_edge_visitor_.visit_vertex_pair(anchor, v, g_);
      put(degree_, anchor, get(degree_, anchor) + 1);
      put(degree_, v, get(degree_, v) + 1);
    }
  }

  Graph& g_;
  VertexIndexMap vm_;
  AddEdgeVisitor add_edge_visitor_;
  std::size_t timestamp_ = 0;
  vertex_vector_t vertices_on_face_;
  vertex_to_size_map_t marked_;
  vertex_to_size_map_t degree_;
};

template <concepts::VertexAndEdgeListGraph G, concepts::LvalueVertexPropertyMap<G> PlanarEmbedding,
          concepts::ReadableVertexIndexMap<G> VertexIndexMap, concepts::ReadableEdgeIndexMap<G> EdgeIndexMap,
          typename AddEdgeVisitor>
void make_maximal_planar(G& g, PlanarEmbedding embedding, VertexIndexMap vm, EdgeIndexMap em, AddEdgeVisitor& vis) {
  planar_face_traversal(g, embedding, triangulation_visitor(g, vm, vis), em);
}

template <concepts::VertexAndEdgeListGraph G, concepts::LvalueVertexPropertyMap<G> PlanarEmbedding,
          concepts::ReadableVertexIndexMap<G> VertexIndexMap, concepts::ReadableEdgeIndexMap<G> EdgeIndexMap>
void make_maximal_planar(G& g, PlanarEmbedding embedding, VertexIndexMap vm, EdgeIndexMap em) {
  default_add_edge_visitor vis;
  make_maximal_planar(g, embedding, vm, em, vis);
}

template <concepts::VertexAndEdgeListGraph G, concepts::LvalueVertexPropertyMap<G> PlanarEmbedding,
          concepts::ReadableVertexIndexMap<G> VertexIndexMap>
void make_maximal_planar(G& g, PlanarEmbedding embedding, VertexIndexMap vm) {
  make_maximal_planar(g, embedding, vm, get(edge_index, g));
}

template <concepts::VertexAndEdgeListGraph G, concepts::LvalueVertexPropertyMap<G> PlanarEmbedding>
void make_maximal_planar(G& g, PlanarEmbedding embedding) {
  make_maximal_planar(g, embedding, get(vertex_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_MAKE_MAXIMAL_PLANAR_H_
