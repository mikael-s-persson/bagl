// Copyright 2008
// Author: Matyas W Egyhazy
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_METRIC_TSP_APPROX_H_
#define BAGL_BAGL_METRIC_TSP_APPROX_H_

// metric_tsp_approx
// Generates an approximate tour solution for the traveling salesperson
// problem in polynomial time. The current algorithm guarantees a tour with a
// length at most as long as 2x optimal solution. The graph should have
// 'natural' (metric) weights such that the triangle inequality is maintained.
// Graphs must be fully interconnected.

// TODO:
// There are a couple of improvements that could be made.
// 1) Change implementation to lower upper bound Christofides heuristic
// 2) Implement a less restrictive TSP heuristic (one that does not rely on
//    triangle inequality).
// 3) Determine if the algorithm can be implemented without creating a new
//    graph.

#include <memory>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/graph_as_tree.h"
#include "bagl/graph_traits.h"
#include "bagl/lookup_edge.h"
#include "bagl/prim_minimum_spanning_tree.h"
#include "bagl/tree_adaptor.h"

namespace bagl {

namespace concepts {
template <typename V, typename G>
concept TSPVertexVisitor = std::copy_constructible<V> &&
    requires(const V& vis, const G& g, graph_vertex_descriptor_t<G> u) {
  vis.visit_vertex(u, g);
};

}  // namespace concepts

// Tree visitor that keeps track of a preorder traversal of a tree
// TODO: Consider migrating this to the graph_as_tree header.
// TODO: Parameterize the underlying stores so it doesn't have to be a vector.
template <typename Vertex, typename Tree>
class preorder_traverser {
 public:
  explicit preorder_traverser(std::vector<Vertex>& p) : path_(p) {}

  void preorder(Vertex n, const Tree&) { path_.push_back(n); }
  void inorder(Vertex, const Tree&) const {}
  void postorder(Vertex, const Tree&) const {}

  auto begin() const { return path_.begin(); }
  auto end() const { return path_.end(); }

 private:
  std::vector<Vertex>& path_;
};

// Forward declarations
template <typename OutItr>
class tsp_tour_visitor;
template <concepts::Graph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          std::output_iterator<graph_vertex_descriptor_t<G>> OutIter, typename Length>
class tsp_tour_len_visitor;

template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutIter>
void metric_tsp_approx_tour(const G& g, OutIter o) {
  metric_tsp_approx_from_vertex(g, *vertices(g).first, get(edge_weight, g), get(vertex_index, g), tsp_tour_visitor{o});
}

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          std::output_iterator<graph_vertex_descriptor_t<G>> OutIter>
void metric_tsp_approx_tour(const G& g, WeightMap w, OutIter o) {
  metric_tsp_approx_from_vertex(g, *vertices(g).first, w, tsp_tour_visitor{o});
}

template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutIter>
void metric_tsp_approx_tour_from_vertex(const G& g, graph_vertex_descriptor_t<G> start, OutIter o) {
  metric_tsp_approx_from_vertex(g, start, get(edge_weight, g), get(vertex_index, g), tsp_tour_visitor{o});
}

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          std::output_iterator<graph_vertex_descriptor_t<G>> OutIter>
void metric_tsp_approx_tour_from_vertex(const G& g, graph_vertex_descriptor_t<G> start, WeightMap w, OutIter o) {
  metric_tsp_approx_from_vertex(g, start, w, get(vertex_index, g), tsp_tour_visitor<OutIter>(o));
}

template <concepts::VertexListGraph G, concepts::TSPVertexVisitor<G> V>
void metric_tsp_approx(const G& g, V vis) {
  metric_tsp_approx_from_vertex(g, *vertices(g).first, get(edge_weight, g), get(vertex_index, g), vis);
}

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> Weightmap,
          concepts::TSPVertexVisitor<G> V>
void metric_tsp_approx(const G& g, Weightmap w, V vis) {
  metric_tsp_approx_from_vertex(g, *vertices(g).first, w, get(vertex_index, g), vis);
}

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> VertexIndexMap, concepts::TSPVertexVisitor<G> V>
void metric_tsp_approx(const G& g, WeightMap w, VertexIndexMap id, V vis) {
  metric_tsp_approx_from_vertex(g, *vertices(g).first, w, id, vis);
}

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::TSPVertexVisitor<G> V>
void metric_tsp_approx_from_vertex(const G& g, graph_vertex_descriptor_t<G> start, WeightMap w, V vis) {
  metric_tsp_approx_from_vertex(g, start, w, get(vertex_index, g), vis);
}

template <concepts::VertexListGraph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> VertexIndexMap, concepts::TSPVertexVisitor<G> V>
void metric_tsp_approx_from_vertex(const G& g, graph_vertex_descriptor_t<G> start, WeightMap weightmap,
                                   VertexIndexMap indexmap, V vis) {
  // Types related to the input graph (GVertex is a template parameter).
  using GVertex = graph_vertex_descriptor_t<G>;

  // We build a custom graph in this algorithm.
  using MSTImpl = adjacency_list<vec_s, vec_s, directed_s, no_property, no_property>;
  using MSTVertex = graph_vertex_descriptor_t<MSTImpl>;

  // And then re-cast it as a tree.

  // A predecessor map.
  auto pred_pmap = make_vector_property_map<GVertex>(num_vertices(g), indexmap);

  // Compute a spanning tree over the in put g.
  prim_minimum_spanning_tree(g, pred_pmap, root_vertex(start).vertex_index_map(indexmap).weight_map(weightmap));

  // Build a MST using the predecessor map from prim mst
  MSTImpl mst{num_vertices(g)};
  std::size_t cnt = 0;
  auto mst_verts = vertices(mst);
  for (auto v : pred_pmap.storage_range()) {
    if (indexmap[v] != cnt) {
      add_edge(*std::next(mst_verts.begin(), indexmap[v]), *std::next(mst_verts.begin(), cnt), mst);
    }
    ++cnt;
  }

  // Build a tree abstraction over the MST.
  auto t = make_graph_as_tree(mst, tree_root(mst),
                              make_vector_property_map<MSTVertex>(num_vertices(mst), get(vertex_index, mst)));
  using Tree = decltype(t);
  using Node = tree_node_descriptor_t<Tree>;

  // Create tour using a preorder traversal of the mst
  std::vector<Node> tour;
  preorder_traverser<Node, Tree> tvis(tour);
  traverse_tree(indexmap[start], t, tvis);

  auto g_verts = vertices(g);
  for (auto curr : tvis) {
    // TODO: This is will be O(n^2) if vertex storage of g != vec_s.
    GVertex v = *std::next(g_verts.begin(), get(vertex_index, mst)[curr]);
    vis.visit_vertex(v, g);
  }

  // Connect back to the start of the tour
  vis.visit_vertex(start, g);
}

// Default tsp tour visitor that puts the tour in an OutputIterator
template <typename OutItr>
class tsp_tour_visitor {
 public:
  explicit tsp_tour_visitor(OutItr itr) : itr_(itr) {}

  template <typename Vertex, typename Graph>
  requires std::output_iterator<OutItr, Vertex>
  void visit_vertex(Vertex v, const Graph&) { *itr_++ = v; }

 private:
  OutItr itr_;
};

// Tsp tour visitor that adds the total tour length.
template <concepts::Graph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          std::output_iterator<graph_vertex_descriptor_t<G>> OutIter, typename Length>
class tsp_tour_len_visitor {
 public:
  using Vertex = graph_vertex_descriptor_t<G>;

  tsp_tour_len_visitor(const G&, OutIter iter, Length& l, WeightMap& map)
      : iter_(iter), tourlen_(l), wmap_(map), previous_(null()) {}

  void visit_vertex(Vertex v, const G& g) {
    // If it is not the start, then there is a
    // previous vertex
    if (previous_ != null()) {
      // NOTE: For non-adjacency matrix graphs g, this bit of code
      // will be linear in the degree of previous_ or v. A better
      // solution would be to visit edges of the graph, but that
      // would require revisiting the core algorithm.
      auto [e, found] = lookup_edge(previous_, v, g);
      if (!found) {
        throw not_complete();
      }
      tourlen_ += wmap_[e];
    }

    previous_ = v;
    *iter_++ = v;
  }

 private:
  OutIter iter_;
  Length& tourlen_;
  WeightMap& wmap_;
  Vertex previous_;

  // Helper function for getting the null vertex.
  Vertex null() { return graph_traits<G>::null_vertex(); }
};

// Object generator(s)
template <typename OutIter>
auto make_tsp_tour_visitor(OutIter iter) {
  return tsp_tour_visitor<OutIter>(iter);
}

template <concepts::Graph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          std::output_iterator<graph_vertex_descriptor_t<G>> OutIter, typename Length>
auto make_tsp_tour_len_visitor(const G& g, OutIter iter, Length& l, WeightMap map) {
  return tsp_tour_len_visitor<G, WeightMap, OutIter, Length>(g, iter, l, map);
}

}  // namespace bagl

#endif  // BAGL_BAGL_METRIC_TSP_APPROX_H_
