// Copyright 2001 Jeremy Siek
// Copyright 2004 Douglas Gregor
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_BICONNECTED_COMPONENTS_H_
#define BAGL_BAGL_BICONNECTED_COMPONENTS_H_

#include <algorithm>  // for std::min and std::max
#include <cassert>
#include <limits>
#include <stack>
#include <vector>

#include "bagl/depth_first_search.h"
#include "bagl/detail/container_generators.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"

namespace bagl {
namespace biconnected_detail {
template <typename ComponentMap, typename DiscoverTimeMap, typename LowPointMap, typename PredecessorMap,
          typename OutputIterator, typename Stack, typename ArticulationVector, typename IndexMap, typename DFSVisitor>
struct biconnected_components_visitor : public dfs_visitor<> {
  biconnected_components_visitor(ComponentMap comp, std::size_t& c, std::size_t& children_of_root, DiscoverTimeMap dtm,
                                 std::size_t& dfs_time, LowPointMap lowpt, PredecessorMap pred, OutputIterator out,
                                 Stack& S, ArticulationVector& is_articulation_point, IndexMap index_map,
                                 DFSVisitor vis)
      : comp_(comp),
        c_(c),
        children_of_root_(children_of_root),
        dtm_(dtm),
        dfs_time_(dfs_time),
        lowpt_(lowpt),
        pred_(pred),
        out_(out),
        stack_(S),
        is_articulation_point_(is_articulation_point),
        index_map_(index_map),
        dfs_vis_(vis) {}

  template <typename Vertex, typename Graph>
  void initialize_vertex(const Vertex& u, Graph& g) {
    put(pred_, u, u);
    dfs_vis_.initialize_vertex(u, g);
  }

  template <typename Vertex, typename Graph>
  void start_vertex(const Vertex& u, Graph& g) {
    children_of_root_ = 0;
    dfs_vis_.start_vertex(u, g);
  }

  template <typename Vertex, typename Graph>
  void discover_vertex(const Vertex& u, Graph& g) {
    put(dtm_, u, ++dfs_time_);
    put(lowpt_, u, get(dtm_, u));
    dfs_vis_.discover_vertex(u, g);
  }

  template <typename Edge, typename Graph>
  void examine_edge(const Edge& e, Graph& g) {
    dfs_vis_.examine_edge(e, g);
  }

  template <typename Edge, typename Graph>
  void tree_edge(const Edge& e, Graph& g) {
    auto src = source(e, g);
    auto tgt = target(e, g);

    stack_.push(e);
    put(pred_, tgt, src);
    if (get(pred_, src) == src) {
      ++children_of_root_;
    }
    dfs_vis_.tree_edge(e, g);
  }

  template <typename Edge, typename Graph>
  void back_edge(const Edge& e, Graph& g) {
    auto src = source(e, g);
    auto tgt = target(e, g);
    if (tgt != get(pred_, src)) {
      stack_.push(e);
      put(lowpt_, src, std::min(get(lowpt_, src), get(dtm_, tgt)));
    }
    dfs_vis_.back_edge(e, g);
  }

  template <typename Edge, typename Graph>
  void forward_or_cross_edge(const Edge& e, Graph& g) {
    dfs_vis_.forward_or_cross_edge(e, g);
  }

  template <typename Vertex, typename Graph>
  void finish_vertex(const Vertex& u, Graph& g) {
    Vertex parent = get(pred_, u);
    if (parent == u) {  // Root of tree is special
      is_articulation_point_[get(index_map_, u)] = (children_of_root_ > 1);
    } else {
      put(lowpt_, parent, std::min(get(lowpt_, parent), get(lowpt_, u)));
      if (get(lowpt_, u) >= get(dtm_, parent)) {
        is_articulation_point_[get(index_map_, parent)] = true;
        while (get(dtm_, source(stack_.top(), g)) >= get(dtm_, u)) {
          put(comp_, stack_.top(), c_);
          stack_.pop();
        }
        assert(source(S.top(), g) == parent);
        assert(target(S.top(), g) == u);
        put(comp_, stack_.top(), c_);
        stack_.pop();
        ++c_;
      }
    }
    if (is_articulation_point_[get(index_map_, u)]) {
      *out_++ = u;
    }
    dfs_vis_.finish_vertex(u, g);
  }

  ComponentMap comp_;
  std::size_t& c_;
  std::size_t& children_of_root_;
  DiscoverTimeMap dtm_;
  std::size_t& dfs_time_;
  LowPointMap lowpt_;
  PredecessorMap pred_;
  OutputIterator out_;
  Stack& stack_;
  ArticulationVector& is_articulation_point_;
  IndexMap index_map_;
  DFSVisitor dfs_vis_;
};

}  // namespace biconnected_detail

template <concepts::IncidenceGraph G, concepts::WritableEdgePropertyMap<G> ComponentMap,
          std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator,
          concepts::ReadableVertexPropertyMap<G> VertexIndexMap,
          concepts::ReadWriteVertexPropertyMap<G> DiscoverTimeMap, concepts::ReadWriteVertexPropertyMap<G> LowPointMap,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap, concepts::DFSVisitor<G> V>
requires concepts::VertexListGraph<G> std::pair<std::size_t, OutputIterator> biconnected_components(
    const G& g, ComponentMap comp, OutputIterator out, VertexIndexMap index_map, DiscoverTimeMap dtm, LowPointMap lowpt,
    PredecessorMap pred, V dfs_vis) {
  using Edge = graph_edge_descriptor_t<G>;

  std::size_t num_components = 0;
  std::size_t children_of_root = 0;
  std::size_t dfs_time = 0;
  std::stack<Edge> edge_stack;
  std::vector<char> is_articulation_point(num_vertices(g));

  biconnected_detail::biconnected_components_visitor vis{
      comp, num_components, children_of_root,      dtm,       dfs_time, lowpt, pred,
      out,  edge_stack,     is_articulation_point, index_map, dfs_vis};

  depth_first_search(g, vis, make_vector_property_map<default_color_type>(index_map));

  return {num_components, vis.out_};
}

template <concepts::IncidenceGraph G, concepts::WritableEdgePropertyMap<G> ComponentMap,
          std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator,
          concepts::ReadWriteVertexPropertyMap<G> DiscoverTimeMap, concepts::ReadWriteVertexPropertyMap<G> LowPointMap>
requires concepts::VertexListGraph<G> std::pair<std::size_t, OutputIterator> biconnected_components(
    const G& g, ComponentMap comp, OutputIterator out, DiscoverTimeMap dtm, LowPointMap lowpt) {
  auto index_map = get(vertex_index, g);
  return biconnected_components(g, comp, out, index_map, dtm, lowpt,
                                make_vector_property_map(num_vertices(g), index_map, graph_traits<G>::null_vertex()),
                                make_dfs_visitor(null_visitor()));
}

template <concepts::IncidenceGraph G, concepts::WritableEdgePropertyMap<G> ComponentMap,
          std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator>
requires concepts::VertexListGraph<G> std::pair<std::size_t, OutputIterator> biconnected_components(
    const G& g, ComponentMap comp, OutputIterator out) {
  auto index_map = get(vertex_index, g);
  auto dtm = make_vector_property_map(num_vertices(g), index_map, std::size_t{0});
  auto lowpt = make_vector_property_map(num_vertices(g), index_map, std::size_t{0});

  return biconnected_components(g, comp, out, index_map, dtm, lowpt);
}

template <concepts::IncidenceGraph G, concepts::WritableEdgePropertyMap<G> ComponentMap>
requires concepts::VertexListGraph<G> std::size_t biconnected_components(const G& g, ComponentMap comp) {
  return biconnected_components(g, comp, container_detail::ignore_output_iter()).first;
}

template <concepts::IncidenceGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator>
requires concepts::VertexListGraph<G> OutputIterator articulation_points(const G& g, OutputIterator out) {
  return biconnected_components(g, dummy_property_map(), out).second;
}

}  // namespace bagl

#endif  // BAGL_BAGL_BICONNECTED_COMPONENTS_H_
