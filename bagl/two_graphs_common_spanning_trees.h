// Copyright (C) 2012, Michele Caini.
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_TWO_GRAPHS_COMMON_SPANNING_TREES_H_
#define BAGL_BAGL_TWO_GRAPHS_COMMON_SPANNING_TREES_H_

#include <map>
#include <stack>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "bagl/connected_components.h"
#include "bagl/filtered_graph.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/undirected_dfs.h"
#include "bagl/vector_property_map.h"

namespace bagl {

//          Two Graphs Common Spanning Trees Algorithm
//      Based on academic article of Mint, Read and Tarjan
//     Efficient Algorithm for Common Spanning Tree Problem
// Electron. Lett., 28 April 1983, Volume 19, Issue 9, p.346-347

namespace common_trees_detail {

template <typename TreeMap, typename PredMap, typename DistMap, typename LowMap, typename Buffer>
struct bridges_visitor : public default_dfs_visitor {
  bridges_visitor(TreeMap tree, PredMap pred, DistMap dist, LowMap low, Buffer& buffer)
      : tree_(tree), pred_(pred), dist_(dist), low_(low), buffer_(buffer) {}

  template <typename Vertex, typename Graph>
  void initialize_vertex(const Vertex& u, const Graph& /*g*/) {
    put(pred_, u, u);
    put(dist_, u, -1);
  }

  template <typename Vertex, typename Graph>
  void discover_vertex(const Vertex& u, const Graph& /*g*/) {
    put(dist_, u, ++num_);
    put(low_, u, get(dist_, u));
  }

  template <typename Edge, typename Graph>
  void tree_edge(const Edge& e, const Graph& g) {
    put(pred_, target(e, g), source(e, g));
    put(tree_, target(e, g), e);
  }

  template <typename Edge, typename Graph>
  void back_edge(const Edge& e, const Graph& g) {
    put(low_, source(e, g), (std::min)(get(low_, source(e, g)), get(dist_, target(e, g))));
  }

  template <typename Vertex, typename Graph>
  void finish_vertex(const Vertex& u, const Graph& /*g*/) {
    Vertex parent = get(pred_, u);
    if (get(low_, u) > get(dist_, parent)) {
      buffer_.push_back(get(tree_, u));
    }
    put(low_, parent, (std::min)(get(low_, parent), get(low_, u)));
  }

  TreeMap tree_;
  PredMap pred_;
  DistMap dist_;
  LowMap low_;
  Buffer& buffer_;
  int num_ = -1;
};

template <typename Buffer>
struct cycle_finder {
  explicit cycle_finder(Buffer* buffer) : buffer_(buffer) {}
  template <typename Edge, typename Graph>
  void back_edge(const Edge& e, const Graph& /*g*/) {
    buffer_->push_back(e);
  }
  Buffer* buffer_ = nullptr;
};

template <typename DeletedMap>
struct deleted_edge_status {
  deleted_edge_status() = default;
  explicit deleted_edge_status(DeletedMap map) : map_(map) {}
  template <typename Edge>
  bool operator()(const Edge& e) const {
    return (!get(map_, e));
  }
  DeletedMap map_;
};

template <typename InLMap>
struct in_l_edge_status {
  in_l_edge_status() = default;
  explicit in_l_edge_status(InLMap map) : map_(map) {}
  template <typename Edge>
  bool operator()(const Edge& e) const {
    return get(map_, e);
  }
  InLMap map_;
};

template <typename Graph, typename VertexIndex, typename EdgeCont, typename EdgeIndex, typename Func, typename Seq,
          typename Map>
void rec_two_graphs_common_spanning_trees(const Graph& i_g, const VertexIndex& i_g_vindex,
                                          const EdgeCont& i_g_id_to_edge, const EdgeIndex& i_g_eindex, Map ai_g_in_l,
                                          Map di_g, const Graph& v_g, const VertexIndex& v_g_vindex,
                                          const EdgeCont& v_g_id_to_edge, const EdgeIndex& v_g_eindex, Map av_g_in_l,
                                          Map dv_g, Func func, Seq in_l) {
  using edge_descriptor = graph_edge_descriptor_t<Graph>;

  //
  //  [ Michele Caini ]
  //
  //  Using the condition (edges != 0) leads to the accidental submission
  //  of sub-graphs ((V-1+1)-fake-tree, named here fat-tree).
  //  Remove this condition is a workaround for the problem of fat-trees.
  //  Please do not add that condition, even if it improves performance.
  //
  //  Here is proposed the previous guard (that was wrong):
  //     for(std::size_t i = 0; (i < in_l.size()) && (edges != 0); ++i)
  //
  int edges_not_in_l = num_vertices(i_g) - 1;
  for (std::size_t i = 0; i < in_l.size(); ++i) {
    if (in_l[i]) {
      --edges_not_in_l;
    }
  }

  if (edges_not_in_l < 0) {
    return;
  }
  if (edges_not_in_l == 0) {
    func(in_l);
    return;
  }

  auto i_g_vertex_color = vector_property_map(num_vertices(i_g), i_g_vindex, default_color_type::white_color);
  auto i_g_edge_color = vector_property_map(num_edges(i_g), i_g_eindex, default_color_type::white_color);
  auto v_g_vertex_color = vector_property_map(num_vertices(v_g), v_g_vindex, default_color_type::white_color);
  auto v_g_edge_color = vector_property_map(num_edges(v_g), v_g_eindex, default_color_type::white_color);

  std::vector<edge_descriptor> i_g_buf;
  std::vector<edge_descriptor> v_g_buf;
  bool found = false;

  std::size_t m = in_l.size();
  for (std::size_t j = 0; j < in_l.size() && !found; ++j) {
    if (in_l[j] || get(di_g, i_g_id_to_edge[j]) || get(dv_g, v_g_id_to_edge[j])) {
      continue;
    }

    put(ai_g_in_l, i_g_id_to_edge[j], true);
    put(av_g_in_l, v_g_id_to_edge[j], true);

    undirected_dfs(filtered_graph(i_g, in_l_edge_status(ai_g_in_l)), make_dfs_visitor(cycle_finder(&i_g_buf)),
                   i_g_vertex_color.ref(), i_g_edge_color.ref());
    undirected_dfs(filtered_graph(v_g, in_l_edge_status(av_g_in_l)), make_dfs_visitor(cycle_finder(&v_g_buf)),
                   v_g_vertex_color.ref(), v_g_edge_color.ref());

    if (i_g_buf.empty() && v_g_buf.empty()) {
      in_l[j] = true;
      found = true;
      m = j;
      break;
    }
    i_g_buf.clear();
    v_g_buf.clear();
    put(ai_g_in_l, i_g_id_to_edge[j], false);
    put(av_g_in_l, v_g_id_to_edge[j], false);
  }

  if (!found) {
    return;
  }

  std::vector<edge_descriptor> i_g_buf_copy;
  std::vector<edge_descriptor> v_g_buf_copy;
  for (std::size_t j = 0; j < in_l.size(); ++j) {
    if (in_l[j] || get(di_g, i_g_id_to_edge[j]) || get(dv_g, v_g_id_to_edge[j])) {
      continue;
    }

    put(ai_g_in_l, i_g_id_to_edge[j], true);
    put(av_g_in_l, v_g_id_to_edge[j], true);

    undirected_dfs(filtered_graph(i_g, in_l_edge_status(ai_g_in_l)), make_dfs_visitor(cycle_finder(&i_g_buf)),
                   i_g_vertex_color.ref(), i_g_edge_color.ref());
    undirected_dfs(filtered_graph(v_g, in_l_edge_status(av_g_in_l)), make_dfs_visitor(cycle_finder(&v_g_buf)),
                   v_g_vertex_color.ref(), v_g_edge_color.ref());

    if (!i_g_buf.empty() || !v_g_buf.empty()) {
      i_g_buf.clear();
      v_g_buf.clear();
      put(di_g, i_g_id_to_edge[j], true);
      put(dv_g, v_g_id_to_edge[j], true);
      i_g_buf_copy.push_back(i_g_id_to_edge[j]);
      v_g_buf_copy.push_back(v_g_id_to_edge[j]);
    }

    put(ai_g_in_l, i_g_id_to_edge[j], false);
    put(av_g_in_l, v_g_id_to_edge[j], false);
  }

  // REC
  rec_two_graphs_common_spanning_trees(i_g, i_g_vindex, i_g_id_to_edge, i_g_eindex, ai_g_in_l.ref(), di_g, v_g,
                                       v_g_vindex, v_g_id_to_edge, v_g_eindex, ai_g_in_l.ref(), dv_g, func, in_l);

  while (!i_g_buf_copy.empty()) {
    put(di_g, i_g_buf_copy.back(), false);
    put(dv_g, v_g_id_to_edge[get(i_g_eindex, i_g_buf_copy.back())], false);
    i_g_buf_copy.pop_back();
  }
  while (!v_g_buf_copy.empty()) {
    put(dv_g, v_g_buf_copy.back(), false);
    put(di_g, i_g_id_to_edge[get(i_g_eindex, v_g_buf_copy.back())], false);
    v_g_buf_copy.pop_back();
  }

  in_l[m] = false;
  put(ai_g_in_l, i_g_id_to_edge[m], false);
  put(av_g_in_l, v_g_id_to_edge[m], false);

  put(di_g, i_g_id_to_edge[m], true);
  put(dv_g, v_g_id_to_edge[m], true);

  auto i_g_tree_map = vector_property_map(num_vertices(i_g), i_g_vindex, edge_descriptor{});
  auto v_g_tree_map = vector_property_map(num_vertices(v_g), v_g_vindex, edge_descriptor{});
  auto i_g_pred_map = vector_property_map(num_vertices(i_g), i_g_vindex, graph_traits<Graph>::null_vertex());
  auto v_g_pred_map = vector_property_map(num_vertices(v_g), v_g_vindex, graph_traits<Graph>::null_vertex());
  auto i_g_dist_map = vector_property_map(num_vertices(i_g), i_g_vindex, int{0});
  auto v_g_dist_map = vector_property_map(num_vertices(v_g), v_g_vindex, int{0});
  auto i_g_low_map = vector_property_map(num_vertices(i_g), i_g_vindex, int{0});
  auto v_g_low_map = vector_property_map(num_vertices(v_g), v_g_vindex, int{0});

  bridges_visitor i_g_vis(i_g_tree_map.ref(), i_g_pred_map.ref(), i_g_dist_map.ref(), i_g_low_map.ref(), i_g_buf);
  bridges_visitor v_g_vis(v_g_tree_map.ref(), v_g_pred_map.ref(), v_g_dist_map.ref(), v_g_low_map.ref(), v_g_buf);

  undirected_dfs(filtered_graph(i_g, deleted_edge_status(di_g)), i_g_vis, i_g_vertex_color.ref(), i_g_edge_color.ref());
  undirected_dfs(filtered_graph(v_g, deleted_edge_status(dv_g)), v_g_vis, v_g_vertex_color.ref(), v_g_edge_color.ref());

  found = false;
  std::vector<edge_descriptor> i_g_buf_tmp;
  std::vector<edge_descriptor> v_g_buf_tmp;
  while (!i_g_buf.empty() && !found) {
    if (in_l[get(i_g_eindex, i_g_buf.back())]) {
      i_g_buf.pop_back();
      continue;
    }

    put(ai_g_in_l, i_g_buf.back(), true);
    put(av_g_in_l, v_g_id_to_edge[get(i_g_eindex, i_g_buf.back())], true);

    undirected_dfs(filtered_graph(i_g, in_l_edge_status(ai_g_in_l.ref())), make_dfs_visitor(cycle_finder(&i_g_buf_tmp)),
                   i_g_vertex_color.ref(), i_g_edge_color.ref());
    undirected_dfs(filtered_graph(v_g, in_l_edge_status(av_g_in_l.ref())), make_dfs_visitor(cycle_finder(&v_g_buf_tmp)),
                   v_g_vertex_color.ref(), v_g_edge_color.ref());

    if (!i_g_buf_tmp.empty() || !v_g_buf_tmp.empty()) {
      found = true;
    } else {
      i_g_buf_tmp.clear();
      v_g_buf_tmp.clear();
      i_g_buf_copy.push_back(i_g_buf.back());
    }

    put(ai_g_in_l, i_g_buf.back(), false);
    put(av_g_in_l, v_g_id_to_edge[get(i_g_eindex, i_g_buf.back())], false);
    i_g_buf.pop_back();
  }
  while (!v_g_buf.empty() && !found) {
    if (in_l[get(v_g_eindex, v_g_buf.back())]) {
      v_g_buf.pop_back();
      continue;
    }

    put(av_g_in_l, v_g_buf.back(), true);
    put(ai_g_in_l, i_g_id_to_edge[get(i_g_eindex, v_g_buf.back())], true);

    undirected_dfs(filtered_graph(i_g, in_l_edge_status(ai_g_in_l.ref())), make_dfs_visitor(cycle_finder(&i_g_buf_tmp)),
                   i_g_vertex_color.ref(), i_g_edge_color.ref());
    undirected_dfs(filtered_graph(v_g, in_l_edge_status(av_g_in_l.ref())), make_dfs_visitor(cycle_finder(&v_g_buf_tmp)),
                   v_g_vertex_color.ref(), v_g_edge_color.ref());

    if (!i_g_buf_tmp.empty() || !v_g_buf_tmp.empty()) {
      found = true;
    } else {
      i_g_buf_tmp.clear();
      v_g_buf_tmp.clear();
      v_g_buf_copy.push_back(v_g_buf.back());
    }

    put(av_g_in_l, v_g_buf.back(), false);
    put(ai_g_in_l, i_g_id_to_edge[get(v_g_eindex, v_g_buf.back())], false);

    v_g_buf.pop_back();
  }

  if (!found) {
    while (!i_g_buf_copy.empty()) {
      in_l[get(i_g_eindex, i_g_buf_copy.back())] = true;
      put(ai_g_in_l, i_g_buf_copy.back(), true);
      put(av_g_in_l, v_g_id_to_edge[get(i_g_eindex, i_g_buf_copy.back())], true);
      i_g_buf.push_back(i_g_buf_copy.back());
      i_g_buf_copy.pop_back();
    }
    while (!v_g_buf_copy.empty()) {
      in_l[get(v_g_eindex, v_g_buf_copy.back())] = true;
      put(av_g_in_l, v_g_buf_copy.back(), true);
      put(ai_g_in_l, i_g_id_to_edge[get(v_g_eindex, v_g_buf_copy.back())], true);
      v_g_buf.push_back(v_g_buf_copy.back());
      v_g_buf_copy.pop_back();
    }

    // REC
    rec_two_graphs_common_spanning_trees(i_g, i_g_vindex, i_g_id_to_edge, i_g_eindex, ai_g_in_l.ref(), di_g, v_g,
                                         v_g_vindex, v_g_id_to_edge, v_g_eindex, ai_g_in_l.ref(), dv_g, func, in_l);

    while (!i_g_buf.empty()) {
      in_l[get(i_g_eindex, i_g_buf.back())] = false;
      put(ai_g_in_l, i_g_buf.back(), false);
      put(av_g_in_l, v_g_id_to_edge[get(i_g_eindex, i_g_buf.back())], false);
      i_g_buf.pop_back();
    }
    while (!v_g_buf.empty()) {
      in_l[get(v_g_eindex, v_g_buf.back())] = false;
      put(av_g_in_l, v_g_buf.back(), false);
      put(ai_g_in_l, i_g_id_to_edge[get(v_g_eindex, v_g_buf.back())], false);
      v_g_buf.pop_back();
    }
  }

  put(di_g, i_g_id_to_edge[m], false);
  put(dv_g, v_g_id_to_edge[m], false);
}

}  // namespace common_trees_detail

template <typename Coll, std::ranges::random_access_range Seq>
struct tree_collector {
 public:
  static_assert(std::is_convertible_v<decltype(std::declval<Seq>()[0]), bool>);

  explicit tree_collector(Coll& seqs) : seqs_(seqs) {}

  void operator()(Seq seq) { seqs_.emplace_back(std::move(seq)); }

 private:
  Coll& seqs_;
};

template <concepts::VertexAndEdgeListGraph Graph, concepts::ReadableVertexIndexMap<Graph> VertexIndex,
          concepts::ReadablePropertyMap<std::size_t> EdgeOrder, concepts::ReadableEdgeIndexMap<Graph> EdgeIndex,
          std::ranges::random_access_range Seq, std::invocable<Seq&> Func>
requires concepts::IncidenceGraph<Graph> && std::ranges::output_range<Seq, bool> &&
    std::convertible_to<property_traits_value_t<EdgeOrder>, graph_edge_descriptor_t<Graph>>
void two_graphs_common_spanning_trees(const Graph& i_g, VertexIndex i_g_vindex, EdgeOrder i_g_id_to_edge,
                                      EdgeIndex i_g_eindex, const Graph& v_g, VertexIndex v_g_vindex,
                                      EdgeOrder v_g_id_to_edge, EdgeIndex v_g_eindex, Func func, Seq in_l) {
  using edge_descriptor = graph_edge_descriptor_t<Graph>;

  static_assert(is_undirected_graph_v<Graph>);

  if (num_vertices(i_g) != num_vertices(v_g)) {
    return;
  }

  if (in_l.size() != num_edges(i_g) || in_l.size() != num_edges(v_g)) {
    return;
  }

  std::vector<edge_descriptor> i_g_buf;
  std::vector<edge_descriptor> v_g_buf;

  auto i_g_tree_map = vector_property_map(num_vertices(i_g), i_g_vindex, edge_descriptor{});
  auto v_g_tree_map = vector_property_map(num_vertices(v_g), v_g_vindex, edge_descriptor{});
  auto i_g_pred_map = vector_property_map(num_vertices(i_g), i_g_vindex, graph_traits<Graph>::null_vertex());
  auto v_g_pred_map = vector_property_map(num_vertices(v_g), v_g_vindex, graph_traits<Graph>::null_vertex());
  auto i_g_dist_map = vector_property_map(num_vertices(i_g), i_g_vindex, int{0});
  auto v_g_dist_map = vector_property_map(num_vertices(v_g), v_g_vindex, int{0});
  auto i_g_low_map = vector_property_map(num_vertices(i_g), i_g_vindex, int{0});
  auto v_g_low_map = vector_property_map(num_vertices(v_g), v_g_vindex, int{0});

  common_trees_detail::bridges_visitor i_g_vis(i_g_tree_map.ref(), i_g_pred_map.ref(), i_g_dist_map.ref(),
                                               i_g_low_map.ref(), i_g_buf);
  common_trees_detail::bridges_visitor v_g_vis(v_g_tree_map.ref(), v_g_pred_map.ref(), v_g_dist_map.ref(),
                                               v_g_low_map.ref(), v_g_buf);

  auto i_g_vertex_color = vector_property_map(num_vertices(i_g), i_g_vindex, default_color_type::white_color);
  auto i_g_edge_color = vector_property_map(num_edges(i_g), i_g_eindex, default_color_type::white_color);
  auto v_g_vertex_color = vector_property_map(num_vertices(v_g), v_g_vindex, default_color_type::white_color);
  auto v_g_edge_color = vector_property_map(num_edges(v_g), v_g_eindex, default_color_type::white_color);

  undirected_dfs(i_g, i_g_vis, i_g_vertex_color.ref(), i_g_edge_color.ref());
  undirected_dfs(v_g, v_g_vis, v_g_vertex_color.ref(), v_g_edge_color.ref());

  while (!i_g_buf.empty()) {
    in_l[get(i_g_eindex, i_g_buf.back())] = true;
    i_g_buf.pop_back();
  }
  while (!v_g_buf.empty()) {
    in_l[get(v_g_eindex, v_g_buf.back())] = true;
    v_g_buf.pop_back();
  }

  auto ai_g_in_l = vector_property_map(num_edges(i_g), i_g_eindex, false);
  auto av_g_in_l = vector_property_map(num_edges(v_g), v_g_eindex, false);

  for (std::size_t i = 0; i < in_l.size(); ++i) {
    if (in_l[i]) {
      put(ai_g_in_l, i_g_id_to_edge[i], true);
      put(av_g_in_l, v_g_id_to_edge[i], true);
    } else {
      put(ai_g_in_l, i_g_id_to_edge[i], false);
      put(av_g_in_l, v_g_id_to_edge[i], false);
    }
  }

  undirected_dfs(filtered_graph(i_g, common_trees_detail::in_l_edge_status(ai_g_in_l.ref())),
                 make_dfs_visitor(common_trees_detail::cycle_finder(&i_g_buf)), i_g_vertex_color.ref(),
                 i_g_edge_color.ref());
  undirected_dfs(filtered_graph(v_g, common_trees_detail::in_l_edge_status(av_g_in_l.ref())),
                 make_dfs_visitor(common_trees_detail::cycle_finder(&v_g_buf)), v_g_vertex_color.ref(),
                 v_g_edge_color.ref());

  if (!i_g_buf.empty() || !v_g_buf.empty()) {
    return;
  }

  auto di_g = vector_property_map(num_edges(i_g), i_g_eindex, false);
  auto dv_g = vector_property_map(num_edges(v_g), v_g_eindex, false);

  for (std::size_t j = 0; j < in_l.size(); ++j) {
    if (in_l[j]) {
      continue;
    }
    put(ai_g_in_l, i_g_id_to_edge[j], true);
    put(av_g_in_l, v_g_id_to_edge[j], true);

    undirected_dfs(filtered_graph(i_g, common_trees_detail::in_l_edge_status(ai_g_in_l.ref())),
                   make_dfs_visitor(common_trees_detail::cycle_finder(&i_g_buf)), i_g_vertex_color.ref(),
                   i_g_edge_color.ref());
    undirected_dfs(filtered_graph(v_g, common_trees_detail::in_l_edge_status(av_g_in_l.ref())),
                   make_dfs_visitor(common_trees_detail::cycle_finder(&v_g_buf)), v_g_vertex_color.ref(),
                   v_g_edge_color.ref());

    if (!i_g_buf.empty() || !v_g_buf.empty()) {
      i_g_buf.clear();
      v_g_buf.clear();
      put(di_g, i_g_id_to_edge[j], true);
      put(dv_g, v_g_id_to_edge[j], true);
    }

    put(ai_g_in_l, i_g_id_to_edge[j], false);
    put(av_g_in_l, v_g_id_to_edge[j], false);
  }

  int cc = 0;

  auto i_g_com_map = vector_property_map(num_vertices(i_g), i_g_vindex, int{0});
  auto v_g_com_map = vector_property_map(num_vertices(v_g), v_g_vindex, int{0});
  cc += connected_components(filtered_graph(i_g, common_trees_detail::deleted_edge_status(di_g.ref())),
                             i_g_com_map.ref());
  cc += connected_components(filtered_graph(v_g, common_trees_detail::deleted_edge_status(dv_g.ref())),
                             v_g_com_map.ref());

  if (cc != 2) {
    return;
  }

  // REC
  common_trees_detail::rec_two_graphs_common_spanning_trees(
      i_g, i_g_vindex, i_g_id_to_edge, i_g_eindex, ai_g_in_l.ref(), di_g.ref(), v_g, v_g_vindex, v_g_id_to_edge,
      v_g_eindex, ai_g_in_l.ref(), dv_g.ref(), func, in_l);
}

template <concepts::EdgeListGraph Graph, concepts::ReadableVertexIndexMap<Graph> VertexIndex,
          concepts::ReadablePropertyMap<std::size_t> EdgeOrder, typename Func, typename Seq>
requires std::convertible_to<property_traits_value_t<EdgeOrder>, graph_edge_descriptor_t<Graph>>
void two_graphs_common_spanning_trees(const Graph& i_g, VertexIndex i_g_vindex, EdgeOrder i_g_id_to_edge,
                                      const Graph& v_g, VertexIndex v_g_vindex, EdgeOrder v_g_id_to_edge, Func func,
                                      Seq in_l) {
  using edge_descriptor = graph_edge_descriptor_t<Graph>;

  std::unordered_map<edge_descriptor, std::size_t> i_g_eindex_store;
  i_g_eindex_store.reserve(num_edges(i_g));
  for (std::size_t i = 0; i < num_edges(i_g); ++i) {
    i_g_eindex_store[get(i_g_id_to_edge, i)] = i;
  }
  std::unordered_map<edge_descriptor, std::size_t> v_g_eindex_store;
  v_g_eindex_store.reserve(num_edges(v_g));
  for (std::size_t i = 0; i < num_edges(v_g); ++i) {
    v_g_eindex_store[get(v_g_id_to_edge, i)] = i;
  }

  for (auto e : edges(i_g)) {
    if (i_g_eindex_store.find(e) == i_g_eindex_store.end()) {
      return;
    }
  }
  for (auto e : edges(v_g)) {
    if (v_g_eindex_store.find(e) == v_g_eindex_store.end()) {
      return;
    }
  }
  auto i_g_eindex = associative_property_map(i_g_eindex_store);
  auto v_g_eindex = associative_property_map(v_g_eindex_store);

  two_graphs_common_spanning_trees(i_g, i_g_vindex, i_g_id_to_edge, i_g_eindex, v_g, v_g_vindex, v_g_id_to_edge,
                                   v_g_eindex, func, in_l);
}

template <concepts::EdgeListGraph Graph, concepts::ReadableVertexIndexMap<Graph> VertexIndex,
          concepts::ReadableEdgeIndexMap<Graph> EdgeIndex, typename Func, typename Seq>
void two_graphs_common_spanning_trees(const Graph& i_g, VertexIndex i_g_vindex, EdgeIndex i_g_eindex, const Graph& v_g,
                                      VertexIndex v_g_vindex, EdgeIndex v_g_eindex, Func func, Seq in_l) {
  using edge_descriptor = graph_edge_descriptor_t<Graph>;

  auto i_g_o = vector_property_map<edge_descriptor>(num_edges(i_g));
  for (auto e : edges(i_g)) {
    put(i_g_o, get(i_g_eindex, e), e);
  }

  auto v_g_o = vector_property_map<edge_descriptor>(num_edges(v_g));
  for (auto e : edges(v_g)) {
    put(v_g_o, get(v_g_eindex, e), e);
  }

  two_graphs_common_spanning_trees(i_g, i_g_vindex, i_g_o.ref(), i_g_eindex, v_g, v_g_vindex, v_g_o.ref(), v_g_eindex,
                                   func, in_l);
}

template <concepts::EdgeListGraph Graph, concepts::ReadableVertexIndexMap<Graph> VertexIndex, typename Func,
          typename Seq>
void two_graphs_common_spanning_trees(const Graph& i_g, VertexIndex i_g_vindex, const Graph& v_g,
                                      VertexIndex v_g_vindex, Func func, Seq in_l) {
  using edge_descriptor = graph_edge_descriptor_t<Graph>;

  std::unordered_map<edge_descriptor, std::size_t> i_g_eindex_store;
  i_g_eindex_store.reserve(num_edges(i_g));
  std::size_t i = 0;
  for (auto e : edges(i_g)) {
    i_g_eindex_store[e] = i++;
  }
  std::unordered_map<edge_descriptor, std::size_t> v_g_eindex_store;
  v_g_eindex_store.reserve(num_edges(v_g));
  i = 0;
  for (auto e : edges(v_g)) {
    v_g_eindex_store[e] = i++;
  }
  auto i_g_eindex = associative_property_map(i_g_eindex_store);
  auto v_g_eindex = associative_property_map(v_g_eindex_store);

  two_graphs_common_spanning_trees(i_g, i_g_vindex, i_g_eindex, v_g, v_g_vindex, v_g_eindex, func, in_l);
}

template <concepts::EdgeListGraph Graph, typename Func, typename Seq>
void two_graphs_common_spanning_trees(const Graph& i_g, const Graph& v_g, Func func, Seq in_l) {
  two_graphs_common_spanning_trees(i_g, get(vertex_index, i_g), v_g, get(vertex_index, v_g), func, in_l);
}

}  // namespace bagl

#endif  // BAGL_BAGL_TWO_GRAPHS_COMMON_SPANNING_TREES_H_
