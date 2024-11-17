// Copyright 2013-2024 Mikael Persson

#ifndef BAGL_BAGL_LTREE_CONTAINERS_H_
#define BAGL_BAGL_LTREE_CONTAINERS_H_

#include <iterator>
#include <memory>
#include <queue>
#include <stack>
#include <type_traits>
#include <utility>
#include <vector>

#include "bagl/detail/adjlist_containers.h"
#include "bagl/detail/adjlist_ranges.h"
#include "bagl/detail/container_generators.h"
#include "bagl/graph_selectors.h"  // for directed_s, undirected_s, bidirectional_s

// Use adjlist_detail namespace due to mixing overloads with general adj-list functions.
namespace bagl::adjlist_detail {

/*************************************************************************
 *        value-types for vertices and edges
 * **********************************************************************/

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
struct ltree_vertex_stored_type;  // forward declaration.

// this is fine because boost-containers allow incomplete types:
template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
struct ltree_vertex_config {
  using stored_type = ltree_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>;
  using container = typename container_gen<VertexListS, stored_type>::type;
  using value_type = typename container::value_type;
  using descriptor = typename container_detail::select_descriptor<container>::type;

  static descriptor null_vertex() { return container_detail::null_desc<descriptor>::value(); }
};

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
struct ltree_edge_config {
  using vertex_descriptor =
      typename ltree_vertex_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>::descriptor;

  using stored_type = std::pair<vertex_descriptor, EdgeProperties>;
  using container = typename container_gen<OutEdgeListS, stored_type>::type;
  using value_type = typename container::value_type;
  using raw_descriptor = typename container_detail::select_descriptor<container>::type;
  using descriptor = container_detail::edge_desc<vertex_descriptor, raw_descriptor>;

  static descriptor null_edge() { return container_detail::null_desc<descriptor>::value(); }
};

/*************************************************************************
 *        vertex values, as stored in the vertex containers.
 * **********************************************************************/

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
struct ltree_vertex_stored_type {
  using self = ltree_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>;
  using directed_tag = DirectedS;
  using Config = ltree_edge_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>;
  using edge_container = typename Config::container;
  using edge_descriptor = typename Config::descriptor;
  using in_edge_iterator = edge_descriptor*;

  VertexProperties data;
  edge_container out_edges;
  edge_descriptor in_edge;

  edge_container& get_out_edges() { return out_edges; }
  const edge_container& get_out_edges() const { return out_edges; }

  auto get_in_edges() { return std::ranges::subrange(&in_edge, &in_edge + 1); }
  auto get_in_edges() const { return std::ranges::subrange(&in_edge, &in_edge + 1); }

  ltree_vertex_stored_type() : data(), out_edges(), in_edge() {}
  template <typename VProp>
  requires std::constructible_from<VertexProperties, VProp&&>
  explicit ltree_vertex_stored_type(VProp&& aData) : data(std::forward<VProp>(aData)), out_edges(), in_edge() {}

  ltree_vertex_stored_type(const ltree_vertex_stored_type&) = delete;
  ltree_vertex_stored_type& operator=(const ltree_vertex_stored_type&) = delete;
  ltree_vertex_stored_type(ltree_vertex_stored_type&&) noexcept = default;
  ltree_vertex_stored_type& operator=(ltree_vertex_stored_type&&) noexcept = default;
  ~ltree_vertex_stored_type() = default;
};

template <typename VertexListS, typename OutEdgeListS, typename VertexProperties, typename EdgeProperties>
struct ltree_vertex_stored_type<VertexListS, OutEdgeListS, directed_s, VertexProperties, EdgeProperties> {
  using self = ltree_vertex_stored_type<VertexListS, OutEdgeListS, directed_s, VertexProperties, EdgeProperties>;
  using directed_tag = directed_s;
  using Config = ltree_edge_config<VertexListS, OutEdgeListS, directed_s, VertexProperties, EdgeProperties>;
  using edge_container = typename Config::container;
  using edge_descriptor = typename Config::descriptor;
  using in_edge_iterator = int*;

  VertexProperties data;
  edge_container out_edges;

  edge_container& get_out_edges() { return out_edges; }
  const edge_container& get_out_edges() const { return out_edges; }

  auto get_in_edges() { return std::ranges::empty_view<edge_descriptor>{}; }
  auto get_in_edges() const { return std::ranges::empty_view<edge_descriptor>{}; }

  ltree_vertex_stored_type() : data(), out_edges() {}
  template <typename VProp>
  requires std::constructible_from<VertexProperties, VProp&&>
  explicit ltree_vertex_stored_type(VProp&& aData) : data(std::forward<VProp>(aData)), out_edges() {}

  ltree_vertex_stored_type(const ltree_vertex_stored_type&) = delete;
  ltree_vertex_stored_type& operator=(const ltree_vertex_stored_type&) = delete;
  ltree_vertex_stored_type(ltree_vertex_stored_type&&) noexcept = default;
  ltree_vertex_stored_type& operator=(ltree_vertex_stored_type&&) noexcept = default;
  ~ltree_vertex_stored_type() = default;
};

// NOTE: ltree_add_edge = adjlist_add_edge
// NOTE: ltree_find_edge_to = adjlist_find_edge_to
// NOTE: ltree_add_vertex = adjlist_add_vertex

/*************************************************************************
 *        helper functions for erasing in-edges / out-edges of a vertex
 * **********************************************************************/
// NOTE: Time complexities:
//               vec_s      pool_s      list_s    (multi)set_s      unordered_(multi)set_s
// directed_s:  O(E/V)      O(E/V)      O(E/V)    O(log(E/V))       O(1)
// bidir:     O((E/V)^2)    O(E/V)      O(E/V)    O(log(E/V))       O(1)

// NOTE: ltree_erase_edges_to = adjlist_erase_edges_to
//  As long as these special versions of update_in_edge_id are used (by ADL):

template <typename VertexListS, typename OutEdgeListS, typename VertexProperties, typename EdgeProperties,
          typename VertexDesc>
void adjlist_update_in_edge_id(
    ltree_vertex_stored_type<VertexListS, OutEdgeListS, directed_s, VertexProperties, EdgeProperties>& vp, VertexDesc v,
    std::size_t old_id, std::size_t new_id) {}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties, typename VertexDesc>
void adjlist_update_in_edge_id(
    ltree_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& vp, VertexDesc v,
    std::size_t old_id, std::size_t new_id) {
  using VertexValue = ltree_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>;
  using EdgeDesc = typename VertexValue::edge_descriptor;

  if (vp.in_edge == EdgeDesc::null_value()) {
    return;
  }
  if ((vp.in_edge.source == v) && (vp.in_edge.edge_id == old_id)) {
    vp.in_edge.edge_id = new_id;
  }
}

/*************************************************************************
 *        functions for erasing an edge
 * **********************************************************************/
// NOTE: Time complexities:
//               vec_s      pool_s      list_s
// directed_s:    O(1)        O(1)        O(1)
// bidir:       O(E/V)      O(E/V)      O(E/V)   (because of in-edge erasure (vector storage))

// for OutEdgeListS = listS, setS, ...
template <typename Container, typename EdgeDesc, typename VertexCont, typename VertexDesc>
void ltree_erase_edge(Container& cont, EdgeDesc e, VertexCont& vcont, VertexDesc v) {
  adjlist_erase_edge(cont, e, vcont, v);
}

// for OutEdgeListS = vecS
template <typename ValueType, typename EdgeDesc, typename VertexCont, typename VertexDesc>
void ltree_erase_edge(std::vector<ValueType>& cont, EdgeDesc e, VertexCont& vertex_cont, VertexDesc v) {
  using std::swap;

  auto it = cont.begin() + e.edge_id;
  auto it_last = cont.end();
  --it_last;
  if (it != it_last) {
    swap(*it, *it_last);
    // If this graph has in-edge references, then they must be updated.
    adjlist_update_in_edge_id(container_detail::desc_to_value(vertex_cont, it->first), v, it_last - cont.begin(),
                              it - cont.begin());
  }
  cont.erase(it_last, cont.end());
}

/*************************************************************************
 *        functions for adding an edge
 * **********************************************************************/
// NOTE: Time complexities:
//               vec_s      pool_s      list_s
// directed_s:    O(1)        O(1)        O(1)
// bidir:         O(1)        O(1)        O(1)

template <typename VertexListS, typename OutEdgeListS, typename VertexProperties, typename EdgeProperties,
          typename EdgeDesc>
void ltree_add_in_edge(
    ltree_vertex_stored_type<VertexListS, OutEdgeListS, directed_s, VertexProperties, EdgeProperties>& vp, EdgeDesc e) {
}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties, typename EdgeDesc>
void ltree_add_in_edge(
    ltree_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& vp, EdgeDesc e) {
  vp.in_edge = e;
}

/*************************************************************************
 *        functions for erasing a vertex (including updating edges of surrounding vertices
 * **********************************************************************/
// NOTE: The function 'ltree_erase_vertex' works for all graph types.
// NOTE: Time complexities:
//               vec_s      pool_s      list_s
// directed_s:     O(E)        O(1)        O(1)
// bidir:     O((E/V)^2)      O(1)        O(1)

template <typename Container, typename VDescContainer>
void ltree_erase_vertices(Container& cont, VDescContainer& v_list) {
  for (auto v : v_list) {
    adjlist_erase_vertex(cont, v);
  }
}

// NOTE: ltree_update_out_edges_impl = adjlist_update_out_edges_impl

// O(E)
template <typename DirectedS, typename ValueType>
std::enable_if_t<std::is_same_v<DirectedS, directed_s>> ltree_update_out_edges(std::vector<ValueType>& cont,
                                                                               std::size_t old_v_id,
                                                                               std::size_t new_v_id) {
  adjlist_update_out_edges<DirectedS>(cont, old_v_id, new_v_id);
}

// O(E/V)
template <typename DirectedS, typename ValueType>
std::enable_if_t<!std::is_same_v<DirectedS, directed_s>> ltree_update_out_edges(std::vector<ValueType>& cont,
                                                                                std::size_t old_v_id,
                                                                                std::size_t new_v_id) {
  using EdgeDesc = typename ValueType::edge_descriptor;
  // first, update in-edge vertices
  if (cont[new_v_id].in_edge != EdgeDesc::null_value()) {
    EdgeDesc& e_in = cont[new_v_id].in_edge;
    ValueType& up = cont[e_in.source];
    adjlist_update_out_edges_impl(up.get_out_edges(), old_v_id, new_v_id, e_in.edge_id);
  }
  // second, update out-edge vertices
  for (auto ei = cont[new_v_id].get_out_edges().begin(); ei != cont[new_v_id].get_out_edges().end(); ++ei) {
    if (!container_detail::is_elem_valid(*ei)) {
      continue;
    }
    auto& [e_target, e_prop] = container_detail::get_value(*ei);
    ValueType& wp = cont[e_target];
    EdgeDesc& e_in = wp.in_edge;
    if (e_in != EdgeDesc::null_value()) {
      if ((e_in.source == old_v_id) &&
          (e_in.edge_id == container_detail::iterator_to_desc(cont[new_v_id].get_out_edges(), ei))) {
        e_in.source = new_v_id;
        break;
      }
    }
  }
}

// O(1)
template <typename DirectedS, typename ValueType>
std::enable_if_t<std::is_same_v<DirectedS, directed_s>> ltree_invalidate_in_edges(std::vector<ValueType>& /*unused*/,
                                                                                  std::size_t /*unused*/) {}

// O(1)
template <typename DirectedS, typename ValueType>
std::enable_if_t<!std::is_same_v<DirectedS, directed_s>> ltree_invalidate_in_edges(std::vector<ValueType>& cont,
                                                                                   std::size_t u) {
  using EdgeDesc = typename ValueType::edge_descriptor;

  for (auto& e : cont[u].get_out_edges()) {
    if (!container_detail::is_elem_valid(e)) {
      continue;
    }
    cont[container_detail::get_value(e).first].in_edge = EdgeDesc::null_value();
  }
}

template <typename ValueType>
void ltree_erase_vertices(std::vector<ValueType>& cont, std::vector<std::size_t>& v_list) {
  using DirectedS = typename ValueType::directed_tag;
  using std::swap;

  for (auto v_it = v_list.begin(); v_it != v_list.end(); ++v_it) {
    // NOTE: First, invalidate the in-edge references of the children, then destroy out-edge list.
    ltree_invalidate_in_edges<DirectedS>(cont, *v_it);
    auto it = cont.begin() + *v_it;
    auto it_last = cont.end();
    --it_last;
    if (it == it_last) {
      cont.erase(it_last);
      continue;
    }
    swap(*it, *it_last);
    std::size_t old_id = it_last - cont.begin();
    std::size_t new_id = it - cont.begin();
    cont.erase(it_last);
    ltree_update_out_edges<DirectedS>(cont, old_id, new_id);
    for (auto v_it2 = v_it + 1; v_it2 != v_list.end(); ++v_it2) {
      if (*v_it2 == old_id) {
        *v_it2 = new_id;
      }
    }
  }
}

template <typename ValueType>
void ltree_erase_vertices(container_detail::pooled_vector<ValueType>& cont, std::vector<std::size_t>& v_list) {
  for (std::size_t& it : v_list) {
    adjlist_erase_vertex(cont, it);
  }
}

// O(1) for bidirectional_s
// O(E) for directed_s
template <typename DirectedS, typename Container, typename VertexValue, typename Vertex>
auto ltree_get_in_edge(Container& vcont, VertexValue& vp, Vertex v) {
  if constexpr (std::is_same_v<DirectedS, directed_s>) {
    using EdgeDesc = typename VertexValue::edge_descriptor;
    for (auto it = vcont.begin(); it != vcont.end(); ++it) {
      if (container_detail::is_elem_valid(*it)) {
        VertexValue& up = container_detail::get_value(*it);
        for (auto ei = up.get_out_edges().begin(); ei != up.get_out_edges().end(); ++ei) {
          if (container_detail::is_elem_valid(*ei) && (container_detail::get_value(*ei).first == v)) {
            return EdgeDesc(container_detail::iterator_to_desc(vcont, it),
                            container_detail::iterator_to_desc(up.get_out_edges(), ei));
          }
        }
      }
    }
    return EdgeDesc::null_value();
  } else {
    return vp.in_edge;
  }
}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
struct ltree_vertex_container {
  using self = ltree_vertex_container<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>;

  using VConfig = ltree_vertex_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>;
  using vertex_container = typename VConfig::container;
  using vertices_size_type = typename vertex_container::size_type;
  using vertex_descriptor = typename VConfig::descriptor;
  using vertex_stored_type = typename VConfig::stored_type;
  using vertex_value_type = typename VConfig::value_type;

  using EConfig = ltree_edge_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>;
  using edge_container = typename EConfig::container;
  using edges_size_type = typename edge_container::size_type;
  using edge_descriptor = typename EConfig::descriptor;
  using edge_stored_type = typename EConfig::stored_type;
  using edge_value_type = typename EConfig::value_type;

  mutable vertex_container m_vertices;
  vertex_descriptor m_root;

  ltree_vertex_container() : m_vertices(), m_root(container_detail::null_desc<vertex_descriptor>::value()) {}
  ~ltree_vertex_container() { clear(); }

  ltree_vertex_container(const ltree_vertex_container&) = delete;
  ltree_vertex_container& operator=(const ltree_vertex_container&) = delete;

  void swap(ltree_vertex_container& rhs) {
    using std::swap;
    m_vertices.swap(rhs.m_vertices);
    // swap(m_vertices, rhs.m_vertices);
    swap(m_root, rhs.m_root);
  }

  ltree_vertex_container(ltree_vertex_container&& rhs) noexcept
      : m_vertices(), m_root(container_detail::null_desc<vertex_descriptor>::value()) {
    swap(rhs);
  }
  ltree_vertex_container& operator=(ltree_vertex_container&& rhs) noexcept {
    swap(rhs);
    return *this;
  }

  /********************************************************************************
   * NOTE NOTE NOTE - FUNCTIONS THAT ARE THE SAME AS ADJ-LIST-BC - NOTE NOTE NOTE *
   * ******************************************************************************/

  std::size_t size() const { return m_vertices.size(); }
  std::size_t capacity() const { return m_vertices.capacity(); }

  std::size_t num_edges() const { return (m_vertices.size() == 0 ? 0 : m_vertices.size() - 1); }

  vertex_stored_type& get_stored_vertex(vertex_descriptor v) const {
    return container_detail::desc_to_value(m_vertices, v);
  }

  const auto& get_stored_edge_property(edge_descriptor e) const {
    return container_detail::desc_to_value(get_stored_vertex(e.source).get_out_edges(), e.edge_id).second;
  }
  auto& get_stored_edge_property(edge_descriptor e) {
    return container_detail::desc_to_value(get_stored_vertex(e.source).get_out_edges(), e.edge_id).second;
  }
  vertex_descriptor get_stored_edge_target(edge_descriptor e) const {
    return container_detail::desc_to_value(get_stored_vertex(e.source).get_out_edges(), e.edge_id).first;
  }

  std::size_t get_out_degree(vertex_descriptor v) const { return get_stored_vertex(v).get_out_edges().size(); }

  std::size_t get_in_degree(vertex_descriptor v) const {
    if ((v != container_detail::null_desc<vertex_descriptor>::value()) && (v != m_root)) {
      return 1;
    }
    return 0;
  }

  // NOTE: This WORKS for ALL vertex container types.
  void clear() {
    m_vertices.clear();
    m_root = container_detail::null_desc<vertex_descriptor>::value();
  }

  // NOTE: This WORKS for ALL vertex container types.
  using VRangeSelect = adjlist_select_vertex_range<VertexListS, vertex_container>;
  auto vertices() const { return VRangeSelect::create_range(m_vertices); }

  // NOTE: This WORKS for ALL vertex container types.
  // NOTE: This WORKS for ALL edge container types.
  using OERangeSelect = adjlist_select_out_edge_range<OutEdgeListS, edge_container, edge_descriptor>;
  auto out_edges(vertex_descriptor u) const {
    return OERangeSelect::create_range(u, get_stored_vertex(u).get_out_edges());
  }

  // NOTE: This WORKS for ALL vertex container types.
  auto edges() const {
    return adjlist_edges_from_out_edges(*this);
  }

  /***************************************************************************************
   * NOTE NOTE NOTE - FUNCTIONS THAT ARE THE DIFFERENT FROM ADJ-LIST-BC - NOTE NOTE NOTE *
   * *************************************************************************************/

  // NOTE: This WORKS for ALL vertex container types.
  auto in_edges(vertex_descriptor v) const {
    using in_edge_range = std::ranges::subrange<const edge_descriptor*>;
    if (get_stored_vertex(v).in_edge == EConfig::null_edge()) {
      return in_edge_range{nullptr, nullptr};
    }
    return in_edge_range{&(get_stored_vertex(v).in_edge), &(get_stored_vertex(v).in_edge) + 1};
  }

  // NOTE: This WORKS for ALL vertex container types.
  vertex_descriptor get_parent(vertex_descriptor v) const {
    return ltree_get_in_edge<DirectedS>(m_vertices, get_stored_vertex(v), v).source;
  }

  // NOTE: This WORKS for ALL vertex container types.
  // NOTE: This WORKS for ALL edge container types.
  std::size_t get_depth(const vertex_descriptor& u) const {
    using TaskType = std::pair<std::size_t, vertex_descriptor>;

    std::size_t max_depth = 0;
    std::stack<TaskType> tasks;
    tasks.push(TaskType(0, u));
    while (!tasks.empty()) {
      TaskType cur = tasks.top();
      tasks.pop();
      ++(cur.first);
      if (cur.first > max_depth) {
        max_depth = cur.first;
      }
      vertex_stored_type& vp = get_stored_vertex(cur.second);
      for (auto e : vp.get_out_edges()) {
        tasks.push(TaskType(cur.first, container_detail::get_value(e).first));
      }
    }
    return max_depth;
  }

  // NOTE: This WORKS for ALL vertex container types.
  // NOTE: This WORKS for ALL edge container types.
  std::pair<edge_descriptor, bool> get_edge(vertex_descriptor u, vertex_descriptor v) const {
    using RawEDesc = typename edge_descriptor::edge_id_type;

    std::pair<RawEDesc, bool> raw_result = adjlist_find_edge_to(get_stored_vertex(u).get_out_edges(), v);

    if (raw_result.second) {
      return {edge_descriptor(u, raw_result.first), true};
    }
    return {edge_descriptor(), false};
  }

  // NOTE: This WORKS for ALL vertex container types.
  // NOTE: This WORKS for ALL edge container types.
  template <typename... VP>
  vertex_descriptor add_new_vertex(VP&&... vp) {
    return adjlist_add_vertex(m_vertices, std::forward<VP>(vp)...);
  }
  template <typename... VP>
  void add_root_vertex(VP&&... vp) {
    m_root = add_new_vertex(std::forward<VP>(vp)...);
  }

  // NOTE: this operation does not invalidate anything.
  // NOTE: This WORKS for ALL vertex container types.
  template <typename EP>
  std::pair<edge_descriptor, bool> add_new_edge(vertex_descriptor u, vertex_descriptor v, EP&& ep) {
    using RawEDesc = typename edge_descriptor::edge_id_type;
    std::pair<RawEDesc, bool> raw_result =
        adjlist_add_edge(get_stored_vertex(u).get_out_edges(), std::forward<EP>(ep), v);
    if (raw_result.second) {
      ltree_add_in_edge(get_stored_vertex(v), edge_descriptor(u, raw_result.first));
      return {edge_descriptor(u, raw_result.first), true};
    }
    return {edge_descriptor(), false};
  }

  // NOTE: this operation does not invalidate anything.
  // NOTE: This WORKS for ALL vertex container types.
  template <typename VP, typename EP>
  std::tuple<vertex_descriptor, edge_descriptor, bool> add_child(vertex_descriptor v, VP&& vp, EP&& ep) {
    vertex_descriptor new_node = add_new_vertex(std::forward<VP>(vp));
    std::pair<edge_descriptor, bool> new_edge = add_new_edge(v, new_node, std::forward<EP>(ep));
    if (!new_edge.second) {
      adjlist_erase_vertex(m_vertices, new_node);
      new_node = container_detail::null_desc<vertex_descriptor>::value();
    }
    return {new_node, new_edge.first, new_edge.second};
  }

  template <typename Vertex_OIter, typename Edge_OIter>
  std::pair<Vertex_OIter, Edge_OIter> clear_children_impl(vertex_descriptor v, Vertex_OIter vit_out,
                                                          Edge_OIter eit_out) {
    std::vector<vertex_descriptor> death_row;
    std::queue<vertex_descriptor> bft_queue;
    bft_queue.push(v);
    // Put all children on death-row:
    while (!bft_queue.empty()) {
      vertex_stored_type& v_value = get_stored_vertex(bft_queue.front());
      bft_queue.pop();
      for (auto e : v_value.get_out_edges()) {
        if (!container_detail::is_elem_valid(e)) {
          continue;
        }
        death_row.push_back(container_detail::get_value(e).first);
        bft_queue.push(container_detail::get_value(e).first);
        vertex_stored_type& u_value = get_stored_vertex(container_detail::get_value(e).first);
        *(vit_out++) = std::move(u_value.data);
        *(eit_out++) = std::move(container_detail::get_value(e).second);
      }
    }

    // Check if we removed the root:
    if (v == m_root) {  // v is the root vertex.
      clear();
    } else {
      // remove the out-edges.
      get_stored_vertex(v).get_out_edges().clear();
      // Execute the death-row vertices!
      ltree_erase_vertices(m_vertices, death_row);
    }

    return {vit_out, eit_out};
  }

  template <typename OutputIter>
  OutputIter clear_children_impl(vertex_descriptor v, OutputIter it_out) {
    return clear_children_impl(v, it_out, container_detail::ignore_output_iter()).first;
  }

  void clear_children_impl(vertex_descriptor v) {
    clear_children_impl(v, container_detail::ignore_output_iter(), container_detail::ignore_output_iter());
  }

  template <typename Vertex_OIter, typename Edge_OIter>
  std::pair<Vertex_OIter, Edge_OIter> remove_branch_impl(vertex_descriptor v, Vertex_OIter vit_out,
                                                         Edge_OIter eit_out) {
    edge_descriptor p_edge = ltree_get_in_edge<DirectedS>(m_vertices, get_stored_vertex(v), v);

    if (p_edge.source == VConfig::null_vertex()) {
      *(eit_out++) = edge_value_type();
    } else {
      *(eit_out++) = std::move(get_stored_edge_property(p_edge));
    }

    std::vector<vertex_descriptor> death_row;
    death_row.push_back(v);
    std::queue<vertex_descriptor> bft_queue;
    bft_queue.push(v);
    // Put all children on death-row:
    while (!bft_queue.empty()) {
      vertex_stored_type& v_value = get_stored_vertex(bft_queue.front());
      *(vit_out++) = std::move(v_value.data);
      bft_queue.pop();
      for (auto& e : v_value.get_out_edges()) {
        if (!container_detail::is_elem_valid(e)) {
          continue;
        }
        death_row.push_back(container_detail::get_value(e).first);
        bft_queue.push(container_detail::get_value(e).first);
        *(eit_out++) = std::move(container_detail::get_value(e).second);
      }
    }

    // Check if we removed the root:
    if (p_edge.source == VConfig::null_vertex()) {
      // v must be the root vertex.
      clear();
      return {vit_out, eit_out};
    }  // remove the edge.
    ltree_erase_edge(get_stored_vertex(p_edge.source).get_out_edges(), p_edge, m_vertices, p_edge.source);

    // Execute the death-row vertices!
    ltree_erase_vertices(m_vertices, death_row);

    return {vit_out, eit_out};
  }

  template <typename OutputIter>
  OutputIter remove_branch_impl(vertex_descriptor v, OutputIter it_out) {
    return remove_branch_impl(v, it_out, container_detail::ignore_output_iter()).first;
  }

  void remove_branch_impl(vertex_descriptor v) {
    remove_branch_impl(v, container_detail::ignore_output_iter(), container_detail::ignore_output_iter());
  }
};

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
void swap(ltree_vertex_container<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& lhs,
          ltree_vertex_container<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& rhs) {
  lhs.swap(rhs);
}

}  // namespace bagl::adjlist_detail

#endif  // BAGL_BAGL_LTREE_CONTAINERS_H_
