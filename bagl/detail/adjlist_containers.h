// Copyright 2013-2024 Mikael Persson.

#ifndef BAGL_BAGL_DETAIL_ADJLIST_CONTAINERS_H_
#define BAGL_BAGL_DETAIL_ADJLIST_CONTAINERS_H_

#include <cstddef>
#include <iterator>
#include <ranges>
#include <type_traits>
#include <utility>

#include "bagl/detail/adjlist_ranges.h"
#include "bagl/detail/container_generators.h"
#include "bagl/edges_from_out_edges.h"
#include "bagl/graph_selectors.h"  // for directed_s, undirected_s, bidirectional_s

namespace bagl::adjlist_detail {

/*************************************************************************
 *        value-types for vertices and edges
 * **********************************************************************/

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
struct adjlist_vertex_stored_type;  // forward declaration.

// this is fine because boost-containers allow incomplete types:
template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
struct adjlist_vertex_config {
  using stored_type =
      adjlist_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>;
  using container = typename container_gen<VertexListS, stored_type>::type;
  using value_type = typename container::value_type;
  using descriptor = typename container_detail::select_descriptor<container>::type;

  static descriptor null_vertex() { return container_detail::null_desc<descriptor>::value(); }
};

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
struct adjlist_edge_config {
  using vertex_descriptor = typename adjlist_vertex_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties,
                                                           EdgeProperties>::descriptor;

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
struct adjlist_vertex_stored_type {
  using self = adjlist_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>;
  using directed_tag = DirectedS;
  using Config = adjlist_edge_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>;
  using edge_container = typename Config::container;
  using edge_descriptor = typename Config::descriptor;
  using in_edge_container = std::vector<edge_descriptor>;
  using in_edge_iterator = typename in_edge_container::iterator;

  VertexProperties data;
  edge_container out_edges;
  in_edge_container in_edges;

  edge_container& get_out_edges() { return out_edges; }
  const edge_container& get_out_edges() const { return out_edges; }

  adjlist_vertex_stored_type() : data(), out_edges(), in_edges() {}
  template <typename VProp>
  requires std::constructible_from<VertexProperties, VProp&&>
  explicit adjlist_vertex_stored_type(VProp&& aData) : data(std::forward<VProp>(aData)), out_edges(), in_edges() {}

  adjlist_vertex_stored_type(const adjlist_vertex_stored_type&) = delete;
  adjlist_vertex_stored_type& operator=(const adjlist_vertex_stored_type&) = delete;
  adjlist_vertex_stored_type(adjlist_vertex_stored_type&&) noexcept = default;
  adjlist_vertex_stored_type& operator=(adjlist_vertex_stored_type&&) noexcept = default;
  ~adjlist_vertex_stored_type() = default;
};

template <typename VertexListS, typename OutEdgeListS, typename VertexProperties, typename EdgeProperties>
struct adjlist_vertex_stored_type<VertexListS, OutEdgeListS, directed_s, VertexProperties, EdgeProperties> {
  using self = adjlist_vertex_stored_type<VertexListS, OutEdgeListS, directed_s, VertexProperties, EdgeProperties>;
  using directed_tag = directed_s;
  using Config = adjlist_edge_config<VertexListS, OutEdgeListS, directed_s, VertexProperties, EdgeProperties>;
  using edge_container = typename Config::container;
  using edge_descriptor = typename Config::descriptor;
  using in_edge_iterator = int*;

  VertexProperties data;
  edge_container out_edges;

  edge_container& get_out_edges() { return out_edges; }
  const edge_container& get_out_edges() const { return out_edges; }

  adjlist_vertex_stored_type() : data(), out_edges() {}
  template <typename VProp>
  requires std::constructible_from<VertexProperties, VProp&&>
  explicit adjlist_vertex_stored_type(VProp&& aData) : data(std::forward<VProp>(aData)), out_edges() {}

  adjlist_vertex_stored_type(const adjlist_vertex_stored_type&) = delete;
  adjlist_vertex_stored_type& operator=(const adjlist_vertex_stored_type&) = delete;
  adjlist_vertex_stored_type(adjlist_vertex_stored_type&&) noexcept = default;
  adjlist_vertex_stored_type& operator=(adjlist_vertex_stored_type&&) noexcept = default;
  ~adjlist_vertex_stored_type() = default;
};

/*************************************************************************
 *        functions for erasing an edge
 * **********************************************************************/
// NOTE: Time complexities:
//               vec_s      pool_s      list_s
// directed_s:    O(1)        O(1)        O(1)
// bidir:       O(E/V)      O(E/V)      O(E/V)   (because of in-edge erasure (vector storage))

template <typename VertexListS, typename OutEdgeListS, typename VertexProperties, typename EdgeProperties,
          typename EdgeDesc>
void adjlist_erase_in_edge(
    adjlist_vertex_stored_type<VertexListS, OutEdgeListS, directed_s, VertexProperties, EdgeProperties>& vp,
    EdgeDesc e) {}

// for vector of edge-desc (in-edges)
template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties, typename EdgeDesc>
void adjlist_erase_in_edge(
    adjlist_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& vp,
    EdgeDesc e) {
  vp.in_edges.erase(std::find(vp.in_edges.begin(), vp.in_edges.end(), e));
}

// for OutEdgeListS = list_s, set_s, ...
template <typename Container, typename EdgeDesc, typename VertexCont, typename VertexDesc>
void adjlist_erase_edge(Container& cont, EdgeDesc e, VertexCont& /*unused*/, VertexDesc /*unused*/) {
  cont.erase(e.edge_id);
}

// for OutEdgeListS = unordered_set_s
template <typename Key, typename Value, typename EdgeDesc, typename VertexCont, typename VertexDesc>
void adjlist_erase_edge(std::unordered_map<Key, Value>& cont, EdgeDesc e, VertexCont& /*unused*/,
                        VertexDesc /*unused*/) {
  auto [ei, ei_end] = cont.equal_range(e.edge_id.first);
  for (; ei != ei_end; ++ei) {
    if (&ei->second == e.edge_id.second) {
      cont.erase(ei);
      break;
    }
  }
}

// for OutEdgeListS = unordered_set_s
template <typename Key, typename Value, typename EdgeDesc, typename VertexCont, typename VertexDesc>
void adjlist_erase_edge(std::unordered_multimap<Key, Value>& cont, EdgeDesc e, VertexCont& /*unused*/,
                        VertexDesc /*unused*/) {
  auto [ei, ei_end] = cont.equal_range(e.edge_id.first);
  for (; ei != ei_end; ++ei) {
    if (&ei->second == e.edge_id.second) {
      cont.erase(ei);
      break;
    }
  }
}

// for OutEdgeListS = vec_s
template <typename ValueType, typename EdgeDesc, typename VertexCont, typename VertexDesc>
void adjlist_erase_edge(std::vector<ValueType>& cont, EdgeDesc e, VertexCont& vertex_cont, VertexDesc v) {
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

// for OutEdgeListS = pool_s
template <typename ValueType, typename EdgeDesc, typename VertexCont, typename VertexDesc>
void adjlist_erase_edge(container_detail::pooled_vector<ValueType>& cont, EdgeDesc e, VertexCont& /*unused*/,
                        VertexDesc /*unused*/) {
  cont.m_data[e.edge_id] = cont.m_first_hole;
  cont.m_first_hole = container_detail::hole_desc(e.edge_id);
  --(cont.m_num_elements);
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
void adjlist_add_in_edge(
    adjlist_vertex_stored_type<VertexListS, OutEdgeListS, directed_s, VertexProperties, EdgeProperties>& vp,
    EdgeDesc e) {}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties, typename EdgeDesc>
void adjlist_add_in_edge(
    adjlist_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& vp,
    EdgeDesc e) {
  vp.in_edges.push_back(e);
}

// for OutEdgeListS = list_s, multiset_s, ...
template <typename Container, typename EdgeProperties, typename VertexDesc>
auto adjlist_add_edge(Container& cont, EdgeProperties&& ep, VertexDesc v) {
  using ValueType = typename Container::value_type;
  return std::pair(cont.insert(cont.end(), ValueType(v, std::forward<EdgeProperties>(ep))), true);
}

// for OutEdgeListS = set_s
template <typename Key, typename ValueType, typename EdgeProperties, typename VertexDesc>
auto adjlist_add_edge(std::map<Key, ValueType>& cont, EdgeProperties&& ep, VertexDesc v) {
  return cont.emplace(v, std::forward<EdgeProperties>(ep));
}

// for OutEdgeListS = multiset_s
template <typename Key, typename ValueType, typename EdgeProperties, typename VertexDesc>
auto adjlist_add_edge(std::multimap<Key, ValueType>& cont, EdgeProperties&& ep, VertexDesc v) {
  return std::pair(cont.emplace(v, std::forward<EdgeProperties>(ep)), true);
}

// for OutEdgeListS = unordered_set_s
template <typename Key, typename ValueType, typename EdgeProperties, typename VertexDesc>
std::pair<std::pair<VertexDesc, ValueType*>, bool> adjlist_add_edge(std::unordered_map<Key, ValueType>& cont,
                                                                    EdgeProperties&& ep, VertexDesc v) {
  auto [it, inserted] = cont.emplace(v, std::forward<EdgeProperties>(ep));
  if (!inserted) {
    return {{v, nullptr}, false};
  }
  return {{v, &it->second}, true};
}

// for OutEdgeListS = unordered_multiset_s
template <typename Key, typename ValueType, typename EdgeProperties, typename VertexDesc>
std::pair<std::pair<VertexDesc, ValueType*>, bool> adjlist_add_edge(std::unordered_multimap<Key, ValueType>& cont,
                                                                    EdgeProperties&& ep, VertexDesc v) {
  auto it = cont.emplace(v, std::forward<EdgeProperties>(ep));
  return {{v, &it->second}, true};
}

// for OutEdgeListS = vec_s
template <typename ValueType, typename EdgeProperties, typename VertexDesc>
std::pair<std::size_t, bool> adjlist_add_edge(std::vector<ValueType>& cont, EdgeProperties&& ep, VertexDesc v) {
  cont.emplace_back(v, std::forward<EdgeProperties>(ep));
  return {cont.size() - 1, true};
}

// for OutEdgeListS = pool_s
template <typename ValueType, typename EdgeProperties, typename VertexDesc>
std::pair<std::size_t, bool> adjlist_add_edge(container_detail::pooled_vector<ValueType>& cont, EdgeProperties&& ep,
                                              VertexDesc v) {
  if (cont.m_first_hole == container_detail::hole_desc()) {
    auto it = cont.m_data.insert(cont.m_data.end(), ValueType(v, std::forward<EdgeProperties>(ep)));
    ++(cont.m_num_elements);
    return {it - cont.m_data.begin(), true};
  }
  auto it = cont.m_data.begin() + cont.m_first_hole.value;
  cont.m_first_hole = std::get<container_detail::hole_desc>(cont.m_data[cont.m_first_hole.value]);
  *it = ValueType(v, std::forward<EdgeProperties>(ep));
  ++(cont.m_num_elements);
  return {it - cont.m_data.begin(), true};
}

/*************************************************************************
 *        helper functions for finding an edgee between two vertices
 * **********************************************************************/
// NOTE: Time complexities:
//               vec_s      pool_s      list_s    (multi)set_s      unordered_(multi)set_s
// any-dir.:    O(E/V)      O(E/V)      O(E/V)    O(log(E/V))       O(1)

// for OutEdgeListS = list_s
template <typename ValueType, typename VertexDesc>
auto adjlist_find_edge_to(std::list<ValueType>& cont, VertexDesc v) {
  auto it = cont.begin();
  for (; it != cont.end(); ++it) {
    if (it->first == v) {
      return std::pair(it, true);
    }
  }
  return std::pair(it, false);
}

// for OutEdgeListS = set_s, multiset_s
template <typename Key, typename Value, typename VertexDesc>
auto adjlist_find_edge_to(std::map<Key, Value>& cont, VertexDesc v) {
  auto it = cont.find(v);
  return std::pair(it, (it != cont.end()));
}
template <typename Key, typename Value, typename VertexDesc>
auto adjlist_find_edge_to(std::multimap<Key, Value>& cont, VertexDesc v) {
  auto it = cont.find(v);
  return std::pair(it, (it != cont.end()));
}

// for OutEdgeListS = unordered_set_s, unordered_multiset_s
template <typename Key, typename Value, typename VertexDesc>
std::pair<std::pair<VertexDesc, Value*>, bool> adjlist_find_edge_to(std::unordered_map<Key, Value>& cont,
                                                                    VertexDesc v) {
  auto it = cont.find(v);
  if (it == cont.end()) {
    return {{v, nullptr}, false};
  }
  return {{v, &it->second}, true};
}
template <typename Key, typename Value, typename VertexDesc>
std::pair<std::pair<VertexDesc, Value*>, bool> adjlist_find_edge_to(std::unordered_multimap<Key, Value>& cont,
                                                                    VertexDesc v) {
  auto it = cont.find(v);
  if (it == cont.end()) {
    return {{v, nullptr}, false};
  }
  return {{v, &it->second}, true};
}

// for OutEdgeListS = vec_s
template <typename ValueType, typename VertexDesc>
std::pair<std::size_t, bool> adjlist_find_edge_to(std::vector<ValueType>& cont, VertexDesc v) {
  for (auto it = cont.begin(); it != cont.end(); ++it) {
    if (it->first == v) {
      return {it - cont.begin(), true};
    }
  }
  return {cont.size(), false};
}

// for OutEdgeListS = pool_s
template <typename ValueType, typename VertexDesc>
std::pair<std::size_t, bool> adjlist_find_edge_to(container_detail::pooled_vector<ValueType>& cont, VertexDesc v) {
  for (auto it = cont.m_data.begin(); it != cont.m_data.end(); ++it) {
    if ((it->index() == 0) && (get<ValueType>(*it).first == v)) {
      return {it - cont.m_data.begin(), true};
    }
  }
  return {cont.m_data.size(), false};
}

/*************************************************************************
 *        helper functions for erasing in-edges / out-edges of a vertex
 * **********************************************************************/
// NOTE: Time complexities:
//               vec_s      pool_s      list_s    (multi)set_s      unordered_(multi)set_s
// directed_s:  O(E/V)      O(E/V)      O(E/V)    O(log(E/V))       O(1)
// bidir:     O((E/V)^2)    O(E/V)      O(E/V)    O(log(E/V))       O(1)

// for OutEdgeListS = list_s, set_s, multiset_s, ...
template <typename Container, typename VertexCont, typename VertexDesc, typename IterPred>
void adjlist_erase_edges_if(Container& cont, VertexCont& /*vcont*/, VertexDesc /*u*/, std::size_t& e_count,
                            IterPred pred) {
  for (auto it = cont.begin(); it != cont.end();) {
    if (pred(it)) {
      it = cont.erase(it);
      --e_count;
    } else {
      ++it;
    }
  }
}

// for OutEdgeListS = list_s
// List container must scan for target.
template <typename ValueType, typename VertexCont, typename VertexDesc>
void adjlist_erase_edges_to(std::list<ValueType>& cont, VertexCont& vcont, VertexDesc u, VertexDesc v,
                            std::size_t& e_count) {
  adjlist_erase_edges_if(cont, vcont, u, e_count, [v](const auto& it) { return it->first == v; });
}

// for OutEdgeListS = set_s, multiset_s, ...
// Associative container can directly look up target.
template <typename Container, typename VertexCont, typename VertexDesc>
void adjlist_erase_edges_to(Container& cont, VertexCont& /*vcont*/, VertexDesc /*u*/, VertexDesc v,
                            std::size_t& e_count) {
  e_count -= cont.erase(v);
}

template <typename VertexListS, typename OutEdgeListS, typename VertexProperties, typename EdgeProperties,
          typename VertexDesc>
void adjlist_update_in_edge_id(
    adjlist_vertex_stored_type<VertexListS, OutEdgeListS, directed_s, VertexProperties, EdgeProperties>& vp,
    VertexDesc v, std::size_t old_id, std::size_t new_id) {}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties, typename VertexDesc>
void adjlist_update_in_edge_id(
    adjlist_vertex_stored_type<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& vp,
    VertexDesc v, std::size_t old_id, std::size_t new_id) {
  for (auto& e : vp.in_edges) {
    if ((e.source == v) && (e.edge_id == old_id)) {
      e.edge_id = new_id;
      break;
    }
  }
}

// for OutEdgeListS = vec_s
template <typename ValueType, typename VertexCont, typename VertexDesc, typename IterPred>
void adjlist_erase_edges_if(std::vector<ValueType>& cont, VertexCont& vcont, VertexDesc u, std::size_t& e_count,
                            IterPred pred) {
  using std::swap;

  auto it_last = cont.end();
  for (auto it = cont.begin(); it != it_last;) {
    if (pred(it)) {
      --it_last;
      if (it != it_last) {
        swap(*it, *it_last);
        // If this graph has in-edge references, then they must be updated.
        adjlist_update_in_edge_id(container_detail::desc_to_value(vcont, it->first), u, it_last - cont.begin(),
                                  it - cont.begin());
      }
      --e_count;
    } else {
      ++it;
    }
  }
  cont.erase(it_last, cont.end());
}

template <typename ValueType, typename VertexCont, typename VertexDesc>
void adjlist_erase_edges_to(std::vector<ValueType>& cont, VertexCont& vcont, VertexDesc u, VertexDesc v,
                            std::size_t& e_count) {
  adjlist_erase_edges_if(cont, vcont, u, e_count,
                         [v](const auto& it) { return container_detail::get_value(*it).first == v; });
}

// for OutEdgeListS = pool_s
template <typename ValueType, typename VertexCont, typename VertexDesc, typename IterPred>
void adjlist_erase_edges_if(container_detail::pooled_vector<ValueType>& cont, VertexCont& /*vcont*/, VertexDesc /*u*/,
                            std::size_t& e_count, IterPred pred) {
  for (auto it = cont.m_data.begin(); it != cont.m_data.end(); ++it) {
    if ((it->index() == 0) && pred(it)) {
      *it = cont.m_first_hole;
      cont.m_first_hole = container_detail::hole_desc(it - cont.m_data.begin());
      --(cont.m_num_elements);
      --e_count;
    }
  }
}

template <typename ValueType, typename VertexCont, typename VertexDesc>
void adjlist_erase_edges_to(container_detail::pooled_vector<ValueType>& cont, VertexCont& vcont, VertexDesc u,
                            VertexDesc v, std::size_t& e_count) {
  adjlist_erase_edges_if(cont, vcont, u, e_count,
                         [v](const auto& it) { return container_detail::get_value(*it).first == v; });
}

/*************************************************************************
 *        functions for clearing the edges of a vertex
 * **********************************************************************/
// NOTE: The function 'adjlist_clear_vertex' works for all graph types.

// NOTE: Time complexities:
//                 vec_s         pool_s         list_s
// directed_s:   O(V*(E/V))     O(V*(E/V))     O(V*(E/V))
// bidir:        O((E/V)^2)     O((E/V)^2)     O((E/V)^2)

template <typename DirectedS, typename VertexCont, typename VertexValue, typename VertexDesc>
std::enable_if_t<std::is_same_v<DirectedS, directed_s>> adjlist_clear_vertex(VertexCont& vcont, VertexValue& vp,
                                                                             VertexDesc v, std::size_t& e_count) {
  // first, just clear the out-going edges. No need to synchronize in-edges (there are none).
  e_count -= vp.get_out_edges().size();
  vp.get_out_edges().clear();

  // now, the stupid part... we have to traverse all other vertices (and their edges)
  // to look for in-edges that lead to "v", and erase them.
  for (auto ui = vcont.begin(); ui != vcont.end(); ++ui) {
    auto u_desc = container_detail::iterator_to_desc(vcont, ui);
    if ((u_desc == v) || !container_detail::is_elem_valid(*ui)) {
      continue;
    }
    VertexValue& up = container_detail::get_value(*ui);
    adjlist_erase_edges_to(up.get_out_edges(), vcont, u_desc, v, e_count);
  }
}

template <typename DirectedS, typename VertexCont, typename VertexValue, typename VertexDesc>
std::enable_if_t<!std::is_same_v<DirectedS, directed_s>> adjlist_clear_vertex(VertexCont& vcont, VertexValue& vp,
                                                                              VertexDesc v, std::size_t& e_count) {
  // first, remove the in-edge references from the adjacent vertices of v.
  for (auto ei = vp.get_out_edges().begin(), ei_end = vp.get_out_edges().end(); ei != ei_end; ++ei) {
    if (!container_detail::is_elem_valid(*ei)) {
      continue;
    }
    auto& [e_target, e_prop] = container_detail::get_value(*ei);
    VertexValue& wp = container_detail::desc_to_value(vcont, e_target);
    for (auto iei = wp.in_edges.begin(); iei != wp.in_edges.end(); ++iei) {
      if ((iei->source == v) && (container_detail::iterator_to_desc(vp.get_out_edges(), ei) == iei->edge_id)) {
        wp.in_edges.erase(iei);
        break;
      }
    }
  }

  // then, clear the out-going edges.
  e_count -= vp.get_out_edges().size();
  vp.get_out_edges().clear();

  // finally, remove the required out-edges of the "parent" vertices of v.
  for (const auto& ie : vp.in_edges) {
    VertexDesc u = ie.source;
    VertexValue& up = container_detail::desc_to_value(vcont, ie.source);
    adjlist_erase_edges_to(up.get_out_edges(), vcont, u, v, e_count);
  }
  vp.in_edges.clear();
}

/*************************************************************************
 *        functions for erasing a vertex (including updating edges of surrounding vertices
 * **********************************************************************/
// NOTE: The function 'adjlist_erase_vertex' works for all graph types.
// NOTE: Time complexities:
//               vec_s      pool_s      list_s
// directed_s:     O(E)        O(1)        O(1)
// bidir:      O((E/V)^2)      O(1)        O(1)

template <typename Container, typename VertexDesc>
void adjlist_erase_vertex(Container& cont, VertexDesc v) {
  cont.erase(v);
}

template <typename Container>
void adjlist_update_out_edges_impl(Container& cont, std::size_t old_v_id, std::size_t new_v_id) {
  for (auto& e : cont) {
    if (container_detail::is_elem_valid(e) && (container_detail::get_value(e).first == old_v_id)) {
      container_detail::get_value(e).first = new_v_id;
    }
  }
}

template <typename Container>
void adjlist_update_assoc_out_edges_impl(Container& cont, std::size_t old_v_id, std::size_t new_v_id) {
  auto [eq_first, eq_last] = cont.equal_range(old_v_id);
  for (auto eq_it = eq_first; eq_it != eq_last;) {
    // Get the pointer to old prop, that won't get invalidated.
    auto* old_e_prop_ptr = &eq_it->second;
    cont.emplace(new_v_id, std::move(eq_it->second));
    // After insert, iterators may be garbage (for unordered map).
    std::tie(eq_first, eq_last) = cont.equal_range(old_v_id);
    eq_it = std::find_if(eq_first, eq_last, [old_e_prop_ptr](const auto& p) { return (&p.second == old_e_prop_ptr); });
    // The old edge must still be there.
    assert(eq_it != eq_last);
    eq_it = cont.erase(eq_it);
  }
}

template <typename Whole>
Whole adjlist_update_assoc_out_edges_first_or_whole(Whole w) {
  return w;
}
template <typename T, typename U>
T adjlist_update_assoc_out_edges_first_or_whole(std::pair<T, U> p) {
  return p.first;
}

template <typename Container, typename EDesc>
void adjlist_update_assoc_out_edges_impl(Container& cont, std::size_t old_v_id, std::size_t new_v_id, EDesc& ie_desc) {
  auto [eq_first, eq_last] = cont.equal_range(old_v_id);
  for (auto eq_it = eq_first; eq_it != eq_last; ++eq_it) {
    if (container_detail::iterator_to_desc(cont, eq_it) != ie_desc) {
      continue;
    }
    // Get the pointer to old prop, that won't get invalidated.
    auto* old_e_prop_ptr = &eq_it->second;
    auto new_eq_it = adjlist_update_assoc_out_edges_first_or_whole(cont.emplace(new_v_id, std::move(eq_it->second)));
    ie_desc = container_detail::iterator_to_desc(cont, new_eq_it);
    // After insert, iterators may be garbage (for unordered map).
    std::tie(eq_first, eq_last) = cont.equal_range(old_v_id);
    eq_it = std::find_if(eq_first, eq_last, [old_e_prop_ptr](const auto& p) { return (&p.second == old_e_prop_ptr); });
    // The old edge must still be there.
    assert(eq_it != eq_last);
    cont.erase(eq_it);
    // We break because this function is called for each in-edge descriptor.
    break;
  }
}

template <typename ValueType>
void adjlist_update_out_edges_impl(std::map<std::size_t, ValueType>& cont, std::size_t old_v_id, std::size_t new_v_id) {
  adjlist_update_assoc_out_edges_impl(cont, old_v_id, new_v_id);
}
template <typename ValueType>
void adjlist_update_out_edges_impl(std::multimap<std::size_t, ValueType>& cont, std::size_t old_v_id,
                                   std::size_t new_v_id) {
  adjlist_update_assoc_out_edges_impl(cont, old_v_id, new_v_id);
}
template <typename ValueType>
void adjlist_update_out_edges_impl(std::unordered_map<std::size_t, ValueType>& cont, std::size_t old_v_id,
                                   std::size_t new_v_id) {
  adjlist_update_assoc_out_edges_impl(cont, old_v_id, new_v_id);
}
template <typename ValueType>
void adjlist_update_out_edges_impl(std::unordered_multimap<std::size_t, ValueType>& cont, std::size_t old_v_id,
                                   std::size_t new_v_id) {
  adjlist_update_assoc_out_edges_impl(cont, old_v_id, new_v_id);
}

template <typename Container, typename EDesc>
void adjlist_update_out_edges_impl(Container& cont, std::size_t old_v_id, std::size_t new_v_id, EDesc e_desc) {
  container_detail::desc_to_value(cont, e_desc).first = new_v_id;
}
template <typename ValueType, typename EDesc>
void adjlist_update_out_edges_impl(std::map<std::size_t, ValueType>& cont, std::size_t old_v_id, std::size_t new_v_id,
                                   EDesc& ie_desc) {
  adjlist_update_assoc_out_edges_impl(cont, old_v_id, new_v_id, ie_desc);
}
template <typename ValueType, typename EDesc>
void adjlist_update_out_edges_impl(std::multimap<std::size_t, ValueType>& cont, std::size_t old_v_id,
                                   std::size_t new_v_id, EDesc& ie_desc) {
  adjlist_update_assoc_out_edges_impl(cont, old_v_id, new_v_id, ie_desc);
}
template <typename ValueType, typename EDesc>
void adjlist_update_out_edges_impl(std::unordered_map<std::size_t, ValueType>& cont, std::size_t old_v_id,
                                   std::size_t new_v_id, EDesc& ie_desc) {
  adjlist_update_assoc_out_edges_impl(cont, old_v_id, new_v_id, ie_desc);
}
template <typename ValueType, typename EDesc>
void adjlist_update_out_edges_impl(std::unordered_multimap<std::size_t, ValueType>& cont, std::size_t old_v_id,
                                   std::size_t new_v_id, EDesc& ie_desc) {
  adjlist_update_assoc_out_edges_impl(cont, old_v_id, new_v_id, ie_desc);
}

// O(E)
template <typename DirectedS, typename ValueType>
std::enable_if_t<std::is_same_v<DirectedS, directed_s>> adjlist_update_out_edges(std::vector<ValueType>& cont,
                                                                                 std::size_t old_v_id,
                                                                                 std::size_t new_v_id) {
  for (auto& v : cont) {
    adjlist_update_out_edges_impl(v.get_out_edges(), old_v_id, new_v_id);
  }
}

// O((E/V)^2)
template <typename DirectedS, typename ValueType>
std::enable_if_t<!std::is_same_v<DirectedS, directed_s>> adjlist_update_out_edges(std::vector<ValueType>& cont,
                                                                                  std::size_t old_v_id,
                                                                                  std::size_t new_v_id) {
  // first, update in-edge vertices
  for (auto iei = cont[new_v_id].in_edges.begin(); iei != cont[new_v_id].in_edges.end(); ++iei) {
    ValueType& up = cont[iei->source];
    adjlist_update_out_edges_impl(up.get_out_edges(), old_v_id, new_v_id, iei->edge_id);
  }

  // second, update out-edge vertices
  for (auto ei = cont[new_v_id].get_out_edges().begin(); ei != cont[new_v_id].get_out_edges().end(); ++ei) {
    if (!container_detail::is_elem_valid(*ei)) {
      continue;
    }
    auto& [e_target, e_prop] = container_detail::get_value(*ei);
    ValueType& wp = cont[e_target];
    for (auto iei = wp.in_edges.begin(); iei != wp.in_edges.end(); ++iei) {
      if ((iei->source == old_v_id) &&
          (iei->edge_id == container_detail::iterator_to_desc(cont[new_v_id].get_out_edges(), ei))) {
        iei->source = new_v_id;
        break;
      }
    }
  }
}

template <typename ValueType>
void adjlist_erase_vertex(std::vector<ValueType>& cont, std::size_t v) {
  using DirectedS = typename ValueType::directed_tag;
  using std::swap;

  auto it = cont.begin() + v;
  auto it_last = cont.end();
  --it_last;
  if (it == it_last) {
    cont.erase(it_last);
    return;
  }
  swap(*it, *it_last);
  std::size_t old_id = it_last - cont.begin();
  std::size_t new_id = it - cont.begin();
  cont.erase(it_last);
  adjlist_update_out_edges<DirectedS>(cont, old_id, new_id);
}

template <typename ValueType>
void adjlist_erase_vertex(container_detail::pooled_vector<ValueType>& cont, std::size_t v) {
  // the first_hole will become v and v will be referring to first_hole (i.e., like a simple linked-list):
  cont.m_data[v] = cont.m_first_hole;
  cont.m_first_hole = container_detail::hole_desc(v);
  --(cont.m_num_elements);
}

/*************************************************************************
 *        functions for adding an edge
 * **********************************************************************/
// NOTE: Time complexities:
//               vec_s      pool_s      list_s
// directed_s:    O(1)        O(1)        O(1)
// bidir:         O(1)        O(1)        O(1)

template <typename Container, typename VertexProperties>
auto adjlist_add_vertex(Container& cont, VertexProperties&& vp) {
  using ValueType = typename Container::value_type;
  return cont.insert(cont.end(), ValueType{std::forward<VertexProperties>(vp)});
}

template <typename ValueType, typename VertexProperties>
std::size_t adjlist_add_vertex(std::vector<ValueType>& cont, VertexProperties&& vp) {
  auto it = cont.insert(cont.end(), ValueType{std::forward<VertexProperties>(vp)});
  return it - cont.begin();
}

template <typename ValueType, typename VertexProperties>
std::size_t adjlist_add_vertex(container_detail::pooled_vector<ValueType>& cont, VertexProperties&& vp) {
  if (cont.m_first_hole == container_detail::hole_desc()) {
    auto it = cont.m_data.insert(cont.m_data.end(), ValueType{std::forward<VertexProperties>(vp)});
    ++(cont.m_num_elements);
    return it - cont.m_data.begin();
  }
  auto it = cont.m_data.begin() + cont.m_first_hole.value;
  cont.m_first_hole = std::get<container_detail::hole_desc>(cont.m_data[cont.m_first_hole.value]);
  *it = ValueType(std::forward<VertexProperties>(vp));
  ++(cont.m_num_elements);
  return it - cont.m_data.begin();
}

template <typename Container, typename VertexProperties>
auto adjlist_add_vertex(Container* cont, VertexProperties&& vp) {
  return adjlist_add_vertex(*cont, std::forward<VertexProperties>(vp));
}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
struct adjlist_vertex_container {
  using self = adjlist_vertex_container<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>;

  using VConfig = adjlist_vertex_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>;
  using vertex_container = typename VConfig::container;
  using vertices_size_type = typename vertex_container::size_type;
  using vertex_descriptor = typename VConfig::descriptor;
  using vertex_stored_type = typename VConfig::stored_type;
  using vertex_value_type = typename VConfig::value_type;

  using EConfig = adjlist_edge_config<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>;
  using edge_container = typename EConfig::container;
  using edges_size_type = typename edge_container::size_type;
  using edge_descriptor = typename EConfig::descriptor;
  using edge_stored_type = typename EConfig::stored_type;
  using edge_value_type = typename EConfig::value_type;

  static vertex_descriptor null_vertex() { return VConfig::null_vertex(); }

  mutable vertex_container m_vertices;
  std::size_t m_num_edges{0};

  adjlist_vertex_container() : m_vertices() {}
  explicit adjlist_vertex_container(std::size_t num_vertices) : m_vertices() {
    for (std::size_t i = 0; i < num_vertices; ++i) {
      this->add_vertex(VertexProperties{});
    }
  }

  ~adjlist_vertex_container() { clear(); }

  adjlist_vertex_container(const adjlist_vertex_container&) = delete;
  adjlist_vertex_container& operator=(const adjlist_vertex_container&) = delete;

  void swap(adjlist_vertex_container& rhs) {
    using std::swap;
    m_vertices.swap(rhs.m_vertices);
    // swap(m_vertices, rhs.m_vertices);
    swap(m_num_edges, rhs.m_num_edges);
  }

  adjlist_vertex_container(adjlist_vertex_container&& rhs) noexcept : m_vertices() { swap(rhs); }
  adjlist_vertex_container& operator=(adjlist_vertex_container&& rhs) noexcept {
    swap(rhs);
    return *this;
  }

  std::size_t size() const { return m_vertices.size(); }
  std::size_t capacity() const { return m_vertices.capacity(); }

  std::size_t num_edges() const { return m_num_edges; }

  vertex_stored_type& get_stored_vertex(vertex_descriptor v) const {
    return container_detail::desc_to_value(m_vertices, v);
  }

  auto& get_stored_edge_property(const edge_descriptor& e) {
    return container_detail::desc_to_value(get_stored_vertex(e.source).get_out_edges(), e.edge_id).second;
  }
  const auto& get_stored_edge_property(const edge_descriptor& e) const {
    return container_detail::desc_to_value(get_stored_vertex(e.source).get_out_edges(), e.edge_id).second;
  }
  vertex_descriptor get_stored_edge_target(const edge_descriptor& e) const {
    return container_detail::desc_to_value(get_stored_vertex(e.source).get_out_edges(), e.edge_id).first;
  }

  std::size_t get_out_degree(vertex_descriptor v) const { return get_stored_vertex(v).get_out_edges().size(); }

  std::size_t get_in_degree(vertex_descriptor v) const { return get_stored_vertex(v).in_edges.size(); }

  // NOTE: this operation does not invalidate anything.
  // NOTE: This WORKS for ALL vertex container types.
  template <typename VProp>
  vertex_descriptor add_vertex(VProp&& vp) { return adjlist_add_vertex(m_vertices, std::forward<VProp>(vp)); }

  // NOTE: this operation only invalidates existing vertex-iterators,
  // and possibly edge-descriptors linked to vertices adjacent to v (if edge-list is vec_s).
  // NOTE: This WORKS for ALL vertex container types.
  void clear_vertex(vertex_descriptor v) {
    adjlist_clear_vertex<DirectedS>(m_vertices, get_stored_vertex(v), v, m_num_edges);
  }

  // NOTE: this operation only invalidates existing vertex-descriptors,
  // and possibly edge-descriptors linked to vertices adjacent to v (if edge-list is vec_s).
  // NOTE: This WORKS for ALL vertex container types.
  void remove_vertex(vertex_descriptor v) {
    clear_vertex(v);
    adjlist_erase_vertex(m_vertices, v);
  }

  // NOTE: this operation does not invalidate anything.
  // NOTE: This WORKS for ALL vertex container types.
  template <typename EProp>
  std::pair<edge_descriptor, bool> add_edge(vertex_descriptor u, vertex_descriptor v, EProp&& ep) {
    using RawEDesc = typename edge_descriptor::edge_id_type;

    std::pair<RawEDesc, bool> raw_result =
        adjlist_add_edge(get_stored_vertex(u).get_out_edges(), std::forward<EProp>(ep), v);

    if (raw_result.second) {
      ++m_num_edges;
      adjlist_add_in_edge(get_stored_vertex(v), edge_descriptor(u, raw_result.first));
      return std::pair<edge_descriptor, bool>(edge_descriptor(u, raw_result.first), true);
    }
    return std::pair<edge_descriptor, bool>(edge_descriptor(), false);
  }

  // NOTE: this operation might invalidate other out-edge iterators/descriptors of the same source vertex.
  // NOTE: This WORKS for ALL vertex container types.
  void remove_edge(const edge_descriptor& e) {
    adjlist_erase_in_edge(get_stored_vertex(get_stored_edge_target(e)), e);
    adjlist_erase_edge(get_stored_vertex(e.source).get_out_edges(), e, m_vertices, e.source);
    --m_num_edges;
  }

  // NOTE: this operation might invalidate other out-edge iterators/descriptors of the same source vertex.
  // NOTE: This WORKS for ALL vertex container types.
  template <typename EdgePred>
  void remove_out_edge_if(vertex_descriptor u, const EdgePred& pred) {
    auto& out_econt = get_stored_vertex(u).get_out_edges();
    adjlist_erase_edges_if(out_econt, m_vertices, u, m_num_edges, [u, &pred, &out_econt](const auto& ei) {
      return pred(edge_descriptor{u, container_detail::iterator_to_desc(out_econt, ei)});
    });
  }

  // NOTE: this operation might invalidate other out-edge iterators/descriptors of the same source vertex.
  // NOTE: This WORKS for ALL vertex container types.
  template <typename EdgePred>
  void remove_edge_if(const EdgePred& pred) {
    for (auto u : VRangeSelect::create_range(m_vertices)) {
      this->remove_out_edge_if(u, pred);
    }
  }

  // NOTE: this operation might invalidate other out-edge iterators/descriptors of the same source vertex.
  // NOTE: This WORKS for ALL vertex container types.
  template <typename EdgePred>
  void remove_in_edge_if(vertex_descriptor v, const EdgePred& pred) {
    static_assert(!std::is_same_v<DirectedS, directed_s>);
    auto& in_econt = get_stored_vertex(v).in_edges;
    for (std::size_t i = 0; i < in_econt.size();) {
      if (pred(in_econt[i])) {
        this->remove_edge(in_econt[i]);
      } else {
        ++i;
      }
    }
  }

  // NOTE: This WORKS for ALL vertex container types.
  void clear() {
    m_vertices.clear();
    m_num_edges = 0;
  }

  // NOTE: This WORKS for ALL vertex container types.
  using VRangeSelect = adjlist_select_vertex_range<VertexListS, vertex_container>;
  auto vertices() const { return VRangeSelect::create_range(m_vertices); }

  // NOTE: This WORKS for ALL vertex container types.
  // NOTE: This WORKS for ALL edge container types.
  using OERangeSelect = adjlist_select_out_edge_range<OutEdgeListS, edge_container, edge_descriptor>;
  auto out_edges(vertex_descriptor u) const {
    if (u == VConfig::null_vertex()) {
      return OERangeSelect::create_empty_range(u);
    }
    return OERangeSelect::create_range(u, get_stored_vertex(u).get_out_edges());
  }

  // NOTE: This WORKS for ALL vertex container types.
  // NOTE: This WORKS for ALL edge container types.
  std::pair<edge_descriptor, bool> get_edge(vertex_descriptor u, vertex_descriptor v) const {
    auto [e_id, e_found] = adjlist_find_edge_to(get_stored_vertex(u).get_out_edges(), v);
    if (e_found) {
      return {edge_descriptor(u, e_id), true};
    }
    return {edge_descriptor(), false};
  }

  // NOTE: This WORKS for ALL vertex container types.
  auto in_edges(vertex_descriptor v) const {
    if constexpr (!std::is_same_v<DirectedS, directed_s>) {
      return std::ranges::subrange(get_stored_vertex(v).in_edges.begin(), get_stored_vertex(v).in_edges.end());
    } else {
      return int{};
    }
  }

  // NOTE: This WORKS for ALL vertex container types.
  auto edges() const { return edges_from_out_edges{*this}; }
};

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
void swap(adjlist_vertex_container<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& lhs,
          adjlist_vertex_container<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& rhs) {
  lhs.swap(rhs);
}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
auto vertices(
    const adjlist_vertex_container<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& g) {
  return g.vertices();
}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
auto out_edges(
    typename adjlist_vertex_container<VertexListS, OutEdgeListS, DirectedS, VertexProperties,
                                      EdgeProperties>::vertex_descriptor u,
    const adjlist_vertex_container<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& g) {
  return g.out_edges(u);
}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
auto in_edges(
    typename adjlist_vertex_container<VertexListS, OutEdgeListS, DirectedS, VertexProperties,
                                      EdgeProperties>::vertex_descriptor u,
    const adjlist_vertex_container<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& g) {
  return g.in_edges(u);
}

template <typename VertexListS, typename OutEdgeListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
auto edges(const adjlist_vertex_container<VertexListS, OutEdgeListS, DirectedS, VertexProperties, EdgeProperties>& g) {
  return g.edges();
}

}  // namespace bagl::adjlist_detail

#endif  // BAGL_BAGL_DETAIL_ADJLIST_CONTAINERS_H_
