// Copyright Michael Drexl 2005, 2006.
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_R_C_SHORTEST_PATHS_H_
#define BAGL_BAGL_R_C_SHORTEST_PATHS_H_

#include <list>
#include <map>
#include <memory>
#include <queue>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"

namespace bagl {

// r_c_shortest_paths_label struct
template <class Graph, class ResourceContainer>
struct r_c_shortest_paths_label {
  explicit r_c_shortest_paths_label(std::size_t n, const ResourceContainer& rc = ResourceContainer(),
                                    std::shared_ptr<r_c_shortest_paths_label<Graph, ResourceContainer>> pl = {},
                                    graph_edge_descriptor_t<Graph> ed = {}, graph_vertex_descriptor_t<Graph> vd = {})
      : num(n), cumulated_resource_consumption(rc), p_pred_label(pl), pred_edge(ed), resident_vertex(vd) {}

  std::size_t num = 0;
  ResourceContainer cumulated_resource_consumption;
  std::shared_ptr<r_c_shortest_paths_label<Graph, ResourceContainer>> p_pred_label;
  graph_edge_descriptor_t<Graph> pred_edge;
  graph_vertex_descriptor_t<Graph> resident_vertex;
  bool b_is_dominated = false;
  bool b_is_processed = false;
};

template <class Graph, class ResourceContainer>
auto operator<=>(const r_c_shortest_paths_label<Graph, ResourceContainer>& l1,
                 const r_c_shortest_paths_label<Graph, ResourceContainer>& l2) {
  return l1.cumulated_resource_consumption <=> l2.cumulated_resource_consumption;
}

namespace concepts {
template <typename F, typename G, typename RC>
concept ResourceExtensionFunction = std::invocable<F, const G&, RC&, const RC&, graph_edge_descriptor_t<G>>;

template <typename F, typename RC>
concept DominanceFunction = std::invocable<F, const RC&, const RC&>;
}

namespace r_c_shortest_paths_detail {

// r_c_shortest_paths_dispatch function (body/implementation)
template <concepts::VertexListGraph G, concepts::ReadableVertexIndexMap<G> VertexIndexMap,
          concepts::ReadableEdgeIndexMap<G> EdgeIndexMap, class ResourceContainer,
          concepts::ResourceExtensionFunction<G, ResourceContainer> RFunc,
          concepts::DominanceFunction<ResourceContainer> DFunc, class Visitor>
void r_c_shortest_paths_dispatch(const G& g, const VertexIndexMap& vertex_index_map,
                                 const EdgeIndexMap& /*edge_index_map*/, graph_vertex_descriptor_t<G> s,
                                 graph_vertex_descriptor_t<G> t,
                                 // each inner vector corresponds to a pareto-optimal path
                                 std::vector<std::vector<graph_edge_descriptor_t<G>>>& pareto_optimal_solutions,
                                 std::vector<ResourceContainer>& pareto_optimal_resource_containers,
                                 bool b_all_pareto_optimal_solutions,
                                 // to initialize the first label/resource container
                                 // and to carry the type information
                                 const ResourceContainer& rc, RFunc& ref, DFunc& dominance,
                                 // to specify the memory management strategy for the labels
                                 Visitor vis) {
  pareto_optimal_resource_containers.clear();
  pareto_optimal_solutions.clear();

  std::size_t i_label_num = 0;
  using LabelType = r_c_shortest_paths_label<G, ResourceContainer>;
  using Splabel = std::shared_ptr<LabelType>;
  std::priority_queue<Splabel, std::vector<Splabel>, std::greater<>> unprocessed_labels;

  bool b_feasible = true;
  Splabel splabel_first_label =
      std::make_shared<LabelType>(i_label_num++, rc, std::shared_ptr<LabelType>(), graph_edge_descriptor_t<G>(), s);

  unprocessed_labels.push(splabel_first_label);
  auto vec_vertex_labels = vector_property_map(num_vertices(g), vertex_index_map, std::list<Splabel>{});
  vec_vertex_labels[s].push_back(splabel_first_label);
  auto vec_last_valid_positions_for_dominance =
      vector_property_map(num_vertices(g), vertex_index_map, typename std::list<Splabel>::iterator{});
  for (auto v : vertices(g)) {
    put(vec_last_valid_positions_for_dominance, v, vec_vertex_labels[v].begin());
  }
  auto vec_last_valid_index_for_dominance = vector_property_map(num_vertices(g), vertex_index_map, std::size_t{0});
  auto b_vec_vertex_already_checked_for_dominance = vector_property_map(num_vertices(g), vertex_index_map, false);

  while (!unprocessed_labels.empty() && vis.on_enter_loop(unprocessed_labels, g)) {
    Splabel cur_label = unprocessed_labels.top();
    unprocessed_labels.pop();
    vis.on_label_popped(*cur_label, g);
    // an Splabel object in unprocessed_labels and the respective
    // Splabel object in the respective list<Splabel> of
    // vec_vertex_labels share their embedded r_c_shortest_paths_label
    // object to avoid memory leaks, dominated r_c_shortest_paths_label
    // objects are marked and deleted when popped from
    // unprocessed_labels, as they can no longer be deleted at the end
    // of the function; only the Splabel object in unprocessed_labels
    // still references the r_c_shortest_paths_label object this is also
    // for efficiency, because the else branch is executed only if there
    // is a chance that extending the label leads to new undominated
    // labels, which in turn is possible only if the label to be
    // extended is undominated
    if (!cur_label->b_is_dominated) {
      auto i_cur_resident_vertex = cur_label->resident_vertex;
      auto& list_labels_cur_vertex = get(vec_vertex_labels, i_cur_resident_vertex);
      if (list_labels_cur_vertex.size() >= 2 &&
          vec_last_valid_index_for_dominance[i_cur_resident_vertex] < list_labels_cur_vertex.size()) {
        auto outer_iter = list_labels_cur_vertex.begin();
        bool b_outer_iter_at_or_beyond_last_valid_pos_for_dominance = false;
        while (outer_iter != list_labels_cur_vertex.end()) {
          Splabel cur_outer_splabel = *outer_iter;
          auto inner_iter = outer_iter;
          if (!b_outer_iter_at_or_beyond_last_valid_pos_for_dominance &&
              outer_iter == get(vec_last_valid_positions_for_dominance, i_cur_resident_vertex)) {
            b_outer_iter_at_or_beyond_last_valid_pos_for_dominance = true;
          }
          if (!get(b_vec_vertex_already_checked_for_dominance, i_cur_resident_vertex) ||
              b_outer_iter_at_or_beyond_last_valid_pos_for_dominance) {
            ++inner_iter;
          } else {
            inner_iter = get(vec_last_valid_positions_for_dominance, i_cur_resident_vertex);
            ++inner_iter;
          }
          bool b_outer_iter_erased = false;
          while (inner_iter != list_labels_cur_vertex.end()) {
            Splabel cur_inner_splabel = *inner_iter;
            if (dominance(cur_outer_splabel->cumulated_resource_consumption,
                          cur_inner_splabel->cumulated_resource_consumption)) {
              auto buf = inner_iter;
              ++inner_iter;
              list_labels_cur_vertex.erase(buf);
              if (cur_inner_splabel->b_is_processed) {
                cur_inner_splabel.reset();
              } else {
                cur_inner_splabel->b_is_dominated = true;
              }
              continue;
            }

            ++inner_iter;

            if (dominance(cur_inner_splabel->cumulated_resource_consumption,
                          cur_outer_splabel->cumulated_resource_consumption)) {
              auto buf = outer_iter;
              ++outer_iter;
              list_labels_cur_vertex.erase(buf);
              b_outer_iter_erased = true;
              if (cur_outer_splabel->b_is_processed) {
                cur_outer_splabel.reset();
              } else {
                cur_outer_splabel->b_is_dominated = true;
              }
              break;
            }
          }
          if (!b_outer_iter_erased) {
            ++outer_iter;
          }
        }
        if (list_labels_cur_vertex.size() > 1) {
          put(vec_last_valid_positions_for_dominance, i_cur_resident_vertex, (--(list_labels_cur_vertex.end())));
        } else {
          put(vec_last_valid_positions_for_dominance, i_cur_resident_vertex, list_labels_cur_vertex.begin());
        }
        put(b_vec_vertex_already_checked_for_dominance, i_cur_resident_vertex, true);
        put(vec_last_valid_index_for_dominance, i_cur_resident_vertex, list_labels_cur_vertex.size() - 1);
      }
    }
    if (!b_all_pareto_optimal_solutions && cur_label->resident_vertex == t) {
      // the devil don't sleep
      if (cur_label->b_is_dominated) {
        cur_label.reset();
      }
      while (unprocessed_labels.size()) {
        Splabel l = unprocessed_labels.top();
        unprocessed_labels.pop();
        // delete only dominated labels, because nondominated labels
        // are deleted at the end of the function
        if (l->b_is_dominated) {
          l.reset();
        }
      }
      break;
    }
    if (!cur_label->b_is_dominated) {
      cur_label->b_is_processed = true;
      vis.on_label_not_dominated(*cur_label, g);
      auto cur_vertex = cur_label->resident_vertex;
      for (auto e : out_edges(cur_vertex, g)) {
        b_feasible = true;
        Splabel new_label = std::make_shared<LabelType>(i_label_num++, cur_label->cumulated_resource_consumption,
                                                        cur_label, e, target(e, g));
        b_feasible = ref(g, new_label->cumulated_resource_consumption,
                         new_label->p_pred_label->cumulated_resource_consumption, new_label->pred_edge);

        if (!b_feasible) {
          vis.on_label_not_feasible(*new_label, g);
          new_label.reset();
        } else {
          vis.on_label_feasible(*new_label, g);
          vec_vertex_labels[new_label->resident_vertex].push_back(new_label);
          unprocessed_labels.push(new_label);
        }
      }
    } else {
      vis.on_label_dominated(*cur_label, g);
      cur_label.reset();
    }
  }
  std::list<Splabel> dsplabels = get(vec_vertex_labels, t);
  if (!b_all_pareto_optimal_solutions) {
    dsplabels.sort();
  }
  // if d could be reached from o
  if (!dsplabels.empty()) {
    for (auto csi : dsplabels) {
      std::vector<graph_edge_descriptor_t<G>> cur_pareto_optimal_path;
      auto p_cur_label = csi;
      pareto_optimal_resource_containers.push_back(p_cur_label->cumulated_resource_consumption);
      while (p_cur_label->num != 0) {
        cur_pareto_optimal_path.push_back(p_cur_label->pred_edge);
        p_cur_label = p_cur_label->p_pred_label;

        // assertion b_is_valid beyond this point is not correct if
        // the domination function requires resource levels to be
        // strictly greater than existing values
        //
        // Example
        // Customers
        // id   min_arrival   max_departure
        //  2             0             974
        //  3             0             972
        //  4             0             964
        //  5           678             801
        //
        // Path A: 2-3-4-5 (times: 0-16-49-84-678)
        // Path B: 3-2-4-5 (times: 0-18-51-62-678)
        // The partial path 3-2-4 dominates the other partial path
        // 2-3-4, though the path 3-2-4-5 does not strictly dominate
        // the path 2-3-4-5
      }
      pareto_optimal_solutions.push_back(cur_pareto_optimal_path);
      if (!b_all_pareto_optimal_solutions) {
        break;
      }
    }
  }

  for (auto i : vertices(g)) {
    for (auto si : vec_vertex_labels[i]) {
      si.reset();
    }
  }
}  // r_c_shortest_paths_dispatch

}  // namespace r_c_shortest_paths_detail

// default_r_c_shortest_paths_visitor struct
struct default_r_c_shortest_paths_visitor {
  template <class Label, class Graph>
  void on_label_popped(const Label& /*unused*/, const Graph& /*unused*/) {}
  template <class Label, class Graph>
  void on_label_feasible(const Label& /*unused*/, const Graph& /*unused*/) {}
  template <class Label, class Graph>
  void on_label_not_feasible(const Label& /*unused*/, const Graph& /*unused*/) {}
  template <class Label, class Graph>
  void on_label_dominated(const Label& /*unused*/, const Graph& /*unused*/) {}
  template <class Label, class Graph>
  void on_label_not_dominated(const Label& /*unused*/, const Graph& /*unused*/) {}
  template <class Queue, class Graph>
  bool on_enter_loop(const Queue& /*queue*/, const Graph& /*graph*/) {
    return true;
  }
};  // default_r_c_shortest_paths_visitor

// default_r_c_shortest_paths_allocator
using default_r_c_shortest_paths_allocator = std::allocator<int>;
// default_r_c_shortest_paths_allocator

// r_c_shortest_paths functions (handle/interface)
// first overload:
// - return all pareto-optimal solutions
// - specify Visitor arguments
template <concepts::VertexListGraph G, concepts::ReadableVertexIndexMap<G> VertexIndexMap,
          concepts::ReadableEdgeIndexMap<G> EdgeIndexMap, class ResourceContainer,
          concepts::ResourceExtensionFunction<G, ResourceContainer> RFunc,
          concepts::DominanceFunction<ResourceContainer> DFunc, class Visitor>
void r_c_shortest_paths(const G& g, const VertexIndexMap& vertex_index_map, const EdgeIndexMap& edge_index_map,
                        graph_vertex_descriptor_t<G> s, graph_vertex_descriptor_t<G> t,
                        // each inner vector corresponds to a pareto-optimal path
                        std::vector<std::vector<graph_edge_descriptor_t<G>>>& pareto_optimal_solutions,
                        std::vector<ResourceContainer>& pareto_optimal_resource_containers,
                        // to initialize the first label/resource container
                        // and to carry the type information
                        const ResourceContainer& rc, const RFunc& ref, const DFunc& dominance,
                        // to specify the memory management strategy for the labels
                        Visitor vis) {
  r_c_shortest_paths_detail::r_c_shortest_paths_dispatch(g, vertex_index_map, edge_index_map, s, t,
                                                         pareto_optimal_solutions, pareto_optimal_resource_containers,
                                                         true, rc, ref, dominance, vis);
}

// second overload:
// - return only one pareto-optimal solution
// - specify Visitor arguments
template <concepts::VertexListGraph G, concepts::ReadableVertexIndexMap<G> VertexIndexMap,
          concepts::ReadableEdgeIndexMap<G> EdgeIndexMap, class ResourceContainer,
          concepts::ResourceExtensionFunction<G, ResourceContainer> RFunc,
          concepts::DominanceFunction<ResourceContainer> DFunc, class Visitor>
void r_c_shortest_paths(const G& g, const VertexIndexMap& vertex_index_map, const EdgeIndexMap& edge_index_map,
                        graph_vertex_descriptor_t<G> s, graph_vertex_descriptor_t<G> t,
                        std::vector<graph_edge_descriptor_t<G>>& pareto_optimal_solution,
                        ResourceContainer& pareto_optimal_resource_container,
                        // to initialize the first label/resource container
                        // and to carry the type information
                        const ResourceContainer& rc, const RFunc& ref, const DFunc& dominance,
                        // to specify the memory management strategy for the labels
                        Visitor vis) {
  // each inner vector corresponds to a pareto-optimal path
  std::vector<std::vector<graph_edge_descriptor_t<G>>> pareto_optimal_solutions;
  std::vector<ResourceContainer> pareto_optimal_resource_containers;
  r_c_shortest_paths_detail::r_c_shortest_paths_dispatch(g, vertex_index_map, edge_index_map, s, t,
                                                         pareto_optimal_solutions, pareto_optimal_resource_containers,
                                                         false, rc, ref, dominance, vis);
  if (!pareto_optimal_solutions.empty()) {
    pareto_optimal_solution = pareto_optimal_solutions[0];
    pareto_optimal_resource_container = pareto_optimal_resource_containers[0];
  }
}

// third overload:
// - return all pareto-optimal solutions
// - use default Visitor
template <concepts::VertexListGraph G, concepts::ReadableVertexIndexMap<G> VertexIndexMap,
          concepts::ReadableEdgeIndexMap<G> EdgeIndexMap, class ResourceContainer,
          concepts::ResourceExtensionFunction<G, ResourceContainer> RFunc,
          concepts::DominanceFunction<ResourceContainer> DFunc>
void r_c_shortest_paths(
    const G& g, const VertexIndexMap& vertex_index_map, const EdgeIndexMap& edge_index_map,
    graph_vertex_descriptor_t<G> s, graph_vertex_descriptor_t<G> t,
    // each inner vector corresponds to a pareto-optimal path
    std::vector<std::vector<graph_edge_descriptor_t<G>>>& pareto_optimal_solutions,
    std::vector<ResourceContainer>& pareto_optimal_resource_containers,
    // to initialize the first label/resource container
    // and to carry the type information
    const ResourceContainer& rc, const RFunc& ref, const DFunc& dominance) {
  r_c_shortest_paths_detail::r_c_shortest_paths_dispatch(
      g, vertex_index_map, edge_index_map, s, t, pareto_optimal_solutions, pareto_optimal_resource_containers, true, rc,
      ref, dominance, default_r_c_shortest_paths_visitor());
}

// fourth overload:
// - return only one pareto-optimal solution
// - use default Visitor
template <concepts::VertexListGraph G, concepts::ReadableVertexIndexMap<G> VertexIndexMap,
          concepts::ReadableEdgeIndexMap<G> EdgeIndexMap, class ResourceContainer,
          concepts::ResourceExtensionFunction<G, ResourceContainer> RFunc,
          concepts::DominanceFunction<ResourceContainer> DFunc>
void r_c_shortest_paths(const G& g, const VertexIndexMap& vertex_index_map, const EdgeIndexMap& edge_index_map,
                        graph_vertex_descriptor_t<G> s,
                        graph_vertex_descriptor_t<G> t,
                        std::vector<graph_edge_descriptor_t<G>>& pareto_optimal_solution,
                        ResourceContainer& pareto_optimal_resource_container,
                        // to initialize the first label/resource container
                        // and to carry the type information
                        const ResourceContainer& rc, const RFunc& ref,
                        const DFunc& dominance) {
  // each inner vector corresponds to a pareto-optimal path
  std::vector<std::vector<graph_edge_descriptor_t<G>>> pareto_optimal_solutions;
  std::vector<ResourceContainer> pareto_optimal_resource_containers;
  r_c_shortest_paths_detail::r_c_shortest_paths_dispatch(
      g, vertex_index_map, edge_index_map, s, t, pareto_optimal_solutions, pareto_optimal_resource_containers, false,
      rc, ref, dominance, default_r_c_shortest_paths_allocator(), default_r_c_shortest_paths_visitor());
  if (!pareto_optimal_solutions.empty()) {
    pareto_optimal_solution = pareto_optimal_solutions[0];
    pareto_optimal_resource_container = pareto_optimal_resource_containers[0];
  }
}
// r_c_shortest_paths

// check_r_c_path function
template <concepts::VertexListGraph G, class ResourceContainer,
          concepts::ResourceExtensionFunction<G, ResourceContainer> RFunc>
void check_r_c_path(const G& g, const std::vector<graph_edge_descriptor_t<G>>& ed_vec_path,
                    const ResourceContainer& initial_resource_levels,
                    // if true, computed accumulated final resource levels must
                    // be equal to desired_final_resource_levels
                    // if false, computed accumulated final resource levels must
                    // be less than or equal to desired_final_resource_levels
                    bool b_result_must_be_equal_to_desired_final_resource_levels,
                    const ResourceContainer& desired_final_resource_levels,
                    ResourceContainer& actual_final_resource_levels, const RFunc& ref,
                    bool& b_is_a_path_at_all, bool& b_feasible, bool& b_correctly_extended,
                    graph_edge_descriptor_t<G>& ed_last_extended_arc) {
  size_t i_size_ed_vec_path = ed_vec_path.size();
  std::vector<graph_edge_descriptor_t<G>> buf_path;
  if (i_size_ed_vec_path == 0) {
    b_feasible = true;
  } else {
    if (i_size_ed_vec_path == 1 || target(ed_vec_path[0], g) == source(ed_vec_path[1], g)) {
      buf_path = ed_vec_path;
    } else {
      for (size_t i = i_size_ed_vec_path; i > 0; --i) {
        buf_path.push_back(ed_vec_path[i - 1]);
      }
    }
    for (size_t i = 0; i < i_size_ed_vec_path - 1; ++i) {
      if (target(buf_path[i], g) != source(buf_path[i + 1], g)) {
        b_is_a_path_at_all = false;
        b_feasible = false;
        b_correctly_extended = false;
        return;
      }
    }
  }
  b_is_a_path_at_all = true;
  b_feasible = true;
  b_correctly_extended = false;
  ResourceContainer current_resource_levels = initial_resource_levels;
  actual_final_resource_levels = current_resource_levels;
  for (size_t i = 0; i < i_size_ed_vec_path; ++i) {
    ed_last_extended_arc = buf_path[i];
    b_feasible = ref(g, actual_final_resource_levels, current_resource_levels, buf_path[i]);
    current_resource_levels = actual_final_resource_levels;
    if (!b_feasible) {
      return;
    }
  }
  if (b_result_must_be_equal_to_desired_final_resource_levels) {
    b_correctly_extended = static_cast<bool>(actual_final_resource_levels == desired_final_resource_levels);
  } else {
    if (actual_final_resource_levels < desired_final_resource_levels ||
        actual_final_resource_levels == desired_final_resource_levels) {
      b_correctly_extended = true;
    }
  }
}  // check_path

}  // namespace bagl

#endif  // BAGL_BAGL_R_C_SHORTEST_PATHS_H_
