// Copyright (c) 2018 Yi Ji
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_MAXIMUM_WEIGHTED_MATCHING_H_
#define BAGL_BAGL_MAXIMUM_WEIGHTED_MATCHING_H_

#include <algorithm>
#include <limits>
#include <memory>
#include <tuple>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/max_cardinality_matching.h"
#include "bagl/partial_range.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

namespace bagl {

template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> MateMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap, concepts::ReadableEdgePropertyMap<Graph> WeightMap>
auto matching_weight_sum(const Graph& g, MateMap mate, VertexIndexMap vm, WeightMap weight) {
  using edge_weight_t = property_traits_value_t<WeightMap>;

  edge_weight_t weight_sum = 0;
  for (auto v : vertices(g)) {
    if (get(mate, v) != graph_traits<Graph>::null_vertex() && get(vm, v) < get(vm, get(mate, v))) {
      weight_sum += get(weight, edge(v, mate[v], g).first);
    }
  }
  return weight_sum;
}

template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> MateMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
auto matching_weight_sum(const Graph& g, MateMap mate, VertexIndexMap vm) {
  return matching_weight_sum(g, mate, vm, get(edge_weight, g));
}

template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> MateMap>
auto matching_weight_sum(const Graph& g, MateMap mate) {
  return matching_weight_sum(g, mate, get(vertex_index, g));
}

template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> MateMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap, concepts::ReadableEdgePropertyMap<Graph> WeightMap>
class weighted_augmenting_path_finder {
 public:
  template <typename T>
  using map_vertex_to = vector_property_map<T, VertexIndexMap>;
  using vertex_state_t = max_cardinality_detail::vertex_state;
  using vertex_descriptor_t = graph_vertex_descriptor_t<Graph>;
  using edge_descriptor_t = graph_edge_descriptor_t<Graph>;
  using edge_weight_t = property_traits_value_t<WeightMap>;
  static constexpr edge_weight_t min_edge_weight = std::numeric_limits<edge_weight_t>::min();
  static constexpr edge_weight_t max_edge_weight = std::numeric_limits<edge_weight_t>::max();
  using vertex_list_t = std::deque<vertex_descriptor_t>;
  using edge_list_t = std::vector<edge_descriptor_t>;
  using vertex_to_vertex_map_t = map_vertex_to<vertex_descriptor_t>;
  using vertex_to_weight_map_t = map_vertex_to<edge_weight_t>;
  using vertex_to_bool_map_t = map_vertex_to<bool>;
  using vertex_to_pair_map_t = map_vertex_to<std::pair<vertex_descriptor_t, vertex_descriptor_t>>;
  using vertex_to_edge_map_t = map_vertex_to<std::pair<edge_descriptor_t, bool>>;
  using vertex_pair_to_edge_map_t = map_vertex_to<vertex_to_edge_map_t>;

  class blossom {
   public:
    using blossom_ptr_t = std::shared_ptr<blossom>;
    std::vector<blossom_ptr_t> sub_blossoms;
    edge_weight_t dual_var;
    blossom_ptr_t father;

    blossom() : dual_var(0), father(blossom_ptr_t()) {}

    // get the base vertex of a blossom by recursively getting
    // its base sub-blossom, which is always the first one in
    // sub_blossoms because of how we create and maintain blossoms
    virtual vertex_descriptor_t get_base() const {
      const blossom* b = this;
      while (!b->sub_blossoms.empty()) {
        b = b->sub_blossoms[0].get();
      }
      return b->get_base();
    }

    // set a sub-blossom as a blossom's base by exchanging it
    // with its first sub-blossom
    void set_base(const blossom_ptr_t& sub) {
      for (auto& bi : sub_blossoms) {
        if (sub == bi) {
          std::swap(sub_blossoms.front(), bi);
          break;
        }
      }
    }

    // get all vertices inside recursively
    virtual std::vector<vertex_descriptor_t> vertices() const {
      std::vector<vertex_descriptor_t> all_vertices;
      for (const auto& bi : sub_blossoms) {
        std::vector<vertex_descriptor_t> some_vertices = bi->vertices();
        all_vertices.insert(all_vertices.end(), some_vertices.begin(), some_vertices.end());
      }
      return all_vertices;
    }
  };

  // a trivial_blossom only has one vertex and no sub-blossom;
  // for each vertex v, in_blossom[v] is the trivial_blossom that contains it
  // directly
  class trivial_blossom : public blossom {
   public:
    explicit trivial_blossom(vertex_descriptor_t v) : trivial_vertex_(v) {}
    virtual vertex_descriptor_t get_base() const { return trivial_vertex_; }

    virtual std::vector<vertex_descriptor_t> vertices() const {
      std::vector<vertex_descriptor_t> all_vertices;
      all_vertices.push_back(trivial_vertex_);
      return all_vertices;
    }

   private:
    vertex_descriptor_t trivial_vertex_;
  };

  using blossom_ptr_t = std::shared_ptr<blossom>;
  using blossom_iterator_t = typename std::vector<blossom_ptr_t>::iterator;
  using vertex_to_blossom_map_t = map_vertex_to<blossom_ptr_t>;

  weighted_augmenting_path_finder(const Graph& arg_g, MateMap arg_mate, VertexIndexMap arg_vm, WeightMap arg_weight)
      : g_(arg_g),
        vm_(arg_vm),
        weight_(arg_weight),
        null_edge_(num_edges(g_) == 0 ? edge_descriptor_t() : *edges(g_).begin(), false),
        mate_(num_vertices(g_), vm_, vm_, graph_traits<Graph>::null_vertex()),
        label_s_(num_vertices(g_), vm_, graph_traits<Graph>::null_vertex()),
        label_t_(num_vertices(g_), vm_, graph_traits<Graph>::null_vertex()),
        outlet_(num_vertices(g_), vm_, graph_traits<Graph>::null_vertex()),
        tau_idx_(num_vertices(g_), vm_, graph_traits<Graph>::null_vertex()),
        dual_var_(num_vertices(g_), vm_, min_edge_weight),
        pi_(num_vertices(g_), vm_, max_edge_weight),
        gamma_(num_vertices(g_), vm_, max_edge_weight),
        tau_(num_vertices(g_), vm_, max_edge_weight),
        in_blossom_(num_vertices(g_), vm_, blossom_ptr_t{}),
        old_label_(num_vertices(g_), vm_,
                   std::pair(graph_traits<Graph>::null_vertex(), graph_traits<Graph>::null_vertex())),
        critical_edge_(
            vertex_pair_to_edge_map_t(num_vertices(g_), vm_, vertex_to_edge_map_t(num_vertices(g_), vm_, null_edge_))) {
    auto max_weight = min_edge_weight;
    for (auto e : edges(g_)) {
      max_weight = std::max(max_weight, get(weight_, e));
    }

    for (auto u : vertices(g_)) {
      mate_[u] = get(arg_mate, u);
      dual_var_[u] = 2 * max_weight;
      in_blossom_[u] = std::make_shared<trivial_blossom>(u);
      outlet_[u] = u;
    }

    init();
  }

  // return the top blossom where v is contained inside
  blossom_ptr_t in_top_blossom(vertex_descriptor_t v) const {
    blossom_ptr_t b = in_blossom_[v];
    while (b->father != blossom_ptr_t()) {
      b = b->father;
    }
    return b;
  }

  // check if vertex v is in blossom b
  bool is_in_blossom(blossom_ptr_t b, vertex_descriptor_t v) const {
    if (v == graph_traits<Graph>::null_vertex()) {
      return false;
    }
    blossom_ptr_t vb = in_blossom_[v]->father;
    while (vb != blossom_ptr_t()) {
      if (vb.get() == b.get()) {
        return true;
      }
      vb = vb->father;
    }
    return false;
  }

  // return the base vertex of the top blossom that contains v
  vertex_descriptor_t base_vertex(vertex_descriptor_t v) const { return in_top_blossom(v)->get_base(); }

  // add an existed top blossom of base vertex v into new top
  // blossom b as its sub-blossom
  void add_sub_blossom(blossom_ptr_t b, vertex_descriptor_t v) {
    blossom_ptr_t sub = in_top_blossom(v);
    sub->father = b;
    b->sub_blossoms.push_back(sub);
    if (sub->sub_blossoms.empty()) {
      return;
    }
    top_blossoms_.erase(std::find(top_blossoms_.begin(), top_blossoms_.end(), sub));
  }

  // when a top blossom is created or its base vertex getting an S-label,
  // add all edges incident to this blossom into even_edges
  void bloom(blossom_ptr_t b) {
    auto vertices_of_b = b->vertices();
    for (auto v : vertices_of_b) {
      for (auto e : out_edges(v, g_)) {
        if (target(e, g_) != v && mate_[v] != target(e, g_)) {
          even_edges_.push_back(e);
        }
      }
    }
  }

  // assigning a T-label to a non S-vertex, along with outlet and updating pi
  // value if updated pi[v] equals zero, augment the matching from its mate
  // vertex
  void put_t_label(vertex_descriptor_t v, vertex_descriptor_t T_label, vertex_descriptor_t outlet_v,
                   edge_weight_t pi_v) {
    if (label_s_[v] != graph_traits<Graph>::null_vertex()) {
      return;
    }

    label_t_[v] = T_label;
    outlet_[v] = outlet_v;
    pi_[v] = pi_v;

    vertex_descriptor_t v_mate = mate_[v];
    if (pi_[v] == 0) {
      label_t_[v_mate] = graph_traits<Graph>::null_vertex();
      label_s_[v_mate] = v;
      bloom(in_top_blossom(v_mate));
    }
  }

  // get the missing T-label for a to-be-expanded base vertex
  // the missing T-label is the last vertex of the path from outlet[v] to v
  std::pair<vertex_descriptor_t, vertex_descriptor_t> missing_label(vertex_descriptor_t b_base) {
    vertex_descriptor_t missing_outlet = outlet_[b_base];

    if (outlet_[b_base] == b_base) {
      return {graph_traits<Graph>::null_vertex(), missing_outlet};
    }

    for (auto v : vertices(g_)) {
      old_label_[v] = {label_t_[v], outlet_[v]};
    }

    auto child = std::pair{outlet_[b_base], vertex_state_t::even};
    blossom_ptr_t b = in_blossom_[child.first];
    for (; b->father->father != nullptr; b = b->father) {
    }
    child.first = b->get_base();

    if (child.first == b_base) {
      return {graph_traits<Graph>::null_vertex(), missing_outlet};
    }

    while (true) {
      auto [child_parent, child_parent_state] = parent(child, true);
      for (b = in_blossom_[child_parent]; b->father->father != nullptr; b = b->father) {
      }
      missing_outlet = child_parent;
      child_parent = b->get_base();
      if (child_parent == b_base) {
        break;
      }
      child = {child_parent, child_parent_state};
    }
    return std::make_pair(child.first, missing_outlet);
  }

  // expand a top blossom, put all its non-trivial sub-blossoms into
  // top_blossoms
  blossom_iterator_t expand_blossom(blossom_iterator_t bi, std::vector<blossom_ptr_t>& new_ones) {
    blossom_ptr_t b = *bi;
    for (const blossom_ptr_t& sub_blossom : b->sub_blossoms) {
      vertex_descriptor_t sub_base = sub_blossom->get_base();
      label_s_[sub_base] = label_t_[sub_base] = graph_traits<Graph>::null_vertex();
      outlet_[sub_base] = sub_base;
      sub_blossom->father = blossom_ptr_t{};
      // new top blossoms cannot be pushed back into top_blossoms
      // immediately, because push_back() may cause reallocation and then
      // invalid iterators
      if (!sub_blossom->sub_blossoms.empty()) {
        new_ones.emplace_back(sub_blossom);
      }
    }
    return top_blossoms_.erase(bi);
  }

  // when expanding a T-blossom with base v, it requires more operations:
  // supply the missing T-labels for new base vertices by picking the minimum
  // tau from vertices of each corresponding new top-blossoms; when label_T[v]
  // is null or we have a smaller tau from missing_label(v), replace T-label
  // and outlet of v (but don't bloom v)
  blossom_iterator_t expand_t_blossom(blossom_iterator_t bi, std::vector<blossom_ptr_t>& new_ones) {
    blossom_ptr_t b = *bi;

    vertex_descriptor_t b_base = b->get_base();
    auto t_and_outlet = missing_label(b_base);

    blossom_iterator_t next_bi = expand_blossom(bi, new_ones);

    for (const blossom_ptr_t& sub_blossom : b->sub_blossoms) {
      vertex_descriptor_t sub_base = sub_blossom->get_base();
      vertex_descriptor_t min_tau_v = graph_traits<Graph>::null_vertex();
      edge_weight_t min_tau = max_edge_weight;

      std::vector<vertex_descriptor_t> sub_vertices = sub_blossom->vertices();
      for (auto v : sub_vertices) {
        if (tau_[v] < min_tau) {
          min_tau = tau_[v];
          min_tau_v = v;
        }
      }

      if (min_tau < max_edge_weight) {
        put_t_label(sub_base, tau_idx_[min_tau_v], min_tau_v, tau_[min_tau_v]);
      }
    }

    if (label_t_[b_base] == graph_traits<Graph>::null_vertex() || tau_[old_label_[b_base].second] < pi_[b_base]) {
      std::tie(label_t_[b_base], outlet_[b_base]) = t_and_outlet;
    }

    return next_bi;
  }

  // when vertices v and w are matched to each other by augmenting,
  // we must set v/w as base vertex of any blossom who contains v/w and
  // is a sub-blossom of their lowest (smallest) common blossom
  void adjust_blossom(vertex_descriptor_t v, vertex_descriptor_t w) {
    blossom_ptr_t vb = in_blossom_[v];
    blossom_ptr_t wb = in_blossom_[w];
    blossom_ptr_t lowest_common_blossom;
    std::vector<blossom_ptr_t> v_ancestors;
    std::vector<blossom_ptr_t> w_ancestors;

    while (vb->father != blossom_ptr_t()) {
      v_ancestors.push_back(vb->father);
      vb = vb->father;
    }
    while (wb->father != blossom_ptr_t()) {
      w_ancestors.push_back(wb->father);
      wb = wb->father;
    }

    auto i = v_ancestors.rbegin();
    auto j = w_ancestors.rbegin();
    while (i != v_ancestors.rend() && j != w_ancestors.rend() && i->get() == j->get()) {
      lowest_common_blossom = *i;
      ++i;
      ++j;
    }

    vb = in_blossom_[v];
    wb = in_blossom_[w];
    while (vb->father != lowest_common_blossom) {
      vb->father->set_base(vb);
      vb = vb->father;
    }
    while (wb->father != lowest_common_blossom) {
      wb->father->set_base(wb);
      wb = wb->father;
    }
  }

  // every edge weight is multiplied by 4 to ensure integer weights
  // throughout the algorithm if all input weights are integers
  edge_weight_t slack(const edge_descriptor_t& e) const {
    auto v = source(e, g_);
    auto w = target(e, g_);
    return dual_var_[v] + dual_var_[w] - 4 * get(weight_, e);
  }

  // backtrace one step on vertex v along the augmenting path
  // by its labels and its vertex state;
  // boolean parameter "use_old" means whether we are updating labels,
  // if we are, then we use old labels to backtrace and also we
  // don't jump to its base vertex when we reach an odd vertex
  std::pair<vertex_descriptor_t, vertex_state_t> parent(std::pair<vertex_descriptor_t, vertex_state_t> v,
                                                        bool use_old = false) const {
    if (v.second == vertex_state_t::even) {
      // a paranoid check: label_S shoule be the same as mate in
      // backtracing
      if (label_s_[v.first] == graph_traits<Graph>::null_vertex()) {
        label_s_[v.first] = mate_[v.first];
      }
      return {label_s_[v.first], vertex_state_t::odd};
    }
    if (v.second == vertex_state_t::odd) {
      vertex_descriptor_t w = use_old ? old_label_[v.first].first : base_vertex(label_t_[v.first]);
      return {w, vertex_state_t::even};
    }
    return {v.first, vertex_state_t::unreached};
  }

  // backtrace from vertices v and w to their free (unmatched) ancesters,
  // return the nearest common ancestor (null_vertex if none) of v and w
  vertex_descriptor_t nearest_common_ancestor(vertex_descriptor_t v, vertex_descriptor_t w,
                                              vertex_descriptor_t& v_free_ancestor,
                                              vertex_descriptor_t& w_free_ancestor) const {
    auto v_up = std::pair{v, vertex_state_t::even};
    auto w_up = std::pair{w, vertex_state_t::even};
    vertex_descriptor_t nca = graph_traits<Graph>::null_vertex();
    w_free_ancestor = graph_traits<Graph>::null_vertex();
    v_free_ancestor = graph_traits<Graph>::null_vertex();

    std::vector<bool> ancestor_of_w_vector(num_vertices(g_), false);
    std::vector<bool> ancestor_of_v_vector(num_vertices(g_), false);
    vertex_to_bool_map_t ancestor_of_w(ancestor_of_w_vector.begin(), vm_);
    vertex_to_bool_map_t ancestor_of_v(ancestor_of_v_vector.begin(), vm_);

    while (nca == graph_traits<Graph>::null_vertex() && (v_free_ancestor == graph_traits<Graph>::null_vertex() ||
                                                         w_free_ancestor == graph_traits<Graph>::null_vertex())) {
      ancestor_of_w[w_up.first] = true;
      ancestor_of_v[v_up.first] = true;

      if (w_free_ancestor == graph_traits<Graph>::null_vertex()) {
        w_up = parent(w_up);
      }
      if (v_free_ancestor == graph_traits<Graph>::null_vertex()) {
        v_up = parent(v_up);
      }

      if (mate_[v_up.first] == graph_traits<Graph>::null_vertex()) {
        v_free_ancestor = v_up.first;
      }
      if (mate_[w_up.first] == graph_traits<Graph>::null_vertex()) {
        w_free_ancestor = w_up.first;
      }

      if (ancestor_of_w[v_up.first] || v_up.first == w_up.first) {
        nca = v_up.first;
      } else if (ancestor_of_v[w_up.first]) {
        nca = w_up.first;
      } else if (v_free_ancestor == w_free_ancestor && v_free_ancestor != graph_traits<Graph>::null_vertex()) {
        nca = v_up.first;
      }
    }

    return nca;
  }

  // when a new top blossom b is created by connecting (v, w), we add
  // sub-blossoms into b along backtracing from v_prime and w_prime to
  // stop_vertex (the base vertex); also, we set labels and outlet for each
  // base vertex we pass by
  void make_blossom(blossom_ptr_t b, vertex_descriptor_t w_prime, vertex_descriptor_t v_prime,
                    vertex_descriptor_t stop_vertex) {
    auto u = std::pair{v_prime, vertex_state_t::odd};
    auto u_up = std::pair{w_prime, vertex_state_t::even};

    for (; u_up.first != stop_vertex; u = u_up, u_up = parent(u)) {
      if (u_up.second == vertex_state_t::even) {
        if (!in_top_blossom(u_up.first)->sub_blossoms.empty()) {
          outlet_[u_up.first] = label_t_[u.first];
        }
        label_t_[u_up.first] = outlet_[u.first];
      } else if (u_up.second == vertex_state_t::odd) {
        label_s_[u_up.first] = u.first;
      }

      add_sub_blossom(b, u_up.first);
    }
  }

  // the design of recursively expanding augmenting path in
  // (reversed_)retrieve_augmenting_path functions is inspired by same
  // functions in max_cardinality_matching.hpp; except that in weighted
  // matching, we use "outlet" vertices instead of "bridge" vertex pairs: if
  // blossom b is the smallest non-trivial blossom that contains its base
  // vertex v, then v and outlet[v] are where augmenting path enters and
  // leaves b
  void retrieve_augmenting_path(vertex_descriptor_t v, vertex_descriptor_t w, vertex_state_t v_state) {
    if (v == w) {
      aug_path_.push_back(v);
    } else if (v_state == vertex_state_t::even) {
      aug_path_.push_back(v);
      retrieve_augmenting_path(label_s_[v], w, vertex_state_t::odd);
    } else if (v_state == vertex_state_t::odd) {
      if (outlet_[v] == v) {
        aug_path_.push_back(v);
      } else {
        reversed_retrieve_augmenting_path(outlet_[v], v, vertex_state_t::even);
      }
      retrieve_augmenting_path(label_t_[v], w, vertex_state_t::even);
    }
  }

  void reversed_retrieve_augmenting_path(vertex_descriptor_t v, vertex_descriptor_t w, vertex_state_t v_state) {
    if (v == w) {
      aug_path_.push_back(v);
    } else if (v_state == vertex_state_t::even) {
      reversed_retrieve_augmenting_path(label_s_[v], w, vertex_state_t::odd);
      aug_path_.push_back(v);
    } else if (v_state == vertex_state_t::odd) {
      reversed_retrieve_augmenting_path(label_t_[v], w, vertex_state_t::even);
      if (outlet_[v] != v) {
        retrieve_augmenting_path(outlet_[v], v, vertex_state_t::even);
      } else {
        aug_path_.push_back(v);
      }
    }
  }

  // correct labels for vertices in the augmenting path
  void relabel(vertex_descriptor_t v) {
    blossom_ptr_t b = in_blossom_[v]->father;

    if (!is_in_blossom(b, mate_[v])) {  // if v is a new base vertex
      auto u = std::pair{v, vertex_state_t::even};
      while (label_s_[u.first] != u.first && is_in_blossom(b, label_s_[u.first])) {
        u = parent(u, true);
      }

      vertex_descriptor_t old_base = u.first;
      if (label_s_[old_base] != old_base) {  // if old base is not exposed
        label_t_[v] = label_s_[old_base];
        outlet_[v] = old_base;
      } else {  // if old base is exposed then new label_T[v] is not in b,
        // we must (i) make b2 the smallest blossom containing v but not
        // as base vertex (ii) backtrace from b2's new base vertex to b
        label_t_[v] = graph_traits<Graph>::null_vertex();
        for (b = b->father; b != nullptr && b->get_base() == v; b = b->father) {
        }
        if (b != nullptr) {
          u = {b->get_base(), vertex_state_t::odd};
          while (!is_in_blossom(in_blossom_[v]->father, old_label_[u.first].first)) {
            u = parent(u, true);
          }
          label_t_[v] = u.first;
          outlet_[v] = old_label_[u.first].first;
        }
      }
    } else if (label_s_[v] == v || !is_in_blossom(b, label_s_[v])) {
      // if v is an old base vertex
      // let u be the new base vertex; backtrace from u's old T-label
      auto u = std::pair{b->get_base(), vertex_state_t::odd};
      while (old_label_[u.first].first != graph_traits<Graph>::null_vertex() && old_label_[u.first].first != v) {
        u = parent(u, true);
      }
      label_t_[v] = old_label_[u.first].second;
      outlet_[v] = v;
    } else {
      // if v is neither a new nor an old base vertex
      label_t_[v] = label_s_[v];
    }
  }

  void augmenting(vertex_descriptor_t v, vertex_descriptor_t v_free_ancestor, vertex_descriptor_t w,
                  vertex_descriptor_t w_free_ancestor) {
    // retrieve the augmenting path and put it in aug_path
    reversed_retrieve_augmenting_path(v, v_free_ancestor, vertex_state_t::even);
    retrieve_augmenting_path(w, w_free_ancestor, vertex_state_t::even);

    // augment the matching along aug_path
    vertex_descriptor_t a;
    vertex_descriptor_t b;
    vertex_list_t reversed_aug_path;
    while (!aug_path_.empty()) {
      a = aug_path_.front();
      aug_path_.pop_front();
      reversed_aug_path.push_back(a);
      b = aug_path_.front();
      aug_path_.pop_front();
      reversed_aug_path.push_back(b);

      mate_[a] = b;
      mate_[b] = a;

      // reset base vertex for every blossom in augment path
      adjust_blossom(a, b);
    }

    for (auto u : vertices(g_)) {
      old_label_[u] = {label_t_[u], outlet_[u]};
    }

    // correct labels for in-blossom vertices along aug_path
    while (!reversed_aug_path.empty()) {
      a = reversed_aug_path.front();
      reversed_aug_path.pop_front();

      if (in_blossom_[a]->father != blossom_ptr_t{}) {
        relabel(a);
      }
    }

    for (auto u : vertices(g_)) {
      if (mate_[u] != graph_traits<Graph>::null_vertex()) {
        label_s_[u] = mate_[u];
      }
    }

    // expand blossoms with zero dual variables
    std::vector<blossom_ptr_t> new_top_blossoms;
    for (auto bi = top_blossoms_.begin(); bi != top_blossoms_.end();) {
      if ((*bi)->dual_var <= 0) {
        bi = expand_blossom(bi, new_top_blossoms);
      } else {
        ++bi;
      }
    }
    top_blossoms_.insert(top_blossoms_.end(), new_top_blossoms.begin(), new_top_blossoms.end());
    init();
  }

  // create a new blossom and set labels for vertices inside
  void blossoming(vertex_descriptor_t v, vertex_descriptor_t v_prime, vertex_descriptor_t w,
                  vertex_descriptor_t w_prime, vertex_descriptor_t nca) {
    vertex_to_bool_map_t is_old_base(num_vertices(g_), vm_, false);
    for (auto i : vertices(g_)) {
      if (i == base_vertex(i)) {
        is_old_base[i] = true;
      }
    }

    blossom_ptr_t b = std::make_shared<blossom>();
    add_sub_blossom(b, nca);

    label_t_[w_prime] = v;
    label_t_[v_prime] = w;
    outlet_[w_prime] = w;
    outlet_[v_prime] = v;

    make_blossom(b, w_prime, v_prime, nca);
    make_blossom(b, v_prime, w_prime, nca);

    label_t_[nca] = graph_traits<Graph>::null_vertex();
    outlet_[nca] = nca;

    top_blossoms_.push_back(b);
    bloom(b);

    // set gamma[b_base] = min_slack{critical_edge(b_base, other_base)}
    // where each critical edge is updated before, by
    // argmin{slack(old_bases_in_b, other_base)};
    std::vector<vertex_descriptor_t> b_vertices = b->vertices();
    std::vector<vertex_descriptor_t> old_base_in_b;
    vertex_descriptor_t b_base = b->get_base();
    for (auto i : b_vertices) {
      if (is_old_base[i]) {
        old_base_in_b.push_back(i);
      }
    }
    std::vector<vertex_descriptor_t> other_base;
    for (auto i : vertices(g_)) {
      if (i != b_base && i == base_vertex(i)) {
        other_base.push_back(i);
      }
    }
    for (auto i : other_base) {
      edge_weight_t min_slack = max_edge_weight;
      std::pair<edge_descriptor_t, bool> b_vi = null_edge_;
      for (auto j : old_base_in_b) {
        if (critical_edge_[j][i] != null_edge_ && min_slack > slack(critical_edge_[j][i].first)) {
          min_slack = slack(critical_edge_[j][i].first);
          b_vi = critical_edge_[j][i];
        }
      }
      critical_edge_[b_base][i] = critical_edge_[i][b_base] = b_vi;
    }
    gamma_[b_base] = max_edge_weight;
    for (auto i : other_base) {
      if (critical_edge_[b_base][i] != null_edge_) {
        gamma_[b_base] = std::min(gamma_[b_base], slack(critical_edge_[b_base][i].first));
      }
    }
  }

  void init() {
    even_edges_.clear();
    for (auto u : vertices(g_)) {
      gamma_[u] = tau_[u] = pi_[u] = max_edge_weight;
      critical_edge_[u] = vertex_to_edge_map_t(num_vertices(g_), vm_, null_edge_);

      if (base_vertex(u) != u) {
        continue;
      }

      label_s_[u] = label_t_[u] = graph_traits<Graph>::null_vertex();
      outlet_[u] = u;

      if (mate_[u] == graph_traits<Graph>::null_vertex()) {
        label_s_[u] = u;
        bloom(in_top_blossom(u));
      }
    }
  }

  bool augment_matching() {
    auto v = graph_traits<Graph>::null_vertex();
    auto w = graph_traits<Graph>::null_vertex();
    auto v_free_ancestor = graph_traits<Graph>::null_vertex();
    auto w_free_ancestor = graph_traits<Graph>::null_vertex();
    bool found_alternating_path = false;

    // note that we only use edges of zero slack value for augmenting
    while (!even_edges_.empty() && !found_alternating_path) {
      // search for augmenting paths depth-first
      edge_descriptor_t current_edge = even_edges_.back();
      even_edges_.pop_back();

      v = source(current_edge, g_);
      w = target(current_edge, g_);

      auto v_prime = base_vertex(v);
      auto w_prime = base_vertex(w);

      // w_prime == v_prime implies that we get an edge that has been
      // shrunk into a blossom
      if (v_prime == w_prime) {
        continue;
      }

      // a paranoid check
      if (label_s_[v_prime] == graph_traits<Graph>::null_vertex()) {
        std::swap(v_prime, w_prime);
        std::swap(v, w);
      }

      // w_prime may be unlabeled or have a T-label; replace the existed
      // T-label if the edge slack is smaller than current pi[w_prime] and
      // update it. Note that a T-label is "deserved" only when pi equals
      // zero. also update tau and tau_idx so that tau_idx becomes T-label
      // when a T-blossom is expanded
      if (label_s_[w_prime] == graph_traits<Graph>::null_vertex()) {
        if (slack(current_edge) < pi_[w_prime]) {
          put_t_label(w_prime, v, w, slack(current_edge));
        }
        if (slack(current_edge) < tau_[w]) {
          if (in_blossom_[w]->father == blossom_ptr_t() || label_t_[w_prime] == v ||
              label_t_[w_prime] == graph_traits<Graph>::null_vertex() ||
              nearest_common_ancestor(v_prime, label_t_[w_prime], v_free_ancestor, w_free_ancestor) ==
                  graph_traits<Graph>::null_vertex()) {
            tau_[w] = slack(current_edge);
            tau_idx_[w] = v;
          }
        }
      }

      else {
        if (slack(current_edge) > 0) {
          // update gamma and critical_edges when we have a smaller
          // edge slack
          gamma_[v_prime] = std::min(gamma_[v_prime], slack(current_edge));
          gamma_[w_prime] = std::min(gamma_[w_prime], slack(current_edge));
          if (critical_edge_[v_prime][w_prime] == null_edge_ ||
              slack(critical_edge_[v_prime][w_prime].first) > slack(current_edge)) {
            critical_edge_[v_prime][w_prime] = {current_edge, true};
            critical_edge_[w_prime][v_prime] = {current_edge, true};
          }
          continue;
        }
        if (slack(current_edge) == 0) {
          // if nca is null_vertex then we have an augmenting path;
          // otherwise we have a new top blossom with nca as its base
          // vertex
          vertex_descriptor_t nca = nearest_common_ancestor(v_prime, w_prime, v_free_ancestor, w_free_ancestor);

          if (nca == graph_traits<Graph>::null_vertex()) {
            // to break out of the loop
            found_alternating_path = true;
          } else {
            blossoming(v, v_prime, w, w_prime, nca);
          }
        }
      }
    }

    if (!found_alternating_path) {
      return false;
    }

    augmenting(v, v_free_ancestor, w, w_free_ancestor);
    return true;
  }

  // slack the vertex and blossom dual variables when there is no augmenting
  // path found according to the primal-dual method
  bool adjust_dual() {
    edge_weight_t delta1 = max_edge_weight;
    edge_weight_t delta2 = max_edge_weight;
    edge_weight_t delta3 = max_edge_weight;
    edge_weight_t delta4 = max_edge_weight;

    for (auto v : vertices(g_)) {
      delta1 = std::min(delta1, dual_var_[v]);
      delta4 = pi_[v] > 0 ? std::min(delta4, pi_[v]) : delta4;
      if (v == base_vertex(v)) {
        delta3 = std::min(delta3, gamma_[v] / 2);
      }
    }

    for (const auto& b : top_blossoms_) {
      vertex_descriptor_t b_base = b->get_base();
      if (label_t_[b_base] != graph_traits<Graph>::null_vertex() && pi_[b_base] == 0) {
        delta2 = std::min(delta2, b->dual_var / 2);
      }
    }

    edge_weight_t delta = std::min(std::min(delta1, delta2), std::min(delta3, delta4));

    // start updating dual variables, note that the order is important

    for (auto v : vertices(g_)) {
      auto v_prime = base_vertex(v);

      if (label_s_[v_prime] != graph_traits<Graph>::null_vertex()) {
        dual_var_[v] -= delta;
      } else if (label_t_[v_prime] != graph_traits<Graph>::null_vertex() && pi_[v_prime] == 0) {
        dual_var_[v] += delta;
      }

      if (v == v_prime) {
        gamma_[v] -= 2 * delta;
      }
    }

    for (auto v : vertices(g_)) {
      vertex_descriptor_t v_prime = base_vertex(v);
      if (pi_[v_prime] > 0) {
        tau_[v] -= delta;
      }
    }

    for (const auto& b : top_blossoms_) {
      vertex_descriptor_t b_base = b->get_base();
      if (label_t_[b_base] != graph_traits<Graph>::null_vertex() && pi_[b_base] == 0) {
        b->dual_var -= 2 * delta;
      }
      if (label_s_[b_base] != graph_traits<Graph>::null_vertex()) {
        b->dual_var += 2 * delta;
      }
    }

    for (auto v : vertices(g_)) {
      if (pi_[v] > 0) {
        pi_[v] -= delta;
      }

      // when some T-vertices have zero pi value, bloom their mates so
      // that matching can be further augmented
      if (label_t_[v] != graph_traits<Graph>::null_vertex() && pi_[v] == 0) {
        put_t_label(v, label_t_[v], outlet_[v], pi_[v]);
      }
    }

    // optimal solution reached, halt
    if (delta == delta1) {
      return false;
    }

    // expand odd blossoms with zero dual variables and zero pi value of
    // their base vertices
    if (delta == delta2 && delta != delta3) {
      std::vector<blossom_ptr_t> new_top_blossoms;
      for (blossom_iterator_t bi = top_blossoms_.begin(); bi != top_blossoms_.end();) {
        const blossom_ptr_t b = *bi;
        vertex_descriptor_t b_base = b->get_base();
        if (b->dual_var == 0 && label_t_[b_base] != graph_traits<Graph>::null_vertex() && pi_[b_base] == 0) {
          bi = expand_t_blossom(bi, new_top_blossoms);
        } else {
          ++bi;
        }
      }
      top_blossoms_.insert(top_blossoms_.end(), new_top_blossoms.begin(), new_top_blossoms.end());
    }

    while (true) {
      // find a zero-slack critical edge (v, w) of zero gamma values
      std::pair<edge_descriptor_t, bool> best_edge = null_edge_;
      std::vector<vertex_descriptor_t> base_nodes;
      for (auto v : vertices(g_)) {
        if (v == base_vertex(v)) {
          base_nodes.push_back(v);
        }
      }
      for (auto u : base_nodes) {
        if (gamma_[u] == 0) {
          for (auto v : base_nodes) {
            if (critical_edge_[u][v] != null_edge_ && slack(critical_edge_[u][v].first) == 0) {
              best_edge = critical_edge_[u][v];
            }
          }
        }
      }

      // if not found, continue finding other augment matching
      if (best_edge == null_edge_) {
        return augment_matching() || delta != delta1;
      }
      // if found, determine either augmenting or blossoming
      auto v = source(best_edge.first, g_);
      auto w = target(best_edge.first, g_);
      auto v_prime = base_vertex(v);
      auto w_prime = base_vertex(w);
      auto v_free_ancestor = v_prime;
      auto w_free_ancestor = w_prime;
      vertex_descriptor_t nca = nearest_common_ancestor(v_prime, w_prime, v_free_ancestor, w_free_ancestor);
      if (nca == graph_traits<Graph>::null_vertex()) {
        augmenting(v, v_free_ancestor, w, w_free_ancestor);
        return true;
      }
      blossoming(v, v_prime, w, w_prime, nca);
    }

    return false;
  }

  template <typename PropertyMap>
  void get_current_matching(PropertyMap pm) {
    for (auto v : vertices(g_)) {
      put(pm, v, mate_[v]);
    }
  }

 private:
  const Graph& g_;
  VertexIndexMap vm_;
  WeightMap weight_;
  const std::pair<edge_descriptor_t, bool> null_edge_;

  // iterator property maps
  vertex_to_vertex_map_t mate_;
  // v has an S-label -> v can be an even vertex, label_S[v] is its mate
  vertex_to_vertex_map_t label_s_;
  // v has a T-label -> v can be an odd vertex, label_T[v] is its
  // predecessor in aug_path
  vertex_to_vertex_map_t label_t_;
  vertex_to_vertex_map_t outlet_;
  vertex_to_vertex_map_t tau_idx_;
  vertex_to_weight_map_t dual_var_;
  vertex_to_weight_map_t pi_;
  vertex_to_weight_map_t gamma_;
  vertex_to_weight_map_t tau_;
  // map any vertex v to the trivial blossom containing v
  vertex_to_blossom_map_t in_blossom_;
  // <old T-label, old outlet> before relabeling or expanding T-blossoms
  vertex_to_pair_map_t old_label_;
  // an not matched edge (v, w) is critical if v and w belongs to different S-blossoms
  vertex_pair_to_edge_map_t critical_edge_;

  vertex_list_t aug_path_;
  edge_list_t even_edges_;
  std::vector<blossom_ptr_t> top_blossoms_;
};

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> MateMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap, concepts::ReadableEdgePropertyMap<Graph> WeightMap>
void maximum_weighted_matching(const Graph& g, MateMap mate, VertexIndexMap vm, WeightMap weight) {
  empty_matching(g, mate);
  weighted_augmenting_path_finder augmentor(g, mate, vm, weight);

  // can have |V| times augmenting at most
  for (std::size_t t = 0; t < num_vertices(g); ++t) {
    bool augmented = false;
    while (!augmented) {
      augmented = augmentor.augment_matching();
      if (!augmented) {
        // halt if adjusting dual variables can't bring potential
        // augment
        if (!augmentor.adjust_dual()) {
          break;
        }
      }
    }
    if (!augmented) {
      break;
    }
  }

  augmentor.get_current_matching(mate);
}

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> MateMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
void maximum_weighted_matching(const Graph& g, MateMap mate, VertexIndexMap vm) {
  maximum_weighted_matching(g, mate, vm, get(edge_weight, g));
}

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> MateMap>
void maximum_weighted_matching(const Graph& g, MateMap mate) {
  maximum_weighted_matching(g, mate, get(vertex_index, g));
}

// brute-force matcher searches all possible combinations of matched edges to
// get the maximum weighted matching which can be used for testing on small
// graphs (within dozens vertices)
template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> MateMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap, concepts::ReadableEdgePropertyMap<Graph> WeightMap>
class brute_force_matching {
 public:
  using vertex_descriptor_t = graph_vertex_descriptor_t<Graph>;
  using vertex_to_vertex_map_t = vector_property_map<vertex_descriptor_t, VertexIndexMap>;

  brute_force_matching(const Graph& g, MateMap mate, VertexIndexMap vm, WeightMap weight)
      : g_(g),
        vm_(vm),
        weight_(weight),
        mate_(num_vertices(g_), vm_, graph_traits<Graph>::null_vertex()),
        best_mate_(num_vertices(g_), vm_, graph_traits<Graph>::null_vertex()) {
    for (auto v : vertices(g_)) {
      best_mate_[v] = mate_[v] = get(mate, v);
    }
  }

  template <typename PropertyMap>
  void find_matching(PropertyMap pm) {
    select_edge(partial_view(edges(g_)));

    for (auto v : vertices(g_)) {
      put(pm, v, best_mate_[v]);
    }
  }

 private:
  const Graph& g_;
  VertexIndexMap vm_;
  WeightMap weight_;
  vertex_to_vertex_map_t mate_;
  vertex_to_vertex_map_t best_mate_;

  template <typename PartialEdgeRange>
  void select_edge(PartialEdgeRange e_rg) {
    if (e_rg.empty()) {
      if (matching_weight_sum(g_, mate_, vm_, weight_) > matching_weight_sum(g_, best_mate_, vm_, weight_)) {
        for (auto v : vertices(g_)) {
          best_mate_[v] = mate_[v];
        }
      }
      return;
    }
    auto v = source(*e_rg.begin(), g_);
    auto w = target(*e_rg.begin(), g_);

    e_rg.move_to_next();
    select_edge(e_rg);

    if (mate_[v] == graph_traits<Graph>::null_vertex() && mate_[w] == graph_traits<Graph>::null_vertex()) {
      mate_[v] = w;
      mate_[w] = v;
      select_edge(e_rg);
      mate_[v] = mate_[w] = graph_traits<Graph>::null_vertex();
    }
  }
};

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> MateMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap, concepts::ReadableEdgePropertyMap<Graph> WeightMap>
void brute_force_maximum_weighted_matching(const Graph& g, MateMap mate, VertexIndexMap vm, WeightMap weight) {
  empty_matching(g, mate);
  brute_force_matching brute_force_matcher(g, mate, vm, weight);
  brute_force_matcher.find_matching(mate);
}

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> MateMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
void brute_force_maximum_weighted_matching(const Graph& g, MateMap mate, VertexIndexMap vm) {
  brute_force_maximum_weighted_matching(g, mate, vm, get(edge_weight, g));
}

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> MateMap>
void brute_force_maximum_weighted_matching(const Graph& g, MateMap mate) {
  brute_force_maximum_weighted_matching(g, mate, get(vertex_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_MAXIMUM_WEIGHTED_MATCHING_H_
