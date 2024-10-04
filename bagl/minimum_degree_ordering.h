// Copyright 1997-2001 University of Notre Dame.
// Authors: Lie-Quan Lee, Jeremy Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_MINIMUM_DEGREE_ORDERING_H_
#define BAGL_BAGL_MINIMUM_DEGREE_ORDERING_H_

#include <cassert>
#include <cstddef>
#include <limits>
#include <vector>

#include "bagl/bucket_sorter.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"

namespace bagl {

namespace ordering_detail {

//
// Given a set of n integers (where the integer values range from
// zero to n-1), we want to keep track of a collection of stacks
// of integers. It so happens that an integer will appear in at
// most one stack at a time, so the stacks form disjoint sets.
// Because of these restrictions, we can use one big array to
// store all the stacks, intertwined with one another.
// No allocation/deallocation happens in the push()/pop() methods
// so this is faster than using std::stack's.
//
class stacks {
  using value_type = std::ptrdiff_t;
  using size_type = std::size_t;
  static constexpr value_type invalid_value = -std::numeric_limits<value_type>::max();

 public:
  explicit stacks(size_type n) : data_(n, invalid_value) {}

  //: stack
  class stack {
    using Iterator = typename std::vector<value_type>::iterator;

   public:
    explicit stack(Iterator data, value_type head = invalid_value) : data_(data), current_(head) {}

    void pop() {
      assert(!empty());
      current_ = data_[current_];
    }
    void push(value_type v) {
      data_[v] = current_;
      current_ = v;
    }
    [[nodiscard]] bool empty() const { return current_ == invalid_value; }
    [[nodiscard]] value_type& top() { return current_; }

   private:
    Iterator data_;
    value_type current_;
  };

  // To return a stack object
  [[nodiscard]] auto make_stack() { return stack(data_.begin()); }

 private:
  std::vector<value_type> data_;
};

// marker class, a generalization of coloring.
//
// This class is to provide a generalization of coloring which has
// complexity of amortized constant time to set all vertices' color
// back to be untagged. It implemented by increasing a tag.
//
// The colors are:
//   not tagged
//   tagged
//   multiple_tagged
//   done
//
template <class Vertex, class VertexIndexMap>
class marker {
  using value_type = std::ptrdiff_t;
  using size_type = std::size_t;

  static constexpr value_type invalid_value = -std::numeric_limits<value_type>::max();
  static constexpr value_type base_tag_value = value_type{1} - std::numeric_limits<value_type>::max();
  static constexpr value_type done_value = std::numeric_limits<value_type>::max() / 2;

 public:
  marker(size_type num, VertexIndexMap index_map) : tag_(base_tag_value), data_(num, invalid_value), id_(index_map) {}

  void mark_done(Vertex node) { data_[get(id_, node)] = done_value; }

  [[nodiscard]] bool is_done(Vertex node) const { return data_[get(id_, node)] == done_value; }

  void mark_tagged(Vertex node) { data_[get(id_, node)] = tag_; }

  void mark_multiple_tagged(Vertex node) { data_[get(id_, node)] = multiple_tag_; }

  [[nodiscard]] bool is_tagged(Vertex node) const { return data_[get(id_, node)] >= tag_; }

  [[nodiscard]] bool is_not_tagged(Vertex node) const { return data_[get(id_, node)] < tag_; }

  [[nodiscard]] bool is_multiple_tagged(Vertex node) const { return data_[get(id_, node)] >= multiple_tag_; }

  void increment_tag() {
    const size_type num = data_.size();
    ++tag_;
    if (tag_ >= done_value) {
      tag_ = base_tag_value;
      for (size_type i = 0; i < num; ++i) {
        if (data_[i] < done_value) {
          data_[i] = invalid_value;
        }
      }
    }
  }

  void set_multiple_tag(value_type mdeg0) {
    const size_type num = data_.size();
    multiple_tag_ = tag_ + mdeg0;

    if (multiple_tag_ >= done_value) {
      tag_ = base_tag_value;

      for (size_type i = 0; i < num; i++) {
        if (data_[i] < done_value) {
          data_[i] = invalid_value;
        }
      }

      multiple_tag_ = tag_ + mdeg0;
    }
  }

  void set_tag_as_multiple_tag() { tag_ = multiple_tag_; }

 private:
  value_type tag_ = base_tag_value;
  value_type multiple_tag_ = invalid_value;
  std::vector<value_type> data_;
  VertexIndexMap id_;
};

template <concepts::MutableLvaluePropertyMap<std::ptrdiff_t> Iterator, class Vertex, class VertexIndexMap,
          int offset = 1>
class numbering {
  using number_type = std::ptrdiff_t;

 public:
  numbering(Iterator data, number_type max_num, VertexIndexMap id) : data_(data), max_num_(max_num), id_(id) {}
  void operator()(Vertex node) { data_[get(id_, node)] = -num_; }
  [[nodiscard]] bool all_done(number_type i = 0) const { return num_ + i > max_num_; }
  void increment(number_type i = 1) { num_ += i; }
  [[nodiscard]] bool is_numbered(Vertex node) const { return data_[get(id_, node)] < 0; }
  void indistinguishable(Vertex i, Vertex j) { data_[get(id_, i)] = -(get(id_, j) + offset); }

 private:
  number_type num_ = 1;  // start from 1 instead of zero
  Iterator data_;
  number_type max_num_;
  VertexIndexMap id_;
};

template <class Vertex, class VertexIndexMap>
class degreelists_marker {
 public:
  using value_type = std::ptrdiff_t;
  using size_type = std::size_t;
  degreelists_marker(size_type n, VertexIndexMap id) : marks_(n, 0), id_(id) {}
  void mark_need_update(Vertex i) { marks_[get(id_, i)] = 1; }
  [[nodiscard]] bool need_update(Vertex i) const { return marks_[get(id_, i)] == 1; }
  [[nodiscard]] bool outmatched_or_done(Vertex i) const { return marks_[get(id_, i)] == -1; }
  void mark(Vertex i) { marks_[get(id_, i)] = -1; }
  void unmark(Vertex i) { marks_[get(id_, i)] = 0; }

 private:
  std::vector<value_type> marks_;
  VertexIndexMap id_;
};

// Helper function object for edge removal
template <class Graph, class MarkerP, class NumberD, class Stack, class VertexIndexMap>
class predicate_remove_edge {
  using vertex_t = graph_vertex_descriptor_t<Graph>;
  using edge_t = graph_edge_descriptor_t<Graph>;

 public:
  predicate_remove_edge(Graph& g, MarkerP& marker, NumberD numbering, Stack& n_e, VertexIndexMap id)
      : g_(&g), marker_(&marker), numbering_(numbering), neighbor_elements_(&n_e), id_(id) {}

  bool operator()(edge_t e) {
    vertex_t dist = target(e, *g_);
    if (marker_->is_tagged(dist)) {
      return true;
    }
    marker_->mark_tagged(dist);
    if (numbering_.is_numbered(dist)) {
      neighbor_elements_->push(get(id_, dist));
      return true;
    }
    return false;
  }

 private:
  Graph* g_;
  MarkerP* marker_;
  NumberD numbering_;
  Stack* neighbor_elements_;
  VertexIndexMap id_;
};

// Helper function object for edge removal
template <class Graph, class MarkerP>
class predicate_remove_tagged_edges {
  using vertex_t = graph_vertex_descriptor_t<Graph>;
  using edge_t = graph_edge_descriptor_t<Graph>;

 public:
  predicate_remove_tagged_edges(Graph& g, MarkerP& marker) : g_(&g), marker_(&marker) {}

  bool operator()(edge_t e) {
    vertex_t dist = target(e, *g_);
    return marker_->is_tagged(dist);
  }

 private:
  Graph* g_;
  MarkerP* marker_;
};

template <concepts::EdgeMutableGraph Graph, concepts::MutableLvalueVertexPropertyMap<Graph> DegreeMap,
          concepts::MutableLvaluePropertyMap<std::ptrdiff_t> InversePermutationMap,
          concepts::MutableLvalueVertexPropertyMap<Graph> SuperNodeMap, class VertexIndexMap>
requires concepts::VertexListGraph<Graph> && concepts::IncidenceGraph<Graph> && concepts::AdjacencyGraph<Graph>
class mmd_impl {
  // Typedefs
  using size_type = std::size_t;
  using diff_t = std::ptrdiff_t;
  using vertex_t = graph_vertex_descriptor_t<Graph>;
  using IndexVertexMap = iterator_property_map<const vertex_t*, identity_property_map>;
  using Workspace = stacks;
  using DegreeLists = bucket_sorter<size_type, vertex_t, DegreeMap, VertexIndexMap>;
  using NumberingD = numbering<InversePermutationMap, vertex_t, VertexIndexMap>;
  using DegreeListsMarker = degreelists_marker<vertex_t, VertexIndexMap>;
  using MarkerP = marker<vertex_t, VertexIndexMap>;

  // Data Members
  bool has_no_edges_ = true;

  // input parameters
  Graph& g_;
  int delta_;
  DegreeMap degree_;
  InversePermutationMap inverse_perm_;
  SuperNodeMap supernode_size_;
  VertexIndexMap vertex_index_map_;

  // internal data-structures
  std::vector<vertex_t> index_vertex_vec_;
  size_type n_;
  IndexVertexMap index_vertex_map_;
  DegreeLists degreelists_;
  NumberingD numbering_;
  DegreeListsMarker degree_lists_marker_;
  MarkerP marker_;
  Workspace work_space_;

 public:
  mmd_impl(Graph& g, size_type n, int delta, DegreeMap degree, InversePermutationMap inverse_perm,
           SuperNodeMap supernode_size, VertexIndexMap id)
      : g_(g),
        delta_(delta),
        degree_(degree),
        inverse_perm_(inverse_perm),
        supernode_size_(supernode_size),
        vertex_index_map_(id),
        index_vertex_vec_(n),
        n_(n),
        index_vertex_map_(IndexVertexMap(index_vertex_vec_.data())),
        degreelists_(n_ + 1, n_, degree, id),
        numbering_(inverse_perm, n_, vertex_index_map_),
        degree_lists_marker_(n_, vertex_index_map_),
        marker_(n_, vertex_index_map_),
        work_space_(n_) {
    size_type vid = 0;
    for (auto v : vertices(g_)) {
      index_vertex_vec_[vid] = v;
      ++vid;
    }

    // Initialize degreelists.  Degreelists organizes the nodes
    // according to their degree.
    for (auto v : vertices(g_)) {
      std::size_t d = out_degree(v, g_);
      put(degree, v, d);
      if (0 < d) {
        has_no_edges_ = false;
      }
      degreelists_.push(v);
    }
  }

  void do_mmd() {
    // Eliminate the isolated nodes -- these are simply the nodes
    // with no neighbors, which are accessible as a list (really, a
    // stack) at location 0.  Since these don't affect any other
    // nodes, we can eliminate them without doing degree updates.
    auto list_isolated = degreelists_[0];
    while (!list_isolated.empty()) {
      vertex_t node = list_isolated.top();
      marker_.mark_done(node);
      numbering_(node);
      numbering_.increment();
      list_isolated.pop();
    }

    if (has_no_edges_) {
      return;
    }

    size_type min_degree = 1;
    auto list_min_degree = degreelists_[min_degree];
    while (list_min_degree.empty()) {
      ++min_degree;
      list_min_degree = degreelists_[min_degree];
    }

    // check if the whole eliminating process is done
    while (!numbering_.all_done()) {
      size_type min_degree_limit = min_degree + delta_;  // WARNING
      auto llist = work_space_.make_stack();

      // multiple elimination
      while (delta_ >= 0) {
        // Find the next non-empty degree
        list_min_degree = degreelists_[min_degree];
        while (list_min_degree.empty() && min_degree <= min_degree_limit) {
          list_min_degree = degreelists_[++min_degree];
        }
        if (min_degree > min_degree_limit) {
          break;
        }

        const vertex_t node = list_min_degree.top();
        const size_type node_id = get(vertex_index_map_, node);
        list_min_degree.pop();
        numbering_(node);

        // check if node is the last one
        if (numbering_.all_done(supernode_size_[node])) {
          numbering_.increment(supernode_size_[node]);
          break;
        }
        marker_.increment_tag();
        marker_.mark_tagged(node);

        eliminate(node);

        numbering_.increment(supernode_size_[node]);
        llist.push(node_id);
      }

      if (numbering_.all_done()) {
        break;
      }

      update(llist, min_degree);
    }
  }

  void eliminate(vertex_t node) {
    auto element_neighbor = work_space_.make_stack();

    // Create two function objects for edge removal
    auto p = predicate_remove_edge(g_, marker_, numbering_, element_neighbor, vertex_index_map_);
    auto p2 = predicate_remove_tagged_edges(g_, marker_);

    // Reconstruct the adjacent node list, push element neighbor in a
    // List.
    remove_out_edge_if(node, p, g_);
    // during removal element neighbors are collected.

    while (!element_neighbor.empty()) {
      // element absorb
      size_type e_id = element_neighbor.top();
      vertex_t element = get(index_vertex_map_, e_id);
      for (auto i_node : adjacent_vertices(element, g_)) {
        if (!marker_.is_tagged(i_node) && !numbering_.is_numbered(i_node)) {
          marker_.mark_tagged(i_node);
          add_edge(node, i_node, g_);
        }
      }
      element_neighbor.pop();
    }
    for (auto v_node : adjacent_vertices(node, g_)) {
      if (!degree_lists_marker_.need_update(v_node) && !degree_lists_marker_.outmatched_or_done(v_node)) {
        degreelists_.remove(v_node);
      }
      // update out edges of v_node
      remove_out_edge_if(v_node, p2, g_);

      if (out_degree(v_node, g_) == 0) {  // indistinguishable nodes
        supernode_size_[node] += supernode_size_[v_node];
        supernode_size_[v_node] = 0;
        numbering_.indistinguishable(v_node, node);
        marker_.mark_done(v_node);
        degree_lists_marker_.mark(v_node);
      } else {  // not indistinguishable nodes
        add_edge(v_node, node, g_);
        degree_lists_marker_.mark_need_update(v_node);
      }
    }
  }  // eliminate()

  template <class Stack>
  void update(Stack llist, size_type& min_degree) {
    size_type min_degree0 = min_degree + delta_ + 1;

    while (!llist.empty()) {
      size_type deg = 0;
      size_type deg0 = 0;

      marker_.set_multiple_tag(min_degree0);
      auto q2list = work_space_.make_stack();
      auto qxlist = work_space_.make_stack();

      vertex_t current = get(index_vertex_map_, llist.top());
      for (auto i_node : adjacent_vertices(current, g_)) {
        const size_type i_id = get(vertex_index_map_, i_node);
        if (supernode_size_[i_node] != 0) {
          deg0 += supernode_size_[i_node];
          marker_.mark_multiple_tagged(i_node);
          if (degree_lists_marker_.need_update(i_node)) {
            if (out_degree(i_node, g_) == 2) {
              q2list.push(i_id);
            } else {
              qxlist.push(i_id);
            }
          }
        }
      }

      while (!q2list.empty()) {
        const size_type u_id = q2list.top();
        vertex_t u_node = get(index_vertex_map_, u_id);
        // if u_id is outmatched by others, no need to update degree
        if (degree_lists_marker_.outmatched_or_done(u_node)) {
          q2list.pop();
          continue;
        }
        marker_.increment_tag();
        deg = deg0;

        auto nu_rg = adjacent_vertices(u_node, g_);
        auto nu = nu_rg.begin();
        vertex_t neighbor = *nu;
        if (neighbor == u_node) {
          ++nu;
          neighbor = *nu;
        }
        if (numbering_.is_numbered(neighbor)) {
          for (auto i_node : adjacent_vertices(neighbor, g_)) {
            if (i_node == u_node || supernode_size_[i_node] == 0) {
              continue;
            }
            if (marker_.is_tagged(i_node)) {
              if (degree_lists_marker_.need_update(i_node)) {
                if (out_degree(i_node, g_) == 2) {  // is indistinguishable
                  supernode_size_[u_node] += supernode_size_[i_node];
                  supernode_size_[i_node] = 0;
                  numbering_.indistinguishable(i_node, u_node);
                  marker_.mark_done(i_node);
                  degree_lists_marker_.mark(i_node);
                } else {  // is outmatched
                  degree_lists_marker_.mark(i_node);
                }
              }
            } else {
              marker_.mark_tagged(i_node);
              deg += supernode_size_[i_node];
            }
          }
        } else {
          deg += supernode_size_[neighbor];
        }

        deg -= supernode_size_[u_node];
        degree_[u_node] = deg;  // update degree
        degreelists_[deg].push(u_node);
        // u_id has been pushed back into degreelists
        degree_lists_marker_.unmark(u_node);
        if (min_degree > deg) {
          min_degree = deg;
        }
        q2list.pop();
      }  // while (!q2list.empty())

      while (!qxlist.empty()) {
        const size_type u_id = qxlist.top();
        const vertex_t u_node = get(index_vertex_map_, u_id);

        // if u_id is outmatched by others, no need to update degree
        if (degree_lists_marker_.outmatched_or_done(u_node)) {
          qxlist.pop();
          continue;
        }
        marker_.increment_tag();
        deg = deg0;
        for (auto i_node : adjacent_vertices(u_node, g_)) {
          if (marker_.is_tagged(i_node)) {
            continue;
          }
          marker_.mark_tagged(i_node);

          if (numbering_.is_numbered(i_node)) {
            for (const vertex_t j_node : adjacent_vertices(i_node, g_)) {
              if (marker_.is_not_tagged(j_node)) {
                marker_.mark_tagged(j_node);
                deg += supernode_size_[j_node];
              }
            }
          } else {
            deg += supernode_size_[i_node];
          }
        }  // for adjacent vertices of u_node
        deg -= supernode_size_[u_node];
        degree_[u_node] = deg;
        degreelists_[deg].push(u_node);
        // u_id has been pushed back into degreelists
        degree_lists_marker_.unmark(u_node);
        if (min_degree > deg) {
          min_degree = deg;
        }
        qxlist.pop();
      }  // while (!qxlist.empty()) {

      marker_.set_tag_as_multiple_tag();
      llist.pop();
    }  // while (! llist.empty())

  }  // update()

  template <concepts::MutableLvaluePropertyMap<diff_t> PermutationMap>
  void build_permutation(InversePermutationMap next, PermutationMap prev) {
    // collect the permutation info
    for (diff_t i = 0; i < n_; ++i) {
      diff_t size = supernode_size_[get(index_vertex_map_, i)];
      if (size <= 0) {
        prev[i] = next[i];
        supernode_size_[get(index_vertex_map_, i)] = next[i] + 1;  // record the supernode info
      } else {
        prev[i] = -next[i];
      }
    }
    for (diff_t i = 1; i < n_ + 1; ++i) {
      if (prev[i - 1] > 0) {
        continue;
      }
      diff_t parent = i;
      while (prev[parent - 1] < 0) {
        parent = -prev[parent - 1];
      }

      diff_t root = parent;
      diff_t num = prev[root - 1] + 1;
      next[i - 1] = -num;
      prev[root - 1] = num;

      parent = i;
      diff_t next_node = -prev[parent - 1];
      while (next_node > 0) {
        prev[parent - 1] = -root;
        parent = next_node;
        next_node = -prev[parent - 1];
      }
    }
    for (diff_t i = 0; i < n_; i++) {
      diff_t num = -next[i] - 1;
      next[i] = num;
      prev[num] = i;
    }
  }  // build_permutation()
};

}  // namespace ordering_detail

// MMD algorithm
//
// The implementation presently includes the enhancements for mass
// elimination, incomplete degree update, multiple elimination, and
// external degree.
//
// Important Note: This implementation requires the BGL graph to be
// directed.  Therefore, nonzero entry (i, j) in a symmetrical matrix
// A coresponds to two directed edges (i->j and j->i).
//
// see Alan George and Joseph W. H. Liu, The Evolution of the Minimum
// Degree Ordering Algorithm, SIAM Review, 31, 1989, Page 1-19
template <concepts::EdgeMutableGraph Graph, concepts::MutableLvalueVertexPropertyMap<Graph> DegreeMap,
          concepts::MutableLvaluePropertyMap<std::ptrdiff_t> InversePermutationMap,
          concepts::MutableLvaluePropertyMap<std::ptrdiff_t> PermutationMap,
          concepts::MutableLvalueVertexPropertyMap<Graph> SuperNodeMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
requires concepts::VertexListGraph<Graph> && concepts::IncidenceGraph<Graph> && concepts::AdjacencyGraph<Graph>
void minimum_degree_ordering(Graph& g, DegreeMap degree, InversePermutationMap inverse_perm, PermutationMap perm,
                             SuperNodeMap supernode_size, int delta, VertexIndexMap vertex_index_map) {
  ordering_detail::mmd_impl impl(g, num_vertices(g), delta, degree, inverse_perm, supernode_size, vertex_index_map);
  impl.do_mmd();
  impl.build_permutation(inverse_perm, perm);
}

}  // namespace bagl

#endif  // BAGL_BAGL_MINIMUM_DEGREE_ORDERING_H_
