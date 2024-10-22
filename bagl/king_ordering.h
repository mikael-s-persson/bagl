// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Copyright 2004, 2005 Trustees of Indiana University
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek,
//          Doug Gregor, D. Kevin McGrath
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_KING_ORDERING_H_
#define BAGL_BAGL_KING_ORDERING_H_

#include <algorithm>
#include <deque>
#include <ranges>
#include <tuple>
#include <vector>

#include "bagl/detail/sparse_ordering.h"
#include "bagl/graph_utility.h"
#include "bagl/two_bit_color_map.h"
#include "bagl/vector_property_map.h"

// King Algorithm for matrix reordering

namespace bagl {
namespace king_detail {
template <typename OutputIterator, typename Buffer, typename Compare, typename PseudoDegreeMap, typename VecMap,
          typename VertexIndexMap>
class bfs_king_visitor : public default_bfs_visitor {
 public:
  bfs_king_visitor(OutputIterator* iter, Buffer* b, Compare compare, PseudoDegreeMap deg, std::vector<int> loc,
                   VecMap color, VertexIndexMap vertices)
      : permutation_(iter),
        q_ptr_(b),
        degree_(deg),
        comp_(compare),
        q_location_(std::move(loc)),
        colors_(color),
        vertex_map_(vertices) {}

  template <typename Vertex, typename Graph>
  void finish_vertex(Vertex /*unused*/, Graph& g) {
    auto rend = q_ptr_->rend() - index_begin_;
    auto rbegin = q_ptr_->rbegin();

    // heap the vertices already there
    std::make_heap(rbegin, rend, [this](Vertex u, Vertex v) { return comp_(v, u); });

    for (std::size_t i = index_begin_; i != q_ptr_->size(); ++i) {
      colors_[get(vertex_map_, (*q_ptr_)[i])] = 1;
      q_location_[get(vertex_map_, (*q_ptr_)[i])] = i;
    }

    for (std::size_t i = 0; rbegin != rend; rend--) {
      percolate_down<Vertex>(i);
      auto w = (*q_ptr_)[index_begin_ + i];
      for (auto e : out_edges(w, g)) {
        auto v = target(e, g);
        put(degree_, v, get(degree_, v) - 1);

        if (colors_[get(vertex_map_, v)] == 1) {
          percolate_up<Vertex>(get(vertex_map_, v), i);
        }
      }

      colors_[get(vertex_map_, w)] = 0;
      i++;
    }
  }

  template <typename Vertex, typename Graph>
  void examine_vertex(Vertex u, const Graph& /*unused*/) {
    *(*permutation_)++ = u;
    index_begin_ = q_ptr_->size();
  }

 protected:
  // this function replaces pop_heap, and tracks state information
  template <typename Vertex>
  void percolate_down(int offset) {
    int heap_last = index_begin_ + offset;
    int heap_first = q_ptr_->size() - 1;

    // pop_heap functionality:
    // swap first, last
    std::swap((*q_ptr_)[heap_last], (*q_ptr_)[heap_first]);

    // swap in the location queue
    std::swap(q_location_[heap_first], q_location_[heap_last]);

    // set drifter, children
    int drifter = heap_first;
    int drifter_heap = q_ptr_->size() - drifter;

    int right_child_heap = drifter_heap * 2 + 1;
    int right_child = q_ptr_->size() - right_child_heap;

    int left_child_heap = drifter_heap * 2;
    int left_child = q_ptr_->size() - left_child_heap;

    // check that we are staying in the heap
    bool valid = right_child >= heap_last;

    // pick smallest child of drifter, and keep in mind there might only
    // be left child
    int smallest_child = (valid && get(degree_, (*q_ptr_)[left_child]) > get(degree_, (*q_ptr_)[right_child]))
                             ? right_child
                             : left_child;

    while (valid && smallest_child < heap_last && comp_((*q_ptr_)[drifter], (*q_ptr_)[smallest_child])) {
      // if smallest child smaller than drifter, swap them
      std::swap((*q_ptr_)[smallest_child], (*q_ptr_)[drifter]);
      std::swap(q_location_[drifter], q_location_[smallest_child]);

      // update the values, run again, as necessary
      drifter = smallest_child;
      drifter_heap = q_ptr_->size() - drifter;

      right_child_heap = drifter_heap * 2 + 1;
      right_child = q_ptr_->size() - right_child_heap;

      left_child_heap = drifter_heap * 2;
      left_child = q_ptr_->size() - left_child_heap;

      valid = right_child >= heap_last;

      smallest_child = (valid && get(degree_, (*q_ptr_)[left_child]) > get(degree_, (*q_ptr_)[right_child]))
                           ? right_child
                           : left_child;
    }
  }

  // this is like percolate down, but we always compare against the
  // parent, as there is only a single choice
  template <typename Vertex>
  void percolate_up(int vertex, int offset) {
    int child_location = q_location_[vertex];
    int heap_child_location = q_ptr_->size() - child_location;
    int heap_parent_location = (heap_child_location / 2);
    int parent_location = q_ptr_->size() - heap_parent_location;

    bool valid =
        (heap_parent_location != 0 && child_location > index_begin_ + offset && parent_location < q_ptr_->size());

    while (valid && comp_((*q_ptr_)[child_location], (*q_ptr_)[parent_location])) {
      // swap in the heap
      std::swap((*q_ptr_)[child_location], (*q_ptr_)[parent_location]);

      // swap in the location queue
      std::swap(q_location_[child_location], q_location_[parent_location]);

      child_location = parent_location;
      heap_child_location = heap_parent_location;
      heap_parent_location = (heap_child_location / 2);
      parent_location = q_ptr_->size() - heap_parent_location;
      valid = (heap_parent_location != 0 && child_location > index_begin_ + offset);
    }
  }

  OutputIterator* permutation_ = nullptr;
  int index_begin_ = 0;
  Buffer* q_ptr_ = nullptr;
  PseudoDegreeMap degree_;
  Compare comp_;
  std::vector<int> q_location_;
  VecMap colors_;
  VertexIndexMap vertex_map_;
};

}  // namespace king_detail

template <concepts::VertexListGraph G, std::ranges::input_range VertexRange,
          std::output_iterator<graph_vertex_descriptor_t<G>> OutputIter,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, concepts::ReadableVertexPropertyMap<G> DegreeMap,
          concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
OutputIter king_ordering(const G& g, const VertexRange& vertex_rg, OutputIter permutation, ColorMap color,
                         DegreeMap degree, VertexIndexMap index_map) {
  using DegreeType = property_traits_value_t<DegreeMap>;
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;
  using Vertex = graph_vertex_descriptor_t<G>;

  // Copy degree to pseudo_degree
  // initialize the color map
  auto pseudo_degree = vector_property_map(num_vertices(g), index_map, DegreeType{});
  for (auto u : vertices(g)) {
    put(pseudo_degree, u, get(degree, u));
    put(color, u, Color::white());
  }

  auto comp = [&pseudo_degree](auto u, auto v) { return get(pseudo_degree, u) < get(pseudo_degree, v); };
  std::vector<int> colors(num_vertices(g));

  for (std::size_t i = 0; i < num_vertices(g); i++) {
    colors[i] = 0;
  }

  std::vector<int> loc(num_vertices(g));

  // create the visitor
  sparse_detail::sparse_ordering_queue<Vertex> q;
  king_detail::bfs_king_visitor vis{&permutation, &q, comp, pseudo_degree.ref(), loc, colors, index_map};

  for (Vertex s : vertex_rg) {
    // call BFS with visitor
    breadth_first_visit(g, s, q, vis, color);
  }

  return permutation;
}

// This is the case where only a single starting vertex is supplied.
template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIter,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, concepts::ReadableVertexPropertyMap<G> DegreeMap,
          concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
OutputIter king_ordering(const G& g, graph_vertex_descriptor_t<G> start, OutputIter permutation, ColorMap color,
                         DegreeMap degree, VertexIndexMap index_map) {
  return king_ordering(g, std::ranges::single_view{start}, permutation, color, degree, index_map);
}

template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIter,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, concepts::ReadableVertexPropertyMap<G> DegreeMap,
          concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
OutputIter king_ordering(const G& g, OutputIter permutation, ColorMap color, DegreeMap degree,
                         VertexIndexMap index_map) {
  if (has_no_vertices(g)) {
    return permutation;
  }

  using Vertex = graph_vertex_descriptor_t<G>;
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;

  // Mark everything white
  for (auto v : vertices(g)) {
    put(color, v, Color::white());
  }

  // Find one vertex from each connected component
  std::vector<Vertex> vertex_queue;
  for (auto v : vertices(g)) {
    if (get(color, v) != Color::white()) {
      continue;
    }
    depth_first_visit(g, v, dfs_visitor<>(), color);
    vertex_queue.push_back(v);
  }

  // Find starting nodes for all vertices
  // TBD: How to do this with a directed graph?
  for (auto& v : vertex_queue) {
    v = find_starting_node(g, v, color, degree);
  }

  return king_ordering(g, vertex_queue, permutation, color, degree, index_map);
}

template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIter,
          concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
OutputIter king_ordering(const G& g, OutputIter permutation, VertexIndexMap index_map) {
  if (has_no_vertices(g)) {
    return permutation;
  }

  return king_ordering(g, permutation, two_bit_color_map(num_vertices(g), index_map).ref(), make_out_degree_map(g),
                       index_map);
}

template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIter>
OutputIter king_ordering(const G& g, OutputIter permutation) {
  return king_ordering(g, permutation, get(vertex_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_KING_ORDERING_H_
