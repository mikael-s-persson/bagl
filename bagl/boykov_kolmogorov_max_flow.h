// Copyright 2006 Stephan Diederich
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_BOYKOV_KOLMOGOROV_MAX_FLOW_H_
#define BAGL_BAGL_BOYKOV_KOLMOGOROV_MAX_FLOW_H_

#include <algorithm>  // for std::min and std::max
#include <cassert>
#include <iosfwd>
#include <limits>
#include <list>
#include <memory>
#include <queue>
#include <utility>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/lookup_edge.h"
#include "bagl/partial_range.h"
#include "bagl/property_map.h"

// The algorithm impelemented here is described in:
//
// Boykov, Y., Kolmogorov, V. "An Experimental Comparison of Min-Cut/Max-Flow
// Algorithms for Energy Minimization in Vision", In IEEE Transactions on
// Pattern Analysis and Machine Intelligence, vol. 26, no. 9, pp. 1124-1137,
// Sep 2004.
//
// For further reading, also see:
//
// Kolmogorov, V. "Graph Based Algorithms for Scene Reconstruction from Two or
// More Views". PhD thesis, Cornell University, Sep 2003.

namespace bagl {

namespace max_flow_detail {

template <class Graph, class EdgeCapacityMap, class ResidualCapacityEdgeMap, class ReverseEdgeMap, class PredecessorMap,
          class ColorMap, class DistanceMap, class IndexMap>
class bk_max_flow {
  using edge_flow_value = property_traits_value_t<EdgeCapacityMap>;
  using vertex_descriptor = graph_vertex_descriptor_t<Graph>;
  using edge_descriptor = graph_edge_descriptor_t<Graph>;
  using out_edge_range = graph_out_edge_range_t<Graph>;
  using vertex_queue = std::queue<vertex_descriptor>;  // queue of vertices, used in adoption-stage
  using color_value = property_traits_value_t<ColorMap>;
  using color = color_traits<color_value>;
  using distance_value = property_traits_value_t<DistanceMap>;

 public:
  bk_max_flow(Graph& g, EdgeCapacityMap cap, ResidualCapacityEdgeMap res, ReverseEdgeMap rev, PredecessorMap pre,
              ColorMap color, DistanceMap dist, IndexMap idx, vertex_descriptor src, vertex_descriptor sink)
      : g_(g),
        index_map_(idx),
        cap_map_(cap),
        res_cap_map_(res),
        rev_edge_map_(rev),
        pre_map_(pre),
        tree_map_(color),
        dist_map_(dist),
        source_(src),
        sink_(sink),
        in_active_list_map_(num_vertices(g), index_map_, false),
        has_parent_map_(num_vertices(g), index_map_, false),
        time_map_(num_vertices(g), index_map_, 0) {
    // initialize the color-map with gray-values
    for (auto v : vertices(g_)) {
      set_tree(v, color::gray());
    }
    // Initialize flow to zero which means initializing
    // the residual capacity equal to the capacity
    for (auto e : edges(g_)) {
      put(res_cap_map_, e, get(cap_map_, e));
      // check if the reverse edge map is build up properly
      assert(get(rev_edge_map_, get(rev_edge_map_, e)) == e);
    }
    // init the search trees with the two terminals
    set_tree(source_, color::black());
    set_tree(sink_, color::white());
    put(time_map_, source_, 1);
    put(time_map_, sink_, 1);
  }

  edge_flow_value max_flow() {
    // augment direct paths from SOURCE->SINK and SOURCE->VERTEX->SINK
    augment_direct_paths();
    // start the main-loop
    while (true) {
      auto [connecting_edge, path_found] = grow();  // find a path from source to sink
      if (!path_found) {
        // we're finished, no more paths were found
        break;
      }
      ++time_;
      augment(connecting_edge);  // augment that path
      adopt();                   // rebuild search tree structure
    }
    return flow_;
  }

  // the complete class is protected, as we want access to members in
  // derived test-class (see test/boykov_kolmogorov_max_flow_test.cpp)
 protected:
  void augment_direct_paths() {
    // in a first step, we augment all direct paths from
    // source->NODE->sink and additionally paths from source->sink. This
    // improves especially graphcuts for segmentation, as most of the
    // nodes have source/sink connects but shouldn't have an impact on
    // other maxflow problems (this is done in grow() anyway)
    for (auto from_source : out_edges(source_, g_)) {
      auto current_node = target(from_source, g_);
      if (current_node == sink_) {
        edge_flow_value cap = get(res_cap_map_, from_source);
        put(res_cap_map_, from_source, 0);
        flow_ += cap;
        continue;
      }
      auto [to_sink, is_there] = lookup_edge(current_node, sink_, g_);
      if (is_there) {
        edge_flow_value cap_from_source = get(res_cap_map_, from_source);
        edge_flow_value cap_to_sink = get(res_cap_map_, to_sink);
        if (cap_from_source > cap_to_sink) {
          set_tree(current_node, color::black());
          add_active_node(current_node);
          set_edge_to_parent(current_node, from_source);
          put(dist_map_, current_node, 1);
          put(time_map_, current_node, 1);
          // add stuff to flow and update residuals. we dont need
          // to update reverse_edges, as incoming/outgoing edges
          // to/from source/sink don't count for max-flow
          put(res_cap_map_, from_source, get(res_cap_map_, from_source) - cap_to_sink);
          put(res_cap_map_, to_sink, 0);
          flow_ += cap_to_sink;
        } else if (cap_to_sink > 0) {
          set_tree(current_node, color::white());
          add_active_node(current_node);
          set_edge_to_parent(current_node, to_sink);
          put(dist_map_, current_node, 1);
          put(time_map_, current_node, 1);
          // add stuff to flow and update residuals. we dont need
          // to update reverse_edges, as incoming/outgoing edges
          // to/from source/sink don't count for max-flow
          put(res_cap_map_, to_sink, get(res_cap_map_, to_sink) - cap_from_source);
          put(res_cap_map_, from_source, 0);
          flow_ += cap_from_source;
        }
      } else if (get(res_cap_map_, from_source)) {
        // there is no sink connect, so we can't augment this path,
        // but to avoid adding source_ to the active nodes, we just
        // activate this node and set the approciate things
        set_tree(current_node, color::black());
        set_edge_to_parent(current_node, from_source);
        put(dist_map_, current_node, 1);
        put(time_map_, current_node, 1);
        add_active_node(current_node);
      }
    }
    for (auto e : out_edges(sink_, g_)) {
      edge_descriptor to_sink = get(rev_edge_map_, e);
      vertex_descriptor current_node = source(to_sink, g_);
      if (get(res_cap_map_, to_sink)) {
        set_tree(current_node, color::white());
        set_edge_to_parent(current_node, to_sink);
        put(dist_map_, current_node, 1);
        put(time_map_, current_node, 1);
        add_active_node(current_node);
      }
    }
  }

  /**
   * Returns a pair of an edge and a boolean. if the bool is true, the
   * edge is a connection of a found path from s->t , read "the link" and
   * source(returnVal, g_) is the end of the path found in the
   * source-tree target(returnVal, g_) is the beginning of the path found
   * in the sink-tree
   */
  std::pair<edge_descriptor, bool> grow() {
    assert(orphans_.empty());
    vertex_descriptor current_node;
    while ((current_node = get_next_active_node()) != graph_traits<Graph>::null_vertex()) {  // if there is one
      assert(get_tree(current_node) != color::gray() &&
             (has_parent(current_node) || current_node == source_ || current_node == sink_));

      if (get_tree(current_node) == color::black()) {
        // source tree growing
        if (current_node != last_grow_vertex_) {
          last_grow_vertex_ = current_node;
          last_grow_edge_range_ = std::make_unique<partial_view<out_edge_range>>(out_edges(current_node, g_));
        }
        for (; !last_grow_edge_range_->empty(); last_grow_edge_range_->move_to_next()) {
          edge_descriptor out_edge = *last_grow_edge_range_->begin();
          if (get(res_cap_map_, out_edge) > 0) {  // check if we have capacity left on this edge
            vertex_descriptor other_node = target(out_edge, g_);
            if (get_tree(other_node) == color::gray()) {  // it's a free node
              set_tree(other_node,
                       color::black());                                      // aquire other node
                                                                             // to our search
                                                                             // tree
              set_edge_to_parent(other_node, out_edge);                      // set us as parent
              put(dist_map_, other_node, get(dist_map_, current_node) + 1);  // and update the
                                                                             // distance-heuristic
              put(time_map_, other_node, get(time_map_, current_node));
              add_active_node(other_node);
            } else if (get_tree(other_node) == color::black()) {
              // we do this to get shorter paths. check if we
              // are nearer to the source as its parent is
              if (is_closer_to_terminal(current_node, other_node)) {
                set_edge_to_parent(other_node, out_edge);
                put(dist_map_, other_node, get(dist_map_, current_node) + 1);
                put(time_map_, other_node, get(time_map_, current_node));
              }
            } else {
              assert(get_tree(other_node) == color::white());
              // kewl, found a path from one to the other
              // search tree, return
              // the connecting edge in src->sink dir
              return {out_edge, true};
            }
          }
        }  // for all out-edges
      }    // source-tree-growing
      else {
        assert(get_tree(current_node) == color::white());
        if (current_node != last_grow_vertex_) {
          last_grow_vertex_ = current_node;
          last_grow_edge_range_ = std::make_unique<partial_view<out_edge_range>>(out_edges(current_node, g_));
        }
        for (; !last_grow_edge_range_->empty(); last_grow_edge_range_->move_to_next()) {
          edge_descriptor in_edge = get(rev_edge_map_, *last_grow_edge_range_->begin());
          if (get(res_cap_map_, in_edge) > 0) {  // check if there is capacity left
            vertex_descriptor other_node = source(in_edge, g_);
            if (get_tree(other_node) == color::gray()) {  // it's a free node
              set_tree(other_node,
                       color::white());                                      // aquire that node
                                                                             // to our search
                                                                             // tree
              set_edge_to_parent(other_node, in_edge);                       // set us as parent
              add_active_node(other_node);                                   // activate that node
              put(dist_map_, other_node, get(dist_map_, current_node) + 1);  // set its distance
              put(time_map_, other_node, get(time_map_, current_node));      // and time
            } else if (get_tree(other_node) == color::white()) {
              if (is_closer_to_terminal(current_node, other_node)) {
                // we are closer to the sink than its parent
                // is, so we "adopt" him
                set_edge_to_parent(other_node, in_edge);
                put(dist_map_, other_node, get(dist_map_, current_node) + 1);
                put(time_map_, other_node, get(time_map_, current_node));
              }
            } else {
              assert(get_tree(other_node) == color::black());
              // kewl, found a path from one to the other
              // search tree,
              // return the connecting edge in src->sink dir
              return {in_edge, true};
            }
          }
        }  // for all out-edges
      }    // sink-tree growing

      // all edges of that node are processed, and no more paths were
      // found.
      // remove if from the front of the active queue
      finish_node(current_node);
    }  // while active_nodes not empty

    // no active nodes anymore and no path found, we're done
    return {edge_descriptor(), false};
  }

  /**
   * augments path from s->t and updates residual graph
   * source(e, g_) is the end of the path found in the source-tree
   * target(e, g_) is the beginning of the path found in the sink-tree
   * this phase generates orphans on satured edges, if the attached verts
   * are from different search-trees orphans are ordered in distance to
   * sink/source. first the farest from the source are front_inserted into
   * the orphans list, and after that the sink-tree-orphans are
   * front_inserted. when going to adoption stage the orphans are
   * popped_front, and so we process the nearest verts to the terminals
   * first
   */
  void augment(edge_descriptor e) {
    assert(get_tree(target(e, g_)) == color::white());
    assert(get_tree(source(e, g_)) == color::black());
    assert(orphans_.empty());

    const edge_flow_value bottleneck = find_bottleneck(e);
    // now we push the found flow through the path
    // for each edge we saturate we have to look for the verts that
    // belong to that edge, one of them becomes an orphans now process
    // the connecting edge
    put(res_cap_map_, e, get(res_cap_map_, e) - bottleneck);
    assert(get(res_cap_map_, e) >= 0);
    put(res_cap_map_, get(rev_edge_map_, e), get(res_cap_map_, get(rev_edge_map_, e)) + bottleneck);

    // now we follow the path back to the source
    vertex_descriptor current_node = source(e, g_);
    while (current_node != source_) {
      edge_descriptor pred = get_edge_to_parent(current_node);
      put(res_cap_map_, pred, get(res_cap_map_, pred) - bottleneck);
      assert(get(res_cap_map_, pred) >= 0);
      put(res_cap_map_, get(rev_edge_map_, pred), get(res_cap_map_, get(rev_edge_map_, pred)) + bottleneck);
      if (get(res_cap_map_, pred) == 0) {
        set_no_parent(current_node);
        orphans_.push_front(current_node);
      }
      current_node = source(pred, g_);
    }
    // then go forward in the sink-tree
    current_node = target(e, g_);
    while (current_node != sink_) {
      edge_descriptor pred = get_edge_to_parent(current_node);
      put(res_cap_map_, pred, get(res_cap_map_, pred) - bottleneck);
      assert(get(res_cap_map_, pred) >= 0);
      put(res_cap_map_, get(rev_edge_map_, pred), get(res_cap_map_, get(rev_edge_map_, pred)) + bottleneck);
      if (get(res_cap_map_, pred) == 0) {
        set_no_parent(current_node);
        orphans_.push_front(current_node);
      }
      current_node = target(pred, g_);
    }
    // and add it to the max-flow
    flow_ += bottleneck;
  }

  /**
   * returns the bottleneck of a s->t path (end_of_path is last vertex in
   * source-tree, begin_of_path is first vertex in sink-tree)
   */
  edge_flow_value find_bottleneck(edge_descriptor e) {
    edge_flow_value minimum_cap = get(res_cap_map_, e);
    vertex_descriptor current_node = source(e, g_);
    // first go back in the source tree
    while (current_node != source_) {
      edge_descriptor pred = get_edge_to_parent(current_node);
      minimum_cap = std::min(minimum_cap, get(res_cap_map_, pred));
      current_node = source(pred, g_);
    }
    // then go forward in the sink-tree
    current_node = target(e, g_);
    while (current_node != sink_) {
      edge_descriptor pred = get_edge_to_parent(current_node);
      minimum_cap = std::min(minimum_cap, get(res_cap_map_, pred));
      current_node = target(pred, g_);
    }
    return minimum_cap;
  }

  /**
   * rebuild search trees
   * empty the queue of orphans, and find new parents for them or just
   * drop them from the search trees
   */
  void adopt() {
    while (!orphans_.empty() || !child_orphans_.empty()) {
      vertex_descriptor current_node;
      if (child_orphans_.empty()) {
        // get the next orphan from the main-queue  and remove it
        current_node = orphans_.front();
        orphans_.pop_front();
      } else {
        current_node = child_orphans_.front();
        child_orphans_.pop();
      }
      if (get_tree(current_node) == color::black()) {
        // we're in the source-tree
        distance_value min_distance = std::numeric_limits<distance_value>::max();
        edge_descriptor new_parent_edge;
        for (auto e : out_edges(current_node, g_)) {
          const edge_descriptor in_edge = get(rev_edge_map_, e);
          // we should be the target of this edge
          assert(target(in_edge, g_) == current_node);
          if (get(res_cap_map_, in_edge) > 0) {
            vertex_descriptor other_node = source(in_edge, g_);
            if (get_tree(other_node) == color::black() && has_source_connect(other_node) &&
                get(dist_map_, other_node) < min_distance) {
              min_distance = get(dist_map_, other_node);
              new_parent_edge = in_edge;
            }
          }
        }
        if (min_distance != std::numeric_limits<distance_value>::max()) {
          set_edge_to_parent(current_node, new_parent_edge);
          put(dist_map_, current_node, min_distance + 1);
          put(time_map_, current_node, time_);
        } else {
          put(time_map_, current_node, 0);
          for (auto e : out_edges(current_node, g_)) {
            edge_descriptor in_edge = get(rev_edge_map_, e);
            vertex_descriptor other_node = source(in_edge, g_);
            if (get_tree(other_node) == color::black() && other_node != source_) {
              if (get(res_cap_map_, in_edge) > 0) {
                add_active_node(other_node);
              }
              if (has_parent(other_node) && source(get_edge_to_parent(other_node), g_) == current_node) {
                // we are the parent of that node
                // it has to find a new parent, too
                set_no_parent(other_node);
                child_orphans_.push(other_node);
              }
            }
          }
          set_tree(current_node, color::gray());
        }  // no parent found
      }    // source-tree-adoption
      else {
        // now we should be in the sink-tree, check that...
        assert(get_tree(current_node) == color::white());
        edge_descriptor new_parent_edge;
        distance_value min_distance = std::numeric_limits<distance_value>::max();
        for (auto out_edge : out_edges(current_node, g_)) {
          if (get(res_cap_map_, out_edge) > 0) {
            const vertex_descriptor other_node = target(out_edge, g_);
            if (get_tree(other_node) == color::white() && has_sink_connect(other_node)) {
              if (get(dist_map_, other_node) < min_distance) {
                min_distance = get(dist_map_, other_node);
                new_parent_edge = out_edge;
              }
            }
          }
        }
        if (min_distance != std::numeric_limits<distance_value>::max()) {
          set_edge_to_parent(current_node, new_parent_edge);
          put(dist_map_, current_node, min_distance + 1);
          put(time_map_, current_node, time_);
        } else {
          put(time_map_, current_node, 0);
          for (auto out_edge : out_edges(current_node, g_)) {
            const vertex_descriptor other_node = target(out_edge, g_);
            if (get_tree(other_node) == color::white() && other_node != sink_) {
              if (get(res_cap_map_, out_edge) > 0) {
                add_active_node(other_node);
              }
              if (has_parent(other_node) && target(get_edge_to_parent(other_node), g_) == current_node) {
                // we were it's parent, so it has to find a
                // new one, too
                set_no_parent(other_node);
                child_orphans_.push(other_node);
              }
            }
          }
          set_tree(current_node, color::gray());
        }  // no parent found
      }    // sink-tree adoption
    }      // while !orphans.empty()
  }        // adopt

  /**
   * return next active vertex if there is one, otherwise a null_vertex
   */
  vertex_descriptor get_next_active_node() {
    while (true) {
      if (active_nodes_.empty()) {
        return graph_traits<Graph>::null_vertex();
      }
      vertex_descriptor v = active_nodes_.front();

      // if it has no parent, this node can't be active (if its not
      // source or sink)
      if (!has_parent(v) && v != source_ && v != sink_) {
        active_nodes_.pop();
        put(in_active_list_map_, v, false);
      } else {
        assert(get_tree(v) == color::black() || get_tree(v) == color::white());
        return v;
      }
    }
  }

  /**
   * adds v as an active vertex, but only if its not in the list already
   */
  void add_active_node(vertex_descriptor v) {
    assert(get_tree(v) != color::gray());
    if (get(in_active_list_map_, v)) {
      if (last_grow_vertex_ == v) {
        last_grow_vertex_ = graph_traits<Graph>::null_vertex();
      }
      return;
    }

    put(in_active_list_map_, v, true);
    active_nodes_.push(v);
  }

  /**
   * finish_node removes a node from the front of the active queue (its
   * called in grow phase, if no more paths can be found using this node)
   */
  void finish_node(vertex_descriptor v) {
    assert(active_nodes_.front() == v);
    active_nodes_.pop();
    put(in_active_list_map_, v, false);
    last_grow_vertex_ = graph_traits<Graph>::null_vertex();
  }

  /**
   * removes a vertex from the queue of active nodes (actually this does
   * nothing, but checks if this node has no parent edge, as this is the
   * criteria for being no more active)
   */
  void remove_active_node(vertex_descriptor v) { assert(!has_parent(v)); }

  /**
   * returns the search tree of v; color_value::black() for source tree,
   * white() for sink tree, gray() for no tree
   */
  color_value get_tree(vertex_descriptor v) const { return get(tree_map_, v); }

  /**
   * sets search tree of v; color_value::black() for source tree, white()
   * for sink tree, gray() for no tree
   */
  void set_tree(vertex_descriptor v, color_value t) { put(tree_map_, v, t); }

  /**
   * returns edge to parent vertex of v;
   */
  edge_descriptor get_edge_to_parent(vertex_descriptor v) const { return get(pre_map_, v); }

  /**
   * returns true if the edge stored in pre_map_[v] is a valid entry
   */
  bool has_parent(vertex_descriptor v) const { return get(has_parent_map_, v); }

  /**
   * sets edge to parent vertex of v;
   */
  void set_edge_to_parent(vertex_descriptor v, edge_descriptor f_edge_to_parent) {
    assert(get(res_cap_map_, f_edge_to_parent) > 0);
    put(pre_map_, v, f_edge_to_parent);
    put(has_parent_map_, v, true);
  }

  /**
   * removes the edge to parent of v (this is done by invalidating the
   * entry an additional map)
   */
  void set_no_parent(vertex_descriptor v) { put(has_parent_map_, v, false); }

  /**
   * checks if vertex v has a connect to the sink-vertex (@var sink_)
   * @param v the vertex which is checked
   * @return true if a path to the sink was found, false if not
   */
  bool has_sink_connect(vertex_descriptor v) {
    distance_value current_distance = 0;
    vertex_descriptor current_vertex = v;
    while (true) {
      if (get(time_map_, current_vertex) == time_) {
        // we found a node which was already checked this round. use
        // it for distance calculations
        current_distance += get(dist_map_, current_vertex);
        break;
      }
      if (current_vertex == sink_) {
        put(time_map_, sink_, time_);
        break;
      }
      if (has_parent(current_vertex)) {
        // it has a parent, so get it
        current_vertex = target(get_edge_to_parent(current_vertex), g_);
        ++current_distance;
      } else {
        // no path found
        return false;
      }
    }
    current_vertex = v;
    while (get(time_map_, current_vertex) != time_) {
      put(dist_map_, current_vertex, current_distance);
      --current_distance;
      put(time_map_, current_vertex, time_);
      current_vertex = target(get_edge_to_parent(current_vertex), g_);
    }
    return true;
  }

  /**
   * checks if vertex v has a connect to the source-vertex (@var source_)
   * @param v the vertex which is checked
   * @return true if a path to the source was found, false if not
   */
  bool has_source_connect(vertex_descriptor v) {
    distance_value current_distance = 0;
    vertex_descriptor current_vertex = v;
    while (true) {
      if (get(time_map_, current_vertex) == time_) {
        // we found a node which was already checked this round. use
        // it for distance calculations
        current_distance += get(dist_map_, current_vertex);
        break;
      }
      if (current_vertex == source_) {
        put(time_map_, source_, time_);
        break;
      }
      if (has_parent(current_vertex)) {
        // it has a parent, so get it
        current_vertex = source(get_edge_to_parent(current_vertex), g_);
        ++current_distance;
      } else {
        // no path found
        return false;
      }
    }
    current_vertex = v;
    while (get(time_map_, current_vertex) != time_) {
      put(dist_map_, current_vertex, current_distance);
      --current_distance;
      put(time_map_, current_vertex, time_);
      current_vertex = source(get_edge_to_parent(current_vertex), g_);
    }
    return true;
  }

  /**
   * returns true, if p is closer to a terminal than q
   */
  bool is_closer_to_terminal(vertex_descriptor p, vertex_descriptor q) {
    // checks the timestamps first, to build no cycles, and after that
    // the real distance
    return (get(time_map_, q) <= get(time_map_, p) && get(dist_map_, q) > get(dist_map_, p) + 1);
  }

  ////////
  // member vars
  ////////
  Graph& g_;
  IndexMap index_map_;
  EdgeCapacityMap cap_map_;
  ResidualCapacityEdgeMap res_cap_map_;
  ReverseEdgeMap rev_edge_map_;
  PredecessorMap pre_map_;  // stores paths found in the growth stage
  ColorMap tree_map_;       // maps each vertex into one of the two search tree
                            // or none (gray())
  DistanceMap dist_map_;    // stores distance to source/sink nodes
  vertex_descriptor source_;
  vertex_descriptor sink_;

  vertex_queue active_nodes_;
  vector_property_map<bool, IndexMap> in_active_list_map_;

  std::list<vertex_descriptor> orphans_;
  vertex_queue child_orphans_;  // we use a second queuqe for child orphans, as
                                // they are FIFO processed

  vector_property_map<bool, IndexMap> has_parent_map_;

  // timestamp of each node, used for
  // sink/source-path calculations
  vector_property_map<std::int64_t, IndexMap> time_map_;
  edge_flow_value flow_ = 0;
  std::int64_t time_ = 1;
  vertex_descriptor last_grow_vertex_ = graph_traits<Graph>::null_vertex();
  std::unique_ptr<partial_view<out_edge_range>> last_grow_edge_range_;
};

}  // namespace max_flow_detail

/**
 * non-named-parameter version, given everything
 * this is the catch all version
 */
template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> CapacityEdgeMap,
          concepts::ReadWriteEdgePropertyMap<G> ResidualCapacityEdgeMap,
          concepts::ReadableEdgePropertyMap<G> ReverseEdgeMap, concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadableVertexPropertyMap<G> IndexMap>
requires concepts::VertexAndEdgeListGraph<G> property_traits_value_t<CapacityEdgeMap> boykov_kolmogorov_max_flow(
    G& g, CapacityEdgeMap cap, ResidualCapacityEdgeMap res_cap, ReverseEdgeMap rev_map, PredecessorMap pre_map,
    ColorMap color, DistanceMap dist, IndexMap idx, graph_vertex_descriptor_t<G> src,
    graph_vertex_descriptor_t<G> sink) {
  assert(num_vertices(g) >= 2);
  assert(src != sink);

  max_flow_detail::bk_max_flow algo{g, cap, res_cap, rev_map, pre_map, color, dist, idx, src, sink};
  return algo.max_flow();
}

/**
 * non-named-parameter version, given capacity, residual_capacity,
 * reverse_edges, and an index map.
 */
template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> CapacityEdgeMap,
          concepts::ReadWriteEdgePropertyMap<G> ResidualCapacityEdgeMap,
          concepts::ReadableEdgePropertyMap<G> ReverseEdgeMap, concepts::ReadableVertexPropertyMap<G> IndexMap>
requires concepts::VertexAndEdgeListGraph<G> property_traits_value_t<CapacityEdgeMap> boykov_kolmogorov_max_flow(
    G& g, CapacityEdgeMap cap, ResidualCapacityEdgeMap res_cap, ReverseEdgeMap rev, IndexMap idx,
    graph_vertex_descriptor_t<G> src, graph_vertex_descriptor_t<G> sink) {
  std::size_t n_verts = num_vertices(g);
  return boykov_kolmogorov_max_flow(g, cap, res_cap, rev,
                                    make_vector_property_map(n_verts, idx, graph_edge_descriptor_t<G>{}),
                                    make_vector_property_map(n_verts, idx, default_color_type::white_color),
                                    make_vector_property_map(n_verts, idx, std::size_t{0}), idx, src, sink);
}

/**
 * non-named-parameter version, some given: capacity, residual_capacity,
 * reverse_edges, color_map and an index map. Use this if you are interested in
 * the minimum cut, as the color map provides that info.
 */
template <concepts::IncidenceGraph G, concepts::ReadableEdgePropertyMap<G> CapacityEdgeMap,
          concepts::ReadWriteEdgePropertyMap<G> ResidualCapacityEdgeMap,
          concepts::ReadableEdgePropertyMap<G> ReverseEdgeMap, concepts::ReadWriteVertexPropertyMap<G> ColorMap,
          concepts::ReadableVertexPropertyMap<G> IndexMap>
requires concepts::VertexAndEdgeListGraph<G> property_traits_value_t<CapacityEdgeMap> boykov_kolmogorov_max_flow(
    G& g, CapacityEdgeMap cap, ResidualCapacityEdgeMap res_cap, ReverseEdgeMap rev, ColorMap color, IndexMap idx,
    graph_vertex_descriptor_t<G> src, graph_vertex_descriptor_t<G> sink) {
  std::size_t n_verts = num_vertices(g);
  return boykov_kolmogorov_max_flow(g, cap, res_cap, rev,
                                    make_vector_property_map(n_verts, idx, graph_edge_descriptor_t<G>{}), color,
                                    make_vector_property_map(n_verts, idx, std::size_t{0}), idx, src, sink);
}

template <concepts::IncidenceGraph G>
requires concepts::VertexAndEdgeListGraph<G>
auto boykov_kolmogorov_max_flow(G& g, graph_vertex_descriptor_t<G> src, graph_vertex_descriptor_t<G> sink) {
  return boykov_kolmogorov_max_flow(g, get(edge_capacity, g), get(edge_residual_capacity, g), get(edge_reverse, g),
                                    get(vertex_predecessor, g), get(vertex_color, g), get(vertex_distance, g),
                                    get(vertex_index, g), src, sink);
}

}  // namespace bagl

#endif  // BAGL_BAGL_BOYKOV_KOLMOGOROV_MAX_FLOW_H_
