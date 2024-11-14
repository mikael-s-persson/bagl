// Copyright (c) Aaron Windsor 2007
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PLANAR_DETAIL_BOYER_MYRVOLD_IMPL_H_
#define BAGL_BAGL_PLANAR_DETAIL_BOYER_MYRVOLD_IMPL_H_

#include <iterator>
#include <list>
#include <memory>
#include <tuple>
#include <vector>

#include "bagl/depth_first_search.h"
#include "bagl/graph_traits.h"
#include "bagl/planar_detail/bucket_sort.h"
#include "bagl/planar_detail/face_handles.h"
#include "bagl/planar_detail/face_iterators.h"
#include "bagl/property_map.h"

namespace bagl {
enum class bm_case_t { bm_no_case_chosen, bm_case_a, bm_case_b, bm_case_c, bm_case_d, bm_case_e };

template <typename LowPointMap, typename DFSParentMap, typename DFSNumberMap, typename LeastAncestorMap,
          typename DFSParentEdgeMap>
struct planar_dfs_visitor : public dfs_visitor<> {
  planar_dfs_visitor(LowPointMap lpm, DFSParentMap dfs_p, DFSNumberMap dfs_n, LeastAncestorMap lam,
                     DFSParentEdgeMap dfs_edge)
      : low_(lpm), parent_(dfs_p), df_number_(dfs_n), least_ancestor_(lam), df_edge_(dfs_edge), count_(0) {}

  template <typename Vertex, typename Graph>
  void start_vertex(const Vertex& u, Graph& /*unused*/) {
    put(parent_, u, u);
    put(least_ancestor_, u, count_);
  }

  template <typename Vertex, typename Graph>
  void discover_vertex(const Vertex& u, Graph& /*unused*/) {
    put(low_, u, count_);
    put(df_number_, u, count_);
    ++count_;
  }

  template <typename Edge, typename Graph>
  void tree_edge(const Edge& e, Graph& g) {
    using vertex_t = typename graph_traits<Graph>::vertex_descriptor;
    vertex_t s(source(e, g));
    vertex_t t(target(e, g));

    put(parent_, t, s);
    put(df_edge_, t, e);
    put(least_ancestor_, t, get(df_number_, s));
  }

  template <typename Edge, typename Graph>
  void back_edge(const Edge& e, Graph& g) {
    using vertex_t = typename graph_traits<Graph>::vertex_descriptor;

    vertex_t s(source(e, g));
    vertex_t t(target(e, g));

    if (t != get(parent_, s)) {
      std::size_t s_low_df_number = get(low_, s);
      std::size_t t_df_number = get(df_number_, t);
      std::size_t s_least_ancestor_df_number = get(least_ancestor_, s);

      put(low_, s, std::min(s_low_df_number, t_df_number));

      put(least_ancestor_, s, std::min(s_least_ancestor_df_number, t_df_number));
    }
  }

  template <typename Vertex, typename Graph>
  void finish_vertex(const Vertex& u, Graph& /*unused*/) {
    Vertex u_parent = get(parent_, u);
    std::size_t u_parent_lowpoint = get(low_, u_parent);
    std::size_t u_lowpoint = get(low_, u);

    if (u_parent != u) {
      put(low_, u_parent, std::min(u_lowpoint, u_parent_lowpoint));
    }
  }

  LowPointMap low_;
  DFSParentMap parent_;
  DFSNumberMap df_number_;
  LeastAncestorMap least_ancestor_;
  DFSParentEdgeMap df_edge_;
  std::size_t count_;
};

template <typename Graph, typename VertexIndexMap, typename StoreOldHandlesPolicy = planar_detail::store_old_handles,
          typename StoreEmbeddingPolicy = planar_detail::recursive_lazy_list>
class boyer_myrvold_impl {
  using vertex_t = typename graph_traits<Graph>::vertex_descriptor;
  using edge_t = typename graph_traits<Graph>::edge_descriptor;
  using face_handle_t = planar_detail::face_handle<Graph, StoreOldHandlesPolicy, StoreEmbeddingPolicy>;
  using vertex_vector_t = std::vector<vertex_t>;
  using edge_vector_t = std::vector<edge_t>;
  using vertex_list_t = std::list<vertex_t>;
  using face_handle_list_t = std::list<face_handle_t>;
  using face_handle_list_ptr_t = std::shared_ptr<face_handle_list_t>;
  using vertex_list_ptr_t = std::shared_ptr<vertex_list_t>;
  using merge_stack_frame_t = std::tuple<vertex_t, bool, bool>;
  using merge_stack_t = std::vector<merge_stack_frame_t>;

  template <typename T>
  using map_vertex_to = vector_property_map<T, VertexIndexMap>;

  using vertex_to_v_size_map_t = map_vertex_to<std::size_t>;
  using vertex_to_vertex_map_t = map_vertex_to<vertex_t>;
  using vertex_to_edge_map_t = map_vertex_to<edge_t>;
  using vertex_to_vertex_list_ptr_map_t = map_vertex_to<vertex_list_ptr_t>;
  using vertex_to_edge_vector_map_t = map_vertex_to<edge_vector_t>;
  using vertex_to_bool_map_t = map_vertex_to<bool>;
  using vertex_to_face_handle_map_t = map_vertex_to<face_handle_t>;
  using vertex_to_face_handle_list_ptr_map_t = map_vertex_to<face_handle_list_ptr_t>;
  using vertex_to_separated_node_map_t = map_vertex_to<typename vertex_list_t::iterator>;

  template <typename BicompSideToTraverse = single_side, typename VisitorType = lead_visitor,
            typename Time = current_iteration>
  struct face_vertex_iterator {
    using type = face_iterator<Graph, property_map_ref<vertex_to_face_handle_map_t>, vertex_t, BicompSideToTraverse,
                               VisitorType, Time>;
  };

  template <typename BicompSideToTraverse = single_side, typename Time = current_iteration>
  struct face_edge_iterator {
    using type = face_iterator<Graph, property_map_ref<vertex_to_face_handle_map_t>, edge_t, BicompSideToTraverse,
                               lead_visitor, Time>;
  };

 public:
  boyer_myrvold_impl(const Graph& g, VertexIndexMap vm)
      : g_(g),
        vm_(vm),

        vertices_by_dfs_num_(num_vertices(g_)),

        low_point_(num_vertices(g_), vm_),
        dfs_parent_(num_vertices(g_), vm_),
        dfs_number_(num_vertices(g_), vm_),
        least_ancestor_(num_vertices(g_), vm_),
        pertinent_roots_(num_vertices(g_), vm_),
        backedge_flag_(num_vertices(g_), vm_, num_vertices(g_) + 1),
        visited_(num_vertices(g_), vm_, num_vertices(g_) + 1),
        face_handles_(num_vertices(g_), vm_),
        dfs_child_handles_(num_vertices(g_), vm_),
        separated_dfs_child_list_(num_vertices(g_), vm_),
        separated_node_in_parent_list_(num_vertices(g_), vm_),
        canonical_dfs_child_(num_vertices(g_), vm_),
        flipped_(num_vertices(g_), vm_, false),
        backedges_(num_vertices(g_), vm_),
        dfs_parent_edge_(num_vertices(g_), vm_)

  {
    planar_dfs_visitor vis(low_point_.ref(), dfs_parent_.ref(), dfs_number_.ref(), least_ancestor_.ref(),
                           dfs_parent_edge_.ref());

    // Perform a depth-first search to find each vertex's low point, least
    // ancestor, and dfs tree information
    depth_first_search(g_, vis, two_bit_color_map(num_vertices(g_), vm_).ref());

    // Sort vertices by their lowpoint - need this later in the constructor
    vertex_vector_t vertices_by_lowpoint(num_vertices(g_));
    auto v_rg = vertices(g_);
    std::copy(v_rg.begin(), v_rg.end(), vertices_by_lowpoint.begin());
    bucket_sort(vertices_by_lowpoint.begin(), vertices_by_lowpoint.end(), low_point_.ref());

    // Sort vertices by their dfs number - need this to iterate by reverse
    // DFS number in the main loop.
    std::copy(v_rg.begin(), v_rg.end(), vertices_by_dfs_num_.begin());
    bucket_sort(vertices_by_dfs_num_.begin(), vertices_by_dfs_num_.end(), dfs_number_.ref());

    // Initialize face handles. A face handle is an abstraction that serves
    // two uses in our implementation - it allows us to efficiently move
    // along the outer face of embedded bicomps in a partially embedded
    // graph, and it provides storage for the planar embedding. Face
    // handles are implemented by a sequence of edges and are associated
    // with a particular vertex - the sequence of edges represents the
    // current embedding of edges around that vertex, and the first and
    // last edges in the sequence represent the pair of edges on the outer
    // face that are adjacent to the associated vertex. This lets us embed
    // edges in the graph by just pushing them on the front or back of the
    // sequence of edges held by the face handles.
    //
    // Our algorithm starts with a DFS tree of edges (where every vertex is
    // an articulation point and every edge is a singleton bicomp) and
    // repeatedly merges bicomps by embedding additional edges. Note that
    // any bicomp at any point in the algorithm can be associated with a
    // unique edge connecting the vertex of that bicomp with the lowest DFS
    // number (which we refer to as the "root" of the bicomp) with its DFS
    // child in the bicomp: the existence of two such edges would contradict
    // the properties of a DFS tree. We refer to the DFS child of the root
    // of a bicomp as the "canonical DFS child" of the bicomp. Note that a
    // vertex can be the root of more than one bicomp.
    //
    // We move around the external faces of a bicomp using a few property
    // maps, which we'll initialize presently:
    //
    // - face_handles_: maps a vertex to a face handle that can be used to
    //   move "up" a bicomp. For a vertex that isn't an articulation point,
    //   this holds the face handles that can be used to move around that
    //   vertex's unique bicomp. For a vertex that is an articulation point,
    //   this holds the face handles associated with the unique bicomp that
    //   the vertex is NOT the root of. These handles can therefore be used
    //   to move from any point on the outer face of the tree of bicomps
    //   around the current outer face towards the root of the DFS tree.
    //
    // - dfs_child_handles_: these are used to hold face handles for
    //   vertices that are articulation points - dfs_child_handles_[v] holds
    //   the face handles corresponding to vertex u in the bicomp with root
    //   u and canonical DFS child v.
    //
    // - canonical_dfs_child_: this property map allows one to determine the
    //   canonical DFS child of a bicomp while traversing the outer face.
    //   This property map is only valid when applied to one of the two
    //   vertices adjacent to the root of the bicomp on the outer face. To
    //   be more precise, if v is the canonical DFS child of a bicomp,
    //   canonical_dfs_child_[dfs_child_handles_[v].first_vertex()] == v and
    //   canonical_dfs_child_[dfs_child_handles_[v].second_vertex()] == v.
    //
    // - pertinent_roots_: given a vertex v, pertinent_roots_[v] contains a
    //   list of face handles pointing to the top of bicomps that need to
    //   be visited by the current walkdown traversal (since they lead to
    //   backedges that need to be embedded). These lists are populated by
    //   the walkup and consumed by the walkdown.

    for (auto v : vertices(g_)) {
      vertex_t parent = dfs_parent_[v];

      if (parent != v) {
        edge_t parent_edge = dfs_parent_edge_[v];
        add_to_embedded_edges(parent_edge, StoreOldHandlesPolicy());
        face_handles_[v] = face_handle_t(v, parent_edge, g_);
        dfs_child_handles_[v] = face_handle_t(parent, parent_edge, g_);
      } else {
        face_handles_[v] = face_handle_t(v);
        dfs_child_handles_[v] = face_handle_t(parent);
      }

      canonical_dfs_child_[v] = v;
      pertinent_roots_[v] = face_handle_list_ptr_t(new face_handle_list_t);
      separated_dfs_child_list_[v] = vertex_list_ptr_t(new vertex_list_t);
    }

    // We need to create a list of not-yet-merged depth-first children for
    // each vertex that will be updated as bicomps get merged. We sort each
    // list by ascending lowpoint, which allows the externally_active
    // function to run in constant time, and we keep a pointer to each
    // vertex's representation in its parent's list, which allows merging
    // in constant time.

    for (auto v : vertices_by_lowpoint) {
      vertex_t parent(dfs_parent_[v]);
      if (v != parent) {
        separated_node_in_parent_list_[v] =
            separated_dfs_child_list_[parent]->insert(separated_dfs_child_list_[parent]->end(), v);
      }
    }

    // The merge stack holds path information during a walkdown iteration
    merge_stack_.reserve(num_vertices(g_));
  }

  bool is_planar() {
    // This is the main algorithm: starting with a DFS tree of embedded
    // edges (which, since it's a tree, is planar), iterate through all
    // vertices by reverse DFS number, attempting to embed all backedges
    // connecting the current vertex to vertices with higher DFS numbers.
    //
    // The walkup is a procedure that examines all such backedges and sets
    // up the required data structures so that they can be searched by the
    // walkdown in linear time. The walkdown does the actual work of
    // embedding edges and flipping bicomps, and can identify when it has
    // come across a kuratowski subgraph.
    //
    // store_old_face_handles caches face handles from the previous
    // iteration - this is used only for the kuratowski subgraph isolation,
    // and is therefore dispatched based on the StoreOldHandlesPolicy.
    //
    // clean_up_embedding does some clean-up and fills in values that have
    // to be computed lazily during the actual execution of the algorithm
    // (for instance, whether or not a bicomp is flipped in the final
    // embedding). It's dispatched on the the StoreEmbeddingPolicy, since
    // it's not needed if an embedding isn't desired.
    for (auto v : std::views::all(vertices_by_dfs_num_) | std::views::reverse) {
      store_old_face_handles(StoreOldHandlesPolicy());

      walkup(v);

      if (!walkdown(v)) {
        return false;
      }
    }

    clean_up_embedding(StoreEmbeddingPolicy());

    return true;
  }

 private:
  void walkup(vertex_t v) {
    // The point of the walkup is to follow all backedges from v to
    // vertices with higher DFS numbers, and update pertinent_roots_
    // for the bicomp roots on the path from backedge endpoints up
    // to v. This will set the stage for the walkdown to efficiently
    // traverse the graph of bicomps down from v.

    using walkup_iterator_t = typename face_vertex_iterator<both_sides>::type;

    for (auto e : out_edges(v, g_)) {
      vertex_t e_source(source(e, g_));
      vertex_t e_target(target(e, g_));

      if (e_source == e_target) {
        self_loops_.push_back(e);
        continue;
      }

      vertex_t w(e_source == v ? e_target : e_source);

      // continue if not a back edge or already embedded
      if (dfs_number_[w] < dfs_number_[v] || e == dfs_parent_edge_[w]) {
        continue;
      }

      backedges_[w].push_back(e);

      std::size_t timestamp = dfs_number_[v];
      backedge_flag_[w] = timestamp;

      walkup_iterator_t walkup_itr(w, face_handles_.ref());
      walkup_iterator_t walkup_end;
      vertex_t lead_vertex = w;

      while (true) {
        // Move to the root of the current bicomp or the first visited
        // vertex on the bicomp by going up each side in parallel

        while (walkup_itr != walkup_end && visited_[*walkup_itr] != timestamp) {
          lead_vertex = *walkup_itr;
          visited_[lead_vertex] = timestamp;
          ++walkup_itr;
        }

        // If we've found the root of a bicomp through a path we haven't
        // seen before, update pertinent_roots_ with a handle to the
        // current bicomp. Otherwise, we've just seen a path we've been
        // up before, so break out of the main while loop.

        if (walkup_itr == walkup_end) {
          vertex_t dfs_child = canonical_dfs_child_[lead_vertex];
          vertex_t parent = dfs_parent_[dfs_child];

          visited_[dfs_child_handles_[dfs_child].first_vertex()] = timestamp;
          visited_[dfs_child_handles_[dfs_child].second_vertex()] = timestamp;

          if (low_point_[dfs_child] < dfs_number_[v] || least_ancestor_[dfs_child] < dfs_number_[v]) {
            pertinent_roots_[parent]->push_back(dfs_child_handles_[dfs_child]);
          } else {
            pertinent_roots_[parent]->push_front(dfs_child_handles_[dfs_child]);
          }

          if (parent != v && visited_[parent] != timestamp) {
            walkup_itr = walkup_iterator_t(parent, face_handles_.ref());
            lead_vertex = parent;
          } else {
            break;
          }
        } else {
          break;
        }
      }
    }
  }

  bool walkdown(vertex_t v) {
    // This procedure is where all of the action is - pertinent_roots_
    // has already been set up by the walkup, so we just need to move
    // down bicomps from v until we find vertices that have been
    // labeled as backedge endpoints. Once we find such a vertex, we
    // embed the corresponding edge and glue together the bicomps on
    // the path connecting the two vertices in the edge. This may
    // involve flipping bicomps along the way.

    vertex_t w;  // the other endpoint of the edge we're embedding

    while (!pertinent_roots_[v]->empty()) {
      face_handle_t root_face_handle = pertinent_roots_[v]->front();
      face_handle_t curr_face_handle = root_face_handle;
      pertinent_roots_[v]->pop_front();

      merge_stack_.clear();

      while (true) {
        typename face_vertex_iterator<>::type first_face_itr;
        typename face_vertex_iterator<>::type second_face_itr;
        typename face_vertex_iterator<>::type face_end;
        vertex_t first_side_vertex = graph_traits<Graph>::null_vertex();
        vertex_t second_side_vertex = graph_traits<Graph>::null_vertex();
        vertex_t first_tail;
        vertex_t second_tail;

        first_tail = second_tail = curr_face_handle.get_anchor();
        first_face_itr = typename face_vertex_iterator<>::type(curr_face_handle, face_handles_.ref(), first_side());
        second_face_itr = typename face_vertex_iterator<>::type(curr_face_handle, face_handles_.ref(), second_side());

        for (; first_face_itr != face_end; ++first_face_itr) {
          vertex_t face_vertex(*first_face_itr);
          if (pertinent(face_vertex, v) || externally_active(face_vertex, v)) {
            first_side_vertex = face_vertex;
            second_side_vertex = face_vertex;
            break;
          }
          first_tail = face_vertex;
        }

        if (first_side_vertex == graph_traits<Graph>::null_vertex() ||
            first_side_vertex == curr_face_handle.get_anchor()) {
          break;
        }

        for (; second_face_itr != face_end; ++second_face_itr) {
          vertex_t face_vertex(*second_face_itr);
          if (pertinent(face_vertex, v) || externally_active(face_vertex, v)) {
            second_side_vertex = face_vertex;
            break;
          }
          second_tail = face_vertex;
        }

        vertex_t chosen;
        bool chose_first_upper_path = false;
        if (internally_active(first_side_vertex, v)) {
          chosen = first_side_vertex;
          chose_first_upper_path = true;
        } else if (internally_active(second_side_vertex, v)) {
          chosen = second_side_vertex;
          chose_first_upper_path = false;
        } else if (pertinent(first_side_vertex, v)) {
          chosen = first_side_vertex;
          chose_first_upper_path = true;
        } else if (pertinent(second_side_vertex, v)) {
          chosen = second_side_vertex;
          chose_first_upper_path = false;
        } else {
          // If there's a pertinent vertex on the lower face
          // between the first_face_itr and the second_face_itr,
          // this graph isn't planar.
          for (; *first_face_itr != second_side_vertex; ++first_face_itr) {
            vertex_t p(*first_face_itr);
            if (pertinent(p, v)) {
              // Found a Kuratowski subgraph
              kuratowski_v_ = v;
              kuratowski_x_ = first_side_vertex;
              kuratowski_y_ = second_side_vertex;
              return false;
            }
          }

          // Otherwise, the fact that we didn't find a pertinent
          // vertex on this face is fine - we should set the
          // short-circuit edges and break out of this loop to
          // start looking at a different pertinent root.

          if (first_side_vertex == second_side_vertex) {
            if (first_tail != v) {
              vertex_t first = face_handles_[first_tail].first_vertex();
              vertex_t second = face_handles_[first_tail].second_vertex();
              std::tie(first_side_vertex, first_tail) =
                  std::tuple(first_tail, first == first_side_vertex ? second : first);
            } else if (second_tail != v) {
              vertex_t first = face_handles_[second_tail].first_vertex();
              vertex_t second = face_handles_[second_tail].second_vertex();
              std::tie(second_side_vertex, second_tail) =
                  std::tuple(second_tail, first == second_side_vertex ? second : first);
            } else {
              break;
            }
          }

          canonical_dfs_child_[first_side_vertex] = canonical_dfs_child_[root_face_handle.first_vertex()];
          canonical_dfs_child_[second_side_vertex] = canonical_dfs_child_[root_face_handle.second_vertex()];
          root_face_handle.set_first_vertex(first_side_vertex);
          root_face_handle.set_second_vertex(second_side_vertex);

          if (face_handles_[first_side_vertex].first_vertex() == first_tail) {
            face_handles_[first_side_vertex].set_first_vertex(v);
          } else {
            face_handles_[first_side_vertex].set_second_vertex(v);
          }

          if (face_handles_[second_side_vertex].first_vertex() == second_tail) {
            face_handles_[second_side_vertex].set_first_vertex(v);
          } else {
            face_handles_[second_side_vertex].set_second_vertex(v);
          }

          break;
        }

        // When we unwind the stack, we need to know which direction
        // we came down from on the top face handle

        bool chose_first_lower_path = (chose_first_upper_path && face_handles_[chosen].first_vertex() == first_tail) ||
                                      (!chose_first_upper_path && face_handles_[chosen].first_vertex() == second_tail);

        // If there's a backedge at the chosen vertex, embed it now
        if (backedge_flag_[chosen] == dfs_number_[v]) {
          w = chosen;

          backedge_flag_[chosen] = num_vertices(g_) + 1;
          add_to_merge_points(chosen, StoreOldHandlesPolicy());

          for (auto e : backedges_[chosen]) {
            add_to_embedded_edges(e, StoreOldHandlesPolicy());

            if (chose_first_lower_path) {
              face_handles_[chosen].push_first(e, g_);
            } else {
              face_handles_[chosen].push_second(e, g_);
            }
          }
        } else {
          merge_stack_.emplace_back(chosen, chose_first_upper_path, chose_first_lower_path);
          curr_face_handle = *pertinent_roots_[chosen]->begin();
          continue;
        }

        // Unwind the merge stack to the root, merging all bicomps

        bool bottom_path_follows_first = false;
        bool top_path_follows_first = false;
        bool next_bottom_follows_first = chose_first_upper_path;

        vertex_t merge_point = chosen;

        while (!merge_stack_.empty()) {
          bottom_path_follows_first = next_bottom_follows_first;
          std::tie(merge_point, next_bottom_follows_first, top_path_follows_first) = merge_stack_.back();
          merge_stack_.pop_back();

          face_handle_t top_handle(face_handles_[merge_point]);
          face_handle_t bottom_handle(*pertinent_roots_[merge_point]->begin());

          vertex_t bottom_dfs_child = canonical_dfs_child_[pertinent_roots_[merge_point]->begin()->first_vertex()];

          remove_vertex_from_separated_dfs_child_list(
              canonical_dfs_child_[pertinent_roots_[merge_point]->begin()->first_vertex()]);

          pertinent_roots_[merge_point]->pop_front();

          add_to_merge_points(top_handle.get_anchor(), StoreOldHandlesPolicy());

          if (top_path_follows_first && bottom_path_follows_first) {
            bottom_handle.flip();
            top_handle.glue_first_to_second(bottom_handle);
          } else if (!top_path_follows_first && bottom_path_follows_first) {
            flipped_[bottom_dfs_child] = true;
            top_handle.glue_second_to_first(bottom_handle);
          } else if (top_path_follows_first && !bottom_path_follows_first) {
            flipped_[bottom_dfs_child] = true;
            top_handle.glue_first_to_second(bottom_handle);
          } else  //! top_path_follows_first &&
                  //! !bottom_path_follows_first
          {
            bottom_handle.flip();
            top_handle.glue_second_to_first(bottom_handle);
          }
        }

        // Finally, embed all edges (v,w) at their upper end points
        canonical_dfs_child_[w] = canonical_dfs_child_[root_face_handle.first_vertex()];

        add_to_merge_points(root_face_handle.get_anchor(), StoreOldHandlesPolicy());

        for (auto e : backedges_[chosen]) {
          if (next_bottom_follows_first) {
            root_face_handle.push_first(e, g_);
          } else {
            root_face_handle.push_second(e, g_);
          }
        }

        backedges_[chosen].clear();
        curr_face_handle = root_face_handle;

      }  // while(true)

    }  // while(!pertinent_roots_[v]->empty())

    return true;
  }

  void store_old_face_handles(planar_detail::no_old_handles /*unused*/) {}

  void store_old_face_handles(planar_detail::store_old_handles /*unused*/) {
    for (auto mp : current_merge_points_) {
      face_handles_[mp].store_old_face_handles();
    }
    current_merge_points_.clear();
  }

  void add_to_merge_points(vertex_t /*unused*/, planar_detail::no_old_handles /*unused*/) {}

  void add_to_merge_points(vertex_t v, planar_detail::store_old_handles /*unused*/) {
    current_merge_points_.push_back(v);
  }

  void add_to_embedded_edges(edge_t /*unused*/, planar_detail::no_old_handles /*unused*/) {}

  void add_to_embedded_edges(edge_t e, planar_detail::store_old_handles /*unused*/) { embedded_edges_.push_back(e); }

  void clean_up_embedding(planar_detail::no_embedding /*unused*/) {}

  void clean_up_embedding(planar_detail::store_embedding /*unused*/) {
    // If the graph isn't biconnected, we'll still have entries
    // in the separated_dfs_child_list_ for some vertices. Since
    // these represent articulation points, we can obtain a
    // planar embedding no matter what order we embed them in.

    for (auto x : vertices(g_)) {
      if (separated_dfs_child_list_[x]->empty()) {
        continue;
      }
      for (auto y : *separated_dfs_child_list_[x]) {
        dfs_child_handles_[y].flip();
        face_handles_[x].glue_first_to_second(dfs_child_handles_[y]);
      }
    }

    // Up until this point, we've flipped bicomps lazily by setting
    // flipped[v] to true if the bicomp rooted at v was flipped (the
    // lazy aspect of this flip is that all descendents of that vertex
    // need to have their orientations reversed as well). Now, we
    // traverse the DFS tree by DFS number and perform the actual
    // flipping as needed
    for (auto v : vertices_by_dfs_num_) {
      bool v_flipped = flipped_[v];
      bool p_flipped = flipped_[dfs_parent_[v]];
      if (v_flipped && !p_flipped) {
        face_handles_[v].flip();
      } else if (p_flipped && !v_flipped) {
        face_handles_[v].flip();
        flipped_[v] = true;
      } else {
        flipped_[v] = false;
      }
    }

    // If there are any self-loops in the graph, they were flagged
    // during the walkup, and we should add them to the embedding now.
    // Adding a self loop anywhere in the embedding could never
    // invalidate the embedding, but they would complicate the traversal
    // if they were added during the walkup/walkdown.
    for (auto e : self_loops_) {
      face_handles_[source(e, g_)].push_second(e, g_);
    }
  }

  bool pertinent(vertex_t w, vertex_t v) {
    // w is pertinent with respect to v if there is a backedge (v,w) or if
    // w is the root of a bicomp that contains a pertinent vertex.

    return backedge_flag_[w] == dfs_number_[v] || !pertinent_roots_[w]->empty();
  }

  bool externally_active(vertex_t w, vertex_t v) {
    // Let a be any proper depth-first search ancestor of v. w is externally
    // active with respect to v if there exists a backedge (a,w) or a
    // backedge (a,w_0) for some w_0 in a descendent bicomp of w.

    std::size_t dfs_number_of_v = dfs_number_[v];
    return (least_ancestor_[w] < dfs_number_of_v) ||
           (!separated_dfs_child_list_[w]->empty() &&
            low_point_[separated_dfs_child_list_[w]->front()] < dfs_number_of_v);
  }

  bool internally_active(vertex_t w, vertex_t v) { return pertinent(w, v) && !externally_active(w, v); }

  void remove_vertex_from_separated_dfs_child_list(vertex_t v) {
    auto to_delete = separated_node_in_parent_list_[v];
    garbage_.splice(garbage_.end(), *separated_dfs_child_list_[dfs_parent_[v]], to_delete, std::next(to_delete));
  }

  // End of the implementation of the basic Boyer-Myrvold Algorithm. The rest
  // of the code below implements the isolation of a Kuratowski subgraph in
  // the case that the input graph is not planar. This is by far the most
  // complicated part of the implementation.

 public:
  template <typename EdgeToBoolPropertyMap, typename EdgeContainer>
  vertex_t kuratowski_walkup(vertex_t v, EdgeToBoolPropertyMap forbidden_edge, EdgeToBoolPropertyMap goal_edge,
                             EdgeToBoolPropertyMap is_embedded, EdgeContainer& path_edges) {
    vertex_t current_endpoint;
    bool seen_goal_edge = false;

    for (auto e : out_edges(v, g_)) {
      forbidden_edge[e] = true;
    }

    for (auto e : out_edges(v, g_)) {
      path_edges.clear();

      current_endpoint = target(e, g_) == v ? source(e, g_) : target(e, g_);

      if (dfs_number_[current_endpoint] < dfs_number_[v] || is_embedded[e] || v == current_endpoint) {
        // Not a backedge
        continue;
      }

      path_edges.push_back(e);
      if (goal_edge[e]) {
        return current_endpoint;
      }

      using walkup_itr_t = typename face_edge_iterator<>::type;

      walkup_itr_t walkup_itr(current_endpoint, face_handles_.ref(), first_side());
      walkup_itr_t walkup_end;

      seen_goal_edge = false;

      while (true) {
        if (walkup_itr != walkup_end && forbidden_edge[*walkup_itr]) {
          break;
        }

        while (walkup_itr != walkup_end && !goal_edge[*walkup_itr] && !forbidden_edge[*walkup_itr]) {
          edge_t f(*walkup_itr);
          forbidden_edge[f] = true;
          path_edges.push_back(f);
          current_endpoint = source(f, g_) == current_endpoint ? target(f, g_) : source(f, g_);
          ++walkup_itr;
        }

        if (walkup_itr != walkup_end && goal_edge[*walkup_itr]) {
          path_edges.push_back(*walkup_itr);
          seen_goal_edge = true;
          break;
        }

        walkup_itr = walkup_itr_t(current_endpoint, face_handles_.ref(), first_side());
      }

      if (seen_goal_edge) {
        break;
      }
    }

    if (seen_goal_edge) {
      return current_endpoint;
    }
    return graph_traits<Graph>::null_vertex();
  }

  template <typename OutputIterator, typename EdgeIndexMap>
  void extract_kuratowski_subgraph(OutputIterator o_itr, EdgeIndexMap em) {
    // If the main algorithm has failed to embed one of the back-edges from
    // a vertex v, we can use the current state of the algorithm to isolate
    // a Kuratowksi subgraph. The isolation process breaks down into five
    // cases, A - E. The general configuration of all five cases is shown in
    //                  figure 1. There is a vertex v from which the planar
    //         v        embedding process could not proceed. This means that
    //         |        there exists some bicomp containing three vertices
    //       -----      x,y, and z as shown such that x and y are externally
    //      |     |     active with respect to v (which means that there are
    //      x     y     two vertices x_0 and y_0 such that (1) both x_0 and
    //      |     |     y_0 are proper depth-first search ancestors of v and
    //       --z--      (2) there are two disjoint paths, one connecting x
    //                  and x_0 and one connecting y and y_0, both
    //                  consisting
    //       fig. 1     entirely of unembedded edges). Furthermore, there
    //                  exists a vertex z_0 such that z is a depth-first
    // search ancestor of z_0 and (v,z_0) is an unembedded back-edge from v.
    // x,y and z all exist on the same bicomp, which consists entirely of
    // embedded edges. The five subcases break down as follows, and are
    // handled by the algorithm logically in the order A-E: First, if v is
    // not on the same bicomp as x,y, and z, a K_3_3 can be isolated - this
    // is case A. So, we'll assume that v is on the same bicomp as x,y, and
    // z. If z_0 is on a different bicomp than x,y, and z, a K_3_3 can also
    // be isolated - this is a case B - so we'll assume from now on that v
    // is on the same bicomp as x, y, and z=z_0. In this case, one can use
    // properties of the Boyer-Myrvold algorithm to show the existence of an
    // "x-y path" connecting some vertex on the "left side" of the x,y,z
    // bicomp with some vertex on the "right side" of the bicomp (where the
    // left and right are split by a line drawn through v and z.If either of
    // the endpoints of the x-y path is above x or y on the bicomp, a K_3_3
    // can be isolated - this is a case C. Otherwise, both endpoints are at
    // or below x and y on the bicomp. If there is a vertex alpha on the x-y
    // path such that alpha is not x or y and there's a path from alpha to v
    // that's disjoint from any of the edges on the bicomp and the x-y path,
    // a K_3_3 can be isolated - this is a case D. Otherwise, properties of
    // the Boyer-Myrvold algorithm can be used to show that another vertex
    // w exists on the lower half of the bicomp such that w is externally
    // active with respect to v. w can then be used to isolate a K_5 - this
    // is the configuration of case E.

    // Clear the short-circuit edges - these are needed for the planar
    // testing/embedding algorithm to run in linear time, but they'll
    // complicate the kuratowski subgraph isolation
    for (auto v : vertices(g_)) {
      face_handles_[v].reset_vertex_cache();
      dfs_child_handles_[v].reset_vertex_cache();
    }

    vertex_t v = kuratowski_v_;
    vertex_t x = kuratowski_x_;
    vertex_t y = kuratowski_y_;

    using edge_to_bool_map_t = iterator_property_map<typename std::vector<bool>::iterator, EdgeIndexMap>;

    std::vector<bool> is_in_subgraph_vector(num_edges(g_), false);
    edge_to_bool_map_t is_in_subgraph(is_in_subgraph_vector.begin(), em);

    std::vector<bool> is_embedded_vector(num_edges(g_), false);
    edge_to_bool_map_t is_embedded(is_embedded_vector.begin(), em);

    for (auto emb : embedded_edges_) {
      is_embedded[emb] = true;
    }

    // upper_face_vertex is true for x,y, and all vertices above x and y in
    // the bicomp
    vertex_to_bool_map_t upper_face_vertex(num_vertices(g_), vm_, false);
    vertex_to_bool_map_t lower_face_vertex(num_vertices(g_), vm_, false);

    // These next few variable declarations are all things that we need
    // to find.
    vertex_t z = graph_traits<Graph>::null_vertex();
    vertex_t bicomp_root;
    vertex_t w = graph_traits<Graph>::null_vertex();
    face_handle_t w_handle;
    face_handle_t v_dfchild_handle;
    vertex_t first_x_y_path_endpoint = graph_traits<Graph>::null_vertex();
    vertex_t second_x_y_path_endpoint = graph_traits<Graph>::null_vertex();
    vertex_t w_ancestor = v;

    bm_case_t chosen_case = bm_case_t::bm_no_case_chosen;

    std::vector<edge_t> x_external_path;
    std::vector<edge_t> y_external_path;
    std::vector<edge_t> case_d_edges;

    std::vector<edge_t> z_v_path;
    std::vector<edge_t> w_path;

    // first, use a walkup to find a path from V that starts with a
    // backedge from V, then goes up until it hits either X or Y
    //(but doesn't find X or Y as the root of a bicomp)

    typename face_vertex_iterator<>::type x_upper_itr(x, face_handles_.ref(), first_side());
    typename face_vertex_iterator<>::type x_lower_itr(x, face_handles_.ref(), second_side());
    typename face_vertex_iterator<>::type face_itr;
    typename face_vertex_iterator<>::type face_end;

    // Don't know which path from x is the upper or lower path -
    // we'll find out here
    for (face_itr = x_upper_itr; face_itr != face_end; ++face_itr) {
      if (*face_itr == y) {
        std::swap(x_upper_itr, x_lower_itr);
        break;
      }
    }

    upper_face_vertex[x] = true;

    vertex_t current_vertex = x;
    vertex_t previous_vertex;
    for (face_itr = x_upper_itr; face_itr != face_end; ++face_itr) {
      previous_vertex = current_vertex;
      current_vertex = *face_itr;
      upper_face_vertex[current_vertex] = true;
    }

    v_dfchild_handle = dfs_child_handles_[canonical_dfs_child_[previous_vertex]];

    for (face_itr = x_lower_itr; *face_itr != y; ++face_itr) {
      vertex_t u(*face_itr);
      lower_face_vertex[u] = true;

      // haven't found a w yet
      if (w == graph_traits<Graph>::null_vertex()) {
        for (auto& root : *pertinent_roots_[u]) {
          if (low_point_[canonical_dfs_child_[root.first_vertex()]] < dfs_number_[u]) {
            w = u;
            w_handle = root;
            break;
          }
        }
      }
    }

    for (; face_itr != face_end; ++face_itr) {
      vertex_t u(*face_itr);
      upper_face_vertex[u] = true;
      bicomp_root = u;
    }

    using walkup_itr_t = typename face_edge_iterator<>::type;

    std::vector<bool> outer_face_edge_vector(num_edges(g_), false);
    edge_to_bool_map_t outer_face_edge(outer_face_edge_vector.begin(), em);

    walkup_itr_t walkup_end;
    for (walkup_itr_t walkup_itr(x, face_handles_.ref(), first_side()); walkup_itr != walkup_end; ++walkup_itr) {
      outer_face_edge[*walkup_itr] = true;
      is_in_subgraph[*walkup_itr] = true;
    }

    for (walkup_itr_t walkup_itr(x, face_handles_.ref(), second_side()); walkup_itr != walkup_end; ++walkup_itr) {
      outer_face_edge[*walkup_itr] = true;
      is_in_subgraph[*walkup_itr] = true;
    }

    std::vector<bool> forbidden_edge_vector(num_edges(g_), false);
    edge_to_bool_map_t forbidden_edge(forbidden_edge_vector.begin(), em);

    std::vector<bool> goal_edge_vector(num_edges(g_), false);
    edge_to_bool_map_t goal_edge(goal_edge_vector.begin(), em);

    // Find external path to x and to y

    for (auto e : edges(g_)) {
      goal_edge[e] = !outer_face_edge[e] && (source(e, g_) == x || target(e, g_) == x);
      forbidden_edge[e] = outer_face_edge[e];
    }

    vertex_t x_ancestor = v;
    vertex_t x_endpoint = graph_traits<Graph>::null_vertex();

    while (x_endpoint == graph_traits<Graph>::null_vertex()) {
      x_ancestor = dfs_parent_[x_ancestor];
      x_endpoint = kuratowski_walkup(x_ancestor, forbidden_edge, goal_edge, is_embedded, x_external_path);
    }

    for (auto e : edges(g_)) {
      goal_edge[e] = !outer_face_edge[e] && (source(e, g_) == y || target(e, g_) == y);
      forbidden_edge[e] = outer_face_edge[e];
    }

    vertex_t y_ancestor = v;
    vertex_t y_endpoint = graph_traits<Graph>::null_vertex();

    while (y_endpoint == graph_traits<Graph>::null_vertex()) {
      y_ancestor = dfs_parent_[y_ancestor];
      y_endpoint = kuratowski_walkup(y_ancestor, forbidden_edge, goal_edge, is_embedded, y_external_path);
    }

    vertex_t parent;
    vertex_t child;

    // If v isn't on the same bicomp as x and y, it's a case A
    if (bicomp_root != v) {
      chosen_case = bm_case_t::bm_case_a;

      for (auto u : vertices(g_)) {
        if (!lower_face_vertex[u]) {
          continue;
        }
        for (auto e : out_edges(u, g_)) {
          if (!outer_face_edge[e]) {
            goal_edge[e] = true;
          }
        }
      }

      for (auto e : edges(g_)) {
        forbidden_edge[e] = outer_face_edge[e];
      }

      z = kuratowski_walkup(v, forbidden_edge, goal_edge, is_embedded, z_v_path);
    } else if (w != graph_traits<Graph>::null_vertex()) {
      chosen_case = bm_case_t::bm_case_b;

      for (auto e : edges(g_)) {
        goal_edge[e] = false;
        forbidden_edge[e] = outer_face_edge[e];
      }

      goal_edge[w_handle.first_edge()] = true;
      goal_edge[w_handle.second_edge()] = true;

      z = kuratowski_walkup(v, forbidden_edge, goal_edge, is_embedded, z_v_path);

      for (auto e : edges(g_)) {
        forbidden_edge[e] = outer_face_edge[e];
      }

      for (auto e : z_v_path) {
        goal_edge[e] = true;
      }

      w_ancestor = v;
      vertex_t w_endpoint = graph_traits<Graph>::null_vertex();

      while (w_endpoint == graph_traits<Graph>::null_vertex()) {
        w_ancestor = dfs_parent_[w_ancestor];
        w_endpoint = kuratowski_walkup(w_ancestor, forbidden_edge, goal_edge, is_embedded, w_path);
      }

      // We really want both the w walkup and the z walkup to finish on
      // exactly the same edge, but for convenience (since we don't have
      // control over which side of a bicomp a walkup moves up) we've
      // defined the walkup to either end at w_handle.first_edge() or
      // w_handle.second_edge(). If both walkups ended at different edges,
      // we'll do a little surgery on the w walkup path to make it follow
      // the other side of the final bicomp.

      if ((w_path.back() == w_handle.first_edge() && z_v_path.back() == w_handle.second_edge()) ||
          (w_path.back() == w_handle.second_edge() && z_v_path.back() == w_handle.first_edge())) {
        walkup_itr_t wi;
        walkup_itr_t wi_end;
        edge_t final_edge = w_path.back();
        vertex_t anchor =
            source(final_edge, g_) == w_handle.get_anchor() ? target(final_edge, g_) : source(final_edge, g_);
        if (face_handles_[anchor].first_edge() == final_edge) {
          wi = walkup_itr_t(anchor, face_handles_.ref(), second_side());
        } else {
          wi = walkup_itr_t(anchor, face_handles_.ref(), first_side());
        }

        w_path.pop_back();

        for (; wi != wi_end; ++wi) {
          edge_t e(*wi);
          if (w_path.back() == e) {
            w_path.pop_back();
          } else {
            w_path.push_back(e);
          }
        }
      }
    } else {
      // We need to find a valid z, since the x-y path re-defines the
      // lower face, and the z we found earlier may now be on the upper
      // face.

      chosen_case = bm_case_t::bm_case_e;

      // The z we've used so far is just an externally active vertex on
      // the lower face path, but may not be the z we need for a case C,
      // D, or E subgraph. the z we need now is any externally active
      // vertex on the lower face path with both old_face_handles edges on
      // the outer face. Since we know an x-y path exists, such a z must
      // also exist.

      // TODO: find this z in the first place.

      // find the new z

      for (face_itr = x_lower_itr; *face_itr != y; ++face_itr) {
        vertex_t possible_z(*face_itr);
        if (pertinent(possible_z, v) && outer_face_edge[face_handles_[possible_z].old_first_edge()] &&
            outer_face_edge[face_handles_[possible_z].old_second_edge()]) {
          z = possible_z;
          break;
        }
      }

      // find x-y path, and a w if one exists.

      if (externally_active(z, v)) {
        w = z;
      }

      using old_face_iterator_t = typename face_edge_iterator<single_side, previous_iteration>::type;

      old_face_iterator_t first_old_face_itr(z, face_handles_.ref(), first_side());
      old_face_iterator_t second_old_face_itr(z, face_handles_.ref(), second_side());
      old_face_iterator_t old_face_end;

      std::vector<old_face_iterator_t> old_face_iterators;
      old_face_iterators.push_back(first_old_face_itr);
      old_face_iterators.push_back(second_old_face_itr);

      vertex_to_bool_map_t x_y_path_vertex(num_vertices(g_), vm_, false);

      for (old_face_iterator_t old_face_itr : old_face_iterators) {
        vertex_t prev_v;
        bool seen_x_or_y = false;
        vertex_t curr_v = z;
        for (; old_face_itr != old_face_end; ++old_face_itr) {
          edge_t e(*old_face_itr);
          prev_v = curr_v;
          curr_v = source(e, g_) == curr_v ? target(e, g_) : source(e, g_);

          if (curr_v == x || curr_v == y) {
            seen_x_or_y = true;
          }

          if (w == graph_traits<Graph>::null_vertex() && externally_active(curr_v, v) && outer_face_edge[e] &&
              outer_face_edge[*std::next(old_face_itr)] && !seen_x_or_y) {
            w = curr_v;
          }

          if (!outer_face_edge[e]) {
            if (!upper_face_vertex[curr_v] && !lower_face_vertex[curr_v]) {
              x_y_path_vertex[curr_v] = true;
            }

            is_in_subgraph[e] = true;
            if (upper_face_vertex[source(e, g_)] || lower_face_vertex[source(e, g_)]) {
              if (first_x_y_path_endpoint == graph_traits<Graph>::null_vertex()) {
                first_x_y_path_endpoint = source(e, g_);
              } else {
                second_x_y_path_endpoint = source(e, g_);
              }
            }
            if (upper_face_vertex[target(e, g_)] || lower_face_vertex[target(e, g_)]) {
              if (first_x_y_path_endpoint == graph_traits<Graph>::null_vertex()) {
                first_x_y_path_endpoint = target(e, g_);
              } else {
                second_x_y_path_endpoint = target(e, g_);
              }
            }
          } else if (prev_v == x || prev_v == y) {
            chosen_case = bm_case_t::bm_case_c;
          }
        }
      }

      // Look for a case D - one of v's embedded edges will connect to the
      // x-y path along an inner face path.

      // First, get a list of all of v's embedded child edges

      for (auto embedded_edge : out_edges(v, g_)) {
        if (!is_embedded[embedded_edge] || embedded_edge == dfs_parent_edge_[v]) {
          continue;
        }

        case_d_edges.push_back(embedded_edge);

        vertex_t curr_v = source(embedded_edge, g_) == v ? target(embedded_edge, g_) : source(embedded_edge, g_);

        typename face_edge_iterator<>::type internal_face_itr;
        typename face_edge_iterator<>::type internal_face_end;
        if (face_handles_[curr_v].first_vertex() == v) {
          internal_face_itr = typename face_edge_iterator<>::type(curr_v, face_handles_.ref(), second_side());
        } else {
          internal_face_itr = typename face_edge_iterator<>::type(curr_v, face_handles_.ref(), first_side());
        }

        while (internal_face_itr != internal_face_end && !outer_face_edge[*internal_face_itr] &&
               !x_y_path_vertex[curr_v]) {
          edge_t e(*internal_face_itr);
          case_d_edges.push_back(e);
          curr_v = source(e, g_) == curr_v ? target(e, g_) : source(e, g_);
          ++internal_face_itr;
        }

        if (x_y_path_vertex[curr_v]) {
          chosen_case = bm_case_t::bm_case_d;
          break;
        }
        case_d_edges.clear();
      }
    }

    if (chosen_case != bm_case_t::bm_case_b && chosen_case != bm_case_t::bm_case_a) {
      // Finding z and w.

      for (auto e : edges(g_)) {
        goal_edge[e] = !outer_face_edge[e] && (source(e, g_) == z || target(e, g_) == z);
        forbidden_edge[e] = outer_face_edge[e];
      }

      kuratowski_walkup(v, forbidden_edge, goal_edge, is_embedded, z_v_path);

      if (chosen_case == bm_case_t::bm_case_e) {
        for (auto e : edges(g_)) {
          forbidden_edge[e] = outer_face_edge[e];
          goal_edge[e] = !outer_face_edge[e] && (source(e, g_) == w || target(e, g_) == w);
        }

        for (auto e : out_edges(w, g_)) {
          if (!outer_face_edge[e]) {
            goal_edge[e] = true;
          }
        }

        for (auto e : z_v_path) {
          goal_edge[e] = true;
        }

        w_ancestor = v;
        vertex_t w_endpoint = graph_traits<Graph>::null_vertex();

        while (w_endpoint == graph_traits<Graph>::null_vertex()) {
          w_ancestor = dfs_parent_[w_ancestor];
          w_endpoint = kuratowski_walkup(w_ancestor, forbidden_edge, goal_edge, is_embedded, w_path);
        }
      }
    }

    // We're done isolating the Kuratowski subgraph at this point -
    // but there's still some cleaning up to do.

    // Update is_in_subgraph with the paths we just found
    for (auto e : x_external_path) {
      is_in_subgraph[e] = true;
    }

    for (auto e : y_external_path) {
      is_in_subgraph[e] = true;
    }

    for (auto e : z_v_path) {
      is_in_subgraph[e] = true;
    }

    for (auto e : case_d_edges) {
      is_in_subgraph[e] = true;
    }

    for (auto e : w_path) {
      is_in_subgraph[e] = true;
    }

    child = bicomp_root;
    parent = dfs_parent_[child];
    while (child != parent) {
      is_in_subgraph[dfs_parent_edge_[child]] = true;
      std::tie(parent, child) = std::make_pair(dfs_parent_[parent], parent);
    }

    // At this point, we've already isolated the Kuratowski subgraph and
    // collected all of the edges that compose it in the is_in_subgraph
    // property map. But we want the verification of such a subgraph to be
    // a deterministic process, and we can simplify the function
    // is_kuratowski_subgraph by cleaning up some edges here.

    if (chosen_case == bm_case_t::bm_case_b) {
      is_in_subgraph[dfs_parent_edge_[v]] = false;
    } else if (chosen_case == bm_case_t::bm_case_c) {
      // In a case C subgraph, at least one of the x-y path endpoints
      // (call it alpha) is above either x or y on the outer face. The
      // other endpoint may be attached at x or y OR above OR below. In
      // any of these three cases, we can form a K_3_3 by removing the
      // edge attached to v on the outer face that is NOT on the path to
      // alpha.

      typename face_vertex_iterator<single_side, follow_visitor>::type follow_itr;
      if (face_handles_[v_dfchild_handle.first_vertex()].first_edge() == v_dfchild_handle.first_edge()) {
        follow_itr = typename face_vertex_iterator<single_side, follow_visitor>::type(
            v_dfchild_handle.first_vertex(), face_handles_.ref(), second_side());
      } else {
        follow_itr = typename face_vertex_iterator<single_side, follow_visitor>::type(
            v_dfchild_handle.first_vertex(), face_handles_.ref(), first_side());
      }

      for (; true; ++follow_itr) {
        vertex_t curr_v(*follow_itr);
        if (curr_v == x || curr_v == y) {
          is_in_subgraph[v_dfchild_handle.first_edge()] = false;
          break;
        }
        if (curr_v == first_x_y_path_endpoint || curr_v == second_x_y_path_endpoint) {
          is_in_subgraph[v_dfchild_handle.second_edge()] = false;
          break;
        }
      }
    } else if (chosen_case == bm_case_t::bm_case_d) {
      // Need to remove both of the edges adjacent to v on the outer face.
      // remove the connecting edges from v to bicomp, then
      // is_kuratowski_subgraph will shrink vertices of degree 1
      // automatically...

      is_in_subgraph[v_dfchild_handle.first_edge()] = false;
      is_in_subgraph[v_dfchild_handle.second_edge()] = false;
    } else if (chosen_case == bm_case_t::bm_case_e) {
      // Similarly to case C, if the endpoints of the x-y path are both
      // below x and y, we should remove an edge to allow the subgraph to
      // contract to a K_3_3.

      if ((first_x_y_path_endpoint != x && first_x_y_path_endpoint != y) ||
          (second_x_y_path_endpoint != x && second_x_y_path_endpoint != y)) {
        is_in_subgraph[dfs_parent_edge_[v]] = false;

        vertex_t deletion_endpoint;
        vertex_t other_endpoint;
        if (lower_face_vertex[first_x_y_path_endpoint]) {
          deletion_endpoint = second_x_y_path_endpoint;
          other_endpoint = first_x_y_path_endpoint;
        } else {
          deletion_endpoint = first_x_y_path_endpoint;
          other_endpoint = second_x_y_path_endpoint;
        }

        typename face_edge_iterator<>::type del_face_itr;
        typename face_edge_iterator<>::type del_face_end;

        bool found_other_endpoint = false;
        for (del_face_itr = typename face_edge_iterator<>::type(deletion_endpoint, face_handles_.ref(), first_side());
             del_face_itr != del_face_end; ++del_face_itr) {
          edge_t e(*del_face_itr);
          if (source(e, g_) == other_endpoint || target(e, g_) == other_endpoint) {
            found_other_endpoint = true;
            break;
          }
        }

        if (found_other_endpoint) {
          is_in_subgraph[face_handles_[deletion_endpoint].first_edge()] = false;
        } else {
          is_in_subgraph[face_handles_[deletion_endpoint].second_edge()] = false;
        }
      }
    }

    for (auto e : edges(g_)) {
      if (is_in_subgraph[e]) {
        *o_itr = e;
      }
    }
  }

  template <typename EdgePermutation>
  void make_edge_permutation(EdgePermutation perm) {
    for (auto v : vertices(g_)) {
      perm[v].clear();
      face_handles_[v].get_list(std::back_inserter(perm[v]));
    }
  }

 private:
  const Graph& g_;
  VertexIndexMap vm_;

  vertex_t kuratowski_v_;
  vertex_t kuratowski_x_;
  vertex_t kuratowski_y_;

  vertex_list_t garbage_;  // we delete items from linked lists by
                           // splicing them into garbage_

  // only need these two for kuratowski subgraph isolation
  std::vector<vertex_t> current_merge_points_;
  std::vector<edge_t> embedded_edges_;

  // property map storage
  edge_vector_t self_loops_;
  vertex_vector_t vertices_by_dfs_num_;

  // property maps
  vertex_to_v_size_map_t low_point_;
  vertex_to_vertex_map_t dfs_parent_;
  vertex_to_v_size_map_t dfs_number_;
  vertex_to_v_size_map_t least_ancestor_;
  vertex_to_face_handle_list_ptr_map_t pertinent_roots_;
  vertex_to_v_size_map_t backedge_flag_;
  vertex_to_v_size_map_t visited_;
  vertex_to_face_handle_map_t face_handles_;
  vertex_to_face_handle_map_t dfs_child_handles_;
  vertex_to_vertex_list_ptr_map_t separated_dfs_child_list_;
  vertex_to_separated_node_map_t separated_node_in_parent_list_;
  vertex_to_vertex_map_t canonical_dfs_child_;
  vertex_to_bool_map_t flipped_;
  vertex_to_edge_vector_map_t backedges_;
  vertex_to_edge_map_t dfs_parent_edge_;  // only need for kuratowski

  merge_stack_t merge_stack_;
};

}  // namespace bagl

#endif  // BAGL_BAGL_PLANAR_DETAIL_BOYER_MYRVOLD_IMPL_H_
