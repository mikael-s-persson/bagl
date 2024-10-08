// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Jeremy G. Siek, Andrew Lumsdaine, Lie-Quan Lee
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_READ_DIMACS_H_
#define BAGL_BAGL_READ_DIMACS_H_

#include <array>
#include <charconv>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"

namespace bagl {

// Reads maximal flow problem in extended DIMACS format.
// This works, but could use some polishing.

namespace dimacs_detail {

template <class Graph, class CapacityMap, class ReverseEdgeMap>
bool read_dimacs_max_flow_internal(Graph& g, CapacityMap capacity, ReverseEdgeMap reverse_edge,
                                   graph_vertex_descriptor_t<Graph>& src, graph_vertex_descriptor_t<Graph>& sink,
                                   std::istream& in, bool require_source_and_sink, const std::string& problem_type) {
  using vertex_descriptor = graph_vertex_descriptor_t<Graph>;

  std::vector<vertex_descriptor> verts;

  // number of edges and nodes
  int m = 0;
  int n = 0;
  int i = 0;
  int head = 0;
  int tail = 0;
  int cap = 0;

  int no_lines = 0;   /* no of current input line */
  int no_plines = 0;  /* no of problem-lines */
  int no_nslines = 0; /* no of node-source-lines */
  int no_nklines = 0; /* no of node-source-lines */
  int no_alines = 0;  /* no of arc-lines */

  // -------------- error messages ----------------
  constexpr std::array<std::string_view, 20> err_message = {
      /* 0*/ "more than one problem line.",
      /* 1*/ "wrong number of parameters in the problem line.",
      /* 2*/ "it is not a Max Flow problem line.",
      /* 3*/ "bad value of a parameter in the problem line.",
      /* 4*/ "can't obtain enough memory to solve this problem.",
      /* 5*/ "more than one line with the problem name.",
      /* 6*/ "can't read problem name.",
      /* 7*/ "problem description must be before node description.",
      /* 8*/ "this parser doesn't support multiply sources and sinks.",
      /* 9*/ "wrong number of parameters in the node line.",
      /*10*/ "wrong value of parameters in the node line.",
      /*11*/ "source and sink descriptions must be before arc descriptions.",
      /*12*/ "too many arcs in the input.",
      /*13*/ "wrong number of parameters in the arc line.",
      /*14*/ "wrong value of parameters in the arc line.",
      /*15*/ "unknown line type in the input.",
      /*16*/ "reading error.",
      /*17*/ "not enough arcs in the input.",
      /*18*/ "source or sink doesn't have incident arcs.",
      /*19*/ "can't read anything from the input file.",
  };
  /* --------------------------------------------------------------- */

  // The main loop:
  // -  reads the line of the input,
  // -  analyses its type,
  // -  checks correctness of parameters,
  // -  puts data to the arrays,
  // -  does service functions

  std::vector<std::string_view> fields;
  std::string in_line;
  while (std::getline(in, in_line)) {
    ++no_lines;

    fields.clear();
    std::string::size_type sp_last_pos = 0;
    auto sp_pos = in_line.find_first_of(' ');
    while (sp_pos != std::string::npos) {
      fields.emplace_back(&in_line[sp_last_pos], sp_pos - sp_last_pos);
      sp_last_pos = sp_pos + 1;
      sp_pos = in_line.find_first_of(' ', sp_last_pos);
    }
    fields.emplace_back(&in_line[sp_last_pos], in_line.size() - sp_last_pos);

    switch (fields[0][0]) {
      case 'c':  // skip lines with comments
        [[fallthrough]];
      case '\n':  // skip empty lines
        [[fallthrough]];
      case '\0':  // skip empty lines at the end of file
        break;

      case 'p':
        // problem description

        if (no_plines > 0) {
          // more than one problem line
          std::cerr << "line " << no_lines << " of input - " << err_message[0] << std::endl;
          return false;
        }
        no_plines = 1;

        // reading problem line: type of problem, no of nodes, no of arcs
        if (fields.size() != 4) {
          // wrong number of parameters in the problem line
          std::cerr << "line " << no_lines << " of input - " << err_message[1] << std::endl;
          return false;
        }

        if (fields[1] != problem_type) {
          // wrong problem type
          std::cerr << "line " << no_lines << " of input - " << err_message[2] << std::endl;
          return false;
        }

        if (std::from_chars(fields[2].data(), fields[2].data() + fields[2].size(), n).ec != std::errc{} ||
            std::from_chars(fields[3].data(), fields[3].data() + fields[3].size(), m).ec != std::errc{} || n <= 0 ||
            m <= 0) {
          // wrong value of no of arcs or nodes
          std::cerr << "line " << no_lines << " of input - " << err_message[3] << std::endl;
          return false;
        }

        for (int vi = 0; vi < n; ++vi) {
          verts.push_back(add_vertex(g));
        }
        break;

      case 'n':
        // source(s) description

        if (no_plines == 0) {
          // there was not problem line above
          std::cerr << "line " << no_lines << " of input - " << err_message[7] << std::endl;
          return false;
        }

        if (fields.size() != 3) {
          // node line is incorrect
          std::cerr << "line " << no_lines << " of input - " << err_message[9] << std::endl;
          return false;
        }

        if (std::from_chars(fields[1].data(), fields[1].data() + fields[1].size(), i).ec != std::errc{} || --i < 0 ||
            i > n) {
          std::cerr << "line " << no_lines << " of input - " << err_message[10] << std::endl;
          return false;
        }

        switch (fields[2][0]) {
          case 's':
            // source line

            if (no_nslines != 0) {
              // more than one source line
              std::cerr << "line " << no_lines << " of input - " << err_message[8] << std::endl;
              return false;
            }

            no_nslines = 1;
            src = verts[i];
            break;

          case 't':
            // sink line

            if (no_nklines != 0) {
              // more than one sink line
              std::cerr << "line " << no_lines << " of input - " << err_message[8] << std::endl;
              return false;
            }

            no_nklines = 1;
            sink = verts[i];
            break;

          default:
            // wrong type of node-line
            std::cerr << "line " << no_lines << " of input - " << err_message[10] << std::endl;
            return false;
        }
        break;

      case 'a':
        // arc description

        if (require_source_and_sink && (no_nslines == 0 || no_nklines == 0)) {
          // there was not source and sink description above
          std::cerr << "line " << no_lines << " of input - " << err_message[11] << std::endl;
          return false;
        }

        if (no_alines >= m) {
          // too many arcs on input
          std::cerr << "line " << no_lines << " of input - " << err_message[12] << std::endl;
          return false;
        }

        if (fields.size() != 4) {
          // arc description is not correct
          std::cerr << "line " << no_lines << " of input - " << err_message[13] << std::endl;
          return false;
        }

        if (std::from_chars(fields[1].data(), fields[1].data() + fields[1].size(), tail).ec != std::errc{} ||
            std::from_chars(fields[1].data(), fields[1].data() + fields[1].size(), head).ec != std::errc{} ||
            std::from_chars(fields[1].data(), fields[1].data() + fields[1].size(), cap).ec != std::errc{} ||
            --tail < 0 || tail > n || --head < 0 || head > n) {
          // wrong value of nodes
          std::cerr << "line " << no_lines << " of input - " << err_message[14] << std::endl;
          return false;
        }

        {
          auto [e1, in1] = add_edge(verts[tail], verts[head], g);
          auto [e2, in2] = add_edge(verts[head], verts[tail], g);
          if (!in1 || !in2) {
            std::cerr << "unable to add edge (" << head << "," << tail << ")" << std::endl;
            return false;
          }
          capacity[e1] = cap;
          capacity[e2] = 0;
          reverse_edge[e1] = e2;
          reverse_edge[e2] = e1;
        }
        ++no_alines;
        break;

      default:
        // unknown type of line
        std::cerr << "line " << no_lines << " of input - " << err_message[15] << std::endl;
        return false;

    }  // end of switch
  }    // end of input loop

  // ----- all is red  or  error while reading -----

  if (in.eof() == 0) {
    // reading error
    std::cerr << "line " << no_lines << " of input - " << err_message[16] << std::endl;
    return false;
  }

  if (no_lines == 0) {
    // empty input
    std::cerr << "line " << no_lines << " of input - " << err_message[19] << std::endl;
    return false;
  }

  if (no_alines < m) {
    // not enough arcs
    std::cerr << "line " << no_lines << " of input - " << err_message[17] << std::endl;
    return false;
  }

  if (require_source_and_sink && (out_degree(src, g) == 0 || out_degree(sink, g) == 0)) {
    // no arc goes out of the source
    std::cerr << "line " << no_lines << " of input - " << err_message[18] << std::endl;
    return false;
  }

  // Thanks God! all is done
  return true;
}

}  // namespace dimacs_detail

template <concepts::MutableGraph Graph, concepts::MutableLvalueEdgePropertyMap<Graph> CapacityMap,
          concepts::MutableLvalueEdgePropertyMap<Graph> ReverseEdgeMap>
bool read_dimacs_max_flow(Graph& g, CapacityMap capacity, ReverseEdgeMap reverse_edge,
                          graph_vertex_descriptor_t<Graph>& src, graph_vertex_descriptor_t<Graph>& sink,
                          std::istream& in = std::cin) {
  return dimacs_detail::read_dimacs_max_flow_internal(g, capacity, reverse_edge, src, sink, in, true, "max");
}

template <concepts::MutableGraph Graph, concepts::MutableLvalueEdgePropertyMap<Graph> CapacityMap,
          concepts::MutableLvalueEdgePropertyMap<Graph> ReverseEdgeMap>
bool read_dimacs_min_cut(Graph& g, CapacityMap capacity, ReverseEdgeMap reverse_edge, std::istream& in = std::cin) {
  graph_vertex_descriptor_t<Graph> dummy_src;
  graph_vertex_descriptor_t<Graph> dummy_sink;  // Not filled in
  return dimacs_detail::read_dimacs_max_flow_internal(g, capacity, reverse_edge, dummy_src, dummy_sink, in, false,
                                                      "cut");
}

}  // namespace bagl

#endif  // BAGL_BAGL_READ_DIMACS_H_
