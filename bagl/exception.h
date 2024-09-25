// Copyright 2002 Indiana University.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_EXCEPTION_H_
#define BAGL_BAGL_EXCEPTION_H_

#include <stdexcept>
#include <string>
#include <utility>

namespace bagl {

struct bad_graph : std::invalid_argument {
  explicit bad_graph(const std::string& what_arg) : std::invalid_argument(what_arg) {}
};

struct not_a_dag : bad_graph {
  not_a_dag() : bad_graph("The graph must be a DAG.") {}
};

struct negative_edge : bad_graph {
  negative_edge() : bad_graph("The graph may not contain an edge with negative weight.") {}
};

struct negative_cycle : bad_graph {
  negative_cycle() : bad_graph("The graph may not contain negative cycles.") {}
};

struct not_connected : bad_graph {
  not_connected() : bad_graph("The graph must be connected.") {}
};

struct not_complete : bad_graph {
  not_complete() : bad_graph("The graph must be complete.") {}
};

struct graph_exception : std::exception {
  [[nodiscard]] const char* what() const noexcept override = 0;
};

struct bad_parallel_edge : graph_exception {
  std::string from;
  std::string to;
  mutable std::string statement;
  bad_parallel_edge(std::string i, std::string j) : from(std::move(i)), to(std::move(j)) {}

  [[nodiscard]] const char* what() const noexcept override {
    if (statement.empty()) {
      statement = std::string("Failed to add parallel edge: (") + from + "," + to + ")\n";
    }

    return statement.c_str();
  }
};

struct directed_graph_error : graph_exception {
  [[nodiscard]] const char* what() const noexcept override {
    return "read_graphviz: Tried to read a directed graph into an undirected graph.";
  }
};

struct undirected_graph_error : graph_exception {
  [[nodiscard]] const char* what() const noexcept override {
    return "read_graphviz: Tried to read an undirected graph into a directed graph.";
  }
};

}  // namespace bagl

#endif  // BAGL_BAGL_EXCEPTION_H_
