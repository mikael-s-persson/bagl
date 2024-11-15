# Born Again Graph Library (BAGL)

The Born Again Graph Library (BAGL) is a port of the
[Boost Graph Library](https://www.boost.org/doc/libs/release/libs/graph/doc/index.html) to C++20
and decoupled from Boost. This library provides an extensive suite of generic graph
algorithms, interfaces and data structures. As mostly a direct port, the documentation of
the [Boost Graph Library](https://www.boost.org/doc/libs/release/libs/graph/doc/index.html) can
be used as documentation of BAGL, with some common changes noted below.

## Building

This library uses the [Bazel](https://bazel.build/) build system and relies on
the [toolchains_llvm](https://github.com/bazel-contrib/toolchains_llvm) to provide a portable
access to a modern toolchain with sufficient C++20 support (LLVM/Clang/libc++ at version 17 or
greater). Building the library works as usual with Bazel:

```bash
$ # Building the library:
$ bazel build //...:all
$ # Running the test suite:
$ bazel test //...:all
```

In addition to the usual compilation modes (e.g., `-c opt`, `-c dbg` and `-c fastbuild`),
configurations are also provided for sanitizer builds, i.e., `--config=asan` and `--config=msan`.

The availability of the LLVM toolchain for a given host platform can be sparse, which means that
one might need to change the version, see the documentation
of [toolchains_llvm](https://github.com/bazel-contrib/toolchains_llvm) for instructions, and
especially the current [list of supported distributions](https://github.com/bazel-contrib/toolchains_llvm/blob/master/toolchain/internal/llvm_distributions.bzl).

Importing this library into your own Bazel project should be as easy as any other native
Bazel module. For importing it into another build system, then it could be as straightforward
as just pointing it to the include path, since most of BAGL is header-only and has no external
dependency. GraphML and GraphViz support do require building the corresponding libraries (and
GraphML uses the tinyxml2 library), and for those, integration steps would depend on the build
system you use, but given how simple these libraries are, it should be trivial to do.

## Motivation

The BGL is a classic (aka, old) graph library that provides an extensive suite of algorithms based
on very well-designed generic interfaces, aka concepts. However, given its complexity and age,
the library has largely become "abandonware". Moreover, its popularity and widespread use, and
support for older C++ versions, also means that most of its components are effectively frozen to
avoid breaking user-code. And finally, as one of the oldest Boost libraries, the BGL has accumulated
a large number of interdependencies on other Boost library components (over 40 Boost dependencies),
most of which have since become part of the C++ Standard Library in one form or another.

Furthermore, C++20 has introduced, among other things, two major new
features: [Ranges](https://en.cppreference.com/w/cpp/ranges)
and [Concepts](https://en.cppreference.com/w/cpp/concepts). The BGL was one of the pioneers
in making extensive use of both of those generic design patterns, of course, in their pre-standard
Boost forms, that is, Boost.Range and Boost.ConceptCheck. Especially for the latter, the integration
of concepts as a core language feature in C++20 is a major leap forward compared to its legacy Boost
counterpart in terms of syntax, readability, compilation, and expressiveness. In short, at this
juncture, it seemed like the perfect time to bring this library into the modern era.

In addition, some personal motivations include wanting to gain experience with the new standard
features of C++20, and it's a labor of love for a library that shouldn't be left to rot.

## Differences compared to the Boost Graph Library

There are many differences between BAGL and the original BGL, mostly modernizing
and decoupling from the Boost "monolith". However, there are also a few other changes or improvements.
It's not practical to list all specific changes, but here is a list of general differences and
specific improvements:

 - The code lives in the `bagl` namespace, rather than `boost` or `boost::graph`.
 - Many modernization changes have been made since the BGL is effectively a C++03 library, which
   means that, in BAGL, a lot of code has been simplified or improved by features such as move
   semantics, return type deductions, `auto`, `<type_traits>`, range-based for-loops, etc.
    - For what its worth, BAGL has about 25% fewer lines of library code and about 50% fewer lines
      of test code, while having (modestly) increased the feature set and test suite.
 - Concepts in the BGL are named something like `boost::FooConcept` and use the Boost.ConceptCheck
   library. In BAGL, the corresponding concepts are named `bagl::concepts::Foo` and are native C++
   `concept` types (are concepts types?).
    - The concepts, as documented in BGL, have been semantically preserved as much as
      possible but some alterations have been made, mostly to make use of move semantics and perfect
      forwarding, for example, in BGL, adding a vertex with properties to a graph would be `add_vertex(vprop, g)`,
      while BAGL uses an "emplace" style `add_vertex(g, vp_args...)` where the underlying
      vertex property is constructed in-place with perfectly forwarded variadic arguments.
    - Most, if not all, uses of overloading and SFINAE in BGL have been replaced with concept-based
      constraints which not only improves compilation times, but also allows for richer overload
      sets. However, note that this is difficult to achieve in practice and may lead to unintended
      ambiguous overloads or substitution failures (that are errors, but shouldn't be).
 - Many BGL interfaces rely on iterator pairs as a sort of primitive version of ranges. In BAGL,
   those have all been changed to use C++ standard ranges and views, and can therefore be used with
   range-based for-loops and transformed with other views.
    - BAGL mostly tries to stick with standard ranges and views, but also has to introduce a number
      of its own views as the standard ranges library is still lacking many useful features (as of
      C++20, and even C++23, in my opinion).
    - Implementations of standard ranges in both `libstdc++` and `libc++` appear to still be
      rather experimental (as of 2024) in my estimation (even if, officially, they are not),
      so one should expect some rough edges.
 - The BGL uses many Boost library components that are now standard, e.g., `<memory>`,
   `<tuple>`, `<random>`, `<regex>`, etc., and obviously, BAGL swaps all of those out for their
   direct standard library equivalents with little to no change in behavior.
 - Anyone familiar with BGL will know that it relies heavily on its sibling the Boost.PropertyMap
   library. BAGL is, effectively, also a port of that (much smaller) library as well, which means
   that utilities usually found in `boost/property_map` can be found in `bagl` directly. And that
   is also true of a few other small bits and pieces of Boost that were brought into BAGL.
    - Property maps have been pretty significantly revamped to benefit from perfect forwarding
      and return type deduction. They function the same but differ a lot under-the-hood.
 - The BGL used Boost.Parameter library in some places to help reduce overload sets through
   named parameters. This has been removed in BAGL to avoid bringing in that library (which is
   also very template-heavy), but somewhat compensates for it with a richer concepts-based
   overload set for many of the algorithms.
 - The main "Swiss Army knife" graph data structure of the BGL, the `adjacency_list`, has been
   completely reimplemented. Historically, this class had to work around limitations of standard
   library containers (e.g., support for incomplete types) which created some very significant
   performance problems. Years ago, I solved these issues using Boost.Container. Now, the fixing
   of the standard library container defects (and the addition of `<unordered_set>`) means that
   a new and improved `adjacency_list` class was in order. The main differences are:
    - Removal of many unnecessary levels of indirection.
    - Support for unordered set for edge lists (for vertex lists, it doesn't make sense).
    - Support for `pool_s` containers which are "vectors with holes", that is, a free-list
      allocator for (reusable) removed elements that provides a useful "third way" compared to a
      contiguous random-access container (`vec_s`) or a persistent node-based container (`list_s`).
      Many graph mutation operations are very expensive (in the big-O sense) with contiguous
      containers (e.g., reallocating, copying, re-indexing), but on the other hand, traversals
      with node-based containers are abysmal, while random-access is rarely needed, and so, that
      is the niche that `pool_s` attempts to fill.
    - Supports perfect forwarding and move semantics.
 - GraphML support (reading and writing graphml files) has been migrated
   to [tinyxml2](https://github.com/leethomason/tinyxml2) from the original Boost-based xml parser
   and ad hoc writer.
 - The BGL test suite has been migrated to use the [googletest](https://github.com/google/googletest)
   library. The BGL test suite was very inconsistent, some used Boost.Test, some not, and some were
   just executables to be run manually. BAGL migrated all of those into a consistent test suite.

