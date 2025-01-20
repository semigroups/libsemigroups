//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2019-2023 Finn Smith
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <cstddef>        // for size_t
#include <stdexcept>      // for runtime_error
#include <unordered_set>  // for unordered_set
#include <vector>         // for vector

#include "catch_amalgamated.hpp"  // for REQUIRE, REQUIRE_THROWS_AS, REQUI...
#include "test-main.hpp"          // for LIBSEMIGROUPS_TEST_CASE
#include "word-graph-test-common.hpp"  // for add_clique etc

#include "libsemigroups/forest.hpp"      // for Forest
#include "libsemigroups/paths.hpp"       // for cbegin_pilo
#include "libsemigroups/word-graph.hpp"  // for WordGraph
#include "libsemigroups/word-range.hpp"  // for literals, WordRange

#include "libsemigroups/detail/string.hpp"  // for detail::to_string

namespace libsemigroups {

  using namespace literals;

  struct LibsemigroupsException;  // forward decl

  namespace {
    void add_path(WordGraph<size_t>& digraph, size_t n) {
      size_t old_nodes = digraph.number_of_nodes();
      digraph.add_nodes(n);
      for (size_t i = old_nodes; i < digraph.number_of_nodes() - 1; ++i) {
        digraph.target(i, 0, i + 1);
      }
    }

    WordGraph<size_t> path(size_t n) {
      WordGraph<size_t> g(0, 1);
      add_path(g, n);
      return g;
    }
  }  // namespace

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "000",
                          "constructor with 1  default arg",
                          "[quick][digraph]") {
    WordGraph<size_t> g;
    REQUIRE(g.number_of_nodes() == 0);
    REQUIRE(g.number_of_edges() == 0);
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "001",
                          "constructor with 0 default args",
                          "[quick][digraph]") {
    for (size_t j = 0; j < 100; ++j) {
      WordGraph<size_t> g(j);
      REQUIRE(g.number_of_nodes() == j);
      REQUIRE(g.number_of_edges() == 0);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph", "002", "add nodes", "[quick][digraph]") {
    WordGraph<size_t> g(3);
    REQUIRE(g.number_of_nodes() == 3);
    REQUIRE(g.number_of_edges() == 0);

    for (size_t i = 1; i < 100; ++i) {
      g.add_nodes(i);
      REQUIRE(g.number_of_nodes() == 3 + i * (i + 1) / 2);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph", "003", "add edges", "[quick][digraph]") {
    WordGraph<size_t> g(17, 31);

    for (size_t i = 0; i < 17; ++i) {
      for (size_t j = 0; j < 31; ++j) {
        g.target(i, j, (7 * i + 23 * j) % 17);
      }
    }

    REQUIRE(g.number_of_edges() == 31 * 17);
    REQUIRE(g.number_of_nodes() == 17);
    REQUIRE_THROWS_AS(g.target(0, 0, 32), LibsemigroupsException);

    for (size_t i = 0; i < 17; ++i) {
      for (size_t j = 0; j < 31; ++j) {
        REQUIRE(g.target(i, j) == (7 * i + 23 * j) % 17);
      }
    }

    g.add_to_out_degree(10);
    REQUIRE(g.out_degree() == 41);
    REQUIRE(g.number_of_nodes() == 17);
    REQUIRE(!word_graph::is_complete(g));

    for (size_t i = 0; i < 17; ++i) {
      for (size_t j = 0; j < 10; ++j) {
        g.target(i, 31 + j, (7 * i + 23 * j) % 17);
      }
    }

    REQUIRE(g.number_of_edges() == 41 * 17);
    REQUIRE(g.number_of_nodes() == 17);
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "004",
                          "exceptions",
                          "[quick][digraph]") {
    WordGraph<size_t> graph(10, 5);
    REQUIRE_THROWS_AS(graph.target(10, 0), LibsemigroupsException);
    REQUIRE(graph.target(0, 1) == UNDEFINED);

    REQUIRE_THROWS_AS(graph.target(0, 10, 0), LibsemigroupsException);
    REQUIRE_THROWS_AS(graph.target(10, 0, 0), LibsemigroupsException);
    for (size_t i = 0; i < 5; ++i) {
      graph.target(0, i, 1);
      graph.target(2, i, 2);
    }
    REQUIRE_NOTHROW(graph.target(0, 0, 1));
    REQUIRE_NOTHROW(graph.target(2, 0, 2));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph", "005", "random", "[quick][digraph]") {
    WordGraph<size_t> graph = WordGraph<size_t>::random(10, 10);
    REQUIRE(graph.number_of_nodes() == 10);
    REQUIRE(graph.number_of_edges() == 100);
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph", "006", "reserve", "[quick][digraph]") {
    WordGraph<size_t> graph;
    graph.reserve(10, 10);
    REQUIRE(graph.number_of_nodes() == 0);
    REQUIRE(graph.number_of_edges() == 0);
    graph.add_nodes(1);
    REQUIRE(graph.number_of_nodes() == 1);
    graph.add_nodes(9);
    REQUIRE(graph.number_of_nodes() == 10);
    REQUIRE(graph.number_of_edges() == 0);
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "007",
                          "default constructors",
                          "[quick][digraph]") {
    auto g1 = WordGraph<size_t>();
    g1.add_to_out_degree(1);
    word_graph::add_cycle(g1, 10);

    // Copy constructor
    auto g2(g1);
    REQUIRE(g2.number_of_edges() == 10);
    REQUIRE(g2.number_of_nodes() == 10);

    // Move constructor
    auto g3(std::move(g2));
    REQUIRE(g3.number_of_edges() == 10);
    REQUIRE(g3.number_of_nodes() == 10);

    // Copy assignment
    g2 = g3;
    REQUIRE(g2.number_of_edges() == 10);
    REQUIRE(g2.number_of_nodes() == 10);
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "008",
                          "iterator to edges",
                          "[quick][digraph]") {
    for (size_t n = 10; n < 512; n *= 4) {
      auto g = clique(n);
      REQUIRE(g.number_of_nodes() == n);
      REQUIRE(g.number_of_edges() == n * n);

      using node_type = decltype(g)::node_type;

      auto expected = std::vector<node_type>(n, 0);
      std::iota(expected.begin(), expected.end(), 0);

      for (auto it = g.cbegin_nodes(); it < g.cend_nodes(); ++it) {
        auto result = std::vector<node_type>(g.cbegin_targets(*it),
                                             g.cend_targets(*it));
        REQUIRE(result == expected);
      }
      REQUIRE_THROWS_AS(g.cbegin_targets(n), LibsemigroupsException);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "009",
                          "reverse node iterator",
                          "[quick]") {
    using node_type = WordGraph<size_t>::node_type;
    WordGraph<size_t> wg;
    wg.add_nodes(10);
    REQUIRE(wg.number_of_nodes() == 10);
    REQUIRE(std::vector<node_type>(wg.cbegin_nodes(), wg.cend_nodes())
            == std::vector<node_type>({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}));

    auto it = wg.cbegin_nodes();
    REQUIRE(*it == 0);
    auto copy(it);
    REQUIRE(*copy == 0);
    it       = wg.cend_nodes();
    auto tmp = it;
    REQUIRE(*--tmp == 9);
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "010",
                          "random/random_acyclic exceptions",
                          "[quick]") {
    // Too few nodes
    REQUIRE_THROWS_AS(word_graph::random_acyclic<size_t>(0, 0),
                      LibsemigroupsException);
    // Out degree too low
    REQUIRE_NOTHROW(WordGraph<size_t>::random(2, 0));

    REQUIRE_THROWS_AS(word_graph::random_acyclic<size_t>(2, 0),
                      LibsemigroupsException);
    REQUIRE_NOTHROW(word_graph::random_acyclic<size_t>(2, 2));
    auto wg = word_graph::random_acyclic<size_t>(2, 2);
    REQUIRE(word_graph::is_acyclic(wg));
    REQUIRE(word_graph::is_connected(wg));
    wg = word_graph::random_acyclic<size_t>(10, 10);
    REQUIRE(word_graph::is_acyclic(wg));
    REQUIRE(word_graph::is_connected(wg));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "011",
                          "unsafe (next) neighbour",
                          "[quick]") {
    auto wg = binary_tree(10);
    REQUIRE(wg.target_no_checks(0, 1) == wg.target(0, 1));
    REQUIRE(wg.next_label_and_target_no_checks(0, 1)
            == wg.next_label_and_target(0, 1));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "012",
                          "number_of_egdes incident to a node",
                          "[quick]") {
    auto wg = binary_tree(10);
    REQUIRE(wg.number_of_nodes() == 1023);
    REQUIRE(
        std::count_if(wg.cbegin_nodes(),
                      wg.cend_nodes(),
                      [&wg](size_t n) { return wg.number_of_edges(n) == 2; })
        == 511);
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "013",
                          "induced_subgraph_no_checks",
                          "[quick]") {
    WordGraph<size_t> wg;
    wg.add_nodes(3);
    wg.add_to_out_degree(2);
    wg.target(0, 0, 1);
    wg.target(1, 0, 0);
    wg.target(2, 0, 0);

    wg.induced_subgraph_no_checks(0, 2);
    static_assert(!IsWordGraph<size_t>);
    static_assert(IsWordGraph<WordGraph<size_t>>);
    REQUIRE(wg == make<WordGraph<size_t>>(2, {{1, UNDEFINED}, {0}}));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "014",
                          "remove_target_no_checks",
                          "[quick]") {
    WordGraph<size_t> wg;
    wg.add_nodes(3);
    wg.add_to_out_degree(2);
    wg.target(0, 0, 1);
    wg.target(1, 0, 0);
    wg.target(2, 0, 0);

    wg.remove_target_no_checks(0, 0);  // remove edge from 0 labelled 0
    REQUIRE(wg
            == make<WordGraph<size_t>>(3, {{UNDEFINED, UNDEFINED}, {0}, {0}}));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "015",
                          "swap_edge_no_checks",
                          "[quick]") {
    WordGraph<size_t> wg;
    wg.add_nodes(3);
    wg.add_to_out_degree(2);
    wg.target(0, 0, 1);
    wg.target(1, 0, 0);
    wg.target(2, 0, 2);

    // swap edge from 0 labelled 0 with edge from 1 labelled 0
    wg.swap_targets_no_checks(0, 1, 0);
    REQUIRE(wg == make<WordGraph<size_t>>(3, {{0, UNDEFINED}, {1}, {2}}));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph", "016", "operator<<", "[quick]") {
    WordGraph<uint64_t> wg;
    wg.add_nodes(3);
    wg.add_to_out_degree(2);
    wg.target(0, 0, 1);
    wg.target(1, 0, 0);
    wg.target(2, 0, 2);

    std::ostringstream oss;
    oss << wg;
    REQUIRE(oss.str()
            == "{3, {{1, 18446744073709551615}, {0, 18446744073709551615}, {2, "
               "18446744073709551615}}}");
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "017",
                          "is_acyclic | 2-cycle",
                          "[quick]") {
    WordGraph<size_t> wg;
    wg.add_nodes(2);
    wg.add_to_out_degree(1);
    wg.target(0, 0, 1);
    wg.target(1, 0, 0);
    REQUIRE(!word_graph::is_acyclic(wg));
    REQUIRE(word_graph::topological_sort(wg).empty());
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "018",
                          "is_acyclic | 1-cycle",
                          "[quick]") {
    WordGraph<size_t> wg;
    wg.add_nodes(1);
    wg.add_to_out_degree(1);
    wg.target(0, 0, 0);
    REQUIRE(!word_graph::is_acyclic(wg));
    REQUIRE(word_graph::topological_sort(wg).empty());
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "019",
                          "is_acyclic | multi-digraph",
                          "[quick]") {
    using node_type = WordGraph<size_t>::node_type;
    WordGraph<size_t> wg;
    wg.add_nodes(2);
    wg.add_to_out_degree(2);
    wg.target(0, 0, 1);
    wg.target(0, 1, 1);
    REQUIRE(word_graph::is_acyclic(wg));
    REQUIRE(word_graph::topological_sort(wg) == std::vector<node_type>({1, 0}));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "020",
                          "is_acyclic | complete digraph 100",
                          "[quick]") {
    WordGraph<size_t> wg;
    size_t const      n = 100;
    wg.add_nodes(n);
    wg.add_to_out_degree(n);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        if (i != j) {
          wg.target(i, j, j);
        }
      }
    }
    REQUIRE(!word_graph::is_acyclic(wg));
    REQUIRE(word_graph::topological_sort(wg).empty());
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "021",
                          "is_acyclic | acyclic digraph with 20000 nodes",
                          "[quick]") {
    WordGraph<size_t> wg;
    size_t const      n = 20000;
    wg.add_nodes(n);
    wg.add_to_out_degree(2);
    for (size_t i = 0; i < (n / 2 - 1); ++i) {
      wg.target(i, 0, i + 1);
    }
    wg.target(n / 2 - 1, 1, n - 1);
    wg.target(n / 2 + 1, 1, (3 * n) / 4 - 1);
    wg.target(n / 2, 1, 0);
    for (size_t i = n / 2; i < n - 1; ++i) {
      wg.target(i, 0, i + 1);
    }
    REQUIRE(word_graph::is_acyclic(wg));
    REQUIRE(word_graph::topological_sort(wg).size() == wg.number_of_nodes());
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "022",
                          "is_acyclic | acyclic digraph with 10 million nodes",
                          "[standard]") {
    WordGraph<size_t> wg;
    size_t const      n = 10000000;
    wg.add_nodes(n);
    wg.add_to_out_degree(2);
    for (size_t i = 0; i < (n / 2 - 1); ++i) {
      wg.target(i, 0, i + 1);
    }
    wg.target(n / 2 - 1, 1, n - 1);
    wg.target(n / 2 + 1, 1, (3 * n) / 4 - 1);
    wg.target(n / 2, 1, 0);
    for (size_t i = n / 2; i < n - 1; ++i) {
      wg.target(i, 0, i + 1);
    }
    REQUIRE(word_graph::is_acyclic(wg));
    REQUIRE(word_graph::topological_sort(wg).size() == n);
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "023",
                          "is_acyclic | for a node",
                          "[quick]") {
    using node_type = WordGraph<size_t>::node_type;
    WordGraph<size_t> wg;
    size_t const      n = 100;
    wg.add_nodes(n);
    wg.add_to_out_degree(2);
    for (size_t i = 0; i < n - 1; ++i) {
      wg.target(i, i % 2, i + 1);
    }
    word_graph::add_cycle(wg, 100);

    REQUIRE(std::all_of(
        wg.cbegin_nodes(), wg.cbegin_nodes() + 100, [&wg](node_type const& v) {
          return word_graph::is_acyclic(wg, v);
        }));

    REQUIRE(!word_graph::is_acyclic(wg));
    REQUIRE(word_graph::is_acyclic(wg, 10, 20));
    REQUIRE(!word_graph::is_acyclic(wg, 100, 120));
    REQUIRE(word_graph::is_acyclic(wg, 10, 120));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "024",
                          "is_acyclic | for a node | 2",
                          "[quick]") {
    WordGraph<size_t> wg;
    using node_type = decltype(wg)::node_type;
    wg.add_nodes(4);
    wg.add_to_out_degree(1);
    wg.target(0, 0, 1);
    wg.target(1, 0, 0);
    wg.target(2, 0, 3);
    REQUIRE(!word_graph::is_acyclic(wg));
    REQUIRE(word_graph::topological_sort(wg).empty());
    REQUIRE(!word_graph::is_acyclic(wg, 0));
    REQUIRE(word_graph::topological_sort(wg, 0).empty());
    REQUIRE(!word_graph::is_acyclic(wg, 1));
    REQUIRE(word_graph::topological_sort(wg, 1).empty());

    REQUIRE(word_graph::is_acyclic(wg, 2));
    REQUIRE(word_graph::topological_sort(wg, 2)
            == std::vector<node_type>({3, 2}));
    REQUIRE(word_graph::is_acyclic(wg, 3));
    REQUIRE(word_graph::topological_sort(wg, 3) == std::vector<node_type>({3}));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "025",
                          "is_reachable | acyclic 20 node digraph",
                          "[quick]") {
    WordGraph<size_t> wg;
    size_t const      n = 20;
    wg.add_nodes(n);
    wg.add_to_out_degree(2);
    for (size_t i = 0; i < (n / 2 - 1); ++i) {
      wg.target(i, 0, i + 1);
    }
    wg.target(n / 2 - 1, 1, n - 1);
    wg.target(n / 2 + 1, 1, (3 * n) / 4 - 1);
    wg.target(n / 2, 1, 0);
    for (size_t i = n / 2; i < n - 1; ++i) {
      wg.target(i, 0, i + 1);
    }
    REQUIRE(!word_graph::is_reachable(wg, 1, 10));
    REQUIRE(word_graph::is_reachable(wg, 10, 1));
    REQUIRE_THROWS_AS(word_graph::is_reachable(wg, 20, 1),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(word_graph::is_reachable(wg, 1, 20),
                      LibsemigroupsException);
    REQUIRE(word_graph::is_reachable(wg, 1, 1));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "026",
                          "is_reachable | 100 node path",
                          "[quick][no-valgrind]") {
    WordGraph<size_t> wg;
    size_t const      n = 100;
    wg.add_nodes(n);
    wg.add_to_out_degree(2);
    for (size_t i = 0; i < n - 1; ++i) {
      wg.target(i, i % 2, i + 1);
    }
    for (auto it1 = wg.cbegin_nodes(); it1 < wg.cend_nodes(); ++it1) {
      for (auto it2 = it1 + 1; it2 < wg.cend_nodes(); ++it2) {
        REQUIRE(word_graph::is_reachable(wg, *it1, *it2));
        REQUIRE(!word_graph::is_reachable(wg, *it2, *it1));
      }
    }
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "027",
                          "is_reachable | 100 node cycle",
                          "[quick][no-valgrind]") {
    WordGraph<size_t> wg;
    wg.add_to_out_degree(1);
    word_graph::add_cycle(wg, 100);
    for (auto it1 = wg.cbegin_nodes(); it1 < wg.cend_nodes(); ++it1) {
      for (auto it2 = it1 + 1; it2 < wg.cend_nodes(); ++it2) {
        REQUIRE(word_graph::is_reachable(wg, *it1, *it2));
        REQUIRE(word_graph::is_reachable(wg, *it2, *it1));
      }
    }
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "028",
                          "is_reachable | 20 node clique",
                          "[quick]") {
    WordGraph<size_t> wg = clique(20);
    for (auto it1 = wg.cbegin_nodes(); it1 < wg.cend_nodes(); ++it1) {
      for (auto it2 = it1 + 1; it2 < wg.cend_nodes(); ++it2) {
        REQUIRE(word_graph::is_reachable(wg, *it1, *it2));
        REQUIRE(word_graph::is_reachable(wg, *it2, *it1));
      }
    }
    REQUIRE(word_graph::is_complete(wg));
    REQUIRE(word_graph::topological_sort(wg).empty());
    REQUIRE(word_graph::topological_sort(wg, 0).empty());
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "029",
                          "follow_path | 20 node path",
                          "[quick]") {
    WordGraph<size_t> wg = path(20);
    for (auto it = cbegin_pilo(wg, 0); it != cend_pilo(wg); ++it) {
      REQUIRE(word_graph::follow_path(wg, 0, *it) == it.target());
      REQUIRE(word_graph::follow_path_no_checks(wg, 0, *it) == it.target());
    }
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "030",
                          "throw_if_label_out_of_bounds | 20 node path",
                          "[quick]") {
    WordGraph<size_t> wg = path(20);
    REQUIRE_THROWS_AS(word_graph::throw_if_label_out_of_bounds(wg, 10),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "031",
                          "last_node_on_path_no_checks | 20 node path",
                          "[quick]") {
    WordGraph<size_t> wg   = path(20);
    word_type         path = {};
    for (size_t i = 0; i < 19; ++i) {
      path.push_back(0);
      REQUIRE(word_graph::last_node_on_path_no_checks(
                  wg, 0, path.cbegin(), path.cend())
                  .first
              == i + 1);
    }
    path.push_back(0);
    auto p = word_graph::last_node_on_path_no_checks(
        wg, 0, path.cbegin(), path.cend());
    REQUIRE(p.first == 19);
    REQUIRE(p.second == path.cend() - 1);
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph", "032", "to_string", "[quick]") {
    WordGraph<uint64_t> wg = path(6);
    REQUIRE(detail::to_string(wg)
            == "{6, {{1}, {2}, {3}, {4}, {5}, {18446744073709551615}}}");
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph", "033", "make<WordGraph>", "[quick]") {
    auto wg = make<WordGraph<uint8_t>>(5, {{0, 0}, {1, 1}, {2}, {3, 3}});
    REQUIRE(detail::to_string(wg)
            == "{5, {{0, 0}, {1, 1}, {2, 255}, {3, 3}, {255, 255}}}");
    REQUIRE_THROWS_AS(
        make<WordGraph<uint8_t>>(5, {{0, 0}, {1, 1, 1}, {2}, {3, 3}}),
        LibsemigroupsException);
    wg = WordGraph<uint8_t>(5, 2);
    REQUIRE(
        detail::to_string(wg)
        == "{5, {{255, 255}, {255, 255}, {255, 255}, {255, 255}, {255, 255}}}");
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph", "034", "is_connected", "[quick]") {
    auto wg = make<WordGraph<size_t>>(5, {{0, 0}, {1, 1}, {2}, {3, 3}});
    REQUIRE(!word_graph::is_connected(wg));
    wg = path(1'000);
    REQUIRE(word_graph::is_connected(wg));
    REQUIRE(wg.number_of_nodes() == 1'000);
    word_graph::add_cycle(wg, 100);
    REQUIRE(wg.number_of_nodes() == 1'100);

    REQUIRE(!word_graph::is_connected(wg));
    wg.add_to_out_degree(1);
    wg.target(0, 1, 1'000);
    REQUIRE(word_graph::is_connected(wg));
    wg = WordGraph<size_t>();
    REQUIRE(word_graph::is_connected(wg));
  }
  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "035",
                          "is_strictly_cyclic",
                          "[quick][no-valgrind]") {
    auto wg = make<WordGraph<size_t>>(5, {{0, 0}, {1, 1}, {2}, {3, 3}});
    REQUIRE(!word_graph::is_strictly_cyclic(wg));
    wg = path(1'000);
    REQUIRE(word_graph::is_strictly_cyclic(wg));
    REQUIRE(wg.number_of_nodes() == 1'000);
    word_graph::add_cycle(wg, 100);
    REQUIRE(wg.number_of_nodes() == 1'100);

    REQUIRE(!word_graph::is_strictly_cyclic(wg));
    wg.add_to_out_degree(1);
    wg.target(0, 1, 1'000);
    REQUIRE(word_graph::is_strictly_cyclic(wg));
    wg = WordGraph<size_t>();
    REQUIRE(word_graph::is_strictly_cyclic(wg));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph", "046", "Joiner x 1", "[quick]") {
    WordGraph<size_t> x(
        make<WordGraph<size_t>>(3, {{0, 1, 2}, {0, 1, 2}, {0, 1, 2}}));
    WordGraph<size_t> y = x;

    Joiner join;

    WordGraph<size_t> xy;
    join(xy, x, y);
    REQUIRE(xy == x);
    join(xy, y, x);
    REQUIRE(xy == x);
    REQUIRE(join(x, y) == x);
    REQUIRE(join(y, x) == x);

    Meeter meet;
    meet(xy, x, y);
    REQUIRE(xy == x);
    REQUIRE(xy == y);

    y.target_no_checks(0, 0, 10);
    REQUIRE_THROWS_AS(join(x, y), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph", "050", "Joiner x 2", "[quick]") {
    WordGraph<size_t> x(
        make<WordGraph<size_t>>(3, {{1, 1, 1}, {2, 2, 2}, {2, 2, 2}}));

    WordGraph<size_t> y(
        make<WordGraph<size_t>>(3, {{1, 1, 2}, {1, 1, 2}, {1, 1, 2}}));

    WordGraph<size_t> xy;

    Joiner join;
    xy = join(x, y);
    REQUIRE(x != y);
    REQUIRE(xy == make<WordGraph<size_t>>(2, {{1, 1, 1}, {1, 1, 1}}));
    REQUIRE(join.is_subrelation(x, xy));
    REQUIRE(join.is_subrelation(y, xy));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph", "036", "Meeter x 1", "[quick]") {
    // These word graphs were taken from the lattice of
    // 2-sided congruences of the free semigroup with 2
    // generators.
    WordGraph<size_t> x(make<WordGraph<size_t>>(3, {{1, 2}, {1, 1}, {2, 2}}));
    WordGraph<size_t> y(make<WordGraph<size_t>>(3, {{1, 2}, {1, 1}, {1, 1}}));

    WordGraph<size_t> xy;

    Meeter meet;
    meet(xy, x, y);

    REQUIRE(xy == make<WordGraph<size_t>>(4, {{1, 2}, {1, 1}, {3, 3}, {3, 3}}));

    y = make<WordGraph<size_t>>(3, {{1, 2}, {2, 2}, {2, 2}});

    meet(xy, x, y);
    REQUIRE(xy == make<WordGraph<size_t>>(4, {{1, 2}, {3, 3}, {2, 2}, {3, 3}}));

    word_graph::standardize(xy);
    REQUIRE(xy == make<WordGraph<size_t>>(4, {{1, 2}, {3, 3}, {2, 2}, {3, 3}}));

    x = xy;
    meet(xy, x, y);
    REQUIRE(xy == make<WordGraph<size_t>>(4, {{1, 2}, {3, 3}, {2, 2}, {3, 3}}));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph", "037", "Meeter x 2", "[quick]") {
    auto x = make<WordGraph<size_t>>(5, {{1, 0}, {1, 2}, {1, 2}});
    auto y = make<WordGraph<size_t>>(5, {{0, 1}, {0, 1}});
    REQUIRE(word_graph::number_of_nodes_reachable_from(x, 0) == 3);
    REQUIRE(word_graph::number_of_nodes_reachable_from(y, 0) == 2);

    Meeter meet;
    auto   xy = meet(x, y);
    REQUIRE(xy == make<WordGraph<size_t>>(4, {{1, 2}, {1, 3}, {1, 2}, {1, 3}}));
    word_graph::standardize(xy);
    REQUIRE(xy == make<WordGraph<size_t>>(4, {{1, 2}, {1, 3}, {1, 2}, {1, 3}}));

    Joiner join;
    join(xy, x, y);
    REQUIRE(xy == make<WordGraph<size_t>>(1, {{0, 0}}));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph", "038", "Joiner incomplete", "[quick]") {
    WordGraph<uint32_t> wg(0, 1);
    word_graph::add_cycle(wg, 5);
    wg.remove_target(0, 0);
    Joiner join;
    REQUIRE(join(wg, wg) == make<WordGraph<uint32_t>>(1, {{UNDEFINED}}));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph", "039", "Meeter incomplete", "[quick]") {
    WordGraph<uint32_t> wg(0, 1);
    word_graph::add_cycle(wg, 5);
    wg.remove_target(0, 0);
    Meeter meet;
    REQUIRE(meet(wg, wg) == make<WordGraph<uint32_t>>(1, {{UNDEFINED}}));
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "040",
                          "WordGraph to_input_string",
                          "[quick]") {
    WordGraph<uint32_t> wg(0, 1);
    word_graph::add_cycle(wg, 5);
    wg.remove_target(0, 0);
    REQUIRE(to_input_string(wg) == "5, {{4294967295}, {2}, {3}, {4}, {0}}");
    REQUIRE(
        to_input_string(wg, "make<WordGraph<uint32_t>>(", "[]", ")")
        == "make<WordGraph<uint32_t>>(5, [[4294967295], [2], [3], [4], [0]])");
  }

  LIBSEMIGROUPS_TEST_CASE("WordGraph",
                          "041",
                          "WordGraph to_input_string",
                          "[quick]") {
    WordGraph<uint32_t> wg(0, 1);
    word_graph::add_cycle(wg, 5);
    wg.remove_target(0, 0);
    REQUIRE_NOTHROW(wg.hash_value());
  }
}  // namespace libsemigroups
