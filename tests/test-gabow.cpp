//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2023 James Mitchell
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

#include <algorithm>  // for all_of
#include <cstddef>    // for size_t
#include <numeric>    // for iota
#include <vector>     // for vector, operator==...

#include "catch_amalgamated.hpp"       // for operator""_catch_sr, Ass...
#include "test-main.hpp"               // for LIBSEMIGROUPS_TEST_CASE
#include "word-graph-test-common.hpp"  // for clique, add_clique

#include "libsemigroups/constants.hpp"   // for UNDEFINED, Undefined, Max
#include "libsemigroups/exception.hpp"   // for LibsemigroupsException
#include "libsemigroups/forest.hpp"      // for make<Forest>, Forest, opera...
#include "libsemigroups/gabow.hpp"       // for Gabow
#include "libsemigroups/ranges.hpp"      // for equal
#include "libsemigroups/word-graph.hpp"  // for WordGraph, to_action...

#include "libsemigroups/ranges.hpp"  // for Inner, operator|, remove...

namespace libsemigroups {

  using namespace rx;

  LIBSEMIGROUPS_TEST_CASE("Gabow", "000", "for a cycle", "[quick][gabow]") {
    WordGraph<size_t> wg;
    wg.add_to_out_degree(1);
    word_graph::add_cycle(wg, 33);
    // REQUIRE(wg.scc_id(0) == 0);
    Gabow scc(wg);
    REQUIRE(scc.id(0) == 0);
    REQUIRE(scc.components()
            == std::vector<std::vector<size_t>>(
                {{32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22,
                  21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11,
                  10, 9,  8,  7,  6,  5,  4,  3,  2,  1,  0}}));
    REQUIRE((wg.nodes() | all_of([&scc](auto v) { return scc.id(v) == 0; })));
  }

  LIBSEMIGROUPS_TEST_CASE("Gabow",
                          "001",
                          "no edges",
                          "[quick][gabow][no-valgrind]") {
    WordGraph<size_t> wg = WordGraph<size_t>(0);
    Gabow             scc(wg);
    for (size_t j = 1; j < 100; ++j) {
      wg.add_nodes(j);
      scc.init(wg);

      for (size_t i = 0; i < j * (j + 1) / 2; ++i) {
        REQUIRE(scc.id(i) == i);
      }
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Gabow", "002", "disjoint cycles", "[quick][gabow]") {
    WordGraph<size_t> wg;
    wg.add_to_out_degree(1);
    Gabow scc(wg);
    for (size_t j = 2; j < 50; ++j) {
      word_graph::add_cycle(wg, j);
      scc.init(wg);
      REQUIRE((wg.nodes()
               | filter([&scc, j](auto v) { return scc.id(v) == j - 2; })
               | count())
              == j);
    }

    REQUIRE(wg.number_of_nodes() == 1'224);
    REQUIRE(wg.number_of_edges() == 1'224);
    REQUIRE(word_graph::is_complete(wg));
  }

  LIBSEMIGROUPS_TEST_CASE("Gabow", "003", "complete graphs", "[quick][gabow]") {
    for (size_t k = 2; k < 50; ++k) {
      auto  wg = clique(k);
      Gabow scc(wg);
      for (size_t i = 0; i < k; ++i) {
        REQUIRE(scc.id(i) == 0);
      }
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Gabow", "004", "exceptions", "[quick][gabow]") {
    WordGraph<size_t> wg(10, 5);
    Gabow             scc(wg);
    REQUIRE_THROWS_AS(scc.id(10), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("Gabow",
                          "005",
                          "spanning forest - complete graphs",
                          "[quick][gabow]") {
    for (size_t k = 2; k < 50; ++k) {
      auto  wg = clique(k);
      Gabow scc(wg);
      REQUIRE(scc.number_of_components() == 1);

      Forest const& f = scc.spanning_forest();
      REQUIRE(f.parent(k - 1) == UNDEFINED);
      REQUIRE_NOTHROW(scc.reverse_spanning_forest());
    }
    auto  wg = clique(3);
    Gabow scc(wg);
    REQUIRE(scc.number_of_components() == 1);
    {
      auto const& f = scc.spanning_forest();
      REQUIRE(f == make<Forest>({2, 2, UNDEFINED}, {0, 1, UNDEFINED}));
    }
    {
      auto const& f = scc.reverse_spanning_forest();
      REQUIRE(f == make<Forest>({2, 2, UNDEFINED}, {2, 2, UNDEFINED}));
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Gabow",
                          "006",
                          "spanning forest - disjoint cycles",
                          "[quick][gabow]") {
    size_t            j = 33;
    WordGraph<size_t> wg;
    wg.add_to_out_degree(1);

    for (size_t k = 0; k < 10; ++k) {
      wg.add_nodes(j);
      for (size_t i = k * j; i < (k + 1) * j - 1; ++i) {
        wg.target(i, 0, i + 1);
      }
      wg.target((k + 1) * j - 1, 0, k * j);
    }

    Gabow scc(wg);
    for (size_t i = 0; i < 10 * j; ++i) {
      REQUIRE(scc.id(i) == i / j);
    }
    Forest const& forest = scc.spanning_forest();

    REQUIRE(forest.parents()
            == std::vector<size_t>(
                {32,  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,
                 10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  20,
                 21,  22,  23,  24,  25,  26,  27,  28,  29,  30,  UNDEFINED,
                 65,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,
                 43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,
                 54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  UNDEFINED,
                 98,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,
                 76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,
                 87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  UNDEFINED,
                 131, 99,  100, 101, 102, 103, 104, 105, 106, 107, 108,
                 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
                 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, UNDEFINED,
                 164, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141,
                 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152,
                 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, UNDEFINED,
                 197, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174,
                 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185,
                 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, UNDEFINED,
                 230, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
                 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218,
                 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, UNDEFINED,
                 263, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240,
                 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251,
                 252, 253, 254, 255, 256, 257, 258, 259, 260, 261, UNDEFINED,
                 296, 264, 265, 266, 267, 268, 269, 270, 271, 272, 273,
                 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284,
                 285, 286, 287, 288, 289, 290, 291, 292, 293, 294, UNDEFINED,
                 329, 297, 298, 299, 300, 301, 302, 303, 304, 305, 306,
                 307, 308, 309, 310, 311, 312, 313, 314, 315, 316, 317,
                 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, UNDEFINED}));
  }

  // TODO(FLS) uncomment this
  //  LIBSEMIGROUPS_TEST_CASE("Gabow",
  //                          "007",
  //                          "scc root paths - complete graphs",
  //                          "[quick][gabow]") {
  //    for (size_t k = 2; k < 50; ++k) {
  //      WordGraph<size_t> graph(k);
  //
  //      for (size_t i = 0; i < k; ++i) {
  //        for (size_t j = 0; j < k; ++j) {
  //          graph.target(i,  j,  j);
  //        }
  //      }
  //
  //      for (size_t i = 0; i < k; ++i) {
  //        size_t pos = i;
  //        for (auto it = graph.cbegin_path_to_root(i); it <
  //        graph.cend_path_to_root(i); ++it) {
  //          pos = graph.target(pos, *it);
  //        }
  //        REQUIRE(pos == graph.cbegin_sccs()[graph.scc_id(i)][0]);
  //      }
  //    }
  //  }

  // LIBSEMIGROUPS_TEST_CASE("Gabow",
  //                         "008",
  //                         "scc root paths - disjoint cycles",
  //                         "[quick][gabow]") {
  //   for (size_t j = 2; j < 35; ++j) {
  //     WordGraph<size_t> graph;

  //     for (size_t k = 0; k < 6; ++k) {
  //       graph.add_nodes(j);
  //       for (size_t i = k * j; i < (k + 1) * j - 1; ++i) {
  //         graph.target(i,  0,  i + 1);
  //       }
  //       graph.target((k + 1) * j - 1,  0,  k * j);
  //     }

  //     for (size_t i = 0; i < graph.number_of_nodes(); ++i) {
  //       size_t pos = i;
  //       for (auto it = graph.cbegin_path_to_root(i);
  //            it < graph.cend_path_to_root(i);
  //            ++it) {
  //         pos = graph.target(pos, *it);
  //       }
  //       REQUIRE(pos == graph.cbegin_sccs()[graph.scc_id(i)][0]);
  //     }
  //   }
  // }

  LIBSEMIGROUPS_TEST_CASE("Gabow",
                          "009",
                          "large cycle",
                          "[quick][gabow][no-valgrind]") {
    WordGraph<size_t> wg;
    wg.add_to_out_degree(1);
    word_graph::add_cycle(wg, 100000);
    using node_type = decltype(wg)::node_type;
    Gabow scc(wg);
    REQUIRE(
        (wg.nodes() | all_of([&scc](node_type i) { return scc.id(i) == 0; })));

    word_graph::add_cycle(wg, 10101);
    scc.init(wg);
    REQUIRE((wg.nodes() | take(100000)
             | all_of([&scc](node_type i) { return scc.id(i) == 0; })));
    REQUIRE((wg.nodes() | skip_n(100000)
             | all_of([&scc](node_type i) { return scc.id(i) == 1; })));
  }

  LIBSEMIGROUPS_TEST_CASE("Gabow", "010", "component", "[quick][gabow]") {
    using node_type = decltype(clique(1))::node_type;

    for (size_t n = 10; n < 512; n *= 4) {
      auto wg = clique(n);
      REQUIRE(wg.number_of_nodes() == n);
      REQUIRE(wg.number_of_edges() == n * n);
      Gabow scc(wg);
      REQUIRE(scc.number_of_components() == 1);

      add_clique(wg, n);

      REQUIRE(wg.number_of_nodes() == 2 * n);
      REQUIRE(wg.number_of_edges() == 2 * n * n);

      scc.init(wg);
      REQUIRE(scc.number_of_components() == 2);

      auto expected = std::vector<node_type>(n, 0);
      std::iota(expected.begin(), expected.end(), 0);
      auto result = scc.component(0);
      std::sort(result.begin(), result.end());
      REQUIRE(result == expected);

      std::iota(expected.begin(), expected.end(), n);
      result = scc.component(1);
      std::sort(result.begin(), result.end());
      REQUIRE(result == expected);
      REQUIRE_THROWS_AS(scc.component(2), LibsemigroupsException);

      result = (scc.roots() | transform([&scc](auto v) { return scc.id(v); })
                | to_vector());
      REQUIRE(result == std::vector<node_type>({0, 1}));
    }
    {
      auto wg = clique(10);
      for (size_t n = 0; n < 99; ++n) {
        add_clique(wg, 10);
      }
      REQUIRE(wg.number_of_nodes() == 1000);
      REQUIRE(wg.number_of_edges() == 10000);

      Gabow scc(wg);
      REQUIRE(scc.number_of_components() == 100);

      auto result
          = (scc.roots() | transform([&scc](auto v) { return scc.id(v); }));

      auto expected = (seq<size_t>() | take(100));
      REQUIRE(equal(result, expected));
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Gabow", "011", "root of scc", "[quick][gabow]") {
    auto wg = clique(10);
    for (size_t n = 0; n < 99; ++n) {
      add_clique(wg, 10);
    }
    REQUIRE(wg.number_of_nodes() == 1000);
    REQUIRE(wg.number_of_edges() == 10000);

    Gabow scc(wg);
    REQUIRE(scc.number_of_components() == 100);

    using node_type = decltype(wg)::node_type;

    for (auto const& comp : scc.components()) {
      REQUIRE(
          std::all_of(comp.cbegin(), comp.cend(), [&scc, &comp](node_type v) {
            return scc.root_of(v) == comp[0];
          }));
    }
    REQUIRE_THROWS_AS(scc.root_of(1000), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("Gabow",
                          "012",
                          "reverse_spanning_forest",
                          "[quick][gabow]") {
    auto wg = make<WordGraph<size_t>>(
        5, {{0, 1, 4, 3}, {2}, {2, 0, 3, 3}, {4, 1}, {1, 0, 2}});
    Gabow scc(wg);
    REQUIRE(scc.number_of_components() == 1);
    REQUIRE(scc.reverse_spanning_forest()
            == make<Forest>({4, 2, 0, 4, UNDEFINED}, {2, 0, 1, 0, UNDEFINED}));
  }

  LIBSEMIGROUPS_TEST_CASE("Gabow",
                          "013",
                          "to_human_readable_repr",
                          "[quick][gabow]") {
    auto wg = make<WordGraph<size_t>>(
        5, {{0, 1, 4, 3}, {2}, {2, 0, 3, 3}, {4, 1}, {1, 0, 2}});
    Gabow scc(wg);
    REQUIRE(to_human_readable_repr(scc)
            == "<Gabow with 5 nodes and components not yet found>");
    REQUIRE(scc.number_of_components() == 1);
    REQUIRE(to_human_readable_repr(scc)
            == "<Gabow with 5 nodes and 1 component>");
  }
}  // namespace libsemigroups
