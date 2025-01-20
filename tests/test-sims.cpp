//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2022-24 James D. Mitchell
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

#define CATCH_CONFIG_ENABLE_PAIR_STRINGMAKER
#define CATCH_CONFIG_ENABLE_TUPLE_STRINGMAKER

#include <algorithm>    // for fill, transform, equal, sort, for_each
#include <array>        // for array
#include <atomic>       // for atomic, atomic_uint64_t, __atomic_b...
#include <cctype>       // for isprint
#include <chrono>       // for milliseconds
#include <cmath>        // for pow
#include <cstddef>      // for size_t
#include <cstdint>      // for uint64_t, uint32_t, uint8_t
#include <iostream>     // for cout
#include <iterator>     // for distance
#include <string>       // for basic_string, operator==, string
#include <thread>       // for thread
#include <tuple>        // for tuple, operator==
#include <type_traits>  // for decay_t
#include <utility>      // for pair, move, swap
#include <vector>       // for vector

#include "catch_amalgamated.hpp"  // for REQUIRE, REQUIRE_THROWS_AS, REQUI...
#include "test-main.hpp"          // for LIBSEMIGROUPS_TEST_CASE

#include "libsemigroups/bipart.hpp"           // for Bipartition
#include "libsemigroups/constants.hpp"        // for UNDEFINED
#include "libsemigroups/exception.hpp"        // for LibsemigroupsException
#include "libsemigroups/fpsemi-examples.hpp"  // for brauer_monoid etc
#include "libsemigroups/froidure-pin.hpp"     // for FroidurePin
#include "libsemigroups/gabow.hpp"            // for Gabow
#include "libsemigroups/knuth-bendix.hpp"     // for redundant_rule
#include "libsemigroups/matrix.hpp"           // for DynamicMatrix, BMat
#include "libsemigroups/order.hpp"            // for shortlex_compare, Order
#include "libsemigroups/presentation.hpp"     // for Presentation
#include "libsemigroups/ranges.hpp"           // for shortlex_compare
#include "libsemigroups/sims.hpp"             // for Sims1
#include "libsemigroups/to-froidure-pin.hpp"  // for make
#include "libsemigroups/to-presentation.hpp"  // for to_presentation
#include "libsemigroups/todd-coxeter.hpp"     // for ToddCoxeter
#include "libsemigroups/transf.hpp"           // for Transf
#include "libsemigroups/types.hpp"            // for word_type
#include "libsemigroups/word-graph.hpp"       // for WordGraph

#include "libsemigroups/detail/eigen.hpp"           // for DenseBase::row
#include "libsemigroups/detail/felsch-graph.hpp"    // for FelschGraph
#include "libsemigroups/detail/fmt.hpp"             // for format, print
#include "libsemigroups/detail/iterator.hpp"        // for operator+
#include "libsemigroups/detail/path-iterators.hpp"  // for const_pstilo_iterator
#include "libsemigroups/detail/report.hpp"          // for ReportGuard
#include "libsemigroups/detail/word-graph-with-sources.hpp"  // for WordGra...

#include "rx/ranges.hpp"  // for operator|, iterator_range

namespace libsemigroups {

  using word_graph_type = typename Sims1::word_graph_type;
  using node_type       = typename word_graph_type::node_type;

  using namespace literals;

  using fpsemigroup::author;
  using fpsemigroup::brauer_monoid;
  using fpsemigroup::chinese_monoid;
  using fpsemigroup::fibonacci_semigroup;
  using fpsemigroup::full_transformation_monoid;
  using fpsemigroup::monogenic_semigroup;
  using fpsemigroup::partition_monoid;
  using fpsemigroup::plactic_monoid;
  using fpsemigroup::rectangular_band;
  using fpsemigroup::stellar_monoid;
  using fpsemigroup::stylic_monoid;
  using fpsemigroup::symmetric_inverse_monoid;
  using fpsemigroup::temperley_lieb_monoid;
  using fpsemigroup::zero_rook_monoid;

  namespace {
    // TODO(2) check_exclude
    template <typename P>
    void check_include(P const& p, std::vector<word_type> const& e, size_t n) {
      auto foo = [&e](auto const& wg) {
        using word_graph::follow_path_no_checks;
        for (auto it = e.cbegin(); it != e.cend(); it += 2) {
          if (follow_path_no_checks(wg, 0, *it)
              != follow_path_no_checks(wg, 0, *(it + 1))) {
            return false;
          }
        }
        return true;
      };
      Sims1 S(p);
      Sims1 T(p);

      for (auto it = e.cbegin(); it != e.cend(); it += 2) {
        sims::add_included_pair(T, *it, *(it + 1));
      }

      REQUIRE(T.included_pairs() == e);

      // auto it = std::find_if(T.cbegin(n), T.cend(n), foo);
      // REQUIRE(*it == make<WordGraph<uint32_t>>(1, {{0, 0}}));
      REQUIRE(std::all_of(T.cbegin(n), T.cend(n), foo));
      REQUIRE(static_cast<uint64_t>(std::count_if(S.cbegin(n), S.cend(n), foo))
              == T.number_of_congruences(n));
    }

    constexpr unsigned int factorial(unsigned int n) {
      return n > 1 ? n * factorial(n - 1) : 1;
    }

    // The following checks whether the return value of generating_pairs is
    // correct by using a ToddCoxeter instance. This works for both left and
    // onesided congruences because we really compute a onesided congruence on
    // the dual semigroup when setting the kind to left. Thus when we get the
    // generating pairs they generate that onesided congruence on the dual,
    // which the function below checks. This seems potentially confusing.
    template <typename Sims1Or2, typename Node>
    void check_right_generating_pairs(Sims1Or2 const&        s,
                                      WordGraph<Node> const& wg) {
      ToddCoxeter tc(congruence_kind::onesided, s.presentation());

      for (auto const& p : sims::right_generating_pairs(wg)) {
        todd_coxeter::add_generating_pair(tc, p.first, p.second);
      }
      tc.run();

      tc.standardize(Order::shortlex);
      auto expected = tc.word_graph();
      auto result   = wg;

      REQUIRE(expected.number_of_nodes() >= result.number_of_active_nodes());

      expected.induced_subgraph_no_checks(0, result.number_of_active_nodes());
      result.induced_subgraph_no_checks(0, result.number_of_active_nodes());
      REQUIRE(result == expected);

      tc.init(congruence_kind::onesided, s.presentation());
      for (auto const& p : sims::right_generating_pairs(s.presentation(), wg)) {
        todd_coxeter::add_generating_pair(tc, p.first, p.second);
      }
      tc.run();
      tc.standardize(Order::shortlex);
      expected = tc.word_graph();
      REQUIRE(expected.number_of_nodes() >= result.number_of_active_nodes());

      expected.induced_subgraph_no_checks(0, result.number_of_active_nodes());
      REQUIRE(result == expected);
    }

    void check_two_sided_generating_pairs(Sims2 const&                  s,
                                          Sims2::word_graph_type const& wg) {
      ToddCoxeter tc(congruence_kind::twosided, s.presentation());

      for (auto const& p : sims::right_generating_pairs(wg)) {
        todd_coxeter::add_generating_pair(tc, p.first, p.second);
      }
      tc.run();

      tc.standardize(Order::shortlex);
      auto expected = tc.word_graph();
      auto result   = wg;

      REQUIRE(expected.number_of_nodes() >= result.number_of_active_nodes());

      expected.induced_subgraph_no_checks(0, result.number_of_active_nodes());
      result.induced_subgraph_no_checks(0, result.number_of_active_nodes());
      REQUIRE(result == expected);

      tc.init(congruence_kind::onesided, s.presentation());
      for (auto const& p : sims::right_generating_pairs(s.presentation(), wg)) {
        todd_coxeter::add_generating_pair(tc, p.first, p.second);
      }
      tc.run();
      tc.standardize(Order::shortlex);
      expected = tc.word_graph();
      REQUIRE(expected.number_of_nodes() >= result.number_of_active_nodes());

      expected.induced_subgraph_no_checks(0, result.number_of_active_nodes());
      REQUIRE(result == expected);
    }

    template <typename Iterator>
    void check_meets_and_joins(Iterator first, Iterator last) {
      using WordGraph_ = std::decay_t<decltype(*first)>;
      std::vector<WordGraph_> graphs(first, last);
      size_t const            n = graphs.size();
      Joiner                  joiner;
      Meeter                  meeter;
      WordGraph_              tmp;
      for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
          REQUIRE(std::tuple(meeter.is_subrelation(graphs[i], graphs[j]),
                             graphs[i],
                             graphs[j])
                  == std::tuple(joiner.is_subrelation(graphs[i], graphs[j]),
                                graphs[i],
                                graphs[j]));
          // TODO(1): FIXME the below doesn't seem to work, but JDM expected
          // it to.
          //
          // joiner(tmp, graphs[i], graphs[j]);
          // REQUIRE(meeter.is_subrelation(graphs[j], tmp));
          // REQUIRE(meeter.is_subrelation(graphs[i], tmp));
          // REQUIRE(joiner.is_subrelation(graphs[j], tmp));
          // REQUIRE(joiner.is_subrelation(graphs[i], tmp));
          // REQUIRE(meeter.is_subrelation(tmp, graphs[j]));
          // REQUIRE(meeter.is_subrelation(tmp, graphs[i]));
          // REQUIRE(joiner.is_subrelation(tmp, graphs[j]));
          // REQUIRE(joiner.is_subrelation(tmp, graphs[i]));
        }
      }
    }

    template <typename Sims1Or2>
    void check_congruence_count_with_free_object(Sims1Or2 const& sims,
                                                 size_t          index,
                                                 size_t          expected) {
      auto                    p = sims.presentation();
      Sims1Or2                SF(sims);
      Presentation<word_type> F;
      F.alphabet(p.alphabet());
      F.contains_empty_word(p.contains_empty_word());
      std::atomic<size_t> count = 0;
      SF.presentation(F);
      SF.for_each(index, [&p, &count](auto const& wg) {
        count += word_graph::is_compatible(wg,
                                           wg.cbegin_nodes(),
                                           wg.cbegin_nodes()
                                               + wg.number_of_active_nodes(),
                                           p.rules.cbegin(),
                                           p.rules.cend());
      });
      REQUIRE(count == expected);
    }

  }  // namespace

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "000",
                          "fp example 1",
                          "[quick][low-index]") {
    auto rg = ReportGuard(false);

    Presentation<word_type> p;
    p.contains_empty_word(true);
    p.alphabet(01_w);
    presentation::add_rule(p, 000_w, 0_w);
    presentation::add_rule(p, 11_w, 1_w);
    presentation::add_rule(p, 0101_w, 0_w);

    {
      Sims1 S;
      REQUIRE(S.presentation(p).number_of_threads(2).number_of_congruences(5)
              == 6);
      REQUIRE_THROWS_AS(S.number_of_congruences(0), LibsemigroupsException);
      REQUIRE_THROWS_AS(S.find_if(0, [](auto) { return false; }),
                        LibsemigroupsException);
      REQUIRE_THROWS_AS(S.for_each(0, [](auto) { return false; }),
                        LibsemigroupsException);
      REQUIRE_THROWS_AS(S.cbegin(0), LibsemigroupsException);
      REQUIRE_THROWS_AS(S.cend(0), LibsemigroupsException);
      REQUIRE(S.number_of_congruences(1) == 1);

      auto it = S.cbegin(1);
      REQUIRE(*it == make<WordGraph<node_type>>(1, {{0, 0}}));

      it = S.cbegin(5);
      REQUIRE(*(it++) == make<WordGraph<node_type>>(5, {{0, 0}}));
      REQUIRE(*(it++) == make<WordGraph<node_type>>(5, {{1, 0}, {1, 1}}));
      REQUIRE(*(it++) == make<WordGraph<node_type>>(5, {{1, 1}, {1, 1}}));
      REQUIRE(*(it++)
              == make<WordGraph<node_type>>(5, {{1, 2}, {1, 1}, {1, 2}}));
      REQUIRE(*(it++)
              == make<WordGraph<node_type>>(5, {{1, 2}, {1, 1}, {2, 2}}));
      REQUIRE(
          *(it++)
          == make<WordGraph<node_type>>(5, {{1, 2}, {1, 1}, {3, 2}, {3, 3}}));
      REQUIRE(*(it++) == WordGraph<node_type>(0, 2));
      REQUIRE(*(it++) == WordGraph<node_type>(0, 2));
      REQUIRE(*(it++) == WordGraph<node_type>(0, 2));

      it = S.cbegin(3);
      REQUIRE(*it == make<WordGraph<node_type>>(3, {{0, 0}}));
      // Note that Catch's REQUIRE macro is not thread safe, see:
      // https://github.com/catchorg/Catch2/issues/99
      // as such we cannot call any function (like
      // check_right_generating_pairs) that uses REQUIRE in multiple threads.
      S.number_of_threads(1).for_each(
          5, [S](auto const& wg) { check_right_generating_pairs(S, wg); });
    }
    // [[[0, 0]],
    // [[1, 2], [1, 1], [3, 2], [3, 3]],
    // [[1, 2], [1, 1], [2, 2]],
    // [[1, 2], [1, 1], [1, 2]],
    // [[1, 1], [1, 1]],
    // [[1, 0], [1, 1]]]
    {
      presentation::reverse(p);
      Sims1 S;
      REQUIRE(S.presentation(p).number_of_congruences(5) == 9);
      for (auto it = S.cbegin(5); it != S.cend(5); ++it) {
        REQUIRE(word_graph::follow_path_no_checks(*it, 0, 1010_w)
                == word_graph::follow_path_no_checks(*it, 0, {0}));
      }
      S.for_each(5,
                 [&S](auto const& wg) { check_right_generating_pairs(S, wg); });
      auto mat = sims::poset(S.cbegin(5), S.cend(5));
      REQUIRE(mat
              == BMat<>({{0, 0, 0, 0, 0, 0, 0, 0, 0},
                         {1, 0, 0, 0, 0, 0, 0, 0, 0},
                         {1, 0, 0, 0, 0, 0, 0, 0, 0},
                         {1, 0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 1, 1, 0, 0, 0, 0, 0},
                         {1, 0, 0, 0, 0, 0, 0, 0, 0},
                         {0, 0, 0, 1, 0, 1, 0, 0, 0},
                         {0, 1, 1, 0, 0, 1, 0, 0, 0},
                         {0, 0, 0, 0, 1, 0, 1, 1, 0}}));
      check_meets_and_joins(S.cbegin(5), S.cend(5));
      // sims::dot_poset("example-000", S.cbegin(5), S.cend(5));
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "001",
                          "fp example 2",
                          "[quick][low-index][no-valgrind]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.contains_empty_word(true);

    p.alphabet(012_w);
    presentation::add_rule(p, 010_w, 00_w);
    presentation::add_rule(p, 22_w, 00_w);
    presentation::add_rule(p, 000_w, 00_w);
    presentation::add_rule(p, 21_w, 12_w);
    presentation::add_rule(p, 20_w, 00_w);
    presentation::add_rule(p, 11_w, 1_w);
    presentation::add_rule(p, 02_w, 00_w);

    Sims1 S;
    S.presentation(p);
    REQUIRE(S.number_of_congruences(1) == 1);
    REQUIRE(S.number_of_congruences(2) == 3);
    REQUIRE(S.number_of_congruences(3) == 13);
    REQUIRE(S.number_of_congruences(4) == 36);
    REQUIRE(S.number_of_congruences(5) == 82);
    REQUIRE(S.number_of_congruences(6) == 135);
    REQUIRE(S.number_of_congruences(7) == 166);
    REQUIRE(S.number_of_congruences(8) == 175);
    REQUIRE(S.number_of_congruences(9) == 176);
    REQUIRE(S.number_of_congruences(10) == 176);

    auto it = S.cbegin(2);
    REQUIRE(*(it++) == make<WordGraph<node_type>>(2, {{0, 0, 0}}));
    REQUIRE(*(it++) == make<WordGraph<node_type>>(2, {{1, 0, 1}, {1, 1, 1}}));
    REQUIRE(*(it++) == make<WordGraph<node_type>>(2, {{1, 1, 1}, {1, 1, 1}}));
    REQUIRE(*(it++) == WordGraph<node_type>(0, 3));
    REQUIRE(*(it++) == WordGraph<node_type>(0, 3));
    // sims::dot_poset("example-001", S.cbegin(4), S.cend(4));

    presentation::reverse(p);
    S.init(p);
    REQUIRE(S.number_of_congruences(11) == 176);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "002",
                          "ToddCoxeter failing example (word_type)",
                          "[quick][low-index][no-valgrind]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.contains_empty_word(false);

    //         aAbBcCe
    p.alphabet(0123456_w);
    presentation::add_identity_rules(p, 6);
    presentation::add_inverse_rules(p, 1032546_w, 6);
    presentation::add_rule(p, 00504_w, 6_w);
    presentation::add_rule(p, 0422152_w, 6_w);
    presentation::add_rule(p, 1302444_w, 6_w);
    Sims1 S;
    S.presentation(p);
    // REQUIRE(S.number_of_congruences(1) == 1);
    // REQUIRE(S.number_of_congruences(3) == 14);
    // REQUIRE(S.number_of_congruences(4) == 14);
    // REQUIRE(S.number_of_congruences(5) == 14);
    // sims::dot_poset("example-002", S.cbegin(3), S.cend(3));

    S.for_each(3,
               [&S](auto const& wg) { check_right_generating_pairs(S, wg); });

    auto it = S.cbegin(3);
    REQUIRE((sims::right_generating_pairs(*it) | rx::to_vector())
            == std::vector<relation_type>({{1_w, 0_w},
                                           {2_w, 0_w},
                                           {3_w, 0_w},
                                           {4_w, 0_w},
                                           {5_w, 0_w},
                                           {6_w, 0_w},
                                           {00_w, 0_w},
                                           {01_w, 0_w},
                                           {02_w, 0_w},
                                           {03_w, 0_w},
                                           {04_w, 0_w},
                                           {05_w, 0_w},
                                           {06_w, 0_w}}));

    REQUIRE((sims::right_generating_pairs(p, *it) | rx::to_vector())
            == std::vector<relation_type>({{1_w, 0_w},
                                           {2_w, 0_w},
                                           {3_w, 0_w},
                                           {4_w, 0_w},
                                           {5_w, 0_w},
                                           {6_w, 0_w}}));

    check_right_generating_pairs(S, *it);

    ++it;
    check_right_generating_pairs(S, *it);
    REQUIRE(
        (sims::right_generating_pairs(*it) | rx::to_vector())
        == std::vector<relation_type>(
            {{1_w, 0_w},  {2_w, 0_w},  {3_w, 0_w},  {6_w, 0_w},  {00_w, 0_w},
             {01_w, 0_w}, {02_w, 0_w}, {03_w, 0_w}, {04_w, 4_w}, {05_w, 5_w},
             {06_w, 0_w}, {40_w, 4_w}, {41_w, 4_w}, {42_w, 4_w}, {43_w, 4_w},
             {44_w, 5_w}, {45_w, 0_w}, {46_w, 4_w}, {50_w, 5_w}, {51_w, 5_w},
             {52_w, 5_w}, {53_w, 5_w}, {54_w, 0_w}, {55_w, 4_w}, {56_w, 5_w}}));
    REQUIRE((sims::right_generating_pairs(p, *it) | rx::to_vector())
            == std::vector<relation_type>({{1_w, 0_w},
                                           {2_w, 0_w},
                                           {3_w, 0_w},
                                           {6_w, 0_w},
                                           {40_w, 4_w},
                                           {42_w, 4_w},
                                           {44_w, 5_w},
                                           {50_w, 5_w}}));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "003",
                          "ToddCoxeter failing example (std::string)",
                          "[quick][low-index]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.contains_empty_word(false);

    p.alphabet("aAbBcCe");
    presentation::add_identity_rules(p, 'e');

    presentation::add_inverse_rules(p, "AaBbCce", 'e');
    presentation::add_rule(p, "aaCac", "e");
    presentation::add_rule(p, "acbbACb", "e");
    presentation::add_rule(p, "ABabccc", "e");
    Sims1 S;
    S.presentation(p);
    REQUIRE(S.number_of_congruences(3) == 14);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "004",
                          "partition_monoid(2) onesided",
                          "[quick][low-index][no-valgrind]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.contains_empty_word(false);

    p.alphabet(0123_w);
    presentation::add_identity_rules(p, 0);
    presentation::add_rule(p, 11_w, 0_w);
    presentation::add_rule(p, 13_w, 3_w);
    presentation::add_rule(p, 22_w, 2_w);
    presentation::add_rule(p, 31_w, 3_w);
    presentation::add_rule(p, 33_w, 3_w);
    presentation::add_rule(p, 232_w, 2_w);
    presentation::add_rule(p, 323_w, 3_w);
    presentation::add_rule(p, 1212_w, 212_w);
    presentation::add_rule(p, 2121_w, 212_w);

    Sims1 S(p);
    REQUIRE(S.number_of_congruences(2) == 4);
    REQUIRE(S.number_of_congruences(3) == 7);
    REQUIRE(S.number_of_congruences(4) == 14);
    REQUIRE(S.number_of_congruences(5) == 23);
    REQUIRE(S.number_of_congruences(6) == 36);
    REQUIRE(S.number_of_congruences(7) == 51);
    REQUIRE(S.number_of_congruences(8) == 62);
    REQUIRE(S.number_of_congruences(9) == 74);
    REQUIRE(S.number_of_congruences(10) == 86);
    REQUIRE(S.number_of_congruences(11) == 95);
    REQUIRE(S.number_of_congruences(12) == 100);
    REQUIRE(S.number_of_congruences(13) == 102);
    REQUIRE(S.number_of_congruences(14) == 104);
    REQUIRE(S.number_of_congruences(15) == 105);
    REQUIRE(S.number_of_congruences(16) == 105);
    REQUIRE(S.number_of_congruences(17) == 105);

    Sims2 T;
    T.presentation(p);
    T.number_of_threads(2);
    REQUIRE(T.number_of_congruences(16) == 13);

    // sims::dot_poset("example-004", S.cbegin(16), S.cend(16));

    MinimalRepOrc orc;
    auto          d = orc.presentation(p)
                 .target_size(15)
                 .number_of_threads(std::thread::hardware_concurrency())
                 .word_graph();

    REQUIRE(d.number_of_nodes() == 7);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "005",
                          "partition_monoid(3)",
                          "[standard][low-index][no-coverage]") {
    auto rg = ReportGuard(false);
    auto p  = partition_monoid(3, author::Machine);
    REQUIRE(!p.contains_empty_word());
    REQUIRE(p.alphabet() == 01234_w);

    Sims1 S;
    S.presentation(p).long_rule_length(11).number_of_threads(
        4);  // This actually helps here!
    REQUIRE(S.number_of_congruences(17) == 1'589);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "006",
                          "full_transformation_monoid(3) onesided",
                          "[quick][low-index][no-valgrind]") {
    auto rg = ReportGuard(false);
    auto S  = to_froidure_pin({make<Transf<3>>({1, 2, 0}),
                               make<Transf<3>>({1, 0, 2}),
                               make<Transf<3>>({0, 1, 0})});
    REQUIRE(S.size() == 27);
    REQUIRE(S.number_of_generators() == 3);
    REQUIRE(S.number_of_rules() == 16);
    auto p = to_presentation<word_type>(S);
    REQUIRE(static_cast<size_t>(std::distance(p.rules.cbegin(), p.rules.cend()))
            == 2 * S.number_of_rules());
    Sims1 C;
    C.presentation(p);
    REQUIRE(C.number_of_congruences(27) == 287);
    // sims::dot_poset("example-006-T-3", C.cbegin(27), C.cend(27));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "007",
                          "full_transformation_monoid(3) left",
                          "[quick][low-index][no-valgrind]") {
    auto rg = ReportGuard(false);
    auto S  = to_froidure_pin(
        {Transf<3>({1, 2, 0}), Transf<3>({1, 0, 2}), Transf<3>({0, 1, 0})});
    REQUIRE(S.size() == 27);
    auto p = to_presentation<word_type>(S);
    presentation::reverse(p);
    Sims1 C(p);
    REQUIRE(C.number_of_congruences(27) == 120);
    // sims::dot_poset("example-006-T-3-left", C.cbegin(27), C.cend(27));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "008",
                          "full_transformation_monoid(4) left",
                          "[fail][low-index][babbage]") {
    auto rg = ReportGuard(true);

    auto p = full_transformation_monoid(4, author::Iwahori);

    REQUIRE(p.alphabet().size() == 4);

    auto w = presentation::longest_subword_reducing_length(p);
    while (!w.empty()) {
      presentation::replace_word_with_new_generator(p, w);
      w = presentation::longest_subword_reducing_length(p);
    }

    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    presentation::remove_duplicate_rules(p);
    presentation::reduce_complements(p);
    presentation::remove_trivial_rules(p);

    do {
      auto it = knuth_bendix::redundant_rule(p, std::chrono::milliseconds(100));
      p.rules.erase(it, it + 2);
    } while (presentation::length(p) > 700);

    Sims1 C;
    C.presentation(p);
    // Takes about 1h31m to run!
    REQUIRE(C.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(256)
            == 22'069'828);
    presentation::reverse(p);
    C.init(p);
    REQUIRE(C.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(256)
            == 120'121);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "009",
                          "symmetric_inverse_monoid(2, author::Gay)",
                          "[quick][low-index]") {
    auto  rg = ReportGuard(false);
    Sims1 C;
    C.presentation(symmetric_inverse_monoid(2, author::Gay));
    REQUIRE(C.number_of_congruences(7) == 10);  // Should be 10
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "010",
                          "symmetric_inverse_monoid(2) from FroidurePin",
                          "[quick][low-index]") {
    auto rg = ReportGuard(false);
    auto S  = to_froidure_pin({PPerm<2>({1, 0}), PPerm<2>({0}, {0}, 2)});
    REQUIRE(S.size() == 7);
    auto p = to_presentation<word_type>(S);
    presentation::reverse(p);
    Sims1 C(p);
    REQUIRE(C.number_of_congruences(7) == 10);
    // sims::dot_poset("example-010-I-2-left", C.cbegin(7), C.cend(7));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "011",
                          "symmetric_inverse_monoid(3)",
                          "[quick][low-index][no-valgrind]") {
    auto rg = ReportGuard(false);
    auto p  = symmetric_inverse_monoid(3, author::Gay);
    presentation::reverse(p);
    Sims1 C(p);
    REQUIRE(C.number_of_congruences(34) == 274);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "012",
                          "symmetric_inverse_monoid(4)",
                          "[extreme][low-index]") {
    auto p = symmetric_inverse_monoid(4, author::Gay);
    presentation::remove_duplicate_rules(p);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);

    REQUIRE(presentation::length(p) == 48);
    REQUIRE(p.alphabet().size() == 4);
    REQUIRE(*presentation::shortest_rule(p) == 00_w);
    REQUIRE(*(presentation::shortest_rule(p) + 1) == ""_w);
    REQUIRE(presentation::longest_rule_length(p) == 8);

    Sims1 C;
    C.presentation(p);
    sims::add_excluded_pair(C, ""_w, 11_w);
    REQUIRE(C.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(209)
            == 0);
    C.clear_excluded_pairs();

    auto rg = ReportGuard(true);
    REQUIRE(C.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(209)
            == 195'709);

    Sims2 S(p);
    REQUIRE(S.number_of_threads(1).number_of_congruences(209) == 11);
    REQUIRE(S.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(209)
            == 11);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "013",
                          "symmetric_inverse_monoid(5)",
                          "[fail][low-index]") {
    // This might take an extremely long time to terminate
    auto rg = ReportGuard(true);
    auto p  = symmetric_inverse_monoid(5, author::Gay);
    presentation::reverse(p);
    Sims1 C(p);
    // NOTE: Never ran to completion, there should be a non-zero number of
    // congruences.
    REQUIRE(C.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(1'546)
            == 0);
    // On 24/08/2022 JDM ran this for approx. 16 hours overnight on his
    // laptop, the last line of output was: #4: Sims1: found 63'968'999
    // congruences in 52156s! #21: Sims1: found 759'256'468 congruences in
    // 244617.546875 #12: Sims1: found 943'233'501 congruences in
    // 321019.531250! #30: Sims1: found 1'005'857'634 congruences in
    // 350411.000000! #45: Sims1: found 1'314'588'296 congruences in
    // 487405.000000! #20: Sims1: found 4'619'043'843 congruences in
    // 2'334'184.500000! #49: Sims1: found 5'582'499'404 congruences in
    // 2'912'877.000000! #10: Sims1: found 6'825'113'083 congruences in
    // 3705611.250000!
    //
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "014",
                          "temperley_lieb_monoid(3) from presentation",
                          "[quick][low-index]") {
    auto  rg = ReportGuard(false);
    auto  p  = temperley_lieb_monoid(3);
    Sims1 S(p);
    REQUIRE(S.number_of_congruences(14) == 9);
    // sims::dot_poset("example-014-TL-3-onesided", S.cbegin(14), S.cend(14));
    presentation::reverse(p);
    S.init(p);
    REQUIRE(S.number_of_congruences(14) == 9);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "015",
                          "temperley_lieb_monoid(4) from presentation",
                          "[quick][low-index]") {
    auto  rg = ReportGuard(false);
    auto  p  = temperley_lieb_monoid(4);
    Sims1 S(p);
    REQUIRE(S.number_of_congruences(14) == 79);
    // sims::dot_poset("example-014-TL-4-onesided", S.cbegin(14), S.cend(14));
    presentation::reverse(p);
    S.init(p);
    REQUIRE(S.number_of_congruences(14) == 79);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "016",
                          "fp semigroup containing given pairs #1",
                          "[quick][low-index]") {
    auto rg = ReportGuard(false);

    Presentation<word_type> p;
    p.contains_empty_word(true);

    p.alphabet({0, 1});
    presentation::add_rule(p, 000_w, 0_w);
    presentation::add_rule(p, 11_w, 1_w);
    presentation::add_rule(p, 0101_w, 0_w);

    Sims1 S;
    S.presentation(p);
    sims::add_included_pair(S, 0_w, 1_w);
    REQUIRE(S.number_of_congruences(5) == 2);
    check_include(p, {0_w, 1_w}, 5);

    sims::add_excluded_pair(S, 0_w, 1_w);
    S.clear_included_pairs();
    REQUIRE(S.pruners().size() == 1);
    REQUIRE(S.number_of_congruences(5) == 4);
    S.clear_excluded_pairs();
    REQUIRE(S.number_of_congruences(5) == 6);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "017",
                          "fp semigroup containing given pairs #2",
                          "[quick][low-index]") {
    auto rg = ReportGuard(false);

    Presentation<word_type> p;
    p.contains_empty_word(true);

    p.alphabet({0, 1});
    presentation::add_rule(p, 000_w, 0_w);
    presentation::add_rule(p, 11_w, 1_w);
    presentation::add_rule(p, 0101_w, 0_w);
    Sims1 T;
    T.presentation(p);
    sims::add_included_pair(T, 01_w, 1_w);
    REQUIRE(T.number_of_congruences(5) == 2);
    check_include(p, T.included_pairs(), 5);
    presentation::reverse(p);
    T.init(p);
    sims::add_included_pair(T, 10_w, 1_w);
    REQUIRE(T.number_of_congruences(5) == 2);
    check_include(p, T.included_pairs(), 5);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "018",
                          "fp semigroup containing given pairs #3",
                          "[quick][low-index]") {
    auto rg = ReportGuard(false);

    Presentation<word_type> p;
    p.contains_empty_word(true);
    p.alphabet({0, 1});
    presentation::add_rule(p, 000_w, 0_w);
    presentation::add_rule(p, 11_w, 1_w);
    presentation::add_rule(p, 0101_w, 0_w);

    {
      Sims1 T;
      T.presentation(p);
      sims::add_included_pair(T, 0101_w, 0_w);
      REQUIRE(T.number_of_congruences(5) == 6);
    }
    {
      Sims1 T;
      presentation::reverse(p);
      T.presentation(p);
      sims::add_included_pair(T, 0101_w, 0_w);
      REQUIRE(T.included_pairs() == std::vector<word_type>({0101_w, 0_w}));
      REQUIRE(T.number_of_congruences(5) == 4);  // Verified with GAP
    }
    check_include(p, {0101_w, 0_w}, 5);
    presentation::reverse(p);
    check_include(p, {1010_w, 0_w}, 5);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "019",
                          "ToddCoxeter failing example",
                          "[quick][low-index]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.contains_empty_word(false);

    p.alphabet("aAbBcCe");
    presentation::add_identity_rules(p, 'e');

    presentation::add_inverse_rules(p, "AaBbCce", 'e');
    presentation::add_rule(p, "aaCac", "e");
    presentation::add_rule(p, "acbbACb", "e");
    presentation::add_rule(p, "ABabccc", "e");

    Sims1  S;
    ToWord to_word(p.alphabet());
    S.presentation(p);
    sims::add_included_pair(S, to_word("a"), to_word("A"));
    sims::add_included_pair(S, to_word("a"), to_word("b"));
    REQUIRE(S.number_of_congruences(3) == 2);

    check_include(S.presentation(), S.included_pairs(), 3);
    presentation::reverse(p);
    S.presentation(p);
    check_include(S.presentation(), S.included_pairs(), 3);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "020",
                          "fp example 2 (check_include)",
                          "[quick][low-index]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.contains_empty_word(true);

    p.alphabet({0, 1, 2});
    presentation::add_rule(p, 010_w, 00_w);
    presentation::add_rule(p, 22_w, 00_w);
    presentation::add_rule(p, 000_w, 00_w);
    presentation::add_rule(p, 21_w, 12_w);
    presentation::add_rule(p, 20_w, 00_w);
    presentation::add_rule(p, 11_w, 1_w);
    presentation::add_rule(p, 02_w, 00_w);

    std::vector<word_type> e = {1_w, 00_w};
    check_include(p, e, 11);
    presentation::reverse(p);
    check_include(p, e, 11);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1", "021", "exceptions", "[quick][low-index]") {
    Presentation<word_type> p;
    p.alphabet({0, 1, 2});
    presentation::add_rule(p, {0, 1, 0}, {0, 0});

    Presentation<word_type> e;
    e.alphabet({0, 1});
    REQUIRE_NOTHROW(Sims1().presentation(p).presentation(e));
    REQUIRE_NOTHROW(Sims1());
    Sims1 S;
    REQUIRE_THROWS_AS(S.number_of_threads(0), LibsemigroupsException);
    RepOrc ro;
    REQUIRE_THROWS_AS(ro.number_of_threads(0), LibsemigroupsException);
    MinimalRepOrc mro;
    REQUIRE_THROWS_AS(mro.number_of_threads(0), LibsemigroupsException);
    REQUIRE_NOTHROW(S.presentation(p));
    REQUIRE_THROWS_AS(S.cbegin_long_rules(S.presentation().rules.cend() + 1),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(S.cbegin_long_rules(S.presentation().rules.cbegin() - 1),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(S.cbegin_long_rules(S.presentation().rules.cbegin() + 1),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(S.idle_thread_restarts(0), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "022",
                          "singular_brauer_monoid(4) (Maltcev-Mazorchuk)",
                          "[standard][sims1][no-coverage]") {
    auto                     rg = ReportGuard(true);
    FroidurePin<Bipartition> S;
    S.add_generator(Bipartition({{1, 2}, {3, -1}, {4, -2}, {-3, -4}}));
    S.add_generator(Bipartition({{1, 2}, {3, -1}, {4, -4}, {-2, -3}}));
    S.add_generator(Bipartition({{1, 2}, {3, -3}, {4, -1}, {-2, -4}}));
    S.add_generator(Bipartition({{1, 2}, {3, -2}, {4, -3}, {-1, -4}}));
    S.add_generator(Bipartition({{1, 2}, {3, -2}, {4, -4}, {-1, -3}}));
    S.add_generator(Bipartition({{1, 3}, {2, -4}, {4, -3}, {-1, -2}}));
    S.add_generator(Bipartition({{1, -4}, {2, 3}, {4, -3}, {-1, -2}}));
    S.add_generator(Bipartition({{1, 4}, {2, -3}, {3, -4}, {-1, -2}}));
    S.add_generator(Bipartition({{1, -3}, {2, 4}, {3, -4}, {-1, -2}}));
    S.add_generator(Bipartition({{1, -3}, {2, -4}, {3, 4}, {-1, -2}}));
    REQUIRE(S.size() == 81);

    /// auto p  = singular_brauer_monoid(4);
    auto p = to_presentation<word_type>(S);
    p.validate();
    REQUIRE(p.alphabet().size() == 10);
    REQUIRE(presentation::length(p) == 719);

    // REQUIRE(*presentation::shortest_rule(p) == 0_w);
    // REQUIRE(*(presentation::shortest_rule(p) + 1) == 3_w);

    // presentation::remove_redundant_generators(p);

    // presentation::remove_duplicate_rules(p);
    // presentation::sort_each_rule(p);
    // presentation::sort_rules(p);

    // REQUIRE(presentation::shortest_rule_length(p) == 3);
    // REQUIRE(*presentation::shortest_rule(p) == 00_w);
    // REQUIRE(*(presentation::shortest_rule(p) + 1) == 0_w);

    // REQUIRE(presentation::longest_rule_length(p) == 6);
    // REQUIRE(*presentation::longest_rule(p) == "048"_w);
    // REQUIRE(*(presentation::longest_rule(p) + 1) == "028"_w);

    // REQUIRE(p.alphabet().size() == 6);
    // REQUIRE(presentation::length(p) == 462);
    // REQUIRE(p.rules.size() == 186);
    // REQUIRE(p.rules == std::vector<word_type>());

    // p.contains_empty_word(true);

    std::vector<word_type> forbid = {{0},
                                     {3, 0},
                                     {0, 0},
                                     {0, 1},
                                     {0, 0},
                                     {0, 2},
                                     {0, 2},
                                     {0, 1},
                                     {0, 0},
                                     {5, 9},
                                     {0, 0},
                                     {6, 9},
                                     {5, 9},
                                     {6, 9}};

    auto filter = [&forbid](auto const& wg) {
      auto first = forbid.cbegin();
      auto last  = forbid.cend();
      for (auto it = first; it != last; it += 2) {
        bool this_rule_compatible = true;
        for (auto n : wg.nodes()) {
          auto l = word_graph::follow_path_no_checks(wg, n, *it);
          auto r = word_graph::follow_path_no_checks(wg, n, *(it + 1));
          if (l != r) {
            this_rule_compatible = false;
            break;
          }
        }
        if (this_rule_compatible) {
          return false;
        }
      }
      return true;
    };

    Sims1 sims;
    sims.presentation(p);
    {
      SimsRefinerFaithful pruno(forbid);

      auto const& wg = sims.number_of_threads(1).add_pruner(pruno).find_if(
          82, [](auto const&) { return true; });
      REQUIRE(wg.number_of_active_nodes() == 18);
      REQUIRE(sims.stats().total_pending_now == 370'719);

      auto const& wg2 = sims.find_if(wg.number_of_active_nodes() - 2,
                                     [](auto const&) { return true; });
      REQUIRE(sims.stats().total_pending_now == 1'930'725);
      REQUIRE(wg2.number_of_active_nodes() == 0);
    }
    {
      sims.clear_pruners();

      auto const& wg = sims.number_of_threads(1).find_if(82, filter);
      REQUIRE(wg.number_of_active_nodes() == 18);
      REQUIRE(sims.stats().total_pending_now == 1'014'357);
      auto const& wg2 = sims.find_if(wg.number_of_active_nodes() - 2, filter);
      REQUIRE(sims.stats().total_pending_now == 3'374'651);
      REQUIRE(wg2.number_of_active_nodes() == 0);
    }
    {
      MinimalRepOrc orc;
      p.contains_empty_word(true);
      auto d
          = orc.presentation(p)
                .target_size(82)
                //                 .number_of_threads(std::thread::hardware_concurrency())
                .word_graph();
      REQUIRE(d.number_of_nodes() == 18);
      REQUIRE(orc.target_size() == 82);
      REQUIRE(orc.stats().total_pending_now == 3'626'612);
    }
    {
      MinimalRepOrc       orc;
      SimsRefinerFaithful pruno(forbid);
      p.contains_empty_word(true);
      auto d
          = orc.add_pruner(pruno)
                .presentation(p)
                .target_size(82)
                //                 .number_of_threads(std::thread::hardware_concurrency())
                .word_graph();
      REQUIRE(d.number_of_nodes() == 18);
      REQUIRE(orc.target_size() == 82);
      REQUIRE(orc.stats().total_pending_now == 2'074'472);
    }

    //   // p.contains_empty_word(false);

    //   // Sims1 C;
    //   // C.presentation(p);
    //   // REQUIRE(C.presentation().rules.size() == 186);

    //   // REQUIRE(C.number_of_threads(std::thread::hardware_concurrency())
    //   //             .number_of_congruences(81)
    //   //         == 601'265);

    //   // Sims2 sim2(p);
    //   // REQUIRE(C.number_of_threads(1).number_of_congruences(81) == 75);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "023",
                          "Brauer(4) from FroidurePin",
                          "[extreme][sims1]") {
    auto                     rg = ReportGuard(true);
    FroidurePin<Bipartition> S;
    S.add_generator(Bipartition({{1, -1}, {2, -2}, {3, -3}, {4, -4}}));
    S.add_generator(Bipartition({{1, -2}, {2, -3}, {3, -4}, {4, -1}}));
    S.add_generator(Bipartition({{1, -2}, {2, -1}, {3, -3}, {4, -4}}));
    S.add_generator(Bipartition({{1, 2}, {3, -3}, {4, -4}, {-1, -2}}));
    REQUIRE(S.size() == 105);

    auto p = to_presentation<word_type>(S);
    REQUIRE(presentation::length(p) == 359);
    presentation::remove_duplicate_rules(p);
    REQUIRE(presentation::length(p) == 359);
    presentation::reduce_complements(p);
    REQUIRE(presentation::length(p) == 359);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    REQUIRE(p.rules.size() == 86);

    Sims1 C;
    C.presentation(p).long_rule_length(12);
    REQUIRE(C.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(105)
            == 103'406);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "024",
                          "brauer_monoid(4) (Kudryavtseva-Mazorchuk)",
                          "[extreme][sims1]") {
    auto rg = ReportGuard(true);
    auto p  = brauer_monoid(4);
    REQUIRE(p.alphabet().size() == 6);
    REQUIRE(presentation::length(p) == 140);
    presentation::remove_duplicate_rules(p);
    REQUIRE(presentation::length(p) == 120);
    presentation::reduce_complements(p);
    REQUIRE(presentation::length(p) == 117);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    REQUIRE(p.rules.size() == 60);

    // auto d = MinimalRepOrc().presentation(p).target_size(105).word_graph();
    // REQUIRE(d.number_of_nodes() == 22);
    // REQUIRE(word_graph::is_strictly_cyclic(d));
    // REQUIRE(
    //     d
    //     == make<WordGraph<uint32_t>>(
    //         22, {{0, 0, 1, 0, 2, 3, 2},        {1, 4, 0, 5, 6, 3, 7},
    //              {2, 2, 2, 2, 2, 2, 2},        {3, 8, 3, 9, 6, 3, 7},
    //              {4, 1, 4, 10, 6, 2, 11},      {5, 10, 5, 1, 12, 2, 7},
    //              {6, 6, 8, 12, 6, 3, 13},      {7, 11, 9, 7, 13, 3, 7},
    //              {8, 3, 6, 14, 6, 3, 11},      {9, 14, 7, 3, 12, 3, 7},
    //              {10, 5, 15, 4, 12, 16, 11},   {11, 7, 17, 11, 13, 16, 11},
    //              {12, 12, 18, 6, 12, 16, 13},  {13, 13, 19, 13, 13, 20,
    //              13}, {14, 9, 21, 8, 12, 20, 11},   {15, 15, 10, 15, 2, 16,
    //              2}, {16, 18, 16, 17, 12, 16, 11}, {17, 21, 11, 16, 6, 16,
    //              11}, {18, 16, 12, 21, 12, 16, 7},  {19, 20, 13, 20, 13,
    //              20, 13}, {20, 19, 20, 19, 13, 20, 13}, {21, 17, 14, 18, 6,
    //              20, 7}}));

    // auto S = to_froidure_pin<Transf<0, node_type>>(d);
    // REQUIRE(S.size() == 105);
    // REQUIRE(S.generator(0) == Transf<0, node_type>::identity(22));
    // REQUIRE(
    //     S.generator(1)
    //     == Transf<0, node_type>({0, 4,  2,  8, 1,  10, 6,  11, 3,  14, 5,
    //                              7, 12, 13, 9, 15, 18, 21, 16, 20, 19,
    //                              17}));
    // REQUIRE(
    //     S.generator(2)
    //     == Transf<0, node_type>({1,  0,  2,  3,  4,  5,  8,  9,  6,  7, 15,
    //                              17, 18, 19, 21, 10, 16, 11, 12, 13, 20,
    //                              14}));
    // REQUIRE(S.generator(3)
    //         == Transf<0, node_type>({0, 5,  2, 9,  10, 1,  12, 7,  14, 3,
    //         4, 11,
    //                                  6, 13, 8, 15, 17, 16, 21, 20, 19,
    //                                  18}));
    // REQUIRE(S.generator(4)
    //         == Transf<0, node_type>({2,  6,  2,  6,  6, 12, 6, 13, 6,  12,
    //         12,
    //                                  13, 12, 13, 12, 2, 12, 6, 12, 13, 13,
    //                                  6}));
    // REQUIRE(
    //     S.generator(5)
    //     == Transf<0, node_type>({3,  3,  2,  3,  2,  2,  3,  3,  3,  3, 16,
    //                              16, 16, 20, 20, 16, 16, 16, 16, 20, 20,
    //                              20}));
    // REQUIRE(
    //     S.generator(6)
    //     == Transf<0, node_type>({2,  7,  2,  7,  11, 7,  13, 7, 11, 7,  11,
    //                              11, 13, 13, 11, 2,  11, 11, 7, 13, 13,
    //                              7}));

    Sims1 C;
    C.presentation(p);
    REQUIRE(C.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(105)
            == 103'406);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "025",
                          "brauer_monoid(5) (Kudryavtseva-Mazorchuk)",
                          "[extreme][sims1]") {
    auto rg = ReportGuard(true);

    // FroidurePin<Bipartition> S;
    // S.add_generator(Bipartition({{1, -1}, {2, -2}, {3, -3}, {4, -4}, {5,
    // -5}})); S.add_generator(Bipartition({{1, -2}, {2, -3}, {3, -4}, {4,
    // -5}, {5, -1}})); S.add_generator(Bipartition({{1, -2}, {2, -1}, {3,
    // -3}, {4, -4}, {5, -5}})); S.add_generator(Bipartition({{1, 2}, {3, -3},
    // {4, -4}, {5, -5}, {-1, -2}})); REQUIRE(S.size() == 945);

    // auto p = to_presentation<word_type>(S);
    auto p = brauer_monoid(5);

    REQUIRE(presentation::length(p) == 240);

    presentation::remove_duplicate_rules(p);
    presentation::reduce_complements(p);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    REQUIRE(presentation::length(p) == 194);

    REQUIRE(presentation::shortest_rule_length(p) == 2);
    REQUIRE(*presentation::shortest_rule(p) == 00_w);
    REQUIRE(*(presentation::shortest_rule(p) + 1) == ""_w);

    REQUIRE(p.alphabet().size() == 8);

    presentation::remove_redundant_generators(p);
    REQUIRE(p.alphabet() == 01234567_w);
    REQUIRE(p.alphabet().size() == 8);
    REQUIRE(presentation::length(p) == 194);

    REQUIRE(*presentation::longest_rule(p) == 101_w);
    REQUIRE(*(presentation::longest_rule(p) + 1) == 010_w);

    REQUIRE(presentation::longest_subword_reducing_length(p) == 76_w);
    presentation::replace_word_with_new_generator(p, 76_w);
    REQUIRE(presentation::length(p) == 193);

    std::vector<word_type> forbid =  // {{3, 1, 2, 1, 1, 3},
                                     //  {1, 1, 3, 1, 2, 1, 1, 3},
                                     //  {1, 1, 2, 1, 3, 1, 1, 3, 1},
                                     //  {1, 1, 2, 1, 3, 1, 1, 3}};
        {{4, 8}, {4, 5, 7}, {4, 7}, {5, 4, 7}};
    Sims1 sims;
    sims.presentation(p);
    // TODO(2) use SimsRefinerFaithful instead
    for (auto it = forbid.cbegin(); it != forbid.cend(); it += 2) {
      sims::add_excluded_pair(sims, *it, *(it + 1));
    }

    Sims1::word_graph_type wg
        = sims.number_of_threads(std::thread::hardware_concurrency())
              .find_if(945, [](auto const&) { return true; });

    // REQUIRE(wg.number_of_active_nodes() == expected[index++]);
    while (wg.number_of_active_nodes() != 0) {
      wg = sims.number_of_threads(std::thread::hardware_concurrency())
               .find_if(wg.number_of_active_nodes() - 1,
                        [](auto const&) { return true; });
      fmt::print("NUMBER OF NODES IS {}\n", wg.number_of_active_nodes());
    }
    REQUIRE(wg.number_of_active_nodes() == 0);

    // auto d = MinimalRepOrc()
    //              .presentation(p)
    //              .target_size(945)
    //              .number_of_threads(8)
    //              .word_graph();
    // // WARNING: the number below is not necessarily the minimal degree of
    // an
    // // action on onesided congruences, only the minimal degree of an action
    // on
    // // onesided congruences containing the pair {0}, {1}.
    // REQUIRE(d.number_of_nodes() == 51);
    // auto S = to_froidure_pin<Transf<0, node_type>>(d);
    // REQUIRE(S.size() == 945);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "026",
                          "uniform_block_bijection_monoid(4) (Fitzgerald)",
                          "[extreme][sims1]") {
    auto rg = ReportGuard(true);
    auto p  = uniform_block_bijection_monoid(4, author::FitzGerald);
    presentation::remove_duplicate_rules(p);
    presentation::reduce_complements(p);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);

    Sims1 C;
    C.presentation(p);
    REQUIRE(C.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(131)
            == 280'455);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "Sims1",
      "027",
      "from https://mathoverflow.net/questions/423541/ (semigroup)",
      "[quick][sims1][no-valgrind]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.contains_empty_word(false);
    p.alphabet("aAbBe");
    presentation::add_identity_rules(p, 'e');
    presentation::add_inverse_rules(p, "AaBbe", 'e');
    presentation::add_rule(p, "aaa", "e");
    presentation::add_rule(p, "baBBBABA", "e");
    Sims1 C;
    C.presentation(p);
    REQUIRE(C.number_of_congruences(10) == 3);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "Sims1",
      "028",
      "from https://mathoverflow.net/questions/423541/ (monoid)",
      "[quick][sims1][no-valgrind]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("aAbB");
    presentation::add_inverse_rules(p, "AaBb");
    presentation::add_rule(p, "aaa", "");
    presentation::add_rule(p, "baBBBABA", "");
    Sims1 C;
    C.presentation(p);
    REQUIRE(C.number_of_congruences(10) == 3);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "029",
                          "fibonacci_semigroup(4, 6)",
                          "[standard][sims1][no-valgrind]") {
    std::cout << "\n";            // So that the reporting looks good
    auto rg = ReportGuard(true);  // for code coverage
    auto p  = fibonacci_semigroup(4, 6);
    presentation::remove_duplicate_rules(p);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    REQUIRE(presentation::length(p) == 30);
    REQUIRE(p.rules.size() == 12);
    REQUIRE(p.rules[0].size() + p.rules[1].size() == 5);

    Sims1 C;
    C.presentation(p);
    REQUIRE(C.number_of_congruences(3) == 5);

    C.number_of_threads(2);
    REQUIRE(C.number_of_congruences(3) == 5);
    REQUIRE_THROWS_AS(C.find_if(0, [](auto const&) { return true; }),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "030",
                          "presentation with one free generator",
                          "[quick][sims1]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(4);
    presentation::add_rule_no_checks(p, 121_w, 11_w);
    presentation::add_rule_no_checks(p, 33_w, 11_w);
    presentation::add_rule_no_checks(p, 111_w, 11_w);
    presentation::add_rule_no_checks(p, 32_w, 23_w);
    presentation::add_rule_no_checks(p, 31_w, 11_w);
    presentation::add_rule_no_checks(p, 22_w, 2_w);
    presentation::add_rule_no_checks(p, 13_w, 11_w);
    p.validate();
    Sims1 C;
    C.presentation(p);
    REQUIRE(C.number_of_congruences(2) == 67);
    // sims::dot_poset("example-030-onesided", C.cbegin(2), C.cend(2));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "031",
                          "presentation with non-zero index generators",
                          "[quick][sims1]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    presentation::add_rule_no_checks(p, 121_w, 11_w);
    presentation::add_rule_no_checks(p, 33_w, 11_w);
    presentation::add_rule_no_checks(p, 111_w, 11_w);
    presentation::add_rule_no_checks(p, 32_w, 23_w);
    presentation::add_rule_no_checks(p, 31_w, 11_w);
    presentation::add_rule_no_checks(p, 22_w, 2_w);
    presentation::add_rule_no_checks(p, 13_w, 11_w);
    p.alphabet_from_rules();
    p.validate();

    Sims1 C;
    C.presentation(p);
    REQUIRE(C.number_of_congruences(2) == 7);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "032",
                          "presentation with empty word",
                          "[quick][sims1][no-valgrind]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.contains_empty_word(true);

    //         aAbBcC
    p.alphabet(012345_w);
    presentation::add_inverse_rules(p, 103254_w);
    presentation::add_rule(p, 00504_w, {});
    presentation::add_rule(p, 0422152_w, {});
    presentation::add_rule(p, 1302444_w, {});
    Sims1 S;
    S.presentation(p);
    REQUIRE(S.number_of_congruences(3) == 14);
    REQUIRE(S.number_of_congruences(4) == 14);
    REQUIRE(S.number_of_congruences(5) == 14);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1", "033", "constructors", "[quick][sims1]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p, pp;
    p.contains_empty_word(true);

    //         aAbBcC
    p.alphabet(012345_w);
    presentation::add_inverse_rules(p, 103254_w);
    presentation::add_rule(p, 00504_w, {});
    presentation::add_rule(p, 0422152_w, {});
    presentation::add_rule(p, 1302444_w, {});
    Sims1 S;
    S.presentation(p);

    Sims1 T(S);
    REQUIRE(S.number_of_congruences(3) == 14);
    REQUIRE(T.number_of_congruences(3) == 14);

    Sims1 U(std::move(S));
    REQUIRE(U.number_of_congruences(3) == 14);
    REQUIRE(T.number_of_congruences(3) == 14);

    S = U;
    REQUIRE(S.number_of_congruences(3) == 14);

    S = std::move(U);
    REQUIRE(S.number_of_congruences(3) == 14);

    Sims1 C;
    REQUIRE_THROWS_AS(sims::add_included_pair(C.presentation(p), 0127_w, 0_w),
                      LibsemigroupsException);

    Sims1 SP(p);
    REQUIRE(SP.number_of_congruences(3) == 14);

    pp = p;
    Sims1 SPP(std::move(pp));
    REQUIRE(SPP.number_of_congruences(3) == 14);

    Sims2 S2;
    S2.presentation(p);

    Sims2 T2(S2);
    REQUIRE(S2.number_of_congruences(3) == 14);
    REQUIRE(T2.number_of_congruences(3) == 14);

    Sims2 U2(std::move(S2));
    REQUIRE(U2.number_of_congruences(3) == 14);
    REQUIRE(T2.number_of_congruences(3) == 14);

    S2 = U2;
    REQUIRE(S2.number_of_congruences(3) == 14);

    S2 = std::move(U2);
    REQUIRE(S2.number_of_congruences(3) == 14);

    Sims2 C2;
    REQUIRE_THROWS_AS(sims::add_included_pair(C2.presentation(p), 0127_w, 0_w),
                      LibsemigroupsException);

    Sims2 SP2(p);
    REQUIRE(SP2.number_of_congruences(3) == 14);

    pp = p;
    Sims2 SPP2(std::move(pp));
    REQUIRE(SPP2.number_of_congruences(3) == 14);

    Presentation<word_type> q;
    q.alphabet({0, 1});
    q.contains_empty_word(true);
    presentation::add_rule(q, 000_w, 0_w);
    presentation::add_rule(q, 111_w, ""_w);
    presentation::add_rule(q, 011_w, 10_w);

    RepOrc Ro;
    REQUIRE(Ro.presentation(q)
                .target_size(9)
                .min_nodes(2)
                .max_nodes(6)
                .number_of_threads(4)
                .word_graph()
                .number_of_active_nodes()
            == 6);
    RepOrc Ro2(Ro);
    REQUIRE(Ro2.word_graph().number_of_active_nodes() == 6);
    RepOrc Ro3;
    Ro3 = Ro2;
    REQUIRE(Ro3.word_graph().number_of_active_nodes() == 6);
    RepOrc Ro4(std::move(Ro3));
    REQUIRE(Ro4.word_graph().number_of_active_nodes() == 6);
    RepOrc Ro5;
    Ro5 = std::move(Ro4);
    REQUIRE(Ro5.word_graph().number_of_active_nodes() == 6);

    MinimalRepOrc Mro;
    REQUIRE(Mro.presentation(q)
                .target_size(9)
                .number_of_threads(4)
                .word_graph()
                .number_of_active_nodes()
            == 6);
    MinimalRepOrc Mro2(Mro);
    REQUIRE(Mro2.word_graph().number_of_active_nodes() == 6);
    MinimalRepOrc Mro3;
    Mro3 = Mro2;
    REQUIRE(Mro3.word_graph().number_of_active_nodes() == 6);
    MinimalRepOrc Mro4(std::move(Mro3));
    REQUIRE(Mro4.word_graph().number_of_active_nodes() == 6);
    MinimalRepOrc Mro5;
    Mro5 = std::move(Mro4);
    REQUIRE(Mro5.word_graph().number_of_active_nodes() == 6);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "034",
                          "cbegin_long_rules",
                          "[standard][sims1]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.contains_empty_word(true);

    //         aAbBcC
    p.alphabet(012345_w);
    presentation::add_inverse_rules(p, 103254_w);
    presentation::add_rule(p, 00504_w, {});
    presentation::add_rule(p, 0422152_w, {});
    presentation::add_rule(p, 1302444_w, {});
    // F := FreeGroup("a", "b", "c");
    // AssignGeneratorVariables(F);
    // R := [ a*a*c^-1*a*c, a*c*b*b*a^-1*c^-1*b, a^-1*b^-1*a*b*c*c*c];
    // G := F/ R;
    // LowIndexSubgroups(G, 3); # returns 14
    Sims1 S;
    S.presentation(p);
    REQUIRE(S.number_of_congruences(3) == 14);
    REQUIRE(S.number_of_congruences(4) == 14);
    // sims::dot_poset("example-034-onesided", S.cbegin(3), S.cend(3));

    REQUIRE_THROWS_AS(S.cbegin_long_rules(p.rules.size() + 1),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(S.cbegin_long_rules(9), LibsemigroupsException);
    REQUIRE_NOTHROW(S.cbegin_long_rules(0));

    REQUIRE(!S.presentation().rules.empty());

    auto const& rules = S.presentation().rules;

    for (size_t i = 0; i <= rules.size() / 2; ++i) {
      S.cbegin_long_rules(2 * i);
      REQUIRE(static_cast<size_t>(
                  std::distance(rules.cbegin(), S.cbegin_long_rules()))
              == 2 * i);
    }
    REQUIRE(S.presentation().rules.size() == p.rules.size());
    for (size_t i = p.rules.size() / 2; i > 0; --i) {
      S.cbegin_long_rules(2 * i);
      REQUIRE(static_cast<size_t>(
                  std::distance(rules.cbegin(), S.cbegin_long_rules()))
              == 2 * i);
    }

    Presentation<word_type> q;
    q.contains_empty_word(p.contains_empty_word()).alphabet(p.alphabet());
    q.rules.insert(q.rules.end(), p.rules.begin(), p.rules.begin() + 8);
    q.validate();
    REQUIRE(q.alphabet() == 012345_w);
    REQUIRE(
        q.rules
        == std::vector<word_type>({01_w, {}, 10_w, {}, 23_w, {}, 32_w, {}}));
    q.validate();

    std::atomic_uint64_t num = 0;
    Sims1                T;
    T.presentation(q);

    REQUIRE(T.number_of_long_rules() == 0);
    T.for_each(3, [&](auto const& wg) {
      num += word_graph::is_compatible_no_checks(
          wg,
          wg.cbegin_nodes(),
          wg.cbegin_nodes() + wg.number_of_active_nodes(),
          p.rules.cbegin(),
          p.rules.cend());
    });
    REQUIRE(num == 14);  // 14 is the correct value

    num = 0;

    S.number_of_threads(1);
    REQUIRE(S.number_of_congruences(3) == 14);
    REQUIRE(rules.size() == 18);
    S.for_each(3, [&](auto const& wg) {
      REQUIRE(wg.out_degree() == 6);
      num += word_graph::is_compatible_no_checks(
                 wg,
                 wg.cbegin_nodes(),
                 wg.cbegin_nodes() + wg.number_of_active_nodes(),
                 rules.cbegin(),
                 S.cbegin_long_rules())
             && word_graph::is_complete(wg,
                                        wg.cbegin_nodes(),
                                        wg.cbegin_nodes()
                                            + wg.number_of_active_nodes());
    });
    REQUIRE(S.presentation().rules == p.rules);
    REQUIRE(num == 14);
    S.clear_long_rules();
    num = 0;
    REQUIRE(rules.size() == 18);
    REQUIRE(S.number_of_threads() == 1);
    REQUIRE(S.presentation().rules == p.rules);
    REQUIRE(S.number_of_congruences(3) == 14);
    S.for_each(3, [&](auto const& wg) {
      REQUIRE(wg.out_degree() == 6);
      num += word_graph::is_compatible_no_checks(wg,
                                                 wg.cbegin_nodes(),
                                                 wg.cend_nodes(),
                                                 rules.cbegin(),
                                                 S.cbegin_long_rules())
             && word_graph::is_complete(wg,
                                        wg.cbegin_nodes(),
                                        wg.cbegin_nodes()
                                            + wg.number_of_active_nodes());
    });
    REQUIRE(num == 14);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1", "035", "stats", "[quick][sims1]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.contains_empty_word(true);

    //         aAbBcC
    p.alphabet(012345_w);
    presentation::add_inverse_rules(p, 103254_w);
    presentation::add_rule(p, 00504_w, {});
    presentation::add_rule(p, 0422152_w, {});
    presentation::add_rule(p, 1302444_w, {});
    Sims1 S;
    S.presentation(p);

    REQUIRE(S.number_of_congruences(2) == 1);
    REQUIRE(S.stats().max_pending != 0);

    auto it = S.cbegin(2);
    ++it;
    REQUIRE(it.stats().max_pending != 0);

    Sims2 S2;
    S2.presentation(p);

    REQUIRE(S2.number_of_congruences(2) == 1);
    REQUIRE(S2.stats().max_pending != 0);

    auto it2 = S2.cbegin(2);
    ++it2;
    REQUIRE(it2.stats().max_pending != 0);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "036",
                          "check iterator requirements",
                          "[quick][sims1]") {
    Presentation<word_type> p;
    p.contains_empty_word(true);
    p.alphabet(01_w);
    presentation::add_rule(p, 000_w, 0_w);
    presentation::add_rule(p, 11_w, 1_w);
    presentation::add_rule(p, 0101_w, 0_w);

    Sims1 S;
    S.presentation(p);
    verify_forward_iterator_requirements(S.cbegin(10));
    auto it = S.cbegin(10);
    REQUIRE(it->number_of_nodes() == 10);
    REQUIRE(&it.sims() == &S);
    REQUIRE(it.maximum_number_of_classes() == 10);

    presentation::reverse(p);
    S.init(p);
    verify_forward_iterator_requirements(S.cbegin(10));
    it = S.cbegin(10);
    REQUIRE(it->number_of_nodes() == 10);
    REQUIRE(&it.sims() == &S);
    REQUIRE(it.maximum_number_of_classes() == 10);
    auto itc(std::move(it));
    REQUIRE(itc->number_of_nodes() == 10);
    REQUIRE(&itc.sims() == &S);
    REQUIRE(itc.maximum_number_of_classes() == 10);
    Sims1::iterator itcc;
    itcc = itc;
    REQUIRE(itcc->number_of_nodes() == 10);
    REQUIRE(&itcc.sims() == &S);
    REQUIRE(itcc.maximum_number_of_classes() == 10);

    Sims2 S2;
    S2.presentation(p);
    verify_forward_iterator_requirements(S2.cbegin(10));
    auto it2 = S2.cbegin(10);
    REQUIRE(it2->number_of_nodes() == 10);
    REQUIRE(&it2.sims() == &S2);
    REQUIRE(it2.maximum_number_of_classes() == 10);

    presentation::reverse(p);
    S2.init(p);
    verify_forward_iterator_requirements(S2.cbegin(10));
    it2 = S2.cbegin(10);
    REQUIRE(it2->number_of_nodes() == 10);
    REQUIRE(&it2.sims() == &S2);
    REQUIRE(it2.maximum_number_of_classes() == 10);
    auto itc2(std::move(it2));
    REQUIRE(itc2->number_of_nodes() == 10);
    REQUIRE(&itc2.sims() == &S2);
    REQUIRE(itc2.maximum_number_of_classes() == 10);
    Sims2::iterator itcc2;
    itcc2 = itc2;
    REQUIRE(itcc2->number_of_nodes() == 10);
    REQUIRE(&itcc2.sims() == &S2);
    REQUIRE(itcc2.maximum_number_of_classes() == 10);
  }

  // Takes about 4s
  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "037",
                          "rectangular_band(9, 2)",
                          "[extreme][sims1]") {
    auto rg = ReportGuard(true);
    auto p  = rectangular_band(9, 2);
    presentation::remove_duplicate_rules(p);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);

    REQUIRE(MinimalRepOrc()
                .presentation(p)
                .target_size(18)
                .number_of_threads(std::thread::hardware_concurrency())
                .word_graph()
                .number_of_nodes()
            == 0);
    p.contains_empty_word(true);
    auto mro
        = MinimalRepOrc().presentation(p).target_size(19).number_of_threads(
            std::thread::hardware_concurrency());
    auto d = mro.word_graph();
    REQUIRE(d.number_of_nodes() == 11);
    REQUIRE(word_graph::is_strictly_cyclic(d));
    auto S = to_froidure_pin<Transf<0, node_type>>(d);
    S.add_generator(one(S.generator(0)));
    REQUIRE(S.size() == 19);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "038",
                          "partition_monoid(3) - minimal o.r.c. rep",
                          "[extreme][sims1]") {
    auto rg = ReportGuard(true);
    auto p  = partition_monoid(3, author::Machine);
    REQUIRE(!p.contains_empty_word());
    REQUIRE(p.alphabet() == 01234_w);

    auto d = RepOrc()
                 .presentation(p)
                 .target_size(203)
                 .min_nodes(1)
                 .max_nodes(22)
                 .number_of_threads(std::thread::hardware_concurrency())
                 .word_graph();
    REQUIRE(d.number_of_nodes() == 22);

    auto mro
        = MinimalRepOrc().presentation(p).target_size(203).number_of_threads(4);
    d = mro.word_graph();

    REQUIRE(word_graph::is_strictly_cyclic(d));
    auto S = to_froidure_pin<Transf<0, node_type>>(d);
    REQUIRE(S.size() == 203);
    // The actual digraph obtained is non-deterministic because we just take
    // whichever one is found first, in multiple threads
    // REQUIRE(
    //     d
    //     == make<WordGraph<node_type>>(
    //         22,
    //         {{0, 1, 0, 2, 0},      {1, 3, 3, 4, 5},      {2, 6, 6, 2, 0},
    //          {3, 0, 1, 2, 5},      {4, 4, 4, 4, 4},      {5, 5, 5, 7, 5},
    //          {6, 8, 2, 4, 0},      {7, 9, 9, 7, 5},      {8, 2, 8, 4, 10},
    //          {9, 10, 7, 11, 5},    {10, 7, 10, 12, 10},  {11, 13, 11, 11,
    //          14}, {12, 11, 13, 12, 10}, {13, 12, 12, 15, 10}, {14, 16, 14,
    //          11, 14}, {15, 15, 15, 15, 17}, {16, 18, 18, 19, 5},  {17, 19,
    //          17, 15, 17}, {18, 14, 16, 12, 5},  {19, 20, 20, 19, 21}, {20,
    //          17, 19, 15, 21}, {21, 21, 21, 19, 21}}));
    REQUIRE(d.number_of_nodes() == 22);

    std::vector<WordGraph<uint32_t>> all;

    auto hook = [&](WordGraph<uint32_t> const& x) {
      auto first = 1;
      auto SS    = to_froidure_pin<Transf<0, node_type>>(
          x, first, x.number_of_active_nodes());
      SuppressReportFor supp("FroidurePin");

      if (SS.size() == 203) {
        return all.push_back(x);
      }
    };

    auto SS = Sims1(p);

    SS.for_each(22, hook);
    REQUIRE(all.size() == 24);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "039",
                          "temperley_lieb_monoid(n) - n = 3 .. 6, minimal rep",
                          "[standard][sims1]") {
    auto rg = ReportGuard(false);

    std::array<uint64_t, 11> const sizes
        = {0, 1, 2, 5, 14, 42, 132, 429, 1'430, 4'862, 16'796};
    std::array<uint64_t, 11> const min_degrees
        = {0, 0, 2, 4, 7, 10, 20, 29, 63, 91, 208};
    // The values 63 and 91 are not verified

    for (size_t n = 3; n <= 6; ++n) {
      auto p = temperley_lieb_monoid(n);
      // There are no relations containing the empty word so we just manually
      // add it.
      p.contains_empty_word(true);
      auto orc
          = MinimalRepOrc().presentation(p).number_of_threads(2).target_size(
              sizes[n]);

      auto d = orc.word_graph();
      REQUIRE(orc.target_size() == sizes[n]);
      REQUIRE(word_graph::is_strictly_cyclic(d));
      auto S = to_froidure_pin<Transf<0, node_type>>(d);
      S.add_generator(one(S.generator(0)));
      REQUIRE(S.size() == sizes[n]);
      REQUIRE(d.number_of_nodes() == min_degrees[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "040",
                          "TransitiveGroup(10, 32) - minimal rep",
                          "[quick][sims1][no-valgrind]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.contains_empty_word(true);
    p.alphabet(01234_w);
    presentation::add_rule(p, 00_w, {});
    presentation::add_rule(p, 11_w, {});
    presentation::add_rule(p, 22_w, {});
    presentation::add_rule(p, 33_w, {});
    presentation::add_rule(p, 44_w, {});
    presentation::add_rule(p, 010101_w, {});
    presentation::add_rule(p, 0202_w, {});
    presentation::add_rule(p, 0303_w, {});
    presentation::add_rule(p, 0404_w, {});
    presentation::add_rule(p, 121212_w, {});
    presentation::add_rule(p, 1313_w, {});
    presentation::add_rule(p, 1414_w, {});
    presentation::add_rule(p, 232323_w, {});
    presentation::add_rule(p, 2424_w, {});
    presentation::add_rule(p, 343434_w, {});
    REQUIRE(MinimalRepOrc()
                .presentation(p)
                .target_size(0)
                .word_graph()
                .number_of_nodes()
            == 0);

    REQUIRE(RepOrc()
                .presentation(p)
                .min_nodes(0)
                .max_nodes(0)
                .target_size(0)
                .word_graph()
                .number_of_nodes()
            == 0);

    auto d = MinimalRepOrc().presentation(p).target_size(720).word_graph();
    REQUIRE(d.number_of_nodes() == 6);
    REQUIRE(word_graph::is_strictly_cyclic(d));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "041",
                          "rectangular_band(4, 4) - minimal o.r.c. rep",
                          "[standard][sims1]") {
    auto rg = ReportGuard(false);
    auto p  = rectangular_band(4, 4);
    p.contains_empty_word(true);
    auto d = MinimalRepOrc()
                 .presentation(p)
                 .number_of_threads(2)
                 .target_size(17)
                 .word_graph();
    REQUIRE(word_graph::is_strictly_cyclic(d));
    auto S = to_froidure_pin<Transf<0, node_type>>(d);
    REQUIRE(S.size() == 16);
    REQUIRE(d.number_of_nodes() == 7);

    p.contains_empty_word(false);
    d = MinimalRepOrc()
            .presentation(p)
            .target_size(16)
            .number_of_threads(2)
            .word_graph();
    REQUIRE(d.number_of_nodes() == 0);
  }

  LIBSEMIGROUPS_TEST_CASE("MinimalRepOrc",
                          "042",
                          "rectangular_band(m, n) - m = 1 .. 5, n = 1 .. 5",
                          "[fail][sims1]") {
    // This doesn't fail it's just very extreme
    std::vector<std::array<size_t, 6>> results = {{0, 0, 0, 0, 0, 0},
                                                  {0, 2, 2, 3, 4, 5},
                                                  {0, 3, 4, 5, 5, 6},
                                                  {0, 4, 5, 6, 6, 7},
                                                  {0, 5, 6, 7, 7, 8},
                                                  {0, 6, 7, 8, 8, 9}};

    auto rg = ReportGuard(true);
    for (size_t m = 1; m <= 5; ++m) {
      for (size_t n = 1; n <= 5; ++n) {
        std::cout << std::string(72, '#') << "\n"
                  << "CASE m, n = " << m << ", " << n << "\n"
                  << std::string(72, '#') << std::endl;

        auto p = rectangular_band(m, n);
        p.contains_empty_word(true);
        auto d = MinimalRepOrc()
                     .presentation(p)
                     .target_size(m * n + 1)
                     .number_of_threads(std::thread::hardware_concurrency())
                     .word_graph();
        REQUIRE(word_graph::is_strictly_cyclic(d));
        auto S = to_froidure_pin<Transf<0, node_type>>(d);
        REQUIRE(S.size() == m * n);
        REQUIRE(d.number_of_nodes() == results[m][n]);
      }
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "043",
                          "rectangular_band(2, 2) - with and without identity",
                          "[quick][sims1]") {
    auto rg = ReportGuard(false);
    auto p  = rectangular_band(2, 2);
    REQUIRE(!p.contains_empty_word());
    Sims1 S;
    S.presentation(p);

    REQUIRE(S.number_of_congruences(4) == 6);

    p.contains_empty_word(true);

    Sims1 T;
    T.presentation(p);
    REQUIRE(T.number_of_congruences(5) == 9);
    // sims::dot_poset("example-043-rectangular-band^1", T.cbegin(5),
    // T.cend(5));

    auto it = S.cbegin(4);

    REQUIRE(*it++
            == make<WordGraph<node_type>>(
                5, {{1, 1, 1, 1}, {1, 1, 1, 1}}));  // Good
    REQUIRE(*it++
            == make<WordGraph<node_type>>(
                5, {{1, 1, 1, 2}, {1, 1, 1, 2}, {1, 1, 1, 2}}));  // Good
    REQUIRE(*it++
            == make<WordGraph<node_type>>(
                5, {{1, 2, 1, 1}, {1, 1, 1, 1}, {2, 2, 2, 2}}));  // Good
    REQUIRE(
        *it++
        == make<WordGraph<node_type>>(
            5,
            {{1, 2, 1, 1}, {1, 1, 1, 1}, {2, 2, 2, 3}, {2, 2, 2, 3}}));  // Good
    REQUIRE(
        *it++
        == make<WordGraph<node_type>>(
            5,
            {{1, 2, 1, 3}, {1, 1, 1, 3}, {2, 2, 2, 2}, {1, 1, 1, 3}}));  // Good
    REQUIRE(*it++
            == make<WordGraph<node_type>>(5,
                                          {{1, 2, 1, 3},
                                           {1, 1, 1, 3},
                                           {2, 2, 2, 4},
                                           {1, 1, 1, 3},
                                           {2, 2, 2, 4}}));  // Good
    REQUIRE(it->number_of_nodes() == 0);

    it = T.cbegin(5);

    REQUIRE(*it++ == make<WordGraph<node_type>>(5, {{0, 0, 0, 0}}));
    REQUIRE(*it++
            == make<WordGraph<node_type>>(5, {{0, 0, 0, 1}, {0, 0, 0, 1}}));
    REQUIRE(*it++
            == make<WordGraph<node_type>>(5, {{1, 1, 1, 0}, {1, 1, 1, 0}}));
    REQUIRE(*it++
            == make<WordGraph<node_type>>(5, {{1, 1, 1, 1}, {1, 1, 1, 1}}));
    REQUIRE(*it++
            == make<WordGraph<node_type>>(
                5, {{1, 1, 1, 2}, {1, 1, 1, 2}, {1, 1, 1, 2}}));
    REQUIRE(*it++
            == make<WordGraph<node_type>>(
                5, {{1, 2, 1, 1}, {1, 1, 1, 1}, {2, 2, 2, 2}}));
    REQUIRE(*it++
            == make<WordGraph<node_type>>(
                5, {{1, 2, 1, 1}, {1, 1, 1, 1}, {2, 2, 2, 3}, {2, 2, 2, 3}}));
    REQUIRE(*it++
            == make<WordGraph<node_type>>(
                5, {{1, 2, 1, 3}, {1, 1, 1, 3}, {2, 2, 2, 2}, {1, 1, 1, 3}}));
    REQUIRE(*it++
            == make<WordGraph<node_type>>(5,
                                          {{1, 2, 1, 3},
                                           {1, 1, 1, 3},
                                           {2, 2, 2, 4},
                                           {1, 1, 1, 3},
                                           {2, 2, 2, 4}}));
    REQUIRE(it->number_of_nodes() == 0);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "044",
                          "trivial group - minimal o.r.c. rep",
                          "[quick][sims1]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("aAbB");
    p.contains_empty_word(true);
    presentation::add_inverse_rules(p, "AaBb");
    presentation::add_rule(p, "ab", "");
    presentation::add_rule(p, "abb", "");

    Sims1 S;
    S.presentation(p);

    REQUIRE(S.number_of_congruences(10) == 1);
    auto d = MinimalRepOrc().presentation(p).target_size(1).word_graph();
    REQUIRE(d.number_of_nodes() == 1);
    REQUIRE(word_graph::is_strictly_cyclic(d));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "045",
                          "onesided zero semigroup - minimal o.r.c. rep",
                          "[quick][sims1]") {
    // This is an example of a semigroup with a strictly cyclic faithful
    // onesided representation.
    auto         rg = ReportGuard(false);
    size_t const n  = 5;
    auto         p  = rectangular_band(1, n);
    auto d = MinimalRepOrc().presentation(p).target_size(n).word_graph();
    REQUIRE(word_graph::is_strictly_cyclic(d));
    auto S = to_froidure_pin<Transf<0, node_type>>(d);
    REQUIRE(S.size() == n);
    REQUIRE(d.number_of_nodes() == 5);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "046",
                          "semigroup with faithful non-strictly cyclic action",
                          "[quick][sims1]") {
    // Found with Smallsemi, this example is minimal wrt size of the
    // semigroup.

    auto rg = ReportGuard(false);

    auto S = to_froidure_pin({make<Transf<6>>({0, 0, 2, 1, 4, 1}),
                              make<Transf<6>>({0, 0, 2, 3, 4, 3}),
                              make<Transf<6>>({0, 2, 2, 0, 4, 4})});

    REQUIRE(S.size() == 5);
    auto p = to_presentation<word_type>(S);
    auto d = MinimalRepOrc().presentation(p).target_size(5).word_graph();
    REQUIRE(word_graph::is_strictly_cyclic(d));
    REQUIRE(d.number_of_nodes() == 4);
    REQUIRE(d
            == make<WordGraph<uint32_t>>(
                4, {{2, 2, 3}, {0, 1, 2}, {2, 2, 2}, {3, 3, 3}}));
    auto T = to_froidure_pin<Transf<4>>(d);
    REQUIRE(T.generator(0) == Transf<4>({2, 0, 2, 3}));
    REQUIRE(T.generator(1) == Transf<4>({2, 1, 2, 3}));
    REQUIRE(T.generator(2) == Transf<4>({3, 2, 2, 3}));
    REQUIRE(T.size() == 5);

    auto dd = make<WordGraph<uint8_t>>(5,
                                       {{0, 0, 0, 0, 0},
                                        {0, 0, 0, 0, 2},
                                        {2, 2, 2, 2, 2},
                                        {0, 1, 2, 3, 0},
                                        {4, 4, 4, 4, 4}});

    REQUIRE(!word_graph::is_strictly_cyclic(dd));
    REQUIRE(dd.number_of_nodes() == 5);
    auto U = to_froidure_pin<Transf<5>>(dd);
    REQUIRE(U.size() == 5);

    Sims1 C;
    C.presentation(p);
    REQUIRE(C.number_of_congruences(5) == 9);
    uint64_t strictly_cyclic_count     = 0;
    uint64_t non_strictly_cyclic_count = 0;

    for (auto it = C.cbegin(5); it != C.cend(5); ++it) {
      auto W = to_froidure_pin<Transf<0, node_type>>(
          *it, 1, it->number_of_active_nodes());
      if (p.contains_empty_word()) {
        auto id = one(W.generator(0));
        if (!W.contains(id)) {
          W.add_generator(id);
        }
      }
      if (W.size() == 5) {
        auto result = *it;
        result.induced_subgraph_no_checks(1, result.number_of_active_nodes());
        result.number_of_active_nodes(result.number_of_active_nodes() - 1);
        if (word_graph::is_strictly_cyclic(result)) {
          strictly_cyclic_count++;
        } else {
          REQUIRE(W.generator(0) == Transf<0, node_type>({3, 0, 2, 3, 4}));
          REQUIRE(W.generator(1) == Transf<0, node_type>({3, 1, 2, 3, 4}));
          REQUIRE(W.generator(2) == Transf<0, node_type>({4, 3, 2, 3, 4}));
          REQUIRE(
              result
              == make<WordGraph<uint32_t>>(
                  5, {{3, 3, 4}, {0, 1, 3}, {2, 2, 2}, {3, 3, 3}, {4, 4, 4}}));
          non_strictly_cyclic_count++;
        }
      }
    }
    REQUIRE(strictly_cyclic_count == 2);
    REQUIRE(non_strictly_cyclic_count == 1);
  }

  // Takes about 3 to 4 minutes
  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "047",
                          "rectangular_band(m, n) - m = 1 .. 5, n = 1 .. 5",
                          "[fail][sims1]") {
    // This doesn't fail it's just very extreme
    // Note: num_congs[n][m] is the number of right congruences of
    // rectangular_band(m, n)
    std::vector<std::array<size_t, 6>> num_congs
        = {{0, 0, 0, 0, 0, 0},
           {0, 0, 0, 0, 0, 0},
           {0, 0, 6, 22, 94, 454},
           {0, 0, 30, 205, 1'555, 12'880},
           {0, 0, 240, 4'065, 72'465, 1'353'390},
           {0, 0, 2'756, 148'772, 8'174'244, 456'876'004}};

    // Seems like the m,n-th entry of the table above is:
    // {m, n} ->  Sum([0 .. n], k -> Bell(m)^k*Stirling2(n, k));

    auto rg = ReportGuard(true);
    for (size_t m = 2; m <= 5; ++m) {
      for (size_t n = 2; n <= 5; ++n) {
        std::cout << std::string(72, '#') << "\n"
                  << "CASE m, n = " << m << ", " << n << "\n"
                  << std::string(72, '#') << std::endl;

        auto  p = rectangular_band(m, n);
        Sims1 S(p);
        REQUIRE(S.number_of_threads(std::thread::hardware_concurrency())
                    .number_of_congruences(m * n)
                == num_congs[n][m]);
        presentation::reverse(p);
        S.init(p);
        REQUIRE(S.number_of_threads(std::thread::hardware_concurrency())
                    .number_of_congruences(m * n)
                == num_congs[m][n]);
      }
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "048",
                          "stellar_monoid(n) n = 3",
                          "[quick][sims1][babbage]") {
    auto   rg = ReportGuard(true);
    size_t n  = 3;
    auto   p  = zero_rook_monoid(n);
    auto   q  = stellar_monoid(n);
    p.rules.insert(p.rules.end(), q.rules.cbegin(), q.rules.cend());
    p.validate();
    REQUIRE(p.alphabet().size() == n);
    presentation::reverse(p);
    Sims1 S;
    S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
    REQUIRE(S.number_of_congruences(16) == 1'550);
    presentation::reverse(p);
    S.presentation(p);
    REQUIRE(S.number_of_congruences(16) == 1'521);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "049",
                          "stylic_monoid(n) n = 3, 4",
                          "[fail][sims1]") {
    // This doesn't fail it's just very extreme
    auto                           rg   = ReportGuard(true);
    std::array<uint64_t, 10> const size = {0, 0, 0, 14, 51};
    //               1505s
    std::array<uint64_t, 10> const num_left  = {0, 0, 0, 1'318, 1'431'795'099};
    std::array<uint64_t, 10> const num_right = {0, 0, 0, 1'318, 1'431'795'099};

    auto p = stylic_monoid(4);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    auto q = p;
    presentation::sort_each_rule(q);
    presentation::sort_rules(q);
    REQUIRE(p == q);

    for (size_t n = 3; n < 5; ++n) {
      auto  p = stylic_monoid(n);
      Sims1 S;
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(size[n]) == num_right[n]);
      presentation::reverse(p);
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(size[n]) == num_left[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "050",
                          "(2, 3, 7)-triangle group - index 50",
                          "[extreme][sims1]") {
    auto                      rg = ReportGuard(true);
    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("xy");
    presentation::add_rule(p, "xx", "");
    presentation::add_rule(p, "yyy", "");
    presentation::add_rule(p, "xyxyxyxy", "yyxyyxyyx");
    Sims1 S;
    S.presentation(p).number_of_threads(1);
    REQUIRE(S.number_of_congruences(50) == 75'971);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "051",
                          "Heineken group - index 10",
                          "[extreme][sims1]") {
    auto                      rg = ReportGuard(true);
    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("xXyY");
    presentation::add_inverse_rules(p, "XxYy");
    presentation::add_rule(p, "yXYYxyYYxyyXYYxyyXyXYYxy", "x");
    presentation::add_rule(p, "YxyyXXYYxyxYxyyXYXyXYYxxyyXYXyXYYxyx", "y");
    presentation::balance_no_checks(p, p.alphabet(), std::string("XxYy"));

    REQUIRE(p.rules
            == std::vector<std::string>({"xX",
                                         "",
                                         "Xx",
                                         "",
                                         "yY",
                                         "",
                                         "Yy",
                                         "",
                                         "yXYYxyYYxyyXY",
                                         "xYXyyxYxYYXy",
                                         "YxyyXXYYxyxYxyyXYXy",
                                         "yXYXyyxYxyxYYXXyyx"}));

    Sims1 S;
    S.presentation(p)
        .number_of_threads(std::thread::hardware_concurrency())
        .long_rule_length(37);
    REQUIRE(S.number_of_congruences(10) == 1);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "052",
                          "temperley_lieb_monoid(n) - n = 3 .. 6",
                          "[extreme][low-index][babbage]") {
    std::array<uint64_t, 10> const size = {0, 0, 0, 5, 14, 42, 132, 429};
    std::array<uint64_t, 10> const num_right
        = {0, 0, 0, 9, 79, 2'157, 4'326'459};

    auto rg = ReportGuard(true);
    for (size_t n = 3; n < 7; ++n) {
      auto p = temperley_lieb_monoid(n);
      p.contains_empty_word(true);
      Sims1 S;
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(size[n]) == num_right[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "053",
                          "partial_transformation_monoid(3)",
                          "[extreme][low-index]") {
    auto  rg = ReportGuard(true);
    auto  p  = partial_transformation_monoid(3, author::Machine);
    Sims1 S;
    S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
    REQUIRE(S.number_of_congruences(64) == 92'703);
    presentation::reverse(p);
    S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
    REQUIRE(S.number_of_congruences(64) == 371);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "054",
                          "partial_transformation_monoid(4) from FroidurePin",
                          "[fail][low-index]") {
    using Transf_ = Transf<5>;
    auto rg       = ReportGuard(true);

    auto S = to_froidure_pin({Transf_({1, 0, 2, 3, 4}),
                              Transf_({3, 0, 1, 2, 4}),
                              Transf_({4, 1, 2, 3, 4}),
                              Transf_({1, 1, 2, 3, 4})});
    REQUIRE(S.size() == 625);
    auto p = to_presentation<word_type>(S);
    presentation::reverse(p);
    Sims1 C(p);
    REQUIRE(presentation::longest_rule_length(p) == 18);
    REQUIRE(presentation::shortest_rule_length(p) == 3);

    presentation::remove_duplicate_rules(p);
    presentation::remove_trivial_rules(p);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);

    word_type w = presentation::longest_subword_reducing_length(p);
    while (!w.empty()) {
      presentation::replace_word_with_new_generator(p, w);
      w = presentation::longest_subword_reducing_length(p);
    }

    REQUIRE(presentation::length(p) == 1419);
    REQUIRE(presentation::longest_rule_length(p) == 6);

    C.presentation(p).long_rule_length(6).number_of_threads(
        std::thread::hardware_concurrency());
    // NOTE: Never ran to completion, there should be a non-zero number of
    // congruences.
    REQUIRE(C.number_of_congruences(625) == 0);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "055",
                          "Plactic semigroup 3 up to index 8",
                          "[extreme][low-index][plactic]") {
    std::array<uint64_t, 9> const num
        = {0, 1, 29, 484, 6'896, 103'204, 1'773'360, 35'874'182, 849'953'461};
    auto rg = ReportGuard(true);
    auto p  = plactic_monoid(3);
    p.contains_empty_word(false);
    Sims1 S;
    for (size_t n = 2; n < 9; ++n) {
      S.init(p);
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
      presentation::reverse(p);
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "056",
                          "Plactic semigroup 4 up to index 6",
                          "[extreme][low-index][plactic]") {
    std::array<uint64_t, 8> const num
        = {0, 1, 67, 2'794, 106'264, 4'795'980, 278'253'841, 20'855'970'290};
    // Last value took 1h34m to compute so is not included.
    auto rg = ReportGuard(true);
    auto p  = plactic_monoid(4);
    p.contains_empty_word(false);
    for (size_t n = 2; n < 7; ++n) {
      Sims1 S;
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
      presentation::reverse(p);
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "057",
                          "Plactic semigroup 5 up to index 5",
                          "[extreme][low-index][plactic]") {
    std::array<uint64_t, 7> const num
        = {0, 1, 145, 14'851, 1'496'113, 198'996'912, 37'585'675'984};
    // Last value took 5h11m to compute
    auto rg = ReportGuard(true);
    auto p  = plactic_monoid(5);
    p.contains_empty_word(false);
    for (size_t n = 3; n < 6; ++n) {
      Sims1 S;
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
      presentation::reverse(p);
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "058",
                          "Plactic semigroup 6 up to index 4",
                          "[extreme][low-index][plactic]") {
    std::array<uint64_t, 6> const num
        = {0, 1, 303, 77'409, 20'526'128, 7'778'840'717};
    // The last value took 4h5m to run and is omitted.
    auto rg = ReportGuard(true);
    auto p  = plactic_monoid(6);
    p.contains_empty_word(false);
    for (size_t n = 2; n < 5; ++n) {
      Sims1 S;
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      presentation::reverse(p);
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "059",
                          "Plactic semigroup 7 up to index 3",
                          "[extreme][low-index][plactic]") {
    std::array<uint64_t, 5> const num = {0, 1, 621, 408'024, 281'600'130};
    // The last value took approx. 12m34s to run and is omitted from the
    // extreme test.
    auto rg = ReportGuard(true);
    auto p  = plactic_monoid(7);
    p.contains_empty_word(false);
    for (size_t n = 2; n < 4; ++n) {
      Sims1 S;
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
      presentation::reverse(p);
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "060",
                          "Plactic semigroup 8 up to index 3",
                          "[extreme][low-index][plactic]") {
    std::array<uint64_t, 4> const num = {0, 1, 1'259, 2'201'564};
    auto                          rg  = ReportGuard(true);
    auto                          p   = plactic_monoid(8);
    p.contains_empty_word(false);
    for (size_t n = 2; n < 4; ++n) {
      Sims1 S;
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
      presentation::reverse(p);
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "061",
                          "Chinese semigroup 3 up to index 8",
                          "[extreme][low-index][chinese]") {
    std::array<uint64_t, 9> const num
        = {0, 1, 31, 559, 8'904, 149'529, 2'860'018, 63'828'938, 1'654'488'307};
    // index 8 is doable and the value is included above, but it took about X
    // minutes, where X could be considered large, so isn't included in the
    // loop below.
    auto rg = ReportGuard(true);
    auto p  = chinese_monoid(3);
    p.contains_empty_word(false);
    for (size_t n = 2; n < 8; ++n) {
      Sims1 S;
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "062",
                          "Chinese semigroup 4 up to index 6",
                          "[extreme][low-index][chinese]") {
    std::array<uint64_t, 8> const num
        = {0, 1, 79, 3'809, 183'995, 10'759'706, 804'802'045, 77'489'765'654};
    // n = 6 took between 3 and 4 minutes
    // n = 7 took 6h16m
    // 7 is omitted
    auto rg = ReportGuard(true);
    auto p  = chinese_monoid(4);
    p.contains_empty_word(false);
    for (size_t n = 3; n < 7; ++n) {
      Sims1 S;
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "063",
                          "Chinese semigroup 5 up to index 5",
                          "[extreme][low-index][chinese]") {
    std::array<uint64_t, 7> const num
        = {0, 1, 191, 23'504, 3'382'921, 685'523'226, 199'011'439'587};
    // The last value took 21h32m and so is omitted
    auto rg = ReportGuard(true);
    auto p  = chinese_monoid(5);
    p.contains_empty_word(false);
    for (size_t n = 3; n < 6; ++n) {
      Sims1 S;
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "064",
                          "Chinese semigroup 6 up to index 4",
                          "[extreme][low-index][chinese]") {
    std::array<uint64_t, 6> const num
        = {0, 1, 447, 137'694, 58'624'384, 40'823'448'867};
    // The last value took 9h54m to compute, and is omitted!
    auto rg = ReportGuard(true);
    auto p  = chinese_monoid(6);
    p.contains_empty_word(false);
    for (size_t n = 3; n < 5; ++n) {
      Sims1 S;
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "065",
                          "Chinese semigroup 7 up to index 4",
                          "[extreme][low-index][chinese]") {
    std::array<uint64_t, 5> const num = {0, 1, 1'023, 786'949, 988'827'143};
    // Last value took about 50m to compute
    auto rg = ReportGuard(true);
    auto p  = chinese_monoid(7);
    p.contains_empty_word(false);
    for (size_t n = 2; n < 4; ++n) {
      Sims1 S;
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "066",
                          "Chinese semigroup 8 up to index 3",
                          "[extreme][low-index][chinese]") {
    std::array<uint64_t, 4> const num = {0, 1, 2'303, 4'459'599};
    auto                          rg  = ReportGuard(true);
    auto                          p   = chinese_monoid(8);
    p.contains_empty_word(false);
    for (size_t n = 2; n < 4; ++n) {
      Sims1 S;
      S.presentation(p).number_of_threads(std::thread::hardware_concurrency());
      REQUIRE(S.number_of_congruences(n) == num[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "067",
                          "FreeSemigroup(n) up to index 3",
                          "[extreme][low-index]") {
    // (27^n - 9^n)/2 - 12^n + 6^n
    std::array<uint64_t, 10> const num
        = {0, 2, 229, 8'022, 243'241, 6'904'866, 190'509'229, 5'192'249'502};
    // 4 ^ n - 2 ^ n
    //  {0, 2, 13, 57, 241, 993, 4033};
    //    = {0, 1, 29, 249, 2'033, 16'353, 131'009};
    // = {0, 1, 830, 81'762, 7'008'614};
    auto rg = ReportGuard(true);
    for (size_t n = 2; n < 8; ++n) {
      Presentation<word_type> p;
      p.contains_empty_word(true);
      p.alphabet(n);
      Sims1 S;
      S.presentation(p);
      REQUIRE(S.number_of_threads(std::thread::hardware_concurrency())
                  .number_of_congruences(3)
              == num[n]);
    }

    // For n >= 1, a(n) is the number of deterministic, completely-defined,
    // initially-connected finite automata with n inputs and 3 unlabeled
    // states. A020522 counts similar automata with n inputs and 2 unlabeled
    // states.

    // According to a comment by Nelma Moreira in A006689 and A006690, the
    // number of such automata with N inputs and M unlabeled states is Sum
    // (Product_{i=1..M-1} i^(f_i - f_{i-1} - 1)) * M^(M*N - f_{M-1} - 1),
    // where the sum is taken over integers f_1, ..., f_{M-1} satisfying 0 <=
    // f_1 < N and f_{i-1} < f_{i} < i*N for i = 2..M-1. (See Theorem 8 in
    // Almeida, Moreira, and Reis (2007). The value of f_0 is not relevant.)
    // For this sequence we have M = 3 unlabeled states, for A020522 we have M
    // = 2 unlabeled states, for A006689 we have N = 2 inputs, and for A006690
    // we have N = 3 inputs.
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "068",
                          "RepOrc",
                          "[quick][low-index][no-valgrind]") {
    auto rg = ReportGuard(true);

    auto   p = temperley_lieb_monoid(9);
    RepOrc orc;
    // Check bad input
    auto d = orc.presentation(p)
                 .min_nodes(100)
                 .max_nodes(90)
                 .target_size(4'862)
                 .word_graph();
    REQUIRE(d.number_of_nodes() == 0);

    d = orc.presentation(p)
            .min_nodes(80)
            .max_nodes(100)
            .target_size(4'862)
            .word_graph();

    auto S = to_froidure_pin<Transf<0, node_type>>(d);
    S.add_generator(one(S.generator(0)));
    REQUIRE(S.size() == 4'862);
    REQUIRE(orc.min_nodes() == 80);
    REQUIRE(orc.max_nodes() == 100);
    REQUIRE(orc.target_size() == 4'862);
    REQUIRE(orc.presentation().rules.size() == 128);
    REQUIRE(orc.number_of_long_rules() == 0);
    REQUIRE(d.number_of_nodes() == 91);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "069",
                          "fp example 1 (settings)",
                          "[quick][low-index]") {
    auto rg = ReportGuard(true);

    Presentation<word_type> p;
    p.contains_empty_word(true);
    p.alphabet(01_w);
    presentation::add_rule(p, 000_w, 0_w);
    presentation::add_rule(p, 11_w, 1_w);
    presentation::add_rule(p, 0101_w, 0_w);

    Sims1 S;
    REQUIRE(S.presentation(p)
                .cbegin_long_rules(4)
                .number_of_threads(1)
                .number_of_congruences(5)
            == 6);
    S.long_rule_length(5);
    REQUIRE(S.number_of_long_rules() == 1);
    REQUIRE(S.presentation().rules.size() == 6);
    S.long_rule_length(4);
    REQUIRE(S.number_of_long_rules() == 2);
    REQUIRE(S.presentation().rules.size() == 6);
    REQUIRE(S.settings().number_of_long_rules() == 2);
    REQUIRE(S.settings().presentation().rules.size() == 6);

    presentation::reverse(p);
    REQUIRE(S.presentation(p)
                .cbegin_long_rules(4)
                .number_of_threads(4)
                .number_of_congruences(5)
            == 9);
    S.clear_long_rules();
    REQUIRE(S.presentation(p)
                .cbegin_long_rules(4)
                .number_of_threads(4)
                .number_of_congruences(5)
            == 9);

    Sims2 S2;
    REQUIRE(S2.presentation(p)
                .cbegin_long_rules(4)
                .number_of_threads(4)
                .number_of_congruences(5)
            == 6);
    S2.long_rule_length(5);
    REQUIRE(S2.number_of_long_rules() == 1);
    REQUIRE(S2.presentation().rules.size() == 6);
    S2.long_rule_length(4);
    REQUIRE(S2.number_of_long_rules() == 2);
    REQUIRE(S2.presentation().rules.size() == 6);
    REQUIRE(S2.settings().number_of_long_rules() == 2);
    REQUIRE(S2.settings().presentation().rules.size() == 6);

    presentation::reverse(p);
    REQUIRE(S2.presentation(p)
                .cbegin_long_rules(4)
                .number_of_threads(4)
                .number_of_congruences(5)
            == 6);
    S2.clear_long_rules();
    REQUIRE(S2.presentation(p)
                .cbegin_long_rules(4)
                .number_of_threads(4)
                .number_of_congruences(5)
            == 6);

    Presentation<word_type> q;

    q.alphabet({0, 1});
    q.contains_empty_word(true);
    presentation::add_rule(q, 000_w, 0_w);
    presentation::add_rule(q, 111_w, ""_w);
    presentation::add_rule(q, 011_w, 10_w);

    std::vector<word_type> forbid = {0_w, 01_w, 00_w, ""_w};
    SimsRefinerFaithful    pruno(forbid);

    RepOrc Ro;
    REQUIRE(Ro.presentation(q)
                .target_size(9)
                .min_nodes(2)
                .max_nodes(6)
                .cbegin_long_rules(4)
                .number_of_threads(4)
                .add_pruner(pruno)
                .word_graph()
                .number_of_active_nodes()
            == 6);
    Ro.long_rule_length(5);
    REQUIRE(Ro.number_of_long_rules() == 1);
    REQUIRE(Ro.presentation().rules.size() == 6);
    Ro.long_rule_length(4);
    REQUIRE(Ro.number_of_long_rules() == 2);
    REQUIRE(Ro.presentation().rules.size() == 6);
    REQUIRE(Ro.settings().number_of_long_rules() == 2);
    REQUIRE(Ro.settings().presentation().rules.size() == 6);
    Ro.clear_long_rules();
    Ro.clear_pruners();
    REQUIRE(Ro.presentation(q)
                .target_size(9)
                .min_nodes(2)
                .max_nodes(6)
                .cbegin_long_rules(4)
                .number_of_threads(4)
                .word_graph()
                .number_of_active_nodes()
            == 6);
    Ro.init();
    REQUIRE(Ro.presentation(q)
                .target_size(9)
                .min_nodes(2)
                .max_nodes(6)
                .cbegin_long_rules(4)
                .number_of_threads(4)
                .add_pruner(pruno)
                .word_graph()
                .number_of_active_nodes()
            == 6);

    MinimalRepOrc Mro;
    REQUIRE(Mro.presentation(q)
                .target_size(9)
                .cbegin_long_rules(4)
                .number_of_threads(4)
                .add_pruner(pruno)
                .word_graph()
                .number_of_active_nodes()
            == 6);
    Mro.long_rule_length(5);
    REQUIRE(Mro.number_of_long_rules() == 1);
    REQUIRE(Mro.presentation().rules.size() == 6);
    Mro.long_rule_length(4);
    REQUIRE(Mro.number_of_long_rules() == 2);
    REQUIRE(Mro.presentation().rules.size() == 6);
    REQUIRE(Mro.settings().number_of_long_rules() == 2);
    REQUIRE(Mro.settings().presentation().rules.size() == 6);
    Mro.clear_long_rules();
    Mro.clear_pruners();
    REQUIRE(Mro.presentation(q)
                .target_size(9)
                .cbegin_long_rules(4)
                .number_of_threads(4)
                .word_graph()
                .number_of_active_nodes()
            == 6);
    Mro.init();
    REQUIRE(Mro.presentation(q)
                .target_size(9)
                .cbegin_long_rules(4)
                .number_of_threads(4)
                .add_pruner(pruno)
                .word_graph()
                .number_of_active_nodes()
            == 6);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "070",
                          "temperley_lieb_monoid(3) - n = minimal rep "
                          "(single-threaded, reporting on)",
                          "[standard][sims1]") {
    std::cout << "\n";  // So that the reporting looks good
    auto rg = ReportGuard(true);

    for (size_t n = 3; n <= 3; ++n) {
      auto p = temperley_lieb_monoid(n);
      // There are no relations containing the empty word so we just manually
      // add it.
      p.contains_empty_word(true);
      auto d = MinimalRepOrc()
                   .presentation(p)
                   .number_of_threads(1)
                   .target_size(5)
                   .word_graph();
      REQUIRE(word_graph::is_strictly_cyclic(d));
      auto S = to_froidure_pin<Transf<0, node_type>>(d);
      S.add_generator(one(S.generator(0)));
      REQUIRE(S.size() == 5);
      REQUIRE(d.number_of_nodes() == 4);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "071",
                          "FreeSemigroup(2) up to index 4",
                          "[quick][low-index]") {
    auto                    rg = ReportGuard(true);
    Presentation<word_type> p;
    p.contains_empty_word(true);
    p.alphabet(2);
    Sims1 S;
    S.presentation(p);
    std::cout << "\n";  // So that the reporting looks good
    REQUIRE(S.number_of_threads(2).number_of_congruences(4) == 5'477);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "072",
                          "symmetric_group(n) for n = 4",
                          "[quick][low-index]") {
    std::array<uint64_t, 10> const num = {0, 0, 0, 6, 30, 156, 1'455};
    auto                           rg  = ReportGuard(false);
    size_t                         n   = 4;
    auto                           p   = symmetric_group(n, author::Carmichael);
    Sims1                          C;
    C.presentation(p).number_of_threads(4);
    REQUIRE(C.number_of_congruences(factorial(n)) == num[n]);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "073",
                          "corner case no generators + no relations",
                          "[quick][low-index]") {
    Presentation<word_type> p;
    p.alphabet(0);
    Sims1 S;
    REQUIRE_THROWS_AS(S.presentation(p), LibsemigroupsException);
    REQUIRE_THROWS_AS(S.number_of_congruences(1), LibsemigroupsException);
    REQUIRE_THROWS_AS(S.cbegin(2), LibsemigroupsException);
    REQUIRE_THROWS_AS(S.cend(2), LibsemigroupsException);
    REQUIRE_THROWS_AS(S.find_if(2, [](auto) { return true; }),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(S.for_each(2, [](auto) {}), LibsemigroupsException);
    REQUIRE_THROWS_AS(sims::add_included_pair(S, 01_w, 10_w),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(sims::add_excluded_pair(S, 01_w, 10_w),
                      LibsemigroupsException);

    p.alphabet(2);
    S.presentation(p);
    REQUIRE_THROWS_AS(sims::add_excluded_pair(S, 01_w, 102_w),
                      LibsemigroupsException);

    p.alphabet(3);
    S.presentation(p);
    sims::add_excluded_pair(S, 01_w, 102_w);
    p.alphabet(2);
    REQUIRE_THROWS_AS(S.presentation(p), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "074",
                          "monogenic_semigroup(m, r) for m, r = 1 .. 10",
                          "[quick][low-index]") {
    auto                                        rg = ReportGuard(false);
    std::vector<std::array<uint64_t, 11>> const num
        = {{1, 2, 2, 3, 2, 4, 2, 4, 3, 4},
           {2, 4, 4, 6, 4, 8, 4, 8, 6, 8},
           {3, 6, 6, 9, 6, 12, 6, 12, 9, 12},
           {4, 8, 8, 12, 8, 16, 8, 16, 12, 16},
           {5, 10, 10, 15, 10, 20, 10, 20, 15, 20},
           {6, 12, 12, 18, 12, 24, 12, 24, 18, 24},
           {7, 14, 14, 21, 14, 28, 14, 28, 21, 28},
           {8, 16, 16, 24, 16, 32, 16, 32, 24, 32},
           {9, 18, 18, 27, 18, 36, 18, 36, 27, 36},
           {10, 20, 20, 30, 20, 40, 20, 40, 30, 40}};

    // m * number of divisors of r

    for (size_t m = 1; m <= 10; ++m) {
      for (size_t r = 1; r <= 10; ++r) {
        // Cyclic groups
        auto p = monogenic_semigroup(m, r);

        Sims1 C;
        C.presentation(p);
        // std::cout << C.number_of_congruences(m + r) << ", ";
        REQUIRE(C.number_of_congruences(m + r) == num[m - 1][r - 1]);
      }
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "075",
                          "partial_transformation_monoid(4)",
                          "[fail][low-index]") {
    auto rg = ReportGuard(true);
    auto p  = partial_transformation_monoid(4, author::Sutov);
    auto w  = presentation::longest_subword_reducing_length(p);
    while (!w.empty()) {
      presentation::replace_word_with_new_generator(
          p, presentation::longest_subword_reducing_length(p));
      w = presentation::longest_subword_reducing_length(p);
    }

    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    presentation::remove_duplicate_rules(p);
    presentation::reduce_complements(p);
    presentation::remove_trivial_rules(p);

    do {
      auto it = knuth_bendix::redundant_rule(p, std::chrono::milliseconds(100));
      p.rules.erase(it, it + 2);
    } while (presentation::length(p) > 800);
    presentation::reverse(p);
    Sims1 C(p);
    C.presentation(p).number_of_threads(std::thread::hardware_concurrency());
    // NOTE: Never ran to completion, there should be a non-zero number of
    // congruences.
    REQUIRE(C.number_of_congruences(624) == 0);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "076",
                          "uninitialized RepOrc",
                          "[quick][low-index]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abc");
    presentation::add_rule(p, "cc", "c");
    presentation::add_rule(p, "abb", "a");
    presentation::add_rule(p, "aca", "aba");

    RepOrc orc;
    orc.presentation(p).number_of_threads(std::thread::hardware_concurrency());
    REQUIRE(orc.min_nodes() == 0);
    REQUIRE(orc.max_nodes() == 0);
    REQUIRE(orc.target_size() == 0);
    REQUIRE(orc.word_graph().number_of_nodes() == 0);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "077",
                          "2x2 simple semigroups over S(4)",
                          "[quick][sims1]") {
    auto rg = ReportGuard(true);

    Presentation<std::string> p;
    p.alphabet("abc");

    // S := ReesMatrixSemigroup(SymmetricGroup(4), [[(1, 2), ()], [(), ()]]);
    presentation::add_rule(p, "cc", "c");
    presentation::add_rule(p, "abb", "a");
    presentation::add_rule(p, "aca", "aba");
    presentation::add_rule(p, "acb", "a");
    presentation::add_rule(p, "bba", "a");
    presentation::add_rule(p, "bbb", "b");
    presentation::add_rule(p, "bca", "a");
    presentation::add_rule(p, "bcb", "b");
    presentation::add_rule(p, "cbc", "c");
    presentation::add_rule(p, "aaaa", "bb");
    presentation::add_rule(p, "baaa", "abab");
    presentation::add_rule(p, "baba", "aaab");
    presentation::add_rule(p, "abaaba", "baab");
    presentation::add_rule(p, "baabaa", "aabaab");
    presentation::add_rule(p, "baabab", "aaabaa");
    presentation::add_rule(p, "aaabaab", "baaba");
    // Minimum rep. o.r.c. 6
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    p.validate();

    MinimalRepOrc orc;
    auto          d = orc.presentation(p)
                 .target_size(96)
                 .number_of_threads(std::thread::hardware_concurrency())
                 .word_graph();
    REQUIRE(d.number_of_nodes() == 6);
    REQUIRE(orc.target_size() == 96);

    // S := ReesMatrixSemigroup(SymmetricGroup(4), [[(1, 2, 3, 4), ()], [(),
    // ()]]);
    p.init();
    p.alphabet("abc");
    presentation::add_rule(p, "cc", "c");
    presentation::add_rule(p, "abb", "a");
    presentation::add_rule(p, "aca", "a");
    presentation::add_rule(p, "acb", "b");
    presentation::add_rule(p, "bba", "a");
    presentation::add_rule(p, "bbb", "b");
    presentation::add_rule(p, "bca", "b");
    presentation::add_rule(p, "bcb", "aba");
    presentation::add_rule(p, "cac", "c");
    presentation::add_rule(p, "aaaa", "bb");
    presentation::add_rule(p, "baaa", "abab");
    presentation::add_rule(p, "baba", "aaab");
    presentation::add_rule(p, "abaaba", "baab");
    presentation::add_rule(p, "baabaa", "aabaab");
    presentation::add_rule(p, "baabab", "aaabaa");
    presentation::add_rule(p, "aaabaab", "baaba");
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    p.validate();
    // Minimum rep. o.r.c. 8
    d = orc.presentation(p)
            .target_size(96)
            .number_of_threads(std::thread::hardware_concurrency())
            .word_graph();
    REQUIRE(d.number_of_nodes() == 8);
    REQUIRE(orc.target_size() == 96);

    // S := ReesMatrixSemigroup(SymmetricGroup(4), [[(1, 2, 3), ()], [(),
    // ()]]);
    p.init();
    p.alphabet("abc");
    presentation::add_rule(p, "cc", "c");
    presentation::add_rule(p, "abb", "a");
    presentation::add_rule(p, "bba", "a");
    presentation::add_rule(p, "bbb", "b");
    presentation::add_rule(p, "bcb", "aca");
    presentation::add_rule(p, "aaaa", "bb");
    presentation::add_rule(p, "aaca", "bab");
    presentation::add_rule(p, "abca", "baa");
    presentation::add_rule(p, "acaa", "aab");
    presentation::add_rule(p, "baaa", "abab");
    presentation::add_rule(p, "baba", "aaab");
    presentation::add_rule(p, "baca", "acba");
    presentation::add_rule(p, "bacb", "acbb");
    presentation::add_rule(p, "bcaa", "bab");
    presentation::add_rule(p, "bcab", "aacb");
    presentation::add_rule(p, "aaaba", "acab");
    presentation::add_rule(p, "aaacb", "baab");
    presentation::add_rule(p, "aabaa", "acbb");
    presentation::add_rule(p, "aabab", "bbca");
    presentation::add_rule(p, "aacba", "acb");
    presentation::add_rule(p, "aacbb", "bca");
    presentation::add_rule(p, "abaab", "acba");
    presentation::add_rule(p, "acaba", "bca");
    presentation::add_rule(p, "acaca", "a");
    presentation::add_rule(p, "acacb", "b");
    presentation::add_rule(p, "acbaa", "baab");
    presentation::add_rule(p, "acbab", "abaa");
    presentation::add_rule(p, "acbca", "aba");
    presentation::add_rule(p, "baaba", "aacb");
    presentation::add_rule(p, "baacb", "aaba");
    presentation::add_rule(p, "bcaca", "b");
    presentation::add_rule(p, "bcacb", "aba");
    presentation::add_rule(p, "cacac", "c");
    presentation::add_rule(p, "acbbca", "aaab");
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    p.validate();
    // min. rep. o.r.c. is 7
    d = orc.presentation(p)
            .target_size(96)
            .number_of_threads(std::thread::hardware_concurrency())
            .word_graph();
    REQUIRE(d.number_of_nodes() == 7);
    REQUIRE(orc.target_size() == 96);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "078",
                          "order_preserving_monoid(5)",
                          "[extreme][sims1]") {
    auto rg = ReportGuard(true);
    auto p  = fpsemigroup::order_preserving_monoid(5);

    REQUIRE(p.rules.size() == 50);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    presentation::remove_duplicate_rules(p);
    presentation::reduce_complements(p);
    presentation::remove_trivial_rules(p);
    REQUIRE(p.rules.size() == 50);

    Sims1 S;
    REQUIRE(S.presentation(p)
                .number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(126)
            == 37'951);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "079",
                          "order_preserving_monoid(6)",
                          "[fail][sims1]") {
    // This doesn't fail it's just very extreme
    auto rg = ReportGuard(true);
    auto p  = fpsemigroup::order_preserving_monoid(6);

    REQUIRE(p.rules.size() == 72);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    presentation::remove_duplicate_rules(p);
    presentation::reduce_complements(p);
    presentation::remove_trivial_rules(p);
    REQUIRE(p.rules.size() == 72);

    Sims1 S;
    // Took 1h38min on RC office computer
    REQUIRE(S.presentation(p)
                .number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(462)
            == 91'304'735);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "080",
                          "fibonacci_group(2, 9) x 1",
                          "[quick][no-valgrind]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abcdefghiABCDEFGHI");
    p.contains_empty_word(true);
    presentation::add_inverse_rules(p, "ABCDEFGHIabcdefghi");
    REQUIRE(p.rules.size() == 36);
    presentation::add_rule(p, "ab", "c");
    presentation::add_rule(p, "bc", "d");
    presentation::add_rule(p, "cd", "e");
    presentation::add_rule(p, "de", "f");
    presentation::add_rule(p, "ef", "g");
    presentation::add_rule(p, "fg", "h");
    presentation::add_rule(p, "gh", "i");
    presentation::add_rule(p, "hi", "a");
    presentation::add_rule(p, "ia", "b");
    Sims2 S;
    S.presentation(p);
    REQUIRE(S.number_of_threads(4).number_of_congruences(12) == 6);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "081",
                          "fibonacci_group(2, 9) x 2",
                          "[extreme][sims1]") {
    using words::pow;
    auto                      rg = ReportGuard(true);
    Presentation<std::string> p;
    p.alphabet("abAB");
    p.contains_empty_word(true);
    presentation::add_inverse_rules(p, "ABab");
    presentation::add_rule(p, "Abababbab", "aBaaBaB");
    presentation::add_rule(p, "babbabbAb", "ABaaBaa");
    presentation::add_rule(p, "abbabbAbA", "BABaaBa");
    presentation::add_rule(p, "bbabbAbAA", "ABABaaB");
    presentation::add_rule(p, "babbAbAAb", "BABABaa");
    presentation::add_rule(p, "abbAbAAbA", "BBABABa");
    presentation::add_rule(p, "bbAbAAbAA", "ABBABAB");
    presentation::add_rule(p, "bAbAAbAAb", "BABBABA");
    presentation::add_rule(p, "AbAAbAAba", "BBABBAB");
    presentation::add_rule(p, "bAAbAAbab", "aBBABBA");
    presentation::add_rule(p, "AAbAAbaba", "BaBBABB");

    presentation::add_rule(p, "AAbababb", "BaaBaBBA");
    presentation::add_rule(p, "Abababba", "aBaaBaBB");
    presentation::add_rule(p, "abbabaaBaaB", "bAbAAbA");
    presentation::add_rule(p, "babaaBaaBaB", "BAbAbAA");

    Sims1 S;
    S.presentation(p);
    REQUIRE(S.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(12)
            == 6);

    Sims2 T;
    T.presentation(p);
    REQUIRE(T.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(12)
            == 6);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1", "082", "trivial group", "[fail][sims1]") {
    // This doesn't fail it's just very extreme
    auto                      rg = ReportGuard();
    Presentation<std::string> p;
    p.alphabet("rstRST");
    p.contains_empty_word(true);
    presentation::add_inverse_rules(p, "RSTrst");
    presentation::add_rule(p, "rt", "trr");
    presentation::add_rule(p, "sr", "rss");
    presentation::add_rule(p, "ts", "stt");

    ToddCoxeter tc(congruence_kind::twosided, p);
    tc.strategy(ToddCoxeter::options::strategy::felsch);
    REQUIRE(tc.number_of_classes() == 1);
    tc.shrink_to_fit();
    REQUIRE(tc.word_graph().number_of_nodes() == 1);

    Sims1 S;
    S.presentation(p);
    // Took 19min11s on RC office computer
    REQUIRE(S.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(20)
            == 1);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1", "083", "M11 x 1", "[extreme][sims1]") {
    using words::pow;
    Presentation<std::string> p;
    p.alphabet("abcABC");
    p.contains_empty_word(true);
    presentation::add_inverse_rules(p, "ABCabc");
    presentation::add_rule(p, pow("a", 6), pow("A", 5));
    presentation::add_rule(p, pow("b", 5), "");
    presentation::add_rule(p, pow("c", 4), "");
    presentation::add_rule(p, "aca", "CAC");
    presentation::add_rule(p, "bc", "cbb");
    presentation::add_rule(p, "Aba", "aab");
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    presentation::balance_no_checks(p, "abcABC", "ABCabc");

    REQUIRE(presentation::longest_subword_reducing_length(p) == "aa");
    presentation::replace_word_with_new_generator(p, "aa");

    ToddCoxeter tc(congruence_kind::twosided, p);
    REQUIRE(tc.number_of_classes() == 7'920);

    Sims1 S;
    S.presentation(p);
    REQUIRE(S.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(12)
            == 24);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1", "084", "M11 x 2", "[extreme][sims1]") {
    using words::pow;
    Presentation<std::string> p;
    p.alphabet("abcABC");
    p.contains_empty_word(true);
    presentation::add_inverse_rules(p, "ABCabc");
    presentation::add_rule(p, pow("b", 5), "");
    presentation::add_rule(p, pow("c", 4), "");
    presentation::add_rule(p, "acacac", "");
    presentation::add_rule(p, "bc", "cbb");
    presentation::add_rule(p, "ba", "aaab");
    presentation::add_rule(p, "aabba", "bb");
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);

    REQUIRE(presentation::longest_subword_reducing_length(p) == "bb");
    presentation::replace_word_with_new_generator(p, "bb");

    ToddCoxeter tc(congruence_kind::twosided, p);
    REQUIRE(tc.number_of_classes() == 7'920);

    Sims1 S;
    S.presentation(p);
    REQUIRE(S.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(16)
            == 24);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "085",
                          "JonesMonoid(4)",
                          "[extreme][sims1]") {
    using words::pow;
    Presentation<std::string> p
        = to_presentation<std::string>(temperley_lieb_monoid(4));
    REQUIRE(p.contains_empty_word());
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);

    Sims1 S;
    S.presentation(p);
    REQUIRE(S.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(10)
            == 69);
    REQUIRE(S.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(11)
            == 74);
  }

  // To keep GCC from complaining
  WordGraph<uint32_t> find_quotient(Presentation<std::string> const& p,
                                    size_t                           skip);

  WordGraph<uint32_t> find_quotient(Presentation<std::string> const& p,
                                    size_t                           skip) {
    SuppressReportFor suppressor("FroidurePin");

    auto   T              = Sims1(p);
    size_t skipped_so_far = 0;

    auto hook = [&](word_graph_type const& x) {
      auto first = (T.presentation().contains_empty_word() ? 0 : 1);
      auto S     = to_froidure_pin<Transf<0, node_type>>(
          x, first, x.number_of_active_nodes());
      if (T.presentation().contains_empty_word()) {
        auto id = one(S.generator(0));

        if (!S.contains(id)) {
          S.add_generator(id);
        }
      }
      if (S.size() == 120) {
        auto scc = Gabow(S.right_cayley_graph());
        if (scc.number_of_components() != 26) {
          return false;
        }
        scc.init(S.left_cayley_graph());
        if (skipped_so_far == skip && scc.number_of_components() == 26) {
          skipped_so_far++;
          return true;
        }
      }
      return false;
    };

    auto result = T.find_if(120, hook);
    return result;
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "087",
                          "2-sylvester monoid",
                          "[fail][sims1]") {
    Presentation<word_type> p;
    p.alphabet(3);
    p.rules = {100_w,       010_w,       200_w,       020_w,       201_w,
               021_w,       211_w,       121_w,       1010_w,      0110_w,
               1020_w,      0120_w,      2020_w,      0220_w,      2021_w,
               0221_w,      2101_w,      1021_w,      2121_w,      1221_w,
               10110_w,     01110_w,     10120_w,     01120_w,     10210_w,
               01210_w,     10220_w,     01220_w,     20220_w,     02220_w,
               20221_w,     02221_w,     21021_w,     10221_w,     21221_w,
               12221_w,     101110_w,    011110_w,    101120_w,    011120_w,
               101210_w,    011210_w,    101220_w,    011220_w,    102120_w,
               012120_w,    102210_w,    012210_w,    102220_w,    012220_w,
               202220_w,    022220_w,    202221_w,    022221_w,    210221_w,
               102221_w,    212221_w,    122221_w,    1011210_w,   0111210_w,
               1012120_w,   0112120_w,   1012210_w,   0112210_w,   1021220_w,
               0121220_w,   1022120_w,   0122120_w,   1022210_w,   0122210_w,
               10112210_w,  01112210_w,  10122120_w,  01122120_w,  10122210_w,
               01122210_w,  10221220_w,  01221220_w,  10222120_w,  01222120_w,
               101122210_w, 011122210_w, 101222120_w, 011222120_w, 101222210_w,
               011222210_w, 102221220_w, 012221220_w};

    auto S = Sims1(p);
    // NOTE: Never ran to completion, there should be a non-zero number of
    // congruences.
    REQUIRE(S.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(31)
            == 0);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1", "088", "Brauer monoid", "[fail][sims1]") {
    // This doesn't fail it's just very extreme
    auto          p = brauer_monoid(5);
    MinimalRepOrc orc;
    auto          d = orc.presentation(brauer_monoid(5))
                 .target_size(945)
                 .number_of_threads(std::thread::hardware_concurrency())
                 // The following are pairs of words in the GAP BrauerMonoid
                 // that generate the minimal 2-sided congruences of
                 // BrauerMonoid(5), the generating sets are not the same
                 // though and so this doesn't work.
                 // sims::add_excluded_pair(d, 201002_w, 00201002_w)
                 // sims::add_excluded_pair(d, 00102002_w, 001020020_w)
                 .word_graph();

    // sigma_i = (i, i + 1)
    // theta_i = Bipartition([[i, i + 1], [-i, -(i + 1)], [j, -j], j neq i]);

    REQUIRE(d.number_of_nodes() == 46);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "089",
                          "partial Brauer monoid",
                          "[fail][sims1]") {
    // This doesn't fail it's just very extreme
    std::array<uint64_t, 6> const sizes       = {0, 2, 10, 76, 764, 9496};
    std::array<uint64_t, 6> const min_degrees = {0, 2, 6, 14, 44, 143};

    for (size_t n = 1; n < 5; ++n) {
      std::cout << std::string(80, '#') << std::endl;
      auto          p = fpsemigroup::partial_brauer_monoid(n);
      MinimalRepOrc orc;
      auto          d = orc.presentation(p)
                   .target_size(sizes[n])
                   .number_of_threads(std::thread::hardware_concurrency())
                   .word_graph();

      // sigma_i = (i, i + 1)
      // theta_i = Bipartition([[i, i + 1], [-i, -(i + 1)], [j, -j], j neq
      // i]);

      REQUIRE(d.number_of_nodes() == min_degrees[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "090",
                          "possible full transf. monoid 8",
                          "[extreme][sims1]") {
    auto                    rg = ReportGuard(true);
    Presentation<word_type> p;
    p.rules = std::vector<word_type>(
        {00_w,         {},           11_w,        {},         22_w,     {},
         33_w,         {},           44_w,        {},         55_w,     {},
         66_w,         {},           101_w,       010_w,      212_w,    121_w,
         323_w,        232_w,        434_w,       343_w,      545_w,    454_w,
         656_w,        565_w,        606_w,       060_w,      2010_w,   0102_w,
         3010_w,       0103_w,       4010_w,      0104_w,     5010_w,   0105_w,
         6010_w,       0106_w,       1210_w,      0121_w,     3121_w,   1213_w,
         4121_w,       1214_w,       5121_w,      1215_w,     6121_w,   1216_w,
         2320_w,       0232_w,       2321_w,      1232_w,     4232_w,   2324_w,
         5232_w,       2325_w,       6232_w,      2326_w,     3430_w,   0343_w,
         3431_w,       1343_w,       3432_w,      2343_w,     5343_w,   3435_w,
         6343_w,       3436_w,       4540_w,      0454_w,     4541_w,   1454_w,
         4542_w,       2454_w,       4543_w,      3454_w,     6454_w,   4546_w,
         5650_w,       0565_w,       5651_w,      1565_w,     5652_w,   2565_w,
         5653_w,       3565_w,       5654_w,      4565_w,     6061_w,   1606_w,
         6062_w,       2606_w,       6063_w,      3606_w,     6064_w,   4606_w,
         6065_w,       5606_w,       071654321_w, 16543217_w, 217121_w, 17171_w,
         7010270102_w, 0102720107_w, 7010701_w,   1070170_w});
    p.alphabet_from_rules();
    auto q = full_transformation_monoid(8);

    std::array<uint64_t, 9> const num = {0, 1, 2, 3, 3, 3, 3, 3, 11};
    Sims1                         s(p);
    for (size_t n = 1; n < num.size(); ++n) {
      s.presentation(q);
      REQUIRE(s.number_of_threads(std::thread::hardware_concurrency())
                  .number_of_congruences(n)
              == num[n]);
      s.presentation(p);
      REQUIRE(s.number_of_threads(std::thread::hardware_concurrency())
                  .number_of_congruences(n)
              == num[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "091",
                          "free semilattice n = 8",
                          "[fail][sims1]") {
    Presentation<std::string> p;
    p.alphabet("abcdef");
    presentation::add_rule(p, "a^2"_p, "a");
    presentation::add_rule(p, "b^2"_p, "b");
    presentation::add_rule(p, "ba"_p, "ab");
    presentation::add_rule(p, "c^2"_p, "c");
    presentation::add_rule(p, "ca"_p, "ac");
    presentation::add_rule(p, "cb"_p, "bc");
    presentation::add_rule(p, "d^2"_p, "d");
    presentation::add_rule(p, "da"_p, "ad");
    presentation::add_rule(p, "db"_p, "bd");
    presentation::add_rule(p, "dc"_p, "cd");
    presentation::add_rule(p, "e^2"_p, "e");
    presentation::add_rule(p, "ea"_p, "ae");
    presentation::add_rule(p, "eb"_p, "be");
    presentation::add_rule(p, "ec"_p, "ce");
    presentation::add_rule(p, "ed"_p, "de");
    presentation::add_rule(p, "f^2"_p, "f");
    presentation::add_rule(p, "fa"_p, "af");
    presentation::add_rule(p, "fb"_p, "bf");
    presentation::add_rule(p, "fc"_p, "cf");
    presentation::add_rule(p, "fd"_p, "df");
    presentation::add_rule(p, "fe"_p, "ef");
    Sims1 s(p);

    // NOTE: Never ran to completion, there should be a non-zero number of
    // congruences.
    REQUIRE(s.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(std::pow(2, 6))
            == 0);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "092",
                          "temperley_lieb_monoid(4) from presentation",
                          "[quick][sims2][low-index]") {
    auto  rg = ReportGuard(false);
    Sims2 S;
    S.presentation(fpsemigroup::temperley_lieb_monoid(4));
    REQUIRE(S.number_of_congruences(14) == 9);

    auto p = fpsemigroup::temperley_lieb_monoid(4);
    presentation::reverse(p);
    Sims1 T(p);
    REQUIRE(T.number_of_congruences(14) == 79);

    std::atomic_size_t count = 0;
    T.for_each(14, [&](auto const& wg) {
      count += sims::is_two_sided_congruence_no_checks(T.presentation(), wg);
    });
    REQUIRE(count == 9);
  }

  // Takes approx. 13.5s in debug mode.
  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "093",
                          "2-sided T_4",
                          "[standard][sims2][no-valgrind][no-coverage]") {
    auto  rg = ReportGuard(false);
    Sims2 S(fpsemigroup::full_transformation_monoid(4, author::Iwahori));

    REQUIRE(S.number_of_congruences(256) == 11);  // Verified with GAP
    // sims::dot_poset("example-093-T-4-2-sided", S.cbegin(256), S.cend(256));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "094",
                          "2-sided T_4 Iwahori presentation",
                          "[standard][sims2][low-index][no-valgrind]") {
    auto  rg = ReportGuard(false);
    Sims2 S(fpsemigroup::full_transformation_monoid(
        4, fpsemigroup::author::Iwahori));
    REQUIRE(S.number_of_congruences(256) == 11);
  }

  // Not sure if the next test case runs to completion or not, JDM ran this
  // for 2m30s and it didn't terminate.
  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "095",
                          "2-sided T_4 Aizenstat presentation",
                          "[fail][sims2][low-index]") {
    auto  rg = ReportGuard(true);
    Sims2 S(fpsemigroup::full_transformation_monoid(
        4, fpsemigroup::author::Aizenstat));
    S.number_of_threads(std::thread::hardware_concurrency());
    // The below test takes too long to terminate
    REQUIRE(S.number_of_congruences(256) == 11);
  }

  // Not sure if the next test case runs to completion or not, JDM ran this
  // for 2m30s and it didn't terminate.
  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "096",
                          "2-sided S_6 Burnside+Miller presentation",
                          "[fail][sims2][low-index]") {
    auto  rg = ReportGuard(true);
    Sims2 S(fpsemigroup::symmetric_group(
        7, fpsemigroup::author::Burnside + fpsemigroup::author::Miller));
    S.number_of_threads(std::thread::hardware_concurrency());
    REQUIRE(S.number_of_congruences(720) == 3);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "097",
                          "2-sided CI_4 Fernandes presentation",
                          "[standard][sims2][low-index]") {
    auto  rg = ReportGuard(false);
    Sims2 S(fpsemigroup::cyclic_inverse_monoid(
        4, fpsemigroup::author::Fernandes, 0));
    REQUIRE(S.number_of_congruences(61) == 14);
    S.presentation(fpsemigroup::cyclic_inverse_monoid(
        4, fpsemigroup::author::Fernandes, 1));
    REQUIRE(S.number_of_congruences(61) == 14);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "098",
                          "2-sided CI_4 Froidure-Pin presentation",
                          "[standard][sims2][low-index]") {
    auto                  rg = ReportGuard(false);
    FroidurePin<PPerm<4>> T;
    T.add_generator(make<PPerm<4>>({1, 2, 3, 0}));
    T.add_generator(make<PPerm<4>>({1, 2, 3}, {1, 2, 3}, 4));
    T.add_generator(make<PPerm<4>>({0, 2, 3}, {0, 2, 3}, 4));
    T.add_generator(make<PPerm<4>>({0, 1, 3}, {0, 1, 3}, 4));
    T.add_generator(make<PPerm<4>>({0, 1, 2}, {0, 1, 2}, 4));
    REQUIRE(T.size() == 61);

    auto p = to_presentation<word_type>(T);

    Sims2 S(p);
    REQUIRE(S.number_of_congruences(61) == 14);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "099",
                          "2-sided (2,3,7) triangle group",
                          "[quick][sims2][low-index][no-valgrind]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("xy");
    presentation::add_rule(p, "xx", "");
    presentation::add_rule(p, "yyy", "");
    presentation::add_rule(p, "xyxyxyxyxyxyxy", "");
    Sims2 S(p);
    // Smallest non-trivial homomorphic image has size 168, see
    // https://mathoverflow.net/questions/180231/
    // for more details
    REQUIRE(S.number_of_congruences(168) == 2);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "100",
                          "2-sided Heineken group",
                          "[extreme][sims2][low-index]") {
    auto                      rg = ReportGuard(true);
    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("xXyY");
    presentation::add_inverse_rules(p, "XxYy");
    presentation::add_rule(p, "yXYYxyYYxyyXYYxyyXyXYYxy", "x");
    presentation::add_rule(p, "YxyyXXYYxyxYxyyXYXyXYYxxyyXYXyXYYxyx", "y");

    Sims2 S(p);
    REQUIRE(S.number_of_threads(1).number_of_congruences(50) == 1);
    REQUIRE(S.number_of_threads(2).number_of_congruences(50) == 1);
    REQUIRE(S.number_of_threads(4).number_of_congruences(50) == 1);
    REQUIRE(S.number_of_threads(8).number_of_congruences(50) == 1);
    REQUIRE(S.number_of_threads(16).number_of_congruences(50) == 1);
    REQUIRE(S.number_of_threads(32).number_of_congruences(50) == 1);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "101",
                          "2-sided Catalan monoid n=4",
                          "[quick][sims2][low-index][no-valgrind]") {
    auto                   rg = ReportGuard(false);
    FroidurePin<Transf<4>> S;
    S.add_generator(make<Transf<4>>({0, 1, 2, 3}));
    S.add_generator(make<Transf<4>>({0, 0, 2, 3}));
    S.add_generator(make<Transf<4>>({0, 1, 1, 3}));
    S.add_generator(make<Transf<4>>({0, 1, 2, 2}));
    REQUIRE(S.size() == 14);
    auto p = to_presentation<word_type>(S);

    Sims2 C(p);
    REQUIRE(C.number_of_threads(1).number_of_congruences(S.size()) == 133);
    REQUIRE(C.number_of_threads(2).number_of_congruences(S.size()) == 133);
    REQUIRE(C.number_of_threads(4).number_of_congruences(S.size()) == 133);
    REQUIRE(C.number_of_threads(8).number_of_congruences(S.size()) == 133);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "102",
                          "2-sided Heineken monoid",
                          "[extreme][sims2][low-index]") {
    auto                      rg = ReportGuard(true);
    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("xyXY");
    presentation::add_rule(p, "yXYYxyYYxyyXYYxyyXyXYYxyX", "");
    presentation::add_rule(p, "YxyyXXYYxyxYxyyXYXyXYYxxyyXYXyXYYxyxY", "");
    // REQUIRE(presentation::to_gap_string(p, "S") == "");
    Sims2 S(p);
    sims::add_included_pair(S, 0_w, 2_w);
    REQUIRE(S.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(8)
            == 9);
    check_congruence_count_with_free_object(S, 8, 9);
    S.clear_included_pairs();
    sims::add_excluded_pair(S, 0_w, 2_w);
    REQUIRE(S.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(8)
            == 63 - 9);
    check_congruence_count_with_free_object(S, 8, 63 - 9);
  }

  // Takes approx. 1 minute
  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "103",
                          "2-sided Fibonacci(2, 9)",
                          "[extreme][sims2][low-index]") {
    auto                      rg = ReportGuard(true);
    Presentation<std::string> p;
    p.alphabet("abAB");
    p.contains_empty_word(true);
    presentation::add_inverse_rules(p, "ABab");
    presentation::add_rule(p, "Abababbab", "aBaaBaB");
    presentation::add_rule(p, "babbabbAb", "ABaaBaa");
    presentation::add_rule(p, "abbabbAbA", "BABaaBa");
    presentation::add_rule(p, "bbabbAbAA", "ABABaaB");
    presentation::add_rule(p, "babbAbAAb", "BABABaa");
    presentation::add_rule(p, "abbAbAAbA", "BBABABa");
    presentation::add_rule(p, "bbAbAAbAA", "ABBABAB");
    presentation::add_rule(p, "bAbAAbAAb", "BABBABA");
    presentation::add_rule(p, "AbAAbAAba", "BBABBAB");
    presentation::add_rule(p, "bAAbAAbab", "aBBABBA");
    presentation::add_rule(p, "AAbAAbaba", "BaBBABB");

    presentation::add_rule(p, "AAbababb", "BaaBaBBA");
    presentation::add_rule(p, "Abababba", "aBaaBaBB");
    presentation::add_rule(p, "abbabaaBaaB", "bAbAAbA");
    presentation::add_rule(p, "babaaBaaBaB", "BAbAbAA");

    // REQUIRE(presentation::to_gap_string(p, "S") == "");
    Sims2 S(p);
    // TODO(2): check correctness
    REQUIRE(S.number_of_threads(std::thread::hardware_concurrency())
                .number_of_congruences(64)
            == 10);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "104",
                          "2-sided one-relation baaabaaa=aba",
                          "[standard][sims2][low-index][no-coverage]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("ab");
    p.contains_empty_word(true);
    presentation::add_rule(p, "baaabaaa", "aba");

    Sims2 S(p);
    REQUIRE(S.number_of_threads(8).number_of_congruences(1) == 1);
    check_congruence_count_with_free_object(S, 1, 1);
    REQUIRE(S.number_of_threads(8).number_of_congruences(2) == 5);
    check_congruence_count_with_free_object(S, 2, 5);
    REQUIRE(S.number_of_threads(8).number_of_congruences(3) == 17);
    check_congruence_count_with_free_object(S, 3, 17);
    REQUIRE(S.number_of_threads(8).number_of_congruences(4) == 52);
    check_congruence_count_with_free_object(S, 4, 52);
    // sims::dot_poset("example-104-baaabaaa=aba", S.cbegin(4), S.cend(4));

    std::atomic_size_t count = 0;
    Sims1              T(p);
    REQUIRE(T.number_of_congruences(4) == 977);
    T.for_each(4, [&](auto const& wg) {
      count += sims::is_two_sided_congruence(T.presentation(), wg);
    });
    REQUIRE(count == 52);
    count = 0;

    presentation::reverse(p);
    T.presentation(p);
    REQUIRE(T.number_of_congruences(4) == 227);
    T.for_each(4, [&](auto const& wg) {
      count += sims::is_two_sided_congruence(T.presentation(), wg);
    });
    REQUIRE(count == 52);

    // Note that Catch's REQUIRE macro is not thread safe, see:
    // https://github.com/catchorg/Catch2/issues/99
    // as such we cannot call any function (like check_right_generating_pairs)
    // that uses REQUIRE in multiple threads.
    S.number_of_threads(1).for_each(
        5, [&S](auto const& wg) { check_two_sided_generating_pairs(S, wg); });
    S.for_each(5,
               [&S](auto const& wg) { check_right_generating_pairs(S, wg); });

    REQUIRE(S.number_of_threads(8).number_of_congruences(5) == 148);
    REQUIRE(S.number_of_threads(8).number_of_congruences(6) == 413);
    REQUIRE(S.number_of_threads(8).number_of_congruences(7) == 1'101);
    REQUIRE(S.number_of_threads(8).number_of_congruences(8) == 2'901);
    REQUIRE(S.number_of_threads(8).number_of_congruences(9) == 7'569);
    REQUIRE(S.number_of_threads(8).number_of_congruences(10) == 19'756);
    REQUIRE(S.number_of_threads(8).number_of_congruences(11) == 50'729);
    REQUIRE(S.number_of_threads(8).number_of_congruences(12) == 129'157);
    REQUIRE(S.number_of_threads(8).number_of_congruences(13) == 330'328);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "Sims2",
      "105",
      "2-sided one-relation baabbaa=a",
      "[extreme][sims2][low-index][no-valgrind][no-coverage]") {
    auto                      rg = ReportGuard(true);
    Presentation<std::string> p;
    p.alphabet("ab");
    p.contains_empty_word(true);
    presentation::add_rule(p, "baabbaa", "a");

    Sims2 S(p);

    // sims::dot_poset("example-105-baabbaa=a-2-sided", S.cbegin(8),
    // S.cend(8));
    // Takes a long time to run, seems like we get all the congruences quite
    // early on, but then spend very long checking that there are no more.
    // Perhaps if we had some sort of upper bound could speed things up?
    size_t num_threads = std::thread::hardware_concurrency();
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(1) == 1);
    check_congruence_count_with_free_object(S, 1, 1);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(2) == 4);
    check_congruence_count_with_free_object(S, 2, 4);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(3) == 13);
    check_congruence_count_with_free_object(S, 3, 13);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(4) == 28);
    check_congruence_count_with_free_object(S, 4, 28);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(5) == 49);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(6) == 86);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(7) == 134);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(8) == 200);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(9) == 284);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(10) == 392);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(11) == 518);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(12) == 693);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(13) == 891);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(14)
            == 1'127);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(15)
            == 1'402);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(16)
            == 1'733);
    REQUIRE(S.number_of_threads(num_threads).number_of_congruences(17)
            == 2'094);
    // This is a little too extreme with the next 3 lines uncommented.
    // REQUIRE(S.number_of_threads(num_threads).number_of_congruences(18) ==
    // 2'531);
    // REQUIRE(S.number_of_threads(num_threads).number_of_congruences(19) ==
    // 3'012);
    // REQUIRE(S.number_of_threads(num_threads).number_of_congruences(20) ==
    // 3'574);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "106",
                          "2-sided full transformation monoid 2",
                          "[quick][sims2]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(2);
    p.contains_empty_word(true);
    presentation::add_rule(p, 00_w, {});
    presentation::add_rule(p, 01_w, 1_w);
    presentation::add_rule(p, 11_w, 1_w);
    Sims2 s(p);
    REQUIRE(s.number_of_congruences(4) == 4);  // Verified with GAP
    auto it = s.cbegin(4);
    REQUIRE(*(it++) == make<WordGraph<node_type>>(4, {{0, 0}}));          // ok
    REQUIRE(*(it++) == make<WordGraph<node_type>>(4, {{0, 1}, {1, 1}}));  // ok
    REQUIRE(*(it++)
            == make<WordGraph<node_type>>(4, {{1, 2}, {0, 2}, {2, 2}}));  // ok

    REQUIRE(*(it++)
            == make<WordGraph<node_type>>(
                4, {{1, 2}, {0, 2}, {3, 2}, {2, 2}}));  // ok
    REQUIRE(it == s.cend(4));

    s.for_each(
        4, [&s](auto const& wg) { check_two_sided_generating_pairs(s, wg); });
    auto not_in_p = rx::filter([&p](auto const& rel) {
      return !presentation::contains_rule(p, rel.first, rel.second);
    });

    it = s.cbegin(4);
    REQUIRE((sims::right_generating_pairs(*it) | rx::to_vector())
            == std::vector<relation_type>({{0_w, {}}, {1_w, {}}}));
    REQUIRE((sims::right_generating_pairs(p, *it) | rx::to_vector())
            == std::vector<relation_type>({{0_w, {}}, {1_w, {}}}));
    ++it;
    REQUIRE(
        (sims::right_generating_pairs(*it) | rx::to_vector())
        == std::vector<relation_type>({{0_w, {}}, {10_w, 1_w}, {11_w, 1_w}}));
    REQUIRE((sims::right_generating_pairs(p, *it) | rx::to_vector())
            == std::vector<relation_type>({{0_w, {}}, {10_w, 1_w}}));
    ++it;
    REQUIRE((sims::right_generating_pairs(*it) | rx::to_vector())
            == std::vector<relation_type>(
                {{00_w, {}}, {01_w, 1_w}, {10_w, 1_w}, {11_w, 1_w}}));
    REQUIRE((sims::right_generating_pairs(p, *it) | not_in_p | rx::to_vector())
            == std::vector<relation_type>({{10_w, 1_w}}));

    // Note that all the rules below follow from the rules in the
    // presentation, and so this congruence is the trivial one.
    ++it;
    REQUIRE((sims::right_generating_pairs(*it) | rx::to_vector())
            == std::vector<relation_type>({{00_w, {}},
                                           {01_w, 1_w},
                                           {11_w, 1_w},
                                           {100_w, 1_w},
                                           {101_w, 1_w}}));
    REQUIRE((sims::right_generating_pairs(p, *it) | not_in_p | rx::to_vector())
            == std::vector<relation_type>({{100_w, 1_w}}));
    REQUIRE(
        (sims::two_sided_generating_pairs(p, *it) | not_in_p | rx::to_vector())
        == std::vector<relation_type>({{100_w, 1_w}}));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2", "107", "2-sided example", "[quick][sims1]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.contains_empty_word(true);
    p.alphabet(01_w);
    presentation::add_rule(p, 000_w, 0_w);
    presentation::add_rule(p, 11_w, 1_w);
    presentation::add_rule(p, 0101_w, 0_w);
    Sims2 s(p);

    REQUIRE(s.number_of_congruences(4) == 6);  // Verified with GAP
    auto it = s.cbegin(5);
    // Verified in 000
    REQUIRE(*(it++) == make<WordGraph<node_type>>(5, {{0, 0}}));
    REQUIRE(*(it++) == make<WordGraph<node_type>>(5, {{1, 0}, {1, 1}}));
    REQUIRE(*(it++) == make<WordGraph<node_type>>(5, {{1, 1}, {1, 1}}));
    REQUIRE(*(it++) == make<WordGraph<node_type>>(5, {{1, 2}, {1, 1}, {1, 2}}));
    REQUIRE(*(it++) == make<WordGraph<node_type>>(5, {{1, 2}, {1, 1}, {2, 2}}));
    REQUIRE(*(it++)
            == make<WordGraph<node_type>>(5, {{1, 2}, {1, 1}, {3, 2}, {3, 3}}));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "108",
                          "2-sided full transf. monoid 3",
                          "[quick][sims2]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abc");
    p.contains_empty_word(true);
    presentation::add_rule(p, "b^2"_p, ""_p);
    presentation::add_rule(p, "bc"_p, "ac"_p);
    presentation::add_rule(p, "c^2"_p, "c"_p);
    presentation::add_rule(p, "a^3"_p, ""_p);
    presentation::add_rule(p, "a^2b"_p, "ba"_p);
    presentation::add_rule(p, "aba"_p, "b"_p);
    presentation::add_rule(p, "baa"_p, "ab"_p);
    presentation::add_rule(p, "bab"_p, "aa"_p);
    presentation::add_rule(p, "bac"_p, "c"_p);
    presentation::add_rule(p, "cac"_p, "cb"_p);
    presentation::add_rule(p, "aca^2c"_p, "ca^2c"_p);
    presentation::add_rule(p, "ca^2cb"_p, "ca^2ca"_p);
    presentation::add_rule(p, "ca^2cab"_p, "ca^2c"_p);
    Sims2 s(p);
    REQUIRE(s.number_of_congruences(27) == 7);  // Verified with GAP

    auto it = s.cbegin(27);

    REQUIRE(*(it++) == make<WordGraph<node_type>>(27, {{0, 0, 0}}));  // ok
    REQUIRE(*(it++)
            == make<WordGraph<node_type>>(27, {{0, 0, 1}, {1, 1, 1}}));  // ok
    REQUIRE(*(it++)
            == make<WordGraph<node_type>>(
                27, {{0, 1, 2}, {1, 0, 2}, {2, 2, 2}}));  // ok
    REQUIRE(*(it++)
            == make<WordGraph<node_type>>(27,
                                          {{1, 2, 3},
                                           {4, 5, 3},
                                           {6, 0, 3},
                                           {3, 3, 3},
                                           {0, 6, 3},
                                           {2, 1, 3},
                                           {5, 4, 3}}));  // ok
    REQUIRE(*(it++)
            == make<WordGraph<node_type>>(27,
                                          {{1, 2, 3},
                                           {4, 5, 6},
                                           {7, 0, 6},
                                           {8, 3, 3},
                                           {0, 7, 9},
                                           {2, 1, 9},
                                           {10, 6, 6},
                                           {5, 4, 3},
                                           {11, 11, 3},
                                           {12, 9, 9},
                                           {13, 13, 6},
                                           {3, 8, 14},
                                           {15, 15, 9},
                                           {6, 10, 14},
                                           {14, 14, 14},
                                           {9, 12, 14}}));  // ok
    REQUIRE(*(it++)
            == make<WordGraph<node_type>>(
                27, {{1, 2, 3},    {4, 5, 6},    {7, 0, 6},    {8, 9, 3},
                     {0, 7, 10},   {2, 1, 10},   {11, 12, 6},  {5, 4, 3},
                     {13, 14, 9},  {15, 3, 9},   {16, 17, 10}, {18, 19, 12},
                     {20, 6, 12},  {3, 15, 21},  {9, 8, 21},   {14, 13, 3},
                     {22, 23, 17}, {24, 10, 17}, {6, 20, 21},  {12, 11, 21},
                     {19, 18, 6},  {21, 21, 21}, {10, 24, 21}, {17, 16, 21},
                     {23, 22, 10}}));  // ok
    REQUIRE(*(it++)
            == make<WordGraph<node_type>>(
                27, {{1, 2, 3},    {4, 5, 6},    {7, 0, 6},    {8, 9, 3},
                     {0, 7, 10},   {2, 1, 10},   {11, 12, 6},  {5, 4, 3},
                     {13, 14, 9},  {15, 3, 9},   {16, 17, 10}, {18, 19, 12},
                     {20, 6, 12},  {3, 15, 21},  {9, 8, 21},   {14, 13, 3},
                     {22, 23, 17}, {24, 10, 17}, {6, 20, 21},  {12, 11, 21},
                     {19, 18, 6},  {25, 25, 21}, {10, 24, 21}, {17, 16, 21},
                     {23, 22, 10}, {26, 21, 25}, {21, 26, 21}}));  // ok
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "109",
                          "2-sided 2-generated free monoid",
                          "[standard][sims2][no-coverage]") {
    auto                      rg = ReportGuard(true);
    Presentation<std::string> p;
    p.alphabet("ab");
    p.contains_empty_word(true);
    Sims2 s(p);
    s.number_of_threads(4);
    REQUIRE(s.number_of_congruences(1) == 1);
    REQUIRE(s.number_of_congruences(2) == 7);    // verified with GAP
    REQUIRE(s.number_of_congruences(3) == 27);   // verified with GAP
    REQUIRE(s.number_of_congruences(4) == 94);   // verified with GAP
    REQUIRE(s.number_of_congruences(5) == 275);  // verified with GAP
    REQUIRE(s.number_of_congruences(6) == 833);
    REQUIRE(s.number_of_congruences(7) == 2'307);
    REQUIRE(s.number_of_congruences(8) == 6'488);
    REQUIRE(s.number_of_congruences(9) == 18'207);
    REQUIRE(s.number_of_congruences(10) == 52'960);
    REQUIRE(s.number_of_congruences(11) == 156'100);
    REQUIRE(s.number_of_congruences(12) == 462'271);
    REQUIRE(s.number_of_congruences(13) == 1'387'117);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "110",
                          "symmetric inverse monoid (Gay)",
                          "[standard][sims2]") {
    auto rg = ReportGuard(true);
    auto p  = symmetric_inverse_monoid(5, fpsemigroup::author::Gay);
    presentation::remove_duplicate_rules(p);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);

    Sims2 C(p);
    REQUIRE(C.number_of_threads(1).number_of_congruences(1'546) == 14);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "111",
                          "2-sided congruence-free monoid n=3",
                          "[quick][sims2][no-valgrind]") {
    auto rg = ReportGuard(false);
    // Presentation taken from
    // Al-Kharousi, F., Cain, A.J., Maltcev, V. et al.
    // A countable family of finitely presented infinite congruence-free
    // monoids https://doi.org/10.14232/actasm-013-028-z
    Presentation<std::string> p;
    p.alphabet("abcdz");
    p.contains_empty_word(true);
    presentation::add_zero_rules(p, 'z');
    presentation::add_rule(p, "ac", "");
    presentation::add_rule(p, "db", "");
    presentation::add_rule(p, "dc", "");
    // n = 3
    presentation::add_rule(p, "dab", "");
    presentation::add_rule(p, "da^2b"_p, "");
    presentation::add_rule(p, "a^3b"_p, "z");

    Sims2 s(p);
    s.number_of_threads(1);
    REQUIRE(s.number_of_congruences(1) == 1);
    REQUIRE(s.number_of_congruences(2) == 1);
    REQUIRE(s.number_of_congruences(3) == 1);
    REQUIRE(s.number_of_congruences(4) == 1);
    REQUIRE(s.number_of_congruences(5) == 1);
    REQUIRE(s.number_of_congruences(10) == 1);
    REQUIRE(s.number_of_congruences(20) == 1);
    REQUIRE(s.number_of_congruences(30) == 1);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "112",
                          "2-sided congruence-free monoid n=8",
                          "[standard][sims2][no-valgrind]") {
    auto rg = ReportGuard(false);
    // Presentation taken from
    // Al-Kharousi, F., Cain, A.J., Maltcev, V. et al.
    // A countable family of finitely presented infinite congruence-free
    // monoids https://doi.org/10.14232/actasm-013-028-z
    Presentation<std::string> p;
    p.alphabet("abcdz");
    p.contains_empty_word(true);
    presentation::add_zero_rules(p, 'z');
    presentation::add_rule(p, "ac", "");
    presentation::add_rule(p, "db", "");
    presentation::add_rule(p, "dc", "");
    // n = 8
    presentation::add_rule(p, "dab", "");
    presentation::add_rule(p, "da^2b"_p, "");
    presentation::add_rule(p, "da^3b"_p, "");
    presentation::add_rule(p, "da^4b"_p, "");
    presentation::add_rule(p, "da^5b"_p, "");
    presentation::add_rule(p, "da^6b"_p, "");
    presentation::add_rule(p, "da^7b"_p, "");
    presentation::add_rule(p, "a^8b"_p, "z");

    Sims2 s(p);
    s.number_of_threads(1);
    REQUIRE(s.number_of_congruences(1) == 1);
    REQUIRE(s.number_of_congruences(2) == 1);
    REQUIRE(s.number_of_congruences(3) == 1);
    REQUIRE(s.number_of_congruences(4) == 1);
    REQUIRE(s.number_of_congruences(5) == 1);
    REQUIRE(s.number_of_congruences(10) == 1);
    REQUIRE(s.number_of_congruences(20) == 1);
    REQUIRE(s.number_of_congruences(30) == 1);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "113",
                          "2-sided bicyclic monoid",
                          "[quick][sims2][no-valgrind]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("ab");
    p.contains_empty_word(true);
    presentation::add_rule(p, "ab"_p, ""_p);
    Sims2 s(p);
    s.number_of_threads(1);
    for (size_t i = 1; i < 50; ++i) {
      REQUIRE(s.number_of_congruences(i) == i);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "114",
                          "2-sided 2-generated free commutative monoid",
                          "[standard][sims2][no-valgrind]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("ab");
    p.contains_empty_word(true);
    presentation::add_rule(p, "ab", "ba");

    Sims2 s(p);
    s.number_of_threads(4);
    REQUIRE(s.number_of_congruences(1) == 1);
    REQUIRE(s.number_of_congruences(2) == 7);  // verified by hand
    REQUIRE(s.number_of_congruences(3) == 25);
    REQUIRE(s.number_of_congruences(4) == 76);
    REQUIRE(s.number_of_congruences(5) == 184);
    REQUIRE(s.number_of_congruences(6) == 432);
    REQUIRE(s.number_of_congruences(7) == 892);

    REQUIRE(s.number_of_congruences(8) == 1'800);
    REQUIRE(s.number_of_congruences(9) == 3'402);
    REQUIRE(s.number_of_congruences(10) == 6'280);
    REQUIRE(s.number_of_congruences(11) == 11'051);

    REQUIRE(s.number_of_congruences(12) == 19'245);
    REQUIRE(s.number_of_congruences(13) == 32'299);

    check_meets_and_joins(s.cbegin(5), s.cend(5));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "115",
                          "free semilattice n = 8",
                          "[standard][sims1][no-coverage]") {
    auto rg = ReportGuard(true);
    // https://oeis.org/A102894
    constexpr std::array<size_t, 6> results = {0, 1, 4, 45, 2'271, 1'373'701};
    for (auto A : {"a", "ab", "abc", "abcd", "abcde"}) {
      Presentation<std::string> p;
      p.alphabet(A);
      presentation::add_commutes_rules_no_checks(p, p.alphabet());
      presentation::add_idempotent_rules_no_checks(p, p.alphabet());
      Sims2 s(p);

      size_t const n = p.alphabet().size();
      REQUIRE(s.number_of_threads(4).number_of_congruences(std::pow(2, n))
              == results[n]);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "116",
                          "2-sided 2-generated free semigroup",
                          "[quick][sims2]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("ab");
    p.contains_empty_word(false);
    Sims2 s(p);
    s.number_of_threads(4);
    // Number of congruences with up to 7 classes given in:
    // A. Bailey, M. Finn-Sell and R. Snocken
    // "SUBSEMIGROUP, IDEAL AND CONGRUENCE GROWTH OF FREE SEMIGROU6S"
    REQUIRE(s.number_of_congruences(1) == 1);
    REQUIRE(s.number_of_congruences(2) == 11);  // From Bailey et al
    REQUIRE(s.number_of_congruences(3) == 51);  // From Bailey et al
    // sims::dot_poset("free-semigroup-2-sided", s.cbegin(3), s.cend(3));
    REQUIRE(s.number_of_congruences(4) == 200);  // From Bailey et al

    // REQUIRE(s.number_of_congruences(5) == 657);    // From Bailey et al
    // REQUIRE(s.number_of_congruences(6) == 2'037);  // From Bailey et al
    // REQUIRE(s.number_of_congruences(7) == 5'977);  // From Bailey et al
    // sims::add_included_pair(s, 0_w, 1_w);
    // REQUIRE(s.number_of_congruences(8) == 36);
    // sims::add_excluded_pair(s, 0_w, 1_w);
    // REQUIRE(s.number_of_congruences(8) == 0);
    // s.clear_included_pairs();
    // REQUIRE(s.number_of_congruences(8) == 17381);
    // s.clear_excluded_pairs();
    // REQUIRE(s.number_of_congruences(8) == 17381 + 36);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "117",
                          "1-sided ideals 2-generated free semigroup",
                          "[quick][sims1]") {
    Presentation<std::string> p;
    p.alphabet("ab");
    p.contains_empty_word(true);

    SimsRefinerIdeals ip(p);

    Sims1 s(p);
    s.number_of_threads(4);
    s.add_pruner(ip);
    REQUIRE(s.number_of_congruences(1) == 1);
    REQUIRE(s.number_of_congruences(2) == 2);

    auto it = s.cbegin(2);
    REQUIRE(*(it++) == make<WordGraph<uint32_t>>(2, {{0, 0}}));
    REQUIRE(*(it++) == make<WordGraph<uint32_t>>(2, {{1, 1}, {1, 1}}));

    REQUIRE(s.number_of_congruences(3) == 4);
    REQUIRE(s.number_of_congruences(4) == 9);
    REQUIRE(s.number_of_congruences(5) == 23);
    REQUIRE(s.number_of_congruences(6) == 65);
    REQUIRE(s.number_of_congruences(7) == 197);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "118",
                          "1-sided ideals partition monoid, n = 2",
                          "[quick][sims1]") {
    Presentation<word_type> p;
    p.contains_empty_word(true);

    p.alphabet(012_w);
    presentation::add_rule(p, 00_w, ""_w);
    presentation::add_rule(p, 02_w, 2_w);
    presentation::add_rule(p, 11_w, 1_w);
    presentation::add_rule(p, 20_w, 2_w);
    presentation::add_rule(p, 22_w, 2_w);
    presentation::add_rule(p, 121_w, 1_w);
    presentation::add_rule(p, 212_w, 2_w);
    presentation::add_rule(p, 0101_w, 101_w);
    presentation::add_rule(p, 1010_w, 101_w);

    // KnuthBendix kb(congruence_kind::twosided, p);
    // REQUIRE(kb.number_of_classes() == 15);

    SimsRefinerIdeals ip(p);

    Sims1 s(p);
    s.add_pruner(ip);
    REQUIRE(s.number_of_congruences(15) == 15);  // correct value is 15
    REQUIRE(s.number_of_threads(2).number_of_congruences(15)
            == 15);  // correct value is 15
    REQUIRE(s.number_of_threads(4).number_of_congruences(15)
            == 15);  // correct value is 15
    REQUIRE(s.number_of_threads(8).number_of_congruences(15)
            == 15);  // correct value is 15

    // p  = partition_monoid(3, author::Machine);
    // ip = SimsRefinerIdeals(to_presentation<std::string>(p));
    // s.init(p).add_pruner(ip);
    // REQUIRE(s.number_of_congruences(203) == 5767);  // checked in GAP
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "119",
                          "2-sided ideals Jura's example",
                          "[quick][sims1][no-valgrind]") {
    Presentation<std::string> p;
    p.alphabet("ab");
    p.contains_empty_word(false);
    presentation::add_rule(p, "aaa", "bb");
    presentation::add_rule(p, "aab", "ba");

    auto q = to_presentation<word_type>(p);

    ToddCoxeter tc(congruence_kind::twosided, q);
    REQUIRE(tc.number_of_classes() == 12);

    Sims2             s(q);
    SimsRefinerIdeals ip(s.presentation());
    s.add_pruner(ip);

    REQUIRE(s.number_of_congruences(1) == 1);   // computed using GAP
    REQUIRE(s.number_of_congruences(2) == 3);   // computed using GAP
    REQUIRE(s.number_of_congruences(3) == 5);   // computed using GAP
    REQUIRE(s.number_of_congruences(4) == 7);   // computed using GAP
    REQUIRE(s.number_of_congruences(5) == 9);   // computed using GAP
    REQUIRE(s.number_of_congruences(6) == 11);  // computed using GAP
    REQUIRE(s.number_of_congruences(7) == 12);  // computed using GAP
    for (size_t nr_classes = 8; nr_classes < 16; ++nr_classes)
      REQUIRE(s.number_of_congruences(nr_classes) == 12);  // computed using GAP
    REQUIRE(s.number_of_threads(2).number_of_congruences(7)
            == 12);  // computed using GAP
    REQUIRE(s.number_of_threads(4).number_of_congruences(7)
            == 12);  // computed using GAP
    REQUIRE(s.number_of_threads(8).number_of_congruences(7)
            == 12);  // computed using GAP
  }

  // about 2 seconds
  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "120",
                          "order_preserving_monoid(5)",
                          "[standard][sims1][no-coverage]") {
    auto rg = ReportGuard(false);
    auto p  = fpsemigroup::order_preserving_monoid(5);
    REQUIRE(p.contains_empty_word());
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    presentation::remove_duplicate_rules(p);
    presentation::reduce_complements(p);
    presentation::remove_trivial_rules(p);

    // presentation::reverse(p);
    // REQUIRE(presentation::length(p) == 88);
    Sims1       s(p);
    auto const& pp = s.presentation();

    // REQUIRE(to_presentation<std::string>(s.presentation()).rules
    //         == std::vector<std::string>());
    s.number_of_threads(std::thread::hardware_concurrency());
    // REQUIRE(s.number_of_congruences(126) == 37'951);

    // auto wg_found = s.find_if(35, [&pp](auto const& wg) {
    //   return !sims::is_right_congruence(pp, wg);
    // });

    // REQUIRE(wg_found.number_of_active_nodes() == 0);
    // auto wg_expected = make<WordGraph<uint32_t>>(6,
    //                                            {{1, 0, 0, 2, 0, 0},
    //                                             {1, 0, 3, 1, 1, 1},
    //                                             {2, 2, 2, 2, 0, 2},
    //                                             {3, 3, 3, 3, 3, 1}});
    //  REQUIRE(sims::is_right_congruence(p, wg_expected));

    // REQUIRE(word_graph::is_complete(
    //     wg_found, wg_found.cbegin_nodes(), wg_found.cbegin_nodes() + 4));
    // auto i = 0;
    // //    static_assert(::libsemigroups::detail::HasLessEqual<, >::value);
    // // REQUIRE(pp.rules == std::vector<word_type>());
    // for (auto it = pp.rules.cbegin(); it < pp.rules.cend(); it += 2) {
    //   fmt::print("i = {}\n", i);
    //   REQUIRE(word_graph::is_compatible_no_checks(wg_found,
    //                                     wg_found.cbegin_nodes(),
    //                                     wg_found.cbegin_nodes() + 4,
    //                                     it,
    //                                     it + 1));
    //   i += 1;
    // }

    std::atomic_uint64_t result = 0;
    s.for_each(125, [&result, &pp](auto const& wg) {
      if (sims::is_maximal_right_congruence(pp, wg)) {
        result++;
      }
    });
    REQUIRE(result == 31);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "121",
                          "order_preserving_monoid(6)",
                          "[fail][sims1]") {
    auto rg = ReportGuard(false);
    auto p  = fpsemigroup::order_preserving_monoid(6);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    presentation::remove_duplicate_rules(p);
    presentation::reduce_complements(p);
    presentation::remove_trivial_rules(p);
    REQUIRE(presentation::length(p) == 128);
    Sims1 s(p);
    s.number_of_threads(std::thread::hardware_concurrency());
    // RC ran this for a while and this is the last line of output before it
    // was terminated: #1: Sims1: total        90,217,309 (cong.)   |
    // 1,370,495,022,905 (nodes) #1: Sims1: diff             +2,237 (cong.) |
    // +63,535,732 (nodes) #1: Sims1: mean              3,633 (cong./s) |
    // 218,936 (node/s) #1: Sims1: time last s.      447µs (/cong.)  | 15ns
    // (/node) #1: Sims1: mean time         275µs (/cong.)  | 18ns (/node) #1:
    // Sims1: time         6h53min48s (total)   |

    // NOTE: Never ran to completion, there should be a non-zero number of
    // congruences.
    REQUIRE(s.number_of_congruences(462) == 0);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "122",
                          "partition_monoid(2)",
                          "[quick][sims1]") {
    Presentation<word_type> p;
    p.contains_empty_word(true);

    p.alphabet(123_w);
    presentation::add_rule(p, 11_w, ""_w);
    presentation::add_rule(p, 13_w, 3_w);
    presentation::add_rule(p, 22_w, 2_w);
    presentation::add_rule(p, 31_w, 3_w);
    presentation::add_rule(p, 33_w, 3_w);
    presentation::add_rule(p, 232_w, 2_w);
    presentation::add_rule(p, 323_w, 3_w);
    presentation::add_rule(p, 1212_w, 212_w);
    presentation::add_rule(p, 2121_w, 212_w);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    presentation::remove_duplicate_rules(p);
    presentation::reduce_complements(p);
    presentation::remove_trivial_rules(p);

    Sims1                s(p);
    auto                 pp     = s.presentation();
    std::atomic_uint64_t result = 0;
    s.for_each(11, [&result, &pp](auto const& wg) {
      if (sims::is_maximal_right_congruence(pp, wg)) {
        fmt::print("Index {}\n", wg.number_of_active_nodes());
        result++;
      }
    });
    REQUIRE(result == 6);
    REQUIRE(s.number_of_congruences(15) == 105);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "123",
                          "Adding and removing pruners",
                          "[quick][low-index]") {
    auto                      rg = ReportGuard(true);
    Presentation<std::string> p;
    p.alphabet("ab");
    p.contains_empty_word(false);
    presentation::add_rule(p, "aaa", "bb");
    presentation::add_rule(p, "aab", "ba");

    Sims2             s(p);
    SimsRefinerIdeals ip(s.presentation());
    s.add_pruner(ip);
    REQUIRE(s.number_of_congruences(12) == 12);  // computed using GAP
    s.clear_pruners();
    REQUIRE(s.number_of_congruences(12) == 41);  // computed using GAP
    s.add_pruner(ip);
    REQUIRE(s.number_of_congruences(12) == 12);  // computed using GAP
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "124",
                          "onesided congruence checking",
                          "[quick][low-index]") {
    Presentation<word_type> p;
    p.alphabet(01_w);
    p.contains_empty_word(true);
    presentation::add_rule(p, 000_w, 11_w);
    presentation::add_rule(p, 001_w, 10_w);

    word_graph_type wg;

    // Wrong alphabet size
    wg = make<WordGraph<node_type>>(3, {{1, 1, 1}, {2, 2, 2}, {2, 2, 2}});
    wg.number_of_active_nodes(3);
    REQUIRE(!sims::is_right_congruence(p, wg));

    // Incomplete
    wg = make<WordGraph<node_type>>(2, {{1, 1}, {1, UNDEFINED}});
    wg.number_of_active_nodes(2);
    REQUIRE(!sims::is_right_congruence(p, wg));

    // Incompatible
    wg = make<WordGraph<node_type>>(2, {{1, 1}, {1, 0}});
    wg.number_of_active_nodes(2);
    REQUIRE(!sims::is_right_congruence(p, wg));
    REQUIRE_THROWS_AS(sims::throw_if_not_right_congruence(p, wg),
                      LibsemigroupsException);

    // Works
    wg = make<WordGraph<node_type>>(4, {{1, 2}, {2, 2}, {3, 3}, {3, 3}});
    wg.number_of_active_nodes(4);
    REQUIRE(sims::is_right_congruence(p, wg));

    // Non maximal
    wg = make<WordGraph<node_type>>(2, {{1, 1}, {1, 0}});
    wg.number_of_active_nodes(2);
    REQUIRE(!sims::is_maximal_right_congruence(p, wg));
    wg = make<WordGraph<node_type>>(4, {{1, 2}, {2, 2}, {3, 3}, {3, 3}});
    wg.number_of_active_nodes(4);
    REQUIRE(!sims::is_maximal_right_congruence(p, wg));
    wg = make<WordGraph<node_type>>(1, {{0, 0}});
    wg.number_of_active_nodes(1);
    REQUIRE(!sims::is_maximal_right_congruence(p, wg));

    // Is maximal
    wg = make<WordGraph<node_type>>(2, {{1, 1}, {1, 1}});
    wg.number_of_active_nodes(2);
    REQUIRE(sims::is_maximal_right_congruence(p, wg));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "125",
                          "Two-sided congruence checking",
                          "[quick][low-index]") {
    Presentation<word_type> p;
    p.alphabet(01_w);
    p.contains_empty_word(true);
    presentation::add_rule(p, 000_w, 11_w);
    presentation::add_rule(p, 001_w, 10_w);

    word_graph_type wg;

    // Wrong alphabet size
    wg = make<WordGraph<node_type>>(3, {{1, 1, 1}, {2, 2, 2}, {2, 2, 2}});
    wg.number_of_active_nodes(3);
    REQUIRE(!sims::is_two_sided_congruence(p, wg));

    // Incomplete
    wg = make<WordGraph<node_type>>(2, {{1, 1}, {1, UNDEFINED}});
    wg.number_of_active_nodes(2);
    REQUIRE(!sims::is_two_sided_congruence(p, wg));

    // Incompatible
    wg = make<WordGraph<node_type>>(2, {{1, 1}, {1, 0}});
    wg.number_of_active_nodes(2);
    REQUIRE(!sims::is_two_sided_congruence(p, wg));
    REQUIRE_THROWS_AS(sims::throw_if_not_two_sided_congruence(p, wg),
                      LibsemigroupsException);

    // Not compatible with X_Gamma
    wg = make<WordGraph<node_type>>(4, {{1, 2}, {2, 2}, {3, 3}, {3, 3}});
    wg.number_of_active_nodes(4);
    REQUIRE(!sims::is_two_sided_congruence(p, wg));
    REQUIRE_THROWS_AS(sims::throw_if_not_two_sided_congruence(p, wg),
                      LibsemigroupsException);

    // Works
    wg = make<WordGraph<node_type>>(2, {{1, 1}, {1, 1}});
    wg.number_of_active_nodes(2);
    REQUIRE(sims::is_two_sided_congruence(p, wg));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "126",
                          "to_human_readable_repr test",
                          "[quick][low-index]") {
    Presentation<word_type> p;
    p.alphabet(01_w);
    p.contains_empty_word(true);
    presentation::add_rule(p, 000_w, 11_w);
    presentation::add_rule(p, 001_w, 10_w);

    SimsStats           sims_stats;
    Sims1               sims1;
    Sims2               sims2;
    RepOrc              rep_orc;
    MinimalRepOrc       minimal_rep_orc;
    SimsRefinerIdeals   sims_refiner_ideals;
    SimsRefinerFaithful sims_refiner_faithful;
    REQUIRE(to_human_readable_repr(sims_stats) == "<SimsStats object>");
    REQUIRE(to_human_readable_repr(sims_refiner_ideals)
            == fmt::format(
                "<SimsRefinerIdeals object over presentation {}>",
                to_human_readable_repr(sims_refiner_ideals.presentation())));
    sims_refiner_ideals.init(p);
    REQUIRE(to_human_readable_repr(sims_refiner_ideals)
            == fmt::format(
                "<SimsRefinerIdeals object over presentation {}>",
                to_human_readable_repr(sims_refiner_ideals.presentation())));
    REQUIRE(to_human_readable_repr(sims_refiner_faithful)
            == "<SimsRefinerFaithful object with 0 forbidden pairs>");
    sims_refiner_faithful.init({01_w, 10_w});
    REQUIRE(to_human_readable_repr(sims_refiner_faithful)
            == "<SimsRefinerFaithful object with 1 forbidden pair>");
    sims_refiner_faithful.init({01_w, 10_w, 100_w, 0_w, 011_w, 111_w});
    REQUIRE(to_human_readable_repr(sims_refiner_faithful)
            == "<SimsRefinerFaithful object with 3 forbidden pairs>");
    REQUIRE(to_human_readable_repr(sims1)
            == fmt::format("<Sims1 object over {} with 1 thread>",
                           to_human_readable_repr(sims1.presentation())));
    REQUIRE(to_human_readable_repr(sims2)
            == fmt::format("<Sims2 object over {} with 1 thread>",
                           to_human_readable_repr(sims2.presentation())));
    REQUIRE(to_human_readable_repr(rep_orc)
            == fmt::format("<RepOrc object over {} with node bounds [0, 0), "
                           "target size 0 and 1 thread>",
                           to_human_readable_repr(rep_orc.presentation())));
    REQUIRE(
        to_human_readable_repr(minimal_rep_orc)
        == fmt::format(
            "<MinimalRepOrc object over {} with target size 0 and 1 thread>",
            to_human_readable_repr(minimal_rep_orc.presentation())));
    sims1.presentation(p);
    sims2.presentation(p);
    rep_orc.presentation(p);
    minimal_rep_orc.presentation(p);
    REQUIRE(to_human_readable_repr(sims1)
            == fmt::format("<Sims1 object over {} with 1 thread>",
                           to_human_readable_repr(p)));
    REQUIRE(to_human_readable_repr(sims2)
            == fmt::format("<Sims2 object over {} with 1 thread>",
                           to_human_readable_repr(p)));
    REQUIRE(to_human_readable_repr(rep_orc)
            == fmt::format("<RepOrc object over {} with node bounds [0, 0), "
                           "target size 0 and 1 thread>",
                           to_human_readable_repr(p)));
    REQUIRE(
        to_human_readable_repr(minimal_rep_orc)
        == fmt::format(
            "<MinimalRepOrc object over {} with target size 0 and 1 thread>",
            to_human_readable_repr(p)));
    sims1.number_of_threads(std::thread::hardware_concurrency());
    sims2.number_of_threads(std::thread::hardware_concurrency());
    rep_orc.number_of_threads(std::thread::hardware_concurrency());
    minimal_rep_orc.number_of_threads(std::thread::hardware_concurrency());
    REQUIRE(to_human_readable_repr(sims1)
            == fmt::format("<Sims1 object over {} with {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(sims2)
            == fmt::format("<Sims2 object over {} with {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(rep_orc)
            == fmt::format("<RepOrc object over {} with node bounds [0, 0), "
                           "target size 0 and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(minimal_rep_orc)
            == fmt::format("<MinimalRepOrc object over {} with target size 0 "
                           "and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    sims::add_included_pair(sims1, 01_w, 10_w);
    sims::add_included_pair(sims1, 010_w, 101_w);
    sims::add_included_pair(sims2, 01_w, 10_w);
    sims::add_included_pair(sims2, 010_w, 101_w);
    sims::add_included_pair(rep_orc, 01_w, 10_w);
    sims::add_included_pair(rep_orc, 010_w, 101_w);
    sims::add_included_pair(minimal_rep_orc, 01_w, 10_w);
    sims::add_included_pair(minimal_rep_orc, 010_w, 101_w);
    REQUIRE(to_human_readable_repr(sims1)
            == fmt::format(
                "<Sims1 object over {} with 2 included pairs and {} threads>",
                to_human_readable_repr(p),
                std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(sims2)
            == fmt::format(
                "<Sims2 object over {} with 2 included pairs and {} threads>",
                to_human_readable_repr(p),
                std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(rep_orc)
            == fmt::format("<RepOrc object over {} with 2 included pairs, "
                           "node bounds [0, 0), "
                           "target size 0 and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(minimal_rep_orc)
            == fmt::format("<MinimalRepOrc object over {} with 2 included "
                           "pairs, target size 0 and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    sims1.clear_included_pairs();
    sims2.clear_included_pairs();
    rep_orc.clear_included_pairs();
    minimal_rep_orc.clear_included_pairs();
    sims::add_excluded_pair(sims1, 11_w, 10_w);
    sims::add_excluded_pair(sims2, 11_w, 10_w);
    sims::add_excluded_pair(rep_orc, 11_w, 10_w);
    sims::add_excluded_pair(minimal_rep_orc, 11_w, 10_w);
    REQUIRE(to_human_readable_repr(sims1)
            == fmt::format(
                "<Sims1 object over {} with 1 excluded pair and {} threads>",
                to_human_readable_repr(p),
                std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(sims2)
            == fmt::format(
                "<Sims2 object over {} with 1 excluded pair and {} threads>",
                to_human_readable_repr(p),
                std::thread::hardware_concurrency()));
    REQUIRE(
        to_human_readable_repr(rep_orc)
        == fmt::format(
            "<RepOrc object over {} with 1 excluded pair, node bounds [0, 0), "
            "target size 0 and {} threads>",
            to_human_readable_repr(p),
            std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(minimal_rep_orc)
            == fmt::format("<MinimalRepOrc object over {} with 1 excluded "
                           "pair, target size 0 and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    sims::add_included_pair(sims1, 01_w, 10_w);
    sims::add_included_pair(sims1, 010_w, 101_w);
    sims::add_included_pair(sims2, 01_w, 10_w);
    sims::add_included_pair(sims2, 010_w, 101_w);
    sims::add_included_pair(rep_orc, 01_w, 10_w);
    sims::add_included_pair(rep_orc, 010_w, 101_w);
    sims::add_included_pair(minimal_rep_orc, 01_w, 10_w);
    sims::add_included_pair(minimal_rep_orc, 010_w, 101_w);
    REQUIRE(
        to_human_readable_repr(sims1)
        == fmt::format("<Sims1 object over {} with 2 included and 1 excluded "
                       "pairs and {} threads>",
                       to_human_readable_repr(p),
                       std::thread::hardware_concurrency()));
    REQUIRE(
        to_human_readable_repr(sims2)
        == fmt::format("<Sims2 object over {} with 2 included and 1 excluded "
                       "pairs and {} threads>",
                       to_human_readable_repr(p),
                       std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(rep_orc)
            == fmt::format("<RepOrc object over {} with 2 included and 1 "
                           "excluded pairs, node bounds [0, 0), "
                           "target size 0 and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(minimal_rep_orc)
            == fmt::format(
                "<MinimalRepOrc object over {} with 2 included and 1 excluded "
                "pairs, target size 0 and {} threads>",
                to_human_readable_repr(p),
                std::thread::hardware_concurrency()));

    std::vector<word_type> forbid = {0_w, 01_w, 00_w, ""_w};
    SimsRefinerFaithful    pruno(forbid);
    sims1.add_pruner(pruno);
    sims2.add_pruner(pruno);
    rep_orc.add_pruner(pruno);
    minimal_rep_orc.add_pruner(pruno);
    REQUIRE(
        to_human_readable_repr(sims1)
        == fmt::format("<Sims1 object over {} with 2 included and 1 excluded "
                       "pairs, 1 pruner and {} threads>",
                       to_human_readable_repr(p),
                       std::thread::hardware_concurrency()));
    REQUIRE(
        to_human_readable_repr(sims2)
        == fmt::format("<Sims2 object over {} with 2 included and 1 excluded "
                       "pairs, 1 pruner and {} threads>",
                       to_human_readable_repr(p),
                       std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(rep_orc)
            == fmt::format("<RepOrc object over {} with 2 included and 1 "
                           "excluded pairs, node bounds [0, 0), "
                           "target size 0, 1 pruner and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(minimal_rep_orc)
            == fmt::format(
                "<MinimalRepOrc object over {} with 2 included and 1 excluded "
                "pairs, target size 0, 1 pruner and {} threads>",
                to_human_readable_repr(p),
                std::thread::hardware_concurrency()));
    SimsRefinerIdeals ideal_pruner(p);
    sims1.add_pruner(ideal_pruner);
    sims2.add_pruner(ideal_pruner);
    rep_orc.add_pruner(ideal_pruner);
    minimal_rep_orc.add_pruner(ideal_pruner);
    REQUIRE(
        to_human_readable_repr(sims1)
        == fmt::format("<Sims1 object over {} with 2 included and 1 excluded "
                       "pairs, 2 pruners and {} threads>",
                       to_human_readable_repr(p),
                       std::thread::hardware_concurrency()));
    REQUIRE(
        to_human_readable_repr(sims2)
        == fmt::format("<Sims2 object over {} with 2 included and 1 excluded "
                       "pairs, 2 pruners and {} threads>",
                       to_human_readable_repr(p),
                       std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(rep_orc)
            == fmt::format("<RepOrc object over {} with 2 included and 1 "
                           "excluded pairs, node bounds [0, 0), "
                           "target size 0, 2 pruners and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(minimal_rep_orc)
            == fmt::format(
                "<MinimalRepOrc object over {} with 2 included and 1 excluded "
                "pairs, target size 0, 2 pruners and {} threads>",
                to_human_readable_repr(p),
                std::thread::hardware_concurrency()));
    sims1.clear_excluded_pairs();
    sims2.clear_excluded_pairs();
    rep_orc.clear_excluded_pairs();
    minimal_rep_orc.clear_excluded_pairs();
    sims1.clear_pruners();
    sims2.clear_pruners();
    rep_orc.clear_pruners();
    minimal_rep_orc.clear_pruners();
    sims1.add_pruner(pruno);
    sims2.add_pruner(pruno);
    rep_orc.add_pruner(pruno);
    minimal_rep_orc.add_pruner(pruno);
    REQUIRE(to_human_readable_repr(sims1)
            == fmt::format("<Sims1 object over {} with 2 included pairs, 1 "
                           "pruner and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(sims2)
            == fmt::format("<Sims2 object over {} with 2 included pairs, 1 "
                           "pruner and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(rep_orc)
            == fmt::format("<RepOrc object over {} with 2 included pairs, "
                           "node bounds [0, 0), "
                           "target size 0, 1 pruner and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(minimal_rep_orc)
            == fmt::format("<MinimalRepOrc object over {} with 2 included "
                           "pairs, target size 0, 1 pruner and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    sims1.add_pruner(ideal_pruner);
    sims2.add_pruner(ideal_pruner);
    rep_orc.add_pruner(ideal_pruner);
    minimal_rep_orc.add_pruner(ideal_pruner);
    REQUIRE(to_human_readable_repr(sims1)
            == fmt::format("<Sims1 object over {} with 2 included pairs, 2 "
                           "pruners and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(sims2)
            == fmt::format("<Sims2 object over {} with 2 included pairs, 2 "
                           "pruners and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(rep_orc)
            == fmt::format("<RepOrc object over {} with 2 included pairs, "
                           "node bounds [0, 0), "
                           "target size 0, 2 pruners and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
    REQUIRE(to_human_readable_repr(minimal_rep_orc)
            == fmt::format("<MinimalRepOrc object over {} with 2 included "
                           "pairs, target size 0, 2 pruners and {} threads>",
                           to_human_readable_repr(p),
                           std::thread::hardware_concurrency()));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "127",
                          "symmetric_inverse_monoid(3)",
                          "[quick][low-index]") {
    auto p = symmetric_inverse_monoid(3, author::Gay);
    presentation::remove_duplicate_rules(p);
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);

    REQUIRE(presentation::length(p) == 32);
    REQUIRE(p.alphabet().size() == 3);
    REQUIRE(*presentation::shortest_rule(p) == 00_w);
    REQUIRE(*(presentation::shortest_rule(p) + 1) == ""_w);
    REQUIRE(presentation::longest_rule_length(p) == 8);

    Sims1 C;
    C.presentation(p);
    sims::add_excluded_pair(C, ""_w, 11_w);
    REQUIRE(C.number_of_threads(2).number_of_congruences(34) == 0);
    C.clear_excluded_pairs();

    auto rg = ReportGuard(true);
    REQUIRE(C.number_of_threads(2).number_of_congruences(34) == 274);

    sims::add_excluded_pair(C, ""_w, 11_w);
    C.clear_pruners();
    REQUIRE(C.number_of_threads(2).number_of_congruences(34) == 0);

    Sims2 S(p);
    S.presentation(p);
    sims::add_excluded_pair(S, ""_w, 11_w);
    REQUIRE(S.number_of_threads(2).number_of_congruences(34) == 0);
    S.clear_excluded_pairs();

    REQUIRE(S.number_of_threads(2).number_of_congruences(34) == 7);

    sims::add_excluded_pair(S, ""_w, 11_w);
    S.clear_pruners();
    REQUIRE(S.number_of_threads(2).number_of_congruences(34) == 0);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "128",
                          "SimsRefinerFaithful test",
                          "[quick][low-index]") {
    auto                    rg = ReportGuard(true);
    Presentation<word_type> p;

    p.alphabet({0, 1});
    p.contains_empty_word(true);
    presentation::add_rule(p, 000_w, 0_w);
    presentation::add_rule(p, 111_w, ""_w);
    presentation::add_rule(p, 011_w, 10_w);

    std::vector<word_type> forbid = {0_w, 01_w, 00_w, ""_w};
    SimsRefinerFaithful    pruno(forbid);

    Sims1 S;
    S.presentation(p);
    S.add_pruner(pruno);
    REQUIRE(S.number_of_threads(2).number_of_congruences(9)
            == 4);  // Verified with GAP

    auto it = S.cbegin(9);
    REQUIRE(*(it++)
            == make<WordGraph<node_type>>(
                9, {{1, 2}, {1, 3}, {4, 5}, {4, 4}, {3, 1}, {3, 0}}));
    REQUIRE(*(it++)
            == make<WordGraph<node_type>>(
                9, {{1, 2}, {3, 3}, {4, 5}, {1, 4}, {4, 1}, {3, 0}}));
    REQUIRE(*(it++)
            == make<WordGraph<node_type>>(
                9, {{1, 2}, {3, 4}, {3, 5}, {1, 1}, {4, 3}, {4, 0}}));
    REQUIRE(*(it++)
            == make<WordGraph<node_type>>(9,
                                          {{1, 2},
                                           {3, 4},
                                           {5, 6},
                                           {1, 7},
                                           {8, 5},
                                           {7, 1},
                                           {4, 0},
                                           {5, 8},
                                           {4, 3}}));
    REQUIRE(*(it++) == WordGraph<node_type>(0, 2));
    REQUIRE(*(it++) == WordGraph<node_type>(0, 2));
    REQUIRE(*(it++) == WordGraph<node_type>(0, 2));
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "129",
                          "Threading tests",
                          "[quick][low-index]") {
    auto                    rg = ReportGuard(true);
    Presentation<word_type> p;
    Sims1                   S;

    p.alphabet({0, 1, 2});
    p.contains_empty_word(true);
    S.presentation(p);

    S.number_of_threads(std::thread::hardware_concurrency() + 1);
    REQUIRE(S.number_of_threads() == std::thread::hardware_concurrency());
    S.idle_thread_restarts(1000);

    auto wg = S.number_of_threads(1).find_if(
        4, [](auto const& wg) { return wg.number_of_active_nodes() == 2; });
    REQUIRE(wg.number_of_active_nodes() == 2);
    wg = S.number_of_threads(1).find_if(3, [](auto const&) { return false; });
    REQUIRE(wg.number_of_active_nodes() == 0);
    sims::add_excluded_pair(S, {0, 0}, {0, 1});
    sims::add_included_pair(S, {0, 0}, {0, 1});
    REQUIRE(S.number_of_congruences(3) == 0);
  }

  LIBSEMIGROUPS_TEST_CASE("Sims1",
                          "130",
                          "MinimalRepOrc test",
                          "[quick][low-index]") {
    auto                    rg = ReportGuard(true);
    Presentation<word_type> p;

    p.alphabet({0, 1});
    p.contains_empty_word(true);
    presentation::add_rule(p, 000_w, 0_w);
    presentation::add_rule(p, 111_w, ""_w);
    presentation::add_rule(p, 011_w, 10_w);

    MinimalRepOrc orc;
    auto          d = orc.presentation(p)
                 .target_size(9)
                 .number_of_threads(std::thread::hardware_concurrency())
                 .word_graph();

    REQUIRE(d.number_of_nodes() == 6);
  }

  // TODO(0): test fails for now, uncomment when we get Sims<Word> to work as
  // per forthcoming changes to congruence interface and friends
  // LIBSEMIGROUPS_TEST_CASE("Sims1",
  //                         "131",
  //                         "fp semigroup containing given pairs, std::string",
  //                         "[quick][low-index]") {
  //   auto rg = ReportGuard(false);
  //
  //   Presentation<std::string> p;
  //   p.contains_empty_word(true);
  //   p.alphabet("ab");
  //   presentation::add_rule(p, "aaa", "a");
  //   presentation::add_rule(p, "bb", "b");
  //   presentation::add_rule(p, "abab", "a");
  //
  //   {
  //     Sims1 T;
  //     T.presentation(p);
  //     sims::add_included_pair(T, "abab", "a");
  //     REQUIRE(T.number_of_congruences(5) == 6);
  //   }
  //   {
  //     Sims1 T;
  //     presentation::reverse(p);
  //     T.presentation(p);
  //     sims::add_included_pair(T, "abab", "a");
  //     REQUIRE(T.included_pairs() == std::vector<word_type>({"abab", "a"}));
  //     REQUIRE(T.number_of_congruences(5) == 4);  // Verified with GAP
  //   }
  //   check_include(p, {"abab", "a"}, 5);
  //   presentation::reverse(p);
  //   check_include(p, {"abab", "a"}, 5);
  // }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "256",
                          "Partition monoid mfrc",
                          "[fail][low-index]") {
    // This doesn't fail it's just very extreme
    auto rg = ReportGuard(true);

    using words::operator+;

    Presentation<word_type> p;
    // plusses indicate that the result was reached but could not be verified.
    std::vector<size_t>  results = {1, 2, 7, 22, +84, +364, +1734, +8943};
    WordGraph<node_type> wg;
    Sims1                sims;
    // The forbid list depends on the presentation used by
    // fpsemigroup::partition_monoid
    // If this changes, then the test will no longer be correct.
    std::vector<word_type> forbid;
    SimsRefinerFaithful    pruno(forbid);
    for (size_t n = 2; n <= 7; ++n) {
      p = partition_monoid(n, author::Halverson + author::Ram);
      presentation::sort_rules(p);
      presentation::remove_duplicate_rules(p);

      sims.clear_pruners();
      forbid.clear();
      word_type alpha = {};
      word_type beta  = {};
      word_type q     = {n - 1};
      for (size_t i = 0; i < n; ++i) {
        if (i != n - 1) {
          alpha.push_back(i + 2 * n - 2);
        }
        beta.push_back(i + 2 * n - 2);
      }
      forbid = {alpha, beta, beta + q, beta, q + beta, beta};

      pruno.init(forbid);
      wg = sims.presentation(p)
               .number_of_threads(std::thread::hardware_concurrency())
               .add_pruner(pruno)
               .find_if(results[n], [](auto const&) { return true; });
      REQUIRE(wg.number_of_nodes() == results[n]);
      wg = sims.find_if(wg.number_of_nodes() - 1,
                        [](auto const&) { return true; });
      REQUIRE(wg.number_of_nodes() == 0);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "257",
                          "Temperley-Lieb monoid mfrc",
                          "[fail][low-index]") {
    // This doesn't fail it's just very extreme
    auto rg = ReportGuard(true);

    using words::operator+;

    Presentation<word_type> p;
    // plusses indicate that the result was reached but could not be verified.
    std::vector<size_t>    results = {1,
                                      1,
                                      2,
                                      4,
                                      7,
                                      10,
                                      20,
                                      29,
                                      +63,
                                      +91,
                                      +208,
                                      +298,
                                      +705,
                                      +1002,
                                      +2432,
                                      +3433,
                                      +8503};
    WordGraph<node_type>   wg;
    Sims1                  sims;
    std::vector<word_type> forbid;
    SimsRefinerFaithful    pruno(forbid);
    for (size_t n = 3; n <= 16; ++n) {
      p = temperley_lieb_monoid(n);
      presentation::sort_rules(p);
      presentation::remove_duplicate_rules(p);

      sims.clear_pruners();
      forbid.clear();
      word_type alpha = {};
      word_type beta  = {};
      word_type q     = {1};
      for (size_t i = 0; i < n - 1; i += 2) {
        if (i != n - 2) {
          alpha.push_back(i);
        }
        beta.push_back(i);
      }
      if (n % 2 == 0) {
        forbid = {alpha, beta, beta + q, beta, q + beta, beta};
      } else {
        forbid = {beta + q, beta, q + beta, beta};
      }

      pruno.init(forbid);
      wg = sims.presentation(p)
               .number_of_threads(std::thread::hardware_concurrency())
               .add_pruner(pruno)
               .find_if(results[n], [](auto const&) { return true; });
      REQUIRE(wg.number_of_nodes() == results[n]);
      wg = sims.find_if(wg.number_of_nodes() - 1,
                        [](auto const&) { return true; });
      REQUIRE(wg.number_of_nodes() == 0);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "258",
                          "Partial Brauer monoid mfrc",
                          "[fail][low-index]") {
    // This doesn't fail it's just very extreme
    auto rg = ReportGuard(true);

    using words::operator+;

    Presentation<word_type> p;
    // plusses indicate that the result was reached but could not be verified.
    std::vector<size_t>  results = {1, 2, 6, 14, +39, +117, +383, +1311, +7080};
    WordGraph<node_type> wg;
    Sims1                sims;
    std::vector<word_type> forbid;
    SimsRefinerFaithful    pruno(forbid);
    for (size_t n = 2; n <= 8; ++n) {
      p = fpsemigroup::partial_brauer_monoid(n);
      presentation::sort_rules(p);
      presentation::remove_duplicate_rules(p);

      sims.clear_pruners();
      forbid.clear();
      word_type alpha = {};
      word_type beta  = {};
      word_type q     = {n - 1};
      for (size_t i = 0; i < n; ++i) {
        if (i != n - 1) {
          alpha.push_back(i + 2 * n - 2);
        }
        beta.push_back(i + 2 * n - 2);
      }
      forbid = {alpha, beta, beta + q, beta, q + beta, beta};

      pruno.init(forbid);
      wg = sims.presentation(p)
               .number_of_threads(std::thread::hardware_concurrency())
               .add_pruner(pruno)
               .find_if(results[n], [](auto const&) { return true; });
      REQUIRE(wg.number_of_nodes() == results[n]);
      wg = sims.find_if(wg.number_of_nodes() - 1,
                        [](auto const&) { return true; });
      REQUIRE(wg.number_of_nodes() == 0);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "259",
                          "Motzkin monoid mfrc",
                          "[fail][low-index]") {
    // This doesn't fail it's just very extreme
    auto rg = ReportGuard(true);

    using words::operator+;

    Presentation<word_type> p;
    // plusses indicate that the result was reached but could not be verified.
    std::vector<size_t> results
        = {1, 2, 6, 13, +31, +77, +197, +513, +1354, +3611, +9714};
    WordGraph<node_type>   wg;
    Sims1                  sims;
    std::vector<word_type> forbid;
    SimsRefinerFaithful    pruno(forbid);
    for (size_t n = 3; n <= 10; ++n) {
      p = fpsemigroup::motzkin_monoid(n);
      presentation::sort_rules(p);
      presentation::remove_duplicate_rules(p);

      sims.clear_pruners();
      forbid.clear();
      word_type alpha = {};
      word_type beta  = {};
      word_type q     = {0};
      for (size_t i = 0; i < n - 1; ++i) {
        if (i != n - 2) {
          alpha.push_back(i + n - 1);
          alpha.push_back(i + n - 1);
        }
        beta.push_back(i + n - 1);
        beta.push_back(i + n - 1);
      }
      forbid = {alpha, beta, beta + q, beta, q + beta, beta};

      pruno.init(forbid);
      wg = sims.presentation(p)
               .number_of_threads(std::thread::hardware_concurrency())
               .add_pruner(pruno)
               .find_if(results[n], [](auto const&) { return true; });
      REQUIRE(wg.number_of_nodes() == results[n]);
      wg = sims.find_if(wg.number_of_nodes() - 1,
                        [](auto const&) { return true; });
      REQUIRE(wg.number_of_nodes() == 0);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Sims2",
                          "260",
                          "Brauer monoid mfrc",
                          "[fail][low-index]") {
    // This doesn't fail it's just very extreme
    auto rg = ReportGuard(true);

    using words::operator+;

    Presentation<word_type> p;
    // plusses indicate that the result was reached but could not be verified.
    std::vector<size_t>    results = {1, 1, 3, 7, 22, 46, +196, +523};
    WordGraph<node_type>   wg;
    Sims1                  sims;
    std::vector<word_type> forbid;
    SimsRefinerFaithful    pruno(forbid);
    for (size_t n = 2; n <= 7; ++n) {
      p = fpsemigroup::brauer_monoid(n);
      presentation::sort_rules(p);
      presentation::remove_duplicate_rules(p);

      sims.clear_pruners();
      forbid.clear();
      word_type alpha = {};
      word_type beta  = {};
      word_type q     = {n};
      word_type s     = {n - 2};
      for (size_t i = 0; i < n - 1; i += 2) {
        if (i != n - 2) {
          alpha.push_back(i + n - 1);
        }
        beta.push_back(i + n - 1);
      }
      if (n % 2 == 0) {
        forbid = {alpha + s, alpha, beta + q, beta, q + beta, beta};
      } else {
        forbid = {beta + q, beta, q + beta, beta};
      }

      pruno.init(forbid);
      wg = sims.presentation(p)
               .number_of_threads(std::thread::hardware_concurrency())
               .add_pruner(pruno)
               .find_if(results[n], [](auto const&) { return true; });
      REQUIRE(wg.number_of_active_nodes() == results[n]);
      wg = sims.find_if(wg.number_of_active_nodes() - 1,
                        [](auto const&) { return true; });
      REQUIRE(wg.number_of_active_nodes() == 0);
    }
  }

}  // namespace libsemigroups
