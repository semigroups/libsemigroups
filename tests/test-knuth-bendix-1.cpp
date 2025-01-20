// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2019-2023 James D. Mitchell
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

// This file is one of six that contains tests for the KnuthBendix classes. In
// a mostly vain attempt to speed up compilation the tests are split across 6
// files as follows:
//
// 1: contains quick tests for KnuthBendix created from rules and all commented
//    out tests.
//
// 2: contains more quick tests for KnuthBendix created from rules
//
// 3: contains yet more quick tests for KnuthBendix created from rules
//
// 4: contains standard and extreme test for KnuthBendix created from rules
//
// 5: contains tests for KnuthBendix created from FroidurePin instances
//
// 6: contains tests for KnuthBendix using word_type presentations

// TODO(later)
// * The other examples from Sims' book (Chap.s 5 and 6) which use
//   reduction orderings different from shortlex
// * Examples from MAF

#define CATCH_CONFIG_ENABLE_ALL_STRINGMAKERS

#include <algorithm>  // for fill
#include <chrono>     // for milliseconds
#include <cstddef>    // for size_t
#include <string>     // for basic_string, operator==
#include <utility>    // for move
#include <vector>     // for vector, operator==

#include "catch_amalgamated.hpp"  // for AssertionHandler, ope...
#include "test-main.hpp"          // for LIBSEMIGROUPS_TEMPLATE_TEST_CASE

#include "libsemigroups/constants.hpp"        // for operator==, operator!=
#include "libsemigroups/exception.hpp"        // for LibsemigroupsException
#include "libsemigroups/knuth-bendix.hpp"     // for KnuthBendix, normal_f...
#include "libsemigroups/obvinf.hpp"           // for is_obviously_infinite
#include "libsemigroups/paths.hpp"            // for Paths
#include "libsemigroups/presentation.hpp"     // for add_rule, Presentation
#include "libsemigroups/ranges.hpp"           // for equal
#include "libsemigroups/to-froidure-pin.hpp"  // for to_froidure_pin
#include "libsemigroups/types.hpp"            // for word_type
#include "libsemigroups/word-graph.hpp"       // for WordGraph

#include "libsemigroups/detail/report.hpp"  // for ReportGuard

namespace libsemigroups {
  using namespace rx;
  using literals::operator""_w;

  congruence_kind constexpr twosided = congruence_kind::twosided;

  namespace {
    // Generate the 'normal forms' defined by an arbitrary WordGraph.
    // If `wg` corresponds to the Gilman graph of some KnuthBendix instance,
    // then the words returned are the normal forms of that KnuthBendix
    // instance. Since the node labels returned by gilman_graph() are
    // implementation dependent, the below function can be used to check that
    // `gilman_graph()` returns something that generates the correct normal
    // forms.
    template <typename Rewriter, typename ReductionOrder, typename WordType>
    [[nodiscard]] inline auto
    normal_forms_from_word_graph(KnuthBendix<Rewriter, ReductionOrder>& kb,
                                 WordGraph<WordType>&                   wg) {
      Paths paths(wg);
      paths.source(0);
      if (!kb.presentation().contains_empty_word()) {
        paths.next();
      }
      return paths;
    }
  }  // namespace

  using RewriteTrie     = detail::RewriteTrie;
  using RewriteFromLeft = detail::RewriteFromLeft;

#define KNUTH_BENDIX_TYPES RewriteTrie, RewriteFromLeft

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "000",
                                   "confluent fp semigroup 1 (infinite)",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abc");
    p.rules = {"ab", "ba", "ac", "ca",  "aa", "a",  "ac", "a",  "ca", "a", "bb",
               "bb", "bc", "cb", "bbb", "b",  "bc", "b",  "cb", "b",  "a", "b"};

    KnuthBendix<TestType> kb(twosided, p);

    // kb.process_pending_rules();

    REQUIRE(kb.number_of_active_rules() == 4);
    REQUIRE(kb.confluent());
    REQUIRE(knuth_bendix::reduce(kb, "ca") == "a");
    REQUIRE(knuth_bendix::reduce(kb, "ac") == "a");
    REQUIRE(knuth_bendix::contains(kb, "ca", "a"));
    REQUIRE(knuth_bendix::contains(kb, "ac", "a"));
    REQUIRE(kb.number_of_classes() == POSITIVE_INFINITY);
    REQUIRE(is_obviously_infinite(kb));

    auto nf = knuth_bendix::normal_forms(kb).min(1).max(5);

    REQUIRE((nf | to_vector())
            == std::vector<std::string>({"a", "c", "cc", "ccc", "cccc"}));
    // REQUIRE(knuth_bendix::is_reduced(kb));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "001",
                                   "confluent fp semigroup 2 (infinite)",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("abc");
    presentation::add_rule_no_checks(p, "ac", "ca");
    presentation::add_rule_no_checks(p, "aa", "a");
    presentation::add_rule_no_checks(p, "ac", "a");
    presentation::add_rule_no_checks(p, "ca", "a");
    presentation::add_rule_no_checks(p, "bb", "bb");
    presentation::add_rule_no_checks(p, "bc", "cb");
    presentation::add_rule_no_checks(p, "bbb", "b");
    presentation::add_rule_no_checks(p, "bc", "b");
    presentation::add_rule_no_checks(p, "cb", "b");
    presentation::add_rule_no_checks(p, "a", "b");

    KnuthBendix<TestType> kb(twosided, p);

    // kb.process_pending_rules();

    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 4);
    REQUIRE(is_obviously_infinite(kb));
    auto nf = knuth_bendix::normal_forms(kb).min(1).max(5);

    REQUIRE((nf | to_vector())
            == std::vector<std::string>({"a", "c", "cc", "ccc", "cccc"}));
    REQUIRE(kb.number_of_classes() == POSITIVE_INFINITY);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "002",
                                   "confluent fp semigroup 3 (infinite)",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("012");
    presentation::add_rule_no_checks(p, "01", "10");
    presentation::add_rule_no_checks(p, "02", "20");
    presentation::add_rule_no_checks(p, "00", "0");
    presentation::add_rule_no_checks(p, "02", "0");
    presentation::add_rule_no_checks(p, "20", "0");
    presentation::add_rule_no_checks(p, "11", "11");
    presentation::add_rule_no_checks(p, "12", "21");
    presentation::add_rule_no_checks(p, "111", "1");
    presentation::add_rule_no_checks(p, "12", "1");
    presentation::add_rule_no_checks(p, "21", "1");
    presentation::add_rule_no_checks(p, "0", "1");

    KnuthBendix<TestType> kb(twosided, p);

    // kb.process_pending_rules();

    REQUIRE(kb.number_of_active_rules() == 4);
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 4);
    REQUIRE(kb.number_of_classes() == POSITIVE_INFINITY);

    auto nf = knuth_bendix::normal_forms(kb);
    REQUIRE((nf.min(1).max(2) | to_vector())
            == std::vector<std::string>({"0", "2"}));

    REQUIRE((nf.min(1).max(12) | to_vector())
            == std::vector<std::string>({"0",
                                         "2",
                                         "22",
                                         "222",
                                         "2222",
                                         "22222",
                                         "222222",
                                         "2222222",
                                         "22222222",
                                         "222222222",
                                         "2222222222",
                                         "22222222222"}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "003",
                                   "non-confluent example wikipedia",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("01");
    presentation::add_rule_no_checks(p, "000", "");
    presentation::add_rule_no_checks(p, "111", "");
    presentation::add_rule_no_checks(p, "010101", "");

    KnuthBendix<TestType> kb(twosided, p);
    REQUIRE(kb.presentation().alphabet() == "01");
    REQUIRE(!kb.confluent());
    kb.run();
    REQUIRE(kb.number_of_active_rules() == 4);
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_classes() == POSITIVE_INFINITY);

    auto nf = knuth_bendix::normal_forms(kb);

    REQUIRE((nf.min(0).max(5) | to_vector())
            == std::vector<std::string>(
                {"",     "0",    "1",    "00",   "01",   "10",   "11",
                 "001",  "010",  "011",  "100",  "101",  "110",  "0010",
                 "0011", "0100", "0101", "0110", "1001", "1011", "1101"}));
    REQUIRE((nf.min(0).max(10) | all_of([&kb](auto const& w) {
               return knuth_bendix::reduce(kb, w) == w;
             })));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "004",
                                   "Example 5.1 in Sims (infinite)",

                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("abcd");
    presentation::add_rule_no_checks(p, "ab", "");
    presentation::add_rule_no_checks(p, "ba", "");
    presentation::add_rule_no_checks(p, "cd", "");
    presentation::add_rule_no_checks(p, "dc", "");
    presentation::add_rule_no_checks(p, "ca", "ac");

    KnuthBendix<TestType> kb(twosided, p);

    REQUIRE(!kb.confluent());
    kb.run();
    REQUIRE(kb.number_of_active_rules() == 8);
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_classes() == POSITIVE_INFINITY);

    auto nf = knuth_bendix::normal_forms(kb);
    REQUIRE((nf.min(0).max(5) | to_vector())
            == std::vector<std::string>(
                {"",     "a",    "b",    "c",    "d",    "aa",   "ac",
                 "ad",   "bb",   "bc",   "bd",   "cc",   "dd",   "aaa",
                 "aac",  "aad",  "acc",  "add",  "bbb",  "bbc",  "bbd",
                 "bcc",  "bdd",  "ccc",  "ddd",  "aaaa", "aaac", "aaad",
                 "aacc", "aadd", "accc", "addd",  // codespell:ignore
                 "bbbb", "bbbc", "bbbd", "bbcc", "bbdd", "bccc", "bddd",
                 "cccc", "dddd"}));
    REQUIRE((nf.min(0).max(6) | all_of([&kb](auto const& w) {
               return knuth_bendix::reduce(kb, w) == w;
             })));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "005",
                                   "Example 5.1 in Sims (infinite) x 2",

                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("aAbB");
    presentation::add_inverse_rules(p, "AaBb");
    presentation::add_rule_no_checks(p, "ba", "ab");

    KnuthBendix<TestType> kb(twosided, p);

    REQUIRE(!kb.confluent());
    kb.run();
    REQUIRE(kb.number_of_active_rules() == 8);
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_classes() == POSITIVE_INFINITY);

    auto nf = knuth_bendix::normal_forms(kb);
    REQUIRE((nf.min(0).max(5) | to_vector())
            == std::vector<std::string>(
                {"",     "a",    "A",    "b",    "B",    "aa",   "ab",
                 "aB",   "AA",   "Ab",   "AB",   "bb",   "BB",   "aaa",
                 "aab",  "aaB",  "abb",  "aBB",  "AAA",  "AAb",  "AAB",
                 "Abb",  "ABB",  "bbb",  "BBB",  "aaaa", "aaab", "aaaB",
                 "aabb", "aaBB", "abbb", "aBBB", "AAAA", "AAAb", "AAAB",
                 "AAbb", "AABB", "Abbb", "ABBB", "bbbb", "BBBB"}));
    REQUIRE((nf.min(0).max(6) | all_of([&kb](auto const& w) {
               return knuth_bendix::reduce(kb, w) == w;
             })));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "006",
                                   "Example 5.3 in Sims",

                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("ab");
    presentation::add_rule_no_checks(p, "aa", "");
    presentation::add_rule_no_checks(p, "bbb", "");
    presentation::add_rule_no_checks(p, "ababab", "");

    KnuthBendix<TestType> kb(twosided, p);

    REQUIRE(!kb.confluent());
    kb.run();
    REQUIRE(kb.number_of_active_rules() == 6);
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_classes() == 12);

    auto nf = knuth_bendix::normal_forms(kb);
    REQUIRE(nf.count() == 12);

    REQUIRE((nf | to_vector())
            == std::vector<std::string>({"",
                                         "a",
                                         "b",
                                         "ab",
                                         "ba",
                                         "bb",
                                         "aba",
                                         "abb",
                                         "bab",
                                         "bba",
                                         "babb",
                                         "bbab"}));
    REQUIRE((nf.min(0).max(6) | all_of([&kb](auto const& w) {
               return knuth_bendix::reduce(kb, w) == w;
             })));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "007",
                                   "Example 5.4 in Sims",

                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("Bab");
    presentation::add_rule_no_checks(p, "aa", "");
    presentation::add_rule_no_checks(p, "bB", "");
    presentation::add_rule_no_checks(p, "bbb", "");
    presentation::add_rule_no_checks(p, "ababab", "");

    KnuthBendix<TestType> kb(twosided, p);

    REQUIRE(!kb.confluent());
    kb.run();
    REQUIRE(kb.number_of_active_rules() == 11);
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_classes() == 12);

    auto nf = knuth_bendix::normal_forms(kb).min(1).max(5);
    REQUIRE(nf.size_hint() == 11);
    REQUIRE((nf | to_vector())
            == std::vector<std::string>({"B",
                                         "a",
                                         "b",
                                         "Ba",
                                         "aB",
                                         "ab",
                                         "ba",
                                         "BaB",
                                         "Bab",
                                         "aBa",
                                         "baB"}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "008",
                                   "Example 6.4 in Sims",

                                   "[quick][knuth-bendix][no-valgrind]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("abc");
    p.contains_empty_word(true);

    presentation::add_rule_no_checks(p, "aa", "");
    presentation::add_rule_no_checks(p, "bc", "");
    presentation::add_rule_no_checks(p, "bbb", "");
    presentation::add_rule_no_checks(p, "ababababababab", "");
    presentation::add_rule_no_checks(p, "abacabacabacabac", "");

    KnuthBendix<TestType> kb(twosided, p);

    REQUIRE(!kb.confluent());
    REQUIRE(!is_obviously_infinite(kb));
    // REQUIRE(!kb.is_obviously_finite());
    kb.run();
    REQUIRE(kb.number_of_active_rules() == 40);
    REQUIRE(kb.confluent());
    REQUIRE(knuth_bendix::reduce(kb, "cc") == "b");
    REQUIRE(knuth_bendix::reduce(kb, "ccc") == "");
    REQUIRE(kb.number_of_classes() == 168);

    auto nf = knuth_bendix::normal_forms(kb).min(1).max(5);
    REQUIRE((nf | to_vector())
            == std::vector<std::string>(
                {"a",    "b",    "c",    "ab",   "ac",   "ba",   "ca",
                 "aba",  "aca",  "bab",  "bac",  "cab",  "cac",  "abab",
                 "abac", "acab", "acac", "baba", "baca", "caba", "caca"}));
    auto S = to_froidure_pin(kb);
    REQUIRE(S.size() == 168);
    REQUIRE(S.generator(2).string(kb) == "c");
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "009",
                                   "random example",

                                   "[quick][knuth-bendix][no-valgrind]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("012");

    presentation::add_rule_no_checks(p, "000", "2");
    presentation::add_rule_no_checks(p, "111", "2");
    presentation::add_rule_no_checks(p, "010101", "2");
    presentation::add_identity_rules(p, '2');

    KnuthBendix<TestType> kb(twosided, p);

    REQUIRE(!kb.confluent());
    kb.run();
    REQUIRE(kb.number_of_active_rules() == 9);
    REQUIRE(kb.confluent());

    auto& wg = kb.gilman_graph();
    REQUIRE(wg.number_of_nodes() == 9);
    REQUIRE(wg.number_of_edges() == 13);
    REQUIRE(!word_graph::is_acyclic(wg));

    auto fp = to_froidure_pin(kb);
    fp.enumerate(100);

    auto expected = froidure_pin::current_normal_forms(fp);

    Paths paths(wg);
    paths.source(0).min(1).max(fp.current_max_word_length() + 1);

    REQUIRE(equal(expected, paths));

    auto nf = knuth_bendix::normal_forms(kb).min(1).max(5);
    REQUIRE((nf | to_vector())
            == std::vector<std::string>(
                {"0",    "1",    "2",    "00",   "01",   "10",   "11",
                 "001",  "010",  "011",  "100",  "101",  "110",  "0010",
                 "0011", "0100", "0101", "0110", "1001", "1011", "1101"}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "010",
                                   "SL(2, 7) from Chap. 3, Prop. 1.5 in NR",
                                   "[quick][knuth-bendix][no-valgrind]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("abAB");
    p.contains_empty_word(true);

    presentation::add_rule_no_checks(p, "aaaaaaa", "");
    presentation::add_rule_no_checks(p, "bb", "ababab");
    presentation::add_rule_no_checks(p, "bb", "aaaabaaaabaaaabaaaab");
    presentation::add_rule_no_checks(p, "aA", "");
    presentation::add_rule_no_checks(p, "Aa", "");
    presentation::add_rule_no_checks(p, "bB", "");
    presentation::add_rule_no_checks(p, "Bb", "");

    KnuthBendix<TestType> kb(twosided, p);

    REQUIRE(!kb.confluent());

    kb.run();
    REQUIRE(kb.number_of_active_rules() == 152);
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_classes() == 336);

    // Test copy constructor
    auto T = to_froidure_pin(kb);
    auto S = froidure_pin::copy_closure(T, {T.generator(0)});

    REQUIRE(S.size() == 336);
    // 5 because S is generated as semigroup by 5 generators, while p is a
    // monoid presentation
    REQUIRE(S.number_of_generators() == 5);

    auto& wg = kb.gilman_graph();
    REQUIRE(wg.number_of_nodes() == 232);
    REQUIRE(wg.number_of_edges() == 265);
    REQUIRE(word_graph::is_acyclic(wg));
    Paths paths(wg);
    paths.source(0).min(0).max(13);
    REQUIRE(paths.count() == 336);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "011",
                                   "F(2, 5) - Chap. 9, Sec. 1 in NR",
                                   "[knuth-bendix][quick]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abcde");

    presentation::add_rule_no_checks(p, "ab", "c");
    presentation::add_rule_no_checks(p, "bc", "d");
    presentation::add_rule_no_checks(p, "cd", "e");
    presentation::add_rule_no_checks(p, "de", "a");
    presentation::add_rule_no_checks(p, "ea", "b");
    KnuthBendix<TestType> kb(twosided, p);

    REQUIRE(!kb.confluent());
    kb.run();
    REQUIRE(kb.number_of_active_rules() == 24);
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_classes() == 11);

    auto& wg = kb.gilman_graph();
    REQUIRE(wg.number_of_nodes() == 8);
    REQUIRE(wg.number_of_edges() == 11);
    REQUIRE(word_graph::is_acyclic(wg));
    Paths paths(wg);
    paths.source(0).min(0).max(5);
    REQUIRE(paths.count() == 12);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "012",
                                   "Reinis example 1",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abc");

    presentation::add_rule_no_checks(p, "a", "abb");
    presentation::add_rule_no_checks(p, "b", "baa");
    KnuthBendix<TestType> kb(twosided, p);

    REQUIRE(!kb.confluent());
    kb.run();
    REQUIRE(kb.number_of_active_rules() == 4);

    auto& wg = kb.gilman_graph();
    REQUIRE(wg.number_of_nodes() == 7);
    REQUIRE(wg.number_of_edges() == 17);
    REQUIRE(!word_graph::is_acyclic(wg));
    Paths paths(wg);
    paths.source(0).min(0).max(10);
    REQUIRE(paths.count() == 13'044);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "013",
                                   "redundant_rule (std::string)",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abc");
    presentation::add_rule(p, "a", "abb");
    presentation::add_rule(p, "b", "baa");
    presentation::add_rule(p, "c", "abbabababaaababababab");

    auto it = knuth_bendix::redundant_rule(p, std::chrono::milliseconds(100));
    REQUIRE(it == p.rules.cend());

    presentation::add_rule(p, "b", "baa");
    it = knuth_bendix::redundant_rule(p, std::chrono::milliseconds(100));
    REQUIRE(it != p.rules.cend());
    REQUIRE(*it == "b");
    REQUIRE(*(it + 1) == "baa");
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "014",
                                   "redundant_rule (word_type)",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    using literals::operator""_w;

    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(3);
    presentation::add_rule(p, 0_w, 011_w);
    presentation::add_rule(p, 1_w, 100_w);
    presentation::add_rule(p, 2_w, 011010101000101010101_w);

    auto it = knuth_bendix::redundant_rule(p, std::chrono::milliseconds(10));
    REQUIRE(it == p.rules.cend());

    presentation::add_rule(p, 1_w, 100_w);
    it = knuth_bendix::redundant_rule(p, std::chrono::milliseconds(10));
    REQUIRE(it != p.rules.cend());
    REQUIRE(*it == 1_w);
    REQUIRE(*(it + 1) == 100_w);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "015",
                                   "constructors/init for finished",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p1;
    p1.contains_empty_word(true);
    p1.alphabet("abcd");
    presentation::add_rule_no_checks(p1, "ab", "");
    presentation::add_rule_no_checks(p1, "ba", "");
    presentation::add_rule_no_checks(p1, "cd", "");
    presentation::add_rule_no_checks(p1, "dc", "");
    presentation::add_rule_no_checks(p1, "ca", "ac");

    Presentation<std::string> p2;
    p2.contains_empty_word(true);
    p2.alphabet("01");
    presentation::add_rule_no_checks(p2, "000", "");
    presentation::add_rule_no_checks(p2, "111", "");
    presentation::add_rule_no_checks(p2, "010101", "");

    KnuthBendix<TestType> kb1(twosided, p1);
    REQUIRE(!kb1.confluent());
    REQUIRE(!kb1.finished());
    kb1.run();
    REQUIRE(kb1.confluent());
    REQUIRE(knuth_bendix::reduce(kb1, "abababbdbcbdbabdbdb") == "bbbbbbddd");

    kb1.init(twosided, p2);
    REQUIRE(!kb1.confluent());
    REQUIRE(!kb1.finished());
    REQUIRE(kb1.presentation() == p2);
    kb1.run();
    REQUIRE(kb1.finished());
    REQUIRE(kb1.confluent());
    REQUIRE(kb1.confluent_known());

    kb1.init(twosided, p1);
    REQUIRE(!kb1.confluent());
    REQUIRE(!kb1.finished());
    REQUIRE(kb1.presentation() == p1);
    kb1.run();
    REQUIRE(kb1.finished());
    REQUIRE(kb1.confluent());
    REQUIRE(kb1.confluent_known());
    REQUIRE(knuth_bendix::reduce(kb1, "abababbdbcbdbabdbdb") == "bbbbbbddd");

    KnuthBendix<TestType> kb2(std::move(kb1));
    REQUIRE(kb2.confluent());
    REQUIRE(kb2.confluent_known());
    REQUIRE(kb2.finished());
    REQUIRE(knuth_bendix::reduce(kb2, "abababbdbcbdbabdbdb") == "bbbbbbddd");

    kb1 = std::move(kb2);
    REQUIRE(kb1.confluent());
    REQUIRE(kb1.confluent_known());
    REQUIRE(kb1.finished());
    REQUIRE(knuth_bendix::reduce(kb1, "abababbdbcbdbabdbdb") == "bbbbbbddd");

    kb1.init(twosided, std::move(p1));
    REQUIRE(!kb1.confluent());
    REQUIRE(!kb1.finished());
    kb1.run();
    REQUIRE(kb1.finished());
    REQUIRE(kb1.confluent());
    REQUIRE(kb1.confluent_known());
    REQUIRE(knuth_bendix::reduce(kb1, "abababbdbcbdbabdbdb") == "bbbbbbddd");
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "016",
                                   "constructors/init for partially run",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    using literals::operator""_w;

    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("abc");

    presentation::add_rule_no_checks(p, "aa", "");
    presentation::add_rule_no_checks(p, "bc", "");
    presentation::add_rule_no_checks(p, "bbb", "");
    presentation::add_rule_no_checks(p, "ababababababab", "");
    presentation::add_rule_no_checks(p, "abacabacabacabacabacabacabacabac", "");

    KnuthBendix<TestType> kb1(twosided, p);
    REQUIRE(!kb1.confluent());
    REQUIRE(!kb1.finished());
    kb1.run_for(std::chrono::milliseconds(10));
    REQUIRE(!kb1.confluent());
    REQUIRE(!kb1.finished());

    kb1.init(twosided, p);
    REQUIRE(!kb1.confluent());
    REQUIRE(!kb1.finished());
    REQUIRE(kb1.presentation() == p);
    kb1.run_for(std::chrono::milliseconds(10));
    REQUIRE(!kb1.confluent());
    REQUIRE(!kb1.finished());

    KnuthBendix<TestType> kb2(kb1);
    REQUIRE(!kb2.confluent());
    REQUIRE(!kb2.finished());
    REQUIRE(kb2.presentation() == p);
    REQUIRE(kb1.number_of_active_rules() == kb2.number_of_active_rules());
    kb2.run_for(std::chrono::milliseconds(10));
    REQUIRE(!kb2.confluent());
    REQUIRE(!kb2.finished());

    size_t const M = kb2.number_of_active_rules();
    kb1            = std::move(kb2);
    REQUIRE(kb1.number_of_active_rules() == M);
    REQUIRE(!kb1.finished());
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "017",
                                   "non-trivial classes",

                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abc");
    presentation::add_rule_no_checks(p, "ab", "ba");
    presentation::add_rule_no_checks(p, "ac", "ca");
    presentation::add_rule_no_checks(p, "aa", "a");
    presentation::add_rule_no_checks(p, "ac", "a");
    presentation::add_rule_no_checks(p, "ca", "a");
    presentation::add_rule_no_checks(p, "bc", "cb");
    presentation::add_rule_no_checks(p, "bbb", "b");
    presentation::add_rule_no_checks(p, "bc", "b");
    presentation::add_rule_no_checks(p, "cb", "b");

    KnuthBendix<TestType> kb1(twosided, p);

    presentation::add_rule_no_checks(p, "a", "b");

    KnuthBendix<TestType> kb2(twosided, p);

    REQUIRE(knuth_bendix::contains(kb2, "a", "b"));
    REQUIRE(knuth_bendix::contains(kb2, "a", "ba"));
    REQUIRE(knuth_bendix::contains(kb2, "a", "bb"));
    REQUIRE(knuth_bendix::contains(kb2, "a", "bab"));

    REQUIRE(knuth_bendix::non_trivial_classes(kb1, kb2)
            == std::vector<std::vector<std::string>>(
                {{"b", "ab", "bb", "abb", "a"}}));
    REQUIRE(knuth_bendix::non_trivial_classes<word_type>(kb1, kb2)
            == std::vector<std::vector<word_type>>(
                {{{98}, {97, 98}, {98, 98}, {97, 98, 98}, {97}}}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "018",
                                   "non-trivial classes x 2",

                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abc");
    presentation::add_rule_no_checks(p, "ab", "ba");
    presentation::add_rule_no_checks(p, "ac", "ca");
    presentation::add_rule_no_checks(p, "aa", "a");
    presentation::add_rule_no_checks(p, "ac", "a");
    presentation::add_rule_no_checks(p, "ca", "a");
    presentation::add_rule_no_checks(p, "bc", "cb");
    presentation::add_rule_no_checks(p, "bbb", "b");
    presentation::add_rule_no_checks(p, "bc", "b");
    presentation::add_rule_no_checks(p, "cb", "b");

    KnuthBendix<TestType> kb1(twosided, p);
    REQUIRE(kb1.number_of_classes() == POSITIVE_INFINITY);

    presentation::add_rule_no_checks(p, "b", "c");

    KnuthBendix<TestType> kb2(twosided, p);
    REQUIRE(kb2.number_of_classes() == 2);

    REQUIRE_THROWS_AS(knuth_bendix::non_trivial_classes(kb1, kb2),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "019",
                                   "non-trivial classes x 3",

                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abc");
    presentation::add_rule_no_checks(p, "ab", "ba");
    presentation::add_rule_no_checks(p, "ac", "ca");
    presentation::add_rule_no_checks(p, "aa", "a");
    presentation::add_rule_no_checks(p, "ac", "a");
    presentation::add_rule_no_checks(p, "ca", "a");
    presentation::add_rule_no_checks(p, "bc", "cb");
    presentation::add_rule_no_checks(p, "bbb", "b");
    presentation::add_rule_no_checks(p, "bc", "b");
    presentation::add_rule_no_checks(p, "cb", "b");

    KnuthBendix<TestType> kb1(twosided, p);

    presentation::add_rule_no_checks(p, "bb", "a");

    KnuthBendix<TestType> kb2(twosided, p);

    REQUIRE(knuth_bendix::non_trivial_classes(kb1, kb2)
            == std::vector<std::vector<std::string>>(
                {{"ab", "b"}, {"bb", "abb", "a"}}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "020",
                                   "non-trivial classes x 4",

                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(4);
    presentation::add_rule_no_checks(p, 01_w, 10_w);
    presentation::add_rule_no_checks(p, 02_w, 20_w);
    presentation::add_rule_no_checks(p, 00_w, 0_w);
    presentation::add_rule_no_checks(p, 02_w, 0_w);
    presentation::add_rule_no_checks(p, 20_w, 0_w);
    presentation::add_rule_no_checks(p, 12_w, 21_w);
    presentation::add_rule_no_checks(p, 111_w, 1_w);
    presentation::add_rule_no_checks(p, 12_w, 1_w);
    presentation::add_rule_no_checks(p, 21_w, 1_w);
    presentation::add_rule_no_checks(p, 03_w, 0_w);
    presentation::add_rule_no_checks(p, 30_w, 0_w);
    presentation::add_rule_no_checks(p, 13_w, 1_w);
    presentation::add_rule_no_checks(p, 31_w, 1_w);
    presentation::add_rule_no_checks(p, 23_w, 2_w);
    presentation::add_rule_no_checks(p, 32_w, 2_w);

    KnuthBendix<TestType> kb1(twosided, p);

    presentation::add_rule_no_checks(p, 0_w, 1_w);

    KnuthBendix<TestType> kb2(twosided, p);
    REQUIRE(knuth_bendix::non_trivial_classes(kb1, kb2)
            == std::vector<std::vector<std::string>>(
                {{{1}, {0, 1}, {1, 1}, {0, 1, 1}, {0}}}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "021",
                                   "non-triv. cong. on infinite fp semigp",
                                   "[quick][knuth-bendix][no-valgrind]",
                                   KNUTH_BENDIX_TYPES) {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(5);
    presentation::add_rule_no_checks(p, 01_w, 0_w);
    presentation::add_rule_no_checks(p, 10_w, 0_w);
    presentation::add_rule_no_checks(p, 02_w, 0_w);
    presentation::add_rule_no_checks(p, 20_w, 0_w);
    presentation::add_rule_no_checks(p, 03_w, 0_w);
    presentation::add_rule_no_checks(p, 30_w, 0_w);
    presentation::add_rule_no_checks(p, 00_w, 0_w);
    presentation::add_rule_no_checks(p, 11_w, 0_w);
    presentation::add_rule_no_checks(p, 22_w, 0_w);
    presentation::add_rule_no_checks(p, 33_w, 0_w);
    presentation::add_rule_no_checks(p, 12_w, 0_w);
    presentation::add_rule_no_checks(p, 21_w, 0_w);
    presentation::add_rule_no_checks(p, 13_w, 0_w);
    presentation::add_rule_no_checks(p, 31_w, 0_w);
    presentation::add_rule_no_checks(p, 23_w, 0_w);
    presentation::add_rule_no_checks(p, 32_w, 0_w);
    presentation::add_rule_no_checks(p, 40_w, 0_w);
    presentation::add_rule_no_checks(p, 41_w, 1_w);
    presentation::add_rule_no_checks(p, 42_w, 2_w);
    presentation::add_rule_no_checks(p, 43_w, 3_w);
    presentation::add_rule_no_checks(p, 04_w, 0_w);
    presentation::add_rule_no_checks(p, 14_w, 1_w);
    presentation::add_rule_no_checks(p, 24_w, 2_w);
    presentation::add_rule_no_checks(p, 34_w, 3_w);

    KnuthBendix<TestType> kb1(twosided, p);

    WordGraph test_wg1 = make<WordGraph<size_t>>(
        6,
        {{1, 2, 3, 4, 5},
         {},
         {},
         {},
         {},
         {UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED, 5}});
    REQUIRE(kb1.number_of_classes() == POSITIVE_INFINITY);

    REQUIRE(
        equal((knuth_bendix::normal_forms<word_type>(kb1) | rx::take(1000)),
              (normal_forms_from_word_graph(kb1, test_wg1) | rx::take(1000))));

    presentation::add_rule_no_checks(p, 1_w, 2_w);
    KnuthBendix<TestType> kb2(twosided, p);

    WordGraph test_wg2 = make<WordGraph<size_t>>(
        5,
        {{1, 2, UNDEFINED, 3, 4},
         {},
         {},
         {},
         {UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED, 4}});

    REQUIRE(kb1.number_of_classes() == POSITIVE_INFINITY);
    REQUIRE(
        equal((knuth_bendix::normal_forms<word_type>(kb2) | rx::take(1000)),
              (normal_forms_from_word_graph(kb2, test_wg2) | rx::take(1000))));

    REQUIRE(knuth_bendix::contains(kb2, 1_w, 2_w));

    auto ntc = knuth_bendix::non_trivial_classes(kb1, kb2);
    REQUIRE(ntc.size() == 1);
    REQUIRE(ntc[0].size() == 2);
    REQUIRE(ntc == decltype(ntc)({{{2}, {1}}}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "022",
                                   "non-triv. cong. on infinite fp semigroup",
                                   "[quick][kbp]",
                                   KNUTH_BENDIX_TYPES) {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(5);
    presentation::add_rule_no_checks(p, 01_w, 0_w);
    presentation::add_rule_no_checks(p, 10_w, 0_w);
    presentation::add_rule_no_checks(p, 02_w, 0_w);
    presentation::add_rule_no_checks(p, 20_w, 0_w);
    presentation::add_rule_no_checks(p, 03_w, 0_w);
    presentation::add_rule_no_checks(p, 30_w, 0_w);
    presentation::add_rule_no_checks(p, 00_w, 0_w);
    presentation::add_rule_no_checks(p, 11_w, 0_w);
    presentation::add_rule_no_checks(p, 22_w, 0_w);
    presentation::add_rule_no_checks(p, 33_w, 0_w);
    presentation::add_rule_no_checks(p, 12_w, 0_w);
    presentation::add_rule_no_checks(p, 21_w, 0_w);
    presentation::add_rule_no_checks(p, 13_w, 0_w);
    presentation::add_rule_no_checks(p, 31_w, 0_w);
    presentation::add_rule_no_checks(p, 23_w, 0_w);
    presentation::add_rule_no_checks(p, 32_w, 0_w);
    presentation::add_rule_no_checks(p, 40_w, 0_w);
    presentation::add_rule_no_checks(p, 41_w, 2_w);
    presentation::add_rule_no_checks(p, 42_w, 3_w);
    presentation::add_rule_no_checks(p, 43_w, 1_w);
    presentation::add_rule_no_checks(p, 04_w, 0_w);
    presentation::add_rule_no_checks(p, 14_w, 2_w);
    presentation::add_rule_no_checks(p, 24_w, 3_w);
    presentation::add_rule_no_checks(p, 34_w, 1_w);

    KnuthBendix<TestType> kb1(twosided, p);

    presentation::add_rule_no_checks(p, 2_w, 3_w);

    KnuthBendix<TestType> kb2(twosided, p);
    auto                  ntc = knuth_bendix::non_trivial_classes(kb1, kb2);
    REQUIRE(ntc.size() == 1);
    REQUIRE(ntc[0].size() == 3);
    REQUIRE(ntc == decltype(ntc)({{{2}, {3}, {1}}}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "023",
                                   "triv. cong. on finite fp semigp",
                                   "[quick][kbp]",
                                   KNUTH_BENDIX_TYPES) {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(2);
    presentation::add_rule_no_checks(p, 001_w, 00_w);
    presentation::add_rule_no_checks(p, 0000_w, 00_w);
    presentation::add_rule_no_checks(p, 0110_w, 00_w);
    presentation::add_rule_no_checks(p, 0111_w, 000_w);
    presentation::add_rule_no_checks(p, 1110_w, 110_w);
    presentation::add_rule_no_checks(p, 1111_w, 111_w);
    presentation::add_rule_no_checks(p, 01000_w, 0101_w);
    presentation::add_rule_no_checks(p, 01010_w, 0100_w);
    presentation::add_rule_no_checks(p, 01011_w, 0101_w);

    KnuthBendix<TestType> kb1(twosided, p);
    KnuthBendix<TestType> kb2(twosided, p);

    REQUIRE(!p.contains_empty_word());
    REQUIRE(kb1.number_of_classes() == 27);
    REQUIRE(kb2.number_of_classes() == 27);
    auto ntc = knuth_bendix::non_trivial_classes(kb1, kb2);
    REQUIRE(ntc.empty());
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "024",
                                   "universal cong. on finite fp semigroup",
                                   "[quick][kbp]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<word_type> p;
    p.alphabet(2);
    presentation::add_rule_no_checks(p, 001_w, 00_w);
    presentation::add_rule_no_checks(p, 0000_w, 00_w);
    presentation::add_rule_no_checks(p, 0110_w, 00_w);
    presentation::add_rule_no_checks(p, 0111_w, 000_w);
    presentation::add_rule_no_checks(p, 1110_w, 110_w);
    presentation::add_rule_no_checks(p, 1111_w, 111_w);
    presentation::add_rule_no_checks(p, 01000_w, 0101_w);
    presentation::add_rule_no_checks(p, 01010_w, 0100_w);
    presentation::add_rule_no_checks(p, 01011_w, 0101_w);

    KnuthBendix<TestType> kb1(twosided, p);

    presentation::add_rule_no_checks(p, 0_w, 1_w);
    presentation::add_rule_no_checks(p, 00_w, 0_w);

    KnuthBendix<TestType> kb2(twosided, p);

    REQUIRE(kb2.number_of_classes() == 1);

    auto ntc = knuth_bendix::non_trivial_classes(kb1, kb2);

    REQUIRE(ntc.size() == 1);
    REQUIRE(ntc[0].size() == 27);
    std::vector<std::string> expected = {{0},
                                         {1},
                                         {0, 0},
                                         {0, 1},
                                         {1, 0},
                                         {1, 1},
                                         {0, 0, 0},
                                         {1, 0, 0},
                                         {0, 1, 0},
                                         {1, 0, 1},
                                         {0, 1, 1},
                                         {1, 1, 0},
                                         {1, 1, 1},
                                         {1, 0, 0, 0},
                                         {0, 1, 0, 0},
                                         {1, 1, 0, 0},
                                         {1, 0, 1, 0},
                                         {0, 1, 0, 1},
                                         {1, 1, 0, 1},
                                         {1, 0, 1, 1},
                                         {1, 1, 0, 0, 0},
                                         {1, 0, 1, 0, 0},
                                         {1, 1, 0, 1, 0},
                                         {1, 0, 1, 0, 1},
                                         {1, 1, 0, 1, 1},
                                         {1, 1, 0, 1, 0, 0},
                                         {1, 1, 0, 1, 0, 1}};
    std::sort(expected.begin(), expected.end());
    std::sort(ntc[0].begin(), ntc[0].end());
    REQUIRE(ntc[0] == expected);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "025",
                                   "finite fp semigroup, size 16",

                                   "[quick][kbp]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<word_type> p;
    p.alphabet(11);
    presentation::add_rule_no_checks(p, {2}, {1});
    presentation::add_rule_no_checks(p, {4}, {3});
    presentation::add_rule_no_checks(p, {5}, {0});
    presentation::add_rule_no_checks(p, {6}, {3});
    presentation::add_rule_no_checks(p, {7}, {1});
    presentation::add_rule_no_checks(p, {8}, {3});
    presentation::add_rule_no_checks(p, {9}, {3});
    presentation::add_rule_no_checks(p, {10}, {0});
    presentation::add_rule_no_checks(p, {0, 2}, {0, 1});
    presentation::add_rule_no_checks(p, {0, 4}, {0, 3});
    presentation::add_rule_no_checks(p, {0, 5}, {0, 0});
    presentation::add_rule_no_checks(p, {0, 6}, {0, 3});
    presentation::add_rule_no_checks(p, {0, 7}, {0, 1});
    presentation::add_rule_no_checks(p, {0, 8}, {0, 3});
    presentation::add_rule_no_checks(p, {0, 9}, {0, 3});

    presentation::add_rule_no_checks(p, {0, 10}, {0, 0});
    presentation::add_rule_no_checks(p, {1, 1}, {1});
    presentation::add_rule_no_checks(p, {1, 2}, {1});
    presentation::add_rule_no_checks(p, {1, 4}, {1, 3});
    presentation::add_rule_no_checks(p, {1, 5}, {1, 0});
    presentation::add_rule_no_checks(p, {1, 6}, {1, 3});
    presentation::add_rule_no_checks(p, {1, 7}, {1});
    presentation::add_rule_no_checks(p, {1, 8}, {1, 3});
    presentation::add_rule_no_checks(p, {1, 9}, {1, 3});
    presentation::add_rule_no_checks(p, {1, 10}, {1, 0});
    presentation::add_rule_no_checks(p, {3, 1}, {3});
    presentation::add_rule_no_checks(p, {3, 2}, {3});
    presentation::add_rule_no_checks(p, {3, 3}, {3});
    presentation::add_rule_no_checks(p, {3, 4}, {3});
    presentation::add_rule_no_checks(p, {3, 5}, {3, 0});
    presentation::add_rule_no_checks(p, {3, 6}, {3});
    presentation::add_rule_no_checks(p, {3, 7}, {3});
    presentation::add_rule_no_checks(p, {3, 8}, {3});
    presentation::add_rule_no_checks(p, {3, 9}, {3});
    presentation::add_rule_no_checks(p, {3, 10}, {3, 0});
    presentation::add_rule_no_checks(p, {0, 0, 0}, {0});
    presentation::add_rule_no_checks(p, {0, 0, 1}, {1});
    presentation::add_rule_no_checks(p, {0, 0, 3}, {3});
    presentation::add_rule_no_checks(p, {0, 1, 3}, {1, 3});
    presentation::add_rule_no_checks(p, {1, 0, 0}, {1});
    presentation::add_rule_no_checks(p, {1, 0, 3}, {0, 3});
    presentation::add_rule_no_checks(p, {3, 0, 0}, {3});
    presentation::add_rule_no_checks(p, {0, 1, 0, 1}, {1, 0, 1});
    presentation::add_rule_no_checks(p, {0, 3, 0, 3}, {3, 0, 3});
    presentation::add_rule_no_checks(p, {1, 0, 1, 0}, {1, 0, 1});
    presentation::add_rule_no_checks(p, {1, 3, 0, 1}, {1, 0, 1});
    presentation::add_rule_no_checks(p, {1, 3, 0, 3}, {3, 0, 3});
    presentation::add_rule_no_checks(p, {3, 0, 1, 0}, {3, 0, 1});
    presentation::add_rule_no_checks(p, {3, 0, 3, 0}, {3, 0, 3});

    KnuthBendix<TestType> kb1(twosided, p);
    REQUIRE(kb1.gilman_graph().number_of_nodes() == 16);

    WordGraph test_wg1
        = make<WordGraph<size_t>>(16,
                                  {{3,
                                    1,
                                    UNDEFINED,
                                    2,
                                    UNDEFINED,
                                    UNDEFINED,
                                    UNDEFINED,
                                    UNDEFINED,
                                    UNDEFINED,
                                    UNDEFINED,
                                    UNDEFINED},
                                   {6, UNDEFINED, UNDEFINED, 12},
                                   {7, UNDEFINED},
                                   {4, 5, UNDEFINED, 9},
                                   {},
                                   {8},
                                   {UNDEFINED, 11},
                                   {UNDEFINED, 14, UNDEFINED, 15},
                                   {},
                                   {10},
                                   {UNDEFINED, 14},
                                   {},
                                   {13},
                                   {UNDEFINED}});
    REQUIRE(equal(knuth_bendix::normal_forms<word_type>(kb1),
                  normal_forms_from_word_graph(kb1, test_wg1)));

    presentation::add_rule_no_checks(p, {1}, {3});
    KnuthBendix<TestType> kb2(twosided, p);

    WordGraph test_wg2 = make<WordGraph<size_t>>(4,
                                                 {{2,
                                                   1,
                                                   UNDEFINED,
                                                   UNDEFINED,
                                                   UNDEFINED,
                                                   UNDEFINED,
                                                   UNDEFINED,
                                                   UNDEFINED,
                                                   UNDEFINED,
                                                   UNDEFINED,
                                                   UNDEFINED},
                                                  {},
                                                  {3}});

    REQUIRE(equal(knuth_bendix::normal_forms<word_type>(kb2),
                  normal_forms_from_word_graph(kb2, test_wg2)));

    auto ntc = knuth_bendix::non_trivial_classes(kb1, kb2);

    std::vector<std::string> expected = {{1},
                                         {3},
                                         {0, 1},
                                         {0, 3},
                                         {1, 0},
                                         {3, 0},
                                         {1, 3},
                                         {0, 1, 0},
                                         {0, 3, 0},
                                         {1, 0, 1},
                                         {3, 0, 1},
                                         {3, 0, 3},
                                         {1, 3, 0},
                                         {0, 3, 0, 1}};
    std::sort(expected.begin(), expected.end());
    std::sort(ntc[0].begin(), ntc[0].end());
    REQUIRE(ntc[0] == expected);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "026",
                                   "non_trivial_classes exceptions",
                                   "[quick][kbp]",
                                   KNUTH_BENDIX_TYPES) {
    Presentation<word_type> p;
    p.alphabet(1);
    KnuthBendix<TestType> kbp(twosided, p);

    {
      Presentation<word_type> q;
      q.alphabet(2);
      KnuthBendix<TestType> kbq(twosided, q);
      REQUIRE_THROWS_AS(knuth_bendix::non_trivial_classes(kbp, kbq),
                        LibsemigroupsException);
      REQUIRE(kbq.number_of_inactive_rules() == 0);
    }
    {
      presentation::add_rule_no_checks(p, 0000_w, 00_w);
      kbp.init(twosided, p);

      Presentation<word_type> q;
      q.alphabet(1);
      presentation::add_rule_no_checks(q, 00_w, 0_w);

      KnuthBendix<TestType> kbq(twosided, q);
      REQUIRE_THROWS_AS(knuth_bendix::non_trivial_classes(kbq, kbp),
                        LibsemigroupsException);
    }
  }

  ////////////////////////////////////////////////////////////////////////
  // Commented out test cases
  ////////////////////////////////////////////////////////////////////////

  // This example verifies the nilpotence of the group using the Sims
  // algorithm. The original presentation was <a,b| [b,a,b], [b,a,a,a,a],
  // [b,a,a,a,b,a,a] >. (where [] mean left-normed commutators). The
  // presentation here was derived by first applying the NQA to find the
  // maximal nilpotent quotient, and then introducing new generators for
  // the PCP generators. It is essential for success that reasonably low
  // values of the maxstoredlen parameter are given.
  // LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
  //           //     "013",
  //     "(from kbmag/standalone/kb_data/verifynilp)",
  //     "000","[quick][knuth-bendix]"[kbmag][recursive],
  //  KNUTH_BENDIX_TYPES){}
  //   KnuthBendix<TestType> kb(new RECURSIVE(), "hHgGfFyYdDcCbBaA");
  //   presentation::add_rule_no_checks(p, "BAba", "c");
  //   presentation::add_rule_no_checks(p, "CAca", "d");
  //   presentation::add_rule_no_checks(p, "DAda", "y");
  //   presentation::add_rule_no_checks(p, "YByb", "f");
  //   presentation::add_rule_no_checks(p, "FAfa", "g");
  //   presentation::add_rule_no_checks(p, "ga", "ag");
  //   presentation::add_rule_no_checks(p, "GBgb", "h");
  //   presentation::add_rule_no_checks(p, "cb", "bc");
  //   presentation::add_rule_no_checks(p, "ya", "ay");
  //   auto rg = ReportGuard(false);
  //
  //   REQUIRE(kb.confluent());
  //
  //   kb.run();
  //   REQUIRE(kb.confluent());
  //   REQUIRE(kb.number_of_active_rules() == 9);
  //
  //   REQUIRE(knuth_bendix::contains(kb, "BAba", "c"));
  //   REQUIRE(knuth_bendix::contains(kb, "CAca", "d"));
  //   REQUIRE(knuth_bendix::contains(kb, "DAda", "y"));
  //   REQUIRE(knuth_bendix::contains(kb, "YByb", "f"));
  //   REQUIRE(knuth_bendix::contains(kb, "FAfa", "g"));
  //   REQUIRE(knuth_bendix::contains(kb, "ga", "ag"));
  //   REQUIRE(knuth_bendix::contains(kb, "GBgb", "h"));
  //   REQUIRE(knuth_bendix::contains(kb, "cb", "bc"));
  //   REQUIRE(knuth_bendix::contains(kb, "ya", "ay"));
  // REQUIRE(kb.active_rules() == std::vector<std::pair<std::string,
  // std::string>>({}));
  // }

  // TODO(later): temporarily commented out to because of change to
  // FpSemigroupInterface that forbids adding rules after started(), and
  // because the copy constructors for KnuthBendix<TestType> et al. don't
  // currently work LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
  //                         "(cong) finite transformation semigroup "
  //                         "congruence (21 classes)",
  // "000","[quick][congruence][knuth-bendix]"[cong], KNUTH_BENDIX_TYPES){
  //   auto rg      = ReportGuard(false);
  //   using Transf = LeastTransf<5>;
  //   FroidurePin<Transf> S({Transf({1, 3, 4, 2, 3}), Transf({3, 2, 1, 3,
  //   3})});

  //   REQUIRE(S.size() == 88);
  //   REQUIRE(S.number_of_rules() == 18);

  //   KnuthBendix<TestType> kb(twosided, S);
  //   auto&       P = kb.quotient_froidure_pin();
  //   REQUIRE(P.size() == 88);
  //   kb.add_generating_pair(S.factorisation(Transf({3, 4, 4, 4, 4})),
  //               S.factorisation(Transf({3, 1, 3, 3, 3})));
  //   // P is now invalid, it's a reference to something that was deleted
  //   in
  //   // kb.

  //   REQUIRE(kb.number_of_classes() == 21);
  //   REQUIRE(kb.number_of_classes() == 21);
  //   auto& Q = kb.quotient_froidure_pin();  // quotient

  //   REQUIRE(Q.size() == 21);
  //   REQUIRE(Q.number_of_idempotents() == 3);

  //   std::vector<word_type>
  //   v(static_cast<FroidurePin<detail::KBE>&>(Q).cbegin(),
  // static_cast<FroidurePin<detail::KBE>&>(Q).cend());
  //   REQUIRE(v
  //           == std::vector<word_type>({{0},
  //                                      {1},
  //                                      {0, 0},
  //                                      {0, 1},
  //                                      {1, 0},
  //                                      {1, 1},
  //                                      {0, 0, 0},
  //                                      {0, 0, 1},
  //                                      {0, 1, 0},
  //                                      {0, 1, 1},
  //                                      {1, 0, 0},
  //                                      {1, 1, 0},
  //                                      {0, 0, 0, 0},
  //                                      {0, 1, 0, 0},
  //                                      {0, 1, 1, 0},
  //                                      {1, 0, 0, 0},
  //                                      {1, 1, 0, 0},
  //                                      {0, 1, 0, 0, 0},
  //                                      {0, 1, 1, 0, 0},
  //                                      {1, 1, 0, 0, 0},
  //                                      {0, 1, 1, 0, 0, 0}}));

  //   REQUIRE(
  //       kb.word_to_class_index(S.factorisation(Transf({1, 3, 1, 3,
  //       3})))
  //       == kb.word_to_class_index(S.factorisation(Transf({4, 2, 4, 4,
  //       2}))));

  //   REQUIRE(kb.number_of_non_trivial_classes() == 1);
  //   REQUIRE(kb.number_of_generators() == 2);
  //   REQUIRE(kb.cbegin_ntc()->size() == 68);
  // }

  //  A nonhopfian group
  // LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
  //                         "(from kbmag/standalone/kb_data/nonhopf)",
  // "027","[quick][knuth-bendix]"[kbmag][recursive], KNUTH_BENDIX_TYPES){
  //   KnuthBendix<TestType> kb(new RECURSIVE(), "aAbB");
  //   presentation::add_rule_no_checks(p, "Baab", "aaa");
  //   auto rg = ReportGuard(false);

  //   REQUIRE(kb.confluent());

  //   kb.run();
  //   REQUIRE(kb.confluent());
  //   REQUIRE(kb.number_of_active_rules() == 1);

  //   REQUIRE(knuth_bendix::contains(kb, "Baab", "aaa"));
  //   REQUIRE(kb.active_rules()
  //           == std::vector<rule_type>({}));
  // }

  // LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
  //                         "(from kbmag/standalone/kb_data/freenilpc3)",
  // "028","[quick][knuth-bendix]"[kbmag][recursive], KNUTH_BENDIX_TYPES){
  //   KnuthBendix<TestType> kb(new RECURSIVE(), "yYdDcCbBaA");
  //   presentation::add_rule_no_checks(p, "BAba", "c");
  //   presentation::add_rule_no_checks(p, "CAca", "d");
  //   presentation::add_rule_no_checks(p, "CBcb", "y");
  //   presentation::add_rule_no_checks(p, "da", "wg");
  //   presentation::add_rule_no_checks(p, "ya", "ay");
  //   presentation::add_rule_no_checks(p, "db", "bd");
  //   presentation::add_rule_no_checks(p, "yb", "by");
  //   auto rg = ReportGuard(false);

  //   REQUIRE(kb.confluent());

  //   kb.run();
  //   REQUIRE(kb.confluent());
  //   REQUIRE(kb.number_of_active_rules() == 7);

  //   REQUIRE(knuth_bendix::contains(kb, "BAba", "c"));
  //   REQUIRE(knuth_bendix::contains(kb, "CAca", "d"));
  //   REQUIRE(knuth_bendix::contains(kb, "CBcb", "y"));
  //   REQUIRE(knuth_bendix::contains(kb, "da", "wg"));
  //   REQUIRE(knuth_bendix::contains(kb, "ya", "ay"));
  //   REQUIRE(knuth_bendix::contains(kb, "db", "bd"));
  //   REQUIRE(knuth_bendix::contains(kb, "yb", "by"));
  //   REQUIRE(kb.active_rules()
  //           == std::vector<rule_type>({}));
  // }

  // TODO(later): temporarily commented out to because of change to
  // FpSemigroupInterface that forbids adding rules after started(), and
  // because the copy constructors for KnuthBendix<TestType> et al. don't
  // currently work LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
  //                         "add_rule after knuth_bendix",
  //                         "029","[quick][knuth-bendix]",
  //  KNUTH_BENDIX_TYPES){
  //   auto        rg = ReportGuard(false);
  //   KnuthBendix<TestType> kb;
  //   Presentation<std::string> p;
  // //p.alphabet("Bab");
  //   presentation::add_rule_no_checks(p, "aa", "");
  //   presentation::add_rule_no_checks(p, "bB", "");
  //   presentation::add_rule_no_checks(p, "bbb", "");
  //   presentation::add_rule_no_checks(p, "ababab", "");

  //   REQUIRE(!kb.confluent());
  //   kb.run_for(FOREVER);
  //   REQUIRE(kb.finished());
  //   // The next line tests what happens when run_for is called when
  //   finished. kb.run_for(FOREVER); REQUIRE(kb.number_of_active_rules()
  //   == 11); REQUIRE(kb.confluent()); REQUIRE(kb.size() == 12);

  //   REQUIRE(knuth_bendix::contains(kb, "aa", ""));
  //   REQUIRE(!knuth_bendix::contains(kb, "a", "b"));

  //   KnuthBendix<TestType> kb2(&kb);
  //   REQUIRE(kb2.number_of_active_rules() == 11);
  //   kb2.add_rule("a", "b");
  //   REQUIRE(kb2.number_of_rules() == 5);
  //   // Adding a rule does not change the number of active rules until
  //   *after*
  //   // kb.run() is called again.
  //   REQUIRE(kb2.number_of_active_rules() == 11);

  //   using rules_type = std::vector<rule_type>;

  //   REQUIRE(rules_type(kb2.cbegin_rules(), kb2.cend_rules())
  //           == rules_type({{"aa", ""},
  //                          {"bB", ""},
  //                          {"bbb", ""},
  //                          {"ababab", ""},
  //                          {"a", "b"}}));

  //   REQUIRE(!kb2.confluent());
  //   REQUIRE(kb2.size() == 1);
  //   REQUIRE(kb2.confluent());
  //   REQUIRE(kb2.number_of_active_rules() == 3);
  //   REQUIRE(kb2.active_rules() == rules_type({{"B", ""}, {"a", ""},
  //   {"b", "a"}}));
  // }

  // Free nilpotent group of rank 2 and class 2
  // LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
  //                         "(from kbmag/standalone/kb_data/nilp2)",
  // "030","[quick][knuth-bendix]"[kbmag][recursive], KNUTH_BENDIX_TYPES){
  //   KnuthBendix<TestType> kb(new RECURSIVE(), "cCbBaA");
  //   presentation::add_rule_no_checks(p, "ba", "abc");
  //   presentation::add_rule_no_checks(p, "ca", "ac");
  //   presentation::add_rule_no_checks(p, "cb", "bc");
  //   auto rg = ReportGuard(false);
  //
  //   REQUIRE(kb.confluent());
  //
  //   kb.run();
  //   REQUIRE(kb.confluent());
  //
  //   REQUIRE(kb.number_of_active_rules() == 3);
  // }

  // monoid presentation of F(2,7) - should produce a monoid of length 30
  // which is the same as the group, together with the empty word. This
  // is a very difficult calculation indeed, however.
  //
  // KBMAG does not terminate when SHORTLEX order is used.
  // LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
  //                         "(from kbmag/standalone/kb_data/f27monoid)",
  //                         "[fail][knuth-bendix]"[kbmag][recursive],
  //  KNUTH_BENDIX_TYPES){
  //   KnuthBendix<TestType> kb(new RECURSIVE(), "abcdefg");
  //   presentation::add_rule_no_checks(p, "ab", "c");
  //   presentation::add_rule_no_checks(p, "bc", "d");
  //   presentation::add_rule_no_checks(p, "cd", "e");
  //   presentation::add_rule_no_checks(p, "de", "f");
  //   presentation::add_rule_no_checks(p, "ef", "g");
  //   presentation::add_rule_no_checks(p, "fg", "a");
  //   presentation::add_rule_no_checks(p, "ga", "b");
  //   auto rg = ReportGuard(false);

  //   REQUIRE(!kb.confluent());

  //   kb.run();
  //   REQUIRE(kb.confluent());
  //   REQUIRE(kb.number_of_active_rules() == 32767);
  // }

  // This example verifies the nilpotence of the group using the Sims
  // algorithm. The original presentation was <a,b| [b,a,a,a],
  // [b^-1,a,a,a], [a,b,b,b], [a^-1,b,b,b], [a,a*b,a*b,a*b],
  // [a^-1,a*b,a*b,a*b] >. (where [] mean left-normed commutators. The
  // presentation here was derived by first applying the NQA to find the
  // maximal nilpotent quotient, and then introducing new generators for
  // the PCP generators. LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
  //           //     "020",
  //     "(from kbmag/standalone/kb_data/heinnilp)",
  //     "[fail][knuth-bendix]"[kbmag][recursive], KNUTH_BENDIX_TYPES){
  //   // TODO(later) fails because internal_rewrite expect rules to be
  //   length reducing KnuthBendix<TestType> kb(new RECURSIVE(),
  //   "fFyYdDcCbBaA"); presentation::add_rule_no_checks(p, "BAba", "c");
  //   presentation::add_rule_no_checks(p, "CAca", "d");
  //   presentation::add_rule_no_checks(p, "CBcb", "y");
  //   presentation::add_rule_no_checks(p, "DBdb", "f");
  //   presentation::add_rule_no_checks(p, "cBCb", "bcBC");
  //   presentation::add_rule_no_checks(p, "babABaBA", "abABaBAb");
  //   presentation::add_rule_no_checks(p, "cBACab", "abcBAC");
  //   presentation::add_rule_no_checks(p, "BabABBAbab", "aabABBAb");
  //   auto rg = ReportGuard(false);

  //   REQUIRE(!kb.confluent());

  //   kb.run();
  //   REQUIRE(kb.confluent());
  //   REQUIRE(kb.number_of_active_rules() == 32767);
  //  }

}  // namespace libsemigroups
