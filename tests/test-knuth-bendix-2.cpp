// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2020-2023 James D. Mitchell
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

// This file is the second of six that contains tests for the KnuthBendix
// classes. In a mostly vain attempt to speed up compilation the tests are
// split across 6 files as follows:
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
// 6: contains tests for KnuthBendix created from word_type presentations

#include <cstddef>      // for size_t
#include <cstdint>      // for uint64_t
#include <iostream>     // for string, ostringstream
#include <string>       // for allocator, basic_string
#include <type_traits>  // for is_same_v
#include <utility>      // for move
#include <vector>       // for vector, operator==

#define CATCH_CONFIG_ENABLE_ALL_STRINGMAKERS
#define CATCH_CONFIG_ENABLE_PAIR_STRINGMAKER

#include "catch_amalgamated.hpp"  // for AssertionHandler, oper...
#include "test-main.hpp"          // for LIBSEMIGROUPS_TEMPLATE_TEST_CASE

#include "libsemigroups/constants.hpp"        // for operator==, operator!=
#include "libsemigroups/exception.hpp"        // for LibsemigroupsException
#include "libsemigroups/fpsemi-examples.hpp"  // for chinese
#include "libsemigroups/knuth-bendix.hpp"     // for KnuthBendix, normal_forms
#include "libsemigroups/order.hpp"            // for shortlex_compare
#include "libsemigroups/paths.hpp"            // for Paths
#include "libsemigroups/presentation.hpp"     // for add_rule, Presentation
#include "libsemigroups/to-froidure-pin.hpp"  // for to_froidure_pin
#include "libsemigroups/word-graph.hpp"       // for WordGraph

#include "libsemigroups/detail/report.hpp"  // for ReportGuard

#include "libsemigroups/ranges.hpp"  // for operator|, Inner, to_v...

namespace libsemigroups {

  congruence_kind constexpr twosided = congruence_kind::twosided;

  using namespace rx;

  struct LibsemigroupsException;

  using rule_type = KnuthBendix<>::rule_type;

  using RewriteTrie     = detail::RewriteTrie;
  using RewriteFromLeft = detail::RewriteFromLeft;

#define KNUTH_BENDIX_TYPES RewriteTrie, RewriteFromLeft

  namespace {
    struct weird_cmp {
      bool operator()(rule_type const& x, rule_type const& y) const noexcept {
        return shortlex_compare(x.first, y.first)
               || (x.first == y.first && shortlex_compare(x.second, y.second));
      }
    };
  }  // namespace

  // Fibonacci group F(2,5) - monoid presentation - has order 12 (group
  // elements + empty word)
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "031",
                                   "kbmag/standalone/kb_data/f25monoid",
                                   "[quick][knuth-bendix][kbmag][shortlex]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("abcde");

    presentation::add_rule(p, "ab", "c");
    presentation::add_rule(p, "bc", "d");
    presentation::add_rule(p, "cd", "e");
    presentation::add_rule(p, "de", "a");
    presentation::add_rule(p, "ea", "b");

    KnuthBendix<TestType> kb(twosided, p);
    REQUIRE(!kb.confluent());

    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 24);

    REQUIRE(knuth_bendix::contains(kb, "ab", "c"));
    REQUIRE(knuth_bendix::contains(kb, "bc", "d"));
    REQUIRE(knuth_bendix::contains(kb, "cd", "e"));
    REQUIRE(knuth_bendix::contains(kb, "de", "a"));
    REQUIRE(knuth_bendix::contains(kb, "ea", "b"));
    REQUIRE(knuth_bendix::contains(kb, "cc", "ad"));
    REQUIRE(knuth_bendix::contains(kb, "dd", "be"));
    REQUIRE(knuth_bendix::contains(kb, "ee", "ca"));
    REQUIRE(knuth_bendix::contains(kb, "ec", "bb"));
    REQUIRE(knuth_bendix::contains(kb, "db", "aa"));
    REQUIRE(knuth_bendix::contains(kb, "aac", "be"));
    REQUIRE(knuth_bendix::contains(kb, "bd", "aa"));
    REQUIRE(knuth_bendix::contains(kb, "bbe", "aad"));
    REQUIRE(knuth_bendix::contains(kb, "aaa", "e"));
    REQUIRE(knuth_bendix::contains(kb, "eb", "be"));
    REQUIRE(knuth_bendix::contains(kb, "ba", "c"));
    REQUIRE(knuth_bendix::contains(kb, "da", "ad"));
    REQUIRE(knuth_bendix::contains(kb, "ca", "ac"));
    REQUIRE(knuth_bendix::contains(kb, "ce", "bb"));
    REQUIRE(knuth_bendix::contains(kb, "cb", "d"));
    REQUIRE(knuth_bendix::contains(kb, "ed", "a"));
    REQUIRE(knuth_bendix::contains(kb, "dc", "e"));
    REQUIRE(knuth_bendix::contains(kb, "ae", "b"));
    REQUIRE(knuth_bendix::contains(kb, "bbb", "a"));

    // REQUIRE(knuth_bendix::is_reduced(kb));
    REQUIRE(knuth_bendix::reduce_no_run(kb, "ca") == "ac");
    REQUIRE((kb.active_rules() | sort(weird_cmp()) | to_vector())
            == std::vector<rule_type>(
                {{"ab", "c"},  {"ae", "b"},   {"ba", "c"},  {"bc", "d"},
                 {"bd", "aa"}, {"ca", "ac"},  {"cb", "d"},  {"cc", "ad"},
                 {"cd", "e"},  {"ce", "bb"},  {"da", "ad"}, {"db", "aa"},
                 {"dc", "e"},  {"dd", "be"},  {"de", "a"},  {"ea", "b"},
                 {"eb", "be"}, {"ec", "bb"},  {"ed", "a"},  {"ee", "ac"},
                 {"aaa", "e"}, {"aac", "be"}, {"bbb", "a"}, {"bbe", "aad"}}));

    auto nf = knuth_bendix::normal_forms(kb);
    REQUIRE(
        (nf.min(1).max(5) | to_vector())
        == std::vector<std::string>(
            {"a", "b", "c", "d", "e", "aa", "ac", "ad", "bb", "be", "aad"}));
    REQUIRE(kb.number_of_classes() == 11);
    REQUIRE(nf.min(1).max(POSITIVE_INFINITY).count() == 11);
  }

  // trivial group - BHN presentation
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
      "KnuthBendix",
      "032",
      "kbmag/standalone/kb_data/degen4a",
      "[quick][knuth-bendix][kbmag][shortlex][no-valgrind]",
      KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("aAbBcC");
    p.contains_empty_word(true);
    presentation::add_inverse_rules(p, "AaBbCc");

    presentation::add_rule(p, "Aba", "bb");
    presentation::add_rule(p, "Bcb", "cc");
    presentation::add_rule(p, "Cac", "aa");

    KnuthBendix<TestType> kb(twosided, p);

    REQUIRE(!kb.confluent());

    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 6);

    REQUIRE(knuth_bendix::contains(kb, "Aba", "bb"));
    REQUIRE(knuth_bendix::contains(kb, "Bcb", "cc"));
    REQUIRE(knuth_bendix::contains(kb, "Cac", "aa"));
    REQUIRE((kb.active_rules() | sort(weird_cmp()) | to_vector())
            == std::vector<rule_type>({{"A", ""},
                                       {"B", ""},
                                       {"C", ""},
                                       {"a", ""},
                                       {"b", ""},
                                       {"c", ""}}));
    REQUIRE(kb.number_of_classes() == 1);
    auto nf = knuth_bendix::normal_forms(kb);
    REQUIRE(nf.count() == 1);
  }

  // Torus group
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "033",
                                   "kbmag/standalone/kb_data/torus",
                                   "[quick][knuth-bendix][kbmag][shortlex]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("aAcCbBdD");
    p.contains_empty_word(true);
    presentation::add_inverse_rules(p, "AaCcBbDd");
    presentation::add_rule(p, "ABab", "DCdc");

    KnuthBendix<TestType> kb(twosided, p);
    REQUIRE(!kb.confluent());
    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 16);

    REQUIRE(knuth_bendix::contains(kb, "DCdc", "ABab"));
    REQUIRE((kb.active_rules() | sort(weird_cmp()) | to_vector())
            == std::vector<rule_type>({{"Aa", ""},
                                       {"Bb", ""},
                                       {"Cc", ""},
                                       {"Dd", ""},
                                       {"aA", ""},
                                       {"bB", ""},
                                       {"cC", ""},
                                       {"dD", ""},
                                       {"BAba", "CDcd"},
                                       {"BabC", "aDCd"},
                                       {"DCdc", "ABab"},
                                       {"DcdA", "cBAb"},
                                       {"bCDc", "AbaD"},
                                       {"baDC", "abCD"},
                                       {"dABa", "CdcB"},
                                       {"dcBA", "cdAB"}}));
    REQUIRE(kb.number_of_classes() == POSITIVE_INFINITY);
    auto nf = knuth_bendix::normal_forms(kb).min(0).max(7);
    REQUIRE(nf.count() == 155'577);
    REQUIRE(
        (nf.min(0).max(3) | to_vector())
        == std::vector<std::string>(
            {"",   "a",  "A",  "c",  "C",  "b",  "B",  "d",  "D",  "aa", "ac",
             "aC", "ab", "aB", "ad", "aD", "AA", "Ac", "AC", "Ab", "AB", "Ad",
             "AD", "ca", "cA", "cc", "cb", "cB", "cd", "cD", "Ca", "CA", "CC",
             "Cb", "CB", "Cd", "CD", "ba", "bA", "bc", "bC", "bb", "bd", "bD",
             "Ba", "BA", "Bc", "BC", "BB", "Bd", "BD", "da", "dA", "dc", "dC",
             "db", "dB", "dd", "Da", "DA", "Dc", "DC", "Db", "DB", "DD"}));
  }

  //  3-fold cover of A_6
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
      "KnuthBendix",
      "034",
      "kbmag/standalone/kb_data/3a6",
      "[quick][knuth-bendix][kbmag][shortlex][no-valgrind]",
      KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("abAB");

    presentation::add_inverse_rules(p, "ABab");

    presentation::add_rule(p, "aaa", "");
    presentation::add_rule(p, "bbb", "");
    presentation::add_rule(p, "abababab", "");
    presentation::add_rule(p, "aBaBaBaBaB", "");

    KnuthBendix<TestType> kb(twosided, p);
    REQUIRE(!kb.confluent());

    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 183);

    REQUIRE(knuth_bendix::contains(kb, "aaa", ""));
    REQUIRE(knuth_bendix::contains(kb, "bbb", ""));
    REQUIRE(knuth_bendix::contains(kb, "BaBaBaBaB", "aa"));
    REQUIRE(knuth_bendix::contains(kb, "bababa", "aabb"));
    REQUIRE(knuth_bendix::contains(kb, "ababab", "bbaa"));
    REQUIRE(knuth_bendix::contains(kb, "aabbaa", "babab"));
    REQUIRE(knuth_bendix::contains(kb, "bbaabb", "ababa"));
    REQUIRE(knuth_bendix::contains(kb, "bababbabab", "aabbabbaa"));
    REQUIRE(knuth_bendix::contains(kb, "ababaababa", "bbaabaabb"));
    REQUIRE(knuth_bendix::contains(kb, "bababbabaababa", "aabbabbaabaabb"));
    REQUIRE(knuth_bendix::contains(kb, "bbaabaabbabbaa", "ababaababbabab"));

    REQUIRE(kb.number_of_classes() == 1080);

    auto nf = knuth_bendix::normal_forms(kb);

    REQUIRE(nf.count() == 1080);
    REQUIRE((nf.min(0).max(3) | to_vector())
            == std::vector<std::string>({"",
                                         "a",
                                         "b",
                                         "A",
                                         "B",
                                         "ab",
                                         "aB",
                                         "ba",
                                         "bA",
                                         "Ab",
                                         "AB",
                                         "Ba",
                                         "BA"}));
  }

  //  Free group on 2 generators
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "035",
                                   "kbmag/standalone/kb_data/f2",
                                   "[quick][knuth-bendix][kbmag][shortlex]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("aAbB");
    p.contains_empty_word(true);
    presentation::add_inverse_rules(p, "AaBb");

    KnuthBendix<TestType> kb(twosided, p);
    // kb.process_pending_rules();
    REQUIRE(kb.confluent());
    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 4);
    REQUIRE(kb.number_of_classes() == POSITIVE_INFINITY);

    auto nf = knuth_bendix::normal_forms(kb).min(0).max(3);

    REQUIRE((nf | to_vector())
            == std::vector<std::string>({"",
                                         "a",
                                         "A",
                                         "b",
                                         "B",
                                         "aa",
                                         "ab",
                                         "aB",
                                         "AA",
                                         "Ab",
                                         "AB",
                                         "ba",
                                         "bA",
                                         "bb",
                                         "Ba",
                                         "BA",
                                         "BB"}));
    REQUIRE(nf.min(0).max(5).count() == 161);
  }

  // Symmetric group S_16
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
      "KnuthBendix",
      "036",
      "kbmag/standalone/kb_data/s16",
      "[quick][knuth-bendix][kbmag][shortlex][no-valgrind]",
      KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abcdefghijklmno");
    p.contains_empty_word(true);

    presentation::add_inverse_rules(p, "abcdefghijklmno");

    presentation::add_rule(p, "bab", "aba");
    presentation::add_rule(p, "ca", "ac");
    presentation::add_rule(p, "da", "ad");
    presentation::add_rule(p, "ea", "ae");
    presentation::add_rule(p, "fa", "af");
    presentation::add_rule(p, "ga", "ag");
    presentation::add_rule(p, "ha", "ah");
    presentation::add_rule(p, "ia", "ai");
    presentation::add_rule(p, "ja", "aj");
    presentation::add_rule(p, "ka", "ak");
    presentation::add_rule(p, "la", "al");
    presentation::add_rule(p, "ma", "am");
    presentation::add_rule(p, "na", "an");
    presentation::add_rule(p, "oa", "ao");
    presentation::add_rule(p, "cbc", "bcb");
    presentation::add_rule(p, "db", "bd");
    presentation::add_rule(p, "eb", "be");
    presentation::add_rule(p, "fb", "bf");
    presentation::add_rule(p, "gb", "bg");
    presentation::add_rule(p, "hb", "bh");
    presentation::add_rule(p, "ib", "bi");
    presentation::add_rule(p, "jb", "bj");
    presentation::add_rule(p, "kb", "bk");
    presentation::add_rule(p, "lb", "bl");
    presentation::add_rule(p, "mb", "bm");
    presentation::add_rule(p, "nb", "bn");
    presentation::add_rule(p, "ob", "bo");
    presentation::add_rule(p, "dcd", "cdc");
    presentation::add_rule(p, "ec", "ce");
    presentation::add_rule(p, "fc", "cf");
    presentation::add_rule(p, "gc", "cg");
    presentation::add_rule(p, "hc", "ch");
    presentation::add_rule(p, "ic", "ci");
    presentation::add_rule(p, "jc", "cj");
    presentation::add_rule(p, "kc", "ck");
    presentation::add_rule(p, "lc", "cl");
    presentation::add_rule(p, "mc", "cm");
    presentation::add_rule(p, "nc", "cn");
    presentation::add_rule(p, "oc", "co");
    presentation::add_rule(p, "ede", "ded");
    presentation::add_rule(p, "fd", "df");
    presentation::add_rule(p, "gd", "dg");
    presentation::add_rule(p, "hd", "dh");
    presentation::add_rule(p, "id", "di");
    presentation::add_rule(p, "jd", "dj");
    presentation::add_rule(p, "kd", "dk");
    presentation::add_rule(p, "ld", "dl");
    presentation::add_rule(p, "md", "dm");
    presentation::add_rule(p, "nd", "dn");
    presentation::add_rule(p, "od", "do");
    presentation::add_rule(p, "fef", "efe");
    presentation::add_rule(p, "ge", "eg");
    presentation::add_rule(p, "he", "eh");
    presentation::add_rule(p, "ie", "ei");
    presentation::add_rule(p, "je", "ej");
    presentation::add_rule(p, "ke", "ek");
    presentation::add_rule(p, "le", "el");
    presentation::add_rule(p, "me", "em");
    presentation::add_rule(p, "ne", "en");
    presentation::add_rule(p, "oe", "eo");
    presentation::add_rule(p, "gfg", "fgf");
    presentation::add_rule(p, "hf", "fh");
    presentation::add_rule(p, "if", "fi");
    presentation::add_rule(p, "jf", "fj");
    presentation::add_rule(p, "kf", "fk");
    presentation::add_rule(p, "lf", "fl");
    presentation::add_rule(p, "mf", "fm");
    presentation::add_rule(p, "nf", "fn");
    presentation::add_rule(p, "of", "fo");
    presentation::add_rule(p, "hgh", "ghg");
    presentation::add_rule(p, "ig", "gi");
    presentation::add_rule(p, "jg", "gj");
    presentation::add_rule(p, "kg", "gk");
    presentation::add_rule(p, "lg", "gl");
    presentation::add_rule(p, "mg", "gm");
    presentation::add_rule(p, "ng", "gn");
    presentation::add_rule(p, "og", "go");
    presentation::add_rule(p, "ihi", "hih");
    presentation::add_rule(p, "jh", "hj");
    presentation::add_rule(p, "kh", "hk");
    presentation::add_rule(p, "lh", "hl");
    presentation::add_rule(p, "mh", "hm");
    presentation::add_rule(p, "nh", "hn");
    presentation::add_rule(p, "oh", "ho");
    presentation::add_rule(p, "jij", "iji");
    presentation::add_rule(p, "ki", "ik");
    presentation::add_rule(p, "li", "il");
    presentation::add_rule(p, "mi", "im");
    presentation::add_rule(p, "ni", "in");
    presentation::add_rule(p, "oi", "io");
    presentation::add_rule(p, "kjk", "jkj");
    presentation::add_rule(p, "lj", "jl");
    presentation::add_rule(p, "mj", "jm");
    presentation::add_rule(p, "nj", "jn");
    presentation::add_rule(p, "oj", "jo");
    presentation::add_rule(p, "lkl", "klk");
    presentation::add_rule(p, "mk", "km");
    presentation::add_rule(p, "nk", "kn");
    presentation::add_rule(p, "ok", "ko");
    presentation::add_rule(p, "mlm", "lml");
    presentation::add_rule(p, "nl", "ln");
    presentation::add_rule(p, "ol", "lo");
    presentation::add_rule(p, "nmn", "mnm");
    presentation::add_rule(p, "om", "mo");
    presentation::add_rule(p, "ono", "non");
    KnuthBendix<TestType> kb(twosided, p);

    REQUIRE(!kb.confluent());

    kb.run();  // faster
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 211);  // verified with KBMAG
    REQUIRE(kb.gilman_graph().number_of_nodes() == 121);
    auto g = kb.gilman_graph_node_labels();
    std::sort(g.begin(), g.end(), [](std::string x, std::string y) {
      return shortlex_compare(x, y);
    });
    REQUIRE(g
            == std::vector<std::string>({"",
                                         "a",
                                         "b",
                                         "c",
                                         "d",
                                         "e",
                                         "f",
                                         "g",
                                         "h",
                                         "i",
                                         "j",
                                         "k",
                                         "l",
                                         "m",
                                         "n",
                                         "o",
                                         "ba",
                                         "cb",
                                         "dc",
                                         "ed",
                                         "fe",
                                         "gf",
                                         "hg",
                                         "ih",
                                         "ji",
                                         "kj",
                                         "lk",
                                         "ml",
                                         "nm",
                                         "on",
                                         "cba",
                                         "dcb",
                                         "edc",
                                         "fed",
                                         "gfe",
                                         "hgf",
                                         "ihg",
                                         "jih",
                                         "kji",
                                         "lkj",
                                         "mlk",
                                         "nml",
                                         "onm",
                                         "dcba",
                                         "edcb",
                                         "fedc",
                                         "gfed",
                                         "hgfe",
                                         "ihgf",
                                         "jihg",
                                         "kjih",
                                         "lkji",
                                         "mlkj",
                                         "nmlk",
                                         "onml",
                                         "edcba",
                                         "fedcb",
                                         "gfedc",
                                         "hgfed",
                                         "ihgfe",
                                         "jihgf",
                                         "kjihg",
                                         "lkjih",
                                         "mlkji",
                                         "nmlkj",
                                         "onmlk",
                                         "fedcba",
                                         "gfedcb",
                                         "hgfedc",
                                         "ihgfed",
                                         "jihgfe",
                                         "kjihgf",
                                         "lkjihg",
                                         "mlkjih",
                                         "nmlkji",
                                         "onmlkj",
                                         "gfedcba",
                                         "hgfedcb",
                                         "ihgfedc",
                                         "jihgfed",
                                         "kjihgfe",
                                         "lkjihgf",
                                         "mlkjihg",
                                         "nmlkjih",
                                         "onmlkji",
                                         "hgfedcba",
                                         "ihgfedcb",
                                         "jihgfedc",
                                         "kjihgfed",
                                         "lkjihgfe",
                                         "mlkjihgf",
                                         "nmlkjihg",
                                         "onmlkjih",
                                         "ihgfedcba",
                                         "jihgfedcb",
                                         "kjihgfedc",
                                         "lkjihgfed",
                                         "mlkjihgfe",
                                         "nmlkjihgf",
                                         "onmlkjihg",
                                         "jihgfedcba",
                                         "kjihgfedcb",
                                         "lkjihgfedc",
                                         "mlkjihgfed",
                                         "nmlkjihgfe",
                                         "onmlkjihgf",
                                         "kjihgfedcba",
                                         "lkjihgfedcb",
                                         "mlkjihgfedc",
                                         "nmlkjihgfed",
                                         "onmlkjihgfe",
                                         "lkjihgfedcba",
                                         "mlkjihgfedcb",
                                         "nmlkjihgfedc",
                                         "onmlkjihgfed",
                                         "mlkjihgfedcba",
                                         "nmlkjihgfedcb",
                                         "onmlkjihgfedc",
                                         "nmlkjihgfedcba",
                                         "onmlkjihgfedcb",
                                         "onmlkjihgfedcba"}));
    REQUIRE(kb.gilman_graph().number_of_edges() == 680);

    // verified with KBMAG
    auto nf = knuth_bendix::normal_forms(kb).min(0).max(7);
    REQUIRE(nf.count() == 49'436);

    // verified with KBMAG
    REQUIRE(nf.min(0).max(11).count() == 2'554'607);
    REQUIRE(std::is_same_v<decltype(nf.size_hint()), uint64_t>);
    REQUIRE(nf.max(POSITIVE_INFINITY).size_hint() == 20'922'789'888'000);
    REQUIRE(kb.number_of_classes() == 20'922'789'888'000);
  }

  // Presentation of group A_4 regarded as monoid presentation - gives
  // infinite monoid.
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "037",
                                   "kbmag/standalone/kb_data/a4monoid",
                                   "[quick][knuth-bendix][kbmag][shortlex]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("abB");

    presentation::add_rule(p, "bb", "B");
    presentation::add_rule(p, "BaB", "aba");

    KnuthBendix<TestType> kb(twosided, p);
    REQUIRE(!kb.confluent());

    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 6);

    REQUIRE(knuth_bendix::contains(kb, "bb", "B"));
    REQUIRE(knuth_bendix::contains(kb, "BaB", "aba"));
    REQUIRE(knuth_bendix::contains(kb, "Bb", "bB"));
    REQUIRE(knuth_bendix::contains(kb, "Baaba", "abaaB"));
    REQUIRE(knuth_bendix::contains(kb, "BabB", "abab"));
    REQUIRE(knuth_bendix::contains(kb, "Bababa", "ababaB"));
    REQUIRE((kb.active_rules() | sort(weird_cmp()) | to_vector())
            == std::vector<rule_type>({{{"Bb", "bB"},
                                        {"bb", "B"},
                                        {"BaB", "aba"},
                                        {"BabB", "abab"},
                                        {"Baaba", "abaaB"},
                                        {"Bababa", "ababaB"}}}));
  }

  // fairly clearly the trivial group
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
      "KnuthBendix",
      "038",
      "kbmag/standalone/kb_data/degen3",
      "[quick][knuth-bendix][kbmag][shortlex][no-valgrind]",
      KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("aAbB");
    p.contains_empty_word(true);
    presentation::add_rule(p, "ab", "");
    presentation::add_rule(p, "abb", "");

    KnuthBendix<TestType> kb(twosided, p);
    // kb.process_pending_rules();
    REQUIRE((kb.active_rules() | sort(weird_cmp()) | to_vector())
            == std::vector<rule_type>({{"a", ""}, {"b", ""}}));
    REQUIRE(kb.number_of_active_rules() == 2);
    REQUIRE(kb.confluent());

    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 2);

    REQUIRE(knuth_bendix::contains(kb, "b", ""));
    REQUIRE(knuth_bendix::contains(kb, "a", ""));
    REQUIRE((kb.active_rules() | sort(weird_cmp()) | to_vector())
            == std::vector<rule_type>({{"a", ""}, {"b", ""}}));
  }

  // infinite cyclic group
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "039",
                                   "kbmag/standalone/kb_data/ab1",
                                   "[quick][knuth-bendix][kbmag][shortlex]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("aA");
    p.contains_empty_word(true);
    presentation::add_inverse_rules(p, "Aa");

    KnuthBendix<TestType> kb(twosided, p);
    // kb.process_pending_rules();
    REQUIRE(kb.confluent());

    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 2);
    REQUIRE(kb.number_of_classes() == POSITIVE_INFINITY);
  }

  // A generator, but trivial.
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "040",
                                   "kbmag/standalone/kb_data/degen2",
                                   "[quick][knuth-bendix][kbmag][shortlex]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("aA");
    p.contains_empty_word(true);
    presentation::add_rule(p, "a", "");

    KnuthBendix<TestType> kb(twosided, p);
    // kb.process_pending_rules();
    REQUIRE(kb.confluent());

    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 1);

    REQUIRE(knuth_bendix::contains(kb, "a", ""));
    REQUIRE((kb.active_rules() | sort(weird_cmp()) | to_vector())
            == std::vector<rule_type>({{"a", ""}}));
  }

  // Fibonacci group F(2,5)
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "041",
                                   "kbmag/standalone/kb_data/f25",
                                   "[quick][knuth-bendix][kbmag][shortlex]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("aAbBcCdDyY");

    presentation::add_rule(p, "ab", "c");
    presentation::add_rule(p, "bc", "d");
    presentation::add_rule(p, "cd", "y");
    presentation::add_rule(p, "dy", "a");
    presentation::add_rule(p, "ya", "b");

    KnuthBendix<TestType> kb(twosided, p);
    REQUIRE(!kb.confluent());

    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 24);

    REQUIRE(knuth_bendix::contains(kb, "ab", "c"));
    REQUIRE(knuth_bendix::contains(kb, "bc", "d"));
    REQUIRE(knuth_bendix::contains(kb, "cd", "y"));
    REQUIRE(knuth_bendix::contains(kb, "dy", "a"));
    REQUIRE(knuth_bendix::contains(kb, "ya", "b"));
    REQUIRE(knuth_bendix::contains(kb, "cc", "ad"));
    REQUIRE(knuth_bendix::contains(kb, "dd", "by"));
    REQUIRE(knuth_bendix::contains(kb, "yy", "ac"));
    REQUIRE(knuth_bendix::contains(kb, "yc", "bb"));
    REQUIRE(knuth_bendix::contains(kb, "db", "aa"));
    REQUIRE(knuth_bendix::contains(kb, "aac", "by"));
    REQUIRE(knuth_bendix::contains(kb, "bd", "aa"));
    REQUIRE(knuth_bendix::contains(kb, "bby", "aad"));
    REQUIRE(knuth_bendix::contains(kb, "aaa", "y"));
    REQUIRE(knuth_bendix::contains(kb, "yb", "by"));
    REQUIRE(knuth_bendix::contains(kb, "ba", "c"));
    REQUIRE(knuth_bendix::contains(kb, "da", "ad"));
    REQUIRE(knuth_bendix::contains(kb, "ca", "ac"));
    REQUIRE(knuth_bendix::contains(kb, "cy", "bb"));
    REQUIRE(knuth_bendix::contains(kb, "cb", "d"));
    REQUIRE(knuth_bendix::contains(kb, "yd", "a"));
    REQUIRE(knuth_bendix::contains(kb, "dc", "y"));
    REQUIRE(knuth_bendix::contains(kb, "ay", "b"));
    REQUIRE(knuth_bendix::contains(kb, "bbb", "a"));
    REQUIRE((kb.active_rules() | sort(weird_cmp()) | to_vector())
            == std::vector<rule_type>(
                {{"ab", "c"},  {"ay", "b"},   {"ba", "c"},  {"bc", "d"},
                 {"bd", "aa"}, {"ca", "ac"},  {"cb", "d"},  {"cc", "ad"},
                 {"cd", "y"},  {"cy", "bb"},  {"da", "ad"}, {"db", "aa"},
                 {"dc", "y"},  {"dd", "by"},  {"dy", "a"},  {"ya", "b"},
                 {"yb", "by"}, {"yc", "bb"},  {"yd", "a"},  {"yy", "ac"},
                 {"aaa", "y"}, {"aac", "by"}, {"bbb", "a"}, {"bby", "aad"}}));
  }

  // Von Dyck (2,3,7) group - infinite hyperbolic
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "042",
                                   "kbmag/standalone/kb_data/237",
                                   "[quick][knuth-bendix][kbmag][shortlex]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("aAbBc");
    p.contains_empty_word(true);

    presentation::add_inverse_rules(p, "AaBbc");

    presentation::add_rule(p, "aaaa", "AAA");
    presentation::add_rule(p, "bb", "B");
    presentation::add_rule(p, "BA", "c");

    KnuthBendix<TestType> kb(twosided, p);
    REQUIRE(!kb.confluent());

    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 32);
    REQUIRE((kb.active_rules() | sort(weird_cmp()) | to_vector())
            == std::vector<rule_type>({{"Aa", ""},
                                       {"Ac", "b"},
                                       {"BA", "c"},
                                       {"BB", "b"},
                                       {"Bb", ""},
                                       {"Bc", "bA"},
                                       {"aA", ""},
                                       {"ab", "c"},
                                       {"bB", ""},
                                       {"ba", "AB"},
                                       {"bb", "B"},
                                       {"bc", "A"},
                                       {"cB", "a"},
                                       {"ca", "B"},
                                       {"cb", "aB"},
                                       {"cc", ""},
                                       {"BaB", "bAb"},
                                       {"bAB", "Ba"},
                                       {"cAB", "aBa"},
                                       {"AAAA", "aaa"},
                                       {"AAAb", "aaac"},
                                       {"aaaa", "AAA"},
                                       {"bAbA", "Bac"},
                                       {"cAAA", "Baaa"},
                                       {"cAbA", "aBac"},
                                       {"ABaaa", "bAAA"},
                                       {"Baaac", "cAAb"},
                                       {"bAABaac", "BacAAb"},
                                       {"cAABaac", "aBacAAb"},
                                       {"BaaaBaaa", "cAAbAAA"},
                                       {"bAABaaBaaa", "BacAAbAAA"},
                                       {"cAABaaBaaa", "aBacAAbAAA"}}));
    REQUIRE(kb.number_of_classes() == POSITIVE_INFINITY);
  }

  // Cyclic group of order 2.
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "043",
                                   "kbmag/standalone/kb_data/c2",
                                   "[quick][knuth-bendix][kbmag][shortlex]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("a");
    p.contains_empty_word(true);
    presentation::add_rule(p, "aa", "");

    KnuthBendix<TestType> kb(twosided, p);
    // kb.process_pending_rules();
    REQUIRE(kb.confluent());

    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 1);

    REQUIRE((kb.active_rules() | sort(weird_cmp()) | to_vector())
            == std::vector<rule_type>({{"aa", ""}}));
  }

  // The group is S_4, and the subgroup H of order 4. There are 30 reduced
  // words - 24 for the group elements, and 6 for the 6 cosets Hg.
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "044",
                                   "kbmag/standalone/kb_data/cosets",
                                   "[quick][knuth-bendix][kbmag][shortlex]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("HaAbB");

    presentation::add_rule(p, "aaa", "");
    presentation::add_rule(p, "bbbb", "");
    presentation::add_rule(p, "abab", "");
    presentation::add_rule(p, "Hb", "H");
    presentation::add_rule(p, "HH", "H");
    presentation::add_rule(p, "aH", "H");
    presentation::add_rule(p, "bH", "H");

    KnuthBendix<TestType> kb(twosided, p);
    REQUIRE(!kb.confluent());

    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 14);

    REQUIRE(knuth_bendix::contains(kb, "aaa", ""));
    REQUIRE(knuth_bendix::contains(kb, "Hb", "H"));
    REQUIRE(knuth_bendix::contains(kb, "HH", "H"));
    REQUIRE(knuth_bendix::contains(kb, "aH", "H"));
    REQUIRE(knuth_bendix::contains(kb, "bH", "H"));
    REQUIRE(knuth_bendix::contains(kb, "bab", "aa"));
    REQUIRE(knuth_bendix::contains(kb, "bbb", "aba"));
    REQUIRE(knuth_bendix::contains(kb, "Hab", "Haa"));
    REQUIRE(knuth_bendix::contains(kb, "abaab", "bbaa"));
    REQUIRE(knuth_bendix::contains(kb, "baaba", "aabb"));
    REQUIRE(knuth_bendix::contains(kb, "Haabb", "Haaba"));
    REQUIRE(knuth_bendix::contains(kb, "bbaabb", "abba"));
    REQUIRE(knuth_bendix::contains(kb, "aabbaa", "baab"));
    REQUIRE(knuth_bendix::contains(kb, "baabba", "abbaab"));
    REQUIRE((kb.active_rules() | sort(weird_cmp()) | to_vector())
            == std::vector<rule_type>({{{"HH", "H"},
                                        {"Hb", "H"},
                                        {"aH", "H"},
                                        {"bH", "H"},
                                        {"Hab", "Haa"},
                                        {"aaa", ""},
                                        {"bab", "aa"},
                                        {"bbb", "aba"},
                                        {"Haabb", "Haaba"},
                                        {"abaab", "bbaa"},
                                        {"baaba", "aabb"},
                                        {"aabbaa", "baab"},
                                        {"baabba", "abbaab"},
                                        {"bbaabb", "abba"}}}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "045",
                                   "Ex. 5.1 in Sims (KnuthBendix 09 again)",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("aAbB");
    p.contains_empty_word(true);

    presentation::add_rule(p, "aA", "");
    presentation::add_rule(p, "Aa", "");
    presentation::add_rule(p, "bB", "");
    presentation::add_rule(p, "Bb", "");
    presentation::add_rule(p, "ba", "ab");

    KnuthBendix<TestType> kb(twosided, p);
    REQUIRE(!kb.confluent());

    kb.run();
    REQUIRE(kb.number_of_active_rules() == 8);
    REQUIRE(kb.confluent());
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "046",
                                   "kbmag/standalone/kb_data/nilp2",
                                   "[quick][knuth-bendix][kbmag][shortlex]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("cCbBaA");
    p.contains_empty_word(true);

    presentation::add_inverse_rules(p, "CcBbAa");

    presentation::add_rule(p, "ba", "abc");
    presentation::add_rule(p, "ca", "ac");
    presentation::add_rule(p, "cb", "bc");

    KnuthBendix<TestType> kb(twosided, p);
    REQUIRE(!kb.confluent());
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "047",
                                   "Ex. 6.4 in Sims",
                                   "[quick][knuth-bendix][no-valgrind]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("abc");
    presentation::add_rule(p, "aa", "");
    presentation::add_rule(p, "bc", "");
    presentation::add_rule(p, "bbb", "");
    presentation::add_rule(p, "ababababababab", "");
    presentation::add_rule(p, "abacabacabacabac", "");

    KnuthBendix<TestType> kb(twosided, p);
    // kb.process_pending_rules();
    REQUIRE(kb.number_of_active_rules() == 5);
    REQUIRE(!kb.confluent());

    kb.max_rules(10);
    kb.run();
    REQUIRE(kb.number_of_active_rules() > 10);
    REQUIRE(!kb.confluent());

    kb.run();
    REQUIRE(kb.number_of_active_rules() > 10);
    REQUIRE(!kb.confluent());

    kb.max_rules(20);
    kb.run();
    REQUIRE(kb.number_of_active_rules() > 20);
    REQUIRE(!kb.confluent());

    kb.max_rules(LIMIT_MAX);
    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 40);
  }

  // Von Dyck (2,3,7) group - infinite hyperbolic
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "KnuthBendix 071 again",
                                   "[no-valgrind]048",
                                   "[quick][knuth-bendix][shortlex]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("aAbBc");
    p.contains_empty_word(true);

    presentation::add_inverse_rules(p, "AaBbc");

    presentation::add_rule(p, "BA", "c");
    presentation::add_rule(p, "Bb", "bB");
    presentation::add_rule(p, "bb", "B");
    presentation::add_rule(p, "AAAa", "aAAA");
    presentation::add_rule(p, "aaaa", "AAA");
    presentation::add_rule(p, "BaAAA", "cAAa");
    presentation::add_rule(p, "BaaAAA", "cAAaa");
    presentation::add_rule(p, "BaAaAAA", "cAAaAa");
    presentation::add_rule(p, "BaaaAAA", "cAAaaa");
    presentation::add_rule(p, "BaAAaAAA", "cAAaAAa");
    presentation::add_rule(p, "BaAaaAAA", "cAAaAaa");
    presentation::add_rule(p, "BaaAaAAA", "cAAaaAa");
    presentation::add_rule(p, "BaAAaaAAA", "cAAaAAaa");
    presentation::add_rule(p, "BaAaAaAAA", "cAAaAaAa");
    presentation::add_rule(p, "BaAaaaAAA", "cAAaAaaa");
    presentation::add_rule(p, "BaaAAaAAA", "cAAaaAAa");
    presentation::add_rule(p, "BaaAaaAAA", "cAAaaAaa");
    presentation::add_rule(p, "BaAAaAaAAA", "cAAaAAaAa");
    presentation::add_rule(p, "BaAAaaaAAA", "cAAaAAaaa");
    presentation::add_rule(p, "BaAaAAaAAA", "cAAaAaAAa");
    presentation::add_rule(p, "BaAaAaaAAA", "cAAaAaAaa");
    presentation::add_rule(p, "BaAaaAaAAA", "cAAaAaaAa");
    presentation::add_rule(p, "BaaAAaaAAA", "cAAaaAAaa");
    presentation::add_rule(p, "BaaAaAaAAA", "cAAaaAaAa");
    presentation::add_rule(p, "BaAAaAAaAAA", "cAAaAAaAAa");
    presentation::add_rule(p, "BaAAaAaaAAA", "cAAaAAaAaa");
    presentation::add_rule(p, "BaAAaaAaAAA", "cAAaAAaaAa");
    presentation::add_rule(p, "BaAaAAaaAAA", "cAAaAaAAaa");
    presentation::add_rule(p, "BaAaAaAaAAA", "cAAaAaAaAa");
    presentation::add_rule(p, "BaAaaAAaAAA", "cAAaAaaAAa");
    presentation::add_rule(p, "BaaAAaAaAAA", "cAAaaAAaAa");
    presentation::add_rule(p, "BaaAaAAaAAA", "cAAaaAaAAa");
    presentation::add_rule(p, "BaAAaAAaaAAA", "cAAaAAaAAaa");
    presentation::add_rule(p, "BaAAaAaAaAAA", "cAAaAAaAaAa");
    presentation::add_rule(p, "BaAAaaAAaAAA", "cAAaAAaaAAa");
    presentation::add_rule(p, "BaAaAAaAaAAA", "cAAaAaAAaAa");
    presentation::add_rule(p, "BaAaAaAAaAAA", "cAAaAaAaAAa");
    presentation::add_rule(p, "BaaAAaAAaAAA", "cAAaaAAaAAa");
    presentation::add_rule(p, "BaAAaAAaAaAAA", "cAAaAAaAAaAa");
    presentation::add_rule(p, "BaAAaAaAAaAAA", "cAAaAAaAaAAa");
    presentation::add_rule(p, "BaAaAAaAAaAAA", "cAAaAaAAaAAa");
    presentation::add_rule(p, "BaAAaAAaAAaAAA", "cAAaAAaAAaAAa");

    KnuthBendix<TestType> kb(twosided, p);
    // kb.process_pending_rules();
    REQUIRE(kb.number_of_active_rules() == 9);
    REQUIRE(!kb.confluent());
    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 32);
    REQUIRE(kb.number_of_classes() == POSITIVE_INFINITY);
    auto nf = knuth_bendix::normal_forms(kb);

    REQUIRE(nf.min(4).max(5).count() == 24);
    nf.min(4);
    REQUIRE((nf | to_vector())
            == std::vector<std::string>({"aaaB", "aaac", "aaBa", "aacA", "aBaa",
                                         "aBac", "acAA", "acAb", "AAAB", "AAbA",
                                         "AABa", "AbAA", "AbAb", "ABaa", "ABac",
                                         "bAAA", "bAAb", "bAAB", "Baaa", "BaaB",
                                         "Baac", "BacA", "cAAb", "cAAB"}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "049",
                                   "Sims Ex. 5.4 - alt. overlap policy",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("Bab");
    p.contains_empty_word(true);
    presentation::add_rule(p, "aa", "");
    presentation::add_rule(p, "bB", "");
    presentation::add_rule(p, "bbb", "");
    presentation::add_rule(p, "ababab", "");

    KnuthBendix<TestType> kb(twosided, p);
    kb.overlap_policy(KnuthBendix<TestType>::options::overlap::AB_BC);

    REQUIRE(!kb.confluent());

    knuth_bendix::by_overlap_length(kb);
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 11);
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_classes() == 12);
    auto nf1 = knuth_bendix::normal_forms(kb).min(4).max(5);
    REQUIRE(nf1.count() == 0);

    REQUIRE((knuth_bendix::normal_forms(kb) | to_vector())
            == std::vector<std::string>({"",
                                         "B",
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
                                   "050",
                                   "Sims - Ex. 5.4 - alt. overlap policy",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("Bab");
    presentation::add_rule(p, "aa", "");
    presentation::add_rule(p, "bB", "");
    presentation::add_rule(p, "bbb", "");
    presentation::add_rule(p, "ababab", "");

    KnuthBendix<TestType> kb(twosided, p);
    kb.overlap_policy(KnuthBendix<TestType>::options::overlap::MAX_AB_BC);
    // The next line tests that we don't delete
    // the old OverlapMeasure.
    kb.overlap_policy(KnuthBendix<TestType>::options::overlap::MAX_AB_BC);

    REQUIRE(!kb.confluent());

    knuth_bendix::by_overlap_length(kb);
    REQUIRE(kb.number_of_active_rules() == 11);
    REQUIRE(kb.confluent());
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "051",
                                   "operator<<",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    std::ostringstream os;

    Presentation<std::string> p;
    p.alphabet("Bab");
    p.contains_empty_word(true);
    presentation::add_rule(p, "aa", "");
    presentation::add_rule(p, "bB", "");
    presentation::add_rule(p, "bbb", "");
    presentation::add_rule(p, "ababab", "");

    KnuthBendix<TestType> kb1(twosided, p);
    os << kb1;  // Does not do anything visible
    p.alphabet("cbaB");
    presentation::add_rule(p, "aa", "");
    presentation::add_rule(p, "bB", "");
    presentation::add_rule(p, "bbb", "");
    presentation::add_rule(p, "ababab", "");
    KnuthBendix<TestType> kb2(twosided, p);
    os << kb2;  // Does not do anything visible
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "052",
                                   "confluence_interval",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("Bab");
    presentation::add_rule(p, "aa", "");
    presentation::add_rule(p, "bB", "");
    presentation::add_rule(p, "bbb", "");
    presentation::add_rule(p, "ababab", "");
    KnuthBendix<TestType> kb(twosided, p);
    kb.check_confluence_interval(LIMIT_MAX);
    kb.check_confluence_interval(10);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "053",
                                   "max_overlap",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("Bab");

    presentation::add_rule(p, "aa", "");
    presentation::add_rule(p, "bB", "");
    presentation::add_rule(p, "bbb", "");
    presentation::add_rule(p, "ababab", "");

    KnuthBendix<TestType> kb(twosided, p);
    kb.max_overlap(10);
    kb.max_overlap(-11);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
      "KnuthBendix",
      "054",
      "kbmag/standalone/kb_data/d22",
      "[quick][knuth-bendix][fpsemi][kbmag][shortlex]",
      KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("ABCDYFabcdyf");
    p.contains_empty_word(true);

    presentation::add_inverse_rules(p, "abcdyfABCDYF");

    presentation::add_rule(p, "aCAd", "");
    presentation::add_rule(p, "bfBY", "");
    presentation::add_rule(p, "cyCD", "");
    presentation::add_rule(p, "dFDa", "");
    presentation::add_rule(p, "ybYA", "");
    presentation::add_rule(p, "fCFB", "");

    KnuthBendix<TestType> kb(twosided, p);
    REQUIRE(!kb.confluent());

    knuth_bendix::by_overlap_length(kb);
    REQUIRE(kb.confluent());
    REQUIRE(kb.presentation().rules.size() / 2 == 18);
    REQUIRE(kb.number_of_classes() == 22);

    REQUIRE((knuth_bendix::normal_forms(kb) | to_vector())
            == std::vector<std::string>({"",    "A",   "B",   "C",  "D",  "Y",
                                         "F",   "AB",  "AC",  "AD", "AY", "AF",
                                         "BA",  "BD",  "BY",  "CY", "DB", "ABA",
                                         "ABD", "ABY", "ACY", "ADB"}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
      "KnuthBendix",
      "055",
      "kbmag/standalone/kb_data/d22",
      "[quick][knuth-bendix][fpsemi][kbmag][shortlex]",
      KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("aAbBcCdDyYfF");
    p.contains_empty_word(true);

    presentation::add_inverse_rules(p, "AaBbCcDdYyFf");

    presentation::add_rule(p, "aCAd", "");
    presentation::add_rule(p, "bfBY", "");
    presentation::add_rule(p, "cyCD", "");
    presentation::add_rule(p, "dFDa", "");
    presentation::add_rule(p, "ybYA", "");
    presentation::add_rule(p, "fCFB", "");
    KnuthBendix<TestType> kb(twosided, p);
    REQUIRE(!kb.confluent());

    knuth_bendix::by_overlap_length(kb);
    REQUIRE(kb.confluent());
    REQUIRE(kb.presentation().rules.size() / 2 == 18);
    REQUIRE(kb.number_of_classes() == 22);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "056",
                                   "small example",
                                   "[quick][knuth-bendix][shortlex]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("ab");
    presentation::add_rule(p, "aaa", "a");
    presentation::add_rule(p, "bbbb", "b");
    presentation::add_rule(p, "ababababab", "aa");
    KnuthBendix<TestType> kb(twosided, p);
    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_classes() == 243);
    auto nf = knuth_bendix::normal_forms(kb).min(1).max(3);
    REQUIRE((nf | to_vector())
            == std::vector<std::string>({"a", "b", "aa", "ab", "ba", "bb"}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "057",
                                   "code coverage",
                                   "[quick]",
                                   KNUTH_BENDIX_TYPES) {
    KnuthBendix<TestType> kb1;
    KnuthBendix<TestType> kb2(kb1);
    REQUIRE(kb1.number_of_classes() == 0);

    Presentation<std::string> p;
    p.alphabet("ab");
    presentation::add_rule(p, "aaa", "a");
    KnuthBendix<TestType> kb3(twosided, p);
    REQUIRE(kb3.presentation().rules.size() / 2 == 1);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "058",
                                   "small overlap 1",
                                   "[quick]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("BCA");
    presentation::add_rule(p, "AABC", "ACBA");

    KnuthBendix<TestType> kb(twosided, p);
    // kb.process_pending_rules();
    REQUIRE(kb.confluent());
    REQUIRE(knuth_bendix::reduce(kb, "CBACBAABCAABCACBACBA")
            == "CBACBACBAACBAACBACBA");
    REQUIRE(knuth_bendix::contains(
        kb, "CBAABCABCAABCAABCABC", "CBACBAABCAABCACBACBA"));
    REQUIRE(knuth_bendix::contains(
        kb, "CBAABCABCAABCAABCABC", "CBACBAABCAABCACBACBA"));
    REQUIRE(knuth_bendix::contains(kb,
                                   "AABCAABCCACAACBBCBCCACBBAABCBA",
                                   "ACBAACBACACAACBBCBCCACBBACBABA"));
    REQUIRE(knuth_bendix::contains(kb,
                                   "CACCBABACCBABACCAAAABCAABCBCAA",
                                   "CACCBABACCBABACCAAACBAACBABCAA"));
    REQUIRE(knuth_bendix::contains(kb,
                                   "CAAACAABCCBABCCBCCBCACABACBBAC",
                                   "CAAACACBACBABCCBCCBCACABACBBAC"));
    REQUIRE(knuth_bendix::contains(kb,
                                   "BABCACBACBCCCCCAACCAAABAABCBCC",
                                   "BABCACBACBCCCCCAACCAAABACBABCC"));
    REQUIRE(kb.number_of_classes() == POSITIVE_INFINITY);
  }

  // Symmetric group S_9
  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "059",
                                   "kbmag/standalone/kb_data/s9",
                                   "[quick][knuth-bendix][kbmag][shortlex]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("abcdefgh");
    p.contains_empty_word(true);

    presentation::add_inverse_rules(p, "abcdefgh");

    presentation::add_rule(p, "bab", "aba");
    presentation::add_rule(p, "ca", "ac");
    presentation::add_rule(p, "da", "ad");
    presentation::add_rule(p, "ea", "ae");
    presentation::add_rule(p, "fa", "af");
    presentation::add_rule(p, "ga", "ag");
    presentation::add_rule(p, "ha", "ah");
    presentation::add_rule(p, "cbc", "bcb");
    presentation::add_rule(p, "db", "bd");
    presentation::add_rule(p, "eb", "be");
    presentation::add_rule(p, "fb", "bf");
    presentation::add_rule(p, "gb", "bg");
    presentation::add_rule(p, "hb", "bh");
    presentation::add_rule(p, "dcd", "cdc");
    presentation::add_rule(p, "ec", "ce");
    presentation::add_rule(p, "fc", "cf");
    presentation::add_rule(p, "gc", "cg");
    presentation::add_rule(p, "hc", "ch");
    presentation::add_rule(p, "ede", "ded");
    presentation::add_rule(p, "fd", "df");
    presentation::add_rule(p, "gd", "dg");
    presentation::add_rule(p, "hd", "dh");
    presentation::add_rule(p, "fef", "efe");
    presentation::add_rule(p, "ge", "eg");
    presentation::add_rule(p, "he", "eh");
    presentation::add_rule(p, "gfg", "fgf");
    presentation::add_rule(p, "hf", "fh");
    presentation::add_rule(p, "hgh", "ghg");

    KnuthBendix<TestType> kb(twosided, p);
    REQUIRE(!kb.confluent());
    kb.run();
    REQUIRE(kb.confluent());
    REQUIRE(kb.number_of_active_rules() == 57);
    REQUIRE(kb.number_of_classes() == 362'880);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "060",
                                   "C(4) monoid",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    Presentation<std::string> p;
    p.alphabet("abcde");
    presentation::add_rule(p, "bceac", "aeebbc");
    presentation::add_rule(p, "aeebbc", "dabcd");

    KnuthBendix<TestType> kb(twosided, p);
    kb.run();
    REQUIRE(kb.confluent());
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "061",
                                   "1-relation hard case",
                                   "[fail][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(true);
    Presentation<std::string> p;
    p.alphabet("ab");
    p.contains_empty_word(true);
    presentation::add_rule(p, "baaababaaa", "aaba");

    KnuthBendix<TestType> kb(twosided, p);
    // knuth_bendix::by_overlap_length(kb);
    REQUIRE(!kb.confluent());
    kb.run();
    REQUIRE(kb.confluent());
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "062",
                                   "1-relation hard case x 2",
                                   "[quick][knuth-bendix][no-valgrind]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("abcd");
    presentation::add_rule(p, "aa", "a");
    presentation::add_rule(p, "ad", "d");
    presentation::add_rule(p, "bb", "b");
    presentation::add_rule(p, "ca", "ac");
    presentation::add_rule(p, "cc", "c");
    presentation::add_rule(p, "da", "d");
    presentation::add_rule(p, "dc", "cd");
    presentation::add_rule(p, "dd", "d");
    presentation::add_rule(p, "aba", "a");
    presentation::add_rule(p, "abd", "d");
    presentation::add_rule(p, "acd", "cd");
    presentation::add_rule(p, "bab", "b");
    presentation::add_rule(p, "bcb", "b");
    presentation::add_rule(p, "bcd", "cd");
    presentation::add_rule(p, "cbc", "c");
    presentation::add_rule(p, "cdb", "cd");
    presentation::add_rule(p, "dba", "d");
    presentation::add_rule(p, "dbd", "d");
    presentation::add_rule(p, "acba", "ac");
    presentation::add_rule(p, "acbd", "cd");
    presentation::add_rule(p, "cbac", "ac");
    auto it = knuth_bendix::redundant_rule(p, std::chrono::milliseconds(100));
    while (it != p.rules.end()) {
      // std::cout << std::endl
      //           << "REMOVING " << *it << " = " << *(it + 1) << std::endl;
      p.rules.erase(it, it + 2);
      it = knuth_bendix::redundant_rule(p, std::chrono::milliseconds(100));
    }
    REQUIRE(p.rules
            == std::vector<std::string>(
                {"aa",  "a", "ad",  "d",  "bb",  "b", "ca",  "ac", "cc",  "c",
                 "da",  "d", "dc",  "cd", "dd",  "d", "aba", "a",  "bab", "b",
                 "bcb", "b", "bcd", "cd", "cbc", "c", "cdb", "cd"}));
    KnuthBendix<TestType> kb(congruence_kind::twosided, p);
    REQUIRE(kb.number_of_classes() == 24);
    REQUIRE(knuth_bendix::reduce(kb, "dcb") == "cd");
    REQUIRE(knuth_bendix::reduce(kb, "dca") == "cd");
    REQUIRE(knuth_bendix::reduce(kb, "da") == "d");
    REQUIRE(knuth_bendix::reduce(kb, "cda") == "cd");
    REQUIRE(knuth_bendix::reduce(kb, "cdb") == "cd");
    REQUIRE(knuth_bendix::reduce(kb, "cdc") == "cd");
    REQUIRE(knuth_bendix::reduce(kb, "cdd") == "cd");
    REQUIRE(knuth_bendix::reduce(kb, "dad") == "d");
    REQUIRE(!knuth_bendix::contains(kb, "bd", "db"));
    REQUIRE(knuth_bendix::reduce(kb, "bd") == "bd");
    REQUIRE(knuth_bendix::reduce(kb, "db") == "db");
    REQUIRE(knuth_bendix::reduce(kb, "cbdcbd") == "cd");
    REQUIRE((knuth_bendix::normal_forms(kb) | to_vector())
            == std::vector<std::string>(
                {"",    "a",   "b",   "c",   "d",    "ab",   "ac",   "ba",
                 "bc",  "bd",  "cb",  "cd",  "db",   "abc",  "acb",  "bac",
                 "bdb", "cba", "cbd", "dbc", "bacb", "bdbc", "cbdb", "cbdbc"}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "063",
                                   "search for a monoid that might not exist",
                                   "[quick][knuth-bendix]",
                                   KNUTH_BENDIX_TYPES) {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.contains_empty_word(true);
    p.alphabet("abcde");
    presentation::add_rule(p, "aa", "a");
    presentation::add_rule(p, "ad", "d");
    presentation::add_rule(p, "bb", "b");
    presentation::add_rule(p, "ca", "ac");
    presentation::add_rule(p, "cc", "c");
    presentation::add_rule(p, "da", "d");
    presentation::add_rule(p, "dc", "cd");
    presentation::add_rule(p, "dd", "d");
    presentation::add_rule(p, "aba", "a");
    presentation::add_rule(p, "bab", "b");
    presentation::add_rule(p, "bcb", "b");
    presentation::add_rule(p, "bcd", "cd");
    presentation::add_rule(p, "cbc", "c");
    presentation::add_rule(p, "cdb", "cd");
    presentation::change_alphabet(p, "cbade");

    presentation::add_rule(p, "ea", "ae");
    presentation::add_rule(p, "be", "eb");
    presentation::add_rule(p, "ee", "e");
    presentation::add_rule(p, "cec", "c");
    presentation::add_rule(p, "ece", "e");

    presentation::add_rule(p, "ead", "ad");
    presentation::add_rule(p, "ade", "ad");
    // presentation::add_rule(p, "de", "ed");
    KnuthBendix<TestType> kb(congruence_kind::twosided, p);
    REQUIRE(kb.number_of_classes() == POSITIVE_INFINITY);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "064",
                                   "Chinese monoid",
                                   "[knuth-bendix][quick][no-valgrind]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);
    // fmt::print(bg(fmt::color::white) | fg(fmt::color::black),
    //            "062",": Chinese monoid STARTING . . .\n");

    std::array<uint64_t, 11> const num
        = {0, 0, 22, 71, 181, 391, 750, 1'317, 2'161, 3'361, 5'006};

    for (size_t n = 2; n < 11; ++n) {
      auto p = fpsemigroup::chinese_monoid(n);
      p.contains_empty_word(true);
      KnuthBendix<TestType> kb(twosided, p);
      kb.run();
      REQUIRE(knuth_bendix::normal_forms(kb).min(0).max(5).count() == num[n]);
    }
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "065",
                                   "hypostylic",
                                   "[todd-coxeter][quick]",
                                   KNUTH_BENDIX_TYPES) {
    using namespace literals;
    using words::operator+;

    auto   rg = ReportGuard(false);
    size_t n  = 2;
    auto   p  = fpsemigroup::hypo_plactic_monoid(n);
    p.contains_empty_word(true);
    presentation::add_idempotent_rules_no_checks(
        p, (seq<size_t>() | take(n) | to_vector()));
    KnuthBendix<TestType> kb(congruence_kind::twosided, p);
    kb.run();
    REQUIRE(kb.presentation().alphabet() == std::string({0, 1}));
    REQUIRE((knuth_bendix::normal_forms<word_type>(kb)
             | filter([&kb](auto const& w) {
                 auto ww = w;
                 ww.insert(ww.begin(), w.begin(), w.end());
                 return knuth_bendix::reduce(kb, ww) == w;
               })
             | to_vector())
            == std::vector<word_type>({{}, 0_w, 1_w, 10_w}));
    REQUIRE((kb.active_rules() | sort(weird_cmp()) | to_vector())
            == std::vector<std::pair<std::string, std::string>>(
                {{{0, 0}, {0}},
                 {{1, 1}, {1}},
                 {{0, 1, 0}, {1, 0}},
                 {{1, 0, 1}, {1, 0}}}));
    // The gilman_graph generated is isomorphic to the word_graph given, but not
    // identical. Since the normal forms are correct (see above) the below check
    // is omitted. REQUIRE(kb.gilman_graph()
    //         == make<WordGraph<size_t>>(5, {{1, 3}, {UNDEFINED, 2}, {},
    //         {4}}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "066",
                                   "Chinese id monoid",
                                   "[todd-coxeter][quick]",
                                   KNUTH_BENDIX_TYPES) {
    auto rg = ReportGuard(false);
    auto n  = 4;
    auto p  = fpsemigroup::chinese_monoid(n);
    p.contains_empty_word(true);
    presentation::add_idempotent_rules_no_checks(p, p.alphabet());

    KnuthBendix<TestType> kb(twosided, to_presentation<std::string>(p));
    kb.run();

    REQUIRE(knuth_bendix::reduce(kb, "cbda") == "bcda");
    REQUIRE(knuth_bendix::reduce(kb, "badc") == "badc");
    REQUIRE(knuth_bendix::reduce(kb, "cadb") == "cadb");
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "067",
                                   "sigma sylvester monoid",
                                   "[todd-coxeter][quick][no-valgrind]",
                                   KNUTH_BENDIX_TYPES) {
    using namespace literals;
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(4);
    p.contains_empty_word(true);
    presentation::add_rule(p, 00_w, 0_w);
    presentation::add_rule(p, 11_w, 1_w);
    presentation::add_rule(p, 22_w, 2_w);
    presentation::add_rule(p, 33_w, 3_w);
    presentation::add_rule(p, 010_w, 01_w);
    presentation::add_rule(p, 020_w, 02_w);
    presentation::add_rule(p, 030_w, 03_w);
    presentation::add_rule(p, 121_w, 12_w);
    presentation::add_rule(p, 131_w, 13_w);
    presentation::add_rule(p, 232_w, 23_w);
    presentation::add_rule(p, 0120_w, 012_w);
    presentation::add_rule(p, 0130_w, 013_w);
    presentation::add_rule(p, 0210_w, 021_w);
    presentation::add_rule(p, 0230_w, 023_w);
    presentation::add_rule(p, 0310_w, 031_w);
    presentation::add_rule(p, 0320_w, 032_w);
    presentation::add_rule(p, 1202_w, 120_w);
    presentation::add_rule(p, 1231_w, 123_w);
    presentation::add_rule(p, 1303_w, 130_w);
    presentation::add_rule(p, 1321_w, 132_w);
    presentation::add_rule(p, 2303_w, 230_w);
    presentation::add_rule(p, 2313_w, 231_w);
    presentation::add_rule(p, 01230_w, 0123_w);
    presentation::add_rule(p, 01320_w, 0132_w);
    presentation::add_rule(p, 02120_w, 0212_w);
    presentation::add_rule(p, 02130_w, 0213_w);
    presentation::add_rule(p, 02310_w, 0231_w);
    presentation::add_rule(p, 03120_w, 0312_w);
    presentation::add_rule(p, 03130_w, 0313_w);
    presentation::add_rule(p, 03210_w, 0321_w);
    presentation::add_rule(p, 03230_w, 0323_w);
    presentation::add_rule(p, 10212_w, 1021_w);
    presentation::add_rule(p, 10313_w, 1031_w);
    presentation::add_rule(p, 12012_w, 1201_w);
    presentation::add_rule(p, 12032_w, 1203_w);
    presentation::add_rule(p, 12302_w, 1230_w);
    presentation::add_rule(p, 13013_w, 1301_w);
    presentation::add_rule(p, 13202_w, 1320_w);
    presentation::add_rule(p, 13231_w, 1323_w);
    presentation::add_rule(p, 20313_w, 2031_w);
    presentation::add_rule(p, 20323_w, 2032_w);
    presentation::add_rule(p, 21323_w, 2132_w);
    presentation::add_rule(p, 23013_w, 2301_w);
    presentation::add_rule(p, 23023_w, 2302_w);
    presentation::add_rule(p, 23103_w, 2310_w);
    presentation::add_rule(p, 23123_w, 2312_w);
    presentation::add_rule(p, 013230_w, 01323_w);
    presentation::add_rule(p, 021230_w, 02123_w);
    presentation::add_rule(p, 021320_w, 02132_w);
    presentation::add_rule(p, 023120_w, 02312_w);
    presentation::add_rule(p, 031230_w, 03123_w);
    presentation::add_rule(p, 031320_w, 03132_w);
    presentation::add_rule(p, 032120_w, 03212_w);
    presentation::add_rule(p, 032130_w, 03213_w);
    presentation::add_rule(p, 032310_w, 03231_w);
    presentation::add_rule(p, 102132_w, 10213_w);
    presentation::add_rule(p, 102312_w, 10231_w);
    presentation::add_rule(p, 103212_w, 10321_w);
    presentation::add_rule(p, 120132_w, 12013_w);
    presentation::add_rule(p, 120312_w, 12031_w);
    presentation::add_rule(p, 123012_w, 12301_w);
    presentation::add_rule(p, 130212_w, 13021_w);
    presentation::add_rule(p, 132012_w, 13201_w);
    presentation::add_rule(p, 132032_w, 13203_w);
    presentation::add_rule(p, 132302_w, 13230_w);
    presentation::add_rule(p, 201323_w, 20132_w);
    presentation::add_rule(p, 203123_w, 20312_w);
    presentation::add_rule(p, 203213_w, 20321_w);
    presentation::add_rule(p, 210323_w, 21032_w);
    presentation::add_rule(p, 213023_w, 21302_w);
    presentation::add_rule(p, 213203_w, 21320_w);
    presentation::add_rule(p, 230123_w, 23012_w);
    presentation::add_rule(p, 230213_w, 23021_w);
    presentation::add_rule(p, 231013_w, 23101_w);
    presentation::add_rule(p, 231023_w, 23102_w);
    presentation::add_rule(p, 231203_w, 23120_w);
    presentation::add_rule(p, 0313230_w, 031323_w);
    presentation::add_rule(p, 0321230_w, 032123_w);
    presentation::add_rule(p, 0321320_w, 032132_w);
    presentation::add_rule(p, 0323120_w, 032312_w);
    presentation::add_rule(p, 1032132_w, 103213_w);
    presentation::add_rule(p, 1032312_w, 103231_w);
    presentation::add_rule(p, 1302132_w, 130213_w);
    presentation::add_rule(p, 1302312_w, 130231_w);
    presentation::add_rule(p, 1320132_w, 132013_w);
    presentation::add_rule(p, 1320312_w, 132031_w);
    presentation::add_rule(p, 1323012_w, 132301_w);
    presentation::add_rule(p, 2032123_w, 203212_w);
    presentation::add_rule(p, 2101323_w, 210132_w);
    presentation::add_rule(p, 2103123_w, 210312_w);
    presentation::add_rule(p, 2103213_w, 210321_w);
    presentation::add_rule(p, 2130123_w, 213012_w);
    presentation::add_rule(p, 2130213_w, 213021_w);
    presentation::add_rule(p, 2132013_w, 213201_w);
    presentation::add_rule(p, 2302123_w, 230212_w);
    presentation::add_rule(p, 2310123_w, 231012_w);
    presentation::add_rule(p, 2310213_w, 231021_w);
    presentation::add_rule(p, 2312013_w, 231201_w);
    REQUIRE(p.rules.size() == 196);
    // {
    //   ReportGuard rg(false);
    //   auto it = knuth_bendix::redundant_rule(p,
    //   std::chrono::milliseconds(100)); while (it != p.rules.cend()) {
    //     p.rules.erase(it, it + 2);
    //     it = knuth_bendix::redundant_rule(p, std::chrono::milliseconds(100));
    //   }

    //   REQUIRE(p.rules.size() == 58);
    //   REQUIRE(
    //       p.rules
    //       == std::vector<word_type>(
    //           {00_w,      0_w,      11_w,      1_w,      22_w,      2_w,
    //            33_w,      3_w,      010_w,     01_w,     020_w,     02_w,
    //            030_w,     03_w,     121_w,     12_w,     131_w,     13_w,
    //            232_w,     23_w,     1202_w,    120_w,    1303_w,    130_w,
    //            2303_w,    230_w,    2313_w,    231_w,    10212_w,   1021_w,
    //            10313_w,   1031_w,   20313_w,   2031_w,   20323_w,   2032_w,
    //            21323_w,   2132_w,   102312_w,  10231_w,  103212_w,  10321_w,
    //            201323_w,  20132_w,  203123_w,  20312_w,  210323_w,  21032_w,
    //            213023_w,  21302_w,  1032312_w, 103231_w, 2101323_w, 210132_w,
    //            2103123_w, 210312_w, 2130123_w, 213012_w}));
    // }
    KnuthBendix<TestType> kb(twosided, p);
    kb.run();
    REQUIRE(kb.number_of_classes() == 312);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "027",
                                   "Reinis MFE",
                                   "[todd-coxeter][quick]",
                                   KNUTH_BENDIX_TYPES) {
    using literals::        operator""_w;
    Presentation<word_type> p;
    p.alphabet(2);
    presentation::add_rule(p, "000"_w, "11"_w);
    presentation::add_rule(p, "001"_w, "10"_w);
    KnuthBendix<> kb(congruence_kind::twosided, p);

    REQUIRE(knuth_bendix::contains(kb, "000"_w, "11"_w));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("KnuthBendix",
                                   "068",
                                   "sigma sylvester monoid x 2",
                                   "[todd-coxeter][quick]",
                                   KNUTH_BENDIX_TYPES) {
    using namespace literals;
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(3);
    p.contains_empty_word(true);
    presentation::add_rule(p, 00_w, 0_w);
    presentation::add_rule(p, 11_w, 1_w);
    presentation::add_rule(p, 22_w, 2_w);
    presentation::add_rule(p, 010_w, 01_w);
    presentation::add_rule(p, 0120_w, 012_w);
    presentation::add_rule(p, 020_w, 02_w);
    presentation::add_rule(p, 0210_w, 021_w);
    presentation::add_rule(p, 02120_w, 0212_w);
    presentation::add_rule(p, 10212_w, 1021_w);
    presentation::add_rule(p, 121_w, 12_w);
    presentation::add_rule(p, 12012_w, 1201_w);
    presentation::add_rule(p, 1202_w, 120_w);

    p.rules.clear();
    p.alphabet(2);
    presentation::add_idempotent_rules_no_checks(p, 01_w);
    using words::operator+;
    WordRange    words;
    words.alphabet_size(2).min(0).max(3);
    size_t n = 2;
    for (size_t a = 0; a < n - 1; ++a) {
      for (size_t b = a; b < n - 1; ++b) {
        for (size_t c = b + 1; c < n; ++c) {
          for (auto& u : words) {
            for (auto& v : words) {
              for (auto& w : words) {
                presentation::add_rule(
                    p, u + a + c + v + b + w, u + c + a + v + b + w);
              }
            }
          }
        }
      }
    }
    presentation::sort_each_rule(p);
    presentation::sort_rules(p);
    presentation::remove_trivial_rules(p);
    p.contains_empty_word(true);

    KnuthBendix<TestType> kb(twosided, p);
    p = to_presentation<word_type>(kb);

    auto S = to_froidure_pin(kb);
    REQUIRE(S.contains_one());
    REQUIRE(S.size() == kb.number_of_classes());
    REQUIRE(S.number_of_idempotents() == 5);
    REQUIRE(kb.number_of_classes() == 6);
  }
}  // namespace libsemigroups
