//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2022-2023 James D. Mitchell
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

#define CATCH_CONFIG_ENABLE_ALL_STRINGMAKERS
// TODO use _w literals

#include <cstddef>   // for size_t
#include <iostream>  // for cout

#include <algorithm>    // for for_each
#include <chrono>       // for duration, microsec...
#include <cstddef>      // for size_t
#include <cstdint>      // for size_t
#include <fstream>      // for ofstream
#include <iostream>     // for string, char_traits
#include <iterator>     // for advance
#include <string>       // for swap, basic_string
#include <type_traits>  // for remove_reference<>...
#include <utility>      // for move
#include <vector>       // for vector, operator==

#include "catch_amalgamated.hpp"  // for REQUIRE, REQUIRE_THROWS_AS, REQUI...
#include "test-main.hpp"          // for LIBSEMIGROUPS_TEST_CASE

#include "libsemigroups/constants.hpp"        // for UNDEFINED
#include "libsemigroups/exception.hpp"        // for LibsemigroupsExcep...
#include "libsemigroups/fpsemi-examples.hpp"  // for make, fibonacci_se...
#include "libsemigroups/order.hpp"            // for LexicographicalCom...
#include "libsemigroups/presentation.hpp"     // for add_rule
#include "libsemigroups/stephen.hpp"          // for Stephen, Stephen::...
#include "libsemigroups/todd-coxeter.hpp"     // for ToddCoxeter
#include "libsemigroups/types.hpp"            // for word_type
#include "libsemigroups/word-graph.hpp"       // for WordGraph, ope...
#include "libsemigroups/word-range.hpp"       // for StringToWord, word...

#include "libsemigroups/detail/report.hpp"  // for ReportGuard

namespace libsemigroups {
  using namespace literals;
  using rx::operator|;

  namespace {
    template <typename PresentationType>
    void check_000(Stephen<PresentationType>& s) {
      s.set_word({0}).run();
      REQUIRE(s.word_graph().number_of_nodes() == 2);
      REQUIRE(
          s.word_graph()
          == make<WordGraph<uint32_t>>(2, {{1, UNDEFINED}, {UNDEFINED, 1}}));
      REQUIRE(stephen::number_of_words_accepted(s) == POSITIVE_INFINITY);
      {
        REQUIRE((stephen::words_accepted(s) | rx::take(10) | rx::to_vector())
                == std::vector<word_type>({0_w,
                                           01_w,
                                           011_w,
                                           0111_w,
                                           01111_w,
                                           011111_w,
                                           0111111_w,
                                           01111111_w,
                                           011111111_w,
                                           0111111111_w}));
      }
      {
        REQUIRE((stephen::left_factors(s) | rx::take(10) | rx::to_vector())
                == std::vector<word_type>({{},
                                           0_w,
                                           01_w,
                                           011_w,
                                           0111_w,
                                           01111_w,
                                           011111_w,
                                           0111111_w,
                                           01111111_w,
                                           011111111_w}));
      }
    }

    void verify_c4_normal_form(Presentation<std::string> const& p,
                               std::string const&               word,
                               std::string const&               nf) {
      ToWord  to_word(p.alphabet());
      Stephen S(p);
      S.set_word(to_word(word)).run();

      REQUIRE((stephen::words_accepted(S) | ToString(p.alphabet())
               | rx::sort(LexicographicalCompare()) | rx::take(1))
                  .get()
              == nf);

      REQUIRE((stephen::words_accepted(S) | rx::all_of([&S](auto const& w) {
                 return stephen::accepts(S, w);
               })));
      REQUIRE(stephen::number_of_words_accepted(S)
              == stephen::words_accepted(S).count());
    }

    void verify_c4_equal_to(Presentation<std::string> const& p,
                            std::string const&               word1,
                            std::string const&               word2) {
      ToWord  to_word(p.alphabet());
      Stephen S(p);
      S.set_word(to_word(word1)).run();
      REQUIRE(stephen::accepts(S, to_word(word2)));
      S.set_word(to_word(word2)).run();
      REQUIRE(stephen::accepts(S, to_word(word1)));
    }

    void verify_c4_not_equal_to(Presentation<std::string> const& p,
                                std::string const&               word1,
                                std::string const&               word2) {
      ToWord  to_word(p.alphabet());
      Stephen S(p);
      S.set_word(to_word(word1)).run();
      REQUIRE(!stephen::accepts(S, to_word(word2)));
      S.set_word(to_word(word2)).run();
      REQUIRE(!stephen::accepts(S, to_word(word1)));
    }

  }  // namespace

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "000",
                          "basic test 1",
                          "[quick][stephen]") {
    auto                    rg = ReportGuard(true);
    Presentation<word_type> p;
    p.alphabet(2);
    presentation::add_rule(p, {0}, {0, 1});
    Stephen s(p);
    check_000(s);
    s.init(p);
    check_000(s);
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "001",
                          "basic test 2",
                          "[quick][stephen]") {
    auto                    rg = ReportGuard(true);
    Presentation<word_type> p;
    p.alphabet(2);
    presentation::add_rule(p, 000_w, 0_w);
    presentation::add_rule(p, 111_w, 1_w);
    presentation::add_rule(p, 0101_w, 00_w);
    Stephen s(p);
    s.set_word(1101_w).run();
    REQUIRE(s.word_graph().number_of_nodes() == 7);
    REQUIRE(s.word_graph()
            == make<WordGraph<uint32_t>>(7,
                                         {{UNDEFINED, 1},
                                          {UNDEFINED, 2},
                                          {3, 1},
                                          {4, 5},
                                          {3, 6},
                                          {6, 3},
                                          {5, 4}}));
    REQUIRE(stephen::number_of_words_accepted(s) == POSITIVE_INFINITY);

    word_type w = 1101_w;

    REQUIRE(word_graph::last_node_on_path_no_checks(
                s.word_graph(), 0, w.begin(), w.end())
                .first
            == 5);
    w = 110010_w;
    REQUIRE(word_graph::last_node_on_path_no_checks(
                s.word_graph(), 0, w.begin(), w.end())
                .first
            == 5);

    REQUIRE(stephen::accepts(s, 110010_w));
    REQUIRE(stephen::accepts(s, 110010_w));
    REQUIRE(!stephen::accepts(s, {}));
    REQUIRE(!stephen::accepts(s, 0000000000_w));
    REQUIRE(!stephen::accepts(s, 111_w));
    {
      REQUIRE((stephen::words_accepted(s) | rx::take(10) | rx::to_vector())
              == std::vector<word_type>({1101_w,
                                         110001_w,
                                         110010_w,
                                         110100_w,
                                         110111_w,
                                         111101_w,
                                         11000001_w,
                                         11000010_w,
                                         11000100_w,
                                         11000111_w}));
    }
    {
      REQUIRE((stephen::left_factors(s) | rx::take(10) | rx::to_vector())
              == std::vector<word_type>({{},
                                         1_w,
                                         11_w,
                                         110_w,
                                         111_w,
                                         1100_w,
                                         1101_w,
                                         1111_w,
                                         11000_w,
                                         11001_w}));
      REQUIRE(stephen::number_of_left_factors(s) == POSITIVE_INFINITY);
      REQUIRE((stephen::left_factors(s) | rx::take(10)
               | rx::all_of([&s](auto const& ww) {
                   return stephen::is_left_factor(s, ww);
                 })));
    }

    s.set_word(00_w).run();
    REQUIRE(s.word_graph().number_of_nodes() == 5);
    REQUIRE(s.word_graph()
            == make<WordGraph<uint32_t>>(
                5, {{1, UNDEFINED}, {2, 3}, {1, 4}, {4, 1}, {3, 2}}));

    p.rules.clear();
    presentation::add_rule(p, 000_w, 0_w);
    presentation::add_rule(p, 111_w, 1_w);
    s.init(p).set_word(00_w).run();
    REQUIRE(s.word_graph().number_of_nodes() == 3);
    REQUIRE(s.word_graph()
            == make<WordGraph<uint32_t>>(
                3, {{1, UNDEFINED}, {2, UNDEFINED}, {1, UNDEFINED}}));
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "002",
                          "full transf monoid",
                          "[quick][stephen][no-valgrind]") {
    using namespace fpsemigroup;
    auto   rg = ReportGuard(true);
    size_t n  = 5;
    auto   p  = full_transformation_monoid(n, author::Iwahori);
    p.validate();

    Stephen s(std::move(p));
    s.set_word(010111020120_w).run();
    REQUIRE(s.word_graph().number_of_nodes() == 120);
    REQUIRE(
        s.word_graph()
        == make<WordGraph<uint32_t>>(
            120,
            {{1, 2, 3, 4, UNDEFINED},        {0, 5, 6, 7, UNDEFINED},
             {8, 0, 9, 10, UNDEFINED},       {11, 12, 0, 13, UNDEFINED},
             {14, 15, 16, 0, UNDEFINED},     {17, 1, 18, 19, UNDEFINED},
             {20, 21, 1, 22, UNDEFINED},     {23, 24, 25, 1, UNDEFINED},
             {2, 17, 26, 27, UNDEFINED},     {28, 29, 2, 30, UNDEFINED},
             {31, 32, 33, 2, UNDEFINED},     {3, 34, 20, 35, UNDEFINED},
             {36, 3, 29, 37, UNDEFINED},     {38, 39, 40, 3, UNDEFINED},
             {4, 41, 42, 23, UNDEFINED},     {43, 4, 44, 32, UNDEFINED},
             {45, 46, 4, 40, UNDEFINED},     {5, 8, 47, 48, UNDEFINED},
             {49, 50, 5, 51, UNDEFINED},     {52, 53, 54, 5, UNDEFINED},
             {6, 55, 11, 56, UNDEFINED},     {57, 6, 50, 58, UNDEFINED},
             {59, 60, 61, 6, UNDEFINED},     {7, 62, 63, 14, UNDEFINED},
             {64, 7, 65, 53, UNDEFINED},     {66, 67, 7, 61, UNDEFINED},
             {55, 57, 8, 68, UNDEFINED},     {62, 64, 69, 8, UNDEFINED},
             {9, 49, 55, 70, UNDEFINED},     {50, 9, 12, 71, UNDEFINED},
             {72, 73, 74, 9, UNDEFINED},     {10, 52, 75, 62, UNDEFINED},
             {53, 10, 76, 15, UNDEFINED},    {77, 78, 10, 74, UNDEFINED},
             {47, 11, 49, 79, UNDEFINED},    {63, 80, 66, 11, UNDEFINED},
             {12, 47, 57, 81, UNDEFINED},    {82, 76, 78, 12, UNDEFINED},
             {13, 83, 59, 63, UNDEFINED},    {84, 13, 73, 76, UNDEFINED},
             {61, 74, 13, 16, UNDEFINED},    {48, 14, 85, 52, UNDEFINED},
             {56, 86, 14, 59, UNDEFINED},    {15, 48, 87, 64, UNDEFINED},
             {88, 71, 15, 73, UNDEFINED},    {16, 89, 56, 66, UNDEFINED},
             {90, 16, 71, 78, UNDEFINED},    {34, 36, 17, 91, UNDEFINED},
             {41, 43, 92, 17, UNDEFINED},    {18, 28, 34, 93, UNDEFINED},
             {29, 18, 21, 94, UNDEFINED},    {95, 96, 97, 18, UNDEFINED},
             {19, 31, 98, 41, UNDEFINED},    {32, 19, 99, 24, UNDEFINED},
             {100, 101, 19, 97, UNDEFINED},  {26, 20, 28, 102, UNDEFINED},
             {42, 103, 45, 20, UNDEFINED},   {21, 26, 36, 104, UNDEFINED},
             {105, 99, 101, 21, UNDEFINED},  {22, 106, 38, 42, UNDEFINED},
             {107, 22, 96, 99, UNDEFINED},   {40, 97, 22, 25, UNDEFINED},
             {27, 23, 108, 31, UNDEFINED},   {35, 109, 23, 38, UNDEFINED},
             {24, 27, 110, 43, UNDEFINED},   {111, 94, 24, 96, UNDEFINED},
             {25, 112, 35, 45, UNDEFINED},   {113, 25, 94, 101, UNDEFINED},
             {106, 107, 114, 26, UNDEFINED}, {112, 113, 27, 114, UNDEFINED},
             {108, 111, 112, 28, UNDEFINED}, {115, 44, 46, 29, UNDEFINED},
             {30, 95, 106, 108, UNDEFINED},  {96, 30, 39, 44, UNDEFINED},
             {114, 40, 30, 33, UNDEFINED},   {102, 105, 31, 106, UNDEFINED},
             {116, 37, 32, 39, UNDEFINED},   {33, 100, 102, 112, UNDEFINED},
             {101, 33, 37, 46, UNDEFINED},   {98, 116, 100, 34, UNDEFINED},
             {110, 35, 111, 116, UNDEFINED}, {109, 110, 113, 36, UNDEFINED},
             {37, 98, 105, 109, UNDEFINED},  {91, 38, 95, 98, UNDEFINED},
             {39, 91, 107, 110, UNDEFINED},  {93, 115, 41, 95, UNDEFINED},
             {104, 42, 115, 105, UNDEFINED}, {103, 104, 43, 107, UNDEFINED},
             {44, 93, 103, 111, UNDEFINED},  {92, 45, 93, 100, UNDEFINED},
             {46, 92, 104, 113, UNDEFINED},  {83, 84, 117, 47, UNDEFINED},
             {89, 90, 48, 117, UNDEFINED},   {85, 88, 89, 49, UNDEFINED},
             {118, 65, 67, 50, UNDEFINED},   {51, 72, 83, 85, UNDEFINED},
             {73, 51, 60, 65, UNDEFINED},    {117, 61, 51, 54, UNDEFINED},
             {79, 82, 52, 83, UNDEFINED},    {119, 58, 53, 60, UNDEFINED},
             {54, 77, 79, 89, UNDEFINED},    {78, 54, 58, 67, UNDEFINED},
             {75, 119, 77, 55, UNDEFINED},   {87, 56, 88, 119, UNDEFINED},
             {86, 87, 90, 57, UNDEFINED},    {58, 75, 82, 86, UNDEFINED},
             {68, 59, 72, 75, UNDEFINED},    {60, 68, 84, 87, UNDEFINED},
             {70, 118, 62, 72, UNDEFINED},   {81, 63, 118, 82, UNDEFINED},
             {80, 81, 64, 84, UNDEFINED},    {65, 70, 80, 88, UNDEFINED},
             {69, 66, 70, 77, UNDEFINED},    {67, 69, 81, 90, UNDEFINED},
             {74, 117, 68, 69, UNDEFINED},   {71, 85, 86, 118, UNDEFINED},
             {76, 79, 119, 80, UNDEFINED},   {97, 114, 91, 92, UNDEFINED},
             {94, 108, 109, 115, UNDEFINED}, {99, 102, 116, 103, UNDEFINED}}));
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "003",
                          "from step_hen 002",
                          "[quick][stephen]") {
    Presentation<std::string> p;
    p.alphabet("ab");
    presentation::add_rule(p, "aaa", "a");
    presentation::add_rule(p, "bbb", "b");
    presentation::add_rule(p, "abab", "aa");

    Stephen S(p);
    S.set_word("bbab"_w);

    REQUIRE(stephen::accepts(S, "bbaaba"_w));
    REQUIRE(!stephen::accepts(S, ""_w));
    REQUIRE(!stephen::accepts(S, "aaaaaaaaaa"_w));
    REQUIRE(!stephen::accepts(S, "bbb"_w));

    S.set_word("bba"_w);
    REQUIRE(stephen::accepts(S, "bbabb"_w));
    REQUIRE(stephen::accepts(S, "bba"_w));
    REQUIRE(!stephen::accepts(S, "bbb"_w));
    REQUIRE(!stephen::accepts(S, "a"_w));
    REQUIRE(!stephen::accepts(S, "ab"_w));

    S.set_word("bbaab"_w);
    REQUIRE(stephen::accepts(S, "bbaba"_w));
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "004",
                          "from step_hen 003",
                          "[quick][stephen]") {
    Presentation<std::string> p;
    p.alphabet("abcdefg");
    presentation::add_rule(p, "aaaeaa", "abcd");
    presentation::add_rule(p, "ef", "dg");

    Stephen S(p);

    S.set_word("abcef"_w).run();
    REQUIRE("abcef"_w == 01245_w);
    REQUIRE(S.word_graph()
            == make<WordGraph<uint32_t>>(
                11,
                {{1,
                  UNDEFINED,
                  UNDEFINED,
                  UNDEFINED,
                  UNDEFINED,
                  UNDEFINED,
                  UNDEFINED},
                 {2, 3},
                 {4},
                 {UNDEFINED, UNDEFINED, 5},
                 {UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED, 6},
                 {UNDEFINED, UNDEFINED, UNDEFINED, 7, 8},
                 {9},
                 {UNDEFINED,
                  UNDEFINED,
                  UNDEFINED,
                  UNDEFINED,
                  UNDEFINED,
                  UNDEFINED,
                  10},
                 {UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED, 10},
                 {7}}));
    ToWord to_word;
    auto   rule = to_word(p.rules[0]);
    auto   m    = word_graph::last_node_on_path(
                 S.word_graph(), 0, rule.cbegin(), rule.cend())
                 .first;
    rule   = to_word(p.rules[1]);
    auto n = word_graph::last_node_on_path(
                 S.word_graph(), 0, rule.cbegin(), rule.cend())
                 .first;
    REQUIRE(m != UNDEFINED);
    REQUIRE(n != UNDEFINED);
    REQUIRE(m == n);
    REQUIRE(S.word_graph().number_of_nodes() == 11);
    REQUIRE(stephen::accepts(S, "aaaeaag"_w));
    REQUIRE(stephen::number_of_words_accepted(S) == 3);
    REQUIRE((stephen::words_accepted(S) | rx::to_vector())
            == std::vector<word_type>({01236_w, 01245_w, 0004006_w}));

    S.set_word("aaaeaaaeaa"_w).run();

    REQUIRE(S.word_graph().number_of_nodes() == 15);
    REQUIRE(stephen::number_of_words_accepted(S) == 3);
    REQUIRE(stephen::accepts(S, "aaaeabcd"_w));
    REQUIRE(std::vector<word_type>(stephen::cbegin_words_accepted(S),
                                   stephen::cend_words_accepted(S))
            == std::vector<word_type>({00040123_w, 01230400_w, 0004000400_w}));

    S.set_word("aaaeaag"_w).run();
    REQUIRE(S.word_graph().number_of_nodes() == 11);
    REQUIRE(stephen::number_of_words_accepted(S) == 3);
    REQUIRE(std::vector<word_type>(stephen::cbegin_words_accepted(S),
                                   stephen::cend_words_accepted(S))
            == std::vector<word_type>({01236_w, 01245_w, 0004006_w}));
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "005",
                          "from step_hen 004",
                          "[quick][stephen]") {
    Presentation<std::string> p;
    p.alphabet("abc");
    presentation::add_rule(p, "ab", "ba");
    presentation::add_rule(p, "ac", "cc");
    presentation::add_rule(p, "ac", "a");
    presentation::add_rule(p, "cc", "a");
    presentation::add_rule(p, "bc", "cc");
    presentation::add_rule(p, "bcc", "b");
    presentation::add_rule(p, "bc", "b");
    presentation::add_rule(p, "cc", "b");
    presentation::add_rule(p, "a", "b");

    Stephen S(p);
    S.set_word("abcc"_w).run();
    REQUIRE(stephen::accepts(S, "baac"_w));
    REQUIRE(S.word_graph().number_of_nodes() == 3);
    REQUIRE(stephen::number_of_words_accepted(S) == POSITIVE_INFINITY);
    REQUIRE((stephen::words_accepted(S) | rx::take(10) | ToString("abc")
             | rx::to_vector())
            == std::vector<std::string>(
                {"a", "b", "aa", "ab", "ac", "ba", "bb", "bc", "ca", "cb"}));
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "006",
                          "from step_hen 005",
                          "[quick][stephen]") {
    Presentation<std::string> p;
    p.alphabet("abcd");
    presentation::add_rule(p, "bb", "c");
    presentation::add_rule(p, "caca", "abab");
    presentation::add_rule(p, "bc", "d");
    presentation::add_rule(p, "cb", "d");
    presentation::add_rule(p, "aa", "d");
    presentation::add_rule(p, "ad", "a");
    presentation::add_rule(p, "da", "a");
    presentation::add_rule(p, "bd", "b");
    presentation::add_rule(p, "db", "b");
    presentation::add_rule(p, "cd", "c");
    presentation::add_rule(p, "dc", "c");

    Stephen S(p);
    S.set_word("dabdaaadabab"_w).run();
    REQUIRE(stephen::accepts(S, "abdadcaca"_w));
    REQUIRE(S.word_graph().number_of_nodes() == 25);
    REQUIRE(stephen::number_of_words_accepted(S) == POSITIVE_INFINITY);
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "007",
                          "Fibonacci(4, 6)",
                          "[stephen][extreme]") {
    using namespace fpsemigroup;
    auto    rg = ReportGuard(true);
    Stephen S(fibonacci_semigroup(4, 6));
    S.set_word(0123_w).run_for(std::chrono::seconds(10));
    REQUIRE(!S.finished());
  }

  LIBSEMIGROUPS_TEST_CASE(
      "Stephen",
      "008",
      "C(4) monoid normal form (test_case_knuth_bendix_055)",
      "[stephen][quick]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abcdefg");
    presentation::add_rule(p, "abcd", "ce");
    presentation::add_rule(p, "df", "dg");

    Stephen S(p);
    S.set_word("dfabcdf"_w).run();

    REQUIRE(S.word_graph().number_of_nodes() == 9);
    REQUIRE(stephen::number_of_words_accepted(S) == 8);

    REQUIRE(
        (stephen::words_accepted(S) | ToString(p.alphabet()) | rx::to_vector())
        == std::vector<std::string>({"dfcef",
                                     "dfceg",
                                     "dgcef",
                                     "dgceg",
                                     "dfabcdf",
                                     "dfabcdg",
                                     "dgabcdf",
                                     "dgabcdg"}));

    REQUIRE((stephen::words_accepted(S) | rx::sort(LexicographicalCompare())
             | rx::take(1) | ToString(p.alphabet()))
                .get()
            == "dfabcdf");

    REQUIRE((stephen::words_accepted(S) | rx::all_of([&S](auto const& w) {
               return stephen::accepts(S, w);
             })));
    REQUIRE(stephen::number_of_words_accepted(S)
            == stephen::words_accepted(S).count());

    S.set_word("abcdfceg"_w).run();
    REQUIRE(stephen::number_of_words_accepted(S) == 16);

    REQUIRE((stephen::words_accepted(S) | ToString(p.alphabet())
             | rx::sort(LexicographicalCompare()) | rx::to_vector())
            == std::vector<std::string>({"abcdfabcdf",
                                         "abcdfabcdg",
                                         "abcdfcef",
                                         "abcdfceg",
                                         "abcdgabcdf",
                                         "abcdgabcdg",
                                         "abcdgcef",
                                         "abcdgceg",
                                         "cefabcdf",
                                         "cefabcdg",
                                         "cefcef",
                                         "cefceg",
                                         "cegabcdf",
                                         "cegabcdg",
                                         "cegcef",
                                         "cegceg"}));

    REQUIRE(ToString(p.alphabet())(stephen::words_accepted(S)
                                   | rx::sort(LexicographicalCompare())
                                   | rx::take(1))
                .get()
            == "abcdfabcdf");
    REQUIRE(stephen::accepts(S, "abcdfabcdf"_w));
  }

  LIBSEMIGROUPS_TEST_CASE(
      "Stephen",
      "009",
      "C(4) monoid normal form (test_case_gap_smalloverlap_85)",
      "[stephen][quick]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("cab");
    presentation::add_rule(p, "aabc", "acba");

    ToWord  word("cab");
    Stephen S(p);
    S.set_word(word("a")).run();
    REQUIRE(!stephen::accepts(S, word("b")));

    S.set_word(word("aabcabc")).run();
    REQUIRE(stephen::accepts(S, word("aabccba")));

    S.set_word(word("aabccba")).run();
    REQUIRE(stephen::accepts(S, word("aabcabc")));

    S.set_word(word("acba")).run();
    auto words = stephen::words_accepted(S) | ToString(p.alphabet());

    REQUIRE((words | rx::to_vector())
            == std::vector<std::string>({"acba", "aabc"}));
    verify_c4_normal_form(p, "acba", "aabc");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "010",
                          "code coverage",
                          "[stephen][quick]") {
    auto                      rg = ReportGuard(true);
    Presentation<std::string> p;
    REQUIRE_THROWS_AS(Stephen(p), LibsemigroupsException);
    p.alphabet("abcdefg");

    Stephen S(p);
    auto    w = "abbbddbcbcbc"_w;
    S.set_word(w);
    S.run();
    REQUIRE(S.finished());
    S.run();
    S.set_word("abbbddbcbcbc"_w);  // resets
    S.report_every(std::chrono::microseconds(10));

    S.run();

    Stephen T(S);
    REQUIRE(stephen::accepts(T, w));
    REQUIRE(!stephen::accepts(T, "abbbd"_w));
    REQUIRE(stephen::number_of_words_accepted(T) == 1);
    REQUIRE(stephen::number_of_left_factors(T) == w.size() + 1);
    REQUIRE(stephen::accepts(S, w));
    REQUIRE(!stephen::accepts(S, "abbbd"_w));
    REQUIRE(stephen::number_of_words_accepted(S) == 1);
    REQUIRE(stephen::number_of_left_factors(S) == w.size() + 1);

    Stephen U(std::move(S));
    REQUIRE(stephen::accepts(U, w));
    REQUIRE(!stephen::accepts(U, "abbbd"_w));
    REQUIRE(stephen::number_of_words_accepted(U) == 1);
    REQUIRE(stephen::number_of_left_factors(U) == w.size() + 1);

    S = T;
    REQUIRE(stephen::accepts(T, w));
    REQUIRE(!stephen::accepts(T, "abbbd"_w));
    REQUIRE(stephen::number_of_words_accepted(T) == 1);
    REQUIRE(stephen::number_of_left_factors(T) == w.size() + 1);
    REQUIRE(stephen::accepts(S, w));
    REQUIRE(!stephen::accepts(S, "abbbd"_w));
    REQUIRE(stephen::number_of_words_accepted(S) == 1);
    REQUIRE(stephen::number_of_left_factors(S) == w.size() + 1);

    decltype(S) V;
    V = std::move(S);
    REQUIRE(stephen::accepts(V, w));
    REQUIRE(!stephen::accepts(V, "abbbd"_w));
    REQUIRE(stephen::number_of_words_accepted(V) == 1);
    REQUIRE(stephen::number_of_left_factors(V) == w.size() + 1);
    REQUIRE(V.word() == w);
    REQUIRE(V.accept_state() == 12);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "Stephen",
      "011",
      "C(4) monoid normal form (test_case_gap_smalloverlap_49)",
      "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcdefgh");

    presentation::add_rule(p, "abcd", "ce");
    presentation::add_rule(p, "df", "hd");
    verify_c4_equal_to(p, "abchd", "abcdf");
    verify_c4_equal_to(p, "abchd", "abchd");
    verify_c4_equal_to(p, "abchdf", "abchhd");
    verify_c4_equal_to(p, "abchd", "cef");
    verify_c4_equal_to(p, "cef", "abchd");
    verify_c4_not_equal_to(p, "abchf", "abcdf");

    verify_c4_equal_to(p, "hdfabce", "dffababcd");

    verify_c4_normal_form(p, "hdfabce", "dffababcd");
  }

  LIBSEMIGROUPS_TEST_CASE(
      "Stephen",
      "012",
      "C(4) monoid normal form (test_case_gap_smalloverlap_63)",
      "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcdefgh");

    presentation::add_rule(p, "afh", "bgh");
    presentation::add_rule(p, "hc", "d");
    verify_c4_equal_to(p, "afd", "bgd");
    verify_c4_equal_to(p, "bghcafhbgd", "afdafhafd");
    verify_c4_normal_form(p, "bghcafhbgd", "afdafhafd");
  }

  LIBSEMIGROUPS_TEST_CASE(
      "Stephen",
      "013",
      "C(4) monoid equal to (test_case_gap_smalloverlap_70)",
      "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcdefghij");

    presentation::add_rule(p, "afh", "bgh");
    presentation::add_rule(p, "hc", "de");
    presentation::add_rule(p, "ei", "j");

    verify_c4_equal_to(p, "afdj", "bgdj");
    verify_c4_not_equal_to(p, "jjjjjjjjjjjjjjjjjjjjjjj", "b");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "014",
                          "C(4) monoid normal form (test_case_ex_3_13_14)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcd");
    presentation::add_rule(p, "abbba", "cdc");
    verify_c4_normal_form(p, "cdcdcabbbabbbabbcd", "abbbadcabbbabbbabbcd");
    verify_c4_equal_to(p, "cdcdcabbbabbbabbcd", "abbbadcabbbabbbabbcd");

    verify_c4_equal_to(p, "abbbadcbbba", "cdabbbcdc");
    verify_c4_equal_to(p, "cdabbbcdc", "cdabbbcdc");
    verify_c4_normal_form(p, "cdabbbcdc", "abbbadcbbba");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "015",
                          "C(4) monoid normal form (test_case_ex_3_15)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcd");
    presentation::add_rule(p, "aabc", "acba");
    std::string original = "cbacbaabcaabcacbacba";
    std::string expected = "cbaabcabcaabcaabcabc";

    verify_c4_equal_to(p, "cbaabcabcaabcaabccba", original);
    verify_c4_equal_to(p, original, expected);
    verify_c4_equal_to(p, expected, original);
    verify_c4_equal_to(p, "cbaabcabcaabcaabccba", expected);

    verify_c4_equal_to(p, original, "cbaabcabcaabcaabccba");

    verify_c4_equal_to(p, expected, "cbaabcabcaabcaabccba");
    verify_c4_normal_form(p, original, expected);
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "016",
                          "C(4) monoid normal form (test_case_ex_3_16)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcd");
    presentation::add_rule(p, "abcd", "acca");
    std::string original = "bbcabcdaccaccabcddd";
    std::string expected = "bbcabcdabcdbcdbcddd";

    verify_c4_equal_to(p, original, expected);
    verify_c4_equal_to(p, expected, original);

    verify_c4_normal_form(p, original, expected);
    verify_c4_normal_form(p, expected, expected);
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "017",
                          "C(4) monoid normal form (test_case_mt_3)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcd");
    presentation::add_rule(p, "abcd", "accca");

    verify_c4_normal_form(p, "bbcabcdaccaccabcddd", "bbcabcdaccaccabcddd");
    verify_c4_equal_to(p, "bbcabcdaccaccabcddd", "bbcabcdaccaccabcddd");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "018",
                          "C(4) monoid normal form (test_case_mt_5)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abc");
    presentation::add_rule(p, "ac", "cbbbbc");

    verify_c4_normal_form(p, "acbbbbc", "aac");
    verify_c4_equal_to(p, "acbbbbc", "aac");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "019",
                          "C(4) monoid normal form (test_case_mt_6)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abc");
    presentation::add_rule(p, "ccab", "cbac");

    verify_c4_normal_form(p, "bacbaccabccabcbacbac", "bacbacbaccbaccbacbac");
    verify_c4_equal_to(p, "bacbaccabccabcbacbac", "bacbacbaccbaccbacbac");
    verify_c4_normal_form(p, "ccabcbaccab", "cbaccbacbac");
    verify_c4_equal_to(p, "ccabcbaccab", "cbaccbacbac");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "020",
                          "C(4) monoid normal form (test_case_mt_10)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcdefghij");
    presentation::add_rule(p, "afh", "bgh");
    presentation::add_rule(p, "hc", "de");
    presentation::add_rule(p, "ei", "j");

    verify_c4_normal_form(p, "bgdj", "afdei");
    verify_c4_equal_to(p, "bgdj", "afdei");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "021",
                          "C(4) monoid normal form (test_case_mt_13)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcd");
    presentation::add_rule(p, "abcd", "dcba");

    verify_c4_normal_form(p, "dcbdcba", "abcdbcd");
    verify_c4_equal_to(p, "dcbdcba", "abcdbcd");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "022",
                          "C(4) monoid normal form (test_case_mt_14)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcd");
    presentation::add_rule(p, "abca", "dcbd");

    verify_c4_normal_form(p, "dcbabca", "abcacbd");
    verify_c4_equal_to(p, "dcbabca", "abcacbd");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "023",
                          "C(4) monoid normal form (test_case_mt_15)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcd");
    presentation::add_rule(p, "abcd", "dcba");
    presentation::add_rule(p, "adda", "dbbd");

    verify_c4_normal_form(p, "dbbabcd", "addacba");
    verify_c4_equal_to(p, "dbbabcd", "addacba");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "024",
                          "C(4) monoid normal form (test_case_mt_16)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcdefg");
    presentation::add_rule(p, "abcd", "acca");
    presentation::add_rule(p, "gf", "ge");

    verify_c4_normal_form(p, "accabcdgf", "abcdbcdge");
    verify_c4_equal_to(p, "accabcdgf", "abcdbcdge");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "025",
                          "C(4) monoid normal form (test_case_mt_17)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcd");
    presentation::add_rule(
        p, "ababbabbbabbbb", "abbbbbabbbbbbabbbbbbbabbbbbbbb");
    presentation::add_rule(
        p, "cdcddcdddcdddd", "cdddddcddddddcdddddddcdddddddd");

    verify_c4_normal_form(
        p, "abbbacdddddcddddddcdddddddcdddddddd", "abbbacdcddcdddcdddd");
    verify_c4_equal_to(
        p, "abbbacdddddcddddddcdddddddcdddddddd", "abbbacdcddcdddcdddd");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "026",
                          "C(4) monoid normal form (test_case_weak_1)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcd");
    presentation::add_rule(p, "acba", "aabc");
    presentation::add_rule(p, "acba", "dbbbd");

    verify_c4_equal_to(p, "aaabc", "adbbbd");
    verify_c4_equal_to(p, "adbbbd", "aaabc");

    verify_c4_equal_to(p, "aaabcadbbbd", "adbbbdadbbbd");
    verify_c4_equal_to(p, "aaabcaaabc", "adbbbdadbbbd");
    verify_c4_equal_to(p, "acba", "dbbbd");
    verify_c4_equal_to(p, "acbabbbd", "aabcbbbd");
    verify_c4_equal_to(p, "aabcbbbd", "acbabbbd");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "027",
                          "C(4) monoid normal form (test_case_weak_2)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcd");
    presentation::add_rule(p, "acba", "aabc");
    presentation::add_rule(p, "acba", "adbd");
    verify_c4_equal_to(p, "acbacba", "aabcabc");
    verify_c4_normal_form(p, "acbacba", "aabcabc");
    verify_c4_equal_to(p, "aabcabc", "acbacba");
    verify_c4_normal_form(p, "aabcabc", "aabcabc");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "028",
                          "C(4) monoid normal form (test_case_weak_3)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcde");
    presentation::add_rule(p, "bceac", "aeebbc");
    presentation::add_rule(p, "aeebbc", "dabcd");
    verify_c4_normal_form(p, "bceacdabcd", "aeebbcaeebbc");
    verify_c4_normal_form(p, "aeebbcaeebbc", "aeebbcaeebbc");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "029",
                          "C(4) monoid normal form (test_case_weak_4)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcd");
    presentation::add_rule(p, "acba", "aabc");
    presentation::add_rule(p, "acba", "dbbd");
    verify_c4_normal_form(p, "bbacbcaaabcbbd", "bbacbcaaabcbbd");
    verify_c4_normal_form(p, "acbacba", "aabcabc");
    verify_c4_normal_form(p, "aabcabc", "aabcabc");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "030",
                          "C(4) monoid normal form (test_case_weak_5)",
                          "[stephen][quick]") {
    Presentation<std::string> p;
    p.alphabet("abcd");
    presentation::add_rule(p, "acba", "aabc");
    presentation::add_rule(p, "acba", "adbd");
    verify_c4_normal_form(p, "acbacba", "aabcabc");
    verify_c4_normal_form(p, "aabcabc", "aabcabc");
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "031",
                          "Test behaviour when uninitialised",
                          "[stephen][quick]") {
    Stephen<Presentation<word_type>> S;

    REQUIRE_THROWS_AS(S.accept_state(), LibsemigroupsException);
    REQUIRE_THROWS_AS(stephen::cbegin_words_accepted(S),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(stephen::cend_words_accepted(S), LibsemigroupsException);
    REQUIRE_THROWS_AS(stephen::cbegin_left_factors(S), LibsemigroupsException);
    REQUIRE_THROWS_AS(stephen::cend_left_factors(S), LibsemigroupsException);

    REQUIRE_THROWS_AS(stephen::is_left_factor(S, {0, 0, 0}),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(stephen::accepts(S, {0, 0, 0}), LibsemigroupsException);
    REQUIRE_THROWS_AS(stephen::number_of_left_factors(S),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(stephen::number_of_words_accepted(S),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(S.run(), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "034",
                          "(inverse presentation) "
                          "step_hen/tests/test_schutzenbergergraph.py:test_001",
                          "[stephen][quick]") {
    ToWord to_word("abcABC");

    InversePresentation<word_type> p;
    p.alphabet(to_word("abcABC"));
    p.inverses_no_checks(to_word("ABCabc"));

    auto S = Stephen(std::move(p));

    S.set_word(to_word("aBcAbC")).run();

    REQUIRE(S.finished());
    REQUIRE(S.word_graph().number_of_nodes() == 7);
    REQUIRE(!stephen::accepts(S, to_word("BaAbaBcAbC")));
    REQUIRE(stephen::accepts(S, to_word("aBcCbBcAbC")));

    S.set_word(to_word("aBcCbBcAbC"));
    REQUIRE(stephen::accepts(S, to_word("aBcAbC")));

    S.set_word(to_word("BaAbaBcAbC"));
    REQUIRE(stephen::accepts(S, to_word("aBcAbC")));
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "035",
                          "(inverse presentation) "
                          "step_hen/tests/test_schutzenbergergraph.py:test_001",
                          "[stephen][quick]") {
    InversePresentation<std::string> p;
    p.alphabet("abcABC");
    p.inverses_no_checks("ABCabc");
    ToWord to_word(p.alphabet());

    Stephen S(p);
    REQUIRE(to_word("aBcAbC") == 042315_w);

    S.set_word(to_word("aBcAbC")).run();

    REQUIRE(S.finished());
    REQUIRE(S.word_graph().number_of_nodes() == 7);
    REQUIRE(!stephen::accepts(S, to_word("BaAbaBcAbC")));
    REQUIRE(stephen::accepts(S, to_word("aBcCbBcAbC")));

    S.set_word(to_word("aBcCbBcAbC"));
    REQUIRE(stephen::accepts(S, to_word("aBcAbC")));

    S.set_word(to_word("BaAbaBcAbC"));
    REQUIRE(stephen::accepts(S, to_word("aBcAbC")));
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "036",
                          "(inverse presentation) "
                          "step_hen/tests/test_schutzenbergergraph.py:test_002",
                          "[stephen][quick]") {
    ToWord                         to_word("abcABC");
    InversePresentation<word_type> p;
    p.alphabet(to_word("abcABC"));
    p.inverses_no_checks(to_word("ABCabc"));

    auto S = Stephen(p);
    S.set_word(to_word("aBbcABAabCc")).run();

    REQUIRE(S.accept_state() == 4);
    REQUIRE(word_graph::follow_path(S.word_graph(), 0, S.word()) == 4);
    REQUIRE(stephen::number_of_words_accepted(S) == POSITIVE_INFINITY);
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "037",
                          "(inverse presentation) "
                          "step_hen/tests/test_schutzenbergergraph.py:test_003",
                          "[stephen][quick]") {
    ToWord                         to_word("xyXY");
    InversePresentation<word_type> p;
    p.alphabet(to_word("xyXY"));
    p.inverses_no_checks(to_word("XYxy"));

    auto S = Stephen(p);
    S.set_word(to_word("xxxyyy")).run();

    REQUIRE(stephen::accepts(S, to_word("xxxyyyYYYXXXxxxyyy")));
    S.set_word(to_word("xxxyyyYYYXXXxxxyyy"));
    REQUIRE(stephen::accepts(S, to_word("xxxyyy")));
    REQUIRE(!stephen::accepts(S, to_word("xxx")));
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "038",
                          "(inverse presentation) "
                          "step_hen/tests/test_schutzenbergergraph.py:test_004",
                          "[stephen][quick]") {
    ToWord                         to_word("xyXY");
    InversePresentation<word_type> p;
    p.alphabet(to_word("xyXY"));
    p.inverses_no_checks(to_word("XYxy"));
    presentation::add_rule(p, to_word("xyXxyX"), to_word("xyX"));

    auto S = Stephen(p);
    // FIXME seems like every word is accepted when we don't set the word
    S.set_word(to_word("xyXyy"));
    std::string ys = "";
    for (size_t i = 0; i < 10; ++i) {
      REQUIRE(stephen::accepts(S, to_word(std::string("x") + ys + "Xyy")));
      ys += "y";
    }

    REQUIRE(!stephen::accepts(S, to_word("xXyx")));
    REQUIRE(!stephen::accepts(S, to_word("xXxx")));
    REQUIRE(!stephen::accepts(S, to_word("xXxy")));
    REQUIRE(!stephen::accepts(S, to_word("xXxX")));
    REQUIRE(!stephen::accepts(S, to_word("xXyY")));
    REQUIRE(stephen::accepts(S, to_word("xyXyy")));
    REQUIRE(stephen::number_of_words_accepted(S) == POSITIVE_INFINITY);
    REQUIRE(S.word_graph().number_of_nodes() == 4);
    REQUIRE(S.word_graph().number_of_edges() == 8);

    REQUIRE(
        S.word_graph()
        == make<WordGraph<uint32_t>>(4,
                                     {{1, 2, UNDEFINED, UNDEFINED},
                                      {UNDEFINED, 1, 0, 1},
                                      {UNDEFINED, 3, UNDEFINED, 0},
                                      {UNDEFINED, UNDEFINED, UNDEFINED, 2}}));
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "039",
                          "(inverse presentation) "
                          "step_hen/tests/test_schutzenbergergraph.py:test_005",
                          "[stephen][quick]") {
    ToWord                         to_word("xyXY");
    InversePresentation<word_type> p;
    p.alphabet(to_word("xyXY"));
    p.inverses_no_checks(to_word("XYxy"));
    presentation::add_rule(p, to_word("xyXxyX"), to_word("xyX"));
    presentation::add_rule(p, to_word("xyxy"), to_word("xy"));

    auto S = Stephen(p);
    S.set_word(to_word("xyXyy"));
    REQUIRE(stephen::accepts(S, to_word("y")));
    REQUIRE(stephen::accepts(S, to_word("xxxxxxxxxxxxx")));
    REQUIRE(stephen::accepts(S, to_word("xyXxyxyxyxyxyXyy")));
    REQUIRE(S.word_graph().number_of_nodes() == 1);
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "040",
                          "(inverse presentation) "
                          "step_hen/tests/test_schutzenbergergraph.py:test_006",
                          "[stephen][quick]") {
    ToWord                         to_word("abcABC");
    InversePresentation<word_type> p;
    p.alphabet(to_word("abcABC"));
    p.inverses_no_checks(to_word("ABCabc"));
    presentation::add_rule(p, to_word("ac"), to_word("ca"));
    presentation::add_rule(p, to_word("ab"), to_word("ba"));
    presentation::add_rule(p, to_word("bc"), to_word("cb"));

    auto S = Stephen(p);
    S.set_word(to_word("BaAbaBcAbC"));
    S.run();
    REQUIRE(S.word_graph().number_of_nodes() == 7);
    REQUIRE(S.word_graph()
            == make<WordGraph<uint32_t>>(
                7,
                {{1, UNDEFINED, 2, UNDEFINED, 3, UNDEFINED},
                 {UNDEFINED, UNDEFINED, UNDEFINED, 0, 4, UNDEFINED},
                 {UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED, 5, 0},
                 {4, 0, 5, UNDEFINED, UNDEFINED, UNDEFINED},
                 {UNDEFINED, 1, 6, 3, UNDEFINED, UNDEFINED},
                 {6, 2, UNDEFINED, UNDEFINED, UNDEFINED, 3},
                 {UNDEFINED, UNDEFINED, UNDEFINED, 5, UNDEFINED, 4}}));
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen", "041", "corner case", "[stephen][quick]") {
    ToWord to_word("x");

    Presentation<word_type> p;
    p.contains_empty_word(true);

    p.alphabet(to_word("x"));
    presentation::add_rule(p, to_word("xxxx"), to_word("xx"));

    Stephen S(p);
    S.set_word(""_w);
    S.run();
    REQUIRE(S.accept_state() == 0);
    REQUIRE(S.word_graph().number_of_nodes() == 1);
    REQUIRE(!stephen::accepts(S, to_word("x")));
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen", "042", "empty word", "[stephen][quick]") {
    auto p = fpsemigroup::symmetric_inverse_monoid(4);
    REQUIRE(p.contains_empty_word());
    REQUIRE(p.alphabet().size() == 4);

    auto s = Stephen(p);
    s.set_word({}).run();
    REQUIRE(s.word_graph().number_of_nodes() == 24);
    s.set_word(0_w).run();
    REQUIRE(s.word_graph().number_of_nodes() == 24);
    s.set_word(1_w).run();
    REQUIRE(s.word_graph().number_of_nodes() == 24);
    s.set_word(2_w).run();
    REQUIRE(s.word_graph().number_of_nodes() == 24);
    s.set_word(3_w).run();
    REQUIRE(s.word_graph().number_of_nodes() == 48);
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen", "043", "shared_ptr", "[stephen][quick]") {
    ToWord                         to_word("abcABC");
    InversePresentation<word_type> p;
    p.alphabet(to_word("abcABC"));
    p.inverses_no_checks(to_word("ABCabc"));
    presentation::add_rule(p, to_word("ac"), to_word("ca"));
    presentation::add_rule(p, to_word("ab"), to_word("ba"));
    presentation::add_rule(p, to_word("bc"), to_word("cb"));

    auto ptr = std::make_shared<decltype(p)>(p);
    auto S   = Stephen<std::shared_ptr<InversePresentation<word_type>>>(ptr);
    S.set_word(to_word("BaAbaBcAbC"));
    S.run();
    REQUIRE(S.word_graph().number_of_nodes() == 7);
    REQUIRE(S.word_graph()
            == make<WordGraph<uint32_t>>(
                7,
                {{1, UNDEFINED, 2, UNDEFINED, 3, UNDEFINED},
                 {UNDEFINED, UNDEFINED, UNDEFINED, 0, 4, UNDEFINED},
                 {UNDEFINED, UNDEFINED, UNDEFINED, UNDEFINED, 5, 0},
                 {4, 0, 5, UNDEFINED, UNDEFINED, UNDEFINED},
                 {UNDEFINED, 1, 6, 3, UNDEFINED, UNDEFINED},
                 {6, 2, UNDEFINED, UNDEFINED, UNDEFINED, 3},
                 {UNDEFINED, UNDEFINED, UNDEFINED, 5, UNDEFINED, 4}}));
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "044",
                          "inverse presentation -- operator==",
                          "[stephen][quick][no-valgrind]") {
    ReportGuard rg(false);
    ToddCoxeter tc;
    {
      auto p = fpsemigroup::symmetric_inverse_monoid(4);
      REQUIRE(p.contains_empty_word());
      REQUIRE(p.alphabet().size() == 4);
      p.validate();
      tc.init(congruence_kind::twosided, p);
    }

    {
      InversePresentation<word_type> p
          = fpsemigroup::symmetric_inverse_monoid(4);
      p.inverses_no_checks(0123_w);
      p.validate();

      Stephen S(p);
      auto    w = 0120112010222201_w;
      S.set_word(w);

      REQUIRE(stephen::number_of_words_accepted(S) == POSITIVE_INFINITY);

      {
        auto const index = todd_coxeter::index_of(tc, w);
        auto       T     = Stephen(p);

        for (auto const& ww : (stephen::words_accepted(S) | rx::take(1024))) {
          REQUIRE(todd_coxeter::index_of(tc, ww) == index);
          REQUIRE(stephen::accepts(T.set_word(ww), ww));
        }
      }
    }
  }

  // This test fails in the fp-inverse-monoids branch and here, not sure if
  // it's supposed to or not should look at it more carefully
  //   LIBSEMIGROUPS_TEST_CASE("Stephen",
  //                           "045",
  //                           "inverse presentation",
  //                           "[stephen][quick]") {
  //     ToWord           to_word("abcABC");
  //     InversePresentation<word_type> p;
  //     p.alphabet(to_word("abcABC"));
  //     p.inverses_no_checks(to_word("ABCabc"));
  //     presentation::add_rule(p, to_word("ac"), to_word("ca"));
  //     presentation::add_rule(p, to_word("ab"), to_word("ba"));
  //     presentation::add_rule(p, to_word("bc"), to_word("cb"));
  //     auto S = Stephen(p);
  //     auto T = Stephen(p);
  //
  //     StringRange strings;
  //     strings.letters("abcABC").first("aaaaa").last("aaaaaaaaa");
  //
  //     for (auto const& w : strings) {
  //       S.set_word(to_word(w));
  //       auto first = ++stephen::cbegin_words_accepted(S);
  //       auto last  = first;
  //       std::advance(last, 1024);
  //       for (auto it = first; it != last; ++it) {
  //         REQUIRE(stephen::accepts(S, *it));
  //         if (!stephen::accepts(T.set_word(*it), S.word())) {
  //           std::cout << "S.word() = " << S.word() << ", T.word() = " <<
  //           T.word()
  //                     << std::endl;
  //         }
  //         REQUIRE(S == T);
  //       }
  //     }
  //   }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "046",
                          "non-inverse presentation -- operator==",
                          "[stephen][quick][no-valgrind]") {
    ReportGuard rg(false);
    auto        p = fpsemigroup::symmetric_inverse_monoid(4);

    ToddCoxeter tc(congruence_kind::twosided, p);

    auto      S = Stephen(p);
    word_type w = 0120112010222201_w;
    S.set_word(w);
    REQUIRE(stephen::number_of_words_accepted(S) == POSITIVE_INFINITY);

    {
      auto const index = todd_coxeter::index_of(tc, w);

      for (auto const& u : (stephen::words_accepted(S) | rx::take(1024))) {
        REQUIRE(todd_coxeter::index_of(tc, u) == index);
      }
    }
    {
      for (auto const& u : todd_coxeter::class_of(tc, w) | rx::take(1024)) {
        REQUIRE(stephen::accepts(S, u));
      }
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "032",
                          "Plactic monoid",
                          "[stephen][quick]") {
    auto p = fpsemigroup::plactic_monoid(4);
    p.contains_empty_word(true);
    Stephen s(p);
    s.set_word(0013122_w).run();
    REQUIRE(!stephen::accepts(s, 0013212_w));
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "033",
                          "Whyte's 4-relation full transf monoid 8",
                          "[stephen][fail]") {
    auto                    rg = ReportGuard(true);
    Presentation<word_type> p;
    p.rules = {00_w,       {},         11_w,         {},           22_w,
               {},         33_w,       {},           44_w,         {},
               55_w,       {},         66_w,         {},           010101_w,
               {},         121212_w,   {},           232323_w,     {},
               343434_w,   {},         454545_w,     {},           565656_w,
               {},         606060_w,   {},           01020102_w,   {},
               01030103_w, {},         01040104_w,   {},           01050105_w,
               {},         01060106_w, {},           12101210_w,   {},
               12131213_w, {},         12141214_w,   {},           12151215_w,
               {},         12161216_w, {},           23202320_w,   {},
               23212321_w, {},         23242324_w,   {},           23252325_w,
               {},         23262326_w, {},           34303430_w,   {},
               34313431_w, {},         34323432_w,   {},           34353435_w,
               {},         34363436_w, {},           45404540_w,   {},
               45414541_w, {},         45424542_w,   {},           45434543_w,
               {},         45464546_w, {},           56505650_w,   {},
               56515651_w, {},         56525652_w,   {},           56535653_w,
               {},         56545654_w, {},           60616061_w,   {},
               60626062_w, {},         60636063_w,   {},           60646064_w,
               {},         60656065_w, {},           071654321_w,  16543217_w,
               217121_w,   17171_w,    0102720107_w, 7010270102_w, 107017_w,
               70107010_w, 1217_w,     7121_w};
    p.alphabet_from_rules();
    presentation::balance_no_checks(p, 0123456_w, 0123456_w);

    Stephen s(p);
    s.set_word(1217_w);

    while (!s.finished()) {
      s.run_for(std::chrono::seconds(1));
    }

    REQUIRE(word_graph::last_node_on_path(s.word_graph(), 0, 1217_w).first
            == word_graph::last_node_on_path(s.word_graph(), 0, 7121_w).first);
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "045",
                          "Munn tree products",
                          "[stephen][quick]") {
    using words::pow;
    ToWord to_word("abcABC");

    InversePresentation<word_type> p;
    p.alphabet(to_word("abcABC"));
    p.inverses_no_checks(to_word("ABCabc"));

    auto S = Stephen(p);
    auto T = Stephen(p);

    S.set_word(to_word("aBbcaABAabCc"));
    T.set_word(to_word("aBbcaABAabCc"));
    S.run();
    REQUIRE(S.word_graph().number_of_nodes() == 7);

    T.run();
    REQUIRE(T.word_graph().number_of_nodes() == 7);
    S *= T;
    REQUIRE(S.word_graph().number_of_nodes() == 14);
    // fmt::print("{}", word_graph::dot(S.word_graph()));
    REQUIRE(!S.finished());
    S.run();
    REQUIRE(S.finished());
    REQUIRE(S.word_graph().number_of_nodes() == 12);
    REQUIRE(stephen::accepts(S, pow(T.word(), 2)));
    T.set_word(pow(T.word(), 2));
    T.run();
    REQUIRE(S == T);
    REQUIRE(stephen::words_accepted(S).get() == 0202_w);
    REQUIRE(stephen::number_of_left_factors(S) == POSITIVE_INFINITY);
    S.set_word(to_word("aBbcaABAabCc"));
    // std::ofstream f("wordgraph.gv");
    // f << stephen::dot(S).to_string();
    T.set_word(to_word("CcBAabaACBbA"));
    S *= T;
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "047",
                          "bicyclic monoid",
                          "[stephen][fail]") {
    ReportGuard                    rg(true);
    ToWord                         to_word("aA");
    InversePresentation<word_type> p;
    p.alphabet(to_word("aA"));
    p.inverses_no_checks(to_word("Aa"));
    p.contains_empty_word(true);
    presentation::add_rule(p, to_word("aA"), to_word(""));

    Stephen S(p);
    S.set_word(""_w).run();
    REQUIRE(S.word_graph().number_of_nodes() == 0);
  }

  LIBSEMIGROUPS_TEST_CASE("Stephen",
                          "048",
                          "plactic monoid",
                          "[stephen][fail]") {
    ReportGuard rg(true);
    auto        p = fpsemigroup::chinese_monoid(3);

    Stephen S(p);
    S.set_word(0110_w).run();
    REQUIRE(S.word_graph().number_of_nodes() == 8);
    // TODO rename words_accepted -> class_of
    REQUIRE((stephen::words_accepted(S) | rx::to_vector())
            == std::vector<word_type>());
  }

  // TODO the examples from Stephen's paper/thesis?
}  // namespace libsemigroups
