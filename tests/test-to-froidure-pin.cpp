//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2022 James D. Mitchell
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

#include <cstdint>  // for uint8_t
#include <vector>   // for vector

#include "catch_amalgamated.hpp"  // for REQUIRE, REQUIRE_THROWS_AS, REQUI...
#include "test-main.hpp"          // for LIBSEMIGROUPS_TEST_CASE

#include "libsemigroups/constants.hpp"        // for UNDEFINED
#include "libsemigroups/exception.hpp"        // for LIBSEMIGROUPS_EXCEPTION
#include "libsemigroups/froidure-pin.hpp"     // for FroidurePin, Froidure...
#include "libsemigroups/to-froidure-pin.hpp"  // for make
#include "libsemigroups/transf.hpp"           // for Transf
#include "libsemigroups/word-graph.hpp"       // for WordGraph

#include "libsemigroups/detail/report.hpp"  // for ReportGuard

namespace libsemigroups {

  LIBSEMIGROUPS_TEST_CASE("to_froidure_pin<Transf>",
                          "000",
                          "from WordGraph",
                          "[quick][to_froidure_pin]") {
    auto rg = ReportGuard(false);
    auto ad = make<WordGraph<uint8_t>>(
        5,
        {{1, 3, 4, 1}, {0, 0, 1, 1}, {2, 1, 2, 2}, {3, 2, 3, 3}, {4, 4, 4, 4}});
    auto T = to_froidure_pin<Transf<5>>(ad);
    REQUIRE(T.size() == 625);
    auto U = to_froidure_pin<Transf<0, uint8_t>>(ad);
    REQUIRE(U.size() == 625);
    auto V = to_froidure_pin<Transf<0, uint8_t>>(ad, 4, 5);
    REQUIRE(V.size() == 1);
    auto W = to_froidure_pin<Transf<0, uint8_t>>(ad, 0, 0);
    REQUIRE(W.size() == 1);
  }

  LIBSEMIGROUPS_TEST_CASE("to_froidure_pin<Transf>",
                          "001",
                          "from WordGraph (exceptions)",
                          "[quick][make]") {
    auto rg = ReportGuard(false);
    auto ad = make<WordGraph<uint8_t>>(
        5,
        {{1, 3, 4, 1}, {0, 0, 1, 1}, {2, 1, 2, 2}, {3, 2, 3, 3}, {4, 4, 4, 4}});
    REQUIRE_THROWS_AS((to_froidure_pin<Transf<0, uint8_t>>(ad, 10, 0)),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS((to_froidure_pin<Transf<0, uint8_t>>(ad, 10, 11)),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS((to_froidure_pin<Transf<0, uint8_t>>(ad, 0, 11)),
                      LibsemigroupsException);
  }

  namespace {
    template <typename Word, typename OtherWord = Word>
    void check_from_ke(Presentation<Word> const& p) {
      using literals::    operator""_w;
      Kambites<OtherWord> k(congruence_kind::twosided, p);
      auto                ptr = to_froidure_pin(k);
      auto& s = static_cast<FroidurePin<detail::KE<OtherWord>>&>(*ptr);
      REQUIRE(s.is_finite() == tril::FALSE);
      s.enumerate(100);
      REQUIRE(s.current_size() == 8'205);
      REQUIRE(s[0].to_word(k) == 0_w);
      REQUIRE(s[1].to_word(k) == 1_w);
      REQUIRE(s[2].to_word(k) == 2_w);
      REQUIRE(s[4].to_word(k) == 4_w);
      REQUIRE(s[8].to_word(k) == 01_w);
      REQUIRE(s[16].to_word(k) == 12_w);
      REQUIRE(s[32].to_word(k) == 34_w);
      REQUIRE(s[64].to_word(k) == 012_w);

      REQUIRE(froidure_pin::factorisation(s, s[0]) == 0_w);
      REQUIRE(froidure_pin::factorisation(s, s[1]) == 1_w);
      REQUIRE(froidure_pin::factorisation(s, s[2]) == 2_w);
      REQUIRE(froidure_pin::factorisation(s, s[4]) == 4_w);
      REQUIRE(froidure_pin::factorisation(s, s[8]) == 01_w);
      REQUIRE(froidure_pin::factorisation(s, s[16]) == 12_w);
      REQUIRE(froidure_pin::factorisation(s, s[32]) == 34_w);
      REQUIRE(froidure_pin::factorisation(s, s[64]) == 012_w);

      REQUIRE(s[0] < s[1]);

      auto t = s[128];
      t      = s[64];
      REQUIRE(froidure_pin::factorisation(s, t) == 012_w);
      auto u = s[128];
      u.swap(t);
      REQUIRE(froidure_pin::factorisation(s, t) == 134_w);
      REQUIRE(froidure_pin::factorisation(s, u) == 012_w);
      t = std::move(u);
      REQUIRE(froidure_pin::factorisation(s, t) == 012_w);

      auto v(std::move(t));
      REQUIRE(froidure_pin::factorisation(s, v) == 012_w);

      if constexpr (!std::is_same_v<Word, word_type>) {
        REQUIRE(s[0].to_string() == "a");
        REQUIRE(s[1].to_string() == "b");
        REQUIRE(s[2].to_string() == "c");
        REQUIRE(s[4].to_string() == "e");
        REQUIRE(s[8].to_string() == "ab");
        REQUIRE(s[16].to_string() == "bc");
        REQUIRE(s[32].to_string() == "de");
        REQUIRE(s[64].to_string() == "abc");
      }

      REQUIRE(Complexity<decltype(v)>()(v) == LIMIT_MAX);
      IncreaseDegree<decltype(v)>()(v);
      REQUIRE(Degree<decltype(v)>()(v) == 0);
      v = One<decltype(v)>()();
      REQUIRE(v.to_word(k) == ""_w);
    }
  }  // namespace

  LIBSEMIGROUPS_TEST_CASE("to_froidure_pin<KE>",
                          "002",
                          "from Kambites (code cov)",
                          "[quick][make][no-valgrind]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abcdefg");
    presentation::add_rule(p, "abcd", "aaaeaa");
    presentation::add_rule(p, "ef", "dg");
    check_from_ke(p);
    check_from_ke<std::string, detail::MultiStringView>(p);
    check_from_ke<word_type>(to_presentation<word_type>(p));
  }

  LIBSEMIGROUPS_TEST_CASE("to_froidure_pin<KE>",
                          "003",
                          "from Kambites (exceptions)",
                          "[quick][make]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("ab");
    presentation::add_rule(p, "bababa", "aba");
    REQUIRE_THROWS_AS(check_from_ke(p), LibsemigroupsException);
    REQUIRE_THROWS_AS((check_from_ke<std::string, detail::MultiStringView>(p)),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS((check_from_ke<word_type>(to_presentation<word_type>(p))),
                      LibsemigroupsException);
  }

  /*
  LIBSEMIGROUPS_TEST_CASE("to_froidure_pin<TCE>",
                          "002",
                          "from Todd-Coxeter",
                          "[quick][make]") {
    using TCE            = detail::TCE;
    using FroidurePinTCE = FroidurePin<TCE>;

    auto rg = ReportGuard(false);

    Presentation<word_type> p;
    p.alphabet(4);
    presentation::add_rule(p, {0, 0}, {0});
    presentation::add_rule(p, {1, 0}, {1});
    presentation::add_rule(p, {0, 1}, {1});
    presentation::add_rule(p, {2, 0}, {2});
    presentation::add_rule(p, {0, 2}, {2});
    presentation::add_rule(p, {3, 0}, {3});
    presentation::add_rule(p, {0, 3}, {3});
    presentation::add_rule(p, {1, 1}, {0});
    presentation::add_rule(p, {2, 3}, {0});
    presentation::add_rule(p, {2, 2, 2}, {0});
    presentation::add_rule(
        p, {1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2}, {0});
    presentation::add_rule(p,
                                     {1, 2, 1, 3, 1, 2, 1, 3, 1, 2, 1,
                                      3, 1, 2, 1, 3, 1, 2, 1, 3, 1, 2,
                                      1, 3, 1, 2, 1, 3, 1, 2, 1, 3},
                                     {0});
    // TODO
    ToddCoxeter tc(twosided, p);

    section_felsch(tc);

    // check_hlt(tc);
    // check_felsch(tc);
    // check_random(tc);
    // check_Rc_style(tc);
    // check_R_over_C_style(tc);
    // check_CR_style(tc);
    // check_Cr_style(tc);
    REQUIRE(tc.number_of_classes() == 10'752);
    check_complete_compatible(tc);
  }
      // Take a copy to test copy constructor
      auto& S = static_cast<FroidurePinTCE&>(*tc.quotient_froidure_pin());
      auto  T = S.copy_closure({S.generator(0)});

      REQUIRE(T.size() == S.size());
      REQUIRE(T.number_of_generators() == S.number_of_generators());

      REQUIRE(S.size() == 10'752);
      REQUIRE(S.number_of_idempotents() == 1);
      for (size_t c = 0; c < tc.number_of_classes(); ++c) {
        REQUIRE(tc.class_index_to_word(c) == froidure_pin::factorisation(S, c));
        REQUIRE(tc.word_to_class_index(tc.class_index_to_word(c)) == c);
      }
      REQUIRE(tc.finished());

      tc.standardize(order::recursive);
      REQUIRE(std::is_sorted(tc.cbegin_normal_forms(),
                             tc.cend_normal_forms(),
                             RecursivePathCompare<word_type>{}));
      REQUIRE(std::vector<word_type>(tc.cbegin_normal_forms(),
                                     tc.cbegin_normal_forms() + 10)
              == std::vector<word_type>({{{0},
                                          {1},
                                          {2},
                                          {2, 1},
                                          {1, 2},
                                          {1, 2, 1},
                                          {2, 2},
                                          {2, 2, 1},
                                          {2, 1, 2},
                                          {2, 1, 2, 1}}}));

      tc.standardize(order::lex);
      for (size_t c = 0; c < tc.number_of_classes(); ++c) {
        REQUIRE(tc.word_to_class_index(tc.class_index_to_word(c)) == c);
      }
      REQUIRE(std::is_sorted(tc.cbegin_normal_forms(),
                             tc.cend_normal_forms(),
                             LexicographicalCompare<word_type>{}));
      REQUIRE(std::vector<word_type>(tc.cbegin_normal_forms(),
                                     tc.cbegin_normal_forms() + 10)
              == std::vector<word_type>({{0},
                                         {0, 1},
                                         {0, 1, 2},
                                         {0, 1, 2, 1},
                                         {0, 1, 2, 1, 2},
                                         {0, 1, 2, 1, 2, 1},
                                         {0, 1, 2, 1, 2, 1, 2},
                                         {0, 1, 2, 1, 2, 1, 2, 1},
                                         {0, 1, 2, 1, 2, 1, 2, 1, 2},
                                         {0, 1, 2, 1, 2, 1, 2, 1, 2, 1}}));
      tc.standardize(order::shortlex);
      for (size_t c = 0; c < tc.number_of_classes(); ++c) {
        REQUIRE(tc.word_to_class_index(tc.class_index_to_word(c)) == c);
      }
      REQUIRE(std::is_sorted(tc.cbegin_normal_forms(),
                             tc.cend_normal_forms(),
                             ShortLexCompare<word_type>{}));
      REQUIRE(std::vector<word_type>(tc.cbegin_normal_forms(),
                                     tc.cbegin_normal_forms() + 10)
              == std::vector<word_type>({{0},
                                         {1},
                                         {2},
                                         {3},
                                         {1, 2},
                                         {1, 3},
                                         {2, 1},
                                         {3, 1},
                                         {1, 2, 1},
                                         {1, 3, 1}}));
    }
*/

}  // namespace libsemigroups
