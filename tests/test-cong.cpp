//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2019 James D. Mitchell
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

#include "catch_amalgamated.hpp"  // for TEST_CASE
#include "test-main.hpp"          // for LIBSEMIGROUPS_TEST_CASE_V3

#include "libsemigroups/bmat-fastest.hpp"     // for BMatFastest
#include "libsemigroups/bmat8.hpp"            // for BMat8
#include "libsemigroups/cong.hpp"             // for Congruence
#include "libsemigroups/fpsemi-examples.hpp"  // for rook_monoid
#include "libsemigroups/froidure-pin.hpp"     // for FroidurePin
#include "libsemigroups/obvinf.hpp"           // for is_obviously_infinite
#include "libsemigroups/pbr.hpp"              // for PBR
#include "libsemigroups/to-froidure-pin.hpp"
#include "libsemigroups/transf.hpp"      // for Transf<>
#include "libsemigroups/types.hpp"       // for word_type
#include "libsemigroups/word-range.hpp"  // for literals

#include "libsemigroups/detail/report.hpp"  // for ReportGuard

namespace libsemigroups {

  using congruence::contains;
  using congruence::non_trivial_classes;
  using congruence::reduce;

  namespace {
    auto sizes_ntc(std::vector<std::vector<word_type>> const& ntc) {
      std::vector<size_t> sizes(ntc.size(), 0);
      std::transform(ntc.cbegin(),
                     ntc.cend(),
                     sizes.begin(),
                     std::mem_fn(&std::vector<word_type>::size));
      std::sort(sizes.begin(), sizes.end());
      return sizes;
    }

    auto unique_sizes_ntc(std::vector<std::vector<word_type>> const& ntc) {
      auto sizes = sizes_ntc(ntc);
      sizes.erase(std::unique(sizes.begin(), sizes.end()), sizes.end());
      return sizes;
    }
  }  // namespace

  // Forward declarations
  struct LibsemigroupsException;

  using literals::operator""_w;

  constexpr congruence_kind twosided = congruence_kind::twosided;
  constexpr congruence_kind onesided = congruence_kind::onesided;

  using fpsemigroup::stellar_monoid;
  using fpsemigroup::zero_rook_monoid;

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "000",
                             "left congruence from presentation",
                             "[quick][cong]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(2);
    presentation::add_rule(p, 000_w, 0_w);
    presentation::add_rule(p, 0_w, 11_w);
    presentation::reverse(p);

    Congruence cong(onesided, p);

    REQUIRE(cong.number_of_classes() == 5);
    REQUIRE(congruence::contains(cong, 100110_w, 100_w));

    REQUIRE(congruence::contains(cong, 100_w, 10000_w));
    REQUIRE(congruence::contains(cong, 100110_w, 100_w));
    REQUIRE(!congruence::contains(cong, 000_w, 100_w));
    REQUIRE(!congruence::contains(cong, 1_w, 0000_w));
    REQUIRE(!congruence::contains(cong, 0000_w, 100_w));

    if (!cong.has<KnuthBendix<>>()) {
      REQUIRE(cong.presentation().rules == p.rules);
    }
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "001",
                             "2-sided congruence from presentation",
                             "[quick][cong]") {
    auto rg = ReportGuard(false);

    Presentation<word_type> p;
    p.alphabet(2);
    presentation::add_rule(p, 000_w, 0_w);
    presentation::add_rule(p, 0_w, 11_w);
    // REQUIRE(!is_obviously_infinite(p));

    Congruence cong(twosided, p);

    REQUIRE(cong.number_of_classes() == 5);

    REQUIRE(congruence::contains(cong, 001_w, 00001_w));
    REQUIRE(congruence::contains(cong, 001_w, 001_w));
    REQUIRE(congruence::contains(cong, 001_w, 00001_w));
    REQUIRE(congruence::contains(cong, 00001_w, 011001_w));
    REQUIRE(!congruence::contains(cong, 000_w, 001_w));
    REQUIRE(!congruence::contains(cong, 000_w, 1_w));
    REQUIRE(!congruence::contains(cong, 1_w, 000_w));
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "007",
                             "2-sided congruence from presentation",
                             "[quick][cong]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(3);
    presentation::add_rule(p, 01_w, 10_w);
    presentation::add_rule(p, 02_w, 22_w);
    presentation::add_rule(p, 02_w, 0_w);
    presentation::add_rule(p, 02_w, 0_w);
    presentation::add_rule(p, 22_w, 0_w);
    presentation::add_rule(p, 12_w, 12_w);
    presentation::add_rule(p, 12_w, 22_w);
    presentation::add_rule(p, 122_w, 1_w);
    presentation::add_rule(p, 12_w, 1_w);
    presentation::add_rule(p, 22_w, 1_w);
    presentation::add_rule(p, 0_w, 1_w);
    presentation::add_rule(p, 0_w, 1_w);

    Congruence cong(twosided, p);

    REQUIRE(cong.number_of_classes() == 2);
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "008",
                             "2-sided congruence from presentation",
                             "[quick][cong]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(3);
    presentation::add_rule(p, 01_w, 10_w);
    presentation::add_rule(p, 02_w, 22_w);
    presentation::add_rule(p, 02_w, 0_w);
    presentation::add_rule(p, 02_w, 0_w);
    presentation::add_rule(p, 22_w, 0_w);
    presentation::add_rule(p, 12_w, 12_w);
    presentation::add_rule(p, 12_w, 22_w);
    presentation::add_rule(p, 122_w, 1_w);
    presentation::add_rule(p, 12_w, 1_w);
    presentation::add_rule(p, 22_w, 1_w);
    presentation::add_rule(p, 0_w, 1_w);

    Congruence cong(twosided, p);

    REQUIRE(congruence::contains(cong, 0_w, 1_w));
    REQUIRE(congruence::contains(cong, 0_w, 10_w));
    REQUIRE(congruence::contains(cong, 0_w, 11_w));
    REQUIRE(congruence::contains(cong, 0_w, 101_w));

    REQUIRE(congruence::contains(cong, 1_w, 11_w));
    REQUIRE(congruence::contains(cong, 101_w, 10_w));
    REQUIRE(cong.number_of_classes() == 2);
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "009",
                             "infinite 2-sided congruence from presentation",
                             "[quick][cong]") {
    auto rg = ReportGuard(false);

    Presentation<word_type> p;
    p.alphabet(3);
    presentation::add_rule(p, 01_w, 10_w);
    presentation::add_rule(p, 02_w, 20_w);
    presentation::add_rule(p, 00_w, 0_w);
    presentation::add_rule(p, 02_w, 0_w);
    presentation::add_rule(p, 20_w, 0_w);
    presentation::add_rule(p, 12_w, 21_w);
    presentation::add_rule(p, 111_w, 1_w);
    presentation::add_rule(p, 12_w, 1_w);
    presentation::add_rule(p, 21_w, 1_w);

    Congruence cong(twosided, p);
    congruence::add_generating_pair(cong, 0_w, 1_w);

    REQUIRE(cong.number_of_classes() == POSITIVE_INFINITY);

    // Used to require KnuthBendix<>CongruenceByPairs to work
    REQUIRE(congruence::contains(cong, 0_w, 1_w));
    REQUIRE(congruence::contains(cong, 0_w, 10_w));
    REQUIRE(congruence::contains(cong, 0_w, 11_w));
    REQUIRE(congruence::contains(cong, 0_w, 101_w));
    REQUIRE(congruence::contains(cong, 1_w, 11_w));
    REQUIRE(congruence::contains(cong, 101_w, 10_w));

    REQUIRE(cong.has<KnuthBendix<>>());

    KnuthBendix<> kb(twosided, p);
    REQUIRE(knuth_bendix::non_trivial_classes(kb, *cong.get<KnuthBendix<>>())
            == std::vector<std::vector<std::string>>(
                {{{1}, {0, 1}, {1, 1}, {0, 1, 1}, {0}}}));

    REQUIRE(non_trivial_classes(cong, p)
            == std::vector<std::vector<word_type>>(
                {{1_w, 01_w, 11_w, 011_w, 0_w}}));

    REQUIRE(non_trivial_classes(
                cong, to_presentation<std::string>(p, [](auto x) { return x; }))
            == std::vector<std::vector<std::string>>(
                {{{1}, {0, 1}, {1, 1}, {0, 1, 1}, {0}}}));
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "010",
                             "2-sided congruence on finite semigroup",
                             "[quick][cong][no-valgrind]") {
    auto rg      = ReportGuard(false);
    using Transf = LeastTransf<8>;
    auto S       = to_froidure_pin({Transf({7, 3, 5, 3, 4, 2, 7, 7}),
                                    Transf({1, 2, 4, 4, 7, 3, 0, 7}),
                                    Transf({0, 6, 4, 2, 2, 6, 6, 4}),
                                    Transf({3, 6, 3, 4, 0, 6, 0, 7})});

    REQUIRE(S.size() == 11804);
    REQUIRE(S.number_of_rules() == 2460);

    Congruence cong(twosided, S, S.right_cayley_graph());
    congruence::add_generating_pair(cong, 0321322_w, 322133_w);

    REQUIRE(cong.number_of_classes() == 525);

    REQUIRE(congruence::contains(cong, 0001_w, 00100_w));
    REQUIRE(congruence::contains(cong, 00101_w, 1101_w));
    REQUIRE(!congruence::contains(cong, 1100_w, 0001_w));
    REQUIRE(!congruence::contains(cong, 003_w, 0001_w));
    REQUIRE(!congruence::contains(cong, 1100_w, 003_w));
    REQUIRE(congruence::contains(cong, 12133212_w, 2133210_w));
    REQUIRE(congruence::contains(cong, 0311132210_w, 03221_w));
    REQUIRE(!congruence::contains(cong, 0321333_w, 003_w));
    REQUIRE(!congruence::contains(cong, 110_w, 1332210_w));

    REQUIRE(congruence::contains(cong, 12133212_w, 2133210_w));
    REQUIRE(!congruence::contains(cong, 110_w, 1332210_w));

    REQUIRE(cong.number_of_classes() == 525);
  }

  // TODO(1) this does not seem to be functioning all that well, one of the
  // threads blocks the others from stopping, extending the time
  // taken for this to run.
  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "011",
                             "congruence on full PBR monoid on 2 points",
                             "[extreme][cong]") {
    auto rg = ReportGuard(true);
    auto S  = to_froidure_pin({PBR({{2}, {3}, {0}, {1}}),
                               PBR({{}, {2}, {1}, {0, 3}}),
                               PBR({{0, 3}, {2}, {1}, {}}),
                               PBR({{1, 2}, {3}, {0}, {1}}),
                               PBR({{2}, {3}, {0}, {1, 3}}),
                               PBR({{3}, {1}, {0}, {1}}),
                               PBR({{3}, {2}, {0}, {0, 1}}),
                               PBR({{3}, {2}, {0}, {1}}),
                               PBR({{3}, {2}, {0}, {3}}),
                               PBR({{3}, {2}, {1}, {0}}),
                               PBR({{3}, {2, 3}, {0}, {1}})});

    // REQUIRE(S.size() == 65536);
    // REQUIRE(S.number_of_rules() == 45416);

    Congruence cong(twosided, S, S.right_cayley_graph());
    congruence::add_generating_pair(
        cong, {7, 10, 9, 3, 6, 9, 4, 7, 9, 10}, {9, 3, 6, 6, 10, 9, 4, 7});
    congruence::add_generating_pair(
        cong, {8, 7, 5, 8, 9, 8}, {6, 3, 8, 6, 1, 2, 4});

    REQUIRE(cong.number_of_classes() == 19'009);
    auto ntc
        = congruence::non_trivial_classes(cong, froidure_pin::normal_forms(S));
    REQUIRE(ntc.size() == 577);

    std::vector<size_t> v(577, 0);
    std::transform(ntc.cbegin(),
                   ntc.cend(),
                   v.begin(),
                   std::mem_fn(&std::vector<word_type>::size));
    REQUIRE(std::count(v.cbegin(), v.cend(), 4) == 384);
    REQUIRE(std::count(v.cbegin(), v.cend(), 16) == 176);
    REQUIRE(std::count(v.cbegin(), v.cend(), 96) == 16);
    REQUIRE(std::count(v.cbegin(), v.cend(), 41216) == 1);
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "012",
                             "2-sided congruence on finite semigroup",
                             "[quick][cong][no-valgrind]") {
    auto rg = ReportGuard(false);

    auto S
        = to_froidure_pin({LeastPPerm<6>({0, 1, 2}, {4, 0, 1}, 6),
                           LeastPPerm<6>({0, 1, 2, 3, 5}, {2, 5, 3, 0, 4}, 6),
                           LeastPPerm<6>({0, 1, 2, 3}, {5, 0, 3, 1}, 6),
                           LeastPPerm<6>({0, 2, 5}, {3, 4, 1}, 6),
                           LeastPPerm<6>({0, 2, 5}, {0, 2, 5}, 6),
                           LeastPPerm<6>({0, 1, 4}, {1, 2, 0}, 6),
                           LeastPPerm<6>({0, 2, 3, 4, 5}, {3, 0, 2, 5, 1}, 6),
                           LeastPPerm<6>({0, 1, 3, 5}, {1, 3, 2, 0}, 6),
                           LeastPPerm<6>({1, 3, 4}, {5, 0, 2}, 6)});

    // REQUIRE(S.size() == 712);
    // REQUIRE(S.number_of_rules() == 1121);

    Congruence cong(twosided, S, S.right_cayley_graph());
    congruence::add_generating_pair(cong, {2, 7}, {1, 6, 6, 1});
    REQUIRE(cong.number_of_classes() == 32);
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "013",
                             "trivial 2-sided congruence on bicyclic monoid",
                             "[quick][cong][no-valgrind]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.contains_empty_word(true);
    p.alphabet(2);
    presentation::add_rule(p, 01_w, {});

    Congruence cong(twosided, p);
    REQUIRE(congruence::contains(cong, {}, 010011_w));
    REQUIRE(congruence::contains(cong, {}, 0101_w));
    REQUIRE(congruence::contains(cong, 10_w, 011001_w));
    REQUIRE(congruence::contains(cong, 10_w, 011001_w));
  }

  LIBSEMIGROUPS_TEST_CASE_V3(
      "Congruence",
      "014",
      "non-trivial 2-sided congruence on bicyclic monoid",
      "[quick][cong]") {
    auto rg = ReportGuard(false);

    Presentation<word_type> p;
    p.contains_empty_word(true).alphabet(2);
    presentation::add_rule(p, 01_w, {});

    REQUIRE(is_obviously_infinite(p));

    Congruence cong(twosided, p);
    congruence::add_generating_pair(cong, 111_w, {});
    REQUIRE(cong.number_of_classes() == 3);

    if (cong.has<KnuthBendix<>>()) {
      KnuthBendix<> kb(twosided, p);
      REQUIRE_THROWS_AS(
          knuth_bendix::non_trivial_classes(kb, *cong.get<KnuthBendix<>>()),
          LibsemigroupsException);
    } else if (cong.has<ToddCoxeter>()) {
      ToddCoxeter tc(twosided, p);
      REQUIRE_THROWS_AS(
          todd_coxeter::non_trivial_classes(*cong.get<ToddCoxeter>(), tc),
          LibsemigroupsException);
    }
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "015",
                             "2-sided congruence on free abelian monoid",
                             "[quick][cong]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(3);
    presentation::add_rule(p, 12_w, 21_w);
    presentation::add_identity_rules(p, 0);

    Congruence cong(twosided, p);
    congruence::add_generating_pair(cong, 11111_w, 1_w);
    congruence::add_generating_pair(cong, 222_w, 2_w);

    REQUIRE(cong.number_of_classes() == 15);
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "016",
                             "example where TC works but KB doesn't",
                             "[quick][cong]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abBe");
    presentation::add_identity_rules(p, 'e');
    presentation::add_rule(p, "aa", "e");
    presentation::add_rule(p, "BB", "b");
    presentation::add_rule(p, "BaBaBaB", "abababa");
    presentation::add_rule(p, "aBabaBabaBabaBab", "BabaBabaBabaBaba");

    Congruence cong(twosided, p);
    congruence::add_generating_pair(cong, "a", "b");

    REQUIRE(cong.number_of_classes() == 4);
    REQUIRE(!to_froidure_pin(cong)->contains_one());
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "017",
                             "2-sided congruence on finite semigroup",
                             "[quick][cong]") {
    auto rg      = ReportGuard(false);
    using Transf = LeastTransf<5>;
    auto S
        = to_froidure_pin({Transf({1, 3, 4, 2, 3}), Transf({3, 2, 1, 3, 3})});

    REQUIRE(S.size() == 88);
    REQUIRE(S.number_of_rules() == 18);

    word_type w1 = froidure_pin::factorisation(S, Transf({3, 4, 4, 4, 4}));
    word_type w2 = froidure_pin::factorisation(S, Transf({3, 4, 4, 4, 4}));

    Congruence cong(twosided, S, S.right_cayley_graph());
    congruence::add_generating_pair(
        cong,
        froidure_pin::factorisation(S, Transf({3, 4, 4, 4, 4})),
        froidure_pin::factorisation(S, Transf({3, 1, 3, 3, 3})));
    REQUIRE(cong.number_of_classes() == 21);

    word_type u = froidure_pin::factorisation(S, Transf({1, 3, 1, 3, 3}));
    word_type v = froidure_pin::factorisation(S, Transf({4, 2, 4, 4, 2}));
    REQUIRE(congruence::contains(cong, u, v));
    REQUIRE(congruence::contains(cong, u, v));
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "018",
                             "infinite fp semigroup from GAP library",
                             "[quick][cong]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(3);
    presentation::add_rule(p, 00_w, 00_w);
    presentation::add_rule(p, 01_w, 10_w);
    presentation::add_rule(p, 02_w, 20_w);
    presentation::add_rule(p, 00_w, 0_w);
    presentation::add_rule(p, 02_w, 0_w);
    presentation::add_rule(p, 20_w, 0_w);
    presentation::add_rule(p, 10_w, 01_w);
    presentation::add_rule(p, 11_w, 11_w);
    presentation::add_rule(p, 12_w, 21_w);
    presentation::add_rule(p, 111_w, 1_w);
    presentation::add_rule(p, 12_w, 1_w);
    presentation::add_rule(p, 21_w, 1_w);

    // REQUIRE(is_obviously_infinite(p));

    Congruence cong(twosided, p);
    congruence::add_generating_pair(cong, 0_w, 1_w);

    REQUIRE(is_obviously_infinite(cong));
    REQUIRE(cong.number_of_classes() == POSITIVE_INFINITY);

    REQUIRE(cong.has<KnuthBendix<>>());

    KnuthBendix<> kb(twosided, p);
    auto          ntc
        = knuth_bendix::non_trivial_classes(kb, *cong.get<KnuthBendix<>>());
    REQUIRE(ntc.size() == 1);
    REQUIRE(ntc[0].size() == 5);
    REQUIRE(ntc[0]
            == std::vector<std::string>({{1}, {0, 1}, {1, 1}, {0, 1, 1}, {0}}));
  }

  LIBSEMIGROUPS_TEST_CASE_V3(
      "Congruence",
      "019",
      "2-sided cong. from presentation with infinite classes ",
      "[quick][cong]") {
    using words::operator+;
    using words::pow;

    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(2);
    presentation::add_rule(p, 01_w, 10_w);
    presentation::add_rule(p, 000_w, 00_w);

    Congruence cong(twosided, p);
    congruence::add_generating_pair(cong, {0}, {1});

    word_type x = "0"_w + pow(1_w, 20);
    word_type y = "00"_w + pow(1_w, 20);

    REQUIRE(congruence::contains(cong, x, y));
    REQUIRE(congruence::contains(cong, y, x));
    REQUIRE(congruence::contains(cong, x, x));
    REQUIRE(congruence::contains(cong, y, y));
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "020",
                             "trivial cong. on an fp semigroup",
                             "[quick][cong][no-valgrind]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("ab");
    presentation::add_rule(p, "ab", "ba");
    presentation::add_rule(p, "a", "b");
    presentation::reverse(p);

    Congruence cong(onesided, p);

    // No generating pairs for the congruence (not the fp semigroup) means no
    // non-trivial classes.
    // REQUIRE(cong.number_of_non_trivial_classes() == 0);
    REQUIRE(cong.number_of_classes() == POSITIVE_INFINITY);
    REQUIRE(cong.has<KnuthBendix<>>());
    REQUIRE(knuth_bendix::non_trivial_classes(*cong.get<KnuthBendix<>>(),
                                              *cong.get<KnuthBendix<>>())
                .empty());
    REQUIRE(cong.finished());
    REQUIRE(cong.started());
    REQUIRE_THROWS_AS(congruence::add_generating_pair(cong, 00_w, 0_w),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "021",
                             "duplicate generators",
                             "[quick][cong]") {
    auto rg      = ReportGuard(false);
    using Transf = LeastTransf<8>;
    auto       S = to_froidure_pin({Transf({7, 3, 5, 3, 4, 2, 7, 7}),
                                    Transf({7, 3, 5, 3, 4, 2, 7, 7}),
                                    Transf({7, 3, 5, 3, 4, 2, 7, 7}),
                                    Transf({3, 6, 3, 4, 0, 6, 0, 7})});
    Congruence cong(twosided, S, S.right_cayley_graph());
    REQUIRE(cong.number_of_classes() == S.size());
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "022",
                             "non-trivial classes",
                             "[quick][cong]") {
    auto rg = ReportGuard(false);

    Presentation<word_type> p;
    p.alphabet(2);
    presentation::add_rule(p, 000_w, 0_w);
    presentation::add_rule(p, 100_w, 10_w);
    presentation::add_rule(p, 10111_w, 10_w);
    presentation::add_rule(p, 11111_w, 11_w);
    presentation::add_rule(p, 110110_w, 101011_w);
    presentation::add_rule(p, 0010110_w, 010110_w);
    presentation::add_rule(p, 0011010_w, 011010_w);
    presentation::add_rule(p, 0101010_w, 101010_w);
    presentation::add_rule(p, 1010101_w, 101010_w);
    presentation::add_rule(p, 1010110_w, 101011_w);
    presentation::add_rule(p, 1011010_w, 101101_w);
    presentation::add_rule(p, 1101010_w, 101010_w);
    presentation::add_rule(p, 1111010_w, 1010_w);
    presentation::add_rule(p, 00111010_w, 111010_w);

    Congruence cong(twosided, p);
    congruence::add_generating_pair(cong, 0_w, 1_w);
    REQUIRE(cong.number_of_classes() == 1);
    if (cong.has<ToddCoxeter>() && cong.get<ToddCoxeter>()->finished()) {
      ToddCoxeter tc(twosided, p);
      REQUIRE(tc.number_of_classes() == 78);
      auto ntc
          = todd_coxeter::non_trivial_classes(*cong.get<ToddCoxeter>(), tc);
      REQUIRE(ntc.size() == 1);
      REQUIRE(ntc[0].size() == 78);
    }
    if (cong.has<KnuthBendix<>>() && cong.get<KnuthBendix<>>()->finished()) {
      KnuthBendix<> kb(twosided, p);
      REQUIRE(kb.number_of_classes() == 78);
      auto ntc
          = knuth_bendix::non_trivial_classes(kb, *cong.get<KnuthBendix<>>());
      REQUIRE(ntc.size() == 1);
      REQUIRE(ntc[0].size() == 78);
    }
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "023",
                             "onesided congruence on finite semigroup",
                             "[quick][cong][no-valgrind]") {
    auto rg      = ReportGuard(false);
    using Transf = LeastTransf<8>;
    auto S       = to_froidure_pin({Transf({0, 1, 2, 3, 4, 5, 6, 7}),
                                    Transf({1, 2, 3, 4, 5, 0, 6, 7}),
                                    Transf({1, 0, 2, 3, 4, 5, 6, 7}),
                                    Transf({0, 1, 2, 3, 4, 0, 6, 7}),
                                    Transf({0, 1, 2, 3, 4, 5, 7, 6})});

    REQUIRE(S.size() == 93'312);
    std::vector<Transf> elms = {Transf({0, 0, 0, 0, 0, 0, 7, 6}),
                                Transf({0, 0, 0, 0, 0, 0, 6, 7}),
                                Transf({0, 0, 0, 0, 0, 0, 6, 7}),
                                Transf({1, 1, 1, 1, 1, 1, 6, 7}),
                                Transf({0, 0, 0, 0, 0, 0, 6, 7}),
                                Transf({2, 2, 2, 2, 2, 2, 6, 7}),
                                Transf({0, 0, 0, 0, 0, 0, 6, 7}),
                                Transf({3, 3, 3, 3, 3, 3, 6, 7}),
                                Transf({0, 0, 0, 0, 0, 0, 6, 7}),
                                Transf({4, 4, 4, 4, 4, 4, 6, 7}),
                                Transf({0, 0, 0, 0, 0, 0, 6, 7}),
                                Transf({5, 5, 5, 5, 5, 5, 6, 7}),
                                Transf({0, 0, 0, 0, 0, 0, 7, 6}),
                                Transf({0, 1, 2, 3, 4, 5, 7, 6})};
    REQUIRE(
        std::all_of(elms.cbegin(), elms.cend(), [&S](Transf const& x) -> bool {
          return S.contains(x);
        }));

    Congruence cong(onesided, S, S.right_cayley_graph());
    word_type  w1, w2;
    for (size_t i = 0; i < elms.size(); i += 2) {
      froidure_pin::factorisation(S, w1, S.position(elms[i]));
      froidure_pin::factorisation(S, w2, S.position(elms[i + 1]));
      congruence::add_generating_pair(cong, w1, w2);
    }
    REQUIRE(cong.number_of_classes() == 1);
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "024",
                             "redundant generating pairs",
                             "[quick][cong]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(1);

    Congruence cong(twosided, p);
    congruence::add_generating_pair(cong, 00_w, 00_w);
    REQUIRE(congruence::contains(cong, 00_w, 00_w));
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "025",
                             "2-sided cong. on free semigroup",
                             "[quick][cong]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("a");
    Congruence cong(twosided, p);
    cong.run();
    REQUIRE(congruence::contains(cong, "aa", "aa"));
    REQUIRE(!congruence::contains(cong, "aa", "a"));
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "026",
                             "is_obviously_(in)finite",
                             "[quick][cong]") {
    auto rg = ReportGuard(false);
    {
      Presentation<word_type> p;
      p.alphabet(3);
      presentation::add_rule(p, {0, 1}, {0});
      Congruence cong(twosided, p);
      congruence::add_generating_pair(cong, {2, 2}, {2});
      REQUIRE(is_obviously_infinite(cong));
      REQUIRE(cong.number_of_classes() == POSITIVE_INFINITY);
    }
    {
      Presentation<word_type> p;
      p.alphabet(3);
      presentation::add_rule(p, {0, 1}, {0});
      presentation::add_rule(p, {0, 0}, {0});
      Congruence cong(twosided, p);
      congruence::add_generating_pair(cong, {1, 1}, {1});
      REQUIRE(is_obviously_infinite(cong));
    }
    {
      Presentation<word_type> p;
      p.alphabet(3);
      presentation::add_rule(p, {0, 1}, {0});
      presentation::add_rule(p, {0, 0}, {0});
      Congruence cong(twosided, p);
      congruence::add_generating_pair(cong, {1, 2}, {1});
      REQUIRE(is_obviously_infinite(cong));
    }
    {
      Presentation<word_type> p;
      p.alphabet(3);
      presentation::add_rule(p, {0, 1}, {0});
      Congruence cong(onesided, p);
      congruence::add_generating_pair(cong, {2, 2}, {2});
      REQUIRE(is_obviously_infinite(cong));
    }
    {
      Presentation<word_type> p;
      p.alphabet(3);
      presentation::add_rule(p, {0, 1}, {0});
      presentation::add_rule(p, {0, 0}, {0});
      Congruence cong(onesided, p);
      congruence::add_generating_pair(cong, {1, 1}, {1});
      REQUIRE(is_obviously_infinite(cong));
    }
    {
      Presentation<word_type> p;
      p.alphabet(3);
      presentation::add_rule(p, {0, 1}, {0});
      presentation::add_rule(p, {0, 0}, {0});
      Congruence cong(onesided, p);
      congruence::add_generating_pair(cong, {1, 2}, {1});
      REQUIRE(is_obviously_infinite(cong));
    }
    {
      Presentation<word_type> p;
      p.alphabet(3);
      presentation::add_rule(p, {0, 1}, {0});
      presentation::reverse(p);
      Congruence cong(onesided, p);
      congruence::add_generating_pair(cong, {2, 2}, {2});
      REQUIRE(is_obviously_infinite(cong));
    }
    {
      Presentation<word_type> p;
      p.alphabet(3);
      presentation::add_rule(p, {0, 1}, {0});
      presentation::add_rule(p, {0, 0}, {0});
      presentation::reverse(p);
      Congruence cong(onesided, p);
      congruence::add_generating_pair(cong, {1, 1}, {1});
      REQUIRE(is_obviously_infinite(cong));
    }
    {
      Presentation<word_type> p;
      p.alphabet(3);
      presentation::add_rule(p, {0, 1}, {0});
      presentation::add_rule(p, {0, 0}, {0});
      presentation::reverse(p);
      Congruence cong(onesided, p);
      congruence::add_generating_pair(cong, {1, 2}, {1});
      REQUIRE(is_obviously_infinite(cong));
    }

    using Transf = LeastTransf<3>;
    auto fp      = to_froidure_pin({Transf({0, 1, 0}), Transf({0, 1, 2})});
    REQUIRE(fp.size() == 2);
    {
      Congruence cong(twosided, fp, fp.right_cayley_graph());
      congruence::add_generating_pair(cong, {1}, {0});
      REQUIRE(!is_obviously_infinite(cong));

      REQUIRE(cong.number_of_classes() == 1);
    }
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "028",
                             "2-sided congruences of BMat8 semigroup",
                             "[quick][cong][no-valgrind]") {
#ifdef LIBSEMIGROUPS_HPCOMBI_ENABLED
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic ignored "-Winline"
#endif
    auto rg    = ReportGuard(false);
    using BMat = BMatFastest<4>;
    std::vector<BMat> gens
        = {BMat({{0, 1, 0, 0}, {1, 0, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}}),
           BMat({{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}),
           BMat({{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {1, 0, 0, 1}}),
           BMat({{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 0}})};
    {
      auto S = to_froidure_pin(gens);

      Congruence cong(twosided, S, S.right_cayley_graph());
      congruence::add_generating_pair(cong, {1}, {0});

      REQUIRE(cong.number_of_classes() == 3);
      REQUIRE(congruence::contains(cong, {1}, {0}));

      auto ntc = non_trivial_classes(cong, froidure_pin::normal_forms(S));
      REQUIRE(ntc.size() == 3);
      REQUIRE(ntc[0].size() == 12);
      REQUIRE(ntc[1].size() == 63'880);
      REQUIRE(ntc[2].size() == 12);
      REQUIRE(ntc[0]
              == std::vector<word_type>({0_w,
                                         1_w,
                                         010_w,
                                         011_w,
                                         101_w,
                                         110_w,
                                         111_w,
                                         01011_w,
                                         01101_w,
                                         10110_w,
                                         10111_w,
                                         11011_w}));
    }
    {
      auto       S = to_froidure_pin({gens[0], gens[2], gens[3]});
      Congruence cong(twosided, S, S.right_cayley_graph());
      congruence::add_generating_pair(cong, {1}, {0});

      REQUIRE(cong.number_of_classes() == 2);
      REQUIRE(congruence::contains(cong, {1}, {0}));

      auto ntc = non_trivial_classes(cong, froidure_pin::normal_forms(S));
      REQUIRE(ntc.size() == 2);
      REQUIRE(ntc[0].size() == 8);
      REQUIRE(ntc[1].size() == 8);

      REQUIRE(ntc[0]
              == std::vector<word_type>(
                  {0_w, 1_w, 00_w, 01_w, 10_w, 010_w, 101_w, 0101_w}));
    }
#ifdef LIBSEMIGROUPS_HPCOMBI_ENABLED
#pragma GCC diagnostic pop
#endif
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "029",
                             "left congruence on finite semigroup",
                             "[quick][cong]") {
    auto                  rg = ReportGuard(false);
    FroidurePin<Transf<>> S;
    S.add_generator(Transf<>({1, 3, 4, 2, 3}));
    S.add_generator(Transf<>({3, 2, 1, 3, 3}));

    REQUIRE(S.size() == 88);
    REQUIRE(S.degree() == 5);

    Congruence cong(onesided, S, S.left_cayley_graph());

    congruence::add_generating_pair(cong, 001100010_w, 10001_w);

    REQUIRE(cong.number_of_classes() == 69);
    REQUIRE(cong.number_of_classes() == 69);

    word_type w3
        = reverse(froidure_pin::factorisation(S, Transf<>({1, 3, 1, 3, 3})));
    word_type w4
        = reverse(froidure_pin::factorisation(S, Transf<>({4, 2, 4, 4, 2})));
    REQUIRE(!congruence::contains(cong, w3, w4));
    REQUIRE(congruence::contains(cong, w3, 10100_w));
    REQUIRE(congruence::contains(cong, 101001_w, 1000100_w));
    REQUIRE(!congruence::contains(cong, 000110_w, 11_w));
    REQUIRE(!congruence::contains(cong, 00010001_w, 1001_w));

    REQUIRE(congruence::contains(cong, 101001_w, 1000100_w));
    REQUIRE(!congruence::contains(cong, 00010001_w, 1001_w));
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "030",
                             "onesided congruence on finite semigroup",
                             "[quick][cong]") {
    auto                  rg = ReportGuard(false);
    FroidurePin<Transf<>> S;
    S.add_generator(Transf<>({1, 3, 4, 2, 3}));
    S.add_generator(Transf<>({3, 2, 1, 3, 3}));

    // REQUIRE(S.size() == 88);
    // REQUIRE(S.degree() == 5);
    Congruence cong(onesided, S, S.right_cayley_graph());
    congruence::add_generating_pair(cong, 010001100_w, 10001_w);

    REQUIRE(cong.number_of_classes() == 72);
    REQUIRE(cong.number_of_classes() == 72);

    word_type w3 = froidure_pin::factorisation(S, Transf<>({1, 3, 1, 3, 3}));
    word_type w4 = froidure_pin::factorisation(S, Transf<>({4, 2, 4, 4, 2}));
    REQUIRE(!congruence::contains(cong, w3, w4));
    REQUIRE(!congruence::contains(cong, w3, 00101_w));
    REQUIRE(!congruence::contains(cong, 100101_w, 0010001_w));
    REQUIRE(!congruence::contains(cong, 011000_w, 11_w));
    REQUIRE(!congruence::contains(cong, 10001000_w, 1001_w));

    REQUIRE(!congruence::contains(cong, 100101_w, 0010001_w));
    REQUIRE(!congruence::contains(cong, 10001000_w, 1001_w));
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "031",
                             "onesided congruence on finite semigroup",
                             "[quick][cong]") {
    auto                  rg = ReportGuard(false);
    FroidurePin<Transf<>> S;
    S.add_generator(Transf<>({1, 3, 4, 2, 3}));
    S.add_generator(Transf<>({3, 2, 1, 3, 3}));

    REQUIRE(S.size() == 88);
    REQUIRE(S.number_of_rules() == 18);
    REQUIRE(S.degree() == 5);
    word_type w1, w2;
    froidure_pin::factorisation(S, w1, S.position(Transf<>({3, 4, 4, 4, 4})));
    froidure_pin::factorisation(S, w2, S.position(Transf<>({3, 1, 3, 3, 3})));
    Congruence cong(onesided, S, S.right_cayley_graph());
    congruence::add_generating_pair(cong, w1, w2);

    REQUIRE(cong.number_of_classes() == 72);
    REQUIRE(cong.number_of_classes() == 72);
    word_type w3, w4, w5, w6;
    froidure_pin::factorisation(S, w3, S.position(Transf<>({1, 3, 3, 3, 3})));
    froidure_pin::factorisation(S, w4, S.position(Transf<>({4, 2, 4, 4, 2})));
    froidure_pin::factorisation(S, w5, S.position(Transf<>({2, 3, 2, 2, 2})));
    froidure_pin::factorisation(S, w6, S.position(Transf<>({2, 3, 3, 3, 3})));
    REQUIRE(!congruence::contains(cong, w3, w4));
    REQUIRE(congruence::contains(cong, w5, w6));
    REQUIRE(!congruence::contains(cong, w3, w6));

    REQUIRE(congruence::contains(cong, w1, w2));
    REQUIRE(congruence::contains(cong, w5, w6));
    REQUIRE(!congruence::contains(cong, w3, w5));
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence", "032", "contains", "[quick][cong]") {
    using namespace std::string_literals;
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("ab");
    REQUIRE(!p.contains_empty_word());
    Congruence cong(twosided, p);
    congruence::add_generating_pair(cong, "aa", "a");
    congruence::add_generating_pair(cong, "ab", "a");
    congruence::add_generating_pair(cong, "ba", "a");

    REQUIRE(cong.has<KnuthBendix<>>());
    auto& kb = *cong.get<KnuthBendix<>>();
    REQUIRE((kb.active_rules() | rx::to_vector() | rx::count()) == 0);

    cong.run();

    REQUIRE(cong.has<KnuthBendix<>>());
    REQUIRE(kb.finished());
    REQUIRE(kb.confluent());
    REQUIRE(kb.presentation().alphabet() == "ab");
    REQUIRE(kb.presentation().rules
            == std::vector<std::string>({"aa", "a", "ab", "a", "ba", "a"}));
    REQUIRE(!kb.presentation().contains_empty_word());
    REQUIRE((kb.active_rules() | rx::to_vector())
            == std::vector<std::decay_t<decltype(kb)>::rule_type>(
                {{"ba", "a"}, {"ab", "a"}, {"aa", "a"}}));

    REQUIRE(knuth_bendix::reduce(kb, "bb") == "bb");
    REQUIRE(knuth_bendix::reduce(kb, "b") == "b");
    REQUIRE(!knuth_bendix::contains(kb, "bb", "b"));
    REQUIRE(congruence::currently_contains_no_checks(cong, "bb", "b")
            == tril::FALSE);
    REQUIRE(congruence::contains(cong, "aa", "a"));
    REQUIRE(congruence::contains(cong, "ab", "a"));
    REQUIRE(congruence::contains(cong, "ba", "a"));
    REQUIRE(congruence::contains(cong, "ba", "abab"));
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "033",
                             "stellar_monoid S2",
                             "[quick][cong]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p  = zero_rook_monoid(2);
    presentation::change_alphabet(p, 10_w);

    REQUIRE(!is_obviously_infinite(p));

    Congruence cong(twosided, p);
    auto       q = stellar_monoid(2);
    presentation::change_alphabet(q, 10_w);
    for (auto it = q.rules.cbegin(); it != q.rules.cend(); it += 2) {
      congruence::add_generating_pair(cong, *it, *(it + 1));
    }

    REQUIRE(!is_obviously_infinite(cong));
    REQUIRE(cong.number_of_classes() == 5);

    ToddCoxeter tc(twosided, p);
    REQUIRE(tc.number_of_classes() == 7);

    auto ntc = non_trivial_classes(cong, todd_coxeter::normal_forms(tc));
    REQUIRE(ntc.size() == 1);
    std::sort(ntc[0].begin(), ntc[0].end());
    REQUIRE(ntc[0] == std::vector<word_type>({010_w, 10_w, 101_w}));
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "034",
                             "stellar_monoid S3",
                             "[quick][cong]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p  = zero_rook_monoid(3);

    REQUIRE(!is_obviously_infinite(p));

    Congruence cong(twosided, p);
    auto       q = stellar_monoid(3);
    for (auto it = q.rules.cbegin(); it != q.rules.cend(); it += 2) {
      congruence::add_generating_pair(cong, *it, *(it + 1));
    }
    REQUIRE(!is_obviously_infinite(cong));
    REQUIRE(cong.number_of_classes() == 16);

    ToddCoxeter tc(twosided, p);
    REQUIRE(tc.number_of_classes() == 34);

    auto ntc = non_trivial_classes(cong, todd_coxeter::normal_forms(tc));
    REQUIRE(ntc.size() == 4);
    std::for_each(ntc.begin(), ntc.end(), [](auto& val) {
      std::sort(val.begin(), val.end());
    });
    std::sort(ntc.begin(), ntc.end());
    REQUIRE(ntc
            == decltype(ntc)({{0102_w,
                               01020_w,
                               010201_w,
                               012012_w,
                               0120120_w,
                               01202_w,
                               102_w,
                               1020_w,
                               10201_w,
                               12012_w,
                               120120_w,
                               1202_w,
                               201202_w},
                              {012_w, 0201_w, 2012_w},
                              {0120_w, 01201_w, 20120_w},
                              {02_w, 020_w, 202_w}}));
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "035",
                             "stellar_monoid S4",
                             "[quick][cong][no-valgrind]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p  = zero_rook_monoid(4);

    REQUIRE(!is_obviously_infinite(p));

    Congruence cong(twosided, p);
    auto       q = stellar_monoid(4);
    for (auto it = q.rules.cbegin(); it != q.rules.cend(); it += 2) {
      congruence::add_generating_pair(cong, *it, *(it + 1));
    }

    REQUIRE(!is_obviously_infinite(cong));
    REQUIRE(cong.number_of_classes() == 65);

    ToddCoxeter tc(twosided, p);
    REQUIRE(tc.number_of_classes() == 209);

    auto ntc = non_trivial_classes(cong, todd_coxeter::normal_forms(tc));
    REQUIRE(ntc.size() == 17);

    auto sizes = sizes_ntc(ntc);
    REQUIRE(sizes
            == std::vector<size_t>(
                {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 13, 13, 13, 13, 73}));
    REQUIRE(std::accumulate(sizes.cbegin(), sizes.cend(), 0)
                + (cong.number_of_classes() - ntc.size())
            == 209);
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "036",
                             "stellar_monoid S5",
                             "[quick][cong][no-valgrind]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p  = zero_rook_monoid(5);

    Congruence cong(twosided, p);
    auto       q = stellar_monoid(5);
    for (auto it = q.rules.cbegin(); it != q.rules.cend(); it += 2) {
      congruence::add_generating_pair(cong, *it, *(it + 1));
    }

    REQUIRE(!is_obviously_infinite(cong));
    REQUIRE(cong.number_of_classes() == 326);

    ToddCoxeter tc(twosided, p);
    REQUIRE(tc.number_of_classes() == 1'546);

    auto ntc = non_trivial_classes(cong, todd_coxeter::normal_forms(tc));
    REQUIRE(ntc.size() == 86);

    auto sizes = sizes_ntc(ntc);
    // REQUIRE(cong.number_of_non_trivial_classes() == 86);

    REQUIRE(std::count(sizes.cbegin(), sizes.cend(), 3) == 60);
    REQUIRE(std::count(sizes.cbegin(), sizes.cend(), 13) == 20);
    REQUIRE(std::count(sizes.cbegin(), sizes.cend(), 73) == 5);
    REQUIRE(std::count(sizes.cbegin(), sizes.cend(), 501) == 1);
    REQUIRE(std::accumulate(sizes.cbegin(), sizes.cend(), 0)
                + (cong.number_of_classes() - ntc.size())
            == tc.number_of_classes());
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "037",
                             "stellar_monoid S6",
                             "[quick][cong][no-valgrind]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p  = zero_rook_monoid(6);

    Congruence cong(twosided, p);
    auto       q = stellar_monoid(6);
    for (auto it = q.rules.cbegin(); it != q.rules.cend(); it += 2) {
      congruence::add_generating_pair(cong, *it, *(it + 1));
    }

    REQUIRE(!is_obviously_infinite(cong));
    REQUIRE(cong.number_of_classes() == 1'957);

    ToddCoxeter tc(twosided, p);
    REQUIRE(tc.number_of_classes() == 13'327);

    auto ntc = non_trivial_classes(cong, todd_coxeter::normal_forms(tc));
    REQUIRE(ntc.size() == 517);

    REQUIRE(unique_sizes_ntc(ntc)
            == std::vector<size_t>({3, 13, 73, 501, 4051}));
    auto sizes = sizes_ntc(ntc);
    REQUIRE(std::accumulate(sizes.cbegin(), sizes.cend(), 0)
                + (cong.number_of_classes() - ntc.size())
            == tc.number_of_classes());
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "038",
                             "stellar_monoid S7",
                             "[quick][cong][no-valgrind]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p  = zero_rook_monoid(7);

    Congruence cong(twosided, p);
    auto       q = stellar_monoid(7);
    for (auto it = q.rules.cbegin(); it != q.rules.cend(); it += 2) {
      congruence::add_generating_pair(cong, *it, *(it + 1));
    }

    REQUIRE(!is_obviously_infinite(cong));
    REQUIRE(cong.number_of_classes() == 13'700);

    ToddCoxeter tc(twosided, p);
    REQUIRE(tc.number_of_classes() == 130'922);

    auto ntc = non_trivial_classes(cong, todd_coxeter::normal_forms(tc));
    REQUIRE(ntc.size() == 3'620);

    REQUIRE(unique_sizes_ntc(ntc)
            == std::vector<size_t>({3, 13, 73, 501, 4'051, 37'633}));
    auto sizes = sizes_ntc(ntc);
    REQUIRE(std::accumulate(sizes.cbegin(), sizes.cend(), 0)
                + (cong.number_of_classes() - ntc.size())
            == tc.number_of_classes());
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "039",
                             "left cong. on an f.p. semigroup",
                             "[quick][cong]") {
    auto rg = ReportGuard(false);

    Presentation<std::string> p;
    p.alphabet("abe");
    presentation::add_identity_rules(p, 'e');
    presentation::add_rule(p, "abb", "bb");
    presentation::add_rule(p, "bbb", "bb");
    presentation::add_rule(p, "aaaa", "a");
    presentation::add_rule(p, "baab", "bb");
    presentation::add_rule(p, "baaab", "b");
    presentation::add_rule(p, "babab", "b");
    presentation::add_rule(p, "bbaaa", "bb");
    presentation::add_rule(p, "bbaba", "bbaa");

    presentation::reverse(p);

    Congruence cong1(onesided, p);
    congruence::add_generating_pair(cong1, "a", "bbb");
    REQUIRE(cong1.number_of_classes() == 11);

    Congruence cong2(onesided, p);
    congruence::add_generating_pair(cong2, "bb", "aaaaaaa");
    REQUIRE(cong1.number_of_classes() == cong2.number_of_classes());
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "040",
                             "2-sided cong. on infinite f.p. semigroup",
                             "[quick][cong]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(3);

    Congruence cong(twosided, p);
    congruence::add_generating_pair(cong, 1_w, 2_w);
    congruence::add_generating_pair(cong, 00_w, 0_w);
    congruence::add_generating_pair(cong, 01_w, 10_w);
    congruence::add_generating_pair(cong, 01_w, 1_w);
    congruence::add_generating_pair(cong, 02_w, 20_w);
    congruence::add_generating_pair(cong, 02_w, 2_w);

    REQUIRE(!congruence::contains(cong, 1_w, 2222222222_w));
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "042",
                             "const_contains",
                             "[quick][cong][no-valgrind]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(2);
    Congruence cong(twosided, p);
    congruence::add_generating_pair(cong, 000_w, 0_w);
    congruence::add_generating_pair(cong, 1111_w, 1_w);
    congruence::add_generating_pair(cong, 01110_w, 00_w);
    congruence::add_generating_pair(cong, 1001_w, 11_w);
    congruence::add_generating_pair(cong, 001010101010_w, 00_w);

    REQUIRE(congruence::currently_contains(cong, 1111_w, 11_w)
            == tril::unknown);
    REQUIRE(congruence::contains(cong, 1111_w, 1_w));
    REQUIRE(!congruence::contains(cong, 1111_w, 11_w));
    if (cong.has<ToddCoxeter>()) {
      REQUIRE_NOTHROW(cong.get<ToddCoxeter>());
    }
    if (cong.has<KnuthBendix<>>()) {
      REQUIRE_NOTHROW(cong.get<KnuthBendix<>>());
    }
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "043",
                             "no winner",
                             "[quick][cong]") {
    auto                    rg = ReportGuard(false);
    Presentation<word_type> p;
    p.alphabet(2);

    for (auto const& knd : {onesided, twosided}) {
      Congruence cong(knd, p);
      // Required in case of using a 1 core computer, otherwise the tests
      // below fail.
      cong.max_threads(2);
      REQUIRE(congruence::contains(cong, 0000_w, 0000_w));
      REQUIRE(!congruence::contains(cong, 0000_w, 0001_w));
      if (knd == twosided) {
        REQUIRE_NOTHROW(to_froidure_pin(cong));
      } else {
        REQUIRE_THROWS_AS(to_froidure_pin(cong), LibsemigroupsException);
      }

      WordRange w;
      w.alphabet_size(2).min(1).max(5);

      REQUIRE(w.count() == 30);

      REQUIRE(non_trivial_classes(cong, w).empty());
    }
  }

  LIBSEMIGROUPS_TEST_CASE_V3("Congruence",
                             "044",
                             "congruence over smalloverlap",
                             "[quick][cong]") {
    auto                      rg = ReportGuard(false);
    Presentation<std::string> p;
    p.alphabet("abcdefg");
    presentation::add_rule(p, "abcd", "aaaeaa");
    Congruence cong(twosided, p);
    congruence::add_generating_pair(cong, "ef", "dg");
    REQUIRE(is_obviously_infinite(cong));
    REQUIRE(cong.number_of_generating_pairs() == 1);
    cong.run();
    REQUIRE(cong.number_of_classes() == POSITIVE_INFINITY);
    REQUIRE(!cong.has<ToddCoxeter>());
    REQUIRE(cong.has<Kambites<word_type>>());
    REQUIRE(cong.get<Kambites<word_type>>()->finished());
    REQUIRE(cong.get<Kambites<word_type>>()->success());

    StringRange w;
    w.alphabet("abcdefg").min(1).max(4);
    REQUIRE(w.count() == 399);
    // REQUIRE(cong.get<Kambites<word_type>>()->presentation().alphabet()
    //         == word_type({0, 1, 2, 3, 4, 5, 6}));
    REQUIRE(non_trivial_classes(cong, w)
            == std::vector<std::vector<std::string>>({{"dg", "ef"},
                                                      {"adg", "aef"},
                                                      {"bdg", "bef"},
                                                      {"cdg", "cef"},
                                                      {"ddg", "def"},
                                                      {"dga", "efa"},
                                                      {"dgb", "efb"},
                                                      {"dgc", "efc"},
                                                      {"dgd", "efd"},
                                                      {"dge", "efe"},
                                                      {"dgf", "eff"},
                                                      {"dgg", "efg"},
                                                      {"edg", "eef"},
                                                      {"fdg", "fef"},
                                                      {"gdg", "gef"}}));
  }

}  // namespace libsemigroups
