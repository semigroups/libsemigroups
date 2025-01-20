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

#include <cstddef>  // for size_t
#include <vector>   // for vector

#include "catch_amalgamated.hpp"  // for REQUIRE
#include "test-main.hpp"          // for LIBSEMIGROUPS_TEST_CASE

#include "libsemigroups/bipart.hpp"        // for Bipartition
#include "libsemigroups/froidure-pin.hpp"  // for FroidurePin
#include "libsemigroups/word-range.hpp"    // for namespace literals

namespace libsemigroups {
  using namespace literals;  // for operator""_w

  // Forward declaration
  struct LibsemigroupsException;

  constexpr bool REPORT = false;

  LIBSEMIGROUPS_TEST_CASE("FroidurePin",
                          "000",
                          "small example 1 .................. Bipartition",
                          "[quick][froidure-pin][bipartition][bipart]") {
    auto        rg = ReportGuard(REPORT);
    FroidurePin S  = to_froidure_pin(
        {make<Bipartition>({{1, 5, 8, -1, -2, -4, -10},
                             {2, 4, 7, -8},
                             {3, 6, 9, 10, -3},
                             {-5, -9},
                             {-6, -7}}),
          make<Bipartition>({{1},
                             {2, 3, 4, 5, -5, -6, -7},
                             {6, 8, -2, -4, -8, -10},
                             {7, -9},
                             {9, -3},
                             {10, -1}}),
          make<Bipartition>({{1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
                              -1, -2, -3, -4, -5, -6, -7, -8, -9, -10}})});

    S.reserve(10);

    REQUIRE(S.size() == 10);
    REQUIRE(S.number_of_idempotents() == 6);

    size_t pos = 0;
    for (auto const& x : froidure_pin::elements(S)) {
      REQUIRE(S.position(x) == pos);
      pos++;
    }

    auto x = make<Bipartition>({{1, -2, -10},
                                {2, 4, 5, 7, -1, -8},
                                {3, 9, -4},
                                {6, 10, -3, -5, -9},
                                {8, -7},
                                {-6}});
    S.add_generator(x);
    REQUIRE(S.number_of_generators() == 4);
    REQUIRE(S.size() == 21);
    froidure_pin::closure(S, {x});
    REQUIRE(S.number_of_generators() == 4);
    REQUIRE(S.size() == 21);
    REQUIRE(froidure_pin::minimal_factorisation(
                S, S.generator(0) * S.generator(3) * S.generator(1))
            == 031_w);
    REQUIRE(froidure_pin::minimal_factorisation(S, 11) == 03_w);
    REQUIRE(S.at(11) == S.generator(0) * S.generator(3));
    REQUIRE_THROWS_AS(froidure_pin::minimal_factorisation(S, 1'000'000'000),
                      LibsemigroupsException);
    pos = 0;
    std::for_each(
        S.cbegin_idempotents(), S.cend_idempotents(), [&pos](auto const& x) {
          REQUIRE(x * x == x);
          ++pos;
        });
    REQUIRE(pos == S.number_of_idempotents());
    REQUIRE(std::is_sorted(S.cbegin_sorted(), S.cend_sorted()));
  }

  LIBSEMIGROUPS_TEST_CASE("FroidurePin",
                          "001",
                          "default constructed .............. Bipartition",
                          "[quick][froidure-pin][bipartition][bipart]") {
    auto rg = ReportGuard(REPORT);

    FroidurePin<Bipartition> S;
    S.add_generator(Bipartition(
        {0, 1, 2, 1, 0, 2, 1, 0, 2, 2, 0, 0, 2, 0, 3, 4, 4, 1, 3, 0}));
    S.add_generator(Bipartition(
        {0, 1, 1, 1, 1, 2, 3, 2, 4, 5, 5, 2, 4, 2, 1, 1, 1, 2, 3, 2}));
    S.add_generator(Bipartition(
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}));

    S.reserve(10);

    REQUIRE(S.size() == 10);
    REQUIRE(S.number_of_idempotents() == 6);

    size_t pos = 0;
    for (auto const& x : froidure_pin::elements(S)) {
      REQUIRE(S.position(x) == pos);
      pos++;
    }

    S.add_generator(Bipartition(
        {0, 1, 2, 1, 1, 3, 1, 4, 2, 3, 1, 0, 3, 2, 3, 5, 4, 1, 3, 0}));

    REQUIRE(S.size() == 21);
    froidure_pin::closure(S, {Bipartition({0, 1, 2, 1, 1, 3, 1, 4, 2, 3,
                                           1, 0, 3, 2, 3, 5, 4, 1, 3, 0})});
    REQUIRE(S.size() == 21);
    REQUIRE(froidure_pin::minimal_factorisation(
                S,
                Bipartition({0, 1, 2, 1, 0, 2, 1, 0, 2, 2,
                             0, 0, 2, 0, 3, 4, 4, 1, 3, 0})
                    * Bipartition({0, 1, 2, 1, 1, 3, 1, 4, 2, 3,
                                   1, 0, 3, 2, 3, 5, 4, 1, 3, 0})
                    * Bipartition({0, 1, 1, 1, 1, 2, 3, 2, 4, 5,
                                   5, 2, 4, 2, 1, 1, 1, 2, 3, 2}))
            == 031_w);
    REQUIRE(froidure_pin::minimal_factorisation(S, 11) == 03_w);
    REQUIRE(S.at(11)
            == Bipartition(
                   {0, 1, 2, 1, 0, 2, 1, 0, 2, 2, 0, 0, 2, 0, 3, 4, 4, 1, 3, 0})
                   * Bipartition({0, 1, 2, 1, 1, 3, 1, 4, 2, 3,
                                  1, 0, 3, 2, 3, 5, 4, 1, 3, 0}));
    REQUIRE_THROWS_AS(froidure_pin::minimal_factorisation(S, 1'000'000'000),
                      LibsemigroupsException);
    pos = 0;
    for (auto const& x : froidure_pin::idempotents(S)) {
      REQUIRE(x * x == x);
      pos++;
    }
    REQUIRE(pos == S.number_of_idempotents());
    for (auto it = S.cbegin_sorted() + 1; it < S.cend_sorted(); ++it) {
      REQUIRE(*(it - 1) < *it);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("FroidurePin",
                          "002",
                          "small example 2 - Bipartition",
                          "[quick][froidure-pin][element]") {
    auto                     rg = ReportGuard(REPORT);
    FroidurePin<Bipartition> S;
    S.add_generator(Bipartition(
        {0, 1, 2, 1, 0, 2, 1, 0, 2, 2, 0, 0, 2, 0, 3, 4, 4, 1, 3, 0}));
    S.add_generator(Bipartition(
        {0, 1, 1, 1, 1, 2, 3, 2, 4, 5, 5, 2, 4, 2, 1, 1, 1, 2, 3, 2}));
    S.add_generator(Bipartition(
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}));

    REQUIRE(S.size() == 10);
    REQUIRE(S.degree() == 10);
    REQUIRE(S.number_of_idempotents() == 6);
    REQUIRE(S.number_of_generators() == 3);
    REQUIRE(S.number_of_rules() == 14);

    REQUIRE(S[0] == S.generator(0));
    REQUIRE(S[1] == S.generator(1));
    REQUIRE(S[2] == S.generator(2));

    REQUIRE(S.position(S.generator(0)) == 0);
    REQUIRE(S.contains(S.generator(0)));

    REQUIRE(S.position(S.generator(1)) == 1);
    REQUIRE(S.contains(S.generator(1)));

    REQUIRE(S.position(S.generator(2)) == 2);
    REQUIRE(S.contains(S.generator(2)));

    Bipartition y({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
    REQUIRE(S.position(y) == 2);
    REQUIRE(S.contains(y));
    y.product_inplace_no_checks(S.generator(0), S.generator(1));
    REQUIRE(S.position(y) == 4);
    REQUIRE(S.contains(y));
    y.product_inplace_no_checks(S.generator(1), S.generator(2));
    REQUIRE(S.position(y) == 7);
    REQUIRE(S.contains(y));
  }

  LIBSEMIGROUPS_TEST_CASE("FroidurePin",
                          "003",
                          "number of idempotents - Bipartition",
                          "[extreme][froidure-pin][element]") {
    auto                     rg = ReportGuard();
    FroidurePin<Bipartition> S;
    S.add_generator(Bipartition({0, 1, 2, 3, 4, 5, 5, 0, 1, 2, 3, 4}));
    S.add_generator(Bipartition({0, 1, 2, 3, 4, 5, 1, 0, 2, 3, 4, 5}));
    S.add_generator(Bipartition({0, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5}));
    S.add_generator(Bipartition({0, 0, 1, 2, 3, 4, 0, 0, 1, 2, 3, 4}));
    S.reserve(4'213'597);
    REQUIRE(S.size() == 4'213'597);
    REQUIRE(S.number_of_idempotents() == 541'254);
  }

  LIBSEMIGROUPS_TEST_CASE("FroidurePin",
                          "004",
                          "exception: is_idempotent - Bipartition",
                          "[quick][froidure-pin][element]") {
    FroidurePin<Bipartition> S;
    S.add_generator(Bipartition(
        {0, 1, 2, 1, 0, 2, 1, 0, 2, 2, 0, 0, 2, 0, 3, 4, 4, 1, 3, 0}));
    S.add_generator(Bipartition(
        {0, 1, 1, 1, 1, 2, 3, 2, 4, 5, 5, 2, 4, 2, 1, 1, 1, 2, 3, 2}));
    S.add_generator(Bipartition(
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}));

    for (size_t i = 0; i < S.size(); ++i) {
      REQUIRE_NOTHROW(S.is_idempotent(i));
    }
    REQUIRE_THROWS_AS(S.is_idempotent(10), LibsemigroupsException);
  }

}  // namespace libsemigroups
