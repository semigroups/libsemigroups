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
#include <cstdint>  // for int64_t

#include "catch_amalgamated.hpp"  // for REQUIRE, REQUIRE_THROWS_AS
#include "test-main.hpp"          // for LIBSEMIGROUPS_TEST_CASE

#include "libsemigroups/constants.hpp"     // for UNDEFINED
#include "libsemigroups/froidure-pin.hpp"  // for FroidurePin
#include "libsemigroups/matrix.hpp"        // for IntMat
#include "libsemigroups/word-range.hpp"    // for namespace literals

namespace libsemigroups {
  using namespace literals;  // for operator""_w

  // Forward declaration
  struct LibsemigroupsException;

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
      "FroidurePin",
      "032",
      "Example 000",
      "[quick][froidure-pin][intmat][no-sanitize-undefined]",
      (IntMat<0, 0, int64_t>),
      (IntMat<2, 2, int64_t>) ) {
    // this test seemingly causes undefined behaviour (multiplication of
    // signed integers that overflows), which is either a bug or a feature
    // depending on your perspective.
    auto rg = ReportGuard(false);

    FroidurePin<TestType> S;
    S.add_generator(TestType({{0, 1}, {0, -1}}));
    S.add_generator(TestType({{0, 1}, {2, 0}}));
    REQUIRE(to_human_readable_repr(S)
            == "<partially enumerated FroidurePin with 2 generators, 2 "
               "elements, Cayley graph ⌀ 1, & 0 rules>");

    REQUIRE(TestType({{0, 1}, {0, -1}}) * TestType({{0, 1}, {2, 0}})
                * TestType({{0, 1}, {2, 0}})
            == S.generator(0) * S.generator(1) * S.generator(0));
    REQUIRE(make<TestType>({{64, 0}, {-64, 0}})
            == S.generator(0) * S.generator(1) * S.generator(0) * S.generator(1)
                   * S.generator(0) * S.generator(1) * S.generator(0)
                   * S.generator(1) * S.generator(0) * S.generator(1)
                   * S.generator(0) * S.generator(1));

    S.reserve(10'000);

    S.enumerate(10'000);
    REQUIRE(S.finished());
    REQUIRE(S.current_size() == 631);
    size_t pos = 0;

    for (auto it = S.cbegin(); it < S.cend(); ++it) {
      REQUIRE(S.position(*it) == pos);
      pos++;
    }
    S.enumerate(1'000'000);
    REQUIRE(S.current_size() == 631);
    REQUIRE(to_human_readable_repr(S)
            == "<fully enumerated FroidurePin with 2 generators, 631 "
               "elements, Cayley graph ⌀ 128, & 7 rules>");
    REQUIRE(froidure_pin::minimal_factorisation(
                S,
                TestType({{0, 1}, {0, -1}}) * TestType({{0, 1}, {2, 0}})
                    * TestType({{0, 1}, {2, 0}}))
            == 010_w);
    REQUIRE(froidure_pin::minimal_factorisation(S, 52) == 010101010101_w);
    REQUIRE(S.at(52) == TestType({{64, 0}, {-64, 0}}));
    REQUIRE_THROWS_AS(froidure_pin::minimal_factorisation(S, 1'000'000'000),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("FroidurePin",
                                   "034",
                                   "Example 001",
                                   "[quick][froidure-pin][intmat]",
                                   (IntMat<0, 0, int64_t>),
                                   (IntMat<2, 2, int64_t>) ) {
    auto                  rg = ReportGuard(false);
    FroidurePin<TestType> S;
    S.add_generator(TestType({{0, 0}, {0, 1}}));
    S.add_generator(TestType({{0, 1}, {-1, 0}}));

    REQUIRE(S.size() == 13);
    REQUIRE(S.degree() == 2);
    REQUIRE(S.number_of_idempotents() == 4);
    REQUIRE(S.number_of_generators() == 2);
    REQUIRE(S.number_of_rules() == 6);
    REQUIRE(S[0] == S.generator(0));
    REQUIRE(S[1] == S.generator(1));

    REQUIRE(S.position(S.generator(0)) == 0);
    REQUIRE(S.contains(S.generator(0)));

    REQUIRE(S.position(S.generator(1)) == 1);
    REQUIRE(S.contains(S.generator(1)));

    TestType x({{-2, 2}, {-1, 0}});
    REQUIRE(S.position(x) == UNDEFINED);
    REQUIRE(!S.contains(x));

    x.product_inplace_no_checks(S.generator(1), S.generator(1));
    REQUIRE(S.position(x) == 4);
    REQUIRE(S.contains(x));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("FroidurePin",
                                   "036",
                                   "exception: current_position",
                                   "[quick][froidure-pin][element]",
                                   IntMat<2>,
                                   IntMat<>) {
    FroidurePin<TestType> T;
    T.add_generator(TestType({{0, 0}, {0, 1}}));
    T.add_generator(TestType({{0, 1}, {-1, 0}}));
    REQUIRE(froidure_pin::current_position(T, {}) == UNDEFINED);
    REQUIRE_NOTHROW(froidure_pin::current_position(T, 0011_w));
    REQUIRE(froidure_pin::current_position(T, 0011_w) == UNDEFINED);
    auto w = froidure_pin::to_element(T, 0011_w);
    REQUIRE(T.current_position(w) == UNDEFINED);
    REQUIRE_THROWS_AS(froidure_pin::current_position(T, 0012_w),
                      LibsemigroupsException);

    REQUIRE(T.size() == 13);
    REQUIRE(froidure_pin::current_position(T, 0011_w) == 6);
    w = froidure_pin::to_element(T, 0011_w);
    REQUIRE(T.current_position(w) == 6);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("FroidurePin",
                                   "037",
                                   "exception: to_element",
                                   "[quick][froidure-pin][element]",
                                   IntMat<2>,
                                   IntMat<>) {
    FroidurePin<TestType> T;
    T.add_generator(TestType({{0, 0}, {0, 1}}));
    T.add_generator(TestType({{0, 1}, {-1, 0}}));

    REQUIRE_THROWS_AS(froidure_pin::to_element(T, {}), LibsemigroupsException);
    REQUIRE_THROWS_AS(froidure_pin::to_element(T, {0, 0, 1, 2}),
                      LibsemigroupsException);

    auto t = froidure_pin::to_element(T, {0, 0, 1, 1});
    REQUIRE(
        t == T.generator(0) * T.generator(0) * T.generator(1) * T.generator(1));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE(
      "FroidurePin",
      "038",
      "exception: prefix, suffix, first_letter",
      "[quick][froidure-pin][element][no-valgrind]",
      IntMat<2>,
      IntMat<>) {
    FroidurePin<TestType> T;
    T.add_generator(TestType({{0, 0}, {0, 1}}));
    T.add_generator(TestType({{0, 1}, {-1, 0}}));

    for (size_t i = 0; i < T.size(); ++i) {
      REQUIRE_NOTHROW(T.prefix(i));
      REQUIRE_THROWS_AS(T.prefix(i + T.size()), LibsemigroupsException);
    }
    for (size_t i = 0; i < T.size(); ++i) {
      REQUIRE_NOTHROW(T.suffix(i));
      REQUIRE_THROWS_AS(T.suffix(i + T.size()), LibsemigroupsException);
    }
    for (size_t i = 0; i < T.size(); ++i) {
      REQUIRE_NOTHROW(T.first_letter(i));
      REQUIRE_THROWS_AS(T.first_letter(i + T.size()), LibsemigroupsException);
    }
    for (size_t i = 0; i < T.size(); ++i) {
      REQUIRE_NOTHROW(T.final_letter(i));
      REQUIRE_THROWS_AS(T.final_letter(i + T.size()), LibsemigroupsException);
    }
    for (size_t i = 0; i < T.size(); ++i) {
      REQUIRE_NOTHROW(T.current_length(i));
      REQUIRE_THROWS_AS(T.current_length(i + T.size()), LibsemigroupsException);
    }
    for (size_t i = 0; i < T.size(); ++i) {
      for (size_t j = 0; j < T.size(); ++j) {
        REQUIRE_NOTHROW(froidure_pin::product_by_reduction(T, i, j));
        REQUIRE_THROWS_AS(
            froidure_pin::product_by_reduction(T, i + T.size(), j),
            LibsemigroupsException);
        REQUIRE_THROWS_AS(
            froidure_pin::product_by_reduction(T, i, j + T.size()),
            LibsemigroupsException);
        REQUIRE_THROWS_AS(
            froidure_pin::product_by_reduction(T, i + T.size(), j + T.size()),
            LibsemigroupsException);
      }
    }
    for (size_t i = 0; i < T.size(); ++i) {
      for (size_t j = 0; j < T.size(); ++j) {
        REQUIRE_NOTHROW(T.fast_product(i, j));
        REQUIRE_THROWS_AS(T.fast_product(i + T.size(), j),
                          LibsemigroupsException);
        REQUIRE_THROWS_AS(T.fast_product(i, j + T.size()),
                          LibsemigroupsException);
        REQUIRE_THROWS_AS(T.fast_product(i + T.size(), j + T.size()),
                          LibsemigroupsException);
      }
    }
  }
}  // namespace libsemigroups
