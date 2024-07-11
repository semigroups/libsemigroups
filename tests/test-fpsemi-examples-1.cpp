// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2022 Murray T. Whyte
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

// This file is the first of three containing tests for the fpsemi-examples
// functions. The presentations here define finite semigroups, and we use
// ToddCoxeter in testing them. Exceptions and defaults are also checked in this
// file.

// #define CATCH_CONFIG_ENABLE_PAIR_STRINGMAKER

#include <cmath>                              // for pow
#include <cstddef>                            // for size_t
#include <unordered_map>                      // for operator!=
#include <vector>                             // for vector, allocator, oper...
                                              //
#include "catch_amalgamated.hpp"              // for StringRef, SourceLineInfo
#include "test-main.hpp"                      // for LIBSEMIGROUPS_TEST_CASE
                                              //
#include "libsemigroups/constants.hpp"        // for operator!=
#include "libsemigroups/exception.hpp"        // for LibsemigroupsException
#include "libsemigroups/fpsemi-examples.hpp"  // for author, symmetric_group
#include "libsemigroups/froidure-pin.hpp"     // for FroidurePin
#include "libsemigroups/obvinf.hpp"           // for is_obviously_infinite
#include "libsemigroups/presentation.hpp"     // for Presentation, operator==
#include "libsemigroups/ranges.hpp"           // for operator|, chain
#include "libsemigroups/to-froidure-pin.hpp"  // for to_froidure_pin
#include "libsemigroups/todd-coxeter.hpp"     // for ToddCoxeter
#include "libsemigroups/types.hpp"            // for congruence_kind, word_type
#include "libsemigroups/words.hpp"            // for operator""_w

#include "libsemigroups/detail/fmt.hpp"  // for format, print                                           //
#include "libsemigroups/detail/report.hpp"  // for ReportGuard

namespace libsemigroups {
  using literals::operator""_w;
  struct LibsemigroupsException;
  constexpr bool REPORT = false;

  using fpsemigroup::author;

  using fpsemigroup::alternating_group;
  using fpsemigroup::brauer_monoid;
  using fpsemigroup::chinese_monoid;
  using fpsemigroup::cyclic_inverse_monoid;
  using fpsemigroup::dual_symmetric_inverse_monoid;
  using fpsemigroup::fibonacci_semigroup;
  using fpsemigroup::full_transformation_monoid;
  using fpsemigroup::hypo_plactic_monoid;
  using fpsemigroup::monogenic_semigroup;
  using fpsemigroup::motzkin_monoid;
  using fpsemigroup::not_renner_type_B_monoid;
  using fpsemigroup::not_renner_type_D_monoid;
  using fpsemigroup::not_symmetric_group;
  using fpsemigroup::order_preserving_cyclic_inverse_monoid;
  using fpsemigroup::order_preserving_monoid;
  using fpsemigroup::orientation_preserving_monoid;
  using fpsemigroup::orientation_preserving_reversing_monoid;
  using fpsemigroup::partial_brauer_monoid;
  using fpsemigroup::partial_isometries_cycle_graph_monoid;
  using fpsemigroup::partial_transformation_monoid;
  using fpsemigroup::partition_monoid;
  using fpsemigroup::plactic_monoid;
  using fpsemigroup::rectangular_band;
  using fpsemigroup::renner_type_B_monoid;
  using fpsemigroup::renner_type_D_monoid;
  using fpsemigroup::sigma_stylic_monoid;
  using fpsemigroup::singular_brauer_monoid;
  using fpsemigroup::special_linear_group_2;
  using fpsemigroup::stellar_monoid;
  using fpsemigroup::stylic_monoid;
  using fpsemigroup::symmetric_group;
  using fpsemigroup::symmetric_inverse_monoid;
  using fpsemigroup::temperley_lieb_monoid;
  using fpsemigroup::uniform_block_bijection_monoid;
  using fpsemigroup::zero_rook_monoid;

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "000",
                          "test default values",
                          "[fpsemi-examples][quick]") {
    // author defaults
    REQUIRE(symmetric_group(4) == symmetric_group(4, author::Carmichael));
    REQUIRE(alternating_group(4) == alternating_group(4, author::Moore));
    REQUIRE(full_transformation_monoid(4)
            == full_transformation_monoid(4, author::Mitchell + author::Whyte));
    REQUIRE(
        partial_transformation_monoid(4)
        == partial_transformation_monoid(4, author::Mitchell + author::Whyte));
    REQUIRE(symmetric_inverse_monoid(4)
            == symmetric_inverse_monoid(4, author::Mitchell + author::Whyte));
    REQUIRE(dual_symmetric_inverse_monoid(4)
            == dual_symmetric_inverse_monoid(
                4, author::Easdown + author::East + author::FitzGerald));
    REQUIRE(uniform_block_bijection_monoid(4)
            == uniform_block_bijection_monoid(4, author::FitzGerald));
    REQUIRE(partition_monoid(4) == partition_monoid(4, author::East));
    REQUIRE(cyclic_inverse_monoid(4)
            == cyclic_inverse_monoid(4, author::Fernandes));

    // author::Any defaults
    REQUIRE(stellar_monoid(4, author::Any) == stellar_monoid(4));
    REQUIRE(singular_brauer_monoid(4, author::Any)
            == singular_brauer_monoid(4));
    REQUIRE(orientation_preserving_monoid(4, author::Any)
            == orientation_preserving_monoid(4));
    REQUIRE(orientation_preserving_reversing_monoid(4, author::Any)
            == orientation_preserving_reversing_monoid(4));
    REQUIRE(temperley_lieb_monoid(4, author::Any) == temperley_lieb_monoid(4));
    REQUIRE(brauer_monoid(4, author::Any) == brauer_monoid(4));
    REQUIRE(partial_brauer_monoid(4, author::Any) == partial_brauer_monoid(4));
    REQUIRE(motzkin_monoid(4, author::Any) == motzkin_monoid(4));
    REQUIRE(fibonacci_semigroup(5, 2, author::Any)
            == fibonacci_semigroup(5, 2));
    REQUIRE(plactic_monoid(4, author::Any) == plactic_monoid(4));
    REQUIRE(stylic_monoid(4, author::Any) == stylic_monoid(4));
    REQUIRE(rectangular_band(5, 3, author::Any) == rectangular_band(5, 3));
    REQUIRE(chinese_monoid(5, author::Any) == chinese_monoid(5));
    REQUIRE(monogenic_semigroup(6, 3, author::Any)
            == monogenic_semigroup(6, 3));
    REQUIRE(order_preserving_monoid(4, author::Any)
            == order_preserving_monoid(4));
    REQUIRE(order_preserving_cyclic_inverse_monoid(4, author::Any)
            == order_preserving_cyclic_inverse_monoid(4));
    REQUIRE(partial_isometries_cycle_graph_monoid(4, author::Any)
            == partial_isometries_cycle_graph_monoid(4));
    REQUIRE(special_linear_group_2(5, author::Any)
            == special_linear_group_2(5));
    REQUIRE(hypo_plactic_monoid(4, author::Any) == hypo_plactic_monoid(4));
    REQUIRE(sigma_stylic_monoid({3, 4}, author::Any)
            == sigma_stylic_monoid({3, 4}));
    REQUIRE(zero_rook_monoid(4, author::Any) == zero_rook_monoid(4));
    REQUIRE(renner_type_B_monoid(4, 1, author::Any)
            == renner_type_B_monoid(4, 1));
    REQUIRE(not_renner_type_B_monoid(4, 1, author::Any)
            == not_renner_type_B_monoid(4, 1));
    REQUIRE(renner_type_D_monoid(4, 1, author::Any)
            == renner_type_D_monoid(4, 1));
    REQUIRE(not_renner_type_D_monoid(4, 1, author::Any)
            == not_renner_type_D_monoid(4, 1));

    // index defaults
    REQUIRE(symmetric_group(4, author::Moore)
            == symmetric_group(4, author::Moore, 0));
    REQUIRE(cyclic_inverse_monoid(4, author::Fernandes)
            == cyclic_inverse_monoid(4, author::Fernandes, 1));
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "102",
                          "test semigroup/monoid status",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE(stellar_monoid(5).contains_empty_word());
    REQUIRE(dual_symmetric_inverse_monoid(5).contains_empty_word());
    REQUIRE(uniform_block_bijection_monoid(5).contains_empty_word());
    REQUIRE(partition_monoid(5, author::East).contains_empty_word());
    REQUIRE(!partition_monoid(3, author::Machine).contains_empty_word());
    REQUIRE(!singular_brauer_monoid(5).contains_empty_word());
    REQUIRE(orientation_preserving_monoid(5).contains_empty_word());
    REQUIRE(orientation_preserving_reversing_monoid(5).contains_empty_word());
    REQUIRE(temperley_lieb_monoid(5).contains_empty_word());
    REQUIRE(brauer_monoid(5).contains_empty_word());
    REQUIRE(partial_brauer_monoid(5).contains_empty_word());
    REQUIRE(motzkin_monoid(5).contains_empty_word());
    REQUIRE(!fibonacci_semigroup(5, 2).contains_empty_word());
    REQUIRE(plactic_monoid(5).contains_empty_word());
    REQUIRE(stylic_monoid(5).contains_empty_word());
    REQUIRE(symmetric_group(5, author::Burnside + author::Miller)
                .contains_empty_word());
    REQUIRE(symmetric_group(5, author::Carmichael).contains_empty_word());
    REQUIRE(symmetric_group(5, author::Coxeter + author::Moser)
                .contains_empty_word());
    REQUIRE(symmetric_group(5, author::Moore).contains_empty_word());
    REQUIRE(alternating_group(5).contains_empty_word());
    REQUIRE(!rectangular_band(5, 5).contains_empty_word());
    REQUIRE(
        full_transformation_monoid(5, author::Iwahori).contains_empty_word());
    REQUIRE(
        full_transformation_monoid(5, author::Aizenstat).contains_empty_word());
    REQUIRE(full_transformation_monoid(5, author::Mitchell + author::Whyte, 0)
                .contains_empty_word());
    REQUIRE(full_transformation_monoid(5, author::Mitchell + author::Whyte, 1)
                .contains_empty_word());

    REQUIRE(
        partial_transformation_monoid(5, author::Sutov).contains_empty_word());
    REQUIRE(partial_transformation_monoid(3, author::Machine)
                .contains_empty_word());
    REQUIRE(symmetric_inverse_monoid(5).contains_empty_word());
    REQUIRE(chinese_monoid(5).contains_empty_word());
    REQUIRE(monogenic_semigroup(0, 5).contains_empty_word());
    REQUIRE(!monogenic_semigroup(2, 6).contains_empty_word());
    REQUIRE(order_preserving_monoid(5).contains_empty_word());
    REQUIRE(cyclic_inverse_monoid(5).contains_empty_word());
    REQUIRE(order_preserving_cyclic_inverse_monoid(5).contains_empty_word());
    REQUIRE(partial_isometries_cycle_graph_monoid(5).contains_empty_word());
    REQUIRE(not_symmetric_group(5).contains_empty_word());
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "001",
                          "full_transformation_monoid author index except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(full_transformation_monoid(5, author::Burnside),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(full_transformation_monoid(5, author::Iwahori, 1),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(
        full_transformation_monoid(5, author::Mitchell + author::Whyte, 2),
        LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "002",
                          "full_transformation_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(full_transformation_monoid(3, author::Iwahori),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(full_transformation_monoid(1), LibsemigroupsException);
    REQUIRE_THROWS_AS(
        full_transformation_monoid(6, author::Mitchell + author::Whyte, 1),
        LibsemigroupsException);
    REQUIRE_THROWS_AS(full_transformation_monoid(1), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "003",
                          "partial_transformation_monoid author except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(partial_transformation_monoid(4, author::Burnside),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "004",
                          "partial_transformation_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(partial_transformation_monoid(3, author::Sutov),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(partial_transformation_monoid(1), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "005",
                          "symmetric_group auth except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(symmetric_group(9, author::Sutov),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "006",
                          "symmetric_group degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(symmetric_group(3, author::Carmichael),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "007",
                          "symmetric_group index except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(symmetric_group(5, author::Moore, 2),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(symmetric_group(5, author::Carmichael, 1),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "008",
                          "dual_symmetric_inverse_monoid auth except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(dual_symmetric_inverse_monoid(9, author::East),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "009",
                          "dual_symmetric_inverse_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(
        dual_symmetric_inverse_monoid(
            2, author::Easdown + author::East + author::FitzGerald),
        LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "010",
                          "uniform_block_bijection_monoid auth except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(uniform_block_bijection_monoid(9, author::Burnside),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "011",
                          "uniform_block_bijection_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(uniform_block_bijection_monoid(2, author::FitzGerald),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "012",
                          "partition_monoid auth except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(partition_monoid(3, author::Burnside),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(partition_monoid(5, author::Machine),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "013",
                          "partition_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(partition_monoid(4, author::Machine),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(partition_monoid(3, author::East),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "014",
                          "alternating_group auth except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(alternating_group(5, author::Sutov),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "015",
                          "alternating_group degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(alternating_group(3, author::Moore),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "016",
                          "chinese_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(chinese_monoid(1), LibsemigroupsException);
    REQUIRE_THROWS_AS(chinese_monoid(0), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "017",
                          "monogenic_semigroup degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(monogenic_semigroup(0, 0), LibsemigroupsException);
    REQUIRE_THROWS_AS(monogenic_semigroup(5, 0), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "018",
                          "rectangular_band degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(rectangular_band(0, 0), LibsemigroupsException);
    REQUIRE_THROWS_AS(rectangular_band(1, 0), LibsemigroupsException);
    REQUIRE_THROWS_AS(rectangular_band(0, 1), LibsemigroupsException);
    REQUIRE_THROWS_AS(rectangular_band(0, 100), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "019",
                          "stellar_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(stellar_monoid(0), LibsemigroupsException);
    REQUIRE_THROWS_AS(stellar_monoid(1), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "103",
                          "stellar_monoid auth except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(stellar_monoid(5, author::Sutov), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "020",
                          "plactic_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(plactic_monoid(0), LibsemigroupsException);
    REQUIRE_NOTHROW(plactic_monoid(1));
    REQUIRE_NOTHROW(plactic_monoid(2));
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "104",
                          "plactic_monoid auth except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(plactic_monoid(5, author::Iwahori),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "021",
                          "stylic_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(stylic_monoid(0), LibsemigroupsException);
    REQUIRE_THROWS_AS(stylic_monoid(1), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "105",
                          "stylic_monoid auth except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(stylic_monoid(5, author::Sutov), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "022",
                          "temperley_lieb_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(temperley_lieb_monoid(0), LibsemigroupsException);
    REQUIRE_THROWS_AS(temperley_lieb_monoid(1), LibsemigroupsException);
    REQUIRE_THROWS_AS(temperley_lieb_monoid(2), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "106",
                          "temperley_lieb_monoid auth except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(temperley_lieb_monoid(5, author::Sutov),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "023",
                          "singular_brauer_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(singular_brauer_monoid(0), LibsemigroupsException);
    REQUIRE_THROWS_AS(singular_brauer_monoid(1), LibsemigroupsException);
    REQUIRE_THROWS_AS(singular_brauer_monoid(2), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "107",
                          "singular_brauer_monoid auth except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(singular_brauer_monoid(5, author::Sutov),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "024",
                          "orientation_preserving_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(orientation_preserving_monoid(0), LibsemigroupsException);
    REQUIRE_THROWS_AS(orientation_preserving_monoid(1), LibsemigroupsException);
    REQUIRE_THROWS_AS(orientation_preserving_monoid(2), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "107",
                          "orientation_preserving_monoid auth except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(orientation_preserving_monoid(5, author::Sutov),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "025",
      "orientation_preserving_reversing_monoid degree except",
      "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(orientation_preserving_reversing_monoid(0),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(orientation_preserving_reversing_monoid(1),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(orientation_preserving_reversing_monoid(2),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "108",
                          "orientation_preserving_reversing_monoid auth except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(orientation_preserving_reversing_monoid(0, author::Sutov),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "026",
                          "order_preserving_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(order_preserving_monoid(0), LibsemigroupsException);
    REQUIRE_THROWS_AS(order_preserving_monoid(1), LibsemigroupsException);
    REQUIRE_THROWS_AS(order_preserving_monoid(2), LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "109",
                          "order_preserving_monoid auth except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(order_preserving_monoid(5, author::Sutov),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "027",
                          "cyclic_inverse_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(cyclic_inverse_monoid(0, author::Fernandes, 0),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(cyclic_inverse_monoid(1, author::Fernandes, 0),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(cyclic_inverse_monoid(2, author::Fernandes, 0),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(cyclic_inverse_monoid(0, author::Fernandes, 1),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(cyclic_inverse_monoid(1, author::Fernandes, 1),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(cyclic_inverse_monoid(2, author::Fernandes, 1),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "028",
                          "cyclic_inverse_monoid author except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(cyclic_inverse_monoid(5, author::Burnside, 0),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(cyclic_inverse_monoid(5, author::Fernandes, 3),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "029",
      "order_preserving_cyclic_inverse_monoid degree except",
      "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(order_preserving_cyclic_inverse_monoid(0),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(order_preserving_cyclic_inverse_monoid(1),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(order_preserving_cyclic_inverse_monoid(2),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "110",
      "order_preserving_cyclic_inverse_monoid author except",
      "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(order_preserving_cyclic_inverse_monoid(2, author::Sutov),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "030",
                          "partial_isometries_cycle_graph_monoid degree except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(partial_isometries_cycle_graph_monoid(0),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(partial_isometries_cycle_graph_monoid(1),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(partial_isometries_cycle_graph_monoid(2),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "111",
                          "partial_isometries_cycle_graph_monoid auth except",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    REQUIRE_THROWS_AS(partial_isometries_cycle_graph_monoid(5, author::Sutov),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "031",
                          "order_preserving_monoid(5)",
                          "[fpsemi-examples][quick]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 5;
    auto   p  = order_preserving_monoid(n);
    p.validate();

    ToddCoxeter tc(congruence_kind::twosided, p);
    REQUIRE(tc.number_of_classes() == 126);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "032",
                          "monogenic_semigroup(4, 9)",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    auto p  = monogenic_semigroup(4, 9);
    p.validate();
    ToddCoxeter tc(congruence_kind::twosided, p);
    REQUIRE(tc.number_of_classes() == 12);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "033",
                          "dual_symmetric_inverse_monoid(5)",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 5;
    auto   p  = dual_symmetric_inverse_monoid(
        n, author::Easdown + author::East + author::FitzGerald);
    p.validate();
    ToddCoxeter tc(congruence_kind::twosided, p);
    REQUIRE(tc.number_of_classes() == 6'721);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "034",
                          "uniform_block_bijection_monoid(5)",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 5;
    auto   p  = uniform_block_bijection_monoid(n, author::FitzGerald);
    p.validate();
    ToddCoxeter tc(congruence_kind::twosided, p);
    REQUIRE(tc.number_of_classes() == 1'496);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "035",
                          "partition_monoid(5)",
                          "[fpsemi-examples][standard]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided, partition_monoid(5));
    REQUIRE(tc.number_of_classes() == 115'975);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "036",
                          "singular_brauer_monoid(6)",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 6;
    ToddCoxeter tc(congruence_kind::twosided, singular_brauer_monoid(n));
    REQUIRE(tc.number_of_classes() == 9'675);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "037",
                          "orientation_preserving_monoid(6)",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 6;
    ToddCoxeter tc(congruence_kind::twosided, orientation_preserving_monoid(n));
    REQUIRE(tc.number_of_classes() == 2'742);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "038",
                          "orientation_preserving_reversing_monoid(5)",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 5;
    ToddCoxeter tc(congruence_kind::twosided,
                   orientation_preserving_reversing_monoid(n));
    REQUIRE(tc.number_of_classes() == 1'015);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "039",
                          "partition_monoid(3)",
                          "[fpsemi-examples][quick]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   partition_monoid(3, author::Machine));
    REQUIRE(tc.number_of_classes() == 203);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "040",
                          "temperley_lieb_monoid(10)",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 10;
    ToddCoxeter tc(congruence_kind::twosided, temperley_lieb_monoid(n));
    REQUIRE(tc.number_of_classes() == 16'796);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "041",
                          "brauer_monoid(5)",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 5;
    ToddCoxeter tc(congruence_kind::twosided, brauer_monoid(n));
    REQUIRE(tc.number_of_classes() == 945);

    tc.init(congruence_kind::twosided, brauer_monoid(1));
    REQUIRE(tc.number_of_classes() == 1);

    tc.init(congruence_kind::twosided, brauer_monoid(2));
    REQUIRE(tc.number_of_classes() == 3);

    tc.init(congruence_kind::twosided, brauer_monoid(3));
    REQUIRE(tc.number_of_classes() == 15);

    tc.init(congruence_kind::twosided, brauer_monoid(4));
    REQUIRE(tc.number_of_classes() == 105);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "042",
                          "rectangular_band(5, 9)",
                          "[fpsemi-examples][quick]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided, rectangular_band(5, 9));
    REQUIRE(tc.number_of_classes() == 45);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "043",
                          "symmetric_group(6) Burnside + Miller",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 6;
    ToddCoxeter tc(congruence_kind::twosided,
                   symmetric_group(n, author::Burnside + author::Miller));
    REQUIRE(tc.number_of_classes() == 720);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "044",
                          "symmetric_group(6) Carmichael",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 6;
    ToddCoxeter tc(congruence_kind::twosided,
                   symmetric_group(n, author::Carmichael));
    REQUIRE(tc.number_of_classes() == 720);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "045",
                          "symmetric_group(6) Moore index 0",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 6;

    ToddCoxeter tc(congruence_kind::twosided,
                   symmetric_group(n, author::Moore, 0));
    REQUIRE(tc.number_of_classes() == 720);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "046",
                          "symmetric_group(7) Moore index 1",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 7;
    ToddCoxeter tc(congruence_kind::twosided,
                   symmetric_group(n, author::Moore, 1));
    REQUIRE(tc.number_of_classes() == 5'040);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "047",
                          "symmetric_group(6) Coxeter + Moser",
                          "[fpsemi-examples][quick]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 6;

    ToddCoxeter tc(congruence_kind::twosided,
                   symmetric_group(n, author::Coxeter + author::Moser));
    REQUIRE(tc.number_of_classes() == 720);
    REQUIRE(presentation::length(
                symmetric_group(n, author::Coxeter + author::Moser))
            == 58);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "048",
                          "symmetric inverse monoid Gay presentation",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    auto p  = symmetric_inverse_monoid(4, author::Gay);

    presentation::sort_each_rule(p);
    presentation::sort_rules(p);

    REQUIRE(p.rules
            == std::vector<word_type>(
                {00_w,  {},    11_w,   {},    22_w,   {},    33_w,  3_w,
                 20_w,  02_w,  31_w,   13_w,  32_w,   23_w,  101_w, 010_w,
                 212_w, 121_w, 0303_w, 303_w, 3030_w, 0303_w}));

    ToddCoxeter tc(congruence_kind::twosided, p);
    REQUIRE(tc.number_of_classes() == 209);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "073",
                          "0-Hecke rook monoid",
                          "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);
    auto p  = zero_rook_monoid(4);

    presentation::sort_each_rule(p);
    presentation::sort_rules(p);

    ToddCoxeter tc(congruence_kind::twosided, p);
    REQUIRE(tc.number_of_classes() == 209);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "049",
                          "full_transformation_monoid(5) Iwahori",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 5;

    ToddCoxeter tc(congruence_kind::twosided,
                   full_transformation_monoid(n, author::Iwahori));
    REQUIRE(tc.number_of_classes() == 3'125);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "117",
      "full_transformation_monoid(n = 2, 3) Mitchell + Whyte",
      "[fpsemi-examples][quick]") {
    auto rg = ReportGuard(REPORT);

    ToddCoxeter tc2(
        congruence_kind::twosided,
        full_transformation_monoid(2, author::Mitchell + author::Whyte));
    REQUIRE(tc2.number_of_classes() == 4);

    ToddCoxeter tc3(
        congruence_kind::twosided,
        full_transformation_monoid(3, author::Mitchell + author::Whyte));
    REQUIRE(tc3.number_of_classes() == 27);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "112",
                          "full_transformation_monoid(5) Mitchell + Whyte",
                          "[fpsemi-examples][quick]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 5;

    ToddCoxeter tc(
        congruence_kind::twosided,
        full_transformation_monoid(n, author::Mitchell + author::Whyte));
    REQUIRE(tc.number_of_classes() == 3'125);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "118",
      "full_transformation_monoid(5) Mitchell + Whyte index 1",
      "[fpsemi-examples]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 5;

    ToddCoxeter tc(
        congruence_kind::twosided,
        full_transformation_monoid(n, author::Mitchell + author::Whyte, 1));
    REQUIRE(tc.number_of_classes() == 3'125);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "113",
                          "full_transformation_monoid(6) Mitchell + Whyte",
                          "[fpsemi-examples]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 6;

    ToddCoxeter tc(
        congruence_kind::twosided,
        full_transformation_monoid(n, author::Mitchell + author::Whyte));
    REQUIRE(tc.number_of_classes() == 46'656);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "050",
                          "partial_transformation_monoid(5) Sutov",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 5;
    ToddCoxeter tc(congruence_kind::twosided,
                   partial_transformation_monoid(n, author::Sutov));
    REQUIRE(tc.number_of_classes() == 7'776);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "116",
                          "partial_transformation_monoid(5) Mitchell + Whyte",
                          "[fpsemi-examples][quick]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 5;
    ToddCoxeter tc(
        congruence_kind::twosided,
        partial_transformation_monoid(n, author::Mitchell + author::Whyte));
    REQUIRE(tc.number_of_classes() == 7'776);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "051",
                          "partial_transformation_monoid(3) Machine",
                          "[fpsemi-examples][quick]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   partial_transformation_monoid(3, author::Machine));
    REQUIRE(tc.number_of_classes() == 64);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "052",
                          "symmetric_inverse_monoid(5)",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 5;

    ToddCoxeter tc(congruence_kind::twosided,
                   symmetric_inverse_monoid(n, author::Sutov));
    REQUIRE(tc.number_of_classes() == 1'546);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "114",
                          "symmetric_inverse_monoid(4) Mitchell + Whyte",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 4;

    ToddCoxeter tc(
        congruence_kind::twosided,
        symmetric_inverse_monoid(n, author::Mitchell + author::Whyte));
    REQUIRE(tc.number_of_classes() == 209);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "115",
                          "symmetric_inverse_monoid(5) Mitchell + Whyte",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 5;

    ToddCoxeter tc(
        congruence_kind::twosided,
        symmetric_inverse_monoid(n, author::Mitchell + author::Whyte));
    REQUIRE(tc.number_of_classes() == 1'546);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "053",
                          "fibonacci_semigroup(2, 5)",
                          "[fpsemi-examples][quick]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided, fibonacci_semigroup(2, 5));
    REQUIRE(tc.number_of_classes() == 11);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "054",
                          "alternating_group(7) Moore",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 7;
    ToddCoxeter tc(congruence_kind::twosided,
                   alternating_group(n, author::Moore));
    REQUIRE(tc.number_of_classes() == 2'520);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "055",
                          "full_transformation_monoid(4) Aizenstat",
                          "[fpsemi-examples][quick]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 4;

    ToddCoxeter tc(congruence_kind::twosided,
                   full_transformation_monoid(n, author::Aizenstat));
    REQUIRE(tc.number_of_classes() == 256);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "056",
                          "order_preserving_monoid(10)",
                          "[fpsemi-examples][standard]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 10;
    ToddCoxeter tc(congruence_kind::twosided, order_preserving_monoid(n));
    REQUIRE(tc.number_of_classes() == 92'378);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "057",
                          "cyclic_inverse_monoid(4) Fernandes 1",
                          "[fpsemi-examples][quick]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 4;

    ToddCoxeter tc(congruence_kind::twosided,
                   cyclic_inverse_monoid(n, author::Fernandes, 1));
    REQUIRE(tc.number_of_classes() == 61);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "058",
                          "cyclic_inverse_monoid(8) Fernandes index 1",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 8;
    ToddCoxeter tc(congruence_kind::twosided,
                   cyclic_inverse_monoid(n, author::Fernandes, 1));
    REQUIRE(tc.number_of_classes() == 2'041);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "059",
                          "cyclic_inverse_monoid Fernandes index 0",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto rg = ReportGuard(REPORT);
    for (size_t n = 3; n < 10; ++n) {
      auto p = fpsemigroup::cyclic_inverse_monoid(n, author::Fernandes, 0);
      REQUIRE(p.rules.size() == (n * n + 3 * n + 4));
      if (n == 3) {
        presentation::sort_each_rule(p);
        presentation::sort_rules(p);

        REQUIRE(p.rules
                == std::vector<word_type>(
                    {000_w, {},   11_w, 1_w,  22_w,   2_w,  33_w, 3_w,
                     10_w,  02_w, 20_w, 03_w, 21_w,   12_w, 30_w, 01_w,
                     31_w,  13_w, 32_w, 23_w, 0123_w, 123_w}));
      }

      ToddCoxeter tc(congruence_kind::twosided, p);
      REQUIRE(tc.number_of_classes() == n * std::pow(2, n) - n + 1);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "060",
                          "order_preserving_cyclic_inverse_monoid(4)",
                          "[fpsemi-examples][quick]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 4;
    ToddCoxeter tc(congruence_kind::twosided,
                   order_preserving_cyclic_inverse_monoid(n));
    REQUIRE(tc.number_of_classes() == 38);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "061",
                          "order_preserving_cyclic_inverse_monoid(11)",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 11;
    ToddCoxeter tc(congruence_kind::twosided,
                   order_preserving_cyclic_inverse_monoid(n));
    REQUIRE(tc.number_of_classes() == 6'120);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "062",
                          "order_preserving_cyclic_inverse_monoid(11)",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 11;
    ToddCoxeter tc(congruence_kind::twosided,
                   order_preserving_cyclic_inverse_monoid(n));
    REQUIRE(tc.number_of_classes() == 6'120);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "063",
                          "partial_isometries_cycle_graph_monoid(5)",
                          "[fpsemi-examples][quick]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 5;
    ToddCoxeter tc(congruence_kind::twosided,
                   partial_isometries_cycle_graph_monoid(n));
    REQUIRE(tc.number_of_classes() == 286);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "064",
                          "partial_isometries_cycle_graph_monoid(4)",
                          "[fpsemi-examples][quick]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 4;
    ToddCoxeter tc(congruence_kind::twosided,
                   partial_isometries_cycle_graph_monoid(n));
    REQUIRE(tc.number_of_classes() == 97);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "065",
      "partial_isometries_cycle_graph_monoid(10)",
      "[fpsemi-examples][quick][no-valgrind][no-coverage]") {
    auto        rg = ReportGuard(REPORT);
    size_t      n  = 10;
    ToddCoxeter tc(congruence_kind::twosided,
                   partial_isometries_cycle_graph_monoid(n));
    REQUIRE(tc.number_of_classes() == 20'311);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "066",
                          "not_symmetric_group(4)",
                          "[fpsemi-examples][quick]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 4;

    ToddCoxeter tc(congruence_kind::twosided, not_symmetric_group(n));
    REQUIRE(tc.number_of_classes() == 72);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "072",
                          "partial_brauer_monoid(4)",
                          "[fpsemi-examples][quick][no-valgrind]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 4;

    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::partial_brauer_monoid(n));
    REQUIRE(tc.number_of_classes() == 764);
    REQUIRE_THROWS_AS(fpsemigroup::partial_brauer_monoid(0),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("fpsemi-examples",
                          "220",
                          "motzkin_monoid(5)",
                          "[fpsemi-examples][quick]") {
    auto   rg = ReportGuard(REPORT);
    size_t n  = 5;

    ToddCoxeter tc(congruence_kind::twosided, fpsemigroup::motzkin_monoid(n));
    REQUIRE(tc.number_of_classes() == 2188);
    REQUIRE_THROWS_AS(fpsemigroup::motzkin_monoid(0), LibsemigroupsException);

    tc.init(congruence_kind::twosided, motzkin_monoid(1));
    REQUIRE(tc.number_of_classes() == 2);

    tc.init(congruence_kind::twosided, motzkin_monoid(2));
    REQUIRE(tc.number_of_classes() == 9);

    tc.init(congruence_kind::twosided, motzkin_monoid(3));
    REQUIRE(tc.number_of_classes() == 51);

    tc.init(congruence_kind::twosided, motzkin_monoid(4));
    REQUIRE(tc.number_of_classes() == 323);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "074",
      "not_renner_type_B_monoid(2, 1) (Godelle presentation)",
      "[quick][fpsemi-examples][hivert]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::not_renner_type_B_monoid(2, 1));
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(!tc.started());
    REQUIRE(!tc.finished());
    REQUIRE(tc.number_of_classes() == 57);
    REQUIRE(tc.started());
    REQUIRE(tc.finished());
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "075",
      "not_renner_type_B_monoid(2, 0) (Godelle presentation)",
      "[quick][fpsemi-examples][hivert]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::not_renner_type_B_monoid(2, 0));
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 57);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "078",
      "renner_type_B_monoid(2, 1) (Gay-Hivert presentation)",
      "[quick][fpsemi-examples][hivert]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_B_monoid(2, 1));
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 57);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "079",
      "renner_type_B_monoid(2, 0) (Gay-Hivert presentation)",
      "[quick][fpsemi-examples][hivert]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_B_monoid(2, 0));
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 57);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "080",
      "renner_type_B_monoid(3, 1) (Gay-Hivert presentation)",
      "[quick][fpsemi-examples][hivert][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_B_monoid(3, 1));
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 757);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "081",
      "renner_type_B_monoid(3, 0) (Gay-Hivert presentation)",
      "[quick][fpsemi-examples][hivert][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_B_monoid(3, 0));
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 757);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "082",
      "renner_type_B_monoid(4, 1) (Gay-Hivert presentation)",
      "[quick][fpsemi-examples][hivert][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_B_monoid(4, 1));
    REQUIRE(tc.presentation().rules.size() == 186);
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 13'889);
    // TODO REQUIRE(to_froidure_pin(tc).run().number_of_rules() == 356);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "083",
      "renner_type_B_monoid(4, 0) (Gay-Hivert presentation)",
      "[quick][fpsemi-examples][hivert][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_B_monoid(4, 0));
    REQUIRE(tc.presentation().rules.size() == 186);
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 13'889);
    // REQUIRE(S.froidure_pin()->number_of_rules() == 356);
  }

  // This appears to be an example where KB + FP is faster than TC
  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "084",
      "renner_type_B_monoid(5, 1) (Gay-Hivert presentation)",
      "[extreme][fpsemi-examples][hivert]") {
    auto        rg = ReportGuard(true);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_B_monoid(5, 1));
    REQUIRE(tc.presentation().rules.size() == 272);
    REQUIRE(!is_obviously_infinite(tc));
    tc.strategy(ToddCoxeter::options::strategy::felsch);
    REQUIRE(tc.number_of_classes() == 322'021);
    // REQUIRE(S.froidure_pin()->number_of_rules() == 1453);
    // {
    //   congruence::ToddCoxeter tc(
    //       twosided,
    //       S.froidure_pin(),
    //       congruence::ToddCoxeter::options::froidure_pin::use_cayley_graph);
    //   REQUIRE(tc.number_of_classes() == 322021);  // Works!
    // }
    // {
    //   fpsemigroup::ToddCoxeter tc(S.froidure_pin());
    //   REQUIRE(tc.number_of_rules() == 1453);
    //   REQUIRE(tc.size() == 322021);
    // }
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "085",
      "renner_type_B_monoid(5, 0) (Gay-Hivert presentation)",
      "[extreme][fpsemi-examples][hivert]") {
    auto        rg = ReportGuard(true);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_B_monoid(5, 0));
    tc.strategy(ToddCoxeter::options::strategy::felsch);

    REQUIRE(tc.presentation().rules.size() == 272);
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 322'021);

    auto S = to_froidure_pin(tc);
    S.run();
    REQUIRE(S.number_of_rules() == 1430);

    //  FIXME the following doesn't seem to work
    // tc.init(congruence_kind::twosided, S.right_cayley_graph());
    // REQUIRE(tc.number_of_classes() == 322'021);  // Works!
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "092",
      "renner_type_D_monoid(2, 1) (Gay-Hivert presentation)",
      "[quick][fpsemi-examples][hivert]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_D_monoid(2, 1));
    REQUIRE(tc.presentation().rules.size() == 68);
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 37);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "093",
      "renner_type_D_monoid(2, 0) (Gay-Hivert presentation)",
      "[quick][fpsemi-examples][hivert]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_D_monoid(2, 0));
    REQUIRE(tc.presentation().rules.size() == 68);
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 37);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "094",
      "renner_type_D_monoid(3, 1) (Gay-Hivert presentation)",
      "[quick][fpsemi-examples][hivert][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_D_monoid(3, 1));
    REQUIRE(tc.presentation().rules.size() == 130);
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 541);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "095",
      "renner_type_D_monoid(3, 0) (Gay-Hivert presentation)",
      "[quick][fpsemi-examples][hivert][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_D_monoid(3, 0));
    REQUIRE(tc.presentation().rules.size() == 130);
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 541);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "096",
      "renner_type_D_monoid(4, 1) (Gay-Hivert presentation)",
      "[quick][fpsemi-examples][hivert][no-valgrind]") {
    auto        rg = ReportGuard(REPORT);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_D_monoid(4, 1));
    REQUIRE(tc.presentation().rules.size() == 208);
    REQUIRE(!is_obviously_infinite(tc));

    REQUIRE(tc.number_of_classes() == 10'625);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "097",
      "renner_type_D_monoid(4, 0) (Gay-Hivert presentation)",
      "[quick][fpsemi-examples][hivert][no-valgrind]") {
    auto        rg = ReportGuard(false);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_D_monoid(4, 0));
    REQUIRE(tc.presentation().rules.size() == 208);
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 10'625);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "098",
      "renner_type_D_monoid(5, 1) (Gay-Hivert presentation)",
      "[extreme][fpsemi-examples][hivert]") {
    auto        rg = ReportGuard(true);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_D_monoid(5, 1));
    tc.strategy(ToddCoxeter::options::strategy::felsch);

    REQUIRE(tc.presentation().rules.size() == 302);
    REQUIRE(!is_obviously_infinite(tc));

    REQUIRE(tc.number_of_classes() == 258'661);
  }

  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "099",
      "renner_type_D_monoid(5, 0) (Gay-Hivert presentation)",
      "[extreme][fpsemi-examples][hivert]") {
    auto        rg = ReportGuard(true);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_D_monoid(5, 0));
    tc.strategy(ToddCoxeter::options::strategy::felsch);

    REQUIRE(tc.presentation().rules.size() == 302);
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 258'661);
  }

  // Takes about 1m46s
  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "100",
      "renner_type_D_monoid(6, 1) (Gay-Hivert presentation)",
      "[extreme][fpsemi-examples][hivert]") {
    auto        rg = ReportGuard(true);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_D_monoid(6, 1));
    tc.strategy(ToddCoxeter::options::strategy::felsch);

    REQUIRE(tc.presentation().rules.size() == 412);
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 7'464'625);
  }

  // Takes about 2m42s
  LIBSEMIGROUPS_TEST_CASE(
      "fpsemi-examples",
      "101",
      "renner_type_D_monoid(6, 0) (Gay-Hivert presentation)",
      "[extreme][fpsemi-examples][hivert]") {
    auto        rg = ReportGuard(true);
    ToddCoxeter tc(congruence_kind::twosided,
                   fpsemigroup::renner_type_D_monoid(6, 0));
    tc.strategy(ToddCoxeter::options::strategy::felsch);

    REQUIRE(tc.presentation().rules.size() == 412);
    REQUIRE(!is_obviously_infinite(tc));
    REQUIRE(tc.number_of_classes() == 7'464'625);
  }
}  // namespace libsemigroups
