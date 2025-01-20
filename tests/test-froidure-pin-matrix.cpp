//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2020 James D. Mitchell
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

#include "catch_amalgamated.hpp"  // REQUIRE
#include "test-main.hpp"          // for LIBSEMIGROUPS_TEST_CASE

#include "libsemigroups/constants.hpp"     // for UNDEFINED
#include "libsemigroups/froidure-pin.hpp"  // for FroidurePin
#include "libsemigroups/matrix.hpp"        // for NTPSemiring, MaxPlusTruncS...

#include "libsemigroups/detail/report.hpp"  // for ReportGuard

namespace libsemigroups {
  struct LibsemigroupsException;
  constexpr bool REPORT = false;

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("FroidurePin",
                                   "040",
                                   "Example 000",
                                   "[quick][froidure-pin][matrix]",
                                   MaxPlusMat<2>,
                                   MaxPlusMat<>) {
    auto                  rg = ReportGuard(REPORT);
    FroidurePin<TestType> S;
    S.add_generator(make<TestType>({{0, -4}, {-4, -1}}));
    S.add_generator(make<TestType>({{0, -3}, {-3, -1}}));

    REQUIRE(S.size() == 26);
    REQUIRE(S.degree() == 2);
    REQUIRE(S.number_of_idempotents() == 4);
    REQUIRE(S.number_of_generators() == 2);
    REQUIRE(S.number_of_rules() == 9);
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
    REQUIRE(S.position(x) == 5);
    REQUIRE(S.contains(x));

    if constexpr (IsDynamicMatrix<TestType>) {
      // If TestType is a static matrix, then the next line leads to out of
      // bounds accesses, since we are constructing a too big matrix without
      // checks.
      x = TestType({{-2, 2, 0}, {-1, 0, 0}, {0, 0, 0}});
      REQUIRE(S.position(x) == UNDEFINED);
      REQUIRE(!S.contains(x));
    }
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("FroidurePin",
                                   "042",
                                   "Example 001",
                                   "[extreme][froidure-pin][matrix]",
                                   (NTPMat<0, 6, 3>),
                                   (NTPMat<0, 6>),
                                   NTPMat<>) {
    auto rg = ReportGuard();

    NTPSemiring<>* sr = nullptr;
    if constexpr (!std::is_same_v<typename TestType::semiring_type, void>) {
      sr = new NTPSemiring<>(0, 6);
    }
    FroidurePin<TestType> S;
    S.add_generator(TestType(sr, {{0, 0, 1}, {0, 1, 0}, {1, 1, 0}}));
    S.add_generator(TestType(sr, {{0, 0, 1}, {0, 1, 0}, {2, 0, 0}}));
    S.add_generator(TestType(sr, {{0, 0, 1}, {0, 1, 1}, {1, 0, 0}}));
    S.add_generator(TestType(sr, {{0, 0, 1}, {0, 1, 0}, {3, 0, 0}}));
    S.reserve(10'077'696);
    REQUIRE(S.size() == 10'077'696);
    REQUIRE(S.number_of_idempotents() == 13'688);
    delete sr;
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("FroidurePin",
                                   "044",
                                   "Example 004",
                                   "[quick][froidure-pin][matrix]",
                                   MinPlusMat<>,
                                   MinPlusMat<2>) {
    auto                  rg = ReportGuard(REPORT);
    FroidurePin<TestType> S;
    S.add_generator(make<TestType>({{1, 0}, {0, POSITIVE_INFINITY}}));

    REQUIRE(S.size() == 3);
    REQUIRE(S.degree() == 2);
    REQUIRE(S.number_of_idempotents() == 1);
    REQUIRE(S.number_of_generators() == 1);
    REQUIRE(S.number_of_rules() == 1);

    REQUIRE(S[0] == S.generator(0));
    REQUIRE(S.position(S.generator(0)) == 0);
    REQUIRE(S.contains(S.generator(0)));

    auto x = TestType({{-2, 2}, {-1, 0}});
    REQUIRE(S.position(x) == UNDEFINED);
    REQUIRE(!S.contains(x));
    x.product_inplace_no_checks(S.generator(0), S.generator(0));
    REQUIRE(S.position(x) == 1);
    REQUIRE(S.contains(x));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("FroidurePin",
                                   "046",
                                   "Example 005",
                                   "[quick][froidure-pin][matrix]",
                                   (MaxPlusTruncMat<33, 3>),
                                   MaxPlusTruncMat<33>,
                                   MaxPlusTruncMat<>) {
    auto rg = ReportGuard(REPORT);

    MaxPlusTruncSemiring<> const* sr = nullptr;
    if constexpr (!std::is_same_v<typename TestType::semiring_type, void>) {
      sr = new MaxPlusTruncSemiring<>(33);
    }
    FroidurePin<TestType> S;
    S.add_generator(make<TestType>(sr, {{22, 21, 0}, {10, 0, 0}, {1, 32, 1}}));
    S.add_generator(make<TestType>(sr, {{0, 0, 0}, {0, 1, 0}, {1, 1, 0}}));

    REQUIRE(S.size() == 119);
    REQUIRE(S.degree() == 3);
    REQUIRE(S.number_of_idempotents() == 1);
    REQUIRE(S.number_of_generators() == 2);
    REQUIRE(S.number_of_rules() == 18);

    REQUIRE(S[0] == S.generator(0));
    REQUIRE(S.position(S.generator(0)) == 0);
    REQUIRE(S.contains(S.generator(0)));

    auto x = make<TestType>(sr, {{2, 2}, {1, 0}});
    REQUIRE(S.position(x) == UNDEFINED);
    delete sr;
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("FroidurePin",
                                   "049",
                                   "Example 006",
                                   "[quick][froidure-pin][matrix]",
                                   (MinPlusTruncMat<11, 3>),
                                   MinPlusTruncMat<11>,
                                   MinPlusTruncMat<>) {
    auto rg = ReportGuard(REPORT);

    MinPlusTruncSemiring<>* sr = nullptr;
    if constexpr (!std::is_same_v<typename TestType::semiring_type, void>) {
      sr = new MinPlusTruncSemiring(11);
    }
    FroidurePin<TestType> S;
    S.add_generator(make<TestType>(sr, {{2, 1, 0}, {10, 0, 0}, {1, 2, 1}}));
    S.add_generator(make<TestType>(sr, {{10, 0, 0}, {0, 1, 0}, {1, 1, 0}}));

    REQUIRE(S.size() == 1039);
    REQUIRE(S.degree() == 3);
    REQUIRE(S.number_of_idempotents() == 5);
    REQUIRE(S.number_of_generators() == 2);
    REQUIRE(S.number_of_rules() == 38);

    REQUIRE(S[0] == S.generator(0));
    REQUIRE(S.position(S.generator(0)) == 0);
    REQUIRE(S.contains(S.generator(0)));

    auto x = make<TestType>(sr, {{2, 2, 0}, {1, 0, 0}, {0, 0, 0}});
    REQUIRE(S.position(x) == UNDEFINED);
    REQUIRE(!S.contains(x));
    x.product_inplace_no_checks(S.generator(0), S.generator(0));
    REQUIRE(S.position(x) == 2);
    REQUIRE(S.contains(x));
    delete sr;
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("FroidurePin",
                                   "052",
                                   "Example 007",
                                   "[quick][froidure-pin][matrix]",
                                   (NTPMat<11, 3>),
                                   NTPMat<>) {
    auto rg = ReportGuard(REPORT);

    NTPSemiring<>* sr = nullptr;
    if constexpr (!std::is_same_v<typename TestType::semiring_type, void>) {
      sr = new NTPSemiring<>(11, 3);
    }
    FroidurePin<TestType> S;
    S.add_generator(make<TestType>(sr, {{2, 1, 0}, {10, 0, 0}, {1, 2, 1}}));
    S.add_generator(make<TestType>(sr, {{10, 0, 0}, {0, 1, 0}, {1, 1, 0}}));

    REQUIRE(S.size() == 86);
    REQUIRE(S.degree() == 3);
    REQUIRE(S.number_of_idempotents() == 10);
    REQUIRE(S.number_of_generators() == 2);
    REQUIRE(S.number_of_rules() == 16);

    REQUIRE(S[0] == S.generator(0));
    REQUIRE(S.position(S.generator(0)) == 0);
    REQUIRE(S.contains(S.generator(0)));

    auto x = make<TestType>(sr, {{2, 2, 0}, {1, 0, 0}, {0, 0, 0}});
    REQUIRE(S.position(x) == UNDEFINED);
    REQUIRE(!S.contains(x));
    x.product_inplace_no_checks(S.generator(1), S.generator(0));
    REQUIRE(S.position(x) == 4);
    REQUIRE(S.contains(x));
    delete sr;
  }

}  // namespace libsemigroups
