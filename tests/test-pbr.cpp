//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2019-21 James D. Mitchell
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

#include "catch_amalgamated.hpp"  // for TEST_CASE
#include "test-main.hpp"          // LIBSEMIGROUPS_TEST_CASE

#include "libsemigroups/exception.hpp"  // for LibsemigroupsException
#include "libsemigroups/pbr.hpp"        // for PBR

#include "libsemigroups/detail/string.hpp"  // for to_string

namespace libsemigroups {
  // Forward decl
  struct LibsemigroupsException;

  LIBSEMIGROUPS_TEST_CASE("PBR",
                          "000",
                          "universal product with convenience constructor",
                          "[quick][pbr]") {
    std::vector<std::vector<int32_t>> const left
        = {{-3, -1}, {-3, -2, -1, 1, 2, 3}, {-3, -2, -1, 1, 3}};

    std::vector<std::vector<int32_t>> const right
        = {{-3, -1, 1, 2, 3}, {-3, 1, 3}, {-3, -2, -1, 2, 3}};
    PBR x(left, right);

    PBR y({{-3, -2, -1, 1}, {-3, -2, 3}, {-3, 2, 3}},
          {{-3, -2, -1, 3}, {-3, -2, -1, 3}, {-2, 2, 3}});

    PBR z({{-3, -2, -1, 1}, {-3, -2, 3}, {-3, 2, 3}},
          {{-3, -2, -1, 3}, {-3, -2, -1, 3}, {-2, 2, 3}});

    PBR xx({{3, 5},
            {0, 1, 2, 3, 4, 5},
            {0, 2, 3, 4, 5},
            {0, 1, 2, 3, 5},
            {0, 2, 5},
            {1, 2, 3, 4, 5}});
    PBR yy({{0, 3, 4, 5},
            {2, 4, 5},
            {1, 2, 5},
            {2, 3, 4, 5},
            {2, 3, 4, 5},
            {1, 2, 4}});

    REQUIRE(x == xx);
    REQUIRE(y == yy);
    REQUIRE(x.degree() == 3);

    z.product_inplace(x, y);

    PBR expected({{0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5}});
    REQUIRE(z == expected);
  }

  LIBSEMIGROUPS_TEST_CASE("PBR", "001", "universal product", "[quick][pbr]") {
    PBR x({{3, 5},
           {0, 1, 2, 3, 4, 5},
           {0, 2, 3, 4, 5},
           {0, 1, 2, 3, 5},
           {0, 2, 5},
           {1, 2, 3, 4, 5}});
    PBR y({{0, 3, 4, 5},
           {2, 4, 5},
           {1, 2, 5},
           {2, 3, 4, 5},
           {2, 3, 4, 5},
           {1, 2, 4}});

    PBR z({{0, 3, 4, 5},
           {2, 4, 5},
           {1, 2, 5},
           {2, 3, 4, 5},
           {2, 3, 4, 5},
           {1, 2, 4}});
    z.product_inplace(x, y);

    PBR expected({{0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5}});
    REQUIRE(z == expected);
  }

  LIBSEMIGROUPS_TEST_CASE("PBR",
                          "002",
                          "product [bigger than previous]",
                          "[quick][pbr]") {
    PBR x({{3, 5},
           {0, 1, 2, 3, 4, 5},
           {0, 2, 3, 4, 5},
           {0, 1, 2, 3, 5},
           {0, 2, 5},
           {1, 2, 3, 4, 5},
           {},
           {}});
    PBR y({{3, 5},
           {0, 1, 2, 3, 4, 5},
           {0, 2, 3, 4, 5},
           {0, 1, 2, 3, 5},
           {0, 2, 5},
           {1, 2, 3, 4, 5},
           {},
           {6}});
    x.product_inplace(y, y);
    PBR expected({{0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {},
                  {6}});

    REQUIRE(x == expected);

    x = PBR({{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {7}});
    y = PBR({{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {7}});

    x.product_inplace(y, y);
    expected = PBR(
        {{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {7}});
    REQUIRE(x == expected);
  }

  LIBSEMIGROUPS_TEST_CASE("PBR", "003", "hash", "[quick][pbr][no-valgrind]") {
    PBR x({{1}, {4}, {3}, {1}, {0, 2}, {0, 3, 4, 5}});
    for (size_t i = 0; i < 1000000; i++) {
      x.hash_value();
    }
  }

  LIBSEMIGROUPS_TEST_CASE("PBR", "004", "delete/copy", "[quick][pbr]") {
    PBR x = make<PBR>({{1}, {4}, {3}, {1}, {0, 2}, {0, 3, 4, 5}});
    PBR y(x);
    REQUIRE(x == y);
    PBR z({{1}, {4}, {3}, {1}, {0, 2}, {0, 3, 4, 5}});
    REQUIRE(y == z);
    PBR yy(y);
    REQUIRE(yy == y);
    PBR zz(yy);
    PBR a({{1}, {4}, {3}, {1}, {0, 2}, {0, 3, 4, 5}});
    REQUIRE(z == a);
    REQUIRE(zz == a);
    PBR tt(std::move(zz));
    REQUIRE(tt == a);
    tt = z;
    REQUIRE(tt == z);
  }

  LIBSEMIGROUPS_TEST_CASE("PBR", "005", "exceptions", "[quick][pbr]") {
    REQUIRE_THROWS_AS(make<PBR>({{1}, {4}, {3}, {10}, {0, 2}, {0, 3, 4, 5}}),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(make<PBR>({{4}, {3}, {0}, {0, 2}, {0, 3, 4, 5}}),
                      LibsemigroupsException);

    REQUIRE_NOTHROW(PBR({{-3, -1}, {-3, -2, -1, 1, 2, 3}, {-3, -2, -1, 1, 3}},
                        {{-3, -1, 1, 2, 3}, {-3, 1, 3}, {-3, -2, -1, 2, 3}}));
    REQUIRE_NOTHROW(
        make<PBR>({{-3, -1}, {-3, -2, -1, 1, 2, 3}, {-3, -2, -1, 1, 3}},
                  {{-3, -1, 1, 2, 3}, {-3, 1, 3}, {-3, -2, -1, 2, 3}}));

    REQUIRE_NOTHROW(PBR({{}, {}}));

    REQUIRE_THROWS_AS(
        make<PBR>({{-4, -1}, {-3, -2, -1, 1, 2, 3}, {-3, -2, -1, 1, 3}},
                  {{-3, -1, 1, 2, 3}, {-3, 1, 3}, {-3, -2, -1, 2, 3}}),
        LibsemigroupsException);

    REQUIRE_THROWS_AS(
        make<PBR>({{-4, -1}, {-3, -2, -1, 1, 2, 3}, {-3, -2, -1, 1, 3}},
                  {{-3, -1, 1, 2, 3}, {-3, 1, 3}, {-3, -2, -1, 2, 3}}),
        LibsemigroupsException);

    REQUIRE_THROWS_AS(
        make<PBR>(
            {{-4, -1}, {-3, -2, -1, 1, 2, 3}, {-3, -2, -1, 1, 3}},
            {{-3, -1, 1, 2, 3}, {-3, 1, 3}, {-3, -2, -1, 2, 3}, {-1, -2}}),
        LibsemigroupsException);
    REQUIRE_THROWS_AS(
        make<PBR>({{-3, -1, 1, 2, 3}, {-3, 1, 3}, {-3, -2, -1, 2, 3}},
                  {{-4, -1}, {-3, -2, -1, 1, 2, 3}, {-3, -2, -1, 1, 3}}),
        LibsemigroupsException);

    REQUIRE_THROWS_AS(make<PBR>({{}, {2}, {1}, {3, 0}}),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("PBR", "006", "product exceptions", "[quick][pbr]") {
    PBR x(3);
    REQUIRE_THROWS(x.product_inplace(PBR(3), PBR(4)));
    REQUIRE_THROWS(x.product_inplace(PBR(4), PBR(4)));
    REQUIRE_THROWS(x.product_inplace(x, PBR(3)));
    REQUIRE_THROWS(
        x.product_inplace(PBR({{4}, {5}, {6}, {7}, {8}, {9}}), PBR(3)));
  }

  LIBSEMIGROUPS_TEST_CASE("PBR", "007", "operators", "[quick][pbr]") {
    PBR x({{3, 5},
           {0, 1, 2, 3, 4, 5},
           {0, 2, 3, 4, 5},
           {0, 1, 2, 3, 5},
           {0, 2, 5},
           {1, 2, 3, 4, 5}});
    REQUIRE(x[0] == std::vector<uint32_t>({3, 5}));
    REQUIRE(x[5] == std::vector<uint32_t>({1, 2, 3, 4, 5}));
    REQUIRE(x[0] == x.at(0));
    REQUIRE(x[5] == x.at(5));
    REQUIRE_THROWS(x.at(6));

    PBR y({{0, 3, 4, 5},
           {2, 4, 5},
           {1, 2, 5},
           {2, 3, 4, 5},
           {2, 3, 4, 5},
           {1, 2, 4}});

    PBR expected({{0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5},
                  {0, 1, 2, 3, 4, 5}});
    REQUIRE(x * y == expected);
    REQUIRE(y * y != expected);
    REQUIRE(y * y > expected);
    REQUIRE(y * y >= expected);
    REQUIRE(expected < y * y);
    REQUIRE(expected <= y * y);
    REQUIRE(x * x >= expected);
    REQUIRE(expected <= x * x);
  }

  LIBSEMIGROUPS_TEST_CASE("PBR", "008", "to_string", "[quick][pbr]") {
    PBR x({{3, 5},
           {0, 1, 2, 3, 4, 5},
           {0, 2, 3, 4, 5},
           {0, 1, 2, 3, 5},
           {0, 2, 5},
           {1, 2, 3, 4, 5}});
    REQUIRE_NOTHROW(detail::to_string(x));
    x = PBR({});
    REQUIRE_NOTHROW(detail::to_string(x));
    std::stringbuf buf;
    std::ostream   os(&buf);
    os << x;  // doesn't do anything visible
  }

  LIBSEMIGROUPS_TEST_CASE("PBR", "009", "one", "[quick][pbr]") {
    PBR x({{3, 5},
           {0, 1, 2, 3, 4, 5},
           {0, 2, 3, 4, 5},
           {0, 1, 2, 3, 5},
           {0, 2, 5},
           {1, 2, 3, 4, 5}});
    REQUIRE(x == x * pbr::one(x));
    REQUIRE(x == pbr::one(x) * x);
    REQUIRE(x == x * pbr::one(3));
    REQUIRE(x == pbr::one(3) * x);
  }

  LIBSEMIGROUPS_TEST_CASE("PBR", "010", "adapters", "[quick][pbr]") {
    PBR x({});
    REQUIRE_NOTHROW(IncreaseDegree<PBR>()(x, 0));
  }

  LIBSEMIGROUPS_TEST_CASE("PBR", "011", "make<PBR>", "[quick][pbr]") {
    PBR x({{-1, 1}, {2}}, {{-2, 1}, {-1, 2}});
    pbr::throw_if_invalid(x);
    REQUIRE(make<PBR>({}) == PBR({}));
    REQUIRE(make<PBR>({{1, 2}, {0, 3}, {2, 3}, {1}})
            == PBR({{1, 2}, {0, 3}, {2, 3}, {1}}));
    REQUIRE(make<PBR>({{-1, 1}, {2}}, {{-2, 1}, {-1, 2}})
            == PBR({{-1, 1}, {2}}, {{-2, 1}, {-1, 2}}));
  }
}  // namespace libsemigroups
