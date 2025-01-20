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

#include <algorithm>    // for lexicographical_compare
#include <array>        // for array
#include <cstddef>      // for size_t
#include <numeric>      // for accumulate
#include <type_traits>  // for move, swap, dec...
#include <utility>      // for operator==, pair
#include <vector>       // for vector

#define CATCH_CONFIG_ENABLE_PAIR_STRINGMAKER

#include "catch_amalgamated.hpp"  // for REQUIRE, REQUIRE_THROWS_AS, REQUIRE_NOTHROW
#include "test-main.hpp"          // for LIBSEMIGROUPS_TEST_CASE

#include "libsemigroups/adapters.hpp"       // for Complexity, Degree
#include "libsemigroups/bmat-fastest.hpp"   // for BMatFastest
#include "libsemigroups/bmat8.hpp"          // for BMat8
#include "libsemigroups/constants.hpp"      // for NEGATIVE_INFINITY
#include "libsemigroups/debug.hpp"          // for LIBSEMIGROUPS_A...
#include "libsemigroups/detail/report.hpp"  // for ReportGuard
#include "libsemigroups/exception.hpp"      // for LibsemigroupsEx...
#include "libsemigroups/matrix.hpp"         // for BMat, NTPMat

#include "libsemigroups/detail/containers.hpp"  // for StaticVector1

namespace libsemigroups {
  template <size_t N>
  class BitSet;

  constexpr bool REPORT = false;

  namespace {

    ////////////////////////////////////////////////////////////////////////
    // For comparison
    ////////////////////////////////////////////////////////////////////////

    template <typename Plus, typename Container>
    struct RowAddition {
      void operator()(Container& x, Container const& y) const {
        LIBSEMIGROUPS_ASSERT(x.size() == y.size());
        for (size_t i = 0; i < x.size(); ++i) {
          x[i] = Plus()(x[i], y[i]);
        }
      }

      void operator()(Container&       res,
                      Container const& x,
                      Container const& y) const {
        LIBSEMIGROUPS_ASSERT(res.size() == x.size());
        LIBSEMIGROUPS_ASSERT(x.size() == y.size());
        for (size_t i = 0; i < x.size(); ++i) {
          res[i] = Plus()(x[i], y[i]);
        }
      }
    };

    template <typename Prod, typename Container>
    Container scalar_row_product(Container                      row,
                                 typename Container::value_type scalar) {
      Container out(row);
      for (size_t i = 0; i < out.size(); ++i) {
        out[i] = Prod()(out[i], scalar);
      }
      return out;
    }

    template <size_t dim, size_t thresh>
    void tropical_max_plus_row_basis(std::vector<std::array<int, dim>>& rows) {
      static thread_local std::vector<std::array<int, dim>> buf;
      buf.clear();
      std::sort(rows.begin(), rows.end());
      for (size_t row = 0; row < rows.size(); ++row) {
        std::array<int, dim> sum;
        sum.fill(NEGATIVE_INFINITY);
        if (row == 0 || rows[row] != rows[row - 1]) {
          for (size_t row2 = 0; row2 < row; ++row2) {
            int max_scalar = thresh;
            for (size_t col = 0; col < dim; ++col) {
              if (rows[row2][col] == NEGATIVE_INFINITY) {
                continue;
              }
              if (rows[row][col] >= rows[row2][col]) {
                if (rows[row][col] != thresh) {
                  max_scalar
                      = std::min(max_scalar, rows[row][col] - rows[row2][col]);
                }
              } else {
                max_scalar = NEGATIVE_INFINITY;
                break;
              }
            }
            if (max_scalar != NEGATIVE_INFINITY) {
              auto scalar_prod
                  = scalar_row_product<MaxPlusTruncProd<thresh, int>,
                                       std::array<int, dim>>(rows[row2],
                                                             max_scalar);
              RowAddition<MaxPlusPlus<int>, std::array<int, dim>>()(
                  sum, scalar_prod);
            }
          }
          if (sum != rows[row]) {
            buf.push_back(rows[row]);
          }
        }
      }
      std::swap(buf, rows);
    }

  }  // namespace

  ////////////////////////////////////////////////////////////////////////
  // Test cases - BMat
  ////////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "000",
                                   "boolean matrix - test 1",
                                   "[quick]",
                                   BMat<2>,
                                   BMat<>) {
    auto rg = ReportGuard(REPORT);
    {
      TestType m = make<TestType>({{0, 1}, {0, 1}});
      REQUIRE_NOTHROW(matrix::throw_if_bad_entry(m));
      REQUIRE(m == TestType({{0, 1}, {0, 1}}));
      REQUIRE(!(m == TestType({{0, 0}, {0, 1}})));
      REQUIRE(m == TestType({{0, 1}, {0, 1}}));
      m.product_inplace_no_checks(TestType({{0, 0}, {0, 0}}),
                                  TestType({{0, 0}, {0, 0}}));
      REQUIRE(m == TestType({{0, 0}, {0, 0}}));
      m.product_inplace_no_checks(TestType({{0, 0}, {0, 0}}),
                                  TestType({{1, 1}, {1, 1}}));
      REQUIRE(m == TestType({{0, 0}, {0, 0}}));
      m.product_inplace_no_checks(TestType({{1, 1}, {1, 1}}),
                                  TestType({{0, 0}, {0, 0}}));
      REQUIRE(m == TestType({{0, 0}, {0, 0}}));

      m.product_inplace_no_checks(TestType({{0, 1}, {1, 0}}),
                                  TestType({{1, 0}, {1, 0}}));
      REQUIRE(m == TestType({{1, 0}, {1, 0}}));
      size_t const M = detail::BitSetCapacity<TestType>::value;
      detail::StaticVector1<BitSet<M>, M> result;
      matrix::bitset_rows(m, result);
      REQUIRE(result.size() == 2);
      REQUIRE(matrix::bitset_rows(m).size() == 2);
      result.clear();
      matrix::bitset_row_basis(m, result);
      REQUIRE(result.size() == 1);
      REQUIRE(matrix::bitset_row_basis(m).size() == 1);
      REQUIRE(std::vector<bool>(m.cbegin(), m.cend())
              == std::vector<bool>({true, false, true, false}));
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
      REQUIRE(std::vector<bool>(m.begin(), m.end())  // ****
              == std::vector<bool>({true, false, true, false}));
#pragma GCC diagnostic pop
    }

    {
      TestType m({{1, 1}, {0, 0}});
      using RowView = typename TestType::RowView;

      auto r = matrix::rows(m);
      REQUIRE(std::vector<bool>(r[0].cbegin(), r[0].cend())
              == std::vector<bool>({true, true}));
      REQUIRE(std::vector<bool>(r[1].cbegin(), r[1].cend())
              == std::vector<bool>({false, false}));
      REQUIRE(r.size() == 2);
      std::sort(r.begin(), r.end(), [](RowView const& rv1, RowView const& rv2) {
        return std::lexicographical_compare(
            rv1.begin(), rv1.end(), rv2.begin(), rv2.end());
      });
      REQUIRE(std::vector<bool>(r[0].cbegin(), r[0].cend())
              == std::vector<bool>({false, false}));
      REQUIRE(std::vector<bool>(r[1].cbegin(), r[1].cend())
              == std::vector<bool>({true, true}));
    }

    {
      using Row = typename TestType::Row;

      TestType A(2, 2);
      std::fill(A.begin(), A.end(), false);
      REQUIRE(A.number_of_rows() == 2);
      REQUIRE(A.number_of_cols() == 2);
      REQUIRE(A == TestType({{false, false}, {false, false}}));

      A(0, 0) = true;
      A(1, 1) = true;
      REQUIRE(A == TestType({{true, false}, {false, true}}));

      TestType B(2, 2);
      B(0, 1) = true;
      B(1, 0) = true;
      B(0, 0) = false;
      B(1, 1) = false;
      REQUIRE(B == TestType({{false, true}, {true, false}}));

      REQUIRE(A + B == TestType({{true, true}, {true, true}}));
      REQUIRE(A * B == B);
      REQUIRE(B * A == B);
      REQUIRE(B * B == A);
      REQUIRE((A + B) * B == TestType({{true, true}, {true, true}}));

      Row C({0, 1});
      REQUIRE(C.number_of_rows() == 1);
      REQUIRE(C.number_of_cols() == 2);

      auto rv = A.row(0);
      Row  D(rv);
      REQUIRE(D.number_of_rows() == 1);
      REQUIRE(D.number_of_cols() == 2);
      REQUIRE(D != C);
      auto views = matrix::rows(A);
      REQUIRE(B < A);
      B.swap(A);
      REQUIRE(A < B);
      std::swap(B, A);
      REQUIRE(B < A);
      REQUIRE(views[0] == Row({true, false}));
      REQUIRE(Row({true, false}) == views[0]);
      REQUIRE(Row({true, true}) != views[0]);
      REQUIRE(Row({false, false}) < views[0]);
      REQUIRE(A.hash_value() != 0);
      A *= false;
      REQUIRE(A == TestType({{false, false}, {false, false}}));
      auto r = Row({true, false});
      views  = matrix::rows(B);
      REQUIRE(views[0].size() == 2);
      r += views[0];
      REQUIRE(r.number_of_cols() == 2);
      REQUIRE(r.number_of_rows() == 1);
      REQUIRE(r == Row({true, true}));

      auto E = TestType::one(2);
      REQUIRE(E.number_of_rows() == 2);
      REQUIRE(E.number_of_cols() == 2);
      auto viewse = matrix::rows(E);
      REQUIRE(viewse.size() == 2);

      std::ostringstream oss;
      oss << E;  // Does not do anything visible

      std::stringbuf buff;
      std::ostream   os(&buff);
      os << E;  // Also does not do anything visible
    }
    {
      TestType m({{0, 0}, {0, 0}});
      using scalar_type = typename TestType::scalar_type;
      auto it           = m.cbegin();
      REQUIRE(m.coords(it) == std::pair<scalar_type, scalar_type>({0, 0}));
      REQUIRE(m.coords(++it) == std::pair<scalar_type, scalar_type>({0, 1}));
      REQUIRE(m.coords(++it) == std::pair<scalar_type, scalar_type>({1, 0}));
      REQUIRE(m.coords(++it) == std::pair<scalar_type, scalar_type>({1, 1}));
    }
    {
      REQUIRE_THROWS_AS(make<TestType>({{0, 0}, {0, 2}}),
                        LibsemigroupsException);
    }
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "002",
                                   "boolean matrix - test 2",
                                   "[quick]",
                                   BMat<3>,
                                   BMat<>) {
    auto     rg = ReportGuard(REPORT);
    TestType m(3, 3);
    m.product_inplace_no_checks(TestType({{1, 1, 0}, {0, 0, 1}, {1, 0, 1}}),
                                TestType({{1, 0, 1}, {0, 0, 1}, {1, 1, 0}}));
    REQUIRE(m == TestType({{1, 0, 1}, {1, 1, 0}, {1, 1, 1}}));
    TestType* A = new TestType();
    delete A;
    typename TestType::Row* B = new typename TestType::Row();
    delete B;
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "003",
                                   "boolean matrix - test 3",
                                   "[quick]",
                                   BMat<2>,
                                   BMat<>) {
    TestType AB(2, 2), A(2, 2), B(2, 2);
    std::fill(A.begin(), A.end(), false);
    std::fill(B.begin(), B.end(), false);
    std::fill(AB.begin(), AB.end(), false);
    A(1, 1) = true;

    AB.product_inplace_no_checks(A, B);
    REQUIRE(AB == B);
    REQUIRE(A.one() == TestType({{true, false}, {false, true}}));
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "004",
                                   "boolean matrix - test 4",
                                   "[quick]",
                                   BMat<3>,
                                   BMat<>) {
    auto x = TestType({{1, 0, 1}, {0, 1, 0}, {0, 1, 0}});
    auto y = TestType({{0, 0, 0}, {0, 0, 0}, {0, 0, 0}});
    auto z = TestType({{0, 0, 0}, {0, 0, 0}, {0, 0, 0}});
    REQUIRE(y == z);
    z.product_inplace_no_checks(x, y);
    REQUIRE(y == z);
    z.product_inplace_no_checks(y, x);
    REQUIRE(y == z);
    REQUIRE(!(y < z));
    auto id = x.one();
    z.product_inplace_no_checks(id, x);
    REQUIRE(z == x);
    z.product_inplace_no_checks(x, id);
    REQUIRE(z == x);
    REQUIRE(x.hash_value() != 0);
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "006",
                                   "boolean matrix - row_basis",
                                   "[quick]",
                                   BMat<3>,
                                   BMat<>) {
    using RowView = typename TestType::RowView;
    auto x        = make<TestType>({{1, 0, 0}, {1, 0, 0}, {1, 0, 0}});
    REQUIRE(matrix::row_basis(x).size() == 1);
    REQUIRE(matrix::row_space_size(x) == 1);
    x = make<TestType>({{1, 0, 0}, {1, 1, 0}, {1, 1, 1}});
    REQUIRE(matrix::row_basis(x).size() == 3);
    REQUIRE_THROWS_AS(x.row(3), LibsemigroupsException);
    std::vector<RowView> v = {x.row(0), x.row(2)};
    REQUIRE(matrix::row_basis<TestType>(v).size() == 2);
    REQUIRE(matrix::row_space_size(x) == 3);
    x = make<TestType>({{1, 0, 0}, {0, 1, 1}, {1, 1, 1}});
    REQUIRE(matrix::row_basis(x).size() == 2);
    REQUIRE(matrix::row_space_size(x) == 3);
    x = make<TestType>({{1, 0, 0}, {0, 0, 1}, {0, 1, 0}});
    REQUIRE(matrix::row_space_size(x) == 7);
    std::vector<typename TestType::RowView> views;
    std::vector<typename TestType::RowView> result;
    matrix::row_basis<TestType, std::vector<typename TestType::RowView>&>(
        views, result);
  }

  ////////////////////////////////////////////////////////////////////////
  // Test cases - IntMat
  ////////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "008",
                                   "integer matrix - test 1",
                                   "[quick]",
                                   IntMat<3>,
                                   IntMat<>) {
    {
      auto x        = TestType({{-2, 2, 0}, {-1, 0, 0}, {1, -3, 1}});
      auto expected = TestType({{-2, 2, 0}, {-1, 0, 0}, {1, -3, 1}});
      REQUIRE(x == expected);

      auto y = TestType({{-100, 0, 0}, {0, 1, 0}, {1, -1, 0}});
      REQUIRE(!(x == y));

      y.product_inplace_no_checks(x, x);
      expected = TestType({{2, -4, 0}, {2, -2, 0}, {2, -1, 1}});
      REQUIRE(y == expected);
      REQUIRE(y.number_of_rows() == 3);

      REQUIRE(x < y);
      REQUIRE(Degree<TestType>()(x) == 3);
      REQUIRE(Degree<TestType>()(y) == 3);
      REQUIRE(Complexity<TestType>()(x) == 27);
      REQUIRE(Complexity<TestType>()(y) == 27);
      auto id = x.one();
      y.product_inplace_no_checks(id, x);
      REQUIRE(y == x);
      y.product_inplace_no_checks(x, id);
      REQUIRE(y == x);
    }
    {
      auto x        = TestType({{-2, 2, 0}, {-1, 0, 0}, {1, -3, 1}});
      auto expected = TestType({{-2, 2, 0}, {-1, 0, 0}, {1, -3, 1}});
      REQUIRE(x == expected);

      auto y = TestType({{-100, 0, 0}, {0, 1, 0}, {1, -1, 0}});
      REQUIRE(!(x == y));

      y.product_inplace_no_checks(x, x);
      expected = TestType({{2, -4, 0}, {2, -2, 0}, {2, -1, 1}});
      REQUIRE(y == expected);

      REQUIRE(x < y);
      auto id = x.one();
      y.product_inplace_no_checks(id, x);
      REQUIRE(y == x);
      y.product_inplace_no_checks(x, id);
      REQUIRE(y == x);
      REQUIRE(Hash<TestType>()(y) != 0);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Matrix",
                          "010",
                          "integer matrix - code cov",
                          "[quick]") {
    IntMat<>* A = new IntMat<>();
    delete A;
    IntMat<3>* B = new IntMat<3>();
    delete B;
    IntMat<3>::Row* C = new IntMat<3>::Row();
    delete C;
    IntMat<2>::Row* D = new IntMat<2>::Row();
    delete D;
    IntMat<2>* E = new IntMat<2>();
    delete E;
    IntMat<> F(3, 3);
    IntMat<> G(4, 4);
    std::swap(F, G);
    REQUIRE(G.number_of_cols() == 3);
    REQUIRE(G.number_of_rows() == 3);
    REQUIRE(F.number_of_cols() == 4);
    REQUIRE(F.number_of_rows() == 4);
  }

  ////////////////////////////////////////////////////////////////////////
  // Test cases - MaxPlusMat
  ////////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "011",
                                   "max-plus matrix - test 1",
                                   "[quick]",
                                   MaxPlusMat<>,
                                   MaxPlusMat<3>) {
    auto x        = TestType({{-2, 2, 0}, {-1, 0, 0}, {1, -3, 1}});
    auto expected = TestType({{-2, 2, 0}, {-1, 0, 0}, {1, -3, 1}});
    REQUIRE(x == expected);

    auto y = TestType{{-100, 0, 0}, {0, 1, 0}, {1, -1, 0}};
    REQUIRE(!(x == y));
    REQUIRE(x != y);

    y.product_inplace_no_checks(x, x);
    expected = TestType({{1, 2, 2}, {1, 1, 1}, {2, 3, 2}});
    REQUIRE(y == expected);

    REQUIRE(x < y);
    REQUIRE(Degree<TestType>()(x) == 3);
    REQUIRE(Degree<TestType>()(y) == 3);
    REQUIRE(Complexity<TestType>()(x) == 27);
    REQUIRE(Complexity<TestType>()(y) == 27);
    auto id = x.one();
    y.product_inplace_no_checks(id, x);
    REQUIRE(y == x);
    y.product_inplace_no_checks(x, id);
    REQUIRE(y == x);
    REQUIRE(Hash<TestType>()(y) != 0);
  }

  LIBSEMIGROUPS_TEST_CASE("Matrix", "013", "MaxPlusMat code cov", "[quick]") {
    MaxPlusMat<3>* B = new MaxPlusMat<3>();
    delete B;
    MaxPlusMat<3>::Row* C = new MaxPlusMat<3>::Row();
    delete C;
  }

  ////////////////////////////////////////////////////////////////////////
  // Test cases - MinPlusMat
  ////////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "014",
                                   "min-plus matrix - test 1",
                                   "[quick]",
                                   MinPlusMat<3>,
                                   MinPlusMat<>) {
    {
      auto x        = TestType({{-2, 2, 0}, {-1, 0, 0}, {1, -3, 1}});
      auto expected = TestType({{-2, 2, 0}, {-1, 0, 0}, {1, -3, 1}});
      // Just testing the below doesn't compile
      // matrix::row_basis(x);
      REQUIRE(x == expected);

      auto y = TestType({{-100, 0, 0}, {0, 1, 0}, {1, -1, 0}});
      REQUIRE(!(x == y));

      y.product_inplace_no_checks(x, x);
      expected = TestType({{-4, -3, -2}, {-3, -3, -1}, {-4, -3, -3}});
      REQUIRE(y == expected);

      REQUIRE(!(x < y));
      REQUIRE(Degree<TestType>()(x) == 3);
      REQUIRE(Degree<TestType>()(y) == 3);
      REQUIRE(Complexity<TestType>()(x) == 27);
      REQUIRE(Complexity<TestType>()(y) == 27);
      auto id = x.one();
      y.product_inplace_no_checks(id, x);
      REQUIRE(y == x);
      y.product_inplace_no_checks(x, id);
      REQUIRE(y == x);
    }
    {
      auto x        = TestType({{22, 21, 0}, {10, 0, 0}, {1, 32, 1}});
      auto expected = TestType({{22, 21, 0}, {10, 0, 0}, {1, 32, 1}});
      REQUIRE(x == expected);

      auto y = TestType({{10, 0, 0}, {0, 1, 0}, {1, 1, 0}});
      REQUIRE(!(x == y));

      y.product_inplace_no_checks(x, x);
      REQUIRE(y == TestType({{1, 21, 1}, {1, 0, 0}, {2, 22, 1}}));

      REQUIRE(x > y);
      REQUIRE(Degree<TestType>()(x) == 3);
      REQUIRE(Degree<TestType>()(y) == 3);
      REQUIRE(Complexity<TestType>()(x) == 27);
      REQUIRE(Complexity<TestType>()(y) == 27);
      auto id = x.one();
      y.product_inplace_no_checks(id, x);
      REQUIRE(y == x);
      y.product_inplace_no_checks(x, id);
      REQUIRE(y == x);
      REQUIRE(Hash<TestType>()(y) != 0);
    }
  }

  LIBSEMIGROUPS_TEST_CASE("Matrix", "016", "MinPlusMat code cov", "[quick]") {
    MinPlusMat<3>* B = new MinPlusMat<3>();
    delete B;
    MinPlusMat<3>::Row* C = new MinPlusMat<3>::Row();
    delete C;
  }

  ////////////////////////////////////////////////////////////////////////
  // Test cases - MaxPlusTruncMat
  ////////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "017",
                                   "max-plus trunc. matrix - test 1",
                                   "[quick]",
                                   (MaxPlusTruncMat<5, 2>),
                                   MaxPlusTruncMat<5>,
                                   MaxPlusTruncMat<>) {
    REQUIRE_THROWS_AS(MaxPlusTruncSemiring<>(-1), LibsemigroupsException);
    // Threshold 5, 2 x 2
    MaxPlusTruncSemiring<>* sr = nullptr;
    if constexpr (!std::is_same_v<typename TestType::semiring_type, void>) {
      sr = new MaxPlusTruncSemiring(5);
    }
    using scalar_type = typename TestType::scalar_type;
    {
      TestType m1(sr, 2, 2);
      std::fill(m1.begin(), m1.end(), NEGATIVE_INFINITY);
      REQUIRE(m1
              == make<TestType>(sr,
                                {{NEGATIVE_INFINITY, NEGATIVE_INFINITY},
                                 {NEGATIVE_INFINITY, NEGATIVE_INFINITY}}));
      TestType m2(sr, 2, 2);
      std::fill(m2.begin(), m2.end(), 4);
      REQUIRE(m1 + m2 == m2);
      REQUIRE(m2(0, 1) == 4);
    }

    auto rg = ReportGuard(REPORT);
    {
      std::vector<std::array<scalar_type, 2>> expected;
      expected.push_back({1, 1});
      expected.push_back({0, 0});
      tropical_max_plus_row_basis<2, 5>(expected);
      REQUIRE(expected.size() == 1);
      REQUIRE(expected.at(0) == std::array<scalar_type, 2>({0, 0}));

      TestType m(sr, {{1, 1}, {0, 0}});
      auto     r = matrix::row_basis(m);
      REQUIRE(r.size() == 1);
      REQUIRE(std::vector<scalar_type>(r[0].cbegin(), r[0].cend())
              == std::vector<scalar_type>({0, 0}));
    }
    {
      TestType m(sr, {{1, 1}, {0, 0}});
      m      = m.one();
      auto r = matrix::row_basis(m);
      REQUIRE(r.size() == 2);
      REQUIRE(std::vector<scalar_type>(r[0].cbegin(), r[0].cend())
              == std::vector<scalar_type>({NEGATIVE_INFINITY, 0}));
      REQUIRE(std::vector<scalar_type>(r[1].cbegin(), r[1].cend())
              == std::vector<scalar_type>({0, NEGATIVE_INFINITY}));
    }
    std::vector<typename TestType::RowView> views;
    std::vector<typename TestType::RowView> result;
    matrix::row_basis<TestType>(views, result);
    delete sr;
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "020",
                                   "max-plus trunc. matrix - test 2",
                                   "[quick]",
                                   (MaxPlusTruncMat<5, 4>),
                                   MaxPlusTruncMat<5>,
                                   MaxPlusTruncMat<>) {
    // Threshold 5, 4 x 4
    using scalar_type = typename TestType::scalar_type;
    using Row         = typename TestType::Row;

    MaxPlusTruncSemiring<>* sr = nullptr;
    if constexpr (!std::is_same_v<typename TestType::semiring_type, void>) {
      sr = new MaxPlusTruncSemiring(5);
    }

    auto m  = make<TestType>(sr,
                            {{2, 2, 0, 1},
                              {0, 0, 1, 3},
                              {1, NEGATIVE_INFINITY, 0, 0},
                              {0, 1, 0, 1}});
    auto rg = ReportGuard(REPORT);
    auto r  = matrix::row_basis(m);
    REQUIRE(r.size() == 4);
    REQUIRE(r[0] == make<Row>(sr, {0, 0, 1, 3}));
    REQUIRE(r[1] == make<Row>(sr, {0, 1, 0, 1}));
    REQUIRE(r[2] == make<Row>(sr, {1, NEGATIVE_INFINITY, 0, 0}));
    REQUIRE(r[3] == make<Row>(sr, {2, 2, 0, 1}));
    m.transpose();
    REQUIRE(m
            == make<TestType>(sr,
                              {{2, 0, 1, 0},
                               {2, 0, NEGATIVE_INFINITY, 1},
                               {0, 1, 0, 0},
                               {1, 3, 0, 1}}));
    m.transpose();
    REQUIRE(m
            == make<TestType>(sr,
                              {{2, 2, 0, 1},
                               {0, 0, 1, 3},
                               {1, NEGATIVE_INFINITY, 0, 0},
                               {0, 1, 0, 1}}));

    std::vector<std::array<scalar_type, 4>> expected;
    expected.push_back({2, 2, 0, 1});
    expected.push_back({0, 0, 1, 3});
    expected.push_back({1, NEGATIVE_INFINITY, 0, 0});
    expected.push_back({0, 1, 0, 1});
    tropical_max_plus_row_basis<4, 5>(expected);
    REQUIRE(expected.size() == 4);
    REQUIRE(m * TestType::one(sr, 4) == m);
    REQUIRE(TestType::one(sr, 4) * m == m);
    delete sr;
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "023",
                                   "max-plus trunc. matrix - test 3",
                                   "[quick]",
                                   (MaxPlusTruncMat<33, 3>),
                                   MaxPlusTruncMat<33>,
                                   MaxPlusTruncMat<>) {
    // Threshold 33, 3 x 3
    MaxPlusTruncSemiring<>* sr = nullptr;
    if constexpr (!std::is_same_v<typename TestType::semiring_type, void>) {
      sr = new MaxPlusTruncSemiring(33);
    }
    auto x        = TestType(sr, {{22, 21, 0}, {10, 0, 0}, {1, 32, 1}});
    auto expected = TestType(sr, {{22, 21, 0}, {10, 0, 0}, {1, 32, 1}});
    REQUIRE(x == expected);

    REQUIRE_THROWS_AS(make<TestType>(sr, {{-100, 0, 0}, {0, 1, 0}, {1, -1, 0}}),
                      LibsemigroupsException);
    auto y = TestType(sr, {{10, 0, 0}, {0, 1, 0}, {1, 1, 0}});
    REQUIRE(!(x == y));

    y.product_inplace_no_checks(x, x);
    expected = TestType(sr, {{33, 33, 22}, {32, 32, 10}, {33, 33, 32}});
    REQUIRE(y == expected);

    REQUIRE(x < y);
    auto id = x.one();
    y.product_inplace_no_checks(id, x);
    REQUIRE(y == x);
    y.product_inplace_no_checks(x, id);
    REQUIRE(y == x);
    REQUIRE(Hash<TestType>()(y) != 0);
    REQUIRE(x * TestType::one(sr, 3) == x);
    REQUIRE(TestType::one(sr, 3) * x == x);
    delete sr;
  }

  LIBSEMIGROUPS_TEST_CASE("Matrix", "026", "MaxPlusMat code cov", "[quick]") {
    MaxPlusTruncMat<33, 3>* B = new MaxPlusTruncMat<33, 3>();
    delete B;
    MaxPlusTruncMat<5, 2>::Row* C = new MaxPlusTruncMat<5, 2>::Row();
    delete C;
    MaxPlusTruncMat<5, 4>::Row* D = new MaxPlusTruncMat<5, 4>::Row();
    delete D;
    auto x = MaxPlusTruncMat<5, 2>();
    REQUIRE(x.number_of_rows() == 2);
    MaxPlusTruncMat<5, 2>* E = new MaxPlusTruncMat<5, 2>();
    delete E;
    MaxPlusTruncMat<5, 4>* F = new MaxPlusTruncMat<5, 4>();
    delete F;
  }

  ////////////////////////////////////////////////////////////////////////
  // Test cases - MinPlusTruncMat
  ////////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "027",
                                   "min-plus trunc. matrix - test 1",
                                   "[quick]",
                                   (MinPlusTruncMat<33, 3>),
                                   MinPlusTruncMat<33>,
                                   MinPlusTruncMat<>) {
    // Threshold 33, 3 x 3
    MinPlusTruncSemiring<>* sr = nullptr;
    if constexpr (!std::is_same_v<typename TestType::semiring_type, void>) {
      sr = new MinPlusTruncSemiring(33);
    }
    auto x = TestType(sr, {{22, 21, 0}, {10, 0, 0}, {1, 32, 1}});

    auto expected = make<TestType>(sr, {{22, 21, 0}, {10, 0, 0}, {1, 32, 1}});
    REQUIRE(x == expected);

    auto y = TestType(sr, {{10, 0, 0}, {0, 1, 0}, {1, 1, 0}});
    REQUIRE(!(x == y));

    y.product_inplace_no_checks(x, x);
    expected = TestType(sr, {{1, 21, 1}, {1, 0, 0}, {2, 22, 1}});
    REQUIRE(y == expected);

    REQUIRE(!(x < y));
    REQUIRE(Degree<TestType>()(x) == 3);
    REQUIRE(Degree<TestType>()(y) == 3);
    REQUIRE(Complexity<TestType>()(x) == 27);
    REQUIRE(Complexity<TestType>()(y) == 27);
    auto id = x.one();
    y.product_inplace_no_checks(id, x);
    REQUIRE(y == x);
    y.product_inplace_no_checks(x, id);
    REQUIRE(y == x);
    REQUIRE(Hash<TestType>()(y) != 0);
    REQUIRE_THROWS_AS(
        make<TestType>(sr, {{-22, 21, 0}, {10, 0, 0}, {1, 32, 1}}),
        LibsemigroupsException);
    REQUIRE(x * TestType::one(sr, 3) == x);
    REQUIRE(TestType::one(sr, 3) * x == x);
    delete sr;
  }
  ////////////////////////////////////////////////////////////////////////
  // Test cases - NTPMat
  ////////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "030",
                                   "3x3 matrix, t = 0, p = 3",
                                   "[quick]",
                                   (NTPMat<0, 3, 3, 3>),
                                   (NTPMat<0, 3, 3>),
                                   NTPMat<>) {
    REQUIRE_THROWS_AS(NTPSemiring<int>(4, -1), LibsemigroupsException);
    REQUIRE_THROWS_AS(NTPSemiring<int>(-1, -1), LibsemigroupsException);

    NTPSemiring<>* sr = nullptr;
    if constexpr (!std::is_same_v<typename TestType::semiring_type, void>) {
      sr = new NTPSemiring<>(0, 3);
    }
    using Row   = typename TestType::Row;
    auto     rg = ReportGuard(REPORT);
    TestType m(sr, 3, 3);
    // REQUIRE(matrix::throw_if_bad_entry(m)); // m might not be valid!
    m.product_inplace_no_checks(
        make<TestType>(sr, {{1, 1, 0}, {0, 0, 1}, {1, 0, 1}}),
        make<TestType>(sr, {{1, 0, 1}, {0, 0, 1}, {1, 1, 0}}));
    REQUIRE(m == make<TestType>(sr, {{1, 0, 2}, {1, 1, 0}, {2, 1, 1}}));
    REQUIRE(m.row(0) == Row(sr, {1, 0, 2}));
    REQUIRE(m.row(0).size() == 3);
    auto r = matrix::rows(m);
    REQUIRE(r[0] == Row(sr, {1, 0, 2}));
    REQUIRE(r[1] == Row(sr, {1, 1, 0}));
    REQUIRE(r[2] == Row(sr, {2, 1, 1}));
    REQUIRE(m * TestType::one(sr, 3) == m);
    REQUIRE(TestType::one(sr, 3) * m == m);
    delete sr;
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "034",
                                   "4x4 matrix, t = 0, p = 10",
                                   "[quick]",
                                   (NTPMat<0, 10>),
                                   NTPMat<>) {
    NTPSemiring<>* sr = nullptr;
    if constexpr (!std::is_same_v<typename TestType::semiring_type, void>) {
      sr = new NTPSemiring<>(0, 10);
    }
    using Row         = typename TestType::Row;
    using RowView     = typename TestType::RowView;
    using scalar_type = typename TestType::scalar_type;

    auto rg = ReportGuard(REPORT);

    TestType m = make<TestType>(
        sr, {{1, 1, 0, 0}, {2, 0, 2, 0}, {1, 2, 3, 9}, {0, 0, 0, 7}});
    REQUIRE(m.number_of_cols() == 4);
    REQUIRE(m.number_of_rows() == 4);
    auto r = matrix::rows(m);
    REQUIRE(r.size() == 4);
    REQUIRE(std::vector<scalar_type>(r[0].cbegin(), r[0].cend())
            == std::vector<scalar_type>({1, 1, 0, 0}));
    r[0] += r[1];
    REQUIRE(std::vector<scalar_type>(r[0].cbegin(), r[0].cend())
            == std::vector<scalar_type>({3, 1, 2, 0}));
    REQUIRE(std::vector<scalar_type>(r[1].cbegin(), r[1].cend())
            == std::vector<scalar_type>({2, 0, 2, 0}));
    REQUIRE(m
            == make<TestType>(
                sr, {{3, 1, 2, 0}, {2, 0, 2, 0}, {1, 2, 3, 9}, {0, 0, 0, 7}}));
    REQUIRE(r[0][0] == 3);
    REQUIRE(r[0](0) == 3);
    REQUIRE(r[2](3) == 9);
    std::sort(r[0].begin(), r[0].end());
    REQUIRE(std::vector<scalar_type>(r[0].cbegin(), r[0].cend())
            == std::vector<scalar_type>({0, 1, 2, 3}));
    REQUIRE(m
            == make<TestType>(
                sr, {{0, 1, 2, 3}, {2, 0, 2, 0}, {1, 2, 3, 9}, {0, 0, 0, 7}}));
    r[0] += 9;
    REQUIRE(std::vector<scalar_type>(r[0].cbegin(), r[0].cend())
            == std::vector<scalar_type>({9, 0, 1, 2}));
    REQUIRE(m
            == make<TestType>(
                sr, {{9, 0, 1, 2}, {2, 0, 2, 0}, {1, 2, 3, 9}, {0, 0, 0, 7}}));
    r[1] *= 3;
    REQUIRE(m
            == make<TestType>(
                sr, {{9, 0, 1, 2}, {6, 0, 6, 0}, {1, 2, 3, 9}, {0, 0, 0, 7}}));
    REQUIRE(std::vector<scalar_type>(r[1].cbegin(), r[1].cend())
            == std::vector<scalar_type>({6, 0, 6, 0}));
    REQUIRE(r[2] < r[1]);
    r[1] = r[2];
    REQUIRE(m
            == make<TestType>(
                sr, {{9, 0, 1, 2}, {6, 0, 6, 0}, {1, 2, 3, 9}, {0, 0, 0, 7}}));
    REQUIRE(r[1] == r[2]);
    REQUIRE(r[1] == make<Row>(sr, {{1, 2, 3, 9}}));

    RowView rv;
    {
      rv = r[0];
      REQUIRE(rv == r[0]);
      REQUIRE(&rv != &r[0]);
    }
    delete sr;
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "036",
                                   "4x4 matrix, t = 0, p = 10",
                                   "[quick]",
                                   (NTPMat<0, 10, 4, 4>),
                                   (NTPMat<0, 10, 4>),
                                   (NTPMat<0, 10>),
                                   NTPMat<>) {
    using Row         = typename TestType::Row;
    NTPSemiring<>* sr = nullptr;
    if constexpr (!std::is_same_v<typename TestType::semiring_type, void>) {
      sr = new NTPSemiring<>(0, 10);
    }

    auto     rg = ReportGuard(REPORT);
    TestType m(sr, {{1, 1, 0, 0}, {2, 0, 2, 0}, {1, 2, 3, 9}, {0, 0, 0, 7}});
    REQUIRE(m.number_of_cols() == 4);
    REQUIRE(m.number_of_rows() == 4);
    auto r = matrix::rows(m);
    REQUIRE(r.size() == 4);
    REQUIRE(r[0] == make<Row>(sr, {{1, 1, 0, 0}}));
    REQUIRE(r[1] == make<Row>(sr, {{2, 0, 2, 0}}));
    REQUIRE(r[0] != make<Row>(sr, {{2, 0, 2, 0}}));
    REQUIRE(r[1] != make<Row>(sr, {{1, 1, 0, 0}}));
    REQUIRE(make<Row>(sr, {{1, 1, 0, 0}}) == r[0]);
    REQUIRE(make<Row>(sr, {{2, 0, 2, 0}}) == r[1]);
    REQUIRE(make<Row>(sr, {{2, 0, 2, 0}}) != r[0]);
    REQUIRE(make<Row>(sr, {{1, 1, 0, 0}}) != r[1]);
    REQUIRE(make<Row>(sr, {{1, 1, 0, 0}}) < Row(sr, {{9, 9, 9, 9}}));
    REQUIRE(r[0] < make<Row>(sr, {{9, 9, 9, 9}}));
    REQUIRE(!(make<Row>(sr, {{9, 9, 9, 9}}) < r[0]));
    Row x(r[3]);
    x *= 3;
    REQUIRE(x == make<Row>(sr, {{0, 0, 0, 1}}));
    REQUIRE(x.number_of_rows() == 1);
    REQUIRE(x.number_of_cols() == 4);
    REQUIRE(r[3] == make<Row>(sr, {{0, 0, 0, 7}}));
    REQUIRE(r[3] != x);
    REQUIRE(x != r[3]);
    REQUIRE(!(x != x));
    delete sr;
  }

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "038",
                                   "3x3 matrix, t = 33, p = 2",
                                   "[quick]",
                                   (NTPMat<33, 2>),
                                   NTPMat<>) {
    NTPSemiring<>* sr = nullptr;
    if constexpr (!std::is_same_v<typename TestType::semiring_type, void>) {
      sr = new NTPSemiring<>(33, 2);
    }
    auto x        = make<TestType>(sr, {{22, 21, 0}, {10, 0, 0}, {1, 32, 1}});
    auto expected = make<TestType>(sr, {{22, 21, 0}, {10, 0, 0}, {1, 32, 1}});
    REQUIRE(x == expected);
    REQUIRE(x.number_of_cols() == 3);
    REQUIRE(x.number_of_rows() == 3);

    auto y = make<TestType>(sr, {{10, 0, 0}, {0, 1, 0}, {1, 1, 0}});
    REQUIRE(!(x == y));

    y.product_inplace_no_checks(x, x);
    expected = make<TestType>(sr, {{34, 34, 0}, {34, 34, 0}, {33, 33, 1}});
    REQUIRE(y == expected);

    REQUIRE(x < y);
    auto id = x.one();
    y.product_inplace_no_checks(id, x);
    REQUIRE(y == x);
    y.product_inplace_no_checks(x, id);
    REQUIRE(y == x);
    REQUIRE(Hash<TestType>()(y) != 0);
    delete sr;
  }

  ////////////////////////////////////////////////////////////////////////
  // Test cases - ProjMaxPlusMat
  ////////////////////////////////////////////////////////////////////////

  LIBSEMIGROUPS_TEMPLATE_TEST_CASE("Matrix",
                                   "042",
                                   "3x3 matrix",
                                   "[quick]",
                                   (ProjMaxPlusMat<3, 3>),
                                   ProjMaxPlusMat<3>,
                                   ProjMaxPlusMat<>) {
    using Row     = typename TestType::Row;
    auto x        = make<TestType>({{-2, 2, 0}, {-1, 0, 0}, {1, -3, 1}});
    auto expected = make<TestType>({{-4, 0, -2}, {-3, -2, -2}, {-1, -5, -1}});
    REQUIRE(x == expected);
    REQUIRE(x.scalar_zero() == NEGATIVE_INFINITY);
    REQUIRE(x.scalar_one() == 0);

    auto y = make<TestType>({{NEGATIVE_INFINITY, 0, 0}, {0, 1, 0}, {1, -1, 0}});
    expected = make<TestType>(
        {{NEGATIVE_INFINITY, -1, -1}, {-1, 0, -1}, {0, -2, -1}});
    REQUIRE(y == expected);
    REQUIRE(!(x == y));

    y.product_inplace_no_checks(x, x);
    expected = TestType({{-2, -1, -1}, {-2, -2, -2}, {-1, 0, -1}});
    REQUIRE(y == expected);

    REQUIRE(x < y);
    REQUIRE(y > x);
    REQUIRE(Degree<TestType>()(x) == 3);
    REQUIRE(Degree<TestType>()(y) == 3);
    REQUIRE(Complexity<TestType>()(x) == 27);
    REQUIRE(Complexity<TestType>()(y) == 27);
    auto id = x.one();
    y.product_inplace_no_checks(id, x);
    REQUIRE(y == x);
    y.product_inplace_no_checks(x, id);
    REQUIRE(y == x);

    REQUIRE(make<TestType>({{-2, 2, 0}, {-1, 0, 0}, {1, -3, 1}}).hash_value()
            != 0);

    y = x;
    REQUIRE(&x != &y);
    REQUIRE(x == y);
    // TODO uncomment or remove
    // REQUIRE(y
    //         == make<TestType>(nullptr, {{-2, 2, 0}, {-1, 0, 0}, {1, -3,
    //         1}}));

    auto yy(y);
    REQUIRE(yy == y);

    std::ostringstream oss;
    oss << y;  // Does not do anything visible
    std::stringbuf buff;
    std::ostream   os(&buff);
    os << y;  // Also does not do anything visible

    REQUIRE(y.row(0) == make<Row>({-4, 0, -2}));
    REQUIRE(y.row(1) == Row({-3, -2, -2}));
    REQUIRE(Row(y.row(0)) == y.row(0));

    auto zz(std::move(y));

    TestType tt;
    REQUIRE(tt != zz);
    REQUIRE(TestType::one(3)
            == TestType({{0, NEGATIVE_INFINITY, NEGATIVE_INFINITY},
                         {NEGATIVE_INFINITY, 0, NEGATIVE_INFINITY},
                         {NEGATIVE_INFINITY, NEGATIVE_INFINITY, 0}}));
    REQUIRE(zz(0, 0) == -4);
    REQUIRE(zz.number_of_cols() == 3);
    zz += zz;
    REQUIRE(zz == TestType({{-2, 2, 0}, {-1, 0, 0}, {1, -3, 1}}));
    zz *= 2;
    REQUIRE(zz == TestType({{-2, 2, 0}, {-1, 0, 0}, {1, -3, 1}}));
    REQUIRE(zz + x == TestType({{-2, 2, 0}, {-1, 0, 0}, {1, -3, 1}}));
    REQUIRE(zz * x == TestType({{-2, -1, -1}, {-2, -2, -2}, {-1, 0, -1}}));
    REQUIRE(std::accumulate(zz.cbegin(), zz.cend(), 0) == -20);
    REQUIRE(std::accumulate(zz.begin(), zz.end(), 0) == -20);
    x.transpose();
    REQUIRE(x == TestType({{-4, -3, -1}, {0, -2, -5}, {-2, -2, -1}}));
    x.swap(zz);
    REQUIRE(zz == TestType({{-4, -3, -1}, {0, -2, -5}, {-2, -2, -1}}));
    REQUIRE(matrix::pow(x, 100)
            == TestType({{-1, 0, -1}, {-2, -1, -2}, {-1, 0, -1}}));
    REQUIRE_THROWS_AS(matrix::pow(x, -100), LibsemigroupsException);
    REQUIRE(matrix::pow(x, 1)
            == TestType({{-4, 0, -2}, {-3, -2, -2}, {-1, -5, -1}}));
    REQUIRE(matrix::pow(x, 0) == TestType::one(3));
  }

  LIBSEMIGROUPS_TEST_CASE("Matrix", "044", "exceptions", "[quick][element]") {
    using Mat               = NTPMat<>;
    using scalar_type       = typename Mat::scalar_type;
    NTPSemiring<> const* sr = new NTPSemiring<>(23, 1);
    REQUIRE(sr->scalar_one() == 1);
    REQUIRE(sr->scalar_zero() == 0);
    auto x = make<Mat>(sr, std::vector<std::vector<scalar_type>>());
    REQUIRE(x.number_of_cols() == x.number_of_rows());
    REQUIRE(x.number_of_cols() == 0);

    REQUIRE_THROWS_AS(make<Mat>(sr, {{2, 2, 0}, {0, 0, 140}, {1, 3, 1}}),
                      LibsemigroupsException);
    delete sr;
  }

  LIBSEMIGROUPS_TEST_CASE("Matrix", "045", "code coverage", "[quick]") {
    {
      BMat<> x(0, 0);
      x.transpose();
      REQUIRE(x == x);
      BMat<> y(2, 1);
      REQUIRE_THROWS_AS(matrix::pow(y, 2), LibsemigroupsException);
    }
    { REQUIRE_THROWS_AS(make<BMat<>>({{0, 1}, {0}}), LibsemigroupsException); }
    {
      BMat<> y(2, 2);
      std::fill(y.begin(), y.end(), 0);
      auto        r   = matrix::rows(y);
      auto const& val = r[0](0);
      REQUIRE(!val);
      REQUIRE(r[0] + r[1] == r[0]);
      std::ostringstream oss;
      oss << r[0];  // Does not do anything visible
      std::stringbuf buff;
      std::ostream   os(&buff);
      os << r[0];  // Also does not do anything visible
    }
    {
      using Row = BMat<2>::Row;
      BMat<2> x;
      std::fill(x.begin(), x.end(), 1);
      auto r = matrix::rows(x);
      Row  y(nullptr, r[0]);
      REQUIRE(y == Row({1, 1}));
    }
    {
      std::vector<std::vector<int>> v = {{true, false}, {true, false}};
      BMat<>                        x(v);
      REQUIRE(x == BMat<>({{true, false}, {true, false}}));
    }
    { BMat<> x; }
    {
      using Mat               = NTPMat<>;
      NTPSemiring<> const* sr = new NTPSemiring<>(23, 1);
      auto                 x  = Mat(sr, {{0, 1}, {0, 2}});
      Mat                  y  = x;
      REQUIRE(x == y);
      REQUIRE(&x != &y);

      Mat z(std::move(x));
      REQUIRE(z == y);
      REQUIRE(&z != &y);
      z(1, 0) = 2;
      y.swap(z);
      REQUIRE(y(1, 0) == 2);
      REQUIRE(z(1, 0) == 0);
      delete sr;
    }
  }

  LIBSEMIGROUPS_TEST_CASE("BMatFastest",
                          "046",
                          "check no throw",
                          "[quick][matrix]") {
    REQUIRE_NOTHROW(BMatFastest<3>({{0, 1}, {0, 1}}));
  }

  LIBSEMIGROUPS_TEST_CASE("Matrix",
                          "047",
                          "to_human_readable_repr",
                          "[quick][matrix]") {
    BMat<3> x({{0, 1, 0}, {0, 1, 0}, {0, 0, 0}});
    REQUIRE(to_human_readable_repr(x, "BMat<3>(")
            == "BMat<3>({{0, 1, 0},\n"
               "         {0, 1, 0},\n"
               "         {0, 0, 0}})");
    MinPlusMat<> y({{-2, 2, 0}, {-1, 0, 0}, {1, -3, POSITIVE_INFINITY}});
    REQUIRE(to_human_readable_repr(y, "MinPlusMat<>(")
            == "MinPlusMat<>({{-2,  2,  0},\n"
               "              {-1,  0,  0},\n"
               "              { 1, -3, +∞}})");
    REQUIRE(
        to_human_readable_repr(y, "Matrix(MatrixKind.MaxPlusMat, ", "", "[]")
        == "Matrix(MatrixKind.MaxPlusMat, [[-2,  2,  0],\n"
           "                               [-1,  0,  0],\n"
           "                               [ 1, -3, +∞]])");
    y = MinPlusMat<>({{-2, 2, 0}, {-1, 0, 0}, {1, -3, 666666}});
    REQUIRE(to_human_readable_repr(y, "MinPlusMat<>(")
            == "MinPlusMat<>({{    -2,      2,      0},\n"
               "              {    -1,      0,      0},\n"
               "              {     1,     -3, 666666}})");

    REQUIRE(to_human_readable_repr(y, "MinPlusMat<>(", "", "{}", 28)
            == "MinPlusMat<>({{-2, 2, 0},\n"
               "              {-1, 0, 0},\n"
               "              {1, -3, 666666}})");
    REQUIRE(
        to_human_readable_repr(y, "MinPlusMat<>(", "max-plus matrix", "{}", 20)
        == "<3x3 max-plus matrix>");
    MinPlusTruncSemiring<> const* sr = new MinPlusTruncSemiring<>(5);
    MinPlusTruncMat<> z(sr, {{0, 2, 0}, {1, 0, 0}, {1, 3, POSITIVE_INFINITY}});
    REQUIRE(to_human_readable_repr(z, "MinPlusTruncMat<>(")
            == "MinPlusTruncMat<>({{ 0,  2,  0},\n"
               "                   { 1,  0,  0},\n"
               "                   { 1,  3, +∞}})");
    delete sr;
  }

}  // namespace libsemigroups
