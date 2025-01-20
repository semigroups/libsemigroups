//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2019-2024 James D. Mitchell
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

// This file contains implementations of the PBR class.

#include "libsemigroups/pbr.hpp"

#include <algorithm>  // for all_of, fill
#include <ostream>    // for operator<<, cha...
#include <string>     // for operator+, char...
#include <thread>     // for thread

#include "libsemigroups/debug.hpp"      // for LIBSEMIGROUPS_A...
#include "libsemigroups/exception.hpp"  // for LIBSEMIGROUPS_E...

#include "libsemigroups/detail/containers.hpp"  // for DynamicArray2
#include "libsemigroups/detail/string.hpp"      // for to_string

namespace libsemigroups {

  namespace {
    void throw_if_invalid_side(PBR::vector_type<int32_t> side,
                               std::string               position) {
      size_t n = side.size();
      for (std::vector<int32_t> const& vec : side) {
        for (int32_t x : vec) {
          if (x == 0 || x < -static_cast<int32_t>(n)
              || x > static_cast<int32_t>(n)) {
            LIBSEMIGROUPS_EXCEPTION(
                "value out of bounds in the {} argument, expected values in "
                "[-{}, -1] or [1, {}] but found {}",
                position,
                n,
                n,
                x);
          }
        }
      }
    }

    void throw_if_invalid_left_right(PBR::vector_type<int32_t> left,
                                     PBR::vector_type<int32_t> right) {
      size_t n = left.size();

      if (n != right.size()) {
        LIBSEMIGROUPS_EXCEPTION("the two vectors must have the same length");
      }
      if (n > 0x40000000) {
        LIBSEMIGROUPS_EXCEPTION("too many points!");
      }

      throw_if_invalid_side(left, std::string("1st"));
      throw_if_invalid_side(right, std::string("2nd"));
    }

    std::vector<std::vector<int32_t>>
    sorted_side(PBR::vector_type<int32_t> side) {
      std::vector<std::vector<int32_t>> out(side);
      for (std::vector<int32_t>& vec : out) {
        if (!std::is_sorted(vec.cbegin(), vec.cend())) {
          std::sort(vec.begin(), vec.end());
        }
      }
      return out;
    }

    void process_side_no_checks(std::vector<std::vector<uint32_t>>& out,
                                PBR::vector_type<int32_t>           side) {
      std::vector<uint32_t> v;
      size_t                n = side.size();

      for (std::vector<int32_t> const& vec : side) {
        v.clear();
        for (int32_t x : vec) {
          if (x > 0) {
            v.push_back(static_cast<uint32_t>(x - 1));
          }
        }
        // We have to go backwards through the vector to add the negative
        // entries, so that users can input those negative entries in the
        // natural order
        for (auto it = vec.rbegin(); it < vec.rend(); ++it) {
          if (*it < 0) {
            v.push_back(static_cast<uint32_t>(n - *it - 1));
          }
        }
        out.push_back(v);
      }
    }

    std::vector<std::vector<uint32_t>>
    process_left_right_no_checks(PBR::vector_type<int32_t> left,
                                 PBR::vector_type<int32_t> right) {
      std::vector<std::vector<uint32_t>> out;
      process_side_no_checks(out, left);
      process_side_no_checks(out, right);
      return out;
    }

    void unite_rows(detail::DynamicArray2<bool>& out,
                    detail::DynamicArray2<bool>& tmp,
                    size_t const&                i,
                    size_t const&                j) {
      for (size_t k = 0; k < out.number_of_cols(); k++) {
        out.set(i, k, (out.get(i, k) || tmp.get(j, k + 1)));
      }
    }

    void y_dfs(std::vector<bool>&,
               std::vector<bool>&,
               detail::DynamicArray2<bool>&,
               uint32_t const&,
               uint32_t const&,
               PBR const* const,
               PBR const* const,
               size_t const&);

    void x_dfs(std::vector<bool>&           x_seen,
               std::vector<bool>&           y_seen,
               detail::DynamicArray2<bool>& tmp,
               uint32_t const&              n,
               uint32_t const&              i,
               PBR const* const             x,
               PBR const* const             y,
               size_t const&                adj) {
      if (!x_seen[i]) {
        x_seen[i] = true;
        for (auto const& j : (*x)[i]) {
          if (j < n) {
            tmp.set(adj, j + 1, true);
          } else {
            y_dfs(x_seen, y_seen, tmp, n, j - n, x, y, adj);
          }
        }
      }
    }

    void y_dfs(std::vector<bool>&           x_seen,
               std::vector<bool>&           y_seen,
               detail::DynamicArray2<bool>& tmp,
               uint32_t const&              n,
               uint32_t const&              i,
               PBR const* const             x,
               PBR const* const             y,
               size_t const&                adj) {
      if (!y_seen[i]) {
        y_seen[i] = true;
        for (auto const& j : (*y)[i]) {
          if (j >= n) {
            tmp.set(adj, j + 1, true);
          } else {
            x_dfs(x_seen, y_seen, tmp, n, j + n, x, y, adj);
          }
        }
      }
    }

  }  // namespace

  std::vector<std::vector<uint32_t>>
  detail::process_left_right(PBR::vector_type<int32_t> left,
                             PBR::vector_type<int32_t> right) {
    throw_if_invalid_left_right(left, right);
    PBR::vector_type<int32_t> sorted_left(sorted_side(left));
    PBR::vector_type<int32_t> sorted_right(sorted_side(right));
    return process_left_right_no_checks(sorted_left, sorted_right);
  }

  PBR operator*(PBR const& x, PBR const& y) {
    PBR xy(x.degree());
    xy.product_inplace_no_checks(x, y);
    return xy;
  }

  PBR pbr::one(size_t n) {
    std::vector<std::vector<uint32_t>> adj;
    adj.reserve(2 * n);
    for (uint32_t i = 0; i < 2 * n; i++) {
      adj.push_back(std::vector<uint32_t>());
    }
    for (uint32_t i = 0; i < n; i++) {
      adj[i].push_back(i + n);
      adj[i + n].push_back(i);
    }
    return PBR(adj);
  }

  PBR pbr::one(PBR const& x) {
    return pbr::one(x.degree());
  }

  void pbr::throw_if_not_even_length(PBR const& x) {
    size_t n(x.number_of_points());
    if (n % 2 == 1) {
      LIBSEMIGROUPS_EXCEPTION("expected argument of even length, found {}",
                              detail::to_string(n));
    }
  }

  void pbr::throw_if_entry_out_of_bounds(PBR const& x) {
    size_t n(x.number_of_points());
    for (size_t u = 0; u < n; ++u) {
      for (auto const& v : x[u]) {
        if (v >= n) {
          LIBSEMIGROUPS_EXCEPTION(
              "entry out of bounds, vertex " + detail::to_string(u)
              + " is adjacent to " + detail::to_string(v)
              + ", should be less than " + detail::to_string(n));
        }
      }
    }
  }

  void pbr::throw_if_adjacencies_unsorted(PBR const& x) {
    size_t n(x.number_of_points());
    for (size_t u = 0; u < n; ++u) {
      if (!std::is_sorted(x[u].cbegin(), x[u].cend())) {
        LIBSEMIGROUPS_EXCEPTION("the adjacencies of vertex {} are unsorted",
                                detail::to_string(u));
      }
    }
  }

  [[nodiscard]] std::string to_human_readable_repr(PBR const& x) {
    // TODO(2) allow different braces
    // TODO(now) Make this better, probably by including some data from
    // x._vector
    return fmt::format("<PBR of degree {}>", x.degree());
  }

  ////////////////////////////////////////////////////////////////////////
  // Partitioned binary relations (PBRs)
  ////////////////////////////////////////////////////////////////////////

  PBR::PBR(PBR::vector_type<uint32_t> vec) : _vector(vec) {}

  PBR::PBR(PBR::initializer_list_type<uint32_t> vec) : _vector(vec) {}

  PBR::PBR(size_t degree)
      : PBR(std::vector<std::vector<uint32_t>>(degree * 2,
                                               std::vector<uint32_t>())) {}

  PBR::PBR(PBR::initializer_list_type<int32_t> left,
           PBR::initializer_list_type<int32_t> right)
      : PBR(process_left_right_no_checks(left, right)) {}

  PBR::PBR(PBR::vector_type<int32_t> left, PBR::vector_type<int32_t> right)
      : PBR(process_left_right_no_checks(left, right)) {}

  std::ostringstream& operator<<(std::ostringstream& os, PBR const& pbr) {
    if (pbr.degree() == 0) {
      os << "{}";
      return os;
    }
    os << "{";
    for (size_t i = 0; i < pbr.degree() * 2 - 1; ++i) {
      os << "{";
      if (!pbr[i].empty()) {
        for (size_t j = 0; j + 1 < pbr[i].size(); ++j) {
          os << pbr[i][j] << ", ";
        }
        os << detail::to_string(pbr[i].back());
      }
      os << "}, ";
    }

    os << "{";
    if (!pbr[2 * pbr.degree() - 1].empty()) {
      for (size_t j = 0; j + 1 < pbr[2 * pbr.degree() - 1].size(); ++j) {
        os << pbr[2 * pbr.degree() - 1][j] << ", ";
      }
      os << detail::to_string(pbr[2 * pbr.degree() - 1].back());
    }
    os << "}}";
    return os;
  }

  std::ostream& operator<<(std::ostream& os, PBR const& pbr) {
    os << detail::to_string(pbr);
    return os;
  }

  size_t PBR::degree() const noexcept {
    return _vector.size() / 2;
  }

  size_t PBR::number_of_points() const noexcept {
    return _vector.size();
  }

  void PBR::product_inplace_no_checks(PBR const& xx,
                                      PBR const& yy,
                                      size_t     thread_id) {
    LIBSEMIGROUPS_ASSERT(xx.degree() == yy.degree());
    LIBSEMIGROUPS_ASSERT(xx.degree() == this->degree());
    LIBSEMIGROUPS_ASSERT(&xx != this && &yy != this);

    static std::vector<std::vector<bool>> _x_seen(
        std::thread::hardware_concurrency() + 1);
    static std::vector<std::vector<bool>> _y_seen(
        std::thread::hardware_concurrency() + 1);
    static std::vector<detail::DynamicArray2<bool>> _out(
        std::thread::hardware_concurrency() + 1);
    static std::vector<detail::DynamicArray2<bool>> _tmp(
        std::thread::hardware_concurrency() + 1);

    PBR const& x = static_cast<PBR const&>(xx);
    PBR const& y = static_cast<PBR const&>(yy);

    uint32_t const n = this->degree();

    std::vector<bool>&           x_seen = _x_seen.at(thread_id);
    std::vector<bool>&           y_seen = _y_seen.at(thread_id);
    detail::DynamicArray2<bool>& tmp    = _tmp.at(thread_id);
    detail::DynamicArray2<bool>& out    = _out.at(thread_id);

    if (x_seen.size() != 2 * n) {
      x_seen.clear();
      x_seen.resize(2 * n, false);
      y_seen.clear();
      y_seen.resize(2 * n, false);
      out.clear();
      out.add_cols(2 * n);
      out.add_rows(2 * n);
      tmp.clear();
      tmp.add_cols(2 * n + 1);
    } else {
      std::fill(x_seen.begin(), x_seen.end(), false);
      std::fill(y_seen.begin(), y_seen.end(), false);
      std::fill(out.begin(), out.end(), false);
      std::fill(tmp.begin(), tmp.end(), false);
    }

    for (size_t i = 0; i < n; i++) {
      for (auto const& j : x[i]) {
        if (j < n) {
          out.set(i, j, true);
        } else if (j < tmp.number_of_rows() && tmp.get(j, 0)) {
          unite_rows(out, tmp, i, j);
        } else {
          if (j >= tmp.number_of_rows()) {
            tmp.add_rows(j - tmp.number_of_rows() + 1);
          }
          tmp.set(j, 0, true);
          x_seen[i] = true;
          y_dfs(x_seen, y_seen, tmp, n, j - n, &x, &y, j);
          unite_rows(out, tmp, i, j);
          std::fill(x_seen.begin(), x_seen.end(), false);
          std::fill(y_seen.begin(), y_seen.end(), false);
        }
        if (std::all_of(out.begin_row(i), out.end_row(i), [](bool val) {
              return val;
            })) {
          break;
        }
      }
    }

    for (size_t i = n; i < 2 * n; i++) {
      for (auto const& j : y[i]) {
        if (j >= n) {
          out.set(i, j, true);
        } else if (j < tmp.number_of_rows() && tmp.get(j, 0)) {
          unite_rows(out, tmp, i, j);
        } else {
          if (j >= tmp.number_of_rows()) {
            tmp.add_rows(j - tmp.number_of_rows() + 1);
          }
          tmp.set(j, 0, true);
          y_seen[i] = true;
          x_dfs(x_seen, y_seen, tmp, n, j + n, &x, &y, j);
          unite_rows(out, tmp, i, j);
          std::fill(x_seen.begin(), x_seen.end(), false);
          std::fill(y_seen.begin(), y_seen.end(), false);
        }
        if (std::all_of(out.begin_row(i), out.end_row(i), [](bool val) {
              return val;
            })) {
          break;
        }
      }
    }

    for (size_t i = 0; i < 2 * n; i++) {
      _vector[i].clear();
      for (size_t j = 0; j < 2 * n; j++) {
        if (out.get(i, j)) {
          _vector[i].push_back(j);
        }
      }
    }
  }

  void PBR::product_inplace(PBR const& xx, PBR const& yy, size_t thread_id) {
    if (xx.degree() != yy.degree()) {
      LIBSEMIGROUPS_EXCEPTION(
          "the degree of the first argument ({}) is not equal to the degree "
          "of the second argument ({})",
          xx.degree(),
          yy.degree());
    }
    if (xx.degree() != this->degree()) {
      LIBSEMIGROUPS_EXCEPTION(
          "the degree of the this ({}) is not equal to the degree "
          "of the first argument ({})",
          this->degree(),
          xx.degree());
    }
    if (&xx == this || &yy == this) {
      LIBSEMIGROUPS_EXCEPTION(
          "the address of this is the same as that of one of the arguments, "
          "expected this to be distinct from the arguments");
    }
    pbr::throw_if_invalid(xx);
    pbr::throw_if_invalid(yy);
    this->product_inplace_no_checks(xx, yy, thread_id);
  }

  std::vector<uint32_t>& PBR::at(size_t i) {
    if (i >= this->number_of_points()) {
      LIBSEMIGROUPS_EXCEPTION(
          "index out of range, expected values in [0, {}) but found {}",
          this->number_of_points(),
          i);
    }
    return this->operator[](i);
  }

  std::vector<uint32_t> const& PBR::at(size_t i) const {
    if (i >= this->number_of_points()) {
      LIBSEMIGROUPS_EXCEPTION(
          "index out of range, expected values in [0, {}) but found {}",
          this->number_of_points(),
          i);
    }
    return this->operator[](i);
  }

}  // namespace libsemigroups
