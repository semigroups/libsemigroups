//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2020 James D. Mitchell
//                    + TODO the other guys
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

#include "libsemigroups/freeband.hpp"

#include <algorithm>  // for max_element
#include <iostream>   // for cout
#include <vector>     // for vector

#include "libsemigroups/constants.hpp"
#include "libsemigroups/libsemigroups-debug.hpp"
#include "libsemigroups/string.hpp"

namespace libsemigroups {
  namespace {
    bool is_standardized(word_type const& w) {
      size_t m = 0;
      for (auto const& l : w) {
        if (l > m + 1) {
          return false;
        } else if (l > m) {
          ++m;
        }
      }
      return true;
    }

    void standardize(word_type& x) {
      if (x.empty()) {
        return;
      }

      size_t              distinct_chars = 0;
      size_t const        N = *std::max_element(x.cbegin(), x.cend()) + 1;
      std::vector<size_t> lookup(N, N);
      lookup[x[0]] = 0;
      x[0]         = 0;
      for (size_t i = 1; i < x.size(); ++i) {
        if (lookup[x[i]] == N) {
          lookup[x[i]] = ++distinct_chars;
        }
        x[i] = lookup[x[i]];
      }
    }

  }  // namespace

    word_type right(word_type const& w, size_t const k) {
      LIBSEMIGROUPS_ASSERT(is_standardized(w));
      word_type out;
      if (w.empty()) {
        return out;
      }
      size_t              j            = 0;
      size_t              content_size = 0;
      size_t const        N = *std::max_element(w.cbegin(), w.cend()) + 1;
      std::vector<size_t> multiplicity(N + 1, 0);
      for (size_t i = 0; i < w.size(); ++i) {
        if (i != 0) {
          LIBSEMIGROUPS_ASSERT(multiplicity[w[i - 1]] != 0);
          --multiplicity[w[i - 1]];
          if (multiplicity[w[i - 1]] == 0) {
            --content_size;
          }
        }
        while (j < w.size()
              && (multiplicity[w[j]] != 0 || content_size < k)) {
          if (multiplicity[w[j]] == 0) {
            ++content_size;
          }
          ++multiplicity[w[j]];
          ++j;
        }
        out.push_back(content_size == k ? j - 1 : UNDEFINED);
      }
      return out;
    }
  bool freeband_equal_to(word_type& x, word_type& y) {
    standardize(x);
    standardize(y);
    right(x, 4);
    // std::cout << "x = " << right(x, 4) << std::endl;
    // std::cout << "y = " << y << std::endl;
    return false;
  }

}  // namespace libsemigroups
