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

#include "libsemigroups/forest.hpp"

#include <algorithm>         // for fill
#include <cstddef>           // for size_t
#include <initializer_list>  // for initializer_list
#include <iosfwd>            // for ostream
#include <vector>            // for vector, allocator, operator==

#include "libsemigroups/bipart.hpp"
#include "libsemigroups/constants.hpp"  // for Max, operator!=, operator==
#include "libsemigroups/exception.hpp"  // for LIBSEMIGROUPS_EXCEPTION
#include "libsemigroups/types.hpp"      // for word_type

#include "libsemigroups/detail/string.hpp"  // for operator<<

// namespace forest {

//! Returns a \ref normal_form_iterator pointing at the first normal
//! form.
//!
//! Returns a const iterator pointing to the normal form of the first
//! class of the congruence represented by an instance of ToddCoxeter.
//! The order of the classes, and the normal form, that is returned are
//! controlled by standardize(order).
//!
//! \returns A value of type \ref normal_form_iterator.
//!
//! \exceptions
//! \no_libsemigroups_except
//!
// TODO(2): update the doc
// TODO(2) This is probably not very efficient, better to have a custom range
//
// where the word_type is cached, and next just deletes the end of the word
// and adds whatever's appropriate, rather than retracing the entire path
// every time
//
// inline auto paths(Forest const& f) {
//   rx::seq<Forest::node_type>() | rx::transform([&f](auto i) {
//     word_type w;
//     while (f.parent_no_checks(i) != UNDEFINED) {
//       w.push_back(f.label_no_checks(i));
//       i = f.parent_no_checks(i);
//     }
//     return w;
//   });
// }
//}  // namespace forest

namespace libsemigroups {
  Forest::~Forest() = default;

  Forest& Forest::init(size_t n) {
    _edge_label.resize(n);
    std::fill(std::begin(_edge_label),
              std::end(_edge_label),
              static_cast<size_t>(UNDEFINED));
    _parent.resize(n);
    std::fill(
        std::begin(_parent), std::end(_parent), static_cast<size_t>(UNDEFINED));
    return *this;
  }

  Forest& Forest::add_nodes(size_t n) {
    size_t const old_nr_nodes = number_of_nodes();
    try {
      _edge_label.insert(_edge_label.cend(), n, static_cast<size_t>(UNDEFINED));
      _parent.insert(_parent.cend(), n, static_cast<size_t>(UNDEFINED));
    } catch (...) {
      _edge_label.resize(old_nr_nodes);
      _parent.resize(old_nr_nodes);
      throw;
    }
    return *this;
  }

  void Forest::path_to_root_no_checks(word_type& w, node_type i) const {
    path_to_root_no_checks(std::back_inserter(w), i);
  }

  void Forest::throw_if_node_out_of_bounds(node_type v) const {
    if (v >= number_of_nodes()) {
      LIBSEMIGROUPS_EXCEPTION("node value out of bounds, expected value in "
                              "the range [0, {}), got {}",
                              number_of_nodes(),
                              v);
    }
  }

  std::string to_human_readable_repr(Forest const& f) {
    size_t const num_roots
        = std::count(f.parents().begin(), f.parents().end(), UNDEFINED);
    return fmt::format("<Forest with {} nodes, {} edges, and {} root{}>",
                       f.number_of_nodes(),
                       f.parents().size(),
                       num_roots,
                       num_roots == 1 ? "" : "s");
  }

}  // namespace libsemigroups
