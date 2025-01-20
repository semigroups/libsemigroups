//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2023 James D. Mitchell
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

// This file contains declarations of some ranges (in the sense of rx-ranges)
// that augment those in rx-ranges.

#ifndef LIBSEMIGROUPS_RANGES_HPP_
#define LIBSEMIGROUPS_RANGES_HPP_

#include <algorithm>    // for is_sorted
#include <cstddef>      // for size_t
#include <functional>   // for less
#include <ostream>      // for operator<<, ostream
#include <random>       // for mt19937, random_device, uniform_int_distribution
#include <type_traits>  // for enable_if_t
#include <utility>      // for forward

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#include <rx/ranges.hpp>  // for count, iterator_range
#pragma GCC diagnostic pop

#include "config.hpp"  // for PARSED_BY_DOXYGEN

#include "detail/fmt.hpp"     // for fmt::join
#include "detail/string.hpp"  // for detail::to_string

namespace libsemigroups {

  //! \defgroup ranges_group Ranges
  //!
  //! This page contains links to some classes and functions in
  //! `libsemigroups` that extend the functionality of range objects from
  //! [rx::ranges][].
  //!
  //! Please consult the [rx::ranges][] [README][] file for more information
  //! about ranges in general and the existing functionality of [rx::ranges][].
  //!
  //! Range objects are used extensively elsewhere in `libsemigroups`, and
  //! are returned by many functions.
  //!
  //! \note So that range classes in `libsemigroups` can be used efficiently
  //! with the functionality of [rx::ranges][], the usual naming conventions in
  //! `libsemigroups` are not used for the member functions of range objects.
  //! In particular, none of the member functions required by [rx::ranges][]
  //! check their arguments, but they do not have the suffix `_no_checks`.
  //!
  //! [rx::ranges]: https://github.com/simonask/rx-ranges/
  //! [README]: https://github.com/simonask/rx-ranges/blob/develop/README.md

  ////////////////////////////////////////////////////////////////////////
  // Custom ranges
  ////////////////////////////////////////////////////////////////////////

  //! \ingroup ranges_group
  //!
  //! \brief Aggregator for choosing a single random item in a range.
  //!
  //! Defined in `ranges.hpp`.
  //!
  //! This struct can be used to select a single random item from an input
  //! range. This can be extremely inefficient if the input range is large,
  //! because every item in the range prior to the randomly chosen one must be
  //! iterated through when `get` is called.
  //!
  //! \par Example
  //! \code
  //! auto wg = make<WordGraph<uint8_t>>(4, [[0, 1], [1, 0], [2, 2]]);
  //! Paths p(wg);
  //! p.source(0).max(10);
  //! p.count();            // returns 1023
  //! (p | Random()).get(); // returns random path in p
  //! \endcode
  // TODO(2) this should be able to emit any number of random items not only
  // one.
  struct Random {
    //! \brief Default constructor.
    //!
    //! Default constructor.
    Random() = default;

    //! \brief Default copy constructor.
    //!
    //! Default copy constructor.
    Random(Random const&) = default;

    //! \brief Default move constructor.
    //!
    //! Default move constructor.
    Random(Random&&) = default;

    //! \brief Default copy assignment operator.
    //!
    //! Default copy assignment operator.
    Random& operator=(Random const&) = default;

    //! \brief Default move assignment operator.
    //!
    //! Default move assignment operator.
    Random& operator=(Random&&) = default;

    template <typename InputRange>
    struct Range {
      using output_type = typename InputRange::output_type;

      // TODO(1) this is probably not correct, should depend on InputRange
      static constexpr bool is_finite = rx::is_finite_v<InputRange>;
      // TODO(1) this is probably not correct, should depend on InputRange
      static constexpr bool is_idempotent = rx::is_idempotent_v<InputRange>;

      bool                _at_end;
      InputRange          _input;
      mutable output_type _val;
      mutable bool        _val_set;

      explicit constexpr Range(InputRange&& input) noexcept
          : _at_end(false), _input(std::move(input)), _val(), _val_set(false) {}

      [[nodiscard]] output_type get() const noexcept;

      constexpr void next() noexcept {
        _at_end = true;
      }

      [[nodiscard]] constexpr bool at_end() const noexcept {
        return _at_end;
      }

      [[nodiscard]] constexpr size_t size_hint() const noexcept {
        if (!at_end()) {
          return 1;
        } else {
          return 0;
        }
      }
    };

#ifndef PARSED_BY_DOXYGEN
    template <typename InputRange>
    [[nodiscard]] constexpr auto operator()(InputRange&& input) const {
      using Inner = rx::get_range_type_t<InputRange>;
      return Range<Inner>(std::forward<InputRange>(input));
    }
#endif  // PARSED_BY_DOXYGEN
  };

  ////////////////////////////////////////////////////////////////////////
  // Comparison functions
  ////////////////////////////////////////////////////////////////////////

  //! \ingroup ranges_group
  //!
  //! \brief Check if a range is sorted according to \p comp.
  //!
  //! Defined in `ranges.hpp`.
  //!
  //! This function checks if a range is sorted according to the comparison
  //! function object \p comp, and is analogous to std::is_sorted.
  //!
  //! \tparam Range the type of the range to check for sortedness.
  //! \tparam Compare the type of the comparison function object.
  //!
  //! \param r the input range.
  //! \param comp the comparison function object which returns \c true if the
  //! first argument is less than the second.
  //!
  //! \returns \c true if the elements in the range are sorted in non-descending
  //! order according to \p comp and \c false otherwise.
  //!
  //! \note This function returns \c true for ranges of size \c 0 and \c 1.
  template <typename Range, typename Compare>
  constexpr bool is_sorted(Range r, Compare&& comp);

  //! \ingroup ranges_group
  //!
  //! \brief Check if a range is sorted according to std::less.
  //!
  //! Defined in `ranges.hpp`.
  //!
  //! This function checks if a range is sorted according to the comparison
  //! function object \p comp, and is analogous to std::is_sorted.
  //!
  //! \tparam Range the type of the range to check for sortedness.
  //!
  //! \param r the input range.
  //!
  //! \returns \c true if the elements in the range are sorted in non-descending
  //! order according to std::less and \c false otherwise.
  //!
  //! \note This function returns \c true for ranges of size \c 0 and \c 1.
  template <typename Range>
  constexpr bool is_sorted(Range r) {
    return std::is_sorted(r, std::less<>());
  }

  //! \ingroup ranges_group
  //!
  //! \brief Check two ranges for equality.
  //!
  //! Defined in `ranges.hpp`.
  //!
  //! This function checks if the ranges \p r1 and \p r2 are equal, in that
  //! they contain the same items (compared using `operator==`) in the same
  //! order. This is analogous to std::equal.
  //!
  //! \tparam Range1 the type of the first range.
  //! \tparam Range2 the type of the second range.
  //!
  //! \param r1 the first range.
  //! \param r2 the second range.
  //!
  //! \returns \c true if the ranges are equal and \c false if not.
  template <typename Range1, typename Range2>
  constexpr bool equal(Range1 r1, Range2 r2);

  //! \ingroup ranges_group
  //!
  //! \brief Check if a range is lexicographically less than another.
  //!
  //! Defined in `ranges.hpp`.
  //!
  //! This function checks if the range \p r1 is lexicographically less than
  //! the range \p r2, using `operator<` on the items in the ranges.
  //! This is analogous to std::lexicographical_compare.
  //!
  //! \tparam Range1 the type of the first range.
  //! \tparam Range2 the type of the second range.
  //!
  //! \param r1 the first range.
  //! \param r2 the second range.
  //!
  //! \returns \c true if the range \p r1 is lexicographically less than the
  //! range \p r2, and \c false it not.
  template <typename Range1, typename Range2>
  constexpr bool lexicographical_compare(Range1 r1, Range2 r2);

  //! \ingroup ranges_group
  //!
  //! \brief Check if a range is shortlex less than another.
  //!
  //! Defined in `ranges.hpp`.
  //!
  //! This function checks if the range \p r1 is shortlex less than
  //! the range \p r2, using `operator<` on the items in the ranges.
  //! This is analogous to \ref shortlex_compare.
  //!
  //! \tparam Range1 the type of the first range.
  //! \tparam Range2 the type of the second range.
  //!
  //! \param r1 the first range.
  //! \param r2 the second range.
  //!
  //! \returns \c true if the range \p r1 is shortlex less than the
  //! range \p r2, and \c false it not.
  template <typename Range1, typename Range2>
  constexpr bool shortlex_compare(Range1 r1, Range2 r2);

  // TODO: recursive_path_compare?

  ////////////////////////////////////////////////////////////////////////
  // Custom combinators
  ////////////////////////////////////////////////////////////////////////

  //! \ingroup ranges_group
  //!
  //! \brief Chain objects (const references)
  //!
  //! Defined in `ranges.hpp`.
  //!
  //! This function is just a wrapper for \c rx::chain but taking objects by
  //! reference rather than using iterators for convenience, and it returns a
  //! range objects where the items in \p x come before any of the items in \p
  //! y.
  //!
  //! \tparam S the type of the first argument.
  //! \tparam T the type of the second argument.
  //!
  //! \param x the first object to chain.
  //! \param y the second object to chain.
  //!
  //! \returns A range object.
  template <typename S, typename T>
  auto chain(S const& x, T const& y) {
    return rx::chain(rx::iterator_range(x.cbegin(), x.cend()),
                     rx::iterator_range(y.cbegin(), y.cend()));
  }

  //! \ingroup ranges_group
  //!
  //! \brief Enumerate an object (by const reference)
  //!
  //! Defined in `ranges.hpp`.
  //!
  //! This function is just a wrapper for \c rx::enumerate but taking objects by
  //! reference rather than using iterators for convenience, and it returns a
  //! a pair consisting of the index of an item and that item. This is similar
  //! to the function of the same name in python.
  //!
  //! \tparam T the type of the argument.
  //!
  //! \param thing the object to enumerate.
  //!
  //! \returns A range object.
  template <typename T>
  auto enumerate(T const& thing) {
    return rx::enumerate(rx::iterator_range(thing));
  }

  ////////////////////////////////////////////////////////////////////////
  // String representation
  ////////////////////////////////////////////////////////////////////////

  //! \ingroup ranges_group
  //!
  //! \brief Insert a range into an output stream.
  //!
  //! Defined in `ranges.hpp`.
  //!
  //! This function inserts a human readable representation of a range object
  //! into a std::ostream by applying std::to_string to each item in the range.
  //!
  //! \tparam Range the type of the range (must satisfy
  //! `rx::is_input_or_sink_v<Range>`).
  //!
  //! \param os the output stream
  //! \param r the range object
  //!
  //! \returns A reference to the first parameter \p os.
  template <typename Range,
            typename = std::enable_if_t<rx::is_input_or_sink_v<Range>>>
  std::ostream& operator<<(std::ostream& os, Range r);

}  // namespace libsemigroups

#include "ranges.tpp"

#endif  // LIBSEMIGROUPS_RANGES_HPP_
