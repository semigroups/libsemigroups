//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2021-2024 James D. Mitchell
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

// This file contains the declaration of the Bipartition and Blocks classes.

#ifndef LIBSEMIGROUPS_BIPART_HPP_
#define LIBSEMIGROUPS_BIPART_HPP_

// TODO(2)
// * benchmarks
// * use Duf/Suf where possible (later?)
// * template like transformations/pperms etc (later?)

#include <algorithm>         // for max
#include <cstddef>           // for size_t
#include <cstdint>           // for uint32_t, int32_t
#include <cstdlib>           // for abs
#include <initializer_list>  // for initializer_list
#include <string_view>       // for string_view
#include <type_traits>       // for decay_t, false_type, is_signed, true_type
#include <unordered_set>     // for unordered_set
#include <vector>            // for vector

#include "adapters.hpp"   // for Hash
#include "debug.hpp"      // for LIBSEMIGROUPS_ASSERT
#include "exception.hpp"  // for LIBSEMIGROUPS_EXCEPTION
#include "types.hpp"      // for enable_if_is_same

#include "detail/fmt.hpp"

namespace libsemigroups {
  //! \defgroup bipart_group Bipartitions
  //!
  //! Defined `bipart.hpp`.
  //!
  //! This page contains an overview of the functionality in `libsemigroups`
  //! for bipartitions and blocks.
  //!
  //! Helper functions for bipartitions and blocks are documented:
  //! * \ref libsemigroups::bipartition "Helper functions for bipartitions"
  //! * \ref libsemigroups::blocks "Helper functions for blocks"

  // Forward decls
  class Bipartition;
  class Blocks;

  //! \ingroup bipart_group
  //!
  //! \brief Namespace for Blocks helper functions.
  //!
  //! This namespace contains helper functions for the Blocks class.
  namespace blocks {
    //! \brief Validate a Blocks object.
    //!
    //! This function validates a Blocks object, and throws an exception if the
    //! object is not valid.
    //!
    //! \param x the blocks object to validate.
    //!
    //! \throws LibsemigroupsException if \p x is invalid.
    void validate(Blocks const& x);

    //! \brief Return the underlying partition of a Blocks object.
    //!
    //! The *underlying partition* of a Blocks object \c x is the
    //! partition of a subset \f$P\f$ of \f$\{-n, \ldots, -1\}\cup \{1, \ldots,
    //! n\}\f$ such that:
    //! * \f$\{|x|\mid x\in P\} = \{1, \ldots, n\}\f$;
    //! * a block of the partition consists of negative numbers if and only if
    //! the corresponding block  of \c x is a transverse block.
    //!
    //! \param x the Blocks object.
    //!
    //! \returns
    //! A vector of vectors of integers.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(n)\f$ where \f$n\f$ is the degree().
    [[nodiscard]] std::vector<std::vector<int32_t>>
    underlying_partition(Blocks const& x);

  }  // namespace blocks

  //! \ingroup bipart_group
  //!
  //! \brief Namespace for Bipartition helper functions.
  //!
  //! This namespace contains helper functions for the Bipartition class.
  namespace bipartition {
    //! \brief Return the identity bipartition with the same degree as the
    //! given bipartition.
    //!
    //! \param f the bipartition.
    //! Returns the identity bipartition of degree equal to `f.degree()`.
    //!
    //! The *identity bipartition* of degree \f$n\f$ has blocks \f$\{i, -i\}\f$
    //! for all \f$i\in \{0, \ldots, n - 1\}\f$. This member function returns a
    //! new identity bipartition of degree equal to the degree of \c this.
    //!
    //! \returns A newly constructed Bipartition.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    [[nodiscard]] Bipartition one(Bipartition const& f);

    //! \brief Return the underlying partition of a Bipartition object.
    //!
    //! The *underlying partition* of a bipartition \c x is the
    //! partition of a subset \f$P\f$ of \f$\{-n, \ldots, -1\}\cup \{1, \ldots,
    //! n\}\f$ such that:
    //! * \f$\{|x|\mid x\in P\} = \{1, \ldots, n\}\f$;
    //! * a block of the partition consists of negative numbers if and only if
    //! the corresponding block of \c x is a transverse block.
    //!
    //! \param x the bipartition.
    //!
    //! \returns
    //! A vector of vectors of integers.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(n)\f$ where \f$n\f$ is the degree().
    [[nodiscard]] std::vector<std::vector<int32_t>>
    underlying_partition(Bipartition const&);

    //! \brief Validates a bipartition.
    //!
    //! This function validates a Bipartition object, and throws an exception if
    //! the object is not valid.
    //!
    //! \param x the bipartition.
    //!
    //! \throws LibsemigroupsException if \p x is invalid.
    void validate(Bipartition const& x);

  }  // namespace bipartition

  namespace detail {

    template <typename T>
    struct IsBipartitionHelper : std::false_type {};

    template <>
    struct IsBipartitionHelper<Bipartition> : std::true_type {};

  }  // namespace detail

  //! \ingroup bipart_group
  //! \brief Helper variable template.
  //!
  //! The value of this variable is \c true if the template parameter \p T
  //! decays to \ref Bipartition.
  //!
  //! \tparam T a type.
  template <typename T>
  static constexpr bool IsBipartition
      = detail::IsBipartitionHelper<std::decay_t<T>>::value;

  namespace detail {

    //! No doc
    template <typename BipartitionOrBlocks, typename Scalar>
    static void validate_args(std::vector<std::vector<Scalar>> const& blocks) {
      static_assert(std::is_same_v<BipartitionOrBlocks, Bipartition>
                    || std::is_same_v<BipartitionOrBlocks, Blocks>);
      static_assert(std::is_signed<Scalar>::value,
                    "the 2nd template parameter Scalar must be signed");
      int32_t offset = 2;
      if (!IsBipartition<BipartitionOrBlocks>) {
        offset = 1;
      }
      int32_t                    m   = 0;
      int32_t                    deg = 0;
      std::unordered_set<Scalar> vals;
      for (size_t i = 0; i < blocks.size(); ++i) {
        auto const& block = blocks[i];
        if (block.empty()) {
          LIBSEMIGROUPS_EXCEPTION("the argument (blocks) is invalid, "
                                  "expected all blocks to be non-empty, but "
                                  "found empty block in position {}",
                                  i);
        }
        bool positive = block[0] >= 0;

        for (size_t j = 0; j < block.size(); ++j) {
          auto x = block[j];
          vals.insert(x);
          if (x == 0) {
            LIBSEMIGROUPS_EXCEPTION("the argument (blocks) is invalid, "
                                    "expected non-zero values but found 0 in "
                                    "position {} of the block with index {}",
                                    j,
                                    i);
          } else if (!IsBipartition<BipartitionOrBlocks> && positive && x < 0) {
            LIBSEMIGROUPS_EXCEPTION(
                "the argument (blocks) is invalid, expected every value in the "
                "block with index {} to be {}tive, but found {} in position {}",
                i,
                positive ? "posi" : "nega",
                x,
                j);
          }
          x = std::abs(x);

          m = std::max(x, m);
          deg++;
        }
      }

      if (m >= static_cast<int32_t>(0x40000000)) {
        LIBSEMIGROUPS_EXCEPTION(
            "too many points, expected at most {}, found {}", 0x40000000, m);
      } else if (deg != offset * m || vals.size() != size_t(deg)) {
        std::string range, prefix;
        if (IsBipartition<BipartitionOrBlocks>) {
          prefix = "the union of";
          range  = fmt::format("{{{}, ..., -1, 1, ..., {}}}", -m, m);
        } else {
          prefix
              = "the set consisting of the absolute values of the entries in ";
          range = fmt::format("[1, {}]", m);
        }
        LIBSEMIGROUPS_EXCEPTION("{} the given blocks is not {},"
                                " only {} values were given",
                                prefix,
                                range,
                                deg);
      }
    }

    //! No doc
    template <typename BipartitionOrBlocks, typename Scalar>
    static void
    validate_args(std::initializer_list<std::vector<Scalar>> const& blocks) {
      std::vector<std::vector<Scalar>> arg(blocks);
      validate_args<BipartitionOrBlocks>(arg);
    }

    //! No doc
    template <typename BipartitionOrBlocks, typename Scalar>
    static void validate_args(std::vector<Scalar> const&) {
      // checks for this argument type are done in validate
    }

    //! No doc
    template <typename BipartitionOrBlocks, typename Scalar>
    static void validate_args(std::initializer_list<Scalar> const&) {
      // checks for this argument type are done in validate
    }
  }  // namespace detail

  //! \ingroup bipart_group
  //!
  //! \brief A Blocks object represents a signed partition of the set
  //! \f$\{0, \ldots, n - 1\}\f$.
  //!
  //! Defined in `bipart.hpp`.
  //!
  //! It is possible to associate to every Bipartition a pair of blocks,
  //! Bipartition::left_blocks() and Bipartition::right_blocks(), which
  //! determine the Green's \f$\mathscr{L}\f$- and \f$\mathscr{R}\f$-classes of
  //! the Bipartition in the monoid of all bipartitions. This is the purpose of
  //! this class.
  //!
  //! The Blocks class is not currently used widely in `libsemigroups`
  //! but are used extensively in the GAP package
  //! [Semigroups package for GAP](https://semigroups.github.io/Semigroups/).
  class Blocks {
   private:
    std::vector<uint32_t> _blocks;
    std::vector<bool>     _lookup;

   public:
    //! \brief Type for const iterators pointing to the transverse block lookup.
    //!
    //! Type for const iterators pointing to the transverse block lookup.
    using lookup_const_iterator = std::vector<bool>::const_iterator;

    //! \brief Type for iterators pointing to the index of the block.
    //!
    //! Type for iterators pointing to the index of the block.
    using iterator = typename std::vector<uint32_t>::iterator;

    //! \brief Type for const iterators pointing to the index of the block.
    //!
    //! Type for const iterators pointing to the index of the block.
    using const_iterator = typename std::vector<uint32_t>::const_iterator;

    //! \brief Constructs a blocks object of size 0.
    //!
    //! Default construct an empty Blocks object.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    Blocks() noexcept = default;

    //! \brief Construct a blocks object from iterators.
    //!
    //! The degree of the blocks object constructed is `last - first / 2`.
    //!
    //! \param first iterator pointing to the index of the block containing the
    //! first point.
    //! \param last iterator pointing one past the index of the block containing
    //! the last point.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Linear in `last - first`.
    //!
    //! \warning
    //! No checks are made on the validity of the arguments to this function.
    //!
    //! \sa validate(Blocks const&)
    Blocks(const_iterator first, const_iterator last);

    //! \brief Constructs a blocks object of given degree.
    //!
    //! \param degree the degree of the blocks object to construct.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Linear in \p degree.
    explicit Blocks(size_t degree) : _blocks(degree), _lookup() {}

    //! \brief Constructs a Blocks object from a vector of vectors of integers.
    //!
    //! This function constructs a Blocks object from a vector of vectors of
    //! (signed) integers, so that the blocks consisting of negative values are
    //! transverse and those consisting of positive values are not.
    //!
    //! \param blocks the blocks.
    //!
    //! \throws LibsemigroupsException if the set consisting of the absolute
    //! values of the entries in \p blocks is not \f${1, \ldots, n\}\f$ where
    //! \f$n\f$ is the maximum such value.
    //! \throws LibsemigroupsException if \c 0 is an item in any block.
    //! \throws LibsemigroupsException if any block is empty.
    //! \throws LibsemigroupsException if any block contains both negative and
    //! positive values.
    //! \throws LibsemigroupsException if the constructed Blocks object is not
    //! valid.
    //!
    //! \complexity Linear in the sum of the sizes of the vectors in \p blocks.
    explicit Blocks(std::vector<std::vector<int32_t>> const& blocks);

    //! \brief Default copy assignment operator.
    //!
    //! Default copy assignment operator.
    Blocks& operator=(Blocks const&) = default;

    //! \brief Default move assignment operator.
    //!
    //! Default move assignment operator.
    Blocks& operator=(Blocks&&) = default;

    //! \brief Default copy constructor.
    //!
    //! Default copy constructor.
    Blocks(Blocks const& copy) = default;

    //! \brief Default move constructor.
    //!
    //! Default move constructor.
    Blocks(Blocks&& copy) = default;

    ~Blocks();

    //! \brief Set whether or not the block containing a point is transverse.
    //!
    //! This function can be used to set whether or not the block containing \p
    //! i is transverse.
    //!
    //! \param i the point.
    //! \param val whether or not the block containing \p i is transverse.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    //!
    //! \warning
    //! No checks are made on the validity of the arguments to this function.
    Blocks& is_transverse_block_no_checks(size_t i, bool val) {
      LIBSEMIGROUPS_ASSERT(i < _lookup.size());
      _lookup[i] = val;
      return *this;
    }

    //! \brief Set whether or not the block containing a point is transverse.
    //!
    //! This function can be used to set whether or not the block containing \p
    //! i is transverse.
    //!
    //! \param i the point.
    //! \param val whether or not the block containing \p i is transverse.
    //!
    //! \throws LibsemigroupsException if \p i is not in the range \f$[0, n)\f$
    //! where \f$n\f$ is the return value of \ref number_of_blocks.
    //!
    //! \complexity
    //! Constant.
    Blocks& is_transverse_block(size_t i, bool val) {
      throw_if_class_index_out_of_range(i);
      return is_transverse_block_no_checks(i, val);
    }

    //! \brief Check if a block is a transverse block.
    //!
    //! This function returns \c true if the block with index \p index is a
    //! transverse (or signed) block and it returns \c false if it is not
    //! transverse (or unsigned).
    //!
    //! \param index the index of a block.
    //!
    //! \returns Whether or not the given block is transverse.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    //!
    //! \warning
    //! No checks are made on the validity of the arguments to this function.
    [[nodiscard]] bool is_transverse_block_no_checks(size_t index) const {
      LIBSEMIGROUPS_ASSERT(index < _lookup.size());
      return _lookup[index];
    }

    //! \brief Check if a block is a transverse block.
    //!
    //! This function returns \c true if the block with index \p index is a
    //! transverse (or signed) block and it returns \c false if it is not
    //! transverse (or unsigned).
    //!
    //! \param index the index of a block.
    //!
    //! \returns Whether or not the given block is transverse.
    //!
    //! \throws LibsemigroupsException if \p i is not in the range \f$[0, n)\f$
    //! where \f$n\f$ is the return value of \ref number_of_blocks.
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] bool is_transverse_block(size_t index) const {
      throw_if_class_index_out_of_range(index);
      return is_transverse_block_no_checks(index);
    }

    //! \brief Set the block that a point belongs to.
    //!
    //! This function can be used to set the block containing the point \p i (to
    //! equal \p val).
    //!
    //! \param i the point.
    //! \param val the block that \p i should belong to.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst linear in `degree()`.
    //!
    //! \warning
    //! No checks are made on the validity of the arguments to this function.
    Blocks& block_no_checks(size_t i, uint32_t val);

    //! \brief Set the block that a point belongs to.
    //!
    //! This function can be used to set the block containing the point \p i (to
    //! equal \p val).
    //!
    //! \param i the point.
    //! \param val the block that \p i should belong to.
    //!
    //! \throws LibsemigroupsException if \p i is not in the range \f$[0, n)\f$
    //! where \f$n\f$ is the return value of \ref number_of_blocks.
    //!
    //! \complexity
    //! At worst linear in `degree()`.
    Blocks& block(size_t i, uint32_t val);

    //! \brief Compare two blocks objects for equality.
    //!
    //! Two Blocks objects are equal if and only if their underlying signed
    //! partitions are equal. It is ok to compare blocks of different
    //! degree with this operator.
    //!
    //! \param that the Blocks instance for comparison.
    //!
    //! \returns Whether or not `*this` equals \p that.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst linear in `degree()`.
    [[nodiscard]] bool operator==(Blocks const& that) const {
      return _blocks == that._blocks && _lookup == that._lookup;
    }

    //! \brief Compare two blocks objects for inequality.
    //!
    //! Two Blocks objects are equal if and only if their underlying signed
    //! partitions are equal. It is ok to compare blocks of different
    //! degree with this operator.
    //!
    //! \param that the Blocks instance for comparison.
    //!
    //! \returns Whether or not  `*this` and \p that are not equal.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst linear in `degree()`.
    [[nodiscard]] bool operator!=(Blocks const& that) const {
      return !(*this == that);
    }

    //! \brief Compare two blocks objects for less.
    //!
    //! This operator defines a total order on the set of all Blocks objects
    //! (including those of different degree).
    //!
    //! \param that the Blocks instance for comparison.
    //!
    //! \returns Whether or not `*this` is less than \p that.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Linear in `degree()`.
    [[nodiscard]] bool operator<(Blocks const& that) const;

    // TODO(2) operator<=, operator>, operator>=

    //! \brief Return the degree of a blocks object.
    //!
    //! The *degree* of a Blocks object is the size of the set of
    //! which it is a partition, or the size of the \p blocks parameter to
    //! Blocks::Blocks.
    //!
    //! \returns The degree of a Blocks object.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] uint32_t degree() const noexcept {
      return _blocks.size();
    }

    //! \brief Return the number of blocks in a Blocks object.
    //!
    //! This function returns the number of parts in the partition that
    //! instances of this class represent.
    //!
    //! \returns The number of blocks.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! At worst \f$O(n)\f$ where \f$n\f$ is degree().
    [[nodiscard]] uint32_t number_of_blocks() const noexcept {
      return _lookup.size();
    }

    //! \brief Return the number of transverse blocks.
    //!
    //! This function returns the number of \c true values in
    //! lookup().
    //!
    //! \returns The number of signed (transverse) blocks in \c this.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At most linear in the number of blocks.
    [[nodiscard]] uint32_t rank() const;

    //! \brief Return a hash value for a Blocks instance.
    //!
    //! This value is recomputed every time this function is called.
    //!
    //! \returns A hash value for a \c this.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Linear in `degree()`.
    [[nodiscard]] size_t hash_value() const noexcept;

    //! \brief Return a const iterator pointing to the first transverse
    //! block lookup.
    //!
    //! The value pointed to is \c true if the \c i th block of \c this is a
    //! transverse block; and \c false otherwise.
    //!
    //! \returns A \ref lookup_const_iterator.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] lookup_const_iterator cbegin_lookup() const noexcept {
      return _lookup.cbegin();
    }

    //! \brief Return a const iterator pointing to the first transverse
    //! block lookup.
    //!
    //! \sa cbegin_lookup.
    [[nodiscard]] lookup_const_iterator cend_lookup() const noexcept {
      return _lookup.cend();
    }

    //! \brief Return a const reference to the transverse blocks lookup.
    //!
    //! The value in position \c i of the returned vector is \c true if the
    //! block with index \c i is transverse and \c false if it is not
    //! transverse.
    //!
    //! \returns A const reference to a `std::vector<bool>`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] std::vector<bool> const& lookup() const noexcept {
      return _lookup;
    }

    //! \brief Return a const iterator pointing to the index of the first
    //! block.
    //!
    //! \returns A value of type \ref const_iterator.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] const_iterator cbegin() const noexcept {
      return _blocks.cbegin();
    }

    //! \brief Return a const iterator pointing one past-the-end of the last
    //! block.
    //!
    //! \returns A value of type \ref const_iterator
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] const_iterator cend() const noexcept {
      return _blocks.cend();
    }

    //! \brief Return a const reference to the index of the block containing a
    //! point.
    //!
    //! \param i the point.
    //!
    //! \returns A value const reference to a value of type \c uint32_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] uint32_t const& operator[](size_t i) const {
      LIBSEMIGROUPS_ASSERT(i < _blocks.size());
      return _blocks[i];
    }

    //! \brief Return a const reference to the index of the block containing a
    //! point.
    //!
    //! \param i the point.
    //!
    //! \returns A value const reference to a value of type \c uint32_t.
    //!
    //! \throws std::out_of_range if \p i is out of range.
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] uint32_t const& at(size_t i) const {
      // TODO(2) better error
      return _blocks.at(i);
    }

   private:
    void throw_if_class_index_out_of_range(size_t index) const;
  };  // class Blocks

  //! \relates Blocks
  //!
  //! \brief Validate the arguments, construct a Blocks object, and validate
  //! it.
  //!
  //! \tparam Container the type of the parameter \p cont.
  //! \tparam Return the return type. Must satisfy
  //! `std::is_same<Return, Blocks>`.
  //!
  //! \param cont container containing a lookup for the blocks.
  //!
  //! \throws LibsemigroupsException if the arguments do not describe a
  //! signed partition.
  //!
  //! \throws LibsemigroupsException if the constructed Blocks object is not
  //! valid.
  template <typename Return, typename Container>
  [[nodiscard]] enable_if_is_same<Return, Blocks> make(Container const& cont) {
    detail::validate_args<Blocks>(cont);
    Blocks result(cont);
    blocks::validate(result);
    return result;
  }

  //! \relates Blocks
  //!
  //! \brief Validate the arguments, construct a Blocks object, and validate
  //! it.
  //!
  //! \tparam Return the return type. Must satisfy
  //! `std::is_same<Return, Blocks>`.
  //!
  //! \param cont container containing a lookup for the blocks.
  //!
  //! \throws LibsemigroupsException if the arguments do not describe a
  //! signed partition.
  //!
  //! \throws LibsemigroupsException if the constructed Blocks object is not
  //! valid.
  template <typename Return>
  [[nodiscard]] enable_if_is_same<Return, Blocks>
  make(std::initializer_list<std::vector<int32_t>> const& cont) {
    return make<Blocks, std::initializer_list<std::vector<int32_t>>>(cont);
  }

  //! \brief Return a human readable representation of a blocks object.
  //!
  //! Return a human readable representation (std::string) of a Blocks object.
  //! The braces to be used in the returns string can be specified using the
  //! argument \p braces. By default the returned string can be used to
  //! reconstruct the bipartition \p x. If the width of this returned string
  //! would be greater than \p max_width, then an abbreviated string is returned
  //! instead.
  //!
  //! \param x the Blocks object.
  //! \param braces the braces to use in the returned string (default: `"{}"`).
  //! \param max_width the maximum width of the returned string (default: \c
  //! 72).
  //!
  //! \returns A std::string representation of \p x.
  //!
  //! \throws LibsemigroupsException if \p braces does not have length \c 2.
  [[nodiscard]] std::string to_human_readable_repr(Blocks const&    x,
                                                   std::string_view braces
                                                   = "{}",
                                                   size_t max_width = 72);

  //! \ingroup bipart_group
  //!
  //! \brief Class for representing bipartitions.
  //!
  //! Defined in `bipart.hpp`.
  //!
  //! A *bipartition* is a partition of the set \f$\{0, ..., 2n - 1\}\f$ for
  //! some non-negative integer \f$n\f$; see the [Semigroups package for GAP
  //! documentation](https://semigroups.github.io/Semigroups/doc/chap3_mj.html)
  //! for more details.  The Bipartition class is more complex (i.e. has more
  //! member functions) than are used in `libsemigroups` because they are
  //! used in the GAP package [Semigroups package for
  //! GAP](https://semigroups.github.io/Semigroups/).
  //!
  //! \sa bipartition::validate(Bipartition const&).
  // TODO(2) add more explanation to the doc here
  class Bipartition {
   private:
    mutable size_t            _nr_blocks;
    mutable size_t            _nr_left_blocks;
    mutable std::vector<bool> _trans_blocks_lookup;
    mutable size_t            _rank;
    std::vector<uint32_t>     _vector;

   public:
    //! \brief Type of iterators pointing to block lookup.
    //!
    //! Type for iterators pointing to the lookup for the blocks of a
    //! bipartition.
    using iterator = std::vector<uint32_t>::iterator;

    //! \brief Type of const iterators pointing to block lookup.
    //!
    //! Type for const iterators pointing to the lookup for the blocks of a
    //! bipartition.
    using const_iterator = std::vector<uint32_t>::const_iterator;

    //! \brief Type of const iterators pointing to transverse blocks lookup.
    //!
    //! Type for  iterators pointing to the lookup for transverse blocks of a
    //! bipartition.
    using lookup_const_iterator = typename std::vector<bool>::const_iterator;

    //! \brief Construct an uninitialised bipartition of degree \c 0.
    //!
    //! Constructs an uninitialised bipartition of degree \c 0.
    Bipartition();

    //! \brief Construct an uninitialised bipartition of given degree.
    //!
    //! Constructs a bipartition of degree \p N that is uninitialised.
    //!
    //! \param N the degree of the bipartition.
    explicit Bipartition(size_t N);

    //! \brief Construct a bipartition from a const reference to blocks lookup.
    //!
    //! The parameter \p blocks:
    //! * is copied;
    //! * must have length \f$2n\f$ for some positive integer \f$n\f$;
    //! * consist of non-negative integers; and
    //! * have the property that if \f$i\f$, \f$i > 0\f$ occurs in \p blocks,
    //! then \f$i - 1\f$ occurs earlier in \p blocks.  The value of `blocks[i]`
    //! should represent the index of the block containing \c i.
    //!
    //! None of these conditions are verified.
    //!
    //! For example, if \p blocks is `{0, 1, 1, 2, 1, 1, 3, 1, 1, 4, 5, 6}`,
    //! then the above conditions are satisfied, but if \p blocks is `{1, 0, 1,
    //! 10}` then they are not.
    //!
    //! \param blocks a lookup for the blocks of the bipartition being
    //! constructed.
    //!
    //! \sa libsemigroups::validate(Bipartition const&).
    explicit Bipartition(std::vector<uint32_t> const& blocks);

    //! \brief Construct a bipartition from an rvalue reference to blocks
    //! lookup.
    //!
    //! \param blocks a lookup for the blocks of the bipartition being
    //! constructed.
    //!
    //! \sa Bipartition(std::vector<uint32_t> const&)
    //!  and libsemigroups::validate(Bipartition const&).
    explicit Bipartition(std::vector<uint32_t>&& blocks);

    //! \brief Construct a bipartition from an initializer list blocks lookup.
    //!
    //! \param blocks a lookup for the blocks of the bipartition being
    //! constructed.
    //!
    //! \sa Bipartition(std::vector<uint32_t> const&)
    //!  and libsemigroups::validate(Bipartition const&).
    Bipartition(std::initializer_list<uint32_t> const& blocks);

    //! \brief Construct a bipartition from a partition.
    //!
    //! The items in \p blocks should be:
    //! * duplicate-free;
    //! * pairwise disjoint; and
    //! * partition the set \f$\{-n, \ldots, -1, 1, \ldots, n\}\f$
    //! for some positive integer \f$n\f$.
    //!
    //! \param blocks the partition.
    //!
    //! \warning None of these conditions is checked by the constructor.
    //!
    //! \sa libsemigroups::validate(Bipartition const&).
    Bipartition(std::initializer_list<std::vector<int32_t>> const& blocks);

    //! \copydoc Bipartition(std::initializer_list<std::vector<int32_t>>
    //! const&)
    explicit Bipartition(std::vector<std::vector<int32_t>> const& blocks);

    //! \brief Default copy constructor.
    //!
    //! Default copy constructor.
    Bipartition(Bipartition const&);

    //! \brief Default move constructor.
    //!
    //! Default move constructor.
    Bipartition(Bipartition&&);

    //! \brief Default copy assignment operator.
    //!
    //! Default copy assignment operator.
    Bipartition& operator=(Bipartition const&);

    //! \brief Default move assignment operator.
    //!
    //! Default move assignment operator.
    Bipartition& operator=(Bipartition&&);

    ~Bipartition();

    //! \brief Compare bipartitions for equality.
    //!
    //! \param that the Bipartition for comparison.
    //!
    //! \returns Whether or not `*this` equals \p that.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst linear in degree().
    [[nodiscard]] bool operator==(Bipartition const& that) const {
      return _vector == that._vector;
    }

    //! \brief Compare bipartitions for less.
    //!
    //! This operator defines a total order on Bipartition objects. It is ok
    //! to compare Bipartition objects of different degrees.
    //!
    //! \param that the Bipartition for comparison.
    //!
    //! \returns Whether or not `*this` is less than \p that.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst linear in degree().
    [[nodiscard]] bool operator<(Bipartition const& that) const {
      return _vector < that._vector;
    }

    // TODO(2) other operators <=, >, >=, !=

    //! \brief Return a hash value.
    //!
    //! Returns a hash value for a Bipartition object.
    //!
    //! \returns
    //! A value of \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Linear in degree().
    // not noexcept because Hash<T>::operator() isn't
    [[nodiscard]] size_t hash_value() const {
      return Hash<std::vector<uint32_t>>()(_vector);
    }

    //! \brief Return the index of the block containing a value.
    //!
    //! Returns the index of the block containing a value.
    //! No bound checks are performed on the parameter \p i.
    //!
    //! \param i an integer
    //!
    //! \returns A reference to the index of the block containing \p i.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] uint32_t& operator[](size_t i) {
      return _vector[i];
    }

    //! \brief Return the index of the block containing a value.
    //!
    //! Returns the index of the block containing a value.
    //!
    //! No bound checks are performed on the parameter \p i.
    //!
    //! \param i an integer
    //!
    //! \returns A const reference to the index of the block containing \p i.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] uint32_t const& operator[](size_t i) const {
      return _vector[i];
    }

    //! \brief Return a reference to the index of the block containing a value.
    //!
    //! Returns a reference to the index of the block containing a value.
    //!
    //! \param i an integer
    //!
    //! \returns A reference to the index of the block containing \p i.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \throws std::out_of_range if the parameter \p i is out of range.
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] uint32_t& at(size_t i) {
      return _vector.at(i);
    }

    //! \brief Return a const reference to the index of the block containing a
    //! value.
    //!
    //! Returns a const reference to the index of the block containing a
    //! value.
    //!
    //! \param i an integer
    //!
    //! \returns A const reference to the index of the block containing \p i.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \throws std::out_of_range if the parameter \p i is out of range.
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] uint32_t const& at(size_t i) const {
      return _vector.at(i);
    }

    //! \brief Return a const iterator pointing to the index of the first
    //! block.
    //!
    //! Returns a const iterator pointing to the index of the first block.
    //!
    //! \returns A value of type \ref const_iterator.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] const_iterator cbegin() const noexcept {
      return _vector.cbegin();
    }

    //! \brief Return a const iterator pointing one beyond the last index of
    //! the last block.
    //!
    //! Returns a const iterator pointing one beyond the last index of the last
    //! block.
    //!
    //! \returns A value of type \ref const_iterator.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] const_iterator cend() const noexcept {
      return _vector.cend();
    }

    //! \brief Const iterator pointing to the index of the first left
    //! block.
    //!
    //! Returns a const iterator pointing to the index of the first left
    //! block.
    //!
    //! \returns A value of type \ref const_iterator.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] const_iterator cbegin_left_blocks() const noexcept {
      return cbegin();
    }

    //! \brief Const iterator pointing one beyond the last index of the
    //! last left block.
    //!
    //! Returns a const iterator pointing one beyond the last index of the
    //! last left block.
    //!
    //! \returns A value of type \ref const_iterator.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] const_iterator cend_left_blocks() const noexcept {
      return cbegin() + degree();
    }

    //! \brief Const iterator pointing to the index of the first right
    //! block.
    //!
    //! Returns a const iterator pointing to the index of the first right
    //! block.
    //!
    //! \returns A value of type \ref const_iterator.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] const_iterator cbegin_right_blocks() const noexcept {
      return cend_left_blocks();
    }

    //! \brief Const iterator pointing one beyond the last index of the
    //! last right block.
    //!
    //! Returns a const iterator pointing one beyond the last index of the
    //! last right block.
    //!
    //! \returns A value of type \ref const_iterator.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] const_iterator cend_right_blocks() const noexcept {
      return cend();
    }

    //! \brief Return the degree of the bipartition.
    //!
    //! Returns the degree of the bipartition.
    //!
    //! A bipartition is of degree \f$n\f$ if it is a partition of
    //! \f$\{0, \ldots, 2n -  1\}\f$.
    //!
    //! \returns A value of type `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] size_t degree() const noexcept;

    //! \brief Return an identity bipartition of given degree.
    //!
    //! Returns an identity bipartition of degree \p n.
    //!
    //! The *identity bipartition* of degree \f$n\f$ has blocks \f$\{i, -i\}\f$
    //! for all \f$i\in \{0, \ldots, n - 1\}\f$. This member function returns a
    //! new identity bipartition of degree equal to \p n.
    //!
    //! \param n the degree of the identity to be returned.
    //!
    //! \returns A newly constructed Bipartition.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    [[nodiscard]] static Bipartition one(size_t n);

    //! \brief Modify the current bipartition in-place to contain the product of
    //! two bipartitions.
    //!
    //! The parameter \p thread_id can be used some temporary storage is
    //! required to find the product of \p x and \p y.
    //!
    //! \param x the first bipartition to multiply
    //! \param y the second bipartition to multiply
    //! \param thread_id the index of the calling thread (defaults to \c 0)
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Quadratic in `x.degree()`.
    //!
    //! \warning
    //! If different threads call this function concurrently with the same
    //! parameter \p thread_id, then bad things will happen.
    //!
    //! \warning This function expects its arguments to have equal degree, but
    //! this is not checked.
    void product_inplace_no_checks(Bipartition const& x,
                                   Bipartition const& y,
                                   size_t             thread_id = 0);

    //! \brief Return the number of transverse blocks.
    //!
    //! The *rank* of a bipartition is the number of blocks containing both
    //! positive and negative values, which are referred to as the *transverse*
    //! blocks.
    //!
    //! \returns The number of transverse blocks.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(n)\f$ where \f$n\f$ is degree().
    [[nodiscard]] size_t rank() const;

    //! \brief Return the number of blocks in a Bipartition.
    //!
    //! This function returns the number of parts in the partition that
    //! instances of this class represent.
    //!
    //! \returns The number of blocks.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst \f$O(n)\f$ where \f$n\f$ is the degree().
    [[nodiscard]] uint32_t number_of_blocks() const;

    //! \brief Return the number of blocks containing a positive integer.
    //!
    //! The *left blocks* of a bipartition is the partition of
    //! \f$\{0, \ldots, n - 1\}\f$ induced by the bipartition. This member
    //! function returns the number of blocks in this partition.
    //!
    //! \returns
    //! A value of type \c uint32_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst \f$O(n)\f$ where \f$n\f$ is the degree().
    [[nodiscard]] uint32_t number_of_left_blocks() const;

    //! \brief Return the number of blocks containing a negative integer.
    //!
    //! The *right blocks* of a bipartition is the partition of
    //! \f$\{n, \ldots, 2n - 1\}\f$ induced by the bipartition. This member
    //! function returns the number of blocks in this partition.
    //!
    //! \returns
    //! A value of type \c uint32_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst \f$O(n)\f$ where \f$n\f$ is the degree().
    [[nodiscard]] uint32_t number_of_right_blocks() const;

    //! \brief Check if a block is a transverse block.
    //!
    //! A block of a biparition is *transverse* if it contains integers less
    //! than and greater than \f$n\f$, which is the degree of the bipartition.
    //!
    //! \param index the index of a block
    //!
    //! \returns Whether or not the given block is transverse.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst \f$O(n)\f$ where \f$n\f$ is the degree().
    //!
    //! \warning This function does no checks on its arguments.
    [[nodiscard]] bool is_transverse_block_no_checks(size_t index) const;

    //! \brief Check if a block is a transverse block.
    //!
    //! A block of a biparition is *transverse* if it contains integers less
    //! than and greater than \f$n\f$, which is the degree of the bipartition.
    //!
    //! \param index the index of a block.
    //!
    //! \returns Whether or not the given block is transverse.
    //!
    //! \throws LibsemigroupsException if \p index is not in the range from \c
    //! 0 to \ref number_of_left_blocks.
    //!
    //! \complexity
    //! At worst \f$O(n)\f$ where \f$n\f$ is the degree().
    [[nodiscard]] bool is_transverse_block(size_t index) const;

    //! \brief Return a pointer to the left blocks of a bipartition.
    //!
    //! The *left blocks* of a bipartition is the partition of
    //! \f$\{0, \ldots, n - 1\}\f$ induced by the bipartition. This member
    //! function returns a Blocks object representing this partition.
    //!
    //! \returns
    //! A pointer to a newly constructed Blocks object.
    //!
    //! \throws LibsemigroupsException if \c this is not valid.
    //!
    //! \complexity
    //! \f$O(n)\f$ where \f$n\f$ is the degree().
    // TODO(2) remove this
    [[nodiscard]] Blocks* left_blocks() const;

    //! \brief Return a pointer to the left blocks of a bipartition.
    //!
    //! The *left blocks* of a bipartition is the partition of
    //! \f$\{0, \ldots, n - 1\}\f$ induced by the bipartition. This member
    //! function returns a Blocks object representing this partition.
    //!
    //! \returns
    //! A pointer to a newly constructed Blocks object.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(n)\f$ where \f$n\f$ is the degree().
    [[nodiscard]] Blocks* left_blocks_no_checks() const;

    //! \brief Return a pointer to the right blocks of a bipartition.
    //!
    //! The *right blocks* of a bipartition is the partition of
    //! \f$\{n, \ldots, 2n - 1\}\f$ induced by the bipartition.
    //!
    //! \returns
    //! A pointer to a newly constructed Blocks object.
    //!
    //! \throws LibsemigroupsException if \c this is not valid.
    //!
    //! \complexity
    //! \f$O(n)\f$ where \f$n\f$ is the degree().
    // TODO(2) remove this
    [[nodiscard]] Blocks* right_blocks() const;

    //! \brief Return a pointer to the right blocks of a bipartition.
    //!
    //! The *right blocks* of a bipartition is the partition of
    //! \f$\{n, \ldots, 2n - 1\}\f$ induced by the bipartition.
    //!
    //! \returns
    //! A pointer to a newly constructed Blocks object.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(n)\f$ where \f$n\f$ is the degree().
    [[nodiscard]] Blocks* right_blocks_no_checks() const;

    //! \brief Set the number of blocks.
    //!
    //! This function sets the number of blocks of \c this to \p n. No checks
    //! are performed.
    //!
    //! \param n the number of blocks.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    void set_number_of_blocks(size_t n) noexcept;

    //! \brief Set the number of left blocks.
    //!
    //! This function sets the number of left blocks of \c this to \p n. No
    //! checks are performed.
    //!
    //! \param n the number of blocks.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    void set_number_of_left_blocks(size_t n) noexcept;

    //! \brief Set the rank.
    //!
    //! This function sets the \c rank of \c this to \p n. No
    //! checks are performed.
    //!
    //! \param n the rank.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    void set_rank(size_t n) noexcept;

    //! \brief Return a const iterator pointing to the first transverse
    //! block lookup.
    //!
    //! The value pointed to is \c true if the \c i th block of \c this is a
    //! transverse block; and \c false otherwise.
    //!
    //! \returns A \ref lookup_const_iterator.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] lookup_const_iterator cbegin_lookup() const noexcept {
      init_trans_blocks_lookup();
      return _trans_blocks_lookup.cbegin();
    }

    //! \brief Return a const iterator pointing to the first transverse
    //! block lookup.
    //!
    //! \sa cbegin_lookup.
    [[nodiscard]] lookup_const_iterator cend_lookup() const noexcept {
      init_trans_blocks_lookup();
      return _trans_blocks_lookup.cend();
    }

    //! \brief Return a const reference to the transverse blocks lookup.
    //!
    //! The value in position \c i of the returned vector is \c true if the
    //! block with index \c i is transverse and \c false if it is not
    //! transverse.
    //!
    //! \returns A const reference to a `std::vector<bool>`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] std::vector<bool> const& lookup() const noexcept {
      return _trans_blocks_lookup;
    }

   private:
    void init_trans_blocks_lookup() const;
  };  // class Bipartition

  //! \relates Bipartition
  //!
  //! \brief Validate the arguments, construct a bipartition, and validate it.
  //!
  //! \tparam Container the type of the parameter \p cont.
  //! \tparam Return the return type. Must satisfy
  //! `std::is_same<Return, Bipartition>`.
  //!
  //! \param cont either a vector providing a lookup for the blocks of the
  //! bipartition or a vector of vectors (or initializer list).
  //!
  //! \throws LibsemigroupsException if the arguments do not describe a
  //! bipartition.
  //!
  //! \throws LibsemigroupsException if the constructed bipartition is not
  //! valid.
  template <typename Return, typename Container>
  [[nodiscard]] enable_if_is_same<Return, Bipartition>
  make(Container const& cont) {
    detail::validate_args<Bipartition>(cont);
    Bipartition result(cont);
    bipartition::validate(result);
    return result;
  }

  //! \relates Bipartition
  //!
  //! \brief Validate the arguments, construct a bipartition, and validate it.
  //!
  //! \copydoc make(Container const&)
  template <typename Return>
  [[nodiscard]] enable_if_is_same<Return, Bipartition>
  make(std::initializer_list<uint32_t> const& cont) {
    return make<Bipartition, std::initializer_list<uint32_t>>(cont);
  }

  //! \relates Bipartition
  //!
  //! \copydoc make(Container const&)
  template <typename Return>
  [[nodiscard]] enable_if_is_same<Return, Bipartition>
  make(std::initializer_list<std::vector<int32_t>> const& cont) {
    return make<Bipartition, std::initializer_list<std::vector<int32_t>>>(cont);
  }

  //! \ingroup bipart_group
  //!
  //! \brief Return a human readable representation of a bipartition.
  //!
  //! Return a human readable representation (std::string) of a bipartition.
  //! The braces to be used in the returns string can be specified using the
  //! argument \p braces. By default the returned string can be used to
  //! reconstruct the bipartition \p x. If the width of this returned string
  //! would be greater than \p max_width, then an abbreviated string is returned
  //! instead.
  //!
  //! \param x the Bipartition object.
  //! \param braces the braces to use in the returned string (default: `"{}"`).
  //! \param max_width the maximum width of the returned string (default: \c
  //! 72).
  //!
  //! \returns A std::string representation of \p x.
  //!
  //! \throws LibsemigroupsException if \p braces does not have length \c 2.
  [[nodiscard]] std::string to_human_readable_repr(Bipartition const& x,
                                                   std::string_view   braces
                                                   = "{}",
                                                   size_t max_width = 72);

  //! \ingroup bipart_group
  //!
  //! \brief Multiply two bipartitions.
  //!
  //! Returns a newly constructed bipartition equal to the product of \p x and
  //! \p y.
  //!
  //! \param x a bipartition
  //! \param y a bipartition
  //!
  //! \returns
  //! A value of type \c Bipartition
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  //!
  //! \complexity
  //! Quadratic in `x.degree()`.
  //!
  //! \warning This function expects its arguments to have equal degree, but
  //! this is not checked.
  [[nodiscard]] Bipartition operator*(Bipartition const& x,
                                      Bipartition const& y);

  //! \ingroup bipart_group
  //!
  //! \brief Check bipartitions for inequality.
  //!
  //! \param x a bipartition
  //! \param y a bipartition
  //!
  //! \returns
  //! A value of type \c bool.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  //!
  //! \complexity
  //! At worst linear in the degree of \p x and \p y.
  [[nodiscard]] inline bool operator!=(Bipartition const& x,
                                       Bipartition const& y) {
    return !(x == y);
  }

  //! \ingroup bipart_group
  //!
  //! \brief Compare bipartitions.
  //!
  //! Convenience function that just calls `operator<` and `operator==`.
  [[nodiscard]] inline bool operator<=(Bipartition const& x,
                                       Bipartition const& y) {
    return x < y || x == y;
  }

  //! \ingroup bipart_group
  //!
  //! \brief Compare bipartitions.
  //!
  //! Convenience function that just calls `operator<` with the arguments
  //! switched.
  [[nodiscard]] inline bool operator>(Bipartition const& x,
                                      Bipartition const& y) {
    return y < x;
  }

  //! \ingroup bipart_group
  //!
  //! \brief Compare bipartitions.
  //!
  //! Convenience function that just calls `operator<=` with the arguments
  //! switched.
  [[nodiscard]] inline bool operator>=(Bipartition const& x,
                                       Bipartition const& y) {
    return y <= x;
  }

  ////////////////////////////////////////////////////////////////////////
  // Adapters
  ////////////////////////////////////////////////////////////////////////

  //! Returns the approximate time complexity of multiplication.
  //!
  //! In the case of a Bipartition of degree *n* the value *2n ^ 2* is
  //! returned.
  template <>
  struct Complexity<Bipartition> {
    //! Call operator.
    //!
    //! \param x a const reference to a bipartition.
    //!
    //! \returns
    //! A value of type `size_t` representing the complexity of multiplying the
    //! parameter \p x by another bipartition of the same degree.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] size_t operator()(Bipartition const& x) const noexcept {
      return x.degree() * x.degree();
    }
  };

  template <>
  struct Degree<Bipartition> {
    [[nodiscard]] size_t operator()(Bipartition const& x) const noexcept {
      return x.degree();
    }
  };

  template <>
  struct Hash<Bipartition> {
    [[nodiscard]] size_t operator()(Bipartition const& x) const {
      return x.hash_value();
    }
  };

  template <>
  struct One<Bipartition> {
    [[nodiscard]] Bipartition operator()(Bipartition const& x) const {
      return (*this)(x.degree());
    }

    [[nodiscard]] Bipartition operator()(size_t N = 0) const {
      return Bipartition::one(N);
    }
  };

  template <>
  struct Product<Bipartition> {
    void operator()(Bipartition&       xy,
                    Bipartition const& x,
                    Bipartition const& y,
                    size_t             thread_id = 0) {
      xy.product_inplace_no_checks(x, y, thread_id);
    }
  };

  template <>
  struct IncreaseDegree<Bipartition> {
    void operator()(Bipartition&, size_t) {}
  };

}  // namespace libsemigroups
#endif  // LIBSEMIGROUPS_BIPART_HPP_
