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

// This file contains the declaration of the PBR class.

#ifndef LIBSEMIGROUPS_PBR_HPP_
#define LIBSEMIGROUPS_PBR_HPP_

#include <cstddef>           // for size_t
#include <cstdint>           // for uint32_t, int32_t
#include <initializer_list>  // for initializer_list
#include <iosfwd>            // for ostream, ostringstream
#include <string>            // for string
#include <type_traits>       // for forward
#include <utility>           // for forward
#include <vector>            // for vector, operator<, operator==, allocator

#include "adapters.hpp"  // for Hash
#include "types.hpp"     // for enable_if_is_same

namespace libsemigroups {

  //! \defgroup pbr_group Partitioned binary relations (PBRs)
  //!
  //! Defined `pbr.hpp`.
  //!
  //! This page contains an overview of the functionality in `libsemigroups`
  //! for partitioned binary relations (PBRs).
  //!
  //! PBRs are a generalisation of bipartitions, which were introduced by
  //! Martin and Mazorchuk in \cite Martin2011aa.
  //!
  //! Helper functions for PBRs are documented at
  //! \ref libsemigroups::pbr "Helper functions for PBRs"
  //!
  //! @{

  //! \brief Class for representing PBRs.
  //!
  //! Defined in `pbr.hpp`.
  //!
  //! *Partitioned binary relations* (PBRs) are a generalisation of a
  //! bipartitions, and were introduced by Martin and Mazorchuk in
  //! \cite Martin2011aa.
  //!
  //! \sa pbr::throw_if_invalid(PBR const&).
  class PBR {
   public:
    //! \brief Type of constructor argument.
    //!
    //! Type of constructor argument.
    template <typename T>
    using vector_type = std::vector<std::vector<T>> const&;

    //! \brief Type of constructor argument.
    //!
    //! Type of constructor argument.
    template <typename T>
    using initializer_list_type = std::initializer_list<std::vector<T>> const&;

    //! \brief Deleted.
    //!
    //! Deleted. To avoid the situation where the underlying container is not
    //! defined, it is not possible to default construct a PBR object.
    PBR() = delete;

    //! \brief Default copy constructor.
    //!
    //! Default copy constructor.
    PBR(PBR const&) = default;

    //! \brief Default move constructor.
    //!
    //! Default move constructor.
    PBR(PBR&&) = default;

    //! \brief Default copy assignment operator.
    //!
    //! Default copy assignment operator.
    PBR& operator=(PBR const&) = default;

    //! \brief Default move assignment operator.
    //!
    PBR& operator=(PBR&&) = default;

    //! \brief Construct from adjacencies \c 0 to `2n - 1`.
    //!
    //! Construct from adjacencies \c 0 to `2n - 1`.
    //!
    //! The parameter \p x must be a container of vectors of \c uint32_t with
    //! size \f$2n\f$ for some integer \f$n\f$, and the vector in position
    //! \f$i\f$ is the list of points adjacent to \f$i\f$ in the PBR
    //! constructed.
    //!
    //! \param  x the container of vectors of adjacencies.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning
    //! No checks whatsoever on the validity of \p x are performed.
    //!
    //! \sa \ref pbr::throw_if_invalid(PBR const&)
    explicit PBR(vector_type<uint32_t> x);

    //! \copydoc PBR(vector_type<uint32_t>)
    explicit PBR(initializer_list_type<uint32_t> x);

    //! \brief Construct empty PBR of given \ref degree.
    //!
    //! Construct empty PBR of given \ref degree.
    //!
    //! \param n the degree
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    explicit PBR(size_t n);

    //! \copydoc PBR(vector_type<int32_t>,vector_type<int32_t>)
    PBR(initializer_list_type<int32_t> left,
        initializer_list_type<int32_t> right);

    //! \brief Construct from adjacencies \c 1 to \c n and \c -1 to \c -n.
    //!
    //! Construct from adjacencies \c 1 to \c n and \c -1 to \c -n.
    //!
    //! The parameters \p left and \p right should be containers of
    //! \f$n\f$ vectors of integer values, so that
    //! the vector in position \f$i\f$ of \p left is the list of points
    //! adjacent to \f$i\f$ in the PBR, and the vector in position \f$i\f$
    //! of \p right is the list of points adjacent to \f$n + i\f$ in the PBR.
    //! A negative value \f$i\f$ corresponds to \f$n - i\f$.
    //!
    //! \param left container of adjacencies of \c 1 to \c n
    //! \param right container of adjacencies of \c n + 1 to \c 2n.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning
    //! No checks whatsoever on the validity of \p left or \p right are
    //! performed.
    //!
    //! \sa \ref pbr::throw_if_invalid(PBR const&) and
    //! \ref make(initializer_list_type<int32_t>,
    //! initializer_list_type<int32_t>)
    PBR(vector_type<int32_t> left, vector_type<int32_t> right);

    //! \brief Returns the degree of a PBR.
    //!
    //! Returns the degree of a PBR, where the *degree* of a PBR is half the
    //! number of points in the PBR.
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    size_t degree() const noexcept;

    //! \brief Returns the number of points of a PBR.
    //!
    //! Returns the number of points of a PBR.
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    size_t number_of_points() const noexcept;

    //! \brief Multiply two PBR objects and store the product in \c this.
    //!
    //! Replaces the contents of \c this by the product of \p x and \p y.
    //!
    //! The parameter \p thread_id is required since some temporary storage is
    //! required to find the product of \p x and \p y.  Note that if different
    //! threads call this member function with the same value of \p thread_id
    //! then bad things will happen.
    //!
    //! \param x a PBR.
    //! \param y a PBR.
    //! \param thread_id the index of the calling thread (defaults to \c 0).
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning
    //! No checks are made on whether or not the parameters are compatible. If
    //! \p x and \p y have different degrees, then bad things will happen.
    void product_inplace_no_checks(PBR const& x,
                                   PBR const& y,
                                   size_t     thread_id = 0);

    //! \brief Multiply two PBR objects and store the product in \c this.
    //!
    //! Replaces the contents of \c this by the product of \p x and \p y.
    //!
    //! The parameter \p thread_id is required since some temporary storage is
    //! required to find the product of \p x and \p y.  Note that if different
    //! threads call this member function with the same value of \p thread_id
    //! then bad things will happen.
    //!
    //! \param x a PBR.
    //! \param y a PBR.
    //! \param thread_id the index of the calling thread (defaults to \c 0).
    //!
    //! \throws LibsemigroupsException if:
    //!   * the \ref degree() of \p x is not the same as the \ref degree() of \p
    //!   y;
    //!   * the \ref degree() of `*this` is not the same as the \ref degree() of
    //!   \p x;
    //!   * one if the addresses `&x` and `&y` is the same as that of `this`; or
    //!   * either \p x or \p y are invalid.
    //!
    //! \sa
    //! \ref pbr::throw_if_invalid.
    void product_inplace(PBR const& x, PBR const& y, size_t thread_id = 0);

    //! \brief Compare two PBRs for equality.
    //!
    //! Compare two PBRs for equality.
    //!
    //! \param that a PBR to compare with.
    //!
    //! \returns \c true if \c this equals \p that, and \c false otherwise.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst linear in degree().
    // TODO(later): a better explanation of what equality means for PBRs
    bool operator==(PBR const& that) const {
      return _vector == that._vector;
    }

    //! \brief Compare for less
    //!
    //! Compare for less.
    //!
    //! \param that a PBR to compare with.
    //!
    //! \returns \c true if \c this is less than \p that, and \c false
    //! otherwise.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst linear in degree().
    bool operator<(PBR const& that) const {
      return _vector < that._vector;
    }

    //! \brief Returns a reference to the points adjacent to a given point.
    //!
    //! Returns a reference to the points adjacent to a given point.
    //!
    //! \param i the point.
    //!
    //! \returns A value reference to a `std::vector<uint32_t>`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    std::vector<uint32_t>& operator[](size_t i) {
      return _vector[i];
    }

    //! \brief Returns a const reference to the points adjacent to a given
    //! point.
    //!
    //! Returns a const reference to the points adjacent to a given point.
    //!
    //! \param i the point.
    //!
    //! \returns A value const reference to a `std::vector<uint32_t>`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    std::vector<uint32_t> const& operator[](size_t i) const {
      return _vector[i];
    }

    //! \brief Returns a reference to the points adjacent to a given point, with
    //! bounds checking.
    //!
    //! Returns a reference to the points adjacent to a given point, with bounds
    //! checking.
    //!
    //! \param i the point.
    //!
    //! \returns A value reference to a `std::vector<uint32_t>`.
    //!
    //! \throws std::out_of_range if \p i > \ref number_of_points().
    //!
    //! \complexity
    //! Constant.
    std::vector<uint32_t>& at(size_t i);

    //! \brief Returns a const reference to the points adjacent to a given
    //! point, with bounds checking.
    //!
    //! Returns a const reference to the points adjacent to a given point, with
    //! bounds checking.
    //!
    //! \param i the point.
    //!
    //! \returns A value const reference to a `std::vector<uint32_t>`.
    //!
    //! \throws std::out_of_range if \p i > \ref number_of_points().
    //!
    //! \complexity
    //! Constant.
    std::vector<uint32_t> const& at(size_t i) const;

    //! \brief Returns a hash value for a PBR.
    //!
    //! Returns a hash value for a PBR.
    //!
    //!
    //! \returns A hash value for a \c this.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Linear in `degree()`.
    //!
    //! \note
    //! This value is recomputed every time this function is called.
    // not noexcept because Hash<T>::operator() isn't
    size_t hash_value() const {
      return Hash<std::vector<std::vector<uint32_t>>>()(_vector);
    }

    //! \brief Insertion operator
    //!
    //! This member function allows PBR objects to be inserted into a
    //! std::ostringstream.
    friend std::ostringstream& operator<<(std::ostringstream&, PBR const&);

    //! \brief Insertion operator
    //!
    //! This member function allows PBR objects to be inserted into a
    //! std::ostream.
    friend std::ostream& operator<<(std::ostream&, PBR const&);

   private:
    std::vector<std::vector<uint32_t>> _vector;
  };

  //! \brief Namespace for PBR helper functions.
  //!
  //! This namespace contains helper functions for the PBR class.
  namespace pbr {
    //! \brief Returns the identity PBR with specified degree.
    //!
    //! This function returns a new PBR with degree equal to \p n where every
    //! value is adjacent to its negative. Equivalently, \f$i\f$ is adjacent
    //! \f$i + n\f$ and vice versa for every \f$i\f$ less than the degree
    //! \f$n\f$.
    //!
    //! \param n the degree.
    //!
    //! \returns
    //! A PBR.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \sa
    //! \ref one(PBR const&)
    PBR one(size_t n);

    //! \brief Returns the identity PBR with degree `x.degree()`.
    //!
    //! This member function returns a new \ref PBR with degree equal to the
    //! degree of \p x, where every value is adjacent to its negative.
    //! Equivalently, \f$i\f$ is adjacent \f$i + n\f$ and vice versa for every
    //! \f$i\f$ less than the degree \f$n\f$.
    //!
    //! \param x A PBR.
    //!
    //! \returns
    //! A PBR.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \sa
    //! \ref one(size_t)
    PBR one(PBR const& x);

    // TODO(later) analogue of bipartition::underlying_partition?

    //! \brief Throws if a PBR has an odd number of points.
    //!
    //! This function throws a LibsemigroupsException if the argument \p x does
    //! not describe a binary relation on an even number of points.
    //!
    //! \param x the PBR to validate.
    //!
    //! \throws LibsemigroupsException \p x does not describe a binary relation
    //! on an even number of points.
    //!
    //! \complexity
    //! Constant.
    void throw_if_not_even_length(PBR const& x);

    //! \brief Throws if a PBR has a point related to a point that is greater
    //! than degree().
    //!
    //! This function throws a LibsemigroupsException if the argument \p x has a
    //! point related to a point that is greater than \ref PBR::degree().
    //!
    //! \param x the PBR to validate.
    //!
    //! \throws LibsemigroupsException if \p x has a point related to a point
    //! that is greater than degree().
    //!
    //! \complexity
    //! Linear in the PBR::degree of \p x.
    void throw_if_entry_out_of_bounds(PBR const& x);

    //! \brief Throws if a PBR has a list of points related to a point that is
    //! not sorted.
    //!
    //! This function throws a LibsemigroupsException if the argument \p x has a
    //! list of points related to a point that is not sorted.
    //!
    //! \param x the PBR to validate.
    //!
    //! \throws LibsemigroupsException if \p x has a list of points related to a
    //! point that is not sorted.
    //!
    //! \complexity
    //! Linear in the PBR::degree of \p x.
    void throw_if_adjacencies_unsorted(PBR const& x);

    //! \brief Throws if a PBR is invalid.
    //!
    //! This function throws a LibsemigroupsException if the argument \p x is
    //! not a valid PBR.
    //!
    //! \param x the PBR to validate.
    //!
    //! \throws LibsemigroupsException if any of the following occur:
    //! * \p x does not describe a binary relation on an even number of points;
    //! * \p x has a point related to a point that is greater than degree();
    //! * a list of points related to a point is not sorted.
    //!
    //! \complexity
    //! Linear in the PBR::degree x.
    //!
    //! \sa
    //! * \ref throw_if_not_even_length(PBR const&);
    //! * \ref throw_if_entry_out_of_bounds(PBR const&);
    //! * \ref throw_if_adjacencies_unsorted(PBR const&).
    void inline throw_if_invalid(PBR const& x) {
      throw_if_not_even_length(x);
      throw_if_entry_out_of_bounds(x);
      throw_if_adjacencies_unsorted(x);
    }

  }  // namespace pbr

  // end pbr_group
  //! @}

  namespace detail {
    std::vector<std::vector<uint32_t>>
    process_left_right(PBR::vector_type<int32_t> left,
                       PBR::vector_type<int32_t> right);
  }

  //! \relates PBR
  //!
  //! \brief Construct and validate a \ref PBR.
  //!
  //! Construct and validate a \ref PBR.
  //!
  //! \tparam Return the return type. Must satisfy `std::is_same<Return, PBR>`.
  //! \tparam T the types of the arguments.
  //!
  //! \param args the arguments to forward to the \ref PBR constructor.
  //!
  //! \returns
  //! A PBR constructed from \p args and validated.
  //!
  //! \throws LibsemigroupsException if libsemigroups::throw_if_invalid(PBR
  //! const&) throws when called with the constructed PBR.
  //!
  //! \warning
  //! No checks are performed on the validity of \p args prior to the
  //! construction of the PBR object.
  //!
  //! \sa
  //! #PBR.
  template <typename Return, typename... T>
  [[nodiscard]] enable_if_is_same<Return, PBR> make(T... args) {
    // TODO(later) validate_args
    PBR result(std::forward<T>(args)...);
    pbr::throw_if_invalid(result);
    return result;
  }

  //! \relates PBR
  //!
  //! \brief Construct and validate a \ref PBR.
  //!
  //! Construct and validate a \ref PBR.
  //!
  //! \tparam Return the return type. Must satisfy `std::is_same<Return, PBR>`.
  //!
  //! \param args the arguments to forward to the constructor.
  //!
  //! \returns
  //! A PBR constructed from \p args and validated.
  //!
  //! \throws LibsemigroupsException if libsemigroups::throw_if_invalid(PBR
  //! const&) throws when called with the constructed PBR.
  //!
  //! \warning
  //! No checks are performed on the validity of \p args prior to the
  //! construction of the PBR object.
  //!
  //! \sa
  //! PBR(initializer_list_type<uint32_t>).
  template <typename Return>
  [[nodiscard]] enable_if_is_same<Return, PBR>
  make(PBR::initializer_list_type<uint32_t> args) {
    return make<PBR, decltype(args)>(args);
  }

  //! \relates PBR
  //!
  //! \copydoc make(PBR::vector_type<int32_t>, PBR::vector_type<int32_t>)
  template <typename Return>
  [[nodiscard]] enable_if_is_same<Return, PBR>
  make(PBR::initializer_list_type<int32_t> left,
       PBR::initializer_list_type<int32_t> right) {
    return PBR(detail::process_left_right(left, right));
  }

  //! \relates PBR
  //!
  //! \brief Construct and validate a \ref PBR.
  //!
  //! Construct and validate a \ref PBR.
  //!
  //! \tparam Return the return type. Must satisfy `std::is_same<Return, PBR>`.
  //!
  //! \param left the 1st argument to forward to the constructor.
  //! \param right the 2nd argument to forward to the constructor.
  //!
  //! \returns
  //! A PBR constructed from \p args and validated.
  //!
  //! \throws LibsemigroupsException if libsemigroups::throw_if_invalid(PBR
  //! const&) throws when called with the constructed PBR.
  //!
  //! \warning
  //! No checks are performed on the validity of \p left or \p right prior to
  //! the construction of the PBR object.
  //!
  //! \sa
  //! PBR(initializer_list_type<int32_t>, initializer_list_type<int32_t>).
  template <typename Return>
  [[nodiscard]] enable_if_is_same<Return, PBR>
  make(PBR::vector_type<int32_t> left, PBR::vector_type<int32_t> right) {
    return PBR(detail::process_left_right(left, right));
  }

  //! \relates PBR
  //!
  //! \brief Return a human readable representation of a PBR.
  //!
  //! Return a human readable representation of a PBR.
  //!
  //! \param x the PBR object.
  //!
  //! \returns
  //! A value of type std::string.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  [[nodiscard]] std::string to_human_readable_repr(PBR const& x);

  //! \relates PBR
  //!
  //! \brief Multiply two PBRs.
  //!
  //! Returns a newly constructed PBR equal to the product of \p x and \p y.
  //!
  //! \param x a PBR.
  //! \param y a PBR.
  //!
  //! \returns
  //! A value of type \ref PBR
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  //!
  //! \complexity
  //! Cubic in degree().
  PBR operator*(PBR const& x, PBR const& y);

  //! \relates PBR
  //!
  //! \brief Compare two PBRs for inequality.
  //!
  //! Compare two PBRs for inequality.
  //!
  //! \param x a PBR.
  //! \param y a PBR.
  //!
  //! \returns
  //! \returns \c true if \c this is not equal to \p that, and \c false
  //! otherwise.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  //!
  //! \complexity
  //! At worst quadratic in the degree of \p x and \p y.
  inline bool operator!=(PBR const& x, PBR const& y) {
    return !(x == y);
  }

  //! \relates PBR
  //!
  //! \brief Convenience function that just calls \ref PBR::operator<
  //! "operator<".
  //!
  //! Convenience function that just calls \ref PBR::operator< "operator<".
  inline bool operator>(PBR const& x, PBR const& y) {
    return y < x;
  }

  //! \relates PBR
  //!
  //! \brief Convenience function that just calls \ref PBR::operator<
  //! "operator<" and \ref PBR::operator== "operator==".
  //!
  //! Convenience function that just calls \ref PBR::operator< "operator<" and
  //! \ref PBR::operator== "operator==".
  inline bool operator<=(PBR const& x, PBR const& y) {
    return x < y || x == y;
  }

  //! \relates PBR
  //!
  //! \brief Convenience function that just calls \ref operator<=(PBR const&,
  //! PBR const&) "operator<=".
  //!
  //! Convenience function that just calls \ref operator<=(PBR const&,
  //! PBR const&) "operator<=".
  inline bool operator>=(PBR const& x, PBR const& y) {
    return y <= x;
  }

  namespace detail {

    template <typename T>
    struct IsPBRHelper : std::false_type {};

    template <>
    struct IsPBRHelper<PBR> : std::true_type {};

  }  // namespace detail

  //! \brief Helper variable template.
  //!
  //! The value of this variable is \c true if the template parameter \p T is
  //! \ref PBR.
  //!
  //! \tparam T a type.
  template <typename T>
  static constexpr bool IsPBR = detail::IsPBRHelper<T>::value;

  ////////////////////////////////////////////////////////////////////////
  // Adapters
  ////////////////////////////////////////////////////////////////////////

  //! \defgroup adapters_pbr_group Adapters for PBR
  //!
  //! This page contains links to the specific specialisations for some of the
  //! adapters on \ref adapters_group "this page" for the \ref PBR class.
  //!
  //! @{

  //! \brief Specialization of the adapter Complexity for instances of PBR.
  //!
  //! Specialization of the adapter Complexity for instances of PBR.
  //!
  //! \sa
  //! \ref Complexity.
  template <>
  struct Complexity<PBR> {
    //! \brief Returns the approximate time complexity of multiplying PBRs.
    //!
    //! Returns the approximate time complexity of multiplying PBRs, which is
    //! \f$2n ^ 3\f$ where \f$n\f$ is the degree.
    //!
    //! \param x a PBR.
    //!
    //! \returns
    //! A value of type `size_t` representing the complexity of multiplying the
    //! parameter \p x by another PBR of the same degree.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    size_t operator()(PBR const& x) const noexcept {
      return 8 * x.degree() * x.degree() * x.degree();
    }
  };

  //! \brief Specialization of the adapter Degree for instances of PBR.
  //!
  //! Specialization of the adapter Degree for instances of PBR.
  //!
  //! \sa
  //! \ref Degree.
  template <>
  struct Degree<PBR> {
    //! \brief Returns the degree of \p x.
    //!
    //! Returns the degree of \p x.
    //!
    //! \param x a PBR.
    //!
    //! \returns
    //! A value of type `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    //!
    //! \sa
    //! \ref PBR::degree
    size_t operator()(PBR const& x) const noexcept {
      return x.degree();
    }
  };

  //! \brief Specialization of the adapter Hash for instances of PBR.
  //!
  //! Specialization of the adapter Hash for instances of PBR.
  //!
  //! \sa
  //! \ref Hash.
  template <>
  struct Hash<PBR> {
    //! \brief Returns a hash value for \p x.
    //!
    //! Returns a hash value for \p x.
    //!
    //! \param x a PBR.
    //!
    //! \returns
    //! A value of `size_t`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Linear in `degree()`.
    //!
    //! \sa
    //! \ref PBR::hash_value
    size_t operator()(PBR const& x) const {
      return x.hash_value();
    }
  };

  //! \brief Specialization of the adapter One for instances of PBR.
  //!
  //! Specialization of the adapter One for instances of PBR.
  //!
  //! \sa
  //! \ref One.
  template <>
  struct One<PBR> {
    //! \brief Returns the identity PBR with degree `x.degree()`.
    //!
    //! This member function returns a new \ref PBR with degree equal to the
    //! degree of \p x, where every value is adjacent to its negative.
    //! Equivalently, \f$i\f$ is adjacent \f$i + n\f$ and vice versa for every
    //! \f$i\f$ less than the degree \f$n\f$.
    //!
    //! \returns
    //! A PBR.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \sa
    //! \ref pbr::one
    PBR operator()(PBR const& x) const {
      return pbr::one(x);
    }

    //! \brief Returns the identity PBR with specified degree.
    //!
    //! This function returns a new PBR with degree equal to \p N where every
    //! value is adjacent to its negative. Equivalently, \f$i\f$ is adjacent
    //! \f$i + N\f$ and vice versa for every \f$i\f$ less than the degree
    //! \f$N\f$.
    //!
    //! \param N the degree.
    //!
    //! \returns
    //! A PBR.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \sa
    //! \ref pbr::one
    PBR operator()(size_t N = 0) const {
      return pbr::one(N);
    }
  };

  //! \brief Specialization of the adapter Product for instances of PBR.
  //!
  //! Specialization of the adapter Product for instances of PBR.
  //!
  //! \sa
  //! \ref Product.
  template <>
  struct Product<PBR> {
    //! \brief Multiply two PBR objects and store the product in a third.
    //!
    //! Replaces the contents of \p xy by the product of \p x and \p y.
    //!
    //! The parameter \p thread_id is required since some temporary storage is
    //! required to find the product of \p x and \p y.  Note that if different
    //! threads call this member function with the same value of \p thread_id
    //! then bad things will happen.
    //!
    //! \param xy a PBR whose contents (if any) will be cleared.
    //! \param x a PBR.
    //! \param y a PBR.
    //! \param thread_id the index of the calling thread (defaults to \c 0).
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning
    //! No checks are made on whether or not the parameters are compatible. If
    //! \p x and \p y have different degrees, then bad things will happen.
    //!
    //! \sa
    //! \ref PBR::product_inplace_no_checks
    void operator()(PBR& xy, PBR const& x, PBR const& y, size_t thread_id = 0) {
      xy.product_inplace_no_checks(x, y, thread_id);
    }
  };

  //! \brief Specialization of the adapter IncreaseDegree for instances of PBR.
  //!
  //! Specialization of the adapter IncreaseDegree for instances of PBR.
  //!
  //! \sa
  //! \ref IncreaseDegree.
  template <>
  struct IncreaseDegree<PBR> {
    //! \brief Do nothing.
    void operator()(PBR&, size_t) {}
  };

  // end adapters_pbr_group
  //! @}

}  // namespace libsemigroups
#endif  // LIBSEMIGROUPS_PBR_HPP_
