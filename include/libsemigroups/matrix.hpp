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

// TODO(1) tpp file
// TODO(1) put the detail stuff into detail/matrix-common.hpp
// TODO(1) there're no complete set of init methods for matrices

#ifndef LIBSEMIGROUPS_MATRIX_HPP_
#define LIBSEMIGROUPS_MATRIX_HPP_

#include <algorithm>         // for min
#include <array>             // for array
#include <bitset>            // for bitset
#include <cstddef>           // for size_t
#include <cstdint>           // for uint64_t
#include <initializer_list>  // for initializer_list
#include <iosfwd>            // for ostringstream
#include <iterator>          // for distance
#include <numeric>           // for inner_product
#include <ostream>           // for operator<<, basic_ostream
#include <string>            // for string
#include <tuple>             // for tie
#include <type_traits>       // for false_type, is_signed, true_type
#include <unordered_map>     // for unordered_map
#include <unordered_set>     // for unordered_set
#include <utility>           // for forward, make_pair, pair
#include <vector>            // for vector

#include "adapters.hpp"   // for Degree
#include "bitset.hpp"     // for BitSet
#include "config.hpp"     // for PARSED_BY_DOXYGEN
#include "constants.hpp"  // for POSITIVE_INFINITY
#include "debug.hpp"      // for LIBSEMIGROUPS_ASSERT
#include "exception.hpp"  // for LIBSEMIGROUPS_EXCEPTION

#include "detail/containers.hpp"  // for StaticVector1
#include "detail/formatters.hpp"  // for formatter of POSITIVE_INFINITY ...
#include "detail/string.hpp"      // for detail::to_string

namespace libsemigroups {

  //! \defgroup matrix_group Matrices
  //!
  //! This page contains links to the documentation of the functionality in
  //! libsemigroups for matrices.
  //!
  //! There are several different classes in `libsemigroups` for representing
  //! matrices over various semirings. There are up to three different
  //! representations for every type of matrix:
  //!
  //! 1. those whose dimension and arithmetic operations can be defined at
  //!    compile time: \ref StaticMatrix;
  //! 2. those whose arithmetic operation can be defined at compile time but
  //!    whose dimensions can be set at run time: \ref
  //!    DynamicMatrixStaticArith
  //!    "DynamicMatrix (with compile-time arithmetic)";
  //! 3. those whose arithmetic operations and dimensions can be set at run
  //!    time: \ref DynamicMatrixDynamicArith
  //!    "DynamicMatrix (with run-time arithmetic)".
  //!
  //! It's unlikely that you will want to use the classes described on this
  //! page directly, but rather through the aliases described on the other
  //! matrix pages (such as, for example, \ref BMat).
  //!
  //! \warning All of the matrices in `libsemigroups` have entries that are
  //! integer types, and the code does not currently attempt to detect
  //! overflows.
  //!
  //! ## Over specific semirings
  //!
  //! The following matrix classes are provided which define matrices over some
  //! specific semirings:
  //!
  //! * \ref bmat_index_group
  //! * \ref intmat_group
  //! * \ref maxplusmat_group
  //! * \ref minplusmat_group
  //! * \ref maxplustruncmat_group
  //! * \ref minplustruncmat_group
  //! * \ref ntpmat_group
  //! * \ref projmaxplus_group
  //!
  //! ## Over arbitrary semirings
  //!
  //! The following general matrix classes are provided which can be used to
  //! define matrices over arbitrary semirings:
  //!
  //! * \ref DynamicMatrixStaticArith
  //! "DynamicMatrix (with compile-time arithmetic)"
  //! * \ref DynamicMatrixDynamicArith
  //! "DynamicMatrix (with run-time arithmetic)"
  //! * \ref StaticMatrix
  //!
  //! ## Row views
  //!
  //! A row view is a lightweight representation of a row of a matrix. The
  //! following row view classes are provided:
  //!
  //! * \ref DynamicRowViewStaticArith
  //! "DynamicRowView (compile-time arithmetic)
  //! * \ref DynamicRowViewDynamicArith
  //! "DynamicRowView (run-time arithmetic)
  //! * \ref StaticRowView
  //!
  //! ## Adapters
  //!
  //! There are various specialisations of the adapters described on
  //! \ref adapters_group "this page" for the matrices described on this page:
  //!
  //! * \ref adapters_matrix_group
  //! * \ref adapters_bmat_group
  //! * \ref adapters_bmat8_group

  ////////////////////////////////////////////////////////////////////////
  // Detail
  ////////////////////////////////////////////////////////////////////////

  namespace detail {

    template <typename T>
    struct IsStdBitSetHelper : std::false_type {};

    template <size_t N>
    struct IsStdBitSetHelper<std::bitset<N>> : std::true_type {};

    template <typename T>
    static constexpr bool IsStdBitSet = IsStdBitSetHelper<T>::value;

    struct MatrixPolymorphicBase {};

    template <typename T>
    struct IsMatrixHelper {
      static constexpr bool value
          = std::is_base_of<detail::MatrixPolymorphicBase, T>::value;
    };
  }  // namespace detail

  //! \ingroup matrix_group
  //!
  //! \brief Helper variable template.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This variable has value \c true if the template parameter \c T is
  //! derived from any of \ref DynamicMatrixStaticArith
  //! "DynamicMatrix (with compile-time arithmetic)",
  //! \ref DynamicMatrixDynamicArith
  //! "DynamicMatrix (with run-time arithmetic)",
  //! or StaticMatrix; and \c false otherwise.
  //!
  //! \tparam T the type to check.
  template <typename T>
  constexpr bool IsMatrix = detail::IsMatrixHelper<T>::value;

  namespace matrix {

    //! \brief Throws if a matrix is not square.
    //!
    //! This function throws a LibsemigroupsException if the matrix \p x is not
    //! square.
    //!
    //! \tparam Mat the type of the argument, must satisfy \ref IsMatrix<Mat>.
    //!
    //! \param x the matrix to check.
    //!
    //! \throws LibsemigroupsException if the number of rows in \p x does not
    //! equal the number of columns.
    template <typename Mat>
    auto throw_if_not_square(Mat const& x) -> std::enable_if_t<IsMatrix<Mat>> {
      if (x.number_of_rows() != x.number_of_cols()) {
        LIBSEMIGROUPS_EXCEPTION("expected a square matrix, but found {}x{}",
                                x.number_of_rows(),
                                x.number_of_cols());
      }
    }

    //! \brief Throws if two matrices do not have the same dimensions.
    //!
    //! This function throws a LibsemigroupsException if the matrices \p x and
    //! \p y do not have equal dimensions.
    //!
    //! \tparam Mat the type of the arguments, must satisfy \ref IsMatrix<Mat>.
    //!
    //! \param x the first matrix to check.
    //! \param y the second matrix to check.
    //!
    //! \throws LibsemigroupsException if the number of rows in \p x does not
    //! equal the number of rows of \p y; or the number of columns of \p x does
    //! not equal the number of columns of \p y.
    template <typename Mat>
    auto throw_if_bad_dim(Mat const& x, Mat const& y)
        -> std::enable_if_t<IsMatrix<Mat>> {
      if (x.number_of_rows() != y.number_of_rows()
          || x.number_of_cols() != y.number_of_cols()) {
        LIBSEMIGROUPS_EXCEPTION(
            "expected matrices with the same dimensions, the 1st argument is a "
            "{}x{} matrix, and the 2nd is a {}x{} matrix",
            x.number_of_rows(),
            x.number_of_cols(),
            y.number_of_rows(),
            y.number_of_cols());
      }
    }

    //! \brief Throws the arguments do not index an entry of a matrix.
    //!
    //! This function throws a LibsemigroupsException if \p r is not less than
    //! the number of rows of \p x; or if \p c is not less than the number of
    //! columns of \p x.
    //!
    //! \tparam Mat the type of the arguments, must satisfy \ref IsMatrix<Mat>.
    //!
    //! \param x the matrix.
    //! \param r the row index.
    //! \param c the column index.
    //!
    //! \throws LibsemigroupsException if `(r, c)` does not index an entry in
    //! the matrix \p x.
    template <typename Mat>
    auto throw_if_bad_coords(Mat const& x, size_t r, size_t c)
        -> std::enable_if_t<IsMatrix<Mat>> {
      if (r >= x.number_of_rows()) {
        LIBSEMIGROUPS_EXCEPTION("invalid row index in ({}, {}), expected "
                                "values in [0, {}) x [0, {})",
                                r,
                                c,
                                x.number_of_rows(),
                                x.number_of_cols(),
                                r);
      }
      if (c >= x.number_of_cols()) {
        LIBSEMIGROUPS_EXCEPTION("invalid column index in ({}, {}), expected "
                                "values in [0, {}) x [0, {})",
                                r,
                                c,
                                x.number_of_rows(),
                                x.number_of_cols(),
                                r);
      }
    }
  }  // namespace matrix

  ////////////////////////////////////////////////////////////////////////
  // Detail
  ////////////////////////////////////////////////////////////////////////
  namespace detail {
    template <typename Container,
              typename Subclass,
              typename TRowView,
              typename Semiring = void>
    class MatrixCommon : MatrixPolymorphicBase {
     public:
      ////////////////////////////////////////////////////////////////////////
      // MatrixCommon - Aliases - public
      ////////////////////////////////////////////////////////////////////////

      using scalar_type            = typename Container::value_type;
      using scalar_reference       = typename Container::reference;
      using scalar_const_reference = typename Container::const_reference;
      using semiring_type          = Semiring;

      using container_type = Container;
      using iterator       = typename Container::iterator;
      using const_iterator = typename Container::const_iterator;

      using RowView = TRowView;

      scalar_type scalar_one() const noexcept {
        return static_cast<Subclass const*>(this)->one_impl();
      }

      scalar_type scalar_zero() const noexcept {
        return static_cast<Subclass const*>(this)->zero_impl();
      }

      Semiring const* semiring() const noexcept {
        return static_cast<Subclass const*>(this)->semiring_impl();
      }

     private:
      ////////////////////////////////////////////////////////////////////////
      // MatrixCommon - Semiring arithmetic - private
      ////////////////////////////////////////////////////////////////////////

      scalar_type plus_no_checks(scalar_type x, scalar_type y) const noexcept {
        return static_cast<Subclass const*>(this)->plus_no_checks_impl(y, x);
      }

      scalar_type product_no_checks(scalar_type x,
                                    scalar_type y) const noexcept {
        return static_cast<Subclass const*>(this)->product_no_checks_impl(y, x);
      }

     protected:
      ////////////////////////////////////////////////////////////////////////
      // MatrixCommon - Container functions - protected
      ////////////////////////////////////////////////////////////////////////

      // TODO(1) use constexpr-if, not SFINAE
      template <typename SFINAE = container_type>
      auto resize(size_t r, size_t c) -> std::enable_if_t<
          std::is_same<SFINAE, std::vector<scalar_type>>::value> {
        _container.resize(r * c);
      }

      template <typename SFINAE = container_type>
      auto resize(size_t, size_t) -> std::enable_if_t<
          !std::is_same<SFINAE, std::vector<scalar_type>>::value> {}

     public:
      ////////////////////////////////////////////////////////////////////////
      // MatrixCommon - Constructors + destructor - public
      ////////////////////////////////////////////////////////////////////////

      // none of the constructors are noexcept because they allocate
      MatrixCommon()                               = default;
      MatrixCommon(MatrixCommon const&)            = default;
      MatrixCommon(MatrixCommon&&)                 = default;
      MatrixCommon& operator=(MatrixCommon const&) = default;
      MatrixCommon& operator=(MatrixCommon&&)      = default;

      explicit MatrixCommon(std::initializer_list<scalar_type> const& c)
          : MatrixCommon() {
        resize(1, c.size());
        std::copy(c.begin(), c.end(), _container.begin());
      }

      explicit MatrixCommon(std::vector<std::vector<scalar_type>> const& m)
          : MatrixCommon() {
        init(m);
      }

      MatrixCommon(
          std::initializer_list<std::initializer_list<scalar_type>> const& m)
          : MatrixCommon() {
        init(m);
      }

     private:
      // not noexcept because resize isn't
      template <typename T>
      void init(T const& m) {
        size_t const R = m.size();
        if (R == 0) {
          return;
        }
        size_t const C = m.begin()->size();
        resize(R, C);
        for (size_t r = 0; r < R; ++r) {
          auto row = m.begin() + r;
          for (size_t c = 0; c < C; ++c) {
            _container[r * C + c] = *(row->begin() + c);
          }
        }
      }

      // not noexcept because init isn't
      void
      init(std::initializer_list<std::initializer_list<scalar_type>> const& m) {
        init<std::initializer_list<std::initializer_list<scalar_type>>>(m);
      }

     public:
      explicit MatrixCommon(RowView const& rv) : MatrixCommon() {
        resize(1, rv.size());
        std::copy(rv.cbegin(), rv.cend(), _container.begin());
      }

      ~MatrixCommon() = default;

      // not noexcept because mem allocate is required
      Subclass one() const {
        size_t const n = number_of_rows();
        Subclass     x(semiring(), n, n);
        std::fill(x.begin(), x.end(), scalar_zero());
        for (size_t r = 0; r < n; ++r) {
          x(r, r) = scalar_one();
        }
        return x;
      }

      ////////////////////////////////////////////////////////////////////////
      // Comparison operators
      ////////////////////////////////////////////////////////////////////////

      // not noexcept because apparently vector::operator== isn't
      bool operator==(MatrixCommon const& that) const {
        return _container == that._container;
      }

      // not noexcept because apparently vector::operator== isn't
      bool operator==(RowView const& that) const {
        return number_of_rows() == 1
               && static_cast<RowView>(*static_cast<Subclass const*>(this))
                      == that;
      }

      // not noexcept because apparently vector::operator< isn't
      bool operator<(MatrixCommon const& that) const {
        return _container < that._container;
      }

      // not noexcept because apparently vector::operator< isn't
      bool operator<(RowView const& that) const {
        return number_of_rows() == 1
               && static_cast<RowView>(*static_cast<Subclass const*>(this))
                      < that;
      }

      // not noexcept because operator== isn't
      template <typename T>
      bool operator!=(T const& that) const {
        static_assert(IsMatrix<T> || std::is_same_v<T, RowView>);
        return !(*this == that);
      }

      // not noexcept because operator< isn't
      template <typename T>
      bool operator>(T const& that) const {
        static_assert(IsMatrix<T> || std::is_same_v<T, RowView>);
        return that < *this;
      }

      // not noexcept because operator< isn't
      template <typename T>
      bool operator>=(T const& that) const {
        static_assert(IsMatrix<T> || std::is_same_v<T, RowView>);
        return that < *this || that == *this;
      }

      // not noexcept because operator< isn't
      template <typename T>
      bool operator<=(T const& that) const {
        static_assert(IsMatrix<T> || std::is_same_v<T, RowView>);
        return *this < that || that == *this;
      }

      ////////////////////////////////////////////////////////////////////////
      // Attributes
      ////////////////////////////////////////////////////////////////////////

      // not noexcept because vector::operator[] isn't, and neither is
      // array::operator[]
      scalar_reference operator()(size_t r, size_t c) {
        return this->_container[r * number_of_cols() + c];
      }

      scalar_reference at(size_t r, size_t c) {
        matrix::throw_if_bad_coords(static_cast<Subclass const&>(*this), r, c);
        return this->operator()(r, c);
      }

      // not noexcept because vector::operator[] isn't, and neither is
      // array::operator[]
      scalar_const_reference operator()(size_t r, size_t c) const {
        return this->_container[r * number_of_cols() + c];
      }

      scalar_const_reference at(size_t r, size_t c) const {
        matrix::throw_if_bad_coords(static_cast<Subclass const&>(*this), r, c);
        return this->operator()(r, c);
      }

      // noexcept because number_of_rows_impl is noexcept
      size_t number_of_rows() const noexcept {
        return static_cast<Subclass const*>(this)->number_of_rows_impl();
      }

      // noexcept because number_of_cols_impl is noexcept
      size_t number_of_cols() const noexcept {
        return static_cast<Subclass const*>(this)->number_of_cols_impl();
      }

      // not noexcept because Hash<T>::operator() isn't
      size_t hash_value() const {
        return Hash<Container>()(_container);
      }

      ////////////////////////////////////////////////////////////////////////
      // Arithmetic operators - in-place
      ////////////////////////////////////////////////////////////////////////

      // not noexcept because memory is allocated
      void product_inplace_no_checks(Subclass const& A, Subclass const& B) {
        LIBSEMIGROUPS_ASSERT(number_of_rows() == number_of_cols());
        LIBSEMIGROUPS_ASSERT(A.number_of_rows() == number_of_rows());
        LIBSEMIGROUPS_ASSERT(B.number_of_rows() == number_of_rows());
        LIBSEMIGROUPS_ASSERT(A.number_of_cols() == number_of_cols());
        LIBSEMIGROUPS_ASSERT(B.number_of_cols() == number_of_cols());
        LIBSEMIGROUPS_ASSERT(&A != this);
        LIBSEMIGROUPS_ASSERT(&B != this);

        // Benchmarking boolean matrix multiplication reveals that using a
        // non-static container_type gives the best performance, when compared
        // to static container_type the performance is more or less the same
        // (but not thread-safe), and there appears to be a performance
        // penalty of about 50% when using static thread_local container_type
        // (when compiling with clang).
        size_t const             N = A.number_of_rows();
        std::vector<scalar_type> tmp(N, 0);

        for (size_t c = 0; c < N; c++) {
          for (size_t i = 0; i < N; i++) {
            tmp[i] = B(i, c);
          }
          for (size_t r = 0; r < N; r++) {
            (*this)(r, c) = std::inner_product(
                A._container.begin() + r * N,
                A._container.begin() + (r + 1) * N,
                tmp.begin(),
                scalar_zero(),
                [this](scalar_type x, scalar_type y) {
                  return this->plus_no_checks(x, y);
                },
                [this](scalar_type x, scalar_type y) {
                  return this->product_no_checks(x, y);
                });
          }
        }
      }

      // not noexcept because iterator increment isn't
      void operator*=(scalar_type a) {
        for (auto it = _container.begin(); it < _container.end(); ++it) {
          *it = product_no_checks(*it, a);
        }
      }

      // not noexcept because vector::operator[] and array::operator[] aren't
      void operator+=(Subclass const& that) {
        LIBSEMIGROUPS_ASSERT(that.number_of_rows() == number_of_rows());
        LIBSEMIGROUPS_ASSERT(that.number_of_cols() == number_of_cols());
        for (size_t i = 0; i < _container.size(); ++i) {
          _container[i] = plus_no_checks(_container[i], that._container[i]);
        }
      }

      void operator+=(RowView const& that) {
        LIBSEMIGROUPS_ASSERT(number_of_rows() == 1);
        RowView(*static_cast<Subclass const*>(this)) += that;
      }

      void operator+=(scalar_type a) {
        for (auto it = _container.begin(); it < _container.end(); ++it) {
          *it = plus_no_checks(*it, a);
        }
      }

      // TODO(2) implement operator*=(Subclass const&)

      ////////////////////////////////////////////////////////////////////////
      // Arithmetic operators - not in-place
      ////////////////////////////////////////////////////////////////////////

      // not noexcept because operator+= isn't
      Subclass operator+(Subclass const& y) const {
        Subclass result(*static_cast<Subclass const*>(this));
        result += y;
        return result;
      }

      // not noexcept because product_inplace_no_checks isn't
      Subclass operator*(Subclass const& y) const {
        Subclass result(*static_cast<Subclass const*>(this));
        result.product_inplace_no_checks(*static_cast<Subclass const*>(this),
                                         y);
        return result;
      }

      Subclass operator*(scalar_type a) const {
        Subclass result(*static_cast<Subclass const*>(this));
        result *= a;
        return result;
      }

      Subclass operator+(scalar_type a) const {
        Subclass result(*static_cast<Subclass const*>(this));
        result += a;
        return result;
      }

      ////////////////////////////////////////////////////////////////////////
      // Iterators
      ////////////////////////////////////////////////////////////////////////

      // noexcept because vector::begin and array::begin are noexcept
      iterator begin() noexcept {
        return _container.begin();
      }

      // noexcept because vector::end and array::end are noexcept
      iterator end() noexcept {
        return _container.end();
      }

      // noexcept because vector::begin and array::begin are noexcept
      const_iterator begin() const noexcept {
        return _container.begin();
      }

      // noexcept because vector::end and array::end are noexcept
      const_iterator end() const noexcept {
        return _container.end();
      }

      // noexcept because vector::cbegin and array::cbegin are noexcept
      const_iterator cbegin() const noexcept {
        return _container.cbegin();
      }

      // noexcept because vector::cend and array::cend are noexcept
      const_iterator cend() const noexcept {
        return _container.cend();
      }

      template <typename U>
      std::pair<scalar_type, scalar_type> coords(U const& it) const {
        static_assert(
            std::is_same<U, iterator>::value
                || std::is_same<U, const_iterator>::value,
            "the parameter it must be of type iterator or const_iterator");
        scalar_type const v = std::distance(_container.begin(), it);
        return std::make_pair(v / number_of_cols(), v % number_of_cols());
      }

      ////////////////////////////////////////////////////////////////////////
      // Modifiers
      ////////////////////////////////////////////////////////////////////////

      // noexcept because vector::swap and array::swap are noexcept
      void swap(MatrixCommon& that) noexcept {
        std::swap(_container, that._container);
      }

      // noexcept because swap is noexcept, and so too are number_of_rows and
      // number_of_cols
      void transpose_no_checks() noexcept {
        LIBSEMIGROUPS_ASSERT(number_of_rows() == number_of_cols());
        if (number_of_rows() == 0) {
          return;
        }
        auto& x = *this;
        for (size_t r = 0; r < number_of_rows() - 1; ++r) {
          for (size_t c = r + 1; c < number_of_cols(); ++c) {
            std::swap(x(r, c), x(c, r));
          }
        }
      }

      void transpose() {
        matrix::throw_if_not_square(static_cast<Subclass&>(*this));
        transpose_no_checks();
      }

      ////////////////////////////////////////////////////////////////////////
      // Rows
      ////////////////////////////////////////////////////////////////////////

      // not noexcept because there's an allocation
      RowView row_no_checks(size_t i) const {
        auto& container = const_cast<Container&>(_container);
        return RowView(static_cast<Subclass const*>(this),
                       container.begin() + i * number_of_cols(),
                       number_of_cols());
      }

      RowView row(size_t i) const {
        if (i >= number_of_rows()) {
          LIBSEMIGROUPS_EXCEPTION(
              "index out of range, expected value in [{}, {}), found {}",
              0,
              number_of_rows(),
              i);
        }
        return row_no_checks(i);
      }

      // not noexcept because there's an allocation
      template <typename T>
      void rows(T& x) const {
        auto& container = const_cast<Container&>(_container);
        for (auto itc = container.begin(); itc != container.end();
             itc += number_of_cols()) {
          x.emplace_back(
              static_cast<Subclass const*>(this), itc, number_of_cols());
        }
        LIBSEMIGROUPS_ASSERT(x.size() == number_of_rows());
      }

      ////////////////////////////////////////////////////////////////////////
      // Friend functions
      ////////////////////////////////////////////////////////////////////////

      friend std::ostream& operator<<(std::ostream& os, MatrixCommon const& x) {
        os << detail::to_string(x);
        return os;
      }

     private:
      ////////////////////////////////////////////////////////////////////////
      // Private data
      ////////////////////////////////////////////////////////////////////////
      container_type _container;
    };

    template <typename Scalar>
    class MatrixDynamicDim {
     public:
      MatrixDynamicDim() : _number_of_cols(0), _number_of_rows(0) {}
      MatrixDynamicDim(MatrixDynamicDim const&)            = default;
      MatrixDynamicDim(MatrixDynamicDim&&)                 = default;
      MatrixDynamicDim& operator=(MatrixDynamicDim const&) = default;
      MatrixDynamicDim& operator=(MatrixDynamicDim&&)      = default;

      MatrixDynamicDim(size_t r, size_t c)
          : _number_of_cols(c), _number_of_rows(r) {}

      ~MatrixDynamicDim() = default;

      void swap(MatrixDynamicDim& that) noexcept {
        std::swap(_number_of_cols, that._number_of_cols);
        std::swap(_number_of_rows, that._number_of_rows);
      }

     protected:
      size_t number_of_rows_impl() const noexcept {
        return _number_of_rows;
      }

      size_t number_of_cols_impl() const noexcept {
        return _number_of_cols;
      }

     private:
      size_t _number_of_cols;
      size_t _number_of_rows;
    };

    template <typename PlusOp,
              typename ProdOp,
              typename ZeroOp,
              typename OneOp,
              typename Scalar>
    struct MatrixStaticArithmetic {
      MatrixStaticArithmetic()                              = default;
      MatrixStaticArithmetic(MatrixStaticArithmetic const&) = default;
      MatrixStaticArithmetic(MatrixStaticArithmetic&&)      = default;
      MatrixStaticArithmetic& operator=(MatrixStaticArithmetic const&)
          = default;
      MatrixStaticArithmetic& operator=(MatrixStaticArithmetic&&) = default;

      // TODO(2) from here to the end of MatrixStaticArithmetic should be
      // private or protected
      using scalar_type = Scalar;

      static constexpr scalar_type plus_no_checks_impl(scalar_type x,
                                                       scalar_type y) noexcept {
        return PlusOp()(x, y);
      }

      static constexpr scalar_type
      product_no_checks_impl(scalar_type x, scalar_type y) noexcept {
        return ProdOp()(x, y);
      }

      static constexpr scalar_type one_impl() noexcept {
        return OneOp()();
      }

      static constexpr scalar_type zero_impl() noexcept {
        return ZeroOp()();
      }

      static constexpr void const* semiring_impl() noexcept {
        return nullptr;
      }
    };

    ////////////////////////////////////////////////////////////////////////
    // RowViews - class for cheaply storing iterators to rows
    ////////////////////////////////////////////////////////////////////////

    template <typename Mat, typename Subclass>
    class RowViewCommon {
      static_assert(IsMatrix<Mat>,
                    "the template parameter Mat must be derived from "
                    "MatrixPolymorphicBase");

     public:
      using const_iterator = typename Mat::const_iterator;
      using iterator       = typename Mat::iterator;

      using scalar_type            = typename Mat::scalar_type;
      using scalar_reference       = typename Mat::scalar_reference;
      using scalar_const_reference = typename Mat::scalar_const_reference;

      using Row         = typename Mat::Row;
      using matrix_type = Mat;

      size_t size() const noexcept {
        return static_cast<Subclass const*>(this)->length_impl();
      }

     private:
      scalar_type plus_no_checks(scalar_type x, scalar_type y) const noexcept {
        return static_cast<Subclass const*>(this)->plus_no_checks_impl(y, x);
      }

      scalar_type product_no_checks(scalar_type x,
                                    scalar_type y) const noexcept {
        return static_cast<Subclass const*>(this)->product_no_checks_impl(y, x);
      }

     public:
      RowViewCommon()                                = default;
      RowViewCommon(RowViewCommon const&)            = default;
      RowViewCommon(RowViewCommon&&)                 = default;
      RowViewCommon& operator=(RowViewCommon const&) = default;
      RowViewCommon& operator=(RowViewCommon&&)      = default;

      explicit RowViewCommon(Row const& r)
          : RowViewCommon(const_cast<Row&>(r).begin()) {}

      // Not noexcept because iterator::operator[] isn't
      scalar_const_reference operator[](size_t i) const {
        return _begin[i];
      }

      // Not noexcept because iterator::operator[] isn't
      scalar_reference operator[](size_t i) {
        return _begin[i];
      }

      // Not noexcept because iterator::operator[] isn't
      scalar_const_reference operator()(size_t i) const {
        return (*this)[i];
      }

      // Not noexcept because iterator::operator[] isn't
      scalar_reference operator()(size_t i) {
        return (*this)[i];
      }

      // noexcept because begin() is
      const_iterator cbegin() const noexcept {
        return _begin;
      }

      // not noexcept because iterator arithmetic isn't
      const_iterator cend() const {
        return _begin + size();
      }

      // noexcept because begin() is
      const_iterator begin() const noexcept {
        return _begin;
      }

      // not noexcept because iterator arithmetic isn't
      const_iterator end() const {
        return _begin + size();
      }

      // noexcept because begin() is
      iterator begin() noexcept {
        return _begin;
      }

      // not noexcept because iterator arithmetic isn't
      iterator end() noexcept {
        return _begin + size();
      }

      ////////////////////////////////////////////////////////////////////////
      // Arithmetic operators
      ////////////////////////////////////////////////////////////////////////

      // not noexcept because operator[] isn't
      void operator+=(RowViewCommon const& x) {
        auto& this_ = *this;
        for (size_t i = 0; i < size(); ++i) {
          this_[i] = plus_no_checks(this_[i], x[i]);
        }
      }

      // not noexcept because iterator arithmeic isn't
      void operator+=(scalar_type a) {
        for (auto& x : *this) {
          x = plus_no_checks(x, a);
        }
      }

      // not noexcept because iterator arithmeic isn't
      void operator*=(scalar_type a) {
        for (auto& x : *this) {
          x = product_no_checks(x, a);
        }
      }

      // not noexcept because operator*= isn't
      Row operator*(scalar_type a) const {
        Row result(*static_cast<Subclass const*>(this));
        result *= a;
        return result;
      }

      // not noexcept because operator+= isn't
      Row operator+(RowViewCommon const& that) const {
        Row result(*static_cast<Subclass const*>(this));
        result += static_cast<Subclass const&>(that);
        return result;
      }

      template <typename U>
      bool operator==(U const& that) const {
        // TODO(1) static assert that U is Row or RowView
        return std::equal(begin(), end(), that.begin());
      }

      template <typename U>
      bool operator!=(U const& that) const {
        return !(*this == that);
      }

      template <typename U>
      bool operator<(U const& that) const {
        return std::lexicographical_compare(
            cbegin(), cend(), that.cbegin(), that.cend());
      }

      template <typename U>
      bool operator>(U const& that) const {
        return that < *this;
      }

      void swap(RowViewCommon& that) noexcept {
        std::swap(that._begin, _begin);
      }

      friend std::ostream& operator<<(std::ostream&        os,
                                      RowViewCommon const& x) {
        os << detail::to_string(x);
        return os;
      }

     protected:
      explicit RowViewCommon(iterator first) : _begin(first) {}

     private:
      iterator _begin;
    };

    template <typename Container>
    void throw_if_any_row_wrong_size(Container const& m) {
      if (m.size() <= 1) {
        return;
      }
      uint64_t const C  = m.begin()->size();
      auto           it = std::find_if_not(
          m.begin() + 1, m.end(), [&C](typename Container::const_reference r) {
            return r.size() == C;
          });
      if (it != m.end()) {
        LIBSEMIGROUPS_EXCEPTION("invalid argument, expected every item to "
                                "have length {}, found {} in entry {}",
                                C,
                                it->size(),
                                std::distance(m.begin(), it));
      }
    }

    template <typename Scalar>
    void throw_if_any_row_wrong_size(
        std::initializer_list<std::initializer_list<Scalar>> m) {
      throw_if_any_row_wrong_size<
          std::initializer_list<std::initializer_list<Scalar>>>(m);
    }

  }  // namespace detail

  ////////////////////////////////////////////////////////////////////////
  // Matrix forward declarations
  ////////////////////////////////////////////////////////////////////////

  template <typename PlusOp,
            typename ProdOp,
            typename ZeroOp,
            typename OneOp,
            size_t R,
            size_t C,
            typename Scalar>
  class StaticMatrix;

  template <typename... Args>
  class DynamicMatrix;

  template <typename PlusOp,
            typename ProdOp,
            typename ZeroOp,
            typename OneOp,
            typename Scalar>
  class DynamicMatrix<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>;

  template <typename Semiring, typename Scalar>
  class DynamicMatrix<Semiring, Scalar>;

  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  // StaticRowViews - static arithmetic
  ////////////////////////////////////////////////////////////////////////

  //! Defined in `matrix.hpp`.
  //!
  //! This is a class for views into a row of a matrix over a semiring. An
  //! instance of this class provides access to the elements in a row of a
  //! matrix and is cheap to create and copy. Addition, scalar multiplication,
  //! and other standard vector operations are defined for row views.
  //!
  //! This class is the type of row views into a StaticMatrix; see
  //! the documentation for StaticMatrix for further details.
  //!
  //! \tparam PlusOp a stateless type with a call operator of signature
  //! `scalar_type operator()(scalar_type, scalar_type)` implementing the
  //! addition of the semiring.
  //!
  //! \tparam ProdOp a stateless type with a call operator of signature
  //! `scalar_type operator()(scalar_type, scalar_type)` implementing the
  //! multiplication of the semiring.
  //!
  //! \tparam ZeroOp a stateless type with a call operator of signature
  //! `scalar_type operator()()` returning the zero of the semiring (the
  //! additive identity element).
  //!
  //! \tparam OneOp a stateless type with a call operator of signature
  //! `scalar_type operator()()` returning the one of the semiring (the
  //! multiplicative identity element).
  //!
  //! \tparam C  the number of columns of the underlying matrix.
  //!
  //! \tparam Scalar the type of the entries in the matrices (the type of
  //! elements in the underlying semiring).
  //!
  //! \warning
  //! If the underlying matrix is destroyed, then any row views for that matrix
  //! are invalidated.
  template <typename PlusOp,
            typename ProdOp,
            typename ZeroOp,
            typename OneOp,
            size_t C,
            typename Scalar>
  class StaticRowView
      : public detail::RowViewCommon<
            StaticMatrix<PlusOp, ProdOp, ZeroOp, OneOp, 1, C, Scalar>,
            StaticRowView<PlusOp, ProdOp, ZeroOp, OneOp, C, Scalar>>,
        public detail::
            MatrixStaticArithmetic<PlusOp, ProdOp, ZeroOp, OneOp, Scalar> {
   private:
    using RowViewCommon = detail::RowViewCommon<
        StaticMatrix<PlusOp, ProdOp, ZeroOp, OneOp, 1, C, Scalar>,
        StaticRowView<PlusOp, ProdOp, ZeroOp, OneOp, C, Scalar>>;
    friend RowViewCommon;

    template <size_t R>
    using StaticMatrix_ = ::libsemigroups::
        StaticMatrix<PlusOp, ProdOp, ZeroOp, OneOp, R, C, Scalar>;

   public:
    //! \copydoc StaticMatrix::const_iterator
    using const_iterator = typename RowViewCommon::const_iterator;

    //! \copydoc StaticMatrix::const_iterator
    using iterator = typename RowViewCommon::iterator;

    //! \copydoc StaticMatrix::scalar_type
    using scalar_type = Scalar;

    //! \copydoc StaticMatrix::scalar_reference
    using scalar_reference = typename RowViewCommon::scalar_reference;

    //! \copydoc StaticMatrix::scalar_const_reference
    using scalar_const_reference =
        typename RowViewCommon::scalar_const_reference;

    //! Alias for the type of the underlying matrix.
    using matrix_type = typename RowViewCommon::matrix_type;

    //! Alias for the type of a row in the underlying matrix.
    using Row = typename matrix_type::Row;

    //! \copydoc StaticMatrix::StaticMatrix()
    StaticRowView() = default;

    //! \copydoc StaticMatrix::StaticMatrix(StaticMatrix const&)
    StaticRowView(StaticRowView const&) = default;

    //! \copydoc StaticMatrix::StaticMatrix(StaticMatrix&&)
    StaticRowView(StaticRowView&&) = default;

    //! \copydoc StaticMatrix::operator=(StaticMatrix const&)
    StaticRowView& operator=(StaticRowView const&) = default;

    //! \copydoc StaticMatrix::operator=(StaticMatrix&&)
    StaticRowView& operator=(StaticRowView&&) = default;

#ifndef PARSED_BY_DOXYGEN
    using RowViewCommon::RowViewCommon;

    // TODO(2) This constructor should be private
    template <size_t R>
    StaticRowView(StaticMatrix_<R> const*,
                  typename RowViewCommon::iterator it,
                  size_t)
        : RowViewCommon(it) {}

    using RowViewCommon::size;
#else
    //! \brief Construct a row view from a \ref Row.
    //!
    //! This function constructs a row view from a \ref Row.
    //!
    //! \param r  the row.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    explicit StaticRowView(Row const& r);

    //! \brief Returns the size of the row.
    //!
    //! Returns the size of the underlying row which equals the template
    //! parameter \p C.
    //!
    //! \returns  The size of the row.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant
    static constexpr size_t size() const noexcept;

    //! \brief Returns a iterator pointing at the first entry.
    //!
    //! This function returns a (random access) iterator pointing at the first
    //! entry in the row.
    //!
    //! \returns  A value of type \ref iterator.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant
    iterator begin() noexcept;

    //! \brief Returns a iterator pointing one beyond the last entry of the row.
    //!
    //! This function returns a (random access) iterator pointing at one beyond
    //! the last entry in the row.
    //!
    //! \returns  A value of type \ref iterator.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant
    iterator end();

    //! \brief Returns a const iterator pointing at the first entry.
    //!
    //! This function returns a const (random access) iterator pointing at the
    //! first entry in the row.
    //!
    //! \returns  A value of type \ref const_iterator.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant
    const_iterator cbegin() const noexcept;

    //! \brief Returns a const iterator pointing one beyond the last entry of
    //! the row.
    //!
    //! This function returns a (random access) const iterator pointing at one
    //! beyond the last entry in the row.
    //!
    //! \returns  A value of type \ref iterator.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant
    iterator cend();

    //! \brief Returns a reference to the specified entry of the row.
    //!
    //! This function returns a reference to the specified entry of the
    //! row.
    //!
    //! \param i  the index of the entry.
    //!
    //! \returns  A value of type \ref scalar_reference.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant
    //!
    //! \warning
    //! No checks on the validity of the parameter \p i is performed.
    scalar_reference operator()(size_t i);

    //! \brief Returns a const reference to the specified entry of the row.
    //!
    //! This function returns a const reference to the specified entry of the
    //! row.
    //!
    //! \param i  the index of the entry.
    //!
    //! \returns  A value of type \ref scalar_const_reference.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant
    //!
    //! \warning
    //! No checks on the validity of the parameter \p i is performed.
    scalar_const_reference operator()(size_t i) const;

    //! \copydoc operator()(size_t)
    scalar_reference operator[](size_t i);

    //! \copydoc operator()(size_t) const
    scalar_const_reference operator[](size_t i) const;

    //! \brief Equality operator.
    //!
    //! Equality operator.
    //!
    //! \tparam U  either \ref Row, \ref StaticRowView,
    //! \ref DynamicRowViewStaticArith
    //! "DynamicRowView (compile-time arithmetic)", or \ref
    //! DynamicRowViewDynamicArith "DynamicRowView (run-time arithmetic)" \param
    //! that  \ref Row, \ref StaticRowView,
    //! \ref DynamicRowViewStaticArith
    //! "DynamicRowView (compile-time arithmetic)", or \ref
    //! DynamicRowViewDynamicArith "DynamicRowView (run-time arithmetic)" object
    //! for comparison.
    //!
    //! \returns
    //! Returns \c true if the first \ref size entries of \p that
    //! equal the entries of `this`, and \c false if not.
    //!
    //! \complexity
    //! At worst \f$O(n)\f$ where \f$n\f$ is the size of the row view.
    template <typename U>
    bool operator==(U const& that) const;

    //! \brief Inequality operator.
    //!
    //! Inequality operator.
    //!
    //! \tparam U  either \ref Row, \ref StaticRowView,
    //! \ref DynamicRowViewStaticArith
    //! "DynamicRowView (compile-time arithmetic)", or \ref
    //! DynamicRowViewDynamicArith "DynamicRowView (run-time arithmetic)" \param
    //! that  \ref Row, \ref StaticRowView,
    //! \ref DynamicRowViewStaticArith
    //! "DynamicRowView (compile-time arithmetic)", or \ref
    //! DynamicRowViewDynamicArith "DynamicRowView (run-time arithmetic)" object
    //! for comparison.
    //!
    //! \returns
    //! Returns \c true if the first \ref size entries of \p that
    //! do not equal the entries of `this`, and \c false otherwise.
    //!
    //! \complexity
    //! At worst \f$O(n)\f$ where \f$n\f$ is the size of the row view.
    template <typename U>
    bool operator!=(U const& that) const;

    //! \brief Lexicographical comparison of rows.
    //!
    //! Lexicographical comparison of rows.
    //!
    //! \tparam U  either \ref Row, \ref StaticRowView, \ref
    //! DynamicRowViewStaticArith "DynamicRowView (compile-time arithmetic)", or
    //! \ref DynamicRowViewDynamicArith "DynamicRowView (run-time arithmetic)".
    //!
    //! \param that \ref Row, \ref StaticRowView, \ref
    //! DynamicRowViewStaticArith "DynamicRowView (compile-time arithmetic)", or
    //! \ref DynamicRowViewDynamicArith "DynamicRowView (run-time arithmetic)"
    //! object for comparison.
    //!
    //! \returns
    //! Returns \c true if `*this` is lex less than that `that` and
    //! \c false otherwise.
    //!
    //! \complexity
    //! At worst \f$O(n)\f$ where \f$n\f$ is the size of the row view.
    template <typename U>
    bool operator<(U const& that) const;

    //! \brief Lexicographical comparison of rows.
    //!
    //! Lexicographical comparison of rows.
    //!
    //! \tparam U  either \ref Row, \ref StaticRowView,
    //! \ref DynamicRowViewStaticArith
    //! "DynamicRowView (compile-time arithmetic)", or \ref
    //! DynamicRowViewDynamicArith "DynamicRowView (run-time arithmetic)"
    //!
    //! \param
    //! that \ref Row, \ref StaticRowView, \ref DynamicRowViewStaticArith
    //! "DynamicRowView (compile-time arithmetic)", or \ref
    //! DynamicRowViewDynamicArith "DynamicRowView (run-time arithmetic)" object
    //! for comparison.
    //!
    //! \returns
    //! Returns \c true if `this` is lex greater than that `that` and
    //! \c false otherwise.
    //!
    //! \complexity
    //! At worst \f$O(n)\f$ where \f$n\f$ is the \ref size of the row view.
    template <typename U>
    bool operator<(U const& that) const;

    //! \brief Sum row views.
    //!
    //! This function sums a row view with another row view and returns a newly
    //! allocated \ref Row.
    //!
    //! \param that  the row view to add.
    //! \returns  A \ref Row containing the sum.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(m)\f$ where \f$m\f$ is \ref size.
    //!
    //! \warning
    //! The two row views must be of the same size, although this is not
    //! verified by the implementation.
    //!
    //! \warning This function does not detect overflows of \ref scalar_type.
    Row operator+(StaticRowView const& that);

    //! \brief Sums a row view with another row view in-place.
    //!
    //! This function redefines a RowView object to be the sum of itself and \p
    //! that.
    //!
    //! \param that  the row view to add.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(m)\f$ where \f$m\f$ is \ref size.
    //!
    //! \warning
    //! The two row views must be of the same size, although this is not
    //! verified by the implementation.
    //!
    //! \warning This function does not detect overflows of \ref scalar_type.
    void operator+=(StaticRowView const& that);

    //! \brief Adds a scalar to every entry of a row in-place.
    //!
    //! This function adds a scalar to every entry of the row in-place.
    //!
    //! \param a the scalar to add.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(m)\f$ where \f$m\f$ is \ref size.
    //!
    //! \warning This function does not detect overflows of \ref scalar_type.
    void operator+=(scalar_type a);

    //! \brief Multiplies every entry of the row by a scalar.
    //!
    //! This function multiplies every entry of the row by a scalar.
    //!
    //! \param a the scalar.
    //!
    //! \returns
    //! A value of type \ref Row containing the original row multiplied by the
    //! scalar \p a.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(m)\f$ where \f$m\f$ is \ref size.
    //!
    //! \warning This function does not detect overflows of \ref scalar_type.
    Row operator*(scalar_type a) const;

    //! \brief Multiplies every entry of the row by a scalar in-place.
    //!
    //! This function multiplies every entry of the row by a scalar in-place.
    //!
    //! \param a  the scalar to multiply by.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(m)\f$ where \f$m\f$ is \ref size.
    //!
    //! \warning This function does not detect overflows of \ref scalar_type.
    void operator*=(scalar_type a);
#endif

   private:
    static constexpr size_t length_impl() noexcept {
      return C;
    }
  };

  ////////////////////////////////////////////////////////////////////////
  // DynamicRowViews - static arithmetic
  ////////////////////////////////////////////////////////////////////////

  template <typename... Args>
  class DynamicRowView;

  //! \anchor DynamicRowViewStaticArith
  //!
  //! \brief Row view class for \ref DynamicMatrixStaticArith
  //! "DynamicMatrix (compile-time arithmetic)".
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This class is the type of row views into a \ref DynamicMatrixStaticArith
  //! "DynamicMatrix (with compile-time arithmetic)" with the same template
  //! parameters.
  //!
  //! An instance of this class provides access to the elements in a row of a
  //! matrix and is cheap to create and copy. Addition, scalar multiplication,
  //! and other usual vector operations are defined for row views.
  //!
  //! \tparam PlusOp a stateless type with a call operator of signature
  //! `scalar_type operator()(scalar_type, scalar_type)` implementing the
  //! addition of the semiring
  //!
  //! \tparam ProdOp a stateless type with a call operator of signature
  //! `scalar_type operator()(scalar_type, scalar_type)` implementing the
  //! multiplication of the semiring
  //!
  //! \tparam ZeroOp a stateless type with a call operator of signature
  //! `scalar_type operator()()` returning the zero of the semiring (the
  //! additive identity element)
  //!
  //! \tparam OneOp a stateless type with a call operator of signature
  //! `scalar_type operator()()` returning the one of the semiring (the
  //! multiplicative identity element)
  //!
  //! \tparam Scalar the type of the entries in the matrices (the type of
  //! elements in the underlying semiring)
  //!
  //! \warning
  //! If the underlying matrix is destroyed, then any row views for that
  //! matrix are invalidated.
  template <typename PlusOp,
            typename ProdOp,
            typename ZeroOp,
            typename OneOp,
            typename Scalar>
  class DynamicRowView<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>
      : public detail::RowViewCommon<
            DynamicMatrix<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>,
            DynamicRowView<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>>,
        public detail::
            MatrixStaticArithmetic<PlusOp, ProdOp, ZeroOp, OneOp, Scalar> {
   private:
    using DynamicMatrix_ = DynamicMatrix<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>;
    using RowViewCommon  = detail::RowViewCommon<
        DynamicMatrix_,
        DynamicRowView<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>>;
    friend RowViewCommon;

   public:
    //! \copydoc StaticMatrix::const_iterator
    using const_iterator = typename RowViewCommon::const_iterator;

    //! \copydoc StaticMatrix::const_iterator
    using iterator = typename RowViewCommon::iterator;

    //! \copydoc StaticMatrix::scalar_type
    using scalar_type = Scalar;

    //! \copydoc StaticMatrix::scalar_reference
    using scalar_reference = typename RowViewCommon::scalar_reference;

    //! \copydoc StaticMatrix::scalar_const_reference
    using scalar_const_reference =
        typename RowViewCommon::scalar_const_reference;

    //! Alias for the type of the underlying matrix.
    using matrix_type = typename RowViewCommon::matrix_type;

    //! Alias for the type of a row in the underlying matrix.
    using Row = typename matrix_type::Row;

    //! \copydoc StaticMatrix::StaticMatrix()
    DynamicRowView() = default;

    //! \copydoc StaticMatrix::StaticMatrix(StaticMatrix const&)
    DynamicRowView(DynamicRowView const&) = default;

    //! \copydoc StaticMatrix::StaticMatrix(StaticMatrix&&)
    DynamicRowView(DynamicRowView&&) = default;

    //! \copydoc StaticMatrix::operator=(StaticMatrix const&)
    DynamicRowView& operator=(DynamicRowView const&) = default;

    //! \copydoc StaticMatrix::operator=(StaticMatrix&&)
    DynamicRowView& operator=(DynamicRowView&&) = default;

    //! \copydoc StaticRowView::StaticRowView(Row const&)
    explicit DynamicRowView(Row const& r)
        : RowViewCommon(r), _length(r.number_of_cols()) {}

#ifndef PARSED_BY_DOXYGEN
    using RowViewCommon::RowViewCommon;

    // TODO(2) This constructor should be private
    DynamicRowView(DynamicMatrix_ const*, iterator const& it, size_t N)
        : RowViewCommon(it), _length(N) {}

    using RowViewCommon::size;
#else
    //! \copydoc StaticRowView::size
    size_t size() const noexcept;

    //! \copydoc StaticRowView::begin
    iterator begin() noexcept;

    //! \copydoc StaticRowView::end
    iterator end();

    //! \copydoc StaticRowView::cbegin
    const_iterator cbegin() const noexcept;

    //! \copydoc StaticRowView::cend
    iterator cend();

    //! \copydoc StaticRowView::operator()(size_t)
    scalar_reference operator()(size_t i);

    //! \copydoc StaticRowView::operator()(size_t) const
    scalar_const_reference operator()(size_t i) const;

    //! \copydoc StaticRowView::operator[](size_t)
    scalar_reference operator[](size_t i);

    //! \copydoc StaticRowView::operator[](size_t) const
    scalar_const_reference operator[](size_t i) const;

    //! \copydoc StaticRowView::operator==
    template <typename U>
    bool operator==(U const& that) const;

    //! \copydoc StaticRowView::operator!=
    template <typename U>
    bool operator!=(U const& that) const;

    //! \copydoc StaticRowView::operator<
    template <typename U>
    bool operator<(U const& that) const;

    //! \copydoc StaticRowView::operator<
    template <typename U>
    bool operator<(U const& that) const;

    //! \copydoc StaticRowView::operator+
    Row operator+(DynamicRowView const& that);

    //! \copydoc StaticRowView::operator+=(StaticRowView const&)
    void operator+=(DynamicRowView const& that);

    //! \copydoc StaticRowView::operator+=(scalar_type)
    void operator+=(scalar_type a);

    //! \copydoc StaticRowView::operator*(scalar_type) const
    Row operator*(scalar_type a) const;

    //! \copydoc StaticRowView::operator*=(scalar_type)
    void operator*=(scalar_type a);
#endif  // PARSED_BY_DOXYGEN

   private:
    size_t length_impl() const noexcept {
      return _length;
    }
    size_t _length;
  };

  ////////////////////////////////////////////////////////////////////////
  // DynamicRowViews - dynamic arithmetic
  ////////////////////////////////////////////////////////////////////////

  //! \anchor DynamicRowViewDynamicArith
  //!
  //! \brief Row view class for \ref DynamicMatrixDynamicArith
  //! "DynamicMatrix (run-time arithmetic)".
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This class is the type of row views into a \ref
  //! DynamicMatrixDynamicArith "DynamicMatrix (run-time arithmetic)".
  //!
  //! \tparam Semiring the type of a semiring object which defines the semiring
  //! arithmetic (see requirements in \ref DynamicMatrixDynamicArith
  //! "DynamicMatrix (run-time arithmetic)".
  //!
  //! \tparam Scalar the type of the entries in the matrices (the type of
  //! elements in the underlying semiring).
  //!
  //! \warning
  //! If the underlying matrix is destroyed, then any row views for that
  //! matrix are invalidated.
  template <typename Semiring, typename Scalar>
  class DynamicRowView<Semiring, Scalar>
      : public detail::RowViewCommon<DynamicMatrix<Semiring, Scalar>,
                                     DynamicRowView<Semiring, Scalar>> {
   private:
    using DynamicMatrix_ = DynamicMatrix<Semiring, Scalar>;
    friend DynamicMatrix_;
    using RowViewCommon
        = detail::RowViewCommon<DynamicMatrix_,
                                DynamicRowView<Semiring, Scalar>>;
    friend RowViewCommon;

   public:
    //! \copydoc StaticMatrix::const_iterator
    using const_iterator = typename RowViewCommon::const_iterator;

    //! \copydoc StaticMatrix::const_iterator
    using iterator = typename RowViewCommon::iterator;

    //! \copydoc StaticMatrix::scalar_type
    using scalar_type = Scalar;

    //! \copydoc StaticMatrix::scalar_reference
    using scalar_reference = typename RowViewCommon::scalar_reference;

    //! \copydoc StaticMatrix::scalar_const_reference
    using scalar_const_reference =
        typename RowViewCommon::scalar_const_reference;

    //! Alias for the type of the underlying matrix.
    using matrix_type = typename RowViewCommon::matrix_type;

    //! Alias for the type of a row in the underlying matrix.
    using Row = typename matrix_type::Row;

    //! \copydoc StaticMatrix::StaticMatrix()
    DynamicRowView() = default;

    //! \copydoc StaticMatrix::StaticMatrix(StaticMatrix const&)
    DynamicRowView(DynamicRowView const&) = default;

    //! \copydoc StaticMatrix::StaticMatrix(StaticMatrix&&)
    DynamicRowView(DynamicRowView&&) = default;

    //! \copydoc StaticMatrix::operator=(StaticMatrix const&)
    DynamicRowView& operator=(DynamicRowView const&) = default;

    //! \copydoc StaticMatrix::operator=(StaticMatrix&&)
    DynamicRowView& operator=(DynamicRowView&&) = default;

    //! \copydoc StaticRowView::StaticRowView(Row const&)
    explicit DynamicRowView(Row const& r) : RowViewCommon(r), _matrix(&r) {}

#ifndef PARSED_BY_DOXYGEN
    using RowViewCommon::RowViewCommon;

    // TODO(2) This constructor should be private
    DynamicRowView(DynamicMatrix_ const* mat, iterator const& it, size_t)
        : RowViewCommon(it), _matrix(mat) {}

    using RowViewCommon::size;
#else
    //! \copydoc StaticRowView::size
    size_t size() const noexcept;

    //! \copydoc StaticRowView::begin
    iterator begin() noexcept;

    //! \copydoc StaticRowView::end
    iterator end();

    //! \copydoc StaticRowView::cbegin
    const_iterator cbegin() const noexcept;

    //! \copydoc StaticRowView::cend
    iterator cend();

    //! \copydoc StaticRowView::operator()(size_t)
    scalar_reference operator()(size_t i);

    //! \copydoc StaticRowView::operator()(size_t) const
    scalar_const_reference operator()(size_t i) const;

    //! \copydoc StaticRowView::operator[](size_t)
    scalar_reference operator[](size_t i);

    //! \copydoc StaticRowView::operator[](size_t) const
    scalar_const_reference operator[](size_t i) const;

    //! \copydoc StaticRowView::operator==
    template <typename U>
    bool operator==(U const& that) const;

    //! \copydoc StaticRowView::operator!=
    template <typename U>
    bool operator!=(U const& that) const;

    //! \copydoc StaticRowView::operator<
    template <typename U>
    bool operator<(U const& that) const;

    //! \copydoc StaticRowView::operator<
    template <typename U>
    bool operator<(U const& that) const;

    //! \copydoc StaticRowView::operator+
    Row operator+(DynamicRowView const& that);

    //! \copydoc StaticRowView::operator+=(StaticRowView const&)
    void operator+=(DynamicRowView const& that);

    //! \copydoc StaticRowView::operator+=(scalar_type)
    void operator+=(scalar_type a);

    //! \copydoc StaticRowView::operator*(scalar_type) const
    Row operator*(scalar_type a) const;

    //! \copydoc StaticRowView::operator*=(scalar_type)
    void operator*=(scalar_type a);
#endif

   private:
    size_t length_impl() const noexcept {
      return _matrix->number_of_cols();
    }

    scalar_type plus_no_checks_impl(scalar_type x,
                                    scalar_type y) const noexcept {
      return _matrix->plus_no_checks_impl(x, y);
    }

    scalar_type product_no_checks_impl(scalar_type x,
                                       scalar_type y) const noexcept {
      return _matrix->product_no_checks_impl(x, y);
    }

    DynamicMatrix_ const* _matrix;
  };

  ////////////////////////////////////////////////////////////////////////
  // StaticMatrix with compile time semiring arithmetic
  ////////////////////////////////////////////////////////////////////////

  //! \brief Static matrix class
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a class for matrices where both the arithmetic operations in
  //! the underlying semiring, and the dimensions of the matrices are known
  //! at compile time.
  //!
  //! \tparam PlusOp a stateless type with a call operator of signature
  //!    `scalar_type operator()(scalar_type, scalar_type)`
  //!    implementing the addition of the semiring.
  //!
  //! \tparam ProdOp a stateless type with a call operator of signature
  //!    `scalar_type operator()(scalar_type, scalar_type)`
  //!    implementing the multiplication of the semiring.
  //!
  //! \tparam ZeroOp a stateless type with a call operator of signature
  //!    `scalar_type operator()()` returning the zero of the semiring
  //!    (the additive identity element).
  //!
  //! \tparam OneOp a stateless type with a call operator of signature
  //!    `scalar_type operator()()` returning the one of the semiring
  //!    (the multiplicative identity element).
  //!
  //! \tparam R  the number of rows of the matrix.
  //! \tparam C  the number of columns of the matrix.
  //!
  //! \tparam Scalar the type of the entries in the matrices (the type of
  //! elements in the underlying semiring)
  //!
  //! \note Certain member functions only work for square matrices and some only
  //!    work for rows.
  template <typename PlusOp,
            typename ProdOp,
            typename ZeroOp,
            typename OneOp,
            size_t R,
            size_t C,
            typename Scalar>
  class StaticMatrix
      : public detail::
            MatrixStaticArithmetic<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>,
        public detail::MatrixCommon<
            std::array<Scalar, R * C>,
            StaticMatrix<PlusOp, ProdOp, ZeroOp, OneOp, R, C, Scalar>,
            StaticRowView<PlusOp, ProdOp, ZeroOp, OneOp, C, Scalar>> {
    ////////////////////////////////////////////////////////////////////////
    // StaticMatrix - Aliases - private
    ////////////////////////////////////////////////////////////////////////

    using MatrixCommon = ::libsemigroups::detail::MatrixCommon<
        std::array<Scalar, R * C>,
        StaticMatrix<PlusOp, ProdOp, ZeroOp, OneOp, R, C, Scalar>,
        StaticRowView<PlusOp, ProdOp, ZeroOp, OneOp, C, Scalar>>;
    friend MatrixCommon;

   public:
    ////////////////////////////////////////////////////////////////////////
    // StaticMatrix - Aliases - public
    ////////////////////////////////////////////////////////////////////////

    //! Alias for the template parameter \p Scalar.
    using scalar_type = typename MatrixCommon::scalar_type;

    //! Alias for references to the template parameter \p Scalar.
    using scalar_reference = typename MatrixCommon::scalar_reference;

    //! Alias for const references to the template parameter \p Scalar.
    using scalar_const_reference =
        typename MatrixCommon::scalar_const_reference;

    //! Alias for the type of the rows of a StaticMatrix.
    using Row = StaticMatrix<PlusOp, ProdOp, ZeroOp, OneOp, 1, C, Scalar>;

    //! Alias for the type of row views of a StaticMatrix.
    using RowView = StaticRowView<PlusOp, ProdOp, ZeroOp, OneOp, C, Scalar>;

    //! Alias for the template parameter \p PlusOp.
    using Plus = PlusOp;

    //! Alias for the template parameter \p ProdOp.
    using Prod = ProdOp;

    //! Alias for the template parameter \p ZeroOp.
    using Zero = ZeroOp;

    //! Alias for the template parameter \p OneOp.
    using One = OneOp;

    //! Alias for iterators pointing at entries of a matrix.
    using iterator = typename MatrixCommon::iterator;

    //! Alias for const iterators pointing at entries of a matrix.
    using const_iterator = typename MatrixCommon::const_iterator;

    static constexpr size_t nr_rows = R;
    static constexpr size_t nr_cols = C;

    ////////////////////////////////////////////////////////////////////////
    // StaticMatrix - Constructors + destructor - public
    ////////////////////////////////////////////////////////////////////////

    //! \brief Construct a row (from std::initializer_list).
    //!
    //! Construct a row from an std::initializer_list.
    //!
    //! \param c  the values to be copied into the row.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //!   \f$O(n)\f$ where \f$n\f$ is the size of the row being
    //!   constructed.
    //!
    //! \warning
    //!   This constructor only works for rows, i.e. when the template
    //!   parameter \c R is \c 1.
    explicit StaticMatrix(std::initializer_list<scalar_type> const& c)
        : MatrixCommon(c) {
      static_assert(R == 1,
                    "cannot construct Matrix from the given initializer list, "
                    "incompatible dimensions");
      LIBSEMIGROUPS_ASSERT(c.size() == C);
    }

    //! \brief Construct a matrix.
    //!
    //! Construct a matrix from an std::initializer_list of
    //! std::initializer_list.
    //!
    //! \param m  the values to be copied into the matrix.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //!   \f$O(mn)\f$ where \f$m\f$ is the template parameter \c R
    //!   and \f$n\f$ is the template parameter \c C.
    //!
    //! \par Example
    //! \code
    //!  Mat m({{1, 1}, {0, 0}});
    //! \endcode
    explicit StaticMatrix(
        std::initializer_list<std::initializer_list<scalar_type>> const& m)
        : MatrixCommon(m) {}

    //! \brief Construct a matrix.
    //!
    //! Construct a matrix from an std::vector of std::vector.
    //!
    //! \param m  the values to be copied into the matrix.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //!   \f$O(mn)\f$ where \f$m\f$ is the template parameter \c R
    //!   and \f$n\f$ is the template parameter \c C.
    explicit StaticMatrix(std::vector<std::vector<scalar_type>> const& m)
        : MatrixCommon(m) {}

    //! \brief Construct a row from a row view.
    //!
    //! Construct a row from a row view.
    //!
    //! \param rv  the row view.
    //!
    //! \exceptions
    //!   this function guarantees not to throw a
    //!   LibsemigroupsException.
    //!
    //! \complexity
    //!   \f$O(n)\f$ where \f$n\f$ is the size of the row being
    //!   constructed.
    //!
    //! \warning
    //!   This constructor will only compile for rows, i.e. when the template
    //!   parameter \c R is \c 1.
    explicit StaticMatrix(RowView const& rv) : MatrixCommon(rv) {
      static_assert(
          R == 1,
          "cannot construct Matrix with more than one row from RowView!");
    }

    //! \brief Default constructor.
    //!
    //! Default constructor.
    StaticMatrix() = default;

    //! \brief Default copy constructor.
    //!
    //! Default copy constructor.
    StaticMatrix(StaticMatrix const&) = default;

    //! \brief Default move constructor.
    //!
    //! Default move constructor.
    StaticMatrix(StaticMatrix&&) = default;

    //! \brief Default copy assignment operator.
    //!
    //! Default copy assignment operator.
    StaticMatrix& operator=(StaticMatrix const&) = default;

    //! \brief Default move assignment operator.
    //!
    //! Default move assignment operator.
    StaticMatrix& operator=(StaticMatrix&&) = default;

#ifndef PARSED_BY_DOXYGEN
    // For uniformity of interface, the args do nothing
    StaticMatrix(size_t r, size_t c) : StaticMatrix() {
      (void) r;
      (void) c;
      LIBSEMIGROUPS_ASSERT(r == number_of_rows());
      LIBSEMIGROUPS_ASSERT(c == number_of_cols());
    }

    // For uniformity of interface, the first arg does nothing
    StaticMatrix(void const* ptr, std::initializer_list<scalar_type> const& c)
        : StaticMatrix(c) {
      (void) ptr;
      LIBSEMIGROUPS_ASSERT(ptr == nullptr);
    }

    // For uniformity of interface, the first arg does nothing
    StaticMatrix(
        void const*                                                      ptr,
        std::initializer_list<std::initializer_list<scalar_type>> const& m)
        : StaticMatrix(m) {
      (void) ptr;
      LIBSEMIGROUPS_ASSERT(ptr == nullptr);
    }

    // For uniformity of interface, the first arg does nothing
    explicit StaticMatrix(void const* ptr, RowView const& rv)
        : StaticMatrix(rv) {
      (void) ptr;
      LIBSEMIGROUPS_ASSERT(ptr == nullptr);
    }

    // For uniformity of interface, no arg used for anything
    StaticMatrix(void const* ptr, size_t r, size_t c) : StaticMatrix(r, c) {
      (void) ptr;
      LIBSEMIGROUPS_ASSERT(ptr == nullptr);
    }
#endif

    ~StaticMatrix() = default;

#ifndef PARSED_BY_DOXYGEN
    static StaticMatrix one(size_t n) {
      // If specified the value of n must equal R or otherwise weirdness will
      // ensue...
      LIBSEMIGROUPS_ASSERT(n == 0 || n == R);
      (void) n;
#if defined(__APPLE__) && defined(__clang__) \
    && (__clang_major__ == 13 || __clang_major__ == 14)
      // With Apple clang version 13.1.6 (clang-1316.0.21.2.5) something goes
      // wrong and the value R is optimized away somehow, meaning that the
      // values on the diagonal aren't actually set. This only occurs when
      // libsemigroups is compiled with -O2 or higher.
      size_t volatile const m = R;
#else
      size_t const m = R;
#endif
      StaticMatrix x(m, m);
      std::fill(x.begin(), x.end(), ZeroOp()());
      for (size_t r = 0; r < m; ++r) {
        x(r, r) = OneOp()();
      }
      return x;
    }

    static StaticMatrix one(void const* ptr, size_t n = 0) {
      (void) ptr;
      LIBSEMIGROUPS_ASSERT(ptr == nullptr);
      LIBSEMIGROUPS_ASSERT(n == 0 || n == R);
      return one(n);
    }
#endif  // PARSED_BY_DOXYGEN

    ////////////////////////////////////////////////////////////////////////
    // StaticMatrix - member function aliases - public
    ////////////////////////////////////////////////////////////////////////
#ifdef PARSED_BY_DOXYGEN
    //! \brief Returns a reference to the specified entry of the matrix.
    //!
    //! Returns a reference to the specified entry of the matrix.
    //!
    //! \param r  the index of the row of the entry.
    //! \param c  the index of the column of the entry.
    //!
    //! \returns  A value of type \ref scalar_reference.
    //!
    //! \exceptions
    //!   this function guarantees not to throw a LibsemigroupsException.
    //!
    //! \complexity
    //!   Constant
    //!
    //! \warning
    //!    No checks on the validity of the parameters \c r and \c c are
    //!    performed.
    scalar_reference operator()(size_t r, size_t c);

    //! \brief Returns a reference to the specified entry of the matrix.
    //!
    //! Returns a reference to the specified entry of the matrix.
    //!
    //! \param r  the index of the row of the entry.
    //! \param c  the index of the column of the entry.
    //!
    //! \returns  A value of type \ref scalar_reference.
    //!
    //! \throws LibsemigroupsException if \p r or \p c is out of bounds.
    //!
    //! \complexity
    //!   Constant
    scalar_reference at(size_t r, size_t c);

    //! \brief Returns a const reference to the specified entry of the matrix.
    //!
    //! Returns a const reference to the specified entry of the matrix.
    //!
    //! \param r  the index of the row of the entry.
    //! \param c  the index of the column of the entry.
    //!
    //! \returns  A value of type \ref scalar_const_reference.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant
    //!
    //! \warning
    //! No checks on the validity of the parameters \c r and \c c are
    //! performed.
    scalar_const_reference operator()(size_t r, size_t c) const;

    //! \brief Returns a const reference to the specified entry of the matrix.
    //!
    //! Returns a const reference to the specified entry of the matrix.
    //!
    //! \param r  the index of the row of the entry.
    //! \param c  the index of the column of the entry.
    //!
    //! \returns  A value of type \ref scalar_const_reference.
    //!
    //! \throws LibsemigroupsException if \p r or \p c is out of bounds.
    //!
    //! \complexity
    //! Constant
    scalar_const_reference at(size_t r, size_t c) const;

    //! \brief Returns an iterator pointing at the first entry.
    //!
    //! This function returns a random access iterator point at the first entry
    //! of the matrix.
    //!
    //! \returns A value of type \ref StaticMatrix::iterator "iterator".
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant
    //!
    //! \warning
    //! The order in which entries in the matrix are iterated over is not
    //! specified.
    iterator begin() noexcept;

    //! \brief Returns an iterator pointing one beyond the last entry in the
    //! matrix.
    //!
    //! \returns A value of type \ref StaticMatrix::iterator "iterator".
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant
    //!
    //! \warning
    //! The order in which entries in the matrix are iterated over is not
    //! specified.
    iterator end();

    //! \brief Returns a const iterator pointing at the first entry.
    //!
    //! This function returns a const (random access) iterator pointing at the
    //! first entry in the matrix.
    //!
    //! \returns  A value of type \ref StaticMatrix::const_iterator
    //! "const_iterator".
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant
    //!
    //! \warning
    //!   The order in which entries in the matrix are iterated over is not
    //!   specified.
    const_iterator cbegin() const noexcept;

    //! \brief Returns a const iterator pointing one beyond the last entry in
    //! the matrix.
    //!
    //! This function returns a const (random access) iterator pointing one
    //! passed the last entry of the matrix.
    //!
    //! \returns  A value of type \ref StaticMatrix::const_iterator
    //! "const_iterator".
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant
    //!
    //! \warning
    //! The order in which entries in the matrix are iterated over is not
    //! specified.
    const_iterator cend();

    //! \brief Equality operator.
    //!
    //! This function returns \c true if `*this` is equal to \p that;
    //! and \c false otherwise.
    //!
    //! \param that  matrix or row view for comparison.
    //!
    //! \returns \c true if `*this` and `that` are equal and \c false if they
    //! are not.
    //!
    //! \complexity
    //! At worst \f$O(mn)\f$ where \f$m\f$ is the number of rows and
    //! \f$n\f$ is the number of columns of the matrix.
    bool operator==(StaticMatrix const& that) const;

    //! \copydoc operator==
    bool operator==(RowView const& that) const;

    //! \brief Inequality operator.
    //!
    //! This function returns \c true if `*this` is not equal to \p that; and
    //! \c false otherwise.
    //!
    //! \param that  the matrix or row view for comparison.
    //!
    //! \returns The negation of `operator==(that)`.
    //!
    //! \complexity
    //! see operator==
    bool operator!=(StaticMatrix const& that) const;

    //! \copydoc operator!=
    bool operator!=(RowView const& that) const;

    //! \brief Less than operator.
    //!
    //! This operator defines a total order on the set of matrices of the
    //! same type.
    //!
    //! \param that the matrix or row view for comparison.
    //!
    //! \returns \c true if `*this` is less than \p that  and \c false if it is
    //! not.
    //!
    //! \complexity
    //! At worst \f$O(mn)\f$ where \f$m\f$ is \ref number_of_rows
    //! and \f$n\f$ is \ref number_of_cols.
    bool operator<(StaticMatrix const& that) const;

    //! \copydoc operator<
    bool operator<(RowView const& that) const;

    //! \brief Greater than operator.
    //!
    //! This operator defines a total order on the set of matrices of the same
    //! type.
    //!
    //! \param that  the matrix for comparison.
    //!
    //! \returns \c true if `*this` is less than `that` and \c false if it is
    //! not.
    //!
    //! \complexity
    //! At worst \f$O(mn)\f$ where \f$m\f$ is \ref number_of_rows
    //! and \f$m\f$ is \ref number_of_cols
    bool operator>(StaticMatrix const& that) const;

    //! \brief Get the coordinates of an iterator.
    //!
    //! This function returns a pair containing the row and columns
    //! corresponding to a \ref StaticMatrix::const_iterator "const_iterator"
    //! pointing into a matrix.
    //!
    //! \param it the \ref StaticMatrix::const_iterator "const_iterator".
    //!
    //! \returns A value of type `std::pair<scalar_type, scalar_type>`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant
    std::pair<scalar_type, scalar_type> coords(const_iterator it) const;

    //! \brief Returns the number of rows of the matrix.
    //!
    //! This function returns the number of rows of the matrix.
    //!
    //! \returns A value of type `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant
    size_t number_of_rows() const noexcept;

    //! \brief Returns the number of columns of the matrix.
    //!
    //! This function returns the number of columns of the matrix.
    //!
    //! \returns A value of type `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant
    size_t number_of_cols() const noexcept;

    //! \brief Returns the sum of two matrices.
    //!
    //! This function returns the sum of two matrices.
    //!
    //! \param that  the matrix to add to `this`.
    //!
    //! \returns The sum of the two matrices.
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is \ref number_of_rows
    //! and \f$m\f$ \ref number_of_cols
    //!
    //! \warning
    //! The matrices must be of the same dimensions, although this is not
    //! verified.
    //!
    //! \warning This function does not detect overflows of \ref scalar_type.
    StaticMatrix operator+(StaticMatrix const& that);

    //! \brief Add a matrix to another matrix in-place.
    //!
    //! This function adds a matrix (or the row represented by a RowView) to
    //! another matrix of the same shape in-place.
    //!
    //! \param that  the matrix or row view to add.
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is \ref number_of_rows and \f$m\f$ is \ref
    //! number_of_cols
    //!
    //! \warning
    //! The matrices must be of the same dimensions, although this is not
    //! checked.
    //!
    //! \warning This function does not detect overflows of \ref scalar_type.
    void operator+=(StaticMatrix const& that);

    //! \copydoc operator+=
    void operator+=(RowView const& that);

    //! \brief Adds a scalar to every entry of the matrix in-place.
    //!
    //! This function adds a scalar to every entry of the matrix in-place.
    //!
    //! \param a the \ref scalar_type to add to a \ref StaticMatrix.
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is \ref number_of_rows and \f$m\f$ is \ref
    //! number_of_cols
    //!
    //! \warning This function does not detect overflows of \ref scalar_type.
    void operator+=(scalar_type a);

    //! \brief Returns the product of two matrices.
    //!
    //! This function returns the product of two matrices.
    //!
    //! \param that the matrix to multiply by `this`.
    //!
    //! \returns The product of the two matrices.
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is \ref number_of_rows and \f$m\f$ is \ref
    //! number_of_cols.
    //!
    //! \warning
    //! The matrices must be of the same dimensions, although this is not
    //! verified.
    //!
    //! \warning This function does not detect overflows of \ref scalar_type.
    StaticMatrix operator*(StaticMatrix const& that);

    //! \brief Multiplies every entry of a matrix by a scalar in-place.
    //!
    //! This function multiplies every entry of a matrix by a scalar in-place.
    //!
    //! \param a the scalar to multiply every entry by.
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is \ref number_of_rows
    //! and \f$m\f$ is \ref number_of_cols.
    //!
    //! \warning This function does not detect overflows of \ref scalar_type.
    void operator*=(scalar_type a);

    //! \brief Multiplies \p x and \p y and stores the result in `this`.
    //!
    //! This function redefines `this` to be the product of \p x and \p y. This
    //! is in-place multiplication to avoid allocation of memory for products
    //! which do not need to be stored for future use.
    //!
    //! \param x  the first matrix to multiply.
    //! \param y  the second matrix to multiply.
    //!
    //! \complexity
    //! \f$O(n ^ 3)\f$ where \f$n\f$ is \ref number_of_rows or \ref
    //! number_of_cols.
    //!
    //! \warning
    //! This function only applies to matrices with the same number of rows
    //! and columns but this isn't verified.
    void product_inplace_no_checks(StaticMatrix const& x,
                                   StaticMatrix const& y);

    //! \brief Returns a view into a row.
    //!
    //! This function returns a \ref RowView into the row with index \p i.
    //!
    //! \param i  the index of the row.
    //!
    //! \returns  A value of type \ref RowView.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant
    //!
    //! \warning No checks are made on the argument \p i.
    RowView row_no_checks(size_t i) const;

    //! \brief Returns a view into a row.
    //!
    //! This function returns a \ref RowView into the row with index \p i.
    //!
    //! \param i the index of the row.
    //!
    //! \returns A value of type \ref RowView.
    //!
    //! \throws LibsemigroupsException if \p i is greater than or equal
    //! to \ref number_of_rows.
    RowView row(size_t i) const;

    //! \brief Add row views for every row in the matrix to a container.
    //!
    //! This function adds a \ref RowView for each row in the matrix to the
    //! container \p x.
    //!
    //! \tparam T the type of the container.
    //! \param x  a container.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(m)\f$ where \f$m\f$ is the template parameter \c R.
    template <typename T>
    void rows(T& x) const;

    //! \brief Swaps the contents of `*this` with the contents of \p that.
    //!
    //! This function swaps the contents of `*this` with the contents of \p
    //! that.
    //!
    //! \param that  the matrix to swap contents with.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant
    void swap(StaticMatrix& that) noexcept;

    //! \brief Transpose a matrix in-place.
    //!
    //! This function transpose a matrix object in-place.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \throws LibsemigroupsException if the matrix isn't square.
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is \ref number_of_rows and
    //! \f$n\f$ is \ref number_of_cols.
    void transpose();

    //! \brief Transpose a matrix in-place.
    //!
    //! This function transpose a matrix object in-place.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is \ref number_of_rows and
    //! \f$n\f$ is \ref number_of_cols.
    //!
    //! \warning
    //! This only works when \ref number_of_rows and \ref number_of_cols are
    //! equal (i.e. for square matrices), but this is not verified.
    void transpose_no_checks();

    //! \brief Returns an identity matrix.
    //!
    //! This function returns an \c R times \c R identity matrix.
    //!
    //! \returns  The identity matrix with \c R rows and \c R columns.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(n ^ 2)\f$ where \f$n\f$ is the template parameter \c R.
    static StaticMatrix one() const;

    //! \brief Return a hash value of a matrix.
    //!
    //! This function returns a hash value for a matrix. The return value is
    //! recomputed every time this function is called.
    //!
    //! \returns A hash value for a \c this.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is the \ref number_of_rows and \f$n\f$ is
    //! \ref number_of_cols.
    size_t hash_value() const;

    //! \brief Less than or equal operator.
    //!
    //! This operator defines a total order on the set of matrices of the
    //! same type.
    //!
    //! \tparam U the type of the argument \p that.
    //! \param that the matrix or row view for comparison.
    //!
    //! \returns \c true if `*this` is less than or equal to \p that  and \c
    //! false if it is not.
    //!
    //! \complexity
    //! At worst \f$O(mn)\f$ where \f$m\f$ is \ref number_of_rows
    //! and \f$n\f$ is \ref number_of_cols.
    template <typename U>
    bool operator<=(U const& that) const;

    //! \brief Greater than or equal operator.
    //!
    //! This operator defines a total order on the set of matrices of the
    //! same type.
    //!
    //! \tparam U the type of the argument \p that.
    //! \param that the matrix or row view for comparison.
    //!
    //! \returns \c true if `*this` is greater than or equal to \p that  and \c
    //! false if it is not.
    //!
    //! \complexity
    //! At worst \f$O(mn)\f$ where \f$m\f$ is \ref number_of_rows
    //! and \f$n\f$ is \ref number_of_cols.
    template <typename U>
    bool operator>=(U const& that) const;

    //! \brief Multiply a matrix by a scalar.
    //!
    //! This function returns the product of the matrix `*this` and the scalar
    //! \p a.
    //!
    //! \param a the scalar to add.
    //! \returns  The product of `*this` and \p a.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning This function does not detect overflows in \ref scalar_type.
    StaticMatrix operator*(scalar_type a);

    //! \brief Add a scalar to a matrix.
    //!
    //! This function returns the sum of the matrix `*this` and the scalar
    //! \p a.
    //!
    //! \param a the scalar to add.
    //! \returns  The sum of `*this` and \p a.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning This function does not detect overflows in \ref
    //! scalar_type.
    StaticMatrix operator+(scalar_type a);

    //! Returns the additive identity of the underlying semiring.
    //!
    //! This function returns the additive identity of the underlying
    //! semiring of a matrix.
    //!
    //! \returns The additive identity of the semiring, a \ref
    //! scalar_type.
    //!
    //! \exceptions
    //! \noexcept
    scalar_type scalar_zero() const noexcept;

    //! Returns the multiplicative identity of the underlying semiring.
    //!
    //! This function returns the multiplicative identity of the underlying
    //! semiring of a matrix.
    //!
    //! \returns The multiplicative identity of the semiring, a \ref
    //! scalar_type.
    //!
    //! \exceptions
    //! \noexcept
    scalar_type scalar_one() const noexcept;

    //! \brief Returns the underlying semiring.
    //!
    //! Returns a const pointer to the underlying semiring (if any).
    //!
    //! \returns  A value of type `Semiring const*`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant
    semiring_type const* semiring() const noexcept;

#else
    using MatrixCommon::at;
    using MatrixCommon::begin;
    using MatrixCommon::cbegin;
    using MatrixCommon::cend;
    using MatrixCommon::coords;
    using MatrixCommon::end;
    using MatrixCommon::hash_value;
    using MatrixCommon::number_of_cols;
    using MatrixCommon::number_of_rows;
    using MatrixCommon::one;
    using MatrixCommon::operator();
    using MatrixCommon::operator==;
    using MatrixCommon::operator<;  // NOLINT(whitespace/operators)
    using MatrixCommon::operator<=;
    using MatrixCommon::operator>;  // NOLINT(whitespace/operators)
    using MatrixCommon::operator>=;
    using MatrixCommon::operator!=;
    using MatrixCommon::operator*=;
    using MatrixCommon::operator+=;
    using MatrixCommon::operator*;
    using MatrixCommon::operator+;
    using MatrixCommon::product_inplace_no_checks;
    using MatrixCommon::row;
    using MatrixCommon::row_no_checks;
    using MatrixCommon::rows;
    using MatrixCommon::scalar_one;
    using MatrixCommon::scalar_zero;
    using MatrixCommon::semiring;
    using MatrixCommon::swap;
    using MatrixCommon::transpose;
    using MatrixCommon::transpose_no_checks;
#endif  // PARSED_BY_DOXYGEN

   private:
    ////////////////////////////////////////////////////////////////////////
    // StaticMatrix - implementation of MatrixCommon requirements - private
    ////////////////////////////////////////////////////////////////////////

    static constexpr size_t number_of_rows_impl() noexcept {
      return R;
    }
    static constexpr size_t number_of_cols_impl() noexcept {
      return C;
    }
  };  // namespace libsemigroups

  ////////////////////////////////////////////////////////////////////////
  // DynamicMatrix with compile time semiring arithmetic
  ////////////////////////////////////////////////////////////////////////

  //! \anchor DynamicMatrixStaticArith
  //!
  //! \brief Class for matrices with compile time arithmetic and run-time
  //! dimensions.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a class for matrices where the arithmetic operations in the
  //! underlying semiring are known at compile time, and the dimensions of the
  //! matrix can be set at run time.
  //!
  //! \tparam PlusOp a stateless type with a call operator of signature
  //! `scalar_type operator()(scalar_type, scalar_type)` implementing the
  //! addition of the semiring
  //!
  //! \tparam ProdOp a stateless type with a call operator of signature
  //! `scalar_type operator()(scalar_type, scalar_type)` implementing the
  //! multiplication of the semiring
  //!
  //! \tparam ZeroOp a stateless type with a call operator of signature
  //! `scalar_type operator()()` returning the zero of the semiring (the
  //! additive identity element)
  //!
  //! \tparam OneOp a stateless type with a call operator of signature
  //! `scalar_type operator()()` returning the one of the semiring (the
  //! multiplicative identity element)
  //!
  //! \tparam Scalar the type of the entries in the matrices (the type of
  //! elements in the underlying semiring)
  //!
  //! \note
  //! Certain member functions only work for square matrices and some only
  //! work for rows.
  template <typename PlusOp,
            typename ProdOp,
            typename ZeroOp,
            typename OneOp,
            typename Scalar>
  class DynamicMatrix<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>
      : public detail::MatrixDynamicDim<Scalar>,
        public detail::MatrixCommon<
            std::vector<Scalar>,
            DynamicMatrix<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>,
            DynamicRowView<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>>,
        public detail::
            MatrixStaticArithmetic<PlusOp, ProdOp, ZeroOp, OneOp, Scalar> {
    using MatrixDynamicDim = ::libsemigroups::detail::MatrixDynamicDim<Scalar>;
    using MatrixCommon     = ::libsemigroups::detail::MatrixCommon<
        std::vector<Scalar>,
        DynamicMatrix<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>,
        DynamicRowView<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>>;
    friend MatrixCommon;

   public:
    //! The type of the entries in the matrix.
    using scalar_type = typename MatrixCommon::scalar_type;

    //! The type of references to the entries in the matrix.
    using scalar_reference = typename MatrixCommon::scalar_reference;

    //! The type of const references to the entries in the matrix.
    using scalar_const_reference =
        typename MatrixCommon::scalar_const_reference;

    //! The type of a row of a DynamicMatrix.
    using Row = DynamicMatrix;

    //! The type of a row view into a DynamicMatrix.
    using RowView = DynamicRowView<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>;

    //! Alias for the template parameter \p PlusOp.
    using Plus = PlusOp;

    //! Alias for the template parameter \p ProdOp.
    using Prod = ProdOp;

    //! Alias for the template parameter \p ZeroOp.
    using Zero = ZeroOp;

    //! Alias for the template parameter \p OneOp.
    using One = OneOp;

    //! \brief Alias for the semiring type (`void`).
    //!
    //! The type of the semiring over which the matrix is defined is `void`
    //! because there's no semiring object, the arithmetic is defined by the
    //! template parameters.
    using semiring_type = void;

    //! \brief Default constructor.
    //!
    //! Default constructor.
    DynamicMatrix() = default;

    //! \brief Default copy constructor.
    //!
    //! Default copy constructor.
    DynamicMatrix(DynamicMatrix const&) = default;

    //! \brief Default move constructor.
    //!
    //! Default move constructor.
    DynamicMatrix(DynamicMatrix&&) = default;

    //! \brief Default copy assignment operator.
    //!
    //! Default copy assignment operator.
    DynamicMatrix& operator=(DynamicMatrix const&) = default;

    //! \brief Default move assignment operator.
    //!
    //! Default move assignment operator.
    DynamicMatrix& operator=(DynamicMatrix&&) = default;

    //! \brief Construct a matrix with given dimensions.
    //!
    //! This function constructs a matrix with the given dimensions.
    //!
    //! \param r the number of rows.
    //! \param c the number of columns.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant
    //!
    //! \par Example
    //! \code
    //! Mat m(2, 3); // construct a 2 x 3 matrix
    //! \endcode
    DynamicMatrix(size_t r, size_t c) : MatrixDynamicDim(r, c), MatrixCommon() {
      resize(number_of_rows(), number_of_cols());
    }

    //! \brief Construct a vector from a std::initializer_list.
    //!
    //! This function constructs a \f$1\times n\f$ matrix from a
    //! std::initializer_list of the entries of the row.
    //!
    //! \param c  the values to be copied into the matrix.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(n)\f$ where \f$n\f$ is the number of columns in the matrix being
    //! constructed.
    //!
    //! \par Example
    //! \code
    //! Mat m({1, 1, 0, 0});
    //! \endcode
    explicit DynamicMatrix(std::initializer_list<scalar_type> const& c)
        : MatrixDynamicDim(1, c.size()), MatrixCommon(c) {}

    //! \brief Construct a matrix from std::initializer_list of
    //! std::initializer_list of scalars.
    //!
    //! This function constructs a matrix from a std::initializer_list of
    //! std::initializer_list of the entries.
    //!
    //! \param m  the values to be copied into the matrix.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is the number of rows and \f$n\f$ is the
    //! number of columns in the matrix being constructed.
    //!
    //! \par Example
    //! \code
    //! Mat m({{1, 1}, {0, 0}});
    //! \endcode
    explicit DynamicMatrix(
        std::initializer_list<std::initializer_list<scalar_type>> const& m)
        : MatrixDynamicDim(m.size(), m.begin()->size()), MatrixCommon(m) {}

    //! \brief Construct a matrix from std::vector of
    //! std::vector of scalars.
    //!
    //! This function constructs a matrix from a std::vector of
    //! std::vector of the entries.
    //!
    //! \param m  the values to be copied into the matrix.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is the number of rows and \f$n\f$ is the
    //! number of columns in the matrix being constructed.
    explicit DynamicMatrix(std::vector<std::vector<scalar_type>> const& m)
        : MatrixDynamicDim(m.size(), m.begin()->size()), MatrixCommon(m) {}

    //! \brief Construct a row from a row view.
    //!
    //! Construct a row from a row view.
    //!
    //! \param rv the row view.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(n)\f$ where \f$n\f$ is the size of the row being constructed.
    explicit DynamicMatrix(RowView const& rv)
        : MatrixDynamicDim(1, rv.size()), MatrixCommon(rv) {}

#ifndef PARSED_BY_DOXYGEN
    // For uniformity of interface, the first arg does nothing
    DynamicMatrix(void const* ptr, size_t r, size_t c) : DynamicMatrix(r, c) {
      (void) ptr;
      LIBSEMIGROUPS_ASSERT(ptr == nullptr);
    }

    // For uniformity of interface, the first arg does nothing
    DynamicMatrix(void const* ptr, std::initializer_list<scalar_type> const& c)
        : DynamicMatrix(c) {
      (void) ptr;
      LIBSEMIGROUPS_ASSERT(ptr == nullptr);
    }

    // For uniformity of interface, the first arg does nothing
    DynamicMatrix(
        void const*                                                      ptr,
        std::initializer_list<std::initializer_list<scalar_type>> const& m)
        : DynamicMatrix(m) {
      (void) ptr;
      LIBSEMIGROUPS_ASSERT(ptr == nullptr);
    }

    static DynamicMatrix one(void const* ptr, size_t n) {
      (void) ptr;
      LIBSEMIGROUPS_ASSERT(ptr == nullptr);
      return one(n);
    }
#endif  // PARSED_BY_DOXYGEN

    ~DynamicMatrix();

    //! \brief Construct the \f$n \times n\f$ identity matrix.
    //!
    //! Construct the \f$n \times n\f$ identity matrix.
    //!
    //! \param n the dimension.
    //! \returns The \f$n \times n\f$ identity matrix.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(n ^ 2)\f$.
    static DynamicMatrix one(size_t n) {
      DynamicMatrix x(n, n);
      std::fill(x.begin(), x.end(), ZeroOp()());
      for (size_t r = 0; r < n; ++r) {
        x(r, r) = OneOp()();
      }
      return x;
    }

#ifdef PARSED_BY_DOXYGEN
    //! \copydoc StaticMatrix::operator()(size_t, size_t)
    scalar_reference operator()(size_t r, size_t c);

    //! \copydoc StaticMatrix::operator()(size_t, size_t) const
    scalar_const_reference operator()(size_t r, size_t c) const;

    //! \copydoc StaticMatrix::begin()
    iterator begin() noexcept;

    //! \copydoc StaticMatrix::end()
    iterator end() noexcept;

    //! \copydoc StaticMatrix::cbegin()
    const_iterator cbegin() noexcept;

    //! \copydoc StaticMatrix::cend()
    const_iterator cend() noexcept;

    //! \copydoc StaticMatrix::operator==
    bool operator==(DynamicMatrix const& that) const;

    //! \copydoc StaticMatrix::operator==
    bool operator==(RowView const& that) const;

    //! \copydoc StaticMatrix::operator!=
    bool operator!=(DynamicMatrix const& that) const;

    //! \copydoc StaticMatrix::operator!=
    bool operator!=(RowView const& that) const;

    //! \copydoc StaticMatrix::operator<
    bool operator<(DynamicMatrix const& that) const;

    //! \copydoc StaticMatrix::operator<
    bool operator<(RowView const& that) const;

    //! \copydoc StaticMatrix::operator>
    bool operator>(DynamicMatrix const& that) const;

    //! \copydoc StaticMatrix::coords
    std::pair<scalar_type, scalar_type> coords(const_iterator it) const;

    //! \copydoc StaticMatrix::number_of_rows
    size_t number_of_rows() const noexcept;

    //! \copydoc StaticMatrix::number_of_cols
    size_t number_of_cols() const noexcept;

    //! \copydoc StaticMatrix::operator+
    DynamicMatrix operator+(DynamicMatrix const& that);

    //! \copydoc StaticMatrix::operator+=
    void operator+=(DynamicMatrix const& that);

    //! \copydoc StaticMatrix::operator+=
    void operator+=(RowView const& that);

    //! \copydoc StaticMatrix::operator+=(scalar_type)
    void operator+=(scalar_type a);

    //! \copydoc StaticMatrix::operator*
    DynamicMatrix operator*(DynamicMatrix const& that);

    //! \copydoc StaticMatrix::operator*=
    void operator*=(scalar_type a);

    //! \copydoc StaticMatrix::product_inplace_no_checks
    void product_inplace_no_checks(DynamicMatrix const& x,
                                   DynamicMatrix const& y);

    //! \copydoc StaticMatrix::row_no_checks
    RowView row_no_checks(size_t i) const;

    //! \copydoc StaticMatrix::row
    RowView row(size_t i) const;

    //! \copydoc StaticMatrix::rows
    template <typename T>
    void rows(T& x) const;

    //! \copydoc StaticMatrix::transpose_no_checks
    void transpose_no_checks();

    //! \copydoc StaticMatrix::transpose
    void transpose();

    //! \copydoc StaticMatrix::hash_value
    size_t hash_value() const;

    //! \copydoc StaticMatrix::operator<=
    bool operator<=(T const& that) const;

    //! \copydoc StaticMatrix::operator>=
    bool operator>=(T const& that) const;

    //! \copydoc StaticMatrix::operator*(scalar_type)
    StaticMatrix operator*(scalar_type a);

    //! \copydoc StaticMatrix::operator+(scalar_type)
    StaticMatrix operator+(scalar_type a);

    //! \copydoc StaticMatrix::scalar_zero
    scalar_type scalar_zero() const noexcept;

    //! \copydoc StaticMatrix::scalar_one
    scalar_type scalar_one() const noexcept;

    //! \copydoc StaticMatrix::semiring
    semiring_type const* semiring() const noexcept;
#else
    using MatrixCommon::at;
    using MatrixCommon::begin;
    using MatrixCommon::cbegin;
    using MatrixCommon::cend;
    using MatrixCommon::coords;
    using MatrixCommon::end;
    using MatrixCommon::hash_value;
    using MatrixCommon::number_of_cols;
    using MatrixCommon::number_of_rows;
    using MatrixCommon::one;
    using MatrixCommon::operator();
    using MatrixCommon::operator==;
    using MatrixCommon::operator<;  // NOLINT(whitespace/operators)
    using MatrixCommon::operator<=;
    using MatrixCommon::operator>;  // NOLINT(whitespace/operators)
    using MatrixCommon::operator>=;
    using MatrixCommon::operator!=;
    using MatrixCommon::operator*=;
    using MatrixCommon::operator+=;
    using MatrixCommon::operator*;
    using MatrixCommon::operator+;
    using MatrixCommon::product_inplace_no_checks;
    using MatrixCommon::row;
    using MatrixCommon::row_no_checks;
    using MatrixCommon::rows;
    using MatrixCommon::scalar_one;
    using MatrixCommon::scalar_zero;
    using MatrixCommon::semiring;
    // using MatrixCommon::swap; don't want this, see below
    using MatrixCommon::transpose;
    using MatrixCommon::transpose_no_checks;
#endif  // PARSED_BY_DOXYGEN

    //! \copydoc StaticMatrix::swap
    void swap(DynamicMatrix& that) noexcept {
      static_cast<MatrixDynamicDim&>(*this).swap(
          static_cast<MatrixDynamicDim&>(that));
      static_cast<MatrixCommon&>(*this).swap(static_cast<MatrixCommon&>(that));
    }

   private:
    using MatrixCommon::resize;
  };

  template <typename PlusOp,
            typename ProdOp,
            typename ZeroOp,
            typename OneOp,
            typename Scalar>
  DynamicMatrix<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>::~DynamicMatrix()
      = default;

  ////////////////////////////////////////////////////////////////////////
  // DynamicMatrix with runtime semiring arithmetic
  ////////////////////////////////////////////////////////////////////////

  //! \anchor DynamicMatrixDynamicArith
  //!
  //! \brief Class for matrices with run-time arithmetic and dimensions.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a class for matrices where both the arithmetic operations in
  //! the underlying semiring and the dimensions of the matrix can be set at
  //! run time.
  //!
  //! \tparam Semiring the type of a semiring object which defines the semiring
  //! arithmetic (see requirements below).
  //!
  //! \tparam Scalar the type of the entries in the matrices (the type of
  //! elements in the underlying semiring)
  //!
  //! \note
  //! Certain member functions only work for square matrices and some only
  //! work for rows.
  //!
  //! ## Semiring requirements
  //!
  //! The template parameter Semiring must have the following
  //! member functions:
  //! * `scalar_type scalar_zero()` that returns the multiplicative zero scalar
  //! in the semiring
  //! * `scalar_type scalar_one()` that returns the multiplicative identity
  //! scalar in the semiring
  //! * `scalar_type plus_no_checks(scalar_type x, scalar_type y)` that returns
  //! the sum in the semiring of the scalars \c x and \c y.
  //! * `scalar_type product_no_checks(scalar_type x, scalar_type y)` that
  //! returns the product in the semiring of the scalars \c x and \c y.
  //!
  //! See, for example, \ref MaxPlusTruncSemiring.
  template <typename Semiring, typename Scalar>
  class DynamicMatrix<Semiring, Scalar>
      : public detail::MatrixDynamicDim<Scalar>,
        public detail::MatrixCommon<std::vector<Scalar>,
                                    DynamicMatrix<Semiring, Scalar>,
                                    DynamicRowView<Semiring, Scalar>,
                                    Semiring> {
    using MatrixCommon = detail::MatrixCommon<std::vector<Scalar>,
                                              DynamicMatrix<Semiring, Scalar>,
                                              DynamicRowView<Semiring, Scalar>,
                                              Semiring>;
    friend MatrixCommon;
    using MatrixDynamicDim = ::libsemigroups::detail::MatrixDynamicDim<Scalar>;

   public:
    //! \copydoc StaticMatrix::scalar_type
    using scalar_type = typename MatrixCommon::scalar_type;
    //! \copydoc StaticMatrix::scalar_reference
    using scalar_reference = typename MatrixCommon::scalar_reference;
    //! \copydoc StaticMatrix::scalar_const_reference
    using scalar_const_reference =
        typename MatrixCommon::scalar_const_reference;

    //! \copydoc StaticMatrix::Row
    using Row = DynamicMatrix;

    //! \copydoc StaticMatrix::RowView
    using RowView = DynamicRowView<Semiring, Scalar>;
    friend RowView;

    //! \brief Alias for the template parameter Semiring.
    using semiring_type = Semiring;

    //! \brief Deleted.
    //!
    //! The default constructor for this variant of DynamicMatrix is deleted
    //! because a valid semiring object is required to define the arithmetic at
    //! run-time.
    DynamicMatrix() = delete;

    //! \copydoc StaticMatrix::StaticMatrix(StaticMatrix const&)
    DynamicMatrix(DynamicMatrix const&) = default;

    //! \copydoc StaticMatrix::StaticMatrix(StaticMatrix&&)
    DynamicMatrix(DynamicMatrix&&) = default;

    //! \copydoc StaticMatrix::operator=(StaticMatrix const&)
    DynamicMatrix& operator=(DynamicMatrix const&) = default;

    //! \copydoc StaticMatrix::operator=(StaticMatrix&&)
    DynamicMatrix& operator=(DynamicMatrix&&) = default;

    //! \brief Construct a matrix over a given semiring with given dimensions.
    //!
    //! Construct a matrix over the semiring `semiring` with the given
    //! dimensions.
    //!
    //! \param semiring  a pointer to const semiring object.
    //! \param r  the number of rows in the matrix being constructed.
    //! \param c  the number of columns in the matrix being constructed.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant
    DynamicMatrix(Semiring const* semiring, size_t r, size_t c)
        : MatrixDynamicDim(r, c), MatrixCommon(), _semiring(semiring) {
      resize(number_of_rows(), number_of_cols());
    }

    //! \brief Construct a matrix over a given semiring (std::initializer_list
    //! of std::initializer_list).
    //!
    //! This function constructs a matrix over a given semiring from an
    //! std::initializer_list of the rows (std::initializer_list).
    //!
    //! \param semiring  a pointer to const semiring object.
    //! \param rows  the values to be copied into the matrix.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is \ref number_of_rows and \f$n\f$ is \ref
    //! number_of_cols.
    explicit DynamicMatrix(
        Semiring const* semiring,
        std::initializer_list<std::initializer_list<scalar_type>> const& rows)
        : MatrixDynamicDim(rows.size(), rows.begin()->size()),
          MatrixCommon(rows),
          _semiring(semiring) {}

    //! \brief Construct a matrix over a given semiring (std::vector
    //! of std::vector).
    //!
    //! This function constructs a matrix over a given semiring from an
    //! std::vector of the rows (std::vector).
    //!
    //! \param semiring  a pointer to const semiring object.
    //! \param rows  the values to be copied into the matrix.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is \ref number_of_rows and \f$n\f$ is \ref
    //! number_of_cols.
    explicit DynamicMatrix(Semiring const* semiring,
                           std::vector<std::vector<scalar_type>> const& rows)
        : MatrixDynamicDim(rows.size(),
                           (rows.empty() ? 0 : rows.begin()->size())),
          MatrixCommon(rows),
          _semiring(semiring) {}

    //! \brief Construct a row over a given semiring (std::initializer_list).
    //!
    //! Construct a row over a given semiring from a std::initializer_list of
    //! the entries.
    //!
    //! \param semiring a pointer to const Semiring object.
    //! \param row the values to be copied into the row.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is \ref number_of_cols.
    explicit DynamicMatrix(Semiring const*                           semiring,
                           std::initializer_list<scalar_type> const& row)
        : MatrixDynamicDim(1, row.size()),
          MatrixCommon(row),
          _semiring(semiring) {}

    //! \brief Construct a row over a given semiring (\ref RowView).
    //!
    //! Construct a row over a given semiring from a \ref RowView.
    //!
    //! \param rv  the row view.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(n)\f$ where \f$n\f$ is the size of the row being constructed.
    explicit DynamicMatrix(RowView const& rv)
        : MatrixDynamicDim(1, rv.size()),
          MatrixCommon(rv),
          _semiring(rv._matrix->semiring()) {}

    //! \brief Construct the \f$n \times n\f$ identity matrix.
    //!
    //! Construct the \f$n \times n\f$ identity matrix.
    //!
    //! \param semiring the semiring.
    //! \param n the dimension.
    //!
    //! \returns The \f$n \times n\f$ identity matrix over \p semiring.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(n ^ 2)\f$.
    // No static DynamicMatrix::one(size_t n) because we need a semiring!
    static DynamicMatrix one(Semiring const* semiring, size_t n) {
      DynamicMatrix x(semiring, n, n);
      std::fill(x.begin(), x.end(), x.scalar_zero());
      for (size_t r = 0; r < n; ++r) {
        x(r, r) = x.scalar_one();
      }
      return x;
    }

    ~DynamicMatrix();

#ifdef PARSED_BY_DOXYGEN

    //! \copydoc StaticMatrix::operator()(size_t, size_t)
    scalar_reference operator()(size_t r, size_t c);

    //! \copydoc StaticMatrix::operator()(size_t, size_t) const
    scalar_const_reference operator()(size_t r, size_t c) const;

    //! \copydoc StaticMatrix::begin()
    iterator begin() noexcept;

    //! \copydoc StaticMatrix::end()
    iterator end() noexcept;

    //! \copydoc StaticMatrix::cbegin()
    const_iterator cbegin() noexcept;

    //! \copydoc StaticMatrix::cend()
    const_iterator cend() noexcept;

    //! \copydoc StaticMatrix::operator==
    bool operator==(DynamicMatrix const& that) const;

    //! \copydoc StaticMatrix::operator==
    bool operator==(RowView const& that) const;

    //! \copydoc StaticMatrix::operator!=
    bool operator!=(DynamicMatrix const& that) const;

    //! \copydoc StaticMatrix::operator!=
    bool operator!=(RowView const& that) const;

    //! \copydoc StaticMatrix::operator<
    bool operator<(DynamicMatrix const& that) const;

    //! \copydoc StaticMatrix::operator<
    bool operator<(RowView const& that) const;

    //! \copydoc StaticMatrix::operator>
    bool operator>(DynamicMatrix const& that) const;

    //! \copydoc StaticMatrix::coords
    std::pair<scalar_type, scalar_type> coords(const_iterator it) const;

    //! \copydoc StaticMatrix::number_of_rows
    size_t number_of_rows() const noexcept;

    //! \copydoc StaticMatrix::number_of_cols
    size_t number_of_cols() const noexcept;

    //! \copydoc StaticMatrix::operator+
    DynamicMatrix operator+(DynamicMatrix const& that);

    //! \copydoc StaticMatrix::operator+=
    void operator+=(DynamicMatrix const& that);

    //! \copydoc StaticMatrix::operator+=
    void operator+=(RowView const& that);

    //! \copydoc StaticMatrix::operator+=(scalar_type)
    void operator+=(scalar_type a);

    //! \copydoc StaticMatrix::operator*
    DynamicMatrix operator*(DynamicMatrix const& that);

    //! \copydoc StaticMatrix::operator*=
    void operator*=(scalar_type a);

    //! \copydoc StaticMatrix::product_inplace_no_checks
    void product_inplace_no_checks(DynamicMatrix const& x,
                                   DynamicMatrix const& y);

    //! \copydoc StaticMatrix::row_no_checks
    RowView row_no_checks(size_t i) const;

    //! \copydoc StaticMatrix::row
    RowView row(size_t i) const;

    //! \copydoc StaticMatrix::rows
    template <typename T>
    void rows(T& x) const;

    //! \copydoc StaticMatrix::transpose_no_checks
    void transpose_no_checks();

    //! \copydoc StaticMatrix::transpose
    void transpose();

    //! \copydoc StaticMatrix::hash_value
    size_t hash_value() const;

    //! \copydoc StaticMatrix::operator<=
    bool operator<=(T const& that) const;

    //! \copydoc StaticMatrix::operator>=
    bool operator>=(T const& that) const;

    //! \copydoc StaticMatrix::operator*(scalar_type)
    StaticMatrix operator*(scalar_type a);

    //! \copydoc StaticMatrix::operator+(scalar_type)
    StaticMatrix operator+(scalar_type a);

    //! \copydoc StaticMatrix::scalar_zero
    scalar_type scalar_zero() const noexcept;

    //! \copydoc StaticMatrix::scalar_one
    scalar_type scalar_one() const noexcept;

    //! \copydoc StaticMatrix::semiring
    semiring_type const* semiring() const noexcept;
#else
    using MatrixCommon::at;
    using MatrixCommon::begin;
    using MatrixCommon::cbegin;
    using MatrixCommon::cend;
    using MatrixCommon::coords;
    using MatrixCommon::end;
    using MatrixCommon::hash_value;
    using MatrixCommon::number_of_cols;
    using MatrixCommon::number_of_rows;
    using MatrixCommon::one;
    using MatrixCommon::operator();
    using MatrixCommon::operator==;
    using MatrixCommon::operator<;  // NOLINT(whitespace/operators)
    using MatrixCommon::operator<=;
    using MatrixCommon::operator>;  // NOLINT(whitespace/operators)
    using MatrixCommon::operator>=;
    using MatrixCommon::operator!=;
    using MatrixCommon::operator*=;
    using MatrixCommon::operator+=;
    using MatrixCommon::operator*;
    using MatrixCommon::operator+;
    using MatrixCommon::product_inplace_no_checks;
    using MatrixCommon::row;
    using MatrixCommon::row_no_checks;
    using MatrixCommon::rows;
    using MatrixCommon::scalar_one;
    using MatrixCommon::scalar_zero;
    using MatrixCommon::semiring;
    // using MatrixCommon::swap; // Don't want this use the one below.
    using MatrixCommon::transpose;
    using MatrixCommon::transpose_no_checks;
#endif  // PARSED_BY_DOXYGEN

    //! \copydoc StaticMatrix::swap
    void swap(DynamicMatrix& that) noexcept {
      static_cast<MatrixDynamicDim&>(*this).swap(
          static_cast<MatrixDynamicDim&>(that));
      static_cast<MatrixCommon&>(*this).swap(static_cast<MatrixCommon&>(that));
      std::swap(_semiring, that._semiring);
    }

   private:
    using MatrixCommon::resize;

    scalar_type plus_no_checks_impl(scalar_type x,
                                    scalar_type y) const noexcept {
      return _semiring->plus_no_checks(x, y);
    }

    scalar_type product_no_checks_impl(scalar_type x,
                                       scalar_type y) const noexcept {
      return _semiring->product_no_checks(x, y);
    }

    scalar_type one_impl() const noexcept {
      return _semiring->scalar_one();
    }

    scalar_type zero_impl() const noexcept {
      return _semiring->scalar_zero();
    }

    Semiring const* semiring_impl() const noexcept {
      return _semiring;
    }

    Semiring const* _semiring;
  };

  template <typename Semiring, typename Scalar>
  DynamicMatrix<Semiring, Scalar>::~DynamicMatrix() = default;

  ////////////////////////////////////////////////////////////////////////
  // Helper structs to check if matrix is static, or has a pointer to a
  // semiring
  ////////////////////////////////////////////////////////////////////////

  namespace detail {
    template <typename T>
    struct IsStaticMatrixHelper : std::false_type {};

    template <typename PlusOp,
              typename ProdOp,
              typename ZeroOp,
              typename OneOp,
              size_t R,
              size_t C,
              typename Scalar>
    struct IsStaticMatrixHelper<
        StaticMatrix<PlusOp, ProdOp, ZeroOp, OneOp, R, C, Scalar>>
        : std::true_type {};

    template <typename T>
    struct IsMatWithSemiringHelper : std::false_type {};

    template <typename Semiring, typename Scalar>
    struct IsMatWithSemiringHelper<DynamicMatrix<Semiring, Scalar>>
        : std::true_type {};

    template <typename S, typename T = void>
    struct IsTruncMatHelper : std::false_type {};
  }  // namespace detail

  //! \brief Helper variable template.
  //!
  //! \ingroup matrix_group
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This variable has value \c true if the template parameter \c T is
  //! StaticMatrix; and \c false otherwise.
  //!
  //! \tparam T the type to check.
  template <typename T>
  constexpr bool IsStaticMatrix = detail::IsStaticMatrixHelper<T>::value;

  //! \ingroup matrix_group
  //!
  //! \brief Helper variable template.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This variable has value \c true if the template parameter \c T is
  //! DynamicMatrix; and \c false otherwise.
  //!
  //! \tparam T the type to check.
  template <typename T>
  constexpr bool IsDynamicMatrix = IsMatrix<T> && !IsStaticMatrix<T>;

  //! \ingroup matrix_group
  //!
  //! \brief Helper variable template.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This variable has value \c true if the template parameter \p T is
  //! DynamicMatrix<Semiring, Scalar>; and \c false otherwise.
  //!
  //! \tparam T the type to check.
  template <typename T>
  static constexpr bool IsMatWithSemiring
      = detail::IsMatWithSemiringHelper<T>::value;

  namespace detail {

    template <typename T>
    static constexpr bool IsTruncMat = IsTruncMatHelper<T>::value;

    template <typename Mat>
    void throw_if_semiring_nullptr(Mat const& m) {
      if (IsMatWithSemiring<Mat> && m.semiring() == nullptr) {
        LIBSEMIGROUPS_EXCEPTION(
            "the matrix's pointer to a semiring is nullptr!")
      }
    }
  }  // namespace detail

  //! \ingroup matrix_group
  //!
  //! \brief Namespace for helper functions for matrices.
  //!
  //! This namespace contains various helper functions for the various matrix
  //! classes in `libsemigroups`. These functions could have been member
  //! functions of the matrix classes but they only use public member functions,
  //! and so they are declared as free functions instead.
  namespace matrix {

    //! \brief Returns the threshold of a matrix.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! The threshold of a matrix that does not have entries in a truncated
    //! semiring is \ref UNDEFINED, and this function returns this value.
    //!
    //! \tparam Mat the type of the matrix.
    //!
    //! \returns The value \ref UNDEFINED.
    //!
    //! \exceptions
    //! \noexcept
    template <typename Mat>
    constexpr auto threshold(Mat const&) noexcept
        -> std::enable_if_t<!detail::IsTruncMat<Mat>,
                            typename Mat::scalar_type> {
      return UNDEFINED;
    }

    //! \brief Returns the threshold of a matrix over a truncated semiring.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function returns the threshold of a matrix over a truncated
    //! semiring.
    //!
    //! \tparam Mat the type of the matrix.
    //!
    //! \returns The threshold of any matrix of type \p Mat.
    //!
    //! \exceptions
    //! \noexcept
    template <typename Mat>
    constexpr auto threshold(Mat const&) noexcept
        -> std::enable_if_t<detail::IsTruncMat<Mat> && !IsMatWithSemiring<Mat>,
                            typename Mat::scalar_type> {
      return detail::IsTruncMatHelper<Mat>::threshold;
    }

    //! \brief Returns the threshold of a matrix over a truncated semiring.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function returns the threshold of a matrix over a truncated
    //! semiring.
    //!
    //! \tparam Mat the type of the matrix.
    //!
    //! \param x the matrix.
    //!
    //! \returns The threshold of \p x.
    //!
    //! \exceptions
    //! \noexcept
    template <typename Mat>
    auto threshold(Mat const& x) noexcept
        -> std::enable_if_t<detail::IsTruncMat<Mat> && IsMatWithSemiring<Mat>,
                            typename Mat::scalar_type> {
      return x.semiring()->threshold();
    }
  }  // namespace matrix

  ////////////////////////////////////////////////////////////////////////
  // Boolean matrices - compile time semiring arithmetic
  ////////////////////////////////////////////////////////////////////////

  //! \defgroup bmat_group BMat
  //!
  //! This page describes the functionality for \f$m \times n\f$ boolean
  //! matrices for arbitrary dimensions \f$m\f$ and \f$n\f$. There are two
  //! types of such boolean matrices those whose dimension is known at
  //! compile-time, and those where it is not.  Both types can be accessed via
  //! the alias template `BMat<R, C>`: if \c R or \c C has value
  //! \c 0, then the dimensions can be set at run time, otherwise the
  //! dimensions are \c R and \c C. The default value of \c R is \c 0 and of
  //! \c C is \c R.
  //!
  //! The alias `BMat<R, C>` is either StaticMatrix, \ref
  //! DynamicMatrixStaticArith "DynamicMatrix (compile-time arithmetic)", or
  //! \ref DynamicMatrixDynamicArith "DynamicMatrix (run-time arithmetic)",
  //! please refer to the documentation of these class templates for more
  //! details. The only substantial difference in the interface of static
  //! and dynamic is that the former can be default constructed and the
  //! latter should be constructed using the dimensions.
  //!
  //! Adapters for \ref BMat objects are documented \ref adapters_bmat_group
  //! "here".
  //!
  //! \par Example
  //! \code
  //! // default construct an uninitialized 3 x 3 static matrix
  //! BMat<3> m;
  //! // construct an uninitialized 4 x 4 dynamic matrix
  //! BMat<> m(4, 4);
  //! \endcode

  //! \ingroup bmat_group
  //!
  //! \brief Function object for addition in the boolean semiring.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `bool operator()(bool x, bool y) const noexcept` which returns `x || y`;
  //! representing addition in the boolean semiring.
  struct BooleanPlus {
    //! \brief Call operator for addition.
    //!
    //! This function returns the sum of its arguments in the boolean semiring.
    //!
    //! \param x the first value.
    //! \param y the second value.
    //!
    //! \returns The sum of \p x and \p y in the boolean semiring.
    //!
    //! \exceptions
    //! \noexcept
    constexpr bool operator()(bool x, bool y) const noexcept {
      return x || y;
    }
  };

  //! \ingroup bmat_group
  //!
  //! \brief Function object for multiplication in the boolean semiring.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `bool operator()(bool x, bool y) const noexcept` which returns `x && y`;
  //! representing multiplication in the boolean semiring.
  struct BooleanProd {
    //! \brief Call operator for multiplication.
    //!
    //! This function returns the product of its arguments in the boolean
    //! semiring.
    //!
    //! \param x the first value.
    //! \param y the second value.
    //!
    //! \returns The product of \p x and \p y in the boolean semiring.
    //!
    //! \exceptions
    //! \noexcept
    constexpr bool operator()(bool x, bool y) const noexcept {
      return x & y;
    }
  };

  //! \ingroup bmat_group
  //!
  //! \brief Function object for returning the multiplicative identity.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `bool operator()() const noexcept` which returns \c true; representing
  //! the multiplicative identity of the boolean semiring.
  struct BooleanOne {
    //! \brief Call operator returning the multiplication identity \c true of
    //! the boolean semiring.
    //!
    //! Call operator returning the multiplication identity \c true of
    //! the boolean semiring.
    //!
    //! \returns The value \c true.
    //!
    //! \exceptions
    //! \noexcept
    constexpr bool operator()() const noexcept {
      return true;
    }
  };

  //! \ingroup bmat_group
  //!
  //! \brief Function object for returning the additive identity.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `bool operator()() const noexcept` which returns \c false;
  //! representing the additive identity of the boolean semiring.
  struct BooleanZero {
    //! \brief Call operator returning the additive identity \c false of
    //! the boolean semiring.
    //!
    //! Call operator returning the additive identity \c false of
    //! the boolean semiring.
    //!
    //! \returns The value \c false.
    //!
    //! \exceptions
    //! \noexcept
    constexpr bool operator()() const noexcept {
      return false;
    }
  };

  //! \ingroup bmat_group
  //!
  //! \brief Alias for dynamic boolean matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for the type of dynamic boolean matrices where the dimensions of the
  //! matrices can be defined at run time.
  // The use of `int` rather than `bool` as the scalar type for dynamic
  // boolean matrices is intentional, because the bit iterators implemented in
  // std::vector<bool> entail a significant performance penalty.
  using DynamicBMat
      = DynamicMatrix<BooleanPlus, BooleanProd, BooleanZero, BooleanOne, int>;

  //! \ingroup bmat_group
  //!
  //! \brief Alias for static boolean matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for the type of static boolean matrices where the dimensions of the
  //! matrices are defined at compile time.
  //!
  //! \tparam R the number of rows.
  //! \tparam C the number of columns.
  template <size_t R, size_t C>
  using StaticBMat = StaticMatrix<BooleanPlus,
                                  BooleanProd,
                                  BooleanZero,
                                  BooleanOne,
                                  R,
                                  C,
                                  int>;

  //! \ingroup bmat_group
  //!
  //! \brief Alias template for boolean matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias template for boolean matrices.
  //!
  //! \tparam R the number of rows of the matrix. A value of \c 0 (the default
  //! value) indicates that the dimensions will be set at run time.
  //!
  //! \tparam C the number of columns of the matrix. A value of \c 0 indicates
  //! that the dimensions will be set at run time (the default value is `\c R`).
  // FLS + JDM considered adding BMat8 and decided it wasn't a good idea.
  template <size_t R = 0, size_t C = R>
  using BMat
      = std::conditional_t<R == 0 || C == 0, DynamicBMat, StaticBMat<R, C>>;

  namespace detail {
    template <typename T>
    struct IsBMatHelper : std::false_type {};

    template <size_t R, size_t C>
    struct IsBMatHelper<StaticBMat<R, C>> : std::true_type {};

    template <>
    struct IsBMatHelper<DynamicBMat> : std::true_type {};

    template <typename T>
    struct BitSetCapacity {
      static constexpr size_t value = BitSet<1>::max_size();
    };

    template <size_t R, size_t C>
    struct BitSetCapacity<StaticBMat<R, C>> {
      static_assert(R == C, "the number of rows and columns must be equal");
      static constexpr size_t value = R;
    };
  }  // namespace detail

  //! \ingroup bmat_group
  //!
  //! \brief Helper to check if a type is \ref BMat.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This variable has value \c true if the template parameter \c T is the
  //! same as `BMat<R, C>` for some values of \c R and \c C.
  //!
  //! \tparam T the type to check.
  template <typename T>
  static constexpr bool IsBMat = detail::IsBMatHelper<T>::value;

  namespace detail {
    // This function is required for exceptions and to_human_readable_repr, so
    // that if we encounter an entry of a matrix (Scalar type), then it can be
    // printed correctly. If we just did fmt::format("{}", val) and val ==
    // POSITIVE_INFINITY, but the type of val is, say, size_t, then this
    // wouldn't use the formatter for PositiveInfinity.
    template <typename Scalar>
    std::string entry_repr(Scalar a) {
      static_assert(std::is_integral_v<Scalar>);

      if constexpr (std::is_signed_v<Scalar>) {
        if (a == NEGATIVE_INFINITY) {
          return u8"-\u221E";
        }
      }
      if (a == static_cast<Scalar>(POSITIVE_INFINITY)) {
        return u8"+\u221E";
      }
      return fmt::format("{}", a);
    }
  }  // namespace detail

  namespace matrix {

    //! \ingroup bmat_group
    //!
    //! \brief Check the entries in a boolean matrix are valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function can be used to validate that a matrix contains values in
    //! the underlying semiring.
    //!
    //! \tparam Mat the type of \p m, must satisfy \ref IsBMat<Mat>.
    //!
    //! \param m the boolean matrix to validate.
    //!
    //! \throws LibsemigroupsException if any of the entries in the matrix are
    //! not \c 0 or \c 1. The values in a boolean matrix are of type \c int, but
    //! a matrix shouldn't contain values except \c 0 and \c 1.
    template <typename Mat>
    std::enable_if_t<IsBMat<Mat>> throw_if_bad_entry(Mat const& m) {
      using scalar_type = typename Mat::scalar_type;
      auto it           = std::find_if_not(
          m.cbegin(), m.cend(), [](scalar_type x) { return x == 0 || x == 1; });
      if (it != m.cend()) {
        auto [r, c] = m.coords(it);
        LIBSEMIGROUPS_EXCEPTION(
            "invalid entry, expected 0 or 1 but found {} in entry ({}, {})",
            detail::entry_repr(*it),
            r,
            c);
      }
    }

    //! \ingroup bmat_group
    //!
    //! \brief Check an entry in a boolean matrix is valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function can be used to validate that a matrix contains values in
    //! the underlying semiring.
    //!
    //! \tparam Mat the type of the 1st argument, must satisfy \ref IsBMat<Mat>.
    //!
    //! \param val the entry to validate.
    //!
    //! \throws LibsemigroupsException if \p val is not \c 0 or \c 1. The
    //! values in a boolean matrix are of type \c int, but a matrix shouldn't
    //! contain values except \c 0 and \c 1.
    template <typename Mat>
    std::enable_if_t<IsBMat<Mat>>
    throw_if_bad_entry(Mat const&, typename Mat::scalar_type val) {
      if (val != 0 && val != 1) {
        LIBSEMIGROUPS_EXCEPTION("invalid entry, expected 0 or 1 but found {}",
                                detail::entry_repr(val));
      }
    }
  }  // namespace matrix

  ////////////////////////////////////////////////////////////////////////
  // Integer matrices - compile time semiring arithmetic
  ////////////////////////////////////////////////////////////////////////

  //! \defgroup intmat_group Integer matrices
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This page describes the functionality for \f$m \times n\f$  matrices of
  //! integers for arbitrary dimensions \f$m\f$ and \f$n\f$. There are two
  //! types of such matrices: those whose dimension is known at compile-time,
  //! and those where it is not.  Both types can be accessed via the alias
  //! template IntMat: if \c R or \c C has value \c 0, then the dimensions
  //! can be set at run time, otherwise \c R and \c C are the dimension. The
  //! default value of \c R is \c 0, and of \c C is \c R.
  //!
  //! The alias \ref IntMat is either StaticMatrix, \ref
  //! DynamicMatrixStaticArith "DynamicMatrix (compile-time arithmetic)", or
  //! \ref DynamicMatrixDynamicArith "DynamicMatrix (run-time arithmetic)",
  //! please refer to the documentation of these class templates for more
  //! details. The only substantial difference in the interface of static
  //! and dynamic is that the former can be default constructed and the
  //! latter should be constructed using the dimensions.
  //!
  //! \par Example
  //! \code
  //! // default construct an uninitialized 3 x 3 static matrix
  //! IntMat<3> m;
  //! // construct an uninitialized 4 x 4 dynamic matrix
  //! IntMat<>  m(4, 4);
  //! \endcode

  //! \ingroup intmat_group
  //!
  //! \brief Function object for addition in the ring of integers.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `Scalar operator()(Scalar x, Scalar y) const noexcept`
  //! which returns the usual sum `x + y` of \c x and \c y; representing
  //! addition in the integer semiring.
  //!
  //! \tparam Scalar the type of the entries in the matrix.
  template <typename Scalar>
  struct IntegerPlus {
    //! \brief Call operator for addition.
    //!
    //! This function returns the sum of its arguments in the ring of integers.
    //!
    //! \param x the first value.
    //! \param y the second value.
    //!
    //! \returns The sum of \p x and \p y in the boolean semiring.
    //!
    //! \exceptions
    //! \noexcept
    constexpr Scalar operator()(Scalar x, Scalar y) const noexcept {
      return x + y;
    }
  };

  //! \ingroup intmat_group
  //!
  //! \brief Function object for multiplication in the ring of integers.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `Scalar operator()(Scalar x, Scalar y) const noexcept`
  //! which returns the usual product `x * y` of \c x and \c y;
  //! representing multiplication in the integer semiring.
  //!
  //! \tparam Scalar the type of the entries in the matrix.
  template <typename Scalar>
  struct IntegerProd {
    //! \brief Call operator for multiplication.
    //!
    //! This function returns the product of its arguments in the ring of
    //! integers.
    //!
    //! \param x the first value.
    //! \param y the second value.
    //!
    //! \returns The product of \p x and \p y in the ring of integers.
    //!
    //! \exceptions
    //! \noexcept
    constexpr Scalar operator()(Scalar x, Scalar y) const noexcept {
      return x * y;
    }
  };

  //! \ingroup intmat_group
  //!
  //! \brief Function object for returning the additive identity.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `Scalar operator()() const noexcept` which returns \c 0; representing
  //! the additive identity of the integer semiring.
  template <typename Scalar>
  struct IntegerZero {
    //! \brief Call operator returning the integer \c 0.
    //!
    //! This function returns the additive identity of the ring of
    //! integers (i.e. the integer \c 0).
    //!
    //! \returns The integer \c 0.
    //!
    //! \exceptions
    //! \noexcept
    constexpr Scalar operator()() const noexcept {
      return 0;
    }
  };

  //! \ingroup intmat_group
  //!
  //! \brief Function object for returning the multiplicative identity.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `Scalar operator()() const noexcept` which returns \c 1; representing
  //! the multiplicative identity of the integer semiring.
  template <typename Scalar>
  struct IntegerOne {
    //! \brief Call operator returning the integer \c 1.
    //!
    //! This function returns the multiplicative identity of the ring of
    //! integers (i.e. the integer \c 1).
    //!
    //! \returns The integer \c 1.
    //!
    //! \exceptions
    //! \noexcept
    constexpr Scalar operator()() const noexcept {
      return 1;
    }
  };

  //! \ingroup intmat_group
  //!
  //! \brief Alias for dynamic integer matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for the type of dynamic integer matrices where the dimensions of the
  //! matrices can be defined at run time.
  //!
  //! \tparam Scalar the type of the entries in the matrix (default: `int`).
  template <typename Scalar>
  using DynamicIntMat = DynamicMatrix<IntegerPlus<Scalar>,
                                      IntegerProd<Scalar>,
                                      IntegerZero<Scalar>,
                                      IntegerOne<Scalar>,
                                      Scalar>;

  //! \ingroup intmat_group
  //!
  //! \brief Alias for static integer matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for the type of static integer matrices where the dimensions of the
  //! matrices can be defined at compile time.
  //!
  //! \tparam R the number of rows of the matrix. A value of \c 0 (the default
  //! value) indicates that the dimensions will be set at run time.
  //!
  //! \tparam C the number of columns of the matrix. A value of \c 0 indicates
  //! that the dimensions will be set at run time (the default value is \c R).
  //!
  //! \tparam Scalar the type of the entries in the matrix (default: `int`).
  template <size_t R, size_t C, typename Scalar>
  using StaticIntMat = StaticMatrix<IntegerPlus<Scalar>,
                                    IntegerProd<Scalar>,
                                    IntegerZero<Scalar>,
                                    IntegerOne<Scalar>,
                                    R,
                                    C,
                                    Scalar>;

  //! \ingroup intmat_group
  //!
  //! \brief Alias template for integer matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias template for integer matrices.
  //!
  //! \tparam R the number of rows of the matrix. A value of \c 0 (the
  //! default value) indicates that the dimensions will be set at run time.
  //!
  //! \tparam C the number of columns of the matrix. A value of \c 0
  //! indicates that the dimensions will be set at run time (the default value
  //! is \c R).
  //!
  //! \tparam Scalar the type of the entries in the matrix (default: `int`).
  template <size_t R = 0, size_t C = R, typename Scalar = int>
  using IntMat = std::conditional_t<R == 0 || C == 0,
                                    DynamicIntMat<Scalar>,
                                    StaticIntMat<R, C, Scalar>>;
  namespace detail {
    template <typename T>
    struct IsIntMatHelper : std::false_type {};

    template <size_t R, size_t C, typename Scalar>
    struct IsIntMatHelper<StaticIntMat<R, C, Scalar>> : std::true_type {};

    template <typename Scalar>
    struct IsIntMatHelper<DynamicIntMat<Scalar>> : std::true_type {};
  }  // namespace detail

  //! \ingroup matrix_group
  //!
  //! \brief Helper variable template.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This variable has value \c true if the template parameter \c T is
  //! \ref IntMat; and \c false otherwise.
  //!
  //! \tparam T the type to check.
  template <typename T>
  static constexpr bool IsIntMat = detail::IsIntMatHelper<T>::value;

  namespace matrix {
    //! \ingroup intmat_group
    //!
    //! \brief Validate that an integer matrix is valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function throws an exception if the entries of an integer matrix
    //! are not valid, which is if and only if any of the entries equal \ref
    //! POSITIVE_INFINITY or \ref NEGATIVE_INFINITY.
    //!
    //! \tparam Mat the type of the argument \p x, must satisfy \ref
    //! IsMatrix<Mat>.
    //!
    //! \param x the matrix to check.
    template <typename Mat>
    constexpr std::enable_if_t<IsIntMat<Mat>> throw_if_bad_entry(Mat const& x) {
      using scalar_type = typename Mat::scalar_type;
      auto it = std::find_if(x.cbegin(), x.cend(), [](scalar_type val) {
        return val == POSITIVE_INFINITY || val == NEGATIVE_INFINITY;
      });
      if (it != x.cend()) {
        auto [r, c] = x.coords(it);
        LIBSEMIGROUPS_EXCEPTION(
            "invalid entry, expected entries to be integers, "
            "but found {} in entry ({}, {})",
            detail::entry_repr(*it),
            r,
            c);
      }
    }

    //! \ingroup intmat_group
    //!
    //! \brief Validate that an entry in an integer matrix is valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function throws an exception if the entry \p val of an integer
    //! matrix are not valid, which is if and only if the entry \p val equals
    //! \ref POSITIVE_INFINITY or \ref NEGATIVE_INFINITY.
    //!
    //! The 1st argument is used for overload resolution.
    //!
    //! \tparam Mat the type of the 1st argument, must satisfy \ref
    //! IsMatrix<Mat>.
    //!
    //! \param val the entry to check.
    template <typename Mat>
    std::enable_if_t<IsIntMat<Mat>>
    throw_if_bad_entry(Mat const&, typename Mat::scalar_type val) {
      if (val == POSITIVE_INFINITY || val == NEGATIVE_INFINITY) {
        LIBSEMIGROUPS_EXCEPTION(
            "invalid entry, expected entries to be integers, "
            "but found {}",
            detail::entry_repr(val));
      }
    }
  }  // namespace matrix

  ////////////////////////////////////////////////////////////////////////
  // Max-plus matrices
  ////////////////////////////////////////////////////////////////////////
  //! \defgroup maxplusmat_group Max-plus matrices
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This page describes the functionality for \f$n \times n\f$  matrices over
  //! the max-plus semiring for arbitrary dimension \f$n\f$. There are two types
  //! of such matrices those whose dimension is known at compile-time, and those
  //! where it is not.  Both types can be accessed via the alias template
  //! \ref MaxPlusMat if \c N has value \c 0, then the dimensions can be set at
  //! run time, otherwise \c N is the dimension. The default value of \c N is \c
  //! 0.
  //!
  //! The alias \ref MaxPlusMat is either StaticMatrix, \ref
  //! DynamicMatrixStaticArith "DynamicMatrix (compile-time arithmetic)", or
  //! \ref DynamicMatrixDynamicArith "DynamicMatrix (run-time arithmetic)",
  //! please refer to the documentation of these class templates for more
  //! details. The only substantial difference in the interface of static
  //! and dynamic matrices is that the former can be default constructed and the
  //! latter should be constructed using the dimensions.
  //!
  //! \par Example
  //! \code
  //! // default construct an uninitialized 3 x 3 static matrix
  //! MaxPlusMat<3> m;
  //! // construct an uninitialized 4 x 4 dynamic matrix
  //! MaxPlusMat<>  m(4, 4);
  //! \endcode

  //! \ingroup maxplusmat_group
  //!
  //! \brief Function object for addition in the max-plus semiring.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! \tparam Scalar the type of the values in the semiring (must be signed
  //! integer type).
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `Scalar operator()(Scalar x, Scalar y) const noexcept`
  //! that returns \f$x \oplus y\f$ which is defined by
  //!
  //! \f[
  //!    x\oplus y =
  //!    \begin{cases}
  //!    \max\{x, y\}   & \text{if } x \neq -\infty\text{ and }y \neq -\infty
  //!    \\ \mbox{}
  //!    -\infty & \text{if } x = -\infty \text{ or }y = -\infty;
  //!    \end{cases}
  //! \f]
  //!
  //! representing addition in the max-plus semiring.
  // Static arithmetic
  template <typename Scalar>
  struct MaxPlusPlus {
    static_assert(std::is_signed<Scalar>::value,
                  "MaxPlus requires a signed integer type as parameter!");
    //! \brief Call operator for addition.
    //!
    //! This function returns the sum of its arguments in the max-plus semiring.
    //!
    //! \param x the first value.
    //! \param y the second value.
    //!
    //! \returns The sum of \p x and \p y in the max-plus semiring.
    //!
    //! \exceptions
    //! \noexcept
    Scalar operator()(Scalar x, Scalar y) const noexcept {
      if (x == NEGATIVE_INFINITY) {
        return y;
      } else if (y == NEGATIVE_INFINITY) {
        return x;
      }
      return std::max(x, y);
    }
  };

  //! \ingroup maxplusmat_group
  //!
  //! \brief Function object for multiplication in the max-plus semiring.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //!  This is a stateless struct with a single call operator of signature:
  //!  `Scalar operator()(Scalar x, Scalar y) const noexcept`
  //!  that returns \f$x \otimes y\f$ which is defined bitset_type
  //!  \f[
  //!     x\otimes y =
  //!     \begin{cases}
  //!     x + y   & \text{if } x \neq -\infty\text{ and }y \neq -\infty
  //!     \\ \mbox{}
  //!     -\infty & \text{if } x = -\infty \text{ or }y = -\infty;
  //!     \end{cases}
  //!  \f]
  //!  representing multiplication in the max-plus semiring.
  //!
  //! \tparam Scalar the type of the values in the semiring (must be signed
  //! integer type).
  template <typename Scalar>
  struct MaxPlusProd {
    static_assert(std::is_signed<Scalar>::value,
                  "MaxPlus requires a signed integer type as parameter!");
    //! \brief Call operator for multiplication.
    //!
    //! This function returns the product of its arguments in the max-plus
    //! semiring.
    //!
    //! \param x the first value.
    //! \param y the second value.
    //!
    //! \returns The product of \p x and \p y in the max-plus semiring.
    //!
    //! \exceptions
    //! \noexcept
    Scalar operator()(Scalar x, Scalar y) const noexcept {
      if (x == NEGATIVE_INFINITY || y == NEGATIVE_INFINITY) {
        return NEGATIVE_INFINITY;
      }
      return x + y;
    }
  };

  //! \ingroup maxplusmat_group
  //!
  //! \brief Function object for returning the additive identity of the max-plus
  //! semiring.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `Scalar operator()() const noexcept` which returns \f$-\infty\f$;
  //! representing the additive identity of the max-plus semiring.
  //!
  //! \tparam Scalar the type of the values in the semiring (must be signed
  //! integer type).
  template <typename Scalar>
  struct MaxPlusZero {
    static_assert(std::is_signed<Scalar>::value,
                  "MaxPlus requires a signed integer type as parameter!");
    //! \brief Call operator for additive identity.
    //!
    //! This function returns the additive identity in the max-plus
    //! semiring.
    //!
    //! \returns The additive identity in the max-plus semiring.
    //!
    //! \exceptions
    //! \noexcept
    constexpr Scalar operator()() const noexcept {
      return NEGATIVE_INFINITY;
    }
  };

  //! \ingroup maxplusmat_group
  //!
  //! \brief Alias for dynamic max-plus matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for the type of dynamic max-plus matrices where the dimensions of
  //! the matrices can be defined at run time.
  //!
  //! \tparam Scalar the type of the entries in the matrix.
  template <typename Scalar>
  using DynamicMaxPlusMat = DynamicMatrix<MaxPlusPlus<Scalar>,
                                          MaxPlusProd<Scalar>,
                                          MaxPlusZero<Scalar>,
                                          IntegerZero<Scalar>,
                                          Scalar>;

  //! \ingroup maxplusmat_group
  //!
  //! \brief Alias for static max-plus matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for static max-plus matrices whose arithmetic and dimensions are
  //! defined at compile-time.
  //!
  //! \tparam R the number of rows.
  //! \tparam C the number of columns.
  //! \tparam Scalar the type of the entries in the matrix.
  template <size_t R, size_t C, typename Scalar>
  using StaticMaxPlusMat = StaticMatrix<MaxPlusPlus<Scalar>,
                                        MaxPlusProd<Scalar>,
                                        MaxPlusZero<Scalar>,
                                        IntegerZero<Scalar>,
                                        R,
                                        C,
                                        Scalar>;

  //! \ingroup maxplusmat_group
  //!
  //! \brief Alias template for max-plus matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias template for max-plus matrices.
  //!
  //! \tparam R the number of rows.  A value of \c 0 indicates that the value
  //! will be set at run time (default: \c 0).
  //!
  //! \tparam C the number of columns.  A value of \c 0 indicates that the value
  //! will be set at run time (default: \c R).
  //!
  //! \tparam Scalar the type of the entries in the matrix (default: `int`).
  template <size_t R = 0, size_t C = R, typename Scalar = int>
  using MaxPlusMat = std::conditional_t<R == 0 || C == 0,
                                        DynamicMaxPlusMat<Scalar>,
                                        StaticMaxPlusMat<R, C, Scalar>>;

  namespace detail {
    template <typename T>
    struct IsMaxPlusMatHelper : std::false_type {};

    template <size_t R, size_t C, typename Scalar>
    struct IsMaxPlusMatHelper<StaticMaxPlusMat<R, C, Scalar>> : std::true_type {
    };

    template <typename Scalar>
    struct IsMaxPlusMatHelper<DynamicMaxPlusMat<Scalar>> : std::true_type {};
  }  // namespace detail

  //! \ingroup matrix_group
  //!
  //! \brief Helper variable template.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This variable has value \c true if the template parameter \c T is
  //! \ref MaxPlusMat; and \c false otherwise.
  //!
  //! \tparam T the type to check.
  template <typename T>
  static constexpr bool IsMaxPlusMat = detail::IsMaxPlusMatHelper<T>::value;

  namespace matrix {
    //! \ingroup maxplusmat_group
    //!
    //! \brief Validate that a max-plus matrix is valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function can be used to validate that a matrix contains values in
    //! the underlying semiring.
    //!
    //! \tparam Mat the type of the matrix, must satisfy \ref IsMaxPlusMat.
    //!
    //! \param x the matrix to check.
    //!
    //! \throws LibsemigroupsException if \p x contains \ref POSITIVE_INFINITY.
    template <typename Mat>
    auto throw_if_bad_entry(Mat const& x)
        -> std::enable_if_t<IsMaxPlusMat<Mat>> {
      using scalar_type = typename Mat::scalar_type;
      auto it = std::find_if(x.cbegin(), x.cend(), [](scalar_type val) {
        return val == POSITIVE_INFINITY;
      });
      if (it != x.cend()) {
        auto [r, c] = x.coords(it);
        LIBSEMIGROUPS_EXCEPTION(
            "invalid entry, expected entries to be integers or {}, "
            "but found {} in entry ({}, {})",
            NEGATIVE_INFINITY,
            POSITIVE_INFINITY,
            r,
            c);
      }
    }

    //! \ingroup maxplusmat_group
    //!
    //! \brief Validate that an entry in a max-plus matrix is valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function can be used to validate that an entry in a matrix belongs
    //! to the underlying semiring.
    //!
    //! \tparam Mat the type of the matrix, must satisfy \ref IsMaxPlusMat<Mat>.
    //!
    //! \param val the entry to check.
    //!
    //! \throws LibsemigroupsException if \p val is \ref POSITIVE_INFINITY.
    template <typename Mat>
    std::enable_if_t<IsMaxPlusMat<Mat>>
    throw_if_bad_entry(Mat const&, typename Mat::scalar_type val) {
      if (val == POSITIVE_INFINITY) {
        LIBSEMIGROUPS_EXCEPTION("invalid entry, expected entries to be "
                                "integers or {} but found {}",
                                NEGATIVE_INFINITY,
                                POSITIVE_INFINITY);
      }
    }
  }  // namespace matrix

  ////////////////////////////////////////////////////////////////////////
  // Min-plus matrices
  ////////////////////////////////////////////////////////////////////////

  //! \defgroup minplusmat_group Min-plus matrices
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This page describes the functionality for \f$n \times n\f$  matrices
  //! over the min-plus semiring for arbitrary dimension \f$n\f$. There are
  //! two types of such matrices those whose dimension is known at
  //! compile-time, and those where it is not.  Both types can be accessed via
  //! the alias template \ref MinPlusMatAlias "MinPlusMat": if \c N has value
  //! \c 0, then the dimensions can be set at run time, otherwise \c N is the
  //! dimension. The default value of \c N is \c 0.
  //!
  //! The alias \ref MinPlusMatAlias "MinPlusMat" is either \ref StaticMatrix,
  //! \ref DynamicMatrixStaticArith "DynamicMatrix (compile-time arithmetic),
  //! or \ref DynamicMatrixDynamicArith "DynamicMatrix (run-time arithmetic)"
  //! please refer to the documentation of these class templates for more
  //! details. The only substantial difference in the interface of static and
  //! dynamic matrices is that the former can be default constructed and the
  //! latter should be constructed using the dimensions.
  //!
  //! \par Example
  //!
  //! \code
  //! // default construct an uninitialized 3 x 3 static matrix
  //! MinPlusMat<3> m;
  //! // construct an uninitialized 4 x 4 dynamic matrix
  //! MinPlusMat<> m(4, 4);
  //! \endcode

  //! \ingroup minplusmat_group
  //!
  //! \brief Function object for addition in the min-plus semiring.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `Scalar operator()(Scalar x, Scalar y) const noexcept` that returns \f$x
  //! \oplus y\f$ which is defined by
  //! \f[
  //!    x\oplus y =
  //!    \begin{cases}
  //!    \min\{x, y\}   & \text{if } x \neq \infty\text{ and }y \neq \infty
  //!    \\ \mbox{}
  //!    \infty & \text{if } x = \infty \text{ or }y = \infty;
  //!    \end{cases}
  //! \f]
  //! representing addition in the min-plus semiring.
  //!
  //! \tparam Scalar the type of the values in the semiring (must be signed
  //! integer type).
  //!
  // Static arithmetic
  template <typename Scalar>
  struct MinPlusPlus {
    static_assert(std::is_signed<Scalar>::value,
                  "MinPlus requires a signed integer type as parameter!");
    //! \brief Call operator for addition.
    //!
    //! This function returns the sum of its arguments in the min-plus
    //! semiring.
    //!
    //! \param x the first value.
    //! \param y the second value.
    //!
    //! \returns The sum of \p x and \p y in the min-plus semiring.
    //!
    //! \exceptions
    //! \noexcept
    Scalar operator()(Scalar x, Scalar y) const noexcept {
      if (x == POSITIVE_INFINITY) {
        return y;
      } else if (y == POSITIVE_INFINITY) {
        return x;
      }
      return std::min(x, y);
    }
  };

  //! \ingroup minplusmat_group
  //!
  //! \brief Function object for multiplication in the min-plus semiring.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `Scalar operator()(Scalar x, Scalar y) const noexcept`
  //! that returns \f$x \otimes y\f$ which is defined by
  //! \f[
  //!     x\otimes y =
  //!     \begin{cases}
  //!     x + y  & \text{if } x \neq \infty\text{ and }y \neq \infty \\ \mbox{}
  //!     \infty & \text{if } x = \infty \text{ or }y = \infty;
  //!     \end{cases}
  //! \f]
  //!  representing multiplication in the min-plus semiring.
  //!
  //! \tparam Scalar the type of the values in the semiring (must be signed
  //! integer type).
  template <typename Scalar>
  struct MinPlusProd {
    static_assert(std::is_signed<Scalar>::value,
                  "MinPlus requires a signed integer type as parameter!");
    //! \brief Call operator for multiplication.
    //!
    //! This function returns the product of its arguments in the min-plus
    //! semiring.
    //!
    //! \param x the first value.
    //! \param y the second value.
    //!
    //! \returns The product of \p x and \p y in the min-plus semiring.
    //!
    //! \exceptions
    //! \noexcept
    Scalar operator()(Scalar x, Scalar y) const noexcept {
      if (x == POSITIVE_INFINITY || y == POSITIVE_INFINITY) {
        return POSITIVE_INFINITY;
      }
      return x + y;
    }
  };

  //! \ingroup minplusmat_group
  //!
  //! \brief Function object for returning the additive identity of the
  //! min-plus semiring.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `Scalar operator()() const noexcept` which returns \f$\infty\f$;
  //! representing the additive identity of the min-plus semiring.
  //!
  //! \tparam Scalar the type of the values in the semiring (must be signed
  //! integer type).
  template <typename Scalar>
  struct MinPlusZero {
    static_assert(std::is_signed<Scalar>::value,
                  "MinPlus requires a signed integer type as parameter!");
    //! \brief Call operator for additive identity.
    //!
    //! This function returns the additive identity in the min-plus
    //! semiring.
    //!
    //! \returns The additive identity in the min-plus semiring.
    //!
    //! \exceptions
    //! \noexcept
    constexpr Scalar operator()() const noexcept {
      return POSITIVE_INFINITY;
    }
  };

  //! \ingroup minplusmat_group
  //!
  //! \brief Alias for dynamic min-plus matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for the type of dynamic min-plus matrices where the dimensions of
  //! the matrices can be defined at run time.
  //!
  //! \tparam Scalar the type of the entries in the matrix.
  template <typename Scalar>
  using DynamicMinPlusMat = DynamicMatrix<MinPlusPlus<Scalar>,
                                          MinPlusProd<Scalar>,
                                          MinPlusZero<Scalar>,
                                          IntegerZero<Scalar>,
                                          Scalar>;

  //! \ingroup minplusmat_group
  //!
  //! \brief Alias for static min-plus matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for static min-plus matrices whose arithmetic and dimensions are
  //! defined at compile-time.
  //!
  //! \tparam R  the number of rows.
  //! \tparam C  the number of columns.
  //! \tparam Scalar the type of the entries in the matrix.
  template <size_t R, size_t C, typename Scalar>
  using StaticMinPlusMat = StaticMatrix<MinPlusPlus<Scalar>,
                                        MinPlusProd<Scalar>,
                                        MinPlusZero<Scalar>,
                                        IntegerZero<Scalar>,
                                        R,
                                        C,
                                        Scalar>;
  //! \anchor MinPlusMatAlias
  //! \ingroup minplusmat_group
  //!
  //! \brief Alias template for min-plus matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias template for min-plus matrices.
  //!
  //! \tparam R the number of rows.  A value of \c 0 indicates that the value
  //! will be set at run time (default: \c 0).
  //!
  //! \tparam C the number of columns.  A value of \c 0 indicates that the
  //! value will be set at run time (default: \c R).
  //!
  //! \tparam Scalar the type of the entries in the matrix (default: `int`).
  template <size_t R = 0, size_t C = R, typename Scalar = int>
  using MinPlusMat = std::conditional_t<R == 0 || C == 0,
                                        DynamicMinPlusMat<Scalar>,
                                        StaticMinPlusMat<R, C, Scalar>>;

  namespace detail {
    template <typename T>
    struct IsMinPlusMatHelper : std::false_type {};

    template <size_t R, size_t C, typename Scalar>
    struct IsMinPlusMatHelper<StaticMinPlusMat<R, C, Scalar>> : std::true_type {
    };

    template <typename Scalar>
    struct IsMinPlusMatHelper<DynamicMinPlusMat<Scalar>> : std::true_type {};
  }  // namespace detail

  //! \ingroup matrix_group
  //!
  //! \brief Helper variable template.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This variable has value \c true if the template parameter \c T is
  //! \ref MinPlusMat; and \c false otherwise.
  //!
  //! \tparam T the type to check.
  template <typename T>
  static constexpr bool IsMinPlusMat = detail::IsMinPlusMatHelper<T>::value;

  namespace matrix {
    //! \ingroup minplusmat_group
    //!
    //! \brief Validate that a min-plus matrix is valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function can be used to validate that a matrix contains values in
    //! the underlying semiring.
    //!
    //! \tparam Mat the type of the matrix, must satisfy \ref IsMinPlusMat.
    //!
    //! \param x the matrix to check.
    //!
    //! \throws LibsemigroupsException if \p x contains \ref NEGATIVE_INFINITY.
    template <typename Mat>
    std::enable_if_t<IsMinPlusMat<Mat>> throw_if_bad_entry(Mat const& x) {
      using scalar_type = typename Mat::scalar_type;
      auto it = std::find_if(x.cbegin(), x.cend(), [](scalar_type val) {
        return val == NEGATIVE_INFINITY;
      });
      if (it != x.cend()) {
        auto [r, c] = x.coords(it);
        LIBSEMIGROUPS_EXCEPTION(
            "invalid entry, expected entries to be integers or {}, "
            "but found {} in entry ({}, {})",
            POSITIVE_INFINITY,
            NEGATIVE_INFINITY,
            r,
            c);
      }
    }

    //! \ingroup minplusmat_group
    //!
    //! \brief Validate that an entry in a min-plus matrix is valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function can be used to validate that an entry in a matrix belongs
    //! to the underlying semiring.
    //!
    //! \tparam Mat the type of the matrix, must satisfy \ref IsMinPlusMat<Mat>.
    //!
    //! \param val the entry to check.
    //!
    //! \throws LibsemigroupsException if \p val is \ref NEGATIVE_INFINITY.
    template <typename Mat>
    std::enable_if_t<IsMinPlusMat<Mat>>
    throw_if_bad_entry(Mat const&, typename Mat::scalar_type val) {
      if (val == NEGATIVE_INFINITY) {
        LIBSEMIGROUPS_EXCEPTION("invalid entry, expected entries to be "
                                "integers or {} but found {}",
                                POSITIVE_INFINITY,
                                NEGATIVE_INFINITY);
      }
    }
  }  // namespace matrix

  ////////////////////////////////////////////////////////////////////////
  // Max-plus matrices with threshold
  ////////////////////////////////////////////////////////////////////////

  //! \defgroup maxplustruncmat_group Truncated max-plus matrices
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This page describes the functionality for \f$n \times n\f$  matrices
  //! over the finite quotient of the max-plus semiring by the congruence
  //! \f$t = t + 1\f$ for arbitrary \f$n\f$ and \f$t\f$. The value
  //! \f$t\f$ is referred to as the *threshold*.
  //!
  //! There are three types of such matrices where:
  //!
  //! 1. the dimension is known at compile-time;
  //! 2. the dimension is to be defined a run time but the arithmetic
  //! operations
  //!    are known at compile-time (i.e. the value of \f$t\f$ is known at
  //!    compile time)
  //! 3. both the dimension and the arithmetic operations (i.e. \f$t\f$) are
  //!    to be defined a run time.
  //!
  //! All three of these types can be accessed via the alias template
  //! \ref MaxPlusTruncMat if \c T has value \c 0,
  //! then the threshold can be set at run time, and if \c R or \c C is \c 0,
  //! then the dimension can be set at run time.  The default value of \c T is
  //! \c 0, \c R is \c 0, and of \c C is \c R.
  //!
  //! The alias \ref MaxPlusTruncMat is either
  //! StaticMatrix, \ref DynamicMatrixStaticArith "DynamicMatrix (compile-time
  //! arithmetic)", or \ref DynamicMatrixDynamicArith "DynamicMatrix (run-time
  //! arithmetic)", please refer to the documentation of these class templates
  //! for more details. The only substantial difference in the interface of
  //! static and dynamic matrices is that the former can be default
  //! constructed and the latter should be constructed using the dimensions.
  //!
  //! \par Example
  //! \code
  //! // construct an uninitialized 3 x 3 static matrix with threshold 11
  //! MaxPlusTruncMat<11, 3> m;
  //! // construct an uninitialized 4 x 4 dynamic matrix with threshold 11
  //! MaxPlusTruncMat<11> m(4, 4);
  //! // construct a truncated max-plus semiring with threshold 11
  //! MaxPlusTruncSemiring sr(11);
  //! // construct an uninitialized 5 x 5 dynamic matrix with threshold 11
  //! // (defined at run time)
  //! MaxPlusTruncMat<>  m(sr, 5, 5);
  //! \endcode

  //! \ingroup maxplustruncmat_group
  //!
  //! \brief Function object for multiplication in truncated max-plus
  //! semirings.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `Scalar operator()(Scalar x, Scalar y) const noexcept`
  //! that returns \f$x \otimes y\f$ which is defined by
  //! \f[
  //!    x\otimes y =
  //!    \begin{cases}
  //!    \min\{x + y, T\}   & \text{if } x \neq -\infty\text{ and }y \neq
  //!    -\infty \\ \mbox{}
  //!    -\infty & \text{if } x = -\infty \text{ or }y = -\infty;
  //!    \end{cases}
  //! \f]
  //! representing multiplication in the quotient of the max-plus semiring by
  //! the congruence \f$T = T + 1\f$.
  //!
  //! \tparam T the threshold (point at which the entries in the max-plus
  //! semiring are truncated).
  //! \tparam Scalar the type of the values in the semiring (must be signed
  //! integer type).
  template <size_t T, typename Scalar>
  struct MaxPlusTruncProd {
    static_assert(std::is_signed<Scalar>::value,
                  "MaxPlus requires a signed integer type as parameter!");
    //! \brief Call operator for multiplication.
    //!
    //! This function returns the product of its arguments in a max-plus
    //! truncated semiring.
    //!
    //! \param x the first value.
    //! \param y the second value.
    //!
    //! \returns The product of \p x and \p y in truncated max-plus semiring.
    //!
    //! \exceptions
    //! \noexcept
    Scalar operator()(Scalar x, Scalar y) const noexcept {
      LIBSEMIGROUPS_ASSERT((x >= 0 && x <= static_cast<Scalar>(T))
                           || x == NEGATIVE_INFINITY);
      LIBSEMIGROUPS_ASSERT((y >= 0 && y <= static_cast<Scalar>(T))
                           || y == NEGATIVE_INFINITY);
      if (x == NEGATIVE_INFINITY || y == NEGATIVE_INFINITY) {
        return NEGATIVE_INFINITY;
      }
      return std::min(x + y, static_cast<Scalar>(T));
    }
  };

  //! \ingroup maxplustruncmat_group
  //!
  //! \brief Class representing a truncated max-plus semiring.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This class represents the **truncated max-plus semiring** consists of
  //! the integers \f$\{0, \ldots , t\}\f$ for some value \f$t\f$ (called the
  //! **threshold** of the semiring) and \f$-\infty\f$. Instances of this
  //! class can be used to define the value of the threshold \f$t\f$ at run
  //! time.
  //!
  //! \tparam Scalar the type of the elements of the semiring. This must be a
  //! signed integer type (defaults to \c int).
  template <typename Scalar = int>
  class MaxPlusTruncSemiring {
    static_assert(std::is_signed<Scalar>::value,
                  "MaxPlus requires a signed integer type as parameter!");

   public:
    //! \brief Deleted default constructor.
    //!
    //! Deleted default constructor.
    MaxPlusTruncSemiring() = delete;

    //! \brief Default copy constructor.
    //!
    //! Default copy constructor.
    MaxPlusTruncSemiring(MaxPlusTruncSemiring const&) noexcept = default;

    //! \brief Default move constructor.
    //!
    //! Default move constructor.
    MaxPlusTruncSemiring(MaxPlusTruncSemiring&&) noexcept = default;

    //! \brief Default copy assignment operator.
    //!
    //! Default copy assignment operator.
    MaxPlusTruncSemiring& operator=(MaxPlusTruncSemiring const&) noexcept
        = default;

    //! \brief Default move assignment operator.
    //!
    //! Default move assignment operator.
    MaxPlusTruncSemiring& operator=(MaxPlusTruncSemiring&&) noexcept = default;

    ~MaxPlusTruncSemiring() = default;

    //! \brief Construct from threshold.
    //!
    //! Construct from threshold.
    //!
    //! \param threshold the threshold.
    //!
    //! \throws LibsemigroupsException if `threshold` is less than zero.
    //!
    //! \complexity
    //! Constant.
    explicit MaxPlusTruncSemiring(Scalar threshold) : _threshold(threshold) {
      if (threshold < 0) {
        LIBSEMIGROUPS_EXCEPTION("expected non-negative value, found {}",
                                threshold);
      }
    }

    //! \brief Get the multiplicative identity.
    //!
    //! This function returns the multiplicative identity in a truncated
    //! max-plus semiring.
    //!
    //! \returns The additive identity in a truncated max-plus semiring (the
    //! value \c 0).
    //!
    //! \exceptions
    //! \noexcept
    static constexpr Scalar scalar_one() noexcept {
      return 0;
    }

    //! \brief Get the additive identity.
    //!
    //! This function returns the additive identity in a truncated max-plus
    //! semiring.
    //!
    //! \returns The additive identity in a truncated max-plus semiring (the
    //! value `NEGATIVE_INFINITY`).
    //!
    //! \exceptions
    //! \noexcept
    static constexpr Scalar scalar_zero() noexcept {
      return NEGATIVE_INFINITY;
    }

    //! \brief Multiplication in a truncated max-plus semiring.
    //!
    //! This function returns \f$x \otimes y\f$ which is defined by
    //! \f[
    //!    x\otimes y =
    //!    \begin{cases}
    //!    \min\{x + y, t\}   & \text{if } x \neq -\infty\text{ and }y \neq
    //!    -\infty \\ \mbox{}
    //!    -\infty & \text{if } x = -\infty \text{ or }y = -\infty;
    //!    \end{cases}
    //! \f]
    //! where \f$t\f$ is the threshold; representing multiplication in the
    //! quotient of the max-plus semiring.
    //!
    //! \param x scalar.
    //! \param y scalar.
    //!
    //! \returns A value of type `Scalar`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    Scalar product_no_checks(Scalar x, Scalar y) const noexcept {
      LIBSEMIGROUPS_ASSERT((x >= 0 && x <= _threshold)
                           || x == NEGATIVE_INFINITY);
      LIBSEMIGROUPS_ASSERT((y >= 0 && y <= _threshold)
                           || y == NEGATIVE_INFINITY);
      if (x == NEGATIVE_INFINITY || y == NEGATIVE_INFINITY) {
        return NEGATIVE_INFINITY;
      }
      return std::min(x + y, _threshold);
    }

    //! \brief Addition in a truncated max-plus semiring.
    //!
    //! Returns \f$x \oplus y\f$ which is defined by
    //! \f[
    //!   x\oplus y =
    //!   \begin{cases}
    //!     \max\{x, y\}   & \text{if } x \neq -\infty\text{ and }y \neq
    //!     -\infty
    //!     \\ \mbox{}
    //!     -\infty & \text{if } x = -\infty \text{ or }y = -\infty;
    //!   \end{cases}
    //! \f]
    //! representing addition in the max-plus semiring (and its quotient).
    //!
    //! \param x scalar.
    //! \param y scalar.
    //!
    //! \returns A value of type `Scalar`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    Scalar plus_no_checks(Scalar x, Scalar y) const noexcept {
      LIBSEMIGROUPS_ASSERT((x >= 0 && x <= _threshold)
                           || x == NEGATIVE_INFINITY);
      LIBSEMIGROUPS_ASSERT((y >= 0 && y <= _threshold)
                           || y == NEGATIVE_INFINITY);
      if (x == NEGATIVE_INFINITY) {
        return y;
      } else if (y == NEGATIVE_INFINITY) {
        return x;
      }
      return std::max(x, y);
    }

    //! \brief Get the threshold.
    //!
    //! Returns the threshold value used to construct \ref
    //! MaxPlusTruncSemiring instance.
    //!
    //! \returns A value of type `Scalar`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    Scalar threshold() const noexcept {
      return _threshold;
    }

   public:
    Scalar const _threshold;
  };

  //! \ingroup maxplustruncmat_group
  //!
  //! \brief Alias for dynamic truncated max-plus matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for the type of dynamic truncated max-plus matrices where the
  //! dimension is defined at run time, but the threshold is defined at
  //! compile-time.
  //!
  //! \tparam T the threshold.
  //! \tparam Scalar the type of the entries in the matrix.
  template <size_t T, typename Scalar>
  using DynamicMaxPlusTruncMat = DynamicMatrix<MaxPlusPlus<Scalar>,
                                               MaxPlusTruncProd<T, Scalar>,
                                               MaxPlusZero<Scalar>,
                                               IntegerZero<Scalar>,
                                               Scalar>;

  //! \ingroup maxplustruncmat_group
  //!
  //! \brief Alias for static truncated max-plus matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for static truncated max-plus matrices where the threshold and
  //! dimensions are defined at compile-time.
  //!
  //! \tparam T  the threshold.
  //! \tparam R  the number of rows.
  //! \tparam C  the number of columns.
  //! \tparam Scalar the type of the entries in the matrix.
  template <size_t T, size_t R, size_t C, typename Scalar>
  using StaticMaxPlusTruncMat = StaticMatrix<MaxPlusPlus<Scalar>,
                                             MaxPlusTruncProd<T, Scalar>,
                                             MaxPlusZero<Scalar>,
                                             IntegerZero<Scalar>,
                                             R,
                                             C,
                                             Scalar>;
  //! \ingroup maxplustruncmat_group
  //!
  //! \brief Alias template for truncated max-plus matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias template for truncated max-plus matrices.
  //!
  //! \tparam T the threshold. A value of \c 0 indicates that the value will
  //! be set at run time (default: \c 0).
  //!
  //! \tparam R the number of rows.  A value of \c 0 indicates that the value
  //! will be set at run time (default: \c 0).
  //!
  //! \tparam C the number of columns.  A value of \c 0 indicates that the
  //! value will be set at run time (default: \c R).
  //!
  //! \tparam Scalar the type of the entries in the matrix (default: `int`).
  template <size_t T = 0, size_t R = 0, size_t C = R, typename Scalar = int>
  using MaxPlusTruncMat = std::conditional_t<
      R == 0 || C == 0,
      std::conditional_t<T == 0,
                         DynamicMatrix<MaxPlusTruncSemiring<Scalar>, Scalar>,
                         DynamicMaxPlusTruncMat<T, Scalar>>,
      StaticMaxPlusTruncMat<T, R, C, Scalar>>;

  namespace detail {
    template <typename T>
    struct IsMaxPlusTruncMatHelper : std::false_type {};

    template <size_t T, size_t R, size_t C, typename Scalar>
    struct IsMaxPlusTruncMatHelper<StaticMaxPlusTruncMat<T, R, C, Scalar>>
        : std::true_type {
      static constexpr Scalar threshold = T;
    };

    template <size_t T, typename Scalar>
    struct IsMaxPlusTruncMatHelper<DynamicMaxPlusTruncMat<T, Scalar>>
        : std::true_type {
      static constexpr Scalar threshold = T;
    };

    template <typename Scalar>
    struct IsMaxPlusTruncMatHelper<
        DynamicMatrix<MaxPlusTruncSemiring<Scalar>, Scalar>> : std::true_type {
      static constexpr Scalar threshold = UNDEFINED;
    };
  }  // namespace detail

  //! \ingroup maxplustruncmat_group
  //!
  //! \brief Helper to check if a type is \ref MaxPlusTruncMat.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This variable has value \c true if the template parameter \c T is the
  //! same as \ref MaxPlusTruncMat for some template parameters; and \c false
  //! if it is not.
  //!
  //! \tparam T the type to check.
  template <typename T>
  static constexpr bool IsMaxPlusTruncMat
      = detail::IsMaxPlusTruncMatHelper<T>::value;

  namespace detail {
    template <typename T>
    struct IsTruncMatHelper<T, std::enable_if_t<IsMaxPlusTruncMat<T>>>
        : std::true_type {
      static constexpr typename T::scalar_type threshold
          = IsMaxPlusTruncMatHelper<T>::threshold;
    };
  }  // namespace detail

  namespace matrix {
    //! \ingroup maxplustruncmat_group
    //!
    //! \brief Validate that a truncated max-plus matrix is valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function can be used to validate that a matrix contains values in
    //! the underlying semiring.
    //!
    //! \tparam Mat the type of \p m, must satisfy \ref
    //! IsMaxPlusTruncMat<Mat>.
    //!
    //! \param m the matrix to validate.
    //!
    //! \throws LibsemigroupsException if any entry in the matrix is not in
    //! the set \f$\{0, 1, \ldots, t, -\infty\}\f$ where \f$t\f$ is the
    //! threshold of the matrix or if the underlying semiring is not defined
    //! (only applies to matrices with run time arithmetic).
    template <typename Mat>
    std::enable_if_t<IsMaxPlusTruncMat<Mat>> throw_if_bad_entry(Mat const& m) {
      // TODO(1) to tpp
      detail::throw_if_semiring_nullptr(m);

      using scalar_type   = typename Mat::scalar_type;
      scalar_type const t = matrix::threshold(m);
      auto it = std::find_if_not(m.cbegin(), m.cend(), [t](scalar_type x) {
        return x == NEGATIVE_INFINITY || (0 <= x && x <= t);
      });
      if (it != m.cend()) {
        auto [r, c] = m.coords(it);
        LIBSEMIGROUPS_EXCEPTION(
            "invalid entry, expected values in {{0, 1, ..., {}, {}}} "
            "but found {} in entry ({}, {})",
            t,
            NEGATIVE_INFINITY,
            detail::entry_repr(*it),
            r,
            c);
      }
    }

    //! \ingroup maxplustruncmat_group
    //!
    //! \brief Validate that an entry in a truncated max-plus matrix is valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function can be used to validate that a matrix contains values in
    //! the underlying semiring.
    //!
    //! \tparam Mat the type the 1st argument, must satisfy \ref
    //! IsMaxPlusTruncMat<Mat>.
    //!
    //! \param m the matrix.
    //! \param val the entry to validate.
    //!
    //! \throws LibsemigroupsException if the entry \p val is not in
    //! the set \f$\{0, 1, \ldots, t, -\infty\}\f$ where \f$t\f$ is the
    //! threshold of the matrix or if the underlying semiring is not defined
    //! (only applies to matrices with run time arithmetic).
    template <typename Mat>
    std::enable_if_t<IsMaxPlusTruncMat<Mat>>
    throw_if_bad_entry(Mat const& m, typename Mat::scalar_type val) {
      detail::throw_if_semiring_nullptr(m);
      using scalar_type   = typename Mat::scalar_type;
      scalar_type const t = matrix::threshold(m);
      if (val == POSITIVE_INFINITY || 0 > val || val > t) {
        LIBSEMIGROUPS_EXCEPTION(
            "invalid entry, expected values in {{0, 1, ..., {}, -{}}} "
            "but found {}",
            t,
            NEGATIVE_INFINITY,
            detail::entry_repr(val));
      }
    }
  }  // namespace matrix

  ////////////////////////////////////////////////////////////////////////
  // Min-plus matrices with threshold
  ////////////////////////////////////////////////////////////////////////

  //! \defgroup minplustruncmat_group Truncated min-plus matrices
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This page describes the functionality for \f$n \times n\f$  matrices
  //! over the finite quotient of the min-plus semiring by the congruence
  //! \f$t = t + 1\f$ for arbitrary \f$n\f$ and \f$t\f$. The value
  //! \f$t\f$ is referred to as the *threshold*.
  //!
  //! There are three types of such matrices where:
  //!
  //! 1. the dimension is known at compile-time;
  //! 2. the dimension is to be defined a run time but the arithmetic
  //! operations
  //!    are known at compile-time (i.e. the value of \f$t\f$ is known at
  //!    compile time)
  //! 3. both the dimension and the arithmetic operations (i.e. \f$t\f$) are
  //!    to be defined a run time.
  //!
  //! All three of these types can be accessed via the alias template \ref
  //! MinPlusTruncMat<T, P, R, C, Scalar> if \c T has value \c 0, then the
  //! threshold can be set at run time, and if \c R or \c C is \c 0, then the
  //! dimension can be set at run time. The default value of \c T is \c 0, \c
  //! R is \c 0, and of \c C is \c R.
  //!
  //! The alias \ref MinPlusTruncMat is either StaticMatrix, \ref
  //! DynamicMatrixStaticArith "DynamicMatrix (compile-time arithmetic)", or
  //! \ref DynamicMatrixDynamicArith "DynamicMatrix (run-time arithmetic)",
  //! please refer to the documentation of these class templates for more
  //! details. The only substantial difference in the interface of static
  //! and dynamic matrices is that the former can be default constructed and the
  //! latter should be constructed using the dimensions.
  //!
  //! \par Example
  //! \code
  //! // construct an uninitialized 3 x 3 static matrix with threshold 11
  //! MinPlusTruncMat<11, 3> m;
  //! // construct an uninitialized 4 x 4 dynamic matrix with threshold 11
  //! MinPlusTruncMat<11> m(4, 4);
  //! // construct a truncated min-plus semiring with threshold 11
  //! MinPlusTruncSemiring sr(11);
  //! // construct an uninitialized 5 x 5 dynamic matrix with threshold 11
  //! // (defined at run time)
  //! MinPlusTruncMat<>  m(sr, 5, 5);
  //! \endcode

  //! \ingroup minplustruncmat_group
  //!
  //! \brief Function object for multiplication in min-plus truncated
  //! semirings.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `Scalar operator()(Scalar x, Scalar y) const noexcept` that returns \f$x
  //! \otimes y\f$ which is defined by
  //! \f[
  //!    x\otimes y =
  //!    \begin{cases}
  //!    \min\{x + y, T\}   & \text{if } x \neq \infty\text{ and }y \neq
  //!    \infty \\ \mbox
  //!    \infty & \text{if } x = \infty \text{ or }y = \infty;
  //!    \end{cases}
  //! \f]
  //! representing multiplication in the quotient of the min-plus semiring
  //! by the congruence \f$T = T + 1\f$.
  //!
  //! \tparam T the threshold (point at which the entries in the min-plus
  //! semiring are truncated).
  //! \tparam Scalar the type of the values in the semiring.
  template <size_t T, typename Scalar>
  struct MinPlusTruncProd {
    //! \brief Call operator for multiplication.
    //!
    //! This function returns the product of its arguments in a truncated
    //! min-plus semiring.
    //!
    //! \param x the first value.
    //! \param y the second value.
    //!
    //! \returns The product of \p x and \p y in truncated min-plus semiring.
    //!
    //! \exceptions
    //! \noexcept
    Scalar operator()(Scalar x, Scalar y) const noexcept {
      LIBSEMIGROUPS_ASSERT((x >= 0 && x <= static_cast<Scalar>(T))
                           || x == POSITIVE_INFINITY);
      LIBSEMIGROUPS_ASSERT((y >= 0 && y <= static_cast<Scalar>(T))
                           || y == POSITIVE_INFINITY);
      if (x == POSITIVE_INFINITY || y == POSITIVE_INFINITY) {
        return POSITIVE_INFINITY;
      }
      return std::min(x + y, static_cast<Scalar>(T));
    }
  };

  //! \ingroup minplustruncmat_group
  //!
  //! \brief Class representing a truncated min-plus semiring.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This class represents the **min-plus truncated semiring** consists of
  //! the integers \f$\{0, \ldots , t\}\f$ for some value \f$t\f$ (called the
  //! **threshold** of the semiring) and \f$\infty\f$. Instances of this class
  //! can be used to define the value of the threshold \f$t\f$ at run time.
  //!
  //! \tparam Scalar the type of the elements of the semiring. This must be an
  //! integral type.
  template <typename Scalar = int>
  class MinPlusTruncSemiring {
    static_assert(std::is_integral<Scalar>::value,
                  "MinPlus requires an integral type as parameter!");

   public:
    //! \brief Deleted default constructor.
    //!
    //! Deleted default constructor.
    MinPlusTruncSemiring() = delete;

    //! \brief Default copy constructor.
    //!
    //! Default copy constructor.
    MinPlusTruncSemiring(MinPlusTruncSemiring const&) noexcept = default;

    //! \brief Default move constructor.
    //!
    //! Default move constructor.
    MinPlusTruncSemiring(MinPlusTruncSemiring&&) noexcept = default;

    //! \brief Default copy assignment operator.
    //!
    //! Default copy assignment operator.
    MinPlusTruncSemiring& operator=(MinPlusTruncSemiring const&) noexcept
        = default;

    //! \brief Default move assignment operator.
    //!
    //! Default move assignment operator.
    MinPlusTruncSemiring& operator=(MinPlusTruncSemiring&&) noexcept = default;

    //! \brief Construct from threshold.
    //!
    //! Construct from threshold.
    //!
    //! \param threshold the threshold.
    //!
    //! \throws LibsemigroupsException if `threshold` is less than zero.
    //!
    //! \complexity
    //! Constant.
    explicit MinPlusTruncSemiring(Scalar threshold) : _threshold(threshold) {
      if (std::is_signed<Scalar>::value && threshold < 0) {
        LIBSEMIGROUPS_EXCEPTION("expected non-negative value, found {}",
                                threshold);
      }
    }

    //! \brief Get the multiplicative identity.
    //!
    //! This function returns the multiplicative identity in a truncated
    //! min-plus semiring.
    //!
    //! \returns The multiplicative identity in a truncated min-plus semiring
    //! (the value \c 0).
    //!
    //! \exceptions
    //! \noexcept
    static constexpr Scalar scalar_one() noexcept {
      return 0;
    }

    //! \brief Get the additive identity.
    //!
    //! This function returns the additive identity in a truncated min-plus
    //! semiring.
    //!
    //! \returns The additive identity in a truncated min-plus semiring (the
    //! value `POSITIVE_INFINITY`).
    //!
    //! \exceptions
    //! \noexcept
    // TODO(1) These mem fns (one and zero) aren't needed?
    static constexpr Scalar scalar_zero() noexcept {
      return POSITIVE_INFINITY;
    }

    //! \brief Multiplication in a truncated min-plus semiring.
    //!
    //! Returns \f$x \otimes y\f$ which is defined by
    //! \f[
    //!    x\otimes y =
    //!    \begin{cases}
    //!    \min\{x + y, t\}   & \text{if } x \neq \infty\text{ and }y \neq
    //!    \infty \\ \mbox{}
    //!    \infty & \text{if } x = \infty \text{ or }y = \infty;
    //!    \end{cases}
    //! \f]
    //! where \f$t\f$ is the threshold; representing multiplication in the
    //! quotient of the min-plus semiring.
    //!
    //! \param x scalar.
    //! \param y scalar.
    //!
    //! \returns A value of type `Scalar`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    Scalar product_no_checks(Scalar x, Scalar y) const noexcept {
      LIBSEMIGROUPS_ASSERT((x >= 0 && x <= _threshold)
                           || x == POSITIVE_INFINITY);
      LIBSEMIGROUPS_ASSERT((y >= 0 && y <= _threshold)
                           || y == POSITIVE_INFINITY);
      if (x == POSITIVE_INFINITY || y == POSITIVE_INFINITY) {
        return POSITIVE_INFINITY;
      }
      return std::min(x + y, _threshold);
    }

    //! \brief Addition in a truncated min-plus semiring.
    //!
    //! Returns \f$x \oplus y\f$ which is defined by
    //! \f[
    //!   x\oplus y =
    //!   \begin{cases}
    //!   \min\{x, y\}   & \text{if } x \neq \infty\text{ and }y \neq \infty
    //!   \\ \mbox{}
    //!   \infty & \text{if } x = \infty \text{ or }y = \infty;
    //!   \end{cases}
    //! \f]
    //! representing addition in the min-plus truncated semiring (and its
    //! quotient).
    //!
    //! \param x scalar.
    //! \param y scalar.
    //!
    //! \returns A value of type `Scalar`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    Scalar plus_no_checks(Scalar x, Scalar y) const noexcept {
      LIBSEMIGROUPS_ASSERT((x >= 0 && x <= _threshold)
                           || x == POSITIVE_INFINITY);
      LIBSEMIGROUPS_ASSERT((y >= 0 && y <= _threshold)
                           || y == POSITIVE_INFINITY);
      if (x == POSITIVE_INFINITY) {
        return y;
      } else if (y == POSITIVE_INFINITY) {
        return x;
      }
      return std::min(x, y);
    }

    //! \brief Get the threshold.
    //!
    //! Returns the threshold value used to construct \ref
    //! MinPlusTruncSemiring instance.
    //!
    //! \returns A value of type `Scalar`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    Scalar threshold() const noexcept {
      return _threshold;
    }

   public:
    Scalar const _threshold;
  };

  //! \ingroup minplustruncmat_group
  //!
  //! \brief Alias for dynamic truncated min-plus matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for the type of dynamic truncated min-plus matrices where the
  //! dimension is defined at run time, but the threshold is defined at
  //! compile-time.
  //!
  //! \tparam T the threshold.
  //! \tparam Scalar the type of the entries in the matrix.
  template <size_t T, typename Scalar>
  using DynamicMinPlusTruncMat = DynamicMatrix<MinPlusPlus<Scalar>,
                                               MinPlusTruncProd<T, Scalar>,
                                               MinPlusZero<Scalar>,
                                               IntegerZero<Scalar>,
                                               Scalar>;

  //! \ingroup minplustruncmat_group
  //!
  //! \brief Alias for static truncated min-plus matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for static truncated min-plus matrices where the threshold and
  //! dimensions are defined at compile-time.
  //!
  //! \tparam T  the threshold.
  //! \tparam R  the number of rows.
  //! \tparam C  the number of columns.
  //! \tparam Scalar the type of the entries in the matrix.
  template <size_t T, size_t R, size_t C, typename Scalar>
  using StaticMinPlusTruncMat = StaticMatrix<MinPlusPlus<Scalar>,
                                             MinPlusTruncProd<T, Scalar>,
                                             MinPlusZero<Scalar>,
                                             IntegerZero<Scalar>,
                                             R,
                                             C,
                                             Scalar>;

  //! \ingroup minplustruncmat_group
  //!
  //! \brief Alias template for truncated min-plus matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias template for truncated min-plus matrices.
  //!
  //! \tparam T the threshold. A value of \c 0 indicates that the value will
  //! be set at run time (default: \c 0).
  //!
  //! \tparam R the number of rows.  A value of \c 0 indicates that the value
  //! will be set at run time (default: \c 0).
  //!
  //! \tparam C the number of columns.  A value of \c 0 indicates that the
  //! value will be set at run time (default: \c R).
  //!
  //! \tparam Scalar the type of the entries in the matrix (default: `int`).
  template <size_t T = 0, size_t R = 0, size_t C = R, typename Scalar = int>
  using MinPlusTruncMat = std::conditional_t<
      R == 0 || C == 0,
      std::conditional_t<T == 0,
                         DynamicMatrix<MinPlusTruncSemiring<Scalar>, Scalar>,
                         DynamicMinPlusTruncMat<T, Scalar>>,
      StaticMinPlusTruncMat<T, R, C, Scalar>>;

  namespace detail {
    template <typename T>
    struct IsMinPlusTruncMatHelper : std::false_type {};

    template <size_t T, size_t R, size_t C, typename Scalar>
    struct IsMinPlusTruncMatHelper<StaticMinPlusTruncMat<T, R, C, Scalar>>
        : std::true_type {
      static constexpr Scalar threshold = T;
    };

    template <size_t T, typename Scalar>
    struct IsMinPlusTruncMatHelper<DynamicMinPlusTruncMat<T, Scalar>>
        : std::true_type {
      static constexpr Scalar threshold = T;
    };

    template <typename Scalar>
    struct IsMinPlusTruncMatHelper<
        DynamicMatrix<MinPlusTruncSemiring<Scalar>, Scalar>> : std::true_type {
      static constexpr Scalar threshold = UNDEFINED;
    };
  }  // namespace detail

  //! \ingroup minplustruncmat_group
  //!
  //! \brief Helper to check if a type is \ref MinPlusTruncMat.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This variable has value \c true if the template parameter \c T is the
  //! same as \ref MinPlusTruncMat for some template parameters; and \c false
  //! if it is not.
  //!
  //! \tparam T the type to check.
  template <typename T>
  static constexpr bool IsMinPlusTruncMat
      = detail::IsMinPlusTruncMatHelper<T>::value;

  namespace detail {
    template <typename T>
    struct IsTruncMatHelper<T, std::enable_if_t<IsMinPlusTruncMat<T>>>
        : std::true_type {
      static constexpr typename T::scalar_type threshold
          = IsMinPlusTruncMatHelper<T>::threshold;
    };
  }  // namespace detail

  namespace matrix {
    //! \ingroup minplustruncmat_group
    //!
    //! \brief Validate that a truncated min-plus matrix is valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function can be used to validate that a matrix contains values in
    //! the underlying semiring.
    //!
    //! \tparam Mat the type of the matrix \p m (must satisfy
    //! \ref IsMinPlusTruncMat<Mat>).
    //!
    //! \param m the matrix to validate.
    //!
    //! \throws LibsemigroupsException if any entry in the matrix is not in
    //! the set \f$\{0, 1, \ldots, t, \infty\}\f$ where \f$t\f$ is the
    //! threshold of the matrix or if the underlying semiring is not defined
    //! (only applies to matrices with run time arithmetic).
    // TODO(1) to tpp
    template <typename Mat>
    std::enable_if_t<IsMinPlusTruncMat<Mat>> throw_if_bad_entry(Mat const& m) {
      // Check that the semiring pointer isn't the nullptr if it shouldn't be
      detail::throw_if_semiring_nullptr(m);

      using scalar_type   = typename Mat::scalar_type;
      scalar_type const t = matrix::threshold(m);
      auto it = std::find_if_not(m.cbegin(), m.cend(), [t](scalar_type x) {
        return x == POSITIVE_INFINITY || (0 <= x && x <= t);
      });
      if (it != m.cend()) {
        uint64_t r, c;
        std::tie(r, c) = m.coords(it);

        LIBSEMIGROUPS_EXCEPTION(
            "invalid entry, expected values in {{0, 1, ..., {}, {}}} "
            "but found {} in entry ({}, {})",
            t,
            u8"\u221E",
            detail::entry_repr(*it),
            r,
            c);
      }
    }

    //! \ingroup minplustruncmat_group
    //!
    //! \brief Validate that an entry in a truncated min-plus matrix is valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function can be used to validate that a matrix contains values in
    //! the underlying semiring.
    //!
    //! \tparam Mat the type the 1st argument, must satisfy \ref
    //! IsMinPlusTruncMat<Mat>.
    //!
    //! \param m the matrix.
    //! \param val the entry to validate.
    //!
    //! \throws LibsemigroupsException if the entry \p val is not in
    //! the set \f$\{0, 1, \ldots, t, \infty\}\f$ where \f$t\f$ is the
    //! threshold of the matrix or if the underlying semiring is not defined
    //! (only applies to matrices with run time arithmetic).
    template <typename Mat>
    std::enable_if_t<IsMinPlusTruncMat<Mat>>
    throw_if_bad_entry(Mat const& m, typename Mat::scalar_type val) {
      detail::throw_if_semiring_nullptr(m);

      using scalar_type   = typename Mat::scalar_type;
      scalar_type const t = matrix::threshold(m);
      if (!(val == POSITIVE_INFINITY || (0 <= val && val <= t))) {
        LIBSEMIGROUPS_EXCEPTION(
            "invalid entry, expected values in {{0, 1, ..., {}, {}}} "
            "but found {}",
            t,
            u8"\u221E",
            detail::entry_repr(val));
      }
    }
  }  // namespace matrix

  ////////////////////////////////////////////////////////////////////////
  // NTP matrices
  ////////////////////////////////////////////////////////////////////////

  //! \defgroup ntpmat_group Matrices over the natural numbers mod (t = t + p)
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This page describes the functionality for \f$n \times n\f$ matrices over
  //! the finite quotient of the usual semiring of natural number by the
  //! congruence \f$t = t + p\f$ for arbitrary \f$n\f$, \f$t\f$, and
  //! \f$p\f$. The value \f$t\f$ is referred to as the *threshold* and
  //! \f$p\f$ is called the *period*. The matrices of this type are referred
  //! to by the acronym **ntp** matrices, for "natural threshold period".  The
  //! NTPSemiring has elements \f$\{0, 1, ..., t, t +  1, ..., t
  //! + p - 1\}\f$ where  \f$t\f$, and \f$p\f$ are the threshold and period,
  //! respectively; addition and multiplication in the NTPSemiring is
  //! defined below.
  //!
  //! There are three types of such matrices where:
  //!
  //! 1. the dimension is known at compile-time;
  //! 2. the dimension is to be defined a run time but the arithmetic
  //!    operations are known at compile-time (i.e. the values of
  //!    \f$t\f$ and \f$p\f$ are known at compile time)
  //! 3. both the dimension and the arithmetic operations (i.e. \f$t\f$ and
  //!    \f$p\f$) are to be defined a run time.
  //!
  //! All three of these types can be accessed via the alias template \ref
  //! NTPMat<T, P, R, C, Scalar> if \c T and \c P have value \c 0, then the
  //! threshold and period can be set at run time, and if \c R or \c C is \c
  //! 0, then the dimension can be set at run time.  The default values of \c
  //! T, \c P, and \c R are \c 0, and the default value of \c C is \c R.
  //!
  //! The alias \ref NTPMat is one of StaticMatrix, \ref
  //! DynamicMatrixStaticArith "DynamicMatrix (compile-time arithmetic)", or
  //! \ref DynamicMatrixDynamicArith "DynamicMatrix (run-time arithmetic)";
  //! please refer to the documentation of these class templates for more
  //! details. The only substantial difference in the interface of static and
  //! dynamic matrices is that the former can be default constructed and the
  //! latter should be constructed using the dimensions.
  //!
  //! \par Example
  //! \code
  //! // construct an uninitialized 3 x 3 static matrix with threshold
  //! // 11, period 2
  //! NTPMat<11, 2, 3> m;
  //! // construct an uninitialized 4 x 4 dynamic matrix with threshold 11,
  //! // period 2
  //! NTPMat<11, 2> m(4, 4);
  //! // construct an ntp semiring with threshold 11, period 2
  //! NTPSemiring sr(11, 2);
  //! // construct an uninitialized 5 x 5 dynamic matrix with threshold 11,
  //! // period 2
  //! NTPMat<> m(sr, 5, 5);
  //! \endcode

  namespace detail {
    template <size_t T, size_t P, typename Scalar>
    constexpr Scalar thresholdperiod(Scalar x) noexcept {
      if (x > T) {
        return T + (x - T) % P;
      }
      return x;
    }

    template <typename Scalar>
    inline Scalar thresholdperiod(Scalar x,
                                  Scalar threshold,
                                  Scalar period) noexcept {
      if (x > threshold) {
        return threshold + (x - threshold) % period;
      }
      return x;
    }
  }  // namespace detail

  //! \ingroup ntpmat_group
  //!
  //! \brief Function object for addition in ntp semirings.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `Scalar operator()(Scalar x, Scalar y) const noexcept`
  //! that returns \f$x \oplus y\f$ which is defined by
  //! \f[
  //!    x\oplus y =
  //!    \begin{cases}
  //!    x + y & \text{if } x + y \leq T \\ \mbox{}
  //!    T + ((x + y) - T \pmod{P}) & \text{if } x + y > T
  //!    \end{cases}
  //! \f]
  //! representing addition in the quotient of the semiring natural numbers by
  //! the congruence \f$(T = T + P)\f$.
  //!
  //! \tparam T the threshold.
  //! \tparam P the period.
  //! \tparam Scalar the type of the values in the semiring.
  // Static arithmetic
  template <size_t T, size_t P, typename Scalar>
  struct NTPPlus {
    //! \brief Call operator for addition.
    //!
    //! This function returns the sum of its arguments in an ntp semiring.
    //!
    //! \param x the first value.
    //! \param y the second value.
    //!
    //! \returns The sum of \p x and \p y in an ntp semiring.
    //!
    //! \exceptions
    //! \noexcept
    constexpr Scalar operator()(Scalar x, Scalar y) const noexcept {
      return detail::thresholdperiod<T, P>(x + y);
    }
  };

  //! \ingroup ntpmat_group
  //!
  //! \brief Function object for multiplication in an ntp semirings.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This is a stateless struct with a single call operator of signature:
  //! `Scalar operator()(Scalar x, Scalar y) const noexcept`
  //! that returns \f$x \otimes y\f$ which is defined by
  //!
  //! \f[
  //!   x\otimes y =
  //!   \begin{cases}
  //!   xy & \text{if } xy \leq T \\ \mbox{}
  //!   T + ((xy - T) \pmod{P}) & \text{if } xy > T
  //!   \end{cases}
  //! \f]
  //! representing multiplication in the quotient of the semiring natural
  //! numbers by the congruence \f$(T = T + P)\f$.
  //!
  //! \tparam T the threshold.
  //! \tparam P the period.
  //! \tparam Scalar the type of the values in the semiring.
  template <size_t T, size_t P, typename Scalar>
  struct NTPProd {
    //! \brief Call operator for multiplication.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function returns the product of its arguments in an ntp semiring.
    //!
    //! \param x the first value.
    //! \param y the second value.
    //!
    //! \returns The product of \p x and \p y in an ntp semiring.
    //!
    //! \exceptions
    //! \noexcept
    constexpr Scalar operator()(Scalar x, Scalar y) const noexcept {
      return detail::thresholdperiod<T, P>(x * y);
    }
  };

  //! \ingroup ntpmat_group
  //!
  //! \brief Class representing an ntp semiring.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This class represents the **ntp semiring** consists of the integers
  //! \f$\{0, 1, ..., t, t +  1, ..., t + p - 1\}\f$ for some  \f$t\f$ and
  //! \f$p\f$ (called the **threshold** and **period**). Instances of this
  //! class can be used to define the value of the threshold \f$t\f$ and
  //! period \f$p\f$ at run time.
  //!
  //! \tparam Scalar the type of the elements of the semiring.
  // Dynamic arithmetic
  template <typename Scalar = size_t>
  class NTPSemiring {
   public:
    //! \brief Deleted default constructor.
    //!
    //! Deleted default constructor.
    // Deleted to avoid uninitialised values of period and threshold.
    NTPSemiring() = delete;

    //! \brief Default copy constructor.
    //!
    //! Default copy constructor.
    NTPSemiring(NTPSemiring const&) = default;

    //! \brief Default move constructor.
    //!
    //! Default move constructor.
    NTPSemiring(NTPSemiring&&) = default;

    //! \brief Default copy assignment operator.
    //!
    //! Default copy assignment operator.
    NTPSemiring& operator=(NTPSemiring const&) = default;

    //! \brief Default move assignment operator.
    //!
    //! Default move assignment operator.
    NTPSemiring& operator=(NTPSemiring&&) = default;

    ~NTPSemiring() = default;

    //! \brief Construct from threshold and period.
    //!
    //! Construct from threshold and period.
    //!
    //! \param t the threshold (\f$t \geq 0\f$).
    //! \param p the period (\f$p > 0\f$).
    //!
    //! \throws LibsemigroupsException if \p t is less than zero.
    //! \throws LibsemigroupsException if \p p is less than or equal to zero.
    //!
    //! \complexity
    //! Constant.
    NTPSemiring(Scalar t, Scalar p) : _period(p), _threshold(t) {
      if constexpr (std::is_signed<Scalar>::value) {
        if (t < 0) {
          LIBSEMIGROUPS_EXCEPTION(
              "expected non-negative value for 1st argument, found {}", t);
        }
      }
      if (p <= 0) {
        LIBSEMIGROUPS_EXCEPTION(
            "expected positive value for 2nd argument, found {}", p);
      }
    }

    //! \brief Get the multiplicative identity.
    //!
    //! This function returns \f$1\f$; representing the multiplicative
    //! identity of the quotient of the semiring of natural numbers.
    //!
    //! \returns The multiplicative identity in an ntp semiring (the
    //! value \c 0).
    //!
    //! \exceptions
    //! \noexcept
    static constexpr Scalar scalar_one() noexcept {
      return 1;
    }

    //! \brief Get the additive identity.
    //!
    //! This function returns \f$0\f$ representing the additive identity of
    //! the quotient of the semiring of natural numbers.
    //!
    //! \returns A value of type `Scalar`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    static constexpr Scalar scalar_zero() noexcept {
      return 0;
    }

    //! \brief Multiplication in an ntp semiring.
    //!
    //! This function returns \f$x \otimes y\f$ which is defined by
    //! \f[
    //!   x\otimes y =
    //!   \begin{cases}
    //!   xy & \text{if } xy \leq t \\ \mbox{}
    //!   t + ((xy - t)\pmod{p}) & \text{if } xy > t
    //!   \end{cases}
    //! \f]
    //! where \f$t\f$ is the \ref threshold and \f$p\f$ is the \ref period;
    //! representing multiplication in the quotient of the semiring of natural
    //! numbers.
    //!
    //! \param x scalar (`0 <= x < threshold() + period()`).
    //! \param y scalar (`0 <= y < threshold() + period()`).
    //!
    //! \returns A value of type `Scalar`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    Scalar product_no_checks(Scalar x, Scalar y) const noexcept {
      LIBSEMIGROUPS_ASSERT(x >= 0 && x <= _period + _threshold - 1);
      LIBSEMIGROUPS_ASSERT(y >= 0 && y <= _period + _threshold - 1);
      return detail::thresholdperiod(x * y, _threshold, _period);
    }

    //! \brief Addition in an ntp semiring.
    //!
    //! This function returns \f$x \oplus y\f$ which is defined by
    //! \f[
    //!   x\oplus y =
    //!   \begin{cases}
    //!   x + y                      & \text{if } x + y \leq t \\ \mbox{}
    //!   t + ((x + y - t) \pmod{p}) & \text{if } x + y > t
    //!   \end{cases}
    //! \f]
    //! where \f$t\f$ is the \ref threshold and \f$p\f$ is the \ref period;
    //! representing multiplication in the quotient of the semiring of natural
    //! numbers.
    //!
    //! \param x scalar (\f$0\leq x < t + p\f$).
    //! \param y scalar (\f$0\leq y < t + p\f$).
    //!
    //! \returns A value of type `Scalar`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    Scalar plus_no_checks(Scalar x, Scalar y) const noexcept {
      LIBSEMIGROUPS_ASSERT(x >= 0 && x <= _period + _threshold - 1);
      LIBSEMIGROUPS_ASSERT(y >= 0 && y <= _period + _threshold - 1);
      return detail::thresholdperiod(x + y, _threshold, _period);
    }

    //! \brief Get the threshold.
    //!
    //! Returns the threshold value used to construct \ref NTPSemiring
    //! instance.
    //!
    //! \returns A value of type `Scalar`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    Scalar threshold() const noexcept {
      return _threshold;
    }

    //! \brief Get the period.
    //!
    //! Returns the period value used to construct \ref NTPSemiring
    //! instance.
    //!
    //! \returns A value of type `Scalar`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    Scalar period() const noexcept {
      return _period;
    }

   private:
    Scalar _period;
    Scalar _threshold;
  };

  //! \ingroup ntpmat_group
  //!
  //! \brief Alias for ntp matrices with dynamic threshold and period.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for ntp matrices with dimensions, threshold, and period defined
  //! at runtime.
  //!
  //! \tparam Scalar the type of the entries in the matrix.
  template <typename Scalar>
  using DynamicNTPMatWithSemiring = DynamicMatrix<NTPSemiring<Scalar>, Scalar>;

  //! \ingroup ntpmat_group
  //!
  //! \brief Alias for ntp matrices with static threshold and period.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for the type of dynamic ntp matrices where the dimension
  //! is defined at run time, but the threshold and period are defined at
  //! compile-time.
  //!
  //! \tparam T       the threshold.
  //! \tparam P       the period.
  //! \tparam Scalar  the type of the entries in the matrix.
  template <size_t T, size_t P, typename Scalar>
  using DynamicNTPMatWithoutSemiring = DynamicMatrix<NTPPlus<T, P, Scalar>,
                                                     NTPProd<T, P, Scalar>,
                                                     IntegerZero<Scalar>,
                                                     IntegerOne<Scalar>,
                                                     Scalar>;

  //! \ingroup ntpmat_group
  //!
  //! \brief Alias for ntp matrices with static threshold and period, and
  //! dimensions.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for static ntp matrices where the threshold, period, and
  //! dimensions are defined at compile-time.
  //!
  //! \tparam T  the threshold.
  //!
  //! \tparam P  the period.
  //!
  //! \tparam R  the number of rows.
  //!
  //! \tparam C  the number of columns.
  //!
  //! \tparam Scalar the type of the entries in the matrix (default: `int`).
  template <size_t T, size_t P, size_t R, size_t C, typename Scalar>
  using StaticNTPMat = StaticMatrix<NTPPlus<T, P, Scalar>,
                                    NTPProd<T, P, Scalar>,
                                    IntegerZero<Scalar>,
                                    IntegerOne<Scalar>,
                                    R,
                                    C,
                                    Scalar>;

  //! \ingroup ntpmat_group
  //!
  //! \brief Alias template for ntp matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias template for ntp matrices.
  //!
  //! \tparam T the threshold. If both \c T and \c P are \c 0, this indicates
  //! that the value will be set at run time (default: \c 0).
  //!
  //! \tparam P the period. If both \c T and \c P are \c 0, this indicates
  //! that the value will be set at run time (default: \c 0).
  //!
  //! \tparam R the number of rows.  A value of \c 0 indicates that the value
  //! will be set at run time (default: \c 0).
  //!
  //! \tparam C the number of columns.  A value of \c 0 indicates that the
  //! value will be set at run time (default: \c R).
  //!
  //! \tparam Scalar the type of the entries in the matrix (default:
  //! `size_t`).
  template <size_t T        = 0,
            size_t P        = 0,
            size_t R        = 0,
            size_t C        = R,
            typename Scalar = size_t>
  using NTPMat = std::conditional_t<
      R == 0 || C == 0,
      std::conditional_t<T == 0 && P == 0,
                         DynamicNTPMatWithSemiring<Scalar>,
                         DynamicNTPMatWithoutSemiring<T, P, Scalar>>,
      StaticNTPMat<T, P, R, C, Scalar>>;

  namespace detail {
    template <typename T>
    struct IsNTPMatHelper : std::false_type {};

    template <typename Scalar>
    struct IsNTPMatHelper<DynamicNTPMatWithSemiring<Scalar>> : std::true_type {
      static constexpr Scalar threshold = UNDEFINED;
      static constexpr Scalar period    = UNDEFINED;
    };

    template <size_t T, size_t P, typename Scalar>
    struct IsNTPMatHelper<DynamicNTPMatWithoutSemiring<T, P, Scalar>>
        : std::true_type {
      static constexpr Scalar threshold = T;
      static constexpr Scalar period    = P;
    };

    template <size_t T, size_t P, size_t R, size_t C, typename Scalar>
    struct IsNTPMatHelper<StaticNTPMat<T, P, R, C, Scalar>> : std::true_type {
      static constexpr Scalar threshold = T;
      static constexpr Scalar period    = P;
    };
  }  // namespace detail

  //! \ingroup ntpmat_group
  //!
  //! \brief Helper to check if a type is \ref NTPMat.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This variable has value \c true if the template parameter \c U is the
  //! same as \ref NTPMat<T, P, R, C, Scalar> for some values of \c T, \c P,
  //! \c R, \c C, and `Scalar`; and \c false if it is not.
  //!
  //! \tparam U the type to check.
  template <typename U>
  static constexpr bool IsNTPMat = detail::IsNTPMatHelper<U>::value;

  namespace detail {
    template <typename T>
    struct IsTruncMatHelper<T, std::enable_if_t<IsNTPMat<T>>> : std::true_type {
      static constexpr typename T::scalar_type threshold
          = IsNTPMatHelper<T>::threshold;
      static constexpr typename T::scalar_type period
          = IsNTPMatHelper<T>::period;
    };

  }  // namespace detail

  namespace matrix {
    //! \brief Returns the period of a static ntp matrix.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function returns the template parameter \p P of a static ntp
    //! matrix.
    //!
    //! \tparam T  the threshold.
    //!
    //! \tparam P  the period.
    //!
    //! \tparam R  the number of rows.
    //!
    //! \tparam C  the number of columns.
    //!
    //! \tparam Scalar the type of the entries in the matrix.
    //!
    //! \returns The template parameter \p P.
    //!
    //! \exceptions
    //! \noexcept
    template <size_t T, size_t P, size_t R, size_t C, typename Scalar>
    constexpr Scalar period(StaticNTPMat<T, P, R, C, Scalar> const&) noexcept {
      return P;
    }

    //! \brief Returns the period of a dynamic ntp matrix.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function returns the template parameter \p P of a dynamic ntp
    //! matrix.
    //!
    //! \tparam T  the threshold.
    //!
    //! \tparam P  the period.
    //!
    //! \tparam Scalar the type of the entries in the matrix.
    //!
    //! \returns The template parameter \p P.
    //!
    //! \exceptions
    //! \noexcept
    template <size_t T, size_t P, typename Scalar>
    constexpr Scalar
    period(DynamicNTPMatWithoutSemiring<T, P, Scalar> const&) noexcept {
      return P;
    }

    //! \brief Returns the period of a dynamic ntp matrix.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function returns the period of the dynamic ntp
    //! matrix \p x using its underlying semiring.
    //!
    //! \tparam Scalar the type of the entries in the matrix.
    //!
    //! \param x the dynamic ntp matrix.
    //!
    //! \returns The period of the matrix \p x.
    //!
    //! \exceptions
    //! \noexcept
    template <typename Scalar>
    Scalar period(DynamicNTPMatWithSemiring<Scalar> const& x) noexcept {
      return x.semiring()->period();
    }
  }  // namespace matrix

  namespace matrix {
    //! \ingroup ntpmat_group
    //!
    //! \brief Check that the entries in an ntp matrix are valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function can be used to validate that a matrix contains values in
    //! the underlying semiring.
    //!
    //! \tparam Mat the type of the matrix \p m (must satisfy
    //! \ref IsNTPMat<Mat>).
    //!
    //! \param m the matrix to check.
    //!
    //! \throws LibsemigroupsException if any entry in the matrix is not in
    //! the set \f$\{0, 1, \ldots, t + p - 1\}\f$ where \f$t\f$ is the
    //! threshold, and \f$p\f$ the period, of the matrix.
    //!
    //! \throws LibsemigroupsException if the underlying semiring is not
    //! defined (only applies to matrices with run time arithmetic).
    template <typename Mat>
    std::enable_if_t<IsNTPMat<Mat>> throw_if_bad_entry(Mat const& m) {
      detail::throw_if_semiring_nullptr(m);

      using scalar_type   = typename Mat::scalar_type;
      scalar_type const t = matrix::threshold(m);
      scalar_type const p = matrix::period(m);
      auto it = std::find_if_not(m.cbegin(), m.cend(), [t, p](scalar_type x) {
        return (0 <= x && x < p + t);
      });
      if (it != m.cend()) {
        uint64_t r, c;
        std::tie(r, c) = m.coords(it);

        LIBSEMIGROUPS_EXCEPTION(
            "invalid entry, expected values in {{0, 1, ..., {}}}, but "
            "found {} in entry ({}, {})",
            p + t,
            detail::entry_repr(*it),
            r,
            c);
      }
    }

    //! \ingroup ntpmat_group
    //!
    //! \brief Check that an entry in an ntp matrix is valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function can be used to validate that a matrix contains values in
    //! the underlying semiring.
    //!
    //! \tparam Mat the type of the 1st argument (must satisfy
    //! \ref IsNTPMat<Mat>).
    //!
    //! \param m the matrix.
    //! \param val the entry to check.
    //!
    //! \throws LibsemigroupsException if the entry \p val is not in
    //! the set \f$\{0, 1, \ldots, t + p - 1\}\f$ where \f$t\f$ is the
    //! threshold, and \f$p\f$ the period, of the matrix.
    //!
    //! \throws LibsemigroupsException if the underlying semiring is not
    //! defined (only applies to matrices with run time arithmetic).
    template <typename Mat>
    std::enable_if_t<IsNTPMat<Mat>>
    throw_if_bad_entry(Mat const& m, typename Mat::scalar_type val) {
      detail::throw_if_semiring_nullptr(m);
      using scalar_type   = typename Mat::scalar_type;
      scalar_type const t = matrix::threshold(m);
      scalar_type const p = matrix::period(m);
      if (val < 0 || val >= p + t) {
        LIBSEMIGROUPS_EXCEPTION(
            "invalid entry, expected values in {{0, 1, ..., {}}}, but "
            "found {}",
            p + t,
            detail::entry_repr(val));
      }
    }
  }  // namespace matrix

  ////////////////////////////////////////////////////////////////////////
  // Projective max-plus matrices
  ////////////////////////////////////////////////////////////////////////

  namespace detail {
    template <typename T>
    class ProjMaxPlusMat : MatrixPolymorphicBase {
     public:
      using scalar_type            = typename T::scalar_type;
      using scalar_reference       = typename T::scalar_reference;
      using scalar_const_reference = typename T::scalar_const_reference;
      using semiring_type          = void;

      using container_type = typename T::container_type;
      using iterator       = typename T::iterator;
      using const_iterator = typename T::const_iterator;

      using underlying_matrix_type = T;

      using RowView = typename T::RowView;

      // Note that Rows are never normalised, and that's why we use the
      // underlying matrix Row type and not 1 x n ProjMaxPlusMat's instead
      // (since these will be normalised according to their entries, and
      // this might not correspond to the normalised entries of the matrix).
      using Row = typename T::Row;

      scalar_type scalar_one() const noexcept {
        return _underlying_mat.scalar_one();
      }

      scalar_type scalar_zero() const noexcept {
        return _underlying_mat.scalar_zero();
      }

      ////////////////////////////////////////////////////////////////////////
      // ProjMaxPlusMat - Constructors + destructor - public
      ////////////////////////////////////////////////////////////////////////

      ProjMaxPlusMat() : _is_normalized(false), _underlying_mat() {}
      ProjMaxPlusMat(ProjMaxPlusMat const&)            = default;
      ProjMaxPlusMat(ProjMaxPlusMat&&)                 = default;
      ProjMaxPlusMat& operator=(ProjMaxPlusMat const&) = default;
      ProjMaxPlusMat& operator=(ProjMaxPlusMat&&)      = default;

      ProjMaxPlusMat(size_t r, size_t c)
          : _is_normalized(false), _underlying_mat(r, c) {}

      // TODO(1) other missing constructors
      ProjMaxPlusMat(
          typename underlying_matrix_type::semiring_type const* semiring,
          size_t                                                r,
          size_t                                                c)
          : _is_normalized(false), _underlying_mat(semiring, r, c) {}

      explicit ProjMaxPlusMat(std::vector<std::vector<scalar_type>> const& m)
          : _is_normalized(false), _underlying_mat(m) {
        normalize();
      }

      ProjMaxPlusMat(
          std::initializer_list<std::initializer_list<scalar_type>> const& m)
          : ProjMaxPlusMat(
              std::vector<std::vector<scalar_type>>(m.begin(), m.end())) {}

      ~ProjMaxPlusMat() = default;

      ProjMaxPlusMat one() const {
        auto result = ProjMaxPlusMat(_underlying_mat.one());
        return result;
      }

      static ProjMaxPlusMat one(size_t n) {
        return ProjMaxPlusMat(T::one(n));
      }

      ////////////////////////////////////////////////////////////////////////
      // Comparison operators
      ////////////////////////////////////////////////////////////////////////

      // TODO template
      bool operator==(ProjMaxPlusMat const& that) const {
        normalize();
        that.normalize();
        return _underlying_mat == that._underlying_mat;
      }

      // TODO template
      bool operator!=(ProjMaxPlusMat const& that) const {
        return !(_underlying_mat == that._underlying_mat);
      }

      // TODO template
      bool operator<(ProjMaxPlusMat const& that) const {
        normalize();
        that.normalize();
        return _underlying_mat < that._underlying_mat;
      }

      // TODO template
      bool operator>(ProjMaxPlusMat const& that) const {
        return that < *this;
      }

      template <typename Thing>
      bool operator>=(Thing const& that) const {
        static_assert(IsMatrix<Thing> || std::is_same_v<Thing, RowView>);
        return that < *this || that == *this;
      }

      // not noexcept because operator< isn't
      template <typename Thing>
      bool operator<=(Thing const& that) const {
        static_assert(IsMatrix<Thing> || std::is_same_v<Thing, RowView>);
        return *this < that || that == *this;
      }

      ////////////////////////////////////////////////////////////////////////
      // Attributes
      ////////////////////////////////////////////////////////////////////////

      scalar_reference operator()(size_t r, size_t c) {
        // to ensure the returned value is normalised
        normalize();
        // to ensure that the matrix is renormalised if the returned scalar is
        // assigned.
        _is_normalized = false;
        return _underlying_mat(r, c);
      }

      scalar_reference at(size_t r, size_t c) {
        matrix::throw_if_bad_coords(*this, r, c);
        return this->operator()(r, c);
      }

      scalar_const_reference operator()(size_t r, size_t c) const {
        normalize();
        return _underlying_mat(r, c);
      }

      scalar_const_reference at(size_t r, size_t c) const {
        matrix::throw_if_bad_coords(*this, r, c);
        return this->operator()(r, c);
      }

      size_t number_of_rows() const noexcept {
        return _underlying_mat.number_of_rows();
      }

      size_t number_of_cols() const noexcept {
        return _underlying_mat.number_of_cols();
      }

      size_t hash_value() const {
        normalize();
        return Hash<T>()(_underlying_mat);
      }

      ////////////////////////////////////////////////////////////////////////
      // Arithmetic operators - in-place
      ////////////////////////////////////////////////////////////////////////

      void product_inplace_no_checks(ProjMaxPlusMat const& A,
                                     ProjMaxPlusMat const& B) {
        _underlying_mat.product_inplace_no_checks(A._underlying_mat,
                                                  B._underlying_mat);
        normalize(true);  // force normalize
      }

      void operator+=(ProjMaxPlusMat const& that) {
        _underlying_mat += that._underlying_mat;
        normalize(true);  // force normalize
      }

      void operator*=(scalar_type a) {
        _underlying_mat *= a;
        normalize(true);  // force normalize
      }

      void operator+=(scalar_type a) {
        _underlying_mat += a;
        normalize(true);  // force normalize
      }

      ProjMaxPlusMat operator*(scalar_type a) const {
        ProjMaxPlusMat result(*this);
        result *= a;
        return result;
      }

      ProjMaxPlusMat operator+(scalar_type a) const {
        ProjMaxPlusMat result(*this);
        result += a;
        return result;
      }

      ////////////////////////////////////////////////////////////////////////
      // Arithmetic operators - not in-place
      ////////////////////////////////////////////////////////////////////////

      ProjMaxPlusMat operator+(ProjMaxPlusMat const& that) const {
        return ProjMaxPlusMat(_underlying_mat + that._underlying_mat);
      }

      ProjMaxPlusMat operator*(ProjMaxPlusMat const& that) const {
        return ProjMaxPlusMat(_underlying_mat * that._underlying_mat);
      }

      ////////////////////////////////////////////////////////////////////////
      // Iterators
      ////////////////////////////////////////////////////////////////////////

      // The following should probably be commented out because I can't
      // currently think how to ensure that the matrix is normalised if it's
      // changed this way.

      iterator begin() noexcept {
        // to ensure the returned value is normalised
        normalize();
        // to ensure that the matrix is renormalised if the returned scalar is
        // assigned.
        _is_normalized = false;
        return _underlying_mat.begin();
      }

      iterator end() noexcept {
        // to ensure the returned value is normalised
        normalize();
        // to ensure that the matrix is renormalised if the returned scalar is
        // assigned.
        _is_normalized = false;
        return _underlying_mat.end();
      }

      const_iterator begin() const noexcept {
        normalize();
        return _underlying_mat.begin();
      }

      const_iterator end() const noexcept {
        normalize();
        return _underlying_mat.end();
      }

      const_iterator cbegin() const noexcept {
        normalize();
        return _underlying_mat.cbegin();
      }

      const_iterator cend() const noexcept {
        normalize();
        return _underlying_mat.cend();
      }

      ////////////////////////////////////////////////////////////////////////
      // Modifiers
      ////////////////////////////////////////////////////////////////////////

      void swap(ProjMaxPlusMat& that) noexcept {
        std::swap(_underlying_mat, that._underlying_mat);
      }

      void transpose() noexcept {
        _underlying_mat.transpose();
      }

      void transpose_no_checks() noexcept {
        _underlying_mat.transpose_no_checks();
      }

      ////////////////////////////////////////////////////////////////////////
      // Rows
      ////////////////////////////////////////////////////////////////////////

      RowView row(size_t i) const {
        normalize();
        return _underlying_mat.row(i);
      }

      template <typename C>
      void rows(C& x) const {
        normalize();
        return _underlying_mat.rows(x);
      }

      ////////////////////////////////////////////////////////////////////////
      // Friend functions
      ////////////////////////////////////////////////////////////////////////

      friend std::ostream& operator<<(std::ostream&         os,
                                      ProjMaxPlusMat const& x) {
        x.normalize();
        os << detail::to_string(x._underlying_mat);
        return os;
      }

      T const& underlying_matrix() const noexcept {
        normalize();
        return _underlying_mat;
      }

     private:
      explicit ProjMaxPlusMat(T&& mat)
          : _is_normalized(false), _underlying_mat(std::move(mat)) {
        normalize();
      }

      void normalize(bool force = false) const {
        if ((_is_normalized && !force)
            || (_underlying_mat.number_of_rows() == 0)
            || (_underlying_mat.number_of_cols() == 0)) {
          _is_normalized = true;
          return;
        }
        scalar_type const n = *std::max_element(_underlying_mat.cbegin(),
                                                _underlying_mat.cend());
        std::for_each(_underlying_mat.begin(),
                      _underlying_mat.end(),
                      [&n](scalar_type& s) {
                        if (s != NEGATIVE_INFINITY) {
                          s -= n;
                        }
                      });
        _is_normalized = true;
      }

      mutable bool _is_normalized;
      mutable T    _underlying_mat;
    };
  }  // namespace detail

  //! \defgroup projmaxplus_group Projective max-plus matrices
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This page describes the functionality for \f$n \times n\f$ projective
  //! matrices over the max-plus semiring. Projective max-plus matrices belong
  //! to the quotient of the monoid of all max-plus matrices by the congruence
  //! where two matrices are related if they differ by a scalar multiple; see
  //! MaxPlusMat.
  //!
  //! There are two types of such matrices those whose dimension is known at
  //! compile-time, and those where it is not.  Both types can be accessed via
  //! the alias template \ref ProjMaxPlusMat : if \c R or \c C has value
  //! \c 0, then the dimensions can be set at run time, otherwise \c R and \c
  //! C are the dimensions.  The default value of \c R is \c 0, \c C is \c R,
  //! and `Scalar` is `int`.
  //!
  //! Matrices in both these classes are modified when constructed to be in a
  //! normal form which is obtained by subtracting the maximum finite entry in
  //! the matrix from the every finite entry.
  //!
  //! A \ref ProjMaxPlusMat is neither a `StaticMatrix` nor a `DynamicMatrix`,
  //! but has the same interface as each of these types. Every instance of
  //! \ref ProjMaxPlusMat wraps a \ref MaxPlusMat.
  //!
  //! \note
  //! The types RowView and Row are the same as those
  //! in the wrapped matrix. This means that a Row object for a
  //! projective max-plus matrix is never normalised, because if they were
  //! they would be normalised according to their entries, and this might not
  //! correspond to the normalised entries of the matrix.
  //!
  //! Please refer to the documentation of these class templates for more
  //! details. The only substantial difference in the interface of
  //! static and dynamic matrices is that the former can be default
  //! constructed and the latter should be constructed using the dimensions.
  //!
  //! \par Example
  //! \code
  //! // default construct an uninitialized 3 x 3 static matrix
  //! ProjMaxPlusMat<3> m;
  //! // construct an uninitialized 4 x 4 dynamic matrix
  //! ProjMaxPlusMat<>  m(4, 4);
  //! \endcode

  //! \ingroup projmaxplus_group
  //!
  //! \brief Alias for static projective max-plus matrices with compile-time
  //! arithmetic and dimensions.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for static projective max-plus matrices whose arithmetic and
  //! dimensions are defined at compile-time.
  //!
  //! \tparam R the number of rows.
  //! \tparam C the number of columns.
  //! \tparam Scalar the type of the entries in the matrix.
  template <size_t R, size_t C, typename Scalar>
  using StaticProjMaxPlusMat
      = detail::ProjMaxPlusMat<StaticMaxPlusMat<R, C, Scalar>>;

  //! \ingroup projmaxplus_group
  //!
  //! \brief Alias for dynamic projective max-plus matrices with run-time
  //! dimensions.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias for the type of dynamic projective max-plus matrices where the
  //! dimensions of the matrices can be defined at run time.
  //!
  //! \tparam Scalar the type of the entries in the matrix.
  template <typename Scalar>
  using DynamicProjMaxPlusMat
      = detail::ProjMaxPlusMat<DynamicMaxPlusMat<Scalar>>;

  //! \ingroup projmaxplus_group
  //! \brief Alias template for projective max-plus matrices.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Alias template for projective max-plus matrices.
  //!
  //! \tparam R the number of rows.  A value of \c 0 indicates that the value
  //! will be set at run time (default: \c 0).
  //!
  //! \tparam C the number of columns.  A value of \c 0 indicates that the
  //! value will be set at run time (default: \c R).
  //!
  //! \tparam Scalar The type of the entries in the matrix (default: `int`).
  template <size_t R = 0, size_t C = R, typename Scalar = int>
  using ProjMaxPlusMat = std::conditional_t<R == 0 || C == 0,
                                            DynamicProjMaxPlusMat<Scalar>,
                                            StaticProjMaxPlusMat<R, C, Scalar>>;

  namespace detail {
    template <typename T>
    struct IsProjMaxPlusMatHelper : std::false_type {};

    template <size_t R, size_t C, typename Scalar>
    struct IsProjMaxPlusMatHelper<StaticProjMaxPlusMat<R, C, Scalar>>
        : std::true_type {};

    template <typename Scalar>
    struct IsProjMaxPlusMatHelper<DynamicProjMaxPlusMat<Scalar>>
        : std::true_type {};
  }  // namespace detail

  //! \ingroup projmaxplus_group
  //!
  //! \brief Helper to check if a type is \ref ProjMaxPlusMat.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This variable has value \c true if the template parameter \c T is the
  //! same as \ref ProjMaxPlusMat<R, C, Scalar> for some values of \c R, \c C,
  //! and \c Scalar; and \c false if it is not.
  //!
  //! \tparam T the type to check.
  template <typename T>
  static constexpr bool IsProjMaxPlusMat
      = detail::IsProjMaxPlusMatHelper<T>::value;

  namespace matrix {
    // \ingroup projmaxplus_group
    //
    //! \brief Check that the entries in a projective max-plus matrix are
    //! valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function can be used to validate that the matrix \p x contains
    //! values in the underlying semiring, by checking the underlying matrix.
    //!
    //! \tparam Mat the type of the parameter (must satisfy \ref
    //! IsProjMaxPlusMat<Mat>).
    //!
    //! \param x the matrix to check.
    //!
    //! \throws LibsemigroupsException if
    //! `throw_if_bad_entry(x.underlying_matrix())` throws.
    template <typename Mat>
    constexpr std::enable_if_t<IsProjMaxPlusMat<Mat>>
    throw_if_bad_entry(Mat const& x) {
      throw_if_bad_entry(x.underlying_matrix());
    }

    // \ingroup projmaxplus_group
    //
    //! \brief Check that an entry in a projective max-plus matrix is valid.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function can be used to validate that the matrix \p x contains
    //! values in the underlying semiring, using the underlying matrix.
    //!
    //! \tparam Mat the type of the parameter (must satisfy \ref
    //! IsProjMaxPlusMat<Mat>).
    //!
    //! \param x the matrix.
    //! \param val the entry.
    //!
    //! \throws LibsemigroupsException if
    //! `throw_if_bad_entry(x.underlying_matrix(), val)` throws.
    template <typename Mat>
    constexpr std::enable_if_t<IsProjMaxPlusMat<Mat>>
    throw_if_bad_entry(Mat const& x, typename Mat::scalar_type val) {
      throw_if_bad_entry(x.underlying_matrix(), val);
    }

    ////////////////////////////////////////////////////////////////////////
    // Matrix helpers - pow
    ////////////////////////////////////////////////////////////////////////

    //! \brief Returns a power of a matrix.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function returns the matrix \p x to the power \p e.
    //! If \p e is \c 0, then the identity matrix is returned; if
    //! \p e is \c 1, then a copy of the parameter \p x is returned.
    //!
    //! \tparam Mat the type of the matrix \p x, must satisfy \ref
    //! IsMatrix<Mat>.
    //!
    //! \param x  the matrix  (must have equal number of rows and columns).
    //!
    //! \param e  the exponent (must be \f$\geq 0\f$).
    //!
    //! \returns The matrix \p x to the power \p e.
    //!
    //! \throws LibsemigroupsException if `Mat::scalar_type` is a signed type
    //! and the parameter \p e is less than \c 0.
    //!
    //! \throws LibsemigroupsException if \p x is not a square matrix.
    //!
    //! \complexity
    //! \f$O(m\log_2(e))\f$ where \f$m\f$ is the dimension of the matrix
    //! \p x and \f$e\f$ is the parameter \p e.
    //!
    //! \par Example
    //!
    //! \code
    //! auto x == ProjMaxPlusMat<>::make({{-2, 2, 0}, {-1, 0, 0}, {1, -3,
    //! 1}}));
    //! // returns {{-1, 0, -1}, {-2, -1, -2}, {-1, 0, -1}}
    //! matrix::pow(x, 100);
    //! \endcode
    // TODO(1) pow_no_checks
    // TODO(2) version that changes x in-place
    // TODO(0) update the code in the examples
    template <typename Mat>
    Mat pow(Mat const& x, typename Mat::scalar_type e) {
      using scalar_type = typename Mat::scalar_type;

      if constexpr (std::is_signed<scalar_type>::value) {
        if (e < 0) {
          LIBSEMIGROUPS_EXCEPTION(
              "negative exponent, expected value >= 0, found {}", e);
        }
      }

      throw_if_not_square(x);

      typename Mat::semiring_type const* sr = nullptr;

      if constexpr (IsMatWithSemiring<Mat>) {
        sr = x.semiring();
      }

      if (e == 0) {
        return x.one();
      }

      auto y = Mat(x);
      if (e == 1) {
        return y;
      }
      auto z = (e % 2 == 0 ? x.one() : y);

      Mat tmp(sr, x.number_of_rows(), x.number_of_cols());
      while (e > 1) {
        tmp.product_inplace_no_checks(y, y);
        std::swap(y, tmp);
        e /= 2;
        if (e % 2 == 1) {
          tmp.product_inplace_no_checks(z, y);
          std::swap(z, tmp);
        }
      }
      return z;
    }

    ////////////////////////////////////////////////////////////////////////
    // Matrix helpers - rows
    ////////////////////////////////////////////////////////////////////////

    //! \brief Returns a std::vector of row views into the rows of a dynamic
    //! matrix.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! Returns a std::vector of row views into the rows of the dynamic
    //! matrix \p x.
    //!
    //! \tparam Mat the type of \p x, must satisfy \ref IsDynamicMatrix<Mat>.
    //!
    //! \param x the matrix.
    //!
    //! \returns A std::vector of `Mat::RowView` of size `x.number_of_rows()`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(m)\f$ where \f$m\f$ is the number of rows in the matrix \p x.
    template <typename Mat, typename = std::enable_if_t<IsDynamicMatrix<Mat>>>
    std::vector<typename Mat::RowView> rows(Mat const& x) {
      std::vector<typename Mat::RowView> container;
      x.rows(container);
      return container;
    }

    //! \brief Returns a static vector of row views into the rows of a static
    //! matrix.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! Returns a static vector of row views into the rows of the static
    //! matrix \p x.
    //!
    //! \tparam Mat the type of \p x, must satisfy \ref IsStaticMatrix<Mat>.
    //!
    //! \param x the matrix.
    //!
    //! \returns A static vector of `Mat::RowView` of size
    //! `x.number_of_rows()`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(m)\f$ where \f$m\f$ is the number of rows in the matrix \p x.
    template <typename Mat, typename = std::enable_if_t<IsStaticMatrix<Mat>>>
    detail::StaticVector1<typename Mat::RowView, Mat::nr_rows>
    rows(Mat const& x) {
      detail::StaticVector1<typename Mat::RowView, Mat::nr_rows> container;
      x.rows(container);
      return container;
    }

    ////////////////////////////////////////////////////////////////////////
    // Matrix helpers - bitset_rows
    ////////////////////////////////////////////////////////////////////////

    // The main function
    //! \brief Converts a container of row views of a boolean matrix to bit
    //! sets, and append them to another container.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function converts the row views (of a boolean matrix) in a
    //! container to bit sets and appends them to \p result.
    //!
    //! \tparam Mat the type of matrix. This must be a type so that
    //! IsBMat<Mat> is \c true.
    //!
    //! \tparam R an upper bound for the number of rows in `views`. This value
    //! must be at most `BitSet<1>::max_size()`.
    //!
    //! \tparam C an upper bound for the number of columns in each row
    //! represented in `views`. This value must be at most
    //! `BitSet<1>::max_size()`.
    //!
    //! \tparam Container the type of the container `views`, should be
    //! `std::vector` or `detail::StaticVector1`.
    //!
    //! \param views   a container of `Mat::RowView` or `std::vector<bool>`.
    //! \param result  a static vector of bit sets to contain the resulting
    //! bit sets.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is the number of rows in `views` and
    //! and \f$n\f$ is the number of columns in any vector in `views`.
    template <typename Mat, size_t R, size_t C, typename Container>
    void bitset_rows(Container&&                          views,
                     detail::StaticVector1<BitSet<C>, R>& result) {
      using RowView    = typename Mat::RowView;
      using value_type = typename std::decay_t<Container>::value_type;
      // std::vector<bool> is used as value_type in the benchmarks
      static_assert(IsBMat<Mat>, "IsBMat<Mat> must be true!");
      static_assert(std::is_same<value_type, RowView>::value
                        || std::is_same<value_type, std::vector<bool>>::value,
                    "Container::value_type must equal Mat::RowView or "
                    "std::vector<bool>!!");
      static_assert(R <= BitSet<1>::max_size(),
                    "R must be at most BitSet<1>::max_size()!");
      static_assert(C <= BitSet<1>::max_size(),
                    "C must be at most BitSet<1>::max_size()!");
      LIBSEMIGROUPS_ASSERT(views.size() <= R);
      LIBSEMIGROUPS_ASSERT(views.empty() || views[0].size() <= C);
      for (auto const& v : views) {
        result.emplace_back(v.cbegin(), v.cend());
      }
    }

    //! \brief Converts a container of row views of a boolean matrix to bit
    //! sets, and return them.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function converts the container of row views \p views of a
    //! boolean matrix to bit sets, and return them.
    //!
    //! \tparam Mat the type of matrix. This must be a type so that
    //! IsBMat<Mat> is \c true.
    //!
    //! \tparam R an upper bound for the number of rows in \p views. This
    //! value must be at most BitSet<1>::max_size().
    //!
    //! \tparam C an upper bound for the number of columns in each row
    //! represented in \p views. This value must be at most
    //! BitSet<1>::max_size().
    //!
    //! \tparam Container the type of the container \p views, should be
    //! `std::vector` or `detail::StaticVector1`.
    //!
    //! \param views  a container of `Mat::RowView` or `std::vector<bool>`.
    //!
    //! \returns  A value of type `detail::StaticVector1<BitSet<C>,
    //! R>`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is the number of rows in \p views and
    //! and \f$n\f$ is the number of columns in any vector in \p views.
    template <typename Mat, size_t R, size_t C, typename Container>
    auto bitset_rows(Container&& views) {
      using RowView    = typename Mat::RowView;
      using value_type = typename std::decay_t<Container>::value_type;
      static_assert(IsBMat<Mat>, "IsBMat<Mat> must be true!");
      static_assert(std::is_same<value_type, RowView>::value
                        || std::is_same<value_type, std::vector<bool>>::value,
                    "Container::value_type must equal Mat::RowView or "
                    "std::vector<bool>!!");
      static_assert(R <= BitSet<1>::max_size(),
                    "R must be at most BitSet<1>::max_size()!");
      static_assert(C <= BitSet<1>::max_size(),
                    "C must be at most BitSet<1>::max_size()!");
      LIBSEMIGROUPS_ASSERT(views.size() <= R);
      LIBSEMIGROUPS_ASSERT(views.empty() || views[0].size() <= C);
      detail::StaticVector1<BitSet<C>, R> result;
      bitset_rows<Mat>(std::forward<Container>(views), result);
      return result;
    }

    // Helper
    //! \brief Computes the rows of a boolean matrix as bit sets and appends
    //! them to a container.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! Computes the rows of the matrix \p x as bit sets and
    //! appends them to \p result.
    //!
    //! \tparam Mat the type of matrix. This must be a type so that
    //! IsBMat<Mat> is \c true.
    //!
    //! \tparam R an upper bound for the number of rows in `views`. This value
    //! must be at most BitSet<1>::max_size().
    //!
    //! \tparam C an upper bound for the number of columns in each row
    //! represented in `views`. This value must be at most
    //! BitSet<1>::max_size().
    //!
    //! \param x the boolean matrix.
    //!
    //! \param result the container to append the rows of \p x to.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is the number of rows in \p x and and
    //! \f$n\f$ is the number of columns in \p x.
    template <typename Mat, size_t R, size_t C>
    void bitset_rows(Mat const&                           x,
                     detail::StaticVector1<BitSet<C>, R>& result) {
      static_assert(IsBMat<Mat>, "IsBMat<Mat> must be true!");
      static_assert(R <= BitSet<1>::max_size(),
                    "R must be at most BitSet<1>::max_size()!");
      static_assert(C <= BitSet<1>::max_size(),
                    "C must be at most BitSet<1>::max_size()!");
      LIBSEMIGROUPS_ASSERT(x.number_of_cols() <= C);
      LIBSEMIGROUPS_ASSERT(x.number_of_rows() <= R);
      bitset_rows<Mat>(std::move(rows(x)), result);
    }

    // Helper
    //! \brief Computes the rows of a boolean matrix as bit sets.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This functions returns the rows of a boolean matrix as bit sets.
    //!
    //! \tparam Mat the type of matrix. Must satisfy \ref IsBMat<Mat>.
    //!
    //! \param x  the boolean matrix.
    //!
    //! \returns  A value of type `detail::StaticVector1<BitSet<C>, R>`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is the number of rows in \p x and
    //! and \f$n\f$ is the number of columns in \p x.
    template <typename Mat>
    auto bitset_rows(Mat const& x) {
      static_assert(IsBMat<Mat>, "IsBMat<Mat> must be true!");
      LIBSEMIGROUPS_ASSERT(x.number_of_rows() <= BitSet<1>::max_size());
      LIBSEMIGROUPS_ASSERT(x.number_of_cols() <= BitSet<1>::max_size());
      size_t const M = detail::BitSetCapacity<Mat>::value;
      return bitset_rows<Mat, M, M>(std::move(rows(x)));
    }

    ////////////////////////////////////////////////////////////////////////
    // Matrix helpers - bitset_row_basis
    ////////////////////////////////////////////////////////////////////////

    //! \brief Appends a basis for the space spanned by some bitsets to a
    //! container.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! Appends a basis for the space spanned by the bitsets in \p rows to
    //! the container \p result.
    //!
    //! \tparam Mat a type such that \ref IsBMat<Mat> is \c true.
    //!
    //! \tparam Container a container type with `Container::value_type` being
    //! `BitSet<M>` or `std::bitset<M>` for some \c M.
    //!
    //! \param rows   container of spanning rows represented by bit sets.
    //! \param result  container for the resulting row basis.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(r ^ 2 c)\f$ where \f$r\f$ is the size of `rows` and
    //! \f$c\f$ is the size of each bitset in `rows`.
    // This works with std::vector and StaticVector1, with value_type equal
    // to std::bitset and BitSet.
    template <typename Mat, typename Container>
    void bitset_row_basis(Container&& rows, std::decay_t<Container>& result) {
      using value_type = typename std::decay_t<Container>::value_type;
      static_assert(IsBMat<Mat>, "IsBMat<Mat> must be true!");
      static_assert(IsBitSet<value_type> || detail::IsStdBitSet<value_type>,
                    "Container::value_type must be BitSet or std::bitset");
      LIBSEMIGROUPS_ASSERT(rows.size() <= BitSet<1>::max_size());
      LIBSEMIGROUPS_ASSERT(rows.empty()
                           || rows[0].size() <= BitSet<1>::max_size());

      std::sort(rows.begin(), rows.end(), detail::LessBitSet());
      // Remove duplicates
      rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
      for (size_t i = 0; i < rows.size(); ++i) {
        value_type cup;
        cup.reset();
        for (size_t j = 0; j < i; ++j) {
          if ((rows[i] & rows[j]) == rows[j]) {
            cup |= rows[j];
          }
        }
        for (size_t j = i + 1; j < rows.size(); ++j) {
          if ((rows[i] & rows[j]) == rows[j]) {
            cup |= rows[j];
          }
        }
        if (cup != rows[i]) {
          result.push_back(std::move(rows[i]));
        }
      }
    }

    //! \brief Returns a basis for the space spanned by some bit sets.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! Returns a basis for the space spanned by the bit sets in \p rows.
    //!
    //! \tparam Mat a type such that \ref IsBMat<Mat> is \c true.
    //!
    //! \tparam Container a container type such that `Container::value_type`
    //! is `BitSet<M>` or `std::bitset<M>` for some \c M.
    //!
    //! \param rows container of spanning rows represented by bit sets.
    //!
    //! \returns
    //! A container of type `std::decay_t<Container>` containing the row
    //! basis consisting of bit sets.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(r ^ 2 c)\f$ where \f$r\f$ is the size of \p rows and
    //! \f$c\f$ is the size of each bitset in \p rows.
    template <typename Mat, typename Container>
    std::decay_t<Container> bitset_row_basis(Container&& rows) {
      using value_type = typename std::decay_t<Container>::value_type;
      static_assert(IsBMat<Mat>, "IsBMat<Mat> must be true!");
      static_assert(IsBitSet<value_type> || detail::IsStdBitSet<value_type>,
                    "Container::value_type must be BitSet or std::bitset");
      LIBSEMIGROUPS_ASSERT(rows.size() <= BitSet<1>::max_size());
      LIBSEMIGROUPS_ASSERT(rows.empty()
                           || rows[0].size() <= BitSet<1>::max_size());
      std::decay_t<Container> result;
      bitset_row_basis<Mat>(std::forward<Container>(rows), result);
      return result;
    }

    //! \brief Returns a basis for the space spanned by the rows of the
    //! boolean matrix.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function returns a basis for the space spanned by the rows of the
    //! boolean matrix \p x.
    //!
    //! \tparam Mat a type such that IsBMat<Mat> is \c true.
    //!
    //! \tparam M an upper bound for the dimensions of the returned container.
    //! If \ref IsStaticMatrix<Mat> is \c true, then \p M is the number of
    //! rows (or columns) in the square matrix \p x. Otherwise, if
    //! IsDynamicMatrix<Mat> is \c true, then \p M is BitSet<1>::max_size().
    //!
    //! \param x  the boolean matrix.
    //!
    //! \returns
    //! A container of type `detail::StaticVector1<BitSet<M>, M>>`
    //! containing the row basis of \p x consisting of bitsets.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(r ^ 2 c)\f$ where \f$r\f$ is the number of rows in \p x and
    //! \f$c\f$ is the number of columns in \p x.
    template <typename Mat, size_t M = detail::BitSetCapacity<Mat>::value>
    detail::StaticVector1<BitSet<M>, M> bitset_row_basis(Mat const& x) {
      static_assert(IsBMat<Mat>, "IsBMat<Mat> must be true!");
      LIBSEMIGROUPS_ASSERT(x.number_of_rows() <= BitSet<1>::max_size());
      LIBSEMIGROUPS_ASSERT(x.number_of_cols() <= BitSet<1>::max_size());
      detail::StaticVector1<BitSet<M>, M> result;
      bitset_row_basis<Mat>(std::move(bitset_rows(x)), result);
      return result;
    }

    //! \brief Appends a basis for the rowspace of a boolean matrix to a
    //! container.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function appends a basis for the rowspace of the boolean matrix
    //! \p x to the container \p result.
    //!
    //! \tparam Mat a type such that \ref IsBMat<Mat> is \c true.
    //!
    //! \tparam Container a container type with `Container::value_type` equal
    //! to `BitSet<M>` or `std::bitset<M>` for some \p M.
    //!
    //! \param x  the boolean matrix.
    //!
    //! \param result  container for the resulting rowbasis.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(r ^ 2 c)\f$ where \f$r\f$ is the number of rows in \p x
    //! and \f$c\f$ is the number of columns in \p x.
    template <typename Mat, typename Container>
    void bitset_row_basis(Mat const& x, Container& result) {
      using value_type = typename Container::value_type;
      static_assert(IsBMat<Mat>, "IsBMat<Mat> must be true!");
      static_assert(IsBitSet<value_type> || detail::IsStdBitSet<value_type>,
                    "Container::value_type must be BitSet or std::bitset");
      LIBSEMIGROUPS_ASSERT(x.number_of_rows() <= BitSet<1>::max_size());
      LIBSEMIGROUPS_ASSERT(x.number_of_cols() <= BitSet<1>::max_size());
      bitset_row_basis<Mat>(std::move(bitset_rows(x)), result);
    }

    ////////////////////////////////////////////////////////////////////////
    // Matrix helpers - row_basis - MaxPlusTruncMat
    ////////////////////////////////////////////////////////////////////////

    //! \brief Appends a basis for a space spanned by row views or bit sets to
    //! a container.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function appends a basis for the space spanned by the row
    //! views or bit sets in \p views to the container \p result.
    //!
    //! \tparam Mat a type satisfying \ref IsMaxPlusTruncMat<Mat> or
    //! \ref IsBMat<Mat>.
    //!
    //! \tparam Container a container type (such as `std::vector`, for
    //! example). The `Container::value_type` must be `Mat::RowView`, or if
    //! `Mat` is \ref BMat<N> for some \c N, then  `Container::value_type` can
    //! additionally be `BitSet<M>` or `std::bitset<M>` where \c M is
    //! greater than or equal to \c N.
    //!
    //! \param views   container of spanning row views or bit sets.
    //! \param result  container for the resulting row basis.
    //!
    //! \exceptions
    //! This function guarantees not to throw a
    //! LibsemigroupsException.
    //!
    //! \complexity
    //! \f$O(r ^ 2 c)\f$ where \f$r\f$ is the size of \p views and
    //! \f$c\f$ is the size of each row view or bit set in \p views.
    template <typename Mat, typename Container>
    std::enable_if_t<IsMaxPlusTruncMat<Mat>>
    row_basis(Container&& views, std::decay_t<Container>& result) {
      using value_type = typename std::decay_t<Container>::value_type;
      static_assert(std::is_same<value_type, typename Mat::RowView>::value,
                    "Container::value_type must be Mat::RowView");
      using scalar_type = typename Mat::scalar_type;
      using Row         = typename Mat::Row;

      if (views.empty()) {
        return;
      }

      LIBSEMIGROUPS_ASSERT(result.empty());

      std::sort(views.begin(), views.end());
      Row tmp1(views[0]);

      for (size_t r1 = 0; r1 < views.size(); ++r1) {
        if (r1 == 0 || views[r1] != views[r1 - 1]) {
          std::fill(tmp1.begin(), tmp1.end(), tmp1.scalar_zero());
          for (size_t r2 = 0; r2 < r1; ++r2) {
            scalar_type max_scalar = matrix::threshold(tmp1);
            for (size_t c = 0; c < tmp1.number_of_cols(); ++c) {
              if (views[r2][c] == tmp1.scalar_zero()) {
                continue;
              }
              if (views[r1][c] >= views[r2][c]) {
                if (views[r1][c] != matrix::threshold(tmp1)) {
                  max_scalar
                      = std::min(max_scalar, views[r1][c] - views[r2][c]);
                }
              } else {
                max_scalar = tmp1.scalar_zero();
                break;
              }
            }
            if (max_scalar != tmp1.scalar_zero()) {
              tmp1 += views[r2] * max_scalar;
            }
          }
          if (tmp1 != views[r1]) {
            result.push_back(views[r1]);
          }
        }
      }
    }

    ////////////////////////////////////////////////////////////////////////
    // Matrix helpers - row_basis - BMat
    ////////////////////////////////////////////////////////////////////////

    // This version of row_basis for BMat's is for used for compatibility
    // with the MatrixCommon framework, i.e. so that BMat's exhibit the same
    // interface/behaviour as other matrices.
    //
    // This version takes a container of row views of BMat, converts it to a
    // container of BitSets, computes the row basis using the BitSets, then
    // selects those row views in views that belong to the computed basis.

    //! \brief Appends a row space basis (consisting of row views) for the row
    //! space spanned by some rows.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function appends a basis (consisting of row views) for the row
    //! space spanned by the rows in \p views to the container \p result.
    //!
    //! \tparam Mat a type satisfying \ref IsBMat<Mat>.
    //!
    //! \tparam Container the type of a container for the rows.
    //!
    //! \param views  the spanning set of rows.
    //! \param result  the container to append the basis to the row space to.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    // TODO(2) complexity
    template <typename Mat, typename Container>
    std::enable_if_t<IsBMat<Mat>> row_basis(Container&&              views,
                                            std::decay_t<Container>& result) {
      using RowView    = typename Mat::RowView;
      using value_type = typename std::decay_t<Container>::value_type;
      // std::vector<bool> is used as value_type in the benchmarks
      static_assert(std::is_same<value_type, RowView>::value
                        || std::is_same<value_type, std::vector<bool>>::value,
                    "Container::value_type must equal Mat::RowView or "
                    "std::vector<bool>!!");

      if (views.empty()) {
        return;
      }

      // Convert RowViews to BitSets
      size_t const M  = detail::BitSetCapacity<Mat>::value;
      auto         br = bitset_rows<Mat, M, M>(std::forward<Container>(views));
      using bitset_type = typename decltype(br)::value_type;

      // Map for converting bitsets back to row views
      std::unordered_map<bitset_type, size_t> lookup;
      LIBSEMIGROUPS_ASSERT(br.size() == views.size());
      for (size_t i = 0; i < br.size(); ++i) {
        lookup.insert({br[i], i});
      }

      // Compute rowbasis using bitsets + convert back to rowviews
      for (auto const& bs : bitset_row_basis<Mat>(br)) {
        auto it = lookup.find(bs);
        LIBSEMIGROUPS_ASSERT(it != lookup.end());
        result.push_back(views[it->second]);
      }
    }

    ////////////////////////////////////////////////////////////////////////
    // Matrix helpers - row_basis - generic helpers
    ////////////////////////////////////////////////////////////////////////

    //! \brief Appends a basis (consisting of row views) for the row space of
    //! a matrix to a container.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function appends a basis (consisting of row views) for the row
    //! space of the matrix \p x to the container \p result.
    //!
    //! \tparam Mat a type satisfying \ref IsMatrix<Mat>.
    //!
    //! \tparam Container a container type (such as std::vector).
    //!
    //! \param x  the matrix.
    //! \param result  the container to append the basis to the row space to.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(r ^ 2 c)\f$ where \f$r\f$ is the number of rows in \p x
    //! and \f$c\f$ is the number of columns in \p x.
    // Row basis of rowspace of matrix <x> appended to <result>
    template <typename Mat,
              typename Container,
              typename = std::enable_if_t<IsMatrix<Mat>>>
    void row_basis(Mat const& x, Container& result) {
      row_basis<Mat>(std::move(rows(x)), result);
    }

    //! \brief Returns a row space basis of a dynamic matrix as a std::vector
    //! of row views.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function returns a row space basis of the dynamic matrix \p x as
    //! a std::vector of row views.
    //!
    //! \tparam Mat a type satisfying IsDynamicMatrix<Mat>.
    //!
    //! \param x  the matrix.
    //! \returns The row space basis of \p x.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(r ^ 2 c)\f$ where \f$r\f$ is the number of rows in \p x
    //! and \f$c\f$ is the number of columns in \p x.
    // Row basis of rowspace of matrix <x>
    template <typename Mat, typename = std::enable_if_t<IsDynamicMatrix<Mat>>>
    std::vector<typename Mat::RowView> row_basis(Mat const& x) {
      std::vector<typename Mat::RowView> container;
      row_basis(x, container);
      return container;
    }

    //! \brief Returns a row space basis of a static matrix as a static vector
    //! of row views.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function returns a row space basis of the static matrix \p x as a
    //! static vector of row views.
    //!
    //! \tparam Mat a type satisfying \ref IsStaticMatrix<Mat>.
    //!
    //! \param x  the matrix.
    //! \returns A static vector of `Mat::RowView`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(r ^ 2 c)\f$ where \f$r\f$ is the number of rows in \p x
    //! and \f$c\f$ is the number of columns in \p x.
    template <typename Mat, typename = std::enable_if_t<IsStaticMatrix<Mat>>>
    detail::StaticVector1<typename Mat::RowView, Mat::nr_rows>
    row_basis(Mat const& x) {
      detail::StaticVector1<typename Mat::RowView, Mat::nr_rows> container;
      row_basis(x, container);
      return container;
    }

    //! \brief Returns a row basis for the space spanned by a container of
    //! row views.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function returns a row space basis for the space spanned by the
    //! container \p rows of row views.
    //!
    //! \tparam Mat a type satisfying \ref IsMatrix<Mat>.
    //!
    //! \tparam Container the type of the container \p rows (must have
    //! `value_type` equal to `Mat::RowView`).
    //!
    //! \param rows  a container whose `value_type` is `Mat::RowView`.
    //! \returns A container of type \p Container of `Mat::RowView`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    // TODO(2) complexity
    template <typename Mat, typename Container>
    std::decay_t<Container> row_basis(Container&& rows) {
      using value_type = typename std::decay_t<Container>::value_type;
      static_assert(IsMatrix<Mat>, "IsMatrix<Mat> must be true!");
      static_assert(std::is_same<value_type, typename Mat::RowView>::value,
                    "Container::value_type must be Mat::RowView");

      std::decay_t<Container> result;
      row_basis<Mat>(std::forward<Container>(rows), result);
      return result;
    }

    ////////////////////////////////////////////////////////////////////////
    // Matrix helpers - row_space_size
    ////////////////////////////////////////////////////////////////////////

    //! \brief Returns the size of the row space of a boolean matrix.
    //!
    //! Defined in `matrix.hpp`.
    //!
    //! This function returns the size of the row space of the boolean matrix
    //! \p x. This is currently only implemented for types of matrix `Mat`
    //! satisfying \ref IsBMat<Mat>.
    //!
    //! \tparam Mat the type of \p x, must satisfy \ref IsBMat<Mat>.
    //! \param x  the matrix.
    //!
    //! \returns  The size of the row space of the matrix \p x.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \f$m\f$ is the size of the row basis of \p x
    //! and \f$n\f$ is the size of the row space.
    //!
    //! \warning
    //! If \f$k\f$ is the number of rows of the matrix \p x, then
    //! \f$n\f$ can be as large as \f$2 ^ k\f$.
    //!
    //! \par Example
    //!
    //! \code
    //! auto x = BMat<>::make({{1, 0, 0}, {0, 0, 1}, {0, 1, 0}});
    //! matrix::row_space_size(x); // returns 7
    //! \endcode
    // TODO(0) Update the code in the examples
    template <typename Mat, typename = std::enable_if_t<IsBMat<Mat>>>
    size_t row_space_size(Mat const& x) {
      size_t const M                 = detail::BitSetCapacity<Mat>::value;
      auto         bitset_row_basis_ = bitset_row_basis<Mat>(
          std::move(bitset_rows<Mat, M, M>(std::move(rows(x)))));

      std::unordered_set<BitSet<M>> st;
      st.insert(bitset_row_basis_.cbegin(), bitset_row_basis_.cend());
      std::vector<BitSet<M>> orb(bitset_row_basis_.cbegin(),
                                 bitset_row_basis_.cend());
      for (size_t i = 0; i < orb.size(); ++i) {
        for (auto& row : bitset_row_basis_) {
          auto cup = orb[i];
          for (size_t j = 0; j < x.number_of_rows(); ++j) {
            cup.set(j, cup[j] || row[j]);
          }
          if (st.insert(cup).second) {
            orb.push_back(std::move(cup));
          }
        }
      }
      return orb.size();
    }

  }  // namespace matrix

  //! \ingroup matrix_group
  //!
  //! \brief Add a scalar to a matrix.
  //!
  //! This function returns the sum of the matrix \p x and the scalar
  //! \p a.
  //!
  //! \tparam Mat the type of the 2nd argument, must satisfy \ref
  //! IsMatrix<Mat>.
  //!
  //! \param a the scalar to add.
  //! \param x the matrix.
  //!
  //! \returns  The sum of \p a and \p x.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  //!
  //! \warning This function does not detect overflows of `Mat::scalar_type`.
  template <typename Mat>
  auto operator+(typename Mat::scalar_type a, Mat const& x)
      -> std::enable_if_t<IsMatrix<Mat>, Mat> {
    return x + a;
  }

  //! \ingroup matrix_group
  //!
  //! \brief Multiply a matrix by a scalar.
  //!
  //! This function returns the product of the matrix \p x and the scalar
  //! \p a.
  //!
  //! \tparam Mat the type of the 2nd argument, must satisfy \ref
  //! IsMatrix<Mat>.
  //!
  //! \param a the scalar to multiply.
  //! \param x the matrix.
  //!
  //! \returns The product of \p a and \p x.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  //!
  //! \warning This function does not detect overflows of `Mat::scalar_type`.
  template <typename Mat>
  auto operator*(typename Mat::scalar_type a, Mat const& x)
      -> std::enable_if_t<IsMatrix<Mat>, Mat> {
    return x * a;
  }

  //! \ingroup matrix_group
  //!
  //! \brief Validates the arguments, constructs a matrix, and validates it.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Validates the arguments, constructs a matrix, and validates it.
  //!
  //! \tparam Mat the type of matrix being constructed, must satisfy
  //! \ref IsMatrix<Mat> and not \ref IsMatWithSemiring<Mat>.
  //!
  //! \param rows the values to be copied into the matrix.
  //! \returns The constructed matrix if valid.
  //!
  //! \throws
  //! LibsemigroupsException if `rows` does not represent a
  //! matrix of the correct dimensions.
  //!
  //! \throws
  //! LibsemigroupsException if the constructed matrix
  //! contains values that do not belong to the underlying semiring.
  //!
  //! \complexity
  //! \f$O(mn)\f$ where \f$m\f$ is the number of rows and
  //! \f$n\f$ is the number of columns of the matrix.
  template <typename Mat,
            typename
            = std::enable_if_t<IsMatrix<Mat> && !IsMatWithSemiring<Mat>>>
  Mat make(std::vector<std::vector<typename Mat::scalar_type>> const& rows) {
    detail::throw_if_any_row_wrong_size(rows);
    Mat m(rows);
    matrix::throw_if_bad_entry(m);
    return m;
  }

  //! \ingroup matrix_group
  //!
  //! \brief Validates the arguments, constructs a matrix, and validates it.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Validates the arguments, constructs a matrix, and validates it.
  //!
  //! \tparam Mat the type of matrix being constructed, must satisfy
  //! \ref IsMatrix<Mat> and not \ref IsMatWithSemiring<Mat>.
  //!
  //! \param rows the values to be copied into the matrix.
  //! \returns The constructed matrix if valid.
  //!
  //! \throws
  //! LibsemigroupsException if `rows` does not represent a
  //! matrix of the correct dimensions.
  //!
  //! \throws
  //! LibsemigroupsException if the constructed matrix
  //! contains values that do not belong to the underlying semiring.
  //!
  //! \complexity
  //! \f$O(mn)\f$ where \f$m\f$ is the number of rows and
  //! \f$n\f$ is the number of columns of the matrix.
  template <typename Mat,
            typename
            = std::enable_if_t<IsMatrix<Mat> && !IsMatWithSemiring<Mat>>>
  Mat make(std::initializer_list<std::vector<typename Mat::scalar_type>> const&
               rows) {
    return make<Mat>(std::vector<std::vector<typename Mat::scalar_type>>(rows));
  }

  //! \ingroup matrix_group
  //!
  //! \brief Constructs a row and validates it.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This function constructs a row from a std::initializer_list and then
  //! calls \ref validate.
  //!
  //! \tparam Mat the type of matrix being constructed, must satisfy
  //! \ref IsMatrix<Mat> and not \ref IsMatWithSemiring<Mat>.
  //!
  //! \param row the values to be copied into the row.
  //!
  //! \returns the constructed row if valid.
  //!
  //! \throws
  //! LibsemigroupsException if the constructed row contains
  //! values that do not belong to the underlying semiring.
  //!
  //! \complexity
  //! \f$O(n)\f$ where \f$n\f$ is the number of columns of the matrix.
  //!
  //! \warning
  //! This function only works for rows, i.e. when the template
  //! parameter \c R is \c 1.
  template <typename Mat,
            typename
            = std::enable_if_t<IsMatrix<Mat> && !IsMatWithSemiring<Mat>>>
  Mat make(std::initializer_list<typename Mat::scalar_type> const& row) {
    Mat m(row);
    matrix::throw_if_bad_entry(m);
    return m;
  }
  // TODO(1) vector version of above

  //! \ingroup matrix_group
  //!
  //! \brief Constructs a matrix (from std::initializer_list) and validates
  //! it.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Validates the arguments, constructs a matrix and validates it.
  //!
  //! \tparam Mat the type of the matrix being constructed (must satisfy
  //! \ref IsMatrix<Mat>).
  //!
  //! \tparam Semiring the type of the semiring where arithmetic is performed.
  //!
  //! \param semiring  a pointer to const semiring object.
  //!
  //! \param rows  the values to be copied into the matrix.
  //!
  //! \returns The constructed matrix.
  //!
  //! \throws
  //! LibsemigroupsException if `rows` does not represent a matrix of the
  //! correct dimensions.
  //!
  //! \throws
  //! LibsemigroupsException if the constructed matrix contains values that do
  //! not belong to the underlying semiring.
  //!
  //! \complexity
  //! \f$O(mn)\f$ where \f$m\f$ is the number of rows and
  //! \f$n\f$ is the number of columns of the matrix.
  template <typename Mat,
            typename Semiring,
            typename = std::enable_if_t<IsMatrix<Mat>>>
  // TODO(1) pass Semiring by reference, this is hard mostly due to the way the
  // tests are written, which is not optimal.
  Mat make(Semiring const* semiring,
           std::initializer_list<
               std::initializer_list<typename Mat::scalar_type>> const& rows) {
    Mat m(semiring, rows);
    matrix::throw_if_bad_entry(m);
    return m;
  }

  //! \ingroup matrix_group
  //!
  //! \brief Constructs a matrix (from std::vector of std::vector) and
  //! validates it.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Validates the arguments, constructs a matrix, and validates it.
  //!
  //! \tparam Mat the type of the matrix being constructed (must satisfy
  //! \ref IsMatrix).
  //!
  //! \tparam Semiring the type of the semiring where arithmetic is performed.
  //!
  //! \param semiring  a pointer to const semiring object.
  //!
  //! \param rows  the rows to be copied into the matrix.
  //!
  //! \returns The constructed matrix.
  //!
  //! \throws
  //! LibsemigroupsException if \p rows does not represent a
  //! matrix of the correct dimensions.
  //!
  //! \throws
  //! LibsemigroupsException if the constructed matrix
  //! contains values that do not belong to the underlying semiring.
  //!
  //! \complexity
  //! \f$O(mn)\f$ where \f$m\f$ is the number of rows and
  //! \f$n\f$ is the number of columns of the matrix.
  template <typename Mat,
            typename Semiring,
            typename = std::enable_if_t<IsMatrix<Mat>>>
  Mat make(Semiring const*                                            semiring,
           std::vector<std::vector<typename Mat::scalar_type>> const& rows) {
    Mat m(semiring, rows);
    matrix::throw_if_bad_entry(m);
    return m;
  }

  //! \ingroup matrix_group
  //!
  //! \brief Constructs a row and validates it.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This function constructs a row and validates it.
  //!
  //! \tparam Semiring the type of the semiring where arithmetic is performed.
  //!
  //! \param semiring  a pointer to const semiring object.
  //!
  //! \param row  the values to be copied into the row.
  //!
  //! \returns  the constructed row if valid.
  //!
  //! \throws
  //! LibsemigroupsException if the constructed row contains
  //! values that do not belong to the underlying semiring.
  //!
  //! \complexity
  //! \f$O(n)\f$ where \f$n\f$ is the number of columns of the matrix.
  template <typename Mat,
            typename Semiring,
            typename = std::enable_if_t<IsMatrix<Mat>>>
  Mat make(Semiring const*                                         semiring,
           std::initializer_list<typename Mat::scalar_type> const& row) {
    Mat m(semiring, row);
    matrix::throw_if_bad_entry(m);
    return m;
  }

  //! \ingroup matrix_group
  //!
  //! \brief Constructs a projective max-plus matrix (from
  //! std::initializer_list) and validates it.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Validates the arguments, constructs a matrix, and validates it.
  //!
  //! \tparam Mat the type of the matrix being constructed (must satisfy
  //! \ref IsProjMaxPlusMat<Mat>).
  //!
  //! \param rows  the values to be copied into the matrix.
  //!
  //! \returns The constructed matrix.
  //!
  //! \throws
  //! LibsemigroupsException if `rows` does not represent a matrix of the
  //! correct dimensions.
  //!
  //! \throws
  //! LibsemigroupsException if the constructed matrix contains values that do
  //! not belong to the underlying semiring.
  //!
  //! \complexity
  //! \f$O(mn)\f$ where \f$m\f$ is the number of rows and \f$n\f$ is the
  //! number of columns of the matrix.
  template <size_t R, size_t C, typename Scalar>
  ProjMaxPlusMat<R, C, Scalar>
  make(std::initializer_list<std::initializer_list<Scalar>> const& rows) {
    return ProjMaxPlusMat<R, C, Scalar>(
        make<ProjMaxPlusMat<R, C, Scalar>::underlying_matrix_type>(rows));
  }

  ////////////////////////////////////////////////////////////////////////
  // Printing etc...
  ////////////////////////////////////////////////////////////////////////

  //! \ingroup matrix_group
  //!
  //! \brief Insertion operator for row views.
  //!
  //! This function inserts a human readable representation of a row view into
  //! the string stream \p os.
  //!
  //! \param os the string stream.
  //! \param x the row view.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  template <typename S, typename T>
  std::ostringstream& operator<<(std::ostringstream&                os,
                                 detail::RowViewCommon<S, T> const& x) {
    os << "{";
    for (auto it = x.cbegin(); it != x.cend(); ++it) {
      os << *it;
      if (it != x.cend() - 1) {
        os << ", ";
      }
    }
    os << "}";
    return os;
  }

  //! \ingroup matrix_group
  //!
  //! \brief Insertion operator for matrices.
  //!
  //! This function inserts a human readable representation of a matrix into
  //! the string stream \p os.
  //!
  //! \tparam Mat the type of the argument \p x, must satisfy \ref
  //! IsMatrix<Mat>
  //!
  //! \param os the string stream.
  //!
  //! \param x the matrix.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  template <typename Mat>
  auto operator<<(std::ostringstream& os, Mat const& x)
      -> std::enable_if_t<IsMatrix<Mat>, std::ostringstream&> {
    size_t n = 0;
    if (x.number_of_rows() != 1) {
      os << "{";
    }
    for (auto&& r : matrix::rows(x)) {
      os << r;
      if (n != x.number_of_rows() - 1) {
        os << ", ";
      }
      n++;
    }
    if (x.number_of_rows() != 1) {
      os << "}";
    }
    return os;
  }

  //! \ingroup matrix_group
  //!
  //! \brief Returns a human readable representation of a matrix.
  //!
  //! This function returns a human readable representation of a matrix.
  //!
  //! \tparam Mat the type of the matrix, must satisfy \ref IsMatrix<Mat>.
  //!
  //! \param x the matrix.
  //! \param prefix the prefix for the returned string.
  //! \param short_name the short name of the type of matrix (default: `""`).
  //! \param braces the braces to use to delineate rows (default: `"{}"`).
  //! \param max_width the maximum width of the returned representation
  //! (default: \c 72).
  //!
  //! \throws LibsemigroupsException if \p braces does not have size \c 2.
  template <typename Mat>
  auto to_human_readable_repr(Mat const&         x,
                              std::string const& prefix,
                              std::string const& short_name = "",
                              std::string const& braces     = "{}",
                              size_t             max_width  = 72)
      -> std::enable_if_t<IsMatrix<Mat>, std::string> {
    if (braces.size() != 2) {
      LIBSEMIGROUPS_EXCEPTION(
          "the 4th argument (braces) must have size 2, found {}",
          braces.size());
    }

    size_t const R = x.number_of_rows();
    size_t const C = x.number_of_cols();

    std::vector<size_t> max_col_widths(C, 0);
    std::vector<size_t> row_widths(C, prefix.size() + 1);
    for (size_t r = 0; r < R; ++r) {
      for (size_t c = 0; c < C; ++c) {
        size_t width
            = detail::unicode_string_length(detail::entry_repr(x(r, c)));
        row_widths[r] += width;
        if (width > max_col_widths[c]) {
          max_col_widths[c] = width;
        }
      }
    }
    auto col_width
        = *std::max_element(max_col_widths.begin(), max_col_widths.end());
    // The total width if we pad the entries according to the widest column.
    auto const total_width = col_width * C + prefix.size() + 1;
    if (total_width > max_width) {
      // Padding according to the widest column is too wide!
      if (*std::max_element(row_widths.begin(), row_widths.end()) > max_width) {
        // If the widest row is too wide, then use the short name
        return fmt::format(
            "<{}x{} {}>", x.number_of_rows(), x.number_of_cols(), short_name);
      }
      // If the widest row is not too wide, then just don't pad the entries
      col_width = 0;
    }

    std::string result = fmt::format("{}", prefix);
    std::string rindent;
    auto const  lbrace = braces[0], rbrace = braces[1];
    if (R != 0 && C != 0) {
      result += lbrace;
      std::string csep = "";
      for (size_t r = 0; r < R; ++r) {
        result += fmt::format("{}{}", rindent, lbrace);
        rindent          = std::string(prefix.size() + 1, ' ');
        std::string csep = "";
        for (size_t c = 0; c < C; ++c) {
          result += fmt::format(
              "{}{:>{}}", csep, detail::entry_repr(x(r, c)), col_width);
          csep = ", ";
        }
        result += fmt::format("{}", rbrace);
        if (r != R - 1) {
          result += ",\n";
        }
      }
      result += rbrace;
    }
    result += ")";
    return result;
  }

  ////////////////////////////////////////////////////////////////////////
  // Adapters
  ////////////////////////////////////////////////////////////////////////

  //! \defgroup adapters_matrix_group Adapters for other matrices
  //!
  //! \brief Documentation for specialization of adapters for matrices
  //! satisfying \ref IsMatrix.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! This page contains the documentation of the functionality in
  //! `libsemigroups` that adapts  objects satisfying \ref IsMatrix for
  //! use with the main algorithms in `libsemigroups`.
  //!
  //! Other adapters specifically for \ref BMat objects are available \ref
  //! adapters_bmat_group "here" and for \ref BMat8 objects \ref
  //! adapters_bmat8_group "here".

  //! \ingroup adapters_matrix_group
  //!
  //! \brief Specialization of \ref Complexity for \ref IsMatrix.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Specialization of the adapter \ref Complexity for types \p Mat
  //! satisfying \ref IsMatrix<Mat>.
  //!
  //! \tparam  Mat the type of matrices.
  template <typename Mat>
  struct Complexity<Mat, std::enable_if_t<IsMatrix<Mat>>> {
    //! \brief Call operator.
    //!
    //! Returns `x.number_of_rows()` cubed.
    //!
    //! \param x  a matrix of type `Mat`.
    //!
    //! \returns  A value of type `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    constexpr size_t operator()(Mat const& x) const noexcept {
      return x.number_of_rows() * x.number_of_rows() * x.number_of_rows();
    }
  };

  //! \ingroup adapters_matrix_group
  //!
  //! \brief Specialization of \ref Degree for \ref IsMatrix.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Specialization of the adapter \ref Degree for types \p Mat satisfying
  //! \ref IsMatrix<Mat>.
  //!
  //! \tparam  Mat the type of matrices.
  template <typename Mat>
  struct Degree<Mat, std::enable_if_t<IsMatrix<Mat>>> {
    //! \brief Call operator
    //!
    //! This function returns `x.number_of_rows()`.
    //!
    //! \param x  a matrix of type `Mat`.
    //! \returns  The number of rows of \p x.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    constexpr size_t operator()(Mat const& x) const noexcept {
      return x.number_of_rows();
    }
  };

  //! \ingroup adapters_matrix_group
  //!
  //! \brief Specialization of \ref Hash for \ref IsMatrix.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Specialization of the adapter \ref Hash for types \p Mat satisfying
  //! \ref IsMatrix<Mat>.
  //!
  //! \tparam  Mat the type of matrices.
  template <typename Mat>
  struct Hash<Mat, std::enable_if_t<IsMatrix<Mat>>> {
    //! \brief Call operator
    //!
    //! This function returns `x.hash_value()`.
    //!
    //! \param x  a matrix of type `Mat`.
    //! \returns  A hash value for \p x.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    constexpr size_t operator()(Mat const& x) const {
      return x.hash_value();
    }
  };

  //! \ingroup adapters_matrix_group
  //!
  //! \brief Specialization of \ref IncreaseDegree for \ref IsMatrix.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Specialization of the adapter \ref IncreaseDegree for types \p Mat
  //! satisfying \ref IsMatrix<Mat>.
  //!
  //! \tparam Mat the type of matrices.
  //!
  //! \warning
  //! It is not possible to increase the degree of any of the types
  //! satisfying \ref IsMatrix, and as such the call operator of this type
  //! does nothing.
  template <typename Mat>
  struct IncreaseDegree<Mat, std::enable_if_t<IsMatrix<Mat>>> {
    //! \brief Call operator.
    //!
    //! This function does nothing and should not be used.
    constexpr void operator()(Mat&, size_t) const noexcept {
      // static_assert(false, "Cannot increase degree for Matrix");
      LIBSEMIGROUPS_ASSERT(false);
    }
  };

  //! \ingroup adapters_matrix_group
  //!
  //! \brief Specialization of \ref One for \ref IsMatrix.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Specialization of the adapter \ref One for types \p Mat satisfying
  //! \ref IsMatrix.
  //!
  //! \tparam Mat the type of matrices.
  template <typename Mat>
  struct One<Mat, std::enable_if_t<IsMatrix<Mat>>> {
    //! \brief Call operator.
    //!
    //! Returns the identity matrix.
    //!
    //! \param x  a matrix of type \p Mat.
    //! \returns  The identity matrix of type \p Mat.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(m ^ 2)\f$ where \f$m\f$ is the number of rows of the
    //! matrix \p x.
    inline Mat operator()(Mat const& x) const {
      return x.one();
    }
  };

  //! \ingroup adapters_matrix_group
  //!
  //! \brief Specialization of \ref Product for \ref IsMatrix.
  //!
  //! Defined in `matrix.hpp`.
  //!
  //! Specialization of the adapter \ref Product for types \p Mat satisfying
  //! \ref IsMatrix<Mat>.
  //!
  //! \tparam Mat the type of matrices.
  template <typename Mat>
  struct Product<Mat, std::enable_if_t<IsMatrix<Mat>>> {
    //! \brief Call operator.
    //!
    //! This function replaces the value of \p xy by the product of the
    //! matrices \p x and \p y.
    //!
    //! \param xy  a reference to a matrix of type \p Mat.
    //! \param x  a const reference to a matrix of type \p Mat.
    //! \param y  a const reference to a matrix of type \p Mat.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(m ^ 2)\f$ where \f$m\f$ is the number of rows of the
    //! matrix \p x.
    //!
    //! \warning
    //! This function only works for square matrices.
    inline void
    operator()(Mat& xy, Mat const& x, Mat const& y, size_t = 0) const {
      xy.product_inplace_no_checks(x, y);
    }
  };
}  // namespace libsemigroups

namespace std {
  template <size_t N,
            typename Mat,
            std::enable_if_t<libsemigroups::IsMatrix<Mat>>>
  inline void swap(Mat& x, Mat& y) noexcept {
    x.swap(y);
  }
}  // namespace std

#endif  // LIBSEMIGROUPS_MATRIX_HPP_
