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

#ifndef LIBSEMIGROUPS_MATRIX_HPP_
#define LIBSEMIGROUPS_MATRIX_HPP_

#include <algorithm>      // for min
#include <array>          // for array
#include <bitset>         // for bitset
#include <cstddef>        // for size_t
#include <iosfwd>         // for ostringstream
#include <numeric>        // for inner_product
#include <type_traits>    // for false_type, is_signed, true_type
#include <unordered_map>  // for unordered_map
#include <unordered_set>  // for unordered_set
#include <vector>         // for vector

#include "adapters.hpp"    // for Degree
#include "bitset.hpp"      // for BitSet
#include "constants.hpp"   // for POSITIVE_INFINITY
#include "containers.hpp"  // for StaticVector1
#include "debug.hpp"       // for LIBSEMIGROUPS_ASSERT
#include "exception.hpp"   // for LIBSEMIGROUPS_EXCEPTION
#include "string.hpp"      // for detail::to_string

namespace libsemigroups {

  ////////////////////////////////////////////////////////////////////////
  // Detail
  ////////////////////////////////////////////////////////////////////////

  namespace detail {

    template <typename T>
    struct IsStdBitSetHelper : std::false_type {};

    template <size_t N>
    struct IsStdBitSetHelper<std::bitset<N>> : std::true_type {};

    struct MatrixPolymorphicBase {};

    template <typename T>
    struct IsMatrixHelper {
      static constexpr bool value
          = std::is_base_of<detail::MatrixPolymorphicBase, T>::value;
    };
  }  // namespace detail

  template <typename T>
  static constexpr bool IsStdBitSet = detail::IsStdBitSetHelper<T>::value;

  template <typename T>
  constexpr bool IsMatrix = detail::IsMatrixHelper<T>::value;

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

      scalar_type one() const noexcept {
        return static_cast<Subclass const*>(this)->one_impl();
      }

      scalar_type zero() const noexcept {
        return static_cast<Subclass const*>(this)->zero_impl();
      }

      Semiring const* semiring() const noexcept {
        return static_cast<Subclass const*>(this)->semiring_impl();
      }

     private:
      ////////////////////////////////////////////////////////////////////////
      // MatrixCommon - Semiring arithmetic - private
      ////////////////////////////////////////////////////////////////////////

      scalar_type plus(scalar_type x, scalar_type y) const noexcept {
        return static_cast<Subclass const*>(this)->plus_impl(y, x);
      }

      scalar_type prod(scalar_type x, scalar_type y) const noexcept {
        return static_cast<Subclass const*>(this)->prod_impl(y, x);
      }

     protected:
      ////////////////////////////////////////////////////////////////////////
      // MatrixCommon - Container functions - protected
      ////////////////////////////////////////////////////////////////////////

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
      MatrixCommon()                    = default;
      MatrixCommon(MatrixCommon const&) = default;
      MatrixCommon(MatrixCommon&&)      = default;
      MatrixCommon& operator=(MatrixCommon const&) = default;
      MatrixCommon& operator=(MatrixCommon&&) = default;

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

      static Subclass make(
          std::initializer_list<std::initializer_list<scalar_type>> const& il) {
        validate_args(il);
        Subclass m(il);
        validate(m);
        return m;
      }

      static Subclass make(std::initializer_list<scalar_type> const& il) {
        Subclass m(il);
        validate(m);
        return m;
      }

      static Subclass make(
          void const*,
          std::initializer_list<std::initializer_list<scalar_type>> const& il) {
        return make(il);
      }

      static Subclass make(void const*,
                           std::initializer_list<scalar_type> const& il) {
        return make(il);
      }

      virtual ~MatrixCommon() = default;

      // not noexcept because mem allocate is required
      Subclass identity() const {
        size_t const n = number_of_rows();
        Subclass     x(semiring(), n, n);
        std::fill(x.begin(), x.end(), zero());
        for (size_t r = 0; r < n; ++r) {
          x(r, r) = one();
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
        return !(*this == that);
      }

      // not noexcept because operator< isn't
      template <typename T>
      bool operator>(T const& that) const {
        return that < *this;
      }

      ////////////////////////////////////////////////////////////////////////
      // Attributes
      ////////////////////////////////////////////////////////////////////////

      // not noexcept because vector::operator[] isn't, and neither is
      // array::operator[]
      scalar_reference operator()(size_t r, size_t c) {
        return this->_container[r * number_of_cols() + c];
      }

      // not noexcept because vector::operator[] isn't, and neither is
      // array::operator[]
      scalar_const_reference operator()(size_t r, size_t c) const {
        return this->_container[r * number_of_cols() + c];
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
      void product_inplace(Subclass const& A, Subclass const& B) {
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
                zero(),
                [this](scalar_type x, scalar_type y) {
                  return this->plus(x, y);
                },
                [this](scalar_type x, scalar_type y) {
                  return this->prod(x, y);
                });
          }
        }
      }

      // not noexcept because iterator increment isn't
      void operator*=(scalar_type a) {
        for (auto it = _container.begin(); it < _container.end(); ++it) {
          *it = prod(*it, a);
        }
      }

      // not noexcept because vector::operator[] and array::operator[] aren't
      void operator+=(Subclass const& that) {
        LIBSEMIGROUPS_ASSERT(that.number_of_rows() == number_of_rows());
        LIBSEMIGROUPS_ASSERT(that.number_of_cols() == number_of_cols());
        for (size_t i = 0; i < _container.size(); ++i) {
          _container[i] = plus(_container[i], that._container[i]);
        }
      }

      void operator+=(RowView const& that) {
        LIBSEMIGROUPS_ASSERT(number_of_rows() == 1);
        RowView(*static_cast<Subclass const*>(this)) += that;
      }

      // TODO(later) uncomment and test (this works, just not tested or used
      // for anything, so because time is short commenting out for now)
      // void operator+=(scalar_type a) {
      //   for (auto it = _container.begin(); it < _container.end(); ++it) {
      //     *it = plus(*it, a);
      //   }
      // }

      // TODO(later) implement operator*=(Subclass const&)

      ////////////////////////////////////////////////////////////////////////
      // Arithmetic operators - not in-place
      ////////////////////////////////////////////////////////////////////////

      // not noexcept because operator+= isn't
      Subclass operator+(Subclass const& y) const {
        Subclass result(*static_cast<Subclass const*>(this));
        result += y;
        return result;
      }

      // not noexcept because product_inplace isn't
      Subclass operator*(Subclass const& y) const {
        Subclass result(*static_cast<Subclass const*>(this));
        result.product_inplace(*static_cast<Subclass const*>(this), y);
        return result;
      }

      // TODO(later) implement operator*(Scalar)
      // TODO(later) implement operator+(Scalar)

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
      void transpose() noexcept {
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

      ////////////////////////////////////////////////////////////////////////
      // Rows
      ////////////////////////////////////////////////////////////////////////

      // not noexcept because there's an allocation
      RowView row(size_t i) const {
        if (i >= number_of_rows()) {
          LIBSEMIGROUPS_EXCEPTION(
              "index out of range, expected value in [%llu, %llu), found %llu",
              static_cast<uint64_t>(0),
              static_cast<uint64_t>(number_of_rows()),
              static_cast<uint64_t>(i));
        }
        auto& container = const_cast<Container&>(_container);
        return RowView(static_cast<Subclass const*>(this),
                       container.begin() + i * number_of_cols(),
                       number_of_cols());
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
      template <typename T>
      static void validate_args(T const& m) {
        if (m.size() <= 1) {
          return;
        }
        uint64_t const C  = m.begin()->size();
        auto           it = std::find_if_not(
            m.begin() + 1, m.end(), [&C](typename T::const_reference r) {
              return r.size() == C;
            });
        if (it != m.end()) {
          LIBSEMIGROUPS_EXCEPTION(
              "invalid argument, expected every item to "
              "have length %llu, found %llu in entry %llu",
              C,
              static_cast<uint64_t>(it->size()),
              static_cast<uint64_t>(std::distance(m.begin(), it)));
        }
      }

      static void validate_args(
          std::initializer_list<std::initializer_list<scalar_type>> m) {
        validate_args<
            std::initializer_list<std::initializer_list<scalar_type>>>(m);
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
      MatrixDynamicDim(MatrixDynamicDim const&) = default;
      MatrixDynamicDim(MatrixDynamicDim&&)      = default;
      MatrixDynamicDim& operator=(MatrixDynamicDim const&) = default;
      MatrixDynamicDim& operator=(MatrixDynamicDim&&) = default;

      MatrixDynamicDim(size_t r, size_t c)
          : _number_of_cols(c), _number_of_rows(r) {}

      virtual ~MatrixDynamicDim() = default;

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

      // TODO(later) from here to the end of MatrixStaticArithmetic should be
      // private or protected
      using scalar_type = Scalar;

      static constexpr scalar_type plus_impl(scalar_type x,
                                             scalar_type y) noexcept {
        return PlusOp()(x, y);
      }

      static constexpr scalar_type prod_impl(scalar_type x,
                                             scalar_type y) noexcept {
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
      scalar_type plus(scalar_type x, scalar_type y) const noexcept {
        return static_cast<Subclass const*>(this)->plus_impl(y, x);
      }

      scalar_type prod(scalar_type x, scalar_type y) const noexcept {
        return static_cast<Subclass const*>(this)->prod_impl(y, x);
      }

     public:
      RowViewCommon()                     = default;
      RowViewCommon(RowViewCommon const&) = default;
      RowViewCommon(RowViewCommon&&)      = default;
      RowViewCommon& operator=(RowViewCommon const&) = default;
      RowViewCommon& operator=(RowViewCommon&&) = default;

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
          this_[i] = plus(this_[i], x[i]);
        }
      }

      // not noexcept because iterator arithmeic isn't
      void operator+=(scalar_type a) {
        for (auto& x : *this) {
          x = plus(x, a);
        }
      }

      // not noexcept because iterator arithmeic isn't
      void operator*=(scalar_type a) {
        for (auto& x : *this) {
          x = prod(x, a);
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
  // RowViews - classes for cheaply storing iterators to rows
  ////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////
  // StaticRowViews - static arithmetic
  ////////////////////////////////////////////////////////////////////////

  template <typename PlusOp,
            typename ProdOp,
            typename ZeroOp,
            typename OneOp,
            size_t C,
            typename Scalar>
  class StaticRowView final
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
    using scalar_type = Scalar;

    using detail::RowViewCommon<
        StaticMatrix_<1>,
        StaticRowView<PlusOp, ProdOp, ZeroOp, OneOp, C, Scalar>>::RowViewCommon;

    StaticRowView() = default;  // gcc-6 seems to require this

    template <size_t R>
    StaticRowView(StaticMatrix_<R> const*,
                  typename RowViewCommon::iterator it,
                  size_t const)
        : RowViewCommon(it) {}

    using RowViewCommon::size;

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

  template <typename PlusOp,
            typename ProdOp,
            typename ZeroOp,
            typename OneOp,
            typename Scalar>
  class DynamicRowView<PlusOp, ProdOp, ZeroOp, OneOp, Scalar> final
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
    using iterator = typename RowViewCommon::iterator;
    using Row      = typename DynamicMatrix_::Row;

    using detail::RowViewCommon<
        DynamicMatrix<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>,
        DynamicRowView<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>>::RowViewCommon;

    DynamicRowView() = default;
    DynamicRowView(DynamicMatrix_ const*, iterator const& it, size_t N)
        : RowViewCommon(it), _length(N) {}

    explicit DynamicRowView(Row const& r)
        : RowViewCommon(r), _length(r.number_of_cols()) {}

   private:
    size_t length_impl() const noexcept {
      return _length;
    }
    size_t _length;
  };

  ////////////////////////////////////////////////////////////////////////
  // DynamicRowViews - dynamic arithmetic
  ////////////////////////////////////////////////////////////////////////

  template <typename Semiring, typename Scalar>
  class DynamicRowView<Semiring, Scalar> final
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
    using scalar_type = typename DynamicMatrix_::scalar_type;
    using iterator    = typename RowViewCommon::iterator;
    using Row         = typename DynamicMatrix_::Row;

    using detail::RowViewCommon<
        DynamicMatrix<Semiring, Scalar>,
        DynamicRowView<Semiring, Scalar>>::RowViewCommon;

    DynamicRowView() = default;

    DynamicRowView(DynamicMatrix_ const* mat, iterator const& it, size_t)
        : RowViewCommon(it), _matrix(mat) {}

    explicit DynamicRowView(Row const& r) : RowViewCommon(r), _matrix(&r) {}

   private:
    size_t length_impl() const noexcept {
      return _matrix->number_of_cols();
    }

    scalar_type plus_impl(scalar_type x, scalar_type y) const noexcept {
      return _matrix->plus_impl(x, y);
    }

    scalar_type prod_impl(scalar_type x, scalar_type y) const noexcept {
      return _matrix->prod_impl(x, y);
    }

    DynamicMatrix_ const* _matrix;
  };

  ////////////////////////////////////////////////////////////////////////
  // StaticMatrix with compile time semiring arithmetic
  ////////////////////////////////////////////////////////////////////////

  template <typename PlusOp,
            typename ProdOp,
            typename ZeroOp,
            typename OneOp,
            size_t R,
            size_t C,
            typename Scalar>
  class StaticMatrix final
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

    using scalar_type      = typename MatrixCommon::scalar_type;
    using scalar_reference = typename MatrixCommon::scalar_reference;
    using scalar_const_reference =
        typename MatrixCommon::scalar_const_reference;

    using Row     = StaticMatrix<PlusOp, ProdOp, ZeroOp, OneOp, 1, C, Scalar>;
    using RowView = StaticRowView<PlusOp, ProdOp, ZeroOp, OneOp, C, Scalar>;

    using Plus = PlusOp;
    using Prod = ProdOp;
    using Zero = ZeroOp;
    using One  = OneOp;

    static constexpr size_t nr_rows = R;
    static constexpr size_t nr_cols = C;

    ////////////////////////////////////////////////////////////////////////
    // StaticMatrix - Constructors + destructor - public
    ////////////////////////////////////////////////////////////////////////

    explicit StaticMatrix(std::initializer_list<scalar_type> const& c)
        : MatrixCommon(c) {
      static_assert(R == 1,
                    "cannot construct Matrix from the given initializer list, "
                    "incompatible dimensions");
      LIBSEMIGROUPS_ASSERT(c.size() == C);
    }

    explicit StaticMatrix(
        std::initializer_list<std::initializer_list<scalar_type>> const& m)
        : MatrixCommon(m) {}

    StaticMatrix(size_t r, size_t c) : StaticMatrix() {
      (void) r;
      (void) c;
      LIBSEMIGROUPS_ASSERT(r == number_of_rows());
      LIBSEMIGROUPS_ASSERT(c == number_of_cols());
    }

    explicit StaticMatrix(std::vector<std::vector<scalar_type>> const& m)
        : MatrixCommon(m) {}

    explicit StaticMatrix(RowView const& rv) : MatrixCommon(rv) {
      static_assert(
          R == 1,
          "cannot construct Matrix with more than one row from RowView!");
    }

    StaticMatrix()                    = default;
    StaticMatrix(StaticMatrix const&) = default;
    StaticMatrix(StaticMatrix&&)      = default;
    StaticMatrix& operator=(StaticMatrix const&) = default;
    StaticMatrix& operator=(StaticMatrix&&) = default;

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

    ~StaticMatrix() = default;

    static StaticMatrix identity(size_t n = 0) {
      (void) n;
      static_assert(C == R, "cannot create non-square identity matrix");
      // If specified the value of n must equal R or otherwise weirdness will
      // ensue...
      LIBSEMIGROUPS_ASSERT(n == 0 || n == R);
      StaticMatrix x(R, R);
      std::fill(x.begin(), x.end(), ZeroOp()());
      for (size_t r = 0; r < R; ++r) {
        x(r, r) = OneOp()();
      }
      return x;
    }

    static StaticMatrix identity(void const* ptr, size_t n = 0) {
      (void) ptr;
      LIBSEMIGROUPS_ASSERT(ptr == nullptr);
      LIBSEMIGROUPS_ASSERT(n == 0 || n == R);
      return identity(n);
    }

    ////////////////////////////////////////////////////////////////////////
    // StaticMatrix - member function aliases - public
    ////////////////////////////////////////////////////////////////////////

    using MatrixCommon::operator();
    using MatrixCommon::begin;
    using MatrixCommon::number_of_cols;
    using MatrixCommon::number_of_rows;
    using MatrixCommon::operator==;
    using MatrixCommon::operator!=;
    using MatrixCommon::swap;

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

  template <typename PlusOp,
            typename ProdOp,
            typename ZeroOp,
            typename OneOp,
            typename Scalar>
  class DynamicMatrix<PlusOp, ProdOp, ZeroOp, OneOp, Scalar> final
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
    using scalar_type      = typename MatrixCommon::scalar_type;
    using scalar_reference = typename MatrixCommon::scalar_reference;
    using scalar_const_reference =
        typename MatrixCommon::scalar_const_reference;

    using Row     = DynamicMatrix;
    using RowView = DynamicRowView<PlusOp, ProdOp, ZeroOp, OneOp, Scalar>;

    DynamicMatrix()                     = default;
    DynamicMatrix(DynamicMatrix const&) = default;
    DynamicMatrix(DynamicMatrix&&)      = default;
    DynamicMatrix& operator=(DynamicMatrix const&) = default;
    DynamicMatrix& operator=(DynamicMatrix&&) = default;

    DynamicMatrix(size_t r, size_t c) : MatrixDynamicDim(r, c), MatrixCommon() {
      resize(number_of_rows(), number_of_cols());
    }

    explicit DynamicMatrix(std::initializer_list<scalar_type> const& c)
        : MatrixDynamicDim(1, c.size()), MatrixCommon(c) {}

    explicit DynamicMatrix(
        std::initializer_list<std::initializer_list<scalar_type>> const& m)
        : MatrixDynamicDim(m.size(), m.begin()->size()), MatrixCommon(m) {}

    explicit DynamicMatrix(std::vector<std::vector<scalar_type>> const& m)
        : MatrixDynamicDim(m.size(), m.begin()->size()), MatrixCommon(m) {}

    explicit DynamicMatrix(RowView const& rv)
        : MatrixDynamicDim(1, rv.size()), MatrixCommon(rv) {}

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

    ~DynamicMatrix();

    static DynamicMatrix identity(void const* ptr, size_t n) {
      (void) ptr;
      LIBSEMIGROUPS_ASSERT(ptr == nullptr);
      return identity(n);
    }

    static DynamicMatrix identity(size_t n) {
      DynamicMatrix x(n, n);
      std::fill(x.begin(), x.end(), ZeroOp()());
      for (size_t r = 0; r < n; ++r) {
        x(r, r) = OneOp()();
      }
      return x;
    }

    using MatrixCommon::begin;
    using MatrixCommon::identity;
    using MatrixCommon::number_of_cols;
    using MatrixCommon::number_of_rows;
    using MatrixCommon::resize;

    void swap(DynamicMatrix& that) noexcept {
      static_cast<MatrixDynamicDim&>(*this).swap(
          static_cast<MatrixDynamicDim&>(that));
      static_cast<MatrixCommon&>(*this).swap(static_cast<MatrixCommon&>(that));
    }
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

  template <typename Semiring, typename Scalar>
  class DynamicMatrix<Semiring, Scalar> final
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
    using scalar_type      = typename MatrixCommon::scalar_type;
    using scalar_reference = typename MatrixCommon::scalar_reference;
    using scalar_const_reference =
        typename MatrixCommon::scalar_const_reference;

    using Row     = DynamicMatrix;
    using RowView = DynamicRowView<Semiring, Scalar>;
    friend RowView;

    DynamicMatrix()                     = delete;
    DynamicMatrix(DynamicMatrix const&) = default;
    DynamicMatrix(DynamicMatrix&&)      = default;
    DynamicMatrix& operator=(DynamicMatrix const&) = default;
    DynamicMatrix& operator=(DynamicMatrix&&) = default;

    DynamicMatrix(Semiring const* semiring, size_t r, size_t c)
        : MatrixDynamicDim(r, c), MatrixCommon(), _semiring(semiring) {
      resize(number_of_rows(), number_of_cols());
    }

    explicit DynamicMatrix(Semiring const*                           semiring,
                           std::initializer_list<scalar_type> const& c)
        : MatrixDynamicDim(1, c.size()), MatrixCommon(c), _semiring(semiring) {}

    explicit DynamicMatrix(
        Semiring const* semiring,
        std::initializer_list<std::initializer_list<scalar_type>> const& m)
        : MatrixDynamicDim(m.size(), m.begin()->size()),
          MatrixCommon(m),
          _semiring(semiring) {}

    explicit DynamicMatrix(Semiring const* semiring,
                           std::vector<std::vector<scalar_type>> const& m)
        : MatrixDynamicDim(m.size(), (m.empty() ? 0 : m.begin()->size())),
          MatrixCommon(m),
          _semiring(semiring) {}

    explicit DynamicMatrix(RowView const& rv)
        : MatrixDynamicDim(1, rv.size()),
          MatrixCommon(rv),
          _semiring(rv._matrix->semiring()) {}

    static DynamicMatrix
    make(Semiring const*                                                 sr,
         std::initializer_list<std::initializer_list<scalar_type>> const il) {
      DynamicMatrix m(sr, il);
      validate(m);
      return m;
    }

    static DynamicMatrix make(Semiring const*                             sr,
                              std::vector<std::vector<scalar_type>> const v) {
      DynamicMatrix m(sr, v);
      validate(m);
      return m;
    }

    static DynamicMatrix make(Semiring const*                          sr,
                              std::initializer_list<scalar_type> const il) {
      DynamicMatrix m(sr, il);
      validate(m);
      return m;
    }

    // No static DynamicMatrix::identity(size_t n) because we need a semiring!
    static DynamicMatrix identity(Semiring const* sr, size_t n) {
      DynamicMatrix x(sr, n, n);
      std::fill(x.begin(), x.end(), x.zero());
      for (size_t r = 0; r < n; ++r) {
        x(r, r) = x.one();
      }
      return x;
    }

    ~DynamicMatrix();

    using MatrixCommon::begin;
    using MatrixCommon::identity;
    using MatrixCommon::number_of_cols;
    using MatrixCommon::number_of_rows;
    using MatrixCommon::resize;

    void swap(DynamicMatrix& that) noexcept {
      static_cast<MatrixDynamicDim&>(*this).swap(
          static_cast<MatrixDynamicDim&>(that));
      static_cast<MatrixCommon&>(*this).swap(static_cast<MatrixCommon&>(that));
      std::swap(_semiring, that._semiring);
    }

   private:
    scalar_type plus_impl(scalar_type x, scalar_type y) const noexcept {
      return _semiring->plus(x, y);
    }

    scalar_type prod_impl(scalar_type x, scalar_type y) const noexcept {
      return _semiring->prod(x, y);
    }

    scalar_type one_impl() const noexcept {
      return _semiring->one();
    }

    scalar_type zero_impl() const noexcept {
      return _semiring->zero();
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

  template <typename T>
  constexpr bool IsStaticMatrix = detail::IsStaticMatrixHelper<T>::value;

  template <typename T>
  constexpr bool IsDynamicMatrix = IsMatrix<T> && !IsStaticMatrix<T>;

  template <typename T>
  static constexpr bool IsMatWithSemiring
      = detail::IsMatWithSemiringHelper<T>::value;

  namespace detail {

    template <typename T>
    static constexpr bool IsTruncMat = IsTruncMatHelper<T>::value;

    template <typename Mat>
    void semiring_validate(Mat const& m) {
      if (IsMatWithSemiring<Mat> && m.semiring() == nullptr) {
        LIBSEMIGROUPS_EXCEPTION("the matrix pointer to semiring is nullptr!")
      }
    }
  }  // namespace detail

  template <typename Mat>
  auto matrix_threshold(Mat const&) noexcept
      -> std::enable_if_t<!detail::IsTruncMat<Mat>, typename Mat::scalar_type> {
    return UNDEFINED;
  }

  template <typename Mat>
  auto matrix_threshold(Mat const&) noexcept
      -> std::enable_if_t<detail::IsTruncMat<Mat> && !IsMatWithSemiring<Mat>,
                          typename Mat::scalar_type> {
    return detail::IsTruncMatHelper<Mat>::threshold;
  }

  template <typename Mat>
  auto matrix_threshold(Mat const& x) noexcept
      -> std::enable_if_t<detail::IsTruncMat<Mat> && IsMatWithSemiring<Mat>,
                          typename Mat::scalar_type> {
    return x.semiring()->threshold();
  }

  ////////////////////////////////////////////////////////////////////////
  // Boolean matrices - compile time semiring arithmetic
  ////////////////////////////////////////////////////////////////////////

  struct BooleanPlus {
    bool operator()(bool x, bool y) const noexcept {
      return x || y;
    }
  };

  struct BooleanProd {
    bool operator()(bool x, bool y) const noexcept {
      return x & y;
    }
  };

  struct BooleanOne {
    bool operator()() const noexcept {
      return true;
    }
  };

  struct BooleanZero {
    bool operator()() const noexcept {
      return false;
    }
  };

  // The use of `int` rather than `bool` as the scalar type for dynamic
  // boolean matrices is intentional, because the bit iterators implemented in
  // std::vector<bool> entail a significant performance penalty.
  using DynamicBMat
      = DynamicMatrix<BooleanPlus, BooleanProd, BooleanZero, BooleanOne, int>;

  template <size_t R, size_t C>
  using StaticBMat = StaticMatrix<BooleanPlus,
                                  BooleanProd,
                                  BooleanZero,
                                  BooleanOne,
                                  R,
                                  C,
                                  int>;

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

  template <typename T>
  static constexpr bool IsBMat = detail::IsBMatHelper<T>::value;

  template <typename Mat>
  auto validate(Mat const& m) -> std::enable_if_t<IsBMat<Mat>> {
    using scalar_type = typename Mat::scalar_type;
    auto it           = std::find_if_not(
        m.cbegin(), m.cend(), [](scalar_type x) { return x == 0 || x == 1; });
    if (it != m.cend()) {
      size_t r, c;
      std::tie(r, c) = m.coords(it);
      LIBSEMIGROUPS_EXCEPTION("invalid entry, expected 0 or 1 "
                              "but found %lld in entry (%llu, %llu)",
                              static_cast<int64_t>(*it),
                              static_cast<uint64_t>(r),
                              static_cast<uint64_t>(c));
    }
  }

  ////////////////////////////////////////////////////////////////////////
  // Integer matrices - compile time semiring arithmetic
  ////////////////////////////////////////////////////////////////////////

  template <typename Scalar>
  struct IntegerPlus {
    Scalar operator()(Scalar x, Scalar y) const noexcept {
      return x + y;
    }
  };

  template <typename Scalar>
  struct IntegerProd {
    Scalar operator()(Scalar x, Scalar y) const noexcept {
      return x * y;
    }
  };

  template <typename Scalar>
  struct IntegerZero {
    Scalar operator()() const noexcept {
      return 0;
    }
  };

  template <typename Scalar>
  struct IntegerOne {
    Scalar operator()() const noexcept {
      return 1;
    }
  };

  template <typename Scalar>
  using DynamicIntMat = DynamicMatrix<IntegerPlus<Scalar>,
                                      IntegerProd<Scalar>,
                                      IntegerZero<Scalar>,
                                      IntegerOne<Scalar>,
                                      Scalar>;

  template <size_t R, size_t C, typename Scalar>
  using StaticIntMat = StaticMatrix<IntegerPlus<Scalar>,
                                    IntegerProd<Scalar>,
                                    IntegerZero<Scalar>,
                                    IntegerOne<Scalar>,
                                    R,
                                    C,
                                    Scalar>;

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

  template <typename T>
  static constexpr bool IsIntMat = detail::IsIntMatHelper<T>::value;

  template <typename Mat>
  auto validate(Mat const&) -> std::enable_if_t<IsIntMat<Mat>> {}

  ////////////////////////////////////////////////////////////////////////
  // Max-plus matrices
  ////////////////////////////////////////////////////////////////////////

  // Static arithmetic
  template <typename Scalar>
  struct MaxPlusPlus {
    static_assert(std::is_signed<Scalar>::value,
                  "MaxPlus requires a signed integer type as parameter!");
    Scalar operator()(Scalar x, Scalar y) const noexcept {
      if (x == NEGATIVE_INFINITY) {
        return y;
      } else if (y == NEGATIVE_INFINITY) {
        return x;
      }
      return std::max(x, y);
    }
  };

  template <typename Scalar>
  struct MaxPlusProd {
    static_assert(std::is_signed<Scalar>::value,
                  "MaxPlus requires a signed integer type as parameter!");
    Scalar operator()(Scalar x, Scalar y) const noexcept {
      if (x == NEGATIVE_INFINITY || y == NEGATIVE_INFINITY) {
        return NEGATIVE_INFINITY;
      }
      return x + y;
    }
  };

  template <typename Scalar>
  struct MaxPlusZero {
    static_assert(std::is_signed<Scalar>::value,
                  "MaxPlus requires a signed integer type as parameter!");
    Scalar operator()() const noexcept {
      return NEGATIVE_INFINITY;
    }
  };

  template <typename Scalar>
  using DynamicMaxPlusMat = DynamicMatrix<MaxPlusPlus<Scalar>,
                                          MaxPlusProd<Scalar>,
                                          MaxPlusZero<Scalar>,
                                          IntegerZero<Scalar>,
                                          Scalar>;

  template <size_t R, size_t C, typename Scalar>
  using StaticMaxPlusMat = StaticMatrix<MaxPlusPlus<Scalar>,
                                        MaxPlusProd<Scalar>,
                                        MaxPlusZero<Scalar>,
                                        IntegerZero<Scalar>,
                                        R,
                                        C,
                                        Scalar>;

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

  template <typename T>
  static constexpr bool IsMaxPlusMat = detail::IsMaxPlusMatHelper<T>::value;

  template <typename Mat>
  auto validate(Mat const&) -> std::enable_if_t<IsMaxPlusMat<Mat>> {}

  ////////////////////////////////////////////////////////////////////////
  // Min-plus matrices
  ////////////////////////////////////////////////////////////////////////

  // Static arithmetic
  template <typename Scalar>
  struct MinPlusPlus {
    static_assert(std::is_signed<Scalar>::value,
                  "MinPlus requires a signed integer type as parameter!");
    Scalar operator()(Scalar x, Scalar y) const noexcept {
      if (x == POSITIVE_INFINITY) {
        return y;
      } else if (y == POSITIVE_INFINITY) {
        return x;
      }
      return std::min(x, y);
    }
  };

  template <typename Scalar>
  struct MinPlusProd {
    static_assert(std::is_signed<Scalar>::value,
                  "MinPlus requires a signed integer type as parameter!");
    Scalar operator()(Scalar x, Scalar y) const noexcept {
      if (x == POSITIVE_INFINITY || y == POSITIVE_INFINITY) {
        return POSITIVE_INFINITY;
      }
      return x + y;
    }
  };

  template <typename Scalar>
  struct MinPlusZero {
    static_assert(std::is_signed<Scalar>::value,
                  "MinPlus requires a signed integer type as parameter!");
    Scalar operator()() const noexcept {
      return POSITIVE_INFINITY;
    }
  };

  template <typename Scalar>
  using DynamicMinPlusMat = DynamicMatrix<MinPlusPlus<Scalar>,
                                          MinPlusProd<Scalar>,
                                          MinPlusZero<Scalar>,
                                          IntegerZero<Scalar>,
                                          Scalar>;

  template <size_t R, size_t C, typename Scalar>
  using StaticMinPlusMat = StaticMatrix<MinPlusPlus<Scalar>,
                                        MinPlusProd<Scalar>,
                                        MinPlusZero<Scalar>,
                                        IntegerZero<Scalar>,
                                        R,
                                        C,
                                        Scalar>;

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

  template <typename T>
  static constexpr bool IsMinPlusMat = detail::IsMinPlusMatHelper<T>::value;

  template <typename Mat>
  auto validate(Mat const&) -> std::enable_if_t<IsMinPlusMat<Mat>> {}

  ////////////////////////////////////////////////////////////////////////
  // Max-plus matrices with threshold
  ////////////////////////////////////////////////////////////////////////

  template <size_t T, typename Scalar>
  struct MaxPlusTruncProd {
    static_assert(std::is_signed<Scalar>::value,
                  "MaxPlus requires a signed integer type as parameter!");
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

  template <typename Scalar = int>
  class MaxPlusTruncSemiring final {
    static_assert(std::is_signed<Scalar>::value,
                  "MaxPlus requires a signed integer type as parameter!");

   public:
    MaxPlusTruncSemiring()                                     = delete;
    MaxPlusTruncSemiring(MaxPlusTruncSemiring const&) noexcept = default;
    MaxPlusTruncSemiring(MaxPlusTruncSemiring&&) noexcept      = default;
    MaxPlusTruncSemiring& operator                             =(
        MaxPlusTruncSemiring const&) noexcept = default;  // NOLINT(whitespace/line_length)
    MaxPlusTruncSemiring& operator=(MaxPlusTruncSemiring&&) noexcept = default;
    ~MaxPlusTruncSemiring()                                          = default;

    explicit MaxPlusTruncSemiring(Scalar threshold) : _threshold(threshold) {
      if (threshold < 0) {
        LIBSEMIGROUPS_EXCEPTION("expected non-negative value, found %lld",
                                static_cast<int64_t>(threshold));
      }
    }

    Scalar one() const noexcept {
      return 0;
    }

    Scalar zero() const noexcept {
      return NEGATIVE_INFINITY;
    }

    Scalar prod(Scalar x, Scalar y) const noexcept {
      LIBSEMIGROUPS_ASSERT((x >= 0 && x <= _threshold)
                           || x == NEGATIVE_INFINITY);
      LIBSEMIGROUPS_ASSERT((y >= 0 && y <= _threshold)
                           || y == NEGATIVE_INFINITY);
      if (x == NEGATIVE_INFINITY || y == NEGATIVE_INFINITY) {
        return NEGATIVE_INFINITY;
      }
      return std::min(x + y, _threshold);
    }

    Scalar plus(Scalar x, Scalar y) const noexcept {
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

    Scalar threshold() const noexcept {
      return _threshold;
    }

   public:
    Scalar const _threshold;
  };

  template <typename Scalar>
  using DynamicMaxPlusTruncMatSR
      = DynamicMatrix<MaxPlusTruncSemiring<Scalar>, Scalar>;

  template <size_t T, typename Scalar>
  using DynamicMaxPlusTruncMat = DynamicMatrix<MaxPlusPlus<Scalar>,
                                               MaxPlusTruncProd<T, Scalar>,
                                               MaxPlusZero<Scalar>,
                                               IntegerZero<Scalar>,
                                               Scalar>;

  template <size_t T, size_t R, size_t C, typename Scalar>
  using StaticMaxPlusTruncMat = StaticMatrix<MaxPlusPlus<Scalar>,
                                             MaxPlusTruncProd<T, Scalar>,
                                             MaxPlusZero<Scalar>,
                                             IntegerZero<Scalar>,
                                             R,
                                             C,
                                             Scalar>;

  template <size_t T = 0, size_t R = 0, size_t C = R, typename Scalar = int>
  using MaxPlusTruncMat = std::conditional_t<
      R == 0 || C == 0,
      std::conditional_t<T == 0,
                         DynamicMaxPlusTruncMatSR<Scalar>,
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
    struct IsMaxPlusTruncMatHelper<DynamicMaxPlusTruncMatSR<Scalar>>
        : std::true_type {
      static constexpr Scalar threshold = UNDEFINED;
    };
  }  // namespace detail

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

  template <typename Mat>
  auto validate(Mat const& m) -> std::enable_if_t<IsMaxPlusTruncMat<Mat>> {
    // Check that the semiring pointer isn't the nullptr if it shouldn't be
    detail::semiring_validate(m);

    using scalar_type   = typename Mat::scalar_type;
    scalar_type const t = matrix_threshold(m);
    auto it = std::find_if_not(m.cbegin(), m.cend(), [t](scalar_type x) {
      return x == NEGATIVE_INFINITY || (0 <= x && x <= t);
    });
    if (it != m.cend()) {
      uint64_t r, c;
      std::tie(r, c) = m.coords(it);
      LIBSEMIGROUPS_EXCEPTION("invalid entry, expected values in [0, %llu] "
                              "%s {-%s} but found %lld in entry (%llu, %llu)",
                              static_cast<uint64_t>(t),
                              u8"\u222A",
                              u8"\u221E",
                              static_cast<int64_t>(*it),
                              r,
                              c);
    }
  }

  ////////////////////////////////////////////////////////////////////////
  // Min-plus matrices with threshold
  ////////////////////////////////////////////////////////////////////////

  template <size_t T, typename Scalar>
  struct MinPlusTruncProd {
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

  template <typename Scalar = int>
  class MinPlusTruncSemiring final {
    static_assert(std::is_integral<Scalar>::value,
                  "MinPlus requires an integral type as parameter!");

   public:
    MinPlusTruncSemiring() noexcept                            = delete;
    MinPlusTruncSemiring(MinPlusTruncSemiring const&) noexcept = default;
    MinPlusTruncSemiring(MinPlusTruncSemiring&&) noexcept      = default;
    MinPlusTruncSemiring& operator                             =(
        MinPlusTruncSemiring const&) noexcept = default;  // NOLINT(whitespace/line_length)
    MinPlusTruncSemiring& operator=(MinPlusTruncSemiring&&) noexcept = default;
    ~MinPlusTruncSemiring()                                          = default;

    explicit MinPlusTruncSemiring(Scalar threshold) : _threshold(threshold) {
      if (std::is_signed<Scalar>::value && threshold < 0) {
        LIBSEMIGROUPS_EXCEPTION("expected non-negative value, found %lld",
                                static_cast<int64_t>(threshold));
      }
    }

    Scalar one() const noexcept {
      return 0;
    }

    // These mem fns (one and zero) aren't needed
    Scalar zero() const noexcept {
      return POSITIVE_INFINITY;
    }

    Scalar prod(Scalar x, Scalar y) const noexcept {
      LIBSEMIGROUPS_ASSERT((x >= 0 && x <= _threshold)
                           || x == POSITIVE_INFINITY);
      LIBSEMIGROUPS_ASSERT((y >= 0 && y <= _threshold)
                           || y == POSITIVE_INFINITY);
      if (x == POSITIVE_INFINITY || y == POSITIVE_INFINITY) {
        return POSITIVE_INFINITY;
      }
      return std::min(x + y, _threshold);
    }

    Scalar plus(Scalar x, Scalar y) const noexcept {
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

    Scalar threshold() const noexcept {
      return _threshold;
    }

   public:
    Scalar const _threshold;
  };

  template <typename Scalar>
  using DynamicMinPlusTruncMatSR
      = DynamicMatrix<MinPlusTruncSemiring<Scalar>, Scalar>;

  template <size_t T, typename Scalar>
  using DynamicMinPlusTruncMat = DynamicMatrix<MinPlusPlus<Scalar>,
                                               MinPlusTruncProd<T, Scalar>,
                                               MinPlusZero<Scalar>,
                                               IntegerZero<Scalar>,
                                               Scalar>;

  template <size_t T, size_t R, size_t C, typename Scalar>
  using StaticMinPlusTruncMat = StaticMatrix<MinPlusPlus<Scalar>,
                                             MinPlusTruncProd<T, Scalar>,
                                             MinPlusZero<Scalar>,
                                             IntegerZero<Scalar>,
                                             R,
                                             C,
                                             Scalar>;

  template <size_t T = 0, size_t R = 0, size_t C = R, typename Scalar = int>
  using MinPlusTruncMat = std::conditional_t<
      R == 0 || C == 0,
      std::conditional_t<T == 0,
                         DynamicMinPlusTruncMatSR<Scalar>,
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
    struct IsMinPlusTruncMatHelper<DynamicMinPlusTruncMatSR<Scalar>>
        : std::true_type {
      static constexpr Scalar threshold = UNDEFINED;
    };
  }  // namespace detail

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

  template <typename Mat>
  auto validate(Mat const& m) -> std::enable_if_t<IsMinPlusTruncMat<Mat>> {
    // Check that the semiring pointer isn't the nullptr if it shouldn't be
    detail::semiring_validate(m);

    using scalar_type   = typename Mat::scalar_type;
    scalar_type const t = matrix_threshold(m);
    auto it = std::find_if_not(m.cbegin(), m.cend(), [t](scalar_type x) {
      return x == POSITIVE_INFINITY || (0 <= x && x <= t);
    });
    if (it != m.cend()) {
      uint64_t r, c;
      std::tie(r, c) = m.coords(it);
      LIBSEMIGROUPS_EXCEPTION("invalid entry, expected values in [0, %llu] "
                              "%s {%s} but found %llu in entry (%llu, %llu)",
                              static_cast<uint64_t>(t),
                              u8"\u222A",
                              u8"\u221E",
                              static_cast<uint64_t>(*it),
                              r,
                              c);
    }
  }

  ////////////////////////////////////////////////////////////////////////
  // NTP matrices
  ////////////////////////////////////////////////////////////////////////

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

  // Static arithmetic
  template <size_t T, size_t P, typename Scalar>
  struct NTPPlus {
    Scalar operator()(Scalar x, Scalar y) const noexcept {
      return detail::thresholdperiod<T, P>(x + y);
    }
  };

  template <size_t T, size_t P, typename Scalar>
  struct NTPProd {
    Scalar operator()(Scalar x, Scalar y) const noexcept {
      return detail::thresholdperiod<T, P>(x * y);
    }
  };

  // Dynamic arithmetic
  template <typename Scalar = size_t>
  class NTPSemiring final {
   public:
    // Deleted to avoid uninitialised values of period and threshold.
    NTPSemiring()                   = delete;
    NTPSemiring(NTPSemiring const&) = default;
    NTPSemiring(NTPSemiring&&)      = default;

    NTPSemiring& operator=(NTPSemiring const&) = default;
    NTPSemiring& operator=(NTPSemiring&&) = default;

    ~NTPSemiring() = default;

    NTPSemiring(Scalar t, Scalar p) : _period(p), _threshold(t) {
      if (std::is_signed<Scalar>::value) {
        if (t < 0) {
          LIBSEMIGROUPS_EXCEPTION(
              "expected non-negative value for 1st argument, found %lld",
              static_cast<int64_t>(t));
        } else if (p <= 0) {
          LIBSEMIGROUPS_EXCEPTION(
              "expected non-negative value for 2nd argument, found %lld",
              static_cast<int64_t>(p));
        }
      }
    }

    Scalar one() const noexcept {
      return 1;
    }

    Scalar zero() const noexcept {
      return 0;
    }

    Scalar prod(Scalar x, Scalar y) const noexcept {
      LIBSEMIGROUPS_ASSERT(x >= 0 && x <= _period + _threshold - 1);
      LIBSEMIGROUPS_ASSERT(y >= 0 && y <= _period + _threshold - 1);
      return detail::thresholdperiod(x * y, _threshold, _period);
    }

    Scalar plus(Scalar x, Scalar y) const noexcept {
      LIBSEMIGROUPS_ASSERT(x >= 0 && x <= _period + _threshold - 1);
      LIBSEMIGROUPS_ASSERT(y >= 0 && y <= _period + _threshold - 1);
      return detail::thresholdperiod(x + y, _threshold, _period);
    }

    Scalar threshold() const noexcept {
      return _threshold;
    }

    Scalar period() const noexcept {
      return _period;
    }

   private:
    Scalar _period;
    Scalar _threshold;
  };

  template <typename Scalar>
  using DynamicNTPMatWithSemiring = DynamicMatrix<NTPSemiring<Scalar>, Scalar>;

  template <size_t T, size_t P, typename Scalar>
  using DynamicNTPMatWithoutSemiring = DynamicMatrix<NTPPlus<T, P, Scalar>,
                                                     NTPProd<T, P, Scalar>,
                                                     IntegerZero<Scalar>,
                                                     IntegerOne<Scalar>,
                                                     Scalar>;

  template <size_t T, size_t P, size_t R, size_t C, typename Scalar>
  using StaticNTPMat = StaticMatrix<NTPPlus<T, P, Scalar>,
                                    NTPProd<T, P, Scalar>,
                                    IntegerZero<Scalar>,
                                    IntegerOne<Scalar>,
                                    R,
                                    C,
                                    Scalar>;

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

  template <typename T>
  static constexpr bool IsNTPMat = detail::IsNTPMatHelper<T>::value;

  namespace detail {
    template <typename T>
    struct IsTruncMatHelper<T, std::enable_if_t<IsNTPMat<T>>> : std::true_type {
      static constexpr typename T::scalar_type threshold
          = IsNTPMatHelper<T>::threshold;
      static constexpr typename T::scalar_type period
          = IsNTPMatHelper<T>::period;
    };

  }  // namespace detail

  template <typename Mat>
  auto matrix_period(Mat const&) noexcept
      -> std::enable_if_t<!IsNTPMat<Mat>, typename Mat::scalar_type> {
    return UNDEFINED;
  }

  template <typename Mat>
  auto matrix_period(Mat const&) noexcept
      -> std::enable_if_t<IsNTPMat<Mat> && !IsMatWithSemiring<Mat>,
                          typename Mat::scalar_type> {
    return detail::IsTruncMatHelper<Mat>::period;
  }

  template <typename Mat>
  auto matrix_period(Mat const& x) noexcept
      -> std::enable_if_t<IsNTPMat<Mat> && IsMatWithSemiring<Mat>,
                          typename Mat::scalar_type> {
    return x.semiring()->period();
  }

  template <typename Mat>
  auto validate(Mat const& m) -> std::enable_if_t<IsNTPMat<Mat>> {
    // Check that the semiring pointer isn't the nullptr if it shouldn't be
    detail::semiring_validate(m);

    using scalar_type   = typename Mat::scalar_type;
    scalar_type const t = matrix_threshold(m);
    scalar_type const p = matrix_period(m);
    auto it = std::find_if_not(m.cbegin(), m.cend(), [t, p](scalar_type x) {
      return (0 <= x && x < p + t);
    });
    if (it != m.cend()) {
      uint64_t r, c;
      std::tie(r, c) = m.coords(it);
      LIBSEMIGROUPS_EXCEPTION("invalid entry, expected values in [0, %llu) "
                              "but found %llu in entry (%llu, %llu)",
                              static_cast<uint64_t>(p + t),
                              static_cast<uint64_t>(*it),
                              r,
                              c);
    }
  }  // namespace libsemigroups

  ////////////////////////////////////////////////////////////////////////
  // Projective max-plus matrices
  ////////////////////////////////////////////////////////////////////////

  namespace detail {
    template <typename T>
    class ProjMaxPlusMat final : MatrixPolymorphicBase {
     public:
      using scalar_type            = typename T::scalar_type;
      using scalar_reference       = typename T::scalar_reference;
      using scalar_const_reference = typename T::scalar_const_reference;
      using semiring_type          = void;

      using container_type = typename T::container_type;
      using iterator       = typename T::iterator;
      using const_iterator = typename T::const_iterator;

      using RowView = typename T::RowView;

      // Note that Rows are never normalised, and that's why we use the
      // underlying matrix Row type and not 1 x n ProjMaxPlusMat's instead
      // (since these will be normalised according to their entries, and
      // this might not correspond to the normalised entries of the matrix).
      using Row = typename T::Row;

      scalar_type one() const noexcept {
        return _underlying_mat.one();
      }

      scalar_type zero() const noexcept {
        return _underlying_mat.zero();
      }

      ////////////////////////////////////////////////////////////////////////
      // ProjMaxPlusMat - Constructors + destructor - public
      ////////////////////////////////////////////////////////////////////////

      ProjMaxPlusMat() : _is_normalized(false), _underlying_mat() {}
      ProjMaxPlusMat(ProjMaxPlusMat const&) = default;
      ProjMaxPlusMat(ProjMaxPlusMat&&)      = default;
      ProjMaxPlusMat& operator=(ProjMaxPlusMat const&) = default;
      ProjMaxPlusMat& operator=(ProjMaxPlusMat&&) = default;

      ProjMaxPlusMat(size_t r, size_t c)
          : _is_normalized(false), _underlying_mat(r, c) {}

      ProjMaxPlusMat(std::vector<std::vector<scalar_type>> const& m)
          : _is_normalized(false), _underlying_mat(m) {
        normalize();
      }

      ProjMaxPlusMat(
          std::initializer_list<std::initializer_list<scalar_type>> const& m)
          : ProjMaxPlusMat(
              std::vector<std::vector<scalar_type>>(m.begin(), m.end())) {}

      static ProjMaxPlusMat make(
          std::initializer_list<std::initializer_list<scalar_type>> const& il) {
        auto result = ProjMaxPlusMat(T::make(il));
        return result;
      }

      static ProjMaxPlusMat make(
          void const*,
          std::initializer_list<std::initializer_list<scalar_type>> const& il) {
        return make(il);
      }

      ~ProjMaxPlusMat() = default;

      ProjMaxPlusMat identity() const {
        auto result = ProjMaxPlusMat(_underlying_mat.identity());
        return result;
      }

      static ProjMaxPlusMat identity(size_t n) {
        return ProjMaxPlusMat(T::identity(n));
      }

      ////////////////////////////////////////////////////////////////////////
      // Comparison operators
      ////////////////////////////////////////////////////////////////////////

      bool operator==(ProjMaxPlusMat const& that) const {
        normalize();
        that.normalize();
        return _underlying_mat == that._underlying_mat;
      }

      bool operator!=(ProjMaxPlusMat const& that) const {
        return !(_underlying_mat == that._underlying_mat);
      }

      bool operator<(ProjMaxPlusMat const& that) const {
        normalize();
        that.normalize();
        return _underlying_mat < that._underlying_mat;
      }

      bool operator>(ProjMaxPlusMat const& that) const {
        return that < *this;
      }

      ////////////////////////////////////////////////////////////////////////
      // Attributes
      ////////////////////////////////////////////////////////////////////////

      // The following is commented out because it can't be used safely,
      // i.e. when is the matrix normalised again? scalar_reference
      scalar_reference operator()(size_t r, size_t c) {
        // to ensure the returned value is normalised
        normalize();
        // to ensure that the matrix is renormalised if the returned scalar is
        // assigned.
        _is_normalized = false;
        return _underlying_mat(r, c);
      }

      scalar_const_reference operator()(size_t r, size_t c) const {
        normalize();
        return _underlying_mat(r, c);
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

      void product_inplace(ProjMaxPlusMat const& A, ProjMaxPlusMat const& B) {
        _underlying_mat.product_inplace(A._underlying_mat, B._underlying_mat);
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
      ProjMaxPlusMat(T&& mat)
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

  template <size_t R, size_t C, typename Scalar>
  using StaticProjMaxPlusMat
      = detail::ProjMaxPlusMat<StaticMaxPlusMat<R, C, Scalar>>;

  template <typename Scalar>
  using DynamicProjMaxPlusMat
      = detail::ProjMaxPlusMat<DynamicMaxPlusMat<Scalar>>;

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

  template <typename T>
  static constexpr bool IsProjMaxPlusMat
      = detail::IsProjMaxPlusMatHelper<T>::value;

  template <typename Mat>
  auto validate(Mat const& m) -> std::enable_if_t<IsProjMaxPlusMat<Mat>> {
    validate(m.underlying_matrix());
  }

  namespace matrix_helpers {

    ////////////////////////////////////////////////////////////////////////
    // Matrix helpers - pow
    ////////////////////////////////////////////////////////////////////////

    // TODO(later) version that changes x in-place
    template <typename Mat>
    Mat pow(Mat const& x, typename Mat::scalar_type e) {
      using scalar_type = typename Mat::scalar_type;

      if (std::is_signed<scalar_type>::value && e < 0) {
        LIBSEMIGROUPS_EXCEPTION(
            "negative exponent, expected value >= 0, found %lld",
            static_cast<int64_t>(e));
      } else if (x.number_of_cols() != x.number_of_rows()) {
        LIBSEMIGROUPS_EXCEPTION("expected a square matrix, found %llux%llu",
                                static_cast<uint64_t>(x.number_of_rows()),
                                static_cast<uint64_t>(x.number_of_cols()));
      }

      if (e == 0) {
        return x.identity();
      }

      auto y = Mat(x);
      if (e == 1) {
        return y;
      }
      auto z = (e % 2 == 0 ? x.identity() : y);

      Mat tmp(x.number_of_rows(), x.number_of_cols());
      while (e > 1) {
        tmp.product_inplace(y, y);
        std::swap(y, tmp);
        e /= 2;
        if (e % 2 == 1) {
          tmp.product_inplace(z, y);
          std::swap(z, tmp);
        }
      }
      return z;
    }

    ////////////////////////////////////////////////////////////////////////
    // Matrix helpers - bitset_rows
    ////////////////////////////////////////////////////////////////////////

    // The main function
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

    // Helper
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

    ////////////////////////////////////////////////////////////////////////
    // Matrix helpers - bitset_row_basis
    ////////////////////////////////////////////////////////////////////////

    // This works with std::vector and StaticVector1, with value_type equal
    // to std::bitset and BitSet.
    template <typename Mat, typename Container>
    void bitset_row_basis(Container&& rows, std::decay_t<Container>& result) {
      using value_type = typename std::decay_t<Container>::value_type;
      static_assert(IsBMat<Mat>, "IsBMat<Mat> must be true!");
      static_assert(IsBitSet<value_type> || IsStdBitSet<value_type>,
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

    template <typename Mat, typename Container>
    std::decay_t<Container> bitset_row_basis(Container&& rows) {
      using value_type = typename std::decay_t<Container>::value_type;
      static_assert(IsBMat<Mat>, "IsBMat<Mat> must be true!");
      static_assert(IsBitSet<value_type> || IsStdBitSet<value_type>,
                    "Container::value_type must be BitSet or std::bitset");
      LIBSEMIGROUPS_ASSERT(rows.size() <= BitSet<1>::max_size());
      LIBSEMIGROUPS_ASSERT(rows.empty()
                           || rows[0].size() <= BitSet<1>::max_size());
      std::decay_t<Container> result;
      bitset_row_basis<Mat>(std::forward<Container>(rows), result);
      return result;
    }

    ////////////////////////////////////////////////////////////////////////
    // Matrix helpers - rows
    ////////////////////////////////////////////////////////////////////////

    template <typename Mat, typename = std::enable_if_t<IsDynamicMatrix<Mat>>>
    std::vector<typename Mat::RowView> rows(Mat const& x) {
      std::vector<typename Mat::RowView> container;
      x.rows(container);
      return container;
    }

    template <typename Mat, typename = std::enable_if_t<IsStaticMatrix<Mat>>>
    detail::StaticVector1<typename Mat::RowView, Mat::nr_rows>
    rows(Mat const& x) {
      detail::StaticVector1<typename Mat::RowView, Mat::nr_rows> container;
      x.rows(container);
      return container;
    }

    // Helper
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
    template <typename Mat>
    auto bitset_rows(Mat const& x) {
      static_assert(IsBMat<Mat>, "IsBMat<Mat> must be true!");
      LIBSEMIGROUPS_ASSERT(x.number_of_rows() <= BitSet<1>::max_size());
      LIBSEMIGROUPS_ASSERT(x.number_of_cols() <= BitSet<1>::max_size());
      size_t const M = detail::BitSetCapacity<Mat>::value;
      return bitset_rows<Mat, M, M>(std::move(rows(x)));
    }

    template <typename Mat, size_t M = detail::BitSetCapacity<Mat>::value>
    detail::StaticVector1<BitSet<M>, M> bitset_row_basis(Mat const& x) {
      static_assert(IsBMat<Mat>, "IsBMat<Mat> must be true!");
      LIBSEMIGROUPS_ASSERT(x.number_of_rows() <= BitSet<1>::max_size());
      LIBSEMIGROUPS_ASSERT(x.number_of_cols() <= BitSet<1>::max_size());
      detail::StaticVector1<BitSet<M>, M> result;
      bitset_row_basis<Mat>(std::move(bitset_rows(x)), result);
      return result;
    }

    template <typename Mat, typename Container>
    void bitset_row_basis(Mat const& x, Container& result) {
      using value_type = typename Container::value_type;
      static_assert(IsBMat<Mat>, "IsBMat<Mat> must be true!");
      static_assert(IsBitSet<value_type> || IsStdBitSet<value_type>,
                    "Container::value_type must be BitSet or std::bitset");
      LIBSEMIGROUPS_ASSERT(x.number_of_rows() <= BitSet<1>::max_size());
      LIBSEMIGROUPS_ASSERT(x.number_of_cols() <= BitSet<1>::max_size());
      bitset_row_basis<Mat>(std::move(bitset_rows(x)), result);
    }

    ////////////////////////////////////////////////////////////////////////
    // Matrix helpers - row_basis - MaxPlusTruncMat
    ////////////////////////////////////////////////////////////////////////

    template <typename Mat, typename Container>
    auto row_basis(Container&& views, std::decay_t<Container>& result)
        -> std::enable_if_t<IsMaxPlusTruncMat<Mat>> {
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
          std::fill(tmp1.begin(), tmp1.end(), tmp1.zero());
          for (size_t r2 = 0; r2 < r1; ++r2) {
            scalar_type max_scalar = matrix_threshold(tmp1);
            for (size_t c = 0; c < tmp1.number_of_cols(); ++c) {
              if (views[r2][c] == tmp1.zero()) {
                continue;
              }
              if (views[r1][c] >= views[r2][c]) {
                if (views[r1][c] != matrix_threshold(tmp1)) {
                  max_scalar
                      = std::min(max_scalar, views[r1][c] - views[r2][c]);
                }
              } else {
                max_scalar = tmp1.zero();
                break;
              }
            }
            if (max_scalar != tmp1.zero()) {
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
    template <typename Mat, typename Container>
    auto row_basis(Container&& views, std::decay_t<Container>& result)
        -> std::enable_if_t<IsBMat<Mat>> {
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

    // Row basis of rowspace of matrix <x> appended to <result>
    template <typename Mat,
              typename Container,
              typename = std::enable_if_t<IsMatrix<Mat>>>
    void row_basis(Mat const& x, Container& result) {
      row_basis<Mat>(std::move(rows(x)), result);
    }

    // Row basis of rowspace of matrix <x>
    template <typename Mat, typename = std::enable_if_t<IsDynamicMatrix<Mat>>>
    std::vector<typename Mat::RowView> row_basis(Mat const& x) {
      std::vector<typename Mat::RowView> container;
      row_basis(x, container);
      return container;
    }

    // Row basis of rowspace of matrix <x>
    template <typename Mat, typename = std::enable_if_t<IsStaticMatrix<Mat>>>
    detail::StaticVector1<typename Mat::RowView, Mat::nr_rows>
    row_basis(Mat const& x) {
      detail::StaticVector1<typename Mat::RowView, Mat::nr_rows> container;
      row_basis(x, container);
      return container;
    }

    // Row basis of rowspace of space spanned by <rows>
    template <typename Mat,
              typename Container,
              typename = std::enable_if_t<!IsMatrix<std::decay_t<Container>>>>
    std::decay_t<Container> row_basis(Container&& rows) {
      using value_type = typename std::decay_t<Container>::value_type;
      static_assert(IsMatrix<Mat>, "IsMatrix<Mat> must be true!");
      static_assert(std::is_same<value_type, typename Mat::RowView>::value,
                    "Container::value_type must be Mat::RowView");

      std::decay_t<Container> result;
      row_basis<Mat>(std::forward<Container>(rows), result);
      return result;
    }

    template <typename Mat,
              typename Container,
              typename = std::enable_if_t<IsMatrix<Mat>>>
    void row_basis_rows(Mat const& x, Container& result) {
      LIBSEMIGROUPS_ASSERT(result.empty());
      for (auto x : row_basis<Mat>(x)) {
        result.push_back(typename Mat::Row(x));
      }
    }

    template <typename Mat,
              typename Container,
              typename = std::enable_if_t<IsStaticMatrix<Mat>>>
    detail::StaticVector1<typename Mat::Row, Mat::nr_rows>
    row_basis_rows(Mat const& x) {
      detail::StaticVector1<typename Mat::Row, Mat::nr_rows> container;
      row_basis_rows(x, container);
      return container;
    }

    template <typename Mat,
              typename Container,
              typename = std::enable_if_t<IsDynamicMatrix<Mat>>>
    std::vector<typename Mat::Row> row_basis_rows(Mat const& x) {
      std::vector<typename Mat::Row> container;
      row_basis_rows(x, container);
      return container;
    }

    // This modifies the argument rows by moving out those in the basis
    template <typename Mat,
              typename Container,
              typename = std::enable_if_t<IsStaticMatrix<Mat>>>
    detail::StaticVector1<typename Mat::Row, Mat::nr_rows>
    row_basis_rows(Container&& rows) {
      using value_type = typename std::decay_t<Container>::value_type;
      static_assert(IsMatrix<Mat>, "IsMatrix<Mat> must be true!");
      static_assert(std::is_same<value_type, typename Mat::Row>::value,
                    "Container::value_type must be Mat::Row");

      detail::StaticVector1<typename Mat::RowView, Mat::nr_rows> rvs;
      std::unordered_map<typename Mat::scalar_type*, size_t>     lookup;

      for (size_t i = 0; i < rows.size(); ++i) {
        auto rv = typename Mat::RowView(rows[i]);
        rvs.push_back(rv);
        lookup.insert({&(*rv.begin()), i});
      }
      std::decay_t<Container> result;
      for (auto rv : row_basis<Mat>(rvs)) {
        // TODO: really???
        result.push_back(std::move(rows[lookup.at(&(*rv.begin()))]));
      }
      return result;
    }

    ////////////////////////////////////////////////////////////////////////
    // Matrix helpers - row_space_size
    ////////////////////////////////////////////////////////////////////////

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

  }  // namespace matrix_helpers

  ////////////////////////////////////////////////////////////////////////
  // Printing etc...
  ////////////////////////////////////////////////////////////////////////

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

  template <typename T>
  auto operator<<(std::ostringstream& os, T const& x)
      -> std::enable_if_t<IsMatrix<T>, std::ostringstream&> {
    size_t n = 0;
    if (x.number_of_rows() != 1) {
      os << "{";
    }
    for (auto&& r : matrix_helpers::rows(x)) {
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

  ////////////////////////////////////////////////////////////////////////
  // Adapters
  ////////////////////////////////////////////////////////////////////////

  template <typename T>
  struct Complexity<T, std::enable_if_t<IsMatrix<T>>> {
    constexpr size_t operator()(T const& x) const noexcept {
      return x.number_of_rows() * x.number_of_rows() * x.number_of_rows();
    }
  };

  template <typename T>
  struct Degree<T, std::enable_if_t<IsMatrix<T>>> {
    constexpr size_t operator()(T const& x) const noexcept {
      return x.number_of_rows();
    }
  };

  template <typename T>
  struct Hash<T, std::enable_if_t<IsMatrix<T>>> {
    constexpr size_t operator()(T const& x) const {
      return x.hash_value();
    }
  };

  template <typename T>
  struct IncreaseDegree<T, std::enable_if_t<IsMatrix<T>>> {
    constexpr void operator()(T&, size_t) const noexcept {
      // static_assert(false, "Cannot increase degree for Matrix");
      LIBSEMIGROUPS_ASSERT(false);
    }
  };

  template <typename T>
  struct One<T, std::enable_if_t<IsMatrix<T>>> {
    inline T operator()(T const& x) const {
      return x.identity();
    }
  };

  template <typename T>
  struct Product<T, std::enable_if_t<IsMatrix<T>>> {
    inline void operator()(T& xy, T const& x, T const& y, size_t = 0) const {
      xy.product_inplace(x, y);
    }
  };

  namespace matrix_helpers {
    template <typename T>
    struct RowSum {
      void operator()(T& res, T const& pt, T const& x) const {
        res = pt;
        res += x;
      }
    };
  }  // namespace matrix_helpers
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
