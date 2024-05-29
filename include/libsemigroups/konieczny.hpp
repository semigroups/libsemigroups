//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2019-20 Finn Smith
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
// This file contains a generic implementation of Konieczny's algorithm,
// originally for computing subsemigroups of the boolean matrix monoid.

// TODO(later):
// * exception safety!
// * expose iterators to relevant things in D-classes, in particular elements
// * tpp file
// TODO(now):
// * rename _NC mem fns

#ifndef LIBSEMIGROUPS_KONIECZNY_HPP_
#define LIBSEMIGROUPS_KONIECZNY_HPP_

#ifndef PARSED_BY_DOXYGEN
#define NOT_PARSED_BY_DOXYGEN
#endif

#include <algorithm>      // for binary_search
#include <cstddef>        // for size_t
#include <set>            // for set
#include <type_traits>    // for is_pointer
#include <unordered_map>  // for unordered_map
#include <unordered_set>  // for unordered_set
#include <utility>        // for pair, make_pair
#include <vector>         // for vector

#include "action.hpp"     // for LeftAction, RightAction
#include "adapters.hpp"   // for Lambda, etc
#include "constants.hpp"  // for UNDEFINED
#include "debug.hpp"      // for LIBSEMIGROUPS_ASSERT
#include "exception.hpp"  // for LIBSEMIGROUPS_EXCEPTION
#include "runner.hpp"     // for Runner

#include "detail/pool.hpp"    // for detail::Pool
#include "detail/report.hpp"  // for report_default
#include "detail/timer.hpp"   // for Timer

#include "detail/bruidhinn-traits.hpp"  // for BruidhinnTraits

namespace libsemigroups {
  //! Defined in ``konieczny.hpp``.
  //!
  //! This is a traits class for use with Konieczny.
  //!
  //! \tparam Element the type of the elements.
  //!
  //! \sa Konieczny
  template <typename Element>
  struct KoniecznyTraits {
    //! The type of the elements of a Konieczny instance with const removed.
    using element_type = typename detail::BruidhinnTraits<Element>::value_type;

    //! The type of const elements of a Konieczny instance.
    using const_element_type =
        typename detail::BruidhinnTraits<Element>::const_value_type;

    //! Alias for \ref LambdaValue with template parameter \ref
    //! element_type.
    using lambda_value_type =
        typename ::libsemigroups::LambdaValue<element_type>::type;

    //! Alias for RhoValue with template parameter \ref element_type.
    using rho_value_type =
        typename ::libsemigroups::RhoValue<element_type>::type;

    //! Alias for \ref RhoValue with template parameter \ref element_type.
    using rank_state_type = typename ::libsemigroups::RankState<element_type>;

    //! The type of the orbit of the lambda values under ImageRightAction.
    //! \sa ImageRightAction and RightAction
    using lambda_orb_type
        = RightAction<element_type,
                      lambda_value_type,
                      ImageRightAction<element_type, lambda_value_type>>;

    //! The type of the orbit of the rho values under ImageLeftAction
    //! \sa ImageLeftAction and LeftAction
    using rho_orb_type
        = LeftAction<element_type,
                     rho_value_type,
                     ImageLeftAction<element_type, rho_value_type>>;

    //! \copydoc libsemigroups::Lambda
    using Lambda = ::libsemigroups::Lambda<element_type, lambda_value_type>;

    //! \copydoc libsemigroups::Rho
    using Rho = ::libsemigroups::Rho<element_type, rho_value_type>;

    //! \copydoc libsemigroups::Product
    using Product = ::libsemigroups::Product<element_type>;

    //! \copydoc libsemigroups::Rank
    using Rank = ::libsemigroups::Rank<element_type, rank_state_type>;

    //! \copydoc libsemigroups::One
    using One = ::libsemigroups::One<element_type>;

    //! \copydoc libsemigroups::Hash
    using ElementHash = ::libsemigroups::Hash<element_type>;

    //! \copydoc libsemigroups::EqualTo
    using EqualTo = ::libsemigroups::EqualTo<element_type>;

    //! \copydoc libsemigroups::Swap
    using Swap = ::libsemigroups::Swap<element_type>;

    //! \copydoc libsemigroups::Less
    using Less = ::libsemigroups::Less<element_type>;

    //! \copydoc libsemigroups::Degree
    using Degree = ::libsemigroups::Degree<element_type>;
  };

  //! Defined in ``konieczny.hpp``.
  //!
  //! The class template Konieczny implements %Konieczny's algorithm as
  //! described in the article 'Green's equivalences in finite semigroups of
  //! binary relations' by Janusz %Konieczny; see [here] for more details.
  //! This algorithm is similar to that of Lallement and McFadden; see [this]
  //! paper for more details. It differs in being applicable to subsemigroups of
  //! a non-regular semigroup, though is essentially the same algorithm for
  //! elements which happen to be regular.
  //!
  //! A Konieczny instance is defined by a generating set, and the main
  //! function is Konieczny::run, which implements
  //! %Konieczny's Algorithm. If Konieczny::run is invoked and
  //! Konieczny::finished returns \c true, then the size, partial order of
  //! \f$\mathscr{D}\f$-classes, and frames for each
  //! \f$\mathscr{D}\f$-class are known.
  //!
  //! \tparam Element the type of the elements of the semigroup (must not
  //! be a pointer).
  //!
  //! \tparam Traits the type of a traits class with the requirements of
  //! KoniecznyTraits.
  //!
  //! \sa KoniecznyTraits and DClass
  //!
  //! [here]: https://link.springer.com/article/10.1007/BF02573672
  //! [this]:
  //! https://www.sciencedirect.com/science/article/pii/S0747717108800570
  template <typename Element, typename Traits = KoniecznyTraits<Element>>
  class Konieczny : public Runner, private detail::BruidhinnTraits<Element> {
    // pointers are not currently supported
    static_assert(!std::is_pointer_v<Element>,
                  "Pointer types are not currently supported by Konieczny");

    ////////////////////////////////////////////////////////////////////////
    // Konieczny - aliases - private
    ////////////////////////////////////////////////////////////////////////

    using internal_element_type =
        typename detail::BruidhinnTraits<Element>::internal_value_type;
    using internal_const_element_type =
        typename detail::BruidhinnTraits<Element>::internal_const_value_type;
    using internal_const_reference =
        typename detail::BruidhinnTraits<Element>::internal_const_reference;
    using internal_reference =
        typename detail::BruidhinnTraits<Element>::internal_reference;

    using PairHash = Hash<std::pair<size_t, size_t>>;

   public:
    ////////////////////////////////////////////////////////////////////////
    // Konieczny - aliases - public
    ////////////////////////////////////////////////////////////////////////

    //! The type of elements.
    using element_type = typename Traits::element_type;

    //! The type of const elements.
    using const_element_type = typename Traits::const_element_type;

    //! Type of element const references.
    using const_reference =
        typename detail::BruidhinnTraits<Element>::const_reference;

    //! Type of indices of \f$\mathscr{D}\f$-classes.
    //!
    //! \sa cbegin_D_classes and cbegin_regular_D_classes.
    using D_class_index_type = size_t;

    //! The type of lambda values.
    using lambda_value_type = typename Traits::lambda_value_type;

    //! The type of the orbit of the lambda values.
    using lambda_orb_type = typename Traits::lambda_orb_type;

    //! The type of rho values.
    using rho_value_type = typename Traits::rho_value_type;

    //! The type of the orbit of the rho values.
    using rho_orb_type = typename Traits::rho_orb_type;

    //! \copydoc libsemigroups::Degree
    using Degree = typename Traits::Degree;

    //! \copydoc libsemigroups::EqualTo
    using EqualTo = typename Traits::EqualTo;

    //! \copydoc libsemigroups::Lambda
    using Lambda = typename Traits::Lambda;

    //! \copydoc libsemigroups::Less
    using Less = typename Traits::Less;

    //! \copydoc libsemigroups::One
    using One = typename Traits::One;

    //! \copydoc libsemigroups::Product
    using Product = typename Traits::Product;

    //! \copydoc libsemigroups::Rank
    using Rank = typename Traits::Rank;

    //! \copydoc libsemigroups::Rho
    using Rho = typename Traits::Rho;

    //! \copydoc libsemigroups::Swap
    using Swap = typename Traits::Swap;

   private:
    ////////////////////////////////////////////////////////////////////////
    // Konieczny - aliases - private
    ////////////////////////////////////////////////////////////////////////

    using lambda_orb_index_type    = typename lambda_orb_type::index_type;
    using rho_orb_index_type       = typename rho_orb_type::index_type;
    using rank_type                = size_t;
    using rank_state_type          = typename Traits::rank_state_type;
    using left_indices_index_type  = size_t;
    using right_indices_index_type = size_t;

    ////////////////////////////////////////////////////////////////////////
    // Konieczny - internal structs - private
    ////////////////////////////////////////////////////////////////////////

    struct InternalHash : private detail::BruidhinnTraits<Element> {
      size_t operator()(internal_const_element_type x) const {
        return Hash<Element>()(this->to_external_const(x));
      }
    };

    struct InternalEqualTo : private detail::BruidhinnTraits<Element> {
      bool operator()(internal_const_reference x,
                      internal_const_reference y) const {
        return EqualTo()(this->to_external_const(x),
                         this->to_external_const(y));
      }
    };

    struct InternalLess : private detail::BruidhinnTraits<Element> {
      bool operator()(internal_const_reference x,
                      internal_const_reference y) const {
        return Less()(this->to_external_const(x), this->to_external_const(y));
      }
    };

    struct InternalVecEqualTo : private detail::BruidhinnTraits<Element> {
      size_t operator()(std::vector<internal_element_type> const& x,
                        std::vector<internal_element_type> const& y) const {
        LIBSEMIGROUPS_ASSERT(x.size() == y.size());
        return std::equal(x.cbegin(), x.cend(), y.cbegin(), InternalEqualTo());
      }
    };

    struct InternalVecFree : private detail::BruidhinnTraits<Element> {
      void operator()(std::vector<internal_element_type> const& x) {
        for (auto it = x.cbegin(); it != x.cend(); ++it) {
          this->internal_free(*it);
        }
      }
    };

    struct OneParamLambda {
      lambda_value_type operator()(const_reference x) const {
        lambda_value_type lval;
        Lambda()(lval, x);
        return lval;
      }
    };

    struct OneParamRho {
      rho_value_type operator()(const_reference x) const {
        rho_value_type rval;
        Rho()(rval, x);
        return rval;
      }
    };

    struct InternalRank {
      template <typename SFINAE = size_t>
      auto operator()(void*, const_reference x)
          -> std::enable_if_t<std::is_void_v<typename rank_state_type::type>,
                              SFINAE> {
        return Rank()(x);
      }

      template <typename SFINAE = size_t>
      auto operator()(rank_state_type* state, const_reference x)
          -> std::enable_if_t<!std::is_void_v<typename rank_state_type::type>,
                              SFINAE> {
        return Rank()(*state, x);
      }
    };

    // RepInfo does not own its internal_element_type member
    struct RepInfo {
      RepInfo(D_class_index_type    D_idx,
              internal_element_type elt,
              lambda_orb_index_type lambda_idx,
              rho_orb_index_type    rho_idx)
          : _D_idx(D_idx),
            _elt(elt),
            _lambda_idx(lambda_idx),
            _rho_idx(rho_idx) {}

      D_class_index_type    _D_idx;
      internal_element_type _elt;
      lambda_orb_index_type _lambda_idx;
      rho_orb_index_type    _rho_idx;
    };

   public:
    ////////////////////////////////////////////////////////////////////////
    // Konieczny - constructor and destructor - public
    ////////////////////////////////////////////////////////////////////////

    //! Default constructor.
    //!
    //! This is the standard constructor for a Konieczny instance with
    //! unspecified generators.
    //!
    //! \parameters
    //! (None)
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \sa add_generator and add_generators
    Konieczny()
        : _adjoined_identity_contained(false),
          _D_classes(),
          _D_rels(),
          _data_initialised(false),
          _degree(UNDEFINED),
          _element_pool(),
          _gens(),
          _group_indices(),
          _group_indices_rev(),
          _lambda_orb(),
          _lambda_to_D_map(),
          _nonregular_reps(),
          _one(),
          _rank_state(nullptr),
          _ranks(),
          _regular_D_classes(),
          _reg_reps(),
          _reps_processed(0),
          _rho_orb(),
          _rho_to_D_map(),
          _run_initialised(false),
          _tmp_lambda_value1(),
          _tmp_lambda_value2(),
          _tmp_rho_value1(),
          _tmp_rho_value2() {
      _lambda_orb.cache_scc_multipliers(true);
      _rho_orb.cache_scc_multipliers(true);
    }

    //! Deleted.
    // TODO undelete
    Konieczny(Konieczny const&) = delete;

    //! Deleted.
    // TODO undelete
    Konieczny(Konieczny&&) = delete;

    //! Deleted.
    // TODO undelete
    Konieczny& operator=(Konieczny const&) = delete;

    //! Deleted.
    // TODO undelete
    Konieczny& operator=(Konieczny&&) = delete;

    //! Construct from generators.
    //!
    //! This function constructs a Konieczny instance generated by the
    //! specified container of generators.  There can be duplicate generators
    //! and although they do not count as distinct elements, they do count as
    //! distinct generators.  In other words, the generators are precisely (a
    //! copy of) \p gens in the same order they occur in \p gens.
    //!
    //! \param gens the generators represented by \c this.
    //!
    //! \throws LibsemigroupsException if any of the following hold:
    //! * \p gens is empty
    //! * Degree`{}(x)` != Degree`{}(y)` for some \c x, \c y in \p gens.
    explicit Konieczny(std::vector<element_type> const& gens) : Konieczny() {
      if (gens.empty()) {
        LIBSEMIGROUPS_EXCEPTION(
            "expected a positive number of generators, but got 0");
      }
      add_generators(gens.cbegin(), gens.cend());
      init_data();
    }

    ~Konieczny();

    ////////////////////////////////////////////////////////////////////////
    // Konieczny - forward class declarations - public/private
    ////////////////////////////////////////////////////////////////////////

    class DClass;

   private:
    class RegularDClass;
    class NonRegularDClass;

   public:
    ////////////////////////////////////////////////////////////////////////
    // Konieczny - member functions - public
    ////////////////////////////////////////////////////////////////////////

    //! Returns the number of generators.
    //!
    //! This  function returns the number of generators given to \c this.
    //! Note that there may be duplicate generators, and so \c this may have
    //! more generators than unique generators.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns A value of type \c size_t.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \sa add_generator and add_generators
    size_t number_of_generators() const noexcept {
      LIBSEMIGROUPS_ASSERT(!_gens.empty());
      return _gens.size() - 1;
    }

    //! Returns a const reference to the generator given by an index.
    //!
    //! This function returns a const reference to the \p pos generators of \c
    //! this.
    //!
    //! \param pos the index of the generator.
    //!
    //! \returns
    //! A value of type \ref const_reference
    //!
    //! \throws LibsemigroupsException if the value of \p pos is greater than
    //! number_of_generators().
    //!
    //! \complexity
    //! Constant.
    //!
    //! \note There may be duplicate generators, and so \c this may
    //! have more generators than unique generators.
    //!
    //! \sa add_generator and add_generators
    const_reference generator(size_t pos) const {
      if (pos >= number_of_generators()) {
        LIBSEMIGROUPS_EXCEPTION(
            "index out of bounds, expected value in [{}, {}) found {}",
            0,
            number_of_generators(),
            pos);
      }
      return this->to_external_const(_gens[pos]);
    }

    //! Returns the number of \f$\mathscr{D}\f$-classes.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers a full enumeration of the frames of every
    //! \f$\mathscr{D}\f$-class.
    size_t number_of_D_classes() {
      run();
      return std::distance(cbegin_D_classes(), cend_D_classes());
    }

    //! Returns the current number of \f$\mathscr{D}\f$-classes.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers no enumeration.
    size_t current_number_of_D_classes() const {
      return std::distance(cbegin_D_classes(), cend_D_classes());
    }

    //! Returns the number of regular \f$\mathscr{D}\f$-classes.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers a full enumeration of the frames of every
    //! \f$\mathscr{D}\f$-class.
    size_t number_of_regular_D_classes() {
      run();
      return current_number_of_regular_D_classes();
    }

    //! Returns the current number of regular \f$\mathscr{D}\f$-classes
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers no enumeration.
    size_t current_number_of_regular_D_classes() const {
      return std::distance(cbegin_regular_D_classes(),
                           cend_regular_D_classes());
    }

    //! Returns the number of \f$\mathscr{L}\f$-classes.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers a full enumeration of the frames of every
    //! \f$\mathscr{D}\f$-class.
    size_t number_of_L_classes() {
      run();
      return current_number_of_L_classes();
    }

    //! Returns the current number of \f$\mathscr{L}\f$-classes.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers no enumeration.
    size_t current_number_of_L_classes() const {
      size_t val = 0;
      std::for_each(
          cbegin_D_classes(), cend_D_classes(), [&val](DClass const& D) {
            val += D.number_of_L_classes();
          });
      return val;
    }

    //! Returns the number of regular \f$\mathscr{L}\f$-classes.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers a full enumeration of the frames of every
    //! \f$\mathscr{D}\f$-class.
    size_t number_of_regular_L_classes() {
      run();
      return current_number_of_regular_L_classes();
    }

    //! Returns the current number of regular \f$\mathscr{L}\f$-classes.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers no enumeration.
    size_t current_number_of_regular_L_classes() const {
      size_t val = 0;
      std::for_each(
          cbegin_regular_D_classes(),
          cend_regular_D_classes(),
          [&val](DClass const& D) { val += D.number_of_L_classes(); });
      return val;
    }

    //! Returns the number of \f$\mathscr{R}\f$-classes.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers a full enumeration of the frames of every
    //! \f$\mathscr{D}\f$-class.
    size_t number_of_R_classes() {
      run();
      return current_number_of_R_classes();
    }

    //! Returns the current number of regular \f$\mathscr{R}\f$-classes.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers no enumeration.
    size_t current_number_of_R_classes() const {
      size_t val = 0;
      std::for_each(
          cbegin_D_classes(), cend_D_classes(), [&val](DClass const& D) {
            val += D.number_of_R_classes();
          });
      return val;
    }

    //! Returns the number of regular \f$\mathscr{R}\f$-classes.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers a full enumeration of the frames of every
    //! \f$\mathscr{D}\f$-class.
    size_t number_of_regular_R_classes() {
      run();
      return current_number_of_regular_R_classes();
    }

    //! Returns the current number of regular \f$\mathscr{R}\f$-classes.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers no enumeration.
    size_t current_number_of_regular_R_classes() const {
      size_t val = 0;
      std::for_each(
          cbegin_regular_D_classes(),
          cend_regular_D_classes(),
          [&val](DClass const& D) { val += D.number_of_R_classes(); });
      return val;
    }

    //! Returns the number of \f$\mathscr{H}\f$-classes.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers a full enumeration of the frames of every
    //! \f$\mathscr{D}\f$-class.
    size_t number_of_H_classes() {
      run();
      return current_number_of_H_classes();
    }

    //! Returns the current number of \f$\mathscr{H}\f$-classes.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers no enumeration.
    size_t current_number_of_H_classes() const {
      size_t val = 0;
      std::for_each(
          cbegin_D_classes(), cend_D_classes(), [&val](DClass const& D) {
            val += (D.number_of_R_classes() * D.number_of_L_classes());
          });
      return val;
    }

    //! Returns the number of idempotents.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers a full enumeration of the frames of every
    //! \f$\mathscr{D}\f$-class.
    size_t number_of_idempotents() {
      run();
      return current_number_of_idempotents();
    }

    //! Returns the current number of idempotents.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers no enumeration.
    size_t current_number_of_idempotents() const {
      size_t val = 0;
      std::for_each(
          cbegin_regular_D_classes(),
          cend_regular_D_classes(),
          [&val](RegularDClass const& D) { val += D.number_of_idempotents(); });
      return val;
    }

    //! Returns the number of regular elements.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers a full enumeration of the frames of every
    //! \f$\mathscr{D}\f$-class.
    size_t number_of_regular_elements() {
      run();
      return current_number_of_regular_elements();
    }

    //! Returns the current number of regular elements.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers no enumeration.
    size_t current_number_of_regular_elements() const {
      size_t val = 0;
      std::for_each(cbegin_regular_D_classes(),
                    cend_regular_D_classes(),
                    [&val](DClass const& D) { val += D.size(); });
      return val;
    }

    //! Returns the size.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers a full enumeration of the frames of every
    //! \f$\mathscr{D}\f$-class.
    //!
    //! \sa current_size
    size_t size() {
      run();
      return current_size();
    }

    //! Returns the current size.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function triggers no enumeration.
    //!
    //! \sa \ref size.
    size_t current_size() const {
      size_t val = 0;
      std::for_each(cbegin_D_classes(),
                    cend_D_classes(),
                    [&val](DClass const& D) { val += D.size(); });
      return val;
    }

    //! Returns the degree of elements.
    //!
    //! All elements of a Konieczny must have the same degree, as computed by
    //! an instance of Degree; this function returns that degree.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \sa Degree
    size_t degree() const noexcept {
      return _degree;
    }

    //! Test membership of an element.
    //!
    //! Returns \c true if \p x belongs to \c this and \c false if it does
    //! not.
    //!
    //! \param x a const reference to a possible element.
    //!
    //! \returns A value of type \c bool.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function may trigger a (partial) enumeration.
    bool contains(const_reference x) {
      return Degree()(x) == degree()
             && get_containing_D_class(this->to_internal_const(x), true)
                    != UNDEFINED;
    }

    //! Returns the \f$\mathscr{D}\f$-class containing an element.
    //!
    //! \param x a const reference to a possible element.
    //!
    //! \returns A reference to DClass.
    //!
    //! \throws LibsemigroupsException if \p x does not belong to \c this.
    //!
    //! \note This involves computing as many frames for
    //! \f$\mathscr{D}\f$-classes  as necessary.
    DClass& D_class_of_element(const_reference x) {
      D_class_index_type i
          = get_containing_D_class(this->to_internal_const(x), true);
      if (i == UNDEFINED) {
        LIBSEMIGROUPS_EXCEPTION(
            "the argument does not belong to this semigroup!");
      }
      return *_D_classes[i];
    }

    //! Test regularity of an element.
    //!
    //! Returns \c true if \p x is a regular element and \c false if it is not.
    //!
    //! \param x a const reference to a possible element.
    //!
    //! \returns A value of type \c bool.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note
    //! This function triggers the computation of the orbits of the Lambda and
    //! Rho values, if they are not already known.
    bool is_regular_element(const_reference x) {
      return contains(x) && is_regular_element_NC(this->to_internal_const(x));
    }

    //! Add a copy of an element to the generators.
    //!
    //! It is possible, if perhaps not desirable,  to add the same generator
    //! multiple times.
    //!
    //! \param gen the generator to add.
    //!
    //! \returns (None)
    //!
    //! \throws LibsemigroupsException if any of the following hold:
    //! * the degree of \p x is incompatible with the existing degree
    //! * \ref started returns \c true
    void add_generator(const_reference gen) {
      add_generators(&gen, &gen + 1);
    }

    //! Add collection of generators from container.
    //!
    //! See \ref add_generator for a detailed description.
    //!
    //! \tparam T the type of the container for generators to add (must be a
    //! non-pointer type).
    //!
    //! \param coll the collection of generators to add.
    //!
    //! \returns (None)
    //!
    //! \throws LibsemigroupsException if any of the following hold:
    //! * the degree of \p x is incompatible with the existing degree.
    //! * \ref started returns \c true
    template <typename T>
    void add_generators(T const& coll) {
      static_assert(!std::is_pointer_v<T>,
                    "the template parameter T must not be a pointer");
      add_generators(coll.begin(), coll.end());
    }

    //! Add collection of generators from initializer list.
    //!
    //! See \ref add_generator for a detailed description.
    //!
    //! \param coll the collection of generators to add.
    //!
    //! \returns (None)
    //!
    //! \throws LibsemigroupsException if any of the following hold:
    //! * the degree of \p x is incompatible with the existing degree.
    //! * \ref started returns \c true
    void add_generators(std::initializer_list<const_element_type> coll) {
      add_generators<std::initializer_list<const_element_type>>(coll);
    }

    template <typename T>
    void add_generators(T const& first, T const& last);

    ////////////////////////////////////////////////////////////////////////
    // Konieczny - iterators - public
    ////////////////////////////////////////////////////////////////////////

    //! A type for const iterators through elements.
    using const_iterator
        = detail::BruidhinnConstIterator<element_type,
                                         std::vector<internal_element_type>>;

    //! Returns a const iterator pointing to the first generator.
    //!
    //! This  function does not trigger any enumeration; the iterator returned
    //! may be invalidated by any call to a non-const member function of the
    //! Konieczny class.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c const_iterator.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \sa cend_generators
    const_iterator cbegin_generators() const noexcept {
      return const_iterator(_gens.cbegin());
    }

    //! Returns a const iterator pointing to one past the last generator.
    //!
    //! This  function does not trigger any enumeration; the iterator returned
    //! may be invalidated by any call to a non-const member function of the
    //! Konieczny class.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c const_iterator.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \sa cbegin_generators
    const_iterator cend_generators() const noexcept {
      return const_iterator(_gens.cend() - 1);
    }

    // forward decl, defined in konieczny.tpp
    template <typename T>
    struct DClassIteratorTraits;

    //! Return type of \ref cbegin_D_classes and \ref cend_D_classes.
    //!
    //! Type for const random access iterators through the
    //! \f$\mathscr{D}\f$-classes, in the order they were enumerated.
    //!
    //! \sa const_regular_d_class_iterator.
    using const_d_class_iterator
        = detail::ConstIteratorStateless<DClassIteratorTraits<DClass>>;

    //! Returns a const iterator referring to a pointer to the first
    //! \f$\mathscr{D}\f$-class.
    //!
    //! This  function does not trigger any enumeration; the iterator returned
    //! may be invalidated by any call to a non-const member function of the
    //! Konieczny class.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c const_d_class_iterator.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    // not noexcept because operator++ isn't necessarily
    const_d_class_iterator cbegin_D_classes() const {
      auto it = _D_classes.cbegin();
      if (_run_initialised) {
        return const_d_class_iterator(it)
               + (_adjoined_identity_contained ? 0 : 1);
      } else {
        return const_d_class_iterator(it);
      }
    }

    //! Returns a const iterator referring to past the pointer to the last
    //! \f$\mathscr{D}\f$-class.
    //!
    //! This  function does not trigger any enumeration; the iterator returned
    //! may be invalidated by any call to a non-const member function of the
    //! Konieczny class.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c const_d_class_iterator.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    const_d_class_iterator cend_D_classes() const noexcept {
      return const_d_class_iterator(_D_classes.cend());
    }

    //! Return type of \ref cbegin_regular_D_classes and \ref
    //! cend_regular_D_classes.
    //!
    //! A type for const random access iterators through the regular
    //! \f$\mathscr{D}\f$-classes, in the order they were enumerated.
    //!
    //! \sa const_d_class_iterator.
    using const_regular_d_class_iterator
        = detail::ConstIteratorStateless<DClassIteratorTraits<RegularDClass>>;

    //! Returns a const iterator referring to a pointer to the first regular
    //! \f$\mathscr{D}\f$-class.
    //!
    //! This  function does not trigger any enumeration; the iterator returned
    //! may be invalidated by any call to a non-const member function of the
    //! Konieczny class.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c const_d_class_iterator.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \sa cbegin_rdc
    // not noexcept because operator++ isn't necessarily
    const_regular_d_class_iterator cbegin_regular_D_classes() const {
      auto it = _regular_D_classes.cbegin();
      if (_run_initialised) {
        return const_regular_d_class_iterator(it)
               + (_adjoined_identity_contained ? 0 : 1);
      } else {
        return const_regular_d_class_iterator(it);
      }
    }

    //! Shorter form of \ref cbegin_regular_D_classes.
    const_regular_d_class_iterator cbegin_rdc() const noexcept {
      return cbegin_regular_D_classes();
    }

    //! Returns a const iterator referring to past the pointer to the last
    //! regular \f$\mathscr{D}\f$-class.
    //!
    //! This  function does not trigger any enumeration; the iterator returned
    //! may be invalidated by any call to a non-const member function of the
    //! Konieczny class.
    //!
    //! \parameters
    //! (None)
    //!
    //! \returns
    //! A value of type \c const_d_class_iterator.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \sa \ref cend_rdc
    const_regular_d_class_iterator cend_regular_D_classes() const noexcept {
      return const_regular_d_class_iterator(_regular_D_classes.cend());
    }

    //! Shorter form of \ref cend_regular_D_classes.
    const_regular_d_class_iterator cend_rdc() const {
      return cend_regular_D_classes();
    }

   private:
    using PoolGuard = detail::PoolGuard<internal_element_type>;

    ////////////////////////////////////////////////////////////////////////
    // Konieczny - utility methods - private
    ////////////////////////////////////////////////////////////////////////

    // assumes its argument has valid lambda/rho values
    bool is_regular_element_NC(internal_const_reference x,
                               lambda_orb_index_type    lpos = UNDEFINED,
                               rho_orb_index_type       rpos = UNDEFINED) {
      LIBSEMIGROUPS_ASSERT(_lambda_orb.finished() && _rho_orb.finished());
      lpos = lpos != UNDEFINED ? lpos : get_lpos(x);
      rpos = rpos != UNDEFINED ? rpos : get_rpos(x);
      return get_lambda_group_index(x, lpos, rpos) != UNDEFINED;
    }

    // TODO use this func
    lambda_orb_index_type get_lpos(internal_const_reference x) const {
      Lambda()(_tmp_lambda_value1, this->to_external_const(x));
      return _lambda_orb.position(_tmp_lambda_value1);
    }

    // TODO use this func
    rho_orb_index_type get_rpos(internal_const_reference x) const {
      Rho()(_tmp_rho_value1, this->to_external_const(x));
      return _rho_orb.position(_tmp_rho_value1);
    }

    // Returns a lambda orb index corresponding to a group H-class in the R-
    // class of \p x.
    // asserts its argument has lambda/rho values in the orbits.
    // modifies _tmp_lambda_value1
    // modifies _tmp_rho_value1
    lambda_orb_index_type get_lambda_group_index(internal_const_reference x,
                                                 lambda_orb_index_type    lpos
                                                 = UNDEFINED,
                                                 rho_orb_index_type rpos
                                                 = UNDEFINED);

    // Finds a group index of a H-class in the L-class of \p x.
    // modifies _tmp_lambda_value1
    // modifies _tmp_rho_value1
    rho_orb_index_type get_rho_group_index(internal_const_reference x,
                                           lambda_orb_index_type    lpos
                                           = UNDEFINED,
                                           rho_orb_index_type rpos = UNDEFINED);

    //! Finds the idempotent in the H-class of \p x. Note that it is assumed
    //! that \p x is in a group H-class.
    // TODO(later): it must be possible to do better than this
    void idem_in_H_class(internal_reference       res,
                         internal_const_reference x) const;

    //! Finds an idempotent in the \f$\mathscr{D}\f$-class of \c x, if \c x is
    //! regular, and modifies \c x in place to be this idempotent
    // modifies _tmp_lambda_value1
    void make_idem(internal_reference x);

    //! Finds the group inverse of \p x in its H-class; i.e. the element \c y
    //! in the H-class of \p x such that <tt> xy = \p id</tt>. Will run
    //! forever if no such element exists.
    void group_inverse(internal_element_type&   res,
                       internal_const_reference id,
                       internal_const_reference x) const;

    //! Determines whether <tt>(x, y)</tt> forms a group index.
    // modifies _tmp_lambda_value and _tmp_rho_value
    bool is_group_index(internal_const_reference x,
                        internal_const_reference y,
                        lambda_orb_index_type    lpos = UNDEFINED,
                        rho_orb_index_type       rpos = UNDEFINED) const;

    // pass full_check = true to use the contains method of the D-classes
    // instead of the contains_NC
    D_class_index_type get_containing_D_class(internal_const_reference x,
                                              bool const full_check = false);

    void add_to_D_maps(D_class_index_type d);

    ////////////////////////////////////////////////////////////////////////
    // Konieczny - accessor member functions - private
    ////////////////////////////////////////////////////////////////////////

    void add_D_class(RegularDClass* D);
    void add_D_class(NonRegularDClass* D);

    auto cbegin_internal_generators() const noexcept {
      return _gens.cbegin();
    }

    auto cend_internal_generators() const noexcept {
      return _gens.cend();
    }

    detail::Pool<internal_element_type>& element_pool() const {
      return _element_pool;
    }

    size_t max_rank() const noexcept {
      if (_ranks.empty()) {
        return UNDEFINED;
      }
      return *_ranks.rbegin();
    }

    ////////////////////////////////////////////////////////////////////////
    // Konieczny - initialisation member functions - private
    ////////////////////////////////////////////////////////////////////////
    void init_run();

    void init_data();

    void init_rank_state_and_rep_vecs();

    void compute_orbs();

    ////////////////////////////////////////////////////////////////////////
    // Konieczny - validation member functions - private
    ////////////////////////////////////////////////////////////////////////

    void validate_element(const_reference x) const {
      size_t const n = Degree()(x);
      if (degree() != UNDEFINED && n != degree()) {
        LIBSEMIGROUPS_EXCEPTION(
            "element has degree {} but should have degree {}", n, degree());
      }
    }

    template <typename Iterator>
    void validate_element_collection(Iterator first, Iterator last) const;

    ////////////////////////////////////////////////////////////////////////
    // Konieczny - Runner methods - private
    ////////////////////////////////////////////////////////////////////////
    bool finished_impl() const override;
    void run_impl() override;
    void run_report();

    ////////////////////////////////////////////////////////////////////////
    // Konieczny - data - private
    ////////////////////////////////////////////////////////////////////////
    bool                                         _adjoined_identity_contained;
    std::vector<DClass*>                         _D_classes;
    std::vector<std::vector<D_class_index_type>> _D_rels;
    bool                                         _data_initialised;
    size_t                                       _degree;
    mutable detail::Pool<internal_element_type>  _element_pool;
    std::vector<internal_element_type>           _gens;
    std::unordered_map<std::pair<rho_orb_index_type, lambda_orb_index_type>,
                       lambda_orb_index_type,
                       PairHash>
        _group_indices;
    std::unordered_map<std::pair<rho_orb_index_type, lambda_orb_index_type>,
                       rho_orb_index_type,
                       PairHash>
                    _group_indices_rev;
    lambda_orb_type _lambda_orb;
    std::unordered_map<lambda_orb_index_type, std::vector<D_class_index_type>>
                                      _lambda_to_D_map;
    std::vector<std::vector<RepInfo>> _nonregular_reps;
    internal_element_type             _one;
    rank_state_type*                  _rank_state;
    std::set<rank_type>               _ranks;
    std::vector<RegularDClass*>       _regular_D_classes;
    std::vector<std::vector<RepInfo>> _reg_reps;
    size_t                            _reps_processed;
    rho_orb_type                      _rho_orb;
    std::unordered_map<rho_orb_index_type, std::vector<D_class_index_type>>
                              _rho_to_D_map;
    bool                      _run_initialised;
    mutable lambda_value_type _tmp_lambda_value1;
    mutable lambda_value_type _tmp_lambda_value2;
    mutable rho_value_type    _tmp_rho_value1;
    mutable rho_value_type    _tmp_rho_value2;
  };

  template <typename Element, typename Traits>
  bool Konieczny<Element, Traits>::finished_impl() const {
    return _ranks.empty() && _run_initialised;
  }

}  // namespace libsemigroups
#include "konieczny.tpp"
#endif  // LIBSEMIGROUPS_KONIECZNY_HPP_
