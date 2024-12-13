//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2019-2023 James D. Mitchell
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

// This file contains a base class for congruence-like classes.

#ifndef LIBSEMIGROUPS_CONG_INTF_HPP_
#define LIBSEMIGROUPS_CONG_INTF_HPP_

#include <cstddef>           // for size_t
#include <initializer_list>  // for initializer_list
#include <iterator>
#include <vector>  // for vector, operator==, vector<>::const_iter...

#include "exception.hpp"  // for LIBSEMIGROUPS_EXCEPTION
#include "ranges.hpp"     // for is_input_or_sink_v
#include "runner.hpp"     // for Runner
#include "types.hpp"      // for word_type, relation_type, letter_type, tril

// TODO(0) add "Defined in \c cong-intf.hpp" everywhere relevant in this file.

namespace libsemigroups {
  // Forward decls
  class Congruence;
  class ToddCoxeter;

  //! \defgroup cong_all_group Congruences
  //!
  //! This page contains links to the documentation for the classes Congruence
  //! and CongruenceInterface, and to helper functions for all of the classes in
  //! ``libsemigroups`` that are derived from CongruenceInterface.
  //! These classes are:
  //! * \ref Congruence
  //! * \ref Kambites
  //! * \ref KnuthBendix
  //! * \ref todd_coxeter_class_group "ToddCoxeter"

  //! \defgroup cong_all_classes_group Classes
  //! \ingroup cong_all_group
  //!
  //! \brief Documentation for the classes CongruenceInterface and Congruence.
  //!
  //! This page contains links to the documentation for the classes \ref
  //! CongruenceInterface and \ref Congruence.

  //! \ingroup cong_all_classes_group
  //!
  //! \brief Class collecting common aspects of classes representing
  //! congruences.
  //!
  //! Defined in ``cong-intf.hpp``.
  //!
  //! Every class for representing a congruence in ``libsemigroups`` is derived
  //! from CongruenceInterface, which holds the member functions and data that
  //! are common to all its derived classes.
  //!
  //! These classes are:
  //! * \ref Congruence
  //! * \ref Kambites
  //! * \ref KnuthBendix
  //! * \ref todd_coxeter_class_group "ToddCoxeter"
  class CongruenceInterface : public Runner {
    /////////////////////////////////////////////////////////////////////////
    // CongruenceInterface - data members - private
    /////////////////////////////////////////////////////////////////////////

    std::vector<word_type> _generating_pairs;
    congruence_kind        _type;

   protected:
    ////////////////////////////////////////////////////////////////////////////
    // CongruenceInterface - constructors + destructor - protected
    ////////////////////////////////////////////////////////////////////////////

    // Constructors + initializers are protected to prevent construction of
    // useless CongruenceInterface objects
    CongruenceInterface() = default;

    CongruenceInterface& init() {
      _generating_pairs.clear();
      Runner::init();
      return *this;
    }

    explicit CongruenceInterface(congruence_kind type) : CongruenceInterface() {
      init(type);
    }

    CongruenceInterface& init(congruence_kind type) {
      init();
      _type = type;
      return *this;
    }

    CongruenceInterface(CongruenceInterface const&)            = default;
    CongruenceInterface(CongruenceInterface&&)                 = default;
    CongruenceInterface& operator=(CongruenceInterface const&) = default;
    CongruenceInterface& operator=(CongruenceInterface&&)      = default;

    ////////////////////////////////////////////////////////////////////////////
    // CongruenceInterface - validation - protected
    ////////////////////////////////////////////////////////////////////////////

    template <typename Subclass, typename Iterator1, typename Iterator2>
    void throw_if_letter_out_of_bounds(Iterator1 first, Iterator2 last) const {
      static_cast<Subclass const*>(this)->throw_if_letter_out_of_bounds(first,
                                                                        last);
    }

   public:
    ~CongruenceInterface();
    ////////////////////////////////////////////////////////////////////////////
    // CongruenceInterface - public member functions
    ////////////////////////////////////////////////////////////////////////////

    //! \brief Returns the number of generating pairs.
    //!
    //! This function returns the number of generating pairs, which is the size
    //! of \ref generating_pairs divided by \f$2\f$.
    //!
    //! \returns
    //! The number of generating pairs.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] size_t number_of_generating_pairs() const noexcept {
      return _generating_pairs.size() / 2;
    }

    //! \brief Returns a const reference to the std::vector of generating
    //! pairs.
    //!
    //! This function returns the generating pairs of the congruence
    //! represented by any derived class of a CongruenceInterface. This is
    //! always a std::vector of \ref word_type, regardless of the
    //! type of the presentation used by the implementation in the derived
    //! class.
    //!
    //! \returns
    //! A const reference to the generating pairs.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    // TODO(0) should implement this in the derived classes, so that the return
    // type is native_word_type, not always word_type.
    [[nodiscard]] std::vector<word_type> const&
    generating_pairs() const noexcept {
      return _generating_pairs;
    }

    //! \brief The kind of the congruence (1- or 2-sided).
    //!
    //! This function returns the kind of the congruence represented by a
    //! derived class of CongruenceInterface. See \ref congruence_kind for
    //! details.
    //!
    //! \return The kind of the congruence (1- or 2-sided).
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] congruence_kind kind() const noexcept {
      return _type;
    }

   protected:
    // This is protected so that it is not possible to change the kind
    // arbitrarily.
    CongruenceInterface& kind(congruence_kind knd) {
      _type = knd;
      return *this;
    }

    ////////////////////////////////////////////////////////////////////////
    // CongruenceInterface - add_generating_pair
    ////////////////////////////////////////////////////////////////////////

    // The functions in this section are used as aliases in the derived
    // classes, and so can be protected here.
    template <typename Subclass,
              typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    Subclass& add_generating_pair_no_checks(Iterator1 first1,
                                            Iterator2 last1,
                                            Iterator3 first2,
                                            Iterator4 last2) {
      LIBSEMIGROUPS_ASSERT(!started());
      _generating_pairs.emplace_back(first1, last1);
      _generating_pairs.emplace_back(first2, last2);
      return static_cast<Subclass&>(*this);
    }

    template <typename Subclass,
              typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    Subclass& add_generating_pair(Iterator1 first1,
                                  Iterator2 last1,
                                  Iterator3 first2,
                                  Iterator4 last2) {
      throw_if_started();
      throw_if_letter_out_of_bounds<Subclass>(first1, last1);
      throw_if_letter_out_of_bounds<Subclass>(first2, last2);
      return static_cast<Subclass&>(*this).add_generating_pair_no_checks(
          first1, last1, first2, last2);
    }

    ////////////////////////////////////////////////////////////////////////
    // CongruenceInterface - contains
    ////////////////////////////////////////////////////////////////////////

    template <typename Subclass,
              typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    [[nodiscard]] tril currently_contains(Iterator1 first1,
                                          Iterator2 last1,
                                          Iterator3 first2,
                                          Iterator4 last2) const {
      throw_if_letter_out_of_bounds<Subclass>(first1, last1);
      throw_if_letter_out_of_bounds<Subclass>(first2, last2);
      return static_cast<Subclass const&>(*this).currently_contains_no_checks(
          first1, last1, first2, last2);
    }

    template <typename Subclass,
              typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    [[nodiscard]] bool contains(Iterator1 first1,
                                Iterator2 last1,
                                Iterator3 first2,
                                Iterator4 last2) {
      throw_if_letter_out_of_bounds<Subclass>(first1, last1);
      throw_if_letter_out_of_bounds<Subclass>(first2, last2);
      return static_cast<Subclass&>(*this).contains_no_checks(
          first1, last1, first2, last2);
    }

   private:
    void throw_if_started() const;
  };

  namespace congruence_interface {

    //! \defgroup cong_intf_helpers_group Generic congruence helpers
    //! \ingroup cong_all_group
    //!
    //! \brief Helper functions for subclasses of \ref CongruenceInterface.
    //!
    //! This page contains documentation for helper functions for the classes
    //! Congruence, Kambites, KnuthBendix, and \ref todd_coxeter_class_group
    //! "ToddCoxeter". The functions documented on this page belong to all of
    //! the namespaces \ref cong_intf_helpers_group "congruence_interface",
    //! ``congruence``, \ref kambites, \ref knuth_bendix, and \ref
    //! todd_coxeter_helpers_group "todd_coxeter".

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - add_generating_pair
    ////////////////////////////////////////////////////////////////////////

    //! \defgroup cong_intf_helpers_add_pair_group Add generating pairs
    //! \ingroup cong_intf_helpers_group
    //!
    //! \brief Add a generating pair using objects instead of iterators.
    //!
    //! This page contains the documentation of the functions
    //! \ref add_generating_pair and \ref add_generating_pair_no_checks
    //! which can be invoked with a variety of different argument types.
    //!
    //! @{

    //! \brief Helper for adding a generating pair of words.
    //!
    //! This function can be used to add a generating pair to the subclass \p
    //! ci of \ref CongruenceInterface using the objects themselves rather than
    //! using iterators.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam Word the type of the second and third parameters.
    //!
    //! \param ci the subclass of \ref CongruenceInterface.
    //! \param u the left hand side of the pair to add.
    //! \param v the right hand side of the pair to add.
    //!
    //! \return A reference to \p ci.
    //!
    //! \cong_intf_warn_assume_letters_in_bounds
    template <typename Subclass, typename Word>
    Subclass& add_generating_pair_no_checks(Subclass&   ci,
                                            Word const& u,
                                            Word const& v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.add_generating_pair_no_checks(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    //! \brief Helper for adding a generating pair of words
    //! (std::initializer_list).
    //!
    //! See \ref
    //! add_generating_pair_no_checks(Subclass&, Word const&, Word const&)
    //! for details.
    template <typename Subclass, typename Int>
    Subclass&
    add_generating_pair_no_checks(Subclass&                         ci,
                                  std::initializer_list<Int> const& u,
                                  std::initializer_list<Int> const& v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return add_generating_pair_no_checks<Subclass, std::vector<Int>>(
          ci, u, v);
    }

    //! \brief Helper for adding a generating pair of words
    //! (string literals).
    //!
    //! See \ref
    //! add_generating_pair_no_checks(Subclass&, Word const&, Word const&)
    //! for details.
    template <typename Subclass>
    Subclass& add_generating_pair_no_checks(Subclass&   ci,
                                            char const* u,
                                            char const* v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.add_generating_pair_no_checks(
          u, u + std::strlen(u), v, v + std::strlen(v));
    }

    // This version of the function catches the cases when u & v are not of the
    // same type but both convertible to string_view
    //! \brief Helper for adding a generating pair of words
    //! (std::string_view).
    //!
    //! See \ref
    //! add_generating_pair_no_checks(Subclass&, Word const&, Word const&)
    //! for details.
    template <typename Subclass>
    Subclass& add_generating_pair_no_checks(Subclass&        ci,
                                            std::string_view u,
                                            std::string_view v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.add_generating_pair_no_checks(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    //! \brief Helper for adding a generating pair of words.
    //!
    //! This function can be used to add a generating pair to the subclass \p
    //! ci of \ref CongruenceInterface using the objects themselves rather than
    //! using iterators.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam Word the type of the second and third parameters.
    //!
    //! \param ci the subclass of \ref CongruenceInterface.
    //! \param u the left hand side of the pair to add.
    //! \param v the right hand side of the pair to add.
    //!
    //! \return A reference to \p ci.
    //!
    //! \cong_intf_throws_if_letters_out_of_bounds
    template <typename Subclass, typename Word>
    Subclass& add_generating_pair(Subclass& ci, Word const& u, Word const& v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.add_generating_pair(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    //! \brief Helper for adding a generating pair of words
    //! (std::initializer_list).
    //!
    //! See \ref
    //! add_generating_pair(Subclass&, Word const&, Word const&)
    //! for details.
    template <typename Subclass, typename Int>
    Subclass& add_generating_pair(Subclass&                         ci,
                                  std::initializer_list<Int> const& u,
                                  std::initializer_list<Int> const& v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return add_generating_pair<Subclass, std::vector<Int>>(ci, u, v);
    }

    //! \brief Helper for adding a generating pair of words
    //! (string literals).
    //!
    //! See \ref
    //! add_generating_pair(Subclass&, Word const&, Word const&)
    //! for details.
    template <typename Subclass>
    Subclass& add_generating_pair(Subclass& ci, char const* u, char const* v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.add_generating_pair(
          u, u + std::strlen(u), v, v + std::strlen(v));
    }

    // This version of the function catches the cases when u & v are not of the
    // same type but both convertible to string_view
    //! \brief Helper for adding a generating pair of words
    //! (std::string_view).
    //!
    //! See \ref
    //! add_generating_pair(Subclass&, Word const&, Word const&)
    //! for details.
    template <typename Subclass>
    Subclass& add_generating_pair(Subclass&        ci,
                                  std::string_view u,
                                  std::string_view v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.add_generating_pair_no_checks(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    //! @}

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - currently_contains_no_checks
    ////////////////////////////////////////////////////////////////////////

    //! \defgroup cong_intf_helpers_contains_group Containment
    //! \ingroup cong_intf_helpers_group
    //!
    //! \brief Check containment of a pair of words in a congruence.
    //!
    //! This page contains the documentation of the functions
    //! \ref currently_contains_no_checks; \ref currently_contains;
    //! \ref contains_no_checks; and \ref contains which can be invoked with a
    //! variety of different argument types.
    //!
    //! Functions with the prefix ``currently_`` do not perform any enumeration
    //! of the \ref CongruenceInterface derived class instances; and those with
    //! the suffix ``_no_checks`` do not check that the input words are valid.
    //!
    //! @{

    //! \brief Check containment of a pair of words.
    //!
    //! This function checks whether or not the words \p u and \p v are already
    //! known to be contained in the congruence represented by a \ref
    //! CongruenceInterface derived class instance \p ci. This function performs
    //! no enumeration of \p ci, so it is possible for the words to be contained
    //! in the congruence, but that this is not currently known.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam Word the type of the second and third parameters.
    //!
    //! \param ci the subclass of \ref CongruenceInterface.
    //! \param u the left hand side of the pair to add.
    //! \param v the right hand side of the pair to add.
    //!
    //! \returns
    //! * tril::TRUE if the words are known to belong to the congruence;
    //! * tril::FALSE if the words are known to not belong to the congruence;
    //! * tril::unknown otherwise.
    //!
    //! \cong_intf_warn_assume_letters_in_bounds
    template <typename Subclass, typename Word>
    [[nodiscard]] tril currently_contains_no_checks(Subclass const& ci,
                                                    Word const&     u,
                                                    Word const&     v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.currently_contains_no_checks(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    // This version of the function catches the cases when u & v are not of the
    // same type but both convertible to string_view
    //! \brief Helper for checking containment of a pair of words
    //! (std::string_view).
    //!
    //! See \ref
    //! currently_contains_no_checks(Subclass const&, Word const&, Word const&)
    //! for details.
    template <typename Subclass>
    [[nodiscard]] tril currently_contains_no_checks(Subclass const&  ci,
                                                    std::string_view u,
                                                    std::string_view v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.currently_contains_no_checks(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    //! \brief Helper for checking containment of a pair of words
    //! (string literal).
    //!
    //! See \ref
    //! currently_contains_no_checks(Subclass const&, Word const&, Word const&)
    //! for details.
    template <typename Subclass>
    [[nodiscard]] tril currently_contains_no_checks(Subclass const& ci,
                                                    char const*     u,
                                                    char const*     v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.currently_contains_no_checks(
          u, u + std::strlen(u), v, v + std::strlen(v));
    }

    //! \brief Helper for checking containment of a pair of words
    //! (std::initializer_list).
    //!
    //! See \ref
    //! currently_contains_no_checks(Subclass const&, Word const&, Word const&)
    //! for details.
    template <typename Subclass,
              typename Int = typename Subclass::native_letter_type>
    [[nodiscard]] tril
    currently_contains_no_checks(Subclass const&                   ci,
                                 std::initializer_list<Int> const& u,
                                 std::initializer_list<Int> const& v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return currently_contains_no_checks<Subclass, std::initializer_list<Int>>(
          ci, u, v);
    }

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - currently_contains
    ////////////////////////////////////////////////////////////////////////

    //! \brief Check containment of a pair of words.
    //!
    //! This function checks whether or not the words \p u and \p v are already
    //! known to be contained in the congruence represented by a \ref
    //! CongruenceInterface derived class instance \p ci. This function performs
    //! no enumeration of \p ci, so it is possible for the words to be contained
    //! in the congruence, but that this is not currently known.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam Word the type of the second and third parameters.
    //!
    //! \param ci the subclass of \ref CongruenceInterface.
    //! \param u the left hand side of the pair to add.
    //! \param v the right hand side of the pair to add.
    //!
    //! \returns
    //! * tril::TRUE if the words are known to belong to the congruence;
    //! * tril::FALSE if the words are known to not belong to the congruence;
    //! * tril::unknown otherwise.
    //!
    //! \cong_intf_throws_if_letters_out_of_bounds
    template <typename Subclass, typename Word>
    [[nodiscard]] tril currently_contains(Subclass const& ci,
                                          Word const&     u,
                                          Word const&     v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.currently_contains(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    // This version of the function catches the cases when u & v are not of the
    // same type but both convertible to string_view
    //! \brief Helper for checking containment of a pair of words
    //! (std::string_view).
    //!
    //! See \ref
    //! currently_contains(Subclass const&, Word const&, Word const&)
    //! for details.
    template <typename Subclass>
    [[nodiscard]] tril currently_contains(Subclass const&  ci,
                                          std::string_view u,
                                          std::string_view v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.currently_contains(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    //! \brief Helper for checking containment of a pair of words
    //! (string literal).
    //!
    //! See \ref
    //! currently_contains(Subclass const&, Word const&, Word const&)
    //! for details.
    template <typename Subclass>
    [[nodiscard]] tril currently_contains(Subclass const& ci,
                                          char const*     u,
                                          char const*     v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.currently_contains_no_checks(
          u, u + std::strlen(u), v, v + std::strlen(v));
    }

    //! \brief Helper for checking containment of a pair of words
    //! (std::initializer_list).
    //!
    //! See \ref
    //! currently_contains(Subclass const&, Word const&, Word const&)
    //! for details.
    template <typename Subclass,
              typename Int = typename Subclass::native_letter_type>
    [[nodiscard]] tril currently_contains(Subclass const&                   ci,
                                          std::initializer_list<Int> const& u,
                                          std::initializer_list<Int> const& v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return currently_contains<Subclass, std::initializer_list<Int>>(ci, u, v);
    }

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - contains_no_checks
    ////////////////////////////////////////////////////////////////////////

    //! \brief Check containment of a pair of words.
    //!
    //! This function checks whether or not the words \p u and \p v are
    //! contained in the congruence represented by the instance \p ci of a
    //! derived class of \ref CongruenceInterface. This function triggers a
    //! full enumeration of \p ci, which may never terminate.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam Word the type of the second and third parameters.
    //!
    //! \param ci the subclass of \ref CongruenceInterface.
    //! \param u the left hand side of the pair to add.
    //! \param v the right hand side of the pair to add.
    //!
    //! \returns Whether or not the pair belongs to the congruence.
    //!
    //! \cong_intf_warn_assume_letters_in_bounds
    template <typename Subclass, typename Word>
    [[nodiscard]] bool contains_no_checks(Subclass&   ci,
                                          Word const& u,
                                          Word const& v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.contains_no_checks(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    // This version of the function catches the cases when u & v are not of the
    // same type but both convertible to string_view
    //! \brief Helper for checking containment of a pair of words
    //! (std::string_view).
    //!
    //! See \ref
    //! contains_no_checks(Subclass&, Word const&, Word const&)
    //! for details.
    template <typename Subclass>
    [[nodiscard]] bool contains_no_checks(Subclass&        ci,
                                          std::string_view u,
                                          std::string_view v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.contains_no_checks(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    //! \brief Helper for checking containment of a pair of words
    //! (string literal).
    //!
    //! See \ref
    //! contains_no_checks(Subclass&, Word const&, Word const&)
    //! for details.
    template <typename Subclass>
    [[nodiscard]] bool contains_no_checks(Subclass&   ci,
                                          char const* u,
                                          char const* v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.contains_no_checks(
          u, u + std::strlen(u), v, v + std::strlen(v));
    }

    //! \brief Helper for checking containment of a pair of words
    //! (std::initializer_list).
    //!
    //! See \ref
    //! contains_no_checks(Subclass&, Word const&, Word const&)
    //! for details.
    template <typename Subclass,
              typename Int = typename Subclass::native_letter_type>
    [[nodiscard]] bool contains_no_checks(Subclass&                         ci,
                                          std::initializer_list<Int> const& u,
                                          std::initializer_list<Int> const& v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return contains_no_checks<Subclass, std::initializer_list<Int>>(ci, u, v);
    }

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - contains
    ////////////////////////////////////////////////////////////////////////

    //! \brief Check containment of a pair of words.
    //!
    //! This function checks whether or not the words \p u and \p v are
    //! contained in the congruence represented by the instance \p ci of a
    //! derived class of \ref CongruenceInterface. This function triggers a
    //! full enumeration of \p ci, which may never terminate.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam Word the type of the second and third parameters.
    //!
    //! \param ci the subclass of \ref CongruenceInterface.
    //! \param u the left hand side of the pair to add.
    //! \param v the right hand side of the pair to add.
    //!
    //! \returns Whether or not the pair belongs to the congruence.
    //!
    //! \cong_intf_throws_if_letters_out_of_bounds
    template <typename Subclass, typename Word>
    [[nodiscard]] bool contains(Subclass& ci, Word const& u, Word const& v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.contains(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    // This version of the function catches the cases when u & v are not of the
    // same type but both convertible to string_view
    //! \brief Helper for checking containment of a pair of words
    //! (std::string_view).
    //!
    //! See \ref
    //! contains(Subclass&, Word const&, Word const&)
    //! for details.
    template <typename Subclass>
    [[nodiscard]] bool contains(Subclass&        ci,
                                std::string_view u,
                                std::string_view v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.contains(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    //! \brief Helper for checking containment of a pair of words
    //! (string literal).
    //!
    //! See \ref
    //! contains(Subclass&, Word const&, Word const&)
    //! for details.
    template <typename Subclass>
    [[nodiscard]] bool contains(Subclass& ci, char const* u, char const* v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return ci.contains(u, u + std::strlen(u), v, v + std::strlen(v));
    }

    //! \brief Helper for checking containment of a pair of words
    //! (std::initializer_list).
    //!
    //! See \ref
    //! contains_no_checks(Subclass&, Word const&, Word const&)
    //! for details.
    template <typename Subclass,
              typename Int = typename Subclass::native_letter_type>
    [[nodiscard]] bool contains(Subclass&                         ci,
                                std::initializer_list<Int> const& u,
                                std::initializer_list<Int> const& v) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return contains<Subclass, std::initializer_list<Int>>(ci, u, v);
    }
    //! @}

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - reduce_no_run_no_checks
    ////////////////////////////////////////////////////////////////////////

    //! \defgroup cong_intf_helpers_reduce_group Reduce a word
    //! \ingroup cong_intf_helpers_group
    //!
    //! \brief Check containment of a pair of words in a congruence.
    //!
    //! This page contains the documentation of the functions \ref
    //! reduce_no_run_no_checks; \ref reduce_no_run; \ref reduce_no_checks; and
    //! \ref reduce which can be invoked with a variety of different argument
    //! types.
    //!
    //! Functions with the suffix ``_no_run`` do not perform any enumeration
    //! of the \ref CongruenceInterface derived class instances; and those with
    //! the suffix ``_no_checks`` do not check that the input words are valid.
    //!
    //! @{

    //! \brief Reduce a word with no enumeration or checks.
    //!
    //! This function returns a reduced word equivalent to the input word \p w
    //! in the congruence represented by an instance of a derived class of \ref
    //! CongruenceInterface. This function triggers no enumeration. The word
    //! output by this function is equivalent to the input word in the
    //! congruence. If \p ci is ``finished``, then the output word is a normal
    //! form for the input word. If the \p ci is not ``finished``, then it
    //! might be that equivalent input words produce different output words.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam InputWord the type of the second parameter.
    //! \tparam OutputWord the type of word to be returned (defaults to \p
    //! InputWord).
    //!
    //! \param ci the subclass of \ref CongruenceInterface.
    //! \param w the word to reduce.
    //!
    //! \returns An irreducible word equivalent to \p w.
    //!
    //! \cong_intf_warn_assume_letters_in_bounds
    template <typename Subclass,
              typename InputWord,
              typename OutputWord = InputWord>
    [[nodiscard]] OutputWord reduce_no_run_no_checks(Subclass const&  ci,
                                                     InputWord const& w) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      OutputWord result;
      ci.reduce_no_run_no_checks(
          std::back_inserter(result), std::begin(w), std::end(w));
      return result;
    }

    // No string_view version is required because there is only a single "word"
    // parameter, and so the first template will catch every case except
    // initializer_list and char const*

    //! \brief Reduce a word (std::initializer_list).
    //!
    //! See \ref
    //! reduce_no_run_no_checks(Subclass const&, InputWord const&)
    //! for details.
    template <typename Subclass,
              typename Int        = size_t,
              typename OutputWord = std::vector<Int>>
    [[nodiscard]] OutputWord
    reduce_no_run_no_checks(Subclass const&                   ci,
                            std::initializer_list<Int> const& w) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      if constexpr (std::is_same_v<Subclass, ToddCoxeter>) {
        // We have a special case here for ToddCoxeter because ToddCoxeter
        // reverses the output words, in its reduce_no_run_no_checks, which does
        // not have the template parameter Subclass, so we call the 2-template
        // parameter overload of this function for ToddCoxeter. I.e. if we
        // called the 3-template param version as in the "else" case below, then
        // we'd just be calling the function above, which doesn't reverse the
        // words, and we'd get incorrect output.
        return reduce_no_run_no_checks<std::initializer_list<Int>, OutputWord>(
            ci, w);
      } else {
        return reduce_no_run_no_checks<Subclass,
                                       std::initializer_list<Int>,
                                       OutputWord>(ci, w);
      }
    }

    //! \brief Reduce a word (string literal).
    //!
    //! See \ref
    //! reduce_no_run_no_checks(Subclass const&, InputWord const&)
    //! for details.
    template <typename Subclass, typename OutputWord = std::string>
    [[nodiscard]] auto reduce_no_run_no_checks(Subclass const& ci,
                                               char const*     w) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      if constexpr (std::is_same_v<Subclass, ToddCoxeter>) {
        return reduce_no_run_no_checks<std::string_view, OutputWord>(ci, w);
      } else {
        return reduce_no_run_no_checks<Subclass, std::string_view, OutputWord>(
            ci, w);
      }
    }

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - reduce_no_run
    ////////////////////////////////////////////////////////////////////////

    //! \brief Reduce a word with no enumeration.
    //!
    //! This function returns a reduced word equivalent to the input word \p w
    //! in the congruence represented by an instance of a derived class of \ref
    //! CongruenceInterface. This function triggers no enumeration. The word
    //! output by this function is equivalent to the input word in the
    //! congruence. If \p ci is ``finished``, then the output word is a normal
    //! form for the input word. If the \p ci is not ``finished``, then it
    //! might be that equivalent input words produce different output words.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam InputWord the type of the second parameter.
    //! \tparam OutputWord the type of word to be returned (defaults to \p
    //! InputWord).
    //!
    //! \param ci the subclass of \ref CongruenceInterface.
    //! \param w the word to reduce.
    //!
    //! \returns An irreducible word equivalent to \p w.
    //!
    //! \cong_intf_throws_if_letters_out_of_bounds
    template <typename Subclass,
              typename InputWord,
              typename OutputWord = InputWord>
    [[nodiscard]] OutputWord reduce_no_run(Subclass const&  ci,
                                           InputWord const& w) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      OutputWord result;
      ci.reduce_no_run(std::back_inserter(result), std::begin(w), std::end(w));
      return result;
    }

    // No string_view version is required because there is only a single "word"
    // parameter, and so the first template will catch every case except
    // initializer_list and char const*

    //! \brief Reduce a word (std::initializer_list).
    //!
    //! See \ref
    //! reduce_no_run(Subclass const&, InputWord const&)
    //! for details.
    template <typename Subclass,
              typename Int        = size_t,
              typename OutputWord = std::vector<Int>>
    [[nodiscard]] OutputWord
    reduce_no_run(Subclass const& ci, std::initializer_list<Int> const& w) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      if constexpr (std::is_same_v<Subclass, ToddCoxeter>) {
        // See the comment above about why there's a special case here.
        return reduce_no_run<std::initializer_list<Int>, OutputWord>(ci, w);
      } else {
        return reduce_no_run<Subclass, std::initializer_list<Int>, OutputWord>(
            ci, w);
      }
    }

    //! \brief Reduce a word (string literal).
    //!
    //! See \ref
    //! reduce_no_run(Subclass const&, InputWord const&)
    //! for details.
    template <typename Subclass, typename OutputWord = std::string>
    [[nodiscard]] OutputWord reduce_no_run(Subclass const& ci, char const* w) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      if constexpr (std::is_same_v<Subclass, ToddCoxeter>) {
        // See the comment above about why there's a special case here.
        return reduce_no_run<std::string_view, OutputWord>(ci, w);
      } else {
        return reduce_no_run<Subclass, std::string_view, OutputWord>(ci, w);
      }
    }

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - reduce_no_checks
    ////////////////////////////////////////////////////////////////////////

    //! \brief Reduce a word with no checks.
    //!
    //! This function returns a reduced word equivalent to the input word \p w
    //! in the congruence represented by an instance of a derived class of \ref
    //! CongruenceInterface. This function triggers a full enumeration. The word
    //! output by this function is equivalent to the input word in the
    //! congruence. If \p ci is ``finished``, then the output word is a normal
    //! form for the input word. If the \p ci is not ``finished``, then it
    //! might be that equivalent input words produce different output words.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam InputWord the type of the second parameter.
    //! \tparam OutputWord the type of word to be returned (defaults to \p
    //! InputWord).
    //!
    //! \param ci the subclass of \ref CongruenceInterface.
    //! \param w the word to reduce.
    //!
    //! \returns An irreducible word equivalent to \p w.
    //!
    //! \cong_intf_warn_assume_letters_in_bounds
    template <typename Subclass,
              typename InputWord,
              typename OutputWord = InputWord>
    [[nodiscard]] OutputWord reduce_no_checks(Subclass&        ci,
                                              InputWord const& w) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      OutputWord result;
      ci.reduce_no_checks(
          std::back_inserter(result), std::begin(w), std::end(w));
      return result;
    }

    // No string_view version is required because there is only a single "word"
    // parameter, and so the first template will catch every case except
    // initializer_list and char const*

    //! \brief Reduce a word (std::initializer_list).
    //!
    //! See \ref
    //! reduce_no_checks(Subclass&, InputWord const&)
    //! for details.
    template <typename Subclass,
              typename Int        = size_t,
              typename OutputWord = std::vector<Int>>
    [[nodiscard]] OutputWord
    reduce_no_checks(Subclass& ci, std::initializer_list<Int> const& w) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      if constexpr (std::is_same_v<Subclass, ToddCoxeter>) {
        // See the comment above about why there's a special case here.
        return reduce_no_checks<std::initializer_list<Int>, OutputWord>(ci, w);
      } else {
        return reduce_no_checks<Subclass,
                                std::initializer_list<Int>,
                                OutputWord>(ci, w);
      }
    }

    //! \brief Reduce a word (string literal).
    //!
    //! See \ref
    //! reduce_no_checks(Subclass&, InputWord const&)
    //! for details.
    template <typename Subclass, typename OutputWord = std::string>
    [[nodiscard]] OutputWord reduce_no_checks(Subclass& ci, char const* w) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      if constexpr (std::is_same_v<Subclass, ToddCoxeter>) {
        // See the comment above about why there's a special case here.
        return reduce_no_checks<std::string_view, OutputWord>(ci, w);
      } else {
        return reduce_no_checks<Subclass, std::string_view, OutputWord>(ci, w);
      }
    }

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - reduce
    ////////////////////////////////////////////////////////////////////////

    //! \brief Reduce a word.
    //!
    //! This function returns a reduced word equivalent to the input word \p w
    //! in the congruence represented by an instance of a derived class of \ref
    //! CongruenceInterface. This function triggers a full enumeration. The word
    //! output by this function is equivalent to the input word in the
    //! congruence. If \p ci is ``finished``, then the output word is a normal
    //! form for the input word. If the \p ci is not ``finished``, then it
    //! might be that equivalent input words produce different output words.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam InputWord the type of the second parameter.
    //! \tparam OutputWord the type of word to be returned (defaults to \p
    //! InputWord).
    //!
    //! \param ci the subclass of \ref CongruenceInterface.
    //! \param w the word to reduce.
    //!
    //! \returns An irreducible word equivalent to \p w.
    //!
    //! \cong_intf_throws_if_letters_out_of_bounds
    template <typename Subclass,
              typename InputWord,
              typename OutputWord = InputWord>
    [[nodiscard]] OutputWord reduce(Subclass& ci, InputWord const& w) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      OutputWord result;
      ci.reduce(std::back_inserter(result), std::begin(w), std::end(w));
      return result;
    }

    // No string_view version is required because there is only a single "word"
    // parameter, and so the first template will catch every case except
    // initializer_list and char const*

    //! \brief Reduce a word (std::initializer_list).
    //!
    //! See \ref
    //! reduce(Subclass&, InputWord const&)
    //! for details.
    template <typename Subclass,
              typename Int        = size_t,
              typename OutputWord = std::vector<Int>>
    [[nodiscard]] OutputWord reduce(Subclass&                         ci,
                                    std::initializer_list<Int> const& w) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      if constexpr (std::is_same_v<Subclass, ToddCoxeter>) {
        // See the comment above about why there's a special case here.
        return reduce<std::initializer_list<Int>, OutputWord>(ci, w);
      } else {
        return reduce<Subclass, std::initializer_list<Int>, OutputWord>(ci, w);
      }
    }

    //! \brief Reduce a word (string literal).
    //!
    //! See \ref
    //! reduce(Subclass&, InputWord const&)
    //! for details.
    template <typename Subclass, typename OutputWord = std::string>
    [[nodiscard]] auto reduce(Subclass& ci, char const* w) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      if constexpr (std::is_same_v<Subclass, ToddCoxeter>) {
        return reduce<std::string_view, OutputWord>(ci, w);
      } else {
        return reduce<Subclass, std::string_view, OutputWord>(ci, w);
      }
    }

    //! @}

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - normal_forms
    ////////////////////////////////////////////////////////////////////////

    // There's nothing in common to implement in this file.

    //! \defgroup cong_intf_helpers_partition_group Partitioning
    //! \ingroup cong_intf_helpers_group
    //!
    //! \brief Partition a range of words by a congruence.
    //!
    //! This page contains the documentation of the functions \ref partition and
    //! \ref non_trivial_classes for partitioning a range of words by a
    //! congruence.
    //!
    //! @{

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - partition
    ////////////////////////////////////////////////////////////////////////

    //! \brief Partition a range of words.
    //!
    //! This function returns the partition of the words in the range \p r
    //! induced by the instance \p ci of a derived class of \ref
    //! CongruenceInterface. This function triggers a full enumeration of \p
    //! ci.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam Range the type of the input range of words.
    //! \tparam OutputWord the type of the words in the output (defaults to the
    //! type of the words in the input range).
    //!
    //! \param ci the derived class of \ref CongruenceInterface.
    //! \param r the input range of words.
    //!
    //! \returns The partition of the input range.
    //!
    //! \throws LibsemigroupsException if the input range of words is infinite.
    template <typename Subclass,
              typename Range,
              typename OutputWord = std::decay_t<typename Range::output_type>,
              typename = std::enable_if_t<rx::is_input_or_sink_v<Range>>>
    [[nodiscard]] std::vector<std::vector<OutputWord>> partition(Subclass& ci,
                                                                 Range     r) {
      // Congruence + ToddCoxeter have their own overloads for this
      static_assert(!std::is_same_v<Subclass, ToddCoxeter>
                    && !std::is_same_v<Subclass, Congruence>);

      if (!r.is_finite) {
        LIBSEMIGROUPS_EXCEPTION("the 2nd argument (a range) must be finite, "
                                "found an infinite range");
      }

      std::vector<std::vector<OutputWord>> result;

      std::unordered_map<OutputWord, size_t> map;
      size_t                                 index = 0;

      while (!r.at_end()) {
        auto next = r.get();
        if (ci.presentation().contains_empty_word() || !next.empty()) {
          auto next_nf        = congruence_interface::reduce(ci, next);
          auto [it, inserted] = map.emplace(next_nf, index);
          if (inserted) {
            result.emplace_back();
            index++;
          }
          size_t index_of_next_nf = it->second;
          result[index_of_next_nf].push_back(next);
        }
        r.next();
      }
      return result;
    }

    //! \brief Partition a range of words (via iterators)
    //!
    //! See \ref partition(Subclass&, Range) for details.
    template <
        typename Subclass,
        typename Iterator1,
        typename Iterator2,
        typename OutputWord = std::decay_t<
            typename rx::iterator_range<Iterator1, Iterator2>::output_type>>
    std::vector<std::vector<OutputWord>> partition(Subclass& ci,
                                                   Iterator1 first,
                                                   Iterator2 last) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      return partition<Subclass,
                       rx::iterator_range<Iterator1, Iterator2>,
                       OutputWord>(ci, rx::iterator_range(first, last));
    }

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - non_trivial_classes
    ////////////////////////////////////////////////////////////////////////

    //! \brief Find the non-trivial classes in the partition of a range of
    //! words.
    //!
    //! This function returns the classes with size at least \f$2\f$ in the
    //! partition of the words in the range \p r according to \p ci.
    //! This function triggers a full enumeration of \p ci.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam Range the type of the input range of words.
    //! \tparam OutputWord the type of the words in the output (defaults to the
    //! type of the words in the input range).
    //!
    //! \param ci a derived class of CongruenceInterface.
    //! \param r the input range of words.
    //!
    //! \returns The partition of the input range.
    //!
    //! \throws LibsemigroupsException if the input range of words is infinite.
    // couldn't get it to compile without copying
    template <typename Subclass,
              typename Range,
              typename OutputWord = std::decay_t<typename Range::output_type>,
              typename = std::enable_if_t<rx::is_input_or_sink_v<Range>>>
    std::vector<std::vector<OutputWord>> non_trivial_classes(Subclass& ci,
                                                             Range     r) {
      static_assert(std::is_base_of_v<CongruenceInterface, Subclass>);
      auto result = partition(ci, r);
      result.erase(
          std::remove_if(result.begin(),
                         result.end(),
                         [](auto const& x) -> bool { return x.size() <= 1; }),
          result.end());
      return result;
    }

    //! \brief Partition a range of words into non-trivial classes (via
    //! iterators).
    //!
    //! See \ref non_trivial_classes(Subclass&, Range) for details.
    template <
        typename Subclass,
        typename Iterator1,
        typename Iterator2,
        typename OutputWord = std::decay_t<
            typename rx::iterator_range<Iterator1, Iterator2>::output_type>>
    std::vector<std::vector<OutputWord>> non_trivial_classes(Subclass& ci,
                                                             Iterator1 first,
                                                             Iterator2 last) {
      return non_trivial_classes<Subclass,
                                 rx::iterator_range<Iterator1, Iterator2>,
                                 OutputWord>(ci,
                                             rx::iterator_range(first, last));
    }

    //! @}
  }  // namespace congruence_interface
}  // namespace libsemigroups
#endif  // LIBSEMIGROUPS_CONG_INTF_HPP_
