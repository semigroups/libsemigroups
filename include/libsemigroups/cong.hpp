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

// This file contains stuff for creating congruence over FroidurePin objects or
// over Presentation objects.

#ifndef LIBSEMIGROUPS_CONG_HPP_
#define LIBSEMIGROUPS_CONG_HPP_

#include <cstddef>  // for size_t
#include <memory>   // for shared_ptr

#include "cong-intf.hpp"        // for CongruenceInterface
#include "exception.hpp"        // for LIBSEMIGROUPS_EXCEPTION
#include "kambites.hpp"         // for Kambites
#include "knuth-bendix.hpp"     // for KnuthBendix
#include "to-todd-coxeter.hpp"  // for to_todd_coxeter
#include "todd-coxeter.hpp"     // for ToddCoxeter
#include "types.hpp"            // for word_type

#include "detail/race.hpp"  // for Race

namespace libsemigroups {
  class FroidurePinBase;  // Forward declaration, constructor parameters

  //! \ingroup cong_all_classes_group
  //!
  //! \brief Class for running Kambites, KnuthBendix, and \ref
  //! todd_coxeter_class_group "ToddCoxeter" in parallel.
  //!
  //! Defined in ``cong.hpp``.
  //!
  //! On this page we describe the functionality relating to the class
  //! Congruence in ``libsemigroups``. This class can be used for computing a
  //! congruence over a semigroup or monoid by running every applicable
  //! algorithm from ``libsemigroups`` (and some variants of the same algorithm)
  //! in parallel. This class is provided for convenience, at present it is not
  //! very customisable, and lacks some of the fine grained control offered by
  //! the classes implementing individual algorithms, such as Kambites,
  //! KnuthBendix, and \ref todd_coxeter_class_group "ToddCoxeter".
  //!
  //! \sa congruence_kind and tril.
  //!
  //! \par Example
  //! \code
  //! Presentation<word_type> p;
  //! p.alphabet(3)
  //! presentation::add_rule(p, {1, 2}, {0});
  //! Congruence cong(congruence_kind::twosided, p);
  //! is_obviously_infinite(cong);  // false
  //! congruence::add_generating_pair(cong, {1, 1, 1}, {0});
  //! cong.number_of_classes(); // 3
  //! \endcode
  // TODO(0):
  // * impl presentation mem fn
  class Congruence : public CongruenceInterface {
    enum class RunnerKind : size_t { TC = 0, KB = 1, K = 2 };

    /////////////////////////////////////////////////////////////////////////
    // Congruence - data - private
    /////////////////////////////////////////////////////////////////////////

    mutable detail::Race    _race;
    mutable bool            _runners_initted;
    std::vector<RunnerKind> _runner_kinds;

   public:
    //////////////////////////////////////////////////////////////////////////
    // Congruence - constructors - public
    //////////////////////////////////////////////////////////////////////////

    //! \brief Default constructor.
    //!
    //! This function default constructs an uninitialised Congruence instance.
    Congruence()
        : CongruenceInterface(), _race(), _runners_initted(), _runner_kinds() {
      init();
    }

    //! \brief Re-initialize a Congruence instance.
    //!
    //! This function puts a Congruence instance back into the state that it
    //! would have been in if it had just been newly default constructed.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    Congruence& init();

    //! Copy constructor.
    Congruence(Congruence const&) = default;

    //! Move constructor.
    Congruence(Congruence&&) = default;

    //! Copy assignment operator.
    Congruence& operator=(Congruence const&) = default;

    //! Move assignment operator.
    Congruence& operator=(Congruence&&) = default;

    ~Congruence() = default;

    //! \brief Construct from \ref congruence_kind and Presentation.
    //!
    //! This function constructs a Congruence instance representing a
    //! congruence of kind \p knd over the semigroup or monoid defined by the
    //! presentation \p p.
    //!
    //! \param knd the kind (onesided or twosided) of the congruence.
    //! \param p the presentation.
    //!
    //! \throws LibsemigroupsException if \p p is not valid.
    // No rvalue ref version because we anyway must copy p multiple times
    Congruence(congruence_kind knd, Presentation<word_type> const& p)
        : Congruence() {
      init(knd, p);
    }

    //! \brief Re-initialize a Congruence instance.
    //!
    //! This function puts a Congruence instance back into the state that it
    //! would have been in if it had just been newly constructed from \p knd and
    //! \p p.
    //!
    //! \param knd the kind (onesided or twosided) of the congruence.
    //! \param p the presentation.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \throws LibsemigroupsException if \p p is not valid.
    // No rvalue ref version because we anyway must copy p multiple times
    Congruence& init(congruence_kind knd, Presentation<word_type> const& p);

    //! \copydoc Congruence(congruence_kind, Presentation<word_type>&)
    // No rvalue ref version because we are not able to use p directly anyway
    template <typename Word>
    Congruence(congruence_kind knd, Presentation<Word> const& p)
        : Congruence(knd, to_presentation<word_type>(p, [](auto const& x) {
                       return x;
                     })) {}

    //! \copydoc init(congruence_kind, Presentation<word_type>&)
    // No rvalue ref version because we are not able to use p directly anyway
    template <typename Word>
    Congruence& init(congruence_kind knd, Presentation<Word> const& p) {
      init(knd, to_presentation<word_type>(p, [](auto const& x) { return x; }));
      return *this;
    }

    //! \brief Construct from congruence_kind, FroidurePin, and WordGraph.
    //!
    //! Constructs a Congruence over the FroidurePin instance \p S
    //! representing a 1- or 2-sided congruence according to \p knd.
    //!
    //! \tparam Node the type of the nodes in the 3rd argument \p wg (word
    //! graph).
    //!
    //! \param knd the kind (onesided or twosided) of the congruence.
    //! \param S a reference to the FroidurePin over which the congruence
    //! is being defined.
    //! \param wg the left or right Cayley graph of S.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning This constructor does not check its arguments.
    template <typename Node>
    Congruence(congruence_kind        knd,
               FroidurePinBase&       S,
               WordGraph<Node> const& wg)
        : Congruence() {
      init(knd, S, wg);
    }

    //! \brief Re-initialize from congruence_kind, FroidurePin, and WordGraph.
    //!
    //! This function re-initializes a Congruence instance as if it had been
    //! newly constructed over the FroidurePin instance \p S representing a 1-
    //! or 2-sided congruence according to \p knd.
    //!
    //! \tparam Node the type of the nodes in the 3rd argument \p wg (word
    //! graph).
    //!
    //! \param knd the kind (onesided or twosided) of the congruence.
    //! \param S a reference to the FroidurePin over which the congruence
    //! is being defined.
    //! \param wg the left or right Cayley graph of S.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning This constructor does not check its arguments.
    template <typename Node>
    Congruence& init(congruence_kind        knd,
                     FroidurePinBase&       S,
                     WordGraph<Node> const& wg);

    //////////////////////////////////////////////////////////////////////////
    // Congruence - interface requirements - add_generating_pair
    //////////////////////////////////////////////////////////////////////////

    //! \brief Add generating pair via iterators.
    //!
    //! This function adds a generating pair to the congruence represented by a
    //! \ref Congruence instance.
    //!
    //! \cong_intf_params_contains
    //!
    //! \returns A reference to `*this`.
    //!
    //! \cong_intf_warn_assume_letters_in_bounds
    //!
    //! \warning It is assumed that \ref started returns \c false. Adding
    //! generating pairs after \ref started is not permitted (but also not
    //! checked by this function).
    template <typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    Congruence& add_generating_pair_no_checks(Iterator1 first1,
                                              Iterator2 last1,
                                              Iterator3 first2,
                                              Iterator4 last2) {
      return CongruenceInterface::add_generating_pair_no_checks<Congruence>(
          first1, last1, first2, last2);
    }

    //! \brief Add generating pair via iterators.
    //!
    //! This function adds a generating pair to the congruence represented by a
    //! Congruence instance.
    //!
    //! \cong_intf_params_contains
    //!
    //! \returns A reference to `*this`.
    //!
    //! \cong_intf_throws_if_letters_out_of_bounds
    //!
    //! \cong_intf_throws_if_started
    template <typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    Congruence& add_generating_pair(Iterator1 first1,
                                    Iterator2 last1,
                                    Iterator3 first2,
                                    Iterator4 last2) {
      return CongruenceInterface::add_generating_pair<Congruence>(
          first1, last1, first2, last2);
    }

    ////////////////////////////////////////////////////////////////////////
    // Congruence - interface requirements - number_of_classes
    ////////////////////////////////////////////////////////////////////////

    //! \brief Compute the number of classes in the congruence.
    //!
    //! This function computes the number of classes in the congruence
    //! represented by a \ref todd_coxeter_class_group "ToddCoxeter" instance by
    //! running the congruence enumeration until it terminates.
    //!
    //! \returns The number of congruences classes of a \ref
    //! todd_coxeter_class_group "ToddCoxeter" instance if this number is
    //! finite, or \ref POSITIVE_INFINITY in some cases if this number is not
    //! finite.
    //!
    //! \warning Determining the number of classes is undecidable in general,
    //! and this may never terminate.
    [[nodiscard]] uint64_t number_of_classes();

    ////////////////////////////////////////////////////////////////////////
    // Congruence - interface requirements - contains
    ////////////////////////////////////////////////////////////////////////
   private:
    using CongruenceInterface::contains;
    using CongruenceInterface::currently_contains;

   public:
    template <typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    [[nodiscard]] tril currently_contains(Iterator1 first1,
                                          Iterator2 last1,
                                          Iterator3 first2,
                                          Iterator4 last2) const {
      return currently_contains<Congruence>(first1, last1, first2, last2);
    }

    template <typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    [[nodiscard]] bool contains(Iterator1 first1,
                                Iterator2 last1,
                                Iterator3 first2,
                                Iterator4 last2) {
      return contains<Congruence>(first1, last1, first2, last2);
    }

    //! \brief Check containment of a pair of words via iterators.
    //!
    //! This function checks whether or not the words represented by the ranges
    //! \p first1 to \p last1 and \p first2 to \p last2 are already known to be
    //! contained in the congruence represented by a \ref Congruence instance.
    //! This function performs no enumeration, so it is possible for the words
    //! to be contained in the congruence, but that this is not currently known.
    //!
    //! \cong_intf_params_contains
    //!
    //! \returns
    //! * tril::TRUE if the words are known to belong to the congruence;
    //! * tril::FALSE if the words are known to not belong to the congruence;
    //! * tril::unknown otherwise.
    //!
    //! \cong_intf_warn_assume_letters_in_bounds
    template <typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    [[nodiscard]] tril currently_contains_no_checks(Iterator1 first1,
                                                    Iterator2 last1,
                                                    Iterator3 first2,
                                                    Iterator4 last2) const;

    //! \brief Check containment of a pair of words via iterators.
    //!
    //! This function checks whether or not the words represented by the ranges
    //! \p first1 to \p last1 and \p first2 to \p last2 are already known to be
    //! contained in the congruence represented by a \ref Congruence instance.
    //! This function performs no enumeration, so it is possible for the words
    //! to be contained in the congruence, but that this is not currently known.
    //!
    //! \cong_intf_params_contains
    //!
    //! \returns
    //! * tril::TRUE if the words are known to belong to the congruence;
    //! * tril::FALSE if the words are known to not belong to the congruence;
    //! * tril::unknown otherwise.
    //!
    //! \cong_intf_throws_if_letters_out_of_bounds
    template <typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    [[nodiscard]] bool contains_no_checks(Iterator1 first1,
                                          Iterator2 last1,
                                          Iterator3 first2,
                                          Iterator4 last2);

    ////////////////////////////////////////////////////////////////////////
    // OLD
    ////////////////////////////////////////////////////////////////////////

    // TODO(0) replace with iterators
    // TODO(0) out of line
    [[nodiscard]] bool contains(word_type const& u, word_type const& v) {
      return contains(std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    // TODO(0) out of line
    template <typename Iterator1, typename Iterator2>
    void throw_if_letter_out_of_bounds(Iterator1 first, Iterator2 last) const {
      if (!_race.empty()) {
        size_t index = (finished() ? _race.winner_index() : 0);

        if (_runner_kinds[index] == RunnerKind::TC) {
          std::static_pointer_cast<ToddCoxeter>(*_race.begin())
              ->throw_if_letter_out_of_bounds(first, last);
        } else if (_runner_kinds[index] == RunnerKind::KB) {
          std::static_pointer_cast<KnuthBendix<>>(*_race.begin())
              ->throw_if_letter_out_of_bounds(first, last);
        } else {
          LIBSEMIGROUPS_ASSERT(_runner_kinds[index] == RunnerKind::K);
          std::static_pointer_cast<Kambites<word_type>>(*_race.begin())
              ->throw_if_letter_out_of_bounds(first, last);
        }
        return;
      }
      LIBSEMIGROUPS_EXCEPTION(
          "No presentation has been set, so cannot validate the word!");
    }

    //////////////////////////////////////////////////////////////////////////
    // Congruence - member functions - public
    //////////////////////////////////////////////////////////////////////////

    template <typename Thing>
    std::shared_ptr<Thing> get() {
      init_runners();
      auto result = _race.find_runner<Thing>();
      if (result == nullptr) {
        LIBSEMIGROUPS_EXCEPTION("There's no Thing to get!");
      }
      return result;
    }

    template <typename Thing>
    [[nodiscard]] bool has() {
      return _race.find_runner<Thing>() != nullptr;
    }

    //! Returns the KnuthBendix instance used to compute the congruence (if
    //! any).
    //!
    //! \returns A std::shared_ptr to a KnuthBendix or \c nullptr.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    //!
    //! \sa has<KnuthBendix>.

    //! Checks if a KnuthBendix instance is being used to compute
    //! the congruence.
    //!
    //! \returns A value of type `bool`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    //!
    //! \sa get<KnuthBendix>.

    //! Get the current maximum number of threads.
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] size_t max_threads() const noexcept {
      return _race.max_threads();
    }

    //! Set the maximum number of threads.
    //!
    //! \param val the number of threads.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    Congruence& max_threads(size_t val) noexcept {
      _race.max_threads(val);
      return *this;
    }

    [[nodiscard]] size_t number_of_runners() const noexcept {
      return _race.number_of_runners();
    }

   private:
    void add_runner(std::shared_ptr<ToddCoxeter> ptr) {
      _race.add_runner(ptr);  // TODO move?
      _runner_kinds.push_back(RunnerKind::TC);
    }

    void add_runner(std::shared_ptr<KnuthBendix<>> ptr) {
      _race.add_runner(ptr);  // TODO move?
      _runner_kinds.push_back(RunnerKind::KB);
    }

    void add_runner(std::shared_ptr<Kambites<word_type>> ptr) {
      _race.add_runner(ptr);  // TODO move?
      _runner_kinds.push_back(RunnerKind::K);
    }

    //////////////////////////////////////////////////////////////////////////
    // Congruence - member functions - private
    //////////////////////////////////////////////////////////////////////////

    void init_runners() const;

    //////////////////////////////////////////////////////////////////////////
    // Runner - pure virtual member functions - private
    //////////////////////////////////////////////////////////////////////////

    void run_impl() override;
    bool finished_impl() const override {
      return _race.finished();
    }
  };  // class Congruence

  ////////////////////////////////////////////////////////////////////////
  // Congruence - out of line implementations
  ////////////////////////////////////////////////////////////////////////

  template <typename Node>
  Congruence& Congruence::init(congruence_kind        knd,
                               FroidurePinBase&       S,
                               WordGraph<Node> const& wg) {
    if (S.is_finite() != tril::FALSE) {
      S.run();
    } else {
      LIBSEMIGROUPS_EXCEPTION(
          "the 2nd argument does not represent a finite semigroup!");
    }
    CongruenceInterface::init(knd);

    // TODO(later) if necessary make a runner that tries to S.run(), then get
    // the Cayley graph and use that in the ToddCoxeter, at present that'll
    // happen here in the constructor, same for the creation of the
    // presentation this could take place in the Runner so that they are done
    // in parallel
    add_runner(std::make_shared<ToddCoxeter>(to_todd_coxeter(knd, S, wg)));

    auto tc = to_todd_coxeter(knd, S, wg);
    tc.strategy(ToddCoxeter::options::strategy::felsch);
    add_runner(std::make_shared<ToddCoxeter>(std::move(tc)));

    tc = ToddCoxeter(knd, to_presentation<word_type>(S));
    add_runner(std::make_shared<ToddCoxeter>(std::move(tc)));

    tc = ToddCoxeter(knd, to_presentation<word_type>(S));
    tc.strategy(ToddCoxeter::options::strategy::felsch);
    add_runner(std::make_shared<ToddCoxeter>(std::move(tc)));

    return *this;
  }

  template <typename Iterator1,
            typename Iterator2,
            typename Iterator3,
            typename Iterator4>
  [[nodiscard]] tril
  Congruence::currently_contains_no_checks(Iterator1 first1,
                                           Iterator2 last1,
                                           Iterator3 first2,
                                           Iterator4 last2) const {
    if (finished()) {
      auto winner_kind = _runner_kinds[_race.winner_index()];
      if (winner_kind == RunnerKind::TC) {
        return std::static_pointer_cast<ToddCoxeter>(_race.winner())
                       ->contains_no_checks(first1, last1, first2, last2)
                   ? tril::TRUE
                   : tril::FALSE;
      } else if (winner_kind == RunnerKind::KB) {
        return std::static_pointer_cast<KnuthBendix<>>(_race.winner())
                       ->contains_no_checks(first1, last1, first2, last2)
                   ? tril::TRUE
                   : tril::FALSE;
      } else {
        LIBSEMIGROUPS_ASSERT(winner_kind == RunnerKind::K);
        return std::static_pointer_cast<Kambites<word_type>>(_race.winner())
                       ->contains_no_checks(first1, last1, first2, last2)
                   ? tril::TRUE
                   : tril::FALSE;
      }
    }
    init_runners();
    tril result = tril::unknown;
    for (auto const& [i, runner] : rx::enumerate(_race)) {
      if (_runner_kinds[i] == RunnerKind::TC) {
        result
            = std::static_pointer_cast<ToddCoxeter>(runner)
                  ->currently_contains_no_checks(first1, last1, first2, last2);
      } else if (_runner_kinds[i] == RunnerKind::KB) {
        result
            = std::static_pointer_cast<KnuthBendix<>>(runner)
                  ->currently_contains_no_checks(first1, last1, first2, last2);
      } else {
        LIBSEMIGROUPS_ASSERT(_runner_kinds[i] == RunnerKind::K);
        result
            = std::static_pointer_cast<Kambites<word_type>>(runner)
                  ->currently_contains_no_checks(first1, last1, first2, last2);
      }
      if (result != tril::unknown) {
        break;
      }
    }
    return result;
  }

  template <typename Iterator1,
            typename Iterator2,
            typename Iterator3,
            typename Iterator4>
  [[nodiscard]] bool Congruence::contains_no_checks(Iterator1 first1,
                                                    Iterator2 last1,
                                                    Iterator3 first2,
                                                    Iterator4 last2) {
    run();
    LIBSEMIGROUPS_ASSERT(_race.winner_index() != UNDEFINED);
    auto winner_kind = _runner_kinds[_race.winner_index()];
    if (winner_kind == RunnerKind::TC) {
      return std::static_pointer_cast<ToddCoxeter>(_race.winner())
          ->contains_no_checks(first1, last1, first2, last2);
    } else if (winner_kind == RunnerKind::KB) {
      return std::static_pointer_cast<KnuthBendix<>>(_race.winner())
          ->contains_no_checks(first1, last1, first2, last2);
    } else {
      LIBSEMIGROUPS_ASSERT(winner_kind == RunnerKind::K);
      return std::static_pointer_cast<Kambites<word_type>>(_race.winner())
          ->contains_no_checks(first1, last1, first2, last2);
    }
  }

  // There's no doc for anything in this section it's covered by the congruence
  // interface helpers.
  namespace congruence {

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - add_generating_pair
    ////////////////////////////////////////////////////////////////////////

    using congruence_interface::add_generating_pair;
    using congruence_interface::add_generating_pair_no_checks;

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - contains
    ////////////////////////////////////////////////////////////////////////

    using congruence_interface::contains;
    using congruence_interface::contains_no_checks;
    using congruence_interface::currently_contains;
    using congruence_interface::currently_contains_no_checks;

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - reduce
    ////////////////////////////////////////////////////////////////////////

    using congruence_interface::reduce;
    using congruence_interface::reduce_no_checks;
    using congruence_interface::reduce_no_run;
    using congruence_interface::reduce_no_run_no_checks;

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - normal_forms
    ////////////////////////////////////////////////////////////////////////

    // The following doesn't work, because the types of the two returned values
    // aren't the same.
    // TODO(0) implement a class containing a variant for this.
    // template <typename Word>
    // auto normal_forms(Congruence& cong) {
    //   cong.run();
    //   if (cong.has<ToddCoxeter>() && cong.get<ToddCoxeter>()->finished()) {
    //     return todd_coxeter::normal_forms(*cong.get<ToddCoxeter>());
    //   } else if (cong.has<KnuthBendix<>>()
    //              && cong.get<KnuthBendix<>>()->finished()) {
    //     return knuth_bendix::normal_forms(*cong.get<KnuthBendix<>>());
    //   }
    //   // There's currently no normal_forms function for Kambites, so can't
    //   // return anything in that case.
    //   LIBSEMIGROUPS_EXCEPTION("Cannot compute the non-trivial classes!");
    // }

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - partition
    ////////////////////////////////////////////////////////////////////////

    template <typename Range>
    std::vector<std::vector<std::decay_t<typename Range::output_type>>>
    partition(Congruence& cong, Range r) {
      cong.run();
      if (cong.has<ToddCoxeter>() && cong.get<ToddCoxeter>()->finished()) {
        return todd_coxeter::partition(*cong.get<ToddCoxeter>(), r);
      } else if (cong.has<KnuthBendix<>>()
                 && cong.get<KnuthBendix<>>()->finished()) {
        return knuth_bendix::partition(*cong.get<KnuthBendix<>>(), r);
      } else if (cong.has<Kambites<word_type>>()
                 && cong.get<Kambites<word_type>>()->success()) {
        return kambites::partition(*cong.get<Kambites<word_type>>(), r);
      }
      LIBSEMIGROUPS_EXCEPTION("Cannot compute the non-trivial classes!");
    }

    ////////////////////////////////////////////////////////////////////////
    // Interface helpers - non_trivial_classes
    ////////////////////////////////////////////////////////////////////////

    using congruence_interface::non_trivial_classes;

    // This is a copy of the function in the congruence_interface namespace,
    // couldn't get it to compile without copying (probably just require some
    // SFINAE)
    template <typename Range,
              typename Word = std::decay_t<typename Range::output_type>,
              typename      = std::enable_if_t<rx::is_input_or_sink_v<Range>>>
    std::vector<std::vector<Word>> non_trivial_classes(Congruence& ci,
                                                       Range       r) {
      auto result = partition(ci, r);
      result.erase(
          std::remove_if(result.begin(),
                         result.end(),
                         [](auto const& x) -> bool { return x.size() <= 1; }),
          result.end());
      return result;
    }

    // template <typename Word>
    // std::vector<std::vector<Word>> non_trivial_classes(Congruence& cong1,
    //                                                   Congruence& cong2);
    // TODO(1) can't really do this until normal_forms is working.

    // Not required.
    // We have to pass the presentation here, because o/w we have no way of
    // knowing over what we should compute the non-trivial classes (i.e. we
    // cannot always recover p from cong).
    // TODO remove?
    std::vector<std::vector<word_type>>
    non_trivial_classes(Congruence& cong, Presentation<word_type> const& p);

    // TODO remove?
    std::vector<std::vector<std::string>>
    non_trivial_classes(Congruence& cong, Presentation<std::string> const& p);

  }  // namespace congruence
}  // namespace libsemigroups

#endif  // LIBSEMIGROUPS_CONG_HPP_
