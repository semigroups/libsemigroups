//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2022-24 James D. Mitchell
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

// This file contains a declaration of a class for performing the "low-index
// congruence" algorithm for 1-sided or 2-sided congruences of semigroups and
// monoids.

// TODO(0):
// * review the function aliases, and remove them if they are unnecessary
// * const
// * noexcept
// * nodiscard
// * python bindings
// * gap bindings

// TODO(2) (later):
// * use SimsRefinerFaithful in RepOrc + MinimalRepOrc
// * a version which allows specifying the word_graph to Sims1 too
// * implement maximum_2_sided_congruence_contained to compute the kernel of the
//   associated homomorphism, which is the largest 2-sided congruence contained
//   in the right congruence.
// * change RepOrc and MinimalRepOrc to compute minimal 2-sided congruences
//   first (by using generating pairs), and then to try and find a right
//   congruence not containing any of the minimal 2-sided congruences.
//
// TODO(3) (far future):
// * Figure out a way of making refining functions possible. A refining function
// differs from a pruner in that it is also allowed to modify the word graph.
// * Implement one and two sided compatibility checking as refiners, following
// the theoretical description in the paper.
// * Modify SimsBase to perform the backtrack on all standard word graphs
// (without explicitly checking compatibility).
// * Make SimsBase -> Sims, i.e. implement all the relevant functions and
// settings of Sims1, Sims2, RepOrc, MinimalRepOrc into a single unified
// class.
// * Implement Sims1, Sim2, RepOrc, MinimalRepOrc just particular constructors
// of the Sims object, which automatically set the relevant pruners etc.

#ifndef LIBSEMIGROUPS_SIMS_HPP_
#define LIBSEMIGROUPS_SIMS_HPP_

#include <algorithm>    // for max, fill, count, for_each
#include <atomic>       // for atomic_uint64_t, __atomic_base
#include <chrono>       // for operator-
#include <cstddef>      // for size_t
#include <cstdint>      // for uint64_t, uint32_t
#include <functional>   // for function
#include <iterator>     // for distance, forward_iterator_tag, pair
#include <list>         // for operator!=
#include <memory>       // for unique_ptr
#include <mutex>        // for mutex
#include <string>       // for operator+, basic_string, string, to_string
#include <thread>       // for thread
#include <type_traits>  // for decay_t
#include <utility>      // for move, pair, swap
#include <vector>       // for vector

#include "constants.hpp"     // for operator!=, operator==, Max, UNDEFINED
#include "debug.hpp"         // for LIBSEMIGROUPS_ASSERT
#include "dot.hpp"           // for Dot
#include "exception.hpp"     // for LIBSEMIGROUPS_EXCEPTION
#include "forest.hpp"        // for Forest
#include "knuth-bendix.hpp"  // for KnuthBendix, to_presentation
#include "matrix.hpp"        // for DynamicMatrix, MatrixCommon, BMat
#include "presentation.hpp"  // for Presentation, longest_rule_length, short...
#include "ranges.hpp"        // for operator|, iterator_range
#include "runner.hpp"        // for delta, Reporter
#include "todd-coxeter.hpp"  // for ToddCoxeter
#include "types.hpp"         // for word_type, relation_type, congruence_kind
#include "word-graph.hpp"    // for Joiner, WordGraph

#include "detail/felsch-graph.hpp"    // for FelschGraph
#include "detail/fmt.hpp"             // for format, print
#include "detail/iterator.hpp"        // for detail/default_postfix_increment
#include "detail/path-iterators.hpp"  // for const_pstilo_iterator
#include "detail/report.hpp"          // for ReportCell
#include "detail/rewriters.hpp"       // for RewriteTrie
#include "detail/word-graph-with-sources.hpp"  // for WordGraphWithSources

//! \defgroup congruences_group Congruences
//!
//! This file contains documentation for the functionality in
//! `libsemigroups` related to congruences of semigroups and monoids.

namespace libsemigroups {

  //! \ingroup congruences_group
  //!
  //! \brief For keeping track of various statistics arising during the runtime
  //! of the low index algorithm.
  //!
  //! Defined in `sims.hpp`.
  //!
  //! On this page we describe the SimsStats class. The purpose of this
  //! class is to collect some statistics related to Sims1 or Sims2 class
  //! template.
  //!
  //! \sa \ref Sims1
  //! \sa \ref Sims2
  class SimsStats {
   public:
    //! Number of congruences found at time of last report.
    //!
    //! This member holds the number of congruences found by the
    //! Sims1 or Sims2 algorithm at the time of the last call
    //! to \ref stats_check_point.
    //!
    //! \sa \ref stats_check_point
    //! \sa \ref count_now
    // TODO(2) might be better to have a mutex here and just lock it in
    // check_point below
    std::atomic_uint64_t count_last;

    //! Number of congruences found up to this point.
    //!
    //! This member tracks the total number of congruences found during the
    //! running of the Sims1 or Sims2 algorithm.
    //! \sa \ref count_last
    // Atomic so as to avoid races between report_progress_from_thread and the
    // threads modifying count_last
    std::atomic_uint64_t count_now;

    //! The maximum number of pending definitions.
    //!
    //! A *pending definition* is just an edge that will be defined at some
    //! point in the future in the WordGraph represented by a Sims1 or Sims2
    //! instance at any given moment.
    //!
    //! This member tracks the maximum number of such pending definitions that
    //! occur during the running of the algorithms in Sims1 or Sims2.
    std::atomic_uint64_t max_pending;

    //! The total number of pending definitions at time of last report.
    //!
    //! A *pending definition* is just an edge that will be defined at some
    //! point in the future in the WordGraph represented by a Sims1 or Sims2
    //! instance at any given moment.
    //!
    //! This member tracks the total number of pending definitions
    //! that occur at the time of the last call to \ref stats_check_point.
    //! This is the same as the number of nodes in the search
    //! tree encounter during the running of Sims1 or Sims2.
    //!
    //! \sa \ref stats_check_point
    //! \sa \ref total_pending_now
    std::atomic_uint64_t total_pending_last;

    //! The total number of pending definitions.
    //!
    //! A *pending definition* is just an edge that will be defined at some
    //! point in the future in the WordGraph represented by a Sims1 or Sims2
    //! instance at any given moment.
    //!
    //! This member tracks the total number of pending definitions that
    //! occur during the running of the algorithms in Sims1 or Sims2. This is
    //! the same as the number of nodes in the search tree encounter during the
    //! running of Sims1 or Sims2.
    //!
    //! \sa \ref total_pending_last
    // Atomic so as to avoid races between report_progress_from_thread and the
    // threads modifying total_pending_now
    std::atomic_uint64_t total_pending_now;

    //! Default constructor.
    //!
    //! Constructs a SimsStats object with all statistics set to zero.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    SimsStats();

    //! \brief Reinitialize an existing SimsStats object.
    //!
    //! This function puts a SimsStats object back into the same state as if
    //! it had been newly default constructed.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    SimsStats& init() {
      stats_zero();
      return *this;
    }

    //! Copy constructor.
    //!
    //! Returns a SimsStats object that is a copy of \p that. The state
    //! of the new SimsStats object is the same as \p that. This triggers an
    //! atomic load on the member variables of \p that.
    //!
    //! \param that the SimsStats to copy.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    SimsStats(SimsStats const& that) : SimsStats() {
      init(that);
    }

    //! Copy assignment operator.
    //!
    //! Assigns a SimsStats object that is a copy of \p that. The state
    //! of the new SimsStats object is the same as \p that. This triggers an
    //! atomic load on the member variables of \p that.
    //!
    //! \param that the SimsStats to move from.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    SimsStats& operator=(SimsStats const& that) {
      return init(that);
    }

    //! Move constructor.
    //!
    //! Returns a SimsStats object that is initialized from \p that. The state
    //! of the new SimsStats object is the same as \p that. This triggers an
    //! atomic load on the member variables of \p that.
    //!
    //! \param that the SimsStats to move from.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    SimsStats(SimsStats&& that) : SimsStats() {
      init(that);
    }

    //! Move assignment operator.
    //!
    //! Assigns a SimsStats object that is initialized from \p that. The state
    //! of the new SimsStats object is the same as \p that. This triggers an
    //! atomic load on the member variables of \p that.
    //!
    //! \param that the SimsStats to move from.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    SimsStats& operator=(SimsStats&& that) {
      return init(that);
    }

    //! Initialize from SimsStats.
    SimsStats& init(SimsStats const& that);

    //! \brief Set all statistics to zero.
    //!
    //! Set all statistics to zero.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    SimsStats& stats_zero();

    //! \brief Store the current statistic values.
    //!
    //! Overwrites the values of \ref count_last and \ref total_pending_last
    //! with \ref count_now and \ref total_pending_now respectively. Triggers an
    //! atomic load on \ref count_now and \ref total_pending_now.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    SimsStats& stats_check_point() {
      count_last         = count_now.load();
      total_pending_last = total_pending_now.load();
      return *this;
    }
  };

  //! \ingroup congruences_group
  //!
  //! \brief For setting the presentation and various runtime parameters of the
  //! Sims low index algorithm.
  //!
  //! Defined in `sims.hpp`.
  //!
  //! On this page we describe the SimsSettings class. The purpose of this class
  //! is to allow us to use the same interface for settings for Sims1, Sims2,
  //! RepOrc, and MinimalRepOrc.
  //!
  //! \sa \ref Sims1
  //! \sa \ref Sims2
  //! \sa \ref RepOrc
  //! \sa \ref MinimalRepOrc
  template <typename Subclass>
  class SimsSettings {
   public:
    // We use WordGraph, even though the iterators produced by this class
    // hold FelschGraph's, none of the features of FelschGraph are useful
    // for the output, only for the implementation
    //! The type of the associated WordGraph objects.
    using word_graph_type = WordGraph<uint32_t>;

    //! Type for the nodes in the associated WordGraph objects.
    using node_type = typename word_graph_type::node_type;

    //! Type for the labels in the associated WordGraph objects.
    using label_type = typename word_graph_type::label_type;

    //! The size_type of the associated WordGraph objects.
    using size_type = typename word_graph_type::size_type;

    //! Type for letters in the underlying presentation.
    using letter_type = typename word_type::value_type;

   private:
    std::vector<word_type>                 _exclude;
    size_t                                 _exclude_pruner_index;
    size_t                                 _idle_thread_restarts;
    std::vector<word_type>                 _include;
    std::vector<word_type>::const_iterator _longs_begin;
    size_t                                 _num_threads;
    Presentation<word_type>                _presentation;
    std::vector<std::function<bool(word_graph_type const&)>> _pruners;
    mutable SimsStats                                        _stats;

   public:
    //! Default constructor.
    //!
    //! Constructs a SimsSettings object. Defaults to a single thread and 64
    //! idle thread restarts, no other settings set. Use
    //! * \ref presentation to set the presentation;
    //! * \ref number_of_threads to set the number of threads;
    //! * \ref included_pairs to set the pairs to be included;
    //! * \ref excluded_pairs to set the pairs to be excluded;
    //! * \ref add_pruner to add a pruner;
    //! * \ref long_rule_length to set the length of long rules;
    //! * \ref idle_thread_restarts to set the number of idle thread restarts.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    SimsSettings();

    //! \brief Reinitialize an existing SimsSettings object.
    //!
    //! This function puts a SimsSettings object back into the same state as if
    //! it had been newly default constructed.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    Subclass& init();

    //! Copy constructor.
    //!
    //! Returns a SimsSettings object that is a copy of \p that. The state
    //! of the new SimsSettings object is the same as \p that.
    //!
    //! \param that the SimsSettings to copy.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    // Copy constructor is explicitly required, the constructor template is not
    // a substitute. If no copy constructor is implemented, then _longs_begin
    // is not properly initialised, and leads to badness.
    SimsSettings(SimsSettings const& that) {
      init(that);
    }

    //! Copy assignment operator.
    //!
    //! Assigns a SimsSettings object that is a copy of \p that. The state
    //! of the new SimsSettings object is the same as \p that.
    //!
    //! \param that the SimsSettings to copy.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    SimsSettings& operator=(SimsSettings const& that) {
      init(that);
      return *this;
    }

    //! Move constructor.
    //!
    //! Returns a SimsSettings object that is initialized from \p that. The
    //! state of the new SimsSettings object is the same as \p that.
    //!
    //! \param that the SimsSettings to copy.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    SimsSettings(SimsSettings&& that) {
      init(that);
    }

    //! Move assignment operator.
    //!
    //! Assigns a SimsSettings object that is initialized from \p that. The
    //! state of the new SimsSettings object is the same as \p that.
    //!
    //! \param that the SimsSettings to copy.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    SimsSettings& operator=(SimsSettings&& that) {
      init(that);
      return *this;
    }

    //! Construct from SimsSettings with different subclass.
    template <typename OtherSubclass>
    SimsSettings(SimsSettings<OtherSubclass> const& that) {
      init(that);
    }

    //! Initialize from SimsSettings with different subclass.
    template <typename OtherSubclass>
    SimsSettings& init(SimsSettings<OtherSubclass> const& that);

    ~SimsSettings();

    //! Returns the settings object of *this.
    //!
    //! The settings object contains all the settings that are common to
    //! Sims1, Sims2, RepOrc, and MinimalRepOrc, which are currently:
    //! * \ref presentation
    //! * \ref long_rule_length
    //! * \ref number_of_threads
    //!
    //! The return value of this function can be used to initialise another
    //! Sims1, Sims2, RepOrc, or MinimalRepOrc with these settings.
    //!
    //! \returns A const reference to SimsSettings.
    //!
    //! \exceptions
    //! \noexcept
    // So that we can access the settings from the derived class T.
    [[nodiscard]] SimsSettings const& settings() const noexcept {
      return *this;
    }

    //! \brief Set the number of threads.
    //!
    //! This function sets the number of threads to be used by Sims1 or Sims2.
    //!
    //! The default value is `1`.
    //!
    //! \param val the maximum number of threads to use.
    //!
    //! \returns A reference to \c this.
    //!
    //! \throws LibsemigroupsException if the argument \p val is 0.
    //!
    //! \note The value of \p val is capped at
    //! `std::thread::hardware_concurrency()`. Trying to set a higher value is
    //! equivalent to setting `std::thread::hardware_concurrency()`.
    Subclass& number_of_threads(size_t val);

    //! \brief Get the number of threads.
    //!
    //! Returns the current number of threads.
    //!
    //! \returns
    //! A `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] size_t number_of_threads() const noexcept {
      return _num_threads;
    }

    //! \brief Set the presentation over which the congruences produced by an
    //! instance are defined.
    //!
    //! \anchor presentation
    //! Set the presentation over which the congruences produced by an
    //! instance are defined.
    //! These are the rules used at every node in the depth first search
    //! conducted by objects of this type.
    //!
    //! If the template parameter \p Word is not `word_type`, then
    //! the parameter \p p is first converted to a value of type
    //! `Presentation<word_type>` and it is this converted value that is used.
    //!
    //! \tparam Word the type of words in the input presentation.
    //!
    //! \param p the presentation.
    //!
    //! \returns A reference to \c this.
    //!
    //! \throws LibsemigroupsException if `to_presentation<word_type>(p)`
    //! throws.
    //!
    //! \throws LibsemigroupsException if `p` is not valid.
    //!
    //! \throws LibsemigroupsException if the alphabet of `p` is non-empty and
    //! not compatible with \ref included_pairs or \ref excluded_pairs.
    //!
    //! \throws LibsemigroupsException if `p` has 0-generators and 0-relations.
    template <typename Word>
    Subclass& presentation(Presentation<Word> const& p);

    //! \brief Get the presentation over which the congruences produced by an
    //! instance are defined.
    //!
    //! \anchor presentation
    //! Returns a const reference to the current relations of the underlying
    //! presentation.
    //!
    //! This function returns the defining presentation of a Sims1 or Sims2
    //! instance. The congruences computed by `cbegin` and `cend` of the
    //! appropriate subclass are defined over the semigroup or monoid defined by
    //! this presentation.
    //!
    //! \returns A const reference to `Presentation<word_type>`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \warning
    //! The value returned by this function is not guaranteed be the same as
    //! the presentation that was used to construct the object! A Sims1 or Sims2
    //! object requires the generators of the defining presentation
    //! \f$\mathcal{P}\f$ to be \f$\{0, \ldots, n - 1\}\f$ where \f$n\f$ is the
    //! size of the alphabet of \f$\mathcal{P}\f$. Every occurrence of every
    //! generator \c a in the presentation \c p used to construct a Sims1 or
    //! Sims2 instance is replaced by `p.index(a)`.
    [[nodiscard]] Presentation<word_type> const& presentation() const noexcept {
      return _presentation;
    }

    //! \brief Set the beginning of the long rules (iterator).
    //!
    //! Set the beginning of the long rules using an iterator.
    //!
    //! The "long rules" are the rules used after a complete deterministic word
    //! graph has been found in the search. If such a word graph is compatible
    //! with the long rules specified by this function, then this word graph is
    //! accepted, and if not it is rejected.
    //!
    //! The purpose of this is to improve the backtrack search by reducing the
    //! time spent processing "long" rules in each node of the search tree, and
    //! to only check them at the leaves.
    //!
    //! \param it iterator pointing at the left hand side of the first long
    //! rule.
    //!
    //! \returns A reference to \c this.
    //!
    //! \throws LibsemigroupsException if \p it does not point into
    //! `presentation().rules`
    //!
    //! \throws LibsemigroupsException if  \p it does not point at the left hand
    //! side of a rule (i.e. if it points at an odd indexed position in
    //! `presentation().rules`).
    Subclass& cbegin_long_rules(std::vector<word_type>::const_iterator it);

    //! \brief Set the beginning of the long rules (position).
    //!
    //! Set the beginning of the long rules using a position in
    //! `presentation().rules`.
    //!
    //! The "long rules" are the rules used after a complete deterministic word
    //! graph has been found in the search. If such a word graph is compatible
    //! with the long rules specified by this function, then this word graph is
    //! accepted, and if not it is rejected.
    //!
    //! The purpose of this is to improve the backtrack search by reducing the
    //! time spent processing "long" rules in each node of the search tree, and
    //! to only check them at the leaves.
    //!
    //! \param pos position of the the left hand side of the first long
    //! rule.
    //!
    //! \returns A reference to \c this.
    //!
    //! \throws LibsemigroupsException if \p pos is not a valid position in
    //! `presentation().rules`
    //!
    //! \throws LibsemigroupsException if the rule at position \p pos is not
    //! the left hand side of a rule (i.e. if \p pos is odd).
    Subclass& cbegin_long_rules(size_t pos);

    //! \brief Get the pointer to the first long rule.
    //!
    //! Returns the pointer to the first long rule.
    //!
    //! \returns A const reference to `Presentation<word_type>`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] std::vector<word_type>::const_iterator
    cbegin_long_rules() const noexcept {
      LIBSEMIGROUPS_ASSERT(_presentation.rules.cbegin() <= _longs_begin);
      LIBSEMIGROUPS_ASSERT(_longs_begin <= _presentation.rules.cend());
      return _longs_begin;
    }

    //! \brief Clear the set of long rules.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    Subclass& clear_long_rules() {
      return cbegin_long_rules(_presentation.rules.cend());
    }

    //! Returns the number of rules marked as long rules.
    [[nodiscard]] size_t number_of_long_rules() const noexcept {
      return std::distance(_longs_begin, _presentation.rules.cend()) / 2;
    }

    //! \brief Set the length of a long rule.
    //!
    //! Define the length of a "long" rule.
    //!
    //! This function modifies \ref presentation so that the rules whose length
    //! (sum of the lengths of both sides) is at least  \p val (if any) occur at
    //! the end of `presentation().rules` and so that `cbegin_long_rules` points
    //! at the such rule.
    //!
    //! The relative orders of the rules within \ref presentation
    //! may not be preserved.
    //!
    //! \param val the value of the long rule length.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    Subclass& long_rule_length(size_t val);

    //! \brief Get all active pruners of the search tree.
    //!
    //! Returns a const reference to the set of pruners.
    //!
    //! A pruner is any function that takes as input a word graph and returns a
    //! boolean. We require that if a pruner returns false for a word graph
    //! `wg`, then it returns false for all word graphs that are descended from
    //! `wg` in the Sims word graph search tree. The pruners are used to refine
    //! the congruence search tree during the execution of the Sims algorithm.
    //! As such, the congruences computed by this instance are only taken among
    //! those whose word graphs are accepted by all pruners returned by
    //! \ref pruners.
    //!
    //! \returns A const reference to
    //! `std::vector<std::function<bool(word_graph_type const&)>>`, the set of
    //! all pruners.
    //!
    //! \exceptions
    //! \noexcept
    auto const& pruners() const noexcept {
      return _pruners;
    }

    //! \brief Add a pruner to the search tree.
    //!
    //! \param func an rvalue reference to a pruner function.
    //!
    //! \returns A reference to \c this.
    //!
    //! \sa \ref pruners
    //!
    //! \warning When running the Sims low-index backtrack with multiple
    //! threads, each added pruner must be guaranteed thread safe. Failing to do
    //! so could cause bad things to happen.
    template <typename Func>
    Subclass& add_pruner(Func&& func) {
      _pruners.emplace_back(func);
      // TODO(2) could return an iterator to the inserted func in case we want
      // to remove it again
      return static_cast<Subclass&>(*this);
    }

    //! \brief Clear the set of pruners.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    Subclass& clear_pruners();

    //! \brief Get the set of pairs that must be included in every
    //! congruence.
    //!
    //! This function returns the set of included pairs.
    //!
    //! The congruences computed by a Sims1 or Sims2 instance always include the
    //! pairs returned by \ref included_pairs. In other words, the congruences
    //! computed by this instance are only taken among those that contains the
    //! pairs of elements of the underlying semigroup (defined by the
    //! presentation returned by \ref presentation and \ref long_rules)
    //! represented by the relations of the presentation returned by \ref
    //! included_pairs.
    //!
    //! \returns A const reference to `std::vector<word_type>`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] std::vector<word_type> const&
    included_pairs() const noexcept {
      return _include;
    }

    //! \brief Add an included pair via iterators.
    //!
    //! This function adds a pair that must be included in every congruence
    //! returned by a \ref Sims1 or \ref Sims2 instance.
    //!
    //! \cong_intf_params_contains
    //!
    //! \returns A reference to \c this.
    //!
    //! \cong_intf_warn_assume_letters_in_bounds
    template <typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    Subclass& add_included_pair_no_checks(Iterator1 first1,
                                          Iterator2 last1,
                                          Iterator3 first2,
                                          Iterator4 last2) {
      return include_exclude_no_checks(first1, last1, first2, last2, _include);
    }

    //! \brief Add an included pair via iterators (no checks).
    //!
    //! This function adds a pair that must be included in every congruence
    //! returned by a \ref Sims1 or \ref Sims2 instance.
    //!
    //! \cong_intf_params_contains
    //!
    //! \returns A reference to \c this.
    //!
    //! \cong_intf_throws_if_letters_out_of_bounds
    template <typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    Subclass& add_included_pair(Iterator1 first1,
                                Iterator2 last1,
                                Iterator3 first2,
                                Iterator4 last2) {
      return include_exclude(first1, last1, first2, last2, _include);
    }

    //! \brief Clear the set of included pairs.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    Subclass& clear_included_pairs() {
      _include.clear();
      return static_cast<Subclass&>(*this);
    }

    //! \brief Get the set of pairs that must be excluded from every
    //! congruence.
    //!
    //! This function returns the set of excluded pairs.
    //!
    //! The congruences computed by a Sims1 or Sims2 instance will never contain
    //! the pairs returned by \ref excluded_pairs. In other words, the
    //! congruences computed by this instance are only taken among those that do
    //! not contain any of the pairs of elements of the underlying semigroup
    //! (defined by the presentation returned by \ref presentation and \ref
    //! long_rules) represented by the relations of the presentation returned by
    //! \ref excluded_pairs.
    //!
    //! \returns A const reference to `std::vector<word_type>`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] std::vector<word_type> const&
    excluded_pairs() const noexcept {
      return _exclude;
    }

    //! \brief Add an excluded pair via iterators (no checks).
    //!
    //! This function adds a pair that must be excluded from every congruence
    //! returned by a \ref Sims1 or \ref Sims2 instance.
    //!
    //! \cong_intf_params_contains
    //!
    //! \returns A reference to \c this.
    //!
    //! \cong_intf_warn_assume_letters_in_bounds
    template <typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    Subclass& add_excluded_pair_no_checks(Iterator1 first1,
                                          Iterator2 last1,
                                          Iterator3 first2,
                                          Iterator4 last2) {
      add_exclude_pruner();
      return include_exclude_no_checks(first1, last1, first2, last2, _exclude);
    }

    //! \brief Add an excluded pair via iterators.
    //!
    //! This function adds a pair that must be excluded from every congruence
    //! returned by a \ref Sims1 or \ref Sims2 instance.
    //!
    //! \cong_intf_params_contains
    //!
    //! \returns A reference to \c this.
    //!
    //! \cong_intf_throws_if_letters_out_of_bounds
    template <typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    Subclass& add_excluded_pair(Iterator1 first1,
                                Iterator2 last1,
                                Iterator3 first2,
                                Iterator4 last2) {
      add_exclude_pruner();
      return include_exclude(first1, last1, first2, last2, _exclude);
    }

    //! \brief Clear the set of excluded words.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    Subclass& clear_excluded_pairs();

    //! \brief Get the current stats object.
    //!
    //! Returns a const reference to the current stats object.
    //!
    //! The value returned by this function is a `SimsStats` object which
    //! contains some statistics related to the current Sims1 or Sims2 instance
    //! and any part of the depth first search already conducted.
    //!
    //! \returns A const reference to `SimsStats`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] SimsStats& stats() const noexcept {
      return _stats;
    }

    //! \brief Get the idle thread restart attempt count.
    //!
    //! Returns the number of times an idle thread will attempt to restart
    //! before yielding during the execution of Sims1 or Sims2.
    //!
    //! \note This setting has no effect if \ref number_of_threads is 1.
    //!
    //! \returns A `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] size_t idle_thread_restarts() const noexcept {
      return _idle_thread_restarts;
    }

    //! \brief Set the idle thread restart attempt count.
    //!
    //! This function sets the idle thread restart attempt count used in Sims1
    //! or Sims2.
    //!
    //! The default value is `64`.
    //!
    //! \note This setting has no effect if \ref number_of_threads is 1.
    //!
    //! \param val the maximum number of times an idle thread will attempt to
    //! restart before yielding.
    //!
    //! \returns A reference to \c this.
    //!
    //! \throws LibsemigroupsException if the argument \p val is 0.
    // Number of times an idle thread will attempt to restart before yielding.
    Subclass& idle_thread_restarts(size_t val);

    //! \brief Throws if any letter in a range is out of bounds.
    //!
    //! This function throws a LibsemigroupsException if any value pointed at
    //! by an iterator in the range \p first to \p last is out of bounds (i.e.
    //! does not belong to the alphabet of the \ref presentation used to
    //! construct the \ref Sims1 or \ref Sims2 instance).
    //!
    //! \tparam Iterator1 the type of first argument \p first.
    //! \tparam Iterator2 the type of second argument \p last.
    //!
    //! \param first iterator pointing at the first letter of the word.
    //! \param last iterator pointing one beyond the last letter of the word.
    //!
    //! \throw LibsemigroupsException if any letter in the range from \p first
    //! to \p last is out of bounds.
    template <typename Iterator1, typename Iterator2>
    void throw_if_letter_out_of_bounds(Iterator1 first, Iterator2 last) const {
      presentation().validate_word(first, last);
    }

   protected:
    Subclass const& stats_copy_from(SimsStats const& stts) const {
      _stats = std::move(stts);
      return static_cast<Subclass const&>(*this);
    }

   private:
    template <typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    Subclass&
    include_exclude_no_checks(Iterator1               first1,
                              Iterator2               last1,
                              Iterator3               first2,
                              Iterator4               last2,
                              std::vector<word_type>& include_or_exclude) {
      include_or_exclude.emplace_back(first1, last1);
      include_or_exclude.emplace_back(first2, last2);
      return static_cast<Subclass&>(*this);
    }

    template <typename Iterator1,
              typename Iterator2,
              typename Iterator3,
              typename Iterator4>
    Subclass& include_exclude(Iterator1               first1,
                              Iterator2               last1,
                              Iterator3               first2,
                              Iterator4               last2,
                              std::vector<word_type>& include_or_exclude) {
      throw_if_letter_out_of_bounds(first1, last1);
      throw_if_letter_out_of_bounds(first2, last2);
      return static_cast<Subclass&>(*this).include_exclude_no_checks(
          first1, last1, first2, last2, include_or_exclude);
    }

    size_t add_exclude_pruner();
  };

  namespace sims {
    //! \brief Helper for adding an included pair of words.
    //!
    //! This function can be used to add an included pair to \p
    //! sims using the objects themselves rather than
    //! using iterators.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam Word the type of the second and third parameters.
    //!
    //! \param sims a \ref Sims1 or \ref Sims2 object.
    //! \param u the left hand side of the pair to add.
    //! \param v the right hand side of the pair to add.
    //!
    //! \return A reference to \p sims.
    //!
    //! \cong_intf_warn_assume_letters_in_bounds
    //!
    //! \sa \ref Sims1::included_pairs
    template <typename Subclass, typename Word>
    Subclass& add_included_pair_no_checks(SimsSettings<Subclass>& sims,
                                          Word const&             u,
                                          Word const&             v) {
      return sims.add_included_pair_no_checks(
          u.cbegin(), u.cend(), v.cbegin(), v.cend());
    }

    //! \brief Helper for adding a generating pair of words
    //! (std::initializer_list).
    //!
    //! See \ref
    //! add_included_pair_no_checks<Subclass, Word>(SimsSettings<Subclass>&,
    //! Word const&, Word const&) for details.
    template <typename Subclass, typename Int>
    Subclass& add_included_pair_no_checks(SimsSettings<Subclass>& sims,
                                          std::initializer_list<Int> const& u,
                                          std::initializer_list<Int> const& v) {
      return add_included_pair_no_checks<Subclass, std::vector<Int>>(
          sims, u, v);
    }

    //! \brief Helper for adding a generating pair of words
    //! (string literals).
    //!
    //! See \ref
    //! add_included_pair_no_checks<Subclass, Word>(SimsSettings<Subclass>&,
    //! Word const&, Word const&) for details.
    template <typename Subclass>
    Subclass& add_included_pair_no_checks(SimsSettings<Subclass>& sims,
                                          char const*             u,
                                          char const*             v) {
      return sims.add_included_pair_no_checks(
          u, u + std::strlen(u), v, v + std::strlen(v));
    }

    // This version of the function catches the cases when u & v are not of the
    // same type but both convertible to string_view
    //! \brief Helper for adding a generating pair of words
    //! (std::string_view).
    //!
    //! See \ref
    //! add_included_pair_no_checks<Subclass, Word>(SimsSettings<Subclass>&,
    //! Word const&, Word const&) for details.
    template <typename Subclass>
    Subclass& add_included_pair_no_checks(SimsSettings<Subclass>& sims,
                                          std::string_view        u,
                                          std::string_view        v) {
      return sims.add_included_pair_no_checks(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    //! \brief Helper for adding an included pair of words.
    //!
    //! This function can be used to add an included pair to \p
    //! sims using the objects themselves rather than
    //! using iterators.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam Word the type of the second and third parameters.
    //!
    //! \param sims a \ref Sims1 or \ref Sims2 object.
    //! \param u the left hand side of the pair to add.
    //! \param v the right hand side of the pair to add.
    //!
    //! \return A reference to \p sims.
    //!
    //! \cong_intf_throws_if_letters_out_of_bounds
    //!
    //! \sa \ref Sims1::included_pairs
    template <typename Subclass, typename Word>
    Subclass& add_included_pair(SimsSettings<Subclass>& sims,
                                Word const&             u,
                                Word const&             v) {
      return sims.add_included_pair(u.cbegin(), u.cend(), v.cbegin(), v.cend());
    }

    //! \brief Helper for adding a included pair of words
    //! (std::initializer_list).
    //!
    //! See \ref
    //! add_included_pair<Subclass, Word>(SimsSettings<Subclass>&, Word const&,
    //! Word const&) for details.
    template <typename Subclass, typename Int>
    Subclass& add_included_pair(SimsSettings<Subclass>&           sims,
                                std::initializer_list<Int> const& u,
                                std::initializer_list<Int> const& v) {
      return add_included_pair<Subclass, std::vector<Int>>(sims, u, v);
    }

    //! \brief Helper for adding a included pair of words
    //! (string literals).
    //!
    //! See \ref
    //! add_included_pair<Subclass, Word>(SimsSettings<Subclass>&, Word const&,
    //! Word const&) for details.
    template <typename Subclass>
    Subclass& add_included_pair(SimsSettings<Subclass>& sims,
                                char const*             u,
                                char const*             v) {
      return sims.add_included_pair(
          u, u + std::strlen(u), v, v + std::strlen(v));
    }

    // This version of the function catches the cases when u & v are not of the
    // same type but both convertible to string_view
    //! \brief Helper for adding a included pair of words
    //! (std::string_view).
    //!
    //! See \ref
    //! add_included_pair<Subclass, Word>(SimsSettings<Subclass>&, Word const&,
    //! Word const&) for details.
    template <typename Subclass>
    Subclass& add_included_pair(SimsSettings<Subclass>& sims,
                                std::string_view        u,
                                std::string_view        v) {
      return sims.add_included_pair(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    //! \brief Helper for adding an excluded pair of words.
    //!
    //! This function can be used to add an excluded pair to \p
    //! sims using the objects themselves rather than
    //! using iterators.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam Word the type of the second and third parameters.
    //!
    //! \param sims a \ref Sims1 or \ref Sims2 object.
    //! \param u the left hand side of the pair to add.
    //! \param v the right hand side of the pair to add.
    //!
    //! \return A reference to \p sims.
    //!
    //! \cong_intf_warn_assume_letters_in_bounds
    //!
    //! \sa \ref Sims1::excluded_pairs
    template <typename Subclass, typename Word>
    Subclass& add_excluded_pair_no_checks(SimsSettings<Subclass>& sims,
                                          Word const&             u,
                                          Word const&             v) {
      return sims.add_excluded_pair_no_checks(
          u.cbegin(), u.cend(), v.cbegin(), v.cend());
    }

    //! \brief Helper for adding a excluded pair of words
    //! (std::initializer_list).
    //!
    //! See \ref
    //! add_excluded_pair_no_checks<Subclass, Word>(SimsSettings<Subclass>&,
    //! Word const&, Word const&) for details.
    template <typename Subclass, typename Int>
    Subclass& add_excluded_pair_no_checks(SimsSettings<Subclass>& sims,
                                          std::initializer_list<Int> const& u,
                                          std::initializer_list<Int> const& v) {
      return add_excluded_pair_no_checks<Subclass, std::vector<Int>>(
          sims, u, v);
    }

    //! \brief Helper for adding a excluded pair of words
    //! (string literals).
    //!
    //! See \ref
    //! add_excluded_pair_no_checks<Subclass, Word>(SimsSettings<Subclass>&,
    //! Word const&, Word const&) for details.
    template <typename Subclass>
    Subclass& add_excluded_pair_no_checks(SimsSettings<Subclass>& sims,
                                          char const*             u,
                                          char const*             v) {
      return sims.add_excluded_pair_no_checks(
          u, u + std::strlen(u), v, v + std::strlen(v));
    }

    // This version of the function catches the cases when u & v are not of the
    // same type but both convertible to string_view
    //! \brief Helper for adding a excluded pair of words
    //! (std::string_view).
    //!
    //! See \ref
    //! add_excluded_pair_no_checks<Subclass, Word>(SimsSettings<Subclass>&,
    //! Word const&, Word const&) for details.
    template <typename Subclass>
    Subclass& add_excluded_pair_no_checks(SimsSettings<Subclass>& sims,
                                          std::string_view        u,
                                          std::string_view        v) {
      return sims.add_excluded_pair_no_checks(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

    //! \brief Helper for adding an excluded pair of words.
    //!
    //! This function can be used to add an excluded pair to \p
    //! sims using the objects themselves rather than
    //! using iterators.
    //!
    //! \tparam Subclass the type of the first parameter.
    //! \tparam Word the type of the second and third parameters.
    //!
    //! \param sims a \ref Sims1 or \ref Sims2 object.
    //! \param u the left hand side of the pair to add.
    //! \param v the right hand side of the pair to add.
    //!
    //! \return A reference to \p sims.
    //!
    //! \cong_intf_throws_if_letters_out_of_bounds
    //!
    //! \sa \ref Sims1::excluded_pairs
    template <typename Subclass, typename Word>
    Subclass& add_excluded_pair(SimsSettings<Subclass>& sims,
                                Word const&             u,
                                Word const&             v) {
      return sims.add_excluded_pair(u.cbegin(), u.cend(), v.cbegin(), v.cend());
    }

    //! \brief Helper for adding a excluded pair of words
    //! (std::initializer_list).
    //!
    //! See \ref
    //! add_excluded_pair<Subclass, Word>(SimsSettings<Subclass>&, Word const&,
    //! Word const&) for details.
    template <typename Subclass, typename Int>
    Subclass& add_excluded_pair(SimsSettings<Subclass>&           sims,
                                std::initializer_list<Int> const& u,
                                std::initializer_list<Int> const& v) {
      return add_excluded_pair<Subclass, std::vector<Int>>(sims, u, v);
    }

    //! \brief Helper for adding a excluded pair of words
    //! (string literals).
    //!
    //! See \ref
    //! add_excluded_pair<Subclass, Word>(SimsSettings<Subclass>&, Word const&,
    //! Word const&) for details.
    template <typename Subclass>
    Subclass& add_excluded_pair(SimsSettings<Subclass>& sims,
                                char const*             u,
                                char const*             v) {
      return sims.add_excluded_pair(
          u, u + std::strlen(u), v, v + std::strlen(v));
    }

    // This version of the function catches the cases when u & v are not of the
    // same type but both convertible to string_view
    //! \brief Helper for adding a excluded pair of words
    //! (std::string_view).
    //!
    //! See \ref
    //! add_excluded_pair<Subclass, Word>(SimsSettings<Subclass>&, Word const&,
    //! Word const&) for details.
    template <typename Subclass>
    Subclass& add_excluded_pair(SimsSettings<Subclass>& sims,
                                std::string_view        u,
                                std::string_view        v) {
      return sims.add_excluded_pair(
          std::begin(u), std::end(u), std::begin(v), std::end(v));
    }

  }  // namespace sims

  ////////////////////////////////////////////////////////////////////////
  // Sims1 - Sims2 - forward decl
  ////////////////////////////////////////////////////////////////////////

  class Sims1;
  class Sims2;

  namespace detail {
    template <typename Sims1or2>
    class SimsBase : public SimsSettings<Sims1or2>, public Reporter {
      static_assert(std::is_same_v<Sims1or2, Sims1>
                    || std::is_same_v<Sims1or2, Sims2>);

     public:
      //! The type of the associated WordGraph objects.
      // We use WordGraph, even though the iterators produced by this class
      // hold FelschGraph's, none of the features of FelschGraph are useful
      // for the output, only for the implementation
      using word_graph_type = typename SimsSettings<Sims1or2>::word_graph_type;

      //! Type for the nodes in the associated WordGraph objects.
      using node_type = typename word_graph_type::node_type;

      //! Type for the labels in the associated WordGraph objects.
      using label_type = typename word_graph_type::label_type;

      //! The size_type of the associated WordGraph objects.
      using size_type = typename word_graph_type::size_type;

      //! Type for letters in the underlying presentation.
      using letter_type = typename word_type::value_type;

      ////////////////////////////////////////////////////////////////////////
      // SimsBase nested classes - protected
      ////////////////////////////////////////////////////////////////////////

      // This class collects some common aspects of the iterator and
      // thread_iterator nested classes. The mutex does nothing for <iterator>
      // and is an actual std::mutex for <thread_iterator>. Also subclassed by
      // Sims2::iterator_base.
      class IteratorBase {
       public:
        using const_reference = word_graph_type const&;
        using const_pointer   = word_graph_type const*;
        using sims_type       = Sims1or2;

       private:
        size_type _max_num_classes;
        size_type _min_target_node;

       protected:
        using PendingDef = typename Sims1or2::PendingDef;
        using Definition = std::pair<node_type, label_type>;

        // TODO(1) (Sims1) ensure that _felsch_graph's settings are
        // properly initialised
        FelschGraph<word_type, node_type, std::vector<Definition>>
            _felsch_graph;

        // This mutex does nothing for iterator, only does something for
        // thread_iterator
        std::mutex              _mtx;
        std::vector<PendingDef> _pending;
        Sims1or2 const*         _sims1or2;

        // Push initial PendingDef's into _pending, see tpp file for
        // explanation.
        void init(size_type n);

        // We could use the copy constructor, but there's no point in copying
        // anything except the FelschGraph and so we only copy that.
        void partial_copy_for_steal_from(IteratorBase const& that) {
          _felsch_graph = that._felsch_graph;
        }

        // Try to pop from _pending into the argument (reference), returns
        // true if successful and false if not.
        [[nodiscard]] bool try_pop(PendingDef& pd);

        // Try to make the definition represented by PendingDef, returns false
        // if it wasn't possible, and true if it was.
        bool try_define(PendingDef const& current);

        // Install any new pending definitions arising from the definition of
        // "current", this should only be called after "try_define(current)",
        // and is in a separate function so that we can define another version
        // of "try_define" in "Sims2".
        bool install_descendents(PendingDef const& current);

        IteratorBase(Sims1or2 const* s, size_type n);

       public:
        // None of the constructors are noexcept because the corresponding
        // constructors for Presentation aren't currently
        IteratorBase();
        IteratorBase(IteratorBase const& that);
        IteratorBase(IteratorBase&& that);
        IteratorBase& operator=(IteratorBase const& that);
        IteratorBase& operator=(IteratorBase&& that);
        ~IteratorBase();

        [[nodiscard]] bool operator==(IteratorBase const& that) const noexcept {
          return _felsch_graph == that._felsch_graph;
        }

        [[nodiscard]] bool operator!=(IteratorBase const& that) const noexcept {
          return !(operator==(that));
        }

        [[nodiscard]] const_reference operator*() const noexcept {
          return _felsch_graph;
        }

        [[nodiscard]] const_pointer operator->() const noexcept {
          return &_felsch_graph;
        }

        void swap(IteratorBase& that) noexcept;

        SimsStats& stats() noexcept {
          return _sims1or2->stats();
        }

        size_type maximum_number_of_classes() const noexcept {
          return _max_num_classes;
        }

        Sims1or2 const& sims() const noexcept {
          return *_sims1or2;
        }
      };  // class IteratorBase

     public:
      ////////////////////////////////////////////////////////////////////////
      // SimsBase nested classes - public
      ////////////////////////////////////////////////////////////////////////

      //! The return type of \ref cbegin and \ref cend.
      //!
      //! This is a forward iterator values of this type are expensive to copy
      //! due to their internal state, and prefix increment should be
      //! preferred to postfix.
      class iterator : public Sims1or2::iterator_base {
        using iterator_base = typename Sims1or2::iterator_base;
        using PendingDef    = typename Sims1or2::PendingDef;

        using iterator_base::init;
        using iterator_base::install_descendents;
        using iterator_base::try_define;
        using iterator_base::try_pop;

       public:
        using const_pointer   = typename iterator_base::const_pointer;
        using const_reference = typename iterator_base::const_reference;
        using size_type = typename std::vector<word_graph_type>::size_type;
        using difference_type =
            typename std::vector<word_graph_type>::difference_type;
        using pointer    = typename std::vector<word_graph_type>::pointer;
        using reference  = typename std::vector<word_graph_type>::reference;
        using value_type = word_graph_type;
        using iterator_category = std::forward_iterator_tag;

        using sims_type = typename iterator_base::sims_type;

        using iterator_base::iterator_base;

       private:
        // Only want SimsBase to be able to use this constructor.
        iterator(Sims1or2 const* s, size_type n);

        // So that we can use the constructor above.
        friend iterator SimsBase::cbegin(SimsBase::size_type) const;
        friend iterator SimsBase::cend(SimsBase::size_type) const;

       public:
        ~iterator();

        iterator() : iterator_base() {}
        iterator(iterator const& that) : iterator_base(that) {}
        iterator(iterator&& that) : iterator_base(std::move(that)) {}
        iterator& operator=(iterator const& that) {
          iterator_base::operator=(that);
          return *this;
        }
        iterator& operator=(iterator&& that) {
          iterator_base::operator=(std::move(that));
          return *this;
        }

        //! No doc
        // prefix
        iterator const& operator++();

        //! No doc
        // postfix
        iterator operator++(int) {
          return detail::default_postfix_increment(*this);
        }

        using iterator_base::swap;
      };  // class iterator

     private:
      class thread_runner;
      class thread_iterator;

      void report_at_start(size_t num_classes) const;
      void report_progress_from_thread() const;
      void report_final() const;

      void throw_if_not_ready(size_type n) const;

     public:
      SimsBase();
      SimsBase(SimsBase const& other)      = default;
      SimsBase(SimsBase&&)                 = default;
      SimsBase& operator=(SimsBase const&) = default;
      SimsBase& operator=(SimsBase&&)      = default;
      ~SimsBase();

      Sims1or2& init();

      // Required because we are inheriting from Reporter and SimsSettings
      using SimsSettings<Sims1or2>::presentation;
      using SimsSettings<Sims1or2>::number_of_threads;
      using SimsSettings<Sims1or2>::stats;
      using SimsSettings<Sims1or2>::included_pairs;
      using SimsSettings<Sims1or2>::excluded_pairs;
      using SimsSettings<Sims1or2>::cbegin_long_rules;

      [[nodiscard]] iterator cbegin(size_type n) const {
        throw_if_not_ready(n);
        return iterator(static_cast<Sims1or2 const*>(this), n);
      }

      [[nodiscard]] iterator cend(size_type n) const {
        throw_if_not_ready(n);
        return iterator(static_cast<Sims1or2 const*>(this), 0);
      }

      // Apply the function pred to every congruence with at
      // most n classes
      void for_each(size_type                                   n,
                    std::function<void(word_graph_type const&)> pred) const;

      word_graph_type
      find_if(size_type                                   n,
              std::function<bool(word_graph_type const&)> pred) const;

      uint64_t number_of_congruences(size_type n) const;
    };  // SimsBase
  }     // namespace detail

  namespace sims {
    class const_cgp_iterator;
    class const_rcgp_iterator;
  }  // namespace sims

  //! \ingroup congruences_group
  //!
  //! \brief For computing finite index right congruences of a finitely
  //! presented semigroup or monoid.
  //!
  //! Defined in `sims.hpp`.
  //!
  //! On this page we describe the functionality relating to the small index
  //! congruence algorithm for 1-sided congruences. The algorithm implemented by
  //! this class is essentially the low index subgroup algorithm for
  //! finitely presented groups described in Section 5.6 of [Computation with
  //! Finitely Presented Groups](https://doi.org/10.1017/CBO9780511574702) by C.
  //! Sims. The low index subgroups algorithm was adapted for semigroups and
  //! monoids by R. Cirpons, J. D. Mitchell, and M. Tsalakou; see
  //! \cite Anagnostopoulou-Merkouri2023aa
  //!
  //! The purpose of this class is to provide the functions \ref cbegin, \ref
  //! cend, \ref for_each, and \ref find_if which permit iterating through the
  //! one-sided congruences of a semigroup or monoid defined by a presentation
  //! containing (a possibly empty) set of pairs and with at most a given number
  //! of classes. An iterator returned by \ref cbegin points at an WordGraph
  //! instance containing the action of the semigroup or monoid on the classes
  //! of a congruence.
  //!
  //! \sa Sims2 for equivalent functionality for 2-sided congruences.
  //! \sa SimsSettings for the various things that can be set in a Sims1 object.
  class Sims1 : public detail::SimsBase<Sims1> {
    // Aliases
    using SimsBase      = detail::SimsBase<Sims1>;
    using iterator_base = IteratorBase;

    // Friends
    // so that SimsBase can access iterator_base, PendingDef, etc
    friend SimsBase;
    // so that Sims2 can access PendingDef
    friend Sims2;

    // Forward decl
    struct PendingDef;

   public:
    // We use WordGraph, even though the iterators produced by this class
    // hold FelschGraph's, none of the features of FelschGraph are useful
    // for the output, only for the implementation
    //! The type of the associated WordGraph objects.
    using word_graph_type = SimsBase::word_graph_type;

    //! Type of the nodes in the associated WordGraph objects.
    using node_type = word_graph_type::node_type;

    //! Type of the edge labels in the associated WordGraph objects.
    using label_type = word_graph_type::label_type;

    //! Type for letters in the underlying presentation.
    using letter_type = word_type::value_type;

    //! The \ref WordGraph::size_type of the associated WordGraph objects.
    using size_type = word_graph_type::size_type;

   private:
    using Definition = std::pair<node_type, label_type>;

    friend class sims::const_rcgp_iterator;
    friend class sims::const_cgp_iterator;

    using felsch_graph_type
        = detail::FelschGraph<word_type, node_type, std::vector<Definition>>;

   public:
    //! Default constructor
    Sims1() = default;

    using SimsBase::init;

    //! \brief Construct from a presentation.
    //!
    //! Constructs an instance from a presentation of any kind.
    //!
    //! The rules of the presentation \p p are used at every node in the depth
    //! first search conducted by an object of this type.
    //!
    //! If the template parameter \p Word is not \ref word_type, then
    //! the parameter \p p is first converted to a value of type
    //! `Presentation<word_type>` and it is this converted value that is used.
    //!
    //! \throws LibsemigroupsException if `to_presentation<word_type>(p)` throws
    //! \throws LibsemigroupsException if `p` is not valid
    //! \throws LibsemigroupsException if `p` has 0-generators and 0-relations.
    //!
    //! \tparam Word the type of the words in the presentation \p p .
    //!
    //! \param p the presentation.
    //!
    //! \sa presentation
    //! \sa init
    template <typename Word>
    explicit Sims1(Presentation<Word> const& p) : Sims1() {
      presentation(p);
    }

    //! \copydoc Sims1::Sims1(Presentation<Word> const&)
    template <typename Word>
    explicit Sims1(Presentation<Word> const&& p) : Sims1() {
      presentation(std::move(p));
    }

    //! Default copy constructor.
    Sims1(Sims1 const& other) = default;

    //! Default move constructor.
    Sims1(Sims1&&) = default;

    //! Default copy assignment operator.
    Sims1& operator=(Sims1 const&) = default;

    //! Default move assignment operator.
    Sims1& operator=(Sims1&&) = default;

    //! No doc
    ~Sims1() = default;

    //! \brief Reinitialize an existing Sims1 object.
    //!
    //! This function puts an object back into the same state as if it had
    //! been newly constructed from the presentation \p p.
    //!
    //! \tparam Word the type of the words in the presentation \p p .
    //!
    //! \param p the presentation.
    //!
    //! \returns A reference to \c *this.
    //!
    //! \throws LibsemigroupsException if `to_presentation<word_type>(p)` throws
    //! \throws LibsemigroupsException if `p` is not valid
    //! \throws LibsemigroupsException if `p` has 0-generators and 0-relations.
    //!
    //! \warning This function has no exception guarantee, the object will be
    //! in the same state as if it was default constructed if an exception is
    //! thrown.
    //!
    //! \sa presentation(Presentation<Word> const&)
    template <typename Word>
    Sims1& init(Presentation<Word> const& p) {
      init();
      presentation(p);
      return *this;
    }

    //! \copydoc Sims1::init(Presentation<Word> const&)
    template <typename Word>
    Sims1& init(Presentation<Word> const&& p) {
      init();
      presentation(std::move(p));
      return *this;
    }

#ifdef PARSED_BY_DOXYGEN
    //! \brief Reinitialize an existing Sims1 object.
    //!
    //! This function puts a Sims1 object back into the same state as if
    //! it had been newly default constructed.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    Sims1& init();

    //! \brief Returns the number of one-sided congruences with up to a given
    //! number of classes.
    //!
    //! This function is similar to `std::distance(begin(n), end(n))` and
    //! exists to:
    //! * provide some feedback on the progress of the computation if it runs
    //!   for more than 1 second.
    //! * allow for the computation of `std::distance(begin(n), end(n))` to be
    //!   performed using \ref number_of_threads in parallel.
    //!
    //! \param n the maximum number of congruence classes.
    //!
    //! \returns A value of type \c uint64_t.
    //!
    //! \throws LibsemigroupsException if \p n is \c 0.
    //! \throws LibsemigroupsException if `presentation()` has 0-generators
    //! and 0-relations (i.e. it has not been initialised).
    uint64_t number_of_congruences(size_t n);

    //! \brief Apply a unary predicate to every one-sided congruence with at
    //! most a given number of classes.
    //!
    //! Apply the function \p pred to every one-sided congruence with at most
    //! \p n classes
    //!
    //!  This function is similar to `std::for_each(begin(n), end(n), pred)` and
    //!  exists to:
    //!  * provide some feedback on the progress of the computation if it runs
    //!  for more than 1 second.
    //!  * allow for the computation of `std::for_each(begin(n), end(n), pred)`
    //!  to be performed using \ref number_of_threads in parallel.
    //!
    //!  \param n the maximum number of congruence classes.
    //!  \param pred the predicate applied to every congruence found.
    //!
    //!  \throws LibsemigroupsException if \p n is \c 0.
    //!  \throws LibsemigroupsException if `presentation()` has 0-generators and
    //!  0-relations (i.e. it has not been initialised).
    //!
    //!  \sa cbegin
    void for_each(size_type                                   n,
                  std::function<void(word_graph_type const&)> pred) const;

    //! \brief Apply a unary predicate to one-sided congruences with at  most a
    //! given number of classes, until it returns \c true.
    //!
    //! Apply the function \p pred to every one-sided congruence with at most \p
    //! n classes, until it returns \c true.
    //!
    //! This function is similar to `std::find_if(begin(n),
    //! end(n), pred)` and exists to:
    //! * provide some feedback on the progress of the computation if it runs
    //!   for more than 1 second.
    //! * allow for the computation of `std::find_if(begin(n), end(n), pred)`
    //!   to be performed using \ref number_of_threads in parallel.
    //!
    //! \param n the maximum number of congruence classes.
    //! \param pred the predicate applied to every congruence found.
    //!
    //! \returns The first WordGraph for which \p pred returns \c true.
    //!
    //! \throws LibsemigroupsException if \p n is \c 0.
    //! \throws LibsemigroupsException if `presentation()` has 0-generators and
    //! 0-relations (i.e. it has not been initialised).
    //!
    //! \sa cbegin
    word_graph_type
    find_if(size_type                                   n,
            std::function<bool(word_graph_type const&)> pred) const;

    //! \brief Returns a forward iterator pointing at the first congruence.
    //!
    //! Returns a forward iterator pointing to the WordGraph representing the
    //! first congruence described by an object of this type with at most \p n
    //! classes.
    //!
    //! If incremented, the iterator will point to the next such congruence.
    //! The order in which the congruences are returned is implementation
    //! specific. Iterators of the type returned by this function are equal
    //! whenever they point to equal objects. The iterator is exhausted if
    //! and only if it points to an WordGraph with zero nodes.
    //!
    //! The meaning of the WordGraph pointed at by the returned iterator depends
    //! on whether the input is a monoid presentation (i.e.
    //! Presentation::contains_empty_word() returns \c true) or a semigroup
    //! presentation. If the input is a monoid presentation for a monoid
    //! \f$M\f$, then the WordGraph pointed to by an iterator of this type has
    //! precisely \p n nodes, and the right action of \f$M\f$ on the nodes of
    //! the word graph is isomorphic to the action of \f$M\f$ on the classes of
    //! a right congruence.
    //!
    //! If the input is a semigroup presentation for a semigroup \f$S\f$,
    //! then the WordGraph has \p n + 1 nodes, and the right action of
    //! \f$S\f$ on the nodes \f$\{1, \ldots, n\}\f$ of the WordGraph is
    //! isomorphic to the action of \f$S\f$ on the classes of a right
    //! congruence. It'd probably be better in this case if node \f$0\f$ was
    //! not included in the output WordGraph, but it is required in the
    //! implementation of the low-index congruence algorithm, and to avoid
    //! unnecessary copies, we've left it in for the time being. \param n
    //! the maximum number of classes in a congruence.
    //!
    //! \returns An iterator \c it of type \c iterator pointing to an WordGraph
    //! with at most \p n or \p n + 1 nodes.
    //!
    //! \throws LibsemigroupsException if \p n is \c 0.
    //! \throws LibsemigroupsException if `presentation()`
    //! has 0-generators and 0-relations (i.e. it has not
    //! been initialised).
    //!
    //! \warning
    //! Copying iterators of this type is expensive.  As a consequence,
    //! prefix incrementing \c ++it the returned  iterator \c it
    //! significantly cheaper than postfix incrementing \c it++.
    //!
    //! \sa cend
    // TODO(2) (Sims1) it'd be good to remove node 0 to avoid confusion. This
    // seems complicated however, and so isn't done at present.
    [[nodiscard]] iterator cbegin(size_type n) const;

    //! \brief Returns a forward iterator pointing one beyond the last
    //! congruence.
    //!
    //! Returns a forward iterator pointing to the empty WordGraph. If
    //! incremented, the returned iterator remains valid and continues to
    //! point at the empty WordGraph.
    //!
    //! \param n the maximum number of classes in a
    //! congruence.
    //!
    //! \returns An iterator \c it of type \c iterator pointing to an WordGraph
    //! with at most \p 0 nodes.
    //!
    //! \throws LibsemigroupsException if \p n is \c 0.
    //! \throws LibsemigroupsException if `presentation()` has 0-generators
    //! and 0-relations (i.e. it has not been initialised).
    //!
    //! \warning
    //! Copying iterators of this type is expensive.  As a consequence,
    //! prefix incrementing \c ++it the returned  iterator \c it
    //! significantly cheaper than postfix incrementing \c it++.
    //!
    //! \sa cbegin
    [[nodiscard]] iterator cend(size_type n) const;
#endif
  };

  //! \ingroup congruences_group
  //!
  //! \brief For computing finite index two-sided congruences of a finitely
  //! presented semigroup or monoid.
  //!
  //! Defined in `sims.hpp`.
  //!
  //! On this page we describe the functionality relating to the small index
  //! congruence algorithm for 2-sided congruences. The algorithm implemented by
  //! this class is described in \cite Anagnostopoulou-Merkouri2023aa.
  //!
  //! The purpose of this class is to provide the functions \ref cbegin, \ref
  //! cend, \ref for_each, and \ref find_if which permit iterating through the
  //! two-sided congruences of a semigroup or monoid defined by a presentation
  //! containing, or not containing, (possibly empty) sets of pairs and with at
  //! most a given number of classes. An iterator returned by \ref cbegin points
  //! at an WordGraph instance containing the action of the semigroup or monoid
  //! on the classes of a congruence.
  //!
  //! \sa Sims1 for equivalent functionality for 1-sided congruences.
  //! \sa SimsSettings for the various things that can be set in a Sims2 object.
  class Sims2 : public detail::SimsBase<Sims2> {
    using SimsBase = detail::SimsBase<Sims2>;
    // so that SimsBase can access iterator_base, PendingDef, etc
    friend SimsBase;

   public:
    //! The type of the associated WordGraph objects.
    using word_graph_type = SimsBase::word_graph_type;

    //! Type of the nodes in the associated WordGraph objects.
    using node_type = word_graph_type::node_type;

    //! Type of the edge labels in the associated WordGraph objects.
    using label_type = word_graph_type::label_type;

    //! Type for letters in the underlying presentation.
    using letter_type = word_type::value_type;

    //! The \ref WordGraph::size_type of the associated WordGraph objects.
    using size_type = word_graph_type::size_type;

    using SimsBase::init;

    //! Default constructor.
    Sims2() = default;
    //! Default copy constructor.
    Sims2(Sims2 const& other) = default;
    //! Default move constructor.
    Sims2(Sims2&&) = default;
    //! Default copy assignment operator.
    Sims2& operator=(Sims2 const&) = default;
    //! Default move assignment operator.
    Sims2& operator=(Sims2&&) = default;

    ~Sims2() = default;

    //! \copydoc Sims1::Sims1(Presentation<Word> const&)
    template <typename Word>
    explicit Sims2(Presentation<Word> const& p) : Sims2() {
      presentation(p);
    }

    //! \copydoc Sims1::Sims1(Presentation<Word> const&)
    template <typename Word>
    explicit Sims2(Presentation<Word> const&& p) : Sims2() {
      presentation(std::move(p));
    }

    //! \brief Reinitialize an existing Sims2 object.
    //!
    //! This function puts an object back into the same state as if it had
    //! been newly constructed from the presentation \p p.
    //!
    //! \tparam Word the type of the words in the presentation \p p .
    //!
    //! \param p the presentation.
    //!
    //! \returns A reference to \c *this.
    //!
    //! \throws LibsemigroupsException if `to_presentation<word_type>(p)` throws
    //! \throws LibsemigroupsException if `p` is not valid
    //! \throws LibsemigroupsException if `p` has 0-generators and 0-relations.
    //!
    //! \warning This function has no exception guarantee, the object will be
    //! in the same state as if it was default constructed if an exception is
    //! thrown.
    //!
    //! \sa presentation(Presentation<Word> const&)
    template <typename Word>
    Sims2& init(Presentation<Word> const& p) {
      init();
      presentation(p);
      return *this;
    }

    //! \copydoc Sims2::init(Presentation<Word> const&)
    template <typename Word>
    Sims2& init(Presentation<Word> const&& p) {
      init();
      presentation(std::move(p));
      return *this;
    }

#ifdef PARSED_BY_DOXYGEN
    //! \brief Reinitialize an existing Sims2 object.
    //!
    //! This function puts a Sims2 object back into the same state as if
    //! it had been newly default constructed.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    Sims2& init();

    //! \copydoc Sims1::number_of_congruences
    uint64_t number_of_congruences(size_t n);

    //! \copydoc Sims1::for_each
    void for_each(size_type                                   n,
                  std::function<void(word_graph_type const&)> pred) const;

    //! \copydoc Sims1::find_if
    word_graph_type
    find_if(size_type                                   n,
            std::function<bool(word_graph_type const&)> pred) const;

    //! \copydoc Sims1::cbegin
    [[nodiscard]] iterator cbegin(size_type n) const;

    //! \copydoc Sims1::cend
    [[nodiscard]] iterator cend(size_type n) const;
#endif

   private:
    struct PendingDef;

    // This class collects some common aspects of the iterator and
    // thread_iterator nested classes. The mutex does nothing for <iterator>
    // and is an actual std::mutex for <thread_iterator>.
    class iterator_base : public SimsBase::IteratorBase {
      class RuleContainer;

     public:
      using const_reference = SimsBase::IteratorBase::const_reference;
      using const_pointer   = SimsBase::IteratorBase::const_pointer;

     private:
      std::unique_ptr<RuleContainer> _2_sided_include;
      std::vector<word_type>         _2_sided_words;

     protected:
      // We could use the copy constructor, but there's no point in copying
      // anything except the FelschGraph and so we only copy that.
      void partial_copy_for_steal_from(iterator_base const& that);

      // Try to make the definition represented by PendingDef, returns false if
      // it wasn't possible, and true if it was.
      //! No doc
      [[nodiscard]] bool try_define(PendingDef const&);

      iterator_base(Sims2 const* s, size_type n);

     public:
      iterator_base();
      iterator_base(iterator_base const& that);
      iterator_base(iterator_base&& that);
      iterator_base& operator=(iterator_base const& that);
      iterator_base& operator=(iterator_base&& that);
      ~iterator_base();

      //! No doc
      void swap(iterator_base& that) noexcept;

      using SimsBase::IteratorBase::stats;
    };  // class iterator_base
  };    // Sims2

  //! \ingroup congruences_group
  //!
  //! \brief For computing small degree transformation representations of a
  //! finite semigroup or monoid.
  //!
  //! Defined in `sims.hpp`.
  //!
  //! This class is a helper for Sims1 calling the word_graph member
  //! function attempts to find a right congruence, represented as an
  //! WordGraph, of the semigroup or monoid defined by the presentation
  //! consisting of its \ref presentation and \ref long_rules with the
  //! following properties:
  //! * the transformation semigroup defined by the WordGraph has size
  //!   \ref target_size;
  //! * the number of nodes in the WordGraph is at least \ref min_nodes
  //!   and at most \ref max_nodes.
  //!
  //! If no such WordGraph can be found, then an empty WordGraph is
  //! returned (with `0` nodes and `0` edges).
  class RepOrc : public SimsSettings<RepOrc> {
   private:
    size_t _min;
    size_t _max;
    size_t _size;

   public:
    //! Default constructor.
    RepOrc() : _min(), _max(), _size() {
      init();
    }

    //! \brief Reinitialize an existing RepOrc object.
    //!
    //! This function puts a RepOrc object back into the same state as if
    //! it had been newly default constructed.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    RepOrc& init();

    //! \brief Construct from Sims1, Sims2 or MinimalRepOrc.
    //!
    //! This constructor creates a new RepOrc instance with
    //! the same SimsSettings as \p s but that is
    //! otherwise uninitialised.
    //!
    //! \tparam S the type of the argument \p s (which is
    //! derived from `SimsSettings<S>`).
    //!
    //! \param s the Sims1, Sims2 or MinimalRepOrc whose settings
    //! should be used.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    template <typename OtherSubclass>
    explicit RepOrc(SimsSettings<OtherSubclass> const& s) : RepOrc() {
      SimsSettings<RepOrc>::init(s);
    }

    //! \brief Initialize an existing RepOrc object from Sims1, Sims2 or
    //! MinimalRepOrc.
    //!
    //! This function reinitializes a RepOrc instance with
    //! the same SimsSettings as \p s .
    //!
    //! \tparam S the type of the argument \p s (which is
    //! derived from `SimsSettings<S>`).
    //!
    //! \param s the Sims1, Sims2 or MinimalRepOrc whose settings
    //! should be used.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    template <typename OtherSubclass>
    RepOrc& init(SimsSettings<OtherSubclass> const& s) {
      SimsSettings<RepOrc>::init(s);
      return *this;
    }

    //! \brief Set the minimum number of nodes.
    //!
    //! This function sets the minimal number of nodes in
    //! the WordGraph that we are seeking.
    //!
    //! \param val the minimum number of nodes
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \noexcept
    RepOrc& min_nodes(size_t val) noexcept {
      _min = val;
      return *this;
    }

    //! \brief Get the current minimum number of nodes.
    //!
    //! This function returns the current value for the minimum number of nodes
    //! in the WordGraph that we are seeking.
    //!
    //! \returns A value of type `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] size_t min_nodes() const noexcept {
      return _min;
    }

    //! \brief Set the maximum number of nodes.
    //!
    //! This function sets the maximum number of nodes in the WordGraph that we
    //! are seeking.
    //!
    //! \param val the maximum number of nodes
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \noexcept
    RepOrc& max_nodes(size_t val) noexcept {
      _max = val;
      return *this;
    }

    //! \brief Get the current maximum number of nodes.
    //!
    //! This function returns the current value for the maximum number of nodes
    //! in the WordGraph that we are seeking.
    //!
    //! \returns A value of type `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] size_t max_nodes() const noexcept {
      return _max;
    }

    //! \brief Set the target size.
    //!
    //! This function sets the target size, i.e. the size of the
    //! transformation semigroup corresponding to the WordGraph returned by the
    //! function \ref word_graph.
    //!
    //! \param val the target size.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \noexcept
    RepOrc& target_size(size_t val) noexcept {
      _size = val;
      return *this;
    }

    //! \brief Get the current target size.
    //!
    //! This function returns the current value for the target size, i.e. the
    //! size of the transformation semigroup corresponding to the
    //! WordGraph returned by the function \ref word_graph.
    //!
    //! \returns A value of type `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] size_t target_size() const noexcept {
      return _size;
    }

    //! \brief Get the word_graph.
    //!
    //! This function attempts to find a right congruence, represented as an
    //! WordGraph, of the semigroup or monoid defined by the presentation
    //! consisting of its \ref presentation and \ref long_rules with the
    //! following properties:
    //! * the transformation semigroup defined by the WordGraph has size \ref
    //!   target_size;
    //! * the number of nodes in the WordGraph is at least \ref min_nodes
    //!   and at most \ref max_nodes.
    //!
    //! If no such WordGraph can be found, then an empty WordGraph is returned
    //! (with `0` nodes and `0` edges).
    //!
    //! \warning The return value of this function is recomputed every time it
    //! is called.
    //!
    //! \warning If the return value of \ref number_of_threads is greater than
    //! \c 1, then the value returned by this function is non-deterministic, and
    //! may vary even for the same parameters.
    //!
    //! \returns A value of type `WordGraph`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    [[nodiscard]] Sims1::word_graph_type word_graph() const;
  };

  //! \ingroup congruences_group
  //!
  //! \brief For computing the minimal degree of a transformation representation
  //! arising from a right congruences of a finite semigroup or monoid.
  //!
  //! Defined in `sims.hpp`.
  //!
  //! This class is a helper for Sims1, calling the word_graph member
  //! function attempts to find a right congruence, represented as an
  //! WordGraph, with the minimum possible number of nodes such that the
  //! action of the semigroup or monoid defined by the presentation consisting
  //! of its \ref presentation on the nodes of the WordGraph corresponds to a
  //! semigroup of size \ref target_size.
  //!
  //! If no such WordGraph can be found, then an empty WordGraph is
  //! returned (with `0` nodes and `0` edges).
  class MinimalRepOrc : public SimsSettings<MinimalRepOrc> {
   private:
    size_t _size;

   public:
    using SimsSettings<MinimalRepOrc>::stats;

    //! Default constructor
    MinimalRepOrc() : _size() {
      init();
    }

    //! \brief Reinitialize an existing MinimalRepOrc object.
    //!
    //! This function puts a MinimalRepOrc object back into the same state as if
    //! it had been newly default constructed.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    MinimalRepOrc& init() {
      _size = 0;
      return *this;
    }

    //! \brief Set the target size.
    //!
    //! This function sets the target size, i.e. the desired size of the
    //! transformation semigroup corresponding to the WordGraph returned by
    //! the function \ref word_graph.
    //!
    //! \param val the target size.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \noexcept
    MinimalRepOrc& target_size(size_t val) noexcept {
      _size = val;
      return *this;
    }

    //! \brief Get the current target size.
    //!
    //! This function returns the current value for the target size, i.e. the
    //! desired size of the transformation semigroup corresponding to the
    //! WordGraph returned by the function \ref word_graph.
    //!
    //! \returns A value of type `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] size_t target_size() const noexcept {
      return _size;
    }

    //! \brief Get the word graph.
    //!
    //! This function attempts to find a right congruence, represented as an
    //! WordGraph, with the minimum possible number of nodes such that the
    //! action of the semigroup or monoid defined by the presentation
    //! consisting of its \ref presentation and \ref long_rules on the nodes
    //! of the WordGraph corresponds to a semigroup of size \ref
    //! target_size.
    //!
    //! If no such WordGraph can be found, then an empty WordGraph is
    //! returned (with `0` nodes and `0` edges).
    //!
    //! The algorithm implemented by this function repeatedly runs:
    //! \code RepOrc(*this)
    //!     .min_nodes(1)
    //!     .max_nodes(best)
    //!     .target_size(target_size())
    //!     .word_graph();
    //! \endcode
    //! where `best` is initially \ref target_size, until the returned
    //! WordGraph is empty, and then the penultimate WordGraph is returned
    //! (if any).
    //!
    //! \warning The return value of this function is recomputed every time
    //! it is called.
    //!
    //! \warning If the return value of \ref number_of_threads is greater
    //! than \c 1, then the value returned by this function is
    //! non-deterministic, and may vary even for the same parameters.
    //!
    //! \returns A value of type `WordGraph<uint32_t>`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    [[nodiscard]] Sims1::word_graph_type word_graph() const;
  };

  namespace sims {
    class const_cgp_iterator;

    //! \ingroup congruences_group
    //!
    //! \brief For iterating over the right congruence generating pairs.
    //!
    //! Defined in `sims.hpp`.
    //!
    //! This class provides a means of iterating over the right congruence
    //! generating pairs (rcgp) of a finite index right congruence defined by a
    //! \ref WordGraph.
    //!
    //! \sa \ref const_cgp_iterator for the two-sided congruence generating pair
    //! iterator.
    // This is similar to FroidurePinBase::const_rule_iterator
    // Right Congruence Generating Pairs (rcgp)
    class const_rcgp_iterator {
     public:
      using size_type = typename std::vector<relation_type>::size_type;
      using difference_type =
          typename std::vector<relation_type>::difference_type;
      using const_pointer = typename std::vector<relation_type>::const_pointer;
      using pointer       = typename std::vector<relation_type>::pointer;
      using const_reference =
          typename std::vector<relation_type>::const_reference;
      using reference         = typename std::vector<relation_type>::reference;
      using value_type        = relation_type;
      using iterator_category = std::forward_iterator_tag;

      using word_graph_type   = Sims1::word_graph_type;
      using felsch_graph_type = Sims1::felsch_graph_type;
      using node_type         = word_graph_type::node_type;
      using label_type        = word_graph_type::label_type;

     protected:
      felsch_graph_type _reconstructed_word_graph;

     private:
      label_type             _gen;
      node_type              _source;
      mutable relation_type  _relation;
      Forest                 _tree;
      word_graph_type const* _word_graph;

      // Set source to ptr->number_of_active_nodes() for cend
      const_rcgp_iterator(Presentation<word_type> const& p,
                          word_graph_type const*         ptr,
                          node_type                      source,
                          label_type                     gen);

      // To allow the use of the above constructor
      template <typename Node>
      friend const_rcgp_iterator
      cbegin_right_generating_pairs_no_checks(WordGraph<Node> const&);

      template <typename Node>
      friend const_rcgp_iterator
      cbegin_right_generating_pairs_no_checks(Presentation<word_type> const&,
                                              WordGraph<Node> const&);

      // To allow the use of the above constructor
      template <typename Node>
      friend const_rcgp_iterator
      cend_right_generating_pairs_no_checks(WordGraph<Node> const&);

      template <typename Node>
      friend const_rcgp_iterator
      cend_right_generating_pairs_no_checks(Presentation<word_type> const&,
                                            WordGraph<Node> const&);

      // To allow the use of the private constructor from pointer
      template <typename Node>
      friend const_cgp_iterator
      cbegin_two_sided_generating_pairs_no_checks(WordGraph<Node> const&);

      template <typename Node>
      friend const_cgp_iterator cbegin_two_sided_generating_pairs_no_checks(
          Presentation<word_type> const&,
          WordGraph<Node> const&);

      // To allow the use of the above constructor
      template <typename Node>
      friend const_cgp_iterator
      cend_two_sided_generating_pairs_no_checks(WordGraph<Node> const&);

      template <typename Node>
      friend const_cgp_iterator
      cend_two_sided_generating_pairs_no_checks(Presentation<word_type> const&,
                                                WordGraph<Node> const&);

     public:
      //! Default constructor
      const_rcgp_iterator() = default;
      //! Default copy constructor.
      const_rcgp_iterator(const_rcgp_iterator const&) = default;
      //! Default move constructor.
      const_rcgp_iterator(const_rcgp_iterator&&) = default;
      //! Default copy assignment operator.
      const_rcgp_iterator& operator=(const_rcgp_iterator const&) = default;
      //! Default move assignment operator.
      const_rcgp_iterator& operator=(const_rcgp_iterator&&) = default;

      ~const_rcgp_iterator();

      //! Check if both iterators point to the same generating pair
      [[nodiscard]] bool
      operator==(const_rcgp_iterator const& that) const noexcept {
        return _gen == that._gen && _source == that._source;
      }

      //! Check if both iterators point to distinct generating pairs
      [[nodiscard]] bool
      operator!=(const_rcgp_iterator const& that) const noexcept {
        return !(this->operator==(that));
      }

      //! Return a const reference to the generating pair pointed to by the
      //! iterator
      [[nodiscard]] const_reference operator*() const {
        populate_relation();
        return _relation;
      }

      //! Return a const pointer to the generating pair pointed to by the
      //! iterator
      [[nodiscard]] const_pointer operator->() const {
        populate_relation();
        return &_relation;
      }

      //! No doc
      // prefix
      const_rcgp_iterator const& operator++();

      //! No doc
      // postfix
      const_rcgp_iterator operator++(int) noexcept {
        return detail::default_postfix_increment(*this);
      }

      //! No doc
      void swap(const_rcgp_iterator& that) noexcept;

     protected:
      [[nodiscard]] bool at_end() const noexcept {
        return _source == _word_graph->number_of_active_nodes();
      }

     private:
      bool populate_relation() const;
    };  // const_rcgp_iterator

    //! \ingroup congruences_group
    //!
    //! \brief For iterating over the two-sided congruence generating pairs.
    //!
    //! Defined in `sims.hpp`.
    //!
    //! This class provides a means of iterating over the two-sided congruence
    //! generating pairs (cgp) of a finite index two-sided congruence defined by
    //! a \ref WordGraph.
    //!
    //! \sa \ref const_rcgp_iterator for the right congruence generating
    //! pair iterator.
    class const_cgp_iterator : public const_rcgp_iterator {
     public:
      using size_type         = const_rcgp_iterator::size_type;
      using difference_type   = const_rcgp_iterator::difference_type;
      using const_pointer     = const_rcgp_iterator::const_pointer;
      using pointer           = const_rcgp_iterator::pointer;
      using const_reference   = const_rcgp_iterator::const_reference;
      using reference         = const_rcgp_iterator::reference;
      using value_type        = relation_type;
      using iterator_category = std::forward_iterator_tag;

      using word_graph_type   = Sims1::word_graph_type;
      using felsch_graph_type = Sims1::felsch_graph_type;
      using node_type         = word_graph_type::node_type;
      using label_type        = word_graph_type::label_type;

      using const_rcgp_iterator::const_rcgp_iterator;

      ~const_cgp_iterator();

      using const_rcgp_iterator::operator==;
      using const_rcgp_iterator::operator!=;
      using const_rcgp_iterator::operator*;
      using const_rcgp_iterator::operator->;

      //! No doc
      // prefix
      const_cgp_iterator const& operator++();

      //! No doc
      // postfix
      const_cgp_iterator operator++(int) noexcept {
        return detail::default_postfix_increment(*this);
      }

      //! No doc
      void swap(const_cgp_iterator& that) noexcept {
        const_rcgp_iterator::swap(that);
      }

     private:
      bool populate_relation() const;
    };  // const_cgp_iterator

    // Forward decl, documented below
    template <typename Node>
    rx::iterator_range<const_rcgp_iterator>
    right_generating_pairs(Presentation<word_type> const& p,
                           WordGraph<Node> const&         wg);

    // Forward decl, documented below
    template <typename Node>
    rx::iterator_range<const_rcgp_iterator>
    right_generating_pairs(WordGraph<Node> const& wg);

    // Forward decl, documented below
    template <typename Node>
    rx::iterator_range<const_cgp_iterator>
    two_sided_generating_pairs(Presentation<word_type> const& p,
                               WordGraph<Node> const&         wg);

    // Forward decl, documented below
    template <typename Node>
    rx::iterator_range<const_cgp_iterator>
    two_sided_generating_pairs(WordGraph<Node> const& wg);

    //! \ingroup congruences_group
    //!
    //! \brief Check if a word graph defines a right congruence on an
    //! f.p. semigroup or monoid.
    //!
    //! \returns `true` if the word graph \p wg defines a right congruence
    //! on the semigroup or monoid defined by \p p and `false` otherwise.
    template <typename Node>
    bool is_right_congruence(Presentation<word_type> const& p,
                             WordGraph<Node> const&         wg);

    //! No doc
    template <typename Node>
    void throw_if_not_right_congruence(Presentation<word_type> const& p,
                                       WordGraph<Node> const&         wg);

    //! \ingroup congruences_group
    //!
    //! \brief Check if a word graph defines a right congruence on the dual of
    //! an f.p. semigroup or monoid.
    //!
    //! This function checks if \p wg defines a right congruence over the dual
    //! of \p p. This is equivalent to checking if the word graph defines a
    //! left congruence in the semigroup or monoid defined by \p p.
    //!
    //! \returns `true` if the word graph \p wg defines a right congruence
    //! on the dual of the semigroup or monoid defined by \p p and `false`
    //! otherwise.
    template <typename Node>
    bool is_right_congruence_of_dual(Presentation<word_type> const& p,
                                     WordGraph<Node> const&         wg) {
      auto p_rev(p);
      presentation::reverse(p_rev);
      return is_right_congruence(p_rev, wg);
    }

    //! \ingroup congruences_group
    //!
    //! \brief Check if a word graph defines a two-sided congruence on
    //! an f.p. semigroup or monoid (no checks).
    //!
    //! \returns `true` if the word graph \p wg defines a two-sided
    //! congruence on the semigroup or monoid defined by \p p and `false`
    //! otherwise.
    //!
    //! \warning
    //! This function does no checks on its arguments whatsoever and assumes
    //! that \p wg already defines a two-sided congruence on the semigroup or
    //! monoid defined by \p p. If this is not the case, then bad things may
    //! happen.
    template <typename Node>
    bool is_two_sided_congruence_no_checks(Presentation<word_type> const& p,
                                           WordGraph<Node> const&         wg);

    //! \ingroup congruences_group
    //!
    //! \brief Check if a word graph defines a two sided congruence on
    //! an f.p. semigroup or monoid.
    //!
    //! \returns `true` if the word graph \p wg defines a two-sided
    //! congruence on the semigroup or monoid defined by \p p and `false`
    //! otherwise.
    template <typename Node>
    bool is_two_sided_congruence(Presentation<word_type> const& p,
                                 WordGraph<Node> const&         wg);

    //! No doc
    template <typename Node>
    void throw_if_not_two_sided_congruence(Presentation<word_type> const& p,
                                           WordGraph<Node> const&         wg);

    //! \ingroup congruences_group
    //!
    //! \brief Get an iterator pointing to the first right congruence generating
    //! pair (no checks).
    //!
    //! \returns an iterator pointing to the first generating pair of the right
    //! congruence defined by the word graph \p wg on the semigroup or monoid
    //! defined by \p p.
    //!
    //! \warning
    //! This function does no checks on its arguments whatsoever and assumes
    //! that \p wg already defines a right congruence on the semigroup or monoid
    //! defined by \p p. If this is not the case, then bad things may happen.
    // Note that this is the generating pairs of the right
    // congruence so defined not the two-sided congruence.
    template <typename Node>
    const_rcgp_iterator
    cbegin_right_generating_pairs_no_checks(Presentation<word_type> const& p,
                                            WordGraph<Node> const&         wg) {
      return const_rcgp_iterator(p, &wg, 0, 0);
    }

    //! \ingroup congruences_group
    //!
    //! \brief Get an iterator pointing to the first right congruence generating
    //! pair.
    //!
    //! \returns an iterator pointing to the first generating pair of the right
    //! congruence defined by the word graph \p wg on the semigroup or monoid
    //! defined by \p p.
    //!
    //! \throws LibsemigroupsException if the argument \p wg does not define a
    //! right congruence on the monoid presented by \p p.
    template <typename Node>
    const_rcgp_iterator
    cbegin_right_generating_pairs(Presentation<word_type> const& p,
                                  WordGraph<Node> const&         wg) {
      throw_if_not_right_congruence(p, wg);
      return cbegin_right_generating_pairs_no_checks(p, wg);
    }

    //! \ingroup congruences_group
    //!
    //! \brief Get an iterator pointing to the first two-sided congruence
    //! generating pair (no checks).
    //!
    //! \returns an iterator pointing to the first generating pair of the
    //! two-sided congruence defined by the word graph \p wg on the semigroup or
    //! monoid defined by \p p.
    //!
    //! \note
    //! The generating pairs of a two-sided congruence \f$\rho\f$ as a two-sided
    //! congruence differ from the generating pairs of \f$\rho\f$ as a right
    //! congruence. This function returns the two-sided congruence generating
    //! pairs, not the right congruence generating pairs of a two-sided
    //! congruence.
    //!
    //! \warning
    //! This function does no checks on its arguments whatsoever and assumes
    //! that \p wg already defines a two-sided congruence on the semigroup or
    //! monoid defined by \p p. If this is not the case, then bad things may
    //! happen.
    // Note that this is the generating pairs of the right
    // congruence so defined not the two-sided congruence.
    template <typename Node>
    const_cgp_iterator cbegin_two_sided_generating_pairs_no_checks(
        Presentation<word_type> const& p,
        WordGraph<Node> const&         wg) {
      return const_cgp_iterator(p, &wg, 0, 0);
    }

    //! \ingroup congruences_group
    //!
    //! \brief Get an iterator pointing to the first two-sided congruence
    //! generating pair.
    //!
    //! \returns an iterator pointing to the first generating pair of the
    //! two-sided congruence defined by the word graph \p wg on the semigroup or
    //! monoid defined by \p p.
    //!
    //! \throws LibsemigroupsException if the argument \p wg does not define a
    //! two-sided congruence on the semigroup or monoid defined by \p p.
    //!
    //! \note
    //! The generating pairs of a two-sided congruence \f$\rho\f$ as a two-sided
    //! congruence differ from the generating pairs of \f$\rho\f$ as a right
    //! congruence. This function returns the two-sided congruence generating
    //! pairs, not the right congruence generating pairs of a two-sided
    //! congruence.
    template <typename Node>
    const_cgp_iterator
    cbegin_two_sided_generating_pairs(Presentation<word_type> const& p,
                                      WordGraph<Node> const&         wg) {
      throw_if_not_two_sided_congruence(p, wg);
      return cbegin_two_sided_generating_pairs_no_checks(p, wg);
    }

    //! \ingroup congruences_group
    //!
    //! \brief Get an iterator pointing to one after the last right congruence
    //! generating pair (no checks).
    //!
    //! \returns an iterator pointing to one after the last generating pair of
    //! the right congruence defined by the word graph \p wg on the semigroup or
    //! monoid defined by \p p.
    //!
    //! \warning
    //! This function does no checks on its arguments whatsoever and assumes
    //! that \p wg already defines a right congruence on the semigroup or
    //! monoid defined by \p p. If this is not the case, then bad things may
    //! happen.
    template <typename Node>
    const_rcgp_iterator
    cend_right_generating_pairs_no_checks(Presentation<word_type> const& p,
                                          WordGraph<Node> const&         wg) {
      return const_rcgp_iterator(p, &wg, wg.number_of_active_nodes(), 0);
    }

    //! \ingroup congruences_group
    //!
    //! \brief Get an iterator pointing to one after the last right congruence
    //! generating pair.
    //!
    //! \returns an iterator pointing to one after the last generating pair of
    //! the right congruence defined by the word graph \p wg on the semigroup or
    //! monoid defined by \p p.
    //!
    //! \throws LibsemigroupsException if the argument \p wg does not define a
    //! right congruence on the semigroup or monoid defined by \p p.
    template <typename Node>
    const_rcgp_iterator
    cend_right_generating_pairs(Presentation<word_type> const& p,
                                WordGraph<Node> const&         wg) {
      throw_if_not_right_congruence(p, wg);
      return cend_right_generating_pairs_no_checks(p, wg);
    }

    //! \ingroup congruences_group
    //!
    //! \brief Get an iterator pointing to one after the last two-sided
    //! congruence generating pair (no checks).
    //!
    //! \returns an iterator pointing to one after the last generating pair of
    //! the two-sided congruence defined by the word graph \p wg on the
    //! semigroup or monoid defined by \p p.
    //!
    //! \note
    //! The generating pairs of a two-sided congruence \f$\rho\f$ as a two-sided
    //! congruence differ from the generating pairs of \f$\rho\f$ as a right
    //! congruence. This function returns the two-sided congruence generating
    //! pairs, not the right congruence generating pairs of a two-sided
    //! congruence.
    //!
    //! \warning
    //! This function does no checks on its arguments whatsoever and assumes
    //! that \p wg already defines a two-sided congruence on the semigroup or
    //! monoid defined by \p p. If this is not the case, then bad things may
    //! happen.
    template <typename Node>
    const_cgp_iterator
    cend_two_sided_generating_pairs_no_checks(Presentation<word_type> const& p,
                                              WordGraph<Node> const& wg) {
      return const_cgp_iterator(p, &wg, wg.number_of_active_nodes(), 0);
    }

    //! \ingroup congruences_group
    //!
    //! \brief Get an iterator pointing to one after the last two-sided
    //! congruence generating pair.
    //!
    //! \returns an iterator pointing to one after the last generating pair of
    //! the two-sided congruence defined by the word graph \p wg on the
    //! semigroup or monoid defined by \p p.
    //!
    //! \note
    //! The generating pairs of a two-sided congruence \f$\rho\f$ as a two-sided
    //! congruence differ from the generating pairs of \f$\rho\f$ as a right
    //! congruence. This function returns the two-sided congruence generating
    //! pairs, not the right congruence generating pairs of a two-sided
    //! congruence.
    //!
    //! \throws LibsemigroupsException if the argument \p wg does not define a
    //! two-sided congruence on the semigroup or monoid defined by \p p.
    template <typename Node>
    const_cgp_iterator
    cend_two_sided_generating_pairs(Presentation<word_type> const& p,
                                    WordGraph<Node> const&         wg) {
      throw_if_not_two_sided_congruence(p, wg);
      return cend_two_sided_generating_pairs_no_checks(p, wg);
    }

    //! \ingroup congruences_group
    //!
    //! \brief Compute the right congruence generating pairs of a word graph on
    //! an f.p. semigroup or monoid (no checks).
    //!
    //! \returns the right congruence generating pairs of the right
    //! congruence defined by the
    //! word graph \p wg on the semigroup or monoid defined by \p p.
    //!
    //! \warning
    //! This function does no checks on its arguments whatsoever and assumes
    //! that \p wg already defines a right congruence on the semigroup or monoid
    //! defined by \p p. If this is not the case, then bad things may happen.
    template <typename Node>
    rx::iterator_range<const_rcgp_iterator>
    right_generating_pairs_no_checks(Presentation<word_type> const& p,
                                     WordGraph<Node> const&         wg) {
      return rx::iterator_range(cbegin_right_generating_pairs_no_checks(p, wg),
                                cend_right_generating_pairs_no_checks(p, wg));
    }

    //! \ingroup congruences_group
    //!
    //! \brief Compute the right congruence generating pairs of a word graph on
    //! an f.p. semigroup or monoid.
    //!
    //! \returns the right congruence generating pairs of the right
    //! congruence defined by the
    //! word graph \p wg on the semigroup or monoid defined by \p p.
    //!
    //! \throws LibsemigroupsException if the argument \p wg does not define a
    //! right congruence on the semigroup or monoid defined by \p p.
    template <typename Node>
    rx::iterator_range<const_rcgp_iterator>
    right_generating_pairs(Presentation<word_type> const& p,
                           WordGraph<Node> const&         wg) {
      throw_if_not_right_congruence(p, wg);
      return right_generating_pairs_no_checks(p, wg);
    }

    //! \ingroup congruences_group
    //!
    //! \brief Compute the right congruence generating pairs of a word graph on
    //! the free monoid (no checks).
    //!
    //! \returns the right congruence generating pairs of the right
    //! congruence defined by the
    //! word graph \p wg on the free monoid.
    //!
    //! \warning
    //! This function does no checks on its arguments whatsoever and assumes
    //! that \p wg already defines a right congruence on the free monoid. If
    //! this is not the case, then bad things may happen.
    template <typename Node>
    rx::iterator_range<const_rcgp_iterator>
    right_generating_pairs_no_checks(WordGraph<Node> const& wg) {
      Presentation<word_type> p;
      p.alphabet(wg.out_degree());
      return right_generating_pairs_no_checks(p, wg);
    }

    //! \ingroup congruences_group
    //!
    //! \brief Compute the right congruence generating pairs of a word graph on
    //! the free monoid.
    //!
    //! \returns the right congruence generating pairs of the right
    //! congruence defined by the
    //! word graph \p wg on the free monoid.
    //!
    //! \throws LibsemigroupsException if the argument \p wg does not define a
    //! right congruence on the free monoid.
    template <typename Node>
    rx::iterator_range<const_rcgp_iterator>
    right_generating_pairs(WordGraph<Node> const& wg);

    //! \ingroup congruences_group
    //!
    //! \brief Compute the two-sided congruence generating pairs of a word graph
    //! on an f.p. semigroup or monoid (no checks).
    //!
    //! \returns the two-sided congruence generating pairs of the two-sided
    //! congruence defined by the
    //! word graph \p wg on the semigroup or monoid defined by \p p.
    //!
    //! \note
    //! The generating pairs of a two-sided congruence \f$\rho\f$ as a two-sided
    //! congruence differ from the generating pairs of \f$\rho\f$ as a right
    //! congruence. This function returns the two-sided congruence generating
    //! pairs, not the right congruence generating pairs of a two-sided
    //! congruence.
    //!
    //! \warning
    //! This function does no checks on its arguments whatsoever and assumes
    //! that \p wg already defines a two-sided congruence on the semigroup or
    //! monoid defined by \p p. If this is not the case, then bad things may
    //! happen.
    template <typename Node>
    rx::iterator_range<const_cgp_iterator>
    two_sided_generating_pairs_no_checks(Presentation<word_type> const& p,
                                         WordGraph<Node> const&         wg) {
      return rx::iterator_range(
          cbegin_two_sided_generating_pairs_no_checks(p, wg),
          cend_two_sided_generating_pairs_no_checks(p, wg));
    }

    //! \ingroup congruences_group
    //!
    //! \brief Compute the two-sided congruence generating pairs of a word graph
    //! on an f.p. semigroup or monoid.
    //!
    //! \returns the two-sided congruence generating pairs of the two-sided
    //! congruence defined by the
    //! word graph \p wg on the semigroup or monoid defined by \p p.
    //!
    //! \note
    //! The generating pairs of a two-sided congruence \f$\rho\f$ as a two-sided
    //! congruence differ from the generating pairs of \f$\rho\f$ as a right
    //! congruence. This function returns the two-sided congruence generating
    //! pairs, not the right congruence generating pairs of a two-sided
    //! congruence.
    //!
    //! \throws LibsemigroupsException if the argument \p wg does not define a
    //! two-sided congruence on the semigroup or monoid presented by \p p.
    template <typename Node>
    rx::iterator_range<const_cgp_iterator>
    two_sided_generating_pairs(Presentation<word_type> const& p,
                               WordGraph<Node> const&         wg) {
      throw_if_not_two_sided_congruence(p, wg);
      return two_sided_generating_pairs_no_checks(p, wg);
    }

    //! \ingroup congruences_group
    //!
    //! \brief Compute the two-sided congruence generating pairs of a word graph
    //! on the free monoid (no checks).
    //!
    //! \returns the two-sided congruence generating pairs of the two-sided
    //! congruence defined by the
    //! word graph \p wg on the free monoid.
    //!
    //! \note
    //! The generating pairs of a two-sided congruence \f$\rho\f$ as a two-sided
    //! congruence differ from the generating pairs of \f$\rho\f$ as a right
    //! congruence. This function returns the two-sided congruence generating
    //! pairs, not the right congruence generating pairs of a two-sided
    //! congruence.
    //!
    //! \warning
    //! This function does no checks on its arguments whatsoever and assumes
    //! that \p wg already defines a two-sided congruence on the free monoid. If
    //! this is not the case, then bad things may happen.
    template <typename Node>
    rx::iterator_range<const_cgp_iterator>
    two_sided_generating_pairs_no_checks(WordGraph<Node> const& wg) {
      Presentation<word_type> p;
      p.alphabet(wg.out_degree());
      return rx::iterator_range(
          cbegin_two_sided_generating_pairs_no_checks(p, wg),
          cend_two_sided_generating_pairs_no_checks(p, wg));
    }

    //! \ingroup congruences_group
    //!
    //! \brief Compute the two-sided congruence generating pairs of a word graph
    //! on the free monoid.
    //!
    //! \returns the two-sided congruence generating pairs of the two-sided
    //! congruence defined by the
    //! word graph \p wg on the free monoid.
    //!
    //! \note
    //! The generating pairs of a two-sided congruence \f$\rho\f$ as a two-sided
    //! congruence differ from the generating pairs of \f$\rho\f$ as a right
    //! congruence. This function returns the two-sided congruence generating
    //! pairs, not the right congruence generating pairs of a two-sided
    //! congruence.
    //!
    //! \throws LibsemigroupsException if the argument \p wg does not define a
    //! two-sided congruence on the free monoid.
    template <typename Node>
    rx::iterator_range<const_cgp_iterator>
    two_sided_generating_pairs(WordGraph<Node> const& wg) {
      Presentation<word_type> p;
      p.alphabet(wg.out_degree());
      throw_if_not_two_sided_congruence(p, wg);
      return two_sided_generating_pairs_no_checks(p, wg);
    }

    //! \ingroup congruences_group
    //!
    //! \brief Check if a word graph defines a maximal right congruence on
    //! an f.p. semigroup or monoid.
    //!
    //! \returns `true` if the word graph \p wg defines a maximal right
    //! congruence on the semigroup or monoid defined by \p p and `false`
    //! otherwise.
    template <typename Node>
    bool is_maximal_right_congruence(Presentation<word_type> const& p,
                                     WordGraph<Node> const&         wg);

    //! \ingroup congruences_group
    //!
    //! \brief Compute the inclusion poset of a collection of congruences
    //! defined by word graphs.
    //!
    //! This function computes a boolean matrix whose `(i, j)`-th entry is 1 if
    //! and only if the congruence defined by the `i`-th word graph in `[first,
    //! last)` is a subrelation of the congruence defined by the `j`-th word
    //! graph in `[first, last)`. If `[first, last)` is the collection of all
    //! congruences of a given semigroup or monoid, then this is equivalent to
    //! computing the congruence lattice of the semigroup or monoid.
    //!
    //! \returns a \ref Bmat object.
    //!
    //! \warning
    //! This function does no checks on its arguments whatsoever and assumes
    //! that every element in the range `[first, last)` is a complete word
    //! graph. If this is not the case, then bad things may happen.
    template <typename Iterator>
    BMat<> poset(Iterator first, Iterator last);

    //! \ingroup congruences_group
    //!
    //! \brief Construct a \ref Dot object representing the inclusion poset of a
    //! collection of word graphs.
    //!
    //! \returns a \ref Dot object.
    // The following produces a self-contained Dot object which doesn't render
    // very well.
    template <typename Iterator>
    Dot dot_poset(Iterator first, Iterator last);

  }  // namespace sims

  //! \ingroup congruences_group
  //!
  //! \brief For pruning the search tree when looking for congruences arising
  //! from right or two-sided congruences representing faithful actions.
  //!
  //! Defined in `sims.hpp`.
  //!
  //! This class provides a pruner for pruning the search tree when looking for
  //! right congruences representing faithful actions. A right congruence
  //! represents a faithful action if and only if it does not contain any
  //! non-trivial two-sided congruence. Equivalently, a word graph of a right
  //! congruence represents a faithful action if and only if there is no
  //! nontrivial pair of words \f$(u, v)\f$ such that every vertex of the word
  //! graph is compatible with \f$(u, v)\f$.
  //!
  //! \sa \ref SimsSettings::pruners
  //! \sa \ref SimsSettings::add_pruner
  // This struct provides an alternative way of doing MinimalRepOrc, when the
  // generating pairs of the minimal 2-sided congruences are known. These pairs
  // should be added to forbid, and then your SimsRefinerFaithful instance
  // should be passed to a Sims1 object via add_pruner.
  class SimsRefinerFaithful {
   private:
    std::vector<word_type> _forbid;

   public:
    //! Default constructor.
    explicit SimsRefinerFaithful() : _forbid() {}

    //! \brief Reinitialize an existing SimsRefinerFaithful object.
    //!
    //! This function puts an object back into the same state as if it had
    //! been newly default constructed.
    //!
    //! \returns A reference to \c *this.
    SimsRefinerFaithful& init() {
      _forbid.clear();
      return *this;
    }

    //! \brief Construct from set of forbidden pairs.
    //!
    //! Constructs a SimsRefinerFaithful pruner with respect to the set of
    //! forbidden relations in \p forbid.
    //!
    //! If \p forbid contains no trivial pairs (i.e. pairs of words that are
    //! equal in the underlying semigroup or monoid), then all word graphs
    //! rejected by SimsRefinerFaithful are guaranteed to not be extendable to a
    //! word graph representing a faithful congruence. Otherwise, the pruner
    //! will incorrectly reject all word graphs.
    //!
    //! If in addition \p forbid is a set of relations
    //! containing all minimal congruence generating pairs of a given semigroup
    //! or monoid, then SimsRefinerFaithful will also correctly determine if a
    //! complete word graph represents a faithful congruence. Otherwise, the
    //! complete word graphs accepted by SimsRefinerFaithful are not guaranteed
    //! to be faithful and must be checked by some other means.
    //!
    //! \warning
    //! This method does not verify if \p forbid contains trivial pairs or not.
    // TODO(2): Construct from std::vector<std::string>
    // TODO(0) construct from iterators, and move word_type version of to helper
    // namespace
    explicit SimsRefinerFaithful(std::vector<word_type> const& forbid)
        : _forbid(forbid) {}

    //! \brief Reinitialize an existing SimsRefinerFaithful object from a set of
    //! forbidden pairs.
    //!
    //! This function puts an object back into the same state as if it had
    //! been newly constructed from set of forbidden pairs \p forbid.
    //!
    //! \returns A reference to \c *this.
    //!
    //! \warning
    //! This method does not verify if \p forbid contains trivial pairs or not.
    // TODO(0) construct from iterators, and move word_type version of to helper
    // namespace
    SimsRefinerFaithful& init(std::vector<word_type> const& forbid) {
      _forbid.clear();
      _forbid = forbid;
      return *this;
    }

    //! \brief Get the forbidden pairs defining the refiner.
    //!
    //! \anchor forbid
    //! Returns a const reference to the current forbidden pairs.
    //!
    //! This function returns the defining forbidden pairs of a
    //! SimsRefinerFaithful instance.
    //!
    //! \returns A const reference to `std::vector<word_type>`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] std::vector<word_type> const& forbid() const noexcept {
      return _forbid;
    }

    //! \brief Check if a word graph can be extended to one defining a faithful
    //! congruence.
    //!
    //! Returns `false` if there is no way of adding edges and nodes to \par wg
    //! which will result in a word graph defining a faithful congruence.
    //! Otherwise returns `true`.
    //!
    bool operator()(Sims1::word_graph_type const& wg);
  };

  //! \ingroup congruences_group
  //!
  //! \brief For pruning the search tree when looking for congruences arising
  //! from right or two-sided ideals.
  //!
  //! Defined in `sims.hpp`.
  //!
  //! This class provides a pruner for pruning the search tree when looking for
  //! right or two-sided congruences arising from ideals (Rees congruences).
  //!
  //! \sa \ref SimsSettings::pruners
  //! \sa \ref SimsSettings::add_pruner
  class SimsRefinerIdeals {
   private:
    using node_type    = uint32_t;
    using KnuthBendix_ = KnuthBendix<>;
    std::vector<KnuthBendix_> _knuth_bendices;
    Presentation<word_type>   _presentation;

   public:
    //! Default constructor.
    explicit SimsRefinerIdeals()
        : _knuth_bendices(std::thread::hardware_concurrency() + 1,
                          KnuthBendix_()),
          _presentation() {
      init();
    }

    //! \brief Reinitialize an existing SimsRefinerIdeals object.
    //!
    //! This function puts an object back into the same state as if it had
    //! been newly default constructed.
    //!
    //! \returns A reference to \c *this.
    SimsRefinerIdeals& init() {
      _presentation.init();
      _knuth_bendices[0].init();
      std::fill(_knuth_bendices.begin() + 1,
                _knuth_bendices.end(),
                _knuth_bendices[0]);
      return *this;
    }

    //! \brief Construct from presentation.
    //!
    //! Constructs a SimsRefinerIdeals pruner for the semigroup or monoid
    //! defined by \p p.
    //!
    //! \tparam Word the type of words in the input presentation.
    //!
    //! \param p the presentation.
    //!
    //! \warning
    //! This method assumes that KnuthBendix terminates on the input
    //! presentation \p p. If this is not the case then th pruner may not
    //! terminate on certain inputs.
    template <typename Word>
    explicit SimsRefinerIdeals(Presentation<Word> const& p)
        : _knuth_bendices(std::thread::hardware_concurrency() + 1,
                          KnuthBendix_()),
          _presentation() {
      init(p);
    }

    //! \brief Reinitialize an existing SimsRefinerIdeals object from a
    //! word_type presentation.
    //!
    //! This function puts an object back into the same state as if it had
    //! been newly constructed from the presentation \p p.
    //!
    //! \tparam Word the type of words in the input presentation.
    //!
    //! \param p the presentation.
    //!
    //! \returns A reference to \c *this.
    //!
    //! \throws LibsemigroupsException if `p` is not valid
    //! \throws LibsemigroupsException if `p` has 0-generators and 0-relations.
    //!
    //! \warning This function has no exception guarantee, the object will be
    //! in the same state as if it was default constructed if an exception is
    //! thrown.
    //!
    //! \warning
    //! This method assumes that KnuthBendix terminates on the input
    //! presentation \p p. If this is not the case then th pruner may not
    //! terminate on certain inputs.
    //!
    //! \sa presentation(Presentation<word_type> const&)
    template <typename Word>
    SimsRefinerIdeals& init(Presentation<Word> const& p);

    //! \brief Get the presentation over which the refiner is defined.
    //!
    //! \anchor presentation
    //! Returns a const reference to the current relations of the underlying
    //! presentation.
    //!
    //! This function returns the defining presentation of a SimsRefinerIdeals
    //! instance.
    //!
    //! \returns A const reference to `Presentation<word_type>`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] Presentation<word_type> const& presentation() const noexcept {
      return _presentation;
    }

    //! \brief Check if a word graph can be extended to one defining a Rees
    //! congruence.
    //!
    //! Returns `false` if there is no way of adding edges and nodes to \par wg
    //! which will result in a word graph defining a Rees congruence. Otherwise
    //! returns `true`.
    //!
    //! \warning
    //! This method assumes that KnuthBendix terminates on the underlying
    //! presentation that was used to construct the SimsRefinerIdeals object. If
    //! this is not the case then th pruner may not terminate on certain inputs.
    bool operator()(Sims1::word_graph_type const& wg);
  };

  //! \ingroup congruence_group
  //!
  //! \brief Return a human readable representation of a SimsStats object.
  //!
  //! Return a human readable representation of a SimsStats object.
  //!
  //! \param x the SimsStats object.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  [[nodiscard]] std::string to_human_readable_repr(SimsStats const& x);

  //! \ingroup congruence_group
  //!
  //! \brief Return a human readable representation of a Sims1 object.
  //!
  //! Return a human readable representation of a Sims1 object.
  //!
  //! \param x the Sims1 object.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  [[nodiscard]] std::string to_human_readable_repr(Sims1 const& x);

  //! \ingroup congruence_group
  //!
  //! \brief Return a human readable representation of a Sims2 object.
  //!
  //! Return a human readable representation of a Sims2 object.
  //!
  //! \param x the Sims2 object.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  [[nodiscard]] std::string to_human_readable_repr(Sims2 const& x);

  //! \ingroup congruence_group
  //!
  //! \brief Return a human readable representation of a RepOrc object.
  //!
  //! Return a human readable representation of a RepOrc object.
  //!
  //! \param x the RepOrc object.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  [[nodiscard]] std::string to_human_readable_repr(RepOrc const& x);

  //! \ingroup congruence_group
  //!
  //! \brief Return a human readable representation of a MinimalRepOrc object.
  //!
  //! Return a human readable representation of a MinimalRepOrc object.
  //!
  //! \param x the MinimalRepOrc object.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  [[nodiscard]] std::string to_human_readable_repr(MinimalRepOrc const& x);

  //! \ingroup congruence_group
  //!
  //! \brief Return a human readable representation of a SimsRefinerIdeals
  //! object.
  //!
  //! Return a human readable representation of a SimsRefinerIdeals object.
  //!
  //! \param x the SimsRefinerIdeals object.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  [[nodiscard]] std::string to_human_readable_repr(SimsRefinerIdeals const& x);

  //! \ingroup congruence_group
  //!
  //! \brief Return a human readable representation of a SimsRefinerFaithful
  //! object.
  //!
  //! Return a human readable representation of a SimsRefinerFaithful object.
  //!
  //! \param x the SimsRefinerIdeals object.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  [[nodiscard]] std::string
  to_human_readable_repr(SimsRefinerFaithful const& x);

}  // namespace libsemigroups

#include "sims.tpp"

#endif  // LIBSEMIGROUPS_SIMS_HPP_
