//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2022 James D. Mitchell
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
// congruence" algorithm for 1-sided congruences of semigroups and monoids.

// TODO(Sims1):
// * implement joins (HopcroftKarp), meets (not sure), containment (find join
//   and check equality)?
// * generating pairs for congruences defined by "word_graph"?
// * is 2-sided congruence method. One approach would be to compute the kernel
//   of the associated homomorphism, which is the largest 2-sided congruence
//   contained in the right congruence. Not sure if this is a good approach.
// * a version which allows specifying the word_graph to Sims1 too
//
// Notes:
// 1. In 2022, when first writing this file, JDM tried templating the word_type
//    used by the presentations in Sims1 (so that we could use StaticVector1
//    for example, using smaller integer types for letters, and the stack to
//    hold the words rather than the heap), but this didn't seem to give any
//    performance improvement, and so I backed out the changes.

#ifndef LIBSEMIGROUPS_SIMS_HPP_
#define LIBSEMIGROUPS_SIMS_HPP_

#include <algorithm>   // for max
#include <cstddef>     // for size_t
#include <cstdint>     // for uint64_t, uint32_t
#include <functional>  // for function
#include <iterator>    // for forward_iterator_tag
#include <mutex>       // for mutex
#include <string>      // for operator+, basic_string
#include <thread>      // for thread, yield
#include <utility>     // for move
#include <vector>      // for vector

#include "debug.hpp"            // for LIBSEMIGROUPS_ASSERT
#include "exception.hpp"        // for LIBSEMIGROUPS_EXCEPTION
#include "felsch-graph.hpp"     // for FelschGraph
#include "presentation.hpp"     // for Presentation, Presentati...
#include "to-presentation.hpp"  // for to_presentation
#include "types.hpp"            // for word_type, congruence_kind
#include "word-graph.hpp"       // for WordGraph

namespace libsemigroups {

  //! Defined in ``sims.hpp``.
  //!
  //! On this page we describe the `SimsStats` class. The purpose of this
  //! class is to collect some statistics related to `Sims1` class template.
  //!
  //! \sa \ref Sims1
  class SimsStats {
   public:
    // TODO(doc)
    // Not atomic because this is only accessed by report_progress_from_thread
    uint64_t count_last;
    // TODO(doc)
    // Atomic so as to avoid races between report_progress_from_thread and the
    // threads modifying count_last
    std::atomic_uint64_t count_now;

    //! The maximum number of pending definitions.
    //!
    //! A *pending definition* is just an edge that will be defined at some
    //! point in the future in the WordGraph represented by a Sims1
    //! instance at any given moment.
    //!
    //! This member tracks the maximum number of such pending definitions that
    //! occur during the running of the algorithms in Sims1.
    std::atomic_uint64_t max_pending;

    //! The total number of pending definitions.
    //!
    //! A *pending definition* is just an edge that will be defined at some
    //! point in the future in the WordGraph represented by a Sims1
    //! instance at any given moment.
    //!
    //! This member tracks the total number of pending definitions that
    //! occur during the running of the algorithms in Sims1. This is the same
    //! as the number of nodes in the search tree encounter during the running
    //! of Sims1.
    // Not atomic because this is only accessed by report_progress_from_thread
    uint64_t total_pending_last;

    // TODO(doc)
    // Atomic so as to avoid races between report_progress_from_thread and the
    // threads modifying total_pending_now
    std::atomic_uint64_t total_pending_now;

    SimsStats();

    SimsStats(SimsStats const& that) : SimsStats() {
      init_from(that);
    }

    SimsStats& operator=(SimsStats const& that) {
      return init_from(that);
    }

    SimsStats(SimsStats&& that) : SimsStats() {
      init_from(that);
    }

    SimsStats& operator=(SimsStats&& that) {
      return init_from(that);
    }

    SimsStats& stats_zero();

    SimsStats& stats_check_point() {
      count_last         = count_now;
      total_pending_last = total_pending_now;
      return *this;
    }

   private:
    SimsStats& init_from(SimsStats const& that);
  };

  //! No doc
  // This class allows us to use the same interface for settings for Sims1,
  // RepOrc, and MinimalRepOrc without duplicating the code.
  template <typename Subclass>
  class SimsSettings {
   protected:
    // These are protected so that Sims1 can reverse them if necessary for
    // left congruences.
    std::vector<word_type>  _exclude;
    std::vector<word_type>  _include;
    Presentation<word_type> _presentation;

   private:
    size_t                                 _idle_thread_restarts;
    std::vector<word_type>::const_iterator _longs_begin;
    size_t                                 _num_threads;
    mutable SimsStats                      _stats;

   public:
    // TODO(doc)
    SimsSettings();

    // TODO(doc)
    // TODO(tests)
    Subclass& init();

    // Copy constructor is explicitly required, the constructor template is not
    // a substitute. If no copy constructor is implemented, then _longs_begin
    // is not properly initialised, and leads to badness.
    SimsSettings(SimsSettings const& that) {
      init(that);
    }

    SimsSettings& operator=(SimsSettings const& that) {
      init(that);
      return *this;
    }

    //! Construct from SimsSettings with different subclass.
    template <typename OtherSubclass>
    SimsSettings(SimsSettings<OtherSubclass> const& that) {
      init(that);
    }

    template <typename OtherSubclass>
    SimsSettings& init(SimsSettings<OtherSubclass> const& that);

    //! Returns the settings object of *this.
    //!
    //! The settings object contains all the settings that are common to
    //! `Sims1`, `RepOrc`, and `MinimalRepOrc`, which are currently:
    //! * \ref presentation
    //! * \ref long_rules
    //! * \ref number_of_threads
    //! * \ref extra
    //!
    //! The return value of this function can be used to
    //! initialise another `Sims1`, `RepOrc`, or
    //! `MinimalRepOrc` with these settings.
    //!
    //! \param (None) this function has no parameters.
    //!
    //! \returns A const reference to `SimsSettings`.
    //!
    //! \exceptions
    //! \noexcept
    // So that we can access the settings from the derived class T.
    [[nodiscard]] SimsSettings const& settings() const noexcept {
      return *this;
    }

    //! Copy the settings from \p that into `this`.
    //!
    //! The settings object contains all the settings that are common to
    //! `Sims1`, `RepOrc`, and `MinimalRepOrc`, which are currently:
    //! * \ref presentation
    //! * \ref long_rules
    //! * \ref number_of_threads
    //! * \ref extra
    //!
    //! The return value of this function can be used to initialise another
    //! `Sims1`, `RepOrc`, or `MinimalRepOrc` with these settings.
    //!
    //! \param that the object to copy the settings from.
    //!
    //! \returns A const reference to `this`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    Subclass& settings_copy_from(SimsSettings const& that) {
      *this = that;
      return static_cast<Subclass&>(*this);
    }

    //! \anchor number_of_threads
    //! Set the number of threads.
    //!
    //! This function sets the number of threads to be used
    //! by `Sims1`.
    //!
    //! The default value is `1`.
    //!
    //! \param val the maximum number of threads to use.
    //!
    //! \returns A reference to \c this.
    //!
    //! \throws LibsemigroupsException if the argument \p
    //! val is 0.
    //!
    //! \warning If \p val exceeds
    //! `std::thread::hardware_concurrency()`, then this is
    //! likely to have a negative impact on the performance
    //! of the algorithms implemented by `Sims1`.
    Subclass& number_of_threads(size_t val);

    //! Returns the current number of threads.
    //!
    //! \param (None) this function has no parameters.
    //!
    //! \returns
    //! A `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] size_t number_of_threads() const noexcept {
      return _num_threads;
    }

    //! Set the short rules.
    //!
    //! These are the rules used at every node in the depth first search
    //! conducted by `Sims1`.
    //!
    //! If the template parameter \p P is not `Presentation<word_type>`, then
    //! the parameter \p p is first converted to a value of type
    //! `Presentation<word_type>` and it is this converted value that is used.
    //!
    //! \tparam P A specific value of the class template `Presentation`, must
    //! be derived from `PresentationBase`. \param p the presentation.
    //!
    //! \returns A reference to \c this.
    //!
    //! \throws LibsemigroupsException if `to_presentation<word_type>(p)`
    //! throws
    //! \throws LibsemigroupsException if `p` is not valid
    //! \throws LibsemigroupsException if the alphabet of `p` is non-empty and
    //! not equal to that of \ref long_rules or \ref extra. \throws
    //! LibsemigroupsException if `p` has 0-generators and 0-relations.
    template <typename PresentationOfSomeKind>
    Subclass& presentation(PresentationOfSomeKind const& p);

    //! \anchor presentation
    //! Returns a const reference to the current short rules.
    //!
    //! This function returns the defining presentation of a `Sims1` instance.
    //! The congruences computed by \ref cbegin and \ref cend are defined over
    //! the semigroup or monoid defined by this presentation.
    //!
    //! Note that it might not be the case that the value returned by this
    //! function and the presentation used to construct the object are the same.
    //! A `Sims1` object requires the generators of the defining presentation
    //! \f$\mathcal{P}\f$ to be \f$\{0, \ldots, n - 1\}\f$ where \f$n\f$ is the
    //! size of the alphabet of \f$\mathcal{P}\f$. Every occurrence of every
    //! generator \c a in the presentation \c p used to construct a `Sims1`
    //! instance is replaced by `p.index(a)`.
    //!
    //! \returns
    //! A const reference to `Presentation<word_type>`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \warning
    //! If \ref split_at or \ref long_rule_length have been
    //! called, then some of the defining relations may
    //! have been moved from \ref presentation to \ref
    //! long_rules.
    [[nodiscard]] Presentation<word_type> const& presentation() const noexcept {
      return _presentation;
    }

    //! Set the long rules.
    //!
    //! These are the rules used after a complete deterministic word graph
    //! compatible with \ref presentation has been found by `Sims1`. If such a
    //! word graph is compatible with the long rules specified by this
    //! function, then this word graph is accepted, and if not it is not
    //! accepted.
    //!
    //! If the template parameter \p P is not `Presentation<word_type>`, then
    //! the parameter \p p is first converted to a value of type
    //! `Presentation<word_type>` and it is this converted value that is used.
    //!
    //! \tparam P A specific value of the class template `Presentation`, must
    //! be derived from `PresentationBase`. \param p the presentation.
    //!
    //! \returns A reference to \c this.
    //!
    //! \throws LibsemigroupsException if `to_presentation<word_type>(p)` throws
    //! \throws LibsemigroupsException if `p` is not valid
    //! \throws LibsemigroupsException if the alphabet of `p` is non-empty and
    //! not equal to that of \ref presentation or \ref extra.
    Subclass& cbegin_long_rules(std::vector<word_type>::const_iterator p);

    // TODO(doc)
    Subclass& cbegin_long_rules(size_t pos);

    // TODO(doc)
    Subclass& clear_long_rules() {
      return cbegin_long_rules(_presentation.rules.cend());
    }

    // TODO(doc)
    [[nodiscard]] size_t number_of_long_rules() const noexcept {
      return std::distance(_longs_begin, _presentation.rules.cend()) / 2;
    }

    //! \anchor long_rules
    //! Returns the current long rules.
    //!
    //! \param (None) this function has no parameters.
    //!
    //! \returns
    //! A const reference to `Presentation<word_type>`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] std::vector<word_type>::const_iterator
    cbegin_long_rules() const noexcept {
      LIBSEMIGROUPS_ASSERT(_presentation.rules.cbegin() <= _longs_begin);
      LIBSEMIGROUPS_ASSERT(_longs_begin <= _presentation.rules.cend());
      return _longs_begin;
    }

    //! \anchor extra
    //! Returns a const reference to the additional defining pairs.
    //!
    //! The congruences computed by a Sims1 instance always contain the
    //! relations of this presentation. In other words, the congruences
    //! computed by this instance are only taken among those that contains the
    //! pairs of elements of the underlying semigroup (defined by the
    //! presentation returned by \ref presentation and \ref long_rules)
    //! represented by the relations of the presentation returned by
    //! `extra()`.
    //!
    //! \param (None) this function has no parameters.
    //!
    //! \returns
    //! A const reference to `Presentation<word_type>`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] std::vector<word_type> const& include() const noexcept {
      return _include;
    }

    //! Set the extra rules.
    //!
    //! The congruences computed by a `Sims1` instance always contain the
    //! relations of this presentation. In other words, the congruences
    //! computed by this instance are only taken among those that contains the
    //! pairs of elements of the underlying semigroup (defined by the
    //! presentation returned by \ref presentation and \ref long_rules)
    //! represented by the relations of the presentation returned by
    //! `extra()`.
    //!
    //! If the template parameter \p P is not `Presentation<word_type>`, then
    //! the parameter \p p is first converted to a value of type
    //! `Presentation<word_type>` and it is this converted value that is used.
    //!
    //! \tparam P A specific value of the class template `Presentation`, must
    //! be derived from `PresentationBase`.
    //!
    //! \param p the presentation.
    //!
    //! \returns A reference to \c this.
    //!
    //! \throws LibsemigroupsException if `to_presentation<word_type>(p)`
    //! throws
    //! \throws LibsemigroupsException if `p` is not valid \throws
    //! LibsemigroupsException if the alphabet of `p` is non-empty and not
    //! equal to that of \ref presentation or \ref long_rules.
    template <typename iterator>
    Subclass& include(iterator first, iterator last);

    // TODO(doc)
    Subclass& include(word_type const& lhs, word_type const& rhs);

    // TODO(doc)
    template <typename Container>
    Subclass& include(Container const& c) {
      include(std::begin(c), std::end(c));
      return static_cast<Subclass&>(*this);
    }

    // TODO(doc)
    Subclass& clear_include() {
      _include.clear();
      return static_cast<Subclass&>(*this);
    }

    // TODO(doc)
    template <typename iterator>
    Subclass& exclude(iterator first, iterator last);

    // TODO(doc)
    template <typename Container>
    Subclass& exclude(Container const& c) {
      exclude(std::begin(c), std::end(c));
      return static_cast<Subclass&>(*this);
    }

    // TODO(doc)
    Subclass& exclude(word_type const& lhs, word_type const& rhs);

    // TODO(doc)
    [[nodiscard]] std::vector<word_type> const& exclude() const noexcept {
      return _exclude;
    }

    // TODO(doc)
    Subclass& clear_exclude() {
      _exclude.clear();
      return static_cast<Subclass&>(*this);
    }

    // TODO(Sims1) ranges version of include/exclude?

    //! Returns a const reference to the current stats object.
    //!
    //! The value returned by this function is a `SimsStats` object which
    //! contains some statistics related to the current `Sims1` instance and
    //! any part of the depth first search already conducted.
    //!
    //! \param (None) this function has no parameters.
    //!
    //! \returns
    //! A const reference to `SimsStats`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] SimsStats& stats() const noexcept {
      return _stats;
    }

    //! \anchor long_rule_length
    //! Define the long rule length.
    //!
    //! This function modifies \ref presentation and \ref long_rules so that
    //! \ref presentation only contains those rules whose length (sum of the
    //! lengths of the two sides of the rules) is less than \p val (if any)
    //! and \ref long_rules only contains those rules of length at least \p
    //! val (if any). The rules contained in the union of \ref presentation
    //! and \ref long_rules is invariant under this function, but the
    //! distribution of the rules between \ref presentation and \ref
    //! long_rules is not.
    //!
    //! The relative orders of the rules within \ref presentation and \ref
    //! long_rules may not be preserved.
    //!
    //! \param val the value of the long rule length.
    //!
    //! \returns
    //! A const reference to `this`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    Subclass& long_rule_length(size_t val);

    // TODO to tpp
    Subclass& idle_thread_restarts(size_t val) {
      if (val == 0) {
        LIBSEMIGROUPS_EXCEPTION(
            "the argument (idle thread restarts) must be non-zero");
      }
      _idle_thread_restarts = val;
      return static_cast<Subclass&>(*this);
    }

    [[nodiscard]] size_t idle_thread_restarts() const noexcept {
      return _idle_thread_restarts;
    }

   protected:
    Subclass const& stats_copy_from(SimsStats const& stts) const {
      _stats = std::move(stts);
      return static_cast<Subclass const&>(*this);
    }

    void reverse(std::vector<word_type>& vec) {
      std::for_each(vec.begin(), vec.end(), [](word_type& w) {
        std::reverse(w.begin(), w.end());
      });
    }

   private:
    template <typename OtherSubclass>
    SimsSettings& init_from(SimsSettings<OtherSubclass> const& that);
  };

  class Sims1;
  class Sims2;

  namespace detail {
    template <typename Sims1or2>
    class SimsBase : public SimsSettings<Sims1or2>, public Reporter {
      static_assert(std::is_same_v<Sims1or2, Sims1>
                    || std::is_same_v<Sims1or2, Sims2>);

     public:
      // We use WordGraph, even though the iterators produced by this class
      // hold FelschGraph's, none of the features of FelschGraph are useful
      // for the output, only for the implementation
      //! The type of the associated WordGraph objects.
      using word_graph_type = WordGraph<uint32_t>;

      //! Type for the nodes in the associated WordGraph objects.
      using node_type = typename word_graph_type::node_type;

      using label_type = typename word_graph_type::label_type;

      //! The size_type of the associated WordGraph objects.
      using size_type = typename word_graph_type::size_type;

      //! Type for letters in the underlying presentation.
      using letter_type = typename word_type::value_type;

     protected:
      ////////////////////////////////////////////////////////////////////////
      // SimsBase nested classes - protected
      ////////////////////////////////////////////////////////////////////////

      struct PendingDefBase {
        PendingDefBase() = default;

        PendingDefBase(node_type   s,
                       letter_type g,
                       node_type   t,
                       size_type   e,
                       size_type   n,
                       bool) noexcept
            : source(s), generator(g), target(t), num_edges(e), num_nodes(n) {}

        node_type   source;
        letter_type generator;
        node_type   target;
        size_type   num_edges;  // Number of edges in the graph when
                                // *this was added to the stack
        size_type num_nodes;    // Number of nodes in the graph
                                // after the definition is made
      };

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
        // TODO(Sims1) ensure that _felsch_graph's settings are
        // properly initialised
        using Definition = std::pair<node_type, label_type>;

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

        // Try to pop from _pending into the argument (reference), returns true
        // if successful and false if not.
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

        // prefix
        iterator const& operator++();

        // postfix
        //! No doc
        iterator operator++(int) {
          iterator copy(*this);
          ++(*this);
          return copy;
        }

        using iterator_base::swap;
      };  // class iterator

     private:
      ////////////////////////////////////////////////////////////////////////
      // SimsBase nested classes - private
      ////////////////////////////////////////////////////////////////////////
      class thread_runner;

      // Note that this class is private, and not really an iterator in the
      // usual sense. It is designed solely to work with thread_runner.
      class thread_iterator : public Sims1or2::iterator_base {
        using PendingDef    = typename Sims1or2::PendingDef;
        using iterator_base = typename Sims1or2::iterator_base;

        friend class thread_runner;

        using iterator_base::partial_copy_for_steal_from;

       public:
        using sims_type = typename iterator_base::sims_type;

        thread_iterator(sims_type const* s, size_type n)
            : iterator_base(s, n) {}

        ~thread_iterator();

        thread_iterator()                                  = delete;
        thread_iterator(thread_iterator const&)            = delete;
        thread_iterator(thread_iterator&&)                 = delete;
        thread_iterator& operator=(thread_iterator const&) = delete;
        thread_iterator& operator=(thread_iterator&&)      = delete;

        using iterator_base::stats;

        // TODO by value really?
        void push(PendingDef pd) {
          iterator_base::_pending.push_back(std::move(pd));
        }

        void steal_from(thread_iterator& that);
        bool try_steal(thread_iterator& that);
      };  // thread_iterator

      class thread_runner {
       private:
        using thread_iterator = typename Sims1or2::thread_iterator;
        using PendingDef      = typename Sims1or2::PendingDef;

        std::atomic_bool                              _done;
        std::vector<std::unique_ptr<thread_iterator>> _theives;
        std::vector<std::thread>                      _threads;
        std::mutex                                    _mtx;
        size_type                                     _num_threads;
        word_graph_type                               _result;
        Sims1or2 const*                               _sims1or2;

        void worker_thread(unsigned                                    my_index,
                           std::function<bool(word_graph_type const&)> hook);

        bool pop_from_local_queue(PendingDef& pd, unsigned my_index) {
          return _theives[my_index]->try_pop(pd);
        }

        bool pop_from_other_thread_queue(PendingDef& pd, unsigned my_index);

       public:
        thread_runner(Sims1or2 const* s, size_type n, size_type num_threads);
        ~thread_runner();

        word_graph_type const& word_graph() const {
          return _result;
        }

        void run(std::function<bool(word_graph_type const&)> hook);
      };  // class thread_runner

     private:
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
      ~SimsBase()                          = default;

      Sims1or2& init();

      using SimsSettings<Sims1or2>::presentation;
      using SimsSettings<Sims1or2>::number_of_threads;
      using SimsSettings<Sims1or2>::stats;
      using SimsSettings<Sims1or2>::include;
      using SimsSettings<Sims1or2>::exclude;
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
    };
  }  // namespace detail

  //! Defined in ``sims.hpp``.
  //!
  //! On this page we describe the functionality relating to the small index
  //! congruence algorithm. The algorithm implemented by this class template
  //! is essentially the low index subgroup algorithm for finitely presented
  //! groups described in Section 5.6 of [Computation with Finitely Presented
  //! Groups](https://doi.org/10.1017/CBO9780511574702) by C. Sims. The low
  //! index subgroups algorithm was adapted for semigroups and monoids by J.
  //! D. Mitchell and M. Tsalakou.
  //!
  //! The purpose of this class is to provide the functions \ref cbegin, \ref
  //! cend, \ref for_each, and \ref find_if which permit iterating through the
  //! one-sided congruences of a semigroup or monoid defined by a presentation
  //! containing (a possibly empty) set of pairs and with at most a given
  //! number of classes. An iterator returned by \ref cbegin points at an
  //! WordGraph instance containing the action of the semigroup or monoid
  //! on the classes of a congruence.
  class Sims1 : public detail::SimsBase<Sims1> {
    using SimsBase = detail::SimsBase<Sims1>;

    friend SimsBase;

    using iterator_base = IteratorBase;
    using PendingDef    = PendingDefBase;

   public:
    //! Type for the nodes in the associated WordGraph
    //! objects.
    using node_type  = uint32_t;
    using label_type = typename WordGraph<node_type>::label_type;

    //! Type for letters in the underlying presentation.
    using letter_type = typename word_type::value_type;

    //! The size_type of the associated WordGraph objects.
    using size_type = typename WordGraph<node_type>::size_type;

    // We use WordGraph, even though the iterators produced by this class
    // hold FelschGraph's, none of the features of FelschGraph are useful
    // for the output, only for the implementation
    //! The type of the associated WordGraph objects.
    using word_graph_type = WordGraph<node_type>;

   private:
    congruence_kind _kind;

   public:
    //! Default constructor
    Sims1() = default;

    // TODO(doc)
    using SimsBase::init;

    //! Construct from \ref congruence_kind.
    //!
    //! \param ck the handedness of the congruences (left or right)
    //!
    //! \throws LibsemigroupsException if \p ck is \ref
    //! congruence_kind::twosided
    //!
    //! \sa \ref cbegin and \ref cend.
    explicit Sims1(congruence_kind ck) : SimsBase(), _kind() {
      kind(ck);
    }

    // TODO(doc)
    template <typename W>
    Sims1(congruence_kind ck, Presentation<W> p) : Sims1(ck) {
      presentation(p);
    }

    //! Default copy constructor.
    Sims1(Sims1 const& other) = default;

    //! Default move constructor.
    Sims1(Sims1&&) = default;

    //! Default copy assignment operator.
    Sims1& operator=(Sims1 const&) = default;

    //! Default move assignment operator.
    Sims1& operator=(Sims1&&) = default;

    // No doc
    ~Sims1() = default;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
    template <typename PresentationOfSomeKind>
    Sims1& presentation(PresentationOfSomeKind const& p) {
      SimsBase::presentation(p);
      if (_kind == congruence_kind::left) {
        presentation::reverse(_presentation);
      }
      return *this;
    }

    // Must accept at least one argument so that we're not calling the 0-arg
    // include() which is const!
    template <typename Arg, typename... Args>
    Sims1& include(Arg arg, Args&&... args) {
      SimsBase::include(arg, std::forward<Args>(args)...);
      if (_kind == congruence_kind::left) {
        SimsBase::reverse(_include);
      }
      return *this;
    }

    // Must accept at least one argument so that we're not calling the 0-arg
    // exclude() which is const!
    template <typename Arg, typename... Args>
    Sims1& exclude(Arg arg, Args&&... args) {
      SimsBase::exclude(arg, std::forward<Args>(args)...);
      if (_kind == congruence_kind::left) {
        SimsBase::reverse(_exclude);
      }
      return *this;
    }

    // This is required because "using SimsBase::include;" tries to
    // drag in all of the include mem fns from SimsBase, which clash with
    // the fns above.
    [[nodiscard]] std::vector<word_type> const& include() const noexcept {
      return SimsBase::include();
    }

    // This is required because "using SimsBase::exclude;" tries to
    // drag in all of the include mem fns from SimsBase, which clash with
    // the fns above.
    [[nodiscard]] std::vector<word_type> const& exclude() const noexcept {
      return SimsBase::exclude();
    }
#endif

    using SimsBase::cbegin_long_rules;
    using SimsBase::number_of_threads;
    using SimsBase::presentation;

    // TODO(doc)
    [[nodiscard]] congruence_kind kind() const noexcept {
      return _kind;
    }

    // TODO(doc)
    Sims1& kind(congruence_kind ck);

    //! Returns the number of one-sided congruences with up to a given
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
    using SimsBase::number_of_congruences;

    //! Apply the function \p pred to every one-sided
    //! congruence with at most \p n classes
    //!
    //! This function is similar to
    //! `std::for_each(begin(n), end(n), pred)` and exists
    //! to:
    //! * provide some feedback on the progress of the
    //!   computation if it runs for more than 1 second.
    //! * allow for the computation of
    //!   `std::for_each(begin(n), end(n), pred)` to be
    //!   performed using \ref number_of_threads in parallel.
    //!
    //! \param n the maximum number of congruence classes.
    //! \param pred the predicate applied to every congruence found.
    //!
    //! \returns (None)
    //!
    //! \throws LibsemigroupsException if \p n is \c 0.
    //! \throws LibsemigroupsException if `presentation()` has 0-generators and
    //! 0-relations (i.e. it has not been initialised).
    // void for_each(size_type                                   n,
    //              std::function<void(word_graph_type const&)> pred) const;
    using SimsBase::for_each;

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
    //! \returns The first congruence whose WordGraph for which \p pred returns
    //! \c true.
    //!
    //! \throws LibsemigroupsException if \p n is \c 0.
    //! \throws LibsemigroupsException if `presentation()` has 0-generators and
    //! 0-relations (i.e. it has not been initialised).
    using SimsBase::find_if;

    //! Returns a forward iterator pointing at the first congruence.
    //!
    //! Returns a forward iterator pointing to the WordGraph representing
    //! the first congruence described by Sims1 object with at most \p n
    //! classes.
    //!
    //! If incremented, the iterator will point to the next such congruence.
    //! The order which the congruences are returned in is implementation
    //! specific. Iterators of the type returned by this function are equal
    //! whenever they point to equal objects. The iterator is exhausted if
    //! and only if it points to an WordGraph with zero nodes.
    //!
    //! The meaning of the WordGraph pointed at by Sims1 iterators depends
    //! on whether the input is a monoid presentation (i.e.
    //! Presentation::contains_empty_word() returns \c true) or a semigroup
    //! presentation. If the input is a monoid presentation for a monoid
    //! \f$M\f$, then the WordGraph pointed to by an iterator of this type
    //! has precisely \p n nodes, and the right action of \f$M\f$ on the
    //! nodes of the word graph is isomorphic to the action of \f$M\f$ on the
    //! classes of a right congruence.
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
    //! \returns
    //! An iterator \c it of type \c iterator pointing to an WordGraph with
    //! at most \p n nodes.
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
    //! \sa
    //! \ref cend
    // TODO(Sims1) it'd be good to remove node 0 to avoid confusion. This
    // seems complicated however, and so isn't done at present.
    using SimsBase::cbegin;

    //! Returns a forward iterator pointing one beyond the last congruence.
    //!
    //! Returns a forward iterator pointing to the empty WordGraph. If
    //! incremented, the returned iterator remains valid and continues to
    //! point at the empty WordGraph.
    //!
    //! \param n the maximum number of classes in a
    //! congruence.
    //!
    //! \returns
    //! An iterator \c it of type \c iterator pointing to an WordGraph with
    //! at most \p 0 nodes.
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
    //! \sa
    //! \ref cbegin
    using SimsBase::cend;
  };

  class Sims2 : public detail::SimsBase<Sims2> {
    using SimsBase = detail::SimsBase<Sims2>;
    friend SimsBase;

   public:
    using node_type       = SimsBase::node_type;
    using label_type      = SimsBase::label_type;
    using letter_type     = SimsBase::letter_type;
    using size_type       = SimsBase::size_type;
    using word_graph_type = SimsBase::word_graph_type;

    Sims2()                        = default;
    Sims2(Sims2 const& other)      = default;
    Sims2(Sims2&&)                 = default;
    Sims2& operator=(Sims2 const&) = default;
    Sims2& operator=(Sims2&&)      = default;
    ~Sims2()                       = default;

    Sims2& init() {
      SimsSettings<Sims2>::init();
      return *this;
    }

    template <typename Word>
    explicit Sims2(Presentation<Word> p) : Sims2() {
      presentation(p);
    }

   private:
    struct PendingDef : public SimsBase::PendingDefBase {
      PendingDef() = default;

      PendingDef(node_type   s,
                 letter_type g,
                 node_type   t,
                 size_type   e,
                 size_type   n,
                 bool        tin) noexcept
          : SimsBase::PendingDefBase(s, g, t, e, n, tin),
            target_is_new_node(tin) {}
      bool target_is_new_node;
    };

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
      using SimsBase::IteratorBase::init;
      using SimsBase::IteratorBase::try_pop;

      // We could use the copy constructor, but there's no point in copying
      // anything except the FelschGraph and so we only copy that.
      void partial_copy_for_steal_from(iterator_base const& that);

      // Try to make the definition represented by PendingDef, returns false if
      // it wasn't possible, and true if it was.
      //! No doc
      [[nodiscard]] bool try_define(PendingDef const&);

      iterator_base(Sims2 const* s, size_type n);

     public:
      iterator_base() = default;

      iterator_base(iterator_base const& that);
      iterator_base(iterator_base&& that);
      iterator_base& operator=(iterator_base const& that);
      iterator_base& operator=(iterator_base&& that);
      ~iterator_base();

      using SimsBase::IteratorBase::operator==;
      using SimsBase::IteratorBase::operator!=;
      using SimsBase::IteratorBase::operator*;
      using SimsBase::IteratorBase::operator->;

      //! No doc
      void swap(iterator_base& that) noexcept;

      using SimsBase::IteratorBase::stats;
    };  // class iterator_base

   public:
    using SimsSettings::cbegin_long_rules;
    using SimsSettings::exclude;
    using SimsSettings::include;
    using SimsSettings::number_of_threads;
    using SimsSettings::presentation;

    using SimsBase::cbegin;
    using SimsBase::cend;
    using SimsBase::find_if;
    using SimsBase::for_each;
    using SimsBase::number_of_congruences;
  };  // Sims2

  //! Defined in ``sims.hpp``.
  //!
  //! This class is a helper for `Sims1` calling the `word_graph` member
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
    // TODO(doc)
    RepOrc& init();

    //! Construct from Sims1 or MinimalRepOrc.
    //!
    //! This constructor creates a new RepOrc instance with
    //! the same SimsSettings as \p s but that is
    //! otherwise uninitialised.
    //!
    //! \tparam S the type of the argument \p s (which is
    //! derived from `SimsSettings<S>`).
    //!
    //! \param s the Sims1 or MinimalRepOrc whose settings
    //! should be used.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    template <typename OtherSubclass>
    explicit RepOrc(SimsSettings<OtherSubclass> const& s) : RepOrc() {
      SimsSettings<RepOrc>::init(s);
    }

    template <typename OtherSubclass>
    RepOrc& init(SimsSettings<OtherSubclass> const& s) {
      SimsSettings<RepOrc>::init(s);
      return *this;
    }

    //! Set the minimum number of nodes.
    //!
    //! This function sets the minimal number of nodes in
    //! the WordGraph that we are seeking.
    //!
    //! \param val the minimum number of nodes
    //!
    //! \returns A reference to `this`.
    //!
    //! \exceptions
    //! \noexcept
    RepOrc& min_nodes(size_t val) noexcept {
      _min = val;
      return *this;
    }

    //! The current minimum number of nodes.
    //!
    //! This function returns the current value for the minimum number of nodes
    //! in the WordGraph that we are seeking.
    //!
    //! \param (None) this function has no parameters.
    //!
    //! \returns A value of type `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] size_t min_nodes() const noexcept {
      return _min;
    }

    //! Set the maximum number of nodes.
    //!
    //! This function sets the maximum number of nodes in the WordGraph that we
    //! are seeking.
    //!
    //! \param val the maximum number of nodes
    //!
    //! \returns A reference to `this`.
    //!
    //! \exceptions
    //! \noexcept
    RepOrc& max_nodes(size_t val) noexcept {
      _max = val;
      return *this;
    }

    //! The current maximum number of nodes.
    //!
    //! This function returns the current value for the maximum number of nodes
    //! in the WordGraph that we are seeking.
    //!
    //! \param (None) this function has no parameters.
    //!
    //! \returns A value of type `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] size_t max_nodes() const noexcept {
      return _max;
    }

    //! Set the target size.
    //!
    //! This function sets the target size, i.e. the desired size of the
    //! transformation semigroup corresponding to the WordGraph returned by the
    //! function \ref word_graph.
    //!
    //! \param val the target size.
    //!
    //! \returns A reference to `this`.
    //!
    //! \exceptions
    //! \noexcept
    RepOrc& target_size(size_t val) noexcept {
      _size = val;
      return *this;
    }

    //! The current target size.
    //!
    //! This function returns the current value for the target size, i.e. the
    //! desired size of the transformation semigroup corresponding to the
    //! WordGraph returned by the function \ref word_graph.
    //!
    //! \param (None) this function has no parameters.
    //!
    //! \returns A value of type `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] size_t target_size() const noexcept {
      return _size;
    }

    //! Get the word_graph.
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
    //! \tparam T the type of the nodes in the returned word graph. \param
    //! (None) this function has no parameters.
    //!
    //! \returns A value of type `WordGraph`.
    //!
    //! \exceptions \no_libsemigroups_except
    [[nodiscard]] Sims1::word_graph_type word_graph() const;

    using SimsSettings<RepOrc>::presentation;
    using SimsSettings<RepOrc>::cbegin_long_rules;
  };

  //! Defined in ``sims.hpp``.
  //!
  //! This class is a helper for `Sims1`, calling the `word_graph` member
  //! function attempts to find a right congruence, represented as an
  //! WordGraph, with the minimum possible number of nodes such that the
  //! action of the semigroup or monoid defined by the presentation consisting
  //! of its \ref presentation and \ref long_rules on the nodes of the
  //! WordGraph corresponds to a semigroup of size \ref target_size.
  //!
  //! If no such WordGraph can be found, then an empty WordGraph is
  //! returned (with `0` nodes and `0` edges).
  class MinimalRepOrc : public SimsSettings<MinimalRepOrc> {
   private:
    size_t _size;

   public:
    using SimsSettings<MinimalRepOrc>::stats;

    //! Default constructor
    MinimalRepOrc() = default;

    MinimalRepOrc& init() {
      _size = 0;
      return *this;
    }

    //! Set the target size.
    //!
    //! This function sets the target size, i.e. the desired size of the
    //! transformation semigroup corresponding to the WordGraph returned by
    //! the function \ref word_graph.
    //!
    //! \param val the target size.
    //!
    //! \returns A reference to `this`.
    //!
    //! \exceptions
    //! \noexcept
    MinimalRepOrc& target_size(size_t val) noexcept {
      _size = val;
      return *this;
    }

    //! The current target size.
    //!
    //! This function returns the current value for the target size, i.e. the
    //! desired size of the transformation semigroup corresponding to the
    //! WordGraph returned by the function \ref word_graph.
    //!
    //! \param (None) this function has no parameters.
    //!
    //! \returns A value of type `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] size_t target_size() const noexcept {
      return _size;
    }

    //! Get the word graph.
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
    //! \param
    //! (None) this function has no parameters.
    //!
    //! \returns A value of type `WordGraph<uint32_t>`.
    //!
    //! \exceptions \no_libsemigroups_except
    [[nodiscard]] Sims1::word_graph_type word_graph() const;
  };

}  // namespace libsemigroups

#include "sims.tpp"

#endif  // LIBSEMIGROUPS_SIMS_HPP_
