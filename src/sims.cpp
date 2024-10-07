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

// This file contains the implementations of the classes for performing the
// "low-index congruence" algorithm (and using it) for semigroups and monoid.

// Notes:
// 1. In 2022, when first writing this file, JDM tried templating the word_type
//    used by the presentations in Sims1 (so that we could use StaticVector1
//    for example, using smaller integer types for letters, and the stack to
//    hold the words rather than the heap), but this didn't seem to give any
//    performance improvement, and so I backed out the changes.

#include "libsemigroups/sims.hpp"

#include <algorithm>   // for max, find_if, fill
#include <chrono>      // for duration, durat...
#include <functional>  // for ref
#include <memory>      // for unique_ptr, mak...
#include <string>      // for operator+, basi...
#include <thread>      // for thread, yield
#include <utility>     // for swap

#include "libsemigroups/constants.hpp"        // for operator!=, ope...
#include "libsemigroups/debug.hpp"            // for LIBSEMIGROUPS_A...
#include "libsemigroups/exception.hpp"        // for LIBSEMIGROUPS_E...
#include "libsemigroups/froidure-pin.hpp"     // for FroidurePin
#include "libsemigroups/presentation.hpp"     // for Presentation
#include "libsemigroups/to-froidure-pin.hpp"  // for to_word_graph
#include "libsemigroups/transf.hpp"           // for Transf, validate
#include "libsemigroups/types.hpp"            // for congruence_kind
#include "libsemigroups/word-graph.hpp"       // for follow_path_no_...

#include "libsemigroups/detail/felsch-graph.hpp"  // for FelschGraph
#include "libsemigroups/detail/report.hpp"        // for report_default
#include "libsemigroups/detail/stl.hpp"           // for JoinThreads
#include "libsemigroups/detail/string.hpp"        // for group_digits
#include "libsemigroups/detail/timer.hpp"         // for Timer

namespace libsemigroups {

  ////////////////////////////////////////////////////////////////////////
  // SimsStats
  ////////////////////////////////////////////////////////////////////////

  SimsStats::SimsStats()
      : count_last(),
        count_now(),
        max_pending(),
        total_pending_last(),
        total_pending_now() {
    stats_zero();
  }

  SimsStats& SimsStats::stats_zero() {
    count_last         = 0;
    count_now          = 0;
    max_pending        = 0;
    total_pending_last = 0;
    total_pending_now  = 0;
    return *this;
  }

  SimsStats& SimsStats::init(SimsStats const& that) {
    count_last         = that.count_last.load();
    count_now          = that.count_now.load();
    max_pending        = that.max_pending.load();
    total_pending_last = that.total_pending_last.load();
    total_pending_now  = that.total_pending_now.load();
    return *this;
  }

  ////////////////////////////////////////////////////////////////////////
  // SimsSettings
  ////////////////////////////////////////////////////////////////////////

  template <typename Subclass>
  SimsSettings<Subclass>::SimsSettings()
      : _exclude(),
        _exclude_pruner_index(UNDEFINED),
        _idle_thread_restarts(64),
        _include(),
        _longs_begin(),
        _num_threads(1),
        _presentation(),
        _pruners(),
        _stats() {
    _longs_begin = _presentation.rules.cend();
  }

  template <typename Subclass>
  Subclass& SimsSettings<Subclass>::init() {
    _exclude.clear();
    _exclude_pruner_index = UNDEFINED;
    _idle_thread_restarts = 64;
    _include.clear();
    _num_threads = 1;
    _presentation.init();
    _pruners.clear();

    // Need to set after presentation is initialized to avoid out of bound
    // access.
    _longs_begin = _presentation.rules.cend();

    _stats.init();

    return static_cast<Subclass&>(*this);
  }

  template <typename Subclass>
  Subclass& SimsSettings<Subclass>::cbegin_long_rules(
      std::vector<word_type>::const_iterator it) {
    auto const& rules = presentation().rules;
    if (!(rules.cbegin() <= it && it <= rules.cend())) {
      LIBSEMIGROUPS_EXCEPTION(
          "expected an iterator pointing into presentation().rules()");
    } else if (std::distance(it, rules.cend()) % 2 != 0) {
      LIBSEMIGROUPS_EXCEPTION(
          "expected an iterator pointing at the left hand side of a rule "
          "(an even distance from the end of the rules), found distance {}",
          std::distance(it, rules.cend()));
    }
    _longs_begin = it;
    return static_cast<Subclass&>(*this);
  }

  template <typename Subclass>
  Subclass& SimsSettings<Subclass>::cbegin_long_rules(size_t pos) {
    return cbegin_long_rules(presentation().rules.cbegin() + pos);
  }

  template <typename Subclass>
  Subclass& SimsSettings<Subclass>::number_of_threads(size_t val) {
    if (val == 0) {
      LIBSEMIGROUPS_EXCEPTION(
          "the argument (number of threads) must be non-zero");
    }

    if ((std::thread::hardware_concurrency() > 0)
        && (val > std::thread::hardware_concurrency())) {
      val = std::thread::hardware_concurrency();
    }
    _num_threads = val;
    return static_cast<Subclass&>(*this);
  }

  template <typename Subclass>
  Subclass& SimsSettings<Subclass>::long_rule_length(size_t val) {
    presentation::sort_rules(_presentation);
    auto& rules = _presentation.rules;
    auto  it    = rules.cbegin();

    for (; it < rules.cend(); it += 2) {
      if (it->size() + (it + 1)->size() >= val) {
        break;
      }
    }

    _longs_begin = it;
    return static_cast<Subclass&>(*this);
  }

  template <typename Subclass>
  Subclass& SimsSettings<Subclass>::idle_thread_restarts(size_t val) {
    if (val == 0) {
      LIBSEMIGROUPS_EXCEPTION(
          "the argument (idle thread restarts) must be non-zero");
    }
    _idle_thread_restarts = val;
    return static_cast<Subclass&>(*this);
  }

  template <typename Subclass>
  size_t SimsSettings<Subclass>::add_exclude_pruner() {
    if (_exclude_pruner_index != UNDEFINED) {
      return _exclude_pruner_index;
    }
    auto pruner = [this](auto const& wg) {
      auto      first = _exclude.cbegin();
      auto      last  = _exclude.cend();
      node_type root  = 0;

      for (auto it = first; it != last; it += 2) {
        auto l = word_graph::follow_path_no_checks(wg, root, *it);
        if (l != UNDEFINED) {
          auto r = word_graph::follow_path_no_checks(wg, root, *(it + 1));
          if (l == r) {
            return false;
          }
        }
      }
      return true;
    };
    // NOTE:: setting _exclude_pruner_index here depends on the implementation
    // of add_pruner, if that changes bad things may happen.
    _exclude_pruner_index = _pruners.size();
    add_pruner(pruner);
    return _exclude_pruner_index;
  }

  template <typename Subclass>
  Subclass& SimsSettings<Subclass>::include_exclude(
      word_type const&        lhs,
      word_type const&        rhs,
      std::vector<word_type>& include_or_exclude) {
    presentation().validate_word(lhs.cbegin(), lhs.cend());
    presentation().validate_word(rhs.cbegin(), rhs.cend());
    include_or_exclude.push_back(lhs);
    include_or_exclude.push_back(rhs);
    return static_cast<Subclass&>(*this);
  }

  template class SimsSettings<Sims1>;
  template class SimsSettings<Sims2>;
  template class SimsSettings<RepOrc>;
  template class SimsSettings<MinimalRepOrc>;

  ////////////////////////////////////////////////////////////////////////
  // Sims1and2::PendingDef
  ////////////////////////////////////////////////////////////////////////

  struct Sims1::PendingDef {
    PendingDef() = default;

    PendingDef(node_type   s,
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

  struct Sims2::PendingDef : public Sims1::PendingDef {
    PendingDef() = default;

    PendingDef(node_type   s,
               letter_type g,
               node_type   t,
               size_type   e,
               size_type   n,
               bool        tin) noexcept
        : Sims1::PendingDef(s, g, t, e, n, tin), target_is_new_node(tin) {}
    bool target_is_new_node;
  };

  namespace detail {
    ////////////////////////////////////////////////////////////////////////
    // SimsBase - public
    ////////////////////////////////////////////////////////////////////////

    template <typename Sims1or2>
    SimsBase<Sims1or2>::~SimsBase() = default;

    template <typename Sims1or2>
    Sims1or2& SimsBase<Sims1or2>::init() {
      if constexpr (std::is_same_v<Sims1or2, Sims1>) {
        report_prefix("Sims1");
      } else {
        report_prefix("Sims2");
      }
      SimsSettings<Sims1or2>::init();
      return static_cast<Sims1or2&>(*this);
    }

    template <typename Sims1or2>
    SimsBase<Sims1or2>::SimsBase() {
      init();
    }

    template <typename Sims1or2>
    void SimsBase<Sims1or2>::for_each(
        size_type                                   n,
        std::function<void(word_graph_type const&)> pred) const {
      throw_if_not_ready(n);

      report_at_start(n);
      if (number_of_threads() == 1) {
        if (!reporting_enabled()) {
          // Don't care about stats in this case
          std::for_each(cbegin(n), cend(n), pred);
        } else {
          stats().stats_zero();
          detail::Ticker t([this]() { report_progress_from_thread(); });
          auto           it   = cbegin(n);
          auto const     last = cend(n);
          for (; it != last; ++it) {
            pred(*it);
          }
          report_final();
        }
      } else {
        thread_runner den(
            static_cast<Sims1or2 const*>(this), n, number_of_threads());
        auto pred_wrapper = [&pred](word_graph_type const& ad) {
          pred(ad);
          return false;
        };
        if (!reporting_enabled()) {
          den.run(pred_wrapper);
        } else {
          stats().stats_zero();
          detail::Ticker t([this]() { report_progress_from_thread(); });
          den.run(pred_wrapper);
          report_final();
        }
      }
    }

    template <typename Sims1or2>
    typename SimsBase<Sims1or2>::word_graph_type SimsBase<Sims1or2>::find_if(
        size_type                                   n,
        std::function<bool(word_graph_type const&)> pred) const {
      throw_if_not_ready(n);
      report_at_start(n);
      if (number_of_threads() == 1) {
        if (!reporting_enabled()) {
          return *std::find_if(cbegin(n), cend(n), pred);
        } else {
          stats().stats_zero();
          detail::Ticker t([this]() { report_progress_from_thread(); });

          auto       it   = cbegin(n);
          auto const last = cend(n);

          for (; it != last; ++it) {
            if (pred(*it)) {
              report_final();
              return *it;
            }
          }
          report_final();
          return *last;  // the empty digraph
        }
      } else {
        thread_runner den(
            static_cast<Sims1or2 const*>(this), n, number_of_threads());
        if (!reporting_enabled()) {
          den.run(pred);
          return den.word_graph();
        } else {
          stats().stats_zero();
          detail::Ticker t([this]() { report_progress_from_thread(); });
          den.run(pred);
          report_final();
          return den.word_graph();
        }
      }
    }

    template <typename Sims1or2>
    uint64_t SimsBase<Sims1or2>::number_of_congruences(size_type n) const {
      if (number_of_threads() == 1) {
        uint64_t result = 0;
        for_each(n, [&result](word_graph_type const&) { ++result; });
        return result;
      } else {
        std::atomic_uint64_t result(0);
        for_each(n, [&result](word_graph_type const&) { ++result; });
        return result;
      }
    }

    ////////////////////////////////////////////////////////////////////////
    // SimsBase - private
    ////////////////////////////////////////////////////////////////////////

    template <typename Sims1or2>
    void SimsBase<Sims1or2>::report_at_start(size_t num_classes) const {
      std::string num_threads = "0";
      if (number_of_threads() != 1) {
        num_threads = fmt::format("{} / {}",
                                  number_of_threads(),
                                  std::thread::hardware_concurrency());
      }
      auto shortest_short = presentation::shortest_rule_length(presentation());
      auto longest_short  = presentation::longest_rule_length(presentation());

      std::string pairs;
      if (!include().empty() && !exclude().empty()) {
        pairs = fmt::format(", including {} + excluding {} pairs",
                            include().size() / 2,
                            exclude().size() / 2);
      } else if (!include().empty()) {
        pairs = fmt::format(", including {} pairs", include().size() / 2);
      } else if (!exclude().empty()) {
        pairs = fmt::format(", excluding {} pairs", exclude().size() / 2);
      }

      report_no_prefix("{:+<80}\n", "");
      report_default("{}: STARTING with {} additional threads . . . \n",
                     report_prefix(),
                     num_threads);
      report_no_prefix("{:+<80}\n", "");
      report_default("{}: \u2264 {} classes{} for \u27E8A|R\u27E9 with:\n",
                     report_prefix(),
                     num_classes,
                     pairs);
      report_default("{}: |A| = {}, |R| = {}, "
                     "|u| + |v| \u2208 [{}, {}], \u2211(|u| + |v|) = {}\n",
                     report_prefix(),
                     presentation().alphabet().size(),
                     presentation().rules.size() / 2,
                     shortest_short,
                     longest_short,
                     presentation::length(presentation()));

      if (cbegin_long_rules() != presentation().rules.cend()) {
        auto first = presentation().rules.cbegin(), last = cbegin_long_rules();

        report_default("{}: {} \"short\" relations with: ",
                       report_prefix(),
                       std::distance(first, last) / 2);
        report_no_prefix(
            "|u| + |v| \u2208 [{}, {}] and \u2211(|u| + |v|) = {}\n",
            presentation::shortest_rule_length(first, last),
            presentation::longest_rule_length(first, last),
            presentation::length(first, last));

        first = cbegin_long_rules(), last = presentation().rules.cend();
        report_default("{}: {} \"long\" relations with: ",
                       report_prefix(),
                       std::distance(first, last) / 2);
        report_no_prefix(
            "|u| + |v| \u2208 [{}, {}] and \u2211(|u| + |v|) = {}\n",
            presentation::shortest_rule_length(first, last),
            presentation::longest_rule_length(first, last),
            presentation::length(first, last));
      }
      Reporter::reset_start_time();
    }

    template <typename Sims1or2>
    void SimsBase<Sims1or2>::report_progress_from_thread() const {
      using namespace detail;       // NOLINT(build/namespaces)
      using namespace std::chrono;  // NOLINT(build/namespaces)

      auto time_total_ns = delta(start_time());
      auto time_diff     = delta(last_report());

      // Stats
      auto count_now          = stats().count_now.load();
      auto count_diff         = count_now - stats().count_last;
      auto total_pending_now  = stats().total_pending_now.load();
      auto total_pending_diff = total_pending_now - stats().total_pending_last;

      constexpr uint64_t billion = 1'000'000'000;
      uint64_t congs_per_sec = (billion * count_now) / time_total_ns.count();
      uint64_t nodes_per_sec
          = (billion * total_pending_now) / time_total_ns.count();

      nanoseconds time_per_cong_last_sec(0);
      if (count_diff != 0) {
        time_per_cong_last_sec = time_diff / count_diff;
      }

      nanoseconds time_per_node_last_sec(0);
      if (total_pending_diff != 0) {
        time_per_node_last_sec = time_diff / total_pending_diff;
      }

      nanoseconds time_per_cong(0);
      if (count_now != 0) {
        time_per_cong = time_total_ns / count_now;
      }

      nanoseconds time_per_node(0);
      if (total_pending_now != 0) {
        time_per_node = time_total_ns / total_pending_now;
      }

      ReportCell<3> rc;
      rc.min_width(0, 4).min_width(1, 7).min_width(2, 11);
      rc("{}: total        {} (cong.)   | {} (nodes) \n",
         report_prefix(),
         group_digits(count_now),
         group_digits(total_pending_now));

      rc("{}: diff         {} (cong.)   | {} (nodes)\n",
         report_prefix(),
         signed_group_digits(count_diff),
         signed_group_digits(total_pending_diff));

      rc("{}: mean         {} (cong./s) | {} (node/s)\n",
         report_prefix(),
         group_digits(congs_per_sec),
         group_digits(nodes_per_sec));

      rc("{}: time last s. {} (/cong.)  | {} (/node)\n",
         report_prefix(),
         string_time(time_per_cong_last_sec),
         string_time(time_per_node_last_sec));

      rc("{}: mean time    {} (/cong.)  | {} (/node)\n",
         report_prefix(),
         string_time(time_per_cong),
         string_time(time_per_node));

      rc("{}: time         {} (total)   |\n",
         report_prefix(),
         string_time(time_total_ns));

      reset_last_report();
      stats().stats_check_point();
    }

    template <typename Sims1or2>
    void SimsBase<Sims1or2>::report_final() const {
      report_progress_from_thread();
      report_no_prefix("{:+<80}\n", "");
      report_default("{}: FINISHED!\n", report_prefix());
      report_no_prefix("{:+<80}\n", "");
    }

    template <typename Sims1or2>
    void SimsBase<Sims1or2>::throw_if_not_ready(size_type n) const {
      if (n == 0) {
        LIBSEMIGROUPS_EXCEPTION(
            "the argument (number of classes) must be non-zero");
      } else if (presentation().rules.empty()
                 && presentation().alphabet().empty()) {
        LIBSEMIGROUPS_EXCEPTION("the presentation() must be defined before "
                                "calling this function");
      }
    }

    ////////////////////////////////////////////////////////////////////////
    // SimsBase::IteratorBase - protected
    ////////////////////////////////////////////////////////////////////////

    // The following function is separated from the constructor
    // so that it isn't called in the constructor of every thread_iterator
    template <typename Sims1or2>
    void SimsBase<Sims1or2>::IteratorBase::init(size_type n) {
      if (n != 0) {
        if (n > 1 || _min_target_node == 1) {
          _pending.emplace_back(0, 0, 1, 0, 2, true);
        }
        if (_min_target_node == 0) {
          _pending.emplace_back(0, 0, 0, 0, 1, false);
        }
      }
    }

    // Try to pop from _pending into the argument (reference), returns true
    // if successful and false if not.
    template <typename Sims1or2>
    bool SimsBase<Sims1or2>::IteratorBase::try_pop(PendingDef& pd) {
      std::lock_guard<std::mutex> lock(_mtx);
      if (_pending.empty()) {
        return false;
      }
      pd = std::move(_pending.back());
      _pending.pop_back();
      return true;
    }

    // Try to make the definition represented by PendingDef, returns false
    // if it wasn't possible, and true if it was.
    template <typename Sims1or2>
    bool
    SimsBase<Sims1or2>::IteratorBase::try_define(PendingDef const& current) {
      // Try to make the definition represented by PendingDef, returns false
      // if it wasn't possible, and true if it was.
      //! No doc
      LIBSEMIGROUPS_ASSERT(current.target < current.num_nodes);
      LIBSEMIGROUPS_ASSERT(current.num_nodes <= _max_num_classes);
      std::lock_guard<std::mutex> lock(_mtx);
      // Backtrack if necessary
      _felsch_graph.reduce_number_of_edges_to(current.num_edges);

      // It might be that current.target is a new node, in which case
      // _felsch_graph.number_of_active_nodes() includes this new node even
      // before the edge current.source -> current.target is defined.
      _felsch_graph.number_of_active_nodes(current.num_nodes);

      LIBSEMIGROUPS_ASSERT(
          _felsch_graph.target_no_checks(current.source, current.generator)
          == UNDEFINED);

      // Don't call number_of_edges because this calls the function in
      // WordGraph
      size_type start = _felsch_graph.definitions().size();

      _felsch_graph.target_no_checks(
          current.source, current.generator, current.target);

      auto first = _sims1or2->include().cbegin();
      auto last  = _sims1or2->include().cend();
      while (start != _felsch_graph.definitions().size()) {
        if (!_felsch_graph.process_definitions(start)) {
          return false;
        }
        start = _felsch_graph.definitions().size();
        if (!felsch_graph::make_compatible<RegisterDefs>(
                _felsch_graph, 0, 1, first, last)) {
          return false;
        }
      }

      return true;
    }

    template <typename Sims1or2>
    bool SimsBase<Sims1or2>::IteratorBase::install_descendents(
        PendingDef const& current) {
      for (auto const& pruner : _sims1or2->pruners()) {
        if (!pruner(_felsch_graph)) {
          return false;
        }
      }

      letter_type     a        = current.generator + 1;
      size_type const M        = _felsch_graph.number_of_active_nodes();
      size_type const N        = _felsch_graph.number_of_edges();
      size_type const num_gens = _felsch_graph.out_degree();
      auto&           stats    = _sims1or2->stats();

      for (node_type next = current.source; next < M; ++next) {
        for (; a < num_gens; ++a) {
          if (_felsch_graph.target_no_checks(next, a) == UNDEFINED) {
            std::lock_guard<std::mutex> lock(_mtx);
            if (M < _max_num_classes) {
              _pending.emplace_back(next, a, M, N, M + 1, true);
            }
            for (node_type b = M; b-- > _min_target_node;) {
              _pending.emplace_back(next, a, b, N, M, false);
            }
            stats.total_pending_now
                += M - _min_target_node + (M < _max_num_classes);

            // Mutex must be locked here so that we can call _pending.size()
            stats.max_pending = std::max(static_cast<uint64_t>(_pending.size()),
                                         stats.max_pending.load());
            return false;
          }
        }
        a = 0;
      }
      // No undefined edges, word graph is complete
      LIBSEMIGROUPS_ASSERT(N == M * num_gens);

      auto first = _sims1or2->cbegin_long_rules();
      auto last  = _sims1or2->presentation().rules.cend();

      bool result = word_graph::is_compatible_no_checks(
          _felsch_graph,
          _felsch_graph.cbegin_nodes(),
          _felsch_graph.cbegin_nodes() + M,
          first,
          last);
      if (result) {
        // stats.count_now is atomic so this is ok
        ++stats.count_now;
      }
      return result;
    }

    template <typename Sims1or2>
    SimsBase<Sims1or2>::IteratorBase::IteratorBase(Sims1or2 const* s,
                                                   size_type       n)
        :  // private
          _max_num_classes(s->presentation().contains_empty_word() ? n : n + 1),
          _min_target_node(s->presentation().contains_empty_word() ? 0 : 1),
          // protected
          _felsch_graph(),
          _mtx(),
          _pending(),
          _sims1or2(s) {
      Presentation<word_type> p = s->presentation();
      size_t m = std::distance(s->presentation().rules.cbegin(),
                               s->cbegin_long_rules());
      p.rules.erase(p.rules.begin() + m, p.rules.end());
      _felsch_graph.init(std::move(p));
      // n == 0 only when the iterator is cend
      _felsch_graph.number_of_active_nodes(n == 0 ? 0 : 1);
      // = 0 indicates iterator is done
      _felsch_graph.add_nodes(n);
    }

    ////////////////////////////////////////////////////////////////////////
    // SimsBase::IteratorBase - public
    ////////////////////////////////////////////////////////////////////////

    template <typename Sims1or2>
    SimsBase<Sims1or2>::IteratorBase::IteratorBase() = default;

    // Intentionally don't copy the mutex, it doesn't compile, wouldn't make
    // sense if the mutex was used here.
    template <typename Sims1or2>
    SimsBase<Sims1or2>::IteratorBase::IteratorBase(IteratorBase const& that)
        :  // private
          _max_num_classes(that._max_num_classes),
          _min_target_node(that._min_target_node),
          // protected
          _felsch_graph(that._felsch_graph),
          _mtx(),
          _pending(that._pending),
          _sims1or2(that._sims1or2) {}

    // Intentionally don't copy the mutex, it doesn't compile, wouldn't make
    // sense if the mutex was used here.
    template <typename Sims1or2>
    SimsBase<Sims1or2>::IteratorBase::IteratorBase(IteratorBase&& that)
        :  // private
          _max_num_classes(std::move(that._max_num_classes)),
          _min_target_node(std::move(that._min_target_node)),
          // protected
          _felsch_graph(std::move(that._felsch_graph)),
          _mtx(),
          _pending(std::move(that._pending)),
          _sims1or2(that._sims1or2) {}

    // Intentionally don't copy the mutex, it doesn't compile, wouldn't make
    // sense if the mutex was used here.
    template <typename Sims1or2>
    typename SimsBase<Sims1or2>::IteratorBase&
    SimsBase<Sims1or2>::IteratorBase::operator=(IteratorBase const& that) {
      // private
      _max_num_classes = that._max_num_classes;
      _min_target_node = that._min_target_node;
      // protected
      _felsch_graph = that._felsch_graph;
      // keep our own _mtx
      _pending  = that._pending;
      _sims1or2 = that._sims1or2;

      return *this;
    }

    // Intentionally don't copy the mutex, it doesn't compile, wouldn't make
    // sense if the mutex was used here.
    template <typename Sims1or2>
    typename SimsBase<Sims1or2>::IteratorBase&
    SimsBase<Sims1or2>::IteratorBase::operator=(IteratorBase&& that) {
      // private
      _max_num_classes = std::move(that._max_num_classes);
      _min_target_node = std::move(that._min_target_node);

      // protected
      _felsch_graph = std::move(that._felsch_graph);
      _pending      = std::move(that._pending);
      // keep our own _mtx
      _sims1or2 = that._sims1or2;
      return *this;
    }

    template <typename Sims1or2>
    SimsBase<Sims1or2>::IteratorBase::~IteratorBase() = default;

    template <typename Sims1or2>
    void SimsBase<Sims1or2>::IteratorBase::swap(IteratorBase& that) noexcept {
      // private
      std::swap(_max_num_classes, that._max_num_classes);
      std::swap(_min_target_node, that._min_target_node);
      // protected
      std::swap(_felsch_graph, that._felsch_graph);
      std::swap(_pending, that._pending);
      std::swap(_sims1or2, that._sims1or2);
    }

    ////////////////////////////////////////////////////////////////////////
    // SimsBase::iterator - private
    ////////////////////////////////////////////////////////////////////////

    template <typename Sims1or2>
    SimsBase<Sims1or2>::iterator::iterator(Sims1or2 const* s, size_type n)
        : iterator_base(s, n) {
      if (this->_felsch_graph.number_of_active_nodes() == 0) {
        return;
      }
      init(n);
      ++(*this);
      // The increment above is required so that when dereferencing any
      // instance of this type we obtain a valid word graph (o/w the value
      // pointed to here is empty).
    }

    ////////////////////////////////////////////////////////////////////////
    // SimsBase::iterator - public
    ////////////////////////////////////////////////////////////////////////

    template <typename Sims1or2>
    SimsBase<Sims1or2>::iterator::~iterator() = default;

    template <typename Sims1or2>
    typename SimsBase<Sims1or2>::iterator const&
    SimsBase<Sims1or2>::iterator::operator++() {
      PendingDef current;
      while (try_pop(current)) {
        if (try_define(current) && install_descendents(current)) {
          return *this;
        }
      }
      this->_felsch_graph.number_of_active_nodes(0);
      // indicates that the iterator is done
      this->_felsch_graph.induced_subgraph_no_checks(0, 0);
      return *this;
    }

    ////////////////////////////////////////////////////////////////////////
    // SimsBase::thread_iterator - public
    ////////////////////////////////////////////////////////////////////////

    // Note that this class is private, and not really an iterator in the
    // usual sense. It is designed solely to work with thread_runner.
    template <typename Sims1or2>
    class SimsBase<Sims1or2>::thread_iterator : public Sims1or2::iterator_base {
      using PendingDef    = typename Sims1or2::PendingDef;
      using iterator_base = typename Sims1or2::iterator_base;

      friend class thread_runner;

      using iterator_base::partial_copy_for_steal_from;

     public:
      using sims_type = typename iterator_base::sims_type;

      thread_iterator(sims_type const* s, size_type n) : iterator_base(s, n) {}

      ~thread_iterator();

      thread_iterator()                                  = delete;
      thread_iterator(thread_iterator const&)            = delete;
      thread_iterator(thread_iterator&&)                 = delete;
      thread_iterator& operator=(thread_iterator const&) = delete;
      thread_iterator& operator=(thread_iterator&&)      = delete;

      using iterator_base::stats;

      void push(PendingDef&& pd) {
        iterator_base::_pending.push_back(std::move(pd));
      }

      void steal_from(thread_iterator& that);
      bool try_steal(thread_iterator& that);
    };  // thread_iterator

    template <typename Sims1or2>
    class SimsBase<Sims1or2>::thread_runner {
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

    template <typename Sims1or2>
    SimsBase<Sims1or2>::thread_iterator::~thread_iterator() = default;

    template <typename Sims1or2>
    void
    SimsBase<Sims1or2>::thread_iterator::steal_from(thread_iterator& that) {
      // WARNING <that> must be locked before calling this function
      std::lock_guard<std::mutex> lock(iterator_base::_mtx);
      LIBSEMIGROUPS_ASSERT(iterator_base::_pending.empty());
      size_t const n = that._pending.size();
      if (n == 1) {
        return;
      }
      partial_copy_for_steal_from(that);

      // Unzip that._pending into _pending and that._pending, this seems to
      // give better performance in the search than splitting that._pending
      // into [begin, begin + size / 2) and [begin + size / 2, end)
      size_t i = 0;
      for (; i < n - 2; i += 2) {
        iterator_base::_pending.push_back(std::move(that._pending[i]));
        that._pending[i / 2] = std::move(that._pending[i + 1]);
      }
      iterator_base::_pending.push_back(std::move(that._pending[i]));
      if (i == n - 2) {
        that._pending[i / 2] = std::move(that._pending[i + 1]);
      }

      that._pending.erase(that._pending.cbegin() + that._pending.size() / 2,
                          that._pending.cend());
    }

    // TODO(0) expand to longer comment about thread sanitizer warnings and how
    // we can gleefully ignore them
    template <typename Sims1or2>
    bool SimsBase<Sims1or2>::thread_iterator::try_steal(thread_iterator& that) {
      std::lock_guard<std::mutex> lock(iterator_base::_mtx);
      if (iterator_base::_pending.empty()) {
        return false;
      }
      // Copy the FelschGraph and half pending from *this into that
      that.steal_from(
          *this);  // Must call steal_from on that, so that that is locked
      return true;
    }

    ////////////////////////////////////////////////////////////////////////
    // SimsBase::thread_runner - private
    ////////////////////////////////////////////////////////////////////////

    template <typename Sims1or2>
    void SimsBase<Sims1or2>::thread_runner::worker_thread(
        unsigned                                    my_index,
        std::function<bool(word_graph_type const&)> hook) {
      PendingDef pd;
      auto const restarts = _sims1or2->idle_thread_restarts();
      for (size_t i = 0; i < restarts; ++i) {
        while ((pop_from_local_queue(pd, my_index)
                || pop_from_other_thread_queue(pd, my_index))
               && !_done) {
          if (_theives[my_index]->try_define(pd)
              && _theives[my_index]->install_descendents(pd)) {
            if (hook(**_theives[my_index])) {
              // hook returns true to indicate that we should stop early
              std::lock_guard<std::mutex> lock(_mtx);
              if (!_done) {
                _done   = true;
                _result = **_theives[my_index];
              }
              return;
            }
          }
        }
        std::this_thread::yield();
        // It's possible to reach here before all of the work is done,
        // because by coincidence there's nothing in the local queue and
        // nothing in any other queue either, this sometimes leads to
        // threads shutting down earlier than desirable. On the other hand,
        // maybe this is a desirable.
      }
    }

    template <typename Sims1or2>
    bool SimsBase<Sims1or2>::thread_runner::pop_from_other_thread_queue(
        PendingDef& pd,
        unsigned    my_index) {
      for (size_t i = 0; i < _theives.size() - 1; ++i) {
        unsigned const index = (my_index + i + 1) % _theives.size();
        // Could always do something different here, like find
        // the largest queue and steal from that? I tried this and it didn't
        // seem to be faster.
        if (_theives[index]->try_steal(*_theives[my_index])) {
          return pop_from_local_queue(pd, my_index);
        }
      }
      return false;
    }

    ////////////////////////////////////////////////////////////////////////
    // SimsBase::thread_runner - public
    ////////////////////////////////////////////////////////////////////////

    template <typename Sims1or2>
    SimsBase<Sims1or2>::thread_runner::thread_runner(Sims1or2 const* s,
                                                     size_type       n,
                                                     size_type num_threads)
        : _done(false),
          _theives(),
          _threads(),
          _mtx(),
          _num_threads(num_threads),
          _result(),
          _sims1or2(s) {
      for (size_t i = 0; i < _num_threads; ++i) {
        _theives.push_back(std::make_unique<thread_iterator>(s, n));
      }
      _theives.front()->init(n);
    }

    template <typename Sims1or2>
    SimsBase<Sims1or2>::thread_runner::~thread_runner() = default;

    template <typename Sims1or2>
    void SimsBase<Sims1or2>::thread_runner::run(
        std::function<bool(word_graph_type const&)> hook) {
      try {
        detail::JoinThreads joiner(_threads);
        for (size_t i = 0; i < _num_threads; ++i) {
          _threads.push_back(std::thread(
              &thread_runner::worker_thread, this, i, std::ref(hook)));
        }
      } catch (...) {
        _done = true;
        throw;
      }
      // TODO(2) Implement a ThreadIdGuard so that thread ids are automatically
      // reset after they exit scope.
      detail::reset_thread_ids();
    }

    ////////////////////////////////////////////////////////////////////////
    // Instantiate the actual classes
    ////////////////////////////////////////////////////////////////////////
    template class SimsBase<Sims1>;
    template class SimsBase<Sims2>;
  }  // namespace detail

  ////////////////////////////////////////////////////////////////////////
  // Sims2
  ////////////////////////////////////////////////////////////////////////

  class Sims2::iterator_base::RuleContainer {
   private:
    std::vector<word_type> _words;
    // _used_slots[i] is the length of _words when we have i edges
    std::vector<size_t> _used_slots;

   public:
    RuleContainer()                                = default;
    RuleContainer(RuleContainer const&)            = default;
    RuleContainer(RuleContainer&&)                 = default;
    RuleContainer& operator=(RuleContainer const&) = default;
    RuleContainer& operator=(RuleContainer&&)      = default;

    ~RuleContainer() = default;

    void resize(size_t m) {
      // TODO(0) should we use resize?
      _used_slots.assign(m, UNDEFINED);
      if (m > 0) {
        _used_slots[0] = 0;
      }
      _words.assign(m, word_type());
    }

    word_type& next_rule_word(size_t num_edges) {
      return _words[used_slots(num_edges)++];
    }

    auto begin(size_t) const noexcept {
      return _words.begin();
    }

    auto end(size_t num_edges) noexcept {
      return _words.begin() + used_slots(num_edges);
    }

    void backtrack(size_t num_edges) {
      std::fill(_used_slots.begin() + num_edges,
                _used_slots.end(),
                size_t(UNDEFINED));
      if (num_edges == 0) {
        _used_slots[0] = 0;
      }
    }

   private:
    size_t& used_slots(size_t num_edges) {
      LIBSEMIGROUPS_ASSERT(num_edges < _used_slots.size());
      if (_used_slots[0] == UNDEFINED) {
        _used_slots[0] = 0;
      }
      size_t i = num_edges;
      while (_used_slots[i] == UNDEFINED) {
        --i;
      }
      LIBSEMIGROUPS_ASSERT(i < _used_slots.size());
      LIBSEMIGROUPS_ASSERT(_used_slots[i] != UNDEFINED);
      for (size_t j = i + 1; j <= num_edges; ++j) {
        _used_slots[j] = _used_slots[i];
      }
      LIBSEMIGROUPS_ASSERT(_used_slots[num_edges] <= _words.size());
      return _used_slots[num_edges];
    }
  };

  ///////////////////////////////////////////////////////////////////////////////
  // Sims2::iterator_base nested class
  ///////////////////////////////////////////////////////////////////////////////

  void Sims2::iterator_base::partial_copy_for_steal_from(
      Sims2::iterator_base const& that) {
    SimsBase::IteratorBase::partial_copy_for_steal_from(that);
    *_2_sided_include = *that._2_sided_include;
    _2_sided_words    = that._2_sided_words;
  }

  Sims2::iterator_base::iterator_base(Sims2 const* s, size_type n)
      : SimsBase::IteratorBase(s, n),
        // protected
        _2_sided_include(new RuleContainer()),
        _2_sided_words() {
    // TODO(2) could be slightly less space allocated here
    size_t const m = SimsBase::IteratorBase::maximum_number_of_classes();
    Presentation<word_type> const& p = this->_felsch_graph.presentation();
    _2_sided_include->resize(2 * m * p.alphabet().size());
    _2_sided_words.assign(n + 1, word_type());
    // +1 to the number of 2-sided words because if we the defining
    // presentation is a semigroup presentation, then we will use 1 more slot
    // here (i.e. the number of nodes in the resulting graphs is n + 1.
  }

  Sims2::iterator_base::iterator_base(Sims2::iterator_base const& that)
      : SimsBase::IteratorBase(that),
        _2_sided_include(new RuleContainer(*that._2_sided_include)),
        _2_sided_words(that._2_sided_words) {}

  Sims2::iterator_base::iterator_base(Sims2::iterator_base&& that)
      : SimsBase::IteratorBase(std::move(that)),  // TODO(0) std::move correct?
        _2_sided_include(std::move(that._2_sided_include)),
        _2_sided_words(std::move(that._2_sided_words)) {}

  typename Sims2::iterator_base&
  Sims2::iterator_base::operator=(Sims2::iterator_base const& that) {
    SimsBase::IteratorBase::operator=(that);
    *_2_sided_include = *that._2_sided_include;
    _2_sided_words    = that._2_sided_words;

    return *this;
  }

  typename Sims2::iterator_base&
  Sims2::iterator_base::operator=(Sims2::iterator_base&& that) {
    SimsBase::IteratorBase::operator=(
        std::move(that));  // TODO(0) std::move ok?
    _2_sided_include = std::move(that._2_sided_include);
    _2_sided_words   = std::move(that._2_sided_words);
    return *this;
  }

  void Sims2::iterator_base::swap(Sims2::iterator_base& that) noexcept {
    SimsBase::IteratorBase::swap(that);
    std::swap(_2_sided_include, that._2_sided_include);
    std::swap(_2_sided_words, that._2_sided_words);
  }

  Sims2::iterator_base::~iterator_base() = default;

  bool Sims2::iterator_base::try_define(PendingDef const& current) {
    LIBSEMIGROUPS_ASSERT(current.target < current.num_nodes);
    LIBSEMIGROUPS_ASSERT(current.num_nodes <= maximum_number_of_classes());

    if (!SimsBase::IteratorBase::try_define(current)) {
      return false;
    }

    std::lock_guard lg(_mtx);

    size_type start     = _felsch_graph.definitions().size();
    size_type num_nodes = _felsch_graph.number_of_active_nodes();
    auto      first     = _sims1or2->include().cbegin();
    auto      last      = _sims1or2->include().cend();

    // That the graph is compatible at 0 with include is checked in
    // Sims1::iterator_base::try_define.
    if (num_nodes > 1
        && (!detail::felsch_graph::make_compatible<detail::RegisterDefs>(
                _felsch_graph, 1, num_nodes, first, last)
            || !_felsch_graph.process_definitions(start))) {
      // Seems to be important to check include() first then
      // process_definitions
      return false;
    }

    _2_sided_include->backtrack(current.num_edges);

    if (current.target_is_new_node) {
      LIBSEMIGROUPS_ASSERT(current.target < _2_sided_words.size());
      LIBSEMIGROUPS_ASSERT(current.source < _2_sided_words.size());
      _2_sided_words[current.target] = _2_sided_words[current.source];
      _2_sided_words[current.target].push_back(current.generator);
    } else {
      auto const& e = _felsch_graph.definitions()[current.num_edges];
      word_type&  u = _2_sided_include->next_rule_word(current.num_edges);
      u.assign(_2_sided_words[e.first].cbegin(),
               _2_sided_words[e.first].cend());
      u.push_back(e.second);
      word_type&       v = _2_sided_include->next_rule_word(current.num_edges);
      word_type const& w
          = _2_sided_words[_felsch_graph.target_no_checks(e.first, e.second)];
      v.assign(w.begin(), w.end());
    }

    start = current.num_edges;
    while (start < _felsch_graph.definitions().size()) {
      auto first = _2_sided_include->begin(current.num_edges);
      auto last  = _2_sided_include->end(current.num_edges);
      start      = _felsch_graph.definitions().size();
      if (!detail::felsch_graph::make_compatible<detail::RegisterDefs>(
              _felsch_graph,
              0,
              _felsch_graph.number_of_active_nodes(),
              first,
              last)
          || !_felsch_graph.process_definitions(start)) {
        return false;
      }
    }
    return true;
  }

  ////////////////////////////////////////////////////////////////////////
  // RepOrc helper class
  ////////////////////////////////////////////////////////////////////////

  RepOrc& RepOrc::init() {
    _min  = 0;
    _max  = 0;
    _size = 0;
    return *this;
  }

  Sims1::word_graph_type RepOrc::word_graph() const {
    using word_graph_type = typename Sims1::word_graph_type;
    using node_type       = typename word_graph_type::node_type;

    report_no_prefix("{:+<80}\n", "");
    report_default(
        "RepOrc: Searching for a faithful rep. o.r.c. on [{}, {}) points\n",
        _min,
        _max + 1);
    if (_min > _max || _max == 0) {
      report_no_prefix("{:+<80}\n", "");
      report_default(
          "RepOrc: No faithful rep. o.r.c. exists in [{}, {}) = \u2205\n",
          _min,
          _max + 1);
      return word_graph_type(0, 0);
    }

    SuppressReportFor suppressor("FroidurePin");

    std::atomic_uint64_t count(0);

    auto hook = [&](word_graph_type const& x) {
      ++count;
      if (x.number_of_active_nodes() >= _min) {
        auto first = (presentation().contains_empty_word() ? 0 : 1);
        auto S     = to_froidure_pin<Transf<0, node_type>>(
            x, first, x.number_of_active_nodes());
        // It'd be nice to reuse S here, but this is tricky because hook
        // maybe called in multiple threads, and so we can't easily do this.
        if (presentation().contains_empty_word()) {
          auto id = one(S.generator(0));
          if (!S.contains(id)) {
            S.add_generator(id);
          }
        }
        LIBSEMIGROUPS_ASSERT(S.size() <= _size);
        if (S.size() == _size) {
          return true;
        }
      }
      return false;
    };

    Sims1 sims1;
    sims1.SimsSettings::init(*this);
    auto result = sims1.find_if(_max, hook);

    if (result.number_of_active_nodes() == 0) {
      report_default(
          "RepOrc: No faithful rep. o.r.c. on [{}, {}) points found\n",
          _min,
          _max + 1);
      result.induced_subgraph_no_checks(0, 0);
    } else {
      report_default("RepOrc: Found a faithful rep. o.r.c. on {} points\n",
                     presentation().contains_empty_word()
                         ? result.number_of_active_nodes()
                         : result.number_of_active_nodes() - 1);
      if (presentation().contains_empty_word()) {
        result.induced_subgraph_no_checks(0, result.number_of_active_nodes());
      } else {
        result.induced_subgraph_no_checks(1, result.number_of_active_nodes());
        result.number_of_active_nodes(result.number_of_active_nodes() - 1);
      }
    }
    stats_copy_from(sims1.stats());
    return result;
  }

  ////////////////////////////////////////////////////////////////////////
  // MinimalRepOrc helper class
  ////////////////////////////////////////////////////////////////////////

  // An alternative to the approach used below would be to do a sort of
  // binary search for a minimal representation. It seems that in most
  // examples that I've tried, this actually finds the minimal rep close to
  // the start of the search. The binary search approach would only really
  // be useful if the rep found at the start of the search was much larger
  // than the minimal one, and that the minimal one would not be found until
  // late in the search through all the reps with [1, best rep so far). It
  // seems that in most examples, binary search will involve checking many
  // of the same graphs repeatedly, i.e. if we check [1, size = 1000) find a
  // rep on 57 points, then check [1, 57 / 2) find nothing, then check
  // [57/2, 57) and find nothing, then the search through [57/2, 57)
  // actually iterates through all the digraphs with [1, 57 / 2) again (just
  // doesn't construct a FroidurePin object for these). So, it seems to be
  // best to just search through the digraphs with [1, 57) nodes once.
  //
  // TODO(2)(later) perhaps find minimal 2-sided congruences first (or try to)
  // and then run MinimalRepOrc for right congruences excluding all the
  // generating pairs from the minimal 2-sided congruences. Also with this
  // approach FroidurePin wouldn't be required in RepOrc. This might not work,
  // given that the minimal rc might contain some pairs from minimal 2-sided
  // congs, just not all of them.
  Sims1::word_graph_type MinimalRepOrc::word_graph() const {
    auto cr = RepOrc(*this);

    size_t hi   = (presentation().contains_empty_word() ? _size : _size + 1);
    auto   best = cr.min_nodes(1).max_nodes(hi).target_size(_size).word_graph();

    if (best.number_of_nodes() < 1) {
      stats_copy_from(cr.stats());
      return best;
    }

    hi        = best.number_of_nodes();
    auto next = cr.max_nodes(hi - 1).word_graph();
    while (next.number_of_nodes() != 0) {
      hi   = next.number_of_nodes();
      best = std::move(next);
      next = cr.max_nodes(hi - 1).word_graph();
    }
    stats_copy_from(cr.stats());
    return best;
  }

  namespace sims {
    const_rcgp_iterator::~const_rcgp_iterator() = default;

    const_rcgp_iterator::const_rcgp_iterator(Presentation<word_type> const& p,
                                             word_graph_type const*         ptr,
                                             node_type  source,
                                             label_type gen)
        :  // protected
          _reconstructed_word_graph(p),
          // private
          _gen(gen),
          _source(source),
          _relation({}, {}),
          _tree(ptr->number_of_active_nodes()),
          _word_graph(ptr) {
      // if (!word_graph::is_complete(*ptr,
      //                              ptr->cbegin_nodes(),
      //                              ptr->cbegin_nodes()
      //                                  + ptr->number_of_active_nodes())) {
      //   LIBSEMIGROUPS_EXCEPTION("The 1st argument must be a pointer to a "
      //                           "complete word graph on the nodes [0, {})",
      //                           ptr->number_of_active_nodes());
      // }
      _reconstructed_word_graph.add_nodes(ptr->number_of_active_nodes());
      ++(*this);
    }

    const_rcgp_iterator const& const_rcgp_iterator::operator++() {
      if (at_end()) {
        return *this;
      }

      _relation.first.clear();
      _relation.second.clear();

      auto const& wg = *_word_graph;

      size_type start = _reconstructed_word_graph.definitions().size();

      while (_source < wg.number_of_active_nodes()) {
        while (_gen < wg.out_degree()) {
          auto target1 = wg.target_no_checks(_source, _gen);
          if (target1 != UNDEFINED) {
            auto target2
                = _reconstructed_word_graph.target_no_checks(_source, _gen);
            LIBSEMIGROUPS_ASSERT(target1 != UNDEFINED);
            if (target2 == UNDEFINED) {
              _reconstructed_word_graph.target_no_checks(
                  _source, _gen, target1);
              if (_tree.parent_no_checks(target1) == UNDEFINED
                  && target1 > _source) {
                // Tree edges
                _tree.set_parent_and_label_no_checks(target1, _source, _gen);
              } else {
                // Non-tree edge, not implied by other edges
                _gen++;
                std::ignore
                    = _reconstructed_word_graph.process_definitions(start);
                return *this;
              }
            }
          }
          _gen++;
        }
        _gen = 0;
        _source++;
      }
      return *this;
    }

    void const_rcgp_iterator::swap(const_rcgp_iterator& that) noexcept {
      std::swap(_word_graph, that._word_graph);
      std::swap(_gen, that._gen);
      std::swap(_source, that._source);
      std::swap(_tree, that._tree);
      std::swap(_reconstructed_word_graph, that._reconstructed_word_graph);
    }

    bool const_rcgp_iterator::populate_relation() const {
      if (_relation.first.empty() && !at_end()) {
        LIBSEMIGROUPS_ASSERT(_relation.second.empty());
        LIBSEMIGROUPS_ASSERT(_gen > 0);
        _tree.path_to_root_no_checks(_relation.first, _source);
        std::reverse(_relation.first.begin(), _relation.first.end());
        _relation.first.push_back(_gen - 1);
        _tree.path_to_root_no_checks(
            _relation.second, _word_graph->target_no_checks(_source, _gen - 1));
        std::reverse(_relation.second.begin(), _relation.second.end());
        return true;
      }
      return false;
    }

    const_cgp_iterator::~const_cgp_iterator() = default;

  }  // namespace sims

  [[nodiscard]] std::string to_human_readable_repr(SimsStats const&) {
    return fmt::format("<SimsStats object>");
  }

  namespace helper {

    template <typename Subclass>
    [[nodiscard]] std::string
    sims_to_human_readable_repr_helper(SimsSettings<Subclass> const& x,
                                       std::string extra_text) {
      std::string result    = "";
      bool        needs_and = false;

      result += fmt::format("object over {} with",
                            to_human_readable_repr(x.presentation()));
      std::string comma = (x.pruners().empty()
                           || (x.pruners().size() == 1 && !x.exclude().empty()))
                                  && extra_text.empty()
                              ? ""
                              : ",";
      if ((x.include().size() > 0) && (x.exclude().size() > 0)) {
        result += fmt::format(" {} include and {} exclude pairs{}",
                              x.include().size() / 2,
                              x.exclude().size() / 2,
                              comma);
        needs_and = true;
      } else if (x.include().size() > 0) {
        result += fmt::format(" {} include pair{}{}",
                              x.include().size() / 2,
                              x.include().size() / 2 == 1 ? "" : "s",
                              comma);
        needs_and = true;
      } else if (x.exclude().size() > 0) {
        result += fmt::format(" {} exclude pair{}{}",
                              x.exclude().size() / 2,
                              x.exclude().size() / 2 == 1 ? "" : "s",
                              comma);
        needs_and = true;
      }

      result += extra_text;
      if (extra_text.size()) {
        needs_and = true;
      }

      if (x.pruners().size() > 1
          || (x.exclude().empty() && !x.pruners().empty())) {
        result += fmt::format(" {} pruner{}",
                              x.exclude().size() == 0 ? x.pruners().size()
                                                      : x.pruners().size() - 1,
                              x.pruners().size() == 1 ? "" : "s");
        needs_and = true;
      }
      result += needs_and ? " and" : "";
      result += fmt::format(" {} thread{}",
                            x.number_of_threads(),
                            x.number_of_threads() == 1 ? "" : "s");
      return result;
    }

  }  // namespace helper

  // template <typename Subclass>
  // [[nodiscard]] std::string
  // to_human_readable_repr(SimsSettings<Subclass> const& x) {
  //   return fmt::format("<SimsSettings {}>",
  //                      helper::sims_to_human_readable_repr_helper(x, ""));
  // }

  [[nodiscard]] std::string to_human_readable_repr(Sims1 const& x) {
    return fmt::format("<Sims1 {}>",
                       helper::sims_to_human_readable_repr_helper(x, ""));
  }

  [[nodiscard]] std::string to_human_readable_repr(Sims2 const& x) {
    return fmt::format("<Sims2 {}>",
                       helper::sims_to_human_readable_repr_helper(x, ""));
  }

  [[nodiscard]] std::string to_human_readable_repr(RepOrc const& x) {
    return fmt::format("<RepOrc {}>",
                       helper::sims_to_human_readable_repr_helper(
                           x,
                           fmt::format(" node bounds [{}, {}), target size {}",
                                       x.min_nodes(),
                                       x.max_nodes(),
                                       x.target_size())));
  }

  [[nodiscard]] std::string to_human_readable_repr(MinimalRepOrc const& x) {
    return fmt::format("<MinimalRepOrc {}>",
                       helper::sims_to_human_readable_repr_helper(
                           x, fmt::format(" target size {}", x.target_size())));
  }

  [[nodiscard]] std::string to_human_readable_repr(SimsRefinerIdeals const& x) {
    return fmt::format("<SimsRefinerIdeals object over presentation {}>",
                       to_human_readable_repr(x.presentation()));
  }

  [[nodiscard]] std::string
  to_human_readable_repr(SimsRefinerFaithful const& x) {
    return fmt::format("<SimsRefinerFaithful object with {} forbidden pair{}>",
                       x.forbid().size() / 2,
                       x.forbid().size() / 2 == 1 ? "" : "s");
  }

}  // namespace libsemigroups
