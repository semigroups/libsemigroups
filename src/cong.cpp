//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2019 James D. Mitchell
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

// This file contains stuff for creating congruence over FroidurePinBase
// objects or over Presentation objects.

#include "libsemigroups/cong.hpp"

#include <memory>  // for shared_ptr

#include "libsemigroups/debug.hpp"              // for LIBSEMIGROUPS_ASSERT
#include "libsemigroups/exception.hpp"          // for LIBSEMIGROUPS_EXCEPTION
#include "libsemigroups/froidure-pin-base.hpp"  // for FroidurePinBase
#include "libsemigroups/kambites.hpp"           // for Kambites
#include "libsemigroups/knuth-bendix.hpp"       // for KnuthBendix_
#include "libsemigroups/to-presentation.hpp"    // for to_presentation
#include "libsemigroups/to-todd-coxeter.hpp"    // for to_todd_coxeter
#include "libsemigroups/todd-coxeter.hpp"       // for ToddCoxeter
#include "libsemigroups/types.hpp"              // for word_type

namespace libsemigroups {

  //////////////////////////////////////////////////////////////////////////
  // Congruence - constructors - public
  //////////////////////////////////////////////////////////////////////////

  Congruence& Congruence::init() {
    CongruenceInterface::init();
    _race.init();
    _runners_initted = false;
    _runner_kinds.clear();
    return *this;
  }

  Congruence& Congruence::init(congruence_kind                type,
                               Presentation<word_type> const& p) {
    CongruenceInterface::init(type);
    _race.max_threads(POSITIVE_INFINITY);
    if (type == congruence_kind::twosided) {
      add_runner(std::make_shared<Kambites<word_type>>(type, p));
    }
    add_runner(std::make_shared<KnuthBendix<>>(type, p));

    add_runner(std::make_shared<ToddCoxeter>(type, p));
    auto tc = std::make_shared<ToddCoxeter>(type, p);
    tc->strategy(ToddCoxeter::options::strategy::felsch);
    add_runner(std::move(tc));
    return *this;
  }

  uint64_t Congruence::number_of_classes() {
    run();
    auto winner_kind = _runner_kinds[_race.winner_index()];
    if (winner_kind == RunnerKind::TC) {
      return std::static_pointer_cast<ToddCoxeter>(_race.winner())
          ->number_of_classes();
    } else if (winner_kind == RunnerKind::KB) {
      return std::static_pointer_cast<KnuthBendix<>>(_race.winner())
          ->number_of_classes();
    } else {
      LIBSEMIGROUPS_ASSERT(winner_kind == RunnerKind::K);
      return std::static_pointer_cast<Kambites<word_type>>(_race.winner())
          ->number_of_classes();
    }
  }

  void Congruence::init_runners() const {
    if (!_runners_initted) {
      _runners_initted = true;
      for (auto const& [i, runner] : rx::enumerate(_race)) {
        auto first = generating_pairs().cbegin();
        auto last  = generating_pairs().cend();
        if (_runner_kinds[i] == RunnerKind::TC) {
          for (auto it = first; it != last; it += 2) {
            // We have to call the specific add_generating_pair_no_checks for
            // each type of runner so that the generating pairs are correctly
            // modified to match the native letters in the alphabet
            std::static_pointer_cast<ToddCoxeter>(runner)
                ->add_generating_pair_no_checks(std::begin(*it),
                                                std::end(*it),
                                                std::begin(*(it + 1)),
                                                std::end(*(it + 1)));
          }
        } else if (_runner_kinds[i] == RunnerKind::KB) {
          for (auto it = first; it != last; it += 2) {
            std::static_pointer_cast<KnuthBendix<>>(runner)
                ->add_generating_pair_no_checks(std::begin(*it),
                                                std::end(*it),
                                                std::begin(*(it + 1)),
                                                std::end(*(it + 1)));
          }
        } else {
          LIBSEMIGROUPS_ASSERT(_runner_kinds[i] == RunnerKind::K);
          for (auto it = first; it != last; it += 2) {
            std::static_pointer_cast<Kambites<word_type>>(runner)
                ->add_generating_pair_no_checks(std::begin(*it),
                                                std::end(*it),
                                                std::begin(*(it + 1)),
                                                std::end(*(it + 1)));
          }
        }
      }
    }
  }

  void Congruence::run_impl() {
    init_runners();
    _race.run_until([this]() { return this->stopped(); });
  }

  Presentation<word_type> const& Congruence::presentation() const {
    if (!_race.empty()) {
      if (finished()) {
        size_t index = _race.winner_index();

        if (_runner_kinds[index] == RunnerKind::TC) {
          return std::static_pointer_cast<ToddCoxeter>(_race.winner())
              ->presentation();
        } else if (_runner_kinds[index] == RunnerKind::K) {
          return std::static_pointer_cast<Kambites<word_type>>(_race.winner())
              ->presentation();
        }
        // KnuthBendix uses a std::string Presentation and this can't be
        // returned by reference here. There are probably better ways of
        // handling this:
        // 1. store a copy in the congruence object itself
        // 2. return by value
        // 3. Something else?
      } else {
        auto it = std::find_if(
            _runner_kinds.begin(), _runner_kinds.end(), [](auto const& val) {
              return val != RunnerKind::KB;
            });
        if (it != _runner_kinds.end()) {
          size_t index = std::distance(_runner_kinds.begin(), it);
          if (_runner_kinds[index] == RunnerKind::TC) {
            return std::static_pointer_cast<ToddCoxeter>(
                       *(_race.begin() + index))
                ->presentation();
          } else if (_runner_kinds[index] == RunnerKind::K) {
            return std::static_pointer_cast<Kambites<word_type>>(
                       *(_race.begin() + index))
                ->presentation();
          }
        }
      }
    }
    LIBSEMIGROUPS_EXCEPTION(
        "No presentation has been set, or it cannot be returned!");
  }

  namespace congruence {
    // We have to pass the presentation here, because o/w we have no way of
    // knowing over what we should compute the non-trivial classes (i.e. we
    // cannot always recover p from cong).
    std::vector<std::vector<word_type>>
    non_trivial_classes(Congruence& cong, Presentation<word_type> const& p) {
      using rx::operator|;
      cong.run();
      if (cong.has<ToddCoxeter>() && cong.get<ToddCoxeter>()->finished()) {
        ToddCoxeter tc(cong.kind(), p);
        return ::libsemigroups::todd_coxeter::non_trivial_classes(
            *cong.get<ToddCoxeter>(), tc);
      } else if (cong.has<KnuthBendix<>>()
                 && cong.get<KnuthBendix<>>()->finished()) {
        KnuthBendix<> kb(cong.kind(), p);
        auto strings = ::libsemigroups::knuth_bendix::non_trivial_classes(
            kb, *cong.get<KnuthBendix<>>());
        std::vector<std::vector<word_type>> result;
        for (auto const& klass : strings) {
          result.push_back(rx::iterator_range(klass.begin(), klass.end())
                           | ToWord(kb.presentation().alphabet())
                           | rx::to_vector());
        }
        return result;
      }
      // If a Kambites object wins the race in Congruence, then we cannot
      // really compute anything here unless the semigroup defined by p is
      // finite, but that'd be better handled explicitly in any code calling
      // this.
      LIBSEMIGROUPS_EXCEPTION("Cannot compute the non-trivial classes!");
    }

    std::vector<std::vector<std::string>>
    non_trivial_classes(Congruence& cong, Presentation<std::string> const& p) {
      using rx::operator|;
      cong.run();
      if (cong.has<KnuthBendix<>>() && cong.get<KnuthBendix<>>()->finished()) {
        KnuthBendix<> kb(cong.kind(), p);
        return ::libsemigroups::knuth_bendix::non_trivial_classes(
            kb, *cong.get<KnuthBendix<>>());
      }
      if (cong.has<ToddCoxeter>() && cong.get<ToddCoxeter>()->finished()) {
        ToddCoxeter tc(cong.kind(), p);
        auto        words = ::libsemigroups::todd_coxeter::non_trivial_classes(
            *cong.get<ToddCoxeter>(), tc);
        std::vector<std::vector<std::string>> result;
        for (auto const& klass : words) {
          result.push_back(rx::iterator_range(klass.begin(), klass.end())
                           | ToString(p.alphabet()) | rx::to_vector());
        }
        return result;
      }

      // If a Kambites object wins the race in Congruence, then we cannot
      // really compute anything here unless the semigroup defined by p is
      // finite, but that'd be better handled explicitly in any code calling
      // this.
      LIBSEMIGROUPS_EXCEPTION("Cannot compute the non-trivial classes!");
    }
  }  // namespace congruence
}  // namespace libsemigroups
