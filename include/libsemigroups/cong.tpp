//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2024 James D. Mitchell
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

namespace libsemigroups {
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

    // auto tc = to_todd_coxeter(knd, S, wg);
    // tc.strategy(ToddCoxeter::options::strategy::felsch);
    // add_runner(std::make_shared<ToddCoxeter>(std::move(tc)));

    auto tc = ToddCoxeter(knd, to_presentation<word_type>(S));
    add_runner(std::make_shared<ToddCoxeter>(std::move(tc)));

    // tc = ToddCoxeter(knd, to_presentation<word_type>(S));
    // tc.strategy(ToddCoxeter::options::strategy::felsch);
    // add_runner(std::make_shared<ToddCoxeter>(std::move(tc)));

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

  template <typename Iterator1, typename Iterator2>
  void Congruence::throw_if_letter_out_of_bounds(Iterator1 first,
                                                 Iterator2 last) const {
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

  namespace congruence {
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

    template <typename Range, typename Word, typename>
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

    // The following doesn't work, because the types of the two returned values
    // aren't the same.
    // TODO(1) implement a class containing a variant for this.
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

  }  // namespace congruence
}  // namespace libsemigroups
