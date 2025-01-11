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

// This file contains the declaration of the Race class template for
// competitively running different functions/methods in different threads, and
// obtaining the winner.

// TODO(later):
// 1. consider if keeping killed off methods has any uses
// 2. update run_until to be similar to Runner::run_until
// 3. add tpp file

#ifndef LIBSEMIGROUPS_DETAIL_RACE_HPP_
#define LIBSEMIGROUPS_DETAIL_RACE_HPP_

#include <chrono>       // for nanoseconds
#include <cstddef>      // for size_t
#include <memory>       // for std::shared_ptr
#include <thread>       // for mutex
#include <type_traits>  // for invoke_result_t
#include <vector>       // for vector

#include "libsemigroups/constants.hpp"  // for UNDEFINED
#include "libsemigroups/debug.hpp"      // for LIBSEMIGROUPS_ASSERT
#include "libsemigroups/exception.hpp"  // for LIBSEMIGROUPS_EXCEPTION
#include "libsemigroups/runner.hpp"     // for Runner

#include "report.hpp"  // for report_default, REPORT_TIME
#include "stl.hpp"     // for IsCallable
#include "timer.hpp"   // for Timer

namespace libsemigroups {
  namespace detail {
    class Race {
      size_t                               _max_threads;
      std::mutex                           _mtx;
      std::vector<std::shared_ptr<Runner>> _runners;
      std::shared_ptr<Runner>              _winner;
      size_t                               _winner_index;

     public:
      // Construct an empty Race object, with maximum number of threads set to
      // std::thread::hardware_concurrency.
      Race();

      // TODO(1) to cpp
      Race& init() {
        _max_threads = std::thread::hardware_concurrency();
        // do nothing to the _mtx
        _runners.clear();
        _winner       = nullptr;
        _winner_index = UNDEFINED;
        return *this;
      }

      Race(Race const& other) : Race() {
        // Can't use = default because std::mutex is non-copyable.
        _max_threads = other._max_threads;
        // do nothing to the _mtx
        _runners      = other._runners;
        _winner       = other._winner;
        _winner_index = other._winner_index;
      }

      // TODO(1) to cpp
      Race(Race&& other)
          : _max_threads(std::move(other._max_threads)),
            _mtx(),
            _runners(std::move(other._runners)),
            _winner(std::move(other._winner)),
            _winner_index(std::move(other._winner_index)) {}

      // TODO(1) to cpp
      Race& operator=(Race const& other) {
        // Can't use = default because std::mutex is non-copyable.
        _max_threads = other._max_threads;
        // do nothing to the _mtx
        _runners      = other._runners;
        _winner       = other._winner;
        _winner_index = other._winner_index;
        return *this;
      }

      // TODO(1) to cpp
      Race& operator=(Race&& other) {
        // Can't use = default because std::mutex is non-copyable.
        _max_threads = std::move(other._max_threads);
        // do nothing to the _mtx
        _runners      = std::move(other._runners);
        _winner       = std::move(other._winner);
        _winner_index = std::move(other._winner_index);
        return *this;
      }

      ~Race();

      // Set the maximum number of threads, throws if try to set to 0.
      Race& max_threads(size_t val) noexcept {
        LIBSEMIGROUPS_ASSERT(val != 0);
        _max_threads = val;
        return *this;
      }

      [[nodiscard]] size_t max_threads() const noexcept {
        return _max_threads;
      }

      // Runs the method Runner::run on every Runner in the Race, and returns
      // the one that finishes first. The losers are deleted.
      [[nodiscard]] std::shared_ptr<Runner> winner() {
        run();
        return _winner;
      }

      [[nodiscard]] size_t winner_index() const {
        return _winner_index;
      }

      [[nodiscard]] size_t winner_index() {
        run();
        return _winner_index;
      }

      [[nodiscard]] bool finished() const noexcept {
        return _winner != nullptr;
      }

      // Adds a Runner to the race, throws if the race is already over.
      void add_runner(std::shared_ptr<Runner>);

      using const_iterator =
          typename std::vector<std::shared_ptr<Runner>>::const_iterator;

      [[nodiscard]] const_iterator begin() const noexcept {
        return _runners.cbegin();
      }

      [[nodiscard]] const_iterator end() const noexcept {
        return _runners.cend();
      }

      // Returns an iterator pointing to the first Runner in the Race.
      [[nodiscard]] const_iterator cbegin() const noexcept {
        return _runners.cbegin();
      }

      // Returns an iterator pointing to one past the last Runner in the Race.
      [[nodiscard]] const_iterator cend() const noexcept {
        return _runners.cend();
      }

      // Returns \c true if there are no Runners in the race, and \c false
      // otherwise.
      // std::vector::empty is noexcept
      [[nodiscard]] bool empty() const noexcept {
        return _runners.empty();
      }

      // std::vector::size is noexcept
      [[nodiscard]] size_t number_of_runners() const noexcept {
        return _runners.size();
      }

      // Runs the race to completion.
      void run();

      // Runs the race for the specified amount of time.
      void run_for(std::chrono::nanoseconds);

      // Runs until \p func returns \c true, or the race is over. This
      // repeatedly calls Race::run_for for \p check_interval, and then checks
      // whether or not \p func() returns true. The object \p func can be any
      // callable object with 0 parameters and that returns a bool.
      // This is definitely tested but doesn't show up in the code coverage for
      // some reason.
      template <typename TCallable>
      void run_until(TCallable const&         func,
                     std::chrono::nanoseconds check_interval
                     = std::chrono::milliseconds(2)) {
        static_assert(detail::IsCallable<TCallable>::value,
                      "the template parameter TCallable must be callable");
        static_assert(
            std::is_same<std::invoke_result_t<TCallable>, bool>::value,
            "the template parameter TCallable must return a bool");
        if (empty()) {
          LIBSEMIGROUPS_EXCEPTION("no runners given, cannot run_until");
        }
        while (!func() && _winner == nullptr) {
          // if _winner != nullptr, then the race is over.
          run_for(check_interval);
          check_interval *= 2;
        }
      }

      template <typename T>
      bool has() const {
        return find_runner<T>() != nullptr;
      }

      template <typename T>
      [[nodiscard]] std::shared_ptr<T> find_runner() const {
        static_assert(std::is_base_of<Runner, T>::value,
                      "the template parameter must be derived from Runner");
        // We use find_if so that this works even if we haven't computed
        // anything at all.
        auto it = std::find_if(_runners.begin(),
                               _runners.end(),
                               [](std::shared_ptr<Runner> const& m) {
                                 auto& r = *(m.get());
                                 return typeid(r) == typeid(T);
                               });
        if (it != _runners.end()) {
          return std::static_pointer_cast<T>(*it);
        } else {
          return nullptr;
        }
      }

      void erase_runners(const_iterator pos) {
        _runners.erase(pos);
      }

      void erase_runners(const_iterator first, const_iterator last) {
        _runners.erase(first, last);
      }

     private:
      // Runs the callable object \p func on every Runner in parallel.
      template <typename TCallable>
      void run_func(TCallable const& func) {
        static_assert(
            std::is_same<
                std::invoke_result_t<TCallable, std::shared_ptr<Runner>>,
                void>::value,
            "the result of calling the template parameter TCallable "
            "must be void");
        LIBSEMIGROUPS_ASSERT(!empty());
        if (_winner == nullptr) {
          size_t nr_threads = std::min(_runners.size(), _max_threads);
          if (nr_threads == 1) {
            report_default("using 0 additional threads\n");
            detail::Timer tmr;
            func(_runners.at(0));
            _winner       = _runners.at(0);
            _winner_index = 0;
            report_elapsed_time("Race: ", tmr);
            return;
          }
          for (size_t i = 0; i < _runners.size(); ++i) {
            if (_runners[i]->success()) {
              report_default("using 0 additional threads\n");
              _winner       = _runners[i];
              _winner_index = i;
              report_default("#{} already finished successfully!\n", i);
              // delete the other runners?
              return;
            }
          }

          std::vector<std::thread::id> tids(_runners.size(),
                                            std::this_thread::get_id());

          report_default("using {} / {} additional threads\n",
                         nr_threads,
                         std::thread::hardware_concurrency());
          detail::Timer tmr;
          LIBSEMIGROUPS_ASSERT(nr_threads != 0);

          auto thread_func = [this, &func, &tids](size_t pos) {
            tids[pos] = std::this_thread::get_id();
            try {
              func(_runners.at(pos));
            } catch (std::exception const& e) {
              size_t tid = thread_id(tids[pos]);
              report_default("exception thrown by #{}:\n{}\n", tid, e.what());
              return;
            }
            // Stop two Runner* objects from killing each other
            {
              std::lock_guard<std::mutex> lg(_mtx);
              if (_runners.at(pos)->success()) {
                for (auto it = _runners.begin(); it < _runners.begin() + pos;
                     it++) {
                  (*it)->kill();
                }
                for (auto it = _runners.begin() + pos + 1; it < _runners.end();
                     it++) {
                  (*it)->kill();
                }
              }
            }
          };
          // TODO remove the next line?
          reset_thread_ids();

          std::vector<std::thread> t;
          for (size_t i = 0; i < nr_threads; ++i) {
            t.push_back(std::thread(thread_func, i));
          }
          for (size_t i = 0; i < nr_threads; ++i) {
            t.at(i).join();
          }
          report_elapsed_time("Race: ", tmr);
          report_no_prefix("\n");
          for (auto method = _runners.begin(); method < _runners.end();
               ++method) {
            if ((*method)->success()) {
              LIBSEMIGROUPS_ASSERT(_winner == nullptr);
              _winner       = *method;
              _winner_index = method - _runners.begin();
              size_t tid    = thread_id(tids.at(method - _runners.begin()));
              report_default("#{} is the winner!\n", tid);
              break;
            }
          }
          if (_winner != nullptr) {
            for (auto rnnr : _runners) {
              if (rnnr != _winner) {
                rnnr.reset();
              }
            }
            _runners.clear();
            _runners.push_back(_winner);
          }
        }
      }
    };
  }  // namespace detail
}  // namespace libsemigroups

#endif  // LIBSEMIGROUPS_DETAIL_RACE_HPP_
