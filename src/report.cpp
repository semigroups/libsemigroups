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

#include "libsemigroups/report.hpp"

#include <ratio>          // for ratio, nano
#include <unordered_map>  // for operator!=

#include "fmt/format.h"  // for everything

#include "libsemigroups/config.hpp"  // for LIBSEMIGROUPS_FMT_ENABLED
#include "libsemigroups/debug.hpp"   // for LIBSEMIGROUPS_ASSERT
#include "libsemigroups/string.hpp"  // for unicode_string_length, wrap

namespace libsemigroups {
  detail::Reporter        REPORTER;
  detail::ThreadIdManager THREAD_ID_MANAGER;

  namespace detail {

    ThreadIdManager::ThreadIdManager() : _mtx(), _next_tid(0), _thread_map() {
      tid(std::this_thread::get_id());
    }

    void ThreadIdManager::reset() {
      // Only do this from the main thread
      LIBSEMIGROUPS_ASSERT(tid(std::this_thread::get_id()) == 0);
      // Delete all thread_ids
      _thread_map.clear();
      _next_tid = 0;
      // Reinsert the main thread's id
      tid(std::this_thread::get_id());
    }

    size_t ThreadIdManager::tid(std::thread::id t) {
      std::lock_guard<std::mutex> lg(_mtx);
      auto                        it = _thread_map.find(t);
      if (it != _thread_map.end()) {
        return (*it).second;
      } else {
        // Don't check the assert below because on a single thread machine
        // (such as those used by appveyor), for an fp-semigroup more than 1
        // thread will be used, and this assertion will fail.
        // LIBSEMIGROUPS_ASSERT(_next_tid <=
        // std::thread::hardware_concurrency());
        _thread_map.emplace(t, _next_tid++);
        return _next_tid - 1;
      }
    }

    Reporter::Reporter(bool report)
        : _last_msg(), _mtx(), _msg(), _options(), _report(report) {}

#ifdef LIBSEMIGROUPS_FMT_ENABLED
    Reporter& Reporter::color(fmt::color c) {
      if (_report) {
        size_t tid = THREAD_ID_MANAGER.tid(std::this_thread::get_id());
        resize(tid + 1);
        _options[tid].color = c;
      }
      return *this;
    }

    Reporter& Reporter::thread_color() {
      if (_report) {
        std::lock_guard<std::mutex> lg(_mtx);
        size_t tid = THREAD_ID_MANAGER.tid(std::this_thread::get_id());
        resize(tid + 1);
        _options[tid].color = thread_colors[tid % thread_colors.size()];
      }
      return *this;
    }
#endif

    Reporter& Reporter::prefix() {
      if (_report) {
        std::lock_guard<std::mutex> lg(_mtx);
        size_t tid = THREAD_ID_MANAGER.tid(std::this_thread::get_id());
        resize(tid + 1);
        _options[tid].prefix = "";
      }
      return *this;
    }

    Reporter& Reporter::flush_right() {
      if (_report) {
        std::lock_guard<std::mutex> lg(_mtx);
        size_t tid = THREAD_ID_MANAGER.tid(std::this_thread::get_id());
        resize(tid + 1);
        _options[tid].flush_right = true;
      }
      return *this;
    }

    void Reporter::flush() {
      if (_report) {
        std::lock_guard<std::mutex> lg(_mtx);
        size_t      tid  = THREAD_ID_MANAGER.tid(std::this_thread::get_id());
        auto const& prfx = _options[tid].prefix;
        auto        pos  = prfx.find(":");
        if (pos != std::string::npos && pos + 2 <= prfx.size() - 1) {
          std::string class_name(prfx.cbegin() + prfx.find(":") + 2,
                                 prfx.cend() - 2);
          if (_suppressions.find(class_name) != _suppressions.cend()) {
            return;
          }
        }

        size_t pad = 0;
        _msg[tid]  = _options[tid].prefix + _msg[tid];
        if (_options[tid].flush_right
            && _last_msg[tid].size() + unicode_string_length(_msg[tid]) < 80) {
          pad = (80 - _last_msg[tid].size()) - unicode_string_length(_msg[tid]);
          _msg[tid] = std::string(pad, ' ') + _msg[tid];
        }
#ifdef LIBSEMIGROUPS_VERBOSE
        if (_msg[tid].back() != '\n') {
          _msg[tid] += "\n";
        }
#endif
        _msg[tid] = wrap(_options[tid].prefix.length(), _msg[tid]);
#ifdef LIBSEMIGROUPS_FMT_ENABLED
        fmt::print(fg(_options[tid].color), _msg[tid]);
#else
        std::cout << _msg[tid];
#endif
        _options[tid] = Options();
      }
    }

    void Reporter::resize(size_t n) {
      if (n > _msg.size()) {
        _last_msg.resize(n);
        _msg.resize(n);
        _options.resize(n);
      }
    }

    bool string_time_incremental(std::string&              result,
                                 std::chrono::nanoseconds& elapsed,
                                 bool                      use_float) {
      using seconds = std::chrono::seconds;
      seconds x     = std::chrono::duration_cast<seconds>(elapsed);
      if (x.count() > 0) {
        if (use_float) {
          double x_float = static_cast<double>(elapsed.count()) / 1'000'000'000;
          result += fmt::format("{:.3f}s", x_float);
        } else {
          result += fmt::format("{}", x);
        }
        elapsed -= std::chrono::nanoseconds(x);
        return true;
      }
      return false;
    }
  }  // namespace detail

  namespace report {
    bool should_report() noexcept {
      return REPORTER.report();
    }

    void suppress(std::string const& class_name) {
      REPORTER.suppress(class_name);
    }

    void clear_suppressions() {
      REPORTER.clear_suppressions();
    }

  }  // namespace report

}  // namespace libsemigroups
