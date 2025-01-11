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

#define CATCH_CONFIG_ENABLE_ALL_STRINGMAKERS
#define CATCH_CONFIG_FAST_COMPILE
#define CATCH_CONFIG_MAIN

#include <algorithm>  // for transform, find_if
#include <cctype>     // for tolower
#include <chrono>     // for nanoseconds
#include <cmath>      // for pow
#include <cstdlib>    // for exit, size_t
#include <string>     // for char_traits, allocator, hash
#include <vector>     // for vector

#include "catch_amalgamated.hpp"  // for Colour, Colour::Code::BrightRed, CATCH_REGISTER...

#include "libsemigroups/detail/fmt.hpp"
#include "libsemigroups/detail/string.hpp"  // for to_string, unicode_string_length
#include "libsemigroups/detail/timer.hpp"  // for Timer

namespace {
  bool find_tag(Catch::TestCaseInfo const& testInfo, std::string tag) {
    std::transform(tag.begin(), tag.end(), tag.begin(), ::tolower);

    return std::find_if(testInfo.tags.cbegin(),
                        testInfo.tags.cend(),
                        [&tag](Catch::Tag const& t) -> bool {
                          auto t_copy = std::string(t.original);
                          std::transform(t_copy.begin(),
                                         t_copy.end(),
                                         t_copy.begin(),
                                         ::tolower);
                          return t_copy == tag;
                        })
           != testInfo.tags.cend();
  }
}  // namespace

struct LibsemigroupsListener : Catch::EventListenerBase {
  using EventListenerBase::EventListenerBase;  // inherit constructor

  // Catch::TestCaseInfo is non-copyable so we have our own version here.
  struct TestCaseInfo {
    std::string category;
    std::string name;
    std::string number;

    TestCaseInfo() = default;

    explicit TestCaseInfo(Catch::TestCaseInfo const& testInfo)
        : TestCaseInfo() {
      *this = testInfo;
    }

    TestCaseInfo& operator=(Catch::TestCaseInfo const& testInfo) {
      name = testInfo.name;
      for (auto const& tag : testInfo.tags) {
        if (tag.original.size() == 3
            && std::all_of(tag.original.begin(),
                           tag.original.end(),
                           [](auto const& c) { return std::isdigit(c); })) {
          number = std::string(tag.original);
          break;
        }
      }

      for (std::string const& cat : {"quick", "standard", "extreme", "fail"}) {
        if (find_tag(testInfo, cat)) {
          category = cat;
          break;
        }
      }
      return *this;
    }
  };

  // Catch::SectionInfo is non-copyable so we have our own version here.
  struct SectionInfo {
    std::string name;

    SectionInfo& operator=(Catch::SectionInfo const& si) {
      name = si.name;
      return *this;
    }
  };

  // Catch::SectionStats is non-copyable so we have our own version here.
  struct SectionStats {
    std::chrono::nanoseconds duration;
    std::string              name;
    // We could store the SectionStats.sectionInfo here instead of "name", but
    // the name we want is the most nested one, and the
    // SectionStats.sectionInfo.name is the least nested one.

    SectionStats(Catch::SectionStats const& ss, SectionInfo const& si)
        : duration(std::chrono::nanoseconds(
              static_cast<uint64_t>(ss.durationInSeconds * std::pow(10, 9)))),
          name(si.name) {}
  };

  void print_extreme_test_divider(std::string_view sv) const {
    auto msg  = fmt::format(_extreme_emph,
                           "[{}]: {} - {}\n",
                           current_test_case_info().number,
                           current_test_case_info().name,
                           sv);
    auto rule = fmt::format(
        _extreme_emph, "{:=>{}}\n", "", std::max(_line_cols, msg.size() - 9));
    fmt::print(rule + msg + rule);
  }

  // Not currently used
  // void check_category(Catch::TestCaseInfo const& testInfo) {
  //     for (std::string const& cat : {"quick", "standard", "extreme", "fail"})
  //     {
  //       if (find_tag(testInfo, cat)) {
  //         category = cat;
  //         break;
  //       }
  //     }
  //   if (!(find_tag(testInfo, "quick") || find_tag(testInfo, "standard")
  //         || find_tag(testInfo, "extreme") || find_tag(testInfo, "fail"))) {
  //     {
  //       Catch::Colour colourGuard(Catch::Colour::BrightRed);
  //       Catch::cerr()
  //           << "Missing category tag: [quick|standard|extreme|fail]!\n"
  //           << "  in test case at " << testInfo.lineInfo.file << ":"
  //           << testInfo.lineInfo.line << "\n";
  //     }
  //     std::exit(1);
  //   }
  // }

  void set_current_test_case_info(Catch::TestCaseInfo const& testInfo) {
    _current_test_case_info = testInfo;
    _section_depth          = 0;
  }

  TestCaseInfo const& current_test_case_info() const noexcept {
    return _current_test_case_info;
  }

  void set_most_recent_section_info(Catch::SectionInfo const& sectionInfo) {
    _most_recent_section_info = sectionInfo;
  }

  SectionInfo const& most_recent_section_info() const {
    return _most_recent_section_info;
  }

  std::string string_current_test_case_info() const {
    auto const prefix = fmt::format("[{}]: ", current_test_case_info().number);
    auto const prefix_pad = _prefix_cols - prefix.size() - 1;

    std::string const trunc_name(
        current_test_case_info().name.begin(),
        current_test_case_info().name.begin()
            + std::min(prefix_pad, current_test_case_info().name.size()));
    // This is the prefix of length (_line_cols - _time_cols)
    return fmt::format("{}{:<{}} ", prefix, trunc_name, prefix_pad);
  }

  void testCaseStarting(Catch::TestCaseInfo const& testInfo) override {
    _current_section_name = testInfo.name;
    set_current_test_case_info(testInfo);
    if (current_test_case_info().category != "extreme") {
      fmt::print("{}", string_current_test_case_info());
    } else {
      print_extreme_test_divider("START");
    }
  }

  void sectionStarting(Catch::SectionInfo const& sectionInfo) override {
    // TODO handle arbitrary depth subsections
    set_most_recent_section_info(sectionInfo);
    if (current_test_case_info().category != "extreme") {
      if (_section_depth == 1 && _current_section_name != sectionInfo.name) {
        _current_section_name             = sectionInfo.name;
        constexpr std::string_view prefix = "\n-- with ";
        auto const prefix_pad             = _prefix_cols - prefix.size() + 1;
        fmt::print("{}{: <{}}", prefix, sectionInfo.name, prefix_pad);
      } else if (_section_depth == 2) {
        constexpr std::string_view prefix = "\n---- ";
        auto const prefix_pad             = _prefix_cols - prefix.size() + 1;
        fmt::print("{}{: <{}}", prefix, sectionInfo.name, prefix_pad);
      }
    } else {
      if (_section_depth > 0) {
        print_extreme_test_divider(sectionInfo.name + " - START");
      }
    }
    _section_depth++;
  }

  void sectionEnded(Catch::SectionStats const& sectionStats) override {
    using libsemigroups::detail::string_time;

    _section_depth--;

    if (_section_depth == 0) {
      _section_stats.emplace_back(sectionStats, most_recent_section_info());
      _total_time += _section_stats.back().duration;
      std::string section_duration
          = string_time(_section_stats.back().duration);
      if (current_test_case_info().category != "extreme") {
        fmt::print("{:.>{}}", section_duration, _time_cols);
      } else {
        if (most_recent_section_info().name != current_test_case_info().name) {
          // In this case the leaf section that was run was a proper
          // subsection, and not just the entire test case, so we print the end
          // of that proper subsection.
          print_extreme_test_divider(most_recent_section_info().name + " - "
                                     + section_duration + " - STOP");
        }
      }
    }
  }

  void testCaseEnded(Catch::TestCaseStats const&) override {
    using libsemigroups::detail::string_time;

    if (current_test_case_info().category == "extreme") {
      auto section_duration = string_time(_total_time);

      print_extreme_test_divider(section_duration + " - STOP");
      if (_section_stats.size() > 1) {
        fmt::print(_extreme_emph, "{:=>{}}\n", "", _line_cols);
        fmt::print(
            _extreme_emph, "Summary for {}", string_current_test_case_info());
        constexpr std::string_view prefix = "\n-- with ";
        for (SectionStats const& ss : _section_stats) {
          auto const prefix_pad = _prefix_cols - prefix.size() + 1;
          fmt::print(_extreme_emph, "{}{: <{}}", prefix, ss.name, prefix_pad);
          fmt::print(
              _extreme_emph, "{:.>{}}", string_time(ss.duration), _time_cols);
        }
        fmt::print(_extreme_emph, "\n{:=>{}}", "", _line_cols);
        fmt::print("\n");
      }
    } else {
      fmt::print("\n");
    }
    _section_stats.clear();
  }

  void testRunEnded(Catch::TestRunStats const&) override {
    using libsemigroups::detail::string_time;

    constexpr std::string_view prefix     = "Total time ";
    auto const                 t          = string_time(_total_time);
    auto const                 prefix_pad = _line_cols - prefix.size();

    fmt::print(_run_end_emph, "{:=>{}}\n", "", _line_cols);
    fmt::print(_run_end_emph, "{}{:.>{}}\n", prefix, t, prefix_pad);
    // The following =s fill in the line printed by catch to make it the same
    // width as _line_cols
    fmt::print(_run_end_emph, "{:=>{}}", "", _line_cols - 79);
  }

  std::chrono::nanoseconds _total_time = std::chrono::nanoseconds(0);
  size_t                   _section_depth;
  std::string              _current_section_name;

  TestCaseInfo              _current_test_case_info;
  SectionInfo               _most_recent_section_info;
  std::vector<SectionStats> _section_stats;

  static constexpr size_t _line_cols   = 90;
  static constexpr size_t _time_cols   = 12;
  static constexpr size_t _prefix_cols = _line_cols - _time_cols;

  static constexpr auto _extreme_emph = fmt::emphasis::bold;
  static constexpr auto _run_end_emph
      = fmt::emphasis::bold | fg(fmt::terminal_color::bright_green);
};

CATCH_REGISTER_LISTENER(LibsemigroupsListener)
