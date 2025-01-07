//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2023 James D. Mitchell
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

// This file includes relevant headers from fmtlib and suppresses the compiler
// warnings emitted from it.

#ifndef LIBSEMIGROUPS_DETAIL_FMT_HPP_
#define LIBSEMIGROUPS_DETAIL_FMT_HPP_
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wall"
#pragma GCC diagnostic ignored "-Wswitch-default"
#pragma GCC diagnostic ignored "-Wswitch-enum"
#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Wtautological-compare"

#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/compile.h>
#include <fmt/format.h>
#include <fmt/printf.h>
#include <fmt/ranges.h>

#pragma GCC diagnostic pop
#endif  // LIBSEMIGROUPS_DETAIL_FMT_HPP_
