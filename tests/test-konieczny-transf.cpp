// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2020 Finn Smith
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
// TODO(later) replace Transf(...) with Transf::make

#include <cstddef>  // for size_t

#include "catch_amalgamated.hpp"  // for REQUIRE
#include "test-main.hpp"          // FOR LIBSEMIGROUPS_TEST_CASE

#include "libsemigroups/konieczny.hpp"  // for Konieczny
#include "libsemigroups/transf.hpp"     // for Transf<>

namespace libsemigroups {

  constexpr bool REPORT = false;

  LIBSEMIGROUPS_TEST_CASE("Konieczny",
                          "031",
                          "transformations",
                          "[quick][transf]") {
    using Transf         = LeastTransf<5>;
    auto              rg = ReportGuard(REPORT);
    Konieczny<Transf> S({Transf({1, 0, 2, 3, 4}),
                         Transf({1, 2, 3, 4, 0}),
                         Transf({0, 0, 2, 3, 4})});
    S.run();
    REQUIRE(S.size() == 3125);

    size_t sum = 0;
    std::for_each(
        S.cbegin_D_classes(),
        S.cend_D_classes(),
        [&sum, &S](Konieczny<Transf>::DClass const& x) {
          sum += S.D_class_of_element(x.rep()).number_of_idempotents();
        });
    REQUIRE(sum == 196);
    REQUIRE(S.number_of_idempotents() == 196);
    REQUIRE(std::vector<Transf>(S.cbegin_generators(), S.cend_generators())
            == std::vector<Transf>({Transf({1, 0, 2, 3, 4}),
                                    Transf({1, 2, 3, 4, 0}),
                                    Transf({0, 0, 2, 3, 4})}));
  }

  LIBSEMIGROUPS_TEST_CASE("Konieczny",
                          "032",
                          "transformations - JDM favourite example",
                          "[quick][no-valgrind][transf]") {
    using Transf         = LeastTransf<8>;
    auto              rg = ReportGuard(REPORT);
    Konieczny<Transf> S({Transf({1, 7, 2, 6, 0, 4, 1, 5}),
                         Transf({2, 4, 6, 1, 4, 5, 2, 7}),
                         Transf({3, 0, 7, 2, 4, 6, 2, 4}),
                         Transf({3, 2, 3, 4, 5, 3, 0, 1}),
                         Transf({4, 3, 7, 7, 4, 5, 0, 4}),
                         Transf({5, 6, 3, 0, 3, 0, 5, 1}),
                         Transf({6, 0, 1, 1, 1, 6, 3, 4}),
                         Transf({7, 7, 4, 0, 6, 4, 1, 7})});
    REQUIRE(S.size() == 597369);
    size_t sum = 0;
    std::for_each(
        S.cbegin_D_classes(),
        S.cend_D_classes(),
        [&sum, &S](Konieczny<Transf>::DClass const& x) {
          sum += S.D_class_of_element(x.rep()).number_of_idempotents();
        });
    REQUIRE(sum == 8194);
    REQUIRE(S.number_of_idempotents() == 8194);
  }

  LIBSEMIGROUPS_TEST_CASE("Konieczny",
                          "033",
                          "transformations - large example",
                          "[quick][no-valgrind][transf]") {
    auto                        rg   = ReportGuard(REPORT);
    std::vector<Transf<>> const gens = {Transf<>({2, 1, 0, 4, 2, 1, 1, 8, 0}),
                                        Transf<>({1, 7, 6, 2, 5, 1, 1, 4, 3}),
                                        Transf<>({1, 0, 7, 2, 1, 3, 1, 3, 7}),
                                        Transf<>({0, 3, 8, 1, 2, 8, 1, 7, 0}),
                                        Transf<>({0, 0, 0, 2, 7, 7, 5, 5, 3})};
    Konieczny<Transf<>>         S(gens);

    for (auto x : gens) {
      REQUIRE(S.contains(x));
    }

    REQUIRE(S.current_size() < 15000);
    REQUIRE(S.current_number_of_regular_elements() < 10000);
    REQUIRE(S.current_number_of_idempotents() < 500);
    REQUIRE(S.current_number_of_D_classes() < 2000);
    REQUIRE(S.current_number_of_L_classes() < 4000);
    REQUIRE(S.current_number_of_R_classes() < 6500);

    REQUIRE(S.size() == 232511);
    REQUIRE(S.current_number_of_D_classes() == 2122);
    REQUIRE(S.current_number_of_L_classes() == 8450);
    REQUIRE(S.current_number_of_R_classes() == 14706);
  }

  LIBSEMIGROUPS_TEST_CASE("Konieczny",
                          "034",
                          "transformations - large example with stop",
                          "[quick][no-valgrind][transf]") {
    auto                rg = ReportGuard(REPORT);
    Konieczny<Transf<>> S({Transf<>({2, 1, 0, 4, 2, 1, 1, 8, 0}),
                           Transf<>({1, 7, 6, 2, 5, 1, 1, 4, 3}),
                           Transf<>({1, 0, 7, 2, 1, 3, 1, 3, 7}),
                           Transf<>({0, 3, 8, 1, 2, 8, 1, 7, 0}),
                           Transf<>({0, 0, 0, 2, 7, 7, 5, 5, 3})});
    S.run_for(std::chrono::milliseconds(100));
    REQUIRE(S.size() == 232511);
  }

  LIBSEMIGROUPS_TEST_CASE("Konieczny",
                          "035",
                          "transformations - large example with run_until",
                          "[quick][no-valgrind][transf]") {
    auto                rg = ReportGuard(REPORT);
    Konieczny<Transf<>> S({Transf<>({2, 1, 0, 4, 2, 1, 1, 8, 0}),
                           Transf<>({1, 7, 6, 2, 5, 1, 1, 4, 3}),
                           Transf<>({1, 0, 7, 2, 1, 3, 1, 3, 7}),
                           Transf<>({0, 3, 8, 1, 2, 8, 1, 7, 0}),
                           Transf<>({0, 0, 0, 2, 7, 7, 5, 5, 3})});
    S.run_until([&S]() -> bool {
      return S.cend_D_classes() - S.cbegin_D_classes() > 20;
    });

    size_t number_of_classes1 = S.cend_D_classes() - S.cbegin_D_classes();
    REQUIRE(number_of_classes1 >= 20);
    S.run();
    size_t number_of_classes2 = S.cend_D_classes() - S.cbegin_D_classes();
    REQUIRE(S.size() == 232511);
    REQUIRE(number_of_classes1 < number_of_classes2);
    REQUIRE(number_of_classes2 == 2122);
  }

  LIBSEMIGROUPS_TEST_CASE("Konieczny",
                          "036",
                          "transformations - large example with stop in Action",
                          "[quick][no-valgrind][transf]") {
    auto                rg = ReportGuard(REPORT);
    Konieczny<Transf<>> S({Transf<>({2, 1, 0, 4, 2, 1, 1, 8, 0}),
                           Transf<>({1, 7, 6, 2, 5, 1, 1, 4, 3}),
                           Transf<>({1, 0, 7, 2, 1, 3, 1, 3, 7}),
                           Transf<>({0, 3, 8, 1, 2, 8, 1, 7, 0}),
                           Transf<>({0, 0, 0, 2, 7, 7, 5, 5, 3})});
    S.run_for(std::chrono::milliseconds(5));
    S.run_for(std::chrono::milliseconds(5));
    S.run_for(std::chrono::milliseconds(5));
    S.run_for(std::chrono::milliseconds(100));
    S.run_for(std::chrono::milliseconds(100));
    S.run();
    S.run_for(std::chrono::milliseconds(100));
    S.run_for(std::chrono::milliseconds(100));
    REQUIRE(S.size() == 232511);
  }

  LIBSEMIGROUPS_TEST_CASE("Konieczny", "037", "exceptions", "[quick][transf]") {
    auto rg          = ReportGuard(REPORT);
    using point_type = typename Transf<>::point_type;
    std::vector<point_type> v(65, 0);
    std::iota(v.begin(), v.end(), 0);
    REQUIRE_THROWS_AS(Konieczny<Transf<>>({Transf<>(v)}),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("Konieczny",
                          "038",
                          "transformations: contains",
                          "[quick][transf]") {
    auto                rg = ReportGuard(REPORT);
    Konieczny<Transf<>> S({Transf<>({1, 0, 2, 3, 4}),
                           Transf<>({1, 2, 3, 4, 0}),
                           Transf<>({0, 0, 2, 3, 4})});
    REQUIRE(S.contains(Transf<>({1, 0, 2, 3, 4})));
    REQUIRE(S.contains(Transf<>({1, 2, 3, 4, 0})));
    REQUIRE(S.contains(Transf<>({0, 0, 2, 3, 4})));
    REQUIRE(!S.contains(Transf<>({1, 0, 2, 3, 4, 5})));
    REQUIRE(!S.contains(Transf<>({1, 2, 3, 4, 0, 5})));
    REQUIRE(!S.contains(Transf<>({0, 0, 2, 3, 4, 1})));

    REQUIRE_THROWS_AS(S.D_class_of_element(Transf<>({1, 0, 2, 3, 4, 5})),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(S.D_class_of_element(Transf<>({1, 2, 3, 4, 0, 5})),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(S.D_class_of_element(Transf<>({0, 0, 2, 3, 4, 1})),
                      LibsemigroupsException);

    Konieczny<Transf<>> T(
        {Transf<>({1, 0, 3, 4, 2}), Transf<>({0, 0, 2, 3, 4})});
    REQUIRE(T.contains(Transf<>({1, 0, 2, 3, 4})));
    REQUIRE(T.contains(Transf<>({0, 0, 2, 3, 4})));
    REQUIRE(!T.contains(Transf<>({1, 2, 3, 4, 0})));
    REQUIRE(!T.contains(Transf<>({1, 2, 3, 0, 4})));
    REQUIRE(!T.contains(Transf<>({1, 2, 3, 4, 0, 5})));
    REQUIRE(!T.contains(Transf<>({0, 2, 3, 4, 1})));

    REQUIRE_THROWS_AS(T.D_class_of_element(Transf<>({1, 2, 3, 4, 0})),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(T.D_class_of_element(Transf<>({1, 2, 3, 4, 0, 5})),
                      LibsemigroupsException);
    REQUIRE_THROWS_AS(T.D_class_of_element(Transf<>({0, 2, 3, 4, 1})),
                      LibsemigroupsException);
  }

  LIBSEMIGROUPS_TEST_CASE("Konieczny",
                          "039",
                          "transformations Hall monoid 5",
                          "[extreme][transf]") {
    auto rg      = ReportGuard();
    using Transf = LeastTransf<31>;
    Konieczny<Transf> K;
    K.add_generator(
        Transf({0, 16, 1, 17, 2,  18, 3,  19, 4,  20, 5,  21, 6,  22, 7, 23,
                8, 24, 9, 25, 10, 26, 11, 27, 12, 28, 13, 29, 14, 30, 15}));
    K.add_generator(
        Transf({0, 1, 2,  3,  4,  5,  6,  7,  16, 17, 18, 19, 20, 21, 22, 23,
                8, 9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30}));
    K.add_generator(
        Transf({0, 16, 8, 24, 4, 20, 12, 28, 2, 18, 10, 26, 6, 22, 14, 30,
                0, 17, 8, 25, 4, 21, 12, 29, 2, 19, 10, 27, 6, 23, 14}));
    K.add_generator(Transf({0, 0, 0, 16, 0, 8, 4, 28, 2, 2, 2, 18, 2, 10, 6, 30,
                            1, 1, 1, 17, 1, 9, 5, 29, 3, 3, 3, 19, 3, 11, 7}));
    K.add_generator(Transf({0, 0, 0, 16, 0, 8, 0, 24, 0, 4, 0, 20, 0, 12, 2, 30,
                            1, 1, 1, 17, 1, 9, 1, 25, 1, 5, 1, 21, 1, 13, 3}));
    K.add_generator(Transf({0, 0, 0, 16, 0, 8, 0, 24, 0, 0, 4, 20, 2, 10, 6, 30,
                            1, 1, 1, 17, 1, 9, 1, 25, 1, 1, 5, 21, 3, 11, 7}));
    K.add_generator(
        Transf({0, 0, 0, 16, 0, 8,  0, 24, 0, 4, 0, 20, 0, 12, 0, 28,
                0, 2, 0, 18, 0, 10, 0, 26, 0, 6, 0, 22, 0, 14, 1}));
    K.add_generator(Transf({0, 0, 0, 16, 0, 8, 0, 24, 0, 4, 0, 20, 0, 12, 0, 28,
                            0, 0, 2, 18, 0, 8, 2, 26, 0, 4, 2, 22, 1, 13, 3}));
    K.add_generator(
        Transf({0, 0, 0, 16, 0, 8,  0, 24, 0, 0, 4, 20, 0, 8,  4, 28,
                0, 0, 0, 16, 2, 10, 2, 26, 0, 1, 4, 21, 2, 11, 6}));
    K.add_generator(Transf({0, 0, 0, 16, 0, 8, 0, 24, 0, 0, 4, 20, 0, 8, 4, 28,
                            0, 0, 0, 16, 0, 8, 2, 26, 0, 1, 4, 21, 0, 9, 6}));
    K.add_generator(Transf({0, 0, 0, 16, 0, 8, 0, 24, 0, 0, 0, 16, 0, 8, 4, 28,
                            0, 0, 0, 16, 0, 8, 2, 26, 0, 1, 0, 17, 0, 9, 6}));
    K.add_generator(
        Transf({0, 0, 0, 16, 0, 8,  0, 24, 0, 0, 4, 20, 0, 8,  4, 28,
                0, 0, 0, 16, 2, 10, 2, 26, 1, 1, 5, 21, 3, 11, 7}));
    REQUIRE(K.size() == 23191071);
  }

  LIBSEMIGROUPS_TEST_CASE("Konieczny",
                          "040",
                          "transformations - destructor coverage",
                          "[quick][transf][no-valgrind]") {
    auto rg      = ReportGuard(REPORT);
    using Transf = LeastTransf<9>;
    Konieczny<Transf> S({Transf({2, 1, 0, 4, 2, 1, 1, 8, 0}),
                         Transf({1, 7, 6, 2, 5, 1, 1, 4, 3}),
                         Transf({1, 0, 7, 2, 1, 3, 1, 3, 7}),
                         Transf({0, 3, 8, 1, 2, 8, 1, 7, 0}),
                         Transf({0, 0, 0, 2, 7, 7, 5, 5, 3})});
    S.run_until(
        [&S]() -> bool { return S.current_number_of_regular_D_classes() > 2; });

    // if these fail, this test won't get the coverage hoped for
    REQUIRE(S.current_number_of_regular_D_classes() < 5);
    REQUIRE(S.current_number_of_D_classes() - S.number_of_regular_D_classes()
            < 2117);
    // now all of the destructor should run
  }

  LIBSEMIGROUPS_TEST_CASE("Konieczny",
                          "041",
                          "current_number_D_classes",
                          "[quick][transf][no-valgrind]") {
    auto rg      = ReportGuard(REPORT);
    using Transf = LeastTransf<9>;
    Konieczny<Transf> S({Transf({2, 1, 0, 4, 2, 1, 1, 8, 0}),
                         Transf({1, 7, 6, 2, 5, 1, 1, 4, 3}),
                         Transf({1, 0, 7, 2, 1, 3, 1, 3, 7}),
                         Transf({0, 3, 8, 1, 2, 8, 1, 7, 0}),
                         Transf({0, 0, 0, 2, 7, 7, 5, 5, 3})});
    REQUIRE(S.current_number_of_regular_D_classes() == 0);
    REQUIRE(S.current_number_of_D_classes() == 0);
    S.run_until(
        [&S]() -> bool { return S.current_number_of_regular_D_classes() > 2; });

    S.run();
    REQUIRE(S.current_number_of_regular_D_classes() == 5);
    REQUIRE(S.current_number_of_D_classes() - S.number_of_regular_D_classes()
            == 2117);
  }

  LIBSEMIGROUPS_TEST_CASE("Konieczny",
                          "042",
                          "add_generator",
                          "[quick][transf][no-valgrind]") {
    auto rg      = ReportGuard(REPORT);
    using Transf = LeastTransf<5>;

    Konieczny<Transf> S({Transf({1, 0, 2, 3, 4})});
    S.add_generator(Transf({1, 2, 3, 4, 0}));
    S.add_generator(Transf({0, 0, 2, 3, 4}));

    REQUIRE(S.generator(0) == Transf({1, 0, 2, 3, 4}));
    REQUIRE(S.generator(1) == Transf({1, 2, 3, 4, 0}));
    REQUIRE(S.generator(2) == Transf({0, 0, 2, 3, 4}));
    REQUIRE(S.number_of_generators() == 3);
    REQUIRE(S.degree() == Degree<Transf>()(Transf({1, 0, 2, 3, 4})));
  }

  LIBSEMIGROUPS_TEST_CASE("Konieczny",
                          "043",
                          "add_generator, init",
                          "[quick][transf][no-valgrind]") {
    auto rg      = ReportGuard(REPORT);
    using Transf = LeastTransf<5>;

    Konieczny<Transf> S;
    S.add_generator(Transf({1, 2, 3, 4, 0}));
    S.add_generator(Transf({0, 0, 2, 3, 4}));

    REQUIRE(S.degree() == Degree<Transf>()(Transf({1, 2, 3, 4, 0})));
    REQUIRE(S.number_of_generators() == 2);

    S.run();
    REQUIRE(S.size() == 610);

    S.init();
    S.add_generator(Transf({1, 2, 3, 4, 0}));
    S.run();
    REQUIRE(S.size() == 5);

    S.init();
    S.add_generator(Transf({1, 2, 3, 4, 0}));
    S.add_generator(Transf({0, 0, 2, 3, 4}));
    REQUIRE(S.size() == 610);
  }

}  // namespace libsemigroups
