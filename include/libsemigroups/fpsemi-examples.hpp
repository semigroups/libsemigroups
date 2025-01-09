//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2022 Murray T. Whyte
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
// This file contains declarations of functions for constructing various
// standard examples of presentations.

#ifndef LIBSEMIGROUPS_FPSEMI_EXAMPLES_HPP_
#define LIBSEMIGROUPS_FPSEMI_EXAMPLES_HPP_

#include <cstddef>  // for size_t
#include <cstdint>  // for uint64_t
#include <vector>   // for vector

#include "presentation.hpp"  // for Presentation
#include "types.hpp"         // for word_type

namespace libsemigroups {

  //! \ingroup presentations_group
  //!
  //! \brief Namespace for presentations of some finitely presented semigroups.
  //!
  //! Defined in `fpsemi-examples.hpp`.
  //!
  //! This namespace contains functions which give presentations for some common
  //! finitely presented semigroups and monoids.
  namespace fpsemigroup {
    //! \brief Enum class containing various authors.
    //!
    //! The values in this enum class are used to specify the authors of a
    //! presentation. Where there are different presentations by different
    //! authors, values of this type can be passed as an argument to
    //! disambiguate which presentation is wanted.
    enum class author : uint64_t {

      //! Encoded as the `uint64_t` value `0`
      Any = 0,

      //! Encoded as the `uint64_t` value `1`
      Machine = 1,

      //! Encoded as the `uint64_t` value `2`
      Aizenstat = 2,

      //! Encoded as the `uint64_t` value `4`
      Burnside = 4,

      //! Encoded as the `uint64_t` value `8`
      Carmichael = 8,

      //! Encoded as the `uint64_t` value `16`
      Coxeter = 16,

      //! Encoded as the `uint64_t` value `32`
      Easdown = 32,

      //! Encoded as the `uint64_t` value `64`
      East = 64,

      //! Encoded as the `uint64_t` value `128`
      Fernandes = 128,

      //! Encoded as the `uint64_t` value `256`
      FitzGerald = 256,

      //! Encoded as the `uint64_t` value `512`
      Gay = 512,

      //! Encoded as the `uint64_t` value `1024`
      Godelle = 1024,

      //! Encoded as the `uint64_t` value `2048`
      Guralnick = 2048,

      //! Encoded as the `uint64_t` value `4096`
      Iwahori = 4096,

      //! Encoded as the `uint64_t` value `8192`
      Kantor = 8192,

      //! Encoded as the `uint64_t` value `16'384`
      Kassabov = 16'384,

      //! Encoded as the `uint64_t` value `32'768`
      Lubotzky = 32'768,

      //! Encoded as the `uint64_t` value `65'536`
      Miller = 65'536,

      //! Encoded as the `uint64_t` value `131'072`
      Mitchell = 131'072,

      //! Encoded as the `uint64_t` value `262'144`
      Moore = 262'144,

      //! Encoded as the `uint64_t` value `524'288`
      Moser = 524'288,

      //! Encoded as the `uint64_t` value `1'048'576`
      Sutov = 1'048'576,

      //! Encoded as the `uint64_t` value `2'097'152`
      Whyte = 2'097'152
    };

    //! \brief Operator to combine \ref author values.
    //!
    //! This operator can be used arbitrarily to combine \ref author values.
    [[nodiscard]] inline author operator+(author auth1, author auth2) {
      return static_cast<author>(static_cast<uint64_t>(auth1)
                                 + static_cast<uint64_t>(auth2));
    }

    //! \brief A presentation for the stellar monoid.
    //!
    //! Returns the stellar monoid with \p l generators, as in Theorem 4.39 of
    //! [10.48550/arXiv.1910.11740][].
    //!
    //! \param l the number of generators.
    //! \param val the author (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `l < 2`.
    //! \throws LibsemigroupsException if \p val is not author::Any.
    //!
    //! [10.48550/arXiv.1910.11740]: https://doi.org/10.48550/arXiv.1910.11740
    [[nodiscard]] Presentation<word_type> stellar_monoid(size_t l,
                                                         author val
                                                         = author::Any);

    //! \brief A presentation for the dual symmetric inverse monoid.
    //!
    //! Returns a monoid presentation defining the dual symmetric inverse monoid
    //! of degree \p n. The argument \p val determines the specific presentation
    //! which is returned. The only option is:
    //! * `author::Easdown + author::East + author::FitzGerald` (from Section 3
    //! of [10.48550/arxiv.0707.2439][]).
    //!
    //! \param n the degree.
    //! \param val the author of the presentation (default: `author::Easdown +
    //! author::East + author::FitzGerald`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 3`.
    //! \throws LibsemigroupsException if \p val is not
    //! `author::Easdown + author::East + author::FitzGerald`.
    //!
    //! [10.48550/arxiv.0707.2439]: https://doi.org/10.48550/arxiv.0707.2439
    [[nodiscard]] Presentation<word_type> dual_symmetric_inverse_monoid(
        size_t n,
        author val = author::Easdown + author::East + author::FitzGerald);

    //! \brief A presentation for the uniform block bijection monoid.
    //!
    //! Returns a monoid presentation defining the uniform block bijection
    //! monoid of degree \p n. The argument \p val determines the specific
    //! presentation which is returned. The only option is:
    //! * `author::FitzGerald` (see [10.1017/s0004972700037692][]).
    //!
    //! \param n the degree.
    //! \param val the author of the presentation (default:
    //! `author::FitzGerald`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 3`.
    //! \throws LibsemigroupsException if \p val is not `author::FitzGerald`.
    //!
    //! [10.1017/s0004972700037692]: https://doi.org/10.1017/s0004972700037692
    [[nodiscard]] Presentation<word_type>
    uniform_block_bijection_monoid(size_t n, author val = author::FitzGerald);

    //! \brief A presentation for the partition monoid.
    //!
    //! Returns a presentation defining the partition monoid of degree \p n. The
    //! argument \p val determines the specific presentation which is returned.
    //! The options are:
    //! * `author::Machine`;
    //! * `author::East` (see Theorem 41 of
    //! [10.1016/j.jalgebra.2011.04.008][]).
    //!
    //! Note that `author::East` returns a monoid presentation, and
    //! `author::Machine` returns a semigroup presentation.
    //!
    //! \param n the degree.
    //! \param val the author of the presentation (default: `author::East`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `val == author::Machine && n != 3`.
    //! \throws LibsemigroupsException if `val == author::East && n < 4`.
    //! \throws LibsemigroupsException if \p val is not one of the options
    //! specified above.
    //!
    //! [10.1016/j.jalgebra.2011.04.008]:
    //! https://doi.org/10.1016/j.jalgebra.2011.04.008
    [[nodiscard]] Presentation<word_type> partition_monoid(size_t n,
                                                           author val
                                                           = author::East);

    //! \brief A presentation for the singular part of the Brauer monoid.
    //!
    //! Returns a monoid presentation for the singular part of the Brauer monoid
    //! of degree \p n, as in Theorem 5 of [10.21136/MB.2007.134125][]).
    //!
    //! \param n the degree.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 3`.
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    //!
    //! [10.21136/MB.2007.134125]: https://doi.org/10.21136/MB.2007.134125
    [[nodiscard]] Presentation<word_type> singular_brauer_monoid(size_t n,
                                                                 author val
                                                                 = author::Any);

    //! \brief A presentation for the monoid of orientation preserving
    //! mappings.
    //!
    //! Returns a monoid presentation defining the monoid of orientation
    //! preserving mappings on a finite chain of order \p n, as described in
    //! [10.1007/s10012-000-0001-1][].
    //!
    //! \param n the order of the chain.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 3`.
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    //!
    //! [10.1007/s10012-000-0001-1]: https://doi.org/10.1007/s10012-000-0001-1
    [[nodiscard]] Presentation<word_type>
    orientation_preserving_monoid(size_t n, author val = author::Any);

    //! \brief A presentation for the monoid of orientation preserving or
    //! reversing mappings.
    //!
    //! Returns a monoid presentation defining the monoid of orientation
    //! preserving or reversing mappings on a finite chain of order \p n, as
    //! described in [10.1007/s10012-000-0001-1][].
    //!
    //! \param n the order of the chain.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 3`.
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    //!
    //! [10.1007/s10012-000-0001-1]: https://doi.org/10.1007/s10012-000-0001-1
    [[nodiscard]] Presentation<word_type>
    orientation_preserving_reversing_monoid(size_t n, author val = author::Any);

    //! \brief A presentation for the Temperley-Lieb monoid.
    //!
    //! Returns a monoid presentation defining the Temperley-Lieb monoid with
    //! \p n generators, as described in Theorem 2.2 of
    //! [10.1093/qmath/haab001][].
    //!
    //! \param n the number of generators.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 3`.
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    //!
    //! [10.1093/qmath/haab001]: https://doi.org/10.1093/qmath/haab001
    [[nodiscard]] Presentation<word_type> temperley_lieb_monoid(size_t n,
                                                                author val
                                                                = author::Any);

    //! \brief A presentation for the Brauer monoid.
    //!
    //! Returns a monoid presentation defining the Brauer monoid of degree \p n,
    //! as described in Theorem 3.1 of the paper [10.2478/s11533-006-0017-6][].
    //!
    //! \param n the degree.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    //!
    //! [10.2478/s11533-006-0017-6]:
    //! https://doi.org/10.2478/s11533-006-0017-6
    [[nodiscard]] Presentation<word_type> brauer_monoid(size_t n,
                                                        author val
                                                        = author::Any);

    // TODO(2): Does this have a reference?
    //! \brief A presentation for the partial Brauer monoid.
    //!
    //! Returns a monoid presentation defining the partial Brauer monoid of
    //! degree \p n..
    //!
    //! \param n the degree.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    [[nodiscard]] Presentation<word_type> partial_brauer_monoid(size_t n,
                                                                author val
                                                                = author::Any);

    //! \brief A presentation for a Fibonacci semigroup.
    //!
    //! Returns a semigroup presentation defining the Fibonacci semigroup
    //! \f$F(r, n)\f$, where \f$r\f$ is \p r and \f$n\f$ is \p n, as described
    //! in [10.1016/0022-4049(94)90005-1][].
    //!
    //! \param r the length of the left hand sides of the relations.
    //! \param n the number of generators.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n == 0`.
    //! \throws LibsemigroupsException if `r == 0`.
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    //!
    //! [10.1016/0022-4049(94)90005-1]:
    //! https://doi.org/10.1016/0022-4049(94)90005-1
    [[nodiscard]] Presentation<word_type>
    fibonacci_semigroup(size_t r, size_t n, author val = author::Any);

    //! \brief A presentation for the plactic monoid.
    //!
    //! Returns a monoid presentation defining the plactic monoid with \p n
    //! generators (see Section 3 of [10.1007/s00233-022-10285-3][]).
    //!
    //! \param n the number of generators.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 2`.
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    //!
    //! [10.1007/s00233-022-10285-3]: https://doi.org/10.1007/s00233-022-10285-3
    [[nodiscard]] Presentation<word_type> plactic_monoid(size_t n,
                                                         author val
                                                         = author::Any);

    //! \brief A presentation for the stylic monoid.
    //!
    //! Returns a monoid presentation defining the stylic monoid with \p n
    //! generators (see Theorem 8.1 of [10.1007/s00233-022-10285-3]).
    //!
    //! \param n the number of generators.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 2`.
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    //!
    //! [10.1007/s00233-022-10285-3]: https://doi.org/10.1007/s00233-022-10285-3
    [[nodiscard]] Presentation<word_type> stylic_monoid(size_t n,
                                                        author val
                                                        = author::Any);

    //! \brief A presentation for the symmetric group.
    //!
    //! Returns a monoid presentation for the symmetric group of degree \p n.
    //! The arguments \p val and \p index determine the specific presentation
    //! which is returned. The options are:
    //!
    // clang-format off
    //!
    //! Author                              | Index | No. generators  | No. relations             | Reference
    //! ----------------------------------- | ----- | --------------- | ------------------------- | ---------------------------------------------------------------
    //! `author::Burnside + author::Miller` | `0`   | \f$n - 1\f$     | \f$n^3 - 5n^2 + 9n - 5\f$ | p.464 of [10.1017/CBO9781139237253][]   <!-- NOLINT -->
    //! `author::Carmichael`                | `0`   | \f$n - 1\f$     | \f$(n - 1)^2\f$           | Comment 9.5.2 of [10.1007/978-1-84800-281-4][]   <!-- NOLINT -->
    //! `author::Coxeter + author::Moser`   | `0`   | \f$n - 1\f$     | \f$n(n + 1)/2\f$          | Ch.3, Prop 1.2 of [hdl.handle.net/10023/2821][]  <!-- NOLINT -->
    //! `author::Moore`                     | `0`   | \f$2\f$         | \f$n + 1\f$               | Ch. 3, Prop 1.1 of [hdl.handle.net/10023/2821][] <!-- NOLINT -->
    //! ^                                   | `1`   | \f$n - 1\f$     | \f$n(n + 1)/2\f$          | Comment 9.5.3 of [10.1007/978-1-84800-281-4][]   <!-- NOLINT -->
    //!
    //! [10.1017/CBO9781139237253]: https://doi.org/10.1017/CBO9781139237253
    //! [10.1007/978-1-84800-281-4]: https://doi.org/10.1007/978-1-84800-281-4
    //! [hdl.handle.net/10023/2821]: http://hdl.handle.net/10023/2821
    //!
    // clang-format on
    //!
    //! The default for `val` is `author::Carmichael`. The default for `index`
    //! is `0`.
    //!
    //! \param n the degree of the symmetric group.
    //! \param val the author of the presentation (default:
    //! `author::Carmichael`).
    //! \param index which of the author's presentations to return, if multiple
    //! exist (default: `0`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 4`.
    //! \throws LibsemigroupsException if the author-index combination is
    //! not present in the table above.
    //!
    [[nodiscard]] Presentation<word_type> symmetric_group(size_t n,
                                                          author val
                                                          = author::Carmichael,
                                                          size_t index = 0);

    //! \brief A presentation for the alternating group.
    //!
    //! Returns a monoid presentation defining the alternating group of degree
    //! \p n. The argument \p val determines the specific presentation which is
    //! returned. The options are:
    //! * `author::Moore` (see Ch. 3, Prop 1.3 of
    //! [hdl.handle.net/10023/2821][]).
    //!
    //! \param n the degree of the alternating group.
    //! \param val the author of the presentation (default: `author::Moore`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 4`.
    //! \throws LibsemigroupsException if `val` is not `author::Moore`.
    //!
    //! [hdl.handle.net/10023/2821]: http://hdl.handle.net/10023/2821
    [[nodiscard]] Presentation<word_type> alternating_group(size_t n,
                                                            author val
                                                            = author::Moore);

    //! \brief A presentation for a rectangular band.
    //!
    //! Returns a semigroup presentation defining the \p m by \p n rectangular
    //! band, as given in Proposition 4.2 of [10.1007/s002339910016][].
    //!
    //! \param m the number of rows.
    //! \param n the number of columns.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `m == 0`.
    //! \throws LibsemigroupsException if `n == 0`.
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    //!
    //! [10.1007/s002339910016]: https://doi.org/10.1007/s002339910016
    [[nodiscard]] Presentation<word_type>
    rectangular_band(size_t m, size_t n, author val = author::Any);

    //! \brief A presentation for the full transformation monoid.
    //!
    //! Returns a monoid presentation defining the full transformation monoid.
    //! The arguments \p val and \p index determine the specific presentation
    //! which is returned. The options are:
    //! * `author::Aizenstat` and `index = 0` (see Ch. 3, Prop 1.7 of
    //! [http://hdl.handle.net/10023/2821][]);
    //! * `author::Iwahori` and `index = 0` (see Theorem 9.3.1 of
    //! [10.1007/978-1-84800-281-4][]);
    //! * `author::Mitchell + author::Whyte` (see Theorem 1.5 of
    //! [10.48550/arXiv.2406.19294][]):
    //!    * `index = 0` for the presentation with five non-symmetric-group
    //!    relations;
    //!    * `index = 1` for the presentation with four non-symmetric-group
    //!    relations, valid for odd degree.
    //!
    //! \param n the degree of the full transformation monoid.
    //! \param val the author of the presentation (default: `author::Mitchell +
    //! author::Whyte`).
    //! \param index which of the author's presentations to
    //! return, if multiple exist (default: `0`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 4`.
    //! \throws LibsemigroupsException if the author-index combination is
    //! not present in the table above.
    //!
    //! [http://hdl.handle.net/10023/2821]: http://hdl.handle.net/10023/2821
    //! [10.1007/978-1-84800-281-4]: https://doi.org/10.1007/978-1-84800-281-4
    //! [10.48550/arXiv.2406.19294]: https://doi.org/10.48550/arXiv.2406.19294
    [[nodiscard]] Presentation<word_type>
    full_transformation_monoid(size_t n,
                               author val   = author::Mitchell + author::Whyte,
                               size_t index = 0);

    //! \brief A presentation for the partial transformation monoid.
    //!
    //! Returns a monoid presentation defining the partial transformation
    //! monoid. The argument \p val determines the specific presentation which
    //! is returned. The options are:
    //! * `author::Machine`;
    //! * `author::Sutov` (see Theorem 9.4.1 of
    //! [10.1007/978-1-84800-281-4][]);
    //! * `author::Mitchell + author::Whyte`  (See Theorem 1.5 of
    //! [10.48550/arXiv.2406.19294]).
    //!
    //! \param n the degree of the partial transformation monoid.
    //! \param val the author of the presentation (default: `author::Mitchell +
    //! author::Whyte`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if \p val is not one of the options
    //! specified above.
    //!
    //! [10.1007/978-1-84800-281-4]: https://doi.org/10.1007/978-1-84800-281-4
    //! [10.48550/arXiv.2406.19294]: https://doi.org/10.48550/arXiv.2406.19294
    [[nodiscard]] Presentation<word_type> partial_transformation_monoid(
        size_t n,
        author val = author::Mitchell + author::Whyte);

    //! \brief A presentation for the symmetric inverse monoid.
    //!
    //! Returns a monoid presentation defining the symmetric inverse monoid. The
    //! argument \brief val determines the specific presentation which is
    //! returned. The options are:
    //! * `author::Gay`;
    //! * `author::Sutov` (see Theorem 9.2.2 of
    //! [10.1007/978-1-84800-281-4][]);
    //! * `author::Mitchell + author::Whyte` (see Theorem 1.4 of
    //! [10.48550/arXiv.2406.19294][]).
    //!
    //!
    //! \param n the degree of the symmetric inverse monoid.
    //! \param val the author of the presentation (default: `author::Mitchell +
    //! author::Whyte`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if \p val is not one of the options
    //! specified above.
    //!
    //! [10.1007/978-1-84800-281-4]: https://doi.org/10.1007/978-1-84800-281-4
    //! [10.48550/arXiv.2406.19294]: https://doi.org/10.48550/arXiv.2406.19294
    //
    // When val == author::Gay, this is just a presentation for the symmetric
    // inverse monoid, a slightly modified version from Solomon (so that
    // contains the Coxeter+Moser presentation for the symmetric group),
    // Example 7.1.2 in Joel gay's thesis (JDM the presentation in Example 7.1.2
    // seems to have 2n - 1 generators whereas this function returns a monoid on
    // n generators. TODO ask Florent again if this reference is correct
    // Maybe should be Solomon:
    // https://www.sciencedirect.com/science/article/pii/S0021869303005933/pdf
    [[nodiscard]] Presentation<word_type>
    symmetric_inverse_monoid(size_t n,
                             author val = author::Mitchell + author::Whyte);

    //! \brief A presentation for the Chinese monoid.
    //!
    //! Returns a monoid presentation defining the Chinese monoid, as described
    //! in [10.1142/S0218196701000425][].
    //!
    //! \param n the number of generators.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 2`.
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    //!
    //! [10.1142/S0218196701000425]: https://doi.org/10.1142/S0218196701000425
    [[nodiscard]] Presentation<word_type> chinese_monoid(size_t n,
                                                         author val
                                                         = author::Any);

    //! \brief A presentation for a monogenic semigroup.
    //!
    //! Returns a presentation defining the monogenic semigroup defined by the
    //! presentation \f$\langle a \mid a^{m + r} = a^m \rangle\f$.
    //!
    //! If \p m = 0, the presentation returned is a monoid presentation;
    //! otherwise, a semigroup presentation is returned.
    //!
    //! \param m the index.
    //! \param r the period.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `m < 0`.
    //! \throws LibsemigroupsException if `r = 0`.
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    [[nodiscard]] Presentation<word_type>
    monogenic_semigroup(size_t m, size_t r, author val = author::Any);

    //! \brief A presentation for the monoid of order-preserving mappings.
    //!
    //! Returns a monoid presentation defining the monoid of order-preserving
    //! transformations of degree \p n, as described in Section 2 of the paper
    //! [10.1007/s10012-000-0001-1][].
    //!
    //! This presentation has \f$2n - 2\f$ generators and \f$n^2\f$ relations.
    //!
    //! \param n the degree.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 3`.
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    //!
    //! [10.1007/s10012-000-0001-1]: https://doi.org/10.1007/s10012-000-0001-1
    [[nodiscard]] Presentation<word_type>
    order_preserving_monoid(size_t n, author val = author::Any);

    //! \brief A presentation for the cyclic inverse monoid.
    //!
    //! Returns a monoid presentation defining the cyclic inverse monoid of
    //! degree \p n.
    //!
    //! The arguments \p val and \p index determine the specific presentation
    //! which is returned. The options are:
    //! * `val = author::Fernandes`:
    //!    * `index = 0` (see Theorem 2.6 of [10.48550/arxiv.2211.02155][]);
    //!    * `index = 1` (see Theorem 2.7 of [10.48550/arxiv.2211.02155][]).
    //!
    //! The presentation with `val = author::Fernandes` and `index = 0` has
    //! \f$n + 1\f$ generators and \f$\frac{1}{2} \left(n^2 + 3n + 4\right)\f$
    //! relations.
    //!
    //! The presentation with `val = author::Fernandes` and `index = 1` has
    //! \f$2\f$ generators and \f$\frac{1}{2}\left(n^2 - n + 6\right)\f$
    //! relations.
    //!
    //! \param n the degree.
    //! \param val the author (default: `author::Fernandes`).
    //! \param index the index (default: `1`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 3`.
    //! \throws LibsemigroupsException if `val` is not `author::Fernandes`.
    //! \throws LibsemigroupsException if `val == author::Fernandes` and \p
    //! index is not `0` or `1`.
    //!
    //! [10.48550/arxiv.2211.02155]: https://doi.org/10.48550/arxiv.2211.02155
    [[nodiscard]] Presentation<word_type>
    cyclic_inverse_monoid(size_t n,
                          author val   = author::Fernandes,
                          size_t index = 1);

    //! \p A presentation for the order-preserving part of the cyclic inverse
    //! monoid.
    //!
    //! Returns a monoid presentation defining the order-preserving part of the
    //! cyclic inverse monoid of degree \p n, as described in Theorem 2.17 of
    //! the paper [10.48550/arxiv.2211.02155][].
    //!
    //! \param n the degree.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 3`.
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    //!
    //! [10.48550/arxiv.2211.02155]: https://doi.org/10.48550/arxiv.2211.02155
    [[nodiscard]] Presentation<word_type>
    order_preserving_cyclic_inverse_monoid(size_t n, author val = author::Any);

    //! \brief A presentation for the monoid of partial isometries of a cycle
    //! graph.
    //!
    //! Returns a monoid presentation defining the monoid of partial isometries
    //! of an \f$n\f$-cycle graph, as described in Theorem 2.8 of
    //! [10.48550/arxiv.2205.02196][].
    //!
    //! \param n the number of vertices of the cycle graph.
    //! \param val the author of the presentation (default: `author::Any`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 3`.
    //! \throws LibsemigroupsException if \p val is not `author::Any`.
    //!
    //! [10.48550/arxiv.2205.02196]: https://doi.org/10.48550/arxiv.2205.02196
    [[nodiscard]] Presentation<word_type>
    partial_isometries_cycle_graph_monoid(size_t n, author val = author::Any);

    //! \brief A non-presentation for the symmetric group.
    //!
    //! Returns a monoid presentation which is claimed to define the symmetric
    //! group of degree \p n, but does not. The argument \p val determines the
    //! specific non-presentation which is returned. The options are:
    //! * `author::Guralnick + author::Kantor + author::Kassabov +
    //! author::Lubotzky` [doi.org/10.1090/S0894-0347-08-00590-0][].
    //!
    //! \param n the claimed degree of the symmetric group.
    //! \param val the author of the presentation (default: `author::Guralnick +
    //! author::Kantor + author::Kassabov + author::Lubotzky`).
    //!
    //! \returns A value of type `Presentation<libsemigroups::word_type>`.
    //!
    //! \throws LibsemigroupsException if `n < 4`.
    //! \throws LibsemigroupsException if `val` is not listed above.
    //!
    //! [doi.org/10.1090/S0894-0347-08-00590-0]:
    //! https://doi.org/10.1090/S0894-0347-08-00590-0
    [[nodiscard]] Presentation<word_type>
    not_symmetric_group(size_t n,
                        author val = author::Guralnick + author::Kantor
                                     + author::Kassabov + author::Lubotzky);

    // TODO (doc)
    // n should be prime for this presentation to actually defined the claimed
    // group.
    [[nodiscard]] Presentation<word_type> special_linear_group_2(size_t n,
                                                                 author val
                                                                 = author::Any);

    // TODO (doc)
    [[nodiscard]] Presentation<word_type> hypo_plactic_monoid(size_t n,
                                                              author val
                                                              = author::Any);

    // TODO (doc)
    [[nodiscard]] Presentation<word_type>
    sigma_stylic_monoid(std::vector<size_t> const& sigma,
                        author                     val = author::Any);

    // TODO add okada_monoid
    // TODO add free_semilattice

    // TODO update doc
    // When q = 0, the a_i^2 = 1 relations from the C+M symmetric group
    // presentation are replaced with a_i^2 = a_i. See Definition 4.1.1 in Joel
    // Gay's thesis
    // https://theses.hal.science/tel-01861199
    //
    //  This could also be called the RennerTypeAMonoid
    [[nodiscard]] Presentation<word_type> zero_rook_monoid(size_t n,
                                                           author val
                                                           = author::Any);

    // TODO rename renner_type_b_monoid:
    // when q = 0: Definition 8.4.1 + Example 8.4.2 in Joel's thesis.
    // when q = 1: same presentation as for q=0 except pi_i ^ 2 = pi_i replaced
    // with pi_i ^ 2 = 1 (usual Renner monoid) (q is the Iwahori-Hecke
    // deformation of the Renner monoid). q = 1From Joel's thesis,
    // Theorem 8.4.19.
    [[nodiscard]] Presentation<word_type>
    renner_type_B_monoid(size_t l, int q, author val = author::Any);

    // not_renner_type_b_monoid
    // when q = 0, Example 7.1.2, of Joel's thesis
    // when q = 1, is Godelle's presentation from: p40 of
    // https://www.cambridge.org/core/services/aop-cambridge-core/content/view/B6BAB75BD3463916FEDEC15BEDA724FF/S0004972710000365a.pdf/presentation_for_renner_monoids.pdf
    // or
    // https://arxiv.org/abs/0904.0926
    [[nodiscard]] Presentation<word_type>
    not_renner_type_B_monoid(size_t l, int q, author val = author::Any);

    // See Theorem 8.4.43 in Joel's thesis q = 1,
    // for q = 0, Definition 8.4.22 (author = Gay)
    [[nodiscard]] Presentation<word_type>
    renner_type_D_monoid(size_t l, int q, author val = author::Any);

    // Godelle's:
    // when q = 1, p41 of
    // https://www.cambridge.org/core/services/aop-cambridge-core/content/view/B6BAB75BD3463916FEDEC15BEDA724FF/S0004972710000365a.pdf/presentation_for_renner_monoids.pdf
    // q = 0, no reference
    [[nodiscard]] Presentation<word_type>
    not_renner_type_D_monoid(size_t l, int q, author val = author::Any);

  }  // namespace fpsemigroup
}  // namespace libsemigroups
#endif  // LIBSEMIGROUPS_FPSEMI_EXAMPLES_HPP_
