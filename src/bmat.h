//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2017 Finn Smith
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

// This file contains a declaration of fast boolean matrices up to dimension 8.

#ifndef LIBSEMIGROUPS_SRC_BMAT_H_
#define LIBSEMIGROUPS_SRC_BMAT_H_

#include <climits>
#include <functional>
#include <iostream>
#include <random>

#include "libsemigroups-debug.h"
#include "timer.h"

namespace libsemigroups {

  //! Class for fast boolean matrices of dimension up to 8x8
  //!
  //! The methods for these small matrices over the boolean semiring
  //! are more optimised than the generic methods for boolean matrices.
  //! Note that all BMat8 are represented internally as an 8x8 matrix;
  //! any entries not defined by the user are taken to be 0. This does
  //! not affect the results of any calculations.
  class BMat8 {
   public:
    //! A default constructor.
    //!
    //! This constructor gives no guarantees on what the matrix will contain.
    BMat8() = default;

    //! A constructor.
    //!
    //! The 8 8-bit chunks of the binary representation of \p mat
    //! form the rows of the matrix initialized by this constructor.
    explicit BMat8(uint64_t mat) : _data(mat) {}

    //! A constructor.
    //!
    //! This constructor initializes a matrix where the rows of the matrix
    //! are the vectors in \p mat.
    explicit BMat8(std::vector<std::vector<size_t>> const &mat);

    //! A constructor.
    //!
    //! This is the copy constructor.
    BMat8(BMat8 const &) = default;

    //! A constructor.
    //!
    //! This is the move constructor.
    BMat8(BMat8 &&) = default;

    //! A constructor.
    //!
    //! This is the copy assignement constructor.
    BMat8 &operator=(BMat8 const &) = default;

    //! A constructor.
    //!
    //! This is the move assignment  constructor.
    BMat8 &operator=(BMat8 &&) = default;

    //! A default destructor.
    ~BMat8() = default;

    //! Returns \c true if \c this equals \p that.
    //!
    //! This method checks the mathematical equality of two BMat8 objects.
    bool operator==(BMat8 const &that) const { return _data == that._data; }

    //! Returns \c true if \c this does not equal \p that
    //!
    //! This method checks the mathematical inequality of two BMat8 objects.
    bool operator!=(BMat8 const &that) const { return _data != that._data; }

    //! Returns \c true if \c this is less than \p that.
    //!
    //! This method checks whether a BMat8 objects is less than another.
    //! We order by the results of to_int() for each matrix.
    bool operator<(BMat8 const &that) const { return _data < that._data; }

    //! Returns \c true if \c this is greater than \p that.
    //!
    //! This method checks whether a BMat8 objects is greater than another.
    //! We order by the results of to_int() for each matrix.
    bool operator>(BMat8 const &that) const { return _data > that._data; }

    //! Returns the entry in the (\p i, \p j)th position.
    //!
    //! This method returns the entry in the (\p i, \p j)th position.
    //! Note that since all matrices are internally represented as 8x8, it
    //! is possible to access entries that you might not believe exist.
    bool operator()(size_t i, size_t j) const {
      LIBSEMIGROUPS_ASSERT(0 <= i && i < 8);
      LIBSEMIGROUPS_ASSERT(0 <= j && j < 8);
      return (_data << (8 * i + j)) >> 63;
    }

    //! Returns the integer representation of \c this.
    //!
    //! Returns an unsigned integer obtained by interpreting an 8x8
    //! BMat8 as a sequence of 64 bits (reading rows left to right,
    //! from top to bottom)and then this sequence as an unsigned int.
    inline uint64_t to_int() const { return _data; }

    //! Returns the transpose of \c this
    //!
    //! Returns the standard matrix transpose of a BMat8.
    inline BMat8 transpose() const {
      uint64_t x = _data;
      uint64_t y = (x ^ (x >> 7)) & 0xAA00AA00AA00AA;
      x = x ^ y ^ (y << 7);
      y = (x ^ (x >> 14)) & 0xCCCC0000CCCC;
      x = x ^ y ^ (y << 14);
      y = (x ^ (x >> 28)) & 0xF0F0F0F0;
      x = x ^ y ^ (y << 28);
      return BMat8(x);
    }

    //! Returns the matrix product of \c this and \p that
    //!
    //! This method returns the standard matrix product (over the
    //! boolean semiring) of two BMat8.
    // https://stackoverflow.com/a/18448513
    inline BMat8 operator*(BMat8 const &that) const {
      uint64_t y = that.transpose()._data;
      uint64_t data = 0;
      uint64_t tmp = 0;
      uint64_t diag = 0x8040201008040201;
      for (int i = 0; i < 8; ++i) {
        tmp = _data & y;
        tmp |= tmp >> 1;
        tmp |= tmp >> 2;
        tmp |= tmp >> 4;
        tmp &= 0x0101010101010101;
        tmp *= 255;
        tmp &= diag;
        data |= tmp;
        y = cyclic_shift(y);
        tmp = 0;
        diag = cyclic_shift(diag);
      }
      return BMat8(data);
    }

    //! Returns the identity BMat8
    //!
    //! This method returns the 8x8 BMat8 with 1s on the main diagonal.
    inline BMat8 one() const { return BMat8(0x8040201008040201); }

    //! Insertion operator
    //!
    //! This method allows BMat8s to be inserted into a stream.
    friend std::ostream &operator<<(std::ostream &os, BMat8 const &bm) {
      os << bm.to_string();
      return os;
    }

    //! Returns a string representation of \c this
    //!
    //! This method returns a string which represents a BMat8, which may
    //! for example be used to display the BMat8.
    std::string to_string() const;

    //! Returns a random BMat8
    //!
    //! This method returns a BMat8 chosen at random.
    static BMat8 random();

    //! Returns a random square BMat8 up to dimension \p dim.
    //!
    //! This method returns a BMat8 chosen at random, where only the
    //! top-left \p dim x \p dim entries may be non-zero.
    static BMat8 random(size_t dim);

   private:
    uint64_t _data;
    static std::random_device _rd;
    static std::mt19937 _gen;
    static std::uniform_int_distribution<size_t> _dist;
    static std::vector<uint64_t> const ROW_MASK;
    static std::vector<uint64_t> const COL_MASK;

    // Cyclically shifts bits to left by 8m
    // https://stackoverflow.com/a/776523
    static inline uint64_t cyclic_shift(uint64_t n, uint64_t m = 1) {
      const unsigned int mask =
          (CHAR_BIT * sizeof(n) - 1);  // assumes width is a power of 2.

      // assert ( (c<=mask) &&"rotate by type width or more");
      unsigned int c = 8 * m;
      c &= mask;
      return (n << c) | (n >> ((-c) & mask));
    }
  };
}  // namespace libsemigroups

namespace std {
  template <> struct hash<libsemigroups::BMat8> {
    size_t operator()(libsemigroups::BMat8 const& bm) const {
      return hash<uint64_t>()(bm.to_int());
    }
  };
}  // namespace std
#endif  // LIBSEMIGROUPS_SRC_BMAT_H_
