//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2019  Jonas Klappert and Fabian Lange
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//====================================================================================

#pragma once

#include <gmpxx.h>
#include <iostream>
#include <string>
#include <vector>

namespace firefly {

  /**
  * @class FFInt
  * @brief A class for finite field integers
  */
  class FFInt {
  public:
    /**
     *    A constructor
     *    @param n_ an integer which is a member of the finite field
     */
    template<typename T>
    FFInt(const T n_);
    /**
     *    A constructor
     *    @param ffint a FFInt object
     */
    FFInt(const FFInt& ffint);
    /**
     *  A constructer for a mpz_class object
     *  @param in the mpz_class object which should be converted to an FFInt
     */
    FFInt(mpz_class in);
    FFInt(const std::string& str, const std::vector<std::pair<std::string, uint64_t>>& replacements);
    /**
     *    Default constructor
     */
    FFInt();
    /**
     *  A static function to define a new prime field
     *  @param prime the defining prime of the field
     */
    static void set_new_prime(uint64_t prime);

    // defining new operators for finite field arithmetic
    FFInt& operator=(const FFInt&) = default;
    FFInt& operator+=(const FFInt&);
    FFInt& operator-=(const FFInt&);
    FFInt& operator*=(const FFInt&);
    FFInt& operator/=(const FFInt&);
    FFInt operator-() const;
    FFInt operator+() const;
    FFInt operator++();
    FFInt operator++(int);
    FFInt operator--();
    FFInt operator--(int);
    bool operator!() const;
    FFInt pow(const FFInt& ffint) const;

    uint64_t n; /**< the integer member of the finite field */
    static uint64_t p; /**< the prime defining the finite field */
    static uint64_t p_inv; /**< the inverse of the defining prime needed for FFTs*/
  private:
    uint64_t parse_longint(const std::string& str) const;
  };

  bool operator<(const FFInt& a, const FFInt& b);
  bool operator<=(const FFInt& a, const FFInt& b);
  bool operator>(const FFInt& a, const FFInt& b);
  bool operator>=(const FFInt& a, const FFInt& b);
  bool operator==(const FFInt& a, const FFInt& b);
  bool operator!=(const FFInt& a, const FFInt& b);
  FFInt operator/(const FFInt& a, const FFInt& b);
  FFInt operator+(const FFInt& a, const FFInt& b);
  FFInt operator-(const FFInt& a, const FFInt& b);
  FFInt operator*(const FFInt& a, const FFInt& b);
  FFInt pow(const FFInt& ffint, const FFInt& power);
  std::ostream& operator<<(std::ostream& out, const FFInt& ffint);

  template<typename T>
  FFInt::FFInt(const T n_) {
    if (n_ >= 0) {
      if (static_cast<uint64_t>(n_) < p) {
        n = n_;
      } else {
        n = n_ % p;
      }
    } else if (n_ < 0) {
      n = p - (-n_) % p;
    }
  }

  extern "C" {
    void firefly_exists(void);
  }
}
