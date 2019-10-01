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

#include "FFInt.hpp"

namespace firefly {

  /**
   * @class FFIntVec
   * @brief A class for finite field integers stored as a vector for vector arithmetic
   */
  class FFIntVec {
  public:
    /**
     *  Constructor
     *  @param in the initializing vector
     */
    FFIntVec(const std::vector<FFInt>& in);
    /**
     *  Constructor
     *  @param in an initializing FFInt casted to an array with this entry
     */
    FFIntVec(const FFInt& in);
    /**
     *  Default constructor
     */
    FFIntVec(){};

    // defining new operators for finite field arithmetic
    FFIntVec& operator=(const FFIntVec&) = default;
    FFIntVec& operator+=(const FFIntVec&);
    FFIntVec& operator-=(const FFIntVec&);
    FFIntVec& operator*=(const FFIntVec&);
    FFIntVec& operator/=(const FFIntVec&);
    FFIntVec pow(const FFIntVec& power) const;
    bool operator!() const;
    std::vector<FFInt> vec; /**< The stored vector for arithmetic */
    static uint32_t size; /**< Sets the size of the vectors */
  };

  inline bool FFIntVec::operator!() const {
    for(const auto& el : vec) {
      if (!el)
        return true;
    }

    return false;
  }

  bool operator==(const FFIntVec& a, const FFIntVec& b);
  bool operator!=(const FFIntVec& a, const FFIntVec& b);
  FFIntVec operator/(const FFIntVec& a, const FFIntVec& b);
  FFIntVec operator+(const FFIntVec& a, const FFIntVec& b);
  FFIntVec operator-(const FFIntVec& a, const FFIntVec& b);
  FFIntVec operator*(const FFIntVec& a, const FFIntVec& b);
  FFIntVec pow(const FFIntVec& a, const FFIntVec& power);
  std::ostream& operator<<(std::ostream& out, const FFIntVec& ffint_vec);
}
