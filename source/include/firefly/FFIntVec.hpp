//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2020  Jonas Klappert and Fabian Lange
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

#include "firefly/config.hpp"
#include "firefly/FFInt.hpp"

#include <array>

namespace firefly {

  /**
   * @class FFIntVec
   * @brief A class for finite field integers stored as an array for vector arithmetic
   */
  template<int N>
  class FFIntVec {
  template<int NN> friend FFIntVec<NN> operator-(const FFIntVec<NN> &);
  template<int NN> friend FFIntVec<NN> operator+(const FFIntVec<NN> &, const FFIntVec<NN> &);
  template<int NN> friend FFIntVec<NN> operator-(const FFIntVec<NN> &, const FFIntVec<NN> &);
  template<int NN> friend FFIntVec<NN> operator*(const FFIntVec<NN> &, const FFIntVec<NN> &);
  template<int NN> friend FFIntVec<NN> operator/(const FFIntVec<NN> &, const FFIntVec<NN> &);
  template<int NN> friend FFIntVec<NN> operator+(const FFIntVec<NN> &, const FFInt &);
  template<int NN> friend FFIntVec<NN> operator-(const FFIntVec<NN> &, const FFInt &);
  template<int NN> friend FFIntVec<NN> operator*(const FFIntVec<NN> &, const FFInt &);
  template<int NN> friend FFIntVec<NN> operator/(const FFIntVec<NN> &, const FFInt &);
  template<int NN> friend FFIntVec<NN> operator+(const FFInt &, const FFIntVec<NN> &);
  template<int NN> friend FFIntVec<NN> operator-(const FFInt &, const FFIntVec<NN> &);
  template<int NN> friend FFIntVec<NN> operator*(const FFInt &, const FFIntVec<NN> &);
  template<int NN> friend FFIntVec<NN> operator/(const FFInt &, const FFIntVec<NN> &);
  template<int NN> friend FFIntVec<NN> pow(const FFIntVec<NN> &, const FFIntVec<NN> &);
  template<int NN> friend bool operator==(const FFIntVec<NN> &, const FFIntVec<NN> &);
  template<int NN> friend bool operator!=(const FFIntVec<NN> &, const FFIntVec<NN> &);
  template<int NN> friend std::ostream & operator<<(std::ostream &, const FFIntVec<NN> &);
  public:
    /**
     *  Constructor
     *  @param in the initializing array
     */
    FFIntVec(const std::array<FFInt, N>& in) : vec(in) {};
    /**
     *  Constructor
     *  @param in the initializing array
     */
    FFIntVec(const FFIntVec& in) : vec(in.vec) {};
    /**
     *  Constructor
     *  @param in fills the array with this FFInt
     */
    FFIntVec(const FFInt& in) {
      vec.fill(in);
    };
    /**
     *  Default constructor
     */
    FFIntVec() {
      vec.fill(0);
    };
    /**
     *  Converts to a negative integer (only for negative exponents)
     */
    int to_neg_int() const noexcept;

    std::array<FFInt, N> vec; /**< The stored vector for arithmetic */
    // defining new operators for finite field arithmetic
    FFIntVec& operator=(const FFIntVec&) = default;
    FFIntVec& operator=(FFIntVec&&) = default;
    FFIntVec& operator+=(const FFIntVec&);
    FFIntVec& operator-=(const FFIntVec&);
    FFIntVec& operator*=(const FFIntVec&);
    FFIntVec& operator/=(const FFIntVec&);
    FFIntVec pow(const FFIntVec& power) const;
    template<class T, typename=typename std::enable_if<(std::is_enum<T>::value || std::is_integral<T>::value)>::type>
    FFIntVec pow(const std::vector<T>& power) const;
    FFIntVec pow(const FFInt& power) const;
    template<class T, typename=typename std::enable_if<(std::is_enum<T>::value || std::is_integral<T>::value)>::type>
    FFIntVec pow(const T& power) const;
    FFInt& operator[](int i) {return vec[i];};
    FFInt operator[](int i) const {return vec[i];};
    bool operator!() const;
    typename std::array<FFInt, N>::iterator begin() noexcept;
    typename std::array<FFInt, N>::const_iterator begin() const noexcept;
    typename std::array<FFInt, N>::iterator end() noexcept;
    typename std::array<FFInt, N>::const_iterator end() const noexcept;
    size_t size() const noexcept;
    inline const FFInt& at(int i) {
      if (i < N - 1)
        return vec[i];
      else
        throw std::out_of_range("Out of range.");
    }
    //static uint32_t size; /**< Sets the size of the vectors */
  };

  template<int N>
  template<class T, typename>
  FFIntVec<N> FFIntVec<N>::pow(const std::vector<T>& power) const {
    if (power.size() != N) {
      throw std::out_of_range("Power of FFIntVec with vector of wrong length.");
    }

    FFIntVec<N> result(vec);
    for (int i = 0; i != N; ++i) {
      result.vec[i] = result.vec[i].pow(power[i]);
    }

    return result;
  }

  template<int N>
  template<class T, typename>
  FFIntVec<N> FFIntVec<N>::pow(const T& power) const {
    FFIntVec<N> result(vec);
    for (int i = 0; i != N; ++i) {
      result.vec[i] = result.vec[i].pow(power);
    }

    return result;
  }

  template<int N>
  inline bool FFIntVec<N>::operator!() const {
    for(const auto& el : vec) {
      if (!el)
        return true;
    }

    return false;
  }

  template<int N>
  bool operator==(const FFIntVec<N>& a, const FFIntVec<N>& b);
  template<int N>
  bool operator!=(const FFIntVec<N>& a, const FFIntVec<N>& b);
  template<int N>
  bool operator>(const FFIntVec<N>& a, const FFIntVec<N>& b);
  template<int N>
  bool operator>=(const FFIntVec<N>& a, const FFIntVec<N>& b);
  template<int N>
  bool operator<(const FFIntVec<N>& a, const FFIntVec<N>& b);
  template<int N>
  bool operator<=(const FFIntVec<N>& a, const FFIntVec<N>& b);
  template<int N>
  FFIntVec<N> operator/(const FFIntVec<N>& a, const FFIntVec<N>& b);
  template<int N>
  FFIntVec<N> operator+(const FFIntVec<N>& a);
  template<int N>
  FFIntVec<N> operator+(const FFIntVec<N>& a, const FFIntVec<N>& b);
  template<int N>
  FFIntVec<N> operator-(const FFIntVec<N>& a);
  template<int N>
  FFIntVec<N> operator-(const FFIntVec<N>& a, const FFIntVec<N>& b);
  template<int N>
  FFIntVec<N> operator*(const FFIntVec<N>& a, const FFIntVec<N>& b);
  template<int N>
  FFIntVec<N> pow(const FFIntVec<N>& a, const FFIntVec<N>& power);
  template<int N>
  std::ostream& operator<<(std::ostream& out, const FFIntVec<N>& ffint_vec);
}
