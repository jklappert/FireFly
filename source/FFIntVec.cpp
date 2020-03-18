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

#include "FFIntVec.hpp"

namespace firefly {

  template<int N>
  typename std::array<FFInt, N>::const_iterator FFIntVec<N>::begin() const noexcept {
    return vec.begin();
  }

  template<int N>
  typename std::array<FFInt, N>::iterator FFIntVec<N>::begin() noexcept {
    return vec.begin();
  }

  template<int N>
  typename std::array<FFInt, N>::const_iterator FFIntVec<N>::end() const noexcept {
    return vec.end();
  }

  template<int N>
  typename std::array<FFInt, N>::iterator FFIntVec<N>::end() noexcept {
    return vec.end();
  }

  template<int N>
  size_t FFIntVec<N>::size() const noexcept {
    return static_cast<size_t>(N);
  }

  template<int N>
  bool operator==(const FFIntVec<N>& a, const FFIntVec<N>& b) {
    return a.vec == b.vec;
  }

  template<int N>
  bool operator!=(const FFIntVec<N>& a, const FFIntVec<N>& b) {
    return a.vec != b.vec;
  }

  // TODO
  template<int N>
  bool operator>(const FFIntVec<N>& a, const FFIntVec<N>& b) {
    return a.vec[0] > b.vec[0];
  }

  // TODO
  template<int N>
  bool operator>=(const FFIntVec<N>& a, const FFIntVec<N>& b) {
    return a.vec[0] >= b.vec[0];
  }

  // TODO
  template<int N>
  bool operator<(const FFIntVec<N>& a, const FFIntVec<N>& b) {
    return a.vec[0] < b.vec[0];
  }

  // TODO
  template<int N>
  bool operator<=(const FFIntVec<N>& a, const FFIntVec<N>& b) {
    return a.vec[0] <= b.vec[0];
  }

  template<int N>
  FFIntVec<N>& FFIntVec<N>::operator+=(const FFIntVec<N>& a) {
    for (int i = 0; i != N; ++i) {
      vec[i] += a.vec[i];
    }

    return *this;
  }

  template<int N>
  FFIntVec<N>& FFIntVec<N>::operator-=(const FFIntVec<N>& a) {
    for (int i = 0; i != N; ++i) {
      vec[i] -= a.vec[i];
    }

    return *this;
  }

  template<int N>
  FFIntVec<N>& FFIntVec<N>::operator*=(const FFIntVec<N>& a) {
    for (int i = 0; i != N; ++i) {
      vec[i] *= a.vec[i];
    }

    return *this;
  }

  template<int N>
  FFIntVec<N>& FFIntVec<N>::operator/=(const FFIntVec<N>& a) {
    for (int i = 0; i != N; ++i) {
      vec[i] /= a.vec[i];
    }

    return *this;
  }

  template<int N>
  FFIntVec<N> FFIntVec<N>::pow(const FFIntVec<N>& power) const {
    FFIntVec<N> result(vec);
    for (int i = 0; i != N; ++i) {
      result.vec[i] = result.vec[i].pow(power.vec[i]);
    }

    return result;
  }

  template<int N>
  FFIntVec<N> FFIntVec<N>::pow(const FFInt& power) const {
    FFIntVec<N> result(vec);
    for (int i = 0; i != N; ++i) {
      result.vec[i] = result.vec[i].pow(power);
    }

    return result;
  }

  template<int N>
  FFIntVec<N> operator+(const FFIntVec<N>& a, const FFIntVec<N>& b) {
    FFIntVec<N> result(a);

    for (int i = 0; i != N; ++i) {
      result.vec[i] += b.vec[i];
    }

    return result;
  }

  template<int N>
  FFIntVec<N> operator+(const FFIntVec<N>& a) {
    FFIntVec<N> result;

    for (int i = 0; i != N; ++i) {
      result.vec[i] = +a.vec[i];
    }

    return result;
  }

  template<int N>
  FFIntVec<N> operator-(const FFIntVec<N>& a, const FFIntVec<N>& b) {
    FFIntVec<N> result(a);

    for (int i = 0; i != N; ++i) {
      result.vec[i] -= b.vec[i];
    }

    return result;
  }

  template<int N>
  FFIntVec<N> operator-(const FFIntVec<N>& a) {
    FFIntVec<N> result;

    for (int i = 0; i != N; ++i) {
      result.vec[i] = -a.vec[i];
    }

    return result;
  }

  template<int N>
  FFIntVec<N> operator*(const FFIntVec<N>& a, const FFIntVec<N>& b) {
    FFIntVec<N> result(a);

    for (int i = 0; i != N; ++i) {
      result.vec[i] *= b.vec[i];
    }

    return result;
  }

  template<int N>
  FFIntVec<N> operator/(const FFIntVec<N>& a, const FFIntVec<N>& b) {
    FFIntVec<N> result(a);

    for (int i = 0; i != N; ++i) {
      result.vec[i] /= b.vec[i];
    }

    return result;
  }

  template<int N>
  FFIntVec<N> operator+(const FFIntVec<N>& a, const FFInt& b) {
    FFIntVec<N> result(a);

    for (int i = 0; i != N; ++i) {
      result.vec[i] += b;
    }

    return result;
  }

  template<int N>
  FFIntVec<N> operator-(const FFIntVec<N>& a, const FFInt& b) {
    FFIntVec<N> result(a);

    for (int i = 0; i != N; ++i) {
      result.vec[i] -= b;
    }

    return result;
  }

  template<int N>
  FFIntVec<N> operator*(const FFIntVec<N>& a, const FFInt& b) {
    FFIntVec<N> result(a);

    for (int i = 0; i != N; ++i) {
      result.vec[i] *= b;
    }

    return result;
  }

  template<int N>
  FFIntVec<N> operator/(const FFIntVec<N>& a, const FFInt& b) {
    FFIntVec<N> result(a);

    for (int i = 0; i != N; ++i) {
      result.vec[i] /= b;
    }

    return result;
  }

  template<int N>
  FFIntVec<N> operator+(const FFInt& b, const FFIntVec<N>& a) {
    FFIntVec<N> result(a);

    for (int i = 0; i != N; ++i) {
      result.vec[i] += b;
    }

    return result;
  }

  template<int N>
  FFIntVec<N> operator-(const FFInt& b, const FFIntVec<N>& a) {
    FFIntVec<N> result(a);

    for (int i = 0; i != N; ++i) {
      result.vec[i] -= b;
    }

    return result;
  }

  template<int N>
  FFIntVec<N> operator*(const FFInt& b, const FFIntVec<N>& a) {
    FFIntVec<N> result(a);

    for (int i = 0; i != N; ++i) {
      result.vec[i] *= b;
    }

    return result;
  }

  template<int N>
  FFIntVec<N> operator/(const FFInt& b, const FFIntVec<N>& a) {
    FFIntVec<N> result(a);

    for (int i = 0; i != N; ++i) {
      result.vec[i] /= b;
    }

    return result;
  }

  template<int N>
  FFIntVec<N> pow(const FFIntVec<N>& a, const FFIntVec<N>& power) {
    FFIntVec<N> result(a);

    for (int i = 0; i != N; ++i) {
      result.vec[i] = result.vec[i].pow(power.vec[i]);
    }

    return result;
  }

  template<int N>
  std::ostream& operator<<(std::ostream& out, const FFIntVec<N>& ffint_vec) {
    out << "(" << ffint_vec.vec.front();

    for (auto it = ffint_vec.vec.begin() + 1; it != ffint_vec.vec.end(); ++ it) {
      out << ", " << *it;
    }

    out << ")";
    return out;
  }

  // 1
  /*template class FFIntVec<1>;
  template FFIntVec<1> operator-(const FFIntVec<1>&);
  template FFIntVec<1> operator+(const FFIntVec<1>&);
  template FFIntVec<1> operator-(const FFIntVec<1>&, const FFIntVec<1>&);
  template FFIntVec<1> operator+(const FFIntVec<1>&, const FFIntVec<1>&);
  template FFIntVec<1> operator*(const FFIntVec<1>&, const FFIntVec<1>&);
  template FFIntVec<1> operator/(const FFIntVec<1>&, const FFIntVec<1>&);
  template FFIntVec<1> operator+(const FFIntVec<1>&, const FFInt&);
  template FFIntVec<1> operator-(const FFIntVec<1>&, const FFInt&);
  template FFIntVec<1> operator*(const FFIntVec<1>&, const FFInt&);
  template FFIntVec<1> operator/(const FFIntVec<1>&, const FFInt&);
  template FFIntVec<1> operator+(const FFInt&, const FFIntVec<1>&);
  template FFIntVec<1> operator-(const FFInt&, const FFIntVec<1>&);
  template FFIntVec<1> operator*(const FFInt&, const FFIntVec<1>&);
  template FFIntVec<1> operator/(const FFInt&, const FFIntVec<1>&);
  template bool operator==(const FFIntVec<1>&, const FFIntVec<1>&);
  template bool operator!=(const FFIntVec<1>&, const FFIntVec<1>&);
  template bool operator>(const FFIntVec<1>&, const FFIntVec<1>&);
  template bool operator>=(const FFIntVec<1>&, const FFIntVec<1>&);
  template bool operator<(const FFIntVec<1>&, const FFIntVec<1>&);
  template bool operator<=(const FFIntVec<1>&, const FFIntVec<1>&);
  template FFIntVec<1> pow(const FFIntVec<1>&, const FFIntVec<1>&);
  template std::ostream& operator<<(std::ostream&, const FFIntVec<1>&);*/
  // 2
  template class FFIntVec<2>;
  template FFIntVec<2> operator-(const FFIntVec<2>&);
  template FFIntVec<2> operator+(const FFIntVec<2>&);
  template FFIntVec<2> operator-(const FFIntVec<2>&, const FFIntVec<2>&);
  template FFIntVec<2> operator+(const FFIntVec<2>&, const FFIntVec<2>&);
  template FFIntVec<2> operator*(const FFIntVec<2>&, const FFIntVec<2>&);
  template FFIntVec<2> operator/(const FFIntVec<2>&, const FFIntVec<2>&);
  template FFIntVec<2> operator+(const FFIntVec<2>&, const FFInt&);
  template FFIntVec<2> operator-(const FFIntVec<2>&, const FFInt&);
  template FFIntVec<2> operator*(const FFIntVec<2>&, const FFInt&);
  template FFIntVec<2> operator/(const FFIntVec<2>&, const FFInt&);
  template FFIntVec<2> operator+(const FFInt&, const FFIntVec<2>&);
  template FFIntVec<2> operator-(const FFInt&, const FFIntVec<2>&);
  template FFIntVec<2> operator*(const FFInt&, const FFIntVec<2>&);
  template FFIntVec<2> operator/(const FFInt&, const FFIntVec<2>&);
  template bool operator==(const FFIntVec<2>&, const FFIntVec<2>&);
  template bool operator!=(const FFIntVec<2>&, const FFIntVec<2>&);
  template bool operator>(const FFIntVec<2>&, const FFIntVec<2>&);
  template bool operator>=(const FFIntVec<2>&, const FFIntVec<2>&);
  template bool operator<(const FFIntVec<2>&, const FFIntVec<2>&);
  template bool operator<=(const FFIntVec<2>&, const FFIntVec<2>&);
  template FFIntVec<2> pow(const FFIntVec<2>&, const FFIntVec<2>&);
  template std::ostream& operator<<(std::ostream&, const FFIntVec<2>&);
  // 4
  template class FFIntVec<4>;
  template FFIntVec<4> operator-(const FFIntVec<4>&);
  template FFIntVec<4> operator+(const FFIntVec<4>&);
  template FFIntVec<4> operator-(const FFIntVec<4>&, const FFIntVec<4>&);
  template FFIntVec<4> operator+(const FFIntVec<4>&, const FFIntVec<4>&);
  template FFIntVec<4> operator*(const FFIntVec<4>&, const FFIntVec<4>&);
  template FFIntVec<4> operator/(const FFIntVec<4>&, const FFIntVec<4>&);
  template FFIntVec<4> operator+(const FFIntVec<4>&, const FFInt&);
  template FFIntVec<4> operator-(const FFIntVec<4>&, const FFInt&);
  template FFIntVec<4> operator*(const FFIntVec<4>&, const FFInt&);
  template FFIntVec<4> operator/(const FFIntVec<4>&, const FFInt&);
  template FFIntVec<4> operator+(const FFInt&, const FFIntVec<4>&);
  template FFIntVec<4> operator-(const FFInt&, const FFIntVec<4>&);
  template FFIntVec<4> operator*(const FFInt&, const FFIntVec<4>&);
  template FFIntVec<4> operator/(const FFInt&, const FFIntVec<4>&);
  template bool operator==(const FFIntVec<4>&, const FFIntVec<4>&);
  template bool operator!=(const FFIntVec<4>&, const FFIntVec<4>&);
  template bool operator>(const FFIntVec<4>&, const FFIntVec<4>&);
  template bool operator>=(const FFIntVec<4>&, const FFIntVec<4>&);
  template bool operator<(const FFIntVec<4>&, const FFIntVec<4>&);
  template bool operator<=(const FFIntVec<4>&, const FFIntVec<4>&);
  template FFIntVec<4> pow(const FFIntVec<4>&, const FFIntVec<4>&);
  template std::ostream& operator<<(std::ostream&, const FFIntVec<4>&);
  // 8
  template class FFIntVec<8>;
  template FFIntVec<8> operator-(const FFIntVec<8>&);
  template FFIntVec<8> operator+(const FFIntVec<8>&);
  template FFIntVec<8> operator-(const FFIntVec<8>&, const FFIntVec<8>&);
  template FFIntVec<8> operator+(const FFIntVec<8>&, const FFIntVec<8>&);
  template FFIntVec<8> operator*(const FFIntVec<8>&, const FFIntVec<8>&);
  template FFIntVec<8> operator/(const FFIntVec<8>&, const FFIntVec<8>&);
  template FFIntVec<8> operator+(const FFIntVec<8>&, const FFInt&);
  template FFIntVec<8> operator-(const FFIntVec<8>&, const FFInt&);
  template FFIntVec<8> operator*(const FFIntVec<8>&, const FFInt&);
  template FFIntVec<8> operator/(const FFIntVec<8>&, const FFInt&);
  template FFIntVec<8> operator+(const FFInt&, const FFIntVec<8>&);
  template FFIntVec<8> operator-(const FFInt&, const FFIntVec<8>&);
  template FFIntVec<8> operator*(const FFInt&, const FFIntVec<8>&);
  template FFIntVec<8> operator/(const FFInt&, const FFIntVec<8>&);
  template bool operator==(const FFIntVec<8>&, const FFIntVec<8>&);
  template bool operator!=(const FFIntVec<8>&, const FFIntVec<8>&);
  template bool operator>(const FFIntVec<8>&, const FFIntVec<8>&);
  template bool operator>=(const FFIntVec<8>&, const FFIntVec<8>&);
  template bool operator<(const FFIntVec<8>&, const FFIntVec<8>&);
  template bool operator<=(const FFIntVec<8>&, const FFIntVec<8>&);
  template FFIntVec<8> pow(const FFIntVec<8>&, const FFIntVec<8>&);
  template std::ostream& operator<<(std::ostream&, const FFIntVec<8>&);
  // 16
  template class FFIntVec<16>;
  template FFIntVec<16> operator-(const FFIntVec<16>&);
  template FFIntVec<16> operator+(const FFIntVec<16>&);
  template FFIntVec<16> operator-(const FFIntVec<16>&, const FFIntVec<16>&);
  template FFIntVec<16> operator+(const FFIntVec<16>&, const FFIntVec<16>&);
  template FFIntVec<16> operator*(const FFIntVec<16>&, const FFIntVec<16>&);
  template FFIntVec<16> operator/(const FFIntVec<16>&, const FFIntVec<16>&);
  template FFIntVec<16> operator+(const FFIntVec<16>&, const FFInt&);
  template FFIntVec<16> operator-(const FFIntVec<16>&, const FFInt&);
  template FFIntVec<16> operator*(const FFIntVec<16>&, const FFInt&);
  template FFIntVec<16> operator/(const FFIntVec<16>&, const FFInt&);
  template FFIntVec<16> operator+(const FFInt&, const FFIntVec<16>&);
  template FFIntVec<16> operator-(const FFInt&, const FFIntVec<16>&);
  template FFIntVec<16> operator*(const FFInt&, const FFIntVec<16>&);
  template FFIntVec<16> operator/(const FFInt&, const FFIntVec<16>&);
  template bool operator==(const FFIntVec<16>&, const FFIntVec<16>&);
  template bool operator!=(const FFIntVec<16>&, const FFIntVec<16>&);
  template bool operator>(const FFIntVec<16>&, const FFIntVec<16>&);
  template bool operator>=(const FFIntVec<16>&, const FFIntVec<16>&);
  template bool operator<(const FFIntVec<16>&, const FFIntVec<16>&);
  template bool operator<=(const FFIntVec<16>&, const FFIntVec<16>&);
  template FFIntVec<16> pow(const FFIntVec<16>&, const FFIntVec<16>&);
  template std::ostream& operator<<(std::ostream&, const FFIntVec<16>&);
  // 32
  template class FFIntVec<32>;
  template FFIntVec<32> operator-(const FFIntVec<32>&);
  template FFIntVec<32> operator+(const FFIntVec<32>&);
  template FFIntVec<32> operator-(const FFIntVec<32>&, const FFIntVec<32>&);
  template FFIntVec<32> operator+(const FFIntVec<32>&, const FFIntVec<32>&);
  template FFIntVec<32> operator*(const FFIntVec<32>&, const FFIntVec<32>&);
  template FFIntVec<32> operator/(const FFIntVec<32>&, const FFIntVec<32>&);
  template FFIntVec<32> operator+(const FFIntVec<32>&, const FFInt&);
  template FFIntVec<32> operator-(const FFIntVec<32>&, const FFInt&);
  template FFIntVec<32> operator*(const FFIntVec<32>&, const FFInt&);
  template FFIntVec<32> operator/(const FFIntVec<32>&, const FFInt&);
  template FFIntVec<32> operator+(const FFInt&, const FFIntVec<32>&);
  template FFIntVec<32> operator-(const FFInt&, const FFIntVec<32>&);
  template FFIntVec<32> operator*(const FFInt&, const FFIntVec<32>&);
  template FFIntVec<32> operator/(const FFInt&, const FFIntVec<32>&);
  template bool operator==(const FFIntVec<32>&, const FFIntVec<32>&);
  template bool operator!=(const FFIntVec<32>&, const FFIntVec<32>&);
  template bool operator>(const FFIntVec<32>&, const FFIntVec<32>&);
  template bool operator>=(const FFIntVec<32>&, const FFIntVec<32>&);
  template bool operator<(const FFIntVec<32>&, const FFIntVec<32>&);
  template bool operator<=(const FFIntVec<32>&, const FFIntVec<32>&);
  template FFIntVec<32> pow(const FFIntVec<32>&, const FFIntVec<32>&);
  template std::ostream& operator<<(std::ostream&, const FFIntVec<32>&);
  // 64
  template class FFIntVec<64>;
  template FFIntVec<64> operator-(const FFIntVec<64>&);
  template FFIntVec<64> operator+(const FFIntVec<64>&);
  template FFIntVec<64> operator-(const FFIntVec<64>&, const FFIntVec<64>&);
  template FFIntVec<64> operator+(const FFIntVec<64>&, const FFIntVec<64>&);
  template FFIntVec<64> operator*(const FFIntVec<64>&, const FFIntVec<64>&);
  template FFIntVec<64> operator/(const FFIntVec<64>&, const FFIntVec<64>&);
  template FFIntVec<64> operator+(const FFIntVec<64>&, const FFInt&);
  template FFIntVec<64> operator-(const FFIntVec<64>&, const FFInt&);
  template FFIntVec<64> operator*(const FFIntVec<64>&, const FFInt&);
  template FFIntVec<64> operator/(const FFIntVec<64>&, const FFInt&);
  template FFIntVec<64> operator+(const FFInt&, const FFIntVec<64>&);
  template FFIntVec<64> operator-(const FFInt&, const FFIntVec<64>&);
  template FFIntVec<64> operator*(const FFInt&, const FFIntVec<64>&);
  template FFIntVec<64> operator/(const FFInt&, const FFIntVec<64>&);
  template bool operator==(const FFIntVec<64>&, const FFIntVec<64>&);
  template bool operator!=(const FFIntVec<64>&, const FFIntVec<64>&);
  template bool operator>(const FFIntVec<64>&, const FFIntVec<64>&);
  template bool operator>=(const FFIntVec<64>&, const FFIntVec<64>&);
  template bool operator<(const FFIntVec<64>&, const FFIntVec<64>&);
  template bool operator<=(const FFIntVec<64>&, const FFIntVec<64>&);
  template FFIntVec<64> pow(const FFIntVec<64>&, const FFIntVec<64>&);
  template std::ostream& operator<<(std::ostream&, const FFIntVec<64>&);
  // 128
  template class FFIntVec<128>;
  template FFIntVec<128> operator-(const FFIntVec<128>&);
  template FFIntVec<128> operator+(const FFIntVec<128>&);
  template FFIntVec<128> operator-(const FFIntVec<128>&, const FFIntVec<128>&);
  template FFIntVec<128> operator+(const FFIntVec<128>&, const FFIntVec<128>&);
  template FFIntVec<128> operator*(const FFIntVec<128>&, const FFIntVec<128>&);
  template FFIntVec<128> operator/(const FFIntVec<128>&, const FFIntVec<128>&);
  template FFIntVec<128> operator+(const FFIntVec<128>&, const FFInt&);
  template FFIntVec<128> operator-(const FFIntVec<128>&, const FFInt&);
  template FFIntVec<128> operator*(const FFIntVec<128>&, const FFInt&);
  template FFIntVec<128> operator/(const FFIntVec<128>&, const FFInt&);
  template FFIntVec<128> operator+(const FFInt&, const FFIntVec<128>&);
  template FFIntVec<128> operator-(const FFInt&, const FFIntVec<128>&);
  template FFIntVec<128> operator*(const FFInt&, const FFIntVec<128>&);
  template FFIntVec<128> operator/(const FFInt&, const FFIntVec<128>&);
  template bool operator==(const FFIntVec<128>&, const FFIntVec<128>&);
  template bool operator!=(const FFIntVec<128>&, const FFIntVec<128>&);
  template bool operator>(const FFIntVec<128>&, const FFIntVec<128>&);
  template bool operator>=(const FFIntVec<128>&, const FFIntVec<128>&);
  template bool operator<(const FFIntVec<128>&, const FFIntVec<128>&);
  template bool operator<=(const FFIntVec<128>&, const FFIntVec<128>&);
  template FFIntVec<128> pow(const FFIntVec<128>&, const FFIntVec<128>&);
  template std::ostream& operator<<(std::ostream&, const FFIntVec<128>&);
  // 256
  /*template class FFIntVec<256>;
  template FFIntVec<256> operator-(const FFIntVec<256>&);
  template FFIntVec<256> operator+(const FFIntVec<256>&);
  template FFIntVec<256> operator-(const FFIntVec<256>&, const FFIntVec<256>&);
  template FFIntVec<256> operator+(const FFIntVec<256>&, const FFIntVec<256>&);
  template FFIntVec<256> operator*(const FFIntVec<256>&, const FFIntVec<256>&);
  template FFIntVec<256> operator/(const FFIntVec<256>&, const FFIntVec<256>&);
  template FFIntVec<256> operator+(const FFIntVec<256>&, const FFInt&);
  template FFIntVec<256> operator-(const FFIntVec<256>&, const FFInt&);
  template FFIntVec<256> operator*(const FFIntVec<256>&, const FFInt&);
  template FFIntVec<256> operator/(const FFIntVec<256>&, const FFInt&);
  template FFIntVec<256> operator+(const FFInt&, const FFIntVec<256>&);
  template FFIntVec<256> operator-(const FFInt&, const FFIntVec<256>&);
  template FFIntVec<256> operator*(const FFInt&, const FFIntVec<256>&);
  template FFIntVec<256> operator/(const FFInt&, const FFIntVec<256>&);
  template bool operator==(const FFIntVec<256>&, const FFIntVec<256>&);
  template bool operator!=(const FFIntVec<256>&, const FFIntVec<256>&);
  template bool operator>(const FFIntVec<256>&, const FFIntVec<256>&);
  template bool operator>=(const FFIntVec<256>&, const FFIntVec<256>&);
  template bool operator<(const FFIntVec<256>&, const FFIntVec<256>&);
  template bool operator<=(const FFIntVec<256>&, const FFIntVec<256>&);
  template FFIntVec<256> pow(const FFIntVec<256>&, const FFIntVec<256>&);
  template std::ostream& operator<<(std::ostream&, const FFIntVec<256>&);*/
}
