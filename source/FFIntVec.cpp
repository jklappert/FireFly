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

#include "FFIntVec.hpp"

namespace firefly {

  FFIntVec::FFIntVec(const std::vector<FFInt>& in) {
    vec = in;
  }

  bool operator==(const FFIntVec& a, const FFIntVec& b) {
    return a.vec == b.vec;
  }

  bool operator!=(const FFIntVec& a, const FFIntVec& b) {
    return a.vec != b.vec;
  }

  FFIntVec& FFIntVec::operator+=(const FFIntVec& a) {
    for (size_t i = 0; i != a.vec.size(); ++i) {
      vec[i] += a.vec[i];
    }

    return *this;
  }

  FFIntVec& FFIntVec::operator-=(const FFIntVec& a) {
    for (size_t i = 0; i != a.vec.size(); ++i) {
      vec[i] -= a.vec[i];
    }

    return *this;
  }

  FFIntVec& FFIntVec::operator*=(const FFIntVec& a) {
    for (size_t i = 0; i != a.vec.size(); ++i) {
      vec[i] *= a.vec[i];
    }

    return *this;
  }

  FFIntVec& FFIntVec::operator/=(const FFIntVec& a) {
    for (size_t i = 0; i != a.vec.size(); ++i) {
      vec[i] /= a.vec[i];
    }

    return *this;
  }

  FFIntVec FFIntVec::pow(const FFIntVec& power) const {
    FFIntVec result(vec);
    for (size_t i = 0; i != power.vec.size(); ++i) {
      result.vec[i] = result.vec[i].pow(power.vec[i]);
    }

    return *this;
  }

  FFIntVec operator+(const FFIntVec& a, const FFIntVec& b) {
    FFIntVec result(a);

    for (size_t i = 0; i != a.vec.size(); ++i) {
      result.vec[i] += b.vec[i];
    }

    return result;
  }

  FFIntVec operator-(const FFIntVec& a, const FFIntVec& b) {
    FFIntVec result(a);

    for (size_t i = 0; i != a.vec.size(); ++i) {
      result.vec[i] -= b.vec[i];
    }

    return result;
  }

  FFIntVec operator*(const FFIntVec& a, const FFIntVec& b) {
    FFIntVec result(a);

    for (size_t i = 0; i != a.vec.size(); ++i) {
      result.vec[i] *= b.vec[i];
    }

    return result;
  }

  FFIntVec operator/(const FFIntVec& a, const FFIntVec& b) {
    FFIntVec result(a);

    for (size_t i = 0; i != a.vec.size(); ++i) {
      result.vec[i] /= b.vec[i];
    }

    return result;
  }

  FFIntVec pow(const FFIntVec& a, const FFIntVec& power) {
    FFIntVec result(a);

    for (size_t i = 0; i != power.vec.size(); ++i) {
      result.vec[i] = result.vec[i].pow(power.vec[i]);
    }

    return result;
  }

  std::ostream& operator<<(std::ostream& out, const FFIntVec& ffint_vec) {
    out << "(" << ffint_vec.vec.front();

    for (auto it = ffint_vec.vec.begin() + 1; it != ffint_vec.vec.end(); ++ it) {
      out << ", " << *it;
    }

    out << ")";
    return out;
  }
}
