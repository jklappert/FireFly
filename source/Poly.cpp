//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2019  Jonas Klappert, Sven Yannick Klein and Fabian Lange
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
//==================================================================================

#include "Poly.hpp"

#include <math.h>
#include <random>

namespace firefly {

  Poly::Poly() {
    coeff = std::vector<FFInt>();
  }

  Poly::Poly(std::vector<FFInt>& coeff_vector) {
    coeff = coeff_vector;
  }

  Poly::Poly(const Poly& old_poly) {
    coeff = old_poly.coeff;
  }

  Poly::~Poly() {
    std::vector<FFInt>().swap(coeff);
  }

  size_t Poly::get_deg() const {
    for (size_t i = coeff.size() - 1; i >= 1; i--) {
      if (coeff.at(i) != 0) {
        return i;
      }
    }

    return 0;
  }

  void Poly::shrink_to_fit() {
    while (coeff.back() == FFInt(0)) {
      coeff.pop_back();
    };

    coeff.shrink_to_fit();
  }

  void Poly::rev() {
    shrink_to_fit();
    std::reverse(std::begin(coeff), std::end(coeff));
  }

  Poly& Poly::operator-=(const Poly& a) {
    if (a.get_deg() > get_deg()) {
      coeff.reserve(a.get_deg() + 1);
      coeff.resize(a.get_deg() + 1);
    };

    for (size_t i = 0; i <= a.get_deg(); i++) {
      if (i <= get_deg() && i <= a.get_deg()) {
        coeff.at(i) = coeff.at(i) - a.coeff.at(i);
      } else if (i > get_deg() && i <= a.get_deg()) {
        coeff.emplace(coeff.begin() + i, - a.coeff.at(i));
      };
    }

    return *this;
  }

  Poly& Poly::operator+=(const Poly& a) {
    if (a.get_deg() > get_deg()) {
      coeff.reserve(a.get_deg() + 1);
      coeff.resize(a.get_deg() + 1);
    };

    for (size_t i = 0; i <= a.get_deg(); i++) {
      if (i <= get_deg() && i <= a.get_deg()) {
        coeff.at(i) = coeff.at(i) + a.coeff.at(i);
      } else if (i > get_deg() && i <= a.get_deg()) {
        coeff.emplace(coeff.begin() + i, a.coeff.at(i));
      };
    }

    return *this;
  }

  Poly& Poly::operator*=(const FFInt& a) {
    for (size_t i = 0; i <= get_deg(); i++) {
      coeff.at(i) *= a;
    };

    return *this;
  }

  Poly& Poly::operator/=(const FFInt& a) {
    for (size_t i = 0; i <= get_deg(); i++) {
      coeff.at(i) /= a;
    };

    return *this;
  }

  Poly operator*(const Poly& a, const FFInt& b) {
    Poly returner(a);
    returner *= b;
    return returner;
  }

  Poly operator/(const Poly& a, const FFInt& b) {
    Poly returner(a);
    returner /= b;
    return returner;
  }

  Poly operator+(const Poly& a, const Poly& b) {
    Poly returner(a);
    returner += b;
    return returner;
  }

  Poly operator-(const Poly& a, const Poly& b) {
    Poly returner(a);
    returner -= b;
    return returner;
  }

  Poly& Poly::operator*=(const Poly& a) {
    std::vector<FFInt> tmp_coeff {};
    tmp_coeff.reserve(((get_deg() + 1) * (a.get_deg() + 1)));
    tmp_coeff.resize(((get_deg() + 1) * (a.get_deg() + 1)));;

    for (size_t i = 0; i < tmp_coeff.size(); i++) {
      FFInt this_coeff(0);

      for (size_t j = 0; j <= i; j++) {
        if (coeff.size() > j && a.coeff.size() > i - j) {
          this_coeff += coeff.at(j) * a.coeff.at(i - j);
        }

        tmp_coeff.at(i) = this_coeff;
      };
    }

    coeff.swap(tmp_coeff);
    return *this;
  }

  Poly operator*(const Poly& a, const Poly& b) {
    Poly returner(a);
    returner *= b;
    return returner;
  }

  std::pair<Poly, Poly> fast_euclidean_division(const Poly& a, const Poly& z) {
    if (a.get_deg() < z.get_deg()) {
      std::vector<FFInt> a0 {FFInt(0)};
      return std::make_pair(Poly(a0), a);
    }

    Poly b(z);
    b.shrink_to_fit();
    b /= b.coeff.back();
    size_t m = a.get_deg() - b.get_deg();
    size_t k = ((size_t)(ceil(log2((double)m + 1))));
    std::vector<FFInt> a1 {FFInt(1)};
    Poly gi(a1);
    Poly f(b);
    f.rev();

    for (size_t i = 1; i <= k; i++) {
      Poly q(gi);
      q *= FFInt(2);
      Poly tmp(gi);
      tmp *= gi;
      tmp *= f;
      gi = q - tmp;

      while (gi.get_deg() >= ((size_t) exp2((double) i))) {
        gi.coeff.pop_back();
      }
    };

    Poly s(a);

    s.rev();

    s *= gi;

    while (s.get_deg() >= m + 1) {
      s.coeff.pop_back();
    }

    size_t exp = m - s.get_deg();
    s.rev();

    for (size_t i = 1; i <= exp; i++) {
      s.coeff.insert(s.coeff.begin(), FFInt(0));
    };

    Poly r(a);

    r -= b * s;

    s /= z.coeff.at(z.get_deg());

    return std::make_pair(s, r);
  }

  Poly operator/(const Poly& a, const Poly& b) {
    return (fast_euclidean_division(a, b).first);
  }

  Poly operator%(const Poly& a, const Poly& b) {
    return (fast_euclidean_division(a, b).second);
  }

  Poly gcd(const Poly& c, const Poly& d) {
    Poly a(c);
    Poly b(d);

    while (b.get_deg() != 0 || b.coeff.at(0) != FFInt(0)) {
      Poly h;
      h = a % b;
      a = b;
      b = h;
    }

    return a;
  }

  std::vector<FFInt> Poly::roots() {
    Poly a;
    a.coeff = coeff;
    a.shrink_to_fit();
    a /= a.coeff.back();
    std::vector<FFInt> roots {};

    while (a.coeff.at(0) == FFInt(0)) {
      a.coeff.erase(a.coeff.begin());
      roots.emplace_back(FFInt(0));
    };

    std::vector<Poly> queue {a};

    while (!queue.empty()) {
      while (queue.back().get_deg() > 1) {
        std::random_device rd;
        std::mt19937_64 eng(rd());//TODO Mersenne-Twister is too slow. Replace with pcg32
        std::uniform_int_distribution<unsigned long long> distr;
        FFInt Delta(distr(eng));
        size_t s = ((FFInt(0).p - 1) / 2);
        std::vector<FFInt> master_vec {FFInt(1)};
        Poly master(master_vec);
        std::vector<FFInt>().swap(master_vec);
        master.coeff.reserve(queue.back().get_deg() + 2);
        master.coeff.resize(queue.back().get_deg() + 2);
        std::vector<FFInt> multi_vec {FFInt(Delta), FFInt(1)};
        Poly multi(multi_vec);
        std::vector<FFInt>().swap(multi_vec);

        for (size_t i = 1; i <= s; i++) {
          master *= multi;
          master = (master % queue.back());
        };

        master.coeff.at(0) -= FFInt(1);

        Poly h;

        h = gcd(queue.back(), master);

        Poly f;

        f = queue.back() / h;

        queue.back() = h;

        queue.back().shrink_to_fit();

        queue.emplace_back(f);

        queue.back().shrink_to_fit();
      };

      if (queue.back().get_deg() == 0) {
        queue.pop_back();
      }

      while (!queue.empty() && queue.back().get_deg() == 1) {
        roots.emplace_back(- queue.back().coeff.at(0) / queue.back().coeff.at(1));
        queue.pop_back();
      }
    }

    return roots;
  }

  std::ostream& operator<<(std::ostream& out, const Poly& a) {
    if (a.coeff.empty()) {
      return out << "0";
    }

    if (a.get_deg() == 0) {
      return out << a.coeff.at(0);
    }

    if (a.get_deg() == 1 && a.coeff.at(0) != FFInt(0)) {
      return out << a.coeff.at(1) << "x+" << a.coeff.at(0);
    };

    if (a.get_deg() == 1 && a.coeff.at(0) == FFInt(0)) {
      return out << a.coeff.at(1) << "x";
    };

    out << a.coeff.at(a.get_deg()) << "x^" << a.get_deg();

    for (size_t i = a.get_deg() - 1; i > 1; i--) {
      if (a.coeff.at(i) != 0) {
        out << "+" << a.coeff.at(i) << "x^" << i;
      };
    };

    if (a.coeff.at(1) != 0) {
      out << "+" << a.coeff.at(1) << "x";
    }

    if (a.coeff.at(0) != 0) {
      out << "+" << a.coeff.at(0);
    }

    return out;
  }
}
