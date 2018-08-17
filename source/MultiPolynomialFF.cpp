#include "MultiPolynomialFF.hpp"

namespace firefly{

  MultiPolynomialFF::MultiPolynomialFF() {}

  MultiPolynomialFF::MultiPolynomialFF(uint n_, const std::map<uint, MultiPolynomialFF>& coef_) : n(n_), coef(coef_) {
    deg = coef.rend()->first;
  }

  MultiPolynomialFF MultiPolynomialFF::operator+(const MultiPolynomialFF& other) {
    if (n == other.n && n > 1) {
      std::map<uint, MultiPolynomialFF> new_coef;

      if (deg >= other.deg) {
        new_coef = coef;
        for (auto element : other.coef) {
          auto insertion = new_coef.insert(element);
          if (!insertion.second) insertion.first->second = insertion.first->second + element.second;
        }
      } else {
        new_coef = other.coef;
        for (auto element : coef) {
          auto insertion = new_coef.insert(element);
          if (!insertion.second) insertion.first->second = insertion.first->second + element.second;
        }
      }

      return MultiPolynomialFF(n, std::move(new_coef));
    } else throw std::runtime_error("Adding multipolynomials of different variables!");
  }

  MultiPolynomialFF MultiPolynomialFF::operator-(const MultiPolynomialFF& other) {
    if (n == other.n && n > 1) {
      std::map<uint, MultiPolynomialFF> new_coef = coef;
      MultiPolynomialFF empty_coef(n - 1, std::map<uint, MultiPolynomialFF>());

      for (auto element : other.coef) {
        auto insertion = new_coef.find(element.first);
        if (insertion != new_coef.end()) {
          insertion->second = insertion->second - element.second;
        } else {
          new_coef.insert(std::make_pair(element.first, empty_coef - element.second));
        }
      }

      return MultiPolynomialFF(n, std::move(new_coef));
    } else throw std::runtime_error("Subtracting multipolynomials of different variables!");
  }

  MultiPolynomialFF MultiPolynomialFF::operator*(const FFInt& factor) {
    std::map<uint, MultiPolynomialFF> new_coef = coef;

    for (auto element : new_coef) {
      element.second = element.second * factor;
    }

    return MultiPolynomialFF(n, std::move(new_coef));
  }

}
