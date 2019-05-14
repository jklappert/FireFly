#ifndef POLY_H_
#define POLY_H_

#pragma once

#include <string>
#include <vector>
#include "FFInt.hpp"

namespace firefly{

	class Poly {

		public:

			Poly();

			Poly(std::vector<FFInt> &);

			Poly(const Poly &);

			~Poly();

			std::vector<FFInt> coeff;

			size_t get_deg() const;

			void shrink_to_fit();

			void rev();

			std::vector<FFInt> roots();

			Poly& operator=(const Poly&) = default;
	    Poly& operator-=(const Poly&);
	    Poly& operator+=(const Poly&);
	    Poly& operator*=(const FFInt&);
	    Poly& operator/=(const FFInt&);
			Poly& operator*=(const Poly&);

		private:
	};

  Poly operator+(const Poly& a, const Poly& b);
  Poly operator-(const Poly& a, const Poly& b);
  Poly operator*(const Poly& a, const FFInt&);
  Poly operator/(const Poly& a, const FFInt&);
	Poly operator*(const Poly& a, const Poly& b);
	Poly operator/(const Poly& a, const Poly& b);
	Poly operator%(const Poly& a, const Poly& b);
	std::ostream& operator<<(std::ostream& out, const Poly& a);

	std::pair<Poly, Poly> fast_euclidean_division(const Poly& a, const Poly& b);

	Poly gcd(const Poly& a, const Poly& b);
};

#endif /* POLY_H_ */
