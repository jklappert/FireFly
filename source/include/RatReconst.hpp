#pragma once

#include <cstdint>
#include <vector>
#include "FFInt.hpp"
#include "Polynomial.hpp"

namespace firefly{

   class RatReconst{
   public:
      /**
       * 	A constructor
       * 	@param n_ the number of parameters
       */
      RatReconst(int n_);
      /**
       * 	Calls the reconstruction algorithm
       * 	@returns a vector with a pair of integers corrisponding to the
       * 	coefficient of the polynom. The vector is ordered in an anscending
       * 	way such that the first coefficient corresponds to to z^0,...
       */
      std::vector<FFInt> reconst();
      /**
       * 	Constructs the canonical form of the rational function recursevly
       */
      void constrCanonical();
      /**
       * 	Returns the canonical form of the rational function, where
       * 	the first entry corresponds to the polyomial in the numerator
       * 	and the second entry is the polynomial in the denominator
       */
      std::pair<Polynomial, Polynomial> canonical;
   private:
      int n; /**< The number of parameters */
      const uint64_t prime;
      FFInt compAi(int i, int ip, const FFInt& num);
      Polynomial iterateCanonical(uint i);
      FFInt iterateCanonicalNum(uint i, uint ip, const FFInt& num);
      FFInt compFyi(int i, FFInt& y);
      std::vector<FFInt> ai {};
      std::vector<FFInt> yi {};
      std::vector<FFInt> fyi {};
      /**
       * 	A numerical black box function which provides the reconstruction
       * 	algorithm with the finite field member f(y)
       * 	@param p the prime number which defines the finite field
       * 	@param y a member of the finite field which should be evaluated in f(y)
       */
      FFInt num(uint64_t p, const FFInt& y);
   };
}