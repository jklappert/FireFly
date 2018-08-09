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
      const uint64_t prime; /**< The prime number defining the finite field */
      /**
       * 	Computes the coefficient a_i recursevly using eq. (3.11) of
       * 	arXiv:1608.01902
       * 	@param i the order of a_i
       * 	@param ip recursion order
       * 	@param num f(y_i)
       * 	@returns a_i
       */
      FFInt compAi(int i, int ip, const FFInt& num);
      Polynomial iterateCanonical(uint i);
      /**
       * 	A numerical implementation of Thiele's interpolation formula
       * 	from arXiv:1608.01902 eq. (3.10)
       * 	@param i the order of the highest a_i
       * 	@param ip recusrion order
       * 	@param num the argument of f(z) which is a finite field member
       * 	@returns a finite field member which corresponds to one recusion
       * 	step
       */
      FFInt iterateCanonicalNum(uint i, uint ip, const FFInt& num);
      /**
       * 	Calculates f(y_i) using  Thiele's interpolation formula
       * 	@param i order of the highest coefficient a_i
       * 	@param y y_i
       * 	@returns f(y_i)
       */
      FFInt compFyi(int i, FFInt& y);
      std::vector<FFInt> ai {}; /**< A vector which holds all coefficients a_i */
      std::vector<FFInt> yi {}; /**< A vector which holds all arguments y_i */
      std::vector<FFInt> fyi {}; /**< A vector which holds all function values f(y_i) */
      /**
       * 	A numerical black box function which provides the reconstruction
       * 	algorithm with the finite field member f(y)
       * 	@param p the prime number which defines the finite field
       * 	@param y a member of the finite field which should be evaluated in f(y)
       */
      FFInt num(uint64_t p, const FFInt& y);
   };
}