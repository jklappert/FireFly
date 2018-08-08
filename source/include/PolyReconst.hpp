#pragma once

#include <cstdint>
#include <vector>
#include "FFInt.hpp"

namespace firefly{
   
   class PolyReconst{
   public:
      /**
       * 	A constructor
       * 	@param n_ the number of parameters as an integer
       */
      PolyReconst(int n_);
      /**
       * 	Calls the reconstruction algorithm
       * 	@returns a vector with a pair of integers corrisponding to the 
       * 	coefficient of the polynom. The vector is ordered in an anscending
       * 	way such that the first coefficient corresponds to to z^0,...
       */
      std::vector<FFInt> reconst();
   private:
      FFInt num(uint64_t p, const FFInt& y);
      FFInt compAi(int i, int ip, const FFInt& num);
      
      int n; /**< The number of parameters */
      std::vector<FFInt> ai;
      std::vector<FFInt> yi;
   };
}