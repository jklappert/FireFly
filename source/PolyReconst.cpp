#include <cstdlib>
#include "PolyReconst.hpp"
#include "ReconstHelper.hpp"
#include "Logger.hpp"

namespace firefly {

PolyReconst::PolyReconst(int n_) : n(n_), prime(primes().at(0)) {
   ai.reserve(5000);
   yi.reserve(5000);
}

std::vector<FFInt> PolyReconst::reconst(){
   uint maxDegree = 5000;
   const int breakCondition = 3;

   yi.emplace_back(FFInt(std::rand() % prime, prime));
   ai.emplace_back(num(prime, yi.back()));

   for(uint i = 1; i < maxDegree; i++){
      yi.emplace_back(FFInt(std::rand() % prime, prime));
      FFInt fyi;

      bool spuriousPole = true;
      while(spuriousPole){
	 try{
	    fyi = num(prime, yi.back());
	    spuriousPole = false;
	 } catch(const std::exception&){
	    yi.pop_back();
	    yi.emplace_back(FFInt(std::rand() % prime, prime));
	 }
      }

      spuriousPole = true;
      while(spuriousPole){
	 try{
	    ai.emplace_back(compAi(i, i, fyi));
	    spuriousPole = false;
	 } catch (const std::exception&){
	    yi.pop_back();
	    yi.emplace_back(FFInt(std::rand() % prime, prime));
	 }
      }

      if(ai.at(i).n == 0) {
	 if(i > breakCondition){
	    bool nonZero = false;
	    for(uint j = ai.size(); j > ai.size() - breakCondition; j--){
	       if(ai.at(j - 1).n != 0){
		  nonZero = true;
		  break;
	       }
	    }
	    if(!nonZero) break;
	 }
      }
      if(i == maxDegree - 1){
	 maxDegree += 5000;
	 yi.resize(maxDegree);
	 yi.resize(maxDegree);
      }
   }

   for(int i = 0; i < breakCondition; i++){
      yi.pop_back();
      ai.pop_back();
   }
   return ai;
}

FFInt PolyReconst::compAi(int i, int ip, const FFInt& num){
   if(ip == 0){
      return num;
   } else{
      if(yi.at(i) == yi.at(ip - 1)) throw std::runtime_error("Divide by 0 error!");
      return (compAi(i, ip - 1, num) - ai.at(ip - 1))/(yi.at(i) - yi.at(ip - 1));
   }
}

void PolyReconst::constrCanonical(){
   if(ai.size() == 0){
      INFO_MSG("Polynomial not yet reconstructed or 0.");
   } else if(ai.size() == 1){
      std::vector<FFInt> coef{ai.at(0)};
      Polynomial poly(coef);
      canonical = poly;
   } else{
      std::vector<FFInt> coef{ai.at(0)};
      Polynomial poly(coef);
      canonical = poly + iterateCanonical(1);
   }
}

Polynomial PolyReconst::iterateCanonical(uint i){
   if(i < ai.size() - 1){
      std::vector<FFInt> coef1{ (FFInt(0, prime) - yi.at(i - 1)) * ai.at(i), ai.at(i) };
      std::vector<FFInt> coef2{ FFInt(0, prime) - yi.at(i - 1), FFInt(1, prime) };
      Polynomial poly1(coef1);
      Polynomial poly2(coef2);
      return poly1 + poly2 * iterateCanonical(i + 1);
   } else{
      std::vector<FFInt> coef1{ (FFInt(0, prime) - yi.at(i - 1)) * ai.at(i), ai.at(i) };
      Polynomial poly1(coef1);
      return poly1;
   }
}

FFInt PolyReconst::num(uint64_t p, const FFInt& y){
   FFInt a0 (3, p);
   FFInt a1 (6, p);
   FFInt a2 (18, p);
   FFInt a3 (25, p);
   FFInt a4 (30, p);
   FFInt a5 (2, p);
   FFInt a6 (7, p);
   FFInt a7 (100, p);
   FFInt a8 (13, p);
   FFInt exp2 (2, p);
   FFInt exp3 (3, p);
   FFInt exp4 (4, p);
   FFInt exp5 (5, p);
   FFInt exp6 (6, p);
   FFInt exp7 (7, p);
   FFInt exp8 (4000, p);

   return a0 + a1*y + a2*y.pow(exp2) + a3*y.pow(exp3) + a4*y.pow(exp4)
      + a5*y.pow(exp5) + a6*y.pow(exp6) + a7*y.pow(exp7) + a8*y.pow(exp8);
}

}
