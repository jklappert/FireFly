#include <iostream>
#include "PolyReconst.hpp"
#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "Logger.hpp"

int main() {
  uint64_t prime = firefly::primes()[0];
  firefly::RatReconst rec_rat(1, prime);
  //firefly::PolyReconst rec_pol(4);

  try {
    int i = 1;
    while (!rec_rat.done) {
      if(rec_rat.new_prime) {
        prime = firefly::primes()[i];
        i++;
      }
      mpz_class test;
      test = "1234567891098987998798709805302432098098743432098";
      test = test % prime;
      firefly::FFInt a7(std::stoull(test.get_str()), prime);
      std::vector<firefly::FFInt> t = {firefly::FFInt(std::rand() % prime, prime)};
      firefly::FFInt a1(1,prime);
      firefly::FFInt a2(18, prime);
      firefly::FFInt a3(25, prime);
      firefly::FFInt a4(300, prime);
      firefly::FFInt a5(2, prime);
      firefly::FFInt a6(7, prime);
      firefly::FFInt num = (a1-a2*t[0].pow(a4))/(a1 - a6/a5*t[0].pow(a6));
      rec_rat.feed(prime, t[0], t, num);
    }

    std::cout << rec_rat.get_result();
    //auto rat_fun = rec_rat.reconst();
    //auto pol_fun = rec_pol.reconst();
    //std::cout << rat_fun;
    //std::cout << "f(x) = " << pol_fun;
  } catch (std::exception &e) {
    ERROR_MSG(e.what());
  }

  return 0;
}

