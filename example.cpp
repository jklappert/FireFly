#include <iostream>
#include "PolyReconst.hpp"
#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "Logger.hpp"

int main() {
  uint n = 3;
  uint64_t prime = firefly::primes()[0];
  firefly::FFInt::p = prime;
  firefly::RatReconst rec_1(n);
  firefly::PolyReconst rec_2(n);
  std::unordered_map<int, firefly::RatReconst> map;

  for (int m = 0; m < 1; m++) {
    firefly::RatReconst tmp_rec(n);
    map.emplace(m, std::move(tmp_rec));
  }

  try {
    /*std::vector<firefly::FFInt> t_yis;
    std::vector<firefly::FFInt> yis;
    yis.reserve(n - 1);
    t_yis.reserve(n - 1);
    firefly::FFInt t;
    std::vector<firefly::FFInt> yis_2;
    yis_2.reserve(n);
    yis_2.emplace_back(firefly::FFInt(std::rand() % prime));
    yis_2.emplace_back(firefly::FFInt(std::rand() % prime));
    yis_2.emplace_back(firefly::FFInt(std::rand() % prime));
    yis_2.emplace_back(firefly::FFInt(std::rand() % prime));
    while(!rec_2.done){
      if (rec.new_prime) {
        prime = firefly::primes()[i];
        firefly::FFInt::p = prime;
        i++;
      }
      yis_2[rec_2.next_zi - 1] = firefly::FFInt(std::rand() % prime);

      mpz_class test;
      test = "1234567891098987998798709805302432022989874343098";
      test = test % prime;
      firefly::FFInt a7(std::stoull(test.get_str()));
      firefly::FFInt a1(1);
      firefly::FFInt a2(18);
      firefly::FFInt a3(25);
      firefly::FFInt a4(10);
      firefly::FFInt a5(2);
      firefly::FFInt a6(7);
      firefly::FFInt num = (a6 + a3 * yis_2[0] + yis_2[0]*a3*yis_2[2] + yis_2[1]*yis_2[1] + yis_2[3]);
      rec_2.feed(yis_2, num);
    }*/
    for (int m = 0; m < 1; m++) {
//      if(m >= 1) std::cout << map.at(m - 1).get_result() << "\n";
/*      if (m == 0) {
        firefly::RatReconst::shift[0] = 1;
      } else if (m == 1) {
        firefly::RatReconst::shift[0] = 2;
      }*/
      int i = 1;
      prime = firefly::primes()[0];
      firefly::FFInt::p = prime;
      std::vector<firefly::FFInt> t_yis;
      std::vector<firefly::FFInt> yis;
      yis.reserve(n - 1);
      t_yis.reserve(n - 1);
      firefly::FFInt t;
      firefly::RatReconst& rec = map.at(m);
      std::cout << "Reconstruction: " << m << " shift " << firefly::RatReconst::shift[3] << " " << firefly::RatReconst::shift[2] << "\n";
      int kk = 0;
      while (!rec.done) {
        if (yis.size() > 0) {
          for (uint j = 2; j <= n; j++) {
            t_yis[j - 2] = yis[j - 2];
          }
        }

        if (rec.new_prime) {
          prime = firefly::primes()[i];
          firefly::FFInt::p = prime;
          i++;
        }

        t = firefly::FFInt(std::rand() % prime);

        if (t_yis.size() == 0) {
          for (uint j = 2; j <= n; j++) {
            yis.emplace_back(firefly::FFInt(std::rand() % prime));
            t_yis.emplace_back(t * yis[j - 2] + firefly::RatReconst::shift[j - 1]);
          }
        }

        if (rec.zi >= 2) {
          for (uint j = 2; j <= n; j++) {
            if (j <= rec.zi) {
              yis[j - 2] = firefly::FFInt(std::rand() % prime);
            }

            t_yis[j - 2] = t * yis[j - 2] + firefly::RatReconst::shift[j - 1];
          }
        } else {
          for (uint j = 2; j <= n; j++) {
            t_yis[j - 2] = t * yis[j - 2] + firefly::RatReconst::shift[j - 1];
          }
        }

        firefly::FFInt z1 = t + firefly::RatReconst::shift[0];

        mpz_class test;
        test = "1234567891098987998798709805302432022989874343098";
        test = test % prime;
        firefly::FFInt a7(std::stoull(test.get_str()));
        firefly::FFInt a1(1);
        firefly::FFInt a2(18);
        firefly::FFInt a3(25);
        firefly::FFInt a4(10);
        firefly::FFInt a5(1);
        firefly::FFInt a6(3);

        firefly::FFInt num = a1;
        firefly::FFInt den = a1;
        for (uint i = 1; i < 5; i++) {
          num += z1.pow(firefly::FFInt(i));

          for (uint j = 0; j < n - 1; j++) {
            num += t_yis[j].pow(firefly::FFInt(i));
            den += t_yis[j].pow(firefly::FFInt(i));
          }
        }

        rec.feed(t, yis, num / den);
        kk++;
      }
      std::cout << "times iterated " << kk << " new primes " << i << "\n";
      std::cout << rec.get_result() << "\n";
    }

    /*std::cout << map.at(0).get_result();

    if (n == 2) {
      std::vector<std::string> symbols = {"z1", "z2"};
      std::cout << rec_1.get_result().string(symbols) << "\n";
    }

    //std::cout << "f(x) = " << rec_2.get_result();*/
  } catch (std::exception& e) {
    ERROR_MSG(e.what());
  }

  return 0;
}
