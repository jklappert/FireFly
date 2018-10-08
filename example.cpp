#include <iostream>
#include "PolyReconst.hpp"
#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "Logger.hpp"
#include <algorithm>

int main() {
  uint n = 1;
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
      std::vector<firefly::FFInt> t_yis;
      std::vector<firefly::FFInt> yis;
      yis.reserve(n - 1);
      t_yis.reserve(n - 1);
      firefly::FFInt t;
      firefly::RatReconst& rec = map.at(m);
      prime = firefly::primes()[rec.prime_number];
      firefly::FFInt::p = prime;
      //std::cout << "Reconstruction: " << m << "\n";
      int count = 0;
      int kk = 0;
      uint primes_used = 0;

      while (!rec.done) {
        if (yis.size() > 0) {
          for (uint j = 2; j <= n; j++) {
            t_yis[j - 2] = yis[j - 2];
          }
        }

        if (primes_used != rec.prime_number) {
          std::cout << "Set new prime. Iterations for last prime: " << kk << ".\n";
          primes_used = std::max(primes_used, rec.prime_number);
          prime = firefly::primes()[rec.prime_number];
          firefly::FFInt::p = prime;

          if (n > 1)
            yis = std::vector<firefly::FFInt> (rec.curr_zi_order.begin(), rec.curr_zi_order.end() - 1);

          kk = 0;

          for (uint j = 2; j <= n; j++) {
            t_yis[j - 2] = t * yis[j - 2] + firefly::RatReconst::shift[j - 1];
          }
        }

        t = firefly::FFInt(std::rand() % (prime - 1)) + firefly::FFInt(1);

        if (t_yis.size() == 0) {
          for (uint j = 2; j <= n; j++) {
            yis = std::vector<firefly::FFInt> (rec.curr_zi_order.begin(), rec.curr_zi_order.end() - 1);
            t_yis.emplace_back(t * yis[j - 2] + firefly::RatReconst::shift[j - 1]);
          }
        }

        if (rec.zi >= 2) {
          yis = std::vector<firefly::FFInt> (rec.curr_zi_order.begin(), rec.curr_zi_order.end() - 1);

          for (uint j = 2; j <= n; j++) {
            t_yis[j - 2] = t * yis[j - 2] + firefly::RatReconst::shift[j - 1];
          }
        } else {
          for (uint j = 2; j <= n; j++) {
            t_yis[j - 2] = t * yis[j - 2] + firefly::RatReconst::shift[j - 1];
          }
        }

        firefly::FFInt z1 = t + firefly::RatReconst::shift[0];

        mpz_class test;
        test = "123456789109898799879870980";
        mpz_class test2;
        test2 = "123456789109898799879";
        firefly::FFInt a7(test);
        firefly::FFInt a8(test2);
        firefly::FFInt a1(1);
        firefly::FFInt a2(18);
        firefly::FFInt a3(25);
        firefly::FFInt a4(10);
        firefly::FFInt a5(2);
        firefly::FFInt a6(3);
        // example for n = 5
        /*firefly::FFInt num = ((z1*z1 -firefly::FFInt(5)*z1+firefly::FFInt(6))*t_yis[0]
          +firefly::FFInt(2)*z1*z1-firefly::FFInt(10)*z1+firefly::FFInt(12));
        firefly::FFInt den = (((firefly::FFInt(2)*z1-firefly::FFInt(8))*t_yis[1])*t_yis[0]*t_yis[0]
          +((-firefly::FFInt(4)*z1+firefly::FFInt(16))*t_yis[1])*t_yis[0]+(firefly::FFInt(2)*z1-firefly::FFInt(8))*t_yis[1]);*/

        // example for n = 1
        firefly::FFInt num = (firefly::FFInt(576)*z1.pow(firefly::FFInt(12)) - firefly::FFInt(35145)*z1.pow(firefly::FFInt(11))
          +firefly::FFInt(946716)*z1.pow(firefly::FFInt(10))-firefly::FFInt(14842335)*z1.pow(firefly::FFInt(9))
          +firefly::FFInt(150236238)*z1.pow(firefly::FFInt(8))-firefly::FFInt(1028892363)*z1.pow(firefly::FFInt(7))
          +firefly::FFInt(4853217576)*z1.pow(firefly::FFInt(6))-firefly::FFInt(15724949577)*z1.pow(firefly::FFInt(5))
          +firefly::FFInt(34208917206)*z1.pow(firefly::FFInt(4))-firefly::FFInt(47506433412)*z1.pow(firefly::FFInt(3))
          +firefly::FFInt(37933483608)*z1.pow(firefly::FFInt(2))-firefly::FFInt(13296184128)*z1+firefly::FFInt(71850240));
        firefly::FFInt den = (firefly::FFInt(16)*z1.pow(firefly::FFInt(12))-firefly::FFInt(960)*z1.pow(firefly::FFInt(11))
          +firefly::FFInt(25456)*z1.pow(firefly::FFInt(10))-firefly::FFInt(393440)*z1.pow(firefly::FFInt(9))
          +firefly::FFInt(3934768)*z1.pow(firefly::FFInt(8))-firefly::FFInt(26714240)*z1.pow(firefly::FFInt(7))
          +firefly::FFInt(125545488)*z1.pow(firefly::FFInt(6))-firefly::FFInt(408157280)*z1.pow(firefly::FFInt(5))
          +firefly::FFInt(899198016)*z1.pow(firefly::FFInt(4))-firefly::FFInt(1278172800)*z1.pow(firefly::FFInt(3))
          +firefly::FFInt(1055033856)*z1.pow(firefly::FFInt(2))-firefly::FFInt(383201280)*z1);

        kk++;
        count++;
        std::vector<uint> tmp_vec;

        if(n > 1)
          tmp_vec = std::vector<uint>(rec.curr_zi_order.begin(), rec.curr_zi_order.end() - 1);
        rec.feed(t, num/den, tmp_vec);
      }
      std::cout << "Total numerical runs: " << count << ", primes used: " << primes_used + 1 << ".\n";
      std::cout << rec.get_result();
    }
  } catch (std::exception& e) {
    ERROR_MSG(e.what());
  }

  return 0;
}

