#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "Tests.hpp"
#include "Logger.hpp"

#include "utils.hpp"

using namespace firefly;

int main() {
  // Example for the reconstruction of a rational function
  uint n = 4;
  FFInt::set_new_prime(primes()[0]);
  RatReconst rec(n);
  //rec.set_tag("gghh");
  //rec.start_from_saved_file("ff_save/gghh_3.txt");

  std::cout << "--------------------------------------------------------------\n";
  std::cout << "Interpolating rational function\n";
  std::cout << "--------------------------------------------------------------\n";
  // Initialize all values. t_yis are the scaled y-values and yis are the
  // unshifted ones
  std::vector<FFInt> t_yis(n);
  t_yis.reserve(n - 1);
  // t is the scaling variable
  FFInt t;

  // Initialize some counters
  int count = 0;
  int kk = 0;
  uint primes_used = 0;

  // Feed loop
  std::vector<FFInt> shift = rec.get_zi_shift_vec();

  while (!rec.is_done()) {
    // If a new prime is needed, set it, generate new random variables
    // and reset counters
    if (primes_used != rec.get_prime()) {
      if (!rec.need_shift()) {
        rec.disable_shift();
        shift = rec.get_zi_shift_vec();
        std::fill(shift.begin(), shift.end(), 0);
      }

      std::cout << "Set new prime. Iterations for last prime: " << kk << ".\n";
      primes_used = std::max(primes_used, rec.get_prime());
      if(primes_used > 4) std::exit(-1);
      FFInt::set_new_prime(primes()[rec.get_prime()]);
      rec.generate_anchor_points();
      kk = 0;
    }

    // Always set the scaling variable to a random value
    t = rec.get_rand();

    // Add the shift to the scaling variable
    FFInt z1 = t + shift[0];

    for (uint j = 2; j <= n; ++j) {
      t_yis[j - 2] = t * rec.get_rand_zi(j, rec.get_zi_order()[j - 2]) + shift[j - 1];
    }

    std::vector<FFInt> yis(n);
    yis[0] = z1;

    for (uint j = 1; j < n; ++j) {
      yis[j] = t_yis[j - 1];
    }

    //FFInt num = singular_solver(yis); // example for n = 4 which uses the singular_solver
    //FFInt num = n_eq_1(z1); // example for n = 1
    //FFInt num = n_eq_4(yis); // example for n = 4 and the usage of the Chinese Remainder Theorem
    FFInt num = gghh(yis); // example for a large interpolation problem augmented with large coefficients
    //FFInt num = ggh(yis); // example for a three loop gg -> h integral coefficient

    ++kk;
    ++count;

    //FFInt num = yis[0].pow(15)*yis[1].pow(15)*yis[2].pow(15)*yis[3].pow(15)*yis[4].pow(15);// + yis[0].pow(20)*yis[1].pow(20)*yis[2].pow(20)*yis[3].pow(20);
    // Feed the algorithm with the current zi_order
    rec.feed(t, num, rec.get_zi_order(), primes_used);
    rec.interpolate();
  }

  std::cout << "Total numerical runs: " << count << ", primes used: " << primes_used + 1 << ".\n";
  //std::cout << rec.get_result();
  std::cout << "--------------------------------------------------------------\n";

  // Example for the reconstruction of a polynomial
  FFInt::set_new_prime(primes()[0]);
  n = 3;
  PolyReconst rec_poly(n, 5);

  std::cout << "Interpolating polynomial\n";
  std::cout << "--------------------------------------------------------------\n";
  std::vector<FFInt> yis(n);
  rec_poly.generate_anchor_points();

  // Initialize some counters
  count = 0;
  kk = 0;
  primes_used = 0;

  // Reconstruction loop
  while (!rec_poly.is_done()) {
    // If a new prime is needed, set it, generate new random variables
    // and reset counters
    if (primes_used != rec_poly.get_prime()) {
      std::cout << "Set new prime. Iterations for last prime: " << kk << ".\n";
      primes_used = std::max(primes_used, rec_poly.get_prime());

      FFInt::set_new_prime(primes()[rec_poly.get_prime()]);
      rec_poly.generate_anchor_points();
      kk = 0;
    }

    yis = rec_poly.get_rand_zi_vec(rec_poly.get_zi_order());

    FFInt num = pol_n_eq_3(yis);

    ++kk;
    ++count;

    rec_poly.feed(num, rec_poly.get_zi_order(), primes_used);
    rec_poly.interpolate();
  }

  std::cout << "Total numerical runs: " << count << ", primes used: " << primes_used + 1 << ".\n";
  //std::cout << rec_poly.get_result();
  std::cout << "--------------------------------------------------------------\n";

/*  FFInt::set_new_prime(primes()[0]);

  mpz_class p = primes()[0];
  FFInt a = -9;
  FFInt b = 1381219840000;

  mpz_class a_mpz = a.n;
  mpz_class b_mpz = b.n;
  mpz_class gcd_(gcd(a_mpz, b_mpz));
  a_mpz = a_mpz / gcd_;
  b_mpz = b_mpz / gcd_;

  std::cout << a << " / " << b << "\n";
  std::cout << a_mpz.get_str() << " / " << b_mpz.get_str() << "\n";

  std::pair<bool, RationalNumber> tmp = get_rational_coef((a/b).n, p);
  if (tmp.first) {
    std::cout << "Wang: " << tmp.second << "\n";
  } else {
    std::cout << "Wang failed\n";
  }

  tmp = get_rational_coef_mqrr((a/b).n, p);
  if (tmp.first) {
    std::cout << "Monagan: " << tmp.second << "\n";
  } else {
    std::cout << "Monagan failed\n";
  }*/

  return 0;
}
