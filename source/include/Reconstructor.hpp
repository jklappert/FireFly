#pragma once

#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "ThreadPool.hpp"

#include <list>
#include <tuple>

namespace firefly {
  typedef std::list<std::tuple<FFInt, std::vector<uint32_t>, std::future<std::vector<FFInt>>>> future_list;

  class Reconstructor{
  public:
    Reconstructor(uint32_t n_, uint32_t thr_n_, uint32_t mode_);
    void scan_for_sparsest_shift();
    void reconstruct();
    std::vector<RationalFunction> get_result_rf();
    std::vector<Polynomial> get_result_p();
    enum reconst_type {POLY, RAT};
  private:
    uint32_t n;
    uint32_t thr_n;
    uint32_t mode;
    std::vector<RatReconst> reconst {};
    bool scan = false;
    uint32_t prime_it = 0;
    ThreadPool tp;
    std::mutex mut;
    std::condition_variable cond;
    // list containing the parameters and the future of the parallel tasks; t, zi_order, future
    future_list probes {};
    uint32_t jobs_finished = 0;
    std::unordered_map<std::vector<uint32_t>, uint32_t, UintHasher> started_probes {};
    uint32_t fed_ones = 0;
    uint32_t probes_for_next_prime = 0;
    uint32_t feeding_jobs = 0;

    void start_probe_jobs(const std::vector<uint32_t> zi_order, const uint32_t start);
    void interpolate_job(RatReconst& reconst, uint32_t i);
  };
}
