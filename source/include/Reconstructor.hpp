#pragma once

#include "RatReconst.hpp"
#include "ReconstHelper.hpp"

namespace firefly {

  class Reconstructor{
  public:
    Reconstructor(uint32_t n_, uint32_t thr_n_, uint32_t mode_);
    void reconstruct();
    void scan_for_sparsest_shift();
    std::vector<RationalFunction> get_result_rf();
    std::vector<Polynomial> get_result_p();
    enum reconst_type {POLY, RAT};
  private:
    uint32_t n;
    uint32_t thr_n;
    uint32_t mode;
  };
}