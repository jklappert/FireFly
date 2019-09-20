//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2019  Jonas Klappert and Fabian Lange
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//==================================================================================

#pragma once

#include "Reconstructor.hpp"

#include <mpi.h>

namespace firefly {
  const int master = 0;
  enum MPI_tags {VALUES, RESULT, NEW_PRIME, SLOTS, TIMING, END};

  class MPIWorker {
  public:
    MPIWorker(uint32_t n_, uint32_t thr_n_, BlackBoxBase& bb_);

  private:
    const uint32_t n;
    const uint32_t thr_n;
    uint32_t total_iterations = 0;
    double average_black_box_time = 0.;
    ThreadPool tp;
    BlackBoxBase& bb;
    std::vector<uint64_t> results;
    std::mutex mut;
    std::condition_variable cond;
    uint64_t tasks = 0;

    void communicate();
    void compute(const uint64_t index, const std::vector<FFInt>& values);
  };
}
