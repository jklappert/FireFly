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

#include "BlackBoxBase.hpp"
#include "FFIntVec.hpp"
#include "ReconstHelper.hpp"
#include "ThreadPool.hpp"
#include "utils.hpp"

#include <mpi.h>

namespace firefly {
  const int master = 0;
  const uint64_t buffer = 2;
  enum MPI_tags {VALUES, RESULT, NEW_PRIME, TIMING, END};

  template<typename BlackBoxTemp>
  class MPIWorker {
  public:
    // TODO
    MPIWorker(uint32_t n_, uint32_t thr_n_, BlackBoxBase<BlackBoxTemp>& bb_);
    MPIWorker(uint32_t n_, uint32_t thr_n_, uint32_t bunch_size_, BlackBoxBase<BlackBoxTemp>& bb_);

  private:
    const uint32_t n;
    const uint32_t thr_n;
    const uint32_t bunch_size = 1;
    uint32_t total_iterations = 0;
    uint64_t tasks = 0;
    double average_black_box_time = 0.;
    ThreadPool tp;
    BlackBoxBase<BlackBoxTemp>& bb;
    std::vector<uint64_t> results;
    std::mutex mut;
    std::condition_variable cond;

    // TODO
    void run();
    void communicate();
    template<uint32_t N>
    void queue_new_job(const std::vector<uint64_t>& values_list, const uint32_t start);
    template<typename FFIntTemp>
    void compute(const std::vector<uint64_t>& index_vec, const std::vector<FFIntTemp>& values_vec);
  };

  template<typename BlackBoxTemp>
  MPIWorker<BlackBoxTemp>::MPIWorker(uint32_t n_, uint32_t thr_n_, BlackBoxBase<BlackBoxTemp>& bb_) : n(n_), thr_n(thr_n_), tp(thr_n_), bb(bb_) {
    run();
  }

  template<typename BlackBoxTemp>
  MPIWorker<BlackBoxTemp>::MPIWorker(uint32_t n_, uint32_t thr_n_, uint32_t bunch_size_, BlackBoxBase<BlackBoxTemp>& bb_) : n(n_), thr_n(thr_n_), bunch_size(bunch_size_), tp(thr_n_), bb(bb_) {
    run();
  }

  template<typename BlackBoxTemp>
  void MPIWorker<BlackBoxTemp>::run() {
    uint32_t prime;
    MPI_Bcast(&prime, 1, MPI_UINT32_T, master, MPI_COMM_WORLD);

    FFInt::set_new_prime(primes()[prime]);

    if (prime != 0) {
      bb.prime_changed_internal();
    }

    //std::cout << "worker " << FFInt::p << "\n";

    communicate();
  }

  template<typename BlackBoxTemp>
  void MPIWorker<BlackBoxTemp>::communicate() {
    while (true) {
      //std::cout << "w loop\n";
      std::unique_lock<std::mutex> lock(mut);
      std::vector<uint64_t> tmp_results;

      MPI_Request request;

      if (tasks == 0 && results.empty()) {
        lock.unlock();

        uint64_t free = buffer * static_cast<uint64_t>(thr_n);
        MPI_Isend(&free, 1, MPI_UINT64_T, master, RESULT, MPI_COMM_WORLD, &request);
      } else {
        while (results.empty()) {
          cond.wait(lock);
        }

        results.emplace_back(buffer * static_cast<uint64_t>(thr_n) - tasks);
        tmp_results = std::move(results);
        results.clear();

        lock.unlock();

        //std::cout << "worker sending " << tmp_results.size() - 1 << " items\n";

        MPI_Isend(&tmp_results[0], static_cast<int>(tmp_results.size()), MPI_UINT64_T, master, RESULT, MPI_COMM_WORLD, &request);
      }

      MPI_Status status;
      MPI_Probe(master, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      if (status.MPI_TAG == VALUES) {
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        int amount;
        MPI_Get_count(&status, MPI_UINT64_T, &amount);

        if (static_cast<uint64_t>(amount) % static_cast<uint64_t>(n + 1) != 0) {
          ERROR_MSG("Corrupted values recieved!");
          std::exit(EXIT_FAILURE);
        }

        std::vector<uint64_t> values_list;
        values_list.reserve(amount);

        MPI_Recv(&values_list[0], amount, MPI_UINT64_T, master, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        //std::cout << "recieved:\n";
        //for (int i = 0; i != amount; ++i) {
        //  std::cout << values_list[i] << "\n";
        //}

        uint32_t new_tasks = static_cast<uint32_t>(amount) / (n + 1);
        uint32_t started = 0;

        //std::cout << "new tasks " << new_tasks << "\n";

        while (started != new_tasks) {
          uint32_t next_bunch = compute_bunch_size(new_tasks - started, thr_n, bunch_size);

          //std::cout << next_bunch << "\n";

          switch(next_bunch) {
            case 1:
              queue_new_job<1>(values_list, started);
              break;
            case 2:
              queue_new_job<2>(values_list, started);
              break;
            case 4:
              queue_new_job<4>(values_list, started);
              break;
            case 8:
              queue_new_job<8>(values_list, started);
              break;
            case 16:
              queue_new_job<16>(values_list, started);
              break;
            case 32:
              queue_new_job<32>(values_list, started);
              break;
            case 64:
              queue_new_job<64>(values_list, started);
              break;
            case 128:
              queue_new_job<128>(values_list, started);
              break;
            case 256:
              queue_new_job<256>(values_list, started);
              break;
          }

          started += next_bunch;

          lock.lock();

          ++tasks;

          lock.unlock();
        }
      } else if (status.MPI_TAG == NEW_PRIME) {
        //std::cout << "worker new prime\n";

        uint64_t prime;
        MPI_Recv(&prime, 1, MPI_UINT64_T, master, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        tp.kill_all();
        results.clear();
        tasks = 0;

        if (prime != 0) {
          FFInt::set_new_prime(primes()[static_cast<uint32_t>(prime)]);
          bb.prime_changed_internal();
        }

        //std::cout << "worker " << FFInt::p << "\n";

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Wait(&request, MPI_STATUS_IGNORE);

        double timing[2] = {average_black_box_time, static_cast<double>(total_iterations)};

        MPI_Send(timing, 2, MPI_DOUBLE, master, TIMING, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
      } else if (status.MPI_TAG == END) {
        uint64_t tmp;
        MPI_Recv(&tmp, 1, MPI_UINT64_T, master, END, MPI_COMM_WORLD, &status);

        tp.kill_all();
        results.clear();
        tasks = 0;

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        MPI_Barrier(MPI_COMM_WORLD);
        break;
      } else {
        ERROR_MSG("Unknown MPI_TAG!");
        std::exit(EXIT_FAILURE);
      }
    }
  }

  template<typename BlackBoxTemp>
  template<uint32_t N>
  void MPIWorker<BlackBoxTemp>::queue_new_job(const std::vector<uint64_t>& values_list, const uint32_t start) {
    std::vector<uint64_t> indices;
    indices.reserve(N);

    std::vector<FFIntVec<N>> values_vec(n);

    for (uint32_t i = 0; i != N; ++i) {
      indices.emplace_back(values_list[(n + 1) * (start + i)]);

      for (uint32_t j = 0; j != n; ++j) {
        values_vec[j][i] = values_list[(n + 1) * (start + i) + j + 1];
      }
    }

    tp.run_task([this, indices = std::move(indices), values_vec = std::move(values_vec)]() {
      compute(indices, values_vec);
    });
  }

  template<typename BlackBoxTemp>
  template<typename FFIntTemp>
  void MPIWorker<BlackBoxTemp>::compute(const std::vector<uint64_t>& index_vec, const std::vector<FFIntTemp>& values_vec) {
    //std::cout << "calc with:\n";
    //for (size_t j = 0; j != values_vec.front().size(); ++j) {
    //  for (size_t i = 0; i != values_vec.size(); ++i) {
    //     std::cout << values_vec[i][j].n << " ";
    //  }
    //  std::cout << "\n";
    //}

    auto time0 = std::chrono::high_resolution_clock::now();

    std::vector<FFIntTemp> result = bb.eval(values_vec);

    auto time1 = std::chrono::high_resolution_clock::now();

    std::vector<uint64_t> result_uint;
    result_uint.reserve(result.front().size() * (1 + result.size()));

    //std::cout << "emplacing " << static_cast<size_t>(bunch_size) * (1 + result.front().size()) << " " << bunch_size << " " << result.front().size() << "\n";

    for (size_t j = 0; j != result.front().size(); ++j) {
      result_uint.emplace_back(index_vec[j]);
      //std::cout << "index " << result_uint.back() << "\n";

      for (size_t i = 0; i != result.size(); ++i) {
        result_uint.emplace_back(result[i][j].n);
        //std::cout << result_uint.back() << "\n";
      }
    }

    std::unique_lock<std::mutex> lock(mut);
    ++total_iterations;

    auto time = std::chrono::duration<double>(time1 - time0).count();
    average_black_box_time = (average_black_box_time * (total_iterations - 1) + time) / total_iterations;
    results.insert(results.end(), result_uint.begin(), result_uint.end());
    //std::cout << "res back " << results.back() << "\n";
    //std::cout << "results size " << results.size() << "\n";
    --tasks;

    cond.notify_one();
  }
}
