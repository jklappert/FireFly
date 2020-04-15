//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2020  Jonas Klappert and Fabian Lange
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
  constexpr int master = 0;
  constexpr uint64_t buffer = 2;
  enum MPI_tags {VALUES, RESULT, NEW_PRIME, TIMING, FACTORS, END};

  /**
   * @class MPIWorker
   * @brief A class to compute probes of the black box BlackBoxTemp on worker processes with MPI
   */
  template<typename BlackBoxTemp>
  class MPIWorker {
  public:
    /**
     *  A constructor for the MPIWorker class
     *  Automatically runs the worker
     *  @param n_ the number of parameters
     *  @param thr_n_ the number of threads
     *  @param bb_ An instance of a BlackBoxBase class
     */
    MPIWorker(uint32_t n_, uint32_t thr_n_, BlackBoxBase<BlackBoxTemp>& bb_);
    /**
     *  A constructor for the MPIWorker class
     *  Automatically runs the worker
     *  @param n_ the number of parameters
     *  @param thr_n_ the number of threads
     *  @param bunch_size_ the maximum bunch size to be used
     *  @param bb_ An instance of a BlackBoxBase class
     */
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
    std::unordered_map<uint32_t, ShuntingYardParser> parsed_factors {};
    std::mutex mut;
    std::condition_variable cond;

    /**
     *  Runs the worker
     */
    void run();
    /**
     *  Communicate with the Reconstructor class
     */
    void communicate();
    /**
     *  Templated function to queue N new probes at the values in values_list
     *  starting at position (n+1) * start
     *  @param values_list the list of values
     *  @param start the position where to start in values_list
     */
    template<uint32_t N>
    void queue_new_job(const std::vector<uint64_t>& values_list, const uint32_t start);
    /**
     *  Computes the probe index at the values values_vec
     *  @param index index of the probe
     *  @param values_vec values at which the probe is evaluated
     */
    void compute(const uint64_t index, const std::vector<FFInt>& values_vec);
    /**
     *  Computes the N probes in index_vec at the values values_vec as an FFIntVec of length N
     *  @param index_vec indices of the probes
     *  @param values_vec N sets of values at which the probes are evaluated
     */
    template<uint32_t N>
    void compute(const std::vector<uint64_t>& index_vec, const std::vector<FFIntVec<N>>& values_vec);
  };

  template<typename BlackBoxTemp>
  MPIWorker<BlackBoxTemp>::MPIWorker(uint32_t n_, uint32_t thr_n_, BlackBoxBase<BlackBoxTemp>& bb_) : n(n_), thr_n(thr_n_), tp(thr_n_), bb(bb_) {
    if (thr_n == 0) {
      ERROR_MSG("Worker started without threads. Abort.");
      std::exit(EXIT_FAILURE);
    }

    run();
  }

  template<typename BlackBoxTemp>
  MPIWorker<BlackBoxTemp>::MPIWorker(uint32_t n_, uint32_t thr_n_, uint32_t bunch_size_, BlackBoxBase<BlackBoxTemp>& bb_) : n(n_), thr_n(thr_n_), bunch_size(bunch_size_), tp(thr_n_), bb(bb_) {
    if (thr_n == 0) {
      ERROR_MSG("Worker started without threads. Abort.");
      std::exit(EXIT_FAILURE);
    }

    run();
  }

  template<typename BlackBoxTemp>
  void MPIWorker<BlackBoxTemp>::run() {
    uint32_t prime;
    MPI_Bcast(&prime, 1, MPI_UINT32_T, master, MPI_COMM_WORLD);

    FFInt::set_new_prime(primes()[prime]);
    bb.prime_changed_internal();

    communicate();
  }

  template<typename BlackBoxTemp>
  void MPIWorker<BlackBoxTemp>::communicate() {
    while (true) {
      std::unique_lock<std::mutex> lock(mut);
      std::vector<uint64_t> tmp_results;

      MPI_Request request;

      if (tasks == 0 && results.empty()) {
        lock.unlock();

        uint64_t free = buffer * static_cast<uint64_t>(thr_n);
        MPI_Isend(&free, 1, MPI_UINT64_T, master, RESULT, MPI_COMM_WORLD, &request);
      } else {
        cond.wait(lock, [this](){return !results.empty();});

        results.emplace_back(buffer * static_cast<uint64_t>(thr_n) - tasks);
        tmp_results = std::move(results);
        results.clear();

        lock.unlock();

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

        uint32_t new_tasks = static_cast<uint32_t>(amount) / (n + 1);
        uint32_t started = 0;

        while (started != new_tasks) {
          uint32_t next_bunch = compute_bunch_size(new_tasks - started, thr_n, bunch_size);

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
            /*case 256:
              queue_new_job<256>(values_list, started);
              break;*/
          }

          started += next_bunch;

          lock.lock();

          ++tasks;

          lock.unlock();
        }
      } else if (status.MPI_TAG == NEW_PRIME) {
        tp.kill_all();
        results.clear();
        tasks = 0;

        // receive the prime signal but not the actual prime
        uint64_t prime;
        MPI_Recv(&prime, 1, MPI_UINT64_T, master, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        MPI_Wait(&request, MPI_STATUS_IGNORE);

        uint64_t last_package_of_prime = 0;
        MPI_Send(&last_package_of_prime, 1, MPI_UINT64_T, master, RESULT, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        double timing[2] = {average_black_box_time, static_cast<double>(total_iterations)};

        MPI_Send(timing, 2, MPI_DOUBLE, master, TIMING, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        total_iterations = 0;
        average_black_box_time = 0.;

        MPI_Status np_or_fac;
        MPI_Probe(master, MPI_ANY_TAG, MPI_COMM_WORLD, &np_or_fac);

        if (np_or_fac.MPI_TAG == FACTORS) {
          uint64_t tmp;
          MPI_Recv(&tmp, 1, MPI_UINT64_T, master, FACTORS, MPI_COMM_WORLD, &np_or_fac);

          parsed_factors.clear();

          std::vector<std::string> vars (n);

          for (size_t i = 0; i != n; ++i) {
            vars[i] = "x" + std::to_string(i + 1);
          }

          // Receive factors
          while (true) {
            int amount;
            int multiple = 1;
            MPI_Bcast(&amount, 1, MPI_INT, master, MPI_COMM_WORLD);

            if (amount == -1) {
              break;
            } else if (amount < -1) {
              multiple = -amount;
            }

            std::string tmp_fac_s = "";

            for (int j = 0; j != multiple; ++j) {
              if (multiple != 1) {
                MPI_Bcast(&amount, 1, MPI_INT, master, MPI_COMM_WORLD);
              }

              char* fac_c = new char[amount];
              MPI_Bcast(fac_c, amount, MPI_CHAR, master, MPI_COMM_WORLD);

              for (int i = 0; i != amount; ++i) {
                tmp_fac_s += fac_c[i];
              }

              delete[] fac_c;
            }

            uint32_t function_number;
            MPI_Bcast(&function_number, 1, MPI_UINT32_T, master, MPI_COMM_WORLD);

            ShuntingYardParser parser = ShuntingYardParser();
            parser.parse_function(tmp_fac_s, vars);
            parsed_factors.emplace(function_number, parser);
          }

          MPI_Barrier(MPI_COMM_WORLD);
        }

        // receive next prime
        MPI_Recv(&prime, 1, MPI_UINT64_T, master, NEW_PRIME, MPI_COMM_WORLD, &np_or_fac);

        FFInt::set_new_prime(primes()[static_cast<uint32_t>(prime)]);
        bb.prime_changed_internal();

        if (!parsed_factors.empty()) {
          for (auto & el : parsed_factors) {
            el.second.precompute_tokens();
          }
        }

        MPI_Barrier(MPI_COMM_WORLD);
      } else if (status.MPI_TAG == END) {
        tp.kill_all();
        results.clear();
        tasks = 0;

        uint64_t tmp;
        MPI_Recv(&tmp, 1, MPI_UINT64_T, master, END, MPI_COMM_WORLD, &status);

        MPI_Wait(&request, MPI_STATUS_IGNORE);

        uint64_t last_package_of_prime = 0;
        MPI_Send(&last_package_of_prime, 1, MPI_UINT64_T, master, RESULT, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        double timing[2] = {average_black_box_time, static_cast<double>(total_iterations)};

        MPI_Send(timing, 2, MPI_DOUBLE, master, TIMING, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);

        total_iterations = 0;
        average_black_box_time = 0.;

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
    if (N != 1) {
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
        compute<N>(indices, values_vec);
      });
    } else {
      std::vector<FFInt> values_vec;
      values_vec.reserve(n);

      uint64_t index = values_list[(n + 1) * start];

      for (uint32_t j = 0; j != n; ++j) {
        values_vec.emplace_back(values_list[(n + 1) * start + j + 1]);
      }

      tp.run_task([this, index, values_vec = std::move(values_vec)]() {
        compute(index, values_vec);
      });
    }
  }

  template<typename BlackBoxTemp>
  void MPIWorker<BlackBoxTemp>::compute(const uint64_t index, const std::vector<FFInt>& values_vec) {
    auto time0 = std::chrono::high_resolution_clock::now();

    std::vector<FFInt> result = bb.eval(values_vec);

    auto time1 = std::chrono::high_resolution_clock::now();

    std::vector<uint64_t> result_uint;
    result_uint.reserve(result.size());

    result_uint.emplace_back(index);

    for (size_t i = 0; i != result.size(); ++i) {
      if (parsed_factors.find(i) != parsed_factors.end()) {
        auto fac = parsed_factors[i].evaluate_pre(values_vec);
        result[i] /= fac[0];
      }

      result_uint.emplace_back(result[i].n);
    }

    auto time = std::chrono::duration<double>(time1 - time0).count();

    std::lock_guard<std::mutex> lock(mut);

    ++total_iterations;
    average_black_box_time = (average_black_box_time * (total_iterations - 1) + time) / total_iterations;
    results.insert(results.end(), result_uint.begin(), result_uint.end());
    --tasks;

    cond.notify_one();
  }

  template<typename BlackBoxTemp>
  template<uint32_t N>
  void MPIWorker<BlackBoxTemp>::compute(const std::vector<uint64_t>& index_vec, const std::vector<FFIntVec<N>>& values_vec) {
    auto time0 = std::chrono::high_resolution_clock::now();

    std::vector<FFIntVec<N>> result = bb.eval(values_vec);

    auto time1 = std::chrono::high_resolution_clock::now();

    std::vector<uint64_t> result_uint;
    result_uint.reserve(result.front().size() * (1 + result.size()));

    if (!parsed_factors.empty()) {
      // Remove factors from bb result
      for (size_t i = 0; i != result.size(); ++i) {
        if (parsed_factors.find(i) != parsed_factors.end()) {
          auto res = parsed_factors[i].evaluate_pre(values_vec);

          for (size_t j = 0; j != N; ++j) {
            result[i][j] /= res[0][j];
          }
        }
      }
    }

    for (size_t j = 0; j != N; ++j) {
      result_uint.emplace_back(index_vec[j]);

      for (size_t i = 0; i != result.size(); ++i) {
        result_uint.emplace_back(result[i][j].n);
      }
    }

    auto time = std::chrono::duration<double>(time1 - time0).count();

    std::lock_guard<std::mutex> lock(mut);

    total_iterations += result.front().size();
    average_black_box_time = (average_black_box_time * (total_iterations - 1) + time) / total_iterations;
    results.insert(results.end(), result_uint.begin(), result_uint.end());
    --tasks;

    cond.notify_one();
  }
}
