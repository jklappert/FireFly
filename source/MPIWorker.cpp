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

#include "MPIWorker.hpp"
#include "ReconstHelper.hpp"

namespace firefly {
  MPIWorker::MPIWorker(uint32_t n_, uint32_t thr_n_, BlackBoxBase& bb_) : n(n_), thr_n(thr_n_), tp(thr_n_), bb(bb_) {
    uint32_t prime;
    MPI_Bcast(&prime, 1, MPI_UINT32_T, master, MPI_COMM_WORLD);

    FFInt::set_new_prime(primes()[prime]);

    if (prime != 0) {
      bb.prime_changed();
    }

    //std::cout << "worker " << FFInt::p << "\n";

    //MPI_Bcast(&bb_size, 1, MPI_INT, master, MPI_COMM_WORLD);

    communicate();
  }

  void MPIWorker::communicate() {
    while (true) {
      //std::cout << "w loop\n";
      std::unique_lock<std::mutex> lock(mut);

      MPI_Request request;

      if (tasks == 0 && results.empty()) {
        lock.unlock();

        uint64_t free = buffer * static_cast<uint64_t>(thr_n);
        MPI_Isend(&free, 1, MPI_UINT64_T, master, RESULT, MPI_COMM_WORLD, &request);
      } else {
        while (results.empty()) {
          cond.wait(lock);
        }

        results.emplace_back(buffer * thr_n - tasks);
        std::vector<uint64_t> tmp_results = std::move(results);
        results.clear();

        lock.unlock();

        //std::cout << "worker sending " << (tmp_results.size() - 1) / (1 + bb_size) << " results\n";

        MPI_Isend(tmp_results.data(), static_cast<int>(tmp_results.size()), MPI_UINT64_T, master, RESULT, MPI_COMM_WORLD, &request);
      }

      MPI_Status status;
      MPI_Probe(master, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      MPI_Wait(&request, MPI_STATUS_IGNORE);

      if (status.MPI_TAG == VALUES) {
        int amount;
        MPI_Get_count(&status, MPI_UINT64_T, &amount);

        if (static_cast<uint64_t>(amount) % static_cast<uint64_t>(n + 1) != 0) {
          ERROR_MSG("Corrupted values recieved!");
          std::exit(EXIT_FAILURE);
        }

        uint64_t* values_list = new uint64_t[amount];

        MPI_Recv(values_list, amount, MPI_UINT64_T, master, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        uint64_t new_tasks = static_cast<uint64_t>(amount) / static_cast<uint64_t>(n + 1);

        //std::cout << "worker starting " << new_tasks << " jobs\n";

        lock.lock();

        tasks += new_tasks;

        lock.unlock();

        for (uint64_t i = 0; i != new_tasks; ++i) {
          uint64_t index = values_list[i * (n + 1)];

          std::vector<FFInt> values;

          for (uint32_t j = 1; j != n + 1; ++j) {
            values.emplace_back(values_list[i * (n + 1) + j]);
          }

          tp.run_task([this, index, values]() {
            compute(index, values);
          });
        }

        delete[] values_list;
      } else if (status.MPI_TAG == NEW_PRIME) {
        //std::cout << "worker new prime\n";

        uint64_t prime;
        MPI_Recv(&prime, 1, MPI_UINT64_T, master, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        tp.kill_all();
        results.clear();
        tasks = 0;

        if (prime != 0) {
          FFInt::set_new_prime(primes()[static_cast<uint32_t>(prime)]);
          bb.prime_changed();
        }

        //std::cout << "worker " << FFInt::p << "\n";

        MPI_Barrier(MPI_COMM_WORLD);

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
        break;
      } else {
        ERROR_MSG("Unknown MPI_TAG!");
        std::exit(EXIT_FAILURE);
      }
    }
  }

  void MPIWorker::compute(const uint64_t index, const std::vector<FFInt>& values) {
    auto time0 = std::chrono::high_resolution_clock::now();

    std::vector<FFInt> result = bb(values);

    auto time1 = std::chrono::high_resolution_clock::now();

    std::vector<uint64_t> result_uint;
    result_uint.reserve(1 + result.size());

    result_uint.emplace_back(index);

    for (size_t i = 0; i != result.size(); ++i) {
      result_uint.emplace_back(result[i].n);
    }

    std::unique_lock<std::mutex> lock(mut);
    ++total_iterations;

    auto time = std::chrono::duration<double>(time1 - time0).count();
    average_black_box_time = (average_black_box_time * (total_iterations - 1) + time) / total_iterations;
    results.insert(results.end(), result_uint.begin(), result_uint.end());
    --tasks;

    cond.notify_one();
  }
}
