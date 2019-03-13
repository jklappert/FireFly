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

#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "ThreadPool.hpp"

#include <list>
#include <tuple>

namespace firefly {
  typedef std::list<std::tuple<FFInt, std::vector<uint32_t>, std::future<std::vector<FFInt>>>> future_list;

  class Reconstructor {
  public:
    Reconstructor(uint32_t n_, uint32_t thr_n_, uint32_t verbosity_ = IMPORTANT);
    void enable_scan();
    void reconstruct();
    std::vector<RationalFunction> get_result();
    void black_box(std::vector<FFInt>& result, const std::vector<FFInt>& values);
    void set_tags(const std::vector<std::string>& tags_);
    void set_tags();
    void resume_from_saved_state(const std::vector<std::string>& file_paths_);
    enum verbosity_levels {SILENT, IMPORTANT, CHATTY};
  private:
    uint32_t n;
    uint32_t thr_n;
    uint32_t verbosity;
    std::vector<RatReconst> reconst {};
    bool scan = false;
    bool save_states = false;
    bool resume_from_state = false;
    std::vector<std::string> tags {};
    std::vector<std::string> file_paths {};
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
    uint32_t items_done = 0;
    uint32_t feeding_jobs = 0;
    uint32_t items = 0;
    uint32_t total_iterations = 0;
    uint32_t iteration = 0;
    RatReconst tmp_rec;
    std::vector<FFInt> shift {};
    /**
    *  Parses a prime number counter from a file
    *  @param file_name the file name
    */
    uint32_t parse_prime_number(std::string& file_name);
    void scan_for_shift();
    void start_first_runs();
    void run_until_done();
    void start_probe_jobs(const std::vector<uint32_t>& zi_order, const uint32_t start);
    void interpolate_job(RatReconst& rec);
  };
}
