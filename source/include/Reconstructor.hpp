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
  /**
   * @class Reconstructor
   * @brief A class to reconstruct functions from its values
   */
  class Reconstructor {
  public:
    /**
     *  A constructor for the Reconstructor class
     *  @param n_ the number of parameters
     *  @param thr_n_ the number of threads being used during the reconstruction
     *  @param verbosity_ the verbosity level which can be chosen as SILENT (no output), IMPORTANT (only important output), and CHATTY (everything)
     */
    Reconstructor(uint32_t n_, uint32_t thr_n_, uint32_t verbosity_ = IMPORTANT);
    /**
     *  Enables the scan for a sparse shift at the beginning of the reconstruction
     */
    void enable_scan();
    /**
     *  Starts the reconstruction
     */
    void reconstruct();
    /**
     *  @return the vector of reconstructed rational functions
     */
    std::vector<RationalFunction> get_result();
    /**
     *  The black box functions which gets called by the Reconstructor class to evaluate probes
     *  @param result a vector to be filled by the user with the probes in an immutable ordering
     *  @param values the parameter point at which the black box should be probed
     */
    void black_box(std::vector<FFInt>& result, const std::vector<FFInt>& values);
    /**
     *  Sets user defined tags for each reconstruction object and saves intermediate results after each prime field
     *  @param tags_ a vector of user defined tags in an immutable ordering
     */
    void set_tags(const std::vector<std::string>& tags_);
    /**
     *  Sets default tags for each reconstruction object and saves intermediate results after each prime field
     */
    void set_tags();
    /**
     *  Resumes the reconstruction of a set of given functions
     *  @param file_paths_ a vector to the absolute paths to the intermediate results of reconstruction objects
     */
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
    std::mutex future_control;
    std::mutex job_control;
    std::mutex feed_control;
    std::mutex print_control;
    std::mutex status_control;
    std::condition_variable condition_future;
    std::condition_variable condition_feed;
    // list containing the parameters and the future of the parallel tasks; t, zi_order, future
    future_list probes {};
    uint32_t jobs_finished = 0;
    std::unordered_map<std::vector<uint32_t>, uint32_t, UintHasher> started_probes {};
    uint32_t fed_ones = 0;
    uint32_t probes_for_next_prime = 0;
    uint32_t items = 0;
    uint32_t items_done = 0;
    uint32_t items_new_prime = 0;
    uint32_t feed_jobs = 0;
    uint32_t interpolate_jobs = 0;
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
    void start_probe_jobs(const std::vector<uint32_t>& zi_order, const uint32_t to_start);
    void feed_job(const std::vector<uint32_t> zi_order, const FFInt t, std::vector<FFInt>* probe);
    void interpolate_job(RatReconst& rec);
  };
}
