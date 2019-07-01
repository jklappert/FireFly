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
  class BlackBoxBase {
  public:
      BlackBoxBase() {};
      virtual std::vector<FFInt> operator()(const std::vector<FFInt> & values) = 0;
      virtual void prime_changed() = 0;
  };

  enum RatReconst_status {DEFAULT, DONE, DELETED};

  typedef std::tuple<uint64_t, std::mutex *, int, RatReconst *> RatReconst_tuple;
  typedef std::list<RatReconst_tuple> RatReconst_list;
  typedef std::list<std::tuple<FFInt, std::vector<uint32_t>, std::future<std::pair<std::vector<FFInt>, double>>>> future_list;
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
    Reconstructor(uint32_t n_, uint32_t thr_n_, BlackBoxBase & bb_, uint32_t verbosity_ = IMPORTANT);
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
     *  @return a vector of the already reconstructed rational functions and its tag;
     *  they are removed from the internal memory afterwards
     */
    std::vector<std::pair<std::string, RationalFunction>> get_early_results();
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
    /**
     *  Activate the safe interpolation mode where the function is completely interpolated in each prime field,
     *  no optimizations are used after the first prime field
     */
    void set_safe_interpolation();
    enum verbosity_levels {SILENT, IMPORTANT, CHATTY};
  private:
    uint32_t n;
    uint32_t thr_n;
    int verbosity;
    RatReconst_list reconst {};
    bool scan = false;
    bool save_states = false;
    bool resume_from_state = false;
    std::vector<std::string> tags {};
    std::vector<std::string> file_paths {};
    bool safe_mode = false;
    uint32_t prime_it = 0;
    BlackBoxBase & bb;
    ThreadPool tp;
    std::mutex future_control;
    std::mutex job_control;
    std::mutex feed_control;
    std::mutex print_control;
    std::mutex status_control;
    std::mutex clean;
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
    bool one_done = false;
    double average_black_box_time;
    RatReconst tmp_rec;
    std::vector<FFInt> shift {};
    /**
    *  Parses a prime number counter from a file
    *  @param file_name the file name
    */
    uint32_t parse_prime_number(std::string& file_name);
    void scan_for_shift();
    /**
     *  Initializes vector of reconstruction objects and starts first probes
     */
    void start_first_runs();
    /**
     *  Starts new jobs until the reconstruction is done
     */
    void run_until_done();
    /**
     *  Queues a number of jobs corresponding to a given zi_order
     *  @param zi_order the order of which a given number of jobs should be started
     *  @param to_start the number of jobs which should be queued
     */
    void start_probe_jobs(const std::vector<uint32_t>& zi_order, const uint32_t to_start);
    /**
     *  Feeds the reconstruction objects
     *  @param zi_order the order at which the black box was probed
     *  @param t the value of the homogenization variable t
     *  @param probe a vector of black box probes in an immutable order
     */
    void feed_job(const std::vector<uint32_t> zi_order, const FFInt t, std::vector<FFInt>* probe);
    /**
     *  Interpolates a RatReconst and queues new jobs if required
     *  @param rec a reference to the RatReconst
     */
    void interpolate_job(RatReconst_tuple& rec);
    /**
     *  Removes all RatReconst from reconst which are flagged as DELETE
     */
    void clean_reconst();
  };
}
