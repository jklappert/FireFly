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
#include "ThreadPool.hpp"

#include <chrono>
#include <tuple>

namespace firefly {
  /**
   * @class BlackBoxBase
   * @brief The base class of the black box
   */
  class BlackBoxBase {
  public:
    /**
     *  A constructor for the BlackBoxBase class
     */
    BlackBoxBase() {}
    /**
     *  A destructor for the BlackBoxBase class
     */
    virtual ~BlackBoxBase() {}
    /**
     *  The evaluation of the black box. This function is called from Reconstructor.
     *  @param values The values to be inserted for the variables
     *  @return The result vector
     */
    virtual std::vector<FFInt> operator()(const std::vector<FFInt>& values) = 0;
    /**
     *  The evaluation of the black box in bunches. This function is called from Reconstructor.
     *  @param values_vec A vector (bunch) which holds several variable tuples at which the black box should be evaluated.
     *  @return The bunched result vector.
     */
    virtual std::vector<std::vector<FFInt>> operator()(const std::vector<std::vector<FFInt>>& values_vec) {
      std::vector<std::vector<FFInt>> results_vec;
      results_vec.reserve(values_vec.size());

      for (const auto & values : values_vec) {
        results_vec.emplace_back(operator()(values));
      }

      return results_vec;
    }
    /**
     *  Update internal variables of the black box when the prime field changes.
     *  This function is called from Reconstructor.
     */
    virtual void prime_changed() {}
  };

  typedef std::tuple<uint64_t, std::mutex*, int, RatReconst*> RatReconst_tuple;
  typedef std::list<RatReconst_tuple> RatReconst_list;
  typedef std::future<std::pair<std::vector<FFInt>, double>> probe_future;
  typedef std::list<std::tuple<FFInt, std::vector<uint32_t>, probe_future>> future_list;
  typedef std::future<std::pair<std::vector<std::vector<FFInt>>, double>> probe_future_bunch;
  typedef std::list<std::tuple<std::vector<FFInt>, std::vector<uint32_t>, probe_future_bunch>> future_list_bunch;

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
    Reconstructor(uint32_t n_, uint32_t thr_n_, BlackBoxBase& bb_, int verbosity_ = IMPORTANT);
    /**
     *  A constructor for the Reconstructor class
     *  @param n_ the number of parameters
     *  @param thr_n_ the number of threads being used during the reconstruction
     *  @param bunch_size_ the bunch size
     *  @param verbosity_ the verbosity level which can be chosen as SILENT (no output), IMPORTANT (only important output), and CHATTY (everything)
     */
    Reconstructor(uint32_t n_, uint32_t thr_n_, uint32_t bunch_size_, BlackBoxBase& bb_, int verbosity_ = IMPORTANT);
    /**
     *  A destructor for the Reconstructor class
     */
    ~Reconstructor();
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
     *  Enables the scan for a sparse shift at the beginning of the reconstruction
     */
    void enable_scan();
    /**
     *  Activate the safe interpolation mode where the function is completely interpolated in each prime field,
     *  no optimizations are used after the first prime field. Note that this mode cannot handle function changes
     *  which lead to coefficients which will become zero in all but one prime field.
     * < b>This option has to be set before calling resume_from_saved_state().< /b>
     */
    void set_safe_interpolation();
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
     *  @deprecated
     *  Resumes the reconstruction of a set of given functions
     *  @param file_paths_ a vector to the absolute paths to the intermediate results of reconstruction objects
     */
    void resume_from_saved_state(const std::vector<std::string>& file_paths_);
    /**
     *  Resumes the reconstruction of a set of functions which are located in a directory.
     *  The corresponding interpolation objects are created in the same order as they were defined in the prior
     *  run, thus requiring the black box to be probed in the same order.
     */
    void resume_from_saved_state();

    enum verbosity_levels {SILENT, IMPORTANT, CHATTY};
    enum RatReconst_status {RECONSTRUCTING, DONE, DELETED};

  private:
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point prime_start = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point last_print_time = std::chrono::high_resolution_clock::now();
    uint32_t n;
    uint32_t thr_n;
    uint32_t bunch_size = 1;
    uint32_t prime_it = 0;
    uint32_t total_iterations = 0;
    uint32_t iteration = 0;
    uint32_t probes_queued = 0;
    uint32_t probes_finished = 0;
    uint32_t fed_ones = 0;
    uint32_t probes_for_next_prime = 0;
    uint32_t items = 0;
    uint32_t items_done = 0;
    uint32_t items_new_prime = 0;
    uint32_t feed_jobs = 0;
    uint32_t interpolate_jobs = 0;
    uint32_t min_prime_keep_shift = 0; // TODO remove?
    BlackBoxBase& bb;
    int verbosity;
    double average_black_box_time = 0;
    double bunch_time;
    std::atomic<bool> scan = {false};
    std::atomic<bool> new_prime = {false};
    std::atomic<bool> done = {false};
    bool save_states = false;
    bool resume_from_state = false;
    bool safe_mode = false;
    bool one_done = false;
    bool one_new_prime = false;
    bool set_anchor_points = false;
    RatReconst_list reconst;
    std::vector<std::string> tags;
    std::vector<std::string> file_paths;
    ThreadPool tp;
    std::mutex future_control;
    std::mutex job_control;
    std::mutex feed_control;
    std::mutex print_control;
    std::mutex status_control;
    std::mutex mutex_probe_queue;
    std::mutex clean;
    std::mutex chosen_mutex;
    std::condition_variable condition_future;
    std::condition_variable condition_feed;
    // list containing the parameters and the future of the parallel tasks; t, zi_order, future
    future_list probes;
    std::queue<future_list::iterator> finished_probes_it;
    future_list_bunch probes_bunch;
    std::queue<future_list_bunch::iterator> finished_probes_bunch_it;
    std::vector<uint32_t> bunch_zi_order;
    std::vector<FFInt> bunch_t;
    std::vector<std::vector<FFInt>> bunch;
    std::unordered_map<std::vector<uint32_t>, uint32_t, UintHasher> started_probes;
    std::unordered_map<std::vector<uint32_t>, std::unordered_set<uint64_t>, UintHasher> chosen_t;
    RatReconst tmp_rec;
    std::vector<FFInt> shift;
    /**
    *  Scan the black-box functions for a sparse shift
    */
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
     * Gets a probe from probes, probes_bunch, or bunch
     * @param t is set to the t of the returned probe
     * @param zi_order is set to the zi_order of the returned probe
     * @param probe is set to point to a probe from probes, probes_bunch, or bunch
     * @param time is set to the time of the returned probe
     */
    void get_probe(FFInt& t, std::vector<uint32_t>& zi_order, std::vector<FFInt>* probe, double& time);
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
#ifdef WITH_MPI
    int world_size;
    double tmp_average_black_box_time = 0.;
    uint32_t total_thread_count = 0;
    uint32_t tmp_total_iterations = 0;
    uint64_t ind = 0;
    bool proceed = false;
    std::unordered_map<uint64_t, std::pair<FFInt, std::vector<uint32_t>>> index_map;
    std::unordered_map<int, uint64_t> nodes;
    std::queue<std::pair<int, uint64_t>> empty_nodes;
    std::condition_variable cond_val;
    std::atomic<bool> new_jobs = {false};
    std::queue<std::vector<uint64_t>> value_queue;
    std::queue<std::pair<uint64_t, std::vector<FFInt>>> results_queue;
    /**
     *  TODO
     */
    void get_a_job();
    void mpi_setup();
    void send_first_jobs();
    void mpi_communicate();
#endif
  };
}
