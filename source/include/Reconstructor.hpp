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
#include "gzstream.hpp"
#include "ParserUtils.hpp"
#include "RatReconst.hpp"
#include "ReconstHelper.hpp"
#include "ThreadPool.hpp"
#include "tinydir.h"
#include "utils.hpp"
#include "version.hpp"

#if WITH_MPI
#include "MPIWorker.hpp"
#endif

#ifdef FLINT
#include <flint/fmpz_poly.h>
#endif

#include <chrono>
#include <tuple>
#include <sys/stat.h>

namespace firefly {
  typedef std::tuple<uint64_t, std::mutex*, int, RatReconst*> RatReconst_tuple;
  typedef std::list<RatReconst_tuple> RatReconst_list;

  /**
   * @class Reconstructor
   * @brief A class to reconstruct functions from its values
   */
  template<typename BlackBoxTemp>
  class Reconstructor {
  public:
    /**
     *  A constructor for the Reconstructor class
     *  @param n_ the number of parameters
     *  @param thr_n_ the number of threads being used during the reconstruction
     *  @param bb_ An instance of a BlackBoxBase class
     *  @param verbosity_ the verbosity level which can be chosen as SILENT (no output), IMPORTANT (only important output), and CHATTY (everything)
     */
    Reconstructor(const uint32_t n_, const uint32_t thr_n_, BlackBoxBase<BlackBoxTemp>& bb_, const int verbosity_ = IMPORTANT);
    /**
     *  A constructor for the Reconstructor class
     *  @param n_ the number of parameters
     *  @param thr_n_ the number of threads being used during the reconstruction
     *  @param bunch_size_ the bunch size
     *  @param bb_ An instance of a BlackBoxBase class
     *  @param verbosity_ the verbosity level which can be chosen as SILENT (no output), IMPORTANT (only important output), and CHATTY (everything)
     */
    Reconstructor(const uint32_t n_, const uint32_t thr_n_, const uint32_t bunch_size_, BlackBoxBase<BlackBoxTemp>& bb_, const int verbosity_ = IMPORTANT);
    /**
     *  A destructor for the Reconstructor class
     */
    ~Reconstructor();
    /**
     *  Default constructor
     */
    Reconstructor() {};
    /**
     *  Starts the reconstruction
     *  @param prime_counter sets how many interpolations have to be performed at most
     */
    void reconstruct(uint32_t prime_counter = 100);
    /**
     *  @return the vector of reconstructed rational functions
     */
    std::vector<RationalFunction> get_result();
    /**
     *  @return the vector of inteprolated rational functions over the last field
     */
    std::vector<RationalFunctionFF> get_result_ff();
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
     *  Enables the scan for factors
     */
    void enable_factor_scan();
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
    /**
     *  Allows to abort the current reconstruction, only valid after prime changes
     */
    void abort();
    /**
     *  Allows to resume the current reconstruction, only valid after prime changes
     */
    void resume();

    enum verbosity_levels {SILENT, IMPORTANT, CHATTY};
    enum RatReconst_status {RECONSTRUCTING, DONE, DELETED};

  private:
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point prime_start = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point last_print_time = std::chrono::high_resolution_clock::now();
    const uint32_t n;
    const uint32_t thr_n;
    const uint32_t bunch_size = 1;
    uint32_t prime_it = 0;
    uint32_t total_iterations = 0;
    uint32_t iteration = 0;
    uint32_t probes_queued = 0;
    uint32_t fed_ones = 0;
    uint32_t probes_for_next_prime = 0;
    uint32_t items = 0;
    uint32_t items_done = 0;
    uint32_t items_new_prime = 0;
    uint32_t feed_jobs = 0;
    uint32_t interpolate_jobs = 0;
    uint64_t ind = 0;
    BlackBoxBase<BlackBoxTemp>& bb;
    int verbosity;
    double average_black_box_time = 0.;
    std::atomic<bool> scan = {false};
    std::atomic<bool> factor_scan = {false};
    std::atomic<bool> aborted = {false};
    std::atomic<bool> resumed = {false};
    std::atomic<bool> new_prime = {false};
    std::atomic<bool> done = {false};
    std::atomic<bool> change_var_order = {false};
    static bool printed_logo;
    bool save_states = false;
    bool resume_from_state = false;
    bool safe_mode = false;
    bool one_done = false;
    bool one_new_prime = false;
    bool set_anchor_points = false;
    RatReconst_list reconst;
    std::vector<std::string> tags;
    std::vector<std::string> file_paths;
    std::string curr_var = "";
    std::vector<FFInt> rand_zi_fac {};
    std::ofstream logger;
    std::vector<uint32_t> max_degs {};
    ThreadPool tp;
    // TODO tidy up the mutexes
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
    std::vector<std::vector<std::string>> factorizations {};
    std::unordered_set<uint32_t> possible_factors_bb_counter {};
    std::unordered_map<uint32_t, std::list<RationalFunction>> factors_rf {};
    std::unordered_map<uint64_t, std::pair<FFInt, std::vector<uint32_t>>> index_map;
    std::unordered_map<std::vector<uint32_t>, uint32_t, UintHasher> started_probes;
    std::deque<std::pair<uint64_t, std::vector<FFInt>>> requested_probes;
    std::queue<std::pair<std::vector<uint64_t>, std::vector<std::vector<FFInt>>>> computed_probes;
    std::unordered_map<std::vector<uint32_t>, std::unordered_set<uint64_t>, UintHasher> chosen_t;
    std::unordered_map<uint32_t, uint32_t> optimal_var_order {}; /**< first is old position, second is new position */
    std::unordered_map<uint32_t, ShuntingYardParser> parsed_factors {};
    RatReconst tmp_rec;
    std::vector<FFInt> shift;
    /**
    *  Scan the black-box functions for a sparse shift
    */
    void scan_for_shift();
    /**
    *  Scan the black-box functions for a sparse shift
    */
    void scan_for_factors();
    /**
     *  Combines results after factor scan in one prime field
     *  @param poly polynomial in FLINT`s notation
     *  @param combined_ci the map of combined coefficients
     *  @param combined_prime previously used combined prime
     *  @return the combined prime
     */
    mpz_class combine_primes(const std::unordered_map<uint32_t, uint64_t>& poly,
                             std::unordered_map<uint32_t, mpz_class>& combined_ci,
                             const mpz_class& combined_prime);
    /**
     *  Initializes vector of reconstruction objects and starts first probes
     */
    void start_first_runs();
    /**
     *  Starts new jobs until the reconstruction is done
     *  @param prime_counter sets how many interpolations have to be performed at most
     */
    void run_until_done(uint32_t prime_counter = 100);
    /**
     *  Queues a number of probes for a given zi_order
     *  @param zi_order the order of which a given number of probes should be queued
     *  @param to_start the number of probes which should be queued
     * TODO
     */
    void queue_probes(const std::vector<uint32_t>& zi_order, const uint32_t to_start, const bool first = false);
    /**
     * Gets a probe from probes, probes_bunch, or bunch
     * @param t is set to the t of the returned probe
     * @param zi_order is set to the zi_order of the returned probe
     * @param probe is set to point to a probe from probes, probes_bunch, or bunch
     * @param time is set to the time of the returned probe
     */
    // TODO
    void get_probe(std::vector<uint64_t>& indices, std::vector<std::vector<FFInt>>& probes);
    /**
     *  Feeds the reconstruction objects
     *  @param zi_order the order at which the black box was probed
     *  @param t the value of the homogenization variable t
     *  @param probe a vector of black box probes in an immutable order
     */
    // TODO
    void feed_job(const std::vector<uint64_t>& indices, const std::vector<std::vector<FFInt>>& probes);
    /**
     *  Interpolates a RatReconst and queues new jobs if required
     *  @param rec a reference to the RatReconst
     */
    void interpolate_job(RatReconst_tuple& rec);
    /**
     *  Removes all RatReconst from reconst which are flagged as DELETE
     */
    void clean_reconst();
    /**
     *  TODO
     */
    void get_job();
    /**
     *  TODO
     */
    template<uint32_t N>
    void start_new_job(std::unique_lock<std::mutex>& lock_probe_queue);
    /**
     *  TODO
     */
    void reset_new_prime();
#if WITH_MPI
    int world_size;
    uint32_t total_thread_count = 0;
    uint32_t iterations_on_this_node = 0;
    bool proceed = false;
    bool continue_communication = false;
    std::unordered_map<int, uint64_t> nodes;
    std::queue<std::pair<int, uint64_t>> empty_nodes;
    std::condition_variable cond_val;
    std::atomic<bool> new_jobs = {false};
    /**
     *  TODO
     */
    inline void mpi_setup();
    void send_first_jobs();
    void mpi_communicate();
#endif
  };

  template<typename BlackBoxTemp>
  bool Reconstructor<BlackBoxTemp>::printed_logo = false;

  template<typename BlackBoxTemp>
  Reconstructor<BlackBoxTemp>::Reconstructor(const uint32_t n_, const uint32_t thr_n_, BlackBoxBase<BlackBoxTemp>& bb_,
#if !WITH_MPI
                               const int verbosity_): n(n_), thr_n(thr_n_), bb(bb_), verbosity(verbosity_), tp(thr_n_) {
#else
                               int verbosity_): n(n_), thr_n(thr_n_ - 1), bb(bb_), verbosity(verbosity_), tp(thr_n), total_thread_count(thr_n) {
#endif
    FFInt::set_new_prime(primes()[prime_it]);
    uint64_t seed = static_cast<uint64_t>(std::time(0));
    BaseReconst().set_seed(seed);
    tmp_rec = RatReconst(n);

    logger.open("firefly.log");

    if (verbosity > SILENT) {
      if(!printed_logo) {
        std::cout << "\nFire\033[1;32mFly\033[0m " << FireFly_VERSION_MAJOR << "."
                  << FireFly_VERSION_MINOR << "." << FireFly_VERSION_RELEASE << "\n\n";
        printed_logo = true;
      }
      INFO_MSG("Launching " << thr_n << " thread(s) with maximum bunch size 1");
      INFO_MSG("Using seed " + std::to_string(seed) + " for random numbers");
      logger << "\nFireFly " << FireFly_VERSION_MAJOR << "."
                << FireFly_VERSION_MINOR << "." << FireFly_VERSION_RELEASE << "\n\n"
                <<"Launching " << thr_n << " thread(s) with maximum bunch size 1\n"
                <<"Using seed " << std::to_string(seed) << " for random numbers\n";
    }
  }

  template<typename BlackBoxTemp>
  Reconstructor<BlackBoxTemp>::Reconstructor(const uint32_t n_, const uint32_t thr_n_, const uint32_t bunch_size_,
#if !WITH_MPI
                               BlackBoxBase<BlackBoxTemp>& bb_, const int verbosity_): n(n_), thr_n(thr_n_), bunch_size(bunch_size_), bb(bb_), verbosity(verbosity_), tp(thr_n_) {
#else
                               BlackBoxBase<BlackBoxTemp>& bb_, int verbosity_): n(n_), thr_n(thr_n_ - 1), bunch_size(bunch_size_), bb(bb_), verbosity(verbosity_), tp(thr_n), total_thread_count(thr_n) {
#endif
    if (bunch_size != 1 && bunch_size != 2 && bunch_size != 4 && bunch_size != 8 && bunch_size != 16 && bunch_size != 32 && bunch_size != 64 && bunch_size != 128 && bunch_size != 256) {
      ERROR_MSG("Maximum bunch size " + std::to_string(bunch_size) + " is no supported power of 2!\nChoose among 1, 2, 4, 8, 16, 32, 64, 128, 256");
      logger << "Maximum bunch size " << std::to_string(bunch_size) << " is no supported power of 2!\nChoose among 1, 2, 4, 8, 16, 32, 64, 128, 256\n";
      std::exit(EXIT_FAILURE);
    }

    FFInt::set_new_prime(primes()[prime_it]);
    uint64_t seed = static_cast<uint64_t>(std::time(0));
    BaseReconst().set_seed(seed);
    tmp_rec = RatReconst(n);

    logger.open("firefly.log");

    if (verbosity > SILENT) {
            if(!printed_logo) {
        std::cout << "\nFire\033[1;32mFly\033[0m " << FireFly_VERSION_MAJOR << "."
                  << FireFly_VERSION_MINOR << "." << FireFly_VERSION_RELEASE << "\n\n";
        printed_logo = true;
      }
      INFO_MSG("Launching " << thr_n << " thread(s) with maximum bunch size " + std::to_string(bunch_size_));
      INFO_MSG("Using seed " + std::to_string(seed) + " for random numbers");
      logger << "\nFireFly " << FireFly_VERSION_MAJOR << "."
                << FireFly_VERSION_MINOR << "." << FireFly_VERSION_RELEASE << "\n\n"
                <<"Launching " << thr_n << " thread(s) with maximum bunch size " << std::to_string(bunch_size_) <<  "\n"
                <<"Using seed " + std::to_string(seed) + " for random numbers\n";
    }
  }

  template<typename BlackBoxTemp>
  Reconstructor<BlackBoxTemp>::~Reconstructor() {
    logger.close();
    tp.kill_all();

    auto it = reconst.begin();

    while (it != reconst.end()) {
      if (std::get<2>(*it) == DELETED) {
        // delete mutex
        delete std::get<1>(*it);

        // remove from list
        it = reconst.erase(it);
      } else {
        // delete mutex
        delete std::get<1>(*it);

        // delete RatReconst
        delete std::get<3>(*it);

        // remove from list
        it = reconst.erase(it);
      }
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::enable_scan() {
    if (n == 1) {
      WARNING_MSG("Scan disabled for a univariate rational function.");
      logger << "Scan disabled for a univariate rational function.\n";
    } else {
      scan = true;
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::enable_factor_scan() {
#ifndef FLINT
    ERROR_MSG("FireFly is not compiled with FLINT. No polynomial factoring possible!");
#endif
    factor_scan = true;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::set_tags() {
    save_states = true;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::set_tags(const std::vector<std::string>& tags_) {
    save_states = true;
    tags = tags_;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::resume_from_saved_state() {
    tinydir_dir dir;
    tinydir_open_sorted(&dir, "ff_save/states");

    std::vector<std::string> files;
    std::vector<std::string> paths;

    for (size_t i = 0; i != dir.n_files; ++i) {
      tinydir_file file;
      tinydir_readfile_n(&dir, &file, i);

      if (!file.is_dir) {
        files.emplace_back(file.name);
      }
    }

    tinydir_close(&dir);

    std::sort(files.begin(), files.end(), [](const std::string & l, const std::string & r) {
      return std::stoi(l.substr(0, l.find("_"))) < std::stoi(r.substr(0, r.find("_")));
    });

    for (const auto & file : files) {
      paths.emplace_back("ff_save/states/" + file);
    }

    if (paths.size() != 0) {
      resume_from_saved_state(paths);
    } else {
      save_states = true;
      WARNING_MSG("Directory './ff_save' does not exist or has no content");
      logger << "Directory './ff_save' does not exist or has no content\n";
      INFO_MSG("Starting new reconstruction and saving states");
      logger << "Starting new reconstruction and saving states\n";
      return;
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::resume_from_saved_state(const std::vector<std::string>& file_paths_) {
    if (verbosity > SILENT) {
      INFO_MSG("Loading saved states");
      logger << "Loading saved states\n";
    }

    set_anchor_points = false;

    std::ifstream v_file;
    igzstream validation_file;
    validation_file.open("ff_save/validation.gz");
    v_file.open("ff_save/validation.gz");
    std::string line;

    if (v_file.is_open()) {
      std::getline(validation_file, line);
      std::vector<FFInt> values = parse_vector_FFInt(line);

      std::vector<FFInt> result = bb.eval(values);
      size_t counter = 0;

      while (std::getline(validation_file, line)) {
        if (std::stoul(line) != result[counter]) {
          ERROR_MSG("Validation failed: Entry " + std::to_string(counter) + " does not match the black-box result!");
          logger << "Validation failed: Entry " + std::to_string(counter) + " does not match the black-box result!\n";
          std::exit(EXIT_FAILURE);
        }

        ++counter;
      }

      if (counter != result.size()) {
        ERROR_MSG("Validation failed: Number of entries does not match the black box!");
        logger << "Validation failed: Number of entries does not match the black box!\n";
        std::exit(EXIT_FAILURE);
      }
    } else {
      ERROR_MSG("Validation file not found!");
      logger << "Validation file not found!\n";
      std::exit(EXIT_FAILURE);
    }

    v_file.close();
    validation_file.close();

    save_states = true;
    resume_from_state = true;
    file_paths = file_paths_;
    items = static_cast<uint32_t>(file_paths.size());
    prime_it = 200; // increase so that the minimum is the mininmum of the files

    for (uint32_t i = 0; i != items; ++i) {
      prime_it = std::min(prime_it, parse_prime_number(file_paths[i]));
    }

    FFInt::set_new_prime(primes()[prime_it]);

    tmp_rec.start_from_saved_file(file_paths[0]);

    // Get probe files
    tinydir_dir dir;
    tinydir_open_sorted(&dir, "ff_save/probes");

    std::vector<std::string> files;
    std::vector<std::string> probe_files;

    for (size_t i = 0; i != dir.n_files; ++i) {
      tinydir_file file;
      tinydir_readfile_n(&dir, &file, i);

      if (!file.is_dir) {
        files.emplace_back(file.name);
      }
    }

    tinydir_close(&dir);

    std::sort(files.begin(), files.end(), [](const std::string & l, const std::string & r) {
      return std::stoi(l.substr(0, l.find("_"))) < std::stoi(r.substr(0, r.find("_")));
    });

    for (const auto & file : files) {
      probe_files.emplace_back("ff_save/probes/" + file);
    }

    if (probe_files.size() != items) {
      ERROR_MSG("Mismatch in number of probe files");
      logger << "Mismatch in number of probe files\n";
      std::exit(EXIT_FAILURE);
    }

    for (uint32_t i = 0; i != items; ++i) {
      RatReconst* rec = new RatReconst(n);
      /*std::pair<bool, uint32_t> shift_prime = */rec->start_from_saved_file(file_paths[i]);

      // Fill in already used ts from prior run
      auto tmp = rec->read_in_probes(probe_files[i]);

      for (const auto & already_chosen_t : tmp) {
        if (chosen_t.find(already_chosen_t.first) != chosen_t.end()) {
          auto set = chosen_t[already_chosen_t.first];

          for (const auto & ts : already_chosen_t.second) {
            if (set.find(ts) == set.end()) {
              set.emplace(ts);
            }
          }

          chosen_t[already_chosen_t.first] = set;
        } else {
          chosen_t.emplace(std::make_pair(already_chosen_t.first, already_chosen_t.second));
        }
      }

      if (safe_mode)
        rec->set_safe_interpolation();

      rec->set_tag(std::to_string(i));

      if (rec->is_done()) {
        ++items_done;
        std::mutex* mut = new std::mutex;

        reconst.emplace_back(std::make_tuple(i, mut, DONE, rec));
      } else {
        if (rec->is_new_prime()) {
          probes_for_next_prime = std::max(probes_for_next_prime, rec->get_num_eqn());
          ++items_new_prime;

          if (rec->get_prime() != prime_it + 1) {
            set_anchor_points = true;
          }
        }

        std::mutex* mut = new std::mutex;

        reconst.emplace_back(std::make_tuple(i, mut, RECONSTRUCTING, rec));
      }
    }

    if (prime_it == 0 && items != items_new_prime + items_done) {
      set_anchor_points = false;
      std::ifstream anchor_point_file;
      anchor_point_file.open("ff_save/anchor_points");

      if (anchor_point_file.is_open()) {
        std::getline(anchor_point_file, line);
        tmp_rec.set_anchor_points(parse_vector_FFInt(line, static_cast<int>(n)));
      } else {
        ERROR_MSG("Anchor point file not found!");
        logger << "Anchor point file not found!\n";
        std::exit(EXIT_FAILURE);
      }

      anchor_point_file.close();

      std::ifstream shift_file;
      shift_file.open("ff_save/shift");

      if (shift_file.is_open()) {
        std::getline(shift_file, line);
        tmp_rec.set_shift(parse_vector_FFInt(line, static_cast<int>(n)));
        shift = tmp_rec.get_zi_shift_vec();
      } else {
        ERROR_MSG("Shift file not found!");
        logger << "Shift file not found!\n";
        std::exit(EXIT_FAILURE);
      }
    }

    if (safe_mode) {
      scan = false;
    }

    if (scan) {
      if (prime_it == 0 && items_new_prime != items) {
        std::ifstream file;
        file.open("ff_save/scan");

        if (file.is_open()) {
          scan = false;
        } else {
          ERROR_MSG("Cannot resume from saved state because the scan was not completed.");
          ERROR_MSG("Please remove the directory 'ff_save' and start from the beginning.");
          logger << "Cannot resume from saved state because the scan was not completed.\n";
          logger << "Please remove the directory 'ff_save' and start from the beginning.\n";
          std::exit(EXIT_FAILURE);
        }

        file.close();
      } else {
        scan = false;
      }
    }

    logger << "All files loaded | Done: " << std::to_string(items_done) << " / " << std::to_string(items) <<
      " | " << "Needs new prime field: " << std::to_string(items_new_prime) << " / " << std::to_string(items - items_done) << "\n";

    if (verbosity > SILENT) {
      INFO_MSG("All files loaded | Done: " + std::to_string(items_done) + " / " + std::to_string(items) +
               " | " + "Needs new prime field: " + std::to_string(items_new_prime) + " / " + std::to_string(items - items_done));
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::set_safe_interpolation() {
    safe_mode = true;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::reconstruct(uint32_t prime_counter) {
    start = std::chrono::high_resolution_clock::now();

    if (!aborted || resumed)
      done = false;

    if (resumed)
      resumed = false;

    if (aborted)
      aborted = false;

#if WITH_MPI
    ThreadPool tp_comm(1);
    tp_comm.run_priority_task([this]() {
      {
        std::unique_lock<std::mutex> lock(mutex_probe_queue);

        while (!continue_communication) {
          cond_val.wait(lock);
        }

        continue_communication = false;
        proceed = true;

        cond_val.notify_one();
      }

      mpi_communicate();
    });
#endif

    if (!resume_from_state) {
      logger << "\n" << "Promote to new prime field: F(" << std::to_string(primes()[prime_it]) << ")\n";

      if (verbosity > SILENT) {
        std::cout << "\n";
        INFO_MSG("Promote to new prime field: F(" + std::to_string(primes()[prime_it]) + ")");
      }

      if (safe_mode) {
        tmp_rec.set_safe_interpolation();

        if (factor_scan) {
          WARNING_MSG("Disabled factor scan in safe mode!");
          logger << "Disabled factor scan in safe mode!\n";
          factor_scan = false;
        }

        if (scan) {
          WARNING_MSG("Disabled shift scan in safe mode!");
          logger << "Disabled shift scan in safe mode!\n";
          scan = false;
        }
      }

      if (factor_scan) {
        RatReconst::reset();
        scan_for_factors();
        tmp_rec = RatReconst(n);
      }

      if (scan) {
        scan_for_shift();
#if !WITH_MPI
        uint32_t to_start = thr_n ;//* bunch_size; // TODO
        queue_probes(std::vector<uint32_t> (n - 1, 1), to_start);
#else
        uint32_t to_start = buffer * total_thread_count; // TODO: start even more? * bunch_size
        queue_probes(std::vector<uint32_t> (n - 1, 1), to_start, true);
#endif

        started_probes.emplace(std::vector<uint32_t> (n - 1, 1), to_start);

#if WITH_MPI
        cond_val.notify_one();

        {
          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          while (!proceed) {
            cond_val.wait(lock_probe_queue);
          }

          proceed = false;
        }

        for (uint32_t j = 0; j != to_start; ++j) {
          tp.run_task([this]() {
            get_job();
          });
        }
#endif
      } else {
        start_first_runs();
      }
    } else {
      scan = false;

      if (items_done == items) {
        done = true;
      }
    }

    if (!done) {
      if (save_states && !set_anchor_points) {
        std::string tmp_str = "";
        std::ofstream file;
        file.open("ff_save/shift");
        tmp_str = "";

        for (const auto & el : tmp_rec.get_zi_shift_vec()) {
          tmp_str += std::to_string(el.n) + " ";
        }

        tmp_str += std::string("\n");
        file << tmp_str;
        file.close();
      }

      run_until_done(prime_counter);
    }
#if WITH_MPI
    else {
      mpi_setup();

      {
        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        continue_communication = true;

        cond_val.notify_one();

        while (!proceed) {
          cond_val.wait(lock_probe_queue);
        }

        proceed = false;
      }
    }
#endif

    if (one_done || one_new_prime) {
      logger << "Probe: " << std::to_string(iteration) <<
                 " | Done: " << std::to_string(items_done) << " / " << std::to_string(items) <<
                 " | " << "Needs new prime field: " << std::to_string(items_new_prime) << " / " << std::to_string(items - items_done) << "\n";
    }

    logger << "Completed reconstruction in "
      << std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count())
      << " s | " << std::to_string(total_iterations) << " probes in total\n"
      << "Needed prime fields: " << std::to_string(prime_it) << " + 1\n"
      << "Average time of the black-box probe: " << std::to_string(average_black_box_time) << " s\n";

    if (verbosity > SILENT) {
      if (one_done || one_new_prime) {
        INFO_MSG("Probe: " + std::to_string(iteration) +
                 " | Done: " + std::to_string(items_done) + " / " + std::to_string(items) +
                 " | " + "Needs new prime field: " + std::to_string(items_new_prime) + " / " + std::to_string(items - items_done) + "\n");
      }

      INFO_MSG("Completed reconstruction in " +
               std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count()) + " s | " + std::to_string(total_iterations) + " probes in total");
      INFO_MSG("Needed prime fields: " + std::to_string(prime_it) + " + 1");
      INFO_MSG("Average time of the black-box probe: " + std::to_string(average_black_box_time) + " s");
    }
  }

  template<typename BlackBoxTemp>
  std::vector<RationalFunction> Reconstructor<BlackBoxTemp>::get_result() {
    std::vector<RationalFunction> result {};
    uint32_t counter = 0;

    for (auto & rec : reconst) {
      if (std::get<2>(rec) == DONE) {
        result.emplace_back(std::get<3>(rec)->get_result());

        if (change_var_order) {
          result.back().set_var_order(optimal_var_order);
        }

        if (factors_rf.find(counter) != factors_rf.end()) {
          for (const auto& factor : factors_rf[counter]) {
            result.back().add_factor(factor);
          }
        }

        ++counter;
      }
    }

    return result;
  }

  template<typename BlackBoxTemp>
  std::vector<RationalFunctionFF> Reconstructor<BlackBoxTemp>::get_result_ff() {
    std::vector<RationalFunctionFF> result {};

    for (auto & rec : reconst) {
      result.emplace_back(std::get<3>(rec)->get_result_ff());
    }

    return result;
  }

  template<typename BlackBoxTemp>
  std::vector<std::pair<std::string, RationalFunction>> Reconstructor<BlackBoxTemp>::get_early_results() {
    if (scan) {
      return std::vector<std::pair<std::string, RationalFunction>> {};
    }

    std::unique_lock<std::mutex> lock_clean(clean);

    std::vector<std::pair<std::string, RationalFunction>> result;

    for (auto & rec : reconst) {
      std::unique_lock<std::mutex> lock_exists(*(std::get<1>(rec)));

      if (std::get<2>(rec) == DONE) {
        if (save_states) {
          result.emplace_back(std::make_pair(std::get<3>(rec)->get_tag_name(), std::get<3>(rec)->get_result()));
        } else {
          result.emplace_back(std::make_pair(std::to_string(std::get<0>(rec)), std::get<3>(rec)->get_result()));
        }

        std::get<2>(rec) = DELETED;
        delete std::get<3>(rec);
      }
    }

    return result;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::scan_for_shift() {
    logger << "Scanning for a sparse shift\n";

    if (verbosity > SILENT)
      INFO_MSG("Scanning for a sparse shift");

    // Generate all possible combinations of shifting variables
    const auto shift_vec = generate_possible_shifts(n);

    bool first = true;
    bool found_shift = false;
    uint32_t counter = 0;
    uint32_t bound = static_cast<uint32_t>(shift_vec.size());

    tmp_rec.scan_for_sparsest_shift();

    start_first_runs();

    uint32_t max_deg_num = 0;
    uint32_t max_deg_den = 0;

    // Run this loop until a proper shift is found
    while (!found_shift && counter != bound) {
      if (!first) {
        tmp_rec.set_zi_shift(shift_vec[counter]);
        shift = tmp_rec.get_zi_shift_vec();

#if !WITH_MPI
        uint32_t to_start = thr_n ;//* bunch_size; // TODO
        queue_probes(std::vector<uint32_t> (n - 1, 1), to_start);
#else
        uint32_t to_start = buffer * total_thread_count; // TODO: start even more? * bunch_size
        queue_probes(std::vector<uint32_t> (n - 1, 1), to_start, true);
#endif

        started_probes.emplace(std::vector<uint32_t> (n - 1, 1), to_start);

#if WITH_MPI
        cond_val.notify_one();

        {
          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          while (!proceed) {
            cond_val.wait(lock_probe_queue);
          }

          proceed = false;
        }

        for (uint32_t j = 0; j != to_start; ++j) {
          tp.run_task([this]() {
            get_job();
          });
        }
#endif
      }

      run_until_done();

      found_shift = true;

      for (auto & rec : reconst) {
        std::get<2>(rec) = RECONSTRUCTING;

        if (!(std::get<3>(rec)->is_shift_working())) {
          found_shift = false;
        }

        if (first) {
          std::pair<uint32_t, uint32_t> degs = std::get<3>(rec)->get_max_deg();

          max_deg_num = std::max(max_deg_num, degs.first);
          max_deg_den = std::max(max_deg_den, degs.second);
        }
      }

      if (first) {
        found_shift = false;
        first = false;

        logger << "Maximum degree of numerator: " << std::to_string(max_deg_num)
          << " | Maximum degree of denominator: " << std::to_string(max_deg_den) << "\n";
        logger.close();
        logger.open("firefly.log", std::ios_base::app);

        if (verbosity > SILENT) {
          INFO_MSG("Maximum degree of numerator: " + std::to_string(max_deg_num) + " | Maximum degree of denominator: " + std::to_string(max_deg_den));
        }
      } else {
        ++counter;
      }

      reset_new_prime();
      items_done = 0;
      done = false;
    }

    if (found_shift) {
      tmp_rec.set_zi_shift(shift_vec[counter - 1]);
    } else {
      tmp_rec.set_zi_shift(std::vector<uint32_t> (n, 1));
    }

    shift = tmp_rec.get_zi_shift_vec();

    for (auto & rec : reconst) {
      std::get<2>(rec) = RECONSTRUCTING;
      std::get<3>(rec)->accept_shift();
    }

    scan = false;

    if (save_states == true) {
      std::ofstream file;
      file.open("ff_save/scan");
      file.close();
    }

    if (found_shift) {
      std::string msg = "";

      for (const auto & el : shift_vec[counter - 1]) {
        msg += std::to_string(el) + ", ";
      }

      msg = msg.substr(0, msg.size() - 2);

      logger << "Found a sparse shift after " << std::to_string(counter + 1) << " scans\n"
        << "Shift the variable tuple (" << msg << ")\n";

      if (verbosity > SILENT) {
        INFO_MSG("Found a sparse shift after " + std::to_string(counter + 1) + " scans");
        INFO_MSG("Shift the variable tuple (" + msg + ")");
      }
    } else {
      logger << "Found no sparse shift after " << std::to_string(counter + 1) << " scans\n";

      if (verbosity > SILENT)
        INFO_MSG("Found no sparse shift after " + std::to_string(counter + 1) + " scans");
    }

    logger << "Completed scan in "
      << std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - prime_start).count())
      << " s | " << std::to_string(total_iterations) << " probes in total\n";

    logger << "Average time of the black-box probe: " << std::to_string(average_black_box_time) << " s\n\n";
    logger << "Proceeding with interpolation over prime field F(" << std::to_string(primes()[prime_it]) << ")\n";
    logger.close();
    logger.open("firefly.log", std::ios_base::app);

    if (verbosity > SILENT) {
      INFO_MSG("Completed scan in " + std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - prime_start).count()) +
               " s | " + std::to_string(total_iterations) + " probes in total");
      INFO_MSG("Average time of the black-box probe: " + std::to_string(average_black_box_time) + " s\n");
      INFO_MSG("Proceeding with interpolation over prime field F(" + std::to_string(primes()[prime_it]) + ")");
    }

    prime_start = std::chrono::high_resolution_clock::now();
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::scan_for_factors() {
#ifdef FLINT
    auto clock_1 = std::chrono::high_resolution_clock::now();
    factorizations.reserve(n);
    int old_verbosity = verbosity;
    verbosity = SILENT;

    if (old_verbosity > SILENT) {
      INFO_MSG("Scanning for factors");
    }

    logger << "Scanning for factors\n";

    uint32_t total_number_of_factors = 0;
    uint32_t number_of_factors = 0;
    max_degs = std::vector<uint32_t> (n, 0);
    shift = std::vector<FFInt> (n, 0);
    rand_zi_fac = std::vector<FFInt> (n, 0);
    std::unordered_map<uint32_t, std::pair<std::unordered_set<std::string>, std::unordered_set<std::string>>> factors_str {};

    // Run this loop until a proper shift is found
    for (int i = 0; i != n; ++i) {
      std::unordered_map<uint32_t, std::unordered_map<uint32_t, mpz_class>> combined_ni {};
      std::unordered_map<uint32_t, std::unordered_map<uint32_t, mpz_class>> combined_di {};
      std::unordered_map<uint32_t, std::pair<uint32_t, uint32_t>> max_deg_map {};
      mpz_class combined_prime = FFInt::p;
      bool fac_done = false;
      curr_var = "x" + std::to_string(i + 1);
      possible_factors_bb_counter.clear();
      size_t tmp_prime_it = 0;

      while (!fac_done) {
        std::unordered_map<uint32_t,std::pair<std::unordered_set<std::string>, std::unordered_set<std::string>>> possible_factors {};
        mpz_class tmp_combined_prime = combined_prime;

        for (int scan_n = 0; scan_n != 2; ++scan_n) {
          number_of_factors = 0;

          for (int j = 0; j != n; ++j) {
            if (j == i) {
              rand_zi_fac[j] = 1;
            } else {
              rand_zi_fac[j] = tmp_rec.get_rand_64();
            }
          }

          start_first_runs();
          run_until_done();
          prime_it = 0;

          uint32_t counter = 0;

          // Get factors
          for (auto& rec : reconst) {
            std::unique_lock<std::mutex> lock_exists(*(std::get<1>(rec)));

            if (scan_n == 0) {
              if (possible_factors_bb_counter.find(counter) != possible_factors_bb_counter.end()) {
                auto tmp_fac_num_den = std::get<3>(rec)->get_factors_ff();
                std::unordered_set<std::string> tmp_fac_num {};
                std::unordered_set<std::string> tmp_fac_den {};
                tmp_fac_num.insert(tmp_fac_num_den.first.begin(), tmp_fac_num_den.first.end());
                tmp_fac_den.insert(tmp_fac_num_den.second.begin(), tmp_fac_num_den.second.end());
                possible_factors.emplace(counter, std::make_pair(tmp_fac_num, tmp_fac_den));
              }

              if (tmp_prime_it == 0) {
                auto tmp_fac_num_den = std::get<3>(rec)->get_factors_ff();
                std::unordered_set<std::string> tmp_fac_num {};
                std::unordered_set<std::string> tmp_fac_den {};
                tmp_fac_num.insert(tmp_fac_num_den.first.begin(), tmp_fac_num_den.first.end());
                tmp_fac_den.insert(tmp_fac_num_den.second.begin(), tmp_fac_num_den.second.end());
                possible_factors.emplace(counter, std::make_pair(tmp_fac_num, tmp_fac_den));

                // Get maximum degrees
                auto tmp_max_degs = std::get<3>(rec)->get_max_deg();
                max_deg_map.emplace(std::make_pair(counter, tmp_max_degs));

                uint32_t max_val = std::max(tmp_max_degs.first, tmp_max_degs.second);
                if (max_val > max_degs[i]) {
                  max_degs[i] = max_val;
                }

                uint32_t tmp_n_fac = tmp_fac_num_den.first.size() + tmp_fac_num_den.second.size();
                number_of_factors += tmp_n_fac;

                if (tmp_n_fac != 0) {
                  possible_factors_bb_counter.emplace(counter);
                }
              }
            } else if (scan_n == 1 && possible_factors_bb_counter.find(counter) != possible_factors_bb_counter.end()) {
              // Rewrite and store in result objects. Compare to previous factorizations
              std::unordered_set<uint32_t> fac_nums_c {};
              std::unordered_set<uint32_t> fac_dens_c {};

              auto tmp_factors = std::get<3>(rec)->get_factors_ff();

              // Numerator
              uint32_t fac_counter = 0;
              for (const auto& tmp_factor : tmp_factors.first) {
                bool found = possible_factors[counter].first.find(tmp_factor) != possible_factors[counter].first.end();

                if (possible_factors[counter].first.find(tmp_factor) != possible_factors[counter].first.end()) {
                  fac_nums_c.emplace(fac_counter);
                  ++number_of_factors;
                }

                ++fac_counter;
              }

              fac_counter = 0;
              // Denominator
              for (const auto& tmp_factor : tmp_factors.second) {
                if (possible_factors[counter].second.find(tmp_factor) != possible_factors[counter].second.end()) {
                  fac_dens_c.emplace(fac_counter);
                  ++number_of_factors;
                }

                ++fac_counter;
              }

              if (fac_nums_c.empty() && fac_dens_c.empty()) {
                possible_factors_bb_counter.erase(counter);
              } else {
                auto canonical_factors = std::get<3>(rec)->get_canonical_factors(std::make_pair(fac_nums_c, fac_dens_c));

                if (tmp_prime_it == 0) {
                  // Init first combinations
                  std::unordered_map<uint32_t, mpz_class> tmp_combined_ni {};
                  std::unordered_map<uint32_t, mpz_class> tmp_combined_di {};
                  uint32_t fac_max_deg_num = 0, fac_max_deg_den = 0;

                  for (const auto& mon : canonical_factors.first) {
                    tmp_combined_ni.emplace(std::make_pair(mon.first, mon.second));
                    fac_max_deg_num = std::max(fac_max_deg_num, mon.first);
                  }

                  for (const auto& mon : canonical_factors.second) {
                    tmp_combined_di.emplace(std::make_pair(mon.first, mon.second));
                    fac_max_deg_den = std::max(fac_max_deg_den, mon.first);
                  }

                  max_deg_map[counter].first -= fac_max_deg_num;
                  max_deg_map[counter].second -= fac_max_deg_den;

                  combined_ni[counter] = tmp_combined_ni;
                  combined_di[counter] = tmp_combined_di;
                } else {
                  // First, check if done, else combine results
                  bool run_test = true;
                  bool combine_results = false;
                  std::unordered_map<std::vector<uint32_t>, RationalNumber, UintHasher> tmp_gni {};
                  std::unordered_map<std::vector<uint32_t>, RationalNumber, UintHasher> tmp_gdi {};

                  // Reconstruct numerator
                  for (const auto& ci : combined_ni[counter]) {
                    mpz_class a = ci.second;
                    auto res = get_rational_coef(a, combined_prime);

                    if (res.first) {
                      tmp_gni.emplace(std::make_pair(std::vector<uint32_t> (1, ci.first), res.second));
                    } else {
                      run_test = false;
                      break;
                    }
                  }

                  // Reconstruct denominator
                  if (run_test) {
                    for (const auto& ci : combined_di[counter]) {
                      mpz_class a = ci.second;
                      auto res = get_rational_coef(a, combined_prime);

                      if (res.first) {
                        tmp_gdi.emplace(std::make_pair(std::vector<uint32_t> (1, ci.first), res.second));
                      } else {
                        run_test = false;
                        combine_results = true;
                        break;
                      }
                    }
                  } else {
                    combine_results = true;
                  }

                  if (run_test) {
                    FFInt tmp_rand = tmp_rec.get_rand_64();

                    if (!tmp_gni.empty()) {
                      ff_map gi_ffi;
                      FFInt num = 0;

                      for (const auto& g_i : tmp_gni) {
                        FFInt n(g_i.second.numerator);
                        FFInt d(g_i.second.denominator);
                        gi_ffi.emplace(std::make_pair(g_i.first, n / d));
                      }

                      for (const auto& coeff : canonical_factors.first) {
                        num += coeff.second*tmp_rand.pow(coeff.first);
                      }

                      combine_results = !(PolynomialFF(1, gi_ffi).calc({tmp_rand}) == num);
                    }

                    if (!combine_results && !tmp_gdi.empty()) {
                      ff_map gi_ffi;
                      FFInt num = 0;

                      for (const auto& g_i : tmp_gdi) {
                        FFInt n(g_i.second.numerator);
                        FFInt d(g_i.second.denominator);
                        gi_ffi.emplace(std::make_pair(g_i.first, n / d));
                      }

                      for (const auto& coeff : canonical_factors.second) {
                        num += coeff.second*tmp_rand.pow(coeff.first);
                      }

                      combine_results = !(PolynomialFF(1, gi_ffi).calc({tmp_rand}) == num);
                    }
                  }

                  // combine results
                  if (combine_results) {
                    // numerator
                    if (!canonical_factors.first.empty()) {
                      combined_prime = combine_primes(canonical_factors.first, combined_ni[counter], tmp_combined_prime);
                    }

                    // denominator
                    if (!canonical_factors.second.empty()) {
                      combined_prime = combine_primes(canonical_factors.second, combined_di[counter], tmp_combined_prime);
                    }
                  } else {
                    possible_factors_bb_counter.erase(counter);
                    combined_ni.erase(counter);
                    combined_di.erase(counter);
                    Polynomial tmp_numerator;
                    Polynomial tmp_denominator;

                    // Rewrite to result
                    if (!tmp_gni.empty()) {
                      tmp_numerator = Polynomial(tmp_gni);
                      factors_str[counter].first.emplace(tmp_numerator.to_string({curr_var}));
                      tmp_numerator.set_var_pos(i);
                    }

                    if (!tmp_gdi.empty()) {
                      tmp_denominator = Polynomial(tmp_gdi);
                      factors_str[counter].second.emplace(tmp_denominator.to_string({curr_var}));
                      tmp_denominator.set_var_pos(i);
                    }

                    factors_rf[counter].emplace_back(RationalFunction(tmp_numerator, tmp_denominator));
                  }
                }
              }
            }

            ++counter;

            std::get<2>(rec) = DELETED;
            delete std::get<3>(rec);
          }

          uint32_t old_max_deg = max_degs[i];
          if (tmp_prime_it == 0 && scan_n == 1) {
            max_degs[i] = 0;

            for (const auto& el : max_deg_map) {
              uint32_t max_val = std::max(el.second.first, el.second.second);

              if (max_val > max_degs[i]) {
                max_degs[i] = max_val;
              }
            }
          }

          if (old_verbosity > SILENT) {
            if (tmp_prime_it == 0 && scan_n == 0) {
              INFO_MSG("Maximum degree of x" + std::to_string(i + 1)
                + ": " + std::to_string(max_degs[i]));
              if (max_degs[i] != 0) {
                INFO_MSG("Possible factors in x"
                  + std::to_string(i + 1) + ": " + std::to_string(number_of_factors));
              } else {
                INFO_MSG("No factors in x"
                  + std::to_string(i + 1) + ": " + std::to_string(number_of_factors));
              }
            } else if (tmp_prime_it == 0 && scan_n == 1) {
              INFO_MSG("Factors in x" + std::to_string(i + 1) + ": "
                + std::to_string(number_of_factors));
              if (max_degs[i] != old_max_deg) {
                INFO_MSG("Maximum degree of x" + std::to_string(i + 1) + " after factoring: "
                  + std::to_string(max_degs[i]));
              }
             INFO_MSG("Starting reconstruction of coefficients");
            }
          }

          if (tmp_prime_it == 0 && scan_n == 0) {
            logger << "Maximum degree of x" << std::to_string(i + 1) << ": "
              << std::to_string(max_degs[i]) << "\n";

            if (max_degs[i] != 0) {
              logger << "Possible factors in x"
                << std::to_string(i + 1) << ": " << std::to_string(number_of_factors)
                << "\n";
            } else {
              logger << "No factors in x"
                << std::to_string(i + 1) << ": " << std::to_string(number_of_factors)
                << "\n";
            }
          } else if (tmp_prime_it == 0 && scan_n == 1) {
            total_number_of_factors += number_of_factors;
            logger << "Factors in x" << std::to_string(i + 1) << ": "
              << std::to_string(number_of_factors) << "\n";

            if (max_degs[i] != old_max_deg) {
              logger << "Maximum degree of x" << std::to_string(i + 1) << " after factoring: "
                << std::to_string(max_degs[i]) << "\n";
            }

            logger << "Starting reconstruction of coefficients\n";
            logger.close();
            logger.open("firefly.log", std::ios_base::app);
          }

          RatReconst::reset(false);
          clean_reconst();
          reconst.clear();

          reset_new_prime();
          items_done = 0;
          done = false;

          if (possible_factors_bb_counter.empty()) {
            fac_done = true;
            break;
          }
        }

        // Promote to new prime field
        if (!fac_done) {
          ++tmp_prime_it;
          FFInt::set_new_prime(primes()[tmp_prime_it]);
          bb.prime_changed_internal();
        } else {
          if (old_verbosity > SILENT) {
            INFO_MSG("Completed factor scan in " + curr_var + " | "
              + std::to_string(total_iterations) + " probes in total\n");
          }

          logger << "Completed factor scan in " << curr_var << " | "
            << total_iterations << " probes in total\n\n";
        }
      }

      // Reset prime
      tmp_prime_it = 0;
      FFInt::set_new_prime(primes()[tmp_prime_it]);
      bb.prime_changed_internal();
    }

    if (save_states) {
      mkdir("ff_save", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      mkdir("ff_save/factors", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    std::vector<std::string> vars (n);

    // Reorder variables with regards to their maximum degree
    std::vector<uint32_t> indices (n);
    std::iota(indices.begin(), indices.end(), 0);

    std::stable_sort(indices.begin(), indices.end(),
                     [&](uint32_t i1, uint32_t i2) {return max_degs[i1] > max_degs[i2];});

    std::sort(max_degs.begin(), max_degs.end(), std::greater<uint32_t>());

    for (size_t i = 0; i != n; ++i) {
      optimal_var_order.emplace(std::make_pair(i, indices[i]));
      vars[i] = "x" + std::to_string(i + 1);
    }

    std::string var_order = "Using optimized variable order: {";

    for (size_t i = 0; i != n; ++i) {
      if (!change_var_order && optimal_var_order[i] != i) {
        change_var_order = true;
      }

      if (i != n - 1) {
        var_order += "x" + std::to_string(optimal_var_order[i] + 1) + ", ";
      } else {
        var_order += "x" + std::to_string(optimal_var_order[i] + 1) + "}";
      }
    }


    for (const auto& tmp_fac : factors_str) {
      std::string tmp_fac_s = "";

      if(!tmp_fac.second.first.empty()){
        tmp_fac_s += "(";
        for (const auto& tmp_fac_num : tmp_fac.second.first) {
          tmp_fac_s += "(" + tmp_fac_num + ")*";
        }

        tmp_fac_s.pop_back();
        tmp_fac_s += ")";
      }

      if (!tmp_fac.second.second.empty()) {
        if (!tmp_fac.second.first.empty()) {
          tmp_fac_s += "/(";
        } else {
          tmp_fac_s += "1/(";
        }

        for (const auto& tmp_fac_den : tmp_fac.second.second) {
          tmp_fac_s += "(" + tmp_fac_den + ")*";
        }

        tmp_fac_s.pop_back();

        ShuntingYardParser parser = ShuntingYardParser();
        parser.parse_function(tmp_fac_s + ")", vars);
        parser.precompute_tokens();
        parsed_factors.emplace(tmp_fac.first, parser);

        tmp_fac_s += ");\n";
      } else {
        ShuntingYardParser parser = ShuntingYardParser();
        parser.parse_function(tmp_fac_s, vars);
        parser.precompute_tokens();
        parsed_factors.emplace(tmp_fac.first, parser);

        tmp_fac_s += ";\n";
      }

      if (save_states) {
        ogzstream file;
        std::string file_name = "ff_save/factors/" + std::to_string(tmp_fac.first) + ".gz";
        file.open(file_name.c_str());
        file << tmp_fac_s;
        file.close();
      }
    }

    verbosity = old_verbosity;

    logger << "Completed factor scan in "
      << std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - clock_1).count())
      << " s | " << std::to_string(total_iterations) << " probes in total\n";

    logger << "Average time of the black-box probe: " << std::to_string(average_black_box_time) << " s\n";

    logger << var_order << "\n\n";
    logger << "Proceeding with interpolation over prime field F(" << std::to_string(primes()[prime_it]) << ")\n";
    logger.close();
    logger.open("firefly.log", std::ios_base::app);

    if (verbosity > SILENT) {
      INFO_MSG("Completed factor scan in " + std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - clock_1).count()) +
               " s | " + std::to_string(total_iterations) + " probes in total");
      INFO_MSG("Average time of the black-box probe: " + std::to_string(average_black_box_time) + " s");
      INFO_MSG(var_order + "\n");
      INFO_MSG("Proceeding with interpolation over prime field F(" + std::to_string(primes()[prime_it]) + ")");
    }

    factor_scan = false;

    prime_start = std::chrono::high_resolution_clock::now();
#endif
  }

  template<typename BlackBoxTemp>
  mpz_class Reconstructor<BlackBoxTemp>::combine_primes(const std::unordered_map<uint32_t, uint64_t>& poly,
                                                        std::unordered_map<uint32_t, mpz_class>& combined_ci,
                                                        const mpz_class& combined_prime) {
    std::pair<mpz_class, mpz_class> p1;
    std::pair<mpz_class, mpz_class> p2;
    std::pair<mpz_class, mpz_class> p3;
    std::unordered_map<uint32_t, mpz_class> tmp_coefs {};

    // Convert poly to mpz
    for (const auto& mon : poly) {
      tmp_coefs.emplace(std::make_pair(mon.first, mon.second));
    }

    for (auto it = tmp_coefs.begin(); it != tmp_coefs.end(); ++it) {
      p2 = std::make_pair(it->second, FFInt::p);
      p1 = std::make_pair(combined_ci[it->first], combined_prime);
      p3 = run_chinese_remainder(p1, p2);
      combined_ci[it->first] = p3.first;
    }

    return p3.second;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::start_first_runs() {
    prime_start = std::chrono::high_resolution_clock::now();

    std::vector<uint32_t> zi_order;

    if (!factor_scan) {
      shift = tmp_rec.get_zi_shift_vec();
      zi_order = std::vector<uint32_t>(n - 1, 1);
    } else {
      zi_order = std::vector<uint32_t>(0, 1);
    }

    uint32_t to_start = thr_n ;//* bunch_size;

#if !WITH_MPI
    queue_probes(zi_order, to_start);
#else
    queue_probes(zi_order, to_start, true);
#endif
    started_probes.emplace(zi_order, to_start);

#if WITH_MPI
    send_first_jobs();

    for (uint32_t j = 0; j != to_start; ++j) {
      tp.run_task([this]() {
        get_job();
      });
    }
#endif

    std::vector<uint64_t> indices;
    std::vector<std::vector<FFInt>> probes;

    get_probe(indices, probes);

    std::vector<FFInt> t_vec;
    t_vec.reserve(indices.size());
    std::vector<std::vector<uint32_t>> zi_order_vec;
    zi_order_vec.reserve(indices.size());

    uint32_t count_ones = 0;

    {// TODO mutex required?
      std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

      for (const auto& index : indices) {
        auto tmp = std::move(index_map[index]);
        index_map.erase(index);
        t_vec.emplace_back(tmp.first);
        zi_order_vec.emplace_back(std::move(tmp.second));

        if ((prime_it == 0 || safe_mode == true) && (zi_order_vec.back() == std::vector<uint32_t>(n - 1, 1) || factor_scan && zi_order_vec.back() == std::vector<uint32_t>(0, 1))) {
          ++count_ones;
        }
      }
    }

    if (count_ones != 0) {
      std::unique_lock<std::mutex> lock_status(job_control);

      fed_ones += count_ones;
    }

    {
      std::unique_lock<std::mutex> lock(future_control);

      if (!factor_scan)
       logger << "Time for the first black-box probe: " << std::to_string(average_black_box_time) << " s\n";

      if (verbosity > SILENT)
        INFO_MSG("Time for the first black-box probe: " + std::to_string(average_black_box_time) + " s");
    }

    items = static_cast<uint32_t>(probes.size());
    size_t tag_size = tags.size();

#if WITH_MPI
    {
      std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

      continue_communication = true;

      cond_val.notify_one();

      while (!proceed) {
        cond_val.wait(lock_probe_queue);
      }

      proceed = false;
    }
#endif

    if (tag_size != 0 && tag_size != items) {
      logger << "Number of tags does not match the black box!\n";
      ERROR_MSG("Number of tags does not match the black box!");
      std::exit(EXIT_FAILURE);
    }

    ogzstream file;

    if (!factor_scan && save_states) {
      mkdir("ff_save", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      mkdir("ff_save/states", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      mkdir("ff_save/probes", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      std::ofstream anchor_file;
      anchor_file.open("ff_save/anchor_points");

      std::string tmp_str = "";

      for (const auto & el : tmp_rec.get_anchor_points()) {
        tmp_str += std::to_string(el.n) + " ";
      }

      tmp_str += std::string("\n");
      anchor_file << tmp_str;
      anchor_file.close();

      file.open("ff_save/validation.gz");

      std::vector<FFInt> rand_zi;
      rand_zi = tmp_rec.get_rand_zi_vec(zi_order, false);

      file << (t_vec[0] + shift[0]).n << " ";

      for (uint32_t i = 1; i != n; ++i) {
        file << (rand_zi[i - 1] * t_vec[0] + shift[i]).n << " ";
      }

      file << "\n";
    }

    for (uint32_t i = 0; i != items; ++i) {
      RatReconst* rec;
      if (!factor_scan) {
        rec = new RatReconst(n);

        if (safe_mode) {
          rec->set_safe_interpolation();
        }

        if (scan) {
          rec->scan_for_sparsest_shift();
        }

        if (save_states) {
          rec->set_tag(std::to_string(i));

          if (tag_size > 0) {
            rec->set_tag_name(tags[i]);
          } else {
            rec->set_tag_name(std::to_string(i));
          }

          file << ((probes)[i][0]).n << "\n";
        }
      } else {
        rec = new RatReconst(1);
        rec->calc_factors(curr_var);

        // Remove functions that are irreducible
        if (!possible_factors_bb_counter.empty() &&  possible_factors_bb_counter.find(i) == possible_factors_bb_counter.end()) {
          rec->set_prime_to_max();
        }
      }

      rec->feed(t_vec, probes[i], zi_order_vec, prime_it);
      std::tuple<bool, bool, uint32_t> interpolated_done_prime = rec->interpolate();

      if (std::get<2>(interpolated_done_prime) > prime_it) {
        ++items_new_prime;
      }

      std::mutex* mut = new std::mutex;

      reconst.emplace_back(std::make_tuple(i, mut, RECONSTRUCTING, rec));
    }

    if (!factor_scan && save_states) {
      file.close();
      tags.clear();
    }

    if (!factor_scan) {
      logger << "Probe: 1 | Done: 0 / " << std::to_string(items)
        << " | Needs new prime field: " << std::to_string(items_new_prime)
        << " / " << std::to_string(items) << "\n";
      logger.close();
      logger.open("firefly.log", std::ios_base::app);
    }

    if (verbosity > SILENT) {
      INFO_MSG("Probe: 1 | Done: 0 / " + std::to_string(items) + " | Needs new prime field: " + std::to_string(items_new_prime) + " / " + std::to_string(items));
    }

    queue_probes(zi_order, static_cast<uint32_t>(probes.front().size()));
    started_probes[zi_order] += static_cast<uint32_t>(probes.front().size());
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::run_until_done(uint32_t prime_counter) {
    std::vector<uint32_t> zi_order;
    if (!factor_scan) {
      zi_order = std::vector<uint32_t> (n - 1, 1);
    } else {
      zi_order = std::vector<uint32_t> (0, 1);
    }

    new_prime = false;

#if WITH_MPI
    bool mpi_first_send = false;
#endif

    if (resume_from_state) {
      if (prime_it == 0 && items != items_new_prime + items_done) {
        logger << "Resuming in prime field: F(" << std::to_string(primes()[prime_it]) << ")\n";
        INFO_MSG("Resuming in prime field: F(" + std::to_string(primes()[prime_it]) + ")");

        {
          std::unique_lock<std::mutex> lock_status(feed_control);

          interpolate_jobs += items;
        }

        uint32_t counter = 0;

        for (auto & rec : reconst) {
          if (std::get<3>(rec)->get_prime() == 0) {
            ++counter;

            tp.run_priority_task([this, &rec]() {
              interpolate_job(rec);
            });
          }
        }

        {
          std::unique_lock<std::mutex> lock_status(feed_control);

          interpolate_jobs -= (items - counter);
        }

#if !WITH_MPI
        // TODO don't start as much ones
        queue_probes(zi_order, thr_n );//* bunch_size);
        started_probes.emplace(zi_order, thr_n );//* bunch_size);
        /*shift = tmp_rec.get_zi_shift_vec();

        uint32_t to_start = thr_n * bunch_size;
        queue_probes(std::vector<uint32_t>(n - 1, 1), to_start);
        started_probes.emplace(std::vector<uint32_t>(n - 1, 1), to_start);*/
#else
        // TODO optimize
        send_first_jobs();

        for (uint32_t j = 0; j != thr_n; ++j) {
          tp.run_task([this]() {
            get_job();
          });
        }

        {
          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          continue_communication = true;

          cond_val.notify_one();

          while (!proceed) {
            cond_val.wait(lock_probe_queue);
          }

          proceed = false;
        }

        uint32_t to_start = buffer * thr_n ;//* bunch_size; // TODO

        queue_probes(zi_order, to_start);
        started_probes[zi_order] += to_start;

#endif
      } else {
        new_prime = true;
#if WITH_MPI
        proceed = true;
        mpi_first_send = true;
        //send_first_jobs(); // TODO sends useless jobs to prepare all variables

        //{
        //  std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        //  new_jobs = true;

        //  cond_val.notify_one();

        //  while (!proceed) {
        //    cond_val.wait(lock_probe_queue);
        //  }

        //  proceed = false;
        //} // TODO end useless
#endif
      }
    }

    while (!done) {
      if (new_prime) {
        //if (prime_it == 1) exit(-1);
#if WITH_MPI
        {
          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          requested_probes = std::deque<std::pair<uint64_t, std::vector<FFInt>>>();
          new_jobs = false;

          while (!proceed) {
            cond_val.wait(lock_probe_queue);
          }

          proceed = false;
        }
#endif

        tp.kill_all();

        clean_reconst();

        if (!factor_scan && save_states) {
          for (uint32_t item = 0; item != items; ++item) {
            std::string probes_file_old = "ff_save/probes/" + std::to_string(item) + "_" + std::to_string(prime_it) + ".gz";
            std::string probes_file_new = "ff_save/probes/" + std::to_string(item) + "_" + std::to_string(prime_it + 1) + ".gz";
            std::rename(probes_file_old.c_str(), probes_file_new.c_str());

            if (prime_it) {
              std::string file_name_old = "ff_save/states/" + std::to_string(item) + "_" + std::to_string(prime_it - 1) + ".gz";
              std::string file_name_new = "ff_save/states/" + std::to_string(item) + "_" + std::to_string(prime_it) + ".gz";
              std::rename(file_name_old.c_str(), file_name_new.c_str());
            }
          }
        }

#if WITH_MPI
        {
          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          continue_communication = true;
          cond_val.notify_one();

          while (!proceed) {
            cond_val.wait(lock_probe_queue);
          }

          proceed = false;
        }
#endif

        total_iterations += iteration;
        ++prime_it;

        if(prime_it >= prime_counter || factor_scan) {
          done = true;
        }

        {
          std::unique_lock<std::mutex> lock_print(print_control);

          if ((one_done || one_new_prime) && !factor_scan) {
            logger << "Probe: " << std::to_string(iteration)
            << " | Done: " << std::to_string(items_done)
            << " / " + std::to_string(items) << " | " << "Needs new prime field: "
            << std::to_string(items_new_prime) << " / " << std::to_string(items - items_done) << "\n";

            if (verbosity > SILENT) {
            INFO_MSG("Probe: " + std::to_string(iteration) +
                     " | Done: " + std::to_string(items_done) + " / " + std::to_string(items) +
                     " | " + "Needs new prime field: " + std::to_string(items_new_prime) + " / " + std::to_string(items - items_done));
            }
          }

          if (!factor_scan) {
            logger << "Completed current prime field in "
            << std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - prime_start).count())
            << " s | " + std::to_string(total_iterations) << " probes in total\n"
            << "Average time of the black-box probe: " << std::to_string(average_black_box_time) + " s\n\n"
            << "Promote to new prime field: F(" << std::to_string(primes()[prime_it]) << ")\n";

            if (verbosity > SILENT) {
              INFO_MSG("Completed current prime field in " +
                     std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - prime_start).count()) +
                     " s | " + std::to_string(total_iterations) + " probes in total");
              INFO_MSG("Average time of the black-box probe: " + std::to_string(average_black_box_time) + " s\n");
              INFO_MSG("Promote to new prime field: F(" + std::to_string(primes()[prime_it]) + ")");
            }
          }
        }

        prime_start = std::chrono::high_resolution_clock::now();

        reset_new_prime();

        if (factor_scan) {
          break;
        }

        FFInt::set_new_prime(primes()[prime_it]);

        bb.prime_changed_internal();

        if (!parsed_factors.empty()) {
          for (auto& el : parsed_factors) {
            el.second.precompute_tokens();
          }
        }

        // if only a small constant is reconstructed it will not ask for new run
        if (probes_for_next_prime == 0) {
#if !WITH_MPI
          probes_for_next_prime = thr_n ;//* bunch_size;
#else
          probes_for_next_prime = buffer * total_thread_count; // TODO start even more? * bunch_size
#endif
        }

        if (!safe_mode && (!save_states || (save_states && !set_anchor_points)) && !tmp_rec.need_shift(prime_it)) {
          if (tmp_rec.get_zi_shift_vec() != std::vector<FFInt> (n, 0)) {
            logger << "Disable shift\n";

            if (verbosity > SILENT)
              INFO_MSG("Disable shift");

            tmp_rec.disable_shift();
          }
        }

        // Set anchor points and the shift to resume from saved probes
        if (save_states && set_anchor_points) {
          set_anchor_points = false;
          std::string line;
          std::ifstream anchor_point_file;
          anchor_point_file.open("ff_save/anchor_points");

          if (anchor_point_file.is_open()) {
            std::getline(anchor_point_file, line);
            tmp_rec.set_anchor_points(parse_vector_FFInt(line, static_cast<int>(n)));
          } else {
            logger << "Anchor point file not found!\n";
            ERROR_MSG("Anchor point file not found!");
            std::exit(EXIT_FAILURE);
          }

          anchor_point_file.close();

          std::ifstream shift_file;
          shift_file.open("ff_save/shift");

          if (shift_file.is_open()) {
            std::getline(shift_file, line);
            tmp_rec.set_shift(parse_vector_FFInt(line, static_cast<int>(n)));
          } else {
            logger << "Shift file not found!\n";
            ERROR_MSG("Shift file not found!");
            std::exit(EXIT_FAILURE);
          }

          for (auto & rec : reconst) {
            if (!std::get<3>(rec)->is_done() && std::get<3>(rec)->get_prime() > prime_it) {
              ++items_new_prime;
            }
          }
        } else {
          tmp_rec.generate_anchor_points();
        }

        shift = tmp_rec.get_zi_shift_vec();

        if (save_states) {
          std::remove("ff_save/anchor_points");
          std::ofstream file;
          file.open("ff_save/anchor_points");
          std::string tmp_str = "";

          for (const auto & el : tmp_rec.get_anchor_points()) {
            tmp_str += std::to_string(el.n) + " ";
          }

          tmp_str += std::string("\n");
          file << tmp_str;
          file.close();

          std::remove("ff_save/shift");
          file.open("ff_save/shift");
          tmp_str = "";

          for (const auto & el : tmp_rec.get_zi_shift_vec()) {
            tmp_str += std::to_string(el.n) + " ";
          }

          tmp_str += std::string("\n");
          file << tmp_str;
          file.close();
        }

        // start only thr_n jobs first, because the reconstruction can be done after the first feed
        uint32_t to_start = 0;

#if !WITH_MPI
        if (probes_for_next_prime > thr_n /** bunch_size*/) {
          to_start = thr_n /** bunch_size*/;

          {
#else
        if (mpi_first_send) {
          mpi_first_send = false;
          send_first_jobs(); // TODO send only start

          queue_probes(zi_order, thr_n); // * bunch_size
          started_probes[zi_order] += thr_n; // * bunch_size

          new_jobs = true;
          continue_communication = true;
          cond_val.notify_one();

          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          while (!proceed) {
            cond_val.wait(lock_probe_queue);
          }

          proceed = false;
        } else {
          if (probes_for_next_prime > buffer * total_thread_count) { // * bunch_size
            to_start = buffer * total_thread_count; // * bunch_size
#endif

            if (verbosity == CHATTY) {
              INFO_MSG("Starting " + std::to_string(to_start) + " jobs now, the remaining " + std::to_string(probes_for_next_prime - to_start) + " jobs will be started later");
            }

#if !WITH_MPI
            queue_probes(zi_order, to_start);
#else
            queue_probes(zi_order, to_start, true);
#endif
            started_probes.emplace(zi_order, to_start);
#if !WITH_MPI
          }
        } else {
          {
#else
          } else {
#endif
            to_start = probes_for_next_prime;

            if (verbosity == CHATTY) {
              INFO_MSG("Starting " + std::to_string(to_start) + " jobs");
            }

#if !WITH_MPI
            queue_probes(zi_order, to_start);
#else
            queue_probes(zi_order, to_start, true);
#endif
            started_probes.emplace(zi_order, to_start);
          }

#if !WITH_MPI
        }
#else
          cond_val.notify_one();

          {
            std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

            while (!proceed) {
              cond_val.wait(lock_probe_queue);
            }

            proceed = false;
          }

          for (uint32_t j = 0; j != to_start; ++j) {
            tp.run_task([this]() {
              get_job();
            });
          }
        }
#endif

        probes_for_next_prime = 0;
      }

      //if (iteration == 1000) exit(-1);

      std::vector<uint64_t> indices;
      std::vector<std::vector<FFInt>> probes;

      get_probe(indices, probes);

      if (verbosity > SILENT && !scan && std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - last_print_time).count() > 2.) {
        std::unique_lock<std::mutex> lock_print(print_control);
        last_print_time = std::chrono::high_resolution_clock::now();
        std::cerr << "\033[1;34mFireFly info:\033[0m Probe: " << iteration << "\r";
      }

      {
        std::unique_lock<std::mutex> lock_feed(feed_control);

        ++feed_jobs;
      }

      tp.run_priority_task([this, indices = std::move(indices), probes = std::move(probes)]() {
        feed_job(indices, probes);
      });

      {
        std::unique_lock<std::mutex> lock_status(status_control);

        if (items_done == items) {
          done = true;
          continue;
        } else if (items_done + items_new_prime == items) {
          new_prime = true;
          continue;
        }
      }

      {
        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        if (probes_queued == 0) {
          lock_probe_queue.unlock();

          {
            std::unique_lock<std::mutex> lock_feed(feed_control);

            while (feed_jobs > 0 || interpolate_jobs > 0) {
              condition_feed.wait(lock_feed);

              lock_feed.unlock();
              lock_probe_queue.lock();

              if (probes_queued != 0) {
                lock_probe_queue.unlock();

                break;
              }

              lock_probe_queue.unlock();
              lock_feed.lock();
            }

            lock_probe_queue.lock();

            if (probes_queued != 0) {
              continue;
            }

            lock_probe_queue.unlock();
          }

          // no jobs are running anymore, check if done or new_prime else throw error
          if (items_done == items) {
            done = true;
            continue;
          } else if (items_done + items_new_prime == items) {
            new_prime = true;
            continue;
          } else {
            throw std::runtime_error("Nothing left to feed: "
                                     + std::to_string(items)
                                     + " " + std::to_string(items_new_prime)
                                     + " " + std::to_string(items_done) + " | "
                                     + std::to_string(feed_jobs) + " "
                                     + std::to_string(interpolate_jobs) + " | "
                                     + std::to_string(iteration) + " "
                                     + std::to_string(fed_ones) + " | "
                                     + std::to_string(probes_queued) + " "
                                     + std::to_string(computed_probes.size()) + " "
                                     + std::to_string(requested_probes.size()) + "\n");
          }
        }
      }
    }

#if WITH_MPI
    {
      std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

      requested_probes = std::deque<std::pair<uint64_t, std::vector<FFInt>>>();
      new_jobs = false;

      while (!proceed) {
        cond_val.wait(lock_probe_queue);
      }

      proceed = false;
    }
#endif

    tp.kill_all();

    if (!factor_scan && !scan && save_states) {
      for (uint32_t item = 0; item != items; ++item) {
        // remove probe files if the interpolation is done
        std::remove(("ff_save/probes/" + std::to_string(item) + "_" + std::to_string(prime_it) + ".gz").c_str());
        ogzstream gzfile;
        std::string file_name = "ff_save/probes/" + std::to_string(item) + "_" + std::to_string(prime_it + 1) + ".gz";
        gzfile.open(file_name.c_str());
        gzfile.close();

        if (prime_it) {
          std::string file_name_old = "ff_save/states/" + std::to_string(item) + "_" + std::to_string(prime_it - 1) + ".gz";
          std::string file_name_new = "ff_save/states/" + std::to_string(item) + "_" + std::to_string(prime_it) + ".gz";
          std::rename(file_name_old.c_str(), file_name_new.c_str());
        }
      }
    }

#if WITH_MPI
    {
      std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

      continue_communication = true;
      cond_val.notify_one();

      while (!proceed) {
        cond_val.wait(lock_probe_queue);
      }

      proceed = false;
    }
#endif

    total_iterations += iteration;
  }

  // TODO optimize for bunch_size 1?
  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::queue_probes(const std::vector<uint32_t>& zi_order, const uint32_t to_start, const bool first) {
    bool ones = false;

    if ((prime_it == 0 || safe_mode == true) && (zi_order == std::vector<uint32_t> (n - 1, 1) || (factor_scan && zi_order == std::vector<uint32_t> (0, 1)))) {
      ones = true;
    }

    std::vector<FFInt> rand_zi;
    rand_zi.reserve(zi_order.size());

    if (factor_scan) {
      rand_zi = rand_zi_fac;
    } else if (!ones && (prime_it != 0 && safe_mode == false)) {
      rand_zi = tmp_rec.get_rand_zi_vec(zi_order, true);
    } else {
      rand_zi = tmp_rec.get_rand_zi_vec(zi_order, false);
    }//TODO elseif for factor

    for (uint32_t j = 0; j != to_start; ++j) {
      std::vector<FFInt> values(n);
      FFInt t = tmp_rec.get_rand_64();

      // check if t was already used for this zi_order
      {
        std::unique_lock<std::mutex> chosen_lock(chosen_mutex);

        auto it = chosen_t.find(zi_order);

        if (it != chosen_t.end()) {
          auto itt = it->second.find(t.n);

          if (itt != it->second.end()) {
            --j;

            //std::unique_lock<std::mutex> lock_print(print_control);

            //WARNING_MSG("Found a duplicate of t, choosing a new one");

            continue;
          } else {
            it->second.emplace(t.n);
          }
        } else {
          chosen_t.emplace(std::make_pair(zi_order, std::unordered_set<uint64_t>( {t.n})));
        }
      }

      if (!factor_scan) {
        if (change_var_order) {
          for (const auto & el : optimal_var_order) {
            if (el.first == 0) {
              values[el.second] = t + shift[0];
            } else {
              values[el.second] = rand_zi[el.first - 1] * t + shift[el.first];
            }
          }
        } else {
          values[0] = t + shift[0];

          for (uint32_t i = 1; i != n; ++i) {
            values[i] = rand_zi[i - 1] * t + shift[i];
          }
        }
      } else {
        for (uint32_t i = 0; i != n; ++i) {
          if (rand_zi_fac[i] == 1) {
            values[i] = t;
          } else {
            values[i] = FFInt(rand_zi_fac[i].n);
          }
        }
      }

      std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

      if (ones) {
        requested_probes.emplace_front(std::make_pair(ind, std::move(values)));
      } else {
        requested_probes.emplace_back(std::make_pair(ind, std::move(values)));
      }

      index_map.emplace(std::make_pair(ind, std::make_pair(t, zi_order)));
      //std::cout << "emplace " << ind << "\n";
      ++ind;

      ++probes_queued;
      //std::cout << "start " << probes_queued << "\n";
    }

#if WITH_MPI
    if (!first) {
#endif
    for (uint32_t j = 0; j != to_start; ++j) {
      tp.run_task([this]() {
        get_job();
      });
    }
#if WITH_MPI
    }
#endif

/*
    if (bunch_size == 1) {
#ifndef WITH_MPI
      std::vector<FFInt> values(n);
#endif
      FFInt t;

      for (uint32_t j = 0; j != to_start; ++j) {
#ifdef WITH_MPI
        std::vector<uint64_t> values(n + 1);
#endif
        t = tmp_rec.get_rand_64();

        // check if t was already used for this zi_order
        {
          std::unique_lock<std::mutex> chosen_lock(chosen_mutex);

          auto it = chosen_t.find(zi_order);

          if (it != chosen_t.end()) {
            auto itt = it->second.find(t.n);

            if (itt != it->second.end()) {
              --j;

              //std::unique_lock<std::mutex> lock_print(print_control);

              //WARNING_MSG("Found a duplicate of t, choosing a new one");

              continue;
            } else {
              it->second.emplace(t.n);
            }
          } else {
            chosen_t.emplace(std::make_pair(zi_order, std::unordered_set<uint64_t>( {t.n})));
          }
        }

#ifndef WITH_MPI
        values[0] = t + shift[0];

        for (uint32_t i = 1; i != n; ++i) {
          values[i] = rand_zi[i - 1] * t + shift[i];
        }

        // TODO: Why is this required already here?
        std::unique_lock<std::mutex> lock(future_control);

        if (ones) {
          probes.emplace_front(std::make_tuple(t, zi_order, probe_future()));
          auto it = probes.begin();

          probe_future future = tp.run_priority_packaged_task([this, values, it]() {
            auto time0 = std::chrono::high_resolution_clock::now();

            std::vector<FFInt> probe = bb(values);

            auto time1 = std::chrono::high_resolution_clock::now();

            std::unique_lock<std::mutex> lock(future_control);

            ++probes_finished;
            finished_probes_it.emplace(it);

            condition_future.notify_one();

            return std::make_pair(std::move(probe), std::chrono::duration<double>(time1 - time0).count());
          });

          std::get<2>(*it) = std::move(future);
        } else {
          probes.emplace_back(std::make_tuple(t, zi_order, probe_future()));
          auto it = --(probes.end());

          probe_future future = tp.run_packaged_task([this, values, it]() {
            auto time0 = std::chrono::high_resolution_clock::now();

            std::vector<FFInt> probe = bb(values);

            auto time1 = std::chrono::high_resolution_clock::now();

            std::unique_lock<std::mutex> lock(future_control);

            ++probes_finished;
            finished_probes_it.emplace(it);

            condition_future.notify_one();

            return std::make_pair(std::move(probe), std::chrono::duration<double>(time1 - time0).count());
          });

          std::get<2>(*it) = std::move(future);
        }
#else
        values[1] = (t + shift[0]).n;

        for (uint32_t i = 1; i != n; ++i) {
          values[i + 1] = (rand_zi[i - 1] * t + shift[i]).n;
        }

        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        values[0] = ind;

        index_map.emplace(std::make_pair(ind, std::make_pair(t, zi_order)));
        ++ind;

        value_queue.emplace(std::move(values));
#endif
      }
    } else {
      for (uint32_t j = 0; j != to_start / bunch_size; ++j) {
#ifndef WITH_MPI
        std::vector<FFInt> t_vec;
        t_vec.reserve(bunch_size);
        std::vector<std::vector<FFInt>> values_vec;
        values_vec.reserve(bunch_size);
#else
        // TODO this queues the jobs in a correct order so that zi_order is fixed for bunch_size entries
        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);
#endif

        for (uint32_t i = 0; i != bunch_size; ++i) {
#ifndef WITH_MPI
          std::vector<FFInt> values(n);
#else
          std::vector<uint64_t> values(n + 1);
#endif

          FFInt t = tmp_rec.get_rand_64();

          // check if t was already used for this zi_order
          {
            std::unique_lock<std::mutex> chosen_lock(chosen_mutex);

            auto it = chosen_t.find(zi_order);

            if (it != chosen_t.end()) {
              auto itt = it->second.find(t.n);

              if (itt != it->second.end()) {
                --i;

                //std::unique_lock<std::mutex> lock_print(print_control);

                //WARNING_MSG("Found a duplicate of t, choosing a new one");

                continue;
              } else {
                it->second.emplace(t.n);
              }
            } else {
              chosen_t.emplace(std::make_pair(zi_order, std::unordered_set<uint64_t>( {t.n})));
            }
          }

#ifndef WITH_MPI
          values[0] = t + shift[0];

          for (uint32_t i = 1; i != n; ++i) {
            values[i] = rand_zi[i - 1] * t + shift[i];
          }

          t_vec.emplace_back(t);
          values_vec.emplace_back(std::move(values));
#else
          values[1] = (t + shift[0]).n;

          for (uint32_t i = 1; i != n; ++i) {
            values[i + 1] = (rand_zi[i - 1] * t + shift[i]).n;
          }

          values[0] = ind;

          index_map.emplace(std::make_pair(ind, std::make_pair(t, zi_order)));
          ++ind;

          value_queue.emplace(std::move(values));
#endif
        }

#ifndef WITH_MPI
        // TODO: Why is this required already here?
        std::unique_lock<std::mutex> lock(future_control);

        if (ones) {
          probes_bunch.emplace_front(std::make_tuple(t_vec, zi_order, probe_future_bunch()));
          auto it = probes_bunch.begin();

          probe_future_bunch future = tp.run_priority_packaged_task([this, values_vec, it]() {
            auto time0 = std::chrono::high_resolution_clock::now();

            std::vector<std::vector<FFInt>> probe_vec = bb(values_vec);

            auto time1 = std::chrono::high_resolution_clock::now();

            std::unique_lock<std::mutex> lock(future_control);

            probes_finished += bunch_size;
            finished_probes_bunch_it.emplace(it);

            condition_future.notify_one();

            return std::make_pair(std::move(probe_vec), std::chrono::duration<double>(time1 - time0).count());
          });

          std::get<2>(*it) = std::move(future);
        } else {
          probes_bunch.emplace_back(std::make_tuple(t_vec, zi_order, probe_future_bunch()));
          auto it = --(probes_bunch.end());

          probe_future_bunch future = tp.run_packaged_task([this, values_vec, it]() {
            auto time0 = std::chrono::high_resolution_clock::now();

            std::vector<std::vector<FFInt>> probe_vec = bb(values_vec);

            auto time1 = std::chrono::high_resolution_clock::now();

            std::unique_lock<std::mutex> lock(future_control);

            probes_finished += bunch_size;
            finished_probes_bunch_it.emplace(it);

            condition_future.notify_one();

            return std::make_pair(std::move(probe_vec), std::chrono::duration<double>(time1 - time0).count());
          });

          std::get<2>(*it) = std::move(future);
        }
#endif
      }
    }

    std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

    probes_queued += to_start;
*/

#if WITH_MPI
    new_jobs = true;
#endif
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::get_probe(std::vector<uint64_t>& indices, std::vector<std::vector<FFInt>>& probes) {
    {
      std::unique_lock<std::mutex> lock_future(future_control);

      while (computed_probes.size() == 0) {
        condition_future.wait(lock_future);
      }

      //std::cout << computed_probes.size() << "\n";

      indices = std::move(computed_probes.front().first);
      probes = std::move(computed_probes.front().second);
      computed_probes.pop();

      //std::cout << "get " << indices.size() << "\n";
    }

    std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

    probes_queued -= static_cast<uint32_t>(indices.size());

    //std::cout << "got " << indices.size() << " " << probes_queued << "\n";
    //std::cout << "np " << items_new_prime << "\n";
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::feed_job(const std::vector<uint64_t>& indices, const std::vector<std::vector<FFInt>>& probes) {
    {
      std::unique_lock<std::mutex> lock(feed_control);

      interpolate_jobs += items;
    }

    std::vector<FFInt> t_vec;
    t_vec.reserve(indices.size());
    std::vector<std::vector<uint32_t>> zi_order_vec;
    zi_order_vec.reserve(indices.size());

    uint32_t count_ones = 0;

    {
      std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

      for (const auto& index : indices) {
        //std::cout << "access index " << index << "\n";
        auto tmp = std::move(index_map[index]);
        index_map.erase(index);
        t_vec.emplace_back(tmp.first);
        zi_order_vec.emplace_back(std::move(tmp.second));

        if ((factor_scan && static_cast<uint32_t>(zi_order_vec.back().size()) != 0) || (!factor_scan && static_cast<uint32_t>(zi_order_vec.back().size()) != n - 1)) {
          logger << "zi_order of probe has wrong length: " << std::to_string(zi_order_vec.back().size()) << "\n";
          ERROR_MSG("zi_order of probe has wrong length: " + std::to_string(zi_order_vec.back().size()));
          std::exit(EXIT_FAILURE);
        }

        if ((prime_it == 0 || safe_mode == true) && (zi_order_vec.back() == std::vector<uint32_t>(n - 1, 1) || (factor_scan && zi_order_vec.back() == std::vector<uint32_t>(0, 1)))) {
          ++count_ones;
        }
      }
    }

    if (count_ones != 0) {
      std::unique_lock<std::mutex> lock_status(job_control);

      fed_ones += count_ones;
    }

    uint32_t counter = 0;

    for (auto & rec : reconst) {
      std::unique_lock<std::mutex> lock_exists(*(std::get<1>(rec)));

      if (std::get<2>(rec) == RECONSTRUCTING) {
        lock_exists.unlock();

        std::pair<bool, uint32_t> done_prime = std::get<3>(rec)->get_done_and_prime();

        if (!done_prime.first) {
          if (done_prime.second == prime_it) {
            auto interpolate_and_write = std::get<3>(rec)->feed(t_vec, probes[std::get<0>(rec)], zi_order_vec, prime_it);

            if (interpolate_and_write.first) {
              ++counter;

              tp.run_priority_task([this, &rec]() {
                interpolate_job(rec);
              });
            }

            if (interpolate_and_write.second) {
              tp.run_priority_task([this, &rec]() {
                std::get<3>(rec)->write_food_to_file();
              });
            }
          }
        }
      }
    }

    {
      std::unique_lock<std::mutex> lock_status(status_control);

      if (!scan && (one_done || one_new_prime)) {
        one_done = false;
        one_new_prime = false;

        std::unique_lock<std::mutex> lock_future(future_control);
        std::unique_lock<std::mutex> lock_print(print_control);
        logger << "Probe: " + std::to_string(iteration)
        << " | Done: " << std::to_string(items_done) << " / " << std::to_string(items)
        << " | " << "Needs new prime field: " << std::to_string(items_new_prime)
        << " / " << std::to_string(items - items_done) << "\n";
        logger.close();
        logger.open("firefly.log", std::ios_base::app);

        if (verbosity > SILENT) {
          INFO_MSG("Probe: " + std::to_string(iteration) +
                   " | Done: " + std::to_string(items_done) + " / " + std::to_string(items) +
                   " | " + "Needs new prime field: " + std::to_string(items_new_prime) + " / " + std::to_string(items - items_done));
        }
      }
    }

    std::unique_lock<std::mutex> lock(feed_control);

    interpolate_jobs -= (items - counter);
    --feed_jobs;
    condition_feed.notify_one();
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::interpolate_job(RatReconst_tuple& rec) {
    std::unique_lock<std::mutex> lock_exists(*(std::get<1>(rec)));

    if (std::get<2>(rec) == RECONSTRUCTING) {
      lock_exists.unlock();

      std::tuple<bool, bool, uint32_t> interpolated_done_prime = std::get<3>(rec)->interpolate();

      if (std::get<0>(interpolated_done_prime)) { // interpolated
        //lock_exists.lock();
        // start new jobs if required
        if (!std::get<1>(interpolated_done_prime)) { // not done
          if (std::get<2>(interpolated_done_prime) > prime_it) { // compare prime counters
            {
              std::unique_lock<std::mutex> lock_status(status_control);

              one_new_prime = true;
              ++items_new_prime;
            }

            std::unique_lock<std::mutex> lock(job_control);

            if (std::get<3>(rec)->get_num_eqn() > probes_for_next_prime) {
              probes_for_next_prime = std::get<3>(rec)->get_num_eqn();
            }

            //lock_exists.unlock();
          } else if (!safe_mode && prime_it != 0) {
            std::vector<std::pair<uint32_t, uint32_t>> all_required_probes = std::get<3>(rec)->get_needed_feed_vec();

            if (!all_required_probes.empty()) {
              uint32_t counter = 1;

              for (const auto & some_probes : all_required_probes) {
                if (some_probes.first == 0) {
                  ++counter;
                  continue;
                }

                uint32_t required_probes = some_probes.second;

                for (uint32_t i = 0; i != some_probes.first; ++i) {
                  std::vector<uint32_t> zi_order;

                  if (!factor_scan) {
                    zi_order = std::vector<uint32_t>(n - 1, counter);
                  } else {
                    zi_order = std::vector<uint32_t>(0, counter);
                  }

                  ++counter;

                  std::unique_lock<std::mutex> lock(job_control);

                  auto it = started_probes.find(zi_order);

                  if (it != started_probes.end()) {
                    if (required_probes > started_probes[zi_order]) {
                      uint32_t to_start = required_probes - started_probes[zi_order];

                      started_probes[zi_order] = required_probes;

                      lock.unlock();

                      if (verbosity == CHATTY) {
                        std::string msg = "Starting zi_order (";

                        for (const auto & ele : zi_order) {
                          msg += std::to_string(ele) + ", ";
                        }

                        msg = msg.substr(0, msg.length() - 2);
                        msg += ") " + std::to_string(to_start) + " time(s)";

                        std::unique_lock<std::mutex> lock_print(print_control);

                        INFO_MSG(msg);
                      }

                      queue_probes(zi_order, to_start);
                    }
                  } else {
                    started_probes.emplace(zi_order, required_probes);

                    lock.unlock();

                    if (verbosity == CHATTY) {
                      std::string msg = "Starting zi_order (";

                      for (const auto & ele : zi_order) {
                        msg += std::to_string(ele) + ", ";
                      }

                      msg = msg.substr(0, msg.length() - 2);
                      msg += ") " + std::to_string(required_probes) + " time(s)";

                      std::unique_lock<std::mutex> lock_print(print_control);

                      INFO_MSG(msg);
                    }

                    queue_probes(zi_order, required_probes);
                  }
                }
              }
            }
          } else {
            std::vector<uint32_t> zi_order = std::get<3>(rec)->get_zi_order();

            if ((prime_it == 0 || safe_mode == true) && (zi_order == std::vector<uint32_t>(n - 1, 1) || (factor_scan && zi_order == std::vector<uint32_t>(0, 1)))) {
              //lock_exists.unlock();
              std::unique_lock<std::mutex> lock(job_control);

#if !WITH_MPI
              if (started_probes[zi_order] - thr_n /** bunch_size*/ <= fed_ones - 1) {
                uint32_t to_start = fed_ones - started_probes[zi_order] + thr_n /** bunch_size*/;
#else
              if (started_probes[zi_order] - total_thread_count <= fed_ones - 1) { // TODO static_cast<uint32_t>(buffer) * * bunch_size
                uint32_t to_start = fed_ones - started_probes[zi_order] + static_cast<uint32_t>(buffer) * total_thread_count; //TODO * bunch_size
#endif
                started_probes[zi_order] += to_start;

                lock.unlock();

                if (verbosity == CHATTY) {
                  std::unique_lock<std::mutex> lock_print(print_control);

                  INFO_MSG("Starting ones: " + std::to_string(to_start));
                }

                queue_probes(zi_order, to_start);
              }
            } else {
              uint32_t required_probes = std::get<3>(rec)->get_num_eqn();

              //lock_exists.unlock();
              std::unique_lock<std::mutex> lock(job_control);

              auto it = started_probes.find(zi_order);

              if (it != started_probes.end()) {
                if (required_probes > started_probes[zi_order]) {
                  uint32_t to_start = required_probes - started_probes[zi_order];

                  started_probes[zi_order] = required_probes;

                  lock.unlock();

                  if (verbosity == CHATTY) {
                    std::string msg = "Starting zi_order (";

                    for (const auto & ele : zi_order) {
                      msg += std::to_string(ele) + ", ";
                    }

                    msg = msg.substr(0, msg.length() - 2);
                    msg += ") " + std::to_string(to_start) + " time(s)";

                    std::unique_lock<std::mutex> lock_print(print_control);

                    INFO_MSG(msg);
                  }

                  queue_probes(zi_order, to_start);
                }
              } else {
                started_probes.emplace(zi_order, required_probes);

                lock.unlock();

                if (verbosity == CHATTY) {
                  std::string msg = "Starting zi_order (";

                  for (const auto & ele : zi_order) {
                    msg += std::to_string(ele) + ", ";
                  }

                  msg = msg.substr(0, msg.length() - 2);
                  msg += ") " + std::to_string(required_probes) + " time(s)";

                  std::unique_lock<std::mutex> lock_print(print_control);

                  INFO_MSG(msg);
                }

                queue_probes(zi_order, required_probes);
              }
            }
          }
        } else {
          lock_exists.lock();

          // to be sure that no other thread does the same
          // TODO: Why is this necessary?
          if (std::get<2>(rec) == RECONSTRUCTING) {
            std::get<2>(rec) = DONE;

            lock_exists.unlock();
            std::unique_lock<std::mutex> lock_status(status_control);

            ++items_done;
            one_done = true;
          } else {
            lock_exists.unlock();
          }
        }
      }
    } else {
      lock_exists.unlock();
    }

    std::unique_lock<std::mutex> lock(feed_control);

    --interpolate_jobs;
    condition_feed.notify_one();
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::clean_reconst() {
    std::unique_lock<std::mutex> lock_clean(clean);

    auto it = reconst.begin();

    while (it != reconst.end()) {
      if (std::get<2>(*it) == DELETED) {
        // delete mutex
        delete std::get<1>(*it);

        // remove from list
        it = reconst.erase(it);
      } else {
        ++it;
      }
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::abort() {
    logger << "Run aborted\n";
    INFO_MSG("Run aborted");
    done = true;
    aborted = true;
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::resume() {
    logger << "Run resumed\n";
    INFO_MSG("Run resumed");
    if(aborted) {
      done = false;
      aborted = false;
      resumed = true;
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::get_job() {
    std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

    if (!requested_probes.empty()) {
      switch(compute_bunch_size(static_cast<uint32_t>(requested_probes.size()), thr_n, bunch_size)) {
        case 1:
          start_new_job<1>(lock_probe_queue);
          break;
        case 2:
          start_new_job<2>(lock_probe_queue);
          break;
        case 4:
          start_new_job<4>(lock_probe_queue);
          break;
        case 8:
          start_new_job<8>(lock_probe_queue);
          break;
        case 16:
          start_new_job<16>(lock_probe_queue);
          break;
        case 32:
          start_new_job<32>(lock_probe_queue);
          break;
/*        case 64:
          start_new_job<64>(lock_probe_queue);
          break;
        case 128:
          start_new_job<128>(lock_probe_queue);
          break;
        case 256:
          start_new_job<256>(lock_probe_queue);
          break;*/
      }

#if WITH_MPI
      if (requested_probes.empty()) {
        new_jobs = false;
      }
#endif
    }
  }

  template<typename BlackBoxTemp>
  template<uint32_t N>
  void Reconstructor<BlackBoxTemp>::start_new_job(std::unique_lock<std::mutex>& lock_probe_queue) {
    std::vector<uint64_t> indices;
    indices.reserve(N);

    //std::cout << "start with " << N << " of " << requested_probes.size() << "\n";

    if (N != 1) {
      std::vector<FFIntVec<N>> values_vec(n); // TODO do not fill

      for (uint32_t i = 0; i != N; ++i) {
        indices.emplace_back(requested_probes.front().first);

        for (uint32_t j = 0; j != n; ++j) {
          values_vec[j][i] = requested_probes.front().second[j];
        }

        requested_probes.pop_front();
      }

      lock_probe_queue.unlock();

      auto time0 = std::chrono::high_resolution_clock::now();

      std::vector<FFIntVec<N>> probe = bb.eval(values_vec);

      auto time1 = std::chrono::high_resolution_clock::now();

      auto time = std::chrono::duration<double>(time1 - time0).count();

      // TODO
      std::vector<std::vector<FFInt>> tmp;//(probe.size(), std::vector<FFInt>(N));
      tmp.reserve(probe.size());

      for (size_t i = 0; i != probe.size(); ++i) {
        //std::move(probe[i].begin(), probe[i].end(), tmp[i].begin());
        // Remove factor from bb result
        if (!factor_scan && parsed_factors.find(i) != parsed_factors.end()) {
          auto res = parsed_factors[i].evaluate_pre(values_vec);

          for (size_t j = 0; j != N; ++j) {
            probe[i][j] /= res[0][j];
          }
        }

        tmp.emplace_back(std::vector<FFInt>(std::make_move_iterator(probe[i].begin()), std::make_move_iterator(probe[i].end())));
      }

      std::unique_lock<std::mutex> lock(future_control);

      iteration += N;
#if !WITH_MPI
      int tmp_iterations = total_iterations + iteration;
      average_black_box_time = (average_black_box_time * (tmp_iterations - N) + time) / tmp_iterations;
#else
      iterations_on_this_node += N;
      int tmp_iterations = total_iterations + iterations_on_this_node;
      average_black_box_time = (average_black_box_time * (tmp_iterations - N) + time) / tmp_iterations;
#endif

      computed_probes.emplace(std::make_pair(std::move(indices), std::move(tmp)));

      condition_future.notify_one();
    } else {
      indices.emplace_back(requested_probes.front().first);

      std::vector<FFInt> values_vec(std::make_move_iterator(requested_probes.front().second.begin()), std::make_move_iterator(requested_probes.front().second.end()));

      //for (uint32_t j = 0; j != n; ++j) {
      //  values_vec.emplace_back(requested_probes.front().second[j]);
      //}

      requested_probes.pop_front();

      lock_probe_queue.unlock();

      auto time0 = std::chrono::high_resolution_clock::now();

      std::vector<FFInt> probe = bb.eval(values_vec);

      auto time1 = std::chrono::high_resolution_clock::now();

      auto time = std::chrono::duration<double>(time1 - time0).count();

      std::vector<std::vector<FFInt>> tmp;
      tmp.reserve(probe.size());

      for (size_t i = 0; i != probe.size(); ++i) {
        if (!factor_scan && parsed_factors.find(i) != parsed_factors.end()) {
          auto res = parsed_factors[i].evaluate_pre(values_vec);

          for (size_t j = 0; j != N; ++j) {
            probe[i] /= res[0];
          }
        }

        tmp.emplace_back(std::vector<FFInt> (1, probe[i]));
      }

      std::unique_lock<std::mutex> lock(future_control);

      iteration += N;
#if !WITH_MPI
      int tmp_iterations = total_iterations + iteration;
      average_black_box_time = (average_black_box_time * (tmp_iterations - N) + time) / tmp_iterations;
#else
      iterations_on_this_node += N;
      int tmp_iterations = total_iterations + iterations_on_this_node;
      average_black_box_time = (average_black_box_time * (tmp_iterations - N) + time) / tmp_iterations;
#endif

      computed_probes.emplace(std::make_pair(std::move(indices), std::move(tmp)));

      condition_future.notify_one();
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::reset_new_prime() {
    iteration = 0;
  #if WITH_MPI
    iterations_on_this_node = 0;
  #endif

    fed_ones = 0;
    probes_queued = 0;
    started_probes.clear();
    index_map.clear();
    ind = 0;
    feed_jobs = 0;
    interpolate_jobs = 0;
    new_prime = false;
    items_new_prime = 0;
    one_done = false;
    one_new_prime = false;

    // Only reset chosen_t when not resuming from a saved state
    if (!set_anchor_points) {
      chosen_t.clear();
    }

    // TODO mutex required here?
    requested_probes = std::deque<std::pair<uint64_t, std::vector<FFInt>>>();
    computed_probes = std::queue<std::pair<std::vector<uint64_t>, std::vector<std::vector<FFInt>>>>();
  }

#if WITH_MPI
  template<typename BlackBoxTemp>
  inline void Reconstructor<BlackBoxTemp>::mpi_setup() {
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Bcast(&prime_it, 1, MPI_UINT32_T, master, MPI_COMM_WORLD);
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::send_first_jobs() {
    mpi_setup();

    std::vector<uint32_t> zi_order = std::vector<uint32_t>(n - 1, 1);

    std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

    if (started_probes.find(zi_order) != started_probes.end()) {
      started_probes.emplace(std::make_pair(zi_order, 0));
    }

    std::vector<FFInt> rand_zi = tmp_rec.get_rand_zi_vec(zi_order, true);

    std::unique_lock<std::mutex> chosen_lock(chosen_mutex); // TODO is this really required?

    for (int i = 1; i != world_size; ++i) {
      uint64_t to_start;
      MPI_Recv(&to_start, 1, MPI_UINT64_T, i, RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      //std::cout << "starting " << to_start << " jobs on worker " << i << "\n";

      nodes.emplace(std::make_pair(i, to_start));
      total_thread_count += static_cast<uint32_t>(to_start / buffer);

      probes_queued += to_start;
      started_probes[zi_order] += to_start;

      //std::cout << "first send " << probes_queued << "\n";

      std::vector<uint64_t> values;
      values.reserve(static_cast<uint32_t>(to_start) * (n + 1));

      for (uint32_t ii = 0; ii != static_cast<uint32_t>(to_start); ++ii) {
        //values[ii * (n + 1)] = ind;
        values.emplace_back(ind);

        FFInt t = tmp_rec.get_rand_64();

        // check if t was already used for this zi_order
        auto it = chosen_t.find(zi_order);

        if (it != chosen_t.end()) {
          auto itt = it->second.find(t.n);

          if (itt != it->second.end()) {
            --ii;

            //std::unique_lock<std::mutex> lock_print(print_control);

            //WARNING_MSG("Found a duplicate of t, choosing a new one");

            continue;
          } else {
            it->second.emplace(t.n);
          }
        } else {
          chosen_t.emplace(std::make_pair(zi_order, std::unordered_set<uint64_t>( {t.n})));
        }

        //values[ii * (n + 1) + 1] = (t + shift[0]).n;
        values.emplace_back((t + shift[0]).n);

        for (uint32_t j = 1; j != n; ++j) {
          //values[ii * (n + 1) + j + 1] = (t * rand_zi[j - 1] + shift[j]).n;
          values.emplace_back((t * rand_zi[j - 1] + shift[j]).n);
        }

        index_map.emplace(std::make_pair(ind, std::make_pair(t, zi_order)));
        //std::cout << "mpi emplace " << ind << "\n";
        ++ind;
      }

      MPI_Send(&values[0], static_cast<int>(static_cast<uint32_t>(to_start) * (n + 1)), MPI_UINT64_T, i, VALUES, MPI_COMM_WORLD);
    }
  }

  template<typename BlackBoxTemp>
  void Reconstructor<BlackBoxTemp>::mpi_communicate() {
    while (true) {
      int flag_ext = 0;
      MPI_Status status;

      bool restart_empty_nodes = false;

      while (!flag_ext) {
        MPI_Iprobe(MPI_ANY_SOURCE, RESULT, MPI_COMM_WORLD, &flag_ext, &status);

        //std::cout << "empty size: " << empty_nodes.size() << " " << new_jobs.load() << "\n";

        if (new_prime || done) {
          break;
        } else if (!empty_nodes.empty() && new_jobs) {
          restart_empty_nodes = true;
          break;
        }
      }

      if (done && !scan) {
        //std::vector<MPI_Request> requests;
        //requests.reserve(world_size - 1);
        MPI_Request* requests = new MPI_Request[world_size - 1];
        uint64_t tmp;

        for (int i = 1; i != world_size; ++i) {
          MPI_Isend(&tmp, 1, MPI_UINT64_T, i, END, MPI_COMM_WORLD, &requests[i - 1]);
        }

        MPI_Waitall(world_size - 1, requests, MPI_STATUSES_IGNORE);

        delete[] requests;

        MPI_Barrier(MPI_COMM_WORLD);

        {
          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          proceed = true;

          cond_val.notify_one();
        }

        int flag = 1;
        MPI_Status status_rec;

        while (flag) {
          MPI_Iprobe(MPI_ANY_SOURCE, RESULT, MPI_COMM_WORLD, &flag, &status_rec);

          if (flag) {
            int amount;
            MPI_Get_count(&status_rec, MPI_UINT64_T, &amount);

            std::vector<uint64_t> tmp;
            tmp.reserve(amount);

            MPI_Recv(&tmp[0], amount, MPI_UINT64_T, status_rec.MPI_SOURCE, status_rec.MPI_TAG, MPI_COMM_WORLD, &status_rec);
          }
        }

        std::vector<double> timings;
        std::vector<double> weights;
        timings.reserve(world_size);
        weights.reserve(world_size);
        double total_weight = 0.;

        for (int i = 1; i != world_size; ++i) {
          MPI_Status status_time;
          double timing[2];
          MPI_Recv(&timing, 2, MPI_DOUBLE, i, TIMING, MPI_COMM_WORLD, &status_time);
          timings.emplace_back(timing[0]);
          weights.emplace_back(timing[1]);
          total_weight += timing[1];
        }

        MPI_Barrier(MPI_COMM_WORLD);

        {
          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          while (!continue_communication) {
            cond_val.wait(lock_probe_queue);
          }

          continue_communication = false;
        }

        iteration = total_weight + iterations_on_this_node;

        timings.emplace_back(average_black_box_time);
        weights.emplace_back(static_cast<double>(total_iterations + iterations_on_this_node));
        total_weight += static_cast<double>(total_iterations + iterations_on_this_node);

        average_black_box_time = 0.;

        for (int i = 0; i != world_size; ++i) {
          average_black_box_time += weights[i] / total_weight * timings[i];
        }

        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        proceed = true;

        cond_val.notify_one();

        break;
      } else if (new_prime || (done && scan)) {
        MPI_Request* requests = new MPI_Request[world_size - 1];
        //std::vector<MPI_Request> requests;
        //requests.reserve(world_size - 1);

        //std::cout << "com np\n";

        uint64_t prime_tmp = static_cast<uint64_t>(prime_it);

        if (!scan) {
          ++prime_tmp;
        }

        for (int i = 1; i != world_size; ++i) {
          MPI_Isend(&prime_tmp, 1, MPI_UINT64_T, i, NEW_PRIME, MPI_COMM_WORLD, &requests[i - 1]);
        }

        MPI_Waitall(world_size - 1, requests, MPI_STATUSES_IGNORE);

        delete[] requests;

        MPI_Barrier(MPI_COMM_WORLD);

        {
          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          proceed = true;

          cond_val.notify_one();
        }

        int flag = 1;
        MPI_Status status_rec;

        while (flag) {
          MPI_Iprobe(MPI_ANY_SOURCE, RESULT, MPI_COMM_WORLD, &flag, &status_rec);

          if (flag) {
            int amount;
            MPI_Get_count(&status_rec, MPI_UINT64_T, &amount);

            //std::cout << "garbage recieved: " << (static_cast<uint32_t>(amount) - 1) / (items + 1) << "\n";

            std::vector<uint64_t> tmp;
            tmp.reserve(amount);

            MPI_Recv(&tmp[0], amount, MPI_UINT64_T, status_rec.MPI_SOURCE, status_rec.MPI_TAG, MPI_COMM_WORLD, &status_rec);
          }
        }

        std::vector<double> timings;
        std::vector<double> weights;
        timings.reserve(world_size);
        weights.reserve(world_size);
        double total_weight = 0.;

        for (int i = 1; i != world_size; ++i) {
          MPI_Status status_time;
          double timing[2];
          MPI_Recv(&timing, 2, MPI_DOUBLE, i, TIMING, MPI_COMM_WORLD, &status_time);
          timings.emplace_back(timing[0]);
          weights.emplace_back(timing[1]);
          total_weight += timing[1];
        }

        MPI_Barrier(MPI_COMM_WORLD);

        {
          std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

          while (!continue_communication) {
            cond_val.wait(lock_probe_queue);
          }

          continue_communication = false;
        }

        iteration = total_weight + iterations_on_this_node;

        timings.emplace_back(average_black_box_time);
        weights.emplace_back(static_cast<double>(total_iterations + iterations_on_this_node));
        total_weight += static_cast<double>(total_iterations + iterations_on_this_node);

        average_black_box_time = 0.;

        for (int i = 0; i != world_size; ++i) {
          average_black_box_time += weights[i] / total_weight * timings[i];
        }

        //std::cout << "com rec\n";

        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        proceed = true;

        cond_val.notify_one();

        while (!new_jobs) {
          cond_val.wait(lock_probe_queue);
        }

        proceed = true;

        cond_val.notify_one();

        empty_nodes = std::queue<std::pair<int, uint64_t>>();

        //std::cout << "com through\n";

        for (int k = 1; k != world_size; ++k) {
          uint64_t free_slots;
          MPI_Recv(&free_slots, 1, MPI_UINT64_T, k, RESULT, MPI_COMM_WORLD, &status);

          if (requested_probes.size() != 0) {
            uint32_t size = compute_job_number(static_cast<uint32_t>(requested_probes.size()), static_cast<uint32_t>(free_slots), thr_n, bunch_size);

            std::vector<uint64_t> values;
            values.reserve(size * (n + 1));

            for (uint64_t i = 0; i != size; ++i) {
              //values[i * (n + 1)] = requested_probes.front().first;
              values.emplace_back(requested_probes.front().first);

              for (uint64_t j = 0; j != n; ++j) {
                //values[i * (n + 1) + 1 + j] = requested_probes.front().second[j].n;
                values.emplace_back(requested_probes.front().second[j].n);
              }

              requested_probes.pop_front();
            }

            if (requested_probes.empty()) {
              new_jobs = false;
            }

            //std::cout << "sending " << size << " jobs\n";

            MPI_Send(&values[0], static_cast<int>(size * (n + 1)), MPI_UINT64_T, status.MPI_SOURCE, VALUES, MPI_COMM_WORLD);
          } else if (free_slots == nodes[status.MPI_SOURCE]) {
            //std::cout << "com np empty nodes " << status.MPI_SOURCE << " " << free_slots << "\n";
            empty_nodes.emplace(std::make_pair(status.MPI_SOURCE, free_slots));
          }
        }
      } else if (restart_empty_nodes) {
        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        while (!empty_nodes.empty() && !requested_probes.empty()) {
          auto node = empty_nodes.front();
          empty_nodes.pop();

          uint32_t size = compute_job_number(static_cast<uint32_t>(requested_probes.size()), static_cast<uint32_t>(node.second), thr_n, bunch_size);

          std::vector<uint64_t> values;
          values.reserve(size * (n + 1));

          for (uint64_t i = 0; i != size; ++i) {
            //values[i * (n + 1)] = requested_probes.front().first;
            values.emplace_back(requested_probes.front().first);

            for (uint64_t j = 0; j != n; ++j) {
              //values[i * (n + 1) + 1 + j] = requested_probes.front().second[j].n;
              values.emplace_back(requested_probes.front().second[j].n);
            }

            requested_probes.pop_front();
          }

          //std::cout << "restart: sending " << size << " jobs\n";

          MPI_Send(&values[0], static_cast<int>(size * (n + 1)), MPI_UINT64_T, status.MPI_SOURCE, VALUES, MPI_COMM_WORLD);
        }

        if (requested_probes.empty()) {
          new_jobs = false;
        }
      } else {
        int amount;
        MPI_Get_count(&status, MPI_UINT64_T, &amount);

        if ((static_cast<uint32_t>(amount) - 1) % (items + 1) != 0) {
          logger << "Corrupted results recieved: " + std::to_string(amount - 1) << "\n";
          ERROR_MSG("Corrupted results recieved: " + std::to_string(amount - 1));
          std::exit(EXIT_FAILURE);
        }

        uint32_t new_results = (static_cast<uint32_t>(amount) - 1) / (items + 1);

        //std::cout << "comm recieving " << new_results << " results \n";

        // TODO optimize the format
        std::vector<uint64_t> results_list;
        results_list.reserve(amount);
        MPI_Recv(&results_list[0], amount, MPI_UINT64_T, status.MPI_SOURCE, RESULT, MPI_COMM_WORLD, &status);

        for (uint32_t i = 0; i != new_results; ++i) {
          uint64_t index = results_list[i * (items + 1)];
          std::vector<std::vector<FFInt>> results(items, std::vector<FFInt>(1)); // TODO do not initialize here
          //results.reserve(items);

          for (uint32_t j = 1; j != items + 1; ++j) {
            results[j - 1][0] = results_list[i * (items + 1) + j];
          }

          std::unique_lock<std::mutex> lock_res(future_control);

          computed_probes.emplace(std::make_pair(std::vector<uint64_t>(1, index), std::move(results)));
        }

        condition_future.notify_one(); // TODO mutex?

        uint64_t free_slots = results_list[amount - 1];

        std::unique_lock<std::mutex> lock_probe_queue(mutex_probe_queue);

        if (requested_probes.size() != 0) {
          uint32_t size = compute_job_number(static_cast<uint32_t>(requested_probes.size()), static_cast<uint32_t>(free_slots), thr_n, bunch_size);

          std::vector<uint64_t> values;
          values.reserve(size * (n + 1));

          for (uint64_t i = 0; i != size; ++i) {
            //values[i * (n + 1)] = requested_probes.front().first;
            values.emplace_back(requested_probes.front().first);

            for (uint64_t j = 0; j != n; ++j) {
              //values[i * (n + 1) + 1 + j] = requested_probes.front().second[j].n;
              values.emplace_back(requested_probes.front().second[j].n);
            }

            requested_probes.pop_front();
          }

          if (requested_probes.empty()) {
            new_jobs = false;
          }

          //std::cout << "sending " << size << " jobs\n";

          MPI_Send(&values[0], static_cast<int>(size * (n + 1)), MPI_UINT64_T, status.MPI_SOURCE, VALUES, MPI_COMM_WORLD);
        } else if (free_slots == nodes[status.MPI_SOURCE]) {
          //std::cout << "empty node " << status.MPI_SOURCE << " " << free_slots << "\n";
          empty_nodes.emplace(std::make_pair(status.MPI_SOURCE, free_slots));
        } else {
          //std::cout << "sending NULL\n";
          MPI_Send(NULL, 0, MPI_UINT64_T, status.MPI_SOURCE, VALUES, MPI_COMM_WORLD);
        }
      }
    }
  }
#endif
}
