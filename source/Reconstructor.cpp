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

#include "gzstream.hpp"
#include "ParserUtils.hpp"
#include "Reconstructor.hpp"
#include "ReconstHelper.hpp"
#include "tinydir.h"
#include "utils.hpp"
#include "version.hpp"

#ifdef WITH_MPI
#include "MPIWorker.hpp"
#endif

#include <sys/stat.h>

namespace firefly {
  Reconstructor::Reconstructor(uint32_t n_, uint32_t thr_n_, BlackBoxBase& bb_, uint32_t verbosity_): n(n_), thr_n(thr_n_), bb(bb_), verbosity(verbosity_), tp(thr_n_) {
    FFInt::set_new_prime(primes()[prime_it]);
    uint64_t seed = static_cast<uint64_t>(std::time(0));
    BaseReconst().set_seed(seed);
    tmp_rec = RatReconst(n);

    if (verbosity > SILENT) {
      std::cout << "\nFire\033[1;32mFly\033[0m " << FireFly_VERSION_MAJOR << "." << FireFly_VERSION_MINOR << "." << FireFly_VERSION_RELEASE << "\n\n";
      INFO_MSG("Launching " << thr_n_ << " thread(s) with bunch size 1");
      INFO_MSG("Using seed " + std::to_string(seed) + " for random numbers");
    }
  }

  Reconstructor::Reconstructor(uint32_t n_, uint32_t thr_n_, uint32_t bunch_size_, BlackBoxBase& bb_, uint32_t verbosity_): n(n_), thr_n(thr_n_), bunch_size(bunch_size_), bb(bb_), verbosity(verbosity_), tp(thr_n_) {
    FFInt::set_new_prime(primes()[prime_it]);
    uint64_t seed = static_cast<uint64_t>(std::time(0));
    BaseReconst().set_seed(seed);
    tmp_rec = RatReconst(n);

    if (verbosity > SILENT) {
      std::cout << "\nFire\033[1;32mFly\033[0m " << FireFly_VERSION_MAJOR << "." << FireFly_VERSION_MINOR << "." << FireFly_VERSION_RELEASE << "\n\n";
      INFO_MSG("Launching " << thr_n_ << " thread(s) with bunch size " + std::to_string(bunch_size_));
      INFO_MSG("Using seed " + std::to_string(seed) + " for random numbers");
    }
  }

  Reconstructor::~Reconstructor() {
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

  void Reconstructor::enable_scan() {
    if (n == 1) {
      WARNING_MSG("Scan disabled for a univariate rational function.");
    } else {
      scan = true;
    }
  }

  void Reconstructor::set_tags() {
    save_states = true;
  }

  void Reconstructor::set_tags(const std::vector<std::string>& tags_) {
    save_states = true;
    tags = tags_;
  }

  void Reconstructor::resume_from_saved_state() {
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
      INFO_MSG("Starting new reconstruction and saving states");
      return;
    }
  }

  void Reconstructor::resume_from_saved_state(const std::vector<std::string>& file_paths_) {
    if (verbosity > SILENT) {
      INFO_MSG("Loading saved states");
    }

    set_anchor_points = true;

    std::ifstream v_file;
    igzstream validation_file;
    validation_file.open("ff_save/validation.gz");
    v_file.open("ff_save/validation.gz");
    std::string line;

    if (v_file.is_open()) {
      std::getline(validation_file, line);
      std::vector<FFInt> values = parse_vector_FFInt(line);

      std::vector<FFInt> result = bb(values);
      size_t counter = 0;

      while (std::getline(validation_file, line)) {
        if (std::stoul(line) != result[counter]) {
          ERROR_MSG("Validation failed: Entry " + std::to_string(counter) + " does not match the black-box result!");
          std::exit(EXIT_FAILURE);
        }

        ++counter;
      }

      if (counter != result.size()) {
        ERROR_MSG("Validation failed: Number of entries does not match the black box!");
        std::exit(EXIT_FAILURE);
      }
    } else {
      ERROR_MSG("Validation file not found!");
      std::exit(EXIT_FAILURE);
    }

    v_file.close();
    validation_file.close();

    save_states = true;
    resume_from_state = true;
    file_paths = file_paths_;
    items = file_paths.size();
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
      std::exit(EXIT_FAILURE);
    }

    for (uint32_t i = 0; i != items; ++i) {
      RatReconst* rec = new RatReconst(n);
      std::pair<bool, uint32_t> shift_prime = rec->start_from_saved_file(file_paths[i]);

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

      if (shift_prime.first && shift_prime.second > min_prime_keep_shift) {
        min_prime_keep_shift = shift_prime.second;
      }

      rec->set_tag(std::to_string(i));

      if (rec->is_done()) {
        ++items_done;
        std::mutex* mut = new std::mutex;

        reconst.emplace_back(std::make_tuple(i, mut, DONE, rec));
      } else {
        if (rec->is_new_prime()) {
          probes_for_next_prime = std::max(probes_for_next_prime, rec->get_num_eqn());
          ++items_new_prime;
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
        tmp_rec.set_anchor_points(parse_vector_FFInt(line, n));
      } else {
        ERROR_MSG("Anchor point file not found!");
        std::exit(EXIT_FAILURE);
      }

      anchor_point_file.close();

      std::ifstream shift_file;
      shift_file.open("ff_save/shift");

      if (shift_file.is_open()) {
        std::getline(shift_file, line);
        tmp_rec.set_shift(parse_vector_FFInt(line, n));
        shift = tmp_rec.get_zi_shift_vec();
      } else {
        ERROR_MSG("Shift file not found!");
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
          std::exit(EXIT_FAILURE);
        }

        file.close();
      } else {
        scan = false;
      }
    }

    if (verbosity > SILENT) {
      INFO_MSG("All files loaded | Done: " + std::to_string(items_done) + " / " + std::to_string(items) +
               " | " + "Needs new prime field: " + std::to_string(items_new_prime) + " / " + std::to_string(items - items_done));
    }
  }

  void Reconstructor::set_safe_interpolation() {
    safe_mode = true;
  }

  void Reconstructor::reconstruct() {
    start = std::chrono::high_resolution_clock::now();

    done = false;

#ifdef WITH_MPI
    mpi_setup();
    ThreadPool tp_comm(1);
    tp_comm.run_priority_task([this]() {
      {
        std::unique_lock<std::mutex> lock(mut_val);

        while (!new_jobs) {
          cond_val.wait(lock);
        }

        if (value_queue.empty()) {
          new_jobs = false;
        }
      }

      mpi_communicate();
    });
#endif

    if (!resume_from_state) {
      if (verbosity > SILENT) {
        std::cout << "\n";
        INFO_MSG("Promote to new prime field: F(" + std::to_string(primes()[prime_it]) + ")");
      }

      if (safe_mode) {
        tmp_rec.set_safe_interpolation();

        if (scan) {
          WARNING_MSG("Disabled shift scan in safe mode!");
          scan = false;
        }
      }

      if (scan) {
        scan_for_shift();

#ifndef WITH_MPI
        uint32_t start = thr_n * bunch_size;
#else
        uint32_t start = 6 * static_cast<uint32_t>(world_size) * bunch_size; // start even more?
#endif

        start_probe_jobs(std::vector<uint32_t> (n - 1, 1), start);
        started_probes.emplace(std::vector<uint32_t> (n - 1, 1), start);

#ifdef WITH_MPI
        cond_val.notify_one();
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

      run_until_done();
    }

    tp.kill_all();

    if (verbosity > SILENT) {
      if (one_done || one_new_prime) {
        INFO_MSG("Probe: " + std::to_string(iteration) +
                 " | Done: " + std::to_string(items_done) + " / " + std::to_string(items) +
                 " | " + "Needs new prime field: " + std::to_string(items_new_prime) + " / " + std::to_string(items - items_done));
      }

      std::cout << "\n";
      INFO_MSG("Completed reconstruction in " +
               std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count()) + " s | " + std::to_string(total_iterations) + " probes in total");
      INFO_MSG("Needed prime fields: " + std::to_string(prime_it) + " + 1");
      INFO_MSG("Average time of the black-box probe: " + std::to_string(average_black_box_time) + " s");
    }
  }

  std::vector<RationalFunction> Reconstructor::get_result() {
    std::vector<RationalFunction> result {};

    for (auto & rec : reconst) {
      if (std::get<2>(rec) == DONE) {
        result.emplace_back(std::get<3>(rec)->get_result());
      }
    }

    return result;
  }

  std::vector<std::pair<std::string, RationalFunction>> Reconstructor::get_early_results() {
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

  void Reconstructor::scan_for_shift() {
    if (verbosity > SILENT)
      INFO_MSG("Scanning for a sparse shift");

    // Generate all possible combinations of shifting variables
    const auto shift_vec = generate_possible_shifts(n);

    bool first = true;
    bool found_shift = false;
    uint32_t counter = 0;
    uint32_t bound = shift_vec.size();

    tmp_rec.scan_for_sparsest_shift();

    start_first_runs();

    uint32_t max_deg_num = 0;
    uint32_t max_deg_den = 0;

    // Run this loop until a proper shift is found
    while (!found_shift && counter != bound) {
      if (!first) {
        tmp_rec.set_zi_shift(shift_vec[counter]);
        shift = tmp_rec.get_zi_shift_vec();

#ifndef WITH_MPI
        uint32_t start = thr_n * bunch_size;
#else
        uint32_t start = 6 * static_cast<uint32_t>(world_size) * bunch_size; // start even more?
#endif

        start_probe_jobs(std::vector<uint32_t> (n - 1, 1), start);
        started_probes.emplace(std::vector<uint32_t> (n - 1, 1), start);

#ifdef WITH_MPI
        cond_val.notify_one();
#endif
      }

      run_until_done();

      // Kill all jobs
      // otherwise it can happen that a RatReconst is fed with old data
#ifdef WITH_MPI
      {
        std::unique_lock<std::mutex> lock_val(mut_val);

        value_queue = std::queue<std::vector<uint64_t>>();
        new_jobs = false;

        index_map.clear();
        ind = 0;

        cond_val.wait(lock_val);
      }
#endif

      tp.kill_all();

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

        if (verbosity > SILENT) {
          INFO_MSG("Maximal degree of numerator: " + std::to_string(max_deg_num) + " | Maximal degree of denominator: " + std::to_string(max_deg_den));
        }
      } else {
        ++counter;
      }

      probes.clear();
      finished_probes_it = std::queue<future_list::iterator>();
      probes_bunch.clear();
      finished_probes_bunch_it = std::queue<future_list_bunch::iterator>();
      bunch_t.clear();
      bunch.clear();
      jobs_finished = 0;
      started_probes.clear();
      fed_ones = 0;
      feed_jobs = 0;
      interpolate_jobs = 0;
      iteration = 0;
      items_done = 0;
      done = false;

#ifdef WITH_MPI
      probes_queued = 0;

      {
        std::unique_lock<std::mutex> lock(future_control);

        results_queue = std::queue<std::pair<uint64_t, std::vector<FFInt>>>();
      }
#endif
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

    if (verbosity > SILENT) {
      if (found_shift) {
        std::string msg = "";

        for (const auto & el : shift_vec[counter - 1]) {
          msg += std::to_string(el) + ", ";
        }

        msg = msg.substr(0, msg.size() - 2);
        INFO_MSG("Found a sparse shift after " + std::to_string(counter + 1) + " scans");
        INFO_MSG("Shift the variable tuple (" + msg + ")");
      } else {
        INFO_MSG("Found no sparse shift after " + std::to_string(counter + 1) + " scans");
      }

      INFO_MSG("Completed scan in " + std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - prime_start).count()) +
               " s | " + std::to_string(total_iterations) + " probes in total");
      INFO_MSG("Average time of the black-box probe: " + std::to_string(average_black_box_time) + " s");
      std::cout << "\n";
      INFO_MSG("Proceeding with interpolation over prime field F(" + std::to_string(primes()[prime_it]) + ")");
    }

    prime_start = std::chrono::high_resolution_clock::now();
  }

  void Reconstructor::start_first_runs() {
    prime_start = std::chrono::high_resolution_clock::now();
    shift = tmp_rec.get_zi_shift_vec();
    std::vector<uint32_t> zi_order(n - 1, 1);

#ifdef WITH_MPI
    send_first_jobs();
#else
    uint32_t start = thr_n * bunch_size;

    start_probe_jobs(zi_order, start);
    started_probes.emplace(zi_order, start);
#endif

    FFInt t = 1;
    std::vector<FFInt>* probe = new std::vector<FFInt>;

#ifdef WITH_MPI
    ++started_probes[zi_order];

    std::vector<FFInt> values(n);

    // check if t was already used for this zi_order
    {
      std::unique_lock<std::mutex> chosen_lock(chosen_mutex);

      auto it = chosen_t.find(zi_order);

      if (it != chosen_t.end()) {
        while (true) {
          t = tmp_rec.get_rand();

          auto itt = it->second.find(t.n);

          if (itt != it->second.end()) {
            std::unique_lock<std::mutex> lock_print(print_control);

            WARNING_MSG("Found a duplicate of t, choosing a new one");
          } else {
            it->second.emplace(t.n);
            break;
          }
        }
      } else {
        t = tmp_rec.get_rand();
        chosen_t.emplace(std::make_pair(zi_order, std::unordered_set<uint64_t>( {t.n})));
      }
    }

    std::vector<FFInt> rand_zi = tmp_rec.get_rand_zi_vec(zi_order, false);

    values[0] = t + shift[0];

    for (uint32_t i = 1; i != n; ++i) {
      values[i] = rand_zi[i - 1] * t + shift[i];
    }

    *probe = bb(values);
#else
    get_probe(t, zi_order, probe, average_black_box_time);
#endif

    ++iteration;

    if (verbosity > SILENT) {
      INFO_MSG("Time for the first black-box probe: " + std::to_string(average_black_box_time) + " s");
    }

    ++fed_ones;
    items = probe->size();
    size_t tag_size = tags.size();

#ifdef WITH_MPI
    {
      std::unique_lock<std::mutex> lock_val(mut_val);

      new_jobs = true;

      cond_val.notify_one();
    }
#endif

    if (tag_size != 0 && tag_size != items) {
      ERROR_MSG("Number of tags does not match the black box!");
      std::exit(EXIT_FAILURE);
    }

    ogzstream file;

    if (save_states) {
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

      file << (t + shift[0]).n << " ";

      for (uint32_t i = 1; i != n; ++i) {
        file << (rand_zi[i - 1] * t + shift[i]).n << " ";
      }

      file << "\n";
    }

    for (uint32_t i = 0; i != items; ++i) {
      RatReconst* rec = new RatReconst(n);

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

        file << ((*probe)[i]).n << "\n";
      }

      rec->feed(t, (*probe)[i], zi_order, prime_it);
      std::tuple<bool, bool, uint32_t> interpolated_done_prime = rec->interpolate();

      if (std::get<2>(interpolated_done_prime) > prime_it) {
        ++items_new_prime;
      }

      std::mutex* mut = new std::mutex;

      reconst.emplace_back(std::make_tuple(i, mut, RECONSTRUCTING, rec));
    }

    if (save_states) {
      file.close();
      tags.clear();
    }

    delete probe;

    if (verbosity > SILENT) {
      INFO_MSG("Probe: 1 | Done: 0 / " + std::to_string(items) + " | Needs new prime field: " + std::to_string(items_new_prime) + " / " + std::to_string(items));
    }

    start_probe_jobs(zi_order, bunch_size);
    started_probes[zi_order] += bunch_size;
  }

  void Reconstructor::run_until_done() {
    FFInt t = 1;
    std::vector<uint32_t> zi_order(n - 1, 1);

    new_prime = false;

    if (resume_from_state) {
      if (prime_it == 0 && /*items_new_prime != items*/items != items_new_prime + items_done) {
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

        // TODO don't start as much ones
        start_probe_jobs(std::vector<uint32_t>(n - 1, 1), thr_n * bunch_size);
        started_probes.emplace(std::vector<uint32_t>(n - 1, 1), thr_n * bunch_size);
        /*shift = tmp_rec.get_zi_shift_vec();

        uint32_t start = thr_n * bunch_size;
        start_probe_jobs(std::vector<uint32_t>(n - 1, 1), start);
        started_probes.emplace(std::vector<uint32_t>(n - 1, 1), start);*/
      } else {
        new_prime = true;
      }
    }

    while (!done) {
      if (new_prime) {
#ifdef WITH_MPI
        {
          std::unique_lock<std::mutex> lock_val(mut_val);

          value_queue = std::queue<std::vector<uint64_t>>();
          new_jobs = false;

          index_map.clear();
          ind = 0;

          cond_val.wait(lock_val);
        }
#endif

        tp.kill_all();

        clean_reconst();

        if (save_states) {
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

        /*if (save_states && prime_it == 0) {
          std::remove("ff_save/scan");
        }*/

        total_iterations += iteration;
        ++prime_it;

        if (verbosity > SILENT) {
          if (one_done || one_new_prime) {
            INFO_MSG("Probe: " + std::to_string(iteration) +
                     " | Done: " + std::to_string(items_done) + " / " + std::to_string(items) +
                     " | " + "Needs new prime field: " + std::to_string(items_new_prime) + " / " + std::to_string(items - items_done));
          }

          INFO_MSG("Completed current prime field in " +
                   std::to_string(std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - prime_start).count()) +
                   " s | " + std::to_string(total_iterations) + " probes in total");
          INFO_MSG("Average time of the black-box probe: " + std::to_string(average_black_box_time) + " s");
          std::cout << "\n";
          INFO_MSG("Promote to new prime field: F(" + std::to_string(primes()[prime_it]) + ")");
        }

        prime_start = std::chrono::high_resolution_clock::now();

        iteration = 0;

        probes.clear();
        finished_probes_it = std::queue<future_list::iterator>();
        probes_bunch.clear();
        finished_probes_bunch_it = std::queue<future_list_bunch::iterator>();
        bunch_t.clear();
        bunch.clear();
        fed_ones = 0;
        jobs_finished = 0;
        started_probes.clear();

        // Only reset chosen_t when not resuming from a saved state
        if (!set_anchor_points) {
          chosen_t.clear();
        }

        feed_jobs = 0;
        interpolate_jobs = 0;
        new_prime = false;
        items_new_prime = 0;
        one_done = false;
        one_new_prime = false;

        FFInt::set_new_prime(primes()[prime_it]);
        bb.prime_changed();

#ifdef WITH_MPI
        probes_queued = 0;

        {
          std::unique_lock<std::mutex> lock(future_control);

          results_queue = std::queue<std::pair<uint64_t, std::vector<FFInt>>>();
        }
#endif

        // if only a small constant is reconstructed it will not ask for new run
        if (probes_for_next_prime == 0) {
#ifdef WITH_MPI
          probes_for_next_prime = static_cast<uint32_t>(world_size) * bunch_size; // start even more?
#else
          probes_for_next_prime = thr_n * bunch_size;
#endif
        }

        bool shift_disabled = false;

        if (!safe_mode && (!save_states || (save_states && !set_anchor_points)) && prime_it >= min_prime_keep_shift && !tmp_rec.need_shift()) {
          if (tmp_rec.get_zi_shift_vec() != std::vector<FFInt> (n, 0)) {
            if (verbosity > SILENT)
              INFO_MSG("Disable shift");

            tmp_rec.disable_shift();
            shift_disabled = true;
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
            tmp_rec.set_anchor_points(parse_vector_FFInt(line, n));
          } else {
            ERROR_MSG("Anchor point file not found!");
            std::exit(EXIT_FAILURE);
          }

          anchor_point_file.close();

          if (!shift_disabled) {
            std::ifstream shift_file;
            shift_file.open("ff_save/shift");

            if (shift_file.is_open()) {
              std::getline(shift_file, line);
              tmp_rec.set_shift(parse_vector_FFInt(line, n));
            } else {
              ERROR_MSG("Shift file not found!");
              std::exit(EXIT_FAILURE);
            }
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
        if (probes_for_next_prime > thr_n * bunch_size) {
          uint32_t start = thr_n * bunch_size;

          if (verbosity == CHATTY) {
            INFO_MSG("Starting " + std::to_string(start) + " jobs now, the remaining " + std::to_string(probes_for_next_prime - start) + " jobs will be started later");
          }

          start_probe_jobs(std::vector<uint32_t>(n - 1, 1), start);
          started_probes.emplace(std::vector<uint32_t>(n - 1, 1), start);
        } else {
          uint32_t start = (probes_for_next_prime + bunch_size - 1) / bunch_size * bunch_size;

          if (verbosity == CHATTY) {
            INFO_MSG("Starting " + std::to_string(start) + " jobs");
          }

          start_probe_jobs(std::vector<uint32_t>(n - 1, 1), start);
          started_probes.emplace(std::vector<uint32_t>(n - 1, 1), start);
        }

#ifdef WITH_MPI
        cond_val.notify_one();
#endif

        probes_for_next_prime = 0;
      }

      std::vector<FFInt>* probe = new std::vector<FFInt>;
      double time;

      get_probe(t, zi_order, probe, time);

      ++iteration;

      average_black_box_time = (average_black_box_time * (total_iterations + iteration - 1) + time) / (total_iterations + iteration);

      if ((prime_it == 0 || safe_mode == true) && zi_order == std::vector<uint32_t>(n - 1, 1)) {
        std::unique_lock<std::mutex> lock_status(job_control);

        ++fed_ones;
      }

      {
        std::unique_lock<std::mutex> lock_feed(feed_control);

        ++feed_jobs;
      }

      tp.run_priority_task([this, zi_order, t, probe]() {
        feed_job(zi_order, t, probe);
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
        std::unique_lock<std::mutex> lock_future(future_control);

        if (jobs_finished == 0 && ((bunch_size == 1 && probes.empty()) || (bunch_size != 1 && bunch.empty() && probes_bunch.empty()))) {
          lock_future.unlock();

          {
            std::unique_lock<std::mutex> lock_feed(feed_control);

            while (feed_jobs > 0 || interpolate_jobs > 0) {
              condition_feed.wait(lock_feed);

              lock_feed.unlock();
              lock_future.lock();

#ifndef WITH_MPI

              if (jobs_finished > 0 || !probes.empty() || !bunch.empty() || !probes_bunch.empty()) {
#else
              std::unique_lock<std::mutex> lck(mut_val);

              if (probes_queued != 0) {
#endif
                lock_future.unlock();

                break;
              }

              lock_future.unlock();
              lock_feed.lock();
            }

            lock_future.lock();

#ifndef WITH_MPI

            if (jobs_finished > 0 || !probes.empty() || !bunch.empty() || !probes_bunch.empty()) {
#else
            std::unique_lock<std::mutex> lck(mut_val);

            if (probes_queued != 0) {
#endif
              continue;
            }

            lock_future.unlock();
          }

          // no jobs are running anymore, check if done or new_prime else throw error
          if (items_done == items) {
            done = true;
            continue;
          } else if (items_done + items_new_prime == items) {
            new_prime = true;
            continue;
          } else {
#ifndef MPI
            throw std::runtime_error("Nothing left to feed: " + std::to_string(iteration) + " | " + std::to_string(items) + " " + std::to_string(items_new_prime) + " " + std::to_string(items_done) + " | " + std::to_string(jobs_finished) + " " + std::to_string(probes.empty()) + " " + std::to_string(bunch.empty()) + " " + std::to_string(probes_bunch.empty()) + " | " + std::to_string(feed_jobs) + " " + std::to_string(interpolate_jobs) + "\n");
#else
            throw std::runtime_error("Nothing left to feed: " + std::to_string(items) + " " + std::to_string(items_new_prime) + " " + std::to_string(items_done) + " | " + std::to_string(jobs_finished) + " " + std::to_string(probes.empty()) + " " + std::to_string(bunch.empty()) + " " + std::to_string(probes_bunch.empty()) + " | " + std::to_string(feed_jobs) + " " + std::to_string(interpolate_jobs) + " | " + std::to_string(iteration) + " | " + std::to_string(probes_queued) + " " + std::to_string(results_queue.size()) + " " + std::to_string(value_queue.size()) + "\n");
#endif
          }
        }
      }
    }

    if (!scan && save_states) {
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

    total_iterations += iteration;
  }

  void Reconstructor::start_probe_jobs(const std::vector<uint32_t>& zi_order, const uint32_t to_start) {
    bool ones = false;

    if ((prime_it == 0 || safe_mode == true) && zi_order == std::vector<uint32_t> (n - 1, 1)) {
      ones = true;
    }

    std::vector<FFInt> rand_zi;
    rand_zi.reserve(zi_order.size());

    if (!ones && (prime_it != 0 && safe_mode == false)) {
      rand_zi = tmp_rec.get_rand_zi_vec(zi_order, true);
    } else {
      rand_zi = tmp_rec.get_rand_zi_vec(zi_order, false);
    }

    if (bunch_size == 1) {
#ifndef WITH_MPI
      std::vector<FFInt> values(n);
#endif
      FFInt t;

      for (uint32_t j = 0; j != to_start; ++j) {
#ifdef WITH_MPI
        std::vector<uint64_t> values(n + 1);
#endif
        t = tmp_rec.get_rand();

        // check if t was already used for this zi_order
        {
          std::unique_lock<std::mutex> chosen_lock(chosen_mutex);

          auto it = chosen_t.find(zi_order);

          if (it != chosen_t.end()) {
            auto itt = it->second.find(t.n);

            if (itt != it->second.end()) {
              --j;

              std::unique_lock<std::mutex> lock_print(print_control);

              //WARNING_MSG("Found a duplicate of t, choosing a new one");

              continue;
            } else {
              it->second.emplace(t.n);
            }
          } else {
            chosen_t.emplace(std::make_pair(zi_order, std::unordered_set<uint64_t>( {t.n})));
          }
        }

#ifdef WITH_MPI
        values[1] = (t + shift[0]).n;

        for (uint32_t i = 1; i != n; ++i) {
          values[i + 1] = (rand_zi[i - 1] * t + shift[i]).n;
        }

        std::unique_lock<std::mutex> lock_val(mut_val);

        values[0] = ind;

        index_map.emplace(std::make_pair(ind, std::make_pair(t, zi_order)));
        ++ind;

        value_queue.emplace(std::move(values));
#else
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

            ++jobs_finished;
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

            ++jobs_finished;
            finished_probes_it.emplace(it);

            condition_future.notify_one();

            return std::make_pair(std::move(probe), std::chrono::duration<double>(time1 - time0).count());
          });

          std::get<2>(*it) = std::move(future);
        }

#endif
      }

#ifdef WITH_MPI
      std::unique_lock<std::mutex> lock_val(mut_val);

      probes_queued += static_cast<uint64_t>(to_start);
      new_jobs = true;
#endif
    } else {
      for (uint32_t j = 0; j != to_start / bunch_size; ++j) {
        std::vector<FFInt> t_vec;
        t_vec.reserve(bunch_size);
        std::vector<std::vector<FFInt>> values_vec;
        values_vec.reserve(bunch_size);

        for (uint32_t i = 0; i != bunch_size; ++i) {
          std::vector<FFInt> values(n);

          FFInt t = tmp_rec.get_rand();

          // check if t was already used for this zi_order
          {
            std::unique_lock<std::mutex> chosen_lock(chosen_mutex);

            auto it = chosen_t.find(zi_order);

            if (it != chosen_t.end()) {
              auto itt = it->second.find(t.n);

              if (itt != it->second.end()) {
                --i;

                std::unique_lock<std::mutex> lock_print(print_control);

                WARNING_MSG("Found a duplicate of t, choosing a new one");

                continue;
              } else {
                it->second.emplace(t.n);
              }
            } else {
              chosen_t.emplace(std::make_pair(zi_order, std::unordered_set<uint64_t>( {t.n})));
            }
          }

          values[0] = t + shift[0];

          for (uint32_t i = 1; i != n; ++i) {
            values[i] = rand_zi[i - 1] * t + shift[i];
          }

          t_vec.emplace_back(t);
          values_vec.emplace_back(std::move(values));
        }

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

            jobs_finished += bunch_size;
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

            jobs_finished += bunch_size;
            finished_probes_bunch_it.emplace(it);

            condition_future.notify_one();

            return std::make_pair(std::move(probe_vec), std::chrono::duration<double>(time1 - time0).count());
          });

          std::get<2>(*it) = std::move(future);
        }
      }
    }
  }

  void Reconstructor::get_probe(FFInt& t, std::vector<uint32_t>& zi_order, std::vector<FFInt>* probe, double& time) {
    std::unique_lock<std::mutex> lock_future(future_control);

    while (jobs_finished == 0) {
      condition_future.wait(lock_future);
    }

    if (bunch_size == 1) {
#ifdef WITH_MPI

      if (!finished_probes_it.empty()) {
#endif
        auto it = finished_probes_it.front();
        finished_probes_it.pop();

        //auto start = std::chrono::high_resolution_clock::now();
        (std::get<2>(*it)).wait();
        //auto end = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double, std::micro> elapsed = end-start;
        //std::cout << "Waited " << elapsed.count() << " micro s\n";

        t = std::get<0>(*it);
        zi_order = std::get<1>(*it);
        std::pair<std::vector<FFInt>, double> tmp = (std::get<2>(*it)).get();
        *probe = std::move(tmp.first);
        time = tmp.second;
        probes.erase(it);
#ifdef WITH_MPI
      } else {
        auto tmp = std::move(results_queue.front());
        results_queue.pop();
        *probe = std::move(tmp.second);

        std::unique_lock<std::mutex> lck(mut_val);

        auto tmp2 = std::move(index_map[tmp.first]);
        t = tmp2.first;
        zi_order = std::move(tmp2.second);
        --probes_queued;
      }

#endif
    } else {
      if (bunch.empty()) {
        auto it = finished_probes_bunch_it.front();
        finished_probes_bunch_it.pop();

        //auto start = std::chrono::high_resolution_clock::now();
        (std::get<2>(*it)).wait();
        //auto end = std::chrono::high_resolution_clock::now();
        //std::chrono::duration<double, std::micro> elapsed = end-start;
        //std::cout << "Waited " << elapsed.count() << " micro s\n";

        bunch_t = std::get<0>(*it);
        bunch_zi_order = std::get<1>(*it);
        std::pair<std::vector<std::vector<FFInt>>, double> tmp = (std::get<2>(*it)).get();
        bunch = std::move(tmp.first);
        bunch_time = tmp.second;
        probes_bunch.erase(it);
      }

      t = bunch_t.back();
      zi_order = bunch_zi_order;
      *probe = std::move(bunch.back());
      time = bunch_time / bunch_size;

      bunch_t.pop_back();
      bunch.pop_back();
    }

    --jobs_finished;
  }

  void Reconstructor::feed_job(const std::vector<uint32_t> zi_order, const firefly::FFInt t, std::vector<FFInt>* probe) {
    {
      std::unique_lock<std::mutex> lock(feed_control);

      interpolate_jobs += items;
    }

    uint32_t counter = 0;

    uint64_t tmp = 0;

    for (auto & rec : reconst) {
      std::unique_lock<std::mutex> lock_exists(*(std::get<1>(rec)));

      if (std::get<2>(rec) == RECONSTRUCTING) {
        lock_exists.unlock();

        std::pair<bool, uint32_t> done_prime = std::get<3>(rec)->get_done_and_prime();

        if (!done_prime.first) {
          if (done_prime.second == prime_it) {
            auto interpolate_and_write = std::get<3>(rec)->feed(t, (*probe)[std::get<0>(rec)], zi_order, prime_it);

            if (interpolate_and_write.first) {
              ++counter;

              tp.run_priority_task([this, &rec]() {
                interpolate_job(rec);
              });
            }

            if (interpolate_and_write.second) {
              tp.run_priority_task([&rec]() {
                std::get<3>(rec)->write_food_to_file();
              });
            }
          }
        }
      }

      ++tmp;
    }

    delete probe;

    {
      std::unique_lock<std::mutex> lock_status(status_control);

      if (!scan && (one_done || one_new_prime)) {
        one_done = false;
        one_new_prime = false;

        std::unique_lock<std::mutex> lock_print(print_control);

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

  void Reconstructor::interpolate_job(RatReconst_tuple& rec) {
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
                  std::vector<uint32_t> zi_order(n - 1, counter);
                  ++counter;

                  std::unique_lock<std::mutex> lock(job_control);

                  auto it = started_probes.find(zi_order);

                  if (it != started_probes.end()) {
                    if (required_probes > started_probes[zi_order]) {
                      uint32_t start = required_probes - started_probes[zi_order];

                      if (bunch_size != 1) {
                        start = (start + bunch_size - 1) / bunch_size * bunch_size;
                      }

                      started_probes[zi_order] = required_probes;

                      lock.unlock();

                      if (verbosity == CHATTY) {
                        std::string msg = "Starting zi_order (";

                        for (const auto & ele : zi_order) {
                          msg += std::to_string(ele) + ", ";
                        }

                        msg = msg.substr(0, msg.length() - 2);
                        msg += ") " + std::to_string(start) + " time(s)";

                        std::unique_lock<std::mutex> lock_print(print_control);

                        INFO_MSG(msg);
                      }

                      start_probe_jobs(zi_order, start);
                    }
                  } else {
                    if (bunch_size != 1) {
                      required_probes = (required_probes + bunch_size - 1) / bunch_size * bunch_size;
                    }

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

                    start_probe_jobs(zi_order, required_probes);
                  }
                }
              }
            }
          } else {
            std::vector<uint32_t> zi_order = std::get<3>(rec)->get_zi_order();

            if ((prime_it == 0 || safe_mode == true) && zi_order == std::vector<uint32_t>(n - 1, 1)) {
              //lock_exists.unlock();
              std::unique_lock<std::mutex> lock(job_control);

              if (started_probes[zi_order] - thr_n * bunch_size <= fed_ones - 1) {
                uint32_t start = fed_ones - started_probes[zi_order] + thr_n * bunch_size;

                if (bunch_size != 1) {
                  start = (start + bunch_size - 1) / bunch_size * bunch_size;
                }

                started_probes[zi_order] += start;

                lock.unlock();

                if (verbosity == CHATTY) {
                  std::unique_lock<std::mutex> lock_print(print_control);

                  INFO_MSG("Starting ones: " + std::to_string(start));
                }

                start_probe_jobs(zi_order, start);
              }
            } else {
              uint32_t required_probes = std::get<3>(rec)->get_num_eqn();

              //lock_exists.unlock();
              std::unique_lock<std::mutex> lock(job_control);

              auto it = started_probes.find(zi_order);

              if (it != started_probes.end()) {
                if (required_probes > started_probes[zi_order]) {
                  uint32_t start = required_probes - started_probes[zi_order];

                  if (bunch_size != 1) {
                    start = (start + bunch_size - 1) / bunch_size * bunch_size;
                  }

                  started_probes[zi_order] = required_probes;

                  lock.unlock();

                  if (verbosity == CHATTY) {
                    std::string msg = "Starting zi_order (";

                    for (const auto & ele : zi_order) {
                      msg += std::to_string(ele) + ", ";
                    }

                    msg = msg.substr(0, msg.length() - 2);
                    msg += ") " + std::to_string(start) + " time(s)";

                    std::unique_lock<std::mutex> lock_print(print_control);

                    INFO_MSG(msg);
                  }

                  start_probe_jobs(zi_order, start);
                }
              } else {
                if (bunch_size != 1) {
                  required_probes = (required_probes + bunch_size - 1) / bunch_size * bunch_size;
                }

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

                start_probe_jobs(zi_order, required_probes);
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

  void Reconstructor::clean_reconst() {
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

#ifdef WITH_MPI
  void Reconstructor::mpi_setup() {
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Bcast(&prime_it, 1, MPI_UINT32_T, master, MPI_COMM_WORLD);
  }

  void Reconstructor::send_first_jobs() {
    std::vector<uint32_t> zi_order = std::vector<uint32_t>(n - 1, 1);
    started_probes.emplace(std::make_pair(zi_order, 0));

    std::vector<FFInt> rand_zi = tmp_rec.get_rand_zi_vec(std::vector<uint32_t>(n - 1, 1), true);

    for (int i = 1; i != world_size; ++i) {
      uint64_t start;
      MPI_Recv(&start, 1, MPI_UINT64_T, i, RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      //std::cout << "starting " << start << " jobs on worker " << i << "\n";

      nodes.emplace(std::make_pair(i, start));

      probes_queued += start;
      started_probes[zi_order] += start;
      uint64_t values[static_cast<uint32_t>(start) * (n + 1)];

      for (uint32_t ii = 0; ii != static_cast<uint32_t>(start); ++ii) {
        values[ii * (n + 1)] = ind;

        FFInt t = tmp_rec.get_rand();

        // check if t was already used for this zi_order
        auto it = chosen_t.find(zi_order);

        if (it != chosen_t.end()) {
          auto itt = it->second.find(t.n);

          if (itt != it->second.end()) {
            --ii;

            std::unique_lock<std::mutex> lock_print(print_control);

            WARNING_MSG("Found a duplicate of t, choosing a new one");

            continue;
          } else {
            it->second.emplace(t.n);
          }
        } else {
          chosen_t.emplace(std::make_pair(zi_order, std::unordered_set<uint64_t>( {t.n})));
        }

        values[ii * (n + 1) + 1] = (t + shift[0]).n;

        for (uint32_t j = 1; j != n; ++j) {
          values[ii * (n + 1) + j + 1] = (t * rand_zi[j - 1] + shift[j]).n;
        }

        index_map.emplace(std::make_pair(ind, std::make_pair(t, std::vector<uint32_t>(n - 1, 1))));
        ++ind;
      }

      MPI_Send(&values, static_cast<int>(static_cast<uint32_t>(start) * (n + 1)), MPI_UINT64_T, i, VALUES, MPI_COMM_WORLD);
    }
  }

  void Reconstructor::mpi_communicate() {
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
        MPI_Request requests[world_size - 1];

        for (int i = 1; i != world_size; ++i) {
          uint64_t tmp;
          MPI_Isend(&tmp, 1, MPI_UINT64_T, i, END, MPI_COMM_WORLD, &requests[i - 1]);
        }

        MPI_Waitall(world_size - 1, requests, MPI_STATUSES_IGNORE);

        int flag = 1;
        MPI_Status status;

        while (flag) {
          MPI_Iprobe(MPI_ANY_SOURCE, RESULT, MPI_COMM_WORLD, &flag, &status);

          if (flag) {
            int amount;
            MPI_Get_count(&status, MPI_UINT64_T, &amount);

            uint64_t tmp[amount];
            MPI_Recv(tmp, amount, MPI_UINT64_T, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
          }
        }

        break;
      } else if (new_prime || done && scan) {
        MPI_Request requests[world_size - 1];

        //std::cout << "com np\n";

        uint64_t prime_tmp = static_cast<uint64_t>(prime_it);

        if (!scan) {
          ++prime_tmp;
        }

        for (int i = 1; i != world_size; ++i) {
          MPI_Isend(&prime_tmp, 1, MPI_UINT64_T, i, NEW_PRIME, MPI_COMM_WORLD, &requests[i - 1]);
        }

        MPI_Waitall(world_size - 1, requests, MPI_STATUSES_IGNORE);

        int flag = 1;
        MPI_Status status;

        while (flag) {
          MPI_Iprobe(MPI_ANY_SOURCE, RESULT, MPI_COMM_WORLD, &flag, &status);

          if (flag) {
            int amount;
            MPI_Get_count(&status, MPI_UINT64_T, &amount);

            //std::cout << "garbage recieved: " << (static_cast<uint32_t>(amount) - 1) / (items + 1) << "\n";

            uint64_t tmp[amount];
            MPI_Recv(tmp, amount, MPI_UINT64_T, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
          }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        //std::cout << "com rec\n";

        cond_val.notify_one();

        std::unique_lock<std::mutex> lock_val(mut_val);

        while (!new_jobs) {
          cond_val.wait(lock_val);
        }

        empty_nodes = std::queue<std::pair<int, uint64_t>>();

        //std::cout << "com through\n";

        for (int i = 1; i != world_size; ++i) {
          uint64_t free_slots;
          MPI_Recv(&free_slots, 1, MPI_UINT64_T, i, RESULT, MPI_COMM_WORLD, &status);

          uint64_t size = std::min(free_slots, static_cast<uint64_t>(value_queue.size()));

          //std::cout << "sending " << size << " jobs\n";

          if (size != 0) {
            uint64_t values[size * (n + 1)];

            for (uint64_t i = 0; i != size; ++i) {
              values[i * (n + 1)] = value_queue.front()[0];

              for (uint64_t j = 0; j != n; ++j) {
                values[i * (n + 1) + 1 + j] = value_queue.front()[1 + j];
              }

              value_queue.pop();
            }

            MPI_Send(&values, static_cast<int>(size * (n + 1)), MPI_UINT64_T, status.MPI_SOURCE, VALUES, MPI_COMM_WORLD);
          } else if (free_slots == nodes[status.MPI_SOURCE]) {
            //std::cout << "com np empty nodes " << status.MPI_SOURCE << " " << free_slots << "\n";
            empty_nodes.emplace(std::make_pair(status.MPI_SOURCE, free_slots));
          }
        }
      } else if (restart_empty_nodes) {
        std::unique_lock<std::mutex> lock_val(mut_val);

        while (!empty_nodes.empty() && !value_queue.empty()) {
          auto node = empty_nodes.front();
          empty_nodes.pop();

          uint64_t available_jobs = static_cast<uint64_t>(value_queue.size());
          uint64_t size = std::min(node.second, available_jobs);
          uint64_t values[size * (n + 1)];

          for (uint64_t i = 0; i != size; ++i) {
            values[i * (n + 1)] = value_queue.front()[0];

            for (uint64_t j = 0; j != n; ++j) {
              values[i * (n + 1) + 1 + j] = value_queue.front()[1 + j];
            }

            value_queue.pop();
          }

          //std::cout << "restart: sending " << size << " jobs\n";

          MPI_Send(&values, static_cast<int>(size * (n + 1)), MPI_UINT64_T, node.first, VALUES, MPI_COMM_WORLD);
        }

        if (value_queue.empty()) {
          new_jobs = false;
        }
      } else {
        int amount;
        MPI_Get_count(&status, MPI_UINT64_T, &amount);

        if ((static_cast<uint32_t>(amount) - 1) % (items + 1) != 0) {
          ERROR_MSG("Corrupted results recieved: " + std::to_string(amount - 1));
          std::exit(EXIT_FAILURE);
        }

        uint32_t new_results = (static_cast<uint32_t>(amount) - 1) / (items + 1);

        //std::cout << "comm recieving " << new_results << " results \n";

        uint64_t results_list[amount];
        MPI_Recv(results_list, amount, MPI_UINT64_T, status.MPI_SOURCE, RESULT, MPI_COMM_WORLD, &status);

        for (uint32_t i = 0; i != new_results; ++i) {
          uint64_t index = results_list[i * (items + 1)];
          std::vector<FFInt> results;
          results.reserve(items);

          for (uint32_t j = 1; j != items + 1; ++j) {
            results.emplace_back(results_list[i * (items + 1) + j]);
          }

          std::unique_lock<std::mutex> lock_res(future_control);

          results_queue.emplace(std::make_pair(index, std::move(results)));
        }

        {
          std::unique_lock<std::mutex> lock_res(future_control);

          jobs_finished += new_results;

          condition_future.notify_one();
        }

        uint64_t free_slots = results_list[amount - 1];

        std::unique_lock<std::mutex> lock_val(mut_val);

        uint64_t size = std::min(free_slots, static_cast<uint64_t>(value_queue.size()));

        if (size != 0) {
          uint64_t values[size * (n + 1)];

          for (uint64_t i = 0; i != size; ++i) {
            values[i * (n + 1)] = value_queue.front()[0];

            for (uint64_t j = 0; j != n; ++j) {
              values[i * (n + 1) + 1 + j] = value_queue.front()[1 + j];
            }

            value_queue.pop();
          }

          if (value_queue.empty()) {
            new_jobs = false;
          }

          lock_val.unlock();

          //std::cout << "sending " << size << " jobs\n";

          MPI_Send(&values, static_cast<int>(size * (n + 1)), MPI_UINT64_T, status.MPI_SOURCE, VALUES, MPI_COMM_WORLD);
        } else if (free_slots == nodes[status.MPI_SOURCE]) {
          //std::cout << "empty node " << status.MPI_SOURCE << " " << free_slots << "\n";
          empty_nodes.emplace(std::make_pair(status.MPI_SOURCE, free_slots));
        } else {
          MPI_Send(NULL, 0, MPI_UINT64_T, status.MPI_SOURCE, VALUES, MPI_COMM_WORLD);
        }
      }
    }
  }
#endif
}
