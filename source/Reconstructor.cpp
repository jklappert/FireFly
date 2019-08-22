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

#include "ParserUtils.hpp"
#include "Reconstructor.hpp"
#include "ReconstHelper.hpp"
#include "tinydir.h"
#include "utils.hpp"
#include "version.hpp"

#include <fstream>
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

  void Reconstructor::resume_from_saved_state(const std::string& directory) {
    tinydir_dir dir;
    tinydir_open_sorted(&dir, directory.c_str());

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
      return std::stoi(l.substr(0, l.find("_"))) < stoi(r.substr(0, r.find("_")));
    });

    for (const auto & file : files) {
      paths.emplace_back(directory + "/" + file);
    }

    if(paths.size() != 0) {
      resume_from_saved_state(paths);
    } else {
      ERROR_MSG("Directory " + directory + " does not exist or has no content.");
      std::exit(EXIT_FAILURE);
    }
  }

  void Reconstructor::resume_from_saved_state(const std::vector<std::string>& file_paths_) {
    INFO_MSG("Loading saved states");

    std::ifstream validation_file;
    validation_file.open("validation");
    std::string line;

    if (validation_file.is_open()) {
      std::getline(validation_file, line);
      std::vector<FFInt> values = parse_vector(line, "64");

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

    for (uint32_t i = 0; i != items; ++i) {
      RatReconst* rec = new RatReconst(n);
      std::pair<bool, uint32_t> shift_prime = rec->start_from_saved_file(file_paths[i]);

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

    if (scan) {
      if (prime_it == 0 && items_new_prime != items) {
        std::ifstream file;
        file.open("scan");

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
  }

  void Reconstructor::set_safe_interpolation() {
    safe_mode = true;
  }

  void Reconstructor::reconstruct() {
    start = std::chrono::high_resolution_clock::now();

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

        uint32_t start = thr_n * bunch_size;

        start_probe_jobs(std::vector<uint32_t> (n - 1, 1), start);
        started_probes.emplace(std::vector<uint32_t> (n - 1, 1), start);
      } else {
        start_first_runs();
      }
    } else {
      scan = false;
    }

    run_until_done();

    tp.kill_all();

    if (save_states) {
      std::remove("validation");
    }

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

        uint32_t start = thr_n * bunch_size;

        start_probe_jobs(std::vector<uint32_t> (n - 1, 1), start);
        started_probes.emplace(std::vector<uint32_t> (n - 1, 1), start);
      }

      run_until_done();

      // Kill all jobs
      // otherwise it can happen that a RatReconst is fed with old data
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
      file.open("scan");
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

    uint32_t start = thr_n * bunch_size;

    start_probe_jobs(zi_order, start);
    started_probes.emplace(zi_order, start);

    FFInt t = 1;
    std::vector<FFInt>* probe = new std::vector<FFInt>;

    get_probe(t, zi_order, probe, average_black_box_time);

    ++iteration;

    if (verbosity > SILENT) {
      INFO_MSG("Time for the first black-box probe: " + std::to_string(average_black_box_time) + " s");
    }

    ++fed_ones;
    items = probe->size();
    size_t tag_size = tags.size();

    if (tag_size != 0 && tag_size != items) {
      ERROR_MSG("Number of tags does not match the black box!");
      std::exit(EXIT_FAILURE);
    }

    std::ofstream file;

    if (save_states) {
      file.open("validation");

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

    bool done = false;
    bool new_prime = false;

    if (resume_from_state) {
      if (prime_it == 0 && items_new_prime != items) {
        INFO_MSG("Resuming in prime field: F(" + std::to_string(primes()[prime_it]) + ")");

        shift = tmp_rec.get_zi_shift_vec();

        uint32_t start = thr_n * bunch_size;
        start_probe_jobs(std::vector<uint32_t>(n - 1, 1), start);
        started_probes.emplace(std::vector<uint32_t>(n - 1, 1), start);
      } else {
        new_prime = true;
      }
    }

    while (!done) {
      if (new_prime) {
        //exit(-1);
        tp.kill_all();

        clean_reconst();

        if (save_states && prime_it) {
          mkdir("ff_save", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

          for (uint32_t item = 0; item != items; ++item) {
            std::string file_name_old = "ff_save/" + std::to_string(item) + "_" + std::to_string(prime_it - 1) + ".txt";
            std::string file_name_new = "ff_save/" + std::to_string(item) + "_" + std::to_string(prime_it) + ".txt";
            std::rename(file_name_old.c_str(), file_name_new.c_str());
          }
        }

        if (save_states && prime_it == 0) {
          std::remove("scan");
        }

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
                   " s | " + std::to_string(total_iterations) + " probes in total.");
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
        feed_jobs = 0;
        interpolate_jobs = 0;
        new_prime = false;
        items_new_prime = 0;
        one_done = false;
        one_new_prime = false;

        FFInt::set_new_prime(primes()[prime_it]);
        bb.prime_changed();

        // if only a small constant is reconstructed it will not ask for new run
        if (probes_for_next_prime == 0) {
          probes_for_next_prime = bunch_size;
        }

        if (!safe_mode && prime_it >= min_prime_keep_shift && !tmp_rec.need_shift()) {
          if (tmp_rec.get_zi_shift_vec() != std::vector<FFInt> (n, 0)) {
            if (verbosity > SILENT)
              INFO_MSG("Disable shift");

            tmp_rec.disable_shift();
          }
        }

        shift = tmp_rec.get_zi_shift_vec();
        tmp_rec.generate_anchor_points();

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

              if (jobs_finished > 0 || !probes.empty() || !bunch.empty() || !probes_bunch.empty()) {
                lock_future.unlock();

                break;
              }

              lock_future.unlock();
              lock_feed.lock();
            }

            lock_future.lock();

            if (jobs_finished > 0 || !probes.empty() || !bunch.empty() || !probes_bunch.empty()) {
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
            throw std::runtime_error("Nothing left to feed.");
          }
        }
      }
    }

    if (save_states && prime_it) {
      mkdir("ff_save", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

      for (uint32_t item = 0; item != items; ++item) {
        std::string file_name_old = "ff_save/" + std::to_string(item) + "_" + std::to_string(prime_it - 1) + ".txt";
        std::string file_name_new = "ff_save/" + std::to_string(item) + "_" + std::to_string(prime_it) + ".txt";
        std::rename(file_name_old.c_str(), file_name_new.c_str());
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
      std::vector<FFInt> values(n);
      FFInt t;

      for (uint32_t j = 0; j != to_start; ++j) {
        t = tmp_rec.get_rand();
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
      }
    } else {
      for (uint32_t j = 0; j != to_start / bunch_size; ++j) {
        std::vector<FFInt> t_vec;
        t_vec.reserve(bunch_size);
        std::vector<std::vector<FFInt>> values_vec;
        values_vec.reserve(bunch_size);

        for (uint32_t i = 0; i != bunch_size; ++i) {
          std::vector<FFInt> values(n);

          FFInt t = tmp_rec.get_rand();
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

    for (auto & rec : reconst) {
      std::unique_lock<std::mutex> lock_exists(*(std::get<1>(rec)));

      if (std::get<2>(rec) == RECONSTRUCTING) {
        lock_exists.unlock();

        std::pair<bool, uint32_t> done_prime = std::get<3>(rec)->get_done_and_prime();

        if (!done_prime.first) {
          if (done_prime.second == prime_it) {
            if (std::get<3>(rec)->feed(t, (*probe)[std::get<0>(rec)], zi_order, prime_it)) {
              ++counter;

              tp.run_task([this, &rec]() {
                interpolate_job(rec);
              });
            }
          }
        }
      }
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
}
