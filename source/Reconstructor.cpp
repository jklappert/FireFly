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

#include "Reconstructor.hpp"
#include "utils.hpp"
#include "version.hpp"

#include <chrono>

namespace firefly {
  Reconstructor::Reconstructor(uint32_t n_, uint32_t thr_n_, black_box_base * bb_, uint32_t verbosity_): n(n_), thr_n(thr_n_), bb(bb_), verbosity(verbosity_), tp(thr_n_) {
    if (verbosity > SILENT) {
      std::cout << "\nFire\033[1;32mFly\033[0m " << FireFly_VERSION_MAJOR << "." << FireFly_VERSION_MINOR << "." << FireFly_VERSION_RELEASE << "\n\n";
      INFO_MSG("Launching " << thr_n_ << " thread(s).");
    }

    FFInt::set_new_prime(primes()[prime_it]);
    uint64_t seed = static_cast<uint64_t> (std::time(0));
    BaseReconst().set_seed(seed);
    tmp_rec = RatReconst(n);
  }

  void Reconstructor::enable_scan() {
    scan = true;
  }

  void Reconstructor::set_tags() {
    save_states = true;
  }

  void Reconstructor::set_tags(const std::vector<std::string>& tags_) {
    save_states = true;
    tags = tags_;
  }

  void Reconstructor::resume_from_saved_state(const std::vector<std::string>& file_paths_) {
    resume_from_state = true;
    file_paths = file_paths_;
    items = file_paths.size();
    uint32_t counter = 0;

    for (uint32_t i = 0; i != items; ++i) {
      uint32_t old_it = prime_it;
      prime_it = std::max(prime_it, parse_prime_number(file_paths[i]));

      if (prime_it > old_it)
        counter = i;
    }

    tmp_rec.start_from_saved_file(file_paths[counter]);

    for (uint32_t i = 0; i != items; ++i) {
      RatReconst* rec = new RatReconst(n);
      rec->start_from_saved_file(file_paths[i]);

      probes_for_next_prime = std::max(probes_for_next_prime, rec->get_num_eqn());

      std::mutex* mut = new std::mutex;

      reconst.emplace_back(std::make_tuple(i, mut, DEFAULT, rec));
    }
  }

  uint32_t Reconstructor::parse_prime_number(std::string& file_name) {
    std::string reverse_file_name = file_name;
    std::reverse(reverse_file_name.begin(), reverse_file_name.end());
    reverse_file_name.erase(0, 4);
    size_t pos = reverse_file_name.find("_");
    return std::stoi(reverse_file_name.substr(0, pos));
  }

  void Reconstructor::set_safe_interpolation() {
    safe_mode = true;
  }

  void Reconstructor::reconstruct() {
    if (!resume_from_state) {
      if (verbosity > SILENT) {
        INFO_MSG("Promote to new prime field: F(" + std::to_string(primes()[prime_it]) + ").");
      }

      if (safe_mode) {
        tmp_rec.set_safe_interpolation();
      }

      if (scan) {
        scan_for_shift();
        start_probe_jobs(std::vector<uint32_t> (n - 1, 1), thr_n);
        started_probes.emplace(std::vector<uint32_t> (n - 1, 1), thr_n);
      } else {
        start_first_runs();
      }
    }

    run_until_done();

    if (verbosity > SILENT) {
      INFO_MSG("Reconstructed all functions successfully.");
      INFO_MSG(std::to_string(total_iterations) + " probes in total.");
      INFO_MSG("Average time of the black-box probes: " + std::to_string(average_black_box_time) + " s.");
    }

    tp.kill_all();
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
    std::unique_lock<std::mutex> lock_clean(clean);

    std::vector<std::pair<std::string, RationalFunction>> result;

    for (auto & rec : reconst) {
      std::unique_lock<std::mutex> lock_exists(*std::get<1>(rec));

      if (std::get<2>(rec) == DONE) {
        if (tags.size() > 0) {
          result.emplace_back(std::make_pair(tags[std::get<0>(rec)], std::get<3>(rec)->get_result()));
        } else {
          result.emplace_back(std::make_pair(std::to_string(std::get<0>(rec)), std::get<3>(rec)->get_result()));
        }

        std::get<2>(rec) == DELETED;
        delete std::get<3>(rec);
      }
    }
  }

  void Reconstructor::scan_for_shift() {
    if (verbosity > SILENT) {
      INFO_MSG("Scanning for a sparse shift.");
    }

    // Generate all possible combinations of shifting variables
    const auto shift_vec = generate_possible_shifts(n);

    bool first = true;
    bool found_shift = false;
    uint32_t counter = 0;
    uint32_t bound = shift_vec.size();

    tmp_rec.scan_for_sparsest_shift();

    start_first_runs();

    // Run this loop until a proper shift is found
    while (!found_shift && counter != bound) {
      if (!first) {
        tmp_rec.set_zi_shift(shift_vec[counter]);
        shift = tmp_rec.get_zi_shift_vec();
        start_probe_jobs(std::vector<uint32_t> (n - 1, 1), thr_n);
        started_probes.emplace(std::vector<uint32_t> (n - 1, 1), thr_n);
      }

      run_until_done();

      // Kill all jobs
      // otherwise it can happen that a RatReconst is fed with old data
      tp.kill_all();

      found_shift = true;

      for (auto & rec : reconst) {
        std::get<2>(rec) = DEFAULT;

        if (!std::get<3>(rec)->is_shift_working()) {
          found_shift = false;
        }
      }

      if (first) {
        found_shift = false;
        first = false;
      } else {
        ++counter;
      }

      probes.clear();
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
      std::get<2>(rec) = DEFAULT;
      std::get<3>(rec)->accept_shift();
    }

    scan = false;

    if (verbosity > SILENT) {
      if (found_shift) {
        std::string msg = "";

        for (const auto & el : shift_vec[counter - 1]) {
          msg += std::to_string(el) + ", ";
        }

        msg = msg.substr(0, msg.size() - 2);
        INFO_MSG("Shift scan completed successfully. Shift the variable tuple: (" + msg + ").");
      } else {
        INFO_MSG("Shift scan found no sparse shift.");
      }

      INFO_MSG("Total black-box probes for scan: " + std::to_string(total_iterations) + ".");
      INFO_MSG("Average time of the black-box probes: " + std::to_string(average_black_box_time) + " s.");
    }
  }

  void Reconstructor::start_first_runs() {
    shift = tmp_rec.get_zi_shift_vec();
    std::vector<uint32_t> zi_order(n - 1, 1);

    start_probe_jobs(zi_order, thr_n);
    started_probes.emplace(zi_order, thr_n);

    FFInt t = 1;
    std::vector<FFInt> probe {};

    {
      std::unique_lock<std::mutex> lock_future(future_control);

      // sometimes the future is not ready even though the solution job returned already
      // probably there is an additional copy operation involved
      while (probe.empty()) {
        while (jobs_finished == 0) {
          condition_future.wait(lock_future);
        }

        for (auto it = probes.begin(); it != probes.end(); ++it) {
          if ((std::get<2>(*it)).wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
            t = std::get<0>(*it);
            zi_order = std::get<1>(*it);
            std::pair<std::vector<FFInt>, double> tmp = (std::get<2>(*it)).get();
            probe = std::move(tmp.first);
            average_black_box_time = std::move(tmp.second);
            it = probes.erase(it);
            --jobs_finished;
            break;
          }
        }
      }
    }

    ++iteration;

    if (verbosity > SILENT) {
      INFO_MSG("Average time of the black-box probes: " + std::to_string(average_black_box_time) + " s.");
    }

    ++fed_ones;
    items = probe.size();
    size_t tag_size = tags.size();

    for (uint32_t i = 0; i != items; ++i) {
      RatReconst* rec = new RatReconst(n);

      if (safe_mode) {
        rec->set_safe_interpolation();
      }

      if (scan) {
        rec->scan_for_sparsest_shift();
      }

      if (save_states) {
        if (tag_size > 0)
          rec->set_tag(tags[i]);
        else
          rec->set_tag(std::to_string(i));
      }

      if (resume_from_state) {
        rec->start_from_saved_file(file_paths[i]);
      }

      rec->feed(t, probe[i], zi_order, prime_it);
      rec->interpolate();

      std::mutex* mut = new std::mutex;

      reconst.emplace_back(std::make_tuple(i, mut, DEFAULT, rec));
    }

    probe.clear();

    if (verbosity == CHATTY) {
      INFO_MSG("Probe: 1 | Done: 0 / " + std::to_string(items) + " | Needs new prime field: 0 / " + std::to_string(items));
    }

    start_probe_jobs(zi_order, 1);
    ++started_probes.at(zi_order);
  }

  void Reconstructor::run_until_done() {
    FFInt t = 1;
    std::vector<uint32_t> zi_order(n - 1, 1);

    bool done = false;
    bool new_prime = false;

    if (resume_from_state) {
      new_prime = true;
      resume_from_state = false;
    }

    while (!done) {
      if (new_prime) {
        tp.kill_all();

        clean_reconst();

        total_iterations += iteration;
        ++prime_it;

        if (verbosity > SILENT) {
          INFO_MSG("Probes for previous prime field: " + std::to_string(iteration) + ". | " + std::to_string(total_iterations) + " probes in total.");
          INFO_MSG("Average time of the black-box probes: " + std::to_string(average_black_box_time) + " s.");
          INFO_MSG("Reconstructed functions: " + std::to_string(items_done) + " / " + std::to_string(items) + ".");
          INFO_MSG("Promote to new prime field: F(" + std::to_string(primes()[prime_it]) + ").");
        }

        iteration = 0;

        probes.clear();
        fed_ones = 0;
        jobs_finished = 0;
        started_probes.clear();
        feed_jobs = 0;
        interpolate_jobs = 0;
        new_prime = false;
        items_new_prime = 0;

        FFInt::set_new_prime(primes()[prime_it]);

        // if only a small constant is reconstructed it will not ask for new run
        if (probes_for_next_prime == 0) {
          probes_for_next_prime = thr_n;
        }

        if (!tmp_rec.need_shift()) {
          if (tmp_rec.get_zi_shift_vec() != std::vector<FFInt> (n, 0)) {
            if (verbosity > SILENT)
              INFO_MSG("Disable shift.");

            tmp_rec.disable_shift();
          }
        }

        shift = tmp_rec.get_zi_shift_vec();
        tmp_rec.generate_anchor_points();

        // start only thr_n jobs first, because the reconstruction can be done after the first feed
        if (probes_for_next_prime > thr_n) {
          if (verbosity == CHATTY) {
            INFO_MSG("Starting " + std::to_string(thr_n) + " jobs now, the remaining " + std::to_string(probes_for_next_prime - thr_n) + " jobs will be started later.");
          }

          start_probe_jobs(std::vector<uint32_t>(n - 1, 1), thr_n);
          started_probes.emplace(std::vector<uint32_t>(n - 1, 1), thr_n);
        } else {
          if (verbosity == CHATTY) {
            INFO_MSG("Starting " + std::to_string(probes_for_next_prime) + " jobs.");
          }

          start_probe_jobs(std::vector<uint32_t>(n - 1, 1), probes_for_next_prime);
          started_probes.emplace(std::vector<uint32_t>(n - 1, 1), probes_for_next_prime);
        }

        probes_for_next_prime = 0;
      }

      std::vector<FFInt>* probe = new std::vector<FFInt>;
      double time;

      {
        std::unique_lock<std::mutex> lock_future(future_control);

        // sometimes the future is not ready even though the solution job returned already
        // probably there is an additional copy operation involved
        while (probe->empty()) {
          while (jobs_finished == 0) {
            condition_future.wait(lock_future);
          }

          for (auto it = probes.begin(); it != probes.end(); ++it) {
            if ((std::get<2>(*it)).wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
              t = std::get<0>(*it);
              zi_order = std::get<1>(*it);
              std::pair<std::vector<FFInt>, double> tmp = (std::get<2>(*it)).get();
              *probe = std::move(tmp.first);
              time = std::move(tmp.second);
              it = probes.erase(it);
              --jobs_finished;
              break;
            }
          }
        }
      }

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

        if (jobs_finished == 0 && probes.size() == 0) {
          lock_future.unlock();

          bool cont = false;

          {
            std::unique_lock<std::mutex> lock_feed(feed_control);

            while (feed_jobs > 0 || interpolate_jobs > 0) {
              condition_feed.wait(lock_feed);

              lock_feed.unlock();
              lock_future.lock();

              if (jobs_finished > 0 || probes.size() > 0) {
                lock_future.unlock();
                cont = true;
                break;
              }

              lock_future.unlock();
              lock_feed.lock();
            }
          }

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

          if (cont) {
            continue;
          }

          // no jobs are running anymore, check if done or new_prime else throw error
          uint32_t items_new_prime_tmp = 0;

          for (auto & rec : reconst) {
            std::unique_lock<std::mutex> lock_exists(*std::get<1>(rec));

            if (std::get<2>(rec) == DEFAULT) {
              if (!std::get<3>(rec)->is_done()) {
                if (std::get<3>(rec)->get_prime() != prime_it) {
                  ++items_new_prime_tmp;
                }
              }
            }
          }

          if (!scan && items_new_prime_tmp > items_new_prime) {
            items_new_prime = items_new_prime_tmp;
          }

          if (items_done == items) {
            done = true;
          } else if (items_done + items_new_prime_tmp == items) {
            new_prime = true;
          } else {
            throw std::runtime_error("No items to feed anymore");
          }
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

    std::vector<FFInt> values(n);

    for (uint32_t j = 0; j != to_start; ++j) {
      FFInt t = tmp_rec.get_rand();
      values[0] = t + shift[0];

      std::vector<firefly::FFInt> rand_zi = tmp_rec.get_rand_zi_vec(zi_order);

      for (uint32_t i = 1; i != n; ++i) {
        values[i] = rand_zi[i - 1] * t + shift[i];
      }

      std::unique_lock<std::mutex> lock(future_control);

      if (ones) {
        auto future = tp.run_priority_packaged_task([this, values]() {
          auto time0 = std::chrono::high_resolution_clock::now();

//           std::vector<FFInt> probe = (*bb)(values);
          std::vector<FFInt> probe;
          (*bb)(probe, values);

          auto time1 = std::chrono::high_resolution_clock::now();

          std::unique_lock<std::mutex> lock(future_control);

          ++jobs_finished;
          condition_future.notify_one();
          return std::make_pair(std::move(probe), std::chrono::duration<double>(time1 - time0).count());
        });

        probes.emplace_back(std::make_tuple(t, zi_order, std::move(future)));
      } else {
        auto future = tp.run_packaged_task([this, values]() {
          auto time0 = std::chrono::high_resolution_clock::now();

                    std::vector<FFInt> probe;
          (*bb)(probe, values);
          //std::vector<FFInt> probe = (*bb)(values);

          auto time1 = std::chrono::high_resolution_clock::now();

          std::unique_lock<std::mutex> lock(future_control);

          ++jobs_finished;
          condition_future.notify_one();
          return std::make_pair(std::move(probe), std::chrono::duration<double>(time1 - time0).count());
        });

        probes.emplace_back(std::make_tuple(t, zi_order, std::move(future)));
      }
    }
  }

  void Reconstructor::feed_job(const std::vector<uint32_t> zi_order, const firefly::FFInt t, std::vector<FFInt>* probe) {
    {
      std::unique_lock<std::mutex> lock(feed_control);

      interpolate_jobs += items;
    }

    uint32_t items_new_prime_tmp = 0;
    uint32_t counter = 0;

    for (auto & rec : reconst) {
      std::unique_lock<std::mutex> lock_exists(*std::get<1>(rec));

      if (std::get<2>(rec) == DEFAULT) {
        lock_exists.unlock();

        if (!std::get<3>(rec)->is_done()) {
          if (std::get<3>(rec)->get_prime() == prime_it) {
            ++counter;

            std::get<3>(rec)->feed(t, (*probe)[std::get<0>(rec)], zi_order, prime_it);

            tp.run_task([this, &rec]() {
              interpolate_job(rec);
            });
          } else {
            ++items_new_prime_tmp;
          }
        }
      }
    }

    delete probe;

    {
      std::unique_lock<std::mutex> lock_status(status_control);

      if (!scan && (one_done || items_new_prime_tmp > items_new_prime)) {
        if (items_new_prime_tmp > items_new_prime) {
          items_new_prime = items_new_prime_tmp;
        }

        one_done = false;

        std::unique_lock<std::mutex> lock_print(print_control);

        if (verbosity == CHATTY) {
          INFO_MSG("Probe: " + std::to_string(iteration) + " | Done: " + std::to_string(items_done) + " / " + std::to_string(items) + " | " + "Needs new prime field: " + std::to_string(items_new_prime) + " / " + std::to_string(items));
        }
      }
    }

    std::unique_lock<std::mutex> lock(feed_control);

    interpolate_jobs -= (items - counter);
    --feed_jobs;
    condition_feed.notify_one();
  }

  void Reconstructor::interpolate_job(RatReconst_tuple& rec) {
    std::unique_lock<std::mutex> lock_exists(*std::get<1>(rec));

    if (std::get<2>(rec) == DEFAULT) {
      lock_exists.unlock();

      if (!std::get<3>(rec)->interpolate()) {
        //lock_exists.lock();

        // start new jobs if required
        if (!std::get<3>(rec)->is_done()) {
          if (std::get<3>(rec)->get_prime() > prime_it) {
            std::unique_lock<std::mutex> lock(job_control);

            if (std::get<3>(rec)->get_num_eqn() > probes_for_next_prime) {
              probes_for_next_prime = std::get<3>(rec)->get_num_eqn();
            }

            //lock_exists.unlock();
          } else {
            std::vector<uint32_t> zi_order = std::get<3>(rec)->get_zi_order();

            if ((prime_it == 0 || safe_mode == true) && zi_order == std::vector<uint32_t>(n - 1, 1)) {
              //lock_exists.unlock();
              std::unique_lock<std::mutex> lock(job_control);

              if (started_probes.at(zi_order) - thr_n <= fed_ones - 1) {
                uint32_t start = fed_ones - started_probes.at(zi_order) + thr_n;
                started_probes.at(zi_order) += start;

                lock.unlock();

                if (verbosity == CHATTY) {
                  std::unique_lock<std::mutex> lock_print(print_control);

                  INFO_MSG("Starting ones: " + std::to_string(start) + ".");
                }

                start_probe_jobs(zi_order, start);
              }
            } else {
              uint32_t required_probes = std::get<3>(rec)->get_num_eqn();

              //lock_exists.unlock();
              std::unique_lock<std::mutex> lock(job_control);

              auto it = started_probes.find(zi_order);

              if (it != started_probes.end()) {
                if (required_probes > started_probes.at(zi_order)) {
                  uint32_t start = required_probes - started_probes.at(zi_order);
                  started_probes.at(zi_order) = required_probes;

                  lock.unlock();

                  if (verbosity == CHATTY) {
                    std::string msg = "Starting zi_order (";

                    for (const auto & ele : zi_order) {
                      msg += std::to_string(ele) + ", ";
                    }

                    msg = msg.substr(0, msg.length() - 2);
                    msg += ") " + std::to_string(start) + " time(s)";

                    std::unique_lock<std::mutex> lock_print(print_control);

                    INFO_MSG(msg + ".");
                  }

                  start_probe_jobs(zi_order, start);
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

                  INFO_MSG(msg + ".");
                }

                start_probe_jobs(zi_order, required_probes);
              }
            }
          }
        } else /*if (std::get<2>(rec) == DEFAULT)*/ { // to be sure that no other thread does the same
          lock_exists.lock();

          std::get<2>(rec) = DONE;

          lock_exists.unlock();
          std::unique_lock<std::mutex> lock_status(status_control);

          ++items_done;
          one_done = true;
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
