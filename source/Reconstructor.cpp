#include "Reconstructor.hpp"
#include "utils.hpp"

namespace firefly {
  Reconstructor::Reconstructor(uint32_t n_, uint32_t thr_n_, uint32_t verbosity_): n(n_), thr_n(thr_n_), verbosity(verbosity_), tp(thr_n_) {
    FFInt::set_new_prime(primes()[prime_it]);
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
    parse_prime_number(file_paths[0]);
    tmp_rec.start_from_saved_file(file_paths[0]);
    items = file_paths.size();

    for (uint32_t i = 0; i != items; ++i) {
      reconst.emplace_back(RatReconst(n));
      reconst[i].start_from_saved_file(file_paths[i]);
      probes_for_next_prime = std::max(probes_for_next_prime, reconst[i].get_num_eqn());
    }
  }

  void Reconstructor::parse_prime_number(std::string& file_name) {
    std::string reverse_file_name = file_name;
    std::reverse(reverse_file_name.begin(), reverse_file_name.end());
    reverse_file_name.erase(0, 4);
    size_t pos = reverse_file_name.find("_");
    prime_it = std::stoi(reverse_file_name.substr(0, pos));
  }

  void Reconstructor::reconstruct() {
    if (!resume_from_state) {
      if (verbosity > SILENT) {
        INFO_MSG("New prime: 1.");
      }

      if (scan) {
        scan_for_shift();
        start_probe_jobs(std::vector<uint32_t>(n - 1, 1), thr_n);
        started_probes.emplace(std::vector<uint32_t>(n - 1, 1), thr_n);
      } else {
        start_first_runs();
      }
    }

    run_until_done();

    if (verbosity > SILENT) {
      INFO_MSG("Done.");
      INFO_MSG("Primes used: " + std::to_string(prime_it + 1) + ".");
      INFO_MSG("Iterations in total: " + std::to_string(total_iterations) + ".");
    }

    tp.kill_all();
  }

  std::vector<RationalFunction> Reconstructor::get_result() {
    std::vector<RationalFunction> result {};

    for (uint32_t i = 0; i != reconst.size(); ++i) {
      result.emplace_back(reconst[i].get_result());
    }

    return result;
  }

  void Reconstructor::scan_for_shift() {
    // Generate all possible combinations of shifting variables
    const auto shift_vec = generate_possible_shifts(n);

    bool first = true;
    bool found_shift = false;
    uint32_t counter = 0;
    uint32_t bound = shift_vec.size();

    tmp_rec.scan_for_sparsest_shift();
    start_first_runs();
    ++total_iterations;

    // Run this loop until a proper shift is found
    while (!found_shift && counter != bound) {
      if (!first) {
        tmp_rec.set_zi_shift(shift_vec[counter]);
        shift = tmp_rec.get_zi_shift_vec();
        start_probe_jobs(std::vector<uint32_t>(n - 1, 1), thr_n);
        started_probes.emplace(std::vector<uint32_t>(n - 1, 1), thr_n);
      }

      run_until_done();

      found_shift = true;

      for (auto & rec : reconst) {
        if (!rec.is_shift_working()) {
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
      feeding_jobs = 0;
    }

    if (found_shift) {
      tmp_rec.set_zi_shift(shift_vec[counter - 1]);
    } else {
      tmp_rec.set_zi_shift(std::vector<uint32_t> (n, 1));
    }

    shift = tmp_rec.get_zi_shift_vec();

    for (auto & rec : reconst) {
      rec.accept_shift();
    }

    scan = false;

    if (verbosity > SILENT) {
      if (found_shift) {
        std::string msg = "";

        for (const auto & el : shift_vec[counter - 1]) {
          msg += " " + std::to_string(el);
        }

        INFO_MSG("Shift scan completed successfully:" + msg + ".");
      } else {
        INFO_MSG("Shift scan found no sparse shift.");
      }

      INFO_MSG("Total black box evaluations for scan: " + std::to_string(total_iterations) + ".");
    }
  }

  void Reconstructor::start_first_runs() {
    shift = tmp_rec.get_zi_shift_vec();

    start_probe_jobs(std::vector<uint32_t>(n - 1, 1), thr_n);
    started_probes.emplace(std::vector<uint32_t>(n - 1, 1), thr_n);

    FFInt t = 1;
    std::vector<uint32_t> zi_order(n - 1, 1);
    std::vector<FFInt> probe {};

    {
      std::unique_lock<std::mutex> lock(mut);
    find:

      while (jobs_finished == 0) {
        cond.wait(lock);
      }

      for (auto it = probes.begin(); it != probes.end(); ++it) {
        if ((std::get<2>(*it)).wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
          t = std::get<0>(*it);
          zi_order = std::get<1>(*it);
          probe = (std::get<2>(*it)).get();
          it = probes.erase(it);
          --jobs_finished;
          break;
        }
      }

      // sometimes the future is not ready even though the solution job returned already
      // probably there is an additional copy operation involved
      if (probe.empty()) {
        goto find;
      }
    }

    ++fed_ones;

    items = probe.size();

    size_t tag_size = tags.size();

    for (uint32_t i = 0; i != items; ++i) {
      reconst.emplace_back(RatReconst(n));

      if (scan) {
        reconst[i].scan_for_sparsest_shift();
      }

      if (save_states) {
        if (tag_size > 0)
          reconst[i].set_tag(tags[i]);
        else
          reconst[i].set_tag(std::to_string(i));
      }

      if (resume_from_state) {
        reconst[i].start_from_saved_file(file_paths[i]);
      }

      reconst[i].feed(t, probe[i], zi_order, prime_it);
      reconst[i].interpolate();
    }

    probe.clear();

    if (verbosity == CHATTY) {
      VERBOSE_MSG("| Iteration: 1 | Done: 0 / " + std::to_string(items) + " | Want new prime: 0 / " + std::to_string(items) + " |");
    }

    start_probe_jobs(std::vector<uint32_t>(n - 1, 1), 1);
    ++started_probes.at(std::vector<uint32_t>(n - 1, 1));
  }

  void Reconstructor::run_until_done() {
    FFInt t = 1;
    std::vector<uint32_t> zi_order(n - 1, 1);
    std::vector<FFInt> probe {};

    uint32_t iteration;

    if (scan) {
      iteration = 0;
    } else {
      iteration = 1;
    }

    bool done = false;
    bool new_prime = false;
    if(resume_from_state){
      new_prime = true;
      resume_from_state = false;
    }

    while (!done) {
      if (new_prime) {
        total_iterations += iteration;
        ++prime_it;

        if (verbosity > SILENT) {
          INFO_MSG("Iterations for last prime: " + std::to_string(iteration) + ".");
          INFO_MSG("Iterations in total: " + std::to_string(total_iterations) + ".");
          INFO_MSG("New prime: " + std::to_string(prime_it + 1) + ".");
        }

        iteration = 0;

        tp.kill_all();

        probes.clear();
        jobs_finished = 0;
        started_probes.clear();
        new_prime = false;
        feeding_jobs = 0;

        FFInt::set_new_prime(primes()[prime_it]);

        // if only a small constant is reconstructed it will not ask for new run
        if (probes_for_next_prime == 0) {
          probes_for_next_prime = 1;
        }

        if (!tmp_rec.need_shift()) {
          if (verbosity > SILENT) {
            INFO_MSG("Disable shift.");
          }

          tmp_rec.disable_shift();
        }

        shift = tmp_rec.get_zi_shift_vec();

        tmp_rec.generate_anchor_points();

        // start only thr_n jobs first, because the reconstruction can be done after the first feed
        if (probes_for_next_prime > thr_n) {
          if (verbosity == CHATTY) {
            VERBOSE_MSG("Starting " + std::to_string(thr_n) + " jobs now, the remaining " + std::to_string(probes_for_next_prime - thr_n) + " jobs will be started later.");
          }

          start_probe_jobs(std::vector<uint32_t>(n - 1, 1), thr_n);
          started_probes.emplace(std::vector<uint32_t>(n - 1, 1), thr_n);
        } else {
          if (verbosity == CHATTY) {
            VERBOSE_MSG("Starting " + std::to_string(probes_for_next_prime) + " jobs.");
          }

          start_probe_jobs(std::vector<uint32_t>(n - 1, 1), probes_for_next_prime);
          started_probes.emplace(std::vector<uint32_t>(n - 1, 1), probes_for_next_prime);
        }

        probes_for_next_prime = 0;
      }

      {
        std::unique_lock<std::mutex> lock(mut);
      find:

        while (jobs_finished == 0) {
          cond.wait(lock);
        }

        for (auto it = probes.begin(); it != probes.end(); ++it) {
          if ((std::get<2>(*it)).wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
            t = std::get<0>(*it);
            zi_order = std::get<1>(*it);
            probe = (std::get<2>(*it)).get();
            it = probes.erase(it);
            --jobs_finished;
            break;
          }
        }

        // sometimes the future is not ready even though the solution job returned already
        // probably there is an additional copy operation involved
        if (probe.empty()) {
          goto find;
        }
      }

      ++iteration;

      uint32_t items_done = 0;
      uint32_t items_new_prime = 0;
      {
        std::unique_lock<std::mutex> lock(mut);

        if (prime_it == 0 && zi_order == std::vector<uint32_t>(n - 1, 1)) {
          ++fed_ones;
        }
      }

      for (uint32_t i = 0; i != items; ++i) {
        if (!reconst[i].is_done()) {
          if (reconst[i].get_prime() == prime_it) {
            {
              std::unique_lock<std::mutex> lock(mut);
              ++feeding_jobs;
            }

            reconst[i].feed(t, probe[i], zi_order, prime_it);

            tp.run_task([this, i]() {
              interpolate_job(reconst[i]);
            });
          } else {
            ++items_new_prime;
          }
        } else {
          ++items_done;
        }
      }

      probe.clear();

      std::unique_lock<std::mutex> lock(mut);

      if (verbosity == CHATTY) {
        VERBOSE_MSG("| Iteration: " + std::to_string(iteration) + " | Done: " + std::to_string(items_done) + " / " + std::to_string(items)
                    + " | " + "Want new prime: " + std::to_string(items_new_prime) + " / " + std::to_string(items) + " |");
      }

      if (items_done == items) {
        done = true;
      } else if (items_done + items_new_prime == items) {
        new_prime = true;
      } else if (probes.size() == 0) {
        bool cont = false;

        while (feeding_jobs > 0) {
          cond.wait(lock);

          if (probes.size() > 0) {
            cont = true;
            break;
          }
        }

        if (cont) {
          continue;
        }

        // no jobs are running anymore, check if done or new_prime else throw error
        items_done = 0;
        items_new_prime = 0;

        for (uint32_t i = 0; i != items; ++i) {
          if (!reconst[i].is_done()) {
            if (reconst[i].get_prime() != prime_it) {
              ++items_new_prime;
            }
          } else {
            ++items_done;
          }
        }

        if (items_done == items) {
          done = true;
        } else if (items_done + items_new_prime == items) {
          new_prime = true;
        } else {
          throw std::runtime_error("No items to feed anymore");
        }
      }
    }

    total_iterations += iteration;

    if (verbosity > SILENT && !scan) {
      INFO_MSG("Iterations for last prime: " + std::to_string(iteration) + ".");
    }
  }

  void Reconstructor::start_probe_jobs(const std::vector<uint32_t>& zi_order, const uint32_t start) {
    std::vector<FFInt> values(n);

    for (uint32_t j = 0; j != start; ++j) {
      FFInt t = tmp_rec.get_rand();
      values[0] = t + shift[0];

      std::vector<firefly::FFInt> rand_zi = tmp_rec.get_rand_zi_vec(zi_order);

      for (uint32_t i = 1; i != n; ++i) {
        values[i] = rand_zi[i - 1] * t + shift[i];
      }

      std::unique_lock<std::mutex> lock(mut);

      auto future = tp.run_packaged_task([this, values]() {
        std::vector<FFInt> probe {};
        black_box(probe, values);

        std::unique_lock<std::mutex> lock(mut);
        ++jobs_finished;
        cond.notify_one();
        return probe;
      });

      probes.emplace_back(std::make_tuple(t, zi_order, std::move(future)));
    }
  }

  void Reconstructor::interpolate_job(RatReconst& rec) {
    if (!rec.interpolate()) {
      // start new jobs if required
      std::unique_lock<std::mutex> lock(mut);

      if (!rec.is_done()) {
        if (rec.get_prime() > prime_it) {
          if (rec.get_num_eqn() > probes_for_next_prime) {
            probes_for_next_prime = rec.get_num_eqn();
          }
        } else {
          std::vector<uint32_t> zi_order = rec.get_zi_order();

          if (prime_it == 0 && zi_order == std::vector<uint32_t>(n - 1, 1)) {
            if (started_probes.at(zi_order) - thr_n <= fed_ones - 1) {
              uint32_t start = fed_ones - started_probes.at(zi_order) + thr_n;

              if (verbosity == CHATTY) {
                VERBOSE_MSG("Starting ones: " + std::to_string(start) + ".");
              }

              started_probes.at(zi_order) += start;
              lock.unlock();
              start_probe_jobs(zi_order, start);
            }
          } else {
            uint32_t required_probes = rec.get_num_eqn();
            auto it = started_probes.find(zi_order);

            if (it != started_probes.end()) {
              if (required_probes > started_probes.at(zi_order)) {
                uint32_t start = required_probes - started_probes.at(zi_order);

                if (verbosity == CHATTY) {
                  std::string msg = "Starting";

                  for (const auto & ele : zi_order) {
                    msg += " " + std::to_string(ele);
                  }

                  msg += " -- " + std::to_string(start);
                  VERBOSE_MSG(msg + ".");
                }

                started_probes.at(zi_order) = required_probes;
                lock.unlock();
                start_probe_jobs(zi_order, start);
              }
            } else {
              if (verbosity == CHATTY) {
                std::string msg = "Starting";

                for (const auto & ele : zi_order) {
                  msg += " " + std::to_string(ele);
                }

                msg += " -- " + std::to_string(required_probes);
                VERBOSE_MSG(msg + ".");
              }

              started_probes.emplace(zi_order, required_probes);
              lock.unlock();
              start_probe_jobs(zi_order, required_probes);
            }
          }
        }
      }
    }

    std::unique_lock<std::mutex> lock(mut);
    --feeding_jobs;
    cond.notify_one();
  }
}
