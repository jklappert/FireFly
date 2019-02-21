#include "Reconstructor.hpp"

namespace firefly {
  Reconstructor::Reconstructor(uint32_t n_, uint32_t thr_n_, uint32_t verbosity_): n(n_), thr_n(thr_n_), verbosity(verbosity_), tp(thr_n_) {}

  void Reconstructor::scan_for_sparsest_shift() {
    scan = true;
  }

  void Reconstructor::reconstruct() {
    if (verbosity > 0) {
      INFO_MSG("New prime: 1");
    }

    FFInt::set_new_prime(primes()[prime_it]);
    RatReconst tmp(n);

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
          probe = std::move((std::get<2>(*it)).get());
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

    uint32_t items = 0;
    ++fed_ones;

    int count = 0;
    for (const auto & value : probe) {
      reconst.emplace_back(RatReconst(n));
      // save intermediate results
      reconst.back().set_tag(std::to_string(count));
      reconst.back().feed(t, value, zi_order, prime_it);
      reconst.back().interpolate();
      ++items;
      ++count;
    }

    probe.clear();

    if (verbosity == 2) {
      VERBOSE_MSG("| Iteration: 1 | Done: 0 / " + std::to_string(items) + " | Want new prime: 0 / " + std::to_string(items) + " |");
    }

    start_probe_jobs(std::vector<uint32_t>(n - 1, 1), 1);
    ++started_probes.at(std::vector<uint32_t>(n - 1, 1));

    uint32_t iteration = 1;
    uint32_t total_iterations = 0;

    bool done = false;
    bool new_prime = false;

    while (!done) {
      if (new_prime) {
        total_iterations += iteration;
        ++prime_it;

        if (verbosity > 0) {
          INFO_MSG("New prime: " + std::to_string(prime_it + 1));
          INFO_MSG("Iterations for last prime: " + std::to_string(iteration));
          INFO_MSG("Iterations in total: " + std::to_string(total_iterations));
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

        if (!tmp.need_shift()) {
          if (verbosity > 0) {
            INFO_MSG("Disable shift.");
          }
          tmp.disable_shift();
        }

        tmp.generate_anchor_points();

        // start only coreNumber jobs first, because the reconstruction can be done after the first feed
        if (probes_for_next_prime > thr_n) {
          if (verbosity == 2) {
            VERBOSE_MSG("Starting " + std::to_string(thr_n) + " jobs now, the remaining " + std::to_string(probes_for_next_prime - thr_n) + " jobs will be started later.");
          }

          start_probe_jobs(std::vector<uint32_t>(n - 1, 1), thr_n);
          started_probes.emplace(std::vector<uint32_t>(n - 1, 1), thr_n);
        } else {
          if (verbosity == 2) {
            VERBOSE_MSG("Starting " + std::to_string(probes_for_next_prime) + " jobs");
          }

          start_probe_jobs(std::vector<uint32_t>(n - 1, 1), probes_for_next_prime);
          started_probes.emplace(std::vector<uint32_t>(n - 1, 1), probes_for_next_prime);
        }

        probes_for_next_prime = 0;
      }

      {
        std::unique_lock<std::mutex> lock(mut);
        find_loop:
        while (jobs_finished == 0) {
          cond.wait(lock);
        }

        for (auto it = probes.begin(); it != probes.end(); ++it) {
          if ((std::get<2>(*it)).wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
            t = std::get<0>(*it);
            zi_order = std::get<1>(*it);
            probe = std::move((std::get<2>(*it)).get());
            it = probes.erase(it);
            --jobs_finished;
            break;
          }
        }

        // sometimes the future is not ready even though the solution job returned already
        // probably there is an additional copy operation involved
        if (probe.empty()) {
          goto find_loop;
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

      for (uint32_t i = 0; i != reconst.size(); ++i) {
        if (!reconst[i].is_done()) {
          if (reconst[i].get_prime() == prime_it) {
            {
              std::unique_lock<std::mutex> lock(mut);
              ++feeding_jobs;
            }

            reconst[i].feed(t, probe[i], zi_order, prime_it);

            tp.run_task([this, i]() {
              interpolate_job(reconst[i], i);
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

      if (verbosity == 2) {
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

        for (uint32_t i = 0; i != reconst.size(); ++i) {
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

    if (verbosity > 0) {
      INFO_MSG("Done.");
      INFO_MSG("Iterations for last prime: " + std::to_string(iteration));
      INFO_MSG("Primes used: " + std::to_string(prime_it + 1));
      INFO_MSG("Iterations in total: " + std::to_string(total_iterations));
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

  void Reconstructor::start_probe_jobs(const std::vector<uint32_t>& zi_order, const uint32_t start, const uint32_t i) {
    RatReconst tmp(n);
    std::vector<firefly::FFInt> values(n);
    std::vector<firefly::FFInt> shift = tmp.get_zi_shift_vec();

    for (uint32_t j = 0; j != start; ++j) {
      FFInt t = tmp.get_rand();
      values[0] = t + shift[0];

      std::vector<firefly::FFInt> rand_zi;
      if (i == 0) {
        rand_zi = tmp.get_rand_zi_vec(zi_order);
      } else {
        rand_zi = reconst[i].get_rand_zi_vec(zi_order);
      }

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
        return std::move(probe);
      });

      probes.emplace_back(std::make_tuple(t, zi_order, std::move(future)));
    }
  }

  void Reconstructor::interpolate_job(RatReconst& reconst, const uint i) {
    if (!reconst.interpolate()) {
      // start new jobs if required
      std::unique_lock<std::mutex> lock(mut);

      if (!reconst.is_done()) {
        if (reconst.get_prime() > prime_it) {
          if (reconst.get_num_eqn() > probes_for_next_prime) {
            probes_for_next_prime = reconst.get_num_eqn();
          }
        } else {
          std::vector<uint32_t> zi_order = reconst.get_zi_order();

          if (prime_it == 0 && zi_order == std::vector<uint32_t>(n - 1, 1)) {
            if (started_probes.at(zi_order) - thr_n <= fed_ones - 1) {
              uint32_t start = fed_ones - started_probes.at(zi_order) + thr_n;

              if (verbosity == 2) {
                VERBOSE_MSG("Starting ones: " + std::to_string(start));
              }

              started_probes.at(zi_order) += start;
              lock.unlock();
              start_probe_jobs(zi_order, start, i);
            }
          } else {
            uint32_t required_probes = reconst.get_num_eqn();
            auto it = started_probes.find(zi_order);

            if (it != started_probes.end()) {
              if (required_probes > started_probes.at(zi_order)) {
                uint32_t start = required_probes - started_probes.at(zi_order);

                if (verbosity == 2) {
                  std::string msg = "Starting";
                  for (const auto & ele : zi_order) {
                    msg += " " + std::to_string(ele);
                  }
                  msg += " -- " + std::to_string(start);
                  VERBOSE_MSG(msg);
                }

                started_probes.at(zi_order) = required_probes;
                lock.unlock();
                start_probe_jobs(zi_order, start, i);
              }
            } else {
              if (verbosity == 2) {
                std::string msg = "Starting";
                for (const auto & ele : zi_order) {
                  msg += " " + std::to_string(ele);
                }
                msg += " -- " + std::to_string(required_probes);
                VERBOSE_MSG(msg);
              }

              started_probes.emplace(zi_order, required_probes);
              lock.unlock();
              start_probe_jobs(zi_order, required_probes, i);
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
