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

#include "RatReconst.hpp"

namespace firefly {
  RatReconst::RatReconst(const RatReconst& other) : BaseReconst(other) {
    std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
    std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
    std::lock(lock_my_status, lock_other_status);

    first_run = other.first_run;
    coef_mat = other.coef_mat;
    curr_zi = other.curr_zi;
    saved_ti = other.saved_ti;
    t_interpolator = other.t_interpolator;
    coef_n = other.coef_n;
    coef_d = other.coef_d;
    non_solved_degs_den = other.non_solved_degs_den;
    non_solved_degs_num = other.non_solved_degs_num;
    saved_num_num = other.saved_num_num;
    saved_num_den = other.saved_num_den;
    max_deg_num = other.max_deg_num;
    max_deg_den = other.max_deg_den;
    curr_deg_num = other.curr_deg_num;
    curr_deg_den = other.curr_deg_den;
    curr_zi_order_num = other.curr_zi_order_num;
    curr_zi_order_den = other.curr_zi_order_den;
    tmp_solved_coefs_num = other.tmp_solved_coefs_num;
    tmp_solved_coefs_den = other.tmp_solved_coefs_den;
    result = other.result;
    g_ni = other.g_ni;
    g_di = other.g_di;
    combined_ni = other.combined_ni;
    combined_di = other.combined_di;
    coef_mat_num = other.coef_mat_num;
    coef_mat_den = other.coef_mat_den;
    solved_num = other.solved_num;
    solved_den = other.solved_den;
    solved_degs_num = other.solved_degs_num;
    solved_degs_den = other.solved_degs_den;
    normalizer_deg = other.normalizer_deg;
    normalizer_den_num = other.normalizer_den_num;
    is_singular_system = other.is_singular_system;
    queue = other.queue;
    const_den = other.const_den;
    tag = other.tag;
    sub_num = other.sub_num;
    sub_den = other.sub_den;
    parsed_variables = other.parsed_variables;
    curr_parsed_variable = other.curr_parsed_variable;
    scan = other.scan;
    all_shift_max_degs = other.all_shift_max_degs;
    shift_works = other.shift_works;
    normalize_to_den = other.normalize_to_den;
    start_deg_num = other.start_deg_num;
    start_deg_den = other.start_deg_den;
    shifted_max_num_eqn = other.shifted_max_num_eqn;
    shifted_degs_num = other.shifted_degs_num;
    shifted_degs_den = other.shifted_degs_den;
    zero_degs_num = other.zero_degs_num;
    zero_degs_den = other.zero_degs_den;
    num_sub_den = other.num_sub_den;
    num_sub_num = other.num_sub_num;
    interpolations = other.interpolations;
    t_interpolator = other.t_interpolator;
    div_by_zero = other.div_by_zero;
    first_feed = other.first_feed;
    zero_counter = other.zero_counter;
    check_interpolation = other.check_interpolation;

    done = other.done;
    new_prime = other.new_prime;
    check = other.check;
    use_chinese_remainder = other.use_chinese_remainder;
    curr_zi_order = other.curr_zi_order;
    prime_number = other.prime_number;
    num_eqn = other.num_eqn;
    n = other.n;
    type = other.type;
    zi = other.zi;
    combined_prime = other.combined_prime;
  }

  RatReconst::RatReconst(RatReconst && other) : BaseReconst(other) {
    std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
    std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
    std::lock(lock_my_status, lock_other_status);

    first_run = std::move(other.first_run);
    coef_mat = std::move(other.coef_mat);
    curr_zi = std::move(other.curr_zi);
    saved_ti = std::move(other.saved_ti);
    t_interpolator = std::move(other.t_interpolator);
    coef_n = std::move(other.coef_n);
    coef_d = std::move(other.coef_d);
    non_solved_degs_den = std::move(other.non_solved_degs_den);
    non_solved_degs_num = std::move(other.non_solved_degs_num);
    saved_num_num = std::move(other.saved_num_num);
    saved_num_den = std::move(other.saved_num_den);
    max_deg_num = std::move(other.max_deg_num);
    max_deg_den = std::move(other.max_deg_den);
    curr_deg_num = std::move(other.curr_deg_num);
    curr_deg_den = std::move(other.curr_deg_den);
    curr_zi_order_num = std::move(other.curr_zi_order_num);
    curr_zi_order_den = std::move(other.curr_zi_order_den);
    tmp_solved_coefs_num = std::move(other.tmp_solved_coefs_num);
    tmp_solved_coefs_den = std::move(other.tmp_solved_coefs_den);
    result = std::move(other.result);
    g_ni = std::move(other.g_ni);
    g_di = std::move(other.g_di);
    combined_ni = std::move(other.combined_ni);
    combined_di = std::move(other.combined_di);
    coef_mat_num = std::move(other.coef_mat_num);
    coef_mat_den = std::move(other.coef_mat_den);
    solved_num = std::move(other.solved_num);
    solved_den = std::move(other.solved_den);
    solved_degs_num = std::move(other.solved_degs_num);
    solved_degs_den = std::move(other.solved_degs_den);
    normalizer_deg = std::move(other.normalizer_deg);
    normalizer_den_num = std::move(other.normalizer_den_num);
    is_singular_system = std::move(other.is_singular_system);
    queue = std::move(other.queue);
    const_den = std::move(other.const_den);
    tag = std::move(other.tag);
    sub_num = std::move(other.sub_num);
    sub_den = std::move(other.sub_den);
    parsed_variables = std::move(other.parsed_variables);
    curr_parsed_variable = std::move(other.curr_parsed_variable);
    scan = std::move(other.scan);
    all_shift_max_degs = std::move(other.all_shift_max_degs);
    shift_works = std::move(other.shift_works);
    normalize_to_den = std::move(other.normalize_to_den);
    start_deg_num = std::move(other.start_deg_num);
    start_deg_den = std::move(other.start_deg_den);
    shifted_max_num_eqn = std::move(other.shifted_max_num_eqn);
    shifted_degs_num = std::move(other.shifted_degs_num);
    shifted_degs_den = std::move(other.shifted_degs_den);
    zero_degs_num = std::move(other.zero_degs_num);
    zero_degs_den = std::move(other.zero_degs_den);
    num_sub_den = std::move(other.num_sub_den);
    num_sub_num = std::move(other.num_sub_num);
    interpolations = std::move(other.interpolations);
    t_interpolator = std::move(other.t_interpolator);
    div_by_zero = std::move(other.div_by_zero);
    first_feed = std::move(other.first_feed);
    zero_counter = std::move(other.zero_counter);
    check_interpolation = std::move(other.check_interpolation);

    done = std::move(other.done);
    new_prime = std::move(other.new_prime);
    check = std::move(other.check);
    use_chinese_remainder = std::move(other.use_chinese_remainder);
    curr_zi_order = std::move(other.curr_zi_order);
    prime_number = std::move(other.prime_number);
    num_eqn = std::move(other.num_eqn);
    n = std::move(other.n);
    type = std::move(other.type);
    zi = std::move(other.zi);
    combined_prime = std::move(other.combined_prime);
  }

  RatReconst& RatReconst::operator=(const RatReconst& other) {
    if (this != &other) {
      std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
      std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
      std::lock(lock_my_status, lock_other_status);

      first_run = other.first_run;
      coef_mat = other.coef_mat;
      curr_zi = other.curr_zi;
      saved_ti = other.saved_ti;
      t_interpolator = other.t_interpolator;
      coef_n = other.coef_n;
      coef_d = other.coef_d;
      saved_num_num = other.saved_num_num;
      saved_num_den = other.saved_num_den;
      non_solved_degs_den = other.non_solved_degs_den;
      non_solved_degs_num = other.non_solved_degs_num;
      max_deg_num = other.max_deg_num;
      max_deg_den = other.max_deg_den;
      curr_deg_num = other.curr_deg_num;
      curr_deg_den = other.curr_deg_den;
      curr_zi_order_num = other.curr_zi_order_num;
      curr_zi_order_den = other.curr_zi_order_den;
      tmp_solved_coefs_num = other.tmp_solved_coefs_num;
      tmp_solved_coefs_den = other.tmp_solved_coefs_den;
      result = other.result;
      g_ni = other.g_ni;
      g_di = other.g_di;
      combined_ni = other.combined_ni;
      combined_di = other.combined_di;
      coef_mat_num = other.coef_mat_num;
      coef_mat_den = other.coef_mat_den;
      solved_num = other.solved_num;
      solved_den = other.solved_den;
      solved_degs_num = other.solved_degs_num;
      solved_degs_den = other.solved_degs_den;
      normalizer_deg = other.normalizer_deg;
      normalizer_den_num = other.normalizer_den_num;
      is_singular_system = other.is_singular_system;
      queue = other.queue;
      const_den = other.const_den;
      tag = other.tag;
      sub_num = other.sub_num;
      sub_den = other.sub_den;
      parsed_variables = other.parsed_variables;
      curr_parsed_variable = other.curr_parsed_variable;
      scan = other.scan;
      all_shift_max_degs = other.all_shift_max_degs;
      shift_works = other.shift_works;
      normalize_to_den = other.normalize_to_den;
      start_deg_num = other.start_deg_num;
      start_deg_den = other.start_deg_den;
      shifted_max_num_eqn = other.shifted_max_num_eqn;
      shifted_degs_num = other.shifted_degs_num;
      shifted_degs_den = other.shifted_degs_den;
      zero_degs_num = other.zero_degs_num;
      zero_degs_den = other.zero_degs_den;
      num_sub_den = other.num_sub_den;
      num_sub_num = other.num_sub_num;
      interpolations = other.interpolations;
      t_interpolator = other.t_interpolator;
      div_by_zero = other.div_by_zero;
      first_feed = other.first_feed;
      zero_counter = other.zero_counter;
      check_interpolation = other.check_interpolation;

      done = other.done;
      new_prime = other.new_prime;
      check = other.check;
      use_chinese_remainder = other.use_chinese_remainder;
      curr_zi_order = other.curr_zi_order;
      prime_number = other.prime_number;
      num_eqn = other.num_eqn;
      n = other.n;
      type = other.type;
      zi = other.zi;
      combined_prime = other.combined_prime;
    }

    return *this;
  }

  RatReconst& RatReconst::operator=(RatReconst && other) {
    if (this != &other) {
      std::unique_lock<std::mutex> lock_my_status(mutex_status, std::defer_lock);
      std::unique_lock<std::mutex> lock_other_status(other.mutex_status, std::defer_lock);
      std::lock(lock_my_status, lock_other_status);

      first_run = std::move(other.first_run);
      coef_mat = std::move(other.coef_mat);
      curr_zi = std::move(other.curr_zi);
      saved_ti = std::move(other.saved_ti);
      t_interpolator = std::move(other.t_interpolator);
      coef_n = std::move(other.coef_n);
      coef_d = std::move(other.coef_d);
      saved_num_num = std::move(other.saved_num_num);
      saved_num_den = std::move(other.saved_num_den);
      non_solved_degs_den = std::move(other.non_solved_degs_den);
      non_solved_degs_num = std::move(other.non_solved_degs_num);
      max_deg_num = std::move(other.max_deg_num);
      max_deg_den = std::move(other.max_deg_den);
      curr_deg_num = std::move(other.curr_deg_num);
      curr_deg_den = std::move(other.curr_deg_den);
      curr_zi_order_num = std::move(other.curr_zi_order_num);
      curr_zi_order_den = std::move(other.curr_zi_order_den);
      tmp_solved_coefs_num = std::move(other.tmp_solved_coefs_num);
      tmp_solved_coefs_den = std::move(other.tmp_solved_coefs_den);
      result = std::move(other.result);
      g_ni = std::move(other.g_ni);
      g_di = std::move(other.g_di);
      combined_ni = std::move(other.combined_ni);
      combined_di = std::move(other.combined_di);
      coef_mat_num = std::move(other.coef_mat_num);
      coef_mat_den = std::move(other.coef_mat_den);
      solved_num = std::move(other.solved_num);
      solved_den = std::move(other.solved_den);
      solved_degs_num = std::move(other.solved_degs_num);
      solved_degs_den = std::move(other.solved_degs_den);
      normalizer_deg = std::move(other.normalizer_deg);
      normalizer_den_num = std::move(other.normalizer_den_num);
      is_singular_system = std::move(other.is_singular_system);
      queue = std::move(other.queue);
      const_den = std::move(other.const_den);
      tag = std::move(other.tag);
      sub_num = std::move(other.sub_num);
      sub_den = std::move(other.sub_den);
      parsed_variables = std::move(other.parsed_variables);
      curr_parsed_variable = std::move(other.curr_parsed_variable);
      scan = std::move(other.scan);
      all_shift_max_degs = std::move(other.all_shift_max_degs);
      shift_works = std::move(other.shift_works);
      normalize_to_den = std::move(other.normalize_to_den);
      start_deg_num = std::move(other.start_deg_num);
      start_deg_den = std::move(other.start_deg_den);
      shifted_max_num_eqn = std::move(other.shifted_max_num_eqn);
      shifted_degs_num = std::move(other.shifted_degs_num);
      shifted_degs_den = std::move(other.shifted_degs_den);
      zero_degs_num = std::move(other.zero_degs_num);
      zero_degs_den = std::move(other.zero_degs_den);
      num_sub_den = std::move(other.num_sub_den);
      num_sub_num = std::move(other.num_sub_num);
      interpolations = std::move(other.interpolations);
      t_interpolator = std::move(other.t_interpolator);
      div_by_zero = std::move(other.div_by_zero);
      first_feed = std::move(other.first_feed);
      zero_counter = std::move(other.zero_counter);
      check_interpolation = std::move(other.check_interpolation);

      done = std::move(other.done);
      new_prime = std::move(other.new_prime);
      check = std::move(other.check);
      use_chinese_remainder = std::move(other.use_chinese_remainder);
      curr_zi_order = std::move(other.curr_zi_order);
      prime_number = std::move(other.prime_number);
      num_eqn = std::move(other.num_eqn);
      n = std::move(other.n);
      type = std::move(other.type);
      zi = std::move(other.zi);
      combined_prime = std::move(other.combined_prime);
    }

    return *this;
  }

}
