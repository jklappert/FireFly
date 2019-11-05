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

#include <vector>

namespace firefly {
  /**
   *  example for n = 4 which uses the singular_solver
   */
  template<typename FFIntTemp>
  FFIntTemp singular_solver(const std::vector<FFIntTemp>& yis) {
    FFIntTemp cr_1(mpz_class("123456789109898799879870980"));
    FFIntTemp test = 17 * yis[0];
    FFIntTemp num = 17 * yis[0] + 7 * yis[1] + yis[2] + yis[0].pow(2) + yis[1].pow(2) + yis[2].pow(2) + yis[3].pow(2) + yis[0] * yis[3].pow(3) + yis[1].pow(4) + yis[3].pow(7) + yis[2].pow(7);
    FFIntTemp den = cr_1 * yis[1] - yis[3] + yis[0] * yis[1] + yis[1] * yis[2] + yis[0] * yis[3] + yis[0].pow(2) * yis[1].pow(2) + yis[2].pow(4);
    return num / den;
  }
  /**
   *  example for n = 1
   */
  template<typename FFIntTemp>
  FFIntTemp n_eq_1(const FFIntTemp& z1) {
    FFIntTemp num = (576 * z1.pow(12) - 35145 * z1.pow(11)
                 + 946716 * z1.pow(10) - 14842335 * z1.pow(9)
                 + 150236238 * z1.pow(8) - 1028892363 * z1.pow(7)
                 + 4853217576 * z1.pow(6) - 15724949577 * z1.pow(5)
                 + 34208917206 * z1.pow(4) - 47506433412 * z1.pow(3)
                 + 37933483608 * z1.pow(2) - 13296184128 * z1 + 71850240);
    FFIntTemp den = (16 * z1.pow(12) - 960 * z1.pow(11)
                 + 25456 * z1.pow(10) - 393440 * z1.pow(9)
                 + 3934768 * z1.pow(8) - 26714240 * z1.pow(7)
                 + 125545488 * z1.pow(6) - 408157280 * z1.pow(5)
                 + 899198016 * z1.pow(4) - 1278172800 * z1.pow(3)
                 + 1055033856 * z1.pow(2) - 383201280 * z1);
    return num / den;
  }
  /**
   *  example for n = 4 and the usage of the Chinese Remainder Theorem
   */
  template<typename FFIntTemp>
  FFIntTemp n_eq_4(const std::vector<FFIntTemp>& yis) {
    mpz_class cr_1_mpz;
    cr_1_mpz = "123456789109898799879870980";
    mpz_class cr_2_mpz;
    cr_2_mpz = "123456789109898799879";
    FFIntTemp cr_1(cr_1_mpz);
    FFIntTemp cr_2(cr_2_mpz);
    FFIntTemp z1 = yis[0];
    std::vector<FFIntTemp> t_yis(3);
    t_yis[0] = yis[1];
    t_yis[1] = yis[2];
    t_yis[2] = yis[3];
    FFIntTemp den = cr_1 * (((z1.pow(3) - 12 * z1.pow(2) + 48 * z1 - 64) * t_yis[1].pow(2))
                        * t_yis[0].pow(5) + ((-3 * z1.pow(3) + 36 * z1.pow(2)
                                              - 144 * z1 + 192) * t_yis[1].pow(2)) * t_yis[0].pow(4) + ((2 * z1.pow(3) - 24 * z1.pow(2)
                                                  + 96 * z1 - 128) * t_yis[1].pow(2)) * t_yis[0].pow(3) + ((2 * z1.pow(3) - 24 * z1.pow(2)
                                                      + 96 * z1 - 128) * t_yis[1].pow(2)) * t_yis[0].pow(2) + ((-3 * z1.pow(3) + 36 * z1.pow(2)
                                                          - 144 * z1 + 192) * t_yis[1].pow(2)) * t_yis[0] + (z1.pow(3) - 12 * z1.pow(2) + 48 * z1 - 64)
                        * t_yis[1].pow(2));
    FFIntTemp num = cr_2 * ((-6 * z1.pow(3) + 54 * z1.pow(2) - 156 * z1 + 144)
                        * t_yis[0].pow(4) + ((-4 * z1.pow(3) + 36 * z1.pow(2) - 104 * z1 + 96) * t_yis[1]
                                             + 9 * z1.pow(3) - 84 * z1.pow(2) + 252 * z1 - 240) * t_yis[0].pow(3) + ((46 * z1.pow(3)
                                                 - 389 * z1.pow(2) + 1074 * z1 - 960) * t_yis[1] - 3 * z1.pow(3) + 30 * z1.pow(2) - 96 * z1 + 96)
                        * t_yis[0].pow(2) + ((-10 * z1.pow(3) + 93 * z1.pow(2) - 278 * z1 + 264)
                                             * t_yis[1]) * t_yis[0]) + z1.pow(15) * t_yis[0].pow(15) * t_yis[1].pow(15) * t_yis[2].pow(15);
    return num / den;
  }
  /**
   *  example for gg->HH large interpolation problem augmented with large coefficients
   */
  template<typename FFIntTemp>
  FFIntTemp gghh(const std::vector<FFIntTemp>& yis) {
    FFIntTemp z1 = yis[0];
    std::vector<FFIntTemp> t_yis(3);
    t_yis[0] = yis[1];
    t_yis[1] = yis[2];
    t_yis[2] = yis[3];
    mpz_class cr_1_mpz;
    cr_1_mpz = "123456789109898799879870980";
    mpz_class cr_2_mpz;
    cr_2_mpz = "123456789109898799879";
    FFIntTemp cr_1(cr_1_mpz);
    FFIntTemp cr_2(cr_2_mpz);
    FFIntTemp num = cr_1 * (t_yis[0] * (-96 + 96 * z1 - 24 * z1 . pow(2)) * t_yis[2] . pow(9) +
                        t_yis[2] . pow(8) * (432 - 384 * z1 + 84 * z1 . pow(2) +
                                             t_yis[0] * t_yis[1] * (432 - 432 * z1 + 108 * z1 . pow(2)) + t_yis[0] * (1552 - 1500 * z1 + 362 * z1 . pow(2)) +
                                             (16 - 104 * z1 + 48 * z1 . pow(2)) * t_yis[0] . pow(2)) +
                        (1024 - 1024 * z1 + 256 * z1 . pow(2)) * t_yis[0] . pow(6) +
                        (320 - 352 * z1 + 96 * z1 . pow(2)) * t_yis[0] . pow(7) +
                        t_yis[1] * ((-256 + 256 * z1 - 64 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                    (4096 - 4032 * z1 + 992 * z1 . pow(2)) * t_yis[0] . pow(5) +
                                    (1728 - 1848 * z1 + 492 * z1 . pow(2)) * t_yis[0] . pow(6) +
                                    (160 - 196 * z1 + 58 * z1 . pow(2)) * t_yis[0] . pow(7)) +
                        (16 - 24 * z1 + 8 * z1 . pow(2)) * t_yis[0] . pow(8) +
                        ((192 - 96 * z1) * t_yis[0] . pow(3) + (6336 - 6128 * z1 + 1480 * z1 . pow(2)) * t_yis[0] . pow(4) +
                         (3752 - 3908 * z1 + 1016 * z1 . pow(2)) * t_yis[0] . pow(5) +
                         (560 - 620 * z1 + 170 * z1 . pow(2)) * t_yis[0] . pow(6)) * t_yis[1] . pow(2) +
                        t_yis[2] . pow(7) * (-2928 + 2616 * z1 + t_yis[0] * (-7128 + 6824 * z1 - 1630 * z1 . pow(2)) -
                                             576 * z1 . pow(2) + (-4464 + 4882 * z1 - 1325 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                             t_yis[1] * (-1968 + 1748 * z1 + t_yis[0] * (-6864 + 6622 * z1 - 1595 * z1 . pow(2)) -
                                                         382 * z1 . pow(2) + (-352 + 608 * z1 - 216 * z1 . pow(2)) * t_yis[0] . pow(2)) +
                                             (440 - 184 * z1 - 18 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                             t_yis[0] * (-792 + 792 * z1 - 198 * z1 . pow(2)) * t_yis[1] . pow(2)) +
                        ((1152 - 960 * z1 + 192 * z1 . pow(2)) * t_yis[0] . pow(2) +
                         (4672 - 4432 * z1 + 1048 * z1 . pow(2)) * t_yis[0] . pow(3) +
                         (4128 - 4192 * z1 + 1064 * z1 . pow(2)) * t_yis[0] . pow(4) +
                         (960 - 1000 * z1 + 260 * z1 . pow(2)) * t_yis[0] . pow(5)) * t_yis[1] . pow(3) +
                        t_yis[2] . pow(6) * (7488 - 6624 * z1 + 1440 * z1 . pow(2) +
                                             t_yis[0] * (12992 - 12152 * z1 + 2828 * z1 . pow(2)) +
                                             (27364 - 27722 * z1 + 7020 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                             (4372 - 5934 * z1 + 1874 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                             t_yis[1] * (12160 - 10880 * z1 + 2400 * z1 . pow(2) +
                                                         t_yis[0] * (31000 - 29628 * z1 + 7064 * z1 . pow(2)) +
                                                         (17024 - 18046 * z1 + 4767 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                                         (-532 - 40 * z1 + 153 * z1 . pow(2)) * t_yis[0] . pow(3)) +
                                             (-644 + 358 * z1 - 18 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                             (3644 - 3234 * z1 + 706 * z1 . pow(2) + t_yis[0] * (12236 - 11784 * z1 + 2833 * z1 . pow(2)) +
                                              (868 - 1154 * z1 + 360 * z1 . pow(2)) * t_yis[0] . pow(2)) * t_yis[1] . pow(2) +
                                             t_yis[0] * (756 - 756 * z1 + 189 * z1 . pow(2)) * t_yis[1] . pow(3)) +
                        (t_yis[0] * (704 - 608 * z1 + 128 * z1 . pow(2)) + (1600 - 1488 * z1 + 344 * z1 . pow(2)) *
                         t_yis[0] . pow(2) + (2352 - 2328 * z1 + 576 * z1 . pow(2)) * t_yis[0] . pow(3) +
                         (880 - 880 * z1 + 220 * z1 . pow(2)) * t_yis[0] . pow(4)) * t_yis[1] . pow(4) +
                        t_yis[2] . pow(5) * (-8448 + 7296 * z1 + t_yis[0] * (-7152 + 5976 * z1 - 1200 * z1 . pow(2)) -
                                             1536 * z1 . pow(2) + (-55784 + 54016 * z1 - 13062 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                             (-39244 + 42208 * z1 - 11293 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                             (-1584 + 3510 * z1 - 1359 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                             t_yis[1] * (-27264 + 24192 * z1 + t_yis[0] * (-57728 + 54124 * z1 - 12630 * z1 . pow(2)) -
                                                         5280 * z1 . pow(2) + (-89256 + 89926 * z1 - 22649 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                                         (-15604 + 18768 * z1 - 5483 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                                         (708 - 246 * z1 - 54 * z1 . pow(2)) * t_yis[0] . pow(4)) +
                                             (368 - 214 * z1 + 15 * z1 . pow(2)) * t_yis[0] . pow(5) +
                                             (-20192 + 18096 * z1 + t_yis[0] * (-52676 + 50294 * z1 - 11978 * z1 . pow(2)) -
                                              4000 * z1 . pow(2) + (-25224 + 26030 * z1 - 6709 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                              (-84 + 546 * z1 - 252 * z1 . pow(2)) * t_yis[0] . pow(3)) * t_yis[1] . pow(2) +
                                             (-3496 + 3100 * z1 + t_yis[0] * (-11204 + 10772 * z1 - 2585 * z1 . pow(2)) -
                                              676 * z1 . pow(2) + (-820 + 974 * z1 - 282 * z1 . pow(2)) * t_yis[0] . pow(2)) *
                                             t_yis[1] . pow(3) + t_yis[0] * (-396 + 396 * z1 - 99 * z1 . pow(2)) * t_yis[1] . pow(4)) +
                        (t_yis[0] * (192 - 176 * z1 + 40 * z1 . pow(2)) + (608 - 584 * z1 + 140 * z1 . pow(2)) *
                         t_yis[0] . pow(2) + (416 - 404 * z1 + 98 * z1 . pow(2)) * t_yis[0] . pow(3)) * t_yis[1] . pow(5) +
                        t_yis[2] . pow(4) * (t_yis[0] * (-7104 + 6880 * z1 - 1664 * z1 . pow(2)) +
                                             (42016 - 37152 * z1 + 8072 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                             (80976 - 81880 * z1 + 20696 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                             (29232 - 33104 * z1 + 9244 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                             (-200 - 964 * z1 + 532 * z1 . pow(2)) * t_yis[0] . pow(5) +
                                             t_yis[1] * (28160 - 24320 * z1 + 5120 * z1 . pow(2) +
                                                         t_yis[0] * (47008 - 41280 * z1 + 8888 * z1 . pow(2)) +
                                                         (152944 - 150168 * z1 + 36848 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                                         (102164 - 107388 * z1 + 28153 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                                         (6336 - 9320 * z1 + 3076 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                                         (-260 + 100 * z1 + 15 * z1 . pow(2)) * t_yis[0] . pow(5)) +
                                             (-92 + 52 * z1 - 3 * z1 . pow(2)) * t_yis[0] . pow(6) +
                                             (37632 - 33536 * z1 + 7360 * z1 . pow(2) +
                                              t_yis[0] * (91480 - 86148 * z1 + 20204 * z1 . pow(2)) +
                                              (115504 - 115030 * z1 + 28639 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                              (18808 - 21062 * z1 + 5829 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                              (-120 - 120 * z1 + 90 * z1 . pow(2)) * t_yis[0] . pow(4)) * t_yis[1] . pow(2) +
                                             (17056 - 15312 * z1 + 3392 * z1 . pow(2) +
                                              t_yis[0] * (44396 - 42362 * z1 + 10082 * z1 . pow(2)) +
                                              (17808 - 18020 * z1 + 4558 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                              (280 - 440 * z1 + 150 * z1 . pow(2)) * t_yis[0] . pow(3)) * t_yis[1] . pow(3) +
                                             (1824 - 1616 * z1 + 352 * z1 . pow(2) + t_yis[0] * (5536 - 5314 * z1 + 1273 * z1 . pow(2)) +
                                              (340 - 380 * z1 + 105 * z1 . pow(2)) * t_yis[0] . pow(2)) * t_yis[1] . pow(4) +
                                             t_yis[0] * (108 - 108 * z1 + 27 * z1 . pow(2)) * t_yis[1] . pow(5)) +
                        (t_yis[0] * (40 - 36 * z1 + 8 * z1 . pow(2)) + (80 - 76 * z1 + 18 * z1 . pow(2)) * t_yis[0] . pow(2)) *
                        t_yis[1] . pow(6) + t_yis[2] . pow(3) * ((17920 - 17152 * z1 + 4096 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                                                 (-64176 + 59480 * z1 - 13696 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                                                 (-57136 + 59748 * z1 - 15590 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                                                 (-12208 + 14490 * z1 - 4193 * z1 . pow(2)) * t_yis[0] . pow(5) +
                                                                 (340 + 36 * z1 - 103 * z1 . pow(2)) * t_yis[0] . pow(6) +
                                                                 t_yis[1] * ((-4480 + 2240 * z1) * t_yis[0] + (-105696 + 96592 * z1 - 21872 * z1 . pow(2)) *
                                                                             t_yis[0] . pow(2) + (-162520 + 165104 * z1 - 41922 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                                                             (-58104 + 63262 * z1 - 17105 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                                                             (-800 + 2092 * z1 - 846 * z1 . pow(2)) * t_yis[0] . pow(5) + (28 - 8 * z1 - 3 * z1 . pow(2)) *
                                                                             t_yis[0] . pow(6)) + (8 - 4 * z1) * t_yis[0] . pow(7) +
                                                                 (-33792 + 29184 * z1 + t_yis[0] * (-73216 + 65632 * z1 - 14512 * z1 . pow(2)) -
                                                                  6144 * z1 . pow(2) + (-163496 + 160748 * z1 - 39500 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                                                  (-97972 + 100656 * z1 - 25835 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                                                  (-5836 + 7398 * z1 - 2240 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                                                  (20 + 20 * z1 - 15 * z1 . pow(2)) * t_yis[0] . pow(5)) * t_yis[1] . pow(2) +
                                                                 (-24192 + 21696 * z1 + t_yis[0] * (-65808 + 62312 * z1 - 14704 * z1 . pow(2)) -
                                                                  4800 * z1 . pow(2) + (-70952 + 69918 * z1 - 17221 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                                                  (-9308 + 10002 * z1 - 2674 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                                                  (-40 + 80 * z1 - 30 * z1 . pow(2)) * t_yis[0] . pow(4)) * t_yis[1] . pow(3) +
                                                                 (-7696 + 6920 * z1 + t_yis[0] * (-19364 + 18466 * z1 - 4392 * z1 . pow(2)) -
                                                                  1536 * z1 . pow(2) + (-6008 + 5998 * z1 - 1497 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                                                  (-80 + 100 * z1 - 30 * z1 . pow(2)) * t_yis[0] . pow(3)) * t_yis[1] . pow(4) +
                                                                 (-488 + 432 * z1 + t_yis[0] * (-1396 + 1338 * z1 - 320 * z1 . pow(2)) - 94 * z1 . pow(2) +
                                                                  (-52 + 56 * z1 - 15 * z1 . pow(2)) * t_yis[0] . pow(2)) * t_yis[1] . pow(5) +
                                                                 t_yis[0] * (-12 + 12 * z1 - 3 * z1 . pow(2)) * t_yis[1] . pow(6)) +
                        t_yis[2] . pow(2) * ((-8000 + 7840 * z1 - 1920 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                             (39168 - 37520 * z1 + 8968 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                             (21368 - 22876 * z1 + 6096 * z1 . pow(2)) * t_yis[0] . pow(5) +
                                             (2800 - 3512 * z1 + 1056 * z1 . pow(2)) * t_yis[0] . pow(6) +
                                             t_yis[1] * ((-10368 + 10944 * z1 - 2880 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                                         (97664 - 92368 * z1 + 21768 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                                         (80432 - 83712 * z1 + 21748 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                                         (17120 - 19314 * z1 + 5377 * z1 . pow(2)) * t_yis[0] . pow(5) +
                                                         (-140 - 124 * z1 + 97 * z1 . pow(2)) * t_yis[0] . pow(6)) +
                                             (-96 + 34 * z1 + 7 * z1 . pow(2)) * t_yis[0] . pow(7) +
                                             (t_yis[0] * (19712 - 16000 * z1 + 3072 * z1 . pow(2)) +
                                              (92512 - 85936 * z1 + 19840 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                              (115064 - 116396 * z1 + 29432 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                              (38496 - 40630 * z1 + 10691 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                              (540 - 970 * z1 + 350 * z1 . pow(2)) * t_yis[0] . pow(5)) * t_yis[1] . pow(2) +
                                             (16896 - 14592 * z1 + 3072 * z1 . pow(2) +
                                              t_yis[0] * (41120 - 37520 * z1 + 8480 * z1 . pow(2)) +
                                              (76064 - 74744 * z1 + 18356 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                              (40820 - 41172 * z1 + 10381 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                              (1640 - 1960 * z1 + 570 * z1 . pow(2)) * t_yis[0] . pow(4)) * t_yis[1] . pow(3) +
                                             (7104 - 6432 * z1 + 1440 * z1 . pow(2) + t_yis[0] * (21824 - 20768 * z1 + 4928 * z1 . pow(2)) +
                                              (20684 - 20196 * z1 + 4927 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                              (1720 - 1810 * z1 + 475 * z1 . pow(2)) * t_yis[0] . pow(3)) * t_yis[1] . pow(4) +
                                             (1760 - 1584 * z1 + 352 * z1 . pow(2) + t_yis[0] * (4092 - 3898 * z1 + 926 * z1 . pow(2)) +
                                              (804 - 796 * z1 + 197 * z1 . pow(2)) * t_yis[0] . pow(2)) * t_yis[1] . pow(5) +
                                             (52 - 46 * z1 + 10 * z1 . pow(2) + t_yis[0] * (140 - 134 * z1 + 32 * z1 . pow(2))) *
                                             t_yis[1] . pow(6)) + t_yis[2] * ((1024 - 1024 * z1 + 256 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                                 (-10496 + 10304 * z1 - 2528 * z1 . pow(2)) * t_yis[0] . pow(5) +
                                                 (-4080 + 4440 * z1 - 1200 * z1 . pow(2)) * t_yis[0] . pow(6) +
                                                 (-328 + 444 * z1 - 140 * z1 . pow(2)) * t_yis[0] . pow(7) +
                                                 t_yis[1] * ((3456 - 3520 * z1 + 896 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                                             (-34752 + 33664 * z1 - 8144 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                                             (-18784 + 19888 * z1 - 5248 * z1 . pow(2)) * t_yis[0] . pow(5) +
                                                             (-2528 + 2972 * z1 - 854 * z1 . pow(2)) * t_yis[0] . pow(6) +
                                                             (28 - 8 * z1 - 3 * z1 . pow(2)) * t_yis[0] . pow(7)) + (8 - 4 * z1) * t_yis[0] . pow(8) +
                                                 ((-3584 + 2560 * z1 - 384 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                                  (-44432 + 42376 * z1 - 10080 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                                  (-34192 + 35260 * z1 - 9082 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                                  (-7152 + 7752 * z1 - 2088 * z1 . pow(2)) * t_yis[0] . pow(5) +
                                                  (20 + 20 * z1 - 15 * z1 . pow(2)) * t_yis[0] . pow(6)) * t_yis[1] . pow(2) +
                                                 (t_yis[0] * (-8832 + 7488 * z1 - 1536 * z1 . pow(2)) +
                                                  (-27360 + 25680 * z1 - 6000 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                                  (-30832 + 30996 * z1 - 7790 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                                  (-9888 + 10168 * z1 - 2612 * z1 . pow(2)) * t_yis[0] . pow(4) +
                                                  (-40 + 80 * z1 - 30 * z1 . pow(2)) * t_yis[0] . pow(5)) * t_yis[1] . pow(3) +
                                                 (-2816 + 2432 * z1 + t_yis[0] * (-7952 + 7368 * z1 - 1696 * z1 . pow(2)) - 512 * z1 . pow(2) +
                                                  (-13984 + 13708 * z1 - 3358 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                                  (-7112 + 7052 * z1 - 1748 * z1 . pow(2)) * t_yis[0] . pow(3) +
                                                  (-80 + 100 * z1 - 30 * z1 . pow(2)) * t_yis[0] . pow(4)) * t_yis[1] . pow(4) +
                                                 (-768 + 704 * z1 + t_yis[0] * (-2800 + 2668 * z1 - 634 * z1 . pow(2)) - 160 * z1 . pow(2) +
                                                  (-2496 + 2412 * z1 - 582 * z1 . pow(2)) * t_yis[0] . pow(2) +
                                                  (-52 + 56 * z1 - 15 * z1 . pow(2)) * t_yis[0] . pow(3)) * t_yis[1] . pow(5) +
                                                 (-160 + 144 * z1 + t_yis[0] * (-320 + 304 * z1 - 72 * z1 . pow(2)) - 32 * z1 . pow(2) +
                                                  (-12 + 12 * z1 - 3 * z1 . pow(2)) * t_yis[0] . pow(2)) * t_yis[1] . pow(6)));

    FFIntTemp den = cr_2 * (t_yis[2] . pow(10) * (192 - 48 * z1 + (288 - 72 * z1) * t_yis[0] + (-192 + 48 * z1) * t_yis[0] . pow(2)) +
                        t_yis[2] . pow(9) * (-2304 + 576 * z1 + (-3712 + 928 * z1) * t_yis[0] + (960 - 240 * z1) * t_yis[0] . pow(2) +
                                             t_yis[1] * (-1024 + 256 * z1 + (-1440 + 360 * z1) * t_yis[0] + (736 - 184 * z1) * t_yis[0] . pow(2)) +
                                             (736 - 184 * z1) * t_yis[0] . pow(3)) + t_yis[1] * ((1024 - 256 * z1) * t_yis[0] . pow(6) +
                                                 (512 - 128 * z1) * t_yis[0] . pow(7) + (64 - 16 * z1) * t_yis[0] . pow(8)) +
                        ((1024 - 256 * z1) * t_yis[0] . pow(4) + (4608 - 1152 * z1) * t_yis[0] . pow(5) +
                         (2624 - 656 * z1) * t_yis[0] . pow(6) + (384 - 96 * z1) * t_yis[0] . pow(7)) * t_yis[1] . pow(2) +
                        t_yis[2] . pow(8) * (9984 - 2496 * z1 + (17856 - 4464 * z1) * t_yis[0] +
                                             (5776 - 1444 * z1) * t_yis[0] . pow(2) + (-6328 + 1582 * z1) * t_yis[0] . pow(3) +
                                             t_yis[1] * (11904 - 2976 * z1 + (18336 - 4584 * z1) * t_yis[0] + (-3120 + 780 * z1) * t_yis[0] . pow(2) +
                                                         (-2272 + 568 * z1) * t_yis[0] . pow(3)) + (-1136 + 284 * z1) * t_yis[0] . pow(4) +
                                             (2288 - 572 * z1 + (3016 - 754 * z1) * t_yis[0] + (-1136 + 284 * z1) * t_yis[0] . pow(2)) *
                                             t_yis[1] . pow(2)) + ((3072 - 768 * z1) * t_yis[0] . pow(3) + (8192 - 2048 * z1) * t_yis[0] . pow(4) +
                                                                   (5440 - 1360 * z1) * t_yis[0] . pow(5) + (960 - 240 * z1) * t_yis[0] . pow(6)) * t_yis[1] . pow(3) +
                        t_yis[2] . pow(7) * (-18432 + 4608 * z1 + (-38656 + 9664 * z1) * t_yis[0] +
                                             (-46976 + 11744 * z1) * t_yis[0] . pow(2) + (9520 - 2380 * z1) * t_yis[0] . pow(3) +
                                             (11160 - 2790 * z1) * t_yis[0] . pow(4) + t_yis[1] * (-49408 + 12352 * z1 + (-87552 + 21888 * z1) * t_yis[0] +
                                                 (-26704 + 6676 * z1) * t_yis[0] . pow(2) + (18160 - 4540 * z1) * t_yis[0] . pow(3) +
                                                 (2712 - 678 * z1) * t_yis[0] . pow(4)) + (904 - 226 * z1) * t_yis[0] . pow(5) +
                                             (-25600 + 6400 * z1 + (-37360 + 9340 * z1) * t_yis[0] + (3576 - 894 * z1) * t_yis[0] . pow(2) +
                                              (2712 - 678 * z1) * t_yis[0] . pow(3)) * t_yis[1] . pow(2) +
                                             (-2768 + 692 * z1 + (-3424 + 856 * z1) * t_yis[0] + (904 - 226 * z1) * t_yis[0] . pow(2)) *
                                             t_yis[1] . pow(3)) + ((3072 - 768 * z1) * t_yis[0] . pow(2) + (7168 - 1792 * z1) * t_yis[0] . pow(3) +
                                                                   (5760 - 1440 * z1) * t_yis[0] . pow(4) + (1280 - 320 * z1) * t_yis[0] . pow(5)) * t_yis[1] . pow(4) +
                        t_yis[2] . pow(6) * (12288 - 3072 * z1 + (34560 - 8640 * z1) * t_yis[0] +
                                             (98560 - 24640 * z1) * t_yis[0] . pow(2) + (41840 - 10460 * z1) * t_yis[0] . pow(3) +
                                             (-30056 + 7514 * z1) * t_yis[0] . pow(4) + (-9512 + 2378 * z1) * t_yis[0] . pow(5) +
                                             t_yis[1] * (86016 - 21504 * z1 + (191104 - 47776 * z1) * t_yis[0] + (189088 - 47272 * z1) *
                                                         t_yis[0] . pow(2) + (-17264 + 4316 * z1) * t_yis[0] . pow(3) +
                                                         (-25128 + 6282 * z1) * t_yis[0] . pow(4) + (-1568 + 392 * z1) * t_yis[0] . pow(5)) +
                                             (-392 + 98 * z1) * t_yis[0] . pow(6) + (100608 - 25152 * z1 + (172784 - 43196 * z1) * t_yis[0] +
                                                 (51384 - 12846 * z1) * t_yis[0] . pow(2) + (-19448 + 4862 * z1) * t_yis[0] . pow(3) +
                                                 (-2352 + 588 * z1) * t_yis[0] . pow(4)) * t_yis[1] . pow(2) +
                                             (29568 - 7392 * z1 + (40544 - 10136 * z1) * t_yis[0] + (-1560 + 390 * z1) * t_yis[0] . pow(2) +
                                              (-1568 + 392 * z1) * t_yis[0] . pow(3)) * t_yis[1] . pow(3) +
                                             (1952 - 488 * z1 + (2272 - 568 * z1) * t_yis[0] + (-392 + 98 * z1) * t_yis[0] . pow(2)) * t_yis[1] . pow(4)) +
                        ((1024 - 256 * z1) * t_yis[0] + (3072 - 768 * z1) * t_yis[0] . pow(2) + (3200 - 800 * z1) * t_yis[0] . pow(3) +
                         (960 - 240 * z1) * t_yis[0] . pow(4)) * t_yis[1] . pow(5) +
                        t_yis[2] . pow(5) * ((-7168 + 1792 * z1) * t_yis[0] + (-58624 + 14656 * z1) * t_yis[0] . pow(2) +
                                             (-122816 + 30704 * z1) * t_yis[0] . pow(3) + (144 - 36 * z1) * t_yis[0] . pow(4) +
                                             (29192 - 7298 * z1) * t_yis[0] . pow(5) + (4432 - 1108 * z1) * t_yis[0] . pow(6) +
                                             t_yis[1] * (-53248 + 13312 * z1 + (-182016 + 45504 * z1) * t_yis[0] + (-379264 + 94816 * z1) *
                                                         t_yis[0] . pow(2) + (-153200 + 38300 * z1) * t_yis[0] . pow(3) +
                                                         (55456 - 13864 * z1) * t_yis[0] . pow(4) + (15944 - 3986 * z1) * t_yis[0] . pow(5) +
                                                         (440 - 110 * z1) * t_yis[0] . pow(6)) + (88 - 22 * z1) * t_yis[0] . pow(7) +
                                             (-161792 + 40448 * z1 + (-364352 + 91088 * z1) * t_yis[0] + (-310224 + 77556 * z1) *
                                              t_yis[0] . pow(2) + (-280 + 70 * z1) * t_yis[0] . pow(3) + (20360 - 5090 * z1) * t_yis[0] . pow(4) +
                                              (880 - 220 * z1) * t_yis[0] . pow(5)) * t_yis[1] . pow(2) +
                                             (-108288 + 27072 * z1 + (-176592 + 44148 * z1) * t_yis[0] + (-50960 + 12740 * z1) *
                                              t_yis[0] . pow(2) + (9736 - 2434 * z1) * t_yis[0] . pow(3) + (880 - 220 * z1) * t_yis[0] . pow(4)) *
                                             t_yis[1] . pow(3) + (-19712 + 4928 * z1 + (-25216 + 6304 * z1) * t_yis[0] +
                                                                  (8 - 2 * z1) * t_yis[0] . pow(2) + (440 - 110 * z1) * t_yis[0] . pow(3)) * t_yis[1] . pow(4) +
                                             (-800 + 200 * z1 + (-880 + 220 * z1) * t_yis[0] + (88 - 22 * z1) * t_yis[0] . pow(2)) * t_yis[1] . pow(5)) +
                        ((512 - 128 * z1) * t_yis[0] + (832 - 208 * z1) * t_yis[0] . pow(2) + (384 - 96 * z1) * t_yis[0] . pow(3)) *
                        t_yis[1] . pow(6) + t_yis[2] . pow(4) * ((-11264 + 2816 * z1) * t_yis[0] . pow(2) +
                                                                 (77056 - 19264 * z1) * t_yis[0] . pow(3) + (68160 - 17040 * z1) * t_yis[0] . pow(4) +
                                                                 (-18960 + 4740 * z1) * t_yis[0] . pow(5) + (-14168 + 3542 * z1) * t_yis[0] . pow(6) +
                                                                 (-1144 + 286 * z1) * t_yis[0] . pow(7) + t_yis[1] * ((54272 - 13568 * z1) * t_yis[0] +
                                                                     (249344 - 62336 * z1) * t_yis[0] . pow(2) + (358336 - 89584 * z1) * t_yis[0] . pow(3) +
                                                                     (44800 - 11200 * z1) * t_yis[0] . pow(4) + (-41384 + 10346 * z1) * t_yis[0] . pow(5) +
                                                                     (-5144 + 1286 * z1) * t_yis[0] . pow(6) + (-48 + 12 * z1) * t_yis[0] . pow(7)) +
                                                                 (-8 + 2 * z1) * t_yis[0] . pow(8) + (90112 - 22528 * z1 + (337920 - 84480 * z1) * t_yis[0] +
                                                                     (567488 - 141872 * z1) * t_yis[0] . pow(2) + (231376 - 57844 * z1) * t_yis[0] . pow(3) +
                                                                     (-30712 + 7678 * z1) * t_yis[0] . pow(4) + (-8952 + 2238 * z1) * t_yis[0] . pow(5) +
                                                                     (-120 + 30 * z1) * t_yis[0] . pow(6)) * t_yis[1] . pow(2) +
                                                                 (155648 - 38912 * z1 + (343104 - 85776 * z1) * t_yis[0] + (260064 - 65016 * z1) * t_yis[0] . pow(2) +
                                                                  (14664 - 3666 * z1) * t_yis[0] . pow(3) + (-7248 + 1812 * z1) * t_yis[0] . pow(4) +
                                                                  (-160 + 40 * z1) * t_yis[0] . pow(5)) * t_yis[1] . pow(3) +
                                                                 (65792 - 16448 * z1 + (100000 - 25000 * z1) * t_yis[0] + (26944 - 6736 * z1) * t_yis[0] . pow(2) +
                                                                  (-2312 + 578 * z1) * t_yis[0] . pow(3) + (-120 + 30 * z1) * t_yis[0] . pow(4)) * t_yis[1] . pow(4) +
                                                                 (7552 - 1888 * z1 + (8960 - 2240 * z1) * t_yis[0] + (168 - 42 * z1) * t_yis[0] . pow(2) +
                                                                  (-48 + 12 * z1) * t_yis[0] . pow(3)) * t_yis[1] . pow(5) +
                                                                 (176 - 44 * z1 + (184 - 46 * z1) * t_yis[0] + (-8 + 2 * z1) * t_yis[0] . pow(2)) * t_yis[1] . pow(6)) +
                        ((64 - 16 * z1) * t_yis[0] + (64 - 16 * z1) * t_yis[0] . pow(2)) * t_yis[1] . pow(7) +
                        t_yis[2] . pow(2) * ((-1024 + 256 * z1) * t_yis[0] . pow(4) + (11264 - 2816 * z1) * t_yis[0] . pow(5) +
                                             (1472 - 368 * z1) * t_yis[0] . pow(6) + (-2336 + 584 * z1) * t_yis[0] . pow(7) +
                                             (-528 + 132 * z1) * t_yis[0] . pow(8) + t_yis[1] * ((1024 - 256 * z1) * t_yis[0] . pow(3) +
                                                 (79872 - 19968 * z1) * t_yis[0] . pow(4) + (46400 - 11600 * z1) * t_yis[0] . pow(5) +
                                                 (-1408 + 352 * z1) * t_yis[0] . pow(6) + (-2456 + 614 * z1) * t_yis[0] . pow(7) +
                                                 (-48 + 12 * z1) * t_yis[0] . pow(8)) + (-8 + 2 * z1) * t_yis[0] . pow(9) +
                                             ((49152 - 12288 * z1) * t_yis[0] . pow(2) + (195328 - 48832 * z1) * t_yis[0] . pow(3) +
                                              (160064 - 40016 * z1) * t_yis[0] . pow(4) + (22000 - 5500 * z1) * t_yis[0] . pow(5) +
                                              (-4232 + 1058 * z1) * t_yis[0] . pow(6) + (-120 + 30 * z1) * t_yis[0] . pow(7)) * t_yis[1] . pow(2) +
                                             ((75776 - 18944 * z1) * t_yis[0] + (216576 - 54144 * z1) * t_yis[0] . pow(2) +
                                              (220352 - 55088 * z1) * t_yis[0] . pow(3) + (55488 - 13872 * z1) * t_yis[0] . pow(4) +
                                              (-2800 + 700 * z1) * t_yis[0] . pow(5) + (-160 + 40 * z1) * t_yis[0] . pow(6)) * t_yis[1] . pow(3) +
                                             (28672 - 7168 * z1 + (110336 - 27584 * z1) * t_yis[0] + (142592 - 35648 * z1) * t_yis[0] . pow(2) +
                                              (55552 - 13888 * z1) * t_yis[0] . pow(3) + (480 - 120 * z1) * t_yis[0] . pow(4) +
                                              (-120 + 30 * z1) * t_yis[0] . pow(5)) * t_yis[1] . pow(4) +
                                             (20480 - 5120 * z1 + (41216 - 10304 * z1) * t_yis[0] + (25792 - 6448 * z1) * t_yis[0] . pow(2) +
                                              (1672 - 418 * z1) * t_yis[0] . pow(3) + (-48 + 12 * z1) * t_yis[0] . pow(4)) * t_yis[1] . pow(5) +
                                             (3840 - 960 * z1 + (4784 - 1196 * z1) * t_yis[0] + (824 - 206 * z1) * t_yis[0] . pow(2) +
                                              (-8 + 2 * z1) * t_yis[0] . pow(3)) * t_yis[1] . pow(6) + (128 - 32 * z1 + (128 - 32 * z1) * t_yis[0]) *
                                             t_yis[1] . pow(7)) + t_yis[2] . pow(3) * ((7168 - 1792 * z1) * t_yis[0] . pow(3) +
                                                 (-44800 + 11200 * z1) * t_yis[0] . pow(4) + (-16704 + 4176 * z1) * t_yis[0] . pow(5) +
                                                 (10544 - 2636 * z1) * t_yis[0] . pow(6) + (3752 - 938 * z1) * t_yis[0] . pow(7) +
                                                 t_yis[1] * ((-14336 + 3584 * z1) * t_yis[0] . pow(2) + (-199936 + 49984 * z1) * t_yis[0] . pow(3) +
                                                             (-172544 + 43136 * z1) * t_yis[0] . pow(4) + (112 - 28 * z1) * t_yis[0] . pow(5) +
                                                             (14320 - 3580 * z1) * t_yis[0] . pow(6) + (808 - 202 * z1) * t_yis[0] . pow(7)) +
                                                 (152 - 38 * z1) * t_yis[0] . pow(8) + ((-104448 + 26112 * z1) * t_yis[0] +
                                                     (-359424 + 89856 * z1) * t_yis[0] . pow(2) + (-418752 + 104688 * z1) * t_yis[0] . pow(3) +
                                                     (-89328 + 22332 * z1) * t_yis[0] . pow(4) + (18144 - 4536 * z1) * t_yis[0] . pow(5) +
                                                     (1744 - 436 * z1) * t_yis[0] . pow(6)) * t_yis[1] . pow(2) +
                                                 (-73728 + 18432 * z1 + (-284160 + 71040 * z1) * t_yis[0] + (-408960 + 102240 * z1) *
                                                  t_yis[0] . pow(2) + (-164752 + 41188 * z1) * t_yis[0] . pow(3) +
                                                  (4400 - 1100 * z1) * t_yis[0] . pow(4) + (1920 - 480 * z1) * t_yis[0] . pow(5)) * t_yis[1] . pow(3) +
                                                 (-79872 + 19968 * z1 + (-168320 + 42080 * z1) * t_yis[0] + (-115328 + 28832 * z1) *
                                                  t_yis[0] . pow(2) + (-8760 + 2190 * z1) * t_yis[0] . pow(3) + (1080 - 270 * z1) * t_yis[0] . pow(4)) *
                                                 t_yis[1] . pow(4) + (-22272 + 5568 * z1 + (-31008 + 7752 * z1) * t_yis[0] +
                                                                      (-7248 + 1812 * z1) * t_yis[0] . pow(2) + (232 - 58 * z1) * t_yis[0] . pow(3)) * t_yis[1] . pow(5) +
                                                 (-1536 + 384 * z1 + (-1680 + 420 * z1) * t_yis[0] + (-32 + 8 * z1) * t_yis[0] . pow(2)) * t_yis[1] . pow(6) +
                                                 (-16 + 4 * z1 + (-16 + 4 * z1) * t_yis[0]) * t_yis[1] . pow(7)) +
                        t_yis[2] * ((-1024 + 256 * z1) * t_yis[0] . pow(6) + (192 - 48 * z1) * t_yis[0] . pow(8) +
                                    t_yis[1] * ((-14848 + 3712 * z1) * t_yis[0] . pow(5) + (-7168 + 1792 * z1) * t_yis[0] . pow(6) +
                                                (-160 + 40 * z1) * t_yis[0] . pow(7) + (176 - 44 * z1) * t_yis[0] . pow(8)) +
                                    (32 - 8 * z1) * t_yis[0] . pow(9) + ((-11264 + 2816 * z1) * t_yis[0] . pow(3) +
                                                                         (-49408 + 12352 * z1) * t_yis[0] . pow(4) + (-31808 + 7952 * z1) * t_yis[0] . pow(5) +
                                                                         (-3952 + 988 * z1) * t_yis[0] . pow(6) + (384 - 96 * z1) * t_yis[0] . pow(7)) * t_yis[1] . pow(2) +
                                    ((-26624 + 6656 * z1) * t_yis[0] . pow(2) + (-71424 + 17856 * z1) * t_yis[0] . pow(3) +
                                     (-55808 + 13952 * z1) * t_yis[0] . pow(4) + (-10640 + 2660 * z1) * t_yis[0] . pow(5) +
                                     (400 - 100 * z1) * t_yis[0] . pow(6)) * t_yis[1] . pow(3) +
                                    ((-19456 + 4864 * z1) * t_yis[0] + (-50944 + 12736 * z1) * t_yis[0] . pow(2) +
                                     (-48512 + 12128 * z1) * t_yis[0] . pow(3) + (-12960 + 3240 * z1) * t_yis[0] . pow(4) +
                                     (160 - 40 * z1) * t_yis[0] . pow(5)) * t_yis[1] . pow(4) +
                                    (-4096 + 1024 * z1 + (-17152 + 4288 * z1) * t_yis[0] + (-21248 + 5312 * z1) * t_yis[0] . pow(2) +
                                     (-8128 + 2032 * z1) * t_yis[0] . pow(3) + (-48 + 12 * z1) * t_yis[0] . pow(4)) * t_yis[1] . pow(5) +
                                    (-2048 + 512 * z1 + (-4160 + 1040 * z1) * t_yis[0] + (-2480 + 620 * z1) * t_yis[0] . pow(2) +
                                     (-64 + 16 * z1) * t_yis[0] . pow(3)) * t_yis[1] . pow(6) +
                                    (-256 + 64 * z1 + (-272 + 68 * z1) * t_yis[0] + (-16 + 4 * z1) * t_yis[0] . pow(2)) * t_yis[1] . pow(7)));

    return num / den;
  }
  /**
   *  examples for the reconstruction of a polynomial with 3 variables
   */
  template<typename FFIntTemp>
  FFIntTemp pol_n_eq_3(const std::vector<FFIntTemp>& yis) {
    return yis[0].pow(5) + yis[0] * yis[1].pow(4) + yis[0] * yis[1] * yis[2].pow(3) + yis[1].pow(5);
  }
  /**
   *  example for the reconstruction of a sparse polynomial with 20 variables
   */
  template<typename FFIntTemp>
  FFIntTemp pol_20_20(const std::vector<FFIntTemp>& yis) {
    FFIntTemp result(0);
    for(size_t i = 0; i < 20; i++){
      result += yis[i].pow(20);
    }
    return result;
  }
  /**
   *  example for reconstruction of a polynomial with 9 variables
   */
  template<typename FFIntTemp>
  FFIntTemp pol_1(const std::vector<FFIntTemp>& yis) {
    return yis[0].pow(2) * yis[2].pow(3) * yis[3] * yis[5] * yis[7] * yis[8].pow(2) + yis[0] * yis[1] * yis[2] * yis[3].pow(2) * yis[4].pow(2) * yis[7] * yis[8] + yis[1] * yis[2] * yis[3] * yis[4].pow(2) * yis[7] * yis[8] + yis[0] * yis[2].pow(3) * yis[3].pow(2) * yis[4].pow(2) * yis[5].pow(2) * yis[6] * yis[7].pow(2) + yis[1] * yis[2] * yis[3] * yis[4].pow(2) * yis[5] * yis[6] * yis[7].pow(2);
  }
  /**
   *  example for reconstruction of a more complex polynomial with 5 variables
   */
  template<typename FFIntTemp>
  FFIntTemp pol_6(const std::vector<FFIntTemp>& yis) {
    FFIntTemp result(0);
    for(size_t i = 0; i < 6; i++){
      result += yis[0].pow(i) * yis[1].pow(i) * yis[2].pow(i) * yis[3].pow(i) * yis[4].pow(i);
    }
    return result;
  }
  /**
   *  example for reconstuction of a polynomial with 5 variables
   */
  template<typename FFIntTemp>
  FFIntTemp pol_7(const std::vector<FFIntTemp>& yis) {
    FFIntTemp result(0);
    for(size_t i = 0; i < 5; i++){
      FFIntTemp a;
      a = yis[0] + yis[1] + yis[2] + yis[3] + yis[4];
      result += a.pow(i);
    }
    return result;
  }
  /**
   *  example for a three loop gg -> h integral coefficient
   */
  template<typename FFIntTemp>
  FFIntTemp ggh(const std::vector<FFIntTemp>& yis) {
    FFIntTemp num = -4444263936000 + 26953049894400 * yis[0] - 68516061805440 * yis[0] . pow(2) +
                89743434192000 * yis[0] . pow(3) - 49560953693040 * yis[0] . pow(4) -
                28729977426000 * yis[0] . pow(5) + 79757115103800 * yis[0] . pow(6) -
                75482789911620 * yis[0] . pow(7) + 42206858356455 * yis[0] . pow(8) -
                13975295550210 * yis[0] . pow(9) + 1432566459740 * yis[0] . pow(10) +
                1196179610030 * yis[0] . pow(11) - 804996323345 * yis[0] . pow(12) +
                282889964900 * yis[0] . pow(13) - 68550339740 * yis[0] . pow(14) +
                12325429400 * yis[0] . pow(15) - 1698275115 * yis[0] . pow(16) + 182775430 * yis[0] . pow(17) -
                15498340 * yis[0] . pow(18) + 1020190 * yis[0] . pow(19) - 48955 * yis[0] . pow(20) +
                1480 * yis[0] . pow(21) - 20 * yis[0] . pow(22) +
                yis[1] * (-530348829696000 + 2062599499900800 * yis[0] - 2706007765381920 * yis[0] . pow(2) +
                          197480772584640 * yis[0] . pow(3) + 3446064469459152 * yis[0] . pow(4) -
                          4255234347788424 * yis[0] . pow(5) + 1973738647850166 * yis[0] . pow(6) +
                          516150155040984 * yis[0] . pow(7) - 1404315632302002 * yis[0] . pow(8) +
                          1081782645683654 * yis[0] . pow(9) - 525982628267780 * yis[0] . pow(10) +
                          184906655524766 * yis[0] . pow(11) - 50139901610736 * yis[0] . pow(12) +
                          11173502324412 * yis[0] . pow(13) - 2229044486286 * yis[0] . pow(14) +
                          430053951128 * yis[0] . pow(15) - 79736298830 * yis[0] . pow(16) +
                          13049440610 * yis[0] . pow(17) - 1728111252 * yis[0] . pow(18) + 174314610 * yis[0] . pow(19) -
                          12747984 * yis[0] . pow(20) + 635348 * yis[0] . pow(21) - 19328 * yis[0] . pow(22) +
                          272 * yis[0] . pow(23)) + (5172065063424000 - 24484758091568640 * yis[0] +
                                                     47559144062609280 * yis[0] . pow(2) - 44670237145355520 * yis[0] . pow(3) +
                                                     9865157601595824 * yis[0] . pow(4) + 26470011649875696 * yis[0] . pow(5) -
                                                     38130147382831848 * yis[0] . pow(6) + 28504663248932100 * yis[0] . pow(7) -
                                                     14339532548611047 * yis[0] . pow(8) + 5215752044491378 * yis[0] . pow(9) -
                                                     1403669534036122 * yis[0] . pow(10) + 272863632469908 * yis[0] . pow(11) -
                                                     29847088222883 * yis[0] . pow(12) - 4216274980770 * yis[0] . pow(13) +
                                                     3890875305834 * yis[0] . pow(14) - 1455307267904 * yis[0] . pow(15) +
                                                     379024311431 * yis[0] . pow(16) - 73449196518 * yis[0] . pow(17) +
                                                     10693070590 * yis[0] . pow(18) - 1161337776 * yis[0] . pow(19) + 92322075 * yis[0] . pow(20) -
                                                     5188586 * yis[0] . pow(21) + 193706 * yis[0] . pow(22) - 4248 * yis[0] . pow(23) +
                                                     40 * yis[0] . pow(24)) * yis[1] . pow(2) + (-13384429922304000 + 66386168799828480 * yis[0] -
                                                         137767432678744320 * yis[0] . pow(2) + 145396530815220864 * yis[0] . pow(3) -
                                                         56062023313623840 * yis[0] . pow(4) - 57664165857341184 * yis[0] . pow(5) +
                                                         111681955232154720 * yis[0] . pow(6) - 97365254494851048 * yis[0] . pow(7) +
                                                         57087566862186414 * yis[0] . pow(8) - 24885777495898250 * yis[0] . pow(9) +
                                                         8428579447891252 * yis[0] . pow(10) - 2260719976619262 * yis[0] . pow(11) +
                                                         474985413678842 * yis[0] . pow(12) - 71389387652756 * yis[0] . pow(13) +
                                                         4179544520084 * yis[0] . pow(14) + 1602746583044 * yis[0] . pow(15) -
                                                         681317919630 * yis[0] . pow(16) + 152671753078 * yis[0] . pow(17) -
                                                         23759672084 * yis[0] . pow(18) + 2699940770 * yis[0] . pow(19) - 224310602 * yis[0] . pow(20) +
                                                         13300392 * yis[0] . pow(21) - 533812 * yis[0] . pow(22) + 12992 * yis[0] . pow(23) -
                                                         144 * yis[0] . pow(24)) * yis[1] . pow(3) +
                (13260837058560000 - 66924650100648960 * yis[0] + 142128374611461120 * yis[0] . pow(2) -
                 155538595935230976 * yis[0] . pow(3) + 67290968387975040 * yis[0] . pow(4) +
                 53260608978113472 * yis[0] . pow(5) - 116440937764968384 * yis[0] . pow(6) +
                 106618064162228640 * yis[0] . pow(7) - 65202141858635160 * yis[0] . pow(8) +
                 29644869763978208 * yis[0] . pow(9) - 10484920438803016 * yis[0] . pow(10) +
                 2943461646643734 * yis[0] . pow(11) - 653585780258358 * yis[0] . pow(12) +
                 109002691281894 * yis[0] . pow(13) - 10727501182698 * yis[0] . pow(14) -
                 628145179492 * yis[0] . pow(15) + 553831157684 * yis[0] . pow(16) -
                 138174546780 * yis[0] . pow(17) + 22432326660 * yis[0] . pow(18) -
                 2618117178 * yis[0] . pow(19) + 222878226 * yis[0] . pow(20) - 13583434 * yis[0] . pow(21) +
                 563918 * yis[0] . pow(22) - 14328 * yis[0] . pow(23) + 168 * yis[0] . pow(24)) * yis[1] . pow(4) +
                (-4378234871808000 + 22194408939632640 * yis[0] - 47399557295324160 * yis[0] . pow(2) +
                 52256914945516800 * yis[0] . pow(3) - 22898884802466816 * yis[0] . pow(4) -
                 18025212963606336 * yis[0] . pow(5) + 40141269737887680 * yis[0] . pow(6) -
                 37367333013042576 * yis[0] . pow(7) + 23282123645804256 * yis[0] . pow(8) -
                 10809082846095108 * yis[0] . pow(9) + 3910859363288104 * yis[0] . pow(10) -
                 1124875529723404 * yis[0] . pow(11) + 256770669143612 * yis[0] . pow(12) -
                 44649611370420 * yis[0] . pow(13) + 4975030595868 * yis[0] . pow(14) +
                 28562694348 * yis[0] . pow(15) - 168339427772 * yis[0] . pow(16) +
                 45516044888 * yis[0] . pow(17) - 7606048708 * yis[0] . pow(18) + 903336000 * yis[0] . pow(19) -
                 78066096 * yis[0] . pow(20) + 4834416 * yis[0] . pow(21) - 204464 * yis[0] . pow(22) +
                 5312 * yis[0] . pow(23) - 64 * yis[0] . pow(24)) * yis[1] . pow(5) +
                (-135444234240000 + 716861362913280 * yis[0] - 1668354406133760 * yis[0] . pow(2) +
                 2165148665137152 * yis[0] . pow(3) - 1529979184949760 * yis[0] . pow(4) +
                 235299216266496 * yis[0] . pow(5) + 681417157977216 * yis[0] . pow(6) -
                 821414012414400 * yis[0] . pow(7) + 530483275717344 * yis[0] . pow(8) -
                 232393058402832 * yis[0] . pow(9) + 73392774099912 * yis[0] . pow(10) -
                 16874228783772 * yis[0] . pow(11) + 2724383055828 * yis[0] . pow(12) -
                 261137653860 * yis[0] . pow(13) - 1578650412 * yis[0] . pow(14) +
                 5660299956 * yis[0] . pow(15) - 1144397988 * yis[0] . pow(16) + 139612932 * yis[0] . pow(17) -
                 11696556 * yis[0] . pow(18) + 666264 * yis[0] . pow(19) - 23424 * yis[0] . pow(20) +
                 384 * yis[0] . pow(21)) * yis[1] . pow(6) +
                (25333088256000 * yis[0] - 89985613824000 * yis[0] . pow(2) + 122002697871360 * yis[0] . pow(3) -
                 70012594114560 * yis[0] . pow(4) - 4686236328960 * yis[0] . pow(5) +
                 33226392545280 * yis[0] . pow(6) - 21784577299200 * yis[0] . pow(7) +
                 5636087527680 * yis[0] . pow(8) + 908612884800 * yis[0] . pow(9) -
                 1353938860800 * yis[0] . pow(10) + 585698020560 * yis[0] . pow(11) -
                 156226500000 * yis[0] . pow(12) + 28765567440 * yis[0] . pow(13) -
                 3751866720 * yis[0] . pow(14) + 342704880 * yis[0] . pow(15) - 20941920 * yis[0] . pow(16) +
                 771120 * yis[0] . pow(17) - 12960 * yis[0] . pow(18)) * yis[1] . pow(7);

    FFIntTemp den =  yis[1] * (-11851370496000 * yis[0] + 47019842150400 * yis[0] . pow(2) -
                           63221373296640 * yis[0] . pow(3) + 12145773588480 * yis[0] . pow(4) +
                           56468249487360 * yis[0] . pow(5) - 65585862673920 * yis[0] . pow(6) +
                           27466213639680 * yis[0] . pow(7) + 2647003092480 * yis[0] . pow(8) -
                           8231281958400 * yis[0] . pow(9) + 3890701102080 * yis[0] . pow(10) -
                           726722104320 * yis[0] . pow(11) - 89881420800 * yis[0] . pow(12) +
                           91522821120 * yis[0] . pow(13) - 27040550400 * yis[0] . pow(14) +
                           4723330560 * yis[0] . pow(15) - 533660160 * yis[0] . pow(16) + 38545920 * yis[0] . pow(17) -
                           1628160 * yis[0] . pow(18) + 30720 * yis[0] . pow(19)) +
                 (23702740992000 * yis[0] - 94039684300800 * yis[0] . pow(2) + 126442746593280 * yis[0] . pow(3) -
                  24291547176960 * yis[0] . pow(4) - 112936498974720 * yis[0] . pow(5) +
                  131171725347840 * yis[0] . pow(6) - 54932427279360 * yis[0] . pow(7) -
                  5294006184960 * yis[0] . pow(8) + 16462563916800 * yis[0] . pow(9) -
                  7781402204160 * yis[0] . pow(10) + 1453444208640 * yis[0] . pow(11) +
                  179762841600 * yis[0] . pow(12) - 183045642240 * yis[0] . pow(13) +
                  54081100800 * yis[0] . pow(14) - 9446661120 * yis[0] . pow(15) +
                  1067320320 * yis[0] . pow(16) - 77091840 * yis[0] . pow(17) + 3256320 * yis[0] . pow(18) -
                  61440 * yis[0] . pow(19)) * yis[1] . pow(2) +
                 (-11851370496000 * yis[0] + 47019842150400 * yis[0] . pow(2) - 63221373296640 * yis[0] . pow(3) +
                  12145773588480 * yis[0] . pow(4) + 56468249487360 * yis[0] . pow(5) -
                  65585862673920 * yis[0] . pow(6) + 27466213639680 * yis[0] . pow(7) +
                  2647003092480 * yis[0] . pow(8) - 8231281958400 * yis[0] . pow(9) +
                  3890701102080 * yis[0] . pow(10) - 726722104320 * yis[0] . pow(11) -
                  89881420800 * yis[0] . pow(12) + 91522821120 * yis[0] . pow(13) -
                  27040550400 * yis[0] . pow(14) + 4723330560 * yis[0] . pow(15) - 533660160 * yis[0] . pow(16) +
                  38545920 * yis[0] . pow(17) - 1628160 * yis[0] . pow(18) + 30720 * yis[0] . pow(19)) *
                 yis[1] . pow(3);

    return num / den;
  }
  /**
   *  a benchmark function with 20 variables
   */
  template<typename FFIntTemp>
  FFIntTemp bench_1(const std::vector<FFIntTemp>& yis) {
    FFIntTemp num = yis[0].pow(20) + yis[1].pow(20) + yis[2].pow(20) + yis[3].pow(20) +
                yis[4].pow(20) + yis[5].pow(20) + yis[6].pow(20) + yis[7].pow(20) +
                yis[8].pow(20) + yis[9].pow(20) + yis[10].pow(20) + yis[11].pow(20) +
                yis[12].pow(20) + yis[13].pow(20) + yis[14].pow(20) + yis[15].pow(20) +
                yis[16].pow(20) + yis[17].pow(20) + yis[18].pow(20) + yis[19].pow(20);

    FFIntTemp den = yis[19].pow(35) * (yis[0] * yis[1] + yis[2] * yis[3] + yis[4] * yis[5] +
                                   (yis[0] * yis[1] + yis[2] * yis[3] + yis[4] * yis[5]).pow(2) +
                                   (yis[0] * yis[1] + yis[2] * yis[3] + yis[4] * yis[5]).pow(3) +
                                   (yis[0] * yis[1] + yis[2] * yis[3] + yis[4] * yis[5]).pow(4) +
                                   (yis[0] * yis[1] + yis[2] * yis[3] + yis[4] * yis[5]).pow(5));

    return num / den;
  }
  /**
   *  a benchmark function with 5 variables and almost complete dense numerator
   */
  template<typename FFIntTemp>
  FFIntTemp bench_2(const std::vector<FFIntTemp>& yis) {
    mpz_class cr_1_mpz;
    cr_1_mpz = "123456789109898799879870980";
    FFIntTemp num = cr_1_mpz * ((1 + yis[0] + yis[1] + yis[2] + yis[3] + yis[4]).pow(20) - 1);
    FFIntTemp den = yis[3] - yis[1] + (yis[0] * yis[1] * yis[2] * yis[3] * yis[4]).pow(10);

    return num / den;
  }
  /**
   *  a benchmark function with 5 variables and almost complete sparse with high degrees
   */
  template<typename FFIntTemp>
  FFIntTemp bench_3(const std::vector<FFIntTemp>& yis) {
    FFIntTemp num = yis[0].pow(100) + yis[1].pow(200) + yis[2].pow(300);
    FFIntTemp den = yis[0] * yis[1] * yis[2] * yis[3] * yis[4] + (yis[0] * yis[1] * yis[2] * yis[3] * yis[4]).pow(4);

    return num / den;
  }
}
