//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2020  Jonas Klappert and Fabian Lange
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

#include "firefly/config.hpp"

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
