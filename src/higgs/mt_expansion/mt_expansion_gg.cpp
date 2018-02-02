#include "mt_expansion_gg.h"
#include "constants.h"
#include <iostream>



//: coefficients taken from http://arxiv.org/pdf/0909.3420v1.pdf
//: note that the sigma(z)/z vs sigma(z) difference in what is expanded
//: between ours and Harlander et al.'s conventions leads to differences
//: in the regular part: us = (har_reg + b_k * Log[1-z]^k ) /z
//: where b_k is the plus coefficient at the same order in a_s and rho
//: where rho = mh^2/mt^2


//:  analytic results from Robert Harlander et al, http://arxiv.org/abs/0912.2104
//:  code generation and checks at harlander_mt_expansion_nnlo.m




HiggsMtExpansion_gg_nnlo_reg::HiggsMtExpansion_gg_nnlo_reg(const double& rho){
    _rho=rho;
    set_dimensions(2);
    _channel = "gg";
    _name="gg_reg_nnlo_mt_exp";
   
    
    
    cL0[0][0][0]=-308.2285739011688;
    cL0[1][0][0]=1.511111111111111;
    cL0[2][0][0]=0.187989417989418;
    cL0[3][0][0]=0.02890025825144873;
    cL0[0][1][0]=8.210732683007046;
    cL0[1][1][0]=-3.022222222222222;
    cL0[2][1][0]=-0.375978835978836;
    cL0[3][1][0]=-0.05780051650289746;
    cL0[0][2][0]=162.5;
    cL0[1][2][0]=0;
    cL0[2][2][0]=0;
    cL0[3][2][0]=0;
    cL0[0][3][0]=-72.;
    cL0[1][3][0]=0;
    cL0[2][3][0]=0;
    cL0[3][3][0]=0;
    cL0[0][0][1]=632.0593500547479;
    cL0[1][0][1]=-5.048388415853759 - 0.675*log(rho);
    cL0[2][0][1]=-0.211833190859022 + 0.04205357142857143*log(rho);
    cL0[3][0][1]=0.009836756964496915 + 0.0117046130952381*log(rho);
    cL0[0][1][1]=632.8678019509789;
    cL0[1][1][1]=11.21666666666667;
    cL0[2][1][1]=1.005198412698413;
    cL0[3][1][1]=0.1375463907785336;
    cL0[0][2][1]=-559.5833333333333;
    cL0[1][2][1]=-1.0125;
    cL0[2][2][1]=0.0590625;
    cL0[3][2][1]=0.01704575892857143;
    cL0[0][3][1]=216.;
    cL0[1][3][1]=0;
    cL0[2][3][1]=0;
    cL0[3][3][1]=0;
    cL0[0][0][2]=-151.4055684138741;
    cL0[1][0][2]=9.222578230668574 + 0.675*log(rho);
    cL0[2][0][2]=0.7859614509917534 + 0.03830357142857143*log(rho);
    cL0[3][0][2]=0.14394143162835 + 0.0107953869047619*log(rho);
    cL0[0][1][2]=-1026.578534633986;
    cL0[1][1][2]=-10.38055555555556;
    cL0[2][1][2]=-0.8839384920634921;
    cL0[3][1][2]=-0.1450344269652305;
    cL0[0][2][2]=652.8333333333333;
    cL0[1][2][2]=1.0125;
    cL0[2][2][2]=0.04941964285714286;
    cL0[3][2][2]=0.01318861607142857;
    cL0[0][3][2]=-144.;
    cL0[1][3][2]=0;
    cL0[2][3][2]=0;
    cL0[3][3][2]=0;
    cL0[0][0][3]=-400.9408717239148;
    cL0[1][0][3]=-3.614506172839506;
    cL0[2][0][3]=0.01824722561512172 + 0.08598214285714286*log(rho);
    cL0[3][0][3]=0.202130209711551 + 0.04088045634920635*log(rho);
    cL0[0][1][3]=988.789267316993;
    cL0[1][1][3]=3.118518518518519;
    cL0[2][1][3]=0.1272475749559083;
    cL0[3][1][3]=-0.08478650951121189;
    cL0[0][2][3]=-417.2083333333333;
    cL0[1][2][3]=0;
    cL0[2][2][3]=0.1131696428571429;
    cL0[3][2][3]=0.05118303571428571;
    cL0[0][3][3]=72.;
    cL0[1][3][3]=0;
    cL0[2][3][3]=0;
    cL0[3][3][3]=0;
    cL0[0][0][4]=306.1662673448907;
    cL0[1][0][4]=0.8828742283950617;
    cL0[2][0][4]=0.5253703100151366 + 0.105*log(rho);
    cL0[3][0][4]=0.4714274589338589 + 0.07378075396825397*log(rho);
    cL0[0][1][4]=-296.1555555555556;
    cL0[1][1][4]=-0.2467592592592593;
    cL0[2][1][4]=-0.2496698633156966;
    cL0[3][1][4]=-0.2500455752131099;
    cL0[0][2][4]=97.65;
    cL0[1][2][4]=0;
    cL0[2][2][4]=0.1359375;
    cL0[3][2][4]=0.0894921875;
    cL0[0][3][4]=0;
    cL0[1][3][4]=0;
    cL0[2][3][4]=0;
    cL0[3][3][4]=0;
    cL0[0][0][5]=44.78165540105065;
    cL0[1][0][5]=0.5918090277777778;
    cL0[2][0][5]=0.4344150338000721 + 0.1240178571428571*log(rho);
    cL0[3][0][5]=0.7057086513364842 + 0.1179893849206349*log(rho);
    cL0[0][1][5]=-65.29925925925926;
    cL0[1][1][5]=-0.4192592592592593;
    cL0[2][1][5]=-0.2275125661375661;
    cL0[3][1][5]=-0.3892577123750735;
    cL0[0][2][5]=46.64722222222222;
    cL0[1][2][5]=0;
    cL0[2][2][5]=0.1587053571428571;
    cL0[3][2][5]=0.1407388392857143;
    cL0[0][3][5]=0;
    cL0[1][3][5]=0;
    cL0[2][3][5]=0;
    cL0[3][3][5]=0;
    cL0[0][0][6]=13.78740480739543;
    cL0[1][0][6]=0.4293898625808348;
    cL0[2][0][6]=0.4127924037787732 + 0.1430357142857143*log(rho);
    cL0[3][0][6]=1.002150166791469 + 0.1735063492063492*log(rho);
    cL0[0][1][6]=-20.80719576719577;
    cL0[1][1][6]=-0.4785493827160494;
    cL0[2][1][6]=-0.2107628600823045;
    cL0[3][1][6]=-0.5622807772388791;
    cL0[0][2][6]=31.70079365079365;
    cL0[1][2][6]=0;
    cL0[2][2][6]=0.1814732142857143;
    cL0[3][2][6]=0.2049229910714286;
    cL0[0][3][6]=0;
    cL0[1][3][6]=0;
    cL0[2][3][6]=0;
    cL0[3][3][6]=0;
    cL0[0][0][7]=6.807004782702873;
    cL0[1][0][7]=0.3148171348786428;
    cL0[2][0][7]=0.4068875961599876 + 0.1620535714285714*log(rho);
    cL0[3][0][7]=1.344974248461683 + 0.2403316468253968*log(rho);
    cL0[0][1][7]=-5.701746031746032;
    cL0[1][1][7]=-0.5047266313932981;
    cL0[2][1][7]=-0.194524444654405;
    cL0[3][1][7]=-0.7654912475329638;
    cL0[0][2][7]=24.53869047619048;
    cL0[1][2][7]=0;
    cL0[2][2][7]=0.2042410714285714;
    cL0[3][2][7]=0.2820446428571429;
    cL0[0][3][7]=0;
    cL0[1][3][7]=0;
    cL0[2][3][7]=0;
    cL0[3][3][7]=0;
    cL0[0][0][8]=4.565383193639303;
    cL0[1][0][8]=0.2278496367682876;
    cL0[2][0][8]=0.4068642770127756 + 0.1810714285714286*log(rho);
    cL0[3][0][8]=1.727265866866938 + 0.3184652777777778*log(rho);
    cL0[0][1][8]=0.9854497354497354;
    cL0[1][1][8]=-0.5181988536155203;
    cL0[2][1][8]=-0.1775945767195767;
    cL0[3][1][8]=-0.9966638525613435;
    cL0[0][2][8]=20.25515873015873;
    cL0[1][2][8]=0;
    cL0[2][2][8]=0.2270089285714286;
    cL0[3][2][8]=0.3721037946428571;
    cL0[0][3][8]=0;
    cL0[1][3][8]=0;
    cL0[2][3][8]=0;
    cL0[3][3][8]=0;
    cL0[0][0][9]=3.607734433514778;
    cL0[1][0][9]=0.158602726177062;
    cL0[2][0][9]=0.4095061051164047 + 0.2000892857142857*log(rho);
    cL0[3][0][9]=2.143907060238031 + 0.4079072420634921*log(rho);
    cL0[0][1][9]=4.435570252792475;
    cL0[1][1][9]=-0.525896531452087;
    cL0[2][1][9]=-0.1595956090255032;
    cL0[3][1][9]=-1.25407618134874;
    cL0[0][2][9]=17.35906084656085;
    cL0[1][2][9]=0;
    cL0[2][2][9]=0.2497767857142857;
    cL0[3][2][9]=0.4751004464285714;
    cL0[0][3][9]=0;
    cL0[1][3][9]=0;
    cL0[2][3][9]=0;
    cL0[3][3][9]=0;
    cL0[0][0][10]=3.088871574053079;
    cL0[1][0][10]=0.1014963577481903;
    cL0[2][0][10]=0.4134202578815912 + 0.2191071428571429*log(rho);
    cL0[3][0][10]=2.590574664651477 + 0.5086575396825397*log(rho);
    cL0[0][1][10]=6.391172037838705;
    cL0[1][1][10]=-0.5306430041152263;
    cL0[2][1][10]=-0.1404063252358557;
    cL0[3][1][10]=-1.536280414313359;
    cL0[0][2][10]=15.24650072150072;
    cL0[1][2][10]=0;
    cL0[2][2][10]=0.2725446428571429;
    cL0[3][2][10]=0.5910345982142857;
    cL0[0][3][10]=0;
    cL0[1][3][10]=0;
    cL0[2][3][10]=0;
    cL0[3][3][10]=0;
    cL0[0][0][11]=2.756928843977676;
    cL0[1][0][11]=0.0531300879904275;
    cL0[2][0][11]=0.4179048648914258 + 0.238125*log(rho);
    cL0[3][0][11]=3.063456760875681 + 0.6207161706349206*log(rho);
    cL0[0][1][11]=7.564248825435694;
    cL0[1][1][11]=-0.5337448559670782;
    cL0[2][1][11]=-0.1200023532280477;
    cL0[3][1][11]=-1.842012598810022;
    cL0[0][2][11]=13.62495791245791;
    cL0[1][2][11]=0;
    cL0[2][2][11]=0.2953125;
    cL0[3][2][11]=0.71990625;
    cL0[0][3][11]=0;
    cL0[1][3][11]=0;
    cL0[2][3][11]=0;
    cL0[3][3][11]=0;
    cL0[0][0][12]=2.523129406994082;
    cL0[1][0][12]=0.01130858010613903;
    cL0[2][0][12]=0.4225744607129024 + 0.2571428571428571*log(rho);
    cL0[3][0][12]=3.559127745429376 + 0.7440831349206349*log(rho);
    cL0[0][1][12]=8.288616972253336;
    cL0[1][1][12]=-0.5358670033670034;
    cL0[2][1][12]=-0.0984000793078174;
    cL0[3][1][12]=-2.170145161101758;
    cL0[0][2][12]=12.33430458430458;
    cL0[1][2][12]=0;
    cL0[2][2][12]=0.3180803571428571;
    cL0[3][2][12]=0.8617154017857143;
    cL0[0][3][12]=0;
    cL0[1][3][12]=0;
    cL0[2][3][12]=0;
    cL0[3][3][12]=0;
    cL0[0][0][13]=2.349983593799101;
    cL0[1][0][13]=-0.02545356402311369;
    cL0[2][0][13]=0.4272066980817023 + 0.2761607142857143*log(rho);
    cL0[3][0][13]=4.074475384124231 + 0.8787584325396825*log(rho);
    cL0[0][1][13]=8.737888599634354;
    cL0[1][1][13]=-0.5373737373737374;
    cL0[2][1][13]=-0.07563370200126153;
    cL0[3][1][13]=-2.519657610394639;
    cL0[0][2][13]=11.2787962037962;
    cL0[1][2][13]=0;
    cL0[2][2][13]=0.3408482142857143;
    cL0[3][2][13]=1.016462053571429;
    cL0[0][3][13]=0;
    cL0[1][3][13]=0;
    cL0[2][3][13]=0;
    cL0[3][3][13]=0;
    cL0[0][0][14]=2.218315895643321;
    cL0[1][0][14]=-0.05820125358020912;
    cL0[2][0][14]=0.4316701268514003 + 0.2951785714285714*log(rho);
    cL0[3][0][14]=4.606650741751972 + 1.024742063492063*log(rho);
    cL0[0][1][14]=9.009954042362467;
    cL0[1][1][14]=-0.538476800976801;
    cL0[2][1][14]=-0.05174483618233618;
    cL0[3][1][14]=-2.889616586871572;
    cL0[0][2][14]=10.39726384726385;
    cL0[1][2][14]=0;
    cL0[2][2][14]=0.3636160714285714;
    cL0[3][2][14]=1.184146205357143;
    cL0[0][3][14]=0;
    cL0[1][3][14]=0;
    cL0[2][3][14]=0;
    cL0[3][3][14]=0;
    cL0[0][0][15]=2.116568119705976;
    cL0[1][0][15]=-0.0876938905111982;
    cL0[2][0][15]=0.4358867188262607 + 0.3141964285714286*log(rho);
    cL0[3][0][15]=5.153030617371611 + 1.182034027777778*log(rho);
    cL0[0][1][15]=9.16363941630883;
    cL0[1][1][15]=-0.5393053859720526;
    cL0[2][1][15]=-0.02677748869663552;
    cL0[3][1][15]=-3.279161394203301;
    cL0[0][2][15]=9.648572954822955;
    cL0[1][2][15]=0;
    cL0[2][2][15]=0.3863839285714286;
    cL0[3][2][15]=1.364767857142857;
    cL0[0][3][15]=0;
    cL0[1][3][15]=0;
    cL0[2][3][15]=0;
    cL0[3][3][15]=0;
    cL0[0][0][16]=2.036973221615945;
    cL0[1][0][16]=-0.1144988086529505;
    cL0[2][0][16]=0.4398110381929171 + 0.3332142857142857*log(rho);
    cL0[3][0][16]=5.711187912393036 + 1.350634325396825*log(rho);
    cL0[0][1][16]=9.236121718306592;
    cL0[1][1][16]=-0.5399415784832451;
    cL0[2][1][16]=-0.0007755479933059298;
    cL0[3][1][16]=-3.687493053308895;
    cL0[0][2][16]=9.003913488472312;
    cL0[1][2][16]=0;
    cL0[2][2][16]=0.4091517857142857;
    cL0[3][2][16]=1.558327008928571;
    cL0[0][3][16]=0;
    cL0[1][3][16]=0;
    cL0[2][3][16]=0;
    cL0[3][3][16]=0;
    
    
    cL1[0][0][0]=-4.355366341503523;
    cL1[1][0][0]=1.511111111111111;
    cL1[2][0][0]=0.187989417989418;
    cL1[3][0][0]=0.02890025825144873;
    cL1[0][1][0]=-149.;
    cL1[1][1][0]=0;
    cL1[2][1][0]=0;
    cL1[3][1][0]=0;
    cL1[0][2][0]=108.;
    cL1[1][2][0]=0;
    cL1[2][2][0]=0;
    cL1[3][2][0]=0;
    cL1[0][3][0]=0;
    cL1[1][3][0]=0;
    cL1[2][3][0]=0;
    cL1[3][3][0]=0;
    cL1[0][0][1]=-302.4339009754894;
    cL1[1][0][1]=-5.720833333333333;
    cL1[2][0][1]=-0.4946974206349206;
    cL1[3][0][1]=-0.0667088352702192;
    cL1[0][1][1]=539.3333333333333;
    cL1[1][1][1]=0.9;
    cL1[2][1][1]=-0.0525;
    cL1[3][1][1]=-0.01515178571428571;
    cL1[0][2][1]=-324.;
    cL1[1][2][1]=0;
    cL1[2][2][1]=0;
    cL1[3][2][1]=0;
    cL1[0][3][1]=0;
    cL1[1][3][1]=0;
    cL1[2][3][1]=0;
    cL1[3][3][1]=0;
    cL1[0][0][2]=499.539267316993;
    cL1[1][0][2]=4.965277777777778;
    cL1[2][0][2]=0.4698263888888889;
    cL1[3][0][2]=0.08066602300642479;
    cL1[0][1][2]=-621.3333333333333;
    cL1[1][1][2]=-0.9;
    cL1[2][1][2]=-0.04392857142857143;
    cL1[3][1][2]=-0.01172321428571429;
    cL1[0][2][2]=216.;
    cL1[1][2][2]=0;
    cL1[2][2][2]=0;
    cL1[3][2][2]=0;
    cL1[0][3][2]=0;
    cL1[1][3][2]=0;
    cL1[2][3][2]=0;
    cL1[3][3][2]=0;
    cL1[0][0][3]=-468.5196336584965;
    cL1[1][0][3]=-1.409259259259259;
    cL1[2][0][3]=-0.01904789462081129;
    cL1[3][0][3]=0.05914213868417737;
    cL1[0][1][3]=404.8333333333333;
    cL1[1][1][3]=0;
    cL1[2][1][3]=-0.1001785714285714;
    cL1[3][1][3]=-0.04551488095238095;
    cL1[0][2][3]=-108.;
    cL1[1][2][3]=0;
    cL1[2][2][3]=0;
    cL1[3][2][3]=0;
    cL1[0][3][3]=0;
    cL1[1][3][3]=0;
    cL1[2][3][3]=0;
    cL1[3][3][3]=0;
    cL1[0][0][4]=139.8277777777778;
    cL1[1][0][4]=0.1796296296296296;
    cL1[2][0][4]=0.2012858245149912;
    cL1[3][0][4]=0.1639429884994121;
    cL1[0][1][4]=-88.2;
    cL1[1][1][4]=0;
    cL1[2][1][4]=-0.12;
    cL1[3][1][4]=-0.07958630952380952;
    cL1[0][2][4]=0;
    cL1[1][2][4]=0;
    cL1[2][2][4]=0;
    cL1[3][2][4]=0;
    cL1[0][3][4]=0;
    cL1[1][3][4]=0;
    cL1[2][3][4]=0;
    cL1[3][3][4]=0;
    cL1[0][0][5]=38.57462962962963;
    cL1[1][0][5]=0.2396296296296296;
    cL1[2][0][5]=0.2190821759259259;
    cL1[3][0][5]=0.2651579484494415;
    cL1[0][1][5]=-41.92222222222222;
    cL1[1][1][5]=0;
    cL1[2][1][5]=-0.1398214285714286;
    cL1[3][1][5]=-0.1251577380952381;
    cL1[0][2][5]=0;
    cL1[1][2][5]=0;
    cL1[2][2][5]=0;
    cL1[3][2][5]=0;
    cL1[0][3][5]=0;
    cL1[1][3][5]=0;
    cL1[2][3][5]=0;
    cL1[3][3][5]=0;
    cL1[0][0][6]=17.70359788359788;
    cL1[1][0][6]=0.2580246913580247;
    cL1[2][0][6]=0.239533215755438;
    cL1[3][0][6]=0.3949586177861062;
    cL1[0][1][6]=-28.42222222222222;
    cL1[1][1][6]=0;
    cL1[2][1][6]=-0.1596428571428571;
    cL1[3][1][6]=-0.1822291666666667;
    cL1[0][2][6]=0;
    cL1[1][2][6]=0;
    cL1[2][2][6]=0;
    cL1[3][2][6]=0;
    cL1[0][3][6]=0;
    cL1[1][3][6]=0;
    cL1[2][3][6]=0;
    cL1[3][3][6]=0;
    cL1[0][0][7]=10.12319444444444;
    cL1[1][0][7]=0.2652204585537919;
    cL1[2][0][7]=0.2607979366129168;
    cL1[3][0][7]=0.5528219843106996;
    cL1[0][1][7]=-21.98333333333333;
    cL1[1][1][7]=0;
    cL1[2][1][7]=-0.1794642857142857;
    cL1[3][1][7]=-0.2508005952380952;
    cL1[0][2][7]=0;
    cL1[1][2][7]=0;
    cL1[2][2][7]=0;
    cL1[3][2][7]=0;
    cL1[0][3][7]=0;
    cL1[1][3][7]=0;
    cL1[2][3][7]=0;
    cL1[3][3][7]=0;
    cL1[0][0][8]=6.476917989417989;
    cL1[1][0][8]=0.2684744268077601;
    cL1[2][0][8]=0.2824391888699924;
    cL1[3][0][8]=0.7386156768121343;
    cL1[0][1][8]=-18.14444444444444;
    cL1[1][1][8]=0;
    cL1[2][1][8]=-0.1992857142857143;
    cL1[3][1][8]=-0.3308720238095238;
    cL1[0][2][8]=0;
    cL1[1][2][8]=0;
    cL1[2][2][8]=0;
    cL1[3][2][8]=0;
    cL1[0][3][8]=0;
    cL1[1][3][8]=0;
    cL1[2][3][8]=0;
    cL1[3][3][8]=0;
    cL1[0][0][9]=4.398702968841858;
    cL1[1][0][9]=0.2700911228689006;
    cL1[2][0][9]=0.3043020562134319;
    cL1[3][0][9]=0.9522910160806435;
    cL1[0][1][9]=-15.5537037037037;
    cL1[1][1][9]=0;
    cL1[2][1][9]=-0.2191071428571429;
    cL1[3][1][9]=-0.4224434523809524;
    cL1[0][2][9]=0;
    cL1[1][2][9]=0;
    cL1[2][2][9]=0;
    cL1[3][2][9]=0;
    cL1[0][3][9]=0;
    cL1[1][3][9]=0;
    cL1[2][3][9]=0;
    cL1[3][3][9]=0;
    cL1[0][0][10]=3.076199695366362;
    cL1[1][0][10]=0.2709465020576132;
    cL1[2][0][10]=0.3263141320056829;
    cL1[3][0][10]=1.193825101808985;
    cL1[0][1][10]=-13.66565656565657;
    cL1[1][1][10]=0;
    cL1[2][1][10]=-0.2389285714285714;
    cL1[3][1][10]=-0.525514880952381;
    cL1[0][2][10]=0;
    cL1[1][2][10]=0;
    cL1[2][2][10]=0;
    cL1[3][2][10]=0;
    cL1[0][3][10]=0;
    cL1[1][3][10]=0;
    cL1[2][3][10]=0;
    cL1[3][3][10]=0;
    cL1[0][0][11]=2.169693769100335;
    cL1[1][0][11]=0.2714178825289936;
    cL1[2][0][11]=0.348434685813158;
    cL1[3][0][11]=1.463205236104059;
    cL1[0][1][11]=-12.21700336700337;
    cL1[1][1][11]=0;
    cL1[2][1][11]=-0.25875;
    cL1[3][1][11]=-0.6400863095238095;
    cL1[0][2][11]=0;
    cL1[1][2][11]=0;
    cL1[2][2][11]=0;
    cL1[3][2][11]=0;
    cL1[0][3][11]=0;
    cL1[1][3][11]=0;
    cL1[2][3][11]=0;
    cL1[3][3][11]=0;
    cL1[0][0][12]=1.515464241146059;
    cL1[1][0][12]=0.2716835016835017;
    cL1[2][0][12]=0.3706379068462402;
    cL1[3][0][12]=1.760423560589359;
    cL1[0][1][12]=-11.06402486402486;
    cL1[1][1][12]=0;
    cL1[2][1][12]=-0.2785714285714286;
    cL1[3][1][12]=-0.7661577380952381;
    cL1[0][2][12]=0;
    cL1[1][2][12]=0;
    cL1[2][2][12]=0;
    cL1[3][2][12]=0;
    cL1[0][3][12]=0;
    cL1[1][3][12]=0;
    cL1[2][3][12]=0;
    cL1[3][3][12]=0;
    cL1[0][0][13]=1.025593127041678;
    cL1[1][0][13]=0.2718337218337218;
    cL1[2][0][13]=0.3929061374582208;
    cL1[3][0][13]=2.0854748344698;
    cL1[0][1][13]=-10.12097902097902;
    cL1[1][1][13]=0;
    cL1[2][1][13]=-0.2983928571428571;
    cL1[3][1][13]=-0.9037291666666667;
    cL1[0][2][13]=0;
    cL1[1][2][13]=0;
    cL1[2][2][13]=0;
    cL1[3][2][13]=0;
    cL1[0][3][13]=0;
    cL1[1][3][13]=0;
    cL1[2][3][13]=0;
    cL1[3][3][13]=0;
    cL1[0][0][14]=0.6487003728533033;
    cL1[1][0][14]=0.2719169719169719;
    cL1[2][0][14]=0.4152266665213094;
    cL1[3][0][14]=2.438355373889228;
    cL1[0][1][14]=-9.333177933177933;
    cL1[1][1][14]=0;
    cL1[2][1][14]=-0.3182142857142857;
    cL1[3][1][14]=-1.052800595238095;
    cL1[0][2][14]=0;
    cL1[1][2][14]=0;
    cL1[2][2][14]=0;
    cL1[3][2][14]=0;
    cL1[0][3][14]=0;
    cL1[1][3][14]=0;
    cL1[2][3][14]=0;
    cL1[3][3][14]=0;
    cL1[0][0][15]=0.3526901123821911;
    cL1[1][0][15]=0.2719603852937186;
    cL1[2][0][15]=0.4375900110875309;
    cL1[3][0][15]=2.819062487661454;
    cL1[0][1][15]=-8.663888888888889;
    cL1[1][1][15]=0;
    cL1[2][1][15]=-0.3380357142857143;
    cL1[3][1][15]=-1.213372023809524;
    cL1[0][2][15]=0;
    cL1[1][2][15]=0;
    cL1[2][2][15]=0;
    cL1[3][2][15]=0;
    cL1[0][3][15]=0;
    cL1[1][3][15]=0;
    cL1[2][3][15]=0;
    cL1[3][3][15]=0;
    cL1[0][0][16]=0.116431345071051;
    cL1[1][0][16]=0.2719797178130511;
    cL1[2][0][16]=0.459988906866759;
    cL1[3][0][16]=3.227594151123026;
    cL1[0][1][16]=-8.087405731523379;
    cL1[1][1][16]=0;
    cL1[2][1][16]=-0.3578571428571429;
    cL1[3][1][16]=-1.385443452380952;
    cL1[0][2][16]=0;
    cL1[1][2][16]=0;
    cL1[2][2][16]=0;
    cL1[3][2][16]=0;
    cL1[0][3][16]=0;
    cL1[1][3][16]=0;
    cL1[2][3][16]=0;
    cL1[3][3][16]=0;
    
    cL2[0][0][0]=23.75;
    cL2[1][0][0]=0;
    cL2[2][0][0]=0;
    cL2[3][0][0]=0;
    cL2[0][1][0]=-36.;
    cL2[1][1][0]=0;
    cL2[2][1][0]=0;
    cL2[3][1][0]=0;
    cL2[0][2][0]=0;
    cL2[1][2][0]=0;
    cL2[2][2][0]=0;
    cL2[3][2][0]=0;
    cL2[0][3][0]=0;
    cL2[1][3][0]=0;
    cL2[2][3][0]=0;
    cL2[3][3][0]=0;
    cL2[0][0][1]=-114.5833333333333;
    cL2[1][0][1]=0;
    cL2[2][0][1]=0;
    cL2[3][0][1]=0;
    cL2[0][1][1]=108.;
    cL2[1][1][1]=0;
    cL2[2][1][1]=0;
    cL2[3][1][1]=0;
    cL2[0][2][1]=0;
    cL2[1][2][1]=0;
    cL2[2][2][1]=0;
    cL2[3][2][1]=0;
    cL2[0][3][1]=0;
    cL2[1][3][1]=0;
    cL2[2][3][1]=0;
    cL2[3][3][1]=0;
    cL2[0][0][2]=123.8333333333333;
    cL2[1][0][2]=0;
    cL2[2][0][2]=0;
    cL2[3][0][2]=0;
    cL2[0][1][2]=-72.;
    cL2[1][1][2]=0;
    cL2[2][1][2]=0;
    cL2[3][1][2]=0;
    cL2[0][2][2]=0;
    cL2[1][2][2]=0;
    cL2[2][2][2]=0;
    cL2[3][2][2]=0;
    cL2[0][3][2]=0;
    cL2[1][3][2]=0;
    cL2[2][3][2]=0;
    cL2[3][3][2]=0;
    cL2[0][0][3]=-80.58333333333333;
    cL2[1][0][3]=0;
    cL2[2][0][3]=0;
    cL2[3][0][3]=0;
    cL2[0][1][3]=36.;
    cL2[1][1][3]=0;
    cL2[2][1][3]=0;
    cL2[3][1][3]=0;
    cL2[0][2][3]=0;
    cL2[1][2][3]=0;
    cL2[2][2][3]=0;
    cL2[3][2][3]=0;
    cL2[0][3][3]=0;
    cL2[1][3][3]=0;
    cL2[2][3][3]=0;
    cL2[3][3][3]=0;
    cL2[0][0][4]=12.6;
    cL2[1][0][4]=0;
    cL2[2][0][4]=0;
    cL2[3][0][4]=0;
    cL2[0][1][4]=0;
    cL2[1][1][4]=0;
    cL2[2][1][4]=0;
    cL2[3][1][4]=0;
    cL2[0][2][4]=0;
    cL2[1][2][4]=0;
    cL2[2][2][4]=0;
    cL2[3][2][4]=0;
    cL2[0][3][4]=0;
    cL2[1][3][4]=0;
    cL2[2][3][4]=0;
    cL2[3][3][4]=0;
    cL2[0][0][5]=5.755555555555556;
    cL2[1][0][5]=0;
    cL2[2][0][5]=0;
    cL2[3][0][5]=0;
    cL2[0][1][5]=0;
    cL2[1][1][5]=0;
    cL2[2][1][5]=0;
    cL2[3][1][5]=0;
    cL2[0][2][5]=0;
    cL2[1][2][5]=0;
    cL2[2][2][5]=0;
    cL2[3][2][5]=0;
    cL2[0][3][5]=0;
    cL2[1][3][5]=0;
    cL2[2][3][5]=0;
    cL2[3][3][5]=0;
    cL2[0][0][6]=3.826984126984127;
    cL2[1][0][6]=0;
    cL2[2][0][6]=0;
    cL2[3][0][6]=0;
    cL2[0][1][6]=0;
    cL2[1][1][6]=0;
    cL2[2][1][6]=0;
    cL2[3][1][6]=0;
    cL2[0][2][6]=0;
    cL2[1][2][6]=0;
    cL2[2][2][6]=0;
    cL2[3][2][6]=0;
    cL2[0][3][6]=0;
    cL2[1][3][6]=0;
    cL2[2][3][6]=0;
    cL2[3][3][6]=0;
    cL2[0][0][7]=2.94047619047619;
    cL2[1][0][7]=0;
    cL2[2][0][7]=0;
    cL2[3][0][7]=0;
    cL2[0][1][7]=0;
    cL2[1][1][7]=0;
    cL2[2][1][7]=0;
    cL2[3][1][7]=0;
    cL2[0][2][7]=0;
    cL2[1][2][7]=0;
    cL2[2][2][7]=0;
    cL2[3][2][7]=0;
    cL2[0][3][7]=0;
    cL2[1][3][7]=0;
    cL2[2][3][7]=0;
    cL2[3][3][7]=0;
    cL2[0][0][8]=2.425396825396825;
    cL2[1][0][8]=0;
    cL2[2][0][8]=0;
    cL2[3][0][8]=0;
    cL2[0][1][8]=0;
    cL2[1][1][8]=0;
    cL2[2][1][8]=0;
    cL2[3][1][8]=0;
    cL2[0][2][8]=0;
    cL2[1][2][8]=0;
    cL2[2][2][8]=0;
    cL2[3][2][8]=0;
    cL2[0][3][8]=0;
    cL2[1][3][8]=0;
    cL2[2][3][8]=0;
    cL2[3][3][8]=0;
    cL2[0][0][9]=2.083068783068783;
    cL2[1][0][9]=0;
    cL2[2][0][9]=0;
    cL2[3][0][9]=0;
    cL2[0][1][9]=0;
    cL2[1][1][9]=0;
    cL2[2][1][9]=0;
    cL2[3][1][9]=0;
    cL2[0][2][9]=0;
    cL2[1][2][9]=0;
    cL2[2][2][9]=0;
    cL2[3][2][9]=0;
    cL2[0][3][9]=0;
    cL2[1][3][9]=0;
    cL2[2][3][9]=0;
    cL2[3][3][9]=0;
    cL2[0][0][10]=1.835569985569986;
    cL2[1][0][10]=0;
    cL2[2][0][10]=0;
    cL2[3][0][10]=0;
    cL2[0][1][10]=0;
    cL2[1][1][10]=0;
    cL2[2][1][10]=0;
    cL2[3][1][10]=0;
    cL2[0][2][10]=0;
    cL2[1][2][10]=0;
    cL2[2][2][10]=0;
    cL2[3][2][10]=0;
    cL2[0][3][10]=0;
    cL2[1][3][10]=0;
    cL2[2][3][10]=0;
    cL2[3][3][10]=0;
    cL2[0][0][11]=1.646296296296296;
    cL2[1][0][11]=0;
    cL2[2][0][11]=0;
    cL2[3][0][11]=0;
    cL2[0][1][11]=0;
    cL2[1][1][11]=0;
    cL2[2][1][11]=0;
    cL2[3][1][11]=0;
    cL2[0][2][11]=0;
    cL2[1][2][11]=0;
    cL2[2][2][11]=0;
    cL2[3][2][11]=0;
    cL2[0][3][11]=0;
    cL2[1][3][11]=0;
    cL2[2][3][11]=0;
    cL2[3][3][11]=0;
    cL2[0][0][12]=1.495726495726496;
    cL2[1][0][12]=0;
    cL2[2][0][12]=0;
    cL2[3][0][12]=0;
    cL2[0][1][12]=0;
    cL2[1][1][12]=0;
    cL2[2][1][12]=0;
    cL2[3][1][12]=0;
    cL2[0][2][12]=0;
    cL2[1][2][12]=0;
    cL2[2][2][12]=0;
    cL2[3][2][12]=0;
    cL2[0][3][12]=0;
    cL2[1][3][12]=0;
    cL2[2][3][12]=0;
    cL2[3][3][12]=0;
    cL2[0][0][13]=1.372427572427572;
    cL2[1][0][13]=0;
    cL2[2][0][13]=0;
    cL2[3][0][13]=0;
    cL2[0][1][13]=0;
    cL2[1][1][13]=0;
    cL2[2][1][13]=0;
    cL2[3][1][13]=0;
    cL2[0][2][13]=0;
    cL2[1][2][13]=0;
    cL2[2][2][13]=0;
    cL2[3][2][13]=0;
    cL2[0][3][13]=0;
    cL2[1][3][13]=0;
    cL2[2][3][13]=0;
    cL2[3][3][13]=0;
    cL2[0][0][14]=1.269208569208569;
    cL2[1][0][14]=0;
    cL2[2][0][14]=0;
    cL2[3][0][14]=0;
    cL2[0][1][14]=0;
    cL2[1][1][14]=0;
    cL2[2][1][14]=0;
    cL2[3][1][14]=0;
    cL2[0][2][14]=0;
    cL2[1][2][14]=0;
    cL2[2][2][14]=0;
    cL2[3][2][14]=0;
    cL2[0][3][14]=0;
    cL2[1][3][14]=0;
    cL2[2][3][14]=0;
    cL2[3][3][14]=0;
    cL2[0][0][15]=1.181288156288156;
    cL2[1][0][15]=0;
    cL2[2][0][15]=0;
    cL2[3][0][15]=0;
    cL2[0][1][15]=0;
    cL2[1][1][15]=0;
    cL2[2][1][15]=0;
    cL2[3][1][15]=0;
    cL2[0][2][15]=0;
    cL2[1][2][15]=0;
    cL2[2][2][15]=0;
    cL2[3][2][15]=0;
    cL2[0][3][15]=0;
    cL2[1][3][15]=0;
    cL2[2][3][15]=0;
    cL2[3][3][15]=0;
    cL2[0][0][16]=1.105343675931911;
    cL2[1][0][16]=0;
    cL2[2][0][16]=0;
    cL2[3][0][16]=0;
    cL2[0][1][16]=0;
    cL2[1][1][16]=0;
    cL2[2][1][16]=0;
    cL2[3][1][16]=0;
    cL2[0][2][16]=0;
    cL2[1][2][16]=0;
    cL2[2][2][16]=0;
    cL2[3][2][16]=0;
    cL2[0][3][16]=0;
    cL2[1][3][16]=0;
    cL2[2][3][16]=0;
    cL2[3][3][16]=0;
    
    
    compute_Agg2();
    compute_Bgg2();
    compute_Bgg1();
    // if 4*mt^2/mh^2>12 the BFKL coeff cannot be computed from the tabulated values
    // so we switch off matching
    if (_Agg2==0) _matching=false;
    else _matching=true;
}


double HiggsMtExpansion_gg_nnlo_reg::z_times_reg_L0(const double&z){
    
    double res=0.0;
    double max_rho_power=4;
    if (_rho<1e-5) max_rho_power=1;
    for (int k = 0; k < 4; k++){
        for (int m = 0; m < max_rho_power; m++) {
        
            for (int zbarpower = 0; zbarpower < 17; zbarpower++){
                res = res +  cL0[m][k][zbarpower]
                        * pow(_rho,m)
                        * pow(log(1-z),k)
                        * pow(1.-z,zbarpower);
            }

        }
    }
    
    
//    for (int m=0;m<max_rho_power;m++){
//        for (int k=0;k<4;k++){
//            
//            res += c[m][k]* pow(_rho,m) * z * pow(log(1-z),k);
//            
//        }
//    }
    
    return res;
    
    
}

double HiggsMtExpansion_gg_nnlo_reg::z_times_reg_L1(const double&z){
    /*
     double c[4][4];
     
     c[0][0]=-414.1395793757631 - 49.76004030182591/z + 1684.356396028981*z - 3732.787686590407*pow(z,2) + 6148.318422758089*pow(z,3) - 8438.512828999413*pow(z,4) + 10093.66137773521*pow(z,5) - 10093.47537790113*pow(z,6) + 8328.666360694271*pow(z,7) - 5617.08182264796*pow(z,8) + 3061.142797443883*pow(z,9) - 1326.472582962748*pow(z,10) + 446.2589579855482*pow(z,11) - 112.3414133869065*pow(z,12) + 19.91081346711228*pow(z,13) - 2.215591633519006*pow(z,14) + 0.116431345071051*pow(z,15);
     c[1][0]=-34.60893899310566 + 2.729104938271605/z + 183.0847070316515*z - 643.7117283950617*pow(z,2) + 1681.084567901235*pow(z,3) - 3364.202592592593*pow(z,4) + 5287.857530864198*pow(z,5) - 6610.58987654321*pow(z,6) + 6611.019885361552*pow(z,7) - 5289.026499118166*pow(z,8) + 3365.831604938272*pow(z,9) - 1682.945761316872*pow(z,10) + 647.2950280583614*pow(z,11) - 184.9431537598204*pow(z,12) + 36.98888888888889*pow(z,13) - 4.623635870302537*pow(z,14) + 0.2719797178130511*pow(z,15);
     c[2][0]=-46.93758625926665 + 4.402609336114172/z + 259.6322165038353*z - 956.6741398809524*pow(z,2) + 2570.165560056584*pow(z,3) - 5258.313233906526*pow(z,4) + 8402.38879039903*pow(z,5) - 10638.61121713593*pow(z,6) + 10747.20761187878*pow(z,7) - 8668.894465309272*pow(z,8) + 5554.308215513546*pow(z,9) - 2793.127185231949*pow(z,10) + 1079.54730988957*pow(z,11) - 309.7468184783324*pow(z,12) + 62.17774565684536*pow(z,13) - 7.797412520955676*pow(z,14) + 0.459988906866759*pow(z,15);
     c[3][0]=-227.3068446369081 + 18.15772856225618/z + 1394.647294339314*z - 5498.083733408832*pow(z,2) + 15457.5751952482*pow(z,3) - 32676.06499501764*pow(z,4) + 53515.31157393231*pow(z,5) - 69072.57108018279*pow(z,6) + 70862.06102365864*pow(z,7) - 57886.44681811469*pow(z,8) + 37482.93518225987*pow(z,9) - 19018.96822891076*pow(z,10) + 7407.656722362486*pow(z,11) - 2139.676735902266*pow(z,12) + 432.0355908235742*pow(z,13) - 54.46056890562988*pow(z,14) + 3.227594151123026*pow(z,15);
     c[0][1]=1586.205944436827 - 113.5447289638466/z - 7847.284615003733*z + 24880.31323529412*pow(z,2) - 60225.82806229718*pow(z,3) + 114554.8386043695*pow(z,4) - 174255.905840075*pow(z,5) + 212832.3830488139*pow(z,6) - 209155.8502845194*pow(z,7) + 165062.4386043695*pow(z,8) - 103896.992877112*pow(z,9) + 51483.69500167588*pow(z,10) - 19653.10381987294*pow(z,11) + 5579.441013071895*pow(z,12) - 1109.780199049317*pow(z,13) + 138.0623805932629*pow(z,14) - 8.087405731523379*pow(z,15);
     c[1][1]=0.9 - 0.9*z;
     c[2][1]=35.10910714285714 - 3.302678571428571/z - 196.1314285714286*z + 729.045*pow(z,2) - 1969.11*pow(z,3) + 4043.352857142857*pow(z,4) - 6477.746785714286*pow(z,5) + 8217.648214285714*pow(z,6) - 8314.02*pow(z,7) + 6714.295714285714*pow(z,8) - 4306.185*pow(z,9) + 2167.245*pow(z,10) - 838.2214285714286*pow(z,11) + 240.6471428571429*pow(z,12) - 48.33160714285714*pow(z,13) + 6.06375*pow(z,14) - 0.3578571428571429*pow(z,15);
     c[3][1]=98.76257738095238 - 7.950583333333333/z - 603.7223273809524*z + 2374.901416666667*pow(z,2) - 6667.220083333333*pow(z,3) + 14079.36802380952*pow(z,4) - 23040.61283333333*pow(z,5) + 29720.72589880952*pow(z,6) - 30476.00918452381*pow(z,7) + 24885.71026190476*pow(z,8) - 16108.84016666667*pow(z,9) + 8171.414083333333*pow(z,10) - 3181.910845238095*pow(z,11) + 918.8953333333333*pow(z,12) - 185.5065952380952*pow(z,13) + 23.38046726190476*pow(z,14) - 1.385443452380952*pow(z,15);
     c[0][2]=216. - 108./z - 108.*z + 108.*pow(z,2);
     c[1][2]=0;
     c[2][2]=0;
     c[3][2]=0;
     c[0][3]=0;
     c[1][3]=0;
     c[2][3]=0;
     c[3][3]=0;
     */
    double max_rho_power=4;
    if (_rho<1e-5) max_rho_power=1;
    double res=0.0;
    for (int m=0;m<max_rho_power;m++){
        for (int k=0;k<4;k++){
            for (int zbarpower=0;zbarpower<17;zbarpower++){
                res += cL1[m][k][zbarpower]* pow(_rho,m)  * pow(log(1-z),k) * pow(1.-z,zbarpower);
            }
            
        }
    }
    
    return res;
}

double HiggsMtExpansion_gg_nnlo_reg::z_times_reg_L2(const double&z){
    
    /*        double c[4][4];
     
     c[0][0]=-177.6361194905313 - 8.045991100402865/z + 1028.383960429549*z - 3349.580462184874*pow(z,2) + 8179.023056518645*pow(z,3) - 15579.34361014802*pow(z,4) + 23725.86750096309*pow(z,5) - 29003.22614983056*pow(z,6) + 28520.90718350277*pow(z,7) - 22519.76503871945*pow(z,8) + 14180.65533165092*pow(z,9) - 7029.235793890206*pow(z,10) + 2684.046866042454*pow(z,11) - 762.1690624734742*pow(z,12) + 151.6297720253603*pow(z,13) - 18.86678697119874*pow(z,14) + 1.105343675931911*pow(z,15);
     c[1][0]=0;
     c[2][0]=0;
     c[3][0]=0;
     c[0][1]=-72. + 36./z + 36.*z - 36.*pow(z,2);
     c[1][1]=0;
     c[2][1]=0;
     c[3][1]=0;
     c[0][2]=0;
     c[1][2]=0;
     c[2][2]=0;
     c[3][2]=0;
     c[0][3]=0;
     c[1][3]=0;
     c[2][3]=0;
     c[3][3]=0;
     */
    double max_rho_power=4;
    if (_rho<1e-5) max_rho_power=1;
    double res=0.0;
    for (int m=0;m<max_rho_power;m++){
        for (int k=0;k<4;k++){
            for (int zbarpower=0;zbarpower<17;zbarpower++){
                res += cL2[m][k][zbarpower]* pow(_rho,m)  * pow(log(1-z),k) * pow(1.-z,zbarpower);
            }
            
        }
    }
    return res;
}

double HiggsMtExpansion_gg_nnlo_reg::matched_L0(const double& z){
    
    const double expanded = z_times_reg_L0(z);
    
    
    if (_matching){
        return 1./z*(expanded
                     + _Agg2*(-log(z)+MTEXP::log_expanded(z))
                     + pow(1-z,16)*(_Bgg2-expanded));
    }
    else{
        return 1./z*expanded;
    }
    
}

double HiggsMtExpansion_gg_nnlo_reg::matched_L1(const double& z){
    
    const double expanded = z_times_reg_L1(z);
    
    if (_matching){
        return 1./z*(expanded
                     + (-6*_Bgg1)*(-log(z)+MTEXP::log_expanded(z))
                     - pow(1.-z,16)*z_times_reg_L1(1e-14)
                     );
    }
    else
        return 1./z*expanded;
    
}

double HiggsMtExpansion_gg_nnlo_reg::matched_L2(const double& z){
    
    const double expanded = z_times_reg_L2(z);
    
    if (_matching){
        return 1./z*(expanded
                     + 18.*(-log(z)+MTEXP::log_expanded(z))
                     - pow(1.-z,16)*z_times_reg_L2(1e-14)
                     );
    }
    else
        return 1./z*expanded;
}

void HiggsMtExpansion_gg_nnlo_reg::compute_Agg2(){
    const double tau=4./_rho;
    //: data from Robert Harlander et al, http://arxiv.org/abs/0912.2104
    double Agg2Table[23] = {33.0465, 35.9907, 44.2884, 53.1336, 61.8029, 70.1088, 78.0127,
        85.5245, 92.6698, 99.4782, 105.979, 112.199, 118.162, 123.89,
        129.402, 134.716, 139.846, 144.806, 149.609, 154.265, 158.783,
        163.173, 167.442};
    
    _Agg2 =  MTEXP::linear_interpolate(Agg2Table,tau);
}

void  HiggsMtExpansion_gg_nnlo_reg::compute_Bgg1(){
    const double tau=4./_rho;
    //: data from Robert Harlander et al, http://arxiv.org/abs/0912.2104
    double Bgg1Table[23] = {-0.8821, 2.9212, 5.0234, 6.5538, 7.765, 8.7693, 9.6279, 10.3781,
        11.0444, 11.6437, 12.1883, 12.6875, 13.1482, 13.576, 13.9752,
        14.3495, 14.7018, 15.0345, 15.3497, 15.6491, 15.9343, 16.2065, 16.467};
    _Bgg1 = MTEXP::linear_interpolate(Bgg1Table,tau);
    
}



