#include "constants.h"
#include<vector>
#include "as_series.h"
#include "higgs_eft.h"
#include "cppchaplin.h"
#include <complex>
using namespace std;


namespace HEFT {

    
    
    
    //: implementation of Bernhard's triple expansion (Nov 2017)
    //: we have here the full gg reg, including logs(mh^2/muf^2)
    double N3LORegEvaluator::n3lo_reg_complete_gg(const double& z, unsigned int log_muf_mh_squared_power) {
        
        //: main forking into three regions of [0,1]
        double res = 0.0;
        if (z<=1./13.) {
            const double x = sqrt(z);
            const double logx = log(z);
            for (int t = 0; t < NumZTerms; t++) {
                for (unsigned int p=0; p < 6; p++) {
                    //: [a_s][initial state][log_muf_mh_squared_power][x power][log(x) power]
                    //: [3: n3lo]  [0: gg] ....
                    res += pow(x,t-2) * pow(logx,p) * ZExp[3][0][log_muf_mh_squared_power][t][p];
                }
            }
        }
        else if (z>1.0/13.0 and z<=0.75) {
            const double x = 0.5-z;
            for (int t = 0; t < NumWTerms; t++) {
                //for (unsigned int p=0; p < 6; p++) {
                    //: [a_s][initial state][log_muf_mh_squared_power][x power][log(x) power]
                    //: [3: n3lo]  [0: gg] ....
                    res += pow(x,t)  * WExp[3][0][log_muf_mh_squared_power][t][0];
                //}
            }
        }
        else if (z>0.75) {
            const double x = 1-z;
            const double logx = log(x);
            for (int t = 0; t < NumZbTerms; t++) {
                for (unsigned int p=0; p < 6; p++) {
                    //: [a_s][initial state][log_muf_mh_squared_power][x power][log(x) power]
                    //: [3: n3lo]  [0: gg] ....
                    res += pow(x,t) * pow(logx,p) * ZbExp[3][0][log_muf_mh_squared_power][t][p];
                }
            }
        }
        else {
            cerr << "src/higgs/effective_theory/eft_gg_reg_n3lo.cpp: n3lo_reg_complete: You shouldn't be here " << endl;
            cerr << "z = " << z << endl;
            exit(1);
        }
        
        return res;
    }
    
    //: n3lo gg reg no log(muf/mh) with truncation
    
    double n3lo_reg_no_Lf(const double& z,int truncation_order,
                          const double& switch_full_logs = 0.0) {
        if (switch_full_logs > 0.0) {
            cout << "Error in eft_gg_reg_n3lo.cpp: n3lo_reg_no_Lf :  "
                 << " switch_full_logs flag > 0.0 but the full ogs are switched off in this version" << endl;
            exit(0);
        }
       
            return
            gg_n3lo_r_lz0_series(z,truncation_order)
            + gg_n3lo_r_lz1_series(z,truncation_order)
            + gg_n3lo_r_lz2_series(z,truncation_order)
            + gg_n3lo_r_lz3_series(z,truncation_order)
            + gg_n3lo_r_lz4_series(z,truncation_order)
            + gg_n3lo_r_lz5_series(z,truncation_order)
            ;

    }
    
    
    
    
    
 
    double gg_n3lo_r_lz5_series(const double&z, int m){
        const int trunc_order = 37;
        if (m>trunc_order){
            cout<<"\nError: the gg_n3lo coefficients are truncated to O("<<trunc_order<<" and you asked for O("<<m<<"))"<<endl;
            exit(EXIT_FAILURE);
        }
        double coeffs[trunc_order+1] = {-216.0000000000000,432.0000000000000,0,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000,216.0000000000000};
        double res=0.;
        for (int i=0;i<m+1;i++){
            res += coeffs[i]*intpow(1.-z,i);
        }
        return res*intpow(log(1.-z),5.);
    }
    
    double gg_n3lo_r_lz4_series(const double&z, int m){
        const int trunc_order = 37;
        if (m>trunc_order){
            cout<<"\nError: the gg_n3lo coefficients are truncated to O("<<trunc_order<<" and you asked for O("<<m<<"))"<<endl;
            exit(EXIT_FAILURE);
        }
        double coeffs[trunc_order+1] = {959.0000000000000,-2519.058641975309,1336.500000000000,-1067.308641975309,-557.0086419753086,-315.1439300411523,-151.3863609641387,-24.78446502057613,79.70461493239271,169.2859910836763,248.0089314663389,318.4022941050719,382.1663851580518,440.5078279744946,494.3186629214407,544.2797745978302,590.9243460303918,634.6787342362914,675.8898754737047,714.8442813458652,751.7815813254195,786.9044119007354,820.3857891195840,852.3747051373648,883.0004446000849,912.3759608645505,940.6005502021229,967.7619939732580,993.9382921824043,1019.199079386859,1043.606790956772,1067.217631153993,1090.082382436958,1112.247086484692,1133.753620765750,1154.640189436428,1174.941743501599,1194.690342202270};
        double res=0.;
        for (int i=0;i<m+1;i++){
            res += coeffs[i]*intpow(1.-z,i);
        }
        return res*intpow(log(1.-z),4.);
    }
    
    double gg_n3lo_r_lz3_series(const double&z, int m){
        const int trunc_order = 37;
        if (m>trunc_order){
            cout<<"\nError: the gg_n3lo coefficients are truncated to O("<<trunc_order<<" and you asked for O("<<m<<"))"<<endl;
            exit(EXIT_FAILURE);
        }
        double coeffs[trunc_order+1] = {1254.029197963407,5659.907653455902,-5554.774691358025,5499.967715616840,1797.583867880214,848.8516387992814,450.9694244155855,250.5501580176683,143.9288974146617,89.88898564862923,68.59808167206168,69.18098562384142,85.08989613823454,112.1026998926647,147.3536780094870,188.8202481691231,235.0309114282467,284.8890591037855,337.5614694960650,392.4049125393021,448.9162383256810,506.6974979659819,565.4310013033452,624.8611255188904,684.7808191841737,745.0214392743344,805.4449964914292,865.9381683280362,926.4076279678638,986.7763650250711,1046.980762437050,1106.968255808748,1166.695445657320,1226.126564878818,1285.232227051339,1343.988398403758,1402.375549134349,1460.377949455786};
        double res=0.;
        for (int i=0;i<m+1;i++){
            res += coeffs[i]*intpow(1.-z,i);
        }
        return res*intpow(log(1.-z),3.);
    }
    
    double gg_n3lo_r_lz2_series(const double&z, int m){
        const int trunc_order = 37;
        if (m>trunc_order){
            cout<<"\nError: the gg_n3lo coefficients are truncated to O("<<trunc_order<<" and you asked for O("<<m<<"))"<<endl;
            exit(EXIT_FAILURE);
        }
        double coeffs[trunc_order+1] = {-11089.32827383219,1520.081363550104,8805.766854133748,-12506.93235145225,-440.3295946232928,1232.087337872434,1646.424924078631,1781.863742199700,1835.655534409699,1861.361191990836,1876.642805654308,1888.264906837485,1899.174892432198,1910.799494164845,1923.879099381453,1938.805311211227,1955.774173340243,1974.864290475342,1996.081030727286,2019.383592868919,2044.702456515913,2071.950959042357,2101.033072265684,2131.848632636985,2164.296832330716,2198.278511902539,2233.697624706056,2270.462129532721,2308.484490194784,2347.681906921503,2387.976366883976,2429.294574877156,2471.567806747992,2514.731715204146,2558.726108526085,2603.494716290826,2648.984951697837,2695.147676906192};
        
        
        double res=0.;
        for (int i=0;i<m+1;i++){
            res += coeffs[i]*intpow(1.-z,i);
        }
        return res*intpow(log(1.-z),2.);
    }
    
    double gg_n3lo_r_lz1_series(const double&z, int m){
        const int trunc_order = 37;
        if (m>trunc_order){
            cout<<"\nError: the gg_n3lo coefficients are truncated to O("<<trunc_order<<" and you asked for O("<<m<<"))"<<endl;
            exit(EXIT_FAILURE);
        }
        double coeffs[trunc_order+1] = {15738.44121219653,-13580.1844621199,1757.56458213909,16078.8844927298,82.94707010075,222.78697455191,947.71318605604,1490.09976660201,1869.96576944240,2145.30184396462,2354.66082361087,2520.81584752866,2657.14370449272,2771.73305778335,2869.69905035360,2954.45051329636,3028.38335966134,3093.26543299187,3150.45537858291,3201.03139853833,3245.87021475392,3285.69781655391,3321.12367859999,3352.66489765910,3380.76389967040,3405.80185425498,3428.10909766921,3447.97339115332,3465.64656378172,3481.34991963036,3495.27868150531,3507.60567220099,3518.48438516134,3528.05156144484,3536.42936428126,3543.72722333125,3550.04340614963,3555.46636306911};
        
        
        double res=0.;
        for (int i=0;i<m+1;i++){
            res += coeffs[i]*intpow(1.-z,i);
        }
        return res*log(1.-z);
    }
    
    double gg_n3lo_r_lz0_series(const double&z, int m){
        const int trunc_order = 37;
        if (m>trunc_order){
            cout<<"\nError: the gg_n3lo coefficients are truncated to O("<<trunc_order<<" and you asked for O("<<m<<"))"<<endl;
            exit(EXIT_FAILURE);
        }
        double coeffs[trunc_order+1]={-5872.588876870346,13334.44000680573,-8488.609005689016,-4281.156843831510,2157.505166572802,907.6324940250037,234.3221065036382,-49.17942750811318,-157.4287152630259,-187.5793092908441,-182.1817424107532,-160.1699994087546,-130.1493242516253,-96.11498685115764,-59.98060165345529,-22.71001572563966,15.17222738673017,53.35066200596565,91.61503328126026,129.8131482558900,167.8289264108613,205.5715357350490,242.9693889491538,279.9663133888160,316.5187631369506,352.5936067174363,388.1662967792010,423.2193347102258,457.7409820105619,491.7241836664278,525.1656737166427,558.0652359932278,590.4250956263226,622.2494197106588,653.5439084299342,684.3154607418020,714.5719013020989,744.3217575804971};
        
        double res=0.;
        for (int i=0;i<m+1;i++){
            res += coeffs[i]*intpow(1.-z,i);
        }
        return res;
    }

    
    
    
    double LEggN3LOregFalko_L1(const double& z){
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2= consts::z2;
        const double z3= consts::z3;
        const double z4= consts::z4;
        complex<double> L1 = (-20*Nc*pow(nf,2)*(-1 + pow(z,2))*(-625 - 4632*z + 3633*pow(z,2) + 1624*pow(z,3) + 72*(2 + 19*z + 25*pow(z,2) - 4*pow(z,3))*z2 - 432*z3 + 864*z*z3 + 216*pow(z,2)*z3) - 20*pow(Nc,3)*pow(nf,2)*(1 + z)*(36*(8 + 28*z + 15*pow(z,2) - 65*pow(z,3) + 12*pow(z,4))*z2 - (-1 + z)*(-2656 + 3394*pow(z,3) + z*(-1373 + 432*z3) + pow(z,2)*(835 + 432*z3))) + 3*nf*(-1 + pow(z,2))*(30*(-168 - 489*z - 738*pow(z,2) + 488*pow(z,3))*z2 + 10*(410 + 144*z3 + 9*z*(115 + 54*z3) - 6*pow(z,2)*(-725 + 369*z3) + pow(z,3)*(-5795 + 1512*z3)) + 270*z*(51 + 71*z)*z4) - 2*pow(Nc,2)*nf*(-1 + pow(z,2))*(90*(-802 + 4201*z - 263*pow(z,2) + 3614*pow(z,3))*z2 + 5*(37351 + 11232*z3 - 36*z*(4526 + 39*z3) - 9*pow(z,2)*(-25787 + 6240*z3) + pow(z,3)*(-111601 + 22464*z3)) + 810*(68 - 23*z + 407*pow(z,2))*z4) + pow(Nc,4)*nf*(90*(3872 - 9535*z + 8792*pow(z,2) - 3933*pow(z,3) - 13000*pow(z,4) + 13452*pow(z,5))*z2 + 10*(pow(z,4)*(397808 - 73494*z3) + z*(369877 - 28674*z3) + 4*(-57731 + 756*z3) - 5*pow(z,3)*(17815 + 1998*z3) + 6*pow(z,2)*(-27814 + 8361*z3) + pow(z,5)*(-280802 + 29592*z3)) + 810*z*(371 - 529*z - 371*pow(z,2) + 529*pow(z,3))*z4) + 2*pow(Nc,5)*(-180*(6116 - 36640*z + 18668*pow(z,2) - 14094*pow(z,3) - 25173*pow(z,4) + 50881*pow(z,5))*z2 - 5*(z*(1446551 - 665712*z3) + pow(z,4)*(1549609 - 604152*z3) + 27*(-46725 + 1312*z3) - 4*pow(z,3)*(5195 + 55728*z3) + 9*pow(z,5)*(-158419 + 93192*z3) + pow(z,2)*(-288034 + 457056*z3)) + 3240*(103 - 163*z + 1086*pow(z,2) + 317*pow(z,3) - 1030*pow(z,4) + 162*pow(z,5))*z4))/(25920.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((-12*nf*pow(1 + z,2)*(1 - 4*z + pow(z,2))*z2 + 18*pow(Nc,2)*nf*pow(1 + z,2)*(8 - 5*z + 8*pow(z,2))*z2 - 12*pow(Nc,4)*nf*(23 + 25*z + 15*pow(z,2) + 19*pow(z,3) + 20*pow(z,4))*z2 + pow(Nc,5)*(6*(1314 + 2614*z + 2415*pow(z,2) + 1984*pow(z,3) + 1023*pow(z,4))*z2 + 72*(72 + 144*z + 216*pow(z,2) + 122*pow(z,3) + 61*pow(z,4))*z3))*chaplin::HPL(-1,z))/(72.*pow(Nc,2)*z*(1 + z)) + ((-192*pow(Nc,3)*pow(nf,2)*(-1 + z)*z*pow(1 + z,2)*z2 + 24*Nc*pow(nf,2)*(2 - 10*z - 9*pow(z,2) + 10*pow(z,3) + 7*pow(z,4))*z2 - nf*(-1 + pow(z,2))*(6*(40 - 69*z - 129*pow(z,2) + 20*pow(z,3))*z2 - 360*z*(1 + z)*z3) + 2*pow(Nc,2)*nf*(-1 + pow(z,2))*(6*(72 - 329*z - 413*pow(z,2) + 8*pow(z,3))*z2 - 36*(-4 + 13*z + 103*pow(z,2))*z3) + pow(Nc,4)*nf*z*(6*(-1489 - 1523*z + 1461*pow(z,2) + 1563*pow(z,3) + 4*pow(z,4))*z2 + 288*(-1 - 25*z + pow(z,2) + 25*pow(z,3))*z3) + 4*pow(Nc,5)*(6*(-1934 + 2498*z + 4379*pow(z,2) - 3895*pow(z,3) - 2500*pow(z,4) + 1430*pow(z,5))*z2 + 144*(-18 + 39*z + 96*pow(z,2) - 42*pow(z,3) - 97*pow(z,4) + 29*pow(z,5))*z3))*chaplin::HPL(1,z))/(288.*pow(Nc,2)*z*(-1 + pow(z,2))) - (2*pow(Nc,3)*(35 + 70*z + 105*pow(z,2) + 58*pow(z,3) + 29*pow(z,4))*z2*chaplin::HPL(-1,-1,z))/(z*(1 + z)) - ((pow(Nc,2)*nf*pow(1 + z,2)*(10 + 29*z - 8*pow(z,2)) + nf*pow(1 + z,2)*(8 + 175*z + 44*pow(z,2)) - 2*pow(Nc,4)*nf*(151 + 95*z - 129*pow(z,2) - 3*pow(z,3) + 110*pow(z,4)) + 2*pow(Nc,5)*(2003 + 3594*z + 1359*pow(z,2) - 500*pow(z,3) - 288*(6 + 12*z + 18*pow(z,2) + 10*pow(z,3) + 5*pow(z,4))*z2))*chaplin::HPL(-1,0,z))/(72.*pow(Nc,2)*z*(1 + z)) - (pow(Nc,3)*(47 + 94*z + 141*pow(z,2) + 88*pow(z,3) + 44*pow(z,4))*z2*chaplin::HPL(-1,1,z))/(z*(1 + z)) + ((-(nf*pow(-1 + z,2)*z*(1 + z)) + 4*pow(Nc,2)*nf*pow(-1 + z,2)*z*(1 + z) - 6*pow(Nc,4)*nf*pow(-1 + z,2)*z*(1 + z) + 4*pow(Nc,5)*(-18 + z - 40*pow(z,2) - pow(z,3) + 37*pow(z,4) + 24*pow(z,5)))*z2*chaplin::HPL(0,-1,z))/(2.*pow(Nc,2)*z*(-1 + pow(z,2))) - ((-(nf*(-1 + z)*z*pow(1 + z,2)) + pow(Nc,4)*nf*z*(-67 + 157*z + 67*pow(z,2) - 157*pow(z,3)) + 2*pow(Nc,2)*nf*(4 + 37*z - 79*pow(z,2) - 37*pow(z,3) + 75*pow(z,4)) + 8*pow(Nc,5)*(69 + 23*z - 53*pow(z,2) - 23*pow(z,3) + 53*pow(z,4) + 13*pow(z,5)))*z2*chaplin::HPL(0,1,z))/(8.*pow(Nc,2)*z*(-1 + pow(z,2))) - (pow(Nc,3)*(47 + 94*z + 141*pow(z,2) + 88*pow(z,3) + 44*pow(z,4))*z2*chaplin::HPL(1,-1,z))/(z*(1 + z)) + ((-8*Nc*pow(nf,2)*z*(35 + 53*z - 4*pow(z,2))*(-1 + pow(z,2)) + 8*pow(Nc,3)*pow(nf,2)*(1 + z)*(4 - 53*z + 3*pow(z,2) + 48*pow(z,3)) + 4*pow(Nc,5)*(24057 + 1899*z - 17271*pow(z,2) + 16927*pow(z,3) - 6665*pow(z,4) - 18705*pow(z,5) + 144*(-10 + 41*z + 127*pow(z,2) - 42*pow(z,3) - 129*pow(z,4) + 16*pow(z,5))*z2) - nf*(-1 + pow(z,2))*(-268 + 24*pow(z,3) + pow(z,2)*(462 - 468*z2) - 3*z*(73 + 156*z2)) - pow(Nc,4)*nf*(1 + z)*(6068 - 2252*pow(z,4) - z*(8461 + 324*z2) - 6*pow(z,3)*(-648 + 1290*z2) + 3*pow(z,2)*(311 + 2688*z2)) + 2*pow(Nc,2)*nf*(-1 + pow(z,2))*(-1598 - 1014*pow(z,3) - z*(2299 + 288*z2) - pow(z,2)*(1777 + 4320*z2)))*chaplin::HPL(1,0,z))/(288.*pow(Nc,2)*z*(-1 + pow(z,2))) - ((23*nf*z*pow(1 + z,2) + 169*pow(Nc,4)*nf*z*pow(1 + z,2) - 4*pow(Nc,2)*nf*(-2 + 48*z + 97*pow(z,2) + 47*pow(z,3)) + 16*pow(Nc,5)*(7 - 180*z - 365*pow(z,2) - 179*pow(z,3) + 5*pow(z,4)))*z2*chaplin::HPL(1,1,z))/(8.*pow(Nc,2)*z*(1 + z)) + ((-2*nf*pow(1 + z,2)*(1 - 4*z + pow(z,2)) + pow(Nc,2)*nf*pow(1 + z,2)*(4 + 5*z + 4*pow(z,2)) - 2*pow(Nc,4)*nf*(7 + 7*z + 3*pow(z,2) + pow(z,3) + 4*pow(z,4)) + pow(Nc,5)*(400 + 838*z + 735*pow(z,2) + 440*pow(z,3) + 209*pow(z,4)))*chaplin::HPL(-1,-1,0,z))/(6.*pow(Nc,2)*z*(1 + z)) - ((-2*nf*pow(1 + z,2)*(1 - 4*z + pow(z,2)) + pow(Nc,4)*nf*(-14 + 21*z + 84*pow(z,2) + 41*pow(z,3) - 4*pow(z,4)) + 26*pow(Nc,2)*nf*(1 + z + pow(z,3) + pow(z,4)) + pow(Nc,5)*(995 + 1962*z + 1545*pow(z,2) + 1260*pow(z,3) + 660*pow(z,4)))*chaplin::HPL(-1,0,0,z))/(12.*pow(Nc,2)*z*(1 + z)) + ((10*nf*(1 + z + pow(z,3) + pow(z,4)) - 2*pow(Nc,2)*nf*(8 + 9*z + 6*pow(z,2) + 9*pow(z,3) + 8*pow(z,4)) + pow(Nc,3)*(457 + 888*z + 840*pow(z,2) + 772*pow(z,3) + 407*pow(z,4)))*chaplin::HPL(-1,1,0,z))/(6.*z*(1 + z)) + ((-(pow(Nc,2)*nf*(-1 + pow(z,2))*(4 + 9*z + 3*pow(z,2) + 16*pow(z,3))) + nf*z*(9 + 6*z - 17*pow(z,2) - 6*pow(z,3) + 8*pow(z,4)) + pow(Nc,4)*nf*z*(41 - 6*z - 61*pow(z,2) + 6*pow(z,3) + 12*pow(z,4)) - 2*pow(Nc,5)*(-105 - 191*z + 24*pow(z,2) - 150*pow(z,3) + 81*pow(z,4) + 319*pow(z,5)))*chaplin::HPL(0,-1,0,z))/(12.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((28*Nc*pow(nf,2)*(-1 + z)*z*pow(1 + z,2) - 28*pow(Nc,3)*pow(nf,2)*(-1 + z)*z*pow(1 + z,2) - pow(Nc,2)*nf*(-1 + pow(z,2))*(80 + 376*z - 275*pow(z,2) + 16*pow(z,3)) + nf*(-8 - 21*z + 29*pow(z,2) + 13*pow(z,3) - 21*pow(z,4) + 8*pow(z,5)) - pow(Nc,4)*nf*(152 + 960*z + 145*pow(z,2) - 964*pow(z,3) - 217*pow(z,4) + 36*pow(z,5)) + 2*pow(Nc,5)*(1678 + 1360*z + 333*pow(z,2) - 1635*pow(z,3) - 1791*pow(z,4) + 363*pow(z,5)))*chaplin::HPL(0,1,0,z))/(24.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((10*nf*(1 + z + pow(z,3) + pow(z,4)) - 2*pow(Nc,2)*nf*(8 + 9*z + 6*pow(z,2) + 9*pow(z,3) + 8*pow(z,4)) + pow(Nc,3)*(457 + 888*z + 840*pow(z,2) + 772*pow(z,3) + 407*pow(z,4)))*chaplin::HPL(1,-1,0,z))/(6.*z*(1 + z)) + ((-96*pow(Nc,3)*pow(nf,2)*(-1 + z)*z*pow(1 + z,2) + 3*nf*(-1 + pow(z,2))*(4 - 21*z - 35*pow(z,2) + 20*pow(z,3)) - 2*pow(Nc,2)*nf*(-1 + pow(z,2))*(-76 + 415*z - 341*pow(z,2) + 32*pow(z,3)) + 4*Nc*pow(nf,2)*(2 - 26*z - 25*pow(z,2) + 26*pow(z,3) + 23*pow(z,4)) + pow(Nc,4)*nf*(12 - 2393*z - 589*pow(z,2) + 2389*pow(z,3) + 561*pow(z,4) + 52*pow(z,5)) - 8*pow(Nc,5)*(-600 - 898*z + 622*pow(z,2) + 194*pow(z,3) - 33*pow(z,4) + 737*pow(z,5)))*chaplin::HPL(1,0,0,z))/(48.*pow(Nc,2)*z*(-1 + pow(z,2))) - ((32*pow(Nc,3)*pow(nf,2)*z*(-1 + pow(z,2)) - 4*Nc*pow(nf,2)*(2 - 12*z + 3*pow(z,2) + 7*pow(z,3)) + 2*pow(Nc,2)*nf*(-1 + z)*(-52 + 329*z + 413*pow(z,2) + 12*pow(z,3)) + nf*(-1 + z)*(40 - 69*z - 129*pow(z,2) + 20*pow(z,3)) + pow(Nc,4)*nf*(64 + 1433*z + 66*pow(z,2) - 1503*pow(z,3) - 68*pow(z,4)) - 4*pow(Nc,5)*(-1477 + 4406*z - 75*pow(z,2) - 3888*pow(z,3) + 1023*pow(z,4)))*chaplin::HPL(1,1,0,z))/(48.*pow(Nc,2)*(-1 + z)*z) - (12*pow(Nc,3)*(4 + 8*z + 12*pow(z,2) + 6*pow(z,3) + 3*pow(z,4))*chaplin::HPL(-1,-1,-1,0,z))/(z*(1 + z)) + (pow(Nc,3)*(61 + 122*z + 183*pow(z,2) + 100*pow(z,3) + 50*pow(z,4))*chaplin::HPL(-1,-1,0,0,z))/(z*(1 + z)) - (2*pow(Nc,3)*(23 + 46*z + 69*pow(z,2) + 40*pow(z,3) + 20*pow(z,4))*chaplin::HPL(-1,-1,1,0,z))/(z*(1 + z)) + (2*pow(Nc,3)*(14 + 28*z + 42*pow(z,2) + 22*pow(z,3) + 11*pow(z,4))*chaplin::HPL(-1,0,-1,0,z))/(z*(1 + z)) - (pow(Nc,3)*(29 + 58*z + 87*pow(z,2) + 48*pow(z,3) + 24*pow(z,4))*chaplin::HPL(-1,0,0,0,z))/(z*(1 + z)) + (pow(Nc,3)*(29 + 58*z + 87*pow(z,2) + 48*pow(z,3) + 24*pow(z,4))*chaplin::HPL(-1,0,1,0,z))/(z*(1 + z)) - (2*pow(Nc,3)*(23 + 46*z + 69*pow(z,2) + 40*pow(z,3) + 20*pow(z,4))*chaplin::HPL(-1,1,-1,0,z))/(z*(1 + z)) + (pow(Nc,3)*(41 + 82*z + 123*pow(z,2) + 72*pow(z,3) + 36*pow(z,4))*chaplin::HPL(-1,1,0,0,z))/(z*(1 + z)) - (24*pow(Nc,3)*pow(1 + z + pow(z,2),2)*chaplin::HPL(-1,1,1,0,z))/(z*(1 + z)) + ((-(nf*pow(-1 + z,2)*z*(1 + z)) - 2*pow(Nc,4)*nf*pow(-1 + z,2)*z*(1 + z) + 4*pow(Nc,5)*(-6 - z - 12*pow(z,2) + pow(z,3) + 10*pow(z,4) + 9*pow(z,5)))*chaplin::HPL(0,-1,-1,0,z))/(pow(Nc,2)*z*(-1 + pow(z,2))) + ((nf*pow(-1 + z,2)*z*(1 + z) - 5*pow(Nc,2)*nf*pow(-1 + z,2)*z*(1 + z) + 7*pow(Nc,4)*nf*pow(-1 + z,2)*z*(1 + z) - 2*pow(Nc,5)*(-30 + 6*z - 70*pow(z,2) - 6*pow(z,3) + 63*pow(z,4) + 41*pow(z,5)))*chaplin::HPL(0,-1,0,0,z))/(2.*pow(Nc,2)*z*(-1 + pow(z,2))) + (2*(nf*pow(-1 + z,2)*z*(1 + z) - pow(Nc,2)*nf*pow(-1 + z,2)*z*(1 + z) + pow(Nc,3)*(-12 + 2*z - 28*pow(z,2) - 2*pow(z,3) + 27*pow(z,4) + 15*pow(z,5)))*chaplin::HPL(0,-1,1,0,z))/(z*(-1 + pow(z,2))) - ((2*pow(Nc,2)*nf*pow(-1 + z,2)*z*(1 + z) + 2*nf*z*(-1 + pow(z,2)) + pow(Nc,4)*nf*z*(-5 + 3*z + 5*pow(z,2) - 3*pow(z,3)) + 4*pow(Nc,5)*(-4 + 7*z - 13*pow(z,2) - 7*pow(z,3) + 12*pow(z,4) + 7*pow(z,5)))*chaplin::HPL(0,0,-1,0,z))/(2.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((nf*(-1 + z)*z*pow(1 + z,2) + pow(Nc,4)*nf*z*(-47 + 17*z + 47*pow(z,2) - 17*pow(z,3)) + 2*pow(Nc,2)*nf*z*(21 - 11*z - 21*pow(z,2) + 11*pow(z,3)) + 8*pow(Nc,5)*(44 + 15*z + 75*pow(z,2) - 15*pow(z,3) - 76*pow(z,4) + 21*pow(z,5)))*chaplin::HPL(0,0,1,0,z))/(8.*pow(Nc,2)*z*(-1 + pow(z,2))) + (2*(nf*pow(-1 + z,2)*z*(1 + z) - pow(Nc,2)*nf*pow(-1 + z,2)*z*(1 + z) + pow(Nc,3)*(-12 + 2*z - 28*pow(z,2) - 2*pow(z,3) + 27*pow(z,4) + 15*pow(z,5)))*chaplin::HPL(0,1,-1,0,z))/(z*(-1 + pow(z,2))) + ((-9*nf*(-1 + z)*z*pow(1 + z,2) + pow(Nc,4)*nf*z*(-99 + 149*z + 99*pow(z,2) - 149*pow(z,3)) + 2*pow(Nc,2)*nf*(-4 + 41*z - 79*pow(z,2) - 41*pow(z,3) + 83*pow(z,4)) + 8*pow(Nc,5)*(43 + 35*z - 73*pow(z,2) - 35*pow(z,3) + 76*pow(z,4) + 10*pow(z,5)))*chaplin::HPL(0,1,0,0,z))/(8.*pow(Nc,2)*z*(-1 + pow(z,2))) - ((-(nf*(-1 + z)*z*pow(1 + z,2)) - 15*pow(Nc,4)*nf*z*(5 - 11*z - 5*pow(z,2) + 11*pow(z,3)) + 2*pow(Nc,2)*nf*(4 + 41*z - 83*pow(z,2) - 41*pow(z,3) + 79*pow(z,4)) + 8*pow(Nc,5)*(57 + 25*z - 81*pow(z,2) - 25*pow(z,3) + 80*pow(z,4) + 28*pow(z,5)))*chaplin::HPL(0,1,1,0,z))/(8.*pow(Nc,2)*z*(-1 + pow(z,2))) - (2*pow(Nc,3)*(23 + 46*z + 69*pow(z,2) + 40*pow(z,3) + 20*pow(z,4))*chaplin::HPL(1,-1,-1,0,z))/(z*(1 + z)) + (pow(Nc,3)*(41 + 82*z + 123*pow(z,2) + 72*pow(z,3) + 36*pow(z,4))*chaplin::HPL(1,-1,0,0,z))/(z*(1 + z)) - (24*pow(Nc,3)*pow(1 + z + pow(z,2),2)*chaplin::HPL(1,-1,1,0,z))/(z*(1 + z)) + (2*(nf*z*(-1 + pow(z,2)) + pow(Nc,2)*nf*(z - pow(z,3)) + pow(Nc,3)*(13 + 10*z + 39*pow(z,2) + 42*pow(z,3) + 14*pow(z,4)))*chaplin::HPL(1,0,-1,0,z))/(z*(1 + z)) + ((-6*nf*(-1 + z)*z*pow(1 + z,2) + pow(Nc,4)*nf*z*(-21 + 103*z + 21*pow(z,2) - 103*pow(z,3)) + 2*pow(Nc,2)*nf*(-2 + 8*z - 55*pow(z,2) - 8*pow(z,3) + 57*pow(z,4)) - 4*pow(Nc,5)*(-10 + 53*z + 223*pow(z,2) - 57*pow(z,3) - 225*pow(z,4) + 24*pow(z,5)))*chaplin::HPL(1,0,0,0,z))/(4.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((3*nf*(-1 + z)*z*pow(1 + z,2) - 24*pow(Nc,2)*nf*z*(3 - 7*z - 3*pow(z,2) + 7*pow(z,3)) + 15*pow(Nc,4)*nf*z*(5 - 11*z - 5*pow(z,2) + 11*pow(z,3)) + 8*pow(Nc,5)*(-8 - 68*z + 136*pow(z,2) + 66*pow(z,3) - 139*pow(z,4) + 19*pow(z,5)))*chaplin::HPL(1,0,1,0,z))/(8.*pow(Nc,2)*z*(-1 + pow(z,2))) - (24*pow(Nc,3)*pow(1 + z + pow(z,2),2)*chaplin::HPL(1,1,-1,0,z))/(z*(1 + z)) + ((19*nf*(-1 + z)*z*pow(1 + z,2) - 6*pow(Nc,2)*nf*z*(-9 - 49*z + 9*pow(z,2) + 49*pow(z,3)) + 5*pow(Nc,4)*nf*z*(-7 - 55*z + 7*pow(z,2) + 55*pow(z,3)) + 8*pow(Nc,5)*(-11 + 191*z + 385*pow(z,2) - 193*pow(z,3) - 387*pow(z,4) + 13*pow(z,5)))*chaplin::HPL(1,1,0,0,z))/(8.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((-23*nf*z*(1 + z) - 169*pow(Nc,4)*nf*z*(1 + z) + 4*pow(Nc,2)*nf*(-2 + 50*z + 47*pow(z,2)) + 16*pow(Nc,5)*(-1 + 193*z + 190*pow(z,2) + pow(z,3)))*chaplin::HPL(1,1,1,0,z))/(8.*pow(Nc,2)*z) + (20*pow(Nc,3)*(-1 + 2*z - pow(z,2) + pow(z,3))*pow(log(1 - z),4))/z + ((8*Nc*pow(nf,2)*(-1 + pow(z,2))*(116 - 236*pow(z,3) + pow(z,2)*(45 - 72*z2) - 12*z*(-31 + 6*z2)) + 4*pow(Nc,3)*pow(nf,2)*(1 + z)*(460 + 27*pow(z,2) + 708*pow(z,4) - 8*z*(-16 + 18*z2) + pow(z,3)*(-1243 + 144*z2)) + 2*pow(Nc,2)*nf*(1 + z)*(8826 - 50679*z + 69510*pow(z,2) - 61517*pow(z,3) + 34184*pow(z,4) - 36*(48 - 64*z + 513*pow(z,2) - 609*pow(z,3) + 112*pow(z,4))*z2 + 1728*z3 - 1728*z*z3 - 11448*pow(z,2)*z3 + 11448*pow(z,3)*z3) - 4*pow(Nc,5)*(-125239 + 275513*z - 101817*pow(z,2) + 37133*pow(z,3) + 218932*pow(z,4) - 320770*pow(z,5) + 36*(494 - 1308*z + 2759*pow(z,2) - 1189*pow(z,3) - 3154*pow(z,4) + 2816*pow(z,5))*z2 + 1728*z3 + 11664*z*z3 + 37152*pow(z,2)*z3 - 11664*pow(z,3)*z3 - 34560*pow(z,4)*z3 + 28944*pow(z,5)*z3) + nf*(-1 + pow(z,2))*(428 - 4748*pow(z,3) + 72*(-4 - 24*z - 63*pow(z,2) + 40*pow(z,3))*z2 + 9*z*(71 - 216*z3) - 3*pow(z,2)*(-389 + 72*z3)) + pow(Nc,4)*nf*(72*(84 - 61*z + 757*pow(z,2) + pow(z,3) - 805*pow(z,4) + 176*pow(z,5))*z2 - (1 + z)*(76208 + 181580*pow(z,4) - 36*pow(z,2)*(-9127 + 540*z3) - z*(233185 + 1512*z3) + pow(z,3)*(-344239 + 20952*z3))))/(1728.*pow(Nc,2)*z*(-1 + pow(z,2))) - (2*pow(Nc,2)*(11*Nc - 2*nf)*pow(1 + z + pow(z,2),2)*chaplin::HPL(-1,0,z))/(3.*z*(1 + z)) + (2*(176*pow(Nc,4) + 11*Nc*nf - 43*pow(Nc,3)*nf - 2*pow(nf,2) + 2*pow(Nc,2)*pow(nf,2))*(1 + z)*chaplin::HPL(1,0,z))/(3.*Nc))*log(z) + ((-16*Nc*pow(nf,2)*(-1 + pow(z,2))*(-8 - 50*z - 41*pow(z,2) + 12*pow(z,3)) + 16*pow(Nc,3)*pow(nf,2)*(1 + z)*(12 + 32*z + 12*pow(z,2) - 70*pow(z,3) + 16*pow(z,4)) + nf*z*(-1 + pow(z,2))*(573 + 1971*z - 304*pow(z,2) - 792*(1 + z)*z2) + 8*pow(Nc,5)*(10292 - 35134*z + 15203*pow(z,2) - 4791*pow(z,3) - 24302*pow(z,4) + 41654*pow(z,5) - 72*(12 - 5*z + 139*pow(z,2) + 5*pow(z,3) - 139*pow(z,4) + 50*pow(z,5))*z2) + 2*pow(Nc,2)*nf*(-1 + pow(z,2))*(-1584 + 7224*pow(z,3) + z*(9293 - 1080*z2) + pow(z,2)*(-3937 + 5256*z2)) + pow(Nc,4)*nf*(-9984 - 25568*pow(z,5) + pow(z,4)*(29863 - 9000*z2) + z*(22743 - 3384*z2) + pow(z,3)*(553 + 3384*z2) + pow(z,2)*(-21511 + 9000*z2)))*pow(log(z),2))/(1152.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((24*Nc*pow(nf,2)*(-1 + z)*z*pow(1 + z,2) - 24*pow(Nc,3)*pow(nf,2)*(-1 + z)*z*pow(1 + z,2) + nf*z*(-3 + 114*z - 68*pow(z,2))*(-1 + pow(z,2)) + 2*pow(Nc,2)*nf*(-1 + pow(z,2))*(-40 + 3*z - 516*pow(z,2) + 60*pow(z,3)) + pow(Nc,4)*nf*(-176 + 35*z - 1482*pow(z,2) + 33*pow(z,3) + 1562*pow(z,4) - 228*pow(z,5)) + 8*pow(Nc,5)*(286 - 331*z + 1158*pow(z,2) - 538*pow(z,3) - 1378*pow(z,4) + 979*pow(z,5)))*pow(log(z),3))/(288.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((11*nf*(-1 + z)*z*pow(1 + z,2) + 11*pow(Nc,4)*nf*z*(11 - 21*z - 11*pow(z,2) + 21*pow(z,3)) - 10*pow(Nc,2)*nf*z*(11 - 25*z - 11*pow(z,2) + 25*pow(z,3)) + 32*pow(Nc,5)*(4 - 11*z + 55*pow(z,2) + 11*pow(z,3) - 54*pow(z,4) + 13*pow(z,5)))*pow(log(z),4))/(384.*pow(Nc,2)*z*(-1 + pow(z,2))) + pow(log(1 - z),3)*(-(pow(Nc,4)*nf*(1956 - 2053*z + 773*pow(z,2) - 1956*pow(z,3)) + 23*nf*(4 + 3*z - 3*pow(z,2) - 4*pow(z,3)) + 192*pow(Nc,2)*nf*(-4 - 3*z + 3*pow(z,2) + 4*pow(z,3)) + 64*pow(Nc,5)*(-671 + 751*z - 641*pow(z,2) + 671*pow(z,3)))/(288.*pow(Nc,2)*z) - ((23*nf*z*(-1 + pow(z,2)) - 192*pow(Nc,2)*nf*z*(-1 + pow(z,2)) + 169*pow(Nc,4)*nf*z*(-1 + pow(z,2)) + 48*pow(Nc,5)*(51 - 38*z + 153*pow(z,2) - 166*pow(z,3) + 51*pow(z,4)))*log(z))/(48.*pow(Nc,2)*(-1 + z)*z)) + pow(log(1 - z),2)*((-16*Nc*pow(nf,2)*(-4 - 7*z + 7*pow(z,3) + 4*pow(z,4)) + 16*pow(Nc,3)*pow(nf,2)*(-8 - 3*z + 4*pow(z,2) + 7*pow(z,3) + 8*pow(z,4)) + 16*pow(Nc,5)*(-9923 + 438*z + 922*pow(z,2) + 484*pow(z,3) + 9923*pow(z,4) + 108*(13 + 16*z + 47*pow(z,2) + 30*pow(z,3) - 15*pow(z,4))*z2) + 2*pow(Nc,2)*nf*(1 + z)*(-2477 + 2909*pow(z,3) + 6*pow(z,2)*(-233 + 288*z2) + 6*z*(161 + 288*z2)) - nf*(1 + z)*(-278 + 404*pow(z,3) + z*(-543 + 414*z2) + pow(z,2)*(417 + 414*z2)) - pow(Nc,4)*nf*(1 + z)*(-13476 + 14214*pow(z,3) + pow(z,2)*(-12317 + 3042*z2) + z*(14155 + 3042*z2)))/(288.*pow(Nc,2)*z*(1 + z)) - (12*pow(Nc,3)*pow(1 + z + pow(z,2),2)*chaplin::HPL(-1,0,z))/(z*(1 + z)) + ((3072*pow(Nc,5) - 23*nf + 192*pow(Nc,2)*nf - 169*pow(Nc,4)*nf)*(1 + z)*chaplin::HPL(1,0,z))/(16.*pow(Nc,2)) + ((32*Nc*pow(nf,2)*z*(-1 + pow(z,2)) - 32*pow(Nc,3)*pow(nf,2)*z*(-1 + pow(z,2)) + nf*(-1 + z)*(32 + 123*z + 75*pow(z,2) - 92*pow(z,3)) + 2*pow(Nc,2)*nf*(-1 + z)*(-288 - 475*z - 214*pow(z,2) + 252*pow(z,3)) + pow(Nc,4)*nf*(-1296 + 261*z - 1782*pow(z,2) + 3229*pow(z,3) - 1164*pow(z,4)) + 8*pow(Nc,5)*(2783 - 4590*z + 7035*pow(z,2) - 8726*pow(z,3) + 4015*pow(z,4)))*log(z))/(96.*pow(Nc,2)*(-1 + z)*z) + ((19*nf*(-1 + z)*z*pow(1 + z,2) - 6*pow(Nc,2)*nf*z*(-9 - 49*z + 9*pow(z,2) + 49*pow(z,3)) + 5*pow(Nc,4)*nf*z*(-7 - 55*z + 7*pow(z,2) + 55*pow(z,3)) + 48*pow(Nc,5)*(23 + 7*z + 89*pow(z,2) - 7*pow(z,3) - 89*pow(z,4) + 27*pow(z,5)))*pow(log(z),2))/(32.*pow(Nc,2)*z*(-1 + pow(z,2)))) + log(1 - z)*((16*Nc*pow(nf,2)*(-1 + pow(z,2))*(-50 + 59*pow(z,3) + 6*pow(z,2)*(12 + 6*z2) + z*(-81 + 36*z2)) - 2*pow(Nc,3)*pow(nf,2)*(-1 + pow(z,2))*(-804 + 876*pow(z,3) + pow(z,2)*(29 + 288*z2) + z*(47 + 288*z2)) + 3*nf*(-1 + pow(z,2))*(-386 - 860*z - 82*pow(z,2) + 1328*pow(z,3) + 6*(16 + 93*z + 105*pow(z,2) - 76*pow(z,3))*z2 + 360*z*z3 + 360*pow(z,2)*z3) - 6*pow(Nc,2)*nf*(-1 + pow(z,2))*(-3269 + 5412*z - 8577*pow(z,2) + 6542*pow(z,3) + 6*(168 + 400*z + 277*pow(z,2) - 172*pow(z,3))*z2 - 288*z3 + 612*z*z3 + 3636*pow(z,2)*z3) + 12*pow(Nc,5)*(-56807 + 62199*z - 3721*pow(z,2) - 4378*pow(z,3) + 60528*pow(z,4) - 57821*pow(z,5) + 12*(1039 - 848*z + 2018*pow(z,2) - 1880*pow(z,3) - 3090*pow(z,4) + 2739*pow(z,5))*z2 + 3384*z3 - 13320*z*z3 + 6840*pow(z,2)*z3 + 360*pow(z,3)*z3 - 19512*pow(z,4)*z3 + 10152*pow(z,5)*z3) + pow(Nc,4)*nf*(-18*(656 + 823*z + 1125*pow(z,2) - 1459*pow(z,3) - 1829*pow(z,4) + 652*pow(z,5))*z2 + 2*(-1 + pow(z,2))*(-45004 + 52996*pow(z,3) + z*(60017 + 432*z3) + pow(z,2)*(-63691 + 10800*z3))))/(864.*pow(Nc,2)*z*(-1 + pow(z,2))) - (pow(Nc,3)*(47 + 94*z + 141*pow(z,2) + 88*pow(z,3) + 44*pow(z,4))*z2*chaplin::HPL(-1,z))/(z*(1 + z)) - ((23*nf*(-1 + z)*z*pow(1 + z,2) + 169*pow(Nc,4)*nf*(-1 + z)*z*pow(1 + z,2) - 4*pow(Nc,2)*nf*(2 - 50*z - 49*pow(z,2) + 50*pow(z,3) + 47*pow(z,4)) + 8*pow(Nc,5)*(-15 + 375*z + 369*pow(z,2) - 373*pow(z,3) - 367*pow(z,4) + 9*pow(z,5)))*z2*chaplin::HPL(1,z))/(8.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((10*nf*(1 + z + pow(z,3) + pow(z,4)) - 2*pow(Nc,2)*nf*(8 + 9*z + 6*pow(z,2) + 9*pow(z,3) + 8*pow(z,4)) + pow(Nc,3)*(457 + 888*z + 840*pow(z,2) + 772*pow(z,3) + 407*pow(z,4)))*chaplin::HPL(-1,0,z))/(6.*z*(1 + z)) - ((-32*Nc*pow(nf,2)*z*(-1 + pow(z,2)) + 32*pow(Nc,3)*pow(nf,2)*z*(-1 + pow(z,2)) + 2*pow(Nc,2)*nf*(-1 + z)*(36 + 262*z + 415*pow(z,2)) + nf*(-1 + z)*(32 - 75*z - 123*pow(z,2) + 28*pow(z,3)) + pow(Nc,4)*nf*(72 + 1417*z + 90*pow(z,2) - 1535*pow(z,3) - 60*pow(z,4)) - 8*pow(Nc,5)*(-685 + 2114*z + 57*pow(z,2) - 2080*pow(z,3) + 583*pow(z,4)))*chaplin::HPL(1,0,z))/(48.*pow(Nc,2)*(-1 + z)*z) - (2*pow(Nc,3)*(23 + 46*z + 69*pow(z,2) + 40*pow(z,3) + 20*pow(z,4))*chaplin::HPL(-1,-1,0,z))/(z*(1 + z)) + (pow(Nc,3)*(41 + 82*z + 123*pow(z,2) + 72*pow(z,3) + 36*pow(z,4))*chaplin::HPL(-1,0,0,z))/(z*(1 + z)) - (24*pow(Nc,3)*pow(1 + z + pow(z,2),2)*chaplin::HPL(-1,1,0,z))/(z*(1 + z)) + (2*(nf*pow(-1 + z,2)*z*(1 + z) - pow(Nc,2)*nf*pow(-1 + z,2)*z*(1 + z) + pow(Nc,3)*(-12 + 2*z - 28*pow(z,2) - 2*pow(z,3) + 27*pow(z,4) + 15*pow(z,5)))*chaplin::HPL(0,-1,0,z))/(z*(-1 + pow(z,2))) - ((-3*nf*(-1 + z)*z*pow(1 + z,2) + 24*pow(Nc,2)*nf*z*(3 - 7*z - 3*pow(z,2) + 7*pow(z,3)) - 15*pow(Nc,4)*nf*z*(5 - 11*z - 5*pow(z,2) + 11*pow(z,3)) + 8*pow(Nc,5)*(54 + 22*z - 90*pow(z,2) - 22*pow(z,3) + 91*pow(z,4) + 27*pow(z,5)))*chaplin::HPL(0,1,0,z))/(8.*pow(Nc,2)*z*(-1 + pow(z,2))) - (24*pow(Nc,3)*pow(1 + z + pow(z,2),2)*chaplin::HPL(1,-1,0,z))/(z*(1 + z)) + ((19*nf*(-1 + z)*z*pow(1 + z,2) + 5*pow(Nc,4)*nf*z*(-7 - 55*z + 7*pow(z,2) + 55*pow(z,3)) + 2*pow(Nc,2)*nf*(4 + 23*z + 145*pow(z,2) - 23*pow(z,3) - 149*pow(z,4)) + 32*pow(Nc,5)*(-1 + 46*z + 98*pow(z,2) - 47*pow(z,3) - 99*pow(z,4) + 5*pow(z,5)))*chaplin::HPL(1,0,0,z))/(8.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((-23*nf*z*(-1 + pow(z,2)) - 169*pow(Nc,4)*nf*z*(-1 + pow(z,2)) + 4*pow(Nc,2)*nf*(2 - 52*z + 3*pow(z,2) + 47*pow(z,3)) + 8*pow(Nc,5)*(3 - 390*z + 9*pow(z,2) + 376*pow(z,3) + 3*pow(z,4)))*chaplin::HPL(1,1,0,z))/(8.*pow(Nc,2)*(-1 + z)*z) + ((8*Nc*pow(nf,2)*(-1 + pow(z,2))*(-12 - 50*z - 41*pow(z,2) + 16*pow(z,3)) - 8*pow(Nc,3)*pow(nf,2)*(1 + z)*(20 + 26*z + 15*pow(z,2) - 77*pow(z,3) + 24*pow(z,4)) - 8*pow(Nc,5)*(14987 - 26842*z + 8910*pow(z,2) - 6544*pow(z,3) - 21779*pow(z,4) + 35504*pow(z,5) - 36*(47 + 15*z + 321*pow(z,2) - 15*pow(z,3) - 323*pow(z,4) + 99*pow(z,5))*z2) + 4*pow(Nc,2)*nf*(-1 + pow(z,2))*(924 - 2563*pow(z,3) + pow(z,2)*(514 - 2160*z2) - 2*z*(1093 + 72*z2)) + nf*(-1 + pow(z,2))*(-136 + 488*pow(z,3) + 6*z*(-95 + 78*z2) + 3*pow(z,2)*(73 + 156*z2)) + pow(Nc,4)*nf*(1 + z)*(12488 + 21492*pow(z,4) + z*(-26606 + 324*z2) - 21*pow(z,2)*(-1721 + 384*z2) + pow(z,3)*(-40571 + 7740*z2)))*log(z))/(288.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((-16*Nc*pow(nf,2)*(-1 + z)*z*pow(1 + z,2) + 16*pow(Nc,3)*pow(nf,2)*(-1 + z)*z*pow(1 + z,2) + nf*(-1 + pow(z,2))*(-4 - 30*z - 54*pow(z,2) + 40*pow(z,3)) - pow(Nc,2)*nf*(-1 + pow(z,2))*(-120 - 136*z - 571*pow(z,2) + 152*pow(z,3)) + pow(Nc,4)*nf*(276 + 109*z + 1038*pow(z,2) - 217*pow(z,3) - 1158*pow(z,4) + 296*pow(z,5)) - 2*pow(Nc,5)*(1980 - 1860*z + 4051*pow(z,2) - 2133*pow(z,3) - 5602*pow(z,4) + 4510*pow(z,5)))*pow(log(z),2))/(48.*pow(Nc,2)*z*(-1 + pow(z,2))) - ((6*nf*(-1 + z)*z*pow(1 + z,2) - 4*pow(Nc,2)*nf*z*(3 - 28*z - 3*pow(z,2) + 28*pow(z,3)) + pow(Nc,4)*nf*z*(21 - 103*z - 21*pow(z,2) + 103*pow(z,3)) + 4*pow(Nc,5)*(40 + 3*z + 273*pow(z,2) - 3*pow(z,3) - 271*pow(z,4) + 74*pow(z,5)))*pow(log(z),3))/(24.*pow(Nc,2)*z*(-1 + pow(z,2))));
        
        check_imaginary_part(L1,__PRETTY_FUNCTION__);
        
        return real(L1);
        
    }
    
    
    double LEggN3LOregFalko_L2(const double& z){
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2= consts::z2;
        const double z3= consts::z3;
        //const double z4= consts::z4;
        complex<double> L2 = (-120*Nc*pow(nf,2)*(-1 + pow(z,2))*(-58 + 76*pow(z,3) + 3*z*(-41 + 18*z2) + 3*pow(z,2)*(35 + 18*z2)) + 120*pow(Nc,3)*pow(nf,2)*(-1 + z)*(1 + z)*(-115 + 133*pow(z,3) + z*(-32 + 54*z2) + pow(z,2)*(34 + 54*z2)) - 90*pow(Nc,2)*nf*(-1 + pow(z,2))*(1613 - 4234*z - 3089*pow(z,3) + 12*(-36 - 81*z - 57*pow(z,2) + 28*pow(z,3))*z2 + pow(z,2)*(5620 - 1728*z3)) + 15*nf*(-1 + pow(z,2))*(214 - 1357*pow(z,3) + 36*(-4 - 15*z - 6*pow(z,2) + 16*pow(z,3))*z2 + pow(z,2)*(831 - 216*z3) - 24*z*(-13 + 9*z3)) - 30*pow(Nc,5)*(-137875 + 72*(418 - 193*z + 881*pow(z,2) - 852*pow(z,3) - 1310*pow(z,4) + 1056*pow(z,5))*z2 + pow(z,4)*(148235 - 51840*z3) + z*(154387 - 25056*z3) + 6912*z3 - 192*pow(z,3)*(86 + 27*z3) + 40*pow(z,2)*(-259 + 540*z3) + pow(z,5)*(-137875 + 19008*z3)) + 15*pow(Nc,4)*nf*(-1 + z)*(36*(-116 - 367*z - 623*pow(z,2) - 276*pow(z,3) + 104*pow(z,4))*z2 - (1 + z)*(-41536 + 49249*pow(z,3) + z*(58060 - 216*z3) + pow(z,2)*(-61193 + 10152*z3))))/(25920.*pow(Nc,2)*z*(-1 + pow(z,2))) + (8*pow(Nc,3)*pow(1 + z + pow(z,2),2)*z2*chaplin::HPL(-1,z))/(z*(1 + z)) + ((1296*pow(Nc,4)*nf*(-1 + z)*z*pow(1 + z,2)*z2 + 12*pow(Nc,2)*nf*(-120*z - 120*pow(z,2))*(-1 + pow(z,2))*z2 - 6*nf*(-24*z - 24*pow(z,2))*(-1 + pow(z,2))*z2 + 576*pow(Nc,5)*(-1 + 39*z + 39*pow(z,2) - 39*pow(z,3) - 39*pow(z,4) + pow(z,5))*z2)*chaplin::HPL(1,z))/(288.*pow(Nc,2)*z*(-1 + pow(z,2))) - ((6*pow(Nc,2)*nf*pow(1 + z,2)*(4 - 7*z + 4*pow(z,2)) - 6*pow(Nc,4)*nf*(8 + 9*z + 6*pow(z,2) + 9*pow(z,3) + 8*pow(z,4)) + 12*pow(Nc,5)*(99 + 182*z + 177*pow(z,2) + 182*pow(z,3) + 99*pow(z,4)))*chaplin::HPL(-1,0,z))/(72.*pow(Nc,2)*z*(1 + z)) + ((-72*Nc*pow(nf,2)*z*(1 + z)*(-1 + pow(z,2)) + 6*nf*(-1 + pow(z,2))*(8 - 6*z - 15*pow(z,2) + 4*pow(z,3)) + 12*pow(Nc,2)*nf*(-1 + pow(z,2))*(12 + 57*z + 75*pow(z,2) + 4*pow(z,3)) + 8*pow(Nc,3)*pow(nf,2)*(1 + z)*(-9*z + 9*pow(z,3)) - 480*pow(Nc,5)*pow(1 + z,2)*(-11 + 56*z - 56*pow(z,2) + 11*pow(z,3)) - 6*pow(Nc,4)*nf*(1 + z)*(-32 - 300*z - 27*pow(z,2) + 347*pow(z,3) + 12*pow(z,4)))*chaplin::HPL(1,0,z))/(288.*pow(Nc,2)*z*(-1 + pow(z,2))) + (8*pow(Nc,3)*pow(1 + z + pow(z,2),2)*chaplin::HPL(-1,-1,0,z))/(z*(1 + z)) - (6*pow(Nc,3)*pow(1 + z + pow(z,2),2)*chaplin::HPL(-1,0,0,z))/(z*(1 + z)) + (4*pow(Nc,3)*pow(1 + z + pow(z,2),2)*chaplin::HPL(-1,1,0,z))/(z*(1 + z)) + ((6*pow(Nc,4)*nf*pow(-1 + z,2)*z*(1 + z) - pow(Nc,2)*nf*(-1 + pow(z,2))*(-6*z + 6*pow(z,2)) - 48*pow(Nc,5)*(-1 + z - 3*pow(z,2) - pow(z,3) + 3*pow(z,4) + pow(z,5)))*chaplin::HPL(0,-1,0,z))/(12.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((-(pow(Nc,2)*nf*(66*z - 114*pow(z,2))*(-1 + pow(z,2))) - pow(Nc,4)*nf*(66*z - 114*pow(z,2) - 66*pow(z,3) + 114*pow(z,4)) + 96*pow(Nc,5)*(3 + 2*z - 5*pow(z,2) - 2*pow(z,3) + 5*pow(z,4) + 2*pow(z,5)))*chaplin::HPL(0,1,0,z))/(24.*pow(Nc,2)*z*(-1 + pow(z,2))) + (4*pow(Nc,3)*pow(1 + z + pow(z,2),2)*chaplin::HPL(1,-1,0,z))/(z*(1 + z)) + ((-2*pow(Nc,2)*nf*(-6*z - 186*pow(z,2))*(-1 + pow(z,2)) + 3*nf*(-6*z - 6*pow(z,2))*(-1 + pow(z,2)) + pow(Nc,4)*nf*(-6*z + 354*pow(z,2) + 6*pow(z,3) - 354*pow(z,4)) - 96*pow(Nc,5)*(-1 + 15*z + 39*pow(z,2) - 15*pow(z,3) - 39*pow(z,4) + pow(z,5)))*chaplin::HPL(1,0,0,z))/(48.*pow(Nc,2)*z*(-1 + pow(z,2))) - ((nf*(-1 + z)*(-24*z - 24*pow(z,2)) + 2*pow(Nc,2)*nf*(-1 + z)*(120*z + 120*pow(z,2)) - 4*pow(Nc,5)*(960*z - 960*pow(z,3)) + pow(Nc,4)*nf*(216*z - 216*pow(z,3)))*chaplin::HPL(1,1,0,z))/(48.*pow(Nc,2)*(-1 + z)*z) - (16*pow(Nc,3)*(-1 + 2*z - pow(z,2) + pow(z,3))*pow(log(1 - z),3))/z + ((8*Nc*pow(nf,2)*(-1 + pow(z,2))*(24 + 93*z + 75*pow(z,2) - 36*pow(z,3)) + 24*pow(Nc,3)*pow(nf,2)*(1 + z)*(12 + 19*z + 6*pow(z,2) - 49*pow(z,3) + 16*pow(z,4)) + nf*(-1 + pow(z,2))*(297*z - 432*pow(z,2) - 144*pow(z,3) + 72*(-9*z - 9*pow(z,2))*z2) - 24*pow(Nc,5)*(-5238 + 11857*z - 3557*pow(z,2) + 2889*pow(z,3) + 8138*pow(z,4) - 15403*pow(z,5) + 576*(1 + 8*pow(z,2) - 8*pow(z,4) + 2*pow(z,5))*z2) - 3*pow(Nc,4)*nf*(1 + z)*(4992 + 9344*pow(z,4) + pow(z,2)*(16833 - 4032*z2) + 3*z*(-4003 + 216*z2) + 12*pow(z,3)*(-1514 + 282*z2)) + 6*pow(Nc,2)*nf*(-1 + z)*(1 + z)*(-792 + 2264*pow(z,3) + z*(3083 - 216*z2) + 4*pow(z,2)*(-157 + 450*z2)))*log(z))/(1728.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((-6*nf*pow(z,2)*(-33 + 32*z)*(-1 + pow(z,2)) - 16*Nc*pow(nf,2)*(-9*z - 9*pow(z,2))*(-1 + pow(z,2)) + 16*pow(Nc,3)*pow(nf,2)*(1 + z)*(9*z - 9*pow(z,3)) + 12*pow(Nc,2)*nf*(-1 + pow(z,2))*(-40 - 273*pow(z,2) + 32*pow(z,3)) - 6*pow(Nc,4)*nf*(176 + 64*z + 913*pow(z,2) - 96*pow(z,3) - 993*pow(z,4) + 160*pow(z,5)) + 96*pow(Nc,5)*(143 - 104*z + 391*pow(z,2) - 204*pow(z,3) - 501*pow(z,4) + 352*pow(z,5)))*pow(log(z),2))/(1152.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((9*nf*z*(1 + z)*(-1 + pow(z,2)) + 2*pow(Nc,2)*nf*(57*z - 129*pow(z,2))*(-1 + pow(z,2)) + pow(Nc,4)*nf*(123*z - 249*pow(z,2) - 123*pow(z,3) + 249*pow(z,4)) + 192*pow(Nc,5)*(1 - 2*z + 11*pow(z,2) + 2*pow(z,3) - 11*pow(z,4) + 2*pow(z,5)))*pow(log(z),3))/(288.*pow(Nc,2)*z*(-1 + pow(z,2))) + pow(log(1 - z),2)*((-12*nf*(1 + z)*(-4 - 3*z + 3*pow(z,2) + 4*pow(z,3)) + 120*pow(Nc,2)*nf*(1 + z)*(-4 - 3*z + 3*pow(z,2) + 4*pow(z,3)) - 108*pow(Nc,4)*nf*(1 + z)*(-12 + 13*z - 5*pow(z,2) + 12*pow(z,3)) + 144*pow(Nc,5)*(-187 + 21*z + 33*pow(z,2) + 12*pow(z,3) + 187*pow(z,4)))/(288.*pow(Nc,2)*z*(1 + z)) + ((2*pow(Nc,2)*nf*(-1 + z)*(-120*z - 120*pow(z,2)) + nf*(-1 + z)*(24*z + 24*pow(z,2)) + pow(Nc,4)*nf*(-216*z + 216*pow(z,3)) + 384*pow(Nc,5)*(7 - 4*z + 21*pow(z,2) - 24*pow(z,3) + 7*pow(z,4)))*log(z))/(96.*pow(Nc,2)*(-1 + z)*z)) + log(1 - z)*((36*Nc*pow(nf,2)*(-1 + pow(z,2))*(-4 - 3*z + 3*pow(z,2) + 4*pow(z,3)) - 12*pow(Nc,3)*pow(nf,2)*(-1 + pow(z,2))*(-20 + 7*z + pow(z,2) + 20*pow(z,3)) + 12*pow(Nc,5)*(-1 + z)*(17149 - 500*z - 1314*pow(z,2) - 814*pow(z,3) - 17149*pow(z,4) + 144*(-17 - 20*z - 59*pow(z,2) - 38*pow(z,3) + 19*pow(z,4))*z2) + 6*nf*(-1 + pow(z,2))*(-34 + 52*pow(z,3) + 3*z*(-29 + 24*z2) + 3*pow(z,2)*(23 + 24*z2)) - 18*pow(Nc,2)*nf*(-1 + pow(z,2))*(-392 + 440*pow(z,3) + pow(z,2)*(-319 + 240*z2) + z*(271 + 240*z2)) + 12*pow(Nc,4)*nf*(-1 + pow(z,2))*(-1517 + 1580*pow(z,3) + 2*z*(824 + 162*z2) + pow(z,2)*(-1463 + 324*z2)))/(864.*pow(Nc,2)*z*(-1 + pow(z,2))) + (4*pow(Nc,3)*pow(1 + z + pow(z,2),2)*chaplin::HPL(-1,0,z))/(z*(1 + z)) - ((nf*(-1 + z)*(-24*z - 24*pow(z,2)) + 2*pow(Nc,2)*nf*(-1 + z)*(120*z + 120*pow(z,2)) - 8*pow(Nc,5)*(480*z - 480*pow(z,3)) + pow(Nc,4)*nf*(216*z - 216*pow(z,3)))*chaplin::HPL(1,0,z))/(48.*pow(Nc,2)*(-1 + z)*z) + ((8*Nc*pow(nf,2)*(-9*z - 9*pow(z,2))*(-1 + pow(z,2)) - 8*pow(Nc,3)*pow(nf,2)*(1 + z)*(9*z - 9*pow(z,3)) + 6*nf*(-1 + pow(z,2))*(-4 - 15*z - 6*pow(z,2) + 16*pow(z,3)) - 12*pow(Nc,2)*nf*(-1 + pow(z,2))*(-60 - 93*z - 39*pow(z,2) + 44*pow(z,3)) + 6*pow(Nc,4)*nf*(1 + z)*(276 - 41*z + 381*pow(z,2) - 688*pow(z,3) + 232*pow(z,4)) - 96*pow(Nc,5)*(275 - 134*z + 254*pow(z,2) - 196*pow(z,3) - 474*pow(z,4) + 385*pow(z,5)))*log(z))/(288.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((-(pow(Nc,2)*nf*(-6*z - 186*pow(z,2))*(-1 + pow(z,2))) + nf*(-9*z - 9*pow(z,2))*(-1 + pow(z,2)) + pow(Nc,4)*nf*(-3*z + 177*pow(z,2) + 3*pow(z,3) - 177*pow(z,4)) - 96*pow(Nc,5)*(5 + 2*z + 25*pow(z,2) - 2*pow(z,3) - 25*pow(z,4) + 6*pow(z,5)))*pow(log(z),2))/(48.*pow(Nc,2)*z*(-1 + pow(z,2))));
        
        check_imaginary_part(L2,__PRETTY_FUNCTION__);
        
        return real(L2);
    }
    
    
    double LEggN3LOregFalko_L3(const double& z){
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2= consts::z2;

        complex<double> L3=(-240*Nc*pow(nf,2)*(-1 + pow(z,2))*(-4 - 3*z + 3*pow(z,2) + 4*pow(z,3)) + 240*pow(Nc,3)*pow(nf,2)*(1 + z)*(6 - 7*z - 5*pow(z,3) + 6*pow(z,4)) - 45*nf*z*(-1 + pow(z,2))*(23*(-1 + z) + 24*(1 + z)*z2) - 120*pow(Nc,5)*(-1 + pow(z,2))*(6085 - 6326*z + 6205*pow(z,2) - 6085*pow(z,3) + 864*(-1 - 3*pow(z,2) + pow(z,3))*z2) + 90*pow(Nc,2)*nf*(-1 + pow(z,2))*(-352 + 352*pow(z,3) + pow(z,2)*(-493 + 144*z2) + z*(493 + 144*z2)) - 15*pow(Nc,4)*nf*(-1 + pow(z,2))*(-4576 + 4576*pow(z,3) + pow(z,2)*(-5107 + 792*z2) + z*(5459 + 792*z2)))/(25920.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((2304*pow(Nc,5)*(-1 + z)*z*pow(1 + z,2) - 12*nf*z*(1 + z)*(-1 + pow(z,2)) + 144*pow(Nc,2)*nf*z*(1 + z)*(-1 + pow(z,2)) - 132*pow(Nc,4)*nf*z*(1 + z)*(-1 + pow(z,2)))*chaplin::HPL(1,0,z))/(288.*pow(Nc,2)*z*(-1 + pow(z,2))) + (4*pow(Nc,3)*(-1 + z + pow(z,2) + pow(z,4))*pow(log(1 - z),2))/(z*(1 + z)) + ((96*Nc*pow(nf,2)*z*(1 + z)*(-1 + pow(z,2)) - 96*pow(Nc,3)*pow(nf,2)*z*(1 + z)*(-1 + pow(z,2)) + nf*(-1 + pow(z,2))*(-36*pow(z,2) - 48*pow(z,3)) + 48*pow(Nc,2)*nf*(1 + z)*(10 + z - 9*pow(z,2) - 4*pow(z,3) + 2*pow(z,4)) - 12*pow(Nc,4)*nf*(88 + 92*z + 115*pow(z,2) - 96*pow(z,3) - 155*pow(z,4) + 52*pow(z,5)) + 288*pow(Nc,5)*(55 - z + 49*pow(z,2) - 43*pow(z,3) - 93*pow(z,4) + 55*pow(z,5)))*log(z))/(1728.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((12*nf*z*(1 + z)*(-1 + pow(z,2)) - 72*pow(Nc,2)*nf*z*(-3 + 7*z)*(-1 + pow(z,2)) + 12*pow(Nc,4)*nf*z*(19 - 41*z - 19*pow(z,2) + 41*pow(z,3)) + 768*pow(Nc,5)*(1 - z + 7*pow(z,2) + pow(z,3) - 7*pow(z,4) + pow(z,5)))*pow(log(z),2))/(1152.*pow(Nc,2)*z*(-1 + pow(z,2))) + log(1 - z)*((-72*pow(Nc,2)*nf*(-1 + pow(z,2))*(-4 - 3*z + 3*pow(z,2) + 4*pow(z,3)) + 3*nf*(-1 + pow(z,2))*(-8 - 6*z + 6*pow(z,2) + 8*pow(z,3)) - 288*pow(Nc,5)*(55 - 58*z - 8*pow(z,2) + 3*pow(z,3) - 47*pow(z,4) + 55*pow(z,5)) + 6*pow(Nc,4)*nf*(140 - 159*z - 77*pow(z,2) + 19*pow(z,3) - 63*pow(z,4) + 140*pow(z,5)))/(864.*pow(Nc,2)*z*(-1 + pow(z,2))) + ((-12*nf*z*(1 + z)*(-1 + pow(z,2)) + 144*pow(Nc,2)*nf*z*(1 + z)*(-1 + pow(z,2)) - 132*pow(Nc,4)*nf*z*(1 + z)*(-1 + pow(z,2)) - 1152*pow(Nc,5)*(1 + z + 3*pow(z,2) - pow(z,3) - 3*pow(z,4) + pow(z,5)))*log(z))/(288.*pow(Nc,2)*z*(-1 + pow(z,2))));
        check_imaginary_part(L3,__PRETTY_FUNCTION__);
        
        return real(L3);
    }
    
    
    
    double LEggN3LOregFalko(const double& z, const double& L)
    {
        if (abs(L)<1e-10) return 0.0;
        
        return   LEggN3LOregFalko_L3(z)*pow(L,3)
        +LEggN3LOregFalko_L2(z)*pow(L,2)
        +LEggN3LOregFalko_L1(z)*L;
    }



}
