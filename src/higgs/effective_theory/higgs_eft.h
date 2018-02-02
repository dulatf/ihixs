#ifndef HIGGS_EFT_H
#define HIGGS_EFT_H

#include "constants.h"
#include<vector>
#include "as_series.h"
#include <complex>

using namespace std;



// namespace HEFT
// declaration of all coefficients that enter the Higgs cross section at the EFT level
// together with a couple of utility functions

namespace HEFT {
    
    complex<double> operator*(int i,const complex<double>& c);

    complex<double> operator*(const complex<double>& c,int i);
    
    void check_imaginary_part(const complex<double> c, const char* func_name);
    
    double one_minus_z(const double& z, const double& L);
    double one(const double& z, const double& L);
    
    double intpow(const double& x,int m);
    
    AsSeries n_delta_at_mh();
    AsSeries n_D0_at_mh();
    AsSeries n_D1_at_mh();
    AsSeries n_D2_at_mh();
    AsSeries n_D3_at_mh();
    AsSeries n_D4_at_mh();
    AsSeries n_D5_at_mh();
    
    AsSeries n_delta_log_muf(const double& L);
    AsSeries n_D0_log_muf(const double& L);
    AsSeries n_D1_log_muf(const double& L);
    AsSeries n_D2_log_muf(const double& L);
    AsSeries n_D3_log_muf(const double& L);
    AsSeries n_D4_log_muf(const double& L);

    double nlo_reg(const double& z, const double& L);
    double nlo_reg(const double& z);
    double nlo_reg_L(const double& z);
    double nlo_r_lz0(const double& z, const double& L);
    double nlo_r_lz1(const double& z, const double& L);
    
    double nnlo_reg(const double& z, const double& L);
    double nnlo_r_lz0(const double& z, const double& L);
    double nnlo_r_lz1(const double& z, const double& L);
    double nnlo_r_lz2(const double& z, const double& L);
    double nnlo_r_lz3(const double& z, const double& L);
    double nnlo_r_lz0_const(const double& z, const double& L);
    double nnlo_r_lz0_logz(const double& z, const double& L);
    double nnlo_r_lz0_logz_sq(const double& z, const double& L);
    double nnlo_r_lz0_logz_cube(const double& z, const double& L);
    
    double gg_n2lo_lzbar0_lz0_no_log(const double& z);
    double gg_n2lo_lzbar0_lz1_no_log(const double& z);
    double gg_n2lo_lzbar0_lz2_no_log(const double& z);
    double gg_n2lo_lzbar0_lz3_no_log(const double& z);
    
    double gg_n2lo_lzbar0_lz0_L_1(const double& z);
    double gg_n2lo_lzbar0_lz0_L_2(const double& z);
    double gg_n2lo_lzbar0_lz0_L_3(const double& z);
    double gg_n2lo_lzbar0_lz1_L_1(const double& z);
    double gg_n2lo_lzbar0_lz1_L_2(const double& z);
    double gg_n2lo_lzbar0_lz2_L_1(const double& z);
    
    double gg_n2lo_lzbar1_L_0(const double& z);
    double gg_n2lo_lzbar1_L_1(const double& z);
    double gg_n2lo_lzbar1_L_2(const double& z);
    double gg_n2lo_lzbar2_L_0(const double& z);
    double gg_n2lo_lzbar2_L_1(const double& z);
    double gg_n2lo_lzbar3_L_0(const double& z);

    
    

    double n3lo_reg(const double& z, const double& L);
    

    double n3lo_reg_no_Lf(const double& z,int truncation_order, const double& switch_full_logs);

    //:bernhard November 2017 (triple expansion, practicaly exact for all coeffs)
    class N3LORegEvaluator{
        public:
        N3LORegEvaluator(){
#include "XSW.txt"
#include "XSZ.txt"
#include "XSZb.txt"
        }
        double n3lo_reg_complete_gg(const double& z, unsigned int log_muf_mh_squared_power);
        double n3lo_reg_complete_qg(const double& z, unsigned int log_muf_mh_squared_power);
        double n3lo_reg_complete_qqbar(const double& z, unsigned int log_muf_mh_squared_power);
        double n3lo_reg_complete_qq(const double& z, unsigned int log_muf_mh_squared_power);
        double n3lo_reg_complete_q1q2(const double& z, unsigned int log_muf_mh_squared_power);

        private:
        static constexpr auto NumWTerms = 201;
        static constexpr auto NumZTerms = 199;
        static constexpr auto NumZbTerms = 51;
        //: arrays of coefficients (found in included files)
        double WExp[4][5][4][NumWTerms][6] = {{{{{0.0}}}}};
        double ZExp[4][5][4][NumZTerms][6] = {{{{{0.0}}}}};
        double ZbExp[4][5][4][NumZbTerms][6] = {{{{{0.0}}}}};
        
    };

    //bernhard October 2015 (log12345 exact , log0 series expansion)
    double sgg_N3LO_reg(double zbar);
    
    double n3lo_reg_bernhard(const double& z,int truncation_order);
    
    double n3lo_reg_gg_w012(const double& z);
    double n3lo_reg_gg_w345(const double& z);
    
    
    

    double gg_n3lo_r_lz0(const double& z, const double& L);
    double gg_n3lo_r_lz1(const double& z, const double& L);
    double gg_n3lo_r_lz2(const double& z, const double& L);
    double gg_n3lo_r_lz3(const double& z, const double& L);
    double gg_n3lo_r_lz4(const double& z, const double& L);
    double gg_n3lo_r_lz5(const double& z, const double& L);

    double gg_n3lo_r_lz5_full_series(const double&z, const double& L);
    double gg_n3lo_r_lz4_full_series(const double&z, const double& L);
    double gg_n3lo_r_lz3_full_series(const double&z, const double& L);
    double gg_n3lo_r_lz2_full_series(const double&z, const double& L);
    double gg_n3lo_r_lz1_full_series(const double&z, const double& L);
    double gg_n3lo_r_lz0_full_series(const double&z, const double& L);
    
    double gg_n3lo_r_lz5_series(const double&z, int m);
    double gg_n3lo_r_lz4_series(const double&z, int m);
    double gg_n3lo_r_lz3_series(const double&z, int m);
    double gg_n3lo_r_lz2_series(const double&z, int m);
    double gg_n3lo_r_lz1_series(const double&z, int m);
    double gg_n3lo_r_lz0_series(const double&z, int m);
    double convergenceSeries(const double& z, const double& L);
    double convergenceFull(const double& z, const double& L);

    double gg_n3lo_r_lz0_NS(const double& z, const double& L);
    double gg_n3lo_r_lz1_NS(const double& z, const double& L);
    double gg_n3lo_r_lz2_NS(const double& z, const double& L);
    double gg_n3lo_r_lz3_NS(const double& z, const double& L);
    double gg_n3lo_r_lz4_NS(const double& z, const double& L);
    double gg_n3lo_r_lz5_NS(const double& z, const double& L);

    // exact log(1-z)^2 by Bernhard, October 2015
    double gg_n3lo_r_lz2_exact(const double& z);
    // exact log(1-z)^1 by Bernhard, October 2015
    double gg_n3lo_r_lz1_exact(const double& z);
    
    double nFull_minus_nS(const double& z, const double& L);
    double qg_nFull_minus_nS(const double& z, const double& L);


    double qg_nlo_reg(const double& z, const double& L);
    
    double qg_nlo_r_L0(const double& z);
    double qg_nlo_r_L1(const double& z);
    
    double qg_nnlo_reg(const double& z, const double& L);
    double qg_n3lo_reg(const double& z, const double& L);
    double qg_n3lo_reg_no_Lf(const double& z, const double& L);
    double qg_n3lo_reg_no_Lf(const double& z);
    double qg_n3lo_reg_no_Lf(const double& z,int truncation_order,
                             const double& log_switch);
    //double LEggNNNLOregSte(const double& z, const double& L);
    double LEggN3LOregFalko(const double& z, const double& L);

    double LEggN3LOregFalko_L1(const double& z);
    double LEggN3LOregFalko_L2(const double& z);
    double LEggN3LOregFalko_L3(const double& z);
    double LEggN3LOregBernhard_L3(const double&z);
    
    //double LEqgNNNLOregSte(const double& z, const double& L);
    double LEqgN3LOregFalko(const double& z, const double& L);

    double LEqgN3LOregFalko_L1(const double& z);
    double LEqgN3LOregFalko_L2(const double& z);
    double LEqgN3LOregFalko_L3(const double& z);
    
    
    //double LEqqbNNNLOregSte(const double& z, const double& L);
    double LEqqbN3LOregFalko(const double& z, const double& L);
    
    //double LEqqNNNLOregSte(const double& z, const double& L);
    double LEqqN3LOregFalko(const double& z, const double& L);
    
    //double LEq1q2NNNLOregSte(const double& z, const double& L);
    double LEq1q2N3LOregFalko(const double& z, const double& L);
    
    double qg_nnlo_r_lz0(const double& z, const double& L);
    double qg_nnlo_r_lz1(const double& z, const double& L);
    double qg_nnlo_r_lz2(const double& z, const double& L);
    double qg_nnlo_r_lz3(const double& z, const double& L);
    double qg_nnlo_r_lz0_const(const double& z, const double& L);
    double qg_nnlo_r_lz0_logz(const double& z, const double& L);
    double qg_nnlo_r_lz0_logz_sq(const double& z, const double& L);
    double qg_nnlo_r_lz0_logz_cube(const double& z, const double& L);
    
    
    double qg_n2lo_r_lz0_L0(const double& z);
    double qg_n2lo_r_lz0_L1(const double& z);
    double qg_n2lo_r_lz0_L2(const double& z);
    double qg_n2lo_r_lz1_L0(const double& z);
    double qg_n2lo_r_lz1_L1(const double& z);
    double qg_n2lo_r_lz1_L2(const double& z);
    double qg_n2lo_r_lz2_L0(const double& z);
    double qg_n2lo_r_lz2_L1(const double& z);
    double qg_n2lo_r_lz3_L0(const double& z);

    
    
    double qg_n3lo_r_lz0(const double& z, const double& L);
    double qg_n3lo_r_lz1(const double& z, const double& L);
    double qg_n3lo_r_lz2(const double& z, const double& L);
    double qg_n3lo_r_lz3(const double& z, const double& L);
    double qg_n3lo_r_lz4(const double& z, const double& L);
    double qg_n3lo_r_lz5(const double& z, const double& L);

    double qg_n3lo_r_lz0_NS(const double& z, const double& L);
    double qg_n3lo_r_lz1_NS(const double& z, const double& L);
    double qg_n3lo_r_lz2_NS(const double& z, const double& L);
    double qg_n3lo_r_lz3_NS(const double& z, const double& L);
    double qg_n3lo_r_lz4_NS(const double& z, const double& L);
    double qg_n3lo_r_lz5_NS(const double& z, const double& L);
    
    
    double qg_n3lo_r_lz5_full_series(const double&z, const double& L);
    double qg_n3lo_r_lz4_full_series(const double&z, const double& L);
    double qg_n3lo_r_lz3_full_series(const double&z, const double& L);
    double qg_n3lo_r_lz2_full_series(const double&z, const double& L);
    double qg_n3lo_r_lz1_full_series(const double&z, const double& L);
    double qg_n3lo_r_lz0_full_series(const double&z, const double& L);
    
    double qg_n3lo_r_lzX_series(int logzbpow,const double& z, int m);

    
    
    double qg_nlo_r_lz0(const double& z, const double& L);
    double qg_nlo_r_lz1(const double& z, const double& L);
    
    // bernhards full log(1-z)^2, october 2015
    double qg_n3lo_r_lz2_exact(const double& z);
    double qg_n3lo_r_lz1_exact(const double& z);

    
    
    
    
    double qqb_nlo_reg(const double& z, const double& L);
    double qqbar_nlo_r_L0(const double& z);
    
    double qqb_nnlo_reg(const double& z, const double& L);
    double qqb_n3lo_reg(const double& z, const double& L);
    double qqb_n3lo_reg_no_Lf(const double& z, const double& L);
    double qqb_n3lo_reg_no_Lf(const double& z, int truncation_order);

    double qqb_nlo_r_lz0(const double& z, const double& L);
    
    double qqb_nnlo_r_lz0(const double& z, const double& L);
    double qqb_nnlo_r_lz1(const double& z, const double& L);
    double qqb_nnlo_r_lz2(const double& z, const double& L);

    double qqb_nnlo_r_lz0_const(const double& z, const double& L);
    double qqb_nnlo_r_lz0_logz(const double& z, const double& L);
    double qqb_nnlo_r_lz0_logz_sq(const double& z, const double& L);
    double qqb_nnlo_r_lz0_logz_cube(const double& z, const double& L);
    
    double qqb_nnlo_r_lz0_L0(const double& z);
    double qqb_nnlo_r_lz0_L1(const double& z);
    double qqb_nnlo_r_lz0_L2(const double& z);

    double qqb_nnlo_r_lz1_L0(const double& z);
    double qqb_nnlo_r_lz1_L1(const double& z);

    double qqb_nnlo_r_lz2_L0(const double& z);
    
    double qqbar_n3lo_r_lz5_full_series(const double&z, const double& L);
    double qqbar_n3lo_r_lz4_full_series(const double&z, const double& L);
    double qqbar_n3lo_r_lz3_full_series(const double&z, const double& L);
    double qqbar_n3lo_r_lz2_full_series(const double&z, const double& L);
    double qqbar_n3lo_r_lz1_full_series(const double&z, const double& L);
    double qqbar_n3lo_r_lz0_full_series(const double&z, const double& L);
    
    double qqbar_n3lo_r_lzX_series(int logzbpow,const double& z, int m);
    
    
    double LEqqbN3LOregFalko_L1(const double& z);
    double LEqqbN3LOregFalko_L2(const double& z);
    double LEqqbN3LOregFalko_L3(const double& z);

    
    
    
    double qq_nnlo_reg(const double& z, const double& L);
    double qq_n3lo_reg(const double& z, const double& L);
    double qq_n3lo_reg_no_Lf(const double& z, const double& L);
    double qq_n3lo_reg_no_Lf(const double& z, int truncation_order);
    
    double qq_nnlo_r_lz0(const double& z, const double& L);
    double qq_nnlo_r_lz1(const double& z, const double& L);
    double qq_nnlo_r_lz2(const double& z, const double& L);
    
    double qq_nnlo_r_lz0_const(const double& z, const double& L);
    double qq_nnlo_r_lz0_logz(const double& z, const double& L);
    double qq_nnlo_r_lz0_logz_sq(const double& z, const double& L);
    double qq_nnlo_r_lz0_logz_cube(const double& z, const double& L);
    
    double qq_n2lo_lz0_L0(const double& z);
    double qq_n2lo_lz0_L1(const double& z);
    double qq_n2lo_lz0_L2(const double& z);
    double qq_n2lo_lz1_L0(const double& z);
    double qq_n2lo_lz1_L1(const double& z);
    double qq_n2lo_lz2_L0(const double& z);
    
    
    double qq_n3lo_r_lz5_full_series(const double&z, const double& L);
    double qq_n3lo_r_lz4_full_series(const double&z, const double& L);
    double qq_n3lo_r_lz3_full_series(const double&z, const double& L);
    double qq_n3lo_r_lz2_full_series(const double&z, const double& L);
    double qq_n3lo_r_lz1_full_series(const double&z, const double& L);
    double qq_n3lo_r_lz0_full_series(const double&z, const double& L);
    
    double qq_n3lo_r_lzX_series(int logzbpow,const double& z, int m);
    
    double LEqqN3LOregFalko_L1(const double& z);
    double LEqqN3LOregFalko_L2(const double& z);
    double LEqqN3LOregFalko_L3(const double& z);

    double q1q2_nnlo_reg(const double& z, const double& L);
    double q1q2_n3lo_reg(const double& z, const double& L);
    double q1q2_n3lo_reg_no_Lf(const double& z, const double& L);
    double q1q2_n3lo_reg_no_Lf(const double& z,int truncation_order);
    
    double q1q2_nnlo_r_lz0(const double& z, const double& L);
    double q1q2_nnlo_r_lz1(const double& z, const double& L);
    double q1q2_nnlo_r_lz2(const double& z, const double& L);
    
    double q1q2_nnlo_r_lz0_const(const double& z, const double& L);
    double q1q2_nnlo_r_lz0_logz(const double& z, const double& L);
    double q1q2_nnlo_r_lz0_logz_sq(const double& z, const double& L);
    double q1q2_nnlo_r_lz0_logz_cube(const double& z, const double& L);
    
    
    double q1q2_n2lo_lz0_L0(const double& z);
    double q1q2_n2lo_lz0_L1(const double& z);
    double q1q2_n2lo_lz0_L2(const double& z);
    double q1q2_n2lo_lz1_L0(const double& z);
    double q1q2_n2lo_lz1_L1(const double& z);
    double q1q2_n2lo_lz2_L0(const double& z);
    
    
    
    
    
    
    double q1q2_n3lo_r_lz5_full_series(const double&z, const double& L);
    double q1q2_n3lo_r_lz4_full_series(const double&z, const double& L);
    double q1q2_n3lo_r_lz3_full_series(const double&z, const double& L);
    double q1q2_n3lo_r_lz2_full_series(const double&z, const double& L);
    double q1q2_n3lo_r_lz1_full_series(const double&z, const double& L);
    double q1q2_n3lo_r_lz0_full_series(const double&z, const double& L);
    
    double q1q2_n3lo_r_lzX_series(int logzbpow,const double& z, int m);
    
    double LEq1q2N3LOregFalko_L1(const double& z);
    double LEq1q2N3LOregFalko_L2(const double& z);
    double LEq1q2N3LOregFalko_L3(const double& z);
    
    double LEq1q2N3LOregFalko_L1_new(const double& z);
    double LEq1q2N3LOregFalko_L2_new(const double& z);
    double LEq1q2N3LOregFalko_L3_new(const double& z);
    
    
    double n_LO_delta();
    
    double n_NLO_delta();
    
    double n_NNLO_delta();
    double n_NNLO_delta_L();
    double n_NNLO_delta_L2();
    
    double n_N3LO_delta();
    double n_N3LO_delta_L();
    double n_N3LO_delta_L2();
    double n_N3LO_delta_L3();

    double n_NLO_D0_L();
    
    double n_NLO_D1();
    
    double n_NNLO_D0();
    double n_NNLO_D0_L();
    double n_NNLO_D0_L2();
    
    double n_NNLO_D1();
    double n_NNLO_D1_L();
    double n_NNLO_D1_L2();
    
    double n_NNLO_D2();
    double n_NNLO_D2_L();
    
    double n_NNLO_D3();
    
    
    double n_N3LO_D0();
    double n_N3LO_D0_L();
    double n_N3LO_D0_L2();
    double n_N3LO_D0_L3();
    
    double n_N3LO_D1();
    double n_N3LO_D1_L();
    double n_N3LO_D1_L2();
    double n_N3LO_D1_L3();
    
    double n_N3LO_D2();
    double n_N3LO_D2_L();
    double n_N3LO_D2_L2();
    
    double n_N3LO_D3();
    double n_N3LO_D3_L();
    double n_N3LO_D3_L2();
    double n_N3LO_D2_L3();
    
    double n_N3LO_D4();
    double n_N3LO_D4_L();
    
    double n_N3LO_D5();
}



#endif
