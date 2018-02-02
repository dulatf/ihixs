#ifndef INCLUSIVE_HIGGS_EFT_H
#define INCLUSIVE_HIGGS_EFT_H

#include "inclusive_process.h"


class InclusiveHiggsEFT: public InclusiveProcess{
public:
    InclusiveHiggsEFT(const UserInterface& UI):InclusiveProcess(UI){_name = "Higgs Production in gluon fusion: Effective Theory (mt-> infinity)";}
    void SetUpContributions();
    void EvaluatePostVegasDependingOnProcess();
    void AddGG();
    void AddQG();
    void AddQQBAR();
    void AddQQ();
    void AddQ1Q2();
    
    double RescalingCoefficient(){return _input._top_only_LO_coefficient;}
    // EW corrections as a series in a_s
    AsSeries EwCorrections();
    
    AsSeries WC(){return _channels[0]->GiveTermPtr(0)->WC();}
    AsSeries Eta();
    
    ResultPair R_LO_x_EFT(int order);
    ResultPair R_LO_x_EFT_N3LO() {
        AsSeries eft_res = QCDCorrections();
        ResultPair eft_n3lo = eft_res.term_of_order(2)+eft_res.term_of_order(3)+eft_res.term_of_order(4)+eft_res.term_of_order(5);
        return RescalingCoefficient() * eft_n3lo ;
    }
    
    ResultPair EFT_N3LO() {
        AsSeries eft_res = QCDCorrections();
        ResultPair eft_n3lo = eft_res.term_of_order(2)+eft_res.term_of_order(3)+eft_res.term_of_order(4)+eft_res.term_of_order(5);
        return eft_n3lo ;
    }
    
};

#include "higgs_eft.h"



///
class HiggsEFTGGDelta: public LuminosityIntegralAbsolutelyGenericDelta
{
public:
    HiggsEFTGGDelta(){set_dimensions(2); // for cuhre
                    _channel = "gg";
                    _name="gg_delta_lo";
        
}
    
    
};

class HiggsEFTGGPlus: public LuminosityIntegralAbsolutelyGenericPlus
{
public:
    HiggsEFTGGPlus(int m){set_dimensions(2);_channel = "gg";_name="gg_D"+to_string(m);
        _log_power=m;}
};



class HiggsEFTGGNLOReg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsEFTGGNLOReg(){set_dimensions(2);_channel = "gg";_name="gg_reg_nlo";
        }
    void compute_matrix_element_Lf_coefficients(const double* xx){
        _me_lf_coeffs.clear();
        _me_lf_coeffs.push_back(HEFT::nlo_reg(xx[1]));
        _me_lf_coeffs.push_back(HEFT::nlo_reg_L(xx[1]));
    }
    
};

class HiggsEFTGGNNLOReg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsEFTGGNNLOReg(){set_dimensions(2);_channel = "gg";_name="gg_reg_n2lo";
        set_rel_accuracy_multiplier(10.);
        }
    
    void compute_matrix_element_Lf_coefficients(const double* xx){
        const double z=xx[1];
        _me_lf_coeffs.clear();
        // L^0
        _me_lf_coeffs.push_back(HEFT::gg_n2lo_lzbar0_lz0_no_log(z)
                                +HEFT::gg_n2lo_lzbar0_lz1_no_log(z)
                                +HEFT::gg_n2lo_lzbar0_lz2_no_log(z)
                                +HEFT::gg_n2lo_lzbar0_lz3_no_log(z)
                                +HEFT::gg_n2lo_lzbar1_L_0(z)
                                +HEFT::gg_n2lo_lzbar2_L_0(z)
                                +HEFT::gg_n2lo_lzbar3_L_0(z));
        // L^1
        _me_lf_coeffs.push_back(+HEFT::gg_n2lo_lzbar0_lz0_L_1(z)
                                +HEFT::gg_n2lo_lzbar0_lz1_L_1(z)
                                +HEFT::gg_n2lo_lzbar0_lz2_L_1(z)
                                +HEFT::gg_n2lo_lzbar1_L_1(z)
                                +HEFT::gg_n2lo_lzbar2_L_1(z));
        // L^2
        _me_lf_coeffs.push_back(
                                +HEFT::gg_n2lo_lzbar1_L_2(z)
                                +HEFT::gg_n2lo_lzbar0_lz0_L_2(z)
                                +HEFT::gg_n2lo_lzbar0_lz1_L_2(z));
    }
    
};


class HiggsEFTGGN3LOReg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsEFTGGN3LOReg(){set_dimensions(2);_channel = "gg";_name="gg_reg_n3lo";
        set_rel_accuracy_multiplier(3.0);
        //set_integration_method("vegas");
        _truncation_order=37;
        _full_log_switch = 0.0;
        }
    void compute_matrix_element_Lf_coefficients(const double* xx){
        const double z=xx[1];
        _me_lf_coeffs.clear();
        // L^0
        _me_lf_coeffs.push_back(evaluator.n3lo_reg_complete_gg(z,0));
        // L^1
        // minus sign for Bernhard's opposite Lf convention
        _me_lf_coeffs.push_back(-evaluator.n3lo_reg_complete_gg(z,1));
        // L^2
        _me_lf_coeffs.push_back(evaluator.n3lo_reg_complete_gg(z,2));
        // L^3
        // minus sign for Bernhard's opposite Lf convention
        _me_lf_coeffs.push_back(-evaluator.n3lo_reg_complete_gg(z,3));
        
    }
    
    void SetTruncationOrder(int i){_truncation_order=i;}
    void SetFullLogSwitch(double x){_full_log_switch=x;}
private:
    int _truncation_order;
    double _full_log_switch;
    HEFT::N3LORegEvaluator evaluator;
};



class HiggsEFTQGNLOReg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsEFTQGNLOReg(){set_dimensions(2);_channel = "qg";_name="qg_reg_nlo";
        }
    void compute_matrix_element_Lf_coefficients(const double* xx){
        _me_lf_coeffs.clear();
        _me_lf_coeffs.push_back(HEFT::qg_nlo_r_L0(xx[1]));
        _me_lf_coeffs.push_back(HEFT::qg_nlo_r_L1(xx[1]));
    }
    
};

class HiggsEFTQGNNLOReg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsEFTQGNNLOReg(){set_dimensions(2);_channel = "qg";_name="qg_reg_n2lo";
        set_rel_accuracy_multiplier(2.5);
        }
    void compute_matrix_element_Lf_coefficients(const double* xx){
        const double z=xx[1];
        _me_lf_coeffs.clear();
        // L^0
        _me_lf_coeffs.push_back( HEFT::qg_n2lo_r_lz0_L0(z)
                                +HEFT::qg_n2lo_r_lz1_L0(z)
                                +HEFT::qg_n2lo_r_lz2_L0(z)
                                +HEFT::qg_n2lo_r_lz3_L0(z));
        // L^1
        _me_lf_coeffs.push_back(+HEFT::qg_n2lo_r_lz0_L1(z)
                                +HEFT::qg_n2lo_r_lz1_L1(z)
                                +HEFT::qg_n2lo_r_lz2_L1(z));
        // L^2
        _me_lf_coeffs.push_back(
                                +HEFT::qg_n2lo_r_lz0_L2(z)
                                +HEFT::qg_n2lo_r_lz1_L2(z)
                                );
    }
    
};

class HiggsEFTQGN3LOReg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsEFTQGN3LOReg(){set_dimensions(2);_channel = "qg";_name="qg_reg_n3lo";
        set_rel_accuracy_multiplier(10.);
        _truncation_order=37;
        _full_log_switch = 0.0;
        }
    void compute_matrix_element_Lf_coefficients(const double* xx){
        const double z=xx[1];
        _me_lf_coeffs.clear();
        // L^0
        _me_lf_coeffs.push_back(evaluator.n3lo_reg_complete_qg(z,0));
        // L^1
        _me_lf_coeffs.push_back(-evaluator.n3lo_reg_complete_qg(z,1));
        // L^2
        _me_lf_coeffs.push_back(evaluator.n3lo_reg_complete_qg(z,2));
        // L^3
        // minus sign for Bernhard's opposite Lf convention
        _me_lf_coeffs.push_back(-evaluator.n3lo_reg_complete_qg(z,3));
        
        
                
    }
    void SetFullLogSwitch(double x){_full_log_switch=x;}
    void SetTruncationOrder(int i){_truncation_order=i;}
private:
    int _truncation_order;
    double _full_log_switch;
    HEFT::N3LORegEvaluator evaluator;
    
};

class HiggsEFTQQBARNLOReg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsEFTQQBARNLOReg(){set_dimensions(2);_channel = "qqbar";_name="qqb_reg_nlo";
        set_rel_accuracy_multiplier(100.);
    }
    void compute_matrix_element_Lf_coefficients(const double* xx){
        _me_lf_coeffs.clear();
        _me_lf_coeffs.push_back(HEFT::qqbar_nlo_r_L0(xx[1]));
    }
    
};

class HiggsEFTQQBARNNLOReg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsEFTQQBARNNLOReg(){set_dimensions(2);_channel = "qqbar";_name="qqb_reg_n2lo";
        set_rel_accuracy_multiplier(100.);
        }
    void compute_matrix_element_Lf_coefficients(const double* xx){
        const double z=xx[1];
        _me_lf_coeffs.clear();
        // L^0
        _me_lf_coeffs.push_back( HEFT::qqb_nnlo_r_lz0_L0(z)
                                +HEFT::qqb_nnlo_r_lz1_L0(z)
                                +HEFT::qqb_nnlo_r_lz2_L0(z)
                                );
        // L^1
        _me_lf_coeffs.push_back(+HEFT::qqb_nnlo_r_lz0_L1(z)
                                +HEFT::qqb_nnlo_r_lz1_L1(z)
                                );
        // L^2
        _me_lf_coeffs.push_back(
                                +HEFT::qqb_nnlo_r_lz0_L2(z)
                                );
    }
    
};

class HiggsEFTQQBARN3LOReg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsEFTQQBARN3LOReg(){set_dimensions(2);_channel = "qqbar";_name="qqb_reg_n3lo";
    
        set_rel_accuracy_multiplier(500.);
        _truncation_order=37;
        }
    void compute_matrix_element_Lf_coefficients(const double* xx){
        const double z=xx[1];
        _me_lf_coeffs.clear();
        // L^0
        _me_lf_coeffs.push_back(evaluator.n3lo_reg_complete_qqbar(z,0));
        // L^1
        // minus sign for Bernhard's opposite Lf convention
        _me_lf_coeffs.push_back(-evaluator.n3lo_reg_complete_qqbar(z,1));
        // L^2
        _me_lf_coeffs.push_back(evaluator.n3lo_reg_complete_qqbar(z,2));
        // L^3
        _me_lf_coeffs.push_back(-evaluator.n3lo_reg_complete_qqbar(z,3));
    }
    void SetTruncationOrder(int i){_truncation_order=i;}

private:
    int _truncation_order;
    HEFT::N3LORegEvaluator evaluator;
};

class HiggsEFTQQNNLOReg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsEFTQQNNLOReg(){set_dimensions(2);_channel = "qq";_name="qq_reg_n2lo";
    
        set_rel_accuracy_multiplier(100.);
        }
    void compute_matrix_element_Lf_coefficients(const double* xx){
        const double z=xx[1];
        _me_lf_coeffs.clear();
        // L^0
        _me_lf_coeffs.push_back( HEFT::qq_n2lo_lz0_L0(z)
                                +HEFT::qq_n2lo_lz1_L0(z)
                                +HEFT::qq_n2lo_lz2_L0(z)
                                );
        // L^1
        _me_lf_coeffs.push_back(+HEFT::qq_n2lo_lz0_L1(z)
                                +HEFT::qq_n2lo_lz1_L1(z)
                                );
        // L^2
        _me_lf_coeffs.push_back(
                                +HEFT::qq_n2lo_lz0_L2(z)
                                );
    }
    
};

class HiggsEFTQQN3LOReg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsEFTQQN3LOReg(){set_dimensions(2);_channel = "qq";_name="qq_reg_n3lo";
        
        set_rel_accuracy_multiplier(200.);
        _truncation_order=37;
        }
    void compute_matrix_element_Lf_coefficients(const double* xx){
        const double z=xx[1];
        _me_lf_coeffs.clear();
        // L^0
        _me_lf_coeffs.push_back(evaluator.n3lo_reg_complete_qq(z,0));
        // L^1
        // minus sign for Bernhard's opposite Lf convention
        _me_lf_coeffs.push_back(-evaluator.n3lo_reg_complete_qq(z,1));
        // L^2
        _me_lf_coeffs.push_back(evaluator.n3lo_reg_complete_qq(z,2));
        // L^3
        _me_lf_coeffs.push_back(-evaluator.n3lo_reg_complete_qq(z,3));
    }
    void SetTruncationOrder(int i){_truncation_order=i;}

private:
    int _truncation_order;
    HEFT::N3LORegEvaluator evaluator;
};


class HiggsEFTQ1Q2NNLOReg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsEFTQ1Q2NNLOReg(){set_dimensions(2);_channel = "q1q2";
        _name="q1q2_reg_n2lo";
        
        set_rel_accuracy_multiplier(50.);
        }
    void compute_matrix_element_Lf_coefficients(const double* xx){
        const double z=xx[1];
        _me_lf_coeffs.clear();
        // L^0
        _me_lf_coeffs.push_back( HEFT::q1q2_n2lo_lz0_L0(z)
                                +HEFT::q1q2_n2lo_lz1_L0(z)
                                +HEFT::q1q2_n2lo_lz2_L0(z)
                                );
        // L^1
        _me_lf_coeffs.push_back(+HEFT::q1q2_n2lo_lz0_L1(z)
                                +HEFT::q1q2_n2lo_lz1_L1(z)
                                );
        // L^2
        _me_lf_coeffs.push_back(
                                +HEFT::q1q2_n2lo_lz0_L2(z)
                                );
        
    }

};


class HiggsEFTQ1Q2N3LOReg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsEFTQ1Q2N3LOReg(){set_dimensions(2);_channel = "q1q2";_name="q1q2_reg_n3lo";
        
        set_rel_accuracy_multiplier(10000.);
        _truncation_order=37;

    }
    void compute_matrix_element_Lf_coefficients(const double* xx){
        const double z=xx[1];
        _me_lf_coeffs.clear();
        // L^0
        _me_lf_coeffs.push_back(evaluator.n3lo_reg_complete_q1q2(z,0));
        // L^1
        // minus sign for Bernhard's opposite Lf convention
        _me_lf_coeffs.push_back(-evaluator.n3lo_reg_complete_q1q2(z,1));
        // L^2
        _me_lf_coeffs.push_back(evaluator.n3lo_reg_complete_q1q2(z,2));
        // L^3
        // minus sign for Bernhard's opposite Lf convention
        _me_lf_coeffs.push_back(-evaluator.n3lo_reg_complete_q1q2(z,3));
    }
    void SetTruncationOrder(int i){_truncation_order=i;}
private:
    int _truncation_order;
    HEFT::N3LORegEvaluator evaluator;
};

#endif
