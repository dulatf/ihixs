#ifndef MT_EXPANSION_GG_H
#define MT_EXPANSION_GG_H


#include "luminosity_integrals.h"
#include "mt_expansion.h"


class HiggsMtExpansion_gg_nnlo_delta: public LuminosityIntegralAbsolutelyGenericDelta
{
public:
    HiggsMtExpansion_gg_nnlo_delta(){
        set_dimensions(1);
        _channel = "gg";
        _name="gg_delta_nnlo_mt_exp";
    }
};

class HiggsMtExpansion_gg_nnlo_plus: public LuminosityIntegralAbsolutelyGenericPlus
{
public:
    HiggsMtExpansion_gg_nnlo_plus(int m){set_dimensions(2);_channel = "gg";_name="gg_D"+to_string(m)+"_nnlo_mt_exp";
        
        _log_power=m;}
};



class HiggsMtExpansion_gg_nnlo_reg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsMtExpansion_gg_nnlo_reg(const double& rho);
    void compute_matrix_element_Lf_coefficients(const double* xx){
        _me_lf_coeffs.clear();
        _me_lf_coeffs.push_back(matched_L0(xx[1]));
        _me_lf_coeffs.push_back(matched_L1(xx[1]));
        _me_lf_coeffs.push_back(matched_L2(xx[1]));
    }
    
    
    double matched_L0(const double& z);
    double matched_L1(const double& z);
    double matched_L2(const double& z);
    
    
    double z_times_reg_L0(const double& z);
    double z_times_reg_L1(const double& z);
    double z_times_reg_L2(const double& z);
    
    
    
private:
    double _rho;
    double _Bgg2;
    double _Agg2;
    double _Bgg1;
    bool _matching;
    double cL0[4][4][17];
    double cL1[4][4][17];
    double cL2[4][4][17];
    
private:
    void compute_Agg2();
    void compute_Bgg2(){_Bgg2 = z_times_reg_L0(1e-14)*1.;}
    void compute_Bgg1();
};










#endif







