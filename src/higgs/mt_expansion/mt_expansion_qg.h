#ifndef MT_EXPANSION_QG_H
#define MT_EXPANSION_QG_H


#include "luminosity_integrals.h"
#include "mt_expansion.h"

class HiggsMtExpansion_qg_nnlo_reg: public LuminosityIntegralAbsolutelyGenericReg
{
public:
    HiggsMtExpansion_qg_nnlo_reg(const double& rho);
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
    double _Bqg2;
    double _Aqg2;
    double _Bgg1;
    bool _matching;
    
    vector<double> cL0;
    double cL1[4][4][17];
    double cL2[4][4][17];
private:
    void compute_Aqg2();
    void compute_Bqg2(){_Bqg2 = z_times_reg_L0(1e-14)*1.0;}
    void compute_Bgg1();
};








#endif







