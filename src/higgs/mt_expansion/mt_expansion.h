#ifndef MT_EXPANSION_H
#define MT_EXPANSION_H


#include "luminosity_integrals.h"


namespace MTEXP {
    double gg_nlo_delta(const double& rho);
    double gg_nlo_D0(const double& rho);
    double gg_nlo_D1(const double& rho);
    double gg_nlo_D0_L(const double& rho);
    double gg_nlo_reg(const double* xx);
    
    
    double gg_nnlo_delta(const double& rho,const double& lmH, const double& lmtop);
    double gg_nnlo_D0(const double& rho,const double& lmH);
    double gg_nnlo_D1(const double& rho,const double& lmH);
    double gg_nnlo_D2(const double& rho,const double& lmH);
    double gg_nnlo_D3(const double& rho,const double& lmH);
    
    
    
    
    double linear_interpolate(double* my_table, const double&tau);
    
    double log_expanded(const double& z);
}





#endif







