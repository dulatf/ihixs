
#include "constants.h"
#include<vector>
#include "as_series.h"
#include "higgs_eft.h"
#include <complex>
using namespace std;

// gluon - gluon reg
namespace HEFT{

    // nlo
    double nlo_reg(const double& z, const double& L){
        return   nlo_r_lz0(z,L)
                +nlo_r_lz1(z,L);
    }

    double nlo_r_lz0(const double& z, const double& L)
    {
        const double zb=1.-z;
        return 1./z * (
                       -(11./2.)*pow(zb,3.)
                       - 6.*pow(z*z-z+1.,2.)*log(z)/zb
                       )
        +
        L * 1./z*6.*(pow(z,3.)-pow(z,2.)+2.*z-1.)
        ;
    }

    double nlo_r_lz1(const double& z, const double& L)
    {
        return -1./z * 12.*(pow(z,3.)-z*z+2.*z-1.)*log(1.-z);
    }
    
    
    double nlo_reg(const double& z){
        return nlo_r_lz0(z,0.)+nlo_r_lz1(z,0.);
    }
    
    double nlo_reg_L(const double& z)
    {
        return 1./z*6.*(pow(z,3.)-pow(z,2.)+2.*z-1.);
    }
    

   
    
   
}
