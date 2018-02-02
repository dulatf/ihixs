
#include "constants.h"
#include<vector>
#include "as_series.h"
#include "higgs_eft.h"
#include "cppchaplin.h"
#include <complex>
using namespace std;

// quark - gluon reg
namespace HEFT{
    
    
    double qg_nlo_reg(const double& z, const double& L){
        return qg_nlo_r_lz0(z,L) + qg_nlo_r_lz1(z,L);
    }
    
    

    
    double qg_nlo_r_lz0(const double& z, const double& L)
    {
        const double zb=1.-z;
        return 1./z * (
                       -1./3.*(pow(zb,2.)+2.*(1.-2*z))
                       -2./3.*(1.+pow(zb,2.))*log(z)
                       -2./3.*(1.+pow(zb,2.))*L
                       )
        ;
        
    }
    
    double qg_nlo_r_lz1(const double& z, const double& L)
    {
        const double zb=1.-z;
        return 1./z * 4./3. * (1.+pow(zb,2.))*log(zb);
    }
    

    
    
    double qg_nlo_r_L0(const double& z){
        const double zb=1.-z;
        return 1./z * (
                       -1./3.*(pow(zb,2.)+2.*(1.-2*z))
                       -2./3.*(1.+pow(zb,2.))*log(z)
                       )
        +1./z * 4./3. * (1.+pow(zb,2.))*log(zb)
        ;
    }

    double qg_nlo_r_L1(const double& z){
        const double zb=1.-z;
        return 1./z * (
                       -2./3.*(1.+pow(zb,2.))
                       )
        ;
    }
    
    
    double qg_nnlo_reg(const double& z, const double& L){
        return qg_nnlo_r_lz0(z,L)
        + qg_nnlo_r_lz1(z,L)
        + qg_nnlo_r_lz2(z,L)
        + qg_nnlo_r_lz3(z,L);
    }
    
    
    double qg_nnlo_r_lz0(const double& z, const double& L)
    {
        return    qg_nnlo_r_lz0_const(z,L)
        + qg_nnlo_r_lz0_logz(z,L)
        + qg_nnlo_r_lz0_logz_sq(z,L)
        + qg_nnlo_r_lz0_logz_cube(z,L) ;
    }
    double qg_nnlo_r_lz0_const(const double& z, const double& L)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2 = consts::z2;
        const double z3 = consts::z3;
        complex<double> res =
        (3*(-1 + pow(Nc,2))*(2 + 2*z + pow(z,2))*z2*chaplin::HPL(-1,z))/(4.*z) - ((-1 + pow(Nc,2))*(2 + 2*z - pow(z,2) + 4*pow(Nc,2)*(3 + 18*z + 8*pow(z,2)))*z2*chaplin::HPL(1,z))/(8.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(-576 - 2592*pow(Nc,2) - 1728*L*pow(Nc,2) - 864*z - 3456*pow(Nc,2)*z - 1728*L*pow(Nc,2)*z - 864*pow(z,2) - 432*pow(Nc,2)*pow(z,2) - 864*L*pow(Nc,2)*pow(z,2) - 576*pow(z,3))*chaplin::HPL(-1,0,z))/(3456.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(1656 + 864*L - 21024*pow(Nc,2) + 4320*L*pow(Nc,2) - 3024*z + 7344*pow(Nc,2)*z + 15552*L*pow(Nc,2)*z - 2268*pow(z,2) + 5940*pow(Nc,2)*pow(z,2) + 7776*L*pow(Nc,2)*pow(z,2) - 288*pow(z,3) - 1152*pow(Nc,2)*pow(z,3))*chaplin::HPL(1,0,z))/(3456.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(3456*pow(Nc,2) + 3456*pow(Nc,2)*z + 1728*pow(Nc,2)*pow(z,2))*chaplin::HPL(-1,-1,0,z))/(3456.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(-1728*pow(Nc,2) - 1728*pow(Nc,2)*z - 864*pow(Nc,2)*pow(z,2))*chaplin::HPL(-1,0,0,z))/(3456.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(3456*pow(Nc,2) + 3456*pow(Nc,2)*z + 1728*pow(Nc,2)*pow(z,2))*chaplin::HPL(-1,1,0,z))/(3456.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(-864 - 14688*pow(Nc,2) - 432*z - 34128*pow(Nc,2)*z + 216*pow(z,2) - 19224*pow(Nc,2)*pow(z,2))*chaplin::HPL(0,1,0,z))/(3456.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(3456*pow(Nc,2) + 3456*pow(Nc,2)*z + 1728*pow(Nc,2)*pow(z,2))*chaplin::HPL(1,-1,0,z))/(3456.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(-3456 - 12960*pow(Nc,2) + 1728*z - 44064*pow(Nc,2)*z - 864*pow(z,2) - 21168*pow(Nc,2)*pow(z,2))*chaplin::HPL(1,0,0,z))/(3456.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(-864 - 6912*pow(Nc,2) - 864*z - 32832*pow(Nc,2)*z + 432*pow(z,2) - 14688*pow(Nc,2)*pow(z,2))*chaplin::HPL(1,1,0,z))/(3456.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(1393 + 636*L - 44471*pow(Nc,2) - 29700*L*pow(Nc,2) - 8280*pow(L,2)*pow(Nc,2) + 2120*Nc*nf + 1392*L*Nc*nf + 288*pow(L,2)*Nc*nf + 780*z - 648*L*z - 216*pow(L,2)*z + 39068*pow(Nc,2)*z + 21336*L*pow(Nc,2)*z + 6984*pow(L,2)*pow(Nc,2)*z - 3344*Nc*nf*z - 1824*L*Nc*nf*z - 288*pow(L,2)*Nc*nf*z - 51*pow(z,2) - 216*L*pow(z,2) + 54*pow(L,2)*pow(z,2) + 2393*pow(Nc,2)*pow(z,2) + 4560*L*pow(Nc,2)*pow(z,2) - 198*pow(L,2)*pow(Nc,2)*pow(z,2) + 1432*Nc*nf*pow(z,2) + 912*L*Nc*nf*pow(z,2) + 144*pow(L,2)*Nc*nf*pow(z,2) - 2608*pow(z,3) - 744*L*pow(z,3) + 1656*pow(Nc,2)*pow(z,3) + 1128*L*pow(Nc,2)*pow(z,3) + 864*pow(L,2)*pow(Nc,2)*pow(z,3) - 2520*z2 - 864*L*z2 + 3816*pow(Nc,2)*z2 + 6048*L*pow(Nc,2)*z2 + 2592*z*z2 + 1728*L*z*z2 - 12096*pow(Nc,2)*z*z2 + 12096*L*pow(Nc,2)*z*z2 - 5508*pow(z,2)*z2 - 864*L*pow(z,2)*z2 + 2916*pow(Nc,2)*pow(z,2)*z2 + 8640*L*pow(Nc,2)*pow(z,2)*z2 - 576*pow(z,3)*z2 - 4608*pow(Nc,2)*pow(z,3)*z2 + 864*z3 + 4320*pow(Nc,2)*z3 - 7776*pow(Nc,2)*z*z3 + 4320*pow(Nc,2)*pow(z,2)*z3))/(3456.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }
    
    double qg_nnlo_r_lz0_logz(const double& z, const double& L)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2 = consts::z2;
        complex<double> res =
        -((-1 + pow(Nc,2))*(-106 + 4910*pow(Nc,2) + 2232*L*pow(Nc,2) + 432*pow(L,2)*pow(Nc,2) - 232*Nc*nf - 96*L*Nc*nf + 621*z + 72*L*z + 36*pow(L,2)*z - 5941*pow(Nc,2)*z - 2760*L*pow(Nc,2)*z + 396*pow(L,2)*pow(Nc,2)*z + 304*Nc*nf*z + 96*L*Nc*nf*z - 642*pow(z,2) - 342*L*pow(z,2) - 18*pow(L,2)*pow(z,2) - 1516*pow(Nc,2)*pow(z,2) + 390*L*pow(Nc,2)*pow(z,2) + 450*pow(L,2)*pow(Nc,2)*pow(z,2) - 152*Nc*nf*pow(z,2) - 48*L*Nc*nf*pow(z,2) + 112*pow(z,3) - 48*L*pow(z,3) - 560*pow(Nc,2)*pow(z,3) - 528*L*pow(Nc,2)*pow(z,3) + 144*z2 - 432*pow(Nc,2)*z2 - 360*z*z2 - 1080*pow(Nc,2)*z*z2 + 180*pow(z,2)*z2 - 1044*pow(Nc,2)*pow(z,2)*z2))/(576.*pow(Nc,2)*z) - ((-1 + pow(Nc,2))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,0,z))/(4.*z) + ((-1 + pow(Nc,2))*(1 + pow(Nc,2)*(5 + 18*z + 9*pow(z,2)))*chaplin::HPL(1,0,z))/(2.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*log(z);
    }
    
    double qg_nnlo_r_lz0_logz_sq(const double& z, const double& L)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        complex<double> res =
        ((-1 + pow(Nc,2))*(16*Nc*nf*(2 - 2*z + pow(z,2)) + z*(84 + 24*L*(-2 + z) + 189*z + 80*pow(z,2)) - pow(Nc,2)*(744 - 1292*z + 229*pow(z,2) - 192*pow(z,3) + 24*L*(12 + 14*z + 15*pow(z,2)))))/(384.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(z),2.);
    }
    
    double qg_nnlo_r_lz0_logz_cube(const double& z, const double& L)
    {
        const double Nc = QCD::Nc;
        complex<double> res =
        -((-1 + pow(Nc,2))*(-5*(-2 + z)*z + pow(Nc,2)*(48 + 62*z + 65*pow(z,2))))/(192.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(z),3.);
    }
    
    double qg_nnlo_r_lz1(const double& z, const double& L)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2 = consts::z2;
        complex<double> res =
        ((-1 + pow(Nc,2))*(-439 - 162*L - 36*pow(L,2) + 4781*pow(Nc,2) + 2766*L*pow(Nc,2) + 252*pow(L,2)*pow(Nc,2) - 156*Nc*nf - 48*L*Nc*nf + 648*z + 288*L*z + 36*pow(L,2)*z - 3360*pow(Nc,2)*z - 2496*L*pow(Nc,2)*z - 252*pow(L,2)*pow(Nc,2)*z + 192*Nc*nf*z + 48*L*Nc*nf*z - 216*pow(z,2) - 126*L*pow(z,2) - 18*pow(L,2)*pow(z,2) - 828*pow(Nc,2)*pow(z,2) + 6*L*pow(Nc,2)*pow(z,2) + 126*pow(L,2)*pow(Nc,2)*pow(z,2) - 108*Nc*nf*pow(z,2) - 24*L*Nc*nf*pow(z,2) + 124*pow(z,3) - 188*pow(Nc,2)*pow(z,3) - 288*L*pow(Nc,2)*pow(z,3) + 216*z2 - 936*pow(Nc,2)*z2 - 360*z*z2 - 2088*pow(Nc,2)*z*z2 + 180*pow(z,2)*z2 - 1404*pow(Nc,2)*pow(z,2)*z2))/(288.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(-27 + 42*z - 66*pow(z,2) - 8*pow(z,3) - 4*Nc*nf*(2 - 2*z + pow(z,2)) + pow(Nc,2)*(373 - 458*z + 46*pow(z,2) - 88*pow(z,3)) + 12*L*(-pow(-1 + z,2) + pow(Nc,2)*(17 + 6*z + 15*pow(z,2))))*chaplin::HPL(0,z))/(48.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(-2 + 6*z - 3*pow(z,2) + pow(Nc,2)*(34 + 18*z + 35*pow(z,2)))*pow(chaplin::HPL(0,z),2))/(16.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,0,z))/(2.*z) + ((-1 + pow(Nc,2))*(-72 - 648*pow(Nc,2) - 72*z - 2664*pow(Nc,2)*z + 36*pow(z,2) - 1260*pow(Nc,2)*pow(z,2))*chaplin::HPL(1,0,z))/(288.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*log(1.-z);
    }
    
    double qg_nnlo_r_lz2(const double& z, const double& L)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        complex<double> res =
        ((-1 + pow(Nc,2))*(180 - 1748*pow(Nc,2) + 8*Nc*nf - 282*z + 1550*pow(Nc,2)*z - 8*Nc*nf*z + 135*pow(z,2) + 59*pow(Nc,2)*pow(z,2) + 4*Nc*nf*pow(z,2) + 192*pow(Nc,2)*pow(z,3) - 36*L*(-1 + 7*pow(Nc,2))*(2 - 2*z + pow(z,2))))/(192.*pow(Nc,2)*z) - ((-1 + pow(Nc,2))*(-6 + 10*z - 5*pow(z,2) + pow(Nc,2)*(72 + 20*z + 62*pow(z,2)))*chaplin::HPL(0,z))/(16.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(1.-z),2.);
    }
    
    double qg_nnlo_r_lz3(const double& z, const double& L)
    {
        const double Nc = QCD::Nc;
        complex<double> res =
        ((13 - 96*pow(Nc,2) + 83*pow(Nc,4))*(2 - 2*z + pow(z,2)))/(96.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(1.-z),3.);
    }
    
    
    // mathematica file for qg_n2lo_lzbar*_lz*_L* : Higgs_falko_qg_n2

    double qg_n2lo_r_lz0_L0(const double& z)
    {
        //const double nf = consts::nf;
        const double Nc = QCD::Nc;
        //const double z2 = consts::z2;
        //const double z3 = consts::z3;
        complex<double> res =-((-1 + Nc)*(1 + Nc)*(-5184*pow(Nc,2)*consts::z2 - 5184*pow(Nc,2)*z*consts::z2 - 2592*pow(Nc,2)*pow(z,2)*consts::z2)*chaplin::HPL(-1,z))/(3456.*pow(Nc,2)*z) - ((-1 + Nc)*(1 + Nc)*(864*consts::z2 + 5184*pow(Nc,2)*consts::z2 + 864*z*consts::z2 + 31104*pow(Nc,2)*z*consts::z2 - 432*pow(z,2)*consts::z2 + 13824*pow(Nc,2)*pow(z,2)*consts::z2)*chaplin::HPL(1,z))/(3456.*pow(Nc,2)*z) - ((-1 + Nc)*(1 + Nc)*(-3456*pow(Nc,2) - 3456*pow(Nc,2)*z - 1728*pow(Nc,2)*pow(z,2))*chaplin::HPL(-1,-1,0,z))/(3456.*pow(Nc,2)*z) - ((-1 + Nc)*(1 + Nc)*(1728*pow(Nc,2) + 1728*pow(Nc,2)*z + 864*pow(Nc,2)*pow(z,2))*chaplin::HPL(-1,0,0,z))/(3456.*pow(Nc,2)*z) - ((-1 + Nc)*(1 + Nc)*(-3456*pow(Nc,2) - 3456*pow(Nc,2)*z - 1728*pow(Nc,2)*pow(z,2))*chaplin::HPL(-1,1,0,z))/(3456.*pow(Nc,2)*z) - ((-1 + Nc)*(1 + Nc)*(864 + 14688*pow(Nc,2) + 432*z + 34128*pow(Nc,2)*z - 216*pow(z,2) + 19224*pow(Nc,2)*pow(z,2))*chaplin::HPL(0,1,0,z))/(3456.*pow(Nc,2)*z) - ((-1 + Nc)*(1 + Nc)*(-3456*pow(Nc,2) - 3456*pow(Nc,2)*z - 1728*pow(Nc,2)*pow(z,2))*chaplin::HPL(1,-1,0,z))/(3456.*pow(Nc,2)*z) - ((-1 + Nc)*(1 + Nc)*(3456 + 12960*pow(Nc,2) - 1728*z + 44064*pow(Nc,2)*z + 864*pow(z,2) + 21168*pow(Nc,2)*pow(z,2))*chaplin::HPL(1,0,0,z))/(3456.*pow(Nc,2)*z) - ((-1 + Nc)*(1 + Nc)*(864 + 6912*pow(Nc,2) + 864*z + 32832*pow(Nc,2)*z - 432*pow(z,2) + 14688*pow(Nc,2)*pow(z,2))*chaplin::HPL(1,1,0,z))/(3456.*pow(Nc,2)*z) - ((-1 + Nc)*(1 + Nc)*chaplin::HPL(1,0,z)*(-1656 + 21024*pow(Nc,2) + 3024*z - 7344*pow(Nc,2)*z + 2268*pow(z,2) - 5940*pow(Nc,2)*pow(z,2) + 288*pow(z,3) + 1152*pow(Nc,2)*pow(z,3) - 1728*log(z) - 8640*pow(Nc,2)*log(z) - 31104*pow(Nc,2)*z*log(z) - 15552*pow(Nc,2)*pow(z,2)*log(z)))/(3456.*pow(Nc,2)*z) - ((-1 + Nc)*(1 + Nc)*chaplin::HPL(-1,0,z)*(576 + 2592*pow(Nc,2) + 864*z + 3456*pow(Nc,2)*z + 864*pow(z,2) + 432*pow(Nc,2)*pow(z,2) + 576*pow(z,3) + 1728*pow(Nc,2)*log(z) + 1728*pow(Nc,2)*z*log(z) + 864*pow(Nc,2)*pow(z,2)*log(z)))/(3456.*pow(Nc,2)*z) - ((-1 + Nc)*(1 + Nc)*(-1393 + 44471*pow(Nc,2) - 2120*Nc*consts::nf - 780*z - 39068*pow(Nc,2)*z + 3344*Nc*consts::nf*z + 51*pow(z,2) - 2393*pow(Nc,2)*pow(z,2) - 1432*Nc*consts::nf*pow(z,2) + 2608*pow(z,3) - 1656*pow(Nc,2)*pow(z,3) + 2520*consts::z2 - 3816*pow(Nc,2)*consts::z2 - 2592*z*consts::z2 + 12096*pow(Nc,2)*z*consts::z2 + 5508*pow(z,2)*consts::z2 - 2916*pow(Nc,2)*pow(z,2)*consts::z2 + 576*pow(z,3)*consts::z2 + 4608*pow(Nc,2)*pow(z,3)*consts::z2 - 636*log(z) + 29460*pow(Nc,2)*log(z) - 1392*Nc*consts::nf*log(z) + 3726*z*log(z) - 35646*pow(Nc,2)*z*log(z) + 1824*Nc*consts::nf*z*log(z) - 3852*pow(z,2)*log(z) - 9096*pow(Nc,2)*pow(z,2)*log(z) - 912*Nc*consts::nf*pow(z,2)*log(z) + 672*pow(z,3)*log(z) - 3360*pow(Nc,2)*pow(z,3)*log(z) + 864*consts::z2*log(z) - 2592*pow(Nc,2)*consts::z2*log(z) - 2160*z*consts::z2*log(z) - 6480*pow(Nc,2)*z*consts::z2*log(z) + 1080*pow(z,2)*consts::z2*log(z) - 6264*pow(Nc,2)*pow(z,2)*consts::z2*log(z) + 6696*pow(Nc,2)*pow(log(z),2) - 288*Nc*consts::nf*pow(log(z),2) - 756*z*pow(log(z),2) - 11628*pow(Nc,2)*z*pow(log(z),2) + 288*Nc*consts::nf*z*pow(log(z),2) - 1701*pow(z,2)*pow(log(z),2) + 2061*pow(Nc,2)*pow(z,2)*pow(log(z),2) - 144*Nc*consts::nf*pow(z,2)*pow(log(z),2) - 720*pow(z,3)*pow(log(z),2) - 1728*pow(Nc,2)*pow(z,3)*pow(log(z),2) + 864*pow(Nc,2)*pow(log(z),3) + 180*z*pow(log(z),3) + 1116*pow(Nc,2)*z*pow(log(z),3) - 90*pow(z,2)*pow(log(z),3) + 1170*pow(Nc,2)*pow(z,2)*pow(log(z),3) - 864*consts::z3 - 4320*pow(Nc,2)*consts::z3 + 7776*pow(Nc,2)*z*consts::z3 - 4320*pow(Nc,2)*pow(z,2)*consts::z3))/(3456.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }

    double qg_n2lo_r_lz0_L1(const double& z)
    {
        //const double nf = consts::nf;
        const double Nc = QCD::Nc;
        //const double z2 = consts::z2;
        //const double z3 = consts::z3;
        complex<double> res = -((-1 + Nc)*(1 + Nc)*(1728*pow(Nc,2) + 1728*pow(Nc,2)*z + 864*pow(Nc,2)*pow(z,2))*chaplin::HPL(-1,0,z))/(3456.*pow(Nc,2)*z) - ((-1 + Nc)*(1 + Nc)*(-864 - 4320*pow(Nc,2) - 15552*pow(Nc,2)*z - 7776*pow(Nc,2)*pow(z,2))*chaplin::HPL(1,0,z))/(3456.*pow(Nc,2)*z) - ((-1 + Nc)*(1 + Nc)*(-636 + 29700*pow(Nc,2) - 1392*Nc*consts::nf + 648*z - 21336*pow(Nc,2)*z + 1824*Nc*consts::nf*z + 216*pow(z,2) - 4560*pow(Nc,2)*pow(z,2) - 912*Nc*consts::nf*pow(z,2) + 744*pow(z,3) - 1128*pow(Nc,2)*pow(z,3) + 864*consts::z2 - 6048*pow(Nc,2)*consts::z2 - 1728*z*consts::z2 - 12096*pow(Nc,2)*z*consts::z2 + 864*pow(z,2)*consts::z2 - 8640*pow(Nc,2)*pow(z,2)*consts::z2 + 13392*pow(Nc,2)*log(z) - 576*Nc*consts::nf*log(z) + 432*z*log(z) - 16560*pow(Nc,2)*z*log(z) + 576*Nc*consts::nf*z*log(z) - 2052*pow(z,2)*log(z) + 2340*pow(Nc,2)*pow(z,2)*log(z) - 288*Nc*consts::nf*pow(z,2)*log(z) - 288*pow(z,3)*log(z) - 3168*pow(Nc,2)*pow(z,3)*log(z) + 2592*pow(Nc,2)*pow(log(z),2) + 432*z*pow(log(z),2) + 3024*pow(Nc,2)*z*pow(log(z),2) - 216*pow(z,2)*pow(log(z),2) + 3240*pow(Nc,2)*pow(z,2)*pow(log(z),2)))/(3456.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }
    
    double qg_n2lo_r_lz0_L2(const double& z)
    {
        
        const double Nc = QCD::Nc;
        
        complex<double> res = -((-1 + Nc)*(1 + Nc)*(8280*pow(Nc,2) - 288*Nc*consts::nf + 216*z - 6984*pow(Nc,2)*z + 288*Nc*consts::nf*z - 54*pow(z,2) + 198*pow(Nc,2)*pow(z,2) - 144*Nc*consts::nf*pow(z,2) - 864*pow(Nc,2)*pow(z,3) + 2592*pow(Nc,2)*log(z) + 216*z*log(z) + 2376*pow(Nc,2)*z*log(z) - 108*pow(z,2)*log(z) + 2700*pow(Nc,2)*pow(z,2)*log(z)))/(3456.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }
    
    double qg_n2lo_r_lz1_L0(const double& z)
    {
        //const double nf = consts::nf;
        const double Nc = QCD::Nc;
        //const double z2 = consts::z2;
        complex<double> res = ((-1 + pow(Nc,2))*(-439 + 4781*pow(Nc,2) - 156*Nc*consts::nf + 648*z - 3360*pow(Nc,2)*z + 192*Nc*consts::nf*z - 216*pow(z,2) - 828*pow(Nc,2)*pow(z,2) - 108*Nc*consts::nf*pow(z,2) + 124*pow(z,3) - 188*pow(Nc,2)*pow(z,3) + 216*consts::z2 - 936*pow(Nc,2)*consts::z2 - 360*z*consts::z2 - 2088*pow(Nc,2)*z*consts::z2 + 180*pow(z,2)*consts::z2 - 1404*pow(Nc,2)*pow(z,2)*consts::z2))/(288.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(-162 + 2238*pow(Nc,2) - 48*Nc*consts::nf + 252*z - 2748*pow(Nc,2)*z + 48*Nc*consts::nf*z - 396*pow(z,2) + 276*pow(Nc,2)*pow(z,2) - 24*Nc*consts::nf*pow(z,2) - 48*pow(z,3) - 528*pow(Nc,2)*pow(z,3))*chaplin::HPL(0,z))/(288.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(-36 + 612*pow(Nc,2) + 108*z + 324*pow(Nc,2)*z - 54*pow(z,2) + 630*pow(Nc,2)*pow(z,2))*pow(chaplin::HPL(0,z),2))/(288.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(288*pow(Nc,2) + 288*pow(Nc,2)*z + 144*pow(Nc,2)*pow(z,2))*chaplin::HPL(-1,0,z))/(288.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(-72 - 648*pow(Nc,2) - 72*z - 2664*pow(Nc,2)*z + 36*pow(z,2) - 1260*pow(Nc,2)*pow(z,2))*chaplin::HPL(1,0,z))/(288.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*log(1.-z);
    }
    
    double qg_n2lo_r_lz1_L1(const double& z)
    {
        //const double nf = consts::nf;
        const double Nc = QCD::Nc;
        //const double z2 = consts::z2;
        complex<double> res =((-1 + pow(Nc,2))*(-162 + 2766*pow(Nc,2) - 48*Nc*consts::nf + 288*z - 2496*pow(Nc,2)*z + 48*Nc*consts::nf*z - 126*pow(z,2) + 6*pow(Nc,2)*pow(z,2) - 24*Nc*consts::nf*pow(z,2) - 288*pow(Nc,2)*pow(z,3)))/(288.*pow(Nc,2)*z) + ((-1 + pow(Nc,2))*(-72 + 1224*pow(Nc,2) + 144*z + 432*pow(Nc,2)*z - 72*pow(z,2) + 1080*pow(Nc,2)*pow(z,2))*chaplin::HPL(0,z))/(288.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*log(1.-z);
    }
    
    double qg_n2lo_r_lz1_L2(const double& z)
    {
        //const double nf = consts::nf;
        const double Nc = QCD::Nc;
        //const double z2 = consts::z2;
        complex<double> res =((-1 + pow(Nc,2))*(-36 + 252*pow(Nc,2) + 36*z - 252*pow(Nc,2)*z - 18*pow(z,2) + 126*pow(Nc,2)*pow(z,2)))/(288.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*log(1.-z);
    }
    
    double qg_n2lo_r_lz2_L0(const double& z)
    {
        //const double nf = consts::nf;
        const double Nc = QCD::Nc;
        complex<double> res =((-1 + pow(Nc,2))*(180 - 1748*pow(Nc,2) + 8*Nc*consts::nf - 282*z + 1550*pow(Nc,2)*z - 8*Nc*consts::nf*z + 135*pow(z,2) + 59*pow(Nc,2)*pow(z,2) + 4*Nc*consts::nf*pow(z,2) + 192*pow(Nc,2)*pow(z,3)))/(192.*pow(Nc,2)*z) - ((-1 + pow(Nc,2))*(-6 + 10*z - 5*pow(z,2) + pow(Nc,2)*(72 + 20*z + 62*pow(z,2)))*chaplin::HPL(0,z))/(16.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(1.-z),2.);
    }

    double qg_n2lo_r_lz2_L1(const double& z)
    {
        //const double nf = consts::nf;
        const double Nc = QCD::Nc;
        complex<double> res = (-3*(-1 + pow(Nc,2))*(-1 + 7*pow(Nc,2))*(2 - 2*z + pow(z,2)))/(16.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(1.-z),2.);
    }
    
    double qg_n2lo_r_lz3_L0(const double& z)
    {
        const double Nc = QCD::Nc;
        complex<double> res =
        ((13 - 96*pow(Nc,2) + 83*pow(Nc,4))*(2 - 2*z + pow(z,2)))/(96.*pow(Nc,2)*z)
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(1.-z),3.);
    }
    
    
    //--------------------------- n3lo
    //: implementation of Bernhard's triple expansion (Nov 2017)
    //: we have here the full qg reg, including logs(mh^2/muf^2)
    double N3LORegEvaluator::n3lo_reg_complete_qg(const double& z, unsigned int log_muf_mh_squared_power) {
        
        //: main forking into three regions of [0,1]
        double res = 0.0;
        if (z<=1./13.) {
            const double x = sqrt(z);
            const double logx = log(z);
            for (int t = 0; t < NumZTerms; t++) {
                for (unsigned int p=0; p < 6; p++) {
                    //: [a_s][initial state][log_muf_mh_squared_power][x power][log(x) power]
                    //: [3: n3lo]  [0: gg] ....
                    res += pow(x,t-2) * pow(logx,p) * ZExp[3][1][log_muf_mh_squared_power][t][p];
                }
            }
        }
        else if (z>1.0/13.0 and z<=0.75) {
            const double x = 0.5-z;
            //const double logx = log(1-z);
            for (int t = 0; t < NumWTerms; t++) {
                //for (unsigned int p=0; p < 6; p++) {
                    //: [a_s][initial state][log_muf_mh_squared_power][x power][log(x) power]
                    //: [3: n3lo]  [0: gg] ....
                    res += pow(x,t) * WExp[3][1][log_muf_mh_squared_power][t][0];
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
                    res += pow(x,t) * pow(logx,p) * ZbExp[3][1][log_muf_mh_squared_power][t][p];
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
    
    
    
    
    
    double qg_n3lo_reg(const double& z, const double& L){
        return qg_n3lo_r_lz0_full_series(z,L)
        + qg_n3lo_r_lz1_full_series(z,L)
        + qg_n3lo_r_lz2_full_series(z,L)
        + qg_n3lo_r_lz3_full_series(z,L)
        + qg_n3lo_r_lz4_full_series(z,L)
        + qg_n3lo_r_lz5_full_series(z,L)
        + LEqgN3LOregFalko(z,L);
    }
    
    double qg_n3lo_reg_no_Lf(const double& z, const double& L){
        return qg_n3lo_r_lz0_full_series(z,L)
        + qg_n3lo_r_lz1_full_series(z,L)
        + qg_n3lo_r_lz2_full_series(z,L)
        + qg_n3lo_r_lz3_full_series(z,L)
        + qg_n3lo_r_lz4_full_series(z,L)
        + qg_n3lo_r_lz5_full_series(z,L);
    }
    
    
    double qg_n3lo_reg_no_Lf(const double& z){
        return qg_n3lo_r_lzX_series(0,z,30)
        +qg_n3lo_r_lzX_series(1,z,30)
        +qg_n3lo_r_lzX_series(2,z,30)
        +qg_n3lo_r_lzX_series(3,z,30)
        +qg_n3lo_r_lzX_series(4,z,30)
        +qg_n3lo_r_lzX_series(5,z,30)
        ;
    }
    
    
    double qg_n3lo_reg_no_Lf(const double& z,int truncation_order,
                             const double& log_switch=0.0) {
        if (log_switch > 0.0) {
            cout << "Error in eft_qg_reg.cpp: qg_n3lo_reg_no_Lf :  "
            << " log_switch flag > 0.0 but the full logs are switched off in this version" << endl;
            exit(0);
        }
        return qg_n3lo_r_lzX_series(0,z,truncation_order)
        +qg_n3lo_r_lzX_series(1,z,truncation_order)
        +qg_n3lo_r_lzX_series(2,z,truncation_order)
        +qg_n3lo_r_lzX_series(3,z,truncation_order)
        +qg_n3lo_r_lzX_series(4,z,truncation_order)
        +qg_n3lo_r_lzX_series(5,z,truncation_order)
        ;

    }
    
    
    double qg_n3lo_r_lz0(const double& z, const double& L)
    {
        // NS approximation only
        return qg_n3lo_r_lz0_NS(z,L);
    }
    
    double qg_n3lo_r_lz1(const double& z, const double& L)
    {
        // NS approximation only
        return qg_n3lo_r_lz1_NS(z,L);
    }
    
    double qg_n3lo_r_lz2(const double& z, const double& L)
    {
        // NS approximation only
        return qg_n3lo_r_lz2_NS(z,L);
    }
    
    
    
    
    double LEqgN3LOregFalko_L3(const double& z){
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2= consts::z2;
        //const double z3= consts::z3;
        //const double z4= consts::z4;
        
        
        complex<double> L3=((-1 + pow(Nc,2))*(4*Nc*nf*(-496 + 1140*z - 663*pow(z,2) + 64*pow(z,3)) - 4*pow(Nc,3)*nf*(-2636 + 2916*z - 771*pow(z,2) + 272*pow(z,3)) + 3*z*(6 - 33*z + 48*(-2 + z)*z2) - 12*pow(Nc,2)*(8*pow(nf,2)*(2 - 2*z + pow(z,2)) + 3*(-40 + 43*z + 20*pow(z,2)) + 24*(8 + 2*z + 11*pow(z,2))*z2) + pow(Nc,4)*(-129064 + 130890*z - 21165*pow(z,2) + 17344*pow(z,3) + 144*(128 + 118*z + 133*pow(z,2))*z2)))/(13824.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-1 - 2*pow(Nc,2)*(-1 + 6*z + 3*pow(z,2)) + 3*pow(Nc,4)*(9 + 32*z + 16*pow(z,2)))*chaplin::HPL(1,0,z))/(48.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(1 - 10*pow(Nc,2) + 37*pow(Nc,4))*(2 - 2*z + pow(z,2))*pow(log(1 - z),2))/(96.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(180*z - 72*pow(z,2) + 24*Nc*nf*(-32 + 10*z + 31*pow(z,2)) - 12*pow(Nc,2)*z*(476 - 43*z + 64*pow(z,2)) + 24*pow(Nc,3)*nf*(136 + 94*z + 73*pow(z,2)) + 12*pow(Nc,4)*(-4112 - 1019*z - 845*pow(z,2) + 288*pow(z,3)))*log(z))/(13824.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(6*(-2 + z)*z + 96*Nc*nf*(-2 + z)*z - 96*pow(Nc,3)*nf*(-2 + z)*z + 2*pow(Nc,2)*(-84*z - 102*pow(z,2)) + 6*pow(Nc,4)*(224 - 194*z + 481*pow(z,2)))*pow(log(z),2))/(4608.*pow(Nc,3)*z) + log(1 - z)*(-((-1 + pow(Nc,2))*(-3*(-4 + z)*z - 10*Nc*nf*(2 - 2*z + pow(z,2)) + 62*pow(Nc,3)*nf*(2 - 2*z + pow(z,2)) + pow(Nc,2)*(358 - 374*z + 49*pow(z,2) - 32*pow(z,3)) + 2*pow(Nc,4)*(-1333 + 1139*z - 82*pow(z,2) + 128*pow(z,3))))/(576.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(3*(-8*z + 4*pow(z,2)) - 24*pow(Nc,2)*(8 + 2*z + 11*pow(z,2)) + 12*pow(Nc,4)*(128 + 118*z + 133*pow(z,2)))*log(z))/(1152.*pow(Nc,3)*z));
        check_imaginary_part(L3,__PRETTY_FUNCTION__);
        
        return real(L3);
        
    }
    
    double LEqgN3LOregFalko_L2(const double& z){
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2= consts::z2;
        const double z3= consts::z3;
        //const double z4= consts::z4;
        
        complex<double> L2= -((-1 + pow(Nc,2))*(-12*pow(Nc,2)*(5041 - 3848*z + 140*pow(z,2) - 2476*pow(z,3) - 8*pow(nf,2)*(29 - 38*z + 19*pow(z,2)) + 24*(-202 + 297*z - 246*pow(z,2) + 36*pow(z,3))*z2 - 432*z3 - 288*z*z3 - 1944*pow(z,2)*z3) - 9*(120*pow(z,3) + 48*(-6 - 13*z + 20*pow(z,2) + 4*pow(z,3))*z2 + 24*(13 + 8*z3) + pow(z,2)*(293 + 48*z3) - 2*z*(439 + 240*z3)) + 4*pow(Nc,3)*nf*(-44290 + 4504*pow(z,3) + 72*(78 + 94*z + 31*pow(z,2))*z2 + 27*pow(z,2)*(-551 + 64*z3) - 6*z*(-8447 + 576*z3)) + 12*Nc*nf*(1446 - 312*pow(z,3) + 72*(-14 + 6*z + 9*pow(z,2))*z2 + pow(z,2)*(4723 - 576*z3) + 2*z*(-2987 + 576*z3)) + pow(Nc,4)*(1417364 - 93272*pow(z,3) + 144*(-2654 - 653*z - 602*pow(z,2) + 468*pow(z,3))*z2 - 82944*z3 - 198*z*(7397 + 144*z3) - 3*pow(z,2)*(-60223 + 68112*z3))))/(27648.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(1 - 4*pow(Nc,2) + 27*pow(Nc,4))*(2 + 2*z + pow(z,2))*z2*chaplin::HPL(-1,z))/(32.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(-16 + 2*z - pow(z,2) - 2*pow(Nc,2)*(-10 + 70*z + 41*pow(z,2)) + pow(Nc,4)*(324 + 1266*z + 623*pow(z,2)))*z2*chaplin::HPL(1,z))/(64.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(6 + 10*z + 4*pow(z,2) - 4*pow(Nc,3)*nf*(2 + 2*z + pow(z,2)) - pow(Nc,2)*(6 + 14*z + 5*pow(z,2)) + pow(Nc,4)*(292 + 240*z - pow(z,2) + 32*pow(z,3)))*chaplin::HPL(-1,0,z))/(64.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(48*Nc*nf*(-2 + 3*pow(z,2)) + 8*pow(Nc,3)*nf*(38 + 54*z + 15*pow(z,2)) + pow(Nc,2)*(-940 + 672*z + 612*pow(z,2) - 16*pow(z,3)) - 3*(-18 + 18*z + 31*pow(z,2) + 8*pow(z,3)) + pow(Nc,4)*(-1378 - 6162*z - 1395*pow(z,2) + 296*pow(z,3)))*chaplin::HPL(1,0,z))/(384.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(1 - 2*pow(Nc,2) + 13*pow(Nc,4))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,-1,0,z))/(16.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(1 - 3*pow(Nc,2) + 20*pow(Nc,4))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,0,0,z))/(32.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-3*pow(Nc,2)*(2 + 2*z + pow(z,2)) + 21*pow(Nc,4)*(2 + 2*z + pow(z,2)))*chaplin::HPL(-1,1,0,z))/(48.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(6*z + 3*pow(z,2) + pow(Nc,2)*(-12*z - 6*pow(z,2)) + 3*pow(Nc,4)*(48 - 46*z + 49*pow(z,2)))*chaplin::HPL(0,-1,0,z))/(96.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-2 + 2*z + 8*Nc*nf*(-2 + z)*z - 8*pow(Nc,3)*nf*(-2 + z)*z - pow(z,2) + pow(Nc,2)*(32 - 22*z + 7*pow(z,2)) + 2*pow(Nc,4)*(-31 - 74*z + 51*pow(z,2)))*chaplin::HPL(0,1,0,z))/(32.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-3*pow(Nc,2)*(2 + 2*z + pow(z,2)) + 21*pow(Nc,4)*(2 + 2*z + pow(z,2)))*chaplin::HPL(1,-1,0,z))/(48.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-8 - 2*z + 16*Nc*nf*(-2 + z)*z - 16*pow(Nc,3)*nf*(-2 + z)*z + pow(z,2) - 4*pow(Nc,2)*(-11 + 26*z + 13*pow(z,2)) + pow(Nc,4)*(212 + 354*z + 647*pow(z,2)))*chaplin::HPL(1,0,0,z))/(64.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(-16 + 2*z - pow(z,2) - 4*pow(Nc,2)*(-4 + 36*z + 21*pow(z,2)) + pow(Nc,4)*(352 + 1294*z + 637*pow(z,2)))*chaplin::HPL(1,1,0,z))/(64.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(3 - 26*pow(Nc,2) + 99*pow(Nc,4))*(2 - 2*z + pow(z,2))*pow(log(1 - z),3))/(64.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(1836*z - 675*pow(z,2) + 2232*pow(z,3) + 72*(-30*z + 9*pow(z,2))*z2 - 12*pow(Nc,2)*(-376 + 3783*z + 126*pow(z,2) + 500*pow(z,3) + 24*pow(nf,2)*(2 - 2*z + pow(z,2)) + 36*(16 + 22*z + 39*pow(z,2))*z2) + 3*pow(Nc,4)*(-118504 + 218568*z - 16071*pow(z,2) + 38248*pow(z,3) + 216*(64 - 30*z + 161*pow(z,2))*z2) - 12*pow(Nc,3)*nf*(-3164 + 400*pow(z,3) + z*(2694 - 576*z2) + 3*pow(z,2)*(-669 + 96*z2)) + 12*Nc*nf*(-344 + 112*pow(z,3) + z*(2310 - 576*z2) + pow(z,2)*(-651 + 288*z2)))*log(z))/(13824.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(9*z*(-18 + 29*z + 16*pow(z,2)) - 12*Nc*nf*(-32 + 42*z + 33*pow(z,2)) - 12*pow(Nc,3)*nf*(136 + 62*z + 95*pow(z,2)) + 2*pow(Nc,2)*(2664*z - 1710*pow(z,2) + 240*pow(z,3)) - 3*pow(Nc,4)*(-7168 + 394*z - 2669*pow(z,2) + 1296*pow(z,3)))*pow(log(z),2))/(4608.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(-3*z*(6 + z) - 6*pow(Nc,2)*z*(26 + 31*z) + 2*Nc*nf*(-96*z + 48*pow(z,2)) - 2*pow(Nc,3)*nf*(-96*z + 48*pow(z,2)) + 3*pow(Nc,4)*(224 - 358*z + 703*pow(z,2)))*pow(log(z),3))/(2304.*pow(Nc,3)*z) + pow(log(1 - z),2)*(((-1 + pow(Nc,2))*(-52*Nc*nf*(2 - 2*z + pow(z,2)) + 268*pow(Nc,3)*nf*(2 - 2*z + pow(z,2)) - 3*(36 - 76*z + 31*pow(z,2)) + pow(Nc,2)*(2976 - 3308*z + 754*pow(z,2) - 208*pow(z,3)) + pow(Nc,4)*(-17044 + 14768*z - 481*pow(z,2) + 1712*pow(z,3))))/(768.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(8 - 22*z + 11*pow(z,2) - 8*pow(Nc,2)*(21 - 5*z + 22*pow(z,2)) + pow(Nc,4)*(1040 + 606*z + 981*pow(z,2)))*log(z))/(128.*pow(Nc,3)*z)) + log(1 - z)*(-((-1 + pow(Nc,2))*(2*pow(Nc,3)*nf*(2609 - 3096*z + 1020*pow(z,2) - 224*pow(z,3)) + 2*Nc*nf*(-655 + 1332*z - 762*pow(z,2) + 64*pow(z,3)) + 3*(-23 + 12*z - 24*pow(z,2) + 62*pow(z,3) + 6*(4 - 18*z + 9*pow(z,2))*z2) - pow(Nc,2)*(-6739 + 5754*z + 318*pow(z,2) + 976*pow(z,3) + 24*pow(nf,2)*(2 - 2*z + pow(z,2)) + 36*(62 - 2*z + 77*pow(z,2))*z2) + pow(Nc,4)*(-80200 + 77142*z - 8466*pow(z,2) + 8662*pow(z,3) + 18*(736 + 854*z + 829*pow(z,2))*z2)))/(1152.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-3*pow(Nc,2)*(2 + 2*z + pow(z,2)) + 21*pow(Nc,4)*(2 + 2*z + pow(z,2)))*chaplin::HPL(-1,0,z))/(48.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(-16 + 2*z - pow(z,2) - 4*pow(Nc,2)*(-4 + 36*z + 21*pow(z,2)) + pow(Nc,4)*(352 + 1294*z + 637*pow(z,2)))*chaplin::HPL(1,0,z))/(64.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-36*Nc*nf*(-14 + 6*z + 9*pow(z,2)) - 12*pow(Nc,3)*nf*(154 + 30*z + 69*pow(z,2)) + 3*(-78*z + 144*pow(z,2) + 24*pow(z,3)) + 6*pow(Nc,2)*(-422 + 642*z - 519*pow(z,2) + 72*pow(z,3)) - 6*pow(Nc,4)*(-5722 + 1251*z - 819*pow(z,2) + 660*pow(z,3)))*log(z))/(1152.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(15*(-2 + z)*z + 24*Nc*nf*(-4*z + 2*pow(z,2)) - 8*pow(Nc,3)*nf*(-12*z + 6*pow(z,2)) - 6*pow(Nc,2)*(16 + 14*z + 45*pow(z,2)) + 3*pow(Nc,4)*(480 + 86*z + 781*pow(z,2)))*pow(log(z),2))/(384.*pow(Nc,3)*z));
        check_imaginary_part(L2,__PRETTY_FUNCTION__);
        
        return real(L2);
        
    }

    double LEqgN3LOregFalko_L1(const double& z){
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double z2= consts::z2;
        const double z3= consts::z3;
        const double z4= consts::z4;
        
        complex<double> L1=((-1 + pow(Nc,2))*(4*pow(Nc,3)*nf*(360*(-1782 + 1204*z - 2333*pow(z,2) + 400*pow(z,3))*z2 - 5*(-946310 + 67400*pow(z,3) + z*(1107462 - 91584*z3) + 5184*z3 + 9*pow(z,2)*(-32399 + 1248*z3)) + 207360*(-2 + z)*z*z4) - 12*Nc*nf*(120*(-880 + 3924*z - 465*pow(z,2) + 192*pow(z,3))*z2 + 5*(34394 - 10856*pow(z,3) + pow(z,2)*(114477 - 8064*z3) + 3456*z3 + 6*z*(-22931 + 2112*z3)) + 1080*(68 - 60*z + 81*pow(z,2))*z4) + 3*(360*(-368 - 1366*z + 269*pow(z,2) + 88*pow(z,3))*z2 + 5*(360*(139 + 24*z3) + 384*pow(z,3)*(17 + 108*z3) - 2*z*(69139 + 28080*z3) + pow(z,2)*(72773 + 55296*z3)) - 1080*(160 - 274*z + 347*pow(z,2))*z4) + 4*pow(Nc,2)*(-240*pow(nf,2)*(253 - 406*z + 203*pow(z,2)) + 720*(-3260 + 5607*z - 1632*pow(z,2) + 827*pow(z,3))*z2 - 5*(-383783 + z*(229842 - 181872*z3) + 18576*z3 + 640*pow(z,3)*(365 + 54*z3) + 216*pow(z,2)*(-77 + 1867*z3)) + 3240*(186 + 215*z + 390*pow(z,2))*z4) - 3*pow(Nc,4)*(120*(-56072 + 277222*z - 12185*pow(z,2) + 59360*pow(z,3))*z2 + 5*(8502148 + pow(z,2)*(400589 - 598176*z3) - 148608*z3 + 192*pow(z,3)*(-2203 + 2640*z3) + z*(-8245630 + 2804832*z3)) + 1080*(1720 - 7814*z + 6479*pow(z,2))*z4)))/(829440.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(6*(65 + 135*z + 54*pow(z,2) - 16*pow(z,3) - 20*pow(Nc,3)*nf*(2 + 2*z + pow(z,2)) + pow(Nc,2)*(-131 - 195*z - 33*pow(z,2) + 40*pow(z,3)) + pow(Nc,4)*(2788 + 2590*z + 125*pow(z,2) + 264*pow(z,3)))*z2 + 36*(7 - 16*pow(Nc,2) + 129*pow(Nc,4))*(2 + 2*z + pow(z,2))*z3)*chaplin::HPL(-1,z))/(576.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(6*(-182 + 51*z + 210*pow(z,2) + 56*pow(z,3) - 4*pow(Nc,3)*nf*(114 + 162*z + 19*pow(z,2)) - 8*Nc*nf*(-19 - 4*z + 32*pow(z,2)) + pow(Nc,4)*(-1934 + 10545*z + 2479*pow(z,2) - 904*pow(z,3)) - pow(Nc,2)*(-1816 + 2240*z + 2087*pow(z,2) + 32*pow(z,3)))*z2 + 18*(-60 + 14*z + 64*Nc*nf*(-2 + z)*z - 64*pow(Nc,3)*nf*(-2 + z)*z - 23*pow(z,2) - 4*pow(Nc,2)*(-47 + 54*z + 31*pow(z,2)) + pow(Nc,4)*(40 + 866*z + 1615*pow(z,2)))*z3)*chaplin::HPL(1,z))/(1152.*pow(Nc,3)*z) - ((-8 + 23*pow(Nc,2) - 140*pow(Nc,4) + 125*pow(Nc,6))*(2 + 2*z + pow(z,2))*z2*chaplin::HPL(-1,-1,z))/(16.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(-12*pow(Nc,3)*nf*(35 + 34*z + 9*pow(z,2)) - 4*Nc*nf*(43 + 18*z - 9*pow(z,2) + 16*pow(z,3)) + 3*(99 + 376*z + 289*pow(z,2) - 36*(2 + 2*z + pow(z,2))*z2) + 3*pow(Nc,2)*(29 - 12*z + 147*pow(z,2) + 128*pow(z,3) + 180*(2 + 2*z + pow(z,2))*z2) - 12*pow(Nc,4)*(-786 - 629*z + 149*pow(z,2) + 50*pow(z,3) + 258*(2 + 2*z + pow(z,2))*z2))*chaplin::HPL(-1,0,z))/(576.*pow(Nc,3)*z) - ((-2 + 13*pow(Nc,2) - 94*pow(Nc,4) + 83*pow(Nc,6))*(2 + 2*z + pow(z,2))*z2*chaplin::HPL(-1,1,z))/(16.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(2 + 11*z + 3*pow(z,2) - pow(Nc,2)*(8 + 19*z + 9*pow(z,2)) + pow(Nc,4)*(182 - 62*z + 161*pow(z,2)))*z2*chaplin::HPL(0,-1,z))/(16.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(-32 + 6*z - 64*pow(Nc,3)*nf*(-2 + z)*z - 11*pow(z,2) + 4*Nc*nf*(-4 - 36*z + 15*pow(z,2)) + 2*pow(Nc,2)*(148 - 78*z + 31*pow(z,2)) + pow(Nc,4)*(-856 - 1050*z + 565*pow(z,2)))*z2*chaplin::HPL(0,1,z))/(64.*pow(Nc,3)*z) - ((-2 + 13*pow(Nc,2) - 94*pow(Nc,4) + 83*pow(Nc,6))*(2 + 2*z + pow(z,2))*z2*chaplin::HPL(1,-1,z))/(16.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-1408 + 930*z + 2685*pow(z,2) - 728*pow(z,3) + 1440*z2 + 504*z*z2 - 108*pow(z,2)*z2 + 4*pow(Nc,2)*(6997 - 908*z - 2201*pow(z,2) + 514*pow(z,3) + 72*(-8 + 43*z + 34*pow(z,2))*z2) - pow(Nc,4)*(192696 + 786*z + 25005*pow(z,2) + 26560*pow(z,3) + 36*(280 + 766*z + 1953*pow(z,2))*z2) + 4*pow(Nc,3)*nf*(1010 + 48*pow(z,3) - 12*z*(227 + 96*z2) + 3*pow(z,2)*(-351 + 192*z2)) - 4*Nc*nf*(862 + 88*pow(z,3) + z*(988 - 1152*z2) + pow(z,2)*(1153 + 576*z2)))*chaplin::HPL(1,0,z))/(2304.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-3*(14 - 2*z + pow(z,2)) - 10*pow(Nc,2)*(-2 + 30*z + 21*pow(z,2)) + pow(Nc,4)*(726 + 3062*z + 1469*pow(z,2)))*z2*chaplin::HPL(1,1,z))/(32.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(37 + 87*z + 42*pow(z,2) - 8*pow(z,3) - 4*pow(Nc,3)*nf*(2 + 2*z + pow(z,2)) + pow(Nc,2)*(-47 - 63*z - 9*pow(z,2) + 16*pow(z,3)) + pow(Nc,4)*(856 + 818*z + 13*pow(z,2) + 56*pow(z,3)))*chaplin::HPL(-1,-1,0,z))/(48.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(61 + 147*z + 66*pow(z,2) - 20*pow(z,3) + 4*pow(Nc,3)*nf*(2 + 2*z + pow(z,2)) + 9*pow(Nc,2)*(-11 - 21*z - 5*pow(z,2) + 4*pow(z,3)) + pow(Nc,4)*(2198 + 2008*z - 61*pow(z,2) + 208*pow(z,3)))*chaplin::HPL(-1,0,0,z))/(96.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-2*pow(1 + z,2)*(-7 + 2*z) - 8*pow(Nc,3)*nf*(2 + 2*z + pow(z,2)) - 3*pow(Nc,2)*(14 + 22*z + 4*pow(z,2) - 4*pow(z,3)) + pow(Nc,4)*(966 + 886*z + 56*pow(z,2) + 104*pow(z,3)))*chaplin::HPL(-1,1,0,z))/(48.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(18 + 138*z + 30*pow(z,2) - 16*pow(z,3) + 4*pow(Nc,3)*nf*(2 + 6*z + pow(z,2)) - 4*Nc*nf*(2 - 6*z + 3*pow(z,2)) + pow(Nc,2)*(-2 - 66*z - 27*pow(z,2) + 32*pow(z,3)) + pow(Nc,4)*(540 + 726*z - 97*pow(z,2) + 176*pow(z,3)))*chaplin::HPL(0,-1,0,z))/(96.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(2*Nc*nf*(-16 - 42*z + 21*pow(z,2)) + 2*pow(Nc,3)*nf*(256 + 722*z + 279*pow(z,2)) + 9*(-6 + 21*z + 5*pow(z,2) + 8*pow(z,3)) + pow(Nc,4)*(-9394 - 7411*z - 2316*pow(z,2) + 80*pow(z,3)) - pow(Nc,2)*(-756 + 312*z + 303*pow(z,2) + 136*pow(z,3)))*chaplin::HPL(0,1,0,z))/(192.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-2*pow(1 + z,2)*(-7 + 2*z) - 8*pow(Nc,3)*nf*(2 + 2*z + pow(z,2)) - 3*pow(Nc,2)*(14 + 22*z + 4*pow(z,2) - 4*pow(z,3)) + pow(Nc,4)*(966 + 886*z + 56*pow(z,2) + 104*pow(z,3)))*chaplin::HPL(1,-1,0,z))/(48.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-148 + 120*z + 447*pow(z,2) + 160*pow(z,3) - 4*Nc*nf*(-162 + 98*z + 83*pow(z,2)) + 4*pow(Nc,3)*nf*(166 + 1014*z + 481*pow(z,2)) - 2*pow(Nc,2)*(-640 + 554*z + 1715*pow(z,2) + 40*pow(z,3)) - pow(Nc,4)*(10972 + 15036*z + 7885*pow(z,2) + 1648*pow(z,3)))*chaplin::HPL(1,0,0,z))/(384.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(-154 + 99*z + 222*pow(z,2) + 48*pow(z,3) - 4*pow(Nc,3)*nf*(122 + 170*z + 23*pow(z,2)) - 8*Nc*nf*(-19 - 4*z + 32*pow(z,2)) + pow(Nc,4)*(-2 + 12317*z + 2591*pow(z,2) - 696*pow(z,3)) - pow(Nc,2)*(-1732 + 2372*z + 2111*pow(z,2) + 8*pow(z,3)))*chaplin::HPL(1,1,0,z))/(192.*pow(Nc,3)*z) - ((-4 + 9*pow(Nc,2) - 48*pow(Nc,4) + 43*pow(Nc,6))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,-1,-1,0,z))/(8.*pow(Nc,3)*z) + ((-2 + 6*pow(Nc,2) - 31*pow(Nc,4) + 27*pow(Nc,6))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,-1,0,0,z))/(4.*pow(Nc,3)*z) - ((-2 + 7*pow(Nc,2) - 46*pow(Nc,4) + 41*pow(Nc,6))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,-1,1,0,z))/(8.*pow(Nc,3)*z) + ((-5 + 12*pow(Nc,2) - 55*pow(Nc,4) + 48*pow(Nc,6))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,0,-1,0,z))/(16.*pow(Nc,3)*z) - ((-5 + 26*pow(Nc,2) - 123*pow(Nc,4) + 102*pow(Nc,6))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,0,0,0,z))/(32.*pow(Nc,3)*z) + ((-1 + 5*pow(Nc,2) - 31*pow(Nc,4) + 27*pow(Nc,6))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,0,1,0,z))/(8.*pow(Nc,3)*z) - ((-2 + 7*pow(Nc,2) - 46*pow(Nc,4) + 41*pow(Nc,6))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,1,-1,0,z))/(8.*pow(Nc,3)*z) + ((-1 + 7*pow(Nc,2) - 43*pow(Nc,4) + 37*pow(Nc,6))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,1,0,0,z))/(8.*pow(Nc,3)*z) - (3*(1 - 8*pow(Nc,2) + 7*pow(Nc,4))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,1,1,0,z))/(8.*Nc*z) + ((-1 + pow(Nc,2))*(2 + 7*z + pow(z,2) - pow(Nc,2)*(4 + 7*z + 3*pow(z,2)) + pow(Nc,4)*(58 - 14*z + 49*pow(z,2)))*chaplin::HPL(0,-1,-1,0,z))/(8.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(2 + 13*z + pow(z,2) - pow(Nc,2)*(6 + 21*z + 7*pow(z,2)) + pow(Nc,4)*(152 - 66*z + 139*pow(z,2)))*chaplin::HPL(0,-1,0,0,z))/(16.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(z*(2 + z) - pow(Nc,2)*(2 + 6*z + 3*pow(z,2)) + pow(Nc,4)*(62 - 24*z + 56*pow(z,2)))*chaplin::HPL(0,-1,1,0,z))/(8.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(-3*pow(Nc,2)*z*(2 + z) + z*(10 + z) + pow(Nc,4)*(48 - 70*z + 51*pow(z,2)))*chaplin::HPL(0,0,-1,0,z))/(16.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(8 - 14*z - 16*Nc*nf*(-2 + z)*z + 16*pow(Nc,3)*nf*(-2 + z)*z + 15*pow(z,2) - 8*pow(Nc,2)*(16 - 9*z + 13*pow(z,2)) + pow(Nc,4)*(824 + 354*z + 591*pow(z,2)))*chaplin::HPL(0,0,1,0,z))/(64.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(z*(2 + z) - pow(Nc,2)*(2 + 6*z + 3*pow(z,2)) + pow(Nc,4)*(62 - 24*z + 56*pow(z,2)))*chaplin::HPL(0,1,-1,0,z))/(8.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-16 - 2*z - 64*pow(Nc,3)*nf*(-2 + z)*z - 5*pow(z,2) + 6*pow(Nc,2)*(36 - 32*z + 5*pow(z,2)) + 4*Nc*nf*(4 - 28*z + 17*pow(z,2)) + pow(Nc,4)*(-536 - 1070*z + 623*pow(z,2)))*chaplin::HPL(0,1,0,0,z))/(64.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(-32 + 14*z - 64*pow(Nc,3)*nf*(-2 + z)*z - 7*pow(z,2) + 4*Nc*nf*(-4 - 36*z + 15*pow(z,2)) + 2*pow(Nc,2)*(144 - 90*z + 25*pow(z,2)) + pow(Nc,4)*(-608 - 1146*z + 789*pow(z,2)))*chaplin::HPL(0,1,1,0,z))/(64.*pow(Nc,3)*z) - ((-2 + 7*pow(Nc,2) - 46*pow(Nc,4) + 41*pow(Nc,6))*(2 + 2*z + pow(z,2))*chaplin::HPL(1,-1,-1,0,z))/(8.*pow(Nc,3)*z) + ((-1 + 7*pow(Nc,2) - 43*pow(Nc,4) + 37*pow(Nc,6))*(2 + 2*z + pow(z,2))*chaplin::HPL(1,-1,0,0,z))/(8.*pow(Nc,3)*z) - (3*(1 - 8*pow(Nc,2) + 7*pow(Nc,4))*(2 + 2*z + pow(z,2))*chaplin::HPL(1,-1,1,0,z))/(8.*Nc*z) + ((-1 + pow(Nc,2))*(6 - 2*z + 5*pow(z,2) - pow(Nc,2)*(6 + 10*z + 7*pow(z,2)) + 4*pow(Nc,4)*(31 - 12*z + 28*pow(z,2)))*chaplin::HPL(1,0,-1,0,z))/(16.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-4 - 14*z + 32*Nc*nf*(-2 + z)*z - 32*pow(Nc,3)*nf*(-2 + z)*z + 3*pow(z,2) - 2*pow(Nc,2)*(-23 + 66*z + 43*pow(z,2)) + 2*pow(Nc,4)*(73 + 132*z + 427*pow(z,2)))*chaplin::HPL(1,0,0,0,z))/(32.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(7*(-2 + z)*z + 64*Nc*nf*(-2 + z)*z - 64*pow(Nc,3)*nf*(-2 + z)*z - 4*pow(Nc,2)*(-34 + 6*z + 7*pow(z,2)) + pow(Nc,4)*(-8 - 1626*z + 1125*pow(z,2)))*chaplin::HPL(1,0,1,0,z))/(64.*pow(Nc,3)*z) - (3*(1 - 8*pow(Nc,2) + 7*pow(Nc,4))*(2 + 2*z + pow(z,2))*chaplin::HPL(1,1,-1,0,z))/(8.*Nc*z) - ((-1 + pow(Nc,2))*(-52 - 10*z + 64*Nc*nf*(-2 + z)*z - 64*pow(Nc,3)*nf*(-2 + z)*z + 5*pow(z,2) - 2*pow(Nc,2)*(-52 + 214*z + 169*pow(z,2)) + pow(Nc,4)*(972 + 2614*z + 3061*pow(z,2)))*chaplin::HPL(1,1,0,0,z))/(64.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-3*(14 - 2*z + pow(z,2)) - 8*pow(Nc,2)*(-1 + 39*z + 27*pow(z,2)) + pow(Nc,4)*(810 + 3146*z + 1511*pow(z,2)))*chaplin::HPL(1,1,1,0,z))/(32.*pow(Nc,3)*z) - ((-27 + 221*pow(Nc,2) - 933*pow(Nc,4) + 739*pow(Nc,6))*(2 - 2*z + pow(z,2))*pow(log(1 - z),4))/(384.*pow(Nc,3)*z) + (((-1 + pow(Nc,2))*(2808 - 5898*z + 6618*pow(z,2) + 15112*pow(z,3) + 72*(-36 - 147*z + 63*pow(z,2) + 56*pow(z,3))*z2 + 1728*z3 - 9072*z*z3 - 2376*pow(z,2)*z3 + 4*Nc*nf*(-3622 + 25498*z - 23387*pow(z,2) + 600*pow(z,3) - 36*(-60 + 86*z + 17*pow(z,2))*z2 - 1728*z3 - 5184*z*z3 + 1296*pow(z,2)*z3) - 4*pow(Nc,3)*nf*(-42330 + 54350*z - 26467*pow(z,2) + 7672*pow(z,3) + 36*(76 - 62*z + 51*pow(z,2))*z2 - 3456*z*z3 + 1728*pow(z,2)*z3) - 2*pow(Nc,2)*(-26482 + 144308*z - 2071*pow(z,2) + 42820*pow(z,3) + 48*pow(nf,2)*(29 - 38*z + 19*pow(z,2)) + 30240*z2 - 54576*z*z2 + 63144*pow(z,2)*z2 - 2592*pow(z,3)*z2 + 864*z3 - 4752*z*z3 + 7992*pow(z,2)*z3) + 2*pow(Nc,4)*(-662738 + 125272*pow(z,3) - 36*(-2172 + 3727*z - 1521*pow(z,2) + 1200*pow(z,3))*z2 + 10368*z3 + 2*pow(z,2)*(-84451 + 32670*z3) + z*(1196185 + 56808*z3))))/(13824.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(11*Nc - 2*nf)*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,0,z))/(12.*z) + ((-1 + pow(Nc,2))*(11*Nc - 2*nf)*(1 + pow(Nc,2)*(5 + 18*z + 9*pow(z,2)))*chaplin::HPL(1,0,z))/(6.*pow(Nc,2)*z))*log(z) - ((-1 + pow(Nc,2))*(-2*z*(1095 + 765*z - 8*pow(z,2) + 18*(-34 + 5*z)*z2) + 2*pow(Nc,2)*(-752 + 10619*z - 2401*pow(z,2) - 720*pow(z,3) + 48*pow(nf,2)*(2 - 2*z + pow(z,2)) + 288*(4 + 8*z + 11*pow(z,2))*z2) - 2*pow(Nc,4)*(2*(-29218 + 79094*z - 6463*pow(z,2) + 12236*pow(z,3)) + 18*(256 - 218*z + 977*pow(z,2))*z2) - 4*pow(Nc,3)*nf*(3164 - 512*pow(z,3) + pow(z,2)*(2374 - 432*z2) + z*(-3845 + 864*z2)) + 4*Nc*nf*(344 - 120*pow(z,3) + pow(z,2)*(740 - 432*z2) + z*(-2941 + 864*z2)))*pow(log(z),2))/(4608.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(2*Nc*nf*(64 - 46*z - 49*pow(z,2)) + pow(Nc,2)*z*(2372 - 1933*z + 16*pow(z,2)) + z*(39 - 3*z + 176*pow(z,2)) - 2*pow(Nc,3)*nf*(272 + 154*z + 227*pow(z,2)) + pow(Nc,4)*(7168 - 1795*z + 4420*pow(z,2) - 1920*pow(z,3)))*pow(log(z),3))/(2304.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(88*Nc*nf*(-2 + z)*z - 88*pow(Nc,3)*nf*(-2 + z)*z - 22*pow(Nc,2)*z*(6 + 7*z) - z*(22 + 9*z) + pow(Nc,4)*(448 - 902*z + 1659*pow(z,2)))*pow(log(z),4))/(3072.*pow(Nc,3)*z) + pow(log(1 - z),3)*(((-1 + pow(Nc,2))*(50*Nc*nf*(2 - 2*z + pow(z,2)) - 198*pow(Nc,3)*nf*(2 - 2*z + pow(z,2)) + 27*(9 - 16*z + 7*pow(z,2)) + pow(Nc,2)*(-4038 + 4558*z - 1241*pow(z,2) + 248*pow(z,3)) - pow(Nc,4)*(-19543 + 16818*z + 30*pow(z,2) + 2056*pow(z,3))))/(576.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(9*pow(-1 + z,2) - 2*pow(Nc,2)*(61 - 23*z + 58*pow(z,2)) + pow(Nc,4)*(673 + 316*z + 613*pow(z,2)))*log(z))/(48.*pow(Nc,3)*z)) + pow(log(1 - z),2)*(-((-1 + pow(Nc,2))*(1507 - 2538*z + 1413*pow(z,2) - 868*pow(z,3) + 4*Nc*nf*(1415 - 2772*z + 1584*pow(z,2) - 128*pow(z,3)) + 4*pow(Nc,3)*nf*(-4585 + 5648*z - 2122*pow(z,2) + 368*pow(z,3)) - 792*z2 + 2088*z*z2 - 1044*pow(z,2)*z2 + 2*pow(Nc,2)*(-23518 + 21414*z - 957*pow(z,2) + 3232*pow(z,3) + 24*pow(nf,2)*(2 - 2*z + pow(z,2)) + 72*(84 - 14*z + 97*pow(z,2))*z2) - pow(Nc,4)*(-343449 + 322538*z - 28057*pow(z,2) + 34036*pow(z,3) + 36*(1658 + 2130*z + 1935*pow(z,2))*z2)))/(2304.*pow(Nc,3)*z) - (3*(1 - 8*pow(Nc,2) + 7*pow(Nc,4))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,0,z))/(16.*Nc*z) + ((-1 + pow(Nc,2))*(-34 - 2*z + pow(z,2) + pow(Nc,2)*(4 - 308*z - 218*pow(z,2)) + pow(Nc,4)*(806 + 3150*z + 1509*pow(z,2)))*chaplin::HPL(1,0,z))/(64.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(108 - 399*z + 432*pow(z,2) + 56*pow(z,3) - 4*Nc*nf*(-90 + 38*z + 53*pow(z,2)) - 4*pow(Nc,3)*nf*(266 + 26*z + 95*pow(z,2)) + pow(Nc,4)*(26024 - 9859*z + 2432*pow(z,2) - 3448*pow(z,3)) + pow(Nc,2)*(-3360 + 3290*z - 3472*pow(z,2) + 304*pow(z,3)))*log(z))/(384.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(16 - 78*z + 64*Nc*nf*(-2 + z)*z - 64*pow(Nc,3)*nf*(-2 + z)*z + 39*pow(z,2) - 6*pow(Nc,2)*(56 - 2*z + 93*pow(z,2)) + pow(Nc,4)*(2864 + 722*z + 4007*pow(z,2)))*pow(log(z),2))/(256.*pow(Nc,3)*z)) + log(1 - z)*(-((-1 + pow(Nc,2))*(956 - 5496*z + 834*pow(z,2) + 3706*pow(z,3) + 612*z2 - 5670*z*z2 + 6480*pow(z,2)*z2 + 1008*pow(z,3)*z2 - 864*z3 - 1620*z*z3 - 54*pow(z,2)*z3 + pow(Nc,2)*(99775 - 81182*z + 4423*pow(z,2) - 25160*pow(z,3) - 24*pow(nf,2)*(29 - 38*z + 19*pow(z,2)) + 6*(-6930 + 6486*z - 9381*pow(z,2) + 768*pow(z,3))*z2 - 4968*z3 + 3456*z*z3 - 14256*pow(z,2)*z3) - 2*Nc*nf*(5633 - 688*pow(z,3) + 144*(-21 + 8*z + 14*pow(z,2))*z2 + pow(z,2)*(15449 - 1728*z3) + 4*z*(-5089 + 864*z3)) - 2*pow(Nc,3)*nf*(-37739 + 3888*pow(z,3) + 36*(158 + 118*z + 41*pow(z,2))*z2 + z*(46880 - 3456*z3) + pow(z,2)*(-15313 + 1728*z3)) + pow(Nc,4)*(-780899 + 50374*pow(z,3) + 6*(34548 - 1581*z + 5883*pow(z,2) - 7272*pow(z,3))*z2 + z*(791734 - 10044*z3) + 58968*z3 + pow(z,2)*(-74957 + 115614*z3))))/(3456.*pow(Nc,3)*z) - ((-2 + 13*pow(Nc,2) - 94*pow(Nc,4) + 83*pow(Nc,6))*(2 + 2*z + pow(z,2))*z2*chaplin::HPL(-1,z))/(16.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-18 - 5*pow(Nc,2)*(-2 + 30*z + 21*pow(z,2)) + pow(Nc,4)*(356 + 1538*z + 731*pow(z,2)))*z2*chaplin::HPL(1,z))/(16.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-2*pow(1 + z,2)*(-7 + 2*z) - 8*pow(Nc,3)*nf*(2 + 2*z + pow(z,2)) - 3*pow(Nc,2)*(14 + 22*z + 4*pow(z,2) - 4*pow(z,3)) + pow(Nc,4)*(966 + 886*z + 56*pow(z,2) + 104*pow(z,3)))*chaplin::HPL(-1,0,z))/(48.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(-190 + 147*z + 201*pow(z,2) + 48*pow(z,3) - 4*pow(Nc,3)*nf*(126 + 166*z + 25*pow(z,2)) - 4*Nc*nf*(-58 + 6*z + 69*pow(z,2)) + pow(Nc,4)*(438 + 12013*z + 2545*pow(z,2) - 760*pow(z,3)) - 2*pow(Nc,2)*(-862 + 1176*z + 1059*pow(z,2) + 4*pow(z,3)))*chaplin::HPL(1,0,z))/(192.*pow(Nc,3)*z) - ((-2 + 7*pow(Nc,2) - 46*pow(Nc,4) + 41*pow(Nc,6))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,-1,0,z))/(8.*pow(Nc,3)*z) + ((-1 + 7*pow(Nc,2) - 43*pow(Nc,4) + 37*pow(Nc,6))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,0,0,z))/(8.*pow(Nc,3)*z) - (3*(1 - 8*pow(Nc,2) + 7*pow(Nc,4))*(2 + 2*z + pow(z,2))*chaplin::HPL(-1,1,0,z))/(8.*Nc*z) + ((-1 + pow(Nc,2))*(z*(2 + z) - pow(Nc,2)*(2 + 6*z + 3*pow(z,2)) + pow(Nc,4)*(62 - 24*z + 56*pow(z,2)))*chaplin::HPL(0,-1,0,z))/(8.*pow(Nc,3)*z) - ((-1 + pow(Nc,2))*(-32 + 18*z + 64*Nc*nf*(-2 + z)*z - 64*pow(Nc,3)*nf*(-2 + z)*z - 9*pow(z,2) + 16*pow(Nc,2)*(18 - 11*z + 3*pow(z,2)) + pow(Nc,4)*(-544 - 1090*z + 857*pow(z,2)))*chaplin::HPL(0,1,0,z))/(64.*pow(Nc,3)*z) - (3*(1 - 8*pow(Nc,2) + 7*pow(Nc,4))*(2 + 2*z + pow(z,2))*chaplin::HPL(1,-1,0,z))/(8.*Nc*z) - ((-1 + pow(Nc,2))*(-40 - 22*z + 64*Nc*nf*(-2 + z)*z - 64*pow(Nc,3)*nf*(-2 + z)*z + 11*pow(z,2) - 2*pow(Nc,2)*(-56 + 218*z + 167*pow(z,2)) + pow(Nc,4)*(1064 + 2522*z + 3107*pow(z,2)))*chaplin::HPL(1,0,0,z))/(64.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-9 + pow(Nc,2)*(2 - 78*z - 54*pow(z,2)) + pow(Nc,4)*(199 + 790*z + 376*pow(z,2)))*chaplin::HPL(1,1,0,z))/(8.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(-2*pow(Nc,2)*(7021 - 11751*z + 2271*pow(z,2) - 2008*pow(z,3) - 24*pow(nf,2)*(2 - 2*z + pow(z,2)) - 1944*z2 - 576*z*z2 - 3708*pow(z,2)*z2) + 3*(46 - 438*z + 697*pow(z,2) - 204*pow(z,3) + 6*(-8 + 62*z - 27*pow(z,2))*z2) + pow(Nc,4)*(146108 - 233036*z + 11551*pow(z,2) - 37844*pow(z,3) - 18*(1040 + 6*z + 2333*pow(z,2))*z2) - 4*pow(Nc,3)*nf*(3025 - 336*pow(z,3) - 8*pow(z,2)*(-239 + 36*z2) + 4*z*(-695 + 144*z2)) - 4*Nc*nf*(-503 + 128*pow(z,3) - 96*z*(-26 + 6*z2) + pow(z,2)*(-645 + 288*z2)))*log(z))/(1152.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(24*Nc*nf*(7 - 8*z - 5*pow(z,2)) + 3*z*(-26 + 89*z + 32*pow(z,2)) - 8*pow(Nc,3)*nf*(77 - 4*z + 43*pow(z,2)) + pow(Nc,4)*(10212 - 5312*z + 2450*pow(z,2) - 2056*pow(z,3)) - pow(Nc,2)*(844 - 1920*z + 2223*pow(z,2) - 88*pow(z,3)))*pow(log(z),2))/(384.*pow(Nc,3)*z) + ((-1 + pow(Nc,2))*(32*Nc*nf*(-2 + z)*z - 32*pow(Nc,3)*nf*(-2 + z)*z + z*(-18 + 5*z) - pow(Nc,2)*(32 + 54*z + 125*pow(z,2)) + pow(Nc,4)*(480 - 70*z + 1021*pow(z,2)))*pow(log(z),3))/(192.*pow(Nc,3)*z));
        check_imaginary_part(L1,__PRETTY_FUNCTION__);
        
        return real(L1);
        
    }
    
    double LEqgN3LOregFalko(const double& z, const double& L)
    {
//        //const double nf = consts::nf;
//        const double Nc = QCD::Nc;
//        const double z2= consts::z2;
//        const double z3= consts::z3;
//        const double z4= consts::z4;
        if (abs(L)<1e-15) return 0.0;
        
        double L3=LEqgN3LOregFalko_L3(z);
        double L2=LEqgN3LOregFalko_L2(z);
        double L1=LEqgN3LOregFalko_L1(z);
        
        return L3*pow(L,3)+L2*pow(L,2)+L1*L;
    }
    
    
    
//    double LEqgNNNLOregSte(const double& z, const double& L)
//    {
//        if (abs(L)<1e-15) return 0.0;
//        
//        double zz = z;
//        double LL = L;
//        return leqgnnnloreg_(&zz,&LL);
//    }
    
    //---- Next to soft approximation (i.e. constants only in front of logs)
    double qg_n3lo_r_lz5_NS(const double& z, const double& L)
    {
        const double Nc = QCD::Nc;
        
        complex<double> res =
        (-27 + 181*pow(Nc,2) - 741*pow(Nc,4) + 587*pow(Nc,6))/(768.*pow(Nc,3))
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(1.-z),5.);
    }
    
    double qg_n3lo_r_lz4_NS(const double& z, const double& L)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        
        complex<double> res =
        -((-1 + pow(Nc,2))*(687 - 15118*pow(Nc,2) + 9155*pow(Nc,4) + 72*L*(27 - 194*pow(Nc,2) + 739*pow(Nc,4)) + 1012*Nc*nf - 3212*pow(Nc,3)*nf))/(27648.*pow(Nc,3))
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(1.-z),4.);
    }
    
    double qg_n3lo_r_lz3_NS(const double& z, const double& L)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double Z2 = consts::z2;
        complex<double> res =
        ((-1 + pow(Nc,2))*(648*pow(L,2)*(3 - 26*pow(Nc,2) + 99*pow(Nc,4)) + 5264*Nc*nf - 29392*pow(Nc,3)*nf + 72*L*Nc*(-473*Nc + 639*pow(Nc,3) + 50*nf - 198*pow(Nc,2)*nf) + pow(Nc,4)*(166903 - 200952*Z2) - 3*(1899 + 4632*Z2) + 2*pow(Nc,2)*(-3085 + 168*pow(nf,2) + 37728*Z2)))/(41472.*pow(Nc,3))
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(1.-z),3.);
    }
    
    double qg_n3lo_r_lz2_NS(const double& z, const double& L)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double Z2 = consts::z2;
        const double Z3 = consts::z3;
        complex<double> res =
        ((-1 + pow(Nc,2))*(-432*pow(L,3)*(1 - 10*pow(Nc,2) + 37*pow(Nc,4)) - 54*pow(L,2)*(-27 - 214*pow(Nc,2) + 1045*pow(Nc,4) + 52*Nc*nf - 268*pow(Nc,3)*nf) + 4*pow(Nc,3)*nf*(6427 - 6660*Z2) + 4*Nc*nf*(-1313 + 1836*Z2) + 36*L*(243 - 198*Nc*nf + 1382*pow(Nc,3)*nf + 504*Z2 - 3*pow(Nc,2)*(57 + 8*pow(nf,2) + 1056*Z2) + 2*pow(Nc,4)*(-3733 + 4212*Z2)) + 9*(435 + 1392*Z2 + 2472*Z3) - 2*pow(Nc,2)*(-9001 + 528*pow(nf,2) + 52164*Z2 + 93636*Z3) + pow(Nc,4)*(-120073 + 124488*Z2 + 728784*Z3)))/(41472.*pow(Nc,3))
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*pow(log(1.-z),2.);
    }
    
    double qg_n3lo_r_lz1_NS(const double& z, const double& L)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double Z2 = consts::z2;
        const double Z3 = consts::z3;
        const double Z4 = consts::z4;
        complex<double> res =
        -((-1 + pow(Nc,2))*(-2160*pow(L,3)*(-9 - pow(Nc,2) + 296*pow(Nc,4) + 10*Nc*nf - 62*pow(Nc,3)*nf) + 3240*pow(L,2)*(27 - 14*Nc*nf + 206*pow(Nc,3)*nf + 60*Z2 - pow(Nc,2)*(103 + 8*pow(nf,2) + 432*Z2) + pow(Nc,4)*(-954 + 1236*Z2)) + 20*Nc*nf*(-20051 + 34632*Z2 + 41040*Z3) - 20*pow(Nc,3)*nf*(-157411 + 33480*Z2 + 54000*Z3) + 360*L*(-8*pow(Nc,3)*nf*(-571 + 180*Z2) + 4*Nc*nf*(-19 + 234*Z2) + 27*(6*Z2 + 72*Z3) - pow(Nc,2)*(2144 + 240*pow(nf,2) + 8658*Z2 + 15228*Z3) + 2*pow(Nc,4)*(-6874 + 2826*Z2 + 30294*Z3)) - 9*(19395 + 89640*Z2 - 93600*Z3 + 32760*Z4) - 4*pow(Nc,2)*(469775 + 20880*pow(nf,2) + 835200*Z2 + 1258200*Z3 + 245430*Z4) + pow(Nc,4)*(-8205065 + 4055400*Z2 + 4890240*Z3 + 11288160*Z4)))/(1.24416e6*pow(Nc,3))
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res)*log(1.-z);
    }
    
    
    
    double qg_n3lo_r_lz0_NS(const double& z, const double& L)
    {
        const double nf = consts::nf;
        const double Nc = QCD::Nc;
        const double Z2 = consts::z2;
        const double Z3 = consts::z3;
        const double Z4 = consts::z4;
        const double Z5 = consts::z5;
        
        complex<double> res =
        ((-1 + pow(Nc,2))*(1080*pow(L,3)*(-27 + 60*Nc*nf + 292*pow(Nc,3)*nf + 48*Z2 - 4*pow(Nc,2)*(69 + 8*pow(nf,2) + 120*Z2) + pow(Nc,4)*(-665 + 1776*Z2)) - 540*pow(L,2)*(-36*Nc*nf*(13 + 24*Z2) + 4*pow(Nc,3)*nf*(-1327 + 96*Z2) + pow(Nc,4)*(13385 + 1632*Z2 - 31104*Z3) + 27*(17 + 24*Z2 - 48*Z3) + 4*pow(Nc,2)*(1143 + 80*pow(nf,2) + 1242*Z2 + 1980*Z3)) + 12*Nc*nf*(79555 + 115200*Z2 + 105120*Z3 + 74160*Z4) + 4*pow(Nc,3)*nf*(410855 - 533760*Z2 - 2165760*Z3 + 416880*Z4) + 18*L*(20*Nc*nf*(-429 + 4416*Z2 + 5904*Z3) - 20*pow(Nc,3)*nf*(-21013 + 6240*Z2 + 9360*Z3) - 9*(4965 + 13680*Z2 - 7920*Z3 + 4680*Z4) - 4*pow(Nc,2)*(105045 + 4000*pow(nf,2) + 105600*Z2 + 92160*Z3 + 38610*Z4) + pow(Nc,4)*(-1170655 + 708720*Z2 + 867600*Z3 + 1218240*Z4)) + pow(Nc,4)*(-7287205 + 49128480*Z3 - 3840*Z2*(-3691 + 13635*Z3) - 1401840*Z4 + 87454080*Z5) - 9*(2880*Z2*(64 + 93*Z3) + 247680*Z4 - 45*(-1699 + 2560*Z3 + 7872*Z5)) - 4*pow(Nc,2)*(320*pow(nf,2)*(125 + 54*Z3) + 3*(585085 - 253080*Z3 - 720*Z2*(-952 + 2361*Z3) + 589860*Z4 + 2088720*Z5))))/(4.97664e6*pow(Nc,3))
        ;
        check_imaginary_part(res,__PRETTY_FUNCTION__);
        return real(res);
    }
    
    //series expansions of the coeffs of log(1-z)^a
    
    
    
    double qg_n3lo_r_lz5_full_series(const double&z, const double& L){
        return qg_n3lo_r_lzX_series(5,z,37);
    }
    double qg_n3lo_r_lz4_full_series(const double&z, const double& L){
        return qg_n3lo_r_lzX_series(4,z,37);
    }
    double qg_n3lo_r_lz3_full_series(const double&z, const double& L){
        return qg_n3lo_r_lzX_series(3,z,37);
    }
    double qg_n3lo_r_lz2_full_series(const double&z, const double& L){
        return qg_n3lo_r_lzX_series(2,z,37);
    }
    double qg_n3lo_r_lz1_full_series(const double&z, const double& L){
        return qg_n3lo_r_lzX_series(1,z,37);
    }
    double qg_n3lo_r_lz0_full_series(const double&z, const double& L){
        return qg_n3lo_r_lzX_series(0,z,37);
    }
    
    class qgregn3locoeffs{
    public:
        static double give(int logpow,int coeff){
         double coeffs_lzbar_5[38]={17.81944444444444,17.81944444444444,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889,35.63888888888889};
         double coeffs_lzbar_4[38]={-2.011959876543210,-79.33860596707819,-118.5869341563786,-66.30692729766804,-17.68167009602195,17.41584362139918,45.60543552812071,69.40403929061336,90.09070032333921,108.4275752825136,124.9149466816905,139.9029117230506,153.6478503064151,166.3440088289625,178.1425219584942,189.1635303396877,199.5042471974185,209.2445216971361,218.4507894322823,227.1789471996127,235.4764891941739,243.3841233773685,250.9370140566481,258.1657506016084,265.0971121573158,271.7546781351294,278.1593205588922,284.3296048143410,290.2821186071326,296.0317440903770,301.5918845925833,306.9746547716116,312.1910410746560,317.2510379155762,322.1637638610329,326.9375612550373,331.5800820425018,336.0983620288373};
        
         double coeffs_lzbar_3[38]={-114.5690210474115,421.4551559073210,381.4115978228725,237.3542132778336,59.06447562351261,-10.79926580309919,-45.24159719313120,-61.90769360700539,-67.89625556985637,-66.87000822611626,-60.99269181910077,-51.65049894096701,-39.78009611993427,-26.03877104202217,-10.90103385636832,5.282870518373797,22.24976190753040,39.79896467797494,57.77559276427757,76.05888444251552,94.55388643235810,113.1854162926278,131.8936084730658,150.6305821036629,169.3579166259137,188.0447178587626,206.6661213532181,225.2021235095870,243.6366610482722,261.9568805404984,280.1525547172483,298.2156130839379,316.1397622428796,333.9201771242113,351.5532486387011,369.0363765032138,386.3677984408155,403.5464488283852};
        
         double coeffs_lzbar_2[38]={513.5629803085088,-754.7879284323460,-280.9749383083636,-2.010140598903302,503.5296667521548,627.8999113675738,691.4555218612026,733.6075279128865,765.1478777982313,790.6630832471588,812.5754745942297,832.3061986213135,850.7348120736729,868.4218379674948,885.7301021009448,902.8958826319175,920.0726192466645,937.3586556546297,954.8152782284386,972.4786681309675,990.3679409360360,1008.490620992915,1026.846405838890,1045.429774655680,1064.231805765861,1083.241447237601,1102.446405950780,1121.833768460228,1141.390432149639,1161.103401545714,1180.959988465812,1200.947943455886,1221.055538135672,1241.271612542261,1261.585597638517,1281.987520345997,1302.467996442352,1323.018215201810};
        
         double coeffs_lzbar_1[38]={-313.9852299862526,807.2802083194414,673.0163240213602,424.9243723359774,-94.52326000972149,-16.19766747490587,53.68992041930396,107.8211489777683,152.2019089557034,190.1122694361395,223.2479850946826,252.5941586129421,278.8051717983267,302.3631953424077,323.6479529990322,342.9701714491425,360.5895957954527,376.7259919756160,391.5666722532899,405.2720903817970,417.9802283806915,429.8101426132216,440.8648788640696,451.2338909801923,460.9950585290657,470.2163760394976,478.9573713023089,487.2702992574010,495.2011495320082,502.7904989261898,510.0742346362070,517.0841694972283,523.8485668248953,530.3925893947198,536.7386846010706,542.9069157864565,548.9152480465371,554.7797954304352};

         
            double coeffs_lzbar_0[38]={204.6207897981724,94.71170901366414,-336.5212650815272,51.21499932306299,240.5837852155830,132.4535263918228,96.83252991517496,88.26348804400830,90.47571630152621,97.70184465966271,107.5995615353205,119.0633675796899,131.4960886021452,144.5366204177273,157.9476790156004,171.5643241707085,185.2678293389744,198.9710043018505,212.6091867574482,226.1342999583811,239.5107192630002,252.7122881830407,265.7201093557246,278.5208809287844,291.1056286335201,303.4687309956332,315.6071649306303,327.5199188126420,339.2075338596950,350.6717444904992,361.9151954515531,372.9412187905186,383.7536576922293,394.3567271606878,404.7549037798985,414.9528385021117,424.9552877277702,434.7670589554991};
        
            switch (logpow){
                case 5: return coeffs_lzbar_5[coeff];break;
                case 4: return coeffs_lzbar_4[coeff];break;
                case 3: return coeffs_lzbar_3[coeff];break;
                case 2: return coeffs_lzbar_2[coeff];break;
                case 1: return coeffs_lzbar_1[coeff];break;
                case 0: return coeffs_lzbar_0[coeff];break;
                default: cout<<"\nerror: in qg_n3lo_r_lzX, logpow is not in 0-5, logpow="<<logpow<<endl;
                    exit(EXIT_FAILURE);
  
            }
        }
    };
    
    double qg_n3lo_r_lzX_series(int logzbpow,const double& z, int m){
        const int trunc_order = 37;
        if (m>trunc_order){
            cout<<"\nError: the qg_n3lo coefficients are truncated to O("<<trunc_order<<" and you asked for O("<<m<<"))"<<endl;
            exit(EXIT_FAILURE);
        }
        double res=0.;
        for (int i=0;i<m+1;i++){
            res += qgregn3locoeffs::give(logzbpow,i)*intpow(1.-z,i);
        }
        return res*intpow(log(1.-z),logzbpow);
    }
    
   
    
    
}
