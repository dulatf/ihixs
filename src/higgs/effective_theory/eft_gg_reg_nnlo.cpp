#include "constants.h"
#include<vector>
#include "as_series.h"
#include "higgs_eft.h"
#include "cppchaplin.h"

#include <complex>
using namespace std;

namespace HEFT{
// nnlo
double nnlo_reg(const double& z, const double& L)
{
    return    nnlo_r_lz0(z,L)
    + nnlo_r_lz1(z,L)
    + nnlo_r_lz2(z,L)
    + nnlo_r_lz3(z,L) ;
}


double nnlo_r_lz0(const double& z, const double& L)
{
    return    nnlo_r_lz0_const(z,L)
    + nnlo_r_lz0_logz(z,L)
    + nnlo_r_lz0_logz_sq(z,L)
    + nnlo_r_lz0_logz_cube(z,L) ;
}

double nnlo_r_lz0_const(const double& z, const double& L)
{
    const double Nc = QCD::Nc;
    complex<double> res =  (pow(L,2)*pow(Nc,2)*(-99 + 94*z - 83*pow(z,2) + 99*pow(z,3)))/(12.*z) + consts::z2*(-((L*pow(Nc,2)*(-4 - 9*z - 24*pow(z,2) - 14*pow(z,3) + 6*pow(z,4)))/(z*(1 + z))) + (pow(Nc,2)*(-60 + 184*z - 320*pow(z,2) + 322*pow(z,3) + 391*pow(z,4) - 495*pow(z,5) + 6*(-11 - 11*z - 11*pow(z,2) + 17*pow(z,3) + 8*pow(z,4) + 8*pow(z,5))*chaplin::HPL(-1,z) + 6*(-7 + 127*z + 121*pow(z,2) - 125*pow(z,3) - 119*pow(z,4) + pow(z,5))*chaplin::HPL(1,z)))/(12.*z*(-1 + pow(z,2)))) - (L*pow(Nc,2)*(144*(-1 + z)*pow(1 + z + pow(z,2),2)*chaplin::HPL(-1,0,z) + (1 + z)*(-2161 + 3942*z - 3294*pow(z,2) + 3674*pow(z,3) - 2161*pow(z,4) + 288*(3 - 2*z + 9*pow(z,2) - 10*pow(z,3) + 3*pow(z,4))*chaplin::HPL(0,1,z) + 288*(3 + 2*z + 9*pow(z,2) - 14*pow(z,3) + 3*pow(z,4))*chaplin::HPL(1,0,z))))/(72.*z*(-1 + pow(z,2))) + (pow(Nc,2)*(18321 - 18523*z + 1430*pow(z,2) - 2642*pow(z,3) - 19751*pow(z,4) + 21165*pow(z,5) + 72*(14 + 12*z - 35*pow(z,2) - 23*pow(z,3) + 21*pow(z,4) + 11*pow(z,5))*chaplin::HPL(-1,0,z) - 144*(94 - 68*z + 165*pow(z,2) - 75*pow(z,3) - 242*pow(z,4) + 160*pow(z,5))*chaplin::HPL(0,1,z) - 3132*chaplin::HPL(1,0,z) + 6012*z*chaplin::HPL(1,0,z) - 38700*pow(z,2)*chaplin::HPL(1,0,z) + 20124*pow(z,3)*chaplin::HPL(1,0,z) + 39780*pow(z,4)*chaplin::HPL(1,0,z) - 28188*pow(z,5)*chaplin::HPL(1,0,z) - 1296*chaplin::HPL(-1,-1,0,z) - 1296*z*chaplin::HPL(-1,-1,0,z) - 1296*pow(z,2)*chaplin::HPL(-1,-1,0,z) + 3888*pow(z,3)*chaplin::HPL(-1,-1,0,z) + 216*chaplin::HPL(-1,0,0,z) + 216*z*chaplin::HPL(-1,0,0,z) + 216*pow(z,2)*chaplin::HPL(-1,0,0,z) - 2376*pow(z,3)*chaplin::HPL(-1,0,0,z) + 864*pow(z,4)*chaplin::HPL(-1,0,0,z) + 864*pow(z,5)*chaplin::HPL(-1,0,0,z) - 1728*chaplin::HPL(-1,1,0,z) - 1728*z*chaplin::HPL(-1,1,0,z) - 1728*pow(z,2)*chaplin::HPL(-1,1,0,z) + 1728*pow(z,3)*chaplin::HPL(-1,1,0,z) + 1728*pow(z,4)*chaplin::HPL(-1,1,0,z) + 1728*pow(z,5)*chaplin::HPL(-1,1,0,z) + 432*pow(z,4)*chaplin::HPL(0,-1,0,z) - 432*pow(z,5)*chaplin::HPL(0,-1,0,z) - 5184*chaplin::HPL(0,0,1,z) - 1728*z*chaplin::HPL(0,0,1,z) - 12096*pow(z,2)*chaplin::HPL(0,0,1,z) + 1728*pow(z,3)*chaplin::HPL(0,0,1,z) + 12096*pow(z,4)*chaplin::HPL(0,0,1,z) - 5184*pow(z,5)*chaplin::HPL(0,0,1,z) - 2592*chaplin::HPL(0,1,0,z) + 12096*z*chaplin::HPL(0,1,0,z) + 4320*pow(z,2)*chaplin::HPL(0,1,0,z) - 12096*pow(z,3)*chaplin::HPL(0,1,0,z) - 4104*pow(z,4)*chaplin::HPL(0,1,0,z) - 4968*pow(z,5)*chaplin::HPL(0,1,0,z) - 13824*chaplin::HPL(0,1,1,z) - 27648*pow(z,2)*chaplin::HPL(0,1,1,z) + 27648*pow(z,4)*chaplin::HPL(0,1,1,z) - 13824*pow(z,5)*chaplin::HPL(0,1,1,z) - 1728*chaplin::HPL(1,-1,0,z) - 1728*z*chaplin::HPL(1,-1,0,z) - 1728*pow(z,2)*chaplin::HPL(1,-1,0,z) + 1728*pow(z,3)*chaplin::HPL(1,-1,0,z) + 1728*pow(z,4)*chaplin::HPL(1,-1,0,z) + 1728*pow(z,5)*chaplin::HPL(1,-1,0,z) - 6048*chaplin::HPL(1,0,0,z) + 18144*z*chaplin::HPL(1,0,0,z) + 4320*pow(z,2)*chaplin::HPL(1,0,0,z) - 17280*pow(z,3)*chaplin::HPL(1,0,0,z) - 3456*pow(z,4)*chaplin::HPL(1,0,0,z) - 7776*pow(z,5)*chaplin::HPL(1,0,0,z) - 13824*chaplin::HPL(1,0,1,z) + 13824*z*chaplin::HPL(1,0,1,z) - 13824*pow(z,2)*chaplin::HPL(1,0,1,z) - 13824*pow(z,3)*chaplin::HPL(1,0,1,z) + 13824*pow(z,4)*chaplin::HPL(1,0,1,z) - 13824*pow(z,5)*chaplin::HPL(1,0,1,z) - 14472*chaplin::HPL(1,1,0,z) + 42120*z*chaplin::HPL(1,1,0,z) + 13176*pow(z,2)*chaplin::HPL(1,1,0,z) - 41688*pow(z,3)*chaplin::HPL(1,1,0,z) - 12744*pow(z,4)*chaplin::HPL(1,1,0,z) - 14472*pow(z,5)*chaplin::HPL(1,1,0,z) - 648*consts::z3 + 11232*z*consts::z3 + 4320*pow(z,2)*consts::z3 - 2808*pow(z,3)*consts::z3 + 3240*pow(z,4)*consts::z3 - 7128*pow(z,5)*consts::z3))/(432.*z*(-1 + pow(z,2))) + consts::nf*((pow(L,2)*(-4 - 3*z + 3*pow(z,2) + 4*pow(z,3) + pow(Nc,2)*(8 - 5*z + pow(z,2) - 8*pow(z,3))))/(24.*Nc*z) + consts::z2*(-((L*(-1 + pow(Nc,2))*(1 + z))/Nc) + (-4 - 11*z - 9*pow(z,2) + 36*pow(z,3) - 12*pow(z,4) + pow(Nc,2)*(8 + 9*z + 12*pow(z,2) - 47*pow(z,3) + 16*pow(z,4)) + 3*(-1 + z)*(2 + 2*(-9 + 8*pow(Nc,2))*z + (-15 + 16*pow(Nc,2))*pow(z,2))*chaplin::HPL(1,z))/(12.*Nc*(-1 + z)*z)) - (L*(16 - 73*pow(Nc,2) + 84*z + 7*pow(Nc,2)*z - 66*pow(z,2) - 5*pow(Nc,2)*pow(z,2) - 34*pow(z,3) + 91*pow(Nc,2)*pow(z,3) + 36*(-1 + pow(Nc,2))*z*(1 + z)*chaplin::HPL(0,1,z) + 72*(-1 + pow(Nc,2))*z*(1 + z)*chaplin::HPL(1,0,z)))/(36.*Nc*z) - (-72*(-4 - 5*z + 9*pow(z,2) + 8*pow(z,3) - 8*pow(z,4) + pow(Nc,2)*(-4 + 21*z - 33*pow(z,2) + 8*pow(z,3)))*chaplin::HPL(0,1,z) + 36*(4 + 23*z + 3*pow(z,2) - 50*pow(z,3) + 20*pow(z,4) + pow(Nc,2)*(16 - 77*z + 66*pow(z,2) + 13*pow(z,3)))*chaplin::HPL(1,0,z) + (-1 + z)*(-152 - 1643*pow(Nc,2) + 3108*z - 14*pow(Nc,2)*z - 2460*pow(z,2) - 149*pow(Nc,2)*pow(z,2) - 496*pow(z,3) + 2030*pow(Nc,2)*pow(z,3) + 432*(-1 + pow(Nc,2))*z*(1 + z)*chaplin::HPL(0,0,1,z) - 216*(-1 + pow(Nc,2))*z*(1 + z)*chaplin::HPL(0,1,0,z) - 864*z*chaplin::HPL(0,1,1,z) + 864*pow(Nc,2)*z*chaplin::HPL(0,1,1,z) - 864*pow(z,2)*chaplin::HPL(0,1,1,z) + 864*pow(Nc,2)*pow(z,2)*chaplin::HPL(0,1,1,z) - 216*chaplin::HPL(1,0,0,z) + 648*z*chaplin::HPL(1,0,0,z) - 432*pow(Nc,2)*z*chaplin::HPL(1,0,0,z) + 324*pow(z,2)*chaplin::HPL(1,0,0,z) - 432*pow(Nc,2)*pow(z,2)*chaplin::HPL(1,0,0,z) - 216*chaplin::HPL(1,1,0,z) + 1944*z*chaplin::HPL(1,1,0,z) - 1728*pow(Nc,2)*z*chaplin::HPL(1,1,0,z) + 1620*pow(z,2)*chaplin::HPL(1,1,0,z) - 1728*pow(Nc,2)*pow(z,2)*chaplin::HPL(1,1,0,z) + 432*consts::z3 - 864*z*consts::z3 + 432*pow(Nc,2)*z*consts::z3 - 216*pow(z,2)*consts::z3 + 432*pow(Nc,2)*pow(z,2)*consts::z3))/(432.*Nc*(-1 + z)*z));
    
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res);
}



double nnlo_r_lz0_logz(const double& z, const double& L)
{
    const double Nc = QCD::Nc;
    complex<double> res =  (2*pow(L,2)*pow(Nc,2)*(1 + 3*pow(z,2) - 4*pow(z,3) + pow(z,4)))/((-1 + z)*z) + (L*pow(Nc,2)*(77 - 179*z + 279*pow(z,2) - 353*pow(z,3) + 187*pow(z,4)))/(6.*(-1 + z)*z) - (pow(Nc,2)*(2 + 6*z + 18*pow(z,2) - 6*pow(z,3) - 19*pow(z,4) + 10*pow(z,5))*consts::z2)/(z*(-1 + pow(z,2))) - (pow(Nc,2)*(144*(-1 + z)*pow(1 + z + pow(z,2),2)*chaplin::HPL(-1,0,z) - (1 + z)*(2069 - 5083*z + 5655*pow(z,2) - 6427*pow(z,3) + 4054*pow(z,4) + 2304*z*(-1 + pow(z,2))*chaplin::HPL(1,0,z))))/(72.*z*(-1 + pow(z,2))) + consts::nf*((pow(L,2)*(-1 + pow(Nc,2))*(1 + z))/(4.*Nc) + (L*(4 + 8*z - 3*pow(z,2) - 17*pow(z,3) + 8*pow(z,4) - pow(Nc,2)*(8 + 4*z + 9*pow(z,2) - 29*pow(z,3) + 12*pow(z,4))))/(12.*Nc*(-1 + z)*z) - (3*(-1 + pow(Nc,2))*(1 + z)*consts::z2)/(2.*Nc) + (pow(Nc,2)*(-146 + 46*z + 3*pow(z,2) + 337*pow(z,3) - 280*pow(z,4)) + 8*(4 + 23*z - 45*pow(z,2) + pow(z,3) + 17*pow(z,4)) - 144*(-1 + pow(Nc,2))*z*(-1 + pow(z,2))*chaplin::HPL(1,0,z))/(72.*Nc*(-1 + z)*z));
    
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    return real(res)*log(z);
}

double nnlo_r_lz0_logz_sq(const double& z, const double& L)
{
    const double Nc = QCD::Nc;
    complex<double> res =  (pow(Nc,2)*(154 - 389*z + 699*pow(z,2) - 827*pow(z,3) + 374*pow(z,4)))/(24.*(-1 + z)*z) + (L*pow(Nc,2)*(2 + 3*z + 8*pow(z,2) - 3*pow(z,3) - 8*pow(z,4) + 3*pow(z,5)))/(z*(-1 + pow(z,2))) + consts::nf*((L*(-1 + pow(Nc,2))*(1 + z))/(2.*Nc) + (8 + 31*z - 18*pow(z,2) - 45*pow(z,3) + 24*pow(z,4) - pow(Nc,2)*(16 + 15*z + 24*pow(z,2) - 83*pow(z,3) + 32*pow(z,4)))/(48.*Nc*(-1 + z)*z));
    
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    return real(res)*pow(log(z),2.);
}

double nnlo_r_lz0_logz_cube(const double& z, const double& L)
{
    const double Nc = QCD::Nc;
    complex<double> res =  (5*(-1 + pow(Nc,2))*consts::nf*(1 + z))/(24.*Nc) + (pow(Nc,2)*(4 + 7*z + 18*pow(z,2) - 7*pow(z,3) - 17*pow(z,4) + 9*pow(z,5)))/(6.*z*(-1 + pow(z,2)));
    
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    return real(res)*pow(log(z),3.);
}

double nnlo_r_lz1(const double& z, const double& L)
{
    const double Nc=QCD::Nc;
    complex<double> res=
    (-4*pow(L,2)*pow(Nc,2)*(-1 + 2*z - pow(z,2) + pow(z,3)))/z + \
    (2*pow(Nc,2)*(-4 - 25*z - 56*pow(z,2) - 30*pow(z,3) + \
                  6*pow(z,4))*consts::z2)/(z*(1 + z)) - (L*pow(Nc,2)*(110 - 237*z + \
                                                                      243*pow(z,2) - 226*pow(z,3) + 110*pow(z,4) + 6*(13 - 10*z + \
                                                                                                                      39*pow(z,2) - 42*pow(z,3) + 13*pow(z,4))*chaplin::HPL(0,z)))/(3.*(-1 + z)*z) - \
    (pow(Nc,2)*(12*(182 - 191*z + 266*pow(z,2) - 128*pow(z,3) - \
                    420*pow(z,4) + 347*pow(z,5))*chaplin::HPL(0,z) + 36*(13 + 5*z + 33*pow(z,2) - \
                                                                5*pow(z,3) - 33*pow(z,4) + 15*pow(z,5))*pow(chaplin::HPL(0,z),2) - (-1 + \
                                                                                                                           z)*(144*pow(1 + z + pow(z,2),2)*chaplin::HPL(-1,0,z) - (1 + z)*(-2295 + 2159*z \
                                                                                                                                                                                  - 1894*pow(z,2) + 2295*pow(z,3) + 2304*z*(1 + \
                                                                                                                                                                                                                            z)*chaplin::HPL(1,0,z)))))/(36.*z*(-1 + pow(z,2))) + consts::nf*((4*(-1 + \
                                                                                                                                                                                                                                                                                        pow(Nc,2))*(1 + z)*consts::z2)/Nc + (L*(4 + 3*z - 3*pow(z,2) - \
                                                                                                                                                                                                                                                                                                                                4*pow(z,3) + pow(Nc,2)*(-8 + 5*z - pow(z,2) + 8*pow(z,3)) - 12*(-1 + \
                                                                                                                                                                                                                                                                                                                                                                                                pow(Nc,2))*z*(1 + z)*chaplin::HPL(0,z)))/(6.*Nc*z) + (6*(-8 - 13*z + \
                                                                                                                                                                                                                                                                                                                                                                                                                                                12*pow(z,2) + 25*pow(z,3) - 16*pow(z,4) + pow(Nc,2)*(4 + 25*z - \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     24*pow(z,2) - 21*pow(z,3) + 12*pow(z,4)))*chaplin::HPL(0,z) - 54*(-1 + \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              pow(Nc,2))*z*(-1 + pow(z,2))*pow(chaplin::HPL(0,z),2) + (-1 + z)*(4*(8 + 42*z \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          - 33*pow(z,2) - 17*pow(z,3)) + pow(Nc,2)*(-146 + 14*z - 13*pow(z,2) + \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    182*pow(z,3)) + 144*(-1 + pow(Nc,2))*z*(1 + \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            z)*chaplin::HPL(1,0,z)))/(36.*Nc*(-1 + z)*z))
    ;
    
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    return real(res)*log(1.-z);
}

double nnlo_r_lz2(const double& z, const double& L)
{
    const double Nc=QCD::Nc;
    complex<double> res=(12*L*pow(Nc,2)*(-1 + 2*z - pow(z,2) + pow(z,3)))/z + (consts::nf*(-4 \
                                                                                           - 3*z + 3*pow(z,2) + 4*pow(z,3) + pow(Nc,2)*(8 - 5*z + pow(z,2) - \
                                                                                                                                        8*pow(z,3)) + 12*(-1 + pow(Nc,2))*z*(1 + z)*chaplin::HPL(0,z)))/(6.*Nc*z) + \
    (pow(Nc,2)*(220 - 474*z + 486*pow(z,2) - 452*pow(z,3) + 220*pow(z,4) \
                + 3*(63 - 62*z + 189*pow(z,2) - 190*pow(z,3) + \
                     63*pow(z,4))*chaplin::HPL(0,z)))/(6.*(-1 + z)*z)
    ;
    
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    return real(res)*pow(log(1.-z),2.);
}

double nnlo_r_lz3(const double& z, const double& L)
{
    const double Nc=QCD::Nc;
    complex<double> res=(-8*pow(Nc,2)*(-1 + 2*z - pow(z,2) + pow(z,3)))/z;
    
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    return real(res)*pow(log(1.-z),3.);
}


// mathematica file for gg_n2lo_lzbar*_lz*_L* : Higgs_falko_gg_n2

double gg_n2lo_lzbar0_lz0_no_log(const double& z){
    complex<double> res =  (pow(QCD::Nc,2)*(-11 - 11*z - 11*pow(z,2) + 17*pow(z,3) + 8*pow(z,4) + 8*pow(z,5))*consts::z2*chaplin::HPL(-1,z))/(2.*(-1 + z)*z*(1 + z)) + ((16*pow(QCD::Nc,2)*consts::nf*(-1 + z)*z*pow(1 + z,2) + consts::nf*(-2 + 18*z + 17*pow(z,2) - 18*pow(z,3) - 15*pow(z,4)) + 2*pow(QCD::Nc,3)*(-7 + 127*z + 121*pow(z,2) - 125*pow(z,3) - 119*pow(z,4) + pow(z,5)))*consts::z2*chaplin::HPL(1,z))/(4.*QCD::Nc*(-1 + z)*z*(1 + z)) + ((1008*pow(QCD::Nc,3) + 864*pow(QCD::Nc,3)*z - 2520*pow(QCD::Nc,3)*pow(z,2) - 1656*pow(QCD::Nc,3)*pow(z,3) + 1512*pow(QCD::Nc,3)*pow(z,4) + 792*pow(QCD::Nc,3)*pow(z,5))*chaplin::HPL(-1,0,z))/(432.*QCD::Nc*(-1 + z)*z*(1 + z)) + ((-13536*pow(QCD::Nc,3) - 288*consts::nf - 288*pow(QCD::Nc,2)*consts::nf + 9792*pow(QCD::Nc,3)*z - 648*consts::nf*z + 1224*pow(QCD::Nc,2)*consts::nf*z - 23760*pow(QCD::Nc,3)*pow(z,2) + 288*consts::nf*pow(z,2) - 864*pow(QCD::Nc,2)*consts::nf*pow(z,2) + 10800*pow(QCD::Nc,3)*pow(z,3) + 1224*consts::nf*pow(z,3) - 1800*pow(QCD::Nc,2)*consts::nf*pow(z,3) + 34848*pow(QCD::Nc,3)*pow(z,4) + 576*pow(QCD::Nc,2)*consts::nf*pow(z,4) - 23040*pow(QCD::Nc,3)*pow(z,5) - 576*consts::nf*pow(z,5))*chaplin::HPL(0,1,z))/(432.*QCD::Nc*(-1 + z)*z*(1 + z)) + ((-3132*pow(QCD::Nc,3) - 144*consts::nf - 576*pow(QCD::Nc,2)*consts::nf + 6012*pow(QCD::Nc,3)*z - 972*consts::nf*z + 2196*pow(QCD::Nc,2)*consts::nf*z - 38700*pow(QCD::Nc,3)*pow(z,2) - 936*consts::nf*pow(z,2) + 396*pow(QCD::Nc,2)*consts::nf*pow(z,2) + 20124*pow(QCD::Nc,3)*pow(z,3) + 1692*consts::nf*pow(z,3) - 2844*pow(QCD::Nc,2)*consts::nf*pow(z,3) + 39780*pow(QCD::Nc,3)*pow(z,4) + 1080*consts::nf*pow(z,4) - 468*pow(QCD::Nc,2)*consts::nf*pow(z,4) - 28188*pow(QCD::Nc,3)*pow(z,5) - 720*consts::nf*pow(z,5))*chaplin::HPL(1,0,z))/(432.*QCD::Nc*(-1 + z)*z*(1 + z)) + ((-1296*pow(QCD::Nc,3) - 1296*pow(QCD::Nc,3)*z - 1296*pow(QCD::Nc,3)*pow(z,2) + 3888*pow(QCD::Nc,3)*pow(z,3))*chaplin::HPL(-1,-1,0,z))/(432.*QCD::Nc*(-1 + z)*z*(1 + z)) + ((216*pow(QCD::Nc,3) + 216*pow(QCD::Nc,3)*z + 216*pow(QCD::Nc,3)*pow(z,2) - 2376*pow(QCD::Nc,3)*pow(z,3) + 864*pow(QCD::Nc,3)*pow(z,4) + 864*pow(QCD::Nc,3)*pow(z,5))*chaplin::HPL(-1,0,0,z))/(432.*QCD::Nc*(-1 + z)*z*(1 + z)) + ((-1728*pow(QCD::Nc,3) - 1728*pow(QCD::Nc,3)*z - 1728*pow(QCD::Nc,3)*pow(z,2) + 1728*pow(QCD::Nc,3)*pow(z,3) + 1728*pow(QCD::Nc,3)*pow(z,4) + 1728*pow(QCD::Nc,3)*pow(z,5))*chaplin::HPL(-1,1,0,z))/(432.*QCD::Nc*(-1 + z)*z*(1 + z)) + ((432*pow(QCD::Nc,3)*pow(z,4) - 432*pow(QCD::Nc,3)*pow(z,5))*chaplin::HPL(0,-1,0,z))/(432.*QCD::Nc*(-1 + z)*z*(1 + z)) + ((-5184*pow(QCD::Nc,3) - 1728*pow(QCD::Nc,3)*z - 432*consts::nf*z + 432*pow(QCD::Nc,2)*consts::nf*z - 12096*pow(QCD::Nc,3)*pow(z,2) - 432*consts::nf*pow(z,2) + 432*pow(QCD::Nc,2)*consts::nf*pow(z,2) + 1728*pow(QCD::Nc,3)*pow(z,3) + 432*consts::nf*pow(z,3) - 432*pow(QCD::Nc,2)*consts::nf*pow(z,3) + 12096*pow(QCD::Nc,3)*pow(z,4) + 432*consts::nf*pow(z,4) - 432*pow(QCD::Nc,2)*consts::nf*pow(z,4) - 5184*pow(QCD::Nc,3)*pow(z,5))*chaplin::HPL(0,0,1,z))/(432.*QCD::Nc*(-1 + z)*z*(1 + z)) + ((-2592*pow(QCD::Nc,3) + 12096*pow(QCD::Nc,3)*z + 216*consts::nf*z - 216*pow(QCD::Nc,2)*consts::nf*z + 4320*pow(QCD::Nc,3)*pow(z,2) + 216*consts::nf*pow(z,2) - 216*pow(QCD::Nc,2)*consts::nf*pow(z,2) - 12096*pow(QCD::Nc,3)*pow(z,3) - 216*consts::nf*pow(z,3) + 216*pow(QCD::Nc,2)*consts::nf*pow(z,3) - 4104*pow(QCD::Nc,3)*pow(z,4) - 216*consts::nf*pow(z,4) + 216*pow(QCD::Nc,2)*consts::nf*pow(z,4) - 4968*pow(QCD::Nc,3)*pow(z,5))*chaplin::HPL(0,1,0,z))/(432.*QCD::Nc*(-1 + z)*z*(1 + z)) + ((-13824*pow(QCD::Nc,3) - 864*consts::nf*z + 864*pow(QCD::Nc,2)*consts::nf*z - 27648*pow(QCD::Nc,3)*pow(z,2) - 864*consts::nf*pow(z,2) + 864*pow(QCD::Nc,2)*consts::nf*pow(z,2) + 864*consts::nf*pow(z,3) - 864*pow(QCD::Nc,2)*consts::nf*pow(z,3) + 27648*pow(QCD::Nc,3)*pow(z,4) + 864*consts::nf*pow(z,4) - 864*pow(QCD::Nc,2)*consts::nf*pow(z,4) - 13824*pow(QCD::Nc,3)*pow(z,5))*chaplin::HPL(0,1,1,z))/(432.*QCD::Nc*(-1 + z)*z*(1 + z)) + ((-1728*pow(QCD::Nc,3) - 1728*pow(QCD::Nc,3)*z - 1728*pow(QCD::Nc,3)*pow(z,2) + 1728*pow(QCD::Nc,3)*pow(z,3) + 1728*pow(QCD::Nc,3)*pow(z,4) + 1728*pow(QCD::Nc,3)*pow(z,5))*chaplin::HPL(1,-1,0,z))/(432.*QCD::Nc*(-1 + z)*z*(1 + z)) + ((-6048*pow(QCD::Nc,3) - 216*consts::nf + 18144*pow(QCD::Nc,3)*z + 648*consts::nf*z - 432*pow(QCD::Nc,2)*consts::nf*z + 4320*pow(QCD::Nc,3)*pow(z,2) + 540*consts::nf*pow(z,2) - 432*pow(QCD::Nc,2)*consts::nf*pow(z,2) - 17280*pow(QCD::Nc,3)*pow(z,3) - 648*consts::nf*pow(z,3) + 432*pow(QCD::Nc,2)*consts::nf*pow(z,3) - 3456*pow(QCD::Nc,3)*pow(z,4) - 324*consts::nf*pow(z,4) + 432*pow(QCD::Nc,2)*consts::nf*pow(z,4) - 7776*pow(QCD::Nc,3)*pow(z,5))*chaplin::HPL(1,0,0,z))/(432.*QCD::Nc*(-1 + z)*z*(1 + z)) + ((-13824*pow(QCD::Nc,3) + 13824*pow(QCD::Nc,3)*z - 13824*pow(QCD::Nc,3)*pow(z,2) - 13824*pow(QCD::Nc,3)*pow(z,3) + 13824*pow(QCD::Nc,3)*pow(z,4) - 13824*pow(QCD::Nc,3)*pow(z,5))*chaplin::HPL(1,0,1,z))/(432.*QCD::Nc*(-1 + z)*z*(1 + z)) + ((-14472*pow(QCD::Nc,3) - 216*consts::nf + 42120*pow(QCD::Nc,3)*z + 1944*consts::nf*z - 1728*pow(QCD::Nc,2)*consts::nf*z + 13176*pow(QCD::Nc,3)*pow(z,2) + 1836*consts::nf*pow(z,2) - 1728*pow(QCD::Nc,2)*consts::nf*pow(z,2) - 41688*pow(QCD::Nc,3)*pow(z,3) - 1944*consts::nf*pow(z,3) + 1728*pow(QCD::Nc,2)*consts::nf*pow(z,3) - 12744*pow(QCD::Nc,3)*pow(z,4) - 1620*consts::nf*pow(z,4) + 1728*pow(QCD::Nc,2)*consts::nf*pow(z,4) - 14472*pow(QCD::Nc,3)*pow(z,5))*chaplin::HPL(1,1,0,z))/(432.*QCD::Nc*(-1 + z)*z*(1 + z)) + (18321*pow(QCD::Nc,3) - 152*consts::nf - 1643*pow(QCD::Nc,2)*consts::nf - 18523*pow(QCD::Nc,3)*z + 3108*consts::nf*z - 14*pow(QCD::Nc,2)*consts::nf*z + 1430*pow(QCD::Nc,3)*pow(z,2) - 2308*consts::nf*pow(z,2) + 1494*pow(QCD::Nc,2)*consts::nf*pow(z,2) - 2642*pow(QCD::Nc,3)*pow(z,3) - 3604*consts::nf*pow(z,3) + 2044*pow(QCD::Nc,2)*consts::nf*pow(z,3) - 19751*pow(QCD::Nc,3)*pow(z,4) + 2460*consts::nf*pow(z,4) + 149*pow(QCD::Nc,2)*consts::nf*pow(z,4) + 21165*pow(QCD::Nc,3)*pow(z,5) + 496*consts::nf*pow(z,5) - 2030*pow(QCD::Nc,2)*consts::nf*pow(z,5) - 2160*pow(QCD::Nc,3)*consts::z2 - 144*consts::nf*consts::z2 + 288*pow(QCD::Nc,2)*consts::nf*consts::z2 + 6624*pow(QCD::Nc,3)*z*consts::z2 - 540*consts::nf*z*consts::z2 + 612*pow(QCD::Nc,2)*consts::nf*z*consts::z2 - 11520*pow(QCD::Nc,3)*pow(z,2)*consts::z2 - 720*consts::nf*pow(z,2)*consts::z2 + 756*pow(QCD::Nc,2)*consts::nf*pow(z,2)*consts::z2 + 11592*pow(QCD::Nc,3)*pow(z,3)*consts::z2 + 972*consts::nf*pow(z,3)*consts::z2 - 1260*pow(QCD::Nc,2)*consts::nf*pow(z,3)*consts::z2 + 14076*pow(QCD::Nc,3)*pow(z,4)*consts::z2 + 864*consts::nf*pow(z,4)*consts::z2 - 1116*pow(QCD::Nc,2)*consts::nf*pow(z,4)*consts::z2 - 17820*pow(QCD::Nc,3)*pow(z,5)*consts::z2 - 432*consts::nf*pow(z,5)*consts::z2 + 576*pow(QCD::Nc,2)*consts::nf*pow(z,5)*consts::z2 - 648*pow(QCD::Nc,3)*consts::z3 + 432*consts::nf*consts::z3 + 11232*pow(QCD::Nc,3)*z*consts::z3 - 864*consts::nf*z*consts::z3 + 432*pow(QCD::Nc,2)*consts::nf*z*consts::z3 + 4320*pow(QCD::Nc,3)*pow(z,2)*consts::z3 - 648*consts::nf*pow(z,2)*consts::z3 + 432*pow(QCD::Nc,2)*consts::nf*pow(z,2)*consts::z3 - 2808*pow(QCD::Nc,3)*pow(z,3)*consts::z3 + 864*consts::nf*pow(z,3)*consts::z3 - 432*pow(QCD::Nc,2)*consts::nf*pow(z,3)*consts::z3 + 3240*pow(QCD::Nc,3)*pow(z,4)*consts::z3 + 216*consts::nf*pow(z,4)*consts::z3 - 432*pow(QCD::Nc,2)*consts::nf*pow(z,4)*consts::z3 - 7128*pow(QCD::Nc,3)*pow(z,5)*consts::z3)/(432.*QCD::Nc*(-1 + z)*z*(1 + z));
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res);
}


double gg_n2lo_lzbar0_lz1_no_log(const double& z){
    const double Nc=QCD::Nc;
    complex<double> res =
    (2069*pow(QCD::Nc,3) + 32*consts::nf - 146*pow(QCD::Nc,2)*consts::nf - 3014*pow(QCD::Nc,3)*z + 216*consts::nf*z - 100*pow(QCD::Nc,2)*consts::nf*z + 572*pow(QCD::Nc,3)*pow(z,2) - 176*consts::nf*pow(z,2) + 49*pow(QCD::Nc,2)*consts::nf*pow(z,2) - 772*pow(QCD::Nc,3)*pow(z,3) - 352*consts::nf*pow(z,3) + 340*pow(QCD::Nc,2)*consts::nf*pow(z,3) - 2373*pow(QCD::Nc,3)*pow(z,4) + 144*consts::nf*pow(z,4) + 57*pow(QCD::Nc,2)*consts::nf*pow(z,4) + 4054*pow(QCD::Nc,3)*pow(z,5) + 136*consts::nf*pow(z,5) - 280*pow(QCD::Nc,2)*consts::nf*pow(z,5) - 144*pow(QCD::Nc,3)*consts::z2 - 432*pow(QCD::Nc,3)*z*consts::z2 - 108*consts::nf*z*consts::z2 + 108*pow(QCD::Nc,2)*consts::nf*z*consts::z2 - 1296*pow(QCD::Nc,3)*pow(z,2)*consts::z2 - 108*consts::nf*pow(z,2)*consts::z2 + 108*pow(QCD::Nc,2)*consts::nf*pow(z,2)*consts::z2 + 432*pow(QCD::Nc,3)*pow(z,3)*consts::z2 + 108*consts::nf*pow(z,3)*consts::z2 - 108*pow(QCD::Nc,2)*consts::nf*pow(z,3)*consts::z2 + 1368*pow(QCD::Nc,3)*pow(z,4)*consts::z2 + 108*consts::nf*pow(z,4)*consts::z2 - 108*pow(QCD::Nc,2)*consts::nf*pow(z,4)*consts::z2 - 720*pow(QCD::Nc,3)*pow(z,5)*consts::z2)/(72.*QCD::Nc*(-1 + z)*z*(1 + z)) - (2*pow(QCD::Nc,2)*pow(1 + z + pow(z,2),2)*chaplin::HPL(-1,0,z))/(z*(1 + z)) + (2*(16*pow(QCD::Nc,3) + consts::nf - pow(QCD::Nc,2)*consts::nf)*(1 + z)*chaplin::HPL(1,0,z))/Nc;
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res)*log(z);
}

double gg_n2lo_lzbar0_lz2_no_log(const double& z){
    const double Nc=QCD::Nc;
    complex<double> res =
    (308*pow(Nc,3) + 8*consts::nf - 16*pow(Nc,2)*consts::nf - 470*pow(Nc,3)*z + 39*consts::nf*z - 31*pow(Nc,2)*consts::nf*z + 620*pow(Nc,3)*pow(z,2) + 13*consts::nf*pow(z,2) - 39*pow(Nc,2)*consts::nf*pow(z,2) - 256*pow(Nc,3)*pow(z,3) - 63*consts::nf*pow(z,3) + 59*pow(Nc,2)*consts::nf*pow(z,3) - 906*pow(Nc,3)*pow(z,4) - 21*consts::nf*pow(z,4) + 51*pow(Nc,2)*consts::nf*pow(z,4) + 748*pow(Nc,3)*pow(z,5) + 24*consts::nf*pow(z,5) - 32*pow(Nc,2)*consts::nf*pow(z,5))/(48.*Nc*(-1 + z)*z*(1 + z));
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res)*pow(log(z),2.);
}

double gg_n2lo_lzbar0_lz3_no_log(const double& z){
    const double Nc=QCD::Nc;
    complex<double> res =
    (-5*consts::nf*(-1 + z)*z*pow(1 + z,2) + 5*pow(Nc,2)*consts::nf*(-1 + z)*z*pow(1 + z,2) + 4*pow(Nc,3)*(4 + 7*z + 18*pow(z,2) - 7*pow(z,3) - 17*pow(z,4) + 9*pow(z,5)))/(24.*Nc*(-1 + z)*z*(1 + z));
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res)*pow(log(z),3.);
}

double gg_n2lo_lzbar0_lz0_L_1(const double& z){
    const double Nc=3.;
    
    complex<double> res = (12966*pow(Nc,3) + 192*consts::nf - 876*pow(Nc,2)*consts::nf - 10686*pow(Nc,3)*z + 1008*consts::nf*z + 84*pow(Nc,2)*consts::nf*z - 3888*pow(Nc,3)*pow(z,2) - 984*consts::nf*pow(z,2) + 816*pow(Nc,2)*consts::nf*pow(z,2) - 2280*pow(Nc,3)*pow(z,3) - 1416*consts::nf*pow(z,3) + 1008*pow(Nc,2)*consts::nf*pow(z,3) - 9078*pow(Nc,3)*pow(z,4) + 792*consts::nf*pow(z,4) + 60*pow(Nc,2)*consts::nf*pow(z,4) + 12966*pow(Nc,3)*pow(z,5) + 408*consts::nf*pow(z,5) - 1092*pow(Nc,2)*consts::nf*pow(z,5) - 1728*pow(Nc,3)*consts::z2 - 2160*pow(Nc,3)*z*consts::z2 - 432*consts::nf*z*consts::z2 + 432*pow(Nc,2)*consts::nf*z*consts::z2 - 6480*pow(Nc,3)*pow(z,2)*consts::z2 - 432*consts::nf*pow(z,2)*consts::z2 + 432*pow(Nc,2)*consts::nf*pow(z,2)*consts::z2 + 4320*pow(Nc,3)*pow(z,3)*consts::z2 + 432*consts::nf*pow(z,3)*consts::z2 - 432*pow(Nc,2)*consts::nf*pow(z,3)*consts::z2 + 8640*pow(Nc,3)*pow(z,4)*consts::z2 + 432*consts::nf*pow(z,4)*consts::z2 - 432*pow(Nc,2)*consts::nf*pow(z,4)*consts::z2 - 2592*pow(Nc,3)*pow(z,5)*consts::z2)/(432.*Nc*(-1 + z)*z*(1 + z)) + ((864*pow(Nc,3) + 864*pow(Nc,3)*z + 864*pow(Nc,3)*pow(z,2) - 864*pow(Nc,3)*pow(z,3) - 864*pow(Nc,3)*pow(z,4) - 864*pow(Nc,3)*pow(z,5))*chaplin::HPL(-1,0,z))/(432.*Nc*(-1 + z)*z*(1 + z)) + ((-5184*pow(Nc,3) - 1728*pow(Nc,3)*z - 432*consts::nf*z + 432*pow(Nc,2)*consts::nf*z - 12096*pow(Nc,3)*pow(z,2) - 432*consts::nf*pow(z,2) + 432*pow(Nc,2)*consts::nf*pow(z,2) + 1728*pow(Nc,3)*pow(z,3) + 432*consts::nf*pow(z,3) - 432*pow(Nc,2)*consts::nf*pow(z,3) + 12096*pow(Nc,3)*pow(z,4) + 432*consts::nf*pow(z,4) - 432*pow(Nc,2)*consts::nf*pow(z,4) - 5184*pow(Nc,3)*pow(z,5))*chaplin::HPL(0,1,z))/(432.*Nc*(-1 + z)*z*(1 + z)) + ((-5184*pow(Nc,3) - 8640*pow(Nc,3)*z - 864*consts::nf*z + 864*pow(Nc,2)*consts::nf*z - 19008*pow(Nc,3)*pow(z,2) - 864*consts::nf*pow(z,2) + 864*pow(Nc,2)*consts::nf*pow(z,2) + 8640*pow(Nc,3)*pow(z,3) + 864*consts::nf*pow(z,3) - 864*pow(Nc,2)*consts::nf*pow(z,3) + 19008*pow(Nc,3)*pow(z,4) + 864*consts::nf*pow(z,4) - 864*pow(Nc,2)*consts::nf*pow(z,4) - 5184*pow(Nc,3)*pow(z,5))*chaplin::HPL(1,0,z))/(432.*Nc*(-1 + z)*z*(1 + z));
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res);
}

double gg_n2lo_lzbar0_lz1_L_1(const double& z){
    const double Nc=QCD::Nc;
    complex<double> res =
    (924*pow(Nc,3) + 24*consts::nf - 48*pow(Nc,2)*consts::nf - 1224*pow(Nc,3)*z + 72*consts::nf*z - 72*pow(Nc,2)*consts::nf*z + 1200*pow(Nc,3)*pow(z,2) + 30*consts::nf*pow(z,2) - 78*pow(Nc,2)*consts::nf*pow(z,2) - 888*pow(Nc,3)*pow(z,3) - 120*consts::nf*pow(z,3) + 120*pow(Nc,2)*consts::nf*pow(z,3) - 1992*pow(Nc,3)*pow(z,4) - 54*consts::nf*pow(z,4) + 102*pow(Nc,2)*consts::nf*pow(z,4) + 2244*pow(Nc,3)*pow(z,5) + 48*consts::nf*pow(z,5) - 72*pow(Nc,2)*consts::nf*pow(z,5))/(72.*Nc*(-1 + z)*z*(1 + z));
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res)*log(z);
}

double gg_n2lo_lzbar0_lz2_L_1(const double& z){
    const double Nc=QCD::Nc;
    complex<double> res =
    (96*pow(Nc,3) + 144*pow(Nc,3)*z + 24*consts::nf*z - 24*pow(Nc,2)*consts::nf*z + 384*pow(Nc,3)*pow(z,2) + 24*consts::nf*pow(z,2) - 24*pow(Nc,2)*consts::nf*pow(z,2) - 144*pow(Nc,3)*pow(z,3) - 24*consts::nf*pow(z,3) + 24*pow(Nc,2)*consts::nf*pow(z,3) - 384*pow(Nc,3)*pow(z,4) - 24*consts::nf*pow(z,4) + 24*pow(Nc,2)*consts::nf*pow(z,4) + 144*pow(Nc,3)*pow(z,5))/(48.*Nc*(-1 + z)*z*(1 + z));
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res)*pow(log(z),2.);
}


double gg_n2lo_lzbar0_lz0_L_2(const double& z){
    const double Nc=3.;
    
    complex<double> res = (3564*pow(Nc,3) + 72*consts::nf - 144*pow(Nc,2)*consts::nf - 3384*pow(Nc,3)*z + 54*consts::nf*z + 90*pow(Nc,2)*consts::nf*z - 576*pow(Nc,3)*pow(z,2) - 126*consts::nf*pow(z,2) + 126*pow(Nc,2)*consts::nf*pow(z,2) - 180*pow(Nc,3)*pow(z,3) - 126*consts::nf*pow(z,3) + 54*pow(Nc,2)*consts::nf*pow(z,3) - 2988*pow(Nc,3)*pow(z,4) + 54*consts::nf*pow(z,4) + 18*pow(Nc,2)*consts::nf*pow(z,4) + 3564*pow(Nc,3)*pow(z,5) + 72*consts::nf*pow(z,5) - 144*pow(Nc,2)*consts::nf*pow(z,5))/(432.*Nc*(-1 + z)*z*(1 + z));
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res);
}

double gg_n2lo_lzbar0_lz1_L_2(const double& z){
    const double Nc=QCD::Nc;
    complex<double> res =
    (144*pow(Nc,3) + 144*pow(Nc,3)*z + 18*consts::nf*z - 18*pow(Nc,2)*consts::nf*z + 432*pow(Nc,3)*pow(z,2) + 18*consts::nf*pow(z,2) - 18*pow(Nc,2)*consts::nf*pow(z,2) - 144*pow(Nc,3)*pow(z,3) - 18*consts::nf*pow(z,3) + 18*pow(Nc,2)*consts::nf*pow(z,3) - 432*pow(Nc,3)*pow(z,4) - 18*consts::nf*pow(z,4) + 18*pow(Nc,2)*consts::nf*pow(z,4) + 144*pow(Nc,3)*pow(z,5))/(72.*Nc*(-1 + z)*z*(1 + z));
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res)*log(z);
}



double gg_n2lo_lzbar1_L_0(const double& z){
    const double Nc=QCD::Nc;
    complex<double> res =
    -(2295*pow(Nc,3) + 32*consts::nf - 146*pow(Nc,2)*consts::nf - 2159*pow(Nc,3)*z + 168*consts::nf*z + 14*pow(Nc,2)*consts::nf*z - 401*pow(Nc,3)*pow(z,2) - 164*consts::nf*pow(z,2) + 133*pow(Nc,2)*consts::nf*pow(z,2) - 136*pow(Nc,3)*pow(z,3) - 236*consts::nf*pow(z,3) + 168*pow(Nc,2)*consts::nf*pow(z,3) - 1894*pow(Nc,3)*pow(z,4) + 132*consts::nf*pow(z,4) + 13*pow(Nc,2)*consts::nf*pow(z,4) + 2295*pow(Nc,3)*pow(z,5) + 68*consts::nf*pow(z,5) - 182*pow(Nc,2)*consts::nf*pow(z,5) - 288*pow(Nc,3)*consts::z2 - 1512*pow(Nc,3)*z*consts::z2 - 144*consts::nf*z*consts::z2 + 144*pow(Nc,2)*consts::nf*z*consts::z2 - 2232*pow(Nc,3)*pow(z,2)*consts::z2 - 144*consts::nf*pow(z,2)*consts::z2 + 144*pow(Nc,2)*consts::nf*pow(z,2)*consts::z2 + 1872*pow(Nc,3)*pow(z,3)*consts::z2 + 144*consts::nf*pow(z,3)*consts::z2 - 144*pow(Nc,2)*consts::nf*pow(z,3)*consts::z2 + 2592*pow(Nc,3)*pow(z,4)*consts::z2 + 144*consts::nf*pow(z,4)*consts::z2 - 144*pow(Nc,2)*consts::nf*pow(z,4)*consts::z2 - 432*pow(Nc,3)*pow(z,5)*consts::z2)/(36.*Nc*(-1 + z)*z*(1 + z)) - ((2184*pow(Nc,3) + 48*consts::nf - 24*pow(Nc,2)*consts::nf - 2292*pow(Nc,3)*z + 126*consts::nf*z - 174*pow(Nc,2)*consts::nf*z + 3192*pow(Nc,3)*pow(z,2) + 6*consts::nf*pow(z,2) - 6*pow(Nc,2)*consts::nf*pow(z,2) - 1536*pow(Nc,3)*pow(z,3) - 222*consts::nf*pow(z,3) + 270*pow(Nc,2)*consts::nf*pow(z,3) - 5040*pow(Nc,3)*pow(z,4) - 54*consts::nf*pow(z,4) + 54*pow(Nc,2)*consts::nf*pow(z,4) + 4164*pow(Nc,3)*pow(z,5) + 96*consts::nf*pow(z,5) - 72*pow(Nc,2)*consts::nf*pow(z,5))*chaplin::HPL(0,z))/(36.*Nc*(-1 + z)*z*(1 + z)) - ((468*pow(Nc,3) + 180*pow(Nc,3)*z + 54*consts::nf*z - 54*pow(Nc,2)*consts::nf*z + 1188*pow(Nc,3)*pow(z,2) + 54*consts::nf*pow(z,2) - 54*pow(Nc,2)*consts::nf*pow(z,2) - 180*pow(Nc,3)*pow(z,3) - 54*consts::nf*pow(z,3) + 54*pow(Nc,2)*consts::nf*pow(z,3) - 1188*pow(Nc,3)*pow(z,4) - 54*consts::nf*pow(z,4) + 54*pow(Nc,2)*consts::nf*pow(z,4) + 540*pow(Nc,3)*pow(z,5))*pow(chaplin::HPL(0,z),2))/(36.*Nc*(-1 + z)*z*(1 + z)) - ((144*pow(Nc,3) + 144*pow(Nc,3)*z + 144*pow(Nc,3)*pow(z,2) - 144*pow(Nc,3)*pow(z,3) - 144*pow(Nc,3)*pow(z,4) - 144*pow(Nc,3)*pow(z,5))*chaplin::HPL(-1,0,z))/(36.*Nc*(-1 + z)*z*(1 + z)) - ((-2304*pow(Nc,3)*z - 144*consts::nf*z + 144*pow(Nc,2)*consts::nf*z - 2304*pow(Nc,3)*pow(z,2) - 144*consts::nf*pow(z,2) + 144*pow(Nc,2)*consts::nf*pow(z,2) + 2304*pow(Nc,3)*pow(z,3) + 144*consts::nf*pow(z,3) - 144*pow(Nc,2)*consts::nf*pow(z,3) + 2304*pow(Nc,3)*pow(z,4) + 144*consts::nf*pow(z,4) - 144*pow(Nc,2)*consts::nf*pow(z,4))*chaplin::HPL(1,0,z))/(36.*Nc*(-1 + z)*z*(1 + z));
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res)*log(1.-z);
}

double gg_n2lo_lzbar1_L_1(const double& z){
    const double Nc=QCD::Nc;
    complex<double> res =
    -(1320*pow(Nc,3) + 24*consts::nf - 48*pow(Nc,2)*consts::nf - 1524*pow(Nc,3)*z + 18*consts::nf*z + 30*pow(Nc,2)*consts::nf*z + 72*pow(Nc,3)*pow(z,2) - 42*consts::nf*pow(z,2) + 42*pow(Nc,2)*consts::nf*pow(z,2) + 204*pow(Nc,3)*pow(z,3) - 42*consts::nf*pow(z,3) + 18*pow(Nc,2)*consts::nf*pow(z,3) - 1392*pow(Nc,3)*pow(z,4) + 18*consts::nf*pow(z,4) + 6*pow(Nc,2)*consts::nf*pow(z,4) + 1320*pow(Nc,3)*pow(z,5) + 24*consts::nf*pow(z,5) - 48*pow(Nc,2)*consts::nf*pow(z,5))/(36.*Nc*(-1 + z)*z*(1 + z)) - ((936*pow(Nc,3) + 216*pow(Nc,3)*z + 72*consts::nf*z - 72*pow(Nc,2)*consts::nf*z + 2088*pow(Nc,3)*pow(z,2) + 72*consts::nf*pow(z,2) - 72*pow(Nc,2)*consts::nf*pow(z,2) - 216*pow(Nc,3)*pow(z,3) - 72*consts::nf*pow(z,3) + 72*pow(Nc,2)*consts::nf*pow(z,3) - 2088*pow(Nc,3)*pow(z,4) - 72*consts::nf*pow(z,4) + 72*pow(Nc,2)*consts::nf*pow(z,4) + 936*pow(Nc,3)*pow(z,5))*chaplin::HPL(0,z))/(36.*Nc*(-1 + z)*z*(1 + z));
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res)*log(1.-z);
}

double gg_n2lo_lzbar1_L_2(const double& z){
    const double Nc=QCD::Nc;
    complex<double> res =
    -(144*pow(Nc,3) - 288*pow(Nc,3)*z + 144*pow(Nc,3)*pow(z,3) - 144*pow(Nc,3)*pow(z,4) + 144*pow(Nc,3)*pow(z,5))/(36.*Nc*(-1 + z)*z*(1 + z));
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res)*log(1.-z);
}

double gg_n2lo_lzbar2_L_0(const double& z){
    const double Nc=QCD::Nc;
    complex<double> res = (220*pow(Nc,3) + 4*consts::nf - 8*pow(Nc,2)*consts::nf - 474*pow(Nc,3)*z - consts::nf*z + 13*pow(Nc,2)*consts::nf*z + 486*pow(Nc,3)*pow(z,2) - 6*consts::nf*pow(z,2) - 6*pow(Nc,2)*consts::nf*pow(z,2) - 452*pow(Nc,3)*pow(z,3) - consts::nf*pow(z,3) + 9*pow(Nc,2)*consts::nf*pow(z,3) + 220*pow(Nc,3)*pow(z,4) + 4*consts::nf*pow(z,4) - 8*pow(Nc,2)*consts::nf*pow(z,4))/(6.*Nc*(-1 + z)*z) + ((189*pow(Nc,3) - 186*pow(Nc,3)*z + 12*consts::nf*z - 12*pow(Nc,2)*consts::nf*z + 567*pow(Nc,3)*pow(z,2) - 570*pow(Nc,3)*pow(z,3) - 12*consts::nf*pow(z,3) + 12*pow(Nc,2)*consts::nf*pow(z,3) + 189*pow(Nc,3)*pow(z,4))*chaplin::HPL(0,z))/(6.*Nc*(-1 + z)*z)
    ;
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res)*pow(log(1.-z),2.);
}

double gg_n2lo_lzbar2_L_1(const double& z){
    const double Nc=QCD::Nc;
    complex<double> res = (72*pow(Nc,3) - 216*pow(Nc,3)*z + 216*pow(Nc,3)*pow(z,2) - 144*pow(Nc,3)*pow(z,3) + 72*pow(Nc,3)*pow(z,4))/(6.*Nc*(-1 + z)*z)
    ;
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res)*pow(log(1.-z),2.);
}

double gg_n2lo_lzbar3_L_0(const double& z){
    const double Nc=QCD::Nc;
    complex<double> res = (-8*pow(Nc,2)*(-1 + 2*z - pow(z,2) + pow(z,3)))/z
    ;
    check_imaginary_part(res,__PRETTY_FUNCTION__);
    
    
    return real(res)*pow(log(1.-z),3.);
}

}
